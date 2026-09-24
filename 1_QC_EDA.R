# ==============================================================================
#  1. QC and exploratory analysis of the GeoMx DSP data
# ------------------------------------------------------------------------------
#  Shows how QC and EDA were done for the paper. The data are not included,
#  so the inputs below are left empty.
#
#  Steps
#    1) Build a SpatialExperiment from the GeoMx count matrix + annotations
#    2) QC
#       2.1) Gene QC  : low-expression genes flagged per AOI (standR)
#       2.2) AOI QC   : raw reads, read processing, saturation, nuclei, size
#       2.3) LOQ      : gene-level filter based on negative-control probes
#    3) EDA
#       3.1) PCA and UMAP coloured by each design variable
#       3.2) RLE plots (raw, logCPM, TMM) to choose the normalisation
# ==============================================================================

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(standR)
  library(edgeR)
  library(scater)
  library(ggplot2)
  library(patchwork)
})

set.seed(100)  # PCA / UMAP

# ---- 0. Inputs ---------------------------------------------------------------
# count_data   : data.frame, first column = gene symbol ("TargetName"),
#                remaining columns = one per AOI (raw probe counts)
# sample_anno  : data.frame, one row per AOI with the GeoMx segment properties
#                (RawReads, AlignedReads, SequencingSaturation, AOINucleiCount,
#                ...) plus the study design variables:
#                  subject    : donor / patient
#                  biomaterial: graft material
#                  tissue     : tissue compartment defined by the
#                                morphology marker (mineralized / connective)
#                  stage      : sampling depth / maturation stage
# feature_anno : data.frame with probe/gene annotation (incl. negative probes)
count_data   <- NULL
sample_anno  <- NULL
feature_anno <- NULL

spe <- readGeoMx(count_data, sample_anno, feature_anno)

# ---- 1. Gene QC --------------------------------------------------------------
# Flags genes with low expression in > 90% of AOIs (default standR behaviour)
spe <- addPerROIQC(spe)
plotGeneQC(spe)
metadata(spe)$lcpm_threshold        # logCPM threshold used
dim(metadata(spe)$genes_rm_rawCount) # genes removed

# ---- 2. AOI (segment) QC -----------------------------------------------------
#   RawReads                          > 1000
#   Aligned / Stitched / Trimmed reads > 80% of RawReads
#   SequencingSaturation              > 50%
#   AOINucleiCount                    > 0  (lenient: small dataset)
#   Library size                      > 5e4
plotROIQC(spe, x_axis = "AOINucleiCount", x_lab = "Nuclei count",
          y_axis = "lib_size", y_lab = "Library size",
          col = tissue, y_threshold = 5e4)

cd <- colData(spe)
qc_raw        <- cd$RawReads > 1000
qc_processed  <- cd$AlignedReads  > cd$RawReads * 0.8 &
                 cd$StitchedReads > cd$RawReads * 0.8 &
                 cd$TrimmedReads  > cd$RawReads * 0.8
qc_saturation <- cd$SequencingSaturation > 50
qc_nuclei     <- cd$AOINucleiCount > 0
qc_libsize    <- cd$lib_size > 5e4

keep_aoi <- qc_raw & qc_processed & qc_saturation & qc_nuclei & qc_libsize
table(keep_aoi)
spe <- spe[, keep_aoi]

# ---- 3. Limit of quantification (LOQ) ----------------------------------------
#   LOQ_i = geomean(NegProbe_i) * geoSD(NegProbe_i)^2   (minimum 2)
#   Keep genes above LOQ in > 10% of AOIs and with < 5 counts in < 90% of AOIs
compute_loq <- function(neg_mat, n = 2, floor = 2) {
  apply(neg_mat, 2, function(x) {
    x <- x[!is.na(x) & x > 0]
    if (length(x) < 2) return(floor)
    max(exp(mean(log(x))) * exp(sd(log(x)))^n, floor)
  })
}

loq  <- compute_loq(as.matrix(metadata(spe)$NegProbes))
expr <- as.matrix(assay(spe, "counts"))

above_loq  <- sweep(expr, 2, loq[colnames(expr)], `>`)
keep_genes <- rowMeans(above_loq) > 0.10 & rowMeans(expr < 5) < 0.90
table(keep_genes)
spe <- spe[keep_genes, ]

# ---- 4. Dimensionality reduction ---------------------------------------------
design_vars <- c("subject", "biomaterial", "tissue", "stage")

spe <- scater::runPCA(spe)
spe <- scater::runUMAP(spe, dimred = "PCA")   # n_neighbors = 15 (default)

pca_plots  <- lapply(design_vars, function(v)
  drawPCA(spe, precomputed = reducedDim(spe, "PCA"), col = !!rlang::sym(v)) +
    labs(colour = v))
umap_plots <- lapply(design_vars, function(v)
  plotDR(spe, dimred = "UMAP", col = !!rlang::sym(v)) + labs(colour = v))

wrap_plots(pca_plots,  ncol = 2)
wrap_plots(umap_plots, ncol = 2)

# ---- 5. Relative log expression (RLE) ----------------------------------------
# Used to compare normalisations; TMM was retained for downstream analyses.
spe_tmm <- geomxNorm(spe, method = "TMM")

rle_raw    <- plotRLExpr(spe,     ordannots = "subject", assay = 1, color = subject) +
                ggtitle("Raw counts")
rle_logcpm <- plotRLExpr(spe,     ordannots = "subject", assay = 2, color = subject) +
                ggtitle("logCPM")
rle_tmm    <- plotRLExpr(spe_tmm, ordannots = "subject", assay = 2, color = subject) +
                ggtitle("TMM")
rle_raw / rle_logcpm / rle_tmm

# The filtered object `spe` (raw counts, QC-passed AOIs and genes) is the
# input of 2_DEA.R.
