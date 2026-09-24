# ==============================================================================
#  2. Differential expression (edgeR), validated with limma-voom
# ------------------------------------------------------------------------------
#  Shows the models and settings used in the paper. The data are not
#  included, so the inputs below are left empty.
#
#  Every comparison follows the same recipe (run_contrast()):
#    a) Keep only "complete" blocks: AOIs from the same subject and location
#       that contain every level being compared (paired design)
#    b) edgeR: TMM normalisation, robust dispersion, GLM with the block as a
#       fixed effect, likelihood-ratio test
#    c) limma-voom with duplicateCorrelation (block as random effect) as an
#       independent validation of the edgeR results
#    d) Post-hoc power estimate (pi1 = 1 - pi0, qvalue) and diagnostics
#       (BCV, p-value / logFC distributions, volcano plot)
#
#  Comparisons
#    1) Biomaterial                        (all AOIs)
#    2) Biomaterial within one tissue      (mineralized compartment)
#    3) Stage contrasts within one tissue  (connective compartment)
#    4) Tissue compartment                 (all AOIs)
#
#  Plus: heatmap of the DEGs from comparison 1, a paired plot for a gene of
#  interest and a co-expression network among DEGs (pairwise linear mixed
#  models with subject as a random intercept).
# ==============================================================================

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(edgeR)
  library(limma)
  library(qvalue)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(seriation)
  library(ggstatsplot)
  library(lme4)
  library(lmerTest)
  library(igraph)
  library(ggraph)
})

FDR_CUTOFF <- 0.05

# ---- 0. Inputs ---------------------------------------------------------------
# spe : filtered SpatialExperiment from 1_QC_EDA.R (raw counts). colData must
#       contain subject, biomaterial, tissue and stage (factors).
spe <- NULL

counts <- assay(spe, "counts")
meta   <- as.data.frame(colData(spe))

# ---- 1. Generic paired contrast ----------------------------------------------
#' @param subset      logical vector selecting the AOIs to use (NULL = all)
#' @param test_var    column of `meta` with the levels being compared
#' @param contrast    contrast between levels of `test_var`, e.g. "B - A"
#' @param block_vars  columns that define a block (paired unit); only blocks
#'                    containing every level of `test_var` are kept
run_contrast <- function(name, test_var, contrast, block_vars, subset = NULL) {

  m <- if (is.null(subset)) meta else meta[subset, , drop = FALSE]
  m <- droplevels(m)

  # a) Paired design: keep complete blocks only
  m$block  <- interaction(m[, block_vars], drop = TRUE, sep = "_")
  n_levels <- nlevels(m[[test_var]])
  complete <- names(which(rowSums(table(m$block, m[[test_var]]) > 0) == n_levels))
  m        <- droplevels(m[m$block %in% complete, , drop = FALSE])
  y        <- counts[, rownames(m)]
  m$test   <- m[[test_var]]
  message(sprintf("[%s] %d AOIs in %d blocks", name, ncol(y), length(complete)))

  # b) edgeR: GLM with block as fixed effect, LRT on the contrast
  design <- model.matrix(~ 0 + test + block, data = m)
  colnames(design) <- make.names(sub("^test", "", colnames(design)))
  con    <- makeContrasts(contrasts = contrast, levels = design)

  dge <- DGEList(y, samples = m, group = m$test)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- glmLRT(fit, contrast = con)

  res <- topTags(lrt, n = Inf)$table
  res$gene      <- rownames(res)
  res$direction <- case_when(res$FDR < FDR_CUTOFF & res$logFC > 0 ~ "Up",
                             res$FDR < FDR_CUTOFF & res$logFC < 0 ~ "Down",
                             TRUE                                 ~ "ns")

  # c) limma-voom validation: block as a random effect
  design_l <- model.matrix(~ 0 + test, data = m)
  colnames(design_l) <- make.names(sub("^test", "", colnames(design_l)))
  v      <- voom(dge, design_l)
  corfit <- duplicateCorrelation(v, design_l, block = m$block)
  fit_l  <- lmFit(v, design_l, block = m$block,
                  correlation = corfit$consensus.correlation)
  fit_l  <- eBayes(contrasts.fit(fit_l, makeContrasts(contrasts = contrast,
                                                      levels = design_l)))
  res_l  <- topTable(fit_l, number = Inf)

  shared <- intersect(res$gene[res$FDR < FDR_CUTOFF],
                      rownames(res_l)[res_l$adj.P.Val < FDR_CUTOFF])
  message(sprintf("[%s] DEGs edgeR = %d | limma = %d | shared = %d | logFC r = %.2f",
                  name, sum(res$FDR < FDR_CUTOFF),
                  sum(res_l$adj.P.Val < FDR_CUTOFF), length(shared),
                  cor(res$logFC, res_l[res$gene, "logFC"])))

  # d) Post-hoc power: estimated proportion of truly DE genes
  pi1 <- 1 - qvalue(res$PValue)$pi0
  message(sprintf("[%s] estimated pi1 = %.3f", name, pi1))

  list(name = name, edgeR = res, limma = res_l, dge = dge, meta = m, pi1 = pi1)
}

# ---- 2. Diagnostics and volcano ----------------------------------------------
plot_diagnostics <- function(r) {
  res <- r$edgeR
  p_fdr <- ggplot(res, aes(FDR)) +
    geom_histogram(bins = 40, fill = "skyblue", colour = "black") +
    labs(x = "FDR", y = "Genes")
  p_lfc <- ggplot(res, aes(logFC)) +
    geom_histogram(bins = 40, fill = "orange", colour = "black") +
    labs(x = "log2 fold change", y = "Genes")
  p_cpm <- ggplot(res, aes(logCPM, FDR, colour = FDR < FDR_CUTOFF)) +
    geom_point(alpha = 0.7) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"),
                        name = sprintf("FDR < %.2f", FDR_CUTOFF))
  (p_fdr | p_lfc) / p_cpm + plot_annotation(title = r$name) &
    theme_minimal(base_size = 12)
}

plot_volcano <- function(r) {
  ggplot(r$edgeR, aes(logFC, -log10(FDR), colour = direction)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
    geom_text_repel(data = subset(r$edgeR, direction != "ns"),
                    aes(label = gene), size = 3, max.overlaps = Inf,
                    show.legend = FALSE) +
    scale_colour_manual(values = c(Up = "red", Down = "blue", ns = "grey60")) +
    labs(title = r$name, x = "log2 fold change", y = "-log10(FDR)") +
    theme_minimal()
}

# ---- 3. Comparisons ----------------------------------------------------------
# Use the factor levels of your own colData here.
results <- list(

  # 1) Biomaterial: blocks: same subject, tissue and stage
  biomaterial = run_contrast(
    "Biomaterial", test_var = "biomaterial", contrast = "MaterialB - MaterialA",
    block_vars = c("subject", "tissue", "stage")),

  # 2) Biomaterial within the mineralized compartment
  biomaterial_mineralized = run_contrast(
    "Biomaterial | mineralized", test_var = "biomaterial",
    contrast = "MaterialB - MaterialA", block_vars = c("subject", "stage"),
    subset = meta$tissue == "mineralized"),

  # 3) Stage contrasts within the connective compartment
  #    blocks: same subject and biomaterial, all stages present
  stage_2_vs_0 = run_contrast(
    "Stage 2 vs 0 | connective", test_var = "stage", contrast = "S2 - S0",
    block_vars = c("subject", "biomaterial"),
    subset = meta$tissue == "connective"),
  stage_1_vs_0 = run_contrast(
    "Stage 1 vs 0 | connective", test_var = "stage", contrast = "S1 - S0",
    block_vars = c("subject", "biomaterial"),
    subset = meta$tissue == "connective"),

  # 4) Tissue compartment: blocks: same subject, biomaterial and stage
  tissue = run_contrast(
    "Tissue", test_var = "tissue", contrast = "connective - mineralized",
    block_vars = c("subject", "biomaterial", "stage"))
)

lapply(results, plot_diagnostics)
lapply(results, plot_volcano)
lapply(results, function(r) plotBCV(r$dge, main = r$name))

summary_tab <- do.call(rbind, lapply(results, function(r) data.frame(
  comparison = r$name,
  n_AOIs     = ncol(r$dge),
  edgeR_DEGs = sum(r$edgeR$FDR < FDR_CUTOFF),
  limma_DEGs = sum(r$limma$adj.P.Val < FDR_CUTOFF),
  pi1        = round(r$pi1, 3))))
summary_tab

# ---- 4. Heatmap of biomaterial DEGs ------------------------------------------
r    <- results$biomaterial
degs <- r$edgeR$gene[r$edgeR$FDR < FDR_CUTOFF]

log_cpm <- cpm(calcNormFactors(DGEList(r$dge$counts), method = "TMM"), log = TRUE)
z       <- t(scale(t(log_cpm[degs, ])))              # row z-scores

anno <- HeatmapAnnotation(df = r$meta[, c("biomaterial", "tissue", "stage")])
o_rows <- seriate(dist(z),    method = "GW")          # seriation-based ordering
o_cols <- seriate(dist(t(z)), method = "GW")

draw(Heatmap(z, name = "Z-score", top_annotation = anno,
             cluster_rows    = as.dendrogram(o_rows[[1]]),
             cluster_columns = as.dendrogram(o_cols[[1]]),
             row_split = 2, column_split = 2,
             show_column_names = FALSE,
             row_names_gp = grid::gpar(fontsize = 7)))

# ---- 5. Paired plot for a gene of interest -----------------------------------
# Expression of one DEG across biomaterials within matched blocks; the
# subtitle reports the edgeR statistics rather than the plot's own test.
GENE_OF_INTEREST <- "GENE"
r <- results$biomaterial_mineralized

tmm_cpm <- cpm(calcNormFactors(DGEList(r$dge$counts), method = "TMM"))
gene_df <- r$meta %>%
  mutate(log_expr = log2(tmm_cpm[GENE_OF_INTEREST, rownames(r$meta)] + 1)) %>%
  arrange(block, biomaterial)                          # pairs in matching order
stats   <- r$edgeR[GENE_OF_INTEREST, ]

ggwithinstats(
  data = gene_df, x = biomaterial, y = log_expr,          # paired by row order
  type = "non-parametric", results.subtitle = FALSE, bf.message = FALSE,
  xlab = "Biomaterial", ylab = "log2(TMM CPM + 1)",
  title = GENE_OF_INTEREST,
  subtitle = sprintf("edgeR logFC = %.2f, FDR = %.3g", stats$logFC, stats$FDR))

# ---- 6. Co-expression network among DEGs ------------------------------------
# Nodes: DEGs of one comparison changing in the same direction.
# For every gene pair (A, B), both orderings are fitted:
#   A ~ B + biomaterial + tissue + stage + (1 | subject)
# If the random effect is singular, subject is entered as a fixed effect.
# A pair is an edge if BH q < 0.05 (on the larger of the two p-values), the
# two slopes share their sign and the mean |beta| >= 0.2.
res_net   <- results$biomaterial$edgeR
net_genes <- res_net$gene[res_net$FDR < FDR_CUTOFF & res_net$logFC < 0]

log_cpm_all <- cpm(calcNormFactors(DGEList(assay(spe, "counts")), method = "TMM"),
               log = TRUE)
df <- as.data.frame(colData(spe))[, c("subject", "biomaterial", "tissue", "stage")]

fit_pair <- function(y_gene, x_gene) {
  df$Y <- log_cpm_all[y_gene, rownames(df)]
  df$X <- log_cpm_all[x_gene, rownames(df)]
  fit <- suppressMessages(lmerTest::lmer(
    Y ~ X + biomaterial + tissue + stage + (1 | subject), data = df,
    control = lmerControl(optimizer = "bobyqa")))
  if (isSingular(fit))
    fit <- lm(Y ~ X + biomaterial + tissue + stage + subject, data = df)
  co <- coef(summary(fit))["X", ]
  data.frame(Y = y_gene, X = x_gene, beta = co[["Estimate"]],
             p = co[[length(co)]])                   # last column = p-value
}

pairs <- combn(net_genes, 2, simplify = FALSE)
fits  <- do.call(rbind, lapply(pairs, function(g)
  rbind(fit_pair(g[1], g[2]), fit_pair(g[2], g[1]))))

edges <- fits %>%
  mutate(from = pmin(Y, X), to = pmax(Y, X)) %>%
  group_by(from, to) %>%
  summarise(sign_consistent = n_distinct(sign(beta)) == 1,
            beta            = mean(beta),
            p               = max(p),
            .groups = "drop") %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  filter(q < 0.05, sign_consistent, abs(beta) >= 0.2)

g <- graph_from_data_frame(edges, vertices = data.frame(name = net_genes),
                           directed = FALSE)
V(g)$degree <- degree(g)

set.seed(7)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_width = abs(beta)), colour = "#d95f02", alpha = 0.7) +
  geom_node_point(aes(size = degree), shape = 21, fill = "grey80") +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold") +
  scale_edge_width_continuous(range = c(0.4, 2.5), name = "|beta|") +
  scale_size_continuous(range = c(3, 12), name = "Degree") +
  theme_graph(base_family = "sans")

# `results` (edgeR tables per comparison) is the input of
# 3_functional_analysis.R.
