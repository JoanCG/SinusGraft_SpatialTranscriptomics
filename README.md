# SinusGraft spatial transcriptomics

R code for the spatial transcriptomics analyses (NanoString GeoMx DSP, Whole
Transcriptome Atlas) in our maxillary sinus grafting study.

The scripts are meant to show how the analyses were done: which filters,
models and thresholds were used. They are not a turnkey pipeline. The data are
not included, so each script starts with an empty inputs section describing
what it expects.

## Scripts

Run them in order; each one uses the objects created by the previous one.

1. `1_QC_EDA.R`: builds the SpatialExperiment object, runs gene, AOI and LOQ
   filtering, and makes the PCA, UMAP and RLE plots.
2. `2_DEA.R`: differential expression with edgeR (paired design, likelihood
   ratio test), checked against limma-voom with `duplicateCorrelation`. Also
   includes diagnostic plots, a post-hoc power estimate, the DEG heatmap and a
   co-expression network among DEGs built from linear mixed models.
3. `3_functional_analysis.R`: GO and KEGG over-representation analysis
   (clusterProfiler), plus annotation of the DEGs with GeneCards SuperPaths
   grouped into osteogenic and host-response processes.

## Design variables

The code uses four variables per AOI:

- `subject`: patient
- `biomaterial`: graft material
- `tissue`: mineralized or connective compartment
- `stage`: graft maturation stage

## Packages

Analyses were run in R 4.5 (Bioconductor 3.21) with edgeR 4.4.2. Other
edgeR versions can give slightly different DEG counts.

- Bioconductor: SpatialExperiment, standR, scater, edgeR, limma, qvalue,
  ComplexHeatmap, clusterProfiler, org.Hs.eg.db
- CRAN: ggplot2, ggrepel, patchwork, dplyr, tidyr, seriation, ggstatsplot,
  lme4, lmerTest, igraph, ggraph

## Contact

Joan Calle (jcalle@go.ugr.es)
