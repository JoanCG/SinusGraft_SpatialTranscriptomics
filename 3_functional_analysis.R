# ==============================================================================
#  3. Functional analysis of the DEGs
# ------------------------------------------------------------------------------
#  Shows how the DEGs were annotated and tested for enrichment. The data are
#  not included, so the inputs below are left empty.
#
#  Parts
#    A) Over-representation analysis (ORA): GO BP, GO MF and KEGG for the
#       up- and down-regulated DEGs of a comparison (clusterProfiler)
#    B) Superpathway annotation: DEGs are mapped to GeneCards SuperPaths,
#       SuperPaths are grouped into curated biological processes
#       (osteogenesis / host response) and each process is summarised by the
#       direction and magnitude of change of its genes
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

FDR_CUTOFF <- 0.05

# ---- 0. Inputs ---------------------------------------------------------------
# results : list from 2_DEA.R; each element holds the edgeR table in $edgeR
results <- NULL

# ==============================================================================
#  A) Over-representation analysis
# ==============================================================================
run_ora <- function(genes) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)$ENTREZID
  if (length(ids) < 5) return(NULL)                  # too few genes for ORA

  params <- list(pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH",
                 minGSSize = 5, maxGSSize = 500)
  list(
    GO_BP = do.call(enrichGO, c(list(gene = ids, OrgDb = org.Hs.eg.db,
                                     ont = "BP", readable = TRUE), params)),
    GO_MF = do.call(enrichGO, c(list(gene = ids, OrgDb = org.Hs.eg.db,
                                     ont = "MF", readable = TRUE), params)),
    KEGG  = setReadable(do.call(enrichKEGG, c(list(gene = ids,
                                                   organism = "hsa"), params)),
                        OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  )
}

ora_dotplots <- function(ora, label, top_n = 10) {
  lapply(names(ora), function(db) {
    if (is.null(ora[[db]]) || nrow(as.data.frame(ora[[db]])) == 0)
      return(ggplot() + labs(title = paste(label, db, "- no enriched terms")) +
               theme_void())
    dotplot(ora[[db]], showCategory = top_n) + labs(title = paste(label, db))
  })
}

# Example: DEGs higher in each group of one comparison
res      <- results$stage_1_vs_0$edgeR
up_genes <- res$gene[res$FDR < FDR_CUTOFF & res$logFC > 0]
dn_genes <- res$gene[res$FDR < FDR_CUTOFF & res$logFC < 0]

ora_up <- run_ora(up_genes)
ora_dn <- run_ora(dn_genes)

# 3 x 2 panel: rows = GO BP / GO MF / KEGG, columns = direction
wrap_plots(c(rbind(ora_dotplots(ora_up, "Up"),
                   ora_dotplots(ora_dn, "Down"))), ncol = 2)

# ==============================================================================
#  B) Superpathway annotation and process-level summary
# ==============================================================================
# genecards   : GeneCards (GeneALaCart) batch query of the DEG symbols,
#               one row per gene x SuperPath  (InputTerm, Symbol, SuperPath)
# sp_classes  : manual classification of every SuperPath into
#                 general_process : "Osteogenesis", "Immune response", "Both"
#                 osteo_process   : specific osteogenic process (or NA)
#                 immune_process  : specific host-response process (or NA)
genecards  <- NULL
sp_classes <- NULL

annot <- genecards %>%
  inner_join(sp_classes, by = "SuperPath") %>%
  inner_join(res %>% select(gene, logFC, FDR) %>% filter(FDR < FDR_CUTOFF),
             by = c("InputTerm" = "gene")) %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down"))

# Per SuperPath: genes involved and net direction (+1 per up gene, -1 per down)
superpath_tab <- annot %>%
  group_by(SuperPath, general_process, osteo_process, immune_process) %>%
  summarise(genes     = paste(unique(Symbol), collapse = ", "),
            net_dir   = sum(ifelse(direction == "Up", 1, -1)),
            .groups   = "drop")

# Per specific process: share of SuperPaths leaning to each direction
summarise_process <- function(col) {
  superpath_tab %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(process = .data[[col]]) %>%
    summarise(n_superpaths = n(),
              pct_up       = 100 * mean(net_dir > 0),
              pct_down     = 100 * mean(net_dir < 0),
              .groups      = "drop")
}
summarise_process("osteo_process")
summarise_process("immune_process")

# Cumulative and mean logFC per process and direction
plot_process <- function(col, title) {
  d <- annot %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(process = .data[[col]], direction) %>%
    summarise(total_logFC = sum(logFC), mean_logFC = mean(logFC),
              .groups = "drop")
  p <- function(y, lab) ggplot(d, aes(process, .data[[y]], fill = direction)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_fill_manual(values = c(Up = "#e41a1c", Down = "#377eb8")) +
    labs(x = NULL, y = lab) + theme_minimal()
  (p("total_logFC", "Cumulative logFC") | p("mean_logFC", "Mean logFC")) +
    plot_annotation(title = title) + plot_layout(guides = "collect")
}
plot_process("osteo_process",  "Osteogenic processes")
plot_process("immune_process", "Host-response processes")
