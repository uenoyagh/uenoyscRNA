suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

VERSION <- "v6.7.4.2"
RES_COL <- "LSECclean_res0.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.4/objects/",
  "Mouse_MASH_LSEC_clean_RPCA_resolution_scan_v6.7.4.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC sample-specific contamination audit\n")
cat("Version:", VERSION, "\n")
cat("Diagnostic resolution: 0.3\n")
cat("====================================================\n\n")

obj <- readRDS(INPUT_RDS)

if (!RES_COL %in% colnames(obj@meta.data)) {
  stop("Missing metadata: ", RES_COL)
}

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

gene_sets <- list(

  LSEC = c(
    "Clec4g","Stab1","Stab2","Lyve1",
    "Fcgr2b","Mrc1","Oit3","Dnase1l3"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1","Apoa2","Fabp1",
    "Ass1","Cps1","Hnf4a","Aldob",
    "Tat","Mat1a","Pck1","Gcgr",
    "Cyp2c29","Cyp3a11"
  ),

  Myeloid = c(
    "Ptprc","Lyz2","Adgre1","Tyrobp"
  ),

  Mesenchymal = c(
    "Col1a1","Col3a1","Pdgfrb","Rgs5"
  ),

  Biliary = c(
    "Krt19","Epcam","Krt8","Krt18"
  ),

  Stress = c(
    "Fos","Fosb","Jun","Junb",
    "Atf3","Egr1","Hspa1a","Hspa1b"
  ),

  Cycling = c(
    "Mki67","Pcna","Mcm5","Cdt1",
    "Rad51","Brca1","Gins2","Dtl"
  )
)

gene_sets <- lapply(
  gene_sets,
  function(x) intersect(x, rownames(counts))
)

total_counts <- Matrix::colSums(counts)

for (nm in names(gene_sets)) {

  genes <- gene_sets[[nm]]

  if (length(genes) == 0) {
    obj[[paste0("pctCounts_", nm, "_v6742")]] <- 0
    next
  }

  x <- Matrix::colSums(
    counts[genes, , drop=FALSE]
  )

  obj[[paste0("pctCounts_", nm, "_v6742")]] <-
    100 * x / pmax(total_counts, 1)
}

obj$Hep_to_LSEC_log2_v6742 <- log2(
  (
    obj$pctCounts_Hepatocyte_v6742 + 0.01
  ) /
  (
    obj$pctCounts_LSEC_v6742 + 0.01
  )
)

obj$cluster_sample_v6742 <- paste0(
  "C",
  obj@meta.data[[RES_COL]],
  "__",
  obj$sample
)

cat("=== CLUSTER BY SAMPLE ===\n")
print(
  table(
    obj@meta.data[[RES_COL]],
    obj$sample
  )
)

# ---------------------------------------------------------
# Cluster/sample quantitative summary
# ---------------------------------------------------------

metrics <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "pctCounts_LSEC_v6742",
  "pctCounts_Hepatocyte_v6742",
  "Hep_to_LSEC_log2_v6742",
  "pctCounts_Myeloid_v6742",
  "pctCounts_Mesenchymal_v6742",
  "pctCounts_Biliary_v6742",
  "pctCounts_Stress_v6742",
  "pctCounts_Cycling_v6742"
)

metrics <- intersect(
  metrics,
  colnames(obj@meta.data)
)

groups <- unique(
  obj$cluster_sample_v6742
)

summary_list <- lapply(
  groups,
  function(g) {

    idx <- obj$cluster_sample_v6742 == g
    md <- obj@meta.data[idx, , drop=FALSE]

    out <- data.frame(
      group=g,
      cluster=as.character(md[[RES_COL]][1]),
      sample=as.character(md$sample[1]),
      condition=as.character(md$condition[1]),
      n_cells=nrow(md),
      stringsAsFactors=FALSE
    )

    for (m in metrics) {

      v <- md[[m]]

      out[[paste0(m, "_median")]] <-
        median(v, na.rm=TRUE)

      out[[paste0(m, "_Q25")]] <-
        quantile(v, 0.25, na.rm=TRUE)

      out[[paste0(m, "_Q75")]] <-
        quantile(v, 0.75, na.rm=TRUE)
    }

    out
  }
)

summary_df <- do.call(
  rbind,
  summary_list
)

write.csv(
  summary_df,
  file.path(
    TABDIR,
    "LSEC_res0.3_cluster_sample_contamination_summary_v6.7.4.2.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Hepatocyte fraction by cluster/sample
# ---------------------------------------------------------

plot_df <- obj@meta.data

plot_df$cluster <- factor(
  plot_df[[RES_COL]],
  levels=as.character(0:6)
)

p_hep <- ggplot(
  plot_df,
  aes(
    x=cluster,
    y=pctCounts_Hepatocyte_v6742
  )
) +
  geom_boxplot(
    outlier.size=0.15
  ) +
  facet_wrap(
    ~sample,
    scales="free_y",
    ncol=3
  ) +
  scale_y_continuous(
    trans="log1p"
  ) +
  theme_classic(base_size=11) +
  labs(
    title="Hepatocyte-transcript fraction in clean LSEC",
    x="res0.3 cluster",
    y="Hepatocyte gene counts (%)"
  )

ggsave(
  file.path(
    FIGDIR,
    "Hepatocyte_count_fraction_by_cluster_sample_v6.7.4.2.pdf"
  ),
  p_hep,
  width=11,
  height=8
)

# ---------------------------------------------------------
# LSEC vs hepatocyte count fraction
# ---------------------------------------------------------

p_scatter <- ggplot(
  plot_df,
  aes(
    x=pctCounts_LSEC_v6742,
    y=pctCounts_Hepatocyte_v6742,
    color=cluster
  )
) +
  geom_point(
    size=0.35,
    alpha=0.45
  ) +
  facet_wrap(
    ~sample,
    ncol=3
  ) +
  theme_classic(base_size=10) +
  labs(
    title="LSEC identity versus hepatocyte transcript fraction",
    x="LSEC gene counts (%)",
    y="Hepatocyte gene counts (%)",
    color="Cluster"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_vs_Hepatocyte_count_fraction_by_sample_v6.7.4.2.pdf"
  ),
  p_scatter,
  width=12,
  height=8
)

# ---------------------------------------------------------
# UMAP metadata FeaturePlots
# ---------------------------------------------------------

p_meta <- FeaturePlot(
  obj,
  features=c(
    "pctCounts_LSEC_v6742",
    "pctCounts_Hepatocyte_v6742",
    "Hep_to_LSEC_log2_v6742",
    "pctCounts_Stress_v6742"
  ),
  reduction="umapLSECclean",
  order=TRUE,
  ncol=2,
  min.cutoff="q05",
  max.cutoff="q95"
)

ggsave(
  file.path(
    FIGDIR,
    "LSEC_contamination_metadata_UMAP_v6.7.4.2.pdf"
  ),
  p_meta,
  width=10,
  height=8
)

# ---------------------------------------------------------
# Cluster × sample DotPlot
# ---------------------------------------------------------

audit_genes <- c(
  "Clec4g","Stab2","Fcgr2b","Mrc1",
  "Wnt2","Wnt9b","Rspo3",
  "Esm1","Icam1","Cxcl9","Cxcl10",

  "Alb","Ttr","Fabp1","Aldob",
  "Tat","Mat1a","Pck1","Gcgr",
  "Cyp2c29","Cyp3a11",

  "H2-Aa","H2-Ab1","Cd74",

  "Cnmd","Cd209b","Ctsj",
  "Bmper","Ccl19","Gjb2",

  "Mki67","Pcna","Rad51"
)

audit_genes <- audit_genes[
  audit_genes %in% rownames(obj)
]

Idents(obj) <- factor(
  obj$cluster_sample_v6742
)

p_dot <- DotPlot(
  obj,
  features=audit_genes,
  assay="RNA",
  dot.scale=5
) +
  RotatedAxis() +
  theme_classic(base_size=7) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  ) +
  ggtitle(
    "LSEC res0.3 cluster x sample audit"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_res0.3_cluster_sample_marker_DotPlot_v6.7.4.2.pdf"
  ),
  p_dot,
  width=18,
  height=14
)

# ---------------------------------------------------------
# Preserve diagnostic object
# ---------------------------------------------------------

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_LSEC_sample_specific_audit_v6.7.4.2.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.4.2.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.4.2 COMPLETE\n")
cat("Cells:", ncol(obj), "\n")
cat("No cells removed\n")
cat("Diagnostic resolution: 0.3\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
