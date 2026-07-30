#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(openxlsx)
  library(scales)
})

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS3_validation"
)

pdf_dir <- file.path(output_dir, "PDF")
csv_dir <- file.path(output_dir, "CSV")
xlsx_dir <- file.path(output_dir, "Excel")

dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(xlsx_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading RDS3...\n")
obj <- readRDS(rds_file)
meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)

first_existing <- function(candidates, columns, label) {
  hit <- candidates[candidates %in% columns]
  if (length(hit) == 0) {
    stop(paste0("No metadata column found for ", label, ": ",
                paste(candidates, collapse = ", ")))
  }
  hit[[1]]
}

sample_col <- first_existing(
  c("sample_for_R8plot_FIXED2", "sample_display_FIXED2",
    "sample_for_R8plot", "sample", "orig.ident"),
  colnames(meta), "sample"
)

condition_col <- first_existing(
  c("condition_FIXED2", "condition", "group"),
  colnames(meta), "condition"
)

cluster_col <- first_existing(
  c("cluster_for_R8plot_FIXED2",
    "integratedRPCA_snn_res.0.8", "seurat_clusters"),
  colnames(meta), "cluster"
)

celltype_col <- first_existing(
  c("celltype_for_R8plot_FIXED2", "celltype_for_R8plot",
    "celltype_auto_annotation", "celltype"),
  colnames(meta), "cell type"
)

reduction_name <- c("umapRPCA", "umap")[
  c("umapRPCA", "umap") %in% Reductions(obj)
][1]

if (is.na(reduction_name)) stop("No UMAP reduction found.")

meta <- meta %>%
  mutate(
    sample_validation = as.character(.data[[sample_col]]),
    condition_validation = as.character(.data[[condition_col]]),
    cluster_validation = as.character(.data[[cluster_col]]),
    celltype_validation = as.character(.data[[celltype_col]])
  )

cluster_levels <- unique(meta$cluster_validation)
cluster_numeric <- suppressWarnings(as.numeric(cluster_levels))
if (all(!is.na(cluster_numeric))) {
  cluster_levels <- as.character(sort(unique(cluster_numeric)))
} else {
  cluster_levels <- sort(unique(cluster_levels))
}
meta$cluster_validation <- factor(meta$cluster_validation, levels = cluster_levels)

preferred_sample_order <- c(
  "STD_rep1", "CDHFD_rep1", "Sham1", "Sham20", "Tx17", "Tx5"
)
sample_levels <- c(
  preferred_sample_order[preferred_sample_order %in% unique(meta$sample_validation)],
  setdiff(unique(meta$sample_validation), preferred_sample_order)
)
meta$sample_validation <- factor(meta$sample_validation, levels = sample_levels)

preferred_condition_order <- c("STD", "CDHFD", "Sham", "Tx")
condition_levels <- c(
  preferred_condition_order[
    preferred_condition_order %in% unique(meta$condition_validation)
  ],
  setdiff(unique(meta$condition_validation), preferred_condition_order)
)
meta$condition_validation <- factor(
  meta$condition_validation,
  levels = condition_levels
)

save_pdf <- function(p, filename, width, height) {
  ggsave(
    file.path(pdf_dir, filename),
    p,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
}

theme_validation <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

# 01 Object summary
object_summary <- tibble(
  item = c(
    "RDS_file", "Cells", "Features_DefaultAssay", "DefaultAssay",
    "Assays", "Reductions", "Sample_column", "Condition_column",
    "Cluster_column", "Celltype_column", "UMAP_reduction"
  ),
  value = c(
    rds_file, ncol(obj), nrow(obj), DefaultAssay(obj),
    paste(Assays(obj), collapse = "; "),
    paste(Reductions(obj), collapse = "; "),
    sample_col, condition_col, cluster_col, celltype_col, reduction_name
  )
)

sample_count <- meta %>%
  count(sample_validation, name = "cells")
condition_count <- meta %>%
  count(condition_validation, name = "cells")
cluster_count <- meta %>%
  count(cluster_validation, name = "cells")
celltype_count <- meta %>%
  count(celltype_validation, name = "cells") %>%
  arrange(desc(cells))

write_csv(object_summary, file.path(csv_dir, "01_ObjectSummary.csv"))
write_csv(sample_count, file.path(csv_dir, "01_SampleCounts.csv"))
write_csv(condition_count, file.path(csv_dir, "01_ConditionCounts.csv"))
write_csv(cluster_count, file.path(csv_dir, "01_ClusterCounts.csv"))
write_csv(celltype_count, file.path(csv_dir, "01_CelltypeCounts.csv"))

wb1 <- createWorkbook()
for (nm in c("ObjectSummary", "SampleCounts", "ConditionCounts",
             "ClusterCounts", "CelltypeCounts")) {
  addWorksheet(wb1, nm)
}
writeData(wb1, "ObjectSummary", object_summary)
writeData(wb1, "SampleCounts", sample_count)
writeData(wb1, "ConditionCounts", condition_count)
writeData(wb1, "ClusterCounts", cluster_count)
writeData(wb1, "CelltypeCounts", celltype_count)
saveWorkbook(
  wb1,
  file.path(xlsx_dir, "01_ObjectSummary.xlsx"),
  overwrite = TRUE
)

# 02 Cluster composition
cluster_by_sample <- meta %>%
  count(cluster_validation, sample_validation, name = "cells") %>%
  group_by(cluster_validation) %>%
  mutate(fraction_within_cluster = cells / sum(cells)) %>%
  ungroup()

cluster_by_condition <- meta %>%
  count(cluster_validation, condition_validation, name = "cells") %>%
  group_by(cluster_validation) %>%
  mutate(fraction_within_cluster = cells / sum(cells)) %>%
  ungroup()

write_csv(cluster_by_sample, file.path(csv_dir, "02_ClusterBySample.csv"))
write_csv(cluster_by_condition, file.path(csv_dir, "02_ClusterByCondition.csv"))

wb2 <- createWorkbook()
addWorksheet(wb2, "ClusterCounts")
addWorksheet(wb2, "ClusterBySample")
addWorksheet(wb2, "ClusterByCondition")
writeData(wb2, "ClusterCounts", cluster_count)
writeData(wb2, "ClusterBySample", cluster_by_sample)
writeData(wb2, "ClusterByCondition", cluster_by_condition)
saveWorkbook(
  wb2,
  file.path(xlsx_dir, "02_ClusterComposition.xlsx"),
  overwrite = TRUE
)

p_cluster_count <- ggplot(
  cluster_count,
  aes(cluster_validation, cells)
) +
  geom_col() +
  labs(
    title = "RDS3 cluster cell counts",
    x = "Cluster",
    y = "Number of cells"
  ) +
  theme_validation()

save_pdf(p_cluster_count, "02A_ClusterCellCounts.pdf", 12, 6)

p_cluster_sample <- ggplot(
  cluster_by_sample,
  aes(cluster_validation, fraction_within_cluster, fill = sample_validation)
) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Sample composition within each cluster",
    x = "Cluster",
    y = "Fraction",
    fill = "Sample"
  ) +
  theme_validation()

save_pdf(p_cluster_sample, "02B_ClusterComposition_BySample.pdf", 14, 7)

p_cluster_condition <- ggplot(
  cluster_by_condition,
  aes(cluster_validation, fraction_within_cluster, fill = condition_validation)
) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Condition composition within each cluster",
    x = "Cluster",
    y = "Fraction",
    fill = "Condition"
  ) +
  theme_validation()

save_pdf(
  p_cluster_condition,
  "02C_ClusterComposition_ByCondition.pdf",
  14,
  7
)

# 03 Cluster-cell type audit
cluster_celltype <- meta %>%
  count(cluster_validation, celltype_validation, name = "cells") %>%
  group_by(cluster_validation) %>%
  mutate(fraction_within_cluster = cells / sum(cells)) %>%
  ungroup()

cluster_purity <- cluster_celltype %>%
  group_by(cluster_validation) %>%
  slice_max(cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    cluster = cluster_validation,
    dominant_celltype = celltype_validation,
    dominant_cells = cells,
    purity = fraction_within_cluster
  ) %>%
  arrange(cluster)

write_csv(
  cluster_celltype,
  file.path(csv_dir, "03_ClusterCelltypeComposition.csv")
)
write_csv(cluster_purity, file.path(csv_dir, "03_ClusterPurity.csv"))

wb3 <- createWorkbook()
addWorksheet(wb3, "ClusterCelltype")
addWorksheet(wb3, "ClusterPurity")
writeData(wb3, "ClusterCelltype", cluster_celltype)
writeData(wb3, "ClusterPurity", cluster_purity)
saveWorkbook(
  wb3,
  file.path(xlsx_dir, "03_ClusterAnnotationAudit.xlsx"),
  overwrite = TRUE
)

p_cluster_celltype <- ggplot(
  cluster_celltype,
  aes(cluster_validation, fraction_within_cluster, fill = celltype_validation)
) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Cell-type composition within each cluster",
    x = "Cluster",
    y = "Fraction",
    fill = "Cell type"
  ) +
  theme_validation(10) +
  theme(legend.text = element_text(size = 8))

save_pdf(
  p_cluster_celltype,
  "03A_ClusterCelltypeComposition.pdf",
  16,
  8
)

p_cluster_purity <- ggplot(
  cluster_purity,
  aes(cluster, purity)
) +
  geom_col() +
  geom_hline(yintercept = c(0.8, 0.9),
             linetype = c("dashed", "dotted")) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(
    title = "Cluster annotation purity",
    subtitle = "Fraction of the dominant cell type in each cluster",
    x = "Cluster",
    y = "Purity"
  ) +
  theme_validation()

save_pdf(p_cluster_purity, "03B_ClusterPurity.pdf", 12, 6)

# 04 UMAP validation
umap_df <- as.data.frame(Embeddings(obj, reduction = reduction_name))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_barcode <- rownames(umap_df)

umap_df <- umap_df %>%
  left_join(
    meta %>%
      select(
        cell_barcode,
        sample_validation,
        condition_validation,
        cluster_validation,
        celltype_validation
      ),
    by = "cell_barcode"
  )

make_umap <- function(color_col, title, legend_title) {
  ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = .data[[color_col]]
    )
  ) +
    geom_point(size = 0.35, alpha = 0.8) +
    coord_equal() +
    labs(title = title, color = legend_title) +
    theme_void(base_size = 12) +
    theme(legend.position = "right")
}

save_pdf(
  make_umap("cluster_validation", "RDS3 UMAP by cluster", "Cluster"),
  "04A_UMAP_Cluster.pdf",
  11,
  9
)

save_pdf(
  make_umap("celltype_validation", "RDS3 UMAP by cell type", "Cell type"),
  "04B_UMAP_Celltype.pdf",
  12,
  9
)

save_pdf(
  make_umap("sample_validation", "RDS3 UMAP by sample", "Sample"),
  "04C_UMAP_Sample.pdf",
  11,
  9
)

save_pdf(
  make_umap("condition_validation", "RDS3 UMAP by condition", "Condition"),
  "04D_UMAP_Condition.pdf",
  11,
  9
)

p_split <- ggplot(
  umap_df,
  aes(UMAP_1, UMAP_2, color = celltype_validation)
) +
  geom_point(size = 0.25, alpha = 0.75) +
  facet_wrap(~sample_validation, ncol = 3) +
  coord_equal() +
  labs(
    title = "RDS3 cell types split by sample",
    color = "Cell type"
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 10)
  )

save_pdf(
  p_split,
  "04E_UMAP_Celltype_SplitBySample.pdf",
  15,
  10
)

# 05 QC summary
qc_columns <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.mt_for_filter"),
  colnames(meta)
)

if (length(qc_columns) > 0) {
  qc_cluster_summary <- meta %>%
    select(cluster_validation, all_of(qc_columns)) %>%
    pivot_longer(
      cols = all_of(qc_columns),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(cluster_validation, metric) %>%
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q05 = quantile(value, 0.05, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      q95 = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(
    qc_cluster_summary,
    file.path(csv_dir, "05_QC_ByCluster.csv")
  )

  wb5 <- createWorkbook()
  addWorksheet(wb5, "QC_ByCluster")
  writeData(wb5, "QC_ByCluster", qc_cluster_summary)
  saveWorkbook(
    wb5,
    file.path(xlsx_dir, "05_QC_ByCluster.xlsx"),
    overwrite = TRUE
  )
}

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

cat("\nRDS3 Phase 1 validation completed.\n")
cat("Output directory:\n", output_dir, "\n")
cat("\nNext phase: marker extraction, DotPlot and FeaturePlot.\n")
