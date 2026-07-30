#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

#============================================================
# 1. Paths
#============================================================

rds1_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_GSE325222_RH251217117_RPCA_integrated_symbolFixed.rds"
)

rds4_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_RH260519343_GSE325222_RPCA_integrated_celltype_annotated.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS4_diff_audit_20260730"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n", output_dir, "\n\n")

#============================================================
# 2. Load objects
#============================================================

cat("Loading RDS1...\n")
obj1 <- readRDS(rds1_file)

cat("Loading RDS4...\n")
obj4 <- readRDS(rds4_file)

meta1 <- obj1@meta.data
meta4 <- obj4@meta.data

meta1$cell_barcode <- rownames(meta1)
meta4$cell_barcode <- rownames(meta4)

cells1 <- Cells(obj1)
cells4 <- Cells(obj4)

#============================================================
# 3. Basic cell-set comparison
#============================================================

common_cells <- intersect(cells1, cells4)
lost_cells <- setdiff(cells1, cells4)
new_cells <- setdiff(cells4, cells1)

basic_summary <- tibble(
  metric = c(
    "RDS1_cells",
    "RDS4_cells",
    "common_cells",
    "lost_from_RDS1",
    "new_in_RDS4",
    "RDS4_fraction_of_RDS1"
  ),
  value = c(
    length(cells1),
    length(cells4),
    length(common_cells),
    length(lost_cells),
    length(new_cells),
    length(common_cells) / length(cells1)
  )
)

write_csv(
  basic_summary,
  file.path(output_dir, "01_basic_cellset_summary.csv")
)

cat("\nBasic comparison:\n")
print(basic_summary)

cat("\nRDS4 is an exact subset of RDS1? ",
    length(new_cells) == 0 &&
      length(common_cells) == length(cells4),
    "\n")

#============================================================
# 4. Export lost-cell metadata
#============================================================

lost_meta <- meta1 %>%
  filter(cell_barcode %in% lost_cells)

common_meta1 <- meta1 %>%
  filter(cell_barcode %in% common_cells)

write_csv(
  lost_meta,
  file.path(output_dir, "02_lost_cells_RDS1_metadata.csv")
)

write_csv(
  tibble(cell_barcode = lost_cells),
  file.path(output_dir, "03_lost_cell_barcodes.csv")
)

write_csv(
  tibble(cell_barcode = new_cells),
  file.path(output_dir, "04_new_cell_barcodes_in_RDS4.csv")
)

#============================================================
# 5. Sample-level audit
#============================================================

sample_summary_rds1 <- meta1 %>%
  count(sample, name = "RDS1_cells")

sample_summary_rds4 <- meta4 %>%
  count(sample, name = "RDS4_cells")

sample_summary_lost <- lost_meta %>%
  count(sample, name = "lost_cells")

sample_audit <- full_join(
  sample_summary_rds1,
  sample_summary_rds4,
  by = "sample"
) %>%
  full_join(
    sample_summary_lost,
    by = "sample"
  ) %>%
  mutate(
    across(
      c(RDS1_cells, RDS4_cells, lost_cells),
      ~replace_na(.x, 0L)
    ),
    retained_fraction = if_else(
      RDS1_cells > 0,
      RDS4_cells / RDS1_cells,
      NA_real_
    ),
    lost_fraction = if_else(
      RDS1_cells > 0,
      lost_cells / RDS1_cells,
      NA_real_
    )
  ) %>%
  arrange(sample)

write_csv(
  sample_audit,
  file.path(output_dir, "05_sample_level_audit.csv")
)

cat("\nSample-level audit:\n")
print(sample_audit)

#============================================================
# 6. Condition-level audit
#============================================================

condition_summary_rds1 <- meta1 %>%
  count(condition, name = "RDS1_cells")

condition_summary_rds4 <- meta4 %>%
  count(condition, name = "RDS4_cells")

condition_summary_lost <- lost_meta %>%
  count(condition, name = "lost_cells")

condition_audit <- full_join(
  condition_summary_rds1,
  condition_summary_rds4,
  by = "condition"
) %>%
  full_join(
    condition_summary_lost,
    by = "condition"
  ) %>%
  mutate(
    across(
      c(RDS1_cells, RDS4_cells, lost_cells),
      ~replace_na(.x, 0L)
    ),
    retained_fraction = if_else(
      RDS1_cells > 0,
      RDS4_cells / RDS1_cells,
      NA_real_
    ),
    lost_fraction = if_else(
      RDS1_cells > 0,
      lost_cells / RDS1_cells,
      NA_real_
    )
  ) %>%
  arrange(condition)

write_csv(
  condition_audit,
  file.path(output_dir, "06_condition_level_audit.csv")
)

cat("\nCondition-level audit:\n")
print(condition_audit)

#============================================================
# 7. Cluster-level audit using RDS1 clusters
#============================================================

cluster_col_candidates <- c(
  "integratedRPCA_snn_res.0.8",
  "seurat_clusters"
)

cluster_col1 <- cluster_col_candidates[
  cluster_col_candidates %in% colnames(meta1)
][1]

if (!is.na(cluster_col1)) {

  cluster_audit <- meta1 %>%
    mutate(
      cluster_RDS1 = as.character(.data[[cluster_col1]]),
      retained_in_RDS4 = cell_barcode %in% common_cells
    ) %>%
    count(
      cluster_RDS1,
      retained_in_RDS4,
      name = "n_cells"
    ) %>%
    pivot_wider(
      names_from = retained_in_RDS4,
      values_from = n_cells,
      values_fill = 0
    )

  if (!"TRUE" %in% colnames(cluster_audit)) {
    cluster_audit$`TRUE` <- 0L
  }

  if (!"FALSE" %in% colnames(cluster_audit)) {
    cluster_audit$`FALSE` <- 0L
  }

  cluster_audit <- cluster_audit %>%
    rename(
      retained_cells = `TRUE`,
      lost_cells = `FALSE`
    ) %>%
    mutate(
      total_cells = retained_cells + lost_cells,
      retained_fraction = retained_cells / total_cells,
      lost_fraction = lost_cells / total_cells
    ) %>%
    arrange(
      desc(lost_fraction),
      desc(lost_cells)
    )

  write_csv(
    cluster_audit,
    file.path(output_dir, "07_RDS1_cluster_loss_audit.csv")
  )

  cat("\nCluster-level loss audit:\n")
  print(cluster_audit, n = Inf)
}

#============================================================
# 8. Sample × cluster audit
#============================================================

if (!is.na(cluster_col1)) {

  sample_cluster_audit <- meta1 %>%
    mutate(
      cluster_RDS1 = as.character(.data[[cluster_col1]]),
      retained_in_RDS4 = cell_barcode %in% common_cells
    ) %>%
    count(
      sample,
      cluster_RDS1,
      retained_in_RDS4,
      name = "n_cells"
    ) %>%
    pivot_wider(
      names_from = retained_in_RDS4,
      values_from = n_cells,
      values_fill = 0
    )

  if (!"TRUE" %in% colnames(sample_cluster_audit)) {
    sample_cluster_audit$`TRUE` <- 0L
  }

  if (!"FALSE" %in% colnames(sample_cluster_audit)) {
    sample_cluster_audit$`FALSE` <- 0L
  }

  sample_cluster_audit <- sample_cluster_audit %>%
    rename(
      retained_cells = `TRUE`,
      lost_cells = `FALSE`
    ) %>%
    mutate(
      total_cells = retained_cells + lost_cells,
      lost_fraction = if_else(
        total_cells > 0,
        lost_cells / total_cells,
        NA_real_
      )
    ) %>%
    arrange(
      sample,
      desc(lost_fraction),
      desc(lost_cells)
    )

  write_csv(
    sample_cluster_audit,
    file.path(output_dir, "08_sample_by_cluster_loss_audit.csv")
  )
}

#============================================================
# 9. QC comparison: retained vs lost
#============================================================

qc_columns <- intersect(
  c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.mt_for_filter"
  ),
  colnames(meta1)
)

qc_audit_long <- meta1 %>%
  mutate(
    status = if_else(
      cell_barcode %in% common_cells,
      "Retained_in_RDS4",
      "Lost_from_RDS4"
    )
  ) %>%
  select(
    cell_barcode,
    sample,
    condition,
    status,
    all_of(qc_columns)
  ) %>%
  pivot_longer(
    cols = all_of(qc_columns),
    names_to = "QC_metric",
    values_to = "value"
  )

qc_summary <- qc_audit_long %>%
  group_by(
    sample,
    condition,
    status,
    QC_metric
  ) %>%
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
  qc_summary,
  file.path(output_dir, "09_QC_retained_vs_lost_summary.csv")
)

#============================================================
# 10. QC violin plots
#============================================================

if (length(qc_columns) > 0) {

  for (metric_name in qc_columns) {

    plot_data <- meta1 %>%
      mutate(
        status = if_else(
          cell_barcode %in% common_cells,
          "Retained_in_RDS4",
          "Lost_from_RDS4"
        )
      ) %>%
      filter(!is.na(.data[[metric_name]]))

    p <- ggplot(
      plot_data,
      aes(
        x = status,
        y = .data[[metric_name]]
      )
    ) +
      geom_violin(
        trim = TRUE,
        scale = "width"
      ) +
      geom_boxplot(
        width = 0.15,
        outlier.shape = NA
      ) +
      facet_wrap(
        ~sample,
        scales = "free_y",
        ncol = 3
      ) +
      labs(
        title = paste0(
          metric_name,
          ": retained versus lost cells"
        ),
        x = NULL,
        y = metric_name
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(
          angle = 25,
          hjust = 1
        )
      )

    ggsave(
      filename = file.path(
        output_dir,
        paste0(
          "10_QC_",
          metric_name,
          "_retained_vs_lost.pdf"
        )
      ),
      plot = p,
      width = 11,
      height = 8
    )
  }
}

#============================================================
# 11. UMAP audit using RDS1 coordinates
#============================================================

umap_candidates <- c(
  "umapRPCA",
  "umap"
)

umap_name1 <- umap_candidates[
  umap_candidates %in% Reductions(obj1)
][1]

if (!is.na(umap_name1)) {

  umap_df <- as.data.frame(
    Embeddings(obj1, reduction = umap_name1)
  )

  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$cell_barcode <- rownames(umap_df)

  umap_df <- umap_df %>%
    left_join(
      meta1 %>%
        select(
          cell_barcode,
          sample,
          condition,
          all_of(cluster_col1)
        ),
      by = "cell_barcode"
    ) %>%
    mutate(
      status = if_else(
        cell_barcode %in% common_cells,
        "Retained_in_RDS4",
        "Lost_from_RDS4"
      )
    )

  # Overall UMAP
  p_umap_status <- ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = status
    )
  ) +
    geom_point(
      size = 0.15,
      alpha = 0.65
    ) +
    coord_equal() +
    labs(
      title = "RDS1 UMAP: cells retained or lost in RDS4",
      color = NULL
    ) +
    theme_void(base_size = 13) +
    theme(
      legend.position = "right"
    )

  ggsave(
    file.path(
      output_dir,
      "11_UMAP_RDS1_retained_vs_lost.pdf"
    ),
    p_umap_status,
    width = 10,
    height = 8
  )

  # Sample-split UMAP
  p_umap_sample <- ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = status
    )
  ) +
    geom_point(
      size = 0.12,
      alpha = 0.65
    ) +
    facet_wrap(
      ~sample,
      ncol = 3
    ) +
    coord_equal() +
    labs(
      title = "RDS1 UMAP by sample: retained or lost in RDS4",
      color = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 11)
    )

  ggsave(
    file.path(
      output_dir,
      "12_UMAP_RDS1_retained_vs_lost_by_sample.pdf"
    ),
    p_umap_sample,
    width = 13,
    height = 9
  )
}

#============================================================
# 12. Lost-cell fraction plots
#============================================================

p_sample_loss <- ggplot(
  sample_audit,
  aes(
    x = reorder(sample, lost_fraction),
    y = lost_fraction
  )
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(
    title = "Fraction of RDS1 cells absent from RDS4",
    x = NULL,
    y = "Lost fraction"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(
    output_dir,
    "13_sample_lost_fraction.pdf"
  ),
  p_sample_loss,
  width = 8,
  height = 5
)

if (exists("cluster_audit")) {

  p_cluster_loss <- ggplot(
    cluster_audit,
    aes(
      x = reorder(cluster_RDS1, lost_fraction),
      y = lost_fraction
    )
  ) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1)
    ) +
    labs(
      title = "Lost fraction by RDS1 cluster",
      x = "RDS1 cluster",
      y = "Lost fraction"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    file.path(
      output_dir,
      "14_RDS1_cluster_lost_fraction.pdf"
    ),
    p_cluster_loss,
    width = 8,
    height = 10
  )
}

#============================================================
# 13. Cell order and embedding identity checks
#============================================================

cell_identity_check <- tibble(
  check = c(
    "RDS4_is_subset_of_RDS1",
    "RDS4_cell_order_matches_RDS1_subset_order",
    "Common_cell_count_equals_RDS4_count"
  ),
  result = c(
    all(cells4 %in% cells1),
    identical(
      cells4,
      cells1[cells1 %in% cells4]
    ),
    length(common_cells) == length(cells4)
  )
)

write_csv(
  cell_identity_check,
  file.path(output_dir, "15_cell_identity_checks.csv")
)

# Compare common-cell PCA coordinates where possible
if (
  "pca" %in% Reductions(obj1) &&
  "pca" %in% Reductions(obj4)
) {

  pca1 <- Embeddings(obj1, "pca")
  pca4 <- Embeddings(obj4, "pca")

  common_pca_cells <- intersect(
    rownames(pca1),
    rownames(pca4)
  )

  n_pc_compare <- min(
    10,
    ncol(pca1),
    ncol(pca4)
  )

  pca_diff <- pca1[
    common_pca_cells,
    seq_len(n_pc_compare),
    drop = FALSE
  ] -
    pca4[
      common_pca_cells,
      seq_len(n_pc_compare),
      drop = FALSE
    ]

  pca_identity <- tibble(
    metric = c(
      "common_cells_compared",
      "PCs_compared",
      "maximum_absolute_difference",
      "mean_absolute_difference",
      "PCA_coordinates_identical"
    ),
    value = c(
      length(common_pca_cells),
      n_pc_compare,
      max(abs(pca_diff), na.rm = TRUE),
      mean(abs(pca_diff), na.rm = TRUE),
      identical(
        pca1[
          common_pca_cells,
          seq_len(n_pc_compare),
          drop = FALSE
        ],
        pca4[
          common_pca_cells,
          seq_len(n_pc_compare),
          drop = FALSE
        ]
      )
    )
  )

  write_csv(
    pca_identity,
    file.path(output_dir, "16_PCA_coordinate_identity.csv")
  )
}

#============================================================
# 14. Session information
#============================================================

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "17_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("RDS4 differential audit completed\n")
cat("============================================================\n")
cat("Output directory:\n")
cat(output_dir, "\n\n")

cat("Most important output files:\n")
cat("01_basic_cellset_summary.csv\n")
cat("05_sample_level_audit.csv\n")
cat("07_RDS1_cluster_loss_audit.csv\n")
cat("08_sample_by_cluster_loss_audit.csv\n")
cat("09_QC_retained_vs_lost_summary.csv\n")
cat("11_UMAP_RDS1_retained_vs_lost.pdf\n")
cat("12_UMAP_RDS1_retained_vs_lost_by_sample.pdf\n")
cat("15_cell_identity_checks.csv\n")
cat("16_PCA_coordinate_identity.csv\n")
