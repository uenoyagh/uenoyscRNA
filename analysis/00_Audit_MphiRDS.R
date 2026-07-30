#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(openxlsx)
  library(patchwork)
})

# ============================================================
# 00_Audit_MphiRDS.R
# Audit a macrophage-only Seurat object before annotation review
# ============================================================

set.seed(260730)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/results/",
  "Mouse_Mphi_manual_annotation_v2_cluster1_40_fraction_WITH_CLUSTER_LABELS/",
  "Mouse_Mphi_manual_annotation_v2_cluster_fraction_WITH_CLUSTER_LABELS.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/results/",
  "Mouse_Mphi_AnnotationValidation/",
  "00_Audit"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_dir <- file.path(output_dir, "PDF")
csv_dir <- file.path(output_dir, "CSV")
log_dir <- file.path(output_dir, "Logs")

for (x in c(pdf_dir, csv_dir, log_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# Candidate metadata names
# ------------------------------------------------------------

cluster_candidates <- c(
  "cluster_label",
  "cluster_for_plot",
  "cluster_for_R8plot",
  "cluster_for_R8plot_FIXED2",
  "seurat_clusters"
)

annotation_candidates <- c(
  "manual_annotation_v2",
  "manual_annotation",
  "celltype_manual",
  "celltype_annotation",
  "celltype_for_R8plot_FIXED2",
  "celltype_for_R8plot",
  "celltype_auto_annotation"
)

sample_candidates <- c(
  "sample_for_R8plot_FIXED2",
  "sample_for_R8plot",
  "sample",
  "orig.ident"
)

condition_candidates <- c(
  "condition_FIXED2",
  "condition",
  "group",
  "treatment"
)

reduction_candidates <- c(
  "umap",
  "integratedRPCA_umap",
  "rpca_umap"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

first_existing <- function(candidates, available) {
  hit <- candidates[candidates %in% available]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

save_pdf <- function(plot, file, width = 11, height = 8.5) {
  ggsave(
    filename = file,
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    limitsize = FALSE
  )
}

write_excel <- function(file, sheets) {
  wb <- createWorkbook()

  for (nm in names(sheets)) {
    sheet_name <- substr(
      gsub("[\\\\/:*?\\[\\]]", "_", nm),
      1,
      31
    )

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, sheets[[nm]], withFilter = TRUE)
    freezePane(wb, sheet_name, firstRow = TRUE)

    if (ncol(sheets[[nm]]) > 0) {
      setColWidths(
        wb,
        sheet = sheet_name,
        cols = seq_len(ncol(sheets[[nm]])),
        widths = "auto"
      )
    }
  }

  saveWorkbook(wb, file, overwrite = TRUE)
}

count_metadata <- function(object, column_name, output_name) {
  if (is.na(column_name)) {
    return(
      tibble(
        variable = output_name,
        value = "NOT_DETECTED",
        cells = NA_integer_
      )
    )
  }

  object@meta.data %>%
    count(
      value = .data[[column_name]],
      name = "cells"
    ) %>%
    mutate(
      variable = output_name,
      value = as.character(value)
    ) %>%
    select(variable, value, cells) %>%
    arrange(desc(cells))
}

make_dimplot <- function(
  object,
  reduction_name,
  group_column,
  title_text,
  label = TRUE
) {
  if (is.na(group_column)) {
    return(NULL)
  }

  DimPlot(
    object = object,
    reduction = reduction_name,
    group.by = group_column,
    label = label,
    repel = TRUE,
    pt.size = 0.35,
    raster = TRUE
  ) +
    ggtitle(title_text) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.text = element_text(size = 8)
    )
}

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

cat("============================================================\n")
cat("Macrophage RDS audit\n")
cat("============================================================\n")
cat("RDS file:\n", rds_file, "\n\n")

if (!file.exists(rds_file)) {
  stop("RDS file not found: ", rds_file)
}

obj <- readRDS(rds_file)

metadata_columns <- colnames(obj@meta.data)
reduction_names <- Reductions(obj)

cluster_col <- first_existing(
  cluster_candidates,
  metadata_columns
)

annotation_col <- first_existing(
  annotation_candidates,
  metadata_columns
)

sample_col <- first_existing(
  sample_candidates,
  metadata_columns
)

condition_col <- first_existing(
  condition_candidates,
  metadata_columns
)

reduction_name <- first_existing(
  reduction_candidates,
  reduction_names
)

if (is.na(reduction_name)) {
  stop(
    "No supported UMAP reduction detected. Available reductions: ",
    paste(reduction_names, collapse = ", ")
  )
}

cat("Detected fields\n")
cat("  Cluster:   ", cluster_col, "\n")
cat("  Annotation:", annotation_col, "\n")
cat("  Sample:    ", sample_col, "\n")
cat("  Condition: ", condition_col, "\n")
cat("  UMAP:      ", reduction_name, "\n\n")

# ------------------------------------------------------------
# Object summary
# ------------------------------------------------------------

assay_summary <- tibble(
  assay = Assays(obj),
  is_default = Assays(obj) == DefaultAssay(obj),
  features = vapply(
    Assays(obj),
    function(a) nrow(obj[[a]]),
    integer(1)
  )
)

reduction_summary <- tibble(
  reduction = reduction_names,
  dimensions = vapply(
    reduction_names,
    function(r) ncol(Embeddings(obj, reduction = r)),
    integer(1)
  )
)

metadata_summary <- tibble(
  metadata_column = metadata_columns,
  class = vapply(
    obj@meta.data,
    function(x) paste(class(x), collapse = ";"),
    character(1)
  ),
  non_missing = vapply(
    obj@meta.data,
    function(x) sum(!is.na(x)),
    integer(1)
  ),
  unique_values = vapply(
    obj@meta.data,
    function(x) dplyr::n_distinct(x, na.rm = TRUE),
    integer(1)
  )
)

object_summary <- tibble(
  item = c(
    "RDS_file",
    "Cells",
    "Features_default_assay",
    "Default_assay",
    "Assays",
    "Reductions",
    "Metadata_columns",
    "Detected_cluster_column",
    "Detected_annotation_column",
    "Detected_sample_column",
    "Detected_condition_column",
    "Detected_UMAP_reduction"
  ),
  value = c(
    rds_file,
    ncol(obj),
    nrow(obj),
    DefaultAssay(obj),
    paste(Assays(obj), collapse = "; "),
    paste(reduction_names, collapse = "; "),
    ncol(obj@meta.data),
    cluster_col,
    annotation_col,
    sample_col,
    condition_col,
    reduction_name
  )
)

# ------------------------------------------------------------
# Frequency tables
# ------------------------------------------------------------

cluster_counts <- count_metadata(
  obj,
  cluster_col,
  "cluster"
)

annotation_counts <- count_metadata(
  obj,
  annotation_col,
  "annotation"
)

sample_counts <- count_metadata(
  obj,
  sample_col,
  "sample"
)

condition_counts <- count_metadata(
  obj,
  condition_col,
  "condition"
)

all_counts <- bind_rows(
  cluster_counts,
  annotation_counts,
  sample_counts,
  condition_counts
)

# Cross-tabulation
cluster_annotation_table <- tibble()

if (!is.na(cluster_col) && !is.na(annotation_col)) {
  cluster_annotation_table <- obj@meta.data %>%
    count(
      cluster = .data[[cluster_col]],
      annotation = .data[[annotation_col]],
      name = "cells"
    ) %>%
    group_by(cluster) %>%
    mutate(
      cluster_total = sum(cells),
      fraction_within_cluster = cells / cluster_total
    ) %>%
    ungroup() %>%
    arrange(cluster, desc(cells))
}

cluster_sample_table <- tibble()

if (!is.na(cluster_col) && !is.na(sample_col)) {
  cluster_sample_table <- obj@meta.data %>%
    count(
      cluster = .data[[cluster_col]],
      sample = .data[[sample_col]],
      name = "cells"
    ) %>%
    group_by(cluster) %>%
    mutate(
      cluster_total = sum(cells),
      fraction_within_cluster = cells / cluster_total
    ) %>%
    ungroup() %>%
    arrange(cluster, desc(cells))
}

# ------------------------------------------------------------
# UMAP plots
# ------------------------------------------------------------

p_cluster <- make_dimplot(
  obj,
  reduction_name,
  cluster_col,
  "Macrophage RDS: cluster"
)

p_annotation <- make_dimplot(
  obj,
  reduction_name,
  annotation_col,
  "Macrophage RDS: existing annotation"
)

p_sample <- make_dimplot(
  obj,
  reduction_name,
  sample_col,
  "Macrophage RDS: sample",
  label = FALSE
)

p_condition <- make_dimplot(
  obj,
  reduction_name,
  condition_col,
  "Macrophage RDS: condition",
  label = FALSE
)

if (!is.null(p_cluster)) {
  save_pdf(
    p_cluster,
    file.path(pdf_dir, "01_UMAP_Cluster.pdf"),
    12,
    9
  )
}

if (!is.null(p_annotation)) {
  save_pdf(
    p_annotation,
    file.path(pdf_dir, "02_UMAP_ExistingAnnotation.pdf"),
    13,
    9
  )
}

if (!is.null(p_sample)) {
  save_pdf(
    p_sample,
    file.path(pdf_dir, "03_UMAP_Sample.pdf"),
    11,
    9
  )
}

if (!is.null(p_condition)) {
  save_pdf(
    p_condition,
    file.path(pdf_dir, "04_UMAP_Condition.pdf"),
    10,
    8
  )
}

plot_list <- Filter(
  Negate(is.null),
  list(
    Cluster = p_cluster,
    Annotation = p_annotation,
    Sample = p_sample,
    Condition = p_condition
  )
)

if (length(plot_list) > 0) {
  p_combined <- wrap_plots(
    plot_list,
    ncol = 2,
    guides = "collect"
  )

  save_pdf(
    p_combined,
    file.path(pdf_dir, "05_UMAP_AuditCombined.pdf"),
    22,
    17
  )
}

# Split-by-sample annotation UMAP
if (!is.na(annotation_col) && !is.na(sample_col)) {
  p_split <- DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = annotation_col,
    split.by = sample_col,
    pt.size = 0.20,
    raster = TRUE,
    ncol = 3
  ) +
    theme_bw(base_size = 8) +
    theme(legend.text = element_text(size = 7))

  save_pdf(
    p_split,
    file.path(pdf_dir, "06_UMAP_Annotation_SplitBySample.pdf"),
    22,
    16
  )
}

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

write_csv(
  object_summary,
  file.path(csv_dir, "01_ObjectSummary.csv")
)

write_csv(
  assay_summary,
  file.path(csv_dir, "02_AssaySummary.csv")
)

write_csv(
  reduction_summary,
  file.path(csv_dir, "03_ReductionSummary.csv")
)

write_csv(
  metadata_summary,
  file.path(csv_dir, "04_MetadataSummary.csv")
)

write_csv(
  all_counts,
  file.path(csv_dir, "05_MetadataCounts.csv")
)

if (nrow(cluster_annotation_table) > 0) {
  write_csv(
    cluster_annotation_table,
    file.path(csv_dir, "06_ClusterAnnotationComposition.csv")
  )
}

if (nrow(cluster_sample_table) > 0) {
  write_csv(
    cluster_sample_table,
    file.path(csv_dir, "07_ClusterSampleComposition.csv")
  )
}

excel_sheets <- list(
  ObjectSummary = object_summary,
  Assays = assay_summary,
  Reductions = reduction_summary,
  Metadata = metadata_summary,
  Counts = all_counts
)

if (nrow(cluster_annotation_table) > 0) {
  excel_sheets$ClusterAnnotation <- cluster_annotation_table
}

if (nrow(cluster_sample_table) > 0) {
  excel_sheets$ClusterSample <- cluster_sample_table
}

write_excel(
  file.path(output_dir, "Macrophage_RDS_Audit.xlsx"),
  excel_sheets
)

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

capture.output(
  warnings(),
  file = file.path(log_dir, "warnings.txt")
)

saveRDS(
  list(
    detected_cluster_column = cluster_col,
    detected_annotation_column = annotation_col,
    detected_sample_column = sample_col,
    detected_condition_column = condition_col,
    detected_reduction = reduction_name
  ),
  file.path(log_dir, "DetectedFields.rds")
)

cat("============================================================\n")
cat("Macrophage RDS audit completed\n")
cat("============================================================\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main output:\n")
cat("  Macrophage_RDS_Audit.xlsx\n")
cat("  PDF/01_UMAP_Cluster.pdf\n")
cat("  PDF/02_UMAP_ExistingAnnotation.pdf\n")
cat("  PDF/05_UMAP_AuditCombined.pdf\n")
cat("  CSV/04_MetadataSummary.csv\n")
cat("  CSV/06_ClusterAnnotationComposition.csv\n")
