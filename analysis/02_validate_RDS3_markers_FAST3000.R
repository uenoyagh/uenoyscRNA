#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(openxlsx)
  library(tibble)
})

# ============================================================
# RDS3 Phase 2-1 FAST3000
# Marker extraction and Top-marker export
# Separate alternative script; does not overwrite FULL-mode outputs
# ============================================================

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS3_validation/",
  "Phase2_Markers_FAST3000"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# This alternative run uses a dedicated directory:
#   Phase2_Markers_FAST3000
# It therefore does not replace or modify Phase2_Markers from FULL mode.
csv_dir <- file.path(output_dir, "CSV")
xlsx_dir <- file.path(output_dir, "Excel")
log_dir <- file.path(output_dir, "Logs")

dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(xlsx_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------

assay_to_use <- "RNA"

cluster_candidates <- c(
  "cluster_for_R8plot_FIXED2",
  "integratedRPCA_snn_res.0.8",
  "seurat_clusters"
)

# FULL: all cells
# FAST: at most 3000 cells per cluster
analysis_mode <- "FAST"

fast_max_cells_per_ident <- 3000

only_positive_markers <- TRUE
minimum_detection_fraction <- 0.10
log2fc_threshold <- 0.25
minimum_cells_group <- 3

# Ranking rule for Top markers:
# 1) adjusted P value
# 2) average log2 fold change
# 3) pct.1 - pct.2
top_n_values <- c(10, 20, 30)

set.seed(260730)

presto_available <- requireNamespace("presto", quietly = TRUE)

if (presto_available) {
  cat("presto detected: Seurat can use the accelerated Wilcoxon implementation.\n")
} else {
  cat(
    "presto is not installed. FAST3000 will still run with Seurat's ",
    "standard Wilcoxon implementation.\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

first_existing <- function(candidates, columns, label) {
  hit <- candidates[candidates %in% columns]

  if (length(hit) == 0) {
    stop(
      paste0(
        "No metadata column found for ", label, ". Candidates: ",
        paste(candidates, collapse = ", ")
      )
    )
  }

  hit[[1]]
}

safe_sheet_name <- function(x) {
  x <- gsub("[\\\\/:*?\\[\\]]", "_", x)
  substr(x, 1, 31)
}

write_dataframe_sheet <- function(wb, sheet_name, x) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, x, withFilter = TRUE)
  freezePane(wb, sheet_name, firstRow = TRUE)

  if (ncol(x) > 0) {
    setColWidths(
      wb,
      sheet = sheet_name,
      cols = seq_len(ncol(x)),
      widths = "auto"
    )
  }
}

# ------------------------------------------------------------
# Load RDS3
# ------------------------------------------------------------

cat("============================================================\n")
cat("RDS3 Phase 2-1 FAST3000: marker extraction\n")
cat("============================================================\n")
cat("RDS file:\n", rds_file, "\n\n")

obj <- readRDS(rds_file)

cluster_col <- first_existing(
  cluster_candidates,
  colnames(obj@meta.data),
  "cluster"
)

cat("Detected cluster column:", cluster_col, "\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n\n")

if (!assay_to_use %in% Assays(obj)) {
  stop(
    paste0(
      "Assay '", assay_to_use,
      "' is not present. Available assays: ",
      paste(Assays(obj), collapse = ", ")
    )
  )
}

DefaultAssay(obj) <- assay_to_use
Idents(obj) <- obj@meta.data[[cluster_col]]

# Natural numeric ordering where possible
ident_levels <- levels(Idents(obj))
ident_numeric <- suppressWarnings(as.numeric(ident_levels))

if (all(!is.na(ident_numeric))) {
  Idents(obj) <- factor(
    as.character(Idents(obj)),
    levels = as.character(sort(unique(ident_numeric)))
  )
}

# ------------------------------------------------------------
# Prepare Seurat v5 layers
# ------------------------------------------------------------

cat("Preparing RNA assay layers...\n")

obj <- tryCatch(
  {
    JoinLayers(obj, assay = assay_to_use)
  },
  error = function(e) {
    cat(
      "JoinLayers was not applied. Continuing with current assay structure.\n",
      "Reason: ", conditionMessage(e), "\n",
      sep = ""
    )
    obj
  }
)

# Check normalized data. Normalize only when necessary.
data_available <- FALSE

data_available <- tryCatch(
  {
    x <- GetAssayData(
      obj,
      assay = assay_to_use,
      layer = "data"
    )
    nrow(x) > 0 && ncol(x) > 0
  },
  error = function(e) FALSE
)

if (!data_available) {
  cat("RNA data layer was unavailable. Running NormalizeData()...\n")

  obj <- NormalizeData(
    obj,
    assay = assay_to_use,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
} else {
  cat("Existing RNA data layer detected. NormalizeData() skipped.\n")
}

# ------------------------------------------------------------
# Cluster summary before DEG
# ------------------------------------------------------------

cluster_summary <- tibble(
  cluster = as.character(Idents(obj))
) %>%
  count(cluster, name = "cells") %>%
  arrange(suppressWarnings(as.numeric(cluster)), cluster)

write_csv(
  cluster_summary,
  file.path(csv_dir, "01_ClusterCellCounts_BeforeMarkers.csv")
)

cat("\nCluster cell counts:\n")
print(cluster_summary, n = Inf)

# ------------------------------------------------------------
# Optional FAST mode
# ------------------------------------------------------------

max_cells_per_ident <- Inf

if (identical(toupper(analysis_mode), "FAST")) {
  max_cells_per_ident <- fast_max_cells_per_ident
}

settings <- tibble(
  setting = c(
    "RDS_file",
    "assay",
    "cluster_column",
    "analysis_mode",
    "max_cells_per_ident",
    "only_positive_markers",
    "minimum_detection_fraction",
    "log2fc_threshold",
    "minimum_cells_group",
    "test_used"
  ),
  value = c(
    rds_file,
    assay_to_use,
    cluster_col,
    analysis_mode,
    as.character(max_cells_per_ident),
    as.character(only_positive_markers),
    as.character(minimum_detection_fraction),
    as.character(log2fc_threshold),
    as.character(minimum_cells_group),
    "wilcox"
  )
)

write_csv(
  settings,
  file.path(csv_dir, "00_MarkerAnalysisSettings.csv")
)

# ------------------------------------------------------------
# Run FindAllMarkers
# ------------------------------------------------------------

cat("\nRunning FindAllMarkers()...\n")
cat("Analysis mode:", analysis_mode, "\n")
cat("max.cells.per.ident:", max_cells_per_ident, "\n")
cat("FAST3000 uses at most 3000 cells per cluster.\n")
cat("FULL-mode output directories are not modified.\n\n")

start_time <- Sys.time()

markers <- FindAllMarkers(
  object = obj,
  assay = assay_to_use,
  only.pos = only_positive_markers,
  test.use = "wilcox",
  min.pct = minimum_detection_fraction,
  logfc.threshold = log2fc_threshold,
  min.cells.group = minimum_cells_group,
  max.cells.per.ident = max_cells_per_ident,
  return.thresh = 1,
  verbose = TRUE
)

end_time <- Sys.time()

elapsed_minutes <- as.numeric(
  difftime(end_time, start_time, units = "mins")
)

cat("\nFindAllMarkers completed.\n")
cat("Elapsed time:", round(elapsed_minutes, 2), "minutes\n")
cat("Marker rows:", nrow(markers), "\n\n")

if (nrow(markers) == 0) {
  stop(
    paste0(
      "FindAllMarkers returned zero rows. ",
      "Check the RNA data layer and cluster identities."
    )
  )
}

# ------------------------------------------------------------
# Standardize output columns
# ------------------------------------------------------------

markers <- markers %>%
  rownames_to_column(var = "row_id")

if (!"gene" %in% colnames(markers)) {
  markers$gene <- markers$row_id
}

markers <- markers %>%
  select(-row_id)

fc_column <- c(
  "avg_log2FC",
  "avg_logFC",
  "avg_diff"
)

fc_column <- fc_column[
  fc_column %in% colnames(markers)
][1]

if (is.na(fc_column)) {
  stop(
    paste0(
      "No fold-change column detected. Columns: ",
      paste(colnames(markers), collapse = ", ")
    )
  )
}

markers <- markers %>%
  mutate(
    cluster = as.character(cluster),
    detection_difference = pct.1 - pct.2
  ) %>%
  arrange(
    cluster,
    p_val_adj,
    desc(.data[[fc_column]]),
    desc(detection_difference)
  )

write_csv(
  markers,
  file.path(csv_dir, "02_AllClusterMarkers.csv")
)

# Significant subset used for ranking
markers_significant <- markers %>%
  filter(
    !is.na(p_val_adj),
    p_val_adj < 0.05,
    .data[[fc_column]] > 0
  )

write_csv(
  markers_significant,
  file.path(csv_dir, "03_SignificantPositiveMarkers.csv")
)

# ------------------------------------------------------------
# Top marker tables
# ------------------------------------------------------------

top_marker_tables <- list()

for (top_n in top_n_values) {

  top_table <- markers_significant %>%
    group_by(cluster) %>%
    arrange(
      p_val_adj,
      desc(.data[[fc_column]]),
      desc(detection_difference),
      .by_group = TRUE
    ) %>%
    slice_head(n = top_n) %>%
    ungroup()

  top_marker_tables[[paste0("Top", top_n)]] <- top_table

  write_csv(
    top_table,
    file.path(
      csv_dir,
      paste0("04_Top", top_n, "Markers_PerCluster.csv")
    )
  )
}

# ------------------------------------------------------------
# Marker counts per cluster
# ------------------------------------------------------------

marker_count_summary <- markers %>%
  group_by(cluster) %>%
  summarise(
    tested_positive_markers = n(),
    significant_positive_markers = sum(
      p_val_adj < 0.05,
      na.rm = TRUE
    ),
    strong_markers_log2FC_ge_0_5 = sum(
      p_val_adj < 0.05 &
        .data[[fc_column]] >= 0.5,
      na.rm = TRUE
    ),
    strong_markers_log2FC_ge_1 = sum(
      p_val_adj < 0.05 &
        .data[[fc_column]] >= 1,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    cluster_summary,
    by = "cluster"
  ) %>%
  arrange(
    suppressWarnings(as.numeric(cluster)),
    cluster
  )

write_csv(
  marker_count_summary,
  file.path(csv_dir, "05_MarkerCounts_PerCluster.csv")
)

# ------------------------------------------------------------
# Excel: all markers and summary
# ------------------------------------------------------------

excel_row_limit <- 1048576L

wb_all <- createWorkbook()

write_dataframe_sheet(
  wb_all,
  "Settings",
  settings
)

write_dataframe_sheet(
  wb_all,
  "ClusterCounts",
  cluster_summary
)

write_dataframe_sheet(
  wb_all,
  "MarkerCounts",
  marker_count_summary
)

if (nrow(markers) <= excel_row_limit - 1L) {
  write_dataframe_sheet(
    wb_all,
    "AllMarkers",
    markers
  )
} else {
  cat(
    "AllMarkers exceeds the Excel row limit. ",
    "Writing cluster-specific sheets instead.\n",
    sep = ""
  )

  for (cl in unique(markers$cluster)) {
    x <- markers %>%
      filter(cluster == cl)

    write_dataframe_sheet(
      wb_all,
      safe_sheet_name(paste0("Cluster_", cl)),
      x
    )
  }
}

saveWorkbook(
  wb_all,
  file.path(xlsx_dir, "02_AllClusterMarkers.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------------------
# Excel: Top 10, 20 and 30
# ------------------------------------------------------------

for (top_name in names(top_marker_tables)) {

  top_table <- top_marker_tables[[top_name]]
  wb_top <- createWorkbook()

  write_dataframe_sheet(
    wb_top,
    top_name,
    top_table
  )

  for (cl in unique(top_table$cluster)) {
    x <- top_table %>%
      filter(cluster == cl)

    write_dataframe_sheet(
      wb_top,
      safe_sheet_name(paste0("Cluster_", cl)),
      x
    )
  }

  saveWorkbook(
    wb_top,
    file.path(
      xlsx_dir,
      paste0("03_", top_name, "Markers_PerCluster.xlsx")
    ),
    overwrite = TRUE
  )
}

# ------------------------------------------------------------
# Compact gene matrix for manual review
# ------------------------------------------------------------

top30 <- top_marker_tables[["Top30"]]

top30_compact <- top30 %>%
  group_by(cluster) %>%
  summarise(
    top30_genes = paste(gene, collapse = "; "),
    .groups = "drop"
  )

write_csv(
  top30_compact,
  file.path(csv_dir, "06_Top30Genes_Compact.csv")
)

wb_compact <- createWorkbook()

write_dataframe_sheet(
  wb_compact,
  "Top30Compact",
  top30_compact
)

write_dataframe_sheet(
  wb_compact,
  "MarkerCounts",
  marker_count_summary
)

saveWorkbook(
  wb_compact,
  file.path(xlsx_dir, "04_Top30Genes_Compact.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------------------
# Runtime and session information
# ------------------------------------------------------------

runtime_summary <- tibble(
  item = c(
    "start_time",
    "end_time",
    "elapsed_minutes",
    "marker_rows",
    "significant_marker_rows",
    "clusters"
  ),
  value = c(
    as.character(start_time),
    as.character(end_time),
    as.character(elapsed_minutes),
    as.character(nrow(markers)),
    as.character(nrow(markers_significant)),
    as.character(length(unique(markers$cluster)))
  )
)

write_csv(
  runtime_summary,
  file.path(csv_dir, "07_RuntimeSummary.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

capture.output(
  warnings(),
  file = file.path(log_dir, "warnings.txt")
)

cat("\n============================================================\n")
cat("RDS3 Phase 2-1 FAST3000 completed\n")
cat("============================================================\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main outputs:\n")
cat("  CSV/02_AllClusterMarkers.csv\n")
cat("  CSV/04_Top30Markers_PerCluster.csv\n")
cat("  CSV/05_MarkerCounts_PerCluster.csv\n")
cat("  CSV/06_Top30Genes_Compact.csv\n")
cat("  Excel/02_AllClusterMarkers.xlsx\n")
cat("  Excel/03_Top10Markers_PerCluster.xlsx\n")
cat("  Excel/03_Top20Markers_PerCluster.xlsx\n")
cat("  Excel/03_Top30Markers_PerCluster.xlsx\n")
cat("  Excel/04_Top30Genes_Compact.xlsx\n")
cat("\nNext phase: marker DotPlot, FeaturePlot and heatmap.\n")
