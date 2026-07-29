# ============================================================
# Cluster-fraction transition configuration
# uenoy scRNAseq Framework v2.6
# ============================================================

# Dataset key used by project_config.R / get_dataset_dir().
# Recommended for the present macrophage analysis:
cluster_fraction_analysis_target <- "Mouse_MASH_Mphi_RDS"

# NULL = analyze all RDS files directly under the dataset directory.
cluster_fraction_selected_files <- NULL

# Metadata overrides. NULL enables automatic detection.
cluster_fraction_cluster_column_override <- NULL
cluster_fraction_annotation_column_override <- NULL
cluster_fraction_condition_column_override <- NULL
cluster_fraction_sample_column_override <- NULL

# Preferred order. Values absent from an RDS are ignored; additional values
# are appended after these levels.
cluster_fraction_condition_order <- c("STD", "CDAHFD", "Sham", "Tx")

# Optional mapping when condition labels are not already present.
# Names are regex patterns matched against sample labels; values are conditions.
cluster_fraction_condition_regex_map <- c(
  "STD" = "STD",
  "CDAHFD|CDHFD" = "CDAHFD",
  "Sham" = "Sham",
  "Tx" = "Tx"
)

# Plot generation switches.
cluster_fraction_make_line_total <- TRUE
cluster_fraction_make_line_within_annotation <- TRUE
cluster_fraction_make_stacked_cluster <- TRUE
cluster_fraction_make_stacked_annotation <- TRUE
cluster_fraction_make_cell_count_line <- TRUE

# Styling.
cluster_fraction_line_width <- 0.55
cluster_fraction_point_size <- 1.8
cluster_fraction_label_clusters <- TRUE
cluster_fraction_label_last_only <- TRUE
cluster_fraction_label_size <- 2.6
cluster_fraction_facet_ncol <- 2
cluster_fraction_free_y <- TRUE

# Dimensions and export.
cluster_fraction_pdf_width <- 14
cluster_fraction_pdf_height <- 10
cluster_fraction_png_dpi <- 300
cluster_fraction_png_background <- "white"
cluster_fraction_export_pdf <- TRUE
cluster_fraction_export_png <- TRUE

# Reuse local_config.R overwrite_existing when available.
if (!exists("overwrite_existing")) overwrite_existing <- FALSE
