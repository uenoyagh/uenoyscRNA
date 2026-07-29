# ============================================================
# Generic cell-fraction configuration
# uenoy scRNAseq Framework v3.0
# ============================================================

# Run one or more of the four active dataset directories.
cell_fraction_analysis_targets <- c(
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_RDS",
  "Human_MASH_RDS",
  "Human_MASH_Mphi_RDS"
)

# Output below get_result_root(dataset_name).
cell_fraction_result_folder <- "09_cell_fraction"
cell_fraction_run_folder <- "cell_fraction_transition_v3.0"

# Plot switches.
cell_fraction_make_line_total <- TRUE
cell_fraction_make_line_within_parent <- TRUE
cell_fraction_make_cell_count_line <- TRUE
cell_fraction_make_stacked_feature <- TRUE
cell_fraction_make_stacked_parent <- TRUE
cell_fraction_make_heatmap <- TRUE

# Plot appearance.
cell_fraction_line_width <- 0.55
cell_fraction_point_size <- 1.8
cell_fraction_label_features <- TRUE
cell_fraction_label_last_only <- TRUE
cell_fraction_label_size <- 2.6
cell_fraction_facet_ncol <- 2
cell_fraction_free_y <- TRUE

cell_fraction_pdf_width <- 14
cell_fraction_pdf_height <- 10
cell_fraction_heatmap_width <- 10
cell_fraction_heatmap_height <- 8
cell_fraction_png_dpi <- 300
cell_fraction_png_background <- "white"
cell_fraction_export_pdf <- TRUE
cell_fraction_export_png <- TRUE

# Heatmap.
cell_fraction_heatmap_value <- "fraction_total_percent"
cell_fraction_low_color <- "#0033FF"
cell_fraction_mid_color <- "#FFFFFF"
cell_fraction_high_color <- "#FF1A1A"
cell_fraction_heatmap_midpoint <- NULL

if (!exists("overwrite_existing")) overwrite_existing <- FALSE

# ------------------------------------------------------------
# Shared candidate columns
# ------------------------------------------------------------

cf_cluster_candidates <- c(
  "integratedRPCA_snn_res.3.0",
  "integratedRPCA_snn_res.3",
  "integrated_snn_res.3.0",
  "RNA_snn_res.3.0",
  "RNA_snn_res.3",
  "seurat_clusters",
  "cluster",
  "Cluster",
  "res3.0",
  "res3"
)

cf_annotation_candidates <- c(
  "annotation_group",
  "annotation_group_final",
  "cell_annotation",
  "celltype",
  "cell_type",
  "CellType",
  "layer2",
  "Layer2",
  "layer1",
  "Layer1",
  "annotation",
  "Annotation"
)

cf_condition_candidates <- c(
  "condition",
  "Condition",
  "group",
  "Group",
  "treatment",
  "Treatment",
  "diet",
  "Diet",
  "orig.ident"
)

cf_sample_candidates <- c(
  "sample",
  "Sample",
  "sample_id",
  "SampleID",
  "orig.ident"
)

# ------------------------------------------------------------
# Dataset-specific profiles
#
# feature:
#   the individual lines/stacked components, e.g. cluster or cell type
# parent:
#   the facet/upper-level category, e.g. macrophage annotation
#
# For whole-cell datasets, feature and parent may intentionally be the
# same cell-type column. This produces one line per cell type and a
# whole-cell composition stacked plot/heatmap.
# ------------------------------------------------------------

cell_fraction_profiles <- list(

  Mouse_MASH_RDS = list(
    profile_name = "all_cell",
    selected_files = NULL,
    feature_column_override = NULL,
    parent_column_override = NULL,
    condition_column_override = NULL,
    sample_column_override = NULL,
    feature_candidates = cf_annotation_candidates,
    parent_candidates = cf_annotation_candidates,
    condition_candidates = cf_condition_candidates,
    sample_candidates = cf_sample_candidates,
    condition_order = c("STD", "CDAHFD", "Sham", "Tx"),
    condition_regex_map = c(
      "STD" = "STD",
      "CDAHFD|CDHFD" = "CDAHFD",
      "Sham" = "Sham",
      "Tx" = "Tx"
    ),
    include_features = NULL,
    exclude_features = NULL,
    include_parents = NULL,
    exclude_parents = NULL,
    feature_label = "cell type",
    parent_label = "cell type",
    total_denominator_label = "total cells"
  ),

  Mouse_MASH_Mphi_RDS = list(
    profile_name = "macrophage_cluster",
    selected_files = NULL,
    feature_column_override = NULL,
    parent_column_override = NULL,
    condition_column_override = NULL,
    sample_column_override = NULL,
    feature_candidates = cf_cluster_candidates,
    parent_candidates = cf_annotation_candidates,
    condition_candidates = cf_condition_candidates,
    sample_candidates = cf_sample_candidates,
    condition_order = c("STD", "CDAHFD", "Sham", "Tx"),
    condition_regex_map = c(
      "STD" = "STD",
      "CDAHFD|CDHFD" = "CDAHFD",
      "Sham" = "Sham",
      "Tx" = "Tx"
    ),
    include_features = NULL,
    exclude_features = NULL,
    include_parents = NULL,
    exclude_parents = NULL,
    feature_label = "cluster",
    parent_label = "macrophage annotation",
    total_denominator_label = "total Mphi"
  ),

  Human_MASH_RDS = list(
    profile_name = "all_cell",
    selected_files = NULL,
    feature_column_override = NULL,
    parent_column_override = NULL,
    condition_column_override = NULL,
    sample_column_override = NULL,
    feature_candidates = cf_annotation_candidates,
    parent_candidates = cf_annotation_candidates,
    condition_candidates = c(
      "NAS",
      "nas",
      "condition",
      "Condition",
      "group",
      "Group",
      "orig.ident"
    ),
    sample_candidates = cf_sample_candidates,
    condition_order = c("0", "1", "2", "3", "4", "5", "6"),
    condition_regex_map = NULL,
    include_features = NULL,
    exclude_features = NULL,
    include_parents = NULL,
    exclude_parents = NULL,
    feature_label = "cell type",
    parent_label = "cell type",
    total_denominator_label = "total cells"
  ),

  Human_MASH_Mphi_RDS = list(
    profile_name = "macrophage_cluster",
    selected_files = NULL,
    feature_column_override = NULL,
    parent_column_override = NULL,
    condition_column_override = NULL,
    sample_column_override = NULL,
    feature_candidates = cf_cluster_candidates,
    parent_candidates = cf_annotation_candidates,
    condition_candidates = c(
      "NAS",
      "nas",
      "condition",
      "Condition",
      "group",
      "Group",
      "orig.ident"
    ),
    sample_candidates = cf_sample_candidates,
    condition_order = c("0", "1", "2", "3", "4", "5", "6"),
    condition_regex_map = NULL,
    include_features = NULL,
    exclude_features = NULL,
    include_parents = NULL,
    exclude_parents = NULL,
    feature_label = "cluster",
    parent_label = "macrophage annotation",
    total_denominator_label = "total Mphi"
  )
)

cell_fraction_default_profile <- cell_fraction_profiles$Mouse_MASH_RDS
