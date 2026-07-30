#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath(
    "analysis/RDS3_annotation_visualization_v4.1.0/99_run_all.R",
    mustWork = FALSE
  )
}

analysis_dir <- dirname(script_path)
project_root <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = FALSE)

required_files <- c(
  file.path(project_root, "config", "config.R"),
  file.path(project_root, "R", "01_utils_v3.2.R"),
  file.path(project_root, "R", "04_plot_export_v3.2.R"),
  file.path(project_root, "R", "04_plot_export_v4.0.R"),
  file.path(project_root, "R", "04_plot_export_v4.0.2.R"),
  file.path(project_root, "R", "05_split_umap_violin_v4.0.3.R"),
  file.path(project_root, "R", "08_saturated_umap_one_marker_violin_v4.0.9.R"),
  file.path(project_root, "R", "09_clear_group_umap_v4.0.9.R"),
  file.path(project_root, "R", "10_dense_umap_singlepage_violin_v4.1.0.R")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Required files were not found:\n", paste(missing_files, collapse = "\n"))
}

invisible(lapply(required_files, source))

assert_packages(c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "openxlsx"
))

external_root <- detect_external_root(CFG$external_ssd_candidates)
output_parent <- file.path(external_root, dirname(CFG$output_relative_dir))

candidate_rds <- c(
  file.path(
    output_parent,
    "RDS3_annotation_visualization_v4.0.9",
    "objects",
    "RDS3_with_visualization_metadata_v4.0.9.rds"
  ),
  file.path(
    output_parent,
    "RDS3_annotation_visualization_v4.0.8",
    "objects",
    "RDS3_with_visualization_metadata_v4.0.8.rds"
  ),
  file.path(
    output_parent,
    "RDS3_annotation_visualization_v4.0.5",
    "objects",
    "RDS3_with_visualization_metadata_v4.0.5.rds"
  )
)

existing_rds <- candidate_rds[file.exists(candidate_rds)]
if (length(existing_rds) == 0L) {
  stop("Input RDS was not found:\n", paste(candidate_rds, collapse = "\n"))
}

rds_path <- existing_rds[[1]]

output_dir <- file.path(
  output_parent,
  "RDS3_annotation_visualization_v4.1.0"
)

fig_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")
log_dir <- file.path(output_dir, "logs")

ensure_dirs(c(output_dir, fig_dir, table_dir, object_dir, log_dir))

message("Input RDS: ", rds_path)
message("Reading RDS...")

object <- readRDS(rds_path)
reduction <- CFG$reduction
created_at <- Sys.time()

if (!"vote_ueno_summary_v40" %in% colnames(object[[]])) {
  object <- derive_summary_annotation_v40(object)
}

if (!"analysis_group_v408" %in% colnames(object[[]])) {
  object <- derive_analysis_group_v408(
    object = object,
    sample_col = CFG$sample_col,
    output_col = "analysis_group_v408"
  )
}

object <- derive_target_cellclass_v410(
  object = object,
  source_col = "vote_ueno_summary_v40",
  output_col = "target_cellclass_v410"
)

# ============================================================
# 1. Denser integrated group UMAP
# Colors and font sizes unchanged from v4.0.9
# ============================================================

group_levels <- levels(object$analysis_group_v408)
group_ncol <- min(2L, max(1L, length(group_levels)))

p_group <- publish_split_umap_clear_v408(
  object = object,
  group_by = "vote_ueno_summary_v40",
  split_by = "analysis_group_v408",
  reduction = reduction,
  palette = ueno_subtype_palette_v407,
  title = "RDS3: integrated UMAP by analysis group",
  pt_size = 1.18,
  raster = TRUE,
  raster_dpi = c(600, 600),
  source_rds = rds_path,
  created_at = created_at,
  ncol = group_ncol,
  legend_text_size = 9.5,
  strip_text_size = 13
)

save_pdf(
  p_group,
  file.path(
    fig_dir,
    "UMAP_summary_integrated_by_group_dense_v4.1.0.pdf"
  ),
  width = max(14, 6.8 * group_ncol),
  height = max(11, 6.2 * ceiling(length(group_levels) / group_ncol))
)

save_each_integrated_group_umap_v408(
  object = object,
  output_pdf = file.path(
    fig_dir,
    "UMAP_summary_each_integrated_group_dense_multipage_v4.1.0.pdf"
  ),
  analysis_group_col = "analysis_group_v408",
  annotation_col = "vote_ueno_summary_v40",
  reduction = reduction,
  palette = ueno_subtype_palette_v407,
  pt_size = 1.30,
  label_size = 4.4,
  source_rds = rds_path,
  created_at = created_at
)

for (group_value in group_levels) {
  p_one <- publish_single_group_umap_v408(
    object = object,
    group_value = group_value,
    analysis_group_col = "analysis_group_v408",
    annotation_col = "vote_ueno_summary_v40",
    reduction = reduction,
    palette = ueno_subtype_palette_v407,
    pt_size = 1.30,
    label_size = 4.4,
    source_rds = rds_path,
    created_at = created_at
  )

  save_pdf(
    p_one,
    file.path(
      fig_dir,
      paste0("UMAP_summary_integrated_", group_value, "_dense_v4.1.0.pdf")
    ),
    width = 12,
    height = 11
  )
}

# ============================================================
# 2. One marker per requested cell class
# ============================================================

assay_to_use <- Seurat::DefaultAssay(object)

marker_table <- find_target_markers_v410(
  object = object,
  identity_col = "target_cellclass_v410",
  assay = assay_to_use,
  slot = "data",
  min_pct = 0.15,
  logfc_threshold = 0.15
)

safe_write_csv(
  marker_table,
  file.path(
    table_dir,
    "One_marker_per_selected_cellclass_v4.1.0.csv"
  )
)

write_excel_report(
  file.path(
    table_dir,
    "One_marker_per_selected_cellclass_v4.1.0.xlsx"
  ),
  list(
    Selected_cellclass_markers = marker_table
  )
)

detected_classes <- levels(object$target_cellclass_v410)
missing_classes <- setdiff(target_cell_order_v410, detected_classes)

safe_write_csv(
  data.frame(
    requested_cellclass = target_cell_order_v410,
    detected = target_cell_order_v410 %in% detected_classes,
    stringsAsFactors = FALSE
  ),
  file.path(
    table_dir,
    "Selected_cellclass_detection_status_v4.1.0.csv"
  )
)

if (length(missing_classes) > 0L) {
  warning(
    "Requested cell classes not detected and omitted from violin plot: ",
    paste(missing_classes, collapse = ", ")
  )
}

save_singlepage_target_violin_v410(
  object = object,
  marker_table = marker_table,
  output_pdf = file.path(
    fig_dir,
    "Violin_selected_cellclasses_onepage_by_group_v4.1.0.pdf"
  ),
  identity_col = "target_cellclass_v410",
  analysis_group_col = "analysis_group_v408",
  assay = assay_to_use,
  slot = "data",
  source_rds = rds_path,
  created_at = created_at,
  group_order = c("STD", "CDHFD", "Sham", "Tx"),
  width = 22,
  height = 12.5
)

saveRDS(
  object,
  file.path(
    object_dir,
    "RDS3_with_visualization_metadata_v4.1.0.rds"
  ),
  compress = FALSE
)

manifest <- data.frame(
  item = c(
    "version",
    "input_rds",
    "integrated_group_point_size",
    "individual_group_point_size",
    "individual_group_label_size",
    "violin_layout",
    "requested_cellclasses",
    "detected_cellclasses",
    "missing_cellclasses",
    "created_at"
  ),
  value = c(
    "v4.1.0",
    normalizePath(rds_path, mustWork = FALSE),
    "1.18",
    "1.30",
    "4.4",
    "single-page 22 x 12.5 inch; rows=genes; columns=STD/CDHFD/Sham/Tx",
    paste(target_cell_order_v410, collapse = ", "),
    paste(detected_classes, collapse = ", "),
    paste(missing_classes, collapse = ", "),
    format(created_at, "%Y-%m-%d %H:%M:%S %Z")
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  manifest,
  file.path(log_dir, "run_manifest_v4.1.0.csv")
)

message("")
message("RDS3 visualization v4.1.0 completed.")
message("Output directory: ", output_dir)
