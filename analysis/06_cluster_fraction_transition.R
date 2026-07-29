# ============================================================
# 06_cluster_fraction_transition.R
# Cluster count/fraction transition module
# uenoy scRNAseq Framework v2.6
# ============================================================

rm(list = ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"

source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "config", "cluster_fraction_config.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "utils.R"))
source(file.path(script_root, "R", "cluster_fraction_engine.R"))
source(file.path(script_root, "R", "cluster_fraction_plot.R"))

required_packages <- c("Seurat", "SeuratObject", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))

analysis_target <- cluster_fraction_analysis_target
target_dir <- get_dataset_dir(analysis_target)
output_root <- file.path(
  get_result_dir(
    analysis_target,
    "module_score",
    create = TRUE
  ),
  "cluster_fraction_transition_v2.6"
)
plot_dirs <- file.path(output_root, c("01_line_total_fraction", "02_line_within_annotation",
                                      "03_line_cell_count", "04_stacked_cluster",
                                      "05_stacked_annotation", "06_tables"))
invisible(lapply(plot_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

rds_files <- list_rds_files(target_dir, recursive = FALSE)
if (!is.null(cluster_fraction_selected_files)) {
  missing_selected <- setdiff(cluster_fraction_selected_files, basename(rds_files))
  if (length(missing_selected) > 0) stop("Selected RDS file(s) not found:\n", paste(missing_selected, collapse = "\n"))
  rds_files <- rds_files[basename(rds_files) %in% cluster_fraction_selected_files]
}

write_csv_safe_local <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path) || isTRUE(overwrite_existing)) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
}

for (i in seq_along(rds_files)) {
  rds_path <- rds_files[[i]]
  rds_name <- basename(rds_path)
  file_stub <- tools::file_path_sans_ext(rds_name)
  file_stub <- gsub("[^A-Za-z0-9._-]+", "_", file_stub)
  cat("[", i, "/", length(rds_files), "] ", rds_name, "\n", sep = "")

  object <- safe_read_rds(rds_path)
  if (inherits(object, "rds_read_error")) { warning("Failed to read: ", rds_name); next }
  if (!inherits(object, "Seurat")) { warning("Not a Seurat object: ", rds_name); next }

  columns <- cf_resolve_columns(
    object,
    cluster_fraction_cluster_column_override,
    cluster_fraction_annotation_column_override,
    cluster_fraction_condition_column_override,
    cluster_fraction_sample_column_override
  )
  message("Resolved metadata: cluster=", columns$cluster,
          "; annotation=", columns$annotation,
          "; condition=", columns$condition,
          "; sample=", columns$sample)

  md <- cf_prepare_metadata(
    object, columns,
    cluster_fraction_condition_order,
    cluster_fraction_condition_regex_map
  )
  tables <- cf_calculate_tables(md)
  cluster_df <- tables$cluster
  annotation_df <- tables$annotation

  # Restore factor ordering after aggregate/merge.
  cluster_df$condition <- cf_order_factor(cluster_df$condition, cluster_fraction_condition_order)
  annotation_df$condition <- cf_order_factor(annotation_df$condition, cluster_fraction_condition_order)
  cluster_df$cluster <- factor(cluster_df$cluster, levels = cf_natural_levels(cluster_df$cluster))
  cluster_df$annotation <- factor(cluster_df$annotation, levels = unique(as.character(md$.cf_annotation)))
  annotation_df$annotation <- factor(annotation_df$annotation, levels = unique(as.character(md$.cf_annotation)))

  write_csv_safe_local(cluster_df, file.path(plot_dirs[[6]], paste0(file_stub, "_cluster_fraction_table.csv")))
  write_csv_safe_local(annotation_df, file.path(plot_dirs[[6]], paste0(file_stub, "_annotation_fraction_table.csv")))
  write_csv_safe_local(data.frame(role = names(columns), column = unlist(columns), row.names = NULL),
                       file.path(plot_dirs[[6]], paste0(file_stub, "_resolved_metadata_columns.csv")))

  common_save <- function(p, dir_index, suffix, width = cluster_fraction_pdf_width, height = cluster_fraction_pdf_height) {
    cf_save_plot_dual(
      p, file.path(plot_dirs[[dir_index]], paste0(file_stub, "_", suffix)),
      width, height,
      cluster_fraction_export_pdf, cluster_fraction_export_png,
      cluster_fraction_png_dpi, cluster_fraction_png_background,
      overwrite_existing
    )
  }

  if (isTRUE(cluster_fraction_make_line_total)) {
    p <- cf_make_cluster_line_plot(
      cluster_df, "fraction_total_percent", "Fraction among total cells (%)",
      "Cluster fraction transition by macrophage annotation",
      "Each line = cluster; y-axis free by annotation group",
      cluster_fraction_line_width, cluster_fraction_point_size,
      cluster_fraction_facet_ncol, cluster_fraction_free_y,
      cluster_fraction_label_clusters, cluster_fraction_label_last_only,
      cluster_fraction_label_size
    )
    common_save(p, 1, "cluster_fraction_total_line_by_annotation")
  }

  if (isTRUE(cluster_fraction_make_line_within_annotation)) {
    p <- cf_make_cluster_line_plot(
      cluster_df, "fraction_within_annotation_percent", "Fraction within annotation (%)",
      "Cluster fraction within each macrophage annotation",
      "Each annotation sums to 100% within each condition",
      cluster_fraction_line_width, cluster_fraction_point_size,
      cluster_fraction_facet_ncol, cluster_fraction_free_y,
      cluster_fraction_label_clusters, cluster_fraction_label_last_only,
      cluster_fraction_label_size
    )
    common_save(p, 2, "cluster_fraction_within_annotation_line")
  }

  if (isTRUE(cluster_fraction_make_cell_count_line)) {
    p <- cf_make_cluster_line_plot(
      cluster_df, "cell_count", "Cell count",
      "Cluster cell-count transition by macrophage annotation",
      "Raw cell counts; interpret with library-size differences in mind",
      cluster_fraction_line_width, cluster_fraction_point_size,
      cluster_fraction_facet_ncol, cluster_fraction_free_y,
      cluster_fraction_label_clusters, cluster_fraction_label_last_only,
      cluster_fraction_label_size
    )
    common_save(p, 3, "cluster_cell_count_line_by_annotation")
  }

  if (isTRUE(cluster_fraction_make_stacked_cluster)) {
    common_save(cf_make_stacked_cluster_plot(cluster_df), 4, "cluster_fraction_stacked_by_annotation")
  }
  if (isTRUE(cluster_fraction_make_stacked_annotation)) {
    common_save(cf_make_stacked_annotation_plot(annotation_df), 5, "annotation_fraction_stacked")
  }

  rm(object, md, tables, cluster_df, annotation_df)
  gc(verbose = FALSE)
}

cat("\nCluster-fraction pipeline completed.\n")
cat("Output root:\n", output_root, "\n")
