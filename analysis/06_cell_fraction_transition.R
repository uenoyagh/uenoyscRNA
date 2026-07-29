# ============================================================
# 06_cell_fraction_transition.R
# Generic cell/cluster fraction transition pipeline
# uenoy scRNAseq Framework v3.0
# ============================================================

rm(list = ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"

source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "config", "cell_fraction_config.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "utils.R"))
source(file.path(script_root, "R", "cell_fraction_engine.R"))
source(file.path(script_root, "R", "cell_fraction_plot.R"))

required_packages <- c("Seurat", "SeuratObject", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

write_csv_safe_local <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path) || isTRUE(overwrite_existing)) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
}

targets <- cell_fraction_analysis_targets
if (length(targets) == 0L) stop("cell_fraction_analysis_targets is empty.")

for (target_index in seq_along(targets)) {
  analysis_target <- targets[[target_index]]
  target_dir <- get_dataset_dir(analysis_target)

  profile <- cf_get_profile(
    analysis_target,
    cell_fraction_profiles,
    default_profile = cell_fraction_default_profile
  )

  output_root <- cf_get_output_root(
    dataset_name = analysis_target,
    result_folder = cell_fraction_result_folder,
    run_folder = cell_fraction_run_folder,
    create = TRUE
  )

  plot_dirs <- file.path(
    output_root,
    c(
      "01_line_total_fraction",
      "02_line_within_parent",
      "03_line_cell_count",
      "04_stacked_feature",
      "05_stacked_parent",
      "06_heatmap",
      "07_tables",
      "08_QC"
    )
  )
  invisible(lapply(
    plot_dirs, dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  rds_files <- list_rds_files(target_dir, recursive = FALSE)
  if (length(rds_files) == 0L) {
    warning("No RDS files found in: ", target_dir)
    next
  }

  selected_files <- profile$selected_files
  if (!is.null(selected_files)) {
    missing_selected <- setdiff(selected_files, basename(rds_files))
    if (length(missing_selected) > 0L) {
      stop(
        "Selected RDS file(s) not found for ", analysis_target, ":\n",
        paste(missing_selected, collapse = "\n")
      )
    }
    rds_files <- rds_files[basename(rds_files) %in% selected_files]
  }

  cat(
    "\n=== Target ", target_index, "/", length(targets),
    ": ", analysis_target, " ===\n", sep = ""
  )

  for (i in seq_along(rds_files)) {
    rds_path <- rds_files[[i]]
    rds_name <- basename(rds_path)
    file_stub <- cf_safe_filename(tools::file_path_sans_ext(rds_name))

    cat("[", i, "/", length(rds_files), "] ", rds_name, "\n", sep = "")

    object <- safe_read_rds(rds_path)
    if (inherits(object, "rds_read_error")) {
      warning("Failed to read: ", rds_name, "\n", object$message)
      next
    }
    if (!inherits(object, "Seurat")) {
      warning("Not a Seurat object: ", rds_name)
      next
    }

    columns <- cf_resolve_columns(
      object = object,
      feature_override = profile$feature_column_override,
      parent_override = profile$parent_column_override,
      condition_override = profile$condition_column_override,
      sample_override = profile$sample_column_override,
      feature_candidates = profile$feature_candidates,
      parent_candidates = profile$parent_candidates,
      condition_candidates = profile$condition_candidates,
      sample_candidates = profile$sample_candidates
    )

    message(
      "Resolved metadata: feature=", columns$feature,
      "; parent=", columns$parent,
      "; condition=", columns$condition,
      "; sample=", columns$sample
    )

    md <- cf_prepare_metadata(
      object = object,
      columns = columns,
      condition_order = profile$condition_order,
      condition_regex_map = profile$condition_regex_map,
      include_features = profile$include_features,
      exclude_features = profile$exclude_features,
      include_parents = profile$include_parents,
      exclude_parents = profile$exclude_parents,
      drop_na = TRUE
    )

    if (nrow(md) == 0L) {
      warning("No cells remained after filtering: ", rds_name)
      next
    }

    tables <- cf_calculate_tables(md)
    feature_df <- tables$feature
    parent_df <- tables$parent

    feature_df$condition <- cf_order_factor(
      feature_df$condition, profile$condition_order
    )
    parent_df$condition <- cf_order_factor(
      parent_df$condition, profile$condition_order
    )
    feature_df$feature <- factor(
      feature_df$feature,
      levels = cf_natural_levels(feature_df$feature)
    )
    parent_levels <- cf_natural_levels(md$.cf_parent)
    feature_df$parent <- factor(feature_df$parent, levels = parent_levels)
    parent_df$parent <- factor(parent_df$parent, levels = parent_levels)

    prefix <- paste0(file_stub, "_", cf_safe_filename(profile$profile_name))

    write_csv_safe_local(
      feature_df,
      file.path(plot_dirs[[7]], paste0(prefix, "_feature_fraction_table.csv"))
    )
    write_csv_safe_local(
      parent_df,
      file.path(plot_dirs[[7]], paste0(prefix, "_parent_fraction_table.csv"))
    )
    write_csv_safe_local(
      data.frame(
        role = names(columns),
        column = unlist(columns, use.names = FALSE),
        stringsAsFactors = FALSE
      ),
      file.path(plot_dirs[[7]], paste0(prefix, "_resolved_metadata_columns.csv"))
    )
    write_csv_safe_local(
      cf_qc_summary(md, rds_name, analysis_target, profile$profile_name),
      file.path(plot_dirs[[8]], paste0(prefix, "_QC_summary.csv"))
    )

    common_save <- function(p, dir_index, suffix,
                            width = cell_fraction_pdf_width,
                            height = cell_fraction_pdf_height) {
      cf_save_plot_dual(
        plot = p,
        base_path = file.path(
          plot_dirs[[dir_index]],
          paste0(prefix, "_", suffix)
        ),
        width = width,
        height = height,
        export_pdf = cell_fraction_export_pdf,
        export_png = cell_fraction_export_png,
        png_dpi = cell_fraction_png_dpi,
        png_background = cell_fraction_png_background,
        overwrite = overwrite_existing
      )
    }

    denominator_label <- profile$total_denominator_label
    parent_label <- profile$parent_label
    feature_label <- profile$feature_label

    if (isTRUE(cell_fraction_make_line_total)) {
      p <- cf_make_feature_line_plot(
        feature_df,
        y_column = "fraction_total_percent",
        y_label = paste0("Fraction among ", denominator_label, " (%)"),
        title = paste0(feature_label, " fraction transition"),
        subtitle = paste0(
          "Each line = ", feature_label,
          "; panels = ", parent_label,
          "; denominator = ", denominator_label
        ),
        line_width = cell_fraction_line_width,
        point_size = cell_fraction_point_size,
        facet_ncol = cell_fraction_facet_ncol,
        free_y = cell_fraction_free_y,
        label_features = cell_fraction_label_features,
        label_last_only = cell_fraction_label_last_only,
        label_size = cell_fraction_label_size
      )
      common_save(p, 1, "feature_fraction_total_line_by_parent")
    }

    if (isTRUE(cell_fraction_make_line_within_parent)) {
      p <- cf_make_feature_line_plot(
        feature_df,
        y_column = "fraction_within_parent_percent",
        y_label = paste0("Fraction within ", parent_label, " (%)"),
        title = paste0(feature_label, " fraction within ", parent_label),
        subtitle = paste0(
          "Each ", parent_label,
          " sums to 100% within each condition"
        ),
        line_width = cell_fraction_line_width,
        point_size = cell_fraction_point_size,
        facet_ncol = cell_fraction_facet_ncol,
        free_y = cell_fraction_free_y,
        label_features = cell_fraction_label_features,
        label_last_only = cell_fraction_label_last_only,
        label_size = cell_fraction_label_size
      )
      common_save(p, 2, "feature_fraction_within_parent_line")
    }

    if (isTRUE(cell_fraction_make_cell_count_line)) {
      p <- cf_make_feature_line_plot(
        feature_df,
        y_column = "cell_count",
        y_label = "Cell count",
        title = paste0(feature_label, " cell-count transition"),
        subtitle = "Raw cell counts; interpret with sample-size differences in mind",
        line_width = cell_fraction_line_width,
        point_size = cell_fraction_point_size,
        facet_ncol = cell_fraction_facet_ncol,
        free_y = cell_fraction_free_y,
        label_features = cell_fraction_label_features,
        label_last_only = cell_fraction_label_last_only,
        label_size = cell_fraction_label_size
      )
      common_save(p, 3, "feature_cell_count_line_by_parent")
    }

    if (isTRUE(cell_fraction_make_stacked_feature)) {
      common_save(
        cf_make_stacked_feature_plot(
          feature_df,
          denominator_label,
          feature_label,
          parent_label
        ),
        4,
        "feature_fraction_stacked_by_parent"
      )
    }

    if (isTRUE(cell_fraction_make_stacked_parent)) {
      common_save(
        cf_make_stacked_parent_plot(
          parent_df,
          denominator_label,
          parent_label
        ),
        5,
        "parent_fraction_stacked"
      )
    }

    if (isTRUE(cell_fraction_make_heatmap)) {
      heatmap_df <- cf_make_heatmap_table(
        parent_df,
        value_column = cell_fraction_heatmap_value
      )
      write_csv_safe_local(
        heatmap_df,
        file.path(plot_dirs[[7]], paste0(prefix, "_heatmap_table.csv"))
      )
      common_save(
        cf_make_fraction_heatmap(
          heatmap_df,
          value_column = cell_fraction_heatmap_value,
          low_color = cell_fraction_low_color,
          mid_color = cell_fraction_mid_color,
          high_color = cell_fraction_high_color,
          midpoint = cell_fraction_heatmap_midpoint,
          title = paste0(parent_label, " composition heatmap")
        ),
        6,
        "parent_fraction_heatmap",
        width = cell_fraction_heatmap_width,
        height = cell_fraction_heatmap_height
      )
    }

    rm(object, md, tables, feature_df, parent_df)
    gc(verbose = FALSE)
  }
}

cat("\nGeneric cell-fraction pipeline completed.\n")
