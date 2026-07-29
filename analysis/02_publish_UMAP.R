# ============================================================
# 02_publish_UMAP.R
# uenoy scRNAseq Framework v2.3
# ============================================================

rm(list = ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"

source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "utils.R"))
source(file.path(script_root, "R", "umap.R"))
source(file.path(script_root, "R", "plotting.R"))

required_packages <- c("Seurat", "SeuratObject", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

target_dir <- get_dataset_dir(analysis_target)
umap_root <- get_result_dir(analysis_target, "umap", create = TRUE)
output_dirs <- make_umap_output_subdirs(umap_root)
rds_files <- list_rds_files(target_dir, recursive = FALSE)

if (!is.null(umap_selected_files)) {
  rds_files <- rds_files[basename(rds_files) %in% umap_selected_files]
}

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(
  output_dirs[["logs"]],
  paste0("publish_UMAP_v2.3_", run_stamp, ".log")
)

sink(log_file, split = TRUE)
on.exit({
  while (sink.number() > 0) sink()
}, add = TRUE)

cat("uenoy scRNAseq Framework: Publication UMAP\n")
cat("Framework version: ", project_config$framework_version, "\n", sep = "")
cat("Started: ", timestamp_string(), "\n", sep = "")
cat("Target: ", analysis_target, "\n\n", sep = "")

resolver_rows <- list()

for (i in seq_along(rds_files)) {
  rds_path <- rds_files[[i]]
  rds_name <- basename(rds_path)
  file_stub <- safe_umap_filename(rds_name)

  cat("[", i, "/", length(rds_files), "] ", rds_name, "\n", sep = "")

  object <- safe_read_rds(rds_path)
  if (inherits(object, "rds_read_error")) next
  if (!inherits(object, "Seurat")) next

  reduction <- resolve_umap_reduction(object, umap_reduction_override)
  sample_column <- resolve_sample_column(object, umap_sample_column_override)
  condition_column <- resolve_condition_column(object, umap_condition_column_override)
  annotation_column <- resolve_annotation_column(object, umap_annotation_column_override)
  cluster_column <- resolve_cluster_column(object, umap_cluster_column_override)

  for (column in c(
    sample_column, condition_column,
    annotation_column, cluster_column
  )) {
    if (!is.na(column)) {
      object <- prepare_metadata_factor(object, column)
    }
  }

  footer_text <- build_umap_footer(
    rds_name,
    project_config$framework_name,
    project_config$framework_version,
    umap_footer_include_rds,
    umap_footer_include_created,
    umap_footer_include_framework
  )

  umap_df <- extract_umap_dataframe(
    object,
    reduction,
    c(
      sample_column, condition_column,
      annotation_column, cluster_column
    )
  )

  resolver_rows[[length(resolver_rows) + 1]] <- data.frame(
    file = rds_name,
    reduction = reduction,
    sample_column = sample_column,
    condition_column = condition_column,
    annotation_column = annotation_column,
    cluster_column = cluster_column,
    stringsAsFactors = FALSE
  )

  palettes <- c()
  if (isTRUE(umap_generate_R8)) palettes <- c(palettes, "R8")
  if (isTRUE(umap_generate_pastel)) palettes <- c(palettes, "Pastel")

  for (palette_style in palettes) {
    prefix <- if (palette_style == "R8") "r8_" else "pastel_"

    make_standard <- function(
        group_column, title_text, out_key,
        suffix, label = FALSE, hide_legend = FALSE) {

      if (is.na(group_column)) return(invisible(FALSE))

      p <- make_dimplot(
        object = object,
        reduction = reduction,
        group_by = group_column,
        palette_style = palette_style,
        title = title_text,
        point_size = umap_point_size,
        raster = umap_raster,
        raster_dpi = umap_raster_dpi,
        shuffle = umap_shuffle,
        label = label,
        label_size = umap_label_size,
        hide_legend = hide_legend,
        legend_point_size = umap_legend_point_size,
        legend_text_size = umap_legend_text_size,
        title_size = umap_title_size,
        axis_text_size = umap_axis_text_size,
        footer_text = footer_text,
        footer_text_size = umap_footer_text_size
      )

      save_umap_pdf(
        p,
        file.path(
          output_dirs[[paste0(prefix, out_key)]],
          paste0(file_stub, "_", suffix, "_", palette_style, ".pdf")
        ),
        umap_pdf_width,
        umap_pdf_height,
        overwrite_existing
      )
    }

    overview_group <- c(
      annotation_column, cluster_column,
      condition_column, sample_column
    )
    overview_group <- overview_group[!is.na(overview_group)]

    if (umap_make_overview && length(overview_group) > 0) {
      g <- overview_group[[1]]
      is_cluster <- identical(g, cluster_column)
      is_annotation <- identical(g, annotation_column)

      make_standard(
        g,
        if (is_cluster) "Cluster" else if (is_annotation) "Cell annotation" else g,
        "overview",
        paste0("UMAP_overview_", g),
        label = (is_cluster && umap_label_cluster) ||
          (is_annotation && umap_label_annotation),
        hide_legend = (is_cluster && umap_hide_cluster_legend) ||
          (is_annotation && umap_hide_annotation_legend)
      )
    }

    if (umap_make_sample && !is.na(sample_column)) {
      make_standard(
        sample_column, "Sample", "sample",
        paste0("UMAP_sample_", sample_column)
      )
    }

    if (umap_make_condition && !is.na(condition_column)) {
      make_standard(
        condition_column, "Condition", "condition",
        paste0("UMAP_condition_", condition_column)
      )
    }

    if (umap_make_annotation && !is.na(annotation_column)) {
      make_standard(
        annotation_column, "Cell annotation", "annotation",
        paste0("UMAP_annotation_", annotation_column),
        label = umap_label_annotation,
        hide_legend = umap_hide_annotation_legend
      )
    }

    if (umap_make_cluster && !is.na(cluster_column)) {
      make_standard(
        cluster_column, "Cluster", "cluster",
        paste0("UMAP_cluster_", cluster_column),
        label = umap_label_cluster,
        hide_legend = umap_hide_cluster_legend
      )
    }

    split_sample_group <- c(
      annotation_column, cluster_column, condition_column
    )
    split_sample_group <- split_sample_group[
      !is.na(split_sample_group) &
        split_sample_group != sample_column
    ]

    if (
      umap_make_split_by_sample &&
      !is.na(sample_column) &&
      length(split_sample_group) > 0
    ) {
      g <- split_sample_group[[1]]
      n_panels <- length(unique(umap_df[[sample_column]]))
      dims <- calculate_split_dimensions(
        n_panels,
        umap_split_panel_width,
        umap_split_panel_height,
        umap_max_split_columns
      )

      p <- make_faceted_split_umap(
        umap_df, g, sample_column,
        palette_style = palette_style,
        title = paste0("Split by ", sample_column),
        point_size = umap_point_size,
        max_columns = umap_max_split_columns,
        legend_point_size = umap_legend_point_size,
        legend_text_size = umap_legend_text_size,
        title_size = umap_title_size,
        axis_text_size = umap_axis_text_size,
        footer_text = footer_text,
        footer_text_size = umap_footer_text_size
      )

      save_umap_pdf(
        p,
        file.path(
          output_dirs[[paste0(prefix, "split_sample")]],
          paste0(
            file_stub, "_UMAP_", g,
            "_split_by_", sample_column,
            "_", palette_style, ".pdf"
          )
        ),
        dims[["width"]],
        dims[["height"]],
        overwrite_existing
      )
    }

    split_condition_group <- c(
      annotation_column, cluster_column, sample_column
    )
    split_condition_group <- split_condition_group[
      !is.na(split_condition_group) &
        split_condition_group != condition_column
    ]

    if (
      umap_make_split_by_condition &&
      !is.na(condition_column) &&
      length(split_condition_group) > 0
    ) {
      g <- split_condition_group[[1]]
      n_panels <- length(unique(umap_df[[condition_column]]))
      dims <- calculate_split_dimensions(
        n_panels,
        umap_split_panel_width,
        umap_split_panel_height,
        umap_max_split_columns
      )

      p <- make_faceted_split_umap(
        umap_df, g, condition_column,
        palette_style = palette_style,
        title = paste0("Split by ", condition_column),
        point_size = umap_point_size,
        max_columns = umap_max_split_columns,
        legend_point_size = umap_legend_point_size,
        legend_text_size = umap_legend_text_size,
        title_size = umap_title_size,
        axis_text_size = umap_axis_text_size,
        footer_text = footer_text,
        footer_text_size = umap_footer_text_size
      )

      save_umap_pdf(
        p,
        file.path(
          output_dirs[[paste0(prefix, "split_condition")]],
          paste0(
            file_stub, "_UMAP_", g,
            "_split_by_", condition_column,
            "_", palette_style, ".pdf"
          )
        ),
        dims[["width"]],
        dims[["height"]],
        overwrite_existing
      )
    }
  }

  if (
    umap_make_each_sample_monochrome &&
    !is.na(sample_column)
  ) {
    for (sample_value in natural_level_order(umap_df[[sample_column]])) {
      p <- make_single_sample_monochrome_umap(
        umap_df,
        sample_column,
        sample_value,
        title = paste0("Sample: ", sample_value),
        point_color = umap_each_sample_color,
        point_size = umap_each_sample_point_size,
        global_limits = umap_each_sample_show_global_limits,
        transparent_background = umap_each_sample_transparent_background,
        title_size = umap_title_size,
        axis_text_size = umap_axis_text_size,
        footer_text = footer_text,
        footer_text_size = umap_footer_text_size
      )

      save_umap_pdf(
        p,
        file.path(
          output_dirs[["each_sample"]],
          paste0(
            file_stub,
            "_UMAP_sample_",
            safe_umap_filename(sample_value),
            "_transparent.pdf"
          )
        ),
        umap_each_sample_pdf_width,
        umap_each_sample_pdf_height,
        overwrite_existing,
        transparent = umap_each_sample_transparent_background
      )
    }
  }

  rm(object, umap_df)
  gc(verbose = FALSE)
}

safe_write_csv(
  do.call(rbind, resolver_rows),
  file.path(
    output_dirs[["tables"]],
    paste0("UMAP_resolved_columns_v2.3_", run_stamp, ".csv")
  )
)

cat("\nCompleted: ", timestamp_string(), "\n", sep = "")
