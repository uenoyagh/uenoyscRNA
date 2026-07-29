# ============================================================
# 05_signature_pipeline.R
# Unified Heatmap / FeaturePlot / Violin / DotPlot pipeline
# uenoy scRNAseq Framework v2.5
# ============================================================

rm(list = ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"

source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "config", "signature_registry.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "utils.R"))
source(file.path(script_root, "R", "umap.R"))
source(file.path(script_root, "R", "signature_engine.R"))

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "patchwork",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

unknown_sets <- setdiff(
  signature_sets_to_run,
  names(signature_registry)
)

if (length(unknown_sets) > 0) {
  stop(
    "Unknown signature set(s): ",
    paste(unknown_sets, collapse = ", ")
  )
}

target_dir <- get_dataset_dir(analysis_target)

output_root <- file.path(
  get_result_dir(
    analysis_target,
    "module_score",
    create = TRUE
  ),
  "signature_pipeline_v2.5"
)

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

rds_files <- list_rds_files(
  target_dir,
  recursive = FALSE
)

if (!is.null(signature_selected_files)) {
  missing_selected <- setdiff(
    signature_selected_files,
    basename(rds_files)
  )

  if (length(missing_selected) > 0) {
    stop(
      "Selected RDS file(s) not found:\n",
      paste(missing_selected, collapse = "\n")
    )
  }

  rds_files <- rds_files[
    basename(rds_files) %in% signature_selected_files
  ]
}

for (rds_index in seq_along(rds_files)) {
  rds_path <- rds_files[[rds_index]]
  rds_name <- basename(rds_path)
  file_stub <- safe_umap_filename(rds_name)

  cat(
    "[", rds_index, "/", length(rds_files), "] ",
    rds_name, "\n",
    sep = ""
  )

  object <- safe_read_rds(rds_path)

  if (inherits(object, "rds_read_error")) {
    warning(
      "Failed to read: ",
      rds_name,
      "\n",
      object$message
    )
    next
  }

  if (!inherits(object, "Seurat")) {
    warning("Not a Seurat object: ", rds_name)
    next
  }

  # ------------------------------------------------------------
  # Resolve assay
  # ------------------------------------------------------------

  assay <- if (is.null(signature_assay_override)) {
    SeuratObject::DefaultAssay(object)
  } else {
    signature_assay_override
  }

  # ------------------------------------------------------------
  # Join split Seurat v5 layers for downstream plotting
  # ------------------------------------------------------------

  available_layers <- SeuratObject::Layers(
    object[[assay]]
  )

  split_data_layers <- available_layers[
    grepl("^data\\.", available_layers)
  ]

  split_count_layers <- available_layers[
    grepl("^counts\\.", available_layers)
  ]

  if (
    length(split_data_layers) > 0 ||
    length(split_count_layers) > 0
  ) {
    message(
      "Joining Seurat v5 layers for assay '",
      assay,
      "': ",
      paste(
        c(split_count_layers, split_data_layers),
        collapse = ", "
      )
    )

    object <- SeuratObject::JoinLayers(
      object = object,
      assay = assay
    )
  }


  reduction <- resolve_umap_reduction(
    object,
    signature_reduction_override
  )

  sample_column <- resolve_sample_column(object)
  condition_column <- resolve_condition_column(object)

  annotation_column <- resolve_annotation_column(
    object,
    signature_annotation_column_override
  )

  cluster_column <- resolve_cluster_column(
    object,
    signature_cluster_column_override
  )

  group_column <- resolve_signature_group_column(
    object = object,
    override = signature_group_column_override,
    condition_column = condition_column,
    sample_column = sample_column,
    annotation_column = annotation_column,
    cluster_column = cluster_column
  )

  footer_text <- build_signature_footer(
    rds_name,
    project_config$framework_name,
    project_config$framework_version
  )

  for (signature_set_name in signature_sets_to_run) {
    signatures <- signature_registry[[signature_set_name]]

    set_root <- file.path(
      output_root,
      signature_set_name,
      file_stub
    )

    dirs <- file.path(
      set_root,
      c(
        "01_heatmap",
        "02_featureplot",
        "03_violin",
        "04_dotplot",
        "05_tables"
      )
    )

    for (d in dirs) {
      dir.create(
        d,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    score_result <- calculate_cell_signature_scores(
      object = object,
      signatures = signatures,
      assay = assay,
      min_genes = signature_min_genes
    )

    score_df <- score_result$scores

    safe_write_csv(
      cbind(
        data.frame(
          cell = rownames(score_df),
          stringsAsFactors = FALSE
        ),
        score_df
      ),
      file.path(
        dirs[[5]],
        paste0(file_stub, "_cell_signature_scores.csv")
      )
    )

    safe_write_csv(
      score_result$gene_report,
      file.path(
        dirs[[5]],
        paste0(file_stub, "_signature_gene_report.csv")
      )
    )

    # Heatmap uses the resolved cluster column.
    if (
      isTRUE(signature_make_heatmap) &&
      !is.na(cluster_column)
    ) {
      summary_df <- calculate_cluster_signature_summary(
        score_df = score_df,
        group_values = object[[cluster_column, drop = TRUE]],
        group_name = cluster_column
      )

      safe_write_csv(
        summary_df,
        file.path(
          dirs[[5]],
          paste0(file_stub, "_cluster_signature_summary.csv")
        )
      )

      p <- make_signature_heatmap_plot(
        summary_df = summary_df,
        group_column = cluster_column,
        clip = signature_heatmap_clip,
        low_color = signature_low_color,
        mid_color = signature_mid_color,
        high_color = signature_high_color,
        title = paste0(
          signature_set_name,
          ": functional signatures by cluster"
        )
      )

      p <- add_signature_footer(
        p,
        footer_text,
        signature_footer_text_size
      )

      save_plot_dual(
        p,
        file.path(
          dirs[[1]],
          paste0(
            file_stub,
            "_",
            signature_set_name,
            "_heatmap"
          )
        ),
        width = signature_heatmap_pdf_width,
        height = signature_heatmap_pdf_height,
        export_pdf = signature_export_pdf,
        export_png = signature_export_png,
        png_dpi = signature_png_dpi,
        png_background = signature_png_background,
        overwrite = overwrite_existing
      )
    }

    if (isTRUE(signature_make_featureplot)) {
      p <- make_signature_featureplot(
        object = object,
        reduction = reduction,
        score_df = score_df,
        low_color = signature_low_color,
        mid_color = signature_mid_color,
        high_color = signature_high_color,
        q_low = signature_feature_quantile_low,
        q_high = signature_feature_quantile_high,
        max_columns = signature_featureplot_max_columns
      )

      p <- add_signature_footer(
        p,
        footer_text,
        signature_footer_text_size
      )

      n_signatures <- ncol(score_df)
      n_columns <- min(
        signature_featureplot_max_columns,
        n_signatures
      )
      n_rows <- ceiling(n_signatures / n_columns)

      save_plot_dual(
        p,
        file.path(
          dirs[[2]],
          paste0(
            file_stub,
            "_",
            signature_set_name,
            "_featureplot"
          )
        ),
        width = signature_featureplot_panel_width * n_columns,
        height = signature_featureplot_panel_height * n_rows,
        export_pdf = signature_export_pdf,
        export_png = signature_export_png,
        png_dpi = signature_png_dpi,
        png_background = signature_png_background,
        overwrite = overwrite_existing
      )
    }

    if (isTRUE(signature_make_violin)) {
      p <- make_signature_violin_plot(
        score_df = score_df,
        group_values = object[[group_column, drop = TRUE]],
        show_points = signature_violin_show_points
      )

      p <- add_signature_footer(
        p,
        footer_text,
        signature_footer_text_size
      )

      save_plot_dual(
        p,
        file.path(
          dirs[[3]],
          paste0(
            file_stub,
            "_",
            signature_set_name,
            "_violin_by_",
            group_column
          )
        ),
        width = signature_violin_width,
        height = signature_violin_height,
        export_pdf = signature_export_pdf,
        export_png = signature_export_png,
        png_dpi = signature_png_dpi,
        png_background = signature_png_background,
        overwrite = overwrite_existing
      )
    }

    if (
      isTRUE(signature_make_dotplot) &&
      !is.na(cluster_column)
    ) {
      p <- make_signature_gene_dotplot(
        object = object,
        signatures = signatures,
        group_column = cluster_column,
        assay = assay,
        max_genes_per_signature =
          signature_dotplot_max_genes_per_signature
      )

      p <- add_signature_footer(
        p,
        footer_text,
        signature_footer_text_size
      )

      save_plot_dual(
        p,
        file.path(
          dirs[[4]],
          paste0(
            file_stub,
            "_",
            signature_set_name,
            "_dotplot_by_",
            cluster_column
          )
        ),
        width = signature_dotplot_width,
        height = signature_dotplot_height,
        export_pdf = signature_export_pdf,
        export_png = signature_export_png,
        png_dpi = signature_png_dpi,
        png_background = signature_png_background,
        overwrite = overwrite_existing
      )
    }

    rm(score_result, score_df)
    gc(verbose = FALSE)
  }

  rm(object)
  gc(verbose = FALSE)
}

cat("\nSignature pipeline completed.\n")
cat("Output root:\n", output_root, "\n")
