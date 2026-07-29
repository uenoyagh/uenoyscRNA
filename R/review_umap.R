#' Create annotation-review UMAP plots
#'
#' Generates UMAP plots grouped by annotation, sample, and condition using
#' automatically detected review settings unless explicit columns are supplied.
#'
#' @param object A Seurat object.
#' @param annotation_column Optional annotation metadata column.
#' @param sample_column Optional sample metadata column.
#' @param condition_column Optional condition metadata column.
#' @param reduction Optional dimensional reduction name.
#' @param assay Optional assay name.
#' @param palette Optional named vector of colors for annotation values.
#' @param sample_palette Optional named vector of colors for sample values.
#' @param condition_palette Optional named vector of colors for condition values.
#' @param point_size Point size for UMAP cells.
#' @param label Logical; label annotation groups.
#' @param repel Logical; repel annotation labels.
#' @param raster Logical; rasterize UMAP points.
#' @param rds_file Optional source RDS path or filename shown in the footer.
#' @param analysis_date Analysis date shown in the footer. Defaults to `Sys.Date()`.
#' @param show_provenance Show source RDS and analysis date in each plot footer.
#' @param provenance_size Footer text size.
#' @param provenance_color Footer text color.
#' @param output_dir Optional output directory. If `NULL`, plots are returned
#'   without writing files.
#' @param width PDF width.
#' @param height PDF height.
#'
#' @return A named list of ggplot objects and written file paths.
#' @export
review_umap <- function(
    object,
    group_by = NULL,
    split_by = NULL,
    ...,
    annotation_column = NULL,
    sample_column = NULL,
    condition_column = NULL,
    reduction = NULL,
    assay = NULL,
    palette = NULL,
    sample_palette = NULL,
    condition_palette = NULL,
    point_size = 0.6,
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    rds_file = NULL,
    analysis_date = Sys.Date(),
    show_provenance = TRUE,
    provenance_size = 6,
    provenance_color = "#666666",
    output_dir = NULL,
    width = 8,
    height = 7
) {
  advanced_args <- list(...)
  use_advanced <- !is.null(group_by) || !is.null(split_by) || length(advanced_args) > 0L

  if (isTRUE(use_advanced)) {
    base_args <- list(
      object = object,
      group_by = group_by,
      split_by = split_by,
      assay = assay,
      reduction = reduction,
      pt_size = point_size,
      label = label,
      repel = repel,
      raster = raster,
      rds_file = rds_file,
      analysis_date = analysis_date,
      show_provenance = show_provenance,
      provenance_size = provenance_size,
      provenance_color = provenance_color,
      output_dir = output_dir,
      width = width,
      height = height
    )
    if (!is.null(palette) && is.null(advanced_args$colors)) {
      base_args$colors <- palette
    }
    advanced_call <- utils::modifyList(base_args, advanced_args)
    return(do.call(.review_umap_advanced, advanced_call))
  }

  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  settings <- detect_review_settings(
    object = object,
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column,
    reduction = reduction,
    assay = assay
  )

  plots <- list()

  plots$annotation <- .make_review_umap(
    object = object,
    group_column = settings$annotation_column,
    reduction = settings$reduction,
    assay = settings$assay,
    palette = palette,
    title = "UMAP by annotation",
    point_size = point_size,
    label = label,
    repel = repel,
    raster = raster
  )

  if (!is.null(settings$sample_column)) {
    plots$sample <- .make_review_umap(
      object = object,
      group_column = settings$sample_column,
      reduction = settings$reduction,
      assay = settings$assay,
      palette = sample_palette,
      title = "UMAP by sample",
      point_size = point_size,
      label = FALSE,
      repel = repel,
      raster = raster
    )
  }

  if (!is.null(settings$condition_column)) {
    plots$condition <- .make_review_umap(
      object = object,
      group_column = settings$condition_column,
      reduction = settings$reduction,
      assay = settings$assay,
      palette = condition_palette,
      title = "UMAP by condition",
      point_size = point_size,
      label = FALSE,
      repel = repel,
      raster = raster
    )
  }

  if (isTRUE(show_provenance)) {
    caption <- .review_provenance_text(
      rds_file = rds_file,
      analysis_date = analysis_date
    )
    plots <- lapply(plots, function(plot) {
      plot +
        ggplot2::labs(caption = caption) +
        ggplot2::theme(
          plot.caption = ggplot2::element_text(
            size = provenance_size,
            colour = provenance_color,
            hjust = 1,
            margin = ggplot2::margin(t = 5, unit = "pt")
          )
        )
    })
  }

  files <- character()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    file_names <- c(
      annotation = "UMAP_by_annotation.pdf",
      sample = "UMAP_by_sample.pdf",
      condition = "UMAP_by_condition.pdf"
    )

    for (name in names(plots)) {
      path <- file.path(output_dir, file_names[[name]])
      ggplot2::ggsave(
        filename = path,
        plot = plots[[name]],
        width = width,
        height = height,
        units = "in",
        device = grDevices::cairo_pdf
      )
      files[[name]] <- path
    }
  }

  structure(
    list(
      plots = plots,
      files = files,
      settings = settings
    ),
    class = c("uenoy_review_umap", "list")
  )
}


.make_review_umap <- function(
    object,
    group_by = NULL,
    split_by = NULL,
    ...,
    group_column,
    reduction,
    assay,
    palette,
    title,
    point_size,
    label,
    repel,
    raster
) {
  metadata_columns <- colnames(object[[]])

  if (!group_column %in% metadata_columns) {
    stop("Metadata column not found: ", group_column, call. = FALSE)
  }

  old_assay <- SeuratObject::DefaultAssay(object)
  on.exit(SeuratObject::DefaultAssay(object) <- old_assay, add = TRUE)
  SeuratObject::DefaultAssay(object) <- assay

  plot <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_column,
    pt.size = point_size,
    label = label,
    repel = repel,
    raster = raster
  ) +
    ggplot2::labs(
      title = title,
      color = group_column
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(0.45, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      )
    )

  if (!is.null(palette)) {
    values <- unique(as.character(object[[]][[group_column]]))
    values <- values[!is.na(values)]

    if (is.null(names(palette))) {
      if (length(palette) < length(values)) {
        stop(
          "The supplied palette has fewer colors than metadata levels.",
          call. = FALSE
        )
      }
      names(palette) <- values
    }

    missing_colors <- setdiff(values, names(palette))
    if (length(missing_colors) > 0L) {
      stop(
        "Palette is missing colors for: ",
        paste(missing_colors, collapse = ", "),
        call. = FALSE
      )
    }

    plot <- plot +
      ggplot2::scale_color_manual(
        values = palette,
        breaks = values,
        drop = FALSE
      )
  }

  plot
}
