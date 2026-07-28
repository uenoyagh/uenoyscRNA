#' Create publication-style feature plots
#'
#' Creates consistently formatted Seurat feature plots with optional
#' quantile cutoffs and a blue-white-red colour scale.
#'
#' @param object A Seurat object.
#' @param features Character vector of feature names.
#' @param reduction Name of the dimensional reduction to use.
#' @param split.by Optional metadata column used to split plots.
#' @param point_size Numeric point size for cells.
#' @param point_size_factor Numeric multiplier applied to `point_size`.
#' @param alpha Numeric point transparency between 0 and 1.
#' @param min.cutoff Lower cutoff passed to [Seurat::FeaturePlot()].
#' @param max.cutoff Upper cutoff passed to [Seurat::FeaturePlot()].
#'   Quantile strings such as `"q05"` and `"q95"` are supported.
#' @param low Colour used for low expression.
#' @param mid Colour used for midpoint expression.
#' @param high Colour used for high expression.
#' @param midpoint Numeric midpoint used in the colour scale.
#' @param keep.scale Scale handling passed to [Seurat::FeaturePlot()].
#' @param order Logical. Whether high-expression cells should be plotted last.
#' @param raster Logical. Whether to rasterise cell points.
#' @param raster.dpi Numeric vector of length two passed to
#'   [Seurat::FeaturePlot()].
#' @param ncol Optional number of columns.
#' @param combine Logical. Whether to combine panels.
#' @param legend_position Legend position.
#' @param show_axes Logical. Whether to display axes.
#' @param title Optional overall title.
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#'
#' @export
publish_feature <- function(
    object,
    features,
    reduction = "umap",
    split.by = NULL,
    point_size = 0.5,
    point_size_factor = 2,
    alpha = 1,
    min.cutoff = "q05",
    max.cutoff = "q95",
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    keep.scale = "feature",
    order = TRUE,
    raster = FALSE,
    raster.dpi = c(512, 512),
    ncol = NULL,
    combine = TRUE,
    legend_position = "right",
    show_axes = FALSE,
    title = NULL,
    base_size = 11,
    base_family = ""
) {
  .validate_seurat_object(object)

  if (!is.character(features) ||
      length(features) < 1L ||
      anyNA(features) ||
      any(!nzchar(features))) {
    stop(
      "`features` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  validate_single_string(reduction, "reduction")
  validate_optional_single_string(split.by, "split.by")
  validate_optional_single_string(title, "title")

  validate_positive_number(point_size, "point_size")
  validate_positive_number(point_size_factor, "point_size_factor")
  validate_positive_number(base_size, "base_size")

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha < 0 ||
      alpha > 1) {
    stop(
      "`alpha` must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(midpoint) ||
      length(midpoint) != 1L ||
      is.na(midpoint)) {
    stop(
      "`midpoint` must be a single numeric value.",
      call. = FALSE
    )
  }

  logical_arguments <- list(
    order = order,
    raster = raster,
    show_axes = show_axes,
    combine = combine
  )

  invalid_logical <- vapply(
    logical_arguments,
    function(x) {
      !is.logical(x) || length(x) != 1L || is.na(x)
    },
    logical(1)
  )

  if (any(invalid_logical)) {
    stop(
      paste0(
        "`",
        names(logical_arguments)[which(invalid_logical)[1]],
        "` must be TRUE or FALSE."
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(raster.dpi) ||
      length(raster.dpi) != 2L ||
      anyNA(raster.dpi) ||
      any(raster.dpi <= 0)) {
    stop(
      "`raster.dpi` must contain two positive numbers.",
      call. = FALSE
    )
  }

  available_reductions <- names(object@reductions)

  if (!reduction %in% available_reductions) {
    stop(
      paste0(
        "Reduction `",
        reduction,
        "` was not found. Available reductions: ",
        paste(available_reductions, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  .validate_metadata_column(
    object = object,
    column = split.by,
    arg = "split.by",
    allow_null = TRUE
  )

  if (isTRUE(raster) &&
      isTRUE(order) &&
      !requireNamespace("ggrastr", quietly = TRUE)) {
    warning(
      "`ggrastr` is not installed. Falling back to `raster = FALSE`.",
      call. = FALSE
    )
    raster <- FALSE
  }

  effective_point_size <- point_size * point_size_factor

  plots <- Seurat::FeaturePlot(
    object = object,
    features = features,
    reduction = reduction,
    split.by = split.by,
    pt.size = effective_point_size,
    alpha = alpha,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    keep.scale = keep.scale,
    order = order,
    raster = raster,
    raster.dpi = raster.dpi,
    ncol = ncol,
    combine = FALSE
  )

  plots <- lapply(
    plots,
    function(plot) {
      plot <- plot +
        ggplot2::scale_colour_gradient2(
          low = low,
          mid = mid,
          high = high,
          midpoint = midpoint
        ) +
        theme_ueno_scRNA(
          base_size = base_size,
          base_family = base_family,
          legend_position = legend_position,
          panel_border = TRUE,
          axis_text = show_axes,
          axis_ticks = show_axes
        )

      if (!isTRUE(show_axes)) {
        plot <- plot +
          ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
          )
      }

      plot
    }
  )

  if (!isTRUE(combine)) {
    return(plots)
  }

  combined_plot <- patchwork::wrap_plots(
    plots,
    ncol = ncol
  )

  if (!is.null(title)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold",
            size = base_size * 1.15,
            family = base_family
          )
        )
      )
  }

  combined_plot
}
