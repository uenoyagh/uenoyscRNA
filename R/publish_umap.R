#' Create a publication-style UMAP plot
#'
#' Creates a consistently formatted dimensional-reduction plot from a
#' Seurat object. The function wraps [Seurat::DimPlot()] and applies
#' [theme_ueno_scRNA()] together with enlarged legend points and optional
#' removal of UMAP axes.
#'
#' @param object A Seurat object.
#' @param reduction Name of the dimensional reduction to use.
#' @param group.by Metadata column used to colour cells. If `NULL`,
#'   the active identity is used.
#' @param split.by Optional metadata column used to split the plot.
#' @param cols Optional vector of colours. Named vectors are recommended.
#' @param label Logical. Whether to label groups on the plot.
#' @param label.size Numeric size of group labels.
#' @param repel Logical. Whether to repel group labels.
#' @param point_size Numeric point size for cells.
#' @param point_size_factor Numeric multiplier applied to `point_size`.
#'   The default of `2` reflects the standard uenoyscRNA display style.
#' @param alpha Numeric point transparency between 0 and 1.
#' @param raster Logical. Whether to rasterise cell points.
#' @param raster.dpi Numeric vector of length two passed to
#'   [Seurat::DimPlot()].
#' @param shuffle Logical. Whether to shuffle plotting order.
#' @param seed Integer random seed used when `shuffle = TRUE`.
#' @param order Optional vector specifying the plotting order of groups.
#' @param legend_position Position of the legend.
#' @param legend_point_size Numeric point size used only in the legend.
#' @param show_axes Logical. Whether to display axis titles, text, and ticks.
#' @param title Optional plot title.
#' @param ncol Optional number of columns when combining multiple panels.
#' @param combine Logical. Whether to combine panels into one patchwork object.
#' @param base_size Base font size.
#' @param base_family Base font family. The default uses the graphics-device
#'   default font for portability.
#'
#' @return A `ggplot` object, a `patchwork` object, or a list of `ggplot`
#'   objects when `combine = FALSE`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "SeuratObject")
#'
#' publish_umap(
#'   object = pbmc_small,
#'   reduction = "pca",
#'   label = TRUE
#' )
#' }
publish_umap <- function(
    object,
    reduction = "umap",
    group.by = NULL,
    split.by = NULL,
    cols = NULL,
    label = FALSE,
    label.size = 4,
    repel = TRUE,
    point_size = 0.5,
    point_size_factor = 2,
    alpha = 1,
    raster = TRUE,
    raster.dpi = c(512, 512),
    shuffle = FALSE,
    seed = 1L,
    order = NULL,
    legend_position = "right",
    legend_point_size = 4,
    show_axes = FALSE,
    title = NULL,
    ncol = NULL,
    combine = TRUE,
    base_size = 11,
    base_family = ""
) {
  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  validate_single_string(reduction, "reduction")

  validate_optional_single_string(group.by, "group.by")
  validate_optional_single_string(split.by, "split.by")
  validate_optional_single_string(title, "title")

  validate_positive_number(label.size, "label.size")
  validate_positive_number(point_size, "point_size")
  validate_positive_number(point_size_factor, "point_size_factor")
  validate_positive_number(legend_point_size, "legend_point_size")
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

  if (!is.numeric(raster.dpi) ||
      length(raster.dpi) != 2L ||
      anyNA(raster.dpi) ||
      any(raster.dpi <= 0)) {
    stop(
      "`raster.dpi` must contain two positive numbers.",
      call. = FALSE
    )
  }

  logical_arguments <- list(
    label = label,
    repel = repel,
    raster = raster,
    shuffle = shuffle,
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

  metadata_columns <- colnames(object[[]])

  if (!is.null(group.by) && !group.by %in% metadata_columns) {
    stop(
      paste0(
        "`group.by = \"",
        group.by,
        "\"` was not found in object metadata."
      ),
      call. = FALSE
    )
  }

  if (!is.null(split.by) && !split.by %in% metadata_columns) {
    stop(
      paste0(
        "`split.by = \"",
        split.by,
        "\"` was not found in object metadata."
      ),
      call. = FALSE
    )
  }

  effective_point_size <- point_size * point_size_factor

  plots <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group.by,
    split.by = split.by,
    cols = cols,
    label = label,
    label.size = label.size,
    repel = repel,
    pt.size = effective_point_size,
    alpha = alpha,
    raster = raster,
    raster.dpi = raster.dpi,
    shuffle = shuffle,
    seed = seed,
    order = order,
    combine = FALSE
  )

  plots <- lapply(
    plots,
    function(plot) {
      plot <- plot +
        theme_ueno_scRNA(
          base_size = base_size,
          base_family = base_family,
          legend_position = legend_position,
          panel_border = TRUE,
          axis_text = show_axes,
          axis_ticks = show_axes
        ) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            override.aes = list(
              size = legend_point_size,
              alpha = 1
            )
          )
        )

      if (isTRUE(show_axes)) {
        plot <- plot +
          ggplot2::labs(
            x = paste0(toupper(reduction), "_1"),
            y = paste0(toupper(reduction), "_2")
          )
      } else {
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

  if (!is.null(title)) {
    if (length(plots) == 1L) {
      plots[[1]] <- plots[[1]] +
        ggplot2::labs(title = title)
    }
  }

  if (!isTRUE(combine)) {
    return(plots)
  }

  combined_plot <- patchwork::wrap_plots(
    plots,
    ncol = ncol
  )

  if (!is.null(title) && length(plots) > 1L) {
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


validate_single_string <- function(x, argument) {
  if (!is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)) {
    stop(
      paste0(
        "`",
        argument,
        "` must be a single non-empty character string."
      ),
      call. = FALSE
    )
  }

  invisible(x)
}


validate_optional_single_string <- function(x, argument) {
  if (is.null(x)) {
    return(invisible(x))
  }

  validate_single_string(x, argument)
}


validate_positive_number <- function(x, argument) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      x <= 0) {
    stop(
      paste0(
        "`",
        argument,
        "` must be a single positive number."
      ),
      call. = FALSE
    )
  }

  invisible(x)
}
