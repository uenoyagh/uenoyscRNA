#' Create a publication-style dot plot
#'
#' Creates a consistently formatted Seurat dot plot in which dot size
#' represents the percentage of expressing cells and colour represents
#' scaled average expression.
#'
#' @param object A Seurat object.
#' @param features Character vector of feature names, or a named list of
#'   feature vectors.
#' @param assay Optional assay name.
#' @param group.by Optional metadata column used to group cells.
#' @param idents Optional identity classes to include.
#' @param low Colour used for low expression.
#' @param mid Colour used for midpoint expression.
#' @param high Colour used for high expression.
#' @param midpoint Numeric midpoint used in the colour scale.
#' @param colour_quantiles Numeric vector of length two defining lower and
#'   upper quantiles for the colour scale. Use `NULL` to disable
#'   quantile-based limits.
#' @param colour_limits Optional numeric vector of length two defining fixed
#'   colour limits. Takes precedence over `colour_quantiles`.
#' @param dot_min Minimum fraction of cells expressing a feature required
#'   for a dot to be shown.
#' @param dot_scale Maximum dot size.
#' @param scale Logical. Whether average expression should be scaled by
#'   feature.
#' @param scale_by Character string passed to [Seurat::DotPlot()].
#'   Usually `"radius"` or `"size"`.
#' @param scale_min Optional minimum dot-size scaling value.
#' @param scale_max Optional maximum dot-size scaling value.
#' @param feature_order Optional character vector defining feature order.
#' @param group_order Optional character vector defining group order.
#' @param rotate_features Logical. Whether feature labels should be rotated.
#' @param feature_angle Numeric angle for feature labels.
#' @param legend_position Legend position.
#' @param title Optional plot title.
#' @param x_title Optional x-axis title.
#' @param y_title Optional y-axis title.
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A `ggplot` object.
#'
#' @export
publish_dotplot <- function(
    object,
    features,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    colour_quantiles = c(0.05, 0.95),
    colour_limits = NULL,
    dot_min = 0,
    dot_scale = 6,
    scale = TRUE,
    scale_by = "radius",
    scale_min = NA,
    scale_max = NA,
    feature_order = NULL,
    group_order = NULL,
    rotate_features = TRUE,
    feature_angle = 45,
    legend_position = "right",
    title = NULL,
    x_title = NULL,
    y_title = NULL,
    base_size = 11,
    base_family = ""
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  valid_features <- is.character(features) || is.list(features)

  if (!valid_features || length(features) < 1L) {
    stop(
      "`features` must be a non-empty character vector or named list.",
      call. = FALSE
    )
  }

  flattened_features <- unlist(
    features,
    use.names = FALSE
  )

  if (!is.character(flattened_features) ||
      length(flattened_features) < 1L ||
      anyNA(flattened_features) ||
      any(!nzchar(flattened_features))) {
    stop(
      "`features` must contain non-empty character feature names.",
      call. = FALSE
    )
  }

  validate_optional_single_string(assay, "assay")
  validate_optional_single_string(group.by, "group.by")
  validate_optional_single_string(title, "title")
  validate_optional_single_string(x_title, "x_title")
  validate_optional_single_string(y_title, "y_title")

  validate_single_string(scale_by, "scale_by")
  validate_positive_number(dot_scale, "dot_scale")
  validate_positive_number(base_size, "base_size")

  if (!scale_by %in% c("radius", "size")) {
    stop(
      "`scale_by` must be either \"radius\" or \"size\".",
      call. = FALSE
    )
  }

  if (!is.numeric(dot_min) ||
      length(dot_min) != 1L ||
      is.na(dot_min) ||
      dot_min < 0 ||
      dot_min > 1) {
    stop(
      "`dot_min` must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(midpoint) ||
      length(midpoint) != 1L ||
      is.na(midpoint) ||
      !is.finite(midpoint)) {
    stop(
      "`midpoint` must be a single finite numeric value.",
      call. = FALSE
    )
  }

  if (!is.numeric(feature_angle) ||
      length(feature_angle) != 1L ||
      is.na(feature_angle) ||
      !is.finite(feature_angle)) {
    stop(
      "`feature_angle` must be a single finite numeric value.",
      call. = FALSE
    )
  }

  logical_arguments <- list(
    scale = scale,
    rotate_features = rotate_features
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

  if (!is.null(assay)) {
    available_assays <- names(object@assays)

    if (!assay %in% available_assays) {
      stop(
        paste0(
          "Assay `",
          assay,
          "` was not found. Available assays: ",
          paste(available_assays, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
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

  if (!is.null(feature_order)) {
    if (!is.character(feature_order) ||
        length(feature_order) < 1L ||
        anyNA(feature_order) ||
        any(!nzchar(feature_order))) {
      stop(
        "`feature_order` must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }

    missing_feature_order <- setdiff(
      feature_order,
      flattened_features
    )

    if (length(missing_feature_order) > 0L) {
      stop(
        paste0(
          "The following `feature_order` values are not present in ",
          "`features`: ",
          paste(missing_feature_order, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(group_order) &&
      (!is.character(group_order) ||
       length(group_order) < 1L ||
       anyNA(group_order) ||
       any(!nzchar(group_order)))) {
    stop(
      "`group_order` must be NULL or a non-empty character vector.",
      call. = FALSE
    )
  }

  if (!is.null(colour_limits)) {
    if (!is.numeric(colour_limits) ||
        length(colour_limits) != 2L ||
        anyNA(colour_limits) ||
        any(!is.finite(colour_limits)) ||
        colour_limits[1] >= colour_limits[2]) {
      stop(
        paste0(
          "`colour_limits` must contain two finite increasing ",
          "numeric values."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(colour_quantiles)) {
    if (!is.numeric(colour_quantiles) ||
        length(colour_quantiles) != 2L ||
        anyNA(colour_quantiles) ||
        any(!is.finite(colour_quantiles)) ||
        any(colour_quantiles < 0) ||
        any(colour_quantiles > 1) ||
        colour_quantiles[1] >= colour_quantiles[2]) {
      stop(
        paste0(
          "`colour_quantiles` must contain two increasing numbers ",
          "between 0 and 1."
        ),
        call. = FALSE
      )
    }
  }

  plot <- Seurat::DotPlot(
    object = object,
    features = features,
    assay = assay,
    dot.min = dot_min,
    dot.scale = dot_scale,
    idents = idents,
    group.by = group.by,
    scale = scale,
    scale.by = scale_by,
    scale.min = scale_min,
    scale.max = scale_max
  )

  if (!is.null(feature_order)) {
    plot$data$features.plot <- factor(
      as.character(plot$data$features.plot),
      levels = feature_order
    )
  } else {
    feature_levels <- unique(flattened_features)

    plot$data$features.plot <- factor(
      as.character(plot$data$features.plot),
      levels = feature_levels
    )
  }

  if (!is.null(group_order)) {
    observed_groups <- unique(
      as.character(plot$data$id)
    )

    missing_groups <- setdiff(
      group_order,
      observed_groups
    )

    if (length(missing_groups) > 0L) {
      stop(
        paste0(
          "The following `group_order` values were not found in the ",
          "plot data: ",
          paste(missing_groups, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    remaining_groups <- setdiff(
      observed_groups,
      group_order
    )

    plot$data$id <- factor(
      as.character(plot$data$id),
      levels = c(group_order, remaining_groups)
    )
  }

  expression_values <- plot$data$avg.exp.scaled
  expression_values <- expression_values[
    is.finite(expression_values)
  ]

  if (length(expression_values) < 1L) {
    stop(
      "No finite scaled expression values were available for plotting.",
      call. = FALSE
    )
  }

  if (!is.null(colour_limits)) {
    effective_colour_limits <- colour_limits
  } else if (!is.null(colour_quantiles)) {
    effective_colour_limits <- as.numeric(
      stats::quantile(
        expression_values,
        probs = colour_quantiles,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
    )

    if (effective_colour_limits[1] ==
        effective_colour_limits[2]) {
      effective_colour_limits <- range(
        expression_values,
        finite = TRUE
      )
    }
  } else {
    effective_colour_limits <- range(
      expression_values,
      finite = TRUE
    )
  }

  if (effective_colour_limits[1] ==
      effective_colour_limits[2]) {
    expansion <- max(
      abs(effective_colour_limits[1]) * 0.01,
      0.01
    )

    effective_colour_limits <- c(
      effective_colour_limits[1] - expansion,
      effective_colour_limits[2] + expansion
    )
  }

  plot <- suppressMessages(
    plot +
      ggplot2::scale_colour_gradient2(
        low = low,
        mid = mid,
        high = high,
        midpoint = midpoint,
        limits = effective_colour_limits,
        oob = scales::squish,
        name = "Average\nexpression"
      )
  )

  plot <- suppressMessages(
    plot +
      ggplot2::scale_size(
        range = c(0, dot_scale),
        limits = c(0, 100),
        name = "Percent\nexpressed"
      )
  )

  plot <- plot +
    theme_ueno_scRNA(
      base_size = base_size,
      base_family = base_family,
      legend_position = legend_position,
      panel_border = TRUE,
      axis_text = TRUE,
      axis_ticks = TRUE
    ) +
    ggplot2::labs(
      title = title,
      x = x_title,
      y = y_title
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey90"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 8)
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size * 1.15,
        family = base_family
      )
    )

  if (isTRUE(rotate_features)) {
    horizontal_justification <- if (feature_angle >= 0) 1 else 0

    plot <- plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = feature_angle,
          hjust = horizontal_justification,
          vjust = 1
        )
      )
  }

  plot
}
