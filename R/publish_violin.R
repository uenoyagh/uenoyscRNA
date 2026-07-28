#' Create publication-style violin plots
#'
#' Creates consistently formatted Seurat violin plots with optional
#' group ordering, split ordering, custom colours, and point rasterisation.
#'
#' @param object A Seurat object.
#' @param features Character vector of features to plot. Features may be
#'   genes, metadata columns, or other values retrievable by
#'   [SeuratObject::FetchData()].
#' @param assay Optional assay name.
#' @param layer Optional assay layer, such as `"data"` or `"counts"`.
#' @param group.by Optional metadata column used to group cells.
#' @param split.by Optional metadata column used to split violins.
#' @param idents Optional identity classes to include.
#' @param group_order Optional character vector defining group order.
#' @param split_order Optional character vector defining split-group order.
#' @param cols Optional character vector of colours passed to
#'   [Seurat::VlnPlot()].
#' @param point_size Point size for individual cells. Use `0` to hide points.
#' @param alpha Numeric point transparency between 0 and 1.
#' @param adjust Numeric bandwidth adjustment for violin density estimation.
#' @param y_max Optional numeric maximum for the y-axis.
#' @param same_y_limits Logical. Whether all panels should use the same
#'   y-axis limits.
#' @param log Logical. Whether to use a logarithmic feature axis.
#' @param sort Logical or character value controlling group sorting.
#' @param split_plot Logical. Whether split groups should share a violin
#'   shape.
#' @param stack Logical. Whether feature plots should be stacked.
#' @param fill_by Character string specifying whether violin fill is based
#'   on `"feature"` or `"ident"`.
#' @param flip Logical. Whether to flip plot orientation.
#' @param add_noise Logical. Whether Seurat should add small plotting noise.
#' @param raster Logical or `NULL`. Whether points should be rasterised.
#' @param raster_dpi Positive numeric raster resolution.
#' @param ncol Optional number of columns when combining panels.
#' @param combine Logical. Whether panels should be combined.
#' @param show_legend Logical. Whether to display the legend.
#' @param legend_position Legend position.
#' @param rotate_groups Logical. Whether group labels should be rotated.
#' @param group_angle Numeric angle for group labels.
#' @param violin_linewidth Numeric violin outline width.
#' @param title Optional overall title.
#' @param x_title Optional x-axis title.
#' @param y_title Optional y-axis title.
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A `ggplot`, `patchwork`, or list of `ggplot` objects.
#'
#' @export
publish_violin <- function(
    object,
    features,
    assay = NULL,
    layer = NULL,
    group.by = NULL,
    split.by = NULL,
    idents = NULL,
    group_order = NULL,
    split_order = NULL,
    cols = NULL,
    point_size = 0,
    alpha = 1,
    adjust = 1,
    y_max = NULL,
    same_y_limits = FALSE,
    log = FALSE,
    sort = FALSE,
    split_plot = FALSE,
    stack = FALSE,
    fill_by = "feature",
    flip = FALSE,
    add_noise = TRUE,
    raster = FALSE,
    raster_dpi = 300,
    ncol = NULL,
    combine = TRUE,
    show_legend = TRUE,
    legend_position = "right",
    rotate_groups = FALSE,
    group_angle = 45,
    violin_linewidth = 0.4,
    title = NULL,
    x_title = NULL,
    y_title = NULL,
    base_size = 11,
    base_family = ""
) {

  # ---------------------------------------------------------------------------
  # Core object and feature validation
  # ---------------------------------------------------------------------------

  .validate_seurat_object(object)

  if (
    !is.character(features) ||
    length(features) < 1L ||
    anyNA(features) ||
    any(!nzchar(features))
  ) {
    stop(
      "`features` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Character argument validation
  # ---------------------------------------------------------------------------

  validate_optional_single_string(
    assay,
    "assay"
  )

  validate_optional_single_string(
    layer,
    "layer"
  )

  validate_optional_single_string(
    group.by,
    "group.by"
  )

  validate_optional_single_string(
    split.by,
    "split.by"
  )

  validate_optional_single_string(
    title,
    "title"
  )

  validate_optional_single_string(
    x_title,
    "x_title"
  )

  validate_optional_single_string(
    y_title,
    "y_title"
  )

  validate_single_string(
    fill_by,
    "fill_by"
  )

  # ---------------------------------------------------------------------------
  # Assay and metadata validation
  # ---------------------------------------------------------------------------

  .validate_assay(
    object = object,
    assay = assay,
    allow_null = TRUE
  )

  .validate_metadata_column(
    object = object,
    column = group.by,
    arg = "group.by",
    allow_null = TRUE
  )

  if (!is.null(split.by) &&
      split.by != "ident") {

    .validate_metadata_column(
      object = object,
      column = split.by,
      arg = "split.by",
      allow_null = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Positive-number validation
  # ---------------------------------------------------------------------------

  validate_positive_number(
    adjust,
    "adjust"
  )

  validate_positive_number(
    raster_dpi,
    "raster_dpi"
  )

  validate_positive_number(
    base_size,
    "base_size"
  )

  validate_positive_number(
    violin_linewidth,
    "violin_linewidth"
  )

  # ---------------------------------------------------------------------------
  # Enumerated argument validation
  # ---------------------------------------------------------------------------

  if (!fill_by %in% c("feature", "ident")) {
    stop(
      "`fill_by` must be either \"feature\" or \"ident\".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Numeric argument validation
  # ---------------------------------------------------------------------------

  if (!is.numeric(point_size) ||
      length(point_size) != 1L ||
      is.na(point_size) ||
      !is.finite(point_size) ||
      point_size < 0) {
    stop(
      "`point_size` must be a single non-negative number.",
      call. = FALSE
    )
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      !is.finite(alpha) ||
      alpha < 0 ||
      alpha > 1) {
    stop(
      "`alpha` must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.null(y_max) &&
      (!is.numeric(y_max) ||
       length(y_max) != 1L ||
       is.na(y_max) ||
       !is.finite(y_max))) {
    stop(
      "`y_max` must be NULL or a single finite numeric value.",
      call. = FALSE
    )
  }

  if (!is.numeric(group_angle) ||
      length(group_angle) != 1L ||
      is.na(group_angle) ||
      !is.finite(group_angle)) {
    stop(
      "`group_angle` must be a single finite numeric value.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Logical argument validation
  # ---------------------------------------------------------------------------

  logical_arguments <- list(
    same_y_limits = same_y_limits,
    log = log,
    split_plot = split_plot,
    stack = stack,
    flip = flip,
    add_noise = add_noise,
    combine = combine,
    show_legend = show_legend,
    rotate_groups = rotate_groups
  )

  invalid_logical <- vapply(
    logical_arguments,
    function(x) {
      !is.logical(x) ||
        length(x) != 1L ||
        is.na(x)
    },
    logical(1)
  )

  if (any(invalid_logical)) {
    invalid_name <- names(
      logical_arguments
    )[which(invalid_logical)[1]]

    stop(
      paste0(
        "`",
        invalid_name,
        "` must be TRUE or FALSE."
      ),
      call. = FALSE
    )
  }

  if (!is.null(raster) &&
      (!is.logical(raster) ||
       length(raster) != 1L ||
       is.na(raster))) {
    stop(
      "`raster` must be TRUE, FALSE, or NULL.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Sort validation
  # ---------------------------------------------------------------------------

  valid_sort <- (
    is.logical(sort) &&
      length(sort) == 1L &&
      !is.na(sort)
  ) || (
    is.character(sort) &&
      length(sort) == 1L &&
      !is.na(sort) &&
      sort %in% c(
        "increasing",
        "decreasing"
      )
  )

  if (!valid_sort) {
    stop(
      paste0(
        "`sort` must be TRUE, FALSE, ",
        "\"increasing\", or \"decreasing\"."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Identity and colour validation
  # ---------------------------------------------------------------------------

  if (!is.null(idents) &&
      (!is.character(idents) &&
       !is.factor(idents) &&
       !is.numeric(idents))) {
    stop(
      paste0(
        "`idents` must be NULL, or a character, ",
        "factor, or numeric vector."
      ),
      call. = FALSE
    )
  }

  if (!is.null(cols) &&
      (!is.character(cols) ||
       length(cols) < 1L ||
       anyNA(cols) ||
       any(!nzchar(cols)))) {
    stop(
      "`cols` must be NULL or a non-empty character vector.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Group-order argument validation
  # ---------------------------------------------------------------------------

  if (!is.null(group_order)) {

    .validate_character_vector(
      x = group_order,
      arg = "group_order"
    )

    if (any(!nzchar(group_order))) {
      stop(
        paste0(
          "`group_order` must be NULL or a non-empty ",
          "character vector."
        ),
        call. = FALSE
      )
    }

    if (anyDuplicated(group_order)) {
      stop(
        "`group_order` must not contain duplicated values.",
        call. = FALSE
      )
    }
  }

  if (!is.null(split_order)) {

    .validate_character_vector(
      x = split_order,
      arg = "split_order"
    )

    if (any(!nzchar(split_order))) {
      stop(
        paste0(
          "`split_order` must be NULL or a non-empty ",
          "character vector."
        ),
        call. = FALSE
      )
    }

    if (anyDuplicated(split_order)) {
      stop(
        "`split_order` must not contain duplicated values.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Layout validation
  # ---------------------------------------------------------------------------

  if (!is.null(ncol) &&
      (!is.numeric(ncol) ||
       length(ncol) != 1L ||
       is.na(ncol) ||
       !is.finite(ncol) ||
       ncol < 1 ||
       ncol != as.integer(ncol))) {
    stop(
      "`ncol` must be NULL or a positive integer.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Apply group ordering
  #
  # This logic is retained explicitly because group.by = NULL uses active
  # identities rather than a metadata column.
  # ---------------------------------------------------------------------------

  if (!is.null(group_order)) {

    if (is.null(group.by)) {
      observed_groups <- unique(
        as.character(
          SeuratObject::Idents(object)
        )
      )
    } else {
      observed_groups <- unique(
        as.character(
          object[[group.by]][, 1]
        )
      )
    }

    observed_groups <- observed_groups[
      !is.na(observed_groups)
    ]

    missing_groups <- setdiff(
      group_order,
      observed_groups
    )

    if (length(missing_groups) > 0L) {
      stop(
        paste0(
          "The following `group_order` values were not found: ",
          paste(
            missing_groups,
            collapse = ", "
          ),
          "."
        ),
        call. = FALSE
      )
    }

    if (is.null(group.by)) {

      SeuratObject::Idents(object) <-
        .apply_group_order(
          x = SeuratObject::Idents(object),
          order = group_order
        )

    } else {

      object[[group.by]] <-
        .apply_group_order(
          x = object[[group.by]][, 1],
          order = group_order
        )
    }
  }

  # ---------------------------------------------------------------------------
  # Apply split ordering
  # ---------------------------------------------------------------------------

  if (!is.null(split_order)) {

    if (is.null(split.by)) {
      stop(
        paste0(
          "`split_order` can only be used when ",
          "`split.by` is specified."
        ),
        call. = FALSE
      )
    }

    if (split.by == "ident") {
      observed_splits <- unique(
        as.character(
          SeuratObject::Idents(object)
        )
      )
    } else {
      observed_splits <- unique(
        as.character(
          object[[split.by]][, 1]
        )
      )
    }

    observed_splits <- observed_splits[
      !is.na(observed_splits)
    ]

    missing_splits <- setdiff(
      split_order,
      observed_splits
    )

    if (length(missing_splits) > 0L) {
      stop(
        paste0(
          "The following `split_order` values were not found: ",
          paste(
            missing_splits,
            collapse = ", "
          ),
          "."
        ),
        call. = FALSE
      )
    }

    if (split.by == "ident") {

      SeuratObject::Idents(object) <-
        .apply_group_order(
          x = SeuratObject::Idents(object),
          order = split_order
        )

    } else {

      object[[split.by]] <-
        .apply_group_order(
          x = object[[split.by]][, 1],
          order = split_order
        )
    }
  }

  # ---------------------------------------------------------------------------
  # Rasterisation fallback
  # ---------------------------------------------------------------------------

  if (isTRUE(raster) &&
      point_size > 0 &&
      !requireNamespace(
        "ggrastr",
        quietly = TRUE
      )) {
    warning(
      paste0(
        "`ggrastr` is not installed. ",
        "Falling back to `raster = FALSE`."
      ),
      call. = FALSE
    )

    raster <- FALSE
  }

  # ---------------------------------------------------------------------------
  # Create violin plots
  # ---------------------------------------------------------------------------

  plots <- Seurat::VlnPlot(
    object = object,
    features = features,
    cols = cols,
    pt.size = point_size,
    alpha = alpha,
    idents = idents,
    sort = sort,
    assay = assay,
    group.by = group.by,
    split.by = split.by,
    adjust = adjust,
    y.max = y_max,
    same.y.lims = same_y_limits,
    log = log,
    ncol = ncol,
    layer = layer,
    split.plot = split_plot,
    stack = stack,
    combine = FALSE,
    fill.by = fill_by,
    flip = flip,
    add.noise = add_noise,
    raster = raster,
    raster.dpi = raster_dpi
  )

  if (!is.list(plots)) {
    plots <- list(plots)
  }

  # ---------------------------------------------------------------------------
  # Apply publication formatting
  # ---------------------------------------------------------------------------

  plots <- lapply(
    plots,
    function(plot) {

      plot <- plot +
        theme_ueno_scRNA(
          base_size = base_size,
          base_family = base_family,
          legend_position = if (isTRUE(show_legend)) {
            legend_position
          } else {
            "none"
          },
          panel_border = TRUE,
          axis_text = TRUE,
          axis_ticks = TRUE
        ) +
        ggplot2::labs(
          x = x_title,
          y = y_title
        ) +
        ggplot2::theme(
          panel.grid.major.x =
            ggplot2::element_blank(),
          panel.grid.minor =
            ggplot2::element_blank(),
          axis.title.x =
            ggplot2::element_text(
              margin = ggplot2::margin(t = 8)
            ),
          axis.title.y =
            ggplot2::element_text(
              margin = ggplot2::margin(r = 8)
            ),
          plot.title =
            ggplot2::element_text(
              hjust = 0.5,
              face = "bold",
              size = base_size * 1.05,
              family = base_family
            ),
          legend.title =
            ggplot2::element_text(
              face = "bold"
            )
        )

      violin_layers <- vapply(
        plot$layers,
        function(layer_object) {
          inherits(
            layer_object$geom,
            "GeomViolin"
          )
        },
        logical(1)
      )

      if (any(violin_layers)) {
        for (layer_index in which(violin_layers)) {
          plot$layers[[layer_index]]$aes_params$linewidth <-
            violin_linewidth
        }
      }

      if (isTRUE(rotate_groups) &&
          !isTRUE(flip)) {

        horizontal_justification <- if (
          group_angle >= 0
        ) {
          1
        } else {
          0
        }

        plot <- plot +
          ggplot2::theme(
            axis.text.x =
              ggplot2::element_text(
                angle = group_angle,
                hjust = horizontal_justification,
                vjust = 1
              )
          )
      }

      plot
    }
  )

  # ---------------------------------------------------------------------------
  # Return separate plots when combine = FALSE
  # ---------------------------------------------------------------------------

  if (!isTRUE(combine)) {
    return(plots)
  }

  # ---------------------------------------------------------------------------
  # Combine plots
  # ---------------------------------------------------------------------------

  combined_plot <- patchwork::wrap_plots(
    plots,
    ncol = ncol
  )

  if (!is.null(title)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title =
            ggplot2::element_text(
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
