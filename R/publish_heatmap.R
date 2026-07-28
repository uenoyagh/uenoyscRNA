#' Create a publication-ready Seurat heatmap
#'
#' A publication-oriented wrapper around [Seurat::DoHeatmap()] with
#' input validation, optional group ordering, consistent formatting,
#' and support for Seurat assay layers.
#'
#' @param object A Seurat object.
#' @param features Character vector of genes or features to display.
#' @param assay Assay to use. If `NULL`, the default assay is used.
#' @param slot Data layer used by `Seurat::DoHeatmap()`. Usually
#'   `"scale.data"`, `"data"`, or `"counts"`.
#' @param cells Optional character vector of cell names to include.
#' @param group.by Metadata column used to group cells. Use `"ident"`
#'   for the active identity class.
#' @param group_order Optional order of groups.
#' @param group.bar Logical; display the coloured group bar.
#' @param group.colors Optional character vector of group colours.
#' @param group.bar.height Height of the group bar.
#' @param disp.min Minimum displayed expression value.
#' @param disp.max Maximum displayed expression value. If `NULL`,
#'   Seurat selects the default value.
#' @param raster Logical; rasterise the heatmap body.
#' @param draw.lines Logical; draw separating lines between groups.
#' @param lines.width Width of group-separating lines.
#' @param label Logical; display group labels.
#' @param label_size Size of group labels.
#' @param label_angle Angle of group labels.
#' @param label_hjust Horizontal justification of group labels.
#' @param feature_labels Logical; display feature names.
#' @param reverse_features Logical; reverse the feature order.
#' @param show_legend Logical; display the legend.
#' @param base_size Base font size.
#' @param title Optional plot title.
#' @param combine Logical; combine plots when multiple grouping variables
#'   are supplied.
#'
#' @return A ggplot or patchwork object. If `combine = FALSE`, a list of
#'   ggplot objects may be returned.
#'
#' @export
publish_heatmap <- function(
    object,
    features,
    assay = NULL,
    slot = "scale.data",
    cells = NULL,
    group.by = "ident",
    group_order = NULL,
    group.bar = TRUE,
    group.colors = NULL,
    group.bar.height = 0.02,
    disp.min = -2.5,
    disp.max = NULL,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    label = TRUE,
    label_size = 5,
    label_angle = 45,
    label_hjust = 0,
    feature_labels = TRUE,
    reverse_features = FALSE,
    show_legend = TRUE,
    base_size = 11,
    title = NULL,
    combine = TRUE
) {

  # ---------------------------------------------------------------------------
  # Validate object
  # ---------------------------------------------------------------------------

  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Validate features
  # ---------------------------------------------------------------------------

  if (
    !is.character(features) ||
    length(features) == 0L ||
    anyNA(features) ||
    any(!nzchar(features))
  ) {
    stop(
      "`features` must be a non-empty character vector without missing values.",
      call. = FALSE
    )
  }

  features <- unique(features)

  # ---------------------------------------------------------------------------
  # Validate assay
  # ---------------------------------------------------------------------------

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  if (
    !is.character(assay) ||
    length(assay) != 1L ||
    is.na(assay) ||
    !nzchar(assay)
  ) {
    stop(
      "`assay` must be NULL or a single assay name.",
      call. = FALSE
    )
  }

  available_assays <- SeuratObject::Assays(object)

  if (!assay %in% available_assays) {
    stop(
      sprintf(
        "Assay `%s` was not found. Available assays: %s.",
        assay,
        paste(available_assays, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  object <- object
  SeuratObject::DefaultAssay(object) <- assay

  # ---------------------------------------------------------------------------
  # Validate slot
  # ---------------------------------------------------------------------------

  allowed_slots <- c(
    "scale.data",
    "data",
    "counts"
  )

  if (
    !is.character(slot) ||
    length(slot) != 1L ||
    is.na(slot) ||
    !slot %in% allowed_slots
  ) {
    stop(
      sprintf(
        "`slot` must be one of: %s.",
        paste(allowed_slots, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Check feature availability
  # ---------------------------------------------------------------------------

  assay_features <- rownames(object[[assay]])

  missing_from_assay <- setdiff(
    features,
    assay_features
  )

  if (length(missing_from_assay) > 0L) {
    stop(
      sprintf(
        paste0(
          "The following features were not found in assay `%s`: %s."
        ),
        assay,
        paste(missing_from_assay, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  layer_features <- tryCatch(
    {
      layer_data <- SeuratObject::LayerData(
        object = object[[assay]],
        layer = slot
      )

      rownames(layer_data)
    },
    error = function(e) {
      NULL
    }
  )

  if (!is.null(layer_features)) {
    missing_from_layer <- setdiff(
      features,
      layer_features
    )

    if (length(missing_from_layer) > 0L) {
      stop(
        sprintf(
          paste0(
            "The following features were not found in `%s` for assay `%s`: ",
            "%s. Run ScaleData() or select another slot."
          ),
          slot,
          assay,
          paste(missing_from_layer, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Validate cells
  # ---------------------------------------------------------------------------

  if (!is.null(cells)) {
    if (
      !is.character(cells) ||
      length(cells) == 0L ||
      anyNA(cells) ||
      any(!nzchar(cells))
    ) {
      stop(
        "`cells` must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }

    missing_cells <- setdiff(
      cells,
      colnames(object)
    )

    if (length(missing_cells) > 0L) {
      stop(
        sprintf(
          "The following cells were not found: %s.",
          paste(missing_cells, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    cells <- unique(cells)
  }

  # ---------------------------------------------------------------------------
  # Validate group.by
  # ---------------------------------------------------------------------------

  if (
    !is.character(group.by) ||
    length(group.by) == 0L ||
    anyNA(group.by) ||
    any(!nzchar(group.by))
  ) {
    stop(
      "`group.by` must be a non-empty character vector.",
      call. = FALSE
    )
  }

  metadata_columns <- colnames(object[[]])

  missing_group_columns <- setdiff(
    group.by,
    c("ident", metadata_columns)
  )

  if (length(missing_group_columns) > 0L) {
    stop(
      sprintf(
        "The following `group.by` columns were not found: %s.",
        paste(missing_group_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Apply group order
  # ---------------------------------------------------------------------------

  if (!is.null(group_order)) {
    if (
      !is.character(group_order) ||
      length(group_order) == 0L ||
      anyNA(group_order) ||
      any(!nzchar(group_order)) ||
      anyDuplicated(group_order)
    ) {
      stop(
        "`group_order` must be a unique, non-empty character vector.",
        call. = FALSE
      )
    }

    if (length(group.by) != 1L) {
      stop(
        paste0(
          "`group_order` can only be used when exactly one ",
          "`group.by` variable is supplied."
        ),
        call. = FALSE
      )
    }

    if (identical(group.by, "ident")) {
      observed_groups <- unique(
        as.character(SeuratObject::Idents(object))
      )

      missing_groups <- setdiff(
        group_order,
        observed_groups
      )

      if (length(missing_groups) > 0L) {
        stop(
          sprintf(
            "The following groups were not found in active identities: %s.",
            paste(missing_groups, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      remaining_groups <- setdiff(
        observed_groups,
        group_order
      )

      full_group_order <- c(
        group_order,
        remaining_groups
      )

      SeuratObject::Idents(object) <- factor(
        as.character(SeuratObject::Idents(object)),
        levels = full_group_order
      )
    } else {
      observed_groups <- unique(
        as.character(object[[group.by, drop = TRUE]])
      )

      missing_groups <- setdiff(
        group_order,
        observed_groups
      )

      if (length(missing_groups) > 0L) {
        stop(
          sprintf(
            "The following groups were not found in `%s`: %s.",
            group.by,
            paste(missing_groups, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      remaining_groups <- setdiff(
        observed_groups,
        group_order
      )

      full_group_order <- c(
        group_order,
        remaining_groups
      )

      object[[group.by]] <- factor(
        as.character(object[[group.by, drop = TRUE]]),
        levels = full_group_order
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Validate colours
  # ---------------------------------------------------------------------------

  if (!is.null(group.colors)) {
    if (
      !is.character(group.colors) ||
      length(group.colors) == 0L ||
      anyNA(group.colors) ||
      any(!nzchar(group.colors))
    ) {
      stop(
        "`group.colors` must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Validate numeric arguments
  # ---------------------------------------------------------------------------

  validate_single_number <- function(
    value,
    name,
    lower = -Inf,
    strictly_greater = FALSE,
    allow_null = FALSE
  ) {
    if (allow_null && is.null(value)) {
      return(invisible(TRUE))
    }

    valid <- is.numeric(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      is.finite(value)

    if (strictly_greater) {
      valid <- valid && value > lower
    } else {
      valid <- valid && value >= lower
    }

    if (!valid) {
      comparator <- if (strictly_greater) "greater than" else "at least"

      stop(
        sprintf(
          "`%s` must be a single finite number %s %s.",
          name,
          comparator,
          lower
        ),
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  validate_single_number(
    group.bar.height,
    "group.bar.height",
    lower = 0,
    strictly_greater = TRUE
  )

  validate_single_number(
    disp.min,
    "disp.min"
  )

  validate_single_number(
    disp.max,
    "disp.max",
    allow_null = TRUE
  )

  validate_single_number(
    label_size,
    "label_size",
    lower = 0,
    strictly_greater = TRUE
  )

  validate_single_number(
    label_angle,
    "label_angle"
  )

  validate_single_number(
    label_hjust,
    "label_hjust"
  )

  validate_single_number(
    base_size,
    "base_size",
    lower = 0,
    strictly_greater = TRUE
  )

  if (!is.null(lines.width)) {
    validate_single_number(
      lines.width,
      "lines.width",
      lower = 0
    )
  }

  if (!is.null(disp.max) && disp.max <= disp.min) {
    stop(
      "`disp.max` must be greater than `disp.min`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Validate logical arguments
  # ---------------------------------------------------------------------------

  logical_arguments <- list(
    group.bar = group.bar,
    raster = raster,
    draw.lines = draw.lines,
    label = label,
    feature_labels = feature_labels,
    reverse_features = reverse_features,
    show_legend = show_legend,
    combine = combine
  )

  for (argument_name in names(logical_arguments)) {
    argument_value <- logical_arguments[[argument_name]]

    if (
      !is.logical(argument_value) ||
      length(argument_value) != 1L ||
      is.na(argument_value)
    ) {
      stop(
        sprintf(
          "`%s` must be TRUE or FALSE.",
          argument_name
        ),
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Validate title
  # ---------------------------------------------------------------------------

  if (
    !is.null(title) &&
    (
      !is.character(title) ||
      length(title) != 1L ||
      is.na(title)
    )
  ) {
    stop(
      "`title` must be NULL or a single character string.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Feature ordering
  # ---------------------------------------------------------------------------

  if (isTRUE(reverse_features)) {
    features <- rev(features)
  }

  # ---------------------------------------------------------------------------
  # Build heatmap
  # ---------------------------------------------------------------------------

  heatmap <- suppressMessages(
    Seurat::DoHeatmap(
      object = object,
      features = features,
      cells = cells,
      group.by = group.by,
      group.bar = group.bar,
      group.colors = group.colors,
      disp.min = disp.min,
      disp.max = disp.max,
      slot = slot,
      assay = assay,
      label = label,
      size = label_size,
      hjust = label_hjust,
      angle = label_angle,
      raster = raster,
      draw.lines = draw.lines,
      lines.width = lines.width,
      group.bar.height = group.bar.height,
      combine = combine
    )
  )

  # ---------------------------------------------------------------------------
  # Apply publication formatting
  # ---------------------------------------------------------------------------

  format_heatmap <- function(plot_object) {
    plot_object <- plot_object +
      ggplot2::theme(
        text = ggplot2::element_text(
          size = base_size
        ),
        axis.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.y = if (isTRUE(feature_labels)) {
          ggplot2::element_text(
            size = base_size,
            colour = "black",
            face = "plain"
          )
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          size = base_size + 2,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = ggplot2::element_text(
          size = base_size
        ),
        legend.text = ggplot2::element_text(
          size = base_size - 1
        ),
        plot.margin = ggplot2::margin(
          t = 5.5,
          r = 5.5,
          b = 5.5,
          l = 5.5
        )
      )

    if (!isTRUE(show_legend)) {
      plot_object <- plot_object +
        ggplot2::theme(
          legend.position = "none"
        )
    }

    if (!is.null(title)) {
      plot_object <- plot_object +
        ggplot2::labs(
          title = title
        )
    }

    plot_object
  }

  if (is.list(heatmap) && !inherits(heatmap, "ggplot")) {
    heatmap <- lapply(
      heatmap,
      format_heatmap
    )

    if (isTRUE(combine)) {
      heatmap <- patchwork::wrap_plots(
        heatmap
      )
    }
  } else {
    heatmap <- format_heatmap(
      heatmap
    )
  }

  heatmap
}
