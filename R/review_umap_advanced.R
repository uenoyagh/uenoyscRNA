#' Review UMAP plots
#'
#' Generate one or more UMAP plots from a Seurat object.
#'
#' When `split_by` is specified, each split level is plotted separately.
#' All panels use one shared color mapping. Panel legends are suppressed and
#' a dedicated legend-only plot is added at the right side, so the final figure
#' always contains exactly one common legend even when some split groups do not
#' contain every cell type.
#'
#' @param object A Seurat object.
#' @param group_by Metadata column used to color cells. If `NULL`, active
#'   identities are used.
#' @param split_by Metadata column used to create separate panels.
#' @param group_order Optional order of grouping levels.
#' @param split_order Optional order of split panels.
#' @param assay Optional assay name.
#' @param reduction Dimensional reduction passed to `Seurat::DimPlot()`.
#' @param cells Optional vector of cell names.
#' @param pt_size Point size.
#' @param alpha Point transparency.
#' @param shuffle Randomize plotting order.
#' @param seed Random seed.
#' @param label Add group labels.
#' @param label_size Label text size.
#' @param repel Repel labels.
#' @param label_box Draw labels inside boxes.
#' @param raster Rasterize points.
#' @param colors Optional color vector. A named vector is recommended.
#' @param cells_highlight Optional cells to highlight.
#' @param cols_highlight Highlight color or colors.
#' @param cols_background Background color when highlighting cells.
#' @param legend_position Legend position. Usually `"right"`.
#' @param legend_ncol Number of legend columns.
#' @param legend_dot_size Legend point size.
#' @param legend_text_size Legend text size.
#' @param legend_title_size Legend title size.
#' @param legend_key_size Legend key size in cm.
#' @param legend_title Optional legend title.
#' @param collect_guides Retained for backward compatibility. When `TRUE`, one
#'   dedicated common legend is shown. When `FALSE`, each panel keeps its own
#'   legend.
#' @param legend_width Relative width assigned to the dedicated legend column.
#' @param panel_expand Multiplicative expansion of x/y limits.
#' @param panel_margin Plot margin in points.
#' @param panel_title_size Panel title size.
#' @param base_size Base theme size.
#' @param fixed_aspect Preserve UMAP aspect ratio.
#' @param ncol Number of panel columns.
#' @param title Optional overall title.
#' @param title_size Overall title size.
#' @param rds_file Optional source RDS file path or filename. Only the basename
#'   is printed in the figure footer.
#' @param analysis_date Analysis date shown in the figure footer. Defaults to
#'   `Sys.Date()`.
#' @param show_provenance Show the RDS filename and analysis date in the footer.
#' @param provenance_size Footer text size.
#' @param provenance_color Footer text color.
#' @param output_dir Optional PDF output directory.
#' @param filename PDF filename.
#' @param width PDF width in inches.
#' @param height PDF height in inches.
#'
#' @return A `review_umap_result` object.
.review_umap_advanced <- function(
    object,
    group_by = NULL,
    split_by = NULL,
    group_order = NULL,
    split_order = NULL,
    assay = NULL,
    reduction = NULL,
    cells = NULL,
    pt_size = NULL,
    alpha = 1,
    shuffle = FALSE,
    seed = 1L,
    label = TRUE,
    label_size = 4,
    repel = TRUE,
    label_box = FALSE,
    raster = NULL,
    colors = NULL,
    cells_highlight = NULL,
    cols_highlight = "#FF1A1A",
    cols_background = "#D9D9D9",
    legend_position = "right",
    legend_ncol = 1L,
    legend_dot_size = 4,
    legend_text_size = 8,
    legend_title_size = 9,
    legend_key_size = 0.32,
    legend_title = NULL,
    collect_guides = TRUE,
    legend_width = 0.24,
    panel_expand = 0.10,
    panel_margin = 10,
    panel_title_size = 11,
    base_size = 11,
    fixed_aspect = TRUE,
    ncol = NULL,
    title = NULL,
    title_size = 14,
    rds_file = NULL,
    analysis_date = Sys.Date(),
    show_provenance = TRUE,
    provenance_size = 6,
    provenance_color = "#666666",
    output_dir = NULL,
    filename = "review_umap.pdf",
    width = 12,
    height = 9
) {

  .review_umap_require_namespace("Seurat")
  .review_umap_require_namespace("ggplot2")
  .review_umap_require_namespace("patchwork")

  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  .review_umap_validate_flag(label, "label")
  .review_umap_validate_flag(repel, "repel")
  .review_umap_validate_flag(label_box, "label_box")
  .review_umap_validate_flag(shuffle, "shuffle")
  .review_umap_validate_flag(collect_guides, "collect_guides")
  .review_umap_validate_flag(fixed_aspect, "fixed_aspect")
  .review_umap_validate_flag(show_provenance, "show_provenance")

  if (!is.null(raster)) {
    .review_umap_validate_flag(raster, "raster")
  }

  .review_umap_validate_scalar_numeric(alpha, "alpha", 0, 1)
  .review_umap_validate_scalar_numeric(label_size, "label_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(legend_dot_size, "legend_dot_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(legend_text_size, "legend_text_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(legend_title_size, "legend_title_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(legend_key_size, "legend_key_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(legend_width, "legend_width", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(panel_expand, "panel_expand", 0, Inf)
  .review_umap_validate_scalar_numeric(panel_margin, "panel_margin", 0, Inf)
  .review_umap_validate_scalar_numeric(panel_title_size, "panel_title_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(base_size, "base_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(title_size, "title_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(provenance_size, "provenance_size", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(width, "width", 0, Inf, FALSE)
  .review_umap_validate_scalar_numeric(height, "height", 0, Inf, FALSE)

  if (!is.null(pt_size)) {
    .review_umap_validate_scalar_numeric(pt_size, "pt_size", 0, Inf, FALSE)
  }

  allowed_legend_positions <- c("right", "left", "top", "bottom", "none")
  if (
    length(legend_position) != 1L ||
      !is.character(legend_position) ||
      is.na(legend_position) ||
      !legend_position %in% allowed_legend_positions
  ) {
    stop(
      "`legend_position` must be one of: ",
      paste(allowed_legend_positions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  legend_ncol <- .review_umap_validate_positive_integer(legend_ncol, "legend_ncol")
  if (!is.null(ncol)) {
    ncol <- .review_umap_validate_positive_integer(ncol, "ncol")
  }

  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed)) {
    stop("`seed` must be one numeric value.", call. = FALSE)
  }
  seed <- as.integer(seed)

  .review_umap_validate_optional_string(group_by, "group_by")
  .review_umap_validate_optional_string(split_by, "split_by")
  .review_umap_validate_optional_string(assay, "assay")
  .review_umap_validate_optional_string(reduction, "reduction")
  .review_umap_validate_optional_string(legend_title, "legend_title", allow_empty = TRUE)
  .review_umap_validate_optional_string(title, "title", allow_empty = TRUE)
  .review_umap_validate_optional_string(rds_file, "rds_file", allow_empty = TRUE)
  .review_umap_validate_required_string(provenance_color, "provenance_color")
  .review_umap_validate_required_string(filename, "filename")

  analysis_date <- .review_umap_normalize_date(analysis_date)

  if (!is.null(assay)) {
    if (!assay %in% names(object@assays)) {
      stop("Assay `", assay, "` was not found.", call. = FALSE)
    }
    Seurat::DefaultAssay(object) <- assay
  }

  all_cells <- colnames(object)
  if (is.null(cells)) {
    selected_cells <- all_cells
  } else {
    cells <- unique(as.character(cells))
    missing_cells <- setdiff(cells, all_cells)
    if (length(missing_cells) > 0L) {
      warning(
        length(missing_cells),
        " cell name(s) were not found and were ignored.",
        call. = FALSE
      )
    }
    selected_cells <- intersect(cells, all_cells)
    if (length(selected_cells) == 0L) {
      stop("None of the requested cells were found.", call. = FALSE)
    }
  }

  metadata <- object[[]]
  internal_group_column <- ".review_umap_group"

  if (is.null(group_by)) {
    group_values <- as.character(Seurat::Idents(object))
    group_source_name <- "Identity"
  } else {
    if (!group_by %in% colnames(metadata)) {
      stop("Metadata column `", group_by, "` was not found.", call. = FALSE)
    }
    group_values <- as.character(metadata[[group_by]])
    group_source_name <- group_by
  }

  names(group_values) <- rownames(metadata)
  selected_group_values <- group_values[selected_cells]

  if (all(is.na(selected_group_values))) {
    stop("All grouping values are missing for the selected cells.", call. = FALSE)
  }

  group_levels_observed <- unique(selected_group_values[!is.na(selected_group_values)])
  group_levels <- .review_umap_resolve_order(
    observed = group_levels_observed,
    requested = group_order,
    argument_name = "group_order"
  )

  object[[internal_group_column]] <- factor(group_values, levels = group_levels)

  if (is.null(split_by)) {
    panel_cells <- list(All = selected_cells)
    panel_names <- "All"
  } else {
    if (!split_by %in% colnames(metadata)) {
      stop("Metadata column `", split_by, "` was not found.", call. = FALSE)
    }

    split_values <- as.character(metadata[[split_by]])
    names(split_values) <- rownames(metadata)
    selected_split_values <- split_values[selected_cells]

    if (all(is.na(selected_split_values))) {
      stop("All split values are missing for the selected cells.", call. = FALSE)
    }

    split_levels_observed <- unique(selected_split_values[!is.na(selected_split_values)])
    split_levels <- .review_umap_resolve_order(
      observed = split_levels_observed,
      requested = split_order,
      argument_name = "split_order"
    )

    panel_cells <- lapply(split_levels, function(level_i) {
      selected_cells[
        !is.na(selected_split_values) & selected_split_values == level_i
      ]
    })
    names(panel_cells) <- split_levels

    nonempty <- lengths(panel_cells) > 0L
    if (any(!nonempty)) {
      warning(
        "Empty split levels were omitted: ",
        paste(names(panel_cells)[!nonempty], collapse = ", "),
        call. = FALSE
      )
      panel_cells <- panel_cells[nonempty]
    }

    if (length(panel_cells) == 0L) {
      stop("No non-empty split panels were available.", call. = FALSE)
    }

    panel_names <- names(panel_cells)
  }

  group_colors <- .review_umap_resolve_colors(colors, group_levels)
  legend_title_final <- if (is.null(legend_title)) group_source_name else legend_title

  common_color_scale <- ggplot2::scale_colour_manual(
    values = group_colors,
    limits = group_levels,
    breaks = group_levels,
    drop = FALSE,
    name = legend_title_final,
    na.value = "#BDBDBD"
  )

  common_fill_scale <- ggplot2::scale_fill_manual(
    values = group_colors,
    limits = group_levels,
    breaks = group_levels,
    drop = FALSE,
    name = legend_title_final,
    na.value = "#BDBDBD"
  )

  plot_list <- vector("list", length(panel_cells))
  names(plot_list) <- panel_names

  for (i in seq_along(panel_cells)) {
    panel_cells_i <- panel_cells[[i]]
    panel_name_i <- panel_names[[i]]

    dimplot_args <- list(
      object = object,
      reduction = reduction,
      cells = panel_cells_i,
      group.by = internal_group_column,
      pt.size = pt_size,
      alpha = alpha,
      shuffle = shuffle,
      seed = seed,
      label = label,
      label.size = label_size,
      repel = repel,
      label.box = label_box,
      raster = raster,
      combine = TRUE
    )

    if (!is.null(cells_highlight)) {
      dimplot_args$cells.highlight <- .review_umap_subset_highlight(
        cells_highlight,
        panel_cells_i
      )
      dimplot_args$cols.highlight <- cols_highlight
      dimplot_args$cols <- cols_background
    }

    dimplot_args <- dimplot_args[
      !vapply(dimplot_args, is.null, logical(1))
    ]

    set.seed(seed)
    p <- do.call(Seurat::DimPlot, dimplot_args)

    if (is.null(cells_highlight)) {
      p <- p + common_color_scale + common_fill_scale
    }

    if (!is.null(split_by)) {
      p <- p + ggplot2::labs(title = panel_name_i)
    } else {
      p <- p + ggplot2::labs(title = NULL)
    }

    p <- p +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = panel_expand)
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = panel_expand)
      )

    if (isTRUE(fixed_aspect)) {
      p <- p + ggplot2::coord_fixed(ratio = 1, clip = "off")
    } else {
      p <- p + ggplot2::coord_cartesian(clip = "off")
    }

    panel_legend_position <- if (isTRUE(collect_guides)) "none" else legend_position

    p <- p +
      .review_umap_theme(
        base_size = base_size,
        panel_margin = panel_margin,
        panel_title_size = panel_title_size,
        legend_position = panel_legend_position,
        legend_text_size = legend_text_size,
        legend_title_size = legend_title_size,
        legend_key_size = legend_key_size
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          ncol = legend_ncol,
          byrow = TRUE,
          override.aes = list(size = legend_dot_size, alpha = 1)
        ),
        fill = ggplot2::guide_legend(
          ncol = legend_ncol,
          byrow = TRUE,
          override.aes = list(size = legend_dot_size, alpha = 1)
        )
      )

    plot_list[[i]] <- p
  }

  number_of_panels <- length(plot_list)
  ncol_final <- if (is.null(ncol)) {
    .review_umap_default_ncol(number_of_panels)
  } else {
    ncol
  }

  panel_grid <- patchwork::wrap_plots(
    plotlist = plot_list,
    ncol = ncol_final,
    guides = "keep"
  )

  if (isTRUE(collect_guides) && legend_position != "none" && is.null(cells_highlight)) {
    legend_plot <- .review_umap_make_legend_plot(
      group_levels = group_levels,
      group_colors = group_colors,
      legend_title = legend_title_final,
      legend_position = legend_position,
      legend_ncol = legend_ncol,
      legend_dot_size = legend_dot_size,
      legend_text_size = legend_text_size,
      legend_title_size = legend_title_size,
      legend_key_size = legend_key_size
    )

    if (legend_position %in% c("right", "left")) {
      if (legend_position == "right") {
        assembled_plot <- (panel_grid | legend_plot) +
          patchwork::plot_layout(widths = c(1, legend_width))
      } else {
        assembled_plot <- (legend_plot | panel_grid) +
          patchwork::plot_layout(widths = c(legend_width, 1))
      }
    } else {
      if (legend_position == "top") {
        assembled_plot <- (legend_plot / panel_grid) +
          patchwork::plot_layout(heights = c(legend_width, 1))
      } else {
        assembled_plot <- (panel_grid / legend_plot) +
          patchwork::plot_layout(heights = c(1, legend_width))
      }
    }
  } else {
    assembled_plot <- panel_grid
  }

  provenance_caption <- if (isTRUE(show_provenance)) {
    .review_umap_provenance_text(
      rds_file = rds_file,
      analysis_date = analysis_date
    )
  } else {
    NULL
  }

  assembled_plot <- assembled_plot +
    patchwork::plot_annotation(
      title = if (!is.null(title) && nzchar(title)) title else NULL,
      caption = provenance_caption,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = title_size,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 8, unit = "pt")
        ),
        plot.caption = ggplot2::element_text(
          size = provenance_size,
          colour = provenance_color,
          hjust = 1,
          margin = ggplot2::margin(t = 5, r = 2, unit = "pt")
        )
      )
    )

  output_files <- character(0)

  if (!is.null(output_dir)) {
    .review_umap_validate_required_string(output_dir, "output_dir")
    output_dir <- path.expand(output_dir)

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(output_dir)) {
      stop("Failed to create output directory: ", output_dir, call. = FALSE)
    }

    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }

    output_path <- normalizePath(
      file.path(output_dir, filename),
      mustWork = FALSE
    )

    grDevices::pdf(
      file = output_path,
      width = width,
      height = height,
      onefile = TRUE,
      bg = "white",
      useDingbats = FALSE
    )
    print(assembled_plot)
    grDevices::dev.off()

    output_files <- output_path
  }

  result <- list(
    plot = assembled_plot,
    plots = plot_list,
    files = output_files,
    group_by = group_source_name,
    split_by = split_by,
    group_levels = group_levels,
    split_levels = if (is.null(split_by)) NULL else panel_names,
    colors = group_colors,
    parameters = list(
      assay = assay,
      reduction = reduction,
      pt_size = pt_size,
      alpha = alpha,
      shuffle = shuffle,
      seed = seed,
      label = label,
      label_size = label_size,
      repel = repel,
      label_box = label_box,
      raster = raster,
      legend_position = legend_position,
      legend_ncol = legend_ncol,
      legend_dot_size = legend_dot_size,
      legend_text_size = legend_text_size,
      legend_title_size = legend_title_size,
      legend_key_size = legend_key_size,
      legend_title = legend_title_final,
      collect_guides = collect_guides,
      legend_width = legend_width,
      panel_expand = panel_expand,
      panel_margin = panel_margin,
      panel_title_size = panel_title_size,
      base_size = base_size,
      fixed_aspect = fixed_aspect,
      ncol = ncol_final,
      title = title,
      title_size = title_size,
      rds_file = rds_file,
      analysis_date = analysis_date,
      show_provenance = show_provenance,
      provenance_size = provenance_size,
      provenance_color = provenance_color,
      width = width,
      height = height
    )
  )

  class(result) <- c("review_umap_result", "list")
  result
}


#' Print a review UMAP result
#'
#' @param x A `review_umap_result` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.review_umap_result <- function(x, ...) {
  cat("\n")
  cat("uenoyscRNA UMAP Review\n")
  cat("----------------------\n")
  cat(sprintf("%-20s: %d\n", "Panels generated", length(x$plots)))
  cat(sprintf("%-20s: %s\n", "Grouping variable", x$group_by))
  cat(sprintf(
    "%-20s: %s\n",
    "Split variable",
    if (is.null(x$split_by)) "None" else x$split_by
  ))
  cat(sprintf("%-20s: %s\n", "Legend position", x$parameters$legend_position))
  cat(sprintf("%-20s: %d\n", "PDF files written", length(x$files)))
  cat("\n")
  invisible(x)
}


.review_umap_make_legend_plot <- function(
    group_levels,
    group_colors,
    legend_title,
    legend_position,
    legend_ncol,
    legend_dot_size,
    legend_text_size,
    legend_title_size,
    legend_key_size
) {
  legend_data <- data.frame(
    x = seq_along(group_levels),
    y = 1,
    group = factor(group_levels, levels = group_levels)
  )

  ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = x, y = y, colour = group)
  ) +
    ggplot2::geom_point(size = 0.01, alpha = 0) +
    ggplot2::scale_colour_manual(
      values = group_colors,
      limits = group_levels,
      breaks = group_levels,
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        byrow = TRUE,
        override.aes = list(size = legend_dot_size, alpha = 1)
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.height = grid::unit(legend_key_size, "cm"),
      legend.key.width = grid::unit(legend_key_size, "cm"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
    )
}


.review_umap_normalize_date <- function(analysis_date) {
  if (inherits(analysis_date, "Date")) {
    if (length(analysis_date) != 1L || is.na(analysis_date)) {
      stop("`analysis_date` must contain one non-missing date.", call. = FALSE)
    }
    return(analysis_date)
  }

  if (length(analysis_date) != 1L || is.na(analysis_date)) {
    stop(
      "`analysis_date` must be one Date value or a date string in YYYY-MM-DD format.",
      call. = FALSE
    )
  }

  parsed_date <- suppressWarnings(as.Date(as.character(analysis_date)))
  if (is.na(parsed_date)) {
    stop(
      "`analysis_date` could not be interpreted as a date. Use YYYY-MM-DD.",
      call. = FALSE
    )
  }

  parsed_date
}


.review_umap_provenance_text <- function(
    rds_file = NULL,
    analysis_date = Sys.Date()
) {
  rds_text <- if (
    is.null(rds_file) ||
      length(rds_file) == 0L ||
      is.na(rds_file) ||
      !nzchar(rds_file)
  ) {
    "Not specified"
  } else {
    basename(path.expand(rds_file))
  }

  paste0(
    "RDS: ",
    rds_text,
    "    |    Analysis date: ",
    format(.review_umap_normalize_date(analysis_date), "%Y-%m-%d")
  )
}


.review_umap_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package `", package, "` is required.", call. = FALSE)
  }
  invisible(TRUE)
}


.review_umap_validate_flag <- function(value, argument_name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop("`", argument_name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}


.review_umap_validate_scalar_numeric <- function(
    value,
    argument_name,
    lower = -Inf,
    upper = Inf,
    lower_inclusive = TRUE,
    upper_inclusive = TRUE
) {
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value)
  ) {
    stop("`", argument_name, "` must be one finite numeric value.", call. = FALSE)
  }

  lower_invalid <- if (isTRUE(lower_inclusive)) value < lower else value <= lower
  upper_invalid <- if (isTRUE(upper_inclusive)) value > upper else value >= upper

  if (lower_invalid || upper_invalid) {
    stop("`", argument_name, "` is outside the allowed range.", call. = FALSE)
  }

  invisible(TRUE)
}


.review_umap_validate_positive_integer <- function(value, argument_name) {
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      value < 1 ||
      value != as.integer(value)
  ) {
    stop("`", argument_name, "` must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}


.review_umap_validate_optional_string <- function(
    value,
    argument_name,
    allow_empty = FALSE
) {
  if (is.null(value)) {
    return(invisible(TRUE))
  }

  invalid <- length(value) != 1L || !is.character(value) || is.na(value)
  if (!allow_empty) {
    invalid <- invalid || !nzchar(value)
  }

  if (invalid) {
    stop("`", argument_name, "` must be NULL or one character value.", call. = FALSE)
  }
  invisible(TRUE)
}


.review_umap_validate_required_string <- function(value, argument_name) {
  if (
    length(value) != 1L ||
      !is.character(value) ||
      is.na(value) ||
      !nzchar(value)
  ) {
    stop("`", argument_name, "` must be one non-empty character value.", call. = FALSE)
  }
  invisible(TRUE)
}


.review_umap_resolve_order <- function(observed, requested = NULL, argument_name) {
  observed <- unique(as.character(observed))
  observed <- observed[!is.na(observed)]

  if (is.null(requested)) {
    return(observed)
  }

  requested <- unique(as.character(requested))
  requested <- requested[!is.na(requested)]

  unknown <- setdiff(requested, observed)
  if (length(unknown) > 0L) {
    stop(
      "`", argument_name, "` contains values not present in selected cells: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  c(requested, setdiff(observed, requested))
}


.review_umap_resolve_colors <- function(colors, group_levels) {
  number_of_groups <- length(group_levels)
  if (number_of_groups < 1L) {
    stop("No grouping levels were available.", call. = FALSE)
  }

  if (is.null(colors)) {
    generated <- grDevices::hcl.colors(number_of_groups, palette = "Dynamic")
    names(generated) <- group_levels
    return(generated)
  }

  colors <- as.character(colors)

  if (!is.null(names(colors))) {
    missing_names <- setdiff(group_levels, names(colors))
    if (length(missing_names) > 0L) {
      stop(
        "The named `colors` vector is missing groups: ",
        paste(missing_names, collapse = ", "),
        call. = FALSE
      )
    }
    resolved <- colors[group_levels]
  } else {
    if (length(colors) < number_of_groups) {
      stop(
        "`colors` contains ", length(colors), " colors, but ",
        number_of_groups, " groups must be colored.",
        call. = FALSE
      )
    }
    resolved <- colors[seq_len(number_of_groups)]
    names(resolved) <- group_levels
  }

  resolved
}


.review_umap_subset_highlight <- function(cells_highlight, panel_cells) {
  if (is.list(cells_highlight)) {
    result <- lapply(cells_highlight, function(x) {
      intersect(as.character(x), panel_cells)
    })
    names(result) <- names(cells_highlight)
    return(result)
  }

  intersect(as.character(cells_highlight), panel_cells)
}


.review_umap_default_ncol <- function(number_of_panels) {
  if (number_of_panels <= 1L) return(1L)
  if (number_of_panels == 2L) return(2L)
  if (number_of_panels <= 4L) return(2L)
  if (number_of_panels <= 9L) return(3L)
  ceiling(sqrt(number_of_panels))
}


.review_umap_theme <- function(
    base_size,
    panel_margin,
    panel_title_size,
    legend_position,
    legend_text_size,
    legend_title_size,
    legend_key_size
) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = panel_title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 5, unit = "pt")
      ),
      plot.margin = ggplot2::margin(
        panel_margin,
        panel_margin,
        panel_margin,
        panel_margin,
        unit = "pt"
      ),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.height = grid::unit(legend_key_size, "cm"),
      legend.key.width = grid::unit(legend_key_size, "cm"),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.box.spacing = grid::unit(4, "pt"),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
