# ============================================================
# Generic cell-fraction plotting functions
# uenoy scRNAseq Framework v3.0
# ============================================================

cf_theme <- function(show_legend = FALSE) {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16, face = "bold", hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10, hjust = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "white", colour = "black", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(
        size = 9, face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = if (isTRUE(show_legend)) "right" else "none"
    )
}

cf_add_feature_labels <- function(
    p,
    df,
    y_column,
    label_last_only = TRUE,
    label_size = 2.6) {

  label_df <- df
  if (isTRUE(label_last_only)) {
    condition_levels <- levels(df$condition)
    last_level <- tail(condition_levels, 1L)
    label_df <- df[
      as.character(df$condition) == last_level &
        is.finite(df[[y_column]]) &
        df[[y_column]] > 0,
      , drop = FALSE
    ]
  } else {
    label_df <- df[
      is.finite(df[[y_column]]) &
        df[[y_column]] > 0,
      , drop = FALSE
    ]
  }

  if (nrow(label_df) == 0L) return(p)

  p +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = feature),
      size = label_size,
      hjust = -0.15,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 5.5)
    )
}

cf_make_feature_line_plot <- function(
    df,
    y_column,
    y_label,
    title,
    subtitle,
    line_width = 0.55,
    point_size = 1.8,
    facet_ncol = 2,
    free_y = TRUE,
    label_features = TRUE,
    label_last_only = TRUE,
    label_size = 2.6) {

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = .data[[y_column]],
      group = feature,
      colour = parent
    )
  ) +
    ggplot2::geom_line(
      linewidth = line_width,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = point_size,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(parent),
      ncol = facet_ncol,
      scales = if (isTRUE(free_y)) "free_y" else "fixed"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_label
    ) +
    cf_theme(show_legend = FALSE)

  if (isTRUE(label_features)) {
    p <- cf_add_feature_labels(
      p,
      df,
      y_column,
      label_last_only,
      label_size
    )
  }
  p
}

cf_make_stacked_feature_plot <- function(
    df,
    denominator_label = "total cells",
    feature_label = "cluster",
    parent_label = "annotation") {

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = fraction_total_percent,
      fill = feature
    )
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::facet_wrap(
      ggplot2::vars(parent),
      ncol = 2
    ) +
    ggplot2::labs(
      title = paste0(feature_label, " composition by ", parent_label),
      subtitle = paste0(
        "Stacked fractions among ", denominator_label
      ),
      x = NULL,
      y = paste0("Fraction among ", denominator_label, " (%)"),
      fill = feature_label
    ) +
    cf_theme(show_legend = TRUE)
}

cf_make_stacked_parent_plot <- function(
    df,
    denominator_label = "total cells",
    parent_label = "cell type") {

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = fraction_total_percent,
      fill = parent
    )
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::labs(
      title = paste0(parent_label, " composition"),
      subtitle = paste0(
        "Fraction of each ", parent_label,
        " among ", denominator_label
      ),
      x = NULL,
      y = paste0("Fraction among ", denominator_label, " (%)"),
      fill = parent_label
    ) +
    cf_theme(show_legend = TRUE)
}

cf_make_fraction_heatmap <- function(
    df,
    value_column = "fraction_total_percent",
    low_color = "#0033FF",
    mid_color = "#FFFFFF",
    high_color = "#FF1A1A",
    midpoint = NULL,
    title = "Cell composition heatmap") {

  values <- df[[value_column]]
  if (is.null(midpoint)) {
    midpoint <- stats::median(values, na.rm = TRUE)
  }
  if (!is.finite(midpoint)) midpoint <- 0

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = parent,
      fill = .data[[value_column]]
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = midpoint,
      name = "Fraction (%)"
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        face = "bold"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

cf_save_plot_dual <- function(
    plot,
    base_path,
    width,
    height,
    export_pdf = TRUE,
    export_png = TRUE,
    png_dpi = 300,
    png_background = "white",
    overwrite = FALSE) {

  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(export_pdf)) {
    path <- paste0(base_path, ".pdf")
    if (!file.exists(path) || isTRUE(overwrite)) {
      grDevices::pdf(
        path,
        width = width,
        height = height,
        useDingbats = FALSE
      )
      tryCatch(
        print(plot),
        finally = grDevices::dev.off()
      )
    }
  }

  if (isTRUE(export_png)) {
    path <- paste0(base_path, ".png")
    if (!file.exists(path) || isTRUE(overwrite)) {
      ggplot2::ggsave(
        filename = path,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = png_dpi,
        bg = png_background,
        limitsize = FALSE
      )
    }
  }

  invisible(base_path)
}
