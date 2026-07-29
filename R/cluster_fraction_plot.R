# ============================================================
# Cluster-fraction plotting functions
# uenoy scRNAseq Framework v2.6
# ============================================================

cf_theme <- function() {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none"
    )
}

cf_add_cluster_labels <- function(p, df, y_column, label_last_only = TRUE, label_size = 2.6) {
  label_df <- df
  if (isTRUE(label_last_only)) {
    last_level <- tail(levels(df$condition), 1)
    label_df <- df[as.character(df$condition) == last_level & df[[y_column]] > 0, , drop = FALSE]
  } else {
    label_df <- df[df[[y_column]] > 0, , drop = FALSE]
  }
  if (nrow(label_df) == 0) return(p)
  p + ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(label = cluster),
    size = label_size,
    hjust = -0.15,
    check_overlap = TRUE,
    show.legend = FALSE
  ) + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 22, 5.5, 5.5))
}

cf_make_cluster_line_plot <- function(df,
                                      y_column,
                                      y_label,
                                      title,
                                      subtitle,
                                      line_width = 0.55,
                                      point_size = 1.8,
                                      facet_ncol = 2,
                                      free_y = TRUE,
                                      label_clusters = TRUE,
                                      label_last_only = TRUE,
                                      label_size = 2.6) {
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = .data[[y_column]],
      group = cluster,
      colour = annotation
    )
  ) +
    ggplot2::geom_line(linewidth = line_width, na.rm = TRUE) +
    ggplot2::geom_point(size = point_size, na.rm = TRUE) +
    ggplot2::facet_wrap(
      ggplot2::vars(annotation),
      ncol = facet_ncol,
      scales = if (isTRUE(free_y)) "free_y" else "fixed"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = y_label) +
    cf_theme()

  if (isTRUE(label_clusters)) {
    p <- cf_add_cluster_labels(p, df, y_column, label_last_only, label_size)
  }
  p
}

cf_make_stacked_cluster_plot <- function(df) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = condition, y = fraction_total_percent, fill = cluster)
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::facet_wrap(ggplot2::vars(annotation), ncol = 2) +
    ggplot2::labs(
      title = "Cluster composition by macrophage annotation",
      subtitle = "Stacked cluster fractions among total analyzed cells",
      x = NULL,
      y = "Fraction among total cells (%)",
      fill = "Cluster"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

cf_make_stacked_annotation_plot <- function(df) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = condition, y = fraction_total_percent, fill = annotation)
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::labs(
      title = "Macrophage annotation composition",
      subtitle = "Fraction of each annotation among total analyzed cells",
      x = NULL,
      y = "Fraction among total cells (%)",
      fill = "Annotation"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

cf_save_plot_dual <- function(plot, base_path, width, height,
                              export_pdf = TRUE, export_png = TRUE,
                              png_dpi = 300, png_background = "white",
                              overwrite = FALSE) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(export_pdf)) {
    path <- paste0(base_path, ".pdf")
    if (!file.exists(path) || isTRUE(overwrite)) {
      grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
      print(plot)
      grDevices::dev.off()
    }
  }
  if (isTRUE(export_png)) {
    path <- paste0(base_path, ".png")
    if (!file.exists(path) || isTRUE(overwrite)) {
      ggplot2::ggsave(path, plot = plot, width = width, height = height,
                      units = "in", dpi = png_dpi, bg = png_background,
                      limitsize = FALSE)
    }
  }
  invisible(base_path)
}
