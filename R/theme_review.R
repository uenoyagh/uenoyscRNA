# Internal helper: standard review plot theme
.review_theme <- function(
    base_size = 11,
    remove_axes = TRUE,
    centered_title = TRUE
) {
  plot_theme <- ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = if (isTRUE(centered_title)) 0.5 else 0
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      )
    )

  if (isTRUE(remove_axes)) {
    plot_theme <- plot_theme +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  plot_theme
}
