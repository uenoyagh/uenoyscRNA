#' Publication theme for scRNA-seq figures
#'
#' Provides a clean publication-style ggplot2 theme for UMAP,
#' violin plots, dot plots, feature plots, and related figures.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @param legend_position Legend position passed to ggplot2.
#' @param panel_border Logical. Whether to draw a panel border.
#' @param axis_text Logical. Whether to show axis text.
#' @param axis_ticks Logical. Whether to show axis ticks.
#'
#' @return A ggplot2 theme object.
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_ueno_scRNA()
theme_ueno_scRNA <- function(
    base_size = 11,
    base_family = "",
    legend_position = "right",
    panel_border = TRUE,
    axis_text = TRUE,
    axis_ticks = TRUE
) {

  validate_positive_number(
    x = base_size,
    argument = "base_size"
  )

  if (
    !is.character(base_family) ||
    length(base_family) != 1L ||
    is.na(base_family)
  ) {
    stop(
      "`base_family` must be a single character string.",
      call. = FALSE
    )
  }

  theme <- ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.15),
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        colour = "black"
      ),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.margin = ggplot2::margin(
        t = 6,
        r = 6,
        b = 6,
        l = 6
      )
    )

  if (isTRUE(panel_border)) {

    theme <- theme +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          fill = NA,
          colour = "black",
          linewidth = 0.6
        )
      )

  } else {

    theme <- theme +
      ggplot2::theme(
        panel.border = ggplot2::element_blank()
      )
  }

  if (!isTRUE(axis_text)) {

    theme <- theme +
      ggplot2::theme(
        axis.text = ggplot2::element_blank()
      )
  }

  if (!isTRUE(axis_ticks)) {

    theme <- theme +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank()
      )
  }

  theme
}
