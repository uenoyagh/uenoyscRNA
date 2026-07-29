# Internal helper: save a plot as PDF
.save_review_pdf <- function(
    plot,
    filename,
    width,
    height,
    bg = "white"
) {
  if (!inherits(plot, c("ggplot", "patchwork"))) {
    stop(
      "`plot` must inherit from class `ggplot` or `patchwork`.",
      call. = FALSE
    )
  }

  if (!is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !nzchar(filename)) {
    stop(
      "`filename` must be one non-empty character value.",
      call. = FALSE
    )
  }

  output_dir <- dirname(filename)

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )

  if (!file.exists(filename)) {
    stop(
      paste0("PDF file was not created: ", filename),
      call. = FALSE
    )
  }

  normalizePath(
    filename,
    winslash = "/",
    mustWork = TRUE
  )
}
