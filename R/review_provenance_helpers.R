# Internal provenance helpers for uenoyscRNA plotting functions.
#
# Add this file under R/ when several plotting functions should share the same
# RDS filename / analysis-date footer implementation.

.review_normalize_analysis_date <- function(analysis_date = Sys.Date()) {
  if (inherits(analysis_date, "Date")) {
    if (length(analysis_date) != 1L || is.na(analysis_date)) {
      stop("`analysis_date` must contain one non-missing date.", call. = FALSE)
    }
    return(analysis_date)
  }

  if (length(analysis_date) != 1L || is.na(analysis_date)) {
    stop(
      "`analysis_date` must be one Date value or a YYYY-MM-DD date string.",
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


.review_provenance_text <- function(
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
    format(.review_normalize_analysis_date(analysis_date), "%Y-%m-%d")
  )
}


.review_add_ggplot_provenance <- function(
    plot,
    rds_file = NULL,
    analysis_date = Sys.Date(),
    show_provenance = TRUE,
    provenance_size = 6,
    provenance_color = "#666666"
) {
  if (!isTRUE(show_provenance)) {
    return(plot)
  }

  plot +
    ggplot2::labs(
      caption = .review_provenance_text(
        rds_file = rds_file,
        analysis_date = analysis_date
      )
    ) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        size = provenance_size,
        colour = provenance_color,
        hjust = 1,
        margin = ggplot2::margin(t = 5, r = 2, unit = "pt")
      )
    )
}
