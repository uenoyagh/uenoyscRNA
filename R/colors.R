#' Blue-white-red color palette
#'
#' Returns a blue-white-red color palette for continuous expression values.
#'
#' @param n Number of colors to generate.
#' @param low Color used for low values.
#' @param mid Color used for midpoint values.
#' @param high Color used for high values.
#'
#' @return A character vector containing hexadecimal colors.
#' @export
#'
#' @examples
#' ueno_blue_white_red()
#' ueno_blue_white_red(11)
ueno_blue_white_red <- function(
    n = 101L,
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A"
) {
  if (
    length(n) != 1L ||
    is.na(n) ||
    n < 2L ||
    n != as.integer(n)
  ) {
    stop(
      "`n` must be a single integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  grDevices::colorRampPalette(
    colors = c(low, mid, high)
  )(as.integer(n))
}
