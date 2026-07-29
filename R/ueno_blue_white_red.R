#' Blue-white-red color palette
#'
#' Generate a blue-white-red color palette with the requested number
#' of colors.
#'
#' @param n A single integer greater than or equal to 2.
#'
#' @return A character vector of hexadecimal color codes.
#'
#' @export
#'
#' @examples
#' ueno_blue_white_red(11)
ueno_blue_white_red <- function(n) {

  if (
    length(n) != 1L ||
    is.na(n) ||
    !is.numeric(n) ||
    !is.finite(n) ||
    n != as.integer(n) ||
    n < 2L
  ) {
    stop(
      "`n` must be a single integer greater than or equal to 2.",
      call. = FALSE
    )
  }

  n <- as.integer(n)

  grDevices::colorRampPalette(
    c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  )(n)
}
