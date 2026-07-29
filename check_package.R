# Run from the package root directory.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools first: install.packages('devtools')")
}

devtools::document()

devtools::check(
  ".",
  cran = TRUE,
  error_on = "warning"
)
