# Run from /Users/uenoya/Projects/uenoyscRNA
# This script does not initialize a new renv project automatically.
# It restores the existing lockfile when present and checks dependencies.

required <- c(
  "Seurat", "SeuratObject", "Matrix", "presto", "future",
  "future.apply", "ggplot2", "patchwork", "dplyr", "tibble",
  "scales", "openxlsx", "here", "renv"
)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

if (file.exists("renv.lock")) {
  message("Existing renv.lock detected.")
  message("Use renv::restore() when package restoration is required.")
} else {
  message("No renv.lock detected. Existing Framework v1.0 policy should determine snapshot timing.")
}

missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  message("Missing packages: ", paste(missing, collapse = ", "))
  message("Install within the active renv environment, then run renv::snapshot().")
} else {
  message("All v3.2 dependencies are available.")
}
