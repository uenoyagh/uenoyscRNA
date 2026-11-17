#!/usr/bin/env Rscript

# Install only the packages needed for the DecontX sensitivity analysis.
# Does NOT run renv::snapshot() and does NOT update unrelated packages.

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

need <- c("decontX", "SingleCellExperiment")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  BiocManager::install(
    missing,
    ask = FALSE,
    update = FALSE
  )
} else {
  message("Required packages are already installed.")
}

message("=== package versions ===")
for (p in need) {
  message(p, ": ", as.character(packageVersion(p)))
}
