# Install Framework v3.0 generic cell-fraction module.
source_dir <- "/path/to/uenoyscRNA_Framework_v3.0_cell_fraction"
target_root <- "/Users/uenoya/Projects/uenoyscRNA"

files <- c(
  "analysis/06_cell_fraction_transition.R",
  "R/cell_fraction_engine.R",
  "R/cell_fraction_plot.R",
  "config/cell_fraction_config.R"
)

for (rel in files) {
  from <- file.path(source_dir, rel)
  to <- file.path(target_root, rel)

  if (!file.exists(from)) {
    stop("Source file not found: ", from)
  }

  dir.create(
    dirname(to),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!file.copy(from, to, overwrite = TRUE)) {
    stop("Failed to copy: ", rel)
  }

  message("Installed: ", to)
}

message("Framework v3.0 cell-fraction module installation completed.")
message("project_config.R was not modified.")
message(
  "Run: source('/Users/uenoya/Projects/uenoyscRNA/",
  "analysis/06_cell_fraction_transition.R')"
)
