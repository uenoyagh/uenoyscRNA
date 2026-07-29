# Install v2.6 cluster-fraction module into the existing framework.
source_dir <- "/path/to/uenoyscRNA_Framework_v2.6_cluster_fraction"
target_root <- "/Users/uenoya/Projects/uenoyscRNA"

files <- c(
  "analysis/06_cluster_fraction_transition.R",
  "R/cluster_fraction_engine.R",
  "R/cluster_fraction_plot.R",
  "config/cluster_fraction_config.R"
)

for (rel in files) {
  from <- file.path(source_dir, rel)
  to <- file.path(target_root, rel)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Failed to copy: ", rel)
  message("Installed: ", to)
}
