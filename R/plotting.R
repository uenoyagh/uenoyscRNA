# ============================================================
# Plotting helpers
# uenoy scRNAseq Framework v2.3
# ============================================================

make_umap_output_subdirs <- function(root) {
  dirs <- file.path(
    root,
    c(
      "R8/00_overview",
      "R8/01_sample",
      "R8/02_condition",
      "R8/03_annotation",
      "R8/04_cluster",
      "R8/05_split_by_sample",
      "R8/06_split_by_condition",
      "Pastel/00_overview",
      "Pastel/01_sample",
      "Pastel/02_condition",
      "Pastel/03_annotation",
      "Pastel/04_cluster",
      "Pastel/05_split_by_sample",
      "Pastel/06_split_by_condition",
      "Monochrome_transparent/07_each_sample",
      "08_tables",
      "09_logs"
    )
  )

  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  names(dirs) <- c(
    "r8_overview", "r8_sample", "r8_condition", "r8_annotation",
    "r8_cluster", "r8_split_sample", "r8_split_condition",
    "pastel_overview", "pastel_sample", "pastel_condition",
    "pastel_annotation", "pastel_cluster",
    "pastel_split_sample", "pastel_split_condition",
    "each_sample", "tables", "logs"
  )

  dirs
}
