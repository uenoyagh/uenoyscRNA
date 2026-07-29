# ============================================================
# Project-wide RDS discovery and registered-dataset inventory
# uenoy scRNAseq Framework v4.0
# ============================================================

rm(list = ls())
gc()

script_root <- Sys.getenv("UENOY_SCRNA_ROOT", unset = "/Users/uenoya/Projects/uenoyscRNA")
if (!dir.exists(script_root) && file.exists(file.path(getwd(), "DESCRIPTION"))) {
  script_root <- normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(script_root, "R", "rds_registry.R"))
source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "seurat_helpers.R"))
source(file.path(script_root, "R", "metadata.R"))
source(file.path(script_root, "R", "utils.R"))

# ------------------------------------------------------------
# Inventory settings: safe defaults
# ------------------------------------------------------------
if (!exists("inventory_max_unique_values", inherits = FALSE) ||
    length(inventory_max_unique_values) != 1L ||
    !is.numeric(inventory_max_unique_values) ||
    is.na(inventory_max_unique_values) ||
    inventory_max_unique_values < 1) {
  inventory_max_unique_values <- 50L
}

inventory_max_unique_values <- as.integer(inventory_max_unique_values)

if (!exists("inventory_save_rds_summary", inherits = FALSE) ||
    length(inventory_save_rds_summary) != 1L ||
    is.na(inventory_save_rds_summary)) {
  inventory_save_rds_summary <- TRUE
}

inventory_save_rds_summary <- isTRUE(inventory_save_rds_summary)

catalog <- discover_rds_files(project_config$external_data_root)
registry_summary <- summarize_dataset_registry(
  dataset_registry, project_config$external_data_root, catalog
)
print_dataset_registry_summary(registry_summary)

output_dir <- file.path(project_config$external_data_root, "analysis_results", "00_project_inventory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(catalog, file.path(output_dir, "RDS_project_catalog.csv"), row.names = FALSE, na = "")
utils::write.csv(registry_summary, file.path(output_dir, "dataset_registry_summary.csv"), row.names = FALSE, na = "")

# Deep inspection is deliberately limited to the resolved registered files,
# avoiding accidental loading of every historical bundle and large RDS.
resolved <- unique(unlist(lapply(
  names(dataset_registry), get_dataset_files, catalog = catalog, strict = FALSE
)))

inventory_rows <- metadata_rows <- layer_rows <- candidate_rows <- list()
object_summaries <- list()

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("Required package is missing: SeuratObject", call. = FALSE)
}

for (i in seq_along(resolved)) {
  path <- resolved[[i]]
  fn <- basename(path)
  cat("[", i, "/", length(resolved), "] ", fn, "\n", sep = "")
  obj <- safe_read_rds(path)
  if (inherits(obj, "rds_read_error")) {
    inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
      file = fn, full_path = path, file_size_gb = file_size_gb(path),
      read_status = "ERROR", object_class = NA, seurat_object = FALSE,
      n_cells = NA, n_features = NA, default_assay = NA,
      assays = NA, reductions = NA, metadata_columns = NA,
      error_message = obj$message, stringsAsFactors = FALSE
    )
    next
  }
  seu <- is_seurat_object(obj)
  meta <- if (seu) obj[[]] else NULL
  assays <- if (seu) get_assay_names(obj) else character()
  reductions <- if (seu) get_reduction_names(obj) else character()
  inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
    file = fn, full_path = path, file_size_gb = file_size_gb(path),
    read_status = "OK", object_class = paste(class(obj), collapse = "/"),
    seurat_object = seu, n_cells = if (seu) ncol(obj) else NA,
    n_features = if (seu) nrow(obj) else NA,
    default_assay = if (seu) get_default_assay_safe(obj) else NA,
    assays = collapse_or_na(assays), reductions = collapse_or_na(reductions),
    metadata_columns = if (seu) ncol(meta) else NA,
    error_message = NA, stringsAsFactors = FALSE
  )
  if (seu) {
    z <- get_layer_names(obj)
    if (nrow(z)) { z$file <- fn; layer_rows[[length(layer_rows)+1L]] <- z[, c("file", "assay", "layer")] }
    m <- metadata_column_summary(meta, inventory_max_unique_values)
    m$file <- fn
    metadata_rows[[length(metadata_rows)+1L]] <- m[, c("file", "column", "class", "n_non_na", "n_na", "n_unique", "example_values")]
    cand <- detect_metadata_candidates(colnames(meta))
    candidate_rows[[length(candidate_rows)+1L]] <- data.frame(
      file = fn,
      sample_candidates = collapse_or_na(cand$sample),
      condition_candidates = collapse_or_na(cand$condition),
      cluster_candidates = collapse_or_na(cand$cluster),
      annotation_candidates = collapse_or_na(cand$annotation),
      reduction_candidates = collapse_or_na(reductions),
      stringsAsFactors = FALSE
    )
    object_summaries[[fn]] <- list(
      file = fn, full_path = path,
      dimensions = c(features = nrow(obj), cells = ncol(obj)),
      default_assay = get_default_assay_safe(obj), assays = assays,
      reductions = reductions, metadata_columns = colnames(meta),
      metadata_candidates = cand
    )
  }
  rm(obj); gc(verbose = FALSE)
}

bind_or_empty <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
safe_write_csv(bind_or_empty(inventory_rows), file.path(output_dir, "RDS_inventory_summary.csv"))
safe_write_csv(bind_or_empty(metadata_rows), file.path(output_dir, "RDS_metadata_columns.csv"))
safe_write_csv(bind_or_empty(layer_rows), file.path(output_dir, "RDS_assay_layers.csv"))
safe_write_csv(bind_or_empty(candidate_rows), file.path(output_dir, "RDS_detected_metadata_candidates.csv"))
saveRDS(object_summaries, file.path(output_dir, "RDS_inventory_object_summaries.rds"))
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Output: ", output_dir, "\n", sep = "")
