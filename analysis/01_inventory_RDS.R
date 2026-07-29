rm(list=ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"
source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "seurat_helpers.R"))
source(file.path(script_root, "R", "metadata.R"))
source(file.path(script_root, "R", "utils.R"))

if (!requireNamespace("SeuratObject", quietly=TRUE)) {
  stop("Required package is missing: SeuratObject")
}

target_dir <- get_dataset_dir(analysis_target)
output_dir <- get_result_dir(analysis_target, "inventory", TRUE)
log_dir <- get_result_dir(analysis_target, "logs", TRUE)
rds_files <- list_rds_files(target_dir)

if (!length(rds_files)) stop("No RDS files found in: ", target_dir)

inventory_rows <- metadata_rows <- layer_rows <- candidate_rows <- list()
object_summaries <- list()

log_file <- file.path(log_dir, paste0("inventory_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split=TRUE)
on.exit(while (sink.number() > 0) sink(), add=TRUE)

cat("Started: ", timestamp_string(), "\n", sep="")
cat("Target: ", analysis_target, "\n", sep="")

for (i in seq_along(rds_files)) {
  path <- rds_files[[i]]
  fn <- basename(path)
  cat("[", i, "/", length(rds_files), "] ", fn, "\n", sep="")
  obj <- safe_read_rds(path)

  if (inherits(obj, "rds_read_error")) {
    inventory_rows[[length(inventory_rows)+1]] <- data.frame(
      file=fn, full_path=path, file_size_gb=file_size_gb(path),
      read_status="ERROR", object_class=NA, seurat_object=FALSE,
      n_cells=NA, n_features=NA, default_assay=NA,
      assays=NA, reductions=NA, metadata_columns=NA,
      error_message=obj$message, stringsAsFactors=FALSE)
    next
  }

  seu <- is_seurat_object(obj)
  meta <- if (seu) obj[[]] else NULL
  assays <- if (seu) get_assay_names(obj) else character(0)
  reductions <- if (seu) get_reduction_names(obj) else character(0)

  inventory_rows[[length(inventory_rows)+1]] <- data.frame(
    file=fn, full_path=path, file_size_gb=file_size_gb(path),
    read_status="OK", object_class=paste(class(obj), collapse="/"),
    seurat_object=seu,
    n_cells=if (seu) ncol(obj) else NA,
    n_features=if (seu) nrow(obj) else NA,
    default_assay=if (seu) get_default_assay_safe(obj) else NA,
    assays=collapse_or_na(assays),
    reductions=collapse_or_na(reductions),
    metadata_columns=if (seu) ncol(meta) else NA,
    error_message=NA, stringsAsFactors=FALSE)

  if (seu) {
    z <- get_layer_names(obj)
    if (nrow(z)) {
      z$file <- fn
      layer_rows[[length(layer_rows)+1]] <- z[,c("file","assay","layer")]
    }

    m <- metadata_column_summary(meta, inventory_max_unique_values)
    m$file <- fn
    metadata_rows[[length(metadata_rows)+1]] <- m[,c(
      "file","column","class","n_non_na","n_na","n_unique","example_values")]

    cand <- detect_metadata_candidates(colnames(meta))
    candidate_rows[[length(candidate_rows)+1]] <- data.frame(
      file=fn,
      sample_candidates=collapse_or_na(cand$sample),
      condition_candidates=collapse_or_na(cand$condition),
      cluster_candidates=collapse_or_na(cand$cluster),
      annotation_candidates=collapse_or_na(cand$annotation),
      reduction_candidates=collapse_or_na(reductions),
      stringsAsFactors=FALSE)

    object_summaries[[fn]] <- list(
      file=fn, dimensions=c(features=nrow(obj), cells=ncol(obj)),
      default_assay=get_default_assay_safe(obj),
      assays=assays, reductions=reductions,
      metadata_columns=colnames(meta),
      metadata_candidates=cand
    )
  }

  rm(obj)
  gc(verbose=FALSE)
}

bind_or_empty <- function(x) if (length(x)) do.call(rbind, x) else data.frame()

safe_write_csv(bind_or_empty(inventory_rows),
               file.path(output_dir, "RDS_inventory_summary.csv"))
safe_write_csv(bind_or_empty(metadata_rows),
               file.path(output_dir, "RDS_metadata_columns.csv"))
safe_write_csv(bind_or_empty(layer_rows),
               file.path(output_dir, "RDS_assay_layers.csv"))
safe_write_csv(bind_or_empty(candidate_rows),
               file.path(output_dir, "RDS_detected_metadata_candidates.csv"))

if (isTRUE(inventory_save_rds_summary)) {
  saveRDS(object_summaries,
          file.path(output_dir, "RDS_inventory_object_summaries.rds"))
}

writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))

cat("Completed: ", timestamp_string(), "\n", sep="")
cat("Output: ", output_dir, "\n", sep="")
