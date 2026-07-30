#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("analysis/RDS3_annotation_validation_v3.2/99_run_all.R", mustWork = FALSE)
}
analysis_dir <- dirname(script_path)
project_root_guess <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = FALSE)

source(file.path(project_root_guess, "config", "config.R"))
source(file.path(project_root_guess, "R", "01_utils_v3.2.R"))
source(file.path(project_root_guess, "R", "02_markers_v3.2.R"))
source(file.path(project_root_guess, "R", "03_analysis_functions_v3.2.R"))
source(file.path(project_root_guess, "R", "04_plot_export_v3.2.R"))

required <- c(
  "Seurat", "SeuratObject", "Matrix", "presto", "future",
  "future.apply", "ggplot2", "patchwork", "dplyr", "tibble",
  "scales", "openxlsx", "here"
)
assert_packages(required)

project_root <- resolve_project_root(CFG$project_root)
external_root <- detect_external_root(CFG$external_ssd_candidates)
rds_path <- file.path(external_root, CFG$rds_relative_path)
output_dir <- file.path(external_root, CFG$output_relative_dir)

if (!file.exists(rds_path)) {
  stop(
    "RDS was not found on the external SSD:\n", rds_path,
    "\nEdit CFG$rds_relative_path in config/config.R."
  )
}

dirs <- list(
  output = output_dir,
  checkpoint = file.path(output_dir, "checkpoints"),
  tables = file.path(output_dir, "tables"),
  figures = file.path(output_dir, "figures"),
  logs = file.path(output_dir, "logs"),
  objects = file.path(output_dir, "objects")
)
ensure_dirs(unname(dirs))

ctx <- list(
  cfg = CFG,
  project_root = project_root,
  external_root = external_root,
  rds_path = rds_path,
  output_dir = output_dir,
  checkpoint_dir = dirs$checkpoint
)

log_file <- file.path(dirs$logs, "annotation_validation_v3.2.log")
cat("", file = log_file)

log_msg("Pipeline: ", CFG$pipeline_version, log_file = log_file)
log_msg("Project root: ", project_root, log_file = log_file)
log_msg("External root: ", external_root, log_file = log_file)
log_msg("Input RDS: ", rds_path, log_file = log_file)
log_msg("Output: ", output_dir, log_file = log_file)

if (isTRUE(CFG$use_future)) {
  options(future.globals.maxSize = CFG$future_max_size_gb * 1024^3)
  future::plan(future::multisession, workers = CFG$workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  log_msg("future multisession workers: ", CFG$workers, log_file = log_file)
} else {
  future::plan(future::sequential)
}

# STEP 1: load and validate
if (should_run(ctx, "01_object_prepared")) {
  log_msg("STEP 1: reading RDS from external SSD.", log_file = log_file)
  object <- readRDS(rds_path)
  if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.")
  object <- prepare_seurat_for_markers(object, CFG, log_file)

  meta <- object[[]]
  CFG$cluster_col <- detect_metadata_column(meta, CFG$cluster_col, CFG$cluster_candidates, "cluster")
  CFG$annotation_col <- detect_metadata_column(meta, CFG$annotation_col, CFG$annotation_candidates, "annotation")
  CFG$sample_col <- tryCatch(
    detect_metadata_column(meta, CFG$sample_col, CFG$sample_candidates, "sample"),
    error = function(e) NULL
  )
  CFG$condition_col <- tryCatch(
    detect_metadata_column(meta, CFG$condition_col, CFG$condition_candidates, "condition"),
    error = function(e) NULL
  )

  if (!CFG$reduction %in% names(object@reductions)) {
    stop("Reduction '", CFG$reduction, "' was not found.")
  }

  object[[CFG$cluster_col]][, 1] <- factor(object[[]][[CFG$cluster_col]])
  save_checkpoint(list(object = object, resolved_cfg = CFG), ctx, "01_object_prepared")
} else {
  tmp <- load_checkpoint(ctx, "01_object_prepared")
  object <- tmp$object
  CFG <- tmp$resolved_cfg
  ctx$cfg <- CFG
  log_msg("STEP 1 restored from checkpoint.", log_file = log_file)
}

safe_write_csv(
  data.frame(
    field = c("rds_path","output_dir","cluster_col","annotation_col","sample_col","condition_col",
              "n_cells","n_features"),
    value = c(rds_path, output_dir, CFG$cluster_col, CFG$annotation_col,
              CFG$sample_col %||% NA, CFG$condition_col %||% NA,
              ncol(object), nrow(object))
  ),
  file.path(dirs$tables, "run_manifest.csv")
)

# STEP 2: markers
if (should_run(ctx, "02_cluster_markers")) {
  log_msg("STEP 2: cluster marker analysis.", log_file = log_file)
  markers <- if (identical(tolower(CFG$marker_engine), "presto")) {
    run_presto_markers(object, CFG$cluster_col, CFG, log_file)
  } else {
    run_seurat_markers(object, CFG$cluster_col, CFG, log_file)
  }
  markers_top <- summarize_top_markers(markers, CFG$top_n_markers)
  save_checkpoint(list(all = markers, top = markers_top), ctx, "02_cluster_markers")
} else {
  marker_ck <- load_checkpoint(ctx, "02_cluster_markers")
  markers <- marker_ck$all
  markers_top <- marker_ck$top
  log_msg("STEP 2 restored from checkpoint.", log_file = log_file)
}
safe_write_csv(markers, file.path(dirs$tables, "cluster_markers_all.csv"))
safe_write_csv(markers_top, file.path(dirs$tables, "cluster_markers_top.csv"))

# STEP 3: marker voting
marker_reference <- get_marker_reference()
safe_write_csv(marker_reference, file.path(dirs$tables, "marker_reference_used.csv"))

if (should_run(ctx, "03_marker_voting")) {
  log_msg("STEP 3: General/Ueno marker voting.", log_file = log_file)
  voting <- run_marker_voting(object, CFG$cluster_col, marker_reference, CFG)
  object <- attach_voting_metadata(object, CFG$cluster_col, voting$winners)
  save_checkpoint(list(voting = voting, object = object), ctx, "03_marker_voting")
} else {
  vote_ck <- load_checkpoint(ctx, "03_marker_voting")
  voting <- vote_ck$voting
  object <- vote_ck$object
  log_msg("STEP 3 restored from checkpoint.", log_file = log_file)
}

safe_write_csv(voting$scores, file.path(dirs$tables, "marker_voting_all_scores.csv"))
safe_write_csv(voting$winners, file.path(dirs$tables, "marker_voting_winners.csv"))

# Annotation agreement table
cluster_annotation <- as.data.frame(table(
  cluster = as.character(object[[]][[CFG$cluster_col]]),
  current_annotation = as.character(object[[]][[CFG$annotation_col]])
))
cluster_annotation <- cluster_annotation[cluster_annotation$Freq > 0, , drop = FALSE]
cluster_annotation <- do.call(rbind, lapply(split(cluster_annotation, cluster_annotation$cluster), function(x) {
  x <- x[order(-x$Freq), , drop = FALSE]
  x$cluster_total <- sum(x$Freq)
  x$fraction <- x$Freq / x$cluster_total
  x
}))
safe_write_csv(cluster_annotation, file.path(dirs$tables, "cluster_current_annotation_composition.csv"))

# STEP 4: plots
log_msg("STEP 4: UMAP, DotPlot and Violin export.", log_file = log_file)

umap_fields <- unique(na.omit(c(
  CFG$cluster_col,
  CFG$annotation_col,
  "vote_general_label",
  "vote_ueno_label",
  CFG$sample_col,
  CFG$condition_col
)))
for (field in umap_fields) {
  if (!field %in% colnames(object[[]])) next
  p <- publish_umap_v32(object, field, CFG, title = paste("RDS3:", field))
  save_pdf(
    p,
    file.path(dirs$figures, paste0("UMAP_", gsub("[^A-Za-z0-9_]+", "_", field), ".pdf")),
    CFG$pdf_width, CFG$pdf_height
  )
}

for (src in c("General", "Ueno")) {
  p_dot <- make_marker_dotplot(object, marker_reference, CFG$cluster_col, CFG, source = src)
  save_pdf(
    p_dot,
    file.path(dirs$figures, paste0("DotPlot_", src, "_markers_by_cluster.pdf")),
    18, 14
  )

  pages <- make_marker_violin_pages(object, marker_reference, CFG$cluster_col, CFG, source = src)
  pdf_path <- file.path(dirs$figures, paste0("Violin_", src, "_markers_by_cluster.pdf"))
  grDevices::cairo_pdf(pdf_path, width = 16, height = 12, onefile = TRUE)
  for (p in pages) print(p)
  grDevices::dev.off()
}

# STEP 5: Excel and final object
log_msg("STEP 5: Excel report and final checkpoint.", log_file = log_file)

write_excel_report(
  file.path(dirs$tables, "RDS3_annotation_validation_v3.2.xlsx"),
  list(
    Run_manifest = data.frame(
      Item = c("Pipeline","Input_RDS","External_root","Cluster_column","Annotation_column",
               "Cells","Features","Marker_engine","Workers"),
      Value = c(CFG$pipeline_version, rds_path, external_root, CFG$cluster_col,
                CFG$annotation_col, ncol(object), nrow(object),
                CFG$marker_engine, CFG$workers)
    ),
    Cluster_markers_top = markers_top,
    Voting_winners = voting$winners,
    Voting_all_scores = voting$scores,
    Current_annotation = cluster_annotation,
    Marker_reference = marker_reference,
    Marker_mean_expression = voting$expression
  )
)

if (isTRUE(CFG$save_intermediate_rds)) {
  saveRDS(
    object,
    file.path(dirs$objects, "RDS3_with_annotation_validation_v3.2.rds"),
    compress = FALSE
  )
}
save_checkpoint(list(completed = TRUE, finished_at = timestamp()), ctx, "99_complete")
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))

log_msg("All v3.2 steps completed successfully.", log_file = log_file)
log_msg("Results: ", output_dir, log_file = log_file)
