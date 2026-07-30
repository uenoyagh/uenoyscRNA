#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("analysis/RDS3_annotation_validation_v3.4/99_run_all.R", mustWork = FALSE)
}
analysis_dir <- dirname(script_path)
project_root_guess <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = FALSE)

source(file.path(project_root_guess, "config", "config.R"))
source(file.path(project_root_guess, "R", "01_utils_v3.2.R"))
source(file.path(project_root_guess, "R", "02_markers_v3.4.R"))
source(file.path(project_root_guess, "R", "03_analysis_functions_v3.2.R"))
source(file.path(project_root_guess, "R", "03_hierarchical_voting_v3.3.R"))
source(file.path(project_root_guess, "R", "03_evidence_voting_v3.4.R"))
source(file.path(project_root_guess, "R", "04_plot_export_v3.2.R"))

CFG$pipeline_version <- "RDS3_annotation_validation_v3.4"
CFG$output_relative_dir <- file.path(
  dirname(CFG$output_relative_dir),
  "RDS3_annotation_validation_v3.4"
)

# Evidence thresholds.
CFG$v34_min_auc <- 0.55
CFG$v34_min_logfc <- 0.15
CFG$v34_min_pct_diff <- 5
CFG$v34_negative_weight <- 1.25
CFG$v34_gate_penalty <- 2.0
CFG$v34_prior_bonus <- 0.20
CFG$v34_unknown_score <- 0.50
CFG$v34_moderate_delta <- 0.35
CFG$v34_high_delta <- 0.80

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
v33_dir <- file.path(dirname(output_dir), "RDS3_annotation_validation_v3.3")

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

log_file <- file.path(dirs$logs, "annotation_validation_v3.4.log")
cat("", file = log_file)
log_msg("Pipeline: ", CFG$pipeline_version, log_file = log_file)
log_msg("Input RDS: ", rds_path, log_file = log_file)
log_msg("Output: ", output_dir, log_file = log_file)

# STEP 1: object
if (should_run(ctx, "01_object_prepared")) {
  object <- readRDS(rds_path)
  if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.")
  object <- prepare_seurat_for_markers(object, CFG, log_file)

  meta <- object[[]]
  CFG$cluster_col <- detect_metadata_column(
    meta, CFG$cluster_col, CFG$cluster_candidates, "cluster"
  )
  CFG$annotation_col <- detect_metadata_column(
    meta, CFG$annotation_col, CFG$annotation_candidates, "annotation"
  )
  CFG$sample_col <- tryCatch(
    detect_metadata_column(meta, CFG$sample_col, CFG$sample_candidates, "sample"),
    error = function(e) NULL
  )
  CFG$condition_col <- tryCatch(
    detect_metadata_column(meta, CFG$condition_col, CFG$condition_candidates, "condition"),
    error = function(e) NULL
  )

  save_checkpoint(list(object = object, resolved_cfg = CFG), ctx, "01_object_prepared")
} else {
  z <- load_checkpoint(ctx, "01_object_prepared")
  object <- z$object
  CFG <- z$resolved_cfg
  ctx$cfg <- CFG
}

# STEP 2: reuse v3.3 presto output when possible.
v33_markers_csv <- file.path(v33_dir, "tables", "cluster_markers_all.csv")
if (file.exists(v33_markers_csv)) {
  log_msg("Reusing v3.3 cluster_markers_all.csv.", log_file = log_file)
  markers <- utils::read.csv(v33_markers_csv, check.names = FALSE)
} else if (should_run(ctx, "02_cluster_markers")) {
  log_msg("v3.3 marker table unavailable; running presto.", log_file = log_file)
  markers <- run_presto_markers(object, CFG$cluster_col, CFG, log_file)
  save_checkpoint(markers, ctx, "02_cluster_markers")
} else {
  markers <- load_checkpoint(ctx, "02_cluster_markers")
}
safe_write_csv(markers, file.path(dirs$tables, "cluster_markers_all.csv"))

# STEP 3: evidence-weighted voting.
marker_reference <- get_marker_reference_v34()
safe_write_csv(
  marker_reference,
  file.path(dirs$tables, "marker_reference_used_v3.4.csv")
)

if (should_run(ctx, "03_evidence_voting")) {
  log_msg("Running evidence-weighted hierarchical voting.", log_file = log_file)
  voting <- run_evidence_voting_v34(
    object,
    CFG$cluster_col,
    CFG$annotation_col,
    markers,
    marker_reference,
    CFG
  )
  object <- attach_evidence_metadata_v34(object, CFG$cluster_col, voting)
  save_checkpoint(
    list(voting = voting, object = object),
    ctx,
    "03_evidence_voting"
  )
} else {
  z <- load_checkpoint(ctx, "03_evidence_voting")
  voting <- z$voting
  object <- z$object
}

audit <- make_annotation_audit_v34(voting)

safe_write_csv(
  voting$marker_evidence,
  file.path(dirs$tables, "marker_evidence_v3.4.csv")
)
safe_write_csv(
  voting$scores,
  file.path(dirs$tables, "evidence_voting_all_scores_v3.4.csv")
)
safe_write_csv(
  voting$winners,
  file.path(dirs$tables, "evidence_voting_winners_v3.4.csv")
)
safe_write_csv(
  audit,
  file.path(dirs$tables, "annotation_audit_v3.4.csv")
)

# STEP 4: UMAPs
umap_fields <- c(
  CFG$annotation_col,
  "vote_ueno_lineage_v34",
  "vote_ueno_celltype_v34",
  "vote_ueno_subtype_v34",
  "vote_ueno_recommended_v34",
  "vote_ueno_lineage_v34_confidence",
  "vote_ueno_celltype_v34_confidence",
  "vote_ueno_subtype_v34_confidence"
)

for (field in umap_fields) {
  if (!field %in% colnames(object[[]])) next
  p <- publish_umap_v32(
    object, field, CFG, title = paste("RDS3 v3.4:", field)
  )
  save_pdf(
    p,
    file.path(
      dirs$figures,
      paste0("UMAP_", gsub("[^A-Za-z0-9_]+", "_", field), ".pdf")
    ),
    CFG$pdf_width,
    CFG$pdf_height
  )
}

# STEP 5: report and object
write_excel_report(
  file.path(dirs$tables, "RDS3_annotation_validation_v3.4.xlsx"),
  list(
    Run_manifest = data.frame(
      Item = c(
        "Pipeline","Input_RDS","Cluster_column","Annotation_column",
        "Cells","Features","Marker_source"
      ),
      Value = c(
        CFG$pipeline_version, rds_path, CFG$cluster_col, CFG$annotation_col,
        ncol(object), nrow(object),
        ifelse(file.exists(v33_markers_csv), "v3.3 reused", "recomputed")
      )
    ),
    Annotation_audit = audit,
    Voting_winners = voting$winners,
    Voting_all_scores = voting$scores,
    Marker_evidence = voting$marker_evidence,
    Current_annotation = voting$current_annotation,
    Marker_reference = marker_reference
  )
)

saveRDS(
  object,
  file.path(dirs$objects, "RDS3_with_annotation_validation_v3.4.rds"),
  compress = FALSE
)

save_checkpoint(list(completed = TRUE, finished_at = timestamp()), ctx, "99_complete")
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))
log_msg("All v3.4 steps completed successfully.", log_file = log_file)
