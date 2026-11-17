#!/usr/bin/env Rscript

# ==============================================================================
# Build_PrimaryClean_Mphi_Sets_v1.0.R
#
# PURPOSE
# -------
# Build analysis-ready Human and Mouse macrophage datasets after the completed
# clustering / annotation / contamination / ambient-RNA audit.
#
# IMPORTANT:
# - Source RDS files are NEVER overwritten.
# - Frozen cluster definitions are preserved.
# - NO PCA / neighbors / clustering / UMAP is recomputed.
# - DecontX-corrected counts remain a sensitivity-analysis resource only.
#
# FINAL PRIMARY QC POLICY
# -----------------------
# Human:
#   keep all frozen Mphi cells EXCEPT cluster 13
#   cluster 13 = coherent neutrophil/non-Mphi contamination
#
# Mouse:
#   keep all frozen Mphi cells EXCEPT cluster 23
#   cluster 23 = Mphi-LSEC heterotypic doublet / mixed-cell candidate
#   cluster 23 remains eligible for sensitivity analysis only.
#
# Frozen resolutions remain:
#   Human res0.8
#   Mouse res2.0
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

VERSION <- "v1.0"

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_IN <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

MOUSE_IN <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

HUMAN_OUTDIR <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Mphi_PrimaryClean_v2.0"
)

MOUSE_OUTDIR <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mouse_Mphi_PrimaryClean_v5.0"
)

CROSS_OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_Mphi",
  "PrimaryClean_Mphi_QC_Finalization_v1.0"
)

for (d in c(
  HUMAN_OUTDIR,
  file.path(HUMAN_OUTDIR, "objects"),
  file.path(HUMAN_OUTDIR, "tables"),
  MOUSE_OUTDIR,
  file.path(MOUSE_OUTDIR, "objects"),
  file.path(MOUSE_OUTDIR, "tables"),
  CROSS_OUTDIR,
  file.path(CROSS_OUTDIR, "tables"),
  file.path(CROSS_OUTDIR, "provenance"),
  file.path(CROSS_OUTDIR, "logs")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOGFILE <- file.path(CROSS_OUTDIR, "logs", "run.log")
cat("", file = LOGFILE)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOGFILE, append = TRUE)
}

write_tsv <- function(x, f) {
  write.table(
    x, f,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

find_first_col <- function(meta, candidates, label) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) {
    stop(
      "Could not find ", label, " column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

summarize_counts <- function(obj, species, cluster_col, sample_col, stage) {
  md <- obj@meta.data
  cl <- as.character(md[[cluster_col]])
  sp <- as.character(md[[sample_col]])

  by_cluster <- as.data.frame(table(cluster = cl), stringsAsFactors = FALSE)
  by_cluster$species <- species
  by_cluster$stage <- stage

  by_sample <- as.data.frame(table(sample = sp), stringsAsFactors = FALSE)
  by_sample$species <- species
  by_sample$stage <- stage

  by_sample_cluster <- as.data.frame(
    table(sample = sp, cluster = cl),
    stringsAsFactors = FALSE
  )
  by_sample_cluster <- by_sample_cluster[by_sample_cluster$Freq > 0, , drop = FALSE]
  by_sample_cluster$species <- species
  by_sample_cluster$stage <- stage

  list(
    cluster = by_cluster[, c("species", "stage", "cluster", "Freq")],
    sample = by_sample[, c("species", "stage", "sample", "Freq")],
    sample_cluster = by_sample_cluster[
      , c("species", "stage", "sample", "cluster", "Freq")
    ]
  )
}

audit_one <- function(
    obj,
    species,
    cluster_col,
    sample_col,
    excluded_cluster,
    exclusion_reason,
    sensitivity_keep_excluded = FALSE) {

  md <- obj@meta.data
  cl <- as.character(md[[cluster_col]])

  if (!excluded_cluster %in% cl) {
    stop(
      species,
      ": excluded cluster ",
      excluded_cluster,
      " not found in ",
      cluster_col
    )
  }

  # Preserve original frozen cluster explicitly.
  md$qc_original_frozen_cluster <- cl

  # Primary analysis inclusion.
  md$qc_primary_include <- cl != excluded_cluster

  # Sensitivity analysis inclusion.
  # Human c13 remains excluded even in standard sensitivity analysis.
  # Mouse c23 can be restored for sensitivity analyses.
  md$qc_sensitivity_include <- md$qc_primary_include
  if (isTRUE(sensitivity_keep_excluded)) {
    md$qc_sensitivity_include[cl == excluded_cluster] <- TRUE
  }

  md$qc_exclusion_reason <- "KEEP"
  md$qc_exclusion_reason[cl == excluded_cluster] <- exclusion_reason

  md$qc_final_disposition <- "PRIMARY_KEEP"
  md$qc_final_disposition[cl == excluded_cluster] <- if (species == "Human") {
    "EXCLUDE_NON_MPHI"
  } else {
    "EXCLUDE_PRIMARY_MIXED_DOUBLET"
  }

  md$qc_audit_version <- "Human_Mouse_Mphi_FinalQC_2026-11-18"
  md$qc_resolution_policy <- if (species == "Human") {
    "retain_res0.8_no_recluster"
  } else {
    "retain_res2.0_no_recluster"
  }

  md$qc_ambient_policy <- paste(
    "Original filtered counts remain primary;",
    "DecontX-corrected whole-liver counts are sensitivity analysis only"
  )

  obj@meta.data <- md

  keep_cells <- rownames(md)[md$qc_primary_include]
  clean <- subset(obj, cells = keep_cells)

  if (ncol(clean) != sum(md$qc_primary_include)) {
    stop(species, ": primary subset cell count mismatch.")
  }

  before <- summarize_counts(
    obj, species, cluster_col, sample_col, "before_QC_exclusion"
  )

  after <- summarize_counts(
    clean, species, cluster_col, sample_col, "primary_clean"
  )

  decision <- data.frame(
    species = species,
    cluster = sort(unique(cl)),
    n_cells = as.integer(table(factor(cl, levels = sort(unique(cl))))),
    primary_include = sort(unique(cl)) != excluded_cluster,
    sensitivity_include = if (sensitivity_keep_excluded) {
      rep(TRUE, length(sort(unique(cl))))
    } else {
      sort(unique(cl)) != excluded_cluster
    },
    final_disposition = ifelse(
      sort(unique(cl)) == excluded_cluster,
      ifelse(
        species == "Human",
        "EXCLUDE_NON_MPHI",
        "EXCLUDE_PRIMARY_MIXED_DOUBLET"
      ),
      "PRIMARY_KEEP"
    ),
    reason = ifelse(
      sort(unique(cl)) == excluded_cluster,
      exclusion_reason,
      "Passed final Mphi QC audit"
    ),
    stringsAsFactors = FALSE
  )

  list(
    flagged = obj,
    clean = clean,
    before = before,
    after = after,
    decision = decision
  )
}

# ------------------------------------------------------------------------------
# Load frozen source objects
# ------------------------------------------------------------------------------

for (f in c(HUMAN_IN, MOUSE_IN)) {
  if (!file.exists(f)) stop("Missing source RDS: ", f)
}

msg("Loading Human frozen Mphi object...")
human <- readRDS(HUMAN_IN)

msg("Loading Mouse frozen Mphi object...")
mouse <- readRDS(MOUSE_IN)

H_CLUSTER <- find_first_col(
  human@meta.data,
  c(
    "res0_8_cluster_frozen_v18",
    "RNA_snn_res.0.8",
    "myeloid_RPCA_res0_8"
  ),
  "Human frozen cluster"
)

M_CLUSTER <- find_first_col(
  mouse@meta.data,
  c(
    "cluster_res2",
    "mphi_rpca_res_2.0"
  ),
  "Mouse frozen cluster"
)

H_SAMPLE <- find_first_col(
  human@meta.data,
  c(
    "sample_v17",
    "orig.ident",
    "sample",
    "sample_display"
  ),
  "Human sample"
)

M_SAMPLE <- find_first_col(
  mouse@meta.data,
  c(
    "sample_display_FIXED2",
    "sample",
    "orig.ident",
    "sample_display"
  ),
  "Mouse sample"
)

msg("Human cluster column: ", H_CLUSTER)
msg("Human sample column:  ", H_SAMPLE)
msg("Mouse cluster column: ", M_CLUSTER)
msg("Mouse sample column:  ", M_SAMPLE)

# ------------------------------------------------------------------------------
# Apply final QC policy
# ------------------------------------------------------------------------------

ha <- audit_one(
  obj = human,
  species = "Human",
  cluster_col = H_CLUSTER,
  sample_col = H_SAMPLE,
  excluded_cluster = "13",
  exclusion_reason = paste(
    "Coherent neutrophil/non-Mphi contamination:",
    "CSF3R/MPO/ELANE/OLFM4 positive; weak macrophage core;",
    "not explained by ambient RNA"
  ),
  sensitivity_keep_excluded = FALSE
)

ma <- audit_one(
  obj = mouse,
  species = "Mouse",
  cluster_col = M_CLUSTER,
  sample_col = M_SAMPLE,
  excluded_cluster = "23",
  exclusion_reason = paste(
    "Mphi-LSEC heterotypic doublet/mixed-cell candidate:",
    "coherent endothelial/LSEC program and high complexity;",
    "excluded from primary Mphi biology but retained for sensitivity analyses"
  ),
  sensitivity_keep_excluded = TRUE
)

# ------------------------------------------------------------------------------
# Hard sanity checks
# ------------------------------------------------------------------------------

expected_human_removed <- sum(
  as.character(human@meta.data[[H_CLUSTER]]) == "13"
)

expected_mouse_removed <- sum(
  as.character(mouse@meta.data[[M_CLUSTER]]) == "23"
)

stopifnot(
  ncol(ha$clean) == ncol(human) - expected_human_removed,
  ncol(ma$clean) == ncol(mouse) - expected_mouse_removed
)

msg(
  "Human: ",
  ncol(human),
  " -> ",
  ncol(ha$clean),
  " (removed ",
  expected_human_removed,
  ")"
)

msg(
  "Mouse: ",
  ncol(mouse),
  " -> ",
  ncol(ma$clean),
  " (removed ",
  expected_mouse_removed,
  ")"
)

# ------------------------------------------------------------------------------
# Save objects
# ------------------------------------------------------------------------------

H_FLAGGED <- file.path(
  HUMAN_OUTDIR,
  "objects",
  "Human_MASH_Mphi_FinalQC_FLAGGED_v2.0.rds"
)

H_CLEAN <- file.path(
  HUMAN_OUTDIR,
  "objects",
  "Human_MASH_Mphi_PrimaryClean_v2.0.rds"
)

M_FLAGGED <- file.path(
  MOUSE_OUTDIR,
  "objects",
  "Mouse_MASH_Mphi_FinalQC_FLAGGED_v5.0.rds"
)

M_CLEAN <- file.path(
  MOUSE_OUTDIR,
  "objects",
  "Mouse_MASH_Mphi_PrimaryClean_v5.0.rds"
)

msg("Saving Human flagged object...")
saveRDS(ha$flagged, H_FLAGGED, compress = "xz")

msg("Saving Human primary-clean object...")
saveRDS(ha$clean, H_CLEAN, compress = "xz")

msg("Saving Mouse flagged object...")
saveRDS(ma$flagged, M_FLAGGED, compress = "xz")

msg("Saving Mouse primary-clean object...")
saveRDS(ma$clean, M_CLEAN, compress = "xz")

# ------------------------------------------------------------------------------
# Save audit tables
# ------------------------------------------------------------------------------

decision_all <- rbind(ha$decision, ma$decision)
write_tsv(
  decision_all,
  file.path(
    CROSS_OUTDIR,
    "tables",
    "01_final_cluster_disposition_v1.0.tsv"
  )
)

cluster_counts <- rbind(
  ha$before$cluster,
  ha$after$cluster,
  ma$before$cluster,
  ma$after$cluster
)

write_tsv(
  cluster_counts,
  file.path(
    CROSS_OUTDIR,
    "tables",
    "02_cluster_counts_before_after_primary_QC_v1.0.tsv"
  )
)

sample_counts <- rbind(
  ha$before$sample,
  ha$after$sample,
  ma$before$sample,
  ma$after$sample
)

write_tsv(
  sample_counts,
  file.path(
    CROSS_OUTDIR,
    "tables",
    "03_sample_counts_before_after_primary_QC_v1.0.tsv"
  )
)

sample_cluster_counts <- rbind(
  ha$before$sample_cluster,
  ha$after$sample_cluster,
  ma$before$sample_cluster,
  ma$after$sample_cluster
)

write_tsv(
  sample_cluster_counts,
  file.path(
    CROSS_OUTDIR,
    "tables",
    "04_sample_cluster_counts_before_after_primary_QC_v1.0.tsv"
  )
)

summary_table <- data.frame(
  species = c("Human", "Mouse"),
  source_cells = c(ncol(human), ncol(mouse)),
  excluded_cluster = c("13", "23"),
  excluded_cells = c(expected_human_removed, expected_mouse_removed),
  primary_clean_cells = c(ncol(ha$clean), ncol(ma$clean)),
  retained_fraction = c(
    ncol(ha$clean) / ncol(human),
    ncol(ma$clean) / ncol(mouse)
  ),
  frozen_resolution = c("res0.8", "res2.0"),
  reclustered = c(FALSE, FALSE),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  file.path(
    CROSS_OUTDIR,
    "tables",
    "05_primary_clean_summary_v1.0.tsv"
  )
)

# ------------------------------------------------------------------------------
# Provenance / policy
# ------------------------------------------------------------------------------

settings <- data.frame(
  item = c(
    "version",
    "human_source",
    "human_flagged_output",
    "human_primary_clean_output",
    "human_excluded_cluster",
    "human_resolution",
    "mouse_source",
    "mouse_flagged_output",
    "mouse_primary_clean_output",
    "mouse_excluded_cluster",
    "mouse_resolution",
    "reclustering_policy",
    "ambient_policy",
    "DecontX_policy"
  ),
  value = c(
    VERSION,
    HUMAN_IN,
    H_FLAGGED,
    H_CLEAN,
    "13",
    "retain res0.8",
    MOUSE_IN,
    M_FLAGGED,
    M_CLEAN,
    "23",
    "retain res2.0",
    "No PCA/neighbors/clustering/UMAP recomputation",
    "Do not exclude Human c3 or Mouse c8/c15/26 based on hepatocyte RNA alone",
    "Original filtered counts primary; DecontX counts sensitivity analysis only"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  settings,
  file.path(
    CROSS_OUTDIR,
    "provenance",
    "analysis_settings_v1.0.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    CROSS_OUTDIR,
    "provenance",
    "sessionInfo.txt"
  ),
  useBytes = TRUE
)

README <- c(
  "Primary clean Human / Mouse macrophage sets - final QC",
  "",
  "FINAL DECISION",
  "- Human frozen resolution: res0.8",
  "- Mouse frozen resolution: res2.0",
  "- Human cluster 13 excluded from primary Mphi analysis as coherent neutrophil/non-Mphi contamination.",
  "- Mouse cluster 23 excluded from primary Mphi analysis as Mphi-LSEC heterotypic doublet/mixed-cell candidate.",
  "- Mouse cluster 23 remains available in the flagged source-derived object for sensitivity analyses.",
  "- Human cluster 3 and Mouse clusters 8/15/26 are retained.",
  "",
  "AMBIENT RNA POLICY",
  "- Original filtered counts remain the primary expression matrix.",
  "- Whole-liver DecontX-corrected counts are sensitivity-analysis data only.",
  "- DecontX showed strong reduction of hepatocyte/LSEC contamination while preserving core Mphi programs.",
  "",
  "NO RECLUSTERING",
  "- Frozen cluster labels are preserved.",
  "- No PCA, neighbors, clustering or UMAP was recomputed in this step.",
  "",
  "SOURCE RDS FILES WERE NOT MODIFIED."
)

writeLines(
  README,
  file.path(CROSS_OUTDIR, "README.txt"),
  useBytes = TRUE
)

msg("================================================================")
msg("STATUS: SUCCESS")
msg("Human primary clean: ", H_CLEAN)
msg("Mouse primary clean: ", M_CLEAN)
msg("Review first:")
msg("  01_final_cluster_disposition_v1.0.tsv")
msg("  03_sample_counts_before_after_primary_QC_v1.0.tsv")
msg("  05_primary_clean_summary_v1.0.tsv")
msg("================================================================")
