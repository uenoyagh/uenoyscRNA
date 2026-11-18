#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_PrimaryClean_Finalization_v1.0.R
#
# PURPOSE
# -------
# Finalize HSC QC dispositions after:
#   - provenance / resolution audit
#   - cluster-level foreign-lineage audit
#   - sample-wise whole-liver filtered-count DecontX sensitivity analysis
#
# FINAL DISPOSITION POLICY
# ------------------------
# Primary analysis:
#   KEEP all frozen res0.8 HSC clusters EXCEPT c14, c15, c16
#
# Sensitivity analysis:
#   c15 is retained as a sensitivity-only HSC-LSEC mixed/doublet candidate
#   c14 and c16 remain excluded
#
# Cluster rationale:
#   c14: coherent B-cell program persists after DecontX
#   c15: coherent LSEC/endothelial program persists after DecontX, but HSC
#        identity is also retained -> sensitivity-only
#   c16: strong striated/skeletal-muscle program persists/increases after
#        DecontX with marked loss of canonical HSC identity
#
# IMPORTANT
# ---------
# - Frozen working resolution remains res0.8.
# - No reclustering.
# - Original filtered counts remain PRIMARY.
# - DecontX-corrected counts remain sensitivity-only and are NOT substituted
#   into the canonical PrimaryClean HSC RDS.
# - Source RDS is NEVER overwritten.
#
# Expected primary-clean cell count:
#   13,166 - 96 - 85 - 44 = 12,941 cells
#
# Version: v1.0.1
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2",
  "RDS",
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

OUTDIR <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0"
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
PROVDIR <- file.path(OUTDIR, "provenance")

for (d in c(OUTDIR, OBJDIR, TABDIR, PROVDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
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

if (!file.exists(HSC_RDS)) {
  stop("Missing source HSC RDS: ", HSC_RDS)
}

message("Loading source HSC object...")
obj <- readRDS(HSC_RDS)
md <- obj@meta.data

required <- c(
  "hsc_cluster_working_v610",
  "hsc_rpca_res_0_8_v610",
  "sample_hsc_v610",
  "HSC_state3_v612",
  "HSC_state2_v612"
)

missing_cols <- setdiff(required, colnames(md))
if (length(missing_cols)) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

working <- as.character(md$hsc_cluster_working_v610)
res08 <- as.character(md$hsc_rpca_res_0_8_v610)

if (!identical(working, res08)) {
  stop("Working HSC cluster is not exactly identical to stored res0.8.")
}

if (ncol(obj) != 13166) {
  stop("Unexpected source cell count: ", ncol(obj), " (expected 13166)")
}

# ------------------------------------------------------------------------------
# Final QC disposition
# ------------------------------------------------------------------------------

PRIMARY_EXCLUDE <- c("14", "15", "16")
SENSITIVITY_EXCLUDE <- c("14", "16")

reason <- rep(NA_character_, ncol(obj))

reason[working == "14"] <- paste(
  "Coherent B-cell program persists after sample-wise whole-liver DecontX;",
  "Ebf1/Cd19/Cd79a/Cd79b/Pax5/Ighm remain strongly enriched;",
  "non-HSC mixed/doublet candidate"
)

reason[working == "15"] <- paste(
  "Coherent LSEC/endothelial program persists after DecontX;",
  "Egfl7/Clec4g/Kdr/Pecam1 and other endothelial markers remain enriched;",
  "HSC identity is also retained; HSC-LSEC mixed/doublet candidate"
)

reason[working == "16"] <- paste(
  "Strong striated/skeletal-muscle program persists or increases after DecontX;",
  "Tnnt3/Ttn/Neb/Actn2/Ryr1/Stac3/Pygm dominate while canonical HSC identity is weak;",
  "non-HSC candidate"
)

primary_include <- !(working %in% PRIMARY_EXCLUDE)
sensitivity_include <- !(working %in% SENSITIVITY_EXCLUDE)

disposition <- rep("KEEP_PRIMARY", ncol(obj))
disposition[working == "14"] <- "EXCLUDE_PRIMARY_AND_SENSITIVITY"
disposition[working == "15"] <- "EXCLUDE_PRIMARY_KEEP_SENSITIVITY"
disposition[working == "16"] <- "EXCLUDE_PRIMARY_AND_SENSITIVITY"

obj$qc_original_frozen_cluster <- working
obj$qc_primary_include <- primary_include
obj$qc_sensitivity_include <- sensitivity_include
obj$qc_exclusion_reason <- reason
obj$qc_final_disposition <- disposition
obj$qc_audit_version <- "Mouse_HSC_FinalQC_v1.0"
obj$qc_resolution_policy <- paste(
  "Frozen HSC working solution = stored res0.8;",
  "hsc_cluster_working_v610 exactly identical to hsc_rpca_res_0_8_v610;",
  "no reclustering after QC"
)
obj$qc_ambient_policy <- paste(
  "Original filtered counts are primary;",
  "sample-wise whole-liver filtered-count DecontX used for sensitivity only;",
  "raw/empty-droplet background unavailable/not used"
)

# ------------------------------------------------------------------------------
# Save all-cell audited object
# ------------------------------------------------------------------------------

AUDITED_RDS <- file.path(
  OBJDIR,
  "Mouse_MASH_HSC_QCAudited_v7.0.rds"
)

saveRDS(obj, AUDITED_RDS, compress = FALSE)

# ------------------------------------------------------------------------------
# Primary Clean HSC
# ------------------------------------------------------------------------------

primary_cells <- colnames(obj)[obj$qc_primary_include %in% TRUE]
primary <- subset(obj, cells = primary_cells)

PRIMARY_RDS <- file.path(
  OBJDIR,
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

saveRDS(primary, PRIMARY_RDS, compress = FALSE)

# ------------------------------------------------------------------------------
# Sensitivity HSC
# ------------------------------------------------------------------------------

sensitivity_cells <- colnames(obj)[obj$qc_sensitivity_include %in% TRUE]
sensitivity <- subset(obj, cells = sensitivity_cells)

SENS_RDS <- file.path(
  OBJDIR,
  "Mouse_MASH_HSC_SensitivityClean_v7.0.rds"
)

saveRDS(sensitivity, SENS_RDS, compress = FALSE)

# ------------------------------------------------------------------------------
# Verification tables
# ------------------------------------------------------------------------------

disposition_summary <- as.data.frame(
  table(
    cluster = working,
    disposition = disposition
  ),
  stringsAsFactors = FALSE
)
disposition_summary <- disposition_summary[
  disposition_summary$Freq > 0,
  ,
  drop = FALSE
]

write_tsv(
  disposition_summary,
  file.path(TABDIR, "01_cluster_final_disposition.tsv")
)

source_sample <- as.character(obj$sample_hsc_v610)

sample_before <- as.data.frame(
  table(sample = source_sample),
  stringsAsFactors = FALSE
)
names(sample_before)[2] <- "before_QC"

sample_primary <- as.data.frame(
  table(sample = as.character(primary$sample_hsc_v610)),
  stringsAsFactors = FALSE
)
names(sample_primary)[2] <- "primary_clean"

sample_sensitivity <- as.data.frame(
  table(sample = as.character(sensitivity$sample_hsc_v610)),
  stringsAsFactors = FALSE
)
names(sample_sensitivity)[2] <- "sensitivity_clean"

sample_counts <- Reduce(
  function(x, y) merge(x, y, by = "sample", all = TRUE),
  list(sample_before, sample_primary, sample_sensitivity)
)

for (cn in setdiff(colnames(sample_counts), "sample")) {
  sample_counts[[cn]][is.na(sample_counts[[cn]])] <- 0
}

sample_counts$primary_removed <- (
  sample_counts$before_QC - sample_counts$primary_clean
)
sample_counts$sensitivity_removed <- (
  sample_counts$before_QC - sample_counts$sensitivity_clean
)

write_tsv(
  sample_counts,
  file.path(TABDIR, "02_sample_counts_before_after_QC.tsv")
)

state3_before <- as.data.frame(
  table(
    sample = as.character(obj$sample_hsc_v610),
    HSC_state3 = as.character(obj$HSC_state3_v612)
  ),
  stringsAsFactors = FALSE
)

state3_primary <- as.data.frame(
  table(
    sample = as.character(primary$sample_hsc_v610),
    HSC_state3 = as.character(primary$HSC_state3_v612)
  ),
  stringsAsFactors = FALSE
)

state3_primary <- state3_primary[state3_primary$Freq > 0, , drop = FALSE]

write_tsv(
  state3_before[state3_before$Freq > 0, , drop = FALSE],
  file.path(TABDIR, "03_HSC3state_before_QC.tsv")
)

write_tsv(
  state3_primary,
  file.path(TABDIR, "04_HSC3state_primary_clean.tsv")
)

state2_primary <- as.data.frame(
  table(
    sample = as.character(primary$sample_hsc_v610),
    HSC_state2 = as.character(primary$HSC_state2_v612)
  ),
  stringsAsFactors = FALSE
)

state2_primary <- state2_primary[state2_primary$Freq > 0, , drop = FALSE]

write_tsv(
  state2_primary,
  file.path(TABDIR, "05_HSC2state_primary_clean.tsv")
)

# Exact flagged-cell inventory
flagged <- data.frame(
  cell = colnames(obj)[working %in% PRIMARY_EXCLUDE],
  sample = as.character(obj$sample_hsc_v610[working %in% PRIMARY_EXCLUDE]),
  cluster = working[working %in% PRIMARY_EXCLUDE],
  HSC_state3 = as.character(obj$HSC_state3_v612[working %in% PRIMARY_EXCLUDE]),
  HSC_state2 = as.character(obj$HSC_state2_v612[working %in% PRIMARY_EXCLUDE]),
  qc_primary_include = primary_include[working %in% PRIMARY_EXCLUDE],
  qc_sensitivity_include = sensitivity_include[working %in% PRIMARY_EXCLUDE],
  qc_final_disposition = disposition[working %in% PRIMARY_EXCLUDE],
  qc_exclusion_reason = reason[working %in% PRIMARY_EXCLUDE],
  stringsAsFactors = FALSE
)

write_tsv(
  flagged,
  file.path(TABDIR, "06_flagged_cells_inventory.tsv")
)

# ------------------------------------------------------------------------------
# Final integrity checks
# ------------------------------------------------------------------------------

EXPECTED_PRIMARY <- 12941L
EXPECTED_SENSITIVITY <- 13026L

checks <- data.frame(
  item = c(
    "source_cells",
    "primary_clean_cells",
    "expected_primary_clean_cells",
    "primary_count_matches_expected",
    "sensitivity_clean_cells",
    "expected_sensitivity_clean_cells",
    "sensitivity_count_matches_expected",
    "working_equals_res0.8",
    "primary_all_qc_primary_include_TRUE",
    "sensitivity_all_qc_sensitivity_include_TRUE",
    "source_RDS_overwritten",
    "DecontX_counts_used_as_primary"
  ),
  value = c(
    ncol(obj),
    ncol(primary),
    EXPECTED_PRIMARY,
    as.character(ncol(primary) == EXPECTED_PRIMARY),
    ncol(sensitivity),
    EXPECTED_SENSITIVITY,
    as.character(ncol(sensitivity) == EXPECTED_SENSITIVITY),
    as.character(identical(working, res08)),
    as.character(all(primary$qc_primary_include %in% TRUE)),
    as.character(all(sensitivity$qc_sensitivity_include %in% TRUE)),
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  checks,
  file.path(TABDIR, "07_final_integrity_checks.tsv")
)

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

settings <- data.frame(
  item = c(
    "version",
    "source_RDS",
    "audited_all_cell_RDS",
    "primary_clean_RDS",
    "sensitivity_clean_RDS",
    "frozen_resolution",
    "primary_excluded_clusters",
    "sensitivity_excluded_clusters",
    "ambient_policy",
    "reclustering_policy"
  ),
  value = c(
    "v1.0.1",
    HSC_RDS,
    AUDITED_RDS,
    PRIMARY_RDS,
    SENS_RDS,
    "res0.8",
    paste(PRIMARY_EXCLUDE, collapse = ";"),
    paste(SENSITIVITY_EXCLUDE, collapse = ";"),
    "Original filtered counts primary; whole-liver DecontX sensitivity only",
    "No reclustering"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  settings,
  file.path(PROVDIR, "analysis_settings_v1.0.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(PROVDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Source cells: ", ncol(obj))
message("Primary Clean cells: ", ncol(primary))
message("Sensitivity Clean cells: ", ncol(sensitivity))
message("Primary excluded clusters: ", paste(PRIMARY_EXCLUDE, collapse = ", "))
message("Sensitivity excluded clusters: ", paste(SENSITIVITY_EXCLUDE, collapse = ", "))
message("Primary RDS:")
message("  ", PRIMARY_RDS)
message("Sensitivity RDS:")
message("  ", SENS_RDS)
message("Review:")
message("  01_cluster_final_disposition.tsv")
message("  02_sample_counts_before_after_QC.tsv")
message("  04_HSC3state_primary_clean.tsv")
message("  07_final_integrity_checks.tsv")
message("============================================================")
