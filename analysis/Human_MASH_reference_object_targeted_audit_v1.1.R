#!/usr/bin/env Rscript

# Human_MASH_reference_object_targeted_audit_v1.1.R
#
# READ-ONLY targeted audit.
# Purpose:
#   1) Resolve the correct sample-identity column in the current 12-sample human object.
#   2) Verify LT1-LT8 NAS mapping.
#   3) Compare LT1-LT8 cells between the older clustered2.0 object and the current
#      12-sample celltype-annotated RPCA object.
#   4) Audit RNA feature retention and annotation metadata.
#
# This script DOES NOT normalize, integrate, recluster, calculate UMAP, perform
# DE/pseudobulk/pathway/cell-cell analysis, modify metadata, or save any RDS.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.1"

OLD_RDS <- "/Users/uenoya/Desktop/clustered2.0_20251003.rds"
NEW_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/rds/",
  "Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples_celltype_annotated.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Reference_Object_Targeted_Audit_", VERSION
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(OLD_RDS))
stopifnot(file.exists(NEW_RDS))

suppressPackageStartupMessages({
  library(SeuratObject)
})

message("====================================================")
message("Human MASH reference-object targeted audit ", VERSION)
message("READ-ONLY")
message("====================================================")
message("OLD: ", OLD_RDS)
message("NEW: ", NEW_RDS)
message("OUT: ", OUTDIR)

old <- readRDS(OLD_RDS)
new <- readRDS(NEW_RDS)

# ---------- helpers ----------

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  u <- unique(x)
  if (!length(u)) return(NA_character_)
  paste(u, collapse = " | ")
}

get_layers <- function(obj, assay) {
  if (!assay %in% Assays(obj)) return(NA_character_)
  paste(Layers(obj[[assay]]), collapse = " | ")
}

audit_col <- function(md, nm) {
  if (!nm %in% colnames(md)) {
    return(data.frame(
      column = nm,
      exists = FALSE,
      n_missing = NA_integer_,
      n_unique_nonmissing = NA_integer_,
      values = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  x <- md[[nm]]
  data.frame(
    column = nm,
    exists = TRUE,
    n_missing = sum(is.na(x) | as.character(x) == ""),
    n_unique_nonmissing = length(unique(as.character(x)[!is.na(x) & as.character(x) != ""])),
    values = collapse_unique(x),
    stringsAsFactors = FALSE
  )
}

# ---------- 01 object comparison ----------

obj_compare <- data.frame(
  object_role = c("older_clustered2", "current_12sample_celltype_annotated"),
  path = c(OLD_RDS, NEW_RDS),
  cells = c(ncol(old), ncol(new)),
  features_default_assay = c(nrow(old), nrow(new)),
  default_assay = c(DefaultAssay(old), DefaultAssay(new)),
  assays = c(paste(Assays(old), collapse = " | "),
             paste(Assays(new), collapse = " | ")),
  RNA_layers = c(get_layers(old, "RNA"), get_layers(new, "RNA")),
  reductions = c(paste(Reductions(old), collapse = " | "),
                 paste(Reductions(new), collapse = " | ")),
  metadata_columns = c(ncol(old[[]]), ncol(new[[]])),
  stringsAsFactors = FALSE
)

write.csv(
  obj_compare,
  file.path(OUTDIR, "01_object_reference_comparison.csv"),
  row.names = FALSE
)

# ---------- 02 sample identity columns in new object ----------

new_md <- new[[]]

sample_cols <- c(
  "orig.ident",
  "sample_id_original",
  "sample_id",
  "sample_original_value",
  "sample",
  "canonical_sample",
  "sample_group"
)

sample_col_audit <- do.call(rbind, lapply(sample_cols, function(nm) {
  audit_col(new_md, nm)
}))

write.csv(
  sample_col_audit,
  file.path(OUTDIR, "02_new_sample_identity_columns.csv"),
  row.names = FALSE
)

if (!"canonical_sample" %in% colnames(new_md)) {
  stop("canonical_sample is absent from the current 12-sample object.")
}
if (any(is.na(new_md$canonical_sample) | as.character(new_md$canonical_sample) == "")) {
  stop("canonical_sample contains missing/blank values.")
}

canonical <- as.character(new_md$canonical_sample)

# Crosswalk each candidate column against canonical_sample.
crosswalk_list <- list()

for (nm in setdiff(sample_cols, c("canonical_sample", "sample_group"))) {
  if (!nm %in% colnames(new_md)) next
  x <- as.character(new_md[[nm]])
  x[is.na(x) | x == ""] <- "<NA/blank>"
  tab <- as.data.frame(
    table(
      canonical_sample = canonical,
      compared_column_value = x,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  tab$compared_column <- nm
  crosswalk_list[[length(crosswalk_list) + 1L]] <- tab[
    , c("compared_column", "canonical_sample", "compared_column_value", "Freq")
  ]
}

sample_crosswalk <- do.call(rbind, crosswalk_list)

write.csv(
  sample_crosswalk,
  file.path(OUTDIR, "03_new_sample_identity_crosswalk.csv"),
  row.names = FALSE
)

sample_counts <- as.data.frame(table(canonical_sample = canonical),
                               stringsAsFactors = FALSE)
colnames(sample_counts)[2] <- "n_cells"
sample_counts <- sample_counts[order(sample_counts$canonical_sample), ]

write.csv(
  sample_counts,
  file.path(OUTDIR, "04_new_sample_counts_canonical.csv"),
  row.names = FALSE
)

# ---------- 05 NAS mapping ----------

nas_cols <- intersect(c("NAS_score", "NAS_group", "sample_group"), colnames(new_md))

nas_rows <- lapply(sort(unique(canonical)), function(s) {
  idx <- canonical == s
  row <- data.frame(
    canonical_sample = s,
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
  for (nm in nas_cols) {
    row[[nm]] <- collapse_unique(new_md[[nm]][idx])
  }
  row
})

nas_map <- do.call(rbind, nas_rows)

write.csv(
  nas_map,
  file.path(OUTDIR, "05_NAS_mapping_by_canonical_sample.csv"),
  row.names = FALSE
)

# ---------- 06 LT1-LT8 cell overlap old vs new ----------

LT_SAMPLES <- paste0("LT", 1:8)

if (!"orig.ident" %in% colnames(old[[]])) {
  stop("Older clustered2.0 object lacks orig.ident.")
}

old_sample <- as.character(old$orig.ident)
new_sample <- canonical

overlap_rows <- lapply(LT_SAMPLES, function(s) {
  old_cells <- colnames(old)[old_sample == s]
  new_cells <- colnames(new)[new_sample == s]

  data.frame(
    sample = s,
    old_cells = length(old_cells),
    new_cells = length(new_cells),
    exact_colname_overlap = length(intersect(old_cells, new_cells)),
    old_only = length(setdiff(old_cells, new_cells)),
    new_only = length(setdiff(new_cells, old_cells)),
    exact_colname_sets_identical = setequal(old_cells, new_cells),
    stringsAsFactors = FALSE
  )
})

lt_overlap <- do.call(rbind, overlap_rows)

write.csv(
  lt_overlap,
  file.path(OUTDIR, "06_LT_cell_overlap_old_vs_new.csv"),
  row.names = FALSE
)

# If exact colnames are not identical, test the preserved original-cell-barcode field.
barcode_overlap <- data.frame()

if ("original_cell_barcode" %in% colnames(new_md)) {
  bc_rows <- lapply(LT_SAMPLES, function(s) {
    old_cells <- colnames(old)[old_sample == s]
    new_barcodes <- as.character(new_md$original_cell_barcode[new_sample == s])

    data.frame(
      sample = s,
      old_cells = length(old_cells),
      new_original_cell_barcodes = length(new_barcodes),
      old_colname_vs_new_original_barcode_overlap =
        length(intersect(old_cells, new_barcodes)),
      sets_identical = setequal(old_cells, new_barcodes),
      stringsAsFactors = FALSE
    )
  })
  barcode_overlap <- do.call(rbind, bc_rows)
}

write.csv(
  barcode_overlap,
  file.path(OUTDIR, "07_LT_original_barcode_overlap.csv"),
  row.names = FALSE
)

# ---------- 08 RNA feature overlap ----------

old_rna_features <- if ("RNA" %in% Assays(old)) rownames(old[["RNA"]]) else character(0)
new_rna_features <- if ("RNA" %in% Assays(new)) rownames(new[["RNA"]]) else character(0)

feature_summary <- data.frame(
  old_RNA_features = length(old_rna_features),
  new_RNA_features = length(new_rna_features),
  overlap = length(intersect(old_rna_features, new_rna_features)),
  old_only = length(setdiff(old_rna_features, new_rna_features)),
  new_only = length(setdiff(new_rna_features, old_rna_features)),
  new_subset_of_old = all(new_rna_features %in% old_rna_features),
  stringsAsFactors = FALSE
)

write.csv(
  feature_summary,
  file.path(OUTDIR, "08_RNA_feature_overlap_summary.csv"),
  row.names = FALSE
)

writeLines(
  sort(setdiff(old_rna_features, new_rna_features)),
  file.path(OUTDIR, "08b_RNA_features_old_only.txt")
)
writeLines(
  sort(setdiff(new_rna_features, old_rna_features)),
  file.path(OUTDIR, "08c_RNA_features_new_only.txt")
)

# ---------- 09 annotation metadata audit ----------

annotation_cols <- c(
  "auto_celltype",
  "auto_celltype_final",
  "celltype",
  "celltype_annotation_basis"
)

annotation_audit <- do.call(rbind, lapply(annotation_cols, function(nm) {
  audit_col(new_md, nm)
}))

write.csv(
  annotation_audit,
  file.path(OUTDIR, "09_annotation_metadata_audit.csv"),
  row.names = FALSE
)

annotation_pair_summary <- data.frame()

if (all(c("celltype", "auto_celltype_final") %in% colnames(new_md))) {
  a <- as.character(new_md$celltype)
  b <- as.character(new_md$auto_celltype_final)
  annotation_pair_summary <- data.frame(
    comparison = "celltype_vs_auto_celltype_final",
    n_cells = nrow(new_md),
    n_equal = sum(a == b, na.rm = TRUE),
    n_different = sum(a != b, na.rm = TRUE),
    n_missing_either = sum(is.na(a) | is.na(b)),
    stringsAsFactors = FALSE
  )
}

write.csv(
  annotation_pair_summary,
  file.path(OUTDIR, "10_annotation_pair_concordance.csv"),
  row.names = FALSE
)

# ---------- summary ----------

sample_id_missing_by_canonical <- character(0)
if ("sample_id" %in% colnames(new_md)) {
  sid <- as.character(new_md$sample_id)
  miss <- is.na(sid) | sid == ""
  if (any(miss)) {
    tt <- table(canonical[miss])
    sample_id_missing_by_canonical <- paste(
      paste0(names(tt), "=", as.integer(tt)),
      collapse = " | "
    )
  }
}

summary_lines <- c(
  paste0("Human MASH reference-object targeted audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no RDS object was modified or saved.",
  "",
  "OBJECTS:",
  paste0("OLD: ", OLD_RDS),
  paste0("NEW: ", NEW_RDS),
  "",
  "CURRENT 12-SAMPLE SAMPLE IDENTITY:",
  paste0("canonical_sample missing: ",
         sum(is.na(new_md$canonical_sample) | canonical == "")),
  paste0("canonical_sample unique: ", length(unique(canonical))),
  paste0("canonical_sample values: ",
         paste(sort(unique(canonical)), collapse = " | ")),
  paste0("sample_id missing by canonical_sample: ",
         ifelse(length(sample_id_missing_by_canonical),
                sample_id_missing_by_canonical, "<none>")),
  "",
  "LT1-LT8 EXACT CELL OVERLAP:",
  paste(capture.output(print(lt_overlap, row.names = FALSE)), collapse = "\n"),
  "",
  "RNA FEATURE OVERLAP:",
  paste(capture.output(print(feature_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "DECISION GATE:",
  "- Do not start biological analysis until canonical sample identity and LT cell mapping are confirmed.",
  "- Do not assume the 9,137-feature annotated RNA assay is the optimal source for pseudobulk until feature provenance is resolved.",
  "",
  "OUTPUT FILES:",
  "01_object_reference_comparison.csv",
  "02_new_sample_identity_columns.csv",
  "03_new_sample_identity_crosswalk.csv",
  "04_new_sample_counts_canonical.csv",
  "05_NAS_mapping_by_canonical_sample.csv",
  "06_LT_cell_overlap_old_vs_new.csv",
  "07_LT_original_barcode_overlap.csv",
  "08_RNA_feature_overlap_summary.csv",
  "08b_RNA_features_old_only.txt",
  "08c_RNA_features_new_only.txt",
  "09_annotation_metadata_audit.csv",
  "10_annotation_pair_concordance.csv",
  "11_audit_summary.txt"
)

writeLines(
  summary_lines,
  file.path(OUTDIR, "11_audit_summary.txt"),
  useBytes = TRUE
)

message("====================================================")
message("TARGETED AUDIT COMPLETE")
message("Output: ", OUTDIR)
message("No RDS object was modified or saved.")
message("====================================================")
