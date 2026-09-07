#!/usr/bin/env Rscript

# Human_MASH_LT_cellmap_feature_provenance_audit_v1.2.R
#
# READ-ONLY targeted audit.
#
# Purpose:
#   1) Resolve LT1-LT8 cell identity between the older clustered2.0 object
#      and the current 12-sample celltype-annotated object despite renamed cells.
#   2) Verify RNA count equality for matched LT cells across the 9,137 common genes.
#   3) Test whether the current 9,137-gene RNA universe is explained by the
#      intersection of the older LT object and the later 2-sample raw object.
#   4) Audit retention/loss of key macrophage/MASH marker genes.
#
# This script DOES NOT:
#   - normalize / SCTransform
#   - integrate
#   - cluster / recluster
#   - calculate UMAP
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#   - modify metadata
#   - save or overwrite any RDS
#
# Seurat v5 convention:
#   GetAssayData(..., layer = "counts") only; no slot= usage.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.2"

OLD_RDS <- "/Users/uenoya/Desktop/clustered2.0_20251003.rds"

CURRENT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/rds/",
  "Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples_celltype_annotated.rds"
)

TX_RAW_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/rds/",
  "Human_RH251217117_2samples_raw_SeuratObject_fromMTX.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/LT_CellMap_FeatureProvenance_Audit_", VERSION
)

for (p in c(OLD_RDS, CURRENT_RDS, TX_RAW_RDS)) {
  if (!file.exists(p)) stop("Required RDS not found:\n", p)
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

message("============================================================")
message("Human MASH LT cell-map + feature-provenance audit ", VERSION)
message("READ-ONLY")
message("============================================================")

extract_barcode_core <- function(x) {
  x <- as.character(x)
  m <- regexpr("[ACGT]+-[0-9]+", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)
  out
}

safe_get_counts <- function(obj) {
  if (!"RNA" %in% Assays(obj)) stop("RNA assay is absent.")
  GetAssayData(obj, assay = "RNA", layer = "counts")
}

collapse_head <- function(x, n = 5L) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  paste(utils::head(x, n), collapse = " | ")
}

LT_SAMPLES <- paste0("LT", 1:8)

# ------------------------------------------------------------------
# Load OLD + CURRENT together because exact matched-count comparison
# requires both objects in memory.
# ------------------------------------------------------------------

message("Reading OLD object:")
message(OLD_RDS)
old <- readRDS(OLD_RDS)

message("Reading CURRENT annotated 12-sample object:")
message(CURRENT_RDS)
cur <- readRDS(CURRENT_RDS)

old_md <- old[[]]
cur_md <- cur[[]]

if (!"orig.ident" %in% colnames(old_md)) {
  stop("OLD object lacks orig.ident.")
}
if (!"canonical_sample" %in% colnames(cur_md)) {
  stop("CURRENT object lacks canonical_sample.")
}
if (!"original_cell_barcode" %in% colnames(cur_md)) {
  stop("CURRENT object lacks original_cell_barcode.")
}

old_sample <- as.character(old_md$orig.ident)
cur_sample <- as.character(cur_md$canonical_sample)

# ---------- 01 naming examples ----------

naming_examples <- do.call(rbind, lapply(LT_SAMPLES, function(s) {
  old_cells <- colnames(old)[old_sample == s]
  cur_cells <- colnames(cur)[cur_sample == s]
  cur_orig <- as.character(cur_md$original_cell_barcode[cur_sample == s])

  n <- min(5L, length(old_cells), length(cur_cells), length(cur_orig))
  data.frame(
    sample = s,
    index = seq_len(n),
    old_colname = old_cells[seq_len(n)],
    old_barcode_core = extract_barcode_core(old_cells[seq_len(n)]),
    current_colname = cur_cells[seq_len(n)],
    current_colname_barcode_core = extract_barcode_core(cur_cells[seq_len(n)]),
    current_original_cell_barcode = cur_orig[seq_len(n)],
    current_original_barcode_core = extract_barcode_core(cur_orig[seq_len(n)]),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  naming_examples,
  file.path(OUTDIR, "01_LT_cell_naming_examples.csv"),
  row.names = FALSE
)

# ---------- 02 barcode-core overlap ----------

barcode_overlap_rows <- list()

for (s in LT_SAMPLES) {
  old_cells <- colnames(old)[old_sample == s]
  cur_cells <- colnames(cur)[cur_sample == s]
  cur_orig <- as.character(cur_md$original_cell_barcode[cur_sample == s])

  old_core <- extract_barcode_core(old_cells)
  cur_col_core <- extract_barcode_core(cur_cells)
  cur_orig_core <- extract_barcode_core(cur_orig)

  barcode_overlap_rows[[length(barcode_overlap_rows) + 1L]] <- data.frame(
    sample = s,
    old_cells = length(old_cells),
    current_cells = length(cur_cells),

    old_core_nonmissing = sum(!is.na(old_core)),
    current_col_core_nonmissing = sum(!is.na(cur_col_core)),
    current_original_core_nonmissing = sum(!is.na(cur_orig_core)),

    old_vs_current_col_core_overlap =
      length(intersect(old_core[!is.na(old_core)], cur_col_core[!is.na(cur_col_core)])),

    old_vs_current_original_core_overlap =
      length(intersect(old_core[!is.na(old_core)], cur_orig_core[!is.na(cur_orig_core)])),

    old_vs_current_original_core_sets_identical =
      setequal(old_core[!is.na(old_core)], cur_orig_core[!is.na(cur_orig_core)]),

    old_core_unique = !anyDuplicated(old_core[!is.na(old_core)]),
    current_original_core_unique = !anyDuplicated(cur_orig_core[!is.na(cur_orig_core)]),

    stringsAsFactors = FALSE
  )
}

barcode_overlap <- do.call(rbind, barcode_overlap_rows)

write.csv(
  barcode_overlap,
  file.path(OUTDIR, "02_LT_barcode_core_overlap.csv"),
  row.names = FALSE
)

all_lt_barcode_resolved <- all(
  barcode_overlap$old_vs_current_original_core_sets_identical &
    barcode_overlap$old_core_nonmissing == barcode_overlap$old_cells &
    barcode_overlap$current_original_core_nonmissing == barcode_overlap$current_cells &
    barcode_overlap$old_core_unique &
    barcode_overlap$current_original_core_unique
)

# ---------- 03 explicit LT cell map ----------

cell_map <- data.frame()

if (all_lt_barcode_resolved) {
  maps <- list()

  for (s in LT_SAMPLES) {
    old_cells <- colnames(old)[old_sample == s]
    cur_cells <- colnames(cur)[cur_sample == s]
    cur_orig <- as.character(cur_md$original_cell_barcode[cur_sample == s])

    old_core <- extract_barcode_core(old_cells)
    cur_core <- extract_barcode_core(cur_orig)

    idx <- match(old_core, cur_core)

    if (anyNA(idx)) {
      stop("Unexpected unmatched barcode after set-equality check in ", s)
    }

    tmp <- data.frame(
      sample = s,
      barcode_core = old_core,
      old_cell = old_cells,
      current_cell = cur_cells[idx],
      current_original_cell_barcode = cur_orig[idx],
      stringsAsFactors = FALSE
    )

    if ("celltype" %in% colnames(cur_md)) {
      tmp$celltype <- as.character(cur_md[tmp$current_cell, "celltype"])
    }
    if ("NAS_score" %in% colnames(cur_md)) {
      tmp$NAS_score <- cur_md[tmp$current_cell, "NAS_score"]
    }

    maps[[length(maps) + 1L]] <- tmp
  }

  cell_map <- do.call(rbind, maps)
}

write.csv(
  cell_map,
  file.path(OUTDIR, "03_LT_old_to_current_cell_map.csv"),
  row.names = FALSE
)

# ---------- 04 common-feature count equality ----------

old_features <- rownames(old[["RNA"]])
cur_features <- rownames(cur[["RNA"]])
common_features <- intersect(old_features, cur_features)

count_compare_rows <- list()

if (all_lt_barcode_resolved) {
  message("Retrieving RNA counts for exact matched-cell comparison...")
  old_counts <- safe_get_counts(old)
  cur_counts <- safe_get_counts(cur)

  # Compare sample-by-sample to limit temporary sparse-matrix size.
  for (s in LT_SAMPLES) {
    message("Count comparison: ", s)

    mp <- cell_map[cell_map$sample == s, , drop = FALSE]

    a <- old_counts[common_features, mp$old_cell, drop = FALSE]
    b <- cur_counts[common_features, mp$current_cell, drop = FALSE]

    # Ignore column names only after explicit barcode-based alignment.
    colnames(a) <- NULL
    colnames(b) <- NULL

    d <- a - b
    n_diff <- Matrix::nnzero(d)
    max_abs_diff <- if (length(d@x)) max(abs(d@x)) else 0

    count_compare_rows[[length(count_compare_rows) + 1L]] <- data.frame(
      sample = s,
      matched_cells = ncol(a),
      common_features = nrow(a),
      differing_matrix_entries = n_diff,
      max_absolute_count_difference = max_abs_diff,
      counts_exactly_identical = (n_diff == 0),
      stringsAsFactors = FALSE
    )

    rm(a, b, d)
    gc(verbose = FALSE)
  }

  rm(old_counts, cur_counts)
  gc(verbose = FALSE)
}

count_compare <- if (length(count_compare_rows)) {
  do.call(rbind, count_compare_rows)
} else {
  data.frame()
}

write.csv(
  count_compare,
  file.path(OUTDIR, "04_LT_common_gene_count_identity.csv"),
  row.names = FALSE
)

# ---------- 05 current-vs-old feature relationship ----------

old_current_feature_summary <- data.frame(
  old_RNA_features = length(old_features),
  current_RNA_features = length(cur_features),
  old_current_overlap = length(intersect(old_features, cur_features)),
  old_only = length(setdiff(old_features, cur_features)),
  current_only = length(setdiff(cur_features, old_features)),
  current_subset_of_old = all(cur_features %in% old_features),
  stringsAsFactors = FALSE
)

write.csv(
  old_current_feature_summary,
  file.path(OUTDIR, "05_old_current_feature_relationship.csv"),
  row.names = FALSE
)

# We no longer need the two large objects for feature-provenance work.
rm(old, cur, old_md, cur_md)
gc(verbose = FALSE)

# ------------------------------------------------------------------
# Load later 2-sample raw source and test feature-set provenance.
# ------------------------------------------------------------------

message("Reading later 2-sample raw object:")
message(TX_RAW_RDS)
txraw <- readRDS(TX_RAW_RDS)

tx_features <- rownames(txraw[["RNA"]])

feature_provenance <- data.frame(
  old_features = length(old_features),
  tx_raw_features = length(tx_features),
  current_features = length(cur_features),

  old_intersect_tx_raw = length(intersect(old_features, tx_features)),
  old_intersect_tx_raw_equals_current =
    setequal(intersect(old_features, tx_features), cur_features),

  current_subset_of_old = all(cur_features %in% old_features),
  current_subset_of_tx_raw = all(cur_features %in% tx_features),

  tx_raw_only_vs_current = length(setdiff(tx_features, cur_features)),
  old_only_vs_current = length(setdiff(old_features, cur_features)),

  stringsAsFactors = FALSE
)

write.csv(
  feature_provenance,
  file.path(OUTDIR, "06_feature_provenance_old_vs_txraw_vs_current.csv"),
  row.names = FALSE
)

writeLines(
  sort(setdiff(tx_features, cur_features)),
  file.path(OUTDIR, "06b_TXraw_features_not_in_current.txt")
)

# ---------- 07 key macrophage/MASH marker retention ----------

key_markers <- c(
  "CD163", "MRC1", "FOLR2", "TIMD4", "CLEC4F", "MARCO", "VSIG4",
  "MERTK", "CD5L", "ADGRE1", "CD68",
  "SPP1", "TREM2", "CD9",
  "IL1B", "TNF", "CCR2", "SPN",
  "F13A1", "TGFB1", "VCAN", "POSTN", "THBS1", "FN1",
  "LGALS3", "ANXA1", "IL10", "IL10RA", "IL10RB",
  "LPL", "FABP5", "APOE",
  "C1QA", "C1QB", "C1QC",
  "MSR1", "OLR1", "CD36"
)

marker_retention <- data.frame(
  gene = key_markers,
  in_old_RNA = key_markers %in% old_features,
  in_TX_raw_RNA = key_markers %in% tx_features,
  in_current_RNA = key_markers %in% cur_features,
  lost_from_old_in_current =
    (key_markers %in% old_features) & !(key_markers %in% cur_features),
  stringsAsFactors = FALSE
)

write.csv(
  marker_retention,
  file.path(OUTDIR, "07_key_macrophage_marker_retention.csv"),
  row.names = FALSE
)

rm(txraw)
gc(verbose = FALSE)

# ---------- summary ----------

mapping_status <- if (all_lt_barcode_resolved) "RESOLVED" else "NOT_RESOLVED"

count_status <- if (nrow(count_compare) &&
                    all(count_compare$counts_exactly_identical)) {
  "EXACTLY_IDENTICAL_FOR_ALL_LT1-LT8_ON_COMMON_FEATURES"
} else if (!nrow(count_compare)) {
  "NOT_TESTED_BECAUSE_CELL_MAPPING_NOT_RESOLVED"
} else {
  "DIFFERENCES_DETECTED"
}

lost_key <- marker_retention$gene[marker_retention$lost_from_old_in_current]

summary_lines <- c(
  paste0("Human MASH LT cell-map + feature-provenance audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no RDS object was modified or saved.",
  "",
  paste0("LT1-LT8 barcode mapping status: ", mapping_status),
  paste0("LT1-LT8 common-gene raw-count status: ", count_status),
  "",
  "BARCODE CORE OVERLAP:",
  paste(capture.output(print(barcode_overlap, row.names = FALSE)), collapse = "\n"),
  "",
  "FEATURE PROVENANCE:",
  paste(capture.output(print(feature_provenance, row.names = FALSE)), collapse = "\n"),
  "",
  "KEY MΦ/MASH MARKERS LOST FROM OLD RNA IN CURRENT RNA:",
  if (length(lost_key)) paste(lost_key, collapse = " | ") else "<none>",
  "",
  "DECISION LOGIC:",
  "- If LT barcode mapping is RESOLVED and common-gene counts are exactly identical,",
  "  current cell annotations can be transferred back to the older LT1-LT8 count universe without ambiguity.",
  "- If old ∩ TX-raw feature set exactly equals the current feature set, the 9,137-gene",
  "  universe is explained by feature intersection during cross-source object construction.",
  "- Do not start pseudobulk/DE from the 9,137-gene current RNA assay if key macrophage",
  "  markers are absent and a wider matched LT count universe is available.",
  "",
  "OUTPUT FILES:",
  "01_LT_cell_naming_examples.csv",
  "02_LT_barcode_core_overlap.csv",
  "03_LT_old_to_current_cell_map.csv",
  "04_LT_common_gene_count_identity.csv",
  "05_old_current_feature_relationship.csv",
  "06_feature_provenance_old_vs_txraw_vs_current.csv",
  "06b_TXraw_features_not_in_current.txt",
  "07_key_macrophage_marker_retention.csv",
  "08_audit_summary.txt"
)

writeLines(
  summary_lines,
  file.path(OUTDIR, "08_audit_summary.txt"),
  useBytes = TRUE
)

message("============================================================")
message("AUDIT COMPLETE")
message("LT barcode mapping: ", mapping_status)
message("LT count identity: ", count_status)
message("Output: ", OUTDIR)
message("No RDS object was modified or saved.")
message("============================================================")
