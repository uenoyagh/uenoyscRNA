#!/usr/bin/env Rscript

# Human_MASH_LT1to8_build_canonical_object_v1.3.R
#
# PURPOSE
# -------
# Build a NEW LT1-LT8 canonical human MASH Seurat object using:
#   - OLD clustered2.0_20251003.rds as the 16,424-gene RNA/count source
#   - CURRENT 12-sample celltype-annotated RPCA object as the source of
#     canonical sample identity, NAS metadata, cell-type metadata, and frozen
#     current umapRPCA coordinates.
#
# SAFETY / NON-PURPOSE
# --------------------
# This script DOES NOT:
#   - modify or overwrite either source RDS
#   - normalize / SCTransform
#   - integrate
#   - cluster / recluster
#   - calculate a new UMAP
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#   - declare the transferred celltype annotation biologically final
#
# The existing 12-sample umapRPCA coordinates are TRANSFERRED EXACTLY
# by barcode mapping; they are not recalculated.
#
# Seurat v5 convention:
#   GetAssayData(..., layer = "counts") only; no slot= usage.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.3"

OLD_RDS <- "/Users/uenoya/Desktop/clustered2.0_20251003.rds"

REF12_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/rds/",
  "Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples_celltype_annotated.rds"
)

BASE_OUT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_", VERSION
)

OBJDIR <- file.path(BASE_OUT, "objects")
TABDIR <- file.path(BASE_OUT, "tables")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_Canonical_", VERSION, ".rds")
)

LT_SAMPLES <- paste0("LT", 1:8)

EXPECTED_NAS <- c(
  LT1 = 6, LT2 = 6,
  LT3 = 4, LT4 = 4,
  LT5 = 2, LT6 = 2,
  LT7 = 0, LT8 = 0
)

# Current-reference metadata that must exist before construction.
REQUIRED_REF12_META <- c(
  "canonical_sample",
  "original_cell_barcode",
  "NAS_score",
  "celltype"
)

# Preserve the full current-reference metadata under a ref12_ prefix.
# This avoids overwriting any metadata already present in clustered2.0.
TRANSFER_ALL_REF12_METADATA <- TRUE

for (p in c(OLD_RDS, REF12_RDS)) {
  if (!file.exists(p)) {
    stop("Required source RDS not found:\n", p)
  }
}

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SeuratObject)
})

message("============================================================")
message("Human MASH LT1-LT8 canonical object build ", VERSION)
message("============================================================")
message("OLD count source: ", OLD_RDS)
message("REF12 metadata/UMAP source: ", REF12_RDS)
message("Output RDS: ", OUT_RDS)

extract_barcode_core <- function(x) {
  x <- as.character(x)
  m <- regexpr("[ACGT]+-[0-9]+", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)
  out
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  u <- unique(x)
  if (!length(u)) return(NA_character_)
  paste(u, collapse = " | ")
}

message("Reading OLD 16,424-gene object...")
old <- readRDS(OLD_RDS)

message("Reading REF12 annotated object...")
ref12 <- readRDS(REF12_RDS)

if (!inherits(old, "Seurat")) stop("OLD_RDS is not a Seurat object.")
if (!inherits(ref12, "Seurat")) stop("REF12_RDS is not a Seurat object.")

if (!"RNA" %in% Assays(old)) stop("OLD object lacks RNA assay.")
if (!"RNA" %in% Assays(ref12)) stop("REF12 object lacks RNA assay.")
if (!"umapRPCA" %in% Reductions(ref12)) {
  stop("REF12 object lacks reduction 'umapRPCA'.")
}

old_md <- old[[]]
ref_md <- ref12[[]]

missing_ref_meta <- setdiff(REQUIRED_REF12_META, colnames(ref_md))
if (length(missing_ref_meta)) {
  stop(
    "REF12 object is missing required metadata: ",
    paste(missing_ref_meta, collapse = ", ")
  )
}

if (!"orig.ident" %in% colnames(old_md)) {
  stop("OLD object lacks orig.ident.")
}

# ------------------------------------------------------------------
# 1. Identify LT1-LT8 cells and reconstruct exact 1:1 barcode map.
# ------------------------------------------------------------------

old_sample <- as.character(old_md$orig.ident)
ref_sample <- as.character(ref_md$canonical_sample)

old_lt_cells <- colnames(old)[old_sample %in% LT_SAMPLES]
ref_lt_cells <- colnames(ref12)[ref_sample %in% LT_SAMPLES]

if (length(old_lt_cells) != 81348L) {
  warning(
    "OLD LT1-LT8 total is ", length(old_lt_cells),
    " rather than the previously audited 81,348. Continuing only if all per-sample checks pass."
  )
}

map_list <- list()
map_summary_list <- list()

for (s in LT_SAMPLES) {
  old_cells_s <- colnames(old)[old_sample == s]
  ref_cells_s <- colnames(ref12)[ref_sample == s]
  ref_orig_s <- as.character(ref_md$original_cell_barcode[ref_sample == s])

  old_core <- extract_barcode_core(old_cells_s)
  ref_core <- extract_barcode_core(ref_orig_s)

  if (anyNA(old_core)) {
    stop("Could not extract barcode core from OLD cells in ", s)
  }
  if (anyNA(ref_core)) {
    stop("Could not extract barcode core from REF12 original_cell_barcode in ", s)
  }
  if (anyDuplicated(old_core)) {
    stop("Duplicated OLD barcode core detected in ", s)
  }
  if (anyDuplicated(ref_core)) {
    stop("Duplicated REF12 barcode core detected in ", s)
  }
  if (!setequal(old_core, ref_core)) {
    stop("Barcode-core sets do not match between OLD and REF12 for ", s)
  }

  idx <- match(old_core, ref_core)
  if (anyNA(idx)) stop("Unexpected unmatched barcode in ", s)

  map_s <- data.frame(
    canonical_sample = s,
    barcode_core = old_core,
    old_cell_id = old_cells_s,
    ref12_cell_id = ref_cells_s[idx],
    ref12_original_cell_barcode = ref_orig_s[idx],
    stringsAsFactors = FALSE
  )

  map_list[[s]] <- map_s

  map_summary_list[[s]] <- data.frame(
    canonical_sample = s,
    old_cells = length(old_cells_s),
    ref12_cells = length(ref_cells_s),
    matched_cells = nrow(map_s),
    mapping_complete = nrow(map_s) == length(old_cells_s) &&
      nrow(map_s) == length(ref_cells_s),
    stringsAsFactors = FALSE
  )
}

cell_map <- do.call(rbind, map_list)
rownames(cell_map) <- NULL

map_summary <- do.call(rbind, map_summary_list)
rownames(map_summary) <- NULL

if (nrow(cell_map) != length(old_lt_cells)) {
  stop("Final LT mapping row count does not equal OLD LT1-LT8 cell count.")
}
if (anyDuplicated(cell_map$old_cell_id)) {
  stop("Duplicate OLD cell ID in final cell map.")
}
if (anyDuplicated(cell_map$ref12_cell_id)) {
  stop("Duplicate REF12 cell ID in final cell map.")
}

write.csv(
  map_summary,
  file.path(TABDIR, "01_cell_mapping_summary_v1.3.csv"),
  row.names = FALSE
)

write.csv(
  cell_map,
  file.path(TABDIR, "02_LT_old_to_ref12_cell_map_v1.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 2. Validate NAS mapping before building.
# ------------------------------------------------------------------

nas_validation_list <- list()

for (s in LT_SAMPLES) {
  ref_cells_s <- cell_map$ref12_cell_id[cell_map$canonical_sample == s]
  nas_values <- unique(ref_md[ref_cells_s, "NAS_score"])

  nas_values_nonmissing <- nas_values[!is.na(nas_values)]
  if (length(nas_values_nonmissing) != 1L) {
    stop(
      "Expected exactly one non-missing NAS_score for ", s,
      "; observed: ", paste(nas_values, collapse = " | ")
    )
  }

  observed_nas <- as.numeric(as.character(nas_values_nonmissing[[1]]))
  expected_nas <- unname(EXPECTED_NAS[[s]])

  if (is.na(observed_nas) || observed_nas != expected_nas) {
    stop(
      "NAS mismatch for ", s,
      ": observed=", observed_nas,
      ", expected=", expected_nas
    )
  }

  nas_validation_list[[s]] <- data.frame(
    canonical_sample = s,
    n_cells = length(ref_cells_s),
    NAS_score = observed_nas,
    expected_NAS = expected_nas,
    NAS_validated = TRUE,
    stringsAsFactors = FALSE
  )
}

nas_validation <- do.call(rbind, nas_validation_list)
rownames(nas_validation) <- NULL

write.csv(
  nas_validation,
  file.path(TABDIR, "03_sample_NAS_validation_v1.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 3. Subset OLD object to LT1-LT8 ONLY.
#    This is the source of the wider 16,424-gene RNA universe.
# ------------------------------------------------------------------

message("Subsetting OLD object to LT1-LT8...")
canonical <- subset(old, cells = old_lt_cells)

if (ncol(canonical) != nrow(cell_map)) {
  stop("Canonical subset cell number mismatch.")
}

old_feature_n <- nrow(old[["RNA"]])
canonical_feature_n <- nrow(canonical[["RNA"]])

if (canonical_feature_n != old_feature_n) {
  stop(
    "Feature loss occurred during cell-only subset: OLD=", old_feature_n,
    ", canonical=", canonical_feature_n
  )
}

if (canonical_feature_n != 16424L) {
  warning(
    "Canonical RNA feature count is ", canonical_feature_n,
    " rather than the previously audited 16,424."
  )
}

# ------------------------------------------------------------------
# 4. Transfer REF12 metadata by exact barcode map.
# ------------------------------------------------------------------

map_for_old <- cell_map[match(colnames(canonical), cell_map$old_cell_id), ]

if (anyNA(map_for_old$ref12_cell_id)) {
  stop("Could not align cell map to canonical object cell order.")
}

ref_meta_aligned <- ref_md[map_for_old$ref12_cell_id, , drop = FALSE]

if (nrow(ref_meta_aligned) != ncol(canonical)) {
  stop("Aligned REF12 metadata row count mismatch.")
}

# Preserve ALL reference metadata under a prefix, without overwriting OLD metadata.
if (TRANSFER_ALL_REF12_METADATA) {
  ref_meta_prefixed <- ref_meta_aligned
  colnames(ref_meta_prefixed) <- paste0("ref12_", colnames(ref_meta_prefixed))
  rownames(ref_meta_prefixed) <- colnames(canonical)

  duplicate_prefixed <- intersect(
    colnames(ref_meta_prefixed),
    colnames(canonical[[]])
  )
  if (length(duplicate_prefixed)) {
    stop(
      "Prefixed REF12 metadata would overwrite existing columns: ",
      paste(duplicate_prefixed, collapse = ", ")
    )
  }

  canonical <- AddMetaData(canonical, metadata = ref_meta_prefixed)
}

# Canonical standardized metadata.
std_meta <- data.frame(
  canonical_sample_v13 =
    as.character(ref_meta_aligned$canonical_sample),

  NAS_score_v13 =
    as.numeric(as.character(ref_meta_aligned$NAS_score)),

  NAS_group_v13 =
    paste0("NAS", as.numeric(as.character(ref_meta_aligned$NAS_score))),

  celltype_ref12_v13 =
    as.character(ref_meta_aligned$celltype),

  barcode_core_v13 =
    map_for_old$barcode_core,

  source_old_cell_id_v13 =
    map_for_old$old_cell_id,

  ref12_cell_id_v13 =
    map_for_old$ref12_cell_id,

  ref12_original_cell_barcode_v13 =
    map_for_old$ref12_original_cell_barcode,

  canonical_object_version_v13 =
    VERSION,

  count_source_v13 =
    "clustered2.0_20251003.rds",

  annotation_source_v13 =
    "Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples_celltype_annotated.rds",

  stringsAsFactors = FALSE,
  row.names = colnames(canonical)
)

if ("auto_celltype_final" %in% colnames(ref_meta_aligned)) {
  std_meta$auto_celltype_final_ref12_v13 <-
    as.character(ref_meta_aligned$auto_celltype_final)
}

if ("celltype_annotation_basis" %in% colnames(ref_meta_aligned)) {
  std_meta$celltype_annotation_basis_ref12_v13 <-
    as.character(ref_meta_aligned$celltype_annotation_basis)
}

canonical <- AddMetaData(canonical, metadata = std_meta)

# Validate standardized sample and NAS columns.
if (anyNA(canonical$canonical_sample_v13)) {
  stop("canonical_sample_v13 contains NA.")
}
if (anyNA(canonical$NAS_score_v13)) {
  stop("NAS_score_v13 contains NA within LT1-LT8.")
}
if (!setequal(unique(as.character(canonical$canonical_sample_v13)), LT_SAMPLES)) {
  stop("canonical_sample_v13 does not contain exactly LT1-LT8.")
}

# ------------------------------------------------------------------
# 5. Transfer the CURRENT frozen umapRPCA coordinates EXACTLY.
#    No UMAP is recomputed.
# ------------------------------------------------------------------

message("Transferring frozen REF12 umapRPCA coordinates...")
ref_umap <- Embeddings(ref12, reduction = "umapRPCA")

if (!all(map_for_old$ref12_cell_id %in% rownames(ref_umap))) {
  stop("Some mapped REF12 cells are absent from umapRPCA embeddings.")
}

umap_aligned <- ref_umap[map_for_old$ref12_cell_id, , drop = FALSE]
rownames(umap_aligned) <- colnames(canonical)

umap_key <- "UMAPREF12_"

canonical[["umapRPCA_ref12_v13"]] <- CreateDimReducObject(
  embeddings = umap_aligned,
  key = umap_key,
  assay = "RNA"
)

# Coordinate identity audit after transfer.
transferred_umap <- Embeddings(
  canonical,
  reduction = "umapRPCA_ref12_v13"
)

if (!isTRUE(all.equal(
  unname(transferred_umap),
  unname(umap_aligned),
  check.attributes = FALSE,
  tolerance = 0
))) {
  stop("Transferred UMAP coordinates are not exactly identical.")
}

# ------------------------------------------------------------------
# 6. Build compact QC tables. NO biological inference is performed.
# ------------------------------------------------------------------

sample_counts <- as.data.frame(
  table(
    canonical_sample = canonical$canonical_sample_v13,
    NAS_score = canonical$NAS_score_v13
  ),
  stringsAsFactors = FALSE
)
sample_counts <- sample_counts[sample_counts$Freq > 0, , drop = FALSE]
colnames(sample_counts)[3] <- "n_cells"

write.csv(
  sample_counts,
  file.path(TABDIR, "04_sample_cell_counts_v1.3.csv"),
  row.names = FALSE
)

celltype_counts <- as.data.frame(
  table(
    canonical_sample = canonical$canonical_sample_v13,
    celltype_ref12 = canonical$celltype_ref12_v13,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
celltype_counts <- celltype_counts[celltype_counts$Freq > 0, , drop = FALSE]
colnames(celltype_counts)[3] <- "n_cells"

write.csv(
  celltype_counts,
  file.path(TABDIR, "05_celltype_counts_by_sample_v1.3.csv"),
  row.names = FALSE
)

feature_manifest <- data.frame(
  assay = "RNA",
  n_features = nrow(canonical[["RNA"]]),
  inherited_from = OLD_RDS,
  feature_filter_applied_in_v1.3 = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  feature_manifest,
  file.path(TABDIR, "06_feature_manifest_v1.3.csv"),
  row.names = FALSE
)

metadata_manifest <- data.frame(
  metadata_group = c(
    "OLD metadata inherited with OLD LT subset",
    "REF12 metadata transferred with ref12_ prefix",
    "Canonical standardized v1.3 metadata"
  ),
  n_columns = c(
    ncol(old_md),
    if (TRANSFER_ALL_REF12_METADATA) ncol(ref_md) else 0L,
    ncol(std_meta)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_manifest,
  file.path(TABDIR, "07_metadata_transfer_manifest_v1.3.csv"),
  row.names = FALSE
)

reduction_manifest <- data.frame(
  reduction = Reductions(canonical),
  source = vapply(
    Reductions(canonical),
    function(r) {
      if (r == "umapRPCA_ref12_v13") {
        "Transferred exactly from REF12 umapRPCA; no recomputation"
      } else {
        "Inherited from OLD LT subset"
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  reduction_manifest,
  file.path(TABDIR, "08_reduction_manifest_v1.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 7. Final pre-save invariants.
# ------------------------------------------------------------------

stopifnot(ncol(canonical) == 81348L)
stopifnot(nrow(canonical[["RNA"]]) == canonical_feature_n)
stopifnot(all(canonical$canonical_sample_v13 %in% LT_SAMPLES))
stopifnot(all(canonical$NAS_score_v13 %in% c(0, 2, 4, 6)))
stopifnot(!anyNA(canonical$celltype_ref12_v13))
stopifnot("umapRPCA_ref12_v13" %in% Reductions(canonical))

# IMPORTANT:
# celltype_ref12_v13 is a transferred reference annotation, NOT a frozen final
# biological annotation. Broad-celltype validation is the next analysis stage.

message("Saving NEW canonical object...")
saveRDS(canonical, OUT_RDS, compress = FALSE)

if (!file.exists(OUT_RDS)) {
  stop("saveRDS completed but output file was not found.")
}

out_info <- file.info(OUT_RDS)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 canonical object build ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  "SOURCE OBJECTS (UNCHANGED):",
  paste0("Count/expression source: ", OLD_RDS),
  paste0("Metadata/annotation/UMAP reference: ", REF12_RDS),
  "",
  "CANONICAL OBJECT:",
  paste0("Output RDS: ", OUT_RDS),
  paste0("Cells: ", ncol(canonical)),
  paste0("RNA features: ", nrow(canonical[["RNA"]])),
  paste0("Assays: ", paste(Assays(canonical), collapse = " | ")),
  paste0("Reductions: ", paste(Reductions(canonical), collapse = " | ")),
  paste0("Metadata columns: ", ncol(canonical[[]])),
  paste0("Output size GB: ", round(out_info$size / 1024^3, 3)),
  "",
  "SAMPLE / NAS:",
  paste(
    capture.output(print(sample_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "PROVENANCE:",
  "- RNA/count universe inherited from OLD clustered2.0 LT1-LT8 subset.",
  "- No feature filtering was performed by v1.3.",
  "- REF12 metadata transferred by exact sample-specific 10x barcode mapping.",
  "- REF12 metadata preserved under ref12_ prefix.",
  "- Standard v1.3 fields added without replacing OLD metadata.",
  "- REF12 umapRPCA coordinates transferred exactly into umapRPCA_ref12_v13.",
  "- No UMAP was recalculated.",
  "",
  "IMPORTANT INTERPRETATION:",
  "- celltype_ref12_v13 is a transferred reference annotation.",
  "- It is NOT yet the biologically validated/frozen human broad-celltype annotation.",
  "- No biological differential analysis was performed in v1.3.",
  "",
  "NEXT DECISION GATE:",
  "Broad whole-liver celltype annotation validation using the 16,424-gene canonical object.",
  "",
  "SOURCE RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "09_build_summary_v1.3.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.3 COMPLETE")
message("Cells: ", ncol(canonical))
message("RNA features: ", nrow(canonical[["RNA"]]))
message("Output: ", OUT_RDS)
message("No source RDS was modified.")
message("============================================================")
