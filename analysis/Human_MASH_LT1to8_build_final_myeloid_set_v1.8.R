#!/usr/bin/env Rscript

# Human_MASH_LT1to8_build_final_myeloid_set_v1.8.R
#
# PURPOSE
# -------
# Build a NEW final clean Human MΦ/Monocyte ENTRY SET from the v1.7 RPCA object,
# using res0.8 as the working resolution and the v1.7.2/v1.7.3 audits.
#
# RETAIN:
#   0,1,2,3,4,5,6,7,8,9,10,11,13,14,16
#
# EXCLUDE:
#   12 = DC-like
#   15 = T/NK-mixed myeloid
#
# IMPORTANT
# ---------
# This script freezes the CLEAN CELL SET and the res0.8 numeric cluster
# membership, but DOES NOT freeze final biological macrophage-state names.
#
# It adds only PRELIMINARY broad candidate classes:
#   Macrophage_candidate
#   Monocyte_candidate
#   Transitional_myeloid_candidate
#   SPP1_TREM2_LAM_candidate
#
# Cluster-specific review flags are retained for small / sample-limited groups.
#
# THIS SCRIPT DOES:
#   - subset explicitly by retained res0.8 clusters
#   - build a fresh Seurat object from RAW RNA counts
#   - preserve metadata
#   - transfer the existing v1.7 myeloid UMAP coordinates exactly for QC
#   - save a NEW RDS
#
# THIS SCRIPT DOES NOT:
#   - modify/overwrite the v1.7 RDS
#   - normalize / integrate / recluster
#   - calculate a new UMAP
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#   - declare preliminary candidate classes to be final biological states
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.8"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_RPCA_v1.7/",
  "objects/Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidSet_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_FinalMyeloidSet_", VERSION, ".rds")
)

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.7 RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

message("============================================================")
message("Human MASH LT1-LT8 final clean myeloid set build ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("Output: ", OUT_RDS)

CLUSTER_COL <- "myeloid_RPCA_res0_8"

RETAIN_CLUSTERS <- c(
  "0","1","2","3","4","5","6","7","8","9","10","11","13","14","16"
)

EXCLUDE_CLUSTERS <- c("12", "15")

EXPECTED_TOTAL <- 6720L

# Preliminary broad candidate classes only.
candidate_class_map <- c(
  "0"  = "Transitional_myeloid_candidate",
  "1"  = "Macrophage_candidate",
  "2"  = "Macrophage_candidate",
  "3"  = "Macrophage_candidate",
  "4"  = "Macrophage_candidate",
  "5"  = "Macrophage_candidate",
  "6"  = "SPP1_TREM2_LAM_candidate",
  "7"  = "Monocyte_candidate",
  "8"  = "Macrophage_candidate",
  "9"  = "Macrophage_candidate",
  "10" = "Transitional_myeloid_candidate",
  "11" = "Macrophage_candidate",
  "13" = "Monocyte_candidate",
  "14" = "Transitional_myeloid_candidate",
  "16" = "Macrophage_candidate"
)

review_flag_map <- c(
  "0"  = "none",
  "1"  = "none",
  "2"  = "none",
  "3"  = "none",
  "4"  = "none",
  "5"  = "none",
  "6"  = "none",
  "7"  = "none",
  "8"  = "none",
  "9"  = "none",
  "10" = "none",
  "11" = "none",
  "13" = "retain_low_RNA_monocyte_candidate",
  "14" = "retain_small_6sample_lipid_inflammatory_candidate",
  "16" = "retain_sample_limited_2sample_not_independent_state"
)

gate_manifest <- data.frame(
  cluster_res0_8 = as.character(0:16),
  action = ifelse(
    as.character(0:16) %in% RETAIN_CLUSTERS,
    "retain",
    "exclude"
  ),
  rationale = c(
    "retain: monocyte/macrophage transitional evidence",           # 0
    "retain: resident macrophage-supported",                       # 1
    "retain: resident macrophage-supported",                       # 2
    "retain: resident macrophage-supported",                       # 3
    "retain: resident macrophage-supported",                       # 4
    "retain: resident/repair macrophage-supported",                # 5
    "retain: SPP1/TREM2-LAM-supported",                            # 6
    "retain: inflammatory-monocyte-supported",                     # 7
    "retain: resident/repair macrophage-supported",                # 8
    "retain: resident/repair macrophage-supported",                # 9
    "retain: inflammatory/transitional myeloid",                   # 10
    "retain: resident/lipid-handling macrophage",                  # 11
    "exclude: coherent DC-like program exceeds myeloid identity", # 12
    "retain: low-RNA monocyte candidate; non-myeloid signal low",  # 13
    "retain: small lipid/inflammatory myeloid candidate",          # 14
    "exclude: coherent T/NK mixed signal within myeloid cluster",  # 15
    "retain: coherent resident/repair macrophage but only 2 samples; not independent state" # 16
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gate_manifest,
  file.path(TABDIR, "01_res08_final_gate_manifest_v1.8.csv"),
  row.names = FALSE
)

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  CLUSTER_COL,
  "sample_v17",
  "NAS_v17",
  "parent_cluster_v17"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")
if (!"umap.myeloid.rpca" %in% Reductions(obj)) {
  stop("umap.myeloid.rpca reduction absent.")
}

cl <- as.character(obj[[CLUSTER_COL, drop = TRUE]])

retain_cells <- colnames(obj)[cl %in% RETAIN_CLUSTERS]

if (length(retain_cells) != EXPECTED_TOTAL) {
  stop(
    "Retained total differs from audited expectation. Observed=",
    length(retain_cells),
    ", expected=", EXPECTED_TOTAL
  )
}

message("Retained final clean myeloid cells: ", length(retain_cells))

# ------------------------------------------------------------------
# Fresh object from raw counts only.
# ------------------------------------------------------------------

raw_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
keep_counts <- raw_counts[, retain_cells, drop = FALSE]
keep_meta <- obj[[]][retain_cells, , drop = FALSE]

clean <- CreateSeuratObject(
  counts = keep_counts,
  assay = "RNA",
  meta.data = keep_meta,
  project = "Human_MASH_LT1to8_FinalMyeloidSet_v1.8",
  min.cells = 0,
  min.features = 0
)

clean$res0_8_cluster_frozen_v18 <- as.character(
  clean[[CLUSTER_COL, drop = TRUE]]
)

clean$myeloid_candidate_class_v18 <- unname(
  candidate_class_map[clean$res0_8_cluster_frozen_v18]
)

clean$review_flag_v18 <- unname(
  review_flag_map[clean$res0_8_cluster_frozen_v18]
)

clean$final_myeloid_set_version_v18 <- VERSION
clean$final_myeloid_set_status_v18 <- "clean_cell_set_frozen_state_labels_not_final"

if (anyNA(clean$myeloid_candidate_class_v18)) {
  stop("Candidate-class mapping produced NA.")
}

# ------------------------------------------------------------------
# Transfer parent v1.7 UMAP coordinates exactly for QC.
# ------------------------------------------------------------------

parent_umap <- Embeddings(obj, reduction = "umap.myeloid.rpca")

if (!all(colnames(clean) %in% rownames(parent_umap))) {
  stop("Some retained cells are absent from parent v1.7 UMAP.")
}

umap_keep <- parent_umap[colnames(clean), , drop = FALSE]

clean[["umap_myeloid_parent_v18"]] <- CreateDimReducObject(
  embeddings = umap_keep,
  key = "MyeloidParentUMAP_",
  assay = "RNA"
)

transferred <- Embeddings(
  clean,
  reduction = "umap_myeloid_parent_v18"
)

if (!isTRUE(all.equal(
  unname(transferred),
  unname(umap_keep),
  check.attributes = FALSE,
  tolerance = 0
))) {
  stop("Transferred UMAP coordinates are not exactly identical.")
}

# ------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(cluster_res0_8 = clean$res0_8_cluster_frozen_v18),
  stringsAsFactors = FALSE
)
colnames(cluster_counts)[2] <- "n_cells"

write.csv(
  cluster_counts,
  file.path(TABDIR, "02_retained_cluster_counts_v1.8.csv"),
  row.names = FALSE
)

class_counts <- as.data.frame(
  table(candidate_class = clean$myeloid_candidate_class_v18),
  stringsAsFactors = FALSE
)
colnames(class_counts)[2] <- "n_cells"

write.csv(
  class_counts,
  file.path(TABDIR, "03_candidate_class_counts_v1.8.csv"),
  row.names = FALSE
)

sample_counts <- as.data.frame(
  table(sample = clean$sample_v17, NAS = clean$NAS_v17),
  stringsAsFactors = FALSE
)
sample_counts <- sample_counts[sample_counts$Freq > 0, , drop = FALSE]
colnames(sample_counts)[3] <- "n_cells"

write.csv(
  sample_counts,
  file.path(TABDIR, "04_sample_counts_v1.8.csv"),
  row.names = FALSE
)

cluster_sample <- as.data.frame(
  table(
    cluster_res0_8 = clean$res0_8_cluster_frozen_v18,
    sample = clean$sample_v17
  ),
  stringsAsFactors = FALSE
)
cluster_sample <- cluster_sample[cluster_sample$Freq > 0, , drop = FALSE]
colnames(cluster_sample)[3] <- "n_cells"

write.csv(
  cluster_sample,
  file.path(TABDIR, "05_cluster_by_sample_v1.8.csv"),
  row.names = FALSE
)

excluded_manifest <- data.frame(
  excluded_cluster_res0_8 = EXCLUDE_CLUSTERS,
  n_cells_parent_v17 = c(
    sum(cl == "12"),
    sum(cl == "15")
  ),
  reason = c(
    "DC-like: Dendritic >=2 markers 78.3%; exclude from final MΦ/Monocyte set",
    "T/NK-mixed: T >=2 25.3%, NK >=2 32.0%; exclude from final MΦ/Monocyte set"
  ),
  preserved_in_source_v17 = TRUE,
  stringsAsFactors = FALSE
)

write.csv(
  excluded_manifest,
  file.path(TABDIR, "06_excluded_cluster_manifest_v1.8.csv"),
  row.names = FALSE
)

feature_manifest <- data.frame(
  assay = "RNA",
  n_features = nrow(clean[["RNA"]]),
  layers = paste(Layers(clean[["RNA"]]), collapse = " | "),
  source = INPUT_RDS,
  additional_feature_filter_applied_v18 = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  feature_manifest,
  file.path(TABDIR, "07_feature_manifest_v1.8.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Final invariants and save NEW object
# ------------------------------------------------------------------

stopifnot(ncol(clean) == EXPECTED_TOTAL)
stopifnot(nrow(clean[["RNA"]]) == 16424L)
stopifnot(
  setequal(
    unique(clean$res0_8_cluster_frozen_v18),
    RETAIN_CLUSTERS
  )
)
stopifnot(!any(clean$res0_8_cluster_frozen_v18 %in% EXCLUDE_CLUSTERS))
stopifnot(!anyNA(clean$sample_v17))
stopifnot(!anyNA(clean$NAS_v17))

message("Saving NEW final clean myeloid-set object...")
saveRDS(clean, OUT_RDS, compress = FALSE)

if (!file.exists(OUT_RDS)) stop("Output RDS was not created.")

info <- file.info(OUT_RDS)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final clean myeloid set build ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Source v1.7 RDS: ", INPUT_RDS),
  paste0("Output v1.8 RDS: ", OUT_RDS),
  "",
  paste0("Cells: ", ncol(clean)),
  paste0("RNA features: ", nrow(clean[["RNA"]])),
  paste0("RNA layers: ", paste(Layers(clean[["RNA"]]), collapse = " | ")),
  paste0("Transferred reduction: umap_myeloid_parent_v18"),
  paste0("Output size GB: ", round(info$size / 1024^3, 3)),
  "",
  "WORKING/FROZEN NUMERIC RESOLUTION:",
  "res0.8 numeric cluster membership is frozen for the clean-cell set.",
  "",
  "RETAINED res0.8 CLUSTERS:",
  paste(RETAIN_CLUSTERS, collapse = " | "),
  "",
  "EXCLUDED res0.8 CLUSTERS:",
  "12 = DC-like",
  "15 = T/NK-mixed myeloid",
  "",
  "PRELIMINARY CANDIDATE CLASS COUNTS:",
  paste(
    capture.output(print(class_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "IMPORTANT:",
  "- The CLEAN CELL SET and numeric res0.8 cluster membership are frozen.",
  "- myeloid_candidate_class_v18 is PRELIMINARY and is NOT a final biological state annotation.",
  "- cluster 13 is retained as low-RNA monocyte candidate.",
  "- cluster 14 is retained as small lipid/inflammatory myeloid candidate.",
  "- cluster 16 is retained because its myeloid/resident signal is coherent, but it is sample-limited",
  "  and should not be frozen as an independent biological state.",
  "- No normalization, reclustering, new UMAP, DE, pseudobulk, or pathway analysis was performed.",
  "",
  "NEXT DECISION GATE:",
  "Perform detailed Human MΦ/Monocyte state annotation on this 6,720-cell clean set,",
  "using cluster marker programs and then disease-stage behavior as validation.",
  "",
  "SOURCE RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "08_build_summary_v1.8.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.8 COMPLETE")
message("Final clean myeloid cells: ", ncol(clean))
message("RNA features: ", nrow(clean[["RNA"]]))
message("Output RDS: ", OUT_RDS)
message("Final biological state labels NOT yet frozen.")
message("============================================================")
