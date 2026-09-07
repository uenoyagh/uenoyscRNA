#!/usr/bin/env Rscript

# Human_MASH_LT1to8_build_clean_myeloid_object_v1.6.R
#
# PURPOSE
# -------
# Build a NEW clean human LT1-LT8 myeloid-entry object after the v1.5.1/v1.5.2
# immune identity audits.
#
# RETAINED res0.6 clusters:
#   1, 2, 3, 4, 9, 10, 18
#
# EXCLUDED FROM THE CLEAN CORE:
#   0,5,6,7,8,12,13,14,16,17,20 = non-myeloid / unresolved immune
#   15 = mesenchymal contamination
#   11 = macrophage + strict-hepatocyte mixed cluster
#   19 = macrophage + mesenchymal mixed cluster
#
# IMPORTANT:
#   cluster 11 and 19 are NOT deleted from source objects.
#   They are excluded only from this clean core object and remain recoverable
#   from the v1.5 immune object for later dedicated review.
#
# THIS SCRIPT DOES:
#   - subset the v1.5 immune object by the explicit retained cluster list
#   - build a FRESH Seurat object from RAW RNA counts
#   - preserve source metadata
#   - transfer the existing parent immune UMAP coordinates exactly for QC
#   - save a NEW RDS
#
# THIS SCRIPT DOES NOT:
#   - modify the v1.5 immune object
#   - normalize / SCTransform
#   - integrate
#   - recluster
#   - calculate a new UMAP
#   - assign final macrophage/monocyte biological states
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5/",
  "objects/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_MyeloidClean_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_MyeloidClean_", VERSION, ".rds")
)

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.5 immune RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

message("============================================================")
message("Human MASH LT1-LT8 clean myeloid object build ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("Output: ", OUT_RDS)

CLUSTER_COL <- "immune_RPCA_res0_6"

RETAIN_CLUSTERS <- c("1", "2", "3", "4", "9", "10", "18")

gate_manifest <- data.frame(
  cluster = as.character(0:20),
  action = c(
    "exclude",  # 0
    "retain",   # 1
    "retain",   # 2
    "retain",   # 3
    "retain",   # 4
    "exclude",  # 5
    "exclude",  # 6
    "exclude",  # 7
    "exclude",  # 8
    "retain",   # 9
    "retain",   # 10
    "exclude",  # 11
    "exclude",  # 12
    "exclude",  # 13
    "exclude",  # 14
    "exclude",  # 15
    "exclude",  # 16
    "exclude",  # 17
    "retain",   # 18
    "exclude",  # 19
    "exclude"   # 20
  ),
  reason = c(
    "T/NK mixed; not clean myeloid core",                         # 0
    "myeloid/macrophage-supported",                               # 1
    "resident-macrophage-supported; retain with contamination QC flag", # 2
    "resident-macrophage-supported",                              # 3
    "myeloid/monocyte-macrophage-supported",                      # 4
    "T-cell-supported",                                           # 5
    "NK-cell-supported",                                          # 6
    "NK-cell-supported",                                          # 7
    "T-cell-supported",                                           # 8
    "macrophage-supported; retain with contamination QC flag",    # 9
    "myeloid/monocyte-macrophage-supported",                      # 10
    "exclude from clean core: macrophage + strict-hepatocyte mixed; elevated RNA burden", # 11
    "T/NK mixed; not clean myeloid core",                         # 12
    "plasma-cell-supported",                                      # 13
    "B-cell-supported",                                           # 14
    "mesenchymal contamination",                                  # 15
    "NK/mixed immune; not clean myeloid core",                    # 16
    "low-confidence/non-core immune; exclude from clean core",    # 17
    "inflammatory-monocyte-supported",                            # 18
    "exclude from clean core: macrophage + mesenchymal mixed",    # 19
    "tiny unresolved cluster; only 5 samples; exclude from clean core" # 20
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gate_manifest,
  file.path(TABDIR, "01_cluster_gate_manifest_v1.6.csv"),
  row.names = FALSE
)

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  CLUSTER_COL,
  "sample_v15",
  "NAS_v15",
  "source_broad_label_v15"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))
if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")
if (!"umap.immune.rpca" %in% Reductions(obj)) {
  stop("Parent immune UMAP reduction 'umap.immune.rpca' absent.")
}

parent_cluster <- as.character(obj[[CLUSTER_COL, drop = TRUE]])

# Audit actual cluster universe before subsetting.
observed_clusters <- sort(unique(parent_cluster))
expected_clusters <- as.character(0:20)

if (!setequal(observed_clusters, expected_clusters)) {
  warning(
    "Observed res0.6 cluster universe differs from 0:20.\nObserved: ",
    paste(observed_clusters, collapse = " | ")
  )
}

retain_cells <- colnames(obj)[parent_cluster %in% RETAIN_CLUSTERS]

if (!length(retain_cells)) stop("No retained cells found.")

message("Retained clean-myeloid cells: ", length(retain_cells))

# Expected from v1.5 summary:
# 1=1562, 2=1554, 3=1247, 4=1223, 9=693, 10=579, 18=98 => total 6956
EXPECTED_TOTAL <- 6956L

if (length(retain_cells) != EXPECTED_TOTAL) {
  stop(
    "Retained cell total changed from audited expectation. Observed=",
    length(retain_cells), ", expected=", EXPECTED_TOTAL
  )
}

# ------------------------------------------------------------------
# Fresh object from raw counts only.
# ------------------------------------------------------------------

parent_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
clean_counts <- parent_counts[, retain_cells, drop = FALSE]
clean_meta <- obj[[]][retain_cells, , drop = FALSE]

clean <- CreateSeuratObject(
  counts = clean_counts,
  assay = "RNA",
  meta.data = clean_meta,
  project = "Human_MASH_LT1to8_MyeloidClean_v1.6",
  min.cells = 0,
  min.features = 0
)

clean$myeloid_gate_v16 <- "retained_clean_core"
clean$parent_immune_cluster_res0_6_v16 <- as.character(
  clean[[CLUSTER_COL, drop = TRUE]]
)
clean$myeloid_clean_version_v16 <- VERSION

# Preserve an intentionally broad entry-class for navigation only.
# This is NOT a final biological annotation.
clean$entry_class_v16 <- ifelse(
  clean$parent_immune_cluster_res0_6_v16 == "18",
  "Monocyte_candidate",
  "Macrophage_or_myeloid_candidate"
)

# Flag the macrophage clusters that had measurable hepatocyte co-detection but
# were retained because they remained strongly macrophage-supported.
clean$contamination_review_flag_v16 <- ifelse(
  clean$parent_immune_cluster_res0_6_v16 == "2",
  "retain_flag_hepatocyte_signal",
  ifelse(
    clean$parent_immune_cluster_res0_6_v16 == "9",
    "retain_flag_low_hepatocyte_signal",
    "none"
  )
)

# ------------------------------------------------------------------
# Transfer parent immune UMAP coordinates exactly; do not recalculate.
# ------------------------------------------------------------------

parent_umap <- Embeddings(obj, reduction = "umap.immune.rpca")
if (!all(colnames(clean) %in% rownames(parent_umap))) {
  stop("Some retained cells are absent from parent immune UMAP.")
}

umap_clean <- parent_umap[colnames(clean), , drop = FALSE]

clean[["umap_parent_immune_v16"]] <- CreateDimReducObject(
  embeddings = umap_clean,
  key = "ParentImmuneUMAP_",
  assay = "RNA"
)

transferred <- Embeddings(clean, reduction = "umap_parent_immune_v16")

if (!isTRUE(all.equal(
  unname(transferred),
  unname(umap_clean),
  check.attributes = FALSE,
  tolerance = 0
))) {
  stop("Transferred parent UMAP coordinates are not exactly identical.")
}

# ------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(
    parent_cluster = clean$parent_immune_cluster_res0_6_v16
  ),
  stringsAsFactors = FALSE
)
colnames(cluster_counts)[2] <- "n_cells"

write.csv(
  cluster_counts,
  file.path(TABDIR, "02_retained_cluster_counts_v1.6.csv"),
  row.names = FALSE
)

sample_counts <- as.data.frame(
  table(
    sample = clean$sample_v15,
    NAS = clean$NAS_v15
  ),
  stringsAsFactors = FALSE
)
sample_counts <- sample_counts[sample_counts$Freq > 0, , drop = FALSE]
colnames(sample_counts)[3] <- "n_cells"

write.csv(
  sample_counts,
  file.path(TABDIR, "03_retained_sample_counts_v1.6.csv"),
  row.names = FALSE
)

cluster_sample <- as.data.frame(
  table(
    parent_cluster = clean$parent_immune_cluster_res0_6_v16,
    sample = clean$sample_v15
  ),
  stringsAsFactors = FALSE
)
cluster_sample <- cluster_sample[cluster_sample$Freq > 0, , drop = FALSE]
colnames(cluster_sample)[3] <- "n_cells"

write.csv(
  cluster_sample,
  file.path(TABDIR, "04_retained_cluster_by_sample_v1.6.csv"),
  row.names = FALSE
)

feature_manifest <- data.frame(
  assay = "RNA",
  n_features = nrow(clean[["RNA"]]),
  source = INPUT_RDS,
  layer = paste(Layers(clean[["RNA"]]), collapse = " | "),
  additional_feature_filter_applied_v16 = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  feature_manifest,
  file.path(TABDIR, "05_feature_manifest_v1.6.csv"),
  row.names = FALSE
)

# Keep an explicit list of mixed clusters that remain recoverable from parent.
mixed_preservation <- data.frame(
  excluded_cluster = c("11", "19"),
  n_cells_parent = c(
    sum(parent_cluster == "11"),
    sum(parent_cluster == "19")
  ),
  status = c(
    "preserved_in_parent_v1.5_not_in_clean_core",
    "preserved_in_parent_v1.5_not_in_clean_core"
  ),
  reason = c(
    "macrophage + strict-hepatocyte co-detection; high nCount/nFeature",
    "macrophage + mesenchymal co-detection"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mixed_preservation,
  file.path(TABDIR, "06_mixed_cluster_preservation_manifest_v1.6.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Final invariants and NEW RDS save
# ------------------------------------------------------------------

stopifnot(ncol(clean) == EXPECTED_TOTAL)
stopifnot(nrow(clean[["RNA"]]) == 16424L)
stopifnot(
  setequal(
    unique(clean$parent_immune_cluster_res0_6_v16),
    RETAIN_CLUSTERS
  )
)
stopifnot(all(clean$sample_v15 %in% paste0("LT", 1:8)))
stopifnot(!anyNA(clean$sample_v15))
stopifnot(!anyNA(clean$NAS_v15))

message("Saving NEW clean myeloid object...")
saveRDS(clean, OUT_RDS, compress = FALSE)

if (!file.exists(OUT_RDS)) stop("Output RDS was not created.")

info <- file.info(OUT_RDS)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 clean myeloid object build ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Source immune RDS: ", INPUT_RDS),
  paste0("Output clean RDS: ", OUT_RDS),
  "",
  paste0("Cells: ", ncol(clean)),
  paste0("RNA features: ", nrow(clean[["RNA"]])),
  paste0("RNA layers: ", paste(Layers(clean[["RNA"]]), collapse = " | ")),
  paste0("Transferred reductions: ", paste(Reductions(clean), collapse = " | ")),
  paste0("Output size GB: ", round(info$size / 1024^3, 3)),
  "",
  "RETAINED PARENT res0.6 CLUSTERS:",
  paste(RETAIN_CLUSTERS, collapse = " | "),
  "",
  "RETAINED CLUSTER COUNTS:",
  paste(
    capture.output(print(cluster_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "MIXED CLUSTERS EXCLUDED FROM CLEAN CORE BUT PRESERVED IN PARENT:",
  "cluster 11: macrophage + strict-hepatocyte mixed",
  "cluster 19: macrophage + mesenchymal mixed",
  "",
  "IMPORTANT:",
  "- This object is an entry-gated clean myeloid core, not a final MΦ/Monocyte annotation.",
  "- No normalization, RPCA integration, reclustering, or new UMAP was performed.",
  "- Parent immune UMAP coordinates were transferred exactly for QC only.",
  "- Cluster 2 and 9 are retained but carry contamination-review flags.",
  "",
  "NEXT DECISION GATE:",
  "Perform sample-aware RPCA reclustering of this 6,956-cell clean myeloid object",
  "and then define Human MΦ / Monocyte states.",
  "",
  "SOURCE RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "07_build_summary_v1.6.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.6 COMPLETE")
message("Clean myeloid cells: ", ncol(clean))
message("RNA features: ", nrow(clean[["RNA"]]))
message("Output RDS: ", OUT_RDS)
message("No source RDS modified.")
message("============================================================")
