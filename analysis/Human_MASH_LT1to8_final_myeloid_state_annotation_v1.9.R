#!/usr/bin/env Rscript

# Human_MASH_LT1to8_final_myeloid_state_annotation_v1.9.R
#
# PURPOSE
# -------
# Freeze FINAL broad Human MΦ/Monocyte state labels on the v1.8 clean-cell set
# using the audited res0.8 numeric clusters and multi-marker programs.
#
# IMPORTANT PRINCIPLE
# -------------------
# Biological labels are assigned from the established cluster-level programs
# defined in v1.7.2/v1.7.3. Disease-stage composition (NAS0/2/4/6) is NOT used
# to assign labels; it is carried only for later validation.
#
# THIS SCRIPT DOES:
#   - read v1.8 final clean myeloid set
#   - assign final broad state labels by explicit res0.8 cluster map
#   - preserve numeric cluster IDs and provenance
#   - generate frozen annotation tables
#   - generate frozen UMAP with final labels using transferred coordinates
#   - generate final marker DotPlot
#   - save a NEW annotated RDS
#
# THIS SCRIPT DOES NOT:
#   - recluster
#   - recompute UMAP
#   - modify/overwrite v1.8
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#
# State framework:
#   Resident/Kupffer-like MΦ
#   Repair/Resolution MΦ
#   Lipid-associated/TREM2 MΦ
#   Inflammatory/Transitional MΦ
#   Inflammatory Monocyte
#
# Small/sample-limited clusters are assigned to the nearest supported broad state,
# not frozen as independent states.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.9"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidSet_v1.8/",
  "objects/Human_MASH_LT1to8_FinalMyeloidSet_v1.8.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidAnnotation_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_FinalMyeloidAnnotation_", VERSION, ".rds")
)

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.8 RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 final myeloid annotation ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("No reclustering; no UMAP recomputation.")

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  "res0_8_cluster_frozen_v18",
  "sample_v17",
  "NAS_v17"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")
if (!"umap_myeloid_parent_v18" %in% Reductions(obj)) {
  stop("Transferred parent myeloid UMAP absent.")
}

cl <- as.character(obj$res0_8_cluster_frozen_v18)

# ------------------------------------------------------------------
# FINAL broad state map.
#
# Rationale:
#  - 6: strongest SPP1/TREM2/LAM program -> Lipid-associated/TREM2 MΦ
#  - 7,13: monocyte-inflammatory dominant -> Inflammatory Monocyte
#  - 0,10,14: mixed monocyte/macrophage inflammatory-transition programs
#              -> Inflammatory/Transitional MΦ
#  - 1,2,3,4,5,8,9,11,16: macrophage-dominant.
#    Within these, clusters with especially strong repair/resolution/lipid programs
#    are grouped as Repair/Resolution MΦ; the remaining macrophage-dominant
#    clusters are Resident/Kupffer-like MΦ.
#
# Small cluster 16 is absorbed into Repair/Resolution MΦ rather than treated as
# an independent state because it is sample-limited.
# ------------------------------------------------------------------

state_map <- c(
  "0"  = "Inflammatory/Transitional_Mphi",
  "1"  = "Resident/Kupffer-like_Mphi",
  "2"  = "Resident/Kupffer-like_Mphi",
  "3"  = "Resident/Kupffer-like_Mphi",
  "4"  = "Repair/Resolution_Mphi",
  "5"  = "Repair/Resolution_Mphi",
  "6"  = "Lipid-associated/TREM2_Mphi",
  "7"  = "Inflammatory_Monocyte",
  "8"  = "Repair/Resolution_Mphi",
  "9"  = "Repair/Resolution_Mphi",
  "10" = "Inflammatory/Transitional_Mphi",
  "11" = "Resident/Kupffer-like_Mphi",
  "13" = "Inflammatory_Monocyte",
  "14" = "Inflammatory/Transitional_Mphi",
  "16" = "Repair/Resolution_Mphi"
)

rationale_map <- c(
  "0"  = "monocyte-inflammatory plus resident-macrophage transitional program",
  "1"  = "strong resident/Kupffer-like macrophage program",
  "2"  = "resident/Kupffer-like macrophage program",
  "3"  = "resident/Kupffer-like macrophage program",
  "4"  = "resident macrophage with strong repair/resolution program",
  "5"  = "resident macrophage with very strong repair/resolution program",
  "6"  = "strongest SPP1/TREM2/LAM program",
  "7"  = "strong inflammatory-monocyte program",
  "8"  = "resident macrophage with strong repair/resolution program",
  "9"  = "resident macrophage with strong repair/resolution program",
  "10" = "inflammatory/transitional myeloid program",
  "11" = "resident/lipid-handling macrophage program",
  "13" = "low-RNA inflammatory monocyte candidate; retained as monocyte",
  "14" = "small lipid-handling plus inflammatory/transitional myeloid program",
  "16" = "sample-limited resident/repair macrophage; absorbed into repair state"
)

audit_note_map <- c(
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
  "13" = "low_RNA_monocyte",
  "14" = "small_6sample_cluster",
  "16" = "sample_limited_2sample_absorbed_into_broad_state"
)

if (!all(unique(cl) %in% names(state_map))) {
  bad <- setdiff(unique(cl), names(state_map))
  stop("Unmapped retained cluster(s): ", paste(bad, collapse = ", "))
}

obj$myeloid_state_FINAL_v19 <- unname(state_map[cl])
obj$myeloid_state_FINAL_rationale_v19 <- unname(rationale_map[cl])
obj$myeloid_state_FINAL_audit_note_v19 <- unname(audit_note_map[cl])
obj$myeloid_state_annotation_version_v19 <- VERSION
obj$myeloid_state_annotation_status_v19 <- "FINAL_broad_state_annotation"

# ------------------------------------------------------------------
# Annotation manifest
# ------------------------------------------------------------------

manifest <- data.frame(
  cluster_res0_8 = names(state_map),
  final_state_v19 = unname(state_map),
  rationale = unname(rationale_map),
  audit_note = unname(audit_note_map),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(TABDIR, "01_final_cluster_to_state_manifest_v1.9.csv"),
  row.names = FALSE
)

state_counts <- as.data.frame(
  table(final_state_v19 = obj$myeloid_state_FINAL_v19),
  stringsAsFactors = FALSE
)
colnames(state_counts)[2] <- "n_cells"

write.csv(
  state_counts,
  file.path(TABDIR, "02_final_state_counts_v1.9.csv"),
  row.names = FALSE
)

state_by_sample <- as.data.frame(
  table(
    final_state_v19 = obj$myeloid_state_FINAL_v19,
    sample = obj$sample_v17
  ),
  stringsAsFactors = FALSE
)
state_by_sample <- state_by_sample[state_by_sample$Freq > 0, , drop = FALSE]
colnames(state_by_sample)[3] <- "n_cells"

write.csv(
  state_by_sample,
  file.path(TABDIR, "03_final_state_by_sample_v1.9.csv"),
  row.names = FALSE
)

state_by_nas <- as.data.frame(
  table(
    final_state_v19 = obj$myeloid_state_FINAL_v19,
    NAS = obj$NAS_v17
  ),
  stringsAsFactors = FALSE
)
state_by_nas <- state_by_nas[state_by_nas$Freq > 0, , drop = FALSE]
colnames(state_by_nas)[3] <- "n_cells"

tot <- aggregate(n_cells ~ final_state_v19, state_by_nas, sum)
colnames(tot)[2] <- "state_total"
state_by_nas <- merge(state_by_nas, tot, by = "final_state_v19")
state_by_nas$fraction_within_state <- state_by_nas$n_cells / state_by_nas$state_total

write.csv(
  state_by_nas,
  file.path(TABDIR, "04_final_state_by_NAS_descriptive_v1.9.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Marker panel for final annotation QC
# ------------------------------------------------------------------

marker_sets <- list(
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Resident_Kupffer_like = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2",
    "TIMD4", "MARCO", "VSIG4", "CD5L"
  ),
  Monocyte_inflammatory = c(
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA",
    "FCGR3A", "IL1B", "VCAN"
  ),
  SPP1_TREM2_LAM = c(
    "SPP1", "TREM2", "GPNMB", "LGALS3", "LPL", "FABP5", "APOE"
  ),
  ECM_inflammatory = c(
    "F13A1", "TGFB1", "VCAN", "THBS1", "FN1", "POSTN"
  ),
  Repair_resolution = c(
    "MERTK", "ANXA1", "LGALS3", "MRC1", "IL10RA", "IL10RB"
  ),
  Lipid_handling = c(
    "CD36", "MSR1", "OLR1", "APOE", "LPL", "FABP5"
  ),
  IFN_response = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44", "CXCL10"
  )
)

features <- rownames(obj[["RNA"]])
marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

# The v1.8 object has raw counts only; do not normalize it here.
# Reuse the final-state assignment for count-based detection QC.
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

states <- sort(unique(obj$myeloid_state_FINAL_v19))

marker_summary_list <- list()

for (st in states) {
  idx <- obj$myeloid_state_FINAL_v19 == st
  pct <- Matrix::rowMeans(counts[marker_order, idx, drop = FALSE] > 0) * 100

  marker_summary_list[[st]] <- data.frame(
    final_state_v19 = st,
    gene = marker_order,
    pct_detected = as.numeric(pct[marker_order]),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, marker_summary_list)
rownames(marker_summary) <- NULL

write.csv(
  marker_summary,
  file.path(TABDIR, "05_final_state_marker_detection_v1.9.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Final UMAP
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umap_myeloid_parent_v18")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  final_state = obj$myeloid_state_FINAL_v19,
  sample = obj$sample_v17,
  stringsAsFactors = FALSE
)

state_levels <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

plot_df$final_state <- factor(
  plot_df$final_state,
  levels = state_levels
)

centers <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ final_state,
  data = plot_df,
  FUN = median
)

# High-saturation palette consistent with project preference.
palette_values <- c(
  "Resident/Kupffer-like_Mphi"      = "#00A6D6",
  "Repair/Resolution_Mphi"          = "#00C853",
  "Lipid-associated/TREM2_Mphi"     = "#F0007A",
  "Inflammatory/Transitional_Mphi"  = "#FF7A00",
  "Inflammatory_Monocyte"           = "#7B1FA2"
)

p_umap <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = final_state)
) +
  geom_point(size = 0.48, alpha = 0.88) +
  geom_text(
    data = centers,
    aes(label = final_state),
    color = "black",
    fontface = "bold",
    size = 3.4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = palette_values, drop = FALSE) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 | FINAL broad myeloid states v1.9",
    "Frozen res0.8 cell set; transferred UMAP coordinates"
  )

ggsave(
  file.path(FIGDIR, "01_FINAL_myeloid_state_UMAP_v1.9.pdf"),
  p_umap,
  width = 10.5,
  height = 7.8
)

# Detection dot plot (size only because v1.8 intentionally contains raw counts only).
dot_df <- marker_summary
dot_df$final_state_v19 <- factor(
  dot_df$final_state_v19,
  levels = state_levels
)
dot_df$gene <- factor(
  dot_df$gene,
  levels = rev(marker_order)
)

p_dot <- ggplot(
  dot_df,
  aes(
    x = final_state_v19,
    y = gene,
    size = pct_detected
  )
) +
  geom_point() +
  scale_size(
    range = c(0.2, 7),
    limits = c(0, 100),
    name = "% detected"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 7.3)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 | FINAL broad myeloid-state marker detection",
    "Raw-count detection only; no new normalization in v1.9"
  )

ggsave(
  file.path(FIGDIR, "02_FINAL_myeloid_state_marker_detection_v1.9.pdf"),
  p_dot,
  width = 11,
  height = max(12, 0.22 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Save NEW annotated object
# ------------------------------------------------------------------

message("Saving NEW final annotated myeloid object...")
saveRDS(obj, OUT_RDS, compress = FALSE)

if (!file.exists(OUT_RDS)) stop("Output RDS was not created.")

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final broad myeloid annotation ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Input v1.8 RDS: ", INPUT_RDS),
  paste0("Output v1.9 RDS: ", OUT_RDS),
  "",
  paste0("Cells: ", ncol(obj)),
  paste0("RNA features: ", nrow(obj[["RNA"]])),
  paste0("Frozen numeric resolution: res0.8"),
  "",
  "FINAL BROAD STATES:",
  paste(
    capture.output(print(state_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "IMPORTANT:",
  "- Final broad state labels are frozen in myeloid_state_FINAL_v19.",
  "- Disease-stage/NAS composition was NOT used to assign labels.",
  "- NAS behavior remains descriptive and will be used only as downstream validation.",
  "- Small/sample-limited cluster 16 was absorbed into Repair/Resolution_Mphi.",
  "- No reclustering, new UMAP, DE, pseudobulk, or pathway analysis was performed.",
  "",
  "NEXT STEP:",
  "Quantify NAS0 -> NAS2 -> NAS4 -> NAS6 changes in final myeloid-state fractions,",
  "then perform sample-level pseudobulk/state-specific expression analyses.",
  "",
  "SOURCE v1.8 RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "06_build_summary_v1.9.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.9 COMPLETE")
message("Final broad state annotation frozen.")
message("Output RDS: ", OUT_RDS)
message("No source RDS modified.")
message("============================================================")
