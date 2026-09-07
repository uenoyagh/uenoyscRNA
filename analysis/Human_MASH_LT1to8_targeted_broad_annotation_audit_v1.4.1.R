#!/usr/bin/env Rscript

# Human_MASH_LT1to8_targeted_broad_annotation_audit_v1.4.1.R
#
# PURPOSE
# -------
# Resolve the remaining uncertainties identified by v1.4 before freezing
# broad whole-liver annotations in the canonical LT1-LT8 object.
#
# Focus:
#   1) Is "Inflammatory_monocyte_macrophage" predominantly myeloid, or does it
#      contain appreciable T/NK/B/plasma-cell populations?
#   2) Does "SPP1_TREM2_macrophage" show coherent macrophage/SPP1-TREM2 evidence,
#      or mixed mesenchymal/hepatocyte evidence?
#   3) What lineage underlies "Cycling_cell"?
#   4) Is "Hepatoblast_or_immature_hepatocyte" clearly hepatocyte-lineage, and
#      is "immature/hepatoblast" specifically supported?
#   5) Are the tiny Cholangiocyte/Fibroblast/Mesothelial labels credible?
#
# This is an AUDIT ONLY.
# It DOES NOT:
#   - relabel cells
#   - modify or save an RDS
#   - normalize / SCTransform
#   - integrate
#   - cluster / recluster
#   - calculate UMAP
#   - perform DE / pseudobulk / pathway / cell-cell analysis
#
# Seurat v5:
#   GetAssayData(..., layer=...), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.4.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/Human_MASH_LT1to8_Canonical_v1.3.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_TargetedBroadAnnotationAudit_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 targeted broad annotation audit ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("NO relabeling; NO RDS save; NO UMAP recomputation.")

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  "celltype_ref12_v13",
  "canonical_sample_v13",
  "NAS_score_v13"
)

miss_meta <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss_meta)) {
  stop("Missing metadata: ", paste(miss_meta, collapse = ", "))
}

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")

rna_layers <- Layers(obj[["RNA"]])
if (!all(c("counts", "data") %in% rna_layers)) {
  stop("RNA counts and/or data layer absent. v1.4.1 will not create them.")
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

ann <- as.character(obj$celltype_ref12_v13)
sample_id <- as.character(obj$canonical_sample_v13)
nas <- as.numeric(obj$NAS_score_v13)
cell_ids <- colnames(obj)

# ------------------------------------------------------------------
# Marker sets
# ------------------------------------------------------------------

lineage_sets <- list(
  Hepatocyte = c("ALB", "TTR", "APOA1", "ASGR1", "CPS1", "CYP3A4", "TDO2"),
  Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "KRT8", "MUC1", "KRT17"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "PLVAP", "RGCC"),
  LSEC = c("CLEC4G", "STAB1", "STAB2", "LYVE1", "FCGR2B", "CLEC4M"),
  Mesenchymal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "RGS5"),
  Myeloid = c("LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"),
  T_cell = c("CD3D", "CD3E", "TRAC", "CD247"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "CD37"),
  Plasma_cell = c("MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHM")
)

macrophage_state_markers <- c(
  "C1QA", "C1QB", "C1QC",
  "CD163", "MRC1", "FOLR2", "TIMD4", "MARCO", "VSIG4", "CD5L",
  "SPP1", "TREM2"
)

monocyte_markers <- c(
  "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "CTSS"
)

cycling_markers <- c("MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1")

immature_hep_markers <- c(
  "AFP", "EPCAM", "KRT19", "KRT7", "KRT8", "KRT17"
)

mature_hep_markers <- c(
  "ALB", "TTR", "APOA1", "ASGR1", "CPS1", "CYP3A4", "TDO2"
)

tiny_review_markers <- unique(c(
  lineage_sets$Hepatocyte,
  lineage_sets$Cholangiocyte,
  lineage_sets$Mesenchymal,
  lineage_sets$Endothelial,
  lineage_sets$LSEC,
  lineage_sets$Myeloid,
  "MSLN", "WT1",
  macrophage_state_markers
))

all_requested <- unique(c(
  unlist(lineage_sets, use.names = FALSE),
  macrophage_state_markers,
  monocyte_markers,
  cycling_markers,
  immature_hep_markers,
  mature_hep_markers,
  tiny_review_markers
))

marker_presence <- data.frame(
  gene = all_requested,
  present = all_requested %in% rownames(obj[["RNA"]]),
  stringsAsFactors = FALSE
)

write.csv(
  marker_presence,
  file.path(TABDIR, "01_marker_presence_v1.4.1.csv"),
  row.names = FALSE
)

present <- function(g) intersect(g, rownames(obj[["RNA"]]))

# ------------------------------------------------------------------
# Per-cell lineage evidence:
#   detected_n = number of marker genes with raw count > 0.
#   mean_data  = mean expression across available markers in existing data layer.
# No arbitrary relabeling is performed.
# ------------------------------------------------------------------

evidence <- data.frame(
  cell_id = cell_ids,
  transferred_annotation = ann,
  canonical_sample = sample_id,
  NAS_score = nas,
  stringsAsFactors = FALSE
)

for (nm in names(lineage_sets)) {
  genes <- present(lineage_sets[[nm]])
  if (!length(genes)) {
    evidence[[paste0(nm, "_detected_n")]] <- NA_integer_
    evidence[[paste0(nm, "_mean_data")]] <- NA_real_
    next
  }

  evidence[[paste0(nm, "_detected_n")]] <-
    as.integer(Matrix::colSums(counts[genes, , drop = FALSE] > 0))

  evidence[[paste0(nm, "_mean_data")]] <-
    as.numeric(Matrix::colMeans(data_mat[genes, , drop = FALSE]))
}

# Do not write all 81k rows unless useful; targeted subsets are written below.

# ------------------------------------------------------------------
# 02 Group-level evidence fractions by transferred annotation
# ------------------------------------------------------------------

ann_levels <- sort(unique(ann))
group_evidence_list <- list()

for (ct in ann_levels) {
  idx <- ann == ct

  for (nm in names(lineage_sets)) {
    det <- evidence[[paste0(nm, "_detected_n")]][idx]
    avg <- evidence[[paste0(nm, "_mean_data")]][idx]

    group_evidence_list[[length(group_evidence_list) + 1L]] <- data.frame(
      transferred_annotation = ct,
      lineage_evidence = nm,
      n_cells = sum(idx),
      n_markers_present = length(present(lineage_sets[[nm]])),
      mean_detected_marker_count = mean(det, na.rm = TRUE),
      fraction_ge1_marker = mean(det >= 1, na.rm = TRUE),
      fraction_ge2_markers = mean(det >= 2, na.rm = TRUE),
      fraction_ge3_markers = mean(det >= 3, na.rm = TRUE),
      mean_marker_data = mean(avg, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

group_evidence <- do.call(rbind, group_evidence_list)

write.csv(
  group_evidence,
  file.path(TABDIR, "02_lineage_evidence_by_transferred_annotation_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 03 Marker-level summaries for targeted groups
# ------------------------------------------------------------------

target_groups <- c(
  "Inflammatory_monocyte_macrophage",
  "SPP1_TREM2_macrophage",
  "Cycling_cell",
  "Hepatoblast_or_immature_hepatocyte",
  "Mesothelial",
  "Cholangiocyte_or_ductal",
  "Fibroblast_mesenchymal"
)

target_markers <- unique(c(
  unlist(lineage_sets, use.names = FALSE),
  macrophage_state_markers,
  monocyte_markers,
  cycling_markers,
  immature_hep_markers
))
target_markers <- present(target_markers)

marker_group_list <- list()

for (ct in target_groups) {
  idx <- ann == ct
  if (!any(idx)) next

  pct <- Matrix::rowMeans(counts[target_markers, idx, drop = FALSE] > 0) * 100
  avg <- Matrix::rowMeans(data_mat[target_markers, idx, drop = FALSE])

  marker_group_list[[ct]] <- data.frame(
    transferred_annotation = ct,
    gene = target_markers,
    pct_detected = as.numeric(pct[target_markers]),
    mean_data = as.numeric(avg[target_markers]),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}

marker_group_summary <- do.call(rbind, marker_group_list)
rownames(marker_group_summary) <- NULL

write.csv(
  marker_group_summary,
  file.path(TABDIR, "03_target_group_marker_summary_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 04 Inflammatory_monocyte_macrophage: immune-lineage heterogeneity
# ------------------------------------------------------------------

inflam_idx <- ann == "Inflammatory_monocyte_macrophage"

inflam_cols <- c(
  "cell_id", "canonical_sample", "NAS_score",
  "Myeloid_detected_n", "T_cell_detected_n", "NK_cell_detected_n",
  "B_cell_detected_n", "Plasma_cell_detected_n",
  "Myeloid_mean_data", "T_cell_mean_data", "NK_cell_mean_data",
  "B_cell_mean_data", "Plasma_cell_mean_data"
)

inflam_evidence <- evidence[inflam_idx, inflam_cols, drop = FALSE]

write.csv(
  inflam_evidence,
  file.path(TABDIR, "04_inflammatory_mm_cell_lineage_evidence_v1.4.1.csv"),
  row.names = FALSE
)

immune_lineages <- c("Myeloid", "T_cell", "NK_cell", "B_cell", "Plasma_cell")

inflam_summary <- do.call(
  rbind,
  lapply(immune_lineages, function(nm) {
    det <- evidence[[paste0(nm, "_detected_n")]][inflam_idx]
    data.frame(
      lineage = nm,
      n_cells = sum(inflam_idx),
      mean_detected_n = mean(det, na.rm = TRUE),
      pct_cells_ge1 = mean(det >= 1, na.rm = TRUE) * 100,
      pct_cells_ge2 = mean(det >= 2, na.rm = TRUE) * 100,
      pct_cells_ge3 = mean(det >= 3, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  inflam_summary,
  file.path(TABDIR, "05_inflammatory_mm_lineage_summary_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 06 SPP1_TREM2 group targeted audit
# ------------------------------------------------------------------

spp_idx <- ann == "SPP1_TREM2_macrophage"
spp_genes <- present(unique(c(
  lineage_sets$Myeloid,
  macrophage_state_markers,
  lineage_sets$Mesenchymal,
  lineage_sets$Hepatocyte
)))

spp_summary <- data.frame()

if (any(spp_idx)) {
  spp_pct <- Matrix::rowMeans(counts[spp_genes, spp_idx, drop = FALSE] > 0) * 100
  spp_avg <- Matrix::rowMeans(data_mat[spp_genes, spp_idx, drop = FALSE])

  spp_summary <- data.frame(
    gene = spp_genes,
    pct_detected = as.numeric(spp_pct[spp_genes]),
    mean_data = as.numeric(spp_avg[spp_genes]),
    stringsAsFactors = FALSE
  )
}

write.csv(
  spp_summary,
  file.path(TABDIR, "06_SPP1_TREM2_group_marker_audit_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 07 Cycling-cell lineage background
# ------------------------------------------------------------------

cyc_idx <- ann == "Cycling_cell"
background_lineages <- c(
  "Hepatocyte", "Cholangiocyte", "Endothelial", "LSEC",
  "Mesenchymal", "Myeloid", "T_cell", "NK_cell", "B_cell", "Plasma_cell"
)

cycling_summary <- do.call(
  rbind,
  lapply(background_lineages, function(nm) {
    det <- evidence[[paste0(nm, "_detected_n")]][cyc_idx]
    data.frame(
      lineage = nm,
      n_cells = sum(cyc_idx),
      mean_detected_n = mean(det, na.rm = TRUE),
      pct_cells_ge2 = mean(det >= 2, na.rm = TRUE) * 100,
      pct_cells_ge3 = mean(det >= 3, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  cycling_summary,
  file.path(TABDIR, "07_cycling_cell_lineage_background_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 08 Hepatoblast/immature-hepatocyte audit
# ------------------------------------------------------------------

hb_idx <- ann == "Hepatoblast_or_immature_hepatocyte"
hep_audit_genes <- present(unique(c(
  mature_hep_markers,
  immature_hep_markers,
  cycling_markers
)))

hep_audit <- data.frame()

if (any(hb_idx)) {
  hb_pct <- Matrix::rowMeans(counts[hep_audit_genes, hb_idx, drop = FALSE] > 0) * 100
  hb_avg <- Matrix::rowMeans(data_mat[hep_audit_genes, hb_idx, drop = FALSE])

  hep_audit <- data.frame(
    gene = hep_audit_genes,
    pct_detected = as.numeric(hb_pct[hep_audit_genes]),
    mean_data = as.numeric(hb_avg[hep_audit_genes]),
    stringsAsFactors = FALSE
  )
}

write.csv(
  hep_audit,
  file.path(TABDIR, "08_immature_hepatocyte_marker_audit_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 09 Tiny labels: per-cell key-marker expression
# ------------------------------------------------------------------

tiny_labels <- c(
  "Mesothelial",
  "Cholangiocyte_or_ductal",
  "Fibroblast_mesenchymal"
)

tiny_idx <- ann %in% tiny_labels
tiny_genes <- present(tiny_review_markers)

tiny_table <- data.frame()

if (any(tiny_idx)) {
  tiny_cells <- cell_ids[tiny_idx]

  expr <- as.matrix(data_mat[tiny_genes, tiny_cells, drop = FALSE])
  expr <- t(expr)

  tiny_table <- cbind(
    data.frame(
      cell_id = tiny_cells,
      transferred_annotation = ann[tiny_idx],
      canonical_sample = sample_id[tiny_idx],
      NAS_score = nas[tiny_idx],
      stringsAsFactors = FALSE
    ),
    as.data.frame(expr, check.names = FALSE)
  )
}

write.csv(
  tiny_table,
  file.path(TABDIR, "09_tiny_label_per_cell_marker_expression_v1.4.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figure 1: lineage evidence heatmap
# fraction of cells with >=2 detected markers
# ------------------------------------------------------------------

heat_df <- group_evidence
heat_df$annotation <- factor(
  heat_df$transferred_annotation,
  levels = rev(ann_levels)
)
heat_df$lineage <- factor(
  heat_df$lineage_evidence,
  levels = names(lineage_sets)
)

p_heat <- ggplot(
  heat_df,
  aes(
    x = lineage,
    y = annotation,
    fill = fraction_ge2_markers
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * fraction_ge2_markers)),
    size = 2.5
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Cells with\n≥2 markers"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 targeted broad-lineage audit",
    "Cell-level evidence: fraction with ≥2 detected lineage markers"
  )

ggsave(
  file.path(FIGDIR, "01_lineage_evidence_heatmap_v1.4.1.pdf"),
  p_heat,
  width = 11.5,
  height = 7.5
)

# ------------------------------------------------------------------
# Figure 2: targeted-group marker dot plot
# ------------------------------------------------------------------

target_gene_order <- unique(c(
  lineage_sets$Hepatocyte,
  immature_hep_markers,
  cycling_markers,
  lineage_sets$Mesenchymal,
  lineage_sets$Endothelial,
  lineage_sets$LSEC,
  lineage_sets$Myeloid,
  monocyte_markers,
  macrophage_state_markers,
  lineage_sets$T_cell,
  lineage_sets$NK_cell,
  lineage_sets$B_cell,
  lineage_sets$Plasma_cell
))
target_gene_order <- target_gene_order[target_gene_order %in% target_markers]

dot_df <- marker_group_summary
dot_df$gene <- factor(dot_df$gene, levels = rev(target_gene_order))
dot_df$transferred_annotation <- factor(
  dot_df$transferred_annotation,
  levels = target_groups
)

color_upper <- as.numeric(
  quantile(dot_df$mean_data, 0.98, na.rm = TRUE)
)
if (!is.finite(color_upper) || color_upper <= 0) color_upper <- 1

p_dot <- ggplot(
  dot_df,
  aes(
    x = transferred_annotation,
    y = gene,
    size = pct_detected,
    color = mean_data
  )
) +
  geom_point() +
  scale_size(
    range = c(0.2, 6),
    limits = c(0, 100),
    name = "% detected"
  ) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(0, color_upper),
    oob = scales::squish,
    name = "Mean RNA data"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 50, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7.5)
  ) +
  ggtitle(
    "Targeted review of unresolved transferred labels",
    "No relabeling performed"
  )

ggsave(
  file.path(FIGDIR, "02_target_group_marker_dotplot_v1.4.1.pdf"),
  p_dot,
  width = 11.5,
  height = max(13, 0.20 * length(target_gene_order) + 3)
)

# ------------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------------

get_pct_ge2 <- function(group, lineage) {
  z <- group_evidence[
    group_evidence$transferred_annotation == group &
      group_evidence$lineage_evidence == lineage,
    "fraction_ge2_markers"
  ]
  if (!length(z)) return(NA_real_)
  100 * z[[1]]
}

get_gene_pct <- function(group, gene) {
  z <- marker_group_summary[
    marker_group_summary$transferred_annotation == group &
      marker_group_summary$gene == gene,
    "pct_detected"
  ]
  if (!length(z)) return(NA_real_)
  z[[1]]
}

summary_lines <- c(
  paste0("Human MASH LT1-LT8 targeted broad annotation audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY / NO RELABELING / NO RDS SAVE",
  "",
  "INFLAMMATORY_MONOCYTE_MACROPHAGE:",
  sprintf("Cells: %d", sum(inflam_idx)),
  sprintf("Myeloid >=2 markers: %.1f%%", get_pct_ge2("Inflammatory_monocyte_macrophage", "Myeloid")),
  sprintf("T-cell >=2 markers: %.1f%%", get_pct_ge2("Inflammatory_monocyte_macrophage", "T_cell")),
  sprintf("NK-cell >=2 markers: %.1f%%", get_pct_ge2("Inflammatory_monocyte_macrophage", "NK_cell")),
  sprintf("B-cell >=2 markers: %.1f%%", get_pct_ge2("Inflammatory_monocyte_macrophage", "B_cell")),
  sprintf("Plasma-cell >=2 markers: %.1f%%", get_pct_ge2("Inflammatory_monocyte_macrophage", "Plasma_cell")),
  "",
  "SPP1_TREM2_MACROPHAGE:",
  sprintf("Cells: %d", sum(spp_idx)),
  sprintf("Myeloid >=2 markers: %.1f%%", get_pct_ge2("SPP1_TREM2_macrophage", "Myeloid")),
  sprintf("Mesenchymal >=2 markers: %.1f%%", get_pct_ge2("SPP1_TREM2_macrophage", "Mesenchymal")),
  sprintf("Hepatocyte >=2 markers: %.1f%%", get_pct_ge2("SPP1_TREM2_macrophage", "Hepatocyte")),
  sprintf("SPP1 detected: %.1f%%", get_gene_pct("SPP1_TREM2_macrophage", "SPP1")),
  sprintf("TREM2 detected: %.1f%%", get_gene_pct("SPP1_TREM2_macrophage", "TREM2")),
  "",
  "CYCLING_CELL:",
  sprintf("Cells: %d", sum(cyc_idx)),
  sprintf("Hepatocyte >=2 markers: %.1f%%", get_pct_ge2("Cycling_cell", "Hepatocyte")),
  sprintf("Mesenchymal >=2 markers: %.1f%%", get_pct_ge2("Cycling_cell", "Mesenchymal")),
  sprintf("LSEC >=2 markers: %.1f%%", get_pct_ge2("Cycling_cell", "LSEC")),
  sprintf("Myeloid >=2 markers: %.1f%%", get_pct_ge2("Cycling_cell", "Myeloid")),
  "",
  "HEPATOBLAST_OR_IMMATURE_HEPATOCYTE:",
  sprintf("Cells: %d", sum(hb_idx)),
  sprintf("Hepatocyte >=2 markers: %.1f%%", get_pct_ge2("Hepatoblast_or_immature_hepatocyte", "Hepatocyte")),
  sprintf("Cholangiocyte >=2 markers: %.1f%%", get_pct_ge2("Hepatoblast_or_immature_hepatocyte", "Cholangiocyte")),
  sprintf("AFP detected: %.1f%%", get_gene_pct("Hepatoblast_or_immature_hepatocyte", "AFP")),
  sprintf("EPCAM detected: %.1f%%", get_gene_pct("Hepatoblast_or_immature_hepatocyte", "EPCAM")),
  "",
  "TINY LABELS:",
  paste0("Mesothelial cells: ", sum(ann == "Mesothelial")),
  paste0("Cholangiocyte_or_ductal cells: ", sum(ann == "Cholangiocyte_or_ductal")),
  paste0("Fibroblast_mesenchymal cells: ", sum(ann == "Fibroblast_mesenchymal")),
  "",
  "DECISION GATE:",
  "- Do not freeze broad annotations solely from v1.4 group-average DotPlot.",
  "- Use cell-level lineage evidence here to decide whether immune/cycling/tiny labels require reannotation.",
  "- No biological downstream analysis should begin until this decision gate is resolved.",
  "",
  "OUTPUTS:",
  "01_marker_presence_v1.4.1.csv",
  "02_lineage_evidence_by_transferred_annotation_v1.4.1.csv",
  "03_target_group_marker_summary_v1.4.1.csv",
  "04_inflammatory_mm_cell_lineage_evidence_v1.4.1.csv",
  "05_inflammatory_mm_lineage_summary_v1.4.1.csv",
  "06_SPP1_TREM2_group_marker_audit_v1.4.1.csv",
  "07_cycling_cell_lineage_background_v1.4.1.csv",
  "08_immature_hepatocyte_marker_audit_v1.4.1.csv",
  "09_tiny_label_per_cell_marker_expression_v1.4.1.csv",
  "figures/01_lineage_evidence_heatmap_v1.4.1.pdf",
  "figures/02_target_group_marker_dotplot_v1.4.1.pdf"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "10_targeted_audit_summary_v1.4.1.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.4.1 COMPLETE")
message("Output: ", OUTDIR)
message("No labels changed; no RDS saved; no UMAP recomputed.")
message("============================================================")
