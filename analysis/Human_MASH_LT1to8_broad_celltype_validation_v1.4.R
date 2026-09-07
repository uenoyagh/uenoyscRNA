#!/usr/bin/env Rscript

# Human_MASH_LT1to8_broad_celltype_validation_v1.4.R
#
# PURPOSE
# -------
# Validate the transferred broad whole-liver annotation in
# Human_MASH_LT1to8_Canonical_v1.3 using the restored 16,424-gene RNA universe.
#
# This is a VALIDATION / AUDIT step only.
#
# It DOES:
#   - inspect annotation counts
#   - audit availability of broad liver lineage markers
#   - summarize marker detection and average expression by transferred celltype
#   - visualize the transferred frozen REF12 UMAP
#   - create a broad-lineage marker dot plot
#
# It DOES NOT:
#   - change any cell label
#   - save a modified RDS
#   - normalize / SCTransform
#   - integrate
#   - cluster / recluster
#   - calculate a new UMAP
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#
# Seurat v5:
#   GetAssayData(..., layer=...), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.4"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/Human_MASH_LT1to8_Canonical_v1.3.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_BroadCelltypeValidation_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) {
  stop("Input canonical RDS not found:\n", INPUT_RDS)
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 broad celltype validation ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("Output: ", OUTDIR)
message("NO label modification; NO new UMAP.")

obj <- readRDS(INPUT_RDS)

if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.")
if (!"RNA" %in% Assays(obj)) stop("RNA assay is absent.")
if (!"celltype_ref12_v13" %in% colnames(obj[[]])) {
  stop("celltype_ref12_v13 metadata is absent.")
}
if (!"canonical_sample_v13" %in% colnames(obj[[]])) {
  stop("canonical_sample_v13 metadata is absent.")
}
if (!"NAS_score_v13" %in% colnames(obj[[]])) {
  stop("NAS_score_v13 metadata is absent.")
}
if (!"umapRPCA_ref12_v13" %in% Reductions(obj)) {
  stop("Transferred frozen UMAP 'umapRPCA_ref12_v13' is absent.")
}

annotation <- as.character(obj$celltype_ref12_v13)
sample_id <- as.character(obj$canonical_sample_v13)
nas <- as.numeric(obj$NAS_score_v13)

if (anyNA(annotation) || any(annotation == "")) {
  stop("Transferred annotation contains missing/blank values.")
}

# ------------------------------------------------------------------
# Broad lineage marker panel.
# This panel is intentionally redundant: validation should rely on marker sets,
# not on one marker alone.
# ------------------------------------------------------------------

marker_sets <- list(
  Hepatocyte = c(
    "ALB", "TTR", "APOA1", "ASGR1", "CPS1", "CYP3A4", "TDO2"
  ),
  Cholangiocyte = c(
    "KRT19", "KRT7", "EPCAM", "KRT8", "KRT18", "MUC1", "KRT17"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "EMCN", "ENG", "PLVAP", "RGCC"
  ),
  LSEC_enriched = c(
    "CLEC4G", "STAB1", "STAB2", "LYVE1", "FCGR2B", "CLEC4M"
  ),
  HSC_Mesenchymal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "RGS5"
  ),
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Macrophage = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2", "TIMD4",
    "MARCO", "VSIG4", "TREM2", "SPP1", "CD5L"
  ),
  Monocyte = c(
    "FCN1", "S100A8", "S100A9", "CCR2", "LILRB1", "CTSS", "MNDA"
  ),
  Neutrophil = c(
    "FCGR3B", "CSF3R", "S100A8", "S100A9", "MPO", "ELANE", "CXCR2"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  T_cell = c(
    "CD3D", "CD3E", "TRAC", "CD247", "IL7R"
  ),
  NK_cell = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1"
  ),
  B_cell = c(
    "MS4A1", "CD79A", "CD79B", "CD74", "CD37"
  ),
  Plasma_cell = c(
    "MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHM"
  ),
  Mesothelial = c(
    "MSLN", "KRT19", "KRT8", "KRT18", "WT1"
  ),
  Platelet = c(
    "PPBP", "PF4", "GP9", "GP6", "ITGA2B"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1"
  )
)

marker_manifest <- do.call(
  rbind,
  lapply(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    data.frame(
      marker_set = set_name,
      gene = genes,
      present = genes %in% rownames(obj[["RNA"]]),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_manifest,
  file.path(TABDIR, "01_marker_presence_manifest_v1.4.csv"),
  row.names = FALSE
)

present_markers <- unique(marker_manifest$gene[marker_manifest$present])
missing_markers <- unique(marker_manifest$gene[!marker_manifest$present])

message("Markers present: ", length(present_markers))
message("Markers missing: ", length(missing_markers))

# ------------------------------------------------------------------
# Annotation inventory.
# ------------------------------------------------------------------

annotation_counts <- as.data.frame(
  table(celltype_ref12_v13 = annotation, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(annotation_counts)[2] <- "n_cells"
annotation_counts <- annotation_counts[
  order(annotation_counts$n_cells, decreasing = TRUE),
  , drop = FALSE
]

write.csv(
  annotation_counts,
  file.path(TABDIR, "02_transferred_annotation_counts_v1.4.csv"),
  row.names = FALSE
)

annotation_by_sample <- as.data.frame(
  table(
    canonical_sample_v13 = sample_id,
    celltype_ref12_v13 = annotation,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
annotation_by_sample <- annotation_by_sample[
  annotation_by_sample$Freq > 0,
  , drop = FALSE
]
colnames(annotation_by_sample)[3] <- "n_cells"

write.csv(
  annotation_by_sample,
  file.path(TABDIR, "03_annotation_counts_by_sample_v1.4.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# RNA layer audit.
# We use an EXISTING normalized 'data' layer if present.
# We never create/normalize a layer in v1.4.
# counts are used independently for percent detected.
# ------------------------------------------------------------------

rna_layers <- Layers(obj[["RNA"]])

layer_audit <- data.frame(
  assay = "RNA",
  layers = paste(rna_layers, collapse = " | "),
  data_layer_present = "data" %in% rna_layers,
  counts_layer_present = "counts" %in% rna_layers,
  stringsAsFactors = FALSE
)

write.csv(
  layer_audit,
  file.path(TABDIR, "04_RNA_layer_audit_v1.4.csv"),
  row.names = FALSE
)

if (!"counts" %in% rna_layers) {
  stop("RNA counts layer is absent.")
}
if (!"data" %in% rna_layers) {
  stop(
    "RNA data layer is absent. v1.4 intentionally does not normalize the object.\n",
    "Stop here rather than creating a new normalized representation."
  )
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

counts_m <- counts[present_markers, , drop = FALSE]
data_m <- data_mat[present_markers, , drop = FALSE]

# ------------------------------------------------------------------
# Marker summary by transferred annotation.
# pct_detected = percentage of cells with raw count > 0.
# avg_data = mean expression in the EXISTING RNA data layer.
# ------------------------------------------------------------------

celltypes <- unique(annotation)
celltypes <- celltypes[order(celltypes)]

summary_list <- vector("list", length(celltypes))

for (i in seq_along(celltypes)) {
  ct <- celltypes[[i]]
  idx <- which(annotation == ct)

  pct <- Matrix::rowMeans(counts_m[, idx, drop = FALSE] > 0) * 100
  avg <- Matrix::rowMeans(data_m[, idx, drop = FALSE])

  summary_list[[i]] <- data.frame(
    celltype_ref12_v13 = ct,
    gene = names(pct),
    pct_detected = as.numeric(pct),
    avg_data = as.numeric(avg[names(pct)]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, summary_list)

gene_to_set <- marker_manifest[
  marker_manifest$present,
  c("gene", "marker_set"),
  drop = FALSE
]

# A gene can occur in more than one conceptual set.
marker_summary_long <- merge(
  marker_summary,
  gene_to_set,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  marker_summary_long,
  file.path(TABDIR, "05_marker_expression_by_annotation_v1.4.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Set-level summary.
# This is descriptive only; NO automated relabeling is performed.
# ------------------------------------------------------------------

set_summary_list <- list()

for (set_name in names(marker_sets)) {
  genes <- intersect(marker_sets[[set_name]], present_markers)
  if (!length(genes)) next

  for (ct in celltypes) {
    z <- marker_summary[
      marker_summary$celltype_ref12_v13 == ct &
        marker_summary$gene %in% genes,
      , drop = FALSE
    ]

    set_summary_list[[length(set_summary_list) + 1L]] <- data.frame(
      celltype_ref12_v13 = ct,
      marker_set = set_name,
      n_markers_present = length(genes),
      mean_pct_detected = mean(z$pct_detected),
      mean_avg_data = mean(z$avg_data),
      stringsAsFactors = FALSE
    )
  }
}

set_summary <- do.call(rbind, set_summary_list)

write.csv(
  set_summary,
  file.path(TABDIR, "06_marker_set_summary_by_annotation_v1.4.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Frozen UMAP visualizations.
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umapRPCA_ref12_v13")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  celltype = annotation,
  sample = sample_id,
  NAS = factor(nas, levels = c(0, 2, 4, 6)),
  stringsAsFactors = FALSE
)

# High-saturation palette; deterministic, annotation-count agnostic.
base_palette <- c(
  "#00A6D6", "#009E73", "#0072B2", "#D55E00", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0007A", "#00C853", "#7B1FA2",
  "#F4511E", "#3949AB", "#00897B", "#C0CA33", "#6D4C41",
  "#EF5350", "#26A69A", "#AB47BC", "#5C6BC0", "#FFA000"
)

if (length(celltypes) > length(base_palette)) {
  palette_values <- grDevices::hcl.colors(
    length(celltypes),
    palette = "Dynamic"
  )
} else {
  palette_values <- base_palette[seq_along(celltypes)]
}
names(palette_values) <- celltypes

p_umap <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = celltype)
) +
  geom_point(size = 0.22, alpha = 0.82) +
  scale_color_manual(values = palette_values) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.key.height = grid::unit(0.42, "cm")
  ) +
  ggtitle(
    "Human MASH LT1-LT8 | transferred broad annotation\n",
    "Frozen REF12 UMAP coordinates; no UMAP recomputation"
  )

ggsave(
  file.path(FIGDIR, "01_whole_UMAP_transferred_annotation_v1.4.pdf"),
  p_umap,
  width = 10.5,
  height = 8.2
)

p_sample <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = celltype)
) +
  geom_point(size = 0.16, alpha = 0.82) +
  scale_color_manual(values = palette_values) +
  coord_equal() +
  facet_wrap(~ sample, ncol = 4) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  ggtitle("Human MASH LT1-LT8 | transferred annotation by biological sample")

ggsave(
  file.path(FIGDIR, "02_UMAP_by_sample_transferred_annotation_v1.4.pdf"),
  p_sample,
  width = 14,
  height = 9.5
)

# ------------------------------------------------------------------
# Broad marker dot plot.
# Color = existing normalized RNA data layer mean.
# Size  = raw-count detection percentage.
# No scaling across groups is applied.
# ------------------------------------------------------------------

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% present_markers]

dot_df <- marker_summary[
  marker_summary$gene %in% marker_order,
  , drop = FALSE
]

dot_df$gene <- factor(dot_df$gene, levels = rev(marker_order))
dot_df$celltype_ref12_v13 <- factor(
  dot_df$celltype_ref12_v13,
  levels = celltypes
)

# Robust upper cutoff only for visualization; raw summary remains in CSV.
color_upper <- as.numeric(
  stats::quantile(dot_df$avg_data, 0.98, na.rm = TRUE)
)
if (!is.finite(color_upper) || color_upper <= 0) {
  color_upper <- max(dot_df$avg_data, na.rm = TRUE)
}
if (!is.finite(color_upper) || color_upper <= 0) {
  color_upper <- 1
}

p_dot <- ggplot(
  dot_df,
  aes(
    x = celltype_ref12_v13,
    y = gene,
    size = pct_detected,
    color = avg_data
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
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 7.5)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 broad-lineage marker validation",
    "Transferred celltype labels | restored 16,424-gene RNA universe"
  )

ggsave(
  file.path(FIGDIR, "03_broad_lineage_marker_dotplot_v1.4.pdf"),
  p_dot,
  width = max(11, 0.65 * length(celltypes) + 5),
  height = max(12, 0.22 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

summary_lines <- c(
  paste0("Human MASH LT1-LT8 broad celltype validation ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "INPUT:",
  INPUT_RDS,
  "",
  "OBJECT:",
  paste0("Cells: ", ncol(obj)),
  paste0("RNA features: ", nrow(obj[["RNA"]])),
  paste0("RNA layers: ", paste(rna_layers, collapse = " | ")),
  paste0("Transferred annotation categories: ", length(celltypes)),
  paste0("Transferred UMAP reduction: umapRPCA_ref12_v13"),
  "",
  "MARKERS:",
  paste0("Marker genes requested: ", length(unique(marker_manifest$gene))),
  paste0("Marker genes present: ", length(present_markers)),
  paste0("Marker genes missing: ", length(missing_markers)),
  paste0(
    "Missing marker genes: ",
    if (length(missing_markers)) {
      paste(missing_markers, collapse = " | ")
    } else {
      "<none>"
    }
  ),
  "",
  "ANNOTATION COUNTS:",
  paste(
    capture.output(print(annotation_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION RULE:",
  "- v1.4 does NOT relabel cells.",
  "- Evaluate whether each transferred celltype is supported by a coherent positive marker set",
  "  and lacks strong incompatible-lineage signal.",
  "- Any ambiguous/mixed group should be flagged for targeted review before annotation freeze.",
  "",
  "NEXT STEP:",
  "Review 03_broad_lineage_marker_dotplot_v1.4.pdf plus marker-set CSVs.",
  "Only after review should a broad-celltype annotation be frozen.",
  "",
  "SOURCE RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "07_validation_summary_v1.4.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.4 COMPLETE")
message("Transferred annotation categories: ", length(celltypes))
message("Markers present: ", length(present_markers))
message("Markers missing: ", length(missing_markers))
message("Output: ", OUTDIR)
message("No labels changed; no RDS saved; no UMAP recomputed.")
message("============================================================")
