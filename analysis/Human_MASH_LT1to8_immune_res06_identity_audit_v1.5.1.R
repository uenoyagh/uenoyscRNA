#!/usr/bin/env Rscript

# Human_MASH_LT1to8_immune_res06_identity_audit_v1.5.1.R
#
# PURPOSE
# -------
# Audit biological identity of the v1.5 immune RPCA clusters at the working
# resolution 0.6 before any annotation is frozen.
#
# This is a READ-ONLY review step:
#   - NO reclustering
#   - NO UMAP recalculation
#   - NO cell relabeling
#   - NO RDS save
#   - NO DE / pseudobulk / pathway analysis
#
# Outputs provide:
#   1) cluster-level lineage evidence
#   2) key-marker detection/expression
#   3) sample balance
#   4) transferred source-label composition
#   5) marker heatmap and dot plot
#
# Important:
# ALB/TTR/APOA1 are highly abundant liver transcripts and may reflect ambient
# RNA in immune droplets. Therefore the "Hepatocyte_strict" contamination set
# deliberately emphasizes ASGR1/CPS1/CYP3A4/TDO2 rather than ALB alone.
#
# Seurat v5:
#   GetAssayData(..., layer=...), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.5.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5/",
  "objects/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Immune_res06_IdentityAudit_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.5 immune RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 immune res0.6 identity audit ", VERSION)
message("============================================================")
message("READ-ONLY: no clustering, no relabeling, no RDS save.")

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "immune_RPCA_res0_6"

required_meta <- c(
  CLUSTER_COL,
  "sample_v15",
  "source_broad_label_v15"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")

rna_layers <- Layers(obj[["RNA"]])
if (!all(c("counts", "data") %in% rna_layers)) {
  stop("RNA counts/data layers are required.")
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

cl <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
sample_id <- as.character(obj$sample_v15)
source_label <- as.character(obj$source_broad_label_v15)

clusters <- sort(unique(cl))
features <- rownames(obj[["RNA"]])

# ------------------------------------------------------------------
# Marker sets
# ------------------------------------------------------------------

marker_sets <- list(
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Monocyte_inflammatory = c(
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "FCGR3A"
  ),
  Resident_macrophage = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2",
    "TIMD4", "MARCO", "VSIG4", "CD5L"
  ),
  SPP1_TREM2_state = c(
    "SPP1", "TREM2", "GPNMB", "LGALS3", "LPL", "FABP5", "APOE"
  ),
  T_cell = c(
    "CD3D", "CD3E", "TRAC", "CD247", "IL7R"
  ),
  NK_cell = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1", "GZMB", "PRF1"
  ),
  B_cell = c(
    "MS4A1", "CD79A", "CD79B", "CD37", "CD74", "CD19"
  ),
  Plasma_cell = c(
    "MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHM"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  Mesenchymal_strict = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  ),
  Hepatocyte_strict = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
  ),
  Ambient_hepatocyte = c(
    "ALB", "TTR", "APOA1"
  )
)

marker_manifest <- do.call(
  rbind,
  lapply(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    data.frame(
      marker_set = set_name,
      gene = genes,
      present = genes %in% features,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_manifest,
  file.path(TABDIR, "01_marker_presence_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Cluster-level marker-set evidence.
#
# For each marker set:
#   mean_pct_detected = average % of cells positive for each marker
#   frac_cells_ge2    = fraction of cells with >=2 markers detected
#   mean_marker_data  = mean normalized RNA-data signal
# ------------------------------------------------------------------

set_evidence_list <- list()

for (cluster_id in clusters) {
  idx <- which(cl == cluster_id)

  for (set_name in names(marker_sets)) {
    genes <- intersect(marker_sets[[set_name]], features)

    if (!length(genes)) next

    sub_counts <- counts[genes, idx, drop = FALSE]
    sub_data <- data_mat[genes, idx, drop = FALSE]

    per_cell_detected <- Matrix::colSums(sub_counts > 0)
    pct_each_gene <- Matrix::rowMeans(sub_counts > 0) * 100

    set_evidence_list[[length(set_evidence_list) + 1L]] <- data.frame(
      cluster = cluster_id,
      marker_set = set_name,
      n_cells = length(idx),
      n_markers_present = length(genes),
      mean_pct_detected = mean(pct_each_gene),
      frac_cells_ge1 = mean(per_cell_detected >= 1),
      frac_cells_ge2 = mean(per_cell_detected >= 2),
      frac_cells_ge3 = mean(per_cell_detected >= 3),
      mean_marker_data = mean(Matrix::rowMeans(sub_data)),
      stringsAsFactors = FALSE
    )
  }
}

set_evidence <- do.call(rbind, set_evidence_list)

write.csv(
  set_evidence,
  file.path(TABDIR, "02_res06_cluster_marker_set_evidence_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Key-marker summary for direct biological inspection.
# ------------------------------------------------------------------

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

marker_summary_list <- list()

for (cluster_id in clusters) {
  idx <- which(cl == cluster_id)

  pct <- Matrix::rowMeans(counts[marker_order, idx, drop = FALSE] > 0) * 100
  avg <- Matrix::rowMeans(data_mat[marker_order, idx, drop = FALSE])

  marker_summary_list[[cluster_id]] <- data.frame(
    cluster = cluster_id,
    gene = marker_order,
    pct_detected = as.numeric(pct[marker_order]),
    mean_data = as.numeric(avg[marker_order]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, marker_summary_list)
rownames(marker_summary) <- NULL

write.csv(
  marker_summary,
  file.path(TABDIR, "03_res06_key_marker_summary_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Sample balance.
# ------------------------------------------------------------------

sample_tab <- as.data.frame(
  table(cluster = cl, sample = sample_id),
  stringsAsFactors = FALSE
)
sample_tab <- sample_tab[sample_tab$Freq > 0, , drop = FALSE]
colnames(sample_tab)[3] <- "n_cells"

write.csv(
  sample_tab,
  file.path(TABDIR, "04_res06_cluster_by_sample_v1.5.1.csv"),
  row.names = FALSE
)

sample_balance_list <- list()

for (cluster_id in clusters) {
  idx <- cl == cluster_id
  tt <- table(sample_id[idx])
  p <- as.numeric(tt) / sum(tt)
  h <- -sum(p * log(p))
  h_norm <- if (length(tt) > 1) h / log(8) else 0

  sample_balance_list[[cluster_id]] <- data.frame(
    cluster = cluster_id,
    n_cells = sum(idx),
    n_samples_present = length(tt),
    dominant_sample = names(tt)[which.max(tt)],
    dominant_sample_fraction = max(p),
    normalized_sample_entropy = h_norm,
    stringsAsFactors = FALSE
  )
}

sample_balance <- do.call(rbind, sample_balance_list)
rownames(sample_balance) <- NULL

write.csv(
  sample_balance,
  file.path(TABDIR, "05_res06_sample_balance_summary_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Source-label composition.
# ------------------------------------------------------------------

source_tab <- as.data.frame(
  table(cluster = cl, source_label = source_label),
  stringsAsFactors = FALSE
)
source_tab <- source_tab[source_tab$Freq > 0, , drop = FALSE]
colnames(source_tab)[3] <- "n_cells"

source_total <- aggregate(n_cells ~ cluster, source_tab, sum)
colnames(source_total)[2] <- "cluster_total"
source_tab <- merge(source_tab, source_total, by = "cluster", all.x = TRUE)
source_tab$fraction <- source_tab$n_cells / source_tab$cluster_total

write.csv(
  source_tab,
  file.path(TABDIR, "06_res06_source_label_composition_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Screening rank for review ONLY.
# It is NOT written back to the Seurat object.
#
# Broad immune lineages use frac_cells_ge2. Myeloid state scores are kept
# separately and do not compete with broad T/NK/B/Plasma/DC identity.
# ------------------------------------------------------------------

BROAD_SETS <- c(
  "Myeloid_core",
  "T_cell",
  "NK_cell",
  "B_cell",
  "Plasma_cell",
  "Dendritic",
  "Mesenchymal_strict",
  "Hepatocyte_strict"
)

screen_list <- list()

for (cluster_id in clusters) {
  z <- set_evidence[
    set_evidence$cluster == cluster_id &
      set_evidence$marker_set %in% BROAD_SETS,
    c("marker_set", "frac_cells_ge2"),
    drop = FALSE
  ]

  z <- z[order(z$frac_cells_ge2, decreasing = TRUE), , drop = FALSE]

  top1 <- z[1, , drop = FALSE]
  top2 <- z[2, , drop = FALSE]

  get_state <- function(set_name) {
    zz <- set_evidence[
      set_evidence$cluster == cluster_id &
        set_evidence$marker_set == set_name,
      "frac_cells_ge2"
    ]
    if (!length(zz)) NA_real_ else zz[[1]]
  }

  screen_list[[cluster_id]] <- data.frame(
    cluster = cluster_id,
    n_cells = sum(cl == cluster_id),
    evidence_rank1 = top1$marker_set,
    rank1_frac_ge2 = top1$frac_cells_ge2,
    evidence_rank2 = top2$marker_set,
    rank2_frac_ge2 = top2$frac_cells_ge2,
    rank_margin = top1$frac_cells_ge2 - top2$frac_cells_ge2,
    resident_macrophage_frac_ge2 = get_state("Resident_macrophage"),
    monocyte_inflammatory_frac_ge2 = get_state("Monocyte_inflammatory"),
    SPP1_TREM2_state_frac_ge2 = get_state("SPP1_TREM2_state"),
    ambient_hepatocyte_frac_ge2 = get_state("Ambient_hepatocyte"),
    stringsAsFactors = FALSE
  )
}

screen <- do.call(rbind, screen_list)
rownames(screen) <- NULL

screen <- merge(
  screen,
  sample_balance,
  by = c("cluster", "n_cells"),
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  screen,
  file.path(TABDIR, "07_res06_identity_screening_rank_v1.5.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

# Heatmap: fraction of cells with >=2 markers.
heat_df <- set_evidence[
  set_evidence$marker_set %in% BROAD_SETS,
  , drop = FALSE
]

heat_df$cluster <- factor(heat_df$cluster, levels = clusters)
heat_df$marker_set <- factor(heat_df$marker_set, levels = BROAD_SETS)

p_heat <- ggplot(
  heat_df,
  aes(x = cluster, y = marker_set, fill = frac_cells_ge2)
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * frac_cells_ge2)),
    size = 2.6
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Cells with\n>=2 markers"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 immune RPCA res0.6 | broad-lineage evidence",
    "Screening only; no cluster labels assigned"
  )

ggsave(
  file.path(FIGDIR, "01_res06_broad_lineage_evidence_heatmap_v1.5.1.pdf"),
  p_heat,
  width = max(11, 0.55 * length(clusters) + 4),
  height = 6.5
)

# Dot plot with broad + myeloid state markers.
dot_df <- marker_summary
dot_df$cluster <- factor(dot_df$cluster, levels = clusters)
dot_df$gene <- factor(dot_df$gene, levels = rev(marker_order))

upper <- as.numeric(quantile(dot_df$mean_data, 0.98, na.rm = TRUE))
if (!is.finite(upper) || upper <= 0) upper <- 1

p_dot <- ggplot(
  dot_df,
  aes(
    x = cluster,
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
    limits = c(0, upper),
    oob = scales::squish,
    name = "Mean RNA data"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7.5)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 immune RPCA res0.6 | identity audit",
    "Ambient hepatocyte markers shown separately; no relabeling"
  )

ggsave(
  file.path(FIGDIR, "02_res06_identity_marker_dotplot_v1.5.1.pdf"),
  p_dot,
  width = max(12, 0.58 * length(clusters) + 5),
  height = max(13, 0.22 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------------

small_or_sample_limited <- sample_balance[
  sample_balance$n_cells < 100 |
    sample_balance$n_samples_present < 8 |
    sample_balance$dominant_sample_fraction > 0.5,
  , drop = FALSE
]

summary_lines <- c(
  paste0("Human MASH LT1-LT8 immune res0.6 identity audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no reclustering, no relabeling, no RDS save.",
  "",
  paste0("Input immune cells: ", ncol(obj)),
  paste0("RNA features: ", nrow(obj[["RNA"]])),
  paste0("res0.6 clusters: ", length(clusters)),
  "",
  "SAMPLE-BALANCE FLAGS:",
  if (nrow(small_or_sample_limited)) {
    paste(
      capture.output(print(small_or_sample_limited, row.names = FALSE)),
      collapse = "\n"
    )
  } else {
    "<none>"
  },
  "",
  "IDENTITY SCREENING RANK:",
  paste(
    capture.output(print(screen, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION NOTES:",
  "- evidence_rank1/rank2 are screening ranks only, not biological labels.",
  "- ALB/TTR/APOA1 are treated as ambient-prone and do not define strict hepatocyte contamination.",
  "- Resident_macrophage / Monocyte_inflammatory / SPP1_TREM2_state are state evidence within myeloid clusters.",
  "- Mesenchymal_strict is included specifically to identify non-immune collagen/DCN/LUM-rich contamination.",
  "",
  "NEXT DECISION GATE:",
  "Use the CSVs and two figures to assign broad immune identities to res0.6 clusters.",
  "Only then create a clean macrophage/monocyte object for detailed human MASH analysis."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "08_identity_audit_summary_v1.5.1.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.5.1 COMPLETE")
message("res0.6 clusters: ", length(clusters))
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
