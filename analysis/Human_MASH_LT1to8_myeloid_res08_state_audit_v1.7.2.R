#!/usr/bin/env Rscript

# Human_MASH_LT1to8_myeloid_res08_state_audit_v1.7.2.R
#
# PURPOSE
# -------
# Use res0.8 as the WORKING (not yet frozen) myeloid resolution and perform
# a read-only preliminary Human MΦ / Monocyte state audit.
#
# WHY res0.8
# ----------
# v1.7.1 showed that res0.8 adds only one net cluster relative to res0.6
# (17 vs 16), while reorganizing the effective res0.6 splits (0, 3, 8)
# into more coherent transcriptional partitions. The goal here is to test
# whether those res0.8 clusters map to interpretable myeloid states.
#
# THIS SCRIPT DOES:
#   - audit all res0.8 clusters
#   - quantify broad myeloid identity and contamination
#   - quantify macrophage/monocyte state marker sets
#   - flag DC-like / cycling / contamination clusters
#   - summarize sample and NAS composition
#   - create numeric-cluster UMAP, marker DotPlot, and state-evidence heatmap
#
# THIS SCRIPT DOES NOT:
#   - recluster
#   - recompute UMAP
#   - modify metadata
#   - save/overwrite an RDS
#   - freeze biological labels
#   - perform condition DE / pseudobulk / pathway / cell-cell analysis
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.7.2"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_RPCA_v1.7/",
  "objects/Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_res08_StateAudit_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.7 RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 myeloid res0.8 state audit ", VERSION)
message("============================================================")
message("READ-ONLY: no reclustering, no relabeling, no RDS save.")

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "myeloid_RPCA_res0_8"

required_meta <- c(
  CLUSTER_COL,
  "sample_v17",
  "NAS_v17",
  "parent_cluster_v17"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) {
  stop("RNA counts/data layers required.")
}
if (!"umap.myeloid.rpca" %in% Reductions(obj)) {
  stop("umap.myeloid.rpca reduction absent.")
}

cl <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
sample_id <- as.character(obj$sample_v17)
nas <- as.numeric(obj$NAS_v17)
parent_cluster <- as.character(obj$parent_cluster_v17)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
features <- rownames(obj[["RNA"]])

clusters <- sort(unique(cl))

# ------------------------------------------------------------------
# Marker sets
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
  ),

  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),

  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1"
  ),

  T_NK_contamination = c(
    "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"
  ),

  B_Plasma_contamination = c(
    "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"
  ),

  Mesenchymal_contamination = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  ),

  Hepatocyte_strict_contamination = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
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
  file.path(TABDIR, "01_marker_presence_v1.7.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Cluster inventory / sample / NAS / parent provenance
# ------------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(cluster_res0_8 = cl),
  stringsAsFactors = FALSE
)
colnames(cluster_counts)[2] <- "n_cells"

write.csv(
  cluster_counts,
  file.path(TABDIR, "02_cluster_counts_v1.7.2.csv"),
  row.names = FALSE
)

cluster_sample <- as.data.frame(
  table(cluster_res0_8 = cl, sample = sample_id),
  stringsAsFactors = FALSE
)
cluster_sample <- cluster_sample[cluster_sample$Freq > 0, , drop = FALSE]
colnames(cluster_sample)[3] <- "n_cells"

write.csv(
  cluster_sample,
  file.path(TABDIR, "03_cluster_by_sample_v1.7.2.csv"),
  row.names = FALSE
)

cluster_nas <- as.data.frame(
  table(cluster_res0_8 = cl, NAS = nas),
  stringsAsFactors = FALSE
)
cluster_nas <- cluster_nas[cluster_nas$Freq > 0, , drop = FALSE]
colnames(cluster_nas)[3] <- "n_cells"

cluster_tot <- aggregate(n_cells ~ cluster_res0_8, cluster_nas, sum)
colnames(cluster_tot)[2] <- "cluster_total"
cluster_nas <- merge(cluster_nas, cluster_tot, by = "cluster_res0_8")
cluster_nas$fraction <- cluster_nas$n_cells / cluster_nas$cluster_total

write.csv(
  cluster_nas,
  file.path(TABDIR, "04_cluster_by_NAS_descriptive_v1.7.2.csv"),
  row.names = FALSE
)

cluster_parent <- as.data.frame(
  table(cluster_res0_8 = cl, parent_v1_5_cluster = parent_cluster),
  stringsAsFactors = FALSE
)
cluster_parent <- cluster_parent[cluster_parent$Freq > 0, , drop = FALSE]
colnames(cluster_parent)[3] <- "n_cells"

write.csv(
  cluster_parent,
  file.path(TABDIR, "05_cluster_by_parent_provenance_v1.7.2.csv"),
  row.names = FALSE
)

# Sample-balance summary
balance_list <- list()

for (cluster_id in clusters) {
  idx <- cl == cluster_id
  tt <- table(sample_id[idx])
  p <- as.numeric(tt) / sum(tt)
  h <- -sum(p * log(p))
  h_norm <- if (length(tt) > 1) h / log(8) else 0

  balance_list[[cluster_id]] <- data.frame(
    cluster_res0_8 = cluster_id,
    n_cells = sum(idx),
    n_samples_present = length(tt),
    dominant_sample = names(tt)[which.max(tt)],
    dominant_sample_fraction = max(p),
    normalized_sample_entropy = h_norm,
    stringsAsFactors = FALSE
  )
}

balance <- do.call(rbind, balance_list)
rownames(balance) <- NULL

write.csv(
  balance,
  file.path(TABDIR, "06_sample_balance_v1.7.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Marker-set evidence by cluster
# ------------------------------------------------------------------

set_summary_list <- list()

for (cluster_id in clusters) {
  idx <- cl == cluster_id

  for (set_name in names(marker_sets)) {
    genes <- intersect(marker_sets[[set_name]], features)
    if (!length(genes)) next

    sub_counts <- counts[genes, idx, drop = FALSE]
    sub_data <- data_mat[genes, idx, drop = FALSE]

    per_cell_n <- Matrix::colSums(sub_counts > 0)
    pct_each_gene <- Matrix::rowMeans(sub_counts > 0) * 100

    set_summary_list[[length(set_summary_list) + 1L]] <- data.frame(
      cluster_res0_8 = cluster_id,
      marker_set = set_name,
      n_cells = sum(idx),
      n_markers_present = length(genes),
      mean_pct_detected = mean(pct_each_gene),
      frac_cells_ge1 = mean(per_cell_n >= 1),
      frac_cells_ge2 = mean(per_cell_n >= 2),
      frac_cells_ge3 = mean(per_cell_n >= 3),
      mean_marker_data = mean(Matrix::rowMeans(sub_data)),
      stringsAsFactors = FALSE
    )
  }
}

set_summary <- do.call(rbind, set_summary_list)

write.csv(
  set_summary,
  file.path(TABDIR, "07_marker_set_evidence_v1.7.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Per-marker summary
# ------------------------------------------------------------------

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

marker_summary_list <- list()

for (cluster_id in clusters) {
  idx <- cl == cluster_id

  pct <- Matrix::rowMeans(
    counts[marker_order, idx, drop = FALSE] > 0
  ) * 100

  avg <- Matrix::rowMeans(
    data_mat[marker_order, idx, drop = FALSE]
  )

  marker_summary_list[[cluster_id]] <- data.frame(
    cluster_res0_8 = cluster_id,
    gene = marker_order,
    pct_detected = as.numeric(pct[marker_order]),
    mean_data = as.numeric(avg[marker_order]),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, marker_summary_list)
rownames(marker_summary) <- NULL

write.csv(
  marker_summary,
  file.path(TABDIR, "08_marker_summary_v1.7.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Screening ranks: broad identity + state evidence
# Screening only; NOT labels.
# ------------------------------------------------------------------

BROAD_SETS <- c(
  "Myeloid_core",
  "Dendritic",
  "Cycling",
  "T_NK_contamination",
  "B_Plasma_contamination",
  "Mesenchymal_contamination",
  "Hepatocyte_strict_contamination"
)

STATE_SETS <- c(
  "Resident_Kupffer_like",
  "Monocyte_inflammatory",
  "SPP1_TREM2_LAM",
  "ECM_inflammatory",
  "Repair_resolution",
  "Lipid_handling",
  "IFN_response"
)

screen_list <- list()

for (cluster_id in clusters) {
  broad <- set_summary[
    set_summary$cluster_res0_8 == cluster_id &
      set_summary$marker_set %in% BROAD_SETS,
    c("marker_set", "frac_cells_ge2"),
    drop = FALSE
  ]
  broad <- broad[order(broad$frac_cells_ge2, decreasing = TRUE), , drop = FALSE]

  state <- set_summary[
    set_summary$cluster_res0_8 == cluster_id &
      set_summary$marker_set %in% STATE_SETS,
    c("marker_set", "frac_cells_ge2"),
    drop = FALSE
  ]
  state <- state[order(state$frac_cells_ge2, decreasing = TRUE), , drop = FALSE]

  screen_list[[cluster_id]] <- data.frame(
    cluster_res0_8 = cluster_id,
    n_cells = sum(cl == cluster_id),

    broad_rank1 = broad$marker_set[1],
    broad_rank1_frac_ge2 = broad$frac_cells_ge2[1],
    broad_rank2 = broad$marker_set[2],
    broad_rank2_frac_ge2 = broad$frac_cells_ge2[2],

    state_rank1 = state$marker_set[1],
    state_rank1_frac_ge2 = state$frac_cells_ge2[1],
    state_rank2 = state$marker_set[2],
    state_rank2_frac_ge2 = state$frac_cells_ge2[2],

    stringsAsFactors = FALSE
  )
}

screen <- do.call(rbind, screen_list)
rownames(screen) <- NULL

screen <- merge(
  screen,
  balance,
  by = c("cluster_res0_8", "n_cells"),
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  screen,
  file.path(TABDIR, "09_preliminary_identity_state_screen_v1.7.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umap.myeloid.rpca")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  cluster = factor(cl, levels = clusters),
  sample = sample_id,
  NAS = factor(nas, levels = c(0, 2, 4, 6)),
  stringsAsFactors = FALSE
)

centers <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cluster,
  data = plot_df,
  FUN = median
)

p_umap <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = cluster)
) +
  geom_point(size = 0.45, alpha = 0.88) +
  geom_text(
    data = centers,
    aes(label = cluster),
    color = "black",
    fontface = "bold",
    size = 4,
    show.legend = FALSE
  ) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 myeloid RPCA | working res0.8",
    "Preliminary state audit; no biological labels frozen"
  )

ggsave(
  file.path(FIGDIR, "01_res08_numeric_UMAP_v1.7.2.pdf"),
  p_umap,
  width = 9,
  height = 7.5
)

# State evidence heatmap
heat_df <- set_summary[
  set_summary$marker_set %in% c(STATE_SETS, "Dendritic", "Cycling"),
  , drop = FALSE
]

heat_df$cluster_res0_8 <- factor(
  heat_df$cluster_res0_8,
  levels = clusters
)
heat_df$marker_set <- factor(
  heat_df$marker_set,
  levels = c(STATE_SETS, "Dendritic", "Cycling")
)

p_heat <- ggplot(
  heat_df,
  aes(
    x = cluster_res0_8,
    y = marker_set,
    fill = frac_cells_ge2
  )
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
    axis.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 myeloid res0.8 | state evidence",
    "Screening only; no state labels assigned"
  )

ggsave(
  file.path(FIGDIR, "02_res08_state_evidence_heatmap_v1.7.2.pdf"),
  p_heat,
  width = max(11, 0.58 * length(clusters) + 4),
  height = 7
)

# Marker DotPlot
dot_df <- marker_summary
dot_df$cluster_res0_8 <- factor(
  dot_df$cluster_res0_8,
  levels = clusters
)
dot_df$gene <- factor(
  dot_df$gene,
  levels = rev(marker_order)
)

upper <- as.numeric(
  quantile(dot_df$mean_data, 0.98, na.rm = TRUE)
)
if (!is.finite(upper) || upper <= 0) upper <- 1

p_dot <- ggplot(
  dot_df,
  aes(
    x = cluster_res0_8,
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
    axis.text.y = element_text(size = 7.3)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 myeloid res0.8 | marker audit",
    "No macrophage/monocyte state labels frozen"
  )

ggsave(
  file.path(FIGDIR, "03_res08_state_marker_dotplot_v1.7.2.pdf"),
  p_dot,
  width = max(12, 0.62 * length(clusters) + 5),
  height = max(14, 0.21 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------------

flags <- screen[
  screen$n_samples_present < 8 |
    screen$dominant_sample_fraction > 0.5 |
    screen$broad_rank1 != "Myeloid_core",
  , drop = FALSE
]

summary_lines <- c(
  paste0("Human MASH LT1-LT8 myeloid res0.8 state audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no reclustering, no relabeling, no RDS save.",
  "",
  paste0("Cells: ", ncol(obj)),
  paste0("RNA features: ", nrow(obj[["RNA"]])),
  paste0("Working resolution: 0.8"),
  paste0("Clusters: ", length(clusters)),
  "",
  "PRELIMINARY SCREEN:",
  paste(
    capture.output(print(screen, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "QC / BROAD-IDENTITY FLAGS:",
  if (nrow(flags)) {
    paste(
      capture.output(print(flags, row.names = FALSE)),
      collapse = "\n"
    )
  } else {
    "<none>"
  },
  "",
  "INTERPRETATION RULES:",
  "- broad_rank/state_rank are screening ranks only, not frozen labels.",
  "- Dendritic-like or non-myeloid clusters must be removed from the final MΦ/Monocyte object.",
  "- Cycling may be retained as a separate cycling-myeloid state only if myeloid evidence is coherent.",
  "- State labels should require coherent multi-marker support, not a single marker.",
  "- NAS composition is descriptive here; no inferential disease-stage test is performed.",
  "",
  "NEXT DECISION GATE:",
  "Assign preliminary cluster identities at res0.8 and define the final clean Human MΦ/Monocyte set."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "10_state_audit_summary_v1.7.2.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.7.2 COMPLETE")
message("Working resolution: 0.8")
message("Clusters: ", length(clusters))
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
