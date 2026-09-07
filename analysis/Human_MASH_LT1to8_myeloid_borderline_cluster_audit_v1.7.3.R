#!/usr/bin/env Rscript

# Human_MASH_LT1to8_myeloid_borderline_cluster_audit_v1.7.3.R
#
# PURPOSE
# -------
# Resolve the remaining borderline res0.8 clusters before defining the final
# Human MΦ / Monocyte set.
#
# TARGETS:
#   12 = DC-like (expected exclusion control)
#   13 = low-myeloid / monocyte-like candidate
#   14 = small lipid-handling + inflammatory myeloid candidate
#   15 = myeloid/resident cluster with T/NK contamination signal
#   16 = tiny, sample-limited resident/repair cluster
#
# CONTROLS:
#   2 = resident macrophage control
#   6 = SPP1/TREM2-LAM control
#   7 = inflammatory monocyte control
#
# READ-ONLY:
#   - NO reclustering
#   - NO UMAP recalculation
#   - NO relabeling
#   - NO RDS save
#   - NO DE / pseudobulk / pathway analysis
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.7.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_RPCA_v1.7/",
  "objects/Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_BorderlineAudit_", VERSION
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
message("Human MASH LT1-LT8 myeloid borderline audit ", VERSION)
message("============================================================")
message("READ-ONLY: no reclustering, no relabeling, no RDS save.")

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "myeloid_RPCA_res0_8"

TARGETS <- c("12", "13", "14", "15", "16")
CONTROLS <- c("2", "6", "7")
AUDIT_CLUSTERS <- c(TARGETS, CONTROLS)

required_meta <- c(
  CLUSTER_COL,
  "sample_v17",
  "NAS_v17",
  "parent_cluster_v17"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")
if (!"umap.myeloid.rpca" %in% Reductions(obj)) stop("UMAP reduction absent.")

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
md <- obj[[]]

cl <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
sample_id <- as.character(obj$sample_v17)
nas <- as.numeric(obj$NAS_v17)
parent_cluster <- as.character(obj$parent_cluster_v17)
features <- rownames(obj[["RNA"]])

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
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  T_cell = c(
    "CD3D", "CD3E", "TRAC", "CD247", "IL7R"
  ),
  NK_cell = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1"
  ),
  B_Plasma = c(
    "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"
  ),
  Mesenchymal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  ),
  Hepatocyte_strict = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
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
      present = genes %in% features,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_manifest,
  file.path(TABDIR, "01_marker_presence_v1.7.3.csv"),
  row.names = FALSE
)

present <- function(x) intersect(x, features)

# ------------------------------------------------------------------
# Cell-level evidence for targets + controls
# ------------------------------------------------------------------

audit_idx <- which(cl %in% AUDIT_CLUSTERS)
audit_cells <- colnames(obj)[audit_idx]

cell_ev <- data.frame(
  cell_id = audit_cells,
  cluster_res0_8 = cl[audit_idx],
  sample = sample_id[audit_idx],
  NAS = nas[audit_idx],
  parent_cluster = parent_cluster[audit_idx],
  stringsAsFactors = FALSE,
  row.names = audit_cells
)

for (set_name in names(marker_sets)) {
  genes <- present(marker_sets[[set_name]])

  if (!length(genes)) {
    cell_ev[[paste0(set_name, "_detected_n")]] <- NA_integer_
  } else {
    cell_ev[[paste0(set_name, "_detected_n")]] <- as.integer(
      Matrix::colSums(counts[genes, audit_cells, drop = FALSE] > 0)
    )
  }
}

for (nm in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
  if (nm %in% colnames(md)) {
    cell_ev[[nm]] <- md[audit_cells, nm]
  }
}

write.csv(
  cell_ev,
  file.path(TABDIR, "02_cell_level_lineage_evidence_v1.7.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Cluster-level evidence
# ------------------------------------------------------------------

cluster_summary_list <- list()

for (cluster_id in AUDIT_CLUSTERS) {
  z <- cell_ev[cell_ev$cluster_res0_8 == cluster_id, , drop = FALSE]

  row <- data.frame(
    cluster_res0_8 = cluster_id,
    n_cells = nrow(z),
    n_samples_present = length(unique(z$sample)),
    dominant_sample = names(which.max(table(z$sample))),
    dominant_sample_fraction = max(table(z$sample)) / nrow(z),
    stringsAsFactors = FALSE
  )

  for (set_name in names(marker_sets)) {
    x <- z[[paste0(set_name, "_detected_n")]]
    row[[paste0(set_name, "_ge1")]] <- mean(x >= 1, na.rm = TRUE)
    row[[paste0(set_name, "_ge2")]] <- mean(x >= 2, na.rm = TRUE)
  }

  if ("nCount_RNA" %in% colnames(z)) {
    row$nCount_RNA_median <- median(z$nCount_RNA, na.rm = TRUE)
    row$nCount_RNA_Q90 <- quantile(z$nCount_RNA, 0.90, na.rm = TRUE)
  }

  if ("nFeature_RNA" %in% colnames(z)) {
    row$nFeature_RNA_median <- median(z$nFeature_RNA, na.rm = TRUE)
    row$nFeature_RNA_Q90 <- quantile(z$nFeature_RNA, 0.90, na.rm = TRUE)
  }

  cluster_summary_list[[cluster_id]] <- row
}

cluster_summary <- do.call(rbind, cluster_summary_list)
rownames(cluster_summary) <- NULL

write.csv(
  cluster_summary,
  file.path(TABDIR, "03_borderline_cluster_summary_v1.7.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Parent provenance and NAS composition
# ------------------------------------------------------------------

parent_tab <- as.data.frame(
  table(
    cluster_res0_8 = cl[cl %in% AUDIT_CLUSTERS],
    parent_v1_5_cluster = parent_cluster[cl %in% AUDIT_CLUSTERS]
  ),
  stringsAsFactors = FALSE
)
parent_tab <- parent_tab[parent_tab$Freq > 0, , drop = FALSE]
colnames(parent_tab)[3] <- "n_cells"

write.csv(
  parent_tab,
  file.path(TABDIR, "04_parent_cluster_provenance_v1.7.3.csv"),
  row.names = FALSE
)

nas_tab <- as.data.frame(
  table(
    cluster_res0_8 = cl[cl %in% AUDIT_CLUSTERS],
    NAS = nas[cl %in% AUDIT_CLUSTERS]
  ),
  stringsAsFactors = FALSE
)
nas_tab <- nas_tab[nas_tab$Freq > 0, , drop = FALSE]
colnames(nas_tab)[3] <- "n_cells"

nas_total <- aggregate(n_cells ~ cluster_res0_8, nas_tab, sum)
colnames(nas_total)[2] <- "cluster_total"
nas_tab <- merge(nas_tab, nas_total, by = "cluster_res0_8")
nas_tab$fraction <- nas_tab$n_cells / nas_tab$cluster_total

write.csv(
  nas_tab,
  file.path(TABDIR, "05_NAS_composition_descriptive_v1.7.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Direct per-marker summary
# ------------------------------------------------------------------

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

marker_summary_list <- list()

for (cluster_id in AUDIT_CLUSTERS) {
  idx <- cl == cluster_id

  pct <- Matrix::rowMeans(counts[marker_order, idx, drop = FALSE] > 0) * 100

  marker_summary_list[[cluster_id]] <- data.frame(
    cluster_res0_8 = cluster_id,
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
  file.path(TABDIR, "06_marker_detection_summary_v1.7.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umap.myeloid.rpca")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  cluster = cl,
  stringsAsFactors = FALSE
)

plot_df$highlight <- ifelse(
  plot_df$cluster %in% TARGETS,
  paste0("target_", plot_df$cluster),
  ifelse(
    plot_df$cluster %in% CONTROLS,
    paste0("control_", plot_df$cluster),
    "other"
  )
)

p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
  geom_point(
    data = plot_df[plot_df$highlight == "other", , drop = FALSE],
    color = "grey88",
    size = 0.22
  ) +
  geom_point(
    data = plot_df[plot_df$highlight != "other", , drop = FALSE],
    aes(color = highlight),
    size = 0.55,
    alpha = 0.9
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
    "Human MASH LT1-LT8 myeloid res0.8 | borderline-cluster audit",
    "Targets 12-16; controls 2, 6, 7"
  )

ggsave(
  file.path(FIGDIR, "01_borderline_cluster_UMAP_v1.7.3.pdf"),
  p_umap,
  width = 10,
  height = 7.5
)

# Heatmap: >=2 markers per cell fraction
heat_rows <- list()

for (cluster_id in AUDIT_CLUSTERS) {
  z <- cluster_summary[
    cluster_summary$cluster_res0_8 == cluster_id,
    , drop = FALSE
  ]

  for (set_name in names(marker_sets)) {
    heat_rows[[length(heat_rows) + 1L]] <- data.frame(
      cluster_res0_8 = cluster_id,
      marker_set = set_name,
      frac_ge2 = z[[paste0(set_name, "_ge2")]],
      stringsAsFactors = FALSE
    )
  }
}

heat_df <- do.call(rbind, heat_rows)
heat_df$cluster_res0_8 <- factor(
  heat_df$cluster_res0_8,
  levels = AUDIT_CLUSTERS
)
heat_df$marker_set <- factor(
  heat_df$marker_set,
  levels = names(marker_sets)
)

p_heat <- ggplot(
  heat_df,
  aes(
    x = cluster_res0_8,
    y = marker_set,
    fill = frac_ge2
  )
) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * frac_ge2)),
    size = 2.7
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Cells with\n>=2 markers"
  ) +
  theme_classic(base_size = 10) +
  theme(axis.title = element_blank()) +
  ggtitle(
    "Human MASH LT1-LT8 myeloid res0.8 | borderline lineage evidence",
    "Targets 12-16 compared with resident/LAM/monocyte controls"
  )

ggsave(
  file.path(FIGDIR, "02_borderline_cluster_evidence_heatmap_v1.7.3.pdf"),
  p_heat,
  width = 10,
  height = 7
)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

summary_lines <- c(
  paste0("Human MASH LT1-LT8 myeloid borderline-cluster audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no reclustering, no relabeling, no RDS save.",
  "",
  "TARGET CLUSTERS:",
  "12 = DC-like candidate",
  "13 = low-myeloid / monocyte-like candidate",
  "14 = lipid-handling + inflammatory myeloid candidate",
  "15 = myeloid/resident with T/NK contamination signal",
  "16 = tiny sample-limited resident/repair candidate",
  "",
  "CONTROLS:",
  "2 = resident macrophage",
  "6 = SPP1/TREM2-LAM",
  "7 = inflammatory monocyte",
  "",
  "CLUSTER SUMMARY:",
  paste(
    capture.output(print(cluster_summary, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "DECISION LOGIC:",
  "- Exclude a cluster from the final MΦ/Monocyte set if non-myeloid lineage evidence",
  "  is coherent and exceeds the myeloid program.",
  "- Retain low-complexity monocyte clusters if monocyte markers are coherent and",
  "  non-myeloid contamination is low.",
  "- Tiny/sample-limited clusters should not be frozen as independent states unless",
  "  multi-marker support is strong and provenance is credible.",
  "",
  "NEXT STEP:",
  "Finalize the res0.8 Human MΦ/Monocyte cluster set and then write the frozen",
  "broad MΦ/Monocyte annotation object."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "07_borderline_audit_summary_v1.7.3.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.7.3 COMPLETE")
message("Targets: ", paste(TARGETS, collapse = " | "))
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
