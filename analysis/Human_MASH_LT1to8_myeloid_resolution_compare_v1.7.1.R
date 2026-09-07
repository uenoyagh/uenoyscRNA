#!/usr/bin/env Rscript

# Human_MASH_LT1to8_myeloid_resolution_compare_v1.7.1.R
#
# PURPOSE
# -------
# Targeted READ-ONLY comparison of myeloid RPCA resolutions 0.6 vs 0.8
# in the v1.7 object before choosing/fixing a working resolution.
#
# WHY 0.6 vs 0.8
# --------------
# v1.7 produced:
#   res0.6 = 16 clusters
#   res0.8 = 17 clusters
#   res1.0 = 19 clusters
#   res1.2 = 21 clusters
# Sample balance remains broadly good at res0.6 and res0.8.
# This script determines whether the single extra split at res0.8 carries
# coherent biological marker information or merely fragments a res0.6 state.
#
# THIS SCRIPT DOES:
#   - exact res0.6 <-> res0.8 cluster crosswalk
#   - parent-to-child purity / split metrics
#   - sample balance for split clusters
#   - marker summaries for BOTH resolutions
#   - targeted marker comparison for clusters actually split at res0.8
#   - UMAPs highlighting the split
#
# THIS SCRIPT DOES NOT:
#   - recluster
#   - recompute UMAP
#   - modify labels
#   - save/overwrite an RDS
#   - perform DE / pseudobulk / pathway analysis
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.7.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_RPCA_v1.7/",
  "objects/Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_ResolutionCompare_", VERSION
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
message("Human MASH LT1-LT8 myeloid resolution comparison ", VERSION)
message("============================================================")
message("READ-ONLY: comparing res0.6 vs res0.8.")

obj <- readRDS(INPUT_RDS)

COL06 <- "myeloid_RPCA_res0_6"
COL08 <- "myeloid_RPCA_res0_8"

required_meta <- c(COL06, COL08, "sample_v17", "NAS_v17")
miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) {
  stop("RNA counts/data layers are required.")
}
if (!"umap.myeloid.rpca" %in% Reductions(obj)) {
  stop("umap.myeloid.rpca reduction absent.")
}

cl06 <- as.character(obj[[COL06, drop = TRUE]])
cl08 <- as.character(obj[[COL08, drop = TRUE]])
sample_id <- as.character(obj$sample_v17)

# ------------------------------------------------------------------
# 01 Exact cluster crosswalk
# ------------------------------------------------------------------

cross <- as.data.frame(
  table(res0_6 = cl06, res0_8 = cl08),
  stringsAsFactors = FALSE
)
cross <- cross[cross$Freq > 0, , drop = FALSE]
colnames(cross)[3] <- "n_cells"

parent_total <- aggregate(n_cells ~ res0_6, cross, sum)
colnames(parent_total)[2] <- "parent_n"

child_total <- aggregate(n_cells ~ res0_8, cross, sum)
colnames(child_total)[2] <- "child_n"

cross <- merge(cross, parent_total, by = "res0_6", all.x = TRUE)
cross <- merge(cross, child_total, by = "res0_8", all.x = TRUE)
cross$fraction_of_res0_6 <- cross$n_cells / cross$parent_n
cross$fraction_of_res0_8 <- cross$n_cells / cross$child_n

cross <- cross[
  order(as.numeric(cross$res0_6), -cross$fraction_of_res0_6),
  , drop = FALSE
]

write.csv(
  cross,
  file.path(TABDIR, "01_res06_res08_cluster_crosswalk_v1.7.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 02 Parent split summary
# ------------------------------------------------------------------

split_summary_list <- list()

for (p in sort(unique(cl06))) {
  z <- cross[cross$res0_6 == p, , drop = FALSE]
  z <- z[order(z$n_cells, decreasing = TRUE), , drop = FALSE]

  split_summary_list[[p]] <- data.frame(
    res0_6_cluster = p,
    n_cells = unique(z$parent_n)[1],
    n_res0_8_children = nrow(z),
    dominant_res0_8_child = z$res0_8[1],
    dominant_child_fraction = z$fraction_of_res0_6[1],
    second_child_fraction = if (nrow(z) >= 2) z$fraction_of_res0_6[2] else 0,
    effective_split = nrow(z) >= 2 && z$fraction_of_res0_6[2] >= 0.10,
    stringsAsFactors = FALSE
  )
}

split_summary <- do.call(rbind, split_summary_list)
rownames(split_summary) <- NULL

write.csv(
  split_summary,
  file.path(TABDIR, "02_res06_parent_split_summary_v1.7.1.csv"),
  row.names = FALSE
)

split_parents <- split_summary$res0_6_cluster[split_summary$effective_split]

# ------------------------------------------------------------------
# 03 Sample balance for res0.8 children of effectively split parents
# ------------------------------------------------------------------

split_child_pairs <- cross[
  cross$res0_6 %in% split_parents &
    cross$fraction_of_res0_6 >= 0.10,
  c("res0_6", "res0_8"),
  drop = FALSE
]

split_balance_list <- list()

if (nrow(split_child_pairs)) {
  for (i in seq_len(nrow(split_child_pairs))) {
    p <- split_child_pairs$res0_6[i]
    c08 <- split_child_pairs$res0_8[i]

    idx <- cl06 == p & cl08 == c08
    tt <- table(sample_id[idx])
    frac <- as.numeric(tt) / sum(tt)
    h <- -sum(frac * log(frac))
    h_norm <- if (length(tt) > 1) h / log(8) else 0

    split_balance_list[[length(split_balance_list) + 1L]] <- data.frame(
      res0_6_parent = p,
      res0_8_child = c08,
      n_cells = sum(idx),
      n_samples_present = length(tt),
      dominant_sample = names(tt)[which.max(tt)],
      dominant_sample_fraction = max(frac),
      normalized_sample_entropy = h_norm,
      stringsAsFactors = FALSE
    )
  }
}

split_balance <- if (length(split_balance_list)) {
  do.call(rbind, split_balance_list)
} else {
  data.frame()
}

write.csv(
  split_balance,
  file.path(TABDIR, "03_effective_split_child_sample_balance_v1.7.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Marker panel
# ------------------------------------------------------------------

marker_sets <- list(
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Resident_macrophage = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2",
    "TIMD4", "MARCO", "VSIG4", "CD5L"
  ),
  Monocyte_inflammatory = c(
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "FCGR3A",
    "IL1B", "VCAN"
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
  Interferon_response = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44", "CXCL10"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  Hepatocyte_strict_contamination = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
  ),
  Mesenchymal_contamination = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1"
  )
)

features <- rownames(obj[["RNA"]])
marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

summarize_resolution <- function(cluster_vec, resolution_label) {
  out <- list()

  for (cluster_id in sort(unique(cluster_vec))) {
    idx <- cluster_vec == cluster_id

    pct <- Matrix::rowMeans(
      counts[marker_order, idx, drop = FALSE] > 0
    ) * 100

    avg <- Matrix::rowMeans(
      data_mat[marker_order, idx, drop = FALSE]
    )

    out[[cluster_id]] <- data.frame(
      resolution = resolution_label,
      cluster = cluster_id,
      gene = marker_order,
      pct_detected = as.numeric(pct[marker_order]),
      mean_data = as.numeric(avg[marker_order]),
      n_cells = sum(idx),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

marker06 <- summarize_resolution(cl06, "0.6")
marker08 <- summarize_resolution(cl08, "0.8")

marker_all <- rbind(marker06, marker08)
rownames(marker_all) <- NULL

write.csv(
  marker_all,
  file.path(TABDIR, "04_marker_summary_res06_res08_v1.7.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 05 Marker-set evidence for effective splits only
# ------------------------------------------------------------------

split_evidence_list <- list()

if (length(split_parents)) {
  for (p in split_parents) {
    children <- unique(
      cross$res0_8[
        cross$res0_6 == p &
          cross$fraction_of_res0_6 >= 0.10
      ]
    )

    # Parent
    idx_parent <- cl06 == p

    for (set_name in names(marker_sets)) {
      genes <- intersect(marker_sets[[set_name]], features)
      if (!length(genes)) next

      det <- Matrix::colSums(counts[genes, idx_parent, drop = FALSE] > 0)

      split_evidence_list[[length(split_evidence_list) + 1L]] <- data.frame(
        res0_6_parent = p,
        level = "parent_res0.6",
        cluster = p,
        marker_set = set_name,
        n_cells = sum(idx_parent),
        frac_cells_ge2 = mean(det >= 2),
        frac_cells_ge3 = mean(det >= 3),
        stringsAsFactors = FALSE
      )
    }

    # Children
    for (c08 in children) {
      idx_child <- cl06 == p & cl08 == c08

      for (set_name in names(marker_sets)) {
        genes <- intersect(marker_sets[[set_name]], features)
        if (!length(genes)) next

        det <- Matrix::colSums(counts[genes, idx_child, drop = FALSE] > 0)

        split_evidence_list[[length(split_evidence_list) + 1L]] <- data.frame(
          res0_6_parent = p,
          level = "child_res0.8",
          cluster = c08,
          marker_set = set_name,
          n_cells = sum(idx_child),
          frac_cells_ge2 = mean(det >= 2),
          frac_cells_ge3 = mean(det >= 3),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

split_evidence <- if (length(split_evidence_list)) {
  do.call(rbind, split_evidence_list)
} else {
  data.frame()
}

write.csv(
  split_evidence,
  file.path(TABDIR, "05_effective_split_marker_set_evidence_v1.7.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umap.myeloid.rpca")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  res0_6 = cl06,
  res0_8 = cl08,
  sample = sample_id,
  stringsAsFactors = FALSE,
  row.names = colnames(obj)
)

# Highlight only effectively split res0.6 parent(s).
if (length(split_parents)) {
  split_plot <- plot_df
  split_plot$highlight <- ifelse(
    split_plot$res0_6 %in% split_parents,
    paste0("res0.6_", split_plot$res0_6, " -> res0.8_", split_plot$res0_8),
    "other"
  )

  p_split <- ggplot(split_plot, aes(UMAP_1, UMAP_2)) +
    geom_point(
      data = split_plot[split_plot$highlight == "other", , drop = FALSE],
      color = "grey88",
      size = 0.22
    ) +
    geom_point(
      data = split_plot[split_plot$highlight != "other", , drop = FALSE],
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
      "Human MASH LT1-LT8 myeloid RPCA | res0.6 to res0.8 effective split",
      "Only parent clusters with >=10% secondary child are highlighted"
    )

  ggsave(
    file.path(FIGDIR, "01_effective_split_UMAP_v1.7.1.pdf"),
    p_split,
    width = 10,
    height = 7.5
  )
}

# DotPlots for both resolutions with identical marker order.
plot_dot <- function(df, title, outfile) {
  df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster)))
  df$gene <- factor(df$gene, levels = rev(marker_order))

  upper <- as.numeric(quantile(df$mean_data, 0.98, na.rm = TRUE))
  if (!is.finite(upper) || upper <= 0) upper <- 1

  p <- ggplot(
    df,
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
      axis.text.y = element_text(size = 7.3)
    ) +
    ggtitle(title)

  ggsave(
    outfile,
    p,
    width = max(11, 0.62 * length(unique(df$cluster)) + 5),
    height = max(14, 0.21 * length(marker_order) + 3)
  )
}

plot_dot(
  marker06,
  "Human MASH LT1-LT8 myeloid | res0.6 marker review",
  file.path(FIGDIR, "02_marker_dotplot_res0.6_v1.7.1.pdf")
)

plot_dot(
  marker08,
  "Human MASH LT1-LT8 myeloid | res0.8 marker review",
  file.path(FIGDIR, "03_marker_dotplot_res0.8_v1.7.1.pdf")
)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

effective_text <- if (length(split_parents)) {
  paste(split_parents, collapse = " | ")
} else {
  "<none>"
}

summary_lines <- c(
  paste0("Human MASH LT1-LT8 myeloid resolution comparison ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no reclustering, no relabeling, no RDS save.",
  "",
  paste0("res0.6 clusters: ", length(unique(cl06))),
  paste0("res0.8 clusters: ", length(unique(cl08))),
  paste0("Effective res0.6 parents split at res0.8: ", effective_text),
  "",
  "PARENT SPLIT SUMMARY:",
  paste(
    capture.output(print(split_summary, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "DECISION RULE:",
  "- Prefer res0.6 if the extra res0.8 split lacks a coherent marker-state difference",
  "  or creates a sample-skewed fragment.",
  "- Prefer res0.8 if the split is sample-balanced and separates coherent, biologically",
  "  interpretable macrophage/monocyte programs.",
  "- Do not choose resolution by cluster count alone.",
  "",
  "NEXT STEP:",
  "Choose/freeze the working myeloid resolution, then perform preliminary state annotation."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "06_resolution_compare_summary_v1.7.1.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.7.1 COMPLETE")
message("Effective split parents: ", effective_text)
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
