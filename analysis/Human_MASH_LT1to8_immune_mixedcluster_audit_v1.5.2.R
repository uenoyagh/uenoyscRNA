#!/usr/bin/env Rscript

# Human_MASH_LT1to8_immune_mixedcluster_audit_v1.5.2.R
#
# PURPOSE
# -------
# Resolve the two res0.6 clusters that materially affect construction of a
# clean human macrophage/monocyte object:
#
#   cluster 11: strong macrophage program + unusually high strict hepatocyte signal
#   cluster 19: macrophage program + strong mesenchymal signal
#
# Controls:
#   clusters 2,3,9  = macrophage-rich controls
#   cluster 18      = inflammatory-monocyte control
#   cluster 15      = mesenchymal-contamination control
#
# This audit asks whether mixed signals occur WITHIN the same cells and whether
# those cells have elevated nCount_RNA/nFeature_RNA or pre-existing doublet
# metadata, which would support doublet/multiplet contamination.
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

VERSION <- "v1.5.2"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5/",
  "objects/Human_MASH_LT1to8_ImmuneCandidate_RPCA_v1.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Immune_MixedClusterAudit_", VERSION
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
message("Human MASH LT1-LT8 immune mixed-cluster audit ", VERSION)
message("============================================================")
message("READ-ONLY: no clustering, no relabeling, no RDS save.")

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "immune_RPCA_res0_6"
TARGET_CLUSTERS <- c("11", "19")
CONTROL_CLUSTERS <- c("2", "3", "9", "18", "15")
AUDIT_CLUSTERS <- c(TARGET_CLUSTERS, CONTROL_CLUSTERS)

required_meta <- c(
  CLUSTER_COL,
  "sample_v15",
  "source_broad_label_v15"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")
if (!"umap.immune.rpca" %in% Reductions(obj)) {
  stop("umap.immune.rpca reduction absent.")
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
md <- obj[[]]

cl <- as.character(obj[[CLUSTER_COL, drop = TRUE]])
sample_id <- as.character(obj$sample_v15)
source_label <- as.character(obj$source_broad_label_v15)
features <- rownames(obj[["RNA"]])

# ------------------------------------------------------------------
# Marker sets
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
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "FCGR3A"
  ),
  SPP1_TREM2_state = c(
    "SPP1", "TREM2", "GPNMB", "LGALS3", "LPL", "FABP5", "APOE"
  ),
  Hepatocyte_strict = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
  ),
  Hepatocyte_ambient_prone = c(
    "ALB", "TTR", "APOA1"
  ),
  Mesenchymal_strict = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  )
)

marker_presence <- do.call(
  rbind,
  lapply(names(marker_sets), function(nm) {
    data.frame(
      marker_set = nm,
      gene = marker_sets[[nm]],
      present = marker_sets[[nm]] %in% features,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_presence,
  file.path(TABDIR, "01_marker_presence_v1.5.2.csv"),
  row.names = FALSE
)

present <- function(x) intersect(x, features)

# ------------------------------------------------------------------
# Cell-level detected-marker counts for all audit clusters
# ------------------------------------------------------------------

audit_idx <- which(cl %in% AUDIT_CLUSTERS)
audit_cells <- colnames(obj)[audit_idx]

cell_ev <- data.frame(
  cell_id = audit_cells,
  cluster = cl[audit_idx],
  sample = sample_id[audit_idx],
  source_label = source_label[audit_idx],
  stringsAsFactors = FALSE,
  row.names = audit_cells
)

for (nm in names(marker_sets)) {
  genes <- present(marker_sets[[nm]])
  if (!length(genes)) {
    cell_ev[[paste0(nm, "_detected_n")]] <- NA_integer_
  } else {
    cell_ev[[paste0(nm, "_detected_n")]] <- as.integer(
      Matrix::colSums(counts[genes, audit_cells, drop = FALSE] > 0)
    )
  }
}

# Standard QC columns.
for (nm in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
  if (nm %in% colnames(md)) {
    cell_ev[[nm]] <- md[audit_cells, nm]
  }
}

write.csv(
  cell_ev,
  file.path(TABDIR, "02_audit_cluster_cell_level_evidence_v1.5.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Mixed-lineage co-detection summary
# ------------------------------------------------------------------

mix_summary_list <- list()

for (cluster_id in AUDIT_CLUSTERS) {
  z <- cell_ev[cell_ev$cluster == cluster_id, , drop = FALSE]

  my2 <- z$Myeloid_core_detected_n >= 2
  rm2 <- z$Resident_macrophage_detected_n >= 2
  mono2 <- z$Monocyte_inflammatory_detected_n >= 2
  hep1 <- z$Hepatocyte_strict_detected_n >= 1
  hep2 <- z$Hepatocyte_strict_detected_n >= 2
  mes1 <- z$Mesenchymal_strict_detected_n >= 1
  mes2 <- z$Mesenchymal_strict_detected_n >= 2

  mix_summary_list[[cluster_id]] <- data.frame(
    cluster = cluster_id,
    n_cells = nrow(z),

    myeloid_ge2 = mean(my2),
    resident_macrophage_ge2 = mean(rm2),
    monocyte_inflammatory_ge2 = mean(mono2),

    hepatocyte_strict_ge1 = mean(hep1),
    hepatocyte_strict_ge2 = mean(hep2),
    mesenchymal_strict_ge1 = mean(mes1),
    mesenchymal_strict_ge2 = mean(mes2),

    myeloid_ge2_AND_hepatocyte_ge1 = mean(my2 & hep1),
    myeloid_ge2_AND_hepatocyte_ge2 = mean(my2 & hep2),

    resident_ge2_AND_hepatocyte_ge1 = mean(rm2 & hep1),
    resident_ge2_AND_hepatocyte_ge2 = mean(rm2 & hep2),

    myeloid_ge2_AND_mesenchymal_ge1 = mean(my2 & mes1),
    myeloid_ge2_AND_mesenchymal_ge2 = mean(my2 & mes2),

    resident_ge2_AND_mesenchymal_ge1 = mean(rm2 & mes1),
    resident_ge2_AND_mesenchymal_ge2 = mean(rm2 & mes2),

    stringsAsFactors = FALSE
  )
}

mix_summary <- do.call(rbind, mix_summary_list)
rownames(mix_summary) <- NULL

write.csv(
  mix_summary,
  file.path(TABDIR, "03_mixed_lineage_codetection_summary_v1.5.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# QC burden summary: nCount/nFeature/percent.mt if available
# ------------------------------------------------------------------

qc_cols <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  colnames(cell_ev)
)

qc_summary <- data.frame()

if (length(qc_cols)) {
  qlist <- list()

  for (cluster_id in AUDIT_CLUSTERS) {
    z <- cell_ev[cell_ev$cluster == cluster_id, , drop = FALSE]

    row <- data.frame(
      cluster = cluster_id,
      n_cells = nrow(z),
      stringsAsFactors = FALSE
    )

    for (nm in qc_cols) {
      x <- as.numeric(z[[nm]])
      row[[paste0(nm, "_median")]] <- median(x, na.rm = TRUE)
      row[[paste0(nm, "_Q25")]] <- quantile(x, 0.25, na.rm = TRUE)
      row[[paste0(nm, "_Q75")]] <- quantile(x, 0.75, na.rm = TRUE)
      row[[paste0(nm, "_Q90")]] <- quantile(x, 0.90, na.rm = TRUE)
    }

    qlist[[cluster_id]] <- row
  }

  qc_summary <- do.call(rbind, qlist)
  rownames(qc_summary) <- NULL
}

write.csv(
  qc_summary,
  file.path(TABDIR, "04_QC_burden_by_cluster_v1.5.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Audit any pre-existing doublet/singlet-related metadata.
# ------------------------------------------------------------------

doublet_pattern <- paste(
  c(
    "doublet", "singlet", "pann", "scdbl", "scrub",
    "dbl", "multiplet"
  ),
  collapse = "|"
)

doublet_cols <- grep(
  doublet_pattern,
  colnames(md),
  ignore.case = TRUE,
  value = TRUE
)

writeLines(
  doublet_cols,
  file.path(TABDIR, "05_doublet_related_metadata_columns_v1.5.2.txt")
)

doublet_summary <- data.frame()

if (length(doublet_cols)) {
  dlist <- list()

  for (nm in doublet_cols) {
    x <- md[audit_cells, nm]

    if (is.numeric(x) || is.integer(x)) {
      tmp <- do.call(
        rbind,
        lapply(AUDIT_CLUSTERS, function(cluster_id) {
          idx <- cell_ev$cluster == cluster_id
          xx <- as.numeric(x[idx])

          data.frame(
            metadata_column = nm,
            cluster = cluster_id,
            type = "numeric",
            level = NA_character_,
            n_cells = sum(idx),
            value = mean(xx, na.rm = TRUE),
            median = median(xx, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
      )
    } else {
      tmp <- as.data.frame(
        table(
          cluster = cell_ev$cluster,
          level = as.character(x),
          useNA = "ifany"
        ),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[tmp$Freq > 0, , drop = FALSE]
      tmp$metadata_column <- nm
      tmp$type <- "categorical"
      tmp$n_cells <- tmp$Freq
      tmp$value <- NA_real_
      tmp$median <- NA_real_
      tmp <- tmp[
        , c(
          "metadata_column", "cluster", "type", "level",
          "n_cells", "value", "median"
        )
      ]
    }

    dlist[[nm]] <- tmp
  }

  doublet_summary <- do.call(rbind, dlist)
  rownames(doublet_summary) <- NULL
}

write.csv(
  doublet_summary,
  file.path(TABDIR, "06_doublet_metadata_summary_v1.5.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Sample and transferred-source composition
# ------------------------------------------------------------------

sample_tab <- as.data.frame(
  table(
    cluster = cl[cl %in% AUDIT_CLUSTERS],
    sample = sample_id[cl %in% AUDIT_CLUSTERS]
  ),
  stringsAsFactors = FALSE
)
sample_tab <- sample_tab[sample_tab$Freq > 0, , drop = FALSE]
colnames(sample_tab)[3] <- "n_cells"

write.csv(
  sample_tab,
  file.path(TABDIR, "07_audit_clusters_by_sample_v1.5.2.csv"),
  row.names = FALSE
)

source_tab <- as.data.frame(
  table(
    cluster = cl[cl %in% AUDIT_CLUSTERS],
    source_label = source_label[cl %in% AUDIT_CLUSTERS]
  ),
  stringsAsFactors = FALSE
)
source_tab <- source_tab[source_tab$Freq > 0, , drop = FALSE]
colnames(source_tab)[3] <- "n_cells"

source_tot <- aggregate(n_cells ~ cluster, source_tab, sum)
colnames(source_tot)[2] <- "cluster_total"
source_tab <- merge(source_tab, source_tot, by = "cluster", all.x = TRUE)
source_tab$fraction <- source_tab$n_cells / source_tab$cluster_total

write.csv(
  source_tab,
  file.path(TABDIR, "08_audit_clusters_source_label_composition_v1.5.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

# UMAP highlight of audit clusters.
um <- Embeddings(obj, reduction = "umap.immune.rpca")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  cluster = cl,
  highlight = ifelse(
    cl %in% TARGET_CLUSTERS,
    paste0("target_", cl),
    ifelse(cl %in% CONTROL_CLUSTERS, paste0("control_", cl), "other")
  ),
  stringsAsFactors = FALSE
)

p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
  geom_point(
    data = plot_df[plot_df$highlight == "other", , drop = FALSE],
    color = "grey85",
    size = 0.20
  ) +
  geom_point(
    data = plot_df[plot_df$highlight != "other", , drop = FALSE],
    aes(color = highlight),
    size = 0.45,
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
    "Human MASH LT1-LT8 immune RPCA | mixed-cluster audit",
    "Targets: cluster 11 and 19; controls shown for comparison"
  )

ggsave(
  file.path(FIGDIR, "01_target_control_UMAP_v1.5.2.pdf"),
  p_umap,
  width = 9,
  height = 7.5
)

# Cell-level marker-count scatter.
scatter_df <- cell_ev

p_hep <- ggplot(
  scatter_df,
  aes(
    x = Resident_macrophage_detected_n,
    y = Hepatocyte_strict_detected_n
  )
) +
  geom_jitter(width = 0.15, height = 0.08, alpha = 0.28, size = 0.55) +
  facet_wrap(~ cluster, ncol = 4) +
  theme_classic(base_size = 10) +
  labs(
    x = "Resident-macrophage markers detected / cell",
    y = "Strict hepatocyte markers detected / cell",
    title = "Within-cell macrophage + hepatocyte co-detection"
  )

ggsave(
  file.path(FIGDIR, "02_macrophage_hepatocyte_codetection_v1.5.2.pdf"),
  p_hep,
  width = 11,
  height = 7.5
)

p_mes <- ggplot(
  scatter_df,
  aes(
    x = Resident_macrophage_detected_n,
    y = Mesenchymal_strict_detected_n
  )
) +
  geom_jitter(width = 0.15, height = 0.08, alpha = 0.28, size = 0.55) +
  facet_wrap(~ cluster, ncol = 4) +
  theme_classic(base_size = 10) +
  labs(
    x = "Resident-macrophage markers detected / cell",
    y = "Strict mesenchymal markers detected / cell",
    title = "Within-cell macrophage + mesenchymal co-detection"
  )

ggsave(
  file.path(FIGDIR, "03_macrophage_mesenchymal_codetection_v1.5.2.pdf"),
  p_mes,
  width = 11,
  height = 7.5
)

# QC distributions if available.
if ("nFeature_RNA" %in% colnames(cell_ev)) {
  p_nf <- ggplot(
    cell_ev,
    aes(x = factor(cluster), y = nFeature_RNA)
  ) +
    geom_boxplot(outlier.size = 0.25) +
    theme_classic(base_size = 10) +
    labs(
      x = "res0.6 cluster",
      y = "nFeature_RNA",
      title = "nFeature_RNA burden in mixed-cluster audit"
    )

  ggsave(
    file.path(FIGDIR, "04_nFeature_RNA_by_cluster_v1.5.2.pdf"),
    p_nf,
    width = 8,
    height = 5.5
  )
}

if ("nCount_RNA" %in% colnames(cell_ev)) {
  p_nc <- ggplot(
    cell_ev,
    aes(x = factor(cluster), y = nCount_RNA)
  ) +
    geom_boxplot(outlier.size = 0.25) +
    theme_classic(base_size = 10) +
    labs(
      x = "res0.6 cluster",
      y = "nCount_RNA",
      title = "nCount_RNA burden in mixed-cluster audit"
    )

  ggsave(
    file.path(FIGDIR, "05_nCount_RNA_by_cluster_v1.5.2.pdf"),
    p_nc,
    width = 8,
    height = 5.5
  )
}

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

summary_rows <- apply(
  mix_summary,
  1,
  function(r) {
    paste0(
      "cluster ", r[["cluster"]],
      " | n=", r[["n_cells"]],
      " | Mφ>=2=", fmt_pct(as.numeric(r[["resident_macrophage_ge2"]])),
      " | hep>=2=", fmt_pct(as.numeric(r[["hepatocyte_strict_ge2"]])),
      " | Mφ>=2 & hep>=2=",
      fmt_pct(as.numeric(r[["resident_ge2_AND_hepatocyte_ge2"]])),
      " | mes>=2=", fmt_pct(as.numeric(r[["mesenchymal_strict_ge2"]])),
      " | Mφ>=2 & mes>=2=",
      fmt_pct(as.numeric(r[["resident_ge2_AND_mesenchymal_ge2"]]))
    )
  }
)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 immune mixed-cluster audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "READ-ONLY: no reclustering, no relabeling, no RDS save.",
  "",
  "TARGETS:",
  "cluster 11 = macrophage-rich + strict-hepatocyte signal",
  "cluster 19 = macrophage/mesenchymal mixed signal",
  "",
  "CONTROLS:",
  "clusters 2,3,9 = macrophage-rich",
  "cluster 18 = inflammatory-monocyte",
  "cluster 15 = mesenchymal contamination",
  "",
  "WITHIN-CELL CO-DETECTION:",
  summary_rows,
  "",
  "DOUBLEt-RELATED METADATA COLUMNS FOUND:",
  if (length(doublet_cols)) paste(doublet_cols, collapse = " | ") else "<none>",
  "",
  "DECISION RULE:",
  "- High within-cell resident-macrophage + strict-hepatocyte/mesenchymal co-detection,",
  "  especially with elevated nCount/nFeature or doublet metadata, supports exclusion",
  "  as mixed/doublet contamination.",
  "- If strict non-myeloid signal is low at the single-cell level despite high ambient",
  "  liver transcripts, retain as macrophage.",
  "",
  "NEXT STEP:",
  "Use this audit to decide cluster 11 and 19 inclusion/exclusion before creating",
  "the clean human macrophage/monocyte object."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "09_mixedcluster_audit_summary_v1.5.2.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.5.2 COMPLETE")
message("Targets audited: 11, 19")
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
