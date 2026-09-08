# ============================================================
# Human MASH v1.3.2B-QC1
# Targeted LSEC-vs-Hepatocyte Audit for clusters 0 and 7
#
# PURPOSE
#   Resolve the identity of hepatobiliary RPCA clusters 0 and 7,
#   which show endothelial/LSEC-associated markers despite being
#   largely derived from the prior Hepatocyte compartment.
#
# IMPORTANT
#   - NO RPCA reintegration.
#   - NO reclustering.
#   - NO FindAllMarkers rerun.
#   - Reads the v1.3.2B pre-plot checkpoint.
#   - Uses validated canonical LSEC and Endothelial cells only as
#     RNA-expression reference populations.
#   - Canonical v1.3 and v1.3.2B objects are NEVER modified.
#   - This is a targeted QC/audit, not a final annotation step.
#
# TARGETS
#   cluster 0
#   cluster 7
#
# REFERENCE GROUPS
#   Other_Hepatocyte_labelled:
#       current REF12 Hepatocyte cells excluding clusters 0/7
#   Canonical_LSEC_reference:
#       canonical celltype_ref12_v13 == "LSEC"
#   Canonical_Endothelial_reference:
#       canonical celltype_ref12_v13 == "Endothelial"
#
# MARKER AXES
#   Hepatocyte:
#       ALB TTR APOA1 ASGR1 HNF4A CYP3A4 CPS1 GLUL
#   LSEC:
#       CLEC4G FCGR2B STAB2 STAB1 EMCN PLVAP RAMP2
#       KDR PECAM1 VWF EGFL7 ACKR3
#   Macrophage contamination:
#       PTPRC LST1 TYROBP FCER1G C1QC
#   Mesenchymal contamination:
#       COL1A1 COL1A2 DCN LUM
#
# OUTPUT
#   Group-level expression QC
#   Target-cluster score summaries
#   DotPlot-like marker validation
#   UMAP highlighting clusters 0/7
#   FeaturePlots on hepatobiliary UMAP
#   Signature-score Violin/boxplots
#   Review flags
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

options(stringsAsFactors = FALSE)

seed <- 20260908
set.seed(seed)

# ------------------------------------------------------------
# 0. Paths / configuration
# ------------------------------------------------------------

canonical_rds <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/",
  "Human_MASH_LT1to8_Canonical_v1.3.rds"
)

parent_outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B_hepatobiliary_refinement"
)

checkpoint_rds <- file.path(
  parent_outdir,
  "objects",
  "Human_MASH_Hepatobiliary_v1.3.2B_preplot_checkpoint.rds"
)

top20_csv <- file.path(
  parent_outdir,
  "15_hepatobiliary_cluster_markers_top20.csv"
)

outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B_QC1_LSEC_vs_Hepatocyte_audit"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

celltype_col <- "celltype_ref12_v13"
sample_col <- "canonical_sample_v13"
cluster_col <- "hepatobiliary_cluster_v132B"
reduction_name <- "umapHepatoRPCA_v132B"

expected_canonical_cells <- 81348L
expected_hepatobiliary_cells <- 47425L
expected_clusters <- as.character(0:16)

target_clusters <- c(
  "0",
  "7"
)

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

get_data_layer <- function(
    object,
    assay = "RNA"
) {

  tryCatch(
    GetAssayData(
      object,
      assay = assay,
      layer = "data"
    ),
    error = function(e) {

      GetAssayData(
        object,
        assay = assay,
        slot = "data"
      )
    }
  )
}

join_rna_if_possible <- function(object) {

  if (
    "JoinLayers" %in%
      getNamespaceExports(
        "SeuratObject"
      )
  ) {

    joined <- try(
      SeuratObject::JoinLayers(
        object[["RNA"]]
      ),
      silent = TRUE
    )

    if (
      !inherits(
        joined,
        "try-error"
      )
    ) {

      object[["RNA"]] <- joined
    }
  }

  object
}

safe_present <- function(
    genes,
    universe
) {

  intersect(
    unique(genes),
    universe
  )
}

row_zscore <- function(mat) {

  if (
    nrow(mat) == 0 ||
    ncol(mat) == 0
  ) {

    return(mat)
  }

  means <- Matrix::rowMeans(mat)

  sq_means <- Matrix::rowMeans(
    mat^2
  )

  sds <- sqrt(
    pmax(
      sq_means -
        means^2,
      0
    )
  )

  sds[
    !is.finite(sds) |
      sds == 0
  ] <- 1

  z <- sweep(
    mat,
    1,
    means,
    "-"
  )

  z <- sweep(
    z,
    1,
    sds,
    "/"
  )

  z
}

signature_score_from_z <- function(
    zmat,
    genes
) {

  genes <- intersect(
    genes,
    rownames(zmat)
  )

  if (
    length(
      genes
    ) == 0
  ) {

    return(
      rep(
        NA_real_,
        ncol(zmat)
      )
    )
  }

  Matrix::colMeans(
    zmat[
      genes,
      ,
      drop = FALSE
    ]
  )
}

# ------------------------------------------------------------
# 2. Marker definitions
# ------------------------------------------------------------

marker_sets <- list(

  Hepatocyte = c(
    "ALB",
    "TTR",
    "APOA1",
    "ASGR1",
    "HNF4A",
    "CYP3A4",
    "CPS1",
    "GLUL"
  ),

  LSEC = c(
    "CLEC4G",
    "FCGR2B",
    "STAB2",
    "STAB1",
    "EMCN",
    "PLVAP",
    "RAMP2",
    "KDR",
    "PECAM1",
    "VWF",
    "EGFL7",
    "ACKR3"
  ),

  Macrophage = c(
    "PTPRC",
    "LST1",
    "TYROBP",
    "FCER1G",
    "C1QC"
  ),

  Mesenchymal = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM"
  )
)

all_marker_genes <- unique(
  unlist(
    marker_sets,
    use.names = FALSE
  )
)

# ------------------------------------------------------------
# 3. Load / validate v1.3.2B checkpoint
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.2B-QC1\n")
cat("LSEC-vs-Hepatocyte targeted audit\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(checkpoint_rds)
)

hepato <- readRDS(
  checkpoint_rds
)

DefaultAssay(
  hepato
) <- "RNA"

required_meta <- c(
  celltype_col,
  sample_col,
  cluster_col,
  "current_hepatobiliary_label_v132B"
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    hepato@meta.data
  )
)

if (
  length(
    missing_meta
  ) > 0
) {

  stop(
    "Missing checkpoint metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  ncol(
    hepato
  ) != expected_hepatobiliary_cells
) {

  stop(
    "Checkpoint cell-count mismatch. Expected ",
    expected_hepatobiliary_cells,
    " but found ",
    ncol(
      hepato
    ),
    "."
  )
}

observed_clusters <- unique(
  as.character(
    hepato[[cluster_col]][, 1]
  )
)

if (
  !setequal(
    observed_clusters,
    expected_clusters
  )
) {

  stop(
    "Unexpected cluster set in checkpoint. Observed: ",
    paste(
      sort(
        observed_clusters
      ),
      collapse = ","
    )
  )
}

if (
  !reduction_name %in%
    Reductions(
      hepato
    )
) {

  stop(
    "Required UMAP reduction missing: ",
    reduction_name
  )
}

cat(
  "Hepatobiliary checkpoint cells:",
  ncol(
    hepato
  ),
  "\n"
)

cat(
  "Clusters:",
  length(
    observed_clusters
  ),
  "\n"
)

# ------------------------------------------------------------
# 4. Define target / hepatocyte reference cells
# ------------------------------------------------------------

cluster_vec <- as.character(
  hepato[[cluster_col]][, 1]
)

current_label_vec <- as.character(
  hepato$current_hepatobiliary_label_v132B
)

target0_cells <- colnames(
  hepato
)[
  cluster_vec == "0"
]

target7_cells <- colnames(
  hepato
)[
  cluster_vec == "7"
]

other_hepatocyte_cells <- colnames(
  hepato
)[
  current_label_vec == "Hepatocyte" &
    !cluster_vec %in%
      target_clusters
]

if (
  length(
    target0_cells
  ) == 0 ||
    length(
      target7_cells
    ) == 0 ||
    length(
      other_hepatocyte_cells
    ) == 0
) {

  stop(
    "One or more required target/reference groups are empty."
  )
}

cat(
  "Cluster 0 cells:",
  length(
    target0_cells
  ),
  "\n"
)

cat(
  "Cluster 7 cells:",
  length(
    target7_cells
  ),
  "\n"
)

cat(
  "Other Hepatocyte-labelled reference cells:",
  length(
    other_hepatocyte_cells
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 5. Load canonical object only to obtain validated
#    LSEC / Endothelial reference populations
# ------------------------------------------------------------

stopifnot(
  file.exists(
    canonical_rds
  )
)

canonical <- readRDS(
  canonical_rds
)

DefaultAssay(
  canonical
) <- "RNA"

if (
  ncol(
    canonical
  ) != expected_canonical_cells
) {

  stop(
    "Canonical cell-count mismatch. Expected ",
    expected_canonical_cells,
    " but found ",
    ncol(
      canonical
    ),
    "."
  )
}

if (
  !celltype_col %in%
    colnames(
      canonical@meta.data
    )
) {

  stop(
    "Canonical celltype metadata missing: ",
    celltype_col
  )
}

canonical_celltype <- as.character(
  canonical[[celltype_col]][, 1]
)

lsec_cells <- colnames(
  canonical
)[
  canonical_celltype == "LSEC"
]

endothelial_cells <- colnames(
  canonical
)[
  canonical_celltype == "Endothelial"
]

if (
  length(
    lsec_cells
  ) == 0 ||
    length(
      endothelial_cells
    ) == 0
) {

  stop(
    "Canonical LSEC or Endothelial reference group is empty."
  )
}

cat(
  "Canonical LSEC reference cells:",
  length(
    lsec_cells
  ),
  "\n"
)

cat(
  "Canonical Endothelial reference cells:",
  length(
    endothelial_cells
  ),
  "\n\n"
)

# Reduce memory: keep only reference cells from canonical
canonical_ref <- subset(
  canonical,
  cells = c(
    lsec_cells,
    endothelial_cells
  )
)

rm(
  canonical
)

gc()

canonical_ref <- join_rna_if_possible(
  canonical_ref
)

hepato <- join_rna_if_possible(
  hepato
)

# ------------------------------------------------------------
# 6. Gene availability audit
# ------------------------------------------------------------

hepato_genes <- rownames(
  hepato
)

ref_genes <- rownames(
  canonical_ref
)

gene_availability <- data.frame(
  gene =
    all_marker_genes,
  present_hepatobiliary =
    all_marker_genes %in%
      hepato_genes,
  present_canonical_reference =
    all_marker_genes %in%
      ref_genes,
  stringsAsFactors = FALSE
)

gene_availability$usable <-
  gene_availability$present_hepatobiliary &
  gene_availability$present_canonical_reference

write.csv(
  gene_availability,
  file.path(
    outdir,
    "01_marker_gene_availability.csv"
  ),
  row.names = FALSE
)

usable_genes <- gene_availability$gene[
  gene_availability$usable
]

if (
  length(
    usable_genes
  ) < 10
) {

  stop(
    "Too few QC marker genes are jointly available: ",
    length(
      usable_genes
    )
  )
}

# ------------------------------------------------------------
# 7. Extract RNA expression matrices
# ------------------------------------------------------------

hepato_data <- get_data_layer(
  hepato,
  assay = "RNA"
)

ref_data <- get_data_layer(
  canonical_ref,
  assay = "RNA"
)

hepato_audit_cells <- c(
  target0_cells,
  target7_cells,
  other_hepatocyte_cells
)

hepato_mat <- hepato_data[
  usable_genes,
  hepato_audit_cells,
  drop = FALSE
]

ref_audit_cells <- c(
  lsec_cells,
  endothelial_cells
)

ref_mat <- ref_data[
  usable_genes,
  ref_audit_cells,
  drop = FALSE
]

# Ensure identical gene order
ref_mat <- ref_mat[
  rownames(
    hepato_mat
  ),
  ,
  drop = FALSE
]

audit_mat <- cbind(
  hepato_mat,
  ref_mat
)

# ------------------------------------------------------------
# 8. Define audit groups
# ------------------------------------------------------------

audit_group <- c(
  rep(
    "Cluster0_target",
    length(
      target0_cells
    )
  ),
  rep(
    "Cluster7_target",
    length(
      target7_cells
    )
  ),
  rep(
    "Other_Hepatocyte_labelled",
    length(
      other_hepatocyte_cells
    )
  ),
  rep(
    "Canonical_LSEC_reference",
    length(
      lsec_cells
    )
  ),
  rep(
    "Canonical_Endothelial_reference",
    length(
      endothelial_cells
    )
  )
)

names(
  audit_group
) <- colnames(
  audit_mat
)

audit_group_levels <- c(
  "Cluster0_target",
  "Cluster7_target",
  "Other_Hepatocyte_labelled",
  "Canonical_LSEC_reference",
  "Canonical_Endothelial_reference"
)

audit_group <- factor(
  audit_group,
  levels =
    audit_group_levels
)

group_counts <- as.data.frame(
  table(
    group =
      audit_group
  ),
  stringsAsFactors = FALSE
)

names(
  group_counts
)[2] <- "n_cells"

write.csv(
  group_counts,
  file.path(
    outdir,
    "02_QC1_group_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Gene-level summary by audit group
# ------------------------------------------------------------

gene_summary_list <- list()

for (
  grp in audit_group_levels
) {

  idx <- which(
    audit_group == grp
  )

  mat_grp <- audit_mat[
    ,
    idx,
    drop = FALSE
  ]

  gene_summary_list[[grp]] <- data.frame(
    group = grp,
    gene =
      rownames(
        mat_grp
      ),
    mean_expression =
      Matrix::rowMeans(
        mat_grp
      ),
    pct_detected =
      100 *
        Matrix::rowMeans(
          mat_grp > 0
        ),
    stringsAsFactors = FALSE
  )
}

gene_summary <- do.call(
  rbind,
  gene_summary_list
)

rownames(
  gene_summary
) <- NULL

write.csv(
  gene_summary,
  file.path(
    outdir,
    "03_marker_expression_summary_by_group.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Gene-wise Z scoring and lineage signature scores
#
# This normalization is ONLY for the QC comparison.
# It allows LSEC vs hepatocyte signatures to be compared on a
# common per-gene scale.
# ------------------------------------------------------------

zmat <- row_zscore(
  audit_mat
)

marker_sets_usable <- lapply(
  marker_sets,
  function(x) {

    intersect(
      x,
      rownames(
        zmat
      )
    )
  }
)

score_df <- data.frame(
  cell =
    colnames(
      zmat
    ),
  group =
    as.character(
      audit_group
    ),
  Hepatocyte_score =
    signature_score_from_z(
      zmat,
      marker_sets_usable$Hepatocyte
    ),
  LSEC_score =
    signature_score_from_z(
      zmat,
      marker_sets_usable$LSEC
    ),
  Macrophage_score =
    signature_score_from_z(
      zmat,
      marker_sets_usable$Macrophage
    ),
  Mesenchymal_score =
    signature_score_from_z(
      zmat,
      marker_sets_usable$Mesenchymal
    ),
  stringsAsFactors = FALSE
)

score_df$LSEC_minus_Hepatocyte <-
  score_df$LSEC_score -
  score_df$Hepatocyte_score

write.csv(
  score_df,
  file.path(
    outdir,
    "04_cell_level_QC_signature_scores.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Group-level signature summary
# ------------------------------------------------------------

score_summary <- do.call(
  rbind,
  lapply(
    audit_group_levels,
    function(grp) {

      x <- score_df[
        score_df$group == grp,
        ,
        drop = FALSE
      ]

      data.frame(
        group = grp,
        n_cells =
          nrow(x),
        median_Hepatocyte_score =
          median(
            x$Hepatocyte_score,
            na.rm = TRUE
          ),
        median_LSEC_score =
          median(
            x$LSEC_score,
            na.rm = TRUE
          ),
        median_LSEC_minus_Hepatocyte =
          median(
            x$LSEC_minus_Hepatocyte,
            na.rm = TRUE
          ),
        pct_cells_LSEC_gt_Hepatocyte =
          100 *
            mean(
              x$LSEC_minus_Hepatocyte > 0,
              na.rm = TRUE
            ),
        median_Macrophage_score =
          median(
            x$Macrophage_score,
            na.rm = TRUE
          ),
        median_Mesenchymal_score =
          median(
            x$Mesenchymal_score,
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  score_summary,
  file.path(
    outdir,
    "05_group_level_QC_signature_summary.csv"
  ),
  row.names = FALSE
)

target_score_summary <- score_summary[
  score_summary$group %in%
    c(
      "Cluster0_target",
      "Cluster7_target"
    ),
  ,
  drop = FALSE
]

write.csv(
  target_score_summary,
  file.path(
    outdir,
    "06_target_cluster_QC_signature_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Export existing top20 markers for clusters 0/7
# ------------------------------------------------------------

if (
  file.exists(
    top20_csv
  )
) {

  top20 <- read.csv(
    top20_csv,
    check.names = FALSE
  )

  if (
    "cluster" %in%
      colnames(
        top20
      )
  ) {

    top20$cluster <- as.character(
      top20$cluster
    )

    target_top20 <- top20[
      top20$cluster %in%
        target_clusters,
      ,
      drop = FALSE
    ]

    write.csv(
      target_top20,
      file.path(
        outdir,
        "07_cluster0_7_existing_top20_markers.csv"
      ),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------
# 13. Fig01 - hepatobiliary UMAP with clusters 0/7 highlighted
# ------------------------------------------------------------

hepato$QC1_highlight <- "Other_hepatobiliary"

hepato$QC1_highlight[
  cluster_vec == "0"
] <- "Cluster_0"

hepato$QC1_highlight[
  cluster_vec == "7"
] <- "Cluster_7"

highlight_palette <- c(
  "Other_hepatobiliary" =
    "#D9D9D9",
  "Cluster_0" =
    "#FF1A1A",
  "Cluster_7" =
    "#0033FF"
)

p01 <- DimPlot(
  hepato,
  reduction =
    reduction_name,
  group.by =
    "QC1_highlight",
  cols =
    highlight_palette,
  pt.size = 0.10,
  raster = TRUE,
  order = FALSE
) +
  ggtitle(
    "v1.3.2B-QC1: target clusters 0 and 7"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig01_TargetClusters0_7_Highlight_UMAP.pdf"
  ),
  p01,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 14. Fig02 - DotPlot-like marker audit across target/reference
# ------------------------------------------------------------

plot_gene_order <- unlist(
  marker_sets,
  use.names = FALSE
)

plot_gene_order <- plot_gene_order[
  plot_gene_order %in%
    usable_genes
]

plot_gene_order <- unique(
  plot_gene_order
)

group_gene_z <- gene_summary

# Per-gene z-score of group mean expression
group_gene_z$mean_z <- ave(
  group_gene_z$mean_expression,
  group_gene_z$gene,
  FUN = function(x) {

    s <- sd(x)

    if (
      !is.finite(s) ||
      s == 0
    ) {

      return(
        rep(
          0,
          length(x)
        )
      )
    }

    (x - mean(x)) / s
  }
)

group_gene_z$group <- factor(
  group_gene_z$group,
  levels =
    rev(
      audit_group_levels
    )
)

group_gene_z$gene <- factor(
  group_gene_z$gene,
  levels =
    plot_gene_order
)

p02 <- ggplot(
  group_gene_z,
  aes(
    x = gene,
    y = group
  )
) +
  geom_point(
    aes(
      size =
        pct_detected,
      color =
        mean_z
    )
  ) +
  scale_size(
    range = c(
      0.5,
      8
    ),
    limits = c(
      0,
      100
    )
  ) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(
      -2,
      2
    ),
    oob =
      scales::squish
  ) +
  labs(
    title =
      "LSEC-vs-Hepatocyte targeted marker audit",
    subtitle =
      "Cluster 0/7 compared with canonical LSEC and Hepatocyte references",
    x = NULL,
    y = NULL,
    size =
      "% detected",
    color =
      "Relative\nmean expression"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig02_LSEC_vs_Hepatocyte_Marker_DotPlot.pdf"
  ),
  p02,
  width = 19,
  height = 6.5
)

write.csv(
  group_gene_z,
  file.path(
    outdir,
    "08_Fig02_marker_plot_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Fig03 - Hepatocyte vs LSEC score distributions
# ------------------------------------------------------------

score_plot_df <- score_df

score_plot_df$group <- factor(
  score_plot_df$group,
  levels =
    audit_group_levels
)

p03a <- ggplot(
  score_plot_df,
  aes(
    x = group,
    y =
      Hepatocyte_score
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title =
      "Hepatocyte signature score",
    x = NULL,
    y =
      "Mean gene-wise Z score"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

p03b <- ggplot(
  score_plot_df,
  aes(
    x = group,
    y =
      LSEC_score
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title =
      "LSEC signature score",
    x = NULL,
    y =
      "Mean gene-wise Z score"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig03A_HepatocyteSignature_Violin.pdf"
  ),
  p03a,
  width = 10,
  height = 6
)

ggsave(
  file.path(
    outdir,
    "Fig03B_LSECSignature_Violin.pdf"
  ),
  p03b,
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# 16. Fig04 - LSEC minus Hepatocyte score delta
# ------------------------------------------------------------

p04 <- ggplot(
  score_plot_df,
  aes(
    x = group,
    y =
      LSEC_minus_Hepatocyte
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  labs(
    title =
      "LSEC minus Hepatocyte signature score",
    subtitle =
      "Positive values indicate an LSEC-like transcriptomic balance",
    x = NULL,
    y =
      "LSEC score - Hepatocyte score"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig04_LSEC_minus_Hepatocyte_ScoreDelta.pdf"
  ),
  p04,
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# 17. Fig05 - targeted FeaturePlots on hepatobiliary UMAP
# ------------------------------------------------------------

feature_genes <- c(
  "ALB",
  "ASGR1",
  "HNF4A",
  "CYP3A4",
  "CLEC4G",
  "FCGR2B",
  "STAB2",
  "STAB1",
  "PLVAP",
  "RAMP2",
  "PECAM1",
  "VWF",
  "EGFL7",
  "ACKR3"
)

feature_present <- safe_present(
  feature_genes,
  rownames(
    hepato
  )
)

feature_availability <- data.frame(
  gene =
    feature_genes,
  present =
    feature_genes %in%
      rownames(
        hepato
      ),
  stringsAsFactors = FALSE
)

write.csv(
  feature_availability,
  file.path(
    outdir,
    "09_featureplot_gene_availability.csv"
  ),
  row.names = FALSE
)

if (
  length(
    feature_present
  ) > 0
) {

  p05 <- FeaturePlot(
    hepato,
    features =
      feature_present,
    reduction =
      reduction_name,
    ncol = 4,
    pt.size = 0.08,
    raster = TRUE,
    order = FALSE,
    cols = c(
      "#F2F2F2",
      "#FF1A1A"
    )
  )

  ggsave(
    file.path(
      outdir,
      "Fig05_LSEC_Hepatocyte_Sentinel_FeaturePlots.pdf"
    ),
    p05,
    width = 16,
    height = 16
  )
}

# ------------------------------------------------------------
# 18. Fig06 - macrophage / mesenchymal contamination check
# ------------------------------------------------------------

contam_df <- score_plot_df[
  score_plot_df$group %in%
    c(
      "Cluster0_target",
      "Cluster7_target",
      "Canonical_LSEC_reference",
      "Canonical_Endothelial_reference"
    ),
  ,
  drop = FALSE
]

contam_long <- rbind(
  data.frame(
    cell =
      contam_df$cell,
    group =
      contam_df$group,
    signature =
      "Macrophage",
    score =
      contam_df$Macrophage_score
  ),
  data.frame(
    cell =
      contam_df$cell,
    group =
      contam_df$group,
    signature =
      "Mesenchymal",
    score =
      contam_df$Mesenchymal_score
  )
)

contam_long$group <- factor(
  contam_long$group,
  levels =
    c(
      "Cluster0_target",
      "Cluster7_target",
      "Canonical_LSEC_reference",
      "Canonical_Endothelial_reference"
    )
)

p06 <- ggplot(
  contam_long,
  aes(
    x = group,
    y = score
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  facet_wrap(
    ~ signature,
    scales = "free_y"
  ) +
  labs(
    title =
      "Alternative-lineage contamination audit",
    x = NULL,
    y =
      "Mean gene-wise Z score"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig06_AlternativeLineage_Contamination_Audit.pdf"
  ),
  p06,
  width = 12,
  height = 6
)

# ------------------------------------------------------------
# 19. Sample contribution for target clusters 0 and 7
# ------------------------------------------------------------

target_sample <- data.frame(
  cluster =
    cluster_vec[
      cluster_vec %in%
        target_clusters
    ],
  sample =
    as.character(
      hepato[[sample_col]][, 1]
    )[
      cluster_vec %in%
        target_clusters
    ],
  stringsAsFactors = FALSE
)

target_sample_tab <- as.data.frame(
  table(
    cluster =
      target_sample$cluster,
    sample =
      target_sample$sample
  ),
  stringsAsFactors = FALSE
)

target_sample_tab$cluster_total <- ave(
  target_sample_tab$Freq,
  target_sample_tab$cluster,
  FUN = sum
)

target_sample_tab$within_cluster_pct <- 100 *
  target_sample_tab$Freq /
  target_sample_tab$cluster_total

write.csv(
  target_sample_tab,
  file.path(
    outdir,
    "10_target_cluster_sample_contribution.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 20. Conservative decision-support classification
#
# This is NOT a final cell-type annotation.
# It only summarizes which reference axis each target is closer
# to on the QC signature score.
# ------------------------------------------------------------

ref_summary <- score_summary[
  score_summary$group %in%
    c(
      "Other_Hepatocyte_labelled",
      "Canonical_LSEC_reference"
    ),
  ,
  drop = FALSE
]

hep_ref_delta <- ref_summary$median_LSEC_minus_Hepatocyte[
  ref_summary$group ==
    "Other_Hepatocyte_labelled"
]

lsec_ref_delta <- ref_summary$median_LSEC_minus_Hepatocyte[
  ref_summary$group ==
    "Canonical_LSEC_reference"
]

target_decision <- target_score_summary

target_decision$distance_to_hepatocyte_reference <-
  abs(
    target_decision$median_LSEC_minus_Hepatocyte -
      hep_ref_delta
  )

target_decision$distance_to_LSEC_reference <-
  abs(
    target_decision$median_LSEC_minus_Hepatocyte -
      lsec_ref_delta
  )

target_decision$QC_axis_closer_to <- ifelse(
  target_decision$distance_to_LSEC_reference <
    target_decision$distance_to_hepatocyte_reference,
  "LSEC_reference",
  "Hepatocyte_reference"
)

write.csv(
  target_decision,
  file.path(
    outdir,
    "11_target_cluster_QC_decision_support.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 21. Review flags
# ------------------------------------------------------------

flags <- c(
  "Human MASH v1.3.2B-QC1",
  "=======================",
  "",
  "Targeted LSEC-vs-Hepatocyte audit for hepatobiliary clusters 0 and 7.",
  "",
  paste0(
    "Cluster 0 cells: ",
    length(
      target0_cells
    )
  ),
  paste0(
    "Cluster 7 cells: ",
    length(
      target7_cells
    )
  ),
  paste0(
    "Other Hepatocyte-labelled reference cells: ",
    length(
      other_hepatocyte_cells
    )
  ),
  paste0(
    "Canonical LSEC reference cells: ",
    length(
      lsec_cells
    )
  ),
  paste0(
    "Canonical Endothelial reference cells: ",
    length(
      endothelial_cells
    )
  ),
  "",
  "TARGET SCORE SUMMARY:",
  paste(
    capture.output(
      print(
        target_decision,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION RULE:",
  "  QC_axis_closer_to is decision support only.",
  "  Final identity must agree with marker DotPlot, FeaturePlots,",
  "  existing top markers, target sample distribution, and contamination audit.",
  "",
  "IMPORTANT:",
  "  No RPCA/reclustering/FindAllMarkers was rerun.",
  "  Canonical v1.3 object was not modified.",
  "  v1.3.2B checkpoint object was not modified.",
  "",
  "NEXT STEP:",
  "  If clusters 0/7 are resolved, proceed to v1.3.2B.1 manual hepatobiliary annotation/cleaning."
)

writeLines(
  flags,
  file.path(
    outdir,
    "18_v1.3.2B_QC1_review_flags.txt"
  )
)

cat(
  paste(
    flags,
    collapse = "\n"
  )
)

cat("\n\n")

# ------------------------------------------------------------
# 22. Manifest / provenance
# ------------------------------------------------------------

manifest <- c(
  "version=v1.3.2B-QC1",
  paste0(
    "checkpoint_input=",
    checkpoint_rds
  ),
  paste0(
    "checkpoint_md5=",
    unname(
      tools::md5sum(
        checkpoint_rds
      )
    )
  ),
  paste0(
    "canonical_reference_input=",
    canonical_rds
  ),
  paste0(
    "canonical_reference_md5=",
    unname(
      tools::md5sum(
        canonical_rds
      )
    )
  ),
  "target_clusters=0,7",
  "RPCA_recomputed=FALSE",
  "clustering_recomputed=FALSE",
  "FindAllMarkers_recomputed=FALSE",
  "canonical_object_modified=FALSE",
  "checkpoint_object_modified=FALSE",
  "status=TARGETED_QC_REVIEW_REQUIRED"
)

writeLines(
  manifest,
  file.path(
    outdir,
    "98_manifest.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)

cat("============================================================\n")
cat("v1.3.2B-QC1 COMPLETE\n")
cat("TARGETED LSEC-vs-HEPATOCYTE AUDIT COMPLETE\n")
cat("NO RPCA / RECLUSTERING / FindAllMarkers RERUN\n")
cat("CANONICAL v1.3 OBJECT NOT MODIFIED\n")
cat("v1.3.2B CHECKPOINT NOT MODIFIED\n")
cat("\n")
cat("Output directory:\n")
cat(outdir, "\n")
cat("============================================================\n")
