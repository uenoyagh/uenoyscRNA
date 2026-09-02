#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6100)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# HSC_Mesenchymal subclustering
#
# Version: v6.1.0
#
# BASE CHECKPOINT
#   mouse-mash-mphi-staining-markers-v5.9.3
#   commit 1381eeae98099a835745313be30c7c6bbc7d1674
#
# INPUT
#   Frozen whole-cell parent v5.1.1
#
# SAMPLES
#   Sham1, Sham20, Tx17, Tx5
#
# PURPOSE
#   1) Extract HSC_Mesenchymal from the frozen whole-cell object.
#   2) Recompute HSC-specific RNA normalization / HVGs / PCA.
#   3) Perform sample-aware RPCA integration.
#   4) Recluster at multiple resolutions.
#   5) Evaluate qHSC / ECM-activated / contractile / proliferating / inflammatory
#      programs without forcing a final annotation.
#   6) Export cluster markers, module-score summaries, UMAPs, DotPlots,
#      sample/condition compositions, and contamination-audit plots.
#
# IMPORTANT
#   - No modification of the frozen parent object.
#   - No final qHSC/aHSC label is written in v6.1.0.
#   - Final HSC annotation will be decided after inspecting v6.1.0 outputs.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(p, file, width, height) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}

resolve_first <- function(x, candidates, what) {
  hit <- candidates[candidates %in% x]

  if (!length(hit)) {
    stop(
      "Could not resolve ",
      what,
      ". Candidates: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[[1]]
}

canonical_sample <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    grepl("^Sham1$", x, ignore.case = TRUE) ~ "Sham1",
    grepl("^Sham20$", x, ignore.case = TRUE) ~ "Sham20",
    grepl("^Tx17$", x, ignore.case = TRUE) ~ "Tx17",
    grepl("^Tx5$", x, ignore.case = TRUE) ~ "Tx5",
    grepl("Sham.?1", x, ignore.case = TRUE) ~ "Sham1",
    grepl("Sham.?20", x, ignore.case = TRUE) ~ "Sham20",
    grepl("Tx.?17", x, ignore.case = TRUE) ~ "Tx17",
    grepl("Tx.?5", x, ignore.case = TRUE) ~ "Tx5",
    TRUE ~ NA_character_
  )
}

canonical_condition <- function(sample) {
  dplyr::case_when(
    grepl("^Sham", sample, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", sample, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )
}

resolve_hsc_label <- function(labels) {

  preferred <- c(
    "HSC_Mesenchymal",
    "HSC/Mesenchymal",
    "HSC-Mesenchymal",
    "HSC Mesenchymal"
  )

  exact <- preferred[preferred %in% labels]

  if (length(exact)) {
    return(exact[[1]])
  }

  hit <- labels[
    grepl(
      "HSC.*Mesench|Mesench.*HSC|Hepatic.?stellate",
      labels,
      ignore.case = TRUE
    )
  ]

  if (!length(hit)) {
    stop(
      "Could not resolve HSC_Mesenchymal label. Available labels:\n",
      paste(sort(labels), collapse = " | ")
    )
  }

  hit[[1]]
}

present_genes <- function(object, genes) {
  genes[genes %in% rownames(object)]
}

safe_join_rna_layers <- function(object) {

  if (!"RNA" %in% Assays(object)) {
    stop("RNA assay not found.")
  }

  current_layers <- Layers(object[["RNA"]])

  msg(
    "RNA layers before JoinLayers: ",
    paste(current_layers, collapse = ", ")
  )

  if (length(current_layers) > 1) {
    object[["RNA"]] <- JoinLayers(object[["RNA"]])
  }

  msg(
    "RNA layers after JoinLayers: ",
    paste(Layers(object[["RNA"]]), collapse = ", ")
  )

  object
}

add_program_score <- function(object, genes, score_name) {

  genes_use <- present_genes(object, genes)

  if (length(genes_use) < 2) {
    warning(
      "Program ",
      score_name,
      " has <2 genes present and will be skipped."
    )
    return(object)
  }

  object <- AddModuleScore(
    object = object,
    features = list(genes_use),
    assay = "RNA",
    name = paste0(score_name, "_tmp"),
    seed = 6100
  )

  generated_col <- paste0(score_name, "_tmp1")
  final_col <- paste0(score_name, "_v610")

  object[[final_col]] <- object@meta.data[[generated_col]]
  object[[generated_col]] <- NULL

  object
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_Subclustering_v6.1.0"
)

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

for (d in c(
  OUT,
  RDS_OUT,
  TAB_OUT,
  FIG_OUT,
  LOG_OUT
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

WHOLE_LAYER1_COL <- "wholecell_layer1_FINAL_v511"

SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "sample_FIXED2",
  "sample",
  "orig.ident"
)

SAMPLE_LEVELS <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

NFEATURES <- 3000
N_PCS <- 40
DIMS_USE <- 1:30

RESOLUTIONS <- c(
  0.4,
  0.6,
  0.8,
  1.0,
  1.2
)

WORKING_RESOLUTION <- 0.8

HSC_UMAP <- "umap.hsc.rpca"

PROGRAMS <- list(

  "Quiescent" = c(
    "Lrat",
    "Rbp1",
    "Reln",
    "Cygb",
    "Des"
  ),

  "ECM_activated" = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Fn1",
    "Lum",
    "Dcn",
    "Bgn",
    "Timp1"
  ),

  "Contractile" = c(
    "Acta2",
    "Tagln",
    "Myl9",
    "Cnn1",
    "Tpm2",
    "Myh11"
  ),

  "Proliferating" = c(
    "Mki67",
    "Top2a",
    "Ccnb1",
    "Ccna2",
    "Cdk1"
  ),

  "Inflammatory_stress" = c(
    "Ccl2",
    "Cxcl1",
    "Il6",
    "Icam1",
    "Vcam1",
    "Nfkbia"
  )
)

CANONICAL_MARKERS <- unique(
  unlist(
    PROGRAMS,
    use.names = FALSE
  )
)

CONTAMINATION_MARKERS <- list(

  "Endothelial" = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Vwf"
  ),

  "Macrophage" = c(
    "Ptprc",
    "Lyz2",
    "Adgre1",
    "C1qa",
    "Cd68"
  ),

  "Hepatocyte" = c(
    "Alb",
    "Apoa1",
    "Ttr",
    "Cyp2f2"
  ),

  "Cholangiocyte" = c(
    "Krt19",
    "Krt8",
    "Krt18",
    "Epcam"
  )
)

CONTAMINATION_GENES <- unique(
  unlist(
    CONTAMINATION_MARKERS,
    use.names = FALSE
  )
)


# ==============================================================================
# 4. Load frozen parent and extract HSC
# ==============================================================================

if (!file.exists(WHOLE_RDS)) {
  stop(
    "Whole-cell frozen RDS not found: ",
    WHOLE_RDS
  )
}

msg("Loading frozen whole-cell parent...")
whole <- readRDS(
  WHOLE_RDS
)

if (!"RNA" %in% Assays(whole)) {
  stop("RNA assay missing from whole-cell parent.")
}

if (!WHOLE_LAYER1_COL %in% colnames(whole@meta.data)) {
  stop(
    "Missing whole-cell annotation column: ",
    WHOLE_LAYER1_COL
  )
}

SAMPLE_COL <- resolve_first(
  colnames(whole@meta.data),
  SAMPLE_CANDIDATES,
  "whole-cell sample column"
)

whole$sample_hsc_v610 <- canonical_sample(
  whole@meta.data[[SAMPLE_COL]]
)

whole$condition_hsc_v610 <- canonical_condition(
  whole$sample_hsc_v610
)

layer1 <- as.character(
  whole@meta.data[[WHOLE_LAYER1_COL]]
)

HSC_SOURCE_LABEL <- resolve_hsc_label(
  unique(layer1)
)

msg(
  "HSC source label: ",
  HSC_SOURCE_LABEL
)

hsc_cells <- colnames(whole)[
  layer1 == HSC_SOURCE_LABEL &
    whole$sample_hsc_v610 %in% SAMPLE_LEVELS
]

if (!length(hsc_cells)) {
  stop(
    "No HSC cells found for Sham1/Sham20/Tx17/Tx5."
  )
}

hsc <- subset(
  whole,
  cells = hsc_cells
)

rm(whole)
gc()

DefaultAssay(hsc) <- "RNA"

hsc$sample_hsc_v610 <- factor(
  hsc$sample_hsc_v610,
  levels = SAMPLE_LEVELS
)

hsc$condition_hsc_v610 <- factor(
  hsc$condition_hsc_v610,
  levels = c(
    "Sham",
    "Tx"
  )
)

if (!"percent.mt" %in% colnames(hsc@meta.data)) {
  hsc[["percent.mt"]] <- PercentageFeatureSet(
    hsc,
    pattern = "^mt-"
  )
}

msg(
  "HSC cells extracted: ",
  ncol(hsc)
)

msg("HSC sample counts:")
print(
  table(
    hsc$sample_hsc_v610,
    useNA = "ifany"
  )
)

write.csv(
  hsc@meta.data %>%
    as_tibble(
      rownames = "cell"
    ) %>%
    count(
      sample_hsc_v610,
      condition_hsc_v610,
      name = "n_cells"
    ),
  file.path(
    TAB_OUT,
    "01_HSC_cells_by_sample_v6.1.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 5. RNA preparation
# ==============================================================================

msg("Joining RNA layers before HSC-specific integration...")
hsc <- safe_join_rna_layers(
  hsc
)

# Split RNA assay by biological sample for RPCA integration.
hsc[["RNA"]] <- split(
  hsc[["RNA"]],
  f = hsc$sample_hsc_v610
)

msg(
  "RNA layers after sample split: ",
  paste(
    Layers(hsc[["RNA"]]),
    collapse = ", "
  )
)


# ==============================================================================
# 6. HSC-specific normalization / HVG / PCA
# ==============================================================================

msg("NormalizeData...")
hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

msg("FindVariableFeatures...")
hsc <- FindVariableFeatures(
  hsc,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = NFEATURES,
  verbose = FALSE
)

msg("ScaleData...")
hsc <- ScaleData(
  hsc,
  assay = "RNA",
  features = VariableFeatures(hsc),
  verbose = FALSE
)

msg("RunPCA...")
hsc <- RunPCA(
  hsc,
  assay = "RNA",
  features = VariableFeatures(hsc),
  npcs = N_PCS,
  reduction.name = "pca.hsc",
  reduction.key = "HSCPCA_",
  verbose = FALSE
)

# Elbow
p_elbow <- ElbowPlot(
  hsc,
  reduction = "pca.hsc",
  ndims = N_PCS
) +
  ggtitle(
    "HSC-specific PCA"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_elbow,
  file.path(
    FIG_OUT,
    "01_HSC_PCA_elbow_v6.1.0.pdf"
  ),
  7,
  5
)


# ==============================================================================
# 7. RPCA integration
# ==============================================================================

msg("IntegrateLayers using RPCAIntegration...")

hsc <- IntegrateLayers(
  object = hsc,
  method = RPCAIntegration,
  orig.reduction = "pca.hsc",
  new.reduction = "integrated.hsc.rpca",
  verbose = FALSE
)

if (!"integrated.hsc.rpca" %in% Reductions(hsc)) {
  stop(
    "RPCA integration did not create integrated.hsc.rpca."
  )
}

# Join RNA layers again for marker testing.
msg("Joining RNA layers after RPCA integration...")
hsc <- safe_join_rna_layers(
  hsc
)


# ==============================================================================
# 8. HSC clustering / UMAP
# ==============================================================================

msg("FindNeighbors...")
hsc <- FindNeighbors(
  hsc,
  reduction = "integrated.hsc.rpca",
  dims = DIMS_USE,
  graph.name = c(
    "hsc_nn_v610",
    "hsc_snn_v610"
  ),
  verbose = FALSE
)

for (res in RESOLUTIONS) {

  res_label <- gsub(
    "\\.",
    "_",
    as.character(res)
  )

  cluster_col <- paste0(
    "hsc_rpca_res_",
    res_label,
    "_v610"
  )

  msg(
    "FindClusters resolution=",
    res
  )

  hsc <- FindClusters(
    hsc,
    graph.name = "hsc_snn_v610",
    resolution = res,
    algorithm = 1,
    cluster.name = cluster_col,
    random.seed = 6100,
    verbose = FALSE
  )
}

msg("RunUMAP...")
hsc <- RunUMAP(
  hsc,
  reduction = "integrated.hsc.rpca",
  dims = DIMS_USE,
  reduction.name = HSC_UMAP,
  reduction.key = "HSCUMAP_",
  seed.use = 6100,
  verbose = FALSE
)

WORKING_CLUSTER_COL <- paste0(
  "hsc_rpca_res_",
  gsub(
    "\\.",
    "_",
    as.character(
      WORKING_RESOLUTION
    )
  ),
  "_v610"
)

if (!WORKING_CLUSTER_COL %in% colnames(hsc@meta.data)) {
  stop(
    "Working cluster column not found: ",
    WORKING_CLUSTER_COL
  )
}

hsc$hsc_cluster_working_v610 <- factor(
  hsc@meta.data[[WORKING_CLUSTER_COL]]
)

Idents(
  hsc
) <- "hsc_cluster_working_v610"


# ==============================================================================
# 9. Resolution comparison UMAP
# ==============================================================================

msg("Creating resolution comparison UMAP...")

resolution_plots <- lapply(
  RESOLUTIONS,
  function(res) {

    col <- paste0(
      "hsc_rpca_res_",
      gsub(
        "\\.",
        "_",
        as.character(res)
      ),
      "_v610"
    )

    DimPlot(
      hsc,
      reduction = HSC_UMAP,
      group.by = col,
      label = TRUE,
      repel = TRUE,
      pt.size = 0.28,
      raster = FALSE
    ) +
      ggtitle(
        paste0(
          "HSC RPCA res ",
          res
        )
      ) +
      theme_classic(
        base_size = 8
      )
  }
)

p_res <- wrap_plots(
  resolution_plots,
  ncol = 2
) +
  plot_annotation(
    title =
      "HSC subclustering | resolution comparison"
  )

save_pdf(
  p_res,
  file.path(
    FIG_OUT,
    "02_HSC_resolution_comparison_UMAP_v6.1.0.pdf"
  ),
  13,
  14
)


# ==============================================================================
# 10. Working-resolution UMAP by cluster / sample / condition
# ==============================================================================

p_cluster <- DimPlot(
  hsc,
  reduction = HSC_UMAP,
  group.by = "hsc_cluster_working_v610",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.32,
  raster = FALSE
) +
  ggtitle(
    paste0(
      "HSC working resolution ",
      WORKING_RESOLUTION
    )
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_cluster,
  file.path(
    FIG_OUT,
    "03_HSC_working_resolution_UMAP_v6.1.0.pdf"
  ),
  8,
  7
)

p_sample <- DimPlot(
  hsc,
  reduction = HSC_UMAP,
  group.by = "hsc_cluster_working_v610",
  split.by = "sample_hsc_v610",
  pt.size = 0.25,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "HSC clusters | Sham1 / Sham20 / Tx17 / Tx5"
  )

save_pdf(
  p_sample,
  file.path(
    FIG_OUT,
    "04_HSC_working_resolution_UMAP_by_sample_v6.1.0.pdf"
  ),
  13,
  10
)

p_condition <- DimPlot(
  hsc,
  reduction = HSC_UMAP,
  group.by = "hsc_cluster_working_v610",
  split.by = "condition_hsc_v610",
  pt.size = 0.25,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "HSC clusters | Sham vs Tx"
  )

save_pdf(
  p_condition,
  file.path(
    FIG_OUT,
    "05_HSC_working_resolution_UMAP_Sham_vs_Tx_v6.1.0.pdf"
  ),
  13,
  6
)


# ==============================================================================
# 11. Program scores
# ==============================================================================

msg("Calculating HSC program scores...")

for (program_name in names(PROGRAMS)) {
  hsc <- add_program_score(
    hsc,
    PROGRAMS[[program_name]],
    program_name
  )
}

PROGRAM_SCORE_COLS <- paste0(
  names(PROGRAMS),
  "_v610"
)

PROGRAM_SCORE_COLS <- PROGRAM_SCORE_COLS[
  PROGRAM_SCORE_COLS %in%
    colnames(hsc@meta.data)
]

program_summary <- hsc@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  group_by(
    hsc_cluster_working_v610
  ) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(PROGRAM_SCORE_COLS),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write.csv(
  program_summary,
  file.path(
    TAB_OUT,
    "02_HSC_cluster_program_scores_v6.1.0.csv"
  ),
  row.names = FALSE
)

# Dominant program is exploratory only.
if (length(PROGRAM_SCORE_COLS)) {

  score_matrix <- as.matrix(
    hsc@meta.data[
      ,
      PROGRAM_SCORE_COLS,
      drop = FALSE
    ]
  )

  dominant_idx <- max.col(
    score_matrix,
    ties.method = "first"
  )

  hsc$dominant_program_v610 <- sub(
    "_v610$",
    "",
    PROGRAM_SCORE_COLS[
      dominant_idx
    ]
  )

  p_program <- DimPlot(
    hsc,
    reduction = HSC_UMAP,
    group.by = "dominant_program_v610",
    pt.size = 0.30,
    raster = FALSE
  ) +
    ggtitle(
      "Exploratory dominant HSC program"
    ) +
    theme_classic(
      base_size = 9
    )

  save_pdf(
    p_program,
    file.path(
      FIG_OUT,
      "06_HSC_dominant_program_UMAP_v6.1.0.pdf"
    ),
    8,
    7
  )
}


# ==============================================================================
# 12. Canonical marker FeaturePlots
# ==============================================================================

canonical_present <- present_genes(
  hsc,
  CANONICAL_MARKERS
)

write.csv(
  tibble(
    requested_gene = CANONICAL_MARKERS,
    present = CANONICAL_MARKERS %in%
      rownames(hsc)
  ),
  file.path(
    TAB_OUT,
    "03_HSC_canonical_marker_gene_audit_v6.1.0.csv"
  ),
  row.names = FALSE
)

if (length(canonical_present)) {

  p_features <- FeaturePlot(
    hsc,
    features = canonical_present,
    reduction = HSC_UMAP,
    ncol = 4,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q98",
    cols = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    ),
    raster = FALSE,
    pt.size = 0.25
  ) &
    theme_classic(
      base_size = 7
    )

  save_pdf(
    p_features,
    file.path(
      FIG_OUT,
      "07_HSC_canonical_marker_FeaturePlots_v6.1.0.pdf"
    ),
    15,
    3.7 * ceiling(
      length(canonical_present) / 4
    )
  )
}


# ==============================================================================
# 13. DotPlot
# ==============================================================================

if (length(canonical_present)) {

  p_dot <- DotPlot(
    hsc,
    features = canonical_present,
    group.by = "hsc_cluster_working_v610",
    assay = "RNA"
  ) +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    RotatedAxis() +
    labs(
      title =
        "HSC canonical programs by cluster",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_dot,
    file.path(
      FIG_OUT,
      "08_HSC_canonical_marker_DotPlot_v6.1.0.pdf"
    ),
    16,
    6
  )
}


# ==============================================================================
# 14. Program-score violin
# ==============================================================================

if (length(PROGRAM_SCORE_COLS)) {

  p_vln <- VlnPlot(
    hsc,
    features = PROGRAM_SCORE_COLS,
    group.by = "hsc_cluster_working_v610",
    pt.size = 0,
    ncol = 2
  ) &
    theme_classic(
      base_size = 8
    ) &
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_pdf(
    p_vln,
    file.path(
      FIG_OUT,
      "09_HSC_program_score_Violin_v6.1.0.pdf"
    ),
    13,
    10
  )
}


# ==============================================================================
# 15. Cluster markers
# ==============================================================================

msg("FindAllMarkers for working resolution...")

DefaultAssay(hsc) <- "RNA"
Idents(hsc) <- "hsc_cluster_working_v610"

markers_all <- FindAllMarkers(
  hsc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

write.csv(
  markers_all,
  file.path(
    TAB_OUT,
    "04_HSC_cluster_markers_ALL_v6.1.0.csv"
  ),
  row.names = FALSE
)

fc_col <- if (
  "avg_log2FC" %in% colnames(markers_all)
) {
  "avg_log2FC"
} else if (
  "avg_logFC" %in% colnames(markers_all)
) {
  "avg_logFC"
} else {
  stop(
    "Could not find average logFC column in FindAllMarkers output."
  )
}

markers_top20 <- markers_all %>%
  group_by(
    cluster
  ) %>%
  arrange(
    desc(
      .data[[fc_col]]
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write.csv(
  markers_top20,
  file.path(
    TAB_OUT,
    "05_HSC_cluster_markers_TOP20_v6.1.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Cluster composition by sample / condition
# ==============================================================================

cluster_by_sample <- hsc@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    sample_hsc_v610,
    condition_hsc_v610,
    hsc_cluster_working_v610,
    name = "n_cells"
  ) %>%
  group_by(
    sample_hsc_v610
  ) %>%
  mutate(
    fraction_of_HSC =
      n_cells / sum(n_cells)
  ) %>%
  ungroup()

write.csv(
  cluster_by_sample,
  file.path(
    TAB_OUT,
    "06_HSC_cluster_fraction_by_sample_v6.1.0.csv"
  ),
  row.names = FALSE
)

p_fraction <- ggplot(
  cluster_by_sample,
  aes(
    x = sample_hsc_v610,
    y = fraction_of_HSC,
    group = hsc_cluster_working_v610,
    linetype = hsc_cluster_working_v610
  )
) +
  geom_line(
    linewidth = 0.55
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    title =
      "HSC cluster fraction by biological sample",
    x = NULL,
    y = "Fraction of HSC",
    linetype = "Cluster"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_pdf(
  p_fraction,
  file.path(
    FIG_OUT,
    "10_HSC_cluster_fraction_by_sample_v6.1.0.pdf"
  ),
  9,
  6
)


# ==============================================================================
# 17. Cluster x sample heatmap
# ==============================================================================

p_heat <- ggplot(
  cluster_by_sample,
  aes(
    x = hsc_cluster_working_v610,
    y = sample_hsc_v610,
    fill = fraction_of_HSC
  )
) +
  geom_tile(
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = n_cells
    ),
    size = 3
  ) +
  scale_fill_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) +
  labs(
    title =
      "HSC cluster composition | cell count labels",
    x = "HSC cluster",
    y = NULL,
    fill = "Fraction"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_pdf(
  p_heat,
  file.path(
    FIG_OUT,
    "11_HSC_cluster_by_sample_heatmap_v6.1.0.pdf"
  ),
  10,
  5
)


# ==============================================================================
# 18. Contamination audit
# ==============================================================================

contamination_present <- present_genes(
  hsc,
  CONTAMINATION_GENES
)

write.csv(
  tibble(
    requested_gene = CONTAMINATION_GENES,
    present = CONTAMINATION_GENES %in%
      rownames(hsc)
  ),
  file.path(
    TAB_OUT,
    "07_HSC_contamination_marker_gene_audit_v6.1.0.csv"
  ),
  row.names = FALSE
)

if (length(contamination_present)) {

  p_contam <- DotPlot(
    hsc,
    features = contamination_present,
    group.by = "hsc_cluster_working_v610",
    assay = "RNA"
  ) +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    RotatedAxis() +
    labs(
      title =
        "HSC subcluster contamination audit",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_contam,
    file.path(
      FIG_OUT,
      "12_HSC_contamination_audit_DotPlot_v6.1.0.pdf"
    ),
    14,
    6
  )
}


# ==============================================================================
# 19. Save HSC subclustered RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_HSC_subclustered_v6.1.0.rds"
)

saveRDS(
  hsc,
  RDS_FILE,
  compress = FALSE
)

msg(
  "Saved HSC subclustered RDS: ",
  RDS_FILE
)


# ==============================================================================
# 20. Metadata / audit
# ==============================================================================

resolution_manifest <- tibble(
  resolution = RESOLUTIONS,
  cluster_column = paste0(
    "hsc_rpca_res_",
    gsub(
      "\\.",
      "_",
      as.character(
        RESOLUTIONS
      )
    ),
    "_v610"
  ),
  working_resolution =
    RESOLUTIONS == WORKING_RESOLUTION
)

write.csv(
  resolution_manifest,
  file.path(
    TAB_OUT,
    "08_HSC_resolution_manifest_v6.1.0.csv"
  ),
  row.names = FALSE
)

analysis_metadata <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "HSC_source_label",
    "sample_column",
    "n_HSC_cells",
    "nfeatures",
    "n_pcs",
    "dims_use",
    "working_resolution",
    "RPCA_reduction",
    "UMAP_reduction"
  ),
  value = c(
    "v6.1.0",
    WHOLE_RDS,
    HSC_SOURCE_LABEL,
    SAMPLE_COL,
    as.character(ncol(hsc)),
    as.character(NFEATURES),
    as.character(N_PCS),
    paste(range(DIMS_USE), collapse = ":"),
    as.character(WORKING_RESOLUTION),
    "integrated.hsc.rpca",
    HSC_UMAP
  )
)

write.csv(
  analysis_metadata,
  file.path(
    LOG_OUT,
    "analysis_metadata_v6.1.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_OUT,
    "sessionInfo_v6.1.0.txt"
  )
)

msg("DONE.")
msg("Output directory: ", OUT)
