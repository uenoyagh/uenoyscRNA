#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6110)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# HSC refined annotation after v6.1.0 subclustering
#
# Version: v6.1.2
#
# INPUT
#   Mouse_MASH_HSC_subclustered_v6.1.0.rds
#
# v6.1.2 CHANGE FROM v6.1.1
#   - Fixes over-stringent HSC purity filtering that could classify every
#     cluster as Review_boundary / Excluded_nonHSC.
#   - HSC identity now requires >=1 detected HSC-core gene per cell
#     (cluster-level evidence remains the main criterion).
#   - Ambient-prone hepatocyte genes Alb/Apoa1/Ttr are NOT used for exclusion.
#   - Broad epithelial Krt8/Krt18 are NOT used for cholangiocyte exclusion.
#   - Contamination exclusion now requires lineage-specific multi-gene evidence
#     AND comparison against HSC-core evidence.
#   - No reclustering / reintegration / UMAP recalculation.
#
# PURPOSE
#   1) Quantitatively identify non-HSC contamination/boundary clusters.
#   2) Retain a conservative genuine-HSC core.
#   3) Assign genuine HSC clusters to:
#        qHSC
#        ECM-activated HSC
#        Contractile HSC
#   4) Simultaneously create a two-state annotation:
#        qHSC
#        aHSC = ECM-activated + Contractile
#   5) Export sample-level cell counts/fractions for future CellChat.
#
# DESIGN PRINCIPLES
#   - No reclustering.
#   - No reintegration.
#   - No new UMAP.
#   - Existing HSC RPCA/UMAP coordinates from v6.1.0 are retained.
#   - Contamination is evaluated by MULTI-GENE co-expression, not by a single
#     ambient RNA-prone marker.
#   - "Review_boundary" cells are excluded from the primary interaction-ready
#     clean HSC object but are preserved in the annotated audit RDS.
#
# HSC STATES
#   qHSC:
#     Lrat, Rbp1, Reln, Cygb, Des
#
#   ECM-activated HSC:
#     Col1a1, Col1a2, Col3a1, Col5a1, Fn1, Lum, Dcn, Bgn, Timp1
#
#   Contractile HSC:
#     Acta2, Tagln, Myl9, Cnn1, Tpm2, Myh11
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

present_genes <- function(object, genes) {
  intersect(
    genes,
    rownames(object)
  )
}

row_z <- function(x) {
  sx <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(sx) ||
    sx == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / sx
}

add_module_score_clean <- function(
  object,
  genes,
  final_name
) {

  genes_use <- present_genes(
    object,
    genes
  )

  if (
    length(
      genes_use
    ) < 2
  ) {
    stop(
      "Too few genes available for module: ",
      final_name
    )
  }

  tmp_name <- paste0(
    final_name,
    "_TMP_"
  )

  object <- AddModuleScore(
    object,
    features = list(
      genes_use
    ),
    assay = "RNA",
    name = tmp_name,
    seed = 6110
  )

  generated <- paste0(
    tmp_name,
    "1"
  )

  object[[final_name]] <-
    object@meta.data[[generated]]

  object[[generated]] <- NULL

  object
}

multi_gene_hit_count <- function(
  counts,
  genes
) {

  genes_use <- intersect(
    genes,
    rownames(counts)
  )

  if (
    !length(
      genes_use
    )
  ) {
    return(
      rep(
        0L,
        ncol(counts)
      )
    )
  }

  Matrix::colSums(
    counts[
      genes_use,
      ,
      drop = FALSE
    ] > 0
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_Subclustering_v6.1.0",
  "RDS",
  "Mouse_MASH_HSC_subclustered_v6.1.0.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2"
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

for (
  d in c(
    OUT,
    RDS_OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

CLUSTER_COL <- "hsc_cluster_working_v610"

SAMPLE_COL <- "sample_hsc_v610"

CONDITION_COL <- "condition_hsc_v610"

HSC_UMAP <- "umap.hsc.rpca"

SAMPLE_LEVELS <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

STATE3_LEVELS <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

STATE2_LEVELS <- c(
  "qHSC",
  "aHSC"
)

# ------------------------------------------------------------------
# HSC programs
# ------------------------------------------------------------------

QHSC_GENES <- c(
  "Lrat",
  "Rbp1",
  "Reln",
  "Cygb",
  "Des"
)

ECM_GENES <- c(
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Col5a1",
  "Fn1",
  "Lum",
  "Dcn",
  "Bgn",
  "Timp1"
)

CONTRACTILE_GENES <- c(
  "Acta2",
  "Tagln",
  "Myl9",
  "Cnn1",
  "Tpm2",
  "Myh11"
)

# Broad HSC identity: requires multi-gene support.
HSC_CORE_GENES <- unique(
  c(
    "Lrat",
    "Rbp1",
    "Reln",
    "Cygb",
    "Des",
    "Pdgfrb",
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Lum",
    "Dcn",
    "Bgn",
    "Acta2",
    "Tagln",
    "Tpm2"
  )
)

# ------------------------------------------------------------------
# Contamination core panels
#
# A cell is counted as lineage-positive only when >=2 genes in that
# lineage are detected, reducing sensitivity to ambient RNA.
# ------------------------------------------------------------------

CONTAMINATION_PANELS <- list(

  "Endothelial" = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Vwf",
    "Erg"
  ),

  "Macrophage" = c(
    "Ptprc",
    "Adgre1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Csf1r",
    "Fcer1g",
    "Tyrobp"
  ),

  # Alb/Apoa1/Ttr are intentionally omitted here because they are highly
  # abundant liver transcripts and may appear through ambient RNA.
  "Hepatocyte" = c(
    "Hnf4a",
    "Asgr1",
    "Slc10a1",
    "Cps1",
    "Cyp2f2",
    "Ass1"
  ),

  # Krt8/Krt18 are intentionally omitted because they are too broad for
  # conservative exclusion.
  "Cholangiocyte" = c(
    "Krt19",
    "Epcam",
    "Krt7",
    "Sox9"
  ),

  "Mesothelial" = c(
    "Msln",
    "Upk3b",
    "Wt1",
    "Calb2"
  )
)

# Conservative purity thresholds.
HSC_CORE_HIT_THRESHOLD <- 1L
CONTAM_HIT_THRESHOLD <- 2L

STRICT_CONTAM_FRACTION <- 0.30
REVIEW_CONTAM_FRACTION <- 0.15

MIN_HSC_CORE_FRACTION <- 0.20

# Contractile HSC should be a clearly enriched program, not every aHSC.
CONTRACTILE_Z_THRESHOLD <- 0.75
CONTRACTILE_MARGIN <- 0.20


# ==============================================================================
# 4. Load
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input v6.1.0 HSC RDS not found: ",
    INPUT_RDS
  )
}

msg(
  "Loading v6.1.0 HSC object..."
)

hsc <- readRDS(
  INPUT_RDS
)

if (
  !"RNA" %in%
    Assays(hsc)
) {
  stop(
    "RNA assay missing."
  )
}

DefaultAssay(
  hsc
) <- "RNA"

required_meta <- c(
  CLUSTER_COL,
  SAMPLE_COL,
  CONDITION_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    hsc@meta.data
  )
)

if (
  length(
    missing_meta
  )
) {
  stop(
    "Missing metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  !HSC_UMAP %in%
    Reductions(hsc)
) {
  stop(
    "HSC UMAP reduction missing: ",
    HSC_UMAP
  )
}

# Ensure one joined RNA matrix for counts/data extraction.
if (
  length(
    Layers(
      hsc[["RNA"]]
    )
  ) > 1
) {
  hsc[["RNA"]] <- JoinLayers(
    hsc[["RNA"]]
  )
}

counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

msg(
  "Cells: ",
  ncol(hsc)
)

msg(
  "Clusters: ",
  paste(
    sort(
      unique(
        as.character(
          hsc@meta.data[[CLUSTER_COL]]
        )
      )
    ),
    collapse = ", "
  )
)


# ==============================================================================
# 5. Recompute clean HSC-state module scores
# ==============================================================================

msg(
  "Computing HSC-state module scores..."
)

hsc <- add_module_score_clean(
  hsc,
  QHSC_GENES,
  "qHSC_score_v612"
)

hsc <- add_module_score_clean(
  hsc,
  ECM_GENES,
  "ECM_score_v612"
)

hsc <- add_module_score_clean(
  hsc,
  CONTRACTILE_GENES,
  "Contractile_score_v612"
)


# ==============================================================================
# 6. Cell-level multi-gene lineage evidence
# ==============================================================================

msg(
  "Computing multi-gene lineage evidence..."
)

hsc_hits <- multi_gene_hit_count(
  counts,
  HSC_CORE_GENES
)

hsc$HSC_core_gene_hits_v612 <-
  as.integer(
    hsc_hits
  )

hsc$HSC_core_positive_v612 <-
  hsc_hits >=
    HSC_CORE_HIT_THRESHOLD

for (
  nm in names(
    CONTAMINATION_PANELS
  )
) {

  hits <- multi_gene_hit_count(
    counts,
    CONTAMINATION_PANELS[[nm]]
  )

  hit_col <- paste0(
    nm,
    "_gene_hits_v612"
  )

  positive_col <- paste0(
    nm,
    "_multiGenePositive_v612"
  )

  hsc[[hit_col]] <-
    as.integer(
      hits
    )

  hsc[[positive_col]] <-
    hits >=
      CONTAM_HIT_THRESHOLD
}


# ==============================================================================
# 7. Cluster-level purity audit
# ==============================================================================

msg(
  "Summarizing cluster purity..."
)

cluster_meta <- hsc@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    cluster = as.character(
      .data[[CLUSTER_COL]]
    )
  )

contam_positive_cols <- paste0(
  names(
    CONTAMINATION_PANELS
  ),
  "_multiGenePositive_v612"
)

purity_summary <- cluster_meta %>%
  group_by(
    cluster
  ) %>%
  summarise(
    n_cells = n(),

    HSC_core_fraction =
      mean(
        HSC_core_positive_v612,
        na.rm = TRUE
      ),

    qHSC_score =
      mean(
        qHSC_score_v612,
        na.rm = TRUE
      ),

    ECM_score =
      mean(
        ECM_score_v612,
        na.rm = TRUE
      ),

    Contractile_score =
      mean(
        Contractile_score_v612,
        na.rm = TRUE
      ),

    across(
      all_of(
        contam_positive_cols
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),

    .groups = "drop"
  )

# Rename lineage fractions for readability.
for (
  nm in names(
    CONTAMINATION_PANELS
  )
) {

  old_col <- paste0(
    nm,
    "_multiGenePositive_v612"
  )

  new_col <- paste0(
    nm,
    "_fraction"
  )

  colnames(
    purity_summary
  )[
    colnames(
      purity_summary
    ) == old_col
  ] <- new_col
}

contam_fraction_cols <- paste0(
  names(
    CONTAMINATION_PANELS
  ),
  "_fraction"
)

purity_summary$max_contamination_fraction <-
  apply(
    purity_summary[
      ,
      contam_fraction_cols,
      drop = FALSE
    ],
    1,
    max,
    na.rm = TRUE
  )

purity_summary$dominant_contamination_lineage <-
  names(
    CONTAMINATION_PANELS
  )[
    max.col(
      as.matrix(
        purity_summary[
          ,
          contam_fraction_cols,
          drop = FALSE
        ]
      ),
      ties.method = "first"
    )
  ]


# ==============================================================================
# 8. Conservative HSC purity class
# ==============================================================================

purity_summary <- purity_summary %>%
  mutate(

    contamination_to_HSC_ratio =
      max_contamination_fraction /
      pmax(
        HSC_core_fraction,
        1e-6
      ),

    purity_class = case_when(

      # Strong non-HSC lineage evidence dominates HSC identity.
      max_contamination_fraction >=
        STRICT_CONTAM_FRACTION &
        contamination_to_HSC_ratio >=
          1.0 ~
        "Excluded_nonHSC",

      # Very weak HSC identity requires manual review.
      HSC_core_fraction <
        MIN_HSC_CORE_FRACTION ~
        "Review_boundary",

      # Moderate contamination is only a boundary call when it is substantial
      # relative to the HSC signal. This avoids ambient-RNA overcalling.
      max_contamination_fraction >=
        REVIEW_CONTAM_FRACTION &
        contamination_to_HSC_ratio >=
          0.50 ~
        "Review_boundary",

      TRUE ~
        "Genuine_HSC"
    )
  )

write.csv(
  purity_summary,
  file.path(
    TAB_OUT,
    "01_HSC_cluster_purity_audit_v6.1.2.csv"
  ),
  row.names = FALSE
)

msg(
  "Purity classification:"
)

print(
  purity_summary %>%
    select(
      cluster,
      n_cells,
      HSC_core_fraction,
      max_contamination_fraction,
      contamination_to_HSC_ratio,
      dominant_contamination_lineage,
      purity_class
    )
)


# ==============================================================================
# 9. HSC-state classification among genuine clusters
# ==============================================================================

genuine_programs <- purity_summary %>%
  filter(
    purity_class ==
      "Genuine_HSC"
  ) %>%
  mutate(
    qHSC_z =
      row_z(
        qHSC_score
      ),

    ECM_z =
      row_z(
        ECM_score
      ),

    Contractile_z =
      row_z(
        Contractile_score
      )
  )

genuine_programs <- genuine_programs %>%
  mutate(

    HSC_state3 = case_when(

      Contractile_z >=
        CONTRACTILE_Z_THRESHOLD &
        Contractile_z >
          ECM_z +
            CONTRACTILE_MARGIN &
        Contractile_z >
          qHSC_z ~
        "Contractile HSC",

      ECM_z >=
        qHSC_z ~
        "ECM-activated HSC",

      TRUE ~
        "qHSC"
    ),

    HSC_state2 = ifelse(
      HSC_state3 ==
        "qHSC",
      "qHSC",
      "aHSC"
    )
  )

cluster_annotation <- purity_summary %>%
  left_join(
    genuine_programs %>%
      select(
        cluster,
        qHSC_z,
        ECM_z,
        Contractile_z,
        HSC_state3,
        HSC_state2
      ),
    by = "cluster"
  ) %>%
  mutate(

    HSC_state3_audit = case_when(
      purity_class ==
        "Excluded_nonHSC" ~
        "Excluded_nonHSC",

      purity_class ==
        "Review_boundary" ~
        "Review_boundary",

      TRUE ~
        HSC_state3
    ),

    HSC_state2_audit = case_when(
      purity_class ==
        "Excluded_nonHSC" ~
        "Excluded_nonHSC",

      purity_class ==
        "Review_boundary" ~
        "Review_boundary",

      TRUE ~
        HSC_state2
    )
  ) %>%
  arrange(
    suppressWarnings(
      as.numeric(
        cluster
      )
    )
  )

write.csv(
  cluster_annotation,
  file.path(
    TAB_OUT,
    "02_HSC_cluster_annotation_v6.1.2.csv"
  ),
  row.names = FALSE
)

msg(
  "Cluster annotation:"
)

print(
  cluster_annotation %>%
    select(
      cluster,
      n_cells,
      purity_class,
      HSC_state3_audit,
      HSC_state2_audit,
      qHSC_z,
      ECM_z,
      Contractile_z
    )
)


# ==============================================================================
# 10. Transfer cluster annotation to cells
# ==============================================================================

annotation_map3 <- setNames(
  cluster_annotation$HSC_state3_audit,
  cluster_annotation$cluster
)

annotation_map2 <- setNames(
  cluster_annotation$HSC_state2_audit,
  cluster_annotation$cluster
)

cell_cluster <- as.character(
  hsc@meta.data[[CLUSTER_COL]]
)

hsc$HSC_state3_v612 <-
  unname(
    annotation_map3[
      cell_cluster
    ]
  )

hsc$HSC_state2_v612 <-
  unname(
    annotation_map2[
      cell_cluster
    ]
  )

hsc$HSC_primary_keep_v612 <-
  hsc$HSC_state3_v612 %in%
    STATE3_LEVELS


# ==============================================================================
# 11. Save annotated audit RDS and clean HSC RDS
# ==============================================================================

ANNOTATED_RDS <- file.path(
  RDS_OUT,
  "Mouse_MASH_HSC_refined_annotated_AUDIT_v6.1.2.rds"
)

saveRDS(
  hsc,
  ANNOTATED_RDS,
  compress = FALSE
)

clean_cells <- colnames(hsc)[
  hsc$HSC_primary_keep_v612
]

if (
  length(
    clean_cells
  ) == 0
) {

  stop(
    "No Genuine_HSC cells remained after v6.1.2 purity filtering. ",
    "The purity audit CSV has been written. ",
    "Do not force annotation; inspect 01_HSC_cluster_purity_audit_v6.1.2.csv."
  )
}

hsc_clean <- subset(
  hsc,
  cells = clean_cells
)

hsc_clean$HSC_state3_v612 <- factor(
  hsc_clean$HSC_state3_v612,
  levels = STATE3_LEVELS
)

hsc_clean$HSC_state2_v612 <- factor(
  hsc_clean$HSC_state2_v612,
  levels = STATE2_LEVELS
)

hsc_clean[[SAMPLE_COL]] <- factor(
  hsc_clean@meta.data[[SAMPLE_COL]],
  levels = SAMPLE_LEVELS
)

CLEAN_RDS <- file.path(
  RDS_OUT,
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

saveRDS(
  hsc_clean,
  CLEAN_RDS,
  compress = FALSE
)

msg(
  "Annotated audit cells: ",
  ncol(hsc)
)

msg(
  "Primary clean HSC cells: ",
  ncol(hsc_clean)
)


# ==============================================================================
# 12. UMAP: purity
# ==============================================================================

PURITY_COLORS <- c(
  "Genuine_HSC" = "#00A8E8",
  "Review_boundary" = "#FFB000",
  "Excluded_nonHSC" = "#D62728"
)

cluster_to_purity <- setNames(
  cluster_annotation$purity_class,
  cluster_annotation$cluster
)

hsc$HSC_purity_v612 <- unname(
  cluster_to_purity[
    cell_cluster
  ]
)

p_purity <- DimPlot(
  hsc,
  reduction = HSC_UMAP,
  group.by = "HSC_purity_v612",
  cols = PURITY_COLORS,
  pt.size = 0.30,
  raster = FALSE
) +
  ggtitle(
    "HSC purity audit"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_purity,
  file.path(
    FIG_OUT,
    "01_HSC_purity_audit_UMAP_v6.1.2.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 13. UMAP: 3-state / 2-state
# ==============================================================================

STATE3_COLORS <- c(
  "qHSC" = "#00BFC4",
  "ECM-activated HSC" = "#FF8C00",
  "Contractile HSC" = "#E7298A"
)

STATE2_COLORS <- c(
  "qHSC" = "#00BFC4",
  "aHSC" = "#E64B35"
)

p_state3 <- DimPlot(
  hsc_clean,
  reduction = HSC_UMAP,
  group.by = "HSC_state3_v612",
  cols = STATE3_COLORS,
  pt.size = 0.32,
  raster = FALSE
) +
  ggtitle(
    "Clean HSC | 3-state annotation"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_state3,
  file.path(
    FIG_OUT,
    "02_HSC_clean_3state_UMAP_v6.1.2.pdf"
  ),
  8,
  7
)

p_state2 <- DimPlot(
  hsc_clean,
  reduction = HSC_UMAP,
  group.by = "HSC_state2_v612",
  cols = STATE2_COLORS,
  pt.size = 0.32,
  raster = FALSE
) +
  ggtitle(
    "Clean HSC | qHSC vs aHSC"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_state2,
  file.path(
    FIG_OUT,
    "03_HSC_clean_2state_UMAP_v6.1.2.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 14. UMAP by sample / Sham vs Tx
# ==============================================================================

p_state3_sample <- DimPlot(
  hsc_clean,
  reduction = HSC_UMAP,
  group.by = "HSC_state3_v612",
  split.by = SAMPLE_COL,
  cols = STATE3_COLORS,
  pt.size = 0.25,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "Clean HSC 3-state | Sham1 / Sham20 / Tx17 / Tx5"
  )

save_pdf(
  p_state3_sample,
  file.path(
    FIG_OUT,
    "04_HSC_clean_3state_UMAP_by_sample_v6.1.2.pdf"
  ),
  13,
  10
)

p_state3_condition <- DimPlot(
  hsc_clean,
  reduction = HSC_UMAP,
  group.by = "HSC_state3_v612",
  split.by = CONDITION_COL,
  cols = STATE3_COLORS,
  pt.size = 0.25,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "Clean HSC 3-state | Sham vs Tx"
  )

save_pdf(
  p_state3_condition,
  file.path(
    FIG_OUT,
    "05_HSC_clean_3state_UMAP_Sham_vs_Tx_v6.1.2.pdf"
  ),
  13,
  6
)


# ==============================================================================
# 15. Marker DotPlot for final states
# ==============================================================================

MARKER_PANEL <- unique(
  c(
    QHSC_GENES,
    ECM_GENES,
    CONTRACTILE_GENES
  )
)

MARKER_PANEL <- present_genes(
  hsc_clean,
  MARKER_PANEL
)

p_dot3 <- DotPlot(
  hsc_clean,
  features = MARKER_PANEL,
  group.by = "HSC_state3_v612",
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
      "Clean HSC 3-state marker validation",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_dot3,
  file.path(
    FIG_OUT,
    "06_HSC_clean_3state_marker_DotPlot_v6.1.2.pdf"
  ),
  15,
  5
)


# ==============================================================================
# 16. Program-score Violin
# ==============================================================================

p_vln <- VlnPlot(
  hsc_clean,
  features = c(
    "qHSC_score_v612",
    "ECM_score_v612",
    "Contractile_score_v612"
  ),
  group.by = "HSC_state3_v612",
  pt.size = 0,
  ncol = 3
) &
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_vln,
  file.path(
    FIG_OUT,
    "07_HSC_clean_3state_program_score_Violin_v6.1.2.pdf"
  ),
  13,
  4.5
)


# ==============================================================================
# 17. Cell counts / fractions by sample
# ==============================================================================

state3_by_sample <- hsc_clean@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    .data[[SAMPLE_COL]],
    .data[[CONDITION_COL]],
    HSC_state3_v612,
    name = "n_cells"
  ) %>%
  complete(
    !!sym(SAMPLE_COL) :=
      factor(
        SAMPLE_LEVELS,
        levels = SAMPLE_LEVELS
      ),
    HSC_state3_v612 :=
      factor(
        STATE3_LEVELS,
        levels = STATE3_LEVELS
      ),
    fill = list(
      n_cells = 0
    )
  ) %>%
  group_by(
    .data[[SAMPLE_COL]]
  ) %>%
  mutate(
    fraction_of_clean_HSC =
      n_cells /
      sum(
        n_cells
      )
  ) %>%
  ungroup()

write.csv(
  state3_by_sample,
  file.path(
    TAB_OUT,
    "03_HSC_3state_counts_fraction_by_sample_v6.1.2.csv"
  ),
  row.names = FALSE
)

state2_by_sample <- hsc_clean@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    .data[[SAMPLE_COL]],
    .data[[CONDITION_COL]],
    HSC_state2_v612,
    name = "n_cells"
  ) %>%
  group_by(
    .data[[SAMPLE_COL]]
  ) %>%
  mutate(
    fraction_of_clean_HSC =
      n_cells /
      sum(
        n_cells
      )
  ) %>%
  ungroup()

write.csv(
  state2_by_sample,
  file.path(
    TAB_OUT,
    "04_HSC_2state_counts_fraction_by_sample_v6.1.2.csv"
  ),
  row.names = FALSE
)

msg(
  "3-state sample counts:"
)

print(
  state3_by_sample
)


# ==============================================================================
# 18. Cell-count heatmap
# ==============================================================================

p_count_heat <- ggplot(
  state3_by_sample,
  aes(
    x = HSC_state3_v612,
    y = .data[[SAMPLE_COL]],
    fill = log10(
      n_cells + 1
    )
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = n_cells
    ),
    size = 3.5
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
      "Clean HSC cells available for interaction analysis",
    x = NULL,
    y = NULL,
    fill = "log10(n+1)"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_pdf(
  p_count_heat,
  file.path(
    FIG_OUT,
    "08_HSC_3state_cellcount_heatmap_v6.1.2.pdf"
  ),
  9,
  5
)


# ==============================================================================
# 19. Purity heatmap
# ==============================================================================

purity_long <- purity_summary %>%
  select(
    cluster,
    HSC_core_fraction,
    all_of(
      contam_fraction_cols
    )
  ) %>%
  pivot_longer(
    cols = -cluster,
    names_to = "program",
    values_to = "fraction"
  )

p_purity_heat <- ggplot(
  purity_long,
  aes(
    x = program,
    y = factor(
      cluster,
      levels = rev(
        sort(
          unique(
            suppressWarnings(
              as.numeric(
                cluster
              )
            )
          )
        )
      )
    ),
    fill = fraction
  )
) +
  geom_tile(
    linewidth = 0.3
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
      "Cluster purity audit | multi-gene co-expression",
    x = NULL,
    y = "HSC cluster",
    fill = "Fraction"
  ) +
  theme_classic(
    base_size = 8
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
  p_purity_heat,
  file.path(
    FIG_OUT,
    "09_HSC_cluster_purity_heatmap_v6.1.2.pdf"
  ),
  11,
  7
)



# ==============================================================================
# 20. Marker-panel audit
# ==============================================================================

panel_audit <- bind_rows(

  tibble(
    panel = "HSC_core",
    gene = HSC_CORE_GENES
  ),

  bind_rows(
    lapply(
      names(
        CONTAMINATION_PANELS
      ),
      function(nm) {
        tibble(
          panel = nm,
          gene = CONTAMINATION_PANELS[[nm]]
        )
      }
    )
  )
) %>%
  mutate(
    present_in_RNA =
      gene %in%
        rownames(hsc)
  )

write.csv(
  panel_audit,
  file.path(
    TAB_OUT,
    "05_HSC_purity_marker_panel_audit_v6.1.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 21. Final manifest / audit

# ==============================================================================

retention_summary <- tibble(
  metric = c(
    "input_HSC_cells",
    "primary_clean_HSC_cells",
    "primary_retention_fraction",
    "genuine_clusters",
    "review_boundary_clusters",
    "excluded_nonHSC_clusters"
  ),
  value = c(
    as.character(
      ncol(hsc)
    ),
    as.character(
      ncol(hsc_clean)
    ),
    as.character(
      ncol(hsc_clean) /
        ncol(hsc)
    ),
    paste(
      cluster_annotation$cluster[
        cluster_annotation$purity_class ==
          "Genuine_HSC"
      ],
      collapse = ","
    ),
    paste(
      cluster_annotation$cluster[
        cluster_annotation$purity_class ==
          "Review_boundary"
      ],
      collapse = ","
    ),
    paste(
      cluster_annotation$cluster[
        cluster_annotation$purity_class ==
          "Excluded_nonHSC"
      ],
      collapse = ","
    )
  )
)

write.csv(
  retention_summary,
  file.path(
    LOG_OUT,
    "retention_summary_v6.1.2.csv"
  ),
  row.names = FALSE
)

write.csv(
  tibble(
    parameter = c(
      "version",
      "input_RDS",
      "cluster_column",
      "UMAP",
      "HSC_core_hit_threshold",
      "contamination_hit_threshold",
      "strict_contam_fraction",
      "review_contam_fraction",
      "min_HSC_core_fraction",
      "contractile_z_threshold",
      "contractile_margin"
    ),
    value = c(
      "v6.1.2",
      INPUT_RDS,
      CLUSTER_COL,
      HSC_UMAP,
      HSC_CORE_HIT_THRESHOLD,
      CONTAM_HIT_THRESHOLD,
      STRICT_CONTAM_FRACTION,
      REVIEW_CONTAM_FRACTION,
      MIN_HSC_CORE_FRACTION,
      CONTRACTILE_Z_THRESHOLD,
      CONTRACTILE_MARGIN
    )
  ),
  file.path(
    LOG_OUT,
    "analysis_parameters_v6.1.2.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_OUT,
    "sessionInfo_v6.1.2.txt"
  )
)

msg(
  "DONE."
)

msg(
  "Annotated audit RDS: ",
  ANNOTATED_RDS
)

msg(
  "Clean HSC RDS: ",
  CLEAN_RDS
)

msg(
  "Output directory: ",
  OUT
)
