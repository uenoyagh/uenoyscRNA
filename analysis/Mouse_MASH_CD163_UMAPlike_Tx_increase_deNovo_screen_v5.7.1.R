#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# Cd163 UMAP-like + Tx-increase broad de novo screen
#
# Version: v5.7.1
#
# v5.7.1 CHANGE FROM v5.7.0
#   SPEED OPTIMIZATION ONLY:
#   Sham/Tx per-gene summary calculations are vectorized.
#   Analysis logic, thresholds, weights, and biological definitions are unchanged.
#
# PURPOSE
#   Search broadly across ALL RNA genes for molecules that:
#
#   A) show a spatial expression pattern on the existing macrophage UMAP
#      similar to Cd163, using STD Kupffer_Macrophage as the reference template;
#
#   B) increase after transplantation (Sham -> Tx) in STRICT-QC-clean
#      Kupffer_Macrophage cells, at both:
#         - sample-level pseudobulk expression
#         - RNA-positive cell fraction
#
#   C) show biological-replicate consistency across:
#         Sham1, Sham20, Tx17, Tx5
#
# IMPORTANT
#   - NO restriction to cell-surface proteins.
#   - Cytoplasmic, lysosomal, nuclear, metabolic, secreted, and membrane proteins
#     are all eligible.
#   - No reintegration, reclustering, or re-UMAP.
#   - Existing umapRPCA coordinates are used.
#   - No cell-level p-values are used for biological inference.
#   - n=2 Sham vs n=2 Tx: treatment statistics are descriptive/exploratory.
#
# SPATIAL SIMILARITY
#   Cd163-like UMAP pattern is evaluated in STD Kupffer_Macrophage cells using:
#     1) cell-level Spearman expression correlation
#     2) Cd163-positive binary Jaccard overlap
#     3) UMAP-grid mean-expression Spearman correlation
#     4) UMAP-grid positive-fraction Spearman correlation
#
# STRICT CONTAMINATION QC
#   Uses the v5.6.4 principle:
#   only lineage-specific CORE markers trigger Monocyte / Neutrophil / B / T / NK
#   contamination flags. Supportive markers do not trigger exclusion.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(5700)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS_CANDIDATES <- c(
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
  ),
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
  )
)

MPHI_RDS <- MPHI_RDS_CANDIDATES[
  file.exists(
    MPHI_RDS_CANDIDATES
  )
][1]

V550_RANKING <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163_Positive_deNovo_Marker_Discovery_v5.5.0",
  "Tables",
  "05_ALL_GENES_Cd163_deNovo_ranking_v5.5.0.csv"
)

OUTPUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163_UMAPlike_TxIncrease_deNovo_v5.7.1"
)

DIR_TABLE <- file.path(
  OUTPUT_DIR,
  "Tables"
)

DIR_UMAP <- file.path(
  OUTPUT_DIR,
  "UMAP"
)

DIR_PLOT <- file.path(
  OUTPUT_DIR,
  "Plots"
)

DIR_LOG <- file.path(
  OUTPUT_DIR,
  "Logs"
)

for (d in c(
  OUTPUT_DIR,
  DIR_TABLE,
  DIR_UMAP,
  DIR_PLOT,
  DIR_LOG
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. Settings
# ==============================================================================

ASSAY_USE <- "RNA"

REFERENCE_GENE <- "Cd163"

DISCOVERY_CONDITION <- "STD"

TREATMENT_GROUPS <- c(
  "Sham",
  "Tx"
)

EXPECTED_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

MACROPHAGE_LINEAGE_VALUE <- "Kupffer_Macrophage"

LINEAGE_COLUMN_CANDIDATES <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_v440",
  "layer1_original",
  "celltype_for_R8plot",
  "celltype_auto_annotation"
)

SAMPLE_COLUMN_CANDIDATES <- c(
  "sample_for_annotation",
  "sample",
  "orig.ident"
)

CONDITION_COLUMN_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group",
  "sample_for_annotation",
  "sample"
)

UMAP_CANDIDATES <- c(
  "umapRPCA",
  "umap"
)

# UMAP grid for spatial-pattern comparison.
GRID_NX <- 30
GRID_NY <- 30
MIN_CELLS_PER_GRID_BIN <- 3

# Expression filter.
# Candidate must be detected in at least 1% of either:
#   - STD cells
#   - Sham cells
#   - Tx cells
MIN_DETECTION_FRACTION <- 0.01

# Ranking weights.
W_UMAP_SPATIAL <- 0.50
W_TX_PSEUDOBULK_INCREASE <- 0.20
W_TX_POSITIVE_FRACTION_INCREASE <- 0.15
W_REPLICATE_CONSISTENCY <- 0.15

# Within UMAP-spatial score.
W_GRID_MEAN <- 0.35
W_GRID_PCT <- 0.30
W_CELL_SPEARMAN <- 0.20
W_BINARY_JACCARD <- 0.15

# Stringent shortlist thresholds.
MIN_SPATIAL_SCORE <- 0.55
MIN_SHAM_PCT_POSITIVE <- 2
MIN_TX_PCT_POSITIVE <- 5
MIN_LOG2FC_TX_VS_SHAM <- 0.25
MIN_DELTA_PCT_TX_MINUS_SHAM <- 3
MIN_PAIRWISE_INCREASE_FRACTION <- 0.75

TOP_N_TABLE <- 200
TOP_N_RANK_PLOT <- 40
TOP_N_FEATURE <- 24

# Technical genes remain in the ALL-gene raw table but are excluded from the
# biological shortlist/top plots.
TECHNICAL_GENE_PATTERN <- "^(mt-|Rpl|Rps|Hba-|Hbb-)"

# Previously interesting genes are annotations only; they DO NOT affect ranking.
REFERENCE_PANEL <- c(
  "Cd163",
  "Timd4",
  "C6",
  "Clec4f",
  "Fcna",
  "Folr2",
  "Stab2",
  "Cdh5",
  "Cd5l",
  "Vsig4",
  "Marco",
  "Apoc1",
  "Slc40a1",
  "Mrc1",
  "Mertk",
  "Hmox1",
  "Arg1",
  "Siglece",
  "Ms4a7",
  "Hexb",
  "Alox5ap",
  "Ltc4s",
  "Ntpcr",
  "Aif1"
)

# ==============================================================================
# 3. Strict contamination marker sets from v5.6.4
# ==============================================================================

CONTAMINATION_CORE_MARKERS <- list(

  Monocyte = c(
    "Ccr2",
    "Ly6c2",
    "Plac8",
    "Sell",
    "Fcgr3",
    "Ms4a4c"
  ),

  Neutrophil = c(
    "Ly6g",
    "Csf3r",
    "Retnlg",
    "Camp",
    "Ngp",
    "Mpo",
    "Elane"
  ),

  B = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd19",
    "Ebf1",
    "Pax5",
    "Cd22"
  ),

  T = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trac",
    "Lck",
    "Lat",
    "Cd247",
    "Il7r"
  ),

  NK = c(
    "Ncr1",
    "Klrd1",
    "Klrk1",
    "Prf1",
    "Gzmb",
    "Ccl5"
  )
)

PRIMARY_MIN_CORE_MARKERS <- 2

# ==============================================================================
# 4. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(
  plot,
  filename,
  width,
  height
) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
}

resolve_first <- function(
  available,
  candidates
) {

  hit <- candidates[
    candidates %in%
      available
  ]

  if (
    !length(
      hit
    )
  ) {
    return(
      NA_character_
    )
  }

  hit[[1]]
}

canonical_condition <- function(x) {

  x <- as.character(x)

  out <- case_when(
    grepl(
      "^STD",
      x,
      ignore.case = TRUE
    ) ~ "STD",

    grepl(
      "CDHFD|CDAHFD",
      x,
      ignore.case = TRUE
    ) ~ "CDAHFD",

    grepl(
      "^Sham",
      x,
      ignore.case = TRUE
    ) ~ "Sham",

    grepl(
      "^Tx",
      x,
      ignore.case = TRUE
    ) ~ "Tx",

    TRUE ~ NA_character_
  )

  factor(
    out,
    levels = c(
      "STD",
      "CDAHFD",
      "Sham",
      "Tx"
    )
  )
}

resolve_condition <- function(obj) {

  meta_cols <- colnames(
    obj@meta.data
  )

  for (
    src in CONDITION_COLUMN_CANDIDATES
  ) {

    if (
      !src %in%
        meta_cols
    ) {
      next
    }

    cond <- canonical_condition(
      obj@meta.data[[src]]
    )

    if (
      !all(
        is.na(
          cond
        )
      )
    ) {

      return(
        list(
          source = src,
          condition = cond
        )
      )
    }
  }

  stop(
    "No usable condition metadata column found."
  )
}

safe_spearman <- function(
  x,
  y
) {

  ok <- is.finite(x) &
    is.finite(y)

  x <- x[ok]
  y <- y[ok]

  if (
    length(x) < 3 ||
    length(
      unique(x)
    ) < 2 ||
    length(
      unique(y)
    ) < 2
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )
}

safe_rescale01 <- function(x) {

  out <- rep(
    NA_real_,
    length(x)
  )

  ok <- is.finite(
    x
  )

  if (
    !any(
      ok
    )
  ) {
    return(
      out
    )
  }

  r <- range(
    x[ok],
    na.rm = TRUE
  )

  if (
    !all(
      is.finite(
        r
      )
    ) ||
    diff(
      r
    ) == 0
  ) {

    out[ok] <- 0.5

    return(
      out
    )
  }

  out[ok] <- scales::rescale(
    x[ok],
    to = c(
      0,
      1
    )
  )

  out
}

binary_jaccard <- function(
  reference_positive,
  candidate_positive
) {

  inter <- sum(
    reference_positive &
      candidate_positive
  )

  union <- sum(
    reference_positive |
      candidate_positive
  )

  if (
    union == 0
  ) {
    return(
      NA_real_
    )
  }

  inter /
    union
}

pairwise_increase_fraction <- function(
  sham_values,
  tx_values
) {

  sham_values <- sham_values[
    is.finite(
      sham_values
    )
  ]

  tx_values <- tx_values[
    is.finite(
      tx_values
    )
  ]

  if (
    !length(
      sham_values
    ) ||
    !length(
      tx_values
    )
  ) {
    return(
      NA_real_
    )
  }

  mean(
    outer(
      tx_values,
      sham_values,
      FUN = ">"
    )
  )
}

# ==============================================================================
# 5. Validate and load
# ==============================================================================

if (
  length(
    MPHI_RDS
  ) == 0L ||
  is.na(
    MPHI_RDS
  ) ||
  !file.exists(
    MPHI_RDS
  )
) {
  stop(
    "Clean-B macrophage RDS not found."
  )
}

msg(
  "Loading Clean-B macrophage RDS..."
)

mphi <- readRDS(
  MPHI_RDS
)

if (
  !ASSAY_USE %in%
    Assays(
      mphi
    )
) {
  stop(
    "RNA assay missing."
  )
}

DefaultAssay(
  mphi
) <- ASSAY_USE

if (
  !REFERENCE_GENE %in%
    rownames(
      mphi
    )
) {
  stop(
    "Cd163 not found in RNA assay."
  )
}

# ==============================================================================
# 6. Resolve metadata
# ==============================================================================

condition_info <- resolve_condition(
  mphi
)

mphi$condition_v571 <-
  condition_info$condition

lineage_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  LINEAGE_COLUMN_CANDIDATES
)

sample_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  SAMPLE_COLUMN_CANDIDATES
)

umap_name <- resolve_first(
  Reductions(
    mphi
  ),
  UMAP_CANDIDATES
)

if (
  is.na(
    lineage_col
  )
) {
  stop(
    "Could not resolve lineage column."
  )
}

if (
  is.na(
    sample_col
  )
) {
  stop(
    "Could not resolve sample column."
  )
}

if (
  is.na(
    umap_name
  )
) {
  stop(
    "Could not resolve UMAP reduction."
  )
}

msg(
  "Condition source: ",
  condition_info$source
)

msg(
  "Lineage column: ",
  lineage_col
)

msg(
  "Sample column: ",
  sample_col
)

msg(
  "UMAP: ",
  umap_name
)

# ==============================================================================
# 7. Extract all Kupffer_Macrophage cells
# ==============================================================================

lineage <- as.character(
  mphi@meta.data[[lineage_col]]
)

kupffer_cells <- colnames(
  mphi
)[
  lineage ==
    MACROPHAGE_LINEAGE_VALUE
]

if (
  length(
    kupffer_cells
  ) < 100
) {
  stop(
    "Too few Kupffer_Macrophage cells."
  )
}

kupffer <- subset(
  mphi,
  cells = kupffer_cells
)

DefaultAssay(
  kupffer
) <- ASSAY_USE

kupffer$sample_v571 <- as.character(
  kupffer@meta.data[[sample_col]]
)

kupffer$condition_v571 <- as.character(
  kupffer$condition_v571
)

msg(
  "Kupffer_Macrophage cells: ",
  ncol(
    kupffer
  )
)

# ==============================================================================
# 8. RNA matrices
# ==============================================================================

counts_all <- GetAssayData(
  kupffer,
  assay = ASSAY_USE,
  layer = "counts"
)

data_all <- GetAssayData(
  kupffer,
  assay = ASSAY_USE,
  layer = "data"
)

all_genes <- rownames(
  counts_all
)

# ==============================================================================
# 9. STD reference cells for Cd163 UMAP template
# ==============================================================================

std_idx <- which(
  kupffer$condition_v571 ==
    DISCOVERY_CONDITION
)

if (
  length(
    std_idx
  ) < 50
) {
  stop(
    "Too few STD Kupffer cells."
  )
}

std_cells <- colnames(
  kupffer
)[
  std_idx
]

std_counts <- counts_all[
  ,
  std_idx,
  drop = FALSE
]

std_data <- data_all[
  ,
  std_idx,
  drop = FALSE
]

std_umap <- Embeddings(
  kupffer,
  reduction = umap_name
)[
  std_cells,
  ,
  drop = FALSE
]

msg(
  "STD reference cells: ",
  length(
    std_cells
  )
)

cd163_std_expr <- as.numeric(
  std_data[
    REFERENCE_GENE,
    ,
    drop = TRUE
  ]
)

cd163_std_positive <- as.numeric(
  std_counts[
    REFERENCE_GENE,
    ,
    drop = TRUE
  ]
) > 0

msg(
  "STD Cd163-positive cells: ",
  sum(
    cd163_std_positive
  ),
  " / ",
  length(
    cd163_std_positive
  ),
  " (",
  round(
    100 *
      mean(
        cd163_std_positive
      ),
    2
  ),
  "%)"
)

# ==============================================================================
# 10. Build fixed UMAP grid
# ==============================================================================

x <- std_umap[
  ,
  1
]

y <- std_umap[
  ,
  2
]

x_breaks <- seq(
  min(
    x
  ) -
    1e-8,
  max(
    x
  ) +
    1e-8,
  length.out =
    GRID_NX +
    1
)

y_breaks <- seq(
  min(
    y
  ) -
    1e-8,
  max(
    y
  ) +
    1e-8,
  length.out =
    GRID_NY +
    1
)

x_bin <- cut(
  x,
  breaks = x_breaks,
  include.lowest = TRUE,
  labels = FALSE
)

y_bin <- cut(
  y,
  breaks = y_breaks,
  include.lowest = TRUE,
  labels = FALSE
)

grid_id_raw <- (
  y_bin -
    1L
) *
  GRID_NX +
  x_bin

occupied_grid <- sort(
  unique(
    grid_id_raw
  )
)

grid_map <- match(
  grid_id_raw,
  occupied_grid
)

n_bins <- length(
  occupied_grid
)

membership <- Matrix::sparseMatrix(
  i = seq_along(
    grid_map
  ),
  j = grid_map,
  x = 1,
  dims = c(
    length(
      grid_map
    ),
    n_bins
  )
)

bin_n <- as.numeric(
  Matrix::colSums(
    membership
  )
)

keep_bins <- which(
  bin_n >=
    MIN_CELLS_PER_GRID_BIN
)

if (
  length(
    keep_bins
  ) < 10
) {
  stop(
    "Too few occupied UMAP bins after MIN_CELLS_PER_GRID_BIN filter."
  )
}

msg(
  "UMAP grid bins retained: ",
  length(
    keep_bins
  )
)

# Gene x grid-bin mean log-normalized expression.
grid_sum_expr <- std_data %*%
  membership

grid_mean_expr <- sweep(
  as.matrix(
    grid_sum_expr[
      ,
      keep_bins,
      drop = FALSE
    ]
  ),
  2,
  bin_n[
    keep_bins
  ],
  FUN = "/"
)

rm(
  grid_sum_expr
)

gc(
  verbose = FALSE
)

# Gene x grid-bin RNA-positive fraction.
std_binary <- std_counts > 0

grid_sum_positive <- std_binary %*%
  membership

grid_pct_positive <- sweep(
  as.matrix(
    grid_sum_positive[
      ,
      keep_bins,
      drop = FALSE
    ]
  ),
  2,
  bin_n[
    keep_bins
  ],
  FUN = "/"
)

rm(
  grid_sum_positive
)

gc(
  verbose = FALSE
)

cd163_grid_mean <- as.numeric(
  grid_mean_expr[
    REFERENCE_GENE,
    ,
    drop = TRUE
  ]
)

cd163_grid_pct <- as.numeric(
  grid_pct_positive[
    REFERENCE_GENE,
    ,
    drop = TRUE
  ]
)

# ==============================================================================
# 11. Expression filter for broad all-gene search
# ==============================================================================

sham_idx_all <- which(
  kupffer$condition_v571 ==
    "Sham"
)

tx_idx_all <- which(
  kupffer$condition_v571 ==
    "Tx"
)

pct_std <- as.numeric(
  Matrix::rowMeans(
    std_counts > 0
  )
)

pct_sham_raw <- as.numeric(
  Matrix::rowMeans(
    counts_all[
      ,
      sham_idx_all,
      drop = FALSE
    ] > 0
  )
)

pct_tx_raw <- as.numeric(
  Matrix::rowMeans(
    counts_all[
      ,
      tx_idx_all,
      drop = FALSE
    ] > 0
  )
)

eligible <- pmax(
  pct_std,
  pct_sham_raw,
  pct_tx_raw
) >=
  MIN_DETECTION_FRACTION

genes_eval <- all_genes[
  eligible
]

genes_eval <- setdiff(
  genes_eval,
  REFERENCE_GENE
)

msg(
  "Genes evaluated: ",
  length(
    genes_eval
  ),
  " / ",
  length(
    all_genes
  )
)

# ==============================================================================
# 12. Direct Cd163 UMAP-spatial similarity
# ==============================================================================

spatial_list <- vector(
  "list",
  length(
    genes_eval
  )
)

names(
  spatial_list
) <- genes_eval

for (
  i in seq_along(
    genes_eval
  )
) {

  gene <- genes_eval[[i]]

  gidx <- match(
    gene,
    all_genes
  )

  expr <- as.numeric(
    std_data[
      gene,
      ,
      drop = TRUE
    ]
  )

  detect <- as.numeric(
    std_counts[
      gene,
      ,
      drop = TRUE
    ]
  ) > 0

  grid_mean_gene <- as.numeric(
    grid_mean_expr[
      gene,
      ,
      drop = TRUE
    ]
  )

  grid_pct_gene <- as.numeric(
    grid_pct_positive[
      gene,
      ,
      drop = TRUE
    ]
  )

  spatial_list[[i]] <- tibble(

    gene = gene,

    STD_pct_positive =
      100 *
      pct_std[
        gidx
      ],

    cell_spearman_with_Cd163 =
      safe_spearman(
        cd163_std_expr,
        expr
      ),

    binary_jaccard_with_Cd163 =
      binary_jaccard(
        cd163_std_positive,
        detect
      ),

    UMAP_grid_mean_spearman =
      safe_spearman(
        cd163_grid_mean,
        grid_mean_gene
      ),

    UMAP_grid_pct_spearman =
      safe_spearman(
        cd163_grid_pct,
        grid_pct_gene
      )
  )
}

spatial_metrics <- bind_rows(
  spatial_list
)

spatial_metrics <- spatial_metrics %>%
  mutate(

    score_grid_mean =
      safe_rescale01(
        pmax(
          UMAP_grid_mean_spearman,
          0
        )
      ),

    score_grid_pct =
      safe_rescale01(
        pmax(
          UMAP_grid_pct_spearman,
          0
        )
      ),

    score_cell_spearman =
      safe_rescale01(
        pmax(
          cell_spearman_with_Cd163,
          0
        )
      ),

    score_binary_jaccard =
      safe_rescale01(
        binary_jaccard_with_Cd163
      ),

    Cd163_UMAP_spatial_score =
      W_GRID_MEAN *
        score_grid_mean +
      W_GRID_PCT *
        score_grid_pct +
      W_CELL_SPEARMAN *
        score_cell_spearman +
      W_BINARY_JACCARD *
        score_binary_jaccard
  ) %>%
  arrange(
    desc(
      Cd163_UMAP_spatial_score
    )
  ) %>%
  mutate(
    spatial_rank =
      row_number()
  )

write.csv(
  spatial_metrics,
  file.path(
    DIR_TABLE,
    "01_ALL_GENES_Cd163_UMAP_spatial_similarity_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 13. Strict contamination QC on Sham/Tx cells
# ==============================================================================

treat_idx <- which(
  kupffer$condition_v571 %in%
    TREATMENT_GROUPS
)

treat <- subset(
  kupffer,
  cells = colnames(
    kupffer
  )[
    treat_idx
  ]
)

DefaultAssay(
  treat
) <- ASSAY_USE

treat_counts <- GetAssayData(
  treat,
  assay = ASSAY_USE,
  layer = "counts"
)

treat_data <- GetAssayData(
  treat,
  assay = ASSAY_USE,
  layer = "data"
)

core_use <- lapply(
  CONTAMINATION_CORE_MARKERS,
  function(x) {
    intersect(
      x,
      rownames(
        treat
      )
    )
  }
)

flag_matrix <- matrix(
  FALSE,
  nrow = ncol(
    treat
  ),
  ncol = length(
    core_use
  ),
  dimnames = list(
    colnames(
      treat
    ),
    names(
      core_use
    )
  )
)

for (
  lineage_name in names(
    core_use
  )
) {

  genes <- core_use[[
    lineage_name]
  ]

  if (
    !length(
      genes
    )
  ) {
    next
  }

  marker_count <- as.numeric(
    Matrix::colSums(
      treat_counts[
        genes,
        ,
        drop = FALSE
      ] > 0
    )
  )

  flag_matrix[
    ,
    lineage_name
  ] <- marker_count >=
    PRIMARY_MIN_CORE_MARKERS
}

strict_flag <- rowSums(
  flag_matrix
) > 0

strict_clean_cells <- colnames(
  treat
)[
  !strict_flag
]

clean <- subset(
  treat,
  cells =
    strict_clean_cells
)

DefaultAssay(
  clean
) <- ASSAY_USE

clean$sample_v571 <- as.character(
  clean@meta.data[[sample_col]]
)

clean$condition_v571 <- as.character(
  clean$condition_v571
)

msg(
  "Sham/Tx cells before strict QC: ",
  ncol(
    treat
  )
)

msg(
  "Sham/Tx cells after strict QC: ",
  ncol(
    clean
  ),
  " (removed ",
  sum(
    strict_flag
  ),
  ")"
)

strict_qc_summary <- tibble(
  metric = c(
    "cells_before_QC",
    "cells_flagged",
    "cells_after_QC",
    "pct_removed"
  ),
  value = c(
    ncol(
      treat
    ),
    sum(
      strict_flag
    ),
    ncol(
      clean
    ),
    100 *
      mean(
        strict_flag
      )
  )
)

write.csv(
  strict_qc_summary,
  file.path(
    DIR_TABLE,
    "02_STRICT_QC_summary_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 14. Sample-level pseudobulk and positive fractions after strict QC
# ==============================================================================

clean_counts <- GetAssayData(
  clean,
  assay = ASSAY_USE,
  layer = "counts"
)

clean_data <- GetAssayData(
  clean,
  assay = ASSAY_USE,
  layer = "data"
)

genes_treat <- intersect(
  genes_eval,
  rownames(
    clean
  )
)

# Include Cd163 as a reference row in treatment tables.
genes_treat_with_ref <- unique(
  c(
    REFERENCE_GENE,
    genes_treat
  )
)

samples_present <- sort(
  unique(
    clean$sample_v571
  )
)

sample_condition_map <- tibble(
  sample = clean$sample_v571,
  condition = clean$condition_v571
) %>%
  distinct()

sample_counts <- tibble(
  sample = clean$sample_v571,
  condition = clean$condition_v571
) %>%
  count(
    condition,
    sample,
    name = "n_cells"
  ) %>%
  arrange(
    condition,
    sample
  )

write.csv(
  sample_counts,
  file.path(
    DIR_TABLE,
    "03_STRICT_QC_sample_cell_counts_v5.7.1.csv"
  ),
  row.names = FALSE
)

pb_counts <- matrix(
  0,
  nrow = length(
    genes_treat_with_ref
  ),
  ncol = length(
    samples_present
  ),
  dimnames = list(
    genes_treat_with_ref,
    samples_present
  )
)

sample_pct <- matrix(
  NA_real_,
  nrow = length(
    genes_treat_with_ref
  ),
  ncol = length(
    samples_present
  ),
  dimnames = list(
    genes_treat_with_ref,
    samples_present
  )
)

sample_mean_data <- matrix(
  NA_real_,
  nrow = length(
    genes_treat_with_ref
  ),
  ncol = length(
    samples_present
  ),
  dimnames = list(
    genes_treat_with_ref,
    samples_present
  )
)

sample_ncells <- setNames(
  integer(
    length(
      samples_present
    )
  ),
  samples_present
)

for (
  s in samples_present
) {

  idx <- which(
    clean$sample_v571 ==
      s
  )

  sample_ncells[[s]] <-
    length(
      idx
    )

  pb_counts[
    ,
    s
  ] <- as.numeric(
    Matrix::rowSums(
      clean_counts[
        genes_treat_with_ref,
        idx,
        drop = FALSE
      ]
    )
  )

  sample_pct[
    ,
    s
  ] <- 100 *
    as.numeric(
      Matrix::rowMeans(
        clean_counts[
          genes_treat_with_ref,
          idx,
          drop = FALSE
        ] > 0
      )
    )

  sample_mean_data[
    ,
    s
  ] <- as.numeric(
    Matrix::rowMeans(
      clean_data[
        genes_treat_with_ref,
        idx,
        drop = FALSE
      ]
    )
  )
}

library_size <- colSums(
  pb_counts
)

pb_cpm <- sweep(
  pb_counts,
  2,
  library_size,
  FUN = "/"
) *
  1e6

pb_log2cpm <- log2(
  pb_cpm +
    1
)

sample_long <- bind_rows(
  lapply(
    samples_present,
    function(s) {

      tibble(
        gene =
          genes_treat_with_ref,

        sample =
          s,

        n_cells =
          sample_ncells[[s]],

        pseudobulk_CPM =
          pb_cpm[
            genes_treat_with_ref,
            s
          ],

        pseudobulk_log2CPM =
          pb_log2cpm[
            genes_treat_with_ref,
            s
          ],

        pct_positive =
          sample_pct[
            genes_treat_with_ref,
            s
          ],

        mean_log_normalized_expression =
          sample_mean_data[
            genes_treat_with_ref,
            s
          ]
      )
    }
  )
) %>%
  left_join(
    sample_condition_map,
    by = "sample"
  ) %>%
  select(
    gene,
    condition,
    sample,
    n_cells,
    pseudobulk_CPM,
    pseudobulk_log2CPM,
    pct_positive,
    mean_log_normalized_expression
  )

write.csv(
  sample_long,
  file.path(
    DIR_TABLE,
    "04_ALL_GENES_samplelevel_metrics_after_STRICT_QC_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 15. Sham -> Tx increase metrics
# ==============================================================================
#
# v5.7.1 SPEED OPTIMIZATION ONLY
#
# v5.7.0 calculated these values by filtering sample_long separately for every
# gene (~12,000 repeated dplyr scans).  v5.7.1 performs the identical
# calculations directly on the already-created gene x sample matrices.
#
# No biological definition, threshold, weight, or scoring rule is changed.
# ==============================================================================

msg(
  "Calculating Sham -> Tx metrics with vectorized matrix operations..."
)

sham_samples <- sample_condition_map %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  pull(
    sample
  )

tx_samples <- sample_condition_map %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  pull(
    sample
  )

sham_samples <- intersect(
  sham_samples,
  colnames(
    pb_cpm
  )
)

tx_samples <- intersect(
  tx_samples,
  colnames(
    pb_cpm
  )
)

if (
  length(
    sham_samples
  ) == 0 ||
  length(
    tx_samples
  ) == 0
) {
  stop(
    "Sham or Tx sample columns are missing from pseudobulk matrices."
  )
}

msg(
  "Sham samples: ",
  paste(
    sham_samples,
    collapse = ", "
  )
)

msg(
  "Tx samples: ",
  paste(
    tx_samples,
    collapse = ", "
  )
)

# ------------------------------------------------------------------
# Group means
# ------------------------------------------------------------------

Sham_mean_CPM_v571 <- rowMeans(
  pb_cpm[
    genes_treat_with_ref,
    sham_samples,
    drop = FALSE
  ],
  na.rm = TRUE
)

Tx_mean_CPM_v571 <- rowMeans(
  pb_cpm[
    genes_treat_with_ref,
    tx_samples,
    drop = FALSE
  ],
  na.rm = TRUE
)

pseudobulk_log2FC_v571 <- log2(
  (
    Tx_mean_CPM_v571 +
      0.5
  ) /
    (
      Sham_mean_CPM_v571 +
        0.5
    )
)

Sham_mean_pct_v571 <- rowMeans(
  sample_pct[
    genes_treat_with_ref,
    sham_samples,
    drop = FALSE
  ],
  na.rm = TRUE
)

Tx_mean_pct_v571 <- rowMeans(
  sample_pct[
    genes_treat_with_ref,
    tx_samples,
    drop = FALSE
  ],
  na.rm = TRUE
)

delta_pct_v571 <-
  Tx_mean_pct_v571 -
  Sham_mean_pct_v571

# ------------------------------------------------------------------
# Pairwise replicate consistency
#
# Same definition as v5.7.0:
#   fraction of all Tx x Sham pairs in which Tx > Sham.
#
# This loops over SAMPLE PAIRS only, never over genes.
# For 2 Sham x 2 Tx this is only four vectorized comparisons.
# ------------------------------------------------------------------

n_gene_v571 <- length(
  genes_treat_with_ref
)

pb_pairwise_sum_v571 <- rep(
  0,
  n_gene_v571
)

pct_pairwise_sum_v571 <- rep(
  0,
  n_gene_v571
)

pb_pairwise_n_v571 <- rep(
  0,
  n_gene_v571
)

pct_pairwise_n_v571 <- rep(
  0,
  n_gene_v571
)

for (
  tx_sample in tx_samples
) {

  for (
    sham_sample in sham_samples
  ) {

    tx_pb <- pb_log2cpm[
      genes_treat_with_ref,
      tx_sample
    ]

    sham_pb <- pb_log2cpm[
      genes_treat_with_ref,
      sham_sample
    ]

    ok_pb <- is.finite(
      tx_pb
    ) &
      is.finite(
        sham_pb
      )

    pb_pairwise_sum_v571[
      ok_pb
    ] <- pb_pairwise_sum_v571[
      ok_pb
    ] +
      as.numeric(
        tx_pb[
          ok_pb
        ] >
          sham_pb[
            ok_pb
          ]
      )

    pb_pairwise_n_v571[
      ok_pb
    ] <- pb_pairwise_n_v571[
      ok_pb
    ] +
      1

    tx_pct <- sample_pct[
      genes_treat_with_ref,
      tx_sample
    ]

    sham_pct <- sample_pct[
      genes_treat_with_ref,
      sham_sample
    ]

    ok_pct <- is.finite(
      tx_pct
    ) &
      is.finite(
        sham_pct
      )

    pct_pairwise_sum_v571[
      ok_pct
    ] <- pct_pairwise_sum_v571[
      ok_pct
    ] +
      as.numeric(
        tx_pct[
          ok_pct
        ] >
          sham_pct[
            ok_pct
          ]
      )

    pct_pairwise_n_v571[
      ok_pct
    ] <- pct_pairwise_n_v571[
      ok_pct
    ] +
      1
  }
}

pseudobulk_pairwise_fraction_v571 <- ifelse(
  pb_pairwise_n_v571 > 0,
  pb_pairwise_sum_v571 /
    pb_pairwise_n_v571,
  NA_real_
)

pct_pairwise_fraction_v571 <- ifelse(
  pct_pairwise_n_v571 > 0,
  pct_pairwise_sum_v571 /
    pct_pairwise_n_v571,
  NA_real_
)

# ------------------------------------------------------------------
# "All Tx above all Sham"
#
# Same condition as v5.7.0:
#   min(Tx) > max(Sham)
# ------------------------------------------------------------------

tx_pb_matrix_v571 <- pb_log2cpm[
  genes_treat_with_ref,
  tx_samples,
  drop = FALSE
]

sham_pb_matrix_v571 <- pb_log2cpm[
  genes_treat_with_ref,
  sham_samples,
  drop = FALSE
]

tx_pct_matrix_v571 <- sample_pct[
  genes_treat_with_ref,
  tx_samples,
  drop = FALSE
]

sham_pct_matrix_v571 <- sample_pct[
  genes_treat_with_ref,
  sham_samples,
  drop = FALSE
]

min_tx_pb_v571 <- apply(
  tx_pb_matrix_v571,
  1,
  min,
  na.rm = TRUE
)

max_sham_pb_v571 <- apply(
  sham_pb_matrix_v571,
  1,
  max,
  na.rm = TRUE
)

min_tx_pct_v571 <- apply(
  tx_pct_matrix_v571,
  1,
  min,
  na.rm = TRUE
)

max_sham_pct_v571 <- apply(
  sham_pct_matrix_v571,
  1,
  max,
  na.rm = TRUE
)

pb_all_tx_above_v571 <-
  min_tx_pb_v571 >
  max_sham_pb_v571

pct_all_tx_above_v571 <-
  min_tx_pct_v571 >
  max_sham_pct_v571

# ------------------------------------------------------------------
# Build the same treatment_metrics table as v5.7.0
# ------------------------------------------------------------------

treatment_metrics <- tibble(

  gene =
    genes_treat_with_ref,

  Sham_mean_CPM =
    as.numeric(
      Sham_mean_CPM_v571
    ),

  Tx_mean_CPM =
    as.numeric(
      Tx_mean_CPM_v571
    ),

  pseudobulk_log2FC_Tx_vs_Sham =
    as.numeric(
      pseudobulk_log2FC_v571
    ),

  Sham_mean_pct_positive =
    as.numeric(
      Sham_mean_pct_v571
    ),

  Tx_mean_pct_positive =
    as.numeric(
      Tx_mean_pct_v571
    ),

  delta_pct_Tx_minus_Sham =
    as.numeric(
      delta_pct_v571
    ),

  pseudobulk_pairwise_increase_fraction =
    as.numeric(
      pseudobulk_pairwise_fraction_v571
    ),

  pct_positive_pairwise_increase_fraction =
    as.numeric(
      pct_pairwise_fraction_v571
    ),

  pseudobulk_all_Tx_above_all_Sham =
    as.logical(
      pb_all_tx_above_v571
    ),

  pct_positive_all_Tx_above_all_Sham =
    as.logical(
      pct_all_tx_above_v571
    )
)

# ------------------------------------------------------------------
# Add explicit sample columns.
#
# Kept equivalent to v5.7.0 output schema.
# ------------------------------------------------------------------

pb_wide <- sample_long %>%
  select(
    gene,
    sample,
    pseudobulk_log2CPM
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = pseudobulk_log2CPM,
    names_prefix = "log2CPM_"
  )

pct_wide <- sample_long %>%
  select(
    gene,
    sample,
    pct_positive
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = pct_positive,
    names_prefix = "pctpos_"
  )

treatment_metrics <- treatment_metrics %>%
  left_join(
    pb_wide,
    by = "gene"
  ) %>%
  left_join(
    pct_wide,
    by = "gene"
  )

write.csv(
  treatment_metrics,
  file.path(
    DIR_TABLE,
    "05_ALL_GENES_Sham_to_Tx_increase_after_STRICT_QC_v5.7.1.csv"
  ),
  row.names = FALSE
)

msg(
  "Vectorized Sham -> Tx metrics completed for ",
  nrow(
    treatment_metrics
  ),
  " genes."
)

# ==============================================================================
# 16. Integrate Cd163 UMAP pattern + Tx increase
# ==============================================================================

ranking <- spatial_metrics %>%
  inner_join(
    treatment_metrics,
    by = "gene"
  ) %>%
  mutate(

    technical_gene =
      grepl(
        TECHNICAL_GENE_PATTERN,
        gene
      ),

    in_reference_panel =
      gene %in%
      REFERENCE_PANEL,

    score_UMAP_spatial =
      safe_rescale01(
        Cd163_UMAP_spatial_score
      ),

    score_Tx_pseudobulk_increase =
      safe_rescale01(
        pmax(
          pseudobulk_log2FC_Tx_vs_Sham,
          0
        )
      ),

    score_Tx_positive_fraction_increase =
      safe_rescale01(
        pmax(
          delta_pct_Tx_minus_Sham,
          0
        )
      ),

    replicate_consistency =
      rowMeans(
        cbind(
          pseudobulk_pairwise_increase_fraction,
          pct_positive_pairwise_increase_fraction
        ),
        na.rm = TRUE
      ),

    score_replicate_consistency =
      replicate_consistency,

    Cd163like_TxIncrease_score =
      W_UMAP_SPATIAL *
        score_UMAP_spatial +
      W_TX_PSEUDOBULK_INCREASE *
        score_Tx_pseudobulk_increase +
      W_TX_POSITIVE_FRACTION_INCREASE *
        score_Tx_positive_fraction_increase +
      W_REPLICATE_CONSISTENCY *
        score_replicate_consistency
  ) %>%
  arrange(
    desc(
      Cd163like_TxIncrease_score
    )
  ) %>%
  mutate(
    final_rank =
      row_number()
  )

# Optional annotation from v5.5.0, if present.
if (
  file.exists(
    V550_RANKING
  )
) {

  v550 <- read.csv(
    V550_RANKING,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  keep_cols <- intersect(
    c(
      "gene",
      "de_novo_rank",
      "de_novo_score",
      "auc_Cd163pos_vs_neg",
      "jaccard_with_Cd163"
    ),
    colnames(
      v550
    )
  )

  v550_small <- v550[
    ,
    keep_cols,
    drop = FALSE
  ]

  ranking <- ranking %>%
    left_join(
      v550_small,
      by = "gene"
    )
}

write.csv(
  ranking,
  file.path(
    DIR_TABLE,
    "06_FINAL_ALL_GENES_Cd163_UMAPlike_TxIncrease_ranking_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 17. Broad biological shortlist
# ==============================================================================

shortlist <- ranking %>%
  filter(
    !technical_gene,

    Cd163_UMAP_spatial_score >=
      MIN_SPATIAL_SCORE,

    Sham_mean_pct_positive >=
      MIN_SHAM_PCT_POSITIVE,

    Tx_mean_pct_positive >=
      MIN_TX_PCT_POSITIVE,

    pseudobulk_log2FC_Tx_vs_Sham >=
      MIN_LOG2FC_TX_VS_SHAM,

    delta_pct_Tx_minus_Sham >=
      MIN_DELTA_PCT_TX_MINUS_SHAM,

    pseudobulk_pairwise_increase_fraction >=
      MIN_PAIRWISE_INCREASE_FRACTION,

    pct_positive_pairwise_increase_fraction >=
      MIN_PAIRWISE_INCREASE_FRACTION
  ) %>%
  arrange(
    desc(
      Cd163like_TxIncrease_score
    )
  ) %>%
  mutate(
    shortlist_rank =
      row_number()
  )

write.csv(
  shortlist,
  file.path(
    DIR_TABLE,
    "07_STRINGENT_Cd163_UMAPlike_TxIncrease_shortlist_v5.7.1.csv"
  ),
  row.names = FALSE
)

top200 <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_TABLE
  )

write.csv(
  top200,
  file.path(
    DIR_TABLE,
    "08_TOP200_broad_candidates_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 18. Reference-panel positions
# ==============================================================================

reference_positions <- ranking %>%
  filter(
    in_reference_panel
  ) %>%
  arrange(
    final_rank
  )

# Add Cd163 itself as treatment reference.
cd163_treatment <- treatment_metrics %>%
  filter(
    gene ==
      REFERENCE_GENE
  )

write.csv(
  reference_positions,
  file.path(
    DIR_TABLE,
    "09_reference_panel_positions_v5.7.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  cd163_treatment,
  file.path(
    DIR_TABLE,
    "10_Cd163_reference_Sham_to_Tx_after_STRICT_QC_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 19. Ranking plot
# ==============================================================================

top_plot <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_RANK_PLOT
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(
        gene
      )
    )
  )

p_rank <- ggplot(
  top_plot,
  aes(
    x =
      Cd163like_TxIncrease_score,
    y =
      gene
  )
) +
  geom_col() +
  labs(
    title =
      "Broad de novo screen: Cd163-like UMAP pattern + Tx increase",
    x =
      "Integrated score",
    y =
      NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title =
      element_text(
        size = 13,
        face = "bold"
      )
  )

save_pdf(
  p_rank,
  file.path(
    DIR_PLOT,
    "01_TOP40_Cd163_UMAPlike_TxIncrease_ranking_v5.7.1.pdf"
  ),
  9,
  11
)

# ==============================================================================
# 20. Ranking metric heatmap
# ==============================================================================

heat_top <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_RANK_PLOT
  )

heat_gene_order <- heat_top$gene

heat_df <- heat_top %>%
  select(
    gene,
    score_UMAP_spatial,
    score_Tx_pseudobulk_increase,
    score_Tx_positive_fraction_increase,
    score_replicate_consistency,
    Cd163like_TxIncrease_score
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "score"
  )

heat_df$gene <- factor(
  heat_df$gene,
  levels = rev(
    heat_gene_order
  )
)

heat_df$metric <- factor(
  heat_df$metric,
  levels = c(
    "score_UMAP_spatial",
    "score_Tx_pseudobulk_increase",
    "score_Tx_positive_fraction_increase",
    "score_replicate_consistency",
    "Cd163like_TxIncrease_score"
  )
)

p_heat <- ggplot(
  heat_df,
  aes(
    x = metric,
    y = gene,
    fill = score
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      1
    )
  ) +
  labs(
    title =
      "Top 40: spatial similarity and treatment-increase components",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    plot.title =
      element_text(
        size = 13,
        face = "bold"
      )
  )

save_pdf(
  p_heat,
  file.path(
    DIR_PLOT,
    "02_TOP40_metric_heatmap_v5.7.1.pdf"
  ),
  9,
  11
)

# ==============================================================================
# 21. Sample-level heatmaps
# ==============================================================================

top_genes <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_RANK_PLOT
  ) %>%
  pull(
    gene
  )

sample_heat <- sample_long %>%
  filter(
    gene %in%
      top_genes
  ) %>%
  mutate(
    gene =
      factor(
        gene,
        levels = rev(
          top_genes
        )
      ),
    sample =
      factor(
        sample,
        levels =
          EXPECTED_SAMPLES
      )
  )

p_pb <- ggplot(
  sample_heat,
  aes(
    x = sample,
    y = gene,
    fill =
      pseudobulk_log2CPM
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint =
      median(
        sample_heat$pseudobulk_log2CPM,
        na.rm = TRUE
      )
  ) +
  labs(
    title =
      "Top 40: sample-level pseudobulk expression after strict QC",
    x = NULL,
    y = NULL,
    fill = "log2 CPM"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_pb,
  file.path(
    DIR_PLOT,
    "03_TOP40_samplelevel_pseudobulk_heatmap_v5.7.1.pdf"
  ),
  7,
  11
)

p_pct <- ggplot(
  sample_heat,
  aes(
    x = sample,
    y = gene,
    fill =
      pct_positive
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  labs(
    title =
      "Top 40: RNA-positive cell fraction after strict QC",
    x = NULL,
    y = NULL,
    fill = "% positive"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_pct,
  file.path(
    DIR_PLOT,
    "04_TOP40_samplelevel_positive_fraction_heatmap_v5.7.1.pdf"
  ),
  7,
  11
)

# ==============================================================================
# 22. UMAP FeaturePlots of top spatial/treatment candidates
# ==============================================================================

top_feature_genes <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_FEATURE
  ) %>%
  pull(
    gene
  )

# STD reference UMAP panel: directly compare candidate spatial pattern with Cd163.
std_obj <- subset(
  kupffer,
  cells =
    std_cells
)

std_plot_genes <- unique(
  c(
    REFERENCE_GENE,
    top_feature_genes
  )
)

p_std <- FeaturePlot(
  std_obj,
  features =
    std_plot_genes,
  reduction =
    umap_name,
  ncol = 5,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  raster = FALSE,
  pt.size = 0.40
) &
  scale_colour_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) &
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_std,
  file.path(
    DIR_UMAP,
    "01_STD_Cd163_plus_TOP24_UMAP_pattern_v5.7.1.pdf"
  ),
  16,
  14
)

# Sham/Tx split UMAP after strict contamination QC.
clean$condition_v571 <- factor(
  clean$condition_v571,
  levels = c(
    "Sham",
    "Tx"
  )
)

for (
  gene in unique(
    c(
      REFERENCE_GENE,
      top_feature_genes
    )
  )
) {

  p_gene <- FeaturePlot(
    clean,
    features = gene,
    reduction =
      umap_name,
    split.by =
      "condition_v571",
    keep.scale = "all",
    ncol = 2,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    raster = FALSE,
    pt.size = 0.35
  ) &
    scale_colour_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) &
    theme_classic(
      base_size = 9
    )

  save_pdf(
    p_gene,
    file.path(
      DIR_UMAP,
      paste0(
        "Sham_vs_Tx_",
        gene,
        "_after_STRICT_QC_v5.7.1.pdf"
      )
    ),
    9,
    4.5
  )
}

# ==============================================================================
# 23. Direct spatial-metric plot
# ==============================================================================

spatial_plot <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n =
      TOP_N_RANK_PLOT
  ) %>%
  select(
    gene,
    UMAP_grid_mean_spearman,
    UMAP_grid_pct_spearman,
    cell_spearman_with_Cd163,
    binary_jaccard_with_Cd163,
    Cd163_UMAP_spatial_score
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "value"
  )

spatial_plot$gene <- factor(
  spatial_plot$gene,
  levels = rev(
    heat_gene_order
  )
)

p_spatial <- ggplot(
  spatial_plot,
  aes(
    x = metric,
    y = gene,
    fill = value
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "Direct Cd163 UMAP-pattern similarity metrics",
    x = NULL,
    y = NULL,
    fill = "Metric"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

save_pdf(
  p_spatial,
  file.path(
    DIR_PLOT,
    "05_TOP40_direct_Cd163_spatial_metrics_v5.7.1.pdf"
  ),
  10,
  11
)

# ==============================================================================
# 24. Metadata and output index
# ==============================================================================

analysis_metadata <- tibble(
  parameter = c(
    "script_version",
    "input_RDS",
    "assay",
    "condition_source",
    "sample_column",
    "lineage_column",
    "lineage_value",
    "UMAP",
    "discovery_condition",
    "GRID_NX",
    "GRID_NY",
    "MIN_CELLS_PER_GRID_BIN",
    "MIN_DETECTION_FRACTION",
    "strict_QC_core_marker_threshold",
    "weight_UMAP_spatial",
    "weight_Tx_pseudobulk_increase",
    "weight_Tx_positive_fraction_increase",
    "weight_replicate_consistency",
    "weight_grid_mean",
    "weight_grid_pct",
    "weight_cell_spearman",
    "weight_binary_jaccard"
  ),
  value = c(
    "v5.7.1",
    MPHI_RDS,
    ASSAY_USE,
    condition_info$source,
    sample_col,
    lineage_col,
    MACROPHAGE_LINEAGE_VALUE,
    umap_name,
    DISCOVERY_CONDITION,
    GRID_NX,
    GRID_NY,
    MIN_CELLS_PER_GRID_BIN,
    MIN_DETECTION_FRACTION,
    PRIMARY_MIN_CORE_MARKERS,
    W_UMAP_SPATIAL,
    W_TX_PSEUDOBULK_INCREASE,
    W_TX_POSITIVE_FRACTION_INCREASE,
    W_REPLICATE_CONSISTENCY,
    W_GRID_MEAN,
    W_GRID_PCT,
    W_CELL_SPEARMAN,
    W_BINARY_JACCARD
  )
)

write.csv(
  analysis_metadata,
  file.path(
    DIR_LOG,
    "analysis_metadata_v5.7.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_LOG,
    "sessionInfo_v5.7.1.txt"
  )
)

index <- tibble(
  item = c(
    "All-gene Cd163 UMAP spatial similarity",
    "Strict QC summary",
    "Strict QC sample counts",
    "All-gene sample metrics",
    "All-gene Sham-to-Tx increase metrics",
    "Final integrated all-gene ranking",
    "Stringent shortlist",
    "Top 200 candidates",
    "Reference-panel positions",
    "Cd163 treatment reference",
    "Top 40 ranking plot",
    "Top 40 metric heatmap",
    "Top 40 pseudobulk heatmap",
    "Top 40 positive-fraction heatmap",
    "STD Cd163 plus top24 UMAP panel",
    "Top 40 direct spatial metrics"
  ),
  path = c(
    file.path(
      DIR_TABLE,
      "01_ALL_GENES_Cd163_UMAP_spatial_similarity_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "02_STRICT_QC_summary_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "03_STRICT_QC_sample_cell_counts_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "04_ALL_GENES_samplelevel_metrics_after_STRICT_QC_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "05_ALL_GENES_Sham_to_Tx_increase_after_STRICT_QC_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "06_FINAL_ALL_GENES_Cd163_UMAPlike_TxIncrease_ranking_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "07_STRINGENT_Cd163_UMAPlike_TxIncrease_shortlist_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "08_TOP200_broad_candidates_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "09_reference_panel_positions_v5.7.1.csv"
    ),
    file.path(
      DIR_TABLE,
      "10_Cd163_reference_Sham_to_Tx_after_STRICT_QC_v5.7.1.csv"
    ),
    file.path(
      DIR_PLOT,
      "01_TOP40_Cd163_UMAPlike_TxIncrease_ranking_v5.7.1.pdf"
    ),
    file.path(
      DIR_PLOT,
      "02_TOP40_metric_heatmap_v5.7.1.pdf"
    ),
    file.path(
      DIR_PLOT,
      "03_TOP40_samplelevel_pseudobulk_heatmap_v5.7.1.pdf"
    ),
    file.path(
      DIR_PLOT,
      "04_TOP40_samplelevel_positive_fraction_heatmap_v5.7.1.pdf"
    ),
    file.path(
      DIR_UMAP,
      "01_STD_Cd163_plus_TOP24_UMAP_pattern_v5.7.1.pdf"
    ),
    file.path(
      DIR_PLOT,
      "05_TOP40_direct_Cd163_spatial_metrics_v5.7.1.pdf"
    )
  )
)

write.csv(
  index,
  file.path(
    OUTPUT_DIR,
    "Cd163_UMAPlike_TxIncrease_INDEX_v5.7.1.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 25. Final terminal report
# ==============================================================================

msg(
  "DONE."
)

msg(
  "Output: ",
  OUTPUT_DIR
)

cat(
  "\n============================================================\n"
)

cat(
  "Cd163 UMAP-like + Tx-increase broad de novo screen v5.7.1\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Cd163 treatment reference after strict QC:\n"
)

print(
  cd163_treatment
)

cat(
  "\nTOP 40 BROAD CANDIDATES:\n"
)

print(
  ranking %>%
    filter(
      !technical_gene
    ) %>%
    select(
      final_rank,
      gene,
      Cd163like_TxIncrease_score,
      spatial_rank,
      Cd163_UMAP_spatial_score,
      UMAP_grid_mean_spearman,
      UMAP_grid_pct_spearman,
      cell_spearman_with_Cd163,
      binary_jaccard_with_Cd163,
      pseudobulk_log2FC_Tx_vs_Sham,
      Sham_mean_pct_positive,
      Tx_mean_pct_positive,
      delta_pct_Tx_minus_Sham,
      pseudobulk_pairwise_increase_fraction,
      pct_positive_pairwise_increase_fraction,
      in_reference_panel
    ) %>%
    slice_head(
      n = 40
    )
)

cat(
  "\nSTRINGENT SHORTLIST: ",
  nrow(
    shortlist
  ),
  " genes\n",
  sep = ""
)

if (
  nrow(
    shortlist
  ) > 0
) {

  print(
    shortlist %>%
      select(
        shortlist_rank,
        gene,
        Cd163like_TxIncrease_score,
        spatial_rank,
        Cd163_UMAP_spatial_score,
        pseudobulk_log2FC_Tx_vs_Sham,
        Sham_mean_pct_positive,
        Tx_mean_pct_positive,
        delta_pct_Tx_minus_Sham,
        pseudobulk_pairwise_increase_fraction,
        pct_positive_pairwise_increase_fraction
      ) %>%
      slice_head(
        n = 40
      )
  )
}

cat(
  "\n============================================================\n"
)
