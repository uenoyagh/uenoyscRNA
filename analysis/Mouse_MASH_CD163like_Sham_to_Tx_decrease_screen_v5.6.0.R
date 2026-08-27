#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# Cd163-like Sham -> Tx decrease screen
#
# Version: v5.6.0
#
# PURPOSE
#   Prioritize Cd163-like macrophage markers that DECREASE from Sham to Tx.
#
# INPUTS
#   1) Frozen Clean-B macrophage RDS (FINAL v4.14.5 / Res2.0)
#   2) v5.5.0 all-gene Cd163 de novo ranking
#
# DESIGN
#   - Restrict to Kupffer_Macrophage cells.
#   - Compare individual biological samples:
#       Sham1, Sham20 vs Tx17, Tx5
#   - Use sample-level pseudobulk RNA counts (CPM/logCPM).
#   - Use sample-level fraction of cells with RNA counts > 0.
#   - Require/score consistency across the 2 Sham and 2 Tx samples.
#   - Combine these Sham->Tx decrease features with the v5.5.0 Cd163 de novo score.
#
# IMPORTANT
#   - No reintegration, reclustering, or re-UMAP.
#   - No cell-level p-values are used for biological inference.
#   - With n=2 Sham and n=2 Tx, statistics are descriptive/exploratory.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(5600)

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

MPHI_RDS <- MPHI_RDS_CANDIDATES[file.exists(MPHI_RDS_CANDIDATES)][1]

DENOVO_V550 <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163_Positive_deNovo_Marker_Discovery_v5.5.0",
  "Tables",
  "05_ALL_GENES_Cd163_deNovo_ranking_v5.5.0.csv"
)

OUTPUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163like_Sham_to_Tx_Decrease_v5.6.0"
)

DIR_TABLE <- file.path(OUTPUT_DIR, "Tables")
DIR_PLOT <- file.path(OUTPUT_DIR, "Plots")
DIR_UMAP <- file.path(OUTPUT_DIR, "UMAP")
DIR_LOG <- file.path(OUTPUT_DIR, "Logs")

for (d in c(
  OUTPUT_DIR,
  DIR_TABLE,
  DIR_PLOT,
  DIR_UMAP,
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

TARGET_GROUPS <- c(
  "Sham",
  "Tx"
)

EXPECTED_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

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

MACROPHAGE_LINEAGE_VALUE <- "Kupffer_Macrophage"

# Ranking weights.
W_CD163_LIKE <- 0.40
W_PSEUDOBULK_DECREASE <- 0.25
W_POSITIVE_FRACTION_DECREASE <- 0.20
W_DIRECTION_CONSISTENCY <- 0.15

# Stringent shortlist criteria.
MAX_DENOVO_RANK <- 300
MIN_SHAM_PCT_POSITIVE <- 10
MAX_PB_LOG2FC_TX_VS_SHAM <- -0.25
MIN_DELTA_PCT_SHAM_MINUS_TX <- 5
MIN_PAIRWISE_DECREASE_FRACTION <- 0.75

TOP_N_TABLE <- 100
TOP_N_PLOT <- 30
TOP_N_FEATURE <- 12

TECHNICAL_GENE_PATTERN <- "^(mt-|Rpl|Rps|Hba-|Hbb-)"

# Previously notable v5.5 genes, annotation only.
V550_NOTABLE <- c(
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
  "Slc40a1"
)

# Histology-oriented surface/membrane candidates.
# This annotation does NOT affect the data-driven score.
HISTOLOGY_SURFACE_PRIORITY <- c(
  "Timd4",
  "Clec4f",
  "Folr2",
  "Vsig4",
  "Marco",
  "Slc40a1",
  "Stab2",
  "Cdh5"
)

# ==============================================================================
# 3. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
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
    candidates %in% available
  ]

  if (!length(hit)) {
    return(NA_character_)
  }

  hit[[1]]
}

canonical_condition <- function(x) {

  x <- as.character(x)

  out <- case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
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

  for (src in CONDITION_COLUMN_CANDIDATES) {

    if (!src %in% meta_cols) {
      next
    }

    cond <- canonical_condition(
      obj@meta.data[[src]]
    )

    if (!all(is.na(cond))) {
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

safe_rescale01 <- function(x) {

  out <- rep(
    NA_real_,
    length(x)
  )

  ok <- is.finite(x)

  if (!any(ok)) {
    return(out)
  }

  r <- range(
    x[ok],
    na.rm = TRUE
  )

  if (
    !all(is.finite(r)) ||
    diff(r) == 0
  ) {
    out[ok] <- 0.5
    return(out)
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

pairwise_decrease_fraction <- function(
  sham_values,
  tx_values
) {

  sham_values <- sham_values[
    is.finite(sham_values)
  ]

  tx_values <- tx_values[
    is.finite(tx_values)
  ]

  if (
    !length(sham_values) ||
    !length(tx_values)
  ) {
    return(NA_real_)
  }

  comparisons <- outer(
    tx_values,
    sham_values,
    FUN = "<"
  )

  mean(
    comparisons
  )
}

# ==============================================================================
# 4. Validate inputs
# ==============================================================================

if (
  length(MPHI_RDS) == 0L ||
  is.na(MPHI_RDS) ||
  !file.exists(MPHI_RDS)
) {
  stop(
    "Clean-B macrophage RDS not found."
  )
}

if (
  !file.exists(DENOVO_V550)
) {
  stop(
    paste0(
      "v5.5.0 de novo ranking not found:\n",
      DENOVO_V550
    )
  )
}

# ==============================================================================
# 5. Load data
# ==============================================================================

msg(
  "Loading Clean-B macrophage RDS..."
)

mphi <- readRDS(
  MPHI_RDS
)

msg(
  "Loading v5.5.0 de novo ranking..."
)

denovo <- read.csv(
  DENOVO_V550,
  check.names = FALSE
)

if (
  !ASSAY_USE %in%
    Assays(mphi)
) {
  stop(
    "RNA assay missing from macrophage object."
  )
}

DefaultAssay(
  mphi
) <- ASSAY_USE

if (
  !all(
    c(
      "gene",
      "de_novo_rank",
      "de_novo_score"
    ) %in%
      colnames(denovo)
  )
) {
  stop(
    "Required columns missing from v5.5.0 de novo ranking."
  )
}

# ==============================================================================
# 6. Resolve metadata
# ==============================================================================

condition_info <- resolve_condition(
  mphi
)

mphi$condition_v560 <- condition_info$condition

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
  Reductions(mphi),
  UMAP_CANDIDATES
)

if (
  is.na(lineage_col)
) {
  stop(
    "Could not resolve lineage column."
  )
}

if (
  is.na(sample_col)
) {
  stop(
    "Could not resolve biological sample column."
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
  ifelse(
    is.na(umap_name),
    "<not found>",
    umap_name
  )
)

# ==============================================================================
# 7. Select Sham/Tx Kupffer macrophages
# ==============================================================================

lineage <- as.character(
  mphi@meta.data[[lineage_col]]
)

sample_id <- as.character(
  mphi@meta.data[[sample_col]]
)

condition <- as.character(
  mphi$condition_v560
)

keep <- lineage ==
  MACROPHAGE_LINEAGE_VALUE &
  condition %in%
  TARGET_GROUPS

cells_use <- colnames(mphi)[
  keep &
    !is.na(keep)
]

if (
  length(cells_use) < 100
) {
  stop(
    paste0(
      "Too few Sham/Tx Kupffer macrophages: ",
      length(cells_use)
    )
  )
}

obj <- subset(
  mphi,
  cells = cells_use
)

DefaultAssay(
  obj
) <- ASSAY_USE

obj$sample_v560 <- as.character(
  obj@meta.data[[sample_col]]
)

obj$condition_v560 <- as.character(
  obj$condition_v560
)

# Keep only Sham/Tx samples.
obj <- subset(
  obj,
  subset =
    condition_v560 %in%
    TARGET_GROUPS
)

samples_present <- sort(
  unique(
    obj$sample_v560
  )
)

msg(
  "Samples present: ",
  paste(
    samples_present,
    collapse = ", "
  )
)

missing_expected <- setdiff(
  EXPECTED_SAMPLES,
  samples_present
)

if (
  length(missing_expected)
) {
  warning(
    paste0(
      "Expected samples missing: ",
      paste(
        missing_expected,
        collapse = ", "
      )
    )
  )
}

sample_counts <- tibble(
  sample =
    obj$sample_v560,
  condition =
    obj$condition_v560
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
    "00_sample_cell_counts_v5.6.0.csv"
  ),
  row.names = FALSE
)

print(
  sample_counts
)

# ==============================================================================
# 8. RNA matrices
# ==============================================================================

counts_mat <- GetAssayData(
  obj,
  assay = ASSAY_USE,
  layer = "counts"
)

data_mat <- GetAssayData(
  obj,
  assay = ASSAY_USE,
  layer = "data"
)

genes_obj <- rownames(
  counts_mat
)

genes_denovo <- intersect(
  denovo$gene,
  genes_obj
)

genes_eval <- unique(
  c(
    REFERENCE_GENE,
    genes_denovo
  )
)

if (
  !REFERENCE_GENE %in%
    genes_obj
) {
  stop(
    "Cd163 is absent from RNA assay."
  )
}

# ==============================================================================
# 9. Sample-level pseudobulk counts
# ==============================================================================

sample_vector <- obj$sample_v560

sample_levels <- sort(
  unique(
    sample_vector
  )
)

pb_counts <- matrix(
  0,
  nrow = length(
    genes_eval
  ),
  ncol = length(
    sample_levels
  ),
  dimnames = list(
    genes_eval,
    sample_levels
  )
)

sample_ncells <- setNames(
  integer(
    length(
      sample_levels
    )
  ),
  sample_levels
)

for (
  s in sample_levels
) {

  idx <- which(
    sample_vector == s
  )

  sample_ncells[[s]] <- length(
    idx
  )

  pb_counts[
    ,
    s
  ] <- as.numeric(
    Matrix::rowSums(
      counts_mat[
        genes_eval,
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
) * 1e6

pb_logcpm <- log2(
  pb_cpm + 1
)

# ==============================================================================
# 10. Sample-level positive fractions and mean log-normalized expression
# ==============================================================================

sample_pct <- matrix(
  NA_real_,
  nrow = length(
    genes_eval
  ),
  ncol = length(
    sample_levels
  ),
  dimnames = list(
    genes_eval,
    sample_levels
  )
)

sample_mean_data <- matrix(
  NA_real_,
  nrow = length(
    genes_eval
  ),
  ncol = length(
    sample_levels
  ),
  dimnames = list(
    genes_eval,
    sample_levels
  )
)

for (
  s in sample_levels
) {

  idx <- which(
    sample_vector == s
  )

  sample_pct[
    ,
    s
  ] <- 100 *
    as.numeric(
      Matrix::rowMeans(
        counts_mat[
          genes_eval,
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
      data_mat[
        genes_eval,
        idx,
        drop = FALSE
      ]
    )
  )
}

# ==============================================================================
# 11. Long-form per-sample table
# ==============================================================================

sample_condition_map <- tibble(
  sample =
    obj$sample_v560,
  condition =
    obj$condition_v560
) %>%
  distinct() %>%
  arrange(
    condition,
    sample
  )

sample_metric_long <- bind_rows(
  lapply(
    sample_levels,
    function(s) {

      tibble(
        gene = genes_eval,
        sample = s,
        pseudobulk_count =
          pb_counts[
            genes_eval,
            s
          ],
        pseudobulk_CPM =
          pb_cpm[
            genes_eval,
            s
          ],
        pseudobulk_log2CPM =
          pb_logcpm[
            genes_eval,
            s
          ],
        pct_positive =
          sample_pct[
            genes_eval,
            s
          ],
        mean_log_normalized_expression =
          sample_mean_data[
            genes_eval,
            s
          ],
        n_cells =
          sample_ncells[[s]]
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
    pseudobulk_count,
    pseudobulk_CPM,
    pseudobulk_log2CPM,
    pct_positive,
    mean_log_normalized_expression
  )

write.csv(
  sample_metric_long,
  file.path(
    DIR_TABLE,
    "01_gene_by_sample_metrics_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 12. Sham vs Tx summary
# ==============================================================================

summary_list <- vector(
  "list",
  length(
    genes_eval
  )
)

names(
  summary_list
) <- genes_eval

for (
  i in seq_along(
    genes_eval
  )
) {

  gene <- genes_eval[[i]]

  df <- sample_metric_long %>%
    filter(
      gene == !!gene
    )

  sham <- df %>%
    filter(
      condition == "Sham"
    )

  tx <- df %>%
    filter(
      condition == "Tx"
    )

  mean_cpm_sham <- mean(
    sham$pseudobulk_CPM,
    na.rm = TRUE
  )

  mean_cpm_tx <- mean(
    tx$pseudobulk_CPM,
    na.rm = TRUE
  )

  pb_log2fc <- log2(
    (
      mean_cpm_tx +
        0.5
    ) /
      (
        mean_cpm_sham +
          0.5
      )
  )

  mean_logcpm_sham <- mean(
    sham$pseudobulk_log2CPM,
    na.rm = TRUE
  )

  mean_logcpm_tx <- mean(
    tx$pseudobulk_log2CPM,
    na.rm = TRUE
  )

  mean_pct_sham <- mean(
    sham$pct_positive,
    na.rm = TRUE
  )

  mean_pct_tx <- mean(
    tx$pct_positive,
    na.rm = TRUE
  )

  delta_pct_sham_minus_tx <-
    mean_pct_sham -
    mean_pct_tx

  pb_pairwise <- pairwise_decrease_fraction(
    sham$pseudobulk_log2CPM,
    tx$pseudobulk_log2CPM
  )

  pct_pairwise <- pairwise_decrease_fraction(
    sham$pct_positive,
    tx$pct_positive
  )

  pb_strict <- if (
    nrow(sham) > 0 &&
    nrow(tx) > 0
  ) {
    max(
      tx$pseudobulk_log2CPM,
      na.rm = TRUE
    ) <
      min(
        sham$pseudobulk_log2CPM,
        na.rm = TRUE
      )
  } else {
    FALSE
  }

  pct_strict <- if (
    nrow(sham) > 0 &&
    nrow(tx) > 0
  ) {
    max(
      tx$pct_positive,
      na.rm = TRUE
    ) <
      min(
        sham$pct_positive,
        na.rm = TRUE
      )
  } else {
    FALSE
  }

  summary_list[[i]] <- tibble(
    gene = gene,

    Sham_mean_CPM =
      mean_cpm_sham,

    Tx_mean_CPM =
      mean_cpm_tx,

    pseudobulk_log2FC_Tx_vs_Sham =
      pb_log2fc,

    Sham_mean_log2CPM =
      mean_logcpm_sham,

    Tx_mean_log2CPM =
      mean_logcpm_tx,

    Sham_mean_pct_positive =
      mean_pct_sham,

    Tx_mean_pct_positive =
      mean_pct_tx,

    delta_pct_Sham_minus_Tx =
      delta_pct_sham_minus_tx,

    pseudobulk_pairwise_decrease_fraction =
      pb_pairwise,

    pct_positive_pairwise_decrease_fraction =
      pct_pairwise,

    pseudobulk_all_Tx_below_all_Sham =
      pb_strict,

    pct_positive_all_Tx_below_all_Sham =
      pct_strict,

    both_metrics_strictly_decrease =
      pb_strict &&
      pct_strict
  )
}

sham_tx_summary <- bind_rows(
  summary_list
)

# Add each sample as explicit columns.
pb_wide <- sample_metric_long %>%
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

pct_wide <- sample_metric_long %>%
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

sham_tx_summary <- sham_tx_summary %>%
  left_join(
    pb_wide,
    by = "gene"
  ) %>%
  left_join(
    pct_wide,
    by = "gene"
  )

write.csv(
  sham_tx_summary,
  file.path(
    DIR_TABLE,
    "02_ALL_GENES_Sham_vs_Tx_samplelevel_summary_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 13. Integrate with v5.5.0 Cd163-like score
# ==============================================================================

ranking <- denovo %>%
  inner_join(
    sham_tx_summary,
    by = "gene"
  ) %>%
  mutate(

    technical_gene =
      grepl(
        TECHNICAL_GENE_PATTERN,
        gene
      ),

    notable_v550 =
      gene %in%
      V550_NOTABLE,

    histology_surface_priority =
      gene %in%
      HISTOLOGY_SURFACE_PRIORITY,

    score_Cd163_like =
      safe_rescale01(
        de_novo_score
      ),

    score_pseudobulk_decrease =
      safe_rescale01(
        pmax(
          -pseudobulk_log2FC_Tx_vs_Sham,
          0
        )
      ),

    score_positive_fraction_decrease =
      safe_rescale01(
        pmax(
          delta_pct_Sham_minus_Tx,
          0
        )
      ),

    direction_consistency =
      rowMeans(
        cbind(
          pseudobulk_pairwise_decrease_fraction,
          pct_positive_pairwise_decrease_fraction
        ),
        na.rm = TRUE
      ),

    score_direction_consistency =
      direction_consistency,

    Sham_to_Tx_decrease_score =
      W_CD163_LIKE *
        score_Cd163_like +
      W_PSEUDOBULK_DECREASE *
        score_pseudobulk_decrease +
      W_POSITIVE_FRACTION_DECREASE *
        score_positive_fraction_decrease +
      W_DIRECTION_CONSISTENCY *
        score_direction_consistency
  ) %>%
  arrange(
    desc(
      Sham_to_Tx_decrease_score
    )
  ) %>%
  mutate(
    decrease_rank =
      row_number()
  )

write.csv(
  ranking,
  file.path(
    DIR_TABLE,
    "03_FINAL_Cd163like_Sham_to_Tx_decrease_ranking_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 14. Stringent shortlist
# ==============================================================================

shortlist <- ranking %>%
  filter(
    !technical_gene,
    de_novo_rank <=
      MAX_DENOVO_RANK,
    Sham_mean_pct_positive >=
      MIN_SHAM_PCT_POSITIVE,
    pseudobulk_log2FC_Tx_vs_Sham <=
      MAX_PB_LOG2FC_TX_VS_SHAM,
    delta_pct_Sham_minus_Tx >=
      MIN_DELTA_PCT_SHAM_MINUS_TX,
    pseudobulk_pairwise_decrease_fraction >=
      MIN_PAIRWISE_DECREASE_FRACTION,
    pct_positive_pairwise_decrease_fraction >=
      MIN_PAIRWISE_DECREASE_FRACTION
  ) %>%
  arrange(
    desc(
      Sham_to_Tx_decrease_score
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
    "04_STRINGENT_Cd163like_Sham_to_Tx_decrease_shortlist_v5.6.0.csv"
  ),
  row.names = FALSE
)

top100 <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_TABLE
  )

write.csv(
  top100,
  file.path(
    DIR_TABLE,
    "05_TOP100_Cd163like_Sham_to_Tx_decrease_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 15. Reference Cd163 sample-level result
# ==============================================================================

cd163_reference <- sham_tx_summary %>%
  filter(
    gene == REFERENCE_GENE
  )

write.csv(
  cd163_reference,
  file.path(
    DIR_TABLE,
    "06_Cd163_reference_Sham_vs_Tx_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 16. Previously notable v5.5 genes
# ==============================================================================

notable_result <- ranking %>%
  filter(
    notable_v550
  ) %>%
  arrange(
    decrease_rank
  )

write.csv(
  notable_result,
  file.path(
    DIR_TABLE,
    "07_v5.5_notable_genes_Sham_to_Tx_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 17. Histology-oriented surface candidates
# ==============================================================================

surface_result <- ranking %>%
  filter(
    histology_surface_priority
  ) %>%
  arrange(
    decrease_rank
  )

write.csv(
  surface_result,
  file.path(
    DIR_TABLE,
    "08_surface_marker_candidates_Sham_to_Tx_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 18. Top ranking plot
# ==============================================================================

top_plot <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_PLOT
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(gene)
    )
  )

p_rank <- ggplot(
  top_plot,
  aes(
    x = Sham_to_Tx_decrease_score,
    y = gene
  )
) +
  geom_col() +
  labs(
    title =
      "Cd163-like macrophage markers decreasing from Sham to Tx",
    x =
      "Integrated Sham-to-Tx decrease score",
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_rank,
  file.path(
    DIR_PLOT,
    "01_TOP30_Cd163like_Sham_to_Tx_decrease_ranking_v5.6.0.pdf"
  ),
  9,
  9
)

# ==============================================================================
# 19. Metric heatmap
# ==============================================================================

heat_df <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_PLOT
  ) %>%
  select(
    gene,
    score_Cd163_like,
    score_pseudobulk_decrease,
    score_positive_fraction_decrease,
    score_direction_consistency,
    Sham_to_Tx_decrease_score
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "score"
  )

heat_gene_order <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_PLOT
  ) %>%
  pull(
    gene
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
    "score_Cd163_like",
    "score_pseudobulk_decrease",
    "score_positive_fraction_decrease",
    "score_direction_consistency",
    "Sham_to_Tx_decrease_score"
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
      "Top Cd163-like Sham-to-Tx decrease candidates",
    x = NULL,
    y = NULL,
    fill = "Score"
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
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_heat,
  file.path(
    DIR_PLOT,
    "02_TOP30_decrease_metric_heatmap_v5.6.0.pdf"
  ),
  9,
  9
)

# ==============================================================================
# 20. Per-sample pseudobulk heatmap
# ==============================================================================

top_genes <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_PLOT
  ) %>%
  pull(
    gene
  )

pb_heat <- sample_metric_long %>%
  filter(
    gene %in%
      top_genes
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(
        top_genes
      )
    ),
    sample = factor(
      sample,
      levels = c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
      )
    )
  )

p_pb <- ggplot(
  pb_heat,
  aes(
    x = sample,
    y = gene,
    fill = pseudobulk_log2CPM
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = median(
      pb_heat$pseudobulk_log2CPM,
      na.rm = TRUE
    )
  ) +
  labs(
    title =
      "Sample-level pseudobulk expression",
    x = NULL,
    y = NULL,
    fill = "log2 CPM"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_pb,
  file.path(
    DIR_PLOT,
    "03_TOP30_samplelevel_pseudobulk_heatmap_v5.6.0.pdf"
  ),
  7,
  9
)

# ==============================================================================
# 21. Per-sample positive fraction heatmap
# ==============================================================================

p_pct <- ggplot(
  pb_heat,
  aes(
    x = sample,
    y = gene,
    fill = pct_positive
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
      "Sample-level RNA-positive cell fraction",
    x = NULL,
    y = NULL,
    fill = "% positive"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_pct,
  file.path(
    DIR_PLOT,
    "04_TOP30_samplelevel_positive_fraction_heatmap_v5.6.0.pdf"
  ),
  7,
  9
)

# ==============================================================================
# 22. Top candidates: individual-sample violins
# ==============================================================================

top_feature_genes <- ranking %>%
  filter(
    !technical_gene
  ) %>%
  slice_head(
    n = TOP_N_FEATURE
  ) %>%
  pull(
    gene
  )

plot_genes <- unique(
  c(
    REFERENCE_GENE,
    top_feature_genes
  )
)

obj$sample_v560 <- factor(
  obj$sample_v560,
  levels = c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )
)

p_vln <- VlnPlot(
  obj,
  features = plot_genes,
  group.by = "sample_v560",
  assay = ASSAY_USE,
  layer = "data",
  pt.size = 0,
  ncol = 4
) &
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_vln,
  file.path(
    DIR_PLOT,
    "05_TOP12_plus_Cd163_violin_by_sample_v5.6.0.pdf"
  ),
  15,
  11
)

# ==============================================================================
# 23. Sham/Tx FeaturePlots
# ==============================================================================

if (
  !is.na(umap_name)
) {

  obj$condition_v560 <- factor(
    obj$condition_v560,
    levels = c(
      "Sham",
      "Tx"
    )
  )

  for (
    gene in plot_genes
  ) {

    p <- FeaturePlot(
      obj,
      features = gene,
      reduction = umap_name,
      split.by = "condition_v560",
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
      p,
      file.path(
        DIR_UMAP,
        paste0(
          "Sham_vs_Tx_",
          gene,
          "_FeaturePlot_v5.6.0.pdf"
        )
      ),
      9,
      4.5
    )
  }
}

# ==============================================================================
# 24. Explicit v5.5 notable-gene heatmap
# ==============================================================================

notable_heat <- sample_metric_long %>%
  filter(
    gene %in%
      V550_NOTABLE
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(
        V550_NOTABLE[
          V550_NOTABLE %in%
            gene
        ]
      )
    ),
    sample = factor(
      sample,
      levels = c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
      )
    )
  )

if (
  nrow(notable_heat) > 0
) {

  p_notable <- ggplot(
    notable_heat,
    aes(
      x = sample,
      y = gene,
      fill = pseudobulk_log2CPM
    )
  ) +
    geom_tile(
      linewidth = 0.2
    ) +
    scale_fill_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = median(
        notable_heat$pseudobulk_log2CPM,
        na.rm = TRUE
      )
    ) +
    labs(
      title =
        "v5.5 Cd163-like markers: Sham vs Tx",
      x = NULL,
      y = NULL,
      fill = "log2 CPM"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold"
      )
    )

  save_pdf(
    p_notable,
    file.path(
      DIR_PLOT,
      "06_v5.5_notable_genes_samplelevel_heatmap_v5.6.0.pdf"
    ),
    7,
    7
  )
}

# ==============================================================================
# 25. Analysis metadata
# ==============================================================================

analysis_metadata <- tibble(
  parameter = c(
    "script_version",
    "input_RDS",
    "input_deNovo_ranking",
    "assay",
    "condition_source",
    "sample_column",
    "lineage_column",
    "lineage_value",
    "target_groups",
    "expected_samples",
    "n_cells",
    "weight_Cd163_like",
    "weight_pseudobulk_decrease",
    "weight_positive_fraction_decrease",
    "weight_direction_consistency",
    "max_deNovo_rank_shortlist",
    "min_Sham_pct_positive",
    "max_pb_log2FC_Tx_vs_Sham",
    "min_delta_pct_Sham_minus_Tx",
    "min_pairwise_decrease_fraction"
  ),
  value = c(
    "v5.6.0",
    MPHI_RDS,
    DENOVO_V550,
    ASSAY_USE,
    condition_info$source,
    sample_col,
    lineage_col,
    MACROPHAGE_LINEAGE_VALUE,
    paste(
      TARGET_GROUPS,
      collapse = ";"
    ),
    paste(
      EXPECTED_SAMPLES,
      collapse = ";"
    ),
    ncol(obj),
    W_CD163_LIKE,
    W_PSEUDOBULK_DECREASE,
    W_POSITIVE_FRACTION_DECREASE,
    W_DIRECTION_CONSISTENCY,
    MAX_DENOVO_RANK,
    MIN_SHAM_PCT_POSITIVE,
    MAX_PB_LOG2FC_TX_VS_SHAM,
    MIN_DELTA_PCT_SHAM_MINUS_TX,
    MIN_PAIRWISE_DECREASE_FRACTION
  )
)

write.csv(
  analysis_metadata,
  file.path(
    DIR_LOG,
    "analysis_metadata_v5.6.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_LOG,
    "sessionInfo_v5.6.0.txt"
  )
)

# ==============================================================================
# 26. Output index
# ==============================================================================

index <- tibble(
  item = c(
    "Sample cell counts",
    "Gene by sample metrics",
    "All-gene Sham vs Tx summary",
    "Final integrated ranking",
    "Stringent shortlist",
    "Top 100",
    "Cd163 reference",
    "v5.5 notable genes",
    "Surface-marker candidates",
    "Top 30 ranking plot",
    "Top 30 metric heatmap",
    "Top 30 pseudobulk heatmap",
    "Top 30 positive-fraction heatmap",
    "Top 12 plus Cd163 violins",
    "v5.5 notable sample heatmap"
  ),
  path = c(
    file.path(
      DIR_TABLE,
      "00_sample_cell_counts_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "01_gene_by_sample_metrics_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "02_ALL_GENES_Sham_vs_Tx_samplelevel_summary_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "03_FINAL_Cd163like_Sham_to_Tx_decrease_ranking_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "04_STRINGENT_Cd163like_Sham_to_Tx_decrease_shortlist_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "05_TOP100_Cd163like_Sham_to_Tx_decrease_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "06_Cd163_reference_Sham_vs_Tx_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "07_v5.5_notable_genes_Sham_to_Tx_v5.6.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "08_surface_marker_candidates_Sham_to_Tx_v5.6.0.csv"
    ),
    file.path(
      DIR_PLOT,
      "01_TOP30_Cd163like_Sham_to_Tx_decrease_ranking_v5.6.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "02_TOP30_decrease_metric_heatmap_v5.6.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "03_TOP30_samplelevel_pseudobulk_heatmap_v5.6.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "04_TOP30_samplelevel_positive_fraction_heatmap_v5.6.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "05_TOP12_plus_Cd163_violin_by_sample_v5.6.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "06_v5.5_notable_genes_samplelevel_heatmap_v5.6.0.pdf"
    )
  )
)

write.csv(
  index,
  file.path(
    OUTPUT_DIR,
    "Cd163like_Sham_to_Tx_INDEX_v5.6.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 27. Final report
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
  "Cd163-like Sham -> Tx decrease screen v5.6.0\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Cd163 reference:\n"
)

print(
  cd163_reference
)

cat(
  "\nTOP 30 CANDIDATES:\n"
)

print(
  ranking %>%
    filter(
      !technical_gene
    ) %>%
    select(
      decrease_rank,
      gene,
      Sham_to_Tx_decrease_score,
      de_novo_rank,
      de_novo_score,
      pseudobulk_log2FC_Tx_vs_Sham,
      Sham_mean_pct_positive,
      Tx_mean_pct_positive,
      delta_pct_Sham_minus_Tx,
      pseudobulk_pairwise_decrease_fraction,
      pct_positive_pairwise_decrease_fraction,
      both_metrics_strictly_decrease,
      histology_surface_priority
    ) %>%
    slice_head(
      n = 30
    )
)

cat(
  "\nSTRINGENT SHORTLIST: ",
  nrow(shortlist),
  " genes\n",
  sep = ""
)

if (
  nrow(shortlist) > 0
) {

  print(
    shortlist %>%
      select(
        shortlist_rank,
        gene,
        Sham_to_Tx_decrease_score,
        de_novo_rank,
        pseudobulk_log2FC_Tx_vs_Sham,
        Sham_mean_pct_positive,
        Tx_mean_pct_positive,
        delta_pct_Sham_minus_Tx,
        pseudobulk_pairwise_decrease_fraction,
        pct_positive_pairwise_decrease_fraction,
        histology_surface_priority
      ) %>%
      slice_head(
        n = 30
      )
  )
}

cat(
  "\n============================================================\n"
)
