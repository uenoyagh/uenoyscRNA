#!/usr/bin/env Rscript

# ============================================================
# Mouse MASH Mphi-HSC Human-axis validation v6.4.1
#
# PURPOSE
#   Validate in Mouse Sham vs Tx the three Human MASH
#   macrophage-to-HSC candidate axes:
#
#     1) Spp1   -- primary mechanistic axis
#     2) Cdcp1  -- strongest novel TREM2-associated candidate
#     3) Tgfb1  -- robust TREM2-state-associated auxiliary axis
#
# PRIMARY QUESTIONS
#   A. Does Tx reduce ligand burden in Lipid-associated/TREM2-Mphi?
#   B. Does Tx reduce HSC activation/fibrogenic output?
#   C. Across the four biological samples, do sender burden and
#      HSC response move in the expected same direction?
#
# IMPORTANT STATISTICAL NOTE
#   Sham n=2 and Tx n=2.
#   Therefore:
#     - no inferential group P-values are used as primary evidence
#     - sample is the biological unit
#     - group effects are descriptive
#     - four-sample Spearman rho values are descriptive only
#
# HSC FRAMEWORK
#   Identity:
#     qHSC / aHSC
#
#   Functional response:
#     IEG / initiatory pseudobulk score
#     ECM / fibrogenic pseudobulk score
#
#   Existing 3-state annotation is also retained:
#     qHSC
#     ECM-activated HSC
#     Contractile HSC
#
# No reclustering.
# No source RDS overwrite.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
})

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG <- file.path(OUT, "figures")
DIR_OBJ <- file.path(OUT, "objects")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")

for (d in c(
  DIR_TABLE,
  DIR_FIG,
  DIR_OBJ,
  DIR_LOG,
  DIR_PROV
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LOG_FILE <- file.path(
  DIR_LOG,
  "Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1.log"
)

zz <- file(
  LOG_FILE,
  open = "wt"
)

sink(
  zz,
  type = "output",
  split = TRUE
)

sink(
  zz,
  type = "message"
)

cat("============================================================\n")
cat("Mouse MASH Mphi-HSC Human-axis validation v6.4.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(MPHI_RDS)
)

stopifnot(
  file.exists(HSC_RDS)
)

# ------------------------------------------------------------
# 1. Constants
# ------------------------------------------------------------

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

sample_to_group <- c(
  "Sham1" = "Sham",
  "Sham20" = "Sham",
  "Tx17" = "Tx",
  "Tx5" = "Tx"
)

MPHI_STATE_COL_PRIMARY <-
  "macrophage_class_Res2_FINAL_v4145_char"

MPHI_STATE_COL_FALLBACK <-
  "macrophage_class_Res2_FINAL_v4145"

MPHI_SAMPLE_CANDIDATES <- c(
  "sample_4group",
  "sample_for_annotation",
  "sample_FIXED2",
  "sample"
)

PRIMARY_MPHI_STATES <- c(
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Inflammatory-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Repair/Resolution-Mphi"
)

TREM2_STATE <-
  "Lipid-associated/TREM2-Mphi"

FOCUS_LIGANDS <- c(
  "Spp1",
  "Cdcp1",
  "Tgfb1"
)

HSC_SAMPLE_COL <-
  "sample_hsc_v610"

HSC_STATE2_COL <-
  "HSC_state2_v612"

HSC_STATE3_COL <-
  "HSC_state3_v612"

HSC_QC_COL <-
  "qc_primary_include"

IEG_GENES <- c(
  "Thbs1",
  "Timp1",
  "Junb",
  "Fos",
  "Fosb"
)

ECM_GENES <- c(
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Col5a1",
  "Col5a2",
  "Fn1",
  "Sparc",
  "Thbs2",
  "Loxl2",
  "Timp1",
  "Tgfbi",
  "Cthrc1",
  "Serpine1",
  "Postn"
)

# Human-supported receptor axes to check in Mouse HSC
RECEPTOR_MAP <- list(
  Spp1 = c(
    "Itga9",
    "Itgb1",
    "Itgav",
    "Itgb5",
    "Itga4",
    "Itga5"
  ),
  Cdcp1 = c(
    "Egfr"
  ),
  Tgfb1 = c(
    "Tgfbr1",
    "Tgfbr2",
    "Tgfbr3",
    "Itgb1",
    "Itgav"
  )
)

# ------------------------------------------------------------
# 2. Load objects
# ------------------------------------------------------------

mphi <- readRDS(
  MPHI_RDS
)

hsc <- readRDS(
  HSC_RDS
)

# Primary-Clean integrity check
EXPECTED_HSC_PRIMARY <- 12941L

if (ncol(hsc) != EXPECTED_HSC_PRIMARY) {
  stop(
    "Unexpected Primary-Clean HSC cell count: expected ",
    EXPECTED_HSC_PRIMARY,
    ", observed ",
    ncol(hsc)
  )
}

if (!"qc_primary_include" %in% colnames(hsc[[]])) {
  stop(
    "Primary-Clean HSC RDS is missing qc_primary_include."
  )
}

qc_check <- hsc[[]][["qc_primary_include"]]

qc_check_logical <- if (is.logical(qc_check)) {
  qc_check
} else if (is.numeric(qc_check)) {
  !is.na(qc_check) & qc_check != 0
} else {
  tolower(trimws(as.character(qc_check))) %in%
    c("true", "t", "1", "yes", "y", "include", "included", "primary", "pass")
}

if (!all(qc_check_logical)) {
  stop(
    "Primary-Clean HSC RDS contains qc_primary_include != TRUE."
  )
}

DefaultAssay(mphi) <- "RNA"
DefaultAssay(hsc) <- "RNA"

mphi <- NormalizeData(
  mphi,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 3. Resolve Mphi metadata columns
# ------------------------------------------------------------

mphi_meta_cols <- colnames(
  mphi[[]]
)

mphi_state_col <- if (
  MPHI_STATE_COL_PRIMARY %in%
    mphi_meta_cols
) {
  MPHI_STATE_COL_PRIMARY
} else if (
  MPHI_STATE_COL_FALLBACK %in%
    mphi_meta_cols
) {
  MPHI_STATE_COL_FALLBACK
} else {
  stop(
    "Final Mphi annotation column not found."
  )
}

mphi_sample_col <- MPHI_SAMPLE_CANDIDATES[
  MPHI_SAMPLE_CANDIDATES %in%
    mphi_meta_cols
][1]

if (
  length(mphi_sample_col) == 0 ||
  is.na(mphi_sample_col)
) {
  stop(
    "No supported Mphi sample column found."
  )
}

mmd <- mphi[[]] %>%
  rownames_to_column(
    "cell"
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[mphi_sample_col]]
      ),
    sender_state =
      as.character(
        .data[[mphi_state_col]]
      )
  ) %>%
  filter(
    sample %in%
      SAMPLES
  ) %>%
  mutate(
    group =
      unname(
        sample_to_group[sample]
      )
  )

# ------------------------------------------------------------
# 4. Validate Mphi states/samples
# ------------------------------------------------------------

mphi_sample_inventory <- mmd %>%
  count(
    sample,
    group,
    sender_state,
    name = "cells"
  ) %>%
  arrange(
    sample,
    sender_state
  )

write_tsv(
  mphi_sample_inventory,
  file.path(
    DIR_TABLE,
    "01_Mphi_sample_state_inventory.tsv"
  )
)

missing_states <- setdiff(
  PRIMARY_MPHI_STATES,
  unique(
    mmd$sender_state
  )
)

if (
  length(missing_states) > 0
) {
  stop(
    "Missing expected primary Mphi states: ",
    paste(
      missing_states,
      collapse = ";"
    )
  )
}

if (
  !all(
    SAMPLES %in%
      unique(
        mmd$sample
      )
  )
) {
  stop(
    "Not all Sham/Tx samples are present in Mphi object."
  )
}

# ------------------------------------------------------------
# 5. Restrict to primary five Mphi states
# ------------------------------------------------------------

mmd_primary <- mmd %>%
  filter(
    sender_state %in%
      PRIMARY_MPHI_STATES
  )

primary_cells <- mmd_primary$cell

mphi_primary <- subset(
  mphi,
  cells = primary_cells
)

# ------------------------------------------------------------
# 6. Focus-ligand availability
# ------------------------------------------------------------

ligand_audit <- tibble(
  ligand =
    FOCUS_LIGANDS,
  measurable =
    FOCUS_LIGANDS %in%
    rownames(
      mphi_primary
    )
)

write_tsv(
  ligand_audit,
  file.path(
    DIR_TABLE,
    "02_focus_ligand_availability.tsv"
  )
)

if (
  !all(
    ligand_audit$measurable
  )
) {
  stop(
    "Missing focus ligand(s): ",
    paste(
      ligand_audit$ligand[
        !ligand_audit$measurable
      ],
      collapse = ";"
    )
  )
}

# ------------------------------------------------------------
# 7. Primary-Mphi ligand burden
# ------------------------------------------------------------

mphi_expr <- FetchData(
  mphi_primary,
  vars = FOCUS_LIGANDS,
  layer = "data"
) %>%
  rownames_to_column(
    "cell"
  )

mphi_count <- FetchData(
  mphi_primary,
  vars = FOCUS_LIGANDS,
  layer = "counts"
) %>%
  rownames_to_column(
    "cell"
  )

mmd_primary2 <- mmd_primary %>%
  filter(
    cell %in%
      colnames(
        mphi_primary
      )
  )

total_primary_n <- mmd_primary2 %>%
  count(
    sample,
    group,
    name = "total_primary_mphi_n"
  )

state_n <- mmd_primary2 %>%
  count(
    sample,
    group,
    sender_state,
    name = "state_n"
  )

expr_long <- mphi_expr %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "expr"
  ) %>%
  left_join(
    mmd_primary2,
    by = "cell"
  )

count_long <- mphi_count %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "count"
  ) %>%
  left_join(
    mmd_primary2,
    by = "cell"
  )

expr_sum <- expr_long %>%
  group_by(
    sample,
    group,
    sender_state,
    ligand
  ) %>%
  summarise(
    ligand_mean =
      mean(
        expr,
        na.rm = TRUE
      ),
    ligand_expr_sum =
      sum(
        expr,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

det_sum <- count_long %>%
  group_by(
    sample,
    group,
    sender_state,
    ligand
  ) %>%
  summarise(
    ligand_detection =
      mean(
        count > 0,
        na.rm = TRUE
      ),
    ligand_detected_n =
      sum(
        count > 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

burden <- expr_sum %>%
  full_join(
    det_sum,
    by = c(
      "sample",
      "group",
      "sender_state",
      "ligand"
    )
  ) %>%
  left_join(
    state_n,
    by = c(
      "sample",
      "group",
      "sender_state"
    )
  ) %>%
  left_join(
    total_primary_n,
    by = c(
      "sample",
      "group"
    )
  ) %>%
  group_by(
    sample,
    group,
    ligand
  ) %>%
  mutate(
    total_ligand_expr_sum =
      sum(
        ligand_expr_sum,
        na.rm = TRUE
      ),
    total_ligand_detected_n =
      sum(
        ligand_detected_n,
        na.rm = TRUE
      )
  ) %>%
  ungroup() %>%
  mutate(
    state_fraction =
      state_n /
      total_primary_mphi_n,

    ligand_burden_expr =
      ligand_expr_sum /
      total_primary_mphi_n,

    ligand_burden_detection =
      ligand_detected_n /
      total_primary_mphi_n,

    ligand_expression_share =
      ifelse(
        total_ligand_expr_sum > 0,
        ligand_expr_sum /
          total_ligand_expr_sum,
        0
      ),

    ligand_positive_cell_share =
      ifelse(
        total_ligand_detected_n > 0,
        ligand_detected_n /
          total_ligand_detected_n,
        0
      )
  ) %>%
  arrange(
    ligand,
    sample,
    sender_state
  )

write_tsv(
  burden,
  file.path(
    DIR_TABLE,
    "03_focus_ligand_burden_by_sample_state.tsv"
  )
)

trem2_burden <- burden %>%
  filter(
    sender_state ==
      TREM2_STATE
  ) %>%
  arrange(
    ligand,
    sample
  )

write_tsv(
  trem2_burden,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_TREM2_focus_ligand_burden.tsv"
  )
)

# ------------------------------------------------------------
# 8. Descriptive Sham vs Tx Mphi effects
# ------------------------------------------------------------

safe_log2fc <- function(
  tx,
  sham,
  pseudocount = 1e-6
) {
  log2(
    (
      tx +
        pseudocount
    ) /
      (
        sham +
          pseudocount
      )
  )
}

MOUSE_SENDER_METRICS <- c(
  "state_fraction",
  "ligand_mean",
  "ligand_detection",
  "ligand_burden_expr",
  "ligand_burden_detection",
  "ligand_expression_share"
)

sender_group_summary <- trem2_burden %>%
  select(
    sample,
    group,
    ligand,
    all_of(
      MOUSE_SENDER_METRICS
    )
  ) %>%
  pivot_longer(
    cols =
      all_of(
        MOUSE_SENDER_METRICS
      ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(
    ligand,
    metric,
    group
  ) %>%
  summarise(
    n_samples = n(),
    mean_value =
      mean(
        value,
        na.rm = TRUE
      ),
    sample1 =
      value[1],
    sample2 =
      value[
        min(
          2,
          length(
            value
          )
        )
      ],
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      n_samples,
      mean_value,
      sample1,
      sample2
    )
  ) %>%
  mutate(
    Tx_minus_Sham =
      mean_value_Tx -
      mean_value_Sham,

    log2_Tx_over_Sham =
      safe_log2fc(
        mean_value_Tx,
        mean_value_Sham
      )
  )

write_tsv(
  sender_group_summary,
  file.path(
    DIR_TABLE,
    "05_PRIMARY_TREM2_ligand_Sham_vs_Tx_descriptive.tsv"
  )
)

# ------------------------------------------------------------
# 9. HSC primary-clean subset
# ------------------------------------------------------------

hsc_meta_cols <- colnames(
  hsc[[]]
)

required_hsc_cols <- c(
  HSC_SAMPLE_COL,
  HSC_STATE2_COL,
  HSC_STATE3_COL,
  HSC_QC_COL
)

missing_hsc_cols <- setdiff(
  required_hsc_cols,
  hsc_meta_cols
)

if (
  length(
    missing_hsc_cols
  ) > 0
) {
  stop(
    "Missing HSC metadata: ",
    paste(
      missing_hsc_cols,
      collapse = ";"
    )
  )
}

hmd_all <- hsc[[]] %>%
  rownames_to_column(
    "cell"
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[HSC_SAMPLE_COL]]
      ),
    state2 =
      as.character(
        .data[[HSC_STATE2_COL]]
      ),
    state3 =
      as.character(
        .data[[HSC_STATE3_COL]]
      ),
    qc_primary_raw =
      .data[[HSC_QC_COL]]
  ) %>%
  filter(
    sample %in%
      SAMPLES
  ) %>%
  mutate(
    group =
      unname(
        sample_to_group[sample]
      )
  )

qc_to_logical <- function(x) {

  if (is.logical(x)) {
    return(
      x
    )
  }

  if (is.numeric(x)) {
    return(
      !is.na(x) &
        x != 0
    )
  }

  x2 <- tolower(
    trimws(
      as.character(x)
    )
  )

  x2 %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y",
    "include",
    "included",
    "primary",
    "pass"
  )
}

hmd_all$qc_primary_include_resolved <-
  qc_to_logical(
    hmd_all$qc_primary_raw
  )

hmd <- hmd_all %>%
  filter(
    qc_primary_include_resolved
  )

if (
  nrow(hmd) == 0
) {
  stop(
    "Primary-clean HSC filter returned zero cells."
  )
}

hsc_primary <- subset(
  hsc,
  cells = hmd$cell
)

hsc_inventory <- hmd %>%
  count(
    sample,
    group,
    state2,
    state3,
    name = "cells"
  ) %>%
  arrange(
    sample,
    state2,
    state3
  )

write_tsv(
  hsc_inventory,
  file.path(
    DIR_TABLE,
    "06_HSC_primary_clean_inventory.tsv"
  )
)

# ------------------------------------------------------------
# 10. HSC state fractions
# ------------------------------------------------------------

state2_frac <- hmd %>%
  count(
    sample,
    group,
    state2,
    name = "cells"
  ) %>%
  group_by(
    sample,
    group
  ) %>%
  mutate(
    fraction =
      cells /
      sum(cells)
  ) %>%
  ungroup()

state3_frac <- hmd %>%
  count(
    sample,
    group,
    state3,
    name = "cells"
  ) %>%
  group_by(
    sample,
    group
  ) %>%
  mutate(
    fraction =
      cells /
      sum(cells)
  ) %>%
  ungroup()

write_tsv(
  state2_frac,
  file.path(
    DIR_TABLE,
    "07_HSC_state2_fraction_by_sample.tsv"
  )
)

write_tsv(
  state3_frac,
  file.path(
    DIR_TABLE,
    "08_HSC_state3_fraction_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 11. Whole-HSC and aHSC pseudobulk
# ------------------------------------------------------------

hsc_counts <- GetAssayData(
  hsc_primary,
  assay = "RNA",
  layer = "counts"
)

score_gene_set_pb <- function(
  counts_mat,
  cell_meta,
  genes,
  samples_use
) {

  available <- intersect(
    genes,
    rownames(
      counts_mat
    )
  )

  if (
    length(
      available
    ) < 3
  ) {
    stop(
      "Insufficient pseudobulk genes: ",
      paste(
        genes,
        collapse = ";"
      )
    )
  }

  pb_list <- lapply(
    samples_use,
    function(s) {

      cells_s <-
        cell_meta$cell[
          cell_meta$sample ==
            s
        ]

      cells_s <- intersect(
        cells_s,
        colnames(
          counts_mat
        )
      )

      if (
        length(
          cells_s
        ) == 0
      ) {
        return(
          rep(
            0,
            nrow(
              counts_mat
            )
          )
        )
      }

      Matrix::rowSums(
        counts_mat[
          ,
          cells_s,
          drop = FALSE
        ]
      )
    }
  )

  pb_counts <- do.call(
    cbind,
    pb_list
  )

  rownames(
    pb_counts
  ) <- rownames(
    counts_mat
  )

  colnames(
    pb_counts
  ) <- samples_use

  lib <- colSums(
    pb_counts
  )

  logcpm <- log2(
    sweep(
      pb_counts,
      2,
      lib,
      "/"
    ) *
      1e6 +
      1
  )

  sample_x_gene <- t(
    logcpm[
      available,
      ,
      drop = FALSE
    ]
  )

  z_fun <- function(x) {

    sx <- sd(
      x,
      na.rm = TRUE
    )

    if (
      !is.finite(
        sx
      ) ||
      sx == 0
    ) {
      return(
        rep(
          0,
          length(
            x
          )
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

  z <- apply(
    sample_x_gene,
    2,
    z_fun
  )

  if (
    is.vector(
      z
    )
  ) {
    z <- matrix(
      z,
      ncol = 1
    )
  }

  tibble(
    sample =
      samples_use,

    score =
      rowMeans(
        z,
        na.rm = TRUE
      ),

    genes_available_n =
      length(
        available
      ),

    genes_available =
      paste(
        available,
        collapse = ";"
      )
  )
}

whole_IEG <- score_gene_set_pb(
  hsc_counts,
  hmd,
  IEG_GENES,
  SAMPLES
) %>%
  rename(
    HSC_whole_IEG_PB =
      score
  )

whole_ECM <- score_gene_set_pb(
  hsc_counts,
  hmd,
  ECM_GENES,
  SAMPLES
) %>%
  rename(
    HSC_whole_ECM_PB =
      score
  )

aHSC_md <- hmd %>%
  filter(
    state2 ==
      "aHSC"
  )

aHSC_IEG <- score_gene_set_pb(
  hsc_counts,
  aHSC_md,
  IEG_GENES,
  SAMPLES
) %>%
  rename(
    aHSC_IEG_PB =
      score
  )

aHSC_ECM <- score_gene_set_pb(
  hsc_counts,
  aHSC_md,
  ECM_GENES,
  SAMPLES
) %>%
  rename(
    aHSC_ECM_PB =
      score
  )

hsc_pb <- tibble(
  sample =
    SAMPLES,
  group =
    unname(
      sample_to_group[
        SAMPLES
      ]
    )
) %>%
  left_join(
    whole_IEG %>%
      select(
        sample,
        HSC_whole_IEG_PB
      ),
    by = "sample"
  ) %>%
  left_join(
    whole_ECM %>%
      select(
        sample,
        HSC_whole_ECM_PB
      ),
    by = "sample"
  ) %>%
  left_join(
    aHSC_IEG %>%
      select(
        sample,
        aHSC_IEG_PB
      ),
    by = "sample"
  ) %>%
  left_join(
    aHSC_ECM %>%
      select(
        sample,
        aHSC_ECM_PB
      ),
    by = "sample"
  )

write_tsv(
  hsc_pb,
  file.path(
    DIR_TABLE,
    "09_HSC_pseudobulk_functional_scores.tsv"
  )
)

# ------------------------------------------------------------
# 12. HSC endpoint table
# ------------------------------------------------------------

state2_wide <- state2_frac %>%
  select(
    sample,
    group,
    state2,
    fraction
  ) %>%
  pivot_wider(
    names_from =
      state2,
    values_from =
      fraction,
    names_prefix =
      "state2_fraction_",
    values_fill = 0
  )

state3_wide <- state3_frac %>%
  select(
    sample,
    group,
    state3,
    fraction
  ) %>%
  pivot_wider(
    names_from =
      state3,
    values_from =
      fraction,
    names_prefix =
      "state3_fraction_",
    values_fill = 0
  )

hsc_endpoints <- hsc_pb %>%
  left_join(
    state2_wide,
    by = c(
      "sample",
      "group"
    )
  ) %>%
  left_join(
    state3_wide,
    by = c(
      "sample",
      "group"
    )
  )

write_tsv(
  hsc_endpoints,
  file.path(
    DIR_TABLE,
    "10_PRIMARY_HSC_endpoints_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 13. Descriptive HSC Sham vs Tx effects
# ------------------------------------------------------------

numeric_endpoints <- setdiff(
  colnames(
    hsc_endpoints
  ),
  c(
    "sample",
    "group"
  )
)

hsc_group_effect <- hsc_endpoints %>%
  pivot_longer(
    cols =
      all_of(
        numeric_endpoints
      ),
    names_to = "endpoint",
    values_to = "value"
  ) %>%
  group_by(
    endpoint,
    group
  ) %>%
  summarise(
    n_samples = n(),
    mean_value =
      mean(
        value,
        na.rm = TRUE
      ),
    sample1 =
      value[1],
    sample2 =
      value[
        min(
          2,
          length(
            value
          )
        )
      ],
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      n_samples,
      mean_value,
      sample1,
      sample2
    )
  ) %>%
  mutate(
    Tx_minus_Sham =
      mean_value_Tx -
      mean_value_Sham
  )

write_tsv(
  hsc_group_effect,
  file.path(
    DIR_TABLE,
    "11_PRIMARY_HSC_Sham_vs_Tx_descriptive.tsv"
  )
)

# ------------------------------------------------------------
# 14. HSC receptor support
# ------------------------------------------------------------

focus_receptors <- unique(
  unlist(
    RECEPTOR_MAP
  )
)

receptor_availability <- tibble(
  receptor =
    focus_receptors,
  measurable =
    focus_receptors %in%
    rownames(
      hsc_primary
    )
)

receptors_measurable <- receptor_availability %>%
  filter(
    measurable
  ) %>%
  pull(
    receptor
  )

receptor_counts <- if (
  length(
    receptors_measurable
  ) > 0
) {
  FetchData(
    hsc_primary,
    vars =
      receptors_measurable,
    layer = "counts"
  ) %>%
    rownames_to_column(
      "cell"
    )
} else {
  tibble(
    cell =
      character()
  )
}

if (
  length(
    receptors_measurable
  ) > 0
) {

  receptor_long <- receptor_counts %>%
    pivot_longer(
      -cell,
      names_to =
        "receptor",
      values_to =
        "count"
    ) %>%
    left_join(
      hmd %>%
        select(
          cell,
          sample,
          group
        ),
      by = "cell"
    )

  receptor_sample <- receptor_long %>%
    group_by(
      sample,
      group,
      receptor
    ) %>%
    summarise(
      HSC_cells =
        n(),
      detection_fraction =
        mean(
          count > 0,
          na.rm = TRUE
        ),
      .groups = "drop"
    )

} else {

  receptor_sample <- tibble(
    sample =
      character(),
    group =
      character(),
    receptor =
      character(),
    HSC_cells =
      integer(),
    detection_fraction =
      numeric()
  )
}

receptor_ligand_map <- bind_rows(
  lapply(
    names(
      RECEPTOR_MAP
    ),
    function(lig) {
      tibble(
        ligand = lig,
        receptor =
          RECEPTOR_MAP[[lig]]
      )
    }
  )
)

receptor_support <- receptor_ligand_map %>%
  left_join(
    receptor_availability,
    by = "receptor"
  ) %>%
  left_join(
    receptor_sample,
    by = "receptor"
  ) %>%
  arrange(
    ligand,
    receptor,
    sample
  )

write_tsv(
  receptor_support,
  file.path(
    DIR_TABLE,
    "12_HSC_receptor_support_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 15. Descriptive sender-receiver correlations across 4 samples
# ------------------------------------------------------------

trem2_primary <- trem2_burden %>%
  select(
    sample,
    group,
    ligand,
    ligand_burden_expr,
    ligand_detection,
    ligand_expression_share,
    state_fraction
  )

corr_input <- trem2_primary %>%
  left_join(
    hsc_endpoints,
    by = c(
      "sample",
      "group"
    )
  )

HSC_CORR_ENDPOINTS <- intersect(
  c(
    "HSC_whole_IEG_PB",
    "HSC_whole_ECM_PB",
    "aHSC_IEG_PB",
    "aHSC_ECM_PB",
    "state2_fraction_aHSC",
    "state2_fraction_qHSC",
    "state3_fraction_ECM-activated HSC",
    "state3_fraction_Contractile HSC"
  ),
  colnames(
    corr_input
  )
)

safe_rho <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  if (
    sum(
      keep
    ) < 3 ||
    length(
      unique(
        x[keep]
      )
    ) < 2 ||
    length(
      unique(
        y[keep]
      )
    ) < 2
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = "spearman"
    )
  )
}

sender_receiver_corr <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      d <- corr_input %>%
        filter(
          ligand ==
            lig
        ) %>%
        arrange(
          sample
        )

      bind_rows(
        lapply(
          HSC_CORR_ENDPOINTS,
          function(ep) {

            tibble(
              ligand = lig,
              sender_metric =
                "TREM2_ligand_burden_expr",
              HSC_endpoint =
                ep,
              n_samples =
                sum(
                  complete.cases(
                    d$ligand_burden_expr,
                    d[[ep]]
                  )
                ),
              rho_descriptive =
                safe_rho(
                  d$ligand_burden_expr,
                  d[[ep]]
                ),
              inference =
                "descriptive_only_n4"
            )
          }
        )
      )
    }
  )
) %>%
  arrange(
    ligand,
    desc(
      abs(
        rho_descriptive
      )
    )
  )

write_tsv(
  sender_receiver_corr,
  file.path(
    DIR_TABLE,
    "13_sender_receiver_correlations_descriptive_n4.tsv"
  )
)

# ------------------------------------------------------------
# 16. Cross-species expected-direction summary
# ------------------------------------------------------------

ligand_group_effect <- sender_group_summary %>%
  filter(
    metric ==
      "ligand_burden_expr"
  ) %>%
  select(
    ligand,
    Sham_mean =
      mean_value_Sham,
    Tx_mean =
      mean_value_Tx,
    Tx_minus_Sham,
    log2_Tx_over_Sham
  ) %>%
  mutate(
    expected_with_effective_Tx =
      "decrease",
    observed_direction =
      case_when(
        Tx_minus_Sham < 0 ~
          "decrease",
        Tx_minus_Sham > 0 ~
          "increase",
        TRUE ~
          "no_change"
      ),
    direction_support =
      observed_direction ==
      expected_with_effective_Tx
  )

write_tsv(
  ligand_group_effect,
  file.path(
    DIR_TABLE,
    "14_PRIMARY_cross_species_ligand_direction_summary.tsv"
  )
)

# ------------------------------------------------------------
# 17. Figures
# ------------------------------------------------------------

plot_ligand <- trem2_burden %>%
  filter(
    ligand %in%
      FOCUS_LIGANDS
  )

p1 <- ggplot(
  plot_ligand,
  aes(
    x = group,
    y = ligand_burden_expr
  )
) +
  geom_point(
    size = 3,
    position =
      position_jitter(
        width = 0.05,
        height = 0
      )
  ) +
  geom_text(
    aes(
      label = sample
    ),
    nudge_y = 0.01,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ ligand,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "TREM2-Mphi ligand burden: Sham vs Tx",
    x = NULL,
    y =
      "Abundance-weighted ligand burden"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TREM2_focus_ligand_burden_Sham_vs_Tx.pdf"
  ),
  p1,
  width = 9,
  height = 4.5
)

pb_plot <- hsc_endpoints %>%
  select(
    sample,
    group,
    HSC_whole_IEG_PB,
    HSC_whole_ECM_PB,
    aHSC_IEG_PB,
    aHSC_ECM_PB
  ) %>%
  pivot_longer(
    cols = -c(
      sample,
      group
    ),
    names_to =
      "endpoint",
    values_to =
      "score"
  )

p2 <- ggplot(
  pb_plot,
  aes(
    x = group,
    y = score
  )
) +
  geom_point(
    size = 3,
    position =
      position_jitter(
        width = 0.05,
        height = 0
      )
  ) +
  geom_text(
    aes(
      label = sample
    ),
    nudge_y = 0.05,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ endpoint,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "HSC pseudobulk functional response: Sham vs Tx",
    x = NULL,
    y =
      "Gene-wise z-aggregated pseudobulk score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_HSC_functional_pseudobulk_Sham_vs_Tx.pdf"
  ),
  p2,
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# 18. Save compact result object
# ------------------------------------------------------------

saveRDS(
  list(
    mphi_sample_inventory =
      mphi_sample_inventory,
    burden =
      burden,
    trem2_burden =
      trem2_burden,
    sender_group_summary =
      sender_group_summary,
    hsc_inventory =
      hsc_inventory,
    state2_frac =
      state2_frac,
    state3_frac =
      state3_frac,
    hsc_endpoints =
      hsc_endpoints,
    hsc_group_effect =
      hsc_group_effect,
    receptor_support =
      receptor_support,
    sender_receiver_corr =
      sender_receiver_corr,
    ligand_group_effect =
      ligand_group_effect
  ),
  file.path(
    DIR_OBJ,
    "Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1_results.rds"
  )
)

# ------------------------------------------------------------
# 19. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "Mphi_source_cells",
    "Mphi_primary_states_n",
    "Mphi_samples_n",
    "focus_ligands_n",
    "HSC_source_cells",
    "HSC_primary_clean_cells",
    "HSC_samples_n",
    "HSC_state2_labels",
    "HSC_state3_labels",
    "group_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      ncol(
        mphi
      )
    ),
    as.character(
      length(
        unique(
          mmd_primary2$sender_state
        )
      )
    ),
    as.character(
      length(
        unique(
          mmd_primary2$sample
        )
      )
    ),
    as.character(
      length(
        FOCUS_LIGANDS
      )
    ),
    as.character(
      ncol(
        hsc
      )
    ),
    as.character(
      nrow(
        hmd
      )
    ),
    as.character(
      length(
        unique(
          hmd$sample
        )
      )
    ),
    paste(
      sort(
        unique(
          hmd$state2
        )
      ),
      collapse = ";"
    ),
    paste(
      sort(
        unique(
          hmd$state3
        )
      ),
      collapse = ";"
    ),
    "descriptive_only_n2_per_group",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v6.4.1.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "Mphi_source_RDS",
    "HSC_source_RDS",
    "Mphi_state_column",
    "Mphi_sample_column",
    "HSC_sample_column",
    "HSC_state2_column",
    "HSC_state3_column",
    "HSC_QC_column",
    "primary_sender_state",
    "focus_ligands",
    "Mphi_burden_denominator",
    "HSC_functional_scoring",
    "group_comparison",
    "cross_sample_correlations"
  ),
  value = c(
    "v6.4.1",
    MPHI_RDS,
    HSC_RDS,
    mphi_state_col,
    mphi_sample_col,
    HSC_SAMPLE_COL,
    HSC_STATE2_COL,
    HSC_STATE3_COL,
    HSC_QC_COL,
    TREM2_STATE,
    paste(
      FOCUS_LIGANDS,
      collapse = ";"
    ),
    "five primary Mphi states; Other excluded",
    "whole-HSC and aHSC pseudobulk logCPM, gene-wise z across 4 biological samples",
    "descriptive Sham n=2 vs Tx n=2; no primary inferential P-values",
    "Spearman rho descriptive only, n=4"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "Mouse_HumanAxis_validation_provenance_v6.4.1.tsv"
  )
)

# ------------------------------------------------------------
# 20. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Mphi SAMPLE / STATE INVENTORY\n")
cat("============================================================\n")
print(
  mphi_sample_inventory
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY TREM2 LIGAND BURDEN\n")
cat("============================================================\n")
print(
  trem2_burden
)

cat("\n")
cat("============================================================\n")
cat("TREM2 LIGAND SHAM vs TX\n")
cat("============================================================\n")
print(
  sender_group_summary
)

cat("\n")
cat("============================================================\n")
cat("HSC ENDPOINTS\n")
cat("============================================================\n")
print(
  hsc_endpoints
)

cat("\n")
cat("============================================================\n")
cat("HSC SHAM vs TX\n")
cat("============================================================\n")
print(
  hsc_group_effect
)

cat("\n")
cat("============================================================\n")
cat("HSC RECEPTOR SUPPORT\n")
cat("============================================================\n")
print(
  receptor_support
)

cat("\n")
cat("============================================================\n")
cat("SENDER-RECEIVER CORRELATIONS (DESCRIPTIVE n=4)\n")
cat("============================================================\n")
print(
  sender_receiver_corr
)

cat("\n")
cat("============================================================\n")
cat("CROSS-SPECIES LIGAND DIRECTION SUMMARY\n")
cat("============================================================\n")
print(
  ligand_group_effect
)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(
  type = "message"
)

sink(
  type = "output"
)

close(zz)
