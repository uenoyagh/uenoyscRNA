#!/usr/bin/env Rscript

# ============================================================
# Mouse MASH HSC 3-axis + total-host-weighted burden v6.4.4
#
# PURPOSE
#   Extend v6.4.3 by separating HSC transcriptional programs into:
#
#     1) Core collagen / matrix-production
#     2) ECM remodeling / fibrosis-response
#     3) IEG-core / immediate-early response
#
#   Then integrate three levels:
#
#     A. HSC abundance among recovered Mouse host cells
#     B. HSC state composition within HSC
#     C. State-specific expression intensity
#
#   Primary host-weighted burden:
#
#     HSC abundance
#       x state fraction
#       x mean linear CP10k program intensity
#
# IMPORTANT INTERPRETATION
#   This is a scRNA-recovered host-cell-weighted transcript burden.
#   It is NOT an absolute tissue collagen-production measurement.
#
# STATISTICS
#   Biological sample is the unit.
#   Sham n=2, Tx n=2.
#   Group comparisons are descriptive.
#
# HSC PRIMARY PARTITION
#   qHSC
#   ECM-activated HSC
#   Contractile HSC
#
# HSC SENSITIVITY PARTITION
#   qHSC
#   aHSC
#
# No reclustering.
# No source-RDS overwrite.
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
# 0. Paths / constants
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

WHOLE_FILENAME <-
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"

WHOLE_CANDIDATES <- c(
  file.path(
    ROOT,
    "rds",
    WHOLE_FILENAME
  ),
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    WHOLE_FILENAME
  ),
  file.path(
    ROOT,
    WHOLE_FILENAME
  )
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4"
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
  "Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4.log"
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
cat("Mouse MASH HSC 3-axis + total-host-weighted burden v6.4.4\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

EXPECTED_HSC_PRIMARY <- 12941L
EXPECTED_WHOLE_CELLS <- 104588L

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

HSC_SAMPLE_COL <- "sample_hsc_v610"
HSC_STATE2_COL <- "HSC_state2_v612"
HSC_STATE3_COL <- "HSC_state3_v612"
HSC_QC_COL <- "qc_primary_include"

STATE3_EXPECTED <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

STATE2_EXPECTED <- c(
  "qHSC",
  "aHSC"
)

# ------------------------------------------------------------
# 1. Primary three-axis definitions
#
# These are intentionally non-overlapping.
# ------------------------------------------------------------

CORE_COLLAGEN <- c(
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Col5a1",
  "Col5a2"
)

ECM_REMODELING <- c(
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

IEG_CORE <- c(
  "Junb",
  "Fos",
  "Fosb"
)

PROGRAMS <- list(
  Core_collagen =
    CORE_COLLAGEN,
  ECM_remodeling =
    ECM_REMODELING,
  IEG_core =
    IEG_CORE
)

# Previous broader IEG panel retained only as a sensitivity audit.
IEG_EXTENDED_SENSITIVITY <- c(
  "Thbs1",
  "Timp1",
  "Junb",
  "Fos",
  "Fosb"
)

# ------------------------------------------------------------
# 2. Resolve whole-Mouse source RDS safely
# ------------------------------------------------------------

WHOLE_RDS <- WHOLE_CANDIDATES[
  file.exists(
    WHOLE_CANDIDATES
  )
][1]

if (
  length(WHOLE_RDS) == 0 ||
  is.na(WHOLE_RDS)
) {

  recursive_hits <- list.files(
    ROOT,
    pattern =
      paste0(
        "^",
        WHOLE_FILENAME,
        "$"
      ),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(recursive_hits) == 1) {

    WHOLE_RDS <- recursive_hits

  } else if (length(recursive_hits) == 0) {

    stop(
      "Whole-Mouse RDS not found: ",
      WHOLE_FILENAME
    )

  } else {

    stop(
      "Multiple whole-Mouse RDS matches found:\n",
      paste(
        recursive_hits,
        collapse = "\n"
      )
    )
  }
}

stopifnot(
  file.exists(HSC_RDS),
  file.exists(WHOLE_RDS)
)

cat("HSC RDS:\n", HSC_RDS, "\n\n")
cat("Whole-Mouse RDS:\n", WHOLE_RDS, "\n\n")

# ------------------------------------------------------------
# 3. Load objects
# ------------------------------------------------------------

hsc <- readRDS(HSC_RDS)
whole <- readRDS(WHOLE_RDS)

DefaultAssay(hsc) <- "RNA"

if (ncol(hsc) != EXPECTED_HSC_PRIMARY) {
  stop(
    "Unexpected HSC cell count: expected ",
    EXPECTED_HSC_PRIMARY,
    ", observed ",
    ncol(hsc)
  )
}

if (ncol(whole) != EXPECTED_WHOLE_CELLS) {
  stop(
    "Unexpected whole-Mouse cell count: expected ",
    EXPECTED_WHOLE_CELLS,
    ", observed ",
    ncol(whole)
  )
}

# ------------------------------------------------------------
# 4. Resolve whole-object sample column
# ------------------------------------------------------------

resolve_sample_column <- function(
  meta,
  required_samples
) {

  preferred <- c(
    "sample_4group",
    "sample_FIXED2",
    "sample",
    "sample_for_annotation",
    "orig.ident"
  )

  for (cc in preferred) {

    if (cc %in% colnames(meta)) {

      vals <- unique(
        as.character(
          meta[[cc]]
        )
      )

      if (all(
        required_samples %in%
          vals
      )) {
        return(cc)
      }
    }
  }

  matching <- colnames(meta)[
    vapply(
      meta,
      function(x) {

        vals <- unique(
          as.character(x)
        )

        all(
          required_samples %in%
            vals
        )
      },
      logical(1)
    )
  ]

  if (length(matching) == 1) {
    return(matching)
  }

  if (length(matching) == 0) {
    stop(
      "Could not resolve whole-object sample column."
    )
  }

  stop(
    "Multiple possible sample columns in whole object: ",
    paste(
      matching,
      collapse = ";"
    )
  )
}

whole_sample_col <- resolve_sample_column(
  whole[[]],
  SAMPLES
)

cat(
  "Resolved whole-object sample column:",
  whole_sample_col,
  "\n\n"
)

# ------------------------------------------------------------
# 5. HSC metadata / QC validation
# ------------------------------------------------------------

required_hsc_meta <- c(
  HSC_SAMPLE_COL,
  HSC_STATE2_COL,
  HSC_STATE3_COL,
  HSC_QC_COL
)

missing_hsc_meta <- setdiff(
  required_hsc_meta,
  colnames(hsc[[]])
)

if (length(missing_hsc_meta) > 0) {
  stop(
    "Missing HSC metadata: ",
    paste(
      missing_hsc_meta,
      collapse = ";"
    )
  )
}

qc_raw <- hsc[[]][[HSC_QC_COL]]

qc_logical <- if (is.logical(qc_raw)) {

  qc_raw

} else if (is.numeric(qc_raw)) {

  !is.na(qc_raw) &
    qc_raw != 0

} else {

  tolower(
    trimws(
      as.character(qc_raw)
    )
  ) %in% c(
    "true", "t", "1",
    "yes", "y",
    "include", "included",
    "primary", "pass"
  )
}

if (!all(qc_logical)) {
  stop(
    "Primary-Clean HSC RDS contains qc_primary_include != TRUE."
  )
}

hmd <- hsc[[]] %>%
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

if (!setequal(
  sort(
    unique(
      hmd$state3
    )
  ),
  sort(
    STATE3_EXPECTED
  )
)) {
  stop(
    "Unexpected state3 labels: ",
    paste(
      sort(
        unique(
          hmd$state3
        )
      ),
      collapse = ";"
    )
  )
}

if (!setequal(
  sort(
    unique(
      hmd$state2
    )
  ),
  sort(
    STATE2_EXPECTED
  )
)) {
  stop(
    "Unexpected state2 labels: ",
    paste(
      sort(
        unique(
          hmd$state2
        )
      ),
      collapse = ";"
    )
  )
}

# ------------------------------------------------------------
# 6. Whole-host denominator inventory
#
# The source object is explicitly the Mouse object.
# We call this "recovered host-cell denominator".
# ------------------------------------------------------------

wmd <- whole[[]] %>%
  rownames_to_column(
    "cell"
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[whole_sample_col]]
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

if (!all(
  SAMPLES %in%
    unique(
      wmd$sample
    )
)) {
  stop(
    "Not all four Sham/Tx samples found in whole-Mouse object."
  )
}

# Strong provenance check: every Primary-Clean HSC barcode should
# be present in the whole Mouse object.
hsc_cell_membership_fraction <-
  mean(
    hmd$cell %in%
      wmd$cell
  )

if (
  hsc_cell_membership_fraction <
    1
) {
  stop(
    "Not all Primary-Clean HSC cells are present in whole-Mouse object. ",
    "Membership fraction = ",
    hsc_cell_membership_fraction
  )
}

whole_counts <- wmd %>%
  count(
    sample,
    group,
    name =
      "recovered_host_cells"
  )

hsc_counts <- hmd %>%
  count(
    sample,
    group,
    name =
      "primary_clean_HSC_cells"
  )

hsc_abundance <- whole_counts %>%
  left_join(
    hsc_counts,
    by = c(
      "sample",
      "group"
    )
  ) %>%
  mutate(
    primary_clean_HSC_cells =
      replace_na(
        primary_clean_HSC_cells,
        0L
      ),

    HSC_fraction_of_recovered_host =
      primary_clean_HSC_cells /
      recovered_host_cells
  ) %>%
  arrange(
    match(
      sample,
      SAMPLES
    )
  )

write_tsv(
  hsc_abundance,
  file.path(
    DIR_TABLE,
    "01_HSC_abundance_among_recovered_host_cells.tsv"
  )
)

# ------------------------------------------------------------
# 7. Program availability
# ------------------------------------------------------------

program_audit <- bind_rows(
  lapply(
    names(PROGRAMS),
    function(p) {

      genes <-
        PROGRAMS[[p]]

      tibble(
        program = p,
        gene = genes,
        measurable =
          genes %in%
          rownames(hsc)
      )
    }
  )
)

write_tsv(
  program_audit,
  file.path(
    DIR_TABLE,
    "02_primary_3axis_gene_availability.tsv"
  )
)

program_coverage <- program_audit %>%
  group_by(program) %>%
  summarise(
    requested_n = n(),
    available_n =
      sum(measurable),
    coverage_fraction =
      mean(measurable),
    available_genes =
      paste(
        gene[measurable],
        collapse = ";"
      ),
    missing_genes =
      paste(
        gene[!measurable],
        collapse = ";"
      ),
    .groups = "drop"
  )

write_tsv(
  program_coverage,
  file.path(
    DIR_TABLE,
    "03_primary_3axis_coverage_summary.tsv"
  )
)

PROGRAMS_AVAIL <- lapply(
  PROGRAMS,
  function(g) {
    intersect(
      g,
      rownames(hsc)
    )
  }
)

if (
  length(PROGRAMS_AVAIL$Core_collagen) <
    4 ||
  length(PROGRAMS_AVAIL$ECM_remodeling) <
    6 ||
  length(PROGRAMS_AVAIL$IEG_core) <
    2
) {
  stop(
    "Insufficient primary 3-axis gene coverage."
  )
}

all_primary_genes <- unique(
  unlist(
    PROGRAMS_AVAIL
  )
)

# ------------------------------------------------------------
# 8. Linear CP10k per HSC cell
# ------------------------------------------------------------

hsc_counts_matrix <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

libsize <- Matrix::colSums(
  hsc_counts_matrix[
    ,
    hmd$cell,
    drop = FALSE
  ]
)

if (any(
  !is.finite(libsize) |
    libsize <= 0
)) {
  stop(
    "Invalid HSC per-cell library size."
  )
}

selected_counts <- hsc_counts_matrix[
  all_primary_genes,
  hmd$cell,
  drop = FALSE
]

cp10k <- t(
  t(selected_counts) /
    libsize
) * 1e4

program_score <- function(
  mat,
  genes
) {

  g <- intersect(
    genes,
    rownames(mat)
  )

  Matrix::colMeans(
    mat[
      g,
      ,
      drop = FALSE
    ]
  )
}

cell_program <- hmd %>%
  mutate(
    Core_collagen =
      as.numeric(
        program_score(
          cp10k,
          PROGRAMS_AVAIL$Core_collagen
        )
      ),

    ECM_remodeling =
      as.numeric(
        program_score(
          cp10k,
          PROGRAMS_AVAIL$ECM_remodeling
        )
      ),

    IEG_core =
      as.numeric(
        program_score(
          cp10k,
          PROGRAMS_AVAIL$IEG_core
        )
      )
  )

# Sensitivity: prior broader IEG panel, if measurable.
ieg_ext_avail <- intersect(
  IEG_EXTENDED_SENSITIVITY,
  rownames(hsc)
)

if (length(ieg_ext_avail) >= 3) {

  ext_counts <- hsc_counts_matrix[
    ieg_ext_avail,
    hmd$cell,
    drop = FALSE
  ]

  ext_cp10k <- t(
    t(ext_counts) /
      libsize
  ) * 1e4

  cell_program$IEG_extended_sensitivity <-
    as.numeric(
      Matrix::colMeans(
        ext_cp10k
      )
    )

} else {

  cell_program$IEG_extended_sensitivity <-
    NA_real_
}

# ------------------------------------------------------------
# 9. Generic partition summarizer
# ------------------------------------------------------------

summarise_partition <- function(
  df,
  state_col,
  partition_label
) {

  total_hsc <- df %>%
    count(
      sample,
      group,
      name =
        "total_HSC_n"
    )

  state_n <- df %>%
    count(
      sample,
      group,
      state =
        .data[[state_col]],
      name =
        "state_n"
    ) %>%
    left_join(
      total_hsc,
      by = c(
        "sample",
        "group"
      )
    ) %>%
    mutate(
      state_fraction_within_HSC =
        state_n /
        total_HSC_n
    )

  program_long <- df %>%
    select(
      sample,
      group,
      state =
        all_of(state_col),
      all_of(
        names(PROGRAMS)
      )
    ) %>%
    pivot_longer(
      cols =
        all_of(
          names(PROGRAMS)
        ),
      names_to =
        "program",
      values_to =
        "cell_program_CP10k"
    ) %>%
    group_by(
      sample,
      group,
      state,
      program
    ) %>%
    summarise(
      mean_intensity_CP10k =
        mean(
          cell_program_CP10k,
          na.rm = TRUE
        ),
      median_intensity_CP10k =
        median(
          cell_program_CP10k,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    left_join(
      state_n,
      by = c(
        "sample",
        "group",
        "state"
      )
    ) %>%
    left_join(
      hsc_abundance %>%
        select(
          sample,
          group,
          recovered_host_cells,
          primary_clean_HSC_cells,
          HSC_fraction_of_recovered_host
        ),
      by = c(
        "sample",
        "group"
      )
    ) %>%
    mutate(
      partition =
        partition_label,

      within_HSC_weighted_burden =
        state_fraction_within_HSC *
        mean_intensity_CP10k,

      state_fraction_of_recovered_host =
        HSC_fraction_of_recovered_host *
        state_fraction_within_HSC,

      host_weighted_burden =
        HSC_fraction_of_recovered_host *
        state_fraction_within_HSC *
        mean_intensity_CP10k,
      .before = 1
    )

  program_long
}

# The summarizer works on cell_program where columns are named
# state2/state3, not the original metadata names.
state3 <- summarise_partition(
  cell_program,
  "state3",
  "state3"
)

state2 <- summarise_partition(
  cell_program,
  "state2",
  "state2"
)

write_tsv(
  state3,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_state3_3axis_host_weighted_burden.tsv"
  )
)

write_tsv(
  state2,
  file.path(
    DIR_TABLE,
    "05_SENSITIVITY_state2_3axis_host_weighted_burden.tsv"
  )
)

# ------------------------------------------------------------
# 10. Whole-HSC and whole-host reconstructed burden
# ------------------------------------------------------------

whole_hsc_axis <- cell_program %>%
  group_by(
    sample,
    group
  ) %>%
  summarise(
    Core_collagen =
      mean(
        Core_collagen,
        na.rm = TRUE
      ),
    ECM_remodeling =
      mean(
        ECM_remodeling,
        na.rm = TRUE
      ),
    IEG_core =
      mean(
        IEG_core,
        na.rm = TRUE
      ),
    IEG_extended_sensitivity =
      mean(
        IEG_extended_sensitivity,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  left_join(
    hsc_abundance,
    by = c(
      "sample",
      "group"
    )
  ) %>%
  pivot_longer(
    cols = c(
      Core_collagen,
      ECM_remodeling,
      IEG_core
    ),
    names_to =
      "program",
    values_to =
      "whole_HSC_mean_intensity_CP10k"
  ) %>%
  mutate(
    whole_host_weighted_HSC_burden =
      HSC_fraction_of_recovered_host *
      whole_HSC_mean_intensity_CP10k
  )

state3_recon <- state3 %>%
  group_by(
    sample,
    group,
    program
  ) %>%
  summarise(
    reconstructed_within_HSC =
      sum(
        within_HSC_weighted_burden,
        na.rm = TRUE
      ),
    reconstructed_host_weighted =
      sum(
        host_weighted_burden,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

recon_check <- whole_hsc_axis %>%
  left_join(
    state3_recon,
    by = c(
      "sample",
      "group",
      "program"
    )
  ) %>%
  mutate(
    within_HSC_abs_error =
      abs(
        whole_HSC_mean_intensity_CP10k -
          reconstructed_within_HSC
      ),

    host_weighted_abs_error =
      abs(
        whole_host_weighted_HSC_burden -
          reconstructed_host_weighted
      )
  )

write_tsv(
  recon_check,
  file.path(
    DIR_TABLE,
    "06_state3_reconstruction_check.tsv"
  )
)

if (
  max(
    recon_check$within_HSC_abs_error,
    na.rm = TRUE
  ) > 1e-8 ||
  max(
    recon_check$host_weighted_abs_error,
    na.rm = TRUE
  ) > 1e-8
) {
  stop(
    "State3 burden reconstruction failed."
  )
}

# ------------------------------------------------------------
# 11. Group summaries
# ------------------------------------------------------------

group_state3 <- state3 %>%
  group_by(
    state,
    program,
    group
  ) %>%
  summarise(
    n_samples = n(),

    mean_HSC_abundance =
      mean(
        HSC_fraction_of_recovered_host
      ),

    mean_state_fraction =
      mean(
        state_fraction_within_HSC
      ),

    mean_intensity =
      mean(
        mean_intensity_CP10k
      ),

    mean_within_HSC_burden =
      mean(
        within_HSC_weighted_burden
      ),

    mean_host_weighted_burden =
      mean(
        host_weighted_burden
      ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      n_samples,
      mean_HSC_abundance,
      mean_state_fraction,
      mean_intensity,
      mean_within_HSC_burden,
      mean_host_weighted_burden
    )
  ) %>%
  mutate(
    delta_HSC_abundance =
      mean_HSC_abundance_Tx -
      mean_HSC_abundance_Sham,

    delta_state_fraction =
      mean_state_fraction_Tx -
      mean_state_fraction_Sham,

    delta_intensity =
      mean_intensity_Tx -
      mean_intensity_Sham,

    delta_within_HSC_burden =
      mean_within_HSC_burden_Tx -
      mean_within_HSC_burden_Sham,

    delta_host_weighted_burden =
      mean_host_weighted_burden_Tx -
      mean_host_weighted_burden_Sham,

    ratio_host_weighted_Tx_over_Sham =
      ifelse(
        mean_host_weighted_burden_Sham > 0,
        mean_host_weighted_burden_Tx /
          mean_host_weighted_burden_Sham,
        NA_real_
      )
  ) %>%
  arrange(
    program,
    state
  )

write_tsv(
  group_state3,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_state3_Sham_vs_Tx_3axis_summary.tsv"
  )
)

# ------------------------------------------------------------
# 12. Exact host-weighted change decomposition
#
# Observed biological-sample burden:
#
#   B_i = A_i * F_i * I_i
#
# where:
#   A = HSC abundance among recovered host cells
#   F = state fraction within HSC
#   I = state program intensity
#
# Group mean:
#   mean(B_i)
#
# Product of group means is not generally equal to mean(B_i).
# Therefore:
#
#   mean(B_i) =
#     mean(A)*mean(F)*mean(I) + residual
#
# We decompose product-of-means Sham→Tx by exact Shapley
# contributions for A/F/I, then add:
#
#   residual_component =
#     residual_Tx - residual_Sham
#
# This exactly reconstructs the observed change in group-mean
# sample-level host-weighted burden.
# ------------------------------------------------------------

shapley_product_3 <- function(
  sham,
  tx
) {

  factors <- c(
    "abundance",
    "composition",
    "intensity"
  )

  perms <- list(
    c("abundance", "composition", "intensity"),
    c("abundance", "intensity", "composition"),
    c("composition", "abundance", "intensity"),
    c("composition", "intensity", "abundance"),
    c("intensity", "abundance", "composition"),
    c("intensity", "composition", "abundance")
  )

  contribution <- setNames(
    rep(
      0,
      length(factors)
    ),
    factors
  )

  value_fun <- function(x) {
    prod(
      unlist(
        x[factors]
      )
    )
  }

  for (perm in perms) {

    current <- sham

    for (fac in perm) {

      before <- value_fun(current)

      current[[fac]] <-
        tx[[fac]]

      after <- value_fun(current)

      contribution[[fac]] <-
        contribution[[fac]] +
        (
          after -
            before
        )
    }
  }

  contribution /
    length(perms)
}

decomp_list <- lapply(
  split(
    state3,
    interaction(
      state3$state,
      state3$program,
      drop = TRUE
    )
  ),
  function(d) {

    state_value <-
      unique(
        d$state
      )

    program_value <-
      unique(
        d$program
      )

    d_sham <- d %>%
      filter(
        group ==
          "Sham"
      )

    d_tx <- d %>%
      filter(
        group ==
          "Tx"
      )

    sham_means <- list(
      abundance =
        mean(
          d_sham$HSC_fraction_of_recovered_host
        ),
      composition =
        mean(
          d_sham$state_fraction_within_HSC
        ),
      intensity =
        mean(
          d_sham$mean_intensity_CP10k
        )
    )

    tx_means <- list(
      abundance =
        mean(
          d_tx$HSC_fraction_of_recovered_host
        ),
      composition =
        mean(
          d_tx$state_fraction_within_HSC
        ),
      intensity =
        mean(
          d_tx$mean_intensity_CP10k
        )
    )

    observed_sham <-
      mean(
        d_sham$host_weighted_burden
      )

    observed_tx <-
      mean(
        d_tx$host_weighted_burden
      )

    product_sham <-
      prod(
        unlist(
          sham_means
        )
      )

    product_tx <-
      prod(
        unlist(
          tx_means
        )
      )

    residual_sham <-
      observed_sham -
      product_sham

    residual_tx <-
      observed_tx -
      product_tx

    phi <- shapley_product_3(
      sham_means,
      tx_means
    )

    residual_component <-
      residual_tx -
      residual_sham

    reconstructed_delta <-
      sum(phi) +
      residual_component

    observed_delta <-
      observed_tx -
      observed_sham

    components <- c(
      HSC_abundance =
        unname(
          phi["abundance"]
        ),
      state_composition =
        unname(
          phi["composition"]
        ),
      expression_intensity =
        unname(
          phi["intensity"]
        ),
      sample_interaction_residual =
        residual_component
    )

    tibble(
      state =
        state_value,
      program =
        program_value,

      observed_Sham =
        observed_sham,
      observed_Tx =
        observed_tx,
      observed_delta =
        observed_delta,

      HSC_abundance_component =
        components[
          "HSC_abundance"
        ],

      state_composition_component =
        components[
          "state_composition"
        ],

      expression_intensity_component =
        components[
          "expression_intensity"
        ],

      sample_interaction_residual =
        components[
          "sample_interaction_residual"
        ],

      reconstructed_delta =
        reconstructed_delta,

      decomposition_error =
        abs(
          reconstructed_delta -
            observed_delta
        ),

      dominant_component =
        names(
          components
        )[
          which.max(
            abs(
              components
            )
          )
        ]
    )
  }
)

decomposition <- bind_rows(
  decomp_list
) %>%
  arrange(
    program,
    state
  )

write_tsv(
  decomposition,
  file.path(
    DIR_TABLE,
    "08_PRIMARY_host_weighted_change_decomposition.tsv"
  )
)

if (
  max(
    decomposition$decomposition_error,
    na.rm = TRUE
  ) > 1e-8
) {
  stop(
    "Host-weighted Shapley decomposition failed."
  )
}

# ------------------------------------------------------------
# 13. Whole-HSC total axis output, both within-HSC and host-weighted
# ------------------------------------------------------------

whole_axis_group <- whole_hsc_axis %>%
  group_by(
    program,
    group
  ) %>%
  summarise(
    n_samples = n(),

    mean_HSC_abundance =
      mean(
        HSC_fraction_of_recovered_host
      ),

    mean_whole_HSC_intensity =
      mean(
        whole_HSC_mean_intensity_CP10k
      ),

    mean_host_weighted_burden =
      mean(
        whole_host_weighted_HSC_burden
      ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      n_samples,
      mean_HSC_abundance,
      mean_whole_HSC_intensity,
      mean_host_weighted_burden
    )
  ) %>%
  mutate(
    delta_whole_HSC_intensity =
      mean_whole_HSC_intensity_Tx -
      mean_whole_HSC_intensity_Sham,

    delta_host_weighted_burden =
      mean_host_weighted_burden_Tx -
      mean_host_weighted_burden_Sham,

    ratio_host_weighted_Tx_over_Sham =
      ifelse(
        mean_host_weighted_burden_Sham > 0,
        mean_host_weighted_burden_Tx /
          mean_host_weighted_burden_Sham,
        NA_real_
      ),

    observed_host_weighted_direction =
      case_when(
        delta_host_weighted_burden < 0 ~
          "decrease",
        delta_host_weighted_burden > 0 ~
          "increase",
        TRUE ~
          "no_change"
      )
  )

write_tsv(
  whole_axis_group,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_whole_HSC_3axis_total_output.tsv"
  )
)

# ------------------------------------------------------------
# 14. Gene-level host-weighted burden
# ------------------------------------------------------------

gene_cell <- as.data.frame(
  t(
    as.matrix(
      cp10k
    )
  )
) %>%
  rownames_to_column(
    "cell"
  ) %>%
  left_join(
    hmd,
    by = "cell"
  ) %>%
  pivot_longer(
    cols =
      all_of(
        all_primary_genes
      ),
    names_to =
      "gene",
    values_to =
      "CP10k"
  )

gene_program_map <- bind_rows(
  lapply(
    names(
      PROGRAMS_AVAIL
    ),
    function(p) {
      tibble(
        program = p,
        gene =
          PROGRAMS_AVAIL[[p]]
      )
    }
  )
)

gene_state3 <- gene_cell %>%
  group_by(
    sample,
    group,
    state3,
    gene
  ) %>%
  summarise(
    state_n = n(),
    mean_gene_CP10k =
      mean(
        CP10k,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  left_join(
    hmd %>%
      count(
        sample,
        group,
        name =
          "total_HSC_n"
      ),
    by = c(
      "sample",
      "group"
    )
  ) %>%
  left_join(
    hsc_abundance %>%
      select(
        sample,
        group,
        HSC_fraction_of_recovered_host
      ),
    by = c(
      "sample",
      "group"
    )
  ) %>%
  mutate(
    state_fraction_within_HSC =
      state_n /
      total_HSC_n,

    within_HSC_gene_burden =
      state_fraction_within_HSC *
      mean_gene_CP10k,

    host_weighted_gene_burden =
      HSC_fraction_of_recovered_host *
      state_fraction_within_HSC *
      mean_gene_CP10k
  ) %>%
  left_join(
    gene_program_map,
    by = "gene"
  )

write_tsv(
  gene_state3,
  file.path(
    DIR_TABLE,
    "10_gene_level_state3_host_weighted_burden.tsv.gz"
  )
)

whole_gene_group <- gene_state3 %>%
  group_by(
    sample,
    group,
    program,
    gene
  ) %>%
  summarise(
    whole_HSC_gene_burden =
      sum(
        within_HSC_gene_burden,
        na.rm = TRUE
      ),
    whole_host_weighted_gene_burden =
      sum(
        host_weighted_gene_burden,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  group_by(
    program,
    gene,
    group
  ) %>%
  summarise(
    mean_whole_HSC_gene_burden =
      mean(
        whole_HSC_gene_burden
      ),

    mean_host_weighted_gene_burden =
      mean(
        whole_host_weighted_gene_burden
      ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      mean_whole_HSC_gene_burden,
      mean_host_weighted_gene_burden
    )
  ) %>%
  mutate(
    delta_whole_HSC_gene =
      mean_whole_HSC_gene_burden_Tx -
      mean_whole_HSC_gene_burden_Sham,

    ratio_whole_HSC_Tx_over_Sham =
      ifelse(
        mean_whole_HSC_gene_burden_Sham > 0,
        mean_whole_HSC_gene_burden_Tx /
          mean_whole_HSC_gene_burden_Sham,
        NA_real_
      ),

    delta_host_weighted_gene =
      mean_host_weighted_gene_burden_Tx -
      mean_host_weighted_gene_burden_Sham,

    ratio_host_weighted_Tx_over_Sham =
      ifelse(
        mean_host_weighted_gene_burden_Sham > 0,
        mean_host_weighted_gene_burden_Tx /
          mean_host_weighted_gene_burden_Sham,
        NA_real_
      )
  ) %>%
  arrange(
    program,
    desc(
      abs(
        delta_host_weighted_gene
      )
    )
  )

write_tsv(
  whole_gene_group,
  file.path(
    DIR_TABLE,
    "11_PRIMARY_gene_level_wholeHSC_and_hostWeighted_Sham_vs_Tx.tsv"
  )
)

# ------------------------------------------------------------
# 15. Biological replicate direction
# ------------------------------------------------------------

replicate_direction <- state3 %>%
  group_by(
    state,
    program
  ) %>%
  summarise(
    Sham1 =
      host_weighted_burden[
        sample ==
          "Sham1"
      ][1],

    Sham20 =
      host_weighted_burden[
        sample ==
          "Sham20"
      ][1],

    Tx17 =
      host_weighted_burden[
        sample ==
          "Tx17"
      ][1],

    Tx5 =
      host_weighted_burden[
        sample ==
          "Tx5"
      ][1],

    Sham_mean =
      mean(
        host_weighted_burden[
          group ==
            "Sham"
        ]
      ),

    Tx_mean =
      mean(
        host_weighted_burden[
          group ==
            "Tx"
        ]
      ),

    Tx_both_below_Sham_mean =
      all(
        host_weighted_burden[
          group ==
            "Tx"
        ] <
          mean(
            host_weighted_burden[
              group ==
                "Sham"
            ]
          )
      ),

    Tx_both_above_Sham_mean =
      all(
        host_weighted_burden[
          group ==
            "Tx"
        ] >
          mean(
            host_weighted_burden[
              group ==
                "Sham"
            ]
          )
      ),

    .groups = "drop"
  )

write_tsv(
  replicate_direction,
  file.path(
    DIR_TABLE,
    "12_PRIMARY_hostWeighted_replicate_direction.tsv"
  )
)

# ------------------------------------------------------------
# 16. HSC abundance Sham vs Tx
# ------------------------------------------------------------

hsc_abundance_group <- hsc_abundance %>%
  group_by(group) %>%
  summarise(
    n_samples = n(),
    mean_recovered_host_cells =
      mean(
        recovered_host_cells
      ),
    mean_primary_clean_HSC_cells =
      mean(
        primary_clean_HSC_cells
      ),
    mean_HSC_fraction =
      mean(
        HSC_fraction_of_recovered_host
      ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      group,
    values_from = c(
      n_samples,
      mean_recovered_host_cells,
      mean_primary_clean_HSC_cells,
      mean_HSC_fraction
    )
  ) %>%
  mutate(
    Tx_minus_Sham_HSC_fraction =
      mean_HSC_fraction_Tx -
      mean_HSC_fraction_Sham,

    Tx_over_Sham_HSC_fraction =
      ifelse(
        mean_HSC_fraction_Sham > 0,
        mean_HSC_fraction_Tx /
          mean_HSC_fraction_Sham,
        NA_real_
      )
  )

write_tsv(
  hsc_abundance_group,
  file.path(
    DIR_TABLE,
    "13_HSC_abundance_Sham_vs_Tx.tsv"
  )
)

# ------------------------------------------------------------
# 17. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  hsc_abundance,
  aes(
    x = group,
    y =
      HSC_fraction_of_recovered_host
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
    nudge_y = 0.005,
    check_overlap = TRUE
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "Primary-Clean HSC abundance among recovered Mouse host cells",
    x = NULL,
    y =
      "HSC / recovered Mouse cells"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_HSC_abundance_Sham_vs_Tx.pdf"
  ),
  p1,
  width = 6,
  height = 5
)

p2 <- ggplot(
  state3,
  aes(
    x = group,
    y =
      host_weighted_burden
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
    nudge_y = 0.002,
    check_overlap = TRUE
  ) +
  facet_grid(
    program ~ state,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 10
  ) +
  labs(
    title =
      "HSC state program burden weighted to recovered Mouse host cells",
    x = NULL,
    y =
      "HSC abundance × state fraction × linear CP10k"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_state3_host_weighted_burden.pdf"
  ),
  p2,
  width = 12,
  height = 8
)

decomp_long <- decomposition %>%
  select(
    state,
    program,
    HSC_abundance_component,
    state_composition_component,
    expression_intensity_component,
    sample_interaction_residual
  ) %>%
  pivot_longer(
    cols = c(
      HSC_abundance_component,
      state_composition_component,
      expression_intensity_component,
      sample_interaction_residual
    ),
    names_to =
      "component",
    values_to =
      "delta_contribution"
  )

p3 <- ggplot(
  decomp_long,
  aes(
    x = state,
    y =
      delta_contribution,
    fill = component
  )
) +
  geom_col(
    position = "stack"
  ) +
  geom_hline(
    yintercept = 0
  ) +
  facet_wrap(
    ~ program,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 30,
        hjust = 1
      )
  ) +
  labs(
    title =
      "Exact decomposition of Sham→Tx host-weighted HSC burden",
    x = NULL,
    y =
      "Contribution to Tx - Sham change"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_host_weighted_change_decomposition.pdf"
  ),
  p3,
  width = 10,
  height = 6
)

p4 <- ggplot(
  whole_hsc_axis,
  aes(
    x = group,
    y =
      whole_host_weighted_HSC_burden
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
    nudge_y = 0.002,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ program,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "Whole-HSC 3-axis burden weighted to recovered Mouse host cells",
    x = NULL,
    y =
      "HSC fraction × whole-HSC mean linear CP10k"
  )

ggsave(
  file.path(
    DIR_FIG,
    "04_whole_HSC_3axis_host_weighted.pdf"
  ),
  p4,
  width = 9,
  height = 5
)

# ------------------------------------------------------------
# 18. Save compact results
# ------------------------------------------------------------

saveRDS(
  list(
    HSC_abundance =
      hsc_abundance,
    program_coverage =
      program_coverage,
    state3 =
      state3,
    state2 =
      state2,
    reconstruction =
      recon_check,
    group_state3 =
      group_state3,
    decomposition =
      decomposition,
    whole_axis_group =
      whole_axis_group,
    gene_group =
      whole_gene_group,
    replicate_direction =
      replicate_direction,
    HSC_abundance_group =
      hsc_abundance_group
  ),
  file.path(
    DIR_OBJ,
    "Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4_results.rds"
  )
)

# ------------------------------------------------------------
# 19. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "whole_Mouse_RDS",
    "whole_Mouse_source_cells",
    "whole_Mouse_expected_cells",
    "whole_Mouse_count_matches",
    "whole_sample_column",
    "HSC_RDS",
    "HSC_source_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "HSC_cell_membership_in_whole_fraction",
    "HSC_samples_n",
    "state3_labels",
    "state2_labels",
    "primary_programs",
    "Core_collagen_genes_available_n",
    "ECM_remodeling_genes_available_n",
    "IEG_core_genes_available_n",
    "expression_scale",
    "host_denominator",
    "state3_reconstruction_passed",
    "host_weighted_decomposition_passed",
    "group_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    WHOLE_RDS,
    as.character(
      ncol(whole)
    ),
    as.character(
      EXPECTED_WHOLE_CELLS
    ),
    as.character(
      ncol(whole) ==
        EXPECTED_WHOLE_CELLS
    ),
    whole_sample_col,
    HSC_RDS,
    as.character(
      ncol(hsc)
    ),
    as.character(
      EXPECTED_HSC_PRIMARY
    ),
    as.character(
      ncol(hsc) ==
        EXPECTED_HSC_PRIMARY
    ),
    as.character(
      hsc_cell_membership_fraction
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
          hmd$state3
        )
      ),
      collapse = ";"
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
      names(
        PROGRAMS
      ),
      collapse = ";"
    ),
    as.character(
      length(
        PROGRAMS_AVAIL$Core_collagen
      )
    ),
    as.character(
      length(
        PROGRAMS_AVAIL$ECM_remodeling
      )
    ),
    as.character(
      length(
        PROGRAMS_AVAIL$IEG_core
      )
    ),
    "linear per-cell CP10k",
    "all recovered cells in Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds",
    "TRUE",
    "TRUE",
    "descriptive_only_n2_per_group",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v6.4.4.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "whole_Mouse_source",
    "HSC_source",
    "host_weighted_definition",
    "primary_partition",
    "sensitivity_partition",
    "Core_collagen_panel",
    "ECM_remodeling_panel",
    "IEG_core_panel",
    "IEG_extended_sensitivity_panel",
    "decomposition_method",
    "statistical_unit",
    "interpretation_limit"
  ),
  value = c(
    "v6.4.4",
    WHOLE_RDS,
    HSC_RDS,
    "HSC fraction among recovered Mouse cells × HSC state fraction × mean linear CP10k program intensity",
    HSC_STATE3_COL,
    HSC_STATE2_COL,
    paste(
      CORE_COLLAGEN,
      collapse = ";"
    ),
    paste(
      ECM_REMODELING,
      collapse = ";"
    ),
    paste(
      IEG_CORE,
      collapse = ";"
    ),
    paste(
      IEG_EXTENDED_SENSITIVITY,
      collapse = ";"
    ),
    "exact Shapley decomposition of product-of-means plus exact sample-interaction residual",
    "biological sample",
    "scRNA-recovered host-cell-weighted transcript burden; not absolute tissue ECM production"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "HSC_3Axis_TotalHostWeighted_provenance_v6.4.4.tsv"
  )
)

# ------------------------------------------------------------
# 20. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("HSC ABUNDANCE AMONG RECOVERED MOUSE CELLS\n")
cat("============================================================\n")
print(
  hsc_abundance
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY 3-AXIS COVERAGE\n")
cat("============================================================\n")
print(
  program_coverage
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY STATE3 SHAM vs TX SUMMARY\n")
cat("============================================================\n")
print(
  group_state3
)

cat("\n")
cat("============================================================\n")
cat("EXACT HOST-WEIGHTED CHANGE DECOMPOSITION\n")
cat("============================================================\n")
print(
  decomposition
)

cat("\n")
cat("============================================================\n")
cat("WHOLE-HSC 3-AXIS TOTAL OUTPUT\n")
cat("============================================================\n")
print(
  whole_axis_group
)

cat("\n")
cat("============================================================\n")
cat("TOP GENE-LEVEL HOST-WEIGHTED CHANGES\n")
cat("============================================================\n")
print(
  whole_gene_group %>%
    group_by(program) %>%
    slice_max(
      order_by =
        abs(
          delta_host_weighted_gene
        ),
      n = 8,
      with_ties = FALSE
    ) %>%
    ungroup()
)

cat("\n")
cat("============================================================\n")
cat("REPLICATE DIRECTION\n")
cat("============================================================\n")
print(
  replicate_direction
)

cat("\n")
cat("============================================================\n")
cat("HSC ABUNDANCE SHAM vs TX\n")
cat("============================================================\n")
print(
  hsc_abundance_group
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
