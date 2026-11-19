#!/usr/bin/env Rscript

# ============================================================
# Mouse MASH Cdcp1 whole-cell source decomposition v6.4.11.1
#
# PURPOSE
#   Correct the whole-host Cdcp1 provenance analysis by using
#   FULL raw RNA library sizes as the CPM denominator.
#
# CRITICAL FIX vs v6.4.10.4
#   For each RNA counts layer:
#     1) restrict to the four Sham/Tx biological samples,
#     2) calculate library size from ALL genes in that full layer,
#     3) only then extract Cdcp1 counts.
#
#   Therefore:
#     CPM = sum(Cdcp1 raw UMI) / sum(ALL raw RNA UMI) * 1e6
#
#   The previous v6.4.10.4 denominator based on 11,913 common /
#   ranking genes is NOT reused for corrected quantitative CPM.
#
# OUTPUTS
#   - corrected Mphi vs non-Mphi source decomposition
#   - corrected Sham vs Tx effects by compartment
#   - whole-liver celltype Cdcp1 source decomposition
#   - Mphi-state Cdcp1 source decomposition
#   - positive-cell UMI audit (fixes prior empty audit)
#   - comparison to v6.4.10.4 restricted-denominator values
#   - explicit provenance / integrity
#
# INFERENCE
#   Mouse n=2/group: descriptive only.
#   This analysis addresses source/provenance, not causality.
#
# No RDS is modified. No reclustering.
# ============================================================

options(stringsAsFactors = FALSE, width = 240)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

VERSION <- "v6.4.11.1"

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

MPHI_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

PB_RANK <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1",
  "tables",
  "08_PRIMARY_Cdcp1_transcriptome_rank_summary.tsv"
)

OLD_ROOT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_WholeHost_Provenance_StressAudit_v6.4.10.4"
)

OLD_PROV <- file.path(
  OLD_ROOT,
  "tables",
  "05_PRIMARY_focus_gene_provenance_and_TREM2_rank.tsv"
)

OLD_COMPARTMENT <- file.path(
  OLD_ROOT,
  "tables",
  "07_focus_gene_compartment_Sham_vs_Tx.tsv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Cdcp1_WholeCell_SourceDecomposition_v6.4.11.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG <- file.path(OUT, "figures")
DIR_OBJ <- file.path(OUT, "objects")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")

for (d in c(DIR_TABLE, DIR_FIG, DIR_OBJ, DIR_LOG, DIR_PROV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(
  DIR_LOG,
  "Mouse_MASH_Cdcp1_WholeCell_SourceDecomposition_v6.4.11.1.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("Mouse MASH Cdcp1 whole-cell source decomposition v6.4.11.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

cat(
  "CRITICAL FIX: full raw RNA library size is calculated BEFORE ",
  "feature subsetting; 11,913-gene restricted denominator is not used.\n\n",
  sep = ""
)

cat(
  "Patch note v6.4.11.1: fixes slice_head() top-celltype plotting error; ",
  "all quantitative source-decomposition definitions are unchanged from v6.4.11.\n\n",
  sep = ""
)

stopifnot(
  file.exists(WHOLE_RDS),
  file.exists(MPHI_RDS),
  file.exists(PB_RANK)
)

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")

sample_to_group <- c(
  Sham1 = "Sham",
  Sham20 = "Sham",
  Tx17 = "Tx",
  Tx5 = "Tx"
)

EXPECTED_WHOLE_4SAMPLE <- 80834L
EXPECTED_MPHI_4SAMPLE <- 15645L

PRIMARY_MPHI_STATES <- c(
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Inflammatory-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Repair/Resolution-Mphi"
)

TREM2_STATE <- "Lipid-associated/TREM2-Mphi"

# ------------------------------------------------------------
# Load objects
# ------------------------------------------------------------

cat("Loading whole-liver RDS:\n", WHOLE_RDS, "\n\n")
whole <- readRDS(WHOLE_RDS)

cat("Loading final Mphi RDS:\n", MPHI_RDS, "\n\n")
mphi <- readRDS(MPHI_RDS)

if (!"RNA" %in% Assays(whole)) {
  stop("RNA assay absent from whole object.")
}

if (!"RNA" %in% Assays(mphi)) {
  stop("RNA assay absent from Mphi object.")
}

# ------------------------------------------------------------
# Resolve whole metadata
# ------------------------------------------------------------

wmd0 <- whole[[]] %>%
  rownames_to_column("cell")

WHOLE_SAMPLE_CANDIDATES <- c(
  "sample_display_FIXED2",
  "sample",
  "sample_for_R8plot_FIXED2",
  "sample_for_R8plot",
  "sample_for_annotation",
  "sample_4group",
  "sample_FIXED2"
)

resolve_sample_column <- function(md, candidates, required_samples) {

  for (cc in candidates) {

    if (cc %in% colnames(md)) {

      vals <- unique(as.character(md[[cc]]))

      if (all(required_samples %in% vals)) {
        return(cc)
      }
    }
  }

  matching <- colnames(md)[vapply(
    md,
    function(x) {
      vals <- unique(as.character(x))
      all(required_samples %in% vals)
    },
    logical(1)
  )]

  if (length(matching) == 1) {
    return(matching)
  }

  if (length(matching) == 0) {
    stop("Could not resolve whole-object sample column.")
  }

  stop(
    "Multiple possible whole sample columns: ",
    paste(matching, collapse = ";")
  )
}

whole_sample_col <- resolve_sample_column(
  wmd0,
  WHOLE_SAMPLE_CANDIDATES,
  SAMPLES
)

CELLTYPE_CANDIDATES <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_for_R8plot",
  "celltype_auto_annotation",
  "celltype"
)

whole_celltype_col <- CELLTYPE_CANDIDATES[
  CELLTYPE_CANDIDATES %in% colnames(wmd0)
][1]

if (length(whole_celltype_col) == 0 || is.na(whole_celltype_col)) {
  stop("Whole-liver celltype column could not be resolved.")
}

cat("Resolved whole sample column:", whole_sample_col, "\n")
cat("Resolved whole celltype column:", whole_celltype_col, "\n\n")

wmd <- wmd0 %>%
  transmute(
    cell,
    sample = as.character(.data[[whole_sample_col]]),
    whole_celltype = as.character(.data[[whole_celltype_col]])
  ) %>%
  filter(sample %in% SAMPLES) %>%
  mutate(
    group = unname(sample_to_group[sample])
  )

if (nrow(wmd) != EXPECTED_WHOLE_4SAMPLE) {
  stop(
    "Unexpected whole-cell count across four samples: expected ",
    EXPECTED_WHOLE_4SAMPLE,
    ", observed ",
    nrow(wmd)
  )
}

if (!all(SAMPLES %in% unique(wmd$sample))) {
  stop("Not all four biological samples present in whole metadata.")
}

# ------------------------------------------------------------
# Resolve Mphi metadata and state map
# ------------------------------------------------------------

mmd0 <- mphi[[]] %>%
  rownames_to_column("cell")

MPHI_STATE_CANDIDATES <- c(
  "macrophage_class_Res2_FINAL_v4145_char",
  "macrophage_class_Res2_FINAL_v4145"
)

mphi_state_col <- MPHI_STATE_CANDIDATES[
  MPHI_STATE_CANDIDATES %in% colnames(mmd0)
][1]

if (length(mphi_state_col) == 0 || is.na(mphi_state_col)) {
  stop("Final Mphi state column could not be resolved.")
}

mphi_state_map <- mmd0 %>%
  transmute(
    cell,
    sender_state = as.character(.data[[mphi_state_col]])
  ) %>%
  filter(cell %in% wmd$cell)

if (nrow(mphi_state_map) != EXPECTED_MPHI_4SAMPLE) {
  stop(
    "Unexpected final-Mphi cell count across four samples: expected ",
    EXPECTED_MPHI_4SAMPLE,
    ", observed ",
    nrow(mphi_state_map)
  )
}

if (anyDuplicated(mphi_state_map$cell)) {
  stop("Duplicated cell IDs in Mphi state map.")
}

if (!all(PRIMARY_MPHI_STATES %in% unique(mphi_state_map$sender_state))) {
  stop(
    "One or more primary Mphi states are missing from four-sample Mphi map."
  )
}

# ------------------------------------------------------------
# FULL-LIBRARY counts extraction
#
# IMPORTANT:
#   libsize is calculated using ALL rows of each relevant counts
#   layer before Cdcp1 is extracted.
# ------------------------------------------------------------

rna_assay <- whole[["RNA"]]

all_layers <- Layers(rna_assay)

count_layers <- grep(
  "^counts($|\\.)",
  all_layers,
  value = TRUE
)

if (length(count_layers) == 0) {
  stop(
    "No RNA counts layers found. Layers: ",
    paste(all_layers, collapse = ";")
  )
}

cat(
  "RNA counts layers before filtering:\n",
  paste(count_layers, collapse = "\n"),
  "\n\n"
)

resolve_cdcp1_feature <- function(features) {

  idx_exact <- which(features == "Cdcp1")

  if (length(idx_exact) == 1L) {
    return(features[idx_exact])
  }

  idx_ci <- which(toupper(features) == "CDCP1")

  if (length(idx_ci) == 1L) {
    return(features[idx_ci])
  }

  if (length(idx_ci) > 1L) {
    stop(
      "Ambiguous case-insensitive Cdcp1 feature match: ",
      paste(features[idx_ci], collapse = ";")
    )
  }

  NA_character_
}

extract_layer_full_library <- function(layer_name) {

  mat <- LayerData(
    rna_assay,
    layer = layer_name
  )

  cells_here <- intersect(
    colnames(mat),
    wmd$cell
  )

  if (length(cells_here) == 0) {
    return(NULL)
  }

  # CRITICAL: use all genes here.
  mat_here <- mat[, cells_here, drop = FALSE]

  full_library_UMI <- Matrix::colSums(mat_here)

  if (any(!is.finite(full_library_UMI) | full_library_UMI <= 0)) {
    stop(
      "Invalid full-library size in layer: ",
      layer_name
    )
  }

  cdcp1_feature <- resolve_cdcp1_feature(
    rownames(mat_here)
  )

  if (is.na(cdcp1_feature)) {

    cdcp1_UMI <- rep(
      0,
      length(cells_here)
    )

  } else {

    cdcp1_UMI <- as.numeric(
      mat_here[
        cdcp1_feature,
        cells_here,
        drop = TRUE
      ]
    )
  }

  tibble(
    cell = cells_here,
    counts_layer = layer_name,
    layer_features_n = nrow(mat_here),
    full_library_UMI = as.numeric(full_library_UMI),
    Cdcp1_UMI = cdcp1_UMI,
    Cdcp1_feature = ifelse(
      is.na(cdcp1_feature),
      NA_character_,
      cdcp1_feature
    )
  )
}

layer_results <- lapply(
  count_layers,
  extract_layer_full_library
)

layer_results <- layer_results[
  !vapply(layer_results, is.null, logical(1))
]

if (length(layer_results) == 0) {
  stop("No counts layer contained Sham/Tx cells.")
}

cell_counts <- bind_rows(layer_results)

if (anyDuplicated(cell_counts$cell)) {
  dup <- unique(
    cell_counts$cell[
      duplicated(cell_counts$cell)
    ]
  )

  stop(
    "Cells appeared in multiple selected counts layers. Example: ",
    paste(head(dup, 10), collapse = ";")
  )
}

if (nrow(cell_counts) != EXPECTED_WHOLE_4SAMPLE) {
  stop(
    "Counts extraction did not recover all four-sample cells: expected ",
    EXPECTED_WHOLE_4SAMPLE,
    ", observed ",
    nrow(cell_counts)
  )
}

if (!setequal(cell_counts$cell, wmd$cell)) {
  stop("Counts-extracted cells do not exactly match four-sample metadata cells.")
}

layers_used <- sort(unique(cell_counts$counts_layer))

cat(
  "Counts layers used:\n",
  paste(layers_used, collapse = "\n"),
  "\n\n"
)

cat(
  "Layer feature counts:\n"
)

print(
  cell_counts %>%
    distinct(
      counts_layer,
      layer_features_n,
      Cdcp1_feature
    )
)

# ------------------------------------------------------------
# Assemble cell-level provenance table
# ------------------------------------------------------------

cell_prov <- wmd %>%
  left_join(
    cell_counts,
    by = "cell"
  ) %>%
  left_join(
    mphi_state_map,
    by = "cell"
  ) %>%
  mutate(
    compartment = if_else(
      !is.na(sender_state),
      "Mphi",
      "nonMphi"
    ),
    Cdcp1_detected = Cdcp1_UMI > 0,
    Cdcp1_cell_CPM_full_library =
      Cdcp1_UMI / full_library_UMI * 1e6
  )

if (sum(cell_prov$compartment == "Mphi") != EXPECTED_MPHI_4SAMPLE) {
  stop("Mphi membership count mismatch after whole-cell counts join.")
}

write_tsv(
  cell_prov %>%
    select(
      cell,
      sample,
      group,
      whole_celltype,
      compartment,
      sender_state,
      counts_layer,
      layer_features_n,
      full_library_UMI,
      Cdcp1_UMI,
      Cdcp1_detected,
      Cdcp1_cell_CPM_full_library
    ),
  file.path(
    DIR_TABLE,
    "01_cell_level_Cdcp1_fullLibrary_provenance.tsv.gz"
  )
)

# ------------------------------------------------------------
# Generic pseudobulk summarizer
# ------------------------------------------------------------

summarise_pb <- function(df, group_cols) {

  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      cells = n(),
      Cdcp1_positive_cells =
        sum(Cdcp1_UMI > 0),
      Cdcp1_detection_fraction =
        mean(Cdcp1_UMI > 0),
      Cdcp1_total_UMI =
        sum(Cdcp1_UMI),
      full_library_total_UMI =
        sum(full_library_UMI),
      Cdcp1_CPM_full_library =
        ifelse(
          full_library_total_UMI > 0,
          Cdcp1_total_UMI /
            full_library_total_UMI *
            1e6,
          NA_real_
        ),
      positive_one_UMI_cells =
        sum(Cdcp1_UMI == 1),
      positive_twoPlus_UMI_cells =
        sum(Cdcp1_UMI >= 2),
      median_UMI_among_positive =
        ifelse(
          any(Cdcp1_UMI > 0),
          median(
            Cdcp1_UMI[Cdcp1_UMI > 0]
          ),
          NA_real_
        ),
      max_UMI =
        ifelse(
          any(Cdcp1_UMI > 0),
          max(Cdcp1_UMI),
          NA_real_
        ),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 1) Mphi vs non-Mphi by sample
# ------------------------------------------------------------

compartment_sample <- summarise_pb(
  cell_prov,
  c(
    "sample",
    "group",
    "compartment"
  )
) %>%
  arrange(
    compartment,
    match(sample, SAMPLES)
  )

write_tsv(
  compartment_sample,
  file.path(
    DIR_TABLE,
    "02_PRIMARY_Mphi_vs_nonMphi_Cdcp1_by_sample_fullLibrary.tsv"
  )
)

# Equal-biological-sample treatment summary.
compartment_group <- compartment_sample %>%
  group_by(
    compartment,
    group
  ) %>%
  summarise(
    n_samples = n(),
    mean_CPM_full_library =
      mean(Cdcp1_CPM_full_library),
    mean_detection_fraction =
      mean(Cdcp1_detection_fraction),
    mean_positive_cells =
      mean(Cdcp1_positive_cells),
    mean_Cdcp1_total_UMI =
      mean(Cdcp1_total_UMI),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(
      n_samples,
      mean_CPM_full_library,
      mean_detection_fraction,
      mean_positive_cells,
      mean_Cdcp1_total_UMI
    )
  ) %>%
  mutate(
    CPM_Tx_minus_Sham =
      mean_CPM_full_library_Tx -
      mean_CPM_full_library_Sham,
    CPM_Tx_over_Sham =
      ifelse(
        mean_CPM_full_library_Sham > 0,
        mean_CPM_full_library_Tx /
          mean_CPM_full_library_Sham,
        NA_real_
      ),
    delta_log2_CPM_plus0.5 =
      log2(
        (
          mean_CPM_full_library_Tx +
            0.5
        ) /
          (
            mean_CPM_full_library_Sham +
              0.5
          )
      ),
    detection_Tx_minus_Sham =
      mean_detection_fraction_Tx -
      mean_detection_fraction_Sham,
    observed_CPM_direction =
      case_when(
        CPM_Tx_minus_Sham < 0 ~ "decrease",
        CPM_Tx_minus_Sham > 0 ~ "increase",
        TRUE ~ "no_change"
      )
  )

write_tsv(
  compartment_group,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_Mphi_vs_nonMphi_Cdcp1_Sham_vs_Tx_fullLibrary.tsv"
  )
)

# ------------------------------------------------------------
# 2) Pooled source attribution
# ------------------------------------------------------------

compartment_pooled <- summarise_pb(
  cell_prov,
  c("compartment")
) %>%
  mutate(
    Cdcp1_UMI_share =
      Cdcp1_total_UMI /
      sum(Cdcp1_total_UMI)
  )

mphi_pool <- compartment_pooled %>%
  filter(compartment == "Mphi")

nonmphi_pool <- compartment_pooled %>%
  filter(compartment == "nonMphi")

if (
  nrow(mphi_pool) != 1 ||
    nrow(nonmphi_pool) != 1
) {
  stop("Could not uniquely resolve pooled Mphi/nonMphi compartments.")
}

mphi_vs_nonmphi <- tibble(
  Mphi_CPM_full_library =
    mphi_pool$Cdcp1_CPM_full_library,
  nonMphi_CPM_full_library =
    nonmphi_pool$Cdcp1_CPM_full_library,
  Mphi_detection_fraction =
    mphi_pool$Cdcp1_detection_fraction,
  nonMphi_detection_fraction =
    nonmphi_pool$Cdcp1_detection_fraction,
  Mphi_Cdcp1_UMI =
    mphi_pool$Cdcp1_total_UMI,
  nonMphi_Cdcp1_UMI =
    nonmphi_pool$Cdcp1_total_UMI,
  Mphi_Cdcp1_UMI_share =
    mphi_pool$Cdcp1_UMI_share,
  nonMphi_Cdcp1_UMI_share =
    nonmphi_pool$Cdcp1_UMI_share,
  log2_Mphi_vs_nonMphi_CPM =
    log2(
      (
        mphi_pool$Cdcp1_CPM_full_library +
          0.5
      ) /
        (
          nonmphi_pool$Cdcp1_CPM_full_library +
            0.5
        )
    ),
  empirical_ambient_risk =
    (
      log2(
        (
          mphi_pool$Cdcp1_CPM_full_library +
            0.5
        ) /
          (
            nonmphi_pool$Cdcp1_CPM_full_library +
              0.5
          )
      ) <= -2
    ) &
      (
        mphi_pool$Cdcp1_detection_fraction <
          0.01
      ),
  source_class = case_when(
    mphi_pool$Cdcp1_total_UMI == 0 ~
      "not_detected_in_Mphi",
    nonmphi_pool$Cdcp1_CPM_full_library >
      mphi_pool$Cdcp1_CPM_full_library ~
      "nonMphi_enriched_but_detected_in_Mphi",
    TRUE ~
      "Mphi_not_depleted_vs_nonMphi"
  ),
  interpretation_boundary =
    "Mphi detection supports a plausible endogenous component, but sparse 1-UMI events and non-Mphi enrichment do not exclude ambient RNA; this is not evidence of Mphi-specific causality."
)

write_tsv(
  compartment_pooled,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_Cdcp1_pooled_compartment_source.tsv"
  )
)

write_tsv(
  mphi_vs_nonmphi,
  file.path(
    DIR_TABLE,
    "05_PRIMARY_Cdcp1_Mphi_vs_nonMphi_source_adjudication.tsv"
  )
)

# ------------------------------------------------------------
# 3) Whole-liver celltype source decomposition
# ------------------------------------------------------------

celltype_sample <- summarise_pb(
  cell_prov,
  c(
    "sample",
    "group",
    "whole_celltype"
  )
)

write_tsv(
  celltype_sample,
  file.path(
    DIR_TABLE,
    "06_Cdcp1_whole_celltype_by_sample_fullLibrary.tsv"
  )
)

celltype_pooled <- summarise_pb(
  cell_prov,
  c("whole_celltype")
) %>%
  left_join(
    cell_prov %>%
      group_by(whole_celltype) %>%
      summarise(
        Mphi_membership_fraction =
          mean(compartment == "Mphi"),
        .groups = "drop"
      ),
    by = "whole_celltype"
  ) %>%
  mutate(
    Cdcp1_UMI_share =
      Cdcp1_total_UMI /
      sum(Cdcp1_total_UMI)
  ) %>%
  arrange(
    desc(Cdcp1_UMI_share)
  )

write_tsv(
  celltype_pooled,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_Cdcp1_whole_celltype_source_ranking.tsv"
  )
)

# ------------------------------------------------------------
# 4) Mphi-state source decomposition
# ------------------------------------------------------------

mphi_cells <- cell_prov %>%
  filter(
    compartment == "Mphi"
  )

mphi_state_sample <- summarise_pb(
  mphi_cells,
  c(
    "sample",
    "group",
    "sender_state"
  )
) %>%
  arrange(
    sender_state,
    match(sample, SAMPLES)
  )

write_tsv(
  mphi_state_sample,
  file.path(
    DIR_TABLE,
    "08_Cdcp1_Mphi_state_by_sample_fullLibrary.tsv"
  )
)

mphi_state_pooled <- summarise_pb(
  mphi_cells,
  c("sender_state")
) %>%
  mutate(
    Cdcp1_UMI_share_within_Mphi =
      Cdcp1_total_UMI /
      sum(Cdcp1_total_UMI),
    primary_state =
      sender_state %in%
        PRIMARY_MPHI_STATES
  ) %>%
  arrange(
    desc(Cdcp1_UMI_share_within_Mphi)
  )

write_tsv(
  mphi_state_pooled,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_Cdcp1_Mphi_state_source_ranking.tsv"
  )
)

# ------------------------------------------------------------
# 5) Positive-cell UMI audit
# ------------------------------------------------------------

positive_compartment <- compartment_sample %>%
  select(
    sample,
    group,
    compartment,
    cells,
    Cdcp1_positive_cells,
    Cdcp1_detection_fraction,
    Cdcp1_total_UMI,
    positive_one_UMI_cells,
    positive_twoPlus_UMI_cells,
    median_UMI_among_positive,
    max_UMI
  )

write_tsv(
  positive_compartment,
  file.path(
    DIR_TABLE,
    "10_Cdcp1_positive_cell_UMI_audit_compartment.tsv"
  )
)

positive_state <- mphi_state_sample %>%
  select(
    sample,
    group,
    sender_state,
    cells,
    Cdcp1_positive_cells,
    Cdcp1_detection_fraction,
    Cdcp1_total_UMI,
    positive_one_UMI_cells,
    positive_twoPlus_UMI_cells,
    median_UMI_among_positive,
    max_UMI
  )

write_tsv(
  positive_state,
  file.path(
    DIR_TABLE,
    "11_Cdcp1_positive_cell_UMI_audit_Mphi_state.tsv"
  )
)

# ------------------------------------------------------------
# 6) Independent TREM2 raw-count pseudobulk support
# ------------------------------------------------------------

pb_rank <- read_tsv(
  PB_RANK,
  show_col_types = FALSE
)

trem2_pb <- pb_rank %>%
  filter(
    gene == "Cdcp1",
    sender_state == TREM2_STATE
  )

if (nrow(trem2_pb) != 1) {
  stop("Expected exactly one TREM2 Cdcp1 row in v6.4.9.1.")
}

write_tsv(
  trem2_pb,
  file.path(
    DIR_TABLE,
    "12_TREM2_Cdcp1_independent_rawCount_pseudobulk_support.tsv"
  )
)

# ------------------------------------------------------------
# 7) Compare corrected denominator to v6.4.10.4
# ------------------------------------------------------------

old_new_comparison <- tibble(
  metric = character(),
  old_restricted_denominator = numeric(),
  corrected_full_library_denominator = numeric()
)

if (file.exists(OLD_PROV)) {

  oldp <- read_tsv(
    OLD_PROV,
    show_col_types = FALSE
  ) %>%
    filter(gene == "Cdcp1")

  if (nrow(oldp) == 1) {

    old_new_comparison <- bind_rows(
      old_new_comparison,
      tibble(
        metric = c(
          "Mphi_CPM",
          "nonMphi_CPM",
          "Mphi_detection_fraction",
          "nonMphi_detection_fraction",
          "log2_Mphi_vs_nonMphi_CPM"
        ),
        old_restricted_denominator = c(
          oldp$Mphi_CPM,
          oldp$nonMphi_CPM,
          oldp$Mphi_detection_fraction,
          oldp$nonMphi_detection_fraction,
          oldp$log2_Mphi_vs_nonMphi_CPM
        ),
        corrected_full_library_denominator = c(
          mphi_vs_nonmphi$Mphi_CPM_full_library,
          mphi_vs_nonmphi$nonMphi_CPM_full_library,
          mphi_vs_nonmphi$Mphi_detection_fraction,
          mphi_vs_nonmphi$nonMphi_detection_fraction,
          mphi_vs_nonmphi$log2_Mphi_vs_nonMphi_CPM
        )
      )
    )
  }
}

write_tsv(
  old_new_comparison,
  file.path(
    DIR_TABLE,
    "13_v6.4.10.4_vs_v6.4.11.1_denominator_comparison.tsv"
  )
)

# ------------------------------------------------------------
# 8) Sample-balanced source share
# ------------------------------------------------------------

sample_source_share <- compartment_sample %>%
  group_by(sample, group) %>%
  mutate(
    sample_Cdcp1_UMI_total =
      sum(Cdcp1_total_UMI),
    sample_Cdcp1_UMI_share =
      ifelse(
        sample_Cdcp1_UMI_total > 0,
        Cdcp1_total_UMI /
          sample_Cdcp1_UMI_total,
        NA_real_
      )
  ) %>%
  ungroup()

sample_balanced_share <- sample_source_share %>%
  group_by(compartment) %>%
  summarise(
    n_samples = sum(
      is.finite(sample_Cdcp1_UMI_share)
    ),
    mean_sample_UMI_share =
      mean(
        sample_Cdcp1_UMI_share,
        na.rm = TRUE
      ),
    median_sample_UMI_share =
      median(
        sample_Cdcp1_UMI_share,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_tsv(
  sample_source_share,
  file.path(
    DIR_TABLE,
    "14_Cdcp1_compartment_source_share_by_sample.tsv"
  )
)

write_tsv(
  sample_balanced_share,
  file.path(
    DIR_TABLE,
    "15_PRIMARY_Cdcp1_sample_balanced_source_share.tsv"
  )
)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

p1 <- compartment_sample %>%
  ggplot(
    aes(
      x = group,
      y = Cdcp1_CPM_full_library
    )
  ) +
  geom_point(
    size = 3,
    position = position_jitter(
      width = 0.05,
      height = 0
    )
  ) +
  geom_text(
    aes(label = sample),
    size = 3,
    nudge_y = 0.05,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ compartment,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title = "Cdcp1 pseudobulk CPM using full raw RNA library denominator",
    x = NULL,
    y = "Cdcp1 CPM (full-library denominator)"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_Cdcp1_Mphi_vs_nonMphi_fullLibrary_CPM.pdf"
  ),
  p1,
  width = 8,
  height = 5
)

top_celltypes <- celltype_pooled %>%
  slice_head(n = min(12L, nrow(celltype_pooled)))

p2 <- ggplot(
  top_celltypes,
  aes(
    x = reorder(
      whole_celltype,
      Cdcp1_UMI_share
    ),
    y = Cdcp1_UMI_share
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(
    base_size = 10
  ) +
  labs(
    title = "Whole-liver Cdcp1 raw-UMI source share",
    x = NULL,
    y = "Share of Cdcp1 raw UMI"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_Cdcp1_whole_celltype_source_share.pdf"
  ),
  p2,
  width = 8,
  height = 6
)

p3 <- mphi_state_sample %>%
  filter(
    sender_state %in%
      PRIMARY_MPHI_STATES
  ) %>%
  ggplot(
    aes(
      x = group,
      y = Cdcp1_CPM_full_library
    )
  ) +
  geom_point(
    size = 2.8,
    position = position_jitter(
      width = 0.05,
      height = 0
    )
  ) +
  geom_text(
    aes(label = sample),
    size = 2.5,
    nudge_y = 0.05,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ sender_state,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 9
  ) +
  labs(
    title = "Cdcp1 across primary Mphi states",
    x = NULL,
    y = "Cdcp1 CPM (full-library denominator)"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_Cdcp1_primary_Mphi_states_fullLibrary_CPM.pdf"
  ),
  p3,
  width = 12,
  height = 7
)

# ------------------------------------------------------------
# Provenance / integrity
# ------------------------------------------------------------

layer_audit <- cell_counts %>%
  group_by(
    counts_layer,
    layer_features_n,
    Cdcp1_feature
  ) %>%
  summarise(
    cells = n(),
    total_full_library_UMI =
      sum(full_library_UMI),
    total_Cdcp1_UMI =
      sum(Cdcp1_UMI),
    .groups = "drop"
  )

write_tsv(
  layer_audit,
  file.path(
    DIR_PROV,
    "counts_layer_fullLibrary_audit_v6.4.11.1.tsv"
  )
)

integrity <- tibble(
  item = c(
    "whole_source",
    "Mphi_source",
    "whole_sample_column",
    "whole_celltype_column",
    "Mphi_state_column",
    "whole_cells_four_samples",
    "expected_whole_cells_four_samples",
    "Mphi_cells_four_samples",
    "expected_Mphi_cells_four_samples",
    "Mphi_barcode_membership_in_whole_four_samples",
    "RNA_counts_layers_all",
    "RNA_counts_layers_used",
    "full_library_denominator_definition",
    "feature_subsetting_before_library_size",
    "v6.4.10.4_11913_gene_denominator_reused",
    "Cdcp1_positive_cell_audit_generated",
    "celltype_source_decomposition_generated",
    "Mphi_state_source_decomposition_generated",
    "statistical_unit",
    "group_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    WHOLE_RDS,
    MPHI_RDS,
    whole_sample_col,
    whole_celltype_col,
    mphi_state_col,
    as.character(nrow(wmd)),
    as.character(EXPECTED_WHOLE_4SAMPLE),
    as.character(nrow(mphi_state_map)),
    as.character(EXPECTED_MPHI_4SAMPLE),
    as.character(
      mean(
        mphi_state_map$cell %in%
          cell_counts$cell
      )
    ),
    paste(count_layers, collapse = ";"),
    paste(layers_used, collapse = ";"),
    "sum of ALL raw RNA counts in each original counts layer after cell restriction and BEFORE Cdcp1 feature extraction",
    "FALSE",
    "FALSE",
    "TRUE",
    "TRUE",
    "TRUE",
    "biological sample",
    "descriptive_only_n2_per_group",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v6.4.11.1.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "primary_question",
    "whole_counts_strategy",
    "Cdcp1_CPM_definition",
    "Mphi_definition",
    "whole_celltype_definition",
    "TREM2_state_definition",
    "prior_Mphi_rawCount_support",
    "old_provenance_comparison",
    "ambient_interpretation_limit"
  ),
  value = c(
    VERSION,
    "Is the Mouse Cdcp1 signal Mphi-specific, non-Mphi dominated, or compatible with a sparse endogenous Mphi component with possible ambient contribution?",
    "process each original RNA counts layer separately; calculate full all-gene library size before extracting Cdcp1",
    "sum Cdcp1 raw UMI / sum all-gene raw RNA library UMI * 1e6",
    "barcode membership in final Mphi v4.14.5 RDS",
    whole_celltype_col,
    TREM2_STATE,
    PB_RANK,
    ifelse(
      file.exists(OLD_PROV),
      OLD_PROV,
      "not available"
    ),
    "scRNA counts alone cannot prove whether sparse 1-UMI Mphi detections are endogenous or ambient; provenance classification is therefore bounded and non-causal"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "Cdcp1_source_decomposition_provenance_v6.4.11.1.tsv"
  )
)

saveRDS(
  list(
    layer_audit = layer_audit,
    compartment_sample = compartment_sample,
    compartment_group = compartment_group,
    compartment_pooled = compartment_pooled,
    source_adjudication = mphi_vs_nonmphi,
    celltype_sample = celltype_sample,
    celltype_pooled = celltype_pooled,
    mphi_state_sample = mphi_state_sample,
    mphi_state_pooled = mphi_state_pooled,
    positive_compartment = positive_compartment,
    positive_state = positive_state,
    TREM2_rawCount_support = trem2_pb,
    old_new_comparison = old_new_comparison,
    sample_balanced_share = sample_balanced_share,
    integrity = integrity,
    provenance = provenance
  ),
  file.path(
    DIR_OBJ,
    "Mouse_MASH_Cdcp1_WholeCell_SourceDecomposition_v6.4.11.1_results.rds"
  )
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PRIMARY: CORRECTED Mphi vs non-Mphi SOURCE ADJUDICATION\n")
cat("============================================================\n")
print(mphi_vs_nonmphi)

cat("\n============================================================\n")
cat("PRIMARY: CORRECTED Mphi vs non-Mphi SHAM vs TX\n")
cat("============================================================\n")
print(compartment_group)

cat("\n============================================================\n")
cat("TOP WHOLE-LIVER CELLTYPE SOURCES\n")
cat("============================================================\n")
print(
  celltype_pooled %>%
    select(
      whole_celltype,
      cells,
      Cdcp1_positive_cells,
      Cdcp1_detection_fraction,
      Cdcp1_total_UMI,
      Cdcp1_CPM_full_library,
      Cdcp1_UMI_share,
      Mphi_membership_fraction
    ) %>%
    slice_head(n = 12)
)

cat("\n============================================================\n")
cat("MPHI STATE SOURCE RANKING\n")
cat("============================================================\n")
print(mphi_state_pooled)

cat("\n============================================================\n")
cat("POSITIVE-CELL UMI AUDIT: COMPARTMENT\n")
cat("============================================================\n")
print(positive_compartment)

cat("\n============================================================\n")
cat("TREM2 Cdcp1 INDEPENDENT RAW-COUNT PSEUDOBULK SUPPORT\n")
cat("============================================================\n")
print(trem2_pb)

cat("\n============================================================\n")
cat("v6.4.10.4 vs v6.4.11.1 DENOMINATOR COMPARISON\n")
cat("============================================================\n")
print(old_new_comparison)

cat("\n============================================================\n")
cat("INTEGRITY\n")
cat("============================================================\n")
print(integrity)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
