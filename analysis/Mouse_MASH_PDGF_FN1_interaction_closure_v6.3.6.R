#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6360)

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
# Integrated closure analysis:
# ECM-activated HSC burden + PDGF/FN1 evidence matrix
#
# Version: v6.3.6
#
# PURPOSE
#
#   1) Quantify whether transplantation reduces the burden of cycling
#      ECM-activated HSC by combining:
#
#        ECM-activated HSC fraction within HSC3
#        x
#        cycling-like fraction within ECM-activated HSC
#
#      The product is verified directly as:
#
#        cycling-like ECM-HSC cells / total HSC3 cells
#
#   2) Summarize the two major macrophage -> HSC signaling axes:
#
#      PDGF / mitogenic axis
#        - Repair/Resolution-Mphi per-cell Pdgfb
#        - CellChat PDGF: Repair/Resolution-Mphi -> ECM-activated HSC
#        - ECM-HSC PDGF mitogenic module
#        - cycling-like ECM-HSC fraction
#        - cycling ECM-HSC burden
#
#      FN1 axis
#        - total Mphi5 population-weighted Fn1
#        - CellChat FN1: all Mphi5 -> ECM-activated HSC
#        - HSC FN1 receptor availability
#        - ECM-HSC FN1 remodeling module
#
# INPUTS
#
#   v6.2.0 interaction-ready object
#   v6.2.2 sample-level CellChat pathway table
#   v6.3.3.1 integrated PDGFB sender/receiver table
#   v6.3.5 PDGF/FN1 rewiring outputs
#
# IMPORTANT
#
#   - No CellChat rerun.
#   - No NicheNet rerun.
#   - No reclustering.
#   - Biological n = 2 Sham vs n = 2 Tx.
#   - Evidence grading is descriptive and replicate-aware.
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

save_pdf <- function(
  p,
  file,
  width,
  height
) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}

canonical_condition <- function(sample) {
  ifelse(
    grepl(
      "^Tx",
      sample
    ),
    "Tx",
    "Sham"
  )
}

pairwise_metrics <- function(
  sham1,
  sham20,
  tx17,
  tx5,
  tol = 1e-12
) {

  vals <- c(
    Sham1 = sham1,
    Sham20 = sham20,
    Tx17 = tx17,
    Tx5 = tx5
  )

  if (
    any(
      !is.finite(
        vals
      )
    )
  ) {
    return(
      tibble(
        Sham_mean = NA_real_,
        Tx_mean = NA_real_,
        Tx_minus_Sham = NA_real_,
        pairwise_up_n = NA_integer_,
        pairwise_down_n = NA_integer_,
        pairwise_tie_n = NA_integer_,
        pairwise_direction = "Missing",
        pairwise_direction_consistency = NA_real_,
        both_Tx_below_both_Sham = NA,
        both_Tx_above_both_Sham = NA,
        evidence_grade = "Missing"
      )
    )
  }

  diffs <- c(
    tx17 - sham1,
    tx17 - sham20,
    tx5 - sham1,
    tx5 - sham20
  )

  up_n <- sum(
    diffs >
      tol
  )

  down_n <- sum(
    diffs <
      -tol
  )

  tie_n <- 4L -
    up_n -
    down_n

  direction <- dplyr::case_when(
    up_n >
      down_n ~
      "Tx_up",
    down_n >
      up_n ~
      "Tx_down",
    TRUE ~
      "Mixed_or_tied"
  )

  consistency <- max(
    up_n,
    down_n
  ) /
    4

  both_down <- max(
    tx17,
    tx5
  ) <
    min(
      sham1,
      sham20
    )

  both_up <- min(
    tx17,
    tx5
  ) >
    max(
      sham1,
      sham20
    )

  delta <- mean(
    c(
      tx17,
      tx5
    )
  ) -
    mean(
      c(
        sham1,
        sham20
      )
    )

  grade <- case_when(
    both_down &
      consistency ==
        1 ~
      "Strong_Tx_down",

    both_up &
      consistency ==
        1 ~
      "Strong_Tx_up",

    delta <
      0 &
      consistency >=
        0.75 ~
      "Moderate_Tx_down",

    delta >
      0 &
      consistency >=
        0.75 ~
      "Moderate_Tx_up",

    TRUE ~
      "Mixed_or_weak"
  )

  tibble(
    Sham_mean =
      mean(
        c(
          sham1,
          sham20
        )
      ),

    Tx_mean =
      mean(
        c(
          tx17,
          tx5
        )
      ),

    Tx_minus_Sham =
      delta,

    pairwise_up_n =
      up_n,

    pairwise_down_n =
      down_n,

    pairwise_tie_n =
      tie_n,

    pairwise_direction =
      direction,

    pairwise_direction_consistency =
      consistency,

    both_Tx_below_both_Sham =
      both_down,

    both_Tx_above_both_Sham =
      both_up,

    evidence_grade =
      grade
  )
}

four_sample_summary <- function(
  df,
  value_col,
  metric_name,
  axis,
  layer
) {

  vals <- df %>%
    select(
      sample,
      value =
        all_of(
          value_col
        )
    ) %>%
    distinct(
      sample,
      .keep_all = TRUE
    ) %>%
    complete(
      sample =
        c(
          "Sham1",
          "Sham20",
          "Tx17",
          "Tx5"
        )
    ) %>%
    pivot_wider(
      names_from =
        sample,
      values_from =
        value
    )

  pm <- pairwise_metrics(
    vals$Sham1[[1]],
    vals$Sham20[[1]],
    vals$Tx17[[1]],
    vals$Tx5[[1]]
  )

  bind_cols(
    tibble(
      axis =
        axis,
      evidence_layer =
        layer,
      metric =
        metric_name,
      Sham1 =
        vals$Sham1[[1]],
      Sham20 =
        vals$Sham20[[1]],
      Tx17 =
        vals$Tx17[[1]],
      Tx5 =
        vals$Tx5[[1]]
    ),
    pm
  )
}

clean_pathway_name <- function(x) {
  toupper(
    trimws(
      as.character(
        x
      )
    )
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_interaction_ready_v6.2.0",
  "RDS",
  "Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds"
)

V622 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_CellChat_refined_v6.2.2"
)

CELLCHAT_SAMPLE_FILE <- file.path(
  V622,
  "Tables",
  "05_pathway_supported_LR_strength_by_sample_v6.2.2.csv"
)

V6331 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_population_weighted_v6.3.3.1"
)

PDGFB_SAMPLE_FILE <- file.path(
  V6331,
  "Tables",
  "09_integrated_PDGFB_sender_receiver_by_sample_v6.3.3.1.csv"
)

V635 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGF_mitogenic_FN1_rewiring_v6.3.5"
)

CYCLING_FILE <- file.path(
  V635,
  "Tables",
  "05_ECM_HSC_cycling_like_fraction_by_sample_v6.3.5.csv"
)

PDGF_MODULE_FILE <- file.path(
  V635,
  "Tables",
  "02_ECM_HSC_PDGF_mitogenic_module_by_sample_v6.3.5.csv"
)

FN1_TOTAL_FILE <- file.path(
  V635,
  "Tables",
  "08_total_Mphi5_weighted_Fn1_by_sample_v6.3.5.csv"
)

FN1_RECEPTOR_FILE <- file.path(
  V635,
  "Tables",
  "09_FN1_receptor_availability_HSC3_v6.3.5.csv"
)

FN1_MODULE_FILE <- file.path(
  V635,
  "Tables",
  "10_ECM_HSC_FN1_remodeling_module_by_sample_v6.3.5.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGF_FN1_interaction_closure_v6.3.6"
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

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

for (
  d in c(
    OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT,
    RDS_OUT
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

GROUP_COL <-
  "interaction_celltype_v620"

SAMPLE_COL <-
  "sample_interaction_v620"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

HSC3 <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

MPHI5 <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

FOCAL_HSC <-
  "ECM-activated HSC"

FOCAL_PDGFB_SENDER <-
  "Repair/Resolution-Mphi"

CONDITION_COLORS <- c(
  "Sham" = "#0072B2",
  "Tx" = "#D55E00"
)


# ==============================================================================
# 4. Preflight / load
# ==============================================================================

required_files <- c(
  INPUT_RDS,
  CELLCHAT_SAMPLE_FILE,
  PDGFB_SAMPLE_FILE,
  CYCLING_FILE,
  PDGF_MODULE_FILE,
  FN1_TOTAL_FILE,
  FN1_RECEPTOR_FILE,
  FN1_MODULE_FILE
)

for (
  f in required_files
) {
  if (
    !file.exists(
      f
    )
  ) {
    stop(
      "Required input missing: ",
      f
    )
  }
}

msg(
  "Loading interaction-ready object..."
)

obj <- readRDS(
  INPUT_RDS
)

meta <- obj@meta.data

required_meta <- c(
  GROUP_COL,
  SAMPLE_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    meta
  )
)

if (
  length(
    missing_meta
  )
) {
  stop(
    "Missing metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

cellchat_sample <- read.csv(
  CELLCHAT_SAMPLE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

pdgfb_sample <- read.csv(
  PDGFB_SAMPLE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

cycling <- read.csv(
  CYCLING_FILE,
  check.names = FALSE
) %>%
  as_tibble()

pdgf_module <- read.csv(
  PDGF_MODULE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

fn1_total <- read.csv(
  FN1_TOTAL_FILE,
  check.names = FALSE
) %>%
  as_tibble()

fn1_receptor <- read.csv(
  FN1_RECEPTOR_FILE,
  check.names = FALSE
) %>%
  as_tibble()

fn1_module <- read.csv(
  FN1_MODULE_FILE,
  check.names = FALSE
) %>%
  as_tibble()


# ==============================================================================
# 5. HSC3 composition by sample
# ==============================================================================

msg(
  "Computing HSC3 composition..."
)

hsc_counts <- meta %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  filter(
    .data[[
      SAMPLE_COL
    ]] %in%
      SAMPLES,
    .data[[
      GROUP_COL
    ]] %in%
      HSC3
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    HSC_state =
      as.character(
        .data[[
          GROUP_COL
        ]]
      )
  ) %>%
  count(
    sample,
    HSC_state,
    name =
      "n_cells"
  ) %>%
  complete(
    sample =
      SAMPLES,
    HSC_state =
      HSC3,
    fill = list(
      n_cells = 0
    )
  ) %>%
  group_by(
    sample
  ) %>%
  mutate(
    total_HSC3 =
      sum(
        n_cells
      ),
    fraction_within_HSC3 =
      n_cells /
        total_HSC3
  ) %>%
  ungroup() %>%
  mutate(
    condition =
      canonical_condition(
        sample
      )
  )

write.csv(
  hsc_counts,
  file.path(
    TAB_OUT,
    "01_HSC3_composition_by_sample_v6.3.6.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. Cycling ECM-HSC burden
# ==============================================================================

msg(
  "Computing cycling ECM-HSC burden..."
)

required_cycling_cols <- c(
  "sample",
  "n_ECM_HSC",
  "n_cycling_like",
  "cycling_like_fraction"
)

missing_cycling <- setdiff(
  required_cycling_cols,
  colnames(
    cycling
  )
)

if (
  length(
    missing_cycling
  )
) {
  stop(
    "Cycling table missing columns: ",
    paste(
      missing_cycling,
      collapse = ", "
    )
  )
}

ecm_fraction <- hsc_counts %>%
  filter(
    HSC_state ==
      FOCAL_HSC
  ) %>%
  select(
    sample,
    condition,
    n_ECM_HSC_from_HSC_table =
      n_cells,
    total_HSC3,
    ECM_HSC_fraction =
      fraction_within_HSC3
  )

burden <- ecm_fraction %>%
  left_join(
    cycling %>%
      select(
        all_of(
          required_cycling_cols
        )
      ),
    by =
      "sample"
  ) %>%
  mutate(
    count_identity_ok =
      n_ECM_HSC_from_HSC_table ==
        n_ECM_HSC,

    cycling_ECM_HSC_burden =
      ECM_HSC_fraction *
        cycling_like_fraction,

    direct_cycling_ECM_HSC_burden =
      n_cycling_like /
        total_HSC3,

    burden_identity_error =
      cycling_ECM_HSC_burden -
        direct_cycling_ECM_HSC_burden,

    cycling_ECM_HSC_per_1000_HSC3 =
      1000 *
        direct_cycling_ECM_HSC_burden
  )

if (
  any(
    !burden$count_identity_ok
  )
) {
  stop(
    "ECM-HSC count mismatch between v6.2.0 object and v6.3.5 cycling table."
  )
}

if (
  any(
    abs(
      burden$burden_identity_error
    ) >
      1e-12
  )
) {
  stop(
    "Cycling ECM-HSC burden identity check failed."
  )
}

write.csv(
  burden,
  file.path(
    TAB_OUT,
    "02_ECM_HSC_cycling_burden_by_sample_v6.3.6.csv"
  ),
  row.names = FALSE
)

burden_metrics <- bind_rows(
  four_sample_summary(
    burden,
    "ECM_HSC_fraction",
    "ECM-activated HSC fraction within HSC3",
    "HSC_state",
    "ECM_HSC_fraction"
  ),
  four_sample_summary(
    burden,
    "cycling_like_fraction",
    "Cycling-like fraction within ECM-activated HSC",
    "HSC_state",
    "cycling_fraction"
  ),
  four_sample_summary(
    burden,
    "direct_cycling_ECM_HSC_burden",
    "Cycling ECM-HSC burden within HSC3",
    "HSC_state",
    "cycling_ECM_HSC_burden"
  )
)

write.csv(
  burden_metrics,
  file.path(
    TAB_OUT,
    "03_ECM_HSC_burden_replicate_summary_v6.3.6.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. PDGF evidence layers
# ==============================================================================

msg(
  "Building PDGF evidence layers..."
)

required_pdgfb_cols <- c(
  "sample",
  "RepairResolution_mean_Pdgfb_CP10k"
)

missing_pdgfb <- setdiff(
  required_pdgfb_cols,
  colnames(
    pdgfb_sample
  )
)

if (
  length(
    missing_pdgfb
  )
) {
  stop(
    "v6.3.3.1 PDGFB table missing columns: ",
    paste(
      missing_pdgfb,
      collapse = ", "
    )
  )
}

pdgf_sender_summary <- four_sample_summary(
  pdgfb_sample %>%
    select(
      sample,
      RepairResolution_mean_Pdgfb_CP10k
    ),
  "RepairResolution_mean_Pdgfb_CP10k",
  "Repair/Resolution-Mphi mean Pdgfb CP10k",
  "PDGF",
  "sender_per_cell_ligand"
)

required_cc_cols <- c(
  "sample",
  "source",
  "target",
  "pathway_name",
  "prob_supported"
)

missing_cc <- setdiff(
  required_cc_cols,
  colnames(
    cellchat_sample
  )
)

if (
  length(
    missing_cc
  )
) {
  stop(
    "CellChat sample table missing columns: ",
    paste(
      missing_cc,
      collapse = ", "
    )
  )
}

pdgf_cc_sample <- cellchat_sample %>%
  mutate(
    pathway_clean =
      clean_pathway_name(
        pathway_name
      )
  ) %>%
  filter(
    source ==
      FOCAL_PDGFB_SENDER,
    target ==
      FOCAL_HSC,
    grepl(
      "^PDGF",
      pathway_clean
    )
  ) %>%
  group_by(
    sample
  ) %>%
  summarise(
    value =
      sum(
        prob_supported,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  complete(
    sample =
      SAMPLES,
    fill = list(
      value = 0
    )
  )

if (
  !nrow(
    pdgf_cc_sample
  )
) {
  warning(
    "No PDGF CellChat pathway rows found."
  )
}

pdgf_cellchat_summary <- four_sample_summary(
  pdgf_cc_sample,
  "value",
  "CellChat PDGF supported strength: Repair/Resolution-Mphi -> ECM-HSC",
  "PDGF",
  "CellChat_pathway"
)

pdgf_receiver_summary <- four_sample_summary(
  pdgf_module %>%
    select(
      sample,
      module_score_z
    ),
  "module_score_z",
  "ECM-HSC PDGF mitogenic module",
  "PDGF",
  "receiver_mitogenic_program"
)

pdgf_cycling_summary <- four_sample_summary(
  burden,
  "cycling_like_fraction",
  "Cycling-like fraction within ECM-HSC",
  "PDGF",
  "receiver_cycling_fraction"
)

pdgf_burden_summary <- four_sample_summary(
  burden,
  "direct_cycling_ECM_HSC_burden",
  "Cycling ECM-HSC burden within HSC3",
  "PDGF",
  "receiver_population_burden"
)


# ==============================================================================
# 8. FN1 evidence layers
# ==============================================================================

msg(
  "Building FN1 evidence layers..."
)

required_fn1_total <- c(
  "sample",
  "total_Mphi5_weighted_Fn1"
)

missing_fn1_total <- setdiff(
  required_fn1_total,
  colnames(
    fn1_total
  )
)

if (
  length(
    missing_fn1_total
  )
) {
  stop(
    "FN1 total table missing columns: ",
    paste(
      missing_fn1_total,
      collapse = ", "
    )
  )
}

fn1_sender_summary <- four_sample_summary(
  fn1_total,
  "total_Mphi5_weighted_Fn1",
  "Total Mphi5 population-weighted Fn1",
  "FN1",
  "sender_population_ligand"
)

fn1_cc_sample <- cellchat_sample %>%
  mutate(
    pathway_clean =
      clean_pathway_name(
        pathway_name
      )
  ) %>%
  filter(
    source %in%
      MPHI5,
    target ==
      FOCAL_HSC,
    grepl(
      "^FN1$",
      pathway_clean
    )
  ) %>%
  group_by(
    sample
  ) %>%
  summarise(
    value =
      sum(
        prob_supported,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  complete(
    sample =
      SAMPLES,
    fill = list(
      value = 0
    )
  )

if (
  !nrow(
    fn1_cc_sample
  )
) {
  warning(
    "No FN1 CellChat pathway rows found."
  )
}

fn1_cellchat_summary <- four_sample_summary(
  fn1_cc_sample,
  "value",
  "CellChat FN1 supported strength: all Mphi5 -> ECM-HSC",
  "FN1",
  "CellChat_pathway"
)

fn1_receiver_summary <- four_sample_summary(
  fn1_module %>%
    select(
      sample,
      module_score_z
    ),
  "module_score_z",
  "ECM-HSC FN1 remodeling module",
  "FN1",
  "receiver_transcriptional_program"
)


# ==============================================================================
# 9. FN1 receptor availability summary
# ==============================================================================

required_receptor_cols <- c(
  "sample",
  "HSC_state",
  "receptor_gene",
  "positive_fraction",
  "mean_logexpr"
)

missing_receptor <- setdiff(
  required_receptor_cols,
  colnames(
    fn1_receptor
  )
)

if (
  length(
    missing_receptor
  )
) {
  stop(
    "FN1 receptor table missing columns: ",
    paste(
      missing_receptor,
      collapse = ", "
    )
  )
}

fn1_receptor_ecm <- fn1_receptor %>%
  filter(
    HSC_state ==
      FOCAL_HSC
  )

fn1_receptor_summary_table <- fn1_receptor_ecm %>%
  group_by(
    receptor_gene
  ) %>%
  summarise(
    Sham1_positive_fraction =
      positive_fraction[
        sample ==
          "Sham1"
      ][1],

    Sham20_positive_fraction =
      positive_fraction[
        sample ==
          "Sham20"
      ][1],

    Tx17_positive_fraction =
      positive_fraction[
        sample ==
          "Tx17"
      ][1],

    Tx5_positive_fraction =
      positive_fraction[
        sample ==
          "Tx5"
      ][1],

    mean_positive_fraction =
      mean(
        positive_fraction,
        na.rm = TRUE
      ),

    mean_logexpr =
      mean(
        mean_logexpr,
        na.rm = TRUE
      ),

    available_all_samples =
      all(
        positive_fraction >
          0
      ),

    .groups = "drop"
  )

write.csv(
  fn1_receptor_summary_table,
  file.path(
    TAB_OUT,
    "04_FN1_receptor_availability_ECM_HSC_summary_v6.3.6.csv"
  ),
  row.names = FALSE
)

fn1_receptor_evidence <- tibble(
  axis =
    "FN1",
  evidence_layer =
    "receiver_receptor_availability",
  metric =
    "FN1 receptors present in ECM-HSC",
  Sham1 =
    NA_real_,
  Sham20 =
    NA_real_,
  Tx17 =
    NA_real_,
  Tx5 =
    NA_real_,
  Sham_mean =
    NA_real_,
  Tx_mean =
    NA_real_,
  Tx_minus_Sham =
    NA_real_,
  pairwise_up_n =
    NA_integer_,
  pairwise_down_n =
    NA_integer_,
  pairwise_tie_n =
    NA_integer_,
  pairwise_direction =
    "Availability",
  pairwise_direction_consistency =
    NA_real_,
  both_Tx_below_both_Sham =
    NA,
  both_Tx_above_both_Sham =
    NA,
  evidence_grade =
    ifelse(
      any(
        fn1_receptor_summary_table$available_all_samples
      ),
      "Receptor_available",
      "Limited_receptor_availability"
    )
)


# ==============================================================================
# 10. Evidence matrix
# ==============================================================================

evidence_matrix <- bind_rows(
  pdgf_sender_summary,
  pdgf_cellchat_summary,
  pdgf_receiver_summary,
  pdgf_cycling_summary,
  pdgf_burden_summary,
  fn1_sender_summary,
  fn1_cellchat_summary,
  fn1_receptor_evidence,
  fn1_receiver_summary
) %>%
  mutate(
    expected_model_direction = case_when(
      axis ==
        "PDGF" &
        evidence_layer !=
          "receiver_receptor_availability" ~
        "Tx_down",

      axis ==
        "FN1" &
        evidence_layer %in%
          c(
            "sender_population_ligand",
            "CellChat_pathway"
          ) ~
        "Tx_up",

      axis ==
        "FN1" &
        evidence_layer ==
          "receiver_transcriptional_program" ~
        "Tx_up_if_transcriptionally_coupled",

      axis ==
        "FN1" &
        evidence_layer ==
          "receiver_receptor_availability" ~
        "Present",

      TRUE ~
        NA_character_
    ),

    model_concordance = case_when(
      expected_model_direction ==
        "Tx_down" &
        evidence_grade %in%
          c(
            "Strong_Tx_down",
            "Moderate_Tx_down"
          ) ~
        "Concordant",

      expected_model_direction ==
        "Tx_up" &
        evidence_grade %in%
          c(
            "Strong_Tx_up",
            "Moderate_Tx_up"
          ) ~
        "Concordant",

      expected_model_direction ==
        "Tx_up_if_transcriptionally_coupled" &
        evidence_grade %in%
          c(
            "Strong_Tx_up",
            "Moderate_Tx_up"
          ) ~
        "Concordant",

      expected_model_direction ==
        "Tx_up_if_transcriptionally_coupled" &
        evidence_grade ==
          "Mixed_or_weak" ~
        "Not_supported",

      expected_model_direction ==
        "Present" &
        evidence_grade ==
          "Receptor_available" ~
        "Concordant",

      evidence_grade ==
        "Mixed_or_weak" ~
        "Mixed",

      TRUE ~
        "Discordant_or_weak"
    )
  )

write.csv(
  evidence_matrix,
  file.path(
    TAB_OUT,
    "05_PDGF_FN1_evidence_matrix_v6.3.6.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Axis-level conclusion table
# ==============================================================================

pdgf_concordant_n <- evidence_matrix %>%
  filter(
    axis ==
      "PDGF",
    model_concordance ==
      "Concordant"
  ) %>%
  nrow()

pdgf_total_n <- evidence_matrix %>%
  filter(
    axis ==
      "PDGF"
  ) %>%
  nrow()

fn1_sender_comm_support <- evidence_matrix %>%
  filter(
    axis ==
      "FN1",
    evidence_layer %in%
      c(
        "sender_population_ligand",
        "CellChat_pathway"
      ),
    model_concordance ==
      "Concordant"
  ) %>%
  nrow()

fn1_transcription_support <- evidence_matrix %>%
  filter(
    axis ==
      "FN1",
    evidence_layer ==
      "receiver_transcriptional_program",
    model_concordance ==
      "Concordant"
  ) %>%
  nrow()

axis_conclusions <- tibble(
  axis = c(
    "PDGF_mitogenic",
    "FN1_communication"
  ),

  conclusion = c(
    ifelse(
      pdgf_concordant_n >=
        4,
      "Coherent_Tx_down_mitogenic_axis",
      "Partial_or_mixed_PDGF_axis"
    ),

    ifelse(
      fn1_sender_comm_support ==
        2 &
        fn1_transcription_support ==
          0,
      "Sender_and_communication_Tx_up_but_receiver_transcription_mixed",
      ifelse(
        fn1_sender_comm_support ==
          2 &
          fn1_transcription_support ==
            1,
        "Coherent_FN1_Tx_up_axis",
        "Partial_or_mixed_FN1_axis"
      )
    )
  ),

  interpretation = c(
    "PDGF-related evidence is evaluated as a selective mitogenic/proliferative arm rather than global HSC signaling suppression.",
    "FN1 is interpreted primarily as increased macrophage-to-HSC communication/adhesion input unless an independent downstream transcriptional program is also reproducibly increased."
  )
)

write.csv(
  axis_conclusions,
  file.path(
    TAB_OUT,
    "06_PDGF_FN1_axis_conclusions_v6.3.6.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Figure 1: HSC3 composition
# ==============================================================================

p1 <- ggplot(
  hsc_counts %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        ),
      HSC_state =
        factor(
          HSC_state,
          levels =
            HSC3
        )
    ),
  aes(
    x =
      sample,
    y =
      100 *
        fraction_within_HSC3,
    fill =
      HSC_state
  )
) +
  geom_col() +
  labs(
    title =
      "HSC3 composition by biological sample",
    x = NULL,
    y =
      "% of HSC3",
    fill =
      "HSC state"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_HSC3_composition_by_sample_v6.3.6.pdf"
  ),
  9,
  6
)


# ==============================================================================
# 13. Figure 2: ECM-HSC fraction, cycling fraction, cycling burden
# ==============================================================================

p2a <- ggplot(
  burden %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      100 *
        ECM_HSC_fraction,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "ECM-activated HSC fraction",
    x = NULL,
    y =
      "% of HSC3",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p2b <- ggplot(
  burden %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      100 *
        cycling_like_fraction,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Cycling-like fraction within ECM-HSC",
    x = NULL,
    y =
      "%",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p2c <- ggplot(
  burden %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      cycling_ECM_HSC_per_1000_HSC3,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Cycling ECM-HSC burden",
    subtitle =
      "Directly calculated as cycling ECM-HSC / total HSC3",
    x = NULL,
    y =
      "Cells per 1,000 HSC3",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p2 <- p2a +
  p2b +
  p2c +
  plot_annotation(
    title =
      "ECM-activated HSC abundance x mitogenic activity | v6.3.6"
  )

save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_ECM_HSC_fraction_cycling_burden_v6.3.6.pdf"
  ),
  15,
  5.5
)


# ==============================================================================
# 14. Figure 3: PDGF evidence
# ==============================================================================

pdgf_plot <- evidence_matrix %>%
  filter(
    axis ==
      "PDGF"
  ) %>%
  select(
    evidence_layer,
    metric,
    Sham1,
    Sham20,
    Tx17,
    Tx5
  ) %>%
  pivot_longer(
    cols = all_of(
      SAMPLES
    ),
    names_to =
      "sample",
    values_to =
      "value"
  ) %>%
  filter(
    is.finite(
      value
    )
  ) %>%
  mutate(
    condition =
      canonical_condition(
        sample
      ),
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p3 <- ggplot(
  pdgf_plot,
  aes(
    x =
      sample,
    y =
      value,
    fill =
      condition
  )
) +
  geom_col() +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "PDGF / mitogenic evidence across biological samples",
    subtitle =
      "Sender ligand, CellChat communication, receiver mitogenic activity and cycling burden",
    x = NULL,
    y = NULL,
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_PDGF_mitogenic_evidence_panel_v6.3.6.pdf"
  ),
  12,
  10
)


# ==============================================================================
# 15. Figure 4: FN1 evidence
# ==============================================================================

fn1_plot <- evidence_matrix %>%
  filter(
    axis ==
      "FN1",
    evidence_layer !=
      "receiver_receptor_availability"
  ) %>%
  select(
    evidence_layer,
    metric,
    Sham1,
    Sham20,
    Tx17,
    Tx5
  ) %>%
  pivot_longer(
    cols = all_of(
      SAMPLES
    ),
    names_to =
      "sample",
    values_to =
      "value"
  ) %>%
  filter(
    is.finite(
      value
    )
  ) %>%
  mutate(
    condition =
      canonical_condition(
        sample
      ),
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p4 <- ggplot(
  fn1_plot,
  aes(
    x =
      sample,
    y =
      value,
    fill =
      condition
  )
) +
  geom_col() +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "FN1 communication evidence across biological samples",
    subtitle =
      "Macrophage weighted Fn1, CellChat FN1 and ECM-HSC remodeling transcription",
    x = NULL,
    y = NULL,
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p4,
  file.path(
    FIG_OUT,
    "04_FN1_communication_evidence_panel_v6.3.6.pdf"
  ),
  12,
  7
)


# ==============================================================================
# 16. Figure 5: evidence matrix
# ==============================================================================

grade_score <- c(
  "Strong_Tx_down" = -2,
  "Moderate_Tx_down" = -1,
  "Mixed_or_weak" = 0,
  "Moderate_Tx_up" = 1,
  "Strong_Tx_up" = 2,
  "Receptor_available" = 0.5,
  "Limited_receptor_availability" = 0
)

evidence_heat <- evidence_matrix %>%
  mutate(
    score =
      unname(
        grade_score[
          evidence_grade
        ]
      ),
    label =
      paste0(
        axis,
        " | ",
        evidence_layer
      ),
    label =
      factor(
        label,
        levels =
          rev(
            unique(
              label
            )
          )
      )
  )

p5 <- ggplot(
  evidence_heat,
  aes(
    x =
      axis,
    y =
      label,
    fill =
      score
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label =
        evidence_grade
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(
      -2,
      2
    )
  ) +
  labs(
    title =
      "PDGF vs FN1 evidence matrix",
    subtitle =
      "Blue = Tx-down evidence; red = Tx-up evidence; receptor availability is qualitative",
    x = NULL,
    y = NULL,
    fill =
      "Direction score"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p5,
  file.path(
    FIG_OUT,
    "05_PDGF_FN1_evidence_matrix_v6.3.6.pdf"
  ),
  10,
  7
)


# ==============================================================================
# 17. Figure 6: integrated closure figure
# ==============================================================================

p6 <- (
  p2 /
    (
      p3 +
        p4
    )
) +
  plot_annotation(
    title =
      "Macrophage -> HSC interaction closure | v6.3.6"
  )

save_pdf(
  p6,
  file.path(
    FIG_OUT,
    "06_Mphi_HSC_interaction_closure_v6.3.6.pdf"
  ),
  15,
  15
)


# ==============================================================================
# 18. Save compact RDS
# ==============================================================================

results <- list(
  HSC3_composition =
    hsc_counts,
  cycling_ECM_HSC_burden =
    burden,
  burden_replicate_summary =
    burden_metrics,
  PDGF_CellChat_sample =
    pdgf_cc_sample,
  FN1_CellChat_sample =
    fn1_cc_sample,
  FN1_receptor_summary =
    fn1_receptor_summary_table,
  evidence_matrix =
    evidence_matrix,
  axis_conclusions =
    axis_conclusions
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_PDGF_FN1_interaction_closure_results_v6.3.6.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 19. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "interaction_ready_input",
    "CellChat_sample_input",
    "PDGFB_sender_input",
    "v6.3.5_input_directory",
    "cycling_burden_definition",
    "PDGF_axis_definition",
    "FN1_axis_definition",
    "biological_replicates",
    "CellChat_recomputed",
    "NicheNet_recomputed",
    "formal_inference_note"
  ),

  value = c(
    "v6.3.6",
    INPUT_RDS,
    CELLCHAT_SAMPLE_FILE,
    PDGFB_SAMPLE_FILE,
    V635,
    "cycling-like ECM-HSC cells / total HSC3 cells",
    "Repair/Resolution per-cell Pdgfb + CellChat PDGF + ECM-HSC mitogenic module + cycling fraction + cycling burden",
    "Mphi5 weighted Fn1 + CellChat FN1 + HSC receptor availability + ECM-HSC remodeling module",
    "Sham n=2; Tx n=2",
    "FALSE",
    "FALSE",
    "Exploratory replicate-aware synthesis; no new formal group-level inference"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.6.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.6.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
