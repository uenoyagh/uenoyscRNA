#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6620)

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Mphi5 <-> Hep5 CellChat replicate-aware refinement
#
# Version: v6.6.2.1
#
# FIX FROM v6.6.2
# -----------------
# v6.6.2 used scalar logical AND (&&) in a dplyr::mutate() expression for
# the vector-valued CellChat p-value column.
# v6.6.2.1 corrects this to vectorized AND (&):
#   supported = !is.na(cellchat_pval) & cellchat_pval < 0.05
#
# The remaining && operators in pairwise_direction_metrics() are intentional:
# they combine scalar summary values, not vectors.
#
#
# INPUT
# -----
# Sample-wise CellChat objects from v6.6.1:
#   CellChat_Sham1_v6.6.1.rds
#   CellChat_Sham20_v6.6.1.rds
#   CellChat_Tx17_v6.6.1.rds
#   CellChat_Tx5_v6.6.1.rds
#
# PURPOSE
# -------
# Re-extract all cross-lineage ligand-receptor interactions from each
# biological sample and separate:
#
#   raw probability
#   supported probability (CellChat pval < 0.05)
#
# Then compare Sham vs Tx at biological-sample level.
#
# PRIMARY FOCUS
# -------------
#   Mphi -> Hepatocyte
#   Hepatocyte -> Mphi
#
# with special emphasis on:
#   PDGF
#   FN1
#   SPP1
#   TNF
#   SEMA4
#   MIF
#   COLLAGEN
#   LAMININ
#
# IMPORTANT
# ---------
# - CellChat inference is NOT rerun.
# - n = 2 Sham vs n = 2 Tx biological samples.
# - No pseudo-replication-based significance claim is made.
# - CellChat pval is used only as within-sample communication support.
# - Final ranking emphasizes:
#       support in biological replicates
#       Tx-Sham effect size
#       four pairwise Tx-vs-Sham direction consistency
# ==============================================================================


# ==============================================================================
# 1. Helpers
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

  print(
    p
  )

  grDevices::dev.off()
}


make_interaction_key <- function(df) {

  if (
    "interaction_name_2" %in%
      colnames(df)
  ) {
    return(
      as.character(
        df$interaction_name_2
      )
    )
  }

  if (
    "interaction_name" %in%
      colnames(df)
  ) {
    return(
      as.character(
        df$interaction_name
      )
    )
  }

  if (
    all(
      c(
        "ligand",
        "receptor"
      ) %in%
        colnames(df)
    )
  ) {
    return(
      paste(
        df$ligand,
        df$receptor,
        sep = " -> "
      )
    )
  }

  paste0(
    "interaction_",
    seq_len(
      nrow(df)
    )
  )
}


priority_family_flag <- function(pathway_name) {

  grepl(
    paste(
      c(
        "^PDGF$",
        "^FN1$",
        "^SPP1$",
        "^TNF$",
        "^SEMA4$",
        "^MIF$",
        "^COLLAGEN$",
        "^LAMININ$",
        "^IL1$",
        "^IL6$",
        "^OSM$",
        "^TGF",
        "TGFB",
        "^GAS$",
        "^IGF$",
        "^HGF$",
        "^VEGF$",
        "^THBS$",
        "^CCL$",
        "^CXCL$",
        "^CSF$",
        "^APP$",
        "^PLAU$",
        "^GRN$",
        "^VCAM$",
        "^ICAM$"
      ),
      collapse = "|"
    ),
    as.character(pathway_name),
    ignore.case = TRUE
  )
}


core_family_flag <- function(pathway_name) {

  grepl(
    paste(
      c(
        "^PDGF$",
        "^FN1$",
        "^SPP1$",
        "^TNF$",
        "^SEMA4$",
        "^MIF$",
        "^COLLAGEN$",
        "^LAMININ$"
      ),
      collapse = "|"
    ),
    as.character(pathway_name),
    ignore.case = TRUE
  )
}


pairwise_direction_metrics <- function(
  sham1,
  sham20,
  tx17,
  tx5
) {

  diffs <- c(
    tx17 - sham1,
    tx17 - sham20,
    tx5 - sham1,
    tx5 - sham20
  )

  up_n <- sum(
    diffs > 0
  )

  down_n <- sum(
    diffs < 0
  )

  tie_n <- sum(
    diffs == 0
  )

  up_fraction <- up_n / 4
  down_fraction <- down_n / 4

  direction <- case_when(
    up_n > down_n ~
      "Tx_up",

    down_n > up_n ~
      "Tx_down",

    TRUE ~
      "Mixed_or_tied"
  )

  consistency <- max(
    up_fraction,
    down_fraction
  )

  both_tx_above <- (
    min(
      tx17,
      tx5
    ) >
      max(
        sham1,
        sham20
      )
  )

  both_tx_below <- (
    max(
      tx17,
      tx5
    ) <
      min(
        sham1,
        sham20
      )
  )

  evidence_grade <- case_when(
    both_tx_above ~
      "Strong_Tx_up",

    both_tx_below ~
      "Strong_Tx_down",

    direction ==
      "Tx_up" &&
      consistency >=
        0.75 ~
      "Moderate_Tx_up",

    direction ==
      "Tx_down" &&
      consistency >=
        0.75 ~
      "Moderate_Tx_down",

    TRUE ~
      "Mixed_or_weak"
  )

  tibble(
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

    both_Tx_above_both_Sham =
      both_tx_above,

    both_Tx_below_both_Sham =
      both_tx_below,

    evidence_grade =
      evidence_grade
  )
}


four_sample_summary <- function(
  df,
  key_cols,
  value_col
) {

  samples <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  keys <- df %>%
    distinct(
      across(
        all_of(
          key_cols
        )
      )
    )

  grid <- tidyr::crossing(
    keys,
    sample =
      samples
  )

  vals <- df %>%
    group_by(
      across(
        all_of(
          key_cols
        )
      ),
      sample
    ) %>%
    summarise(
      value =
        sum(
          .data[[
            value_col
          ]],
          na.rm =
            TRUE
        ),
      .groups =
        "drop"
    )

  x <- grid %>%
    left_join(
      vals,
      by = c(
        key_cols,
        "sample"
      )
    ) %>%
    mutate(
      value =
        replace_na(
          value,
          0
        )
    ) %>%
    pivot_wider(
      names_from =
        sample,
      values_from =
        value,
      values_fill =
        0
    )

  metrics <- lapply(
    seq_len(
      nrow(x)
    ),
    function(i) {

      pairwise_direction_metrics(
        x$Sham1[[i]],
        x$Sham20[[i]],
        x$Tx17[[i]],
        x$Tx5[[i]]
      )
    }
  ) %>%
    bind_rows()

  bind_cols(
    x,
    metrics
  ) %>%
    mutate(
      Sham_mean =
        (
          Sham1 +
            Sham20
        ) /
          2,

      Tx_mean =
        (
          Tx17 +
            Tx5
        ) /
          2,

      Tx_minus_Sham =
        Tx_mean -
          Sham_mean,

      detected_Sham_replicates =
        as.integer(
          Sham1 > 0
        ) +
          as.integer(
            Sham20 > 0
          ),

      detected_Tx_replicates =
        as.integer(
          Tx17 > 0
        ) +
          as.integer(
            Tx5 > 0
          ),

      detected_samples =
        detected_Sham_replicates +
          detected_Tx_replicates
    )
}


extract_all_crosslineage <- function(
  cellchat,
  sample_name,
  condition_name,
  mphi_states,
  hep_states
) {

  extract_one <- function(
    sources,
    targets,
    direction_name
  ) {

    df <- subsetCommunication(
      cellchat,
      sources.use =
        sources,
      targets.use =
        targets,
      thresh =
        1
    )

    if (
      is.null(df) ||
      !nrow(df)
    ) {
      return(
        tibble()
      )
    }

    df <- as_tibble(
      df
    )

    if (
      !"prob" %in%
        colnames(df)
    ) {
      stop(
        "CellChat subsetCommunication output lacks 'prob'."
      )
    }

    if (
      !"pval" %in%
        colnames(df)
    ) {
      stop(
        "CellChat subsetCommunication output lacks 'pval'."
      )
    }

    if (
      !"pathway_name" %in%
        colnames(df)
    ) {
      stop(
        "CellChat subsetCommunication output lacks 'pathway_name'."
      )
    }

    if (
      !"ligand" %in%
        colnames(df)
    ) {
      df$ligand <-
        NA_character_
    }

    if (
      !"receptor" %in%
        colnames(df)
    ) {
      df$receptor <-
        NA_character_
    }

    df$interaction_key <-
      make_interaction_key(
        df
      )

    df %>%
      mutate(
        sample =
          sample_name,

        condition =
          condition_name,

        direction =
          direction_name,

        raw_prob =
          as.numeric(
            prob
          ),

        cellchat_pval =
          as.numeric(
            pval
          ),

        supported =
          !is.na(
            cellchat_pval
          ) &
            cellchat_pval <
              0.05,

        supported_prob =
          ifelse(
            !is.na(
              cellchat_pval
            ) &
              cellchat_pval <
                0.05,
            raw_prob,
            0
          )
      )
  }

  bind_rows(
    extract_one(
      mphi_states,
      hep_states,
      "Mphi_to_Hep"
    ),

    extract_one(
      hep_states,
      mphi_states,
      "Hep_to_Mphi"
    )
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"


IN_DIR <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_v6.6.1",
  "RDS"
)


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_refine_v6.6.2.1"
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
# 3. Biological states
# ==============================================================================

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)


MPHI_STATES <- c(
  "Mphi_Anti-inflammatory",
  "Mphi_Inflammatory",
  "Mphi_ECM-associated-inflammatory",
  "Mphi_Repair-Resolution",
  "Mphi_Lipid-associated-TREM2"
)


HEP_STATES <- c(
  "Hep_Periportal",
  "Hep_Pericentral",
  "Hep_Injury-inflammatory",
  "Hep_Intermediate",
  "Hep_Cycling"
)


# ==============================================================================
# 4. Load sample-wise CellChat objects and extract raw/supported LR
# ==============================================================================

lr_rows <- list()


for (
  sample_name in SAMPLES
) {

  file <- file.path(
    IN_DIR,
    paste0(
      "CellChat_",
      sample_name,
      "_v6.6.1.rds"
    )
  )

  if (
    !file.exists(file)
  ) {
    stop(
      "Missing CellChat sample RDS: ",
      file
    )
  }

  msg(
    "Loading ",
    sample_name
  )

  cc <- readRDS(
    file
  )

  condition_name <- ifelse(
    grepl(
      "^Sham",
      sample_name
    ),
    "Sham",
    "Tx"
  )

  tmp <- extract_all_crosslineage(
    cellchat =
      cc,
    sample_name =
      sample_name,
    condition_name =
      condition_name,
    mphi_states =
      MPHI_STATES,
    hep_states =
      HEP_STATES
  )

  lr_rows[[
    sample_name
  ]] <- tmp

  rm(cc)
  gc()
}


lr_all <- bind_rows(
  lr_rows
)


if (
  !nrow(
    lr_all
  )
) {
  stop(
    "No cross-lineage LR rows extracted."
  )
}


write.csv(
  lr_all,
  file.path(
    TAB_OUT,
    "01_crosslineage_LR_raw_supported_all_samples_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 5. Basic support audit
# ==============================================================================

support_audit <- lr_all %>%
  group_by(
    sample,
    condition,
    direction
  ) %>%
  summarise(
    n_raw_rows =
      n(),

    n_supported_rows =
      sum(
        supported,
        na.rm =
          TRUE
      ),

    supported_fraction =
      n_supported_rows /
        n_raw_rows,

    raw_prob_sum =
      sum(
        raw_prob,
        na.rm =
          TRUE
      ),

    supported_prob_sum =
      sum(
        supported_prob,
        na.rm =
          TRUE
      ),

    .groups =
      "drop"
  )


write.csv(
  support_audit,
  file.path(
    TAB_OUT,
    "02_sample_direction_support_audit_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. LR-level raw and supported comparisons
# ==============================================================================

LR_KEYS <- c(
  "direction",
  "source",
  "target",
  "pathway_name",
  "interaction_key",
  "ligand",
  "receptor"
)


raw_lr_summary <- four_sample_summary(
  lr_all,
  key_cols =
    LR_KEYS,
  value_col =
    "raw_prob"
) %>%
  rename_with(
    ~ paste0(
      "raw_",
      .x
    ),
    c(
      "Sham1",
      "Sham20",
      "Tx17",
      "Tx5",
      "pairwise_up_n",
      "pairwise_down_n",
      "pairwise_tie_n",
      "pairwise_direction",
      "pairwise_direction_consistency",
      "both_Tx_above_both_Sham",
      "both_Tx_below_both_Sham",
      "evidence_grade",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "detected_Sham_replicates",
      "detected_Tx_replicates",
      "detected_samples"
    )
  )


supported_lr_summary <- four_sample_summary(
  lr_all,
  key_cols =
    LR_KEYS,
  value_col =
    "supported_prob"
) %>%
  rename_with(
    ~ paste0(
      "supported_",
      .x
    ),
    c(
      "Sham1",
      "Sham20",
      "Tx17",
      "Tx5",
      "pairwise_up_n",
      "pairwise_down_n",
      "pairwise_tie_n",
      "pairwise_direction",
      "pairwise_direction_consistency",
      "both_Tx_above_both_Sham",
      "both_Tx_below_both_Sham",
      "evidence_grade",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "detected_Sham_replicates",
      "detected_Tx_replicates",
      "detected_samples"
    )
  )


support_counts_lr <- lr_all %>%
  group_by(
    across(
      all_of(
        LR_KEYS
      )
    )
  ) %>%
  summarise(
    Sham_support_replicates =
      n_distinct(
        sample[
          condition ==
            "Sham" &
            supported
        ]
      ),

    Tx_support_replicates =
      n_distinct(
        sample[
          condition ==
            "Tx" &
            supported
        ]
      ),

    total_support_samples =
      n_distinct(
        sample[
          supported
        ]
      ),

    min_supported_pval =
      suppressWarnings(
        min(
          cellchat_pval[
            supported
          ],
          na.rm =
            TRUE
        )
      ),

    .groups =
      "drop"
  ) %>%
  mutate(
    min_supported_pval =
      ifelse(
        is.infinite(
          min_supported_pval
        ),
        NA_real_,
        min_supported_pval
      )
  )


lr_refined <- raw_lr_summary %>%
  left_join(
    supported_lr_summary,
    by =
      LR_KEYS
  ) %>%
  left_join(
    support_counts_lr,
    by =
      LR_KEYS
  ) %>%
  mutate(
    priority_family =
      priority_family_flag(
        pathway_name
      ),

    core_family =
      core_family_flag(
        pathway_name
      ),

    support_pattern =
      case_when(
        Sham_support_replicates ==
          2 &
          Tx_support_replicates ==
            2 ~
          "Supported_2v2",

        Sham_support_replicates ==
          2 &
          Tx_support_replicates ==
            0 ~
          "Sham_only_2v0",

        Sham_support_replicates ==
          0 &
          Tx_support_replicates ==
            2 ~
          "Tx_only_0v2",

        Sham_support_replicates >=
          1 &
          Tx_support_replicates >=
            1 ~
          "Supported_both_conditions",

        Sham_support_replicates >
          0 &
          Tx_support_replicates ==
            0 ~
          "Sham_only_partial",

        Sham_support_replicates ==
          0 &
          Tx_support_replicates >
            0 ~
          "Tx_only_partial",

        TRUE ~
          "Unsupported"
      )
  )


write.csv(
  lr_refined,
  file.path(
    TAB_OUT,
    "03_LR_refined_replicate_comparison_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Pathway sample-level raw/supported values
# ==============================================================================

pathway_sample <- lr_all %>%
  group_by(
    sample,
    condition,
    direction,
    source,
    target,
    pathway_name
  ) %>%
  summarise(
    raw_pathway_strength =
      sum(
        raw_prob,
        na.rm =
          TRUE
      ),

    supported_pathway_strength =
      sum(
        supported_prob,
        na.rm =
          TRUE
      ),

    n_raw_LR =
      n_distinct(
        interaction_key
      ),

    n_supported_LR =
      n_distinct(
        interaction_key[
          supported
        ]
      ),

    .groups =
      "drop"
  )


write.csv(
  pathway_sample,
  file.path(
    TAB_OUT,
    "04_pathway_raw_supported_by_sample_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Pathway-level refined comparison
# ==============================================================================

PATH_KEYS <- c(
  "direction",
  "source",
  "target",
  "pathway_name"
)


raw_path_summary <- four_sample_summary(
  pathway_sample,
  key_cols =
    PATH_KEYS,
  value_col =
    "raw_pathway_strength"
) %>%
  rename_with(
    ~ paste0(
      "raw_",
      .x
    ),
    c(
      "Sham1",
      "Sham20",
      "Tx17",
      "Tx5",
      "pairwise_up_n",
      "pairwise_down_n",
      "pairwise_tie_n",
      "pairwise_direction",
      "pairwise_direction_consistency",
      "both_Tx_above_both_Sham",
      "both_Tx_below_both_Sham",
      "evidence_grade",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "detected_Sham_replicates",
      "detected_Tx_replicates",
      "detected_samples"
    )
  )


supported_path_summary <- four_sample_summary(
  pathway_sample,
  key_cols =
    PATH_KEYS,
  value_col =
    "supported_pathway_strength"
) %>%
  rename_with(
    ~ paste0(
      "supported_",
      .x
    ),
    c(
      "Sham1",
      "Sham20",
      "Tx17",
      "Tx5",
      "pairwise_up_n",
      "pairwise_down_n",
      "pairwise_tie_n",
      "pairwise_direction",
      "pairwise_direction_consistency",
      "both_Tx_above_both_Sham",
      "both_Tx_below_both_Sham",
      "evidence_grade",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "detected_Sham_replicates",
      "detected_Tx_replicates",
      "detected_samples"
    )
  )


path_support_counts <- pathway_sample %>%
  group_by(
    across(
      all_of(
        PATH_KEYS
      )
    )
  ) %>%
  summarise(
    Sham_support_replicates =
      sum(
        condition ==
          "Sham" &
          supported_pathway_strength >
            0
      ),

    Tx_support_replicates =
      sum(
        condition ==
          "Tx" &
          supported_pathway_strength >
            0
      ),

    total_support_samples =
      sum(
        supported_pathway_strength >
          0
      ),

    .groups =
      "drop"
  )


pathway_refined <- raw_path_summary %>%
  left_join(
    supported_path_summary,
    by =
      PATH_KEYS
  ) %>%
  left_join(
    path_support_counts,
    by =
      PATH_KEYS
  ) %>%
  mutate(
    priority_family =
      priority_family_flag(
        pathway_name
      ),

    core_family =
      core_family_flag(
        pathway_name
      ),

    support_pattern =
      case_when(
        Sham_support_replicates ==
          2 &
          Tx_support_replicates ==
            2 ~
          "Supported_2v2",

        Sham_support_replicates ==
          2 &
          Tx_support_replicates ==
            0 ~
          "Sham_only_2v0",

        Sham_support_replicates ==
          0 &
          Tx_support_replicates ==
            2 ~
          "Tx_only_0v2",

        Sham_support_replicates >=
          1 &
          Tx_support_replicates >=
            1 ~
          "Supported_both_conditions",

        Sham_support_replicates >
          0 &
          Tx_support_replicates ==
            0 ~
          "Sham_only_partial",

        Sham_support_replicates ==
          0 &
          Tx_support_replicates >
            0 ~
          "Tx_only_partial",

        TRUE ~
          "Unsupported"
      )
  )


write.csv(
  pathway_refined,
  file.path(
    TAB_OUT,
    "05_pathway_refined_replicate_comparison_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Priority family tables
# ==============================================================================

priority_lr <- lr_refined %>%
  filter(
    priority_family
  ) %>%
  arrange(
    direction,
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


priority_pathway <- pathway_refined %>%
  filter(
    priority_family
  ) %>%
  arrange(
    direction,
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


write.csv(
  priority_lr,
  file.path(
    TAB_OUT,
    "06_priority_LR_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


write.csv(
  priority_pathway,
  file.path(
    TAB_OUT,
    "07_priority_pathway_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Mphi -> Injury/inflammatory Hepatocyte focus
# ==============================================================================

mphi_to_injury_lr <- lr_refined %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


mphi_to_injury_path <- pathway_refined %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


write.csv(
  mphi_to_injury_lr,
  file.path(
    TAB_OUT,
    "08_Mphi5_to_InjuryHep_LR_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


write.csv(
  mphi_to_injury_path,
  file.path(
    TAB_OUT,
    "09_Mphi5_to_InjuryHep_pathway_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Injury/inflammatory Hepatocyte -> Mphi focus
# ==============================================================================

injury_to_mphi_lr <- lr_refined %>%
  filter(
    direction ==
      "Hep_to_Mphi",
    source ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


injury_to_mphi_path <- pathway_refined %>%
  filter(
    direction ==
      "Hep_to_Mphi",
    source ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


write.csv(
  injury_to_mphi_lr,
  file.path(
    TAB_OUT,
    "10_InjuryHep_to_Mphi5_LR_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


write.csv(
  injury_to_mphi_path,
  file.path(
    TAB_OUT,
    "11_InjuryHep_to_Mphi5_pathway_refined_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Core pathway panel
# ==============================================================================

CORE_FAMILIES <- c(
  "PDGF",
  "FN1",
  "SPP1",
  "TNF",
  "SEMA4",
  "MIF",
  "COLLAGEN",
  "LAMININ"
)


core_pathway_panel <- pathway_refined %>%
  filter(
    toupper(
      pathway_name
    ) %in%
      toupper(
        CORE_FAMILIES
      )
  ) %>%
  arrange(
    direction,
    pathway_name,
    source,
    target
  )


write.csv(
  core_pathway_panel,
  file.path(
    TAB_OUT,
    "12_core_pathway_panel_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. High-confidence shortlist
# ==============================================================================

# High-confidence rule:
#   - core or priority biological family
#   - support in at least 3/4 biological samples
#   - supported direction consistency >= 0.75
#   - non-zero supported effect
#
# "Strong" additionally requires both Tx above both Sham OR both Tx below both
# Sham on supported probability.

high_conf_pathway <- pathway_refined %>%
  filter(
    priority_family,
    total_support_samples >=
      3,
    supported_pairwise_direction_consistency >=
      0.75,
    supported_Tx_minus_Sham !=
      0
  ) %>%
  mutate(
    confidence_class =
      case_when(
        supported_both_Tx_above_both_Sham |
          supported_both_Tx_below_both_Sham ~
          "High_Strong",

        TRUE ~
          "High_Moderate"
      )
  ) %>%
  arrange(
    desc(
      core_family
    ),
    confidence_class,
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


high_conf_lr <- lr_refined %>%
  filter(
    priority_family,
    total_support_samples >=
      3,
    supported_pairwise_direction_consistency >=
      0.75,
    supported_Tx_minus_Sham !=
      0
  ) %>%
  mutate(
    confidence_class =
      case_when(
        supported_both_Tx_above_both_Sham |
          supported_both_Tx_below_both_Sham ~
          "High_Strong",

        TRUE ~
          "High_Moderate"
      )
  ) %>%
  arrange(
    desc(
      core_family
    ),
    confidence_class,
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


write.csv(
  high_conf_pathway,
  file.path(
    TAB_OUT,
    "13_high_confidence_pathway_shortlist_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


write.csv(
  high_conf_lr,
  file.path(
    TAB_OUT,
    "14_high_confidence_LR_shortlist_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. NicheNet candidate shortlist
# ==============================================================================

nichenet_candidates <- high_conf_lr %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  select(
    pathway_name,
    source,
    target,
    ligand,
    receptor,
    interaction_key,
    Sham_support_replicates,
    Tx_support_replicates,
    total_support_samples,
    supported_Sham_mean,
    supported_Tx_mean,
    supported_Tx_minus_Sham,
    supported_pairwise_direction,
    supported_pairwise_direction_consistency,
    supported_evidence_grade,
    confidence_class,
    core_family
  ) %>%
  arrange(
    desc(
      core_family
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    ),
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    )
  )


write.csv(
  nichenet_candidates,
  file.path(
    TAB_OUT,
    "15_NicheNet_candidate_Mphi_to_InjuryHep_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. PDGF / FN1 / SPP1 / TNF / SEMA4 dedicated table
# ==============================================================================

MECHANISTIC_FAMILIES <- c(
  "PDGF",
  "FN1",
  "SPP1",
  "TNF",
  "SEMA4"
)


mechanistic_core <- pathway_refined %>%
  filter(
    toupper(
      pathway_name
    ) %in%
      toupper(
        MECHANISTIC_FAMILIES
      )
  ) %>%
  arrange(
    pathway_name,
    direction,
    source,
    target
  )


write.csv(
  mechanistic_core,
  file.path(
    TAB_OUT,
    "16_PDGFFN1SPP1TNFSEMA4_pathway_matrix_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Figure 1: Mphi -> Injury-Hep refined pathways
# ==============================================================================

plot_mphi_injury <- mphi_to_injury_path %>%
  filter(
    priority_family,
    total_support_samples >=
      2
  ) %>%
  mutate(
    label =
      paste0(
        source,
        " | ",
        pathway_name
      )
  ) %>%
  arrange(
    supported_Tx_minus_Sham
  )


if (
  nrow(
    plot_mphi_injury
  )
) {

  p1 <- ggplot(
    plot_mphi_injury,
    aes(
      x =
        supported_Tx_minus_Sham,
      y =
        reorder(
          label,
          supported_Tx_minus_Sham
        ),
      size =
        total_support_samples
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.35
    ) +
    labs(
      title =
        "Mphi -> Injury/inflammatory Hepatocyte",
      subtitle =
        "Supported CellChat probability only (p < 0.05)",
      x =
        "Tx mean - Sham mean supported strength",
      y =
        NULL,
      size =
        "Supported samples"
    ) +
    theme_classic(
      base_size =
        7.5
    )

  save_pdf(
    p1,
    file.path(
      FIG_OUT,
      "01_Mphi5_to_InjuryHep_supported_pathways_v6.6.2.1.pdf"
    ),
    10,
    13
  )
}


# ==============================================================================
# 17. Figure 2: Injury-Hep -> Mphi refined pathways
# ==============================================================================

plot_injury_mphi <- injury_to_mphi_path %>%
  filter(
    priority_family,
    total_support_samples >=
      2
  ) %>%
  mutate(
    label =
      paste0(
        pathway_name,
        " | ",
        target
      )
  ) %>%
  arrange(
    supported_Tx_minus_Sham
  )


if (
  nrow(
    plot_injury_mphi
  )
) {

  p2 <- ggplot(
    plot_injury_mphi,
    aes(
      x =
        supported_Tx_minus_Sham,
      y =
        reorder(
          label,
          supported_Tx_minus_Sham
        ),
      size =
        total_support_samples
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.35
    ) +
    labs(
      title =
        "Injury/inflammatory Hepatocyte -> Mphi",
      subtitle =
        "Supported CellChat probability only (p < 0.05)",
      x =
        "Tx mean - Sham mean supported strength",
      y =
        NULL,
      size =
        "Supported samples"
    ) +
    theme_classic(
      base_size =
        7.5
    )

  save_pdf(
    p2,
    file.path(
      FIG_OUT,
      "02_InjuryHep_to_Mphi5_supported_pathways_v6.6.2.1.pdf"
    ),
    10,
    14
  )
}


# ==============================================================================
# 18. Figure 3: core pathway sample heatmap
# ==============================================================================

core_sample <- pathway_sample %>%
  filter(
    toupper(
      pathway_name
    ) %in%
      toupper(
        CORE_FAMILIES
      )
  ) %>%
  mutate(
    interaction =
      paste0(
        direction,
        " | ",
        source,
        " -> ",
        target,
        " | ",
        pathway_name
      )
  )


top_interactions <- core_sample %>%
  group_by(
    interaction
  ) %>%
  summarise(
    max_supported =
      max(
        supported_pathway_strength,
        na.rm =
          TRUE
      ),
    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      max_supported
    )
  ) %>%
  slice_head(
    n =
      60
  ) %>%
  pull(
    interaction
  )


heat_df <- core_sample %>%
  filter(
    interaction %in%
      top_interactions
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    interaction =
      factor(
        interaction,
        levels =
          rev(
            top_interactions
          )
      )
  )


if (
  nrow(
    heat_df
  )
) {

  p3 <- ggplot(
    heat_df,
    aes(
      x =
        sample,
      y =
        interaction,
      fill =
        supported_pathway_strength
    )
  ) +
    geom_tile(
      linewidth =
        0.2
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
        "Core Mphi <-> Hepatocyte pathways",
      subtitle =
        "Sample-wise supported CellChat strength",
      x =
        NULL,
      y =
        NULL,
      fill =
        "Supported\nstrength"
    ) +
    theme_classic(
      base_size =
        6.2
    )

  save_pdf(
    p3,
    file.path(
      FIG_OUT,
      "03_core_pathways_supported_by_sample_v6.6.2.1.pdf"
    ),
    10,
    16
  )
}


# ==============================================================================
# 19. Figure 4: mechanistic core Tx-Sham matrix
# ==============================================================================

matrix_df <- mechanistic_core %>%
  filter(
    total_support_samples >=
      2
  ) %>%
  mutate(
    interaction =
      paste0(
        source,
        " -> ",
        target
      ),
    pathway_name =
      factor(
        pathway_name,
        levels =
          MECHANISTIC_FAMILIES
      )
  )


if (
  nrow(
    matrix_df
  )
) {

  p4 <- ggplot(
    matrix_df,
    aes(
      x =
        interaction,
      y =
        pathway_name,
      fill =
        supported_Tx_minus_Sham
    )
  ) +
    geom_tile(
      linewidth =
        0.25
    ) +
    scale_fill_gradient2(
      low =
        "#0033FF",
      mid =
        "#FFFFFF",
      high =
        "#FF1A1A",
      midpoint =
        0
    ) +
    labs(
      title =
        "Mechanistic core communication rewiring",
      subtitle =
        "Supported CellChat Tx - Sham",
      x =
        NULL,
      y =
        NULL,
      fill =
        "Tx-Sham"
    ) +
    theme_classic(
      base_size =
        6.5
    ) +
    theme(
      axis.text.x =
        element_text(
          angle =
            60,
          hjust =
            1
        )
    )

  save_pdf(
    p4,
    file.path(
      FIG_OUT,
      "04_PDGFFN1SPP1TNFSEMA4_supported_TxSham_matrix_v6.6.2.1.pdf"
    ),
    18,
    6
  )
}


# ==============================================================================
# 20. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_CellChat_directory",
    "samples",
    "CellChat_recomputed",
    "support_threshold",
    "population_size_setting_in_source_CellChat",
    "source_CellChat_estimator",
    "primary_output",
    "high_confidence_rule",
    "biological_replicates"
  ),
  value = c(
    "v6.6.2.1",
    IN_DIR,
    paste(
      SAMPLES,
      collapse = ","
    ),
    "FALSE",
    "CellChat pval < 0.05",
    "FALSE",
    "triMean",
    "Supported probability + biological-replicate consistency",
    "Priority family; >=3/4 supported samples; pairwise consistency >=0.75; non-zero effect",
    "Sham n=2; Tx n=2"
  )
)


write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.2.1.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.2.1.txt"
    )
)


# ==============================================================================
# 21. Console summary
# ==============================================================================

msg(
  "Sample-level support audit:"
)

print(
  support_audit
)


msg(
  "High-confidence Mphi -> Injury-Hep pathways:"
)

print(
  high_conf_pathway %>%
    filter(
      direction ==
        "Mphi_to_Hep",
      target ==
        "Hep_Injury-inflammatory"
    ) %>%
    select(
      source,
      target,
      pathway_name,
      Sham_support_replicates,
      Tx_support_replicates,
      supported_Sham_mean,
      supported_Tx_mean,
      supported_Tx_minus_Sham,
      supported_pairwise_direction_consistency,
      supported_evidence_grade,
      confidence_class
    )
)


msg(
  "High-confidence Injury-Hep -> Mphi pathways:"
)

print(
  high_conf_pathway %>%
    filter(
      direction ==
        "Hep_to_Mphi",
      source ==
        "Hep_Injury-inflammatory"
    ) %>%
    select(
      source,
      target,
      pathway_name,
      Sham_support_replicates,
      Tx_support_replicates,
      supported_Sham_mean,
      supported_Tx_mean,
      supported_Tx_minus_Sham,
      supported_pairwise_direction_consistency,
      supported_evidence_grade,
      confidence_class
    )
)


msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
