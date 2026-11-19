#!/usr/bin/env Rscript

# ============================================================
# Mouse MASH Mphi -> HSC integrated mechanistic closure v6.4.5
#
# PURPOSE
#   Integrate:
#     - Mouse macrophage ligand burden from v6.4.1
#     - HSC 3-axis host-weighted output from v6.4.4
#
#   Focus ligands:
#     Spp1
#     Cdcp1
#     Tgfb1
#
#   Sender resolution:
#     all five primary macrophage states
#     + total primary-myeloid burden
#
#   Receiver endpoints:
#     Whole-HSC host-weighted:
#       Core_collagen
#       ECM_remodeling
#       IEG_core
#
#     ECM-activated HSC host-weighted:
#       Core_collagen
#       ECM_remodeling
#       IEG_core
#
#     qHSC host-weighted:
#       Core_collagen
#       ECM_remodeling
#       IEG_core
#
# IMPORTANT
#   Sham n=2, Tx n=2.
#   All group effects and n=4 correlations are descriptive.
#   No causal inference.
#   No reclustering.
#   No source RDS overwrite.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RESULT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1",
  "objects",
  "Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1_results.rds"
)

HSC_RESULT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4",
  "objects",
  "Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4_results.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_HSC_IntegratedClosure_v6.4.5"
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
  "Mouse_MASH_Mphi_HSC_IntegratedClosure_v6.4.5.log"
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
cat("Mouse MASH Mphi -> HSC integrated closure v6.4.5\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(MPHI_RESULT_RDS),
  file.exists(HSC_RESULT_RDS)
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

FOCUS_LIGANDS <- c(
  "Spp1",
  "Cdcp1",
  "Tgfb1"
)

PRIMARY_MPHI_STATES <- c(
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Inflammatory-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Repair/Resolution-Mphi"
)

TREM2_STATE <- "Lipid-associated/TREM2-Mphi"

HSC_PROGRAMS <- c(
  "Core_collagen",
  "ECM_remodeling",
  "IEG_core"
)

HSC_STATES_FOCUS <- c(
  "ECM-activated HSC",
  "qHSC"
)

# ------------------------------------------------------------
# 2. Load prior validated result objects
# ------------------------------------------------------------

mres <- readRDS(
  MPHI_RESULT_RDS
)

hres <- readRDS(
  HSC_RESULT_RDS
)

required_mres <- c(
  "burden",
  "mphi_sample_inventory"
)

required_hres <- c(
  "state3",
  "HSC_abundance"
)

if (!all(
  required_mres %in%
    names(mres)
)) {
  stop(
    "v6.4.1 results object missing required elements."
  )
}

if (!all(
  required_hres %in%
    names(hres)
)) {
  stop(
    "v6.4.4 results object missing required elements."
  )
}

burden <- mres$burden
state3 <- hres$state3
hsc_abundance <- hres$HSC_abundance

# ------------------------------------------------------------
# 3. Integrity checks
# ------------------------------------------------------------

if (!all(
  FOCUS_LIGANDS %in%
    unique(
      burden$ligand
    )
)) {
  stop(
    "Not all focus ligands are present in v6.4.1 burden table."
  )
}

if (!all(
  PRIMARY_MPHI_STATES %in%
    unique(
      burden$sender_state
    )
)) {
  stop(
    "Not all five primary macrophage states are present."
  )
}

if (!all(
  HSC_PROGRAMS %in%
    unique(
      state3$program
    )
)) {
  stop(
    "Not all three HSC programs are present in v6.4.4."
  )
}

if (!all(
  SAMPLES %in%
    unique(
      burden$sample
    )
) ||
  !all(
    SAMPLES %in%
      unique(
        state3$sample
      )
  )) {
  stop(
    "Not all four biological samples are present."
  )
}

# ------------------------------------------------------------
# 4. Sender table: all five Mphi states
# ------------------------------------------------------------

sender_all <- burden %>%
  filter(
    ligand %in%
      FOCUS_LIGANDS,
    sender_state %in%
      PRIMARY_MPHI_STATES
  ) %>%
  select(
    sample,
    group,
    sender_state,
    ligand,
    state_n,
    total_primary_mphi_n,
    state_fraction,
    ligand_mean,
    ligand_detection,
    ligand_burden_expr,
    ligand_burden_detection,
    ligand_expression_share,
    ligand_positive_cell_share
  ) %>%
  arrange(
    ligand,
    sender_state,
    sample
  )

write_tsv(
  sender_all,
  file.path(
    DIR_TABLE,
    "01_sender_all5_focus_ligand_burden_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 5. Total-primary-myeloid ligand burden
#
# Sum state-specific expression across five primary states,
# divided by total primary Mphi denominator.
# ------------------------------------------------------------

sender_total <- burden %>%
  filter(
    ligand %in%
      FOCUS_LIGANDS,
    sender_state %in%
      PRIMARY_MPHI_STATES
  ) %>%
  group_by(
    sample,
    group,
    ligand
  ) %>%
  summarise(
    total_primary_mphi_n =
      first(
        total_primary_mphi_n
      ),

    total_ligand_expr_sum =
      sum(
        ligand_expr_sum,
        na.rm = TRUE
      ),

    total_ligand_detected_n =
      sum(
        ligand_detected_n,
        na.rm = TRUE
      ),

    total_myeloid_burden_expr =
      total_ligand_expr_sum /
      total_primary_mphi_n,

    total_myeloid_burden_detection =
      total_ligand_detected_n /
      total_primary_mphi_n,

    .groups = "drop"
  ) %>%
  arrange(
    ligand,
    sample
  )

write_tsv(
  sender_total,
  file.path(
    DIR_TABLE,
    "02_total_primary_myeloid_focus_ligand_burden.tsv"
  )
)

# ------------------------------------------------------------
# 6. Sender Sham vs Tx effects
# ------------------------------------------------------------

group_effect <- function(
  df,
  value_col,
  id_cols
) {

  df %>%
    group_by(
      across(
        all_of(
          c(
            id_cols,
            "group"
          )
        )
      )
    ) %>%
    summarise(
      n_samples = n(),
      mean_value =
        mean(
          .data[[value_col]],
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from =
        group,
      values_from = c(
        n_samples,
        mean_value
      )
    ) %>%
    mutate(
      Tx_minus_Sham =
        mean_value_Tx -
        mean_value_Sham,

      Tx_over_Sham =
        ifelse(
          mean_value_Sham > 0,
          mean_value_Tx /
            mean_value_Sham,
          NA_real_
        ),

      log2_Tx_over_Sham =
        ifelse(
          mean_value_Sham > 0 &
            mean_value_Tx > 0,
          log2(
            mean_value_Tx /
              mean_value_Sham
          ),
          NA_real_
        ),

      observed_direction =
        case_when(
          Tx_minus_Sham < 0 ~
            "decrease",
          Tx_minus_Sham > 0 ~
            "increase",
          TRUE ~
            "no_change"
        )
    )
}

sender_state_effect <- group_effect(
  sender_all,
  "ligand_burden_expr",
  c(
    "ligand",
    "sender_state"
  )
) %>%
  arrange(
    ligand,
    Tx_over_Sham
  )

sender_total_effect <- group_effect(
  sender_total,
  "total_myeloid_burden_expr",
  c(
    "ligand"
  )
) %>%
  arrange(
    ligand
  )

write_tsv(
  sender_state_effect,
  file.path(
    DIR_TABLE,
    "03_sender_state_ligand_burden_Sham_vs_Tx.tsv"
  )
)

write_tsv(
  sender_total_effect,
  file.path(
    DIR_TABLE,
    "04_total_myeloid_ligand_burden_Sham_vs_Tx.tsv"
  )
)

# ------------------------------------------------------------
# 7. Receiver endpoints per biological sample
# ------------------------------------------------------------

whole_hsc <- state3 %>%
  group_by(
    sample,
    group,
    program
  ) %>%
  summarise(
    whole_HSC_host_weighted_burden =
      sum(
        host_weighted_burden,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    receiver_endpoint =
      paste0(
        "WholeHSC_",
        program
      )
  ) %>%
  select(
    sample,
    group,
    receiver_endpoint,
    receiver_value =
      whole_HSC_host_weighted_burden
  )

state_hsc <- state3 %>%
  filter(
    state %in%
      HSC_STATES_FOCUS
  ) %>%
  transmute(
    sample,
    group,
    receiver_endpoint =
      paste0(
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          state
        ),
        "_",
        program
      ),
    receiver_value =
      host_weighted_burden
  )

abundance_endpoint <- hsc_abundance %>%
  transmute(
    sample,
    group,
    receiver_endpoint =
      "HSC_abundance",
    receiver_value =
      HSC_fraction_of_recovered_host
  )

receiver_long <- bind_rows(
  whole_hsc,
  state_hsc,
  abundance_endpoint
) %>%
  arrange(
    receiver_endpoint,
    sample
  )

write_tsv(
  receiver_long,
  file.path(
    DIR_TABLE,
    "05_receiver_host_weighted_endpoints_by_sample.tsv"
  )
)

receiver_effect <- group_effect(
  receiver_long,
  "receiver_value",
  c(
    "receiver_endpoint"
  )
) %>%
  arrange(
    receiver_endpoint
  )

write_tsv(
  receiver_effect,
  file.path(
    DIR_TABLE,
    "06_receiver_endpoints_Sham_vs_Tx.tsv"
  )
)

# ------------------------------------------------------------
# 8. Correlation helper
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 9. All-state sender -> receiver descriptive correlations
# ------------------------------------------------------------

sender_receiver_state <- sender_all %>%
  select(
    sample,
    group,
    ligand,
    sender_state,
    sender_value =
      ligand_burden_expr
  ) %>%
  inner_join(
    receiver_long,
    by = c(
      "sample",
      "group"
    ),
    relationship =
      "many-to-many"
  ) %>%
  group_by(
    ligand,
    sender_state,
    receiver_endpoint
  ) %>%
  summarise(
    n_samples =
      sum(
        complete.cases(
          sender_value,
          receiver_value
        )
      ),

    rho_descriptive =
      safe_rho(
        sender_value,
        receiver_value
      ),

    inference =
      "descriptive_only_n4",

    .groups = "drop"
  ) %>%
  arrange(
    receiver_endpoint,
    desc(
      abs(
        rho_descriptive
      )
    ),
    ligand,
    sender_state
  )

write_tsv(
  sender_receiver_state,
  file.path(
    DIR_TABLE,
    "07_all5_sender_to_receiver_correlations_descriptive_n4.tsv"
  )
)

# ------------------------------------------------------------
# 10. Total-myeloid sender -> receiver correlations
# ------------------------------------------------------------

sender_receiver_total <- sender_total %>%
  select(
    sample,
    group,
    ligand,
    sender_value =
      total_myeloid_burden_expr
  ) %>%
  inner_join(
    receiver_long,
    by = c(
      "sample",
      "group"
    ),
    relationship =
      "many-to-many"
  ) %>%
  group_by(
    ligand,
    receiver_endpoint
  ) %>%
  summarise(
    n_samples =
      sum(
        complete.cases(
          sender_value,
          receiver_value
        )
      ),

    rho_descriptive =
      safe_rho(
        sender_value,
        receiver_value
      ),

    inference =
      "descriptive_only_n4",

    .groups = "drop"
  ) %>%
  arrange(
    receiver_endpoint,
    desc(
      abs(
        rho_descriptive
      )
    ),
    ligand
  )

write_tsv(
  sender_receiver_total,
  file.path(
    DIR_TABLE,
    "08_total_myeloid_to_receiver_correlations_descriptive_n4.tsv"
  )
)

# ------------------------------------------------------------
# 11. Focused mechanistic matrix
#
# Primary receiver endpoint:
#   ECM-activated HSC Core_collagen
#
# Secondary:
#   WholeHSC Core_collagen
#   ECM-activated HSC ECM_remodeling
#   WholeHSC ECM_remodeling
# ------------------------------------------------------------

PRIMARY_RECEIVERS <- c(
  "ECM_activated_HSC_Core_collagen",
  "WholeHSC_Core_collagen",
  "ECM_activated_HSC_ECM_remodeling",
  "WholeHSC_ECM_remodeling"
)

focused_state <- sender_state_effect %>%
  select(
    ligand,
    sender_state,
    sender_Sham =
      mean_value_Sham,
    sender_Tx =
      mean_value_Tx,
    sender_Tx_minus_Sham =
      Tx_minus_Sham,
    sender_Tx_over_Sham =
      Tx_over_Sham,
    sender_direction =
      observed_direction
  ) %>%
  crossing(
    receiver_endpoint =
      PRIMARY_RECEIVERS
  ) %>%
  left_join(
    receiver_effect %>%
      filter(
        receiver_endpoint %in%
          PRIMARY_RECEIVERS
      ) %>%
      select(
        receiver_endpoint,
        receiver_Sham =
          mean_value_Sham,
        receiver_Tx =
          mean_value_Tx,
        receiver_Tx_minus_Sham =
          Tx_minus_Sham,
        receiver_Tx_over_Sham =
          Tx_over_Sham,
        receiver_direction =
          observed_direction
      ),
    by =
      "receiver_endpoint"
  ) %>%
  left_join(
    sender_receiver_state %>%
      filter(
        receiver_endpoint %in%
          PRIMARY_RECEIVERS
      ) %>%
      select(
        ligand,
        sender_state,
        receiver_endpoint,
        rho_descriptive
      ),
    by = c(
      "ligand",
      "sender_state",
      "receiver_endpoint"
    )
  ) %>%
  mutate(
    expected_mechanistic_direction =
      "sender_decrease_and_receiver_decrease",

    both_decrease =
      sender_direction ==
        "decrease" &
      receiver_direction ==
        "decrease",

    positive_cross_sample_concordance =
      !is.na(
        rho_descriptive
      ) &
      rho_descriptive > 0,

    closure_support =
      case_when(
        both_decrease &
          positive_cross_sample_concordance ~
          "direction_plus_positive_rho",

        both_decrease ~
          "direction_only",

        TRUE ~
          "not_supported"
      )
  ) %>%
  arrange(
    receiver_endpoint,
    desc(
      closure_support
    ),
    desc(
      rho_descriptive
    )
  )

write_tsv(
  focused_state,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_focused_mechanistic_matrix.tsv"
  )
)

# Total-myeloid version
focused_total <- sender_total_effect %>%
  select(
    ligand,
    sender_Sham =
      mean_value_Sham,
    sender_Tx =
      mean_value_Tx,
    sender_Tx_minus_Sham =
      Tx_minus_Sham,
    sender_Tx_over_Sham =
      Tx_over_Sham,
    sender_direction =
      observed_direction
  ) %>%
  crossing(
    receiver_endpoint =
      PRIMARY_RECEIVERS
  ) %>%
  left_join(
    receiver_effect %>%
      filter(
        receiver_endpoint %in%
          PRIMARY_RECEIVERS
      ) %>%
      select(
        receiver_endpoint,
        receiver_Sham =
          mean_value_Sham,
        receiver_Tx =
          mean_value_Tx,
        receiver_Tx_minus_Sham =
          Tx_minus_Sham,
        receiver_Tx_over_Sham =
          Tx_over_Sham,
        receiver_direction =
          observed_direction
      ),
    by =
      "receiver_endpoint"
  ) %>%
  left_join(
    sender_receiver_total %>%
      filter(
        receiver_endpoint %in%
          PRIMARY_RECEIVERS
      ) %>%
      select(
        ligand,
        receiver_endpoint,
        rho_descriptive
      ),
    by = c(
      "ligand",
      "receiver_endpoint"
    )
  ) %>%
  mutate(
    both_decrease =
      sender_direction ==
        "decrease" &
      receiver_direction ==
        "decrease",

    positive_cross_sample_concordance =
      !is.na(
        rho_descriptive
      ) &
      rho_descriptive > 0,

    closure_support =
      case_when(
        both_decrease &
          positive_cross_sample_concordance ~
          "direction_plus_positive_rho",

        both_decrease ~
          "direction_only",

        TRUE ~
          "not_supported"
      )
  ) %>%
  arrange(
    receiver_endpoint,
    desc(
      closure_support
    ),
    desc(
      rho_descriptive
    )
  )

write_tsv(
  focused_total,
  file.path(
    DIR_TABLE,
    "10_TOTAL_MYELOID_focused_mechanistic_matrix.tsv"
  )
)

# ------------------------------------------------------------
# 12. Rank treatment-decreased sender axes
#
# This is not an inferential ranking.
# It identifies which ligand/state combinations move in the
# treatment-expected direction most strongly.
# ------------------------------------------------------------

sender_decrease_rank <- sender_state_effect %>%
  mutate(
    pct_change =
      (
        Tx_over_Sham -
          1
      ) *
      100
  ) %>%
  arrange(
    Tx_over_Sham
  )

write_tsv(
  sender_decrease_rank,
  file.path(
    DIR_TABLE,
    "11_sender_state_treatment_direction_rank.tsv"
  )
)

# ------------------------------------------------------------
# 13. Human-prioritized TREM2 axis summary
# ------------------------------------------------------------

trem2_summary <- focused_state %>%
  filter(
    sender_state ==
      TREM2_STATE
  ) %>%
  arrange(
    ligand,
    receiver_endpoint
  )

write_tsv(
  trem2_summary,
  file.path(
    DIR_TABLE,
    "12_HUMAN_PRIOR_TREM2_axis_summary.tsv"
  )
)

# ------------------------------------------------------------
# 14. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  sender_all,
  aes(
    x = group,
    y =
      ligand_burden_expr
  )
) +
  geom_point(
    size = 2.5,
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
    size = 2.5,
    nudge_y = 0.005,
    check_overlap = TRUE
  ) +
  facet_grid(
    ligand ~ sender_state,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 8
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
      "Focus-ligand burden across five primary macrophage states",
    x = NULL,
    y =
      "Abundance-weighted ligand burden"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_focus_ligand_burden_all5_states.pdf"
  ),
  p1,
  width = 15,
  height = 8
)

receiver_plot <- receiver_long %>%
  filter(
    receiver_endpoint %in%
      c(
        "WholeHSC_Core_collagen",
        "ECM_activated_HSC_Core_collagen",
        "WholeHSC_ECM_remodeling",
        "ECM_activated_HSC_ECM_remodeling",
        "WholeHSC_IEG_core",
        "ECM_activated_HSC_IEG_core"
      )
  )

p2 <- ggplot(
  receiver_plot,
  aes(
    x = group,
    y =
      receiver_value
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
    size = 3,
    nudge_y = 0.002,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ receiver_endpoint,
    scales = "free_y"
  ) +
  theme_classic(
    base_size = 10
  ) +
  labs(
    title =
      "Host-weighted HSC receiver endpoints",
    x = NULL,
    y =
      "Host-weighted burden"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_HSC_receiver_endpoints.pdf"
  ),
  p2,
  width = 11,
  height = 7
)

# ------------------------------------------------------------
# 15. Save compact results
# ------------------------------------------------------------

saveRDS(
  list(
    sender_all =
      sender_all,
    sender_total =
      sender_total,
    sender_state_effect =
      sender_state_effect,
    sender_total_effect =
      sender_total_effect,
    receiver_long =
      receiver_long,
    receiver_effect =
      receiver_effect,
    sender_receiver_state =
      sender_receiver_state,
    sender_receiver_total =
      sender_receiver_total,
    focused_state =
      focused_state,
    focused_total =
      focused_total,
    sender_decrease_rank =
      sender_decrease_rank,
    trem2_summary =
      trem2_summary
  ),
  file.path(
    DIR_OBJ,
    "Mouse_MASH_Mphi_HSC_IntegratedClosure_v6.4.5_results.rds"
  )
)

# ------------------------------------------------------------
# 16. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "Mphi_input_result",
    "HSC_input_result",
    "focus_ligands",
    "primary_Mphi_states_n",
    "samples_n",
    "receiver_programs",
    "receiver_primary_partition",
    "total_myeloid_burden_calculated",
    "all5_sender_correlations_calculated",
    "group_inference",
    "correlation_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    MPHI_RESULT_RDS,
    HSC_RESULT_RDS,
    paste(
      FOCUS_LIGANDS,
      collapse = ";"
    ),
    as.character(
      length(
        PRIMARY_MPHI_STATES
      )
    ),
    as.character(
      length(
        SAMPLES
      )
    ),
    paste(
      HSC_PROGRAMS,
      collapse = ";"
    ),
    "state3 host-weighted burden",
    "TRUE",
    "TRUE",
    "descriptive_only_n2_per_group",
    "descriptive_only_n4",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v6.4.5.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "sender_source",
    "receiver_source",
    "sender_metric_primary",
    "total_myeloid_metric",
    "receiver_metric",
    "primary_receiver_endpoint",
    "interpretation"
  ),
  value = c(
    "v6.4.5",
    "validated v6.4.1 macrophage ligand-burden results",
    "validated v6.4.4 HSC host-weighted 3-axis results",
    "state-specific ligand_burden_expr / total primary-Mphi denominator",
    "sum ligand expression across five primary Mphi states / total primary-Mphi denominator",
    "HSC abundance × HSC state fraction × linear CP10k program intensity",
    "ECM-activated HSC Core_collagen host-weighted burden",
    "cross-species mechanistic closure is descriptive, not causal"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "Mphi_HSC_IntegratedClosure_provenance_v6.4.5.tsv"
  )
)

# ------------------------------------------------------------
# 17. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TOTAL MYELOID LIGAND BURDEN: SHAM vs TX\n")
cat("============================================================\n")
print(
  sender_total_effect
)

cat("\n")
cat("============================================================\n")
cat("TOP TREATMENT-DECREASED SENDER AXES\n")
cat("============================================================\n")
print(
  sender_decrease_rank %>%
    select(
      ligand,
      sender_state,
      mean_value_Sham,
      mean_value_Tx,
      Tx_over_Sham,
      observed_direction
    ) %>%
    slice_head(
      n = 12
    )
)

cat("\n")
cat("============================================================\n")
cat("RECEIVER ENDPOINTS: SHAM vs TX\n")
cat("============================================================\n")
print(
  receiver_effect
)

cat("\n")
cat("============================================================\n")
cat("HUMAN-PRIOR TREM2 AXIS SUMMARY\n")
cat("============================================================\n")
print(
  trem2_summary
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY FOCUSED MECHANISTIC MATRIX\n")
cat("============================================================\n")
print(
  focused_state %>%
    filter(
      closure_support !=
        "not_supported"
    )
)

cat("\n")
cat("============================================================\n")
cat("TOTAL-MYELOID FOCUSED MECHANISTIC MATRIX\n")
cat("============================================================\n")
print(
  focused_total
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
