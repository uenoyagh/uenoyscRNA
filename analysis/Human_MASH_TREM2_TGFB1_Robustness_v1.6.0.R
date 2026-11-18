#!/usr/bin/env Rscript

# ============================================================
# Human MASH TREM2-Mphi TGFB1 robustness audit v1.6.0
#
# PURPOSE
#   Robustness-test the focused association identified in v1.5:
#
#     Lipid-associated/TREM2_Mphi TGFB1 burden
#                 ->
#     HSC TGF-beta direct-response score
#
#   Analyses:
#     1. Reproduce primary Spearman association
#     2. Exact permutation test (8! = 40,320 permutations)
#     3. Leave-one-sample-out analysis
#     4. Leave-one-NAS-group-out analysis
#     5. Explicit NAS6-excluded sensitivity analysis
#     6. Partial Spearman controlling for NAS
#     7. Decompose burden into:
#          - state abundance
#          - within-state TGFB1 mean expression
#          - abundance-weighted TGFB1 burden
#     8. Compare TREM2 state with other myeloid states
#     9. Compare against total-myeloid TGFB1 burden
#
# IMPORTANT
#   - Biological sample is the unit of association.
#   - No causal inference is claimed.
#   - No reclustering.
#   - No source RDS is overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

V15 <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5"
)

SENDER_TABLE <- file.path(
  V15,
  "tables",
  "02_TGFB1_sender_state_burden_by_sample.tsv"
)

HSC_TABLE <- file.path(
  V15,
  "tables",
  "01_HSC_endpoints_by_sample.tsv"
)

TOTAL_TGFB1_TABLE <- file.path(
  V15,
  "tables",
  "04_total_myeloid_TGFB1_by_sample.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_TGFB1_Robustness_v1.6"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV,  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_TREM2_TGFB1_Robustness_v1.6.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH TREM2 TGFB1 robustness audit v1.6.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(SENDER_TABLE))
stopifnot(file.exists(HSC_TABLE))
stopifnot(file.exists(TOTAL_TGFB1_TABLE))

# ------------------------------------------------------------
# 1. Load v1.5 derived tables
# ------------------------------------------------------------

sender <- read_tsv(
  SENDER_TABLE,
  show_col_types = FALSE
)

hsc <- read_tsv(
  HSC_TABLE,
  show_col_types = FALSE
)

total_tgfb1 <- read_tsv(
  TOTAL_TGFB1_TABLE,
  show_col_types = FALSE
)

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

required_sender <- c(
  "sample",
  "NAS",
  "sender_state",
  "state_fraction",
  "TGFB1_mean",
  "TGFB1_detection",
  "TGFB1_burden_expr",
  "TGFB1_burden_detection",
  "TGFB1_expression_share"
)

required_hsc <- c(
  "sample",
  "NAS",
  "TGF_direct_mean",
  "TGF_direct_high_fraction",
  "TGF_fibrogenic_mean",
  "ECM_high_fraction",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "ECM_remodeling_mean"
)

missing_sender <- setdiff(
  required_sender,
  colnames(sender)
)

missing_hsc <- setdiff(
  required_hsc,
  colnames(hsc)
)

if (length(missing_sender) > 0) {
  stop(
    "Missing sender columns: ",
    paste(missing_sender, collapse = ";")
  )
}

if (length(missing_hsc) > 0) {
  stop(
    "Missing HSC columns: ",
    paste(missing_hsc, collapse = ";")
  )
}

# ------------------------------------------------------------
# 2. Focus dataset
# ------------------------------------------------------------

trem2 <- sender %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  inner_join(
    hsc,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  arrange(sample)

if (nrow(trem2) != 8) {
  stop(
    "Expected 8 TREM2 sample rows; observed ",
    nrow(trem2)
  )
}

write_tsv(
  trem2,
  file.path(
    DIR_TABLE,
    "01_TREM2_TGFB1_HSC_focus_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

safe_spearman <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  if (
    sum(keep) < 4 ||
    length(unique(x[keep])) < 2 ||
    length(unique(y[keep])) < 2
  ) {
    return(
      tibble(
        n = sum(keep),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n = sum(keep),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

# Exact all-permutations generator for n <= 8
all_permutations <- function(x) {

  if (length(x) == 1) {
    return(
      matrix(
        x,
        nrow = 1
      )
    )
  }

  out <- vector(
    "list",
    length(x)
  )

  for (i in seq_along(x)) {

    rest <- x[-i]

    p_rest <- all_permutations(rest)

    out[[i]] <- cbind(
      x[i],
      p_rest
    )
  }

  do.call(
    rbind,
    out
  )
}

exact_permutation_spearman <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (n > 8) {
    stop(
      "Exact permutation implemented only for n <= 8."
    )
  }

  rho_obs <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  perm_idx <- all_permutations(
    seq_len(n)
  )

  rho_perm <- apply(
    perm_idx,
    1,
    function(idx) {
      suppressWarnings(
        cor(
          x,
          y[idx],
          method = "spearman"
        )
      )
    }
  )

  # exact two-sided permutation P
  p_exact <- mean(
    abs(rho_perm) >=
      abs(rho_obs) -
      1e-12
  )

  tibble(
    n = n,
    rho_observed = rho_obs,
    permutations = length(rho_perm),
    p_exact_two_sided = p_exact,
    perm_rho_mean = mean(rho_perm),
    perm_rho_sd = sd(rho_perm)
  )
}

partial_spearman_one_covariate <- function(
  x,
  y,
  z
) {

  keep <- complete.cases(
    x,
    y,
    z
  )

  x <- x[keep]
  y <- y[keep]
  z <- z[keep]

  n <- length(x)

  if (
    n < 5 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2 ||
    length(unique(z)) < 2
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        t_value = NA_real_,
        df = NA_integer_,
        p_value = NA_real_
      )
    )
  }

  r_xy <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  r_xz <- suppressWarnings(
    cor(
      x,
      z,
      method = "spearman"
    )
  )

  r_yz <- suppressWarnings(
    cor(
      y,
      z,
      method = "spearman"
    )
  )

  denom <- sqrt(
    (1 - r_xz^2) *
      (1 - r_yz^2)
  )

  if (
    !is.finite(denom) ||
    denom <= 0
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        t_value = NA_real_,
        df = n - 3L,
        p_value = NA_real_
      )
    )
  }

  r_partial <- (
    r_xy -
      r_xz *
      r_yz
  ) / denom

  r_partial <- max(
    min(
      r_partial,
      0.999999
    ),
    -0.999999
  )

  df <- n - 3L

  t_value <- r_partial *
    sqrt(
      df /
        (1 - r_partial^2)
    )

  p_value <- 2 *
    pt(
      -abs(t_value),
      df = df
    )

  tibble(
    n = n,
    partial_rho = r_partial,
    t_value = t_value,
    df = df,
    p_value = p_value
  )
}

# ------------------------------------------------------------
# 4. Primary association
# ------------------------------------------------------------

primary <- safe_spearman(
  trem2$TGFB1_burden_expr,
  trem2$TGF_direct_mean
) %>%
  mutate(
    sender_metric = "TREM2_TGFB1_burden_expr",
    receiver_metric = "TGF_direct_mean",
    .before = 1
  )

write_tsv(
  primary,
  file.path(
    DIR_TABLE,
    "02_PRIMARY_Spearman.tsv"
  )
)

# ------------------------------------------------------------
# 5. Exact permutation test
# ------------------------------------------------------------

perm_primary <- exact_permutation_spearman(
  trem2$TGFB1_burden_expr,
  trem2$TGF_direct_mean
)

write_tsv(
  perm_primary,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_exact_permutation.tsv"
  )
)

# ------------------------------------------------------------
# 6. Leave-one-sample-out robustness
# ------------------------------------------------------------

loo <- bind_rows(
  lapply(
    seq_len(nrow(trem2)),
    function(i) {

      sub <- trem2[-i, ]

      safe_spearman(
        sub$TGFB1_burden_expr,
        sub$TGF_direct_mean
      ) %>%
        mutate(
          omitted_sample =
            trem2$sample[i],
          omitted_NAS =
            trem2$NAS[i],
          .before = 1
        )
    }
  )
)

write_tsv(
  loo,
  file.path(
    DIR_TABLE,
    "04_leave_one_sample_out.tsv"
  )
)

loo_summary <- tibble(
  tests = nrow(loo),
  all_positive =
    all(
      loo$rho > 0,
      na.rm = TRUE
    ),
  min_rho =
    min(
      loo$rho,
      na.rm = TRUE
    ),
  median_rho =
    median(
      loo$rho,
      na.rm = TRUE
    ),
  max_rho =
    max(
      loo$rho,
      na.rm = TRUE
    ),
  p_lt_0_05_n =
    sum(
      loo$p_value < 0.05,
      na.rm = TRUE
    ),
  p_lt_0_10_n =
    sum(
      loo$p_value < 0.10,
      na.rm = TRUE
    )
)

write_tsv(
  loo_summary,
  file.path(
    DIR_TABLE,
    "05_leave_one_sample_out_summary.tsv"
  )
)

# ------------------------------------------------------------
# 7. Leave-one-NAS-group-out
# ------------------------------------------------------------

nas_groups <- sort(
  unique(
    trem2$NAS
  )
)

logo_nas <- bind_rows(
  lapply(
    nas_groups,
    function(nas_drop) {

      sub <- trem2 %>%
        filter(
          NAS != nas_drop
        )

      safe_spearman(
        sub$TGFB1_burden_expr,
        sub$TGF_direct_mean
      ) %>%
        mutate(
          omitted_NAS =
            nas_drop,
          samples_remaining =
            paste(
              sub$sample,
              collapse = ";"
            ),
          .before = 1
        )
    }
  )
)

write_tsv(
  logo_nas,
  file.path(
    DIR_TABLE,
    "06_leave_one_NAS_group_out.tsv"
  )
)

# Explicit NAS6 exclusion
nas6_excluded <- trem2 %>%
  filter(
    NAS != 6
  )

nas6_result <- safe_spearman(
  nas6_excluded$TGFB1_burden_expr,
  nas6_excluded$TGF_direct_mean
) %>%
  mutate(
    analysis = "NAS6_excluded",
    samples =
      paste(
        nas6_excluded$sample,
        collapse = ";"
      ),
    .before = 1
  )

write_tsv(
  nas6_result,
  file.path(
    DIR_TABLE,
    "07_NAS6_excluded_Spearman.tsv"
  )
)

# ------------------------------------------------------------
# 8. Partial Spearman controlling for NAS
# ------------------------------------------------------------

partial_nas <- partial_spearman_one_covariate(
  trem2$TGFB1_burden_expr,
  trem2$TGF_direct_mean,
  trem2$NAS
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "TGF_direct_mean",
    covariate =
      "NAS",
    .before = 1
  )

write_tsv(
  partial_nas,
  file.path(
    DIR_TABLE,
    "08_partial_Spearman_controlling_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 9. Burden decomposition
# ------------------------------------------------------------

decomp_metrics <- c(
  "state_fraction",
  "TGFB1_mean",
  "TGFB1_detection",
  "TGFB1_burden_expr",
  "TGFB1_burden_detection",
  "TGFB1_expression_share"
)

decomp <- bind_rows(
  lapply(
    decomp_metrics,
    function(metric) {

      safe_spearman(
        trem2[[metric]],
        trem2$TGF_direct_mean
      ) %>%
        mutate(
          sender_metric = metric,
          receiver_metric =
            "TGF_direct_mean",
          .before = 1
        )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR,
    p_value
  )

write_tsv(
  decomp,
  file.path(
    DIR_TABLE,
    "09_TREM2_burden_decomposition.tsv"
  )
)

# ------------------------------------------------------------
# 10. All sender-state comparison for TGF direct response
# ------------------------------------------------------------

all_states_tgf <- sender %>%
  inner_join(
    hsc %>%
      select(
        sample,
        NAS,
        TGF_direct_mean
      ),
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  group_by(
    sender_state
  ) %>%
  group_modify(
    ~ safe_spearman(
      .x$TGFB1_burden_expr,
      .x$TGF_direct_mean
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    desc(rho)
  )

write_tsv(
  all_states_tgf,
  file.path(
    DIR_TABLE,
    "10_all_sender_states_TGFB1_burden_vs_TGF_direct.tsv"
  )
)

# ------------------------------------------------------------
# 11. Leave-one-out winner robustness across sender states
# ------------------------------------------------------------

samples <- sort(
  unique(
    hsc$sample
  )
)

loo_winner <- bind_rows(
  lapply(
    samples,
    function(drop_sample) {

      sub_sender <- sender %>%
        filter(
          sample != drop_sample
        )

      sub_hsc <- hsc %>%
        filter(
          sample != drop_sample
        ) %>%
        select(
          sample,
          NAS,
          TGF_direct_mean
        )

      joined <- sub_sender %>%
        inner_join(
          sub_hsc,
          by = c(
            "sample",
            "NAS"
          )
        )

      state_res <- joined %>%
        group_by(
          sender_state
        ) %>%
        group_modify(
          ~ safe_spearman(
            .x$TGFB1_burden_expr,
            .x$TGF_direct_mean
          )
        ) %>%
        ungroup() %>%
        arrange(
          desc(rho)
        )

      winner <- state_res %>%
        slice(1)

      tibble(
        omitted_sample = drop_sample,
        winner_state =
          winner$sender_state,
        winner_rho =
          winner$rho,
        TREM2_rho =
          state_res %>%
            filter(
              sender_state ==
                TREM2_STATE
            ) %>%
            pull(rho) %>%
            {
              if (
                length(.) == 0
              ) {
                NA_real_
              } else {
                .
              }
            },
        TREM2_rank =
          match(
            TREM2_STATE,
            state_res$sender_state
          )
      )
    }
  )
)

write_tsv(
  loo_winner,
  file.path(
    DIR_TABLE,
    "11_leave_one_out_sender_winner.tsv"
  )
)

winner_summary <- tibble(
  LOO_tests = nrow(loo_winner),
  TREM2_winner_n =
    sum(
      loo_winner$winner_state ==
        TREM2_STATE
    ),
  TREM2_top2_n =
    sum(
      loo_winner$TREM2_rank <= 2,
      na.rm = TRUE
    ),
  median_TREM2_rank =
    median(
      loo_winner$TREM2_rank,
      na.rm = TRUE
    )
)

write_tsv(
  winner_summary,
  file.path(
    DIR_TABLE,
    "12_leave_one_out_sender_winner_summary.tsv"
  )
)

# ------------------------------------------------------------
# 12. Total-myeloid TGFB1 burden comparison
# ------------------------------------------------------------

total_hsc <- total_tgfb1 %>%
  inner_join(
    hsc,
    by = c(
      "sample",
      "NAS"
    )
  )

total_assoc <- bind_rows(
  safe_spearman(
    total_hsc$total_TGFB1_burden_expr,
    total_hsc$TGF_direct_mean
  ) %>%
    mutate(
      sender_metric =
        "total_myeloid_TGFB1_burden_expr",
      receiver_metric =
        "TGF_direct_mean",
      .before = 1
    ),

  safe_spearman(
    total_hsc$total_TGFB1_detection_fraction,
    total_hsc$TGF_direct_mean
  ) %>%
    mutate(
      sender_metric =
        "total_myeloid_TGFB1_detection_fraction",
      receiver_metric =
        "TGF_direct_mean",
      .before = 1
    )
)

write_tsv(
  total_assoc,
  file.path(
    DIR_TABLE,
    "13_total_myeloid_TGFB1_vs_TGF_direct.tsv"
  )
)

# ------------------------------------------------------------
# 13. Additional endpoints for focused TREM2 burden
# ------------------------------------------------------------

secondary_endpoints <- c(
  "TGF_direct_high_fraction",
  "TGF_fibrogenic_mean",
  "ECM_high_fraction",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "ECM_remodeling_mean",
  "NAS"
)

secondary <- bind_rows(
  lapply(
    secondary_endpoints,
    function(metric) {

      safe_spearman(
        trem2$TGFB1_burden_expr,
        trem2[[metric]]
      ) %>%
        mutate(
          sender_metric =
            "TREM2_TGFB1_burden_expr",
          receiver_metric =
            metric,
          .before = 1
        )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  )

write_tsv(
  secondary,
  file.path(
    DIR_TABLE,
    "14_TREM2_TGFB1_secondary_endpoints.tsv"
  )
)

# ------------------------------------------------------------
# 14. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  trem2,
  aes(
    x = TGFB1_burden_expr,
    y = TGF_direct_mean,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "TREM2-Mphi TGFB1 burden vs HSC TGF-beta response",
    x = "TREM2-Mphi TGFB1 burden",
    y = "HSC TGF-beta direct-response mean"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TREM2_TGFB1_burden_vs_TGF_direct.pdf"
  ),
  p1,
  width = 6.5,
  height = 5.5
)

p2 <- ggplot(
  loo,
  aes(
    x = omitted_sample,
    y = rho
  )
) +
  geom_col() +
  theme_classic(base_size = 14) +
  labs(
    title = "Leave-one-sample-out robustness",
    x = "Omitted sample",
    y = "Spearman rho"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_leave_one_sample_out_rho.pdf"
  ),
  p2,
  width = 7.5,
  height = 5
)

p3 <- ggplot(
  all_states_tgf,
  aes(
    x = reorder(
      sender_state,
      rho
    ),
    y = rho
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "TGFB1 burden vs HSC TGF-beta response by sender state",
    x = NULL,
    y = "Spearman rho"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_sender_state_comparison_TGF_direct.pdf"
  ),
  p3,
  width = 7.5,
  height = 5.5
)

# ------------------------------------------------------------
# 15. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "source_sender_table",
    "source_HSC_table",
    "source_total_TGFB1_table",
    "primary_sender_state",
    "primary_sender_metric",
    "primary_receiver_metric",
    "exact_permutation_n",
    "partial_correlation_covariate",
    "association_unit",
    "reclustering",
    "source_overwritten"
  ),
  value = c(
    "v1.6.0",
    SENDER_TABLE,
    HSC_TABLE,
    TOTAL_TGFB1_TABLE,
    TREM2_STATE,
    "TGFB1_burden_expr",
    "TGF_direct_mean",
    "40320_for_n8",
    "NAS",
    "biological_sample",
    "none",
    "FALSE"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "robustness_provenance_v1.6.tsv"
  )
)

# ------------------------------------------------------------
# 16. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "focus_samples",
    "expected_focus_samples",
    "focus_sample_count_matches",
    "all_4_NAS_groups_present",
    "exact_permutation_completed",
    "LOO_tests",
    "leave_NAS_group_tests",
    "TREM2_state_present",
    "total_myeloid_table_present",
    "source_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      nrow(trem2)
    ),
    "8",
    as.character(
      nrow(trem2) == 8
    ),
    as.character(
      length(
        unique(
          trem2$NAS
        )
      ) == 4
    ),
    as.character(
      perm_primary$permutations ==
        40320
    ),
    as.character(
      nrow(loo)
    ),
    as.character(
      nrow(logo_nas)
    ),
    as.character(
      TREM2_STATE %in%
        sender$sender_state
    ),
    as.character(
      file.exists(
        TOTAL_TGFB1_TABLE
      )
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.6.tsv"
  )
)

# ------------------------------------------------------------
# 17. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PRIMARY ASSOCIATION\n")
cat("============================================================\n")
print(primary)

cat("\n")
cat("============================================================\n")
cat("EXACT PERMUTATION\n")
cat("============================================================\n")
print(perm_primary)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-SAMPLE-OUT\n")
cat("============================================================\n")
print(loo)
print(loo_summary)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-NAS-GROUP-OUT\n")
cat("============================================================\n")
print(logo_nas)

cat("\n")
cat("============================================================\n")
cat("NAS6 EXCLUDED\n")
cat("============================================================\n")
print(nas6_result)

cat("\n")
cat("============================================================\n")
cat("PARTIAL SPEARMAN CONTROLLING NAS\n")
cat("============================================================\n")
print(partial_nas)

cat("\n")
cat("============================================================\n")
cat("TREM2 BURDEN DECOMPOSITION\n")
cat("============================================================\n")
print(decomp)

cat("\n")
cat("============================================================\n")
cat("ALL SENDER STATES\n")
cat("============================================================\n")
print(all_states_tgf)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT SENDER WINNER\n")
cat("============================================================\n")
print(loo_winner)
print(winner_summary)

cat("\n")
cat("============================================================\n")
cat("TOTAL MYELOID TGFB1\n")
cat("============================================================\n")
print(total_assoc)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
