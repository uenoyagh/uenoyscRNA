#!/usr/bin/env Rscript

# ============================================================
# Human MASH PDGF robustness adjudication v1.1.0
#
# PURPOSE
#   Adjudicate the Human PDGFA/PDGFB -> HSC signal after the
#   corrected cross-species v1.0.2 analysis.
#
# PRIMARY QUESTIONS
#   1. Is total-myeloid PDGFA -> HSC PDGFRA association robust?
#   2. Is the apparent TREM2-PDGFA signal driven by sparsity /
#      shared zero patterns?
#   3. Does the PDGFA association remain after controlling for NAS?
#   4. Is the signal robust to leave-one-sample-out analysis?
#   5. Which myeloid state is the dominant PDGFA / PDGFB source?
#   6. Does HSC PDGFRA expression alone explain the association?
#
# STATISTICS
#   n = 8 Human samples (LT1-LT8).
#   Exploratory sample-level analysis.
#   No cell-level pseudoreplication.
#   Partial Spearman = Pearson correlation of rank residuals after
#   rank(NAS) regression.
#
# INTERPRETATION
#   - Sparse predictors (<5/8 nonzero samples) are not treated as
#     continuous dose-response evidence.
#   - LOO sign consistency is descriptive robustness only.
#   - No causal claim from correlation.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

SRC <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2"
)

SAMPLE_AXIS_TSV <- file.path(
  SRC,
  "tables",
  "06_Human_sample_PDGF_axis_summary.tsv"
)

LIGAND_STATE_TSV <- file.path(
  SRC,
  "tables",
  "01_Human_myeloid_state_PDGF_ligand_burden.tsv"
)

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_Receiver_Primary_v1.0",
  "objects",
  "Human_MASH_HSC_Receiver_Primary_v1.0.rds"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGF_Robustness_Adjudication_v1.1.0"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_LOG <- file.path(OUT, "logs")
DIR_OBJ <- file.path(OUT, "objects")
DIR_PROV <- file.path(OUT, "provenance")

for (d in c(DIR_TABLE, DIR_LOG, DIR_OBJ, DIR_PROV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(
  DIR_LOG,
  "Human_MASH_PDGF_Robustness_Adjudication_v1.1.0.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH PDGF robustness adjudication v1.1.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(SAMPLE_AXIS_TSV),
  file.exists(LIGAND_STATE_TSV),
  file.exists(HSC_RDS)
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  }

  tt <- suppressWarnings(
    cor.test(
      x[ok],
      y[ok],
      method = "spearman",
      exact = FALSE
    )
  )

  c(
    n = sum(ok),
    rho = unname(tt$estimate),
    p = tt$p.value
  )
}

partial_spearman_rank_residual <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)

  if (sum(ok) < 4) {
    return(NA_real_)
  }

  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  zr <- rank(z[ok], ties.method = "average")

  rx <- residuals(lm(xr ~ zr))
  ry <- residuals(lm(yr ~ zr))

  suppressWarnings(
    cor(
      rx,
      ry,
      method = "pearson"
    )
  )
}

loo_spearman <- function(x, y, sample_names) {
  out <- lapply(
    seq_along(x),
    function(i) {
      keep <- seq_along(x) != i
      rr <- suppressWarnings(
        cor(
          x[keep],
          y[keep],
          method = "spearman",
          use = "complete.obs"
        )
      )

      tibble(
        left_out = sample_names[i],
        rho = rr
      )
    }
  )

  bind_rows(out)
}

direction_sign <- function(x) {
  ifelse(
    is.na(x),
    NA_integer_,
    ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
  )
}

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

axis <- read_tsv(
  SAMPLE_AXIS_TSV,
  show_col_types = FALSE
)

ligand_state <- read_tsv(
  LIGAND_STATE_TSV,
  show_col_types = FALSE
)

hsc <- readRDS(HSC_RDS)
hmd <- hsc[[]] %>%
  rownames_to_column("cell")

if (!all(c("sample", "NAS_score") %in% colnames(hmd))) {
  stop("Human HSC RDS lacks sample and/or NAS_score metadata.")
}

nas_by_sample <- hmd %>%
  transmute(
    sample = as.character(sample),
    NAS = as.numeric(NAS_score)
  ) %>%
  filter(
    grepl("^LT[1-8]$", sample)
  ) %>%
  group_by(sample) %>%
  summarise(
    NAS = median(NAS, na.rm = TRUE),
    NAS_unique_n = n_distinct(NAS[is.finite(NAS)]),
    .groups = "drop"
  )

if (nrow(nas_by_sample) != 8) {
  stop("Expected 8 LT samples in NAS table.")
}

axis <- axis %>%
  left_join(
    nas_by_sample %>%
      select(sample, NAS),
    by = "sample"
  ) %>%
  arrange(sample)

if (any(!is.finite(axis$NAS))) {
  stop("NAS missing for one or more Human samples.")
}

write_tsv(
  axis,
  file.path(
    DIR_TABLE,
    "01_Human_PDGF_sample_axis_with_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 1. Myeloid-state source decomposition
# ------------------------------------------------------------

source_summary <- ligand_state %>%
  group_by(
    gene,
    sender_state
  ) %>%
  summarise(
    samples_n = n_distinct(sample),
    samples_nonzero_burden =
      sum(
        ligand_burden > 0,
        na.rm = TRUE
      ),
    mean_ligand_burden =
      mean(
        ligand_burden,
        na.rm = TRUE
      ),
    median_ligand_burden =
      median(
        ligand_burden,
        na.rm = TRUE
      ),
    sum_ligand_burden =
      sum(
        ligand_burden,
        na.rm = TRUE
      ),
    mean_detection_fraction =
      mean(
        detection_fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    share_of_sum_burden =
      sum_ligand_burden /
      sum(
        sum_ligand_burden,
        na.rm = TRUE
      )
  ) %>%
  ungroup() %>%
  arrange(
    gene,
    desc(share_of_sum_burden)
  )

write_tsv(
  source_summary,
  file.path(
    DIR_TABLE,
    "02_PDGF_myeloid_state_source_decomposition.tsv"
  )
)

# ------------------------------------------------------------
# 2. Predictor / endpoint definitions
# ------------------------------------------------------------

predictors <- c(
  "total_PDGFA",
  "total_PDGFB",
  "TREM2_PDGFA",
  "TREM2_PDGFB",
  "HSC_PDGFRA",
  "HSC_PDGFRB",
  "total_PDGFA_PDGFRA_product",
  "total_PDGFB_PDGFRB_product",
  "TREM2_PDGFA_PDGFRA_product",
  "TREM2_PDGFB_PDGFRB_product"
)

endpoints <- c(
  "aHSC_fraction",
  "aHSC_core_collagen",
  "aHSC_ECM_remodeling",
  "aHSC_cycling",
  "wholeHSC_core_collagen",
  "wholeHSC_ECM_remodeling",
  "wholeHSC_cycling"
)

stopifnot(
  all(predictors %in% colnames(axis)),
  all(endpoints %in% colnames(axis))
)

# ------------------------------------------------------------
# 3. NAS associations
# ------------------------------------------------------------

predictor_nas <- bind_rows(
  lapply(
    predictors,
    function(pred) {
      ss <- safe_spearman(
        axis[[pred]],
        axis$NAS
      )

      tibble(
        variable = pred,
        variable_type = "predictor",
        n = ss[["n"]],
        nonzero_n =
          sum(
            is.finite(axis[[pred]]) &
              axis[[pred]] > 0
          ),
        rho_vs_NAS = ss[["rho"]],
        p_nominal = ss[["p"]]
      )
    }
  )
)

endpoint_nas <- bind_rows(
  lapply(
    endpoints,
    function(endp) {
      ss <- safe_spearman(
        axis[[endp]],
        axis$NAS
      )

      tibble(
        variable = endp,
        variable_type = "endpoint",
        n = ss[["n"]],
        nonzero_n =
          sum(
            is.finite(axis[[endp]]) &
              axis[[endp]] > 0
          ),
        rho_vs_NAS = ss[["rho"]],
        p_nominal = ss[["p"]]
      )
    }
  )
)

nas_audit <- bind_rows(
  predictor_nas,
  endpoint_nas
) %>%
  mutate(
    p_BH =
      p.adjust(
        p_nominal,
        method = "BH"
      )
  )

write_tsv(
  nas_audit,
  file.path(
    DIR_TABLE,
    "03_PDGF_and_HSC_variables_vs_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 4. Primary robustness matrix
# ------------------------------------------------------------

robustness <- list()
loo_all <- list()
k <- 0L
kk <- 0L

for (pred in predictors) {
  for (endp in endpoints) {

    x <- axis[[pred]]
    y <- axis[[endp]]

    raw <- safe_spearman(x, y)

    partial_rho <- partial_spearman_rank_residual(
      x,
      y,
      axis$NAS
    )

    loo <- loo_spearman(
      x,
      y,
      axis$sample
    )

    raw_sign <- direction_sign(
      raw[["rho"]]
    )

    loo_sign_consistency <-
      if (is.na(raw_sign) || raw_sign == 0) {
        NA_real_
      } else {
        mean(
          direction_sign(
            loo$rho
          ) == raw_sign,
          na.rm = TRUE
        )
      }

    predictor_nonzero_n <-
      sum(
        is.finite(x) &
          x > 0
      )

    endpoint_nonzero_n <-
      sum(
        is.finite(y) &
          y > 0
      )

    both_positive_n <-
      sum(
        is.finite(x) &
          is.finite(y) &
          x > 0 &
          y > 0
      )

    both_zero_n <-
      sum(
        is.finite(x) &
          is.finite(y) &
          x == 0 &
          y == 0
      )

    k <- k + 1L

    robustness[[k]] <- tibble(
      predictor = pred,
      endpoint = endp,
      n = raw[["n"]],
      predictor_nonzero_n = predictor_nonzero_n,
      endpoint_nonzero_n = endpoint_nonzero_n,
      both_positive_n = both_positive_n,
      both_zero_n = both_zero_n,
      rho_raw = raw[["rho"]],
      p_nominal = raw[["p"]],
      rho_partial_NAS = partial_rho,
      loo_rho_median =
        median(
          loo$rho,
          na.rm = TRUE
        ),
      loo_rho_min =
        min(
          loo$rho,
          na.rm = TRUE
        ),
      loo_rho_max =
        max(
          loo$rho,
          na.rm = TRUE
        ),
      loo_sign_consistency =
        loo_sign_consistency
    )

    kk <- kk + 1L
    loo_all[[kk]] <- loo %>%
      mutate(
        predictor = pred,
        endpoint = endp
      ) %>%
      select(
        predictor,
        endpoint,
        left_out,
        rho
      )
  }
}

robustness <- bind_rows(
  robustness
) %>%
  mutate(
    p_BH_all70 =
      p.adjust(
        p_nominal,
        method = "BH"
      ),
    predictor_sparse =
      predictor_nonzero_n < 5,
    endpoint_sparse =
      endpoint_nonzero_n < 5,
    partial_same_direction =
      direction_sign(rho_raw) ==
      direction_sign(rho_partial_NAS),
    interpretation_class =
      case_when(
        predictor_sparse &
          endpoint_sparse ~
          "sparse_predictor_and_endpoint",

        predictor_sparse ~
          "sparse_predictor",

        abs(rho_raw) >= 0.7 &
          abs(rho_partial_NAS) >= 0.5 &
          partial_same_direction &
          loo_sign_consistency >= 0.875 ~
          "robust_association_beyond_NAS",

        abs(rho_raw) >= 0.7 &
          (
            is.na(rho_partial_NAS) |
              abs(rho_partial_NAS) < 0.4
          ) ~
          "possible_NAS_confounding",

        abs(rho_raw) >= 0.6 &
          loo_sign_consistency >= 0.75 ~
          "moderately_robust",

        TRUE ~
          "weak_or_unresolved"
      )
  ) %>%
  arrange(
    desc(
      abs(rho_raw)
    )
  )

write_tsv(
  robustness,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_PDGF_HSC_robustness_adjudication.tsv"
  )
)

write_tsv(
  bind_rows(loo_all),
  file.path(
    DIR_TABLE,
    "05_leave_one_out_correlations.tsv"
  )
)

# ------------------------------------------------------------
# 5. Focused sparse-pattern audit
# ------------------------------------------------------------

focus_pairs <- tribble(
  ~predictor, ~endpoint,
  "TREM2_PDGFA", "aHSC_fraction",
  "TREM2_PDGFA", "aHSC_cycling",
  "TREM2_PDGFA_PDGFRA_product", "aHSC_fraction",
  "TREM2_PDGFA_PDGFRA_product", "aHSC_cycling",
  "total_PDGFA", "aHSC_ECM_remodeling",
  "total_PDGFA_PDGFRA_product", "aHSC_ECM_remodeling",
  "total_PDGFA_PDGFRA_product", "wholeHSC_ECM_remodeling",
  "TREM2_PDGFB_PDGFRB_product", "wholeHSC_core_collagen",
  "TREM2_PDGFB_PDGFRB_product", "wholeHSC_ECM_remodeling"
)

zero_pattern <- bind_rows(
  lapply(
    seq_len(nrow(focus_pairs)),
    function(i) {

      pred <- focus_pairs$predictor[i]
      endp <- focus_pairs$endpoint[i]

      tibble(
        sample = axis$sample,
        NAS = axis$NAS,
        predictor = pred,
        endpoint = endp,
        predictor_value = axis[[pred]],
        endpoint_value = axis[[endp]],
        predictor_positive =
          axis[[pred]] > 0,
        endpoint_positive =
          axis[[endp]] > 0
      )
    }
  )
)

write_tsv(
  zero_pattern,
  file.path(
    DIR_TABLE,
    "06_focus_pair_zero_pattern_by_sample.tsv"
  )
)

zero_summary <- zero_pattern %>%
  group_by(
    predictor,
    endpoint
  ) %>%
  summarise(
    predictor_positive_n =
      sum(
        predictor_positive,
        na.rm = TRUE
      ),
    endpoint_positive_n =
      sum(
        endpoint_positive,
        na.rm = TRUE
      ),
    both_positive_n =
      sum(
        predictor_positive &
          endpoint_positive,
        na.rm = TRUE
      ),
    both_zero_n =
      sum(
        !predictor_positive &
          !endpoint_positive,
        na.rm = TRUE
      ),
    discordant_binary_n =
      sum(
        xor(
          predictor_positive,
          endpoint_positive
        ),
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_tsv(
  zero_summary,
  file.path(
    DIR_TABLE,
    "07_focus_pair_zero_pattern_summary.tsv"
  )
)

# ------------------------------------------------------------
# 6. PDGFRA-alone audit
# ------------------------------------------------------------

receptor_only <- robustness %>%
  filter(
    predictor %in%
      c(
        "HSC_PDGFRA",
        "HSC_PDGFRB"
      )
  )

write_tsv(
  receptor_only,
  file.path(
    DIR_TABLE,
    "08_HSC_receptor_alone_vs_HSC_endpoints.tsv"
  )
)

# ------------------------------------------------------------
# 7. Primary concise interpretation table
# ------------------------------------------------------------

primary_pairs <- robustness %>%
  filter(
    predictor %in%
      c(
        "total_PDGFA",
        "total_PDGFA_PDGFRA_product",
        "TREM2_PDGFA",
        "TREM2_PDGFA_PDGFRA_product",
        "total_PDGFB",
        "total_PDGFB_PDGFRB_product",
        "TREM2_PDGFB",
        "TREM2_PDGFB_PDGFRB_product"
      ),
    endpoint %in%
      c(
        "aHSC_fraction",
        "aHSC_ECM_remodeling",
        "aHSC_cycling",
        "wholeHSC_core_collagen",
        "wholeHSC_ECM_remodeling",
        "wholeHSC_cycling"
      )
  ) %>%
  mutate(
    priority_note =
      case_when(
        predictor_sparse ~
          "Do not interpret as continuous dose-response; sparse predictor.",
        endpoint_sparse ~
          "Endpoint is sparse; treat correlation as sensitivity only.",
        interpretation_class ==
          "robust_association_beyond_NAS" ~
          "Primary robustness candidate.",
        TRUE ~
          "Exploratory."
      )
  )

write_tsv(
  primary_pairs,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_PDGF_interpretation_table.tsv"
  )
)

# ------------------------------------------------------------
# 8. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "source_cross_species_output",
    "HSC_source",
    "samples_n",
    "samples",
    "NAS_source",
    "partial_correlation_method",
    "LOO_method",
    "sparse_predictor_definition",
    "inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    SRC,
    HSC_RDS,
    as.character(nrow(axis)),
    paste(axis$sample, collapse = ";"),
    "Human HSC receiver metadata NAS_score; sample median",
    "Spearman rank residual correlation after rank(NAS) regression",
    "leave-one-biological-sample-out Spearman",
    "nonzero in <5/8 samples",
    "exploratory_sample_level_n8",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.1.0.tsv"
  )
)

saveRDS(
  list(
    axis = axis,
    source_summary = source_summary,
    nas_audit = nas_audit,
    robustness = robustness,
    zero_summary = zero_summary,
    receptor_only = receptor_only,
    primary_pairs = primary_pairs
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_PDGF_Robustness_Adjudication_v1.1.0_results.rds"
  )
)

cat("\n============================================================\n")
cat("PDGF MYELOID SOURCE DECOMPOSITION\n")
cat("============================================================\n")
print(source_summary)

cat("\n============================================================\n")
cat("VARIABLES vs NAS\n")
cat("============================================================\n")
print(nas_audit)

cat("\n============================================================\n")
cat("TOP ROBUSTNESS RESULTS\n")
cat("============================================================\n")
print(
  robustness %>%
    select(
      predictor,
      endpoint,
      predictor_nonzero_n,
      endpoint_nonzero_n,
      rho_raw,
      p_nominal,
      rho_partial_NAS,
      loo_rho_median,
      loo_rho_min,
      loo_rho_max,
      loo_sign_consistency,
      interpretation_class
    ) %>%
    slice_head(n = 30)
)

cat("\n============================================================\n")
cat("ZERO-PATTERN AUDIT\n")
cat("============================================================\n")
print(zero_summary)

cat("\n============================================================\n")
cat("HSC RECEPTOR-ONLY AUDIT\n")
cat("============================================================\n")
print(receptor_only)

cat("\n============================================================\n")
cat("PRIMARY INTERPRETATION TABLE\n")
cat("============================================================\n")
print(primary_pairs)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
