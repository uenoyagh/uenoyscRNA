#!/usr/bin/env Rscript

# ============================================================
# Human MASH PDGFA -> validated HSC phenotype audit v1.3.2
#
# PURPOSE
#   Validate the cross-state myeloid PDGFA signal against the
#   previously validated Human HSC ECM/proliferation phenotype
#   framework (Human_MASH_HSC_ECM_Proliferation_v1.2).
#
# PRIMARY ENDPOINTS
#   - ECM_high_fraction
#   - ECM_remodeling_mean
#
# SECONDARY
#   - ECM_total_mean
#   - ECM_collagen_mean
#   - Proliferation_score_mean
#   - Cycling_strict_fraction
#   - Cycling_sensitive_fraction
#   - Cycling_marker_only_fraction
#
# PREDICTORS
#   - total_PDGFA
#   - total_PDGFA_PDGFRA_product
#   - Resident/Kupffer-like PDGFA burden
#   - TREM2-associated PDGFA burden
#   - total_PDGFB and PDGFB-PDGFRB controls
#
# ROBUSTNESS
#   - sample-level Spearman
#   - partial Spearman after rank(NAS) regression
#   - leave-one-sample-out
#   - sparse predictor/endpoint audit
#
# n=8 Human samples. Exploratory only.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

ROBUST_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGF_Robustness_Adjudication_v1.1.0"
)

CROSS_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2"
)

HSC_ECM_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2"
)

AXIS_TSV <- file.path(
  ROBUST_ROOT,
  "tables",
  "01_Human_PDGF_sample_axis_with_NAS.tsv"
)

LIGAND_TSV <- file.path(
  CROSS_ROOT,
  "tables",
  "01_Human_myeloid_state_PDGF_ligand_burden.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.2"
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
  "Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.2.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH PDGFA -> validated HSC phenotype v1.3.2\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")
cat("Patch note: v1.3.2 retains strict canonical HSC summary selection and fixes endpoint sparsity classification: signed continuous scores are judged by finite/unique values, while fraction endpoints retain positive-count zero-inflation audit.\n\n")

stopifnot(
  file.exists(AXIS_TSV),
  file.exists(LIGAND_TSV),
  dir.exists(file.path(HSC_ECM_ROOT, "tables"))
)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(tibble(n = sum(ok), rho = NA_real_, p_nominal = NA_real_))
  }

  tt <- suppressWarnings(
    cor.test(
      x[ok],
      y[ok],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n = sum(ok),
    rho = unname(tt$estimate),
    p_nominal = tt$p.value
  )
}

partial_rank <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)

  if (
    sum(ok) < 4 ||
    length(unique(x[ok])) < 2 ||
    length(unique(y[ok])) < 2 ||
    length(unique(z[ok])) < 2
  ) {
    return(NA_real_)
  }

  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  zr <- rank(z[ok], ties.method = "average")

  rx <- residuals(lm(xr ~ zr))
  ry <- residuals(lm(yr ~ zr))

  suppressWarnings(cor(rx, ry))
}

loo_stats <- function(x, y) {
  vals <- vapply(
    seq_along(x),
    function(i) {
      keep <- seq_along(x) != i
      xx <- x[keep]
      yy <- y[keep]
      ok <- is.finite(xx) & is.finite(yy)

      if (
        sum(ok) < 3 ||
        length(unique(xx[ok])) < 2 ||
        length(unique(yy[ok])) < 2
      ) {
        return(NA_real_)
      }

      suppressWarnings(
        cor(
          xx[ok],
          yy[ok],
          method = "spearman"
        )
      )
    },
    numeric(1)
  )

  raw <- safe_cor(x, y)$rho
  raw_sign <- sign(raw)

  tibble(
    loo_median = median(vals, na.rm = TRUE),
    loo_min = ifelse(all(is.na(vals)), NA_real_, min(vals, na.rm = TRUE)),
    loo_max = ifelse(all(is.na(vals)), NA_real_, max(vals, na.rm = TRUE)),
    loo_sign_consistency =
      ifelse(
        !is.finite(raw) || raw_sign == 0,
        NA_real_,
        mean(sign(vals) == raw_sign, na.rm = TRUE)
      )
  )
}

# ------------------------------------------------------------
# Find the validated HSC sample-summary TSV automatically
# ------------------------------------------------------------

tsvs <- list.files(
  file.path(HSC_ECM_ROOT, "tables"),
  pattern = "\\.tsv$",
  full.names = TRUE
)

required_cols <- c(
  "sample",
  "NAS",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "ECM_remodeling_mean",
  "ECM_high_fraction",
  "Proliferation_score_mean",
  "Cycling_strict_fraction",
  "Cycling_sensitive_fraction",
  "Cycling_marker_only_fraction"
)

# Strong signature of the canonical 8-row HSC sample-summary table.
sample_summary_signature <- c(
  "HSC_n",
  "ECM_high_n",
  "Cycling_strict_in_ECM_high_n",
  "Cycling_strict_fraction_in_ECM_high",
  "Cycling_sensitive_fraction_in_ECM_high"
)

expected_samples <- paste0("LT", 1:8)

candidate_info <- lapply(
  tsvs,
  function(f) {

    dat <- tryCatch(
      read_tsv(
        f,
        show_col_types = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(dat)) return(NULL)

    sample_values <-
      if ("sample" %in% colnames(dat)) {
        sort(
          unique(
            as.character(
              dat$sample
            )
          )
        )
      } else {
        character(0)
      }

    tibble(
      file = f,
      basename = basename(f),
      has_required =
        all(
          required_cols %in%
            colnames(dat)
        ),
      has_sample_summary_signature =
        all(
          sample_summary_signature %in%
            colnames(dat)
        ),
      rows = nrow(dat),
      unique_samples =
        length(
          sample_values
        ),
      exact_LT1_to_LT8 =
        identical(
          sample_values,
          expected_samples
        )
    )
  }
)

candidate_info <- bind_rows(
  candidate_info
)

write_tsv(
  candidate_info,
  file.path(
    DIR_PROV,
    "01_HSC_sample_summary_file_discovery.tsv"
  )
)

valid_files <- candidate_info %>%
  filter(
    has_required,
    has_sample_summary_signature,
    rows == 8,
    unique_samples == 8,
    exact_LT1_to_LT8
  ) %>%
  pull(file)

if (length(valid_files) != 1) {

  cat(
    "\nValidated HSC sample-summary discovery candidates:\n"
  )

  print(
    candidate_info %>%
      filter(
        has_required |
          has_sample_summary_signature
      )
  )

  stop(
    "Expected exactly one canonical 8-row HSC sample-summary TSV after strict signature audit; found ",
    length(valid_files),
    ". Inspect provenance/01_HSC_sample_summary_file_discovery.tsv"
  )
}

HSC_SAMPLE_TSV <- valid_files[[1]]

cat("Validated HSC sample-summary source:\n")
cat(HSC_SAMPLE_TSV, "\n\n")

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

axis <- read_tsv(
  AXIS_TSV,
  show_col_types = FALSE
) %>%
  arrange(sample)

ligand <- read_tsv(
  LIGAND_TSV,
  show_col_types = FALSE
)

hsc <- read_tsv(
  HSC_SAMPLE_TSV,
  show_col_types = FALSE
) %>%
  filter(
    grepl("^LT[1-8]$", sample)
  ) %>%
  arrange(sample)

if (nrow(hsc) != 8 || n_distinct(hsc$sample) != 8) {
  stop("Validated HSC sample summary does not contain exactly LT1-LT8.")
}

# Exact state-specific PDGFA burden
state_wide <- ligand %>%
  filter(
    gene == "PDGFA",
    sender_state %in%
      c(
        "Resident/Kupffer-like_Mphi",
        "Lipid-associated/TREM2_Mphi"
      )
  ) %>%
  select(
    sample,
    sender_state,
    ligand_burden
  ) %>%
  mutate(
    sender_state = case_when(
      sender_state == "Resident/Kupffer-like_Mphi" ~ "Resident_PDGFA",
      sender_state == "Lipid-associated/TREM2_Mphi" ~ "TREM2_PDGFA",
      TRUE ~ sender_state
    )
  ) %>%
  pivot_wider(
    names_from = sender_state,
    values_from = ligand_burden,
    values_fill = 0
  )

dat <- axis %>%
  select(
    sample,
    NAS,
    total_PDGFA,
    total_PDGFB,
    total_PDGFA_PDGFRA_product,
    total_PDGFB_PDGFRB_product,
    HSC_PDGFRA,
    HSC_PDGFRB
  ) %>%
  left_join(
    state_wide,
    by = "sample"
  ) %>%
  left_join(
    hsc %>%
      select(
        all_of(required_cols)
      ),
    by = "sample",
    suffix = c("_axis", "_validated")
  ) %>%
  mutate(
    NAS =
      ifelse(
        is.finite(NAS_axis),
        NAS_axis,
        NAS_validated
      )
  ) %>%
  select(
    sample,
    NAS,
    total_PDGFA,
    Resident_PDGFA,
    TREM2_PDGFA,
    total_PDGFA_PDGFRA_product,
    HSC_PDGFRA,
    total_PDGFB,
    total_PDGFB_PDGFRB_product,
    HSC_PDGFRB,
    all_of(
      setdiff(
        required_cols,
        c("sample", "NAS")
      )
    )
  )

write_tsv(
  dat,
  file.path(
    DIR_TABLE,
    "01_integrated_PDGF_validated_HSC_sample_matrix.tsv"
  )
)

# ------------------------------------------------------------
# Associations
# ------------------------------------------------------------

predictors <- c(
  "total_PDGFA",
  "Resident_PDGFA",
  "TREM2_PDGFA",
  "total_PDGFA_PDGFRA_product",
  "HSC_PDGFRA",
  "total_PDGFB",
  "total_PDGFB_PDGFRB_product",
  "HSC_PDGFRB"
)

endpoints <- c(
  "ECM_high_fraction",
  "ECM_remodeling_mean",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "Proliferation_score_mean",
  "Cycling_strict_fraction",
  "Cycling_sensitive_fraction",
  "Cycling_marker_only_fraction"
)

rows <- list()
k <- 0L

for (pred in predictors) {
  for (endp in endpoints) {

    x <- dat[[pred]]
    y <- dat[[endp]]

    raw <- safe_cor(x, y)
    pr <- partial_rank(x, y, dat$NAS)
    loo <- loo_stats(x, y)

    nonzero_pred <- sum(is.finite(x) & x > 0)

    # Endpoint sparsity must depend on endpoint type.
    # Continuous signed scores (ECM/proliferation means) are informative
    # even when <= 0, whereas fraction endpoints can be truly zero-inflated.
    endpoint_finite_n <- sum(is.finite(y))
    endpoint_unique_n <- length(unique(y[is.finite(y)]))
    endpoint_positive_n <- sum(is.finite(y) & y > 0)

    endpoint_is_fraction <- grepl(
      "fraction$",
      endp
    )

    endpoint_informative_n <-
      if (endpoint_is_fraction) {
        endpoint_positive_n
      } else {
        endpoint_finite_n
      }

    k <- k + 1L

    rows[[k]] <- tibble(
      predictor = pred,
      endpoint = endp,
      n = raw$n,
      predictor_nonzero_n = nonzero_pred,
      endpoint_finite_n = endpoint_finite_n,
      endpoint_unique_n = endpoint_unique_n,
      endpoint_positive_n = endpoint_positive_n,
      endpoint_informative_n = endpoint_informative_n,
      rho_raw = raw$rho,
      p_nominal = raw$p_nominal,
      rho_partial_NAS = pr,
      loo_median = loo$loo_median,
      loo_min = loo$loo_min,
      loo_max = loo$loo_max,
      loo_sign_consistency = loo$loo_sign_consistency
    )
  }
}

res <- bind_rows(rows) %>%
  mutate(
    p_BH_all =
      p.adjust(
        p_nominal,
        method = "BH"
      ),
    predictor_sparse =
      predictor_nonzero_n < 5,
    endpoint_sparse =
      endpoint_informative_n < 5 |
      endpoint_unique_n < 3,
    partial_same_direction =
      sign(rho_raw) ==
      sign(rho_partial_NAS),
    robustness_class =
      case_when(
        predictor_sparse &
          endpoint_sparse ~
          "sparse_predictor_and_endpoint",

        predictor_sparse ~
          "sparse_predictor",

        endpoint_sparse ~
          "sparse_endpoint",

        abs(rho_raw) >= 0.7 &
          abs(rho_partial_NAS) >= 0.5 &
          partial_same_direction &
          loo_sign_consistency >= 0.875 ~
          "robust_beyond_NAS",

        abs(rho_raw) >= 0.6 &
          abs(rho_partial_NAS) >= 0.4 &
          partial_same_direction &
          loo_sign_consistency >= 0.875 ~
          "moderately_robust",

        TRUE ~
          "weak_or_unresolved"
      )
  ) %>%
  arrange(
    desc(abs(rho_raw))
  )

write_tsv(
  res,
  file.path(
    DIR_TABLE,
    "02_ALL_PDGF_validated_HSC_associations.tsv"
  )
)

primary <- res %>%
  filter(
    predictor %in%
      c(
        "total_PDGFA",
        "Resident_PDGFA",
        "TREM2_PDGFA",
        "total_PDGFA_PDGFRA_product",
        "HSC_PDGFRA"
      ),
    endpoint %in%
      c(
        "ECM_high_fraction",
        "ECM_remodeling_mean",
        "Proliferation_score_mean",
        "Cycling_strict_fraction",
        "Cycling_sensitive_fraction"
      )
  ) %>%
  arrange(
    endpoint,
    desc(abs(rho_raw))
  )

write_tsv(
  primary,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_PDGFA_validated_HSC_phenotype.tsv"
  )
)

# ------------------------------------------------------------
# Adjudication summary
# ------------------------------------------------------------

adjudication <- primary %>%
  mutate(
    interpretation =
      case_when(
        predictor == "TREM2_PDGFA" &
          predictor_sparse ~
          "sparse_TREM2_signal_sensitivity_only",

        endpoint %in%
          c(
            "Proliferation_score_mean",
            "Cycling_strict_fraction",
            "Cycling_sensitive_fraction"
          ) &
          endpoint_sparse ~
          "proliferation_endpoint_sparse",

        robustness_class == "robust_beyond_NAS" ~
          "primary_support",

        robustness_class == "moderately_robust" ~
          "secondary_support",

        TRUE ~
          "not_robust"
      )
  )

write_tsv(
  adjudication,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_interpretation.tsv"
  )
)

# ------------------------------------------------------------
# Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "PDGF_axis_source",
    "PDGF_state_source",
    "validated_HSC_source",
    "samples_n",
    "samples",
    "primary_endpoints",
    "endpoint_sparsity_rule",
    "inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    AXIS_TSV,
    LIGAND_TSV,
    HSC_SAMPLE_TSV,
    as.character(nrow(dat)),
    paste(dat$sample, collapse = ";"),
    "ECM_high_fraction;ECM_remodeling_mean",
    "fraction endpoints: positive_n>=5 and unique_n>=3; continuous signed scores: finite_n>=5 and unique_n>=3",
    "exploratory_sample_level_n8",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.3.2.tsv"
  )
)

saveRDS(
  list(
    sample_matrix = dat,
    associations = res,
    primary = primary,
    adjudication = adjudication
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.2_results.rds"
  )
)

cat("\n============================================================\n")
cat("PRIMARY PDGFA -> VALIDATED HSC PHENOTYPE\n")
cat("============================================================\n")
print(primary)

cat("\n============================================================\n")
cat("PRIMARY INTERPRETATION\n")
cat("============================================================\n")
print(adjudication)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
