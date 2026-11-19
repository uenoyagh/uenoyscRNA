#!/usr/bin/env Rscript

# ============================================================
# Human MASH PDGFA sender-state attribution v1.2.0
#
# PURPOSE
#   Resolve which Human myeloid state(s) drive the robust
#   total-myeloid PDGFA association with HSC remodeling.
#
# INPUTS
#   - v1.0.2 state-level ligand burden table
#   - v1.1.0 sample axis with NAS
#
# COMPONENTS
#   For each myeloid state and PDGFA/PDGFB:
#     1) ligand_burden = state_fraction x mean_CP10k
#     2) state_fraction
#     3) mean_CP10k
#     4) detection_fraction
#
# ENDPOINTS
#   Primary:
#     aHSC_ECM_remodeling
#     wholeHSC_ECM_remodeling
#   Secondary:
#     aHSC_core_collagen
#     wholeHSC_core_collagen
#     aHSC_fraction
#   Sensitivity:
#     aHSC_cycling
#     wholeHSC_cycling
#
# ROBUSTNESS
#   - raw sample-level Spearman
#   - partial Spearman after rank(NAS) regression
#   - leave-one-sample-out Spearman
#   - nonzero sample audit
#
# INFERENCE
#   n=8; exploratory sample-level only.
#   No cell-level pseudoreplication.
#   No source-RDS modification or reclustering.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

CROSS_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2"
)

ROBUST_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGF_Robustness_Adjudication_v1.1.0"
)

LIGAND_STATE_TSV <- file.path(
  CROSS_ROOT,
  "tables",
  "01_Human_myeloid_state_PDGF_ligand_burden.tsv"
)

AXIS_NAS_TSV <- file.path(
  ROBUST_ROOT,
  "tables",
  "01_Human_PDGF_sample_axis_with_NAS.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGFA_SenderState_Attribution_v1.2.0"
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
  "Human_MASH_PDGFA_SenderState_Attribution_v1.2.0.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH PDGFA sender-state attribution v1.2.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(LIGAND_STATE_TSV),
  file.exists(AXIS_NAS_TSV)
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) {
    return(
      tibble(
        n = sum(ok),
        rho = NA_real_,
        p_nominal = NA_real_
      )
    )
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

loo_summary <- function(x, y, sample_names) {

  rhos <- vapply(
    seq_along(x),
    function(i) {

      keep <- seq_along(x) != i

      suppressWarnings(
        cor(
          x[keep],
          y[keep],
          method = "spearman",
          use = "complete.obs"
        )
      )
    },
    numeric(1)
  )

  raw_rho <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman",
      use = "complete.obs"
    )
  )

  raw_sign <- sign(raw_rho)

  sign_consistency <-
    if (
      !is.finite(raw_rho) ||
      raw_sign == 0
    ) {
      NA_real_
    } else {
      mean(
        sign(rhos) == raw_sign,
        na.rm = TRUE
      )
    }

  tibble(
    loo_rho_median = median(rhos, na.rm = TRUE),
    loo_rho_min = min(rhos, na.rm = TRUE),
    loo_rho_max = max(rhos, na.rm = TRUE),
    loo_sign_consistency = sign_consistency
  )
}

component_driver_class <- function(
  burden_row,
  fraction_row,
  expression_row
) {

  robust_component <- function(row) {
    if (nrow(row) == 0) return(FALSE)

    isTRUE(
      row$nonzero_n >= 5 &&
      abs(row$rho_raw) >= 0.6 &&
      abs(row$rho_partial_NAS) >= 0.4 &&
      row$partial_same_direction &&
      row$loo_sign_consistency >= 0.875
    )
  }

  fr <- robust_component(fraction_row)
  ex <- robust_component(expression_row)

  case_when(
    fr & ex ~ "mixed_composition_and_expression",
    fr & !ex ~ "composition_dominant",
    !fr & ex ~ "within_state_expression_dominant",
    TRUE ~ "unresolved_component_driver"
  )
}

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

ligand_state <- read_tsv(
  LIGAND_STATE_TSV,
  show_col_types = FALSE
)

axis <- read_tsv(
  AXIS_NAS_TSV,
  show_col_types = FALSE
) %>%
  arrange(sample)

samples <- axis$sample

if (nrow(axis) != 8) {
  stop("Expected 8 Human samples.")
}

needed_ligand_cols <- c(
  "sample",
  "sender_state",
  "gene",
  "state_fraction",
  "mean_CP10k",
  "detection_fraction",
  "ligand_burden"
)

if (!all(needed_ligand_cols %in% colnames(ligand_state))) {
  stop("Ligand-state input lacks required columns.")
}

# ------------------------------------------------------------
# 1. Sample-level source share
# ------------------------------------------------------------

source_share <- ligand_state %>%
  group_by(
    sample,
    gene
  ) %>%
  mutate(
    total_gene_burden =
      sum(
        ligand_burden,
        na.rm = TRUE
      ),
    source_share =
      ifelse(
        total_gene_burden > 0,
        ligand_burden / total_gene_burden,
        NA_real_
      )
  ) %>%
  ungroup() %>%
  arrange(
    gene,
    sample,
    desc(source_share)
  )

write_tsv(
  source_share,
  file.path(
    DIR_TABLE,
    "01_PDGF_source_share_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 2. Build state-component x endpoint association matrix
# ------------------------------------------------------------

components <- c(
  "ligand_burden",
  "state_fraction",
  "mean_CP10k",
  "detection_fraction"
)

endpoints <- c(
  "aHSC_ECM_remodeling",
  "wholeHSC_ECM_remodeling",
  "aHSC_core_collagen",
  "wholeHSC_core_collagen",
  "aHSC_fraction",
  "aHSC_cycling",
  "wholeHSC_cycling"
)

association_rows <- list()
k <- 0L

for (gene_name in c("PDGFA", "PDGFB")) {

  states <- sort(
    unique(
      ligand_state$sender_state[
        ligand_state$gene ==
          gene_name
      ]
    )
  )

  for (state_name in states) {

    dat_state <- ligand_state %>%
      filter(
        gene == gene_name,
        sender_state == state_name
      ) %>%
      select(
        sample,
        all_of(components)
      ) %>%
      right_join(
        axis %>%
          select(
            sample,
            NAS,
            all_of(endpoints)
          ),
        by = "sample"
      ) %>%
      arrange(sample)

    for (component_name in components) {

      x <- dat_state[[component_name]]

      for (endpoint_name in endpoints) {

        y <- dat_state[[endpoint_name]]

        ss <- safe_spearman(
          x,
          y
        )

        partial_rho <- partial_spearman_rank_residual(
          x,
          y,
          dat_state$NAS
        )

        loo <- loo_summary(
          x,
          y,
          dat_state$sample
        )

        k <- k + 1L

        association_rows[[k]] <- tibble(
          gene = gene_name,
          sender_state = state_name,
          component = component_name,
          endpoint = endpoint_name,
          n = ss$n,
          nonzero_n =
            sum(
              is.finite(x) &
              x > 0
            ),
          rho_raw = ss$rho,
          p_nominal = ss$p_nominal,
          rho_partial_NAS = partial_rho,
          loo_rho_median = loo$loo_rho_median,
          loo_rho_min = loo$loo_rho_min,
          loo_rho_max = loo$loo_rho_max,
          loo_sign_consistency = loo$loo_sign_consistency
        )
      }
    }
  }
}

assoc <- bind_rows(
  association_rows
) %>%
  mutate(
    p_BH_all =
      p.adjust(
        p_nominal,
        method = "BH"
      ),
    sparse_component =
      nonzero_n < 5,
    partial_same_direction =
      sign(rho_raw) ==
      sign(rho_partial_NAS),
    robustness_class =
      case_when(
        sparse_component ~
          "sparse_component",

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
    gene,
    endpoint,
    desc(abs(rho_raw))
  )

write_tsv(
  assoc,
  file.path(
    DIR_TABLE,
    "02_state_component_HSC_associations.tsv"
  )
)

# ------------------------------------------------------------
# 3. Primary PDGFA state attribution
# ------------------------------------------------------------

primary_pdgfa <- assoc %>%
  filter(
    gene == "PDGFA",
    component == "ligand_burden",
    endpoint %in%
      c(
        "aHSC_ECM_remodeling",
        "wholeHSC_ECM_remodeling"
      )
  ) %>%
  arrange(
    endpoint,
    desc(abs(rho_raw))
  )

write_tsv(
  primary_pdgfa,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_PDGFA_sender_state_attribution.tsv"
  )
)

# ------------------------------------------------------------
# 4. Burden-component driver classification
# ------------------------------------------------------------

driver_rows <- list()
kk <- 0L

for (gene_name in c("PDGFA", "PDGFB")) {

  for (
    state_name in unique(
      assoc$sender_state[
        assoc$gene ==
          gene_name
      ]
    )
  ) {

    for (
      endpoint_name in c(
        "aHSC_ECM_remodeling",
        "wholeHSC_ECM_remodeling",
        "aHSC_core_collagen",
        "wholeHSC_core_collagen"
      )
    ) {

      b <- assoc %>%
        filter(
          gene == gene_name,
          sender_state == state_name,
          component == "ligand_burden",
          endpoint == endpoint_name
        )

      f <- assoc %>%
        filter(
          gene == gene_name,
          sender_state == state_name,
          component == "state_fraction",
          endpoint == endpoint_name
        )

      e <- assoc %>%
        filter(
          gene == gene_name,
          sender_state == state_name,
          component == "mean_CP10k",
          endpoint == endpoint_name
        )

      kk <- kk + 1L

      driver_rows[[kk]] <- tibble(
        gene = gene_name,
        sender_state = state_name,
        endpoint = endpoint_name,
        burden_rho = ifelse(nrow(b) == 1, b$rho_raw, NA_real_),
        burden_partial_NAS = ifelse(nrow(b) == 1, b$rho_partial_NAS, NA_real_),
        burden_nonzero_n = ifelse(nrow(b) == 1, b$nonzero_n, NA_integer_),
        fraction_rho = ifelse(nrow(f) == 1, f$rho_raw, NA_real_),
        fraction_partial_NAS = ifelse(nrow(f) == 1, f$rho_partial_NAS, NA_real_),
        expression_rho = ifelse(nrow(e) == 1, e$rho_raw, NA_real_),
        expression_partial_NAS = ifelse(nrow(e) == 1, e$rho_partial_NAS, NA_real_),
        component_driver =
          component_driver_class(
            b,
            f,
            e
          )
      )
    }
  }
}

driver_table <- bind_rows(
  driver_rows
) %>%
  arrange(
    gene,
    endpoint,
    desc(abs(burden_rho))
  )

write_tsv(
  driver_table,
  file.path(
    DIR_TABLE,
    "04_burden_composition_expression_driver.tsv"
  )
)

# ------------------------------------------------------------
# 5. PDGFA sample-level source matrix
# ------------------------------------------------------------

pdgfa_sample_matrix <- source_share %>%
  filter(
    gene == "PDGFA"
  ) %>%
  select(
    sample,
    sender_state,
    ligand_burden,
    state_fraction,
    mean_CP10k,
    detection_fraction,
    source_share
  ) %>%
  left_join(
    axis %>%
      select(
        sample,
        NAS,
        aHSC_ECM_remodeling,
        wholeHSC_ECM_remodeling,
        aHSC_core_collagen,
        wholeHSC_core_collagen,
        aHSC_fraction
      ),
    by = "sample"
  ) %>%
  arrange(
    sample,
    desc(source_share)
  )

write_tsv(
  pdgfa_sample_matrix,
  file.path(
    DIR_TABLE,
    "05_PDGFA_sample_source_matrix.tsv"
  )
)

# ------------------------------------------------------------
# 6. Primary rank summary
# ------------------------------------------------------------

primary_rank <- primary_pdgfa %>%
  group_by(
    endpoint
  ) %>%
  arrange(
    sparse_component,
    desc(abs(rho_partial_NAS)),
    desc(abs(rho_raw)),
    .by_group = TRUE
  ) %>%
  mutate(
    sender_rank =
      row_number()
  ) %>%
  ungroup()

write_tsv(
  primary_rank,
  file.path(
    DIR_TABLE,
    "06_PRIMARY_PDGFA_sender_rank.tsv"
  )
)

# ------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------

plot_df <- pdgfa_sample_matrix %>%
  filter(
    sender_state %in%
      c(
        "Resident/Kupffer-like_Mphi",
        "Lipid-associated/TREM2_Mphi",
        "Repair/Resolution_Mphi",
        "Inflammatory_Monocyte",
        "Inflammatory/Transitional_Mphi"
      )
  )

p1 <- ggplot(
  plot_df,
  aes(
    x = ligand_burden,
    y = aHSC_ECM_remodeling,
    label = sample
  )
) +
  geom_point(size = 2.5) +
  geom_text(
    nudge_y = 0.12,
    check_overlap = TRUE,
    size = 3
  ) +
  facet_wrap(
    ~ sender_state,
    scales = "free_x"
  ) +
  theme_classic(base_size = 10) +
  labs(
    title =
      "Human myeloid-state PDGFA burden vs aHSC ECM remodeling",
    x =
      "PDGFA burden (state fraction x mean CP10k)",
    y =
      "aHSC ECM remodeling score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_PDGFA_state_burden_vs_aHSC_ECM_remodeling.pdf"
  ),
  p1,
  width = 11,
  height = 7
)

# ------------------------------------------------------------
# 8. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "ligand_state_source",
    "sample_axis_source",
    "samples_n",
    "samples",
    "primary_gene",
    "primary_endpoints",
    "components_tested",
    "partial_method",
    "LOO_method",
    "sparse_definition",
    "inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    LIGAND_STATE_TSV,
    AXIS_NAS_TSV,
    as.character(nrow(axis)),
    paste(samples, collapse = ";"),
    "PDGFA",
    "aHSC_ECM_remodeling;wholeHSC_ECM_remodeling",
    paste(components, collapse = ";"),
    "rank residual partial Spearman after NAS",
    "leave-one-biological-sample-out Spearman",
    "component nonzero in <5/8 samples",
    "exploratory_sample_level_n8",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.2.0.tsv"
  )
)

saveRDS(
  list(
    source_share = source_share,
    associations = assoc,
    primary_pdgfa = primary_pdgfa,
    driver_table = driver_table,
    pdgfa_sample_matrix = pdgfa_sample_matrix,
    primary_rank = primary_rank
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_PDGFA_SenderState_Attribution_v1.2.0_results.rds"
  )
)

cat("\n============================================================\n")
cat("PRIMARY PDGFA SENDER-STATE ATTRIBUTION\n")
cat("============================================================\n")
print(primary_pdgfa)

cat("\n============================================================\n")
cat("BURDEN DRIVER CLASSIFICATION\n")
cat("============================================================\n")
print(
  driver_table %>%
    filter(
      gene == "PDGFA",
      endpoint %in%
        c(
          "aHSC_ECM_remodeling",
          "wholeHSC_ECM_remodeling"
        )
    )
)

cat("\n============================================================\n")
cat("PRIMARY PDGFA SENDER RANK\n")
cat("============================================================\n")
print(primary_rank)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
