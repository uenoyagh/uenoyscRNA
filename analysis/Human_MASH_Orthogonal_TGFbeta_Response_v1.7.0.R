#!/usr/bin/env Rscript

# ============================================================
# Human MASH orthogonal TGF-beta response audit v1.7.0
#
# PURPOSE
#   Remove circularity between the HSC ECM score and the
#   TGF-beta response score.
#
# PRIMARY orthogonal TGF-beta response genes:
#   SMAD7, PMEPA1, THBS1
#
# These genes are measurable in the HSC Flex panel and are NOT
# members of the ECM-core score used to define ECM-high HSC.
#
# SECONDARY extended direct-response genes:
#   SERPINE1, SMAD7, PMEPA1, TGFBI, THBS1
#
# SERPINE1 and TGFBI overlap the ECM-core score and therefore
# are retained only as secondary/sensitivity readouts.
#
# Analyses:
#   1. Build orthogonal 3-gene TGF-beta response score
#   2. Test TREM2-Mphi TGFB1 burden vs orthogonal response
#   3. Exact permutation test (8! = 40,320)
#   4. Leave-one-sample-out
#   5. Partial Spearman controlling for NAS
#   6. Compare all 5 myeloid sender states
#   7. Gene-level target concordance for SMAD7/PMEPA1/THBS1
#   8. ECM-high vs ECM-nonhigh paired comparison using only
#      orthogonal response genes
#
# IMPORTANT
#   - Sample is the statistical unit.
#   - No reclustering.
#   - Source RDS is never overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
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

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2",
  "objects",
  "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
)

SENDER_TABLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5",
  "tables",
  "02_TGFB1_sender_state_burden_by_sample.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Orthogonal_TGFbeta_Response_v1.7"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_OBJ   <- file.path(OUT, "objects")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV,  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_Orthogonal_TGFbeta_Response_v1.7.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH orthogonal TGF-beta response v1.7.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(HSC_RDS))
stopifnot(file.exists(SENDER_TABLE))

# ------------------------------------------------------------
# 1. Load
# ------------------------------------------------------------

hsc <- readRDS(HSC_RDS)

DefaultAssay(hsc) <- "RNA"

sender <- read_tsv(
  SENDER_TABLE,
  show_col_types = FALSE
)

cat("HSC cells:", ncol(hsc), "\n")
cat("HSC features:", nrow(hsc), "\n\n")

if (ncol(hsc) != 7640) {
  warning(
    "Expected 7640 canonical HSC cells; observed ",
    ncol(hsc)
  )
}

required_hsc <- c(
  "sample",
  "NAS_score",
  "HSC_state_v1.2",
  "ECM_high_primary_v1.2",
  "ECM_total_score_v1.2"
)

missing_hsc <- setdiff(
  required_hsc,
  colnames(hsc[[]])
)

if (length(missing_hsc) > 0) {
  stop(
    "Missing HSC metadata: ",
    paste(missing_hsc, collapse = ";")
  )
}

required_sender <- c(
  "sample",
  "NAS",
  "sender_state",
  "TGFB1_burden_expr"
)

missing_sender <- setdiff(
  required_sender,
  colnames(sender)
)

if (length(missing_sender) > 0) {
  stop(
    "Missing sender columns: ",
    paste(missing_sender, collapse = ";")
  )
}

# ------------------------------------------------------------
# 2. Normalize RNA
# ------------------------------------------------------------

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 3. Define orthogonal vs extended panels
# ------------------------------------------------------------

ORTHOGONAL_TGF_RESPONSE <- c(
  "SMAD7",
  "PMEPA1",
  "THBS1"
)

EXTENDED_TGF_RESPONSE <- c(
  "SERPINE1",
  "SMAD7",
  "PMEPA1",
  "TGFBI",
  "THBS1"
)

ECM_OVERLAP_GENES <- c(
  "SERPINE1",
  "TGFBI"
)

availability <- bind_rows(
  tibble(
    panel = "ORTHOGONAL_PRIMARY",
    gene = ORTHOGONAL_TGF_RESPONSE,
    overlaps_ECM_core = FALSE
  ),
  tibble(
    panel = "EXTENDED_SECONDARY",
    gene = EXTENDED_TGF_RESPONSE,
    overlaps_ECM_core =
      EXTENDED_TGF_RESPONSE %in%
      ECM_OVERLAP_GENES
  )
) %>%
  mutate(
    available_in_HSC =
      gene %in%
      rownames(hsc)
  )

write_tsv(
  availability,
  file.path(
    DIR_TABLE,
    "01_TGFbeta_panel_availability_and_overlap.tsv"
  )
)

orthogonal_genes <- availability %>%
  filter(
    panel == "ORTHOGONAL_PRIMARY",
    available_in_HSC
  ) %>%
  pull(gene)

extended_genes <- availability %>%
  filter(
    panel == "EXTENDED_SECONDARY",
    available_in_HSC
  ) %>%
  pull(gene)

if (length(orthogonal_genes) != 3) {
  stop(
    "Expected all 3 orthogonal genes to be measurable; found ",
    length(orthogonal_genes)
  )
}

cat("Orthogonal genes:\n")
print(orthogonal_genes)
cat("\n")

# ------------------------------------------------------------
# 4. Scoring helper
# ------------------------------------------------------------

z_fun <- function(x) {

  sx <- sd(
    x,
    na.rm = TRUE
  )

  if (
    is.na(sx) ||
    sx == 0
  ) {
    return(
      rep(0, length(x))
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

calc_score <- function(
  object,
  genes
) {

  x <- FetchData(
    object,
    vars = genes,
    layer = "data"
  )

  xz <- as.data.frame(
    lapply(
      x,
      z_fun
    )
  )

  rowMeans(
    xz,
    na.rm = TRUE
  )
}

# ------------------------------------------------------------
# 5. Compute orthogonal and extended scores
# ------------------------------------------------------------

hsc$TGF_orthogonal_score_v1.7 <-
  calc_score(
    hsc,
    orthogonal_genes
  )

hsc$TGF_extended_score_v1.7 <-
  calc_score(
    hsc,
    extended_genes
  )

# ------------------------------------------------------------
# 6. Cell-level expression table for direct target genes
# ------------------------------------------------------------

target_expr <- FetchData(
  hsc,
  vars = unique(
    c(
      orthogonal_genes,
      extended_genes
    )
  ),
  layer = "data"
)

hmd <- hsc[[]]

cell_df <- hmd %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample),
    NAS = suppressWarnings(
      as.numeric(NAS_score)
    ),
    HSC_state = as.character(
      HSC_state_v1.2
    ),
    ECM_high = as.logical(
      ECM_high_primary_v1.2
    ),
    ECM_total_score = as.numeric(
      ECM_total_score_v1.2
    ),
    TGF_orthogonal_score = as.numeric(
      TGF_orthogonal_score_v1.7
    ),
    TGF_extended_score = as.numeric(
      TGF_extended_score_v1.7
    )
  )

for (g in colnames(target_expr)) {
  cell_df[[g]] <- target_expr[[g]]
}

# ------------------------------------------------------------
# 7. Sample-level HSC summaries
# ------------------------------------------------------------

sample_hsc <- cell_df %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    HSC_n = n(),

    ECM_high_fraction =
      mean(
        ECM_high,
        na.rm = TRUE
      ),

    ECM_total_mean =
      mean(
        ECM_total_score,
        na.rm = TRUE
      ),

    TGF_orthogonal_mean =
      mean(
        TGF_orthogonal_score,
        na.rm = TRUE
      ),

    TGF_extended_mean =
      mean(
        TGF_extended_score,
        na.rm = TRUE
      ),

    across(
      all_of(
        unique(
          c(
            orthogonal_genes,
            extended_genes
          )
        )
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      ),
      .names = "{.col}_mean"
    ),

    .groups = "drop"
  ) %>%
  arrange(sample)

write_tsv(
  sample_hsc,
  file.path(
    DIR_TABLE,
    "02_HSC_TGFbeta_targets_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 8. ECM-high vs ECM-nonhigh summaries
# ------------------------------------------------------------

ecm_group <- cell_df %>%
  mutate(
    ECM_group = ifelse(
      ECM_high,
      "ECM_high",
      "ECM_nonhigh"
    )
  ) %>%
  group_by(
    sample,
    NAS,
    ECM_group
  ) %>%
  summarise(
    cells = n(),

    TGF_orthogonal_mean =
      mean(
        TGF_orthogonal_score,
        na.rm = TRUE
      ),

    TGF_extended_mean =
      mean(
        TGF_extended_score,
        na.rm = TRUE
      ),

    across(
      all_of(
        unique(
          c(
            orthogonal_genes,
            extended_genes
          )
        )
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      ),
      .names = "{.col}_mean"
    ),

    .groups = "drop"
  )

write_tsv(
  ecm_group,
  file.path(
    DIR_TABLE,
    "03_ECMhigh_vs_nonhigh_TGFbeta_targets_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 9. Paired ECM-high vs nonhigh tests
# ------------------------------------------------------------

paired_wide <- ecm_group %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = c(
      cells,
      TGF_orthogonal_mean,
      TGF_extended_mean,
      all_of(
        paste0(
          unique(
            c(
              orthogonal_genes,
              extended_genes
            )
          ),
          "_mean"
        )
      )
    )
  )

paired_metric_names <- c(
  "TGF_orthogonal_mean",
  paste0(
    orthogonal_genes,
    "_mean"
  )
)

paired_results <- bind_rows(
  lapply(
    paired_metric_names,
    function(metric) {

      high_col <- paste0(
        metric,
        "_ECM_high"
      )

      non_col <- paste0(
        metric,
        "_ECM_nonhigh"
      )

      x <- paired_wide[[high_col]]
      y <- paired_wide[[non_col]]

      keep <- complete.cases(
        x,
        y
      )

      if (
        sum(keep) < 4 ||
        all(
          x[keep] ==
            y[keep]
        )
      ) {
        p <- NA_real_
      } else {
        p <- suppressWarnings(
          wilcox.test(
            x[keep],
            y[keep],
            paired = TRUE,
            exact = FALSE
          )$p.value
        )
      }

      tibble(
        metric = metric,
        n_samples = sum(keep),
        median_delta_ECMhigh_minus_nonhigh =
          median(
            x[keep] -
              y[keep],
            na.rm = TRUE
          ),
        positive_delta_samples =
          sum(
            x[keep] >
              y[keep],
            na.rm = TRUE
          ),
        p_value = p
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
  paired_results,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_orthogonal_ECMhigh_paired_tests.tsv"
  )
)

# ------------------------------------------------------------
# 10. Statistical helpers
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
    rho = unname(
      ct$estimate
    ),
    p_value = ct$p.value
  )
}

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

    p_rest <- all_permutations(
      x[-i]
    )

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

exact_perm <- function(
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
      "Exact permutation restricted to n <= 8."
    )
  }

  rho_obs <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  perms <- all_permutations(
    seq_len(n)
  )

  rhos <- apply(
    perms,
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

  tibble(
    n = n,
    rho = rho_obs,
    permutations = length(rhos),
    p_exact =
      mean(
        abs(rhos) >=
          abs(rho_obs) -
          1e-12
      )
  )
}

partial_spearman <- function(
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
    denom <= 0 ||
    n < 5
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  rp <- (
    r_xy -
      r_xz *
      r_yz
  ) / denom

  rp <- max(
    min(
      rp,
      0.999999
    ),
    -0.999999
  )

  df <- n - 3

  tval <- rp *
    sqrt(
      df /
        (1 - rp^2)
    )

  p <- 2 *
    pt(
      -abs(tval),
      df = df
    )

  tibble(
    n = n,
    partial_rho = rp,
    p_value = p
  )
}

# ------------------------------------------------------------
# 11. TREM2 sender merge
# ------------------------------------------------------------

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

trem2 <- sender %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  select(
    sample,
    NAS,
    TGFB1_burden_expr,
    TGFB1_expression_share,
    TGFB1_mean,
    state_fraction
  ) %>%
  inner_join(
    sample_hsc,
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

# ------------------------------------------------------------
# 12. PRIMARY orthogonal module association
# ------------------------------------------------------------

primary_spearman <- safe_spearman(
  trem2$TGFB1_burden_expr,
  trem2$TGF_orthogonal_mean
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "TGF_orthogonal_mean",
    .before = 1
  )

primary_exact <- exact_perm(
  trem2$TGFB1_burden_expr,
  trem2$TGF_orthogonal_mean
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "TGF_orthogonal_mean",
    .before = 1
  )

primary_partial <- partial_spearman(
  trem2$TGFB1_burden_expr,
  trem2$TGF_orthogonal_mean,
  trem2$NAS
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "TGF_orthogonal_mean",
    covariate =
      "NAS",
    .before = 1
  )

write_tsv(
  primary_spearman,
  file.path(
    DIR_TABLE,
    "05_PRIMARY_TREM2_TGFB1_vs_orthogonal_TGF_Spearman.tsv"
  )
)

write_tsv(
  primary_exact,
  file.path(
    DIR_TABLE,
    "06_PRIMARY_TREM2_TGFB1_vs_orthogonal_TGF_exact_permutation.tsv"
  )
)

write_tsv(
  primary_partial,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_partial_Spearman_controlling_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 13. Leave-one-out for orthogonal module
# ------------------------------------------------------------

loo <- bind_rows(
  lapply(
    seq_len(nrow(trem2)),
    function(i) {

      sub <- trem2[-i, ]

      safe_spearman(
        sub$TGFB1_burden_expr,
        sub$TGF_orthogonal_mean
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
    "08_leave_one_sample_out_orthogonal_score.tsv"
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
    "09_leave_one_sample_out_orthogonal_summary.tsv"
  )
)

# ------------------------------------------------------------
# 14. Gene-level orthogonal target concordance
# ------------------------------------------------------------

gene_results <- bind_rows(
  lapply(
    orthogonal_genes,
    function(g) {

      metric <- paste0(
        g,
        "_mean"
      )

      asym <- safe_spearman(
        trem2$TGFB1_burden_expr,
        trem2[[metric]]
      )

      ex <- exact_perm(
        trem2$TGFB1_burden_expr,
        trem2[[metric]]
      )

      pa <- partial_spearman(
        trem2$TGFB1_burden_expr,
        trem2[[metric]],
        trem2$NAS
      )

      tibble(
        gene = g,
        n = asym$n,
        rho = asym$rho,
        p_value = asym$p_value,
        p_exact = ex$p_exact,
        partial_rho_NAS = pa$partial_rho,
        partial_p_NAS = pa$p_value
      )
    }
  )
) %>%
  mutate(
    FDR_asymptotic =
      p.adjust(
        p_value,
        method = "BH"
      ),

    FDR_exact =
      p.adjust(
        p_exact,
        method = "BH"
      ),

    FDR_partial_NAS =
      p.adjust(
        partial_p_NAS,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR_exact,
    p_exact
  )

write_tsv(
  gene_results,
  file.path(
    DIR_TABLE,
    "10_PRIMARY_orthogonal_target_gene_concordance.tsv"
  )
)

# ------------------------------------------------------------
# 15. Compare all sender states using orthogonal score
# ------------------------------------------------------------

all_sender <- sender %>%
  inner_join(
    sample_hsc %>%
      select(
        sample,
        NAS,
        TGF_orthogonal_mean
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
      .x$TGF_orthogonal_mean
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
  all_sender,
  file.path(
    DIR_TABLE,
    "11_all_sender_states_vs_orthogonal_TGF_score.tsv"
  )
)

# ------------------------------------------------------------
# 16. Sensitivity: extended 5-gene score
# ------------------------------------------------------------

extended_sensitivity <- safe_spearman(
  trem2$TGFB1_burden_expr,
  trem2$TGF_extended_mean
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "TGF_extended_mean",
    note =
      "secondary_only_SERPINE1_and_TGFBI_overlap_ECM_core",
    .before = 1
  )

write_tsv(
  extended_sensitivity,
  file.path(
    DIR_TABLE,
    "12_SECONDARY_extended_TGF_score.tsv"
  )
)

# ------------------------------------------------------------
# 17. NAS group descriptive summary
# ------------------------------------------------------------

NAS_summary <- trem2 %>%
  group_by(NAS) %>%
  summarise(
    samples = n(),

    TREM2_TGFB1_burden =
      mean(
        TGFB1_burden_expr,
        na.rm = TRUE
      ),

    TREM2_TGFB1_share =
      mean(
        TGFB1_expression_share,
        na.rm = TRUE
      ),

    TGF_orthogonal_mean =
      mean(
        TGF_orthogonal_mean,
        na.rm = TRUE
      ),

    TGF_extended_mean =
      mean(
        TGF_extended_mean,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high_fraction,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  NAS_summary,
  file.path(
    DIR_TABLE,
    "13_NAS_mean_orthogonal_TGF_summary.tsv"
  )
)

# ------------------------------------------------------------
# 18. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  trem2,
  aes(
    x = TGFB1_burden_expr,
    y = TGF_orthogonal_mean,
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
    title = "TREM2-Mphi TGFB1 burden vs orthogonal HSC TGF-beta response",
    subtitle = "Orthogonal score: SMAD7 + PMEPA1 + THBS1",
    x = "TREM2-Mphi TGFB1 burden",
    y = "HSC orthogonal TGF-beta response score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TREM2_TGFB1_vs_orthogonal_TGF_response.pdf"
  ),
  p1,
  width = 7,
  height = 5.5
)

p2 <- ggplot(
  ecm_group,
  aes(
    x = ECM_group,
    y = TGF_orthogonal_mean,
    group = sample
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Orthogonal TGF-beta response in ECM-high HSC",
    subtitle = "No score genes overlap the ECM-core definition",
    x = NULL,
    y = "Orthogonal TGF-beta response score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_orthogonal_TGF_response_ECMhigh_vs_nonhigh.pdf"
  ),
  p2,
  width = 6.5,
  height = 5.5
)

p3 <- ggplot(
  all_sender,
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
    title = "TGFB1 burden vs orthogonal HSC TGF-beta response",
    x = NULL,
    y = "Spearman rho"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_sender_state_comparison_orthogonal_TGF_response.pdf"
  ),
  p3,
  width = 8,
  height = 5.5
)

# ------------------------------------------------------------
# 19. Save derived object
# ------------------------------------------------------------

saveRDS(
  list(
    availability = availability,
    sample_hsc = sample_hsc,
    ecm_group = ecm_group,
    paired_results = paired_results,
    trem2 = trem2,
    primary_spearman = primary_spearman,
    primary_exact = primary_exact,
    primary_partial = primary_partial,
    loo = loo,
    gene_results = gene_results,
    all_sender = all_sender,
    NAS_summary = NAS_summary
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_Orthogonal_TGFbeta_Response_v1.7_results.rds"
  )
)

# ------------------------------------------------------------
# 20. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "HSC_source_RDS",
    "sender_source_table",
    "HSC_cells",
    "primary_TGF_response_genes",
    "primary_ECM_overlap",
    "secondary_extended_genes",
    "secondary_ECM_overlap_genes",
    "primary_sender_state",
    "primary_sender_metric",
    "primary_association_unit",
    "exact_permutations",
    "reclustering",
    "source_overwritten"
  ),
  value = c(
    "v1.7.0",
    HSC_RDS,
    SENDER_TABLE,
    as.character(
      ncol(hsc)
    ),
    paste(
      orthogonal_genes,
      collapse = ";"
    ),
    "none",
    paste(
      extended_genes,
      collapse = ";"
    ),
    paste(
      intersect(
        extended_genes,
        ECM_OVERLAP_GENES
      ),
      collapse = ";"
    ),
    TREM2_STATE,
    "TGFB1_burden_expr",
    "biological_sample",
    "40320_for_n8",
    "none",
    "FALSE"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "orthogonal_TGFbeta_provenance_v1.7.tsv"
  )
)

# ------------------------------------------------------------
# 21. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "orthogonal_genes_n",
    "orthogonal_genes_overlap_ECM_core",
    "TREM2_samples",
    "all_8_samples_present",
    "exact_permutation_completed",
    "LOO_tests",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      ncol(hsc)
    ),
    "7640",
    as.character(
      ncol(hsc) == 7640
    ),
    as.character(
      length(
        orthogonal_genes
      )
    ),
    "FALSE",
    as.character(
      nrow(trem2)
    ),
    as.character(
      nrow(trem2) == 8
    ),
    as.character(
      primary_exact$permutations ==
        40320
    ),
    as.character(
      nrow(loo)
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.7.tsv"
  )
)

# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("ORTHOGONAL PANEL AUDIT\n")
cat("============================================================\n")
print(availability)

cat("\n")
cat("============================================================\n")
cat("PRIMARY ORTHOGONAL ASSOCIATION\n")
cat("============================================================\n")
print(primary_spearman)
print(primary_exact)
print(primary_partial)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT\n")
cat("============================================================\n")
print(loo)
print(loo_summary)

cat("\n")
cat("============================================================\n")
cat("ORTHOGONAL TARGET GENES\n")
cat("============================================================\n")
print(gene_results)

cat("\n")
cat("============================================================\n")
cat("ECM-HIGH vs ECM-NONHIGH: ORTHOGONAL TARGETS\n")
cat("============================================================\n")
print(paired_results)

cat("\n")
cat("============================================================\n")
cat("ALL SENDER STATES\n")
cat("============================================================\n")
print(all_sender)

cat("\n")
cat("============================================================\n")
cat("NAS SUMMARY\n")
cat("============================================================\n")
print(NAS_summary)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
