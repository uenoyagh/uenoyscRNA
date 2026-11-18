#!/usr/bin/env Rscript

# ============================================================
# Human MASH HSC activation-state matrix v1.4.0
#
# PURPOSE
#   Refine HSC functional-state structure after v1.3.2 showed:
#     - qHSC identity is supported mainly by quiescent-marker enrichment
#     - aHSC_mesenchymal shows consistent IEG/initiatory activation
#     - myofibroblastic identity is NOT uniformly elevated in aHSC
#     - true cycling HSCs are extremely rare
#
# PRIMARY MODEL
#   Keep identity annotation unchanged:
#     qHSC / aHSC_mesenchymal
#
#   Overlay two independent functional axes:
#     1) IEG / initiatory activation
#     2) ECM / fibrogenic output
#
#   Functional states:
#     Low
#     IEG_high_only
#     ECM_high_only
#     IEG_ECM_dual_high
#
# PRIMARY THRESHOLD
#   Within each biological sample:
#     score genes by within-sample z-score
#     define IEG-high and ECM-high using the sample-matched qHSC p90
#
# SENSITIVITY
#   Global gene-z scores + pooled qHSC p90 thresholds
#
# IMPORTANT
#   These are functional states, not replacement HSC subtypes.
#   No reclustering.
#   No source RDS overwrite.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

set.seed(20261118)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2",
  "objects",
  "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_Activation_State_Matrix_v1.4"
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
  "Human_MASH_HSC_Activation_State_Matrix_v1.4.0.log"
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
cat("Human MASH HSC activation-state matrix v1.4.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(HSC_RDS)
)

hsc <- readRDS(HSC_RDS)

DefaultAssay(hsc) <- "RNA"

if (ncol(hsc) != 7640) {
  warning(
    "Expected canonical HSC n=7640; observed ",
    ncol(hsc)
  )
}

required_meta <- c(
  "sample",
  "NAS_score",
  "HSC_state_v1.2"
)

missing_meta <- setdiff(
  required_meta,
  colnames(hsc[[]])
)

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata: ",
    paste(
      missing_meta,
      collapse = ";"
    )
  )
}

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 1. Functional gene panels
# ------------------------------------------------------------

IEG_GENES <- c(
  "THBS1",
  "TIMP1",
  "JUNB",
  "FOS",
  "FOSB"
)

ECM_GENES <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2",
  "FN1",
  "SPARC",
  "THBS2",
  "LOXL2",
  "TIMP1",
  "TGFBI",
  "CTHRC1",
  "SERPINE1",
  "POSTN"
)

QUIESCENT_GENES <- c(
  "RBP1",
  "RELN",
  "CYGB",
  "NGFR",
  "LHX2",
  "COLEC10"
)

panel_audit <- bind_rows(
  tibble(
    panel = "IEG_Initiatory",
    gene = IEG_GENES,
    available =
      IEG_GENES %in%
      rownames(hsc)
  ),
  tibble(
    panel = "ECM_Fibrogenic",
    gene = ECM_GENES,
    available =
      ECM_GENES %in%
      rownames(hsc)
  ),
  tibble(
    panel = "Quiescent",
    gene = QUIESCENT_GENES,
    available =
      QUIESCENT_GENES %in%
      rownames(hsc)
  )
)

write_tsv(
  panel_audit,
  file.path(
    DIR_TABLE,
    "01_panel_gene_availability.tsv"
  )
)

panel_summary <- panel_audit %>%
  group_by(panel) %>%
  summarise(
    requested_n = n(),
    available_n =
      sum(available),
    coverage_fraction =
      mean(available),
    available_genes =
      paste(
        gene[available],
        collapse = ";"
      ),
    missing_genes =
      paste(
        gene[!available],
        collapse = ";"
      ),
    .groups = "drop"
  )

write_tsv(
  panel_summary,
  file.path(
    DIR_TABLE,
    "02_panel_coverage_summary.tsv"
  )
)

IEG_AVAIL <- intersect(
  IEG_GENES,
  rownames(hsc)
)

ECM_AVAIL <- intersect(
  ECM_GENES,
  rownames(hsc)
)

Q_AVAIL <- intersect(
  QUIESCENT_GENES,
  rownames(hsc)
)

if (
  length(IEG_AVAIL) < 4 ||
  length(ECM_AVAIL) < 8 ||
  length(Q_AVAIL) < 4
) {
  stop(
    "Insufficient panel coverage."
  )
}

# ------------------------------------------------------------
# 2. Expression matrix
# ------------------------------------------------------------

ALL_SCORE_GENES <- unique(
  c(
    IEG_AVAIL,
    ECM_AVAIL,
    Q_AVAIL
  )
)

expr <- FetchData(
  hsc,
  vars = ALL_SCORE_GENES,
  layer = "data"
)

md <- hsc[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample =
      as.character(sample),
    NAS =
      suppressWarnings(
        as.numeric(
          NAS_score
        )
      ),
    current_HSC_state =
      as.character(
        HSC_state_v1.2
      )
  )

stopifnot(
  identical(
    rownames(expr),
    md$cell
  )
)

# ------------------------------------------------------------
# 3. Scoring helpers
# ------------------------------------------------------------

safe_z <- function(x) {

  sx <- sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(sx) ||
    sx == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
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

score_from_matrix <- function(
  mat,
  genes
) {

  x <- mat[
    ,
    genes,
    drop = FALSE
  ]

  rowMeans(
    x,
    na.rm = TRUE
  )
}

# ------------------------------------------------------------
# 4. PRIMARY: within-sample gene z-scores
# ------------------------------------------------------------

within_sample_scored <- bind_rows(
  lapply(
    unique(
      md$sample
    ),
    function(s) {

      idx <- which(
        md$sample == s
      )

      x <- expr[
        idx,
        ALL_SCORE_GENES,
        drop = FALSE
      ]

      xz <- as.data.frame(
        lapply(
          x,
          safe_z
        )
      )

      tibble(
        cell =
          rownames(x),

        sample = s,

        IEG_score_ws =
          score_from_matrix(
            xz,
            IEG_AVAIL
          ),

        ECM_score_ws =
          score_from_matrix(
            xz,
            ECM_AVAIL
          ),

        Quiescent_score_ws =
          score_from_matrix(
            xz,
            Q_AVAIL
          )
      )
    }
  )
)

# ------------------------------------------------------------
# 5. Sample-matched qHSC p90 thresholds
# ------------------------------------------------------------

primary_thresholds <- within_sample_scored %>%
  left_join(
    md %>%
      select(
        cell,
        current_HSC_state,
        NAS
      ),
    by = "cell"
  ) %>%
  filter(
    current_HSC_state ==
      "qHSC"
  ) %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    qHSC_n =
      n(),

    IEG_qHSC_p90 =
      quantile(
        IEG_score_ws,
        0.90,
        na.rm = TRUE
      ),

    ECM_qHSC_p90 =
      quantile(
        ECM_score_ws,
        0.90,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  primary_thresholds,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_sample_matched_qHSC_p90_thresholds.tsv"
  )
)

primary_cells <- within_sample_scored %>%
  left_join(
    md,
    by = c(
      "cell",
      "sample"
    )
  ) %>%
  left_join(
    primary_thresholds,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  mutate(
    IEG_high_primary =
      IEG_score_ws >
      IEG_qHSC_p90,

    ECM_high_primary =
      ECM_score_ws >
      ECM_qHSC_p90,

    functional_state_primary =
      case_when(
        IEG_high_primary &
          ECM_high_primary ~
          "IEG_ECM_dual_high",

        IEG_high_primary &
          !ECM_high_primary ~
          "IEG_high_only",

        !IEG_high_primary &
          ECM_high_primary ~
          "ECM_high_only",

        TRUE ~
          "Low"
      )
  )

# ------------------------------------------------------------
# 6. SENSITIVITY: global gene-z + pooled qHSC p90
# ------------------------------------------------------------

global_z <- as.data.frame(
  lapply(
    expr,
    safe_z
  )
)

sensitivity_cells <- md %>%
  mutate(
    IEG_score_global =
      score_from_matrix(
        global_z,
        IEG_AVAIL
      ),

    ECM_score_global =
      score_from_matrix(
        global_z,
        ECM_AVAIL
      ),

    Quiescent_score_global =
      score_from_matrix(
        global_z,
        Q_AVAIL
      )
  )

pooled_q <- sensitivity_cells %>%
  filter(
    current_HSC_state ==
      "qHSC"
  )

IEG_q_p90_global <- quantile(
  pooled_q$IEG_score_global,
  0.90,
  na.rm = TRUE
)

ECM_q_p90_global <- quantile(
  pooled_q$ECM_score_global,
  0.90,
  na.rm = TRUE
)

sensitivity_cells <- sensitivity_cells %>%
  mutate(
    IEG_high_global =
      IEG_score_global >
      IEG_q_p90_global,

    ECM_high_global =
      ECM_score_global >
      ECM_q_p90_global,

    functional_state_global =
      case_when(
        IEG_high_global &
          ECM_high_global ~
          "IEG_ECM_dual_high",

        IEG_high_global &
          !ECM_high_global ~
          "IEG_high_only",

        !IEG_high_global &
          ECM_high_global ~
          "ECM_high_only",

        TRUE ~
          "Low"
      )
  )

# ------------------------------------------------------------
# 7. Combine cell-level calls
# ------------------------------------------------------------

cell_calls <- primary_cells %>%
  left_join(
    sensitivity_cells %>%
      select(
        cell,
        IEG_score_global,
        ECM_score_global,
        Quiescent_score_global,
        IEG_high_global,
        ECM_high_global,
        functional_state_global
      ),
    by = "cell"
  ) %>%
  mutate(
    primary_vs_global_state_agree =
      functional_state_primary ==
      functional_state_global
  )

write_tsv(
  cell_calls,
  file.path(
    DIR_TABLE,
    "04_cell_level_activation_state_calls.tsv.gz"
  )
)

# Save metadata into derived object only.
hsc$HSC_IEG_score_ws_v1.4 <-
  cell_calls[
    colnames(hsc),
    "IEG_score_ws",
    drop = TRUE
  ]

hsc$HSC_ECM_score_ws_v1.4 <-
  cell_calls[
    colnames(hsc),
    "ECM_score_ws",
    drop = TRUE
  ]

hsc$HSC_functional_state_v1.4 <-
  cell_calls[
    colnames(hsc),
    "functional_state_primary",
    drop = TRUE
  ]

# ------------------------------------------------------------
# 8. State fractions by identity and sample
# ------------------------------------------------------------

state_by_identity <- cell_calls %>%
  count(
    current_HSC_state,
    functional_state_primary,
    name = "cells"
  ) %>%
  group_by(
    current_HSC_state
  ) %>%
  mutate(
    fraction =
      cells /
      sum(cells)
  ) %>%
  ungroup()

write_tsv(
  state_by_identity,
  file.path(
    DIR_TABLE,
    "05_functional_state_by_current_identity.tsv"
  )
)

sample_state <- cell_calls %>%
  group_by(
    sample,
    NAS,
    current_HSC_state
  ) %>%
  summarise(
    cells = n(),

    Low_fraction =
      mean(
        functional_state_primary ==
          "Low"
      ),

    IEG_high_only_fraction =
      mean(
        functional_state_primary ==
          "IEG_high_only"
      ),

    ECM_high_only_fraction =
      mean(
        functional_state_primary ==
          "ECM_high_only"
      ),

    Dual_high_fraction =
      mean(
        functional_state_primary ==
          "IEG_ECM_dual_high"
      ),

    Any_IEG_high_fraction =
      mean(
        IEG_high_primary
      ),

    Any_ECM_high_fraction =
      mean(
        ECM_high_primary
      ),

    mean_IEG_score =
      mean(
        IEG_score_ws,
        na.rm = TRUE
      ),

    mean_ECM_score =
      mean(
        ECM_score_ws,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  sample_state,
  file.path(
    DIR_TABLE,
    "06_sample_by_identity_functional_state_fractions.tsv"
  )
)

# ------------------------------------------------------------
# 9. PRIMARY focus: aHSC functional states
# ------------------------------------------------------------

aHSC_state <- sample_state %>%
  filter(
    current_HSC_state ==
      "aHSC_mesenchymal"
  ) %>%
  arrange(sample)

write_tsv(
  aHSC_state,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_aHSC_functional_state_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 10. Exact Spearman permutation helper
# ------------------------------------------------------------

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

    out[[i]] <- cbind(
      x[i],
      all_permutations(
        x[-i]
      )
    )
  }

  do.call(
    rbind,
    out
  )
}

exact_spearman <- function(
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
    p_exact =
      mean(
        abs(rhos) >=
          abs(rho_obs) -
          1e-12
      )
  )
}

# ------------------------------------------------------------
# 11. NAS associations for aHSC functional state
# ------------------------------------------------------------

NAS_METRICS <- c(
  "Low_fraction",
  "IEG_high_only_fraction",
  "ECM_high_only_fraction",
  "Dual_high_fraction",
  "Any_IEG_high_fraction",
  "Any_ECM_high_fraction",
  "mean_IEG_score",
  "mean_ECM_score"
)

nas_assoc <- bind_rows(
  lapply(
    NAS_METRICS,
    function(metric) {

      exact_spearman(
        aHSC_state$NAS,
        aHSC_state[[metric]]
      ) %>%
        mutate(
          metric = metric,
          .before = 1
        )
    }
  )
) %>%
  mutate(
    FDR =
      p.adjust(
        p_exact,
        method = "BH"
      )
  ) %>%
  arrange(FDR)

write_tsv(
  nas_assoc,
  file.path(
    DIR_TABLE,
    "08_PRIMARY_aHSC_functional_state_vs_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 12. Sample-matched aHSC vs qHSC paired deltas
# ------------------------------------------------------------

paired_wide <- sample_state %>%
  select(
    sample,
    NAS,
    current_HSC_state,
    Any_IEG_high_fraction,
    Any_ECM_high_fraction,
    Dual_high_fraction,
    mean_IEG_score,
    mean_ECM_score
  ) %>%
  pivot_wider(
    names_from =
      current_HSC_state,
    values_from = c(
      Any_IEG_high_fraction,
      Any_ECM_high_fraction,
      Dual_high_fraction,
      mean_IEG_score,
      mean_ECM_score
    )
  )

paired_deltas <- paired_wide %>%
  transmute(
    sample,
    NAS,

    IEG_high_fraction_delta =
      Any_IEG_high_fraction_aHSC_mesenchymal -
      Any_IEG_high_fraction_qHSC,

    ECM_high_fraction_delta =
      Any_ECM_high_fraction_aHSC_mesenchymal -
      Any_ECM_high_fraction_qHSC,

    Dual_high_fraction_delta =
      Dual_high_fraction_aHSC_mesenchymal -
      Dual_high_fraction_qHSC,

    IEG_score_delta =
      mean_IEG_score_aHSC_mesenchymal -
      mean_IEG_score_qHSC,

    ECM_score_delta =
      mean_ECM_score_aHSC_mesenchymal -
      mean_ECM_score_qHSC
  )

write_tsv(
  paired_deltas,
  file.path(
    DIR_TABLE,
    "09_sample_matched_aHSC_minus_qHSC_deltas.tsv"
  )
)

paired_test <- function(x) {

  wt <- suppressWarnings(
    wilcox.test(
      x,
      mu = 0,
      exact = FALSE,
      alternative = "two.sided"
    )
  )

  tibble(
    n_samples =
      sum(
        is.finite(x)
      ),

    median_delta =
      median(
        x,
        na.rm = TRUE
      ),

    positive_samples =
      sum(
        x > 0,
        na.rm = TRUE
      ),

    negative_samples =
      sum(
        x < 0,
        na.rm = TRUE
      ),

    p_value =
      wt$p.value
  )
}

paired_summary <- bind_rows(
  lapply(
    setdiff(
      colnames(
        paired_deltas
      ),
      c(
        "sample",
        "NAS"
      )
    ),
    function(metric) {

      paired_test(
        paired_deltas[[metric]]
      ) %>%
        mutate(
          metric = metric,
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
  paired_summary,
  file.path(
    DIR_TABLE,
    "10_PRIMARY_aHSC_vs_qHSC_functional_state_tests.tsv"
  )
)

# ------------------------------------------------------------
# 13. Primary vs sensitivity agreement
# ------------------------------------------------------------

agreement <- cell_calls %>%
  group_by(
    current_HSC_state
  ) %>%
  summarise(
    cells = n(),

    exact_state_agreement =
      mean(
        primary_vs_global_state_agree
      ),

    IEG_call_agreement =
      mean(
        IEG_high_primary ==
          IEG_high_global
      ),

    ECM_call_agreement =
      mean(
        ECM_high_primary ==
          ECM_high_global
      ),

    .groups = "drop"
  )

write_tsv(
  agreement,
  file.path(
    DIR_TABLE,
    "11_primary_vs_global_threshold_agreement.tsv"
  )
)

# ------------------------------------------------------------
# 14. Figures
# ------------------------------------------------------------

STATE_LEVELS <- c(
  "Low",
  "IEG_high_only",
  "ECM_high_only",
  "IEG_ECM_dual_high"
)

plot_sample <- aHSC_state %>%
  select(
    sample,
    NAS,
    Low_fraction,
    IEG_high_only_fraction,
    ECM_high_only_fraction,
    Dual_high_fraction
  ) %>%
  pivot_longer(
    cols = c(
      Low_fraction,
      IEG_high_only_fraction,
      ECM_high_only_fraction,
      Dual_high_fraction
    ),
    names_to = "state",
    values_to = "fraction"
  ) %>%
  mutate(
    state =
      recode(
        state,
        Low_fraction =
          "Low",
        IEG_high_only_fraction =
          "IEG_high_only",
        ECM_high_only_fraction =
          "ECM_high_only",
        Dual_high_fraction =
          "IEG_ECM_dual_high"
      ),
    state =
      factor(
        state,
        levels = STATE_LEVELS
      ),
    sample =
      factor(
        sample,
        levels =
          aHSC_state %>%
          arrange(
            NAS,
            sample
          ) %>%
          pull(sample)
      )
  )

p1 <- ggplot(
  plot_sample,
  aes(
    x = sample,
    y = fraction,
    fill = state
  )
) +
  geom_col() +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "aHSC functional-state composition",
    subtitle =
      "Sample-matched qHSC p90 thresholds",
    x = "Sample",
    y = "Fraction within aHSC",
    fill = "Functional state"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_aHSC_functional_state_composition.pdf"
  ),
  p1,
  width = 8.5,
  height = 5.5
)

p2 <- ggplot(
  cell_calls,
  aes(
    x = IEG_score_ws,
    y = ECM_score_ws
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.35
  ) +
  facet_wrap(
    ~ current_HSC_state
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "Independent HSC activation axes",
    x =
      "IEG / initiatory score",
    y =
      "ECM / fibrogenic score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_HSC_IEG_vs_ECM_axes.pdf"
  ),
  p2,
  width = 9,
  height = 5
)

# ------------------------------------------------------------
# 15. Save derived object
# ------------------------------------------------------------

saveRDS(
  hsc,
  file.path(
    DIR_OBJ,
    "Human_MASH_HSC_Activation_State_Matrix_v1.4.rds"
  )
)

# ------------------------------------------------------------
# 16. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "identity_annotation_changed",
    "functional_axes_n",
    "primary_threshold_reference",
    "sample_matched_thresholds_n",
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
    "FALSE",
    "2",
    "sample-matched qHSC p90",
    as.character(
      nrow(
        primary_thresholds
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
    "integrity_summary_v1.4.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "source_RDS",
    "identity_annotation",
    "functional_axis_1",
    "functional_axis_2",
    "primary_threshold",
    "sensitivity_threshold",
    "cycling_interpretation",
    "state_interpretation"
  ),
  value = c(
    "v1.4.0",
    HSC_RDS,
    "unchanged: qHSC / aHSC_mesenchymal",
    "IEG/initiatory",
    "ECM/fibrogenic",
    "sample-specific qHSC p90 after within-sample gene z-scoring",
    "pooled qHSC p90 after global gene z-scoring",
    "not included as HSC subtype because direct MKI67/TOP2A evidence is extremely rare",
    "functional state overlay; not replacement subtype annotation"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "HSC_activation_state_provenance_v1.4.tsv"
  )
)

# ------------------------------------------------------------
# 17. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PANEL COVERAGE\n")
cat("============================================================\n")
print(
  panel_summary
)

cat("\n")
cat("============================================================\n")
cat("FUNCTIONAL STATE BY CURRENT IDENTITY\n")
cat("============================================================\n")
print(
  state_by_identity
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY aHSC STATE BY SAMPLE\n")
cat("============================================================\n")
print(
  aHSC_state
)

cat("\n")
cat("============================================================\n")
cat("aHSC FUNCTIONAL STATE vs NAS\n")
cat("============================================================\n")
print(
  nas_assoc
)

cat("\n")
cat("============================================================\n")
cat("aHSC vs qHSC PAIRED TESTS\n")
cat("============================================================\n")
print(
  paired_summary
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY vs SENSITIVITY AGREEMENT\n")
cat("============================================================\n")
print(
  agreement
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
