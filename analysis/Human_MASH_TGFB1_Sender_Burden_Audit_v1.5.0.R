#!/usr/bin/env Rscript

# ============================================================
# Human MASH TGFB1 sender-source burden audit v1.5.0
#
# PURPOSE
#   1. Re-evaluate which myeloid state contributes most to the
#      sample-level TGFB1 signal.
#   2. Separate:
#        a) state abundance
#        b) within-state TGFB1 expression
#        c) abundance-weighted TGFB1 burden
#   3. Test each sender state against independent HSC endpoints:
#        - TGF-beta direct-response score
#        - TGF-beta direct-high fraction
#        - fibrogenic-output score
#        - ECM-high HSC fraction
#        - ECM total/remodeling/collagen scores
#   4. Apply BH correction across sender-state x endpoint tests.
#
# IMPORTANT
#   - This is a sample-level source-attribution audit.
#   - It does NOT claim causality.
#   - No reclustering.
#   - No source RDS is overwritten.
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

MYELOID_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

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
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5"
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
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH TGFB1 sender burden audit v1.5.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(MYELOID_RDS))
stopifnot(file.exists(HSC_RDS))

# ------------------------------------------------------------
# 1. Load
# ------------------------------------------------------------

mphi <- readRDS(MYELOID_RDS)
hsc  <- readRDS(HSC_RDS)

DefaultAssay(mphi) <- "RNA"
DefaultAssay(hsc)  <- "RNA"

cat("Myeloid cells:", ncol(mphi), "\n")
cat("HSC cells:", ncol(hsc), "\n\n")

if (ncol(hsc) != 7640) {
  warning(
    "Expected 7640 HSC cells, observed ",
    ncol(hsc)
  )
}

required_mphi <- c(
  "sample_v17",
  "NAS_v17",
  "myeloid_state_FINAL_v19"
)

required_hsc <- c(
  "sample",
  "NAS_score",
  "HSC_state_v1.2",
  "ECM_high_primary_v1.2",
  "ECM_total_score_v1.2",
  "ECM_collagen_score_v1.2",
  "ECM_remodeling_score_v1.2"
)

missing_mphi <- setdiff(
  required_mphi,
  colnames(mphi[[]])
)

missing_hsc <- setdiff(
  required_hsc,
  colnames(hsc[[]])
)

if (length(missing_mphi) > 0) {
  stop(
    "Missing myeloid metadata: ",
    paste(missing_mphi, collapse = ";")
  )
}

if (length(missing_hsc) > 0) {
  stop(
    "Missing HSC metadata: ",
    paste(missing_hsc, collapse = ";")
  )
}

if (!("TGFB1" %in% rownames(mphi))) {
  stop("TGFB1 is not measurable in the myeloid object.")
}

# ------------------------------------------------------------
# 2. Normalize
# ------------------------------------------------------------

mphi <- NormalizeData(
  mphi,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 3. Canonical HSC TGF-beta scores
#    Same conceptual definitions as v1.4
# ------------------------------------------------------------

TGF_DIRECT_RESPONSE <- c(
  "SERPINE1",
  "SMAD7",
  "PMEPA1",
  "KLF10",
  "SKIL",
  "TGFBI",
  "THBS1"
)

TGF_FIBROGENIC_OUTPUT <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "FN1",
  "THBS2",
  "POSTN",
  "CTHRC1",
  "LOXL2"
)

direct_genes <- intersect(
  TGF_DIRECT_RESPONSE,
  rownames(hsc)
)

fibro_genes <- intersect(
  TGF_FIBROGENIC_OUTPUT,
  rownames(hsc)
)

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

  genes <- intersect(
    genes,
    rownames(object)
  )

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

hsc$TGF_direct_response_score_v1.5 <-
  calc_score(
    hsc,
    direct_genes
  )

hsc$TGF_fibrogenic_output_score_v1.5 <-
  calc_score(
    hsc,
    fibro_genes
  )

hmd <- hsc[[]]

nas_hsc <- suppressWarnings(
  as.numeric(
    hmd$NAS_score
  )
)

ref_idx <-
  as.character(
    hmd$HSC_state_v1.2
  ) == "qHSC" &
  nas_hsc == 0

if (sum(ref_idx) < 100) {
  ref_idx <-
    as.character(
      hmd$HSC_state_v1.2
    ) == "qHSC"
}

direct_p90 <- as.numeric(
  quantile(
    hsc$TGF_direct_response_score_v1.5[
      ref_idx
    ],
    probs = 0.90,
    na.rm = TRUE
  )
)

hsc$TGF_direct_high_v1.5 <-
  hsc$TGF_direct_response_score_v1.5 >=
  direct_p90

# ------------------------------------------------------------
# 4. HSC sample-level independent endpoints
# ------------------------------------------------------------

hsc_sample <- hsc[[]] %>%
  transmute(
    sample = as.character(sample),
    NAS = suppressWarnings(
      as.numeric(NAS_score)
    ),
    ECM_high = as.logical(
      ECM_high_primary_v1.2
    ),
    ECM_total = as.numeric(
      ECM_total_score_v1.2
    ),
    ECM_collagen = as.numeric(
      ECM_collagen_score_v1.2
    ),
    ECM_remodeling = as.numeric(
      ECM_remodeling_score_v1.2
    ),
    TGF_direct = as.numeric(
      TGF_direct_response_score_v1.5
    ),
    TGF_direct_high = as.logical(
      TGF_direct_high_v1.5
    ),
    TGF_fibrogenic = as.numeric(
      TGF_fibrogenic_output_score_v1.5
    )
  ) %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    HSC_n = n(),

    TGF_direct_mean =
      mean(
        TGF_direct,
        na.rm = TRUE
      ),

    TGF_direct_high_fraction =
      mean(
        TGF_direct_high,
        na.rm = TRUE
      ),

    TGF_fibrogenic_mean =
      mean(
        TGF_fibrogenic,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high,
        na.rm = TRUE
      ),

    ECM_total_mean =
      mean(
        ECM_total,
        na.rm = TRUE
      ),

    ECM_collagen_mean =
      mean(
        ECM_collagen,
        na.rm = TRUE
      ),

    ECM_remodeling_mean =
      mean(
        ECM_remodeling,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(sample)

write_tsv(
  hsc_sample,
  file.path(
    DIR_TABLE,
    "01_HSC_endpoints_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 5. Myeloid TGFB1 cell-level data
# ------------------------------------------------------------

tgfb1_data <- FetchData(
  mphi,
  vars = "TGFB1",
  layer = "data"
)

tgfb1_counts <- FetchData(
  mphi,
  vars = "TGFB1",
  layer = "counts"
)

mphi_df <- mphi[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(
      sample_v17
    ),
    NAS = suppressWarnings(
      as.numeric(NAS_v17)
    ),
    sender_state = as.character(
      myeloid_state_FINAL_v19
    ),
    TGFB1_expr = as.numeric(
      tgfb1_data$TGFB1
    ),
    TGFB1_count = as.numeric(
      tgfb1_counts$TGFB1
    ),
    TGFB1_detected =
      tgfb1_counts$TGFB1 > 0
  )

# ------------------------------------------------------------
# 6. Sample total myeloid counts
# ------------------------------------------------------------

sample_totals <- mphi_df %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    total_myeloid_n = n(),
    total_TGFB1_expr_sum =
      sum(
        TGFB1_expr,
        na.rm = TRUE
      ),
    total_TGFB1_detected_n =
      sum(
        TGFB1_detected,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. State-specific TGFB1 source burden
#
# burden_expr =
#   state fraction x within-state mean normalized TGFB1
#
# This equals:
#   sum(normalized TGFB1 in the state) / total myeloid cells
#
# burden_detection =
#   fraction of ALL myeloid cells that are in this state and
#   TGFB1-positive.
# ------------------------------------------------------------

sender_state <- mphi_df %>%
  group_by(
    sample,
    NAS,
    sender_state
  ) %>%
  summarise(
    state_n = n(),

    TGFB1_mean =
      mean(
        TGFB1_expr,
        na.rm = TRUE
      ),

    TGFB1_median =
      median(
        TGFB1_expr,
        na.rm = TRUE
      ),

    TGFB1_detection =
      mean(
        TGFB1_detected,
        na.rm = TRUE
      ),

    TGFB1_expr_sum =
      sum(
        TGFB1_expr,
        na.rm = TRUE
      ),

    TGFB1_detected_n =
      sum(
        TGFB1_detected,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  left_join(
    sample_totals,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  mutate(
    state_fraction =
      state_n /
      total_myeloid_n,

    TGFB1_burden_expr =
      TGFB1_expr_sum /
      total_myeloid_n,

    TGFB1_burden_detection =
      TGFB1_detected_n /
      total_myeloid_n,

    TGFB1_expression_share =
      ifelse(
        total_TGFB1_expr_sum > 0,
        TGFB1_expr_sum /
          total_TGFB1_expr_sum,
        NA_real_
      ),

    TGFB1_positive_cell_share =
      ifelse(
        total_TGFB1_detected_n > 0,
        TGFB1_detected_n /
          total_TGFB1_detected_n,
        NA_real_
      )
  ) %>%
  arrange(
    sample,
    sender_state
  )

write_tsv(
  sender_state,
  file.path(
    DIR_TABLE,
    "02_TGFB1_sender_state_burden_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 8. Check contribution shares sum to ~1 per sample
# ------------------------------------------------------------

share_audit <- sender_state %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    sum_expression_share =
      sum(
        TGFB1_expression_share,
        na.rm = TRUE
      ),
    sum_positive_cell_share =
      sum(
        TGFB1_positive_cell_share,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_tsv(
  share_audit,
  file.path(
    DIR_TABLE,
    "03_TGFB1_share_integrity_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 9. Total myeloid TGFB1 burden
# ------------------------------------------------------------

total_myeloid <- sample_totals %>%
  mutate(
    total_TGFB1_burden_expr =
      total_TGFB1_expr_sum /
      total_myeloid_n,

    total_TGFB1_detection_fraction =
      total_TGFB1_detected_n /
      total_myeloid_n
  ) %>%
  arrange(sample)

write_tsv(
  total_myeloid,
  file.path(
    DIR_TABLE,
    "04_total_myeloid_TGFB1_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 10. Merge sender metrics with HSC endpoints
# ------------------------------------------------------------

sender_hsc <- sender_state %>%
  inner_join(
    hsc_sample,
    by = c(
      "sample",
      "NAS"
    )
  )

write_tsv(
  sender_hsc,
  file.path(
    DIR_TABLE,
    "05_sender_state_TGFB1_HSC_endpoints.tsv"
  )
)

# ------------------------------------------------------------
# 11. Spearman helper
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
      tibble(
        n_samples = sum(keep),
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
    n_samples = sum(keep),
    rho = unname(
      ct$estimate
    ),
    p_value = ct$p.value
  )
}

# ------------------------------------------------------------
# 12. Sender metric x HSC endpoint association grid
# ------------------------------------------------------------

sender_metrics <- c(
  "state_fraction",
  "TGFB1_mean",
  "TGFB1_detection",
  "TGFB1_burden_expr",
  "TGFB1_burden_detection",
  "TGFB1_expression_share",
  "TGFB1_positive_cell_share"
)

receiver_metrics <- c(
  "TGF_direct_mean",
  "TGF_direct_high_fraction",
  "TGF_fibrogenic_mean",
  "ECM_high_fraction",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "ECM_remodeling_mean"
)

states <- sort(
  unique(
    sender_hsc$sender_state
  )
)

association_list <- list()
k <- 1L

for (st in states) {

  sub <- sender_hsc %>%
    filter(
      sender_state == st
    )

  for (sm in sender_metrics) {

    for (rm in receiver_metrics) {

      tmp <- safe_spearman(
        sub[[sm]],
        sub[[rm]]
      )

      association_list[[k]] <- tmp %>%
        mutate(
          sender_state = st,
          sender_metric = sm,
          receiver_metric = rm,
          .before = 1
        )

      k <- k + 1L
    }
  }
}

assoc <- bind_rows(
  association_list
) %>%
  mutate(
    FDR_global =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  group_by(
    sender_metric
  ) %>%
  mutate(
    FDR_within_sender_metric =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  arrange(
    FDR_global,
    p_value
  )

write_tsv(
  assoc,
  file.path(
    DIR_TABLE,
    "06_sender_TGFB1_burden_HSC_associations.tsv"
  )
)

# ------------------------------------------------------------
# 13. Focused primary source-attribution family
#
# Primary sender metric:
#   TGFB1_burden_expr
#
# Primary receiver endpoints:
#   TGF_direct_mean
#   TGF_fibrogenic_mean
#   ECM_high_fraction
#
# BH correction is performed only across this prespecified
# 5-state x 3-endpoint family.
# ------------------------------------------------------------

primary_assoc <- assoc %>%
  filter(
    sender_metric ==
      "TGFB1_burden_expr",
    receiver_metric %in%
      c(
        "TGF_direct_mean",
        "TGF_fibrogenic_mean",
        "ECM_high_fraction"
      )
  ) %>%
  mutate(
    FDR_primary_family =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR_primary_family,
    p_value
  )

write_tsv(
  primary_assoc,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_TGFB1_burden_source_attribution.tsv"
  )
)

# ------------------------------------------------------------
# 14. NAS association of sender burden
# ------------------------------------------------------------

sender_NAS <- bind_rows(
  lapply(
    states,
    function(st) {

      sub <- sender_state %>%
        filter(
          sender_state == st
        )

      bind_rows(
        lapply(
          c(
            "state_fraction",
            "TGFB1_mean",
            "TGFB1_burden_expr",
            "TGFB1_expression_share"
          ),
          function(metric) {

            safe_spearman(
              sub[[metric]],
              sub$NAS
            ) %>%
              mutate(
                sender_state = st,
                metric = metric,
                .before = 1
              )
          }
        )
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
  sender_NAS,
  file.path(
    DIR_TABLE,
    "08_sender_state_TGFB1_NAS_associations.tsv"
  )
)

# ------------------------------------------------------------
# 15. NAS-level descriptive source composition
# ------------------------------------------------------------

NAS_sender <- sender_state %>%
  group_by(
    NAS,
    sender_state
  ) %>%
  summarise(
    samples = n(),

    mean_state_fraction =
      mean(
        state_fraction,
        na.rm = TRUE
      ),

    mean_TGFB1_mean =
      mean(
        TGFB1_mean,
        na.rm = TRUE
      ),

    mean_TGFB1_burden_expr =
      mean(
        TGFB1_burden_expr,
        na.rm = TRUE
      ),

    mean_TGFB1_expression_share =
      mean(
        TGFB1_expression_share,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  NAS_sender,
  file.path(
    DIR_TABLE,
    "09_NAS_mean_TGFB1_sender_burden.tsv"
  )
)

# ------------------------------------------------------------
# 16. Rank sender states by TGFB1 contribution per sample
# ------------------------------------------------------------

ranked_sender <- sender_state %>%
  group_by(
    sample,
    NAS
  ) %>%
  arrange(
    desc(TGFB1_burden_expr),
    .by_group = TRUE
  ) %>%
  mutate(
    burden_rank = row_number()
  ) %>%
  ungroup()

write_tsv(
  ranked_sender,
  file.path(
    DIR_TABLE,
    "10_sender_TGFB1_burden_rank_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 17. Figures
# ------------------------------------------------------------

# TGFB1 burden by sender state and sample
p1 <- ggplot(
  sender_state,
  aes(
    x = factor(
      sample,
      levels = paste0(
        "LT",
        1:8
      )
    ),
    y = TGFB1_burden_expr,
    group = sender_state,
    linetype = sender_state
  )
) +
  geom_line() +
  geom_point(size = 2.5) +
  theme_classic(base_size = 11) +
  labs(
    title = "Abundance-weighted TGFB1 burden by myeloid sender state",
    x = NULL,
    y = "TGFB1 burden index",
    linetype = "Sender state"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TGFB1_burden_by_sender_state_sample.pdf"
  ),
  p1,
  width = 9,
  height = 6
)

# TGFB1 expression-share composition
p2 <- ggplot(
  sender_state,
  aes(
    x = factor(
      sample,
      levels = paste0(
        "LT",
        1:8
      )
    ),
    y = TGFB1_expression_share,
    fill = sender_state
  )
) +
  geom_col() +
  theme_classic(base_size = 11) +
  labs(
    title = "Myeloid-state contribution to total TGFB1 expression",
    x = NULL,
    y = "Fraction of myeloid TGFB1 signal",
    fill = "Sender state"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_TGFB1_expression_share_by_sample.pdf"
  ),
  p2,
  width = 9,
  height = 6
)

# Primary associations
p3_df <- primary_assoc %>%
  filter(
    !is.na(rho)
  ) %>%
  mutate(
    label = paste(
      sender_state,
      receiver_metric,
      sep = " | "
    ),
    label = factor(
      label,
      levels = rev(
        unique(label)
      )
    )
  )

if (nrow(p3_df) > 0) {

  p3 <- ggplot(
    p3_df,
    aes(
      x = label,
      y = rho
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 10) +
    labs(
      title = "Primary TGFB1 sender-burden associations",
      x = NULL,
      y = "Spearman rho"
    )

  ggsave(
    file.path(
      DIR_FIG,
      "03_PRIMARY_TGFB1_burden_associations.pdf"
    ),
    p3,
    width = 9,
    height = 7
  )
}

# ------------------------------------------------------------
# 18. Save results object
# ------------------------------------------------------------

saveRDS(
  list(
    hsc_sample = hsc_sample,
    sender_state = sender_state,
    total_myeloid = total_myeloid,
    assoc = assoc,
    primary_assoc = primary_assoc,
    sender_NAS = sender_NAS,
    NAS_sender = NAS_sender,
    ranked_sender = ranked_sender
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5_results.rds"
  )
)

# ------------------------------------------------------------
# 19. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "myeloid_source_RDS",
    "HSC_source_RDS",
    "HSC_cells",
    "TGFB1_sender_metric_primary",
    "TGFB1_burden_expr_definition",
    "primary_receiver_endpoints",
    "multiple_testing_primary",
    "association_unit",
    "TGF_direct_genes_available",
    "TGF_fibrogenic_genes_available",
    "reclustering",
    "source_overwritten"
  ),
  value = c(
    "v1.5.0",
    MYELOID_RDS,
    HSC_RDS,
    as.character(
      ncol(hsc)
    ),
    "TGFB1_burden_expr",
    "sum_within_state_log_normalized_TGFB1_divided_by_total_myeloid_cells",
    "TGF_direct_mean;TGF_fibrogenic_mean;ECM_high_fraction",
    "BH_across_5_states_x_3_endpoints",
    "biological_sample",
    paste(
      direct_genes,
      collapse = ";"
    ),
    paste(
      fibro_genes,
      collapse = ";"
    ),
    "none",
    "FALSE"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "sender_burden_provenance_v1.5.tsv"
  )
)

# ------------------------------------------------------------
# 20. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "all_8_HSC_samples_present",
    "all_8_myeloid_samples_present",
    "TGFB1_measurable",
    "sender_states_n",
    "direct_response_genes_n",
    "fibrogenic_output_genes_n",
    "expression_share_sums_to_1",
    "positive_cell_share_sums_to_1",
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
        unique(
          hsc_sample$sample
        )
      ) == 8
    ),
    as.character(
      length(
        unique(
          sample_totals$sample
        )
      ) == 8
    ),
    as.character(
      "TGFB1" %in%
        rownames(mphi)
    ),
    as.character(
      length(states)
    ),
    as.character(
      length(
        direct_genes
      )
    ),
    as.character(
      length(
        fibro_genes
      )
    ),
    as.character(
      all(
        abs(
          share_audit$sum_expression_share -
            1
        ) < 1e-8 |
          is.na(
            share_audit$sum_expression_share
          )
      )
    ),
    as.character(
      all(
        abs(
          share_audit$sum_positive_cell_share -
            1
        ) < 1e-8 |
          is.na(
            share_audit$sum_positive_cell_share
          )
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
    "integrity_summary_v1.5.tsv"
  )
)

# ------------------------------------------------------------
# 21. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TGFB1 SENDER BURDEN BY SAMPLE\n")
cat("============================================================\n")
print(
  sender_state %>%
    select(
      sample,
      NAS,
      sender_state,
      state_n,
      state_fraction,
      TGFB1_mean,
      TGFB1_detection,
      TGFB1_burden_expr,
      TGFB1_expression_share
    )
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY TGFB1 SOURCE ATTRIBUTION\n")
cat("============================================================\n")
print(primary_assoc)

cat("\n")
cat("============================================================\n")
cat("TOP GLOBAL ASSOCIATIONS\n")
cat("============================================================\n")
print(
  assoc %>%
    head(30)
)

cat("\n")
cat("============================================================\n")
cat("TGFB1 SENDER BURDEN vs NAS\n")
cat("============================================================\n")
print(sender_NAS)

cat("\n")
cat("============================================================\n")
cat("NAS MEAN SENDER BURDEN\n")
cat("============================================================\n")
print(NAS_sender)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
