#!/usr/bin/env Rscript

# ============================================================
# Human MASH TREM2-Mphi -> TGFB1/TGFBR1 -> HSC mechanistic closure
# v1.4.0
#
# PURPOSE
#   1. Apply BH multiple-testing correction to v1.3.1 LR associations
#   2. Quantify sender specificity of TGFB1 across canonical myeloid states
#   3. Build an HSC TGF-beta response target score
#   4. Test whether the target score is enriched in ECM-high HSC
#   5. Correlate TREM2-Mphi TGFB1 with:
#        - HSC TGF-beta response score
#        - HSC ECM total score
#        - ECM-high HSC fraction
#        - TGFBR1/TGFBR2 expression
#   6. Audit the TGFB1/TGFBR1 candidate without reclustering
#
# IMPORTANT
#   - Sample is the unit of association.
#   - No causal inference is claimed from correlation alone.
#   - Canonical HSC receiver remains fixed at 7,640 cells.
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

LR_TABLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Mphi_to_HSC_LR_Audit_v1.3.1",
  "tables",
  "10_candidate_LR_sample_level_associations.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_TGFB1_TGFBR1_Closure_v1.4"
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
  "Human_MASH_TREM2_TGFB1_TGFBR1_Closure_v1.4.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH TREM2/TGFB1/TGFBR1 closure v1.4.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(MYELOID_RDS))
stopifnot(file.exists(HSC_RDS))
stopifnot(file.exists(LR_TABLE))

# ------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------

mphi <- readRDS(MYELOID_RDS)
hsc  <- readRDS(HSC_RDS)

DefaultAssay(mphi) <- "RNA"
DefaultAssay(hsc)  <- "RNA"

cat("Myeloid cells:", ncol(mphi), "\n")
cat("HSC cells:", ncol(hsc), "\n")
cat("HSC features:", nrow(hsc), "\n\n")

if (ncol(hsc) != 7640) {
  warning(
    "Expected 7640 canonical HSC cells; observed ",
    ncol(hsc)
  )
}

# ------------------------------------------------------------
# 2. Metadata requirements
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 3. Normalize RNA
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
# 4. Multiple-testing correction of v1.3.1
# ------------------------------------------------------------

lr <- read_tsv(
  LR_TABLE,
  show_col_types = FALSE
)

required_lr_cols <- c(
  "pair_id",
  "class",
  "ligand",
  "receptor",
  "sender_state",
  "p_ECM_high_fraction",
  "p_NAS"
)

missing_lr <- setdiff(
  required_lr_cols,
  colnames(lr)
)

if (length(missing_lr) > 0) {
  stop(
    "v1.3.1 LR table lacks required columns: ",
    paste(missing_lr, collapse = ";")
  )
}

lr_fdr <- lr %>%
  mutate(
    FDR_ECM_high_fraction =
      p.adjust(
        p_ECM_high_fraction,
        method = "BH"
      ),
    FDR_NAS =
      p.adjust(
        p_NAS,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR_ECM_high_fraction,
    p_ECM_high_fraction
  )

write_tsv(
  lr_fdr,
  file.path(
    DIR_TABLE,
    "01_v1.3.1_LR_associations_with_BH_FDR.tsv"
  )
)

# ------------------------------------------------------------
# 5. TGF-beta response panels
#
# DIRECT_RESPONSE:
#   relatively proximal / canonical transcriptional response.
#
# FIBROGENIC_OUTPUT:
#   downstream extracellular-matrix output.
#
# We keep them separate to avoid treating ECM output itself as
# proof of proximal TGF-beta pathway activation.
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

TGF_RECEPTORS <- c(
  "TGFBR1",
  "TGFBR2"
)

availability <- bind_rows(
  tibble(
    panel = "TGF_DIRECT_RESPONSE",
    gene = TGF_DIRECT_RESPONSE
  ),
  tibble(
    panel = "TGF_FIBROGENIC_OUTPUT",
    gene = TGF_FIBROGENIC_OUTPUT
  ),
  tibble(
    panel = "TGF_RECEPTORS",
    gene = TGF_RECEPTORS
  )
) %>%
  mutate(
    available_in_HSC = gene %in% rownames(hsc)
  )

write_tsv(
  availability,
  file.path(
    DIR_TABLE,
    "02_TGFbeta_feature_availability.tsv"
  )
)

direct_genes <- availability %>%
  filter(
    panel == "TGF_DIRECT_RESPONSE",
    available_in_HSC
  ) %>%
  pull(gene)

fibro_genes <- availability %>%
  filter(
    panel == "TGF_FIBROGENIC_OUTPUT",
    available_in_HSC
  ) %>%
  pull(gene)

receptor_genes <- availability %>%
  filter(
    panel == "TGF_RECEPTORS",
    available_in_HSC
  ) %>%
  pull(gene)

cat("TGF direct-response genes available:\n")
print(direct_genes)
cat("\n")

cat("TGF fibrogenic-output genes available:\n")
print(fibro_genes)
cat("\n")

if (length(direct_genes) < 3) {
  warning(
    "Fewer than 3 direct-response genes are measurable. ",
    "Interpret direct-response score cautiously."
  )
}

if (length(fibro_genes) < 3) {
  stop(
    "Too few fibrogenic-output genes are measurable."
  )
}

if (!all(c("TGFBR1", "TGFBR2") %in% rownames(hsc))) {
  warning(
    "One or more canonical TGF-beta receptors are missing."
  )
}

if (!("TGFB1" %in% rownames(mphi))) {
  stop("TGFB1 is not measurable in the myeloid object.")
}

# ------------------------------------------------------------
# 6. Equal-weight z-score helper
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

  genes <- intersect(
    genes,
    rownames(object)
  )

  if (length(genes) == 0) {
    return(
      rep(
        NA_real_,
        ncol(object)
      )
    )
  }

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
# 7. HSC TGF-beta response scores
# ------------------------------------------------------------

hsc$TGF_direct_response_score_v1.4 <-
  calc_score(
    hsc,
    direct_genes
  )

hsc$TGF_fibrogenic_output_score_v1.4 <-
  calc_score(
    hsc,
    fibro_genes
  )

# NAS0 qHSC reference threshold
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
  warning(
    "NAS0 qHSC reference has fewer than 100 cells; ",
    "falling back to pooled qHSC."
  )

  ref_idx <-
    as.character(
      hmd$HSC_state_v1.2
    ) == "qHSC"
}

TGF_direct_p90 <- as.numeric(
  quantile(
    hsc$TGF_direct_response_score_v1.4[
      ref_idx
    ],
    probs = 0.90,
    na.rm = TRUE
  )
)

hsc$TGF_direct_high_v1.4 <-
  hsc$TGF_direct_response_score_v1.4 >=
  TGF_direct_p90

# ------------------------------------------------------------
# 8. HSC receptor expression
# ------------------------------------------------------------

receptor_data <- FetchData(
  hsc,
  vars = receptor_genes,
  layer = "data"
)

receptor_counts <- FetchData(
  hsc,
  vars = receptor_genes,
  layer = "counts"
)

# ------------------------------------------------------------
# 9. HSC sample-level summary
# ------------------------------------------------------------

hsc_df <- hsc[[]] %>%
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
    ECM_collagen_score = as.numeric(
      ECM_collagen_score_v1.2
    ),
    ECM_remodeling_score = as.numeric(
      ECM_remodeling_score_v1.2
    ),
    TGF_direct_score = as.numeric(
      TGF_direct_response_score_v1.4
    ),
    TGF_fibrogenic_score = as.numeric(
      TGF_fibrogenic_output_score_v1.4
    ),
    TGF_direct_high = as.logical(
      TGF_direct_high_v1.4
    )
  )

for (g in receptor_genes) {

  hsc_df[[paste0(g, "_expr")]] <-
    receptor_data[[g]]

  hsc_df[[paste0(g, "_detected")]] <-
    receptor_counts[[g]] > 0
}

hsc_sample <- hsc_df %>%
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

    ECM_collagen_mean =
      mean(
        ECM_collagen_score,
        na.rm = TRUE
      ),

    ECM_remodeling_mean =
      mean(
        ECM_remodeling_score,
        na.rm = TRUE
      ),

    TGF_direct_mean =
      mean(
        TGF_direct_score,
        na.rm = TRUE
      ),

    TGF_direct_high_fraction =
      mean(
        TGF_direct_high,
        na.rm = TRUE
      ),

    TGF_fibrogenic_mean =
      mean(
        TGF_fibrogenic_score,
        na.rm = TRUE
      ),

    TGFBR1_mean =
      if ("TGFBR1" %in% receptor_genes) {
        mean(
          .data$TGFBR1_expr,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR1_detection =
      if ("TGFBR1" %in% receptor_genes) {
        mean(
          .data$TGFBR1_detected,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR2_mean =
      if ("TGFBR2" %in% receptor_genes) {
        mean(
          .data$TGFBR2_expr,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR2_detection =
      if ("TGFBR2" %in% receptor_genes) {
        mean(
          .data$TGFBR2_detected,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    .groups = "drop"
  ) %>%
  arrange(sample)

write_tsv(
  hsc_sample,
  file.path(
    DIR_TABLE,
    "03_HSC_TGFbeta_ECM_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 10. ECM-high vs ECM-nonhigh HSC:
#     target response + receptors
# ------------------------------------------------------------

hsc_ecm_group <- hsc_df %>%
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

    TGF_direct_mean =
      mean(
        TGF_direct_score,
        na.rm = TRUE
      ),

    TGF_direct_high_fraction =
      mean(
        TGF_direct_high,
        na.rm = TRUE
      ),

    TGF_fibrogenic_mean =
      mean(
        TGF_fibrogenic_score,
        na.rm = TRUE
      ),

    TGFBR1_mean =
      if ("TGFBR1" %in% receptor_genes) {
        mean(
          .data$TGFBR1_expr,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR1_detection =
      if ("TGFBR1" %in% receptor_genes) {
        mean(
          .data$TGFBR1_detected,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR2_mean =
      if ("TGFBR2" %in% receptor_genes) {
        mean(
          .data$TGFBR2_expr,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    TGFBR2_detection =
      if ("TGFBR2" %in% receptor_genes) {
        mean(
          .data$TGFBR2_detected,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    .groups = "drop"
  )

write_tsv(
  hsc_ecm_group,
  file.path(
    DIR_TABLE,
    "04_HSC_TGFbeta_by_ECM_group.tsv"
  )
)

hsc_ecm_wide <- hsc_ecm_group %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = c(
      cells,
      TGF_direct_mean,
      TGF_direct_high_fraction,
      TGF_fibrogenic_mean,
      TGFBR1_mean,
      TGFBR1_detection,
      TGFBR2_mean,
      TGFBR2_detection
    )
  ) %>%
  mutate(
    TGF_direct_delta =
      TGF_direct_mean_ECM_high -
      TGF_direct_mean_ECM_nonhigh,

    TGF_fibrogenic_delta =
      TGF_fibrogenic_mean_ECM_high -
      TGF_fibrogenic_mean_ECM_nonhigh,

    TGFBR1_delta =
      TGFBR1_mean_ECM_high -
      TGFBR1_mean_ECM_nonhigh,

    TGFBR2_delta =
      TGFBR2_mean_ECM_high -
      TGFBR2_mean_ECM_nonhigh
  )

write_tsv(
  hsc_ecm_wide,
  file.path(
    DIR_TABLE,
    "05_HSC_ECMhigh_TGFbeta_deltas_by_sample.tsv"
  )
)

paired_test <- function(
  x_high,
  x_nonhigh
) {

  keep <- complete.cases(
    x_high,
    x_nonhigh
  )

  if (
    sum(keep) < 4 ||
    all(
      x_high[keep] ==
        x_nonhigh[keep]
    )
  ) {
    return(
      tibble(
        n_samples = sum(keep),
        median_delta = median(
          x_high[keep] -
            x_nonhigh[keep],
          na.rm = TRUE
        ),
        p_value = NA_real_
      )
    )
  }

  wt <- suppressWarnings(
    wilcox.test(
      x_high[keep],
      x_nonhigh[keep],
      paired = TRUE,
      exact = FALSE
    )
  )

  tibble(
    n_samples = sum(keep),
    median_delta = median(
      x_high[keep] -
        x_nonhigh[keep],
      na.rm = TRUE
    ),
    p_value = wt$p.value
  )
}

ecm_group_tests <- bind_rows(
  paired_test(
    hsc_ecm_wide$TGF_direct_mean_ECM_high,
    hsc_ecm_wide$TGF_direct_mean_ECM_nonhigh
  ) %>%
    mutate(metric = "TGF_direct_mean"),

  paired_test(
    hsc_ecm_wide$TGF_fibrogenic_mean_ECM_high,
    hsc_ecm_wide$TGF_fibrogenic_mean_ECM_nonhigh
  ) %>%
    mutate(metric = "TGF_fibrogenic_mean"),

  paired_test(
    hsc_ecm_wide$TGFBR1_mean_ECM_high,
    hsc_ecm_wide$TGFBR1_mean_ECM_nonhigh
  ) %>%
    mutate(metric = "TGFBR1_mean"),

  paired_test(
    hsc_ecm_wide$TGFBR2_mean_ECM_high,
    hsc_ecm_wide$TGFBR2_mean_ECM_nonhigh
  ) %>%
    mutate(metric = "TGFBR2_mean")
) %>%
  select(
    metric,
    everything()
  ) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  ecm_group_tests,
  file.path(
    DIR_TABLE,
    "06_ECMhigh_vs_nonhigh_paired_tests.tsv"
  )
)

# ------------------------------------------------------------
# 11. TGFB1 sender expression by sample x myeloid state
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
    TGFB1_detected =
      tgfb1_counts$TGFB1 > 0
  )

tgfb1_by_state <- mphi_df %>%
  group_by(
    sample,
    NAS,
    sender_state
  ) %>%
  summarise(
    cells = n(),
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
    .groups = "drop"
  ) %>%
  arrange(
    sample,
    sender_state
  )

write_tsv(
  tgfb1_by_state,
  file.path(
    DIR_TABLE,
    "07_TGFB1_by_sample_myeloid_state.tsv"
  )
)

# ------------------------------------------------------------
# 12. Sender specificity of TGFB1
# ------------------------------------------------------------

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

trem2_tgfb1 <- tgfb1_by_state %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  rename(
    TREM2_cells = cells,
    TREM2_TGFB1_mean = TGFB1_mean,
    TREM2_TGFB1_median = TGFB1_median,
    TREM2_TGFB1_detection = TGFB1_detection
  ) %>%
  select(
    -sender_state
  )

other_max <- tgfb1_by_state %>%
  filter(
    sender_state != TREM2_STATE
  ) %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    max_other_TGFB1_mean =
      max(
        TGFB1_mean,
        na.rm = TRUE
      ),
    max_other_TGFB1_detection =
      max(
        TGFB1_detection,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

trem2_specificity <- left_join(
  trem2_tgfb1,
  other_max,
  by = c(
    "sample",
    "NAS"
  )
) %>%
  mutate(
    TREM2_minus_max_other_mean =
      TREM2_TGFB1_mean -
      max_other_TGFB1_mean,

    TREM2_minus_max_other_detection =
      TREM2_TGFB1_detection -
      max_other_TGFB1_detection
  )

write_tsv(
  trem2_specificity,
  file.path(
    DIR_TABLE,
    "08_TREM2_TGFB1_sender_specificity_by_sample.tsv"
  )
)

# Pairwise state comparisons against TREM2
other_states <- sort(
  setdiff(
    unique(
      tgfb1_by_state$sender_state
    ),
    TREM2_STATE
  )
)

state_tests <- bind_rows(
  lapply(
    other_states,
    function(st) {

      x <- tgfb1_by_state %>%
        filter(
          sender_state %in%
            c(
              TREM2_STATE,
              st
            )
        ) %>%
        select(
          sample,
          sender_state,
          TGFB1_mean
        ) %>%
        pivot_wider(
          names_from = sender_state,
          values_from = TGFB1_mean
        )

      if (
        !all(
          c(
            TREM2_STATE,
            st
          ) %in%
            colnames(x)
        )
      ) {
        return(
          tibble(
            comparator = st,
            n_samples = 0,
            median_delta_TREM2_minus_other = NA_real_,
            p_value = NA_real_
          )
        )
      }

      a <- x[[TREM2_STATE]]
      b <- x[[st]]

      keep <- complete.cases(
        a,
        b
      )

      if (
        sum(keep) < 4 ||
        all(a[keep] == b[keep])
      ) {
        p <- NA_real_
      } else {
        p <- suppressWarnings(
          wilcox.test(
            a[keep],
            b[keep],
            paired = TRUE,
            exact = FALSE
          )$p.value
        )
      }

      tibble(
        comparator = st,
        n_samples = sum(keep),
        median_delta_TREM2_minus_other =
          median(
            a[keep] -
              b[keep],
            na.rm = TRUE
          ),
        p_value = p
      )
    }
  )
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  state_tests,
  file.path(
    DIR_TABLE,
    "09_TREM2_TGFB1_pairwise_sender_tests.tsv"
  )
)

# ------------------------------------------------------------
# 13. Sample-level TREM2 TGFB1 -> HSC associations
# ------------------------------------------------------------

closure_df <- inner_join(
  trem2_tgfb1,
  hsc_sample,
  by = c(
    "sample",
    "NAS"
  )
) %>%
  arrange(sample)

write_tsv(
  closure_df,
  file.path(
    DIR_TABLE,
    "10_TREM2_TGFB1_HSC_closure_by_sample.tsv"
  )
)

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

association_targets <- c(
  "NAS",
  "TGF_direct_mean",
  "TGF_direct_high_fraction",
  "TGF_fibrogenic_mean",
  "ECM_total_mean",
  "ECM_collagen_mean",
  "ECM_remodeling_mean",
  "ECM_high_fraction",
  "TGFBR1_mean",
  "TGFBR1_detection",
  "TGFBR2_mean",
  "TGFBR2_detection"
)

closure_associations <- bind_rows(
  lapply(
    association_targets,
    function(metric) {

      safe_spearman(
        closure_df$TREM2_TGFB1_mean,
        closure_df[[metric]]
      ) %>%
        mutate(
          sender_metric =
            "TREM2_TGFB1_mean",
          receiver_metric =
            metric,
          .before = 1
        )
    }
  )
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    FDR,
    p_value
  )

write_tsv(
  closure_associations,
  file.path(
    DIR_TABLE,
    "11_TREM2_TGFB1_HSC_sample_level_associations.tsv"
  )
)

# Sensitivity: ligand detection fraction instead of mean expression
closure_detection_associations <- bind_rows(
  lapply(
    association_targets,
    function(metric) {

      safe_spearman(
        closure_df$TREM2_TGFB1_detection,
        closure_df[[metric]]
      ) %>%
        mutate(
          sender_metric =
            "TREM2_TGFB1_detection",
          receiver_metric =
            metric,
          .before = 1
        )
    }
  )
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    FDR,
    p_value
  )

write_tsv(
  closure_detection_associations,
  file.path(
    DIR_TABLE,
    "12_TREM2_TGFB1_detection_HSC_associations.tsv"
  )
)

# ------------------------------------------------------------
# 14. NAS group descriptive summary
# ------------------------------------------------------------

NAS_summary <- closure_df %>%
  group_by(NAS) %>%
  summarise(
    samples = n(),

    TREM2_TGFB1_mean =
      mean(
        TREM2_TGFB1_mean,
        na.rm = TRUE
      ),

    TREM2_TGFB1_detection =
      mean(
        TREM2_TGFB1_detection,
        na.rm = TRUE
      ),

    TGF_direct_mean =
      mean(
        TGF_direct_mean,
        na.rm = TRUE
      ),

    TGF_direct_high_fraction =
      mean(
        TGF_direct_high_fraction,
        na.rm = TRUE
      ),

    TGF_fibrogenic_mean =
      mean(
        TGF_fibrogenic_mean,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high_fraction,
        na.rm = TRUE
      ),

    TGFBR1_mean =
      mean(
        TGFBR1_mean,
        na.rm = TRUE
      ),

    TGFBR2_mean =
      mean(
        TGFBR2_mean,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  NAS_summary,
  file.path(
    DIR_TABLE,
    "13_NAS_mean_TREM2_TGFB1_HSC_closure.tsv"
  )
)

# ------------------------------------------------------------
# 15. Figures
# ------------------------------------------------------------

# TGFB1 by sender state
p1 <- ggplot(
  tgfb1_by_state,
  aes(
    x = sender_state,
    y = TGFB1_mean,
    group = sample
  )
) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "TGFB1 expression across myeloid sender states",
    x = NULL,
    y = "Mean log-normalized TGFB1"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TGFB1_sender_state_by_sample.pdf"
  ),
  p1,
  width = 8,
  height = 6
)

# TREM2 TGFB1 vs ECM-high fraction
p2 <- ggplot(
  closure_df,
  aes(
    x = TREM2_TGFB1_mean,
    y = ECM_high_fraction,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.01,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "TREM2-Mphi TGFB1 vs ECM-high HSC fraction",
    x = "TREM2-Mphi mean TGFB1",
    y = "ECM-high HSC fraction"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_TREM2_TGFB1_vs_ECM_high_fraction.pdf"
  ),
  p2,
  width = 6.5,
  height = 5.5
)

# TREM2 TGFB1 vs direct TGF response
p3 <- ggplot(
  closure_df,
  aes(
    x = TREM2_TGFB1_mean,
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
    title = "TREM2-Mphi TGFB1 vs HSC TGF-beta response",
    x = "TREM2-Mphi mean TGFB1",
    y = "Mean HSC TGF-beta direct-response score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_TREM2_TGFB1_vs_HSC_TGF_direct_response.pdf"
  ),
  p3,
  width = 6.5,
  height = 5.5
)

# ECM-high vs nonhigh direct response
p4 <- ggplot(
  hsc_ecm_group,
  aes(
    x = ECM_group,
    y = TGF_direct_mean,
    group = sample
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "HSC TGF-beta response by ECM functional state",
    x = NULL,
    y = "Mean TGF-beta direct-response score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "04_HSC_TGF_direct_response_ECM_high_vs_nonhigh.pdf"
  ),
  p4,
  width = 6.5,
  height = 5.5
)

# ------------------------------------------------------------
# 16. Save derived object/results
# ------------------------------------------------------------

saveRDS(
  list(
    lr_fdr = lr_fdr,
    availability = availability,
    hsc_sample = hsc_sample,
    hsc_ecm_group = hsc_ecm_group,
    tgfb1_by_state = tgfb1_by_state,
    trem2_specificity = trem2_specificity,
    state_tests = state_tests,
    closure_df = closure_df,
    closure_associations = closure_associations,
    NAS_summary = NAS_summary
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_TREM2_TGFB1_TGFBR1_Closure_v1.4_results.rds"
  )
)

# ------------------------------------------------------------
# 17. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "myeloid_source_RDS",
    "HSC_source_RDS",
    "LR_source_table",
    "HSC_cells",
    "sender_focus",
    "ligand_focus",
    "receiver_focus",
    "functional_receiver",
    "TGF_direct_genes_available",
    "TGF_fibrogenic_genes_available",
    "multiple_testing",
    "association_unit",
    "reclustering",
    "source_overwritten"
  ),
  value = c(
    "v1.4.0",
    MYELOID_RDS,
    HSC_RDS,
    LR_TABLE,
    as.character(ncol(hsc)),
    TREM2_STATE,
    "TGFB1",
    "TGFBR1;TGFBR2",
    "ECM_high_primary_v1.2",
    paste(
      direct_genes,
      collapse = ";"
    ),
    paste(
      fibro_genes,
      collapse = ";"
    ),
    "Benjamini-Hochberg",
    "biological_sample",
    "none",
    "FALSE"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "closure_provenance_v1.4.tsv"
  )
)

# ------------------------------------------------------------
# 18. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "all_8_HSC_samples_present",
    "all_8_myeloid_samples_present",
    "TGFB1_measurable",
    "TGFBR1_measurable",
    "TGFBR2_measurable",
    "direct_response_genes_n",
    "fibrogenic_output_genes_n",
    "v1.3.1_FDR_table_created",
    "closure_sample_table_created",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(ncol(hsc)),
    "7640",
    as.character(ncol(hsc) == 7640),
    as.character(
      length(
        unique(
          hsc_df$sample
        )
      ) == 8
    ),
    as.character(
      length(
        unique(
          mphi_df$sample
        )
      ) == 8
    ),
    as.character(
      "TGFB1" %in%
        rownames(mphi)
    ),
    as.character(
      "TGFBR1" %in%
        rownames(hsc)
    ),
    as.character(
      "TGFBR2" %in%
        rownames(hsc)
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
      nrow(
        lr_fdr
      ) > 0
    ),
    as.character(
      nrow(
        closure_df
      ) > 0
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

# ------------------------------------------------------------
# 19. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TOP v1.3.1 LR RESULTS AFTER BH FDR\n")
cat("============================================================\n")
print(
  lr_fdr %>%
    select(
      pair_id,
      class,
      sender_state,
      rho_ECM_high_fraction,
      p_ECM_high_fraction,
      FDR_ECM_high_fraction,
      rho_NAS,
      p_NAS,
      FDR_NAS
    ) %>%
    head(20)
)

cat("\n")
cat("============================================================\n")
cat("TGF-BETA FEATURE AVAILABILITY\n")
cat("============================================================\n")
print(availability)

cat("\n")
cat("============================================================\n")
cat("TREM2 TGFB1 SENDER SPECIFICITY\n")
cat("============================================================\n")
print(trem2_specificity)

cat("\n")
cat("============================================================\n")
cat("TREM2 TGFB1 PAIRWISE SENDER TESTS\n")
cat("============================================================\n")
print(state_tests)

cat("\n")
cat("============================================================\n")
cat("ECM-HIGH vs ECM-NONHIGH HSC\n")
cat("============================================================\n")
print(ecm_group_tests)

cat("\n")
cat("============================================================\n")
cat("TREM2 TGFB1 -> HSC CLOSURE ASSOCIATIONS\n")
cat("============================================================\n")
print(closure_associations)

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
