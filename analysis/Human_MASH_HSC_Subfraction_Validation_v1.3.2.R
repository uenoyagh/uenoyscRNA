#!/usr/bin/env Rscript

# ============================================================
# Human MASH HSC subfraction validation v1.3.2
#
# PURPOSE
#   Validate the current qHSC / aHSC_mesenchymal annotation
#   without reclustering and without forcing every cell into
#   a finer subtype.
#
# DESIGN PRINCIPLE
#   Identity axis:
#     Quiescent <---- mixed/intermediate ----> Myofibroblastic
#
#   Functional overlays (NOT identity subtypes):
#     - Initiatory / IEG-like program
#     - Cycling evidence
#     - Pericyte/VSMC sensitivity
#     - Portal-fibroblast sensitivity
#
# WHY v1.3.2
#   v1.3.1 showed directional concordance but:
#     1) its discordance code used obsolete label names
#     2) dominant-program assignment forced all cells into one
#        of four programs
#     3) Cycling coverage was only 3/6 genes
#     4) Inflammatory coverage lacked CCL2/IL6/CXCL8
#
#   Therefore v1.3.2 validates the existing 2-state identity
#   and treats finer programs as overlays only.
#
# No reclustering.
# No source RDS overwrite.
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
  "Human_MASH_HSC_Subfraction_Validation_v1.3.2"
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
  "Human_MASH_HSC_Subfraction_Validation_v1.3.2.log"
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
cat("Human MASH HSC subfraction validation v1.3.2\n")
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
    "Missing metadata: ",
    paste(
      missing_meta,
      collapse = ";"
    )
  )
}

current_labels <- sort(
  unique(
    as.character(
      hsc$HSC_state_v1.2
    )
  )
)

expected_labels <- c(
  "aHSC_mesenchymal",
  "qHSC"
)

if (!setequal(
  current_labels,
  expected_labels
)) {
  stop(
    "Unexpected current HSC labels: ",
    paste(
      current_labels,
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
# 1. Marker panels
# ------------------------------------------------------------

PANELS <- list(

  Quiescent = c(
    "RBP1",
    "LRAT",
    "RELN",
    "CYGB",
    "NGFR",
    "SYNM",
    "LHX2",
    "ECM1",
    "COLEC10",
    "COLEC11"
  ),

  Myofibroblastic = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL5A1",
    "COL5A2",
    "ACTA2",
    "TAGLN",
    "CTHRC1",
    "POSTN",
    "THBS2",
    "LOXL2",
    "TGFBI",
    "TIMP1"
  ),

  Initiatory_IEG_overlay = c(
    "CCL2",
    "THBS1",
    "TIMP1",
    "JUNB",
    "FOS",
    "FOSB",
    "IL6",
    "CXCL8"
  ),

  Cycling_overlay = c(
    "MKI67",
    "TOP2A",
    "CENPF",
    "BIRC5",
    "UBE2C",
    "TYMS"
  ),

  Pericyte_VSMC_sensitivity = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "NOTCH3",
    "MYH11",
    "MYL9",
    "PDGFRB",
    "DES"
  ),

  PortalFibroblast_sensitivity = c(
    "PI16",
    "COL15A1",
    "FBLN1",
    "DPT",
    "SFRP2",
    "DCN",
    "LUM"
  )
)

INTERPRETATION <- tibble(
  panel = names(PANELS),
  role = c(
    "PRIMARY_IDENTITY",
    "PRIMARY_IDENTITY",
    "FUNCTIONAL_OVERLAY_ONLY",
    "FUNCTIONAL_OVERLAY_ONLY",
    "SENSITIVITY_ONLY",
    "SENSITIVITY_ONLY"
  ),
  interpretation = c(
    "Supports qHSC identity; LRAT absence in Flex panel limits completeness",
    "Strong activated/myofibroblastic identity support",
    "Do not call inflammatory HSC subtype if CCL2/IL6/CXCL8 are unavailable",
    "Do not call cycling subtype from dominant score alone; require direct cycling-gene evidence",
    "Potential pericyte/VSMC-like admixture only; overlapping mesenchymal genes prevent automatic exclusion",
    "Potential portal-fibroblast-like program only; PI16/SFRP2 absence limits subtype assignment"
  )
)

# ------------------------------------------------------------
# 2. Marker availability
# ------------------------------------------------------------

marker_audit <- bind_rows(
  lapply(
    names(PANELS),
    function(.panel_name) {

      genes <- PANELS[[.panel_name]]

      tibble(
        panel = rep(
          .panel_name,
          length(genes)
        ),
        gene = genes,
        available =
          genes %in%
          rownames(hsc)
      )
    }
  )
)

coverage <- marker_audit %>%
  group_by(panel) %>%
  summarise(
    requested_n = n(),
    available_n = sum(available),
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
  ) %>%
  left_join(
    INTERPRETATION,
    by = "panel"
  )

write_tsv(
  marker_audit,
  file.path(
    DIR_TABLE,
    "01_marker_availability.tsv"
  )
)

write_tsv(
  coverage,
  file.path(
    DIR_TABLE,
    "02_panel_coverage_and_interpretation.tsv"
  )
)

# ------------------------------------------------------------
# 3. Equal-weight z-score modules
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

score_panel <- function(
  object,
  genes
) {

  g <- intersect(
    genes,
    rownames(object)
  )

  if (length(g) < 2) {
    return(
      rep(
        NA_real_,
        ncol(object)
      )
    )
  }

  x <- FetchData(
    object,
    vars = g,
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

for (p in names(PANELS)) {

  nm <- paste0(
    "HSCv132_",
    p
  )

  hsc[[nm]] <- score_panel(
    hsc,
    PANELS[[p]]
  )
}

# ------------------------------------------------------------
# 4. Identity support axis
#
# delta_identity =
#   Myofibroblastic score - Quiescent score
#
# At margin 0.30:
#   <= -0.30 : Quiescent_supported
#   >= +0.30 : Activated_supported
#   otherwise: Mixed_intermediate
#
# Threshold sensitivity is reported separately.
# ------------------------------------------------------------

md <- hsc[[]]

identity_delta <-
  md$HSCv132_Myofibroblastic -
  md$HSCv132_Quiescent

classify_margin <- function(
  delta,
  margin
) {

  case_when(
    delta <= -margin ~
      "Quiescent_supported",

    delta >= margin ~
      "Activated_supported",

    TRUE ~
      "Mixed_intermediate"
  )
}

hsc$HSC_identity_delta_v1.3.2 <-
  identity_delta

hsc$HSC_identity_support_v1.3.2 <-
  classify_margin(
    identity_delta,
    0.30
  )

md <- hsc[[]]

# ------------------------------------------------------------
# 5. Correct current-label concordance
# ------------------------------------------------------------

cell_audit <- md %>%
  rownames_to_column(
    "cell"
  ) %>%
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
      ),

    Quiescent_score =
      HSCv132_Quiescent,

    Myofibroblastic_score =
      HSCv132_Myofibroblastic,

    identity_delta =
      HSC_identity_delta_v1.3.2,

    identity_support =
      as.character(
        HSC_identity_support_v1.3.2
      ),

    Initiatory_IEG_score =
      HSCv132_Initiatory_IEG_overlay,

    Cycling_score =
      HSCv132_Cycling_overlay,

    Pericyte_VSMC_score =
      HSCv132_Pericyte_VSMC_sensitivity,

    PortalFibroblast_score =
      HSCv132_PortalFibroblast_sensitivity
  ) %>%
  mutate(
    high_confidence_concordant =
      case_when(
        current_HSC_state ==
          "qHSC" &
          identity_support ==
            "Quiescent_supported" ~ TRUE,

        current_HSC_state ==
          "aHSC_mesenchymal" &
          identity_support ==
            "Activated_supported" ~ TRUE,

        TRUE ~ FALSE
      ),

    high_confidence_discordant =
      case_when(
        current_HSC_state ==
          "qHSC" &
          identity_support ==
            "Activated_supported" ~ TRUE,

        current_HSC_state ==
          "aHSC_mesenchymal" &
          identity_support ==
            "Quiescent_supported" ~ TRUE,

        TRUE ~ FALSE
      )
  )

write_tsv(
  cell_audit,
  file.path(
    DIR_TABLE,
    "03_cell_level_identity_audit.tsv.gz"
  )
)

concordance_summary <- cell_audit %>%
  group_by(
    current_HSC_state
  ) %>%
  summarise(
    cells = n(),

    concordant_n =
      sum(
        high_confidence_concordant
      ),

    concordant_fraction =
      mean(
        high_confidence_concordant
      ),

    mixed_n =
      sum(
        identity_support ==
          "Mixed_intermediate"
      ),

    mixed_fraction =
      mean(
        identity_support ==
          "Mixed_intermediate"
      ),

    discordant_n =
      sum(
        high_confidence_discordant
      ),

    discordant_fraction =
      mean(
        high_confidence_discordant
      ),

    .groups = "drop"
  )

write_tsv(
  concordance_summary,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_current_label_concordance.tsv"
  )
)

cross_tab <- cell_audit %>%
  count(
    current_HSC_state,
    identity_support,
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
  cross_tab,
  file.path(
    DIR_TABLE,
    "05_current_label_vs_identity_support.tsv"
  )
)

# ------------------------------------------------------------
# 6. Margin sensitivity
# ------------------------------------------------------------

MARGINS <- c(
  0,
  0.10,
  0.25,
  0.30,
  0.50
)

margin_sensitivity <- bind_rows(
  lapply(
    MARGINS,
    function(margin) {

      cls <- classify_margin(
        identity_delta,
        margin
      )

      tibble(
        cell = colnames(hsc),
        current_HSC_state =
          as.character(
            hsc$HSC_state_v1.2
          ),
        identity_support = cls
      ) %>%
        group_by(
          current_HSC_state
        ) %>%
        summarise(
          cells = n(),

          current_label_supported_fraction =
            ifelse(
              unique(
                current_HSC_state
              ) ==
                "qHSC",

              mean(
                identity_support ==
                  "Quiescent_supported"
              ),

              mean(
                identity_support ==
                  "Activated_supported"
              )
            ),

          opposite_label_fraction =
            ifelse(
              unique(
                current_HSC_state
              ) ==
                "qHSC",

              mean(
                identity_support ==
                  "Activated_supported"
              ),

              mean(
                identity_support ==
                  "Quiescent_supported"
              )
            ),

          mixed_fraction =
            mean(
              identity_support ==
                "Mixed_intermediate"
            ),

          .groups = "drop"
        ) %>%
        mutate(
          margin = margin,
          .before = 1
        )
    }
  )
)

write_tsv(
  margin_sensitivity,
  file.path(
    DIR_TABLE,
    "06_identity_margin_sensitivity.tsv"
  )
)

# ------------------------------------------------------------
# 7. Sample-level current-label validation
#
# Biological sample is the unit of inference.
#
# Expected:
#   qHSC has higher Quiescent score than aHSC
#   aHSC has higher Myofibroblastic score than qHSC
# ------------------------------------------------------------

sample_label_scores <- cell_audit %>%
  group_by(
    sample,
    NAS,
    current_HSC_state
  ) %>%
  summarise(
    cells = n(),

    mean_Quiescent =
      mean(
        Quiescent_score,
        na.rm = TRUE
      ),

    mean_Myofibroblastic =
      mean(
        Myofibroblastic_score,
        na.rm = TRUE
      ),

    mean_identity_delta =
      mean(
        identity_delta,
        na.rm = TRUE
      ),

    mean_IEG =
      mean(
        Initiatory_IEG_score,
        na.rm = TRUE
      ),

    mean_Cycling =
      mean(
        Cycling_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  sample_label_scores,
  file.path(
    DIR_TABLE,
    "07_sample_by_current_label_program_scores.tsv"
  )
)

wide <- sample_label_scores %>%
  select(
    sample,
    NAS,
    current_HSC_state,
    mean_Quiescent,
    mean_Myofibroblastic,
    mean_identity_delta,
    mean_IEG,
    mean_Cycling
  ) %>%
  pivot_wider(
    names_from =
      current_HSC_state,

    values_from = c(
      mean_Quiescent,
      mean_Myofibroblastic,
      mean_identity_delta,
      mean_IEG,
      mean_Cycling
    )
  )

paired_deltas <- wide %>%
  transmute(
    sample,
    NAS,

    q_marker_delta_q_minus_a =
      mean_Quiescent_qHSC -
      mean_Quiescent_aHSC_mesenchymal,

    myofib_marker_delta_a_minus_q =
      mean_Myofibroblastic_aHSC_mesenchymal -
      mean_Myofibroblastic_qHSC,

    identity_delta_a_minus_q =
      mean_identity_delta_aHSC_mesenchymal -
      mean_identity_delta_qHSC,

    IEG_delta_a_minus_q =
      mean_IEG_aHSC_mesenchymal -
      mean_IEG_qHSC,

    Cycling_delta_a_minus_q =
      mean_Cycling_aHSC_mesenchymal -
      mean_Cycling_qHSC
  )

write_tsv(
  paired_deltas,
  file.path(
    DIR_TABLE,
    "08_PRIMARY_sample_paired_label_deltas.tsv"
  )
)

paired_test <- function(
  x,
  expected_direction
) {

  keep <- is.finite(x)
  x <- x[keep]

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
      length(x),

    median_delta =
      median(
        x,
        na.rm = TRUE
      ),

    positive_samples =
      sum(
        x > 0
      ),

    negative_samples =
      sum(
        x < 0
      ),

    p_value =
      wt$p.value,

    expected_direction =
      expected_direction
  )
}

paired_test_summary <- bind_rows(

  paired_test(
    paired_deltas$q_marker_delta_q_minus_a,
    "positive"
  ) %>%
    mutate(
      metric =
        "Quiescent_score_qHSC_minus_aHSC",
      .before = 1
    ),

  paired_test(
    paired_deltas$myofib_marker_delta_a_minus_q,
    "positive"
  ) %>%
    mutate(
      metric =
        "Myofibroblastic_score_aHSC_minus_qHSC",
      .before = 1
    ),

  paired_test(
    paired_deltas$identity_delta_a_minus_q,
    "positive"
  ) %>%
    mutate(
      metric =
        "Identity_delta_aHSC_minus_qHSC",
      .before = 1
    ),

  paired_test(
    paired_deltas$IEG_delta_a_minus_q,
    "descriptive_only"
  ) %>%
    mutate(
      metric =
        "IEG_overlay_aHSC_minus_qHSC",
      .before = 1
    ),

  paired_test(
    paired_deltas$Cycling_delta_a_minus_q,
    "descriptive_only"
  ) %>%
    mutate(
      metric =
        "Cycling_overlay_aHSC_minus_qHSC",
      .before = 1
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
  paired_test_summary,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_sample_level_label_validation.tsv"
  )
)

# ------------------------------------------------------------
# 8. Direct cycling-gene evidence
#
# Do not infer cycling subtype from module dominance.
# Report direct detection of the available core markers.
# ------------------------------------------------------------

cycling_available <- intersect(
  PANELS$Cycling_overlay,
  rownames(hsc)
)

cycling_count <- FetchData(
  hsc,
  vars = cycling_available,
  layer = "counts"
)

cycling_detected_n <- rowSums(
  cycling_count > 0
)

cycling_core <- intersect(
  c(
    "MKI67",
    "TOP2A"
  ),
  colnames(
    cycling_count
  )
)

if (length(cycling_core) > 0) {

  cycling_core_any <- rowSums(
    cycling_count[
      ,
      cycling_core,
      drop = FALSE
    ] > 0
  ) > 0

} else {

  cycling_core_any <- rep(
    FALSE,
    nrow(
      cycling_count
    )
  )
}

cycling_evidence <- cell_audit %>%
  select(
    cell,
    sample,
    NAS,
    current_HSC_state
  ) %>%
  mutate(
    cycling_genes_available_n =
      length(
        cycling_available
      ),

    cycling_detected_genes_n =
      cycling_detected_n,

    MKI67_or_TOP2A_detected =
      cycling_core_any,

    two_or_more_cycling_genes_detected =
      cycling_detected_n >= 2
  )

write_tsv(
  cycling_evidence,
  file.path(
    DIR_TABLE,
    "10_cell_level_cycling_evidence.tsv.gz"
  )
)

cycling_summary <- cycling_evidence %>%
  group_by(
    current_HSC_state
  ) %>%
  summarise(
    cells = n(),

    MKI67_or_TOP2A_positive_n =
      sum(
        MKI67_or_TOP2A_detected
      ),

    MKI67_or_TOP2A_positive_fraction =
      mean(
        MKI67_or_TOP2A_detected
      ),

    two_or_more_cycling_genes_n =
      sum(
        two_or_more_cycling_genes_detected
      ),

    two_or_more_cycling_genes_fraction =
      mean(
        two_or_more_cycling_genes_detected
      ),

    .groups = "drop"
  )

write_tsv(
  cycling_summary,
  file.path(
    DIR_TABLE,
    "11_cycling_evidence_summary.tsv"
  )
)

# ------------------------------------------------------------
# 9. Sample-level identity support + NAS correlation
# ------------------------------------------------------------

sample_identity <- cell_audit %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    HSC_n = n(),

    current_qHSC_fraction =
      mean(
        current_HSC_state ==
          "qHSC"
      ),

    current_aHSC_fraction =
      mean(
        current_HSC_state ==
          "aHSC_mesenchymal"
      ),

    Quiescent_supported_fraction =
      mean(
        identity_support ==
          "Quiescent_supported"
      ),

    Activated_supported_fraction =
      mean(
        identity_support ==
          "Activated_supported"
      ),

    Mixed_intermediate_fraction =
      mean(
        identity_support ==
          "Mixed_intermediate"
      ),

    mean_Quiescent_score =
      mean(
        Quiescent_score,
        na.rm = TRUE
      ),

    mean_Myofibroblastic_score =
      mean(
        Myofibroblastic_score,
        na.rm = TRUE
      ),

    mean_identity_delta =
      mean(
        identity_delta,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  sample_identity,
  file.path(
    DIR_TABLE,
    "12_sample_level_identity_summary.tsv"
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

  ct <- suppressWarnings(
    cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n =
      sum(keep),

    rho =
      unname(
        ct$estimate
      ),

    p_value =
      ct$p.value
  )
}

nas_assoc <- bind_rows(

  safe_spearman(
    sample_identity$NAS,
    sample_identity$current_qHSC_fraction
  ) %>%
    mutate(
      metric =
        "current_qHSC_fraction",
      .before = 1
    ),

  safe_spearman(
    sample_identity$NAS,
    sample_identity$current_aHSC_fraction
  ) %>%
    mutate(
      metric =
        "current_aHSC_fraction",
      .before = 1
    ),

  safe_spearman(
    sample_identity$NAS,
    sample_identity$Quiescent_supported_fraction
  ) %>%
    mutate(
      metric =
        "Quiescent_supported_fraction",
      .before = 1
    ),

  safe_spearman(
    sample_identity$NAS,
    sample_identity$Activated_supported_fraction
  ) %>%
    mutate(
      metric =
        "Activated_supported_fraction",
      .before = 1
    ),

  safe_spearman(
    sample_identity$NAS,
    sample_identity$mean_identity_delta
  ) %>%
    mutate(
      metric =
        "mean_identity_delta",
      .before = 1
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
  nas_assoc,
  file.path(
    DIR_TABLE,
    "13_sample_level_NAS_association.tsv"
  )
)

# ------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------

plot_df <- cell_audit %>%
  select(
    current_HSC_state,
    Quiescent_score,
    Myofibroblastic_score,
    identity_delta,
    Initiatory_IEG_score,
    Cycling_score
  ) %>%
  pivot_longer(
    -current_HSC_state,
    names_to = "program",
    values_to = "score"
  )

p1 <- ggplot(
  plot_df,
  aes(
    x = current_HSC_state,
    y = score
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  facet_wrap(
    ~ program,
    scales = "free_y",
    ncol = 3
  ) +
  theme_classic(
    base_size = 10
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
      "Human MASH HSC identity validation",
    x = NULL,
    y = "Score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_HSC_identity_validation_violin.pdf"
  ),
  p1,
  width = 11,
  height = 7
)

p2 <- ggplot(
  cell_audit,
  aes(
    x = current_HSC_state,
    y = identity_delta
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA
  ) +
  geom_hline(
    yintercept = c(
      -0.30,
      0.30
    ),
    linetype = 2
  ) +
  theme_classic(
    base_size = 11
  ) +
  labs(
    title =
      "HSC identity axis",
    subtitle =
      "Myofibroblastic score - Quiescent score",
    x = NULL,
    y =
      "Identity delta"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_HSC_identity_delta_by_current_label.pdf"
  ),
  p2,
  width = 7,
  height = 5.5
)

# ------------------------------------------------------------
# 11. Save derived audit object
# ------------------------------------------------------------

saveRDS(
  hsc,
  file.path(
    DIR_OBJ,
    "Human_MASH_HSC_Subfraction_Validation_v1.3.2.rds"
  )
)

# ------------------------------------------------------------
# 12. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "current_labels",
    "identity_axis_used",
    "forced_multiclass_assignment",
    "cycling_called_by_dominant_module",
    "reclustering_performed",
    "source_RDS_overwritten",
    "derived_RDS_saved"
  ),
  value = c(
    as.character(
      ncol(hsc)
    ),
    "7640",
    as.character(
      ncol(hsc) == 7640
    ),
    paste(
      current_labels,
      collapse = ";"
    ),
    "Quiescent_vs_Myofibroblastic",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "TRUE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.3.2.tsv"
  )
)

provenance <- tibble(
  item = c(
    "analysis_version",
    "source_RDS",
    "current_annotation",
    "primary_identity_validation",
    "identity_support_margin",
    "IEG_interpretation",
    "Cycling_interpretation",
    "Pericyte_VSMC_interpretation",
    "PortalFibroblast_interpretation",
    "ECM_high_interpretation"
  ),
  value = c(
    "v1.3.2",
    HSC_RDS,
    "HSC_state_v1.2: qHSC / aHSC_mesenchymal",
    "Quiescent score versus Myofibroblastic score",
    "0.30; sensitivity margins also reported",
    "functional overlay only; insufficient canonical inflammatory genes for formal subtype",
    "functional overlay only; direct MKI67/TOP2A detection reported separately",
    "sensitivity only",
    "sensitivity only",
    "functional state, not HSC identity subtype"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "HSC_subfraction_validation_provenance_v1.3.2.tsv"
  )
)

# ------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PANEL COVERAGE / INTERPRETATION\n")
cat("============================================================\n")
print(coverage)

cat("\n")
cat("============================================================\n")
cat("PRIMARY CURRENT-LABEL CONCORDANCE\n")
cat("============================================================\n")
print(concordance_summary)

cat("\n")
cat("============================================================\n")
cat("CURRENT LABEL vs IDENTITY SUPPORT\n")
cat("============================================================\n")
print(cross_tab)

cat("\n")
cat("============================================================\n")
cat("SAMPLE-LEVEL LABEL VALIDATION\n")
cat("============================================================\n")
print(paired_test_summary)

cat("\n")
cat("============================================================\n")
cat("CYCLING DIRECT EVIDENCE\n")
cat("============================================================\n")
print(cycling_summary)

cat("\n")
cat("============================================================\n")
cat("NAS ASSOCIATION\n")
cat("============================================================\n")
print(nas_assoc)

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
