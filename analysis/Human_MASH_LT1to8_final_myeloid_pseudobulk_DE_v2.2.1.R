#!/usr/bin/env Rscript

# Human_MASH_LT1to8_final_myeloid_pseudobulk_DE_v2.2.1.R
#
# PURPOSE
# -------
# Corrected and fibrosis-aware interpretation of sample-level state-specific
# pseudobulk DE for the FINAL Human MΦ/Monocyte states.
#
# IMPORTANT CLINICAL DESIGN
# -------------------------
# LT7 = NAS0 / F0
# LT8 = NAS0 / F0
# LT5 = NAS2 / F0
# LT6 = NAS2 / F0
# LT3 = NAS4 / F0
# LT4 = NAS4 / F1
# LT1 = NAS6 / F4
# LT2 = NAS6 / F3
#
# Therefore NAS and fibrosis are strongly confounded at advanced disease.
# This script DOES NOT fit NAS and fibrosis simultaneously as covariates.
#
# TESTS
# -----
# PRIMARY:
#   Overall disease-severity trend across all 8 biological samples.
#   Numerically encoded as NAS_step = NAS/2, but interpreted as an overall
#   disease-severity axis, NOT a pure NAS effect.
#
# SECONDARY A:
#   Advanced disease: NAS6/F3-4 vs NAS0/F0.
#
# SECONDARY B:
#   Early fibrosis-matched activity contrast: NAS2/F0 vs NAS0/F0.
#
# All tests use biological samples as replicates.
#
# FIXES vs v2.2
# --------------
# - edgeR::makeContrasts() removed (not exported by edgeR).
# - Explicit numeric contrast vectors are used.
# - edgeR::normLibSizes() is used instead of deprecated calcNormFactors().
# - Fibrosis metadata and interpretation are carried explicitly.
#
# THIS SCRIPT DOES NOT:
#   - perform pathway enrichment
#   - modify/save Seurat objects
#   - treat cells as statistical replicates
#   - claim that the overall trend is a pure NAS effect
#
# Requires Bioconductor package: edgeR

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v2.2.1"

PB_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_PseudobulkQC_v2.1/",
  "objects/Human_MASH_LT1to8_FinalMyeloid_pseudobulk_counts_v2.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_PseudobulkDE_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PB_RDS)) {
  stop("v2.1 pseudobulk count RDS not found:\n", PB_RDS)
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "Bioconductor package 'edgeR' is required but not installed.\n",
    "Do not continue until edgeR is available."
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 final myeloid pseudobulk DE ", VERSION)
message("============================================================")
message("Statistical unit: biological sample only.")
message("Advanced fibrosis is confounded with NAS6; interpret accordingly.")

pb_list <- readRDS(PB_RDS)

SAMPLE_ORDER <- c("LT7", "LT8", "LT5", "LT6", "LT3", "LT4", "LT1", "LT2")

sample_info <- data.frame(
  sample = SAMPLE_ORDER,
  NAS = c(0, 0, 2, 2, 4, 4, 6, 6),
  fibrosis = factor(
    c("F0", "F0", "F0", "F0", "F0", "F1", "F4", "F3"),
    levels = c("F0", "F1", "F2", "F3", "F4")
  ),
  stringsAsFactors = FALSE
)

sample_info$NAS_step <- sample_info$NAS / 2
sample_info$NAS_stage <- factor(
  sample_info$NAS,
  levels = c(0, 2, 4, 6)
)

write.csv(
  sample_info,
  file.path(TABDIR, "00_sample_clinical_design_v2.2.1.csv"),
  row.names = FALSE
)

STATE_ORDER <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

state_tier <- c(
  "Resident/Kupffer-like_Mphi" = "PRIMARY_robust",
  "Repair/Resolution_Mphi" = "PRIMARY_robust",
  "Lipid-associated/TREM2_Mphi" = "EXPLORATORY_lower_abundance",
  "Inflammatory/Transitional_Mphi" = "PRIMARY_robust",
  "Inflammatory_Monocyte" = "EXPLORATORY_lower_abundance"
)

if (!all(STATE_ORDER %in% names(pb_list))) {
  stop(
    "Pseudobulk RDS is missing expected state(s): ",
    paste(setdiff(STATE_ORDER, names(pb_list)), collapse = " | ")
  )
}

tier_manifest <- data.frame(
  final_state = STATE_ORDER,
  analysis_tier = unname(state_tier[STATE_ORDER]),
  rationale = c(
    "v2.1 min 159 cells/sample; min library 429,223",
    "v2.1 min 161 cells/sample; min library 251,062",
    "v2.1 one sample with 19 cells; min library 31,410",
    "v2.1 min 47 cells/sample; min library 58,537",
    "v2.1 one sample with 22 cells; min library 29,635"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  tier_manifest,
  file.path(TABDIR, "01_state_analysis_tier_manifest_v2.2.1.csv"),
  row.names = FALSE
)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

make_volcano <- function(tab, logfc_col, title, outfile) {
  df <- tab
  df$minus_log10_FDR <- -log10(pmax(df$FDR, 1e-300))
  df$status <- "NS"
  df$status[df$FDR < 0.05 & df[[logfc_col]] >= 1] <- "Up_FDR05_logFC1"
  df$status[df$FDR < 0.05 & df[[logfc_col]] <= -1] <- "Down_FDR05_logFC1"
  df$status[df$FDR < 0.05 & abs(df[[logfc_col]]) < 1] <- "FDR05_small_effect"

  p <- ggplot(
    df,
    aes(
      x = .data[[logfc_col]],
      y = minus_log10_FDR,
      color = status
    )
  ) +
    geom_point(size = 0.8, alpha = 0.65) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    theme_classic(base_size = 10) +
    theme(legend.title = element_blank()) +
    labs(
      x = logfc_col,
      y = "-log10(FDR)"
    ) +
    ggtitle(title)

  ggsave(outfile, p, width = 7.5, height = 6)
}

all_trend <- list()
all_adv <- list()
all_early <- list()
summary_rows <- list()
top_rows <- list()
logcpm_list <- list()

for (st in STATE_ORDER) {
  message("------------------------------------------------------------")
  message("State: ", st)

  counts <- pb_list[[st]]

  if (!all(SAMPLE_ORDER %in% colnames(counts))) {
    stop("State missing one or more samples: ", st)
  }

  counts <- counts[, SAMPLE_ORDER, drop = FALSE]

  y <- edgeR::DGEList(
    counts = counts,
    samples = sample_info
  )

  keep <- edgeR::filterByExpr(
    y,
    group = sample_info$NAS_stage
  )

  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::normLibSizes(y, method = "TMM")

  if (nrow(y) < 100) {
    stop("Too few genes remain after filterByExpr for state: ", st)
  }

  logcpm <- edgeR::cpm(
    y,
    log = TRUE,
    prior.count = 1,
    normalized.lib.sizes = TRUE
  )

  logcpm_list[[st]] <- logcpm

  # ------------------------------------------------------------
  # PRIMARY: overall disease-severity trend
  # ------------------------------------------------------------

  design_trend <- model.matrix(
    ~ NAS_step,
    data = sample_info
  )
  rownames(design_trend) <- sample_info$sample

  y_trend <- edgeR::estimateDisp(
    y,
    design = design_trend,
    robust = TRUE
  )

  fit_trend <- edgeR::glmQLFit(
    y_trend,
    design = design_trend,
    robust = TRUE
  )

  qlf_trend <- edgeR::glmQLFTest(
    fit_trend,
    coef = "NAS_step"
  )

  trend <- edgeR::topTags(
    qlf_trend,
    n = Inf,
    sort.by = "none"
  )$table

  trend$gene <- rownames(trend)
  rownames(trend) <- NULL
  colnames(trend)[colnames(trend) == "logFC"] <-
    "logFC_per_severity_step_NAS2"

  trend$final_state <- st
  trend$analysis_tier <- unname(state_tier[[st]])
  trend$interpretation_axis <-
    "overall_disease_severity_NAS_plus_fibrosis_confounded"

  trend <- trend[
    , c(
      "final_state",
      "analysis_tier",
      "interpretation_axis",
      "gene",
      "logFC_per_severity_step_NAS2",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]

  all_trend[[st]] <- trend

  write.csv(
    trend,
    file.path(
      TABDIR,
      paste0("trend_severity_", safe_name(st), "_v2.2.1.csv")
    ),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # Stage-factor design for endpoint contrasts
  # ------------------------------------------------------------

  design_stage <- model.matrix(
    ~ 0 + NAS_stage,
    data = sample_info
  )
  rownames(design_stage) <- sample_info$sample
  colnames(design_stage) <- c("NAS0", "NAS2", "NAS4", "NAS6")

  y_stage <- edgeR::estimateDisp(
    y,
    design = design_stage,
    robust = TRUE
  )

  fit_stage <- edgeR::glmQLFit(
    y_stage,
    design = design_stage,
    robust = TRUE
  )

  # Explicit contrast vectors.
  contrast_60 <- c(
    NAS0 = -1,
    NAS2 = 0,
    NAS4 = 0,
    NAS6 = 1
  )

  contrast_20 <- c(
    NAS0 = -1,
    NAS2 = 1,
    NAS4 = 0,
    NAS6 = 0
  )

  # ------------------------------------------------------------
  # SECONDARY A: Advanced disease NAS6/F3-4 vs NAS0/F0
  # ------------------------------------------------------------

  qlf_adv <- edgeR::glmQLFTest(
    fit_stage,
    contrast = contrast_60
  )

  adv <- edgeR::topTags(
    qlf_adv,
    n = Inf,
    sort.by = "none"
  )$table

  adv$gene <- rownames(adv)
  rownames(adv) <- NULL
  colnames(adv)[colnames(adv) == "logFC"] <-
    "logFC_advanced_NAS6F3_4_vs_NAS0F0"

  adv$final_state <- st
  adv$analysis_tier <- unname(state_tier[[st]])

  adv <- adv[
    , c(
      "final_state",
      "analysis_tier",
      "gene",
      "logFC_advanced_NAS6F3_4_vs_NAS0F0",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]

  all_adv[[st]] <- adv

  write.csv(
    adv,
    file.path(
      TABDIR,
      paste0("advanced_NAS6F3_4_vs_NAS0F0_", safe_name(st), "_v2.2.1.csv")
    ),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # SECONDARY B: Early fibrosis-matched NAS2/F0 vs NAS0/F0
  # ------------------------------------------------------------

  qlf_early <- edgeR::glmQLFTest(
    fit_stage,
    contrast = contrast_20
  )

  early <- edgeR::topTags(
    qlf_early,
    n = Inf,
    sort.by = "none"
  )$table

  early$gene <- rownames(early)
  rownames(early) <- NULL
  colnames(early)[colnames(early) == "logFC"] <-
    "logFC_early_NAS2F0_vs_NAS0F0"

  early$final_state <- st
  early$analysis_tier <- unname(state_tier[[st]])

  early <- early[
    , c(
      "final_state",
      "analysis_tier",
      "gene",
      "logFC_early_NAS2F0_vs_NAS0F0",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]

  all_early[[st]] <- early

  write.csv(
    early,
    file.path(
      TABDIR,
      paste0("early_NAS2F0_vs_NAS0F0_", safe_name(st), "_v2.2.1.csv")
    ),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # Summary counts
  # ------------------------------------------------------------

  summary_rows[[st]] <- data.frame(
    final_state = st,
    analysis_tier = unname(state_tier[[st]]),
    genes_tested = nrow(trend),

    severity_trend_FDR05 = sum(trend$FDR < 0.05, na.rm = TRUE),
    severity_trend_FDR05_abslogFC1 = sum(
      trend$FDR < 0.05 &
        abs(trend$logFC_per_severity_step_NAS2) >= 1,
      na.rm = TRUE
    ),

    advanced_FDR05 = sum(adv$FDR < 0.05, na.rm = TRUE),
    advanced_FDR05_abslogFC1 = sum(
      adv$FDR < 0.05 &
        abs(adv$logFC_advanced_NAS6F3_4_vs_NAS0F0) >= 1,
      na.rm = TRUE
    ),

    early_F0matched_FDR05 = sum(early$FDR < 0.05, na.rm = TRUE),
    early_F0matched_FDR05_abslogFC1 = sum(
      early$FDR < 0.05 &
        abs(early$logFC_early_NAS2F0_vs_NAS0F0) >= 1,
      na.rm = TRUE
    ),

    common_dispersion_trend = y_trend$common.dispersion,
    common_dispersion_stage = y_stage$common.dispersion,

    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------
  # Top genes compact table
  # ------------------------------------------------------------

  pick_top <- function(tab, effect_col, test_name) {
    z <- tab[order(tab$FDR, tab$PValue), , drop = FALSE]
    z <- utils::head(z, 15)

    data.frame(
      final_state = st,
      analysis_tier = unname(state_tier[[st]]),
      test = test_name,
      gene = z$gene,
      logFC = z[[effect_col]],
      PValue = z$PValue,
      FDR = z$FDR,
      stringsAsFactors = FALSE
    )
  }

  top_rows[[paste0(st, "_trend")]] <- pick_top(
    trend,
    "logFC_per_severity_step_NAS2",
    "overall_disease_severity_trend"
  )

  top_rows[[paste0(st, "_advanced")]] <- pick_top(
    adv,
    "logFC_advanced_NAS6F3_4_vs_NAS0F0",
    "advanced_NAS6F3_4_vs_NAS0F0"
  )

  top_rows[[paste0(st, "_early")]] <- pick_top(
    early,
    "logFC_early_NAS2F0_vs_NAS0F0",
    "early_NAS2F0_vs_NAS0F0"
  )

  # ------------------------------------------------------------
  # Volcano plots
  # ------------------------------------------------------------

  make_volcano(
    trend,
    "logFC_per_severity_step_NAS2",
    paste0(st, " | overall disease-severity trend"),
    file.path(
      FIGDIR,
      paste0("volcano_severity_", safe_name(st), "_v2.2.1.pdf")
    )
  )

  make_volcano(
    adv,
    "logFC_advanced_NAS6F3_4_vs_NAS0F0",
    paste0(st, " | advanced NAS6/F3-4 vs NAS0/F0"),
    file.path(
      FIGDIR,
      paste0("volcano_advanced_", safe_name(st), "_v2.2.1.pdf")
    )
  )

  make_volcano(
    early,
    "logFC_early_NAS2F0_vs_NAS0F0",
    paste0(st, " | early NAS2/F0 vs NAS0/F0"),
    file.path(
      FIGDIR,
      paste0("volcano_early_F0matched_", safe_name(st), "_v2.2.1.pdf")
    )
  )

  # ------------------------------------------------------------
  # Top disease-severity trend gene trajectories
  # ------------------------------------------------------------

  candidate <- trend[
    trend$FDR < 0.10,
    ,
    drop = FALSE
  ]

  if (nrow(candidate) < 6) {
    candidate <- trend[
      order(trend$FDR, trend$PValue),
      ,
      drop = FALSE
    ]
  }

  genes_plot <- utils::head(candidate$gene, 12)
  genes_plot <- genes_plot[genes_plot %in% rownames(logcpm)]

  if (length(genes_plot)) {
    traj <- do.call(
      rbind,
      lapply(genes_plot, function(g) {
        data.frame(
          gene = g,
          sample = SAMPLE_ORDER,
          NAS = sample_info$NAS,
          fibrosis = as.character(sample_info$fibrosis),
          logCPM = as.numeric(logcpm[g, SAMPLE_ORDER]),
          stringsAsFactors = FALSE
        )
      })
    )

    traj$gene <- factor(traj$gene, levels = genes_plot)

    p_traj <- ggplot(
      traj,
      aes(
        x = NAS,
        y = logCPM,
        label = paste0(sample, "/", fibrosis)
      )
    ) +
      geom_point(size = 2) +
      geom_text(
        vjust = -0.7,
        size = 2.3
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linewidth = 0.55
      ) +
      facet_wrap(
        ~ gene,
        scales = "free_y",
        ncol = 3
      ) +
      scale_x_continuous(
        breaks = c(0, 2, 4, 6)
      ) +
      theme_classic(base_size = 9) +
      labs(
        x = "NAS (fibrosis shown in sample labels)",
        y = "TMM-normalized logCPM"
      ) +
      ggtitle(
        paste0(st, " | top overall disease-severity genes"),
        "Advanced NAS6 samples are F3/F4; trend is not a pure NAS effect"
      )

    ggsave(
      file.path(
        FIGDIR,
        paste0(
          "trajectory_top_severity_genes_",
          safe_name(st),
          "_v2.2.1.pdf"
        )
      ),
      p_traj,
      width = 10,
      height = 10
    )
  }
}

# ------------------------------------------------------------------
# Combined outputs
# ------------------------------------------------------------------

trend_all <- do.call(rbind, all_trend)
adv_all <- do.call(rbind, all_adv)
early_all <- do.call(rbind, all_early)
de_summary <- do.call(rbind, summary_rows)
top_all <- do.call(rbind, top_rows)

rownames(trend_all) <- NULL
rownames(adv_all) <- NULL
rownames(early_all) <- NULL
rownames(de_summary) <- NULL
rownames(top_all) <- NULL

write.csv(
  trend_all,
  file.path(TABDIR, "02_all_states_overall_severity_trend_v2.2.1.csv"),
  row.names = FALSE
)

write.csv(
  adv_all,
  file.path(TABDIR, "03_all_states_advanced_NAS6F3_4_vs_NAS0F0_v2.2.1.csv"),
  row.names = FALSE
)

write.csv(
  early_all,
  file.path(TABDIR, "04_all_states_early_NAS2F0_vs_NAS0F0_v2.2.1.csv"),
  row.names = FALSE
)

write.csv(
  de_summary,
  file.path(TABDIR, "05_DE_state_summary_v2.2.1.csv"),
  row.names = FALSE
)

write.csv(
  top_all,
  file.path(TABDIR, "06_top15_genes_each_state_test_v2.2.1.csv"),
  row.names = FALSE
)

saveRDS(
  logcpm_list,
  file.path(
    OBJDIR,
    "Human_MASH_LT1to8_FinalMyeloid_TMM_logCPM_v2.2.1.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------------

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final myeloid pseudobulk DE ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  "CLINICAL DESIGN:",
  "LT7 NAS0/F0",
  "LT8 NAS0/F0",
  "LT5 NAS2/F0",
  "LT6 NAS2/F0",
  "LT3 NAS4/F0",
  "LT4 NAS4/F1",
  "LT1 NAS6/F4",
  "LT2 NAS6/F3",
  "",
  "CRITICAL INTERPRETATION:",
  "NAS and fibrosis are strongly confounded at advanced disease.",
  "The continuous trend is therefore labeled OVERALL DISEASE-SEVERITY TREND,",
  "not a pure NAS-specific effect.",
  "",
  "PRIMARY TEST:",
  "Overall disease-severity trend across all 8 biological samples.",
  "Effect is reported per +2 NAS points numerically, but advanced samples are also F3/F4.",
  "",
  "SECONDARY A:",
  "Advanced disease: NAS6/F3-4 vs NAS0/F0.",
  "",
  "SECONDARY B:",
  "Early fibrosis-matched activity: NAS2/F0 vs NAS0/F0.",
  "",
  "STATE TIERS:",
  paste(
    capture.output(print(tier_manifest, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "DE SUMMARY:",
  paste(
    capture.output(print(de_summary, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION:",
  "- FDR < 0.05 is the primary multiple-testing threshold.",
  "- LAM and Inflammatory Monocyte remain exploratory because of lower pseudobulk support.",
  "- Do not fit NAS and fibrosis together as covariates in the current 8-sample dataset;",
  "  the design is too collinear.",
  "- Additional F3/F4 samples will improve reproducibility/power, but only samples that",
  "  break the NAS-fibrosis correlation will help disentangle activity from fibrosis.",
  "",
  "NEXT DECISION GATE:",
  "Review DE counts, top genes, volcanoes, and sample-level trajectories before pathway analysis."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "07_run_summary_v2.2.1.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v2.2.1 COMPLETE")
message("NO pathway analysis performed.")
message("Output: ", OUTDIR)
message("============================================================")
