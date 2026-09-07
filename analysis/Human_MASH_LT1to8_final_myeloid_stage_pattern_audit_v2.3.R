#!/usr/bin/env Rscript

# Human_MASH_LT1to8_final_myeloid_stage_pattern_audit_v2.3.R
#
# PURPOSE
# -------
# Classify disease-associated myeloid gene trajectories into descriptive
# temporal/severity patterns before pathway analysis or additional F3/F4 samples.
#
# CLINICAL DESIGN
# ---------------
# LT7 = NAS0 / F0
# LT8 = NAS0 / F0
# LT5 = NAS2 / F0
# LT6 = NAS2 / F0
# LT3 = NAS4 / F0
# LT4 = NAS4 / F1
# LT1 = NAS6 / F4
# LT2 = NAS6 / F3
#
# KEY QUESTION
# ------------
# For each final myeloid state, does a candidate gene show:
#   1) Early F0-associated change (already visible by NAS2/F0 and maintained)
#   2) NAS4-onset/intermediate change
#   3) Advanced F3/F4-enriched shift
#   4) Gradual monotonic disease-severity change
#   5) Non-monotonic / sample-sensitive behavior
#
# IMPORTANT
# ---------
# These are DESCRIPTIVE pattern classes, not inferential causal labels.
# NAS and fibrosis remain confounded in advanced disease.
# Same-NAS fibrosis contrasts are n=1 vs n=1 and are reported only as diagnostics.
#
# INPUTS
# ------
# - TMM logCPM matrices from v2.2.1
# - overall disease-severity trend DE table from v2.2.1
#
# THIS SCRIPT DOES:
#   - calculate stage means and adjacent changes
#   - calculate F0-only early trajectory diagnostics
#   - calculate same-NAS fibrosis diagnostic contrasts
#   - assign transparent descriptive pattern classes
#   - rank candidate genes within each final myeloid state
#   - generate state-level trajectory and stage-mean heatmap figures
#
# THIS SCRIPT DOES NOT:
#   - perform new DE
#   - perform pathway enrichment
#   - modify/save Seurat objects
#   - infer a pure fibrosis effect
#   - infer a pure NAS effect
#
# Thresholds are intentionally explicit and configurable below.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v2.3"

DE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_PseudobulkDE_v2.2.1"
)

LOGCPM_RDS <- file.path(
  DE_DIR,
  "objects/Human_MASH_LT1to8_FinalMyeloid_TMM_logCPM_v2.2.1.rds"
)

TREND_CSV <- file.path(
  DE_DIR,
  "tables/02_all_states_overall_severity_trend_v2.2.1.csv"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_StagePatternAudit_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

for (p in c(LOGCPM_RDS, TREND_CSV)) {
  if (!file.exists(p)) stop("Required input not found:\n", p)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 final myeloid stage-pattern audit ", VERSION)
message("============================================================")
message("DESCRIPTIVE pattern classification only.")

# ------------------------------------------------------------------
# Explicit pattern thresholds
# ------------------------------------------------------------------

MIN_MEANINGFUL_DELTA <- 0.50
MIN_ADVANCED_DELTA <- 0.75
MIN_TOTAL_ENDPOINT <- 0.75
MIN_STEP_FOR_MONOTONIC <- 0.10

# ------------------------------------------------------------------
# Sample design
# ------------------------------------------------------------------

SAMPLE_ORDER <- c("LT7", "LT8", "LT5", "LT6", "LT3", "LT4", "LT1", "LT2")

sample_info <- data.frame(
  sample = SAMPLE_ORDER,
  NAS = c(0, 0, 2, 2, 4, 4, 6, 6),
  fibrosis = c("F0", "F0", "F0", "F0", "F0", "F1", "F4", "F3"),
  stringsAsFactors = FALSE
)

STATE_ORDER <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

logcpm_list <- readRDS(LOGCPM_RDS)
trend <- read.csv(TREND_CSV, check.names = FALSE)

required_trend <- c(
  "final_state",
  "gene",
  "logFC_per_severity_step_NAS2",
  "PValue",
  "FDR"
)

miss <- setdiff(required_trend, colnames(trend))
if (length(miss)) {
  stop("Trend table missing columns: ", paste(miss, collapse = ", "))
}

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

same_sign <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || a == 0 || b == 0) return(FALSE)
  sign(a) == sign(b)
}

both_advanced_concordant <- function(v1, v2, ref, direction) {
  if (direction > 0) {
    return(v1 > ref && v2 > ref)
  } else if (direction < 0) {
    return(v1 < ref && v2 < ref)
  }
  FALSE
}

classify_pattern <- function(
  m0, m2, m4, m6,
  f0_nas4,
  lt1, lt2
) {
  d02 <- m2 - m0
  d24 <- m4 - m2
  d46 <- m6 - m4
  dend <- m6 - m0
  d04 <- m4 - m0

  if (!is.finite(dend) || abs(dend) < MIN_TOTAL_ENDPOINT) {
    return("Low_total_effect")
  }

  direction <- sign(dend)

  # Advanced F3/F4-enriched shift:
  # most change appears after NAS4 and BOTH advanced samples move beyond NAS4 mean.
  if (
    abs(d04) < MIN_ADVANCED_DELTA &&
    abs(d46) >= MIN_ADVANCED_DELTA &&
    same_sign(d46, dend) &&
    both_advanced_concordant(lt1, lt2, m4, direction)
  ) {
    return("Advanced_F3F4_enriched_shift")
  }

  # Early F0-associated:
  # already changed by NAS2/F0 and the F0 NAS4 sample remains in same direction.
  if (
    abs(d02) >= MIN_MEANINGFUL_DELTA &&
    same_sign(d02, dend) &&
    same_sign(f0_nas4 - m0, dend) &&
    abs(f0_nas4 - m0) >= MIN_MEANINGFUL_DELTA
  ) {
    return("Early_F0_associated_sustained")
  }

  # NAS4 onset:
  # little early change, then clear NAS4 movement in endpoint direction.
  if (
    abs(d02) < MIN_MEANINGFUL_DELTA &&
    abs(d24) >= MIN_MEANINGFUL_DELTA &&
    same_sign(d24, dend)
  ) {
    return("NAS4_onset_or_intermediate_shift")
  }

  # Gradual monotonic:
  # all adjacent means move in the endpoint direction.
  if (
    same_sign(d02, dend) &&
    same_sign(d24, dend) &&
    same_sign(d46, dend) &&
    abs(d02) >= MIN_STEP_FOR_MONOTONIC &&
    abs(d24) >= MIN_STEP_FOR_MONOTONIC &&
    abs(d46) >= MIN_STEP_FOR_MONOTONIC
  ) {
    return("Gradual_monotonic_severity_shift")
  }

  "Nonmonotonic_or_sample_sensitive"
}

priority_class <- function(fdr, p, effect) {
  if (is.finite(fdr) && fdr < 0.05) return("Tier1_FDR05")
  if (is.finite(fdr) && fdr < 0.10) return("Tier2_FDR10")
  if (is.finite(p) && p < 0.01 && abs(effect) >= 0.75) return("Tier3_P01_effect")
  if (is.finite(p) && p < 0.05 && abs(effect) >= 0.50) return("Tier4_P05_effect")
  "Background"
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ------------------------------------------------------------------
# Build gene-level stage-pattern table
# ------------------------------------------------------------------

rows <- list()

for (st in STATE_ORDER) {
  if (!st %in% names(logcpm_list)) {
    stop("Missing state in logCPM RDS: ", st)
  }

  mat <- logcpm_list[[st]]

  if (!all(SAMPLE_ORDER %in% colnames(mat))) {
    stop("Missing one or more samples in state: ", st)
  }

  tr <- trend[trend$final_state == st, , drop = FALSE]
  tr <- tr[tr$gene %in% rownames(mat), , drop = FALSE]

  for (i in seq_len(nrow(tr))) {
    g <- tr$gene[i]
    x <- as.numeric(mat[g, SAMPLE_ORDER])
    names(x) <- SAMPLE_ORDER

    m0 <- mean(x[c("LT7", "LT8")])
    m2 <- mean(x[c("LT5", "LT6")])
    m4 <- mean(x[c("LT3", "LT4")])
    m6 <- mean(x[c("LT1", "LT2")])

    d02 <- m2 - m0
    d24 <- m4 - m2
    d46 <- m6 - m4
    dend <- m6 - m0

    # F0-only early diagnostic: NAS0/F0 -> NAS2/F0 -> NAS4/F0(LT3 only).
    f0_nas4 <- x[["LT3"]]

    # Same-NAS fibrosis diagnostics, n=1 vs n=1, descriptive only.
    nas4_F1_minus_F0 <- x[["LT4"]] - x[["LT3"]]
    nas6_F4_minus_F3 <- x[["LT1"]] - x[["LT2"]]

    pattern <- classify_pattern(
      m0 = m0,
      m2 = m2,
      m4 = m4,
      m6 = m6,
      f0_nas4 = f0_nas4,
      lt1 = x[["LT1"]],
      lt2 = x[["LT2"]]
    )

    rows[[length(rows) + 1L]] <- data.frame(
      final_state = st,
      gene = g,

      mean_NAS0_F0 = m0,
      mean_NAS2_F0 = m2,
      mean_NAS4_mixed_F0_F1 = m4,
      mean_NAS6_F3_F4 = m6,

      delta_NAS2_vs_NAS0 = d02,
      delta_NAS4_vs_NAS2 = d24,
      delta_NAS6_vs_NAS4 = d46,
      delta_NAS6_vs_NAS0 = dend,

      LT3_NAS4_F0_logCPM = x[["LT3"]],
      LT4_NAS4_F1_logCPM = x[["LT4"]],
      NAS4_F1_minus_F0_descriptive = nas4_F1_minus_F0,

      LT2_NAS6_F3_logCPM = x[["LT2"]],
      LT1_NAS6_F4_logCPM = x[["LT1"]],
      NAS6_F4_minus_F3_descriptive = nas6_F4_minus_F3,

      logFC_per_severity_step_NAS2 =
        tr$logFC_per_severity_step_NAS2[i],
      PValue = tr$PValue[i],
      FDR = tr$FDR[i],

      priority = priority_class(
        tr$FDR[i],
        tr$PValue[i],
        tr$logFC_per_severity_step_NAS2[i]
      ),

      pattern_class = pattern,

      stringsAsFactors = FALSE
    )
  }
}

pattern_table <- do.call(rbind, rows)
rownames(pattern_table) <- NULL

write.csv(
  pattern_table,
  file.path(TABDIR, "01_all_gene_stage_pattern_table_v2.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Candidate-focused table
# ------------------------------------------------------------------

candidate <- pattern_table[
  pattern_table$priority != "Background" &
    pattern_table$pattern_class != "Low_total_effect",
  ,
  drop = FALSE
]

priority_levels <- c(
  "Tier1_FDR05",
  "Tier2_FDR10",
  "Tier3_P01_effect",
  "Tier4_P05_effect",
  "Background"
)

candidate$priority <- factor(
  candidate$priority,
  levels = priority_levels
)

candidate <- candidate[
  order(
    match(candidate$final_state, STATE_ORDER),
    candidate$priority,
    candidate$FDR,
    candidate$PValue
  ),
  ,
  drop = FALSE
]

candidate$priority <- as.character(candidate$priority)

write.csv(
  candidate,
  file.path(TABDIR, "02_priority_candidate_stage_patterns_v2.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Pattern summaries
# ------------------------------------------------------------------

pattern_counts <- as.data.frame(
  table(
    final_state = pattern_table$final_state,
    pattern_class = pattern_table$pattern_class
  ),
  stringsAsFactors = FALSE
)
pattern_counts <- pattern_counts[pattern_counts$Freq > 0, , drop = FALSE]
colnames(pattern_counts)[3] <- "n_genes"

write.csv(
  pattern_counts,
  file.path(TABDIR, "03_pattern_counts_all_tested_genes_v2.3.csv"),
  row.names = FALSE
)

candidate_pattern_counts <- as.data.frame(
  table(
    final_state = candidate$final_state,
    pattern_class = candidate$pattern_class,
    priority = candidate$priority
  ),
  stringsAsFactors = FALSE
)
candidate_pattern_counts <- candidate_pattern_counts[
  candidate_pattern_counts$Freq > 0,
  ,
  drop = FALSE
]
colnames(candidate_pattern_counts)[4] <- "n_genes"

write.csv(
  candidate_pattern_counts,
  file.path(TABDIR, "04_priority_candidate_pattern_counts_v2.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Compact top candidate table per state
# ------------------------------------------------------------------

top_list <- list()

for (st in STATE_ORDER) {
  z <- candidate[candidate$final_state == st, , drop = FALSE]

  if (!nrow(z)) {
    z <- pattern_table[
      pattern_table$final_state == st,
      ,
      drop = FALSE
    ]
    z <- z[order(z$FDR, z$PValue), , drop = FALSE]
  }

  top_list[[st]] <- utils::head(z, 20)
}

top20 <- do.call(rbind, top_list)
rownames(top20) <- NULL

write.csv(
  top20,
  file.path(TABDIR, "05_top20_stage_pattern_candidates_by_state_v2.3.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures: top trajectory candidates per state
# ------------------------------------------------------------------

for (st in STATE_ORDER) {
  z <- top20[top20$final_state == st, , drop = FALSE]
  if (!nrow(z)) next

  genes <- unique(utils::head(z$gene, 12))
  mat <- logcpm_list[[st]]

  genes <- genes[genes %in% rownames(mat)]
  if (!length(genes)) next

  traj <- do.call(
    rbind,
    lapply(genes, function(g) {
      zz <- pattern_table[
        pattern_table$final_state == st &
          pattern_table$gene == g,
        ,
        drop = FALSE
      ]

      data.frame(
        gene = g,
        sample = SAMPLE_ORDER,
        NAS = sample_info$NAS,
        fibrosis = sample_info$fibrosis,
        logCPM = as.numeric(mat[g, SAMPLE_ORDER]),
        pattern_class = zz$pattern_class[1],
        priority = zz$priority[1],
        stringsAsFactors = FALSE
      )
    })
  )

  traj$gene <- factor(traj$gene, levels = genes)

  p <- ggplot(
    traj,
    aes(
      x = NAS,
      y = logCPM,
      label = paste0(sample, "/", fibrosis)
    )
  ) +
    geom_point(size = 2.1) +
    geom_text(vjust = -0.7, size = 2.25) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.5
    ) +
    facet_wrap(
      ~ gene,
      scales = "free_y",
      ncol = 3
    ) +
    scale_x_continuous(breaks = c(0, 2, 4, 6)) +
    theme_classic(base_size = 9) +
    labs(
      x = "NAS (fibrosis shown in labels)",
      y = "TMM-normalized logCPM"
    ) +
    ggtitle(
      paste0(st, " | top stage-pattern candidates"),
      "Descriptive only; advanced NAS6 samples are F3/F4"
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "01_top_stage_pattern_trajectories_",
        safe_name(st),
        "_v2.3.pdf"
      )
    ),
    p,
    width = 10,
    height = 10
  )
}

# ------------------------------------------------------------------
# Figures: stage-mean heatmap for top candidates
# ------------------------------------------------------------------

for (st in STATE_ORDER) {
  z <- top20[top20$final_state == st, , drop = FALSE]
  if (!nrow(z)) next

  z <- utils::head(z, 15)

  heat <- data.frame(
    gene = rep(z$gene, each = 4),
    stage = rep(
      c("NAS0_F0", "NAS2_F0", "NAS4_F0_F1", "NAS6_F3_F4"),
      times = nrow(z)
    ),
    mean_logCPM = as.vector(
      t(
        as.matrix(
          z[
            , c(
              "mean_NAS0_F0",
              "mean_NAS2_F0",
              "mean_NAS4_mixed_F0_F1",
              "mean_NAS6_F3_F4"
            )
          ]
        )
      )
    ),
    stringsAsFactors = FALSE
  )

  heat$gene <- factor(heat$gene, levels = rev(z$gene))
  heat$stage <- factor(
    heat$stage,
    levels = c("NAS0_F0", "NAS2_F0", "NAS4_F0_F1", "NAS6_F3_F4")
  )

  p_heat <- ggplot(
    heat,
    aes(
      x = stage,
      y = gene,
      fill = mean_logCPM
    )
  ) +
    geom_tile() +
    geom_text(
      aes(label = sprintf("%.1f", mean_logCPM)),
      size = 2.6
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      name = "Mean\nlogCPM"
    ) +
    theme_classic(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    ggtitle(
      paste0(st, " | top stage-pattern candidates"),
      "Stage means are descriptive; NAS4 and NAS6 include fibrosis heterogeneity"
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "02_stage_mean_heatmap_",
        safe_name(st),
        "_v2.3.pdf"
      )
    ),
    p_heat,
    width = 7.8,
    height = max(6, 0.34 * nrow(z) + 2.5)
  )
}

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

tier1 <- candidate[candidate$priority == "Tier1_FDR05", , drop = FALSE]
tier2 <- candidate[candidate$priority == "Tier2_FDR10", , drop = FALSE]

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final myeloid stage-pattern audit ", VERSION),
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
  "DESCRIPTIVE PATTERN CLASSES:",
  "Early_F0_associated_sustained",
  "NAS4_onset_or_intermediate_shift",
  "Advanced_F3F4_enriched_shift",
  "Gradual_monotonic_severity_shift",
  "Nonmonotonic_or_sample_sensitive",
  "Low_total_effect",
  "",
  "THRESHOLDS:",
  paste0("MIN_MEANINGFUL_DELTA = ", MIN_MEANINGFUL_DELTA, " log2 units"),
  paste0("MIN_ADVANCED_DELTA = ", MIN_ADVANCED_DELTA, " log2 units"),
  paste0("MIN_TOTAL_ENDPOINT = ", MIN_TOTAL_ENDPOINT, " log2 units"),
  paste0("MIN_STEP_FOR_MONOTONIC = ", MIN_STEP_FOR_MONOTONIC, " log2 units"),
  "",
  "PRIORITY CANDIDATES:",
  paste0("Tier1 FDR<0.05: ", nrow(tier1)),
  paste0("Tier2 FDR<0.10: ", nrow(tier2)),
  paste0("All non-background priority candidates: ", nrow(candidate)),
  "",
  "TIER1 DETAILS:",
  if (nrow(tier1)) {
    paste(
      capture.output(
        print(
          tier1[
            , c(
              "final_state",
              "gene",
              "logFC_per_severity_step_NAS2",
              "FDR",
              "pattern_class"
            )
          ],
          row.names = FALSE
        )
      ),
      collapse = "\n"
    )
  } else {
    "<none>"
  },
  "",
  "CRITICAL INTERPRETATION:",
  "- Pattern labels are descriptive trajectory categories, not causal labels.",
  "- Advanced_F3F4_enriched_shift indicates that most observed change appears after NAS4",
  "  and both NAS6/F3-F4 samples move in the same endpoint direction.",
  "- Early_F0_associated_sustained requires evidence already within F0 samples.",
  "- NAS4_F1_minus_F0 and NAS6_F4_minus_F3 are n=1 vs n=1 diagnostics only.",
  "- These outputs are intended to prioritize which signals should be retested when new",
  "  F3/F4 and NAS/fibrosis-overlap samples are added.",
  "",
  "NEXT DECISION GATE:",
  "Review the Tier1/Tier2 pattern assignments and top trajectories.",
  "Then decide which gene programs proceed to pathway analysis now and which should wait",
  "for additional fibrosis-balanced samples."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "06_stage_pattern_summary_v2.3.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v2.3 COMPLETE")
message("Output: ", OUTDIR)
message("No new DE or pathway analysis performed.")
message("============================================================")
