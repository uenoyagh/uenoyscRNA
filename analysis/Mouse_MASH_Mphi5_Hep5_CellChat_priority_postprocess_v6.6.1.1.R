#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Mphi5 <-> Hep5 CellChat priority post-processing
#
# Version: v6.6.1.1
#
# PURPOSE
# -------
# Re-use the completed v6.6.1 CellChat summary tables WITHOUT re-running
# CellChat, and correct the priority-pathway definition by explicitly adding:
#
#   PDGF
#
# The script:
#   1) reads v6.6.1 LR and pathway summary CSVs
#   2) re-flags priority pathways with PDGF included
#   3) exports updated priority tables
#   4) exports PDGF-only pathway/LR tables
#   5) makes compact priority pathway figures
#
# NO CellChat inference is performed.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}


save_pdf <- function(
  p,
  file,
  width,
  height
) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}


priority_pathway_flag_v6611 <- function(pathway_name) {

  grepl(
    paste(
      c(
        "^TNF$",
        "^IL1$",
        "^IL6$",
        "^OSM$",
        "^TGF",
        "TGFB",
        "^PDGF$",       # FIX in v6.6.1.1
        "SPP1",
        "FN1",
        "GAS",
        "AXL",
        "MERTK",
        "EGF",
        "IGF",
        "FGF",
        "CCL",
        "CXCL",
        "CSF",
        "MIF",
        "ANNEXIN",
        "HEDGEHOG",
        "SHH",
        "NOTCH",
        "VEGF",
        "COLLAGEN",
        "LAMININ",
        "THBS"
      ),
      collapse = "|"
    ),
    as.character(pathway_name),
    ignore.case = TRUE
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

IN_DIR <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_v6.6.1",
  "Tables"
)

LR_FILE <- file.path(
  IN_DIR,
  "06_bidirectional_LR_Sham_vs_Tx_replicate_summary_v6.6.1.csv"
)

PATHWAY_FILE <- file.path(
  IN_DIR,
  "10_bidirectional_pathway_Sham_vs_Tx_summary_v6.6.1.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_priority_postprocess_v6.6.1.1"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

for (
  d in c(
    OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Input checks
# ==============================================================================

for (
  f in c(
    LR_FILE,
    PATHWAY_FILE
  )
) {
  if (
    !file.exists(f)
  ) {
    stop(
      "Required v6.6.1 table missing: ",
      f
    )
  }
}

msg(
  "Reading completed v6.6.1 summary tables..."
)

lr <- read.csv(
  LR_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pathway <- read.csv(
  PATHWAY_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 4. Required-column audit
# ==============================================================================

required_lr <- c(
  "direction",
  "source",
  "target",
  "pathway_name",
  "interaction_key",
  "ligand",
  "receptor",
  "Sham_mean",
  "Tx_mean",
  "Tx_minus_Sham",
  "pairwise_direction_consistency",
  "evidence_grade"
)

required_pathway <- c(
  "direction",
  "source",
  "target",
  "pathway_name",
  "Sham_mean",
  "Tx_mean",
  "Tx_minus_Sham",
  "pairwise_direction_consistency",
  "evidence_grade",
  "detected_samples"
)

missing_lr <- setdiff(
  required_lr,
  colnames(lr)
)

missing_pathway <- setdiff(
  required_pathway,
  colnames(pathway)
)

if (
  length(missing_lr)
) {
  stop(
    "LR summary missing required columns: ",
    paste(
      missing_lr,
      collapse = ", "
    )
  )
}

if (
  length(missing_pathway)
) {
  stop(
    "Pathway summary missing required columns: ",
    paste(
      missing_pathway,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 5. Re-flag priority pathways with PDGF explicitly included
# ==============================================================================

lr <- lr %>%
  mutate(
    priority_pathway_v6611 =
      priority_pathway_flag_v6611(
        pathway_name
      )
  )

pathway <- pathway %>%
  mutate(
    priority_pathway_v6611 =
      priority_pathway_flag_v6611(
        pathway_name
      )
  )


# ==============================================================================
# 6. Save complete re-flagged tables
# ==============================================================================

write.csv(
  lr,
  file.path(
    TAB_OUT,
    "01_LR_summary_priority_reflagged_v6.6.1.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  pathway,
  file.path(
    TAB_OUT,
    "02_pathway_summary_priority_reflagged_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Updated priority tables
# ==============================================================================

priority_pathway <- pathway %>%
  filter(
    priority_pathway_v6611
  ) %>%
  arrange(
    direction,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

priority_lr <- lr %>%
  filter(
    priority_pathway_v6611
  ) %>%
  arrange(
    direction,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

write.csv(
  priority_pathway,
  file.path(
    TAB_OUT,
    "03_priority_crosslineage_pathways_v6.6.1.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  priority_lr,
  file.path(
    TAB_OUT,
    "04_priority_crosslineage_LR_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Biological focus tables
# ==============================================================================

mphi_to_injury <- priority_pathway %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

injury_to_mphi <- priority_pathway %>%
  filter(
    direction ==
      "Hep_to_Mphi",
    source ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

write.csv(
  mphi_to_injury,
  file.path(
    TAB_OUT,
    "05_Mphi5_to_InjuryHep_priority_pathways_v6.6.1.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  injury_to_mphi,
  file.path(
    TAB_OUT,
    "06_InjuryHep_to_Mphi5_priority_pathways_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. PDGF-specific extraction
# ==============================================================================

pdgf_pathway <- pathway %>%
  filter(
    grepl(
      "^PDGF$",
      pathway_name,
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    direction,
    source,
    target,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

pdgf_lr <- lr %>%
  filter(
    grepl(
      "^PDGF$",
      pathway_name,
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    direction,
    source,
    target,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

write.csv(
  pdgf_pathway,
  file.path(
    TAB_OUT,
    "07_PDGF_pathway_all_crosslineage_v6.6.1.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  pdgf_lr,
  file.path(
    TAB_OUT,
    "08_PDGF_LR_all_crosslineage_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. PDGF focus: Mphi -> Injury/inflammatory Hepatocyte
# ==============================================================================

pdgf_mphi_to_injury <- pdgf_pathway %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  )

write.csv(
  pdgf_mphi_to_injury,
  file.path(
    TAB_OUT,
    "09_PDGF_Mphi5_to_InjuryHep_focus_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Compact figure: priority pathways to Injury-Hep
# ==============================================================================

plot_mphi_to_injury <- mphi_to_injury %>%
  filter(
    detected_samples >=
      2
  ) %>%
  mutate(
    sender_pathway =
      paste0(
        source,
        " | ",
        pathway_name
      )
  ) %>%
  arrange(
    Tx_minus_Sham
  )

if (
  nrow(
    plot_mphi_to_injury
  )
) {

  p1 <- ggplot(
    plot_mphi_to_injury,
    aes(
      x =
        Tx_minus_Sham,
      y =
        reorder(
          sender_pathway,
          Tx_minus_Sham
        )
    )
  ) +
    geom_col() +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.35
    ) +
    labs(
      title =
        "Priority Mphi -> Injury/inflammatory Hepatocyte pathways",
      subtitle =
        "v6.6.1 CellChat results re-flagged with PDGF included",
      x =
        "Tx mean - Sham mean",
      y =
        NULL
    ) +
    theme_classic(
      base_size =
        8
    )

  save_pdf(
    p1,
    file.path(
      FIG_OUT,
      "01_Mphi5_to_InjuryHep_priority_pathways_v6.6.1.1.pdf"
    ),
    9,
    11
  )
}


# ==============================================================================
# 12. Compact figure: PDGF pathway
# ==============================================================================

plot_pdgf <- pdgf_mphi_to_injury %>%
  mutate(
    source =
      factor(
        source,
        levels =
          unique(
            source
          )
      )
  )

if (
  nrow(
    plot_pdgf
  )
) {

  pdgf_long <- bind_rows(
    plot_pdgf %>%
      transmute(
        source,
        condition =
          "Sham",
        strength =
          Sham_mean
      ),
    plot_pdgf %>%
      transmute(
        source,
        condition =
          "Tx",
        strength =
          Tx_mean
      )
  )

  p2 <- ggplot(
    pdgf_long,
    aes(
      x =
        source,
      y =
        strength,
      fill =
        condition
    )
  ) +
    geom_col(
      position =
        position_dodge(
          width =
            0.8
        ),
      width =
        0.72
    ) +
    labs(
      title =
        "PDGF: Mphi -> Injury/inflammatory Hepatocyte",
      subtitle =
        "Mean sample-wise CellChat pathway strength",
      x =
        NULL,
      y =
        "CellChat pathway strength",
      fill =
        NULL
    ) +
    theme_classic(
      base_size =
        8
    ) +
    theme(
      axis.text.x =
        element_text(
          angle =
            45,
          hjust =
            1
        )
    )

  save_pdf(
    p2,
    file.path(
      FIG_OUT,
      "02_PDGF_Mphi5_to_InjuryHep_v6.6.1.1.pdf"
    ),
    9,
    6
  )
}


# ==============================================================================
# 13. Audit: what changed relative to old priority flag
# ==============================================================================

if (
  "priority_pathway" %in%
    colnames(
      pathway
    )
) {

  priority_change <- pathway %>%
    mutate(
      old_priority =
        as.logical(
          priority_pathway
        )
    ) %>%
    filter(
      old_priority !=
        priority_pathway_v6611
    ) %>%
    select(
      direction,
      source,
      target,
      pathway_name,
      Sham_mean,
      Tx_mean,
      Tx_minus_Sham,
      evidence_grade,
      old_priority,
      priority_pathway_v6611
    )

} else {

  priority_change <- pathway %>%
    filter(
      priority_pathway_v6611,
      grepl(
        "^PDGF$",
        pathway_name,
        ignore.case = TRUE
      )
    ) %>%
    mutate(
      old_priority =
        NA
    ) %>%
    select(
      direction,
      source,
      target,
      pathway_name,
      Sham_mean,
      Tx_mean,
      Tx_minus_Sham,
      evidence_grade,
      old_priority,
      priority_pathway_v6611
    )
}

write.csv(
  priority_change,
  file.path(
    TAB_OUT,
    "10_priority_flag_changes_v6.6.1_to_v6.6.1.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Manifest / console summary
# ==============================================================================

manifest <- data.frame(
  parameter = c(
    "version",
    "source_LR_table",
    "source_pathway_table",
    "CellChat_recomputed",
    "change_from_v6.6.1",
    "new_priority_family"
  ),
  value = c(
    "v6.6.1.1",
    LR_FILE,
    PATHWAY_FILE,
    "FALSE",
    "Post-processing priority re-flag only",
    "PDGF"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.1.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.1.1.txt"
    )
)

msg(
  "Priority flag changes from v6.6.1:"
)

print(
  priority_change
)

msg(
  "PDGF Mphi -> Injury-Hep rows:"
)

print(
  pdgf_mphi_to_injury
)

msg(
  "DONE."
)

msg(
  "No CellChat inference was rerun."
)

msg(
  "Output directory: ",
  OUT
)
