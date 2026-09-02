#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Strict 3-way concordance:
# NicheNet + sender pseudobulk DE + replicate-aware CellChat
#
# Version: v6.3.1
#
# INPUTS
#   v6.3.0:
#     12_integrated_NicheNet_CellChat_sender_evidence_v6.3.0.csv
#     15_integrated_ligand_target_links_v6.3.0.csv
#
#   v6.2.2:
#     03_LR_refined_replicate_comparison_v6.2.2.csv
#
# PURPOSE
#   Convert the broad v6.3.0 candidate set into an explicitly evidence-graded
#   mechanistic shortlist.
#
# STRICT EVIDENCE CRITERIA
#
#   A. NicheNet
#      nichenet_rank <= 20
#
#   B. Sender expression
#      ligand expressed in >=10% of sender cells
#
#   C. Sender pseudobulk DE direction
#      Tx_up_program:
#        sender ligand logFC > 0
#      Tx_down_program:
#        sender ligand logFC < 0
#
#      "sender_DE_supported" additionally requires P <= 0.10.
#
#   D. CellChat direction
#      Tx_up_program:
#        delta_prob_supported > 0
#      Tx_down_program:
#        delta_prob_supported < 0
#
#   E. CellChat replicate support
#      total_support_replicates >= 2
#
#   F. CellChat pairwise consistency
#      pairwise_direction_consistency >= 0.75
#
# PRIMARY CLASSIFICATION
#
#   3_way_supported
#     NicheNet top20
#     + sender expressed
#     + sender DE supported in matching direction
#     + CellChat replicate-supported in matching direction
#
#   NicheNet_plus_CellChat
#     NicheNet top20 + matching replicate-supported CellChat
#     but no supported matching sender DE
#
#   NicheNet_plus_senderDE
#     NicheNet top20 + matching sender DE
#     but no matching replicate-supported CellChat
#
#   CellChat_plus_senderDE
#     matching CellChat + matching sender DE
#     but NicheNet rank >20
#
#   CellChat_only
#     matching replicate-supported CellChat only
#
#   NicheNet_only
#     NicheNet top20 only
#
#   senderDE_only
#     supported matching sender DE only
#
#   discordant
#     sender DE and/or CellChat provides directional evidence opposite to the
#     receiver program.
#
# IMPORTANT
#   - No CellChat rerun.
#   - No NicheNet rerun.
#   - No cell-level statistical testing.
#   - Biological n=2 Sham vs n=2 Tx remains exploratory.
#   - This script is designed to determine whether:
#       FN1 is CellChat-dominant,
#       SEMA4D is truly 3-way supported,
#       PDGFB is directionally discordant.
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

require_columns <- function(
  df,
  cols,
  label
) {
  missing <- setdiff(
    cols,
    colnames(df)
  )

  if (
    length(
      missing
    )
  ) {
    stop(
      label,
      " is missing required columns: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(
      x
    )
  )
}

program_expected_sign <- function(
  receiver_program
) {

  ifelse(
    receiver_program ==
      "Tx_up_program",
    1,
    ifelse(
      receiver_program ==
        "Tx_down_program",
      -1,
      NA_real_
    )
  )
}

sign_direction <- function(
  x,
  tol = 0
) {

  case_when(
    is.na(x) ~
      "missing",

    x >
      tol ~
      "Tx_up",

    x <
      -tol ~
      "Tx_down",

    TRUE ~
      "zero"
  )
}

priority_family <- function(gene) {

  case_when(
    gene ==
      "Fn1" ~
      "FN1",

    gene %in%
      c(
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Pdgfd"
      ) ~
      "PDGF",

    gene %in%
      c(
        "Sema4a",
        "Sema4b",
        "Sema4c",
        "Sema4d",
        "Sema4f",
        "Sema4g"
      ) ~
      "SEMA4",

    gene ==
      "App" ~
      "APP",

    gene ==
      "Thbs1" ~
      "THBS",

    gene ==
      "Spp1" ~
      "SPP1",

    gene %in%
      c(
        "Tgfb1",
        "Tgfb2",
        "Tgfb3"
      ) ~
      "TGFB",

    gene ==
      "Tnf" ~
      "TNF",

    TRUE ~
      "Other"
  )
}

grade_class <- function(x) {

  case_when(
    x ==
      "3_way_supported" ~
      1L,

    x %in%
      c(
        "NicheNet_plus_CellChat",
        "NicheNet_plus_senderDE",
        "CellChat_plus_senderDE"
      ) ~
      2L,

    x %in%
      c(
        "CellChat_only",
        "NicheNet_only",
        "senderDE_only"
      ) ~
      3L,

    x ==
      "discordant" ~
      4L,

    TRUE ~
      5L
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

V630 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_v6.3.0"
)

V622 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_CellChat_refined_v6.2.2"
)

INTEGRATED_FILE <- file.path(
  V630,
  "Tables",
  "12_integrated_NicheNet_CellChat_sender_evidence_v6.3.0.csv"
)

TARGET_FILE <- file.path(
  V630,
  "Tables",
  "15_integrated_ligand_target_links_v6.3.0.csv"
)

CELLCHAT_FILE <- file.path(
  V622,
  "Tables",
  "03_LR_refined_replicate_comparison_v6.2.2.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_concordance_v6.3.1"
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

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

for (
  d in c(
    OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT,
    RDS_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

NICHE_TOP_N <- 20L

MIN_SENDER_PCT <- 0.10

SENDER_P_THRESHOLD <- 0.10

CELLCHAT_MIN_SUPPORT_REPS <- 2L

CELLCHAT_MIN_PAIRWISE_CONSISTENCY <- 0.75

PRIORITY_LIGANDS <- c(
  "Fn1",
  "Sema4d",
  "Pdgfb",
  "App",
  "Thbs1",
  "Spp1",
  "Tgfb1",
  "Tnf"
)

CLASS_LEVELS <- c(
  "3_way_supported",
  "NicheNet_plus_CellChat",
  "NicheNet_plus_senderDE",
  "CellChat_plus_senderDE",
  "CellChat_only",
  "NicheNet_only",
  "senderDE_only",
  "discordant",
  "insufficient_evidence"
)


# ==============================================================================
# 4. Preflight
# ==============================================================================

for (
  f in c(
    INTEGRATED_FILE,
    TARGET_FILE,
    CELLCHAT_FILE
  )
) {
  if (
    !file.exists(
      f
    )
  ) {
    stop(
      "Required input not found: ",
      f
    )
  }
}

msg(
  "Reading v6.3.0 integrated evidence..."
)

integrated <- read.csv(
  INTEGRATED_FILE,
  check.names = FALSE
) %>%
  as_tibble()

msg(
  "Reading v6.3.0 ligand-target links..."
)

target_links <- read.csv(
  TARGET_FILE,
  check.names = FALSE
) %>%
  as_tibble()

msg(
  "Reading v6.2.2 full refined CellChat LR table..."
)

cellchat <- read.csv(
  CELLCHAT_FILE,
  check.names = FALSE
) %>%
  as_tibble()


# ==============================================================================
# 5. Validate columns
# ==============================================================================

require_columns(
  integrated,
  c(
    "sender",
    "receiver",
    "receiver_program",
    "ligand",
    "nichenet_rank",
    "nichenet_percentile",
    "sender_ligand_logFC",
    "sender_ligand_PValue",
    "sender_ligand_FDR",
    "sender_ligand_pct_expressed"
  ),
  "v6.3.0 integrated evidence"
)

require_columns(
  target_links,
  c(
    "sender",
    "receiver",
    "receiver_program",
    "ligand",
    "target",
    "regulatory_potential",
    "receiver_logFC",
    "receiver_PValue",
    "receiver_FDR"
  ),
  "v6.3.0 ligand-target links"
)

require_columns(
  cellchat,
  c(
    "source",
    "target",
    "ligand",
    "receptor",
    "delta_prob_supported",
    "total_support_replicates",
    "pairwise_direction_consistency"
  ),
  "v6.2.2 refined CellChat"
)


# ==============================================================================
# 6. Clean numeric columns
# ==============================================================================

for (
  nm in c(
    "nichenet_rank",
    "nichenet_percentile",
    "sender_ligand_logFC",
    "sender_ligand_PValue",
    "sender_ligand_FDR",
    "sender_ligand_pct_expressed"
  )
) {
  integrated[[nm]] <-
    safe_num(
      integrated[[nm]]
    )
}

for (
  nm in c(
    "delta_prob_supported",
    "total_support_replicates",
    "pairwise_direction_consistency"
  )
) {
  cellchat[[nm]] <-
    safe_num(
      cellchat[[nm]]
    )
}

for (
  nm in c(
    "regulatory_potential",
    "receiver_logFC",
    "receiver_PValue",
    "receiver_FDR"
  )
) {
  target_links[[nm]] <-
    safe_num(
      target_links[[nm]]
    )
}


# ==============================================================================
# 7. Rebuild CellChat ligand-level evidence from FULL v6.2.2 table
#
# Important:
#   v6.3.0 integrated table used the high-confidence shortlist.
#   Here we return to the full refined LR table so that absence from the
#   high-confidence shortlist is not mistaken for absence of CellChat evidence.
# ==============================================================================

cellchat <- cellchat %>%
  mutate(
    source =
      as.character(
        source
      ),
    target =
      as.character(
        target
      ),
    ligand =
      as.character(
        ligand
      ),
    receptor =
      as.character(
        receptor
      )
  )

cellchat_ligand <- cellchat %>%
  group_by(
    source,
    target,
    ligand
  ) %>%
  arrange(
    desc(
      abs(
        delta_prob_supported
      )
    ),
    desc(
      total_support_replicates
    ),
    desc(
      pairwise_direction_consistency
    ),
    .by_group = TRUE
  ) %>%
  summarise(
    CellChat_present =
      TRUE,

    CellChat_best_delta_prob =
      first(
        delta_prob_supported
      ),

    CellChat_support_replicates =
      max(
        total_support_replicates,
        na.rm = TRUE
      ),

    CellChat_pairwise_consistency =
      max(
        pairwise_direction_consistency,
        na.rm = TRUE
      ),

    CellChat_receptors =
      paste(
        sort(
          unique(
            receptor[
              !is.na(
                receptor
              ) &
                receptor !=
                  ""
            ]
          )
        ),
        collapse = "|"
      ),

    CellChat_pathways =
      "",

    .groups = "drop"
  ) %>%
  rename(
    sender =
      source,
    receiver =
      target
  )

# Rebuild pathways safely without relying on cur_data_all() for newer dplyr.
if (
  "pathway_name" %in%
    colnames(
      cellchat
    )
) {

  cc_pathways <- cellchat %>%
    group_by(
      source,
      target,
      ligand
    ) %>%
    summarise(
      CellChat_pathways_v2 =
        paste(
          sort(
            unique(
              pathway_name[
                !is.na(
                  pathway_name
                ) &
                  pathway_name !=
                    ""
              ]
            )
          ),
          collapse = "|"
        ),
      .groups = "drop"
    ) %>%
    rename(
      sender =
        source,
      receiver =
        target
    )

  cellchat_ligand <- cellchat_ligand %>%
    select(
      -CellChat_pathways
    ) %>%
    left_join(
      cc_pathways,
      by = c(
        "sender",
        "receiver",
        "ligand"
      )
    ) %>%
    rename(
      CellChat_pathways =
        CellChat_pathways_v2
    )
}


# ==============================================================================
# 8. One row per sender x receiver x program x ligand
# ==============================================================================

candidate_base <- integrated %>%
  select(
    sender,
    receiver,
    receiver_program,
    ligand,
    everything()
  ) %>%
  arrange(
    sender,
    receiver,
    receiver_program,
    ligand,
    nichenet_rank
  ) %>%
  group_by(
    sender,
    receiver,
    receiver_program,
    ligand
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup()

candidate <- candidate_base %>%
  select(
    -any_of(
      c(
        "CellChat_supported",
        "CellChat_best_delta_prob",
        "CellChat_direction_consistency",
        "CellChat_support_replicates",
        "CellChat_confidence_class",
        "CellChat_receptors",
        "CellChat_pathways",
        "CellChat_direction_match",
        "CellChat_replicate_supported",
        "evidence_count"
      )
    )
  ) %>%
  left_join(
    cellchat_ligand,
    by = c(
      "sender",
      "receiver",
      "ligand"
    )
  )


# ==============================================================================
# 9. Strict evidence flags
# ==============================================================================

candidate <- candidate %>%
  mutate(
    expected_sign =
      program_expected_sign(
        receiver_program
      ),

    expected_direction =
      ifelse(
        expected_sign >
          0,
        "Tx_up",
        "Tx_down"
      ),

    NicheNet_top20 =
      is.finite(
        nichenet_rank
      ) &
        nichenet_rank <=
          NICHE_TOP_N,

    sender_expressed =
      is.finite(
        sender_ligand_pct_expressed
      ) &
        sender_ligand_pct_expressed >=
          MIN_SENDER_PCT,

    sender_direction =
      sign_direction(
        sender_ligand_logFC
      ),

    sender_direction_match =
      is.finite(
        sender_ligand_logFC
      ) &
        sign(
          sender_ligand_logFC
        ) ==
          expected_sign,

    sender_DE_supported =
      sender_expressed &
        sender_direction_match &
        is.finite(
          sender_ligand_PValue
        ) &
        sender_ligand_PValue <=
          SENDER_P_THRESHOLD,

    sender_direction_opposite =
      is.finite(
        sender_ligand_logFC
      ) &
        sign(
          sender_ligand_logFC
        ) ==
          -expected_sign,

    CellChat_present =
      replace_na(
        CellChat_present,
        FALSE
      ),

    CellChat_direction =
      sign_direction(
        CellChat_best_delta_prob
      ),

    CellChat_direction_match =
      CellChat_present &
        is.finite(
          CellChat_best_delta_prob
        ) &
        sign(
          CellChat_best_delta_prob
        ) ==
          expected_sign,

    CellChat_direction_opposite =
      CellChat_present &
        is.finite(
          CellChat_best_delta_prob
        ) &
        sign(
          CellChat_best_delta_prob
        ) ==
          -expected_sign,

    CellChat_replicate_supported =
      CellChat_present &
        is.finite(
          CellChat_support_replicates
        ) &
        CellChat_support_replicates >=
          CELLCHAT_MIN_SUPPORT_REPS &
        is.finite(
          CellChat_pairwise_consistency
        ) &
        CellChat_pairwise_consistency >=
          CELLCHAT_MIN_PAIRWISE_CONSISTENCY,

    CellChat_supported_matching =
      CellChat_replicate_supported &
        CellChat_direction_match,

    CellChat_supported_opposite =
      CellChat_replicate_supported &
        CellChat_direction_opposite,

    directional_discordance =
      CellChat_supported_opposite |
        (
          sender_expressed &
            sender_direction_opposite &
            is.finite(
              sender_ligand_PValue
            ) &
            sender_ligand_PValue <=
              SENDER_P_THRESHOLD
        ) |
        (
          CellChat_replicate_supported &
            sender_expressed &
            is.finite(
              sender_ligand_logFC
            ) &
            is.finite(
              CellChat_best_delta_prob
            ) &
            sign(
              sender_ligand_logFC
            ) !=
              sign(
                CellChat_best_delta_prob
              ) &
            sender_ligand_logFC !=
              0 &
            CellChat_best_delta_prob !=
              0
        ),

    strict_evidence_count =
      as.integer(
        NicheNet_top20
      ) +
        as.integer(
          sender_DE_supported
        ) +
        as.integer(
          CellChat_supported_matching
        ),

    evidence_class = case_when(

      directional_discordance ~
        "discordant",

      NicheNet_top20 &
        sender_DE_supported &
        CellChat_supported_matching ~
        "3_way_supported",

      NicheNet_top20 &
        CellChat_supported_matching &
        !sender_DE_supported ~
        "NicheNet_plus_CellChat",

      NicheNet_top20 &
        sender_DE_supported &
        !CellChat_supported_matching ~
        "NicheNet_plus_senderDE",

      !NicheNet_top20 &
        sender_DE_supported &
        CellChat_supported_matching ~
        "CellChat_plus_senderDE",

      !NicheNet_top20 &
        !sender_DE_supported &
        CellChat_supported_matching ~
        "CellChat_only",

      NicheNet_top20 &
        !sender_DE_supported &
        !CellChat_supported_matching ~
        "NicheNet_only",

      !NicheNet_top20 &
        sender_DE_supported &
        !CellChat_supported_matching ~
        "senderDE_only",

      TRUE ~
        "insufficient_evidence"
    ),

    evidence_grade =
      grade_class(
        evidence_class
      ),

    mechanism_family =
      priority_family(
        ligand
      ),

    priority_axis =
      ligand %in%
        PRIORITY_LIGANDS,

    concordance_score =
      100 *
        as.integer(
          evidence_class ==
            "3_way_supported"
        ) +
        20 *
          strict_evidence_count +
        10 *
          as.integer(
            NicheNet_top20
          ) +
        5 *
          as.integer(
            CellChat_supported_matching
          ) +
        5 *
          as.integer(
            sender_DE_supported
          ) +
        pmax(
          0,
          1 -
            pmin(
              nichenet_rank,
              100
            ) /
              100
        )
  ) %>%
  arrange(
    evidence_grade,
    desc(
      concordance_score
    ),
    nichenet_rank
  )


# ==============================================================================
# 10. Main tables
# ==============================================================================

write.csv(
  candidate,
  file.path(
    TAB_OUT,
    "01_all_ligand_3way_concordance_v6.3.1.csv"
  ),
  row.names = FALSE
)

three_way <- candidate %>%
  filter(
    evidence_class ==
      "3_way_supported"
  )

write.csv(
  three_way,
  file.path(
    TAB_OUT,
    "02_three_way_supported_ligands_v6.3.1.csv"
  ),
  row.names = FALSE
)

nn_cc <- candidate %>%
  filter(
    evidence_class ==
      "NicheNet_plus_CellChat"
  )

write.csv(
  nn_cc,
  file.path(
    TAB_OUT,
    "03_NicheNet_plus_CellChat_v6.3.1.csv"
  ),
  row.names = FALSE
)

nn_de <- candidate %>%
  filter(
    evidence_class ==
      "NicheNet_plus_senderDE"
  )

write.csv(
  nn_de,
  file.path(
    TAB_OUT,
    "04_NicheNet_plus_senderDE_v6.3.1.csv"
  ),
  row.names = FALSE
)

cc_only <- candidate %>%
  filter(
    evidence_class %in%
      c(
        "CellChat_plus_senderDE",
        "CellChat_only"
      )
  )

write.csv(
  cc_only,
  file.path(
    TAB_OUT,
    "05_CellChat_dominant_candidates_v6.3.1.csv"
  ),
  row.names = FALSE
)

discordant <- candidate %>%
  filter(
    evidence_class ==
      "discordant"
  )

write.csv(
  discordant,
  file.path(
    TAB_OUT,
    "06_directionally_discordant_candidates_v6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Priority-axis audit
# ==============================================================================

priority_audit <- candidate %>%
  filter(
    ligand %in%
      PRIORITY_LIGANDS |
      mechanism_family %in%
        c(
          "FN1",
          "PDGF",
          "SEMA4",
          "APP",
          "THBS",
          "SPP1",
          "TGFB",
          "TNF"
        )
  ) %>%
  select(
    sender,
    receiver,
    receiver_program,
    ligand,
    mechanism_family,
    evidence_class,
    strict_evidence_count,
    NicheNet_top20,
    nichenet_rank,
    nichenet_percentile,
    sender_expressed,
    sender_ligand_pct_expressed,
    sender_ligand_logFC,
    sender_ligand_PValue,
    sender_ligand_FDR,
    sender_direction,
    sender_direction_match,
    sender_DE_supported,
    CellChat_present,
    CellChat_best_delta_prob,
    CellChat_direction,
    CellChat_direction_match,
    CellChat_support_replicates,
    CellChat_pairwise_consistency,
    CellChat_replicate_supported,
    CellChat_supported_matching,
    CellChat_receptors,
    CellChat_pathways,
    directional_discordance
  ) %>%
  arrange(
    mechanism_family,
    ligand,
    evidence_grade =
      grade_class(
        evidence_class
      ),
    sender,
    receiver,
    receiver_program
  )

write.csv(
  priority_audit,
  file.path(
    TAB_OUT,
    "07_priority_axis_audit_FN1_SEMA4D_PDGFB_APP_etc_v6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Dedicated FN1 / SEMA4D / PDGFB status table
# ==============================================================================

key_axis_status <- candidate %>%
  filter(
    ligand %in%
      c(
        "Fn1",
        "Sema4d",
        "Pdgfb"
      )
  ) %>%
  mutate(
    key_question = case_when(
      ligand ==
        "Fn1" ~
        "Is_FN1_CellChat_dominant",

      ligand ==
        "Sema4d" ~
        "Is_SEMA4D_3way_supported",

      ligand ==
        "Pdgfb" ~
        "Is_PDGFB_directionally_discordant",

      TRUE ~
        ""
    )
  ) %>%
  select(
    key_question,
    sender,
    receiver,
    receiver_program,
    ligand,
    evidence_class,
    strict_evidence_count,
    NicheNet_top20,
    nichenet_rank,
    sender_ligand_logFC,
    sender_ligand_PValue,
    sender_DE_supported,
    CellChat_best_delta_prob,
    CellChat_support_replicates,
    CellChat_pairwise_consistency,
    CellChat_supported_matching,
    directional_discordance,
    CellChat_receptors
  ) %>%
  arrange(
    ligand,
    grade =
      grade_class(
        evidence_class
      ),
    sender,
    receiver,
    receiver_program
  )

write.csv(
  key_axis_status,
  file.path(
    TAB_OUT,
    "08_FN1_SEMA4D_PDGFB_status_v6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Evidence class counts
# ==============================================================================

class_counts <- candidate %>%
  count(
    receiver,
    receiver_program,
    evidence_class,
    name =
      "n_candidates"
  ) %>%
  complete(
    receiver,
    receiver_program,
    evidence_class =
      CLASS_LEVELS,
    fill = list(
      n_candidates = 0
    )
  )

write.csv(
  class_counts,
  file.path(
    TAB_OUT,
    "09_evidence_class_counts_v6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Build exact ligand-receptor-target tables
#
# Receptors come from the full refined v6.2.2 CellChat LR table.
# Targets come from v6.3.0 NicheNet regulatory-potential links.
# ==============================================================================

cc_exact <- cellchat %>%
  filter(
    !is.na(
      ligand
    ),
    ligand !=
      "",
    !is.na(
      receptor
    ),
    receptor !=
      ""
  ) %>%
  transmute(
    sender =
      source,
    receiver =
      target,
    ligand,
    receptor,
    pathway_name =
      if (
        "pathway_name" %in%
          colnames(
            cellchat
          )
      ) {
        as.character(
          pathway_name
        )
      } else {
        NA_character_
      },
    interaction_key =
      if (
        "interaction_key" %in%
          colnames(
            cellchat
          )
      ) {
        as.character(
          interaction_key
        )
      } else {
        NA_character_
      },
    delta_prob_supported,
    total_support_replicates,
    pairwise_direction_consistency
  )

lrt <- target_links %>%
  inner_join(
    candidate %>%
      select(
        sender,
        receiver,
        receiver_program,
        ligand,
        evidence_class,
        evidence_grade,
        strict_evidence_count,
        NicheNet_top20,
        nichenet_rank,
        nichenet_percentile,
        sender_ligand_pct_expressed,
        sender_ligand_logFC,
        sender_ligand_PValue,
        sender_DE_supported,
        CellChat_supported_matching,
        directional_discordance,
        mechanism_family,
        concordance_score
      ),
    by = c(
      "sender",
      "receiver",
      "receiver_program",
      "ligand"
    )
  ) %>%
  left_join(
    cc_exact,
    by = c(
      "sender",
      "receiver",
      "ligand"
    )
  ) %>%
  mutate(
    receptor =
      replace_na(
        receptor,
        "No_supported_CellChat_receptor"
      ),

    LRT_score =
      concordance_score +
        10 *
          pmin(
            pmax(
              regulatory_potential,
              0
            ),
            1
          ) +
        2 *
          abs(
            replace_na(
              receiver_logFC,
              0
            )
          )
  ) %>%
  arrange(
    evidence_grade,
    desc(
      LRT_score
    ),
    desc(
      regulatory_potential
    )
  )

# Top links from 3-way supported ligands.
lrt_three_way <- lrt %>%
  filter(
    evidence_class ==
      "3_way_supported"
  ) %>%
  group_by(
    sender,
    receiver,
    receiver_program,
    ligand
  ) %>%
  slice_max(
    order_by =
      LRT_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  lrt_three_way,
  file.path(
    TAB_OUT,
    "10_three_way_ligand_receptor_target_links_v6.3.1.csv"
  ),
  row.names = FALSE
)

# Best 20 L-R-target chains overall from evidence class 1 or 2.
lrt_top20 <- lrt %>%
  filter(
    evidence_grade <=
      2
  ) %>%
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    receptor,
    target,
    .keep_all = TRUE
  ) %>%
  slice_max(
    order_by =
      LRT_score,
    n = 20,
    with_ties = FALSE
  )

write.csv(
  lrt_top20,
  file.path(
    TAB_OUT,
    "11_top20_mechanistic_ligand_receptor_target_chains_v6.3.1.csv"
  ),
  row.names = FALSE
)

# Broader top 50 for manuscript / mechanistic review.
lrt_top50 <- lrt %>%
  filter(
    evidence_grade <=
      3
  ) %>%
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    receptor,
    target,
    .keep_all = TRUE
  ) %>%
  slice_max(
    order_by =
      LRT_score,
    n = 50,
    with_ties = FALSE
  )

write.csv(
  lrt_top50,
  file.path(
    TAB_OUT,
    "12_top50_supported_ligand_receptor_target_chains_v6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Figure 1: evidence-class counts
# ==============================================================================

class_counts_plot <- class_counts %>%
  mutate(
    receiver_program_label =
      paste0(
        receiver,
        " | ",
        receiver_program
      ),

    evidence_class =
      factor(
        evidence_class,
        levels =
          CLASS_LEVELS
      )
  )

p1 <- ggplot(
  class_counts_plot,
  aes(
    x =
      receiver_program_label,
    y =
      n_candidates,
    fill =
      evidence_class
  )
) +
  geom_col(
    position = "stack"
  ) +
  labs(
    title =
      "Mphi -> HSC ligand evidence classes",
    subtitle =
      "Strict integration of NicheNet, sender pseudobulk DE, and replicate-aware CellChat",
    x = NULL,
    y =
      "Candidate ligands",
    fill =
      "Evidence class"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_evidence_class_counts_v6.3.1.pdf"
  ),
  11,
  6
)


# ==============================================================================
# 16. Figure 2: strict concordance candidate dot plot
# ==============================================================================

strict_plot_df <- candidate %>%
  filter(
    evidence_grade <=
      3,
    evidence_class !=
      "discordant"
  ) %>%
  arrange(
    evidence_grade,
    desc(
      concordance_score
    )
  ) %>%
  group_by(
    receiver,
    receiver_program
  ) %>%
  slice_head(
    n = 25
  ) %>%
  ungroup() %>%
  mutate(
    label =
      paste0(
        sender,
        " | ",
        ligand
      ),

    receiver_program_label =
      paste0(
        receiver,
        " | ",
        receiver_program
      )
  )

if (
  nrow(
    strict_plot_df
  )
) {

  p2 <- ggplot(
    strict_plot_df,
    aes(
      x =
        receiver_program_label,
      y =
        label,
      size =
        nichenet_percentile,
      fill =
        strict_evidence_count
    )
  ) +
    geom_point(
      shape = 21
    ) +
    scale_fill_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    labs(
      title =
        "Strict concordant Mphi ligand -> HSC program candidates",
      subtitle =
        "Point size = NicheNet percentile; fill = number of matching evidence layers",
      x = NULL,
      y = NULL,
      size =
        "NicheNet percentile",
      fill =
        "Evidence layers"
    ) +
    theme_classic(
      base_size = 7
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_pdf(
    p2,
    file.path(
      FIG_OUT,
      "02_strict_concordant_candidates_v6.3.1.pdf"
    ),
    12,
    max(
      8,
      min(
        18,
        4 +
          0.16 *
            n_distinct(
              strict_plot_df$label
            )
      )
    )
  )
}


# ==============================================================================
# 17. Figure 3: FN1 / SEMA4D / PDGFB audit
# ==============================================================================

key_plot_df <- candidate %>%
  filter(
    ligand %in%
      c(
        "Fn1",
        "Sema4d",
        "Pdgfb"
      )
  ) %>%
  mutate(
    label =
      paste0(
        sender,
        " -> ",
        receiver,
        " | ",
        receiver_program
      ),

    plot_x =
      replace_na(
        sender_ligand_logFC,
        0
      ),

    plot_y =
      replace_na(
        CellChat_best_delta_prob,
        0
      )
  )

if (
  nrow(
    key_plot_df
  )
) {

  p3 <- ggplot(
    key_plot_df,
    aes(
      x =
        plot_x,
      y =
        plot_y,
      size =
        nichenet_percentile,
      fill =
        evidence_class
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3
    ) +
    geom_point(
      shape = 21
    ) +
    facet_wrap(
      ~ ligand,
      scales = "free",
      ncol = 3
    ) +
    labs(
      title =
        "FN1 / SEMA4D / PDGFB three-way concordance audit",
      subtitle =
        "x = sender ligand Tx/Sham logFC; y = CellChat Tx-Sham supported probability",
      x =
        "Sender ligand log2FC Tx vs Sham",
      y =
        "CellChat delta probability",
      size =
        "NicheNet percentile",
      fill =
        "Evidence class"
    ) +
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p3,
    file.path(
      FIG_OUT,
      "03_FN1_SEMA4D_PDGFB_concordance_audit_v6.3.1.pdf"
    ),
    13,
    7
  )
}


# ==============================================================================
# 18. Figure 4: top mechanistic L-R-target chains
# ==============================================================================

if (
  nrow(
    lrt_top20
  )
) {

  top20_plot <- lrt_top20 %>%
    mutate(
      chain =
        paste0(
          sender,
          " | ",
          ligand,
          " -> ",
          receptor,
          " -> ",
          target
        ),

      chain =
        factor(
          chain,
          levels =
            rev(
              unique(
                chain
              )
            )
        )
    )

  p4 <- ggplot(
    top20_plot,
    aes(
      x =
        regulatory_potential,
      y =
        chain,
      size =
        abs(
          receiver_logFC
        ),
      fill =
        strict_evidence_count
    )
  ) +
    geom_point(
      shape = 21
    ) +
    scale_fill_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    labs(
      title =
        "Top mechanistic Mphi ligand -> HSC receptor -> target chains",
      subtitle =
        "Evidence-grade 1-2 candidates only",
      x =
        "NicheNet regulatory potential",
      y = NULL,
      size =
        "|HSC log2FC|",
      fill =
        "Evidence layers"
    ) +
    theme_classic(
      base_size = 7
    )

  save_pdf(
    p4,
    file.path(
      FIG_OUT,
      "04_top20_ligand_receptor_target_chains_v6.3.1.pdf"
    ),
    13,
    9
  )
}


# ==============================================================================
# 19. Save compact RDS
# ==============================================================================

results <- list(
  all_candidates =
    candidate,
  three_way_supported =
    three_way,
  NicheNet_plus_CellChat =
    nn_cc,
  NicheNet_plus_senderDE =
    nn_de,
  CellChat_dominant =
    cc_only,
  discordant =
    discordant,
  priority_axis_audit =
    priority_audit,
  key_axis_status =
    key_axis_status,
  top20_LRT =
    lrt_top20,
  top50_LRT =
    lrt_top50
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_Mphi5_HSC3_concordance_results_v6.3.1.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 20. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "v630_integrated_input",
    "v630_target_input",
    "v622_CellChat_input",
    "NicheNet_top_N",
    "sender_expression_threshold",
    "sender_DE_P_threshold",
    "CellChat_min_support_replicates",
    "CellChat_min_pairwise_consistency",
    "formal_inference_note"
  ),
  value = c(
    "v6.3.1",
    INTEGRATED_FILE,
    TARGET_FILE,
    CELLCHAT_FILE,
    as.character(
      NICHE_TOP_N
    ),
    as.character(
      MIN_SENDER_PCT
    ),
    as.character(
      SENDER_P_THRESHOLD
    ),
    as.character(
      CELLCHAT_MIN_SUPPORT_REPS
    ),
    as.character(
      CELLCHAT_MIN_PAIRWISE_CONSISTENCY
    ),
    "Exploratory: biological n=2 Sham vs n=2 Tx"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.1.txt"
    )
)

msg(
  "3-way supported ligands: ",
  nrow(
    three_way
  )
)

msg(
  "Directionally discordant ligands: ",
  nrow(
    discordant
  )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
