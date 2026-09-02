#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6631)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Mphi5 -> Hepatocyte strict 3-way concordance post-processing
#
# Version: v6.6.3.1
#
# PURPOSE
# -------
# Post-process the completed v6.6.3 NicheNet analysis WITHOUT rerunning:
#   - NicheNet
#   - edgeR
#   - CellChat
#
# FIXES / IMPROVEMENTS
# --------------------
# 1) Correct NicheNet percentile calculation.
#    v6.6.3 used base::ifelse() with a scalar test inside mutate(), which caused
#    nichenet_percentile to become 1 for all rows.
#
# 2) Rebuild sender expansion explicitly.
#    NicheNet activity is receiver/program/ligand-level, while a ligand can be
#    expressed by multiple macrophage sender states. This is an EXPECTED
#    one-ligand-to-many-senders expansion. It is made explicit here and reduced
#    to one row per:
#       sender x receiver x receiver_program x ligand
#
# 3) Rebuild ligand-level CellChat evidence from the FULL v6.6.2.1 refined LR
#    table, rather than treating absence from the high-confidence shortlist as
#    absence of CellChat evidence.
#
# 4) Apply a strict concordance classification analogous to the prior
#    Mphi-HSC v6.3.1 workflow.
#
# PRIMARY RECEIVER
# ----------------
#   Hep_Injury-inflammatory
#
# KEY AXES
# --------
#   Pdgfb
#   Sema4d
#   Fn1
#   Tnf
#   Plau
#
# Secondary audit:
#   Spp1
#
# STRICT EVIDENCE DEFINITIONS
# ---------------------------
# A. NicheNet
#    nichenet_rank <= 20
#
# B. Sender expression
#    ligand expressed in >=10% of sender cells
#
# C. Sender pseudobulk direction
#    Tx_up_program:   sender ligand logFC > 0
#    Tx_down_program: sender ligand logFC < 0
#
#    sender_DE_supported additionally requires nominal P <= 0.10.
#    FDR is retained but is NOT required because biological n=2/group.
#
# D. CellChat direction
#    Tx_up_program:   supported Tx-Sham probability > 0
#    Tx_down_program: supported Tx-Sham probability < 0
#
# E. CellChat replicate support
#    total_support_samples >= 3/4
#
# F. CellChat pairwise consistency
#    supported_pairwise_direction_consistency >= 0.75
#
# EVIDENCE CLASSES
# ----------------
#   3_way_supported
#     NicheNet top20 + senderDE matching + CellChat matching
#
#   NicheNet_plus_CellChat
#   NicheNet_plus_senderDE
#   CellChat_plus_senderDE
#   CellChat_only
#   NicheNet_only
#   senderDE_only
#   discordant
#   insufficient_evidence
#
# IMPORTANT INTERPRETATION
# ------------------------
# NicheNet ranks ligand ability to explain the selected receiver gene program.
# It does NOT establish the actual in-vivo direction of ligand activity.
#
# Therefore a high NicheNet rank with sender/CellChat moving in the opposite
# direction is explicitly labeled "discordant", not mechanistically supported.
#
# Biological n:
#   Sham n=2
#   Tx   n=2
#
# All inference remains exploratory.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
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
    length(missing)
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
    as.numeric(x)
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


mechanism_family <- function(gene) {

  gene_chr <- as.character(gene)

  case_when(
    gene_chr ==
      "Fn1" ~
      "FN1",

    gene_chr %in%
      c(
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Pdgfd"
      ) ~
      "PDGF",

    gene_chr ==
      "Spp1" ~
      "SPP1",

    grepl(
      "^Sema4",
      gene_chr
    ) ~
      "SEMA4",

    gene_chr ==
      "Tnf" ~
      "TNF",

    gene_chr ==
      "Plau" ~
      "PLAU",

    gene_chr ==
      "Mif" ~
      "MIF",

    gene_chr ==
      "Tgfb1" ~
      "TGFB",

    gene_chr ==
      "Gas6" ~
      "GAS6",

    gene_chr ==
      "Hgf" ~
      "HGF",

    gene_chr ==
      "App" ~
      "APP",

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


V663 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_NicheNet_v6.6.3"
)


V6621 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_refine_v6.6.2.1"
)


RESULT_RDS <- file.path(
  V663,
  "RDS",
  "Mouse_MASH_Mphi5_Hep5_NicheNet_results_v6.6.3.rds"
)


CELLCHAT_FILE <- file.path(
  V6621,
  "Tables",
  "03_LR_refined_replicate_comparison_v6.6.2.1.csv"
)


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_NicheNet_concordance_v6.6.3.1"
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

PRIMARY_RECEIVER <-
  "Hep_Injury-inflammatory"


NICHE_TOP_N <-
  20L


MIN_SENDER_PCT <-
  0.10


SENDER_P_THRESHOLD <-
  0.10


CELLCHAT_MIN_SUPPORT_SAMPLES <-
  3L


CELLCHAT_MIN_PAIRWISE_CONSISTENCY <-
  0.75


KEY_LIGANDS <- c(
  "Pdgfb",
  "Sema4d",
  "Fn1",
  "Tnf",
  "Plau"
)


SECONDARY_LIGANDS <- c(
  "Spp1"
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
    RESULT_RDS,
    CELLCHAT_FILE
  )
) {

  if (
    !file.exists(f)
  ) {

    stop(
      "Required input not found: ",
      f
    )
  }
}


msg(
  "Reading v6.6.3 saved NicheNet results..."
)


res <- readRDS(
  RESULT_RDS
)


required_result_names <- c(
  "nichenet_activity",
  "sender_DE",
  "expression_fraction",
  "sender_receiver_ligand_map",
  "target_links"
)


missing_result_names <- setdiff(
  required_result_names,
  names(res)
)


if (
  length(missing_result_names)
) {

  stop(
    "v6.6.3 result RDS is missing components: ",
    paste(
      missing_result_names,
      collapse = ", "
    )
  )
}


activity <- as_tibble(
  res$nichenet_activity
)


sender_de <- as_tibble(
  res$sender_DE
)


expr_fraction <- as_tibble(
  res$expression_fraction
)


sender_map <- as_tibble(
  res$sender_receiver_ligand_map
)


target_links <- as_tibble(
  res$target_links
)


msg(
  "Reading full v6.6.2.1 refined CellChat LR table..."
)


cellchat <- read.csv(
  CELLCHAT_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


# ==============================================================================
# 5. Validate inputs
# ==============================================================================

require_columns(
  activity,
  c(
    "test_ligand",
    "receiver",
    "receiver_program",
    "nichenet_rank",
    "aupr_corrected"
  ),
  "v6.6.3 nichenet_activity"
)


require_columns(
  sender_de,
  c(
    "sender",
    "gene",
    "logFC",
    "logCPM",
    "PValue",
    "FDR"
  ),
  "v6.6.3 sender_DE"
)


require_columns(
  expr_fraction,
  c(
    "gene",
    "celltype",
    "pct_expressed"
  ),
  "v6.6.3 expression_fraction"
)


require_columns(
  sender_map,
  c(
    "sender",
    "receiver",
    "ligand"
  ),
  "v6.6.3 sender_receiver_ligand_map"
)


require_columns(
  target_links,
  c(
    "ligand",
    "receiver",
    "receiver_program",
    "target",
    "regulatory_potential"
  ),
  "v6.6.3 target_links"
)


require_columns(
  cellchat,
  c(
    "direction",
    "source",
    "target",
    "ligand",
    "receptor",
    "pathway_name",
    "supported_Tx_minus_Sham",
    "total_support_samples",
    "supported_pairwise_direction_consistency"
  ),
  "v6.6.2.1 full refined CellChat"
)


# ==============================================================================
# 6. Correct NicheNet percentile
# ==============================================================================

activity_corrected <- activity %>%
  mutate(
    ligand =
      as.character(
        test_ligand
      ),
    nichenet_rank =
      safe_num(
        nichenet_rank
      ),
    aupr_corrected =
      safe_num(
        aupr_corrected
      )
  ) %>%
  group_by(
    receiver,
    receiver_program
  ) %>%
  arrange(
    nichenet_rank,
    .by_group = TRUE
  ) %>%
  mutate(
    nichenet_n_candidates_corrected =
      n(),

    nichenet_percentile_corrected =
      case_when(
        nichenet_n_candidates_corrected >
          1 ~
          1 -
            (
              nichenet_rank -
                1
            ) /
              (
                nichenet_n_candidates_corrected -
                  1
              ),

        TRUE ~
          1
      )
  ) %>%
  ungroup()


write.csv(
  activity_corrected,
  file.path(
    TAB_OUT,
    "01_NicheNet_activity_percentile_corrected_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


percentile_audit <- activity_corrected %>%
  group_by(
    receiver,
    receiver_program
  ) %>%
  summarise(
    n_candidates =
      n(),

    min_rank =
      min(
        nichenet_rank,
        na.rm = TRUE
      ),

    max_rank =
      max(
        nichenet_rank,
        na.rm = TRUE
      ),

    max_percentile =
      max(
        nichenet_percentile_corrected,
        na.rm = TRUE
      ),

    min_percentile =
      min(
        nichenet_percentile_corrected,
        na.rm = TRUE
      ),

    n_unique_percentiles =
      n_distinct(
        nichenet_percentile_corrected
      ),

    .groups = "drop"
  )


write.csv(
  percentile_audit,
  file.path(
    TAB_OUT,
    "02_NicheNet_percentile_fix_audit_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Explicit sender expansion
#
# NicheNet activity has one row per receiver/program/ligand.
# sender_map can contain the same receiver/ligand in multiple Mphi senders.
# That is biologically expected.
#
# We deliberately expand activity to sender-specific rows and then enforce
# uniqueness at:
#   sender x receiver x receiver_program x ligand
# ==============================================================================

sender_map_unique <- sender_map %>%
  transmute(
    sender =
      as.character(
        sender
      ),
    receiver =
      as.character(
        receiver
      ),
    ligand =
      as.character(
        ligand
      )
  ) %>%
  distinct(
    sender,
    receiver,
    ligand
  )


activity_sender <- activity_corrected %>%
  select(
    -any_of(
      c(
        "nichenet_percentile"
      )
    )
  ) %>%
  inner_join(
    sender_map_unique,
    by = c(
      "receiver",
      "ligand"
    ),
    relationship =
      "many-to-many"
  ) %>%
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    .keep_all = TRUE
  )


sender_expansion_audit <- activity_sender %>%
  count(
    receiver,
    receiver_program,
    ligand,
    name =
      "n_sender_states"
  ) %>%
  arrange(
    desc(
      n_sender_states
    ),
    receiver,
    receiver_program,
    ligand
  )


write.csv(
  sender_expansion_audit,
  file.path(
    TAB_OUT,
    "03_sender_expansion_audit_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Add sender pseudobulk DE and sender expression
# ==============================================================================

sender_de_join <- sender_de %>%
  transmute(
    sender =
      as.character(
        sender
      ),
    ligand =
      as.character(
        gene
      ),
    sender_ligand_logFC =
      safe_num(
        logFC
      ),
    sender_ligand_logCPM =
      safe_num(
        logCPM
      ),
    sender_ligand_PValue =
      safe_num(
        PValue
      ),
    sender_ligand_FDR =
      safe_num(
        FDR
      )
  ) %>%
  distinct(
    sender,
    ligand,
    .keep_all = TRUE
  )


sender_expr_join <- expr_fraction %>%
  filter(
    grepl(
      "^Mphi_",
      celltype
    )
  ) %>%
  transmute(
    sender =
      as.character(
        celltype
      ),
    ligand =
      as.character(
        gene
      ),
    sender_ligand_pct_expressed =
      safe_num(
        pct_expressed
      )
  ) %>%
  distinct(
    sender,
    ligand,
    .keep_all = TRUE
  )


candidate_base <- activity_sender %>%
  left_join(
    sender_de_join,
    by = c(
      "sender",
      "ligand"
    ),
    relationship =
      "many-to-one"
  ) %>%
  left_join(
    sender_expr_join,
    by = c(
      "sender",
      "ligand"
    ),
    relationship =
      "many-to-one"
  )


# ==============================================================================
# 9. Rebuild ligand-level CellChat evidence from FULL v6.6.2.1 table
#
# Strategy:
#   - only Mphi -> Hep rows
#   - select LR rows that have:
#       total_support_samples >= 3
#       pairwise consistency >= 0.75
#   - among supported rows for each sender/receiver/ligand, select the row with
#     largest absolute supported Tx-Sham delta as the representative LR
#
# Also retain all receptor/pathway names for audit.
# ==============================================================================

cellchat_mh <- cellchat %>%
  filter(
    direction ==
      "Mphi_to_Hep"
  ) %>%
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
      ),
    pathway_name =
      as.character(
        pathway_name
      ),
    supported_Tx_minus_Sham =
      safe_num(
        supported_Tx_minus_Sham
      ),
    total_support_samples =
      safe_num(
        total_support_samples
      ),
    supported_pairwise_direction_consistency =
      safe_num(
        supported_pairwise_direction_consistency
      )
  )


cellchat_lr_supported <- cellchat_mh %>%
  filter(
    !is.na(
      ligand
    ),
    ligand !=
      "",
    total_support_samples >=
      CELLCHAT_MIN_SUPPORT_SAMPLES,
    supported_pairwise_direction_consistency >=
      CELLCHAT_MIN_PAIRWISE_CONSISTENCY,
    !is.na(
      supported_Tx_minus_Sham
    ),
    supported_Tx_minus_Sham !=
      0
  )


cellchat_best <- cellchat_lr_supported %>%
  arrange(
    desc(
      abs(
        supported_Tx_minus_Sham
      )
    ),
    desc(
      total_support_samples
    ),
    desc(
      supported_pairwise_direction_consistency
    )
  ) %>%
  group_by(
    source,
    target,
    ligand
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  transmute(
    sender =
      source,
    receiver =
      target,
    ligand,
    CellChat_present =
      TRUE,
    CellChat_best_supported_delta =
      supported_Tx_minus_Sham,
    CellChat_support_samples =
      total_support_samples,
    CellChat_pairwise_consistency =
      supported_pairwise_direction_consistency,
    CellChat_best_receptor =
      receptor,
    CellChat_best_pathway =
      pathway_name
  )


cellchat_all_names <- cellchat_mh %>%
  filter(
    !is.na(
      ligand
    ),
    ligand !=
      ""
  ) %>%
  group_by(
    source,
    target,
    ligand
  ) %>%
  summarise(
    CellChat_all_receptors =
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

    CellChat_all_pathways =
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


cellchat_ligand <- cellchat_all_names %>%
  left_join(
    cellchat_best,
    by = c(
      "sender",
      "receiver",
      "ligand"
    ),
    relationship =
      "one-to-one"
  ) %>%
  mutate(
    CellChat_present =
      replace_na(
        CellChat_present,
        FALSE
      )
  )


write.csv(
  cellchat_ligand,
  file.path(
    TAB_OUT,
    "04_CellChat_ligand_level_rebuilt_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Merge the three evidence layers
# ==============================================================================

candidate <- candidate_base %>%
  left_join(
    cellchat_ligand,
    by = c(
      "sender",
      "receiver",
      "ligand"
    ),
    relationship =
      "many-to-one"
  ) %>%
  mutate(
    CellChat_present =
      replace_na(
        CellChat_present,
        FALSE
      ),

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

    sender_direction_opposite =
      is.finite(
        sender_ligand_logFC
      ) &
      sign(
        sender_ligand_logFC
      ) ==
        -expected_sign,

    sender_DE_supported =
      sender_expressed &
      sender_direction_match &
      is.finite(
        sender_ligand_PValue
      ) &
      sender_ligand_PValue <=
        SENDER_P_THRESHOLD,

    sender_DE_opposite_supported =
      sender_expressed &
      sender_direction_opposite &
      is.finite(
        sender_ligand_PValue
      ) &
      sender_ligand_PValue <=
        SENDER_P_THRESHOLD,

    CellChat_direction =
      sign_direction(
        CellChat_best_supported_delta
      ),

    CellChat_direction_match =
      CellChat_present &
      is.finite(
        CellChat_best_supported_delta
      ) &
      sign(
        CellChat_best_supported_delta
      ) ==
        expected_sign,

    CellChat_direction_opposite =
      CellChat_present &
      is.finite(
        CellChat_best_supported_delta
      ) &
      sign(
        CellChat_best_supported_delta
      ) ==
        -expected_sign,

    CellChat_replicate_supported =
      CellChat_present &
      is.finite(
        CellChat_support_samples
      ) &
      CellChat_support_samples >=
        CELLCHAT_MIN_SUPPORT_SAMPLES &
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

    cross_layer_sender_CellChat_opposite =
      CellChat_replicate_supported &
      sender_expressed &
      is.finite(
        sender_ligand_logFC
      ) &
      is.finite(
        CellChat_best_supported_delta
      ) &
      sender_ligand_logFC !=
        0 &
      CellChat_best_supported_delta !=
        0 &
      sign(
        sender_ligand_logFC
      ) !=
        sign(
          CellChat_best_supported_delta
        ),

    directional_discordance =
      sender_DE_opposite_supported |
      CellChat_supported_opposite |
      cross_layer_sender_CellChat_opposite,

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

    evidence_class =
      case_when(
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
      mechanism_family(
        ligand
      ),

    key_axis =
      ligand %in%
        KEY_LIGANDS,

    secondary_axis =
      ligand %in%
        SECONDARY_LIGANDS,

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
          sender_DE_supported
        ) +
      5 *
        as.integer(
          CellChat_supported_matching
        ) +
      pmax(
        0,
        replace_na(
          nichenet_percentile_corrected,
          0
        )
      )
  ) %>%
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    .keep_all = TRUE
  ) %>%
  arrange(
    evidence_grade,
    desc(
      concordance_score
    ),
    nichenet_rank
  )


# ==============================================================================
# 11. Main concordance tables
# ==============================================================================

write.csv(
  candidate,
  file.path(
    TAB_OUT,
    "05_all_Mphi_Hep_ligand_3way_concordance_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


strict_supported <- candidate %>%
  filter(
    evidence_grade <=
      3,
    evidence_class !=
      "discordant"
  )


write.csv(
  strict_supported,
  file.path(
    TAB_OUT,
    "06_supported_candidate_shortlist_v6.6.3.1.csv"
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
    "07_three_way_supported_ligands_v6.6.3.1.csv"
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
    "08_directionally_discordant_ligands_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Primary Injury-Hepatocyte concordance
# ==============================================================================

injury_candidate <- candidate %>%
  filter(
    receiver ==
      PRIMARY_RECEIVER
  ) %>%
  arrange(
    receiver_program,
    evidence_grade,
    desc(
      concordance_score
    ),
    nichenet_rank
  )


write.csv(
  injury_candidate,
  file.path(
    TAB_OUT,
    "09_InjuryHep_all_ligand_concordance_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


injury_strict <- injury_candidate %>%
  filter(
    evidence_grade <=
      3,
    evidence_class !=
      "discordant"
  )


write.csv(
  injury_strict,
  file.path(
    TAB_OUT,
    "10_InjuryHep_strict_supported_candidates_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Direct key-axis comparison
# ==============================================================================

key_axis_status <- injury_candidate %>%
  filter(
    ligand %in%
      c(
        KEY_LIGANDS,
        SECONDARY_LIGANDS
      )
  ) %>%
  mutate(
    key_question =
      case_when(
        ligand ==
          "Pdgfb" ~
          "PDGFB_down_concordance",

        ligand ==
          "Sema4d" ~
          "SEMA4D_down_concordance",

        ligand ==
          "Fn1" ~
          "FN1_directional_discordance",

        ligand ==
          "Tnf" ~
          "TNF_loss_or_down_support",

        ligand ==
          "Plau" ~
          "PLAU_down_concordance",

        ligand ==
          "Spp1" ~
          "SPP1_secondary_audit",

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
    mechanism_family,
    evidence_class,
    strict_evidence_count,
    NicheNet_top20,
    nichenet_rank,
    nichenet_percentile_corrected,
    aupr_corrected,
    sender_expressed,
    sender_ligand_pct_expressed,
    sender_ligand_logFC,
    sender_ligand_PValue,
    sender_ligand_FDR,
    sender_direction,
    sender_DE_supported,
    CellChat_present,
    CellChat_best_supported_delta,
    CellChat_best_receptor,
    CellChat_best_pathway,
    CellChat_support_samples,
    CellChat_pairwise_consistency,
    CellChat_direction,
    CellChat_supported_matching,
    directional_discordance
  ) %>%
  arrange(
    ligand,
    evidence_grade =
      grade_class(
        evidence_class
      ),
    receiver_program,
    sender
  )


write.csv(
  key_axis_status,
  file.path(
    TAB_OUT,
    "11_InjuryHep_PDGFB_SEMA4D_FN1_TNF_PLAU_status_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Focus specifically on Repair/Resolution-Mphi
# ==============================================================================

repair_key <- key_axis_status %>%
  filter(
    sender ==
      "Mphi_Repair-Resolution"
  )


write.csv(
  repair_key,
  file.path(
    TAB_OUT,
    "12_RepairResolutionMphi_to_InjuryHep_key_axes_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Mechanistic interpretation matrix
#
# This is a machine-readable summary only. It does not claim causality.
# ==============================================================================

mechanistic_matrix <- key_axis_status %>%
  mutate(
    sender_change =
      sign_direction(
        sender_ligand_logFC
      ),

    communication_change =
      sign_direction(
        CellChat_best_supported_delta
      ),

    receiver_program_direction =
      ifelse(
        receiver_program ==
          "Tx_up_program",
        "Tx_up",
        "Tx_down"
      ),

    interpretation_class =
      case_when(
        evidence_class ==
          "3_way_supported" ~
          "NicheNet_sender_CellChat_concordant",

        evidence_class ==
          "CellChat_plus_senderDE" ~
          "Sender_and_CellChat_concordant_NicheNet_not_top20",

        evidence_class ==
          "NicheNet_plus_CellChat" ~
          "NicheNet_and_CellChat_concordant_senderDE_not_supported",

        evidence_class ==
          "NicheNet_plus_senderDE" ~
          "NicheNet_and_senderDE_concordant_CellChat_not_supported",

        evidence_class ==
          "discordant" ~
          "Directional_discordance",

        TRUE ~
          "Partial_or_insufficient_support"
      )
  ) %>%
  select(
    ligand,
    sender,
    receiver_program,
    NicheNet_rank =
      nichenet_rank,
    NicheNet_percentile =
      nichenet_percentile_corrected,
    sender_logFC =
      sender_ligand_logFC,
    sender_PValue =
      sender_ligand_PValue,
    CellChat_delta =
      CellChat_best_supported_delta,
    CellChat_support_samples,
    CellChat_pairwise_consistency,
    sender_change,
    communication_change,
    receiver_program_direction,
    evidence_class,
    interpretation_class
  )


write.csv(
  mechanistic_matrix,
  file.path(
    TAB_OUT,
    "13_InjuryHep_key_axis_mechanistic_matrix_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Rebuild ligand-receptor-target chains without inflated many-to-many joins
# ==============================================================================

target_links_clean <- target_links %>%
  transmute(
    receiver =
      as.character(
        receiver
      ),
    receiver_program =
      as.character(
        receiver_program
      ),
    ligand =
      as.character(
        ligand
      ),
    target =
      as.character(
        target
      ),
    regulatory_potential =
      safe_num(
        regulatory_potential
      ),
    receiver_logFC =
      if (
        "receiver_logFC" %in%
          colnames(
            target_links
          )
      ) {
        safe_num(
          receiver_logFC
        )
      } else {
        NA_real_
      },
    receiver_PValue =
      if (
        "receiver_PValue" %in%
          colnames(
            target_links
          )
      ) {
        safe_num(
          receiver_PValue
        )
      } else {
        NA_real_
      },
    receiver_FDR =
      if (
        "receiver_FDR" %in%
          colnames(
            target_links
          )
      ) {
        safe_num(
          receiver_FDR
        )
      } else {
        NA_real_
      }
  ) %>%
  distinct(
    receiver,
    receiver_program,
    ligand,
    target,
    .keep_all = TRUE
  )


candidate_for_target <- injury_candidate %>%
  select(
    sender,
    receiver,
    receiver_program,
    ligand,
    evidence_class,
    evidence_grade,
    strict_evidence_count,
    nichenet_rank,
    nichenet_percentile_corrected,
    sender_ligand_logFC,
    sender_ligand_PValue,
    CellChat_best_supported_delta,
    CellChat_best_receptor,
    CellChat_best_pathway,
    CellChat_support_samples,
    CellChat_pairwise_consistency,
    mechanism_family,
    concordance_score
  ) %>%
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    .keep_all = TRUE
  )


lrt <- candidate_for_target %>%
  inner_join(
    target_links_clean,
    by = c(
      "receiver",
      "receiver_program",
      "ligand"
    ),
    relationship =
      "many-to-many"
  ) %>%
  mutate(
    receptor =
      replace_na(
        CellChat_best_receptor,
        "No_replicate_supported_CellChat_receptor"
      ),

    LRT_score =
      concordance_score +
      10 *
        pmin(
          pmax(
            replace_na(
              regulatory_potential,
              0
            ),
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
  distinct(
    sender,
    receiver,
    receiver_program,
    ligand,
    receptor,
    target,
    .keep_all = TRUE
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


lrt_key <- lrt %>%
  filter(
    ligand %in%
      c(
        KEY_LIGANDS,
        SECONDARY_LIGANDS
      )
  )


write.csv(
  lrt_key,
  file.path(
    TAB_OUT,
    "14_InjuryHep_key_ligand_receptor_target_links_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


lrt_top30 <- lrt %>%
  filter(
    evidence_grade <=
      3,
    evidence_class !=
      "discordant"
  ) %>%
  slice_max(
    order_by =
      LRT_score,
    n =
      30,
    with_ties =
      FALSE
  )


write.csv(
  lrt_top30,
  file.path(
    TAB_OUT,
    "15_InjuryHep_top30_supported_LR_target_chains_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Evidence class counts
# ==============================================================================

class_counts <- injury_candidate %>%
  count(
    receiver_program,
    evidence_class,
    name =
      "n_candidates"
  ) %>%
  complete(
    receiver_program,
    evidence_class =
      CLASS_LEVELS,
    fill =
      list(
        n_candidates =
          0
      )
  )


write.csv(
  class_counts,
  file.path(
    TAB_OUT,
    "16_InjuryHep_evidence_class_counts_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Figure 1: corrected NicheNet percentile audit
# ==============================================================================

p_percentile <- ggplot(
  activity_corrected %>%
    filter(
      receiver ==
        PRIMARY_RECEIVER
    ),
  aes(
    x =
      nichenet_rank,
    y =
      nichenet_percentile_corrected
  )
) +
  geom_point(
    size =
      1
  ) +
  facet_wrap(
    ~ receiver_program,
    nrow =
      1
  ) +
  labs(
    title =
      "Corrected NicheNet percentile | Injury-Hepatocyte",
    subtitle =
      "Percentile now decreases monotonically with ligand rank",
    x =
      "NicheNet rank",
    y =
      "Corrected percentile"
  ) +
  theme_classic(
    base_size =
      8
  )


save_pdf(
  p_percentile,
  file.path(
    FIG_OUT,
    "01_corrected_NicheNet_percentile_v6.6.3.1.pdf"
  ),
  9,
  4.5
)


# ==============================================================================
# 19. Figure 2: strict supported Injury-Hep candidates
# ==============================================================================

strict_plot_df <- injury_candidate %>%
  filter(
    evidence_grade <=
      3,
    evidence_class !=
      "discordant"
  ) %>%
  group_by(
    receiver_program
  ) %>%
  arrange(
    evidence_grade,
    desc(
      concordance_score
    ),
    .by_group =
      TRUE
  ) %>%
  slice_head(
    n =
      25
  ) %>%
  ungroup() %>%
  mutate(
    label =
      paste0(
        sender,
        " | ",
        ligand
      )
  )


if (
  nrow(
    strict_plot_df
  )
) {

  p_strict <- ggplot(
    strict_plot_df,
    aes(
      x =
        receiver_program,
      y =
        label,
      size =
        nichenet_percentile_corrected,
      fill =
        strict_evidence_count
    )
  ) +
    geom_point(
      shape =
        21
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
        "Strict Mphi ligand -> Injury-Hep concordance",
      subtitle =
        "NicheNet + sender pseudobulk DE + replicate-supported CellChat",
      x =
        NULL,
      y =
        NULL,
      size =
        "Corrected\nNicheNet percentile",
      fill =
        "Matching\nevidence layers"
    ) +
    theme_classic(
      base_size =
        7
    )


  save_pdf(
    p_strict,
    file.path(
      FIG_OUT,
      "02_InjuryHep_strict_concordant_candidates_v6.6.3.1.pdf"
    ),
    10,
    max(
      8,
      min(
        18,
        5 +
          0.16 *
            n_distinct(
              strict_plot_df$label
            )
      )
    )
  )
}


# ==============================================================================
# 20. Figure 3: direct key-axis audit
# ==============================================================================

key_plot_df <- key_axis_status %>%
  filter(
    ligand %in%
      KEY_LIGANDS
  ) %>%
  mutate(
    sender_short =
      gsub(
        "^Mphi_",
        "",
        sender
      ),

    plot_x =
      replace_na(
        sender_ligand_logFC,
        0
      ),

    plot_y =
      replace_na(
        CellChat_best_supported_delta,
        0
      )
  )


if (
  nrow(
    key_plot_df
  )
) {

  p_key <- ggplot(
    key_plot_df,
    aes(
      x =
        plot_x,
      y =
        plot_y,
      size =
        nichenet_percentile_corrected,
      fill =
        evidence_class
    )
  ) +
    geom_hline(
      yintercept =
        0,
      linewidth =
        0.3
    ) +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.3
    ) +
    geom_point(
      shape =
        21
    ) +
    facet_grid(
      receiver_program ~ ligand,
      scales =
        "free"
    ) +
    labs(
      title =
        "PDGFB / SEMA4D / FN1 / TNF / PLAU concordance audit",
      subtitle =
        "x = Mphi sender log2FC Tx/Sham; y = replicate-supported CellChat Tx-Sham",
      x =
        "Sender ligand log2FC Tx vs Sham",
      y =
        "CellChat supported Tx - Sham",
      size =
        "Corrected\nNicheNet percentile",
      fill =
        "Evidence class"
    ) +
    theme_classic(
      base_size =
        7
    )


  save_pdf(
    p_key,
    file.path(
      FIG_OUT,
      "03_PDGFB_SEMA4D_FN1_TNF_PLAU_concordance_audit_v6.6.3.1.pdf"
    ),
    16,
    8
  )
}


# ==============================================================================
# 21. Figure 4: Repair/Resolution-Mphi key axes
# ==============================================================================

repair_plot_df <- repair_key %>%
  filter(
    ligand %in%
      KEY_LIGANDS
  ) %>%
  select(
    receiver_program,
    ligand,
    sender_ligand_logFC,
    CellChat_best_supported_delta,
    nichenet_rank,
    nichenet_percentile_corrected,
    evidence_class
  ) %>%
  pivot_longer(
    cols = c(
      sender_ligand_logFC,
      CellChat_best_supported_delta
    ),
    names_to =
      "evidence_layer",
    values_to =
      "effect"
  ) %>%
  mutate(
    evidence_layer =
      recode(
        evidence_layer,
        sender_ligand_logFC =
          "Sender pseudobulk log2FC",
        CellChat_best_supported_delta =
          "CellChat supported Tx-Sham"
      )
  )


if (
  nrow(
    repair_plot_df
  )
) {

  p_repair <- ggplot(
    repair_plot_df,
    aes(
      x =
        ligand,
      y =
        effect
    )
  ) +
    geom_hline(
      yintercept =
        0,
      linewidth =
        0.3
    ) +
    geom_point(
      size =
        2
    ) +
    facet_grid(
      receiver_program ~ evidence_layer,
      scales =
        "free_y"
    ) +
    labs(
      title =
        "Repair/Resolution-Mphi -> Injury-Hep key-axis direction",
      subtitle =
        "Direct comparison of sender expression and supported CellChat direction",
      x =
        NULL,
      y =
        "Effect"
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
    p_repair,
    file.path(
      FIG_OUT,
      "04_RepairResolutionMphi_key_axis_direction_v6.6.3.1.pdf"
    ),
    11,
    7
  )
}


# ==============================================================================
# 22. Figure 5: top L-R-target chains
# ==============================================================================

if (
  nrow(
    lrt_top30
  )
) {

  lrt_plot <- lrt_top30 %>%
    slice_head(
      n =
        20
    ) %>%
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


  p_lrt <- ggplot(
    lrt_plot,
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
        evidence_class
    )
  ) +
    geom_point(
      shape =
        21
    ) +
    labs(
      title =
        "Top supported Mphi -> Injury-Hep ligand-receptor-target chains",
      subtitle =
        "Receiver target links from NicheNet; receptor from replicate-supported CellChat",
      x =
        "NicheNet regulatory potential",
      y =
        NULL,
      size =
        "|Receiver log2FC|",
      fill =
        "Evidence class"
    ) +
    theme_classic(
      base_size =
        7
    )


  save_pdf(
    p_lrt,
    file.path(
      FIG_OUT,
      "05_top20_InjuryHep_LR_target_chains_v6.6.3.1.pdf"
    ),
    12,
    10
  )
}


# ==============================================================================
# 23. Save post-processing result RDS
# ==============================================================================

out_object <- list(
  version =
    "v6.6.3.1",

  activity_corrected =
    activity_corrected,

  percentile_audit =
    percentile_audit,

  sender_expansion_audit =
    sender_expansion_audit,

  cellchat_ligand =
    cellchat_ligand,

  candidate =
    candidate,

  injury_candidate =
    injury_candidate,

  injury_strict =
    injury_strict,

  key_axis_status =
    key_axis_status,

  repair_key =
    repair_key,

  mechanistic_matrix =
    mechanistic_matrix,

  lrt_key =
    lrt_key,

  lrt_top30 =
    lrt_top30
)


saveRDS(
  out_object,
  file.path(
    RDS_OUT,
    "Mouse_MASH_Mphi5_Hep5_NicheNet_concordance_v6.6.3.1.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 24. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "source_NicheNet_result_RDS",
    "source_full_refined_CellChat",
    "NicheNet_rerun",
    "edgeR_rerun",
    "CellChat_rerun",
    "primary_receiver",
    "NicheNet_topN",
    "sender_expression_threshold",
    "sender_nominal_P_threshold",
    "CellChat_min_support_samples",
    "CellChat_min_pairwise_consistency",
    "key_axes",
    "biological_replicates",
    "percentile_fix"
  ),
  value = c(
    "v6.6.3.1",
    RESULT_RDS,
    CELLCHAT_FILE,
    "FALSE",
    "FALSE",
    "FALSE",
    PRIMARY_RECEIVER,
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
      CELLCHAT_MIN_SUPPORT_SAMPLES
    ),
    as.character(
      CELLCHAT_MIN_PAIRWISE_CONSISTENCY
    ),
    paste(
      KEY_LIGANDS,
      collapse = " | "
    ),
    "Sham n=2; Tx n=2",
    "Corrected within receiver x receiver_program using rank and candidate count"
  )
)


write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.3.1.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.3.1.txt"
    )
)


# ==============================================================================
# 25. Console summary
# ==============================================================================

msg(
  "NicheNet percentile correction audit:"
)


print(
  percentile_audit %>%
    filter(
      receiver ==
        PRIMARY_RECEIVER
    )
)


msg(
  "Three-way supported Injury-Hep candidates:"
)


print(
  injury_candidate %>%
    filter(
      evidence_class ==
        "3_way_supported"
    ) %>%
    select(
      sender,
      receiver_program,
      ligand,
      mechanism_family,
      nichenet_rank,
      nichenet_percentile_corrected,
      sender_ligand_logFC,
      sender_ligand_PValue,
      CellChat_best_supported_delta,
      CellChat_support_samples,
      CellChat_pairwise_consistency,
      evidence_class
    )
)


msg(
  "Key-axis status:"
)


print(
  key_axis_status %>%
    select(
      sender,
      receiver_program,
      ligand,
      evidence_class,
      nichenet_rank,
      nichenet_percentile_corrected,
      sender_ligand_logFC,
      sender_ligand_PValue,
      CellChat_best_supported_delta,
      CellChat_support_samples,
      CellChat_pairwise_consistency,
      directional_discordance
    )
)


msg(
  "DONE."
)


msg(
  "No NicheNet, edgeR, or CellChat computation was rerun."
)


msg(
  "Output directory: ",
  OUT
)
