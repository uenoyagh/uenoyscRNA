#!/usr/bin/env Rscript

# ============================================================
# Human-Mouse MASH Mphi -> HSC 4-axis cross-species hierarchy v1.0
#
# PURPOSE
#   Integrate the canonical Human Mphi->HSC four-axis hierarchy
#   with the canonical Mouse Sham->Tx Mphi/HSC response.
#
# HUMAN AXES
#   1. SPP1 -> HSC integrins
#   2. CDCP1-associated -> HSC EGFR-supported axis
#   3. TGFB1 -> TGF-beta receptors
#   4. PDGFA -> PDGFRA
#
# CROSS-SPECIES PRINCIPLE
#   Compare FUNCTIONAL AXES, not only one-to-one gene identity.
#   In particular, Human PDGFA-PDGFRA and Mouse Pdgfb-Pdgfrb
#   are evaluated as possible context/species-specific usage
#   within a conserved PDGF -> HSC axis.
#
# INFERENCE BOUNDARIES
#   Human: n=8, exploratory sample-level.
#   Mouse: n=2 Sham + n=2 Tx; treatment comparisons are descriptive.
#   Mouse n=4 correlations are descriptive only.
#   Mouse v6.4.4/v6.4.5 host-weighted HSC burden is used as the
#   canonical receiver treatment endpoint.
#
# IMPORTANT Cdcp1 BOUNDARY
#   v6.4.10.4 provenance information is NOT used quantitatively
#   for cross-species scoring because the full-library CPM
#   denominator correction remains pending. Cdcp1 source
#   provenance therefore remains an explicit mechanistic caveat.
#
# No RDS is modified. No reclustering.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

# ------------------------------------------------------------
# Canonical Human checkpoint
# ------------------------------------------------------------

HUMAN_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1"
)

HUMAN_HIERARCHY <- file.path(
  HUMAN_ROOT,
  "tables",
  "02_PRIMARY_Human_Mphi_HSC_4Axis_mechanistic_hierarchy.tsv"
)

# ------------------------------------------------------------
# Canonical Mouse focused Mphi-HSC closure
# ------------------------------------------------------------

MOUSE_CLOSURE_ROOT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_HSC_IntegratedClosure_v6.4.5"
)

MOUSE_TOTAL_MATRIX <- file.path(
  MOUSE_CLOSURE_ROOT,
  "tables",
  "10_TOTAL_MYELOID_focused_mechanistic_matrix.tsv"
)

MOUSE_TREM2_MATRIX <- file.path(
  MOUSE_CLOSURE_ROOT,
  "tables",
  "12_HUMAN_PRIOR_TREM2_axis_summary.tsv"
)

MOUSE_RECEIVER <- file.path(
  MOUSE_CLOSURE_ROOT,
  "tables",
  "06_receiver_endpoints_Sham_vs_Tx.tsv"
)

# ------------------------------------------------------------
# Canonical Mouse raw-count Cdcp1 support
# ------------------------------------------------------------

MOUSE_PB_ROOT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1"
)

MOUSE_CDCP1_RANK <- file.path(
  MOUSE_PB_ROOT,
  "tables",
  "08_PRIMARY_Cdcp1_transcriptome_rank_summary.tsv"
)

# ------------------------------------------------------------
# Human-Mouse PDGF analysis
# ------------------------------------------------------------

PDGF_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2"
)

MOUSE_PDGF_EFFECT <- file.path(
  PDGF_ROOT,
  "tables",
  "13_Mouse_PDGF_Sham_vs_Tx_effects.tsv"
)

# ------------------------------------------------------------
# Optional Mouse provenance boundary
# ------------------------------------------------------------

MOUSE_PROV_ROOT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_WholeHost_Provenance_StressAudit_v6.4.10.4"
)

MOUSE_PROV_TABLE <- file.path(
  MOUSE_PROV_ROOT,
  "tables",
  "05_PRIMARY_focus_gene_provenance_and_TREM2_rank.tsv"
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_DOC <- file.path(OUT, "docs")
DIR_PROV <- file.path(OUT, "provenance")
DIR_LOG <- file.path(OUT, "logs")
DIR_OBJ <- file.path(OUT, "objects")

for (d in c(DIR_TABLE, DIR_DOC, DIR_PROV, DIR_LOG, DIR_OBJ)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(
  DIR_LOG,
  "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human-Mouse MASH Mphi-HSC 4-axis cross-species hierarchy v1.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Input checks
# ------------------------------------------------------------

required_files <- c(
  HUMAN_HIERARCHY,
  MOUSE_TOTAL_MATRIX,
  MOUSE_TREM2_MATRIX,
  MOUSE_RECEIVER,
  MOUSE_CDCP1_RANK,
  MOUSE_PDGF_EFFECT
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Human hierarchy:\n", HUMAN_HIERARCHY, "\n\n")
cat("Mouse total-myeloid matrix:\n", MOUSE_TOTAL_MATRIX, "\n\n")
cat("Mouse TREM2 matrix:\n", MOUSE_TREM2_MATRIX, "\n\n")
cat("Mouse receiver endpoints:\n", MOUSE_RECEIVER, "\n\n")
cat("Mouse Cdcp1 raw-count support:\n", MOUSE_CDCP1_RANK, "\n\n")
cat("Mouse PDGF effects:\n", MOUSE_PDGF_EFFECT, "\n\n")

# ------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------

human <- read_tsv(
  HUMAN_HIERARCHY,
  show_col_types = FALSE
)

mouse_total_matrix <- read_tsv(
  MOUSE_TOTAL_MATRIX,
  show_col_types = FALSE
)

mouse_trem2_matrix <- read_tsv(
  MOUSE_TREM2_MATRIX,
  show_col_types = FALSE
)

mouse_receiver <- read_tsv(
  MOUSE_RECEIVER,
  show_col_types = FALSE
)

mouse_cdcp1_rank <- read_tsv(
  MOUSE_CDCP1_RANK,
  show_col_types = FALSE
)

mouse_pdgf <- read_tsv(
  MOUSE_PDGF_EFFECT,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Human integrity
# ------------------------------------------------------------

expected_human <- c(
  "SPP1",
  "CDCP1",
  "TGFB1",
  "PDGFA"
)

if (!all(expected_human %in% human$ligand)) {
  stop(
    "Human canonical hierarchy does not contain all four expected ligands."
  )
}

if (nrow(human) != 4) {
  stop("Human canonical hierarchy should contain exactly four rows.")
}

human_core <- human %>%
  select(
    hierarchy_order,
    hierarchy_tier,
    axis,
    ligand,
    mechanistic_role,
    sender_scope,
    primary_receiver_endpoint,
    n,
    rho,
    p_value_primary,
    FDR_primary,
    partial_rho_NAS,
    partial_p_NAS,
    receiver_support,
    NicheNet_rank,
    overall_evidence_grade,
    safe_statement,
    do_not_claim
  ) %>%
  arrange(hierarchy_order)

write_tsv(
  human_core,
  file.path(
    DIR_TABLE,
    "01_Human_canonical_4Axis_input.tsv"
  )
)

# ------------------------------------------------------------
# Mouse canonical receiver treatment effect
# Use host-weighted v6.4.4/v6.4.5 endpoints, NOT the simple
# per-cell HSC intensity from the PDGF exploratory script.
# ------------------------------------------------------------

needed_receiver <- c(
  "ECM_activated_HSC_Core_collagen",
  "ECM_activated_HSC_ECM_remodeling",
  "HSC_abundance",
  "WholeHSC_Core_collagen",
  "WholeHSC_ECM_remodeling"
)

if (!all(needed_receiver %in% mouse_receiver$receiver_endpoint)) {
  stop(
    "Mouse receiver table lacks required canonical endpoints."
  )
}

receiver_focus <- mouse_receiver %>%
  filter(
    receiver_endpoint %in% needed_receiver
  ) %>%
  select(
    receiver_endpoint,
    n_samples_Sham,
    n_samples_Tx,
    mean_value_Sham,
    mean_value_Tx,
    Tx_minus_Sham,
    Tx_over_Sham,
    observed_direction
  )

write_tsv(
  receiver_focus,
  file.path(
    DIR_TABLE,
    "02_Mouse_canonical_HSC_receiver_Tx_effects.tsv"
  )
)

core_receiver <- receiver_focus %>%
  filter(
    receiver_endpoint == "ECM_activated_HSC_Core_collagen"
  )

if (nrow(core_receiver) != 1) {
  stop("Expected one ECM_activated_HSC_Core_collagen row.")
}

# ------------------------------------------------------------
# Mouse Spp1 / Cdcp1 / Tgfb1 sender effects
# Deduplicate total-myeloid matrix because sender values repeat
# across receiver endpoints.
# ------------------------------------------------------------

needed_total_cols <- c(
  "ligand",
  "sender_Sham",
  "sender_Tx",
  "sender_Tx_minus_Sham",
  "sender_Tx_over_Sham",
  "sender_direction"
)

if (!all(needed_total_cols %in% colnames(mouse_total_matrix))) {
  stop(
    "Mouse total-myeloid matrix lacks required sender columns."
  )
}

mouse_total_sender <- mouse_total_matrix %>%
  filter(
    ligand %in% c("Spp1", "Cdcp1", "Tgfb1")
  ) %>%
  select(
    all_of(needed_total_cols)
  ) %>%
  distinct() %>%
  arrange(ligand)

if (nrow(mouse_total_sender) != 3) {
  stop(
    "Expected exactly one unique total-myeloid sender row for Spp1, Cdcp1, Tgfb1."
  )
}

# TREM2 values also repeat across receiver endpoints.
needed_trem2_cols <- c(
  "ligand",
  "sender_state",
  "sender_Sham",
  "sender_Tx",
  "sender_Tx_minus_Sham",
  "sender_Tx_over_Sham",
  "sender_direction"
)

if (!all(needed_trem2_cols %in% colnames(mouse_trem2_matrix))) {
  stop(
    "Mouse TREM2 summary lacks required sender columns."
  )
}

mouse_trem2_sender <- mouse_trem2_matrix %>%
  filter(
    ligand %in% c("Spp1", "Cdcp1", "Tgfb1"),
    sender_state == "Lipid-associated/TREM2-Mphi"
  ) %>%
  select(
    all_of(needed_trem2_cols)
  ) %>%
  distinct() %>%
  arrange(ligand)

if (nrow(mouse_trem2_sender) != 3) {
  stop(
    "Expected exactly one unique TREM2 sender row for Spp1, Cdcp1, Tgfb1."
  )
}

mouse_sender <- mouse_total_sender %>%
  rename(
    total_Sham = sender_Sham,
    total_Tx = sender_Tx,
    total_Tx_minus_Sham = sender_Tx_minus_Sham,
    total_Tx_over_Sham = sender_Tx_over_Sham,
    total_direction = sender_direction
  ) %>%
  left_join(
    mouse_trem2_sender %>%
      select(
        ligand,
        TREM2_Sham = sender_Sham,
        TREM2_Tx = sender_Tx,
        TREM2_Tx_minus_Sham = sender_Tx_minus_Sham,
        TREM2_Tx_over_Sham = sender_Tx_over_Sham,
        TREM2_direction = sender_direction
      ),
    by = "ligand"
  )

write_tsv(
  mouse_sender,
  file.path(
    DIR_TABLE,
    "03_Mouse_Spp1_Cdcp1_Tgfb1_sender_Tx_effects.tsv"
  )
)

# ------------------------------------------------------------
# Mouse Cdcp1 independent raw-count support
# ------------------------------------------------------------

cdcp1_trem2_pb <- mouse_cdcp1_rank %>%
  filter(
    gene == "Cdcp1",
    sender_state == "Lipid-associated/TREM2-Mphi"
  )

if (nrow(cdcp1_trem2_pb) != 1) {
  stop(
    "Expected one TREM2 Cdcp1 raw-count pseudobulk rank row."
  )
}

cdcp1_pb_support <- cdcp1_trem2_pb %>%
  transmute(
    gene,
    sender_state,
    total_counts,
    retained,
    delta_log2CPM,
    rank_down,
    percentile_down,
    support_statement =
      "Independent raw-count pseudobulk supports TREM2-state Cdcp1 downregulation; Mouse n=2/group remains descriptive."
  )

write_tsv(
  cdcp1_pb_support,
  file.path(
    DIR_TABLE,
    "04_Mouse_Cdcp1_rawCount_support.tsv"
  )
)

# ------------------------------------------------------------
# Optional Cdcp1 provenance boundary
# v6.4.10.4 is NOT used for quantitative scoring.
# ------------------------------------------------------------

if (file.exists(MOUSE_PROV_TABLE)) {

  prov_tmp <- read_tsv(
    MOUSE_PROV_TABLE,
    show_col_types = FALSE
  )

  cdcp1_prov <- prov_tmp %>%
    filter(gene == "Cdcp1") %>%
    mutate(
      quantitative_use_in_this_analysis = FALSE,
      denominator_boundary =
        "v6.4.10.4 full-library CPM denominator correction pending; provenance retained as qualitative boundary only."
    )

} else {

  cdcp1_prov <- tibble(
    gene = "Cdcp1",
    provenance_class = NA_character_,
    quantitative_use_in_this_analysis = FALSE,
    denominator_boundary =
      "v6.4.10.4 provenance table not found; full-library CPM denominator correction remains a required boundary before quantitative provenance use."
  )
}

write_tsv(
  cdcp1_prov,
  file.path(
    DIR_TABLE,
    "05_Mouse_Cdcp1_provenance_boundary.tsv"
  )
)

# ------------------------------------------------------------
# Mouse PDGF treatment effects
# ------------------------------------------------------------

needed_pdgf_metrics <- c(
  "total_Pdgfa",
  "total_Pdgfb",
  "TREM2_Pdgfa",
  "TREM2_Pdgfb",
  "HSC_Pdgfra",
  "HSC_Pdgfrb",
  "total_Pdgfa_Pdgfra_product",
  "total_Pdgfb_Pdgfrb_product",
  "TREM2_Pdgfa_Pdgfra_product",
  "TREM2_Pdgfb_Pdgfrb_product"
)

if (!all(needed_pdgf_metrics %in% mouse_pdgf$metric)) {
  stop(
    "Mouse PDGF effect table lacks required metrics: ",
    paste(
      setdiff(
        needed_pdgf_metrics,
        mouse_pdgf$metric
      ),
      collapse = ";"
    )
  )
}

mouse_pdgf_focus <- mouse_pdgf %>%
  filter(
    metric %in% needed_pdgf_metrics
  ) %>%
  select(
    metric,
    mean_value_Sham,
    mean_value_Tx,
    n_samples_Sham,
    n_samples_Tx,
    Tx_minus_Sham,
    Tx_over_Sham,
    inference
  ) %>%
  arrange(
    match(metric, needed_pdgf_metrics)
  )

write_tsv(
  mouse_pdgf_focus,
  file.path(
    DIR_TABLE,
    "06_Mouse_PDGF_family_Tx_effects.tsv"
  )
)

get_metric_ratio <- function(metric_name) {
  x <- mouse_pdgf_focus %>%
    filter(metric == metric_name)

  if (nrow(x) != 1) {
    stop("Metric not uniquely resolved: ", metric_name)
  }

  x$Tx_over_Sham[[1]]
}

# ------------------------------------------------------------
# Build axis-specific Mouse evidence
# ------------------------------------------------------------

sender_lookup <- function(mouse_gene) {
  x <- mouse_sender %>%
    filter(ligand == mouse_gene)

  if (nrow(x) != 1) {
    stop("Mouse sender not uniquely resolved: ", mouse_gene)
  }

  x
}

spp1_m <- sender_lookup("Spp1")
cdcp1_m <- sender_lookup("Cdcp1")
tgfb1_m <- sender_lookup("Tgfb1")

mouse_axis <- tribble(
  ~human_ligand, ~mouse_counterpart, ~mouse_functional_axis,
  ~mouse_total_sender_ratio, ~mouse_TREM2_sender_ratio,
  ~mouse_primary_LR_ratio, ~mouse_secondary_LR_ratio,
  ~mouse_receiver_endpoint, ~mouse_receiver_ratio,
  ~mouse_treatment_pattern,

  "SPP1", "Spp1",
  "Spp1 -> HSC integrin/fibrotic-niche axis",
  spp1_m$total_Tx_over_Sham[[1]],
  spp1_m$TREM2_Tx_over_Sham[[1]],
  NA_real_, NA_real_,
  core_receiver$receiver_endpoint[[1]],
  core_receiver$Tx_over_Sham[[1]],
  "Spp1 burden is not attenuated by Tx despite reduced host-weighted ECM-HSC core-collagen output.",

  "CDCP1", "Cdcp1",
  "Cdcp1-associated HSC axis",
  cdcp1_m$total_Tx_over_Sham[[1]],
  cdcp1_m$TREM2_Tx_over_Sham[[1]],
  NA_real_, NA_real_,
  core_receiver$receiver_endpoint[[1]],
  core_receiver$Tx_over_Sham[[1]],
  "Cdcp1 burden decreases with Tx, including TREM2-associated Mphi, in the same direction as reduced host-weighted ECM-HSC core-collagen output; source provenance remains unresolved.",

  "TGFB1", "Tgfb1",
  "Tgfb1 -> TGF-beta receptor axis",
  tgfb1_m$total_Tx_over_Sham[[1]],
  tgfb1_m$TREM2_Tx_over_Sham[[1]],
  NA_real_, NA_real_,
  core_receiver$receiver_endpoint[[1]],
  core_receiver$Tx_over_Sham[[1]],
  "Tgfb1 burden is not attenuated by Tx despite reduced host-weighted ECM-HSC core-collagen output.",

  "PDGFA", "Pdgfa/Pdgfb",
  "PDGF-family -> HSC PDGF receptor axis",
  get_metric_ratio("total_Pdgfa"),
  get_metric_ratio("TREM2_Pdgfa"),
  get_metric_ratio("total_Pdgfa_Pdgfra_product"),
  get_metric_ratio("total_Pdgfb_Pdgfrb_product"),
  core_receiver$receiver_endpoint[[1]],
  core_receiver$Tx_over_Sham[[1]],
  "Mouse Pdgfa-Pdgfra is approximately unchanged, whereas Pdgfb-Pdgfrb is attenuated with Tx, supporting context-specific paralog usage within a conserved PDGF-family HSC axis."
)

# Add TREM2 Pdgfb ratio separately for transparency.
mouse_axis <- mouse_axis %>%
  mutate(
    mouse_TREM2_Pdgfb_ratio =
      ifelse(
        human_ligand == "PDGFA",
        get_metric_ratio(
          "TREM2_Pdgfb_Pdgfrb_product"
        ),
        NA_real_
      )
  )

write_tsv(
  mouse_axis,
  file.path(
    DIR_TABLE,
    "07_Mouse_4Axis_functional_mapping.tsv"
  )
)

# ------------------------------------------------------------
# Cross-species adjudication
#
# IMPORTANT:
#   "conserved" here means preservation of a FUNCTIONAL signaling
#   concept / receiver architecture. It does not mean identical
#   treatment direction or identical ligand paralog usage.
# ------------------------------------------------------------

cross_meta <- tribble(
  ~human_ligand, ~cross_species_class, ~bridge_grade,
  ~conserved_component, ~divergent_component,
  ~cross_species_safe_statement, ~cross_species_do_not_claim,

  "SPP1",
  "conserved_axis_candidate_but_Tx_response_not_attenuated",
  "MODERATE",
  "SPP1/integrin fibrotic-niche biology and HSC receiver architecture are compatible across species.",
  "Mouse Tx does not reduce total or TREM2-associated Spp1 burden while host-weighted ECM-HSC core-collagen output decreases.",
  "SPP1-integrin is a conserved profibrotic candidate axis, but the Mouse therapeutic response is not explained by Spp1 attenuation.",
  "Do not claim that Spp1 suppression mediates the Mouse Tx anti-fibrotic effect.",

  "CDCP1",
  "directionally_concordant_bridge_with_mouse_source_uncertainty",
  "MODERATE-HIGH",
  "Human CDCP1-associated HSC signal and Mouse Tx-associated Cdcp1 reduction converge directionally with reduced fibrogenic HSC burden.",
  "Mouse Cdcp1 is sparse and much more abundant outside Mphi; ambient/global contribution is not excluded.",
  "CDCP1/Cdcp1 is the strongest directionally concordant Human-Mouse candidate bridge, but Mouse cellular source attribution remains unresolved.",
  "Do not claim a proven Mphi-specific Cdcp1 mechanism or direct CDCP1-EGFR causality.",

  "TGFB1",
  "canonical_pathway_conserved_but_Tx_response_not_attenuated",
  "MODERATE",
  "Canonical TGF-beta profibrotic receiver biology is conserved.",
  "Mouse Tx does not reduce total or TREM2-associated Tgfb1 burden while fibrogenic HSC burden decreases.",
  "TGFB1/TGF-beta signaling is a conserved auxiliary profibrotic pathway, but Tgfb1 attenuation is not supported as the principal Mouse Tx mechanism.",
  "Do not claim that reduced Tgfb1 mediates the Mouse Tx effect or that TREM2 Mphi is the dominant absolute Tgfb1 source.",

  "PDGFA",
  "functional_axis_conserved_with_paralog_usage_divergence",
  "MODERATE-HIGH",
  "PDGF-family signaling to HSC is supported in both species.",
  "Human signal is PDGFA-PDGFRA/ECM-high-associated, whereas Mouse Tx preferentially attenuates Pdgfb-Pdgfrb rather than Pdgfa-Pdgfra.",
  "The PDGF-to-HSC functional axis is cross-species conserved with context/species-specific ligand-receptor paralog usage.",
  "Do not require one-to-one PDGFA-PDGFRA conservation; do not claim that Human PDGFA and Mouse Pdgfb are biologically identical signals."
)

cross <- human_core %>%
  left_join(
    mouse_axis,
    by = c("ligand" = "human_ligand")
  ) %>%
  left_join(
    cross_meta,
    by = c("ligand" = "human_ligand")
  ) %>%
  mutate(
    human_inference = "exploratory_sample_level_n8",
    mouse_inference = "descriptive_only_n2_per_group",
    mouse_common_receiver_reduction =
      mouse_receiver_ratio < 1,
    mouse_sender_attenuated =
      case_when(
        ligand == "PDGFA" ~
          mouse_secondary_LR_ratio < 1,
        TRUE ~
          mouse_total_sender_ratio < 1
      )
  ) %>%
  select(
    human_priority = hierarchy_order,
    human_tier = hierarchy_tier,
    human_axis = axis,
    human_ligand = ligand,
    human_endpoint = primary_receiver_endpoint,
    human_rho = rho,
    human_partial_rho_NAS = partial_rho_NAS,
    human_evidence_grade = overall_evidence_grade,
    mouse_counterpart,
    mouse_functional_axis,
    mouse_total_sender_ratio,
    mouse_TREM2_sender_ratio,
    mouse_primary_LR_ratio,
    mouse_secondary_LR_ratio,
    mouse_TREM2_Pdgfb_ratio,
    mouse_receiver_endpoint,
    mouse_receiver_ratio,
    mouse_sender_attenuated,
    mouse_common_receiver_reduction,
    mouse_treatment_pattern,
    cross_species_class,
    bridge_grade,
    conserved_component,
    divergent_component,
    cross_species_safe_statement,
    cross_species_do_not_claim,
    human_inference,
    mouse_inference
  ) %>%
  arrange(human_priority)

write_tsv(
  cross,
  file.path(
    DIR_TABLE,
    "08_PRIMARY_Human_Mouse_4Axis_cross_species_hierarchy.tsv"
  )
)

# ------------------------------------------------------------
# Concise manuscript / slide table
# ------------------------------------------------------------

concise <- cross %>%
  transmute(
    Human_priority = human_priority,
    Functional_axis = case_when(
      human_ligand == "SPP1" ~ "SPP1/Spp1-integrin",
      human_ligand == "CDCP1" ~ "CDCP1/Cdcp1-associated",
      human_ligand == "TGFB1" ~ "TGFB1/Tgfb1-TGF-betaR",
      human_ligand == "PDGFA" ~ "PDGF-family -> HSC PDGFR",
      TRUE ~ human_axis
    ),
    Human_evidence = human_evidence_grade,
    Human_rho = round(human_rho, 3),
    Human_NAS_adjusted_rho =
      round(human_partial_rho_NAS, 3),
    Mouse_Tx_axis_ratio = case_when(
      human_ligand == "PDGFA" ~
        mouse_secondary_LR_ratio,
      TRUE ~
        mouse_total_sender_ratio
    ),
    Mouse_ECM_HSC_coreCollagen_ratio =
      mouse_receiver_ratio,
    Cross_species_class = cross_species_class,
    Bridge_grade = bridge_grade,
    Interpretation =
      cross_species_safe_statement
  )

write_tsv(
  concise,
  file.path(
    DIR_TABLE,
    "09_CONCISE_Human_Mouse_4Axis_cross_species_hierarchy.tsv"
  )
)

# ------------------------------------------------------------
# Evidence-boundary table
# ------------------------------------------------------------

boundaries <- tribble(
  ~topic, ~boundary,

  "Human inference",
  "Human associations are exploratory at n=8; correlations and ligand-target predictions do not prove causality.",

  "Mouse inference",
  "Mouse Sham vs Tx comparisons are descriptive only (n=2/group); n=4 sample correlations are descriptive only.",

  "Mouse HSC receiver",
  "Canonical treatment-response endpoint is host-weighted ECM-activated HSC core-collagen burden from v6.4.4/v6.4.5, not simple per-cell intensity.",

  "Cdcp1 provenance",
  "Mouse Cdcp1 shows a reproducible TREM2-state decrease, but is sparse and non-Mphi enriched; endogenous Mphi signal is plausible while ambient/global contribution is not excluded.",

  "v6.4.10.4 denominator",
  "Whole-host provenance v6.4.10.4 is retained as a qualitative boundary only until the full-library CPM denominator correction is completed.",

  "SPP1",
  "Conserved profibrotic candidate biology does not imply that Mouse Tx acts through Spp1 suppression.",

  "TGFB1",
  "Conserved canonical TGF-beta biology does not imply that Mouse Tx acts through Tgfb1 suppression.",

  "PDGF",
  "Cross-species conservation is adjudicated at the PDGF-family -> HSC functional-axis level; paralog identity may differ between Human and Mouse/context."
)

write_tsv(
  boundaries,
  file.path(
    DIR_TABLE,
    "10_cross_species_interpretation_boundaries.tsv"
  )
)

# ------------------------------------------------------------
# Cross-species hierarchy checkpoint draft
# ------------------------------------------------------------

md <- c(
  "# Human-Mouse MASH MΦ→HSC 4-axis cross-species hierarchy v1.0",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Core principle",
  "",
  "Cross-species conservation is evaluated at the functional signaling-axis level, not by requiring identical ligand paralogs or identical treatment-direction behavior.",
  "",
  "## Integrated hierarchy",
  "",
  "### 1. SPP1/Spp1-integrin",
  "- Human: primary mechanistic axis.",
  "- Mouse: Spp1 burden is not attenuated by Tx despite lower host-weighted ECM-HSC core-collagen burden.",
  "- Classification: **conserved profibrotic candidate biology, but not a Tx-response attenuation mechanism**.",
  "",
  "### 2. CDCP1/Cdcp1-associated",
  "- Human: strongest novel TREM2-linked candidate.",
  "- Mouse: total and TREM2-associated Cdcp1 burden decrease with Tx; independent raw-count pseudobulk supports the TREM2-state decrease.",
  "- Critical boundary: Cdcp1 is sparse and non-MΦ enriched; ambient/global contribution is not excluded.",
  "- Classification: **strongest directionally concordant Human-Mouse candidate bridge, with unresolved Mouse source attribution**.",
  "",
  "### 3. TGFB1/Tgfb1-TGFβR",
  "- Human: robust auxiliary profibrotic pathway.",
  "- Mouse: Tgfb1 burden is not attenuated by Tx despite lower fibrogenic HSC burden.",
  "- Classification: **canonical pathway conservation, but no evidence that Tgfb1 attenuation explains the Mouse Tx effect**.",
  "",
  "### 4. PDGF-family -> HSC PDGFR",
  "- Human: PDGFA-PDGFRA is associated with the ECM-high HSC phenotype.",
  "- Mouse: Pdgfa-Pdgfra is approximately unchanged, whereas Pdgfb-Pdgfrb is attenuated with Tx.",
  "- Classification: **functional-axis conservation with context/species-specific paralog usage divergence**.",
  "",
  "## Working cross-species model",
  "",
  "Human and Mouse data support a conserved multi-axis myeloid-to-HSC fibrogenic network, but individual ligand usage and treatment responsiveness differ by species/context. CDCP1/Cdcp1 provides the strongest directionally concordant candidate bridge, whereas the PDGF family shows conserved receiver-level biology with a Human PDGFA-PDGFRA versus Mouse Pdgfb-Pdgfrb divergence. SPP1 and TGFB1 remain biologically conserved profibrotic axes but are not attenuated in the Mouse Tx response.",
  "",
  "## Inference boundaries",
  "",
  "- Human: exploratory n=8.",
  "- Mouse: descriptive n=2/group.",
  "- No cross-species correlation or treatment-direction match establishes a causal intracellular chain.",
  "- Mouse Cdcp1 provenance remains provisional until the full-library denominator correction is completed."
)

writeLines(
  md,
  file.path(
    DIR_DOC,
    "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0.md"
  )
)

# ------------------------------------------------------------
# Provenance / integrity
# ------------------------------------------------------------

provenance <- tibble(
  component = c(
    "Human canonical 4-axis hierarchy",
    "Mouse focused total-myeloid closure",
    "Mouse TREM2-focused closure",
    "Mouse canonical HSC receiver effects",
    "Mouse Cdcp1 raw-count pseudobulk",
    "Mouse PDGF family effects",
    "Mouse Cdcp1 whole-host provenance boundary"
  ),
  source = c(
    HUMAN_HIERARCHY,
    MOUSE_TOTAL_MATRIX,
    MOUSE_TREM2_MATRIX,
    MOUSE_RECEIVER,
    MOUSE_CDCP1_RANK,
    MOUSE_PDGF_EFFECT,
    ifelse(
      file.exists(MOUSE_PROV_TABLE),
      MOUSE_PROV_TABLE,
      paste0(
        MOUSE_PROV_TABLE,
        " [NOT FOUND]"
      )
    )
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "source_provenance_v1.0.tsv"
  )
)

integrity <- tibble(
  item = c(
    "Human_axes_n",
    "Human_axes",
    "Human_inference",
    "Mouse_samples_n",
    "Mouse_group_inference",
    "Mouse_receiver_primary",
    "Mouse_Cdcp1_rawCount_support",
    "Mouse_v6.4.10.4_quantitative_use",
    "Cross_species_unit",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    "4",
    "SPP1;CDCP1;TGFB1;PDGFA",
    "exploratory_sample_level_n8",
    "4",
    "descriptive_only_n2_per_group",
    "host-weighted ECM-activated HSC core-collagen burden",
    "TRUE",
    "FALSE - qualitative provenance boundary only; full-library denominator correction pending",
    "functional_axis_not_gene_identity",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.0.tsv"
  )
)

# Basic integrity assertions
if (!all(cross$human_priority == 1:4)) {
  stop("Cross-species hierarchy lost canonical Human priority order.")
}

if (!all(is.finite(cross$mouse_receiver_ratio))) {
  stop("Mouse receiver ratio missing.")
}

if (!all(cross$mouse_common_receiver_reduction)) {
  stop(
    "Expected canonical Mouse host-weighted ECM-HSC core-collagen reduction was not reproduced."
  )
}

saveRDS(
  list(
    human = human_core,
    mouse_receiver = receiver_focus,
    mouse_sender = mouse_sender,
    cdcp1_rawCount_support = cdcp1_pb_support,
    cdcp1_provenance_boundary = cdcp1_prov,
    mouse_pdgf = mouse_pdgf_focus,
    mouse_axis = mouse_axis,
    cross_species = cross,
    concise = concise,
    boundaries = boundaries,
    provenance = provenance,
    integrity = integrity
  ),
  file.path(
    DIR_OBJ,
    "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0_results.rds"
  )
)

cat("\n============================================================\n")
cat("PRIMARY CROSS-SPECIES 4-AXIS HIERARCHY\n")
cat("============================================================\n")
print(cross)

cat("\n============================================================\n")
cat("CONCISE CROSS-SPECIES HIERARCHY\n")
cat("============================================================\n")
print(concise)

cat("\n============================================================\n")
cat("MOUSE Cdcp1 RAW-COUNT SUPPORT\n")
cat("============================================================\n")
print(cdcp1_pb_support)

cat("\n============================================================\n")
cat("INTEGRITY\n")
cat("============================================================\n")
print(integrity)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
