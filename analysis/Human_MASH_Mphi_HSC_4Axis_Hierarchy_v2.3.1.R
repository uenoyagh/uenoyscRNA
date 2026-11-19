#!/usr/bin/env Rscript

# ============================================================
# Human MASH Mphi -> HSC four-axis mechanistic hierarchy v2.3.1
#
# Integrates:
#   1. SPP1 -> HSC integrins
#   2. CDCP1-associated -> EGFR-supported HSC axis
#   3. TGFB1 -> TGF-beta receptor axis
#   4. PDGFA -> PDGFRA axis
#
# PRINCIPLES
#   - Do not force all four axes into the same evidence type.
#   - SPP1/CDCP1/TGFB1 inherit the validated Focused Ligand Closure v2.2.
#   - PDGFA inherits the validated HSC phenotype analysis v1.3.2.
#   - Preserve explicit mechanistic boundaries.
#   - Human n=8; exploratory sample-level inference.
#   - No source-RDS modification or reclustering.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

FOCUS_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Focused_Ligand_Closure_v2.2"
)

PDGF_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.2"
)

SENDER_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_PDGFA_SenderState_Attribution_v1.2.0"
)

FOCUS_FINAL <- file.path(
  FOCUS_ROOT,
  "tables",
  "11_FINAL_focus_ligand_comparison.tsv"
)

FOCUS_RECEPTOR <- file.path(
  FOCUS_ROOT,
  "tables",
  "07_focus_ligand_HSC_receptor_support.tsv"
)

FOCUS_TARGET <- file.path(
  FOCUS_ROOT,
  "tables",
  "09_focus_ligand_target_support_summary.tsv"
)

PDGF_PRIMARY <- file.path(
  PDGF_ROOT,
  "tables",
  "03_PRIMARY_PDGFA_validated_HSC_phenotype.tsv"
)

PDGF_INTERP <- file.path(
  PDGF_ROOT,
  "tables",
  "04_PRIMARY_interpretation.tsv"
)

PDGF_SENDER <- file.path(
  SENDER_ROOT,
  "tables",
  "02_state_component_HSC_associations.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_DOC <- file.path(OUT, "docs")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")
DIR_OBJ <- file.path(OUT, "objects")

for (d in c(DIR_TABLE, DIR_DOC, DIR_LOG, DIR_PROV, DIR_OBJ)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG <- file.path(
  DIR_LOG,
  "Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH Mphi -> HSC four-axis hierarchy v2.3.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")
cat("Patch note: v2.3.1 removes the no-finite-value sender warning, preserves raw receptor-prior support separately, uses canonical receptor labels in the hierarchy, and adds numeric integrity checks for n and NicheNet ranks.\n\n")

stopifnot(
  file.exists(FOCUS_FINAL),
  file.exists(FOCUS_RECEPTOR),
  file.exists(FOCUS_TARGET),
  file.exists(PDGF_PRIMARY),
  file.exists(PDGF_INTERP),
  file.exists(PDGF_SENDER)
)

focus_final <- read_tsv(
  FOCUS_FINAL,
  show_col_types = FALSE
)

focus_receptor <- read_tsv(
  FOCUS_RECEPTOR,
  show_col_types = FALSE
)

focus_target <- read_tsv(
  FOCUS_TARGET,
  show_col_types = FALSE
)

pdgf_primary <- read_tsv(
  PDGF_PRIMARY,
  show_col_types = FALSE
)

pdgf_interp <- read_tsv(
  PDGF_INTERP,
  show_col_types = FALSE
)

pdgf_sender <- read_tsv(
  PDGF_SENDER,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Integrity checks
# ------------------------------------------------------------

needed_focus <- c(
  "ligand",
  "rho",
  "p_exact",
  "FDR_exact_within_endpoint",
  "partial_rho_NAS",
  "partial_p_NAS",
  "NicheNet_rank",
  "NicheNet_aupr_corrected",
  "NicheNet_targets_tested",
  "targets_positive_rho",
  "targets_exact_p_lt_0.05",
  "targets_exact_FDR_lt_0.05",
  "targets_partial_p_lt_0.05",
  "best_target",
  "best_target_exact_p"
)

if (!all(needed_focus %in% colnames(focus_final))) {
  stop(
    "Focused Ligand Closure final comparison lacks required columns: ",
    paste(setdiff(needed_focus, colnames(focus_final)), collapse = ", ")
  )
}

needed_pdgf <- c(
  "predictor",
  "endpoint",
  "n",
  "predictor_nonzero_n",
  "rho_raw",
  "p_nominal",
  "rho_partial_NAS",
  "loo_sign_consistency",
  "p_BH_all",
  "predictor_sparse",
  "endpoint_sparse",
  "robustness_class"
)

if (!all(needed_pdgf %in% colnames(pdgf_primary))) {
  stop(
    "PDGFA primary table lacks required columns: ",
    paste(setdiff(needed_pdgf, colnames(pdgf_primary)), collapse = ", ")
  )
}

# ------------------------------------------------------------
# Receptor support summarization for established 3 axes
# ------------------------------------------------------------

receptor_summary <- focus_receptor %>%
  filter(
    from %in% c("SPP1", "CDCP1", "TGFB1"),
    receptor_measurable %in% TRUE
  ) %>%
  group_by(from) %>%
  summarise(
    measurable_receptors =
      paste(
        unique(to),
        collapse = ";"
      ),
    primary_expressed_receptors =
      paste(
        unique(
          to[
            receptor_primary_expressed %in% TRUE
          ]
        ),
        collapse = ";"
      ),
    max_mean_detection_fraction =
      suppressWarnings(
        max(
          mean_detection_fraction,
          na.rm = TRUE
        )
      ),
    .groups = "drop"
  ) %>%
  mutate(
    max_mean_detection_fraction =
      ifelse(
        is.infinite(max_mean_detection_fraction),
        NA_real_,
        max_mean_detection_fraction
      )
  )

# ------------------------------------------------------------
# Established three axes: extract quantitative evidence
# ------------------------------------------------------------

established <- focus_final %>%
  filter(
    ligand %in% c("SPP1", "CDCP1", "TGFB1")
  ) %>%
  left_join(
    receptor_summary,
    by = c("ligand" = "from")
  ) %>%
  left_join(
    focus_target %>%
      select(
        ligand,
        NicheNet_targets_tested,
        targets_positive_rho,
        targets_exact_p_lt_0.05,
        targets_exact_FDR_lt_0.05,
        targets_partial_p_lt_0.05,
        best_target,
        best_target_exact_p
      ) %>%
      rename_with(
        ~ paste0("target_", .x),
        -ligand
      ),
    by = "ligand"
  )

# Prefer final-comparison target fields; target summary is retained as audit.
established_table <- established %>%
  transmute(
    axis = case_when(
      ligand == "SPP1" ~ "SPP1 -> HSC integrins",
      ligand == "CDCP1" ~ "CDCP1-associated -> HSC EGFR-supported axis",
      ligand == "TGFB1" ~ "TGFB1 -> TGF-beta receptor axis",
      TRUE ~ ligand
    ),
    ligand = ligand,
    primary_receiver_endpoint =
      "Hallmark_whole_HSC_pseudobulk_score",
    n = 8L,
    rho = rho,
    p_value_primary = p_exact,
    FDR_primary = FDR_exact_within_endpoint,
    partial_rho_NAS = partial_rho_NAS,
    partial_p_NAS = partial_p_NAS,
    LOO_sign_consistency = NA_real_,
    predictor_nonzero_n = NA_integer_,
    receiver_support_raw = primary_expressed_receptors,
    receiver_support = case_when(
      ligand == "SPP1" ~
        "ITGA4;ITGA5;ITGA9;ITGAV;ITGB1;ITGB5",
      ligand == "CDCP1" ~
        "EGFR (expression support; direct CDCP1-EGFR LR not established)",
      ligand == "TGFB1" ~
        "TGFBR1;TGFBR2;TGFBR3",
      TRUE ~
        primary_expressed_receptors
    ),
    receptor_max_mean_detection_fraction =
      max_mean_detection_fraction,
    NicheNet_rank = NicheNet_rank,
    NicheNet_AUPR_corrected = NicheNet_aupr_corrected,
    NicheNet_targets_tested = NicheNet_targets_tested,
    NicheNet_targets_positive_rho = targets_positive_rho,
    NicheNet_targets_exact_p_lt_0.05 =
      targets_exact_p_lt_0.05,
    NicheNet_targets_FDR_lt_0.05 =
      targets_exact_FDR_lt_0.05,
    NicheNet_targets_partial_p_lt_0.05 =
      targets_partial_p_lt_0.05,
    best_target = best_target,
    best_target_p_exact = best_target_exact_p
  )

# ------------------------------------------------------------
# PDGFA: extract the validated primary Human phenotype result
# ------------------------------------------------------------

pdgfa_row <- pdgf_primary %>%
  filter(
    predictor == "total_PDGFA_PDGFRA_product",
    endpoint == "ECM_high_fraction"
  )

if (nrow(pdgfa_row) != 1) {
  stop(
    "Expected exactly one PDGFA primary row for total_PDGFA_PDGFRA_product -> ECM_high_fraction."
  )
}

pdgfa_table <- pdgfa_row %>%
  transmute(
    axis = "PDGFA -> PDGFRA",
    ligand = "PDGFA",
    primary_receiver_endpoint = endpoint,
    n = n,
    rho = rho_raw,
    p_value_primary = p_nominal,
    FDR_primary = p_BH_all,
    partial_rho_NAS = rho_partial_NAS,
    partial_p_NAS = NA_real_,
    LOO_sign_consistency = loo_sign_consistency,
    predictor_nonzero_n = predictor_nonzero_n,
    receiver_support_raw = "PDGFRA",
    receiver_support = "PDGFRA",
    receptor_max_mean_detection_fraction = NA_real_,
    NicheNet_rank = NA_real_,
    NicheNet_AUPR_corrected = NA_real_,
    NicheNet_targets_tested = NA_integer_,
    NicheNet_targets_positive_rho = NA_integer_,
    NicheNet_targets_exact_p_lt_0.05 = NA_integer_,
    NicheNet_targets_FDR_lt_0.05 = NA_integer_,
    NicheNet_targets_partial_p_lt_0.05 = NA_integer_,
    best_target = NA_character_,
    best_target_p_exact = NA_real_
  )

quantitative <- bind_rows(
  established_table,
  pdgfa_table
)

# ------------------------------------------------------------
# PDGFA sender-state audit
# ------------------------------------------------------------

pdgfa_sender_summary <- pdgf_sender %>%
  filter(
    gene == "PDGFA",
    component == "ligand_burden",
    endpoint %in%
      c(
        "aHSC_ECM_remodeling",
        "wholeHSC_ECM_remodeling"
      )
  ) %>%
  group_by(sender_state) %>%
  summarise(
    nonzero_n_max =
      ifelse(
        any(is.finite(nonzero_n)),
        max(nonzero_n[is.finite(nonzero_n)]),
        NA_integer_
      ),
    max_abs_rho_raw =
      ifelse(
        any(is.finite(rho_raw)),
        max(abs(rho_raw[is.finite(rho_raw)])),
        NA_real_
      ),
    max_abs_partial_NAS =
      ifelse(
        any(is.finite(rho_partial_NAS)),
        max(abs(rho_partial_NAS[is.finite(rho_partial_NAS)])),
        NA_real_
      ),
    any_sparse =
      any(sparse_component %in% TRUE),
    .groups = "drop"
  )

write_tsv(
  pdgfa_sender_summary,
  file.path(
    DIR_TABLE,
    "01_PDGFA_sender_state_audit.tsv"
  )
)

# ------------------------------------------------------------
# Mechanistic hierarchy
# ------------------------------------------------------------

hierarchy_meta <- tribble(
  ~ligand, ~hierarchy_order, ~hierarchy_tier,
  ~mechanistic_role, ~sender_scope,
  ~statistical_position, ~safe_statement, ~do_not_claim,

  "SPP1", 1L, "Tier 1 - Primary mechanistic axis",
  "Fibrotic niche / adhesion / HSC activation",
  "LAM/TREM2-associated myeloid context; not exclusive to TREM2 state",
  "Strong receptor + NicheNet + target support; sample-level association persists but NAS adjustment attenuates",
  "Myeloid SPP1 -> HSC integrin signaling is the primary Human MASH Mphi-HSC profibrotic axis.",
  "Do not call this a TREM2-specific direct axis; do not imply one integrin is uniquely causal.",

  "CDCP1", 2L, "Tier 1B - Strongest novel candidate",
  "Novel HSC-associated signaling candidate",
  "Strong TREM2-state specificity in prior focused closure",
  "Strongest sample-level association and NAS-adjusted robustness among focused ligands; NicheNet rank lower than SPP1",
  "CDCP1-associated signaling is the strongest novel TREM2-linked candidate and is supported by HSC EGFR expression.",
  "Do not state that CDCP1 is a proven direct EGFR ligand or that CDCP1->EGFR causality is established.",

  "TGFB1", 3L, "Tier 2 - Robust auxiliary profibrotic axis",
  "Canonical profibrotic / TGF-beta response support",
  "TREM2-state-associated signal; absolute TGFB1 source is not necessarily TREM2-dominant",
  "Robust sample-level and NAS-adjusted association with broad target support; direct-driver specificity is weaker",
  "TGFB1/TGF-beta signaling is a robust TREM2-state-associated auxiliary profibrotic signal.",
  "Do not call TGFB1 the dominant direct ligand or imply TREM2 is the dominant absolute TGFB1 source.",

  "PDGFA", 4L, "Tier 3 - Robust exploratory ECM-high axis",
  "PDGFRA-positive HSC / ECM-high fibrogenic phenotype association",
  "Cross-state myeloid; Resident/Kupffer-like and TREM2-associated cells are major candidate sources, but TREM2 PDGFA is sparse",
  "Robust ECM-high association after NAS adjustment and LOO; equivalent NicheNet target-closure has not yet been performed",
  "Cross-state myeloid PDGFA-associated signaling is linked to a PDGFRA-positive ECM-high HSC phenotype in Human MASH.",
  "Do not claim PDGFA directly drives HSC proliferation; do not label the axis TREM2-specific; do not equate association with biochemical PDGFRA activation."
)

hierarchy <- quantitative %>%
  left_join(
    hierarchy_meta,
    by = "ligand"
  ) %>%
  arrange(
    hierarchy_order
  )

# ------------------------------------------------------------
# Evidence-grade adjudication
# ------------------------------------------------------------

hierarchy <- hierarchy %>%
  mutate(
    sample_level_support = case_when(
      ligand == "CDCP1" &
        is.finite(partial_p_NAS) &
        partial_p_NAS < 0.01 ~
        "Very strong",

      ligand == "TGFB1" &
        is.finite(partial_p_NAS) &
        partial_p_NAS < 0.05 ~
        "Strong",

      ligand == "SPP1" &
        is.finite(p_value_primary) &
        p_value_primary < 0.05 ~
        "Moderate-strong",

      ligand == "PDGFA" &
        rho >= 0.7 &
        partial_rho_NAS >= 0.5 &
        LOO_sign_consistency == 1 ~
        "Moderate-strong exploratory",

      TRUE ~
        "Exploratory"
    ),

    receptor_support_grade = case_when(
      ligand == "SPP1" ~ "Strong - multiple measurable integrins",
      ligand == "CDCP1" ~ "Present - EGFR measurable; ligand-receptor directness unresolved",
      ligand == "TGFB1" ~ "Strong - TGFBR1/2/3 and integrin-related support",
      ligand == "PDGFA" ~ "Present - PDGFRA measurable; ECM-high enrichment is modest",
      TRUE ~ NA_character_
    ),

    target_support_grade = case_when(
      ligand == "SPP1" ~ "Strong",
      ligand == "CDCP1" ~ "Moderate",
      ligand == "TGFB1" ~ "Strong/broad",
      ligand == "PDGFA" ~ "Not yet tested in equivalent focused NicheNet closure",
      TRUE ~ NA_character_
    ),

    overall_evidence_grade = case_when(
      ligand == "SPP1" ~ "HIGH - primary mechanistic axis",
      ligand == "CDCP1" ~ "HIGH for association / MODERATE for receptor mechanism",
      ligand == "TGFB1" ~ "HIGH as auxiliary pathway",
      ligand == "PDGFA" ~ "MODERATE - robust exploratory axis",
      TRUE ~ NA_character_
    )
  )

# Reorder for readable output.
hierarchy_out <- hierarchy %>%
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
    LOO_sign_consistency,
    predictor_nonzero_n,
    receiver_support,
    receiver_support_raw,
    receptor_support_grade,
    receptor_max_mean_detection_fraction,
    NicheNet_rank,
    NicheNet_AUPR_corrected,
    NicheNet_targets_tested,
    NicheNet_targets_positive_rho,
    NicheNet_targets_exact_p_lt_0.05,
    NicheNet_targets_FDR_lt_0.05,
    NicheNet_targets_partial_p_lt_0.05,
    best_target,
    best_target_p_exact,
    sample_level_support,
    target_support_grade,
    overall_evidence_grade,
    statistical_position,
    safe_statement,
    do_not_claim
  )

write_tsv(
  hierarchy_out,
  file.path(
    DIR_TABLE,
    "02_PRIMARY_Human_Mphi_HSC_4Axis_mechanistic_hierarchy.tsv"
  )
)

# ------------------------------------------------------------
# Concise manuscript-ready table
# ------------------------------------------------------------

concise <- hierarchy_out %>%
  transmute(
    Priority = hierarchy_order,
    Axis = axis,
    Role = mechanistic_role,
    Sender_scope = sender_scope,
    Primary_HSC_endpoint = primary_receiver_endpoint,
    Sample_level_rho = round(rho, 3),
    NAS_adjusted_rho = round(partial_rho_NAS, 3),
    NicheNet_rank = NicheNet_rank,
    Receptor_support = receptor_support_grade,
    Target_support = target_support_grade,
    Evidence_grade = overall_evidence_grade,
    Interpretation = safe_statement
  )

write_tsv(
  concise,
  file.path(
    DIR_TABLE,
    "03_CONCISE_Human_Mphi_HSC_4Axis_hierarchy.tsv"
  )
)

# ------------------------------------------------------------
# Markdown checkpoint
# ------------------------------------------------------------

md <- c(
  "# Human MASH Mphi-HSC four-axis mechanistic hierarchy v2.3.1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Updated hierarchy",
  "",
  "1. **SPP1 -> HSC integrins** — primary mechanistic axis.",
  "2. **CDCP1-associated signaling** — strongest novel TREM2-linked candidate; receptor mechanism remains unresolved.",
  "3. **TGFB1 -> TGF-beta receptor signaling** — robust auxiliary profibrotic pathway.",
  "4. **PDGFA -> PDGFRA** — newly integrated cross-state myeloid axis associated with ECM-high HSC phenotype; robust exploratory evidence.",
  "",
  "## Critical interpretation boundaries",
  "",
  "- SPP1: use `myeloid SPP1 -> HSC integrin axis`; do not call it TREM2-specific direct signaling.",
  "- CDCP1: do not claim proven direct CDCP1 -> EGFR ligand-receptor causality.",
  "- TGFB1: do not call it the dominant direct ligand; TREM2 association is not identical to dominant absolute source.",
  "- PDGFA: do not claim direct HSC proliferation from the current Human scRNA-seq data; TREM2-specificity is not established.",
  "",
  "## Overall working model",
  "",
  "Human MASH HSC fibrogenic activation is best represented as a multi-axis myeloid signaling network rather than a single ligand pathway. SPP1-integrin signaling remains the primary fibrotic-niche axis; CDCP1-associated signaling provides the strongest novel TREM2-linked candidate; TGFB1 provides a canonical auxiliary profibrotic program; and PDGFA-PDGFRA adds a cross-state myeloid axis linked specifically to the ECM-high HSC phenotype.",
  "",
  "## Inference boundary",
  "",
  "All Human sample-level associations are exploratory (n=8). Correlation, receptor transcript detection, and ligand-target prediction do not establish biochemical receptor activation or causal ligand dependence."
)

writeLines(
  md,
  file.path(
    DIR_DOC,
    "Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1.md"
  )
)

# ------------------------------------------------------------
# Numeric integrity audit
# ------------------------------------------------------------

numeric_audit <- hierarchy_out %>%
  transmute(
    ligand,
    n,
    rho,
    partial_rho_NAS,
    NicheNet_rank,
    expected_n = 8,
    n_ok = n == 8,
    expected_rank = case_when(
      ligand == "SPP1" ~ 4,
      ligand == "CDCP1" ~ 20,
      ligand == "TGFB1" ~ 206,
      ligand == "PDGFA" ~ NA_real_,
      TRUE ~ NA_real_
    ),
    rank_ok = case_when(
      is.na(expected_rank) ~ TRUE,
      TRUE ~ NicheNet_rank == expected_rank
    )
  )

if (!all(numeric_audit$n_ok)) {
  stop("Numeric integrity audit failed: one or more axes do not have n=8.")
}

if (!all(numeric_audit$rank_ok)) {
  stop("Numeric integrity audit failed: NicheNet rank mismatch.")
}

write_tsv(
  numeric_audit,
  file.path(
    DIR_PROV,
    "numeric_integrity_audit_v2.3.1.tsv"
  )
)

# ------------------------------------------------------------
# Provenance
# ------------------------------------------------------------

provenance <- tibble(
  component = c(
    "SPP1_CDCP1_TGFB1_final_comparison",
    "SPP1_CDCP1_TGFB1_receptor_support",
    "SPP1_CDCP1_TGFB1_target_support",
    "PDGFA_validated_HSC_primary",
    "PDGFA_validated_HSC_interpretation",
    "PDGFA_sender_state_attribution"
  ),
  source = c(
    FOCUS_FINAL,
    FOCUS_RECEPTOR,
    FOCUS_TARGET,
    PDGF_PRIMARY,
    PDGF_INTERP,
    PDGF_SENDER
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "source_provenance_v2.3.1.tsv"
  )
)

integrity <- tibble(
  item = c(
    "axes_n",
    "axes",
    "Human_sample_inference",
    "PDGFA_equivalent_NicheNet_closure",
    "canonical_receiver_display",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    "4",
    "SPP1;CDCP1;TGFB1;PDGFA",
    "exploratory_n8",
    "NOT_YET_PERFORMED",
    "SPP1=integrins;CDCP1=EGFR-expression-support-only;TGFB1=TGFBR1/2/3;PDGFA=PDGFRA",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v2.3.1.tsv"
  )
)

saveRDS(
  list(
    hierarchy = hierarchy_out,
    concise = concise,
    pdgfa_sender_summary = pdgfa_sender_summary,
    provenance = provenance,
    integrity = integrity
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1_results.rds"
  )
)

cat("\n============================================================\n")
cat("PRIMARY 4-AXIS MECHANISTIC HIERARCHY\n")
cat("============================================================\n")
print(hierarchy_out)

cat("\n============================================================\n")
cat("CONCISE HIERARCHY\n")
cat("============================================================\n")
print(concise)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
