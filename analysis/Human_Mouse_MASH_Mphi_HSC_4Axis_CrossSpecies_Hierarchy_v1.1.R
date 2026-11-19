#!/usr/bin/env Rscript

# ============================================================
# Human-Mouse MASH Mphi -> HSC 4-axis cross-species hierarchy v1.1
#
# PURPOSE
#   Update the validated v1.0 cross-species hierarchy with the
#   corrected Mouse Cdcp1 whole-cell source decomposition v6.4.11.1.
#
# KEY UPDATE
#   Mouse Cdcp1 is directionally concordant with Tx response,
#   but its absolute whole-liver source is predominantly non-Mphi.
#   Therefore CDCP1/Cdcp1 is retained as a directionally concordant
#   gene-level bridge but is no longer described as a conserved
#   Mphi-derived mechanistic bridge.
#
#   PDGF-family -> HSC PDGFR is elevated to the strongest
#   cross-species Mphi->HSC FUNCTIONAL therapeutic-response bridge,
#   while preserving the Human PDGFA-PDGFRA vs Mouse Pdgfb-Pdgfrb
#   paralog/context divergence.
#
# INFERENCE BOUNDARIES
#   Human: n=8 exploratory sample-level.
#   Mouse: n=2/group descriptive.
#   No causal inference.
#
# No RDS is modified. No reclustering.
# ============================================================

options(stringsAsFactors = FALSE, width = 240)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

V10_ROOT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0"
)

V10_PRIMARY <- file.path(
  V10_ROOT,
  "tables",
  "08_PRIMARY_Human_Mouse_4Axis_cross_species_hierarchy.tsv"
)

CDCP1_ROOT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Cdcp1_WholeCell_SourceDecomposition_v6.4.11.1"
)

CDCP1_SOURCE <- file.path(
  CDCP1_ROOT,
  "tables",
  "05_PRIMARY_Cdcp1_Mphi_vs_nonMphi_source_adjudication.tsv"
)

CDCP1_TX <- file.path(
  CDCP1_ROOT,
  "tables",
  "03_PRIMARY_Mphi_vs_nonMphi_Cdcp1_Sham_vs_Tx_fullLibrary.tsv"
)

CDCP1_CELLTYPE <- file.path(
  CDCP1_ROOT,
  "tables",
  "07_PRIMARY_Cdcp1_whole_celltype_source_ranking.tsv"
)

CDCP1_MPHI_STATE <- file.path(
  CDCP1_ROOT,
  "tables",
  "09_PRIMARY_Cdcp1_Mphi_state_source_ranking.tsv"
)

CDCP1_INTEGRITY <- file.path(
  CDCP1_ROOT,
  "provenance",
  "integrity_summary_v6.4.11.1.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.1"
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
  "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.1.log"
)

zz <- file(LOG, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("Human-Mouse MASH Mphi-HSC 4-axis cross-species hierarchy v1.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

required <- c(
  V10_PRIMARY,
  CDCP1_SOURCE,
  CDCP1_TX,
  CDCP1_CELLTYPE,
  CDCP1_MPHI_STATE,
  CDCP1_INTEGRITY
)

missing <- required[!file.exists(required)]

if (length(missing) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing, collapse = "\n")
  )
}

# ------------------------------------------------------------
# Read validated inputs
# ------------------------------------------------------------

v10 <- read_tsv(
  V10_PRIMARY,
  show_col_types = FALSE
)

cdcp1_source <- read_tsv(
  CDCP1_SOURCE,
  show_col_types = FALSE
)

cdcp1_tx <- read_tsv(
  CDCP1_TX,
  show_col_types = FALSE
)

cdcp1_celltype <- read_tsv(
  CDCP1_CELLTYPE,
  show_col_types = FALSE
)

cdcp1_mphi_state <- read_tsv(
  CDCP1_MPHI_STATE,
  show_col_types = FALSE
)

cdcp1_integrity <- read_tsv(
  CDCP1_INTEGRITY,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Input integrity
# ------------------------------------------------------------

stopifnot(
  nrow(v10) == 4,
  all(v10$human_ligand == c("SPP1", "CDCP1", "TGFB1", "PDGFA")),
  nrow(cdcp1_source) == 1
)

get_integrity <- function(key) {
  x <- cdcp1_integrity %>%
    filter(item == key) %>%
    pull(value)

  if (length(x) != 1) {
    stop("Integrity item not uniquely resolved: ", key)
  }

  x[[1]]
}

stopifnot(
  get_integrity("feature_subsetting_before_library_size") == "FALSE",
  get_integrity("v6.4.10.4_11913_gene_denominator_reused") == "FALSE",
  get_integrity("source_RDS_overwritten") == "FALSE",
  get_integrity("reclustering_performed") == "FALSE"
)

# ------------------------------------------------------------
# Extract corrected Cdcp1 source metrics
# ------------------------------------------------------------

mphi_tx <- cdcp1_tx %>%
  filter(compartment == "Mphi")

nonmphi_tx <- cdcp1_tx %>%
  filter(compartment == "nonMphi")

if (nrow(mphi_tx) != 1 || nrow(nonmphi_tx) != 1) {
  stop("Mphi/nonMphi Tx rows not uniquely resolved.")
}

trem2_state <- cdcp1_mphi_state %>%
  filter(sender_state == "Lipid-associated/TREM2-Mphi")

if (nrow(trem2_state) != 1) {
  stop("TREM2 Mphi Cdcp1 source row not uniquely resolved.")
}

top_sources <- cdcp1_celltype %>%
  arrange(desc(Cdcp1_UMI_share))

chol <- top_sources %>%
  filter(whole_celltype == "Cholangiocyte")

hep <- top_sources %>%
  filter(whole_celltype == "Hepatocyte")

kup <- top_sources %>%
  filter(whole_celltype == "Kupffer_Macrophage")

mono <- top_sources %>%
  filter(whole_celltype == "Monocyte")

if (
  nrow(chol) != 1 ||
  nrow(hep) != 1 ||
  nrow(kup) != 1 ||
  nrow(mono) != 1
) {
  stop("Expected whole-celltype source rows were not uniquely resolved.")
}

cdcp1_update <- tibble(
  Mouse_Mphi_CPM_full_library =
    cdcp1_source$Mphi_CPM_full_library,
  Mouse_nonMphi_CPM_full_library =
    cdcp1_source$nonMphi_CPM_full_library,
  Mouse_Mphi_detection_fraction =
    cdcp1_source$Mphi_detection_fraction,
  Mouse_nonMphi_detection_fraction =
    cdcp1_source$nonMphi_detection_fraction,
  Mouse_Mphi_Cdcp1_UMI_share =
    cdcp1_source$Mphi_Cdcp1_UMI_share,
  Mouse_nonMphi_Cdcp1_UMI_share =
    cdcp1_source$nonMphi_Cdcp1_UMI_share,
  Mouse_log2_Mphi_vs_nonMphi_CPM =
    cdcp1_source$log2_Mphi_vs_nonMphi_CPM,
  Mouse_Mphi_CPM_Tx_over_Sham =
    mphi_tx$CPM_Tx_over_Sham,
  Mouse_nonMphi_CPM_Tx_over_Sham =
    nonmphi_tx$CPM_Tx_over_Sham,
  Mouse_TREM2_share_within_Mphi =
    trem2_state$Cdcp1_UMI_share_within_Mphi,
  Mouse_Cholangiocyte_Cdcp1_UMI_share =
    chol$Cdcp1_UMI_share,
  Mouse_Hepatocyte_Cdcp1_UMI_share =
    hep$Cdcp1_UMI_share,
  Mouse_Kupffer_Cdcp1_UMI_share =
    kup$Cdcp1_UMI_share,
  Mouse_Monocyte_Cdcp1_UMI_share =
    mono$Cdcp1_UMI_share,
  Mouse_Cholangiocyte_plus_Hepatocyte_share =
    chol$Cdcp1_UMI_share +
    hep$Cdcp1_UMI_share
)

write_tsv(
  cdcp1_update,
  file.path(
    DIR_TABLE,
    "01_Mouse_Cdcp1_corrected_source_summary.tsv"
  )
)

# ------------------------------------------------------------
# Revised cross-species adjudication
#
# Keep Human priority separate from cross-species therapeutic bridge
# priority.
# ------------------------------------------------------------

revised <- v10 %>%
  mutate(
    cross_species_bridge_priority = case_when(
      human_ligand == "PDGFA" ~ 1L,
      human_ligand == "CDCP1" ~ 2L,
      human_ligand == "SPP1" ~ 3L,
      human_ligand == "TGFB1" ~ 4L,
      TRUE ~ NA_integer_
    ),

    cross_species_bridge_type = case_when(
      human_ligand == "PDGFA" ~
        "functional_Mphi_to_HSC_bridge_with_paralog_divergence",

      human_ligand == "CDCP1" ~
        "directionally_concordant_gene_level_bridge_with_source_divergence",

      human_ligand == "SPP1" ~
        "conserved_profibrotic_axis_not_Tx_attenuated",

      human_ligand == "TGFB1" ~
        "conserved_canonical_axis_not_Tx_attenuated",

      TRUE ~ NA_character_
    ),

    bridge_grade_v1.1 = case_when(
      human_ligand == "PDGFA" ~
        "MODERATE-HIGH",

      human_ligand == "CDCP1" ~
        "MODERATE",

      human_ligand == "SPP1" ~
        "MODERATE",

      human_ligand == "TGFB1" ~
        "MODERATE",

      TRUE ~ NA_character_
    ),

    source_conservation = case_when(
      human_ligand == "PDGFA" ~
        "partially conserved at myeloid-family signaling level; ligand/receptor paralog usage differs",

      human_ligand == "CDCP1" ~
        "NOT conserved as a dominant myeloid source in Mouse; Mouse whole-liver Cdcp1 is predominantly non-Mphi",

      human_ligand == "SPP1" ~
        "compatible with myeloid-source biology, but treatment attenuation absent in Mouse",

      human_ligand == "TGFB1" ~
        "myeloid/TREM2-associated signal exists, but absolute source dominance is not established",

      TRUE ~ NA_character_
    ),

    cross_species_safe_statement_v1.1 = case_when(
      human_ligand == "PDGFA" ~
        paste0(
          "The strongest cross-species Mphi-to-HSC functional therapeutic-response bridge is the PDGF-family axis: ",
          "Human MASH supports a myeloid PDGFA-PDGFRA/ECM-high HSC association, whereas Mouse Tx preferentially attenuates ",
          "Pdgfb-Pdgfrb. This supports functional-axis conservation with context/species-specific paralog usage."
        ),

      human_ligand == "CDCP1" ~
        paste0(
          "CDCP1/Cdcp1 remains strongly directionally concordant at the gene level: Human CDCP1 is the strongest novel ",
          "TREM2-linked HSC-associated candidate, and the sparse Mouse Mphi Cdcp1 component decreases with Tx. ",
          "However, Mouse whole-liver Cdcp1 is predominantly non-Mphi-derived, so conservation of a dominant Mphi-derived ",
          "CDCP1/Cdcp1 axis is not supported."
        ),

      human_ligand == "SPP1" ~
        "SPP1/Spp1-integrin biology is a conserved profibrotic candidate axis, but Mouse Tx does not attenuate Spp1 burden.",

      human_ligand == "TGFB1" ~
        "TGFB1/Tgfb1-TGF-beta signaling is a conserved canonical auxiliary profibrotic axis, but Mouse Tx does not attenuate Tgfb1 burden.",

      TRUE ~ NA_character_
    ),

    cross_species_do_not_claim_v1.1 = case_when(
      human_ligand == "PDGFA" ~
        "Do not claim one-to-one PDGFA-PDGFRA conservation or causal PDGF dependence.",

      human_ligand == "CDCP1" ~
        paste0(
          "Do not call Mouse Cdcp1 Mphi-derived, Mphi-specific, or a conserved Mphi-source mechanism; ",
          "do not claim direct CDCP1-EGFR causality."
        ),

      human_ligand == "SPP1" ~
        "Do not claim Spp1 suppression mediates the Mouse Tx anti-fibrotic effect.",

      human_ligand == "TGFB1" ~
        "Do not claim Tgfb1 suppression mediates the Mouse Tx effect.",

      TRUE ~ NA_character_
    )
  )

# Add corrected Cdcp1 metrics only to the CDCP1 row.
revised <- revised %>%
  mutate(
    Mouse_Mphi_Cdcp1_UMI_share =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_Mphi_Cdcp1_UMI_share,
        NA_real_
      ),

    Mouse_nonMphi_Cdcp1_UMI_share =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_nonMphi_Cdcp1_UMI_share,
        NA_real_
      ),

    Mouse_Cholangiocyte_plus_Hepatocyte_Cdcp1_share =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_Cholangiocyte_plus_Hepatocyte_share,
        NA_real_
      ),

    Mouse_TREM2_Cdcp1_share_within_Mphi =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_TREM2_share_within_Mphi,
        NA_real_
      ),

    Mouse_Mphi_Cdcp1_CPM_Tx_over_Sham =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_Mphi_CPM_Tx_over_Sham,
        NA_real_
      ),

    Mouse_nonMphi_Cdcp1_CPM_Tx_over_Sham =
      ifelse(
        human_ligand == "CDCP1",
        cdcp1_update$Mouse_nonMphi_CPM_Tx_over_Sham,
        NA_real_
      )
  )

write_tsv(
  revised,
  file.path(
    DIR_TABLE,
    "02_PRIMARY_Human_Mouse_4Axis_cross_species_hierarchy_v1.1.tsv"
  )
)

# ------------------------------------------------------------
# Concise bridge hierarchy
# ------------------------------------------------------------

concise_bridge <- revised %>%
  arrange(cross_species_bridge_priority) %>%
  transmute(
    Cross_species_bridge_priority =
      cross_species_bridge_priority,

    Functional_axis = case_when(
      human_ligand == "PDGFA" ~
        "PDGF-family -> HSC PDGFR",

      human_ligand == "CDCP1" ~
        "CDCP1/Cdcp1-associated",

      human_ligand == "SPP1" ~
        "SPP1/Spp1-integrin",

      human_ligand == "TGFB1" ~
        "TGFB1/Tgfb1-TGF-betaR",

      TRUE ~ human_axis
    ),

    Human_priority =
      human_priority,

    Human_evidence =
      human_evidence_grade,

    Human_rho =
      round(human_rho, 3),

    Human_NAS_adjusted_rho =
      round(human_partial_rho_NAS, 3),

    Mouse_Tx_axis_ratio = case_when(
      human_ligand == "PDGFA" ~
        mouse_secondary_LR_ratio,

      human_ligand == "CDCP1" ~
        Mouse_Mphi_Cdcp1_CPM_Tx_over_Sham,

      TRUE ~
        mouse_total_sender_ratio
    ),

    Mouse_ECM_HSC_coreCollagen_ratio =
      mouse_receiver_ratio,

    Source_conservation =
      source_conservation,

    Bridge_grade =
      bridge_grade_v1.1,

    Interpretation =
      cross_species_safe_statement_v1.1
  )

write_tsv(
  concise_bridge,
  file.path(
    DIR_TABLE,
    "03_CONCISE_cross_species_therapeutic_bridge_hierarchy_v1.1.tsv"
  )
)

# ------------------------------------------------------------
# Source divergence table
# ------------------------------------------------------------

source_divergence <- tibble(
  axis = "CDCP1/Cdcp1-associated",
  Mouse_Mphi_Cdcp1_UMI_share =
    cdcp1_update$Mouse_Mphi_Cdcp1_UMI_share,
  Mouse_nonMphi_Cdcp1_UMI_share =
    cdcp1_update$Mouse_nonMphi_Cdcp1_UMI_share,
  Mouse_Cholangiocyte_Cdcp1_UMI_share =
    cdcp1_update$Mouse_Cholangiocyte_Cdcp1_UMI_share,
  Mouse_Hepatocyte_Cdcp1_UMI_share =
    cdcp1_update$Mouse_Hepatocyte_Cdcp1_UMI_share,
  Mouse_Cholangiocyte_plus_Hepatocyte_share =
    cdcp1_update$Mouse_Cholangiocyte_plus_Hepatocyte_share,
  Mouse_TREM2_share_within_Mphi =
    cdcp1_update$Mouse_TREM2_share_within_Mphi,
  Mouse_Mphi_CPM_Tx_over_Sham =
    cdcp1_update$Mouse_Mphi_CPM_Tx_over_Sham,
  Mouse_nonMphi_CPM_Tx_over_Sham =
    cdcp1_update$Mouse_nonMphi_CPM_Tx_over_Sham,
  source_adjudication =
    "Mouse Cdcp1 is non-Mphi-dominant at whole-liver level; sparse Mphi/TREM2-associated Cdcp1 is present and decreases with Tx.",
  mechanistic_consequence =
    "Retain directional cross-species Cdcp1 concordance, but do not interpret Mouse Cdcp1 as a conserved dominant Mphi-derived ligand source."
)

write_tsv(
  source_divergence,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_CDCP1_Cdcp1_source_divergence_adjudication.tsv"
  )
)

# ------------------------------------------------------------
# Final model
# ------------------------------------------------------------

model <- tibble(
  model_component = c(
    "Strongest cross-species Mphi-to-HSC functional therapeutic bridge",
    "Strongest directionally concordant gene-level bridge",
    "Primary Human fibrotic-niche axis",
    "Canonical auxiliary profibrotic axis"
  ),
  axis = c(
    "PDGF-family -> HSC PDGFR",
    "CDCP1/Cdcp1-associated",
    "SPP1/Spp1-integrin",
    "TGFB1/Tgfb1-TGF-betaR"
  ),
  interpretation = c(
    "Human PDGFA-PDGFRA/ECM-high association and Mouse Tx attenuation of Pdgfb-Pdgfrb support conserved PDGF-to-HSC function with paralog/context divergence.",
    "Human CDCP1 is strongly TREM2-linked and Mouse Mphi Cdcp1 decreases with Tx, but Mouse whole-liver Cdcp1 is overwhelmingly non-Mphi-derived; source conservation is not supported.",
    "Strong Human primary axis and conserved fibrotic-niche biology, but Spp1 is not attenuated by Mouse Tx.",
    "Canonical conserved pathway, but Tgfb1 is not attenuated by Mouse Tx."
  )
)

write_tsv(
  model,
  file.path(
    DIR_TABLE,
    "05_FINAL_cross_species_working_model_v1.1.tsv"
  )
)

# ------------------------------------------------------------
# Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "input_v1.0_rows",
    "output_axes_n",
    "Human_priority_preserved",
    "cross_species_bridge_priority",
    "Cdcp1_full_library_denominator_verified",
    "Cdcp1_Mphi_UMI_share",
    "Cdcp1_nonMphi_UMI_share",
    "Cdcp1_cholangiocyte_plus_hepatocyte_share",
    "Cdcp1_TREM2_share_within_Mphi",
    "Cdcp1_Mphi_Tx_over_Sham",
    "Cdcp1_nonMphi_Tx_over_Sham",
    "Human_inference",
    "Mouse_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(nrow(v10)),
    as.character(nrow(revised)),
    as.character(all(revised$human_priority == 1:4)),
    paste(
      concise_bridge$Cross_species_bridge_priority,
      concise_bridge$Functional_axis,
      sep = ":",
      collapse = ";"
    ),
    "TRUE",
    as.character(cdcp1_update$Mouse_Mphi_Cdcp1_UMI_share),
    as.character(cdcp1_update$Mouse_nonMphi_Cdcp1_UMI_share),
    as.character(cdcp1_update$Mouse_Cholangiocyte_plus_Hepatocyte_share),
    as.character(cdcp1_update$Mouse_TREM2_share_within_Mphi),
    as.character(cdcp1_update$Mouse_Mphi_CPM_Tx_over_Sham),
    as.character(cdcp1_update$Mouse_nonMphi_CPM_Tx_over_Sham),
    "exploratory_sample_level_n8",
    "descriptive_only_n2_per_group",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.1.tsv"
  )
)

provenance <- tibble(
  component = c(
    "validated cross-species hierarchy v1.0",
    "corrected Mouse Cdcp1 source adjudication",
    "corrected Mouse Cdcp1 treatment effect",
    "Mouse Cdcp1 celltype source ranking",
    "Mouse Cdcp1 Mphi-state source ranking",
    "Mouse Cdcp1 full-library integrity"
  ),
  source = c(
    V10_PRIMARY,
    CDCP1_SOURCE,
    CDCP1_TX,
    CDCP1_CELLTYPE,
    CDCP1_MPHI_STATE,
    CDCP1_INTEGRITY
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "source_provenance_v1.1.tsv"
  )
)

# Assertions
stopifnot(
  nrow(revised) == 4,
  all(revised$human_priority == 1:4),
  setequal(
    revised$cross_species_bridge_priority,
    1:4
  ),
  cdcp1_update$Mouse_Mphi_Cdcp1_UMI_share < 0.05,
  cdcp1_update$Mouse_nonMphi_Cdcp1_UMI_share > 0.95,
  cdcp1_update$Mouse_Cholangiocyte_plus_Hepatocyte_share > 0.80,
  cdcp1_update$Mouse_Mphi_CPM_Tx_over_Sham < 1,
  cdcp1_update$Mouse_nonMphi_CPM_Tx_over_Sham < 1
)

saveRDS(
  list(
    v1.0 = v10,
    Cdcp1_source_update = cdcp1_update,
    revised_hierarchy = revised,
    concise_bridge = concise_bridge,
    source_divergence = source_divergence,
    final_model = model,
    integrity = integrity,
    provenance = provenance
  ),
  file.path(
    DIR_OBJ,
    "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.1_results.rds"
  )
)

# Markdown
md <- c(
  "# Human–Mouse MASH MΦ→HSC cross-species hierarchy v1.1",
  "",
  "## Key revision",
  "",
  "Corrected full-library source decomposition shows that Mouse whole-liver Cdcp1 is overwhelmingly non-macrophage-derived. Therefore CDCP1/Cdcp1 remains a directionally concordant gene-level bridge, but not a conserved dominant MΦ-derived source mechanism.",
  "",
  "## Cross-species therapeutic-response bridge hierarchy",
  "",
  "1. **PDGF-family → HSC PDGFR** — strongest cross-species MΦ→HSC functional therapeutic-response bridge; Human PDGFA-PDGFRA versus Mouse Pdgfb-Pdgfrb paralog/context divergence.",
  "2. **CDCP1/Cdcp1-associated** — strongest directionally concordant gene-level bridge, but major Mouse source divergence; sparse TREM2-MΦ Cdcp1 decreases with Tx.",
  "3. **SPP1/Spp1-integrin** — conserved profibrotic/fibrotic-niche biology, but not attenuated by Mouse Tx.",
  "4. **TGFB1/Tgfb1-TGFβR** — conserved canonical auxiliary pathway, but not attenuated by Mouse Tx.",
  "",
  "## Mouse Cdcp1 source boundary",
  "",
  paste0(
    "- MΦ Cdcp1 raw-UMI share: ",
    round(100 * cdcp1_update$Mouse_Mphi_Cdcp1_UMI_share, 2),
    "%."
  ),
  paste0(
    "- non-MΦ Cdcp1 raw-UMI share: ",
    round(100 * cdcp1_update$Mouse_nonMphi_Cdcp1_UMI_share, 2),
    "%."
  ),
  paste0(
    "- Cholangiocyte + Hepatocyte share: ",
    round(100 * cdcp1_update$Mouse_Cholangiocyte_plus_Hepatocyte_share, 2),
    "%."
  ),
  paste0(
    "- TREM2-associated MΦ share within the MΦ Cdcp1 pool: ",
    round(100 * cdcp1_update$Mouse_TREM2_share_within_Mphi, 2),
    "%."
  ),
  paste0(
    "- MΦ Cdcp1 CPM Tx/Sham: ",
    round(cdcp1_update$Mouse_Mphi_CPM_Tx_over_Sham, 3),
    "."
  ),
  paste0(
    "- non-MΦ Cdcp1 CPM Tx/Sham: ",
    round(cdcp1_update$Mouse_nonMphi_CPM_Tx_over_Sham, 3),
    "."
  ),
  "",
  "## Inference boundary",
  "",
  "- Human n=8 associations are exploratory.",
  "- Mouse n=2/group comparisons are descriptive.",
  "- Directional concordance and receptor expression do not prove causal ligand dependence.",
  "- Do not describe Mouse Cdcp1 as MΦ-derived or MΦ-specific."
)

writeLines(
  md,
  file.path(
    DIR_DOC,
    "Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.1.md"
  )
)

cat("\n============================================================\n")
cat("PRIMARY v1.1 CROSS-SPECIES HIERARCHY\n")
cat("============================================================\n")
print(
  revised %>%
    select(
      human_priority,
      cross_species_bridge_priority,
      human_ligand,
      human_rho,
      human_partial_rho_NAS,
      mouse_total_sender_ratio,
      mouse_secondary_LR_ratio,
      Mouse_Mphi_Cdcp1_UMI_share,
      Mouse_nonMphi_Cdcp1_UMI_share,
      cross_species_bridge_type,
      bridge_grade_v1.1
    )
)

cat("\n============================================================\n")
cat("CONCISE THERAPEUTIC BRIDGE HIERARCHY\n")
cat("============================================================\n")
print(concise_bridge)

cat("\n============================================================\n")
cat("CDCP1/Cdcp1 SOURCE DIVERGENCE\n")
cat("============================================================\n")
print(source_divergence)

cat("\n============================================================\n")
cat("FINAL WORKING MODEL\n")
cat("============================================================\n")
print(model)

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
