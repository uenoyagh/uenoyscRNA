suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v6.9.10"

BASE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS"
)

V6972_DIR <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.7.2",
  "tables"
)

V698_DIR <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.8",
  "tables"
)

V6991_DIR <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.9.1",
  "tables"
)

OUTDIR <- file.path(
  BASE_DIR,
  paste0("Mouse_MASH_Monocyte_", VERSION)
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte mechanistic closure\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# =========================================================
# Inputs
# =========================================================

GLOBAL_REF_FILE <- file.path(
  V6972_DIR,
  "Monocyte_Hallmark_GSEA_PRIMARY_v6.9.7.2.csv"
)

GLOBAL_STRICT_FILE <- file.path(
  V6972_DIR,
  "Monocyte_Hallmark_GSEA_PRIMARY_strict_sensitivity_v6.9.7.2.csv"
)

STATE_SUPPORT_FILE <- file.path(
  V698_DIR,
  "Monocyte_state_specific_support_policy_v6.9.8.csv"
)

STATE_CONCORDANCE_FILE <- file.path(
  V698_DIR,
  "Monocyte_state_specific_Hallmark_concordance_v6.9.8.csv"
)

CORE_DRIVER_FILE <- file.path(
  V6991_DIR,
  "Monocyte_CORE_ROBUST_leading_edge_drivers_v6.9.9.1.csv"
)

CROSS_CORE_FILE <- file.path(
  V6991_DIR,
  "Monocyte_cross_core_pathway_recurrent_genes_v6.9.9.1.csv"
)

input_files <- c(
  GLOBAL_REF_FILE,
  GLOBAL_STRICT_FILE,
  STATE_SUPPORT_FILE,
  STATE_CONCORDANCE_FILE,
  CORE_DRIVER_FILE,
  CROSS_CORE_FILE
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse="\n")
  )
}

global_ref <- read.csv(
  GLOBAL_REF_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

global_strict <- read.csv(
  GLOBAL_STRICT_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

state_support <- read.csv(
  STATE_SUPPORT_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

state_concordance <- read.csv(
  STATE_CONCORDANCE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

core_drivers <- read.csv(
  CORE_DRIVER_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cross_core <- read.csv(
  CROSS_CORE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

# =========================================================
# Core pathway closure
# =========================================================

core_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_P53_PATHWAY"
)

get_global <- function(df, pathway, value_col) {

  x <- df[
    df$pathway == pathway,
    value_col,
    drop=TRUE
  ]

  if (length(x) == 0) {
    return(NA)
  }

  x[[1]]
}

get_state <- function(
  df,
  pathway,
  filter_name,
  value_col
) {

  x <- df[
    df$pathway == pathway &
    df$reference_filter == filter_name,
    value_col,
    drop=TRUE
  ]

  if (length(x) == 0) {
    return(NA)
  }

  x[[1]]
}

core_summary_rows <- list()

for (pw in core_pathways) {

  core_summary_rows[[pw]] <- data.frame(
    pathway=pw,

    global_reference_NES=
      get_global(
        global_ref,
        pw,
        "NES"
      ),

    global_reference_FDR=
      get_global(
        global_ref,
        pw,
        "padj"
      ),

    global_strict_NES=
      get_global(
        global_strict,
        pw,
        "NES"
      ),

    global_strict_FDR=
      get_global(
        global_strict,
        pw,
        "padj"
      ),

    state_reference_n_states=
      get_state(
        state_concordance,
        pw,
        "REFERENCE_AWARE",
        "n_states"
      ),

    state_reference_n_Tx_enriched=
      get_state(
        state_concordance,
        pw,
        "REFERENCE_AWARE",
        "n_Tx_enriched"
      ),

    state_reference_n_Tx_FDRlt0.10=
      get_state(
        state_concordance,
        pw,
        "REFERENCE_AWARE",
        "n_Tx_FDRlt0.10"
      ),

    state_reference_n_Tx_FDRlt0.25=
      get_state(
        state_concordance,
        pw,
        "REFERENCE_AWARE",
        "n_Tx_FDRlt0.25"
      ),

    state_strict_n_Tx_enriched=
      get_state(
        state_concordance,
        pw,
        "STRICT",
        "n_Tx_enriched"
      ),

    state_strict_n_Tx_FDRlt0.10=
      get_state(
        state_concordance,
        pw,
        "STRICT",
        "n_Tx_FDRlt0.10"
      ),

    state_strict_n_Tx_FDRlt0.25=
      get_state(
        state_concordance,
        pw,
        "STRICT",
        "n_Tx_FDRlt0.25"
      ),

    stringsAsFactors=FALSE
  )
}

core_summary <- do.call(
  rbind,
  core_summary_rows
)

rownames(core_summary) <- NULL

core_summary$closure_class <- ifelse(
  core_summary$global_reference_FDR < 0.10 &
  core_summary$global_strict_FDR < 0.10 &
  core_summary$state_reference_n_Tx_enriched ==
    core_summary$state_reference_n_states &
  core_summary$state_strict_n_Tx_enriched >= 5,
  "CORE_ROBUST_TX_PROGRAM",
  "SUPPORTING_PROGRAM"
)

write.csv(
  core_summary,
  file.path(
    TABDIR,
    "Monocyte_mechanistic_closure_core_pathways_v6.9.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Eligible-state summary
# =========================================================

eligible_states <- state_support[
  state_support$state_specific_eligible,
  ,
  drop=FALSE
]

write.csv(
  eligible_states,
  file.path(
    TABDIR,
    "Monocyte_mechanistic_closure_eligible_states_v6.9.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Recurrent leading-edge network
#
# This is a pathway-level recurrent network.
# Individual genes are NOT declared statistically significant
# simply because they recur in GSEA leading edges.
# =========================================================

cross_core$network_tier <- ifelse(
  cross_core$n_core_pathways >= 3,
  "CROSS_CORE_3_PATHWAYS",
  ifelse(
    cross_core$n_core_pathways == 2,
    "CROSS_CORE_2_PATHWAYS",
    "SINGLE_CORE_PATHWAY"
  )
)

cross_core$replicate_support <- ifelse(
  cross_core$PRIMARY_CORE_replicate_direction ==
    "Tx_both_higher",
  "Tx_both_higher",
  "not_concordant"
)

cross_core <- cross_core[
  order(
    factor(
      cross_core$network_tier,
      levels=c(
        "CROSS_CORE_3_PATHWAYS",
        "CROSS_CORE_2_PATHWAYS",
        "SINGLE_CORE_PATHWAY"
      )
    ),
    -cross_core$max_states_in_core_pathway,
    -abs(cross_core$PRIMARY_CORE_logFC)
  ),
  ,
  drop=FALSE
]

write.csv(
  cross_core,
  file.path(
    TABDIR,
    "Monocyte_mechanistic_closure_recurrent_network_v6.9.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Interpretation tiers
# =========================================================

interpretation <- data.frame(
  level=c(
    "COMPOSITION",
    "CORE_TRANSCRIPTIONAL_PROGRAM",
    "CORE_TRANSCRIPTIONAL_PROGRAM",
    "CORE_TRANSCRIPTIONAL_PROGRAM",
    "RECURRENT_LEADING_EDGE_NETWORK",
    "SUPPORTING_PROGRAMS",
    "SENSITIVITY_LIMIT",
    "MYELOID_CONTEXT",
    "NOT_ESTABLISHED"
  ),

  finding=c(
    "Tx shows modest Monocyte state redistribution, including enrichment of the Cd300e/Pglyrp1/Cd36/S1pr5 activated state; this is not a broad anti-inflammatory normalization.",
    "TNFA/NFKB signaling is strongly Tx-enriched globally and within essentially all adequately sampled primary Monocyte states.",
    "Hypoxia is strongly Tx-enriched globally and across all adequately sampled primary Monocyte states.",
    "p53/stress signaling is Tx-enriched globally and directionally concordant across all adequately sampled primary Monocyte states.",
    "Repeated leading-edge genes link TNFA/NFKB, hypoxia and p53/stress programs; Fos and Atf3 span all three core pathways, while Dusp1, Zfp36, Nfil3, Sat1, Ddit4, Btg2, Foxo3 and Klf6 recur across two core pathways.",
    "UV/stress response, apoptosis, IL2/STAT5 and KRAS-up support a broader reactive/adaptive program. TGF-beta, EMT-like remodeling, adipogenesis and metabolic pathways are secondary because they are more reference-filter sensitive.",
    "Most individual genes are not significant after gene-level FDR correction with n=2/group; recurrent leading-edge membership is therefore pathway-level support, not proof of a causal driver gene.",
    "The Monocyte response differs from the mature macrophage compartment, where Tx is associated with anti-inflammatory and repair/resolution remodeling. The combined result is compatible with heterogeneous myeloid reorganization rather than uniform anti-inflammatory conversion.",
    "Direct Monocyte-to-macrophage differentiation, causal ligand-receptor signaling, and a causal role for individual leading-edge genes are not established by these analyses."
  ),

  stringsAsFactors=FALSE
)

write.csv(
  interpretation,
  file.path(
    TABDIR,
    "Monocyte_mechanistic_closure_interpretation_v6.9.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Compact figure: global pathway strength + state robustness
# =========================================================

plot_df <- core_summary

plot_df$pathway_display <- c(
  "TNF/NF-kB",
  "Hypoxia",
  "p53 / stress"
)[
  match(
    plot_df$pathway,
    core_pathways
  )
]

plot_df$pathway_display <- factor(
  plot_df$pathway_display,
  levels=rev(
    c(
      "TNF/NF-kB",
      "Hypoxia",
      "p53 / stress"
    )
  )
)

p <- ggplot(
  plot_df,
  aes(
    x=global_reference_NES,
    y=pathway_display,
    size=state_reference_n_Tx_FDRlt0.10
  )
) +
  geom_point() +
  geom_vline(
    xintercept=0,
    linetype=2
  ) +
  theme_classic(
    base_size=11
  ) +
  labs(
    title="Mouse Monocyte Tx core transcriptional programs",
    subtitle="Global PRIMARY_CORE GSEA; point size = states with Tx enrichment at FDR < 0.10",
    x="Global Hallmark NES (Tx vs Sham)",
    y=NULL,
    size="States FDR<0.10"
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_mechanistic_closure_core_pathways_v6.9.10.pdf"
  ),
  p,
  width=7,
  height=4.5
)

# =========================================================
# Plain-text closure report
# =========================================================

report_file <- file.path(
  OUTDIR,
  "Monocyte_mechanistic_closure_v6.9.10.txt"
)

report_lines <- c(
  "Mouse MASH Monocyte mechanistic closure v6.9.10",
  "",
  "Primary conclusion:",
  "Tx does not produce a simple anti-inflammatory normalization of the Monocyte compartment.",
  "Instead, Tx is associated with modest state redistribution plus broad within-state reactive/stress-adaptive transcriptional reprogramming.",
  "",
  "Core robust programs:",
  "- TNFA/NFKB",
  "- Hypoxia",
  "- p53/stress",
  "",
  "Recurrent leading-edge network:",
  "- Cross-core 3 pathways: Fos, Atf3",
  "- Cross-core recurrent genes include Dusp1, Zfp36, Nfil3, Sat1, Ddit4, Btg2, Foxo3 and Klf6.",
  "",
  "Interpretive limit:",
  "Individual leading-edge genes are not treated as independently significant causal drivers because gene-level FDR is weak with Sham/Tx n=2/group.",
  "",
  "Myeloid context:",
  "The Monocyte response is distinct from the mature macrophage repair/resolution shift and supports heterogeneous Tx-associated myeloid remodeling.",
  "",
  "Not established:",
  "- Direct Monocyte-to-macrophage differentiation",
  "- Direct ligand-receptor causality",
  "- Causal action of individual leading-edge genes"
)

writeLines(
  report_lines,
  con=report_file
)

# =========================================================
# Terminal output
# =========================================================

cat("\n=== CORE PATHWAY CLOSURE ===\n")
print(
  core_summary,
  row.names=FALSE
)

cat("\n=== RECURRENT CORE NETWORK ===\n")

network_cols <- c(
  "gene",
  "n_core_pathways",
  "max_states_in_core_pathway",
  "pathway",
  "PRIMARY_CORE_logFC",
  "PRIMARY_CORE_FDR",
  "PRIMARY_CORE_replicate_direction",
  "network_tier"
)

print(
  head(
    cross_core[
      ,
      intersect(
        network_cols,
        colnames(cross_core)
      ),
      drop=FALSE
    ],
    30
  ),
  row.names=FALSE
)

cat("\n=== FINAL MONOCYTE INTERPRETATION ===\n")

for (i in seq_len(nrow(interpretation))) {
  cat(
    "[",
    interpretation$level[[i]],
    "] ",
    interpretation$finding[[i]],
    "\n",
    sep=""
  )
}

cat("\n====================================================\n")
cat("v6.9.10 COMPLETE\n")
cat("Mouse Monocyte mechanistic closure complete\n")
cat("No new DE or GSEA testing performed\n")
cat("No source object modified\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.10.txt"
  )
)
