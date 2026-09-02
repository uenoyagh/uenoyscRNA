suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v6.8.10"

BASE <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/"
)

FILE_STATE <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_v6.8.7/tables/",
  "Cholangiocyte_state_axis_summary_NO_QCWATCH_v6.8.7.csv"
)

FILE_MODULE <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_v6.8.7/tables/",
  "Cholangiocyte_module_axis_NO_QCWATCH_v6.8.7.csv"
)

FILE_DE <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_v6.8.8.1/tables/",
  "PRIMARY_CORE_FDRlt0.10_reference_annotated_v6.8.8.1.csv"
)

FILE_GLOBAL_GSEA <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_v6.8.9/tables/",
  "Hallmark_GSEA_PRIMARY_CORE_STRICT_REFERENCE_SENSITIVITY_v6.8.9.csv"
)

FILE_STATE_GSEA <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_v6.8.9.1/tables/",
  "State_specific_Hallmark_GSEA_ALL_RESULTS_v6.8.9.1.csv"
)

OUTDIR <- paste0(
  BASE,
  "Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte mechanistic closure\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

files <- c(
  FILE_STATE,
  FILE_MODULE,
  FILE_DE,
  FILE_GLOBAL_GSEA,
  FILE_STATE_GSEA
)

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "Missing inputs:\n",
    paste(missing_files, collapse="\n")
  )
}

state <- read.csv(
  FILE_STATE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

module <- read.csv(
  FILE_MODULE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

de <- read.csv(
  FILE_DE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

global_gsea <- read.csv(
  FILE_GLOBAL_GSEA,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

state_gsea <- read.csv(
  FILE_STATE_GSEA,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

# =========================================================
# 1. Key state-composition effects
# =========================================================

state_key <- state[
  order(
    -abs(state$Tx_minus_Sham)
  ),
  ,
  drop=FALSE
]

write.csv(
  state_key,
  file.path(
    TABDIR,
    "01_state_composition_Tx_vs_Sham_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 2. Module effects
# =========================================================

module_key <- module[
  order(
    -abs(module$Tx_minus_Sham)
  ),
  ,
  drop=FALSE
]

write.csv(
  module_key,
  file.path(
    TABDIR,
    "02_module_Tx_vs_Sham_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 3. Ambient-aware DE candidates
# =========================================================

de_primary <- de[
  de$Hepatocyte_specificity_class !=
    "Hepatocyte_dominant",
  ,
  drop=FALSE
]

de_primary <- de_primary[
  order(
    de_primary$FDR,
    -abs(de_primary$logFC)
  ),
  ,
  drop=FALSE
]

write.csv(
  de_primary,
  file.path(
    TABDIR,
    "03_PRIMARY_CORE_FDRlt0.10_ambient_aware_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 4. Global strict-reference Hallmark
# =========================================================

selected_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)

global_key <- global_gsea[
  global_gsea$pathway %in% selected_pathways,
  ,
  drop=FALSE
]

global_key <- global_key[
  order(
    global_key$NES,
    decreasing=TRUE
  ),
  ,
  drop=FALSE
]

write.csv(
  global_key,
  file.path(
    TABDIR,
    "04_global_strict_reference_key_Hallmarks_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 5. State-specific strict-reference Hallmark
# =========================================================

primary_states <- c(
  "Homeostatic_combined",
  "Inflammatory_reactive",
  "Ductular_like",
  "Krt20_Cdh17_reactive",
  "Cycling",
  "Ciliated"
)

state_key_gsea <- state_gsea[
  state_gsea$filter_mode ==
    "STRICT_REFERENCE_SENSITIVITY" &
  state_gsea$state %in%
    primary_states &
  state_gsea$pathway %in%
    selected_pathways,
  ,
  drop=FALSE
]

write.csv(
  state_key_gsea,
  file.path(
    TABDIR,
    "05_state_specific_strict_reference_key_Hallmarks_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 6. Final state-specific heatmap
# =========================================================

state_key_gsea$state <- factor(
  state_key_gsea$state,
  levels=primary_states
)

pathway_labels <- sub(
  "^HALLMARK_",
  "",
  selected_pathways
)

pathway_labels <- gsub(
  "_",
  " ",
  pathway_labels
)

state_key_gsea$pathway_label <- sub(
  "^HALLMARK_",
  "",
  state_key_gsea$pathway
)

state_key_gsea$pathway_label <- gsub(
  "_",
  " ",
  state_key_gsea$pathway_label
)

state_key_gsea$pathway_label <- factor(
  state_key_gsea$pathway_label,
  levels=rev(pathway_labels)
)

p <- ggplot(
  state_key_gsea,
  aes(
    x=state,
    y=pathway_label,
    fill=NES
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low="#0033FF",
    mid="#FFFFFF",
    high="#FF1A1A",
    midpoint=0
  ) +
  theme_classic(base_size=9) +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    )
  ) +
  labs(
    title=
      "Mouse MASH Cholangiocyte mechanistic closure",
    subtitle=
      "Tx vs Sham | strict reference sensitivity | state-specific pseudobulk GSEA",
    x=NULL,
    y=NULL,
    fill="NES"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_FINAL_state_specific_mechanistic_heatmap_v6.8.10.pdf"
  ),
  p,
  width=12,
  height=8
)

# =========================================================
# 7. Evidence table
# =========================================================

evidence <- data.frame(
  evidence_level=c(
    "Lineage cleanup",
    "Annotation",
    "QC sensitivity",
    "State composition",
    "Module score",
    "Pseudobulk DE",
    "Global pathway",
    "State-specific pathway",
    "Mechanistic interpretation"
  ),

  conclusion=c(
    "15,755 lineage-clean Cholangiocytes retained after conservative cleanup.",
    "11 frozen res0.4 states; rare/sample-biased/QC-watch states explicitly flagged.",
    "Excluding high-complexity QC-watch cells does not materially change the Tx-axis conclusions.",
    "Tx increases selected reactive/ductular states rather than producing uniform homeostatic normalization.",
    "Tx increases Stress/IEG, biliary identity, inflammatory-adhesion and Krt20/Cdh17-reactive programs.",
    "Ccn1, Jun and Actg1 remain Tx-up after Hepatocyte ambient-risk review; Slc25a47 and Gdf15 are Hepatocyte-dominant and excluded from primary interpretation.",
    "TNFA/NFKB, hypoxia, p53/apoptosis and EMT/remodeling are Tx-enriched after strict reference filtering; E2F/cell-cycle programs decrease.",
    "TNFA/NFKB and hypoxia responses recur across multiple Cholangiocyte states, including Homeostatic_combined.",
    "Tx is associated with broad Cholangiocyte epithelial reprogramming/adaptive-reactive remodeling rather than simple normalization."
  ),

  limitation=c(
    "Conservative parent definition; Cycling-derived rescue rejected.",
    "Rare Tuft-like and Sham20-biased Dmbt1/Duox2 states remain exploratory.",
    "QC-watch state is retained only with explicit sensitivity analysis.",
    "Sham/Tx n=2/group; composition differences are effect-size observations.",
    "AddModuleScore is supportive and not an independent statistical test.",
    "edgeR FDR interpreted cautiously at n=2/group.",
    "GSEA FDR does not overcome the small biological-sample number.",
    "State-specific analyses have variable cell numbers; IEG-stress remains supplemental.",
    "Association does not establish whether the remodeling is beneficial, harmful, or directly caused by graft-derived signaling."
  ),

  stringsAsFactors=FALSE
)

write.csv(
  evidence,
  file.path(
    TABDIR,
    "06_mechanistic_closure_evidence_table_v6.8.10.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 8. Final closure document
# =========================================================

summary_lines <- c(
  "# Mouse MASH Cholangiocyte mechanistic closure v6.8.10",
  "",
  "## Frozen analytical baseline",
  "- Lineage-clean Cholangiocytes: 15,755 cells.",
  "- Final annotation scaffold: RPCA res0.4.",
  "- Frozen states: 11.",
  "- Disease-enriched high-complexity cluster retained as QC-watch with exclusion sensitivity analysis.",
  "",
  "## Disease axis",
  "- STD vs CDHFD is descriptive only because each group contains one biological sample.",
  "- Disease-associated differences are therefore not used for inferential claims.",
  "",
  "## Sham to Tx composition",
  "- Tx does not produce uniform return to a homeostatic Cholangiocyte composition.",
  "- Selected reactive/ductular epithelial populations increase while Cycling/Ciliated and some other states decrease.",
  "- QC-watch exclusion does not materially alter the Tx-axis interpretation.",
  "",
  "## Sham to Tx transcriptional response",
  "- Stress/IEG, inflammatory-adhesion and reactive epithelial programs increase with Tx.",
  "- Biliary identity also increases modestly, arguing against simple loss of Cholangiocyte identity.",
  "- Ambient-aware pseudobulk supports Ccn1, Jun and Actg1 as Tx-up candidates.",
  "- Slc25a47 and Gdf15 are Hepatocyte-dominant and are not used as primary Cholangiocyte evidence.",
  "",
  "## Pathway-level response",
  "- Strict reference-aware GSEA identifies strong Tx-associated TNFA/NFKB and hypoxia activation.",
  "- p53/apoptosis and epithelial-mesenchymal/remodeling programs are also increased.",
  "- E2F and related proliferative programs are reduced globally.",
  "- These findings remain after Hepatocyte and broader reference-lineage sensitivity filtering.",
  "",
  "## State-specific validation",
  "- TNFA/NFKB and hypoxia activation recur across multiple frozen Cholangiocyte states.",
  "- Importantly, Homeostatic_combined Cholangiocytes show the same response.",
  "- Therefore the global pathway signal cannot be explained solely by redistribution toward reactive states.",
  "",
  "## Mechanistic closure",
  "Tx is associated with broad Cholangiocyte epithelial reprogramming rather than simple normalization.",
  "The response combines preserved/increased biliary identity with NF-kB/IEG, hypoxic-stress, p53/apoptotic and EMT/remodeling programs.",
  "The most appropriate interpretation is adaptive/reactive Cholangiocyte remodeling after transplantation.",
  "",
  "## Limits",
  "- Sham vs Tx: n=2 biological samples/group.",
  "- STD vs CDHFD: n=1/group and descriptive only.",
  "- Pathway FDR values should not be interpreted as overcoming the small biological-sample number.",
  "- The analysis does not establish whether Cholangiocyte remodeling is beneficial or detrimental.",
  "- Direct graft-to-Cholangiocyte signaling requires later ligand-receptor/cross-species analysis."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Mouse_MASH_Cholangiocyte_mechanistic_closure_v6.8.10.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.10.txt"
  )
)

cat("\n=== FINAL AMBIENT-AWARE DE ===\n")
print(
  de_primary[
    ,
    intersect(
      c(
        "gene",
        "logFC",
        "FDR",
        "replicate_pattern",
        "Hepatocyte_specificity_class"
      ),
      colnames(de_primary)
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== FINAL GLOBAL HALLMARKS ===\n")
print(
  global_key[
    ,
    c(
      "pathway",
      "NES",
      "padj"
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.8.10 COMPLETE\n")
cat("Cholangiocyte mechanistic closure complete\n")
cat("No cells changed\n")
cat("No annotation changed\n")
cat("No new inferential test added\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
