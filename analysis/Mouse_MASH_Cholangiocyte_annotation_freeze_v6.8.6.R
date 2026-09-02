suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

VERSION <- "v6.8.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.5/",
  "Mouse_MASH_Cholangiocyte_res0.4_state_audit_v6.8.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte annotation freeze\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

if (!"Chol_res04_v685" %in% colnames(obj@meta.data)) {
  stop("Missing Chol_res04_v685.")
}

cluster <- as.character(obj$Chol_res04_v685)

# ---------------------------------------------------------
# Frozen cluster annotation
# ---------------------------------------------------------

annotation_table <- data.frame(
  cluster=as.character(0:10),

  fine_state=c(
    "Pcp4l1_Serpina6_homeostatic_like",
    "Kptn_Pacs2_homeostatic_like",
    "Ccl2_Vcam1_inflammatory_reactive",
    "Msln_Aqp5_ductular_like",
    "Disease_enriched_high_complexity_QC_watch",
    "Cycling_cholangiocyte",
    "Ciliated_cholangiocyte",
    "IEG_stress_response",
    "Krt20_Cdh17_reactive_epithelial",
    "Dmbt1_Duox2_reactive_epithelial",
    "Tuft_like_cholangiocyte"
  ),

  broad_state=c(
    "Homeostatic_like",
    "Homeostatic_like",
    "Inflammatory_reactive",
    "Ductular_reactive",
    "QC_watch",
    "Cycling",
    "Ciliated",
    "Stress_response",
    "Reactive_epithelial",
    "Reactive_epithelial",
    "Tuft_like"
  ),

  analysis_class=c(
    "primary",
    "primary",
    "primary",
    "primary",
    "QC_watch_sensitivity",
    "primary",
    "primary",
    "primary",
    "exploratory_primary",
    "sample_biased_exploratory",
    "rare_exploratory"
  ),

  rationale=c(
    "Pcp4l1/Serpina6/Prox1/Acsm3 dominant conventional biliary state",
    "Kptn/Pacs2/Fam13a dominant conventional biliary state",
    "Ccl2/Vcam1/Cldn4/Ubd inflammatory-reactive program",
    "Msln/Aqp5/Anxa3/Lamb3/Gprc5a ductular epithelial program",
    "Strong STD/CDHFD enrichment; high RNA complexity; 58% derived from precleanup QC-watch cluster",
    "Dtl/Exo1/E2f1/Rad51/Mcm5 proliferation program",
    "Cfap73/Cfap44/Lrrc23/Drc1/Mns1 ciliary program",
    "Fosl1/Egr3/Atf3/Nr4a1/Fosb immediate-early/stress program",
    "Krt20/Cdh17/Inhba/Ccl2 reactive epithelial program",
    "Dmbt1/Duox2/Duoxa2/Slc26a9/Gcnt3 epithelial-reactive program; Sham20 biased",
    "Pou2f3/Trpm5 tuft-like program; extremely rare"
  ),

  stringsAsFactors=FALSE
)

write.csv(
  annotation_table,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_frozen_annotation_v6.8.6.csv"
  ),
  row.names=FALSE
)

state_map <- setNames(
  annotation_table$fine_state,
  annotation_table$cluster
)

broad_map <- setNames(
  annotation_table$broad_state,
  annotation_table$cluster
)

class_map <- setNames(
  annotation_table$analysis_class,
  annotation_table$cluster
)

obj$Chol_state_v686 <- unname(
  state_map[cluster]
)

obj$Chol_broad_state_v686 <- unname(
  broad_map[cluster]
)

obj$Chol_analysis_class_v686 <- unname(
  class_map[cluster]
)

if (anyNA(obj$Chol_state_v686)) {
  stop("NA state annotation detected.")
}

# ---------------------------------------------------------
# Explicit flags
# ---------------------------------------------------------

obj$Chol_primary_state_v686 <-
  obj$Chol_analysis_class_v686 %in%
  c(
    "primary",
    "exploratory_primary"
  )

obj$Chol_QCwatch_v686 <-
  obj$Chol_analysis_class_v686 ==
  "QC_watch_sensitivity"

obj$Chol_sample_biased_v686 <-
  obj$Chol_analysis_class_v686 ==
  "sample_biased_exploratory"

obj$Chol_rare_v686 <-
  obj$Chol_analysis_class_v686 ==
  "rare_exploratory"

# ---------------------------------------------------------
# Counts
# ---------------------------------------------------------

cat("=== FROZEN FINE STATES ===\n")
print(
  table(obj$Chol_state_v686)
)

cat("\n=== FROZEN BROAD STATES ===\n")
print(
  table(obj$Chol_broad_state_v686)
)

cat("\n=== ANALYSIS CLASS ===\n")
print(
  table(obj$Chol_analysis_class_v686)
)

# ---------------------------------------------------------
# State x sample
# ---------------------------------------------------------

state_sample <- table(
  state=obj$Chol_state_v686,
  sample=obj$sample
)

write.csv(
  as.data.frame.matrix(state_sample),
  file.path(
    TABDIR,
    "Cholangiocyte_frozen_state_by_sample_v6.8.6.csv"
  )
)

cat("\n=== FROZEN STATE x SAMPLE ===\n")
print(state_sample)

# ---------------------------------------------------------
# Within-sample state fractions
# ---------------------------------------------------------

state_sample_fraction <- prop.table(
  state_sample,
  margin=2
)

write.csv(
  as.data.frame.matrix(state_sample_fraction),
  file.path(
    TABDIR,
    "Cholangiocyte_frozen_state_fraction_within_sample_v6.8.6.csv"
  )
)

# ---------------------------------------------------------
# UMAP fine states
# ---------------------------------------------------------

fine_levels <- annotation_table$fine_state

obj$Chol_state_v686 <- factor(
  obj$Chol_state_v686,
  levels=fine_levels
)

p_fine <- DimPlot(
  obj,
  reduction="umap",
  group.by="Chol_state_v686",
  label=TRUE,
  repel=TRUE,
  pt.size=0.20
) +
  theme_classic() +
  ggtitle(
    "Mouse MASH Cholangiocyte frozen states v6.8.6"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_frozen_fine_state_UMAP_v6.8.6.pdf"
  ),
  p_fine,
  width=11,
  height=8
)

# ---------------------------------------------------------
# UMAP broad states
# ---------------------------------------------------------

p_broad <- DimPlot(
  obj,
  reduction="umap",
  group.by="Chol_broad_state_v686",
  label=TRUE,
  repel=TRUE,
  pt.size=0.20
) +
  theme_classic() +
  ggtitle(
    "Mouse MASH Cholangiocyte broad states v6.8.6"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_frozen_broad_state_UMAP_v6.8.6.pdf"
  ),
  p_broad,
  width=10,
  height=8
)

# ---------------------------------------------------------
# QC-watch location
# ---------------------------------------------------------

obj$QCwatch_display_v686 <- ifelse(
  obj$Chol_QCwatch_v686,
  "QC_watch_state",
  "Primary_or_exploratory_state"
)

p_qc <- DimPlot(
  obj,
  reduction="umap",
  group.by="QCwatch_display_v686",
  pt.size=0.20
) +
  theme_classic() +
  ggtitle(
    "Cholangiocyte QC-watch state v6.8.6"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_QCwatch_state_UMAP_v6.8.6.pdf"
  ),
  p_qc,
  width=8,
  height=7
)

# ---------------------------------------------------------
# Save frozen object
# ---------------------------------------------------------

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_annotation_frozen_v6.8.6.rds"
  ),
  compress=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte annotation freeze v6.8.6",
  "",
  paste0("- Cells: ", ncol(obj)),
  "- Clean lineage baseline: v6.8.3.2",
  "- RPCA resolution scan: v6.8.4",
  "- Frozen annotation scaffold: res0.4",
  "- Total frozen states: 11",
  "",
  "## Interpretation policy",
  "- Cluster 4 is retained as Disease_enriched_high_complexity_QC_watch.",
  "- Cluster 4 is not treated as a definitive disease-specific biological state.",
  "- Primary analyses should include a sensitivity analysis excluding the QC-watch state.",
  "- Dmbt1_Duox2_reactive_epithelial is sample-biased exploratory because Sham20 contributes strongly.",
  "- Tuft_like_cholangiocyte is rare exploratory.",
  "- STD vs CDHFD remains descriptive because each condition has one biological sample.",
  "- Sham vs Tx uses biological-sample-level inference (n=2 per group)."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_annotation_freeze_summary_v6.8.6.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.6.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.6 COMPLETE\n")
cat("Cells:", ncol(obj), "\n")
cat("Frozen resolution: res0.4\n")
cat("Frozen states:", nrow(annotation_table), "\n")
cat("QC-watch state retained with sensitivity flag\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
