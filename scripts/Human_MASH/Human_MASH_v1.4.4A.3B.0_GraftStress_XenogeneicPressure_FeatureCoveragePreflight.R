#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.3B.0_GraftStress_XenogeneicPressure_FeatureCoveragePreflight.R
#
# Purpose
#   Feature-space preflight only.
#   Determine which candidate human graft-stress / xenogeneic-pressure genes are
#   actually measurable in the 10x Flex 9,137-feature RNA space.
#
# IMPORTANT
#   - No module scoring.
#   - No differential expression.
#   - No biological absence inferred from genes absent from the Flex feature set.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

options(stringsAsFactors = FALSE)

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2A.1/objects/",
  "Human_MASH_iPSC_derived_v1.4.2A.1_hepatocyte_state_interpreted.rds"
)

OUT_ROOT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.0/",
  "analysis/v1.4.4A.3B.0_graft_stress_xenogeneic_pressure_feature_preflight"
)

TAB_DIR <- file.path(OUT_ROOT, "tables")
LOG_DIR <- file.path(OUT_ROOT, "logs")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "v1.4.4A.3B.0_run_log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================\n")
cat("Human MASH v1.4.4A.3B.0\n")
cat("Graft stress / residual xenogeneic immune-pressure\n")
cat("Feature coverage preflight\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input:", INPUT_RDS, "\n")
cat("Output:", OUT_ROOT, "\n\n")

if (!file.exists(INPUT_RDS)) stop("Input RDS not found: ", INPUT_RDS)

obj <- readRDS(INPUT_RDS)

cat("Object class:", paste(class(obj), collapse = ", "), "\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n\n")

if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")
if (!"RNA" %in% Assays(obj)) stop("RNA assay is missing.")

expected_cells <- 43415L
expected_features <- 9137L

if (ncol(obj) != expected_cells) {
  stop("Unexpected cell count: ", ncol(obj), " (expected ", expected_cells, ")")
}

rna_features <- rownames(obj[["RNA"]])

if (length(rna_features) != expected_features) {
  stop(
    "Unexpected RNA feature count: ", length(rna_features),
    " (expected ", expected_features, ")"
  )
}

if (anyDuplicated(rna_features)) stop("Duplicated feature names detected in RNA assay.")

required_meta <- c(
  "sample_v141",
  "broad_lineage_manual_v142A",
  "manual_state_v142A",
  "annotation_status_v142A"
)

missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

sample_counts <- table(obj$sample_v141)
expected_samples <- c("SfLB1", "SfLB2", "Tx17", "Tx5")

if (!setequal(names(sample_counts), expected_samples)) {
  stop("Unexpected sample set: ", paste(names(sample_counts), collapse = ", "))
}

cat("Frozen iPSC cell count: PASS\n")
cat("RNA feature-space size: PASS (", length(rna_features), " genes)\n", sep = "")
cat("\n=== SAMPLE COUNTS ===\n")
print(sample_counts)
cat("\n")

panels <- list(
  IFN_response = c(
    "STAT1","STAT2","IRF1","IRF7",
    "ISG15","IFIT1","IFIT2","IFIT3",
    "IFITM1","IFITM2","IFITM3",
    "MX1","MX2","OAS1","OAS2","OAS3",
    "IFI6","IFI27","IFI35","IFI44","IFI44L",
    "USP18","BST2","CXCL10"
  ),

  Antigen_presentation = c(
    "B2M","HLA-A","HLA-B","HLA-C",
    "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1",
    "TAP1","TAP2","TAPBP",
    "PSMB8","PSMB9","PSMB10",
    "CIITA","NLRC5"
  ),

  TNF_NFkB_response = c(
    "TNF","TNFRSF1A","TNFRSF1B",
    "NFKB1","NFKB2","RELA","RELB",
    "NFKBIA","NFKBIE","TNFAIP3",
    "BIRC3","TRAF1","ICAM1","JUN","FOS"
  ),

  Cellular_oxidative_stress = c(
    "NFE2L2","HMOX1","NQO1",
    "GCLC","GCLM","SOD1","SOD2",
    "TXN","TXNRD1","GPX1","GPX4",
    "ATF3","DDIT3",
    "HSPA1A","HSPA1B","HSP90AA1","DNAJB1",
    "JUN","FOS"
  ),

  Apoptosis_death_survival = c(
    "FAS","FASLG","CASP3","CASP8","CASP9",
    "BAX","BAK1","BID",
    "BCL2","BCL2L1","MCL1",
    "PMAIP1","BBC3"
  ),

  Hypoxia_ischemia_adaptation = c(
    "HIF1A","EPAS1","VEGFA",
    "SLC2A1","LDHA","PDK1",
    "BNIP3","BNIP3L","CA9",
    "EGLN1","EGLN3","ADM","NDRG1",
    "PGK1","ENO1","HK2"
  ),

  Graft_host_innate_interface = c(
    "CD47","B2M","CD274","LGALS9",
    "ICAM1","VCAM1",
    "CCL2","CXCL8","CX3CL1","CSF1",
    "C3","CFB","C1R","C1S","LGALS3"
  )
)

panel_long <- do.call(
  rbind,
  lapply(names(panels), function(panel_name) {
    genes <- unique(panels[[panel_name]])
    data.frame(
      panel = panel_name,
      gene = genes,
      present_in_Flex_RNA = genes %in% rna_features,
      stringsAsFactors = FALSE
    )
  })
)

panel_long$status <- ifelse(
  panel_long$present_in_Flex_RNA,
  "MEASURED_IN_FLEX",
  "NOT_IN_FLEX_FEATURE_SPACE"
)

write.csv(
  panel_long,
  file.path(TAB_DIR, "01_candidate_gene_feature_coverage_long_v1.4.4A.3B.0.csv"),
  row.names = FALSE
)

coverage_summary <- do.call(
  rbind,
  lapply(names(panels), function(panel_name) {
    z <- panel_long[panel_long$panel == panel_name, ]
    present <- z$gene[z$present_in_Flex_RNA]
    absent  <- z$gene[!z$present_in_Flex_RNA]

    data.frame(
      panel = panel_name,
      n_candidate = nrow(z),
      n_present = sum(z$present_in_Flex_RNA),
      n_absent = sum(!z$present_in_Flex_RNA),
      fraction_present = mean(z$present_in_Flex_RNA),
      present_genes = paste(present, collapse = ";"),
      absent_genes = paste(absent, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  coverage_summary,
  file.path(TAB_DIR, "02_panel_feature_coverage_summary_v1.4.4A.3B.0.csv"),
  row.names = FALSE
)

cat("=== PANEL FEATURE COVERAGE ===\n")
print(
  coverage_summary[, c(
    "panel","n_candidate","n_present","n_absent","fraction_present"
  )],
  row.names = FALSE
)

measured_only <- panel_long[
  panel_long$present_in_Flex_RNA,
  c("panel","gene")
]

write.csv(
  measured_only,
  file.path(TAB_DIR, "03_measured_only_candidate_panel_members_v1.4.4A.3B.0.csv"),
  row.names = FALSE
)

module_lines <- unlist(
  lapply(names(panels), function(panel_name) {
    present <- panel_long$gene[
      panel_long$panel == panel_name &
      panel_long$present_in_Flex_RNA
    ]
    c(
      paste0("[", panel_name, "]"),
      paste(present, collapse = ", "),
      ""
    )
  })
)

writeLines(
  module_lines,
  file.path(TAB_DIR, "04_measured_only_candidate_panels_v1.4.4A.3B.0.txt"),
  useBytes = TRUE
)

all_candidates <- sort(unique(panel_long$gene))
all_present <- sort(intersect(all_candidates, rna_features))
all_absent <- sort(setdiff(all_candidates, rna_features))

global_summary <- data.frame(
  metric = c(
    "RNA_features",
    "unique_candidate_genes",
    "unique_candidate_genes_present",
    "unique_candidate_genes_absent",
    "overall_candidate_fraction_present"
  ),
  value = c(
    length(rna_features),
    length(all_candidates),
    length(all_present),
    length(all_absent),
    length(all_present) / length(all_candidates)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  global_summary,
  file.path(TAB_DIR, "05_global_feature_coverage_summary_v1.4.4A.3B.0.csv"),
  row.names = FALSE
)

writeLines(
  all_absent,
  file.path(TAB_DIR, "06_candidate_genes_absent_from_Flex_feature_space_v1.4.4A.3B.0.txt"),
  useBytes = TRUE
)

cat("\n=== GLOBAL FEATURE COVERAGE ===\n")
print(global_summary, row.names = FALSE)

cat("\n=== CANDIDATE GENES ABSENT FROM FLEX FEATURE SPACE ===\n")
if (length(all_absent) == 0) {
  cat("None\n")
} else {
  cat(paste(all_absent, collapse = ", "), "\n")
}

note <- c(
  "Human MASH v1.4.4A.3B.0 interpretation notes",
  "",
  "Purpose:",
  "- Feature-space preflight only.",
  "- No module score, DE, pathway test, or causal inference is performed.",
  "",
  "Critical Flex rule:",
  "- Genes absent from the 9,137-feature RNA space are NOT interpreted as biologically unexpressed.",
  "- Only measured genes can be considered for subsequent transcriptomic modules.",
  "",
  "Candidate program families:",
  "- IFN response",
  "- antigen presentation",
  "- TNF/NF-kB response",
  "- cellular/oxidative stress",
  "- apoptosis/death/survival",
  "- hypoxia/ischemia adaptation",
  "- graft-host innate interface",
  "",
  "Causal boundary:",
  "- SfLB -> Tx is NOT an IL10-only perturbation.",
  "- The contrast combines transgene status, transplantation/engraftment, host microenvironment,",
  "  post-engraftment selection/state transition, and possible residual innate xenogeneic",
  "  immune pressure in NOD/scid recipients.",
  "",
  "Next step:",
  "- Review panel coverage.",
  "- Freeze a measured-gene module set for v1.4.4A.3B.",
  "- Then evaluate SfLB1/SfLB2 versus Tx17/Tx5 at sample, broad-lineage, and state levels.",
  "- Preserve replicate-level results because n=2 per group."
)

writeLines(
  note,
  file.path(OUT_ROOT, "README_interpretation_v1.4.4A.3B.0.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("v1.4.4A.3B.0 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_ROOT, "\n")
cat("Tables/files:", length(list.files(TAB_DIR)), "\n")
cat("No scoring or DE performed.\n")
cat("Genes absent from Flex feature space are not biological negatives.\n")
