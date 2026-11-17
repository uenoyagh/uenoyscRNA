#!/usr/bin/env Rscript

# ==============================================================================
# PostDecontX_ForeignLineage_FinalAudit_v1.1.R
#
# PURPOSE
# -------
# Final targeted audit after whole-liver DecontX sensitivity analysis.
#
# This script DOES NOT modify source objects and DOES NOT recluster.
# It compares original vs DecontX-corrected counts at individual-gene level for:
#
# Human:
#   cluster 3  : hepatocyte / macrophage / phagocytosis
#   cluster 13 : neutrophil / macrophage
#
# Mouse:
#   cluster 8  : hepatocyte / HSC-fibroblast / macrophage
#   cluster 15 : hepatocyte / HSC-fibroblast / macrophage
#   cluster 23 : endothelial/LSEC / macrophage
#   cluster 26 : T/NK / M1 / IFN / macrophage
#
# Outputs are intended to support KEEP / EXCLUDE / REVIEW decisions.
#
# Version: v1.1
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_MPHI <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

MOUSE_MPHI <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

DECONTX_DIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_Mphi",
  "WholeLiver_DecontX_Mphi_Sensitivity_v1.0.1",
  "corrected_Mphi_matrices"
)

HUMAN_CORR <- file.path(
  DECONTX_DIR,
  "Human_Mphi_decontX_corrected_counts_ALL_v1.0.1.rds"
)

MOUSE_CORR <- file.path(
  DECONTX_DIR,
  "Mouse_Mphi_decontX_corrected_counts_ALL_v1.0.1.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_Mphi",
  "PostDecontX_ForeignLineage_FinalAudit_v1.1"
)

TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOGDIR, "run.log")
cat("", file = LOGFILE)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOGFILE, append = TRUE)
}

write_tsv <- function(x, f) {
  write.table(
    x, f,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

for (f in c(HUMAN_MPHI, MOUSE_MPHI, HUMAN_CORR, MOUSE_CORR)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

extract_counts <- function(obj, cells) {
  layers <- Layers(obj[["RNA"]])
  clayers <- grep("^counts($|\\.)", layers, value = TRUE)
  if (!length(clayers)) stop("No RNA counts layer found.")

  pieces <- list()
  for (ly in clayers) {
    m <- LayerData(obj, assay = "RNA", layer = ly)
    keep <- intersect(cells, colnames(m))
    if (length(keep)) pieces[[ly]] <- m[, keep, drop = FALSE]
  }

  if (!length(pieces)) stop("No requested cells in RNA count layers.")

  common_genes <- Reduce(intersect, lapply(pieces, rownames))
  pieces <- lapply(pieces, function(m) m[common_genes, , drop = FALSE])
  x <- do.call(cbind, pieces)

  if (anyDuplicated(colnames(x))) stop("Duplicated cell barcodes in original counts.")
  miss <- setdiff(cells, colnames(x))
  if (length(miss)) stop("Missing original cells: ", length(miss))

  as(x[, cells, drop = FALSE], "dgCMatrix")
}

PAN_H <- list(
  Hepatocyte_ambient = c("ALB", "APOA1", "APOA2", "TTR"),
  Hepatocyte_strict = c("ASGR1", "HNF4A", "CPS1", "CYP3A4", "ALDOB", "APOB"),
  Neutrophil = c("FCGR3B", "CSF3R", "MPO", "ELANE", "CEACAM8", "OLFM4"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "CLDN5", "EGFL7", "PLVAP"),
  HSC_Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Myeloid_core = c("LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "CSF1R", "CD68"),
  Phagocytosis = c(
    "MERTK", "AXL", "TIMD4", "MARCO", "STAB1", "MSR1", "CD36",
    "C1QA", "C1QB", "C1QC", "CTSB", "CTSD", "LAMP1"
  )
)

PAN_M <- list(
  Hepatocyte_ambient = c("Alb", "Apoa1", "Apoa2", "Ttr"),
  Hepatocyte_strict = c("Asgr1", "Hnf4a", "Cps1", "Cyp3a11", "Aldob", "Apob"),
  HSC_Fibroblast = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  Endothelial_core = c("Pecam1", "Vwf", "Kdr", "Emcn", "Ramp2", "Cldn5", "Egfl7", "Plvap"),
  LSEC_extended = c(
    "Oit3", "Dll4", "Clec4g", "Cldn5", "Flt4",
    "Kdr", "F8", "Gpr182", "Ptprb", "Ramp2"
  ),
  T_cell = c("Cd3d", "Cd3e", "Trac"),
  NK = c("Klrd1", "Nkg7", "Prf1", "Gzma"),
  Myeloid_core = c("Lst1", "Tyrobp", "Fcer1g", "Ctss", "Lyz2", "Csf1r", "Cd68"),
  M1 = c("Il1b", "Tnf", "Ccl3", "Ccl4", "Cxcl2"),
  IFN = c("Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"),
  Phagocytosis = c(
    "Mertk", "Axl", "Timd4", "Marco", "Stab1", "Msr1", "Cd36",
    "C1qa", "C1qb", "C1qc", "Ctsb", "Ctsd", "Lamp1"
  )
)

gene_stats <- function(orig, corr, cluster, panels, species, focus_clusters) {

  # common ordering
  common_cells <- intersect(colnames(orig), colnames(corr))
  orig <- orig[, common_cells, drop = FALSE]
  corr <- corr[, common_cells, drop = FALSE]

  cluster <- as.character(cluster[common_cells])

  lib_o <- Matrix::colSums(orig)
  lib_c <- Matrix::colSums(corr)
  lib_o[lib_o <= 0] <- 1
  lib_c[lib_c <= 0] <- 1

  rows <- list()
  k <- 0L

  for (cl in focus_clusters) {
    ix <- which(cluster == as.character(cl))
    if (!length(ix)) next

    for (module in names(panels)) {
      genes <- intersect(panels[[module]], intersect(rownames(orig), rownames(corr)))

      for (g in genes) {
        o <- as.numeric(orig[g, ix])
        c <- as.numeric(corr[g, ix])

        o_log <- log1p(o / lib_o[ix] * 1e4)
        c_log <- log1p(c / lib_c[ix] * 1e4)

        k <- k + 1L
        rows[[k]] <- data.frame(
          species = species,
          cluster = cl,
          n_cells = length(ix),
          module = module,
          gene = g,
          pct_detected_original = mean(o > 0) * 100,
          pct_detected_corrected = mean(c > 0) * 100,
          delta_pct_detected = (mean(c > 0) - mean(o > 0)) * 100,
          mean_logCP10k_original = mean(o_log),
          mean_logCP10k_corrected = mean(c_log),
          delta_mean_logCP10k = mean(c_log) - mean(o_log),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows)
}

module_stats <- function(gene_df) {
  spl <- split(
    gene_df,
    interaction(gene_df$species, gene_df$cluster, gene_df$module, drop = TRUE)
  )

  rr <- lapply(spl, function(x) {
    data.frame(
      species = x$species[1],
      cluster = x$cluster[1],
      module = x$module[1],
      n_genes = nrow(x),
      mean_pct_detected_original = mean(x$pct_detected_original, na.rm = TRUE),
      mean_pct_detected_corrected = mean(x$pct_detected_corrected, na.rm = TRUE),
      delta_mean_pct_detected = mean(x$delta_pct_detected, na.rm = TRUE),
      mean_expression_original = mean(x$mean_logCP10k_original, na.rm = TRUE),
      mean_expression_corrected = mean(x$mean_logCP10k_corrected, na.rm = TRUE),
      delta_mean_expression = mean(x$delta_mean_logCP10k, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rr)
}

msg("Loading Human...")
human <- readRDS(HUMAN_MPHI)
hcorr <- readRDS(HUMAN_CORR)
horig <- extract_counts(human, colnames(human))

if (!all(colnames(human) %in% colnames(hcorr))) {
  stop("Human corrected matrix does not contain all frozen Mphi cells.")
}
hcorr <- hcorr[, colnames(human), drop = FALSE]

hcluster <- human@meta.data[["res0_8_cluster_frozen_v18"]]
names(hcluster) <- colnames(human)

hgene <- gene_stats(
  horig, hcorr, hcluster, PAN_H,
  species = "Human",
  focus_clusters = c("3", "13")
)

write_tsv(
  hgene,
  file.path(TABDIR, "01_Human_focus_gene_original_vs_DecontX_v1.1.tsv")
)

write_tsv(
  module_stats(hgene),
  file.path(TABDIR, "02_Human_focus_foreign_module_original_vs_DecontX_v1.1.tsv")
)

msg("Loading Mouse...")
mouse <- readRDS(MOUSE_MPHI)
mcorr <- readRDS(MOUSE_CORR)
morig <- extract_counts(mouse, colnames(mouse))

if (!all(colnames(mouse) %in% colnames(mcorr))) {
  stop("Mouse corrected matrix does not contain all frozen Mphi cells.")
}
mcorr <- mcorr[, colnames(mouse), drop = FALSE]

mcluster <- mouse@meta.data[["cluster_res2"]]
names(mcluster) <- colnames(mouse)

mgene <- gene_stats(
  morig, mcorr, mcluster, PAN_M,
  species = "Mouse",
  focus_clusters = c("8", "15", "23", "26")
)

write_tsv(
  mgene,
  file.path(TABDIR, "03_Mouse_focus_gene_original_vs_DecontX_v1.1.tsv")
)

write_tsv(
  module_stats(mgene),
  file.path(TABDIR, "04_Mouse_focus_foreign_module_original_vs_DecontX_v1.1.tsv")
)

# Compact decision-support table based only on measured modules.
# This table deliberately leaves final disposition to review.
decision <- data.frame(
  species = c("Human", "Human", "Mouse", "Mouse", "Mouse", "Mouse"),
  cluster = c("3", "13", "8", "15", "23", "26"),
  pre_audit_status = c(
    "KEEP_review_ambient_passenger",
    "EXCLUDE_candidate_neutrophil",
    "KEEP_review_hepatocyte_HSC",
    "KEEP_review_hepatocyte",
    "EXCLUDE_or_REVIEW_endothelial_doublet",
    "KEEP_M1_IFN_state"
  ),
  key_question = c(
    "Do hepatocyte genes fall while macrophage/phagocytosis genes persist?",
    "Do neutrophil genes persist after correction while macrophage core remains low?",
    "Do hepatocyte/HSC genes fall while macrophage genes persist?",
    "Do hepatocyte genes fall while macrophage genes persist?",
    "Do endothelial/LSEC genes persist after correction together with macrophage genes?",
    "Do T/NK genes remain low while M1/IFN/macrophage genes persist?"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  decision,
  file.path(TABDIR, "05_final_decision_questions_v1.1.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt")
)

msg("STATUS: SUCCESS")
msg("Review first:")
msg("  01_Human_focus_gene_original_vs_DecontX_v1.1.tsv")
msg("  02_Human_focus_foreign_module_original_vs_DecontX_v1.1.tsv")
msg("  03_Mouse_focus_gene_original_vs_DecontX_v1.1.tsv")
msg("  04_Mouse_focus_foreign_module_original_vs_DecontX_v1.1.tsv")
msg("  05_final_decision_questions_v1.1.tsv")
