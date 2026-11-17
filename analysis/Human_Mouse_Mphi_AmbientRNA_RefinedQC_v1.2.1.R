#!/usr/bin/env Rscript

# ==============================================================================
# Human_Mouse_Mphi_AmbientRNA_RefinedQC_v1.2.1.R
#
# Refines v1.2.0 ambient/passenger-RNA QC:
# - Uses sample-wise robust z-scores to reduce between-sample ambient differences.
# - Requires >=2 ambient-prone hepatocyte genes AND elevated sample-wise score.
# - Produces gene-level detection/mean-expression tables for focus clusters.
# - Does not modify source RDS objects.
#
# Focus:
# Human clusters 3, 13
# Mouse clusters 8, 15, 23, 26
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)
set.seed(20261118)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_SOURCE <- file.path(
  ROOT, "Human_MASH_RDS", "Human_MASH_LT1to8_Myeloid_RPCA_v1.7",
  "objects", "Human_MASH_LT1to8_Myeloid_RPCA_v1.7.rds"
)

HUMAN_FINAL <- file.path(
  ROOT, "Human_MASH_RDS", "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects", "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

MOUSE_FINAL <- file.path(
  ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS", "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTDIR <- file.path(
  ROOT, "CrossSpecies_MASH_Mphi",
  "Human_Mouse_Mphi_AmbientRNA_RefinedQC_v1.2.1"
)

TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(LOGDIR, "run.log")
cat("", file = LOG)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOG, append = TRUE)
}

write_tsv <- function(x, f) {
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

first_existing <- function(x, candidates) {
  z <- candidates[candidates %in% x]
  if (length(z)) z[1] else NA_character_
}

normalize_mouse_sample <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("STD", "STD_rep1")] <- "STD"
  x[x %in% c("CDHFD", "CDHFD_rep1", "CDAHFD", "CDAHFD_rep1")] <- "CDAHFD"
  x
}

z_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s
}

sample_z <- function(x, sample) {
  out <- rep(NA_real_, length(x))
  for (s in unique(sample)) {
    ix <- which(sample == s)
    out[ix] <- z_safe(x[ix])
  }
  out
}

safe_mad_high <- function(x, k = 3) {
  med <- median(x, na.rm = TRUE)
  md <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(md) || md == 0) rep(FALSE, length(x)) else x > med + k * md
}

for (f in c(HUMAN_SOURCE, HUMAN_FINAL, MOUSE_FINAL)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

human <- readRDS(HUMAN_FINAL)
hsrc <- readRDS(HUMAN_SOURCE)
mouse <- readRDS(MOUSE_FINAL)

HCL <- "res0_8_cluster_frozen_v18"
MCL <- "cluster_res2"

human_sample <- as.character(human$sample_v17)

mouse_sample_col <- first_existing(
  colnames(mouse@meta.data),
  c("sample_display_FIXED2", "sample_for_R8plot_FIXED2",
    "sample_for_annotation", "sample", "sample_internal",
    "sample_4group", "orig.ident")
)

if (is.na(mouse_sample_col)) stop("Mouse sample column not found.")
mouse_sample <- normalize_mouse_sample(mouse@meta.data[[mouse_sample_col]])

PANELS_H <- list(
  Ambient = c("ALB", "APOA1", "APOA2", "TTR"),
  Passenger = c("ALB", "APOA1", "CYP3A4"),
  HepLineage = c("ASGR1", "HNF4A", "CPS1", "ALDOB", "APOB", "CYP3A4"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "CSF1R", "CD68"),
  Phago = c("MERTK", "AXL", "TIMD4", "MARCO", "STAB1", "MSR1", "CD36",
            "C1QA", "C1QB", "C1QC", "CTSB", "CTSD", "LAMP1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "CLDN5", "EGFL7", "PLVAP"),
  HSC_Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Neutrophil = c("FCGR3B", "CSF3R", "MPO", "ELANE", "CEACAM8", "OLFM4"),
  T_cell = c("CD3D", "CD3E", "TRAC")
)

PANELS_M <- list(
  Ambient = c("Alb", "Apoa1", "Apoa2", "Ttr"),
  Passenger = c("Alb", "Apoa1", "Cyp3a11"),
  HepLineage = c("Asgr1", "Hnf4a", "Cps1", "Aldob", "Apob", "Cyp3a11"),
  Myeloid = c("Lst1", "Tyrobp", "Fcer1g", "Ctss", "Lyz2", "Csf1r", "Cd68"),
  Phago = c("Mertk", "Axl", "Timd4", "Marco", "Stab1", "Msr1", "Cd36",
            "C1qa", "C1qb", "C1qc", "Ctsb", "Ctsd", "Lamp1"),
  Endothelial = c("Pecam1", "Vwf", "Kdr", "Emcn", "Ramp2", "Cldn5", "Egfl7", "Plvap"),
  HSC_Fibroblast = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  Neutrophil = c("Csf3r", "Mpo", "Elane", "Camp", "Ltf", "Ngp", "Retnlg"),
  T_cell = c("Cd3d", "Cd3e", "Trac")
)

score_and_gene_stats <- function(obj, sample, cluster, panels, species, focus_clusters) {

  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  genes <- unique(unlist(panels, use.names = FALSE))
  genes <- intersect(genes, rownames(counts))

  sub <- counts[genes, , drop = FALSE]
  lib <- Matrix::colSums(counts)
  lib[lib <= 0] <- 1

  logcp10k <- log1p(t(t(as.matrix(sub)) / lib) * 1e4)

  score_panel <- function(gs) {
    g <- intersect(gs, rownames(logcp10k))
    list(
      score = colMeans(logcp10k[g, , drop = FALSE]),
      detected = colSums(as.matrix(sub[g, , drop = FALSE]) > 0)
    )
  }

  scores <- data.frame(
    cell = colnames(obj),
    sample = sample,
    cluster = as.character(cluster),
    stringsAsFactors = FALSE
  )

  for (nm in names(panels)) {
    z <- score_panel(panels[[nm]])
    scores[[paste0(nm, "_score")]] <- z$score
    scores[[paste0(nm, "_detected")]] <- z$detected
    scores[[paste0(nm, "_z_sample")]] <- sample_z(z$score, sample)
  }

  ncount <- if ("nCount_RNA" %in% colnames(obj@meta.data)) {
    as.numeric(obj$nCount_RNA)
  } else {
    as.numeric(Matrix::colSums(counts))
  }

  nfeat <- if ("nFeature_RNA" %in% colnames(obj@meta.data)) {
    as.numeric(obj$nFeature_RNA)
  } else {
    as.numeric(Matrix::colSums(counts > 0))
  }

  scores$high_nCount <- FALSE
  scores$high_nFeature <- FALSE

  for (s in unique(sample)) {
    ix <- which(sample == s)
    scores$high_nCount[ix] <- safe_mad_high(ncount[ix])
    scores$high_nFeature[ix] <- safe_mad_high(nfeat[ix])
  }

  # Refined flags: sample-relative elevation is mandatory.
  scores$ambient_elevated_refined <- (
    scores$Ambient_detected >= 2 &
      scores$Ambient_z_sample >= 1 &
      scores$Myeloid_z_sample >= -0.5
  )

  scores$passenger_RNA_compatible_refined <- (
    scores$Passenger_detected >= 2 &
      scores$Passenger_z_sample >= 0.5 &
      scores$Phago_z_sample >= 0.5 &
      scores$Myeloid_z_sample >= -0.5 &
      !scores$high_nCount &
      !scores$high_nFeature
  )

  scores$hep_lineage_elevated_refined <- (
    scores$HepLineage_detected >= 2 &
      scores$HepLineage_z_sample >= 1.5
  )

  scores$possible_hepatocyte_doublet_refined <- (
    scores$hep_lineage_elevated_refined &
      (
        scores$Myeloid_z_sample < -0.5 |
          scores$high_nCount |
          scores$high_nFeature
      )
  )

  scores$endothelial_strong_refined <- (
    scores$Endothelial_detected >= 2 &
      scores$Endothelial_z_sample >= 1.5
  )

  scores$neutrophil_strong_refined <- (
    scores$Neutrophil_detected >= 2 &
      scores$Neutrophil_z_sample >= 1.5
  )

  scores$hsc_strong_refined <- (
    scores$HSC_Fibroblast_detected >= 2 &
      scores$HSC_Fibroblast_z_sample >= 1.5
  )

  scores$tcell_strong_refined <- (
    scores$T_cell_detected >= 2 &
      scores$T_cell_z_sample >= 1.5
  )

  scores$likely_nonMphi_refined <- (
    (
      scores$endothelial_strong_refined |
        scores$neutrophil_strong_refined |
        scores$hsc_strong_refined |
        scores$tcell_strong_refined
    ) &
      scores$Myeloid_z_sample < -0.5
  )

  # Cluster summary
  rr <- lapply(sort(unique(scores$cluster)), function(cl) {
    ix <- which(scores$cluster == cl)
    data.frame(
      species = species,
      cluster = cl,
      n_cells = length(ix),
      ambient_refined_fraction = mean(scores$ambient_elevated_refined[ix]),
      passenger_refined_fraction = mean(scores$passenger_RNA_compatible_refined[ix]),
      hep_lineage_refined_fraction = mean(scores$hep_lineage_elevated_refined[ix]),
      hep_doublet_refined_fraction = mean(scores$possible_hepatocyte_doublet_refined[ix]),
      endothelial_refined_fraction = mean(scores$endothelial_strong_refined[ix]),
      neutrophil_refined_fraction = mean(scores$neutrophil_strong_refined[ix]),
      hsc_refined_fraction = mean(scores$hsc_strong_refined[ix]),
      tcell_refined_fraction = mean(scores$tcell_strong_refined[ix]),
      likely_nonMphi_refined_fraction = mean(scores$likely_nonMphi_refined[ix]),
      mean_Myeloid_z_sample = mean(scores$Myeloid_z_sample[ix]),
      mean_Phago_z_sample = mean(scores$Phago_z_sample[ix]),
      stringsAsFactors = FALSE
    )
  })

  cluster_summary <- do.call(rbind, rr)

  # Gene-level focus report
  gene_rows <- list()
  k <- 0L

  focus_plus_rest <- unique(c(as.character(focus_clusters), "REST"))

  for (grp in focus_plus_rest) {

    ix <- if (grp == "REST") {
      which(!scores$cluster %in% as.character(focus_clusters))
    } else {
      which(scores$cluster == grp)
    }

    if (!length(ix)) next

    for (g in genes) {
      k <- k + 1L
      vals <- logcp10k[g, ix]
      raw <- as.numeric(sub[g, ix])

      module <- names(Filter(function(v) g %in% v, panels))
      module <- if (length(module)) module[1] else "Other"

      gene_rows[[k]] <- data.frame(
        species = species,
        group = grp,
        module = module,
        gene = g,
        n_cells = length(ix),
        pct_detected = mean(raw > 0) * 100,
        mean_logCP10k = mean(vals),
        median_logCP10k = median(vals),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    scores = scores,
    cluster_summary = cluster_summary,
    gene_focus = do.call(rbind, gene_rows)
  )
}

ha <- score_and_gene_stats(
  human,
  human_sample,
  human@meta.data[[HCL]],
  PANELS_H,
  "Human",
  c("3", "13")
)

ma <- score_and_gene_stats(
  mouse,
  mouse_sample,
  mouse@meta.data[[MCL]],
  PANELS_M,
  "Mouse",
  c("8", "15", "23", "26")
)

write_tsv(
  ha$cluster_summary,
  file.path(TABDIR, "01_Human_refined_cluster_QC_v1.2.1.tsv")
)

write_tsv(
  ma$cluster_summary,
  file.path(TABDIR, "02_Mouse_refined_cluster_QC_v1.2.1.tsv")
)

write_tsv(
  ha$gene_focus,
  file.path(TABDIR, "03_Human_focus_gene_level_QC_v1.2.1.tsv")
)

write_tsv(
  ma$gene_focus,
  file.path(TABDIR, "04_Mouse_focus_gene_level_QC_v1.2.1.tsv")
)

write_tsv(
  ha$scores[, c(
    "cell", "sample", "cluster",
    "ambient_elevated_refined",
    "passenger_RNA_compatible_refined",
    "hep_lineage_elevated_refined",
    "possible_hepatocyte_doublet_refined",
    "endothelial_strong_refined",
    "neutrophil_strong_refined",
    "hsc_strong_refined",
    "tcell_strong_refined",
    "likely_nonMphi_refined",
    "Myeloid_z_sample",
    "Phago_z_sample"
  )],
  file.path(TABDIR, "05_Human_refined_cell_flags_v1.2.1.tsv")
)

write_tsv(
  ma$scores[, c(
    "cell", "sample", "cluster",
    "ambient_elevated_refined",
    "passenger_RNA_compatible_refined",
    "hep_lineage_elevated_refined",
    "possible_hepatocyte_doublet_refined",
    "endothelial_strong_refined",
    "neutrophil_strong_refined",
    "hsc_strong_refined",
    "tcell_strong_refined",
    "likely_nonMphi_refined",
    "Myeloid_z_sample",
    "Phago_z_sample"
  )],
  file.path(TABDIR, "06_Mouse_refined_cell_flags_v1.2.1.tsv")
)

msg("STATUS: SUCCESS")
msg("Review:")
msg("  01_Human_refined_cluster_QC_v1.2.1.tsv")
msg("  02_Mouse_refined_cluster_QC_v1.2.1.tsv")
msg("  03_Human_focus_gene_level_QC_v1.2.1.tsv")
msg("  04_Mouse_focus_gene_level_QC_v1.2.1.tsv")
