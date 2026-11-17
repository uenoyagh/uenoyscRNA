#!/usr/bin/env Rscript

# ==============================================================================
# WholeLiver_DecontX_Mphi_Sensitivity_v1.0.1.R
#
# PURPOSE
# -------
# Whole-liver, sample-wise DecontX sensitivity analysis for Human and Mouse MΦ.
#
# Principles:
# - Source RDS files are NEVER overwritten.
# - DecontX is run sample-by-sample on whole-liver filtered counts.
# - Broad lineage labels are supplied to DecontX to avoid treating MΦ
#   functional states (M1/M2/LAM/etc.) as separate contamination sources.
# - Only decontaminated counts for the already-frozen MΦ barcodes are retained.
# - Original vs corrected MΦ programs are compared quantitatively.
# - This is a SENSITIVITY ANALYSIS because raw/empty-droplet background is absent.
#
# Human:
#   Whole liver:  Human_MASH_LT1to8_Canonical_v1.3.3.rds
#   Sample:       orig.ident
#   Cell type:    celltype_final_v133
#   Frozen MΦ:    res0_8_cluster_frozen_v18
#
# Mouse:
#   Whole liver:  Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds
#   Sample:       sample_display_FIXED2
#   Cell type:    celltype_for_R8plot_FIXED2
#   Frozen MΦ:    cluster_res2
#
# Focus clusters:
#   Human: 3, 13
#   Mouse: 8, 15, 23, 26
#
# Version: v1.0.1
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)
set.seed(20261118)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(decontX)
})

VERSION <- "v1.0.1"

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_WHOLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_Canonical_v1.3.3",
  "objects",
  "Human_MASH_LT1to8_Canonical_v1.3.3.rds"
)

HUMAN_MPHI <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

MOUSE_WHOLE <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

MOUSE_MPHI <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_Mphi",
  "WholeLiver_DecontX_Mphi_Sensitivity_v1.0.1"
)

TABDIR <- file.path(OUTDIR, "tables")
MATDIR <- file.path(OUTDIR, "corrected_Mphi_matrices")
LOGDIR <- file.path(OUTDIR, "logs")
PROVDIR <- file.path(OUTDIR, "provenance")

for (d in c(OUTDIR, TABDIR, MATDIR, LOGDIR, PROVDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOGFILE <- file.path(LOGDIR, "run.log")
cat("", file = LOGFILE)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOGFILE, append = TRUE)
}

write_tsv <- function(x, f) {
  write.table(
    x,
    f,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

for (f in c(HUMAN_WHOLE, HUMAN_MPHI, MOUSE_WHOLE, MOUSE_MPHI)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

# ------------------------------------------------------------------------------
# Broad lineage mapping
# ------------------------------------------------------------------------------

collapse_broad_lineage <- function(x) {
  x <- as.character(x)
  out <- rep("Other", length(x))

  out[grepl("Ambiguous", x, ignore.case = TRUE)] <- "Ambiguous"
  out[grepl("Cycl", x, ignore.case = TRUE)] <- "Cycling"

  out[grepl("Neutro", x, ignore.case = TRUE)] <- "Neutrophil"

  out[
    grepl("Hepatocyte|Hepatoblast", x, ignore.case = TRUE)
  ] <- "Hepatocyte"

  out[
    grepl("LSEC|Endothelial", x, ignore.case = TRUE)
  ] <- "Endothelial"

  out[
    grepl("HSC|Fibroblast|Mesenchymal|Perivascular", x, ignore.case = TRUE)
  ] <- "Mesenchymal"

  out[
    grepl("Cholangiocyte|Ductular", x, ignore.case = TRUE)
  ] <- "Cholangiocyte"

  out[
    grepl("Macrophage|Monocyte|Dendritic", x, ignore.case = TRUE)
  ] <- "Myeloid"

  out[
    grepl("Plasma|B_cell|Bcell", x, ignore.case = TRUE)
  ] <- "B_Plasma"

  out[
    grepl("T_cell|Tcell|NK", x, ignore.case = TRUE)
  ] <- "T_NK"

  out[
    grepl("RBC|Eryth", x, ignore.case = TRUE)
  ] <- "Erythroid"

  out[
    grepl("Platelet", x, ignore.case = TRUE)
  ] <- "Platelet"

  out[
    grepl("Mesothelial", x, ignore.case = TRUE)
  ] <- "Mesothelial"

  out
}

# ------------------------------------------------------------------------------
# Seurat v5 count-layer extraction
# ------------------------------------------------------------------------------

extract_counts_for_cells <- function(obj, cells, assay = "RNA") {

  cells <- as.character(cells)

  lyr <- Layers(obj[[assay]])
  count_layers <- grep("^counts($|\\.)", lyr, value = TRUE)

  if (!length(count_layers)) {
    stop("No counts layer found in assay ", assay)
  }

  pieces <- list()

  for (ly in count_layers) {
    m <- LayerData(obj, assay = assay, layer = ly)
    keep <- intersect(colnames(m), cells)

    if (length(keep)) {
      pieces[[ly]] <- m[, keep, drop = FALSE]
    }
  }

  if (!length(pieces)) {
    stop("No requested cells found in count layers.")
  }

  common_genes <- Reduce(
    intersect,
    lapply(pieces, rownames)
  )

  if (!length(common_genes)) {
    stop("No common genes across count layers.")
  }

  pieces <- lapply(
    pieces,
    function(m) m[common_genes, , drop = FALSE]
  )

  mat <- do.call(cbind, pieces)

  if (anyDuplicated(colnames(mat))) {
    stop("Duplicated cell barcodes after count-layer assembly.")
  }

  missing_cells <- setdiff(cells, colnames(mat))

  if (length(missing_cells)) {
    stop(
      "Missing ",
      length(missing_cells),
      " requested cells from assembled counts. First: ",
      paste(head(missing_cells, 10), collapse = ", ")
    )
  }

  mat <- mat[, cells, drop = FALSE]
  mat <- as(mat, "dgCMatrix")

  mat
}

# ------------------------------------------------------------------------------
# Functional panels
# ------------------------------------------------------------------------------

PAN_H <- list(
  Ambient_hepatocyte = c("ALB", "APOA1", "APOA2", "TTR"),
  Hepatocyte_lineage = c("ASGR1", "HNF4A", "CPS1", "ALDOB", "APOB", "CYP3A4"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "LYZ", "CSF1R", "CD68"),
  M1 = c("IL1B", "TNF", "CCL3", "CCL4", "CXCL2"),
  M2 = c("MRC1", "CD163", "HMOX1", "MAF"),
  Resident = c("MARCO", "VSIG4", "FOLR2", "TIMD4", "CD5L", "C1QA", "C1QB", "C1QC"),
  Recruited = c("CCR2", "S100A8", "S100A9", "VCAN", "FCN1", "SELL"),
  LAM = c("TREM2", "GPNMB", "CD9", "APOE", "LPL", "FABP5", "SPP1"),
  IFN = c("STAT1", "CXCL10", "ISG15", "IFIT1", "IFIT2", "IFIT3"),
  Phagocytosis = c(
    "MERTK", "AXL", "TIMD4", "MARCO", "STAB1", "MSR1", "CD36",
    "C1QA", "C1QB", "C1QC", "CTSB", "CTSD", "LAMP1"
  )
)

PAN_M <- list(
  Ambient_hepatocyte = c("Alb", "Apoa1", "Apoa2", "Ttr"),
  Hepatocyte_lineage = c("Asgr1", "Hnf4a", "Cps1", "Aldob", "Apob", "Cyp3a11"),
  Myeloid = c("Lst1", "Tyrobp", "Fcer1g", "Ctss", "Lyz2", "Csf1r", "Cd68"),
  M1 = c("Il1b", "Tnf", "Ccl3", "Ccl4", "Cxcl2"),
  M2 = c("Mrc1", "Cd163", "Hmox1", "Maf"),
  Resident = c("Marco", "Vsig4", "Folr2", "Timd4", "Cd5l", "C1qa", "C1qb", "C1qc"),
  Recruited = c("Ccr2", "S100a8", "S100a9", "Vcan", "Fcna", "Sell"),
  LAM = c("Trem2", "Gpnmb", "Cd9", "Apoe", "Lpl", "Fabp5", "Spp1"),
  IFN = c("Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"),
  Phagocytosis = c(
    "Mertk", "Axl", "Timd4", "Marco", "Stab1", "Msr1", "Cd36",
    "C1qa", "C1qb", "C1qc", "Ctsb", "Ctsd", "Lamp1"
  )
)

# ------------------------------------------------------------------------------
# Program scoring
# ------------------------------------------------------------------------------

score_modules <- function(counts, panels) {

  lib <- Matrix::colSums(counts)
  lib[!is.finite(lib) | lib <= 0] <- 1

  genes_needed <- unique(unlist(panels, use.names = FALSE))
  genes_present <- intersect(genes_needed, rownames(counts))

  x <- counts[genes_present, , drop = FALSE]
  x <- as.matrix(x)
  logcp10k <- log1p(t(t(x) / lib) * 10000)

  out <- data.frame(
    cell = colnames(counts),
    stringsAsFactors = FALSE
  )

  for (nm in names(panels)) {
    g <- intersect(panels[[nm]], rownames(logcp10k))

    if (!length(g)) {
      out[[nm]] <- NA_real_
    } else {
      out[[nm]] <- colMeans(
        logcp10k[g, , drop = FALSE],
        na.rm = TRUE
      )
    }
  }

  out
}

# ------------------------------------------------------------------------------
# Run one species
# ------------------------------------------------------------------------------

run_species <- function(
    species,
    whole,
    mphi,
    sample_col,
    celltype_col,
    cluster_col,
    panels,
    focus_clusters) {

  msg("================================================================")
  msg(species, " | starting")

  stopifnot(sample_col %in% colnames(whole@meta.data))
  stopifnot(celltype_col %in% colnames(whole@meta.data))
  stopifnot(cluster_col %in% colnames(mphi@meta.data))

  samples <- as.character(whole@meta.data[[sample_col]])
  names(samples) <- colnames(whole)

  original_ct <- as.character(whole@meta.data[[celltype_col]])
  broad <- collapse_broad_lineage(original_ct)
  names(broad) <- colnames(whole)

  cluster <- as.character(mphi@meta.data[[cluster_col]])
  names(cluster) <- colnames(mphi)

  # Audit broad lineage mapping.
  broad_audit <- as.data.frame(
    table(
      original_celltype = original_ct,
      broad_lineage = broad
    ),
    stringsAsFactors = FALSE
  )
  broad_audit <- broad_audit[broad_audit$Freq > 0, , drop = FALSE]
  broad_audit$species <- species

  write_tsv(
    broad_audit,
    file.path(
      TABDIR,
      paste0(species, "_broad_lineage_mapping_v1.0.1.tsv")
    )
  )

  sample_levels <- unique(samples)

  contamination_rows <- list()
  sample_summary_rows <- list()
  label_count_rows <- list()

  corrected_mphi_list <- list()
  original_mphi_list <- list()

  for (s in sample_levels) {

    msg(species, " | sample ", s, " | extracting counts")

    sample_cells <- names(samples)[samples == s]
    counts <- extract_counts_for_cells(
      whole,
      sample_cells,
      assay = "RNA"
    )

    z0 <- broad[colnames(counts)]
    z0[is.na(z0) | z0 == ""] <- "Other"

    # Merge very rare broad labels within each sample so DecontX does not
    # estimate contamination from tiny pseudo-clusters.
    ztab <- table(z0)
    rare <- names(ztab)[ztab < 20]

    z <- z0
    if (length(rare)) {
      z[z %in% rare] <- "Other_rare"
    }

    if (length(unique(z)) < 2) {
      stop(
        species, " / ", s,
        ": fewer than two broad lineages after mapping; cannot run background-free DecontX."
      )
    }

    labtab <- as.data.frame(table(broad_lineage = z), stringsAsFactors = FALSE)
    labtab$species <- species
    labtab$sample <- s
    label_count_rows[[length(label_count_rows) + 1L]] <- labtab

    msg(
      species, " | sample ", s,
      " | n=", ncol(counts),
      " | broad lineages=", length(unique(z)),
      " | running DecontX"
    )

    sce <- SingleCellExperiment(
      assays = list(counts = counts)
    )

    sce <- decontX(
      sce,
      z = as.integer(factor(z))
    )

    corrected <- decontXcounts(sce)
    corrected <- as(corrected, "dgCMatrix")

    contamination <- as.numeric(
      colData(sce)$decontX_contamination
    )

    cr <- data.frame(
      species = species,
      sample = s,
      cell = colnames(sce),
      broad_lineage = z,
      original_celltype = original_ct[colnames(sce)],
      decontX_contamination = contamination,
      stringsAsFactors = FALSE
    )

    contamination_rows[[length(contamination_rows) + 1L]] <- cr

    sample_summary_rows[[length(sample_summary_rows) + 1L]] <- data.frame(
      species = species,
      sample = s,
      n_cells = ncol(sce),
      n_broad_lineages = length(unique(z)),
      median_contamination = median(contamination, na.rm = TRUE),
      mean_contamination = mean(contamination, na.rm = TRUE),
      p90_contamination = as.numeric(
        quantile(contamination, 0.90, na.rm = TRUE)
      ),
      fraction_gt_0.20 = mean(contamination > 0.20, na.rm = TRUE),
      fraction_gt_0.50 = mean(contamination > 0.50, na.rm = TRUE),
      fraction_gt_0.70 = mean(contamination > 0.70, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    mphi_cells <- intersect(
      colnames(mphi),
      colnames(sce)
    )

    if (length(mphi_cells)) {

      corrected_mphi <- corrected[, mphi_cells, drop = FALSE]
      original_mphi <- counts[, mphi_cells, drop = FALSE]

      corrected_mphi_list[[s]] <- corrected_mphi
      original_mphi_list[[s]] <- original_mphi

      saveRDS(
        corrected_mphi,
        file.path(
          MATDIR,
          paste0(
            species,
            "_",
            gsub("[^A-Za-z0-9_.-]", "_", s),
            "_Mphi_decontX_corrected_counts_v1.0.1.rds"
          )
        ),
        compress = "xz"
      )
    }

    rm(sce, corrected, counts)
    gc(verbose = FALSE)
  }

  contamination_df <- do.call(
    rbind,
    contamination_rows
  )

  sample_summary_df <- do.call(
    rbind,
    sample_summary_rows
  )

  label_count_df <- do.call(
    rbind,
    label_count_rows
  )

  write_tsv(
    contamination_df,
    file.path(
      TABDIR,
      paste0(species, "_whole_liver_decontX_per_cell_v1.0.1.tsv")
    )
  )

  write_tsv(
    sample_summary_df,
    file.path(
      TABDIR,
      paste0(species, "_whole_liver_decontX_sample_summary_v1.0.1.tsv")
    )
  )

  write_tsv(
    label_count_df,
    file.path(
      TABDIR,
      paste0(species, "_DecontX_broad_lineage_counts_by_sample_v1.0.1.tsv")
    )
  )

  # Assemble only the MΦ matrices.
  common_genes <- Reduce(
    intersect,
    lapply(corrected_mphi_list, rownames)
  )

  corrected_mphi_list <- lapply(
    corrected_mphi_list,
    function(m) m[common_genes, , drop = FALSE]
  )

  original_mphi_list <- lapply(
    original_mphi_list,
    function(m) m[common_genes, , drop = FALSE]
  )

  corr_mphi <- do.call(
    cbind,
    corrected_mphi_list
  )

  orig_mphi <- do.call(
    cbind,
    original_mphi_list
  )

  if (anyDuplicated(colnames(corr_mphi))) {
    stop(species, ": duplicated MΦ barcodes after assembly.")
  }

  if (!all(colnames(mphi) %in% colnames(corr_mphi))) {
    missing <- setdiff(
      colnames(mphi),
      colnames(corr_mphi)
    )
    stop(
      species,
      ": corrected MΦ matrix missing ",
      length(missing),
      " frozen MΦ cells."
    )
  }

  corr_mphi <- corr_mphi[, colnames(mphi), drop = FALSE]
  orig_mphi <- orig_mphi[, colnames(mphi), drop = FALSE]

  saveRDS(
    corr_mphi,
    file.path(
      MATDIR,
      paste0(
        species,
        "_Mphi_decontX_corrected_counts_ALL_v1.0.1.rds"
      )
    ),
    compress = "xz"
  )

  # --------------------------------------------------------------------------
  # Original vs corrected program scores
  # --------------------------------------------------------------------------

  orig_score <- score_modules(
    orig_mphi,
    panels
  )

  corr_score <- score_modules(
    corr_mphi,
    panels
  )

  stopifnot(
    identical(
      orig_score$cell,
      corr_score$cell
    )
  )

  orig_score$cluster <- cluster[orig_score$cell]
  corr_score$cluster <- cluster[corr_score$cell]

  # Overall module concordance
  concord <- lapply(
    names(panels),
    function(mod) {
      data.frame(
        species = species,
        module = mod,
        spearman_original_vs_corrected = safe_cor(
          orig_score[[mod]],
          corr_score[[mod]],
          "spearman"
        ),
        mean_original = mean(
          orig_score[[mod]],
          na.rm = TRUE
        ),
        mean_corrected = mean(
          corr_score[[mod]],
          na.rm = TRUE
        ),
        mean_delta_corrected_minus_original = mean(
          corr_score[[mod]] - orig_score[[mod]],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  concord <- do.call(
    rbind,
    concord
  )

  write_tsv(
    concord,
    file.path(
      TABDIR,
      paste0(
        species,
        "_Mphi_module_original_vs_DecontX_v1.0.1.tsv"
      )
    )
  )

  # Cluster-level program comparison
  clusters <- sort(unique(cluster))
  cluster_rows <- list()

  for (cl in clusters) {
    ix <- which(orig_score$cluster == cl)

    for (mod in names(panels)) {
      cluster_rows[[length(cluster_rows) + 1L]] <- data.frame(
        species = species,
        cluster = cl,
        n_cells = length(ix),
        module = mod,
        mean_original = mean(
          orig_score[[mod]][ix],
          na.rm = TRUE
        ),
        mean_corrected = mean(
          corr_score[[mod]][ix],
          na.rm = TRUE
        ),
        delta_corrected_minus_original = mean(
          corr_score[[mod]][ix] -
            orig_score[[mod]][ix],
          na.rm = TRUE
        ),
        spearman_original_vs_corrected = safe_cor(
          orig_score[[mod]][ix],
          corr_score[[mod]][ix],
          "spearman"
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  cluster_compare <- do.call(
    rbind,
    cluster_rows
  )

  write_tsv(
    cluster_compare,
    file.path(
      TABDIR,
      paste0(
        species,
        "_Mphi_cluster_module_original_vs_DecontX_v1.0.1.tsv"
      )
    )
  )

  focus_compare <- cluster_compare[
    cluster_compare$cluster %in% as.character(focus_clusters),
    ,
    drop = FALSE
  ]

  write_tsv(
    focus_compare,
    file.path(
      TABDIR,
      paste0(
        species,
        "_FOCUS_clusters_original_vs_DecontX_v1.0.1.tsv"
      )
    )
  )

  # --------------------------------------------------------------------------
  # DecontX contamination estimates inside frozen MΦ clusters
  # --------------------------------------------------------------------------

  mphi_contam <- contamination_df[
    contamination_df$cell %in% colnames(mphi),
    ,
    drop = FALSE
  ]

  mphi_contam$cluster <- cluster[mphi_contam$cell]

  cc <- lapply(
    sort(unique(mphi_contam$cluster)),
    function(cl) {
      ix <- which(mphi_contam$cluster == cl)

      data.frame(
        species = species,
        cluster = cl,
        n_cells = length(ix),
        median_decontX_contamination = median(
          mphi_contam$decontX_contamination[ix],
          na.rm = TRUE
        ),
        mean_decontX_contamination = mean(
          mphi_contam$decontX_contamination[ix],
          na.rm = TRUE
        ),
        p90_decontX_contamination = as.numeric(
          quantile(
            mphi_contam$decontX_contamination[ix],
            0.90,
            na.rm = TRUE
          )
        ),
        fraction_gt_0.20 = mean(
          mphi_contam$decontX_contamination[ix] > 0.20,
          na.rm = TRUE
        ),
        fraction_gt_0.50 = mean(
          mphi_contam$decontX_contamination[ix] > 0.50,
          na.rm = TRUE
        ),
        fraction_gt_0.70 = mean(
          mphi_contam$decontX_contamination[ix] > 0.70,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  cc <- do.call(
    rbind,
    cc
  )

  write_tsv(
    cc,
    file.path(
      TABDIR,
      paste0(
        species,
        "_Mphi_cluster_decontX_contamination_summary_v1.0.1.tsv"
      )
    )
  )

  msg(species, " | completed")

  invisible(
    list(
      sample_summary = sample_summary_df,
      concordance = concord,
      cluster_compare = cluster_compare,
      contamination_cluster = cc
    )
  )
}

# ------------------------------------------------------------------------------
# Load objects
# ------------------------------------------------------------------------------

msg("Loading Human whole-liver...")
human_whole <- readRDS(HUMAN_WHOLE)

msg("Loading Human frozen MΦ...")
human_mphi <- readRDS(HUMAN_MPHI)

msg("Loading Mouse whole-liver...")
mouse_whole <- readRDS(MOUSE_WHOLE)

msg("Loading Mouse frozen MΦ...")
mouse_mphi <- readRDS(MOUSE_MPHI)

# Hard barcode-overlap check.
if (length(intersect(colnames(human_whole), colnames(human_mphi))) != ncol(human_mphi)) {
  stop("Human frozen MΦ barcodes are not fully present in whole-liver object.")
}

if (length(intersect(colnames(mouse_whole), colnames(mouse_mphi))) != ncol(mouse_mphi)) {
  stop("Mouse frozen MΦ barcodes are not fully present in whole-liver object.")
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

human_result <- run_species(
  species = "Human",
  whole = human_whole,
  mphi = human_mphi,
  sample_col = "orig.ident",
  celltype_col = "celltype_final_v133",
  cluster_col = "res0_8_cluster_frozen_v18",
  panels = PAN_H,
  focus_clusters = c("3", "13")
)

rm(human_whole)
gc(verbose = FALSE)

mouse_result <- run_species(
  species = "Mouse",
  whole = mouse_whole,
  mphi = mouse_mphi,
  sample_col = "sample_display_FIXED2",
  celltype_col = "celltype_for_R8plot_FIXED2",
  cluster_col = "cluster_res2",
  panels = PAN_M,
  focus_clusters = c("8", "15", "23", "26")
)

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

settings <- data.frame(
  item = c(
    "version",
    "human_whole",
    "human_mphi",
    "human_sample_col",
    "human_celltype_col",
    "human_cluster_col",
    "mouse_whole",
    "mouse_mphi",
    "mouse_sample_col",
    "mouse_celltype_col",
    "mouse_cluster_col",
    "decontX_mode",
    "background_matrix",
    "source_object_policy"
  ),
  value = c(
    VERSION,
    HUMAN_WHOLE,
    HUMAN_MPHI,
    "orig.ident",
    "celltype_final_v133",
    "res0_8_cluster_frozen_v18",
    MOUSE_WHOLE,
    MOUSE_MPHI,
    "sample_display_FIXED2",
    "celltype_for_R8plot_FIXED2",
    "cluster_res2",
    "whole-liver sample-wise with supplied broad lineage z; decontX 1.10.0 release API (no embedding argument)",
    "NOT AVAILABLE; sensitivity analysis only",
    "Source RDS files are never overwritten"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  settings,
  file.path(
    PROVDIR,
    "analysis_settings_v1.0.1.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    PROVDIR,
    "sessionInfo.txt"
  ),
  useBytes = TRUE
)

msg("================================================================")
msg("STATUS: SUCCESS")
msg("Primary review files:")
msg("  Human_whole_liver_decontX_sample_summary_v1.0.1.tsv")
msg("  Mouse_whole_liver_decontX_sample_summary_v1.0.1.tsv")
msg("  Human_Mphi_cluster_decontX_contamination_summary_v1.0.1.tsv")
msg("  Mouse_Mphi_cluster_decontX_contamination_summary_v1.0.1.tsv")
msg("  Human_FOCUS_clusters_original_vs_DecontX_v1.0.1.tsv")
msg("  Mouse_FOCUS_clusters_original_vs_DecontX_v1.0.1.tsv")
msg("  Human_Mphi_module_original_vs_DecontX_v1.0.1.tsv")
msg("  Mouse_Mphi_module_original_vs_DecontX_v1.0.1.tsv")
msg("================================================================")
