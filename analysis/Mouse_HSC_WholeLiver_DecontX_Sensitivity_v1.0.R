#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_WholeLiver_DecontX_Sensitivity_v1.0.R
#
# PURPOSE
# -------
# Sample-wise DecontX sensitivity analysis for the refined Mouse HSC object
# using the matched whole-liver filtered-count object as the broad cellular
# background.
#
# This follows the same policy as the previously successful Mouse MΦ DecontX
# audit:
#   - original filtered counts remain PRIMARY
#   - whole-liver filtered-count DecontX is SENSITIVITY ANALYSIS ONLY
#   - no raw/empty-droplet background is supplied
#   - source RDS files are NEVER overwritten
#   - no reclustering is performed
#
# Target HSC working solution:
#   res0.8, hsc_cluster_working_v610
#
# Focus clusters:
#   c14 = B-cell-like mixed cluster candidate
#   c15 = LSEC/endothelial mixed cluster candidate
#   c16 = striated/skeletal-muscle-like non-HSC candidate
#
# DecontX API:
#   decontX(sce, z = as.integer(factor(z)))
# Compatible with the installed decontX/celda workflow already validated
# in this project. Do NOT use the unsupported embedding= argument.
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(celda)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2",
  "RDS",
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_HSC",
  "Mouse_HSC_WholeLiver_DecontX_Sensitivity_v1.0"
)

TABDIR <- file.path(OUTDIR, "tables")
MATDIR <- file.path(OUTDIR, "corrected_HSC_matrices")
PROVDIR <- file.path(OUTDIR, "provenance")

for (d in c(OUTDIR, TABDIR, MATDIR, PROVDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
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

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

pct_change <- function(before, after) {
  if (!is.finite(before) || before == 0) return(NA_real_)
  100 * (after - before) / abs(before)
}

# ------------------------------------------------------------------------------
# Load
# ------------------------------------------------------------------------------

if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)
if (!file.exists(WHOLE_RDS)) stop("Missing whole-liver RDS: ", WHOLE_RDS)

message("Loading HSC object...")
hsc <- readRDS(HSC_RDS)

message("Loading whole-liver object...")
whole <- readRDS(WHOLE_RDS)

hmd <- hsc@meta.data
wmd <- whole@meta.data

required_h <- c(
  "sample_hsc_v610",
  "hsc_cluster_working_v610",
  "HSC_state3_v612",
  "HSC_state2_v612"
)

missing_h <- setdiff(required_h, colnames(hmd))
if (length(missing_h)) {
  stop("Missing HSC metadata: ", paste(missing_h, collapse = ", "))
}

required_w <- c(
  "sample_display_FIXED2",
  "celltype_for_R8plot_FIXED2"
)

missing_w <- setdiff(required_w, colnames(wmd))
if (length(missing_w)) {
  stop("Missing whole-liver metadata: ", paste(missing_w, collapse = ", "))
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
FOCUS_CLUSTERS <- c("14", "15", "16")

# Strict barcode/sample integrity check
overlap <- intersect(colnames(hsc), colnames(whole))

if (length(overlap) != ncol(hsc)) {
  stop(
    "Not all HSC barcodes are present in whole-liver object: ",
    length(overlap), " / ", ncol(hsc)
  )
}

hidx <- match(colnames(hsc), colnames(hsc))
widx <- match(colnames(hsc), colnames(whole))

hs <- as.character(hmd$sample_hsc_v610[hidx])
ws <- as.character(wmd$sample_display_FIXED2[widx])

if (!all(hs == ws)) {
  stop("HSC vs whole-liver sample labels do not match for all HSC barcodes.")
}

# ------------------------------------------------------------------------------
# Helper for Seurat v5 counts
# ------------------------------------------------------------------------------

get_counts_joined <- function(obj) {
  assay <- obj[["RNA"]]

  if (inherits(assay, "Assay5")) {
    layers <- Layers(assay)
    count_layers <- grep("^counts($|\\.)", layers, value = TRUE)

    if (!length(count_layers)) {
      stop("No RNA counts layer found.")
    }

    # If split layers exist, join only within the already subsetted object.
    if (length(count_layers) > 1 || !("counts" %in% count_layers)) {
      obj <- JoinLayers(obj, assay = "RNA")
    }

    mat <- LayerData(obj[["RNA"]], layer = "counts")
  } else {
    mat <- GetAssayData(obj, assay = "RNA", slot = "counts")
  }

  if (!inherits(mat, "Matrix")) {
    mat <- Matrix(mat, sparse = TRUE)
  }

  mat
}

# ------------------------------------------------------------------------------
# Panels
# ------------------------------------------------------------------------------

PANELS <- list(
  HSC_core = c(
    "Lrat", "Rbp1", "Reln", "Cygb", "Des", "Pdgfra"
  ),

  qHSC = c(
    "Lrat", "Rbp1", "Reln", "Cygb", "Colec11", "Vipr1"
  ),

  ECM_activated_HSC = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1",
    "Fn1", "Lum", "Dcn", "Bgn", "Timp1"
  ),

  Contractile_HSC = c(
    "Acta2", "Tagln", "Myl9", "Tpm2", "Myh11"
  ),

  Pericyte_VSMC = c(
    "Rgs5", "Cspg4", "Pdgfrb", "Des",
    "Acta2", "Tagln", "Myh11", "Notch3"
  ),

  Hepatocyte_ambient = c(
    "Alb", "Apoa1", "Apoa2", "Ttr"
  ),

  Hepatocyte_strict = c(
    "Asgr1", "Hnf4a", "Cps1", "Cyp3a11", "Aldob", "Apob"
  ),

  Endothelial_LSEC = c(
    "Pecam1", "Vwf", "Kdr", "Emcn", "Ramp2",
    "Cldn5", "Egfl7", "Plvap", "Clec4g", "Oit3"
  ),

  B_cell = c(
    "Ms4a1", "Cd79a", "Cd79b", "Cd19", "Pax5", "Ighm", "Igkc", "Ebf1"
  ),

  Myeloid = c(
    "Lst1", "Tyrobp", "Fcer1g", "Ctss", "Lyz2", "Csf1r", "Cd68"
  ),

  Neutrophil = c(
    "Csf3r", "Mpo", "Elane", "Camp", "Ltf", "Ngp", "Retnlg"
  ),

  Striated_muscle = c(
    "Tnnt3", "Ttn", "Actn2", "Neb", "Tnni2", "Ryr1",
    "Stac3", "Pygm", "Obscn", "Ldb3", "Tcap", "Atp2a1"
  )
)

ALL_DIRECT_GENES <- unique(unlist(PANELS))

# ------------------------------------------------------------------------------
# Score helper using log1p(CP10k)
# ------------------------------------------------------------------------------

normalize_logcp10k <- function(counts) {
  lib <- Matrix::colSums(counts)

  if (any(lib <= 0)) {
    stop("Zero library size encountered.")
  }

  cp10k <- t(t(counts) / lib * 1e4)
  log1p(cp10k)
}

score_panel <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[g, , drop = FALSE]))
}

# ------------------------------------------------------------------------------
# Run sample-wise DecontX
# ------------------------------------------------------------------------------

contam_rows <- list()
sample_rows <- list()
cluster_panel_rows <- list()
cluster_gene_rows <- list()
corrected_list <- list()

kc <- 0L
ks <- 0L
kp <- 0L
kg <- 0L

set.seed(20261118)

for (s in SAMPLES) {
  message("============================================================")
  message("Sample: ", s)

  # Whole-liver cells for this sample
  whole_cells <- rownames(wmd)[
    as.character(wmd$sample_display_FIXED2) == s
  ]

  if (!length(whole_cells)) {
    stop("No whole-liver cells for sample: ", s)
  }

  whole_s <- subset(whole, cells = whole_cells)

  counts_whole <- get_counts_joined(whole_s)

  # Ensure broad-lineage labels are in the same order as counts columns
  z <- as.character(
    whole_s@meta.data[
      colnames(counts_whole),
      "celltype_for_R8plot_FIXED2"
    ]
  )

  z[is.na(z) | z == ""] <- "Unknown"

  message("Whole-liver cells: ", ncol(counts_whole))
  message("Broad-lineage classes: ", length(unique(z)))

  sce <- SingleCellExperiment(
    assays = list(counts = counts_whole)
  )

  message("Running DecontX...")
  sce <- decontX(
    sce,
    z = as.integer(factor(z))
  )

  corr_whole <- decontXcounts(sce)
  contamination <- as.numeric(colData(sce)$decontX_contamination)

  names(contamination) <- colnames(corr_whole)

  # HSC barcodes for this sample
  hsc_cells <- rownames(hmd)[
    as.character(hmd$sample_hsc_v610) == s
  ]

  missing_hsc <- setdiff(hsc_cells, colnames(corr_whole))

  if (length(missing_hsc)) {
    stop(
      "HSC barcodes missing from DecontX output for ", s, ": ",
      length(missing_hsc)
    )
  }

  # Original HSC counts from source HSC object
  # Subset first so JoinLayers remains lightweight.
  hsc_s <- subset(hsc, cells = hsc_cells)
  orig_hsc <- get_counts_joined(hsc_s)

  # Corrected HSC counts from whole-liver DecontX result
  corr_hsc <- corr_whole[, hsc_cells, drop = FALSE]

  # Gene ordering harmonization
  common_genes <- intersect(rownames(orig_hsc), rownames(corr_hsc))

  orig_hsc <- orig_hsc[common_genes, hsc_cells, drop = FALSE]
  corr_hsc <- corr_hsc[common_genes, hsc_cells, drop = FALSE]

  corrected_list[[s]] <- corr_hsc

  # Per-cell contamination
  cvec <- contamination[hsc_cells]

  for (cell in hsc_cells) {
    kc <- kc + 1L

    contam_rows[[kc]] <- data.frame(
      cell = cell,
      sample = s,
      cluster = as.character(hmd[cell, "hsc_cluster_working_v610"]),
      HSC_state3 = as.character(hmd[cell, "HSC_state3_v612"]),
      HSC_state2 = as.character(hmd[cell, "HSC_state2_v612"]),
      decontX_contamination = cvec[cell],
      stringsAsFactors = FALSE
    )
  }

  ks <- ks + 1L
  sample_rows[[ks]] <- data.frame(
    sample = s,
    n_whole_liver_cells = ncol(counts_whole),
    n_HSC_cells = length(hsc_cells),
    median_HSC_contamination = safe_median(cvec),
    mean_HSC_contamination = safe_mean(cvec),
    stringsAsFactors = FALSE
  )

  orig_log <- normalize_logcp10k(orig_hsc)
  corr_log <- normalize_logcp10k(corr_hsc)

  cluster_s <- as.character(
    hmd[hsc_cells, "hsc_cluster_working_v610"]
  )

  # Cluster × panel pre/post
  for (cl in sort(unique(cluster_s))) {
    ix <- which(cluster_s == cl)

    for (pn in names(PANELS)) {
      before <- safe_mean(score_panel(orig_log[, ix, drop = FALSE], PANELS[[pn]]))
      after  <- safe_mean(score_panel(corr_log[, ix, drop = FALSE], PANELS[[pn]]))

      kp <- kp + 1L
      cluster_panel_rows[[kp]] <- data.frame(
        sample = s,
        cluster = cl,
        n_cells = length(ix),
        panel = pn,
        original_mean_score = before,
        decontX_mean_score = after,
        absolute_change = after - before,
        percent_change = pct_change(before, after),
        stringsAsFactors = FALSE
      )
    }
  }

  # Direct gene pre/post for focus clusters 14/15/16 and ALL_HSC
  direct_genes <- intersect(ALL_DIRECT_GENES, common_genes)

  for (scope in c("ALL_HSC", FOCUS_CLUSTERS)) {
    if (scope == "ALL_HSC") {
      ix <- seq_along(hsc_cells)
    } else {
      ix <- which(cluster_s == scope)
      if (!length(ix)) next
    }

    for (g in direct_genes) {
      orig_cnt <- as.numeric(orig_hsc[g, ix])
      corr_cnt <- as.numeric(corr_hsc[g, ix])

      orig_expr <- as.numeric(orig_log[g, ix])
      corr_expr <- as.numeric(corr_log[g, ix])

      kg <- kg + 1L
      cluster_gene_rows[[kg]] <- data.frame(
        sample = s,
        cluster_scope = scope,
        gene = g,
        n_cells = length(ix),
        original_pct_detected = 100 * mean(orig_cnt > 0),
        decontX_pct_detected = 100 * mean(corr_cnt > 0),
        original_mean_logCP10k = safe_mean(orig_expr),
        decontX_mean_logCP10k = safe_mean(corr_expr),
        percent_change_mean_logCP10k = pct_change(
          safe_mean(orig_expr),
          safe_mean(corr_expr)
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  rm(
    whole_s,
    counts_whole,
    sce,
    corr_whole,
    hsc_s,
    orig_hsc,
    corr_hsc,
    orig_log,
    corr_log
  )
  invisible(gc())
}

# ------------------------------------------------------------------------------
# Save corrected HSC matrices
# ------------------------------------------------------------------------------

saveRDS(
  corrected_list,
  file.path(
    MATDIR,
    "Mouse_HSC_decontX_corrected_counts_by_sample_v1.0.rds"
  )
)

# Combined corrected matrix in HSC source-object cell order
common_genes_all <- Reduce(
  intersect,
  lapply(corrected_list, rownames)
)

corrected_list2 <- lapply(
  corrected_list,
  function(m) m[common_genes_all, , drop = FALSE]
)

corr_all <- do.call(cbind, corrected_list2)
corr_all <- corr_all[, colnames(hsc), drop = FALSE]

saveRDS(
  corr_all,
  file.path(
    MATDIR,
    "Mouse_HSC_decontX_corrected_counts_ALL_v1.0.rds"
  )
)

# ------------------------------------------------------------------------------
# Write raw summaries
# ------------------------------------------------------------------------------

contam_df <- do.call(rbind, contam_rows)
sample_df <- do.call(rbind, sample_rows)
cluster_panel_df <- do.call(rbind, cluster_panel_rows)
cluster_gene_df <- do.call(rbind, cluster_gene_rows)

write_tsv(
  contam_df,
  file.path(TABDIR, "01_HSC_cell_contamination.tsv")
)

write_tsv(
  sample_df,
  file.path(TABDIR, "02_HSC_sample_contamination_summary.tsv")
)

write_tsv(
  cluster_panel_df,
  file.path(TABDIR, "03_cluster_panel_pre_post_by_sample.tsv")
)

write_tsv(
  cluster_gene_df,
  file.path(TABDIR, "04_focus_direct_gene_pre_post_by_sample.tsv")
)

# ------------------------------------------------------------------------------
# Aggregate cluster × panel across biological samples
# ------------------------------------------------------------------------------

agg_panel_rows <- list()
k <- 0L

keys <- unique(cluster_panel_df[, c("cluster", "panel")])

for (i in seq_len(nrow(keys))) {
  cl <- keys$cluster[i]
  pn <- keys$panel[i]

  d <- cluster_panel_df[
    cluster_panel_df$cluster == cl &
      cluster_panel_df$panel == pn,
    ,
    drop = FALSE
  ]

  k <- k + 1L
  agg_panel_rows[[k]] <- data.frame(
    cluster = cl,
    panel = pn,
    n_samples = nrow(d),
    total_cells = sum(d$n_cells),
    original_mean_of_sample_means = safe_mean(d$original_mean_score),
    decontX_mean_of_sample_means = safe_mean(d$decontX_mean_score),
    percent_change = pct_change(
      safe_mean(d$original_mean_score),
      safe_mean(d$decontX_mean_score)
    ),
    stringsAsFactors = FALSE
  )
}

agg_panel <- do.call(rbind, agg_panel_rows)

write_tsv(
  agg_panel,
  file.path(TABDIR, "05_cluster_panel_pre_post_aggregated.tsv")
)

# Focus-cluster compact table
focus_panel <- agg_panel[
  agg_panel$cluster %in% FOCUS_CLUSTERS,
  ,
  drop = FALSE
]

write_tsv(
  focus_panel,
  file.path(TABDIR, "06_FOCUS_cluster14_15_16_panel_summary.tsv")
)

# ------------------------------------------------------------------------------
# Aggregate direct genes across samples
# ------------------------------------------------------------------------------

agg_gene_rows <- list()
k <- 0L

keys <- unique(cluster_gene_df[, c("cluster_scope", "gene")])

for (i in seq_len(nrow(keys))) {
  sc <- keys$cluster_scope[i]
  g <- keys$gene[i]

  d <- cluster_gene_df[
    cluster_gene_df$cluster_scope == sc &
      cluster_gene_df$gene == g,
    ,
    drop = FALSE
  ]

  # Weight detection by cell count, but keep expression as sample-balanced mean.
  total_n <- sum(d$n_cells)

  orig_det <- sum(
    d$original_pct_detected / 100 * d$n_cells
  ) / total_n * 100

  corr_det <- sum(
    d$decontX_pct_detected / 100 * d$n_cells
  ) / total_n * 100

  orig_expr <- safe_mean(d$original_mean_logCP10k)
  corr_expr <- safe_mean(d$decontX_mean_logCP10k)

  k <- k + 1L
  agg_gene_rows[[k]] <- data.frame(
    cluster_scope = sc,
    gene = g,
    n_samples = nrow(d),
    total_cells = total_n,
    original_pct_detected = orig_det,
    decontX_pct_detected = corr_det,
    original_mean_logCP10k = orig_expr,
    decontX_mean_logCP10k = corr_expr,
    percent_change_mean_logCP10k = pct_change(orig_expr, corr_expr),
    stringsAsFactors = FALSE
  )
}

agg_gene <- do.call(rbind, agg_gene_rows)

write_tsv(
  agg_gene,
  file.path(TABDIR, "07_focus_direct_gene_pre_post_aggregated.tsv")
)

# ------------------------------------------------------------------------------
# Focus-cluster contamination summary
# ------------------------------------------------------------------------------

contam_cluster_rows <- list()
k <- 0L

for (cl in sort(unique(contam_df$cluster))) {
  d <- contam_df[contam_df$cluster == cl, , drop = FALSE]

  k <- k + 1L
  contam_cluster_rows[[k]] <- data.frame(
    cluster = cl,
    n_cells = nrow(d),
    n_samples = length(unique(d$sample)),
    median_contamination = safe_median(d$decontX_contamination),
    mean_contamination = safe_mean(d$decontX_contamination),
    stringsAsFactors = FALSE
  )
}

contam_cluster <- do.call(rbind, contam_cluster_rows)

write_tsv(
  contam_cluster,
  file.path(TABDIR, "08_cluster_contamination_summary.tsv")
)

# ------------------------------------------------------------------------------
# QC decision scaffold
# ------------------------------------------------------------------------------

decision <- data.frame(
  cluster = c("14", "15", "16"),
  pre_DecontX_concern = c(
    "Coherent B-cell program; HSC-B mixed/doublet candidate",
    "Coherent LSEC/endothelial program; HSC-LSEC mixed/doublet candidate",
    "Strong striated/skeletal-muscle program with marked loss of HSC identity"
  ),
  pre_DecontX_status = c(
    "PRIMARY_EXCLUSION_CANDIDATE",
    "PRIMARY_EXCLUSION_CANDIDATE_SENSITIVITY_KEEP_CANDIDATE",
    "PRIMARY_EXCLUSION_CANDIDATE"
  ),
  final_status = "PENDING_POST_DECONTX_REVIEW",
  stringsAsFactors = FALSE
)

write_tsv(
  decision,
  file.path(TABDIR, "09_focus_cluster_decision_scaffold.tsv")
)

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

settings <- data.frame(
  item = c(
    "version",
    "HSC_RDS",
    "whole_liver_RDS",
    "HSC_cells",
    "whole_liver_cells",
    "sample_column_HSC",
    "sample_column_whole_liver",
    "broad_lineage_z_column",
    "DecontX_call",
    "raw_empty_droplet_background",
    "primary_data_policy",
    "reclustering",
    "source_RDS_overwrite"
  ),
  value = c(
    "v1.0",
    HSC_RDS,
    WHOLE_RDS,
    ncol(hsc),
    ncol(whole),
    "sample_hsc_v610",
    "sample_display_FIXED2",
    "celltype_for_R8plot_FIXED2",
    "decontX(sce, z = as.integer(factor(z)))",
    "NOT AVAILABLE / NOT USED",
    "Original filtered counts remain primary; DecontX is sensitivity only",
    "NONE",
    "NEVER"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  settings,
  file.path(PROVDIR, "analysis_settings_v1.0.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(PROVDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Corrected matrices:")
message("  Mouse_HSC_decontX_corrected_counts_by_sample_v1.0.rds")
message("  Mouse_HSC_decontX_corrected_counts_ALL_v1.0.rds")
message("Review first:")
message("  02_HSC_sample_contamination_summary.tsv")
message("  06_FOCUS_cluster14_15_16_panel_summary.tsv")
message("  07_focus_direct_gene_pre_post_aggregated.tsv")
message("  08_cluster_contamination_summary.tsv")
message("  09_focus_cluster_decision_scaffold.tsv")
message("============================================================")
