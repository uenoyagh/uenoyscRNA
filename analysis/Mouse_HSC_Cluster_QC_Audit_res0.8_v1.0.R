#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_Cluster_QC_Audit_res0.8_v1.0.R
#
# PURPOSE
# -------
# Formal cluster-level QC audit of the refined Mouse HSC working solution
# (stored res0.8), with special attention to sample-specific clusters 14 and 16.
#
# This script does NOT recluster and does NOT change annotations.
#
# Questions:
#   1) Do all res0.8 clusters retain coherent HSC identity?
#   2) Are clusters 14 / 16 bona fide rare HSC states or foreign-lineage/doublet
#      contaminants?
#   3) Are they unusually high-complexity relative to the HSC population?
#   4) What 3-state / 2-state annotations do they carry?
#   5) Which direct marker genes define them?
#
# POLICY
# ------
# Sample-specificity alone is NOT an exclusion criterion.
# Candidate exclusion requires coherent foreign-lineage signal and/or
# mixed-lineage/high-complexity evidence with loss of HSC identity.
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2",
  "RDS",
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_HSC",
  "Mouse_HSC_Cluster_QC_Audit_res0.8_v1.0"
)

TABDIR <- file.path(OUTDIR, "tables")
PROVDIR <- file.path(OUTDIR, "provenance")

for (d in c(OUTDIR, TABDIR, PROVDIR)) {
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

if (!file.exists(HSC_RDS)) stop("HSC RDS not found: ", HSC_RDS)

message("Loading HSC object...")
obj <- readRDS(HSC_RDS)
md <- obj@meta.data

CLUSTER_COL <- "hsc_cluster_working_v610"
RES08_COL <- "hsc_rpca_res_0_8_v610"
SAMPLE_COL <- "sample_hsc_v610"
STATE3_COL <- "HSC_state3_v612"
STATE2_COL <- "HSC_state2_v612"

required <- c(CLUSTER_COL, RES08_COL, SAMPLE_COL, STATE3_COL, STATE2_COL)
missing_cols <- setdiff(required, colnames(md))

if (length(missing_cols)) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

cluster <- as.character(md[[CLUSTER_COL]])
sample <- as.character(md[[SAMPLE_COL]])
state3 <- as.character(md[[STATE3_COL]])
state2 <- as.character(md[[STATE2_COL]])

if (!identical(cluster, as.character(md[[RES08_COL]]))) {
  stop("Working cluster is not exactly identical to stored res0.8.")
}

# ------------------------------------------------------------------------------
# RNA data
# ------------------------------------------------------------------------------

counts <- LayerData(obj[["RNA"]], layer = "counts")
data <- LayerData(obj[["RNA"]], layer = "data")

if (!inherits(counts, "Matrix")) counts <- Matrix(counts, sparse = TRUE)
if (!inherits(data, "Matrix")) data <- Matrix(data, sparse = TRUE)

if (!identical(colnames(counts), colnames(obj))) {
  counts <- counts[, colnames(obj), drop = FALSE]
}
if (!identical(colnames(data), colnames(obj))) {
  data <- data[, colnames(obj), drop = FALSE]
}

# ------------------------------------------------------------------------------
# QC / identity marker panels
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

  Cycling = c(
    "Mki67", "Top2a", "Cenpf", "Ube2c", "Ccnb1", "Ccnb2"
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

  Myeloid = c(
    "Lst1", "Tyrobp", "Fcer1g", "Ctss", "Lyz2", "Csf1r", "Cd68"
  ),

  Cholangiocyte = c(
    "Krt19", "Krt7", "Epcam", "Sox9"
  ),

  T_cell = c(
    "Cd3d", "Cd3e", "Trac"
  ),

  NK = c(
    "Nkg7", "Klrd1", "Prf1"
  ),

  B_cell = c(
    "Ms4a1", "Cd79a", "Cd79b", "Cd19"
  ),

  Neutrophil = c(
    "Csf3r", "Mpo", "Elane", "Camp", "Ltf", "Ngp", "Retnlg"
  ),

  Erythroid = c(
    "Hbb-bs", "Hbb-bt", "Hba-a1", "Alas2", "Ahsp"
  ),

  Platelet = c(
    "Ppbp", "Pf4", "Gp9", "Itga2b", "Gp1ba"
  )
)

HSC_PANELS <- c(
  "HSC_core",
  "qHSC",
  "ECM_activated_HSC",
  "Contractile_HSC",
  "Pericyte_VSMC",
  "Cycling"
)

FOREIGN_PANELS <- setdiff(names(PANELS), HSC_PANELS)

coverage <- do.call(
  rbind,
  lapply(names(PANELS), function(pn) {
    g <- PANELS[[pn]]
    present <- intersect(g, rownames(data))

    data.frame(
      panel = pn,
      n_requested = length(g),
      n_present = length(present),
      coverage_fraction = length(present) / length(g),
      genes_present = paste(present, collapse = ";"),
      genes_missing = paste(setdiff(g, present), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  coverage,
  file.path(TABDIR, "01_marker_panel_coverage.tsv")
)

score_panel <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[g, , drop = FALSE]))
}

# Cell-level panel scores
score_df <- data.frame(
  cell = colnames(obj),
  cluster = cluster,
  sample = sample,
  HSC_state3 = state3,
  HSC_state2 = state2,
  stringsAsFactors = FALSE
)

for (pn in names(PANELS)) {
  score_df[[pn]] <- score_panel(data, PANELS[[pn]])
}

# ------------------------------------------------------------------------------
# Complexity metadata
# ------------------------------------------------------------------------------

detect_meta <- function(candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit)) hit[1] else NA_character_
}

NCOUNT_COL <- detect_meta(c(
  "nCount_RNA",
  "nCount_originalexp",
  "nCount"
))

NFEATURE_COL <- detect_meta(c(
  "nFeature_RNA",
  "nFeature_originalexp",
  "nFeature"
))

if (is.na(NCOUNT_COL)) {
  nCount_vec <- as.numeric(Matrix::colSums(counts))
} else {
  nCount_vec <- as.numeric(md[[NCOUNT_COL]])
}

if (is.na(NFEATURE_COL)) {
  nFeature_vec <- as.numeric(Matrix::colSums(counts > 0))
} else {
  nFeature_vec <- as.numeric(md[[NFEATURE_COL]])
}

score_df$nCount <- nCount_vec
score_df$nFeature <- nFeature_vec

# ------------------------------------------------------------------------------
# Cluster-level panel / complexity summary
# ------------------------------------------------------------------------------

clusters <- sort(unique(cluster))

cluster_summary_rows <- list()
k <- 0L

for (cl in clusters) {
  ix <- which(cluster == cl)
  stab <- table(sample[ix])

  k <- k + 1L
  row <- data.frame(
    cluster = cl,
    n_cells = length(ix),
    n_samples_present = sum(stab > 0),
    dominant_sample = names(stab)[which.max(stab)],
    dominant_sample_fraction = max(stab) / sum(stab),
    median_nCount = safe_median(nCount_vec[ix]),
    median_nFeature = safe_median(nFeature_vec[ix]),
    mean_nCount = safe_mean(nCount_vec[ix]),
    mean_nFeature = safe_mean(nFeature_vec[ix]),
    stringsAsFactors = FALSE
  )

  for (pn in names(PANELS)) {
    row[[paste0("mean_", pn)]] <- safe_mean(score_df[[pn]][ix])
  }

  cluster_summary_rows[[k]] <- row
}

cluster_summary <- do.call(rbind, cluster_summary_rows)

# Add ranks for HSC identity and foreign burden
cluster_summary$HSC_core_rank_high <- rank(
  -cluster_summary$mean_HSC_core,
  ties.method = "min"
)

foreign_mean_cols <- paste0("mean_", FOREIGN_PANELS)

cluster_summary$max_foreign_panel <- apply(
  cluster_summary[, foreign_mean_cols, drop = FALSE],
  1,
  function(x) {
    FOREIGN_PANELS[which.max(x)]
  }
)

cluster_summary$max_foreign_score <- apply(
  cluster_summary[, foreign_mean_cols, drop = FALSE],
  1,
  max,
  na.rm = TRUE
)

cluster_summary$foreign_rank_high <- rank(
  -cluster_summary$max_foreign_score,
  ties.method = "min"
)

write_tsv(
  cluster_summary,
  file.path(TABDIR, "02_cluster_QC_summary.tsv")
)

# ------------------------------------------------------------------------------
# State × cluster / sample × cluster
# ------------------------------------------------------------------------------

state3_tab <- as.data.frame(
  table(
    cluster = cluster,
    HSC_state3 = state3
  ),
  stringsAsFactors = FALSE
)
state3_tab <- state3_tab[state3_tab$Freq > 0, , drop = FALSE]

state3_tot <- aggregate(Freq ~ cluster, data = state3_tab, FUN = sum)
names(state3_tot)[2] <- "cluster_total"

state3_tab <- merge(
  state3_tab,
  state3_tot,
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)
state3_tab$fraction <- state3_tab$Freq / state3_tab$cluster_total

write_tsv(
  state3_tab,
  file.path(TABDIR, "03_cluster_x_HSC3state.tsv")
)

state2_tab <- as.data.frame(
  table(
    cluster = cluster,
    HSC_state2 = state2
  ),
  stringsAsFactors = FALSE
)
state2_tab <- state2_tab[state2_tab$Freq > 0, , drop = FALSE]

state2_tot <- aggregate(Freq ~ cluster, data = state2_tab, FUN = sum)
names(state2_tot)[2] <- "cluster_total"

state2_tab <- merge(
  state2_tab,
  state2_tot,
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)
state2_tab$fraction <- state2_tab$Freq / state2_tab$cluster_total

write_tsv(
  state2_tab,
  file.path(TABDIR, "04_cluster_x_HSC2state.tsv")
)

sample_tab <- as.data.frame(
  table(
    cluster = cluster,
    sample = sample
  ),
  stringsAsFactors = FALSE
)
sample_tab <- sample_tab[sample_tab$Freq > 0, , drop = FALSE]

write_tsv(
  sample_tab,
  file.path(TABDIR, "05_cluster_x_sample_counts.tsv")
)

# ------------------------------------------------------------------------------
# Direct gene detection / expression by cluster
# ------------------------------------------------------------------------------

ALL_GENES <- unique(unlist(PANELS))
ALL_GENES <- intersect(ALL_GENES, rownames(data))

gene_rows <- list()
k <- 0L

for (cl in clusters) {
  ix <- which(cluster == cl)

  for (g in ALL_GENES) {
    cnt <- as.numeric(counts[g, ix])
    expr <- as.numeric(data[g, ix])

    k <- k + 1L
    gene_rows[[k]] <- data.frame(
      cluster = cl,
      gene = g,
      n_cells = length(ix),
      pct_detected = 100 * mean(cnt > 0),
      mean_logexpr = safe_mean(expr),
      stringsAsFactors = FALSE
    )
  }
}

gene_cluster <- do.call(rbind, gene_rows)

write_tsv(
  gene_cluster,
  file.path(TABDIR, "06_direct_gene_expression_by_cluster.tsv")
)

# Focused direct-gene table for suspicious clusters + reference all-HSC mean
TARGET_CLUSTERS <- c("14", "16", "15")

target_genes <- unique(c(
  PANELS$HSC_core,
  PANELS$qHSC,
  PANELS$ECM_activated_HSC,
  PANELS$Contractile_HSC,
  PANELS$Pericyte_VSMC,
  PANELS$Hepatocyte_ambient,
  PANELS$Hepatocyte_strict,
  PANELS$Endothelial_LSEC,
  PANELS$Myeloid,
  PANELS$Cholangiocyte,
  PANELS$Neutrophil
))

target_genes <- intersect(target_genes, rownames(data))

target_rows <- list()
k <- 0L

for (g in target_genes) {
  # overall
  k <- k + 1L
  target_rows[[k]] <- data.frame(
    cluster = "ALL_HSC",
    gene = g,
    n_cells = ncol(obj),
    pct_detected = 100 * mean(as.numeric(counts[g, ]) > 0),
    mean_logexpr = safe_mean(as.numeric(data[g, ])),
    stringsAsFactors = FALSE
  )

  for (cl in TARGET_CLUSTERS[TARGET_CLUSTERS %in% clusters]) {
    ix <- which(cluster == cl)

    k <- k + 1L
    target_rows[[k]] <- data.frame(
      cluster = cl,
      gene = g,
      n_cells = length(ix),
      pct_detected = 100 * mean(as.numeric(counts[g, ix]) > 0),
      mean_logexpr = safe_mean(as.numeric(data[g, ix])),
      stringsAsFactors = FALSE
    )
  }
}

target_direct <- do.call(rbind, target_rows)

write_tsv(
  target_direct,
  file.path(TABDIR, "07_TARGET_cluster14_16_15_direct_markers.tsv")
)

# ------------------------------------------------------------------------------
# Top markers of each working cluster
# ------------------------------------------------------------------------------

gene_detected <- Matrix::rowSums(data > 0)
marker_genes <- rownames(data)[
  gene_detected >= max(20L, round(0.002 * ncol(obj)))
]

mat <- data[marker_genes, , drop = FALSE]
global_sum <- Matrix::rowSums(mat)
n_total <- ncol(mat)

top_rows <- list()
k <- 0L

for (cl in clusters) {
  ix <- which(cluster == cl)
  n_in <- length(ix)
  n_out <- n_total - n_in

  mean_in <- Matrix::rowMeans(mat[, ix, drop = FALSE])
  sum_in <- mean_in * n_in
  mean_out <- (global_sum - sum_in) / n_out

  delta <- as.numeric(mean_in - mean_out)
  names(delta) <- rownames(mat)

  ord <- order(delta, decreasing = TRUE)
  ord <- ord[seq_len(min(30L, length(ord)))]

  k <- k + 1L
  top_rows[[k]] <- data.frame(
    cluster = cl,
    rank = seq_along(ord),
    gene = names(delta)[ord],
    mean_logexpr_delta = delta[ord],
    stringsAsFactors = FALSE
  )
}

top_markers <- do.call(rbind, top_rows)

write_tsv(
  top_markers,
  file.path(TABDIR, "08_top30_markers_working_res0.8.tsv")
)

# ------------------------------------------------------------------------------
# RPCA centroid distance: identify nearest reference cluster for targets
# ------------------------------------------------------------------------------

REDUCTION <- "integrated.hsc.rpca"
emb <- Embeddings(obj, reduction = REDUCTION)
ndims <- min(30L, ncol(emb))
emb <- emb[, seq_len(ndims), drop = FALSE]

centroids <- do.call(
  rbind,
  lapply(clusters, function(cl) {
    colMeans(emb[cluster == cl, , drop = FALSE])
  })
)
rownames(centroids) <- clusters

dist_mat <- as.matrix(dist(centroids))

centroid_rows <- list()
k <- 0L

for (cl in TARGET_CLUSTERS[TARGET_CLUSTERS %in% clusters]) {
  d <- dist_mat[cl, ]
  d <- d[names(d) != cl]
  ord <- order(d)

  for (j in seq_len(min(5L, length(ord)))) {
    k <- k + 1L
    centroid_rows[[k]] <- data.frame(
      target_cluster = cl,
      neighbor_rank = j,
      neighbor_cluster = names(d)[ord[j]],
      centroid_distance = d[ord[j]],
      stringsAsFactors = FALSE
    )
  }
}

centroid_neighbors <- if (length(centroid_rows)) {
  do.call(rbind, centroid_rows)
} else {
  data.frame(
    target_cluster = character(0),
    neighbor_rank = integer(0),
    neighbor_cluster = character(0),
    centroid_distance = numeric(0)
  )
}

write_tsv(
  centroid_neighbors,
  file.path(TABDIR, "09_TARGET_cluster_centroid_neighbors.tsv")
)

# ------------------------------------------------------------------------------
# Compact target summary
# ------------------------------------------------------------------------------

target_summary <- cluster_summary[
  cluster_summary$cluster %in% TARGET_CLUSTERS,
  ,
  drop = FALSE
]

write_tsv(
  target_summary,
  file.path(TABDIR, "10_TARGET_cluster14_16_15_summary.tsv")
)

# ------------------------------------------------------------------------------
# Readiness / interpretation policy
# ------------------------------------------------------------------------------

policy <- data.frame(
  item = c(
    "working_resolution",
    "working_cluster_column",
    "target_clusters",
    "sample_specificity_policy",
    "exclusion_policy",
    "ambient_policy",
    "source_RDS_policy"
  ),
  value = c(
    "0.8",
    CLUSTER_COL,
    paste(TARGET_CLUSTERS, collapse = ";"),
    "Sample-specificity alone is not an exclusion criterion",
    "Require coherent foreign-lineage signal and/or mixed-lineage high-complexity evidence with weak HSC identity",
    "Ambient RNA will be tested later using whole-liver background; not corrected here",
    "Never overwrite source RDS"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  policy,
  file.path(PROVDIR, "analysis_policy_v1.0.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(PROVDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Working res0.8 clusters: ", length(clusters))
message("Target clusters: ", paste(TARGET_CLUSTERS, collapse = ", "))
message("Review first:")
message("  02_cluster_QC_summary.tsv")
message("  03_cluster_x_HSC3state.tsv")
message("  06_direct_gene_expression_by_cluster.tsv")
message("  07_TARGET_cluster14_16_15_direct_markers.tsv")
message("  08_top30_markers_working_res0.8.tsv")
message("  09_TARGET_cluster_centroid_neighbors.tsv")
message("  10_TARGET_cluster14_16_15_summary.tsv")
message("============================================================")
