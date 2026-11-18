#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_Resolution_Provenance_Audit_v1.0.R
#
# PURPOSE
# -------
# Provenance-anchored audit of the refined Mouse HSC clustering before any
# contamination / ambient-RNA cleanup.
#
# Audits the STORED HSC-specific RPCA clustering solutions:
#   res0.4, res0.6, res0.8, res1.0, res1.2
#
# Key questions:
#   1) Which stored resolution is the working HSC clustering?
#   2) Is the working clustering exactly identical to res0.8 cell-by-cell?
#   3) How do neighboring resolutions compare in:
#        - cluster number / tiny clusters
#        - PCA/RPCA silhouette
#        - sample mixing / single-sample dominance
#        - marker separability
#        - transitions relative to the working solution
#
# IMPORTANT
# ---------
# - Source RDS is NEVER overwritten.
# - No reclustering is performed.
# - Stored cluster labels are audited as-is.
# - Resolution is NOT selected based on Sham/Tx biological direction.
# - State labels (qHSC/ECM-activated/Contractile) are NOT used to choose the
#   resolution.
#
# Version: v1.0.1
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(cluster)
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
  "Mouse_HSC_Resolution_Provenance_Audit_v1.0.1"
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

if (!file.exists(HSC_RDS)) {
  stop("HSC RDS not found: ", HSC_RDS)
}

message("Loading HSC object...")
obj <- readRDS(HSC_RDS)
md <- obj@meta.data

RES_COLS <- c(
  "0.4" = "hsc_rpca_res_0_4_v610",
  "0.6" = "hsc_rpca_res_0_6_v610",
  "0.8" = "hsc_rpca_res_0_8_v610",
  "1.0" = "hsc_rpca_res_1_v610",
  "1.2" = "hsc_rpca_res_1_2_v610"
)

WORKING_COL <- "hsc_cluster_working_v610"
SAMPLE_COL <- "sample_hsc_v610"

required <- c(unname(RES_COLS), WORKING_COL, SAMPLE_COL)

missing_cols <- setdiff(required, colnames(md))
if (length(missing_cols)) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

REDUCTION <- if ("integrated.hsc.rpca" %in% names(obj@reductions)) {
  "integrated.hsc.rpca"
} else if ("pca.hsc" %in% names(obj@reductions)) {
  "pca.hsc"
} else {
  stop("No HSC-specific RPCA/PCA reduction found.")
}

GRAPH <- if ("hsc_snn_v610" %in% names(obj@graphs)) {
  "hsc_snn_v610"
} else {
  NA_character_
}

# ------------------------------------------------------------------------------
# Adjusted Rand Index (dependency-free implementation)
# ------------------------------------------------------------------------------

adjusted_rand_index <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)

  if (length(x) != length(y)) stop("ARI vectors differ in length.")

  tab <- table(x, y)
  n <- sum(tab)

  choose2 <- function(z) z * (z - 1) / 2

  sum_ij <- sum(choose2(tab))
  sum_i <- sum(choose2(rowSums(tab)))
  sum_j <- sum(choose2(colSums(tab)))
  total <- choose2(n)

  if (total == 0) return(NA_real_)

  expected <- (sum_i * sum_j) / total
  max_index <- 0.5 * (sum_i + sum_j)

  denom <- max_index - expected

  if (denom == 0) {
    return(ifelse(sum_ij == expected, 1, NA_real_))
  }

  (sum_ij - expected) / denom
}

# ------------------------------------------------------------------------------
# Provenance identity checks
# ------------------------------------------------------------------------------

working <- as.character(md[[WORKING_COL]])
res08 <- as.character(md[[RES_COLS["0.8"]]])
seurat_clusters <- if ("seurat_clusters" %in% colnames(md)) {
  as.character(md$seurat_clusters)
} else {
  rep(NA_character_, nrow(md))
}
res12 <- as.character(md[[RES_COLS["1.2"]]])

identity_checks <- data.frame(
  comparison = c(
    "working_vs_res0.8",
    "seurat_clusters_vs_res1.2"
  ),
  exact_cellwise_identity = c(
    identical(working, res08),
    if (all(is.na(seurat_clusters))) NA else identical(seurat_clusters, res12)
  ),
  ARI = c(
    adjusted_rand_index(working, res08),
    if (all(is.na(seurat_clusters))) NA_real_ else adjusted_rand_index(seurat_clusters, res12)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  identity_checks,
  file.path(TABDIR, "01_provenance_identity_checks.tsv")
)

# ------------------------------------------------------------------------------
# Common embedding for silhouette
# ------------------------------------------------------------------------------

emb <- Embeddings(obj, reduction = REDUCTION)

NDIMS <- min(30L, ncol(emb))
emb <- emb[, seq_len(NDIMS), drop = FALSE]

# Use the SAME stratified sample of cells for all resolutions.
# Full 13k silhouette would require a large O(n^2) distance matrix.
set.seed(20261118)

MAX_SIL_CELLS <- 5000L
all_cells <- rownames(md)

if (length(all_cells) <= MAX_SIL_CELLS) {
  sil_cells <- all_cells
} else {
  # Stratify by working cluster so that small clusters are represented.
  spl <- split(all_cells, working)
  total_n <- length(all_cells)

  target <- vapply(
    spl,
    function(v) as.integer(max(20L, round(MAX_SIL_CELLS * length(v) / total_n))),
    integer(1)
  )

  # Cap at cluster size.
  target <- pmin(target, vapply(spl, length, integer(1)))

  selected <- unlist(
    Map(
      function(v, n) {
        if (length(v) <= n) v else sample(v, n)
      },
      spl,
      target
    ),
    use.names = FALSE
  )

  if (length(selected) > MAX_SIL_CELLS) {
    selected <- sample(selected, MAX_SIL_CELLS)
  }

  sil_cells <- selected
}

emb_sil <- emb[sil_cells, , drop = FALSE]
dist_sil <- dist(emb_sil)

# ------------------------------------------------------------------------------
# RNA normalized-data matrix for marker separability
# ------------------------------------------------------------------------------

rna_data <- LayerData(obj[["RNA"]], layer = "data")

if (!inherits(rna_data, "Matrix")) {
  rna_data <- Matrix(rna_data, sparse = TRUE)
}

if (!identical(colnames(rna_data), colnames(obj))) {
  rna_data <- rna_data[, colnames(obj), drop = FALSE]
}

# Remove extremely sparse genes from marker ranking to limit noise.
gene_detected <- Matrix::rowSums(rna_data > 0)
marker_genes <- rownames(rna_data)[gene_detected >= max(20L, round(0.002 * ncol(obj)))]

rna_marker <- rna_data[marker_genes, , drop = FALSE]

# ------------------------------------------------------------------------------
# Helper metrics
# ------------------------------------------------------------------------------

normalized_entropy <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(NA_real_)
  p <- x / sum(x)
  h <- -sum(p * log(p))
  hmax <- log(length(x))
  if (!is.finite(hmax) || hmax == 0) return(0)
  h / hmax
}

cluster_sample_metrics <- function(labels, sample) {
  labs <- sort(unique(as.character(labels)))

  rows <- lapply(labs, function(cl) {
    ix <- which(as.character(labels) == cl)
    tab <- table(as.character(sample[ix]))

    data.frame(
      cluster = cl,
      n_cells = length(ix),
      n_samples_present = sum(tab > 0),
      max_sample_fraction = max(tab) / sum(tab),
      sample_entropy = normalized_entropy(as.numeric(tab)),
      dominant_sample = names(tab)[which.max(tab)],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

top_marker_delta <- function(labels, mat, top_n = 10L) {
  labels <- as.character(labels)
  labs <- sort(unique(labels))

  global_sum <- Matrix::rowSums(mat)
  n_total <- ncol(mat)

  rows <- list()
  k <- 0L

  for (cl in labs) {
    ix <- which(labels == cl)
    n_in <- length(ix)
    n_out <- n_total - n_in

    if (n_in < 1 || n_out < 1) next

    mean_in <- Matrix::rowMeans(mat[, ix, drop = FALSE])
    sum_in <- mean_in * n_in
    mean_out <- (global_sum - sum_in) / n_out

    delta <- as.numeric(mean_in - mean_out)
    names(delta) <- rownames(mat)

    ord <- order(delta, decreasing = TRUE)
    ord <- ord[seq_len(min(top_n, length(ord)))]

    k <- k + 1L
    rows[[k]] <- data.frame(
      cluster = cl,
      rank = seq_along(ord),
      gene = names(delta)[ord],
      mean_logexpr_delta = delta[ord],
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) {
    return(data.frame(
      cluster = character(0),
      rank = integer(0),
      gene = character(0),
      mean_logexpr_delta = numeric(0)
    ))
  }

  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# Resolution summaries
# ------------------------------------------------------------------------------

summary_rows <- list()
sample_metric_all <- list()
marker_all <- list()

sample_vec <- as.character(md[[SAMPLE_COL]])

for (res in names(RES_COLS)) {
  cn <- RES_COLS[[res]]
  lab <- as.character(md[[cn]])

  sizes <- table(lab)
  cs <- cluster_sample_metrics(lab, sample_vec)

  # Silhouette on the common sampled cells.
  sil_lab <- factor(lab[match(sil_cells, rownames(md))])

  if (nlevels(sil_lab) >= 2) {
    sil <- silhouette(as.integer(sil_lab), dist_sil)
    sil_mean <- mean(sil[, "sil_width"])
    sil_median <- median(sil[, "sil_width"])
    frac_negative <- mean(sil[, "sil_width"] < 0)
  } else {
    sil_mean <- NA_real_
    sil_median <- NA_real_
    frac_negative <- NA_real_
  }

  # Marker separability.
  tm <- top_marker_delta(lab, rna_marker, top_n = 10L)
  tm$resolution <- res
  marker_all[[res]] <- tm

  cluster_marker_means <- aggregate(
    mean_logexpr_delta ~ cluster,
    data = tm,
    FUN = mean
  )

  summary_rows[[res]] <- data.frame(
    resolution = as.numeric(res),
    cluster_column = cn,
    n_clusters = length(sizes),
    min_cluster_size = min(sizes),
    median_cluster_size = median(as.numeric(sizes)),
    max_cluster_size = max(sizes),
    n_clusters_lt_50 = sum(sizes < 50),
    n_clusters_lt_100 = sum(sizes < 100),
    silhouette_mean = sil_mean,
    silhouette_median = sil_median,
    fraction_negative_silhouette = frac_negative,
    ARI_vs_working = adjusted_rand_index(lab, working),
    median_sample_entropy = median(cs$sample_entropy, na.rm = TRUE),
    n_clusters_sample_dominance_gt_0.90 = sum(cs$max_sample_fraction > 0.90),
    median_top10_marker_delta = median(
      cluster_marker_means$mean_logexpr_delta,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  cs$resolution <- as.numeric(res)
  sample_metric_all[[res]] <- cs
}

resolution_summary <- do.call(rbind, summary_rows)
sample_metric_all <- do.call(rbind, sample_metric_all)
marker_all <- do.call(rbind, marker_all)

write_tsv(
  resolution_summary,
  file.path(TABDIR, "02_resolution_summary.tsv")
)

write_tsv(
  sample_metric_all,
  file.path(TABDIR, "03_cluster_sample_reproducibility.tsv")
)

write_tsv(
  marker_all,
  file.path(TABDIR, "04_top10_cluster_markers_by_resolution.tsv")
)

# ------------------------------------------------------------------------------
# Transition tables between neighboring resolutions
# ------------------------------------------------------------------------------

transition_pairs <- list(
  c("0.4", "0.6"),
  c("0.6", "0.8"),
  c("0.8", "1.0"),
  c("1.0", "1.2")
)

transition_rows <- list()
k <- 0L

for (pair in transition_pairs) {
  a <- pair[1]
  b <- pair[2]

  la <- as.character(md[[RES_COLS[[a]]]])
  lb <- as.character(md[[RES_COLS[[b]]]])

  tab <- as.data.frame(
    table(parent_cluster = la, child_cluster = lb),
    stringsAsFactors = FALSE
  )

  tab <- tab[tab$Freq > 0, , drop = FALSE]

  parent_total <- aggregate(
    Freq ~ parent_cluster,
    data = tab,
    FUN = sum
  )
  names(parent_total)[2] <- "parent_n"

  tab <- merge(
    tab,
    parent_total,
    by = "parent_cluster",
    all.x = TRUE,
    sort = FALSE
  )

  tab$fraction_of_parent <- tab$Freq / tab$parent_n
  tab$from_resolution <- as.numeric(a)
  tab$to_resolution <- as.numeric(b)

  k <- k + 1L
  transition_rows[[k]] <- tab[
    ,
    c(
      "from_resolution",
      "to_resolution",
      "parent_cluster",
      "child_cluster",
      "Freq",
      "parent_n",
      "fraction_of_parent"
    )
  ]
}

transitions <- do.call(rbind, transition_rows)

write_tsv(
  transitions,
  file.path(TABDIR, "05_neighbor_resolution_transitions.tsv")
)

# Summary of fragmentation from each parent
fragmentation <- do.call(
  rbind,
  lapply(
    split(
      transitions,
      interaction(
        transitions$from_resolution,
        transitions$to_resolution,
        transitions$parent_cluster,
        drop = TRUE
      )
    ),
    function(d) {
      data.frame(
        from_resolution = d$from_resolution[1],
        to_resolution = d$to_resolution[1],
        parent_cluster = d$parent_cluster[1],
        parent_n = d$parent_n[1],
        n_child_clusters = nrow(d),
        largest_child_fraction = max(d$fraction_of_parent),
        second_child_fraction = if (nrow(d) >= 2) {
          sort(d$fraction_of_parent, decreasing = TRUE)[2]
        } else {
          0
        },
        stringsAsFactors = FALSE
      )
    }
  )
)

write_tsv(
  fragmentation,
  file.path(TABDIR, "06_neighbor_resolution_fragmentation.tsv")
)

# ------------------------------------------------------------------------------
# Working-cluster sample distribution
# ------------------------------------------------------------------------------

working_sample <- cluster_sample_metrics(working, sample_vec)

write_tsv(
  working_sample,
  file.path(TABDIR, "07_working_res0.8_sample_distribution.tsv")
)

# ------------------------------------------------------------------------------
# Current HSC state × working cluster cross-tab
# (For interpretation only; NOT used in resolution selection.)
# ------------------------------------------------------------------------------

if ("HSC_state3_v612" %in% colnames(md)) {
  cross3 <- as.data.frame(
    table(
      working_cluster = working,
      HSC_state3 = as.character(md$HSC_state3_v612)
    ),
    stringsAsFactors = FALSE
  )
  cross3 <- cross3[cross3$Freq > 0, , drop = FALSE]

  write_tsv(
    cross3,
    file.path(TABDIR, "08_working_cluster_x_HSC3state.tsv")
  )
}

if ("HSC_state2_v612" %in% colnames(md)) {
  cross2 <- as.data.frame(
    table(
      working_cluster = working,
      HSC_state2 = as.character(md$HSC_state2_v612)
    ),
    stringsAsFactors = FALSE
  )
  cross2 <- cross2[cross2$Freq > 0, , drop = FALSE]

  write_tsv(
    cross2,
    file.path(TABDIR, "09_working_cluster_x_HSC2state.tsv")
  )
}

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

settings <- data.frame(
  item = c(
    "version",
    "input_RDS",
    "HSC_specific_reduction",
    "HSC_specific_graph",
    "working_cluster_column",
    "stored_resolution_columns",
    "silhouette_dimensions",
    "silhouette_sample_cells",
    "resolution_selection_policy",
    "state_label_policy",
    "source_object_policy"
  ),
  value = c(
    "v1.0.1",
    HSC_RDS,
    REDUCTION,
    GRAPH,
    WORKING_COL,
    paste(names(RES_COLS), unname(RES_COLS), sep = "=", collapse = ";"),
    NDIMS,
    length(sil_cells),
    "Topology/stability/sample reproducibility/marker separability; not Sham/Tx direction",
    "qHSC/ECM-activated/Contractile labels are not used to choose resolution",
    "Source RDS is never overwritten"
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
message("Reduction: ", REDUCTION)
message("Graph: ", GRAPH)
message("Working vs res0.8 exact: ", identity_checks$exact_cellwise_identity[1])
message("Working vs res0.8 ARI: ", identity_checks$ARI[1])
message("Review first:")
message("  01_provenance_identity_checks.tsv")
message("  02_resolution_summary.tsv")
message("  03_cluster_sample_reproducibility.tsv")
message("  04_top10_cluster_markers_by_resolution.tsv")
message("  05_neighbor_resolution_transitions.tsv")
message("  06_neighbor_resolution_fragmentation.tsv")
message("  07_working_res0.8_sample_distribution.tsv")
message("============================================================")
