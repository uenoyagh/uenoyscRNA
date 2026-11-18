#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_ShamsTx_GeneLevel_Decomposition_v1.0.R
#
# PURPOSE
# -------
# Resolve the apparent divergence observed in the module-level analysis:
#
#   - Anti-inflammatory / resolution score: increased in Tx
#   - Efferocytosis / resolution score: increased in Tx
#   - MASH-pathogenic inflammatory score: decreased in Tx
#   - HSC-activating ligand module: ~unchanged
#   - "Fibrosis-associated MΦ" and ECM-remodeling modules: increased in Tx
#   - IL10 combined axis: not increased at pseudobulk level
#
# This script decomposes those modules to INDIVIDUAL GENES and identifies
# which genes / frozen clusters / FINAL macrophage classes drive the changes.
#
# IMPORTANT:
# - This is descriptive sample-level analysis (Sham n=2; Tx n=2).
# - No cell-level inferential p-values.
# - No formal DE p-values in this version.
# - Original filtered RNA counts are primary.
# - DecontX remains sensitivity-only.
# - The term "fibrosis-associated" is treated as a state label, not proof of
#   net profibrotic activity.
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mouse_Mphi_PrimaryClean_v5.0",
  "objects",
  "Mouse_MASH_Mphi_PrimaryClean_v5.0.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_Mphi",
  "Mouse_PrimaryClean_ShamsTx_GeneLevel_Decomposition_v1.0"
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
PROVDIR <- file.path(OUTDIR, "provenance")

for (d in c(OUTDIR, TABDIR, FIGDIR, PROVDIR)) {
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

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

direction_label <- function(sh, tx) {
  sh <- sh[is.finite(sh)]
  tx <- tx[is.finite(tx)]

  if (!length(sh) || !length(tx)) return("insufficient")
  if (min(tx) > max(sh)) return("Tx_above_all_Sham")
  if (max(tx) < min(sh)) return("Tx_below_all_Sham")

  d <- safe_mean(tx) - safe_mean(sh)
  if (d > 0) return("Tx_higher_overlap")
  if (d < 0) return("Tx_lower_overlap")
  "No_change"
}

# ------------------------------------------------------------------------------
# Load
# ------------------------------------------------------------------------------

if (!file.exists(RDS)) stop("RDS not found: ", RDS)

message("Loading: ", RDS)
obj <- readRDS(RDS)

md <- obj@meta.data

required_meta <- c(
  "sample_display_FIXED2",
  "cluster_res2",
  "qc_primary_include",
  "macrophage_class_Res2_FINAL_v4145"
)

missing_meta <- setdiff(required_meta, colnames(md))
if (length(missing_meta)) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

if (!all(md$qc_primary_include %in% TRUE)) {
  stop("Primary-Clean object contains qc_primary_include != TRUE.")
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(
  Sham1 = "Sham",
  Sham20 = "Sham",
  Tx17 = "Tx",
  Tx5 = "Tx"
)

sample_all <- as.character(md$sample_display_FIXED2)
use_cells <- rownames(md)[sample_all %in% SAMPLES]

obj_st <- subset(obj, cells = use_cells)
md_st <- obj_st@meta.data

sample <- as.character(md_st$sample_display_FIXED2)
cluster <- as.character(md_st$cluster_res2)
final_class <- as.character(md_st$macrophage_class_Res2_FINAL_v4145)

# ------------------------------------------------------------------------------
# Gene sets to decompose
# ------------------------------------------------------------------------------

GENE_GROUPS <- list(
  IL10_receptor_signaling = c(
    "Il10ra", "Il10rb", "Jak1", "Tyk2", "Stat3"
  ),
  IL10_downstream_response = c(
    "Socs3", "Bcl3", "Il1rn"
  ),
  Anti_inflammatory_resolution = c(
    "Il10", "Il1rn", "Areg", "Gas6", "Mfge8",
    "Igf1", "Hgf", "Apoe", "Gdf15", "Edil3", "Serpina1"
  ),
  Fibrosis_associated_state = c(
    "Spp1", "Trem2", "Gpnmb", "Lgals3",
    "Timp1", "Fn1", "Tgfb1", "Pdgfb"
  ),
  Direct_HSC_activating_candidates = c(
    "Tgfb1", "Pdgfa", "Pdgfb", "Spp1"
  ),
  ECM_remodeling = c(
    "Spp1", "Timp1", "Lgals3", "Mmp9",
    "Mmp12", "Ctsb", "Ctsd", "Ctsk"
  ),
  MASH_pathogenic_inflammatory = c(
    "F13a1", "Olr1", "Selplg"
  ),
  IFN_inflammatory = c(
    "Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"
  ),
  Efferocytosis_resolution = c(
    "Mertk", "Axl", "Gas6", "Mfge8"
  )
)

gene_map <- do.call(
  rbind,
  lapply(names(GENE_GROUPS), function(gp) {
    data.frame(
      mechanism = gp,
      gene = GENE_GROUPS[[gp]],
      stringsAsFactors = FALSE
    )
  })
)

gene_map <- unique(gene_map)

# ------------------------------------------------------------------------------
# Counts / CP10k / logCP10k
# ------------------------------------------------------------------------------

counts <- LayerData(obj_st[["RNA"]], layer = "counts")

if (!inherits(counts, "Matrix")) {
  counts <- Matrix(counts, sparse = TRUE)
}

if (!identical(colnames(counts), colnames(obj_st))) {
  counts <- counts[, colnames(obj_st), drop = FALSE]
}

present_genes <- intersect(unique(gene_map$gene), rownames(counts))
missing_genes <- setdiff(unique(gene_map$gene), rownames(counts))

coverage <- gene_map
coverage$present <- coverage$gene %in% rownames(counts)

write_tsv(
  coverage,
  file.path(TABDIR, "01_gene_coverage.tsv")
)

libsize <- Matrix::colSums(counts)
if (any(libsize <= 0)) stop("Zero library size cell.")

cp10k <- t(t(counts) / libsize * 1e4)
logcp10k <- log1p(cp10k)

# ------------------------------------------------------------------------------
# A. Sample-level pseudobulk gene expression
# ------------------------------------------------------------------------------

pb_rows <- list()
k <- 0L

for (s in SAMPLES) {
  cells_s <- colnames(obj_st)[sample == s]
  if (!length(cells_s)) next

  m <- counts[, cells_s, drop = FALSE]
  pb <- Matrix::rowSums(m)
  pb_cp10k <- pb / sum(pb) * 1e4

  for (g in present_genes) {
    k <- k + 1L
    pb_rows[[k]] <- data.frame(
      sample = s,
      group = unname(GROUP[s]),
      n_cells = length(cells_s),
      gene = g,
      CP10k = as.numeric(pb_cp10k[g]),
      log1p_CP10k = log1p(as.numeric(pb_cp10k[g])),
      stringsAsFactors = FALSE
    )
  }
}

pb_gene <- do.call(rbind, pb_rows)

write_tsv(
  pb_gene,
  file.path(TABDIR, "02_sample_pseudobulk_gene_expression.tsv")
)

# Gene-level group effects
gene_effect_rows <- list()
k <- 0L

for (g in present_genes) {
  d <- pb_gene[pb_gene$gene == g, , drop = FALSE]

  sh <- d$CP10k[d$group == "Sham"]
  tx <- d$CP10k[d$group == "Tx"]

  sham_mean <- safe_mean(sh)
  tx_mean <- safe_mean(tx)

  log2fc <- log2((tx_mean + 0.1) / (sham_mean + 0.1))

  k <- k + 1L
  gene_effect_rows[[k]] <- data.frame(
    gene = g,
    sham_mean_CP10k = sham_mean,
    sham_sd = safe_sd(sh),
    tx_mean_CP10k = tx_mean,
    tx_sd = safe_sd(tx),
    tx_minus_sham_CP10k = tx_mean - sham_mean,
    tx_over_sham = ifelse(sham_mean > 0, tx_mean / sham_mean, NA_real_),
    log2FC_Tx_vs_Sham_pseudo0.1 = log2fc,
    directional_consistency = direction_label(sh, tx),
    stringsAsFactors = FALSE
  )
}

gene_effect <- do.call(rbind, gene_effect_rows)

# Attach mechanism membership; a gene may appear in multiple mechanisms.
gene_effect_long <- merge(
  gene_map[gene_map$gene %in% present_genes, , drop = FALSE],
  gene_effect,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

gene_effect_long <- gene_effect_long[
  order(
    gene_effect_long$mechanism,
    -abs(gene_effect_long$log2FC_Tx_vs_Sham_pseudo0.1)
  ),
  ,
  drop = FALSE
]

write_tsv(
  gene_effect_long,
  file.path(TABDIR, "03_gene_level_pseudobulk_effects_by_mechanism.tsv")
)

# ------------------------------------------------------------------------------
# B. Cell-level detection / mean logCP10k by sample
# ------------------------------------------------------------------------------
cell_gene_rows <- list()
k <- 0L

for (s in SAMPLES) {
  ix <- which(sample == s)
  if (!length(ix)) next

  for (g in present_genes) {
    v <- as.numeric(logcp10k[g, ix])
    cnt <- as.numeric(counts[g, ix])

    k <- k + 1L
    cell_gene_rows[[k]] <- data.frame(
      sample = s,
      group = unname(GROUP[s]),
      gene = g,
      n_cells = length(ix),
      pct_detected = 100 * mean(cnt > 0),
      mean_log1p_CP10k = safe_mean(v),
      median_log1p_CP10k = median(v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

cell_gene <- do.call(rbind, cell_gene_rows)

write_tsv(
  cell_gene,
  file.path(TABDIR, "04_sample_cell_gene_expression.tsv")
)

# ------------------------------------------------------------------------------
# C. FINAL-class gene contribution
#
# Contribution = class fraction within sample × class mean log1p(CP10k)
# This separates composition and expression effects in a transparent way.
# ------------------------------------------------------------------------------

class_rows <- list()
k <- 0L

for (s in SAMPLES) {
  ix_s <- which(sample == s)
  n_s <- length(ix_s)

  classes <- sort(unique(final_class[ix_s]))

  for (cl in classes) {
    ix <- which(sample == s & final_class == cl)
    frac <- length(ix) / n_s

    for (g in present_genes) {
      score <- safe_mean(as.numeric(logcp10k[g, ix]))

      k <- k + 1L
      class_rows[[k]] <- data.frame(
        sample = s,
        group = unname(GROUP[s]),
        final_class = cl,
        gene = g,
        n_cells = length(ix),
        class_fraction = frac,
        class_mean_log1p_CP10k = score,
        weighted_contribution = frac * score,
        stringsAsFactors = FALSE
      )
    }
  }
}

class_gene <- do.call(rbind, class_rows)

write_tsv(
  class_gene,
  file.path(TABDIR, "05_finalclass_sample_gene_contribution.tsv")
)

class_effect_rows <- list()
k <- 0L

for (g in present_genes) {
  dg <- class_gene[class_gene$gene == g, , drop = FALSE]

  for (cl in sort(unique(dg$final_class))) {
    d <- dg[dg$final_class == cl, , drop = FALSE]

    full <- data.frame(
      sample = SAMPLES,
      group = unname(GROUP[SAMPLES]),
      stringsAsFactors = FALSE
    )

    full <- merge(
      full,
      d[, c("sample", "weighted_contribution")],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )

    full$weighted_contribution[is.na(full$weighted_contribution)] <- 0

    sh <- full$weighted_contribution[full$group == "Sham"]
    tx <- full$weighted_contribution[full$group == "Tx"]

    k <- k + 1L
    class_effect_rows[[k]] <- data.frame(
      gene = g,
      final_class = cl,
      sham_mean_contribution = safe_mean(sh),
      tx_mean_contribution = safe_mean(tx),
      tx_minus_sham_contribution = safe_mean(tx) - safe_mean(sh),
      abs_delta_contribution = abs(safe_mean(tx) - safe_mean(sh)),
      directional_consistency = direction_label(sh, tx),
      stringsAsFactors = FALSE
    )
  }
}

class_effect <- do.call(rbind, class_effect_rows)

write_tsv(
  class_effect,
  file.path(TABDIR, "06_finalclass_gene_attribution_effects.tsv")
)

# ------------------------------------------------------------------------------
# D. Frozen-cluster gene contribution
# ------------------------------------------------------------------------------

cluster_rows <- list()
k <- 0L

for (s in SAMPLES) {
  ix_s <- which(sample == s)
  n_s <- length(ix_s)

  clusters <- sort(unique(cluster[ix_s]))

  for (cl in clusters) {
    ix <- which(sample == s & cluster == cl)
    frac <- length(ix) / n_s

    for (g in present_genes) {
      score <- safe_mean(as.numeric(logcp10k[g, ix]))

      k <- k + 1L
      cluster_rows[[k]] <- data.frame(
        sample = s,
        group = unname(GROUP[s]),
        cluster = cl,
        gene = g,
        n_cells = length(ix),
        cluster_fraction = frac,
        cluster_mean_log1p_CP10k = score,
        weighted_contribution = frac * score,
        stringsAsFactors = FALSE
      )
    }
  }
}

cluster_gene <- do.call(rbind, cluster_rows)

write_tsv(
  cluster_gene,
  file.path(TABDIR, "07_cluster_sample_gene_contribution.tsv")
)

cluster_effect_rows <- list()
k <- 0L

for (g in present_genes) {
  dg <- cluster_gene[cluster_gene$gene == g, , drop = FALSE]

  for (cl in sort(unique(dg$cluster))) {
    d <- dg[dg$cluster == cl, , drop = FALSE]

    full <- data.frame(
      sample = SAMPLES,
      group = unname(GROUP[SAMPLES]),
      stringsAsFactors = FALSE
    )

    full <- merge(
      full,
      d[, c("sample", "weighted_contribution")],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )

    full$weighted_contribution[is.na(full$weighted_contribution)] <- 0

    sh <- full$weighted_contribution[full$group == "Sham"]
    tx <- full$weighted_contribution[full$group == "Tx"]

    k <- k + 1L
    cluster_effect_rows[[k]] <- data.frame(
      gene = g,
      cluster = cl,
      sham_mean_contribution = safe_mean(sh),
      tx_mean_contribution = safe_mean(tx),
      tx_minus_sham_contribution = safe_mean(tx) - safe_mean(sh),
      abs_delta_contribution = abs(safe_mean(tx) - safe_mean(sh)),
      directional_consistency = direction_label(sh, tx),
      stringsAsFactors = FALSE
    )
  }
}

cluster_effect <- do.call(rbind, cluster_effect_rows)

write_tsv(
  cluster_effect,
  file.path(TABDIR, "08_cluster_gene_attribution_effects.tsv")
)

# Top cluster driver per gene
top_cluster_gene <- do.call(
  rbind,
  lapply(split(cluster_effect, cluster_effect$gene), function(d) {
    head(d[order(-d$abs_delta_contribution), , drop = FALSE], 5)
  })
)

write_tsv(
  top_cluster_gene,
  file.path(TABDIR, "09_TOP_cluster_drivers_by_gene.tsv")
)

# ------------------------------------------------------------------------------
# E. Compact mechanistic decomposition tables
# ------------------------------------------------------------------------------

KEY_MECHANISMS <- c(
  "IL10_receptor_signaling",
  "IL10_downstream_response",
  "Anti_inflammatory_resolution",
  "Fibrosis_associated_state",
  "Direct_HSC_activating_candidates",
  "ECM_remodeling",
  "MASH_pathogenic_inflammatory",
  "Efferocytosis_resolution"
)

key_gene <- gene_effect_long[
  gene_effect_long$mechanism %in% KEY_MECHANISMS,
  ,
  drop = FALSE
]

write_tsv(
  key_gene,
  file.path(TABDIR, "10_KEY_gene_decomposition.tsv")
)

# For each mechanism, list strongest positive and negative gene changes.
mechanism_extremes <- do.call(
  rbind,
  lapply(split(key_gene, key_gene$mechanism), function(d) {
    d <- d[order(d$log2FC_Tx_vs_Sham_pseudo0.1), , drop = FALSE]

    lo <- head(d, 3)
    hi <- tail(d, 3)

    out <- rbind(lo, hi)
    out$extreme_type <- c(
      rep("most_negative", nrow(lo)),
      rep("most_positive", nrow(hi))
    )
    out
  })
)

write_tsv(
  mechanism_extremes,
  file.path(TABDIR, "11_mechanism_gene_extremes.tsv")
)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

plot_gene_group <- function(mechanism_name, filename) {
  d <- gene_effect_long[
    gene_effect_long$mechanism == mechanism_name,
    ,
    drop = FALSE
  ]

  d <- unique(d[, c(
    "gene",
    "log2FC_Tx_vs_Sham_pseudo0.1",
    "directional_consistency"
  )])

  d$gene <- factor(
    d$gene,
    levels = d$gene[order(d$log2FC_Tx_vs_Sham_pseudo0.1)]
  )

  p <- ggplot(
    d,
    aes(
      x = log2FC_Tx_vs_Sham_pseudo0.1,
      y = gene
    )
  ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_col() +
    labs(
      x = "log2FC Tx vs Sham (sample pseudobulk; pseudocount 0.1)",
      y = NULL,
      title = mechanism_name,
      subtitle = "Descriptive gene-level decomposition; n=2/group"
    ) +
    theme_classic(base_size = 10)

  ggsave(
    file.path(FIGDIR, filename),
    p,
    width = 7.5,
    height = max(4.5, 0.35 * nrow(d) + 2)
  )
}

plot_gene_group(
  "IL10_receptor_signaling",
  "01_IL10_receptor_gene_decomposition.pdf"
)

plot_gene_group(
  "IL10_downstream_response",
  "02_IL10_downstream_gene_decomposition.pdf"
)

plot_gene_group(
  "Anti_inflammatory_resolution",
  "03_Anti_inflammatory_resolution_gene_decomposition.pdf"
)

plot_gene_group(
  "Fibrosis_associated_state",
  "04_Fibrosis_associated_state_gene_decomposition.pdf"
)

plot_gene_group(
  "Direct_HSC_activating_candidates",
  "05_HSC_activating_gene_decomposition.pdf"
)

plot_gene_group(
  "ECM_remodeling",
  "06_ECM_remodeling_gene_decomposition.pdf"
)

plot_gene_group(
  "MASH_pathogenic_inflammatory",
  "07_MASH_pathogenic_gene_decomposition.pdf"
)

plot_gene_group(
  "Efferocytosis_resolution",
  "08_Efferocytosis_gene_decomposition.pdf"
)

# Sample × gene pseudobulk heatmap-like plot for key genes
heat <- pb_gene[
  pb_gene$gene %in% unique(key_gene$gene),
  ,
  drop = FALSE
]

heat$sample <- factor(heat$sample, levels = SAMPLES)

# z-score per gene across the four samples
heat$z <- ave(
  heat$log1p_CP10k,
  heat$gene,
  FUN = function(x) {
    sx <- sd(x)
    if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
    (x - mean(x)) / sx
  }
)

p_heat <- ggplot(
  heat,
  aes(
    x = sample,
    y = gene,
    fill = z
  )
) +
  geom_tile() +
  labs(
    x = NULL,
    y = NULL,
    fill = "within-gene\nz score",
    title = "Key mechanistic genes | sample pseudobulk",
    subtitle = "log1p(CP10k), standardized within each gene"
  ) +
  theme_classic(base_size = 8)

ggsave(
  file.path(FIGDIR, "09_KEY_gene_sample_heatmap.pdf"),
  p_heat,
  width = 7,
  height = max(8, 0.22 * length(unique(heat$gene)) + 2)
)

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

design <- data.frame(
  item = c(
    "input_RDS",
    "comparison",
    "replicate_unit",
    "inferential_policy",
    "expression_source",
    "gene_effect_definition",
    "class_contribution_definition",
    "cluster_contribution_definition",
    "fibrosis_label_caution",
    "IL10_label_caution"
  ),
  value = c(
    RDS,
    "Sham1/Sham20 vs Tx17/Tx5",
    "biological sample",
    "descriptive only; n=2/group; no formal p-values",
    "original filtered RNA counts in Primary-Clean MΦ object",
    "sample pseudobulk CP10k; descriptive log2FC with pseudocount 0.1",
    "class fraction × class mean log1p(CP10k)",
    "cluster fraction × cluster mean log1p(CP10k)",
    "fibrosis-associated state is not equivalent to net profibrotic activity",
    "IL10-axis-compatible expression is not proof of causal IL10 signaling"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  design,
  file.path(PROVDIR, "analysis_design_v1.0.tsv")
)

write_tsv(
  gene_map,
  file.path(PROVDIR, "gene_groups_v1.0.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(PROVDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  03_gene_level_pseudobulk_effects_by_mechanism.tsv")
message("  06_finalclass_gene_attribution_effects.tsv")
message("  09_TOP_cluster_drivers_by_gene.tsv")
message("  10_KEY_gene_decomposition.tsv")
message("  11_mechanism_gene_extremes.tsv")
message("Figures:")
message("  04_Fibrosis_associated_state_gene_decomposition.pdf")
message("  05_HSC_activating_gene_decomposition.pdf")
message("  06_ECM_remodeling_gene_decomposition.pdf")
message("  09_KEY_gene_sample_heatmap.pdf")
message("============================================================")
