#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_ShamsTx_IL10_Profibrotic_Attribution_v1.0.R
#
# PURPOSE
# -------
# Follow-up mechanistic analysis using the finalized Mouse Primary-Clean MΦ RDS.
#
# Primary questions:
#   1) Does Tx increase an IL10-axis-compatible response program?
#   2) Does Tx suppress or enhance profibrotic / ECM-promoting MΦ programs?
#   3) Which frozen clusters and FINAL macrophage classes contribute most to
#      the observed Sham -> Tx changes?
#
# IMPORTANT INTERPRETATION
# ------------------------
# - "IL10 response" here means IL10-axis-compatible transcriptional activity.
#   It is NOT proof of causal IL10 signaling because STAT3/SOCS3/IL1RN etc.
#   can be regulated by other pathways.
# - M1/M2 are not used as identities.
# - Cells are NOT treated as independent biological replicates.
# - Sham n=2 and Tx n=2; no formal inferential p-values are reported.
# - Primary results are sample-level effect sizes and directional consistency.
# - Original filtered RNA counts are primary; DecontX remains sensitivity-only.
#
# INPUT
# -----
# Mouse_MASH_Mphi_PrimaryClean_v5.0.rds
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
  "Mouse_PrimaryClean_ShamsTx_IL10_Profibrotic_Attribution_v1.0"
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

if (!file.exists(RDS)) {
  stop("Primary-clean Mouse MΦ RDS not found: ", RDS)
}

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
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

if (!all(md$qc_primary_include %in% TRUE)) {
  stop("Primary-Clean object contains qc_primary_include != TRUE.")
}

SAMPLE_ORDER <- c("Sham1", "Sham20", "Tx17", "Tx5")

sample_all <- as.character(md$sample_display_FIXED2)
use_cells <- rownames(md)[sample_all %in% SAMPLE_ORDER]

obj_st <- subset(obj, cells = use_cells)
md_st <- obj_st@meta.data

md_st$sample <- as.character(md_st$sample_display_FIXED2)
md_st$group <- ifelse(
  md_st$sample %in% c("Sham1", "Sham20"),
  "Sham",
  "Tx"
)

md_st$group <- factor(md_st$group, levels = c("Sham", "Tx"))
md_st$sample <- factor(md_st$sample, levels = SAMPLE_ORDER)

# ------------------------------------------------------------------------------
# Mechanistic panels
# ------------------------------------------------------------------------------

PANELS <- list(
  # IL10 signaling machinery / competence.
  IL10_receptor_signaling = c(
    "Il10ra", "Il10rb", "Jak1", "Tyk2", "Stat3"
  ),

  # More proximal downstream response. Still not IL10-specific.
  IL10_downstream_response = c(
    "Socs3", "Bcl3", "Il1rn"
  ),

  # Combined axis for compact screening.
  IL10_axis_combined = c(
    "Il10ra", "Il10rb", "Jak1", "Tyk2", "Stat3",
    "Socs3", "Bcl3", "Il1rn"
  ),

  # Resolution / anti-inflammatory mediators relevant to macrophage repair.
  Anti_inflammatory_resolution = c(
    "Il10", "Il1rn", "Areg", "Gas6", "Mfge8",
    "Igf1", "Hgf", "Apoe", "Gdf15", "Edil3", "Serpina1"
  ),

  # Fibrosis-associated macrophage state.
  Fibrosis_associated_Mphi = c(
    "Spp1", "Trem2", "Gpnmb", "Lgals3",
    "Timp1", "Fn1", "Tgfb1", "Pdgfb"
  ),

  # Candidate direct HSC-activating ligands.
  HSC_activating_ligands = c(
    "Tgfb1", "Pdgfa", "Pdgfb", "Spp1"
  ),

  # ECM remodeling / matrix turnover. This is not identical to net fibrosis.
  ECM_remodeling = c(
    "Spp1", "Timp1", "Lgals3", "Mmp9",
    "Mmp12", "Ctsb", "Ctsd", "Ctsk"
  ),

  # MASH inflammatory / pathogenic axis from the frozen project panel.
  MASH_pathogenic_inflammatory = c(
    "F13a1", "Olr1", "Selplg"
  ),

  IFN_inflammatory = c(
    "Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"
  ),

  # Efferocytosis / resolution axis retained for linkage to the previous step.
  Efferocytosis_resolution = c(
    "Mertk", "Axl", "Gas6", "Mfge8"
  )
)

# ------------------------------------------------------------------------------
# Primary RNA counts and log1p(CP10k)
# ------------------------------------------------------------------------------

counts <- LayerData(obj_st[["RNA"]], layer = "counts")

if (!inherits(counts, "Matrix")) {
  counts <- Matrix(counts, sparse = TRUE)
}

if (!identical(colnames(counts), colnames(obj_st))) {
  counts <- counts[, colnames(obj_st), drop = FALSE]
}

libsize <- Matrix::colSums(counts)

if (any(libsize <= 0)) {
  stop("Zero-library-size cell found.")
}

cp10k <- t(t(counts) / libsize * 1e4)
logcp10k <- log1p(cp10k)

# ------------------------------------------------------------------------------
# Coverage audit
# ------------------------------------------------------------------------------

coverage <- do.call(
  rbind,
  lapply(names(PANELS), function(pn) {
    g <- PANELS[[pn]]
    present <- intersect(g, rownames(counts))
    missing <- setdiff(g, present)

    data.frame(
      panel = pn,
      n_requested = length(g),
      n_present = length(present),
      coverage_fraction = length(present) / length(g),
      genes_present = paste(present, collapse = ";"),
      genes_missing = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  coverage,
  file.path(TABDIR, "01_panel_coverage.tsv")
)

# ------------------------------------------------------------------------------
# Cell-level scores
# ------------------------------------------------------------------------------

score_panel <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[g, , drop = FALSE]))
}

cell <- data.frame(
  cell = colnames(obj_st),
  sample = as.character(md_st$sample),
  group = as.character(md_st$group),
  cluster = as.character(md_st$cluster_res2),
  final_class = as.character(md_st$macrophage_class_Res2_FINAL_v4145),
  stringsAsFactors = FALSE
)

for (pn in names(PANELS)) {
  cell[[pn]] <- score_panel(logcp10k, PANELS[[pn]])
}

write_tsv(
  cell,
  file.path(TABDIR, "02_cell_mechanistic_scores.tsv")
)

# ------------------------------------------------------------------------------
# Sample-level means / medians
# ------------------------------------------------------------------------------

sample_rows <- list()
k <- 0L

for (s in SAMPLE_ORDER) {
  ix <- which(cell$sample == s)
  if (!length(ix)) next

  gp <- unique(cell$group[ix])[1]

  for (pn in names(PANELS)) {
    x <- cell[[pn]][ix]

    k <- k + 1L
    sample_rows[[k]] <- data.frame(
      sample = s,
      group = gp,
      n_cells = length(ix),
      metric = pn,
      mean_cell_score = safe_mean(x),
      median_cell_score = median(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

sample_score <- do.call(rbind, sample_rows)

write_tsv(
  sample_score,
  file.path(TABDIR, "03_sample_mechanistic_scores.tsv")
)

# ------------------------------------------------------------------------------
# Group-level descriptive effects
# ------------------------------------------------------------------------------

effect_rows <- list()
k <- 0L

for (pn in names(PANELS)) {
  d <- sample_score[sample_score$metric == pn, , drop = FALSE]

  sh <- d$mean_cell_score[d$group == "Sham"]
  tx <- d$mean_cell_score[d$group == "Tx"]

  k <- k + 1L
  effect_rows[[k]] <- data.frame(
    metric = pn,
    n_sham_samples = length(sh),
    n_tx_samples = length(tx),
    sham_mean = safe_mean(sh),
    sham_sd = safe_sd(sh),
    tx_mean = safe_mean(tx),
    tx_sd = safe_sd(tx),
    tx_minus_sham = safe_mean(tx) - safe_mean(sh),
    percent_change_vs_sham = 100 * (
      safe_mean(tx) - safe_mean(sh)
    ) / abs(safe_mean(sh)),
    directional_consistency = direction_label(sh, tx),
    inferential_p_value = NA_real_,
    statistical_note = "Descriptive sample-level effect only; n=2/group",
    stringsAsFactors = FALSE
  )
}

group_effect <- do.call(rbind, effect_rows)

write_tsv(
  group_effect,
  file.path(TABDIR, "04_group_mechanistic_effects.tsv")
)

# ------------------------------------------------------------------------------
# Sample pseudobulk panel expression
# ------------------------------------------------------------------------------

pb_rows <- list()
k <- 0L

for (s in SAMPLE_ORDER) {
  cells_s <- colnames(obj_st)[as.character(md_st$sample) == s]
  if (!length(cells_s)) next

  m <- counts[, cells_s, drop = FALSE]
  pb <- Matrix::rowSums(m)
  pb_cp10k <- pb / sum(pb) * 1e4

  for (pn in names(PANELS)) {
    g <- intersect(PANELS[[pn]], names(pb_cp10k))
    val <- if (length(g)) mean(pb_cp10k[g]) else NA_real_

    k <- k + 1L
    pb_rows[[k]] <- data.frame(
      sample = s,
      group = ifelse(s %in% c("Sham1", "Sham20"), "Sham", "Tx"),
      n_cells = length(cells_s),
      metric = pn,
      mean_gene_CP10k = val,
      stringsAsFactors = FALSE
    )
  }
}

pb_score <- do.call(rbind, pb_rows)

write_tsv(
  pb_score,
  file.path(TABDIR, "05_sample_pseudobulk_mechanistic_scores.tsv")
)

pb_effect_rows <- list()
k <- 0L

for (pn in names(PANELS)) {
  d <- pb_score[pb_score$metric == pn, , drop = FALSE]

  sh <- d$mean_gene_CP10k[d$group == "Sham"]
  tx <- d$mean_gene_CP10k[d$group == "Tx"]

  k <- k + 1L
  pb_effect_rows[[k]] <- data.frame(
    metric = pn,
    sham_mean_CP10k = safe_mean(sh),
    sham_sd = safe_sd(sh),
    tx_mean_CP10k = safe_mean(tx),
    tx_sd = safe_sd(tx),
    tx_minus_sham_CP10k = safe_mean(tx) - safe_mean(sh),
    tx_over_sham = ifelse(
      is.finite(safe_mean(sh)) && safe_mean(sh) != 0,
      safe_mean(tx) / safe_mean(sh),
      NA_real_
    ),
    directional_consistency = direction_label(sh, tx),
    inferential_p_value = NA_real_,
    statistical_note = "Descriptive pseudobulk sample-level effect only; n=2/group",
    stringsAsFactors = FALSE
  )
}

pb_effect <- do.call(rbind, pb_effect_rows)

write_tsv(
  pb_effect,
  file.path(TABDIR, "06_group_pseudobulk_mechanistic_effects.tsv")
)

# ------------------------------------------------------------------------------
# Cluster attribution
#
# For each sample × cluster:
#   contribution = cluster fraction within sample × cluster mean score
#
# This decomposes the sample-wide mean score into cluster-weighted components.
# ------------------------------------------------------------------------------
cluster_rows <- list()
k <- 0L

for (s in SAMPLE_ORDER) {
  ix_s <- which(cell$sample == s)
  n_s <- length(ix_s)
  if (!n_s) next

  cls <- sort(unique(cell$cluster[ix_s]))

  for (cl in cls) {
    ix <- which(cell$sample == s & cell$cluster == cl)
    frac <- length(ix) / n_s

    for (pn in names(PANELS)) {
      mscore <- safe_mean(cell[[pn]][ix])

      k <- k + 1L
      cluster_rows[[k]] <- data.frame(
        sample = s,
        group = unique(cell$group[ix])[1],
        cluster = cl,
        n_cells = length(ix),
        sample_total_cells = n_s,
        cluster_fraction = frac,
        metric = pn,
        cluster_mean_score = mscore,
        weighted_contribution = frac * mscore,
        stringsAsFactors = FALSE
      )
    }
  }
}

cluster_sample <- do.call(rbind, cluster_rows)

write_tsv(
  cluster_sample,
  file.path(TABDIR, "07_cluster_sample_mechanistic_contribution.tsv")
)

cluster_effect_rows <- list()
k <- 0L

for (pn in names(PANELS)) {
  d0 <- cluster_sample[cluster_sample$metric == pn, , drop = FALSE]

  for (cl in sort(unique(d0$cluster))) {
    d <- d0[d0$cluster == cl, , drop = FALSE]

    # Explicitly add zeros for a cluster absent from a sample:
    full <- data.frame(
      sample = SAMPLE_ORDER,
      group = c("Sham", "Sham", "Tx", "Tx"),
      stringsAsFactors = FALSE
    )

    full <- merge(
      full,
      d[, c(
        "sample",
        "cluster_fraction",
        "cluster_mean_score",
        "weighted_contribution"
      )],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )

    full$cluster_fraction[is.na(full$cluster_fraction)] <- 0
    full$weighted_contribution[is.na(full$weighted_contribution)] <- 0

    sh <- full$weighted_contribution[full$group == "Sham"]
    tx <- full$weighted_contribution[full$group == "Tx"]

    k <- k + 1L
    cluster_effect_rows[[k]] <- data.frame(
      metric = pn,
      cluster = cl,
      sham_mean_weighted_contribution = safe_mean(sh),
      tx_mean_weighted_contribution = safe_mean(tx),
      tx_minus_sham_contribution = safe_mean(tx) - safe_mean(sh),
      abs_delta_contribution = abs(safe_mean(tx) - safe_mean(sh)),
      directional_consistency = direction_label(sh, tx),
      stringsAsFactors = FALSE
    )
  }
}

cluster_effect <- do.call(rbind, cluster_effect_rows)

cluster_effect <- cluster_effect[
  order(
    cluster_effect$metric,
    -cluster_effect$abs_delta_contribution
  ),
  ,
  drop = FALSE
]

write_tsv(
  cluster_effect,
  file.path(TABDIR, "08_cluster_attribution_effects.tsv")
)

# Top cluster drivers per metric
top_cluster <- do.call(
  rbind,
  lapply(split(cluster_effect, cluster_effect$metric), function(d) {
    head(
      d[order(-d$abs_delta_contribution), , drop = FALSE],
      8
    )
  })
)

write_tsv(
  top_cluster,
  file.path(TABDIR, "09_TOP_cluster_drivers.tsv")
)

# ------------------------------------------------------------------------------
# FINAL macrophage-class attribution
# ------------------------------------------------------------------------------
class_rows <- list()
k <- 0L

for (s in SAMPLE_ORDER) {
  ix_s <- which(cell$sample == s)
  n_s <- length(ix_s)
  if (!n_s) next

  classes <- sort(unique(cell$final_class[ix_s]))

  for (cl in classes) {
    ix <- which(cell$sample == s & cell$final_class == cl)
    frac <- length(ix) / n_s

    for (pn in names(PANELS)) {
      mscore <- safe_mean(cell[[pn]][ix])

      k <- k + 1L
      class_rows[[k]] <- data.frame(
        sample = s,
        group = unique(cell$group[ix])[1],
        final_class = cl,
        n_cells = length(ix),
        sample_total_cells = n_s,
        class_fraction = frac,
        metric = pn,
        class_mean_score = mscore,
        weighted_contribution = frac * mscore,
        stringsAsFactors = FALSE
      )
    }
  }
}

class_sample <- do.call(rbind, class_rows)

write_tsv(
  class_sample,
  file.path(TABDIR, "10_finalclass_sample_mechanistic_contribution.tsv")
)

class_effect_rows <- list()
k <- 0L

for (pn in names(PANELS)) {
  d0 <- class_sample[class_sample$metric == pn, , drop = FALSE]

  for (cl in sort(unique(d0$final_class))) {
    d <- d0[d0$final_class == cl, , drop = FALSE]

    full <- data.frame(
      sample = SAMPLE_ORDER,
      group = c("Sham", "Sham", "Tx", "Tx"),
      stringsAsFactors = FALSE
    )

    full <- merge(
      full,
      d[, c(
        "sample",
        "class_fraction",
        "class_mean_score",
        "weighted_contribution"
      )],
      by = "sample",
      all.x = TRUE,
      sort = FALSE
    )

    full$class_fraction[is.na(full$class_fraction)] <- 0
    full$weighted_contribution[is.na(full$weighted_contribution)] <- 0

    sh <- full$weighted_contribution[full$group == "Sham"]
    tx <- full$weighted_contribution[full$group == "Tx"]

    k <- k + 1L
    class_effect_rows[[k]] <- data.frame(
      metric = pn,
      final_class = cl,
      sham_mean_weighted_contribution = safe_mean(sh),
      tx_mean_weighted_contribution = safe_mean(tx),
      tx_minus_sham_contribution = safe_mean(tx) - safe_mean(sh),
      abs_delta_contribution = abs(safe_mean(tx) - safe_mean(sh)),
      directional_consistency = direction_label(sh, tx),
      stringsAsFactors = FALSE
    )
  }
}

class_effect <- do.call(rbind, class_effect_rows)

class_effect <- class_effect[
  order(
    class_effect$metric,
    -class_effect$abs_delta_contribution
  ),
  ,
  drop = FALSE
]

write_tsv(
  class_effect,
  file.path(TABDIR, "11_finalclass_attribution_effects.tsv")
)

# ------------------------------------------------------------------------------
# Key combined interpretation table
# ------------------------------------------------------------------------------
KEY <- c(
  "IL10_receptor_signaling",
  "IL10_downstream_response",
  "IL10_axis_combined",
  "Anti_inflammatory_resolution",
  "Fibrosis_associated_Mphi",
  "HSC_activating_ligands",
  "ECM_remodeling",
  "MASH_pathogenic_inflammatory",
  "IFN_inflammatory",
  "Efferocytosis_resolution"
)

key_cell <- group_effect[
  match(KEY, group_effect$metric),
  c(
    "metric",
    "sham_mean",
    "tx_mean",
    "tx_minus_sham",
    "percent_change_vs_sham",
    "directional_consistency"
  ),
  drop = FALSE
]

key_pb <- pb_effect[
  match(KEY, pb_effect$metric),
  c(
    "metric",
    "sham_mean_CP10k",
    "tx_mean_CP10k",
    "tx_minus_sham_CP10k",
    "tx_over_sham",
    "directional_consistency"
  ),
  drop = FALSE
]

names(key_pb)[names(key_pb) == "directional_consistency"] <- (
  "pseudobulk_directional_consistency"
)

key <- merge(
  key_cell,
  key_pb,
  by = "metric",
  all = TRUE,
  sort = FALSE
)

key <- key[match(KEY, key$metric), , drop = FALSE]

write_tsv(
  key,
  file.path(TABDIR, "12_KEY_mechanistic_summary.tsv")
)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

plot_metric <- function(metric_name, filename, ylab = "Mean cell score") {
  d <- sample_score[
    sample_score$metric == metric_name,
    ,
    drop = FALSE
  ]

  d$group <- factor(d$group, levels = c("Sham", "Tx"))

  p <- ggplot(
    d,
    aes(
      x = group,
      y = mean_cell_score
    )
  ) +
    geom_point(
      size = 3,
      position = position_jitter(width = 0.06, height = 0)
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.45,
      linewidth = 0.5
    ) +
    geom_text(
      aes(label = sample),
      nudge_x = 0.12,
      size = 3,
      check_overlap = TRUE
    ) +
    labs(
      x = NULL,
      y = ylab,
      title = paste0(metric_name, " | sample-level"),
      subtitle = "Each point = one biological sample; descriptive only (n=2/group)"
    ) +
    theme_classic(base_size = 11)

  ggsave(
    file.path(FIGDIR, filename),
    p,
    width = 6,
    height = 5
  )
}

plot_metric(
  "IL10_axis_combined",
  "01_IL10_axis_combined_sample_level.pdf",
  "Mean IL10-axis-compatible score"
)

plot_metric(
  "IL10_downstream_response",
  "02_IL10_downstream_response_sample_level.pdf",
  "Mean downstream IL10-compatible response score"
)

plot_metric(
  "Anti_inflammatory_resolution",
  "03_Anti_inflammatory_resolution_sample_level.pdf",
  "Mean anti-inflammatory / resolution score"
)

plot_metric(
  "Fibrosis_associated_Mphi",
  "04_Fibrosis_associated_Mphi_sample_level.pdf",
  "Mean fibrosis-associated MΦ score"
)

plot_metric(
  "HSC_activating_ligands",
  "05_HSC_activating_ligands_sample_level.pdf",
  "Mean HSC-activating ligand score"
)

plot_metric(
  "MASH_pathogenic_inflammatory",
  "06_MASH_pathogenic_inflammatory_sample_level.pdf",
  "Mean MASH-pathogenic inflammatory score"
)

plot_metric(
  "Efferocytosis_resolution",
  "07_Efferocytosis_resolution_sample_level.pdf",
  "Mean efferocytosis / resolution score"
)

# Top cluster attribution for selected metrics
for (pn in c(
  "IL10_axis_combined",
  "Anti_inflammatory_resolution",
  "Fibrosis_associated_Mphi",
  "HSC_activating_ligands",
  "MASH_pathogenic_inflammatory",
  "Efferocytosis_resolution"
)) {
  d <- cluster_effect[
    cluster_effect$metric == pn,
    ,
    drop = FALSE
  ]

  d <- head(
    d[order(-d$abs_delta_contribution), , drop = FALSE],
    10
  )

  d$cluster <- factor(
    d$cluster,
    levels = rev(d$cluster)
  )

  p <- ggplot(
    d,
    aes(
      x = tx_minus_sham_contribution,
      y = cluster
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_col() +
    labs(
      x = "Tx - Sham weighted contribution",
      y = "Frozen cluster",
      title = paste0("Top cluster drivers | ", pn),
      subtitle = "Contribution = sample cluster fraction × cluster mean score"
    ) +
    theme_classic(base_size = 10)

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "cluster_driver_",
        pn,
        ".pdf"
      )
    ),
    p,
    width = 7,
    height = 5
  )
}

# FINAL class attribution heatmap-like tile plot
class_plot <- class_effect[
  class_effect$metric %in% c(
    "IL10_axis_combined",
    "Anti_inflammatory_resolution",
    "Fibrosis_associated_Mphi",
    "HSC_activating_ligands",
    "MASH_pathogenic_inflammatory",
    "Efferocytosis_resolution"
  ),
  ,
  drop = FALSE
]

p_class <- ggplot(
  class_plot,
  aes(
    x = metric,
    y = final_class,
    fill = tx_minus_sham_contribution
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf("%.3f", tx_minus_sham_contribution)
    ),
    size = 2.6
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Tx - Sham\ncontribution",
    title = "FINAL macrophage-class attribution",
    subtitle = "Weighted contribution = class fraction × class mean score"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(FIGDIR, "08_FINAL_class_attribution.pdf"),
  p_class,
  width = 11,
  height = 6
)

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

design <- data.frame(
  item = c(
    "input_RDS",
    "primary_comparison",
    "Sham_samples",
    "Tx_samples",
    "replicate_unit",
    "inferential_policy",
    "IL10_interpretation",
    "profibrotic_interpretation",
    "expression_source",
    "DecontX_policy",
    "cluster_attribution_definition",
    "finalclass_attribution_definition",
    "reclustering"
  ),
  value = c(
    RDS,
    "Sham vs Tx",
    "Sham1;Sham20",
    "Tx17;Tx5",
    "biological sample",
    "No formal p-values; n=2/group; sample-level effect sizes and directional consistency",
    "IL10-axis-compatible transcriptional activity; not proof of causal IL10 signaling",
    "Fibrosis-associated and HSC-activating panels are mechanistic signatures, not direct fibrosis measurements",
    "original filtered RNA counts in Primary-Clean MΦ object",
    "DecontX corrected counts reserved for sensitivity analysis",
    "sample cluster fraction × cluster mean score",
    "sample FINAL-class fraction × class mean score",
    "none"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  design,
  file.path(PROVDIR, "analysis_design_v1.0.tsv")
)

panel_table <- do.call(
  rbind,
  lapply(names(PANELS), function(pn) {
    data.frame(
      panel = pn,
      genes = paste(PANELS[[pn]], collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  panel_table,
  file.path(PROVDIR, "marker_panels_v1.0.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(PROVDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  01_panel_coverage.tsv")
message("  03_sample_mechanistic_scores.tsv")
message("  04_group_mechanistic_effects.tsv")
message("  06_group_pseudobulk_mechanistic_effects.tsv")
message("  09_TOP_cluster_drivers.tsv")
message("  11_finalclass_attribution_effects.tsv")
message("  12_KEY_mechanistic_summary.tsv")
message("Figures:")
message("  01_IL10_axis_combined_sample_level.pdf")
message("  03_Anti_inflammatory_resolution_sample_level.pdf")
message("  04_Fibrosis_associated_Mphi_sample_level.pdf")
message("  05_HSC_activating_ligands_sample_level.pdf")
message("  06_MASH_pathogenic_inflammatory_sample_level.pdf")
message("  07_Efferocytosis_resolution_sample_level.pdf")
message("  08_FINAL_class_attribution.pdf")
message("============================================================")
