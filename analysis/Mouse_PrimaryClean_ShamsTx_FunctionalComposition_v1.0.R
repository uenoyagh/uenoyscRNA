#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_ShamsTx_FunctionalComposition_v1.0.R
#
# PURPOSE
# -------
# First primary downstream analysis using the finalized Mouse MΦ Primary-Clean
# RDS:
#
#   A. Sham vs Tx macrophage functional axes
#      - M1 core
#      - IFN / inflammatory context
#      - M2 core
#      - efferocytosis / resolution
#   B. sample-level pseudobulk M2/M1 ratio
#   C. Resident / Recruited / LAM identity-axis scores
#   D. composition of the frozen FINAL macrophage classes
#
# STATISTICAL POLICY
# ------------------
# Sham has 2 biological samples (Sham1, Sham20).
# Tx has 2 biological samples (Tx17, Tx5).
#
# Therefore:
# - CELLS ARE NOT TREATED AS INDEPENDENT REPLICATES.
# - NO cell-level p-values are reported.
# - NO formal inferential p-values are reported in v1.0 because n=2/group.
# - Primary outputs are sample-level estimates, group means, and effect sizes.
# - Pseudobulk counts are constructed independently for each sample.
#
# M1 / M2 POLICY
# --------------
# M1 and M2 are continuous functional axes, NOT macrophage identities.
#
# PRIMARY EXPRESSION POLICY
# -------------------------
# Original filtered RNA counts in Mouse_MASH_Mphi_PrimaryClean_v5.0.rds are
# primary. DecontX-corrected counts remain sensitivity-analysis data and are
# not substituted here.
#
# This script does NOT modify the input RDS and does NOT recluster.
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
  "Mouse_PrimaryClean_ShamsTx_FunctionalComposition_v1.0"
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
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

# ------------------------------------------------------------------------------
# Load finalized Primary-Clean MΦ object
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

sample <- as.character(md$sample_display_FIXED2)
group <- ifelse(
  sample %in% c("Sham1", "Sham20"),
  "Sham",
  ifelse(sample %in% c("Tx17", "Tx5"), "Tx", "Other")
)

use_cells <- rownames(md)[group %in% c("Sham", "Tx")]

if (!length(use_cells)) stop("No Sham / Tx cells found.")

obj_st <- subset(obj, cells = use_cells)
md_st <- obj_st@meta.data

md_st$analysis_group <- ifelse(
  as.character(md_st$sample_display_FIXED2) %in% c("Sham1", "Sham20"),
  "Sham",
  "Tx"
)

md_st$analysis_group <- factor(
  md_st$analysis_group,
  levels = c("Sham", "Tx")
)

SAMPLE_ORDER <- c("Sham1", "Sham20", "Tx17", "Tx5")

md_st$sample_display_FIXED2 <- factor(
  as.character(md_st$sample_display_FIXED2),
  levels = SAMPLE_ORDER
)

# ------------------------------------------------------------------------------
# Frozen marker panels
# ------------------------------------------------------------------------------

PANELS <- list(
  M1_primary = c("Il1b", "Tnf", "Ccl3", "Ccl4", "Cxcl2"),

  M1_supportive = c(
    "Nfkbia", "Tnfaip3", "Nfkb2", "Nlrp3", "Relb", "Gpr84"
  ),

  M1_IFN = c(
    "Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"
  ),

  M1_MASH_pathogenic = c(
    "F13a1", "Olr1", "Selplg"
  ),

  M1_context = c(
    "Ccl2", "Irf5", "Ptgs2", "Cd86", "Cxcl3"
  ),

  M2_primary = c(
    "Mrc1", "Cd163", "Hmox1", "Maf"
  ),

  M2_supportive = c(
    "Dab2", "Msr1", "Stab1", "Il1rn"
  ),

  M2_resolution_efferocytosis = c(
    "Mertk", "Axl", "Gas6", "Mfge8"
  ),

  M2_regulatory = c(
    "Klf4", "Pparg", "Itga4"
  ),

  M2_context = c(
    "Igf1", "Anxa1"
  ),

  Resident_Kupffer = c(
    "Marco", "Vsig4", "Folr2", "Timd4",
    "Cd5l", "C1qa", "C1qb", "C1qc"
  ),

  Recruited_Monocyte = c(
    "Ccr2", "S100a8", "S100a9", "Vcan", "Fcna", "Sell"
  ),

  LAM = c(
    "Trem2", "Gpnmb", "Cd9", "Apoe", "Lpl", "Fabp5", "Spp1"
  )
)

PANELS$M1_core <- unique(c(
  PANELS$M1_primary,
  PANELS$M1_supportive
))

PANELS$M2_core <- unique(c(
  PANELS$M2_primary,
  PANELS$M2_supportive
))

# ------------------------------------------------------------------------------
# Extract primary RNA counts
# ------------------------------------------------------------------------------

counts <- LayerData(
  obj_st[["RNA"]],
  layer = "counts"
)

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

# log1p(CP10k), sparse-safe
cp10k <- t(t(counts) / libsize * 1e4)
logcp10k <- log1p(cp10k)

# ------------------------------------------------------------------------------
# Cell-level functional scores
#
# Score = mean per-gene log1p(CP10k) within each pre-defined panel.
# These are functional axes only.
# ------------------------------------------------------------------------------

score_panel <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[g, , drop = FALSE]))
}

cell_scores <- data.frame(
  cell = colnames(obj_st),
  sample = as.character(md_st$sample_display_FIXED2),
  group = as.character(md_st$analysis_group),
  cluster = as.character(md_st$cluster_res2),
  final_class = as.character(md_st$macrophage_class_Res2_FINAL_v4145),
  stringsAsFactors = FALSE
)

for (pn in names(PANELS)) {
  cell_scores[[pn]] <- score_panel(logcp10k, PANELS[[pn]])
}

# A difference on the same logCP10k scale is the preferred continuous balance.
cell_scores$M2_minus_M1_core <- (
  cell_scores$M2_core - cell_scores$M1_core
)

write_tsv(
  cell_scores,
  file.path(TABDIR, "01_cell_functional_scores.tsv")
)

# ------------------------------------------------------------------------------
# Sample-level functional summaries
# ------------------------------------------------------------------------------

SCORE_COLS <- c(
  "M1_primary",
  "M1_supportive",
  "M1_core",
  "M1_IFN",
  "M1_MASH_pathogenic",
  "M1_context",
  "M2_primary",
  "M2_supportive",
  "M2_core",
  "M2_resolution_efferocytosis",
  "M2_regulatory",
  "M2_context",
  "M2_minus_M1_core",
  "Resident_Kupffer",
  "Recruited_Monocyte",
  "LAM"
)

sample_rows <- list()
kk <- 0L

for (s in SAMPLE_ORDER) {
  ix <- which(cell_scores$sample == s)
  if (!length(ix)) next

  for (sc in SCORE_COLS) {
    kk <- kk + 1L
    x <- cell_scores[[sc]][ix]

    sample_rows[[kk]] <- data.frame(
      sample = s,
      group = unique(cell_scores$group[ix])[1],
      n_cells = length(ix),
      metric = sc,
      mean_cell_score = safe_mean(x),
      median_cell_score = median(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

sample_function <- do.call(rbind, sample_rows)

write_tsv(
  sample_function,
  file.path(TABDIR, "02_sample_functional_scores.tsv")
)

# Group-level descriptive effect: sample is the replicate
effect_rows <- list()
kk <- 0L

for (sc in SCORE_COLS) {
  x <- sample_function[sample_function$metric == sc, , drop = FALSE]

  sh <- x$mean_cell_score[x$group == "Sham"]
  tx <- x$mean_cell_score[x$group == "Tx"]

  kk <- kk + 1L
  effect_rows[[kk]] <- data.frame(
    metric = sc,
    n_sham_samples = length(sh),
    n_tx_samples = length(tx),
    sham_mean_of_sample_means = safe_mean(sh),
    sham_sd = safe_sd(sh),
    tx_mean_of_sample_means = safe_mean(tx),
    tx_sd = safe_sd(tx),
    tx_minus_sham = safe_mean(tx) - safe_mean(sh),
    direction = ifelse(
      safe_mean(tx) > safe_mean(sh),
      "Tx_higher",
      ifelse(
        safe_mean(tx) < safe_mean(sh),
        "Tx_lower",
        "No_change"
      )
    ),
    inferential_p_value = NA_real_,
    statistical_note = "Descriptive sample-level effect only; n=2/group; cells are not independent replicates",
    stringsAsFactors = FALSE
  )
}

functional_effects <- do.call(rbind, effect_rows)

write_tsv(
  functional_effects,
  file.path(TABDIR, "03_group_functional_effects.tsv")
)

# ------------------------------------------------------------------------------
# Sample pseudobulk and M2/M1 ratio
#
# Ratio is computed from sample-level pseudobulk CP10k, NOT from cell ratios.
# Panel expression = mean linear CP10k across panel genes.
# pseudocount = 0.1 CP10k only for ratio stability.
# ------------------------------------------------------------------------------

PB_PSEUDOCOUNT <- 0.1

pb_rows <- list()
kk <- 0L

for (s in SAMPLE_ORDER) {
  cells_s <- colnames(obj_st)[
    as.character(md_st$sample_display_FIXED2) == s
  ]

  if (!length(cells_s)) next

  m <- counts[, cells_s, drop = FALSE]
  pb <- Matrix::rowSums(m)
  pb_cp10k <- pb / sum(pb) * 1e4

  panel_linear <- function(genes) {
    g <- intersect(genes, names(pb_cp10k))
    if (!length(g)) return(NA_real_)
    mean(pb_cp10k[g])
  }

  m1 <- panel_linear(PANELS$M1_core)
  m2 <- panel_linear(PANELS$M2_core)

  kk <- kk + 1L
  pb_rows[[kk]] <- data.frame(
    sample = s,
    group = ifelse(s %in% c("Sham1", "Sham20"), "Sham", "Tx"),
    n_cells = length(cells_s),
    M1_core_mean_CP10k = m1,
    M2_core_mean_CP10k = m2,
    M2_to_M1_ratio = (m2 + PB_PSEUDOCOUNT) / (m1 + PB_PSEUDOCOUNT),
    log2_M2_to_M1_ratio = log2(
      (m2 + PB_PSEUDOCOUNT) / (m1 + PB_PSEUDOCOUNT)
    ),
    stringsAsFactors = FALSE
  )
}

pb_ratio <- do.call(rbind, pb_rows)

write_tsv(
  pb_ratio,
  file.path(TABDIR, "04_sample_pseudobulk_M2_M1_ratio.tsv")
)

sh_ratio <- pb_ratio$M2_to_M1_ratio[pb_ratio$group == "Sham"]
tx_ratio <- pb_ratio$M2_to_M1_ratio[pb_ratio$group == "Tx"]

pb_ratio_effect <- data.frame(
  metric = "sample_pseudobulk_M2_to_M1_ratio",
  n_sham_samples = length(sh_ratio),
  n_tx_samples = length(tx_ratio),
  sham_mean = safe_mean(sh_ratio),
  sham_sd = safe_sd(sh_ratio),
  tx_mean = safe_mean(tx_ratio),
  tx_sd = safe_sd(tx_ratio),
  tx_minus_sham = safe_mean(tx_ratio) - safe_mean(sh_ratio),
  tx_over_sham = safe_mean(tx_ratio) / safe_mean(sh_ratio),
  inferential_p_value = NA_real_,
  statistical_note = "Descriptive pseudobulk sample-level effect only; n=2/group",
  stringsAsFactors = FALSE
)

write_tsv(
  pb_ratio_effect,
  file.path(TABDIR, "05_group_pseudobulk_M2_M1_effect.tsv")
)

# ------------------------------------------------------------------------------
# Frozen FINAL macrophage-class composition
# ------------------------------------------------------------------------------

comp <- as.data.frame(
  table(
    sample = as.character(md_st$sample_display_FIXED2),
    final_class = as.character(md_st$macrophage_class_Res2_FINAL_v4145)
  ),
  stringsAsFactors = FALSE
)

comp <- comp[comp$sample %in% SAMPLE_ORDER, , drop = FALSE]

sample_totals <- aggregate(
  Freq ~ sample,
  data = comp,
  FUN = sum
)

names(sample_totals)[2] <- "sample_total"

comp <- merge(
  comp,
  sample_totals,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

comp$group <- ifelse(
  comp$sample %in% c("Sham1", "Sham20"),
  "Sham",
  "Tx"
)

comp$fraction <- comp$Freq / comp$sample_total

write_tsv(
  comp,
  file.path(TABDIR, "06_sample_finalclass_composition.tsv")
)

classes <- sort(unique(comp$final_class))

comp_effect_rows <- list()
kk <- 0L

for (cl in classes) {
  x <- comp[comp$final_class == cl, , drop = FALSE]
  sh <- x$fraction[x$group == "Sham"]
  tx <- x$fraction[x$group == "Tx"]

  kk <- kk + 1L
  comp_effect_rows[[kk]] <- data.frame(
    final_class = cl,
    sham_mean_fraction = safe_mean(sh),
    sham_sd = safe_sd(sh),
    tx_mean_fraction = safe_mean(tx),
    tx_sd = safe_sd(tx),
    tx_minus_sham_fraction = safe_mean(tx) - safe_mean(sh),
    tx_minus_sham_percentage_points = 100 * (
      safe_mean(tx) - safe_mean(sh)
    ),
    direction = ifelse(
      safe_mean(tx) > safe_mean(sh),
      "Tx_higher",
      ifelse(
        safe_mean(tx) < safe_mean(sh),
        "Tx_lower",
        "No_change"
      )
    ),
    inferential_p_value = NA_real_,
    statistical_note = "Descriptive sample-level composition only; n=2/group",
    stringsAsFactors = FALSE
  )
}

comp_effect <- do.call(rbind, comp_effect_rows)

write_tsv(
  comp_effect,
  file.path(TABDIR, "07_group_finalclass_composition_effect.tsv")
)

# ------------------------------------------------------------------------------
# Resident / Recruited / LAM continuous identity-axis summaries
# ------------------------------------------------------------------------------

identity_metrics <- c(
  "Resident_Kupffer",
  "Recruited_Monocyte",
  "LAM"
)

identity_sample <- sample_function[
  sample_function$metric %in% identity_metrics,
  ,
  drop = FALSE
]

write_tsv(
  identity_sample,
  file.path(TABDIR, "08_sample_identity_axis_scores.tsv")
)

identity_effect <- functional_effects[
  functional_effects$metric %in% identity_metrics,
  ,
  drop = FALSE
]

write_tsv(
  identity_effect,
  file.path(TABDIR, "09_group_identity_axis_effect.tsv")
)

# ------------------------------------------------------------------------------
# Cluster × group descriptive functional summaries
# ------------------------------------------------------------------------------

cluster_metrics <- c(
  "M1_core",
  "M1_IFN",
  "M1_MASH_pathogenic",
  "M2_core",
  "M2_resolution_efferocytosis",
  "M2_minus_M1_core",
  "Resident_Kupffer",
  "Recruited_Monocyte",
  "LAM"
)

cluster_rows <- list()
kk <- 0L

for (cl in sort(unique(cell_scores$cluster))) {
  for (gp in c("Sham", "Tx")) {
    ix <- which(
      cell_scores$cluster == cl &
      cell_scores$group == gp
    )

    if (!length(ix)) next

    for (sc in cluster_metrics) {
      kk <- kk + 1L
      cluster_rows[[kk]] <- data.frame(
        cluster = cl,
        group = gp,
        n_cells = length(ix),
        metric = sc,
        mean_cell_score = safe_mean(cell_scores[[sc]][ix]),
        median_cell_score = median(
          cell_scores[[sc]][ix],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}

cluster_group_scores <- do.call(rbind, cluster_rows)

write_tsv(
  cluster_group_scores,
  file.path(TABDIR, "10_cluster_by_group_functional_scores.tsv")
)

# Cluster × sample counts/fractions for audit
cluster_comp <- as.data.frame(
  table(
    sample = cell_scores$sample,
    cluster = cell_scores$cluster
  ),
  stringsAsFactors = FALSE
)

cluster_comp <- cluster_comp[
  cluster_comp$sample %in% SAMPLE_ORDER,
  ,
  drop = FALSE
]

ctot <- aggregate(
  Freq ~ sample,
  data = cluster_comp,
  FUN = sum
)

names(ctot)[2] <- "sample_total"

cluster_comp <- merge(
  cluster_comp,
  ctot,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

cluster_comp$fraction <- (
  cluster_comp$Freq / cluster_comp$sample_total
)

write_tsv(
  cluster_comp,
  file.path(TABDIR, "11_cluster_by_sample_composition.tsv")
)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

plot_sample_metric <- function(metric_name, ylab, filename) {
  d <- sample_function[
    sample_function$metric == metric_name,
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
      subtitle = "Each point = one biological sample; no cell-level inferential testing"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(FIGDIR, filename),
    p,
    width = 6,
    height = 5
  )
}

plot_sample_metric(
  "M1_core",
  "Mean cell M1-core score\n(mean log1p CP10k)",
  "01_M1_core_sample_level.pdf"
)

plot_sample_metric(
  "M2_core",
  "Mean cell M2-core score\n(mean log1p CP10k)",
  "02_M2_core_sample_level.pdf"
)

plot_sample_metric(
  "M2_minus_M1_core",
  "M2-core minus M1-core\n(logCP10k score difference)",
  "03_M2_minus_M1_balance_sample_level.pdf"
)

plot_sample_metric(
  "M1_IFN",
  "Mean cell IFN/inflammatory score",
  "04_M1_IFN_sample_level.pdf"
)

plot_sample_metric(
  "M2_resolution_efferocytosis",
  "Mean cell resolution/efferocytosis score",
  "05_M2_resolution_efferocytosis_sample_level.pdf"
)

# Pseudobulk ratio
pb_ratio$group <- factor(
  pb_ratio$group,
  levels = c("Sham", "Tx")
)

p_ratio <- ggplot(
  pb_ratio,
  aes(
    x = group,
    y = M2_to_M1_ratio
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2
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
    y = "Sample pseudobulk M2/M1 ratio",
    title = "M2/M1 balance | sample pseudobulk",
    subtitle = "M2-core / M1-core mean CP10k; each point = one biological sample"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(FIGDIR, "06_pseudobulk_M2_M1_ratio.pdf"),
  p_ratio,
  width = 6,
  height = 5
)

# Frozen class composition
comp$sample <- factor(
  comp$sample,
  levels = SAMPLE_ORDER
)

p_comp <- ggplot(
  comp,
  aes(
    x = sample,
    y = fraction,
    fill = final_class
  )
) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL,
    y = "Fraction of Primary-Clean MΦ",
    fill = "FINAL MΦ class",
    title = "Frozen FINAL macrophage-class composition",
    subtitle = "Primary-Clean Mouse MΦ | Sham vs Tx"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(FIGDIR, "07_finalclass_composition_by_sample.pdf"),
  p_comp,
  width = 8.5,
  height = 5.5
)

# Identity axes
id_plot <- identity_sample
id_plot$group <- factor(
  id_plot$group,
  levels = c("Sham", "Tx")
)

p_id <- ggplot(
  id_plot,
  aes(
    x = group,
    y = mean_cell_score
  )
) +
  geom_point(
    size = 2.5,
    position = position_jitter(width = 0.06, height = 0)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.45
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = "Mean cell score",
    title = "Resident / Recruited / LAM identity axes",
    subtitle = "Continuous identity scores; not forced cell identities"
  ) +
  theme_classic(base_size = 10)

ggsave(
  file.path(FIGDIR, "08_identity_axes_sample_level.pdf"),
  p_id,
  width = 10,
  height = 4.8
)

# ------------------------------------------------------------------------------
# Analysis design / provenance
# ------------------------------------------------------------------------------

design <- data.frame(
  item = c(
    "input_RDS",
    "analysis",
    "primary_comparison",
    "Sham_samples",
    "Tx_samples",
    "replicate_unit",
    "inferential_policy",
    "M1_M2_policy",
    "M2_M1_primary_balance",
    "M2_M1_ratio_definition",
    "composition_annotation",
    "identity_axis_policy",
    "expression_source",
    "DecontX_policy",
    "reclustering"
  ),
  value = c(
    RDS,
    "Mouse Primary-Clean MΦ Sham vs Tx functional/composition analysis",
    "Sham vs Tx",
    "Sham1;Sham20",
    "Tx17;Tx5",
    "biological sample",
    "No formal p-values in v1.0 because n=2/group; descriptive sample-level effects only",
    "M1 and M2 are continuous functional axes, not identities",
    "M2_core minus M1_core on cell log1p(CP10k) scale",
    "sample pseudobulk mean M2-core CP10k / mean M1-core CP10k; pseudocount 0.1",
    "macrophage_class_Res2_FINAL_v4145",
    "Resident/Recruited/LAM are continuous identity axes; no forced reassignment",
    "original RNA counts from finalized Primary-Clean object",
    "DecontX-corrected counts reserved for later sensitivity analysis",
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
message("Primary review tables:")
message("  02_sample_functional_scores.tsv")
message("  03_group_functional_effects.tsv")
message("  04_sample_pseudobulk_M2_M1_ratio.tsv")
message("  05_group_pseudobulk_M2_M1_effect.tsv")
message("  06_sample_finalclass_composition.tsv")
message("  07_group_finalclass_composition_effect.tsv")
message("  08_sample_identity_axis_scores.tsv")
message("  09_group_identity_axis_effect.tsv")
message("Primary review figures:")
message("  03_M2_minus_M1_balance_sample_level.pdf")
message("  06_pseudobulk_M2_M1_ratio.pdf")
message("  07_finalclass_composition_by_sample.pdf")
message("  08_identity_axes_sample_level.pdf")
message("============================================================")
