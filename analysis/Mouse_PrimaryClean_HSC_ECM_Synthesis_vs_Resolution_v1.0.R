#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_HSC_ECM_Synthesis_vs_Resolution_v1.0.R
#
# PURPOSE
# -------
# Decompose the Tx-associated HSC program into:
#   1) collagen synthesis / maturation
#   2) ECM assembly
#   3) ECM degradation / turnover
#   4) TIMP / antiprotease restraint
#   5) myofibroblast contractility
#   6) quiescence / deactivation
#
# Primary-Clean HSC only.
# Biological replicate = SAMPLE (Sham1/Sham20 vs Tx17/Tx5).
# No cell-level inferential p-values.
#
# Key derived index:
#   Matrix turnover index =
#     mean(log1p CP10k of MMP/ADAMTS genes)
#     - mean(log1p CP10k of collagen synthesis/maturation genes)
#
# This is a descriptive transcriptional index, not a direct measurement
# of collagen degradation or extracellular protease activity.
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
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_HSC",
  "Mouse_PrimaryClean_HSC_ECM_Synthesis_vs_Resolution_v1.0"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, f) {
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

direction_rule <- function(sham, tx) {
  sham <- sham[is.finite(sham)]
  tx <- tx[is.finite(tx)]

  if (length(sham) < 2 || length(tx) < 2) return("insufficient")
  if (min(tx) > max(sham)) return("Tx_above_all_Sham")
  if (max(tx) < min(sham)) return("Tx_below_all_Sham")
  "overlap"
}

get_joined_counts <- function(obj) {
  lyr <- Layers(obj[["RNA"]])
  count_layers <- grep("^counts($|\\.)", lyr, value = TRUE)

  if (!length(count_layers)) stop("No RNA counts layer found.")

  if (length(count_layers) > 1 || !("counts" %in% count_layers)) {
    obj <- JoinLayers(obj, assay = "RNA")
  }

  mat <- LayerData(obj[["RNA"]], layer = "counts")

  if (!inherits(mat, "Matrix")) {
    mat <- Matrix(mat, sparse = TRUE)
  }

  mat
}

if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)

message("Loading Primary-Clean HSC...")
obj <- readRDS(HSC_RDS)

if (ncol(obj) != 12941) {
  stop("Unexpected Primary-Clean HSC cell count: ", ncol(obj))
}

md <- obj@meta.data

required <- c(
  "sample_hsc_v610",
  "HSC_state3_v612",
  "qc_primary_include"
)

missing <- setdiff(required, colnames(md))
if (length(missing)) {
  stop("Missing metadata: ", paste(missing, collapse = ", "))
}

if (!all(md$qc_primary_include %in% TRUE)) {
  stop("Primary-Clean HSC contains qc_primary_include != TRUE.")
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(
  Sham1 = "Sham",
  Sham20 = "Sham",
  Tx17 = "Tx",
  Tx5 = "Tx"
)

sample <- as.character(md$sample_hsc_v610)
state3 <- as.character(md$HSC_state3_v612)

keep <- sample %in% SAMPLES
obj <- subset(obj, cells = rownames(md)[keep])
md <- obj@meta.data

sample <- as.character(md$sample_hsc_v610)
state3 <- as.character(md$HSC_state3_v612)

message("Joining RNA counts...")
counts <- get_joined_counts(obj)
counts <- counts[, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# Mechanistic gene panels
# ------------------------------------------------------------------------------

PANELS <- list(
  Collagen_synthesis_core = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1"
  ),

  Collagen_maturation_crosslinking = c(
    "P4ha1", "P4ha2", "P4hb",
    "Plod1", "Plod2",
    "Serpinh1",
    "Pcolce", "Pcolce2",
    "Lox", "Loxl1", "Loxl2"
  ),

  ECM_assembly = c(
    "Fn1", "Sparc", "Dcn", "Lum", "Bgn", "Thbs1", "Thbs2"
  ),

  MMP_degradation = c(
    "Mmp2", "Mmp3", "Mmp9", "Mmp13", "Mmp14", "Mmp19"
  ),

  ADAMTS_degradation = c(
    "Adamts2", "Adamts4", "Adamts5", "Adamts12", "Adamts14"
  ),

  TIMP_antiprotease = c(
    "Timp1", "Timp2", "Timp3"
  ),

  Fibrogenic_stress = c(
    "Ctgf", "Serpine1", "Tgfb1", "Tgfbi"
  ),

  Myofibroblast_contractile = c(
    "Acta2", "Tagln", "Myl9", "Tpm2", "Myh11"
  ),

  Quiescence_deactivation = c(
    "Lrat", "Rbp1", "Reln", "Cygb", "Colec11", "Vipr1", "Pparg", "Bambi", "Smad7"
  )
)

coverage <- do.call(
  rbind,
  lapply(names(PANELS), function(pn) {
    g <- PANELS[[pn]]
    present <- intersect(g, rownames(counts))
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
  file.path(OUTDIR, "01_panel_coverage.tsv")
)

# ------------------------------------------------------------------------------
# Sample-level ALL-HSC pseudobulk CP10k
# ------------------------------------------------------------------------------

all_genes <- unique(unlist(PANELS))
all_genes <- intersect(all_genes, rownames(counts))

gene_rows <- list()
panel_rows <- list()
index_rows <- list()

kg <- 0L
kp <- 0L
ki <- 0L

for (s in SAMPLES) {

  ix <- which(sample == s)
  if (!length(ix)) next

  sums <- Matrix::rowSums(counts[, ix, drop = FALSE])
  lib <- sum(sums)

  cp <- as.numeric(sums) / lib * 1e4
  names(cp) <- rownames(counts)

  for (g in all_genes) {
    kg <- kg + 1L

    gene_rows[[kg]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      gene = g,
      cp10k = cp[g],
      log1p_cp10k = log1p(cp[g]),
      stringsAsFactors = FALSE
    )
  }

  panel_score <- list()

  for (pn in names(PANELS)) {
    gg <- intersect(PANELS[[pn]], rownames(counts))
    vals <- cp[gg]

    score <- mean(log1p(vals))
    panel_score[[pn]] <- score

    kp <- kp + 1L

    panel_rows[[kp]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      panel = pn,
      mean_log1p_CP10k = score,
      sum_CP10k = sum(vals),
      stringsAsFactors = FALSE
    )
  }

  degradation_genes <- unique(c(
    intersect(PANELS$MMP_degradation, rownames(counts)),
    intersect(PANELS$ADAMTS_degradation, rownames(counts))
  ))

  synthesis_genes <- unique(c(
    intersect(PANELS$Collagen_synthesis_core, rownames(counts)),
    intersect(PANELS$Collagen_maturation_crosslinking, rownames(counts))
  ))

  degradation_score <- mean(log1p(cp[degradation_genes]))
  synthesis_score <- mean(log1p(cp[synthesis_genes]))

  ki <- ki + 1L

  index_rows[[ki]] <- data.frame(
    sample = s,
    group = GROUP[[s]],
    degradation_score = degradation_score,
    collagen_synthesis_maturation_score = synthesis_score,
    matrix_turnover_index = degradation_score - synthesis_score,
    quiescence_minus_contractile =
      panel_score$Quiescence_deactivation -
      panel_score$Myofibroblast_contractile,
    stringsAsFactors = FALSE
  )
}

gene_pb <- do.call(rbind, gene_rows)
panel_pb <- do.call(rbind, panel_rows)
index_df <- do.call(rbind, index_rows)

write_tsv(
  gene_pb,
  file.path(OUTDIR, "02_ALLHSC_gene_pseudobulk.tsv")
)

write_tsv(
  panel_pb,
  file.path(OUTDIR, "03_ALLHSC_panel_pseudobulk.tsv")
)

write_tsv(
  index_df,
  file.path(OUTDIR, "04_matrix_turnover_indices_by_sample.tsv")
)

# ------------------------------------------------------------------------------
# ALL-HSC panel group effects
# ------------------------------------------------------------------------------

panel_effect_rows <- list()
k <- 0L

for (pn in unique(panel_pb$panel)) {

  d <- panel_pb[panel_pb$panel == pn, , drop = FALSE]

  sham <- d$mean_log1p_CP10k[d$group == "Sham"]
  tx <- d$mean_log1p_CP10k[d$group == "Tx"]

  sham_sum <- d$sum_CP10k[d$group == "Sham"]
  tx_sum <- d$sum_CP10k[d$group == "Tx"]

  k <- k + 1L

  panel_effect_rows[[k]] <- data.frame(
    panel = pn,
    Sham_mean_log_score = safe_mean(sham),
    Tx_mean_log_score = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    log_score_direction = direction_rule(sham, tx),
    Sham_mean_sum_CP10k = safe_mean(sham_sum),
    Tx_mean_sum_CP10k = safe_mean(tx_sum),
    sum_CP10k_ratio_Tx_over_Sham =
      (safe_mean(tx_sum) + 0.01) / (safe_mean(sham_sum) + 0.01),
    sum_CP10k_log2_ratio = log2(
      (safe_mean(tx_sum) + 0.01) /
        (safe_mean(sham_sum) + 0.01)
    ),
    sum_CP10k_direction = direction_rule(sham_sum, tx_sum),
    stringsAsFactors = FALSE
  )
}

panel_effects <- do.call(rbind, panel_effect_rows)

write_tsv(
  panel_effects,
  file.path(OUTDIR, "05_ALLHSC_panel_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Gene-level effects
# ------------------------------------------------------------------------------

gene_effect_rows <- list()
k <- 0L

for (g in unique(gene_pb$gene)) {

  d <- gene_pb[gene_pb$gene == g, , drop = FALSE]

  sham <- d$cp10k[d$group == "Sham"]
  tx <- d$cp10k[d$group == "Tx"]

  sm <- safe_mean(sham)
  tm <- safe_mean(tx)

  k <- k + 1L

  gene_effect_rows[[k]] <- data.frame(
    gene = g,
    Sham_mean_CP10k = sm,
    Tx_mean_CP10k = tm,
    ratio_Tx_over_Sham = (tm + 0.01) / (sm + 0.01),
    log2_ratio_Tx_over_Sham = log2((tm + 0.01) / (sm + 0.01)),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

gene_effects <- do.call(rbind, gene_effect_rows)

gene_effects <- gene_effects[
  order(-abs(gene_effects$log2_ratio_Tx_over_Sham)),
  ,
  drop = FALSE
]

write_tsv(
  gene_effects,
  file.path(OUTDIR, "06_ALLHSC_gene_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Derived index group effects
# ------------------------------------------------------------------------------

index_effect_rows <- list()
k <- 0L

for (metric in c(
  "degradation_score",
  "collagen_synthesis_maturation_score",
  "matrix_turnover_index",
  "quiescence_minus_contractile"
)) {

  sham <- index_df[[metric]][index_df$group == "Sham"]
  tx <- index_df[[metric]][index_df$group == "Tx"]

  k <- k + 1L

  index_effect_rows[[k]] <- data.frame(
    metric = metric,
    Sham_mean = safe_mean(sham),
    Tx_mean = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

index_effects <- do.call(rbind, index_effect_rows)

write_tsv(
  index_effects,
  file.path(OUTDIR, "07_derived_index_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# State-specific panel scores
# ------------------------------------------------------------------------------

state_panel_rows <- list()
k <- 0L

for (s in SAMPLES) {

  for (st in sort(unique(state3[sample == s]))) {

    ix <- which(sample == s & state3 == st)
    if (!length(ix)) next

    sums <- Matrix::rowSums(counts[, ix, drop = FALSE])
    lib <- sum(sums)

    cp <- as.numeric(sums) / lib * 1e4
    names(cp) <- rownames(counts)

    for (pn in names(PANELS)) {
      gg <- intersect(PANELS[[pn]], rownames(counts))
      vals <- cp[gg]

      k <- k + 1L

      state_panel_rows[[k]] <- data.frame(
        sample = s,
        group = GROUP[[s]],
        HSC_state3 = st,
        panel = pn,
        mean_log1p_CP10k = mean(log1p(vals)),
        sum_CP10k = sum(vals),
        stringsAsFactors = FALSE
      )
    }
  }
}

state_panel <- do.call(rbind, state_panel_rows)

write_tsv(
  state_panel,
  file.path(OUTDIR, "08_state_panel_pseudobulk.tsv")
)

# State panel group effects
keys <- unique(state_panel[, c("HSC_state3", "panel")])

state_effect_rows <- list()
k <- 0L

for (i in seq_len(nrow(keys))) {

  st <- keys$HSC_state3[i]
  pn <- keys$panel[i]

  d <- state_panel[
    state_panel$HSC_state3 == st &
      state_panel$panel == pn,
    ,
    drop = FALSE
  ]

  sham <- d$mean_log1p_CP10k[d$group == "Sham"]
  tx <- d$mean_log1p_CP10k[d$group == "Tx"]

  k <- k + 1L

  state_effect_rows[[k]] <- data.frame(
    HSC_state3 = st,
    panel = pn,
    Sham_mean_log_score = safe_mean(sham),
    Tx_mean_log_score = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

state_effects <- do.call(rbind, state_effect_rows)

write_tsv(
  state_effects,
  file.path(OUTDIR, "09_state_panel_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Compact mechanistic summary
# ------------------------------------------------------------------------------

compact <- panel_effects[
  panel_effects$panel %in% c(
    "Collagen_synthesis_core",
    "Collagen_maturation_crosslinking",
    "ECM_assembly",
    "MMP_degradation",
    "ADAMTS_degradation",
    "TIMP_antiprotease",
    "Fibrogenic_stress",
    "Myofibroblast_contractile",
    "Quiescence_deactivation"
  ),
  ,
  drop = FALSE
]

write_tsv(
  compact,
  file.path(OUTDIR, "10_KEY_ECM_synthesis_vs_resolution_summary.tsv")
)

# ------------------------------------------------------------------------------
# Integrity
# ------------------------------------------------------------------------------

checks <- data.frame(
  item = c(
    "PrimaryClean_HSC_cells",
    "all_qc_primary_include_TRUE",
    "all_four_samples_present",
    "biological_replicate",
    "cell_level_inferential_pvalues_used",
    "matrix_turnover_index_definition"
  ),
  value = c(
    ncol(obj),
    as.character(all(obj$qc_primary_include %in% TRUE)),
    as.character(all(SAMPLES %in% unique(sample))),
    "sample",
    "FALSE",
    "mean log1p(CP10k) MMP+ADAMTS minus mean log1p(CP10k) collagen synthesis+maturation"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  checks,
  file.path(OUTDIR, "11_analysis_integrity_summary.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  05_ALLHSC_panel_group_effects.tsv")
message("  06_ALLHSC_gene_group_effects.tsv")
message("  07_derived_index_group_effects.tsv")
message("  09_state_panel_group_effects.tsv")
message("  10_KEY_ECM_synthesis_vs_resolution_summary.tsv")
message("  11_analysis_integrity_summary.tsv")
message("============================================================")
