#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_HSC_ShamsTx_FunctionalResponse_v1.0.R
#
# PURPOSE
# -------
# Receiver-side validation after MΦ -> HSC LR compatibility analysis.
#
# Primary-Clean HSC only.
# Contrast: Sham1/Sham20 vs Tx17/Tx5.
#
# Biological replicate = SAMPLE (n=2/group).
# No cell-level inferential p-values.
#
# Outputs:
#   - sample-level HSC state composition
#   - qHSC / aHSC ratio
#   - sample pseudobulk expression of HSC functional panels
#   - gene-level ECM/fibrosis output
#   - receptor-side changes relevant to MΦ -> HSC signaling
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
  "Mouse_PrimaryClean_HSC_ShamsTx_FunctionalResponse_v1.0"
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

if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)

message("Loading Primary-Clean HSC...")
obj <- readRDS(HSC_RDS)
md <- obj@meta.data

required <- c(
  "sample_hsc_v610",
  "HSC_state3_v612",
  "HSC_state2_v612",
  "qc_primary_include"
)

missing <- setdiff(required, colnames(md))
if (length(missing)) stop("Missing metadata: ", paste(missing, collapse = ", "))

if (ncol(obj) != 12941) {
  stop("Unexpected Primary-Clean HSC cell count: ", ncol(obj))
}

if (!all(md$qc_primary_include %in% TRUE)) {
  stop("Primary-Clean HSC contains qc_primary_include != TRUE.")
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(Sham1 = "Sham", Sham20 = "Sham", Tx17 = "Tx", Tx5 = "Tx")

sample <- as.character(md$sample_hsc_v610)
state3 <- as.character(md$HSC_state3_v612)
state2 <- as.character(md$HSC_state2_v612)

keep <- sample %in% SAMPLES
obj <- subset(obj, cells = rownames(md)[keep])
md <- obj@meta.data

sample <- as.character(md$sample_hsc_v610)
state3 <- as.character(md$HSC_state3_v612)
state2 <- as.character(md$HSC_state2_v612)

# ------------------------------------------------------------------------------
# HSC state composition
# ------------------------------------------------------------------------------

tab3 <- as.data.frame(
  table(sample = sample, HSC_state3 = state3),
  stringsAsFactors = FALSE
)

tot <- aggregate(Freq ~ sample, data = tab3, FUN = sum)
names(tot)[2] <- "sample_total"
tab3 <- merge(tab3, tot, by = "sample", all.x = TRUE, sort = FALSE)
tab3$fraction <- tab3$Freq / tab3$sample_total
tab3$group <- GROUP[tab3$sample]

write_tsv(tab3, file.path(OUTDIR, "01_HSC3state_sample_composition.tsv"))

tab2 <- as.data.frame(
  table(sample = sample, HSC_state2 = state2),
  stringsAsFactors = FALSE
)

tot2 <- aggregate(Freq ~ sample, data = tab2, FUN = sum)
names(tot2)[2] <- "sample_total"
tab2 <- merge(tab2, tot2, by = "sample", all.x = TRUE, sort = FALSE)
tab2$fraction <- tab2$Freq / tab2$sample_total
tab2$group <- GROUP[tab2$sample]

write_tsv(tab2, file.path(OUTDIR, "02_HSC2state_sample_composition.tsv"))

ratio_rows <- list()
k <- 0L

for (s in SAMPLES) {
  d <- tab2[tab2$sample == s, , drop = FALSE]
  q <- d$Freq[d$HSC_state2 == "qHSC"]
  a <- d$Freq[d$HSC_state2 == "aHSC"]

  if (!length(q)) q <- 0
  if (!length(a)) a <- 0

  k <- k + 1L
  ratio_rows[[k]] <- data.frame(
    sample = s,
    group = GROUP[[s]],
    qHSC_n = q,
    aHSC_n = a,
    qHSC_fraction = q / (q + a),
    aHSC_fraction = a / (q + a),
    qHSC_over_aHSC = (q + 0.5) / (a + 0.5),
    stringsAsFactors = FALSE
  )
}

ratio_df <- do.call(rbind, ratio_rows)
write_tsv(ratio_df, file.path(OUTDIR, "03_qHSC_aHSC_sample_ratio.tsv"))

# Group-level state effects
state_effect_rows <- list()
k <- 0L

for (st in unique(tab3$HSC_state3)) {
  d <- tab3[tab3$HSC_state3 == st, , drop = FALSE]
  sham <- d$fraction[d$group == "Sham"]
  tx <- d$fraction[d$group == "Tx"]

  k <- k + 1L
  state_effect_rows[[k]] <- data.frame(
    metric = paste0("fraction_", st),
    Sham_mean = safe_mean(sham),
    Tx_mean = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

for (metric in c("qHSC_fraction", "aHSC_fraction", "qHSC_over_aHSC")) {
  sham <- ratio_df[[metric]][ratio_df$group == "Sham"]
  tx <- ratio_df[[metric]][ratio_df$group == "Tx"]

  k <- k + 1L
  state_effect_rows[[k]] <- data.frame(
    metric = metric,
    Sham_mean = safe_mean(sham),
    Tx_mean = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

state_effects <- do.call(rbind, state_effect_rows)
write_tsv(state_effects, file.path(OUTDIR, "04_HSC_state_group_effects.tsv"))

# ------------------------------------------------------------------------------
# RNA counts
# ------------------------------------------------------------------------------

get_joined_counts <- function(obj) {
  lyr <- Layers(obj[["RNA"]])
  count_layers <- grep("^counts($|\\.)", lyr, value = TRUE)

  if (!length(count_layers)) stop("No RNA counts layer found.")

  if (length(count_layers) > 1 || !("counts" %in% count_layers)) {
    obj <- JoinLayers(obj, assay = "RNA")
  }

  mat <- LayerData(obj[["RNA"]], layer = "counts")
  if (!inherits(mat, "Matrix")) mat <- Matrix(mat, sparse = TRUE)
  mat
}

message("Joining RNA counts...")
counts <- get_joined_counts(obj)
counts <- counts[, colnames(obj), drop = FALSE]

# ------------------------------------------------------------------------------
# Functional panels
# ------------------------------------------------------------------------------

PANELS <- list(
  Quiescence_identity = c(
    "Lrat", "Rbp1", "Reln", "Cygb", "Colec11", "Vipr1"
  ),

  ECM_synthesis = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1",
    "Fn1", "Lum", "Dcn", "Bgn"
  ),

  Activated_fibrogenic_output = c(
    "Col1a1", "Col1a2", "Col3a1", "Fn1",
    "Timp1", "Ctgf", "Serpine1"
  ),

  Contractile = c(
    "Acta2", "Tagln", "Myl9", "Tpm2", "Myh11"
  ),

  TGF_receptors = c(
    "Tgfbr1", "Tgfbr2"
  ),

  PDGF_receptors = c(
    "Pdgfra", "Pdgfrb"
  ),

  SPP1_receptors = c(
    "Cd44", "Itgav", "Itgb1", "Itgb3"
  ),

  FN1_integrin_receptors = c(
    "Itga5", "Itgb1"
  ),

  EGFR_repair_axis = c(
    "Egfr"
  ),

  IL10_receptors = c(
    "Il10ra", "Il10rb"
  ),

  HGF_receptor = c(
    "Met"
  ),

  BMP7_receptors = c(
    "Bmpr1a", "Bmpr2"
  ),

  GAS6_receptors = c(
    "Axl", "Mertk", "Tyro3"
  ),

  APOE_receptor = c(
    "Lrp1"
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

write_tsv(coverage, file.path(OUTDIR, "05_panel_coverage.tsv"))

# ------------------------------------------------------------------------------
# Sample pseudobulk by state
# ------------------------------------------------------------------------------

pb_cp10k <- function(mat, groups) {
  lev <- unique(groups)

  res <- lapply(lev, function(g) {
    ix <- which(groups == g)
    s <- Matrix::rowSums(mat[, ix, drop = FALSE])
    lib <- sum(s)
    cp <- as.numeric(s) / lib * 1e4
    names(cp) <- rownames(mat)
    cp
  })

  out <- do.call(cbind, res)
  rownames(out) <- rownames(mat)
  colnames(out) <- lev
  out
}

pb_group <- paste(sample, state3, sep = "|||")
pb <- pb_cp10k(counts, pb_group)

# Gene-level table for all panel genes
genes <- unique(unlist(PANELS))
genes <- intersect(genes, rownames(pb))

gene_rows <- list()
k <- 0L

for (grp in colnames(pb)) {
  spl <- strsplit(grp, "\\|\\|\\|")[[1]]
  s <- spl[1]
  st <- paste(spl[-1], collapse = "|||")

  for (g in genes) {
    k <- k + 1L
    gene_rows[[k]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      HSC_state3 = st,
      gene = g,
      cp10k = pb[g, grp],
      log1p_cp10k = log1p(pb[g, grp]),
      stringsAsFactors = FALSE
    )
  }
}

gene_pb <- do.call(rbind, gene_rows)
write_tsv(gene_pb, file.path(OUTDIR, "06_gene_pseudobulk_by_HSC_state.tsv"))

# Panel score = mean log1p(CP10k) of available genes
panel_rows <- list()
k <- 0L

for (grp in colnames(pb)) {
  spl <- strsplit(grp, "\\|\\|\\|")[[1]]
  s <- spl[1]
  st <- paste(spl[-1], collapse = "|||")

  for (pn in names(PANELS)) {
    g <- intersect(PANELS[[pn]], rownames(pb))
    val <- mean(log1p(pb[g, grp]))

    k <- k + 1L
    panel_rows[[k]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      HSC_state3 = st,
      panel = pn,
      panel_mean_log1p_CP10k = val,
      stringsAsFactors = FALSE
    )
  }
}

panel_pb <- do.call(rbind, panel_rows)
write_tsv(panel_pb, file.path(OUTDIR, "07_panel_pseudobulk_by_HSC_state.tsv"))

# ------------------------------------------------------------------------------
# Group effects: panels
# ------------------------------------------------------------------------------

keys <- unique(panel_pb[, c("HSC_state3", "panel")])
panel_effect_rows <- list()
k <- 0L

for (i in seq_len(nrow(keys))) {
  st <- keys$HSC_state3[i]
  pn <- keys$panel[i]

  d <- panel_pb[
    panel_pb$HSC_state3 == st &
      panel_pb$panel == pn,
    ,
    drop = FALSE
  ]

  sham <- d$panel_mean_log1p_CP10k[d$group == "Sham"]
  tx <- d$panel_mean_log1p_CP10k[d$group == "Tx"]

  k <- k + 1L
  panel_effect_rows[[k]] <- data.frame(
    HSC_state3 = st,
    panel = pn,
    Sham_mean = safe_mean(sham),
    Tx_mean = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    percent_change_vs_abs_Sham = ifelse(
      safe_mean(sham) == 0,
      NA_real_,
      100 * (safe_mean(tx) - safe_mean(sham)) / abs(safe_mean(sham))
    ),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

panel_effects <- do.call(rbind, panel_effect_rows)
write_tsv(panel_effects, file.path(OUTDIR, "08_panel_group_effects.tsv"))

# ------------------------------------------------------------------------------
# Group effects: direct genes
# ------------------------------------------------------------------------------

keys <- unique(gene_pb[, c("HSC_state3", "gene")])
gene_effect_rows <- list()
k <- 0L

for (i in seq_len(nrow(keys))) {
  st <- keys$HSC_state3[i]
  g <- keys$gene[i]

  d <- gene_pb[
    gene_pb$HSC_state3 == st &
      gene_pb$gene == g,
    ,
    drop = FALSE
  ]

  sham <- d$cp10k[d$group == "Sham"]
  tx <- d$cp10k[d$group == "Tx"]

  sham_mean <- safe_mean(sham)
  tx_mean <- safe_mean(tx)

  k <- k + 1L
  gene_effect_rows[[k]] <- data.frame(
    HSC_state3 = st,
    gene = g,
    Sham_mean_CP10k = sham_mean,
    Tx_mean_CP10k = tx_mean,
    delta_Tx_minus_Sham = tx_mean - sham_mean,
    ratio_Tx_over_Sham = (tx_mean + 0.01) / (sham_mean + 0.01),
    log2_ratio_Tx_over_Sham = log2(
      (tx_mean + 0.01) / (sham_mean + 0.01)
    ),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

gene_effects <- do.call(rbind, gene_effect_rows)

write_tsv(
  gene_effects,
  file.path(OUTDIR, "09_gene_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Key mechanistic readout
# ------------------------------------------------------------------------------

KEY_GENES <- c(
  "Lrat", "Rbp1", "Cygb",
  "Col1a1", "Col1a2", "Col3a1", "Fn1", "Timp1", "Ctgf", "Serpine1",
  "Acta2", "Tagln",
  "Tgfbr1", "Tgfbr2",
  "Pdgfra", "Pdgfrb",
  "Cd44", "Itgav", "Itgb1", "Itgb3",
  "Itga5",
  "Egfr",
  "Il10ra", "Il10rb",
  "Met",
  "Bmpr1a", "Bmpr2",
  "Axl", "Mertk", "Tyro3",
  "Lrp1"
)

key <- gene_effects[
  gene_effects$gene %in% KEY_GENES,
  ,
  drop = FALSE
]

key <- key[
  order(
    key$HSC_state3,
    -abs(key$log2_ratio_Tx_over_Sham)
  ),
  ,
  drop = FALSE
]

write_tsv(
  key,
  file.path(OUTDIR, "10_KEY_HSC_response.tsv")
)

# ------------------------------------------------------------------------------
# Integrity
# ------------------------------------------------------------------------------

checks <- data.frame(
  item = c(
    "PrimaryClean_HSC_cells",
    "all_qc_primary_include_TRUE",
    "all_four_samples_present",
    "cell_level_inferential_pvalues_used",
    "biological_replicate",
    "interpretation"
  ),
  value = c(
    ncol(obj),
    as.character(all(obj$qc_primary_include %in% TRUE)),
    as.character(all(SAMPLES %in% unique(sample))),
    "FALSE",
    "sample",
    "Receiver-side descriptive validation of HSC state and pseudobulk expression"
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
message("  04_HSC_state_group_effects.tsv")
message("  08_panel_group_effects.tsv")
message("  10_KEY_HSC_response.tsv")
message("  11_analysis_integrity_summary.tsv")
message("============================================================")
