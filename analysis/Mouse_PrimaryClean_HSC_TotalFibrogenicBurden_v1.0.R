#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_HSC_TotalFibrogenicBurden_v1.0.R
#
# PURPOSE
# -------
# Resolve the apparent contrast between:
#   1) favorable HSC composition shift in Tx
#      (qHSC up, ECM-activated HSC down)
#   2) increased within-state ECM/remodeling transcripts in some HSC states.
#
# The analysis quantifies:
#   - ALL-HSC sample pseudobulk transcript burden
#   - state contribution to the total HSC transcript burden
#   - PrimaryClean vs SensitivityClean robustness
#
# Biological replicate = SAMPLE (Sham1, Sham20, Tx17, Tx5).
# No cell-level inferential p-values.
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

PRIMARY_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

SENS_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_SensitivityClean_v7.0.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_HSC",
  "Mouse_PrimaryClean_HSC_TotalFibrogenicBurden_v1.0"
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

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(Sham1 = "Sham", Sham20 = "Sham", Tx17 = "Tx", Tx5 = "Tx")

PANELS <- list(
  Collagen_core = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1"
  ),

  ECM_structural = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1",
    "Lum", "Dcn", "Bgn"
  ),

  Matrix_remodeling = c(
    "Fn1", "Timp1", "Serpine1", "Ctgf"
  ),

  Fibrogenic_combined = c(
    "Col1a1", "Col1a2", "Col3a1", "Col5a1",
    "Fn1", "Timp1", "Serpine1", "Ctgf"
  ),

  Contractile = c(
    "Acta2", "Tagln", "Myl9", "Tpm2", "Myh11"
  ),

  Quiescence_identity = c(
    "Lrat", "Rbp1", "Reln", "Cygb", "Colec11", "Vipr1"
  )
)

KEY_GENES <- unique(unlist(PANELS))

analyze_object <- function(obj, dataset_label) {

  md <- obj@meta.data

  required <- c(
    "sample_hsc_v610",
    "HSC_state3_v612"
  )

  missing <- setdiff(required, colnames(md))
  if (length(missing)) {
    stop(dataset_label, " missing metadata: ", paste(missing, collapse = ", "))
  }

  sample <- as.character(md$sample_hsc_v610)
  state3 <- as.character(md$HSC_state3_v612)

  keep <- sample %in% SAMPLES

  obj <- subset(obj, cells = rownames(md)[keep])
  md <- obj@meta.data

  sample <- as.character(md$sample_hsc_v610)
  state3 <- as.character(md$HSC_state3_v612)

  counts <- get_joined_counts(obj)
  counts <- counts[, colnames(obj), drop = FALSE]

  genes <- intersect(KEY_GENES, rownames(counts))

  # --------------------------------------------------------------------------
  # Cell counts / state composition
  # --------------------------------------------------------------------------

  comp <- as.data.frame(
    table(sample = sample, HSC_state3 = state3),
    stringsAsFactors = FALSE
  )

  tot <- aggregate(Freq ~ sample, data = comp, FUN = sum)
  names(tot)[2] <- "sample_total"

  comp <- merge(comp, tot, by = "sample", all.x = TRUE, sort = FALSE)
  comp$fraction <- comp$Freq / comp$sample_total
  comp$group <- GROUP[comp$sample]
  comp$dataset <- dataset_label

  # --------------------------------------------------------------------------
  # ALL-HSC sample pseudobulk
  # --------------------------------------------------------------------------

  all_gene_rows <- list()
  all_panel_rows <- list()
  state_gene_rows <- list()
  state_panel_rows <- list()

  kg <- 0L
  kp <- 0L
  ksg <- 0L
  ksp <- 0L

  for (s in SAMPLES) {

    ix_all <- which(sample == s)
    if (!length(ix_all)) next

    total_counts <- Matrix::rowSums(counts[, ix_all, drop = FALSE])
    total_lib <- sum(total_counts)

    cp10k_all <- as.numeric(total_counts) / total_lib * 1e4
    names(cp10k_all) <- rownames(counts)

    for (g in genes) {
      kg <- kg + 1L
      all_gene_rows[[kg]] <- data.frame(
        dataset = dataset_label,
        sample = s,
        group = GROUP[[s]],
        gene = g,
        cp10k = cp10k_all[g],
        log1p_cp10k = log1p(cp10k_all[g]),
        stringsAsFactors = FALSE
      )
    }

    for (pn in names(PANELS)) {
      gg <- intersect(PANELS[[pn]], rownames(counts))
      vals <- cp10k_all[gg]

      kp <- kp + 1L
      all_panel_rows[[kp]] <- data.frame(
        dataset = dataset_label,
        sample = s,
        group = GROUP[[s]],
        panel = pn,
        panel_sum_CP10k = sum(vals),
        panel_mean_log1p_CP10k = mean(log1p(vals)),
        stringsAsFactors = FALSE
      )
    }

    # ------------------------------------------------------------------------
    # State contribution to TOTAL HSC transcript pool
    #
    # State contribution CP10k uses the TOTAL HSC library denominator.
    # Therefore contributions across HSC states add to ALL-HSC CP10k.
    # ------------------------------------------------------------------------

    for (st in sort(unique(state3[ix_all]))) {

      ix_st <- which(sample == s & state3 == st)
      if (!length(ix_st)) next

      st_counts <- Matrix::rowSums(counts[, ix_st, drop = FALSE])
      st_contrib_cp10k <- as.numeric(st_counts) / total_lib * 1e4
      names(st_contrib_cp10k) <- rownames(counts)

      state_frac <- length(ix_st) / length(ix_all)

      for (g in genes) {
        ksg <- ksg + 1L
        state_gene_rows[[ksg]] <- data.frame(
          dataset = dataset_label,
          sample = s,
          group = GROUP[[s]],
          HSC_state3 = st,
          state_cell_fraction = state_frac,
          gene = g,
          contribution_CP10k = st_contrib_cp10k[g],
          stringsAsFactors = FALSE
        )
      }

      for (pn in names(PANELS)) {
        gg <- intersect(PANELS[[pn]], rownames(counts))
        vals <- st_contrib_cp10k[gg]

        ksp <- ksp + 1L
        state_panel_rows[[ksp]] <- data.frame(
          dataset = dataset_label,
          sample = s,
          group = GROUP[[s]],
          HSC_state3 = st,
          state_cell_fraction = state_frac,
          panel = pn,
          contribution_sum_CP10k = sum(vals),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    composition = comp,
    all_gene = do.call(rbind, all_gene_rows),
    all_panel = do.call(rbind, all_panel_rows),
    state_gene = do.call(rbind, state_gene_rows),
    state_panel = do.call(rbind, state_panel_rows)
  )
}

if (!file.exists(PRIMARY_RDS)) stop("Missing Primary RDS: ", PRIMARY_RDS)
if (!file.exists(SENS_RDS)) stop("Missing Sensitivity RDS: ", SENS_RDS)

message("Loading PrimaryClean HSC...")
primary <- readRDS(PRIMARY_RDS)

message("Loading SensitivityClean HSC...")
sens <- readRDS(SENS_RDS)

message("Analyzing PrimaryClean...")
p <- analyze_object(primary, "PrimaryClean")

message("Analyzing SensitivityClean...")
s <- analyze_object(sens, "SensitivityClean")

composition <- rbind(p$composition, s$composition)
all_gene <- rbind(p$all_gene, s$all_gene)
all_panel <- rbind(p$all_panel, s$all_panel)
state_gene <- rbind(p$state_gene, s$state_gene)
state_panel <- rbind(p$state_panel, s$state_panel)

write_tsv(
  composition,
  file.path(OUTDIR, "01_HSC_state_composition_primary_sensitivity.tsv")
)

write_tsv(
  all_gene,
  file.path(OUTDIR, "02_ALLHSC_gene_pseudobulk.tsv")
)

write_tsv(
  all_panel,
  file.path(OUTDIR, "03_ALLHSC_panel_pseudobulk.tsv")
)

write_tsv(
  state_gene,
  file.path(OUTDIR, "04_state_gene_contribution_to_totalHSC.tsv")
)

write_tsv(
  state_panel,
  file.path(OUTDIR, "05_state_panel_contribution_to_totalHSC.tsv")
)

# ------------------------------------------------------------------------------
# Group effects: ALL-HSC panels
# ------------------------------------------------------------------------------

panel_effect_rows <- list()
k <- 0L

for (ds in unique(all_panel$dataset)) {
  for (pn in unique(all_panel$panel)) {

    d <- all_panel[
      all_panel$dataset == ds &
        all_panel$panel == pn,
      ,
      drop = FALSE
    ]

    sham_sum <- d$panel_sum_CP10k[d$group == "Sham"]
    tx_sum <- d$panel_sum_CP10k[d$group == "Tx"]

    sham_log <- d$panel_mean_log1p_CP10k[d$group == "Sham"]
    tx_log <- d$panel_mean_log1p_CP10k[d$group == "Tx"]

    k <- k + 1L
    panel_effect_rows[[k]] <- data.frame(
      dataset = ds,
      panel = pn,
      Sham_mean_sum_CP10k = safe_mean(sham_sum),
      Tx_mean_sum_CP10k = safe_mean(tx_sum),
      sum_CP10k_ratio_Tx_over_Sham =
        (safe_mean(tx_sum) + 0.01) / (safe_mean(sham_sum) + 0.01),
      sum_CP10k_log2_ratio = log2(
        (safe_mean(tx_sum) + 0.01) / (safe_mean(sham_sum) + 0.01)
      ),
      sum_CP10k_direction = direction_rule(sham_sum, tx_sum),
      Sham_mean_panel_log = safe_mean(sham_log),
      Tx_mean_panel_log = safe_mean(tx_log),
      panel_log_delta = safe_mean(tx_log) - safe_mean(sham_log),
      panel_log_direction = direction_rule(sham_log, tx_log),
      stringsAsFactors = FALSE
    )
  }
}

panel_effects <- do.call(rbind, panel_effect_rows)

write_tsv(
  panel_effects,
  file.path(OUTDIR, "06_ALLHSC_panel_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Group effects: ALL-HSC key genes
# ------------------------------------------------------------------------------

gene_effect_rows <- list()
k <- 0L

for (ds in unique(all_gene$dataset)) {
  for (g in unique(all_gene$gene)) {

    d <- all_gene[
      all_gene$dataset == ds &
        all_gene$gene == g,
      ,
      drop = FALSE
    ]

    sham <- d$cp10k[d$group == "Sham"]
    tx <- d$cp10k[d$group == "Tx"]

    sm <- safe_mean(sham)
    tm <- safe_mean(tx)

    k <- k + 1L
    gene_effect_rows[[k]] <- data.frame(
      dataset = ds,
      gene = g,
      Sham_mean_CP10k = sm,
      Tx_mean_CP10k = tm,
      ratio_Tx_over_Sham = (tm + 0.01) / (sm + 0.01),
      log2_ratio_Tx_over_Sham = log2((tm + 0.01) / (sm + 0.01)),
      direction_concordance = direction_rule(sham, tx),
      stringsAsFactors = FALSE
    )
  }
}

gene_effects <- do.call(rbind, gene_effect_rows)

write_tsv(
  gene_effects,
  file.path(OUTDIR, "07_ALLHSC_gene_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Group effects: state contribution to total HSC panel burden
# ------------------------------------------------------------------------------

state_effect_rows <- list()
k <- 0L

keys <- unique(
  state_panel[, c("dataset", "HSC_state3", "panel")]
)

for (i in seq_len(nrow(keys))) {

  ds <- keys$dataset[i]
  st <- keys$HSC_state3[i]
  pn <- keys$panel[i]

  d <- state_panel[
    state_panel$dataset == ds &
      state_panel$HSC_state3 == st &
      state_panel$panel == pn,
    ,
    drop = FALSE
  ]

  sham <- d$contribution_sum_CP10k[d$group == "Sham"]
  tx <- d$contribution_sum_CP10k[d$group == "Tx"]

  sm <- safe_mean(sham)
  tm <- safe_mean(tx)

  k <- k + 1L
  state_effect_rows[[k]] <- data.frame(
    dataset = ds,
    HSC_state3 = st,
    panel = pn,
    Sham_mean_contribution_sum_CP10k = sm,
    Tx_mean_contribution_sum_CP10k = tm,
    ratio_Tx_over_Sham = (tm + 0.01) / (sm + 0.01),
    log2_ratio_Tx_over_Sham = log2((tm + 0.01) / (sm + 0.01)),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

state_effects <- do.call(rbind, state_effect_rows)

write_tsv(
  state_effects,
  file.path(OUTDIR, "08_state_contribution_group_effects.tsv")
)

# ------------------------------------------------------------------------------
# Primary vs sensitivity directional robustness
# ------------------------------------------------------------------------------

pfx <- panel_effects[panel_effects$dataset == "PrimaryClean", , drop = FALSE]
sfx <- panel_effects[panel_effects$dataset == "SensitivityClean", , drop = FALSE]

robust <- merge(
  pfx,
  sfx,
  by = "panel",
  suffixes = c("_Primary", "_Sensitivity"),
  all = TRUE
)

robust$same_log2_direction <- sign(robust$sum_CP10k_log2_ratio_Primary) ==
  sign(robust$sum_CP10k_log2_ratio_Sensitivity)

robust$same_concordance_class <- robust$sum_CP10k_direction_Primary ==
  robust$sum_CP10k_direction_Sensitivity

write_tsv(
  robust,
  file.path(OUTDIR, "09_Primary_vs_Sensitivity_panel_robustness.tsv")
)

# ------------------------------------------------------------------------------
# Compact key output
# ------------------------------------------------------------------------------

KEY_PANEL_ORDER <- c(
  "Collagen_core",
  "ECM_structural",
  "Matrix_remodeling",
  "Fibrogenic_combined",
  "Contractile",
  "Quiescence_identity"
)

compact <- panel_effects[
  panel_effects$panel %in% KEY_PANEL_ORDER,
  ,
  drop = FALSE
]

compact$panel <- factor(compact$panel, levels = KEY_PANEL_ORDER)

compact <- compact[
  order(compact$dataset, compact$panel),
  ,
  drop = FALSE
]

compact$panel <- as.character(compact$panel)

write_tsv(
  compact,
  file.path(OUTDIR, "10_KEY_total_HSC_burden_summary.tsv")
)

# ------------------------------------------------------------------------------
# Integrity
# ------------------------------------------------------------------------------

checks <- data.frame(
  item = c(
    "PrimaryClean_cells",
    "SensitivityClean_cells",
    "Primary_expected_cells",
    "Sensitivity_expected_cells",
    "all_four_samples_primary",
    "all_four_samples_sensitivity",
    "biological_replicate",
    "cell_level_inferential_pvalues_used",
    "burden_definition"
  ),
  value = c(
    ncol(primary),
    ncol(sens),
    12941,
    13026,
    as.character(all(SAMPLES %in% unique(primary$sample_hsc_v610))),
    as.character(all(SAMPLES %in% unique(sens$sample_hsc_v610))),
    "sample",
    "FALSE",
    "ALL-HSC sample pseudobulk and state contribution using total-HSC library denominator"
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
message("  06_ALLHSC_panel_group_effects.tsv")
message("  07_ALLHSC_gene_group_effects.tsv")
message("  08_state_contribution_group_effects.tsv")
message("  09_Primary_vs_Sensitivity_panel_robustness.tsv")
message("  10_KEY_total_HSC_burden_summary.tsv")
message("  11_analysis_integrity_summary.tsv")
message("============================================================")
