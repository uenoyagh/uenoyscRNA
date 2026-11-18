#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS <- file.path(
  ROOT, "Mouse_MASH_Mphi_RDS", "Mouse_Mphi_PrimaryClean_v5.0",
  "objects", "Mouse_MASH_Mphi_PrimaryClean_v5.0.rds"
)

HSC_RDS <- file.path(
  ROOT, "Mouse_MASH_RDS", "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects", "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

OUTDIR <- file.path(
  ROOT, "Mouse_MASH_Interaction",
  "Mouse_CrossLineage_Mphi_HSC_ECM_MechanisticClosure_v1.0"
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
  if (!inherits(mat, "Matrix")) mat <- Matrix(mat, sparse = TRUE)
  mat
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(Sham1 = "Sham", Sham20 = "Sham", Tx17 = "Tx", Tx5 = "Tx")

if (!file.exists(MPHI_RDS)) stop("Missing Mphi RDS: ", MPHI_RDS)
if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)

message("Loading Primary-Clean Mphi...")
mphi <- readRDS(MPHI_RDS)

message("Loading Primary-Clean HSC...")
hsc <- readRDS(HSC_RDS)

if (ncol(hsc) != 12941) stop("Unexpected Primary-Clean HSC cell count: ", ncol(hsc))

MPHI_PANELS <- list(
  Mphi_matrix_degradation = c("Mmp12","Mmp9","Mmp13","Mmp14","Mmp19"),
  Mphi_phagolysosomal_remodeling = c("Ctsb","Ctsd","Ctsk","Lgals3"),
  Mphi_direct_HSC_activating_ligands = c("Tgfb1","Pdgfa","Pdgfb","Spp1"),
  Mphi_matrix_adhesion_remodeling = c("Fn1","Lgals3","Gpnmb"),
  Mphi_resolution_repair = c("Il1rn","Areg","Gdf15","Apoe","Igf1","Hbegf")
)

HSC_PANELS <- list(
  HSC_quiescence = c("Lrat","Rbp1","Reln","Cygb","Colec11","Vipr1"),
  HSC_contractile = c("Acta2","Tagln","Myl9","Tpm2","Myh11"),
  HSC_collagen_synthesis = c("Col1a1","Col1a2","Col3a1","Col5a1"),
  HSC_collagen_maturation = c(
    "P4ha1","P4ha2","P4hb","Plod1","Plod2",
    "Serpinh1","Pcolce","Pcolce2","Lox","Loxl1","Loxl2"
  ),
  HSC_fibrogenic_stress = c("Ctgf","Serpine1","Tgfb1","Tgfbi"),
  HSC_matrix_degradation = c(
    "Mmp2","Mmp3","Mmp9","Mmp13","Mmp14","Mmp19",
    "Adamts2","Adamts4","Adamts5","Adamts12","Adamts14"
  )
)

analyze_lineage <- function(obj, sample_col, panels, lineage_label) {
  if (!(sample_col %in% colnames(obj@meta.data))) {
    stop(lineage_label, ": missing sample column ", sample_col)
  }

  smp <- as.character(obj@meta.data[[sample_col]])
  keep <- smp %in% SAMPLES
  obj <- subset(obj, cells = colnames(obj)[keep])
  smp <- as.character(obj@meta.data[[sample_col]])

  counts <- get_joined_counts(obj)
  counts <- counts[, colnames(obj), drop = FALSE]

  coverage <- do.call(rbind, lapply(names(panels), function(pn) {
    g <- panels[[pn]]
    present <- intersect(g, rownames(counts))
    data.frame(
      lineage = lineage_label,
      panel = pn,
      n_requested = length(g),
      n_present = length(present),
      coverage_fraction = length(present) / length(g),
      genes_present = paste(present, collapse = ";"),
      genes_missing = paste(setdiff(g, present), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  sample_rows <- list()
  k <- 0L

  for (s in SAMPLES) {
    ix <- which(smp == s)
    if (!length(ix)) next

    sums <- Matrix::rowSums(counts[, ix, drop = FALSE])
    lib <- sum(sums)
    cp <- as.numeric(sums) / lib * 1e4
    names(cp) <- rownames(counts)

    for (pn in names(panels)) {
      gg <- intersect(panels[[pn]], rownames(counts))
      vals <- cp[gg]
      k <- k + 1L
      sample_rows[[k]] <- data.frame(
        lineage = lineage_label,
        sample = s,
        group = GROUP[[s]],
        panel = pn,
        mean_log1p_CP10k = if (length(vals)) mean(log1p(vals)) else NA_real_,
        sum_CP10k = if (length(vals)) sum(vals) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  sample_df <- do.call(rbind, sample_rows)

  effect_rows <- list()
  k <- 0L

  for (pn in unique(sample_df$panel)) {
    d <- sample_df[sample_df$panel == pn, , drop = FALSE]
    sham <- d$mean_log1p_CP10k[d$group == "Sham"]
    tx <- d$mean_log1p_CP10k[d$group == "Tx"]
    sham_sum <- d$sum_CP10k[d$group == "Sham"]
    tx_sum <- d$sum_CP10k[d$group == "Tx"]

    k <- k + 1L
    effect_rows[[k]] <- data.frame(
      lineage = lineage_label,
      panel = pn,
      Sham_mean_log_score = safe_mean(sham),
      Tx_mean_log_score = safe_mean(tx),
      delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
      direction_concordance = direction_rule(sham, tx),
      Sham_mean_sum_CP10k = safe_mean(sham_sum),
      Tx_mean_sum_CP10k = safe_mean(tx_sum),
      sum_CP10k_ratio_Tx_over_Sham =
        (safe_mean(tx_sum) + 0.01) / (safe_mean(sham_sum) + 0.01),
      sum_CP10k_log2_ratio = log2(
        (safe_mean(tx_sum) + 0.01) / (safe_mean(sham_sum) + 0.01)
      ),
      sum_CP10k_direction = direction_rule(sham_sum, tx_sum),
      stringsAsFactors = FALSE
    )
  }

  list(
    coverage = coverage,
    sample = sample_df,
    effects = do.call(rbind, effect_rows),
    n_cells = ncol(obj),
    sample_values = unique(smp)
  )
}

mphi_res <- analyze_lineage(
  mphi, "sample_display_FIXED2", MPHI_PANELS, "Mphi"
)

hsc_res <- analyze_lineage(
  hsc, "sample_hsc_v610", HSC_PANELS, "HSC"
)

write_tsv(rbind(mphi_res$coverage, hsc_res$coverage),
          file.path(OUTDIR, "01_panel_coverage.tsv"))
write_tsv(mphi_res$sample, file.path(OUTDIR, "02_Mphi_sample_panel_scores.tsv"))
write_tsv(hsc_res$sample, file.path(OUTDIR, "03_HSC_sample_panel_scores.tsv"))
write_tsv(mphi_res$effects, file.path(OUTDIR, "04_Mphi_group_effects.tsv"))
write_tsv(hsc_res$effects, file.path(OUTDIR, "05_HSC_group_effects.tsv"))

hmd <- hsc@meta.data
hs <- as.character(hmd$sample_hsc_v610)
hst <- as.character(hmd$HSC_state3_v612)
keep <- hs %in% SAMPLES
hs <- hs[keep]
hst <- hst[keep]

state_tab <- as.data.frame(table(sample = hs, HSC_state3 = hst),
                           stringsAsFactors = FALSE)
state_tot <- aggregate(Freq ~ sample, data = state_tab, FUN = sum)
names(state_tot)[2] <- "sample_total"
state_tab <- merge(state_tab, state_tot, by = "sample", all.x = TRUE, sort = FALSE)
state_tab$fraction <- state_tab$Freq / state_tab$sample_total
state_tab$group <- GROUP[state_tab$sample]

write_tsv(state_tab, file.path(OUTDIR, "06_HSC_state_composition_by_sample.tsv"))

state_effect_rows <- list()
k <- 0L
for (st in unique(state_tab$HSC_state3)) {
  d <- state_tab[state_tab$HSC_state3 == st, , drop = FALSE]
  sham <- d$fraction[d$group == "Sham"]
  tx <- d$fraction[d$group == "Tx"]

  k <- k + 1L
  state_effect_rows[[k]] <- data.frame(
    component = paste0("HSC_fraction_", st),
    lineage = "HSC",
    Sham_mean = safe_mean(sham),
    Tx_mean = safe_mean(tx),
    delta_Tx_minus_Sham = safe_mean(tx) - safe_mean(sham),
    direction_concordance = direction_rule(sham, tx),
    stringsAsFactors = FALSE
  )
}

state_effects <- do.call(rbind, state_effect_rows)
write_tsv(state_effects,
          file.path(OUTDIR, "07_HSC_state_composition_group_effects.tsv"))

integrated <- rbind(
  data.frame(
    component = mphi_res$effects$panel,
    lineage = "Mphi",
    Sham_mean_log_score = mphi_res$effects$Sham_mean_log_score,
    Tx_mean_log_score = mphi_res$effects$Tx_mean_log_score,
    delta_Tx_minus_Sham = mphi_res$effects$delta_Tx_minus_Sham,
    direction_concordance = mphi_res$effects$direction_concordance,
    stringsAsFactors = FALSE
  ),
  data.frame(
    component = hsc_res$effects$panel,
    lineage = "HSC",
    Sham_mean_log_score = hsc_res$effects$Sham_mean_log_score,
    Tx_mean_log_score = hsc_res$effects$Tx_mean_log_score,
    delta_Tx_minus_Sham = hsc_res$effects$delta_Tx_minus_Sham,
    direction_concordance = hsc_res$effects$direction_concordance,
    stringsAsFactors = FALSE
  )
)

write_tsv(integrated,
          file.path(OUTDIR, "08_integrated_Mphi_HSC_mechanistic_closure.tsv"))

checks <- data.frame(
  item = c(
    "Mphi_cells_ShamsTx",
    "HSC_cells_ShamsTx",
    "all_four_samples_Mphi",
    "all_four_samples_HSC",
    "biological_replicate",
    "cell_level_inferential_pvalues_used",
    "cross_lineage_absolute_CP10k_compared",
    "interpretation"
  ),
  value = c(
    mphi_res$n_cells,
    hsc_res$n_cells,
    as.character(all(SAMPLES %in% mphi_res$sample_values)),
    as.character(all(SAMPLES %in% hsc_res$sample_values)),
    "sample",
    "FALSE",
    "FALSE",
    "Within-lineage Sham-vs-Tx effects integrated across Mphi and HSC; no synthetic causal score"
  ),
  stringsAsFactors = FALSE
)

write_tsv(checks, file.path(OUTDIR, "09_analysis_integrity_summary.tsv"))
writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "sessionInfo.txt"), useBytes = TRUE)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  04_Mphi_group_effects.tsv")
message("  05_HSC_group_effects.tsv")
message("  07_HSC_state_composition_group_effects.tsv")
message("  08_integrated_Mphi_HSC_mechanistic_closure.tsv")
message("  09_analysis_integrity_summary.tsv")
message("============================================================")
