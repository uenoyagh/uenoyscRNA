#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_Mphi_to_HSC_LR_Compatibility_v1.0.R
#
# PURPOSE
# -------
# Sample-level, descriptive ligand-receptor compatibility analysis:
#
#   Sender   = Mouse Primary-Clean macrophages
#   Receiver = Mouse Primary-Clean HSC
#   Contrast = Sham1/Sham20 vs Tx17/Tx5
#
# Statistical policy:
#   - biological replicate = SAMPLE, not cell
#   - n = 2 samples/group
#   - no cell-level inferential p-values
#   - report sample-level effect sizes and direction concordance
#
# IMPORTANT
# ---------
# This is NOT a claim of physical ligand-receptor binding or signaling flux.
# The "LR compatibility score" is a descriptive co-expression proxy:
#
#   sqrt(sender ligand CP10k * receiver receptor CP10k)
#
# calculated from sample-level pseudobulk counts.
#
# Canonical inputs:
#   MΦ  : Mouse_Mphi_PrimaryClean_v5.0.rds
#   HSC : Mouse_MASH_HSC_PrimaryClean_v7.0.rds
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

MPHI_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mouse_Mphi_PrimaryClean_v5.0",
  "objects",
  "Mouse_MASH_Mphi_PrimaryClean_v5.0.rds"
)

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0",
  "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

OUTDIR <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mouse_PrimaryClean_Mphi_to_HSC_LR_Compatibility_v1.0"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

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

direction_rule <- function(sham, tx) {
  sham <- sham[is.finite(sham)]
  tx <- tx[is.finite(tx)]

  if (length(sham) < 2 || length(tx) < 2) return("insufficient")

  if (min(tx) > max(sham)) return("Tx_above_all_Sham")
  if (max(tx) < min(sham)) return("Tx_below_all_Sham")
  "overlap"
}

# ------------------------------------------------------------------------------
# Load canonical Primary-Clean objects
# ------------------------------------------------------------------------------

if (!file.exists(MPHI_RDS)) stop("Missing MΦ RDS: ", MPHI_RDS)
if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)

message("Loading Primary-Clean MΦ...")
mphi <- readRDS(MPHI_RDS)

message("Loading Primary-Clean HSC...")
hsc <- readRDS(HSC_RDS)

mmd <- mphi@meta.data
hmd <- hsc@meta.data

# ------------------------------------------------------------------------------
# Detect metadata columns
# ------------------------------------------------------------------------------

detect_sample_col <- function(md, preferred) {
  for (cn in preferred) {
    if (cn %in% colnames(md)) {
      vals <- unique(as.character(md[[cn]]))
      if (all(c("Sham1", "Sham20", "Tx17", "Tx5") %in% vals)) return(cn)
    }
  }

  for (cn in colnames(md)) {
    vals <- unique(as.character(md[[cn]]))
    if (all(c("Sham1", "Sham20", "Tx17", "Tx5") %in% vals)) return(cn)
  }

  NA_character_
}

detect_mphi_class_col <- function(md) {
  patterns <- c(
    "Anti-inflammatory",
    "Repair/Resolution",
    "Inflammatory",
    "ECM-associated",
    "LAM/TREM2"
  )

  score_col <- function(cn) {
    vals <- unique(as.character(md[[cn]]))
    sum(vapply(
      patterns,
      function(p) any(grepl(p, vals, fixed = TRUE)),
      logical(1)
    ))
  }

  scores <- vapply(colnames(md), score_col, integer(1))
  if (!length(scores) || max(scores) < 3) return(NA_character_)

  names(which.max(scores))
}

M_SAMPLE_COL <- detect_sample_col(
  mmd,
  c("sample_display_FIXED2", "sample", "orig.ident", "sample_display")
)

H_SAMPLE_COL <- detect_sample_col(
  hmd,
  c("sample_hsc_v610", "sample", "sample_display_FIXED2", "orig.ident")
)

M_CLASS_COL <- detect_mphi_class_col(mmd)

H_STATE_COL <- if ("HSC_state3_v612" %in% colnames(hmd)) {
  "HSC_state3_v612"
} else {
  stop("HSC_state3_v612 not found.")
}

if (is.na(M_SAMPLE_COL)) stop("Could not detect MΦ sample column.")
if (is.na(H_SAMPLE_COL)) stop("Could not detect HSC sample column.")
if (is.na(M_CLASS_COL)) {
  stop(
    "Could not detect final MΦ class column. ",
    "Inspect metadata before LR analysis."
  )
}

meta_detect <- data.frame(
  item = c(
    "Mphi_RDS",
    "HSC_RDS",
    "Mphi_cells_total",
    "HSC_cells_total",
    "Mphi_sample_col",
    "Mphi_final_class_col",
    "HSC_sample_col",
    "HSC_state_col"
  ),
  value = c(
    MPHI_RDS,
    HSC_RDS,
    ncol(mphi),
    ncol(hsc),
    M_SAMPLE_COL,
    M_CLASS_COL,
    H_SAMPLE_COL,
    H_STATE_COL
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  meta_detect,
  file.path(OUTDIR, "01_metadata_detection.tsv")
)

# ------------------------------------------------------------------------------
# Restrict to Sham / Tx
# ------------------------------------------------------------------------------

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
GROUP <- c(
  Sham1 = "Sham",
  Sham20 = "Sham",
  Tx17 = "Tx",
  Tx5 = "Tx"
)

m_sample <- as.character(mmd[[M_SAMPLE_COL]])
h_sample <- as.character(hmd[[H_SAMPLE_COL]])

m_cells <- rownames(mmd)[m_sample %in% SAMPLES]
h_cells <- rownames(hmd)[h_sample %in% SAMPLES]

mphi <- subset(mphi, cells = m_cells)
hsc <- subset(hsc, cells = h_cells)

mmd <- mphi@meta.data
hmd <- hsc@meta.data

m_sample <- as.character(mmd[[M_SAMPLE_COL]])
h_sample <- as.character(hmd[[H_SAMPLE_COL]])
m_class <- as.character(mmd[[M_CLASS_COL]])
h_state <- as.character(hmd[[H_STATE_COL]])

# ------------------------------------------------------------------------------
# Integrity checks
# ------------------------------------------------------------------------------

if ("qc_primary_include" %in% colnames(mmd)) {
  if (!all(mmd$qc_primary_include %in% TRUE)) {
    stop("MΦ Primary-Clean object contains qc_primary_include != TRUE.")
  }
}

if ("qc_primary_include" %in% colnames(hmd)) {
  if (!all(hmd$qc_primary_include %in% TRUE)) {
    stop("HSC Primary-Clean object contains qc_primary_include != TRUE.")
  }
}

m_counts_tab <- as.data.frame(
  table(sample = m_sample, sender_class = m_class),
  stringsAsFactors = FALSE
)
m_counts_tab <- m_counts_tab[m_counts_tab$Freq > 0, , drop = FALSE]

h_counts_tab <- as.data.frame(
  table(sample = h_sample, HSC_state3 = h_state),
  stringsAsFactors = FALSE
)
h_counts_tab <- h_counts_tab[h_counts_tab$Freq > 0, , drop = FALSE]

write_tsv(
  m_counts_tab,
  file.path(OUTDIR, "02_sender_Mphi_cell_counts.tsv")
)

write_tsv(
  h_counts_tab,
  file.path(OUTDIR, "03_receiver_HSC_cell_counts.tsv")
)

# ------------------------------------------------------------------------------
# Seurat v5 counts helper
# ------------------------------------------------------------------------------

get_joined_counts <- function(obj) {
  if (!("RNA" %in% Assays(obj))) stop("RNA assay missing.")

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

message("Joining MΦ RNA counts...")
m_counts <- get_joined_counts(mphi)

message("Joining HSC RNA counts...")
h_counts <- get_joined_counts(hsc)

m_counts <- m_counts[, colnames(mphi), drop = FALSE]
h_counts <- h_counts[, colnames(hsc), drop = FALSE]

# ------------------------------------------------------------------------------
# Candidate ligand-receptor pairs
# ------------------------------------------------------------------------------

LR <- data.frame(
  axis = c(
    "PROFIBROTIC", "PROFIBROTIC",
    "PROFIBROTIC", "PROFIBROTIC",
    "PROFIBROTIC", "PROFIBROTIC",
    "SPP1_ADHESION", "SPP1_ADHESION", "SPP1_ADHESION", "SPP1_ADHESION",
    "MATRIX_INTEGRIN", "MATRIX_INTEGRIN",
    "REPAIR_RESOLUTION", "REPAIR_RESOLUTION",
    "REPAIR_RESOLUTION",
    "REPAIR_RESOLUTION", "REPAIR_RESOLUTION", "REPAIR_RESOLUTION",
    "REPAIR_RESOLUTION", "REPAIR_RESOLUTION", "REPAIR_RESOLUTION",
    "REPAIR_RESOLUTION",
    "ANTI_INFLAMMATORY", "ANTI_INFLAMMATORY",
    "HB_EGF",
    "BMP7", "BMP7", "BMP7",
    "CHEMOKINE"
  ),
  ligand = c(
    "Tgfb1", "Tgfb1",
    "Pdgfa", "Pdgfa",
    "Pdgfb", "Pdgfb",
    "Spp1", "Spp1", "Spp1", "Spp1",
    "Fn1", "Fn1",
    "Areg", "Hgf",
    "Igf1",
    "Gas6", "Gas6", "Gas6",
    "Mfge8", "Mfge8", "Mfge8",
    "Apoe",
    "Il10", "Il10",
    "Hbegf",
    "Bmp7", "Bmp7", "Bmp7",
    "Cxcl12"
  ),
  receptor = c(
    "Tgfbr1", "Tgfbr2",
    "Pdgfra", "Pdgfrb",
    "Pdgfra", "Pdgfrb",
    "Cd44", "Itgav", "Itgb1", "Itgb3",
    "Itga5", "Itgb1",
    "Egfr", "Met",
    "Igf1r",
    "Axl", "Mertk", "Tyro3",
    "Itgav", "Itgb3", "Itgb5",
    "Lrp1",
    "Il10ra", "Il10rb",
    "Egfr",
    "Bmpr1a", "Bmpr1b", "Bmpr2",
    "Cxcr4"
  ),
  stringsAsFactors = FALSE
)

LR$ligand_in_Mphi <- LR$ligand %in% rownames(m_counts)
LR$receptor_in_HSC <- LR$receptor %in% rownames(h_counts)
LR$pair_measurable <- LR$ligand_in_Mphi & LR$receptor_in_HSC

write_tsv(
  LR,
  file.path(OUTDIR, "04_LR_pair_coverage.tsv")
)

LRm <- LR[LR$pair_measurable, , drop = FALSE]

if (!nrow(LRm)) stop("No measurable ligand-receptor pairs.")

# ------------------------------------------------------------------------------
# Pseudobulk helper
# ------------------------------------------------------------------------------

pb_cp10k <- function(counts_mat, groups) {
  groups <- as.character(groups)
  lev <- unique(groups)

  res <- lapply(lev, function(g) {
    ix <- which(groups == g)
    s <- Matrix::rowSums(counts_mat[, ix, drop = FALSE])
    lib <- sum(s)

    if (lib <= 0) {
      cp10k <- rep(0, length(s))
    } else {
      cp10k <- as.numeric(s) / lib * 1e4
    }

    names(cp10k) <- rownames(counts_mat)
    cp10k
  })

  mat <- do.call(cbind, res)
  colnames(mat) <- lev
  rownames(mat) <- rownames(counts_mat)
  mat
}

# ------------------------------------------------------------------------------
# Sender pseudobulk: sample × final class + sample × ALL_MPHI
# ------------------------------------------------------------------------------

sender_group <- paste(m_sample, m_class, sep = "|||")
m_pb_class <- pb_cp10k(m_counts, sender_group)

sender_rows <- list()
k <- 0L

ligands <- unique(LRm$ligand)

for (grp in colnames(m_pb_class)) {
  spl <- strsplit(grp, "\\|\\|\\|", fixed = FALSE)[[1]]
  s <- spl[1]
  cl <- paste(spl[-1], collapse = "|||")

  for (g in ligands) {
    k <- k + 1L
    sender_rows[[k]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      sender_class = cl,
      ligand = g,
      ligand_cp10k = m_pb_class[g, grp],
      ligand_log1p_cp10k = log1p(m_pb_class[g, grp]),
      stringsAsFactors = FALSE
    )
  }
}

# Add pooled ALL_MPHI by sample
for (s in SAMPLES) {
  ix <- which(m_sample == s)
  if (!length(ix)) next

  sums <- Matrix::rowSums(m_counts[, ix, drop = FALSE])
  lib <- sum(sums)
  cp <- as.numeric(sums) / lib * 1e4
  names(cp) <- rownames(m_counts)

  for (g in ligands) {
    k <- k + 1L
    sender_rows[[k]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      sender_class = "ALL_MPHI",
      ligand = g,
      ligand_cp10k = cp[g],
      ligand_log1p_cp10k = log1p(cp[g]),
      stringsAsFactors = FALSE
    )
  }
}

sender_pb <- do.call(rbind, sender_rows)

write_tsv(
  sender_pb,
  file.path(OUTDIR, "05_sender_ligand_pseudobulk.tsv")
)

# ------------------------------------------------------------------------------
# Receiver pseudobulk: sample × HSC state
# ------------------------------------------------------------------------------

receiver_group <- paste(h_sample, h_state, sep = "|||")
h_pb <- pb_cp10k(h_counts, receiver_group)

receiver_rows <- list()
k <- 0L

receptors <- unique(LRm$receptor)

for (grp in colnames(h_pb)) {
  spl <- strsplit(grp, "\\|\\|\\|", fixed = FALSE)[[1]]
  s <- spl[1]
  st <- paste(spl[-1], collapse = "|||")

  for (g in receptors) {
    k <- k + 1L
    receiver_rows[[k]] <- data.frame(
      sample = s,
      group = GROUP[[s]],
      HSC_state3 = st,
      receptor = g,
      receptor_cp10k = h_pb[g, grp],
      receptor_log1p_cp10k = log1p(h_pb[g, grp]),
      stringsAsFactors = FALSE
    )
  }
}

receiver_pb <- do.call(rbind, receiver_rows)

write_tsv(
  receiver_pb,
  file.path(OUTDIR, "06_receiver_receptor_pseudobulk.tsv")
)

# ------------------------------------------------------------------------------
# Sample-level LR compatibility
# ------------------------------------------------------------------------------

sender_classes <- unique(sender_pb$sender_class)
receiver_states <- unique(receiver_pb$HSC_state3)

lr_rows <- list()
k <- 0L

for (s in SAMPLES) {
  for (sc in sender_classes) {
    sd <- sender_pb[
      sender_pb$sample == s &
        sender_pb$sender_class == sc,
      ,
      drop = FALSE
    ]

    if (!nrow(sd)) next

    for (st in receiver_states) {
      rd <- receiver_pb[
        receiver_pb$sample == s &
          receiver_pb$HSC_state3 == st,
        ,
        drop = FALSE
      ]

      if (!nrow(rd)) next

      for (i in seq_len(nrow(LRm))) {
        lig <- LRm$ligand[i]
        rec <- LRm$receptor[i]

        lv <- sd$ligand_cp10k[sd$ligand == lig]
        rv <- rd$receptor_cp10k[rd$receptor == rec]

        if (!length(lv) || !length(rv)) next

        score <- sqrt(lv * rv)

        k <- k + 1L
        lr_rows[[k]] <- data.frame(
          sample = s,
          group = GROUP[[s]],
          sender_class = sc,
          HSC_state3 = st,
          axis = LRm$axis[i],
          ligand = lig,
          receptor = rec,
          ligand_cp10k = lv,
          receptor_cp10k = rv,
          LR_compatibility_geomean_cp10k = score,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

lr_sample <- do.call(rbind, lr_rows)

write_tsv(
  lr_sample,
  file.path(OUTDIR, "07_sample_LR_compatibility.tsv")
)

# ------------------------------------------------------------------------------
# Group effects: sample-balanced descriptive contrasts
# ------------------------------------------------------------------------------

keys <- unique(
  lr_sample[
    ,
    c("sender_class", "HSC_state3", "axis", "ligand", "receptor")
  ]
)

group_rows <- list()
k <- 0L

for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]

  d <- lr_sample[
    lr_sample$sender_class == key$sender_class &
      lr_sample$HSC_state3 == key$HSC_state3 &
      lr_sample$axis == key$axis &
      lr_sample$ligand == key$ligand &
      lr_sample$receptor == key$receptor,
    ,
    drop = FALSE
  ]

  sham <- d$LR_compatibility_geomean_cp10k[d$group == "Sham"]
  tx <- d$LR_compatibility_geomean_cp10k[d$group == "Tx"]

  sham_lig <- d$ligand_cp10k[d$group == "Sham"]
  tx_lig <- d$ligand_cp10k[d$group == "Tx"]

  sham_rec <- d$receptor_cp10k[d$group == "Sham"]
  tx_rec <- d$receptor_cp10k[d$group == "Tx"]

  sham_mean <- safe_mean(sham)
  tx_mean <- safe_mean(tx)

  k <- k + 1L
  group_rows[[k]] <- data.frame(
    sender_class = key$sender_class,
    HSC_state3 = key$HSC_state3,
    axis = key$axis,
    ligand = key$ligand,
    receptor = key$receptor,
    Sham_mean_LR = sham_mean,
    Tx_mean_LR = tx_mean,
    delta_Tx_minus_Sham = tx_mean - sham_mean,
    ratio_Tx_over_Sham = (tx_mean + 0.01) / (sham_mean + 0.01),
    log2_ratio_Tx_over_Sham = log2(
      (tx_mean + 0.01) / (sham_mean + 0.01)
    ),
    direction_concordance = direction_rule(sham, tx),
    Sham_mean_ligand_cp10k = safe_mean(sham_lig),
    Tx_mean_ligand_cp10k = safe_mean(tx_lig),
    ligand_log2_ratio = log2(
      (safe_mean(tx_lig) + 0.01) /
        (safe_mean(sham_lig) + 0.01)
    ),
    Sham_mean_receptor_cp10k = safe_mean(sham_rec),
    Tx_mean_receptor_cp10k = safe_mean(tx_rec),
    receptor_log2_ratio = log2(
      (safe_mean(tx_rec) + 0.01) /
        (safe_mean(sham_rec) + 0.01)
    ),
    stringsAsFactors = FALSE
  )
}

group_effects <- do.call(rbind, group_rows)

write_tsv(
  group_effects,
  file.path(OUTDIR, "08_group_LR_effects.tsv")
)

# ------------------------------------------------------------------------------
# Compact mechanistic summaries
# ------------------------------------------------------------------------------

# Priority = pooled MΦ plus anti-inflammatory / repair / inflammatory classes
priority_sender <- c(
  "ALL_MPHI",
  unique(m_class[
    grepl(
      "Anti-inflammatory|Repair/Resolution|Inflammatory|ECM-associated|LAM/TREM2",
      m_class
    )
  ])
)

priority_sender <- unique(priority_sender)

key_effects <- group_effects[
  group_effects$sender_class %in% priority_sender &
    group_effects$axis %in% c(
      "PROFIBROTIC",
      "SPP1_ADHESION",
      "MATRIX_INTEGRIN",
      "REPAIR_RESOLUTION",
      "ANTI_INFLAMMATORY",
      "HB_EGF",
      "BMP7"
    ),
  ,
  drop = FALSE
]

key_effects <- key_effects[
  order(
    key_effects$HSC_state3,
    key_effects$axis,
    -abs(key_effects$log2_ratio_Tx_over_Sham)
  ),
  ,
  drop = FALSE
]

write_tsv(
  key_effects,
  file.path(OUTDIR, "09_KEY_Mphi_HSC_LR_effects.tsv")
)

# Axis-level summary using mean of pair-level LR compatibility changes
axis_keys <- unique(
  group_effects[, c("sender_class", "HSC_state3", "axis")]
)

axis_rows <- list()
k <- 0L

for (i in seq_len(nrow(axis_keys))) {
  key <- axis_keys[i, ]

  d <- group_effects[
    group_effects$sender_class == key$sender_class &
      group_effects$HSC_state3 == key$HSC_state3 &
      group_effects$axis == key$axis,
    ,
    drop = FALSE
  ]

  k <- k + 1L
  axis_rows[[k]] <- data.frame(
    sender_class = key$sender_class,
    HSC_state3 = key$HSC_state3,
    axis = key$axis,
    n_pairs = nrow(d),
    mean_pair_log2_ratio = safe_mean(d$log2_ratio_Tx_over_Sham),
    median_pair_log2_ratio = median(
      d$log2_ratio_Tx_over_Sham,
      na.rm = TRUE
    ),
    n_Tx_above_all_Sham = sum(
      d$direction_concordance == "Tx_above_all_Sham"
    ),
    n_Tx_below_all_Sham = sum(
      d$direction_concordance == "Tx_below_all_Sham"
    ),
    stringsAsFactors = FALSE
  )
}

axis_summary <- do.call(rbind, axis_rows)

write_tsv(
  axis_summary,
  file.path(OUTDIR, "10_axis_level_summary.tsv")
)

# Top changes
top_changes <- group_effects[
  order(-abs(group_effects$log2_ratio_Tx_over_Sham)),
  ,
  drop = FALSE
]

write_tsv(
  head(top_changes, 100),
  file.path(OUTDIR, "11_TOP100_LR_changes.tsv")
)

# ------------------------------------------------------------------------------
# Final readiness / provenance
# ------------------------------------------------------------------------------

checks <- data.frame(
  item = c(
    "Mphi_primary_clean_cells_in_ShamsTx",
    "HSC_primary_clean_cells",
    "Mphi_sample_col",
    "Mphi_final_class_col",
    "HSC_sample_col",
    "HSC_state_col",
    "all_four_samples_present_Mphi",
    "all_four_samples_present_HSC",
    "measurable_LR_pairs",
    "cell_level_inferential_pvalues_used",
    "LR_score_interpretation"
  ),
  value = c(
    ncol(mphi),
    ncol(hsc),
    M_SAMPLE_COL,
    M_CLASS_COL,
    H_SAMPLE_COL,
    H_STATE_COL,
    as.character(all(SAMPLES %in% unique(m_sample))),
    as.character(all(SAMPLES %in% unique(h_sample))),
    nrow(LRm),
    "FALSE",
    "Descriptive pseudobulk co-expression compatibility; not proof of signaling"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  checks,
  file.path(OUTDIR, "12_analysis_integrity_summary.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("MΦ sample column: ", M_SAMPLE_COL)
message("MΦ final class column: ", M_CLASS_COL)
message("HSC sample column: ", H_SAMPLE_COL)
message("HSC state column: ", H_STATE_COL)
message("Measurable LR pairs: ", nrow(LRm), " / ", nrow(LR))
message("Review first:")
message("  02_sender_Mphi_cell_counts.tsv")
message("  03_receiver_HSC_cell_counts.tsv")
message("  04_LR_pair_coverage.tsv")
message("  08_group_LR_effects.tsv")
message("  09_KEY_Mphi_HSC_LR_effects.tsv")
message("  10_axis_level_summary.tsv")
message("  12_analysis_integrity_summary.tsv")
message("============================================================")
