#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_FinalQC_Preflight_v1.0.R
#
# PURPOSE
# -------
# Preflight before formal QC of the refined Mouse HSC object.
#
# Audits:
#   1) exact object provenance / assays / layers / reductions / graphs
#   2) sample structure
#   3) current 3-state and 2-state annotation columns
#   4) cluster / resolution metadata candidates
#   5) HSC identity-marker coverage
#   6) foreign-lineage / doublet-marker coverage
#   7) nearby HSC RDS and repo script provenance
#   8) candidate whole-liver background RDS paths for later DecontX sensitivity
#
# IMPORTANT
# ---------
# - Source RDS is NEVER overwritten.
# - No reclustering.
# - No annotation changes.
# - No DecontX is run here.
# - This script only inventories readiness for the next formal audit.
#
# Canonical refined HSC object:
#   v6.1.2
#   13,166 cells expected
#   3-state column expected: HSC_state3_v612
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
REPO <- "/Users/uenoya/Projects/uenoyscRNA"

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
  "Mouse_HSC_FinalQC_Preflight_v1.0"
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

if (!file.exists(HSC_RDS)) {
  stop("HSC RDS not found: ", HSC_RDS)
}

message("Loading HSC object...")
obj <- readRDS(HSC_RDS)
md <- obj@meta.data

# ------------------------------------------------------------------------------
# 1. Object inventory
# ------------------------------------------------------------------------------

assays <- Assays(obj)
reductions <- names(obj@reductions)
graphs <- names(obj@graphs)

assay_layers <- do.call(
  rbind,
  lapply(assays, function(a) {
    lyr <- tryCatch(
      Layers(obj[[a]]),
      error = function(e) character(0)
    )
    if (!length(lyr)) lyr <- NA_character_

    data.frame(
      assay = a,
      layer = lyr,
      stringsAsFactors = FALSE
    )
  })
)

inventory <- data.frame(
  item = c(
    "HSC_RDS",
    "n_cells",
    "n_features_default_assay",
    "default_assay",
    "assays",
    "reductions",
    "graphs"
  ),
  value = c(
    HSC_RDS,
    ncol(obj),
    nrow(obj),
    DefaultAssay(obj),
    paste(assays, collapse = ";"),
    paste(reductions, collapse = ";"),
    paste(graphs, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  inventory,
  file.path(OUTDIR, "01_object_inventory.tsv")
)

write_tsv(
  assay_layers,
  file.path(OUTDIR, "02_assay_layers.tsv")
)

# ------------------------------------------------------------------------------
# 2. Metadata inventory
# ------------------------------------------------------------------------------

meta_inventory <- do.call(
  rbind,
  lapply(colnames(md), function(cn) {
    vals <- unique(as.character(md[[cn]]))

    data.frame(
      column = cn,
      class = paste(class(md[[cn]]), collapse = ";"),
      n_unique = length(vals),
      example_values = paste(head(vals, 40), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  meta_inventory,
  file.path(OUTDIR, "03_metadata_inventory.tsv")
)

# ------------------------------------------------------------------------------
# 3. Detect sample / state / cluster-resolution columns
# ------------------------------------------------------------------------------

find_first <- function(candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit)) hit[1] else NA_character_
}

find_by_values <- function(target_values, preferred = character(0)) {
  for (cn in preferred) {
    if (cn %in% colnames(md)) {
      vals <- unique(as.character(md[[cn]]))
      if (any(target_values %in% vals)) return(cn)
    }
  }

  hits <- character(0)

  for (cn in colnames(md)) {
    vals <- unique(as.character(md[[cn]]))
    if (any(target_values %in% vals)) {
      hits <- c(hits, cn)
    }
  }

  if (!length(hits)) return(NA_character_)
  hits[1]
}

SAMPLE_COL <- find_first(c(
  "sample",
  "sample_display_FIXED2",
  "orig.ident",
  "sample_display"
))

STATE3_COL <- find_first(c(
  "HSC_state3_v612",
  "HSC_state_3state",
  "hsc_state_3state",
  "state_3state"
))

if (is.na(STATE3_COL)) {
  STATE3_COL <- find_by_values(
    c("qHSC", "ECM-activated HSC", "Contractile HSC"),
    preferred = c("HSC_state3_v612")
  )
}

STATE2_COL <- find_first(c(
  "HSC_state2_v612",
  "HSC_state_2state",
  "hsc_state_2state",
  "state_2state"
))

if (is.na(STATE2_COL)) {
  # Search for columns containing both qHSC and aHSC when possible.
  candidate <- character(0)

  for (cn in colnames(md)) {
    vals <- unique(as.character(md[[cn]]))
    if (all(c("qHSC", "aHSC") %in% vals)) {
      candidate <- c(candidate, cn)
    }
  }

  if (length(candidate)) STATE2_COL <- candidate[1]
}

# Cluster / resolution metadata candidates
cluster_candidate_cols <- grep(
  "cluster|res\\.|res_|snn|seurat_clusters|hsc.*res|subcluster",
  colnames(md),
  value = TRUE,
  ignore.case = TRUE
)

cluster_candidate_values <- if (length(cluster_candidate_cols)) {
  do.call(
    rbind,
    lapply(cluster_candidate_cols, function(cn) {
      tab <- sort(
        table(as.character(md[[cn]]), useNA = "ifany"),
        decreasing = TRUE
      )

      data.frame(
        column = cn,
        value = names(tab),
        n = as.integer(tab),
        stringsAsFactors = FALSE
      )
    })
  )
} else {
  data.frame(
    column = character(0),
    value = character(0),
    n = integer(0)
  )
}

write_tsv(
  cluster_candidate_values,
  file.path(OUTDIR, "04_cluster_resolution_metadata_candidates.tsv")
)

detected_cols <- data.frame(
  item = c(
    "sample_col",
    "state3_col",
    "state2_col",
    "n_cluster_resolution_candidate_cols",
    "cluster_resolution_candidate_cols"
  ),
  value = c(
    SAMPLE_COL,
    STATE3_COL,
    STATE2_COL,
    length(cluster_candidate_cols),
    paste(cluster_candidate_cols, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  detected_cols,
  file.path(OUTDIR, "05_detected_core_metadata_columns.tsv")
)

# ------------------------------------------------------------------------------
# 4. Sample / state counts
# ------------------------------------------------------------------------------

if (!is.na(SAMPLE_COL)) {
  sample <- as.character(md[[SAMPLE_COL]])

  sample_counts <- as.data.frame(
    table(sample = sample),
    stringsAsFactors = FALSE
  )

  write_tsv(
    sample_counts,
    file.path(OUTDIR, "06_sample_counts.tsv")
  )
} else {
  sample <- rep(NA_character_, nrow(md))
}

if (!is.na(STATE3_COL)) {
  state3 <- as.character(md[[STATE3_COL]])

  state3_counts <- as.data.frame(
    table(HSC_state3 = state3),
    stringsAsFactors = FALSE
  )

  write_tsv(
    state3_counts,
    file.path(OUTDIR, "07_HSC_3state_counts.tsv")
  )

  if (!is.na(SAMPLE_COL)) {
    state3_sample <- as.data.frame(
      table(
        sample = sample,
        HSC_state3 = state3
      ),
      stringsAsFactors = FALSE
    )

    state3_sample <- state3_sample[
      state3_sample$Freq > 0,
      ,
      drop = FALSE
    ]

    write_tsv(
      state3_sample,
      file.path(OUTDIR, "08_HSC_3state_by_sample_counts.tsv")
    )
  }
}

if (!is.na(STATE2_COL)) {
  state2 <- as.character(md[[STATE2_COL]])

  state2_counts <- as.data.frame(
    table(HSC_state2 = state2),
    stringsAsFactors = FALSE
  )

  write_tsv(
    state2_counts,
    file.path(OUTDIR, "09_HSC_2state_counts.tsv")
  )

  if (!is.na(SAMPLE_COL)) {
    state2_sample <- as.data.frame(
      table(
        sample = sample,
        HSC_state2 = state2
      ),
      stringsAsFactors = FALSE
    )

    state2_sample <- state2_sample[
      state2_sample$Freq > 0,
      ,
      drop = FALSE
    ]

    write_tsv(
      state2_sample,
      file.path(OUTDIR, "10_HSC_2state_by_sample_counts.tsv")
    )
  }
}

# ------------------------------------------------------------------------------
# 5. HSC identity / state marker panels
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

  # Closely related mesenchymal / mural-lineage context
  Pericyte_VSMC = c(
    "Rgs5", "Cspg4", "Pdgfrb", "Des",
    "Acta2", "Tagln", "Myh11", "Notch3"
  ),

  Cycling = c(
    "Mki67", "Top2a", "Cenpf", "Ube2c", "Ccnb1", "Ccnb2"
  )
)

# Foreign-lineage panels for contamination/doublet screening
FOREIGN <- list(
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

all_panels <- c(PANELS, FOREIGN)

rna_features <- rownames(obj[["RNA"]])

coverage <- do.call(
  rbind,
  lapply(names(all_panels), function(pn) {
    g <- all_panels[[pn]]
    present <- intersect(g, rna_features)
    missing <- setdiff(g, present)

    data.frame(
      panel = pn,
      panel_type = ifelse(
        pn %in% names(PANELS),
        "HSC_identity_state",
        "foreign_lineage"
      ),
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
  file.path(OUTDIR, "11_marker_panel_coverage.tsv")
)

# ------------------------------------------------------------------------------
# 6. Nearby RDS / script provenance inventory
# ------------------------------------------------------------------------------

hsc_dir <- dirname(dirname(HSC_RDS))

nearby_rds <- list.files(
  hsc_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

nearby_rds_df <- data.frame(
  path = nearby_rds,
  size_bytes = if (length(nearby_rds)) {
    file.info(nearby_rds)$size
  } else {
    numeric(0)
  },
  stringsAsFactors = FALSE
)

write_tsv(
  nearby_rds_df,
  file.path(OUTDIR, "12_nearby_HSC_RDS_inventory.tsv")
)

repo_files <- list.files(
  REPO,
  pattern = "HSC.*v6\\.1|v6\\.1.*HSC|HSC.*refin|HSC.*annot",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

repo_files_df <- data.frame(
  path = repo_files,
  stringsAsFactors = FALSE
)

write_tsv(
  repo_files_df,
  file.path(OUTDIR, "13_repo_HSC_provenance_candidates.tsv")
)

# Candidate whole-liver background RDS paths.
whole_candidates <- unique(c(
  list.files(
    ROOT,
    pattern = "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  ),
  list.files(
    ROOT,
    pattern = "RH260519343_mo4sample_clustered_annotated\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
))

whole_df <- data.frame(
  path = whole_candidates,
  size_bytes = if (length(whole_candidates)) {
    file.info(whole_candidates)$size
  } else {
    numeric(0)
  },
  stringsAsFactors = FALSE
)

write_tsv(
  whole_df,
  file.path(OUTDIR, "14_candidate_whole_liver_background_RDS.tsv")
)

# ------------------------------------------------------------------------------
# 7. Readiness summary
# ------------------------------------------------------------------------------

state3_ok <- (
  !is.na(STATE3_COL) &&
  all(
    c(
      "qHSC",
      "ECM-activated HSC",
      "Contractile HSC"
    ) %in% unique(as.character(md[[STATE3_COL]]))
  )
)

state2_ok <- (
  !is.na(STATE2_COL) &&
  all(
    c("qHSC", "aHSC") %in%
      unique(as.character(md[[STATE2_COL]]))
  )
)

readiness <- data.frame(
  item = c(
    "expected_cell_count_13166",
    "actual_cell_count",
    "sample_column_found",
    "state3_column_found",
    "state3_expected_labels_present",
    "state2_column_found",
    "state2_qHSC_aHSC_present",
    "n_cluster_resolution_candidate_cols",
    "graphs_present",
    "reductions_present",
    "HSC_core_marker_coverage",
    "foreign_panels_available_for_QC",
    "nearby_HSC_RDS_candidates",
    "repo_HSC_provenance_candidates",
    "whole_liver_background_candidates"
  ),
  value = c(
    as.character(ncol(obj) == 13166),
    ncol(obj),
    as.character(!is.na(SAMPLE_COL)),
    as.character(!is.na(STATE3_COL)),
    as.character(state3_ok),
    as.character(!is.na(STATE2_COL)),
    as.character(state2_ok),
    length(cluster_candidate_cols),
    as.character(length(graphs) > 0),
    as.character(length(reductions) > 0),
    coverage$coverage_fraction[coverage$panel == "HSC_core"][1],
    as.character(all(
      coverage$coverage_fraction[
        coverage$panel_type == "foreign_lineage"
      ] > 0
    )),
    length(nearby_rds),
    length(repo_files),
    length(whole_candidates)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  readiness,
  file.path(OUTDIR, "15_analysis_readiness_summary.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("HSC cells: ", ncol(obj))
message("Sample column: ", SAMPLE_COL)
message("3-state column: ", STATE3_COL)
message("2-state column: ", STATE2_COL)
message("Graphs: ", paste(graphs, collapse = ";"))
message("Reductions: ", paste(reductions, collapse = ";"))
message("Review first:")
message("  01_object_inventory.tsv")
message("  04_cluster_resolution_metadata_candidates.tsv")
message("  05_detected_core_metadata_columns.tsv")
message("  06_sample_counts.tsv")
message("  07_HSC_3state_counts.tsv")
message("  08_HSC_3state_by_sample_counts.tsv")
message("  09_HSC_2state_counts.tsv")
message("  11_marker_panel_coverage.tsv")
message("  13_repo_HSC_provenance_candidates.tsv")
message("  14_candidate_whole_liver_background_RDS.tsv")
message("  15_analysis_readiness_summary.tsv")
message("============================================================")
