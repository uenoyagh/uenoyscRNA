#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_PrimaryClean_ShamsTx_Preflight_v1.0.R
#
# PURPOSE
# -------
# Preflight for the first downstream analysis using the finalized Mouse
# macrophage Primary-Clean RDS.
#
# Planned downstream analysis:
#   1) Sham vs Tx M1 / M2 functional axes
#   2) M2/M1 balance
#   3) Resident / Recruited / LAM composition
#
# This script:
# - does NOT modify the RDS
# - does NOT recluster
# - does NOT run statistical testing
# - inventories metadata, sample structure, cluster composition,
#   assay/layers, and marker availability.
#
# Version: v1.0.1
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
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
  "Mouse_PrimaryClean_ShamsTx_Preflight_v1.0.1"
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

if (!file.exists(RDS)) {
  stop("Primary-clean Mouse RDS not found: ", RDS)
}

message("Loading: ", RDS)
obj <- readRDS(RDS)

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
    "RDS",
    "n_cells",
    "n_features_default_assay",
    "default_assay",
    "assays",
    "reductions",
    "graphs"
  ),
  value = c(
    RDS,
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

md <- obj@meta.data

meta_inventory <- do.call(
  rbind,
  lapply(colnames(md), function(cn) {
    x <- md[[cn]]
    ux <- unique(as.character(x))
    data.frame(
      column = cn,
      class = paste(class(x), collapse = ";"),
      n_unique = length(ux),
      example_values = paste(head(ux, 20), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  meta_inventory,
  file.path(OUTDIR, "03_metadata_inventory.tsv")
)

# Candidate annotation / identity columns
candidate_pattern <- paste(
  c(
    "sample", "cluster", "annot", "state", "identity", "mphi", "macro",
    "harmon", "layer", "celltype", "class", "resident", "kupffer",
    "mono", "lam", "trem2", "spp1", "m1", "m2", "qc_"
  ),
  collapse = "|"
)

candidate_cols <- grep(
  candidate_pattern,
  colnames(md),
  value = TRUE,
  ignore.case = TRUE
)

candidate_details <- do.call(
  rbind,
  lapply(candidate_cols, function(cn) {
    tab <- sort(table(as.character(md[[cn]])), decreasing = TRUE)
    data.frame(
      column = cn,
      value = names(tab),
      n = as.integer(tab),
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(
  candidate_details,
  file.path(OUTDIR, "04_candidate_metadata_value_counts.tsv")
)

# ------------------------------------------------------------------------------
# 3. Detect canonical sample / frozen cluster columns
# ------------------------------------------------------------------------------

find_first <- function(candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit)) hit[1] else NA_character_
}

SAMPLE_COL <- find_first(c(
  "sample_display_FIXED2",
  "sample",
  "orig.ident",
  "sample_display"
))

CLUSTER_COL <- find_first(c(
  "cluster_res2",
  "mphi_rpca_res_2.0"
))

if (is.na(SAMPLE_COL)) {
  stop("Could not identify sample column.")
}
if (is.na(CLUSTER_COL)) {
  stop("Could not identify frozen res2 cluster column.")
}

message("Sample column: ", SAMPLE_COL)
message("Cluster column: ", CLUSTER_COL)

sample <- as.character(md[[SAMPLE_COL]])
cluster <- as.character(md[[CLUSTER_COL]])

sample_counts <- as.data.frame(
  table(sample = sample),
  stringsAsFactors = FALSE
)

sample_counts$analysis_group <- ifelse(
  grepl("^Sham", sample_counts$sample, ignore.case = TRUE),
  "Sham",
  ifelse(
    grepl("^Tx", sample_counts$sample, ignore.case = TRUE),
    "Tx",
    ifelse(
      grepl("STD", sample_counts$sample, ignore.case = TRUE),
      "STD",
      ifelse(
        grepl("CDAHFD|CDHFD", sample_counts$sample, ignore.case = TRUE),
        "CDAHFD",
        "Other"
      )
    )
  )
)

write_tsv(
  sample_counts,
  file.path(OUTDIR, "05_sample_counts_and_analysis_groups.tsv")
)

cluster_sample <- as.data.frame(
  table(cluster = cluster, sample = sample),
  stringsAsFactors = FALSE
)
cluster_sample <- cluster_sample[cluster_sample$Freq > 0, , drop = FALSE]

write_tsv(
  cluster_sample,
  file.path(OUTDIR, "06_cluster_x_sample_counts.tsv")
)

# ------------------------------------------------------------------------------
# 4. QC flags sanity check
# ------------------------------------------------------------------------------

qc_cols <- intersect(
  c(
    "qc_original_frozen_cluster",
    "qc_primary_include",
    "qc_sensitivity_include",
    "qc_exclusion_reason",
    "qc_final_disposition",
    "qc_audit_version",
    "qc_resolution_policy",
    "qc_ambient_policy"
  ),
  colnames(md)
)

qc_summary <- if (length(qc_cols)) {
  do.call(
    rbind,
    lapply(qc_cols, function(cn) {
      tab <- sort(table(as.character(md[[cn]]), useNA = "ifany"), decreasing = TRUE)
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
  qc_summary,
  file.path(OUTDIR, "07_qc_flag_summary.tsv")
)

# ------------------------------------------------------------------------------
# 5. Frozen functional / identity marker panels and coverage
# ------------------------------------------------------------------------------

PANELS <- list(
  M1_primary = c("Il1b", "Tnf", "Ccl3", "Ccl4", "Cxcl2"),
  M1_supportive = c("Nfkbia", "Tnfaip3", "Nfkb2", "Nlrp3", "Relb", "Gpr84"),
  M1_IFN = c("Stat1", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3"),
  M1_MASH_pathogenic = c("F13a1", "Olr1", "Selplg"),
  M1_context = c("Ccl2", "Irf5", "Ptgs2", "Cd86", "Cxcl3"),

  M2_primary = c("Mrc1", "Cd163", "Hmox1", "Maf"),
  M2_supportive = c("Dab2", "Msr1", "Stab1", "Il1rn"),
  M2_resolution_efferocytosis = c("Mertk", "Axl", "Gas6", "Mfge8"),
  M2_regulatory = c("Klf4", "Pparg", "Itga4"),
  M2_context = c("Igf1", "Anxa1"),

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

rna_features <- rownames(obj[["RNA"]])

coverage <- do.call(
  rbind,
  lapply(names(PANELS), function(pn) {
    g <- PANELS[[pn]]
    present <- g[g %in% rna_features]
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
  file.path(OUTDIR, "08_marker_panel_coverage.tsv")
)

# ------------------------------------------------------------------------------
# 6. Exact gene-level availability
# ------------------------------------------------------------------------------

all_genes <- unique(unlist(PANELS, use.names = FALSE))

gene_availability <- data.frame(
  gene = all_genes,
  present_in_RNA = all_genes %in% rna_features,
  stringsAsFactors = FALSE
)

write_tsv(
  gene_availability,
  file.path(OUTDIR, "09_marker_gene_availability.tsv")
)

# ------------------------------------------------------------------------------
# 7. Compact analysis-readiness summary
# ------------------------------------------------------------------------------

sham_samples <- unique(sample[grepl("^Sham", sample, ignore.case = TRUE)])
tx_samples <- unique(sample[grepl("^Tx", sample, ignore.case = TRUE)])

readiness <- data.frame(
  item = c(
    "sample_column",
    "cluster_column",
    "n_cells",
    "n_sham_samples",
    "sham_samples",
    "n_tx_samples",
    "tx_samples",
    "n_frozen_clusters",
    "primary_qc_all_TRUE",
    "marker_panels_all_ge_0.5_coverage"
  ),
  value = c(
    SAMPLE_COL,
    CLUSTER_COL,
    ncol(obj),
    length(sham_samples),
    paste(sort(sham_samples), collapse = ";"),
    length(tx_samples),
    paste(sort(tx_samples), collapse = ";"),
    length(unique(cluster)),
    if ("qc_primary_include" %in% colnames(md)) {
      as.character(all(md$qc_primary_include %in% TRUE))
    } else {
      "NA"
    },
    as.character(all(coverage$coverage_fraction >= 0.5))
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  readiness,
  file.path(OUTDIR, "10_analysis_readiness_summary.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Review first:")
message("  04_candidate_metadata_value_counts.tsv")
message("  05_sample_counts_and_analysis_groups.tsv")
message("  06_cluster_x_sample_counts.tsv")
message("  08_marker_panel_coverage.tsv")
message("  10_analysis_readiness_summary.tsv")
message("============================================================")
