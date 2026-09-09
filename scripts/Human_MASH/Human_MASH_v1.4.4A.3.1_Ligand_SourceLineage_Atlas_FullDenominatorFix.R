#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.3.1_Ligand_SourceLineage_Atlas_FullDenominatorFix.R
#
# Purpose
#   Corrected source-lineage atlas for SfLB1/SfLB2 vs Tx17/Tx5.
#
# Critical technical correction relative to v1.4.4A.3
#   In B.2C.3 lineage_reconstruction_long.tsv, n_cells is the number of canonical
#   cells that had >=1 selected molecule among the prioritized ligand panel
#   ("candidate-signal cell universe"), NOT the full frozen lineage cell count.
#
#   Therefore:
#     - UMI numerator               : B.2C.3 molecule reconstruction (valid)
#     - positive_cells numerator    : B.2C.3 molecule reconstruction (valid)
#     - lineage denominator n_cells : v1.4.4A.1.1 frozen full lineage cell count
#
#   This restores true zero-expressing cells to the denominator and avoids
#   inflation of UMI/1000 cells and positive fraction.
#
# Biological groups
#   SfLB1/SfLB2 = pre-transplant, non-engineered cultured organoids
#   Tx17/Tx5    = engrafted human cells derived from IL10-engineered organoids
#
# Primary evidence
#   - Non-IL10: SHARED_FILTERED_TRUE_ALL4
#   - IL10: ALL_PROBES for primary lineage-source localization because it matches
#     standard MEX; SHARED_FILTERED_TRUE_ALL4 retained as stringent sensitivity.
#
# Causal boundary
#   SfLB -> Tx is NOT an IL10-only perturbation.
#   It combines transgene status, transplantation/engraftment, host microenvironment,
#   post-engraftment selection/state transition, and possible residual innate
#   xenogeneic immune pressure in NOD/scid recipients.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

B2C3_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3B.2C.3/",
  "analysis/v1.4.3B.2C.3_prioritized_ligand_molecule_reconstruction/tables"
)

LINEAGE_TSV <- file.path(
  B2C3_DIR,
  "v1.4.3B.2C.3_lineage_reconstruction_long.tsv"
)

BARCODE_AUDIT_TSV <- file.path(
  B2C3_DIR,
  "v1.4.3B.2C.3_barcode_mapping_audit.tsv"
)

A11_TSV <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.1.1/",
  "analysis/v1.4.4A.1.1_pretransplant_SfLB_vs_engraftedTx_lineage_composition/",
  "tables/02_broad_lineage_counts_fractions_by_sample_v1.4.4A.1.1.csv"
)

A2_CLASS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.2/",
  "analysis/v1.4.4A.2_probe_aware_ligand_transition_audit/",
  "tables/03_ligand_transition_classification_v1.4.4A.2.csv"
)

OUT_ROOT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3.1/",
  "analysis/v1.4.4A.3.1_ligand_source_lineage_atlas_full_denominator_fix"
)

TAB_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR <- file.path(OUT_ROOT, "figures")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "v1.4.4A.3.1_run_log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================\n")
cat("Human MASH v1.4.4A.3.1\n")
cat("Ligand source-lineage atlas — full denominator correction\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Lineage reconstruction:", LINEAGE_TSV, "\n")
cat("Barcode audit:", BARCODE_AUDIT_TSV, "\n")
cat("A.1.1 composition:", A11_TSV, "\n")
cat("A.2 classification:", A2_CLASS, "\n")
cat("Output:", OUT_ROOT, "\n\n")

for (f in c(LINEAGE_TSV, BARCODE_AUDIT_TSV, A11_TSV, A2_CLASS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

# ------------------------------------------------------------------------------
# 1. Read inputs
# ------------------------------------------------------------------------------

lin <- read.delim(
  LINEAGE_TSV,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expected_cols <- c(
  "sample", "gene", "mode", "group", "n_cells",
  "UMI", "positive_cells", "positive_fraction", "UMI_per_1000_cells"
)

if (!identical(colnames(lin), expected_cols)) {
  stop(
    "Unexpected lineage TSV schema.\nObserved: ",
    paste(colnames(lin), collapse = " | ")
  )
}

if (nrow(lin) != 1620L) {
  stop("Unexpected lineage reconstruction row count: ", nrow(lin), " (expected 1620)")
}

mode_tab <- table(lin$mode)
if (length(mode_tab) != 6L || any(mode_tab != 270L)) {
  stop("Unexpected mode distribution: ",
       paste(names(mode_tab), as.integer(mode_tab), collapse = "; "))
}

audit <- read.delim(
  BARCODE_AUDIT_TSV,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

a11 <- read.csv(A11_TSV, stringsAsFactors = FALSE, check.names = FALSE)
a2  <- read.csv(A2_CLASS, stringsAsFactors = FALSE, check.names = FALSE)

candidate_order <- c(
  "IL10",
  "FGF1",
  "TNFSF12",
  "VEGFA",
  "CXCL12",
  "TGFB3",
  "BMP5",
  "BMP2",
  "HBEGF",
  "BTC",
  "FST",
  "IGFBP2",
  "IGFBP4",
  "INHBA",
  "VEGFB"
)

sample_order <- c("SfLB1", "SfLB2", "Tx17", "Tx5")

all_lineages <- c(
  "Hepatocyte_like",
  "Mesenchymal_like",
  "Endothelial_like",
  "Hepatobiliary_like",
  "Myeloid_like",
  "Unresolved"
)

if (!setequal(unique(lin$gene), candidate_order)) {
  stop("Candidate gene set mismatch in lineage reconstruction.")
}
if (!setequal(unique(lin$sample), sample_order)) {
  stop("Sample set mismatch in lineage reconstruction.")
}

needed_a11 <- c("sample", "category", "n_cells", "sample_total", "fraction")
if (!all(needed_a11 %in% colnames(a11))) {
  stop("A.1.1 table is missing required columns.")
}

a11$sample <- as.character(a11$sample)
a11$lineage <- as.character(a11$category)

# ------------------------------------------------------------------------------
# 2. Prove the B.2C.3 denominator provenance
# ------------------------------------------------------------------------------

# B.2C.3 n_cells is invariant across genes/modes for a sample-lineage pair.
ncheck <- unique(lin[, c("sample", "group", "n_cells")])
names(ncheck)[2] <- "lineage"
names(ncheck)[3] <- "B2C3_candidate_signal_n_cells"

pair_key <- paste(ncheck$sample, ncheck$lineage, sep = "||")
if (anyDuplicated(pair_key)) {
  stop("B.2C.3 contains >1 n_cells value for the same sample-lineage pair.")
}

# Sum B.2C.3 lineage n_cells by sample.
b2_sum <- aggregate(
  B2C3_candidate_signal_n_cells ~ sample,
  ncheck,
  sum
)

needed_audit <- c(
  "sample",
  "canonical_cells_expected",
  "canonical_barcode_indices_with_candidate_signal"
)
if (!all(needed_audit %in% colnames(audit))) {
  stop("Barcode audit missing required denominator-provenance columns.")
}

audit_small <- audit[, needed_audit]

prov <- merge(
  audit_small,
  b2_sum,
  by = "sample",
  all = TRUE,
  sort = FALSE
)

prov <- prov[match(sample_order, prov$sample), ]

prov$B2C3_lineage_sum_matches_candidate_signal_universe <-
  prov$B2C3_candidate_signal_n_cells ==
  prov$canonical_barcode_indices_with_candidate_signal

prov$candidate_signal_fraction_of_full_canonical <-
  prov$canonical_barcode_indices_with_candidate_signal /
  prov$canonical_cells_expected

if (!all(prov$B2C3_lineage_sum_matches_candidate_signal_universe)) {
  stop(
    "DENOMINATOR PROVENANCE FAIL: B.2C.3 lineage n_cells do not sum to ",
    "canonical_barcode_indices_with_candidate_signal."
  )
}

cat("DENOMINATOR PROVENANCE: PASS\n")
cat("B.2C.3 lineage n_cells == candidate-signal cell universe, not full canonical universe.\n")
print(prov, row.names = FALSE)

write.csv(
  prov,
  file.path(TAB_DIR, "00_denominator_provenance_validation_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Build full frozen denominator from A.1.1
# ------------------------------------------------------------------------------

grid_sl <- expand.grid(
  sample = sample_order,
  lineage = all_lineages,
  stringsAsFactors = FALSE
)

a11_small <- a11[, c("sample", "lineage", "n_cells", "sample_total", "fraction")]

full_den <- merge(
  grid_sl,
  a11_small,
  by = c("sample", "lineage"),
  all.x = TRUE,
  sort = FALSE
)

# A.1.1 is zero-completed, but keep this defensive.
full_den$n_cells[is.na(full_den$n_cells)] <- 0
full_den$fraction[is.na(full_den$fraction)] <- 0

sample_totals <- unique(a11[, c("sample", "sample_total")])

# Validate full sample totals.
full_sum <- aggregate(n_cells ~ sample, full_den, sum)
full_sum <- merge(full_sum, sample_totals, by = "sample", all = TRUE, sort = FALSE)
full_sum <- full_sum[match(sample_order, full_sum$sample), ]

if (any(full_sum$n_cells != full_sum$sample_total)) {
  stop("A.1.1 full lineage denominators do not sum to sample totals.")
}

cat("\nA.1.1 FULL FROZEN DENOMINATOR: PASS\n")
print(full_sum, row.names = FALSE)

# Reconciliation table
recon <- merge(
  full_den[, c("sample", "lineage", "n_cells")],
  ncheck,
  by = c("sample", "lineage"),
  all.x = TRUE,
  sort = FALSE
)

recon$B2C3_candidate_signal_n_cells[
  is.na(recon$B2C3_candidate_signal_n_cells)
] <- 0

names(recon)[names(recon) == "n_cells"] <- "A11_full_n_cells"
recon$zero_ligand_panel_cells_restored <-
  recon$A11_full_n_cells - recon$B2C3_candidate_signal_n_cells

recon$B2C3_candidate_signal_fraction_of_full <- ifelse(
  recon$A11_full_n_cells > 0,
  recon$B2C3_candidate_signal_n_cells / recon$A11_full_n_cells,
  NA_real_
)

write.csv(
  recon,
  file.path(TAB_DIR, "01_lineage_denominator_reconciliation_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Lineage composition status
# ------------------------------------------------------------------------------

full_den$bio_group <- ifelse(
  full_den$sample %in% c("SfLB1", "SfLB2"),
  "SfLB",
  "Tx"
)

lineage_group_mean <- aggregate(
  fraction ~ lineage + bio_group,
  full_den,
  mean
)

get_group_mean <- function(lineage_name, grp) {
  z <- lineage_group_mean[
    lineage_group_mean$lineage == lineage_name &
      lineage_group_mean$bio_group == grp,
    "fraction"
  ]
  if (length(z) == 0) return(0)
  z[1]
}

lineage_status <- do.call(
  rbind,
  lapply(all_lineages, function(ll) {
    s <- get_group_mean(ll, "SfLB")
    t <- get_group_mean(ll, "Tx")

    status <- if (s >= 0.05 && t >= 0.05) {
      "Comparable_major_lineage"
    } else if (s >= 0.05 && t < 0.001) {
      "SfLB_enriched_essentially_lost_after_engraftment"
    } else if (s == 0 && t >= 0.05) {
      "Tx_emergent_lineage"
    } else {
      "Minor_or_partial_lineage"
    }

    data.frame(
      lineage = ll,
      SfLB_mean_fraction = s,
      Tx_mean_fraction = t,
      Tx_minus_SfLB_fraction = t - s,
      lineage_composition_status = status,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  lineage_status,
  file.path(TAB_DIR, "02_lineage_composition_status_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

cat("\n=== LINEAGE COMPOSITION STATUS ===\n")
print(lineage_status, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Primary non-IL10 strict shared-TRUE source atlas
# ------------------------------------------------------------------------------

strict_mode <- "SHARED_FILTERED_TRUE_ALL4"

strict <- lin[
  lin$mode == strict_mode & lin$gene != "IL10",
  c("sample", "gene", "group", "UMI", "positive_cells"),
  drop = FALSE
]
names(strict)[names(strict) == "group"] <- "lineage"

# Zero-complete all non-IL10 gene x sample x lineage combinations.
grid <- expand.grid(
  gene = setdiff(candidate_order, "IL10"),
  sample = sample_order,
  lineage = all_lineages,
  stringsAsFactors = FALSE
)

atlas <- merge(
  grid,
  strict,
  by = c("gene", "sample", "lineage"),
  all.x = TRUE,
  sort = FALSE
)

atlas$UMI[is.na(atlas$UMI)] <- 0
atlas$positive_cells[is.na(atlas$positive_cells)] <- 0

# IMPORTANT FIX: full frozen denominator from A.1.1
atlas <- merge(
  atlas,
  full_den[, c("sample", "lineage", "n_cells", "fraction", "sample_total")],
  by = c("sample", "lineage"),
  all.x = TRUE,
  sort = FALSE
)

atlas$full_lineage_n_cells <- atlas$n_cells
atlas$lineage_fraction_of_sample <- atlas$fraction
atlas$n_cells <- NULL
atlas$fraction <- NULL

atlas$positive_fraction_full_denominator <- ifelse(
  atlas$full_lineage_n_cells > 0,
  atlas$positive_cells / atlas$full_lineage_n_cells,
  NA_real_
)

atlas$UMI_per_1000_full_lineage_cells <- ifelse(
  atlas$full_lineage_n_cells > 0,
  1000 * atlas$UMI / atlas$full_lineage_n_cells,
  NA_real_
)

atlas$bio_group <- ifelse(
  atlas$sample %in% c("SfLB1", "SfLB2"),
  "SfLB",
  "Tx"
)

atlas <- merge(
  atlas,
  lineage_status[, c("lineage", "lineage_composition_status")],
  by = "lineage",
  all.x = TRUE,
  sort = FALSE
)

# Gene-sample total UMI is denominator-independent.
gene_sample_tot <- aggregate(
  UMI ~ gene + sample,
  atlas,
  sum
)
names(gene_sample_tot)[3] <- "gene_sample_total_UMI"

atlas <- merge(
  atlas,
  gene_sample_tot,
  by = c("gene", "sample"),
  all.x = TRUE,
  sort = FALSE
)

atlas$source_UMI_fraction <- ifelse(
  atlas$gene_sample_total_UMI > 0,
  atlas$UMI / atlas$gene_sample_total_UMI,
  NA_real_
)

write.csv(
  atlas,
  file.path(TAB_DIR, "03_nonIL10_strict_shared_TRUE_source_lineage_atlas_FULL_DENOM_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Dominant source / source shift
# ------------------------------------------------------------------------------

mean_source <- aggregate(
  source_UMI_fraction ~ gene + bio_group + lineage,
  atlas,
  mean,
  na.rm = TRUE
)
mean_source$source_UMI_fraction[is.nan(mean_source$source_UMI_fraction)] <- NA_real_

dominant_source_for <- function(g, grp) {
  z <- mean_source[mean_source$gene == g & mean_source$bio_group == grp, ]
  z <- z[!is.na(z$source_UMI_fraction), ]
  if (nrow(z) == 0) return(c(lineage = NA, fraction = NA))
  z <- z[order(z$source_UMI_fraction, decreasing = TRUE), ]
  c(lineage = z$lineage[1], fraction = z$source_UMI_fraction[1])
}

emergent_lineages <- lineage_status$lineage[
  lineage_status$lineage_composition_status == "Tx_emergent_lineage"
]
lost_lineages <- lineage_status$lineage[
  lineage_status$lineage_composition_status ==
    "SfLB_enriched_essentially_lost_after_engraftment"
]
comparable_lineages <- lineage_status$lineage[
  lineage_status$lineage_composition_status == "Comparable_major_lineage"
]

gene_source_summary <- do.call(
  rbind,
  lapply(setdiff(candidate_order, "IL10"), function(g) {
    ds <- dominant_source_for(g, "SfLB")
    dt <- dominant_source_for(g, "Tx")

    ms <- mean_source[mean_source$gene == g, ]

    tx_emergent_share <- sum(
      ms$source_UMI_fraction[
        ms$bio_group == "Tx" & ms$lineage %in% emergent_lineages
      ],
      na.rm = TRUE
    )

    sflb_lost_share <- sum(
      ms$source_UMI_fraction[
        ms$bio_group == "SfLB" & ms$lineage %in% lost_lineages
      ],
      na.rm = TRUE
    )

    data.frame(
      gene = g,
      SfLB_dominant_source = as.character(ds["lineage"]),
      SfLB_dominant_source_fraction = as.numeric(ds["fraction"]),
      Tx_dominant_source = as.character(dt["lineage"]),
      Tx_dominant_source_fraction = as.numeric(dt["fraction"]),
      dominant_source_switch =
        !is.na(ds["lineage"]) &&
        !is.na(dt["lineage"]) &&
        as.character(ds["lineage"]) != as.character(dt["lineage"]),
      Tx_mean_source_fraction_from_Tx_emergent_lineages = tx_emergent_share,
      SfLB_mean_source_fraction_from_lineages_lost_after_engraftment = sflb_lost_share,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  gene_source_summary,
  file.path(TAB_DIR, "04_nonIL10_dominant_source_and_source_shift_summary_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Within-lineage replicate-aware transition using FULL denominator
# ------------------------------------------------------------------------------

compare_one_gene_lineage <- function(g, ll) {
  z <- atlas[atlas$gene == g & atlas$lineage == ll, ]

  rr <- setNames(
    z$UMI_per_1000_full_lineage_cells,
    z$sample
  )

  s <- rr[c("SfLB1", "SfLB2")]
  t <- rr[c("Tx17", "Tx5")]

  if (anyNA(c(s, t))) {
    rel <- "Not_comparable_due_to_absent_lineage"
  } else if (min(t) > max(s)) {
    rel <- "Tx_higher_complete_replicate_range_separation"
  } else if (max(t) < min(s)) {
    rel <- "Tx_lower_complete_replicate_range_separation"
  } else {
    rel <- "Replicate_ranges_overlap"
  }

  smean <- mean(s, na.rm = TRUE)
  tmean <- mean(t, na.rm = TRUE)

  data.frame(
    gene = g,
    lineage = ll,
    SfLB1_UMI_per_1000_FULL = s["SfLB1"],
    SfLB2_UMI_per_1000_FULL = s["SfLB2"],
    Tx17_UMI_per_1000_FULL = t["Tx17"],
    Tx5_UMI_per_1000_FULL = t["Tx5"],
    SfLB_mean_UMI_per_1000_FULL = smean,
    Tx_mean_UMI_per_1000_FULL = tmean,
    Tx_minus_SfLB_mean_UMI_per_1000_FULL = tmean - smean,
    Tx_vs_SfLB_log2FC_pseudocount =
      log2((tmean + 0.01) / (smean + 0.01)),
    replicate_relation = rel,
    denominator_note = "FULL frozen lineage cells from v1.4.4A.1.1",
    stringsAsFactors = FALSE
  )
}

within_lineage <- do.call(
  rbind,
  lapply(setdiff(candidate_order, "IL10"), function(g) {
    do.call(
      rbind,
      lapply(comparable_lineages, function(ll) {
        compare_one_gene_lineage(g, ll)
      })
    )
  })
)

write.csv(
  within_lineage,
  file.path(TAB_DIR, "05_nonIL10_within_comparable_lineage_transition_FULL_DENOM_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Integrated source-lineage interpretation
# ------------------------------------------------------------------------------

a2_non <- a2[
  a2$gene != "IL10",
  c(
    "gene",
    "Tx_vs_SfLB_log2FC_per_cell_descriptive",
    "replicate_relation",
    "transition_call"
  )
]

integrated <- merge(
  a2_non,
  gene_source_summary,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

wl_counts <- do.call(
  rbind,
  lapply(setdiff(candidate_order, "IL10"), function(g) {
    z <- within_lineage[within_lineage$gene == g, ]
    data.frame(
      gene = g,
      n_comparable_major_lineages = nrow(z),
      n_comparable_lineages_Tx_higher_complete =
        sum(z$replicate_relation ==
              "Tx_higher_complete_replicate_range_separation"),
      n_comparable_lineages_Tx_lower_complete =
        sum(z$replicate_relation ==
              "Tx_lower_complete_replicate_range_separation"),
      n_comparable_lineages_overlap =
        sum(z$replicate_relation == "Replicate_ranges_overlap"),
      comparable_lineage_detail =
        paste(
          paste0(z$lineage, ":", z$replicate_relation),
          collapse = ";"
        ),
      stringsAsFactors = FALSE
    )
  })
)

integrated <- merge(
  integrated,
  wl_counts,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

integrated$mechanistic_interpretation <- apply(
  integrated,
  1,
  function(r) {
    whole_rel <- r[["replicate_relation"]]
    n_up <- as.numeric(r[["n_comparable_lineages_Tx_higher_complete"]])
    emerg <- as.numeric(r[["Tx_mean_source_fraction_from_Tx_emergent_lineages"]])
    sw <- as.character(r[["dominant_source_switch"]]) == "TRUE"

    if (whole_rel == "Tx_higher_complete_separation" && n_up >= 1 && emerg >= 0.20) {
      "Whole-graft Tx-higher with both within-lineage and Tx-emergent-source contributions"
    } else if (whole_rel == "Tx_higher_complete_separation" && n_up >= 1) {
      "Whole-graft Tx-higher with supporting within-lineage Tx-higher evidence"
    } else if (whole_rel == "Tx_higher_complete_separation" && emerg >= 0.20) {
      "Whole-graft Tx-higher strongly influenced by Tx-emergent source lineages"
    } else if (whole_rel == "Tx_higher_complete_separation" && sw) {
      "Whole-graft Tx-higher with dominant source-lineage switch"
    } else if (whole_rel == "Tx_higher_complete_separation") {
      "Whole-graft Tx-higher; source decomposition mixed"
    } else if (whole_rel == "Replicate_ranges_overlap" && n_up >= 1) {
      "Whole-graft ranges overlap despite supporting within-lineage Tx-higher evidence"
    } else {
      "No strong whole-graft Tx-higher closure; inspect source atlas"
    }
  }
)

write.csv(
  integrated,
  file.path(TAB_DIR, "06_nonIL10_integrated_source_lineage_interpretation_FULL_DENOM_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

cat("\n=== INTEGRATED NON-IL10 SOURCE-LINEAGE INTERPRETATION ===\n")
print(
  integrated[, c(
    "gene",
    "Tx_vs_SfLB_log2FC_per_cell_descriptive",
    "replicate_relation",
    "SfLB_dominant_source",
    "Tx_dominant_source",
    "dominant_source_switch",
    "Tx_mean_source_fraction_from_Tx_emergent_lineages",
    "n_comparable_lineages_Tx_higher_complete",
    "mechanistic_interpretation"
  )],
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. IL10 source-lineage resolution with FULL denominator
# ------------------------------------------------------------------------------

make_il10_mode <- function(mode_name) {
  z <- lin[
    lin$gene == "IL10" & lin$mode == mode_name,
    c("sample", "group", "UMI", "positive_cells"),
    drop = FALSE
  ]
  names(z)[names(z) == "group"] <- "lineage"

  grid_il <- expand.grid(
    sample = sample_order,
    lineage = all_lineages,
    stringsAsFactors = FALSE
  )

  z <- merge(
    grid_il,
    z,
    by = c("sample", "lineage"),
    all.x = TRUE,
    sort = FALSE
  )

  z$UMI[is.na(z$UMI)] <- 0
  z$positive_cells[is.na(z$positive_cells)] <- 0

  z <- merge(
    z,
    full_den[, c("sample", "lineage", "n_cells")],
    by = c("sample", "lineage"),
    all.x = TRUE,
    sort = FALSE
  )

  names(z)[names(z) == "n_cells"] <- "full_lineage_n_cells"

  z$positive_fraction_full_denominator <- ifelse(
    z$full_lineage_n_cells > 0,
    z$positive_cells / z$full_lineage_n_cells,
    NA_real_
  )

  z$UMI_per_1000_full_lineage_cells <- ifelse(
    z$full_lineage_n_cells > 0,
    1000 * z$UMI / z$full_lineage_n_cells,
    NA_real_
  )

  z$mode <- mode_name

  tot <- aggregate(UMI ~ sample, z, sum)
  names(tot)[2] <- "sample_total_IL10_UMI"

  z <- merge(z, tot, by = "sample", all.x = TRUE, sort = FALSE)

  z$source_UMI_fraction <- ifelse(
    z$sample_total_IL10_UMI > 0,
    z$UMI / z$sample_total_IL10_UMI,
    NA_real_
  )

  z
}

il10_all <- make_il10_mode("ALL_PROBES")
il10_strict <- make_il10_mode("SHARED_FILTERED_TRUE_ALL4")
il10_source <- rbind(il10_all, il10_strict)

write.csv(
  il10_source,
  file.path(TAB_DIR, "07_IL10_source_lineage_FULL_DENOM_ALL_PROBES_and_strict_v1.4.4A.3.1.csv"),
  row.names = FALSE
)

cat("\n=== IL10 SOURCE-LINEAGE, ALL_PROBES, FULL DENOMINATOR ===\n")
print(
  il10_all[
    il10_all$UMI > 0,
    c(
      "sample", "lineage", "full_lineage_n_cells",
      "UMI", "positive_cells",
      "UMI_per_1000_full_lineage_cells",
      "positive_fraction_full_denominator",
      "source_UMI_fraction"
    )
  ],
  row.names = FALSE
)

cat("\n=== IL10 SOURCE-LINEAGE, STRICT SHARED-TRUE, FULL DENOMINATOR ===\n")
print(
  il10_strict[
    il10_strict$UMI > 0,
    c(
      "sample", "lineage", "full_lineage_n_cells",
      "UMI", "positive_cells",
      "UMI_per_1000_full_lineage_cells",
      "positive_fraction_full_denominator",
      "source_UMI_fraction"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------------------------

# Fig 1: source UMI fraction atlas (denominator-independent source fraction)
sf <- atlas
sf$sample <- factor(sf$sample, levels = sample_order)
sf$lineage <- factor(sf$lineage, levels = all_lineages)

gene_plot_order <- a2_non$gene[
  order(a2_non$Tx_vs_SfLB_log2FC_per_cell_descriptive, decreasing = TRUE)
]
sf$gene <- factor(sf$gene, levels = rev(gene_plot_order))
sf$xlab <- interaction(sf$sample, sf$lineage, sep = " | ", lex.order = TRUE)

p1 <- ggplot(sf, aes(x = xlab, y = gene, fill = source_UMI_fraction)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A",
    limits = c(0, 1),
    na.value = "#F0F0F0"
  ) +
  labs(
    title = "Ligand source-lineage atlas",
    subtitle = "Strict shared-TRUE mode; fill = fraction of each gene's sample-level UMI contributed by each lineage",
    x = "Sample | lineage",
    y = NULL,
    fill = "Source UMI fraction"
  ) +
  theme_minimal(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7.2),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig01_nonIL10_source_fraction_atlas_v1.4.4A.3.1.pdf"),
  p1, width = 15.5, height = 7.5, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig01_nonIL10_source_fraction_atlas_v1.4.4A.3.1.png"),
  p1, width = 15.5, height = 7.5, dpi = 300
)

# Fig 2: full-denominator within-lineage abundance in comparable major lineages
wl_plot <- atlas[atlas$lineage %in% comparable_lineages, ]
wl_plot$sample <- factor(wl_plot$sample, levels = sample_order)
wl_plot$gene <- factor(wl_plot$gene, levels = rev(gene_plot_order))

wl_plot$display_group <- factor(
  ifelse(
    wl_plot$sample %in% c("SfLB1", "SfLB2"),
    "Pre-transplant SfLB",
    "Engrafted Tx"
  ),
  levels = c("Pre-transplant SfLB", "Engrafted Tx")
)

p2 <- ggplot(
  wl_plot,
  aes(x = display_group, y = UMI_per_1000_full_lineage_cells)
) +
  geom_point(
    aes(shape = sample),
    size = 2.1,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.38,
    linewidth = 0.45
  ) +
  facet_grid(gene ~ lineage, scales = "free_y") +
  labs(
    title = "Within-lineage ligand abundance using full frozen denominators",
    subtitle = "Strict shared-TRUE UMI per 1,000 lineage cells; descriptive n=2 vs n=2",
    x = NULL,
    y = "UMI per 1,000 full lineage cells",
    shape = "Sample"
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 7.5),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig02_nonIL10_within_lineage_FULL_DENOM_v1.4.4A.3.1.pdf"),
  p2, width = 10.5, height = 20.0, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig02_nonIL10_within_lineage_FULL_DENOM_v1.4.4A.3.1.png"),
  p2, width = 10.5, height = 20.0, dpi = 250
)

# Fig 3: mean source composition
ms <- mean_source
ms$bio_group <- factor(ms$bio_group, levels = c("SfLB", "Tx"))
ms$lineage <- factor(ms$lineage, levels = all_lineages)
ms$gene <- factor(ms$gene, levels = rev(gene_plot_order))

p3 <- ggplot(
  ms,
  aes(x = bio_group, y = source_UMI_fraction, fill = lineage)
) +
  geom_col(width = 0.72) +
  facet_wrap(~ gene, ncol = 4) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = function(x) paste0(round(100*x), "%")
  ) +
  labs(
    title = "Mean ligand source composition before and after engraftment",
    subtitle = "Mean of two biological samples per group; source fractions are denominator-independent",
    x = NULL,
    y = "Mean source UMI fraction",
    fill = "Lineage"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 7.5),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig03_nonIL10_mean_source_composition_v1.4.4A.3.1.pdf"),
  p3, width = 12.5, height = 10.5, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig03_nonIL10_mean_source_composition_v1.4.4A.3.1.png"),
  p3, width = 12.5, height = 10.5, dpi = 300
)

# Fig 4: IL10 source-lineage, ALL_PROBES, full denominator
ilp <- il10_all
ilp$sample <- factor(ilp$sample, levels = sample_order)
ilp$lineage <- factor(ilp$lineage, levels = all_lineages)

p4 <- ggplot(
  ilp,
  aes(x = sample, y = UMI, fill = lineage)
) +
  geom_col(width = 0.72) +
  labs(
    title = "IL10 source-lineage resolution",
    subtitle = "ALL_PROBES mode matching standard MEX; full frozen lineage denominators used for rate calculations",
    x = NULL,
    y = "IL10 UMI",
    fill = "Lineage"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig04_IL10_source_lineage_FULL_DENOM_v1.4.4A.3.1.pdf"),
  p4, width = 8.0, height = 5.8, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig04_IL10_source_lineage_FULL_DENOM_v1.4.4A.3.1.png"),
  p4, width = 8.0, height = 5.8, dpi = 300
)

# ------------------------------------------------------------------------------
# 11. Interpretation note
# ------------------------------------------------------------------------------

note <- c(
  "Human MASH v1.4.4A.3.1 interpretation notes",
  "",
  "Technical correction:",
  "- B.2C.3 lineage n_cells equals the canonical cells with >=1 selected molecule",
  "  among the prioritized ligand panel (candidate-signal cell universe).",
  "- It is therefore not the correct denominator for lineage-level expression rates.",
  "- v1.4.4A.3.1 restores all zero-expressing frozen lineage cells using A.1.1 denominators.",
  "- B.2C.3 UMI and positive-cell numerators remain valid.",
  "",
  "Evidence hierarchy:",
  "- Non-IL10: SHARED_FILTERED_TRUE_ALL4 is primary.",
  "- IL10: ALL_PROBES is primary for lineage-source localization because it reproduces",
  "  standard MEX; SHARED_FILTERED_TRUE_ALL4 is a stringent sensitivity analysis.",
  "",
  "Critical causal boundary:",
  "- SfLB -> Tx is NOT an IL10-only perturbation.",
  "- It combines transgene status, transplantation/engraftment, host microenvironment,",
  "  post-engraftment selection/state transition, and possible residual innate",
  "  xenogeneic immune pressure.",
  "- Within-lineage Tx-higher abundance does not prove direct IL10 transcriptional induction.",
  "",
  "Paper-safe language:",
  "- 'Engraftment-associated source-lineage remodeling' is preferred.",
  "- 'Within-lineage Tx-higher abundance' is preferred over 'IL10-induced upregulation'.",
  "- Tx-emergent Myeloid-like/Unresolved sources require cautious biological interpretation.",
  "",
  "Next step:",
  "- Review A.3.1 integrated source-lineage table.",
  "- Then perform v1.4.4A.3B graft stress / residual xenogeneic immune-pressure audit.",
  "- Afterwards proceed to measurable hepatocyte-function audit and LT7/LT8 reference comparison."
)

writeLines(
  note,
  file.path(OUT_ROOT, "README_interpretation_v1.4.4A.3.1.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("v1.4.4A.3.1 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_ROOT, "\n")
cat("Tables:", length(list.files(TAB_DIR)), "\n")
cat("Figures:", length(list.files(FIG_DIR)), "\n")
cat("DENOMINATOR PROVENANCE: PASS\n")
cat("Full frozen lineage denominator restored from v1.4.4A.1.1.\n")
cat("B.2C.3 UMI / positive-cell numerators retained.\n")
cat("No inferential test performed (n=2 per group).\n")
