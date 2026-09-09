#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.3B.1P.3_StrictSharedProbe_GeneDriverConcordanceAudit.R
#
# Purpose
# -------
# Lightweight post-reconstruction audit. No HDF5 scan is performed.
#
# Uses:
#   - v1.4.4A.3B.0 measured-gene panel definitions
#   - v1.4.4A.3B.1P.1 strict cross-sample probe comparability
#   - v1.4.4A.3B.1P.2 full-denominator strict shared-probe reconstruction
#
# For Hepatocyte_like and Mesenchymal_like, this script:
#   1) identifies which strict-comparable genes drive each program,
#   2) counts Tx-higher / Tx-lower / overlapping 2-vs-2 replicate relations,
#   3) compares the direction from standard gene-level matrix metrics with
#      the strict shared-probe molecule reconstruction,
#   4) reports antigen-presentation gene behavior separately.
#
# IMPORTANT
# -------
# - No new biological module score is created yet.
# - No inferential test is performed (n=2/group).
# - Program-level direction is not automatically assigned from an arbitrary
#   majority threshold; gene-level evidence is shown explicitly.
# ==============================================================================

options(stringsAsFactors = FALSE)

BASE <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

A3B0 <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.0",
  "analysis/v1.4.4A.3B.0_graft_stress_xenogeneic_pressure_feature_preflight",
  "tables/01_candidate_gene_feature_coverage_long_v1.4.4A.3B.0.csv"
)

P1 <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.1P.1",
  "analysis/v1.4.4A.3B.1P.1_measured_gene_probe_status_audit",
  "v1.4.4A.3B.1P.1_measured_gene_crosssample_probe_summary.tsv"
)

P2 <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.1P.2",
  "analysis/v1.4.4A.3B.1P.2_strict_shared_probe_molecule_reconstruction",
  "tables/v1.4.4A.3B.1P.2_strict_lineage_reconstruction_FULL_DENOM.tsv"
)

OUT <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.1P.3",
  "analysis/v1.4.4A.3B.1P.3_strict_shared_probe_gene_driver_concordance_audit"
)

TAB <- file.path(OUT, "tables")
LOG <- file.path(OUT, "logs")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

sink(file.path(LOG, "v1.4.4A.3B.1P.3_run_log.txt"), split = TRUE)
on.exit(sink(), add = TRUE)

cat("====================================================================\n")
cat("Human MASH v1.4.4A.3B.1P.3\n")
cat("Strict shared-probe gene-driver / concordance audit\n")
cat("====================================================================\n\n")

for (f in c(A3B0, P1, P2)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

cov <- read.csv(A3B0, stringsAsFactors = FALSE, check.names = FALSE)
p1  <- read.delim(P1, stringsAsFactors = FALSE, check.names = FALSE)
p2  <- read.delim(P2, stringsAsFactors = FALSE, check.names = FALSE)

req_cov <- c("panel", "gene", "present_in_Flex_RNA")
req_p1 <- c("gene", "technical_class")
req_p2 <- c(
  "sample", "gene", "lineage", "full_lineage_n_cells",
  "strict_shared_UMI",
  "strict_shared_UMI_per_1000_full_lineage_cells",
  "standard_matrix_UMI",
  "standard_matrix_UMI_per_1000_cells",
  "n_strict_shared_probe_ids",
  "strict_shared_probe_ids"
)

if (!all(req_cov %in% names(cov))) stop("A.3B.0 schema mismatch.")
if (!all(req_p1 %in% names(p1))) stop("P.1 schema mismatch.")
if (!all(req_p2 %in% names(p2))) stop("P.2 schema mismatch.")

as_bool <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

cov$gene <- toupper(trimws(cov$gene))
p1$gene <- toupper(trimws(p1$gene))
p2$gene <- toupper(trimws(p2$gene))

strict_genes <- sort(unique(
  p1$gene[p1$technical_class == "STRICT_SHARED_PROBE_COMPARABLE"]
))

if (length(strict_genes) != 49L) {
  stop("Expected 49 strict-comparable genes; found ", length(strict_genes))
}

measured_cov <- cov[as_bool(cov$present_in_Flex_RNA), , drop = FALSE]
strict_membership <- measured_cov[
  measured_cov$gene %in% strict_genes,
  c("panel", "gene"),
  drop = FALSE
]
strict_membership <- unique(strict_membership)
strict_membership <- strict_membership[
  order(strict_membership$panel, strict_membership$gene),
]

write.csv(
  strict_membership,
  file.path(TAB, "00_strict_program_membership_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

sample_order <- c("SfLB1", "SfLB2", "Tx17", "Tx5")
lineages <- c("Hepatocyte_like", "Mesenchymal_like")

classify_relation <- function(sf1, sf2, tx17, tx5) {
  z <- c(sf1, sf2, tx17, tx5)
  if (any(!is.finite(z))) return("Not_comparable_missing_value")

  s <- c(sf1, sf2)
  t <- c(tx17, tx5)

  if (min(t) > max(s)) {
    return("Tx_higher_complete_replicate_range_separation")
  }
  if (max(t) < min(s)) {
    return("Tx_lower_complete_replicate_range_separation")
  }
  "Replicate_ranges_overlap"
}

mean_direction <- function(sf_mean, tx_mean, tol = 1e-12) {
  d <- tx_mean - sf_mean
  if (!is.finite(d)) return("NA")
  if (d > tol) return("Tx_higher_mean")
  if (d < -tol) return("Tx_lower_mean")
  "No_mean_difference"
}

# Ensure numeric conversion is explicit.
num_cols <- c(
  "full_lineage_n_cells",
  "strict_shared_UMI",
  "strict_shared_UMI_per_1000_full_lineage_cells",
  "standard_matrix_UMI",
  "standard_matrix_UMI_per_1000_cells",
  "n_strict_shared_probe_ids"
)

for (cc in num_cols) {
  p2[[cc]] <- suppressWarnings(as.numeric(p2[[cc]]))
}

p2c <- p2[
  p2$gene %in% strict_genes &
    p2$lineage %in% lineages &
    p2$sample %in% sample_order,
  ,
  drop = FALSE
]

# Validate one row/gene/lineage/sample.
key <- paste(p2c$gene, p2c$lineage, p2c$sample, sep = "||")
if (anyDuplicated(key)) {
  stop("Duplicated P.2 gene-lineage-sample row detected.")
}

driver_rows <- list()
k <- 1L

for (g in strict_genes) {
  for (ll in lineages) {
    z <- p2c[p2c$gene == g & p2c$lineage == ll, , drop = FALSE]

    if (!all(sample_order %in% z$sample)) {
      stop("Incomplete 4-sample rows for ", g, " / ", ll)
    }

    z <- z[match(sample_order, z$sample), , drop = FALSE]

    strict_vals <- setNames(
      z$strict_shared_UMI_per_1000_full_lineage_cells,
      z$sample
    )

    std_vals <- setNames(
      z$standard_matrix_UMI_per_1000_cells,
      z$sample
    )

    strict_sf <- mean(strict_vals[c("SfLB1", "SfLB2")])
    strict_tx <- mean(strict_vals[c("Tx17", "Tx5")])

    std_sf <- mean(std_vals[c("SfLB1", "SfLB2")])
    std_tx <- mean(std_vals[c("Tx17", "Tx5")])

    strict_rel <- classify_relation(
      strict_vals["SfLB1"],
      strict_vals["SfLB2"],
      strict_vals["Tx17"],
      strict_vals["Tx5"]
    )

    std_rel <- classify_relation(
      std_vals["SfLB1"],
      std_vals["SfLB2"],
      std_vals["Tx17"],
      std_vals["Tx5"]
    )

    strict_dir <- mean_direction(strict_sf, strict_tx)
    std_dir <- mean_direction(std_sf, std_tx)

    driver_rows[[k]] <- data.frame(
      gene = g,
      lineage = ll,

      strict_SfLB1 = strict_vals["SfLB1"],
      strict_SfLB2 = strict_vals["SfLB2"],
      strict_Tx17 = strict_vals["Tx17"],
      strict_Tx5 = strict_vals["Tx5"],
      strict_SfLB_mean = strict_sf,
      strict_Tx_mean = strict_tx,
      strict_Tx_minus_SfLB = strict_tx - strict_sf,
      strict_replicate_relation = strict_rel,
      strict_mean_direction = strict_dir,

      standard_SfLB1 = std_vals["SfLB1"],
      standard_SfLB2 = std_vals["SfLB2"],
      standard_Tx17 = std_vals["Tx17"],
      standard_Tx5 = std_vals["Tx5"],
      standard_SfLB_mean = std_sf,
      standard_Tx_mean = std_tx,
      standard_Tx_minus_SfLB = std_tx - std_sf,
      standard_replicate_relation = std_rel,
      standard_mean_direction = std_dir,

      mean_direction_concordant = strict_dir == std_dir,
      replicate_relation_exact_match = strict_rel == std_rel,

      n_strict_shared_probe_ids = unique(z$n_strict_shared_probe_ids)[1],
      strict_shared_probe_ids = unique(z$strict_shared_probe_ids)[1],
      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }
}

drivers <- do.call(rbind, driver_rows)

write.csv(
  drivers,
  file.path(TAB, "01_gene_level_standard_vs_strict_transition_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

# Expand to program membership: a gene may belong to >1 program.
detail <- merge(
  strict_membership,
  drivers,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

detail <- detail[order(detail$panel, detail$lineage, detail$gene), ]

write.csv(
  detail,
  file.path(TAB, "02_program_gene_driver_detail_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

collapse_genes <- function(x) {
  x <- sort(unique(x))
  if (length(x) == 0) "" else paste(x, collapse = ";")
}

summary_rows <- list()
k <- 1L

for (pp in sort(unique(detail$panel))) {
  for (ll in lineages) {
    z <- detail[
      detail$panel == pp & detail$lineage == ll,
      ,
      drop = FALSE
    ]

    if (nrow(z) == 0) next

    high <- z$gene[
      z$strict_replicate_relation ==
        "Tx_higher_complete_replicate_range_separation"
    ]
    low <- z$gene[
      z$strict_replicate_relation ==
        "Tx_lower_complete_replicate_range_separation"
    ]
    overlap <- z$gene[
      z$strict_replicate_relation == "Replicate_ranges_overlap"
    ]

    dir_conc <- z$gene[z$mean_direction_concordant]
    dir_disc <- z$gene[!z$mean_direction_concordant]

    rel_match <- z$gene[z$replicate_relation_exact_match]
    rel_mismatch <- z$gene[!z$replicate_relation_exact_match]

    summary_rows[[k]] <- data.frame(
      panel = pp,
      lineage = ll,
      n_strict_genes = nrow(z),

      n_Tx_higher_complete = length(high),
      Tx_higher_complete_genes = collapse_genes(high),

      n_Tx_lower_complete = length(low),
      Tx_lower_complete_genes = collapse_genes(low),

      n_replicate_ranges_overlap = length(overlap),
      replicate_ranges_overlap_genes = collapse_genes(overlap),

      n_standard_vs_strict_mean_direction_concordant = length(dir_conc),
      fraction_standard_vs_strict_mean_direction_concordant =
        length(dir_conc) / nrow(z),
      direction_concordant_genes = collapse_genes(dir_conc),

      n_standard_vs_strict_mean_direction_discordant = length(dir_disc),
      direction_discordant_genes = collapse_genes(dir_disc),

      n_exact_replicate_relation_match = length(rel_match),
      fraction_exact_replicate_relation_match =
        length(rel_match) / nrow(z),
      exact_relation_match_genes = collapse_genes(rel_match),
      relation_mismatch_genes = collapse_genes(rel_mismatch),

      stringsAsFactors = FALSE
    )

    k <- k + 1L
  }
}

summary_tbl <- do.call(rbind, summary_rows)

write.csv(
  summary_tbl,
  file.path(TAB, "03_program_direction_and_concordance_summary_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

# Dedicated antigen-presentation detail.
antigen_genes <- c("NLRC5", "PSMB10", "TAP2", "TAPBP")

antigen <- drivers[drivers$gene %in% antigen_genes, , drop = FALSE]

# Add TAP2 exclusion row explicitly so the table documents the boundary.
if (!"TAP2" %in% antigen$gene) {
  tap2_rows <- data.frame(
    gene = "TAP2",
    lineage = lineages,
    strict_SfLB1 = NA_real_,
    strict_SfLB2 = NA_real_,
    strict_Tx17 = NA_real_,
    strict_Tx5 = NA_real_,
    strict_SfLB_mean = NA_real_,
    strict_Tx_mean = NA_real_,
    strict_Tx_minus_SfLB = NA_real_,
    strict_replicate_relation =
      "EXCLUDED_NO_SHARED_INCLUDED_FILTERED_TRUE_PROBE_ALL4",
    strict_mean_direction = "NA",
    standard_SfLB1 = NA_real_,
    standard_SfLB2 = NA_real_,
    standard_Tx17 = NA_real_,
    standard_Tx5 = NA_real_,
    standard_SfLB_mean = NA_real_,
    standard_Tx_mean = NA_real_,
    standard_Tx_minus_SfLB = NA_real_,
    standard_replicate_relation = "NOT_USED_FOR_FINAL_PROBE_AWARE_CALL",
    standard_mean_direction = "NA",
    mean_direction_concordant = NA,
    replicate_relation_exact_match = NA,
    n_strict_shared_probe_ids = 0,
    strict_shared_probe_ids = "",
    stringsAsFactors = FALSE
  )
  antigen <- rbind(antigen, tap2_rows)
}

antigen <- antigen[order(antigen$gene, antigen$lineage), ]

write.csv(
  antigen,
  file.path(TAB, "04_antigen4_probe_aware_gene_driver_detail_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

# Global concordance summary across 49 x 2 = 98 gene-lineage comparisons.
global <- data.frame(
  metric = c(
    "strict_gene_lineage_comparisons",
    "mean_direction_concordant_n",
    "mean_direction_concordant_fraction",
    "exact_replicate_relation_match_n",
    "exact_replicate_relation_match_fraction",
    "strict_Tx_higher_complete_n",
    "strict_Tx_lower_complete_n",
    "strict_overlap_n"
  ),
  value = c(
    nrow(drivers),
    sum(drivers$mean_direction_concordant),
    mean(drivers$mean_direction_concordant),
    sum(drivers$replicate_relation_exact_match),
    mean(drivers$replicate_relation_exact_match),
    sum(
      drivers$strict_replicate_relation ==
        "Tx_higher_complete_replicate_range_separation"
    ),
    sum(
      drivers$strict_replicate_relation ==
        "Tx_lower_complete_replicate_range_separation"
    ),
    sum(
      drivers$strict_replicate_relation ==
        "Replicate_ranges_overlap"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  global,
  file.path(TAB, "05_global_standard_vs_strict_concordance_v1.4.4A.3B.1P.3.csv"),
  row.names = FALSE
)

cat("=== GLOBAL STANDARD vs STRICT CONCORDANCE ===\n")
print(global, row.names = FALSE)

cat("\n=== PROGRAM / LINEAGE GENE-DIRECTION SUMMARY ===\n")
print(
  summary_tbl[, c(
    "panel",
    "lineage",
    "n_strict_genes",
    "n_Tx_higher_complete",
    "n_Tx_lower_complete",
    "n_replicate_ranges_overlap",
    "fraction_standard_vs_strict_mean_direction_concordant"
  )],
  row.names = FALSE
)

cat("\n=== PROGRAM DRIVERS ===\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  cat("\n", r$panel, " | ", r$lineage, "\n", sep = "")
  cat("  Tx higher complete: ",
      ifelse(r$Tx_higher_complete_genes == "", "(none)",
             r$Tx_higher_complete_genes), "\n", sep = "")
  cat("  Tx lower complete : ",
      ifelse(r$Tx_lower_complete_genes == "", "(none)",
             r$Tx_lower_complete_genes), "\n", sep = "")
  cat("  Overlap           : ",
      ifelse(r$replicate_ranges_overlap_genes == "", "(none)",
             r$replicate_ranges_overlap_genes), "\n", sep = "")
  if (r$direction_discordant_genes != "") {
    cat("  Standard/strict mean-direction discordant: ",
        r$direction_discordant_genes, "\n", sep = "")
  }
}

cat("\n=== ANTIGEN-PRESENTATION BOUNDARY ===\n")
for (g in antigen_genes) {
  z <- antigen[antigen$gene == g, , drop = FALSE]
  cat("\n", g, "\n", sep = "")
  for (j in seq_len(nrow(z))) {
    cat(
      "  ", z$lineage[j], ": ",
      z$strict_replicate_relation[j],
      "\n",
      sep = ""
    )
  }
}

note <- c(
  "Human MASH v1.4.4A.3B.1P.3 interpretation boundary",
  "",
  "- This is a gene-driver and technical-concordance audit.",
  "- No new aggregate program score is created.",
  "- No inferential test is performed (n=2/group).",
  "- Complete replicate-range separation is descriptive evidence only.",
  "- Standard-vs-strict concordance evaluates whether the original gene-level",
  "  direction survives restriction to four-sample shared included+filtered TRUE probes.",
  "- TAP2 is excluded from final probe-aware interpretation because no strict",
  "  shared included+filtered TRUE probe exists across all four samples.",
  "- Antigen-presentation behavior must remain gene-specific because only 3/4",
  "  measured genes are strict-comparable and the larger canonical HLA/B2M panel",
  "  is outside the Flex feature space.",
  "- Even technically robust SfLB-vs-Tx differences remain engraftment-associated,",
  "  not IL10-specific, because transgene status, transplantation, environment,",
  "  selection/state transition, and possible residual xenogeneic pressure are confounded."
)

writeLines(
  note,
  file.path(OUT, "README_interpretation_v1.4.4A.3B.1P.3.txt"),
  useBytes = TRUE
)

cat("\n====================================================================\n")
cat("v1.4.4A.3B.1P.3 COMPLETE\n")
cat("====================================================================\n")
cat("Output:", OUT, "\n")
cat("Strict genes:", length(strict_genes), "\n")
cat("Comparable gene-lineage comparisons:", nrow(drivers), "\n")
cat("No HDF5 scan performed.\n")
cat("No aggregate program score rebuilt yet.\n")
