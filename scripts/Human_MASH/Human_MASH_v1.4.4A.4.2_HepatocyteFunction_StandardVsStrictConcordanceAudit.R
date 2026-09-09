#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.4.2_HepatocyteFunction_StandardVsStrictConcordanceAudit.R
#
# Purpose
# -------
# SSD-only concordance audit after v1.4.4A.4.1.
#
# Compare, within Hepatocyte_like cells:
#   A) standard RNA-matrix expression normalized per cell as log1p(CP10k)
#   B) strict shared-probe molecule reconstruction from v1.4.4A.4.1
#
# for the 16 Flex-measurable hepatic function/maturity genes.
#
# This step checks whether the gene-level SfLB-vs-Tx direction survives
# probe-version control and per-cell library-size normalization.
#
# No HDF5 scan. External RH241011956/RH251217117 drives are NOT required.
# No inferential test (n=2/group).
# No final aggregate hepatic-function score is created.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

VERSION <- "v1.4.4A.4.2"

BASE <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2A.1/objects",
  "Human_MASH_iPSC_derived_v1.4.2A.1_hepatocyte_state_interpreted.rds"
)

STRICT_FILE <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.4.1",
  "analysis/v1.4.4A.4.1_hepatocyte_function_strict_shared_probe_reconstruction",
  "tables/v1.4.4A.4.1_Hepatocyte_like_strict_gene_transition.tsv"
)

OUT <- file.path(
  BASE,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.4.2",
  "analysis/v1.4.4A.4.2_hepatocyte_function_standard_vs_strict_concordance"
)

TAB <- file.path(OUT, "tables")
LOG <- file.path(OUT, "logs")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

sink(file.path(LOG, paste0(VERSION, "_run_log.txt")), split = TRUE)
on.exit(sink(), add = TRUE)

cat("====================================================================\n")
cat("Human MASH ", VERSION, "\n", sep = "")
cat("Hepatocyte function: standard vs strict shared-probe concordance\n")
cat("====================================================================\n\n")
cat("SSD-only analysis; external delivery HDDs are not required.\n\n")

for (f in c(INPUT_RDS, STRICT_FILE)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

PANELS <- list(
  Secretory_synthetic = c("TTR", "APOA1", "FGB", "FGG"),
  Urea_nitrogen = c("OTC", "ASS1", "ASL", "NAGS"),
  Bile_transport = c("ABCB11"),
  Lipid_metabolic = c("MTTP", "FABP1"),
  Hepatocyte_identity = c("HNF4A", "CEBPA"),
  Epithelial_immaturity_transition = c("EPCAM", "KRT7", "SOX9")
)

GENES <- sort(unique(unlist(PANELS, use.names = FALSE)))
SAMPLES <- c("SfLB1", "SfLB2", "Tx17", "Tx5")
LINEAGE <- "Hepatocyte_like"

obj <- readRDS(INPUT_RDS)

if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.")
if (!"RNA" %in% Assays(obj)) stop("RNA assay missing.")
if (ncol(obj) != 43415L) stop("Unexpected cell count: ", ncol(obj))

required_meta <- c("sample_v141", "broad_lineage_manual_v142A")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

rna_features <- rownames(obj[["RNA"]])
missing_genes <- setdiff(GENES, rna_features)
if (length(missing_genes) > 0) {
  stop("Expected hepatic genes missing from RNA feature space: ",
       paste(missing_genes, collapse = ", "))
}

counts <- tryCatch(
  LayerData(obj, assay = "RNA", layer = "counts"),
  error = function(e1) {
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e2) {
        GetAssayData(obj, assay = "RNA", slot = "counts")
      }
    )
  }
)

if (ncol(counts) != ncol(obj)) stop("Counts/cell dimension mismatch.")

libsize <- Matrix::colSums(counts)
if (any(!is.finite(libsize)) || any(libsize <= 0)) {
  stop("Invalid RNA library size detected.")
}

md <- obj@meta.data
keep <- which(
  as.character(md$broad_lineage_manual_v142A) == LINEAGE &
    as.character(md$sample_v141) %in% SAMPLES
)

if (length(keep) == 0) stop("No Hepatocyte_like cells found.")

sel_counts <- counts[GENES, keep, drop = FALSE]
sel_dense <- as.matrix(sel_counts)
sel_lib <- libsize[keep]

log_cp10k <- log1p(
  sweep(sel_dense, 2, sel_lib, "/") * 10000
)

sample_vec <- as.character(md$sample_v141[keep])

cat("Hepatocyte_like cells by sample:\n")
print(table(sample_vec))
cat("\n")

standard_rows <- list()
k <- 1L

for (ss in SAMPLES) {
  ii <- which(sample_vec == ss)
  if (length(ii) == 0) stop("No Hepatocyte_like cells for ", ss)

  for (gg in GENES) {
    gi <- match(gg, rownames(sel_dense))
    x_count <- sel_dense[gi, ii]
    x_log <- log_cp10k[gi, ii]

    standard_rows[[k]] <- data.frame(
      sample = ss,
      gene = gg,
      lineage = LINEAGE,
      n_cells = length(ii),
      total_UMI = sum(x_count),
      UMI_per_1000_cells = sum(x_count) / length(ii) * 1000,
      positive_cells = sum(x_count > 0),
      positive_fraction = mean(x_count > 0),
      mean_logCP10k = mean(x_log),
      median_logCP10k = median(x_log),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

standard <- do.call(rbind, standard_rows)

write.csv(
  standard,
  file.path(TAB, paste0(VERSION, "_standard_matrix_Hepatocyte_like_gene_metrics.csv")),
  row.names = FALSE
)

strict <- read.delim(
  STRICT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

req_strict <- c(
  "gene",
  "SfLB1_strict_UMI_per_1000",
  "SfLB2_strict_UMI_per_1000",
  "Tx17_strict_UMI_per_1000",
  "Tx5_strict_UMI_per_1000",
  "replicate_relation"
)

if (!all(req_strict %in% names(strict))) {
  stop("Strict A.4.1 transition table schema mismatch.")
}

strict$gene <- toupper(trimws(strict$gene))

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

get_metric <- function(gene, metric) {
  z <- standard[standard$gene == gene, , drop = FALSE]
  if (!all(SAMPLES %in% z$sample)) {
    stop("Incomplete standard rows for ", gene)
  }
  z <- z[match(SAMPLES, z$sample), , drop = FALSE]
  setNames(as.numeric(z[[metric]]), z$sample)
}

compare_rows <- list()
k <- 1L

for (gg in GENES) {
  srow <- strict[strict$gene == gg, , drop = FALSE]
  if (nrow(srow) != 1L) stop("Expected one strict row for ", gg)

  std_log <- get_metric(gg, "mean_logCP10k")
  std_pos <- get_metric(gg, "positive_fraction")
  std_umi <- get_metric(gg, "UMI_per_1000_cells")

  strict_vals <- c(
    SfLB1 = as.numeric(srow$SfLB1_strict_UMI_per_1000),
    SfLB2 = as.numeric(srow$SfLB2_strict_UMI_per_1000),
    Tx17  = as.numeric(srow$Tx17_strict_UMI_per_1000),
    Tx5   = as.numeric(srow$Tx5_strict_UMI_per_1000)
  )

  strict_sf_mean <- mean(strict_vals[c("SfLB1", "SfLB2")])
  strict_tx_mean <- mean(strict_vals[c("Tx17", "Tx5")])

  std_log_sf_mean <- mean(std_log[c("SfLB1", "SfLB2")])
  std_log_tx_mean <- mean(std_log[c("Tx17", "Tx5")])

  std_pos_sf_mean <- mean(std_pos[c("SfLB1", "SfLB2")])
  std_pos_tx_mean <- mean(std_pos[c("Tx17", "Tx5")])

  strict_dir <- mean_direction(strict_sf_mean, strict_tx_mean)
  std_log_dir <- mean_direction(std_log_sf_mean, std_log_tx_mean)

  strict_rel <- as.character(srow$replicate_relation)
  std_log_rel <- classify_relation(
    std_log["SfLB1"], std_log["SfLB2"],
    std_log["Tx17"], std_log["Tx5"]
  )

  compare_rows[[k]] <- data.frame(
    gene = gg,
    panel = paste(
      names(PANELS)[vapply(PANELS, function(x) gg %in% x, logical(1))],
      collapse = ";"
    ),

    strict_relation = strict_rel,
    strict_mean_direction = strict_dir,
    strict_SfLB_mean_UMI_per_1000 = strict_sf_mean,
    strict_Tx_mean_UMI_per_1000 = strict_tx_mean,

    standard_logCP10k_relation = std_log_rel,
    standard_logCP10k_mean_direction = std_log_dir,
    standard_SfLB_mean_logCP10k = std_log_sf_mean,
    standard_Tx_mean_logCP10k = std_log_tx_mean,
    standard_Tx_minus_SfLB_mean_logCP10k =
      std_log_tx_mean - std_log_sf_mean,

    standard_SfLB_mean_positive_fraction = std_pos_sf_mean,
    standard_Tx_mean_positive_fraction = std_pos_tx_mean,
    standard_Tx_minus_SfLB_positive_fraction =
      std_pos_tx_mean - std_pos_sf_mean,

    standard_SfLB1_mean_logCP10k = std_log["SfLB1"],
    standard_SfLB2_mean_logCP10k = std_log["SfLB2"],
    standard_Tx17_mean_logCP10k = std_log["Tx17"],
    standard_Tx5_mean_logCP10k = std_log["Tx5"],

    strict_vs_standard_mean_direction_concordant =
      strict_dir == std_log_dir,

    strict_vs_standard_exact_relation_match =
      strict_rel == std_log_rel,

    stringsAsFactors = FALSE
  )

  k <- k + 1L
}

cmp <- do.call(rbind, compare_rows)

write.csv(
  cmp,
  file.path(TAB, paste0(VERSION, "_strict_vs_standard_Hepatocyte_like_concordance.csv")),
  row.names = FALSE
)

panel_rows <- list()
k <- 1L

for (pp in names(PANELS)) {
  z <- cmp[cmp$gene %in% PANELS[[pp]], , drop = FALSE]

  panel_rows[[k]] <- data.frame(
    panel = pp,
    n_genes = nrow(z),

    n_strict_Tx_higher_complete = sum(
      z$strict_relation ==
        "Tx_higher_complete_replicate_range_separation"
    ),
    strict_Tx_higher_complete_genes = paste(
      z$gene[
        z$strict_relation ==
          "Tx_higher_complete_replicate_range_separation"
      ],
      collapse = ";"
    ),

    n_strict_Tx_lower_complete = sum(
      z$strict_relation ==
        "Tx_lower_complete_replicate_range_separation"
    ),
    strict_Tx_lower_complete_genes = paste(
      z$gene[
        z$strict_relation ==
          "Tx_lower_complete_replicate_range_separation"
      ],
      collapse = ";"
    ),

    n_strict_overlap = sum(
      z$strict_relation == "Replicate_ranges_overlap"
    ),

    n_direction_concordant = sum(
      z$strict_vs_standard_mean_direction_concordant
    ),
    fraction_direction_concordant = mean(
      z$strict_vs_standard_mean_direction_concordant
    ),

    n_exact_relation_match = sum(
      z$strict_vs_standard_exact_relation_match
    ),
    fraction_exact_relation_match = mean(
      z$strict_vs_standard_exact_relation_match
    ),

    discordant_genes = paste(
      z$gene[
        !z$strict_vs_standard_mean_direction_concordant
      ],
      collapse = ";"
    ),

    stringsAsFactors = FALSE
  )

  k <- k + 1L
}

panel_summary <- do.call(rbind, panel_rows)

write.csv(
  panel_summary,
  file.path(TAB, paste0(VERSION, "_panel_concordance_summary.csv")),
  row.names = FALSE
)

global <- data.frame(
  metric = c(
    "genes",
    "strict_vs_standard_mean_direction_concordant_n",
    "strict_vs_standard_mean_direction_concordant_fraction",
    "strict_vs_standard_exact_relation_match_n",
    "strict_vs_standard_exact_relation_match_fraction"
  ),
  value = c(
    nrow(cmp),
    sum(cmp$strict_vs_standard_mean_direction_concordant),
    mean(cmp$strict_vs_standard_mean_direction_concordant),
    sum(cmp$strict_vs_standard_exact_relation_match),
    mean(cmp$strict_vs_standard_exact_relation_match)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  global,
  file.path(TAB, paste0(VERSION, "_global_concordance_summary.csv")),
  row.names = FALSE
)

cat("=== GLOBAL CONCORDANCE ===\n")
print(global, row.names = FALSE)

cat("\n=== GENE-LEVEL RESULT ===\n")
print(
  cmp[, c(
    "gene",
    "panel",
    "strict_relation",
    "standard_logCP10k_relation",
    "strict_vs_standard_mean_direction_concordant",
    "standard_Tx_minus_SfLB_mean_logCP10k",
    "standard_Tx_minus_SfLB_positive_fraction"
  )],
  row.names = FALSE
)

cat("\n=== PANEL SUMMARY ===\n")
print(
  panel_summary[, c(
    "panel",
    "n_genes",
    "n_strict_Tx_higher_complete",
    "n_strict_Tx_lower_complete",
    "n_strict_overlap",
    "fraction_direction_concordant",
    "fraction_exact_relation_match",
    "discordant_genes"
  )],
  row.names = FALSE
)

note <- c(
  "Human MASH v1.4.4A.4.2 interpretation boundary",
  "",
  "- This is an SSD-only standard-vs-strict technical concordance audit.",
  "- Standard RNA expression is normalized per cell as log1p(CP10k).",
  "- Strict reconstruction uses only four-sample shared included+filtered TRUE probes.",
  "- A gene-level SfLB-vs-Tx direction is more robust when both approaches agree.",
  "- No inferential test is performed (n=2/group).",
  "- Transcript levels do not establish hepatocyte biochemical function.",
  "- Because ALB, CPS1, ARG1, CYP3A4 and several other canonical hepatic genes",
  "  are outside the 9,137-feature Flex space, conclusions must be framed as",
  "  probe-measurable hepatic transcriptomic programs, not comprehensive function.",
  "- SfLB->Tx remains an engraftment-associated comparison and is not IL10-specific."
)

writeLines(
  note,
  file.path(OUT, paste0("README_interpretation_", VERSION, ".txt")),
  useBytes = TRUE
)

cat("\n====================================================================\n")
cat(VERSION, " COMPLETE\n")
cat("====================================================================\n")
cat("Output:", OUT, "\n")
cat("No external delivery HDD required.\n")
cat("No aggregate hepatic-function score created.\n")
