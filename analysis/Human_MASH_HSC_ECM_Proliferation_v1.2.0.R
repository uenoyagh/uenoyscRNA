#!/usr/bin/env Rscript

# ============================================================
# Human MASH HSC ECM two-axis validation + strict proliferation
# v1.2.0
#
# Canonical receiver:
#   Human_MASH_HSC_Receiver_Primary_v1.0.rds
#
# PURPOSE
#   1. Reproduce the v1.1 total ECM score
#   2. Split ECM activation into:
#        A. Collagen-production axis
#        B. ECM-remodeling / fibrogenic axis
#   3. Re-evaluate proliferation without CellCycleScoring()
#   4. Define strict Cycling-high HSC using direct proliferation genes
#   5. Quantify all endpoints by sample, HSC state, and NAS
#
# IMPORTANT
#   - No reclustering
#   - Source RDS is never overwritten
#   - Canonical HSC receiver remains fixed at 7,640 cells
#   - RNA LogNormalize data are used
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_Receiver_Primary_v1.0",
  "objects",
  "Human_MASH_HSC_Receiver_Primary_v1.0.rds"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_OBJ   <- file.path(OUT, "objects")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV,  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_HSC_ECM_Proliferation_v1.2.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH HSC ECM / Proliferation v1.2.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(INPUT_RDS))

# ------------------------------------------------------------
# 1. Load canonical HSC receiver
# ------------------------------------------------------------

obj <- readRDS(INPUT_RDS)

cat("Input RDS:\n", INPUT_RDS, "\n\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n\n")

if (ncol(obj) != 7640) {
  warning(
    "Canonical receiver expected 7640 cells, observed ",
    ncol(obj)
  )
}

DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------
# 2. Metadata discovery
# ------------------------------------------------------------

md <- obj[[]]

find_first_column <- function(candidates, metadata) {
  hit <- candidates[candidates %in% colnames(metadata)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

sample_col <- find_first_column(
  c("sample", "Sample", "sample_id", "sampleID", "orig.ident"),
  md
)

state_col <- find_first_column(
  c(
    "state",
    "state2",
    "HSC_state2",
    "receiver_state2",
    "primary_state2",
    "HSC_state",
    "HSC_state_v1.1"
  ),
  md
)

nas_col <- find_first_column(
  c("NAS", "nas", "NAS_score", "NASscore"),
  md
)

cat("Detected metadata columns\n")
cat("-------------------------\n")
cat("sample:", sample_col, "\n")
cat("state :", state_col, "\n")
cat("NAS   :", nas_col, "\n\n")

if (is.na(sample_col)) {
  stop("Sample metadata column could not be identified.")
}

# ------------------------------------------------------------
# 3. Recover canonical qHSC / aHSC state if necessary
# ------------------------------------------------------------

if (!is.na(state_col)) {

  obj$HSC_state_v1.2 <- as.character(
    obj[[]][[state_col]]
  )

} else {

  label_col <- find_first_column(
    c(
      "celltype",
      "cell_type",
      "annotation",
      "CellType",
      "celltype_final",
      "final_annotation"
    ),
    obj[[]]
  )

  if (is.na(label_col)) {
    stop(
      "Neither HSC state nor source annotation column could be identified."
    )
  }

  lab <- as.character(obj[[]][[label_col]])

  state <- rep(NA_character_, length(lab))

  state[
    grepl(
      "HSC_quiescent",
      lab,
      ignore.case = TRUE
    )
  ] <- "qHSC"

  state[
    grepl(
      "HSC_fibroblast_mesenchymal",
      lab,
      ignore.case = TRUE
    )
  ] <- "aHSC_mesenchymal"

  obj$HSC_state_v1.2 <- state

  cat(
    "State reconstructed from annotation column:",
    label_col,
    "\n"
  )
}

obj$HSC_state_v1.2[
  grepl(
    "^qHSC$|quiescent",
    obj$HSC_state_v1.2,
    ignore.case = TRUE
  )
] <- "qHSC"

obj$HSC_state_v1.2[
  grepl(
    "aHSC|fibroblast|mesenchymal",
    obj$HSC_state_v1.2,
    ignore.case = TRUE
  )
] <- "aHSC_mesenchymal"

if (any(is.na(obj$HSC_state_v1.2))) {
  stop("Some HSC cells could not be assigned to qHSC/aHSC.")
}

cat("\nHSC state counts\n")
print(table(obj$HSC_state_v1.2))
cat("\n")

# ------------------------------------------------------------
# 4. Normalize RNA
# ------------------------------------------------------------

obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 5. Gene panels
# ------------------------------------------------------------

ECM_TOTAL <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2",
  "FN1",
  "SPARC",
  "THBS2",
  "LOXL2",
  "TIMP1",
  "TGFBI",
  "CTHRC1",
  "SERPINE1",
  "POSTN"
)

ECM_COLLAGEN <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2"
)

ECM_REMODELING <- c(
  "FN1",
  "SPARC",
  "THBS2",
  "LOXL2",
  "TIMP1",
  "TGFBI",
  "CTHRC1",
  "SERPINE1",
  "POSTN"
)

PROLIFERATION_PANEL <- c(
  "MKI67",
  "TOP2A",
  "UBE2C",
  "CCNB1",
  "CCNB2",
  "CDC20",
  "BIRC5",
  "CENPF",
  "TYMS",
  "PCNA",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "MCM7"
)

PROLIFERATION_CORE <- c(
  "MKI67",
  "TOP2A",
  "UBE2C",
  "CCNB1",
  "CCNB2",
  "CDC20",
  "BIRC5",
  "CENPF"
)

# ------------------------------------------------------------
# 6. Helper functions
# ------------------------------------------------------------

available_panel <- function(panel, object) {
  intersect(panel, rownames(object))
}

missing_panel <- function(panel, object) {
  setdiff(panel, rownames(object))
}

z_fun <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sx
}

calc_equal_weight_zscore <- function(object, genes) {

  genes_use <- intersect(genes, rownames(object))

  if (length(genes_use) == 0) {
    stop("No genes available for requested score.")
  }

  x <- FetchData(
    object,
    vars = genes_use,
    layer = "data"
  )

  xz <- as.data.frame(
    lapply(x, z_fun)
  )

  rowMeans(
    xz,
    na.rm = TRUE
  )
}

count_detected <- function(object, genes) {

  genes_use <- intersect(genes, rownames(object))

  if (length(genes_use) == 0) {
    return(rep(0L, ncol(object)))
  }

  x <- FetchData(
    object,
    vars = genes_use,
    layer = "data"
  )

  rowSums(
    x > 0,
    na.rm = TRUE
  )
}

# ------------------------------------------------------------
# 7. Panel availability audit
# ------------------------------------------------------------

panel_audit <- bind_rows(
  tibble(
    panel = "ECM_TOTAL",
    requested = length(ECM_TOTAL),
    available = length(available_panel(ECM_TOTAL, obj)),
    available_genes = paste(available_panel(ECM_TOTAL, obj), collapse = ";"),
    missing_genes = paste(missing_panel(ECM_TOTAL, obj), collapse = ";")
  ),
  tibble(
    panel = "ECM_COLLAGEN",
    requested = length(ECM_COLLAGEN),
    available = length(available_panel(ECM_COLLAGEN, obj)),
    available_genes = paste(available_panel(ECM_COLLAGEN, obj), collapse = ";"),
    missing_genes = paste(missing_panel(ECM_COLLAGEN, obj), collapse = ";")
  ),
  tibble(
    panel = "ECM_REMODELING",
    requested = length(ECM_REMODELING),
    available = length(available_panel(ECM_REMODELING, obj)),
    available_genes = paste(available_panel(ECM_REMODELING, obj), collapse = ";"),
    missing_genes = paste(missing_panel(ECM_REMODELING, obj), collapse = ";")
  ),
  tibble(
    panel = "PROLIFERATION_PANEL",
    requested = length(PROLIFERATION_PANEL),
    available = length(available_panel(PROLIFERATION_PANEL, obj)),
    available_genes = paste(available_panel(PROLIFERATION_PANEL, obj), collapse = ";"),
    missing_genes = paste(missing_panel(PROLIFERATION_PANEL, obj), collapse = ";")
  ),
  tibble(
    panel = "PROLIFERATION_CORE",
    requested = length(PROLIFERATION_CORE),
    available = length(available_panel(PROLIFERATION_CORE, obj)),
    available_genes = paste(available_panel(PROLIFERATION_CORE, obj), collapse = ";"),
    missing_genes = paste(missing_panel(PROLIFERATION_CORE, obj), collapse = ";")
  )
)

write_tsv(
  panel_audit,
  file.path(
    DIR_TABLE,
    "01_gene_panel_availability.tsv"
  )
)

cat("Gene panel availability\n")
print(panel_audit)
cat("\n")

if (length(available_panel(ECM_TOTAL, obj)) < 5) {
  stop("Too few ECM genes are available.")
}

if (length(available_panel(PROLIFERATION_PANEL, obj)) < 5) {
  stop("Too few proliferation genes are available.")
}

# ------------------------------------------------------------
# 8. ECM scores
# ------------------------------------------------------------

obj$ECM_total_score_v1.2 <-
  calc_equal_weight_zscore(
    obj,
    ECM_TOTAL
  )

obj$ECM_collagen_score_v1.2 <-
  calc_equal_weight_zscore(
    obj,
    ECM_COLLAGEN
  )

obj$ECM_remodeling_score_v1.2 <-
  calc_equal_weight_zscore(
    obj,
    ECM_REMODELING
  )

# Reproduce the v1.1 primary ECM-high definition:
# pooled qHSC 90th percentile of total ECM score

is_qHSC <- obj$HSC_state_v1.2 == "qHSC"

ECM_primary_threshold <- as.numeric(
  quantile(
    obj$ECM_total_score_v1.2[is_qHSC],
    probs = 0.90,
    na.rm = TRUE
  )
)

obj$ECM_high_primary_v1.2 <-
  obj$ECM_total_score_v1.2 >= ECM_primary_threshold

# ------------------------------------------------------------
# 9. Proliferation score and direct-marker counts
# ------------------------------------------------------------

obj$Proliferation_score_v1.2 <-
  calc_equal_weight_zscore(
    obj,
    PROLIFERATION_PANEL
  )

obj$Proliferation_core_detected_n_v1.2 <-
  count_detected(
    obj,
    PROLIFERATION_CORE
  )

# Individual direct-marker detection
direct_markers <- available_panel(
  PROLIFERATION_CORE,
  obj
)

direct_expr <- FetchData(
  obj,
  vars = direct_markers,
  layer = "data"
)

for (g in direct_markers) {
  obj[[paste0(g, "_detected_v1.2")]] <-
    as.logical(direct_expr[[g]] > 0)
}

# ------------------------------------------------------------
# 10. Define strict proliferation reference
#
# Preferred reference:
#   NAS 0 qHSC
#
# Fallback:
#   pooled qHSC
#
# Primary strict definition:
#   proliferation score >= reference 95th percentile
#   AND >= 2 direct proliferation-core genes detected
#
# Sensitivity definition:
#   proliferation score >= reference 90th percentile
#   AND >= 1 direct proliferation-core gene detected
# ------------------------------------------------------------

md <- obj[[]]

if (!is.na(nas_col)) {

  nas_vector <- suppressWarnings(
    as.numeric(md[[nas_col]])
  )

} else {

  nas_vector <- rep(NA_real_, nrow(md))
}

ref_primary <-
  obj$HSC_state_v1.2 == "qHSC" &
  !is.na(nas_vector) &
  nas_vector == 0

if (sum(ref_primary) >= 100) {

  proliferation_reference <- "NAS0_qHSC"
  ref_idx <- ref_primary

} else {

  proliferation_reference <- "pooled_qHSC"
  ref_idx <- obj$HSC_state_v1.2 == "qHSC"
}

prolif_threshold_p90 <- as.numeric(
  quantile(
    obj$Proliferation_score_v1.2[ref_idx],
    probs = 0.90,
    na.rm = TRUE
  )
)

prolif_threshold_p95 <- as.numeric(
  quantile(
    obj$Proliferation_score_v1.2[ref_idx],
    probs = 0.95,
    na.rm = TRUE
  )
)

obj$Cycling_strict_PRIMARY_v1.2 <-
  obj$Proliferation_score_v1.2 >= prolif_threshold_p95 &
  obj$Proliferation_core_detected_n_v1.2 >= 2

obj$Cycling_sensitive_v1.2 <-
  obj$Proliferation_score_v1.2 >= prolif_threshold_p90 &
  obj$Proliferation_core_detected_n_v1.2 >= 1

# Orthogonal marker-only classification
# This does NOT use a percentile threshold.
# It is intended as a robustness check.
marker_only <- rep(FALSE, ncol(obj))

if ("MKI67" %in% direct_markers) {
  marker_only <- marker_only | (direct_expr$MKI67 > 0)
}

if ("TOP2A" %in% direct_markers) {
  marker_only <- marker_only | (direct_expr$TOP2A > 0)
}

other_core <- intersect(
  c(
    "UBE2C",
    "CCNB1",
    "CCNB2",
    "CDC20",
    "BIRC5",
    "CENPF"
  ),
  direct_markers
)

if (length(other_core) > 0) {

  other_detected_n <- rowSums(
    direct_expr[, other_core, drop = FALSE] > 0,
    na.rm = TRUE
  )

  marker_only <- marker_only | (other_detected_n >= 2)
}

obj$Cycling_marker_only_v1.2 <- marker_only

# ------------------------------------------------------------
# 11. Threshold table
# ------------------------------------------------------------

threshold_table <- tibble(
  item = c(
    "ECM_primary_threshold",
    "ECM_reference",
    "proliferation_reference",
    "proliferation_reference_cells",
    "proliferation_p90",
    "proliferation_p95_PRIMARY",
    "PRIMARY_core_detected_min",
    "SENSITIVE_core_detected_min"
  ),
  value = c(
    as.character(ECM_primary_threshold),
    "pooled_qHSC_90th_percentile",
    proliferation_reference,
    as.character(sum(ref_idx)),
    as.character(prolif_threshold_p90),
    as.character(prolif_threshold_p95),
    "2",
    "1"
  )
)

write_tsv(
  threshold_table,
  file.path(
    DIR_TABLE,
    "02_thresholds_v1.2.tsv"
  )
)

cat("Thresholds\n")
print(threshold_table)
cat("\n")

# ------------------------------------------------------------
# 12. Standardized cell-level table
# ------------------------------------------------------------

md <- obj[[]]

analysis_df <- tibble(
  cell = rownames(md),

  sample = as.character(
    md[[sample_col]]
  ),

  state = as.character(
    md$HSC_state_v1.2
  ),

  NAS = nas_vector,

  ECM_total_score = as.numeric(
    md$ECM_total_score_v1.2
  ),

  ECM_collagen_score = as.numeric(
    md$ECM_collagen_score_v1.2
  ),

  ECM_remodeling_score = as.numeric(
    md$ECM_remodeling_score_v1.2
  ),

  ECM_high_primary = as.logical(
    md$ECM_high_primary_v1.2
  ),

  Proliferation_score = as.numeric(
    md$Proliferation_score_v1.2
  ),

  Proliferation_core_detected_n = as.integer(
    md$Proliferation_core_detected_n_v1.2
  ),

  Cycling_strict_PRIMARY = as.logical(
    md$Cycling_strict_PRIMARY_v1.2
  ),

  Cycling_sensitive = as.logical(
    md$Cycling_sensitive_v1.2
  ),

  Cycling_marker_only = as.logical(
    md$Cycling_marker_only_v1.2
  )
)

write_tsv(
  analysis_df,
  file.path(
    DIR_TABLE,
    "03_HSC_cell_level_ECM_proliferation.tsv.gz"
  )
)

# ------------------------------------------------------------
# 13. Sample-level summary
# ------------------------------------------------------------

sample_summary <- analysis_df %>%
  group_by(sample) %>%
  summarise(
    HSC_n = n(),

    NAS =
      if (
        all(is.na(NAS))
      ) {
        NA_real_
      } else {
        unique(
          NAS[!is.na(NAS)]
        )[1]
      },

    ECM_total_mean =
      mean(
        ECM_total_score,
        na.rm = TRUE
      ),

    ECM_collagen_mean =
      mean(
        ECM_collagen_score,
        na.rm = TRUE
      ),

    ECM_remodeling_mean =
      mean(
        ECM_remodeling_score,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high_primary,
        na.rm = TRUE
      ),

    Proliferation_score_mean =
      mean(
        Proliferation_score,
        na.rm = TRUE
      ),

    Cycling_strict_fraction =
      mean(
        Cycling_strict_PRIMARY,
        na.rm = TRUE
      ),

    Cycling_sensitive_fraction =
      mean(
        Cycling_sensitive,
        na.rm = TRUE
      ),

    Cycling_marker_only_fraction =
      mean(
        Cycling_marker_only,
        na.rm = TRUE
      ),

    ECM_high_n =
      sum(
        ECM_high_primary,
        na.rm = TRUE
      ),

    Cycling_strict_in_ECM_high_n =
      sum(
        ECM_high_primary &
          Cycling_strict_PRIMARY,
        na.rm = TRUE
      ),

    Cycling_strict_fraction_in_ECM_high =
      ifelse(
        sum(
          ECM_high_primary,
          na.rm = TRUE
        ) > 0,

        mean(
          Cycling_strict_PRIMARY[
            ECM_high_primary
          ],
          na.rm = TRUE
        ),

        NA_real_
      ),

    Cycling_sensitive_fraction_in_ECM_high =
      ifelse(
        sum(
          ECM_high_primary,
          na.rm = TRUE
        ) > 0,

        mean(
          Cycling_sensitive[
            ECM_high_primary
          ],
          na.rm = TRUE
        ),

        NA_real_
      ),

    .groups = "drop"
  ) %>%
  arrange(sample)

write_tsv(
  sample_summary,
  file.path(
    DIR_TABLE,
    "04_sample_ECM_proliferation_summary.tsv"
  )
)

# ------------------------------------------------------------
# 14. Sample x HSC state summary
# ------------------------------------------------------------

sample_state_summary <- analysis_df %>%
  group_by(
    sample,
    state
  ) %>%
  summarise(
    cells = n(),

    NAS =
      if (
        all(is.na(NAS))
      ) {
        NA_real_
      } else {
        unique(
          NAS[!is.na(NAS)]
        )[1]
      },

    ECM_total_mean =
      mean(
        ECM_total_score,
        na.rm = TRUE
      ),

    ECM_collagen_mean =
      mean(
        ECM_collagen_score,
        na.rm = TRUE
      ),

    ECM_remodeling_mean =
      mean(
        ECM_remodeling_score,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high_primary,
        na.rm = TRUE
      ),

    Proliferation_score_mean =
      mean(
        Proliferation_score,
        na.rm = TRUE
      ),

    Cycling_strict_fraction =
      mean(
        Cycling_strict_PRIMARY,
        na.rm = TRUE
      ),

    Cycling_sensitive_fraction =
      mean(
        Cycling_sensitive,
        na.rm = TRUE
      ),

    Cycling_marker_only_fraction =
      mean(
        Cycling_marker_only,
        na.rm = TRUE
      ),

    Cycling_strict_fraction_in_ECM_high =
      ifelse(
        sum(
          ECM_high_primary,
          na.rm = TRUE
        ) > 0,

        mean(
          Cycling_strict_PRIMARY[
            ECM_high_primary
          ],
          na.rm = TRUE
        ),

        NA_real_
      ),

    .groups = "drop"
  )

write_tsv(
  sample_state_summary,
  file.path(
    DIR_TABLE,
    "05_sample_state_ECM_proliferation_summary.tsv"
  )
)

# ------------------------------------------------------------
# 15. NAS-level descriptive summary
#
# IMPORTANT:
# Average sample-level quantities so each biological sample
# contributes equally.
# ------------------------------------------------------------

if (!all(is.na(sample_summary$NAS))) {

  NAS_summary <- sample_summary %>%
    filter(!is.na(NAS)) %>%
    group_by(NAS) %>%
    summarise(
      samples = n(),

      ECM_total_mean =
        mean(
          ECM_total_mean,
          na.rm = TRUE
        ),

      ECM_collagen_mean =
        mean(
          ECM_collagen_mean,
          na.rm = TRUE
        ),

      ECM_remodeling_mean =
        mean(
          ECM_remodeling_mean,
          na.rm = TRUE
        ),

      ECM_high_fraction =
        mean(
          ECM_high_fraction,
          na.rm = TRUE
        ),

      Proliferation_score_mean =
        mean(
          Proliferation_score_mean,
          na.rm = TRUE
        ),

      Cycling_strict_fraction =
        mean(
          Cycling_strict_fraction,
          na.rm = TRUE
        ),

      Cycling_sensitive_fraction =
        mean(
          Cycling_sensitive_fraction,
          na.rm = TRUE
        ),

      Cycling_marker_only_fraction =
        mean(
          Cycling_marker_only_fraction,
          na.rm = TRUE
        ),

      Cycling_strict_fraction_in_ECM_high =
        mean(
          Cycling_strict_fraction_in_ECM_high,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  write_tsv(
    NAS_summary,
    file.path(
      DIR_TABLE,
      "06_NAS_mean_ECM_proliferation_summary.tsv"
    )
  )
}

# ------------------------------------------------------------
# 16. Sample-level Spearman tests vs NAS
# ------------------------------------------------------------

if (!all(is.na(sample_summary$NAS))) {

  metrics <- c(
    "ECM_total_mean",
    "ECM_collagen_mean",
    "ECM_remodeling_mean",
    "ECM_high_fraction",
    "Proliferation_score_mean",
    "Cycling_strict_fraction",
    "Cycling_sensitive_fraction",
    "Cycling_marker_only_fraction",
    "Cycling_strict_fraction_in_ECM_high"
  )

  spearman_list <- lapply(
    metrics,
    function(metric) {

      x <- sample_summary$NAS
      y <- sample_summary[[metric]]

      keep <- complete.cases(x, y)

      if (
        sum(keep) < 4 ||
        length(unique(x[keep])) < 2 ||
        length(unique(y[keep])) < 2
      ) {

        return(
          tibble(
            metric = metric,
            n_samples = sum(keep),
            rho = NA_real_,
            p_value = NA_real_
          )
        )
      }

      ct <- suppressWarnings(
        cor.test(
          x[keep],
          y[keep],
          method = "spearman",
          exact = FALSE
        )
      )

      tibble(
        metric = metric,
        n_samples = sum(keep),
        rho = unname(ct$estimate),
        p_value = ct$p.value
      )
    }
  )

  spearman <- bind_rows(
    spearman_list
  )

  write_tsv(
    spearman,
    file.path(
      DIR_TABLE,
      "07_NAS_sample_level_Spearman.tsv"
    )
  )
}

# ------------------------------------------------------------
# 17. ECM two-axis concordance
# ------------------------------------------------------------

axis_tests <- list(
  total_vs_collagen = c(
    "ECM_total_mean",
    "ECM_collagen_mean"
  ),
  total_vs_remodeling = c(
    "ECM_total_mean",
    "ECM_remodeling_mean"
  ),
  collagen_vs_remodeling = c(
    "ECM_collagen_mean",
    "ECM_remodeling_mean"
  )
)

axis_concordance <- bind_rows(
  lapply(
    names(axis_tests),
    function(nm) {

      vars <- axis_tests[[nm]]

      x <- sample_summary[[vars[1]]]
      y <- sample_summary[[vars[2]]]

      keep <- complete.cases(x, y)

      ct <- suppressWarnings(
        cor.test(
          x[keep],
          y[keep],
          method = "spearman",
          exact = FALSE
        )
      )

      tibble(
        comparison = nm,
        variable_x = vars[1],
        variable_y = vars[2],
        n_samples = sum(keep),
        rho = unname(ct$estimate),
        p_value = ct$p.value
      )
    }
  )
)

write_tsv(
  axis_concordance,
  file.path(
    DIR_TABLE,
    "08_ECM_axis_concordance.tsv"
  )
)

# ------------------------------------------------------------
# 18. Direct proliferation-marker detection by sample
# ------------------------------------------------------------

marker_detection_df <- tibble(
  sample = analysis_df$sample
)

for (g in direct_markers) {

  marker_detection_df[[g]] <-
    as.logical(
      obj[[]][[paste0(g, "_detected_v1.2")]]
    )
}

marker_detection_by_sample <- marker_detection_df %>%
  group_by(sample) %>%
  summarise(
    across(
      all_of(direct_markers),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_fraction"
    ),
    .groups = "drop"
  )

write_tsv(
  marker_detection_by_sample,
  file.path(
    DIR_TABLE,
    "09_direct_proliferation_marker_detection_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 19. Figures
# ------------------------------------------------------------

sample_levels <- paste0("LT", 1:8)

sample_summary$sample <- factor(
  sample_summary$sample,
  levels = sample_levels
)

analysis_df$sample <- factor(
  analysis_df$sample,
  levels = sample_levels
)

# Figure 1: ECM two axes by sample
axis_long <- sample_summary %>%
  select(
    sample,
    ECM_collagen_mean,
    ECM_remodeling_mean
  ) %>%
  pivot_longer(
    cols = c(
      ECM_collagen_mean,
      ECM_remodeling_mean
    ),
    names_to = "axis",
    values_to = "score"
  )

p1 <- ggplot(
  axis_long,
  aes(
    x = sample,
    y = score,
    group = axis,
    linetype = axis
  )
) +
  geom_point(size = 2.5) +
  geom_line() +
  theme_classic(base_size = 14) +
  labs(
    title = "Human MASH HSC ECM two-axis validation",
    x = NULL,
    y = "Mean z-score",
    linetype = "ECM axis"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_ECM_two_axes_by_sample.pdf"
  ),
  p1,
  width = 8,
  height = 5
)

# Figure 2: ECM-high fraction
p2 <- ggplot(
  sample_summary,
  aes(
    x = sample,
    y = ECM_high_fraction
  )
) +
  geom_col() +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "ECM-high HSC fraction",
    subtitle = "Primary definition: pooled qHSC 90th percentile",
    x = NULL,
    y = "ECM-high fraction"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_ECM_high_fraction_by_sample.pdf"
  ),
  p2,
  width = 7.5,
  height = 5
)

# Figure 3: strict cycling fraction
p3 <- ggplot(
  sample_summary,
  aes(
    x = sample,
    y = Cycling_strict_fraction
  )
) +
  geom_col() +
  theme_classic(base_size = 14) +
  labs(
    title = "Strict Cycling-high HSC fraction",
    subtitle = paste0(
      proliferation_reference,
      " p95 + >=2 direct proliferation markers"
    ),
    x = NULL,
    y = "Strict Cycling-high fraction"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_Strict_Cycling_fraction_by_sample.pdf"
  ),
  p3,
  width = 7.5,
  height = 5
)

# Figure 4: strict cycling within ECM-high HSC
p4 <- ggplot(
  sample_summary,
  aes(
    x = sample,
    y = Cycling_strict_fraction_in_ECM_high
  )
) +
  geom_col() +
  theme_classic(base_size = 14) +
  labs(
    title = "Strict Cycling-high fraction within ECM-high HSC",
    x = NULL,
    y = "Strict Cycling-high fraction"
  )

ggsave(
  file.path(
    DIR_FIG,
    "04_Strict_Cycling_within_ECM_high_by_sample.pdf"
  ),
  p4,
  width = 7.5,
  height = 5
)

# Figure 5: ECM total vs proliferation score
p5 <- ggplot(
  sample_summary,
  aes(
    x = ECM_total_mean,
    y = Proliferation_score_mean,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "ECM activation vs proliferation",
    x = "Mean ECM total score",
    y = "Mean proliferation score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "05_ECM_vs_Proliferation_sample_level.pdf"
  ),
  p5,
  width = 6.5,
  height = 5.5
)

# ------------------------------------------------------------
# 20. Save annotated analysis object
# ------------------------------------------------------------

saveRDS(
  obj,
  file.path(
    DIR_OBJ,
    "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
  )
)

# ------------------------------------------------------------
# 21. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "source_RDS",
    "source_cells",
    "reclustering",
    "source_overwritten",
    "expression_assay",
    "normalization",
    "ECM_total_score_method",
    "ECM_primary_threshold",
    "ECM_collagen_axis",
    "ECM_remodeling_axis",
    "proliferation_score_method",
    "proliferation_reference",
    "strict_cycling_definition",
    "sensitive_cycling_definition",
    "marker_only_definition"
  ),

  value = c(
    "v1.2.0",
    INPUT_RDS,
    as.character(ncol(obj)),
    "none",
    "FALSE",
    "RNA",
    "LogNormalize_10000",
    "mean_per_gene_z_score",
    "pooled_qHSC_90th_percentile",
    paste(available_panel(ECM_COLLAGEN, obj), collapse = ";"),
    paste(available_panel(ECM_REMODELING, obj), collapse = ";"),
    "mean_per_gene_z_score",
    proliferation_reference,
    "reference_p95_AND_at_least_2_direct_core_genes_detected",
    "reference_p90_AND_at_least_1_direct_core_gene_detected",
    "MKI67_detected_OR_TOP2A_detected_OR_at_least_2_other_core_genes_detected"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "ECM_Proliferation_provenance_v1.2.tsv"
  )
)

# ------------------------------------------------------------
# 22. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "input_cells",
    "expected_cells",
    "cell_count_matches",
    "all_cells_have_state",
    "qHSC_present",
    "aHSC_present",
    "ECM_total_score_complete",
    "ECM_collagen_score_complete",
    "ECM_remodeling_score_complete",
    "proliferation_score_complete",
    "strict_cycling_defined_all_cells",
    "sensitive_cycling_defined_all_cells",
    "marker_only_defined_all_cells",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),

  value = c(
    as.character(ncol(obj)),
    "7640",
    as.character(ncol(obj) == 7640),
    as.character(all(!is.na(obj$HSC_state_v1.2))),
    as.character(any(obj$HSC_state_v1.2 == "qHSC")),
    as.character(any(obj$HSC_state_v1.2 == "aHSC_mesenchymal")),
    as.character(all(is.finite(obj$ECM_total_score_v1.2))),
    as.character(all(is.finite(obj$ECM_collagen_score_v1.2))),
    as.character(all(is.finite(obj$ECM_remodeling_score_v1.2))),
    as.character(all(is.finite(obj$Proliferation_score_v1.2))),
    as.character(all(!is.na(obj$Cycling_strict_PRIMARY_v1.2))),
    as.character(all(!is.na(obj$Cycling_sensitive_v1.2))),
    as.character(all(!is.na(obj$Cycling_marker_only_v1.2))),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.2.tsv"
  )
)

# ------------------------------------------------------------
# 23. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SAMPLE SUMMARY\n")
cat("============================================================\n")
print(sample_summary)

cat("\n")
cat("============================================================\n")
cat("STATE SUMMARY\n")
cat("============================================================\n")
print(sample_state_summary)

if (exists("NAS_summary")) {
  cat("\n")
  cat("============================================================\n")
  cat("NAS SUMMARY\n")
  cat("============================================================\n")
  print(NAS_summary)
}

if (exists("spearman")) {
  cat("\n")
  cat("============================================================\n")
  cat("NAS SPEARMAN\n")
  cat("============================================================\n")
  print(spearman)
}

cat("\n")
cat("============================================================\n")
cat("ECM AXIS CONCORDANCE\n")
cat("============================================================\n")
print(axis_concordance)

cat("\n")
cat("============================================================\n")
cat("DIRECT PROLIFERATION MARKER DETECTION BY SAMPLE\n")
cat("============================================================\n")
print(marker_detection_by_sample)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
