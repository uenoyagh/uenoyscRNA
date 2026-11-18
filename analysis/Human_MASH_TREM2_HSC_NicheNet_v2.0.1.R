#!/usr/bin/env Rscript

# ============================================================
# Human MASH TREM2-Mphi -> HSC NicheNet ligand-target analysis
# v2.0.1
#
# PURPOSE
#   Use the robust sample-level association established in v1.9:
#
#     TREM2-Mphi TGFB1 burden
#              ->
#     whole-HSC pseudobulk transcriptional response
#
#   to perform a sender-focused NicheNet analysis.
#
# DESIGN
#   - Sender: Lipid-associated/TREM2_Mphi
#   - Receiver: all canonical HSC cells (7,640 cells)
#   - Biological unit: sample (LT1-LT8)
#   - Receiver response: whole-HSC pseudobulk genes positively
#     associated with TREM2-Mphi TGFB1 burden across 8 samples
#   - NicheNet ligand ranking metric: aupr_corrected
#   - Primary question:
#       Is TGFB1 prioritized among TREM2-Mphi ligands, and which
#       HSC target genes are predicted downstream?
#
# IMPORTANT INTERPRETATION
#   This is hypothesis prioritization, not proof of causality.
#   The receiver gene set is adapted to a continuous sample-level
#   response rather than a binary condition DE contrast.
#
# NICHENET MODEL
#   Current human NicheNet v2 prior model resources from Zenodo:
#     lr_network_human_21122021.rds
#     ligand_target_matrix_nsga2r_final.rds
#
#   Official NicheNet repository:
#     https://github.com/saeyslab/nichenetr
#
# INPUTS
#   1. Canonical HSC receiver RDS
#   2. Canonical final myeloid annotation RDS
#   3. v1.5 sender burden table
#
# OUTPUT
#   SSD-only results; source RDS files are never overwritten.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(Matrix)
})

if (!requireNamespace("nichenetr", quietly = TRUE)) {
  stop(
    paste0(
      "\nPackage 'nichenetr' is not installed in this project.\n",
      "Install it from the uenoyscRNA project directory with:\n\n",
      "  R -q -e 'renv::install(\"saeyslab/nichenetr\")'\n\n",
      "Then rerun this script. Do not run renv::snapshot() until ",
      "the analysis completes successfully.\n"
    )
  )
}

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2",
  "objects",
  "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
)

MYELOID_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

SENDER_BURDEN_TABLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5",
  "tables",
  "02_TGFB1_sender_state_burden_by_sample.tsv"
)

RESOURCE_DIR <- file.path(
  ROOT,
  "resources",
  "NicheNet_v2_human"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_OBJ   <- file.path(OUT, "objects")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(RESOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH TREM2-Mphi -> HSC NicheNet v2.0.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(HSC_RDS))
stopifnot(file.exists(MYELOID_RDS))
stopifnot(file.exists(SENDER_BURDEN_TABLE))

# ------------------------------------------------------------
# 1. Download/cache official NicheNet v2 human resources
# ------------------------------------------------------------

LR_URL <- paste0(
  "https://zenodo.org/record/7074291/files/",
  "lr_network_human_21122021.rds"
)

LT_URL <- paste0(
  "https://zenodo.org/record/7074291/files/",
  "ligand_target_matrix_nsga2r_final.rds"
)

LR_FILE <- file.path(
  RESOURCE_DIR,
  "lr_network_human_21122021.rds"
)

LT_FILE <- file.path(
  RESOURCE_DIR,
  "ligand_target_matrix_nsga2r_final.rds"
)

validate_rds <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  ok <- tryCatch(
    {
      x <- readRDS(path)
      rm(x)
      gc(FALSE)
      TRUE
    },
    error = function(e) {
      message(
        "RDS validation failed for ",
        basename(path),
        ": ",
        conditionMessage(e)
      )
      FALSE
    }
  )

  ok
}

download_if_missing <- function(url, dest) {

  # Reuse only a file that can actually be read as RDS.
  if (validate_rds(dest)) {
    cat("Validated cached resource:\n", dest, "\n")
    return(invisible(TRUE))
  }

  # Keep an incomplete download as .part and allow curl to resume it.
  part <- paste0(dest, ".part")

  # If the old destination is incomplete, move it to .part so curl -C -
  # can resume from the existing byte count.
  if (file.exists(dest) && !validate_rds(dest)) {
    if (file.exists(part)) {
      file.remove(part)
    }
    ok_rename <- file.rename(dest, part)
    if (!ok_rename) {
      stop("Could not move incomplete resource to .part: ", dest)
    }
  }

  cat("Downloading/resuming:\n", url, "\n")
  cat("Temporary file:\n", part, "\n")

  curl_bin <- Sys.which("curl")

  if (!nzchar(curl_bin)) {
    # Fallback to R download.file with a long timeout.
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(3600, old_timeout))

    ok <- tryCatch(
      {
        download.file(
          url,
          destfile = part,
          mode = "wb",
          quiet = FALSE,
          method = "libcurl"
        )
        TRUE
      },
      error = function(e) {
        message(conditionMessage(e))
        FALSE
      }
    )
  } else {
    # macOS curl: resume partial file, follow redirects, retry transient failures,
    # and do not impose a total transfer timeout.
    args <- c(
      "--location",
      "--fail",
      "--retry", "20",
      "--retry-delay", "5",
      "--retry-all-errors",
      "--connect-timeout", "30",
      "--max-time", "0",
      "--continue-at", "-",
      "--output", part,
      url
    )

    status <- system2(
      curl_bin,
      args = args
    )

    ok <- identical(status, 0L)
  }

  if (!ok) {
    stop(
      "Download did not complete successfully: ",
      url,
      "\nPartial file retained for resume: ",
      part
    )
  }

  if (!validate_rds(part)) {
    stop(
      "Downloaded file is not a valid RDS: ",
      part
    )
  }

  if (file.exists(dest)) {
    file.remove(dest)
  }

  if (!file.rename(part, dest)) {
    stop(
      "Validated resource could not be moved into place: ",
      dest
    )
  }

  cat("Validated resource saved:\n", dest, "\n")
  invisible(TRUE)
}

download_if_missing(
  LR_URL,
  LR_FILE
)

download_if_missing(
  LT_URL,
  LT_FILE
)

lr_network <- readRDS(LR_FILE) %>%
  distinct(from, to)

ligand_target_matrix <- readRDS(LT_FILE)

cat("NicheNet LR interactions:", nrow(lr_network), "\n")
cat(
  "Ligand-target matrix:",
  nrow(ligand_target_matrix),
  "targets x",
  ncol(ligand_target_matrix),
  "ligands\n\n"
)

# ------------------------------------------------------------
# 2. Load biological objects
# ------------------------------------------------------------

hsc <- readRDS(HSC_RDS)
mphi <- readRDS(MYELOID_RDS)

DefaultAssay(hsc) <- "RNA"
DefaultAssay(mphi) <- "RNA"

if (ncol(hsc) != 7640) {
  warning(
    "Expected 7640 canonical HSC cells; observed ",
    ncol(hsc)
  )
}

required_hsc <- c(
  "sample",
  "NAS_score"
)

required_mphi <- c(
  "sample_v17",
  "NAS_v17",
  "myeloid_state_FINAL_v19"
)

if (length(setdiff(required_hsc, colnames(hsc[[]]))) > 0) {
  stop("Required HSC metadata missing.")
}

if (length(setdiff(required_mphi, colnames(mphi[[]]))) > 0) {
  stop("Required myeloid metadata missing.")
}

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

mphi <- NormalizeData(
  mphi,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 3. TREM2-Mphi TGFB1 burden by sample
# ------------------------------------------------------------

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

sender_burden <- read_tsv(
  SENDER_BURDEN_TABLE,
  show_col_types = FALSE
) %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  select(
    sample,
    NAS,
    TGFB1_burden_expr,
    TGFB1_expression_share,
    state_fraction,
    TGFB1_mean,
    TGFB1_detection
  ) %>%
  arrange(sample)

if (nrow(sender_burden) != 8) {
  stop(
    "Expected 8 TREM2-Mphi sample rows; observed ",
    nrow(sender_burden)
  )
}

# ------------------------------------------------------------
# 4. Whole-HSC pseudobulk counts
# ------------------------------------------------------------

hsc_counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

hsc_md <- hsc[[]]

samples <- sort(
  unique(
    as.character(
      hsc_md$sample
    )
  )
)

pb_counts_list <- lapply(
  samples,
  function(s) {

    idx <- which(
      as.character(
        hsc_md$sample
      ) == s
    )

    Matrix::rowSums(
      hsc_counts[
        ,
        idx,
        drop = FALSE
      ]
    )
  }
)

pb_counts <- do.call(
  cbind,
  pb_counts_list
)

rownames(pb_counts) <- rownames(hsc_counts)
colnames(pb_counts) <- samples

pb_lib <- colSums(pb_counts)

pb_cpm <- sweep(
  pb_counts,
  2,
  pb_lib,
  "/"
) * 1e6

pb_logcpm <- log2(
  pb_cpm + 1
)

pb_meta <- tibble(
  sample = samples
) %>%
  left_join(
    hsc_md %>%
      transmute(
        sample = as.character(sample),
        NAS = suppressWarnings(
          as.numeric(NAS_score)
        )
      ) %>%
      distinct(),
    by = "sample"
  ) %>%
  left_join(
    sender_burden,
    by = c(
      "sample",
      "NAS"
    )
  )

if (nrow(pb_meta) != 8) {
  stop("Whole-HSC pseudobulk metadata does not contain 8 samples.")
}

# ------------------------------------------------------------
# 5. Receiver expressed-gene background
#
# Primary background:
#   CPM >= 1 in at least 4 of 8 samples
#   AND target gene exists in NicheNet matrix
# ------------------------------------------------------------

receiver_expressed <- rownames(pb_cpm)[
  rowSums(
    pb_cpm >= 1,
    na.rm = TRUE
  ) >= 4
]

background_expressed_genes <- intersect(
  receiver_expressed,
  rownames(
    ligand_target_matrix
  )
)

cat(
  "Receiver expressed genes:",
  length(receiver_expressed),
  "\n"
)

cat(
  "NicheNet-compatible background genes:",
  length(background_expressed_genes),
  "\n\n"
)

# ------------------------------------------------------------
# 6. Unbiased whole-HSC gene association with TREM2 TGFB1 burden
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  keep <- complete.cases(x, y)

  if (
    sum(keep) < 4 ||
    length(unique(x[keep])) < 2 ||
    length(unique(y[keep])) < 2
  ) {
    return(
      c(
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

  c(
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

receiver_assoc_mat <- t(
  vapply(
    background_expressed_genes,
    function(g) {

      safe_spearman(
        pb_logcpm[
          g,
          pb_meta$sample
        ],
        pb_meta$TGFB1_burden_expr
      )
    },
    FUN.VALUE = c(
      rho = 0,
      p_value = 0
    )
  )
)

receiver_assoc <- as_tibble(
  receiver_assoc_mat,
  rownames = "gene"
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    desc(rho),
    p_value
  )

write_tsv(
  receiver_assoc,
  file.path(
    DIR_TABLE,
    "01_wholeHSC_gene_association_with_TREM2_TGFB1.tsv"
  )
)

# ------------------------------------------------------------
# 7. Define NicheNet receiver gene set of interest
#
# PRIMARY:
#   rho >= 0.60 AND nominal P <= 0.10
#
# This threshold is prespecified for n=8 exploratory sample-level
# response discovery and is NOT claimed as genome-wide significance.
#
# If >500 genes satisfy it, keep the strongest 500.
# Require at least 20 genes to run NicheNet.
# ------------------------------------------------------------

geneset_primary <- receiver_assoc %>%
  filter(
    rho >= 0.60,
    p_value <= 0.10
  ) %>%
  arrange(
    desc(rho),
    p_value
  ) %>%
  slice_head(n = 500) %>%
  pull(gene)

geneset_stringent <- receiver_assoc %>%
  filter(
    rho >= 0.70,
    p_value <= 0.05
  ) %>%
  arrange(
    desc(rho),
    p_value
  ) %>%
  slice_head(n = 500) %>%
  pull(gene)

# Sensitivity: top 10% positively correlated genes
positive_assoc <- receiver_assoc %>%
  filter(
    is.finite(rho),
    rho > 0
  ) %>%
  arrange(
    desc(rho)
  )

n_broad <- max(
  20,
  ceiling(
    0.10 *
      nrow(
        positive_assoc
      )
  )
)

geneset_broad <- positive_assoc %>%
  slice_head(
    n = min(
      n_broad,
      500
    )
  ) %>%
  pull(gene)

geneset_summary <- tibble(
  geneset = c(
    "PRIMARY_rho0.60_p0.10",
    "STRINGENT_rho0.70_p0.05",
    "SENSITIVITY_top10pct_positive_max500"
  ),
  n_genes = c(
    length(geneset_primary),
    length(geneset_stringent),
    length(geneset_broad)
  ),
  definition = c(
    "rho>=0.60 & nominal p<=0.10; max500",
    "rho>=0.70 & nominal p<=0.05; max500",
    "top10pct positive rho; min20 max500"
  )
)

write_tsv(
  geneset_summary,
  file.path(
    DIR_TABLE,
    "02_receiver_geneset_definitions.tsv"
  )
)

write_tsv(
  tibble(
    gene = geneset_primary
  ),
  file.path(
    DIR_TABLE,
    "03_PRIMARY_receiver_geneset.tsv"
  )
)

write_tsv(
  tibble(
    gene = geneset_stringent
  ),
  file.path(
    DIR_TABLE,
    "04_STRINGENT_receiver_geneset.tsv"
  )
)

write_tsv(
  tibble(
    gene = geneset_broad
  ),
  file.path(
    DIR_TABLE,
    "05_SENSITIVITY_receiver_geneset.tsv"
  )
)

if (length(geneset_primary) < 20) {
  stop(
    "Primary receiver gene set has only ",
    length(geneset_primary),
    " genes. Do not force NicheNet; inspect table 01 first."
  )
}

# ------------------------------------------------------------
# 8. Sender ligand expression in TREM2-Mphi
#
# Sample-aware criterion:
#   detection fraction >= 0.05 in at least 2 of 8 samples
# ------------------------------------------------------------

mphi_md <- mphi[[]]

trem2_cells <- rownames(
  mphi_md
)[
  as.character(
    mphi_md$myeloid_state_FINAL_v19
  ) == TREM2_STATE
]

if (length(trem2_cells) == 0) {
  stop("No TREM2-Mphi cells found.")
}

ligand_universe <- intersect(
  unique(
    lr_network$from
  ),
  rownames(mphi)
)

sender_ligand_expr <- FetchData(
  mphi,
  vars = ligand_universe,
  cells = trem2_cells,
  layer = "data"
)

sender_ligand_counts <- FetchData(
  mphi,
  vars = ligand_universe,
  cells = trem2_cells,
  layer = "counts"
)

sender_cell_meta <- mphi_md[
  trem2_cells,
  ,
  drop = FALSE
] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(
      sample_v17
    )
  )

sender_detect_long <- as_tibble(
  sender_ligand_counts,
  rownames = "cell"
) %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "count"
  ) %>%
  left_join(
    sender_cell_meta,
    by = "cell"
  ) %>%
  group_by(
    sample,
    ligand
  ) %>%
  summarise(
    sender_cells = n(),
    detection_fraction =
      mean(
        count > 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

sender_expression_audit <- sender_detect_long %>%
  group_by(ligand) %>%
  summarise(
    samples_present = n(),
    samples_detection_ge_0.05 =
      sum(
        detection_fraction >= 0.05,
        na.rm = TRUE
      ),
    mean_detection_fraction =
      mean(
        detection_fraction,
        na.rm = TRUE
      ),
    max_detection_fraction =
      max(
        detection_fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

sender_expressed_ligands <- sender_expression_audit %>%
  filter(
    samples_detection_ge_0.05 >= 2
  ) %>%
  pull(ligand)

# ------------------------------------------------------------
# 9. Receiver receptor expression
#
# Receptor detection >=5% in at least 4 of 8 HSC samples
# ------------------------------------------------------------

receptor_universe <- intersect(
  unique(
    lr_network$to
  ),
  rownames(hsc)
)

hsc_receptor_counts <- FetchData(
  hsc,
  vars = receptor_universe,
  layer = "counts"
)

hsc_cell_meta <- hsc[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample)
  )

receiver_detect_long <- as_tibble(
  hsc_receptor_counts,
  rownames = "cell"
) %>%
  pivot_longer(
    -cell,
    names_to = "receptor",
    values_to = "count"
  ) %>%
  left_join(
    hsc_cell_meta,
    by = "cell"
  ) %>%
  group_by(
    sample,
    receptor
  ) %>%
  summarise(
    receiver_cells = n(),
    detection_fraction =
      mean(
        count > 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

receiver_expression_audit <- receiver_detect_long %>%
  group_by(receptor) %>%
  summarise(
    samples_present = n(),
    samples_detection_ge_0.05 =
      sum(
        detection_fraction >= 0.05,
        na.rm = TRUE
      ),
    mean_detection_fraction =
      mean(
        detection_fraction,
        na.rm = TRUE
      ),
    max_detection_fraction =
      max(
        detection_fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

receiver_expressed_receptors <- receiver_expression_audit %>%
  filter(
    samples_detection_ge_0.05 >= 4
  ) %>%
  pull(receptor)

# ------------------------------------------------------------
# 10. Potential ligands:
#     sender-expressed ligand + receiver-expressed receptor
# ------------------------------------------------------------

lr_expressed <- lr_network %>%
  filter(
    from %in% sender_expressed_ligands,
    to %in% receiver_expressed_receptors
  ) %>%
  filter(
    from %in%
      colnames(
        ligand_target_matrix
      )
  )

potential_ligands <- sort(
  unique(
    lr_expressed$from
  )
)

write_tsv(
  sender_expression_audit %>%
    arrange(
      desc(
        samples_detection_ge_0.05
      ),
      desc(
        mean_detection_fraction
      )
    ),
  file.path(
    DIR_TABLE,
    "06_TREM2_sender_ligand_expression_audit.tsv"
  )
)

write_tsv(
  receiver_expression_audit %>%
    arrange(
      desc(
        samples_detection_ge_0.05
      ),
      desc(
        mean_detection_fraction
      )
    ),
  file.path(
    DIR_TABLE,
    "07_HSC_receptor_expression_audit.tsv"
  )
)

write_tsv(
  lr_expressed,
  file.path(
    DIR_TABLE,
    "08_expressed_ligand_receptor_pairs.tsv"
  )
)

write_tsv(
  tibble(
    ligand = potential_ligands
  ),
  file.path(
    DIR_TABLE,
    "09_potential_ligands.tsv"
  )
)

cat(
  "Sender-expressed ligands:",
  length(sender_expressed_ligands),
  "\n"
)

cat(
  "Receiver-expressed receptors:",
  length(receiver_expressed_receptors),
  "\n"
)

cat(
  "Potential NicheNet ligands:",
  length(potential_ligands),
  "\n\n"
)

if (!("TGFB1" %in% potential_ligands)) {
  warning(
    "TGFB1 is NOT in the potential ligand set under the primary ",
    "expression thresholds. Interpret before changing thresholds."
  )
}

if (length(potential_ligands) < 5) {
  stop(
    "Too few potential ligands for meaningful NicheNet ranking."
  )
}

# ------------------------------------------------------------
# 11. Run NicheNet ligand activity
# ------------------------------------------------------------

run_nichenet <- function(
  geneset,
  label
) {

  geneset <- intersect(
    geneset,
    background_expressed_genes
  )

  out <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes =
      background_expressed_genes,
    ligand_target_matrix =
      ligand_target_matrix,
    potential_ligands =
      potential_ligands
  ) %>%
    arrange(
      desc(
        aupr_corrected
      )
    ) %>%
    mutate(
      rank = row_number(),
      geneset = label
    )

  out
}

activities_primary <- run_nichenet(
  geneset_primary,
  "PRIMARY_rho0.60_p0.10"
)

activities_stringent <- if (
  length(geneset_stringent) >= 20
) {
  run_nichenet(
    geneset_stringent,
    "STRINGENT_rho0.70_p0.05"
  )
} else {
  tibble()
}

activities_broad <- run_nichenet(
  geneset_broad,
  "SENSITIVITY_top10pct_positive"
)

activities_all <- bind_rows(
  activities_primary,
  activities_stringent,
  activities_broad
)

write_tsv(
  activities_all,
  file.path(
    DIR_TABLE,
    "10_NicheNet_ligand_activities_all_genesets.tsv"
  )
)

# ------------------------------------------------------------
# 12. TGFB1 summary and rank stability
# ------------------------------------------------------------

tgfb1_summary <- activities_all %>%
  filter(
    test_ligand == "TGFB1"
  ) %>%
  select(
    geneset,
    test_ligand,
    auroc,
    aupr,
    aupr_corrected,
    pearson,
    rank
  )

write_tsv(
  tgfb1_summary,
  file.path(
    DIR_TABLE,
    "11_TGFB1_NicheNet_activity_summary.tsv"
  )
)

# ------------------------------------------------------------
# 13. Top ligands, targets, and TGFB1 targets
# ------------------------------------------------------------

top_ligands_primary <- activities_primary %>%
  slice_head(n = 30) %>%
  pull(test_ligand)

top_targets_list <- lapply(
  top_ligands_primary,
  function(lig) {

    x <- nichenetr::get_weighted_ligand_target_links(
      ligand = lig,
      geneset = geneset_primary,
      ligand_target_matrix =
        ligand_target_matrix,
      n = 250
    )

    x %>%
      mutate(
        ligand_rank =
          match(
            lig,
            top_ligands_primary
          )
      )
  }
)

top_ligand_targets <- bind_rows(
  top_targets_list
) %>%
  arrange(
    ligand_rank,
    desc(weight)
  )

write_tsv(
  activities_primary %>%
    slice_head(n = 30),
  file.path(
    DIR_TABLE,
    "12_PRIMARY_top30_ligands.tsv"
  )
)

write_tsv(
  top_ligand_targets,
  file.path(
    DIR_TABLE,
    "13_PRIMARY_top30_ligand_target_links.tsv"
  )
)

if ("TGFB1" %in% potential_ligands) {

  tgfb1_targets <- nichenetr::get_weighted_ligand_target_links(
    ligand = "TGFB1",
    geneset = geneset_primary,
    ligand_target_matrix =
      ligand_target_matrix,
    n = 250
  ) %>%
    arrange(
      desc(weight)
    ) %>%
    left_join(
      receiver_assoc,
      by = c(
        "target" = "gene"
      )
    )

} else {

  tgfb1_targets <- tibble(
    ligand = character(),
    target = character(),
    weight = numeric(),
    rho = numeric(),
    p_value = numeric(),
    FDR = numeric()
  )
}

write_tsv(
  tgfb1_targets,
  file.path(
    DIR_TABLE,
    "14_PRIMARY_TGFB1_predicted_targets.tsv"
  )
)

# ------------------------------------------------------------
# 14. TGFB1 receptor support in HSC
# ------------------------------------------------------------

tgfb1_lr <- lr_network %>%
  filter(
    from == "TGFB1"
  ) %>%
  distinct(
    from,
    to
  ) %>%
  left_join(
    receiver_expression_audit,
    by = c(
      "to" = "receptor"
    )
  ) %>%
  mutate(
    receptor_primary_expressed =
      to %in%
      receiver_expressed_receptors
  ) %>%
  arrange(
    desc(
      receptor_primary_expressed
    ),
    desc(
      mean_detection_fraction
    )
  )

write_tsv(
  tgfb1_lr,
  file.path(
    DIR_TABLE,
    "15_TGFB1_receptor_support_in_HSC.tsv"
  )
)

# ------------------------------------------------------------
# 15. Activity-rank figure
# ------------------------------------------------------------

plot_df <- activities_primary %>%
  slice_head(n = 30) %>%
  mutate(
    ligand = factor(
      test_ligand,
      levels = rev(
        test_ligand
      )
    )
  )

p1 <- ggplot(
  plot_df,
  aes(
    x = ligand,
    y = aupr_corrected
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "NicheNet ligand activity: TREM2-Mphi -> whole HSC",
    subtitle = "Primary receiver response: genes positively associated with TREM2 TGFB1 burden",
    x = NULL,
    y = "AUPR corrected"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_PRIMARY_top30_NicheNet_ligand_activity.pdf"
  ),
  p1,
  width = 8,
  height = 8
)

# TGFB1 target regulatory potential
if (nrow(tgfb1_targets) > 0) {

  p2_df <- tgfb1_targets %>%
    filter(
      !is.na(target),
      !is.na(weight)
    ) %>%
    slice_head(n = 30) %>%
    mutate(
      target = factor(
        target,
        levels = rev(
          target
        )
      )
    )

  if (nrow(p2_df) > 0) {

    p2 <- ggplot(
      p2_df,
      aes(
        x = target,
        y = weight
      )
    ) +
      geom_col() +
      coord_flip() +
      theme_classic(base_size = 11) +
      labs(
        title = "Predicted TGFB1 target genes in HSC response",
        x = NULL,
        y = "NicheNet regulatory potential"
      )

    ggsave(
      file.path(
        DIR_FIG,
        "02_PRIMARY_TGFB1_top_targets.pdf"
      ),
      p2,
      width = 8,
      height = 8
    )
  }
}

# ------------------------------------------------------------
# 16. Save compact result object
# ------------------------------------------------------------

saveRDS(
  list(
    pb_meta = pb_meta,
    receiver_assoc = receiver_assoc,
    geneset_summary = geneset_summary,
    geneset_primary = geneset_primary,
    geneset_stringent = geneset_stringent,
    geneset_broad = geneset_broad,
    sender_expression_audit = sender_expression_audit,
    receiver_expression_audit = receiver_expression_audit,
    lr_expressed = lr_expressed,
    potential_ligands = potential_ligands,
    activities_primary = activities_primary,
    activities_stringent = activities_stringent,
    activities_broad = activities_broad,
    tgfb1_summary = tgfb1_summary,
    tgfb1_targets = tgfb1_targets,
    tgfb1_lr = tgfb1_lr
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_TREM2_HSC_NicheNet_v2.0.1_results.rds"
  )
)

# ------------------------------------------------------------
# 17. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "HSC_source_RDS",
    "myeloid_source_RDS",
    "sender_burden_table",
    "HSC_cells",
    "sender_state",
    "receiver_population",
    "receiver_response_definition",
    "receiver_background_definition",
    "sender_ligand_expression_definition",
    "receiver_receptor_expression_definition",
    "NicheNet_package_version",
    "NicheNet_lr_resource",
    "NicheNet_ligand_target_resource",
    "ligand_ranking_metric",
    "association_unit",
    "source_overwritten",
    "reclustering"
  ),
  value = c(
    "v2.0.1",
    HSC_RDS,
    MYELOID_RDS,
    SENDER_BURDEN_TABLE,
    as.character(
      ncol(hsc)
    ),
    TREM2_STATE,
    "all canonical HSC",
    "whole-HSC pseudobulk genes with Spearman rho>=0.60 and nominal p<=0.10 vs TREM2 TGFB1 burden; max500",
    "CPM>=1 in at least 4/8 whole-HSC pseudobulk samples and present in NicheNet ligand-target matrix",
    "TREM2-Mphi ligand detection>=5% in at least 2/8 samples",
    "HSC receptor detection>=5% in at least 4/8 samples",
    as.character(
      packageVersion(
        "nichenetr"
      )
    ),
    LR_FILE,
    LT_FILE,
    "aupr_corrected",
    "biological_sample",
    "FALSE",
    "none"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "NicheNet_provenance_v2.0.1.tsv"
  )
)

# ------------------------------------------------------------
# 18. Integrity
# ------------------------------------------------------------

tgfb1_rank_primary <- activities_primary %>%
  filter(
    test_ligand == "TGFB1"
  ) %>%
  pull(rank)

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "all_8_samples_present",
    "receiver_background_genes_n",
    "primary_receiver_geneset_n",
    "stringent_receiver_geneset_n",
    "potential_ligands_n",
    "TGFB1_in_potential_ligands",
    "TGFB1_primary_rank",
    "NicheNet_completed",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      ncol(hsc)
    ),
    "7640",
    as.character(
      ncol(hsc) == 7640
    ),
    as.character(
      nrow(pb_meta) == 8
    ),
    as.character(
      length(
        background_expressed_genes
      )
    ),
    as.character(
      length(
        geneset_primary
      )
    ),
    as.character(
      length(
        geneset_stringent
      )
    ),
    as.character(
      length(
        potential_ligands
      )
    ),
    as.character(
      "TGFB1" %in%
        potential_ligands
    ),
    ifelse(
      length(tgfb1_rank_primary) == 1,
      as.character(
        tgfb1_rank_primary
      ),
      "NA"
    ),
    "TRUE",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v2.0.1.tsv"
  )
)

# ------------------------------------------------------------
# 19. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("RECEIVER GENE SETS\n")
cat("============================================================\n")
print(geneset_summary)

cat("\n")
cat("============================================================\n")
cat("POTENTIAL LIGANDS\n")
cat("============================================================\n")
cat("N =", length(potential_ligands), "\n")
cat(
  "TGFB1 included:",
  "TGFB1" %in% potential_ligands,
  "\n"
)

cat("\n")
cat("============================================================\n")
cat("PRIMARY NICHE-NET TOP 30 LIGANDS\n")
cat("============================================================\n")
print(
  activities_primary %>%
    slice_head(n = 30)
)

cat("\n")
cat("============================================================\n")
cat("TGFB1 ACTIVITY SUMMARY\n")
cat("============================================================\n")
print(tgfb1_summary)

cat("\n")
cat("============================================================\n")
cat("TGFB1 RECEPTOR SUPPORT\n")
cat("============================================================\n")
print(tgfb1_lr)

cat("\n")
cat("============================================================\n")
cat("TOP TGFB1 TARGETS\n")
cat("============================================================\n")
print(
  tgfb1_targets %>%
    head(40)
)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
