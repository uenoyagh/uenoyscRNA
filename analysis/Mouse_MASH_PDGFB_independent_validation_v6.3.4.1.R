#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6340)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Independent PDGF-response validation in HSC
#
# Version: v6.3.4.1.1
#
# PURPOSE
#   Independently validate the PDGFB hypothesis without reusing the
#   NicheNet-selected target genes that were derived from the same receiver DE.
#
# PRIMARY RECEIVER
#   ECM-activated HSC
#
# BIOLOGICAL SAMPLES
#   Sham1, Sham20, Tx17, Tx5
#
# INPUTS
#   1) v6.2.0 interaction-ready Seurat object
#   2) v6.3.0 receiver pseudobulk DE table
#   3) v6.3.2 NicheNet-derived PDGFB target list
#      -- used ONLY as an exclusion/audit list, never as a positive gene set.
#   4) v6.3.3.1 sample-level weighted Pdgfb output table
#
# INDEPENDENT A-PRIORI MODULES
#
#   PDGF_IE_MAPK
#     Immediate-early / MAPK-responsive transcriptional program.
#
#   PDGF_GROWTH_CELL_CYCLE
#     Growth and cell-cycle response downstream of mitogenic signaling.
#
#   HSC_ACTIVATION_ECM
#     Fibrogenic / activated-HSC endpoint program.
#
#   HSC_CONTRACTILE
#     Contractile / myofibroblastic endpoint program.
#
# IMPORTANT
#   - Any gene appearing in the v6.3.2 NicheNet target list is removed.
#   - Module scores are computed from sample-level pseudobulk expression.
#   - No cell-level pseudoreplication.
#   - No NicheNet rerun.
#   - No CellChat rerun.
#   - Biological n = 2 Sham vs n = 2 Tx.
#   - Results remain exploratory.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(
  p,
  file,
  width,
  height
) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}

safe_join_rna <- function(object) {

  if (
    length(
      Layers(
        object[["RNA"]]
      )
    ) > 1
  ) {
    object[["RNA"]] <- JoinLayers(
      object[["RNA"]]
    )
  }

  object
}

canonical_condition <- function(sample) {

  ifelse(
    grepl(
      "^Tx",
      sample
    ),
    "Tx",
    "Sham"
  )
}

row_z <- function(x) {

  sx <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      sx
    ) ||
    sx ==
      0
  ) {
    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / sx
}

sample_pseudobulk <- function(
  counts,
  meta,
  sample_col,
  group_col,
  group_value,
  samples
) {

  mats <- lapply(
    samples,
    function(smp) {

      cells <- rownames(
        meta
      )[
        as.character(
          meta[[
            sample_col
          ]]
        ) ==
          smp &
          as.character(
            meta[[
              group_col
            ]]
          ) ==
            group_value
      ]

      if (
        !length(
          cells
        )
      ) {
        stop(
          "No cells for ",
          group_value,
          " / ",
          smp
        )
      }

      Matrix::rowSums(
        counts[
          ,
          cells,
          drop = FALSE
        ]
      )
    }
  )

  mat <- do.call(
    cbind,
    mats
  )

  rownames(
    mat
  ) <- rownames(
    counts
  )

  colnames(
    mat
  ) <- samples

  mat
}

pairwise_direction <- function(
  sham1,
  sham20,
  tx17,
  tx5,
  tol = 1e-12
) {

  diffs <- c(
    tx17 - sham1,
    tx17 - sham20,
    tx5 - sham1,
    tx5 - sham20
  )

  up_n <- sum(
    diffs >
      tol
  )

  down_n <- sum(
    diffs <
      -tol
  )

  tie_n <- 4L -
    up_n -
    down_n

  direction <- dplyr::case_when(
    up_n >
      down_n ~
      "Tx_up",

    down_n >
      up_n ~
      "Tx_down",

    TRUE ~
      "Mixed_or_tied"
  )

  tibble(
    Tx17_minus_Sham1 =
      diffs[[1]],
    Tx17_minus_Sham20 =
      diffs[[2]],
    Tx5_minus_Sham1 =
      diffs[[3]],
    Tx5_minus_Sham20 =
      diffs[[4]],
    pairwise_up_n =
      up_n,
    pairwise_down_n =
      down_n,
    pairwise_tie_n =
      tie_n,
    pairwise_direction =
      direction,
    pairwise_direction_consistency =
      max(
        up_n,
        down_n
      ) / 4
  )
}

safe_spearman <- function(
  x,
  y
) {

  ok <- is.finite(
    x
  ) &
    is.finite(
      y
    )

  if (
    sum(
      ok
    ) <
      3
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    stats::cor(
      x[
        ok
      ],
      y[
        ok
      ],
      method = "spearman"
    )
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_interaction_ready_v6.2.0",
  "RDS",
  "Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds"
)

V630 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_v6.3.0"
)

RECEIVER_DE_FILE <- file.path(
  V630,
  "Tables",
  "04_receiver_DE_all_HSC_states_v6.3.0.csv"
)

V632 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_axis_validation_v6.3.2"
)

NICHE_TARGET_FILE <- file.path(
  V632,
  "Tables",
  "10_PDGFB_target_gene_set_v6.3.2.csv"
)

V6331 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_population_weighted_v6.3.3.1"
)

WEIGHTED_FILE <- file.path(
  V6331,
  "Tables",
  "09_integrated_PDGFB_sender_receiver_by_sample_v6.3.3.1.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_independent_validation_v6.3.4.1"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

for (
  d in c(
    OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT,
    RDS_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

GROUP_COL <-
  "interaction_celltype_v620"

SAMPLE_COL <-
  "sample_interaction_v620"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

HSC3 <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

FOCAL_RECEIVER <-
  "ECM-activated HSC"

CONDITION_COLORS <- c(
  "Sham" = "#0072B2",
  "Tx" = "#D55E00"
)


# ==============================================================================
# 4. A-priori independent gene sets
#
# These genes are defined before looking at the v6.3.4.1 result.
# NicheNet-selected targets are explicitly excluded below.
# ==============================================================================

GENESETS_RAW <- list(

  PDGF_IE_MAPK = c(
    "Egr1",
    "Egr2",
    "Egr3",
    "Fos",
    "Fosb",
    "Jun",
    "Junb",
    "Jund",
    "Fosl1",
    "Dusp1",
    "Dusp4",
    "Dusp5",
    "Dusp6",
    "Spry1",
    "Spry2",
    "Spry4",
    "Etv4",
    "Etv5",
    "Myc"
  ),

  PDGF_GROWTH_CELL_CYCLE = c(
    "Myc",
    "Ccnd1",
    "Ccnd2",
    "Cdk4",
    "Cdk6",
    "E2f1",
    "E2f2",
    "Pcna",
    "Mcm2",
    "Mcm3",
    "Mcm4",
    "Mcm5",
    "Mcm6",
    "Mcm7",
    "Mki67",
    "Top2a"
  ),

  HSC_ACTIVATION_ECM = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Col5a1",
    "Col5a2",
    "Fn1",
    "Lox",
    "Loxl2",
    "Timp1",
    "Ctgf",
    "Serpine1",
    "Acta2",
    "Tagln",
    "Tpm2",
    "Myl9"
  ),

  HSC_CONTRACTILE = c(
    "Acta2",
    "Tagln",
    "Myl9",
    "Tpm2",
    "Cnn1",
    "Myh11",
    "Des",
    "Vim"
  )
)


# ==============================================================================
# 5. Preflight / load
# ==============================================================================

for (
  f in c(
    INPUT_RDS,
    RECEIVER_DE_FILE,
    NICHE_TARGET_FILE,
    WEIGHTED_FILE
  )
) {
  if (
    !file.exists(
      f
    )
  ) {
    stop(
      "Required input missing: ",
      f
    )
  }
}

msg(
  "Loading interaction-ready object..."
)

obj <- readRDS(
  INPUT_RDS
)

DefaultAssay(
  obj
) <- "RNA"

obj <- safe_join_rna(
  obj
)

required_meta <- c(
  GROUP_COL,
  SAMPLE_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    obj@meta.data
  )
)

if (
  length(
    missing_meta
  )
) {
  stop(
    "Missing metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

meta <- obj@meta.data

receiver_de <- read.csv(
  RECEIVER_DE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

niche_targets <- read.csv(
  NICHE_TARGET_FILE,
  check.names = FALSE
) %>%
  as_tibble()

weighted_axis <- read.csv(
  WEIGHTED_FILE,
  check.names = FALSE
) %>%
  as_tibble()

if (
  !"target" %in%
    colnames(
      niche_targets
    )
) {
  stop(
    "NicheNet target file lacks 'target' column."
  )
}

NICHE_EXCLUSION <- unique(
  as.character(
    niche_targets$target
  )
)


# ==============================================================================
# 6. Independence audit and final gene sets
# ==============================================================================

msg(
  "Auditing gene-set independence from NicheNet targets..."
)

gene_audit <- bind_rows(
  lapply(
    names(
      GENESETS_RAW
    ),
    function(gs) {

      genes <- unique(
        GENESETS_RAW[[
          gs
        ]]
      )

      tibble(
        module =
          gs,
        gene =
          genes,
        overlaps_NicheNet_target =
          gene %in%
            NICHE_EXCLUSION,
        present_in_RNA =
          gene %in%
            rownames(
              counts
            )
      )
    }
  )
)

write.csv(
  gene_audit,
  file.path(
    TAB_OUT,
    "01_independent_gene_set_audit_v6.3.4.1.csv"
  ),
  row.names = FALSE
)

GENESETS <- lapply(
  names(
    GENESETS_RAW
  ),
  function(gs) {

    setdiff(
      unique(
        GENESETS_RAW[[
          gs
        ]]
      ),
      NICHE_EXCLUSION
    ) %>%
      intersect(
        rownames(
          counts
        )
      )
  }
)

names(
  GENESETS
) <- names(
  GENESETS_RAW
)

geneset_summary <- tibble(
  module =
    names(
      GENESETS
    ),
  n_raw =
    vapply(
      GENESETS_RAW,
      length,
      integer(
        1
      )
    ),
  n_overlap_removed =
    vapply(
      names(
        GENESETS_RAW
      ),
      function(gs) {
        length(
          intersect(
            GENESETS_RAW[[
              gs
            ]],
            NICHE_EXCLUSION
          )
        )
      },
      integer(
        1
      )
    ),
  n_present_final =
    vapply(
      GENESETS,
      length,
      integer(
        1
      )
    )
)

write.csv(
  geneset_summary,
  file.path(
    TAB_OUT,
    "02_independent_gene_set_summary_v6.3.4.1.csv"
  ),
  row.names = FALSE
)

if (
  any(
    geneset_summary$n_present_final <
      4
  )
) {
  warning(
    "At least one independent module has fewer than 4 genes present."
  )
}

final_gene_table <- bind_rows(
  lapply(
    names(
      GENESETS
    ),
    function(gs) {

      tibble(
        module =
          gs,
        gene =
          GENESETS[[
            gs
          ]]
      )
    }
  )
)

write.csv(
  final_gene_table,
  file.path(
    TAB_OUT,
    "03_final_independent_gene_sets_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Sample-level pseudobulk for all three HSC states
# ==============================================================================

msg(
  "Computing HSC-state pseudobulk..."
)

all_module_genes <- unique(
  unlist(
    GENESETS,
    use.names = FALSE
  )
)

pb_values_list <- list()

for (
  hsc_state in HSC3
) {

  pb <- sample_pseudobulk(
    counts = counts,
    meta = meta,
    sample_col = SAMPLE_COL,
    group_col = GROUP_COL,
    group_value = hsc_state,
    samples = SAMPLES
  )

  lib <- Matrix::colSums(
    pb
  )

  genes_present <- intersect(
    all_module_genes,
    rownames(
      pb
    )
  )

  cpm <- sapply(
    SAMPLES,
    function(smp) {

      1e6 *
        as.numeric(
          pb[
            genes_present,
            smp,
            drop = TRUE
          ]
        ) /
        max(
          lib[[
            smp
          ]],
          1
        )
    }
  )

  if (
    is.null(
      dim(
        cpm
      )
    )
  ) {
    cpm <- matrix(
      cpm,
      nrow =
        length(
          genes_present
        ),
      dimnames = list(
        genes_present,
        SAMPLES
      )
    )
  } else {
    rownames(
      cpm
    ) <- genes_present
    colnames(
      cpm
    ) <- SAMPLES
  }

  logcpm <- log2(
    cpm +
      1
  )

  zmat <- t(
    apply(
      logcpm,
      1,
      row_z
    )
  )

  # apply(..., MARGIN = 1) does not reliably preserve the sample names
  # after transpose. Re-attach dimnames explicitly before pivot_longer().
  rownames(
    zmat
  ) <- genes_present

  colnames(
    zmat
  ) <- SAMPLES

  if (
    !all(
      SAMPLES %in%
        colnames(
          logcpm
        )
    )
  ) {
    stop(
      "logcpm sample columns missing for HSC state: ",
      hsc_state,
      ". Found: ",
      paste(
        colnames(
          logcpm
        ),
        collapse = ", "
      )
    )
  }

  if (
    !all(
      SAMPLES %in%
        colnames(
          zmat
        )
    )
  ) {
    stop(
      "zmat sample columns missing for HSC state: ",
      hsc_state,
      ". Found: ",
      paste(
        colnames(
          zmat
        ),
        collapse = ", "
      )
    )
  }

  logcpm_long <- as.data.frame(
    logcpm,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    pivot_longer(
      cols =
        all_of(
          SAMPLES
        ),
      names_to =
        "sample",
      values_to =
        "log2_CPM1"
    )

  zmat_long <- as.data.frame(
    zmat,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    pivot_longer(
      cols =
        all_of(
          SAMPLES
        ),
      names_to =
        "sample",
      values_to =
        "gene_z"
    )

  pb_values_list[[
    hsc_state
  ]] <- logcpm_long %>%
    left_join(
      zmat_long,
      by = c(
        "gene",
        "sample"
      )
    ) %>%
    mutate(
      HSC_state =
        hsc_state,
      condition =
        canonical_condition(
          sample
        )
    )

  msg(
    "  ",
    hsc_state,
    ": ",
    length(
      genes_present
    ),
    " genes x ",
    length(
      SAMPLES
    ),
    " samples; z-score matrix dimnames OK"
  )
}

pb_values <- bind_rows(
  pb_values_list
)

write.csv(
  pb_values,
  file.path(
    TAB_OUT,
    "04_independent_module_gene_pseudobulk_all_HSC_states_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Module scores
#
# Primary score = mean gene-wise z-score.
# Each gene therefore contributes equally regardless of absolute abundance.
# ==============================================================================

module_score_list <- list()

for (
  gs in names(
    GENESETS
  )
) {

  genes <- GENESETS[[
    gs
  ]]

  tmp <- pb_values %>%
    filter(
      gene %in%
        genes
    ) %>%
    group_by(
      HSC_state,
      sample,
      condition
    ) %>%
    summarise(
      module_score_z =
        mean(
          gene_z,
          na.rm = TRUE
        ),
      mean_log2_CPM1 =
        mean(
          log2_CPM1,
          na.rm = TRUE
        ),
      n_genes =
        n_distinct(
          gene
        ),
      .groups = "drop"
    ) %>%
    mutate(
      module =
        gs
    )

  module_score_list[[
    gs
  ]] <- tmp
}

module_scores <- bind_rows(
  module_score_list
)

write.csv(
  module_scores,
  file.path(
    TAB_OUT,
    "05_independent_module_scores_by_HSC_state_sample_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. ECM-activated HSC replicate-aware module comparison
# ==============================================================================

ecm_scores <- module_scores %>%
  filter(
    HSC_state ==
      FOCAL_RECEIVER
  )

ecm_wide <- ecm_scores %>%
  select(
    module,
    sample,
    module_score_z
  ) %>%
  pivot_wider(
    names_from =
      sample,
    values_from =
      module_score_z
  )

module_metrics <- lapply(
  seq_len(
    nrow(
      ecm_wide
    )
  ),
  function(i) {

    pairwise_direction(
      sham1 =
        ecm_wide$Sham1[[i]],
      sham20 =
        ecm_wide$Sham20[[i]],
      tx17 =
        ecm_wide$Tx17[[i]],
      tx5 =
        ecm_wide$Tx5[[i]]
    )
  }
) %>%
  bind_rows()

ecm_comparison <- bind_cols(
  ecm_wide,
  module_metrics
) %>%
  mutate(
    Sham_mean =
      (
        Sham1 +
          Sham20
      ) / 2,
    Tx_mean =
      (
        Tx17 +
          Tx5
      ) / 2,
    Tx_minus_Sham =
      Tx_mean -
        Sham_mean,
    both_Tx_below_both_Sham =
      max(
        Tx17,
        Tx5
      ) <
        min(
          Sham1,
          Sham20
        ),
    both_Tx_above_both_Sham =
      min(
        Tx17,
        Tx5
      ) >
        max(
          Sham1,
          Sham20
        )
  ) %>%
  arrange(
    Tx_minus_Sham
  )

write.csv(
  ecm_comparison,
  file.path(
    TAB_OUT,
    "06_ECM_HSC_independent_module_Sham_vs_Tx_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Independent gene-level DE audit from v6.3.0 receiver DE
# ==============================================================================

de_audit <- receiver_de %>%
  filter(
    receiver ==
      FOCAL_RECEIVER,
    gene %in%
      all_module_genes
  ) %>%
  inner_join(
    final_gene_table,
    by =
      "gene"
  ) %>%
  arrange(
    module,
    PValue
  )

write.csv(
  de_audit,
  file.path(
    TAB_OUT,
    "07_ECM_HSC_independent_gene_DE_audit_v6.3.4.1.csv"
  ),
  row.names = FALSE
)

de_module_summary <- de_audit %>%
  group_by(
    module
  ) %>%
  summarise(
    n_genes =
      n(),
    n_logFC_negative =
      sum(
        logFC <
          0,
        na.rm = TRUE
      ),
    n_logFC_positive =
      sum(
        logFC >
          0,
        na.rm = TRUE
      ),
    fraction_negative =
      mean(
        logFC <
          0,
        na.rm = TRUE
      ),
    mean_logFC =
      mean(
        logFC,
        na.rm = TRUE
      ),
    median_logFC =
      stats::median(
        logFC,
        na.rm = TRUE
      ),
    n_P_lt_0_05 =
      sum(
        PValue <
          0.05,
        na.rm = TRUE
      ),
    n_FDR_lt_0_10 =
      sum(
        FDR <
          0.10,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  de_module_summary,
  file.path(
    TAB_OUT,
    "08_ECM_HSC_independent_module_DE_direction_summary_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Integrate v6.3.3.1 sender signal with independent modules
# ==============================================================================

required_weighted <- c(
  "sample",
  "condition",
  "total_Mphi5_weighted_Pdgfb",
  "RepairResolution_weighted_Pdgfb",
  "RepairResolution_mean_Pdgfb_CP10k"
)

missing_weighted <- setdiff(
  required_weighted,
  colnames(
    weighted_axis
  )
)

if (
  length(
    missing_weighted
  )
) {
  stop(
    "v6.3.3.1 weighted table missing columns: ",
    paste(
      missing_weighted,
      collapse = ", "
    )
  )
}

ecm_sender_module <- ecm_scores %>%
  select(
    sample,
    condition,
    module,
    module_score_z
  ) %>%
  left_join(
    weighted_axis %>%
      select(
        all_of(
          required_weighted
        )
      ),
    by = c(
      "sample",
      "condition"
    )
  )

write.csv(
  ecm_sender_module,
  file.path(
    TAB_OUT,
    "09_PDGFB_sender_vs_independent_ECM_HSC_modules_v6.3.4.1.csv"
  ),
  row.names = FALSE
)

correlation_summary <- ecm_sender_module %>%
  group_by(
    module
  ) %>%
  summarise(
    rho_total_Mphi5_weighted =
      safe_spearman(
        total_Mphi5_weighted_Pdgfb,
        module_score_z
      ),
    rho_RepairResolution_weighted =
      safe_spearman(
        RepairResolution_weighted_Pdgfb,
        module_score_z
      ),
    rho_RepairResolution_percell =
      safe_spearman(
        RepairResolution_mean_Pdgfb_CP10k,
        module_score_z
      ),
    n_samples =
      n_distinct(
        sample
      ),
    note =
      "Descriptive only; n=4",
    .groups = "drop"
  )

write.csv(
  correlation_summary,
  file.path(
    TAB_OUT,
    "10_descriptive_sender_independent_module_correlations_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Figure 1:
#     Independent module scores in ECM-activated HSC
# ==============================================================================

ecm_plot <- ecm_scores %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p1 <- ggplot(
  ecm_plot,
  aes(
    x =
      sample,
    y =
      module_score_z,
    fill =
      condition
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  facet_wrap(
    ~ module,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Independent a-priori programs in ECM-activated HSC",
    subtitle =
      "NicheNet-derived PDGFB targets excluded; biological samples shown explicitly",
    x = NULL,
    y =
      "Mean gene-wise z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_ECM_HSC_independent_PDGFB_response_modules_v6.3.4.1.pdf"
  ),
  11,
  8
)


# ==============================================================================
# 13. Figure 2:
#     All HSC states module heatmap
# ==============================================================================

module_heat <- module_scores %>%
  mutate(
    HSC_sample =
      paste0(
        HSC_state,
        " | ",
        sample
      ),
    HSC_sample =
      factor(
        HSC_sample,
        levels =
          unlist(
            lapply(
              HSC3,
              function(h) {
                paste0(
                  h,
                  " | ",
                  SAMPLES
                )
              }
            )
          )
      )
  )

p2 <- ggplot(
  module_heat,
  aes(
    x =
      HSC_sample,
    y =
      module,
    fill =
      module_score_z
  )
) +
  geom_tile(
    linewidth = 0.25
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "Independent PDGF-response / HSC-state programs",
    subtitle =
      "Sample-level pseudobulk module scores across qHSC, ECM-activated HSC and Contractile HSC",
    x = NULL,
    y = NULL,
    fill =
      "Module z-score"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1
      )
  )

save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_independent_modules_all_HSC_states_heatmap_v6.3.4.1.pdf"
  ),
  15,
  5
)


# ==============================================================================
# 14. Figure 3:
#     Gene-level ECM-HSC heatmap for independent genes
# ==============================================================================

gene_heat <- pb_values %>%
  filter(
    HSC_state ==
      FOCAL_RECEIVER,
    gene %in%
      all_module_genes
  ) %>%
  inner_join(
    final_gene_table,
    by =
      "gene"
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    gene =
      factor(
        gene,
        levels =
          rev(
            unique(
              final_gene_table$gene
            )
          )
      )
  )

p3 <- ggplot(
  gene_heat,
  aes(
    x =
      sample,
    y =
      gene,
    fill =
      gene_z
  )
) +
  geom_tile(
    linewidth = 0.18
  ) +
  facet_grid(
    module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "ECM-activated HSC | independent gene-level programs",
    subtitle =
      "NicheNet-derived target genes excluded before scoring",
    x = NULL,
    y = NULL,
    fill =
      "Gene z-score"
  ) +
  theme_classic(
    base_size = 7
  )

save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_ECM_HSC_independent_gene_heatmap_v6.3.4.1.pdf"
  ),
  8,
  13
)


# ==============================================================================
# 15. Figure 4:
#     Module effect summary
# ==============================================================================

p4 <- ggplot(
  ecm_comparison,
  aes(
    x =
      reorder(
        module,
        Tx_minus_Sham
      ),
    y =
      Tx_minus_Sham
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  coord_flip() +
  labs(
    title =
      "Independent ECM-HSC program change | Tx vs Sham",
    subtitle =
      "Negative = lower in Tx; based on biological-sample pseudobulk module scores",
    x = NULL,
    y =
      "Tx mean - Sham mean"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p4,
  file.path(
    FIG_OUT,
    "04_ECM_HSC_independent_module_Tx_minus_Sham_v6.3.4.1.pdf"
  ),
  8,
  5
)


# ==============================================================================
# 16. Figure 5:
#     Sender-receiver descriptive coupling for independent modules
# ==============================================================================

coupling_plot <- ecm_sender_module %>%
  mutate(
    label =
      sample
  )

p5a <- ggplot(
  coupling_plot,
  aes(
    x =
      RepairResolution_mean_Pdgfb_CP10k,
    y =
      module_score_z,
    fill =
      condition
  )
) +
  geom_point(
    shape = 21,
    size = 3.5
  ) +
  geom_text(
    aes(
      label =
        label
    ),
    nudge_y = 0.06,
    size = 2.6
  ) +
  facet_wrap(
    ~ module,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Repair/Resolution-Mphi per-cell Pdgfb vs independent ECM-HSC programs",
    subtitle =
      "Four biological samples; descriptive only",
    x =
      "Repair/Resolution-Mphi mean Pdgfb CP10k",
    y =
      "Independent module z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p5a,
  file.path(
    FIG_OUT,
    "05_RepairResolution_percell_Pdgfb_vs_independent_HSC_modules_v6.3.4.1.pdf"
  ),
  11,
  8
)


# ==============================================================================
# 17. Figure 6:
#     Integrated independent-validation panel
# ==============================================================================

p6 <- (
  p1 /
    p4
) +
  plot_annotation(
    title =
      "Independent validation of the macrophage PDGFB -> HSC hypothesis | v6.3.4.1"
  )

save_pdf(
  p6,
  file.path(
    FIG_OUT,
    "06_PDGFB_independent_validation_summary_v6.3.4.1.pdf"
  ),
  11,
  12
)


# ==============================================================================
# 18. Final interpretation table
# ==============================================================================

interpretation <- ecm_comparison %>%
  left_join(
    de_module_summary,
    by =
      "module"
  ) %>%
  left_join(
    correlation_summary,
    by =
      "module"
  ) %>%
  mutate(
    independent_support_class = case_when(

      both_Tx_below_both_Sham &
        pairwise_direction_consistency ==
          1 &
        mean_logFC <
          0 ~
        "Strong_independent_Tx_down",

      Tx_minus_Sham <
        0 &
        pairwise_direction_consistency >=
          0.75 &
        mean_logFC <
          0 ~
        "Moderate_independent_Tx_down",

      Tx_minus_Sham >
        0 &
        pairwise_direction_consistency >=
          0.75 &
        mean_logFC >
          0 ~
        "Independent_Tx_up",

      TRUE ~
        "Mixed_or_weak"
    )
  ) %>%
  arrange(
    Tx_minus_Sham
  )

write.csv(
  interpretation,
  file.path(
    TAB_OUT,
    "11_independent_PDGFB_hypothesis_interpretation_v6.3.4.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 19. Save compact RDS
# ==============================================================================

results <- list(
  gene_set_audit =
    gene_audit,
  final_gene_sets =
    GENESETS,
  pseudobulk_gene_values =
    pb_values,
  module_scores =
    module_scores,
  ECM_module_comparison =
    ecm_comparison,
  independent_gene_DE =
    de_audit,
  independent_module_DE =
    de_module_summary,
  sender_module_table =
    ecm_sender_module,
  correlations =
    correlation_summary,
  interpretation =
    interpretation
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_PDGFB_independent_validation_results_v6.3.4.1.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 20. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "receiver_DE_input",
    "NicheNet_target_exclusion_file",
    "weighted_sender_input",
    "primary_receiver",
    "modules",
    "NicheNet_targets_used_as_positive_genes",
    "score_definition",
    "biological_replicates",
    "formal_inference_note"
  ),
  value = c(
    "v6.3.4.1",
    INPUT_RDS,
    RECEIVER_DE_FILE,
    NICHE_TARGET_FILE,
    WEIGHTED_FILE,
    FOCAL_RECEIVER,
    paste(
      names(
        GENESETS
      ),
      collapse = " | "
    ),
    "FALSE",
    "Mean gene-wise z-score from sample-level pseudobulk log2(CPM+1)",
    "Sham n=2; Tx n=2",
    "Exploratory; module gene sets were defined independently of NicheNet target selection"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.4.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.4.1.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
