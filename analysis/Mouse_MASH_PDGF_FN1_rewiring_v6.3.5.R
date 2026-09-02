#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6350)

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
# PDGF mitogenic arm vs FN1 remodeling arm
#
# Version: v6.3.5
#
# PURPOSE
#   Test the rewiring model suggested by v6.3.4.1:
#
#     PDGF / mitogenic arm       -> lower in Tx
#     FN1 / adhesion-remodeling  -> higher or preserved in Tx
#
# PRIMARY RECEIVER
#   ECM-activated HSC
#
# BIOLOGICAL SAMPLES
#   Sham1, Sham20, Tx17, Tx5
#
# INPUT
#   v6.2.0 interaction-ready Seurat object
#
# DESIGN
#
#   A. PDGF-mitogenic arm in ECM-activated HSC
#
#      1) Sample-level pseudobulk module score
#         genes:
#           Mki67, Top2a, Pcna,
#           Mcm2-7,
#           E2f1, E2f2,
#           Ccna2, Ccnb1, Cdk1,
#           Ccnd1, Ccnd2
#
#      2) Positive-cell fraction for each gene
#
#      3) Cycling-like cell fraction
#         A cell is called cycling-like if:
#           - Mki67 or Top2a is detected
#             OR
#           - at least 2 core proliferation genes are detected
#
#   B. FN1 sender arm in macrophages
#
#      1) Fn1 expression by each of the 5 Mphi subtypes
#      2) Population-weighted Fn1 output:
#
#           subtype fraction x mean per-cell Fn1 CP10k
#
#   C. FN1 receptor availability in HSC
#
#      Sdc4
#      Itga5 + Itgb1
#      Itga3 + Itgb1
#      Itga8 + Itgb1
#
#   D. Independent FN1-associated adhesion/remodeling module in HSC
#
#      Focal adhesion / matrix remodeling genes selected independently
#      of NicheNet target selection:
#
#        Ptk2, Src, Paxillin/Pxn, Vcl, Tln1, Kindlin2/Fermt2,
#        Cdc42, Rac1, Rhoa,
#        Mmp2, Mmp14, Timp1,
#        Vcam1, Icam1,
#        Lox, Loxl2
#
# IMPORTANT
#   - No CellChat rerun.
#   - No NicheNet rerun.
#   - No reclustering.
#   - Biological n = 2 Sham vs n = 2 Tx.
#   - Sample-level outputs are primary.
#   - Cell-level fractions are descriptive and summarized per biological sample.
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

module_from_pb <- function(
  pb,
  genes,
  samples
) {

  genes <- intersect(
    genes,
    rownames(
      pb
    )
  )

  if (
    !length(
      genes
    )
  ) {
    stop(
      "No module genes present."
    )
  }

  lib <- Matrix::colSums(
    pb
  )

  cpm <- sapply(
    samples,
    function(smp) {
      1e6 *
        as.numeric(
          pb[
            genes,
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
          genes
        ),
      dimnames = list(
        genes,
        samples
      )
    )
  } else {
    rownames(
      cpm
    ) <- genes
    colnames(
      cpm
    ) <- samples
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

  rownames(
    zmat
  ) <- genes

  colnames(
    zmat
  ) <- samples

  module <- tibble(
    sample =
      samples,
    condition =
      canonical_condition(
        samples
      ),
    module_score_z =
      colMeans(
        zmat,
        na.rm = TRUE
      ),
    mean_log2_CPM1 =
      colMeans(
        logcpm,
        na.rm = TRUE
      ),
    n_genes =
      length(
        genes
      )
  )

  gene_long <- as.data.frame(
    logcpm,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    pivot_longer(
      cols =
        all_of(
          samples
        ),
      names_to =
        "sample",
      values_to =
        "log2_CPM1"
    ) %>%
    left_join(
      as.data.frame(
        zmat,
        check.names = FALSE
      ) %>%
        rownames_to_column(
          "gene"
        ) %>%
        pivot_longer(
          cols =
            all_of(
              samples
            ),
          names_to =
            "sample",
          values_to =
            "gene_z"
        ),
      by = c(
        "gene",
        "sample"
      )
    ) %>%
    mutate(
      condition =
        canonical_condition(
          sample
        )
    )

  list(
    module =
      module,
    gene_long =
      gene_long,
    genes =
      genes
  )
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
      case_when(
        up_n >
          down_n ~
          "Tx_up",
        down_n >
          up_n ~
          "Tx_down",
        TRUE ~
          "Mixed_or_tied"
      ),
    pairwise_direction_consistency =
      max(
        up_n,
        down_n
      ) /
        4
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

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGF_mitogenic_FN1_rewiring_v6.3.5"
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

MPHI5 <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

HSC3 <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

FOCAL_HSC <-
  "ECM-activated HSC"

FOCAL_MPHI <-
  "Repair/Resolution-Mphi"

PDGF_MITOGENIC_GENES <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Mcm2",
  "Mcm3",
  "Mcm4",
  "Mcm5",
  "Mcm6",
  "Mcm7",
  "E2f1",
  "E2f2",
  "Ccna2",
  "Ccnb1",
  "Cdk1",
  "Ccnd1",
  "Ccnd2"
)

CORE_CYCLING_GENES <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Mcm2",
  "Mcm3",
  "Mcm4",
  "Mcm5",
  "Mcm6",
  "Mcm7",
  "Ccna2",
  "Ccnb1",
  "Cdk1"
)

FN1_LIGAND <-
  "Fn1"

FN1_RECEPTOR_GENES <- c(
  "Sdc4",
  "Itga5",
  "Itga3",
  "Itga8",
  "Itgb1"
)

FN1_REMODELING_GENES <- c(
  "Ptk2",
  "Src",
  "Pxn",
  "Vcl",
  "Tln1",
  "Fermt2",
  "Cdc42",
  "Rac1",
  "Rhoa",
  "Mmp2",
  "Mmp14",
  "Timp1",
  "Vcam1",
  "Icam1",
  "Lox",
  "Loxl2"
)

CONDITION_COLORS <- c(
  "Sham" = "#0072B2",
  "Tx" = "#D55E00"
)


# ==============================================================================
# 4. Load object
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input RDS missing: ",
    INPUT_RDS
  )
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

data_mat <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

meta <- obj@meta.data


# ==============================================================================
# 5. Gene availability audit
# ==============================================================================

gene_audit <- bind_rows(
  tibble(
    analysis =
      "PDGF_mitogenic",
    gene =
      PDGF_MITOGENIC_GENES
  ),
  tibble(
    analysis =
      "core_cycling",
    gene =
      CORE_CYCLING_GENES
  ),
  tibble(
    analysis =
      "FN1_receptor",
    gene =
      FN1_RECEPTOR_GENES
  ),
  tibble(
    analysis =
      "FN1_remodeling",
    gene =
      FN1_REMODELING_GENES
  ),
  tibble(
    analysis =
      "FN1_ligand",
    gene =
      FN1_LIGAND
  )
) %>%
  mutate(
    present =
      gene %in%
        rownames(
          counts
        )
  )

write.csv(
  gene_audit,
  file.path(
    TAB_OUT,
    "01_gene_availability_audit_v6.3.5.csv"
  ),
  row.names = FALSE
)

PDGF_MITOGENIC_PRESENT <- intersect(
  PDGF_MITOGENIC_GENES,
  rownames(
    counts
  )
)

CORE_CYCLING_PRESENT <- intersect(
  CORE_CYCLING_GENES,
  rownames(
    counts
  )
)

FN1_RECEPTOR_PRESENT <- intersect(
  FN1_RECEPTOR_GENES,
  rownames(
    counts
  )
)

FN1_REMODELING_PRESENT <- intersect(
  FN1_REMODELING_GENES,
  rownames(
    counts
  )
)

if (
  !FN1_LIGAND %in%
    rownames(
      counts
    )
) {
  stop(
    "Fn1 not found in RNA assay."
  )
}

if (
  length(
    PDGF_MITOGENIC_PRESENT
  ) <
    5
) {
  stop(
    "Too few PDGF mitogenic genes present."
  )
}

if (
  length(
    CORE_CYCLING_PRESENT
  ) <
    4
) {
  stop(
    "Too few core cycling genes present."
  )
}


# ==============================================================================
# 6. ECM-activated HSC pseudobulk PDGF-mitogenic module
# ==============================================================================

msg(
  "Computing ECM-activated HSC PDGF-mitogenic pseudobulk module..."
)

pb_ecm_hsc <- sample_pseudobulk(
  counts = counts,
  meta = meta,
  sample_col = SAMPLE_COL,
  group_col = GROUP_COL,
  group_value = FOCAL_HSC,
  samples = SAMPLES
)

pdgf_module_res <- module_from_pb(
  pb = pb_ecm_hsc,
  genes = PDGF_MITOGENIC_PRESENT,
  samples = SAMPLES
)

pdgf_module <- pdgf_module_res$module %>%
  mutate(
    module =
      "PDGF_MITOGENIC"
  )

pdgf_gene_values <- pdgf_module_res$gene_long %>%
  mutate(
    module =
      "PDGF_MITOGENIC"
  )

write.csv(
  pdgf_module,
  file.path(
    TAB_OUT,
    "02_ECM_HSC_PDGF_mitogenic_module_by_sample_v6.3.5.csv"
  ),
  row.names = FALSE
)

write.csv(
  pdgf_gene_values,
  file.path(
    TAB_OUT,
    "03_ECM_HSC_PDGF_mitogenic_gene_pseudobulk_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. ECM-activated HSC gene-positive fractions
# ==============================================================================

msg(
  "Computing ECM-activated HSC proliferation marker fractions..."
)

positive_fraction_list <- list()

for (
  smp in SAMPLES
) {

  cells <- rownames(
    meta
  )[
    as.character(
      meta[[
        SAMPLE_COL
      ]]
    ) ==
      smp &
      as.character(
        meta[[
          GROUP_COL
        ]]
      ) ==
        FOCAL_HSC
  ]

  x <- counts[
    PDGF_MITOGENIC_PRESENT,
    cells,
    drop = FALSE
  ]

  pcts <- Matrix::rowMeans(
    x >
      0
  )

  positive_fraction_list[[
    smp
  ]] <- tibble(
    sample =
      smp,
    condition =
      canonical_condition(
        smp
      ),
    gene =
      names(
        pcts
      ),
    positive_fraction =
      as.numeric(
        pcts
      ),
    n_cells =
      length(
        cells
      )
  )
}

positive_fraction <- bind_rows(
  positive_fraction_list
)

write.csv(
  positive_fraction,
  file.path(
    TAB_OUT,
    "04_ECM_HSC_proliferation_gene_positive_fraction_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Cycling-like cell fraction
# ==============================================================================

msg(
  "Calling cycling-like ECM-activated HSC cells..."
)

cycling_summary_list <- list()
cycling_cell_list <- list()

for (
  smp in SAMPLES
) {

  cells <- rownames(
    meta
  )[
    as.character(
      meta[[
        SAMPLE_COL
      ]]
    ) ==
      smp &
      as.character(
        meta[[
          GROUP_COL
        ]]
      ) ==
        FOCAL_HSC
  ]

  x <- counts[
    CORE_CYCLING_PRESENT,
    cells,
    drop = FALSE
  ]

  detected_n <- Matrix::colSums(
    x >
      0
  )

  mki67_detected <- if (
    "Mki67" %in%
      rownames(
        x
      )
  ) {
    as.numeric(
      x[
        "Mki67",
        ,
        drop = TRUE
      ] >
        0
    )
  } else {
    rep(
      0,
      length(
        cells
      )
    )
  }

  top2a_detected <- if (
    "Top2a" %in%
      rownames(
        x
      )
  ) {
    as.numeric(
      x[
        "Top2a",
        ,
        drop = TRUE
      ] >
        0
    )
  } else {
    rep(
      0,
      length(
        cells
      )
    )
  }

  cycling_like <- (
    mki67_detected >
      0
  ) |
    (
      top2a_detected >
        0
    ) |
    (
      detected_n >=
        2
    )

  cycling_cell_list[[
    smp
  ]] <- tibble(
    cell =
      cells,
    sample =
      smp,
    condition =
      canonical_condition(
        smp
      ),
    n_core_cycling_genes_detected =
      as.numeric(
        detected_n
      ),
    Mki67_detected =
      as.logical(
        mki67_detected
      ),
    Top2a_detected =
      as.logical(
        top2a_detected
      ),
    cycling_like =
      cycling_like
  )

  cycling_summary_list[[
    smp
  ]] <- tibble(
    sample =
      smp,
    condition =
      canonical_condition(
        smp
      ),
    n_ECM_HSC =
      length(
        cells
      ),
    n_cycling_like =
      sum(
        cycling_like
      ),
    cycling_like_fraction =
      mean(
        cycling_like
      ),
    Mki67_positive_fraction =
      mean(
        mki67_detected >
          0
      ),
    Top2a_positive_fraction =
      mean(
        top2a_detected >
          0
      ),
    mean_core_cycling_genes_detected =
      mean(
        detected_n
      )
  )
}

cycling_cells <- bind_rows(
  cycling_cell_list
)

cycling_summary <- bind_rows(
  cycling_summary_list
)

write.csv(
  cycling_summary,
  file.path(
    TAB_OUT,
    "05_ECM_HSC_cycling_like_fraction_by_sample_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Replicate-aware summary of PDGF mitogenic arm
# ==============================================================================

pdgf_wide <- pdgf_module %>%
  select(
    sample,
    module_score_z
  ) %>%
  pivot_wider(
    names_from =
      sample,
    values_from =
      module_score_z
  )

pdgf_pairwise <- pairwise_direction(
  pdgf_wide$Sham1[[1]],
  pdgf_wide$Sham20[[1]],
  pdgf_wide$Tx17[[1]],
  pdgf_wide$Tx5[[1]]
)

pdgf_summary <- bind_cols(
  tibble(
    Sham1 =
      pdgf_wide$Sham1[[1]],
    Sham20 =
      pdgf_wide$Sham20[[1]],
    Tx17 =
      pdgf_wide$Tx17[[1]],
    Tx5 =
      pdgf_wide$Tx5[[1]]
  ),
  pdgf_pairwise
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
        )
  )

write.csv(
  pdgf_summary,
  file.path(
    TAB_OUT,
    "06_PDGF_mitogenic_replicate_summary_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Fn1 macrophage expression and weighted output
# ==============================================================================

msg(
  "Computing macrophage Fn1 sender output..."
)

cell_lib <- Matrix::colSums(
  counts
)

fn1_raw <- as.numeric(
  counts[
    FN1_LIGAND,
    ,
    drop = TRUE
  ]
)

names(
  fn1_raw
) <- colnames(
  counts
)

fn1_cp10k <- 10000 *
  fn1_raw /
  pmax(
    as.numeric(
      cell_lib
    ),
    1
  )

names(
  fn1_cp10k
) <- colnames(
  counts
)

mphi_fn1_list <- list()

for (
  smp in SAMPLES
) {

  total_mphi_cells <- rownames(
    meta
  )[
    as.character(
      meta[[
        SAMPLE_COL
      ]]
    ) ==
      smp &
      as.character(
        meta[[
          GROUP_COL
        ]]
      ) %in%
        MPHI5
  ]

  total_n <- length(
    total_mphi_cells
  )

  for (
    ct in MPHI5
  ) {

    cells <- rownames(
      meta
    )[
      as.character(
        meta[[
          SAMPLE_COL
        ]]
      ) ==
        smp &
        as.character(
          meta[[
            GROUP_COL
          ]]
        ) ==
          ct
    ]

    if (
      !length(
        cells
      )
    ) {
      next
    }

    vals <- fn1_cp10k[
      cells
    ]

    raw_vals <- fn1_raw[
      cells
    ]

    frac <- length(
      cells
    ) /
      total_n

    mphi_fn1_list[[
      paste(
        smp,
        ct,
        sep = "__"
      )
    ]] <- tibble(
      sample =
        smp,
      condition =
        canonical_condition(
          smp
        ),
      Mphi_subtype =
        ct,
      n_cells =
        length(
          cells
        ),
      fraction_within_Mphi5 =
        frac,
      Fn1_pct_positive =
        mean(
          raw_vals >
            0
        ),
      mean_Fn1_CP10k =
        mean(
          vals
        ),
      weighted_Fn1_output =
        frac *
          mean(
            vals
          )
    )
  }
}

mphi_fn1 <- bind_rows(
  mphi_fn1_list
) %>%
  group_by(
    sample
  ) %>%
  mutate(
    total_Mphi5_weighted_Fn1 =
      sum(
        weighted_Fn1_output
      ),
    contribution_fraction =
      ifelse(
        total_Mphi5_weighted_Fn1 >
          0,
        weighted_Fn1_output /
          total_Mphi5_weighted_Fn1,
        0
      )
  ) %>%
  ungroup()

write.csv(
  mphi_fn1,
  file.path(
    TAB_OUT,
    "07_Mphi5_Fn1_population_weighted_output_v6.3.5.csv"
  ),
  row.names = FALSE
)

total_fn1 <- mphi_fn1 %>%
  group_by(
    sample,
    condition
  ) %>%
  summarise(
    total_Mphi5_weighted_Fn1 =
      first(
        total_Mphi5_weighted_Fn1
      ),
    .groups = "drop"
  )

write.csv(
  total_fn1,
  file.path(
    TAB_OUT,
    "08_total_Mphi5_weighted_Fn1_by_sample_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. FN1 receptor availability in HSC states
# ==============================================================================

msg(
  "Computing FN1 receptor availability in HSC..."
)

fn1_receptor_list <- list()

for (
  smp in SAMPLES
) {

  for (
    hsc in HSC3
  ) {

    cells <- rownames(
      meta
    )[
      as.character(
        meta[[
          SAMPLE_COL
        ]]
      ) ==
        smp &
        as.character(
          meta[[
            GROUP_COL
          ]]
        ) ==
          hsc
    ]

    if (
      !length(
        cells
      )
    ) {
      next
    }

    x <- counts[
      FN1_RECEPTOR_PRESENT,
      cells,
      drop = FALSE
    ]

    d <- data_mat[
      FN1_RECEPTOR_PRESENT,
      cells,
      drop = FALSE
    ]

    fn1_receptor_list[[
      paste(
        smp,
        hsc,
        sep = "__"
      )
    ]] <- tibble(
      sample =
        smp,
      condition =
        canonical_condition(
          smp
        ),
      HSC_state =
        hsc,
      receptor_gene =
        FN1_RECEPTOR_PRESENT,
      positive_fraction =
        as.numeric(
          Matrix::rowMeans(
            x >
              0
          )
        ),
      mean_logexpr =
        as.numeric(
          Matrix::rowMeans(
            d
          )
        ),
      n_cells =
        length(
          cells
        )
    )
  }
}

fn1_receptors <- bind_rows(
  fn1_receptor_list
)

write.csv(
  fn1_receptors,
  file.path(
    TAB_OUT,
    "09_FN1_receptor_availability_HSC3_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Independent FN1 remodeling module in ECM-HSC
# ==============================================================================

msg(
  "Computing independent FN1 remodeling module in ECM-activated HSC..."
)

fn1_module_res <- module_from_pb(
  pb = pb_ecm_hsc,
  genes = FN1_REMODELING_PRESENT,
  samples = SAMPLES
)

fn1_module <- fn1_module_res$module %>%
  mutate(
    module =
      "FN1_REMODELING"
  )

fn1_gene_values <- fn1_module_res$gene_long %>%
  mutate(
    module =
      "FN1_REMODELING"
  )

write.csv(
  fn1_module,
  file.path(
    TAB_OUT,
    "10_ECM_HSC_FN1_remodeling_module_by_sample_v6.3.5.csv"
  ),
  row.names = FALSE
)

write.csv(
  fn1_gene_values,
  file.path(
    TAB_OUT,
    "11_ECM_HSC_FN1_remodeling_gene_pseudobulk_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Replicate-aware summary of FN1 remodeling arm
# ==============================================================================

fn1_wide <- fn1_module %>%
  select(
    sample,
    module_score_z
  ) %>%
  pivot_wider(
    names_from =
      sample,
    values_from =
      module_score_z
  )

fn1_pairwise <- pairwise_direction(
  fn1_wide$Sham1[[1]],
  fn1_wide$Sham20[[1]],
  fn1_wide$Tx17[[1]],
  fn1_wide$Tx5[[1]]
)

fn1_summary <- bind_cols(
  tibble(
    Sham1 =
      fn1_wide$Sham1[[1]],
    Sham20 =
      fn1_wide$Sham20[[1]],
    Tx17 =
      fn1_wide$Tx17[[1]],
    Tx5 =
      fn1_wide$Tx5[[1]]
  ),
  fn1_pairwise
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
    both_Tx_above_both_Sham =
      min(
        Tx17,
        Tx5
      ) >
        max(
          Sham1,
          Sham20
        )
  )

write.csv(
  fn1_summary,
  file.path(
    TAB_OUT,
    "12_FN1_remodeling_replicate_summary_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Integrated sample-level rewiring table
# ==============================================================================

integrated <- pdgf_module %>%
  select(
    sample,
    condition,
    PDGF_mitogenic_z =
      module_score_z
  ) %>%
  left_join(
    cycling_summary %>%
      select(
        sample,
        cycling_like_fraction,
        Mki67_positive_fraction,
        Top2a_positive_fraction
      ),
    by =
      "sample"
  ) %>%
  left_join(
    total_fn1,
    by = c(
      "sample",
      "condition"
    )
  ) %>%
  left_join(
    fn1_module %>%
      select(
        sample,
        FN1_remodeling_z =
          module_score_z
      ),
    by =
      "sample"
  )

write.csv(
  integrated,
  file.path(
    TAB_OUT,
    "13_integrated_PDGF_FN1_rewiring_by_sample_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Figure 1: PDGF mitogenic arm
# ==============================================================================

p1a <- ggplot(
  pdgf_module %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      module_score_z,
    fill =
      condition
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "ECM-HSC PDGF mitogenic module",
    x = NULL,
    y =
      "Mean gene-wise z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p1b <- ggplot(
  cycling_summary %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      100 *
        cycling_like_fraction,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Cycling-like ECM-HSC fraction",
    x = NULL,
    y =
      "% cycling-like cells",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p1 <- p1a +
  p1b +
  plot_annotation(
    title =
      "PDGF / mitogenic arm in ECM-activated HSC"
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_PDGF_mitogenic_arm_ECM_HSC_v6.3.5.pdf"
  ),
  12,
  5
)


# ==============================================================================
# 16. Figure 2: proliferation gene heatmap
# ==============================================================================

pdgf_gene_heat <- pdgf_gene_values %>%
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
            PDGF_MITOGENIC_PRESENT
          )
      )
  )

p2 <- ggplot(
  pdgf_gene_heat,
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
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "ECM-activated HSC proliferation / cell-cycle genes",
    subtitle =
      "Sample-level pseudobulk; gene-wise z-score",
    x = NULL,
    y = NULL,
    fill =
      "Gene z-score"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_ECM_HSC_proliferation_gene_heatmap_v6.3.5.pdf"
  ),
  7,
  8
)


# ==============================================================================
# 17. Figure 3: macrophage Fn1 sender arm
# ==============================================================================

p3a <- ggplot(
  mphi_fn1 %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        ),
      Mphi_subtype =
        factor(
          Mphi_subtype,
          levels =
            MPHI5
        )
    ),
  aes(
    x =
      sample,
    y =
      weighted_Fn1_output,
    fill =
      Mphi_subtype
  )
) +
  geom_col() +
  labs(
    title =
      "Population-weighted macrophage Fn1 output",
    x = NULL,
    y =
      "Weighted Fn1 CP10k",
    fill =
      "Mphi subtype"
  ) +
  theme_classic(
    base_size = 9
  )

p3b <- ggplot(
  total_fn1 %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      total_Mphi5_weighted_Fn1,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Total Mphi5 weighted Fn1",
    x = NULL,
    y =
      "Weighted Fn1 CP10k",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p3 <- p3a +
  p3b +
  plot_annotation(
    title =
      "FN1 sender arm in macrophages"
  )

save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_FN1_sender_arm_Mphi5_v6.3.5.pdf"
  ),
  13,
  5
)


# ==============================================================================
# 18. Figure 4: FN1 receptor availability
# ==============================================================================

p4 <- ggplot(
  fn1_receptors %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        ),
      HSC_state =
        factor(
          HSC_state,
          levels =
            HSC3
        )
    ),
  aes(
    x =
      sample,
    y =
      HSC_state,
    size =
      positive_fraction,
    fill =
      mean_logexpr
  )
) +
  geom_point(
    shape = 21
  ) +
  facet_wrap(
    ~ receptor_gene,
    ncol = 3
  ) +
  scale_fill_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) +
  scale_size(
    range = c(
      1,
      9
    )
  ) +
  labs(
    title =
      "FN1 receptor availability across HSC states",
    subtitle =
      "Sdc4 and integrin components",
    x = NULL,
    y = NULL,
    size =
      "Fraction expressing",
    fill =
      "Mean expression"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p4,
  file.path(
    FIG_OUT,
    "04_FN1_receptor_availability_HSC3_v6.3.5.pdf"
  ),
  13,
  7
)


# ==============================================================================
# 19. Figure 5: FN1 remodeling arm
# ==============================================================================

p5a <- ggplot(
  fn1_module %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      module_score_z,
    fill =
      condition
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "ECM-HSC FN1 remodeling module",
    x = NULL,
    y =
      "Mean gene-wise z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

fn1_gene_heat <- fn1_gene_values %>%
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
            FN1_REMODELING_PRESENT
          )
      )
  )

p5b <- ggplot(
  fn1_gene_heat,
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
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "FN1-associated remodeling genes",
    x = NULL,
    y = NULL,
    fill =
      "Gene z-score"
  ) +
  theme_classic(
    base_size = 8
  )

p5 <- p5a /
  p5b +
  plot_annotation(
    title =
      "FN1 / adhesion-remodeling arm in ECM-activated HSC"
  )

save_pdf(
  p5,
  file.path(
    FIG_OUT,
    "05_FN1_remodeling_arm_ECM_HSC_v6.3.5.pdf"
  ),
  9,
  11
)


# ==============================================================================
# 20. Figure 6: integrated rewiring panel
# ==============================================================================

p6a <- ggplot(
  integrated %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      PDGF_mitogenic_z,
    fill =
      condition
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "PDGF mitogenic arm",
    x = NULL,
    y =
      "Module z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

p6b <- ggplot(
  integrated %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      100 *
        cycling_like_fraction,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Cycling-like ECM-HSC",
    x = NULL,
    y =
      "%",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

p6c <- ggplot(
  integrated %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      total_Mphi5_weighted_Fn1,
    fill =
      condition
  )
) +
  geom_col() +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Mphi5 weighted Fn1",
    x = NULL,
    y =
      "Weighted CP10k",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

p6d <- ggplot(
  integrated %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      FN1_remodeling_z,
    fill =
      condition
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "FN1 remodeling arm",
    x = NULL,
    y =
      "Module z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

p6 <- (
  p6a +
    p6b
) /
  (
    p6c +
      p6d
  ) +
  plot_annotation(
    title =
      "PDGF mitogenic arm vs FN1 remodeling arm | v6.3.5"
  )

save_pdf(
  p6,
  file.path(
    FIG_OUT,
    "06_PDGF_FN1_rewiring_integrated_panel_v6.3.5.pdf"
  ),
  12,
  9
)


# ==============================================================================
# 21. Final interpretation summary
# ==============================================================================

condition_summary <- integrated %>%
  group_by(
    condition
  ) %>%
  summarise(
    mean_PDGF_mitogenic_z =
      mean(
        PDGF_mitogenic_z
      ),
    mean_cycling_like_fraction =
      mean(
        cycling_like_fraction
      ),
    mean_total_Mphi5_weighted_Fn1 =
      mean(
        total_Mphi5_weighted_Fn1
      ),
    mean_FN1_remodeling_z =
      mean(
        FN1_remodeling_z
      ),
    .groups = "drop"
  )

sham_vals <- condition_summary %>%
  filter(
    condition ==
      "Sham"
  )

tx_vals <- condition_summary %>%
  filter(
    condition ==
      "Tx"
  )

interpretation <- tibble(
  metric = c(
    "PDGF_mitogenic_module_Tx_minus_Sham",
    "cycling_like_fraction_Tx_minus_Sham",
    "total_Mphi5_weighted_Fn1_Tx_minus_Sham",
    "FN1_remodeling_module_Tx_minus_Sham",
    "PDGF_pairwise_direction_consistency",
    "FN1_pairwise_direction_consistency",
    "PDGF_both_Tx_below_both_Sham",
    "FN1_both_Tx_above_both_Sham"
  ),
  value = c(
    tx_vals$mean_PDGF_mitogenic_z -
      sham_vals$mean_PDGF_mitogenic_z,
    tx_vals$mean_cycling_like_fraction -
      sham_vals$mean_cycling_like_fraction,
    tx_vals$mean_total_Mphi5_weighted_Fn1 -
      sham_vals$mean_total_Mphi5_weighted_Fn1,
    tx_vals$mean_FN1_remodeling_z -
      sham_vals$mean_FN1_remodeling_z,
    pdgf_summary$pairwise_direction_consistency,
    fn1_summary$pairwise_direction_consistency,
    as.numeric(
      pdgf_summary$both_Tx_below_both_Sham
    ),
    as.numeric(
      fn1_summary$both_Tx_above_both_Sham
    )
  )
)

write.csv(
  interpretation,
  file.path(
    TAB_OUT,
    "14_PDGF_FN1_rewiring_interpretation_summary_v6.3.5.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 22. Save compact RDS
# ==============================================================================

results <- list(
  gene_audit =
    gene_audit,
  PDGF_module =
    pdgf_module,
  PDGF_gene_values =
    pdgf_gene_values,
  proliferation_positive_fraction =
    positive_fraction,
  cycling_summary =
    cycling_summary,
  PDGF_replicate_summary =
    pdgf_summary,
  Mphi5_Fn1 =
    mphi_fn1,
  total_Fn1 =
    total_fn1,
  FN1_receptors =
    fn1_receptors,
  FN1_module =
    fn1_module,
  FN1_gene_values =
    fn1_gene_values,
  FN1_replicate_summary =
    fn1_summary,
  integrated =
    integrated,
  interpretation =
    interpretation
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_PDGF_FN1_rewiring_results_v6.3.5.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 23. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "primary_HSC",
    "PDGF_mitogenic_genes",
    "cycling_like_definition",
    "FN1_ligand",
    "FN1_receptor_genes",
    "FN1_remodeling_genes",
    "biological_replicates",
    "CellChat_recomputed",
    "NicheNet_recomputed",
    "formal_inference_note"
  ),
  value = c(
    "v6.3.5",
    INPUT_RDS,
    FOCAL_HSC,
    paste(
      PDGF_MITOGENIC_PRESENT,
      collapse = ","
    ),
    "Mki67+ OR Top2a+ OR >=2 core proliferation genes detected",
    FN1_LIGAND,
    paste(
      FN1_RECEPTOR_PRESENT,
      collapse = ","
    ),
    paste(
      FN1_REMODELING_PRESENT,
      collapse = ","
    ),
    "Sham n=2; Tx n=2",
    "FALSE",
    "FALSE",
    "Exploratory; primary inference based on biological-sample patterns"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.5.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.5.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
