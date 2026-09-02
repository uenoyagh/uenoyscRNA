#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6530)

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

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "Package 'edgeR' is required for v6.5.3.1. ",
    "Install it inside the project renv before running."
  )
}

# ==============================================================================
# Mouse MASH scRNA-seq
# Ballooning-focused hepatocyte analysis
#
# Version: v6.5.3.1
#
# FIX FROM v6.5.3:
#   Avoid tidy-evaluation/name-shadowing collisions between loop variables
#   and tibble/dplyr columns named scope, program, and gene.
#
# INPUT:
#   Mouse_MASH_Hepatocyte_FINAL_annotated_v6.5.2.rds
#
# PURPOSE:
#   Test whether scRNA-seq contains molecular/QC findings that are concordant
#   with the histological observation that hepatocyte ballooning is lower
#   after transplantation.
#
# IMPORTANT:
#   This script DOES NOT claim to diagnose ballooning from scRNA-seq.
#   Histological ballooning is a morphological phenotype. The present analysis
#   evaluates orthogonal molecular and QC proxies that may be concordant with
#   lower hepatocyte injury/ballooning burden.
#
# PRIMARY EVIDENCE LAYERS:
#   1) MT-high/QC hepatocyte burden
#   2) Cytoskeletal / proteostasis response
#   3) ER / lipotoxic stress response
#   4) Oxidative / mitochondrial stress response
#   5) Damage / cell-cycle-arrest response
#   6) Acute-phase response
#   7) Direct candidate genes including Shh, Krt8/Krt18, Sqstm1
#
# ANALYSIS SCOPES:
#   - All_primary_Hepatocytes:
#       all final hepatocyte states except MT_high_QC_Hepatocyte
#   - Injury_inflammatory_Hepatocyte:
#       the v6.5.2 injury/inflammatory state
#   - MT_high_QC_Hepatocyte:
#       analyzed separately as a QC/injury-associated state
#
# STATISTICAL NOTE:
#   n = 2 Sham vs n = 2 Tx.
#   edgeR pseudobulk results are exploratory.
#   Replicate direction consistency and effect size are emphasized.
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
    "RNA" %in%
      Assays(object) &&
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

pairwise_direction_summary <- function(
  sham_values,
  tx_values,
  tolerance = 1e-12
) {

  diffs <- as.vector(
    outer(
      tx_values,
      sham_values,
      "-"
    )
  )

  up_n <- sum(
    diffs > tolerance,
    na.rm = TRUE
  )

  down_n <- sum(
    diffs < -tolerance,
    na.rm = TRUE
  )

  tie_n <- sum(
    abs(diffs) <= tolerance,
    na.rm = TRUE
  )

  n_comp <- sum(
    is.finite(diffs)
  )

  direction <- if (
    up_n > down_n
  ) {
    "Tx_up"
  } else if (
    down_n > up_n
  ) {
    "Tx_down"
  } else {
    "Mixed_or_tied"
  }

  consistency <- if (
    n_comp > 0
  ) {
    max(
      up_n,
      down_n
    ) / n_comp
  } else {
    NA_real_
  }

  both_tx_above <- if (
    all(
      is.finite(
        c(
          sham_values,
          tx_values
        )
      )
    )
  ) {
    min(tx_values) >
      max(sham_values)
  } else {
    NA
  }

  both_tx_below <- if (
    all(
      is.finite(
        c(
          sham_values,
          tx_values
        )
      )
    )
  ) {
    max(tx_values) <
      min(sham_values)
  } else {
    NA
  }

  grade <- if (
    isTRUE(
      both_tx_above
    )
  ) {
    "Strong_Tx_up"
  } else if (
    isTRUE(
      both_tx_below
    )
  ) {
    "Strong_Tx_down"
  } else if (
    direction ==
      "Tx_up" &&
    consistency >=
      0.75
  ) {
    "Moderate_Tx_up"
  } else if (
    direction ==
      "Tx_down" &&
    consistency >=
      0.75
  ) {
    "Moderate_Tx_down"
  } else {
    "Mixed_or_weak"
  }

  tibble(
    pairwise_up_n =
      up_n,
    pairwise_down_n =
      down_n,
    pairwise_tie_n =
      tie_n,
    pairwise_direction =
      direction,
    pairwise_direction_consistency =
      consistency,
    both_Tx_above_both_Sham =
      both_tx_above,
    both_Tx_below_both_Sham =
      both_tx_below,
    evidence_grade =
      grade
  )
}

four_sample_summary <- function(
  df,
  value_col,
  metric_name = value_col
) {

  dat <- df %>%
    transmute(
      sample =
        as.character(
          sample
        ),
      condition =
        as.character(
          condition
        ),
      value =
        .data[[
          value_col
        ]]
    )

  sample_values <- setNames(
    dat$value,
    dat$sample
  )

  required_samples <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  missing_samples <- setdiff(
    required_samples,
    names(
      sample_values
    )
  )

  if (
    length(
      missing_samples
    )
  ) {
    stop(
      "Missing sample(s) in four_sample_summary: ",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  }

  sham <- c(
    sample_values[
      "Sham1"
    ],
    sample_values[
      "Sham20"
    ]
  )

  tx <- c(
    sample_values[
      "Tx17"
    ],
    sample_values[
      "Tx5"
    ]
  )

  direction <- pairwise_direction_summary(
    sham,
    tx
  )

  bind_cols(
    tibble(
      metric =
        metric_name,
      Sham1 =
        unname(
          sample_values[
            "Sham1"
          ]
        ),
      Sham20 =
        unname(
          sample_values[
            "Sham20"
          ]
        ),
      Tx17 =
        unname(
          sample_values[
            "Tx17"
          ]
        ),
      Tx5 =
        unname(
          sample_values[
            "Tx5"
          ]
        ),
      Sham_mean =
        mean(
          sham,
          na.rm = TRUE
        ),
      Tx_mean =
        mean(
          tx,
          na.rm = TRUE
        ),
      Tx_minus_Sham =
        mean(
          tx,
          na.rm = TRUE
        ) -
          mean(
            sham,
            na.rm = TRUE
          )
    ),
    direction
  )
}

aggregate_counts_by_group <- function(
  counts,
  groups
) {

  groups <- as.character(
    groups
  )

  group_levels <- unique(
    groups
  )

  result <- lapply(
    group_levels,
    function(g) {

      idx <- which(
        groups ==
          g
      )

      if (
        length(
          idx
        ) ==
          1
      ) {
        as.numeric(
          counts[
            ,
            idx
          ]
        )
      } else {
        Matrix::rowSums(
          counts[
            ,
            idx,
            drop = FALSE
          ]
        )
      }
    }
  )

  mat <- do.call(
    cbind,
    result
  )

  rownames(
    mat
  ) <- rownames(
    counts
  )

  colnames(
    mat
  ) <- group_levels

  mat
}

zscore_rows <- function(mat) {

  t(
    apply(
      mat,
      1,
      function(x) {

        s <- stats::sd(
          x,
          na.rm = TRUE
        )

        if (
          !is.finite(
            s
          ) ||
          s ==
            0
        ) {
          rep(
            0,
            length(
              x
            )
          )
        } else {
          (
            x -
              mean(
                x,
                na.rm = TRUE
              )
          ) /
            s
        }
      }
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
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_FINAL_state_pseudobulk_v6.5.2",
  "RDS",
  "Mouse_MASH_Hepatocyte_FINAL_annotated_v6.5.2.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_ballooning_focused_v6.5.3.1"
)

RDS_OUT <- file.path(
  OUT,
  "RDS"
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

for (
  d in c(
    OUT,
    RDS_OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Metadata and scope settings
# ==============================================================================

STATE_COL <-
  "hepatocyte_state_FINAL_v652"

SAMPLE_COL <-
  "sample_hep_v650"

CONDITION_COL <-
  "condition_hep_v650"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

MT_STATE <-
  "MT_high_QC_Hepatocyte"

INJURY_STATE <-
  "Injury_inflammatory_Hepatocyte"

SCOPES <- c(
  "All_primary_Hepatocytes",
  "Injury_inflammatory_Hepatocyte",
  "MT_high_QC_Hepatocyte"
)

MIN_CELLS_PB <-
  10


# ==============================================================================
# 4. Ballooning-related proxy programs
# ==============================================================================

BALLOONING_PROGRAMS <- list(

  Cytoskeletal_keratin_response = c(
    "Krt8",
    "Krt18",
    "Krt19",
    "Krt20",
    "Vim"
  ),

  Proteostasis_aggresome_response = c(
    "Sqstm1",
    "Ubb",
    "Ubc",
    "Hspa1a",
    "Hspa1b",
    "Hspa5",
    "Hsp90aa1",
    "Dnajb1",
    "Bag3",
    "Vcp"
  ),

  ER_lipotoxic_stress = c(
    "Ddit3",
    "Atf3",
    "Atf4",
    "Xbp1",
    "Hspa5",
    "Herpud1",
    "Trib3",
    "Nupr1"
  ),

  Oxidative_mito_stress = c(
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Txnrd1",
    "Srxn1",
    "Sod2",
    "Prdx1",
    "Prdx3"
  ),

  Damage_cell_cycle_arrest = c(
    "Cdkn1a",
    "Gadd45a",
    "Gadd45b",
    "Bbc3",
    "Pmaip1",
    "Bax",
    "Trp53inp1"
  ),

  Acute_phase_response = c(
    "Saa1",
    "Saa2",
    "Orm1",
    "Orm2",
    "Lcn2",
    "Serpina3n"
  ),

  Lipid_storage_lipotoxic_context = c(
    "Plin2",
    "Cd36",
    "Cidec",
    "Mogat1",
    "Dgat2",
    "Scd1",
    "Fasn",
    "Srebf1",
    "Acaca"
  )
)

DIRECT_CANDIDATE_GENES <- unique(
  c(
    # Ballooning-associated / orthogonal candidates
    "Shh",
    "Ihh",

    # Cytoskeleton / proteostasis
    "Krt8",
    "Krt18",
    "Sqstm1",
    "Ubb",
    "Ubc",
    "Dnajb1",
    "Bag3",

    # ER / lipotoxic injury
    "Ddit3",
    "Atf3",
    "Atf4",
    "Xbp1",
    "Hspa5",
    "Trib3",
    "Nupr1",

    # Oxidative / mitochondrial response
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",

    # Damage / arrest
    "Cdkn1a",
    "Gadd45a",
    "Gadd45b",
    "Bbc3",
    "Pmaip1",

    # Acute phase
    "Saa1",
    "Saa2",
    "Orm1",
    "Orm2",
    "Lcn2",

    # Lipotoxic context
    "Plin2",
    "Cd36",
    "Cidec",
    "Mogat1",
    "Dgat2"
  )
)


# ==============================================================================
# 5. Load final v6.5.2 object
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input v6.5.2 RDS missing: ",
    INPUT_RDS
  )
}

msg(
  "Loading v6.5.2 final Hepatocyte object..."
)

hep <- readRDS(
  INPUT_RDS
)

DefaultAssay(
  hep
) <- "RNA"

hep <- safe_join_rna(
  hep
)

required_meta <- c(
  STATE_COL,
  SAMPLE_COL,
  CONDITION_COL,
  "percent.mt",
  "nCount_RNA",
  "nFeature_RNA"
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    hep@meta.data
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

state_chr <- as.character(
  hep@meta.data[[
    STATE_COL
  ]]
)

if (
  !all(
    c(
      MT_STATE,
      INJURY_STATE
    ) %in%
      unique(
        state_chr
      )
  )
) {
  stop(
    "Required final states are missing."
  )
}


# ==============================================================================
# 6. Define analysis scopes
# ==============================================================================

scope_label <- rep(
  NA_character_,
  ncol(
    hep
  )
)

# Scope memberships are saved separately rather than as one exclusive label.
hep$balloon_scope_primary_v6531 <-
  state_chr !=
    MT_STATE

hep$balloon_scope_injury_v6531 <-
  state_chr ==
    INJURY_STATE

hep$balloon_scope_MT_v6531 <-
  state_chr ==
    MT_STATE

scope_membership <- tibble(
  cell =
    colnames(
      hep
    ),
  sample =
    as.character(
      hep@meta.data[[
        SAMPLE_COL
      ]]
    ),
  condition =
    as.character(
      hep@meta.data[[
        CONDITION_COL
      ]]
    ),
  state =
    state_chr,
  All_primary_Hepatocytes =
    hep$balloon_scope_primary_v6531,
  Injury_inflammatory_Hepatocyte =
    hep$balloon_scope_injury_v6531,
  MT_high_QC_Hepatocyte =
    hep$balloon_scope_MT_v6531
)


# ==============================================================================
# 7. Scope cell counts
# ==============================================================================

scope_count_rows <- list()

for (
  scope_name in SCOPES
) {

  tab <- scope_membership %>%
    filter(
      .data[[
        scope_name
      ]]
    ) %>%
    count(
      sample,
      condition,
      name =
        "n_cells"
    ) %>%
    complete(
      sample =
        SAMPLES,
      fill =
        list(
          n_cells =
            0
        )
    ) %>%
    mutate(
      condition =
        ifelse(
          grepl(
            "^Tx",
            sample
          ),
          "Tx",
          "Sham"
        ),
      scope =
        scope_name
    )

  scope_count_rows[[
    scope_name
  ]] <- tab
}

scope_counts <- bind_rows(
  scope_count_rows
)

write.csv(
  scope_counts,
  file.path(
    TAB_OUT,
    "01_ballooning_scope_cell_counts_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Gene audit
# ==============================================================================

gene_audit <- bind_rows(
  lapply(
    names(
      BALLOONING_PROGRAMS
    ),
    function(program_name) {

      genes_this <- BALLOONING_PROGRAMS[[
        program_name
      ]]

      tibble(
        category =
          "program",
        program =
          program_name,
        gene =
          genes_this,
        present =
          genes_this %in%
            rownames(
              hep
            )
      )
    }
  ),
  tibble(
    category =
      "direct_candidate",
    program =
      "Direct_candidate_gene",
    gene =
      DIRECT_CANDIDATE_GENES,
    present =
      DIRECT_CANDIDATE_GENES %in%
        rownames(
          hep
        )
  )
)

write.csv(
  gene_audit,
  file.path(
    TAB_OUT,
    "02_ballooning_proxy_gene_audit_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. MT-high/QC burden and QC metrics
# ==============================================================================

qc_by_sample <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    condition =
      as.character(
        .data[[
          CONDITION_COL
        ]]
      ),
    state =
      as.character(
        .data[[
          STATE_COL
        ]]
      ),
    percent.mt,
    nCount_RNA,
    nFeature_RNA
  ) %>%
  group_by(
    sample,
    condition
  ) %>%
  summarise(
    n_all_clean =
      n(),
    n_MT_high =
      sum(
        state ==
          MT_STATE
      ),
    MT_high_fraction =
      mean(
        state ==
          MT_STATE
      ),
    median_percent_mt_all_clean =
      median(
        percent.mt,
        na.rm = TRUE
      ),
    pct_all_clean_mt_gt_10 =
      mean(
        percent.mt >
          10,
        na.rm = TRUE
      ),
    pct_all_clean_mt_gt_20 =
      mean(
        percent.mt >
          20,
        na.rm = TRUE
      ),
    median_percent_mt_MT_high =
      median(
        percent.mt[
          state ==
            MT_STATE
        ],
        na.rm = TRUE
      ),
    median_nFeature_MT_high =
      median(
        nFeature_RNA[
          state ==
            MT_STATE
        ],
        na.rm = TRUE
      ),
    .groups =
      "drop"
  )

write.csv(
  qc_by_sample,
  file.path(
    TAB_OUT,
    "03_MT_high_QC_burden_by_sample_v6.5.3.1.csv"
  ),
  row.names = FALSE
)

qc_metric_names <- c(
  "MT_high_fraction",
  "median_percent_mt_all_clean",
  "pct_all_clean_mt_gt_10",
  "pct_all_clean_mt_gt_20",
  "median_percent_mt_MT_high",
  "median_nFeature_MT_high"
)

qc_summary <- bind_rows(
  lapply(
    qc_metric_names,
    function(m) {
      four_sample_summary(
        qc_by_sample,
        m,
        m
      )
    }
  )
)

write.csv(
  qc_summary,
  file.path(
    TAB_OUT,
    "04_MT_high_QC_replicate_summary_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Direct candidate-positive fractions by scope
# ==============================================================================

counts <- GetAssayData(
  hep,
  assay =
    "RNA",
  layer =
    "counts"
)

direct_present <- intersect(
  DIRECT_CANDIDATE_GENES,
  rownames(
    counts
  )
)

positive_fraction_rows <- list()

for (
  scope_name in SCOPES
) {

  scope_cells <- scope_membership$cell[
    scope_membership[[
      scope_name
    ]]
  ]

  for (
    sample_name in SAMPLES
  ) {

    cells <- intersect(
      scope_cells,
      scope_membership$cell[
        scope_membership$sample ==
          sample_name
      ]
    )

    if (
      !length(
        cells
      )
    ) {
      next
    }

    idx <- match(
      cells,
      colnames(
        counts
      )
    )

    mat <- counts[
      direct_present,
      idx,
      drop =
        FALSE
    ]

    frac <- Matrix::rowMeans(
      mat >
        0
    )

    positive_fraction_rows[[
      paste(
        scope_name,
        sample_name,
        sep =
          "|||"
      )
    ]] <- tibble(
      scope =
        scope_name,
      sample =
        sample_name,
      condition =
        ifelse(
          grepl(
            "^Tx",
            sample_name
          ),
          "Tx",
          "Sham"
        ),
      gene =
        names(
          frac
        ),
      pct_expressed =
        as.numeric(
          frac
        ),
      n_cells =
        length(
          cells
        )
    )
  }
}

positive_fraction <- bind_rows(
  positive_fraction_rows
)

write.csv(
  positive_fraction,
  file.path(
    TAB_OUT,
    "05_direct_candidate_positive_fraction_by_scope_sample_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Pseudobulk counts for each analysis scope
# ==============================================================================

pb_count_list <- list()
pb_meta_list <- list()

for (
  scope_name in SCOPES
) {

  use <- scope_membership[[
    scope_name
  ]]

  cells <- scope_membership$cell[
    use
  ]

  idx <- match(
    cells,
    colnames(
      counts
    )
  )

  groups <- scope_membership$sample[
    use
  ]

  mat_scope <- counts[
    ,
    idx,
    drop =
      FALSE
  ]

  pb <- aggregate_counts_by_group(
    mat_scope,
    groups
  )

  missing_pb_samples <- setdiff(
    SAMPLES,
    colnames(
      pb
    )
  )

  if (
    length(
      missing_pb_samples
    )
  ) {
    stop(
      "Scope ",
      scope_name,
      " is missing sample pseudobulk(s): ",
      paste(
        missing_pb_samples,
        collapse = ", "
      )
    )
  }

  # Force canonical sample order.
  pb <- pb[
    ,
    SAMPLES,
    drop =
      FALSE
  ]

  pb_count_list[[
    scope_name
  ]] <- pb

  n_by_sample <- table(
    factor(
      groups,
      levels =
        SAMPLES
    )
  )

  pb_meta_list[[
    scope_name
  ]] <- tibble(
    scope =
      scope_name,
    sample =
      SAMPLES,
    condition =
      c(
        "Sham",
        "Sham",
        "Tx",
        "Tx"
      ),
    n_cells =
      as.integer(
        n_by_sample
      )
  )
}

pb_meta <- bind_rows(
  pb_meta_list
)

write.csv(
  pb_meta,
  file.path(
    TAB_OUT,
    "06_ballooning_pseudobulk_group_cell_counts_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Scope-wise TMM logCPM and exploratory edgeR
# ==============================================================================

logcpm_rows <- list()
de_rows <- list()

for (
  scope_name in SCOPES
) {

  msg(
    "Pseudobulk analysis: ",
    scope_name
  )

  mat <- pb_count_list[[
    scope_name
  ]]

  meta <- pb_meta_list[[
    scope_name
  ]] %>%
    mutate(
      condition =
        factor(
          condition,
          levels = c(
            "Sham",
            "Tx"
          )
        )
    )

  if (
    any(
      meta$n_cells <
        MIN_CELLS_PB
    )
  ) {
    warning(
      scope_name,
      " contains sample pseudobulk(s) with < ",
      MIN_CELLS_PB,
      " cells."
    )
  }

  # A DGEList using all genes is used for direct CPM reporting.
  y_all <- edgeR::DGEList(
    counts =
      mat
  )

  y_all <- edgeR::calcNormFactors(
    y_all,
    method =
      "TMM"
  )

  direct_cpm <- edgeR::cpm(
    y_all,
    log =
      TRUE,
    prior.count =
      2
  )

  direct_df <- as.data.frame(
    direct_cpm
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    as_tibble() %>%
    pivot_longer(
      cols =
        all_of(
          SAMPLES
        ),
      names_to =
        "sample",
      values_to =
        "logCPM"
    ) %>%
    mutate(
      scope =
        scope_name,
      condition =
        ifelse(
          grepl(
            "^Tx",
            sample
          ),
          "Tx",
          "Sham"
        )
    )

  logcpm_rows[[
    scope_name
  ]] <- direct_df

  # Exploratory edgeR DE uses expression filtering.
  y <- edgeR::DGEList(
    counts =
      mat
  )

  keep <- edgeR::filterByExpr(
    y,
    group =
      meta$condition,
    min.count =
      5
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes =
      FALSE
  ]

  y <- edgeR::calcNormFactors(
    y,
    method =
      "TMM"
  )

  design <- model.matrix(
    ~ condition,
    data =
      meta
  )

  y <- edgeR::estimateDisp(
    y,
    design,
    robust =
      TRUE
  )

  fit <- edgeR::glmQLFit(
    y,
    design,
    robust =
      TRUE
  )

  qlf <- edgeR::glmQLFTest(
    fit,
    coef =
      "conditionTx"
  )

  tab <- edgeR::topTags(
    qlf,
    n =
      Inf,
    sort.by =
      "PValue"
  )$table %>%
    rownames_to_column(
      "gene"
    ) %>%
    as_tibble() %>%
    mutate(
      scope =
        scope_name,
      comparison =
        "Tx_vs_Sham",
      interpretation =
        "Exploratory sample-level pseudobulk; n=2 vs n=2"
    )

  de_rows[[
    scope_name
  ]] <- tab
}

logcpm_all <- bind_rows(
  logcpm_rows
)

de_all <- bind_rows(
  de_rows
)

write.csv(
  logcpm_all,
  file.path(
    TAB_OUT,
    "07_ballooning_pseudobulk_logCPM_all_genes_v6.5.3.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  de_all,
  file.path(
    TAB_OUT,
    "08_ballooning_pseudobulk_DE_all_genes_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Direct candidate gene replicate summaries
# ==============================================================================

direct_summary_rows <- list()

direct_logcpm <- logcpm_all %>%
  filter(
    gene %in%
      direct_present
  )

direct_keys <- direct_logcpm %>%
  distinct(
    scope,
    gene
  )

for (
  i in seq_len(
    nrow(
      direct_keys
    )
  )
) {

  scope_name <- direct_keys$scope[
    i
  ]

  gene_name <- direct_keys$gene[
    i
  ]

  dat <- direct_logcpm %>%
    filter(
      .data$scope ==
        .env$scope_name,
      .data$gene ==
        .env$gene_name
    )

  s <- four_sample_summary(
    dat,
    "logCPM",
    gene_name
  ) %>%
    mutate(
      scope =
        scope_name,
      gene =
        gene_name
    )

  direct_summary_rows[[
    paste(
      scope_name,
      gene_name,
      sep =
        "|||"
    )
  ]] <- s
}

direct_summary <- bind_rows(
  direct_summary_rows
) %>%
  select(
    scope,
    gene,
    everything()
  )

de_focused <- de_all %>%
  filter(
    gene %in%
      direct_present
  ) %>%
  select(
    scope,
    gene,
    logFC,
    logCPM,
    PValue,
    FDR
  )

direct_summary <- direct_summary %>%
  left_join(
    de_focused,
    by = c(
      "scope",
      "gene"
    )
  )

write.csv(
  direct_summary,
  file.path(
    TAB_OUT,
    "09_ballooning_direct_candidate_gene_replicate_summary_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Pseudobulk program scoring by scope
# ==============================================================================

program_gene_audit_rows <- list()
program_sample_rows <- list()

for (
  scope_name in SCOPES
) {

  scope_log <- logcpm_all %>%
    filter(
      .data$scope ==
        .env$scope_name
    ) %>%
    select(
      gene,
      sample,
      logCPM
    ) %>%
    pivot_wider(
      names_from =
        sample,
      values_from =
        logCPM
    )

  mat <- as.matrix(
    scope_log[
      ,
      SAMPLES,
      drop =
        FALSE
    ]
  )

  rownames(
    mat
  ) <- scope_log$gene

  for (
    program_name in names(
      BALLOONING_PROGRAMS
    )
  ) {

    genes_this <- BALLOONING_PROGRAMS[[
      program_name
    ]]

    genes <- intersect(
      genes_this,
      rownames(
        mat
      )
    )

    program_gene_audit_rows[[
      paste(
        scope_name,
        program_name,
        sep =
          "|||"
      )
    ]] <- tibble(
      scope =
        scope_name,
      program =
        program_name,
      gene =
        genes_this,
      present =
        genes_this %in%
          genes
    )

    if (
      length(
        genes
      ) <
        2
    ) {
      next
    }

    zmat <- zscore_rows(
      mat[
        genes,
        ,
        drop =
          FALSE
      ]
    )

    score <- colMeans(
      zmat,
      na.rm =
        TRUE
    )

    program_sample_rows[[
      paste(
        scope_name,
        program_name,
        sep =
          "|||"
      )
    ]] <- tibble(
      scope =
        scope_name,
      program =
        program_name,
      sample =
        names(
          score
        ),
      condition =
        ifelse(
          grepl(
            "^Tx",
            names(
              score
            )
          ),
          "Tx",
          "Sham"
        ),
      score =
        as.numeric(
          score
        ),
      n_genes =
        length(
          genes
        )
    )
  }
}

program_gene_audit <- bind_rows(
  program_gene_audit_rows
)

program_sample <- bind_rows(
  program_sample_rows
)

write.csv(
  program_gene_audit,
  file.path(
    TAB_OUT,
    "10_ballooning_program_gene_audit_by_scope_v6.5.3.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  program_sample,
  file.path(
    TAB_OUT,
    "11_ballooning_program_scores_by_scope_sample_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Program replicate summaries
# ==============================================================================

program_summary_rows <- list()

program_keys <- program_sample %>%
  distinct(
    scope,
    program
  )

for (
  i in seq_len(
    nrow(
      program_keys
    )
  )
) {

  scope_name <- program_keys$scope[
    i
  ]

  program_name <- program_keys$program[
    i
  ]

  dat <- program_sample %>%
    filter(
      .data$scope ==
        .env$scope_name,
      .data$program ==
        .env$program_name
    )

  s <- four_sample_summary(
    dat,
    "score",
    program_name
  ) %>%
    mutate(
      scope =
        scope_name,
      program =
        program_name,
      n_genes =
        first(
          dat$n_genes
        )
    )

  program_summary_rows[[
    paste(
      scope_name,
      program_name,
      sep =
        "|||"
    )
  ]] <- s
}

program_summary <- bind_rows(
  program_summary_rows
) %>%
  select(
    scope,
    program,
    n_genes,
    everything()
  )

write.csv(
  program_summary,
  file.path(
    TAB_OUT,
    "12_ballooning_program_replicate_summary_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Orthogonal candidate gene positive-fraction summaries
# ==============================================================================

positive_summary_rows <- list()

positive_keys <- positive_fraction %>%
  distinct(
    scope,
    gene
  )

for (
  i in seq_len(
    nrow(
      positive_keys
    )
  )
) {

  scope_name <- positive_keys$scope[
    i
  ]

  gene_name <- positive_keys$gene[
    i
  ]

  dat <- positive_fraction %>%
    filter(
      .data$scope ==
        .env$scope_name,
      .data$gene ==
        .env$gene_name
    )

  s <- four_sample_summary(
    dat,
    "pct_expressed",
    paste0(
      gene_name,
      "_positive_fraction"
    )
  ) %>%
    mutate(
      scope =
        scope_name,
      gene =
        gene_name
    )

  positive_summary_rows[[
    paste(
      scope_name,
      gene_name,
      sep =
        "|||"
    )
  ]] <- s
}

positive_summary <- bind_rows(
  positive_summary_rows
) %>%
  select(
    scope,
    gene,
    everything()
  )

write.csv(
  positive_summary,
  file.path(
    TAB_OUT,
    "13_ballooning_candidate_positive_fraction_replicate_summary_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Compact evidence table for histology concordance
# ==============================================================================

# Primary stress/injury proxies: a Tx decrease is considered directionally
# concordant with lower ballooning burden. This is not a diagnostic criterion.

program_evidence <- program_summary %>%
  filter(
    scope %in%
      c(
        "All_primary_Hepatocytes",
        "Injury_inflammatory_Hepatocyte"
      )
  ) %>%
  mutate(
    evidence_type =
      "Molecular_program",
    histology_concordant_direction =
      "Tx_down",
    supports_lower_ballooning =
      case_when(
        evidence_grade %in%
          c(
            "Strong_Tx_down",
            "Moderate_Tx_down"
          ) ~
          "Supportive",
        evidence_grade %in%
          c(
            "Strong_Tx_up",
            "Moderate_Tx_up"
          ) ~
          "Opposite_direction",
        TRUE ~
          "Indeterminate"
      )
  ) %>%
  select(
    scope,
    evidence_type,
    feature =
      program,
    Sham_mean,
    Tx_mean,
    Tx_minus_Sham,
    evidence_grade,
    supports_lower_ballooning
  )

qc_evidence <- qc_summary %>%
  filter(
    metric %in%
      c(
        "MT_high_fraction",
        "pct_all_clean_mt_gt_20"
      )
  ) %>%
  mutate(
    scope =
      "All_clean_Hepatocytes",
    evidence_type =
      "QC_burden",
    feature =
      metric,
    supports_lower_ballooning =
      case_when(
        evidence_grade %in%
          c(
            "Strong_Tx_down",
            "Moderate_Tx_down"
          ) ~
          "Supportive",
        evidence_grade %in%
          c(
            "Strong_Tx_up",
            "Moderate_Tx_up"
          ) ~
          "Opposite_direction",
        TRUE ~
          "Indeterminate"
      )
  ) %>%
  select(
    scope,
    evidence_type,
    feature,
    Sham_mean,
    Tx_mean,
    Tx_minus_Sham,
    evidence_grade,
    supports_lower_ballooning
  )

# Shh is treated separately because it is a direct literature-linked candidate,
# not part of a broad stress module. Krt8/Krt18/Sqstm1 are reported as
# orthogonal candidates but their direction is not used as a diagnostic rule.
shh_evidence <- direct_summary %>%
  filter(
    gene ==
      "Shh",
    scope %in%
      c(
        "All_primary_Hepatocytes",
        "Injury_inflammatory_Hepatocyte"
      )
  ) %>%
  mutate(
    evidence_type =
      "Direct_candidate_gene",
    feature =
      "Shh_logCPM",
    supports_lower_ballooning =
      case_when(
        evidence_grade %in%
          c(
            "Strong_Tx_down",
            "Moderate_Tx_down"
          ) ~
          "Supportive",
        evidence_grade %in%
          c(
            "Strong_Tx_up",
            "Moderate_Tx_up"
          ) ~
          "Opposite_direction",
        TRUE ~
          "Indeterminate"
      )
  ) %>%
  select(
    scope,
    evidence_type,
    feature,
    Sham_mean,
    Tx_mean,
    Tx_minus_Sham,
    evidence_grade,
    supports_lower_ballooning
  )

evidence_table <- bind_rows(
  qc_evidence,
  program_evidence,
  shh_evidence
)

write.csv(
  evidence_table,
  file.path(
    TAB_OUT,
    "14_ballooning_histology_concordance_evidence_table_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Key orthogonal marker table
# ==============================================================================

ORTHOGONAL_GENES <- intersect(
  c(
    "Shh",
    "Ihh",
    "Krt8",
    "Krt18",
    "Sqstm1",
    "Ddit3",
    "Atf3",
    "Cdkn1a",
    "Gadd45a",
    "Saa1",
    "Saa2",
    "Plin2",
    "Cd36"
  ),
  direct_present
)

orthogonal_table <- direct_summary %>%
  filter(
    gene %in%
      ORTHOGONAL_GENES,
    scope %in%
      c(
        "All_primary_Hepatocytes",
        "Injury_inflammatory_Hepatocyte"
      )
  ) %>%
  arrange(
    scope,
    gene
  )

write.csv(
  orthogonal_table,
  file.path(
    TAB_OUT,
    "15_ballooning_key_orthogonal_gene_table_v6.5.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 19. Figures: MT-high/QC fraction
# ==============================================================================

qc_plot_df <- qc_by_sample %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p_mt_fraction <- ggplot(
  qc_plot_df,
  aes(
    x =
      sample,
    y =
      MT_high_fraction,
    group =
      1
  )
) +
  geom_line(
    linewidth =
      0.6
  ) +
  geom_point(
    size =
      2.3
  ) +
  labs(
    title =
      "MT-high/QC Hepatocyte burden",
    subtitle =
      "Sample-level fraction among clean hepatocytes",
    x =
      NULL,
    y =
      "MT-high/QC fraction"
  ) +
  theme_classic(
    base_size =
      9
  )

save_pdf(
  p_mt_fraction,
  file.path(
    FIG_OUT,
    "01_MT_high_QC_fraction_by_sample_v6.5.3.1.pdf"
  ),
  7,
  5
)


# ==============================================================================
# 20. Figures: ballooning proxy programs
# ==============================================================================

program_plot_primary <- program_sample %>%
  filter(
    scope ==
      "All_primary_Hepatocytes"
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p_primary <- ggplot(
  program_plot_primary,
  aes(
    x =
      sample,
    y =
      score,
    group =
      1
  )
) +
  geom_hline(
    yintercept =
      0,
    linewidth =
      0.25
  ) +
  geom_line(
    linewidth =
      0.55
  ) +
  geom_point(
    size =
      2
  ) +
  facet_wrap(
    ~ program,
    scales =
      "free_y",
    ncol =
      3
  ) +
  labs(
    title =
      "Ballooning-related proxy programs | all primary hepatocytes",
    subtitle =
      "Sample-level pseudobulk program scores",
    x =
      NULL,
    y =
      "Program score"
  ) +
  theme_classic(
    base_size =
      8
  )

save_pdf(
  p_primary,
  file.path(
    FIG_OUT,
    "02_ballooning_proxy_programs_all_primary_Hepatocytes_v6.5.3.1.pdf"
  ),
  13,
  9
)

program_plot_injury <- program_sample %>%
  filter(
    scope ==
      "Injury_inflammatory_Hepatocyte"
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p_injury <- ggplot(
  program_plot_injury,
  aes(
    x =
      sample,
    y =
      score,
    group =
      1
  )
) +
  geom_hline(
    yintercept =
      0,
    linewidth =
      0.25
  ) +
  geom_line(
    linewidth =
      0.55
  ) +
  geom_point(
    size =
      2
  ) +
  facet_wrap(
    ~ program,
    scales =
      "free_y",
    ncol =
      3
  ) +
  labs(
    title =
      "Ballooning-related proxy programs | injury/inflammatory hepatocytes",
    subtitle =
      "Sample-level pseudobulk program scores",
    x =
      NULL,
    y =
      "Program score"
  ) +
  theme_classic(
    base_size =
      8
  )

save_pdf(
  p_injury,
  file.path(
    FIG_OUT,
    "03_ballooning_proxy_programs_Injury_Hepatocytes_v6.5.3.1.pdf"
  ),
  13,
  9
)


# ==============================================================================
# 21. Figures: key direct candidates
# ==============================================================================

key_gene_plot <- orthogonal_table %>%
  select(
    scope,
    gene,
    Sham1,
    Sham20,
    Tx17,
    Tx5
  ) %>%
  pivot_longer(
    cols =
      all_of(
        SAMPLES
      ),
    names_to =
      "sample",
    values_to =
      "logCPM"
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

if (
  nrow(
    key_gene_plot
  )
) {

  p_key_gene <- ggplot(
    key_gene_plot,
    aes(
      x =
        sample,
      y =
        logCPM,
      group =
        1
    )
  ) +
    geom_line(
      linewidth =
        0.5
    ) +
    geom_point(
      size =
        1.8
    ) +
    facet_grid(
      scope ~ gene,
      scales =
        "free_y"
    ) +
    labs(
      title =
        "Ballooning-focused orthogonal candidate genes",
      subtitle =
        "TMM-normalized pseudobulk logCPM",
      x =
        NULL,
      y =
        "logCPM"
    ) +
    theme_classic(
      base_size =
        7
    )

  save_pdf(
    p_key_gene,
    file.path(
      FIG_OUT,
      "04_ballooning_key_orthogonal_genes_v6.5.3.1.pdf"
    ),
    18,
    7
  )
}


# ==============================================================================
# 22. Figures: evidence heatmap
# ==============================================================================

evidence_plot <- evidence_table %>%
  mutate(
    feature =
      factor(
        feature,
        levels =
          unique(
            feature
          )
      )
  )

p_evidence <- ggplot(
  evidence_plot,
  aes(
    x =
      feature,
    y =
      scope,
    fill =
      Tx_minus_Sham
  )
) +
  geom_tile(
    linewidth =
      0.3
  ) +
  scale_fill_gradient2(
    low =
      "#0033FF",
    mid =
      "#FFFFFF",
    high =
      "#FF1A1A",
    midpoint =
      0
  ) +
  labs(
    title =
      "Ballooning-focused evidence summary | Tx - Sham",
    subtitle =
      "Blue = lower in Tx; red = higher in Tx",
    x =
      NULL,
    y =
      NULL,
    fill =
      "Tx-Sham"
  ) +
  theme_classic(
    base_size =
      8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )

save_pdf(
  p_evidence,
  file.path(
    FIG_OUT,
    "05_ballooning_histology_concordance_heatmap_v6.5.3.1.pdf"
  ),
  14,
  6
)


# ==============================================================================
# 23. Save augmented RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_Hepatocyte_ballooning_focused_v6.5.3.1.rds"
)

saveRDS(
  hep,
  RDS_FILE,
  compress =
    FALSE
)

msg(
  "Saved ballooning-focused Hepatocyte RDS: ",
  RDS_FILE
)


# ==============================================================================
# 24. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "final_state_column",
    "primary_scope",
    "injury_scope",
    "QC_scope",
    "primary_interpretation",
    "pseudobulk_method",
    "comparison",
    "replicates",
    "ballooning_claim"
  ),
  value = c(
    "v6.5.3.1",
    INPUT_RDS,
    STATE_COL,
    "All final Hepatocyte states except MT_high_QC_Hepatocyte",
    INJURY_STATE,
    MT_STATE,
    "Evaluate molecular/QC concordance with histological ballooning reduction",
    "Sample-level raw-count aggregation + edgeR TMM; gene-wise z-score programs",
    "Tx vs Sham",
    "Sham1, Sham20 vs Tx17, Tx5",
    "scRNA-seq proxies only; histology remains the direct ballooning endpoint"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.5.3.1.csv"
  ),
  row.names =
    FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.5.3.1.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
