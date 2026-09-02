#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6330)

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
# Population-weighted PDGFB output and HSC target-module validation
#
# Version: v6.3.3.1.1
#
# PRIMARY QUESTION
#   Repair/Resolution-Mphi increases after transplantation, whereas Pdgfb
#   expression per Repair/Resolution-Mphi cell decreases.
#
#   Which effect dominates the population-level macrophage Pdgfb output?
#
# INPUTS
#   1) v6.2.0 interaction-ready Seurat object
#   2) v6.3.2 PDGFB target-gene set
#
# DESIGN
#
#   A. Mphi abundance
#      subtype fraction among the five macrophage populations
#
#   B. Per-cell ligand output
#      Pdgfb count normalized per cell:
#
#        Pdgfb_CP10k =
#          10,000 * Pdgfb raw count / total RNA counts in that cell
#
#      This is linear-scale normalized expression and is therefore preferable
#      to multiplying a cell fraction by log-normalized expression.
#
#   C. Population-weighted Pdgfb output proxy
#
#        weighted_output(subtype, sample) =
#          subtype_fraction_within_Mphi5 *
#          mean_Pdgfb_CP10k_within_subtype
#
#      Sum over all five macrophage subtypes =
#          mean Pdgfb CP10k per macrophage cell in that sample.
#
#   D. Exact decomposition of the observed Sham -> Tx mean change
#
#      Let f = subtype fraction and e = per-cell Pdgfb CP10k.
#
#      composition_effect =
#          (mean(f)_Tx - mean(f)_Sham) *
#          (mean(e)_Tx + mean(e)_Sham) / 2
#
#      expression_effect =
#          (mean(e)_Tx - mean(e)_Sham) *
#          (mean(f)_Tx + mean(f)_Sham) / 2
#
#      The above two terms exactly reconstruct the change in the PRODUCT OF
#      CONDITION MEANS, but not necessarily the change in mean(f*e), because
#      the two biological samples within a condition can show covariance
#      between abundance and per-cell expression.
#
#      covariance_effect =
#          [mean(f*e)_Tx - mean(f)_Tx*mean(e)_Tx] -
#          [mean(f*e)_Sham - mean(f)_Sham*mean(e)_Sham]
#
#      Therefore:
#
#          observed change in mean(f*e)
#            = composition_effect
#            + expression_effect
#            + covariance_effect
#
#      With n=2 per condition, the covariance term is descriptive only.
#
#   E. ECM-activated HSC PDGFB-target module
#      For each biological sample:
#        - pseudobulk CPM for each v6.3.2 target gene
#        - log2(CPM+1)
#        - gene-wise z-score across the four samples
#        - mean z-score across targets = target-module score
#
#   F. Descriptive sender-receiver coupling
#      Across the four biological samples only:
#        total Mphi5 weighted Pdgfb output vs HSC target-module score
#        Repair/Resolution weighted Pdgfb output vs target-module score
#
# IMPORTANT
#   - No CellChat rerun.
#   - No NicheNet rerun.
#   - No reclustering / reintegration.
#   - No cell-level pseudoreplication.
#   - Biological n = 2 Sham vs n = 2 Tx.
#   - Correlations with n = 4 are descriptive only.
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

V632 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_axis_validation_v6.3.2"
)

TARGET_FILE <- file.path(
  V632,
  "Tables",
  "10_PDGFB_target_gene_set_v6.3.2.csv"
)

V632_SUMMARY_FILE <- file.path(
  V632,
  "Tables",
  "15_PDGFB_axis_validation_summary_v6.3.2.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_population_weighted_v6.3.3.1"
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

CONDITION_COL <-
  "condition_interaction_v620"

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

FOCAL_SENDER <-
  "Repair/Resolution-Mphi"

FOCAL_RECEIVER <-
  "ECM-activated HSC"

LIGAND <-
  "Pdgfb"

DEFAULT_TARGETS <- c(
  "Cd24a",
  "Aff3",
  "Sfrp1",
  "Ypel1",
  "Atp1b1",
  "Kcnma1",
  "Itga4",
  "Slfn4",
  "Has2",
  "Il7r"
)

FOCUSED_TARGETS <- c(
  "Sfrp1",
  "Kcnma1",
  "Has2",
  "Itga4",
  "Atp1b1"
)

CONDITION_COLORS <- c(
  "Sham" = "#0072B2",
  "Tx" = "#D55E00"
)


# ==============================================================================
# 4. Preflight / load
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

if (
  !file.exists(
    TARGET_FILE
  )
) {
  stop(
    "v6.3.2 target file missing: ",
    TARGET_FILE
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
  SAMPLE_COL,
  CONDITION_COL
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

if (
  !LIGAND %in%
    rownames(
      counts
    )
) {
  stop(
    "Pdgfb is absent from RNA counts."
  )
}

target_info <- read.csv(
  TARGET_FILE,
  check.names = FALSE
) %>%
  as_tibble()

if (
  !"target" %in%
    colnames(
      target_info
    )
) {
  stop(
    "Target file has no 'target' column."
  )
}

TARGETS <- unique(
  target_info$target
)

TARGETS <- TARGETS[
  TARGETS %in%
    rownames(
      counts
    )
]

if (
  !length(
    TARGETS
  )
) {
  TARGETS <- intersect(
    DEFAULT_TARGETS,
    rownames(
      counts
    )
  )
}

if (
  !length(
    TARGETS
  )
) {
  stop(
    "No PDGFB target genes found in RNA assay."
  )
}

FOCUSED_PRESENT <- intersect(
  FOCUSED_TARGETS,
  TARGETS
)

msg(
  "PDGFB target genes: ",
  paste(
    TARGETS,
    collapse = ", "
  )
)


# ==============================================================================
# 5. Per-cell linear normalized Pdgfb expression
#
# CP10k = 10,000 * raw Pdgfb count / cell total RNA counts
# ==============================================================================

msg(
  "Computing per-cell linear normalized Pdgfb..."
)

cell_library_size <- Matrix::colSums(
  counts
)

pdgfb_raw <- as.numeric(
  counts[
    LIGAND,
    ,
    drop = TRUE
  ]
)

names(
  pdgfb_raw
) <- colnames(
  counts
)

pdgfb_cp10k <- 10000 *
  pdgfb_raw /
  pmax(
    as.numeric(
      cell_library_size
    ),
    1
  )

names(
  pdgfb_cp10k
) <- colnames(
  counts
)

cell_pdgfb <- tibble(
  cell =
    colnames(
      counts
    ),
  Pdgfb_raw_count =
    pdgfb_raw,
  cell_RNA_counts =
    as.numeric(
      cell_library_size
    ),
  Pdgfb_CP10k =
    pdgfb_cp10k
) %>%
  left_join(
    meta %>%
      as_tibble(
        rownames = "cell"
      ) %>%
      select(
        cell,
        all_of(
          c(
            GROUP_COL,
            SAMPLE_COL,
            CONDITION_COL
          )
        )
      ),
    by =
      "cell"
  )

write.csv(
  cell_pdgfb %>%
    filter(
      .data[[
        GROUP_COL
      ]] %in%
        MPHI5
    ),
  file.path(
    TAB_OUT,
    "01_Mphi5_celllevel_Pdgfb_CP10k_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. Sample x Mphi subtype abundance and per-cell expression
# ==============================================================================

msg(
  "Aggregating Mphi subtype abundance and per-cell Pdgfb..."
)

mphi_cell <- cell_pdgfb %>%
  filter(
    .data[[
      GROUP_COL
    ]] %in%
      MPHI5,
    .data[[
      SAMPLE_COL
    ]] %in%
      SAMPLES
  ) %>%
  mutate(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    condition =
      canonical_condition(
        sample
      ),
    Mphi_subtype =
      as.character(
        .data[[
          GROUP_COL
        ]]
      )
  )

mphi_summary <- mphi_cell %>%
  group_by(
    sample,
    condition,
    Mphi_subtype
  ) %>%
  summarise(
    n_cells =
      n(),

    Pdgfb_positive_cells =
      sum(
        Pdgfb_raw_count >
          0
      ),

    Pdgfb_pct_positive =
      mean(
        Pdgfb_raw_count >
          0
      ),

    mean_Pdgfb_CP10k =
      mean(
        Pdgfb_CP10k
      ),

    median_Pdgfb_CP10k =
      stats::median(
        Pdgfb_CP10k
      ),

    sum_Pdgfb_raw_counts =
      sum(
        Pdgfb_raw_count
      ),

    .groups = "drop"
  ) %>%
  group_by(
    sample
  ) %>%
  mutate(
    total_Mphi5_cells =
      sum(
        n_cells
      ),

    fraction_within_Mphi5 =
      n_cells /
        total_Mphi5_cells,

    weighted_Pdgfb_output =
      fraction_within_Mphi5 *
        mean_Pdgfb_CP10k
  ) %>%
  ungroup()

# Contribution to the sample's total weighted Pdgfb output.
mphi_summary <- mphi_summary %>%
  group_by(
    sample
  ) %>%
  mutate(
    total_weighted_Pdgfb_output =
      sum(
        weighted_Pdgfb_output
      ),

    contribution_fraction =
      ifelse(
        total_weighted_Pdgfb_output >
          0,
        weighted_Pdgfb_output /
          total_weighted_Pdgfb_output,
        0
      )
  ) %>%
  ungroup()

write.csv(
  mphi_summary,
  file.path(
    TAB_OUT,
    "02_Mphi5_population_weighted_Pdgfb_by_sample_subtype_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Sample-level total macrophage Pdgfb output
# ==============================================================================

total_mphi_output <- mphi_summary %>%
  group_by(
    sample,
    condition
  ) %>%
  summarise(
    total_Mphi5_cells =
      first(
        total_Mphi5_cells
      ),

    total_weighted_Pdgfb_output =
      first(
        total_weighted_Pdgfb_output
      ),

    direct_mean_Pdgfb_CP10k_all_Mphi5 =
      weighted.mean(
        mean_Pdgfb_CP10k,
        w =
          n_cells
      ),

    total_Pdgfb_raw_counts_Mphi5 =
      sum(
        sum_Pdgfb_raw_counts
      ),

    .groups = "drop"
  ) %>%
  mutate(
    identity_check_difference =
      total_weighted_Pdgfb_output -
        direct_mean_Pdgfb_CP10k_all_Mphi5
  )

if (
  any(
    abs(
      total_mphi_output$identity_check_difference
    ) >
      1e-10
  )
) {
  stop(
    "Weighted-output identity check failed."
  )
}

write.csv(
  total_mphi_output,
  file.path(
    TAB_OUT,
    "03_total_Mphi5_population_weighted_Pdgfb_by_sample_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Focal Repair/Resolution-Mphi sample values
# ==============================================================================

repair_summary <- mphi_summary %>%
  filter(
    Mphi_subtype ==
      FOCAL_SENDER
  )

write.csv(
  repair_summary,
  file.path(
    TAB_OUT,
    "04_RepairResolution_population_weighted_Pdgfb_by_sample_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Sham-vs-Tx mean decomposition by Mphi subtype
# ==============================================================================

msg(
  "Decomposing composition and per-cell expression effects..."
)

mphi_condition_mean <- mphi_summary %>%
  group_by(
    condition,
    Mphi_subtype
  ) %>%
  summarise(
    mean_fraction =
      mean(
        fraction_within_Mphi5
      ),

    mean_percell_Pdgfb_CP10k =
      mean(
        mean_Pdgfb_CP10k
      ),

    mean_weighted_output =
      mean(
        weighted_Pdgfb_output
      ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from =
      condition,
    values_from = c(
      mean_fraction,
      mean_percell_Pdgfb_CP10k,
      mean_weighted_output
    )
  )

decomposition <- mphi_condition_mean %>%
  mutate(
    delta_fraction =
      mean_fraction_Tx -
        mean_fraction_Sham,

    delta_percell_expression =
      mean_percell_Pdgfb_CP10k_Tx -
        mean_percell_Pdgfb_CP10k_Sham,

    observed_delta_weighted_output =
      mean_weighted_output_Tx -
        mean_weighted_output_Sham,

    composition_effect =
      delta_fraction *
        (
          mean_percell_Pdgfb_CP10k_Tx +
            mean_percell_Pdgfb_CP10k_Sham
        ) /
          2,

    expression_effect =
      delta_percell_expression *
        (
          mean_fraction_Tx +
            mean_fraction_Sham
        ) /
          2,

    product_of_means_Sham =
      mean_fraction_Sham *
        mean_percell_Pdgfb_CP10k_Sham,

    product_of_means_Tx =
      mean_fraction_Tx *
        mean_percell_Pdgfb_CP10k_Tx,

    covariance_component_Sham =
      mean_weighted_output_Sham -
        product_of_means_Sham,

    covariance_component_Tx =
      mean_weighted_output_Tx -
        product_of_means_Tx,

    covariance_effect =
      covariance_component_Tx -
        covariance_component_Sham,

    reconstructed_delta =
      composition_effect +
        expression_effect +
        covariance_effect,

    decomposition_error =
      observed_delta_weighted_output -
        reconstructed_delta
  )

if (
  any(
    abs(
      decomposition$decomposition_error
    ) >
      1e-10
  )
) {
  stop(
    "Three-component decomposition identity check failed."
  )
}

write.csv(
  decomposition,
  file.path(
    TAB_OUT,
    "05_Mphi5_Pdgfb_composition_vs_expression_decomposition_v6.3.3.1.csv"
  ),
  row.names = FALSE
)

decomposition_long <- decomposition %>%
  select(
    Mphi_subtype,
    composition_effect,
    expression_effect,
    covariance_effect
  ) %>%
  pivot_longer(
    cols = c(
      composition_effect,
      expression_effect,
      covariance_effect
    ),
    names_to =
      "effect_type",
    values_to =
      "delta_weighted_output"
  )

write.csv(
  decomposition_long,
  file.path(
    TAB_OUT,
    "06_Mphi5_Pdgfb_decomposition_long_v6.3.3.1.csv"
  ),
  row.names = FALSE
)

decomposition_diagnostic <- decomposition %>%
  select(
    Mphi_subtype,
    mean_fraction_Sham,
    mean_fraction_Tx,
    mean_percell_Pdgfb_CP10k_Sham,
    mean_percell_Pdgfb_CP10k_Tx,
    mean_weighted_output_Sham,
    mean_weighted_output_Tx,
    product_of_means_Sham,
    product_of_means_Tx,
    covariance_component_Sham,
    covariance_component_Tx,
    composition_effect,
    expression_effect,
    covariance_effect,
    observed_delta_weighted_output,
    reconstructed_delta,
    decomposition_error
  )

write.csv(
  decomposition_diagnostic,
  file.path(
    TAB_OUT,
    "06b_Mphi5_Pdgfb_decomposition_diagnostic_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. HSC target pseudobulk CPM
# ==============================================================================

msg(
  "Computing ECM-activated HSC PDGFB target pseudobulk..."
)

pb_hsc <- sample_pseudobulk(
  counts = counts,
  meta = meta,
  sample_col = SAMPLE_COL,
  group_col = GROUP_COL,
  group_value = FOCAL_RECEIVER,
  samples = SAMPLES
)

hsc_lib <- Matrix::colSums(
  pb_hsc
)

target_cpm_mat <- sapply(
  SAMPLES,
  function(smp) {

    1e6 *
      as.numeric(
        pb_hsc[
          TARGETS,
          smp,
          drop = TRUE
        ]
      ) /
      max(
        hsc_lib[[
          smp
        ]],
        1
      )
  }
)

if (
  is.null(
    dim(
      target_cpm_mat
    )
  )
) {
  target_cpm_mat <- matrix(
    target_cpm_mat,
    nrow =
      length(
        TARGETS
      ),
    dimnames = list(
      TARGETS,
      SAMPLES
    )
  )
} else {
  rownames(
    target_cpm_mat
  ) <- TARGETS
  colnames(
    target_cpm_mat
  ) <- SAMPLES
}

target_logcpm_mat <- log2(
  target_cpm_mat +
    1
)

target_z_mat <- t(
  apply(
    target_logcpm_mat,
    1,
    row_z
  )
)

target_long <- as.data.frame(
  target_logcpm_mat
) %>%
  rownames_to_column(
    "target"
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

target_z_long <- as.data.frame(
  target_z_mat
) %>%
  rownames_to_column(
    "target"
  ) %>%
  pivot_longer(
    cols =
      all_of(
        SAMPLES
      ),
    names_to =
      "sample",
    values_to =
      "z"
  )

target_values <- target_long %>%
  left_join(
    target_z_long,
    by = c(
      "target",
      "sample"
    )
  ) %>%
  mutate(
    condition =
      canonical_condition(
        sample
      ),
    focused =
      target %in%
        FOCUSED_PRESENT
  )

write.csv(
  target_values,
  file.path(
    TAB_OUT,
    "07_ECM_HSC_PDGFB_target_pseudobulk_values_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. HSC target module scores
# ==============================================================================

all_target_module <- target_values %>%
  group_by(
    sample,
    condition
  ) %>%
  summarise(
    PDGFB_target_module_z_all =
      mean(
        z,
        na.rm = TRUE
      ),

    mean_target_log2CPM1_all =
      mean(
        log2_CPM1,
        na.rm = TRUE
      ),

    n_targets =
      n_distinct(
        target
      ),

    .groups = "drop"
  )

if (
  length(
    FOCUSED_PRESENT
  )
) {

  focused_target_module <- target_values %>%
    filter(
      focused
    ) %>%
    group_by(
      sample,
      condition
    ) %>%
    summarise(
      PDGFB_target_module_z_focused =
        mean(
          z,
          na.rm = TRUE
        ),

      mean_target_log2CPM1_focused =
        mean(
          log2_CPM1,
          na.rm = TRUE
        ),

      n_focused_targets =
        n_distinct(
          target
        ),

      .groups = "drop"
    )

} else {

  focused_target_module <- all_target_module %>%
    transmute(
      sample,
      condition,
      PDGFB_target_module_z_focused =
        NA_real_,
      mean_target_log2CPM1_focused =
        NA_real_,
      n_focused_targets =
        0L
    )
}

target_module <- all_target_module %>%
  left_join(
    focused_target_module,
    by = c(
      "sample",
      "condition"
    )
  )

write.csv(
  target_module,
  file.path(
    TAB_OUT,
    "08_ECM_HSC_PDGFB_target_module_by_sample_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Integrated sample-level sender-receiver table
# ==============================================================================

sample_axis <- total_mphi_output %>%
  select(
    sample,
    condition,
    total_Mphi5_cells,
    total_Mphi5_weighted_Pdgfb =
      total_weighted_Pdgfb_output
  ) %>%
  left_join(
    repair_summary %>%
      select(
        sample,
        RepairResolution_n_cells =
          n_cells,
        RepairResolution_fraction =
          fraction_within_Mphi5,
        RepairResolution_mean_Pdgfb_CP10k =
          mean_Pdgfb_CP10k,
        RepairResolution_Pdgfb_pct_positive =
          Pdgfb_pct_positive,
        RepairResolution_weighted_Pdgfb =
          weighted_Pdgfb_output,
        RepairResolution_contribution_fraction =
          contribution_fraction
      ),
    by =
      "sample"
  ) %>%
  left_join(
    target_module,
    by = c(
      "sample",
      "condition"
    )
  )

write.csv(
  sample_axis,
  file.path(
    TAB_OUT,
    "09_integrated_PDGFB_sender_receiver_by_sample_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Descriptive correlations across four samples
# ==============================================================================

correlation_summary <- tibble(
  comparison = c(
    "Total_Mphi5_weighted_Pdgfb_vs_all_target_module",
    "RepairResolution_weighted_Pdgfb_vs_all_target_module",
    "RepairResolution_percell_Pdgfb_vs_all_target_module",
    "RepairResolution_fraction_vs_all_target_module",
    "Total_Mphi5_weighted_Pdgfb_vs_focused_target_module",
    "RepairResolution_weighted_Pdgfb_vs_focused_target_module"
  ),

  spearman_rho = c(
    safe_spearman(
      sample_axis$total_Mphi5_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_all
    ),

    safe_spearman(
      sample_axis$RepairResolution_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_all
    ),

    safe_spearman(
      sample_axis$RepairResolution_mean_Pdgfb_CP10k,
      sample_axis$PDGFB_target_module_z_all
    ),

    safe_spearman(
      sample_axis$RepairResolution_fraction,
      sample_axis$PDGFB_target_module_z_all
    ),

    safe_spearman(
      sample_axis$total_Mphi5_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_focused
    ),

    safe_spearman(
      sample_axis$RepairResolution_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_focused
    )
  ),

  n_samples =
    4L,

  interpretation =
    "Descriptive only; no inferential claim with n=4"
)

write.csv(
  correlation_summary,
  file.path(
    TAB_OUT,
    "10_descriptive_PDGFB_sender_receiver_correlations_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Condition means / summary
# ==============================================================================

condition_axis <- sample_axis %>%
  group_by(
    condition
  ) %>%
  summarise(
    n_samples =
      n(),

    mean_total_Mphi5_weighted_Pdgfb =
      mean(
        total_Mphi5_weighted_Pdgfb
      ),

    mean_RepairResolution_fraction =
      mean(
        RepairResolution_fraction
      ),

    mean_RepairResolution_percell_Pdgfb =
      mean(
        RepairResolution_mean_Pdgfb_CP10k
      ),

    mean_RepairResolution_weighted_Pdgfb =
      mean(
        RepairResolution_weighted_Pdgfb
      ),

    mean_target_module_all =
      mean(
        PDGFB_target_module_z_all
      ),

    mean_target_module_focused =
      mean(
        PDGFB_target_module_z_focused,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  condition_axis,
  file.path(
    TAB_OUT,
    "11_PDGFB_axis_condition_means_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Figure 1:
#     Mphi subtype contributions to total Pdgfb output
# ==============================================================================

contribution_plot_df <- mphi_summary %>%
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
  )

p1 <- ggplot(
  contribution_plot_df,
  aes(
    x =
      sample,
    y =
      weighted_Pdgfb_output,
    fill =
      Mphi_subtype
  )
) +
  geom_col(
    width = 0.75
  ) +
  labs(
    title =
      "Population-weighted macrophage Pdgfb output",
    subtitle =
      "Subtype fraction x mean per-cell Pdgfb CP10k; stacked sum = total Mphi5 output",
    x = NULL,
    y =
      "Weighted Pdgfb CP10k",
    fill =
      "Mphi subtype"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_Mphi5_population_weighted_Pdgfb_contributions_v6.3.3.1.pdf"
  ),
  10,
  7
)


# ==============================================================================
# 16. Figure 2:
#     Composition vs expression decomposition
# ==============================================================================

decomp_plot_df <- decomposition_long %>%
  mutate(
    Mphi_subtype =
      factor(
        Mphi_subtype,
        levels =
          MPHI5
      ),

    effect_type =
      factor(
        effect_type,
        levels = c(
          "composition_effect",
          "expression_effect",
          "covariance_effect"
        )
      )
  )

p2 <- ggplot(
  decomp_plot_df,
  aes(
    x =
      Mphi_subtype,
    y =
      delta_weighted_output,
    fill =
      effect_type
  )
) +
  geom_col(
    position = "stack"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35
  ) +
  labs(
    title =
      "Sham -> Tx change in weighted Pdgfb output",
    subtitle =
      "Exact decomposition into composition, per-cell expression, and within-condition covariance effects",
    x = NULL,
    y =
      "Delta weighted Pdgfb CP10k",
    fill =
      "Effect"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_Pdgfb_composition_vs_expression_decomposition_v6.3.3.1.pdf"
  ),
  11,
  7
)


# ==============================================================================
# 17. Figure 3:
#     Repair/Resolution paradox
# ==============================================================================

repair_plot <- sample_axis %>%
  select(
    sample,
    condition,
    RepairResolution_fraction,
    RepairResolution_mean_Pdgfb_CP10k,
    RepairResolution_weighted_Pdgfb
  )

p3a <- ggplot(
  repair_plot,
  aes(
    x =
      sample,
    y =
      100 *
        RepairResolution_fraction,
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
      "Repair/Resolution-Mphi abundance",
    x = NULL,
    y =
      "% of Mphi5",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p3b <- ggplot(
  repair_plot,
  aes(
    x =
      sample,
    y =
      RepairResolution_mean_Pdgfb_CP10k,
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
      "Per-cell Pdgfb output",
    x = NULL,
    y =
      "Mean Pdgfb CP10k",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p3c <- ggplot(
  repair_plot,
  aes(
    x =
      sample,
    y =
      RepairResolution_weighted_Pdgfb,
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
      "Population-weighted Repair/Resolution Pdgfb",
    x = NULL,
    y =
      "Weighted Pdgfb CP10k",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p3 <- p3a +
  p3b +
  p3c +
  plot_annotation(
    title =
      "Repair/Resolution-Mphi: abundance vs per-cell Pdgfb vs net weighted output"
  )

save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_RepairResolution_abundance_expression_weighted_output_v6.3.3.1.pdf"
  ),
  15,
  5
)


# ==============================================================================
# 18. Figure 4:
#     ECM-HSC target module
# ==============================================================================

module_plot_df <- sample_axis %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      )
  )

p4a <- ggplot(
  module_plot_df,
  aes(
    x =
      sample,
    y =
      PDGFB_target_module_z_all,
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
      "All predicted PDGFB targets",
    x = NULL,
    y =
      "Mean gene-wise z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p4b <- ggplot(
  module_plot_df,
  aes(
    x =
      sample,
    y =
      PDGFB_target_module_z_focused,
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
      "Focused HSC-relevant PDGFB targets",
    x = NULL,
    y =
      "Mean gene-wise z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p4 <- p4a +
  p4b +
  plot_annotation(
    title =
      "ECM-activated HSC PDGFB-target module"
  )

save_pdf(
  p4,
  file.path(
    FIG_OUT,
    "04_ECM_HSC_PDGFB_target_module_by_sample_v6.3.3.1.pdf"
  ),
  12,
  5
)


# ==============================================================================
# 19. Figure 5:
#     Descriptive sender-receiver coupling
# ==============================================================================

p5a <- ggplot(
  sample_axis,
  aes(
    x =
      total_Mphi5_weighted_Pdgfb,
    y =
      PDGFB_target_module_z_all,
    fill =
      condition
  )
) +
  geom_point(
    shape = 21,
    size = 4
  ) +
  geom_text(
    aes(
      label =
        sample
    ),
    nudge_y = 0.08,
    size = 3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Total Mphi5 Pdgfb vs ECM-HSC target module",
    subtitle =
      "Four biological samples; descriptive only",
    x =
      "Total Mphi5 weighted Pdgfb CP10k",
    y =
      "PDGFB target module z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p5b <- ggplot(
  sample_axis,
  aes(
    x =
      RepairResolution_weighted_Pdgfb,
    y =
      PDGFB_target_module_z_all,
    fill =
      condition
  )
) +
  geom_point(
    shape = 21,
    size = 4
  ) +
  geom_text(
    aes(
      label =
        sample
    ),
    nudge_y = 0.08,
    size = 3
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Repair/Resolution Pdgfb vs ECM-HSC target module",
    subtitle =
      "Four biological samples; descriptive only",
    x =
      "Repair/Resolution weighted Pdgfb CP10k",
    y =
      "PDGFB target module z-score",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p5 <- p5a +
  p5b

save_pdf(
  p5,
  file.path(
    FIG_OUT,
    "05_PDGFB_sender_receiver_descriptive_coupling_v6.3.3.1.pdf"
  ),
  12,
  5.5
)


# ==============================================================================
# 20. Figure 6:
#     Integrated validation panel
# ==============================================================================

p6a <- ggplot(
  sample_axis,
  aes(
    x =
      sample,
    y =
      total_Mphi5_weighted_Pdgfb,
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
      "Total Mphi5 weighted Pdgfb",
    x = NULL,
    y =
      "Weighted CP10k",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 8
  )

p6b <- p3a +
  theme(
    legend.position = "none"
  )

p6c <- p3b +
  theme(
    legend.position = "none"
  )

p6d <- p3c +
  theme(
    legend.position = "none"
  )

p6e <- p4a +
  theme(
    legend.position = "none"
  )

p6f <- p5b +
  theme(
    legend.position = "none"
  )

p6 <- (
  p6a +
    p6b +
    p6c
) /
  (
    p6d +
      p6e +
      p6f
  ) +
  plot_annotation(
    title =
      "PDGFB axis integrated validation | v6.3.3.1"
  )

save_pdf(
  p6,
  file.path(
    FIG_OUT,
    "06_PDGFB_axis_integrated_validation_panel_v6.3.3.1.pdf"
  ),
  16,
  10
)


# ==============================================================================
# 21. Final interpretation summary
# ==============================================================================

repair_decomp <- decomposition %>%
  filter(
    Mphi_subtype ==
      FOCAL_SENDER
  )

sham_total <- condition_axis %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  pull(
    mean_total_Mphi5_weighted_Pdgfb
  )

tx_total <- condition_axis %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  pull(
    mean_total_Mphi5_weighted_Pdgfb
  )

sham_repair_weighted <- condition_axis %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  pull(
    mean_RepairResolution_weighted_Pdgfb
  )

tx_repair_weighted <- condition_axis %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  pull(
    mean_RepairResolution_weighted_Pdgfb
  )

sham_target <- condition_axis %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  pull(
    mean_target_module_all
  )

tx_target <- condition_axis %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  pull(
    mean_target_module_all
  )

interpretation_summary <- tibble(
  metric = c(
    "Mean_total_Mphi5_weighted_Pdgfb_Sham",
    "Mean_total_Mphi5_weighted_Pdgfb_Tx",
    "Delta_total_Mphi5_weighted_Pdgfb_Tx_minus_Sham",
    "Mean_RepairResolution_weighted_Pdgfb_Sham",
    "Mean_RepairResolution_weighted_Pdgfb_Tx",
    "Delta_RepairResolution_weighted_Pdgfb_Tx_minus_Sham",
    "RepairResolution_composition_effect",
    "RepairResolution_expression_effect",
    "RepairResolution_covariance_effect",
    "RepairResolution_observed_delta_weighted_output",
    "Mean_ECM_HSC_target_module_Sham",
    "Mean_ECM_HSC_target_module_Tx",
    "Delta_ECM_HSC_target_module_Tx_minus_Sham",
    "Spearman_total_Mphi5_Pdgfb_vs_target_module",
    "Spearman_RepairResolution_Pdgfb_vs_target_module"
  ),

  value = c(
    sham_total,
    tx_total,
    tx_total -
      sham_total,
    sham_repair_weighted,
    tx_repair_weighted,
    tx_repair_weighted -
      sham_repair_weighted,
    if (
      nrow(
        repair_decomp
      )
    ) {
      repair_decomp$composition_effect[[1]]
    } else {
      NA_real_
    },
    if (
      nrow(
        repair_decomp
      )
    ) {
      repair_decomp$expression_effect[[1]]
    } else {
      NA_real_
    },
    if (
      nrow(
        repair_decomp
      )
    ) {
      repair_decomp$covariance_effect[[1]]
    } else {
      NA_real_
    },
    if (
      nrow(
        repair_decomp
      )
    ) {
      repair_decomp$observed_delta_weighted_output[[1]]
    } else {
      NA_real_
    },
    sham_target,
    tx_target,
    tx_target -
      sham_target,
    safe_spearman(
      sample_axis$total_Mphi5_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_all
    ),
    safe_spearman(
      sample_axis$RepairResolution_weighted_Pdgfb,
      sample_axis$PDGFB_target_module_z_all
    )
  ),

  note = c(
    rep(
      "Descriptive biological-sample summary",
      13
    ),
    rep(
      "Spearman rho with n=4; descriptive only",
      2
    )
  )
)

write.csv(
  interpretation_summary,
  file.path(
    TAB_OUT,
    "12_PDGFB_axis_interpretation_summary_v6.3.3.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 22. Save compact results RDS
# ==============================================================================

results <- list(
  Mphi5_population_weighted =
    mphi_summary,
  total_Mphi5_output =
    total_mphi_output,
  RepairResolution =
    repair_summary,
  decomposition =
    decomposition,
  target_values =
    target_values,
  target_module =
    target_module,
  sample_axis =
    sample_axis,
  correlations =
    correlation_summary,
  condition_means =
    condition_axis,
  interpretation_summary =
    interpretation_summary
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_PDGFB_population_weighted_results_v6.3.3.1.rds"
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
    "target_file",
    "ligand",
    "focal_sender",
    "focal_receiver",
    "per_cell_normalization",
    "weighted_output_definition",
    "decomposition",
    "target_module",
    "biological_replicates",
    "CellChat_recomputed",
    "NicheNet_recomputed",
    "formal_inference_note"
  ),

  value = c(
    "v6.3.3.1",
    INPUT_RDS,
    TARGET_FILE,
    LIGAND,
    FOCAL_SENDER,
    FOCAL_RECEIVER,
    "Pdgfb CP10k = 10000 * raw Pdgfb count / total RNA counts per cell",
    "Mphi subtype fraction x mean per-cell Pdgfb CP10k",
    "Exact decomposition into composition, expression, and within-condition covariance effects",
    "Mean gene-wise z-score of ECM-HSC target pseudobulk log2(CPM+1)",
    "Sham n=2; Tx n=2",
    "FALSE",
    "FALSE",
    "Exploratory; all correlations n=4 descriptive only"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.3.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.3.1.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
