#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6520)

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

# edgeR is used only for sample-level pseudobulk DE.
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "Package 'edgeR' is required for v6.5.2. ",
    "Install it inside the project renv before running."
  )
}

# ==============================================================================
# Mouse MASH scRNA-seq
# Final clean Hepatocyte state annotation + sample-level pseudobulk analysis
#
# Version: v6.5.2
#
# INPUT:
#   Mouse_MASH_Hepatocyte_highconf_clean_reclustered_v6.5.1.rds
#
# FINAL WORKING ANNOTATION (clean res0.4):
#   0 = Periportal_Hepatocyte_1
#   1 = Injury_inflammatory_Hepatocyte
#   2 = Pericentral_Hepatocyte
#   3 = MT_high_QC_Hepatocyte
#   4 = Periportal_Hepatocyte_2
#   5 = Intermediate_Hepatocyte
#   6 = Cycling_G2M_Hepatocyte
#   7 = Cycling_S_Hepatocyte
#
# IMPORTANT INTERPRETATION RULES:
#   - MT_high_QC_Hepatocyte is retained in the RDS for transparency but is
#     excluded from PRIMARY biological-state denominator and primary
#     mechanistic interpretation.
#   - With n=2 Sham and n=2 Tx, pseudobulk DE is exploratory. Effect direction,
#     replicate consistency, and biological coherence are prioritized.
#   - No cell-level p-value is used as the primary biological inference.
#
# PRIMARY QUESTIONS:
#   1) Which hepatocyte states change in abundance with transplantation?
#   2) Within Injury/inflammatory hepatocytes, do macrophage-recruiting /
#      inflammatory output genes decrease in Tx?
#   3) Are ER/UPR, oxidative, inflammatory, acute-phase, and regenerative
#      programs altered within comparable hepatocyte states?
#   4) Does hepatocyte biology support an upstream mechanism for the observed
#      macrophage shift and fibrosis suppression?
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
    "RNA" %in% Assays(object) &&
    length(Layers(object[["RNA"]])) > 1
  ) {
    object[["RNA"]] <- JoinLayers(
      object[["RNA"]]
    )
  }
  object
}

present_genes <- function(object, genes) {
  intersect(
    unique(genes),
    rownames(object)
  )
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
    consistency >= 0.75
  ) {
    "Moderate_Tx_up"
  } else if (
    direction ==
      "Tx_down" &&
    consistency >= 0.75
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
    select(
      sample,
      condition,
      value =
        all_of(
          value_col
        )
    ) %>%
    mutate(
      sample =
        as.character(
          sample
        ),
      condition =
        as.character(
          condition
        )
    )

  sample_values <- setNames(
    dat$value,
    dat$sample
  )

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

      cells <- which(
        groups ==
          g
      )

      if (
        length(
          cells
        ) ==
          1
      ) {
        as.numeric(
          counts[
            ,
            cells
          ]
        )
      } else {
        Matrix::rowSums(
          counts[
            ,
            cells,
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
  "Hepatocyte_cleanup_lineage_audit_v6.5.1",
  "RDS",
  "Mouse_MASH_Hepatocyte_highconf_clean_reclustered_v6.5.1.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_FINAL_state_pseudobulk_v6.5.2"
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
# 3. Settings
# ==============================================================================

CLUSTER_COL <-
  "hepclean_rpca_res_0.4"

SAMPLE_COL <-
  "sample_hep_v650"

CONDITION_COL <-
  "condition_hep_v650"

UMAP_NAME <-
  "umap.hep.clean.rpca"

FINAL_STATE_COL <-
  "hepatocyte_state_FINAL_v652"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

STATE_MAP <- c(
  "0" =
    "Periportal_Hepatocyte_1",
  "1" =
    "Injury_inflammatory_Hepatocyte",
  "2" =
    "Pericentral_Hepatocyte",
  "3" =
    "MT_high_QC_Hepatocyte",
  "4" =
    "Periportal_Hepatocyte_2",
  "5" =
    "Intermediate_Hepatocyte",
  "6" =
    "Cycling_G2M_Hepatocyte",
  "7" =
    "Cycling_S_Hepatocyte"
)

STATE_LEVELS <- unname(
  STATE_MAP
)

PRIMARY_BIO_STATES <- setdiff(
  STATE_LEVELS,
  "MT_high_QC_Hepatocyte"
)

MIN_CELLS_PB <-
  10

FDR_REPORT <-
  0.10


# ==============================================================================
# 4. A-priori program definitions
# ==============================================================================

PROGRAMS <- list(

  Hepatocyte_identity = c(
    "Alb",
    "Ttr",
    "Apoa1",
    "Apoa2",
    "Hnf4a",
    "Cps1",
    "Ass1"
  ),

  Periportal_Z1 = c(
    "Cps1",
    "Ass1",
    "Asl",
    "Arg1",
    "Pck1",
    "G6pc",
    "Gls2",
    "Hal"
  ),

  Pericentral_Z3 = c(
    "Glul",
    "Cyp2e1",
    "Cyp1a2",
    "Oat",
    "Axin2",
    "Slc1a2",
    "Cyp2a5"
  ),

  Injury_inflammatory = c(
    "Areg",
    "Ccl2",
    "Cxcl10",
    "Rsad2",
    "Lif",
    "Nupr1",
    "Sprr1a",
    "Socs3",
    "Nfkbia"
  ),

  Myeloid_recruitment_output = c(
    "Ccl2",
    "Cxcl1",
    "Cxcl10",
    "Csf1",
    "Il6",
    "Icam1",
    "Vcam1"
  ),

  Acute_phase = c(
    "Saa1",
    "Saa2",
    "Orm1",
    "Orm2",
    "Lcn2",
    "Serpina3n"
  ),

  ER_UPR_stress = c(
    "Ddit3",
    "Atf3",
    "Atf4",
    "Hspa5",
    "Xbp1",
    "Hsp90b1",
    "Herpud1",
    "Trib3"
  ),

  Oxidative_stress = c(
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Txnrd1",
    "Srxn1"
  ),

  Apoptosis_damage_response = c(
    "Bax",
    "Bbc3",
    "Pmaip1",
    "Cdkn1a",
    "Gadd45a",
    "Gadd45b",
    "Trp53inp1"
  ),

  Regenerative_cycling = c(
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
)

FOCUSED_GENES <- unique(
  c(
    # inflammatory / macrophage recruitment
    "Ccl2",
    "Cxcl1",
    "Cxcl10",
    "Csf1",
    "Il6",
    "Icam1",
    "Vcam1",

    # injury / acute phase
    "Areg",
    "Lif",
    "Nupr1",
    "Sprr1a",
    "Rsad2",
    "Socs3",
    "Nfkbia",
    "Saa1",
    "Saa2",
    "Orm1",
    "Orm2",
    "Lcn2",

    # ER / oxidative
    "Ddit3",
    "Atf3",
    "Atf4",
    "Hspa5",
    "Xbp1",
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",

    # death / damage
    "Bax",
    "Bbc3",
    "Pmaip1",
    "Cdkn1a",
    "Gadd45a",

    # cycling
    "Mki67",
    "Top2a",
    "Pcna",
    "Mcm5",
    "Ccna2",
    "Ccnb1",
    "Cdk1"
  )
)


# ==============================================================================
# 5. Load clean v6.5.1 object
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input v6.5.1 RDS missing: ",
    INPUT_RDS
  )
}

msg(
  "Loading v6.5.1 clean Hepatocyte RDS..."
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
  CLUSTER_COL,
  SAMPLE_COL,
  CONDITION_COL
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

clusters_present <- sort(
  unique(
    as.character(
      hep@meta.data[[
        CLUSTER_COL
      ]]
    )
  )
)

if (
  !setequal(
    clusters_present,
    names(
      STATE_MAP
    )
  )
) {
  stop(
    "Unexpected clean res0.4 clusters. Present: ",
    paste(
      clusters_present,
      collapse = ", "
    ),
    " ; expected: ",
    paste(
      names(
        STATE_MAP
      ),
      collapse = ", "
    )
  )
}


# ==============================================================================
# 6. Freeze final state annotation
# ==============================================================================

cluster_chr <- as.character(
  hep@meta.data[[
    CLUSTER_COL
  ]]
)

hep[[
  FINAL_STATE_COL
]] <- factor(
  unname(
    STATE_MAP[
      cluster_chr
    ]
  ),
  levels =
    STATE_LEVELS
)

hep$hepatocyte_primary_biology_v652 <-
  as.character(
    hep@meta.data[[
      FINAL_STATE_COL
    ]]
  ) %in%
  PRIMARY_BIO_STATES

annotation_manifest <- tibble(
  clean_cluster =
    names(
      STATE_MAP
    ),
  final_state =
    unname(
      STATE_MAP
    ),
  primary_biology =
    unname(
      STATE_MAP
    ) %in%
      PRIMARY_BIO_STATES,
  interpretation = c(
    "Periportal-like hepatocyte; Cyp2f2/Sds/Serpina1a/Hao2",
    "Areg/Ccl2/Sprr1a/Rsad2/Plat-enriched hepatocyte",
    "Pericentral-like hepatocyte; Glul/Slc1a2/Slc22a3",
    "High mitochondrial RNA / lower-complexity QC-associated hepatocyte; retained for transparency",
    "Second periportal-like metabolic hepatocyte; Gls2/Hal-enriched",
    "Intermediate zonation hepatocyte; no strong pathological label assigned",
    "Cycling hepatocyte enriched for Mki67/Cenpf/Ckap2/Gtse1",
    "Cycling hepatocyte enriched for Dtl/Uhrf1/Mcm5/Pole/Rad51"
  )
)

write.csv(
  annotation_manifest,
  file.path(
    TAB_OUT,
    "01_FINAL_Hepatocyte_annotation_manifest_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Source -> clean provenance table
# ==============================================================================

if (
  "source_cluster_v650" %in%
    colnames(
      hep@meta.data
    )
) {

  provenance <- hep@meta.data %>%
    as_tibble(
      rownames = "cell"
    ) %>%
    transmute(
      source_cluster =
        as.character(
          source_cluster_v650
        ),
      clean_cluster =
        as.character(
          .data[[
            CLUSTER_COL
          ]]
        ),
      final_state =
        as.character(
          .data[[
            FINAL_STATE_COL
          ]]
        )
    ) %>%
    count(
      source_cluster,
      clean_cluster,
      final_state,
      name =
        "n_cells"
    ) %>%
    group_by(
      clean_cluster
    ) %>%
    mutate(
      fraction_within_clean_cluster =
        n_cells /
          sum(
            n_cells
          )
    ) %>%
    ungroup() %>%
    arrange(
      as.numeric(
        clean_cluster
      ),
      desc(
        n_cells
      )
    )

  write.csv(
    provenance,
    file.path(
      TAB_OUT,
      "02_source_to_clean_cluster_provenance_v6.5.2.csv"
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 8. Final state composition by sample
# ==============================================================================

composition_all <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    cell,
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
          FINAL_STATE_COL
        ]]
      )
  ) %>%
  count(
    sample,
    condition,
    state,
    name =
      "n_cells"
  ) %>%
  complete(
    sample =
      SAMPLES,
    state =
      STATE_LEVELS,
    fill =
      list(
        n_cells = 0
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
      )
  ) %>%
  group_by(
    sample
  ) %>%
  mutate(
    denominator_all_clean =
      sum(
        n_cells
      ),
    fraction_all_clean =
      n_cells /
        denominator_all_clean
  ) %>%
  ungroup()

primary_denominator <- composition_all %>%
  filter(
    state %in%
      PRIMARY_BIO_STATES
  ) %>%
  group_by(
    sample
  ) %>%
  summarise(
    denominator_primary_biology =
      sum(
        n_cells
      ),
    .groups = "drop"
  )

composition <- composition_all %>%
  left_join(
    primary_denominator,
    by =
      "sample"
  ) %>%
  mutate(
    primary_biology =
      state %in%
        PRIMARY_BIO_STATES,
    fraction_primary_biology =
      ifelse(
        primary_biology,
        n_cells /
          denominator_primary_biology,
        NA_real_
      )
  )

write.csv(
  composition,
  file.path(
    TAB_OUT,
    "03_FINAL_Hepatocyte_state_composition_by_sample_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Replicate-aware state composition summary
# ==============================================================================

composition_summaries <- list()

for (
  st in STATE_LEVELS
) {

  dat <- composition %>%
    filter(
      state ==
        st
    )

  s_all <- four_sample_summary(
    dat %>%
      rename(
        value =
          fraction_all_clean
      ),
    "value",
    metric_name =
      "fraction_all_clean"
  ) %>%
    mutate(
      state =
        st,
      denominator =
        "all_clean_hepatocytes"
    )

  composition_summaries[[
    paste0(
      st,
      "_all"
    )
  ]] <- s_all

  if (
    st %in%
      PRIMARY_BIO_STATES
  ) {

    s_primary <- four_sample_summary(
      dat %>%
        rename(
          value =
            fraction_primary_biology
        ),
      "value",
      metric_name =
        "fraction_primary_biology"
    ) %>%
      mutate(
        state =
          st,
        denominator =
          "primary_biology_excluding_MT_high_QC"
      )

    composition_summaries[[
      paste0(
        st,
        "_primary"
      )
    ]] <- s_primary
  }
}

composition_summary <- bind_rows(
  composition_summaries
) %>%
  select(
    state,
    denominator,
    everything()
  )

write.csv(
  composition_summary,
  file.path(
    TAB_OUT,
    "04_FINAL_Hepatocyte_state_composition_replicate_summary_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Combined cycling burden
# ==============================================================================

cycling_states <- c(
  "Cycling_G2M_Hepatocyte",
  "Cycling_S_Hepatocyte"
)

cycling_by_sample <- composition %>%
  filter(
    state %in%
      cycling_states
  ) %>%
  group_by(
    sample,
    condition
  ) %>%
  summarise(
    n_cycling =
      sum(
        n_cells
      ),
    denominator_all_clean =
      first(
        denominator_all_clean
      ),
    denominator_primary_biology =
      first(
        denominator_primary_biology
      ),
    cycling_fraction_all_clean =
      n_cycling /
        denominator_all_clean,
    cycling_fraction_primary_biology =
      n_cycling /
        denominator_primary_biology,
    .groups = "drop"
  )

write.csv(
  cycling_by_sample,
  file.path(
    TAB_OUT,
    "05_combined_cycling_Hepatocyte_fraction_by_sample_v6.5.2.csv"
  ),
  row.names = FALSE
)

cycling_summary <- bind_rows(
  four_sample_summary(
    cycling_by_sample,
    "cycling_fraction_all_clean",
    "combined_cycling_fraction_all_clean"
  ),
  four_sample_summary(
    cycling_by_sample,
    "cycling_fraction_primary_biology",
    "combined_cycling_fraction_primary_biology"
  )
)

write.csv(
  cycling_summary,
  file.path(
    TAB_OUT,
    "06_combined_cycling_Hepatocyte_replicate_summary_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Build state x sample pseudobulk counts
# ==============================================================================

counts <- GetAssayData(
  hep,
  assay =
    "RNA",
  layer =
    "counts"
)

meta_pb <- hep@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  transmute(
    cell,
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
          FINAL_STATE_COL
        ]]
      )
  ) %>%
  mutate(
    pb_group =
      paste(
        state,
        sample,
        sep =
          "|||"
      )
  )

stopifnot(
  all(
    colnames(
      counts
    ) ==
      meta_pb$cell
  )
)

pb_counts_all <- aggregate_counts_by_group(
  counts,
  meta_pb$pb_group
)

pb_meta_all <- tibble(
  pb_group =
    colnames(
      pb_counts_all
    )
) %>%
  separate(
    pb_group,
    into = c(
      "state",
      "sample"
    ),
    sep =
      "\\|\\|\\|",
    remove =
      FALSE
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
      )
  )

cell_n_pb <- meta_pb %>%
  count(
    pb_group,
    name =
      "n_cells"
  )

pb_meta_all <- pb_meta_all %>%
  left_join(
    cell_n_pb,
    by =
      "pb_group"
  )

write.csv(
  pb_meta_all,
  file.path(
    TAB_OUT,
    "07_pseudobulk_group_cell_counts_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. edgeR state-wise pseudobulk DE
# ==============================================================================

de_results <- list()
logcpm_results <- list()

for (
  st in STATE_LEVELS
) {

  msg(
    "Pseudobulk DE: ",
    st
  )

  groups_needed <- paste(
    st,
    SAMPLES,
    sep =
      "|||"
  )

  meta_st <- pb_meta_all %>%
    filter(
      pb_group %in%
        groups_needed
    ) %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        ),
      condition =
        factor(
          condition,
          levels = c(
            "Sham",
            "Tx"
          )
        )
    ) %>%
    arrange(
      sample
    )

  if (
    nrow(
      meta_st
    ) !=
      4
  ) {
    warning(
      "Skipping ",
      st,
      ": not all 4 sample pseudobulks are present."
    )
    next
  }

  if (
    any(
      meta_st$n_cells <
        MIN_CELLS_PB
    )
  ) {
    warning(
      "State ",
      st,
      " has a sample with fewer than ",
      MIN_CELLS_PB,
      " cells. DE retained but flagged."
    )
  }

  mat <- pb_counts_all[
    ,
    meta_st$pb_group,
    drop =
      FALSE
  ]

  colnames(
    mat
  ) <- as.character(
    meta_st$sample
  )

  y <- edgeR::DGEList(
    counts =
      mat
  )

  keep <- edgeR::filterByExpr(
    y,
    group =
      meta_st$condition,
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
      meta_st
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
      state =
        st,
      comparison =
        "Tx_vs_Sham",
      n_Shams =
        2L,
      n_Tx =
        2L,
      min_cells_across_samples =
        min(
          meta_st$n_cells
        ),
      interpretation =
        "Exploratory sample-level pseudobulk; n=2 vs n=2"
    )

  de_results[[
    st
  ]] <- tab

  logcpm <- edgeR::cpm(
    y,
    log =
      TRUE,
    prior.count =
      2
  )

  logcpm_df <- as.data.frame(
    logcpm
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
      state =
        st,
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

  logcpm_results[[
    st
  ]] <- logcpm_df
}

de_all <- bind_rows(
  de_results
)

logcpm_all <- bind_rows(
  logcpm_results
)

write.csv(
  de_all,
  file.path(
    TAB_OUT,
    "08_statewise_pseudobulk_DE_Tx_vs_Sham_all_genes_v6.5.2.csv"
  ),
  row.names = FALSE
)

write.csv(
  logcpm_all,
  file.path(
    TAB_OUT,
    "09_statewise_pseudobulk_logCPM_all_genes_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Focused gene pseudobulk table
# ==============================================================================

focused_present <- intersect(
  FOCUSED_GENES,
  unique(
    logcpm_all$gene
  )
)

focused_logcpm <- logcpm_all %>%
  filter(
    gene %in%
      focused_present
  )

focused_wide <- focused_logcpm %>%
  select(
    state,
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

focused_summary_list <- list()

for (
  i in seq_len(
    nrow(
      focused_wide
    )
  )
) {

  row <- focused_wide[
    i,
    ,
    drop =
      FALSE
  ]

  sham <- c(
    row$Sham1,
    row$Sham20
  )

  tx <- c(
    row$Tx17,
    row$Tx5
  )

  d <- pairwise_direction_summary(
    sham,
    tx
  )

  focused_summary_list[[
    i
  ]] <- bind_cols(
    tibble(
      state =
        row$state,
      gene =
        row$gene,
      Sham1 =
        row$Sham1,
      Sham20 =
        row$Sham20,
      Tx17 =
        row$Tx17,
      Tx5 =
        row$Tx5,
      Sham_mean =
        mean(
          sham,
          na.rm =
            TRUE
        ),
      Tx_mean =
        mean(
          tx,
          na.rm =
            TRUE
        ),
      Tx_minus_Sham =
        mean(
          tx,
          na.rm =
            TRUE
        ) -
          mean(
            sham,
            na.rm =
              TRUE
          )
    ),
    d
  )
}

focused_summary <- bind_rows(
  focused_summary_list
)

focused_de <- de_all %>%
  filter(
    gene %in%
      FOCUSED_GENES
  ) %>%
  select(
    state,
    gene,
    logFC,
    logCPM,
    PValue,
    FDR
  )

focused_summary <- focused_summary %>%
  left_join(
    focused_de,
    by = c(
      "state",
      "gene"
    )
  ) %>%
  arrange(
    state,
    gene
  )

write.csv(
  focused_summary,
  file.path(
    TAB_OUT,
    "10_focused_Hepatocyte_signals_pseudobulk_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Program scores from pseudobulk logCPM
# ==============================================================================

program_gene_audit <- list()
program_sample_rows <- list()

for (
  st in STATE_LEVELS
) {

  state_log <- logcpm_all %>%
    filter(
      state ==
        st
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

  if (
    !all(
      SAMPLES %in%
        colnames(
          state_log
        )
    )
  ) {
    next
  }

  mat <- as.matrix(
    state_log[
      ,
      SAMPLES,
      drop =
        FALSE
    ]
  )

  rownames(
    mat
  ) <- state_log$gene

  for (
    prog in names(
      PROGRAMS
    )
  ) {

    genes <- intersect(
      PROGRAMS[[
        prog
      ]],
      rownames(
        mat
      )
    )

    program_gene_audit[[
      paste(
        st,
        prog,
        sep =
          "|||"
      )
    ]] <- tibble(
      state =
        st,
      program =
        prog,
      gene =
        PROGRAMS[[
          prog
        ]],
      present_in_filtered_pseudobulk =
        PROGRAMS[[
          prog
        ]] %in%
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

    sample_score <- colMeans(
      zmat,
      na.rm =
        TRUE
    )

    program_sample_rows[[
      paste(
        st,
        prog,
        sep =
          "|||"
      )
    ]] <- tibble(
      state =
        st,
      program =
        prog,
      sample =
        names(
          sample_score
        ),
      condition =
        ifelse(
          grepl(
            "^Tx",
            names(
              sample_score
            )
          ),
          "Tx",
          "Sham"
        ),
      score =
        as.numeric(
          sample_score
        ),
      n_genes =
        length(
          genes
        )
    )
  }
}

program_gene_audit <- bind_rows(
  program_gene_audit
)

program_sample <- bind_rows(
  program_sample_rows
)

write.csv(
  program_gene_audit,
  file.path(
    TAB_OUT,
    "11_pseudobulk_program_gene_audit_v6.5.2.csv"
  ),
  row.names = FALSE
)

write.csv(
  program_sample,
  file.path(
    TAB_OUT,
    "12_pseudobulk_program_scores_by_state_sample_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Replicate-aware pseudobulk program summary
# ==============================================================================

program_summary_rows <- list()

keys <- program_sample %>%
  distinct(
    state,
    program
  )

for (
  i in seq_len(
    nrow(
      keys
    )
  )
) {

  st <- keys$state[
    i
  ]

  prog <- keys$program[
    i
  ]

  dat <- program_sample %>%
    filter(
      state ==
        st,
      program ==
        prog
    )

  s <- four_sample_summary(
    dat,
    "score",
    prog
  ) %>%
    mutate(
      state =
        st,
      program =
        prog,
      n_genes =
        first(
          dat$n_genes
        )
    )

  program_summary_rows[[
    paste(
      st,
      prog,
      sep =
        "|||"
    )
  ]] <- s
}

program_summary <- bind_rows(
  program_summary_rows
) %>%
  select(
    state,
    program,
    n_genes,
    everything()
  )

write.csv(
  program_summary,
  file.path(
    TAB_OUT,
    "13_pseudobulk_program_replicate_summary_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Injury/inflammatory state focused interpretation table
# ==============================================================================

injury_state <-
  "Injury_inflammatory_Hepatocyte"

injury_gene_table <- focused_summary %>%
  filter(
    state ==
      injury_state
  ) %>%
  arrange(
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )

write.csv(
  injury_gene_table,
  file.path(
    TAB_OUT,
    "14_Injury_inflammatory_Hepatocyte_focused_genes_v6.5.2.csv"
  ),
  row.names = FALSE
)

injury_program_table <- program_summary %>%
  filter(
    state ==
      injury_state
  ) %>%
  arrange(
    program
  )

write.csv(
  injury_program_table,
  file.path(
    TAB_OUT,
    "15_Injury_inflammatory_Hepatocyte_programs_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Exploratory DE shortlist
# ==============================================================================

de_shortlist <- de_all %>%
  filter(
    state %in%
      PRIMARY_BIO_STATES,
    FDR <=
      FDR_REPORT
  ) %>%
  mutate(
    abs_logFC =
      abs(
        logFC
      ),
    direction =
      ifelse(
        logFC >
          0,
        "Tx_up",
        "Tx_down"
      )
  ) %>%
  arrange(
    state,
    FDR,
    desc(
      abs_logFC
    )
  )

write.csv(
  de_shortlist,
  file.path(
    TAB_OUT,
    "16_exploratory_statewise_DE_FDR0.10_v6.5.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Figures: final annotated UMAP
# ==============================================================================

p_umap <- DimPlot(
  hep,
  reduction =
    UMAP_NAME,
  group.by =
    FINAL_STATE_COL,
  label =
    TRUE,
  repel =
    TRUE,
  pt.size =
    0.28,
  raster =
    FALSE
) +
  ggtitle(
    "Final clean Hepatocyte states | v6.5.2"
  )

save_pdf(
  p_umap,
  file.path(
    FIG_OUT,
    "01_FINAL_Hepatocyte_state_UMAP_v6.5.2.pdf"
  ),
  9,
  7
)

p_umap_sample <- DimPlot(
  hep,
  reduction =
    UMAP_NAME,
  group.by =
    FINAL_STATE_COL,
  split.by =
    SAMPLE_COL,
  ncol =
    2,
  pt.size =
    0.24,
  raster =
    FALSE
) +
  plot_annotation(
    title =
      "Final Hepatocyte states by biological sample | v6.5.2"
  )

save_pdf(
  p_umap_sample,
  file.path(
    FIG_OUT,
    "02_FINAL_Hepatocyte_state_UMAP_by_sample_v6.5.2.pdf"
  ),
  14,
  11
)


# ==============================================================================
# 19. Figures: sample-level state fractions
# ==============================================================================

plot_comp <- composition %>%
  filter(
    state %in%
      PRIMARY_BIO_STATES
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    state =
      factor(
        state,
        levels =
          PRIMARY_BIO_STATES
      )
  )

p_comp <- ggplot(
  plot_comp,
  aes(
    x =
      sample,
    y =
      fraction_primary_biology,
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
      2
  ) +
  facet_wrap(
    ~ state,
    scales =
      "free_y",
    ncol =
      3
  ) +
  labs(
    title =
      "Final Hepatocyte state fractions",
    subtitle =
      "Denominator excludes MT-high/QC hepatocytes",
    x = NULL,
    y =
      "Fraction of primary biological hepatocytes"
  ) +
  theme_classic(
    base_size =
      8
  )

save_pdf(
  p_comp,
  file.path(
    FIG_OUT,
    "03_FINAL_Hepatocyte_state_fractions_by_sample_v6.5.2.pdf"
  ),
  13,
  9
)


# ==============================================================================
# 20. Figures: focused injury-state genes
# ==============================================================================

injury_plot_genes <- intersect(
  c(
    "Ccl2",
    "Cxcl1",
    "Cxcl10",
    "Csf1",
    "Il6",
    "Areg",
    "Lif",
    "Nupr1",
    "Sprr1a",
    "Rsad2",
    "Socs3",
    "Saa1",
    "Saa2",
    "Ddit3",
    "Atf3",
    "Hspa5",
    "Hmox1",
    "Nqo1"
  ),
  injury_gene_table$gene
)

injury_plot_df <- injury_gene_table %>%
  filter(
    gene %in%
      injury_plot_genes
  ) %>%
  select(
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
    injury_plot_df
  )
) {

  p_injury_genes <- ggplot(
    injury_plot_df,
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
    facet_wrap(
      ~ gene,
      scales =
        "free_y",
      ncol =
        4
    ) +
    labs(
      title =
        "Injury/inflammatory Hepatocyte | focused pseudobulk genes",
      subtitle =
        "Biological replicate values; Sham1, Sham20, Tx17, Tx5",
      x = NULL,
      y =
        "TMM-normalized logCPM"
    ) +
    theme_classic(
      base_size =
        8
    )

  save_pdf(
    p_injury_genes,
    file.path(
      FIG_OUT,
      "04_Injury_inflammatory_Hepatocyte_focused_genes_v6.5.2.pdf"
    ),
    14,
    10
  )
}


# ==============================================================================
# 21. Figures: injury-state program panel
# ==============================================================================

injury_program_plot <- program_sample %>%
  filter(
    state ==
      injury_state
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
    injury_program_plot
  )
) {

  p_injury_program <- ggplot(
    injury_program_plot,
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
        "Injury/inflammatory Hepatocyte | pseudobulk programs",
      subtitle =
        "Gene-wise z-scores averaged within each a-priori program",
      x = NULL,
      y =
        "Program score"
    ) +
    theme_classic(
      base_size =
        8
    )

  save_pdf(
    p_injury_program,
    file.path(
      FIG_OUT,
      "05_Injury_inflammatory_Hepatocyte_programs_v6.5.2.pdf"
    ),
    13,
    10
  )
}


# ==============================================================================
# 22. Figures: all-state program heatmap
# ==============================================================================

program_heat <- program_summary %>%
  select(
    state,
    program,
    Tx_minus_Sham,
    evidence_grade
  ) %>%
  mutate(
    state =
      factor(
        state,
        levels =
          STATE_LEVELS
      )
  )

p_program_heat <- ggplot(
  program_heat,
  aes(
    x =
      program,
    y =
      state,
    fill =
      Tx_minus_Sham
  )
) +
  geom_tile(
    linewidth =
      0.25
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
      "Hepatocyte pseudobulk program changes | Tx - Sham",
    subtitle =
      "Sample-level effect direction; n=2 vs n=2",
    x = NULL,
    y = NULL,
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
  p_program_heat,
  file.path(
    FIG_OUT,
    "06_Hepatocyte_program_Tx_minus_Sham_heatmap_v6.5.2.pdf"
  ),
  11,
  7
)


# ==============================================================================
# 23. Save final annotated RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_Hepatocyte_FINAL_annotated_v6.5.2.rds"
)

saveRDS(
  hep,
  RDS_FILE,
  compress =
    FALSE
)

msg(
  "Saved final annotated Hepatocyte RDS: ",
  RDS_FILE
)


# ==============================================================================
# 24. Analysis manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "final_annotation_column",
    "working_cluster_column",
    "final_states",
    "MT_high_QC_primary_exclusion",
    "pseudobulk_method",
    "comparison",
    "biological_replicates",
    "minimum_cells_per_state_sample",
    "formal_inference_note"
  ),
  value = c(
    "v6.5.2",
    INPUT_RDS,
    FINAL_STATE_COL,
    CLUSTER_COL,
    paste(
      STATE_LEVELS,
      collapse =
        ","
    ),
    "TRUE; retained in object but excluded from primary biological-state denominator",
    "sample-by-state raw-count aggregation + edgeR TMM QL",
    "Tx vs Sham",
    "Sham1,Sham20 vs Tx17,Tx5",
    as.character(
      MIN_CELLS_PB
    ),
    "n=2 vs n=2; exploratory. Prioritize effect direction and replicate consistency."
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.5.2.csv"
  ),
  row.names =
    FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.5.2.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
