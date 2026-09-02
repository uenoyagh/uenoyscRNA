#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6643)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Mphi5 -> Hepatocyte mechanistic closure correction
#
# Version: v6.6.4.4
#
# FIX FROM v6.6.4.3
# -------------------
# The Pdgfrb-positive / Pdgfrb-negative marker comparison contained an invalid
# chained comparison. v6.6.4.4 separates the Boolean positivity flag from the
# equality test:
#
#   pdgfrb_flag <- as.numeric(...) > 0
#   cells_group <- cells_state[pdgfrb_flag == pdgfrb_status]
#
# No biological design, thresholds, gene sets, or upstream results are changed.
#
#
# PURPOSE
# -------
# Correct and finalize v6.6.4.2 WITHOUT rerunning:
#   - CellChat
#   - NicheNet
#   - edgeR
#   - clustering
#
# v6.6.4.2 successfully produced sender and receptor validation, but all
# independent receiver-program values were NA because a 2D pseudobulk matrix
# was indexed with one-dimensional character indexing:
#
#     pb_counts[genes_use]
#
# v6.6.4.4 uses a named pseudobulk COUNT VECTOR and indexes it directly:
#
#     pb_counts[genes_use]
#
# This preserves gene names and produces valid CP10k program values.
#
# v6.6.4.4 also tightens interpretation of Pdgfrb-positive Hepatocytes.
# The previous heuristic only labeled a cell suspicious if it had:
#   >=2 mesenchymal genes AND <=2 hepatocyte-identity genes.
#
# Strong hepatocyte identity does not exclude double-positive / mixed
# transcriptional profiles. Therefore this correction explicitly compares:
#   Pdgfrb-positive vs Pdgfrb-negative Hepatocytes
# for hepatocyte and HSC/mesenchymal marker burden.
#
# IMPORTANT BIOLOGICAL DISTINCTION
# --------------------------------
# Sender evidence is retained at three levels:
#
# 1) per-cell ligand expression
# 2) Repair/Resolution-Mphi population-weighted output
# 3) total Mphi5 population-weighted output
#
# These are NOT interchangeable.
#
# CellChat v6.6.1/6.6.2.1 used population.size = FALSE and therefore most
# closely reflects per-cell state-dependent communication rather than
# abundance-weighted ligand mass.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
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


sample_condition <- function(sample_name) {

  ifelse(
    grepl(
      "^Sham",
      sample_name
    ),
    "Sham",
    ifelse(
      grepl(
        "^Tx",
        sample_name
      ),
      "Tx",
      NA_character_
    )
  )
}


pairwise_direction_metrics <- function(
  sham1,
  sham20,
  tx17,
  tx5
) {

  vals <- c(
    sham1,
    sham20,
    tx17,
    tx5
  )

  if (
    any(
      !is.finite(
        vals
      )
    )
  ) {

    return(
      tibble(
        pairwise_up_n =
          NA_integer_,
        pairwise_down_n =
          NA_integer_,
        pairwise_tie_n =
          NA_integer_,
        pairwise_direction =
          NA_character_,
        pairwise_direction_consistency =
          NA_real_,
        both_Tx_above_both_Sham =
          NA,
        both_Tx_below_both_Sham =
          NA,
        evidence_grade =
          "Missing"
      )
    )
  }

  diffs <- c(
    tx17 - sham1,
    tx17 - sham20,
    tx5 - sham1,
    tx5 - sham20
  )

  up_n <- sum(
    diffs > 0
  )

  down_n <- sum(
    diffs < 0
  )

  tie_n <- sum(
    diffs == 0
  )

  up_fraction <- up_n / 4
  down_fraction <- down_n / 4

  direction <- case_when(
    up_n > down_n ~
      "Tx_up",

    down_n > up_n ~
      "Tx_down",

    TRUE ~
      "Mixed_or_tied"
  )

  consistency <- max(
    up_fraction,
    down_fraction
  )

  both_tx_above <- (
    min(
      tx17,
      tx5
    ) >
      max(
        sham1,
        sham20
      )
  )

  both_tx_below <- (
    max(
      tx17,
      tx5
    ) <
      min(
        sham1,
        sham20
      )
  )

  grade <- case_when(
    both_tx_above ~
      "Strong_Tx_up",

    both_tx_below ~
      "Strong_Tx_down",

    direction ==
      "Tx_up" &&
    consistency >=
      0.75 ~
      "Moderate_Tx_up",

    direction ==
      "Tx_down" &&
    consistency >=
      0.75 ~
      "Moderate_Tx_down",

    TRUE ~
      "Mixed_or_weak"
  )

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


summarize_four_samples <- function(
  df,
  key_cols,
  sample_col,
  value_col
) {

  samples <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  keys <- df %>%
    distinct(
      across(
        all_of(
          key_cols
        )
      )
    )

  grid <- tidyr::crossing(
    keys,
    sample =
      samples
  )

  vals <- df %>%
    transmute(
      across(
        all_of(
          key_cols
        )
      ),
      sample =
        as.character(
          .data[[
            sample_col
          ]]
        ),
      value =
        as.numeric(
          .data[[
            value_col
          ]]
        )
    )

  wide <- grid %>%
    left_join(
      vals,
      by = c(
        key_cols,
        "sample"
      )
    ) %>%
    pivot_wider(
      names_from =
        sample,
      values_from =
        value
    )

  metrics <- lapply(
    seq_len(
      nrow(
        wide
      )
    ),
    function(i) {

      pairwise_direction_metrics(
        wide$Sham1[[i]],
        wide$Sham20[[i]],
        wide$Tx17[[i]],
        wide$Tx5[[i]]
      )
    }
  ) %>%
    bind_rows()

  bind_cols(
    wide,
    metrics
  ) %>%
    mutate(
      Sham_mean =
        (
          Sham1 +
          Sham20
        ) /
        2,

      Tx_mean =
        (
          Tx17 +
          Tx5
        ) /
        2,

      Tx_minus_Sham =
        Tx_mean -
        Sham_mean
    )
}


pseudobulk_count_vector <- function(
  counts,
  cells
) {

  if (
    !length(
      cells
    )
  ) {
    stop(
      "pseudobulk_count_vector received zero cells."
    )
  }

  x <- Matrix::rowSums(
    counts[
      ,
      cells,
      drop = FALSE
    ]
  )

  names(x) <-
    rownames(
      counts
    )

  x
}


program_value_from_named_vector <- function(
  pb_counts,
  genes
) {

  genes_use <- intersect(
    genes,
    names(
      pb_counts
    )
  )

  if (
    !length(
      genes_use
    )
  ) {

    return(
      list(
        value =
          NA_real_,
        n_present =
          0L,
        present_genes =
          ""
      )
    )
  }

  lib <- sum(
    pb_counts,
    na.rm = TRUE
  )

  if (
    !is.finite(
      lib
    ) ||
    lib <= 0
  ) {

    return(
      list(
        value =
          NA_real_,
        n_present =
          length(
            genes_use
          ),
        present_genes =
          paste(
            genes_use,
            collapse = "|"
          )
      )
    )
  }

  cp10k <- as.numeric(
    pb_counts[
      genes_use
    ]
  ) /
    lib *
    10000

  value <- mean(
    log1p(
      cp10k
    ),
    na.rm = TRUE
  )

  list(
    value =
      value,
    n_present =
      length(
        genes_use
      ),
    present_genes =
      paste(
        genes_use,
        collapse = "|"
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
  "Mphi5_Hep5_interaction_ready_v6.6.0",
  "RDS",
  "Mouse_MASH_Mphi5_Hep5_interaction_ready_v6.6.0.rds"
)


V6642 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_mechanistic_closure_v6.6.4.2",
  "Tables"
)


OLD_CLOSURE <- file.path(
  V6642,
  "17_mechanistic_closure_summary_v6.6.4.2.csv"
)


OLD_PDGFRB <- file.path(
  V6642,
  "16_Pdgfrb_receiver_closure_v6.6.4.2.csv"
)


OLD_REPAIR_PERCELL <- file.path(
  V6642,
  "08_Mphi_sender_ligand_percell_summary_v6.6.4.2.csv"
)


OLD_REPAIR_WEIGHTED <- file.path(
  V6642,
  "10_RepairResolutionMphi_population_weighted_summary_v6.6.4.2.csv"
)


OLD_TOTAL_WEIGHTED <- file.path(
  V6642,
  "12_total_Mphi5_population_weighted_ligand_summary_v6.6.4.2.csv"
)


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_mechanistic_closure_v6.6.4.4"
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
  "interaction_celltype_v660"


LINEAGE_COL <-
  "interaction_lineage_v660"


SAMPLE_COL <-
  "sample_interaction_v660"


SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)


PRIMARY_HEP_STATES <- c(
  "Hep_Injury-inflammatory",
  "Hep_Pericentral"
)


HEP_IDENTITY_GENES <- c(
  "Alb",
  "Ttr",
  "Apoa1",
  "Apoa2",
  "Hnf4a",
  "Cps1",
  "Ass1"
)


HSC_MESENCHYMAL_GENES <- c(
  "Lrat",
  "Rbp1",
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Dcn",
  "Lum",
  "Pdgfra",
  "Rgs5",
  "Des"
)


PROGRAMS <- list(

  PDGF_IE_MAPK_proxy = c(
    "Fos",
    "Fosb",
    "Jun",
    "Junb",
    "Egr1",
    "Egr2",
    "Dusp1",
    "Dusp5",
    "Spry2",
    "Myc"
  ),

  Proliferation_cell_cycle = c(
    "Mki67",
    "Top2a",
    "Ccna2",
    "Ccnb1",
    "Cdk1",
    "Pcna",
    "Mcm2",
    "Mcm5"
  ),

  TNF_NFkB_response = c(
    "Nfkbia",
    "Tnfaip3",
    "Relb",
    "Icam1",
    "Vcam1",
    "Ccl2",
    "Cxcl10",
    "Cxcl1"
  ),

  Injury_stress = c(
    "Atf3",
    "Gadd45a",
    "Cdkn1a",
    "Ddit3",
    "Nupr1",
    "Hmox1",
    "Sqstm1",
    "Plin2"
  ),

  ECM_contact_remodeling = c(
    "Sdc1",
    "Sdc4",
    "Cd44",
    "Itga3",
    "Itgb1",
    "Itgav",
    "Itgb5",
    "Ptk2",
    "Vcl",
    "Pxn",
    "Ccn2"
  )
)


# ==============================================================================
# 4. Preflight and load
# ==============================================================================

for (
  f in c(
    INPUT_RDS,
    OLD_CLOSURE,
    OLD_PDGFRB,
    OLD_REPAIR_PERCELL,
    OLD_REPAIR_WEIGHTED,
    OLD_TOTAL_WEIGHTED
  )
) {

  if (
    !file.exists(f)
  ) {

    stop(
      "Required input missing: ",
      f
    )
  }
}


msg(
  "Loading v6.6.0 interaction-ready object..."
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


counts <- GetAssayData(
  obj,
  assay =
    "RNA",
  layer =
    "counts"
)


meta <- obj@meta.data


for (
  col in c(
    GROUP_COL,
    LINEAGE_COL,
    SAMPLE_COL
  )
) {

  if (
    !col %in%
      colnames(
        meta
      )
  ) {
    stop(
      "Missing metadata column: ",
      col
    )
  }
}


meta$closure_group <-
  as.character(
    meta[[
      GROUP_COL
    ]]
  )


meta$closure_lineage <-
  as.character(
    meta[[
      LINEAGE_COL
    ]]
  )


meta$closure_sample <-
  as.character(
    meta[[
      SAMPLE_COL
    ]]
  )


old_closure <- read.csv(
  OLD_CLOSURE,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


old_pdgfrb <- read.csv(
  OLD_PDGFRB,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


repair_percell <- read.csv(
  OLD_REPAIR_PERCELL,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


repair_weighted <- read.csv(
  OLD_REPAIR_WEIGHTED,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


total_weighted <- read.csv(
  OLD_TOTAL_WEIGHTED,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  as_tibble()


# ==============================================================================
# 5. Correct independent receiver programs
# ==============================================================================

msg(
  "Correcting receiver program calculation..."
)


program_rows <- list()


for (
  hep_state in PRIMARY_HEP_STATES
) {

  for (
    sample_name in SAMPLES
  ) {

    cells <- rownames(meta)[
      meta$closure_group ==
        hep_state &
      meta$closure_sample ==
        sample_name
    ]

    if (
      !length(
        cells
      )
    ) {
      stop(
        "No cells for ",
        hep_state,
        " / ",
        sample_name
      )
    }

    pb <- pseudobulk_count_vector(
      counts,
      cells
    )

    for (
      program_name in names(
        PROGRAMS
      )
    ) {

      ans <- program_value_from_named_vector(
        pb,
        PROGRAMS[[
          program_name
        ]]
      )

      program_rows[[
        length(
          program_rows
        ) + 1
      ]] <- tibble(
        sample =
          sample_name,
        condition =
          sample_condition(
            sample_name
          ),
        hepatocyte_state =
          hep_state,
        program =
          program_name,
        module_mean_log1p_CP10k =
          ans$value,
        n_present_genes =
          ans$n_present,
        n_requested_genes =
          length(
            PROGRAMS[[
              program_name
            ]]
          ),
        present_genes =
          ans$present_genes
      )
    }
  }
}


program_sample <- bind_rows(
  program_rows
)


if (
  any(
    !is.finite(
      program_sample$module_mean_log1p_CP10k
    )
  )
) {

  bad <- program_sample %>%
    filter(
      !is.finite(
        module_mean_log1p_CP10k
      )
    )

  print(
    bad
  )

  stop(
    "Corrected receiver program calculation still contains non-finite values."
  )
}


write.csv(
  program_sample,
  file.path(
    TAB_OUT,
    "01_Hep_receiver_programs_by_sample_CORRECTED_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


program_summary <- summarize_four_samples(
  program_sample,
  key_cols = c(
    "hepatocyte_state",
    "program"
  ),
  sample_col =
    "sample",
  value_col =
    "module_mean_log1p_CP10k"
)


write.csv(
  program_summary,
  file.path(
    TAB_OUT,
    "02_Hep_receiver_program_replicate_summary_CORRECTED_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. Program gene-availability audit
# ==============================================================================

program_gene_audit <- lapply(
  names(
    PROGRAMS
  ),
  function(program_name) {

    genes <- PROGRAMS[[
      program_name
    ]]

    tibble(
      program =
        program_name,
      gene =
        genes,
      present =
        genes %in%
          rownames(
            counts
          )
    )
  }
) %>%
  bind_rows()


write.csv(
  program_gene_audit,
  file.path(
    TAB_OUT,
    "03_receiver_program_gene_availability_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Stricter Pdgfrb-positive Hepatocyte audit
# ==============================================================================

msg(
  "Running stricter Pdgfrb-positive Hepatocyte audit..."
)


hep_cells <- rownames(meta)[
  meta$closure_lineage ==
    "Hepatocyte"
]


if (
  !"Pdgfrb" %in%
    rownames(
      counts
    )
) {

  stop(
    "Pdgfrb is absent from RNA counts."
  )
}


hep_identity_use <- intersect(
  HEP_IDENTITY_GENES,
  rownames(
    counts
  )
)


hsc_use <- intersect(
  HSC_MESENCHYMAL_GENES,
  rownames(
    counts
  )
)


hep_id_detect <- counts[
  hep_identity_use,
  hep_cells,
  drop = FALSE
] > 0


hsc_detect <- counts[
  hsc_use,
  hep_cells,
  drop = FALSE
] > 0


pdgfrb_positive <- as.numeric(
  counts[
    "Pdgfrb",
    hep_cells,
    drop = TRUE
  ]
) > 0


cell_audit <- tibble(
  cell =
    hep_cells,
  sample =
    meta[
      hep_cells,
      "closure_sample"
    ],
  hepatocyte_state =
    meta[
      hep_cells,
      "closure_group"
    ],
  Pdgfrb_positive =
    pdgfrb_positive,
  hepatocyte_identity_n =
    as.integer(
      Matrix::colSums(
        hep_id_detect
      )
    ),
  HSC_mesenchymal_n =
    as.integer(
      Matrix::colSums(
        hsc_detect
      )
    )
)


write.csv(
  cell_audit,
  file.path(
    TAB_OUT,
    "04_Pdgfrb_Hep_celllevel_identity_mesenchymal_burden_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


pdgfrb_posneg_summary <- cell_audit %>%
  filter(
    hepatocyte_state %in%
      PRIMARY_HEP_STATES
  ) %>%
  group_by(
    hepatocyte_state,
    Pdgfrb_positive
  ) %>%
  summarise(
    n_cells =
      n(),

    median_hepatocyte_identity_n =
      median(
        hepatocyte_identity_n
      ),

    mean_hepatocyte_identity_n =
      mean(
        hepatocyte_identity_n
      ),

    median_HSC_mesenchymal_n =
      median(
        HSC_mesenchymal_n
      ),

    mean_HSC_mesenchymal_n =
      mean(
        HSC_mesenchymal_n
      ),

    pct_HSC_mesenchymal_ge2 =
      mean(
        HSC_mesenchymal_n >=
          2
      ),

    pct_HSC_mesenchymal_ge3 =
      mean(
        HSC_mesenchymal_n >=
          3
      ),

    pct_HSC_mesenchymal_ge4 =
      mean(
        HSC_mesenchymal_n >=
          4
      ),

    pct_Hep_identity_ge5 =
      mean(
        hepatocyte_identity_n >=
          5
      ),

    pct_Hep_identity_ge6 =
      mean(
        hepatocyte_identity_n >=
          6
      ),

    .groups =
      "drop"
  )


write.csv(
  pdgfrb_posneg_summary,
  file.path(
    TAB_OUT,
    "05_Pdgfrb_positive_vs_negative_Hep_summary_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# Marker-level Pdgfrb+ vs Pdgfrb- co-expression.
marker_use <- unique(
  c(
    hep_identity_use,
    hsc_use
  )
)


pdgfrb_marker_rows <- list()


for (
  hep_state in PRIMARY_HEP_STATES
) {

  cells_state <- intersect(
    hep_cells,
    rownames(meta)[
      meta$closure_group ==
        hep_state
    ]
  )

  for (
    pdgfrb_status in c(
      FALSE,
      TRUE
    )
  ) {

    pdgfrb_flag <- as.numeric(
      counts[
        "Pdgfrb",
        cells_state,
        drop = TRUE
      ]
    ) > 0

    cells_group <- cells_state[
      pdgfrb_flag ==
        pdgfrb_status
    ]

    if (
      !length(
        cells_group
      )
    ) {
      next
    }

    det <- counts[
      marker_use,
      cells_group,
      drop = FALSE
    ] > 0

    pdgfrb_marker_rows[[
      length(
        pdgfrb_marker_rows
      ) + 1
    ]] <- tibble(
      hepatocyte_state =
        hep_state,
      Pdgfrb_positive =
        pdgfrb_status,
      gene =
        marker_use,
      n_cells =
        length(
          cells_group
        ),
      pct_expressed =
        as.numeric(
          Matrix::rowMeans(
            det
          )
        ),
      marker_class =
        ifelse(
          marker_use %in%
            hep_identity_use,
          "Hepatocyte_identity",
          "HSC_mesenchymal"
        )
    )
  }
}


pdgfrb_marker_comparison <- bind_rows(
  pdgfrb_marker_rows
)


write.csv(
  pdgfrb_marker_comparison,
  file.path(
    TAB_OUT,
    "06_Pdgfrb_positive_vs_negative_marker_coexpression_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Revised Pdgfrb receiver interpretation
# ==============================================================================

pdgfrb_context <- pdgfrb_posneg_summary %>%
  select(
    hepatocyte_state,
    Pdgfrb_positive,
    n_cells,
    median_HSC_mesenchymal_n,
    pct_HSC_mesenchymal_ge3,
    pct_Hep_identity_ge5
  ) %>%
  pivot_wider(
    names_from =
      Pdgfrb_positive,
    values_from = c(
      n_cells,
      median_HSC_mesenchymal_n,
      pct_HSC_mesenchymal_ge3,
      pct_Hep_identity_ge5
    ),
    names_glue =
      "{.value}_Pdgfrb_{Pdgfrb_positive}"
  )


pdgfrb_revised <- old_pdgfrb %>%
  left_join(
    pdgfrb_context,
    by =
      "hepatocyte_state"
  ) %>%
  mutate(
    revised_Pdgfrb_interpretation =
      case_when(

        receptor_detectability ==
          "Clearly_detectable" &
        !is.na(
          median_HSC_mesenchymal_n_Pdgfrb_TRUE
        ) &
        median_HSC_mesenchymal_n_Pdgfrb_TRUE >=
          3 ~
          "Pdgfrb_detectable_but_Pdgfrb_positive_Hep_show_mesenchymal_coexpression_interpret_cautiously",

        receptor_detectability ==
          "Clearly_detectable" &
        !is.na(
          pct_Hep_identity_ge5_Pdgfrb_TRUE
        ) &
        pct_Hep_identity_ge5_Pdgfrb_TRUE >=
          0.80 ~
          "Pdgfrb_detectable_in_strong_Hep_identity_cells",

        receptor_detectability ==
          "Low_but_detectable" ~
          "Pdgfrb_low_but_detectable",

        TRUE ~
          "Pdgfrb_sparse_or_uncertain"
      )
  )


write.csv(
  pdgfrb_revised,
  file.path(
    TAB_OUT,
    "07_Pdgfrb_receiver_closure_REVISED_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Reintegrate corrected receiver programs into closure summary
# ==============================================================================

program_lookup <- program_summary %>%
  transmute(
    receiver =
      hepatocyte_state,
    receiver_program_proxy =
      program,
    receiver_program_Sham_mean_CORRECTED =
      Sham_mean,
    receiver_program_Tx_mean_CORRECTED =
      Tx_mean,
    receiver_program_delta_CORRECTED =
      Tx_minus_Sham,
    receiver_program_pairwise_consistency_CORRECTED =
      pairwise_direction_consistency,
    receiver_program_grade_CORRECTED =
      evidence_grade
  )


closure_corrected <- old_closure %>%
  select(
    -any_of(
      c(
        "receiver_program_Sham_mean",
        "receiver_program_Tx_mean",
        "receiver_program_delta",
        "receiver_program_grade"
      )
    )
  ) %>%
  left_join(
    program_lookup,
    by = c(
      "receiver",
      "receiver_program_proxy"
    ),
    relationship =
      "many-to-one"
  )


# Add revised Pdgfrb context where relevant.
pdgfrb_lookup <- pdgfrb_revised %>%
  transmute(
    receiver =
      hepatocyte_state,
    receptor_audit_gene =
      receptor,
    revised_Pdgfrb_interpretation
  )


closure_corrected <- closure_corrected %>%
  left_join(
    pdgfrb_lookup,
    by = c(
      "receiver",
      "receptor_audit_gene"
    ),
    relationship =
      "many-to-one"
  ) %>%
  mutate(
    sender_level_note =
      case_when(
        ligand ==
          "Pdgfb" &
        sender_percell_delta <
          0 &
        Repair_weighted_delta >
          0 ~
          "Per_cell_down_but_Repair_population_weighted_up_due_to_subtype_expansion",

        ligand ==
          "Sema4d" &
        sender_percell_delta <
          0 &
        Repair_weighted_delta >
          0 ~
          "Per_cell_down_but_Repair_population_weighted_up_due_to_subtype_expansion",

        ligand ==
          "Fn1" &
        sender_percell_delta >
          0 &
        Repair_weighted_delta >
          0 ~
          "Per_cell_and_Repair_population_weighted_both_up",

        ligand ==
          "Plau" &
        sender_percell_delta <
          0 &
        Repair_weighted_delta <
          0 ~
          "Per_cell_and_Repair_population_weighted_both_down",

        ligand ==
          "Tnf" &
        sender_percell_delta <
          0 &
        Repair_weighted_delta <
          0 ~
          "Per_cell_and_Repair_population_weighted_both_down_or_mixed",

        TRUE ~
          "Mixed_sender_level_pattern"
      ),

    receiver_program_direction_match_to_CellChat =
      case_when(
        is.na(
          CellChat_best_supported_delta
        ) |
        is.na(
          receiver_program_delta_CORRECTED
        ) ~
          NA,

        CellChat_best_supported_delta >
          0 &
        receiver_program_delta_CORRECTED >
          0 ~
          TRUE,

        CellChat_best_supported_delta <
          0 &
        receiver_program_delta_CORRECTED <
          0 ~
          TRUE,

        TRUE ~
          FALSE
      ),

    closure_status_v6643 =
      case_when(
        ligand ==
          "Pdgfb" &
        sender_percell_delta <
          0 &
        !is.na(
          CellChat_best_supported_delta
        ) &
        CellChat_best_supported_delta <
          0 &
        receptor_detectability ==
          "Clearly_detectable" &
        receiver_program_delta_CORRECTED <
          0 ~
          "PDGFB_percell_CellChat_receptor_receiverProgram_concordant_Tx_down",

        ligand ==
          "Pdgfb" &
        sender_percell_delta <
          0 &
        !is.na(
          CellChat_best_supported_delta
        ) &
        CellChat_best_supported_delta <
          0 &
        receptor_detectability ==
          "Clearly_detectable" &
        receiver_program_delta_CORRECTED >=
          0 ~
          "PDGFB_sender_CellChat_down_but_receiverProgram_not_down",

        ligand ==
          "Sema4d" &
        sender_percell_delta <
          0 &
        !is.na(
          CellChat_best_supported_delta
        ) &
        CellChat_best_supported_delta <
          0 &
        receiver_program_delta_CORRECTED <
          0 ~
          "SEMA4D_sender_CellChat_receiverProgram_concordant_Tx_down",

        ligand ==
          "Fn1" &
        sender_percell_delta >
          0 &
        !is.na(
          CellChat_best_supported_delta
        ) &
        CellChat_best_supported_delta >
          0 &
        receiver_program_delta_CORRECTED >
          0 ~
          "FN1_sender_CellChat_receiverProgram_concordant_Tx_up",

        ligand ==
          "Plau" &
        sender_percell_delta <
          0 &
        !is.na(
          CellChat_best_supported_delta
        ) &
        CellChat_best_supported_delta <
          0 &
        receiver_program_delta_CORRECTED <
          0 ~
          "PLAU_sender_CellChat_receiverProgram_concordant_Tx_down",

        ligand ==
          "Tnf" &
        sender_percell_delta <
          0 &
        receiver_program_delta_CORRECTED <
          0 ~
          "TNF_sender_and_receiverProgram_down_CellChat_not_replicate_supported",

        TRUE ~
          "Partial_or_directionally_mixed"
      )
  )


write.csv(
  closure_corrected,
  file.path(
    TAB_OUT,
    "08_FINAL_mechanistic_closure_CORRECTED_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Compact key-axis table
# ==============================================================================

key_axis_compact <- closure_corrected %>%
  select(
    sender,
    receiver,
    receiver_program,
    ligand,
    receptor_audit_gene,
    sender_percell_Sham_mean,
    sender_percell_Tx_mean,
    sender_percell_delta,
    sender_percell_grade,
    Repair_weighted_Sham_mean,
    Repair_weighted_Tx_mean,
    Repair_weighted_delta,
    Repair_weighted_grade,
    total_Mphi5_weighted_Sham_mean,
    total_Mphi5_weighted_Tx_mean,
    total_Mphi5_weighted_delta,
    total_Mphi5_weighted_grade,
    CellChat_best_supported_delta,
    CellChat_support_samples,
    CellChat_pairwise_consistency,
    nichenet_rank,
    nichenet_percentile_corrected,
    evidence_class,
    receptor_Sham_mean_pct,
    receptor_Tx_mean_pct,
    receptor_detectability,
    revised_Pdgfrb_interpretation,
    receiver_program_proxy,
    receiver_program_Sham_mean_CORRECTED,
    receiver_program_Tx_mean_CORRECTED,
    receiver_program_delta_CORRECTED,
    receiver_program_grade_CORRECTED,
    sender_level_note,
    closure_status_v6643
  )


write.csv(
  key_axis_compact,
  file.path(
    TAB_OUT,
    "09_KEY_AXIS_mechanistic_closure_compact_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Repair-vs-total sender level audit
# ==============================================================================

sender_level_audit <- repair_percell %>%
  filter(
    sender ==
      "Mphi_Repair-Resolution",
    ligand %in%
      c(
        "Pdgfb",
        "Sema4d",
        "Plau",
        "Fn1",
        "Tnf"
      )
  ) %>%
  transmute(
    ligand,
    percell_Sham_mean =
      Sham_mean,
    percell_Tx_mean =
      Tx_mean,
    percell_delta =
      Tx_minus_Sham,
    percell_grade =
      evidence_grade
  ) %>%
  left_join(
    repair_weighted %>%
      filter(
        sender ==
          "Mphi_Repair-Resolution"
      ) %>%
      transmute(
        ligand,
        Repair_weighted_Sham_mean =
          Sham_mean,
        Repair_weighted_Tx_mean =
          Tx_mean,
        Repair_weighted_delta =
          Tx_minus_Sham,
        Repair_weighted_grade =
          evidence_grade
      ),
    by =
      "ligand"
  ) %>%
  left_join(
    total_weighted %>%
      transmute(
        ligand,
        total_Mphi5_weighted_Sham_mean =
          Sham_mean,
        total_Mphi5_weighted_Tx_mean =
          Tx_mean,
        total_Mphi5_weighted_delta =
          Tx_minus_Sham,
        total_Mphi5_weighted_grade =
          evidence_grade
      ),
    by =
      "ligand"
  ) %>%
  mutate(
    biological_level_interpretation =
      case_when(
        percell_delta <
          0 &
        Repair_weighted_delta >
          0 &
        total_Mphi5_weighted_delta <
          0 ~
          "Repair_percell_down_Repair_weighted_up_but_total_Mphi5_weighted_down",

        percell_delta <
          0 &
        Repair_weighted_delta <
          0 &
        total_Mphi5_weighted_delta <
          0 ~
          "Down_at_percell_Repair_weighted_and_total_Mphi5_levels",

        percell_delta >
          0 &
        Repair_weighted_delta >
          0 &
        total_Mphi5_weighted_delta >
          0 ~
          "Up_at_percell_Repair_weighted_and_total_Mphi5_levels",

        TRUE ~
          "Mixed_across_sender_levels"
      )
  )


write.csv(
  sender_level_audit,
  file.path(
    TAB_OUT,
    "10_RepairMphi_percell_vs_population_weighted_audit_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Figure: corrected receiver programs
# ==============================================================================

plot_program <- program_sample %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    hepatocyte_state =
      factor(
        hepatocyte_state,
        levels =
          PRIMARY_HEP_STATES
      )
  )


p1 <- ggplot(
  plot_program,
  aes(
    x =
      sample,
    y =
      module_mean_log1p_CP10k,
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
    hepatocyte_state ~ program,
    scales =
      "free_y"
  ) +
  labs(
    title =
      "Corrected independent Hepatocyte receiver programs",
    subtitle =
      "Sample-level pseudobulk mean log1p(CP10k)",
    x =
      NULL,
    y =
      "Program intensity"
  ) +
  theme_classic(
    base_size =
      7
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
  p1,
  file.path(
    FIG_OUT,
    "01_Hep_receiver_programs_CORRECTED_v6.6.4.4.pdf"
  ),
  16,
  7
)


# ==============================================================================
# 13. Figure: Pdgfrb+ vs Pdgfrb- mesenchymal burden
# ==============================================================================

p2_df <- pdgfrb_posneg_summary %>%
  mutate(
    Pdgfrb_status =
      ifelse(
        Pdgfrb_positive,
        "Pdgfrb+",
        "Pdgfrb-"
      )
  )


p2 <- ggplot(
  p2_df,
  aes(
    x =
      Pdgfrb_status,
    y =
      median_HSC_mesenchymal_n
  )
) +
  geom_col() +
  facet_wrap(
    ~ hepatocyte_state,
    nrow =
      1
  ) +
  labs(
    title =
      "Pdgfrb-positive Hepatocyte mesenchymal-marker burden",
    subtitle =
      "Median number of detected HSC/mesenchymal markers per cell",
    x =
      NULL,
    y =
      "Median HSC/mesenchymal genes detected"
  ) +
  theme_classic(
    base_size =
      8
  )


save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_Pdgfrb_positive_vs_negative_mesenchymal_burden_v6.6.4.4.pdf"
  ),
  8,
  4.5
)


# ==============================================================================
# 14. Figure: three sender levels
# ==============================================================================

sender_plot <- sender_level_audit %>%
  select(
    ligand,
    percell_delta,
    Repair_weighted_delta,
    total_Mphi5_weighted_delta
  ) %>%
  pivot_longer(
    cols = c(
      percell_delta,
      Repair_weighted_delta,
      total_Mphi5_weighted_delta
    ),
    names_to =
      "sender_level",
    values_to =
      "Tx_minus_Sham"
  ) %>%
  mutate(
    sender_level =
      recode(
        sender_level,
        percell_delta =
          "Repair-Mphi per-cell",
        Repair_weighted_delta =
          "Repair-Mphi abundance-weighted",
        total_Mphi5_weighted_delta =
          "Total Mphi5 abundance-weighted"
      )
  )


p3 <- ggplot(
  sender_plot,
  aes(
    x =
      sender_level,
    y =
      ligand,
    fill =
      sign(
        Tx_minus_Sham
      )
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
      0,
    limits = c(
      -1,
      1
    )
  ) +
  labs(
    title =
      "Sender-level direction audit",
    subtitle =
      "Blue = Tx down; red = Tx up",
    x =
      NULL,
    y =
      NULL,
    fill =
      "Direction"
  ) +
  theme_classic(
    base_size =
      8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle =
          35,
        hjust =
          1
      )
  )


save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_sender_level_direction_audit_v6.6.4.4.pdf"
  ),
  9,
  5
)


# ==============================================================================
# 15. Save correction object
# ==============================================================================

out_obj <- list(
  version =
    "v6.6.4.4",

  program_sample =
    program_sample,

  program_summary =
    program_summary,

  program_gene_audit =
    program_gene_audit,

  Pdgfrb_cell_audit =
    cell_audit,

  Pdgfrb_posneg_summary =
    pdgfrb_posneg_summary,

  Pdgfrb_marker_comparison =
    pdgfrb_marker_comparison,

  Pdgfrb_revised =
    pdgfrb_revised,

  closure_corrected =
    closure_corrected,

  key_axis_compact =
    key_axis_compact,

  sender_level_audit =
    sender_level_audit
)


saveRDS(
  out_obj,
  file.path(
    RDS_OUT,
    "Mouse_MASH_Mphi5_Hep5_mechanistic_closure_correction_v6.6.4.4.rds"
  ),
  compress =
    FALSE
)


# ==============================================================================
# 16. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "source_v6.6.0_RDS",
    "source_v6.6.4.2_closure",
    "CellChat_rerun",
    "NicheNet_rerun",
    "edgeR_rerun",
    "clustering_rerun",
    "receiver_program_fix",
    "Pdgfrb_audit_change",
    "CellChat_population_size_context",
    "biological_replicates"
  ),
  value = c(
    "v6.6.4.4",
    INPUT_RDS,
    OLD_CLOSURE,
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "Use named pseudobulk count vector instead of 2D one-dimensional character indexing",
    "Compare Pdgfrb-positive vs Pdgfrb-negative Hepatocytes; do not require loss of Hep identity to flag mesenchymal coexpression",
    "Source CellChat used population.size=FALSE; interpret primarily as per-cell communication state",
    "Sham n=2; Tx n=2"
  )
)


write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.4.4.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.4.4.txt"
    )
)


# ==============================================================================
# 17. Console summary
# ==============================================================================

msg(
  "Corrected receiver programs:"
)


print(
  program_summary %>%
    select(
      hepatocyte_state,
      program,
      Sham1,
      Sham20,
      Tx17,
      Tx5,
      Sham_mean,
      Tx_mean,
      Tx_minus_Sham,
      pairwise_direction_consistency,
      evidence_grade
    )
)


msg(
  "Revised Pdgfrb receiver audit:"
)


print(
  pdgfrb_revised %>%
    select(
      hepatocyte_state,
      receptor_detectability,
      Sham_mean_pct,
      Tx_mean_pct,
      Sham_mean_CP10k,
      Tx_mean_CP10k,
      median_HSC_mesenchymal_n_Pdgfrb_FALSE,
      median_HSC_mesenchymal_n_Pdgfrb_TRUE,
      pct_HSC_mesenchymal_ge3_Pdgfrb_FALSE,
      pct_HSC_mesenchymal_ge3_Pdgfrb_TRUE,
      pct_Hep_identity_ge5_Pdgfrb_TRUE,
      revised_Pdgfrb_interpretation
    )
)


msg(
  "Sender-level audit:"
)


print(
  sender_level_audit
)


msg(
  "Corrected final mechanistic closure:"
)


print(
  key_axis_compact %>%
    select(
      sender,
      receiver,
      ligand,
      receptor_audit_gene,
      sender_percell_delta,
      Repair_weighted_delta,
      total_Mphi5_weighted_delta,
      CellChat_best_supported_delta,
      nichenet_rank,
      receptor_detectability,
      receiver_program_proxy,
      receiver_program_delta_CORRECTED,
      receiver_program_grade_CORRECTED,
      sender_level_note,
      closure_status_v6643
    )
)


msg(
  "DONE."
)


msg(
  "No CellChat, NicheNet, edgeR, or clustering was rerun."
)


msg(
  "Output directory: ",
  OUT
)
