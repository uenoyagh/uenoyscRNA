#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6610)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Sample-wise CellChat:
#   5 macrophage subtypes <-> 5 Hepatocyte states
#
# Version: v6.6.1
#
# INPUT
# -----
#   Mouse_MASH_Mphi5_Hep5_interaction_ready_v6.6.0.rds
#
# DESIGN
# ------
# Biological samples are analyzed independently:
#   Sham1
#   Sham20
#   Tx17
#   Tx5
#
# Macrophage populations:
#   Mphi_Anti-inflammatory
#   Mphi_Inflammatory
#   Mphi_ECM-associated-inflammatory
#   Mphi_Repair-Resolution
#   Mphi_Lipid-associated-TREM2
#
# Hepatocyte populations:
#   Hep_Periportal
#   Hep_Pericentral
#   Hep_Injury-inflammatory
#   Hep_Intermediate
#   Hep_Cycling
#
# BOTH CROSS-LINEAGE DIRECTIONS ARE EXTRACTED:
#   A) Mphi -> Hepatocyte
#   B) Hepatocyte -> Mphi
#
# PRIMARY PRINCIPLE
# -----------------
# CellChat is run with:
#   population.size = FALSE
#   type = "triMean"
#
# Therefore population abundance is NOT multiplied directly into communication
# probability. Abundance changes remain a separate result.
#
# REPLICATE-AWARE COMPARISON
# --------------------------
# No pseudo-replication group-level test is performed.
# Each biological sample is retained independently.
#
# For each LR interaction and pathway:
#   Sham_mean = mean(Sham1, Sham20)
#   Tx_mean   = mean(Tx17, Tx5)
#
# Direction consistency is assessed using all four pairwise comparisons:
#   Tx17 - Sham1
#   Tx17 - Sham20
#   Tx5  - Sham1
#   Tx5  - Sham20
#
# IMPORTANT
# ---------
# This v6.6.1 is the primary sample-wise CellChat run.
# A later refinement script should explicitly incorporate CellChat p-values /
# supported probabilities before mechanistic ranking, as was done previously
# for the Mphi-HSC analysis.
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

  print(
    p
  )

  grDevices::dev.off()
}


safe_join_rna <- function(object) {

  if (
    "RNA" %in%
      Assays(
        object
      ) &&
    length(
      Layers(
        object[["RNA"]]
      )
    ) >
      1
  ) {

    object[["RNA"]] <- JoinLayers(
      object[["RNA"]]
    )
  }

  object
}


make_interaction_key <- function(df) {

  if (
    "interaction_name_2" %in%
      colnames(
        df
      )
  ) {

    return(
      as.character(
        df$interaction_name_2
      )
    )
  }

  if (
    "interaction_name" %in%
      colnames(
        df
      )
  ) {

    return(
      as.character(
        df$interaction_name
      )
    )
  }

  if (
    all(
      c(
        "ligand",
        "receptor"
      ) %in%
        colnames(
          df
        )
    )
  ) {

    return(
      paste(
        df$ligand,
        df$receptor,
        sep = " -> "
      )
    )
  }

  paste0(
    "interaction_",
    seq_len(
      nrow(
        df
      )
    )
  )
}


extract_lr_table <- function(
  cellchat,
  sample_name,
  condition_name,
  sources,
  targets,
  direction_name
) {

  df <- subsetCommunication(
    cellchat,
    sources.use = sources,
    targets.use = targets,
    thresh = 1
  )

  if (
    is.null(
      df
    ) ||
    !nrow(
      df
    )
  ) {

    return(
      tibble()
    )
  }

  df <- as_tibble(
    df
  )

  df$interaction_key <-
    make_interaction_key(
      df
    )

  df$sample <-
    sample_name

  df$condition <-
    condition_name

  df$direction <-
    direction_name

  df
}


pairwise_direction_metrics <- function(
  sham1,
  sham20,
  tx17,
  tx5
) {

  diffs <- c(
    tx17 -
      sham1,
    tx17 -
      sham20,
    tx5 -
      sham1,
    tx5 -
      sham20
  )

  up_fraction <- mean(
    diffs >
      0
  )

  down_fraction <- mean(
    diffs <
      0
  )

  direction <- dplyr::case_when(
    up_fraction >
      down_fraction ~
      "Tx_up",

    down_fraction >
      up_fraction ~
      "Tx_down",

    TRUE ~
      "Mixed"
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

  evidence_grade <- dplyr::case_when(
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
    pairwise_up_fraction =
      up_fraction,
    pairwise_down_fraction =
      down_fraction,
    pairwise_direction =
      direction,
    pairwise_direction_consistency =
      consistency,
    both_Tx_above_both_Sham =
      both_tx_above,
    both_Tx_below_both_Sham =
      both_tx_below,
    evidence_grade =
      evidence_grade
  )
}


summarize_four_samples <- function(
  df,
  key_cols,
  value_col = "prob"
) {

  sample_levels <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  all_keys <- df %>%
    distinct(
      across(
        all_of(
          key_cols
        )
      )
    )

  grid <- tidyr::crossing(
    all_keys,
    sample =
      sample_levels
  )

  value_df <- df %>%
    group_by(
      across(
        all_of(
          key_cols
        )
      ),
      sample
    ) %>%
    summarise(
      value =
        sum(
          .data[[
            value_col
          ]],
          na.rm = TRUE
        ),
      .groups =
        "drop"
    )

  complete_df <- grid %>%
    left_join(
      value_df,
      by = c(
        key_cols,
        "sample"
      )
    ) %>%
    mutate(
      value =
        replace_na(
          value,
          0
        )
    )

  wide <- complete_df %>%
    pivot_wider(
      names_from =
        sample,
      values_from =
        value,
      values_fill =
        0
    )

  required_samples <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  missing_samples <- setdiff(
    required_samples,
    colnames(
      wide
    )
  )

  if (
    length(
      missing_samples
    )
  ) {

    for (
      sample_name in missing_samples
    ) {

      wide[[
        sample_name
      ]] <- 0
    }
  }

  metrics <- lapply(
    seq_len(
      nrow(
        wide
      )
    ),
    function(i) {

      pairwise_direction_metrics(
        sham1 =
          wide$Sham1[[
            i
          ]],
        sham20 =
          wide$Sham20[[
            i
          ]],
        tx17 =
          wide$Tx17[[
            i
          ]],
        tx5 =
          wide$Tx5[[
            i
          ]]
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
          Sham_mean,

      log2_Tx_vs_Sham =
        log2(
          (
            Tx_mean +
              1e-12
          ) /
            (
              Sham_mean +
                1e-12
            )
        ),

      detected_samples =
        (
          Sham1 >
            0
        ) +
          (
            Sham20 >
              0
          ) +
          (
            Tx17 >
              0
          ) +
          (
            Tx5 >
              0
          ),

      detected_Sham_replicates =
        (
          Sham1 >
            0
        ) +
          (
            Sham20 >
              0
          ),

      detected_Tx_replicates =
        (
          Tx17 >
            0
        ) +
          (
            Tx5 >
              0
          )
    )
}


priority_pathway_flag <- function(pathway_name) {

  grepl(
    paste(
      c(
        "^TNF$",
        "^IL1$",
        "^IL6$",
        "^OSM$",
        "^TGF",
        "TGFB",
        "SPP1",
        "FN1",
        "GAS",
        "AXL",
        "MERTK",
        "EGF",
        "IGF",
        "FGF",
        "CCL",
        "CXCL",
        "CSF",
        "MIF",
        "ANNEXIN",
        "HEDGEHOG",
        "SHH",
        "NOTCH",
        "VEGF",
        "COLLAGEN",
        "LAMININ",
        "THBS"
      ),
      collapse = "|"
    ),
    pathway_name,
    ignore.case = TRUE
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


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_v6.6.1"
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

GROUP_COL <-
  "interaction_celltype_v660"


SAMPLE_COL <-
  "sample_interaction_v660"


CONDITION_COL <-
  "condition_interaction_v660"


SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)


MPHI_STATES <- c(
  "Mphi_Anti-inflammatory",
  "Mphi_Inflammatory",
  "Mphi_ECM-associated-inflammatory",
  "Mphi_Repair-Resolution",
  "Mphi_Lipid-associated-TREM2"
)


HEP_STATES <- c(
  "Hep_Periportal",
  "Hep_Pericentral",
  "Hep_Injury-inflammatory",
  "Hep_Intermediate",
  "Hep_Cycling"
)


INTERACTION_LEVELS <- c(
  MPHI_STATES,
  HEP_STATES
)


MIN_CELLS <-
  10


# ==============================================================================
# 4. Preflight
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {

  stop(
    "Input interaction-ready RDS not found: ",
    INPUT_RDS
  )
}


cellchat_version <- as.character(
  utils::packageVersion(
    "CellChat"
  )
)


seurat_version <- as.character(
  utils::packageVersion(
    "Seurat"
  )
)


msg(
  "CellChat version: ",
  cellchat_version
)


msg(
  "Seurat version: ",
  seurat_version
)


msg(
  "Loading interaction-ready object..."
)


obj <- readRDS(
  INPUT_RDS
)


if (
  !"RNA" %in%
    Assays(
      obj
    )
) {

  stop(
    "RNA assay missing."
  )
}


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
    "Missing metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


available_groups <- unique(
  as.character(
    obj@meta.data[[
      GROUP_COL
    ]]
  )
)


missing_groups <- setdiff(
  INTERACTION_LEVELS,
  available_groups
)


if (
  length(
    missing_groups
  )
) {

  stop(
    "Missing required interaction populations: ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 5. Cell-count audit
# ==============================================================================

count_audit <- obj@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    interaction_celltype =
      as.character(
        .data[[
          GROUP_COL
        ]]
      )
  ) %>%
  count(
    sample,
    interaction_celltype,
    name =
      "n_cells"
  ) %>%
  complete(
    sample =
      SAMPLES,
    interaction_celltype =
      INTERACTION_LEVELS,
    fill =
      list(
        n_cells =
          0
      )
  ) %>%
  arrange(
    factor(
      sample,
      levels =
        SAMPLES
    ),
    factor(
      interaction_celltype,
      levels =
        INTERACTION_LEVELS
    )
  )


write.csv(
  count_audit,
  file.path(
    TAB_OUT,
    "01_sample_by_celltype_cellcount_audit_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


if (
  any(
    count_audit$n_cells <
      MIN_CELLS
  )
) {

  bad <- count_audit %>%
    filter(
      n_cells <
        MIN_CELLS
    )

  print(
    bad
  )

  stop(
    "At least one sample x celltype group has fewer than ",
    MIN_CELLS,
    " cells. CellChat aborted."
  )
}


# ==============================================================================
# 6. CellChat database audit
# ==============================================================================

DB_USE <-
  CellChatDB.mouse


db_categories <- DB_USE$interaction %>%
  count(
    annotation,
    name =
      "n_interactions"
  )


write.csv(
  db_categories,
  file.path(
    TAB_OUT,
    "02_CellChatDB_mouse_category_audit_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 7. Sample-wise CellChat
# ==============================================================================

cellchat_objects <- list()

mphi_to_hep_lr_tables <- list()

hep_to_mphi_lr_tables <- list()


for (
  sample_name in SAMPLES
) {

  msg(
    "============================================================"
  )

  msg(
    "CellChat sample: ",
    sample_name
  )


  sample_cells <- colnames(
    obj
  )[
    as.character(
      obj@meta.data[[
        SAMPLE_COL
      ]]
    ) ==
      sample_name
  ]


  sample_obj <- subset(
    obj,
    cells =
      sample_cells
  )


  sample_obj <- safe_join_rna(
    sample_obj
  )


  DefaultAssay(
    sample_obj
  ) <- "RNA"


  sample_meta <- sample_obj@meta.data %>%
    as.data.frame()


  sample_meta[[
    GROUP_COL
  ]] <- factor(
    as.character(
      sample_meta[[
        GROUP_COL
      ]]
    ),
    levels =
      INTERACTION_LEVELS
  )


  group_count <- table(
    sample_meta[[
      GROUP_COL
    ]]
  )


  msg(
    "Group counts: ",
    paste(
      names(
        group_count
      ),
      as.integer(
        group_count
      ),
      sep =
        "=",
      collapse =
        "; "
    )
  )


  data_input <- GetAssayData(
    sample_obj,
    assay =
      "RNA",
    layer =
      "data"
  )


  meta_input <- sample_meta[
    ,
    c(
      GROUP_COL
    ),
    drop =
      FALSE
  ]


  cellchat <- createCellChat(
    object =
      data_input,
    meta =
      meta_input,
    group.by =
      GROUP_COL
  )


  cellchat@DB <-
    DB_USE


  msg(
    sample_name,
    ": subsetData"
  )


  cellchat <- subsetData(
    cellchat
  )


  msg(
    sample_name,
    ": identifyOverExpressedGenes"
  )


  cellchat <- identifyOverExpressedGenes(
    cellchat
  )


  msg(
    sample_name,
    ": identifyOverExpressedInteractions"
  )


  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )


  msg(
    sample_name,
    ": computeCommunProb population.size=FALSE"
  )


  cellchat <- computeCommunProb(
    cellchat,
    type =
      "triMean",
    population.size =
      FALSE
  )


  cellchat <- filterCommunication(
    cellchat,
    min.cells =
      MIN_CELLS
  )


  msg(
    sample_name,
    ": computeCommunProbPathway"
  )


  cellchat <- computeCommunProbPathway(
    cellchat
  )


  msg(
    sample_name,
    ": aggregateNet"
  )


  cellchat <- aggregateNet(
    cellchat
  )


  condition_name <- ifelse(
    grepl(
      "^Sham",
      sample_name
    ),
    "Sham",
    "Tx"
  )


  # --------------------------------------------------------------------------
  # Mphi -> Hepatocyte
  # --------------------------------------------------------------------------

  lr_mphi_to_hep <- extract_lr_table(
    cellchat =
      cellchat,
    sample_name =
      sample_name,
    condition_name =
      condition_name,
    sources =
      MPHI_STATES,
    targets =
      HEP_STATES,
    direction_name =
      "Mphi_to_Hep"
  )


  if (
    nrow(
      lr_mphi_to_hep
    )
  ) {

    mphi_to_hep_lr_tables[[
      sample_name
    ]] <- lr_mphi_to_hep

  } else {

    warning(
      "No Mphi -> Hepatocyte LR interactions extracted for ",
      sample_name
    )
  }


  # --------------------------------------------------------------------------
  # Hepatocyte -> Mphi
  # --------------------------------------------------------------------------

  lr_hep_to_mphi <- extract_lr_table(
    cellchat =
      cellchat,
    sample_name =
      sample_name,
    condition_name =
      condition_name,
    sources =
      HEP_STATES,
    targets =
      MPHI_STATES,
    direction_name =
      "Hep_to_Mphi"
  )


  if (
    nrow(
      lr_hep_to_mphi
    )
  ) {

    hep_to_mphi_lr_tables[[
      sample_name
    ]] <- lr_hep_to_mphi

  } else {

    warning(
      "No Hepatocyte -> Mphi LR interactions extracted for ",
      sample_name
    )
  }


  cellchat_objects[[
    sample_name
  ]] <- cellchat


  saveRDS(
    cellchat,
    file.path(
      RDS_OUT,
      paste0(
        "CellChat_",
        sample_name,
        "_v6.6.1.rds"
      )
    ),
    compress =
      FALSE
  )


  rm(
    sample_obj,
    data_input,
    meta_input,
    cellchat
  )


  gc()
}


# ==============================================================================
# 8. Combine LR tables
# ==============================================================================

if (
  !length(
    mphi_to_hep_lr_tables
  )
) {

  stop(
    "No Mphi -> Hepatocyte LR tables were generated."
  )
}


if (
  !length(
    hep_to_mphi_lr_tables
  )
) {

  stop(
    "No Hepatocyte -> Mphi LR tables were generated."
  )
}


lr_mphi_to_hep_all <- bind_rows(
  mphi_to_hep_lr_tables
)


lr_hep_to_mphi_all <- bind_rows(
  hep_to_mphi_lr_tables
)


lr_cross_all <- bind_rows(
  lr_mphi_to_hep_all,
  lr_hep_to_mphi_all
)


write.csv(
  lr_mphi_to_hep_all,
  file.path(
    TAB_OUT,
    "03_Mphi5_to_Hep5_LR_all_samples_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  lr_hep_to_mphi_all,
  file.path(
    TAB_OUT,
    "04_Hep5_to_Mphi5_LR_all_samples_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  lr_cross_all,
  file.path(
    TAB_OUT,
    "05_bidirectional_crosslineage_LR_all_samples_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


msg(
  "Mphi -> Hep sample-level LR rows: ",
  nrow(
    lr_mphi_to_hep_all
  )
)


msg(
  "Hep -> Mphi sample-level LR rows: ",
  nrow(
    lr_hep_to_mphi_all
  )
)


# ==============================================================================
# 9. Standardize LR columns
# ==============================================================================

required_lr_cols <- c(
  "source",
  "target",
  "prob",
  "pathway_name",
  "interaction_key",
  "sample",
  "condition",
  "direction"
)


missing_lr_cols <- setdiff(
  required_lr_cols,
  colnames(
    lr_cross_all
  )
)


if (
  length(
    missing_lr_cols
  )
) {

  stop(
    "Required LR columns missing from CellChat subsetCommunication output: ",
    paste(
      missing_lr_cols,
      collapse = ", "
    )
  )
}


if (
  !"ligand" %in%
    colnames(
      lr_cross_all
    )
) {

  lr_cross_all$ligand <-
    NA_character_
}


if (
  !"receptor" %in%
    colnames(
      lr_cross_all
    )
) {

  lr_cross_all$receptor <-
    NA_character_
}


# ==============================================================================
# 10. LR replicate-aware Sham vs Tx summary
# ==============================================================================

lr_key_cols <- c(
  "direction",
  "source",
  "target",
  "pathway_name",
  "interaction_key",
  "ligand",
  "receptor"
)


lr_compare <- summarize_four_samples(
  df =
    lr_cross_all,
  key_cols =
    lr_key_cols,
  value_col =
    "prob"
) %>%
  mutate(
    priority_pathway =
      priority_pathway_flag(
        pathway_name
      )
  ) %>%
  arrange(
    direction,
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  lr_compare,
  file.path(
    TAB_OUT,
    "06_bidirectional_LR_Sham_vs_Tx_replicate_summary_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  lr_compare %>%
    filter(
      direction ==
        "Mphi_to_Hep"
    ),
  file.path(
    TAB_OUT,
    "07_Mphi_to_Hep_LR_Sham_vs_Tx_summary_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  lr_compare %>%
    filter(
      direction ==
        "Hep_to_Mphi"
    ),
  file.path(
    TAB_OUT,
    "08_Hep_to_Mphi_LR_Sham_vs_Tx_summary_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 11. Pathway sample-level summaries
# ==============================================================================

pathway_sample <- lr_cross_all %>%
  group_by(
    sample,
    condition,
    direction,
    source,
    target,
    pathway_name
  ) %>%
  summarise(
    pathway_strength =
      sum(
        prob,
        na.rm = TRUE
      ),

    n_LR =
      n_distinct(
        interaction_key
      ),

    .groups =
      "drop"
  )


write.csv(
  pathway_sample,
  file.path(
    TAB_OUT,
    "09_bidirectional_pathway_strength_by_sample_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 12. Pathway replicate-aware Sham vs Tx summary
# ==============================================================================

pathway_compare <- summarize_four_samples(
  df =
    pathway_sample %>%
    rename(
      prob =
        pathway_strength
    ),
  key_cols = c(
    "direction",
    "source",
    "target",
    "pathway_name"
  ),
  value_col =
    "prob"
) %>%
  mutate(
    priority_pathway =
      priority_pathway_flag(
        pathway_name
      )
  ) %>%
  arrange(
    direction,
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    ),
    desc(
      Tx_mean +
        Sham_mean
    )
  )


write.csv(
  pathway_compare,
  file.path(
    TAB_OUT,
    "10_bidirectional_pathway_Sham_vs_Tx_summary_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  pathway_compare %>%
    filter(
      priority_pathway
    ),
  file.path(
    TAB_OUT,
    "11_priority_crosslineage_pathways_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 13. Biological focus A:
#     Mphi -> Injury/inflammatory Hepatocyte
# ==============================================================================

mphi_to_injury_hep_lr <- lr_compare %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  mphi_to_injury_hep_lr,
  file.path(
    TAB_OUT,
    "12_Mphi5_to_InjuryHep_LR_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


mphi_to_injury_hep_pathway <- pathway_compare %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    target ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  mphi_to_injury_hep_pathway,
  file.path(
    TAB_OUT,
    "13_Mphi5_to_InjuryHep_pathway_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 14. Biological focus B:
#     Injury/inflammatory Hepatocyte -> Mphi5
# ==============================================================================

injury_hep_to_mphi_lr <- lr_compare %>%
  filter(
    direction ==
      "Hep_to_Mphi",
    source ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  injury_hep_to_mphi_lr,
  file.path(
    TAB_OUT,
    "14_InjuryHep_to_Mphi5_LR_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


injury_hep_to_mphi_pathway <- pathway_compare %>%
  filter(
    direction ==
      "Hep_to_Mphi",
    source ==
      "Hep_Injury-inflammatory"
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  injury_hep_to_mphi_pathway,
  file.path(
    TAB_OUT,
    "15_InjuryHep_to_Mphi5_pathway_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 15. Biological focus C:
#     Anti-inflammatory / Repair-Resolution Mphi -> Hepatocyte
# ==============================================================================

protective_mphi_to_hep <- lr_compare %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    source %in%
      c(
        "Mphi_Anti-inflammatory",
        "Mphi_Repair-Resolution"
      )
  ) %>%
  arrange(
    source,
    target,
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      Tx_minus_Sham
    )
  )


write.csv(
  protective_mphi_to_hep,
  file.path(
    TAB_OUT,
    "16_AntiInflammatory_RepairMphi_to_Hep_LR_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 16. Biological focus D:
#     Inflammatory / ECM-inflammatory Mphi -> Hepatocyte
# ==============================================================================

inflammatory_mphi_to_hep <- lr_compare %>%
  filter(
    direction ==
      "Mphi_to_Hep",
    source %in%
      c(
        "Mphi_Inflammatory",
        "Mphi_ECM-associated-inflammatory"
      )
  ) %>%
  arrange(
    source,
    target,
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  inflammatory_mphi_to_hep,
  file.path(
    TAB_OUT,
    "17_Inflammatory_ECMMphi_to_Hep_LR_focus_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 17. Candidate ligand/receptor focused table
# ==============================================================================

CANDIDATE_REGEX <- paste(
  c(
    "^Tnf$",
    "^Il1",
    "^Il6$",
    "^Osm$",
    "^Tgfb1$",
    "^Spp1$",
    "^Fn1$",
    "^Gas6$",
    "^Areg$",
    "^Hbegf$",
    "^Igf1$",
    "^Ccl2$",
    "^Cxcl1$",
    "^Cxcl10$",
    "^Csf1$",
    "^Saa1$",
    "^Saa2$",
    "^Shh$",
    "^Ihh$"
  ),
  collapse =
    "|"
)


candidate_lr <- lr_compare %>%
  filter(
    grepl(
      CANDIDATE_REGEX,
      ligand,
      ignore.case =
        TRUE
    )
  ) %>%
  arrange(
    direction,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        Tx_minus_Sham
      )
    )
  )


write.csv(
  candidate_lr,
  file.path(
    TAB_OUT,
    "18_candidate_ligand_LR_summary_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


# ==============================================================================
# 18. Figure: bidirectional pathway Tx-vs-Sham heatmaps
# ==============================================================================

for (
  direction_name in c(
    "Mphi_to_Hep",
    "Hep_to_Mphi"
  )
) {

  plot_pathways <- pathway_compare %>%
    filter(
      direction ==
        direction_name,
      detected_samples >=
        2
    ) %>%
    group_by(
      pathway_name
    ) %>%
    mutate(
      pathway_max_strength =
        max(
          Tx_mean +
            Sham_mean,
          na.rm =
            TRUE
        )
    ) %>%
    ungroup() %>%
    arrange(
      desc(
        priority_pathway
      ),
      desc(
        pathway_max_strength
      )
    )


  pathway_keep <- plot_pathways %>%
    distinct(
      pathway_name,
      priority_pathway,
      pathway_max_strength
    ) %>%
    arrange(
      desc(
        priority_pathway
      ),
      desc(
        pathway_max_strength
      )
    ) %>%
    slice_head(
      n =
        30
    ) %>%
    pull(
      pathway_name
    )


  heat_df <- plot_pathways %>%
    filter(
      pathway_name %in%
        pathway_keep
    ) %>%
    mutate(
      sender_receiver =
        paste0(
          source,
          " -> ",
          target
        ),

      pathway_name =
        factor(
          pathway_name,
          levels =
            rev(
              pathway_keep
            )
        )
    )


  if (
    nrow(
      heat_df
    )
  ) {

    p_heat <- ggplot(
      heat_df,
      aes(
        x =
          sender_receiver,
        y =
          pathway_name,
        fill =
          log2_Tx_vs_Sham
      )
    ) +
      geom_tile(
        linewidth =
          0.20
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
          -5,
          5
        ),
        oob =
          scales::squish
      ) +
      labs(
        title =
          paste0(
            direction_name,
            " pathway change | Tx vs Sham"
          ),
        subtitle =
          "Red = higher in Tx; blue = lower in Tx; sample-wise CellChat",
        x =
          NULL,
        y =
          NULL,
        fill =
          "log2 Tx/Sham"
      ) +
      theme_classic(
        base_size =
          6.5
      ) +
      theme(
        axis.text.x =
          element_text(
            angle =
              60,
            hjust =
              1
          ),
        plot.title =
          element_text(
            face =
              "bold",
            hjust =
              0.5
          )
      )


    file_name <- if (
      direction_name ==
        "Mphi_to_Hep"
    ) {
      "01_Mphi5_to_Hep5_pathway_Tx_vs_Sham_heatmap_v6.6.1.pdf"
    } else {
      "02_Hep5_to_Mphi5_pathway_Tx_vs_Sham_heatmap_v6.6.1.pdf"
    }


    save_pdf(
      p_heat,
      file.path(
        FIG_OUT,
        file_name
      ),
      18,
      10
    )
  }
}


# ==============================================================================
# 19. Figure: Mphi -> Injury Hep top LR
# ==============================================================================

injury_plot_df <- mphi_to_injury_hep_lr %>%
  filter(
    detected_samples >=
      2
  ) %>%
  mutate(
    max_strength =
      pmax(
        Sham_mean,
        Tx_mean
      )
  ) %>%
  group_by(
    source
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      max_strength
    ),
    .by_group =
      TRUE
  ) %>%
  slice_head(
    n =
      15
  ) %>%
  ungroup() %>%
  mutate(
    interaction_label =
      ifelse(
        !is.na(
          interaction_key
        ),
        interaction_key,
        paste0(
          ligand,
          " -> ",
          receptor
        )
      )
  )


if (
  nrow(
    injury_plot_df
  )
) {

  p_injury <- ggplot(
    injury_plot_df,
    aes(
      x =
        Tx_minus_Sham,
      y =
        reorder(
          interaction_label,
          Tx_minus_Sham
        ),
      size =
        max_strength
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.4
    ) +
    facet_wrap(
      ~ source,
      scales =
        "free_y",
      ncol =
        1
    ) +
    labs(
      title =
        "Mphi -> Injury/inflammatory Hepatocyte",
      subtitle =
        "Top CellChat LR interactions; absolute Tx-Sham probability difference",
      x =
        "Tx mean - Sham mean",
      y =
        NULL,
      size =
        "Max mean strength"
    ) +
    theme_classic(
      base_size =
        7.5
    ) +
    theme(
      plot.title =
        element_text(
          face =
            "bold",
          hjust =
            0.5
        )
    )


  save_pdf(
    p_injury,
    file.path(
      FIG_OUT,
      "03_Mphi5_to_InjuryHep_top_LR_v6.6.1.pdf"
    ),
    10,
    16
  )
}


# ==============================================================================
# 20. Figure: Injury Hep -> Mphi top LR
# ==============================================================================

reverse_plot_df <- injury_hep_to_mphi_lr %>%
  filter(
    detected_samples >=
      2
  ) %>%
  mutate(
    max_strength =
      pmax(
        Sham_mean,
        Tx_mean
      )
  ) %>%
  group_by(
    target
  ) %>%
  arrange(
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    desc(
      max_strength
    ),
    .by_group =
      TRUE
  ) %>%
  slice_head(
    n =
      15
  ) %>%
  ungroup() %>%
  mutate(
    interaction_label =
      ifelse(
        !is.na(
          interaction_key
        ),
        interaction_key,
        paste0(
          ligand,
          " -> ",
          receptor
        )
      )
  )


if (
  nrow(
    reverse_plot_df
  )
) {

  p_reverse <- ggplot(
    reverse_plot_df,
    aes(
      x =
        Tx_minus_Sham,
      y =
        reorder(
          interaction_label,
          Tx_minus_Sham
        ),
      size =
        max_strength
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept =
        0,
      linewidth =
        0.4
    ) +
    facet_wrap(
      ~ target,
      scales =
        "free_y",
      ncol =
        1
    ) +
    labs(
      title =
        "Injury/inflammatory Hepatocyte -> Mphi",
      subtitle =
        "Top CellChat LR interactions; absolute Tx-Sham probability difference",
      x =
        "Tx mean - Sham mean",
      y =
        NULL,
      size =
        "Max mean strength"
    ) +
    theme_classic(
      base_size =
        7.5
    ) +
    theme(
      plot.title =
        element_text(
          face =
            "bold",
          hjust =
            0.5
        )
    )


  save_pdf(
    p_reverse,
    file.path(
      FIG_OUT,
      "04_InjuryHep_to_Mphi5_top_LR_v6.6.1.pdf"
    ),
    10,
    16
  )
}


# ==============================================================================
# 21. Figure: sample-level priority pathway heatmap
# ==============================================================================

priority_sample <- pathway_sample %>%
  filter(
    priority_pathway_flag(
      pathway_name
    )
  ) %>%
  mutate(
    sender_receiver_pathway =
      paste0(
        direction,
        " | ",
        source,
        " -> ",
        target,
        " | ",
        pathway_name
      )
  )


top_priority_keys <- priority_sample %>%
  group_by(
    sender_receiver_pathway
  ) %>%
  summarise(
    max_strength =
      max(
        pathway_strength,
        na.rm =
          TRUE
      ),
    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      max_strength
    )
  ) %>%
  slice_head(
    n =
      50
  ) %>%
  pull(
    sender_receiver_pathway
  )


priority_sample_plot <- priority_sample %>%
  filter(
    sender_receiver_pathway %in%
      top_priority_keys
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),

    sender_receiver_pathway =
      factor(
        sender_receiver_pathway,
        levels =
          rev(
            top_priority_keys
          )
      )
  )


if (
  nrow(
    priority_sample_plot
  )
) {

  p_sample_heat <- ggplot(
    priority_sample_plot,
    aes(
      x =
        sample,
      y =
        sender_receiver_pathway,
      fill =
        pathway_strength
    )
  ) +
    geom_tile(
      linewidth =
        0.20
    ) +
    scale_fill_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    labs(
      title =
        "Priority Mphi <-> Hepatocyte pathways by biological sample",
      x =
        NULL,
      y =
        NULL,
      fill =
        "CellChat strength"
    ) +
    theme_classic(
      base_size =
        6.5
    ) +
    theme(
      plot.title =
        element_text(
          face =
            "bold",
          hjust =
            0.5
        )
    )


  save_pdf(
    p_sample_heat,
    file.path(
      FIG_OUT,
      "05_priority_bidirectional_pathways_by_sample_v6.6.1.pdf"
    ),
    10,
    15
  )
}


# ==============================================================================
# 22. Save CellChat collection
# ==============================================================================

saveRDS(
  cellchat_objects,
  file.path(
    RDS_OUT,
    "CellChat_samplewise_list_v6.6.1.rds"
  ),
  compress =
    FALSE
)


merged_cellchat <- tryCatch(
  mergeCellChat(
    cellchat_objects,
    add.names =
      names(
        cellchat_objects
      ),
    cell.prefix =
      TRUE
  ),
  error =
    function(e) {

      warning(
        "mergeCellChat failed: ",
        conditionMessage(
          e
        ),
        ". Sample-wise objects and summary tables remain valid."
      )

      NULL
    }
)


if (
  !is.null(
    merged_cellchat
  )
) {

  saveRDS(
    merged_cellchat,
    file.path(
      RDS_OUT,
      "CellChat_merged_4samples_v6.6.1.rds"
    ),
    compress =
      FALSE
  )
}


# ==============================================================================
# 23. Analysis metadata
# ==============================================================================

analysis_metadata <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "CellChat_version",
    "Seurat_version",
    "CellChat_database",
    "population_size",
    "computeCommunProb_type",
    "min_cells",
    "samples",
    "macrophage_states",
    "hepatocyte_states",
    "directions_analyzed",
    "group_level_statistics"
  ),
  value = c(
    "v6.6.1",
    INPUT_RDS,
    cellchat_version,
    seurat_version,
    "CellChatDB.mouse full database",
    "FALSE",
    "triMean",
    as.character(
      MIN_CELLS
    ),
    paste(
      SAMPLES,
      collapse =
        ","
    ),
    paste(
      MPHI_STATES,
      collapse =
        " | "
    ),
    paste(
      HEP_STATES,
      collapse =
        " | "
    ),
    "Mphi_to_Hep | Hep_to_Mphi",
    "None; exploratory sample-level n=2 Sham vs n=2 Tx"
  )
)


write.csv(
  analysis_metadata,
  file.path(
    LOG_OUT,
    "analysis_metadata_v6.6.1.csv"
  ),
  row.names =
    FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.1.txt"
    )
)


# ==============================================================================
# 24. Console summary
# ==============================================================================

msg(
  "Top Mphi -> Hep pathway changes:"
)


print(
  pathway_compare %>%
    filter(
      direction ==
        "Mphi_to_Hep",
      priority_pathway
    ) %>%
    select(
      source,
      target,
      pathway_name,
      Sham_mean,
      Tx_mean,
      Tx_minus_Sham,
      pairwise_direction_consistency,
      evidence_grade
    ) %>%
    arrange(
      desc(
        pairwise_direction_consistency
      ),
      desc(
        abs(
          Tx_minus_Sham
        )
      )
    ) %>%
    head(
      20
    )
)


msg(
  "Top Hep -> Mphi pathway changes:"
)


print(
  pathway_compare %>%
    filter(
      direction ==
        "Hep_to_Mphi",
      priority_pathway
    ) %>%
    select(
      source,
      target,
      pathway_name,
      Sham_mean,
      Tx_mean,
      Tx_minus_Sham,
      pairwise_direction_consistency,
      evidence_grade
    ) %>%
    arrange(
      desc(
        pairwise_direction_consistency
      ),
      desc(
        abs(
          Tx_minus_Sham
        )
      )
    ) %>%
    head(
      20
    )
)


msg(
  "DONE."
)


msg(
  "Output directory: ",
  OUT
)
