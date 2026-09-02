#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6210)

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
# Sample-wise CellChat: 5 Mphi subtypes -> 3 HSC states
#
# Version: v6.2.1
#
# INPUT
#   Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds
#
# DESIGN
#   Biological samples are analyzed independently:
#     Sham1
#     Sham20
#     Tx17
#     Tx5
#
#   Sender populations:
#     Anti-inflammatory-Mphi
#     Inflammatory-Mphi
#     ECM-associated inflammatory-Mphi
#     Repair/Resolution-Mphi
#     Lipid-associated/TREM2-Mphi
#
#   Receiver populations:
#     qHSC
#     ECM-activated HSC
#     Contractile HSC
#
# PRIMARY PRINCIPLE
#   CellChat is run with population.size = FALSE.
#   Therefore changes in population abundance are NOT directly multiplied into
#   communication probability. Population-size changes remain a separate result.
#
# REPLICATE-AWARE COMPARISON
#   No pseudo-replication statistical test is performed.
#   For each LR interaction and pathway, sample-level probabilities are retained.
#   Sham mean and Tx mean are calculated from the two biological replicates.
#   Direction consistency is assessed from all 2 x 2 Tx-vs-Sham sample pairs:
#
#     Tx17 - Sham1
#     Tx17 - Sham20
#     Tx5  - Sham1
#     Tx5  - Sham20
#
#   consistency = 1.0 means all four pairwise comparisons have the same direction.
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

make_interaction_key <- function(df) {

  if (
    "interaction_name_2" %in%
      colnames(df)
  ) {
    return(
      as.character(
        df$interaction_name_2
      )
    )
  }

  if (
    "interaction_name" %in%
      colnames(df)
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
        colnames(df)
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
      nrow(df)
    )
  )
}

extract_lr_table <- function(
  cellchat,
  sample_name,
  condition_name,
  sources,
  targets
) {

  df <- subsetCommunication(
    cellchat,
    sources.use = sources,
    targets.use = targets,
    thresh = 1
  )

  if (
    is.null(df) ||
    !nrow(df)
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

  df
}

pairwise_direction_metrics <- function(
  sham1,
  sham20,
  tx17,
  tx5
) {

  diffs <- c(
    tx17 - sham1,
    tx17 - sham20,
    tx5 - sham1,
    tx5 - sham20
  )

  up_fraction <-
    mean(
      diffs > 0
    )

  down_fraction <-
    mean(
      diffs < 0
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

  consistency <-
    max(
      up_fraction,
      down_fraction
    )

  tibble(
    pairwise_up_fraction =
      up_fraction,
    pairwise_down_fraction =
      down_fraction,
    pairwise_direction =
      direction,
    pairwise_direction_consistency =
      consistency
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
    sample = sample_levels
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
      .groups = "drop"
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
      names_from = sample,
      values_from = value,
      values_fill = 0
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
      smp in missing_samples
    ) {
      wide[[smp]] <- 0
    }
  }

  metrics <- lapply(
    seq_len(
      nrow(wide)
    ),
    function(i) {

      pairwise_direction_metrics(
        sham1 =
          wide$Sham1[[i]],
        sham20 =
          wide$Sham20[[i]],
        tx17 =
          wide$Tx17[[i]],
        tx5 =
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
        ) / 2,

      Tx_mean =
        (
          Tx17 +
            Tx5
        ) / 2,

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
          Sham1 > 0
        ) +
        (
          Sham20 > 0
        ) +
        (
          Tx17 > 0
        ) +
        (
          Tx5 > 0
        ),

      detected_Sham_replicates =
        (
          Sham1 > 0
        ) +
        (
          Sham20 > 0
        ),

      detected_Tx_replicates =
        (
          Tx17 > 0
        ) +
        (
          Tx5 > 0
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

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_CellChat_v6.2.1"
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

SENDERS <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

RECEIVERS <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

INTERACTION_LEVELS <- c(
  SENDERS,
  RECEIVERS
)

MIN_CELLS <- 10

PRIORITY_PATHWAY_REGEX <- paste(
  c(
    "^TGF",
    "THBS",
    "SPP1",
    "PDGF",
    "GAS",
    "AXL",
    "MERTK",
    "CCL",
    "CXCL",
    "IL10",
    "TNF",
    "VEGF",
    "FGF",
    "IGF",
    "COLLAGEN",
    "LAMININ",
    "FN1"
  ),
  collapse = "|"
)


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

cellchat_version <-
  as.character(
    utils::packageVersion(
      "CellChat"
    )
  )

seurat_version <-
  as.character(
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
    Assays(obj)
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

available_groups <-
  unique(
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
    rownames = "cell"
  ) %>%
  count(
    .data[[
      SAMPLE_COL
    ]],
    .data[[
      GROUP_COL
    ]],
    name = "n_cells"
  ) %>%
  complete(
    !!sym(
      SAMPLE_COL
    ) :=
      SAMPLES,

    !!sym(
      GROUP_COL
    ) :=
      INTERACTION_LEVELS,

    fill = list(
      n_cells = 0
    )
  )

write.csv(
  count_audit,
  file.path(
    TAB_OUT,
    "01_sample_by_celltype_cellcount_audit_v6.2.1.csv"
  ),
  row.names = FALSE
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
# 6. CellChat database
# ==============================================================================

DB_USE <- CellChatDB.mouse

db_categories <- DB_USE$interaction %>%
  count(
    annotation,
    name = "n_interactions"
  )

write.csv(
  db_categories,
  file.path(
    TAB_OUT,
    "02_CellChatDB_mouse_category_audit_v6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Sample-wise CellChat
# ==============================================================================

cellchat_objects <- list()
lr_tables <- list()

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
    cells = sample_cells
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
  ]] <-
    factor(
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
      sep = "=",
      collapse = "; "
    )
  )

  data_input <- GetAssayData(
    sample_obj,
    assay = "RNA",
    layer = "data"
  )

  meta_input <- sample_meta[
    ,
    c(
      GROUP_COL
    ),
    drop = FALSE
  ]

  cellchat <- createCellChat(
    object = data_input,
    meta = meta_input,
    group.by = GROUP_COL
  )

  cellchat@DB <- DB_USE

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
    type = "triMean",
    population.size = FALSE
  )

  cellchat <- filterCommunication(
    cellchat,
    min.cells = MIN_CELLS
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

  lr <- extract_lr_table(
    cellchat = cellchat,
    sample_name = sample_name,
    condition_name = condition_name,
    sources = SENDERS,
    targets = RECEIVERS
  )

  if (
    nrow(
      lr
    )
  ) {
    lr_tables[[
      sample_name
    ]] <- lr
  } else {
    warning(
      "No Mphi->HSC LR interactions extracted for ",
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
        "_v6.2.1.rds"
      )
    ),
    compress = FALSE
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
    lr_tables
  )
) {
  stop(
    "No LR tables were generated."
  )
}

lr_all <- bind_rows(
  lr_tables
)

write.csv(
  lr_all,
  file.path(
    TAB_OUT,
    "03_Mphi5_to_HSC3_LR_all_samples_v6.2.1.csv"
  ),
  row.names = FALSE
)

msg(
  "Total sample-level LR rows: ",
  nrow(
    lr_all
  )
)


# ==============================================================================
# 9. Standardize LR columns for comparison
# ==============================================================================

required_lr_cols <- c(
  "source",
  "target",
  "prob",
  "pathway_name",
  "interaction_key",
  "sample",
  "condition"
)

missing_lr_cols <- setdiff(
  required_lr_cols,
  colnames(
    lr_all
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
      lr_all
    )
) {
  lr_all$ligand <-
    NA_character_
}

if (
  !"receptor" %in%
    colnames(
      lr_all
    )
) {
  lr_all$receptor <-
    NA_character_
}


# ==============================================================================
# 10. LR replicate-aware Sham vs Tx summary
# ==============================================================================

lr_key_cols <- c(
  "source",
  "target",
  "pathway_name",
  "interaction_key",
  "ligand",
  "receptor"
)

lr_compare <- summarize_four_samples(
  df = lr_all,
  key_cols = lr_key_cols,
  value_col = "prob"
) %>%
  arrange(
    desc(
      pairwise_direction_consistency
    ),
    desc(
      abs(
        log2_Tx_vs_Sham
      )
    ),
    desc(
      Tx_mean +
        Sham_mean
    )
  )

write.csv(
  lr_compare,
  file.path(
    TAB_OUT,
    "04_LR_Sham_vs_Tx_replicate_summary_v6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Pathway sample-level summary
# ==============================================================================

pathway_sample <- lr_all %>%
  group_by(
    sample,
    condition,
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

    .groups = "drop"
  )

write.csv(
  pathway_sample,
  file.path(
    TAB_OUT,
    "05_pathway_strength_by_sample_v6.2.1.csv"
  ),
  row.names = FALSE
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
    "source",
    "target",
    "pathway_name"
  ),
  value_col = "prob"
) %>%
  mutate(
    priority_pathway =
      grepl(
        PRIORITY_PATHWAY_REGEX,
        pathway_name,
        ignore.case = TRUE
      )
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
        log2_Tx_vs_Sham
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
    "06_pathway_Sham_vs_Tx_replicate_summary_v6.2.1.csv"
  ),
  row.names = FALSE
)

priority_pathways <- pathway_compare %>%
  filter(
    priority_pathway
  )

write.csv(
  priority_pathways,
  file.path(
    TAB_OUT,
    "07_priority_pathways_Sham_vs_Tx_v6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Primary biological focus:
#     ECM-associated inflammatory-Mphi -> HSC3
# ==============================================================================

ecm_inflammatory_focus <- lr_compare %>%
  filter(
    source ==
      "ECM-associated inflammatory-Mphi",
    target %in%
      RECEIVERS
  ) %>%
  arrange(
    target,
    desc(
      pairwise_direction_consistency
    ),
    log2_Tx_vs_Sham
  )

write.csv(
  ecm_inflammatory_focus,
  file.path(
    TAB_OUT,
    "08_ECM_inflammatory_Mphi_to_HSC3_LR_focus_v6.2.1.csv"
  ),
  row.names = FALSE
)

ecm_inflammatory_pathway_focus <-
  pathway_compare %>%
  filter(
    source ==
      "ECM-associated inflammatory-Mphi",
    target %in%
      RECEIVERS
  ) %>%
  arrange(
    target,
    desc(
      priority_pathway
    ),
    desc(
      pairwise_direction_consistency
    ),
    log2_Tx_vs_Sham
  )

write.csv(
  ecm_inflammatory_pathway_focus,
  file.path(
    TAB_OUT,
    "09_ECM_inflammatory_Mphi_to_HSC3_pathway_focus_v6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Secondary biological focus:
#     Anti-inflammatory / Repair-Resolution -> HSC3
# ==============================================================================

resolution_focus <- lr_compare %>%
  filter(
    source %in%
      c(
        "Anti-inflammatory-Mphi",
        "Repair/Resolution-Mphi"
      ),
    target %in%
      RECEIVERS
  ) %>%
  arrange(
    source,
    target,
    desc(
      pairwise_direction_consistency
    ),
    desc(
      log2_Tx_vs_Sham
    )
  )

write.csv(
  resolution_focus,
  file.path(
    TAB_OUT,
    "10_AntiInflammatory_Repair_Mphi_to_HSC3_LR_focus_v6.2.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Figure: pathway change heatmap
# ==============================================================================

plot_pathways <- pathway_compare %>%
  filter(
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
        na.rm = TRUE
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
    n = 30
  ) %>%
  pull(
    pathway_name
  )

heat_df <- pathway_compare %>%
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
      linewidth = 0.20
    ) +
    scale_fill_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(
        -5,
        5
      ),
      oob = scales::squish
    ) +
    labs(
      title =
        "Mphi -> HSC pathway change | Tx vs Sham",
      subtitle =
        "Red = higher in Tx; blue = lower in Tx; sample-wise CellChat",
      x = NULL,
      y = NULL,
      fill =
        "log2 Tx/Sham"
    ) +
    theme_classic(
      base_size = 7
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 60,
          hjust = 1
        ),
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )

  save_pdf(
    p_heat,
    file.path(
      FIG_OUT,
      "01_Mphi5_to_HSC3_pathway_Tx_vs_Sham_heatmap_v6.2.1.pdf"
    ),
    18,
    10
  )
}


# ==============================================================================
# 16. Figure: ECM-inflammatory Mphi -> HSC3 top LR
# ==============================================================================

ecm_plot_df <- ecm_inflammatory_focus %>%
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
      pairwise_direction_consistency
    ),
    desc(
      max_strength
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
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
    ecm_plot_df
  )
) {

  p_ecm <- ggplot(
    ecm_plot_df,
    aes(
      x =
        log2_Tx_vs_Sham,
      y =
        reorder(
          interaction_label,
          log2_Tx_vs_Sham
        ),
      size =
        max_strength
    )
  ) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +
    facet_wrap(
      ~ target,
      scales = "free_y",
      ncol = 1
    ) +
    labs(
      title =
        "ECM-associated inflammatory-Mphi -> HSC",
      subtitle =
        "Top replicate-supported CellChat LR interactions",
      x =
        "log2 Tx/Sham communication probability",
      y = NULL,
      size =
        "Max mean strength"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )

  save_pdf(
    p_ecm,
    file.path(
      FIG_OUT,
      "02_ECM_inflammatory_Mphi_to_HSC3_top_LR_v6.2.1.pdf"
    ),
    10,
    12
  )
}


# ==============================================================================
# 17. Figure: sample-level priority pathway heatmap
# ==============================================================================

priority_sample <- pathway_sample %>%
  filter(
    grepl(
      PRIORITY_PATHWAY_REGEX,
      pathway_name,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    sender_receiver_pathway =
      paste0(
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
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      max_strength
    )
  ) %>%
  slice_head(
    n = 40
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
        levels = SAMPLES
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
      linewidth = 0.20
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
        "Priority Mphi -> HSC pathways by biological sample",
      x = NULL,
      y = NULL,
      fill =
        "CellChat strength"
    ) +
    theme_classic(
      base_size = 7
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        )
    )

  save_pdf(
    p_sample_heat,
    file.path(
      FIG_OUT,
      "03_priority_pathways_by_sample_heatmap_v6.2.1.pdf"
    ),
    9,
    12
  )
}


# ==============================================================================
# 18. Save merged CellChat collection
# ==============================================================================

saveRDS(
  cellchat_objects,
  file.path(
    RDS_OUT,
    "CellChat_samplewise_list_v6.2.1.rds"
  ),
  compress = FALSE
)

merged_cellchat <- tryCatch(
  mergeCellChat(
    cellchat_objects,
    add.names = names(
      cellchat_objects
    ),
    cell.prefix = TRUE
  ),
  error = function(e) {

    warning(
      "mergeCellChat failed: ",
      conditionMessage(e),
      ". Sample-wise objects and all summary tables are still valid."
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
      "CellChat_merged_4samples_v6.2.1.rds"
    ),
    compress = FALSE
  )
}


# ==============================================================================
# 19. Analysis metadata
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
    "senders",
    "receivers"
  ),
  value = c(
    "v6.2.1",
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
      collapse = ","
    ),
    paste(
      SENDERS,
      collapse = " | "
    ),
    paste(
      RECEIVERS,
      collapse = " | "
    )
  )
)

write.csv(
  analysis_metadata,
  file.path(
    LOG_OUT,
    "analysis_metadata_v6.2.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.2.1.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
