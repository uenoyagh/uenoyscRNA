#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6200)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# 5 macrophage subtypes x 3 HSC states interaction-ready object
#
# Version: v6.2.0
#
# BASE CHECKPOINTS
#   Mphi:
#     mouse-mash-mphi-staining-markers-v5.9.3
#     commit 1381eeae98099a835745313be30c7c6bbc7d1674
#
#   HSC:
#     HSC refined annotation v6.1.2
#
# INPUT OBJECTS
#   1) Frozen whole-cell parent v5.1.1
#   2) Clean-B macrophage FINAL v4.14.5
#   3) Clean HSC 3-state/2-state v6.1.2
#
# PURPOSE
#   Create a single whole-cell-derived Seurat object containing:
#
#   Sender macrophages:
#     1) Anti-inflammatory-Mphi
#     2) Inflammatory-Mphi
#     3) ECM-associated inflammatory-Mphi
#     4) Repair/Resolution-Mphi
#     5) Lipid-associated/TREM2-Mphi
#
#   Receiver HSC states:
#     6) qHSC
#     7) ECM-activated HSC
#     8) Contractile HSC
#
# SAMPLES
#   Sham1
#   Sham20
#   Tx17
#   Tx5
#
# IMPORTANT DESIGN
#   - No reclustering
#   - No reintegration
#   - No new UMAP
#   - All interaction cells are selected from the frozen whole-cell parent.
#   - Mphi and HSC labels are transferred by exact cell-barcode matching.
#   - HSC Review_boundary / Excluded_nonHSC cells are not included because
#     input HSC RDS is already the primary clean v6.1.2 object.
#   - The HSC 2-state annotation (qHSC/aHSC) is preserved as metadata.
#   - This script prepares and audits the object; it does NOT run CellChat.
#
# FUTURE WHOLE-CELL DISPLAY
#   A barcode map for the 3 HSC states is exported so that the full whole-cell
#   Layer1 display can later replace HSC_Mesenchymal with:
#     qHSC / ECM-activated HSC / Contractile HSC
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

find_first_existing <- function(paths) {
  hit <- paths[
    file.exists(paths)
  ]

  if (!length(hit)) {
    stop(
      "None of the expected input files exists:\n",
      paste(
        paths,
        collapse = "\n"
      )
    )
  }

  hit[[1]]
}

resolve_first <- function(
  x,
  candidates,
  what
) {
  hit <- candidates[
    candidates %in% x
  ]

  if (!length(hit)) {
    stop(
      "Could not resolve ",
      what,
      ". Candidates: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  hit[[1]]
}

canonical_sample <- function(x) {
  x <- as.character(x)

  dplyr::case_when(
    grepl("^Sham1$", x, ignore.case = TRUE) ~ "Sham1",
    grepl("^Sham20$", x, ignore.case = TRUE) ~ "Sham20",
    grepl("^Tx17$", x, ignore.case = TRUE) ~ "Tx17",
    grepl("^Tx5$", x, ignore.case = TRUE) ~ "Tx5",
    grepl("Sham.?1", x, ignore.case = TRUE) ~ "Sham1",
    grepl("Sham.?20", x, ignore.case = TRUE) ~ "Sham20",
    grepl("Tx.?17", x, ignore.case = TRUE) ~ "Tx17",
    grepl("Tx.?5", x, ignore.case = TRUE) ~ "Tx5",
    TRUE ~ NA_character_
  )
}

canonical_condition <- function(sample) {
  sample <- as.character(sample)

  dplyr::case_when(
    grepl("^Sham", sample, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", sample, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

MPHI_RDS <- find_first_existing(
  c(
    file.path(
      ROOT,
      "Mouse_MASH_Mphi_RDS",
      "Mphi_Res2_CleanB_FINAL_v4.14.5",
      "RDS",
      "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
    ),
    file.path(
      ROOT,
      "Mouse_MASH_Mphi_RDS",
      "Mphi_Res2_CleanB_FINAL_v4.14.5",
      "RDS",
      "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
    )
  )
)

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2",
  "RDS",
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_interaction_ready_v6.2.0"
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
# 3. Frozen metadata definitions
# ==============================================================================

WHOLE_LAYER1_COL <-
  "wholecell_layer1_FINAL_v511"

MPHI_CLASS_COL <-
  "macrophage_class_Res2_FINAL_v4145_char"

HSC_STATE3_COL <-
  "HSC_state3_v612"

HSC_STATE2_COL <-
  "HSC_state2_v612"

WHOLE_SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "sample_FIXED2",
  "sample",
  "orig.ident"
)

MPHI_SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "sample_FIXED2",
  "sample",
  "orig.ident"
)

HSC_SAMPLE_CANDIDATES <- c(
  "sample_hsc_v610",
  "sample_for_annotation",
  "sample_FIXED2",
  "sample",
  "orig.ident"
)

WHOLE_UMAP <- "umapRPCA"

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

INTERACTION_LEVELS <- c(
  MPHI5,
  HSC3
)

SAMPLE_LEVELS <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

CONDITION_LEVELS <- c(
  "Sham",
  "Tx"
)


# ==============================================================================
# 4. Load
# ==============================================================================

if (
  !file.exists(
    WHOLE_RDS
  )
) {
  stop(
    "Whole-cell frozen RDS not found: ",
    WHOLE_RDS
  )
}

if (
  !file.exists(
    HSC_RDS
  )
) {
  stop(
    "Clean HSC v6.1.2 RDS not found: ",
    HSC_RDS
  )
}

msg(
  "Loading whole-cell frozen parent..."
)

whole <- readRDS(
  WHOLE_RDS
)

msg(
  "Loading Clean-B macrophage FINAL..."
)

mphi <- readRDS(
  MPHI_RDS
)

msg(
  "Loading clean HSC v6.1.2..."
)

hsc <- readRDS(
  HSC_RDS
)

for (
  object_name in c(
    "whole",
    "mphi",
    "hsc"
  )
) {

  object_now <- get(
    object_name
  )

  if (
    !"RNA" %in%
      Assays(
        object_now
      )
  ) {
    stop(
      "RNA assay missing from ",
      object_name
    )
  }
}

DefaultAssay(
  whole
) <- "RNA"

DefaultAssay(
  mphi
) <- "RNA"

DefaultAssay(
  hsc
) <- "RNA"

if (
  !WHOLE_LAYER1_COL %in%
    colnames(
      whole@meta.data
    )
) {
  stop(
    "Missing whole-cell Layer1 column: ",
    WHOLE_LAYER1_COL
  )
}

if (
  !MPHI_CLASS_COL %in%
    colnames(
      mphi@meta.data
    )
) {
  stop(
    "Missing macrophage subtype column: ",
    MPHI_CLASS_COL
  )
}

for (
  col in c(
    HSC_STATE3_COL,
    HSC_STATE2_COL
  )
) {
  if (
    !col %in%
      colnames(
        hsc@meta.data
      )
  ) {
    stop(
      "Missing HSC annotation column: ",
      col
    )
  }
}

if (
  !WHOLE_UMAP %in%
    Reductions(
      whole
    )
) {
  stop(
    "Frozen whole-cell UMAP not found: ",
    WHOLE_UMAP
  )
}

WHOLE_SAMPLE_COL <- resolve_first(
  colnames(
    whole@meta.data
  ),
  WHOLE_SAMPLE_CANDIDATES,
  "whole-cell sample column"
)

MPHI_SAMPLE_COL <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  MPHI_SAMPLE_CANDIDATES,
  "macrophage sample column"
)

HSC_SAMPLE_COL <- resolve_first(
  colnames(
    hsc@meta.data
  ),
  HSC_SAMPLE_CANDIDATES,
  "HSC sample column"
)

msg(
  "Whole sample source: ",
  WHOLE_SAMPLE_COL
)

msg(
  "Mphi sample source: ",
  MPHI_SAMPLE_COL
)

msg(
  "HSC sample source: ",
  HSC_SAMPLE_COL
)


# ==============================================================================
# 5. Canonical sample labels
# ==============================================================================

whole$sample_interaction_v620 <-
  canonical_sample(
    whole@meta.data[[
      WHOLE_SAMPLE_COL
    ]]
  )

whole$condition_interaction_v620 <-
  canonical_condition(
    whole$sample_interaction_v620
  )

mphi$sample_interaction_v620 <-
  canonical_sample(
    mphi@meta.data[[
      MPHI_SAMPLE_COL
    ]]
  )

mphi$condition_interaction_v620 <-
  canonical_condition(
    mphi$sample_interaction_v620
  )

hsc$sample_interaction_v620 <-
  canonical_sample(
    hsc@meta.data[[
      HSC_SAMPLE_COL
    ]]
  )

hsc$condition_interaction_v620 <-
  canonical_condition(
    hsc$sample_interaction_v620
  )

msg(
  "Whole sample counts:"
)

print(
  table(
    whole$sample_interaction_v620,
    useNA = "ifany"
  )
)

msg(
  "Mphi sample counts:"
)

print(
  table(
    mphi$sample_interaction_v620,
    useNA = "ifany"
  )
)

msg(
  "HSC sample counts:"
)

print(
  table(
    hsc$sample_interaction_v620,
    useNA = "ifany"
  )
)


# ==============================================================================
# 6. Build transfer maps
# ==============================================================================

mphi_map <- tibble(
  cell = colnames(
    mphi
  ),
  interaction_celltype =
    as.character(
      mphi@meta.data[[
        MPHI_CLASS_COL
      ]]
    ),
  sample =
    mphi$sample_interaction_v620,
  condition =
    mphi$condition_interaction_v620,
  source =
    "CleanB_Mphi_v4.14.5"
) %>%
  filter(
    sample %in%
      SAMPLE_LEVELS,
    interaction_celltype %in%
      MPHI5
  )

hsc_map <- tibble(
  cell = colnames(
    hsc
  ),
  interaction_celltype =
    as.character(
      hsc@meta.data[[
        HSC_STATE3_COL
      ]]
    ),
  HSC_state2 =
    as.character(
      hsc@meta.data[[
        HSC_STATE2_COL
      ]]
    ),
  sample =
    hsc$sample_interaction_v620,
  condition =
    hsc$condition_interaction_v620,
  source =
    "HSC_clean_v6.1.2"
) %>%
  filter(
    sample %in%
      SAMPLE_LEVELS,
    interaction_celltype %in%
      HSC3
  )

if (
  anyDuplicated(
    mphi_map$cell
  )
) {
  stop(
    "Duplicated barcodes in Mphi transfer map."
  )
}

if (
  anyDuplicated(
    hsc_map$cell
  )
) {
  stop(
    "Duplicated barcodes in HSC transfer map."
  )
}

overlap_cells <- intersect(
  mphi_map$cell,
  hsc_map$cell
)

if (
  length(
    overlap_cells
  ) > 0
) {
  stop(
    "Mphi and HSC transfer maps overlap for ",
    length(
      overlap_cells
    ),
    " cells. Review annotation sources before proceeding."
  )
}


# ==============================================================================
# 7. Exact barcode mapping audit
# ==============================================================================

whole_cells <- colnames(
  whole
)

mphi_map <- mphi_map %>%
  mutate(
    exact_match_in_whole =
      cell %in%
        whole_cells
  )

hsc_map <- hsc_map %>%
  mutate(
    exact_match_in_whole =
      cell %in%
        whole_cells
  )

mapping_global <- tibble(
  source = c(
    "Mphi5",
    "HSC3"
  ),
  n_input_cells = c(
    nrow(
      mphi_map
    ),
    nrow(
      hsc_map
    )
  ),
  n_exact_matches = c(
    sum(
      mphi_map$exact_match_in_whole
    ),
    sum(
      hsc_map$exact_match_in_whole
    )
  )
) %>%
  mutate(
    exact_match_fraction =
      n_exact_matches /
      n_input_cells
  )

write.csv(
  mapping_global,
  file.path(
    TAB_OUT,
    "01_barcode_mapping_global_audit_v6.2.0.csv"
  ),
  row.names = FALSE
)

print(
  mapping_global
)

if (
  any(
    !is.finite(
      mapping_global$exact_match_fraction
    )
  ) ||
  any(
    mapping_global$exact_match_fraction <
      0.95
  )
) {
  stop(
    "Exact barcode mapping below 95% for Mphi and/or HSC. ",
    "Review 01_barcode_mapping_global_audit_v6.2.0.csv."
  )
}

write.csv(
  mphi_map,
  file.path(
    TAB_OUT,
    "02_Mphi5_barcode_transfer_map_v6.2.0.csv"
  ),
  row.names = FALSE
)

write.csv(
  hsc_map,
  file.path(
    TAB_OUT,
    "03_HSC3_barcode_transfer_map_v6.2.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Sample consistency audit
# ==============================================================================

mphi_idx <- match(
  mphi_map$cell,
  whole_cells
)

hsc_idx <- match(
  hsc_map$cell,
  whole_cells
)

mphi_whole_sample <-
  whole$sample_interaction_v620[
    mphi_idx
  ]

hsc_whole_sample <-
  whole$sample_interaction_v620[
    hsc_idx
  ]

sample_audit <- bind_rows(

  tibble(
    source = "Mphi5",
    cell = mphi_map$cell,
    source_sample = mphi_map$sample,
    whole_sample =
      mphi_whole_sample
  ),

  tibble(
    source = "HSC3",
    cell = hsc_map$cell,
    source_sample = hsc_map$sample,
    whole_sample =
      hsc_whole_sample
  )
) %>%
  mutate(
    same_sample =
      source_sample ==
      whole_sample
  )

write.csv(
  sample_audit,
  file.path(
    TAB_OUT,
    "04_sample_consistency_audit_v6.2.0.csv"
  ),
  row.names = FALSE
)

if (
  any(
    sample_audit$same_sample %in%
      FALSE,
    na.rm = TRUE
  )
) {
  bad_n <- sum(
    sample_audit$same_sample %in%
      FALSE,
    na.rm = TRUE
  )

  stop(
    "Sample mismatch between source objects and whole-cell parent for ",
    bad_n,
    " cells."
  )
}


# ==============================================================================
# 9. Transfer labels to frozen whole-cell parent
# ==============================================================================

whole$interaction_celltype_v620 <-
  NA_character_

whole$interaction_role_v620 <-
  NA_character_

whole$HSC_state2_v620 <-
  NA_character_

# Mphi
whole$interaction_celltype_v620[
  mphi_idx
] <-
  mphi_map$interaction_celltype

whole$interaction_role_v620[
  mphi_idx
] <-
  "Sender_Mphi"

# HSC
whole$interaction_celltype_v620[
  hsc_idx
] <-
  hsc_map$interaction_celltype

whole$interaction_role_v620[
  hsc_idx
] <-
  "Receiver_HSC"

whole$HSC_state2_v620[
  hsc_idx
] <-
  hsc_map$HSC_state2


# ==============================================================================
# 10. Build interaction-ready object
# ==============================================================================

keep_cells <- colnames(
  whole
)[
  !is.na(
    whole$interaction_celltype_v620
  ) &
    whole$sample_interaction_v620 %in%
      SAMPLE_LEVELS
]

interaction <- subset(
  whole,
  cells = keep_cells
)

interaction$interaction_celltype_v620 <-
  factor(
    interaction$interaction_celltype_v620,
    levels = INTERACTION_LEVELS
  )

interaction$sample_interaction_v620 <-
  factor(
    interaction$sample_interaction_v620,
    levels = SAMPLE_LEVELS
  )

interaction$condition_interaction_v620 <-
  factor(
    interaction$condition_interaction_v620,
    levels = CONDITION_LEVELS
  )

interaction$interaction_role_v620 <-
  factor(
    interaction$interaction_role_v620,
    levels = c(
      "Sender_Mphi",
      "Receiver_HSC"
    )
  )

interaction$HSC_state2_v620 <-
  factor(
    interaction$HSC_state2_v620,
    levels = c(
      "qHSC",
      "aHSC"
    )
  )

Idents(
  interaction
) <-
  "interaction_celltype_v620"


# ==============================================================================
# 11. Count/fraction tables
# ==============================================================================

interaction_meta <- interaction@meta.data %>%
  as_tibble(
    rownames = "cell"
  )

count_table <- interaction_meta %>%
  count(
    sample_interaction_v620,
    condition_interaction_v620,
    interaction_celltype_v620,
    interaction_role_v620,
    name = "n_cells"
  ) %>%
  complete(
    sample_interaction_v620 =
      factor(
        SAMPLE_LEVELS,
        levels = SAMPLE_LEVELS
      ),
    interaction_celltype_v620 =
      factor(
        INTERACTION_LEVELS,
        levels = INTERACTION_LEVELS
      ),
    fill = list(
      n_cells = 0
    )
  ) %>%
  group_by(
    sample_interaction_v620
  ) %>%
  mutate(
    fraction_within_interaction_object =
      n_cells /
      sum(
        n_cells
      )
  ) %>%
  ungroup()

write.csv(
  count_table,
  file.path(
    TAB_OUT,
    "05_sample_by_interaction_celltype_counts_v6.2.0.csv"
  ),
  row.names = FALSE
)

msg(
  "Sample x interaction-celltype counts:"
)

print(
  count_table
)

# Separate receiver table.
hsc_count_table <- interaction_meta %>%
  filter(
    interaction_role_v620 ==
      "Receiver_HSC"
  ) %>%
  count(
    sample_interaction_v620,
    condition_interaction_v620,
    interaction_celltype_v620,
    HSC_state2_v620,
    name = "n_cells"
  )

write.csv(
  hsc_count_table,
  file.path(
    TAB_OUT,
    "06_HSC3_receiver_counts_by_sample_v6.2.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Minimum-cell audit
# ==============================================================================

minimum_cell_audit <- count_table %>%
  mutate(
    below_10 =
      n_cells < 10,
    below_20 =
      n_cells < 20,
    below_30 =
      n_cells < 30,
    below_50 =
      n_cells < 50
  )

write.csv(
  minimum_cell_audit,
  file.path(
    TAB_OUT,
    "07_CellChat_minimum_cell_audit_v6.2.0.csv"
  ),
  row.names = FALSE
)

if (
  any(
    minimum_cell_audit$n_cells <
      10
  )
) {
  warning(
    "At least one sample x interaction-celltype group has <10 cells. ",
    "Review 07_CellChat_minimum_cell_audit_v6.2.0.csv before CellChat."
  )
}


# ==============================================================================
# 13. Full metadata / HSC3 future whole-cell mapping table
# ==============================================================================

write.csv(
  interaction_meta,
  file.path(
    TAB_OUT,
    "08_interaction_ready_cell_metadata_v6.2.0.csv"
  ),
  row.names = FALSE
)

# This compact table is intentionally saved for the future whole-cell
# Layer1 update where HSC_Mesenchymal will be replaced by three HSC states.
wholecell_hsc3_map <- hsc_map %>%
  filter(
    exact_match_in_whole
  ) %>%
  select(
    cell,
    sample,
    condition,
    interaction_celltype,
    HSC_state2
  ) %>%
  rename(
    HSC_state3_v620 =
      interaction_celltype,
    HSC_state2_v620 =
      HSC_state2
  )

write.csv(
  wholecell_hsc3_map,
  file.path(
    TAB_OUT,
    "09_wholecell_HSC3_barcode_map_for_future_Layer1_update_v6.2.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. RNA layer audit
# ==============================================================================

rna_layers <- Layers(
  interaction[["RNA"]]
)

write.csv(
  tibble(
    assay = "RNA",
    layer = rna_layers
  ),
  file.path(
    TAB_OUT,
    "10_RNA_layer_audit_v6.2.0.csv"
  ),
  row.names = FALSE
)

msg(
  "RNA layers: ",
  paste(
    rna_layers,
    collapse = ", "
  )
)


# ==============================================================================
# 15. UMAP QC
# ==============================================================================

PLOT_COLORS <- c(
  "Anti-inflammatory-Mphi" =
    "#00C8FF",
  "Inflammatory-Mphi" =
    "#FF3131",
  "ECM-associated inflammatory-Mphi" =
    "#FF8C00",
  "Repair/Resolution-Mphi" =
    "#A020F0",
  "Lipid-associated/TREM2-Mphi" =
    "#0066FF",
  "qHSC" =
    "#00BFC4",
  "ECM-activated HSC" =
    "#F39C12",
  "Contractile HSC" =
    "#E7298A"
)

p_all <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by =
    "interaction_celltype_v620",
  cols = PLOT_COLORS,
  pt.size = 0.35,
  raster = FALSE
) +
  ggtitle(
    "Interaction-ready object | 5 Mphi subtypes + 3 HSC states"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

save_pdf(
  p_all,
  file.path(
    FIG_OUT,
    "01_interaction_ready_UMAP_all_samples_v6.2.0.pdf"
  ),
  10,
  7
)

p_sample <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by =
    "interaction_celltype_v620",
  split.by =
    "sample_interaction_v620",
  cols = PLOT_COLORS,
  pt.size = 0.28,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "Interaction-ready object | biological samples"
  )

save_pdf(
  p_sample,
  file.path(
    FIG_OUT,
    "02_interaction_ready_UMAP_by_sample_v6.2.0.pdf"
  ),
  14,
  10
)

p_condition <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by =
    "interaction_celltype_v620",
  split.by =
    "condition_interaction_v620",
  cols = PLOT_COLORS,
  pt.size = 0.28,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(
    title =
      "Interaction-ready object | Sham vs Tx"
  )

save_pdf(
  p_condition,
  file.path(
    FIG_OUT,
    "03_interaction_ready_UMAP_Sham_vs_Tx_v6.2.0.pdf"
  ),
  14,
  6
)

# HSC-only view using frozen whole-cell coordinates.
hsc_interaction_cells <- colnames(
  interaction
)[
  interaction$interaction_role_v620 ==
    "Receiver_HSC"
]

hsc_view <- subset(
  interaction,
  cells = hsc_interaction_cells
)

p_hsc <- DimPlot(
  hsc_view,
  reduction = WHOLE_UMAP,
  group.by =
    "interaction_celltype_v620",
  cols = PLOT_COLORS[
    HSC3
  ],
  pt.size = 0.45,
  raster = FALSE
) +
  ggtitle(
    "HSC receiver states | frozen whole-cell UMAP coordinates"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_hsc,
  file.path(
    FIG_OUT,
    "04_HSC3_receiver_UMAP_wholecell_coordinates_v6.2.0.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 16. Cell-count heatmap
# ==============================================================================

p_count_heat <- ggplot(
  count_table,
  aes(
    x =
      interaction_celltype_v620,
    y =
      sample_interaction_v620,
    fill =
      log10(
        n_cells + 1
      )
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = n_cells
    ),
    size = 3
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
      "Cells available for per-sample Mphi-to-HSC interaction analysis",
    x = NULL,
    y = NULL,
    fill =
      "log10(n+1)"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

save_pdf(
  p_count_heat,
  file.path(
    FIG_OUT,
    "05_sample_by_Mphi5_HSC3_cellcount_heatmap_v6.2.0.pdf"
  ),
  13,
  5
)


# ==============================================================================
# 17. Save interaction-ready RDS
# ==============================================================================

INTERACTION_RDS <- file.path(
  RDS_OUT,
  "Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds"
)

saveRDS(
  interaction,
  INTERACTION_RDS,
  compress = FALSE
)

msg(
  "Saved interaction-ready RDS: ",
  INTERACTION_RDS
)


# ==============================================================================
# 18. Final audit
# ==============================================================================

final_audit <- tibble(
  metric = c(
    "wholecell_cells",
    "Mphi5_input_cells",
    "HSC3_input_cells",
    "Mphi5_exact_match_fraction",
    "HSC3_exact_match_fraction",
    "interaction_ready_cells",
    "whole_sample_column",
    "Mphi_sample_column",
    "HSC_sample_column",
    "Mphi_annotation_column",
    "HSC_3state_column",
    "HSC_2state_column",
    "wholecell_UMAP"
  ),
  value = c(
    as.character(
      ncol(
        whole
      )
    ),
    as.character(
      nrow(
        mphi_map
      )
    ),
    as.character(
      nrow(
        hsc_map
      )
    ),
    as.character(
      mapping_global$exact_match_fraction[
        mapping_global$source ==
          "Mphi5"
      ]
    ),
    as.character(
      mapping_global$exact_match_fraction[
        mapping_global$source ==
          "HSC3"
      ]
    ),
    as.character(
      ncol(
        interaction
      )
    ),
    WHOLE_SAMPLE_COL,
    MPHI_SAMPLE_COL,
    HSC_SAMPLE_COL,
    MPHI_CLASS_COL,
    HSC_STATE3_COL,
    HSC_STATE2_COL,
    WHOLE_UMAP
  )
)

write.csv(
  final_audit,
  file.path(
    LOG_OUT,
    "final_audit_v6.2.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.2.0.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
