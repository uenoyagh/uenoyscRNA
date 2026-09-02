#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(7031)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Sham vs Tx balanced BROAD R8 whole-cell UMAP
#
# Version: v7.0.3.1
#
# FIX FROM v7.0.3
#   v7.0.3 incorrectly expected a saved "refined whole-cell Seurat RDS".
#   v7.0.1 was run with SAVE_UPDATED_WHOLE_RDS = FALSE, so that RDS does not
#   exist by design.
#
#   v7.0.3.1 therefore uses the actual validated sources:
#
#     1) UMAP coordinates / whole-cell object:
#        WholeCell_Layer1_ParentFreeze_v5.1.1
#
#     2) Current refined annotation:
#        WholeCell_RefinedAnnotation_UMAP_v7.0.1/Tables/
#        05_wholecell_refined_metadata_v7.0.1.csv
#
#     3) R8 refined palette manifest:
#        WholeCell_RefinedAnnotation_UMAP_v7.0.1/Tables/
#        08_R8_refined_palette_manifest_v7.0.1.csv
#
#     4) Biological sample identity:
#        recovered directly from exact cell-name prefixes
#
# PURPOSE
#   Produce a publication/QC-friendly BROAD whole-liver comparison:
#
#       Sham vs Tx
#
#   after equal downsampling of:
#       Sham1 / Sham20 / Tx17 / Tx5
#
# IMPORTANT
#   - NO reintegration
#   - NO reclustering
#   - NO new UMAP
#   - SAME frozen umapRPCA coordinates
#   - CURRENT v7.0.1 refined lineage ownership is collapsed to broad lineage
#   - Balanced display is visualization-only, not inferential statistics
# ==============================================================================


# ==============================================================================
# 1. Paths
# ==============================================================================

VERSION <- "v7.0.3.1"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

V701_DIR <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_RefinedAnnotation_UMAP_v7.0.1"
)

REFINED_METADATA_CSV <- file.path(
  V701_DIR,
  "Tables",
  "05_wholecell_refined_metadata_v7.0.1.csv"
)

R8_REFINED_PALETTE_CSV <- file.path(
  V701_DIR,
  "Tables",
  "08_R8_refined_palette_manifest_v7.0.1.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_RefinedAnnotation_UMAP_v7.0.3.1_balanced_broad_R8"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

for (d in c(OUT, FIG_OUT, TAB_OUT, LOG_OUT)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

WHOLE_UMAP <- "umapRPCA"

ALL_SAMPLE_LEVELS <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

TARGET_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

TARGET_CONDITIONS <- c(
  "Sham",
  "Tx"
)

EXPECTED_CONDITION_COUNTS <- c(
  "STD" = 12211L,
  "CDAHFD" = 11543L,
  "Sham" = 38382L,
  "Tx" = 42452L
)


# ==============================================================================
# 2. Helpers
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


save_png <- function(
  p,
  file,
  width,
  height,
  dpi = 400
) {

  ggsave(
    filename = file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )
}


sample_from_cellname <- function(cells) {

  cells <- as.character(cells)

  out <- rep(
    NA_character_,
    length(cells)
  )

  out[
    grepl("^STD_rep1_", cells)
  ] <- "STD_rep1"

  out[
    grepl("^CDHFD_rep1_", cells)
  ] <- "CDHFD_rep1"

  out[
    grepl("^Sham20_", cells)
  ] <- "Sham20"

  out[
    grepl("^Sham1_", cells)
  ] <- "Sham1"

  out[
    grepl("^Tx17_", cells)
  ] <- "Tx17"

  out[
    grepl("^Tx5_", cells)
  ] <- "Tx5"

  out
}


sample_to_condition <- function(x) {

  x <- as.character(x)

  out <- rep(
    NA_character_,
    length(x)
  )

  out[
    x == "STD_rep1"
  ] <- "STD"

  out[
    x == "CDHFD_rep1"
  ] <- "CDAHFD"

  out[
    x %in%
      c(
        "Sham1",
        "Sham20"
      )
  ] <- "Sham"

  out[
    x %in%
      c(
        "Tx17",
        "Tx5"
      )
  ] <- "Tx"

  out
}


collapse_refined_to_broad <- function(
  refined,
  broad_original
) {

  refined <- as.character(refined)
  broad_original <- as.character(broad_original)

  # Start from the current broad label carried in the v7.0.1 metadata.
  out <- broad_original

  # Then explicitly impose the CURRENT refined lineage ownership.
  # This is important for cells re-owned by the dedicated HSC/Hep/LSEC/
  # Cholangiocyte/Monocyte objects in v7.0.1.
  out[
    grepl("^Mphi \\|", refined)
  ] <- "Kupffer/Macrophage"

  out[
    grepl("^Mono \\|", refined)
  ] <- "Monocyte"

  out[
    grepl("^HSC \\|", refined)
  ] <- "HSC/Mesenchymal"

  out[
    grepl("^Hep \\|", refined)
  ] <- "Hepatocyte"

  out[
    grepl("^LSEC \\|", refined)
  ] <- "LSEC"

  out[
    grepl("^Chol \\|", refined)
  ] <- "Cholangiocyte"

  # Standardize broad names for display.
  out[
    out %in%
      c(
        "Kupffer_Macrophage",
        "Kupffer/Macrophage"
      )
  ] <- "Kupffer/Macrophage"

  out[
    out %in%
      c(
        "HSC_Mesenchymal",
        "HSC/Mesenchymal"
      )
  ] <- "HSC/Mesenchymal"

  out[
    out %in%
      c(
        "Vascular_endothelial",
        "Vascular endothelial"
      )
  ] <- "Vascular endothelial"

  out[
    out %in%
      c(
        "Dendritic_cell",
        "Dendritic cell"
      )
  ] <- "Dendritic cell"

  out[
    out %in%
      c(
        "NK_cell",
        "NK cell"
      )
  ] <- "NK cell"

  out[
    out %in%
      c(
        "T_cell",
        "T cell"
      )
  ] <- "T cell"

  out[
    out %in%
      c(
        "B_cell",
        "B cell"
      )
  ] <- "B cell"

  out[
    out %in%
      c(
        "Plasma_cell",
        "Plasma cell"
      )
  ] <- "Plasma cell"

  out[
    out %in%
      c(
        "Monocyte/Neutrophil_boundary",
        "Monocyte/Neutrophil boundary"
      )
  ] <- "Monocyte/Neutrophil boundary"

  out[
    out %in%
      c(
        "NK/T_boundary",
        "NK/T boundary"
      )
  ] <- "NK/T boundary"

  out
}


plot_theme <- function(
  base_size = 10
) {

  theme_classic(
    base_size = base_size
  ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        ),
      plot.subtitle =
        element_text(
          hjust = 0.5
        ),
      strip.background =
        element_rect(
          fill = "white",
          colour = "black"
        ),
      strip.text =
        element_text(
          face = "bold"
        ),
      axis.text =
        element_text(
          colour = "black"
        ),
      legend.title =
        element_blank(),
      legend.text =
        element_text(
          size = 8
        ),
      legend.key.height =
        grid::unit(
          0.38,
          "cm"
        ),
      legend.key.width =
        grid::unit(
          0.38,
          "cm"
        )
    )
}


# ==============================================================================
# 3. Broad R8 palette
#
# Preserves the established high-saturation whole-cell direction:
#   LSEC = vivid cyan
#   Hepatocyte = blue-green
#   Cholangiocyte = vivid green
#   HSC/Mesenchymal = vivid pink
# ==============================================================================

BROAD_LEVELS <- c(
  "Hepatocyte",
  "Cholangiocyte",
  "HSC/Mesenchymal",
  "LSEC",
  "Vascular endothelial",
  "Kupffer/Macrophage",
  "Monocyte",
  "Neutrophil",
  "Dendritic cell",
  "NK cell",
  "T cell",
  "B cell",
  "Plasma cell",
  "Mesothelial",
  "Platelet",
  "Cycling",
  "RBC",
  "Monocyte/Neutrophil boundary",
  "NK/T boundary",
  "Other"
)

BROAD_COLORS <- c(
  "Hepatocyte" = "#00A087",
  "Cholangiocyte" = "#00C853",
  "HSC/Mesenchymal" = "#FF4FA3",
  "LSEC" = "#00C8FF",
  "Vascular endothelial" = "#FF5C8A",
  "Kupffer/Macrophage" = "#F04444",
  "Monocyte" = "#2F65FF",
  "Neutrophil" = "#0066FF",
  "Dendritic cell" = "#A6D854",
  "NK cell" = "#5E3CFF",
  "T cell" = "#FF1493",
  "B cell" = "#FF6B6B",
  "Plasma cell" = "#A020F0",
  "Mesothelial" = "#00BFA5",
  "Platelet" = "#C77CFF",
  "Cycling" = "#FFB000",
  "RBC" = "#E83E9B",
  "Monocyte/Neutrophil boundary" = "#8A8A8A",
  "NK/T boundary" = "#66C2D7",
  "Other" = "#B8B8B8"
)


# ==============================================================================
# 4. Input checks
# ==============================================================================

for (
  f in c(
    WHOLE_RDS,
    REFINED_METADATA_CSV,
    R8_REFINED_PALETTE_CSV
  )
) {

  if (!file.exists(f)) {
    stop(
      "Required input not found: ",
      f
    )
  }
}


# ==============================================================================
# 5. Load frozen whole-cell coordinates
# ==============================================================================

msg(
  "Loading frozen whole-cell coordinate source..."
)

whole <- readRDS(
  WHOLE_RDS
)

if (
  !inherits(
    whole,
    "Seurat"
  )
) {
  stop(
    "WHOLE_RDS is not a Seurat object."
  )
}

if (
  !WHOLE_UMAP %in%
    Reductions(
      whole
    )
) {
  stop(
    "Required UMAP reduction not found: ",
    WHOLE_UMAP
  )
}

whole_cells <- colnames(
  whole
)

msg(
  "Whole cells: ",
  length(
    whole_cells
  )
)


# ==============================================================================
# 6. Recover sample identity from cell-name prefix
# ==============================================================================

whole$sample_v7031 <-
  sample_from_cellname(
    whole_cells
  )

whole$condition_v7031 <-
  sample_to_condition(
    whole$sample_v7031
  )

if (
  anyNA(
    whole$sample_v7031
  )
) {

  stop(
    "Unresolved sample prefixes: ",
    sum(
      is.na(
        whole$sample_v7031
      )
    )
  )
}


sample_counts <- as.data.frame(
  table(
    sample =
      factor(
        whole$sample_v7031,
        levels =
          ALL_SAMPLE_LEVELS
      )
  ),
  stringsAsFactors = FALSE
)

names(
  sample_counts
) <- c(
  "sample",
  "n_cells"
)


condition_counts <- as.data.frame(
  table(
    condition =
      factor(
        whole$condition_v7031,
        levels =
          names(
            EXPECTED_CONDITION_COUNTS
          )
      )
  ),
  stringsAsFactors = FALSE
)

names(
  condition_counts
) <- c(
  "condition",
  "n_cells"
)

condition_counts$expected_n_cells <-
  unname(
    EXPECTED_CONDITION_COUNTS[
      condition_counts$condition
    ]
  )

condition_counts$matches_expected <-
  condition_counts$n_cells ==
    condition_counts$expected_n_cells


msg(
  "Resolved sample counts:"
)

print(
  sample_counts
)

msg(
  "Condition-count validation:"
)

print(
  condition_counts
)


if (
  !all(
    condition_counts$matches_expected
  )
) {
  stop(
    "Barcode-prefix sample recovery failed validated condition totals."
  )
}


# ==============================================================================
# 7. Load current v7.0.1 refined metadata
# ==============================================================================

msg(
  "Loading v7.0.1 refined metadata..."
)

meta701 <- read.csv(
  REFINED_METADATA_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "cell",
  "broad_annotation",
  "refined_annotation"
)

missing_cols <- setdiff(
  required_cols,
  colnames(
    meta701
  )
)

if (length(missing_cols)) {

  stop(
    "v7.0.1 metadata missing columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

if (
  anyDuplicated(
    meta701$cell
  )
) {
  stop(
    "Duplicated barcodes in v7.0.1 refined metadata."
  )
}

idx_meta <- match(
  whole_cells,
  meta701$cell
)

if (
  anyNA(
    idx_meta
  )
) {
  stop(
    "v7.0.1 metadata missing whole-cell barcodes: ",
    sum(
      is.na(
        idx_meta
      )
    )
  )
}


whole$broad_original_v7031 <-
  meta701$broad_annotation[
    idx_meta
  ]

whole$refined_annotation_v7031 <-
  meta701$refined_annotation[
    idx_meta
  ]

whole$broad_current_v7031 <-
  collapse_refined_to_broad(
    refined =
      whole$refined_annotation_v7031,
    broad_original =
      whole$broad_original_v7031
  )


broad_present <- sort(
  unique(
    whole$broad_current_v7031
  )
)

unknown_broad <- setdiff(
  broad_present,
  names(
    BROAD_COLORS
  )
)

if (length(unknown_broad)) {

  stop(
    "Broad labels without R8 color mapping: ",
    paste(
      unknown_broad,
      collapse = ", "
    )
  )
}


whole$broad_current_v7031 <-
  factor(
    whole$broad_current_v7031,
    levels =
      BROAD_LEVELS
  )


# ==============================================================================
# 8. Select Sham/Tx and balance by biological sample
# ==============================================================================

target_cells <- whole_cells[
  whole$sample_v7031 %in%
    TARGET_SAMPLES
]

target <- subset(
  whole,
  cells =
    target_cells
)

target$sample_v7031 <-
  factor(
    target$sample_v7031,
    levels =
      TARGET_SAMPLES
  )

target$condition_v7031 <-
  factor(
    target$condition_v7031,
    levels =
      TARGET_CONDITIONS
  )

target$broad_current_v7031 <-
  factor(
    target$broad_current_v7031,
    levels =
      BROAD_LEVELS
  )


target_sample_counts <- as.data.frame(
  table(
    sample =
      target$sample_v7031
  ),
  stringsAsFactors = FALSE
)

names(
  target_sample_counts
) <- c(
  "sample",
  "n_cells"
)

n_per_sample <- min(
  target_sample_counts$n_cells
)

msg(
  "Balanced cells per biological sample: ",
  n_per_sample
)


set.seed(
  70310
)

balanced_cells <- unlist(
  lapply(
    TARGET_SAMPLES,
    function(s) {

      cells_now <- colnames(
        target
      )[
        target$sample_v7031 ==
          s
      ]

      sample(
        cells_now,
        size =
          n_per_sample,
        replace =
          FALSE
      )
    }
  ),
  use.names =
    FALSE
)

balanced <- subset(
  target,
  cells =
    balanced_cells
)

balanced$sample_v7031 <-
  factor(
    balanced$sample_v7031,
    levels =
      TARGET_SAMPLES
  )

balanced$condition_v7031 <-
  factor(
    balanced$condition_v7031,
    levels =
      TARGET_CONDITIONS
  )

balanced$broad_current_v7031 <-
  factor(
    balanced$broad_current_v7031,
    levels =
      BROAD_LEVELS
  )


# ==============================================================================
# 9. Audit tables
# ==============================================================================

balanced_sample_counts <- as.data.frame(
  table(
    sample =
      balanced$sample_v7031
  ),
  stringsAsFactors = FALSE
)

names(
  balanced_sample_counts
) <- c(
  "sample",
  "n_cells"
)


balanced_condition_counts <- as.data.frame(
  table(
    condition =
      balanced$condition_v7031
  ),
  stringsAsFactors = FALSE
)

names(
  balanced_condition_counts
) <- c(
  "condition",
  "n_cells"
)


broad_by_condition <- as.data.frame(
  table(
    condition =
      balanced$condition_v7031,
    broad =
      balanced$broad_current_v7031
  ),
  stringsAsFactors = FALSE
)

names(
  broad_by_condition
) <- c(
  "condition",
  "broad",
  "n_cells"
)

broad_by_condition <- broad_by_condition[
  broad_by_condition$n_cells > 0,
  ,
  drop = FALSE
]


broad_by_sample <- as.data.frame(
  table(
    sample =
      balanced$sample_v7031,
    broad =
      balanced$broad_current_v7031
  ),
  stringsAsFactors = FALSE
)

names(
  broad_by_sample
) <- c(
  "sample",
  "broad",
  "n_cells"
)

broad_by_sample <- broad_by_sample[
  broad_by_sample$n_cells > 0,
  ,
  drop = FALSE
]


write.csv(
  sample_counts,
  file.path(
    TAB_OUT,
    paste0(
      "01_all_sample_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  condition_counts,
  file.path(
    TAB_OUT,
    paste0(
      "02_condition_validation_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  balanced_sample_counts,
  file.path(
    TAB_OUT,
    paste0(
      "03_balanced_sample_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  balanced_condition_counts,
  file.path(
    TAB_OUT,
    paste0(
      "04_balanced_condition_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  broad_by_condition,
  file.path(
    TAB_OUT,
    paste0(
      "05_balanced_broad_by_condition_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  broad_by_sample,
  file.path(
    TAB_OUT,
    paste0(
      "06_balanced_broad_by_sample_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Figure 01: balanced broad Sham vs Tx
# ==============================================================================

msg(
  "Drawing balanced broad Sham vs Tx R8 UMAP..."
)

p_condition <- DimPlot(
  balanced,
  reduction =
    WHOLE_UMAP,
  group.by =
    "broad_current_v7031",
  split.by =
    "condition_v7031",
  cols =
    BROAD_COLORS,
  pt.size =
    0.42,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70311
) +
  patchwork::plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | Sham vs Tx | balanced broad R8",
    subtitle =
      paste0(
        n_per_sample,
        " cells per biological sample | ",
        2 * n_per_sample,
        " cells per condition"
      )
  )

p_condition <-
  p_condition &
  plot_theme(
    base_size = 10
  )


save_pdf(
  p_condition,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_broad_Sham_vs_Tx_BALANCED_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 16,
  height = 9
)

save_png(
  p_condition,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_broad_Sham_vs_Tx_BALANCED_R8_",
      VERSION,
      ".png"
    )
  ),
  width = 16,
  height = 9
)


# ==============================================================================
# 11. Figure 02: four balanced biological samples
# ==============================================================================

msg(
  "Drawing balanced broad biological-replicate R8 UMAP..."
)

p_sample <- DimPlot(
  balanced,
  reduction =
    WHOLE_UMAP,
  group.by =
    "broad_current_v7031",
  split.by =
    "sample_v7031",
  cols =
    BROAD_COLORS,
  pt.size =
    0.36,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70312
) +
  patchwork::plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | Sham1 / Sham20 / Tx17 / Tx5 | balanced broad R8",
    subtitle =
      paste0(
        n_per_sample,
        " cells per biological sample"
      )
  )

p_sample <-
  p_sample &
  plot_theme(
    base_size = 9
  )


save_pdf(
  p_sample,
  file.path(
    FIG_OUT,
    paste0(
      "02_wholecell_broad_Sham1_Sham20_Tx17_Tx5_BALANCED_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 16,
  height = 12
)

save_png(
  p_sample,
  file.path(
    FIG_OUT,
    paste0(
      "02_wholecell_broad_Sham1_Sham20_Tx17_Tx5_BALANCED_R8_",
      VERSION,
      ".png"
    )
  ),
  width = 16,
  height = 12
)


# ==============================================================================
# 12. Palette manifest
# ==============================================================================

palette_manifest <- data.frame(
  broad =
    names(
      BROAD_COLORS
    ),
  hex =
    unname(
      BROAD_COLORS
    ),
  present =
    names(
      BROAD_COLORS
    ) %in%
      broad_present,
  stringsAsFactors = FALSE
)

write.csv(
  palette_manifest,
  file.path(
    TAB_OUT,
    paste0(
      "07_broad_R8_palette_manifest_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. README / session info
# ==============================================================================

readme <- c(
  paste0(
    "Mouse MASH balanced broad R8 Sham vs Tx ",
    VERSION
  ),
  "============================================================",
  "",
  paste0(
    "Whole-cell coordinate source: ",
    WHOLE_RDS
  ),
  paste0(
    "Current refined metadata source: ",
    REFINED_METADATA_CSV
  ),
  "",
  "Sample identity:",
  "- derived directly from exact cell-name prefix",
  "",
  "Broad-lineage policy:",
  "- start from v7.0.1 broad annotation",
  "- override broad lineage using current v7.0.1 refined-lineage ownership",
  "- therefore Monocyte/HSC/Hep/LSEC/Chol dedicated re-ownership is retained",
  "",
  paste0(
    "Balanced cells per biological sample: ",
    n_per_sample
  ),
  paste0(
    "Balanced cells per condition: ",
    2 * n_per_sample
  ),
  "",
  "Balanced display is visualization-only.",
  "No reintegration, reclustering, or UMAP recomputation."
)

writeLines(
  readme,
  file.path(
    OUT,
    paste0(
      "README_",
      VERSION,
      ".txt"
    )
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_OUT,
    paste0(
      "sessionInfo_",
      VERSION,
      ".txt"
    )
  )
)


cat("\n====================================================\n")
cat("BALANCED BROAD R8 UMAP COMPLETE\n")
cat("Version:", VERSION, "\n")

cat("\n=== ALL SAMPLE COUNTS ===\n")
print(
  sample_counts
)

cat("\n=== CONDITION VALIDATION ===\n")
print(
  condition_counts
)

cat("\n=== BALANCED SAMPLE COUNTS ===\n")
print(
  balanced_sample_counts
)

cat("\n=== BALANCED CONDITION COUNTS ===\n")
print(
  balanced_condition_counts
)

cat(
  "\nBalanced cells per biological sample:",
  n_per_sample,
  "\n"
)

cat(
  "\nOutput:",
  OUT,
  "\n"
)

cat("====================================================\n")
