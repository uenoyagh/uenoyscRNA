#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(7023)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Sham vs Tx whole-cell refined UMAP using sample names encoded in cell barcodes
#
# Version: v7.0.2.4
#
# FIX FROM v7.0.2.3
#   Sample recovery from cell-name prefixes was successful and validated:
#     STD    = 12,211
#     CDAHFD = 11,543
#     Sham   = 38,382
#     Tx     = 42,452
#
#   v7.0.2.3 then stopped only at plotting because `guides()` had been
#   incorrectly added to a ggplot2 theme object inside common_plot_theme().
#   v7.0.2.4 removes that invalid theme+guides combination.
#
# FIX
#   v7.0.2 / v7.0.2.1 / v7.0.2.2 tried to recover Sham/Tx sample identity
#   from metadata columns, but the available whole-cell objects have usable
#   sample metadata only for STD/CDHFD in those columns.
#
#   The whole-cell barcode names themselves retain the biological-sample prefix:
#
#     STD_rep1_...
#     CDHFD_rep1_...
#     Sham1_...
#     Sham20_...
#     Tx17_...
#     Tx5_...
#
#   Therefore v7.0.2.3 derives biological sample identity DIRECTLY from the
#   exact cell-name prefix. No external sample-metadata RDS is required.
#
# PURPOSE
#   1) Sham vs Tx whole-cell UMAP, all cells
#   2) Sham vs Tx balanced-display UMAP
#   3) Sham1 / Sham20 / Tx17 / Tx5 replicate UMAP
#   4) Sham vs Tx Mphi + Monocyte highlight UMAP
#
# IMPORTANT
#   - Same frozen whole-cell umapRPCA coordinates
#   - Same v7.0.1 refined annotations
#   - Same v7.0.1 R8 palette
#   - No reintegration
#   - No reclustering
#   - No UMAP recomputation
#   - Balanced plot is visualization-only
# ==============================================================================


# ==============================================================================
# 1. Paths
# ==============================================================================

VERSION <- "v7.0.2.4"

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

R8_PALETTE_CSV <- file.path(
  V701_DIR,
  "Tables",
  "08_R8_refined_palette_manifest_v7.0.1.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_RefinedAnnotation_UMAP_v7.0.2.4_Sham_vs_Tx"
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

PT_ALL <- 0.42
PT_SAMPLE <- 0.36
PT_MYEL <- 0.46


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
    grepl(
      "^STD_rep1_",
      cells
    )
  ] <- "STD_rep1"

  out[
    grepl(
      "^CDHFD_rep1_",
      cells
    )
  ] <- "CDHFD_rep1"

  out[
    grepl(
      "^Sham20_",
      cells
    )
  ] <- "Sham20"

  out[
    grepl(
      "^Sham1_",
      cells
    )
  ] <- "Sham1"

  out[
    grepl(
      "^Tx17_",
      cells
    )
  ] <- "Tx17"

  out[
    grepl(
      "^Tx5_",
      cells
    )
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


common_plot_theme <- function(
  legend_ncol = 3,
  base_size = 9
) {

  # IMPORTANT:
  # Return a theme object ONLY.
  #
  # In current ggplot2, `guides()` is not a theme and cannot be added
  # inside `theme_classic() + theme(...)`.  v7.0.2.3 stopped here.
  # Legend guide formatting is intentionally left to DimPlot/patchwork.
  theme_classic(
    base_size =
      base_size
  ) +
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        ),
      legend.title =
        element_blank(),
      legend.text =
        element_text(
          size = 6.6
        ),
      legend.key.height =
        grid::unit(
          0.34,
          "cm"
        ),
      legend.key.width =
        grid::unit(
          0.34,
          "cm"
        )
    )
}


# ==============================================================================
# 3. Input checks
# ==============================================================================

for (
  f in c(
    WHOLE_RDS,
    REFINED_METADATA_CSV,
    R8_PALETTE_CSV
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
# 4. Load frozen whole-cell coordinate source
# ==============================================================================

msg(
  "Loading whole-cell frozen parent..."
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
# 5. Recover biological sample directly from cell-name prefix
# ==============================================================================

msg(
  "Recovering biological sample from cell-name prefix..."
)

whole$sample_v7023 <-
  sample_from_cellname(
    whole_cells
  )

whole$condition_v7023 <-
  sample_to_condition(
    whole$sample_v7023
  )

n_unresolved <- sum(
  is.na(
    whole$sample_v7023
  )
)

if (n_unresolved > 0) {

  unresolved_examples <- head(
    whole_cells[
      is.na(
        whole$sample_v7023
      )
    ],
    30
  )

  stop(
    "Unresolved biological-sample prefixes: ",
    n_unresolved,
    " cells.\nExamples:\n",
    paste(
      unresolved_examples,
      collapse = "\n"
    )
  )
}


sample_counts <- as.data.frame(
  table(
    sample =
      factor(
        whole$sample_v7023,
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

write.csv(
  sample_counts,
  file.path(
    TAB_OUT,
    paste0(
      "01_cellname_prefix_sample_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


condition_counts <- as.data.frame(
  table(
    condition =
      factor(
        whole$condition_v7023,
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

write.csv(
  condition_counts,
  file.path(
    TAB_OUT,
    paste0(
      "02_cellname_prefix_condition_validation_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


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
  sum(
    sample_counts$n_cells
  ) !=
    ncol(
      whole
    )
) {
  stop(
    "Sample counts do not sum to whole-cell total."
  )
}


if (
  !all(
    condition_counts$matches_expected
  )
) {
  stop(
    "Condition totals recovered from barcode prefixes do not match the ",
    "validated whole-cell totals. Review before plotting."
  )
}


target_counts <- sample_counts[
  sample_counts$sample %in%
    TARGET_SAMPLES,
  ,
  drop = FALSE
]

if (
  nrow(target_counts) != 4 ||
  any(
    target_counts$n_cells <= 0
  )
) {
  stop(
    "Sham1/Sham20/Tx17/Tx5 were not all recovered."
  )
}


# ==============================================================================
# 6. Load v7.0.1 refined annotation by exact barcode
# ==============================================================================

meta701 <- read.csv(
  REFINED_METADATA_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta_cols <- c(
  "cell",
  "broad_annotation",
  "refined_annotation"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(
    meta701
  )
)

if (length(missing_meta_cols)) {
  stop(
    "v7.0.1 metadata missing columns: ",
    paste(
      missing_meta_cols,
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
    "Duplicated cell barcodes in v7.0.1 refined metadata."
  )
}

idx_meta <- match(
  whole_cells,
  meta701$cell
)

if (anyNA(idx_meta)) {

  stop(
    "v7.0.1 refined metadata does not cover all whole-cell barcodes. Missing: ",
    sum(
      is.na(
        idx_meta
      )
    )
  )
}

whole$broad_annotation_v7023 <-
  meta701$broad_annotation[
    idx_meta
  ]

whole$refined_annotation_v7023 <-
  meta701$refined_annotation[
    idx_meta
  ]


# ==============================================================================
# 7. Load exact v7.0.1 R8 palette
# ==============================================================================

palette_df <- read.csv(
  R8_PALETTE_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  !all(
    c(
      "annotation",
      "hex"
    ) %in%
      colnames(
        palette_df
      )
  )
) {
  stop(
    "R8 palette manifest must contain annotation and hex."
  )
}

R8_REFINED <- setNames(
  palette_df$hex,
  palette_df$annotation
)

refined_present <- unique(
  whole$refined_annotation_v7023
)

missing_colors <- setdiff(
  refined_present,
  names(
    R8_REFINED
  )
)

if (length(missing_colors)) {

  fallback <- grDevices::hcl.colors(
    length(
      missing_colors
    ),
    palette =
      "Dynamic"
  )

  names(
    fallback
  ) <- missing_colors

  R8_REFINED <- c(
    R8_REFINED,
    fallback
  )

  warning(
    "Fallback colors assigned to: ",
    paste(
      missing_colors,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 8. Build Sham/Tx object
# ==============================================================================

target_cells <- whole_cells[
  whole$sample_v7023 %in%
    TARGET_SAMPLES
]

shamtx <- subset(
  whole,
  cells =
    target_cells
)

shamtx$sample_v7023 <-
  factor(
    shamtx$sample_v7023,
    levels =
      TARGET_SAMPLES
  )

shamtx$condition_v7023 <-
  factor(
    shamtx$condition_v7023,
    levels =
      TARGET_CONDITIONS
  )

refined_order <- names(
  R8_REFINED
)

refined_order <- refined_order[
  refined_order %in%
    unique(
      shamtx$refined_annotation_v7023
    )
]

shamtx$refined_annotation_v7023 <-
  factor(
    shamtx$refined_annotation_v7023,
    levels =
      refined_order
  )


# ==============================================================================
# 9. Figure 01: Sham vs Tx full display
# ==============================================================================

msg(
  "Drawing Sham vs Tx full refined UMAP..."
)

p_full <- DimPlot(
  shamtx,
  reduction =
    WHOLE_UMAP,
  group.by =
    "refined_annotation_v7023",
  split.by =
    "condition_v7023",
  cols =
    R8_REFINED,
  pt.size =
    PT_ALL,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70231
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | Sham vs Tx | current refined annotations | R8"
  )

p_full <-
  p_full &
  common_plot_theme(
    legend_ncol = 3,
    base_size = 8.5
  )

save_pdf(
  p_full,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_refined_Sham_vs_Tx_FULL_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 20,
  height = 11
)

save_png(
  p_full,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_refined_Sham_vs_Tx_FULL_R8_",
      VERSION,
      ".png"
    )
  ),
  width = 20,
  height = 11
)


# ==============================================================================
# 10. Figure 02: balanced display
#
# Equal number of cells from EACH biological sample.
# Therefore total displayed Sham and Tx cell numbers are also equal.
# Visualization only.
# ==============================================================================

n_per_sample <- min(
  target_counts$n_cells
)

msg(
  "Balanced display cells per biological sample: ",
  n_per_sample
)

set.seed(
  70232
)

balanced_cells <- unlist(
  lapply(
    TARGET_SAMPLES,
    function(s) {

      cells_now <- colnames(
        shamtx
      )[
        shamtx$sample_v7023 ==
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
  shamtx,
  cells =
    balanced_cells
)

balanced$sample_v7023 <-
  factor(
    balanced$sample_v7023,
    levels =
      TARGET_SAMPLES
  )

balanced$condition_v7023 <-
  factor(
    balanced$condition_v7023,
    levels =
      TARGET_CONDITIONS
  )

balanced$refined_annotation_v7023 <-
  factor(
    balanced$refined_annotation_v7023,
    levels =
      refined_order
  )

balanced_counts <- as.data.frame(
  table(
    sample =
      balanced$sample_v7023,
    condition =
      balanced$condition_v7023
  ),
  stringsAsFactors = FALSE
)

balanced_counts <- balanced_counts[
  balanced_counts$Freq > 0,
  ,
  drop = FALSE
]

names(
  balanced_counts
) <- c(
  "sample",
  "condition",
  "n_cells"
)

write.csv(
  balanced_counts,
  file.path(
    TAB_OUT,
    paste0(
      "03_balanced_display_cell_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


p_balanced <- DimPlot(
  balanced,
  reduction =
    WHOLE_UMAP,
  group.by =
    "refined_annotation_v7023",
  split.by =
    "condition_v7023",
  cols =
    R8_REFINED,
  pt.size =
    PT_ALL,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70233
) +
  plot_annotation(
    title = paste0(
      "Mouse liver whole-cell UMAP | Sham vs Tx | balanced display | ",
      n_per_sample,
      " cells / biological sample | R8"
    )
  )

p_balanced <-
  p_balanced &
  common_plot_theme(
    legend_ncol = 3,
    base_size = 8.5
  )

save_pdf(
  p_balanced,
  file.path(
    FIG_OUT,
    paste0(
      "02_wholecell_refined_Sham_vs_Tx_BALANCED_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 20,
  height = 11
)

save_png(
  p_balanced,
  file.path(
    FIG_OUT,
    paste0(
      "02_wholecell_refined_Sham_vs_Tx_BALANCED_R8_",
      VERSION,
      ".png"
    )
  ),
  width = 20,
  height = 11
)


# ==============================================================================
# 11. Figure 03: biological-replicate UMAPs
# ==============================================================================

msg(
  "Drawing Sham1 / Sham20 / Tx17 / Tx5 UMAPs..."
)

p_samples <- DimPlot(
  shamtx,
  reduction =
    WHOLE_UMAP,
  group.by =
    "refined_annotation_v7023",
  split.by =
    "sample_v7023",
  cols =
    R8_REFINED,
  pt.size =
    PT_SAMPLE,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70234
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | Sham1 / Sham20 / Tx17 / Tx5 | R8"
  )

p_samples <-
  p_samples &
  common_plot_theme(
    legend_ncol = 3,
    base_size = 8
  )

save_pdf(
  p_samples,
  file.path(
    FIG_OUT,
    paste0(
      "03_wholecell_refined_Sham1_Sham20_Tx17_Tx5_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 20,
  height = 15
)


# ==============================================================================
# 12. Figure 04: Mphi + Monocyte focused UMAP
# ==============================================================================

msg(
  "Drawing Sham vs Tx Mphi + Monocyte highlight..."
)

annotation_char <- as.character(
  shamtx$refined_annotation_v7023
)

is_myeloid_refined <- grepl(
  "^Mphi \\||^Mono \\|",
  annotation_char
)

myeloid_display <- ifelse(
  is_myeloid_refined,
  annotation_char,
  "Other cells"
)

myeloid_levels <- c(
  "Other cells",
  refined_order[
    grepl(
      "^Mphi \\||^Mono \\|",
      refined_order
    )
  ]
)

shamtx$myeloid_display_v7023 <-
  factor(
    myeloid_display,
    levels =
      myeloid_levels
  )

myeloid_palette <- c(
  "Other cells" =
    "#E6E6E6",
  R8_REFINED[
    names(
      R8_REFINED
    ) %in%
      myeloid_levels
  ]
)

p_myeloid <- DimPlot(
  shamtx,
  reduction =
    WHOLE_UMAP,
  group.by =
    "myeloid_display_v7023",
  split.by =
    "condition_v7023",
  cols =
    myeloid_palette,
  pt.size =
    PT_MYEL,
  raster =
    FALSE,
  ncol = 2,
  shuffle =
    TRUE,
  seed =
    70235
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | Sham vs Tx | Mphi + Monocyte focus | R8"
  )

p_myeloid <-
  p_myeloid &
  common_plot_theme(
    legend_ncol = 2,
    base_size = 9
  )

save_pdf(
  p_myeloid,
  file.path(
    FIG_OUT,
    paste0(
      "04_wholecell_Sham_vs_Tx_Mphi_Monocyte_highlight_R8_",
      VERSION,
      ".pdf"
    )
  ),
  width = 17,
  height = 10
)

save_png(
  p_myeloid,
  file.path(
    FIG_OUT,
    paste0(
      "04_wholecell_Sham_vs_Tx_Mphi_Monocyte_highlight_R8_",
      VERSION,
      ".png"
    )
  ),
  width = 17,
  height = 10
)


# ==============================================================================
# 13. Refined annotation count table by biological sample
# ==============================================================================

annotation_sample_counts <- as.data.frame(
  table(
    sample =
      shamtx$sample_v7023,
    condition =
      shamtx$condition_v7023,
    refined_annotation =
      shamtx$refined_annotation_v7023
  ),
  stringsAsFactors = FALSE
)

annotation_sample_counts <-
  annotation_sample_counts[
    annotation_sample_counts$Freq > 0,
    ,
    drop = FALSE
  ]

names(
  annotation_sample_counts
) <- c(
  "sample",
  "condition",
  "refined_annotation",
  "n_cells"
)

write.csv(
  annotation_sample_counts,
  file.path(
    TAB_OUT,
    paste0(
      "04_ShTx_refined_annotation_by_sample_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. README / session info
# ==============================================================================

readme <- c(
  paste0(
    "Mouse MASH Sham vs Tx refined whole-cell UMAP ",
    VERSION
  ),
  "============================================================",
  "",
  paste0(
    "Whole-cell / UMAP source: ",
    WHOLE_RDS
  ),
  paste0(
    "Refined annotation source: ",
    REFINED_METADATA_CSV
  ),
  paste0(
    "R8 palette source: ",
    R8_PALETTE_CSV
  ),
  "",
  "Biological sample is recovered directly from exact cell-name prefix:",
  "- STD_rep1_",
  "- CDHFD_rep1_",
  "- Sham1_",
  "- Sham20_",
  "- Tx17_",
  "- Tx5_",
  "",
  "No external sample metadata source is used.",
  "No reintegration, reclustering, or UMAP recomputation.",
  "",
  paste0(
    "Balanced display uses ",
    n_per_sample,
    " cells from each biological sample."
  ),
  "Balanced display is visualization-only."
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
cat("SHAM vs TX WHOLE-CELL UMAP COMPLETE\n")
cat("Version:", VERSION, "\n")

cat("\n=== CELL-NAME PREFIX SAMPLE COUNTS ===\n")
print(
  sample_counts
)

cat("\n=== CONDITION COUNT VALIDATION ===\n")
print(
  condition_counts
)

cat("\n=== TARGET SHAM/TX COUNTS ===\n")
print(
  target_counts
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
