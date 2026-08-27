#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Cd163-related candidate UMAP panel
# Visualization-improved version
#
# Version: v5.7.3
#
# CHANGE FROM v5.7.2
#   VISUALIZATION ONLY.
#   No reclustering, reintegration, re-UMAP, or biological reclassification.
#
# Improvements:
#   1) full Clean-B Mphi UMAP context
#   2) Kupffer_Macrophage auto-zoom
#   3) condition panels = 2 x 2
#   4) sample panels = 3 x 2
#   5) same expression scale across panels for each gene
#   6) positive cells are drawn last and larger
#   7) faint full-Kupffer background silhouette for split plots
#
# Target genes:
#   Cd163, Sema6d, Vsig4, Clec4f, Slc40a1, Colec12
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(plot_obj, filename, width, height) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot_obj)
  grDevices::dev.off()
}

resolve_first <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

canonical_condition <- function(x) {
  x <- as.character(x)

  out <- dplyr::case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )

  factor(
    out,
    levels = c("STD", "CDAHFD", "Sham", "Tx")
  )
}

resolve_gene_symbols <- function(requested, available) {
  available_lower <- tolower(available)

  resolved <- vapply(
    requested,
    FUN.VALUE = character(1),
    FUN = function(g) {
      idx <- match(tolower(g), available_lower)
      if (is.na(idx)) NA_character_ else available[[idx]]
    }
  )

  tibble(
    requested_gene = requested,
    resolved_gene = resolved,
    present = !is.na(resolved)
  )
}

calc_zoom_limits <- function(coords, pad_fraction = 0.06) {
  xr <- range(coords[, 1], finite = TRUE)
  yr <- range(coords[, 2], finite = TRUE)

  xpad <- max(diff(xr) * pad_fraction, 0.25)
  ypad <- max(diff(yr) * pad_fraction, 0.25)

  list(
    xlim = c(xr[1] - xpad, xr[2] + xpad),
    ylim = c(yr[1] - ypad, yr[2] + ypad)
  )
}

get_expression_limits <- function(expr) {
  pos <- expr[is.finite(expr) & expr > 0]

  if (!length(pos)) {
    return(c(0, 1))
  }

  lo <- as.numeric(quantile(pos, 0.02, na.rm = TRUE))
  hi <- as.numeric(quantile(pos, 0.98, na.rm = TRUE))

  if (!is.finite(lo)) lo <- min(pos, na.rm = TRUE)
  if (!is.finite(hi)) hi <- max(pos, na.rm = TRUE)

  if (hi <= lo) {
    hi <- max(pos, na.rm = TRUE)
  }

  if (hi <= lo) {
    hi <- lo + 1e-6
  }

  c(lo, hi)
}

make_plot_df <- function(obj, gene, reduction_use) {
  coords <- Embeddings(
    obj,
    reduction = reduction_use
  )

  expr <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "data"
  )[gene, , drop = TRUE]

  tibble(
    cell = colnames(obj),
    UMAP_1 = coords[colnames(obj), 1],
    UMAP_2 = coords[colnames(obj), 2],
    expression = as.numeric(expr[colnames(obj)])
  )
}

make_context_umap <- function(
  df,
  gene,
  limits_xy,
  expression_limits,
  title_text,
  pt_zero = 0.20,
  pt_pos = 0.55
) {

  zero_df <- df %>%
    filter(
      !is.finite(expression) |
        expression <= 0
    )

  pos_df <- df %>%
    filter(
      is.finite(expression),
      expression > 0
    ) %>%
    arrange(expression)

  midpoint_use <- mean(expression_limits)

  ggplot() +
    geom_point(
      data = zero_df,
      aes(
        x = UMAP_1,
        y = UMAP_2
      ),
      color = "grey88",
      size = pt_zero,
      alpha = 0.75
    ) +
    geom_point(
      data = pos_df,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        color = expression
      ),
      size = pt_pos,
      alpha = 0.95
    ) +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = midpoint_use,
      limits = expression_limits,
      oob = scales::squish,
      name = gene
    ) +
    coord_cartesian(
      xlim = limits_xy$xlim,
      ylim = limits_xy$ylim,
      expand = FALSE
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      legend.title = element_text(
        size = 9
      ),
      legend.text = element_text(
        size = 8
      ),
      plot.margin = margin(
        4,
        4,
        4,
        4
      )
    )
}

make_split_panel <- function(
  full_df,
  group_values,
  group_levels,
  gene,
  limits_xy,
  expression_limits,
  ncol,
  title_text,
  background_pt = 0.18,
  zero_pt = 0.32,
  positive_pt = 0.85
) {

  plot_list <- lapply(
    group_levels,
    function(group_now) {

      idx <- which(
        group_values ==
          group_now
      )

      group_df <- full_df[
        idx,
        ,
        drop = FALSE
      ]

      zero_df <- group_df %>%
        filter(
          !is.finite(expression) |
            expression <= 0
        )

      pos_df <- group_df %>%
        filter(
          is.finite(expression),
          expression > 0
        ) %>%
        arrange(expression)

      midpoint_use <- mean(expression_limits)

      ggplot() +
        geom_point(
          data = full_df,
          aes(
            x = UMAP_1,
            y = UMAP_2
          ),
          color = "grey94",
          size = background_pt,
          alpha = 0.45
        ) +
        geom_point(
          data = zero_df,
          aes(
            x = UMAP_1,
            y = UMAP_2
          ),
          color = "grey70",
          size = zero_pt,
          alpha = 0.72
        ) +
        geom_point(
          data = pos_df,
          aes(
            x = UMAP_1,
            y = UMAP_2,
            color = expression
          ),
          size = positive_pt,
          alpha = 1
        ) +
        scale_color_gradient2(
          low = "#0033FF",
          mid = "#FFFFFF",
          high = "#FF1A1A",
          midpoint = midpoint_use,
          limits = expression_limits,
          oob = scales::squish,
          name = gene
        ) +
        coord_cartesian(
          xlim = limits_xy$xlim,
          ylim = limits_xy$ylim,
          expand = FALSE
        ) +
        labs(
          title = group_now,
          x = NULL,
          y = NULL
        ) +
        theme_void(base_size = 10) +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 13,
            hjust = 0.5
          ),
          legend.title = element_text(
            size = 9
          ),
          legend.text = element_text(
            size = 8
          ),
          plot.margin = margin(
            3,
            3,
            3,
            3
          )
        )
    }
  )

  wrap_plots(
    plot_list,
    ncol = ncol,
    guides = "collect"
  ) +
    plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 15,
          hjust = 0.5
        )
      )
    ) &
    theme(
      legend.position = "right"
    )
}

# ==============================================================================
# 1. Paths and settings
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS_CANDIDATES <- c(
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
  ),
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
  )
)

MPHI_RDS <- MPHI_RDS_CANDIDATES[
  file.exists(MPHI_RDS_CANDIDATES)
][1]

OUTPUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "CD163_candidate_UMAP_panel_v5.7.3"
)

FIG_DIR <- file.path(
  OUTPUT_DIR,
  "Figures"
)

TAB_DIR <- file.path(
  OUTPUT_DIR,
  "Tables"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TAB_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

ASSAY_USE <- "RNA"

UMAP_CANDIDATES <- c(
  "umapRPCA",
  "umap",
  "UMAP"
)

CONDITION_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group"
)

SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "sample_display_FIXED2",
  "sample",
  "orig.ident"
)

LINEAGE_CANDIDATES <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_v440",
  "layer1_original",
  "celltype_for_R8plot",
  "vote_ueno_celltype_v34",
  "celltype_auto_annotation"
)

TARGET_LINEAGE <- "Kupffer_Macrophage"

REQUESTED_GENES <- c(
  "Cd163",
  "Sema6d",
  "Vsig4",
  "Clec4f",
  "Slc40a1",
  "Colec12"
)

CONDITION_LEVELS <- c(
  "STD",
  "CDAHFD",
  "Sham",
  "Tx"
)

SAMPLE_ORDER_PREFERRED <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "CDAHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

# ==============================================================================
# 2. Load
# ==============================================================================

if (
  length(MPHI_RDS) == 0L ||
  is.na(MPHI_RDS) ||
  !file.exists(MPHI_RDS)
) {
  stop(
    "Clean-B macrophage RDS not found."
  )
}

msg(
  "Loading Clean-B macrophage RDS..."
)

mphi <- readRDS(
  MPHI_RDS
)

if (
  !ASSAY_USE %in%
    Assays(mphi)
) {
  stop(
    "RNA assay not found."
  )
}

DefaultAssay(
  mphi
) <- ASSAY_USE

condition_col <- resolve_first(
  colnames(mphi@meta.data),
  CONDITION_CANDIDATES
)

sample_col <- resolve_first(
  colnames(mphi@meta.data),
  SAMPLE_CANDIDATES
)

lineage_col <- resolve_first(
  colnames(mphi@meta.data),
  LINEAGE_CANDIDATES
)

umap_use <- resolve_first(
  Reductions(mphi),
  UMAP_CANDIDATES
)

if (is.na(condition_col)) {
  stop(
    "Could not resolve condition column."
  )
}

if (is.na(sample_col)) {
  stop(
    "Could not resolve sample column."
  )
}

if (is.na(lineage_col)) {
  stop(
    "Could not resolve lineage column."
  )
}

if (is.na(umap_use)) {
  stop(
    "Could not resolve UMAP reduction."
  )
}

msg(
  "Condition source: ",
  condition_col
)

msg(
  "Sample source: ",
  sample_col
)

msg(
  "Lineage source: ",
  lineage_col
)

msg(
  "UMAP reduction: ",
  umap_use
)

mphi$condition_v573 <- canonical_condition(
  mphi@meta.data[[condition_col]]
)

mphi$sample_v573 <- as.character(
  mphi@meta.data[[sample_col]]
)

# ==============================================================================
# 3. Gene audit
# ==============================================================================

gene_audit <- resolve_gene_symbols(
  REQUESTED_GENES,
  rownames(mphi)
)

write.csv(
  gene_audit,
  file.path(
    TAB_DIR,
    "00_candidate_gene_audit_v5.7.3.csv"
  ),
  row.names = FALSE
)

print(
  gene_audit
)

genes_use <- gene_audit %>%
  filter(present) %>%
  pull(resolved_gene) %>%
  unique()

if (!length(genes_use)) {
  stop(
    "None of the requested genes were detected."
  )
}

msg(
  "Genes: ",
  paste(
    genes_use,
    collapse = ", "
  )
)

# ==============================================================================
# 4. Full-Mphi and Kupffer objects
# ==============================================================================

lineage_values <- as.character(
  mphi@meta.data[[lineage_col]]
)

kupffer_cells <- colnames(mphi)[
  lineage_values ==
    TARGET_LINEAGE
]

if (
  length(kupffer_cells) < 50
) {
  stop(
    "Too few Kupffer_Macrophage cells."
  )
}

kupffer <- subset(
  mphi,
  cells = kupffer_cells
)

DefaultAssay(
  kupffer
) <- ASSAY_USE

msg(
  "All Clean-B Mphi cells: ",
  ncol(mphi)
)

msg(
  "Kupffer_Macrophage cells: ",
  ncol(kupffer)
)

full_coords <- Embeddings(
  mphi,
  reduction = umap_use
)

kupffer_coords <- Embeddings(
  kupffer,
  reduction = umap_use
)

full_limits <- calc_zoom_limits(
  full_coords,
  pad_fraction = 0.03
)

kupffer_limits <- calc_zoom_limits(
  kupffer_coords,
  pad_fraction = 0.05
)

write.csv(
  tibble(
    region = c(
      "Full_Mphi",
      "Kupffer_zoom"
    ),
    xmin = c(
      full_limits$xlim[1],
      kupffer_limits$xlim[1]
    ),
    xmax = c(
      full_limits$xlim[2],
      kupffer_limits$xlim[2]
    ),
    ymin = c(
      full_limits$ylim[1],
      kupffer_limits$ylim[1]
    ),
    ymax = c(
      full_limits$ylim[2],
      kupffer_limits$ylim[2]
    )
  ),
  file.path(
    TAB_DIR,
    "01_UMAP_plot_limits_v5.7.3.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 5. Integrated full-Mphi overview
# ==============================================================================

msg(
  "Generating full-Mphi overview..."
)

full_plots <- list()

for (gene in genes_use) {

  df_gene <- make_plot_df(
    mphi,
    gene,
    umap_use
  )

  expression_limits <- get_expression_limits(
    df_gene$expression
  )

  full_plots[[gene]] <- make_context_umap(
    df = df_gene,
    gene = gene,
    limits_xy = full_limits,
    expression_limits = expression_limits,
    title_text = gene,
    pt_zero = 0.16,
    pt_pos = 0.48
  )
}

p_full <- wrap_plots(
  full_plots,
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(
    title =
      "Cd163-related candidates | full Clean-B macrophage UMAP"
  ) &
  theme(
    legend.position = "right"
  )

save_pdf(
  p_full,
  file.path(
    FIG_DIR,
    "01_FULL_MPHI_candidate_UMAP_overview_v5.7.3.pdf"
  ),
  width = 14,
  height = 9
)

# ==============================================================================
# 6. Kupffer zoom integrated overview
# ==============================================================================

msg(
  "Generating Kupffer zoom overview..."
)

kupffer_plots <- list()

for (gene in genes_use) {

  df_gene <- make_plot_df(
    kupffer,
    gene,
    umap_use
  )

  expression_limits <- get_expression_limits(
    df_gene$expression
  )

  kupffer_plots[[gene]] <- make_context_umap(
    df = df_gene,
    gene = gene,
    limits_xy = kupffer_limits,
    expression_limits = expression_limits,
    title_text = gene,
    pt_zero = 0.24,
    pt_pos = 0.70
  )
}

p_kupffer <- wrap_plots(
  kupffer_plots,
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(
    title =
      "Cd163-related candidates | Kupffer_Macrophage zoom"
  ) &
  theme(
    legend.position = "right"
  )

save_pdf(
  p_kupffer,
  file.path(
    FIG_DIR,
    "02_KUPFFER_ZOOM_candidate_UMAP_overview_v5.7.3.pdf"
  ),
  width = 14,
  height = 9
)

# ==============================================================================
# 7. Condition split: 2 x 2
# ==============================================================================

msg(
  "Generating condition split UMAPs (2 x 2)..."
)

condition_values <- as.character(
  kupffer$condition_v573
)

condition_levels_use <- CONDITION_LEVELS[
  CONDITION_LEVELS %in%
    unique(condition_values)
]

for (gene in genes_use) {

  df_gene <- make_plot_df(
    kupffer,
    gene,
    umap_use
  )

  expression_limits <- get_expression_limits(
    df_gene$expression
  )

  p_condition <- make_split_panel(
    full_df = df_gene,
    group_values = condition_values,
    group_levels = condition_levels_use,
    gene = gene,
    limits_xy = kupffer_limits,
    expression_limits = expression_limits,
    ncol = 2,
    title_text = paste0(
      gene,
      " | condition split | Kupffer zoom"
    ),
    background_pt = 0.18,
    zero_pt = 0.36,
    positive_pt = 0.95
  )

  save_pdf(
    p_condition,
    file.path(
      FIG_DIR,
      paste0(
        "03_CONDITION_2x2_KUPFFER_ZOOM_",
        gene,
        "_v5.7.3.pdf"
      )
    ),
    width = 11,
    height = 10
  )
}

# ==============================================================================
# 8. Sample split: 3 x 2
# ==============================================================================

msg(
  "Generating sample split UMAPs (3 x 2)..."
)

sample_values <- as.character(
  kupffer$sample_v573
)

sample_present <- unique(
  sample_values
)

sample_levels_use <- c(
  SAMPLE_ORDER_PREFERRED[
    SAMPLE_ORDER_PREFERRED %in%
      sample_present
  ],
  sort(
    setdiff(
      sample_present,
      SAMPLE_ORDER_PREFERRED
    )
  )
)

for (gene in genes_use) {

  df_gene <- make_plot_df(
    kupffer,
    gene,
    umap_use
  )

  expression_limits <- get_expression_limits(
    df_gene$expression
  )

  p_sample <- make_split_panel(
    full_df = df_gene,
    group_values = sample_values,
    group_levels = sample_levels_use,
    gene = gene,
    limits_xy = kupffer_limits,
    expression_limits = expression_limits,
    ncol = 3,
    title_text = paste0(
      gene,
      " | biological samples | Kupffer zoom"
    ),
    background_pt = 0.16,
    zero_pt = 0.34,
    positive_pt = 0.90
  )

  save_pdf(
    p_sample,
    file.path(
      FIG_DIR,
      paste0(
        "04_SAMPLE_3x2_KUPFFER_ZOOM_",
        gene,
        "_v5.7.3.pdf"
      )
    ),
    width = 13,
    height = 9
  )
}

# ==============================================================================
# 9. Sham vs Tx only: large 1 x 2
# ==============================================================================

msg(
  "Generating large Sham-vs-Tx comparison UMAPs..."
)

sham_tx_levels <- c(
  "Sham",
  "Tx"
)

for (gene in genes_use) {

  df_gene <- make_plot_df(
    kupffer,
    gene,
    umap_use
  )

  expression_limits <- get_expression_limits(
    df_gene$expression
  )

  p_sham_tx <- make_split_panel(
    full_df = df_gene,
    group_values = condition_values,
    group_levels = sham_tx_levels,
    gene = gene,
    limits_xy = kupffer_limits,
    expression_limits = expression_limits,
    ncol = 2,
    title_text = paste0(
      gene,
      " | Sham vs Tx | Kupffer zoom"
    ),
    background_pt = 0.20,
    zero_pt = 0.42,
    positive_pt = 1.15
  )

  save_pdf(
    p_sham_tx,
    file.path(
      FIG_DIR,
      paste0(
        "05_SHAM_vs_TX_LARGE_KUPFFER_ZOOM_",
        gene,
        "_v5.7.3.pdf"
      )
    ),
    width = 11,
    height = 5.8
  )
}

# ==============================================================================
# 10. Output index
# ==============================================================================

index <- tibble(
  category = c(
    "full macrophage overview",
    "Kupffer zoom overview",
    rep(
      "condition 2x2",
      length(genes_use)
    ),
    rep(
      "sample 3x2",
      length(genes_use)
    ),
    rep(
      "Sham vs Tx large",
      length(genes_use)
    )
  ),
  file = c(
    "01_FULL_MPHI_candidate_UMAP_overview_v5.7.3.pdf",
    "02_KUPFFER_ZOOM_candidate_UMAP_overview_v5.7.3.pdf",
    paste0(
      "03_CONDITION_2x2_KUPFFER_ZOOM_",
      genes_use,
      "_v5.7.3.pdf"
    ),
    paste0(
      "04_SAMPLE_3x2_KUPFFER_ZOOM_",
      genes_use,
      "_v5.7.3.pdf"
    ),
    paste0(
      "05_SHAM_vs_TX_LARGE_KUPFFER_ZOOM_",
      genes_use,
      "_v5.7.3.pdf"
    )
  )
)

write.csv(
  index,
  file.path(
    OUTPUT_DIR,
    "FIGURE_INDEX_v5.7.3.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo_v5.7.3.txt"
  )
)

msg(
  "DONE."
)

msg(
  "Output: ",
  OUTPUT_DIR
)

msg(
  "Figures: ",
  FIG_DIR
)
