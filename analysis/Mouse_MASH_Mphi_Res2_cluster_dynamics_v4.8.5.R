#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# Res2 manual annotation v4.8.4
# Cluster transition + group-wise cluster summary
# v4.8.5
# ==============================================================================
#
# INPUT:
#   Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds
#
# PURPOSE:
#   1) Show cluster trajectories across:
#        STD -> CDAHFD -> Sham -> Tx
#      within each manually assigned MΦ class.
#
#   2) Show group-wise cluster summaries for each condition:
#        - cluster fraction of all MΦ
#        - cluster cell count
#
#   3) Keep v4.8.4 manual classification fixed.
#
# IMPORTANT:
#   RPCA / UMAP / Res2 clustering are NOT recalculated.
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4850)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

INPUT_RDS <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_manual_annotation_v4.8.4",
  "RDS",
  "Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_cluster_transition_v4.8.5"
)

FIG_OUT_DIR <- file.path(OUTPUT_DIR, "Figures")
CSV_OUT_DIR <- file.path(OUTPUT_DIR, "Tables")

dir.create(FIG_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSV_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Metadata
# ------------------------------------------------------------------------------

RES2_CLUSTER_COL <- "cluster_res2"
CLASS_COL <- "macrophage_class_Res2_v484"
SAMPLE_COL <- "sample_4group"
CONDITION_COL <- "condition_4group"

CONDITION_ORDER <- c(
  "STD",
  "CDAHFD",
  "Sham",
  "Tx"
)

CLASS_ORDER <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "Fibrogenic-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Other"
)

CLASS_LABELS <- c(
  "Inflammatory-Mphi" = "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi" = "Anti-inflammatory-Mphi",
  "Fibrogenic-Mphi" = "Fibrogenic-Mphi",
  "Repair/Resolution-Mphi" = "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-Mphi",
  "Other" = "Other"
)

CONDITION_COLORS <- c(
  "STD"    = "#2F65FF",
  "CDAHFD" = "#F04444",
  "Sham"   = "#777777",
  "Tx"     = "#F28C18"
)

CLASS_COLORS <- c(
  "Inflammatory-Mphi"           = "#E31A1C",
  "Anti-inflammatory-Mphi"      = "#1478FF",
  "Fibrogenic-Mphi"             = "#B218B2",
  "Repair/Resolution-Mphi"      = "#00A65A",
  "Lipid-associated/TREM2-Mphi" = "#F28C18",
  "Other"                       = "#B5B5B5"
)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------

message_time <- function(...) {
  message(
    format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
    ...
  )
}

numeric_cluster_levels <- function(x) {
  z <- unique(as.character(x))
  suppressWarnings(n <- as.numeric(z))

  if (all(!is.na(n))) {
    z[order(n)]
  } else {
    sort(z)
  }
}

save_pdf <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(FIG_OUT_DIR, filename),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

save_png_600 <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(FIG_OUT_DIR, filename),
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    limitsize = FALSE
  )
}

save_jpeg_600 <- function(
  filename,
  plot,
  width,
  height,
  quality = 95
) {

  outfile <- file.path(FIG_OUT_DIR, filename)

  if (requireNamespace("ragg", quietly = TRUE)) {

    ragg::agg_jpeg(
      filename = outfile,
      width = width,
      height = height,
      units = "in",
      res = 600,
      quality = quality,
      background = "white"
    )

    print(plot)
    grDevices::dev.off()

  } else {

    grDevices::jpeg(
      filename = outfile,
      width = width,
      height = height,
      units = "in",
      res = 600,
      quality = quality,
      bg = "white"
    )

    print(plot)
    grDevices::dev.off()
  }
}

# R8tone-like 28 cluster palette
R8TONE_BASE <- c(
  "#F04444",
  "#F28C18",
  "#D8B400",
  "#55A600",
  "#00B85A",
  "#00BFC4",
  "#2F65FF",
  "#B84BE8"
)

make_r8tone_palette <- function(levels_vec) {

  levels_vec <- as.character(levels_vec)
  n <- length(levels_vec)

  tone2 <- c(
    "#D62828",
    "#D96D00",
    "#B49700",
    "#3F8F00",
    "#00994B",
    "#009DA6",
    "#174ED1",
    "#9634C9"
  )

  tone3 <- c(
    "#FF6B6B",
    "#FFA340",
    "#E8C83C",
    "#78C928",
    "#29C978",
    "#28D3D8",
    "#5B82FF",
    "#CD6AF0"
  )

  tone4 <- c(
    "#A91E1E",
    "#B75400",
    "#8C7600",
    "#307000",
    "#00793B",
    "#007A82",
    "#103AA8",
    "#74269F"
  )

  pool <- c(
    R8TONE_BASE,
    tone2,
    tone3,
    tone4
  )

  cols <- pool[seq_len(n)]
  names(cols) <- levels_vec

  cols
}

# ------------------------------------------------------------------------------
# 4. Load object
# ------------------------------------------------------------------------------

message_time("Loading: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
  stop("Input object is not a Seurat object.")
}

required_cols <- c(
  RES2_CLUSTER_COL,
  CLASS_COL,
  SAMPLE_COL,
  CONDITION_COL
)

missing_cols <- setdiff(
  required_cols,
  colnames(mphi@meta.data)
)

if (length(missing_cols) > 0L) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

cluster_levels <- numeric_cluster_levels(
  mphi@meta.data[[RES2_CLUSTER_COL]]
)

mphi@meta.data[[RES2_CLUSTER_COL]] <- factor(
  as.character(
    mphi@meta.data[[RES2_CLUSTER_COL]]
  ),
  levels = cluster_levels
)

mphi@meta.data[[CLASS_COL]] <- factor(
  as.character(
    mphi@meta.data[[CLASS_COL]]
  ),
  levels = CLASS_ORDER
)

mphi@meta.data[[CONDITION_COL]] <- factor(
  as.character(
    mphi@meta.data[[CONDITION_COL]]
  ),
  levels = CONDITION_ORDER
)

message_time(
  "Cells: ",
  ncol(mphi)
)

# ------------------------------------------------------------------------------
# 5. Cell-level metadata
# ------------------------------------------------------------------------------

cell_meta <- mphi@meta.data %>%
  transmute(
    sample = as.character(.data[[SAMPLE_COL]]),

    condition = factor(
      as.character(.data[[CONDITION_COL]]),
      levels = CONDITION_ORDER
    ),

    cluster = factor(
      as.character(.data[[RES2_CLUSTER_COL]]),
      levels = cluster_levels
    ),

    macrophage_class = factor(
      as.character(.data[[CLASS_COL]]),
      levels = CLASS_ORDER
    )
  )

sample_condition <- cell_meta %>%
  distinct(
    sample,
    condition
  )

cluster_class_map <- cell_meta %>%
  distinct(
    cluster,
    macrophage_class
  )

write.csv(
  cluster_class_map,
  file.path(
    CSV_OUT_DIR,
    "00_cluster_class_map_v4.8.5.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Per-sample cluster counts
# ------------------------------------------------------------------------------

raw_counts <- cell_meta %>%
  count(
    sample,
    condition,
    cluster,
    macrophage_class,
    name = "n_cells"
  )

sample_cluster <- tidyr::expand_grid(
  sample = unique(cell_meta$sample),
  cluster = factor(
    cluster_levels,
    levels = cluster_levels
  )
) %>%
  left_join(
    sample_condition,
    by = "sample"
  ) %>%
  left_join(
    cluster_class_map,
    by = "cluster"
  ) %>%
  left_join(
    raw_counts,
    by = c(
      "sample",
      "condition",
      "cluster",
      "macrophage_class"
    )
  ) %>%
  mutate(
    n_cells = tidyr::replace_na(
      n_cells,
      0L
    )
  ) %>%
  group_by(sample) %>%
  mutate(
    total_mphi = sum(n_cells),

    fraction_all_mphi = ifelse(
      total_mphi > 0,
      n_cells / total_mphi,
      0
    ),

    percent_all_mphi =
      100 * fraction_all_mphi
  ) %>%
  ungroup() %>%
  group_by(
    sample,
    macrophage_class
  ) %>%
  mutate(
    total_class = sum(n_cells),

    fraction_within_class = ifelse(
      total_class > 0,
      n_cells / total_class,
      0
    ),

    percent_within_class =
      100 * fraction_within_class
  ) %>%
  ungroup()

write.csv(
  sample_cluster,
  file.path(
    CSV_OUT_DIR,
    "01_cluster_abundance_by_sample_v4.8.5.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Condition-level means
# ------------------------------------------------------------------------------

condition_cluster <- sample_cluster %>%
  group_by(
    condition,
    macrophage_class,
    cluster
  ) %>%
  summarise(
    n_samples = n_distinct(sample),

    mean_cells = mean(n_cells),

    mean_fraction_all_mphi =
      mean(fraction_all_mphi),

    mean_percent_all_mphi =
      mean(percent_all_mphi),

    mean_fraction_within_class =
      mean(fraction_within_class),

    mean_percent_within_class =
      mean(percent_within_class),

    .groups = "drop"
  )

write.csv(
  condition_cluster,
  file.path(
    CSV_OUT_DIR,
    "02_cluster_condition_summary_v4.8.5.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Cluster palette
# ------------------------------------------------------------------------------

cluster_palette <- make_r8tone_palette(
  cluster_levels
)

# ------------------------------------------------------------------------------
# 9. Class-faceted trajectory plot
#    cluster fraction of ALL MΦ
# ------------------------------------------------------------------------------

line_label_all <- condition_cluster %>%
  filter(
    condition == "Tx"
  ) %>%
  mutate(
    cluster_label = paste0(
      "#",
      as.character(cluster)
    )
  )

p_slope_all <- ggplot(
  condition_cluster,
  aes(
    x = condition,
    y = mean_percent_all_mphi,
    group = cluster,
    color = cluster
  )
) +
  geom_line(
    linewidth = 0.95,
    alpha = 0.90
  ) +
  geom_point(
    size = 2.8
  ) +
  ggrepel::geom_text_repel(
    data = line_label_all,
    aes(
      x = condition,
      y = mean_percent_all_mphi,
      label = cluster_label,
      color = cluster
    ),
    inherit.aes = FALSE,
    direction = "y",
    hjust = 0,
    nudge_x = 0.17,
    box.padding = 0.20,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.size = 0.28,
    seed = 4850,
    size = 3.7,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  facet_wrap(
    ~ macrophage_class,
    scales = "free_y",
    ncol = 2,
    labeller = as_labeller(
      CLASS_LABELS
    )
  ) +
  scale_color_manual(
    values = cluster_palette
  ) +
  scale_x_discrete(
    expand = expansion(
      mult = c(
        0.05,
        0.24
      )
    )
  ) +
  labs(
    title =
      "Res2 cluster changes within each MΦ class",
    subtitle =
      "Y-axis = cluster fraction of ALL MΦ cells | cluster number shown at Tx endpoint",
    x = NULL,
    y = "Fraction of total MΦ (%)"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold"
      ),
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      ),
    legend.position =
      "none",
    plot.margin =
      margin(
        5.5,
        35,
        5.5,
        5.5
      )
  ) +
  coord_cartesian(
    clip = "off"
  )

save_pdf(
  "03_cluster_trajectory_allMphi_v4.8.5.pdf",
  p_slope_all,
  14,
  12
)

save_png_600(
  "03_cluster_trajectory_allMphi_v4.8.5_600dpi.png",
  p_slope_all,
  14,
  12
)

# ------------------------------------------------------------------------------
# 10. Class-faceted trajectory plot
#     cluster fraction WITHIN CLASS
# ------------------------------------------------------------------------------

line_label_class <- condition_cluster %>%
  filter(
    condition == "Tx"
  ) %>%
  mutate(
    cluster_label = paste0(
      "#",
      as.character(cluster)
    )
  )

p_slope_class <- ggplot(
  condition_cluster,
  aes(
    x = condition,
    y = mean_percent_within_class,
    group = cluster,
    color = cluster
  )
) +
  geom_line(
    linewidth = 0.95,
    alpha = 0.90
  ) +
  geom_point(
    size = 2.8
  ) +
  ggrepel::geom_text_repel(
    data = line_label_class,
    aes(
      x = condition,
      y = mean_percent_within_class,
      label = cluster_label,
      color = cluster
    ),
    inherit.aes = FALSE,
    direction = "y",
    hjust = 0,
    nudge_x = 0.17,
    box.padding = 0.20,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.size = 0.28,
    seed = 4850,
    size = 3.7,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  facet_wrap(
    ~ macrophage_class,
    scales = "free_y",
    ncol = 2,
    labeller = as_labeller(
      CLASS_LABELS
    )
  ) +
  scale_color_manual(
    values = cluster_palette
  ) +
  scale_x_discrete(
    expand = expansion(
      mult = c(
        0.05,
        0.24
      )
    )
  ) +
  labs(
    title =
      "Res2 cluster composition within each MΦ class",
    subtitle =
      "Y-axis = cluster fraction WITHIN assigned class | cluster number shown at Tx endpoint",
    x = NULL,
    y = "Fraction within MΦ class (%)"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold"
      ),
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      ),
    legend.position =
      "none",
    plot.margin =
      margin(
        5.5,
        35,
        5.5,
        5.5
      )
  ) +
  coord_cartesian(
    clip = "off"
  )

save_pdf(
  "04_cluster_trajectory_withinClass_v4.8.5.pdf",
  p_slope_class,
  14,
  12
)

# ------------------------------------------------------------------------------
# 11. Group-wise cluster fraction plot
#     Each condition shown separately
# ------------------------------------------------------------------------------

p_group_fraction <- ggplot(
  condition_cluster,
  aes(
    x = cluster,
    y = mean_percent_all_mphi,
    fill = macrophage_class
  )
) +
  geom_col(
    width = 0.78
  ) +
  facet_wrap(
    ~ condition,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = CLASS_COLORS,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  labs(
    title =
      "Res2 cluster distribution by condition",
    subtitle =
      "Mean cluster fraction of ALL MΦ cells",
    x = "Res2 cluster",
    y = "Fraction of total MΦ (%)",
    fill = "MΦ class"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
    legend.position =
      "right"
  )

save_pdf(
  "05_cluster_fraction_by_condition_v4.8.5.pdf",
  p_group_fraction,
  15,
  10
)

save_png_600(
  "05_cluster_fraction_by_condition_v4.8.5_600dpi.png",
  p_group_fraction,
  15,
  10
)

# ------------------------------------------------------------------------------
# 12. Group-wise cluster CELL COUNT plot
# ------------------------------------------------------------------------------

p_group_count <- ggplot(
  condition_cluster,
  aes(
    x = cluster,
    y = mean_cells,
    fill = macrophage_class
  )
) +
  geom_col(
    width = 0.78
  ) +
  facet_wrap(
    ~ condition,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = CLASS_COLORS,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(
    title =
      "Res2 cluster cell counts by condition",
    subtitle =
      "Mean recovered cells per biological sample",
    x = "Res2 cluster",
    y = "Cell count",
    fill = "MΦ class"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
    legend.position =
      "right"
  )

save_pdf(
  "06_cluster_cellcount_by_condition_v4.8.5.pdf",
  p_group_count,
  15,
  10
)

# ------------------------------------------------------------------------------
# 13. One-panel grouped bar plot
#     STD / CDAHFD / Sham / Tx side-by-side per cluster
# ------------------------------------------------------------------------------

p_cluster_condition_dodge <- ggplot(
  condition_cluster,
  aes(
    x = cluster,
    y = mean_percent_all_mphi,
    fill = condition
  )
) +
  geom_col(
    position =
      position_dodge(
        width = 0.82
      ),
    width = 0.74
  ) +
  scale_fill_manual(
    values = CONDITION_COLORS
  ) +
  labs(
    title =
      "Res2 clusters: STD / CDAHFD / Sham / Tx comparison",
    subtitle =
      "Mean fraction of ALL MΦ cells",
    x = "Res2 cluster",
    y = "Fraction of total MΦ (%)",
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    legend.position =
      "top"
  )

save_pdf(
  "07_cluster_fraction_grouped_4conditions_v4.8.5.pdf",
  p_cluster_condition_dodge,
  16,
  7
)

# ------------------------------------------------------------------------------
# 14. Biological sample-level cluster fraction
# ------------------------------------------------------------------------------

p_sample_fraction <- ggplot(
  sample_cluster,
  aes(
    x = cluster,
    y = percent_all_mphi,
    group = sample,
    color = condition
  )
) +
  geom_point(
    size = 1.9,
    alpha = 0.80,
    position =
      position_jitter(
        width = 0.08,
        height = 0
      )
  ) +
  facet_wrap(
    ~ condition,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = CONDITION_COLORS
  ) +
  labs(
    title =
      "Biological sample-level Res2 cluster fractions",
    subtitle =
      "Each point = biological sample",
    x = "Res2 cluster",
    y = "Fraction of total MΦ (%)",
    color = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold"
      ),
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7
      ),
    legend.position =
      "none"
  )

save_pdf(
  "08_cluster_fraction_biological_samples_v4.8.5.pdf",
  p_sample_fraction,
  15,
  10
)

# ------------------------------------------------------------------------------
# 15. Disease / Treatment delta table
# ------------------------------------------------------------------------------

wide_condition <- condition_cluster %>%
  select(
    macrophage_class,
    cluster,
    condition,
    mean_percent_all_mphi,
    mean_cells
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(
      mean_percent_all_mphi,
      mean_cells
    ),
    values_fill = 0
  ) %>%
  mutate(
    delta_Disease_percent =
      mean_percent_all_mphi_CDAHFD -
      mean_percent_all_mphi_STD,

    delta_Treatment_percent =
      mean_percent_all_mphi_Tx -
      mean_percent_all_mphi_Sham,

    delta_Disease_cells =
      mean_cells_CDAHFD -
      mean_cells_STD,

    delta_Treatment_cells =
      mean_cells_Tx -
      mean_cells_Sham
  )

write.csv(
  wide_condition,
  file.path(
    CSV_OUT_DIR,
    "03_cluster_Disease_Treatment_changes_v4.8.5.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 16. Disease / treatment delta plot
# ------------------------------------------------------------------------------

delta_long <- bind_rows(

  wide_condition %>%
    transmute(
      macrophage_class,
      cluster,
      comparison = "STD -> CDAHFD",
      delta_percent =
        delta_Disease_percent
    ),

  wide_condition %>%
    transmute(
      macrophage_class,
      cluster,
      comparison = "Sham -> Tx",
      delta_percent =
        delta_Treatment_percent
    )
)

delta_long$comparison <- factor(
  delta_long$comparison,
  levels = c(
    "STD -> CDAHFD",
    "Sham -> Tx"
  )
)

p_delta <- ggplot(
  delta_long,
  aes(
    x = cluster,
    y = delta_percent,
    fill = comparison
  )
) +
  geom_col(
    position =
      position_dodge(
        width = 0.80
      ),
    width = 0.72
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.50
  ) +
  facet_wrap(
    ~ macrophage_class,
    ncol = 2,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "STD -> CDAHFD" = "#B2182B",
      "Sham -> Tx" = "#2166AC"
    )
  ) +
  labs(
    title =
      "Cluster-level disease and treatment changes",
    subtitle =
      "Positive = increase | negative = decrease",
    x = "Res2 cluster",
    y =
      "Δ fraction of total MΦ (percentage points)",
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    strip.text =
      element_text(
        face = "bold"
      ),
    legend.position =
      "top"
  )

save_pdf(
  "09_cluster_delta_Disease_Treatment_v4.8.5.pdf",
  p_delta,
  14,
  10
)

# ------------------------------------------------------------------------------
# 17. Summary figure
# ------------------------------------------------------------------------------

p_summary <- (
  p_slope_all /
    p_group_fraction /
    p_group_count /
    p_delta
) +
  patchwork::plot_layout(
    heights = c(
      1.4,
      1.0,
      1.0,
      1.0
    )
  ) +
  patchwork::plot_annotation(
    title =
      "Mouse MASH MΦ Res2 cluster dynamics v4.8.5",
    subtitle =
      "Manual macrophage classes v4.8.4 fixed",
    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 19
          ),
        plot.subtitle =
          element_text(
            size = 12
          )
      )
  )

save_pdf(
  "10_Mphi_Res2_cluster_dynamics_summary_v4.8.5.pdf",
  p_summary,
  17,
  30
)

save_jpeg_600(
  "10_Mphi_Res2_cluster_dynamics_summary_v4.8.5_600dpi.jpg",
  p_summary,
  17,
  30
)

# ------------------------------------------------------------------------------
# 18. README
# ------------------------------------------------------------------------------

summary_lines <- c(
  "Mouse MASH MΦ Res2 cluster dynamics v4.8.5",
  "",
  paste0(
    "Input: ",
    INPUT_RDS
  ),
  "",
  "v4.8.4 manual annotation is fixed.",
  "",
  "Main plots:",
  "  03 class-faceted cluster trajectories, all-MΦ denominator",
  "  04 class-faceted trajectories, within-class denominator",
  "  05 group-wise cluster fraction",
  "  06 group-wise cluster cell count",
  "  07 four-condition grouped cluster fraction",
  "  08 biological-sample cluster fraction",
  "  09 disease/treatment delta",
  "  10 summary",
  "",
  "Comparisons:",
  "  Disease   = STD -> CDAHFD",
  "  Treatment = Sham -> Tx"
)

writeLines(
  summary_lines,
  file.path(
    OUTPUT_DIR,
    "README_Mphi_Res2_cluster_dynamics_v4.8.5.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo_v4.8.5.txt"
  )
)

message_time("DONE.")
message_time(
  "Output: ",
  OUTPUT_DIR
)
