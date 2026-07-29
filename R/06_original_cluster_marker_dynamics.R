# ============================================================
# Original High-resolution Cluster Marker Dynamics
# uenoy scRNAseq Framework
#
# Purpose:
#   1. Res3.0 cluster 1-40 の検体間構成比変化
#   2. Res3.0 clusterごとのマーカー平均発現変化
#   3. 発現細胞率変化
#   4. Marker module score変化
#   5. Cluster番号を線上に直接表示
#   6. 集約群版と個別検体版の両方を出力
#
# Expected sample order:
#   Individual:
#     STD -> CDAHFD -> Sham1 -> Sham20 -> Tx5 -> Tx17
#
#   Grouped:
#     STD -> CDAHFD -> Sham -> Tx
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(scales)
})

# ============================================================
# 1. User settings
# ============================================================

project_dir <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

## Automatically detect the validated RDS
rds_candidates <- list.files(
  file.path(project_dir, "results"),
  pattern = "Mouse_Mphi_cluster1_40_original_annotation_validated.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(rds_candidates) == 0) {
  stop(
    "Mouse_Mphi_cluster1_40_original_annotation_validated.rds was not found under:\n",
    file.path(project_dir, "results")
  )
}

if (length(rds_candidates) > 1) {

  message("Multiple candidate RDS files detected:")

  for (i in seq_along(rds_candidates)) {
    message(sprintf("[%d] %s", i, rds_candidates[i]))
  }

  message("\nUsing the first candidate.")
}

rds_path <- rds_candidates[1]

message("Using RDS:")
message(rds_path)

cat("\n==============================\n")
cat("RDS auto detection succeeded\n")
cat("==============================\n")
print(rds_path)

out_dir <- file.path(
  project_dir,
  "Original_HighResolution_Cluster_Marker_Dynamics"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_dir <- file.path(out_dir, "PDF")
csv_dir <- file.path(out_dir, "CSV")

dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ============================================================
# 2. Plot settings
# ============================================================

base_family <- "Helvetica"

theme_ueno <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = base_size,
        hjust = 0.5
      ),
      strip.text = element_text(
        face = "bold",
        size = base_size - 1
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

annotation_colors <- c(
  "Resident Kupffer-like Mphi" = "#20B8C0",
  "Monocyte-like Mphi" = "#18B83E",
  "Inflammatory M1-like Mphi" = "#66A63C",
  "Pro-resolution M2-like Mphi" = "#6B9EF2",
  "SPP1/TREM2 MASH-associated Mphi" = "#E56CEB",
  "Fibrosis-associated macrophage" = "#E67C9D",
  "Other macrophage" = "#68B8B6"
)

sample_colors <- c(
  "STD" = "#4472C4",
  "CDAHFD" = "#ED7D31",
  "Sham1" = "#70AD47",
  "Sham20" = "#A5A5A5",
  "Tx5" = "#FFC000",
  "Tx17" = "#5B9BD5",
  "Sham" = "#70AD47",
  "Tx" = "#FFC000"
)

individual_order <- c(
  "STD",
  "CDAHFD",
  "Sham1",
  "Sham20",
  "Tx5",
  "Tx17"
)

grouped_order <- c(
  "STD",
  "CDAHFD",
  "Sham",
  "Tx"
)

# ============================================================
# 3. Marker definitions
# ============================================================

marker_sets <- list(

  Pan_macrophage = c(
    "Adgre1", "Csf1r", "Cd68", "Lyz2",
    "Fcgr1", "C1qa", "C1qb", "C1qc"
  ),

  Resident_Kupffer_identity = c(
    "Clec4f", "Timd4", "Vsig4", "Marco",
    "Cd5l", "Folr2", "Spic", "Id3"
  ),

  Monocyte_identity = c(
    "Ly6c2", "Ccr2", "S100a8", "S100a9",
    "Plac8", "Lgals3", "Chil3"
  ),

  Inflammatory_M1 = c(
    "Il1b", "Tnf", "Nfkbia", "Ccl2",
    "Ccl3", "Ccl4", "Cxcl9", "Cxcl10",
    "Cd80", "Cd86", "Ptgs2", "Nos2"
  ),

  Pro_resolution_M2 = c(
    "Mrc1", "Cd163", "Folr2", "Maf",
    "Mertk", "Il10", "Arg1", "Retnla",
    "Chil3", "C1qa", "C1qb", "C1qc"
  ),

  SPP1_TREM2_MASH = c(
    "Spp1", "Trem2", "Gpnmb", "Lgals3",
    "Cd9", "Fabp5", "Ctsb", "Ctsd",
    "Lpl", "Apoe", "Itgax"
  ),

  Efferocytosis_phagocytosis = c(
    "Mertk", "Axl", "Tyro3", "Gas6",
    "Mfge8", "Lrp1", "Abca1", "Cd36",
    "Fcgr1", "Ctsb", "Ctsd", "Ctsl"
  ),

  Fibrosis_associated = c(
    "Tgfb1", "Spp1", "Pdgfb", "Mmp9",
    "Mmp12", "Mmp14", "Timp1", "Ccl2",
    "Ccl7"
  ),

  IL10_response = c(
    "Il10ra", "Il10rb", "Stat3", "Socs3",
    "Bcl3", "Dusp1", "Klf4", "Maf"
  ),

  Lipid_lysosome = c(
    "Lpl", "Apoe", "Fabp5", "Cd36",
    "Abca1", "Ctsb", "Ctsd", "Ctsl",
    "Lamp1", "Lamp2"
  )
)

priority_genes <- unique(c(
  marker_sets$Resident_Kupffer_identity,
  marker_sets$Monocyte_identity,
  marker_sets$Inflammatory_M1,
  marker_sets$Pro_resolution_M2,
  marker_sets$SPP1_TREM2_MASH,
  marker_sets$Fibrosis_associated
))

# ============================================================
# 4. Utility functions
# ============================================================

detect_metadata_column <- function(
    meta,
    candidates,
    required = TRUE,
    label = "metadata column"
) {

  exact_hit <- candidates[candidates %in% colnames(meta)]

  if (length(exact_hit) > 0) {
    return(exact_hit[[1]])
  }

  lower_cols <- tolower(colnames(meta))
  lower_candidates <- tolower(candidates)

  idx <- match(lower_candidates, lower_cols)
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    return(colnames(meta)[idx[[1]]])
  }

  if (required) {
    stop(
      "Could not detect ", label, ".\n",
      "Candidates:\n  ",
      paste(candidates, collapse = "\n  "),
      "\nAvailable metadata columns:\n  ",
      paste(colnames(meta), collapse = "\n  ")
    )
  }

  return(NULL)
}

standardize_sample_name <- function(x) {

  x <- as.character(x)

  case_when(
    str_detect(x, regex("STD", ignore_case = TRUE)) ~ "STD",
    str_detect(x, regex("CDHFD|CDAHFD", ignore_case = TRUE)) ~ "CDAHFD",
    str_detect(x, regex("Sham.?1$", ignore_case = TRUE)) ~ "Sham1",
    str_detect(x, regex("Sham.?20$", ignore_case = TRUE)) ~ "Sham20",
    str_detect(x, regex("Tx.?5$", ignore_case = TRUE)) ~ "Tx5",
    str_detect(x, regex("Tx.?17$", ignore_case = TRUE)) ~ "Tx17",
    TRUE ~ x
  )
}

make_grouped_sample <- function(x) {

  case_when(
    x == "STD" ~ "STD",
    x == "CDAHFD" ~ "CDAHFD",
    x %in% c("Sham1", "Sham20") ~ "Sham",
    x %in% c("Tx5", "Tx17") ~ "Tx",
    TRUE ~ x
  )
}

normalize_layer2_name <- function(x) {

  x <- as.character(x)

  case_when(
    str_detect(x, regex("Resident", ignore_case = TRUE)) ~
      "Resident Kupffer-like Mphi",

    str_detect(x, regex("Monocyte", ignore_case = TRUE)) ~
      "Monocyte-like Mphi",

    str_detect(x, regex("Inflammatory|M1", ignore_case = TRUE)) ~
      "Inflammatory M1-like Mphi",

    str_detect(x, regex("Pro-resolution|M2", ignore_case = TRUE)) ~
      "Pro-resolution M2-like Mphi",

    str_detect(x, regex("SPP1|TREM2|MASH", ignore_case = TRUE)) ~
      "SPP1/TREM2 MASH-associated Mphi",

    str_detect(x, regex("Fibrosis", ignore_case = TRUE)) ~
      "Fibrosis-associated macrophage",

    TRUE ~ "Other macrophage"
  )
}

clean_cluster_number <- function(x) {

  x <- as.character(x)

  suppressWarnings({
    numeric_value <- as.numeric(
      str_extract(x, "[0-9]+")
    )
  })

  numeric_value
}

safe_scale <- function(x) {

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  s <- sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
}

get_assay_data_safe <- function(object, assay, layer = "data") {

  DefaultAssay(object) <- assay

  out <- tryCatch(
    SeuratObject::LayerData(
      object = object,
      assay = assay,
      layer = layer
    ),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    return(out)
  }

  out <- tryCatch(
    Seurat::GetAssayData(
      object = object,
      assay = assay,
      slot = layer
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    stop(
      "Could not retrieve assay data from assay = ",
      assay,
      ", layer/slot = ",
      layer
    )
  }

  out
}

label_at_last_sample <- function(df, sample_levels) {

  df %>%
    filter(sample_plot %in% sample_levels) %>%
    group_by(cluster_num) %>%
    arrange(
      factor(sample_plot, levels = sample_levels),
      .by_group = TRUE
    ) %>%
    slice_tail(n = 1) %>%
    ungroup()
}

save_plot_pdf <- function(plot, filename, width, height) {

  ggsave(
    filename = file.path(pdf_dir, filename),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

# ============================================================
# 5. Load object
# ============================================================

stopifnot(file.exists(rds_path))

obj <- readRDS(rds_path)

message("Loaded object: ", rds_path)
message("Cells: ", ncol(obj))
message("Features: ", nrow(obj))

meta <- obj@meta.data

# ============================================================
# 6. Detect metadata columns
# ============================================================

sample_col <- detect_metadata_column(
  meta,
  candidates = c(
    "sample",
    "Sample",
    "sample_id",
    "sample_name",
    "orig.ident",
    "dataset",
    "Dataset"
  ),
  label = "sample column"
)

cluster_col <- detect_metadata_column(
  meta,
  candidates = c(
    "RNA_snn_res.3",
    "RNA_snn_res.3.0",
    "SCT_snn_res.3",
    "SCT_snn_res.3.0",
    "integrated_snn_res.3",
    "integrated_snn_res.3.0",
    "seurat_clusters_res3.0",
    "cluster_res3.0",
    "res3.0",
    "res3",
    "seurat_clusters"
  ),
  label = "high-resolution cluster column"
)

layer2_col <- detect_metadata_column(
  meta,
  candidates = c(
    "layer2",
    "Layer2",
    "annotation_layer2",
    "celltype_layer2",
    "layer2_annotation",
    "Mphi_layer2",
    "annotation"
  ),
  label = "Layer2 annotation column"
)

message("Detected sample column: ", sample_col)
message("Detected cluster column: ", cluster_col)
message("Detected Layer2 column: ", layer2_col)

# ============================================================
# 7. Prepare metadata
# ============================================================

obj$sample_standard <- standardize_sample_name(
  obj@meta.data[[sample_col]]
)

obj$sample_grouped <- make_grouped_sample(
  obj$sample_standard
)

obj$cluster_res3 <- as.character(
  obj@meta.data[[cluster_col]]
)

obj$cluster_num <- clean_cluster_number(
  obj$cluster_res3
)

obj$layer2_standard <- normalize_layer2_name(
  obj@meta.data[[layer2_col]]
)

obj$sample_standard <- factor(
  obj$sample_standard,
  levels = individual_order
)

obj$sample_grouped <- factor(
  obj$sample_grouped,
  levels = grouped_order
)

obj$layer2_standard <- factor(
  obj$layer2_standard,
  levels = names(annotation_colors)
)

# Keep cluster 1-40 only
obj <- subset(
  obj,
  subset = !is.na(cluster_num) &
    cluster_num >= 1 &
    cluster_num <= 40
)

message("Cells after cluster 1-40 filtering: ", ncol(obj))

# ============================================================
# 8. Determine assay and available genes
# ============================================================

assay_candidates <- c("RNA", "SCT", DefaultAssay(obj))
assay_candidates <- unique(
  assay_candidates[assay_candidates %in% Assays(obj)]
)

if (length(assay_candidates) == 0) {
  stop("No suitable assay found.")
}

assay_use <- assay_candidates[[1]]
DefaultAssay(obj) <- assay_use

expr_mat <- get_assay_data_safe(
  object = obj,
  assay = assay_use,
  layer = "data"
)

available_genes <- rownames(expr_mat)

marker_sets_present <- map(
  marker_sets,
  ~ intersect(.x, available_genes)
)

marker_sets_missing <- map(
  marker_sets,
  ~ setdiff(.x, available_genes)
)

marker_availability <- tibble(
  marker_set = names(marker_sets),
  n_requested = map_int(marker_sets, length),
  n_present = map_int(marker_sets_present, length),
  present_genes = map_chr(
    marker_sets_present,
    ~ paste(.x, collapse = ";")
  ),
  missing_genes = map_chr(
    marker_sets_missing,
    ~ paste(.x, collapse = ";")
  )
)

write_csv(
  marker_availability,
  file.path(csv_dir, "Marker_availability.csv")
)

priority_genes_present <- intersect(
  priority_genes,
  available_genes
)

message(
  "Priority genes available: ",
  length(priority_genes_present),
  " / ",
  length(priority_genes)
)

# ============================================================
# 9. Add module scores
# ============================================================

set.seed(1234)

module_score_names <- character(0)

for (set_name in names(marker_sets_present)) {

  genes_use <- marker_sets_present[[set_name]]

  if (length(genes_use) < 2) {
    warning(
      "Skipping module score: ",
      set_name,
      " because fewer than 2 genes are available."
    )
    next
  }

  prefix <- paste0("MS_", set_name, "_")

  obj <- AddModuleScore(
    object = obj,
    features = list(genes_use),
    assay = assay_use,
    name = prefix,
    search = FALSE,
    seed = 1234
  )

  generated_col <- paste0(prefix, "1")

  module_score_names[set_name] <- generated_col
}

message(
  "Module scores generated: ",
  paste(names(module_score_names), collapse = ", ")
)

# ============================================================
# 10. Cell count and fraction summaries
# ============================================================

meta2 <- obj@meta.data %>%
  mutate(
    cell_barcode = rownames(obj@meta.data),
    sample_standard = as.character(sample_standard),
    sample_grouped = as.character(sample_grouped),
    cluster_num = as.numeric(cluster_num),
    cluster_label = paste0("C", cluster_num),
    layer2_standard = as.character(layer2_standard)
  )

cluster_annotation_map <- meta2 %>%
  count(
    cluster_num,
    cluster_label,
    layer2_standard,
    name = "n"
  ) %>%
  group_by(cluster_num, cluster_label) %>%
  mutate(
    fraction = n / sum(n)
  ) %>%
  slice_max(
    order_by = n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    cluster_num,
    cluster_label,
    dominant_layer2 = layer2_standard,
    dominant_fraction = fraction
  )

write_csv(
  cluster_annotation_map,
  file.path(csv_dir, "Cluster_to_Layer2_mapping.csv")
)

fraction_individual <- meta2 %>%
  count(
    sample_standard,
    cluster_num,
    cluster_label,
    layer2_standard,
    name = "n_cluster"
  ) %>%
  group_by(sample_standard) %>%
  mutate(
    n_total_mphi = sum(n_cluster),
    fraction_total_mphi = 100 * n_cluster / n_total_mphi
  ) %>%
  ungroup()

fraction_grouped <- meta2 %>%
  count(
    sample_grouped,
    cluster_num,
    cluster_label,
    layer2_standard,
    name = "n_cluster"
  ) %>%
  group_by(sample_grouped) %>%
  mutate(
    n_total_mphi = sum(n_cluster),
    fraction_total_mphi = 100 * n_cluster / n_total_mphi
  ) %>%
  ungroup()

write_csv(
  fraction_individual,
  file.path(csv_dir, "Cluster_fraction_individual_samples.csv")
)

write_csv(
  fraction_grouped,
  file.path(csv_dir, "Cluster_fraction_grouped_samples.csv")
)

# ============================================================
# 11. Generic line plot function
# ============================================================

plot_cluster_dynamics <- function(
    data,
    sample_col_name,
    value_col,
    sample_levels,
    title,
    subtitle,
    y_label,
    facet_scales = "free_y",
    label_clusters = TRUE,
    line_alpha = 0.80,
    point_size = 2.0,
    line_width = 0.7
) {

  df <- data %>%
    mutate(
      sample_plot = factor(
        .data[[sample_col_name]],
        levels = sample_levels
      ),
      cluster_label = paste0("C", cluster_num)
    ) %>%
    filter(!is.na(sample_plot))

  label_df <- label_at_last_sample(
    df = df,
    sample_levels = sample_levels
  )

  p <- ggplot(
    df,
    aes(
      x = sample_plot,
      y = .data[[value_col]],
      group = cluster_num,
      color = layer2_standard
    )
  ) +
    geom_line(
      linewidth = line_width,
      alpha = line_alpha,
      na.rm = TRUE
    ) +
    geom_point(
      size = point_size,
      alpha = 0.95,
      na.rm = TRUE
    ) +
    facet_wrap(
      vars(layer2_standard),
      scales = facet_scales,
      ncol = 2
    ) +
    scale_color_manual(
      values = annotation_colors,
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(mult = c(0.05, 0.20))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_label
    ) +
    theme_ueno()

  if (label_clusters) {

    p <- p +
      geom_text_repel(
        data = label_df,
        aes(
          label = cluster_label
        ),
        size = 2.8,
        direction = "y",
        hjust = 0,
        nudge_x = 0.25,
        segment.size = 0.25,
        segment.alpha = 0.6,
        min.segment.length = 0,
        box.padding = 0.20,
        point.padding = 0.15,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }

  p
}

# ============================================================
# 12. Existing-style fraction plots
# ============================================================

p_fraction_grouped <- plot_cluster_dynamics(
  data = fraction_grouped,
  sample_col_name = "sample_grouped",
  value_col = "fraction_total_mphi",
  sample_levels = grouped_order,
  title = "Cluster fraction transition by macrophage annotation",
  subtitle = "Each line = Res3.0 cluster; y-axis free by annotation group",
  y_label = "Fraction among total Mphi (%)",
  facet_scales = "free_y"
)

save_plot_pdf(
  p_fraction_grouped,
  "01_Cluster1_40_fraction_total_Mphi_grouped_with_labels.pdf",
  width = 14,
  height = 10
)

p_fraction_individual <- plot_cluster_dynamics(
  data = fraction_individual,
  sample_col_name = "sample_standard",
  value_col = "fraction_total_mphi",
  sample_levels = individual_order,
  title = "Cluster fraction transition by individual sample",
  subtitle = "Each line = Res3.0 cluster; y-axis free by annotation group",
  y_label = "Fraction among total Mphi (%)",
  facet_scales = "free_y"
)

save_plot_pdf(
  p_fraction_individual,
  "02_Cluster1_40_fraction_total_Mphi_individual_with_labels.pdf",
  width = 16,
  height = 10
)

# ============================================================
# 13. Average expression and percent expressing
# ============================================================

group_vars <- c(
  "sample_standard",
  "sample_grouped",
  "cluster_num",
  "layer2_standard"
)

cell_meta_expr <- meta2 %>%
  select(
    cell_barcode,
    all_of(group_vars)
  )

summarize_gene_expression <- function(
    gene,
    expr_mat,
    cell_meta
) {

  gene_values <- as.numeric(
    expr_mat[gene, cell_meta$cell_barcode, drop = TRUE]
  )

  cell_meta %>%
    mutate(
      expression = gene_values,
      expressed = expression > 0
    ) %>%
    group_by(
      sample_standard,
      sample_grouped,
      cluster_num,
      layer2_standard
    ) %>%
    summarise(
      gene = gene,
      n_cells = n(),
      avg_expression = mean(expression, na.rm = TRUE),
      pct_expressing = 100 * mean(expressed, na.rm = TRUE),
      median_expression = median(expression, na.rm = TRUE),
      .groups = "drop"
    )
}

message("Summarizing individual genes...")

gene_summary <- map_dfr(
  priority_genes_present,
  summarize_gene_expression,
  expr_mat = expr_mat,
  cell_meta = cell_meta_expr
)

write_csv(
  gene_summary,
  file.path(csv_dir, "Original_cluster_gene_expression_summary.csv")
)

# ============================================================
# 14. Per-gene line plots
# ============================================================

gene_plot_dir <- file.path(pdf_dir, "Individual_genes")
dir.create(gene_plot_dir, recursive = TRUE, showWarnings = FALSE)

for (gene_i in priority_genes_present) {

  gene_df <- gene_summary %>%
    filter(gene == gene_i)

  p_avg_grouped <- plot_cluster_dynamics(
    data = gene_df,
    sample_col_name = "sample_grouped",
    value_col = "avg_expression",
    sample_levels = grouped_order,
    title = paste0(gene_i, ": average expression"),
    subtitle = "Each line = Res3.0 cluster",
    y_label = "Average normalized expression",
    facet_scales = "free_y"
  )

  p_pct_grouped <- plot_cluster_dynamics(
    data = gene_df,
    sample_col_name = "sample_grouped",
    value_col = "pct_expressing",
    sample_levels = grouped_order,
    title = paste0(gene_i, ": percent expressing"),
    subtitle = "Each line = Res3.0 cluster",
    y_label = "Percent expressing (%)",
    facet_scales = "free_y"
  )

  p_avg_individual <- plot_cluster_dynamics(
    data = gene_df,
    sample_col_name = "sample_standard",
    value_col = "avg_expression",
    sample_levels = individual_order,
    title = paste0(gene_i, ": average expression by individual sample"),
    subtitle = "Each line = Res3.0 cluster",
    y_label = "Average normalized expression",
    facet_scales = "free_y"
  )

  p_pct_individual <- plot_cluster_dynamics(
    data = gene_df,
    sample_col_name = "sample_standard",
    value_col = "pct_expressing",
    sample_levels = individual_order,
    title = paste0(gene_i, ": percent expressing by individual sample"),
    subtitle = "Each line = Res3.0 cluster",
    y_label = "Percent expressing (%)",
    facet_scales = "free_y"
  )

  combined_grouped <- p_avg_grouped / p_pct_grouped

  combined_individual <- p_avg_individual / p_pct_individual

  ggsave(
    filename = file.path(
      gene_plot_dir,
      paste0(
        "Gene_",
        gene_i,
        "_grouped_avg_and_pct.pdf"
      )
    ),
    plot = combined_grouped,
    device = cairo_pdf,
    width = 14,
    height = 18,
    units = "in",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      gene_plot_dir,
      paste0(
        "Gene_",
        gene_i,
        "_individual_avg_and_pct.pdf"
      )
    ),
    plot = combined_individual,
    device = cairo_pdf,
    width = 16,
    height = 18,
    units = "in",
    limitsize = FALSE
  )
}

# ============================================================
# 15. Module score summaries
# ============================================================

if (length(module_score_names) > 0) {

  module_meta <- obj@meta.data %>%
    mutate(
      cell_barcode = rownames(obj@meta.data),
      sample_standard = as.character(sample_standard),
      sample_grouped = as.character(sample_grouped),
      cluster_num = as.numeric(cluster_num),
      layer2_standard = as.character(layer2_standard)
    )

  module_summary <- map_dfr(
    names(module_score_names),
    function(module_name) {

      score_col <- module_score_names[[module_name]]

      module_meta %>%
        group_by(
          sample_standard,
          sample_grouped,
          cluster_num,
          layer2_standard
        ) %>%
        summarise(
          module = module_name,
          n_cells = n(),
          mean_score = mean(
            .data[[score_col]],
            na.rm = TRUE
          ),
          median_score = median(
            .data[[score_col]],
            na.rm = TRUE
          ),
          .groups = "drop"
        )
    }
  )

  write_csv(
    module_summary,
    file.path(csv_dir, "Original_cluster_module_score_summary.csv")
  )

  module_plot_dir <- file.path(
    pdf_dir,
    "Module_scores"
  )

  dir.create(
    module_plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (module_i in unique(module_summary$module)) {

    module_df <- module_summary %>%
      filter(module == module_i)

    p_module_grouped <- plot_cluster_dynamics(
      data = module_df,
      sample_col_name = "sample_grouped",
      value_col = "mean_score",
      sample_levels = grouped_order,
      title = paste0(
        module_i,
        " module score transition"
      ),
      subtitle = "Each line = Res3.0 cluster",
      y_label = "Mean module score",
      facet_scales = "free_y"
    )

    p_module_individual <- plot_cluster_dynamics(
      data = module_df,
      sample_col_name = "sample_standard",
      value_col = "mean_score",
      sample_levels = individual_order,
      title = paste0(
        module_i,
        " module score by individual sample"
      ),
      subtitle = "Each line = Res3.0 cluster",
      y_label = "Mean module score",
      facet_scales = "free_y"
    )

    save_plot_pdf(
      p_module_grouped,
      paste0(
        "Module_",
        module_i,
        "_grouped.pdf"
      ),
      width = 14,
      height = 10
    )

    save_plot_pdf(
      p_module_individual,
      paste0(
        "Module_",
        module_i,
        "_individual.pdf"
      ),
      width = 16,
      height = 10
    )
  }
}

# ============================================================
# 16. Heatmap-style summaries
# ============================================================

gene_heatmap_data <- gene_summary %>%
  group_by(
    sample_standard,
    cluster_num,
    layer2_standard,
    gene
  ) %>%
  summarise(
    avg_expression = weighted.mean(
      avg_expression,
      w = pmax(n_cells, 1),
      na.rm = TRUE
    ),
    pct_expressing = weighted.mean(
      pct_expressing,
      w = pmax(n_cells, 1),
      na.rm = TRUE
    ),
    n_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    avg_expression_scaled = safe_scale(avg_expression)
  ) %>%
  ungroup() %>%
  mutate(
    sample_cluster = paste0(
      sample_standard,
      "_C",
      cluster_num
    ),
    sample_standard = factor(
      sample_standard,
      levels = individual_order
    )
  )

p_heatmap_avg <- ggplot(
  gene_heatmap_data,
  aes(
    x = sample_cluster,
    y = gene,
    fill = avg_expression_scaled
  )
) +
  geom_tile() +
  facet_grid(
    cols = vars(sample_standard),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Original high-resolution cluster marker expression",
    subtitle = "Gene-wise scaled average expression",
    x = "Sample and Res3.0 cluster",
    y = "Marker gene",
    fill = "Scaled\nexpression"
  ) +
  theme_ueno(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    axis.text.y = element_text(
      size = 7
    ),
    legend.position = "right"
  )

save_plot_pdf(
  p_heatmap_avg,
  "03_Original_cluster_marker_average_expression_heatmap.pdf",
  width = 24,
  height = 14
)

p_heatmap_pct <- ggplot(
  gene_heatmap_data,
  aes(
    x = sample_cluster,
    y = gene,
    fill = pct_expressing
  )
) +
  geom_tile() +
  facet_grid(
    cols = vars(sample_standard),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 100),
    oob = squish
  ) +
  labs(
    title = "Original high-resolution cluster marker detection rate",
    subtitle = "Percent expressing cells",
    x = "Sample and Res3.0 cluster",
    y = "Marker gene",
    fill = "Percent\nexpressing"
  ) +
  theme_ueno(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    axis.text.y = element_text(
      size = 7
    ),
    legend.position = "right"
  )

save_plot_pdf(
  p_heatmap_pct,
  "04_Original_cluster_marker_percent_expressing_heatmap.pdf",
  width = 24,
  height = 14
)

# ============================================================
# 17. Cluster marker stability
# ============================================================

marker_stability <- gene_summary %>%
  filter(
    sample_standard %in% individual_order
  ) %>%
  group_by(
    cluster_num,
    layer2_standard,
    gene
  ) %>%
  summarise(
    n_samples = n_distinct(sample_standard),
    mean_avg_expression = mean(
      avg_expression,
      na.rm = TRUE
    ),
    sd_avg_expression = sd(
      avg_expression,
      na.rm = TRUE
    ),
    cv_avg_expression = if_else(
      abs(mean_avg_expression) > 1e-8,
      sd_avg_expression / abs(mean_avg_expression),
      NA_real_
    ),
    mean_pct_expressing = mean(
      pct_expressing,
      na.rm = TRUE
    ),
    sd_pct_expressing = sd(
      pct_expressing,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  marker_stability,
  file.path(csv_dir, "Original_cluster_marker_stability.csv")
)

# ============================================================
# 18. Disease-response pattern classification
# ============================================================

grouped_gene_summary <- gene_summary %>%
  group_by(
    sample_grouped,
    cluster_num,
    layer2_standard,
    gene
  ) %>%
  summarise(
    avg_expression = weighted.mean(
      avg_expression,
      w = pmax(n_cells, 1),
      na.rm = TRUE
    ),
    pct_expressing = weighted.mean(
      pct_expressing,
      w = pmax(n_cells, 1),
      na.rm = TRUE
    ),
    n_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample_grouped,
    values_from = c(
      avg_expression,
      pct_expressing,
      n_cells
    )
  )

classify_pattern <- function(std, cdahfd, sham, tx) {

  vals <- c(std, cdahfd, sham, tx)

  if (sum(is.finite(vals)) < 3) {
    return("Insufficient data")
  }

  disease_up <- is.finite(std) &&
    is.finite(cdahfd) &&
    cdahfd > std

  disease_down <- is.finite(std) &&
    is.finite(cdahfd) &&
    cdahfd < std

  tx_restored_down <- disease_up &&
    is.finite(tx) &&
    tx < cdahfd

  tx_restored_up <- disease_down &&
    is.finite(tx) &&
    tx > cdahfd

  sham_specific <- is.finite(sham) &&
    sham > max(
      c(std, cdahfd, tx),
      na.rm = TRUE
    )

  tx_specific <- is.finite(tx) &&
    tx > max(
      c(std, cdahfd, sham),
      na.rm = TRUE
    )

  range_val <- diff(
    range(vals, na.rm = TRUE)
  )

  mean_abs <- mean(
    abs(vals),
    na.rm = TRUE
  )

  stable <- is.finite(range_val) &&
    range_val <= 0.20 * pmax(mean_abs, 0.1)

  case_when(
    tx_restored_down ~
      "CDAHFD-up / Tx-down",

    tx_restored_up ~
      "CDAHFD-down / Tx-restored",

    sham_specific ~
      "Sham-specific",

    tx_specific ~
      "Tx-specific",

    stable ~
      "Stable",

    TRUE ~
      "Inconsistent or complex"
  )
}

required_pattern_cols <- c(
  "avg_expression_STD",
  "avg_expression_CDAHFD",
  "avg_expression_Sham",
  "avg_expression_Tx"
)

for (nm in required_pattern_cols) {
  if (!nm %in% colnames(grouped_gene_summary)) {
    grouped_gene_summary[[nm]] <- NA_real_
  }
}

disease_pattern <- grouped_gene_summary %>%
  rowwise() %>%
  mutate(
    expression_pattern = classify_pattern(
      avg_expression_STD,
      avg_expression_CDAHFD,
      avg_expression_Sham,
      avg_expression_Tx
    )
  ) %>%
  ungroup()

write_csv(
  disease_pattern,
  file.path(csv_dir, "Original_cluster_gene_disease_patterns.csv")
)

# ============================================================
# 19. Cluster-level annotation confidence support table
# ============================================================

identity_modules <- c(
  "Resident_Kupffer_identity",
  "Monocyte_identity",
  "Inflammatory_M1",
  "Pro_resolution_M2",
  "SPP1_TREM2_MASH",
  "Fibrosis_associated"
)

if (exists("module_summary")) {

  cluster_module_consistency <- module_summary %>%
    filter(
      module %in% identity_modules
    ) %>%
    group_by(
      cluster_num,
      layer2_standard,
      module
    ) %>%
    summarise(
      n_samples = n_distinct(sample_standard),
      mean_score = mean(
        mean_score,
        na.rm = TRUE
      ),
      sd_score = sd(
        mean_score,
        na.rm = TRUE
      ),
      cv_score = if_else(
        abs(mean_score) > 1e-8,
        sd_score / abs(mean_score),
        NA_real_
      ),
      .groups = "drop"
    )

  write_csv(
    cluster_module_consistency,
    file.path(
      csv_dir,
      "Cluster_module_consistency_for_annotation.csv"
    )
  )
}

# ============================================================
# 20. Session information
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

message("================================================")
message("Analysis completed.")
message("Output:")
message(out_dir)
message("================================================")
