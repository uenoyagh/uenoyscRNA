#!/usr/bin/env Rscript

# ==============================================================================
# uenoyscRNA
# Mouse macrophage data-driven analysis Step 3
# Five broad macrophage classes + independent stress/QC axis
# Version 1.0.0
# ==============================================================================
#
# Primary method:
#   RPCA
#
# Sensitivity method:
#   Harmony
#
# Five broad classes:
#   1. M1-like Mphi
#   2. Homeostatic / pro-resolution Mphi
#   3. Fibrosis-associated Mphi
#   4. Monocyte / myeloid-like cells
#   5. Cycling Mphi
#
# Important:
#   - Low_quality_stress is NOT used as a biological class.
#   - It is stored independently as qc_state_internal.
#   - Neutrophil-like and dendritic-like states are retained in the fine-level
#     annotation and grouped provisionally into Monocyte / myeloid-like cells.
#   - These cells should later be reviewed for true contamination/doublets.
#
# INPUT:
#   Step 2 characterized RPCA/Harmony RDS files
#
# OUTPUT:
#   - RDS files with:
#       layer2_internal
#       broad_mphi_class_internal
#       qc_state_internal
#   - Broad-class UMAPs
#   - Sample-split UMAPs
#   - Broad-class composition plots
#   - Sham vs Tx abundance summaries
#   - Broad-class DotPlots and module-score summaries
#   - Excel workbook
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
BASE_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  "DataDriven_ShamTx_RPCA_vs_Harmony_v1.0.0"
)

INPUT_DIR <- file.path(
  BASE_DIR,
  "Step2_ClusterCharacterization_v1.0.0",
  "RDS"
)

RPCA_RDS <- file.path(
  INPUT_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_RPCA_cluster_characterized_v1.0.0.rds"
)

HARMONY_RDS <- file.path(
  INPUT_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_Harmony_cluster_characterized_v1.0.0.rds"
)

OUTPUT_DIR <- file.path(
  BASE_DIR,
  "Step3_FiveBroadClasses_v1.0.0"
)

FIGURE_DIR <- file.path(OUTPUT_DIR, "Figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "Tables")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

for (d in c(OUTPUT_DIR, FIGURE_DIR, TABLE_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ASSAY_USE <- "RNA"
RPCA_CLUSTER_COL <- "mphi_rpca_res_1_2"
HARMONY_CLUSTER_COL <- "mphi_harmony_res_1_2"

RPCA_REDUCTION <- "mphi.umap.rpca"
HARMONY_REDUCTION <- "mphi.umap.harmony"

FINE_ANNOTATION_COL <- "layer2_internal_provisional"

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# ------------------------------------------------------------------------------
# 2. Plot colors and factor order
# ------------------------------------------------------------------------------
broad_class_levels <- c(
  "M1-like Mphi",
  "Homeostatic / pro-resolution Mphi",
  "Fibrosis-associated Mphi",
  "Monocyte / myeloid-like cells",
  "Cycling Mphi"
)

broad_class_colors <- c(
  "M1-like Mphi" = "#E53935",
  "Homeostatic / pro-resolution Mphi" = "#00A6A6",
  "Fibrosis-associated Mphi" = "#D81B60",
  "Monocyte / myeloid-like cells" = "#3F51B5",
  "Cycling Mphi" = "#8E24AA"
)

condition_levels <- c("Sham", "Tx")
condition_colors <- c(
  "Sham" = "#6F6F6F",
  "Tx" = "#FF8C1A"
)

sample_levels <- c("Sham1", "Sham20", "Tx17", "Tx5")
sample_colors <- c(
  "Sham1" = "#4D4D4D",
  "Sham20" = "#9A9A9A",
  "Tx17" = "#FF7A00",
  "Tx5" = "#FFB347"
)

qc_state_levels <- c(
  "Stress-low",
  "Stress-intermediate",
  "Stress-high"
)

qc_state_colors <- c(
  "Stress-low" = "#BDBDBD",
  "Stress-intermediate" = "#FFA726",
  "Stress-high" = "#D32F2F"
)

# ------------------------------------------------------------------------------
# 3. Fine-to-broad mapping
# ------------------------------------------------------------------------------
# These labels must match the Step 2 provisional annotation labels exactly.
fine_to_broad_map <- c(
  "Inflammatory M1-like M." = "M1-like Mphi",

  "Efferocytosis-high M." = "Homeostatic / pro-resolution Mphi",
  "Pro-resolution M2-like M." = "Homeostatic / pro-resolution Mphi",
  "Resident Kupffer-like M." = "Homeostatic / pro-resolution Mphi",

  "Fibrosis-associated M." = "Fibrosis-associated Mphi",
  "SPP1/TREM2 MASH-associated M." = "Fibrosis-associated Mphi",

  "Monocyte-like M." = "Monocyte / myeloid-like cells",
  "Neutrophil-like contamination" = "Monocyte / myeloid-like cells",
  "Dendritic-like contamination" = "Monocyte / myeloid-like cells",

  "Cycling M." = "Cycling Mphi"
)

# "Other macrophage" remains unassigned by the automatic mapping.
# It will be assigned by its highest biological module score.
unassigned_fine_labels <- c(
  NA_character_,
  "",
  "Other macrophage"
)

# Biological score priority used only for unresolved "Other macrophage" cells.
# Cycling and contamination states are handled conservatively.
score_to_broad_map <- c(
  "Score_Inflammatory_M1_like" = "M1-like Mphi",
  "Score_Efferocytosis_high" = "Homeostatic / pro-resolution Mphi",
  "Score_Pro_resolution_M2_like" = "Homeostatic / pro-resolution Mphi",
  "Score_Resident_Kupffer" = "Homeostatic / pro-resolution Mphi",
  "Score_Fibrosis_associated" = "Fibrosis-associated Mphi",
  "Score_SPP1_TREM2_MASH" = "Fibrosis-associated Mphi",
  "Score_Monocyte_like" = "Monocyte / myeloid-like cells",
  "Score_Neutrophil_like" = "Monocyte / myeloid-like cells",
  "Score_Dendritic_like" = "Monocyte / myeloid-like cells",
  "Score_Cycling" = "Cycling Mphi"
)

# ------------------------------------------------------------------------------
# 4. Marker panels for figures
# ------------------------------------------------------------------------------
marker_panels <- list(
  M1_like = c(
    "Il1b", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Cxcl2",
    "Ptgs2", "Cd80", "Cd86", "Stat1", "Irf1"
  ),
  Homeostatic_ProResolution = c(
    "Clec4f", "Timd4", "Vsig4", "Marco", "Cd5l",
    "Mertk", "Axl", "Gas6", "Maf", "Mafb",
    "Il10", "Mrc1", "Cd163", "Folr2"
  ),
  Fibrosis_associated = c(
    "Spp1", "Trem2", "Gpnmb", "Cd9", "Lgals3",
    "Tgfb1", "Mmp12", "Mmp14", "Timp1", "Fn1", "Thbs1"
  ),
  Monocyte_Myeloid = c(
    "Ly6c2", "Ccr2", "Sell", "Plac8", "S100a8", "S100a9",
    "Fcgr3", "Ms4a7", "H2-Ab1", "Cd74", "Itgax",
    "Ly6g", "Csf3r", "Mpo"
  ),
  Cycling = c(
    "Mki67", "Top2a", "Stmn1", "Cenpf", "Cdk1", "Pcna"
  ),
  Stress_QC = c(
    "Fos", "Jun", "Junb", "Atf3", "Ddit3", "Hspa1a", "Hspa1b"
  )
)

# ------------------------------------------------------------------------------
# 5. Helpers
# ------------------------------------------------------------------------------
message_time <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

save_pdf <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(FIGURE_DIR, filename),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

safe_join_layers <- function(object, assay = "RNA") {
  DefaultAssay(object) <- assay
  layer_names <- Layers(object[[assay]])

  if (length(grep("^counts\\.", layer_names)) > 0 ||
      length(grep("^data\\.", layer_names)) > 0) {
    object <- JoinLayers(object, assay = assay)
  }

  object
}

existing_genes <- function(object, genes, assay = "RNA") {
  intersect(genes, rownames(object[[assay]]))
}

validate_metadata <- function(object, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(object@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

write_excel_sheets <- function(sheets, filename) {
  wb <- createWorkbook()

  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
    freezePane(wb, nm, firstRow = TRUE)

    if (ncol(sheets[[nm]]) > 0) {
      setColWidths(
        wb,
        nm,
        cols = seq_len(ncol(sheets[[nm]])),
        widths = "auto"
      )
    }
  }

  saveWorkbook(wb, filename, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 6. Assign five broad classes
# ------------------------------------------------------------------------------
assign_broad_class <- function(object) {

  validate_metadata(
    object,
    c(
      FINE_ANNOTATION_COL,
      "condition_internal",
      "sample_internal"
    )
  )

  fine <- as.character(object@meta.data[[FINE_ANNOTATION_COL]])
  broad <- unname(fine_to_broad_map[fine])

  unresolved <- is.na(broad) | fine %in% unassigned_fine_labels

  required_score_cols <- names(score_to_broad_map)
  available_score_cols <- intersect(
    required_score_cols,
    colnames(object@meta.data)
  )

  if (length(available_score_cols) == 0) {
    stop(
      "No biological module-score columns were found for unresolved cells."
    )
  }

  if (any(unresolved)) {
    score_matrix <- as.matrix(
      object@meta.data[unresolved, available_score_cols, drop = FALSE]
    )

    score_matrix[!is.finite(score_matrix)] <- NA_real_

    max_index <- apply(score_matrix, 1, function(x) {
      if (all(is.na(x))) return(NA_integer_)
      which.max(x)
    })

    inferred_score <- rep(NA_character_, length(max_index))
    valid <- !is.na(max_index)
    inferred_score[valid] <- colnames(score_matrix)[max_index[valid]]

    broad[unresolved] <- unname(score_to_broad_map[inferred_score])
  }

  # Last-resort conservative fallback
  broad[is.na(broad)] <- "Monocyte / myeloid-like cells"

  object$broad_mphi_class_internal <- factor(
    broad,
    levels = broad_class_levels
  )

  object
}

# ------------------------------------------------------------------------------
# 7. Independent stress/QC state
# ------------------------------------------------------------------------------
assign_stress_qc <- function(object) {

  stress_col <- "Score_Low_quality_stress"

  if (!(stress_col %in% colnames(object@meta.data))) {
    stop("Stress module score column not found: ", stress_col)
  }

  stress <- object@meta.data[[stress_col]]

  q <- quantile(
    stress,
    probs = c(0.33, 0.67),
    na.rm = TRUE,
    names = FALSE
  )

  qc_state <- case_when(
    is.na(stress) ~ NA_character_,
    stress <= q[1] ~ "Stress-low",
    stress <= q[2] ~ "Stress-intermediate",
    TRUE ~ "Stress-high"
  )

  object$qc_state_internal <- factor(
    qc_state,
    levels = qc_state_levels
  )

  object$stress_score_internal <- stress

  object
}

# ------------------------------------------------------------------------------
# 8. Summary tables
# ------------------------------------------------------------------------------
make_cell_annotation_table <- function(
    object,
    cluster_col,
    method_label) {

  object@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
      method = method_label,
      cell,
      condition_internal = as.character(condition_internal),
      sample_internal = as.character(sample_internal),
      cluster_internal = as.character(.data[[cluster_col]]),
      layer2_internal = as.character(.data[[FINE_ANNOTATION_COL]]),
      broad_mphi_class_internal = as.character(
        broad_mphi_class_internal
      ),
      qc_state_internal = as.character(qc_state_internal),
      stress_score_internal
    )
}

make_composition_table <- function(
    object,
    method_label) {

  object@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = method_label,
      condition_internal,
      sample_internal,
      broad_mphi_class_internal,
      name = "n_cells",
      .drop = FALSE
    ) %>%
    group_by(method, sample_internal) %>%
    mutate(
      sample_total = sum(n_cells),
      fraction_of_sample = if_else(
        sample_total > 0,
        n_cells / sample_total,
        NA_real_
      )
    ) %>%
    ungroup()
}

make_condition_summary <- function(composition_df) {
  composition_df %>%
    group_by(
      method,
      condition_internal,
      broad_mphi_class_internal
    ) %>%
    summarise(
      mean_fraction = mean(fraction_of_sample, na.rm = TRUE),
      sd_fraction = sd(fraction_of_sample, na.rm = TRUE),
      n_samples = n_distinct(sample_internal),
      .groups = "drop"
    )
}

make_tx_vs_sham_summary <- function(composition_df) {
  composition_df %>%
    group_by(
      method,
      condition_internal,
      broad_mphi_class_internal
    ) %>%
    summarise(
      mean_fraction = mean(fraction_of_sample, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = condition_internal,
      values_from = mean_fraction,
      values_fill = 0
    ) %>%
    mutate(
      log2FC_Tx_vs_Sham = log2((Tx + 1e-4) / (Sham + 1e-4)),
      direction = case_when(
        log2FC_Tx_vs_Sham > 0.25 ~ "Higher in Tx",
        log2FC_Tx_vs_Sham < -0.25 ~ "Lower in Tx",
        TRUE ~ "Similar"
      )
    )
}

make_cluster_to_broad_table <- function(
    object,
    cluster_col,
    method_label) {

  object@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = method_label,
      cluster_internal = .data[[cluster_col]],
      layer2_internal = .data[[FINE_ANNOTATION_COL]],
      broad_mphi_class_internal,
      qc_state_internal,
      name = "n_cells"
    ) %>%
    group_by(method, cluster_internal) %>%
    mutate(
      cluster_total = sum(n_cells),
      fraction_within_cluster = n_cells / cluster_total
    ) %>%
    ungroup()
}

# ------------------------------------------------------------------------------
# 9. Plotting helpers
# ------------------------------------------------------------------------------
plot_broad_umap <- function(
    object,
    reduction,
    method_label) {

  DimPlot(
    object,
    reduction = reduction,
    group.by = "broad_mphi_class_internal",
    cols = broad_class_colors,
    pt.size = 0.60,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    labs(
      title = paste0(
        method_label,
        ": five broad macrophage classes"
      ),
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 4)
      )
    )
}

plot_broad_umap_split <- function(
    object,
    reduction,
    method_label) {

  DimPlot(
    object,
    reduction = reduction,
    group.by = "broad_mphi_class_internal",
    split.by = "sample_internal",
    cols = broad_class_colors,
    pt.size = 0.42,
    raster = FALSE,
    ncol = 4
  ) +
    labs(
      title = paste0(
        method_label,
        ": broad classes split by sample"
      ),
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3.5)
      )
    )
}

plot_qc_umap <- function(
    object,
    reduction,
    method_label) {

  DimPlot(
    object,
    reduction = reduction,
    group.by = "qc_state_internal",
    cols = qc_state_colors,
    pt.size = 0.55,
    raster = FALSE
  ) +
    labs(
      title = paste0(
        method_label,
        ": independent stress/QC state"
      ),
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

# ------------------------------------------------------------------------------
# 10. Read objects
# ------------------------------------------------------------------------------
for (f in c(RPCA_RDS, HARMONY_RDS)) {
  if (!file.exists(f)) {
    stop("Input RDS not found: ", f)
  }
}

message_time("Reading Step 2 RPCA object.")
rpca <- readRDS(RPCA_RDS)

message_time("Reading Step 2 Harmony object.")
harmony <- readRDS(HARMONY_RDS)

rpca <- safe_join_layers(rpca, ASSAY_USE)
harmony <- safe_join_layers(harmony, ASSAY_USE)

validate_metadata(
  rpca,
  c(
    RPCA_CLUSTER_COL,
    FINE_ANNOTATION_COL,
    "condition_internal",
    "sample_internal",
    "Score_Low_quality_stress"
  )
)

validate_metadata(
  harmony,
  c(
    HARMONY_CLUSTER_COL,
    FINE_ANNOTATION_COL,
    "condition_internal",
    "sample_internal",
    "Score_Low_quality_stress"
  )
)

# ------------------------------------------------------------------------------
# 11. Assign broad classes and stress/QC state
# ------------------------------------------------------------------------------
message_time("Assigning RPCA broad classes.")
rpca <- assign_broad_class(rpca)
rpca <- assign_stress_qc(rpca)

message_time("Assigning Harmony broad classes.")
harmony <- assign_broad_class(harmony)
harmony <- assign_stress_qc(harmony)

rpca$condition_internal <- factor(
  as.character(rpca$condition_internal),
  levels = condition_levels
)
harmony$condition_internal <- factor(
  as.character(harmony$condition_internal),
  levels = condition_levels
)

rpca$sample_internal <- factor(
  as.character(rpca$sample_internal),
  levels = sample_levels
)
harmony$sample_internal <- factor(
  as.character(harmony$sample_internal),
  levels = sample_levels
)

# ------------------------------------------------------------------------------
# 12. UMAPs
# ------------------------------------------------------------------------------
p_rpca_broad <- plot_broad_umap(
  rpca,
  RPCA_REDUCTION,
  "RPCA"
)

p_harmony_broad <- plot_broad_umap(
  harmony,
  HARMONY_REDUCTION,
  "Harmony"
)

save_pdf(
  "01A_RPCA_UMAP_five_broad_classes.pdf",
  p_rpca_broad,
  11,
  8
)

save_pdf(
  "01B_Harmony_UMAP_five_broad_classes.pdf",
  p_harmony_broad,
  11,
  8
)

p_rpca_split <- plot_broad_umap_split(
  rpca,
  RPCA_REDUCTION,
  "RPCA"
)

p_harmony_split <- plot_broad_umap_split(
  harmony,
  HARMONY_REDUCTION,
  "Harmony"
)

save_pdf(
  "02A_RPCA_UMAP_five_broad_classes_split_by_sample.pdf",
  p_rpca_split,
  16,
  5.5
)

save_pdf(
  "02B_Harmony_UMAP_five_broad_classes_split_by_sample.pdf",
  p_harmony_split,
  16,
  5.5
)

p_rpca_qc <- plot_qc_umap(
  rpca,
  RPCA_REDUCTION,
  "RPCA"
)

p_harmony_qc <- plot_qc_umap(
  harmony,
  HARMONY_REDUCTION,
  "Harmony"
)

save_pdf(
  "03A_RPCA_UMAP_stress_QC_state.pdf",
  p_rpca_qc,
  9,
  7
)

save_pdf(
  "03B_Harmony_UMAP_stress_QC_state.pdf",
  p_harmony_qc,
  9,
  7
)

# ------------------------------------------------------------------------------
# 13. Composition tables and plots
# ------------------------------------------------------------------------------
rpca_composition <- make_composition_table(
  rpca,
  "RPCA"
)

harmony_composition <- make_composition_table(
  harmony,
  "Harmony"
)

composition_all <- bind_rows(
  rpca_composition,
  harmony_composition
)

condition_summary_all <- make_condition_summary(
  composition_all
)

tx_vs_sham_all <- make_tx_vs_sham_summary(
  composition_all
)

write.csv(
  composition_all,
  file.path(
    TABLE_DIR,
    "01_broad_class_fraction_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  condition_summary_all,
  file.path(
    TABLE_DIR,
    "02_broad_class_fraction_condition_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  tx_vs_sham_all,
  file.path(
    TABLE_DIR,
    "03_broad_class_log2FC_Tx_vs_Sham.csv"
  ),
  row.names = FALSE
)

p_composition <- ggplot(
  composition_all,
  aes(
    x = sample_internal,
    y = fraction_of_sample,
    fill = broad_mphi_class_internal
  )
) +
  geom_col(width = 0.78) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(
    values = broad_class_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Five broad macrophage classes by sample",
    x = NULL,
    y = "Fraction of macrophage/monocyte cells",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "right"
  )

save_pdf(
  "04_RPCA_vs_Harmony_broad_class_composition_by_sample.pdf",
  p_composition,
  13,
  6.5
)

p_fraction_line <- ggplot(
  composition_all,
  aes(
    x = broad_mphi_class_internal,
    y = fraction_of_sample,
    group = sample_internal,
    color = sample_internal
  )
) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.85
  ) +
  geom_point(size = 2.2) +
  facet_wrap(~ method, nrow = 1) +
  scale_color_manual(
    values = sample_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1)
  ) +
  labs(
    title = "Broad-class fraction in each biological sample",
    x = NULL,
    y = "Fraction of macrophage/monocyte cells",
    color = "Sample"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_pdf(
  "05_broad_class_fraction_each_sample.pdf",
  p_fraction_line,
  14,
  6.8
)

p_fc <- ggplot(
  tx_vs_sham_all,
  aes(
    x = broad_mphi_class_internal,
    y = log2FC_Tx_vs_Sham,
    fill = broad_mphi_class_internal
  )
) +
  geom_col(width = 0.72) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ method,
    nrow = 1
  ) +
  scale_fill_manual(
    values = broad_class_colors,
    drop = FALSE
  ) +
  labs(
    title = "Exploratory broad-class abundance change: Tx vs Sham",
    x = NULL,
    y = "log2 mean-fraction ratio",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "none"
  )

save_pdf(
  "06_broad_class_abundance_log2FC_Tx_vs_Sham.pdf",
  p_fc,
  13,
  6.5
)

# ------------------------------------------------------------------------------
# 14. DotPlots
# ------------------------------------------------------------------------------
all_marker_genes <- unique(unlist(marker_panels))

rpca_dot_genes <- existing_genes(
  rpca,
  all_marker_genes,
  ASSAY_USE
)

harmony_dot_genes <- existing_genes(
  harmony,
  all_marker_genes,
  ASSAY_USE
)

Idents(rpca) <- rpca$broad_mphi_class_internal
Idents(harmony) <- harmony$broad_mphi_class_internal

p_rpca_dot <- DotPlot(
  rpca,
  features = rpca_dot_genes,
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "RPCA five-class marker profile",
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 7)
  )

p_harmony_dot <- DotPlot(
  harmony,
  features = harmony_dot_genes,
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "Harmony five-class marker profile",
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 7)
  )

save_pdf(
  "07A_RPCA_five_class_marker_DotPlot.pdf",
  p_rpca_dot,
  20,
  6.5
)

save_pdf(
  "07B_Harmony_five_class_marker_DotPlot.pdf",
  p_harmony_dot,
  20,
  6.5
)

# ------------------------------------------------------------------------------
# 15. Stress/QC summaries
# ------------------------------------------------------------------------------
stress_summary <- bind_rows(
  rpca@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = "RPCA",
      broad_mphi_class_internal,
      qc_state_internal,
      name = "n_cells"
    ),
  harmony@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = "Harmony",
      broad_mphi_class_internal,
      qc_state_internal,
      name = "n_cells"
    )
) %>%
  group_by(method, broad_mphi_class_internal) %>%
  mutate(
    class_total = sum(n_cells),
    fraction_within_class = n_cells / class_total
  ) %>%
  ungroup()

write.csv(
  stress_summary,
  file.path(
    TABLE_DIR,
    "04_stress_QC_state_by_broad_class.csv"
  ),
  row.names = FALSE
)

p_stress <- ggplot(
  stress_summary,
  aes(
    x = broad_mphi_class_internal,
    y = fraction_within_class,
    fill = qc_state_internal
  )
) +
  geom_col(width = 0.76) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(
    values = qc_state_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Stress/QC state within each broad macrophage class",
    x = NULL,
    y = "Fraction within broad class",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

save_pdf(
  "08_stress_QC_state_by_broad_class.pdf",
  p_stress,
  14,
  6.5
)

# ------------------------------------------------------------------------------
# 16. Cluster-to-broad mapping and cell-level tables
# ------------------------------------------------------------------------------
rpca_mapping <- make_cluster_to_broad_table(
  rpca,
  RPCA_CLUSTER_COL,
  "RPCA"
)

harmony_mapping <- make_cluster_to_broad_table(
  harmony,
  HARMONY_CLUSTER_COL,
  "Harmony"
)

cluster_mapping_all <- bind_rows(
  rpca_mapping,
  harmony_mapping
)

write.csv(
  cluster_mapping_all,
  file.path(
    TABLE_DIR,
    "05_cluster_fine_broad_QC_mapping.csv"
  ),
  row.names = FALSE
)

rpca_cells <- make_cell_annotation_table(
  rpca,
  RPCA_CLUSTER_COL,
  "RPCA"
)

harmony_cells <- make_cell_annotation_table(
  harmony,
  HARMONY_CLUSTER_COL,
  "Harmony"
)

write.csv(
  rpca_cells,
  file.path(
    TABLE_DIR,
    "06_RPCA_cell_level_five_class_annotation.csv"
  ),
  row.names = FALSE
)

write.csv(
  harmony_cells,
  file.path(
    TABLE_DIR,
    "07_Harmony_cell_level_five_class_annotation.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 17. Broad-class pseudobulk-ready cell count table
# ------------------------------------------------------------------------------
pseudobulk_ready_counts <- bind_rows(
  rpca@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = "RPCA",
      condition_internal,
      sample_internal,
      broad_mphi_class_internal,
      name = "n_cells"
    ),
  harmony@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      method = "Harmony",
      condition_internal,
      sample_internal,
      broad_mphi_class_internal,
      name = "n_cells"
    )
)

write.csv(
  pseudobulk_ready_counts,
  file.path(
    TABLE_DIR,
    "08_pseudobulk_ready_cell_counts_by_broad_class.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 18. Excel workbook
# ------------------------------------------------------------------------------
write_excel_sheets(
  sheets = list(
    Broad_fraction_sample = composition_all,
    Broad_condition_summary = condition_summary_all,
    Tx_vs_Sham = tx_vs_sham_all,
    Stress_QC = stress_summary,
    Cluster_mapping = cluster_mapping_all,
    RPCA_cell_annotation = rpca_cells,
    Harmony_cell_annotation = harmony_cells,
    Pseudobulk_cell_counts = pseudobulk_ready_counts
  ),
  filename = file.path(
    TABLE_DIR,
    "Mouse_Mphi_five_broad_classes_v1.0.0.xlsx"
  )
)

# ------------------------------------------------------------------------------
# 19. Save RDS
# ------------------------------------------------------------------------------
saveRDS(
  rpca,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_RPCA_five_broad_classes_v1.0.0.rds"
  ),
  compress = FALSE
)

saveRDS(
  harmony,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_Harmony_five_broad_classes_v1.0.0.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 20. Run log
# ------------------------------------------------------------------------------
log_lines <- c(
  "analysis_name: Mouse macrophage five broad classes",
  "version: 1.0.0",
  paste0("run_time: ", Sys.time()),
  paste0("RPCA_input: ", RPCA_RDS),
  paste0("Harmony_input: ", HARMONY_RDS),
  "Primary method: RPCA",
  "Sensitivity method: Harmony",
  "Broad classes:",
  paste0("  - ", broad_class_levels),
  "Low_quality_stress handling: independent QC axis",
  "Neutrophil-like handling: provisionally grouped under Monocyte / myeloid-like cells",
  "Dendritic-like handling: provisionally grouped under Monocyte / myeloid-like cells",
  "Final exclusion decisions require marker-level manual review.",
  "sessionInfo:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)

writeLines(
  log_lines,
  file.path(
    LOG_DIR,
    "run_manifest_v1.0.0.txt"
  )
)

message_time("Step 3 completed.")
message_time("Output directory: ", OUTPUT_DIR)
