#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6000)

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
# 5 macrophage subtypes x HSC_Mesenchymal interaction-ready object
#
# Version: v6.0.0
#
# BASE CHECKPOINT
#   mouse-mash-mphi-staining-markers-v5.9.3
#   commit 1381eeae98099a835745313be30c7c6bbc7d1674
#
# PURPOSE
#   Prepare and validate a single whole-cell-derived Seurat object containing:
#     Anti-inflammatory-Mphi
#     Inflammatory-Mphi
#     ECM-associated inflammatory-Mphi
#     Repair/Resolution-Mphi
#     Lipid-associated/TREM2-Mphi
#     HSC_Mesenchymal
#
#   Samples:
#     Sham1, Sham20, Tx17, Tx5
#
# IMPORTANT
#   - No reclustering
#   - No reintegration
#   - No new UMAP
#   - Clean-B FINAL v4.14.5 macrophage labels are transferred to the frozen
#     whole-cell parent v5.1.1 by exact cell-barcode matching.
#   - "Other" macrophages are excluded.
#   - This script performs preparation/QC only.
# ==============================================================================

msg <- function(...) {
  message(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    paste0(...)
  )
}

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("None of the expected files exists:\n", paste(paths, collapse = "\n"))
  }
  hit[[1]]
}

resolve_first <- function(x, candidates, what) {
  hit <- candidates[candidates %in% x]
  if (!length(hit)) {
    stop("Could not resolve ", what, ". Candidates: ", paste(candidates, collapse = ", "))
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
  dplyr::case_when(
    grepl("^Sham", sample, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", sample, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )
}

resolve_hsc_label <- function(labels) {
  preferred <- c("HSC_Mesenchymal", "HSC/Mesenchymal", "HSC-Mesenchymal", "HSC Mesenchymal")
  exact <- preferred[preferred %in% labels]
  if (length(exact)) return(exact[[1]])

  hit <- labels[
    grepl("HSC.*Mesench|Mesench.*HSC|Hepatic.?stellate", labels, ignore.case = TRUE)
  ]

  if (!length(hit)) {
    stop("Could not resolve HSC_Mesenchymal label. Available labels:\n",
         paste(sort(labels), collapse = " | "))
  }
  hit[[1]]
}

save_pdf <- function(p, file, width, height) {
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
}

# Paths
ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

MPHI_RDS <- find_first_existing(c(
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
))

OUT <- file.path(ROOT, "Mouse_MASH_Interaction", "Mphi5_HSC_interaction_ready_v6.0.0")
RDS_OUT <- file.path(OUT, "RDS")
TAB_OUT <- file.path(OUT, "Tables")
FIG_OUT <- file.path(OUT, "Figures")
LOG_OUT <- file.path(OUT, "Logs")

for (d in c(OUT, RDS_OUT, TAB_OUT, FIG_OUT, LOG_OUT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Metadata definitions
WHOLE_LAYER1_COL <- "wholecell_layer1_FINAL_v511"
MPHI_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145_char"

WHOLE_SAMPLE_CANDIDATES <- c("sample_for_annotation", "sample_FIXED2", "sample", "orig.ident")
MPHI_SAMPLE_CANDIDATES <- c("sample_for_annotation", "sample_FIXED2", "sample", "orig.ident")

WHOLE_UMAP <- "umapRPCA"

MPHI5 <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

INTERACTION_LEVELS <- c(MPHI5, "HSC_Mesenchymal")
SAMPLE_LEVELS <- c("Sham1", "Sham20", "Tx17", "Tx5")

# Load
if (!file.exists(WHOLE_RDS)) stop("Whole-cell frozen RDS not found: ", WHOLE_RDS)

msg("Loading whole-cell frozen parent...")
whole <- readRDS(WHOLE_RDS)

msg("Loading Clean-B macrophage FINAL...")
mphi <- readRDS(MPHI_RDS)

if (!WHOLE_LAYER1_COL %in% colnames(whole@meta.data)) {
  stop("Missing whole-cell Layer1 column: ", WHOLE_LAYER1_COL)
}
if (!MPHI_CLASS_COL %in% colnames(mphi@meta.data)) {
  stop("Missing macrophage class column: ", MPHI_CLASS_COL)
}
if (!WHOLE_UMAP %in% Reductions(whole)) {
  stop("Missing frozen whole-cell UMAP: ", WHOLE_UMAP)
}
if (!"RNA" %in% Assays(whole)) stop("RNA assay missing in whole-cell object.")
if (!"RNA" %in% Assays(mphi)) stop("RNA assay missing in macrophage object.")

DefaultAssay(whole) <- "RNA"
DefaultAssay(mphi) <- "RNA"

WHOLE_SAMPLE_COL <- resolve_first(
  colnames(whole@meta.data),
  WHOLE_SAMPLE_CANDIDATES,
  "whole-cell sample column"
)

MPHI_SAMPLE_COL <- resolve_first(
  colnames(mphi@meta.data),
  MPHI_SAMPLE_CANDIDATES,
  "macrophage sample column"
)

msg("Whole-cell sample source: ", WHOLE_SAMPLE_COL)
msg("Mphi sample source: ", MPHI_SAMPLE_COL)

# Canonical sample labels
whole$sample_interaction_v600 <- canonical_sample(
  whole@meta.data[[WHOLE_SAMPLE_COL]]
)
whole$condition_interaction_v600 <- canonical_condition(
  whole$sample_interaction_v600
)

mphi$sample_interaction_v600 <- canonical_sample(
  mphi@meta.data[[MPHI_SAMPLE_COL]]
)
mphi$condition_interaction_v600 <- canonical_condition(
  mphi$sample_interaction_v600
)

msg("Whole-cell canonical sample counts:")
print(table(whole$sample_interaction_v600, useNA = "ifany"))

msg("Mphi canonical sample counts:")
print(table(mphi$sample_interaction_v600, useNA = "ifany"))

# Resolve HSC population
whole_layer1 <- as.character(whole@meta.data[[WHOLE_LAYER1_COL]])
HSC_LABEL_SOURCE <- resolve_hsc_label(unique(whole_layer1))
msg("Resolved HSC label: ", HSC_LABEL_SOURCE)

write.csv(
  tibble(
    role = "HSC receiver population",
    source_column = WHOLE_LAYER1_COL,
    resolved_label = HSC_LABEL_SOURCE,
    final_interaction_label = "HSC_Mesenchymal"
  ),
  file.path(TAB_OUT, "01_HSC_label_resolution_v6.0.0.csv"),
  row.names = FALSE
)

# Exact barcode mapping
mphi_meta <- tibble(
  cell = colnames(mphi),
  mphi_class = as.character(mphi@meta.data[[MPHI_CLASS_COL]]),
  mphi_sample = mphi$sample_interaction_v600,
  mphi_condition = mphi$condition_interaction_v600
)

mphi_meta_shamtx <- mphi_meta %>%
  filter(mphi_sample %in% SAMPLE_LEVELS) %>%
  mutate(exact_match_in_whole = cell %in% colnames(whole))

mapping_summary <- mphi_meta_shamtx %>%
  count(mphi_sample, mphi_class, exact_match_in_whole, name = "n_cells") %>%
  arrange(mphi_sample, mphi_class, desc(exact_match_in_whole))

write.csv(
  mapping_summary,
  file.path(TAB_OUT, "02_barcode_mapping_summary_v6.0.0.csv"),
  row.names = FALSE
)

n_mphi_shamtx <- nrow(mphi_meta_shamtx)
n_exact <- sum(mphi_meta_shamtx$exact_match_in_whole)
mapping_fraction <- n_exact / n_mphi_shamtx

msg(
  "Exact barcode match: ",
  n_exact, " / ", n_mphi_shamtx, " = ",
  sprintf("%.2f%%", 100 * mapping_fraction)
)

write.csv(
  tibble(
    metric = c(
      "Mphi Sham/Tx cells",
      "Exact matches in whole-cell object",
      "Exact-match fraction"
    ),
    value = c(
      n_mphi_shamtx,
      n_exact,
      mapping_fraction
    )
  ),
  file.path(TAB_OUT, "03_barcode_mapping_global_audit_v6.0.0.csv"),
  row.names = FALSE
)

if (!is.finite(mapping_fraction) || mapping_fraction < 0.95) {
  stop(
    "Exact barcode mapping is below 95% (",
    sprintf("%.2f%%", 100 * mapping_fraction),
    "). Do not proceed until barcode naming is resolved."
  )
}

# Transfer labels
whole$interaction_celltype_v600 <- NA_character_

hsc_idx <- which(
  whole_layer1 == HSC_LABEL_SOURCE &
    whole$sample_interaction_v600 %in% SAMPLE_LEVELS
)

whole$interaction_celltype_v600[hsc_idx] <- "HSC_Mesenchymal"

mphi5_map <- mphi_meta_shamtx %>%
  filter(
    exact_match_in_whole,
    mphi_class %in% MPHI5
  ) %>%
  select(cell, mphi_class, mphi_sample, mphi_condition)

if (anyDuplicated(mphi5_map$cell)) {
  stop("Duplicated macrophage cell barcodes detected in transfer map.")
}

transfer_idx <- match(mphi5_map$cell, colnames(whole))
whole$interaction_celltype_v600[transfer_idx] <- mphi5_map$mphi_class

# Sample consistency
whole_sample_at_mphi <- whole$sample_interaction_v600[transfer_idx]

sample_consistency <- tibble(
  cell = mphi5_map$cell,
  mphi_sample = mphi5_map$mphi_sample,
  whole_sample = whole_sample_at_mphi,
  same_sample = mphi5_map$mphi_sample == whole_sample_at_mphi
)

write.csv(
  sample_consistency,
  file.path(TAB_OUT, "04_Mphi_wholecell_sample_consistency_v6.0.0.csv"),
  row.names = FALSE
)

if (any(sample_consistency$same_sample %in% FALSE, na.rm = TRUE)) {
  bad_n <- sum(sample_consistency$same_sample %in% FALSE, na.rm = TRUE)
  stop(
    "Sample mismatch between Clean-B Mphi and whole-cell metadata for ",
    bad_n, " mapped cells."
  )
}

# Build interaction-ready object
keep_cells <- colnames(whole)[
  !is.na(whole$interaction_celltype_v600) &
    whole$sample_interaction_v600 %in% SAMPLE_LEVELS
]

interaction <- subset(whole, cells = keep_cells)

interaction$interaction_celltype_v600 <- factor(
  interaction$interaction_celltype_v600,
  levels = INTERACTION_LEVELS
)

interaction$sample_interaction_v600 <- factor(
  interaction$sample_interaction_v600,
  levels = SAMPLE_LEVELS
)

interaction$condition_interaction_v600 <- factor(
  interaction$condition_interaction_v600,
  levels = c("Sham", "Tx")
)

Idents(interaction) <- "interaction_celltype_v600"

# QC tables
count_table <- interaction@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(
    sample_interaction_v600,
    condition_interaction_v600,
    interaction_celltype_v600,
    name = "n_cells"
  ) %>%
  complete(
    sample_interaction_v600 = factor(SAMPLE_LEVELS, levels = SAMPLE_LEVELS),
    interaction_celltype_v600 = factor(INTERACTION_LEVELS, levels = INTERACTION_LEVELS),
    fill = list(n_cells = 0)
  ) %>%
  group_by(sample_interaction_v600) %>%
  mutate(
    fraction_within_interaction_object = n_cells / sum(n_cells)
  ) %>%
  ungroup()

write.csv(
  count_table,
  file.path(TAB_OUT, "05_sample_by_interaction_celltype_counts_v6.0.0.csv"),
  row.names = FALSE
)

msg("Sample x interaction-celltype counts:")
print(count_table)

min_cell_table <- count_table %>%
  mutate(
    below_10_cells = n_cells < 10,
    below_20_cells = n_cells < 20,
    below_30_cells = n_cells < 30
  )

write.csv(
  min_cell_table,
  file.path(TAB_OUT, "06_CellChat_minimum_cell_audit_v6.0.0.csv"),
  row.names = FALSE
)

if (any(min_cell_table$n_cells < 10)) {
  warning(
    "At least one sample x celltype group has <10 cells. Review ",
    "06_CellChat_minimum_cell_audit_v6.0.0.csv before CellChat."
  )
}

metadata_out <- interaction@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(
    cell,
    sample_interaction_v600,
    condition_interaction_v600,
    interaction_celltype_v600,
    everything()
  )

write.csv(
  metadata_out,
  file.path(TAB_OUT, "07_interaction_ready_cell_metadata_v6.0.0.csv"),
  row.names = FALSE
)

# RNA layer audit
rna_layers <- Layers(interaction[["RNA"]])

write.csv(
  tibble(assay = "RNA", layer = rna_layers),
  file.path(TAB_OUT, "08_RNA_layer_audit_v6.0.0.csv"),
  row.names = FALSE
)

msg("RNA layers: ", paste(rna_layers, collapse = ", "))

# UMAP QC
PLOT_COLORS <- c(
  "Anti-inflammatory-Mphi" = "#00C8FF",
  "Inflammatory-Mphi" = "#FF3131",
  "ECM-associated inflammatory-Mphi" = "#FF8C00",
  "Repair/Resolution-Mphi" = "#A020F0",
  "Lipid-associated/TREM2-Mphi" = "#0066FF",
  "HSC_Mesenchymal" = "#FF3FA4"
)

p_all <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by = "interaction_celltype_v600",
  cols = PLOT_COLORS,
  pt.size = 0.35,
  raster = FALSE
) +
  ggtitle("Interaction-ready object | 5 Mphi subtypes + HSC_Mesenchymal") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_pdf(
  p_all,
  file.path(FIG_OUT, "01_interaction_ready_UMAP_all_samples_v6.0.0.pdf"),
  9, 7
)

p_sample <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by = "interaction_celltype_v600",
  split.by = "sample_interaction_v600",
  cols = PLOT_COLORS,
  pt.size = 0.30,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(title = "Interaction-ready object | biological samples")

save_pdf(
  p_sample,
  file.path(FIG_OUT, "02_interaction_ready_UMAP_by_sample_v6.0.0.pdf"),
  13, 10
)

p_condition <- DimPlot(
  interaction,
  reduction = WHOLE_UMAP,
  group.by = "interaction_celltype_v600",
  split.by = "condition_interaction_v600",
  cols = PLOT_COLORS,
  pt.size = 0.30,
  raster = FALSE,
  ncol = 2
) +
  plot_annotation(title = "Interaction-ready object | Sham vs Tx")

save_pdf(
  p_condition,
  file.path(FIG_OUT, "03_interaction_ready_UMAP_Sham_vs_Tx_v6.0.0.pdf"),
  13, 6
)

# Cell count heatmap
p_count_heat <- ggplot(
  count_table,
  aes(
    x = interaction_celltype_v600,
    y = sample_interaction_v600,
    fill = log10(n_cells + 1)
  )
) +
  geom_tile(linewidth = 0.4) +
  geom_text(aes(label = n_cells), size = 3) +
  scale_fill_gradientn(
    colours = c("#0033FF", "#FFFFFF", "#FF1A1A")
  ) +
  labs(
    title = "Cells available for per-sample interaction analysis",
    x = NULL,
    y = NULL,
    fill = "log10(n+1)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_pdf(
  p_count_heat,
  file.path(FIG_OUT, "04_sample_by_celltype_cellcount_heatmap_v6.0.0.pdf"),
  11, 5
)

# Save interaction-ready RDS
INTERACTION_RDS <- file.path(
  RDS_OUT,
  "Mouse_MASH_Mphi5_HSC_interaction_ready_v6.0.0.rds"
)

saveRDS(
  interaction,
  INTERACTION_RDS,
  compress = FALSE
)

msg("Saved interaction-ready RDS: ", INTERACTION_RDS)

# Final audit
final_audit <- tibble(
  metric = c(
    "wholecell_cells",
    "CleanB_Mphi_cells",
    "CleanB_Mphi_ShamTx_cells",
    "exact_barcode_matches",
    "exact_barcode_match_fraction",
    "interaction_ready_cells",
    "HSC_source_label",
    "wholecell_sample_column",
    "Mphi_sample_column",
    "wholecell_UMAP"
  ),
  value = c(
    as.character(ncol(whole)),
    as.character(ncol(mphi)),
    as.character(n_mphi_shamtx),
    as.character(n_exact),
    as.character(mapping_fraction),
    as.character(ncol(interaction)),
    HSC_LABEL_SOURCE,
    WHOLE_SAMPLE_COL,
    MPHI_SAMPLE_COL,
    WHOLE_UMAP
  )
)

write.csv(
  final_audit,
  file.path(LOG_OUT, "final_audit_v6.0.0.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(LOG_OUT, "sessionInfo_v6.0.0.txt")
)

msg("DONE.")
msg("Output directory: ", OUT)
