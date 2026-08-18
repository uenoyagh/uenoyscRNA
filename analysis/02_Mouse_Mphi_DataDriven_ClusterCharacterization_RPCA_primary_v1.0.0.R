#!/usr/bin/env Rscript

# ==============================================================================
# uenoyscRNA
# Mouse macrophage data-driven analysis Step 2
# RPCA-centered cluster characterization with Harmony sensitivity analysis
# Version 1.0.0
# ==============================================================================
#
# INPUT:
#   Step 1 RPCA and Harmony RDS files
#
# PRIMARY ANALYSIS:
#   RPCA resolution 1.2 clusters
#
# SENSITIVITY ANALYSIS:
#   Harmony resolution 1.2 clusters
#
# OUTPUTS:
#   - Cluster markers
#   - Canonical marker DotPlots
#   - QC/contamination DotPlots
#   - Sample bias tables
#   - Sham vs Tx cell fraction summaries
#   - Preliminary annotation tables
#   - RPCA vs Harmony annotation concordance support
#
# NOTE:
#   This script does not perform final biological annotation automatically.
#   It produces evidence tables and a conservative provisional annotation.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
INPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  "DataDriven_ShamTx_RPCA_vs_Harmony_v1.0.0",
  "RDS"
)

RPCA_RDS <- file.path(
  INPUT_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_RPCA_v1.0.0.rds"
)

HARMONY_RDS <- file.path(
  INPUT_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_Harmony_v1.0.0.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  "DataDriven_ShamTx_RPCA_vs_Harmony_v1.0.0",
  "Step2_ClusterCharacterization_v1.0.0"
)

FIGURE_DIR <- file.path(OUTPUT_DIR, "Figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "Tables")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

for (d in c(OUTPUT_DIR, FIGURE_DIR, TABLE_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RPCA_CLUSTER_COL <- "mphi_rpca_res_1_2"
HARMONY_CLUSTER_COL <- "mphi_harmony_res_1_2"

ASSAY_USE <- "RNA"

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
# 2. Helpers
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

write_excel_sheets <- function(sheets, filename) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
    freezePane(wb, nm, firstRow = TRUE)
    setColWidths(wb, nm, cols = 1:ncol(sheets[[nm]]), widths = "auto")
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}

normalized_entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) <= 1) return(0)
  p <- x / sum(x)
  -sum(p * log(p)) / log(length(p))
}

dominant_sample_metrics <- function(meta, cluster_col) {
  meta %>%
    rownames_to_column("cell") %>%
    count(cluster = .data[[cluster_col]], sample_internal, name = "n_cells") %>%
    group_by(cluster) %>%
    mutate(
      cluster_total = sum(n_cells),
      fraction = n_cells / cluster_total
    ) %>%
    summarise(
      n_cells = first(cluster_total),
      n_samples_present = sum(n_cells > 0),
      dominant_sample = sample_internal[which.max(fraction)],
      dominant_sample_fraction = max(fraction),
      sample_entropy = normalized_entropy(n_cells),
      .groups = "drop"
    )
}

sample_fraction_table <- function(meta, cluster_col, method_label) {
  meta %>%
    rownames_to_column("cell") %>%
    count(
      method = method_label,
      condition_internal,
      sample_internal,
      cluster = .data[[cluster_col]],
      name = "n_cells"
    ) %>%
    group_by(method, condition_internal, sample_internal) %>%
    mutate(
      sample_total = sum(n_cells),
      fraction_of_sample = n_cells / sample_total
    ) %>%
    ungroup()
}

condition_summary <- function(sample_fraction_df) {
  sample_fraction_df %>%
    group_by(method, condition_internal, cluster) %>%
    summarise(
      mean_fraction = mean(fraction_of_sample),
      sd_fraction = sd(fraction_of_sample),
      n_samples = n_distinct(sample_internal),
      .groups = "drop"
    )
}

sham_tx_fold_change <- function(sample_fraction_df) {
  wide <- sample_fraction_df %>%
    group_by(method, condition_internal, cluster) %>%
    summarise(
      mean_fraction = mean(fraction_of_sample),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = condition_internal,
      values_from = mean_fraction,
      values_fill = 0
    )

  if (!all(c("Sham", "Tx") %in% colnames(wide))) {
    stop("Sham and Tx columns were not found in composition table.")
  }

  wide %>%
    mutate(
      log2FC_Tx_vs_Sham = log2((Tx + 1e-4) / (Sham + 1e-4)),
      direction = case_when(
        log2FC_Tx_vs_Sham > 0.25 ~ "Higher in Tx",
        log2FC_Tx_vs_Sham < -0.25 ~ "Lower in Tx",
        TRUE ~ "Similar"
      )
    )
}

# ------------------------------------------------------------------------------
# 3. Marker sets
# ------------------------------------------------------------------------------
marker_sets <- list(
  Resident_Kupffer = c(
    "Clec4f", "Timd4", "Vsig4", "Marco", "Cd5l", "C1qa", "C1qb", "C1qc",
    "Fcrls", "Folr2", "Slc40a1", "Hpgd"
  ),
  Monocyte_like = c(
    "Ly6c2", "Ccr2", "Sell", "Plac8", "S100a8", "S100a9", "Lyz2",
    "Fcgr3", "Ms4a7", "Ctss"
  ),
  Inflammatory_M1_like = c(
    "Il1b", "Tnf", "Nfkbia", "Ptgs2", "Cxcl2", "Cxcl3", "Ccl2",
    "Ccl3", "Ccl4", "Cd80", "Cd86", "Nos2", "Stat1", "Irf1"
  ),
  Pro_resolution_M2_like = c(
    "Il10", "Mertk", "Axl", "Gas6", "Maf", "Mafb", "Cd163", "Mrc1",
    "Retnla", "Chil3", "Arg1", "Clec10a", "Folr2"
  ),
  Efferocytosis_high = c(
    "Mertk", "Axl", "Tyro3", "Gas6", "Mfge8", "Cd36", "Itgav", "Itgb5",
    "Lrp1", "Abca1", "Abcg1", "Nr1h3", "Pparg"
  ),
  SPP1_TREM2_MASH = c(
    "Spp1", "Trem2", "Gpnmb", "Cd9", "Lgals3", "Fabp5", "Lpl",
    "Ctsb", "Ctsd", "Ctsl", "Apoe", "Tyrobp"
  ),
  IL10_response = c(
    "Il10ra", "Il10rb", "Stat3", "Socs3", "Bcl3", "Dusp1", "Dusp5",
    "Tnfaip3", "Klf4", "Maf"
  ),
  Fibrosis_associated = c(
    "Spp1", "Tgfb1", "Pdgfb", "Mmp12", "Mmp14", "Timp1", "Lgals3",
    "Fn1", "Thbs1", "Ccl2"
  ),
  Cycling = c(
    "Mki67", "Top2a", "Tubb5", "Stmn1", "Cenpf", "Cdk1", "Pcna"
  ),
  Dendritic_like = c(
    "Flt3", "Zbtb46", "Clec10a", "Cd209a", "Itgax", "H2-Ab1",
    "H2-Aa", "H2-Eb1", "Ccr7"
  ),
  Neutrophil_like = c(
    "S100a8", "S100a9", "Ly6g", "Csf3r", "Mpo", "Elane", "Camp",
    "Retnlg", "Ngp"
  ),
  Low_quality_stress = c(
    "Fos", "Jun", "Junb", "Ddit3", "Hspa1a", "Hspa1b", "Atf3"
  )
)

# Conservative rule-based provisional annotation.
# This is only a first-pass label based on average module score rank.
annotation_priority <- c(
  "Neutrophil_like",
  "Cycling",
  "Dendritic_like",
  "SPP1_TREM2_MASH",
  "Resident_Kupffer",
  "Monocyte_like",
  "Inflammatory_M1_like",
  "Pro_resolution_M2_like",
  "Efferocytosis_high"
)

annotation_label_map <- c(
  Resident_Kupffer = "Resident Kupffer-like M.",
  Monocyte_like = "Monocyte-like M.",
  Inflammatory_M1_like = "Inflammatory M1-like M.",
  Pro_resolution_M2_like = "Pro-resolution M2-like M.",
  Efferocytosis_high = "Efferocytosis-high M.",
  SPP1_TREM2_MASH = "SPP1/TREM2 MASH-associated M.",
  Cycling = "Cycling M.",
  Dendritic_like = "Dendritic-like contamination",
  Neutrophil_like = "Neutrophil-like contamination"
)

# ------------------------------------------------------------------------------
# 4. Read objects
# ------------------------------------------------------------------------------
for (f in c(RPCA_RDS, HARMONY_RDS)) {
  if (!file.exists(f)) stop("Input RDS not found: ", f)
}

message_time("Reading RPCA object.")
rpca <- readRDS(RPCA_RDS)
message_time("Reading Harmony object.")
harmony <- readRDS(HARMONY_RDS)

rpca <- safe_join_layers(rpca, ASSAY_USE)
harmony <- safe_join_layers(harmony, ASSAY_USE)

for (x in list(rpca = rpca, harmony = harmony)) {
  if (!all(c("condition_internal", "sample_internal") %in% colnames(x@meta.data))) {
    stop("Required internal metadata columns are missing.")
  }
}

if (!(RPCA_CLUSTER_COL %in% colnames(rpca@meta.data))) {
  stop("RPCA cluster column missing: ", RPCA_CLUSTER_COL)
}
if (!(HARMONY_CLUSTER_COL %in% colnames(harmony@meta.data))) {
  stop("Harmony cluster column missing: ", HARMONY_CLUSTER_COL)
}

# ------------------------------------------------------------------------------
# 5. Module scores
# ------------------------------------------------------------------------------
add_scores <- function(object, assay = "RNA") {
  DefaultAssay(object) <- assay

  for (nm in names(marker_sets)) {
    genes <- existing_genes(object, marker_sets[[nm]], assay)
    if (length(genes) >= 3) {
      object <- AddModuleScore(
        object,
        features = list(genes),
        name = paste0("Score_", nm, "_"),
        assay = assay,
        seed = 20260731,
        search = FALSE
      )
      created <- paste0("Score_", nm, "_1")
      names(object@meta.data)[names(object@meta.data) == created] <-
        paste0("Score_", nm)
    } else {
      warning("Insufficient genes for module score: ", nm)
      object[[paste0("Score_", nm)]] <- NA_real_
    }
  }

  object
}

message_time("Calculating RPCA module scores.")
rpca <- add_scores(rpca)

message_time("Calculating Harmony module scores.")
harmony <- add_scores(harmony)

# ------------------------------------------------------------------------------
# 6. Cluster markers
# ------------------------------------------------------------------------------
find_cluster_markers <- function(object, cluster_col, method_label) {
  Idents(object) <- object@meta.data[[cluster_col]]
  DefaultAssay(object) <- ASSAY_USE

  message_time("Finding ", method_label, " cluster markers.")
  markers <- FindAllMarkers(
    object,
    assay = ASSAY_USE,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.10,
    logfc.threshold = 0.20,
    return.thresh = 0.05,
    verbose = FALSE
  )

  markers %>%
    mutate(method = method_label) %>%
    relocate(method)
}

rpca_markers <- find_cluster_markers(rpca, RPCA_CLUSTER_COL, "RPCA")
harmony_markers <- find_cluster_markers(harmony, HARMONY_CLUSTER_COL, "Harmony")

write.csv(
  rpca_markers,
  file.path(TABLE_DIR, "01_RPCA_all_positive_cluster_markers.csv"),
  row.names = FALSE
)
write.csv(
  harmony_markers,
  file.path(TABLE_DIR, "02_Harmony_all_positive_cluster_markers.csv"),
  row.names = FALSE
)

top20_markers <- function(markers) {
  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC"),
    colnames(markers)
  )[1]

  if (is.na(fc_col)) stop("No average log fold-change column found.")

  markers %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = 20, with_ties = FALSE) %>%
    ungroup()
}

rpca_top20 <- top20_markers(rpca_markers)
harmony_top20 <- top20_markers(harmony_markers)

# ------------------------------------------------------------------------------
# 7. DotPlots
# ------------------------------------------------------------------------------
all_canonical_genes <- unique(unlist(marker_sets))
rpca_dot_genes <- existing_genes(rpca, all_canonical_genes, ASSAY_USE)
harmony_dot_genes <- existing_genes(harmony, all_canonical_genes, ASSAY_USE)

Idents(rpca) <- rpca@meta.data[[RPCA_CLUSTER_COL]]
Idents(harmony) <- harmony@meta.data[[HARMONY_CLUSTER_COL]]

p_rpca_dot <- DotPlot(
  rpca,
  features = rpca_dot_genes,
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "RPCA cluster canonical-marker profile",
    x = NULL,
    y = "RPCA cluster"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 7)
  )

p_harmony_dot <- DotPlot(
  harmony,
  features = harmony_dot_genes,
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "Harmony cluster canonical-marker profile",
    x = NULL,
    y = "Harmony cluster"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 7)
  )

save_pdf("01_RPCA_canonical_marker_DotPlot.pdf", p_rpca_dot, 22, 8.5)
save_pdf("02_Harmony_canonical_marker_DotPlot.pdf", p_harmony_dot, 22, 8.5)

qc_genes <- unique(c(
  marker_sets$Neutrophil_like,
  marker_sets$Dendritic_like,
  marker_sets$Cycling,
  marker_sets$Low_quality_stress
))

p_rpca_qc <- DotPlot(
  rpca,
  features = existing_genes(rpca, qc_genes),
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "RPCA contamination and QC marker profile",
    x = NULL,
    y = "RPCA cluster"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_harmony_qc <- DotPlot(
  harmony,
  features = existing_genes(harmony, qc_genes),
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "Harmony contamination and QC marker profile",
    x = NULL,
    y = "Harmony cluster"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_pdf("03_RPCA_QC_contamination_DotPlot.pdf", p_rpca_qc, 15, 7.5)
save_pdf("04_Harmony_QC_contamination_DotPlot.pdf", p_harmony_qc, 15, 7.5)

# ------------------------------------------------------------------------------
# 8. Cluster-level module-score summaries
# ------------------------------------------------------------------------------
score_cols <- grep("^Score_", colnames(rpca@meta.data), value = TRUE)

cluster_score_summary <- function(object, cluster_col, method_label) {
  object@meta.data %>%
    rownames_to_column("cell") %>%
    group_by(cluster = .data[[cluster_col]]) %>%
    summarise(
      across(
        all_of(score_cols),
        list(mean = ~ mean(.x, na.rm = TRUE),
             median = ~ median(.x, na.rm = TRUE))
      ),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    mutate(method = method_label) %>%
    relocate(method)
}

rpca_score_summary <- cluster_score_summary(
  rpca, RPCA_CLUSTER_COL, "RPCA"
)
harmony_score_summary <- cluster_score_summary(
  harmony, HARMONY_CLUSTER_COL, "Harmony"
)

write.csv(
  rpca_score_summary,
  file.path(TABLE_DIR, "03_RPCA_cluster_module_score_summary.csv"),
  row.names = FALSE
)
write.csv(
  harmony_score_summary,
  file.path(TABLE_DIR, "04_Harmony_cluster_module_score_summary.csv"),
  row.names = FALSE
)

# Heatmap-like tile plot for mean scores
score_mean_cols <- grep("_mean$", colnames(rpca_score_summary), value = TRUE)

score_long_plot <- function(summary_df, method_label) {
  summary_df %>%
    select(method, cluster, all_of(score_mean_cols)) %>%
    pivot_longer(
      cols = all_of(score_mean_cols),
      names_to = "module",
      values_to = "mean_score"
    ) %>%
    mutate(
      module = sub("^Score_", "", module),
      module = sub("_mean$", "", module),
      z_score = ave(
        mean_score,
        module,
        FUN = function(x) as.numeric(scale(x))
      )
    ) %>%
    ggplot(aes(x = module, y = factor(cluster), fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title = paste0(method_label, " cluster module-score profile"),
      x = NULL,
      y = paste0(method_label, " cluster"),
      fill = "Cluster-z"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 55, hjust = 1)
    )
}

p_rpca_score <- score_long_plot(rpca_score_summary, "RPCA")
p_harmony_score <- score_long_plot(harmony_score_summary, "Harmony")

save_pdf("05_RPCA_cluster_module_score_heatmap.pdf", p_rpca_score, 12, 8)
save_pdf("06_Harmony_cluster_module_score_heatmap.pdf", p_harmony_score, 12, 8)

# ------------------------------------------------------------------------------
# 9. Sample bias and composition
# ------------------------------------------------------------------------------
rpca_bias <- dominant_sample_metrics(rpca@meta.data, RPCA_CLUSTER_COL) %>%
  mutate(method = "RPCA")
harmony_bias <- dominant_sample_metrics(harmony@meta.data, HARMONY_CLUSTER_COL) %>%
  mutate(method = "Harmony")

bias_all <- bind_rows(rpca_bias, harmony_bias)
write.csv(
  bias_all,
  file.path(TABLE_DIR, "05_cluster_sample_bias_metrics.csv"),
  row.names = FALSE
)

rpca_fraction <- sample_fraction_table(
  rpca@meta.data,
  RPCA_CLUSTER_COL,
  "RPCA"
)
harmony_fraction <- sample_fraction_table(
  harmony@meta.data,
  HARMONY_CLUSTER_COL,
  "Harmony"
)

fraction_all <- bind_rows(rpca_fraction, harmony_fraction)

condition_all <- condition_summary(fraction_all)
fc_all <- sham_tx_fold_change(fraction_all)

write.csv(
  fraction_all,
  file.path(TABLE_DIR, "06_cluster_fraction_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  condition_all,
  file.path(TABLE_DIR, "07_cluster_fraction_condition_summary.csv"),
  row.names = FALSE
)
write.csv(
  fc_all,
  file.path(TABLE_DIR, "08_cluster_fraction_log2FC_Tx_vs_Sham.csv"),
  row.names = FALSE
)

p_fraction <- ggplot(
  fraction_all,
  aes(
    x = factor(cluster),
    y = fraction_of_sample,
    group = sample_internal,
    color = sample_internal
  )
) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 1.8) +
  facet_grid(method ~ condition_internal, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Cluster fraction by biological sample",
    x = "Cluster",
    y = "Fraction within macrophage/monocyte cells",
    color = "Sample"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_pdf("07_cluster_fraction_by_sample_RPCA_vs_Harmony.pdf", p_fraction, 14, 8)

p_fc <- ggplot(
  fc_all,
  aes(
    x = factor(cluster),
    y = log2FC_Tx_vs_Sham,
    fill = direction
  )
) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ method, nrow = 1, scales = "free_x") +
  labs(
    title = "Exploratory cluster abundance change: Tx vs Sham",
    x = "Cluster",
    y = "log2 mean-fraction ratio"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

save_pdf("08_cluster_abundance_log2FC_Tx_vs_Sham.pdf", p_fc, 12, 6)

# ------------------------------------------------------------------------------
# 10. Provisional annotation
# ------------------------------------------------------------------------------
provisional_annotation <- function(score_summary, method_label) {
  mean_cols <- grep("_mean$", colnames(score_summary), value = TRUE)

  score_long <- score_summary %>%
    select(cluster, n_cells, all_of(mean_cols)) %>%
    pivot_longer(
      cols = all_of(mean_cols),
      names_to = "module",
      values_to = "mean_score"
    ) %>%
    mutate(
      module = sub("^Score_", "", module),
      module = sub("_mean$", "", module)
    ) %>%
    group_by(module) %>%
    mutate(z_score = as.numeric(scale(mean_score))) %>%
    ungroup()

  top_module <- score_long %>%
    filter(module %in% annotation_priority) %>%
    group_by(cluster) %>%
    arrange(desc(z_score), match(module, annotation_priority)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      provisional_annotation = unname(annotation_label_map[module]),
      provisional_annotation = if_else(
        is.na(provisional_annotation) | z_score < 0.5,
        "Other macrophage",
        provisional_annotation
      ),
      method = method_label
    ) %>%
    select(
      method,
      cluster,
      n_cells,
      top_module = module,
      top_module_z = z_score,
      provisional_annotation
    )

  top_module
}

rpca_annotation <- provisional_annotation(rpca_score_summary, "RPCA")
harmony_annotation <- provisional_annotation(harmony_score_summary, "Harmony")

# Add sample-bias warning
rpca_annotation <- rpca_annotation %>%
  left_join(
    rpca_bias %>% select(cluster, dominant_sample, dominant_sample_fraction, sample_entropy),
    by = "cluster"
  ) %>%
  mutate(
    sample_bias_warning = dominant_sample_fraction >= 0.70
  )

harmony_annotation <- harmony_annotation %>%
  left_join(
    harmony_bias %>% select(cluster, dominant_sample, dominant_sample_fraction, sample_entropy),
    by = "cluster"
  ) %>%
  mutate(
    sample_bias_warning = dominant_sample_fraction >= 0.70
  )

write.csv(
  rpca_annotation,
  file.path(TABLE_DIR, "09_RPCA_provisional_annotation.csv"),
  row.names = FALSE
)
write.csv(
  harmony_annotation,
  file.path(TABLE_DIR, "10_Harmony_provisional_annotation.csv"),
  row.names = FALSE
)

# Add provisional labels to objects
rpca_map <- setNames(
  rpca_annotation$provisional_annotation,
  rpca_annotation$cluster
)
harmony_map <- setNames(
  harmony_annotation$provisional_annotation,
  harmony_annotation$cluster
)

rpca$layer2_internal_provisional <- unname(
  rpca_map[as.character(rpca@meta.data[[RPCA_CLUSTER_COL]])]
)
harmony$layer2_internal_provisional <- unname(
  harmony_map[as.character(harmony@meta.data[[HARMONY_CLUSTER_COL]])]
)

# ------------------------------------------------------------------------------
# 11. Annotated UMAPs
# ------------------------------------------------------------------------------
p_rpca_annotation <- DimPlot(
  rpca,
  reduction = "mphi.umap.rpca",
  group.by = "layer2_internal_provisional",
  pt.size = 0.55,
  label = TRUE,
  repel = TRUE
) +
  labs(title = "RPCA provisional macrophage annotation") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

p_harmony_annotation <- DimPlot(
  harmony,
  reduction = "mphi.umap.harmony",
  group.by = "layer2_internal_provisional",
  pt.size = 0.55,
  label = TRUE,
  repel = TRUE
) +
  labs(title = "Harmony provisional macrophage annotation") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))

save_pdf("09_RPCA_provisional_annotation_UMAP.pdf", p_rpca_annotation, 11, 8)
save_pdf("10_Harmony_provisional_annotation_UMAP.pdf", p_harmony_annotation, 11, 8)

# ------------------------------------------------------------------------------
# 12. Excel workbook
# ------------------------------------------------------------------------------
write_excel_sheets(
  sheets = list(
    RPCA_annotation = rpca_annotation,
    Harmony_annotation = harmony_annotation,
    RPCA_top20_markers = rpca_top20,
    Harmony_top20_markers = harmony_top20,
    RPCA_module_scores = rpca_score_summary,
    Harmony_module_scores = harmony_score_summary,
    Sample_bias = bias_all,
    Fraction_by_sample = fraction_all,
    Tx_vs_Sham_fraction = fc_all
  ),
  filename = file.path(
    TABLE_DIR,
    "Mouse_Mphi_RPCA_Harmony_cluster_characterization_v1.0.0.xlsx"
  )
)

# ------------------------------------------------------------------------------
# 13. Save annotated objects
# ------------------------------------------------------------------------------
saveRDS(
  rpca,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_DataDriven_RPCA_cluster_characterized_v1.0.0.rds"
  ),
  compress = FALSE
)

saveRDS(
  harmony,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_DataDriven_Harmony_cluster_characterized_v1.0.0.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 14. Run log
# ------------------------------------------------------------------------------
log_lines <- c(
  "analysis_name: Mouse macrophage Step 2 cluster characterization",
  "version: 1.0.0",
  paste0("run_time: ", Sys.time()),
  paste0("RPCA_RDS: ", RPCA_RDS),
  paste0("Harmony_RDS: ", HARMONY_RDS),
  paste0("RPCA_cluster_column: ", RPCA_CLUSTER_COL),
  paste0("Harmony_cluster_column: ", HARMONY_CLUSTER_COL),
  paste0("RPCA_cells: ", ncol(rpca)),
  paste0("Harmony_cells: ", ncol(harmony)),
  "Primary method: RPCA",
  "Sensitivity method: Harmony",
  "Annotation status: provisional; requires manual marker review",
  "sessionInfo:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)

writeLines(
  log_lines,
  file.path(LOG_DIR, "run_manifest_v1.0.0.txt")
)

message_time("Step 2 completed.")
message_time("Output directory: ", OUTPUT_DIR)
