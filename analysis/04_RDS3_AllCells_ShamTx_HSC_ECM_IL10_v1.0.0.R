#!/usr/bin/env Rscript

# ==============================================================================
# uenoyscRNA
# Mouse whole-cell RDS3 analysis
# Sham vs Tx: HSC abundance, ECM production, IL10 downstream signaling,
# and whole-cell IL10/IL10 receptor UMAPs
# Version 1.0.0
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ------------------------------------------------------------------------------
# 0. Paths and settings
# ------------------------------------------------------------------------------
INPUT_RDS <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS",
  "RDS3_ShamTx_HSC_IL10_ECM_v1.0.0"
)

FIGURE_DIR <- file.path(OUTPUT_DIR, "Figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "Tables")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

for (d in c(OUTPUT_DIR, FIGURE_DIR, TABLE_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

CELLTYPE_COL <- "celltype_for_R8plot_FIXED2"
CONDITION_COL <- "condition"
SAMPLE_COL <- "sample"
ASSAY_USE <- "RNA"

TARGET_CONDITIONS <- c("Sham", "Tx")
SAMPLE_ORDER <- c("Sham1", "Sham20", "Tx17", "Tx5")

# HSC labels are auto-detected from these aliases.
AHSC_ALIASES <- c(
  "aHSC", "Activated_HSC", "Activated HSC", "Activated-HSC",
  "HSC_activated", "Activated hepatic stellate cell",
  "Activated hepatic stellate cells"
)

QHSC_ALIASES <- c(
  "qHSC", "Quiescent_HSC", "Quiescent HSC", "Quiescent-HSC",
  "HSC_quiescent", "Quiescent hepatic stellate cell",
  "Quiescent hepatic stellate cells"
)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------
required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "scales", "openxlsx"
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

canonical_token <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", trimws(as.character(x))))
}

match_aliases <- function(values, aliases) {
  canonical_token(values) %in% canonical_token(aliases)
}

safe_join_layers <- function(object, assay = "RNA") {
  DefaultAssay(object) <- assay
  lyr <- Layers(object[[assay]])
  if (any(grepl("^counts\\.", lyr)) || any(grepl("^data\\.", lyr))) {
    object <- JoinLayers(object, assay = assay)
  }
  object
}

existing_genes <- function(object, genes, assay = "RNA") {
  intersect(genes, rownames(object[[assay]]))
}

detect_umap <- function(object) {
  reductions <- Reductions(object)
  preferred <- c("umap", "wnn.umap", "integrated.umap", "rpca.umap")
  hit <- preferred[preferred %in% reductions]
  if (length(hit) > 0L) return(hit[1])

  umap_like <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  if (length(umap_like) == 0L) {
    stop("No UMAP reduction found. Available reductions: ",
         paste(reductions, collapse = ", "))
  }
  umap_like[1]
}

write_excel_sheets <- function(sheets, filename) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
    freezePane(wb, nm, firstRow = TRUE)
    if (ncol(sheets[[nm]]) > 0) {
      setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
    }
  }
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Marker/module definitions
# ------------------------------------------------------------------------------
ecm_production_genes <- c(
  "Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2",
  "Col6a1", "Col6a2", "Col6a3", "Fn1", "Postn",
  "Sparc", "Thbs1", "Timp1", "Timp2"
)

ecm_crosslinking_genes <- c(
  "Lox", "Loxl1", "Loxl2", "Plod2", "P4ha1", "P4ha2"
)

ahsc_genes <- c(
  "Acta2", "Tagln", "Pdgfrb", "Pdgfra", "Col1a1", "Col1a2",
  "Col3a1", "Timp1", "Lox", "Loxl1"
)

qhsc_genes <- c(
  "Lrat", "Reln", "Rgs5", "Dcn", "Igfbp3", "Cygb",
  "Adirf", "Pparg", "Rbp1", "Gpx3"
)

il10_receptor_genes <- c("Il10ra", "Il10rb")

# Conservative IL10/STAT3 downstream-response module.
# Il10 itself and receptor genes are excluded from the downstream score.
il10_downstream_genes <- c(
  "Stat3", "Socs3", "Bcl3", "Tnfaip3", "Dusp1", "Dusp5",
  "Klf4", "Maf", "Mafb", "Nfkbia", "Zfp36", "Btg2"
)

il10_antifibrotic_response_genes <- c(
  "Socs3", "Tnfaip3", "Dusp1", "Klf4", "Maf",
  "Pparg", "Nr1h3", "Abca1"
)

feature_panel <- unique(c(
  "Il10", "Il10ra", "Il10rb",
  ahsc_genes, qhsc_genes,
  ecm_production_genes, ecm_crosslinking_genes,
  il10_downstream_genes
))

# ------------------------------------------------------------------------------
# 4. Load and subset whole-cell RDS3
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_RDS)) {
  stop("Input RDS not found: ", INPUT_RDS)
}

message_time("Reading whole-cell RDS3: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)
obj <- safe_join_layers(obj, ASSAY_USE)

required_meta <- c(CELLTYPE_COL, CONDITION_COL, SAMPLE_COL)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) {
  stop("Missing metadata column(s): ", paste(missing_meta, collapse = ", "))
}

obj$condition_internal <- as.character(obj@meta.data[[CONDITION_COL]])
obj$sample_internal <- as.character(obj@meta.data[[SAMPLE_COL]])
obj$celltype_internal <- as.character(obj@meta.data[[CELLTYPE_COL]])

keep <- obj$condition_internal %in% TARGET_CONDITIONS
obj_st <- subset(obj, cells = colnames(obj)[keep])
rm(obj)
invisible(gc())

obj_st$condition_internal <- factor(
  as.character(obj_st$condition_internal),
  levels = TARGET_CONDITIONS
)
obj_st$sample_internal <- factor(
  as.character(obj_st$sample_internal),
  levels = SAMPLE_ORDER
)

umap_reduction <- detect_umap(obj_st)
message_time("Using UMAP reduction: ", umap_reduction)

# Save available cell type labels for audit.
celltype_values <- obj_st@meta.data %>%
  count(celltype_internal, name = "n_cells", sort = TRUE)

write.csv(
  celltype_values,
  file.path(TABLE_DIR, "00_available_celltype_labels_ShamTx.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Detect aHSC/qHSC labels and create HSC subset
# ------------------------------------------------------------------------------
is_ahsc <- match_aliases(obj_st$celltype_internal, AHSC_ALIASES)
is_qhsc <- match_aliases(obj_st$celltype_internal, QHSC_ALIASES)

if (!any(is_ahsc) || !any(is_qhsc)) {
  stop(
    "Could not detect both aHSC and qHSC labels.\n",
    "Review: ",
    file.path(TABLE_DIR, "00_available_celltype_labels_ShamTx.csv")
  )
}

obj_st$hsc_state_internal <- NA_character_
obj_st$hsc_state_internal[is_ahsc] <- "aHSC"
obj_st$hsc_state_internal[is_qhsc] <- "qHSC"

hsc_cells <- colnames(obj_st)[!is.na(obj_st$hsc_state_internal)]
hsc <- subset(obj_st, cells = hsc_cells)
hsc$hsc_state_internal <- factor(
  as.character(hsc$hsc_state_internal),
  levels = c("qHSC", "aHSC")
)

# ------------------------------------------------------------------------------
# 6. Add HSC module scores
# ------------------------------------------------------------------------------
add_named_module <- function(object, genes, output_name) {
  genes_use <- existing_genes(object, genes, ASSAY_USE)
  if (length(genes_use) < 3L) {
    warning("Insufficient genes for ", output_name, ": ",
            paste(genes_use, collapse = ", "))
    object[[output_name]] <- NA_real_
    return(object)
  }

  object <- AddModuleScore(
    object,
    features = list(genes_use),
    assay = ASSAY_USE,
    name = paste0(output_name, "_tmp"),
    seed = 20260731,
    search = FALSE
  )

  created <- paste0(output_name, "_tmp1")
  object[[output_name]] <- object@meta.data[[created]]
  object[[created]] <- NULL
  object
}

hsc <- add_named_module(hsc, ahsc_genes, "aHSC_score")
hsc <- add_named_module(hsc, qhsc_genes, "qHSC_score")
hsc <- add_named_module(hsc, ecm_production_genes, "ECM_production_score")
hsc <- add_named_module(hsc, ecm_crosslinking_genes, "ECM_crosslinking_score")
hsc <- add_named_module(hsc, il10_downstream_genes, "IL10_downstream_score")
hsc <- add_named_module(
  hsc,
  il10_antifibrotic_response_genes,
  "IL10_antifibrotic_response_score"
)

# Whole-cell receptor composite score
obj_st <- add_named_module(
  obj_st,
  il10_receptor_genes,
  "IL10_receptor_score"
)

# ------------------------------------------------------------------------------
# 7. HSC abundance analyses
# ------------------------------------------------------------------------------
hsc_counts <- hsc@meta.data %>%
  rownames_to_column("cell") %>%
  count(
    condition_internal,
    sample_internal,
    hsc_state_internal,
    name = "n_cells",
    .drop = FALSE
  ) %>%
  group_by(sample_internal) %>%
  mutate(
    total_HSC = sum(n_cells),
    fraction_within_HSC = n_cells / total_HSC
  ) %>%
  ungroup()

whole_counts <- obj_st@meta.data %>%
  rownames_to_column("cell") %>%
  count(condition_internal, sample_internal, name = "n_all_cells")

hsc_fraction_all_cells <- hsc_counts %>%
  left_join(whole_counts, by = c("condition_internal", "sample_internal")) %>%
  mutate(fraction_of_all_cells = n_cells / n_all_cells)

write.csv(
  hsc_counts,
  file.path(TABLE_DIR, "01_HSC_counts_and_fraction_within_HSC.csv"),
  row.names = FALSE
)
write.csv(
  hsc_fraction_all_cells,
  file.path(TABLE_DIR, "02_HSC_fraction_of_all_cells.csv"),
  row.names = FALSE
)

p_hsc_composition <- ggplot(
  hsc_counts,
  aes(
    x = sample_internal,
    y = fraction_within_HSC,
    fill = hsc_state_internal
  )
) +
  geom_col(width = 0.75) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "qHSC/aHSC composition in Sham and Tx samples",
    x = NULL,
    y = "Fraction within HSC",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_pdf("01_HSC_qHSC_aHSC_composition_by_sample.pdf", p_hsc_composition, 8.5, 6)

# ------------------------------------------------------------------------------
# 8. HSC module-score summaries and plots
# ------------------------------------------------------------------------------
score_cols <- c(
  "aHSC_score",
  "qHSC_score",
  "ECM_production_score",
  "ECM_crosslinking_score",
  "IL10_downstream_score",
  "IL10_antifibrotic_response_score"
)

hsc_score_sample_summary <- hsc@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(
    condition_internal,
    sample_internal,
    hsc_state_internal
  ) %>%
  summarise(
    across(
      all_of(score_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      )
    ),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  hsc_score_sample_summary,
  file.path(TABLE_DIR, "03_HSC_module_score_by_sample.csv"),
  row.names = FALSE
)

make_violin <- function(feature, title, ylab) {
  FetchData(
    hsc,
    vars = c(feature, "condition_internal", "hsc_state_internal")
  ) %>%
    ggplot(
      aes(
        x = condition_internal,
        y = .data[[feature]],
        fill = condition_internal
      )
    ) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.75) +
    facet_wrap(~ hsc_state_internal, nrow = 1) +
    scale_fill_manual(values = c("Sham" = "#6F6F6F", "Tx" = "#FF8C1A")) +
    labs(title = title, x = NULL, y = ylab, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

p_ecm <- make_violin(
  "ECM_production_score",
  "HSC ECM-production score: Sham vs Tx",
  "ECM-production module score"
)
p_cross <- make_violin(
  "ECM_crosslinking_score",
  "HSC ECM-crosslinking score: Sham vs Tx",
  "ECM-crosslinking module score"
)
p_il10 <- make_violin(
  "IL10_downstream_score",
  "HSC IL10-downstream response: Sham vs Tx",
  "IL10-downstream module score"
)
p_il10_af <- make_violin(
  "IL10_antifibrotic_response_score",
  "HSC IL10-associated antifibrotic response: Sham vs Tx",
  "IL10-associated antifibrotic score"
)

save_pdf("02_HSC_ECM_production_violin.pdf", p_ecm, 9, 6)
save_pdf("03_HSC_ECM_crosslinking_violin.pdf", p_cross, 9, 6)
save_pdf("04_HSC_IL10_downstream_violin.pdf", p_il10, 9, 6)
save_pdf("05_HSC_IL10_antifibrotic_response_violin.pdf", p_il10_af, 9, 6)

# Sample-level means, emphasized because biological n = 2 per condition
score_long <- hsc_score_sample_summary %>%
  select(
    condition_internal,
    sample_internal,
    hsc_state_internal,
    ends_with("_mean")
  ) %>%
  pivot_longer(
    cols = ends_with("_mean"),
    names_to = "module",
    values_to = "sample_mean"
  ) %>%
  mutate(module = sub("_mean$", "", module))

p_sample_means <- ggplot(
  score_long,
  aes(
    x = condition_internal,
    y = sample_mean,
    group = sample_internal,
    color = condition_internal
  )
) +
  geom_point(size = 3) +
  facet_grid(hsc_state_internal ~ module, scales = "free_y") +
  scale_color_manual(values = c("Sham" = "#6F6F6F", "Tx" = "#FF8C1A")) +
  labs(
    title = "Sample-level HSC module-score means",
    x = NULL,
    y = "Mean module score",
    color = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_pdf("06_HSC_module_scores_sample_level.pdf", p_sample_means, 16, 7.5)

# ------------------------------------------------------------------------------
# 9. HSC gene-level DotPlots
# ------------------------------------------------------------------------------
hsc_dot_genes <- existing_genes(
  hsc,
  unique(c(
    ahsc_genes,
    qhsc_genes,
    ecm_production_genes,
    ecm_crosslinking_genes,
    "Il10ra", "Il10rb",
    il10_downstream_genes
  )),
  ASSAY_USE
)

hsc$group_internal <- interaction(
  hsc$hsc_state_internal,
  hsc$condition_internal,
  sep = "_",
  drop = TRUE
)
Idents(hsc) <- hsc$group_internal

p_hsc_dot <- DotPlot(
  hsc,
  features = hsc_dot_genes,
  assay = ASSAY_USE,
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "HSC state, ECM and IL10-response markers",
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 7)
  )

save_pdf("07_HSC_state_ECM_IL10_DotPlot.pdf", p_hsc_dot, 22, 6.5)

# ------------------------------------------------------------------------------
# 10. Whole-cell UMAPs for Il10 and receptors
# ------------------------------------------------------------------------------
whole_features <- existing_genes(
  obj_st,
  c("Il10", "Il10ra", "Il10rb"),
  ASSAY_USE
)

if (length(whole_features) == 0L) {
  stop("Il10/Il10ra/Il10rb were not found in the RNA assay.")
}

feature_plots <- FeaturePlot(
  obj_st,
  features = whole_features,
  reduction = umap_reduction,
  split.by = NULL,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  pt.size = 0.45,
  combine = FALSE
)

for (i in seq_along(feature_plots)) {
  feature_plots[[i]] <- feature_plots[[i]] +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p_il10_all <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(
    title = "Whole-cell Il10 and IL10-receptor expression",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )

save_pdf("08_whole_cell_Il10_Il10ra_Il10rb_UMAP.pdf", p_il10_all, 15, 5.5)

# Separate condition-specific UMAP panels
split_feature_plots <- FeaturePlot(
  obj_st,
  features = whole_features,
  reduction = umap_reduction,
  split.by = "condition_internal",
  keep.scale = "all",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  pt.size = 0.40,
  combine = FALSE
)

p_split <- wrap_plots(split_feature_plots, ncol = 2) +
  plot_annotation(
    title = "Whole-cell Il10/Il10ra/Il10rb expression: Sham vs Tx",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )

save_pdf(
  "09_whole_cell_Il10_Il10ra_Il10rb_UMAP_split_Sham_Tx.pdf",
  p_split,
  12,
  14
)

# IL10 receptor module score UMAP
p_receptor_score <- FeaturePlot(
  obj_st,
  features = "IL10_receptor_score",
  reduction = umap_reduction,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q99",
  cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
  pt.size = 0.45
) +
  labs(title = "Whole-cell IL10 receptor composite score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_pdf("10_whole_cell_IL10_receptor_score_UMAP.pdf", p_receptor_score, 8.5, 7)

# Cell type UMAP for reference
p_celltype <- DimPlot(
  obj_st,
  reduction = umap_reduction,
  group.by = "celltype_internal",
  pt.size = 0.40,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  labs(title = "Whole-cell RDS3 cell-type annotation: Sham and Tx") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

save_pdf("11_whole_cell_celltype_UMAP_Sham_Tx.pdf", p_celltype, 13, 9)

# ------------------------------------------------------------------------------
# 11. Whole-cell expression summary by cell type
# ------------------------------------------------------------------------------
DefaultAssay(obj_st) <- ASSAY_USE

expr_summary <- FetchData(
  obj_st,
  vars = c(
    "Il10", "Il10ra", "Il10rb", "IL10_receptor_score",
    "condition_internal", "sample_internal", "celltype_internal"
  )
) %>%
  group_by(
    condition_internal,
    sample_internal,
    celltype_internal
  ) %>%
  summarise(
    n_cells = n(),
    Il10_mean = mean(Il10, na.rm = TRUE),
    Il10_pct_positive = mean(Il10 > 0, na.rm = TRUE),
    Il10ra_mean = mean(Il10ra, na.rm = TRUE),
    Il10ra_pct_positive = mean(Il10ra > 0, na.rm = TRUE),
    Il10rb_mean = mean(Il10rb, na.rm = TRUE),
    Il10rb_pct_positive = mean(Il10rb > 0, na.rm = TRUE),
    IL10_receptor_score_mean = mean(IL10_receptor_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  expr_summary,
  file.path(TABLE_DIR, "04_whole_cell_IL10_receptor_expression_by_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Exploratory cell-level HSC comparisons
# ------------------------------------------------------------------------------
celllevel_tests <- bind_rows(lapply(score_cols, function(score_name) {
  bind_rows(lapply(c("qHSC", "aHSC"), function(state) {
    df <- hsc@meta.data %>%
      filter(hsc_state_internal == state)

    sham <- df[df$condition_internal == "Sham", score_name, drop = TRUE]
    tx <- df[df$condition_internal == "Tx", score_name, drop = TRUE]

    wt <- tryCatch(
      wilcox.test(tx, sham, exact = FALSE),
      error = function(e) NULL
    )

    tibble(
      hsc_state = state,
      feature = score_name,
      Sham_mean = mean(sham, na.rm = TRUE),
      Tx_mean = mean(tx, na.rm = TRUE),
      difference_Tx_minus_Sham =
        mean(tx, na.rm = TRUE) - mean(sham, na.rm = TRUE),
      p_value_celllevel = if (is.null(wt)) NA_real_ else wt$p.value,
      note = "Exploratory cell-level test; biological replication is n=2 per condition."
    )
  }))
}))

celllevel_tests <- celllevel_tests %>%
  mutate(FDR_celllevel = p.adjust(p_value_celllevel, method = "BH"))

write.csv(
  celllevel_tests,
  file.path(TABLE_DIR, "05_exploratory_HSC_celllevel_module_tests.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Save RDS and Excel
# ------------------------------------------------------------------------------
saveRDS(
  obj_st,
  file.path(RDS_DIR, "RDS3_allcells_ShamTx_IL10_receptor_scored_v1.0.0.rds"),
  compress = FALSE
)

saveRDS(
  hsc,
  file.path(RDS_DIR, "RDS3_HSC_ShamTx_ECM_IL10_scored_v1.0.0.rds"),
  compress = FALSE
)

write_excel_sheets(
  list(
    Celltype_labels = celltype_values,
    HSC_counts = hsc_counts,
    HSC_fraction_allcells = hsc_fraction_all_cells,
    HSC_module_sample = hsc_score_sample_summary,
    Wholecell_IL10 = expr_summary,
    HSC_celllevel_tests = celllevel_tests
  ),
  file.path(TABLE_DIR, "RDS3_ShamTx_HSC_ECM_IL10_summary_v1.0.0.xlsx")
)

# ------------------------------------------------------------------------------
# 14. Log
# ------------------------------------------------------------------------------
log_lines <- c(
  "analysis_name: RDS3 whole-cell Sham vs Tx HSC/ECM/IL10 analysis",
  "version: 1.0.0",
  paste0("run_time: ", Sys.time()),
  paste0("input_rds: ", INPUT_RDS),
  paste0("output_dir: ", OUTPUT_DIR),
  paste0("celltype_col: ", CELLTYPE_COL),
  paste0("condition_col: ", CONDITION_COL),
  paste0("sample_col: ", SAMPLE_COL),
  paste0("umap_reduction: ", umap_reduction),
  paste0("all_ShamTx_cells: ", ncol(obj_st)),
  paste0("HSC_cells: ", ncol(hsc)),
  "primary comparison: Sham1/Sham20 vs Tx17/Tx5",
  "cell-level p-values are exploratory only",
  "sample-level summaries should guide biological interpretation",
  "sessionInfo:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)

writeLines(
  log_lines,
  file.path(LOG_DIR, "run_manifest_v1.0.0.txt")
)

message_time("Completed RDS3 Sham/Tx HSC, ECM and IL10 analysis.")
message_time("Output directory: ", OUTPUT_DIR)
