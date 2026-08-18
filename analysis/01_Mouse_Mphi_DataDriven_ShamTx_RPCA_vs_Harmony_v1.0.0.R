#!/usr/bin/env Rscript

# ==============================================================================
# uenoyscRNA
# Mouse macrophage data-driven analysis Step 1
# Sham vs Tx reference construction: RPCA versus Harmony
# Version 1.0.0
# ==============================================================================
#
# PRIMARY ANALYSIS:
#   Sham1, Sham20 versus Tx17, Tx5
#
# REFERENCE DATA:
#   STD_rep1 versus CDHFD_rep1
#   These cells are NOT used to construct the Sham/Tx reference.
#   They are saved separately for later FindTransferAnchors()/MapQuery().
#
# Design:
#   1. Read the whole-liver Seurat RDS.
#   2. Extract Layer1 Kupffer_Macrophage and Monocyte cells.
#   3. Split into:
#        - primary reference: Sham/Tx
#        - external query: STD/CDAHFD
#   4. Run identical preprocessing on the Sham/Tx reference.
#   5. Run both:
#        - Seurat RPCAIntegration
#        - Seurat HarmonyIntegration
#   6. Recluster independently at several resolutions.
#   7. Save UMAPs, cluster/sample audits, cross-method agreement tables,
#      and both annotated Seurat objects.
#
# Important:
#   - condition is never used as a batch variable.
#   - sample_internal is the integration batch.
#   - raw/count and normalized RNA expression are retained for downstream DEG.
#   - RPCA/Harmony embeddings are used only for neighbors, clustering, and UMAP.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260731)

# ------------------------------------------------------------------------------
# 0. Paths and settings
# ------------------------------------------------------------------------------
PROJECT_ROOT <- "/Users/uenoya/Projects/uenoyscRNA"

INPUT_RDS <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  "DataDriven_ShamTx_RPCA_vs_Harmony_v1.0.0"
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
FIGURE_DIR <- file.path(OUTPUT_DIR, "Figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

for (d in c(RDS_DIR, FIGURE_DIR, TABLE_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LAYER1_SOURCE_COL <- "celltype_for_R8plot_FIXED2"
CONDITION_SOURCE_COL <- "condition"
SAMPLE_SOURCE_COL <- "sample"

TARGET_LAYER1 <- c("Kupffer_Macrophage", "Monocyte")
PRIMARY_CONDITIONS <- c("Sham", "Tx")
REFERENCE_CONDITIONS <- c("STD", "CDAHFD")
CONDITION_LEVELS <- c("STD", "CDAHFD", "Sham", "Tx")

CONDITION_ALIASES <- list(
  STD = c("STD", "STANDARD", "CONTROL", "CTRL"),
  CDAHFD = c("CDAHFD", "CDHFD", "CDA-HFD", "CDA_HFD", "CDA HFD"),
  Sham = c("SHAM", "SHAM1", "SHAM20"),
  Tx = c("TX", "TX5", "TX17", "TREATMENT", "TRANSPLANT")
)

PREFERRED_ASSAY <- "RNA"

NFEATURES <- 3000L
NPCS <- 50L
DIMS_USE <- 1:30
K_PARAM <- 20L
RESOLUTIONS <- c(0.8, 1.0, 1.2, 1.5)
PRIMARY_RESOLUTION <- 1.2
UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.30

PCA_REDUCTION <- "mphi.pca"
RPCA_REDUCTION <- "mphi.integrated.rpca"
HARMONY_REDUCTION <- "mphi.integrated.harmony"
RPCA_UMAP <- "mphi.umap.rpca"
HARMONY_UMAP <- "mphi.umap.harmony"

SAMPLE_COLORS <- c(
  "Sham1" = "#5B5B5B",
  "Sham20" = "#A0A0A0",
  "Tx17" = "#FF8C1A",
  "Tx5" = "#FFB45A"
)

CONDITION_COLORS <- c(
  "Sham" = "#7A7A7A",
  "Tx" = "#FF8C1A"
)

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
  "harmony"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
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
})

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------
message_time <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

save_pdf <- function(filename, plot, width, height) {
  ggplot2::ggsave(
    filename = file.path(FIGURE_DIR, filename),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

canonicalize_token <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]", "", x)
}

canonicalize_condition <- function(x) {
  token <- canonicalize_token(x)
  out <- rep(NA_character_, length(token))

  for (canonical_name in names(CONDITION_ALIASES)) {
    aliases <- unique(canonicalize_token(CONDITION_ALIASES[[canonical_name]]))
    out[token %in% aliases] <- canonical_name
  }
  out
}

validate_columns <- function(object, columns) {
  missing_cols <- setdiff(columns, colnames(object@meta.data))
  if (length(missing_cols) > 0L) {
    stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
  }
}

join_layers_safe <- function(object, assay) {
  layer_names <- SeuratObject::Layers(object[[assay]])
  split_prefixes <- sub("\\..*$", "", layer_names)

  if (any(duplicated(split_prefixes)) || length(grep("^counts\\.", layer_names)) > 0L) {
    message_time("Joining existing ", assay, " layers before re-splitting.")
    object <- SeuratObject::JoinLayers(object, assay = assay)
  }
  object
}

prepare_layered_object <- function(object, assay, batch_col) {
  DefaultAssay(object) <- assay
  object <- join_layers_safe(object, assay)

  batch <- as.factor(object@meta.data[[batch_col]])
  if (anyNA(batch) || any(!nzchar(as.character(batch)))) {
    stop("Missing batch labels in ", batch_col)
  }
  if (nlevels(batch) < 2L) {
    stop("At least two samples are required for integration.")
  }

  object[[assay]] <- split(object[[assay]], f = batch)

  object <- NormalizeData(
    object,
    assay = assay,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  object <- FindVariableFeatures(
    object,
    assay = assay,
    selection.method = "vst",
    nfeatures = NFEATURES,
    verbose = FALSE
  )

  object <- ScaleData(
    object,
    assay = assay,
    features = VariableFeatures(object, assay = assay),
    verbose = FALSE
  )

  object <- RunPCA(
    object,
    assay = assay,
    features = VariableFeatures(object, assay = assay),
    npcs = NPCS,
    reduction.name = PCA_REDUCTION,
    reduction.key = "MphiPC_",
    verbose = FALSE
  )

  object
}

resolution_tag <- function(x) {
  gsub("\\.", "_", format(x, trim = TRUE, scientific = FALSE))
}

run_clustering_grid <- function(
    object,
    reduction,
    umap_name,
    graph_prefix,
    cluster_prefix,
    method_label) {

  dims_available <- ncol(Embeddings(object, reduction = reduction))
  dims_use <- DIMS_USE[DIMS_USE <= dims_available]

  if (length(dims_use) < 5L) {
    stop(
      method_label,
      " reduction has too few dimensions: ",
      dims_available
    )
  }

  object <- FindNeighbors(
    object,
    reduction = reduction,
    dims = dims_use,
    k.param = K_PARAM,
    graph.name = c(
      paste0(graph_prefix, "_nn"),
      paste0(graph_prefix, "_snn")
    ),
    verbose = FALSE
  )

  for (res in RESOLUTIONS) {
    col_name <- paste0(cluster_prefix, "_res_", resolution_tag(res))

    object <- FindClusters(
      object,
      graph.name = paste0(graph_prefix, "_snn"),
      resolution = res,
      algorithm = 1,
      random.seed = 20260731,
      cluster.name = col_name,
      verbose = FALSE
    )
  }

  primary_col <- paste0(
    cluster_prefix,
    "_res_",
    resolution_tag(PRIMARY_RESOLUTION)
  )

  object[[paste0(cluster_prefix, "_primary")]] <-
    object@meta.data[[primary_col]]

  object <- RunUMAP(
    object,
    reduction = reduction,
    dims = dims_use,
    reduction.name = umap_name,
    reduction.key = if (method_label == "RPCA") "MphiRPCAUMAP_" else "MphiHarmonyUMAP_",
    n.neighbors = UMAP_N_NEIGHBORS,
    min.dist = UMAP_MIN_DIST,
    seed.use = 20260731,
    verbose = FALSE
  )

  object
}

make_umap_df <- function(object, reduction, cluster_col, method_label) {
  emb <- Embeddings(object, reduction = reduction)

  object@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      UMAP_1 = emb[cell, 1],
      UMAP_2 = emb[cell, 2],
      cluster = factor(.data[[cluster_col]]),
      method = method_label,
      condition_internal = factor(
        condition_internal,
        levels = PRIMARY_CONDITIONS
      )
    )
}

make_umap_plot <- function(
    df,
    color_col,
    title,
    colors = NULL,
    split_col = NULL,
    point_size = 0.42) {

  p <- ggplot(
    df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = .data[[color_col]]
    )
  ) +
    geom_point(size = point_size, alpha = 0.90) +
    coord_equal() +
    labs(
      title = title,
      x = "UMAP_1",
      y = "UMAP_2",
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.key.height = grid::unit(0.45, "cm")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3.2, alpha = 1)
      )
    )

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors, drop = FALSE)
  }

  if (!is.null(split_col)) {
    p <- p + facet_wrap(
      stats::as.formula(paste("~", split_col)),
      nrow = 1
    )
  }

  p
}

cluster_sample_audit <- function(object, cluster_col, method_label) {
  object@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      condition_internal,
      sample_internal,
      cluster = .data[[cluster_col]],
      name = "n_cells",
      .drop = FALSE
    ) %>%
    group_by(cluster) %>%
    mutate(
      cluster_total = sum(n_cells),
      fraction_within_cluster = if_else(
        cluster_total > 0,
        n_cells / cluster_total,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(method = method_label)
}

cluster_condition_composition <- function(object, cluster_col, method_label) {
  object@meta.data %>%
    rownames_to_column("cell") %>%
    count(
      condition_internal,
      sample_internal,
      cluster = .data[[cluster_col]],
      name = "n_cells",
      .drop = FALSE
    ) %>%
    group_by(condition_internal, sample_internal) %>%
    mutate(
      sample_total = sum(n_cells),
      fraction_of_sample = if_else(
        sample_total > 0,
        n_cells / sample_total,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(method = method_label)
}

entropy_normalized <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) <= 1L) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p)) / log(length(p))
}

sample_mixing_summary <- function(audit_df) {
  audit_df %>%
    group_by(method, cluster) %>%
    summarise(
      n_cells = sum(n_cells),
      n_samples_present = sum(n_cells > 0),
      dominant_sample_fraction = max(
        if_else(sum(n_cells) > 0, n_cells / sum(n_cells), NA_real_),
        na.rm = TRUE
      ),
      normalized_sample_entropy = entropy_normalized(n_cells),
      .groups = "drop"
    )
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)

  choose2 <- function(z) z * (z - 1) / 2

  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  n <- sum(tab)
  total_pairs <- choose2(n)

  if (total_pairs == 0) return(NA_real_)

  expected <- (sum_ai * sum_bj) / total_pairs
  max_index <- 0.5 * (sum_ai + sum_bj)
  denominator <- max_index - expected

  if (denominator == 0) return(NA_real_)
  (sum_nij - expected) / denominator
}

best_jaccard_mapping <- function(rpca_cluster, harmony_cluster) {
  tab <- table(
    RPCA = as.character(rpca_cluster),
    Harmony = as.character(harmony_cluster)
  )

  out <- list()
  k <- 1L

  for (r in rownames(tab)) {
    for (h in colnames(tab)) {
      intersection <- tab[r, h]
      union <- sum(tab[r, ]) + sum(tab[, h]) - intersection
      out[[k]] <- data.frame(
        rpca_cluster = r,
        harmony_cluster = h,
        intersection_n = as.integer(intersection),
        jaccard = if (union > 0) intersection / union else NA_real_
      )
      k <- k + 1L
    }
  }

  all_pairs <- bind_rows(out)

  best_rpca_to_harmony <- all_pairs %>%
    group_by(rpca_cluster) %>%
    slice_max(jaccard, n = 1, with_ties = FALSE) %>%
    ungroup()

  best_harmony_to_rpca <- all_pairs %>%
    group_by(harmony_cluster) %>%
    slice_max(jaccard, n = 1, with_ties = FALSE) %>%
    ungroup()

  list(
    all_pairs = all_pairs,
    rpca_to_harmony = best_rpca_to_harmony,
    harmony_to_rpca = best_harmony_to_rpca
  )
}

# ------------------------------------------------------------------------------
# 3. Read and harmonize metadata
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_RDS)) {
  stop("Input RDS not found: ", INPUT_RDS)
}

message_time("Reading input RDS: ", INPUT_RDS)
whole <- readRDS(INPUT_RDS)
message_time("Input cells: ", ncol(whole))

validate_columns(
  whole,
  c(
    LAYER1_SOURCE_COL,
    CONDITION_SOURCE_COL,
    SAMPLE_SOURCE_COL
  )
)

source_condition_table <- whole@meta.data %>%
  count(
    source_condition = .data[[CONDITION_SOURCE_COL]],
    name = "n_cells",
    sort = TRUE
  )

source_layer1_table <- whole@meta.data %>%
  count(
    source_layer1 = .data[[LAYER1_SOURCE_COL]],
    name = "n_cells",
    sort = TRUE
  )

write.csv(
  source_condition_table,
  file.path(TABLE_DIR, "00_source_condition_values.csv"),
  row.names = FALSE
)
write.csv(
  source_layer1_table,
  file.path(TABLE_DIR, "00_source_layer1_values.csv"),
  row.names = FALSE
)

whole$condition_internal <- canonicalize_condition(
  whole@meta.data[[CONDITION_SOURCE_COL]]
)
whole$sample_internal <- as.character(
  whole@meta.data[[SAMPLE_SOURCE_COL]]
)
whole$layer1_internal <- as.character(
  whole@meta.data[[LAYER1_SOURCE_COL]]
)

keep <- (
  whole$layer1_internal %in% TARGET_LAYER1 &
  whole$condition_internal %in% c(
    PRIMARY_CONDITIONS,
    REFERENCE_CONDITIONS
  )
)

mphi_all <- subset(whole, cells = colnames(whole)[keep])
rm(whole)
invisible(gc())

mphi_all$condition_internal <- factor(
  mphi_all$condition_internal,
  levels = CONDITION_LEVELS
)

message_time(
  "Extracted macrophage/monocyte cells: ",
  ncol(mphi_all)
)
message_time(
  "Counts by condition: ",
  paste(
    names(table(mphi_all$condition_internal)),
    as.integer(table(mphi_all$condition_internal)),
    collapse = "; "
  )
)

primary_cells <- colnames(mphi_all)[
  as.character(mphi_all$condition_internal) %in% PRIMARY_CONDITIONS
]
reference_cells <- colnames(mphi_all)[
  as.character(mphi_all$condition_internal) %in% REFERENCE_CONDITIONS
]

primary_base <- subset(mphi_all, cells = primary_cells)
reference_query <- subset(mphi_all, cells = reference_cells)
rm(mphi_all)
invisible(gc())

primary_base$condition_internal <- factor(
  as.character(primary_base$condition_internal),
  levels = PRIMARY_CONDITIONS
)

sample_order <- c("Sham1", "Sham20", "Tx17", "Tx5")
observed_primary_samples <- unique(as.character(primary_base$sample_internal))
missing_primary_samples <- setdiff(sample_order, observed_primary_samples)

if (length(missing_primary_samples) > 0L) {
  stop(
    "Missing expected Sham/Tx sample(s): ",
    paste(missing_primary_samples, collapse = ", "),
    "\nObserved samples: ",
    paste(sort(observed_primary_samples), collapse = ", ")
  )
}

primary_base$sample_internal <- factor(
  as.character(primary_base$sample_internal),
  levels = sample_order
)

reference_query$condition_internal <- factor(
  as.character(reference_query$condition_internal),
  levels = REFERENCE_CONDITIONS
)

write.csv(
  primary_base@meta.data %>%
    rownames_to_column("cell") %>%
    count(condition_internal, sample_internal, name = "n_cells"),
  file.path(TABLE_DIR, "01_primary_Sham_Tx_cell_counts.csv"),
  row.names = FALSE
)

write.csv(
  reference_query@meta.data %>%
    rownames_to_column("cell") %>%
    count(condition_internal, sample_internal, name = "n_cells"),
  file.path(TABLE_DIR, "01_reference_STD_CDAHFD_cell_counts.csv"),
  row.names = FALSE
)

saveRDS(
  reference_query,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_STD_CDAHFD_query_before_reference_mapping_v1.0.0.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 4. Shared Sham/Tx preprocessing
# ------------------------------------------------------------------------------
assay_used <- if (PREFERRED_ASSAY %in% Assays(primary_base)) {
  PREFERRED_ASSAY
} else {
  DefaultAssay(primary_base)
}

message_time("Using assay: ", assay_used)
message_time(
  "Primary Sham/Tx samples: ",
  paste(levels(primary_base$sample_internal), collapse = ", ")
)

message_time("Running shared Sham/Tx preprocessing.")
preprocessed <- prepare_layered_object(
  primary_base,
  assay = assay_used,
  batch_col = "sample_internal"
)
rm(primary_base)
invisible(gc())

saveRDS(
  preprocessed,
  file.path(RDS_DIR, "Mouse_Mphi_Sham_Tx_shared_preprocessed_v1.0.0.rds"),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 5. RPCA integration and clustering
# ------------------------------------------------------------------------------
message_time("Running Seurat RPCAIntegration.")
rpca_obj <- IntegrateLayers(
  object = preprocessed,
  method = RPCAIntegration,
  orig.reduction = PCA_REDUCTION,
  new.reduction = RPCA_REDUCTION,
  assay = assay_used,
  dims = DIMS_USE,
  k.anchor = 5,
  verbose = FALSE
)

rpca_obj <- run_clustering_grid(
  object = rpca_obj,
  reduction = RPCA_REDUCTION,
  umap_name = RPCA_UMAP,
  graph_prefix = "mphi_rpca",
  cluster_prefix = "mphi_rpca",
  method_label = "RPCA"
)

# ------------------------------------------------------------------------------
# 6. Harmony integration and clustering
# ------------------------------------------------------------------------------
message_time("Running Seurat HarmonyIntegration.")
harmony_obj <- IntegrateLayers(
  object = preprocessed,
  method = HarmonyIntegration,
  orig.reduction = PCA_REDUCTION,
  new.reduction = HARMONY_REDUCTION,
  assay = assay_used,
  theta = 2,
  verbose = FALSE
)

harmony_obj <- run_clustering_grid(
  object = harmony_obj,
  reduction = HARMONY_REDUCTION,
  umap_name = HARMONY_UMAP,
  graph_prefix = "mphi_harmony",
  cluster_prefix = "mphi_harmony",
  method_label = "Harmony"
)

rm(preprocessed)
invisible(gc())

rpca_primary_col <- paste0(
  "mphi_rpca_res_",
  resolution_tag(PRIMARY_RESOLUTION)
)
harmony_primary_col <- paste0(
  "mphi_harmony_res_",
  resolution_tag(PRIMARY_RESOLUTION)
)

# ------------------------------------------------------------------------------
# 7. UMAP figures
# ------------------------------------------------------------------------------
rpca_df <- make_umap_df(
  rpca_obj,
  reduction = RPCA_UMAP,
  cluster_col = rpca_primary_col,
  method_label = "RPCA"
)

harmony_df <- make_umap_df(
  harmony_obj,
  reduction = HARMONY_UMAP,
  cluster_col = harmony_primary_col,
  method_label = "Harmony"
)

p_rpca_cluster <- make_umap_plot(
  rpca_df,
  color_col = "cluster",
  title = paste0(
    "Sham/Tx macrophage-only RPCA clusters (resolution ",
    PRIMARY_RESOLUTION,
    ")"
  )
)

p_harmony_cluster <- make_umap_plot(
  harmony_df,
  color_col = "cluster",
  title = paste0(
    "Sham/Tx macrophage-only Harmony clusters (resolution ",
    PRIMARY_RESOLUTION,
    ")"
  )
)

save_pdf("01A_RPCA_UMAP_clusters.pdf", p_rpca_cluster, 9.0, 7.2)
save_pdf("01B_Harmony_UMAP_clusters.pdf", p_harmony_cluster, 9.0, 7.2)

p_rpca_sample <- make_umap_plot(
  rpca_df,
  color_col = "sample_internal",
  title = "RPCA UMAP colored by sample",
  colors = SAMPLE_COLORS
)

p_harmony_sample <- make_umap_plot(
  harmony_df,
  color_col = "sample_internal",
  title = "Harmony UMAP colored by sample",
  colors = SAMPLE_COLORS
)

save_pdf("02A_RPCA_UMAP_by_sample.pdf", p_rpca_sample, 9.0, 7.2)
save_pdf("02B_Harmony_UMAP_by_sample.pdf", p_harmony_sample, 9.0, 7.2)

p_rpca_condition <- make_umap_plot(
  rpca_df,
  color_col = "condition_internal",
  title = "RPCA UMAP: Sham vs Tx",
  colors = CONDITION_COLORS
)

p_harmony_condition <- make_umap_plot(
  harmony_df,
  color_col = "condition_internal",
  title = "Harmony UMAP: Sham vs Tx",
  colors = CONDITION_COLORS
)

save_pdf("03A_RPCA_UMAP_Sham_vs_Tx.pdf", p_rpca_condition, 9.0, 7.2)
save_pdf("03B_Harmony_UMAP_Sham_vs_Tx.pdf", p_harmony_condition, 9.0, 7.2)

p_rpca_split_sample <- make_umap_plot(
  rpca_df,
  color_col = "cluster",
  title = "RPCA clusters split by Sham/Tx sample",
  split_col = "sample_internal",
  point_size = 0.34
)

p_harmony_split_sample <- make_umap_plot(
  harmony_df,
  color_col = "cluster",
  title = "Harmony clusters split by Sham/Tx sample",
  split_col = "sample_internal",
  point_size = 0.34
)

save_pdf(
  "04A_RPCA_UMAP_clusters_split_by_sample.pdf",
  p_rpca_split_sample,
  15.5,
  4.8
)
save_pdf(
  "04B_Harmony_UMAP_clusters_split_by_sample.pdf",
  p_harmony_split_sample,
  15.5,
  4.8
)

# ------------------------------------------------------------------------------
# 8. Cluster/sample audits
# ------------------------------------------------------------------------------
rpca_audit <- cluster_sample_audit(
  rpca_obj,
  rpca_primary_col,
  "RPCA"
)

harmony_audit <- cluster_sample_audit(
  harmony_obj,
  harmony_primary_col,
  "Harmony"
)

cluster_sample_all <- bind_rows(rpca_audit, harmony_audit)
write.csv(
  cluster_sample_all,
  file.path(TABLE_DIR, "02_cluster_by_sample_audit_primary_resolution.csv"),
  row.names = FALSE
)

mixing_summary <- sample_mixing_summary(cluster_sample_all)
write.csv(
  mixing_summary,
  file.path(TABLE_DIR, "03_cluster_sample_mixing_summary.csv"),
  row.names = FALSE
)

rpca_composition <- cluster_condition_composition(
  rpca_obj,
  rpca_primary_col,
  "RPCA"
)
harmony_composition <- cluster_condition_composition(
  harmony_obj,
  harmony_primary_col,
  "Harmony"
)

composition_all <- bind_rows(rpca_composition, harmony_composition)
write.csv(
  composition_all,
  file.path(TABLE_DIR, "04_cluster_fraction_by_sample.csv"),
  row.names = FALSE
)

p_composition <- ggplot(
  composition_all,
  aes(
    x = sample_internal,
    y = fraction_of_sample,
    fill = factor(cluster)
  )
) +
  geom_col(width = 0.78) +
  facet_wrap(~ method, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Cluster composition by Sham/Tx sample",
    x = NULL,
    y = "Fraction of macrophage/monocyte cells",
    fill = "Cluster"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

save_pdf(
  "05_RPCA_vs_Harmony_cluster_composition_by_sample.pdf",
  p_composition,
  12.5,
  6.5
)

# ------------------------------------------------------------------------------
# 9. Resolution audit
# ------------------------------------------------------------------------------
resolution_audit_one <- function(object, prefix, method) {
  bind_rows(lapply(RESOLUTIONS, function(res) {
    col_name <- paste0(prefix, "_res_", resolution_tag(res))

    object@meta.data %>%
      rownames_to_column("cell") %>%
      count(
        condition_internal,
        sample_internal,
        cluster = .data[[col_name]],
        name = "n_cells"
      ) %>%
      mutate(
        method = method,
        resolution = res
      )
  }))
}

resolution_audit <- bind_rows(
  resolution_audit_one(rpca_obj, "mphi_rpca", "RPCA"),
  resolution_audit_one(harmony_obj, "mphi_harmony", "Harmony")
)

write.csv(
  resolution_audit,
  file.path(TABLE_DIR, "05_cluster_counts_resolution_grid.csv"),
  row.names = FALSE
)

cluster_number_summary <- resolution_audit %>%
  distinct(method, resolution, cluster) %>%
  count(method, resolution, name = "n_clusters")

write.csv(
  cluster_number_summary,
  file.path(TABLE_DIR, "06_number_of_clusters_by_resolution.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. RPCA versus Harmony cluster agreement
# ------------------------------------------------------------------------------
common_cells <- intersect(colnames(rpca_obj), colnames(harmony_obj))

rpca_labels <- as.character(
  rpca_obj@meta.data[common_cells, rpca_primary_col]
)
harmony_labels <- as.character(
  harmony_obj@meta.data[common_cells, harmony_primary_col]
)

ari_value <- adjusted_rand_index(rpca_labels, harmony_labels)

agreement_summary <- data.frame(
  primary_resolution = PRIMARY_RESOLUTION,
  n_common_cells = length(common_cells),
  adjusted_rand_index = ari_value
)

write.csv(
  agreement_summary,
  file.path(TABLE_DIR, "07_RPCA_Harmony_adjusted_rand_index.csv"),
  row.names = FALSE
)

cross_tab <- as.data.frame.matrix(
  table(
    RPCA_cluster = rpca_labels,
    Harmony_cluster = harmony_labels
  )
)
cross_tab <- rownames_to_column(cross_tab, "RPCA_cluster")

write.csv(
  cross_tab,
  file.path(TABLE_DIR, "08_RPCA_Harmony_cluster_contingency.csv"),
  row.names = FALSE
)

jaccard <- best_jaccard_mapping(rpca_labels, harmony_labels)

write.csv(
  jaccard$all_pairs,
  file.path(TABLE_DIR, "09_RPCA_Harmony_all_cluster_Jaccard_pairs.csv"),
  row.names = FALSE
)
write.csv(
  jaccard$rpca_to_harmony,
  file.path(TABLE_DIR, "10_best_Harmony_match_for_each_RPCA_cluster.csv"),
  row.names = FALSE
)
write.csv(
  jaccard$harmony_to_rpca,
  file.path(TABLE_DIR, "11_best_RPCA_match_for_each_Harmony_cluster.csv"),
  row.names = FALSE
)

p_jaccard <- ggplot(
  jaccard$all_pairs,
  aes(
    x = factor(harmony_cluster),
    y = factor(rpca_cluster),
    fill = jaccard
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A",
    limits = c(0, 1)
  ) +
  labs(
    title = paste0(
      "RPCA-Harmony cluster correspondence (ARI = ",
      round(ari_value, 3),
      ")"
    ),
    x = "Harmony cluster",
    y = "RPCA cluster",
    fill = "Jaccard"
  ) +
  coord_fixed() +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_pdf(
  "06_RPCA_Harmony_cluster_Jaccard_heatmap.pdf",
  p_jaccard,
  8.5,
  7.5
)

# ------------------------------------------------------------------------------
# 11. Save objects
# ------------------------------------------------------------------------------
saveRDS(
  rpca_obj,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_DataDriven_RPCA_v1.0.0.rds"
  ),
  compress = FALSE
)

saveRDS(
  harmony_obj,
  file.path(
    RDS_DIR,
    "Mouse_Mphi_Sham_Tx_DataDriven_Harmony_v1.0.0.rds"
  ),
  compress = FALSE
)

# Cell-level cross-method table
cell_comparison <- tibble(
  cell = common_cells,
  condition_internal = as.character(
    rpca_obj@meta.data[common_cells, "condition_internal"]
  ),
  sample_internal = as.character(
    rpca_obj@meta.data[common_cells, "sample_internal"]
  ),
  rpca_cluster = rpca_labels,
  harmony_cluster = harmony_labels
)

write.csv(
  cell_comparison,
  file.path(TABLE_DIR, "12_cell_level_RPCA_Harmony_cluster_comparison.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Manifest
# ------------------------------------------------------------------------------
manifest <- c(
  "analysis_name: DataDriven_ShamTx_RPCA_vs_Harmony",
  "script_version: 1.0.0",
  paste0("run_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("input_rds: ", INPUT_RDS),
  paste0("output_dir: ", OUTPUT_DIR),
  paste0("assay_used: ", assay_used),
  "primary_conditions: Sham, Tx",
  "primary_samples: Sham1, Sham20, Tx17, Tx5",
  "reference_conditions_saved_for_mapping: STD, CDAHFD",
  "condition_used_as_integration_batch: FALSE",
  "integration_batch: sample_internal",
  paste0("nfeatures: ", NFEATURES),
  paste0("npcs: ", NPCS),
  paste0("dims_used: ", paste(range(DIMS_USE), collapse = "-")),
  paste0("k_param: ", K_PARAM),
  paste0("resolutions: ", paste(RESOLUTIONS, collapse = ",")),
  paste0("primary_resolution: ", PRIMARY_RESOLUTION),
  paste0("RPCA_Harmony_ARI_primary_resolution: ", ari_value),
  paste0("RPCA_cells: ", ncol(rpca_obj)),
  paste0("Harmony_cells: ", ncol(harmony_obj)),
  paste0("reference_query_cells: ", ncol(reference_query)),
  "DEG_policy: use RNA counts/pseudobulk; do not use integrated embeddings as expression",
  "sessionInfo:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)

writeLines(
  manifest,
  file.path(LOG_DIR, "run_manifest_v1.0.0.txt")
)

message_time("Completed Sham/Tx RPCA versus Harmony Step 1.")
message_time("RPCA object: ", file.path(
  RDS_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_RPCA_v1.0.0.rds"
))
message_time("Harmony object: ", file.path(
  RDS_DIR,
  "Mouse_Mphi_Sham_Tx_DataDriven_Harmony_v1.0.0.rds"
))
message_time("STD/CDAHFD query object: ", file.path(
  RDS_DIR,
  "Mouse_Mphi_STD_CDAHFD_query_before_reference_mapping_v1.0.0.rds"
))
message_time("All outputs: ", OUTPUT_DIR)
