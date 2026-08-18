#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH macrophage analysis v4.0.0
# ==============================================================================
# Purpose:
#   1) Extract STD/CDAHFD Kupffer_Macrophage and Monocyte cells from Layer1.
#   2) Reconstruct Layer2 annotation from zero with clear Resident/Monocyte split.
#   3) Apply a strict M1 definition that minimizes Monocyte contamination.
#   4) Exclude the isolated lower-left UMAP population from DISPLAY ONLY.
#   5) Calculate strict M1/M2 proportions using strict M1 + strict M2 as denominator.
#
# Input metadata:
#   layer1    = celltype_for_R8plot_FIXED2
#   condition = condition
#   sample    = sample
#
# Internal metadata names created by this script:
#   condition_internal
#   sample_internal
#   cluster_internal
#   layer2_internal
#
# Layer2 labels:
#   Resident Kupffer-like
#   Monocyte-like
#   Inflammatory M1-like
#   Pro-resolution M2-like
#   SPP1/TREM2 MASH-associated
#   Other
#
# Important:
#   - display_excluded_left_bottom is used only for UMAP plotting.
#   - All quantitative summaries retain these cells.
#   - Cell-level labels are rule-based and intentionally conservative.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4000)

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------
INPUT_RDS <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  "STD_vs_CDAHFD_Layer2_v4.0.0"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source metadata names
LAYER1_SOURCE_COL <- "celltype_for_R8plot_FIXED2"
CONDITION_SOURCE_COL <- "condition"
SAMPLE_SOURCE_COL <- "sample"

# Analysis targets
TARGET_CONDITIONS <- c("STD", "CDAHFD")
TARGET_LAYER1 <- c("Kupffer_Macrophage", "Monocyte")

# Preferred assay and UMAP reduction
PREFERRED_ASSAY <- "RNA"
PREFERRED_UMAP_REDUCTION <- "umap"

# Existing cluster column candidates, in priority order
CLUSTER_COLUMN_CANDIDATES <- c(
  "seurat_clusters",
  "RNA_snn_res.3",
  "RNA_snn_res.2.5",
  "integrated_snn_res.3",
  "integrated_snn_res.2.5"
)

# Lower-left isolated UMAP population: DISPLAY ONLY
# mode = "auto", "manual", or "none"
DISPLAY_EXCLUDE_MODE <- "auto"

# Manual mode excludes cells satisfying BOTH thresholds below.
# Set after inspecting 01_UMAP_display_exclusion_diagnostic.pdf if needed.
MANUAL_UMAP1_MAX <- -Inf
MANUAL_UMAP2_MAX <- -Inf

# Auto mode uses k-means only to identify a small, isolated lower-left component.
AUTO_KMEANS_CENTERS <- 10L
AUTO_MAX_EXCLUDED_FRACTION <- 0.10
AUTO_MIN_EXCLUDED_CELLS <- 20L
AUTO_MIN_CENTROID_SEPARATION_Z <- 1.10

# Strict classification parameters
SCORE_MARGIN_RESIDENT_MONOCYTE <- 0.12
SCORE_MARGIN_M1_M2 <- 0.10
MIN_RESIDENT_MARKERS <- 3L
MIN_MONOCYTE_CORE_MARKERS <- 2L
MIN_M1_MARKERS <- 3L
MIN_M2_MARKERS <- 3L
MIN_SPP1_TREM2_MARKERS <- 3L

# R8-like high-saturation Layer2 colors
LAYER2_COLORS <- c(
  "Resident Kupffer-like" = "#00BFC4",
  "Monocyte-like" = "#2F65FF",
  "Inflammatory M1-like" = "#F04444",
  "Pro-resolution M2-like" = "#00B85A",
  "SPP1/TREM2 MASH-associated" = "#E83E9B",
  "Other" = "#B8B8B8"
)

CONDITION_COLORS <- c(
  "STD" = "#2F65FF",
  "CDAHFD" = "#F04444"
)

LAYER2_LEVELS <- names(LAYER2_COLORS)

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
  "scales"
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
# 2. Helper functions
# ------------------------------------------------------------------------------
message_time <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

save_pdf <- function(filename, plot, width, height) {
  ggplot2::ggsave(
    filename = file.path(OUTPUT_DIR, filename),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

validate_metadata <- function(object, columns) {
  missing_cols <- setdiff(columns, colnames(object@meta.data))
  if (length(missing_cols) > 0L) {
    stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
  }
}

choose_cluster_column <- function(object, candidates) {
  hit <- candidates[candidates %in% colnames(object@meta.data)]
  if (length(hit) > 0L) {
    return(hit[[1]])
  }
  return(NULL)
}

choose_umap_reduction <- function(object, preferred = "umap") {
  reductions <- Reductions(object)
  if (preferred %in% reductions) {
    return(preferred)
  }
  umap_like <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  if (length(umap_like) > 0L) {
    return(umap_like[[1]])
  }
  stop("No UMAP reduction was found in the Seurat object.")
}

prepare_rna_assay <- function(object, preferred_assay = "RNA") {
  if (preferred_assay %in% Assays(object)) {
    DefaultAssay(object) <- preferred_assay
  }

  assay_name <- DefaultAssay(object)
  message_time("Using assay: ", assay_name)

  # Seurat v5 objects can contain multiple assay layers. Join when possible.
  if (exists("Layers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    assay_layers <- tryCatch(
      SeuratObject::Layers(object[[assay_name]]),
      error = function(e) character(0)
    )

    data_like_layers <- assay_layers[grepl("^data($|\\.)", assay_layers)]
    count_like_layers <- assay_layers[grepl("^counts($|\\.)", assay_layers)]

    if (length(data_like_layers) > 1L || length(count_like_layers) > 1L) {
      if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
        message_time("Joining split assay layers for ", assay_name, ".")
        object <- SeuratObject::JoinLayers(object, assay = assay_name)
      } else if (exists("JoinLayers", where = asNamespace("Seurat"), inherits = FALSE)) {
        message_time("Joining split assay layers for ", assay_name, ".")
        object <- Seurat::JoinLayers(object, assay = assay_name)
      }
    }
  }

  # Ensure normalized data exists.
  normalized_ok <- tryCatch({
    x <- GetAssayData(object, assay = assay_name, slot = "data")
    nrow(x) > 0L && ncol(x) > 0L
  }, error = function(e) FALSE)

  if (!normalized_ok) {
    message_time("Normalized data layer not found; running NormalizeData().")
    object <- NormalizeData(object, assay = assay_name, verbose = FALSE)
  }

  object
}

get_normalized_matrix <- function(object, assay = DefaultAssay(object)) {
  mat <- tryCatch(
    GetAssayData(object, assay = assay, slot = "data"),
    error = function(e) NULL
  )

  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    stop("Could not retrieve normalized expression matrix from assay: ", assay)
  }
  mat
}

match_genes_case_insensitive <- function(requested_genes, available_genes) {
  lookup <- setNames(available_genes, toupper(available_genes))
  matched <- unname(lookup[toupper(requested_genes)])
  unique(stats::na.omit(matched))
}

module_mean <- function(expr, genes) {
  genes_present <- match_genes_case_insensitive(genes, rownames(expr))
  if (length(genes_present) == 0L) {
    return(rep(0, ncol(expr)))
  }
  Matrix::colMeans(expr[genes_present, , drop = FALSE])
}

module_detected_count <- function(expr, genes) {
  genes_present <- match_genes_case_insensitive(genes, rownames(expr))
  if (length(genes_present) == 0L) {
    return(rep(0L, ncol(expr)))
  }
  as.integer(Matrix::colSums(expr[genes_present, , drop = FALSE] > 0))
}

any_gene_detected <- function(expr, genes) {
  module_detected_count(expr, genes) >= 1L
}

zscore_safe <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

identify_lower_left_display_island <- function(
    embedding,
    mode = c("auto", "manual", "none"),
    manual_x_max = -Inf,
    manual_y_max = -Inf,
    k_centers = 10L,
    max_fraction = 0.10,
    min_cells = 20L,
    min_separation_z = 1.10) {

  mode <- match.arg(mode)
  n <- nrow(embedding)
  excluded <- rep(FALSE, n)
  names(excluded) <- rownames(embedding)

  if (mode == "none") {
    return(list(excluded = excluded, diagnostics = NULL, reason = "none"))
  }

  x <- embedding[, 1]
  y <- embedding[, 2]

  if (mode == "manual") {
    excluded <- x <= manual_x_max & y <= manual_y_max
    names(excluded) <- rownames(embedding)
    return(list(
      excluded = excluded,
      diagnostics = NULL,
      reason = paste0("manual: UMAP1 <= ", manual_x_max, " and UMAP2 <= ", manual_y_max)
    ))
  }

  # Automatic conservative identification:
  # choose the small k-means component with the lowest standardized centroid sum,
  # but only if it is sufficiently separated from all other centroids.
  if (n < max(50L, k_centers * 5L)) {
    return(list(excluded = excluded, diagnostics = NULL, reason = "auto skipped: too few cells"))
  }

  z <- cbind(zscore_safe(x), zscore_safe(y))
  colnames(z) <- c("UMAP1_z", "UMAP2_z")
  k_use <- max(2L, min(as.integer(k_centers), floor(n / 20L)))

  km <- stats::kmeans(z, centers = k_use, nstart = 50, iter.max = 200)
  cluster_id <- km$cluster

  cluster_table <- tibble(
    kmeans_cluster = seq_len(k_use),
    n_cells = as.integer(tabulate(cluster_id, nbins = k_use)),
    fraction = n_cells / n,
    centroid_x_z = km$centers[, 1],
    centroid_y_z = km$centers[, 2],
    lower_left_score = centroid_x_z + centroid_y_z
  ) %>%
    arrange(lower_left_score)

  candidate <- cluster_table$kmeans_cluster[[1]]
  candidate_n <- cluster_table$n_cells[cluster_table$kmeans_cluster == candidate]
  candidate_fraction <- candidate_n / n

  candidate_center <- km$centers[candidate, , drop = FALSE]
  other_centers <- km$centers[-candidate, , drop = FALSE]
  centroid_distances <- sqrt(rowSums((other_centers - matrix(
    candidate_center,
    nrow = nrow(other_centers),
    ncol = 2,
    byrow = TRUE
  ))^2))
  min_separation <- min(centroid_distances)

  accepted <- candidate_n >= min_cells &&
    candidate_fraction <= max_fraction &&
    is.finite(min_separation) &&
    min_separation >= min_separation_z

  if (accepted) {
    excluded <- cluster_id == candidate
  }
  names(excluded) <- rownames(embedding)

  diagnostics <- cluster_table %>%
    mutate(
      selected_candidate = kmeans_cluster == candidate,
      accepted_for_display_exclusion = selected_candidate & accepted,
      minimum_centroid_separation_z = if_else(selected_candidate, min_separation, NA_real_)
    )

  reason <- if (accepted) {
    paste0(
      "auto accepted k-means component ", candidate,
      " (n=", candidate_n,
      ", fraction=", round(candidate_fraction, 4),
      ", separation_z=", round(min_separation, 3), ")"
    )
  } else {
    paste0(
      "auto candidate rejected; no display exclusion applied",
      " (candidate n=", candidate_n,
      ", fraction=", round(candidate_fraction, 4),
      ", separation_z=", round(min_separation, 3), ")"
    )
  }

  list(excluded = excluded, diagnostics = diagnostics, reason = reason)
}

make_umap_plot <- function(
    plot_df,
    color_col,
    palette,
    title,
    split_col = NULL,
    point_size = 0.45,
    alpha = 0.90) {

  p <- ggplot(
    plot_df,
    aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])
  ) +
    geom_point(size = point_size, alpha = alpha, stroke = 0, raster = FALSE) +
    scale_color_manual(values = palette, drop = FALSE, na.value = "#B8B8B8") +
    coord_equal() +
    labs(title = title, color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.key.height = grid::unit(0.48, "cm"),
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.4, alpha = 1)))

  if (!is.null(split_col)) {
    p <- p + facet_wrap(vars(.data[[split_col]]), nrow = 1)
  }

  p
}

# ------------------------------------------------------------------------------
# 3. Load and validate input
# ------------------------------------------------------------------------------
message_time("Reading input RDS: ", INPUT_RDS)
if (!file.exists(INPUT_RDS)) {
  stop("Input RDS does not exist: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) {
  stop("Input object is not a Seurat object.")
}

validate_metadata(
  obj,
  c(LAYER1_SOURCE_COL, CONDITION_SOURCE_COL, SAMPLE_SOURCE_COL)
)

message_time("Input cells: ", ncol(obj))

# ------------------------------------------------------------------------------
# 4. Extract STD/CDAHFD Kupffer_Macrophage + Monocyte
# ------------------------------------------------------------------------------
meta0 <- obj@meta.data
keep_cells <- rownames(meta0)[
  as.character(meta0[[CONDITION_SOURCE_COL]]) %in% TARGET_CONDITIONS &
    as.character(meta0[[LAYER1_SOURCE_COL]]) %in% TARGET_LAYER1
]

if (length(keep_cells) == 0L) {
  stop(
    "No cells matched conditions [", paste(TARGET_CONDITIONS, collapse = ", "),
    "] and Layer1 [", paste(TARGET_LAYER1, collapse = ", "), "]."
  )
}

mphi <- subset(obj, cells = keep_cells)
rm(obj)
gc(verbose = FALSE)

mphi$condition_internal <- factor(
  as.character(mphi@meta.data[[CONDITION_SOURCE_COL]]),
  levels = TARGET_CONDITIONS
)
mphi$sample_internal <- as.character(mphi@meta.data[[SAMPLE_SOURCE_COL]])

cluster_source_col <- choose_cluster_column(mphi, CLUSTER_COLUMN_CANDIDATES)
if (!is.null(cluster_source_col)) {
  mphi$cluster_internal <- as.character(mphi@meta.data[[cluster_source_col]])
  message_time("cluster_internal copied from: ", cluster_source_col)
} else {
  mphi$cluster_internal <- as.character(Idents(mphi))
  message_time("cluster_internal copied from active identities.")
}

message_time("Extracted macrophage/monocyte cells: ", ncol(mphi))
message_time(
  "Condition counts: ",
  paste(names(table(mphi$condition_internal)), table(mphi$condition_internal), collapse = "; ")
)

# ------------------------------------------------------------------------------
# 5. Prepare normalized expression
# ------------------------------------------------------------------------------
mphi <- prepare_rna_assay(mphi, PREFERRED_ASSAY)
assay_used <- DefaultAssay(mphi)
expr <- get_normalized_matrix(mphi, assay = assay_used)

# Ensure matrix columns follow object cell order.
expr <- expr[, colnames(mphi), drop = FALSE]

# ------------------------------------------------------------------------------
# 6. Marker definitions
# ------------------------------------------------------------------------------
# Resident markers emphasize bona fide Kupffer identity rather than generic Mphi.
resident_core <- c("Clec4f", "Timd4", "Vsig4", "Marco")
resident_extended <- c(
  "Clec4f", "Timd4", "Vsig4", "Marco", "Cd5l", "Adgre1",
  "C1qa", "C1qb", "C1qc", "Fcgrt", "Slc40a1", "Hmox1"
)

# Monocyte definition emphasizes circulating/recently recruited monocyte markers.
monocyte_core <- c("Ly6c2", "Ccr2", "S100a8", "S100a9", "Plac8")
monocyte_extended <- c(
  "Ly6c2", "Ccr2", "S100a8", "S100a9", "Plac8", "Chil3",
  "Sell", "Lair1", "Ctss", "Lyz2", "Ms4a7", "Fcgr1"
)

# Strict M1 markers: inflammatory effector program. Immediate-early genes alone
# are not sufficient. Strict M1 is allowed only after excluding Monocyte-like cells.
m1_markers <- c(
  "Il1b", "Tnf", "Cxcl2", "Ccl3", "Ccl4", "Ptgs2",
  "Cd80", "Cd86", "Nos2", "Il6", "Nfkbia", "Ier3"
)
m1_anchor <- c("Il1b", "Tnf", "Ptgs2", "Cd80", "Cd86", "Nos2", "Il6")

# Pro-resolution/M2 markers emphasize anti-inflammatory, reparative, and
# efferocytosis-associated macrophage states.
m2_markers <- c(
  "Mrc1", "Cd163", "Folr2", "Il10", "Retnla", "Arg1",
  "Mertk", "Axl", "Maf", "Klf4", "C1qa", "C1qb", "C1qc"
)
m2_anchor <- c("Mrc1", "Cd163", "Folr2", "Il10", "Retnla", "Arg1", "Mertk")

# MASH-associated lipid/scar macrophage program.
spp1_trem2_markers <- c(
  "Spp1", "Trem2", "Gpnmb", "Lpl", "Cd9", "Lgals3",
  "Apoe", "Fabp5", "Ctsb", "Ctsd", "Ctsk", "Tyrobp"
)
spp1_trem2_anchor <- c("Spp1", "Trem2", "Gpnmb", "Cd9")

# Mature macrophage support markers used to prevent a monocyte-only inflammatory
# program from being called strict M1.
mature_macrophage_support <- c(
  "Adgre1", "Csf1r", "Fcgr1", "Apoe", "C1qa", "C1qb", "C1qc", "Lgals3"
)

# ------------------------------------------------------------------------------
# 7. Calculate cell-level scores and strict gates
# ------------------------------------------------------------------------------
score_df <- tibble(
  cell = colnames(mphi),
  resident_score = module_mean(expr, resident_extended),
  monocyte_score = module_mean(expr, monocyte_extended),
  m1_score = module_mean(expr, m1_markers),
  m2_score = module_mean(expr, m2_markers),
  spp1_trem2_score = module_mean(expr, spp1_trem2_markers),
  resident_n = module_detected_count(expr, resident_extended),
  resident_core_n = module_detected_count(expr, resident_core),
  monocyte_n = module_detected_count(expr, monocyte_extended),
  monocyte_core_n = module_detected_count(expr, monocyte_core),
  m1_n = module_detected_count(expr, m1_markers),
  m1_anchor_n = module_detected_count(expr, m1_anchor),
  m2_n = module_detected_count(expr, m2_markers),
  m2_anchor_n = module_detected_count(expr, m2_anchor),
  spp1_trem2_n = module_detected_count(expr, spp1_trem2_markers),
  spp1_trem2_anchor_n = module_detected_count(expr, spp1_trem2_anchor),
  mature_mphi_n = module_detected_count(expr, mature_macrophage_support),
  spp1_or_trem2_detected = any_gene_detected(expr, c("Spp1", "Trem2"))
)

# Standardize module means before cross-program comparisons.
score_df <- score_df %>%
  mutate(
    resident_score_z = zscore_safe(resident_score),
    monocyte_score_z = zscore_safe(monocyte_score),
    m1_score_z = zscore_safe(m1_score),
    m2_score_z = zscore_safe(m2_score),
    spp1_trem2_score_z = zscore_safe(spp1_trem2_score)
  )

# Clear Resident vs Monocyte separation.
score_df <- score_df %>%
  mutate(
    resident_strict =
      resident_n >= MIN_RESIDENT_MARKERS &
      resident_core_n >= 1L &
      monocyte_core_n <= 1L &
      resident_score_z >= monocyte_score_z + SCORE_MARGIN_RESIDENT_MONOCYTE,

    monocyte_strict =
      monocyte_core_n >= MIN_MONOCYTE_CORE_MARKERS &
      resident_core_n == 0L &
      monocyte_score_z >= resident_score_z + SCORE_MARGIN_RESIDENT_MONOCYTE
  )

# Strict M1 intentionally excludes cells meeting the Monocyte gate and requires:
#   - multiple inflammatory genes,
#   - at least one inflammatory anchor,
#   - evidence of macrophage maturation,
#   - dominance over the M2 program.
score_df <- score_df %>%
  mutate(
    m1_strict =
      !resident_strict &
      !monocyte_strict &
      monocyte_core_n <= 1L &
      m1_n >= MIN_M1_MARKERS &
      m1_anchor_n >= 1L &
      mature_mphi_n >= 2L &
      m1_score_z >= m2_score_z + SCORE_MARGIN_M1_M2,

    m2_strict =
      !resident_strict &
      !monocyte_strict &
      m2_n >= MIN_M2_MARKERS &
      m2_anchor_n >= 1L &
      mature_mphi_n >= 2L &
      m2_score_z >= m1_score_z + SCORE_MARGIN_M1_M2,

    spp1_trem2_strict =
      !resident_strict &
      !monocyte_strict &
      !m1_strict &
      !m2_strict &
      spp1_or_trem2_detected &
      spp1_trem2_anchor_n >= 1L &
      spp1_trem2_n >= MIN_SPP1_TREM2_MARKERS
  )

# Annotation precedence is explicit and non-overlapping.
score_df <- score_df %>%
  mutate(
    layer2_internal = case_when(
      resident_strict ~ "Resident Kupffer-like",
      monocyte_strict ~ "Monocyte-like",
      m1_strict ~ "Inflammatory M1-like",
      m2_strict ~ "Pro-resolution M2-like",
      spp1_trem2_strict ~ "SPP1/TREM2 MASH-associated",
      TRUE ~ "Other"
    ),
    layer2_internal = factor(layer2_internal, levels = LAYER2_LEVELS)
  )

# Add scores and labels to the Seurat object in exact cell order.
score_df_ordered <- score_df[match(colnames(mphi), score_df$cell), ]
stopifnot(identical(score_df_ordered$cell, colnames(mphi)))

metadata_to_add <- score_df_ordered %>%
  select(-cell) %>%
  as.data.frame()
rownames(metadata_to_add) <- colnames(mphi)

mphi <- AddMetaData(mphi, metadata = metadata_to_add)
mphi$layer2_internal <- factor(mphi$layer2_internal, levels = LAYER2_LEVELS)

# ------------------------------------------------------------------------------
# 8. UMAP display-only exclusion
# ------------------------------------------------------------------------------
umap_reduction <- choose_umap_reduction(mphi, PREFERRED_UMAP_REDUCTION)
umap_embedding <- Embeddings(mphi, reduction = umap_reduction)[, 1:2, drop = FALSE]
colnames(umap_embedding) <- c("UMAP_1", "UMAP_2")

island_result <- identify_lower_left_display_island(
  embedding = umap_embedding,
  mode = DISPLAY_EXCLUDE_MODE,
  manual_x_max = MANUAL_UMAP1_MAX,
  manual_y_max = MANUAL_UMAP2_MAX,
  k_centers = AUTO_KMEANS_CENTERS,
  max_fraction = AUTO_MAX_EXCLUDED_FRACTION,
  min_cells = AUTO_MIN_EXCLUDED_CELLS,
  min_separation_z = AUTO_MIN_CENTROID_SEPARATION_Z
)

mphi$display_excluded_left_bottom <- island_result$excluded[colnames(mphi)]
message_time("Display exclusion: ", island_result$reason)
message_time(
  "Cells excluded from UMAP display only: ",
  sum(mphi$display_excluded_left_bottom),
  " / ", ncol(mphi)
)

if (!is.null(island_result$diagnostics)) {
  write.csv(
    island_result$diagnostics,
    file.path(OUTPUT_DIR, "UMAP_display_exclusion_auto_diagnostics.csv"),
    row.names = FALSE
  )
}

plot_df_all <- bind_cols(
  tibble(cell = rownames(umap_embedding)),
  as_tibble(umap_embedding),
  mphi@meta.data[rownames(umap_embedding), , drop = FALSE] %>%
    rownames_to_column("cell_meta")
) %>%
  mutate(
    condition_internal = factor(condition_internal, levels = TARGET_CONDITIONS),
    layer2_internal = factor(layer2_internal, levels = LAYER2_LEVELS),
    display_status = if_else(
      display_excluded_left_bottom,
      "Display-excluded lower-left island",
      "Displayed"
    )
  )

plot_df <- plot_df_all %>% filter(!display_excluded_left_bottom)

# Diagnostic plot always shows all cells and highlights any display exclusion.
diagnostic_colors <- c(
  "Displayed" = "#B8B8B8",
  "Display-excluded lower-left island" = "#FF1A1A"
)

p_exclusion <- make_umap_plot(
  plot_df_all,
  color_col = "display_status",
  palette = diagnostic_colors,
  title = paste0("UMAP display exclusion diagnostic\n", island_result$reason),
  point_size = 0.48
)
save_pdf("01_UMAP_display_exclusion_diagnostic.pdf", p_exclusion, 8.5, 7.0)

# ------------------------------------------------------------------------------
# 9. UMAP outputs
# ------------------------------------------------------------------------------
p_layer2 <- make_umap_plot(
  plot_df,
  color_col = "layer2_internal",
  palette = LAYER2_COLORS,
  title = "Mouse macrophage Layer2 annotation v4.0",
  point_size = 0.52
)
save_pdf("02_UMAP_Layer2_v4.0.pdf", p_layer2, 9.2, 7.2)

p_layer2_split <- make_umap_plot(
  plot_df,
  color_col = "layer2_internal",
  palette = LAYER2_COLORS,
  title = "Mouse macrophage Layer2 annotation: STD vs CDAHFD",
  split_col = "condition_internal",
  point_size = 0.48
)
save_pdf("03_UMAP_Layer2_split_STD_CDAHFD_v4.0.pdf", p_layer2_split, 13.0, 6.3)

p_condition <- make_umap_plot(
  plot_df,
  color_col = "condition_internal",
  palette = CONDITION_COLORS,
  title = "Mouse macrophages: STD vs CDAHFD",
  point_size = 0.50
)
save_pdf("04_UMAP_condition_STD_CDAHFD_v4.0.pdf", p_condition, 8.2, 6.8)

# ------------------------------------------------------------------------------
# 10. Quantitative summaries: ALL cells retained
# ------------------------------------------------------------------------------
meta_out <- mphi@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    condition_internal = factor(condition_internal, levels = TARGET_CONDITIONS),
    layer2_internal = factor(layer2_internal, levels = LAYER2_LEVELS)
  )

# Cell-level audit table
cell_audit <- meta_out %>%
  select(
    cell,
    condition_internal,
    sample_internal,
    cluster_internal,
    layer2_internal,
    display_excluded_left_bottom,
    resident_strict,
    monocyte_strict,
    m1_strict,
    m2_strict,
    spp1_trem2_strict,
    resident_score,
    monocyte_score,
    m1_score,
    m2_score,
    spp1_trem2_score,
    resident_n,
    resident_core_n,
    monocyte_n,
    monocyte_core_n,
    m1_n,
    m1_anchor_n,
    m2_n,
    m2_anchor_n,
    spp1_trem2_n,
    spp1_trem2_anchor_n,
    mature_mphi_n
  )

write.csv(
  cell_audit,
  file.path(OUTPUT_DIR, "Cell_level_Layer2_annotation_audit_v4.0.csv"),
  row.names = FALSE
)

# Layer2 counts and proportions by sample
layer2_by_sample <- meta_out %>%
  count(condition_internal, sample_internal, layer2_internal, name = "n_cells", .drop = FALSE) %>%
  group_by(condition_internal, sample_internal) %>%
  mutate(
    total_mphi_monocyte = sum(n_cells),
    proportion = if_else(total_mphi_monocyte > 0, n_cells / total_mphi_monocyte, NA_real_)
  ) %>%
  ungroup()

write.csv(
  layer2_by_sample,
  file.path(OUTPUT_DIR, "Layer2_counts_proportions_by_sample_v4.0.csv"),
  row.names = FALSE
)

# Layer2 counts and proportions by condition
layer2_by_condition <- meta_out %>%
  count(condition_internal, layer2_internal, name = "n_cells", .drop = FALSE) %>%
  group_by(condition_internal) %>%
  mutate(
    total_mphi_monocyte = sum(n_cells),
    proportion = if_else(total_mphi_monocyte > 0, n_cells / total_mphi_monocyte, NA_real_)
  ) %>%
  ungroup()

write.csv(
  layer2_by_condition,
  file.path(OUTPUT_DIR, "Layer2_counts_proportions_by_condition_v4.0.csv"),
  row.names = FALSE
)

# Strict M1/M2 denominator is STRICT M1 + STRICT M2 only.
strict_m1_m2_by_sample <- meta_out %>%
  group_by(condition_internal, sample_internal) %>%
  summarise(
    strict_M1_n = sum(layer2_internal == "Inflammatory M1-like", na.rm = TRUE),
    strict_M2_n = sum(layer2_internal == "Pro-resolution M2-like", na.rm = TRUE),
    strict_M1_plus_M2_n = strict_M1_n + strict_M2_n,
    strict_M1_fraction_within_M1_M2 = if_else(
      strict_M1_plus_M2_n > 0,
      strict_M1_n / strict_M1_plus_M2_n,
      NA_real_
    ),
    strict_M2_fraction_within_M1_M2 = if_else(
      strict_M1_plus_M2_n > 0,
      strict_M2_n / strict_M1_plus_M2_n,
      NA_real_
    ),
    strict_M2_to_M1_ratio = if_else(
      strict_M1_n > 0,
      strict_M2_n / strict_M1_n,
      NA_real_
    ),
    log2_strict_M2_to_M1_ratio_pseudocount_0.5 = log2(
      (strict_M2_n + 0.5) / (strict_M1_n + 0.5)
    ),
    .groups = "drop"
  )

write.csv(
  strict_m1_m2_by_sample,
  file.path(OUTPUT_DIR, "Strict_M1_M2_ratio_by_sample_v4.0.csv"),
  row.names = FALSE
)

strict_m1_m2_by_condition <- meta_out %>%
  group_by(condition_internal) %>%
  summarise(
    strict_M1_n = sum(layer2_internal == "Inflammatory M1-like", na.rm = TRUE),
    strict_M2_n = sum(layer2_internal == "Pro-resolution M2-like", na.rm = TRUE),
    strict_M1_plus_M2_n = strict_M1_n + strict_M2_n,
    strict_M1_fraction_within_M1_M2 = if_else(
      strict_M1_plus_M2_n > 0,
      strict_M1_n / strict_M1_plus_M2_n,
      NA_real_
    ),
    strict_M2_fraction_within_M1_M2 = if_else(
      strict_M1_plus_M2_n > 0,
      strict_M2_n / strict_M1_plus_M2_n,
      NA_real_
    ),
    strict_M2_to_M1_ratio = if_else(
      strict_M1_n > 0,
      strict_M2_n / strict_M1_n,
      NA_real_
    ),
    log2_strict_M2_to_M1_ratio_pseudocount_0.5 = log2(
      (strict_M2_n + 0.5) / (strict_M1_n + 0.5)
    ),
    .groups = "drop"
  )

write.csv(
  strict_m1_m2_by_condition,
  file.path(OUTPUT_DIR, "Strict_M1_M2_ratio_by_condition_v4.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Composition plots
# ------------------------------------------------------------------------------
p_composition_sample <- ggplot(
  layer2_by_sample,
  aes(x = sample_internal, y = proportion, fill = layer2_internal)
) +
  geom_col(width = 0.78) +
  facet_grid(. ~ condition_internal, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = LAYER2_COLORS, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    x = NULL,
    y = "Fraction of Kupffer_Macrophage + Monocyte cells",
    fill = NULL,
    title = "Layer2 composition by sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

save_pdf("05_Layer2_composition_by_sample_v4.0.pdf", p_composition_sample, 11.5, 6.5)

strict_long <- strict_m1_m2_by_sample %>%
  select(
    condition_internal,
    sample_internal,
    strict_M1_fraction_within_M1_M2,
    strict_M2_fraction_within_M1_M2
  ) %>%
  pivot_longer(
    cols = starts_with("strict_"),
    names_to = "strict_class",
    values_to = "fraction"
  ) %>%
  mutate(
    strict_class = recode(
      strict_class,
      "strict_M1_fraction_within_M1_M2" = "Strict M1",
      "strict_M2_fraction_within_M1_M2" = "Strict M2"
    ),
    strict_class = factor(strict_class, levels = c("Strict M1", "Strict M2"))
  )

strict_colors <- c("Strict M1" = LAYER2_COLORS[["Inflammatory M1-like"]],
                   "Strict M2" = LAYER2_COLORS[["Pro-resolution M2-like"]])

p_strict_m1_m2 <- ggplot(
  strict_long,
  aes(x = sample_internal, y = fraction, fill = strict_class)
) +
  geom_col(width = 0.72) +
  facet_grid(. ~ condition_internal, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = strict_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Fraction within strict M1 + strict M2",
    fill = NULL,
    title = "Strict M1/M2 composition\nDenominator = strict M1 + strict M2 only"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

save_pdf("06_Strict_M1_M2_fraction_by_sample_v4.0.pdf", p_strict_m1_m2, 10.5, 6.5)

# ------------------------------------------------------------------------------
# 12. Cluster-level audit
# ------------------------------------------------------------------------------
cluster_audit <- meta_out %>%
  count(
    condition_internal,
    cluster_internal,
    layer2_internal,
    name = "n_cells",
    .drop = FALSE
  ) %>%
  group_by(condition_internal, cluster_internal) %>%
  mutate(
    cluster_total = sum(n_cells),
    fraction_within_cluster = if_else(cluster_total > 0, n_cells / cluster_total, NA_real_)
  ) %>%
  ungroup()

write.csv(
  cluster_audit,
  file.path(OUTPUT_DIR, "Layer2_by_cluster_and_condition_audit_v4.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Save annotated Seurat object and run manifest
# ------------------------------------------------------------------------------
output_rds <- file.path(
  OUTPUT_DIR,
  "Mouse_STD_CDAHFD_Mphi_Layer2_annotated_v4.0.0.rds"
)
saveRDS(mphi, output_rds, compress = FALSE)

marker_manifest <- tibble(
  marker_set = c(
    "resident_core",
    "resident_extended",
    "monocyte_core",
    "monocyte_extended",
    "m1_markers",
    "m1_anchor",
    "m2_markers",
    "m2_anchor",
    "spp1_trem2_markers",
    "spp1_trem2_anchor",
    "mature_macrophage_support"
  ),
  genes = c(
    paste(resident_core, collapse = ";"),
    paste(resident_extended, collapse = ";"),
    paste(monocyte_core, collapse = ";"),
    paste(monocyte_extended, collapse = ";"),
    paste(m1_markers, collapse = ";"),
    paste(m1_anchor, collapse = ";"),
    paste(m2_markers, collapse = ";"),
    paste(m2_anchor, collapse = ";"),
    paste(spp1_trem2_markers, collapse = ";"),
    paste(spp1_trem2_anchor, collapse = ";"),
    paste(mature_macrophage_support, collapse = ";")
  )
)
write.csv(
  marker_manifest,
  file.path(OUTPUT_DIR, "Layer2_marker_manifest_v4.0.csv"),
  row.names = FALSE
)

run_manifest <- c(
  paste0("script_version: 4.0.0"),
  paste0("run_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("input_rds: ", INPUT_RDS),
  paste0("output_dir: ", OUTPUT_DIR),
  paste0("assay_used: ", assay_used),
  paste0("umap_reduction_used: ", umap_reduction),
  paste0("cluster_source: ", ifelse(is.null(cluster_source_col), "active_ident", cluster_source_col)),
  paste0("n_cells_extracted: ", ncol(mphi)),
  paste0("display_exclusion_mode: ", DISPLAY_EXCLUDE_MODE),
  paste0("display_exclusion_reason: ", island_result$reason),
  paste0("n_cells_display_excluded_only: ", sum(mphi$display_excluded_left_bottom)),
  paste0("strict_M1_M2_denominator: strict M1 + strict M2 only"),
  "sessionInfo:",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(run_manifest, file.path(OUTPUT_DIR, "run_manifest_v4.0.txt"))

message_time("Completed v4.0.0 analysis.")
message_time("Annotated RDS: ", output_rds)
message_time("All outputs: ", OUTPUT_DIR)
