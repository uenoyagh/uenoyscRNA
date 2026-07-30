#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(openxlsx)
  library(patchwork)
})

# ============================================================
# RDS3 annotation comparison
#
# Compares:
#   1. Existing annotation
#   2. General-marker classification
#   3. Ueno-marker classification
#
# Independent of FindAllMarkers().
# Uses RNA expression and AddModuleScore() directly on RDS3.
#
# Outputs:
#   - Annotation UMAPs
#   - Cluster-level prediction tables
#   - Cell-level confidence/margin summaries
#   - Marker violin plot collections
#   - Comparison and recommended annotation tables
#   - RDS containing all comparison metadata
# ============================================================

set.seed(260730)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS3_validation/",
  "Phase3_AnnotationComparison"
)

dir_umap <- file.path(output_dir, "UMAP")
dir_violin <- file.path(output_dir, "Violin")
dir_tables <- file.path(output_dir, "Tables")
dir_rds <- file.path(output_dir, "RDS")
dir_logs <- file.path(output_dir, "Logs")

for (x in c(output_dir, dir_umap, dir_violin, dir_tables, dir_rds, dir_logs)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

umap_reduction_candidates <- c("umap", "integratedRPCA_umap")
cluster_col_candidates <- c(
  "cluster_for_R8plot_FIXED2",
  "seurat_clusters"
)

existing_annotation_candidates <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_for_R8plot",
  "celltype_auto_annotation"
)

sample_col_candidates <- c(
  "sample_for_R8plot_FIXED2",
  "sample_for_R8plot",
  "sample"
)

# Classification is assigned at cluster level from mean module scores.
# A cell-level top score and margin are also retained as confidence evidence.
minimum_available_markers <- 2
minimum_cluster_top_markers <- 2
high_margin <- 0.20
moderate_margin <- 0.08

# Violin plots are split into multiple PDF pages.
genes_per_violin_page <- 8

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

first_existing <- function(candidates, available, label) {
  hit <- candidates[candidates %in% available]

  if (length(hit) == 0) {
    stop(
      paste0(
        "No suitable ", label, " found. Candidates: ",
        paste(candidates, collapse = ", ")
      )
    )
  }

  hit[[1]]
}

safe_filename <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  str_replace_all(x, "^_|_$", "")
}

natural_levels <- function(x) {
  z <- unique(as.character(x))
  zn <- suppressWarnings(as.numeric(z))

  if (all(!is.na(zn))) {
    as.character(sort(unique(zn)))
  } else {
    sort(z)
  }
}

save_pdf <- function(plot, file, width, height) {
  ggsave(
    filename = file,
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    limitsize = FALSE
  )
}

write_excel <- function(file, sheets) {
  wb <- createWorkbook()

  for (nm in names(sheets)) {
    sheet_nm <- substr(
      str_replace_all(nm, "[\\\\/:*?\\[\\]]", "_"),
      1,
      31
    )

    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, sheets[[nm]], withFilter = TRUE)
    freezePane(wb, sheet_nm, firstRow = TRUE)

    if (ncol(sheets[[nm]]) > 0) {
      setColWidths(
        wb,
        sheet = sheet_nm,
        cols = seq_len(ncol(sheets[[nm]])),
        widths = "auto"
      )
    }
  }

  saveWorkbook(wb, file, overwrite = TRUE)
}

filter_marker_sets <- function(marker_sets, available_features) {
  out <- lapply(marker_sets, intersect, y = available_features)
  out[lengths(out) >= minimum_available_markers]
}

marker_set_table <- function(marker_sets, source_name) {
  bind_rows(
    lapply(names(marker_sets), function(label) {
      tibble(
        source = source_name,
        label = label,
        gene = marker_sets[[label]]
      )
    })
  )
}

# ------------------------------------------------------------
# Marker sets
# ------------------------------------------------------------

general_markers <- list(
  Hepatocyte = c("Alb", "Apoa1", "Ttr", "Fabp1", "Ass1", "Cps1"),
  Cholangiocyte = c("Krt7", "Krt8", "Krt18", "Krt19", "Epcam", "Sox9"),
  LSEC = c("Kdr", "Klf2", "Stab1", "Stab2", "Clec4g", "Pecam1"),
  Vascular_endothelial = c("Pecam1", "Cdh5", "Kdr", "Erg", "Emcn", "Vwf"),
  Kupffer_Macrophage = c(
    "Adgre1", "C1qa", "C1qb", "C1qc",
    "Cd68", "Clec4f", "Timd4", "Marco"
  ),
  Monocyte = c("Lyz2", "Ly6c2", "Ccr2", "S100a8", "S100a9", "Ctss"),
  Neutrophil = c("S100a8", "S100a9", "Ly6g", "Csf3r", "Mpo", "Elane"),
  Dendritic_cell = c("Flt3", "Itgax", "H2-Ab1", "Cd74", "Clec10a", "Xcr1"),
  T_cell = c("Cd3d", "Cd3e", "Trac", "Cd247", "Lck"),
  NK_cell = c("Nkg7", "Klrd1", "Klrk1", "Prf1", "Gzmb"),
  B_cell = c("Cd79a", "Cd79b", "Ms4a1", "Cd74", "H2-Ab1"),
  Plasma_cell = c("Jchain", "Mzb1", "Sdc1", "Xbp1", "Igha"),
  HSC_Mesenchymal = c(
    "Dcn", "Col1a1", "Col1a2", "Col3a1",
    "Lrat", "Rgs5", "Acta2"
  ),
  Mesothelial = c("Msln", "Krt19", "Upk3b", "Wt1", "Krt8"),
  Platelet = c("Pf4", "Ppbp", "Itga2b", "Gp9"),
  RBC = c("Hba-a1", "Hba-a2", "Hbb-bs", "Alas2"),
  Cycling = c("Mki67", "Top2a", "Pcna", "Cdk1", "Ube2c")
)

ueno_markers <- list(
  NK_cell = c("Klrd1", "Nkg7", "Gzmb", "Klrf1"),
  B_cell = c("Cd19", "Cd79a", "Ighm", "Ighd"),
  T_cell = c("Cd8a", "Cd3d", "Cd3e", "Trac", "Trbc1", "Trbc2"),
  T_NKT_candidate = c("Cd4", "Nkg7", "Klrd1"),
  Effector_T = c("Cd3d", "Cd3e", "Gzmb", "Nkg7"),
  Naive_T = c("Sell", "Il7r", "Ccr7", "Ltb"),
  Monocyte = c("S100a8", "S100a9", "Ccr2", "Lyz2"),
  Macrophage_general = c("Aif1", "Csf1r", "Adgre1", "Cd68"),
  Macrophage_resident = c("Timd4", "Siglec1", "Clec4f", "Marco"),
  Macrophage_M1_like = c("Itgax", "Il1b", "Tnf", "Cd80", "Cd86"),
  Macrophage_M2a_like = c("Spic", "Cd209a", "Il10", "Vsig4"),
  Macrophage_M2c_like = c("Marco", "C1qb", "Mertk", "Cd163"),
  Hepatocyte_common = c(
    "Hnf4a", "Alb", "Fgb", "Ttr",
    "Ahsg", "Atf5", "Fabp1"
  ),
  Hepatoblast_like = c("Prox1", "Cdh1", "Afp", "Dlk1"),
  Hepatocyte_zone1 = c("Gls2", "Pck1", "Cps1", "Hal", "Ass1"),
  Hepatocyte_zone3 = c("Glul", "Cyp1a2", "Cyp2e1", "Lect2"),
  Cholangiocyte = c("Krt7", "Cldn3", "Cldn4", "Spp1", "Krt19"),
  Endothelial_general = c("Pecam1", "Erg", "Cdh5", "Cldn5"),
  LSEC = c("Fcgr2b", "Lyve1", "Stab1", "Stab2", "Clec4g"),
  Mesenchymal_general = c("Dcn", "Bgn"),
  Stellate_general = c("Zeb2", "Lrat", "Rgs5"),
  Stellate_nonactivated = c("Lrat", "Lhx2", "Reln", "Rbp1"),
  Stellate_activated = c("Acta2", "Col1a1", "Col3a1", "Tagln"),
  STM1 = c("Ptch1", "Dcn", "Col1a1"),
  STM2 = c("Pcolce", "Col3a1", "Dcn"),
  Meso1 = c("Tbx18", "Wt1", "Msln"),
  Meso2 = c("Foxf1", "Col15a1", "Dcn")
)

# Coarse mapping is required because Ueno labels include subtypes.
ueno_to_general <- c(
  NK_cell = "NK_cell",
  B_cell = "B_cell",
  T_cell = "T_cell",
  T_NKT_candidate = "T_cell",
  Effector_T = "T_cell",
  Naive_T = "T_cell",
  Monocyte = "Monocyte",
  Macrophage_general = "Kupffer_Macrophage",
  Macrophage_resident = "Kupffer_Macrophage",
  Macrophage_M1_like = "Kupffer_Macrophage",
  Macrophage_M2a_like = "Kupffer_Macrophage",
  Macrophage_M2c_like = "Kupffer_Macrophage",
  Hepatocyte_common = "Hepatocyte",
  Hepatoblast_like = "Hepatocyte",
  Hepatocyte_zone1 = "Hepatocyte",
  Hepatocyte_zone3 = "Hepatocyte",
  Cholangiocyte = "Cholangiocyte",
  Endothelial_general = "Vascular_endothelial",
  LSEC = "LSEC",
  Mesenchymal_general = "HSC_Mesenchymal",
  Stellate_general = "HSC_Mesenchymal",
  Stellate_nonactivated = "HSC_Mesenchymal",
  Stellate_activated = "HSC_Mesenchymal",
  STM1 = "HSC_Mesenchymal",
  STM2 = "HSC_Mesenchymal",
  Meso1 = "Mesothelial",
  Meso2 = "Mesothelial"
)

# ------------------------------------------------------------
# Load RDS3
# ------------------------------------------------------------

cat("Loading RDS3...\n")
obj <- readRDS(rds_file)

metadata_cols <- colnames(obj@meta.data)
reduction_name <- first_existing(
  umap_reduction_candidates,
  Reductions(obj),
  "UMAP reduction"
)
cluster_col <- first_existing(
  cluster_col_candidates,
  metadata_cols,
  "cluster column"
)
existing_col <- first_existing(
  existing_annotation_candidates,
  metadata_cols,
  "existing annotation column"
)
sample_col <- first_existing(
  sample_col_candidates,
  metadata_cols,
  "sample column"
)

DefaultAssay(obj) <- "RNA"

obj <- tryCatch(
  JoinLayers(obj, assay = "RNA"),
  error = function(e) obj
)

rna_data_exists <- tryCatch(
  {
    x <- GetAssayData(obj, assay = "RNA", layer = "data")
    nrow(x) > 0 && ncol(x) > 0
  },
  error = function(e) FALSE
)

if (!rna_data_exists) {
  cat("RNA data layer absent. Running NormalizeData()...\n")
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
}

features_present <- rownames(obj[["RNA"]])
cluster_levels <- natural_levels(obj@meta.data[[cluster_col]])

obj$annotation_existing_compare <- as.character(
  obj@meta.data[[existing_col]]
)
obj$cluster_compare <- factor(
  as.character(obj@meta.data[[cluster_col]]),
  levels = cluster_levels
)
obj$sample_compare <- as.character(
  obj@meta.data[[sample_col]]
)

# ------------------------------------------------------------
# Marker availability
# ------------------------------------------------------------

general_available <- filter_marker_sets(
  general_markers,
  features_present
)
ueno_available <- filter_marker_sets(
  ueno_markers,
  features_present
)

availability <- bind_rows(
  marker_set_table(general_markers, "General"),
  marker_set_table(ueno_markers, "Ueno_Custom")
) %>%
  mutate(available_in_RDS = gene %in% features_present)

availability_summary <- availability %>%
  group_by(source, label) %>%
  summarise(
    requested_markers = n(),
    available_markers = sum(available_in_RDS),
    availability_fraction = mean(available_in_RDS),
    available_gene_list = paste(gene[available_in_RDS], collapse = "; "),
    missing_gene_list = paste(gene[!available_in_RDS], collapse = "; "),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Module scoring
# ------------------------------------------------------------

add_marker_scores <- function(object, marker_sets, prefix) {
  score_map <- tibble(
    label = names(marker_sets),
    score_column = NA_character_
  )

  for (i in seq_along(marker_sets)) {
    generated_prefix <- paste0(prefix, sprintf("%02d", i), "_")

    object <- AddModuleScore(
      object = object,
      features = list(marker_sets[[i]]),
      assay = "RNA",
      name = generated_prefix,
      seed = 260730 + i,
      search = FALSE
    )

    score_map$score_column[i] <- paste0(generated_prefix, "1")
  }

  list(object = object, score_map = score_map)
}

cat("Calculating general-marker scores...\n")
tmp <- add_marker_scores(
  obj,
  general_available,
  "GeneralCompare"
)
obj <- tmp$object
general_score_map <- tmp$score_map

cat("Calculating Ueno-marker scores...\n")
tmp <- add_marker_scores(
  obj,
  ueno_available,
  "UenoCompare"
)
obj <- tmp$object
ueno_score_map <- tmp$score_map

# ------------------------------------------------------------
# Classification engine
# ------------------------------------------------------------

classify_from_scores <- function(
  object,
  score_map,
  source_name,
  cluster_col_name
) {
  meta <- object@meta.data %>%
    rownames_to_column("cell_barcode") %>%
    mutate(
      cluster_compare_internal = as.character(
        .data[[cluster_col_name]]
      )
    )

  score_long <- meta %>%
    select(
      cell_barcode,
      cluster_compare_internal,
      all_of(score_map$score_column)
    ) %>%
    pivot_longer(
      cols = all_of(score_map$score_column),
      names_to = "score_column",
      values_to = "marker_score"
    ) %>%
    left_join(score_map, by = "score_column")

  cell_prediction <- score_long %>%
    group_by(cell_barcode) %>%
    arrange(desc(marker_score), .by_group = TRUE) %>%
    summarise(
      predicted_label_cell = label[1],
      top_score_cell = marker_score[1],
      second_label_cell = ifelse(n() >= 2, label[2], NA_character_),
      second_score_cell = ifelse(n() >= 2, marker_score[2], NA_real_),
      score_margin_cell = ifelse(
        n() >= 2,
        marker_score[1] - marker_score[2],
        NA_real_
      ),
      .groups = "drop"
    )

  cluster_scores <- score_long %>%
    group_by(cluster_compare_internal, label) %>%
    summarise(
      mean_score = mean(marker_score, na.rm = TRUE),
      median_score = median(marker_score, na.rm = TRUE),
      positive_fraction = mean(marker_score > 0, na.rm = TRUE),
      .groups = "drop"
    )

  cluster_prediction <- cluster_scores %>%
    group_by(cluster_compare_internal) %>%
    arrange(desc(mean_score), .by_group = TRUE) %>%
    summarise(
      predicted_label_cluster = label[1],
      top_score_cluster = mean_score[1],
      second_label_cluster = ifelse(n() >= 2, label[2], NA_character_),
      second_score_cluster = ifelse(n() >= 2, mean_score[2], NA_real_),
      score_margin_cluster = ifelse(
        n() >= 2,
        mean_score[1] - mean_score[2],
        NA_real_
      ),
      confidence_cluster = case_when(
        score_margin_cluster >= high_margin ~ "High",
        score_margin_cluster >= moderate_margin ~ "Moderate",
        TRUE ~ "Low"
      ),
      .groups = "drop"
    )

  cell_prediction <- cell_prediction %>%
    left_join(
      meta %>%
        select(cell_barcode, cluster_compare_internal),
      by = "cell_barcode"
    ) %>%
    left_join(
      cluster_prediction %>%
        select(
          cluster_compare_internal,
          predicted_label_cluster,
          score_margin_cluster,
          confidence_cluster
        ),
      by = "cluster_compare_internal"
    ) %>%
    mutate(source = source_name)

  cluster_cell_agreement <- cell_prediction %>%
    group_by(
      cluster_compare_internal,
      predicted_label_cluster,
      confidence_cluster
    ) %>%
    summarise(
      cells = n(),
      cell_cluster_agreement = mean(
        predicted_label_cell == predicted_label_cluster
      ),
      median_cell_margin = median(
        score_margin_cell,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  cluster_prediction <- cluster_prediction %>%
    left_join(
      cluster_cell_agreement,
      by = c(
        "cluster_compare_internal",
        "predicted_label_cluster",
        "confidence_cluster"
      )
    )

  list(
    cell_prediction = cell_prediction,
    cluster_scores = cluster_scores,
    cluster_prediction = cluster_prediction
  )
}

general_result <- classify_from_scores(
  obj,
  general_score_map,
  "General",
  cluster_col
)

ueno_result <- classify_from_scores(
  obj,
  ueno_score_map,
  "Ueno_Custom",
  cluster_col
)

# ------------------------------------------------------------
# Add predictions to metadata
# ------------------------------------------------------------

general_cell_meta <- general_result$cell_prediction %>%
  select(
    cell_barcode,
    general_prediction_cell = predicted_label_cell,
    general_score_margin_cell = score_margin_cell
  )

ueno_cell_meta <- ueno_result$cell_prediction %>%
  select(
    cell_barcode,
    ueno_prediction_cell = predicted_label_cell,
    ueno_score_margin_cell = score_margin_cell
  )

cluster_annotation <- tibble(
  cluster_compare_internal = cluster_levels
) %>%
  left_join(
    obj@meta.data %>%
      rownames_to_column("cell_barcode") %>%
      transmute(
        cluster_compare_internal = as.character(
          .data[[cluster_col]]
        ),
        existing_annotation = as.character(
          .data[[existing_col]]
        )
      ) %>%
      count(
        cluster_compare_internal,
        existing_annotation,
        name = "cells"
      ) %>%
      group_by(cluster_compare_internal) %>%
      slice_max(cells, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(
        cluster_compare_internal,
        existing_annotation
      ),
    by = "cluster_compare_internal"
  ) %>%
  left_join(
    general_result$cluster_prediction %>%
      transmute(
        cluster_compare_internal,
        general_annotation = predicted_label_cluster,
        general_margin = score_margin_cluster,
        general_confidence = confidence_cluster,
        general_cell_cluster_agreement = cell_cluster_agreement
      ),
    by = "cluster_compare_internal"
  ) %>%
  left_join(
    ueno_result$cluster_prediction %>%
      transmute(
        cluster_compare_internal,
        ueno_annotation = predicted_label_cluster,
        ueno_margin = score_margin_cluster,
        ueno_confidence = confidence_cluster,
        ueno_cell_cluster_agreement = cell_cluster_agreement
      ),
    by = "cluster_compare_internal"
  ) %>%
  mutate(
    ueno_generalized = unname(
      ueno_to_general[ueno_annotation]
    ),
    general_ueno_coarse_concordance =
      general_annotation == ueno_generalized,
    recommendation = case_when(
      general_ueno_coarse_concordance &
        general_confidence == "High" &
        ueno_confidence %in% c("High", "Moderate") ~
        "General_supported; Ueno subtype may be retained",
      general_ueno_coarse_concordance &
        general_confidence %in% c("High", "Moderate") ~
        "Concordant; manual marker review",
      !general_ueno_coarse_concordance &
        general_confidence == "High" &
        ueno_confidence == "Low" ~
        "Prefer general classification",
      !general_ueno_coarse_concordance &
        general_confidence == "Low" &
        ueno_confidence == "High" ~
        "Prefer Ueno classification after marker review",
      TRUE ~
        "Manual review required"
    )
  ) %>%
  arrange(match(cluster_compare_internal, cluster_levels))

comparison_meta <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  left_join(general_cell_meta, by = "cell_barcode") %>%
  left_join(ueno_cell_meta, by = "cell_barcode") %>%
  mutate(
    cluster_compare_internal = as.character(
      .data[[cluster_col]]
    )
  ) %>%
  left_join(
    cluster_annotation %>%
      select(
        cluster_compare_internal,
        general_annotation,
        ueno_annotation,
        recommendation
      ),
    by = "cluster_compare_internal"
  )

rownames(comparison_meta) <- comparison_meta$cell_barcode
comparison_meta$cell_barcode <- NULL

obj@meta.data <- comparison_meta
obj$annotation_general_compare <- obj$general_annotation
obj$annotation_ueno_compare <- obj$ueno_annotation
obj$annotation_recommendation_compare <- obj$recommendation

# ------------------------------------------------------------
# Annotation UMAPs
# ------------------------------------------------------------

umap_existing <- DimPlot(
  obj,
  reduction = reduction_name,
  group.by = "annotation_existing_compare",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.20,
  raster = TRUE
) +
  ggtitle("RDS3: existing annotation") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 8))

umap_general <- DimPlot(
  obj,
  reduction = reduction_name,
  group.by = "annotation_general_compare",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.20,
  raster = TRUE
) +
  ggtitle("RDS3: general-marker classification") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 8))

umap_ueno <- DimPlot(
  obj,
  reduction = reduction_name,
  group.by = "annotation_ueno_compare",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.20,
  raster = TRUE
) +
  ggtitle("RDS3: Ueno-marker classification") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 8))

save_pdf(
  umap_existing,
  file.path(dir_umap, "01_UMAP_ExistingAnnotation.pdf"),
  12,
  9
)
save_pdf(
  umap_general,
  file.path(dir_umap, "02_UMAP_GeneralClassification.pdf"),
  12,
  9
)
save_pdf(
  umap_ueno,
  file.path(dir_umap, "03_UMAP_UenoClassification.pdf"),
  13,
  9
)

combined_umap <- (
  umap_existing |
  umap_general |
  umap_ueno
) +
  plot_layout(guides = "collect")

save_pdf(
  combined_umap,
  file.path(dir_umap, "04_UMAP_ThreeWayComparison.pdf"),
  30,
  9
)

# Split by sample
for (annotation_col in c(
  "annotation_existing_compare",
  "annotation_general_compare",
  "annotation_ueno_compare"
)) {
  p <- DimPlot(
    obj,
    reduction = reduction_name,
    group.by = annotation_col,
    split.by = "sample_compare",
    pt.size = 0.12,
    raster = TRUE,
    ncol = 3
  ) +
    theme_bw(base_size = 8) +
    theme(legend.text = element_text(size = 7))

  save_pdf(
    p,
    file.path(
      dir_umap,
      paste0(
        "SplitBySample_",
        safe_filename(annotation_col),
        ".pdf"
      )
    ),
    22,
    16
  )
}

# ------------------------------------------------------------
# Violin plots
# ------------------------------------------------------------

make_violin_pages <- function(
  object,
  marker_sets,
  group_column,
  prefix,
  title_prefix
) {
  for (set_name in names(marker_sets)) {
    genes <- marker_sets[[set_name]]
    genes <- intersect(genes, features_present)

    if (length(genes) == 0) {
      next
    }

    page_id <- ceiling(
      seq_along(genes) / genes_per_violin_page
    )

    for (page in unique(page_id)) {
      page_genes <- genes[page_id == page]

      p <- VlnPlot(
        object = object,
        features = page_genes,
        group.by = group_column,
        assay = "RNA",
        pt.size = 0,
        stack = FALSE,
        ncol = 2
      ) &
        theme_bw(base_size = 8) &
        theme(
          axis.text.x = element_text(
            angle = 60,
            hjust = 1,
            vjust = 1
          )
        )

      p <- p +
        plot_annotation(
          title = paste0(
            title_prefix,
            ": ",
            set_name,
            " (page ",
            page,
            ")"
          )
        )

      file_name <- paste0(
        prefix,
        "_",
        safe_filename(set_name),
        "_page",
        page,
        ".pdf"
      )

      save_pdf(
        p,
        file.path(dir_violin, file_name),
        15,
        max(8, 3.8 * ceiling(length(page_genes) / 2))
      )
    }
  }
}

# Existing annotation, evaluated using both marker libraries
make_violin_pages(
  obj,
  general_available,
  "annotation_existing_compare",
  "Existing_GeneralMarkers",
  "Existing annotation / general markers"
)

make_violin_pages(
  obj,
  ueno_available,
  "annotation_existing_compare",
  "Existing_UenoMarkers",
  "Existing annotation / Ueno markers"
)

# General classification, evaluated using general markers
make_violin_pages(
  obj,
  general_available,
  "annotation_general_compare",
  "GeneralClassification",
  "General classification"
)

# Ueno classification, evaluated using Ueno markers
make_violin_pages(
  obj,
  ueno_available,
  "annotation_ueno_compare",
  "UenoClassification",
  "Ueno classification"
)

# ------------------------------------------------------------
# Comparison tables
# ------------------------------------------------------------

existing_cluster_purity <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  transmute(
    cluster = as.character(.data[[cluster_col]]),
    annotation = annotation_existing_compare
  ) %>%
  count(cluster, annotation, name = "cells") %>%
  group_by(cluster) %>%
  mutate(
    cluster_cells = sum(cells),
    fraction = cells / cluster_cells
  ) %>%
  slice_max(fraction, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    existing_annotation = annotation,
    existing_cluster_purity = fraction
  )

cluster_comparison <- cluster_annotation %>%
  left_join(
    existing_cluster_purity %>%
      select(
        cluster,
        existing_cluster_purity
      ),
    by = c(
      "cluster_compare_internal" = "cluster"
    )
  ) %>%
  mutate(
    general_quality_score =
      general_margin *
      general_cell_cluster_agreement,
    ueno_quality_score =
      ueno_margin *
      ueno_cell_cluster_agreement,
    suggested_primary_annotation = case_when(
      general_confidence == "High" &
        general_cell_cluster_agreement >= 0.60 ~
        general_annotation,
      general_ueno_coarse_concordance &
        ueno_confidence %in% c("High", "Moderate") ~
        general_annotation,
      ueno_confidence == "High" &
        ueno_cell_cluster_agreement >= 0.60 ~
        ueno_annotation,
      TRUE ~
        existing_annotation
    ),
    suggested_subtype_annotation = case_when(
      general_ueno_coarse_concordance &
        ueno_confidence %in% c("High", "Moderate") ~
        ueno_annotation,
      TRUE ~
        NA_character_
    )
  )

annotation_counts <- bind_rows(
  obj@meta.data %>%
    count(annotation_existing_compare, name = "cells") %>%
    transmute(
      classification = "Existing",
      annotation = annotation_existing_compare,
      cells
    ),
  obj@meta.data %>%
    count(annotation_general_compare, name = "cells") %>%
    transmute(
      classification = "General",
      annotation = annotation_general_compare,
      cells
    ),
  obj@meta.data %>%
    count(annotation_ueno_compare, name = "cells") %>%
    transmute(
      classification = "Ueno_Custom",
      annotation = annotation_ueno_compare,
      cells
    )
)

write_csv(
  availability,
  file.path(dir_tables, "MarkerAvailability_All.csv")
)
write_csv(
  availability_summary,
  file.path(dir_tables, "MarkerAvailability_Summary.csv")
)
write_csv(
  general_result$cluster_scores,
  file.path(dir_tables, "General_ClusterScores.csv")
)
write_csv(
  general_result$cluster_prediction,
  file.path(dir_tables, "General_ClusterPrediction.csv")
)
write_csv(
  ueno_result$cluster_scores,
  file.path(dir_tables, "Ueno_ClusterScores.csv")
)
write_csv(
  ueno_result$cluster_prediction,
  file.path(dir_tables, "Ueno_ClusterPrediction.csv")
)
write_csv(
  cluster_comparison,
  file.path(dir_tables, "ThreeWay_ClusterComparison.csv")
)
write_csv(
  annotation_counts,
  file.path(dir_tables, "Annotation_CellCounts.csv")
)

write_excel(
  file.path(dir_tables, "RDS3_AnnotationComparison.xlsx"),
  list(
    ThreeWayComparison = cluster_comparison,
    GeneralPrediction = general_result$cluster_prediction,
    UenoPrediction = ueno_result$cluster_prediction,
    GeneralScores = general_result$cluster_scores,
    UenoScores = ueno_result$cluster_scores,
    AnnotationCounts = annotation_counts,
    MarkerAvailability = availability_summary
  )
)

# ------------------------------------------------------------
# Save annotated comparison RDS
# ------------------------------------------------------------

saveRDS(
  obj,
  file.path(
    dir_rds,
    "Mouse_RDS3_with_three_annotation_comparisons.rds"
  ),
  compress = FALSE
)

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------

settings <- tibble(
  item = c(
    "rds_file",
    "reduction",
    "cluster_column",
    "existing_annotation_column",
    "sample_column",
    "minimum_available_markers",
    "high_margin",
    "moderate_margin",
    "genes_per_violin_page"
  ),
  value = c(
    rds_file,
    reduction_name,
    cluster_col,
    existing_col,
    sample_col,
    minimum_available_markers,
    high_margin,
    moderate_margin,
    genes_per_violin_page
  )
)

write_csv(
  settings,
  file.path(dir_logs, "Settings.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(dir_logs, "sessionInfo.txt")
)

capture.output(
  warnings(),
  file = file.path(dir_logs, "warnings.txt")
)

cat("\n============================================================\n")
cat("RDS3 annotation comparison completed\n")
cat("============================================================\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main outputs:\n")
cat("  UMAP/01_UMAP_ExistingAnnotation.pdf\n")
cat("  UMAP/02_UMAP_GeneralClassification.pdf\n")
cat("  UMAP/03_UMAP_UenoClassification.pdf\n")
cat("  UMAP/04_UMAP_ThreeWayComparison.pdf\n")
cat("  Violin/*.pdf\n")
cat("  Tables/RDS3_AnnotationComparison.xlsx\n")
cat("  RDS/Mouse_RDS3_with_three_annotation_comparisons.rds\n")
