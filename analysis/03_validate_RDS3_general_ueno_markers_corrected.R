#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(openxlsx)
  library(stringr)
  library(scales)
  library(tibble)
})

# ============================================================
# RDS3 Phase 2-2 (corrected: tibble::rownames_to_column support)
# General markers + Ueno custom markers
# DotPlot, module scores, cluster-level predictions
# ============================================================

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS3_validation/",
  "Phase2_MarkerValidation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

general_dir <- file.path(output_dir, "General")
ueno_dir <- file.path(output_dir, "Ueno_Custom")
comparison_dir <- file.path(output_dir, "Comparison")
log_dir <- file.path(output_dir, "Logs")

for (x in c(general_dir, ueno_dir, comparison_dir, log_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

set.seed(260730)

# ============================================================
# Helpers
# ============================================================

first_existing <- function(candidates, columns, label) {
  hit <- candidates[candidates %in% columns]
  if (length(hit) == 0) {
    stop(
      paste0(
        "No metadata column found for ", label, ". Candidates: ",
        paste(candidates, collapse = ", ")
      )
    )
  }
  hit[[1]]
}

natural_cluster_levels <- function(x) {
  z <- unique(as.character(x))
  zn <- suppressWarnings(as.numeric(z))
  if (all(!is.na(zn))) {
    as.character(sort(unique(zn)))
  } else {
    sort(z)
  }
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

save_pdf <- function(plot, path, width, height) {
  ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
}

write_excel <- function(path, sheets) {
  wb <- createWorkbook()

  for (nm in names(sheets)) {
    sheet_name <- substr(sanitize_name(nm), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, sheets[[nm]], withFilter = TRUE)
    freezePane(wb, sheet_name, firstRow = TRUE)

    if (ncol(sheets[[nm]]) > 0) {
      setColWidths(
        wb,
        sheet = sheet_name,
        cols = seq_len(ncol(sheets[[nm]])),
        widths = "auto"
      )
    }
  }

  saveWorkbook(wb, path, overwrite = TRUE)
}

# ============================================================
# Load object
# ============================================================

cat("Loading RDS3...\n")
obj <- readRDS(rds_file)

cluster_col <- first_existing(
  c(
    "cluster_for_R8plot_FIXED2",
    "integratedRPCA_snn_res.0.8",
    "seurat_clusters"
  ),
  colnames(obj@meta.data),
  "cluster"
)

existing_celltype_col <- first_existing(
  c(
    "celltype_for_R8plot_FIXED2",
    "celltype_for_R8plot",
    "celltype_auto_annotation",
    "celltype"
  ),
  colnames(obj@meta.data),
  "existing cell type"
)

DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj@meta.data[[cluster_col]]

obj <- tryCatch(
  JoinLayers(obj, assay = "RNA"),
  error = function(e) obj
)

data_available <- tryCatch(
  {
    x <- GetAssayData(obj, assay = "RNA", layer = "data")
    nrow(x) > 0 && ncol(x) > 0
  },
  error = function(e) FALSE
)

if (!data_available) {
  cat("RNA data layer not available. Running NormalizeData()...\n")
  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
}

features_present <- rownames(obj[["RNA"]])

cluster_levels <- natural_cluster_levels(obj@meta.data[[cluster_col]])
Idents(obj) <- factor(
  as.character(obj@meta.data[[cluster_col]]),
  levels = cluster_levels
)

# ============================================================
# Marker sets
# ============================================================

general_markers <- list(
  Hepatocyte = c("Alb", "Apoa1", "Ttr", "Fabp1", "Ass1", "Cps1"),
  Cholangiocyte = c("Krt7", "Krt8", "Krt18", "Krt19", "Epcam", "Sox9"),
  LSEC = c("Kdr", "Klf2", "Stab1", "Stab2", "Clec4g", "Pecam1"),
  Vascular_endothelial = c("Pecam1", "Cdh5", "Kdr", "Erg", "Emcn", "Vwf"),
  Kupffer_Macrophage = c("Adgre1", "C1qa", "C1qb", "C1qc", "Cd68", "Clec4f", "Timd4", "Marco"),
  Monocyte = c("Lyz2", "Ly6c2", "Ccr2", "S100a8", "S100a9", "Ctss"),
  Neutrophil = c("S100a8", "S100a9", "Ly6g", "Csf3r", "Mpo", "Elane"),
  Dendritic_cell = c("Flt3", "Itgax", "H2-Ab1", "Cd74", "Clec10a", "Xcr1"),
  T_cell = c("Cd3d", "Cd3e", "Trac", "Cd247", "Lck"),
  NK_cell = c("Nkg7", "Klrd1", "Klrk1", "Prf1", "Gzmb"),
  B_cell = c("Cd79a", "Cd79b", "Ms4a1", "Cd74", "H2-Ab1"),
  Plasma_cell = c("Jchain", "Mzb1", "Sdc1", "Xbp1", "Igha"),
  HSC_Mesenchymal = c("Dcn", "Col1a1", "Col1a2", "Col3a1", "Lrat", "Rgs5", "Acta2"),
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
  Hepatocyte_common = c("Hnf4a", "Alb", "Fgb", "Ttr", "Ahsg", "Atf5", "Fabp1"),
  Hepatoblast_like = c("Prox1", "Cdh1", "Afp", "Dlk1"),
  Hepatocyte_zone1 = c("Gls2", "Pck1", "Cps1", "Hal", "Ass1"),
  Hepatocyte_zone3 = c("Glul", "Cyp1a2", "Cyp2e1", "Lect2"),
  Cholangiocyte = c("Krt7", "Cldn3", "Cldn4", "Spp1", "Cxcl6", "Krt19"),
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

# Human-only / unavailable / ambiguous entries are retained in audit notes
ueno_notes <- tibble(
  original_gene = c("GNLY", "LYZ", "CD209", "ILR2"),
  mouse_handling = c(
    "No direct mouse ortholog routinely used; excluded",
    "Mapped to Lyz2",
    "Mapped provisionally to Cd209a",
    "Ambiguous symbol; excluded from automatic scoring"
  )
)

# ============================================================
# Availability audit
# ============================================================

audit_marker_sets <- function(marker_sets, source_name) {
  bind_rows(
    lapply(names(marker_sets), function(set_name) {
      tibble(
        source = source_name,
        marker_set = set_name,
        gene = marker_sets[[set_name]],
        available_in_RDS = marker_sets[[set_name]] %in% features_present
      )
    })
  )
}

availability_general <- audit_marker_sets(general_markers, "General")
availability_ueno <- audit_marker_sets(ueno_markers, "Ueno_Custom")

availability_all <- bind_rows(
  availability_general,
  availability_ueno
)

availability_summary <- availability_all %>%
  group_by(source, marker_set) %>%
  summarise(
    requested_genes = n(),
    available_genes = sum(available_in_RDS),
    availability_fraction = mean(available_in_RDS),
    available_gene_list = paste(gene[available_in_RDS], collapse = "; "),
    missing_gene_list = paste(gene[!available_in_RDS], collapse = "; "),
    .groups = "drop"
  )

write_csv(
  availability_all,
  file.path(output_dir, "MarkerAvailability_All.csv")
)

write_csv(
  availability_summary,
  file.path(output_dir, "MarkerAvailability_Summary.csv")
)

write_csv(
  ueno_notes,
  file.path(ueno_dir, "Ueno_MarkerNotes.csv")
)

write_excel(
  file.path(output_dir, "MarkerAvailability.xlsx"),
  list(
    Availability_All = availability_all,
    Availability_Summary = availability_summary,
    Ueno_Notes = ueno_notes
  )
)

# Keep only available genes
filter_available_sets <- function(marker_sets) {
  out <- lapply(marker_sets, function(x) {
    intersect(x, features_present)
  })
  out[lengths(out) >= 2]
}

general_sets_available <- filter_available_sets(general_markers)
ueno_sets_available <- filter_available_sets(ueno_markers)

# ============================================================
# DotPlots
# ============================================================

make_dotplot <- function(obj, marker_sets, title_text) {
  genes_ordered <- unique(unlist(marker_sets, use.names = FALSE))
  genes_ordered <- genes_ordered[genes_ordered %in% features_present]

  DotPlot(
    object = obj,
    features = genes_ordered,
    assay = "RNA",
    group.by = cluster_col,
    scale = TRUE,
    dot.scale = 6
  ) +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title = title_text,
      x = "Marker gene",
      y = "Cluster"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_text(face = "bold")
    )
}

p_general <- make_dotplot(
  obj,
  general_sets_available,
  "RDS3 general marker validation"
)

save_pdf(
  p_general,
  file.path(general_dir, "General_Marker_DotPlot.pdf"),
  width = 22,
  height = 12
)

p_ueno <- make_dotplot(
  obj,
  ueno_sets_available,
  "RDS3 Ueno custom marker validation"
)

save_pdf(
  p_ueno,
  file.path(ueno_dir, "Ueno_Custom_Marker_DotPlot.pdf"),
  width = 24,
  height = 13
)

# ============================================================
# Module scores
# ============================================================

add_scores <- function(obj, marker_sets, prefix) {
  score_map <- tibble(
    marker_set = names(marker_sets),
    score_column = character(length(marker_sets))
  )

  for (i in seq_along(marker_sets)) {
    set_name <- names(marker_sets)[i]
    prefix_i <- paste0(prefix, "_", sprintf("%02d", i), "_")

    obj <- AddModuleScore(
      object = obj,
      features = list(marker_sets[[i]]),
      assay = "RNA",
      name = prefix_i,
      search = FALSE,
      seed = 260730 + i
    )

    generated_col <- paste0(prefix_i, "1")
    score_map$score_column[i] <- generated_col
  }

  list(
    object = obj,
    score_map = score_map
  )
}

cat("Calculating general module scores...\n")
general_scored <- add_scores(
  obj,
  general_sets_available,
  "GeneralScore"
)

obj <- general_scored$object
general_score_map <- general_scored$score_map

cat("Calculating Ueno custom module scores...\n")
ueno_scored <- add_scores(
  obj,
  ueno_sets_available,
  "UenoScore"
)

obj <- ueno_scored$object
ueno_score_map <- ueno_scored$score_map

meta <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    cluster_validation = as.character(.data[[cluster_col]]),
    existing_celltype = as.character(.data[[existing_celltype_col]])
  )

summarize_scores <- function(meta, score_map, source_name) {
  score_columns <- score_map$score_column

  long <- meta %>%
    select(
      cell_barcode,
      cluster_validation,
      existing_celltype,
      all_of(score_columns)
    ) %>%
    pivot_longer(
      cols = all_of(score_columns),
      names_to = "score_column",
      values_to = "module_score"
    ) %>%
    left_join(score_map, by = "score_column") %>%
    mutate(source = source_name)

  cluster_summary <- long %>%
    group_by(
      source,
      cluster_validation,
      existing_celltype,
      marker_set
    ) %>%
    summarise(
      cells = n(),
      mean_score = mean(module_score, na.rm = TRUE),
      median_score = median(module_score, na.rm = TRUE),
      q25 = quantile(module_score, 0.25, na.rm = TRUE),
      q75 = quantile(module_score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    cell_long = long,
    cluster_summary = cluster_summary
  )
}

general_summary <- summarize_scores(
  meta,
  general_score_map,
  "General"
)

ueno_summary <- summarize_scores(
  meta,
  ueno_score_map,
  "Ueno_Custom"
)

write_csv(
  general_summary$cluster_summary,
  file.path(general_dir, "General_ModuleScores_ByCluster.csv")
)

write_csv(
  ueno_summary$cluster_summary,
  file.path(ueno_dir, "Ueno_ModuleScores_ByCluster.csv")
)

# ============================================================
# Cluster-level predictions
# ============================================================

predict_from_scores <- function(cluster_summary, source_name) {
  cluster_summary %>%
    group_by(cluster_validation) %>%
    arrange(desc(mean_score), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    summarise(
      source = source_name,
      existing_celltype = first(existing_celltype),
      predicted_label = marker_set[1],
      top_score = mean_score[1],
      second_label = ifelse(n() >= 2, marker_set[2], NA_character_),
      second_score = ifelse(n() >= 2, mean_score[2], NA_real_),
      score_margin = ifelse(n() >= 2, mean_score[1] - mean_score[2], NA_real_),
      .groups = "drop"
    )
}

general_prediction <- predict_from_scores(
  general_summary$cluster_summary,
  "General"
)

ueno_prediction <- predict_from_scores(
  ueno_summary$cluster_summary,
  "Ueno_Custom"
)

comparison <- general_prediction %>%
  select(
    cluster_validation,
    existing_celltype,
    general_prediction = predicted_label,
    general_top_score = top_score,
    general_second_label = second_label,
    general_score_margin = score_margin
  ) %>%
  left_join(
    ueno_prediction %>%
      select(
        cluster_validation,
        ueno_prediction = predicted_label,
        ueno_top_score = top_score,
        ueno_second_label = second_label,
        ueno_score_margin = score_margin
      ),
    by = "cluster_validation"
  ) %>%
  mutate(
    general_matches_existing =
      str_detect(
        str_to_lower(existing_celltype),
        fixed(str_to_lower(general_prediction))
      ) |
      str_detect(
        str_to_lower(general_prediction),
        fixed(str_to_lower(existing_celltype))
      ),
    review_status = case_when(
      is.na(general_prediction) ~ "Review",
      general_score_margin < 0.05 ~ "Low_margin",
      general_matches_existing ~ "Concordant",
      TRUE ~ "Discordant"
    )
  ) %>%
  arrange(
    suppressWarnings(as.numeric(cluster_validation)),
    cluster_validation
  )

write_csv(
  general_prediction,
  file.path(general_dir, "General_ClusterPrediction.csv")
)

write_csv(
  ueno_prediction,
  file.path(ueno_dir, "Ueno_ClusterPrediction.csv")
)

write_csv(
  comparison,
  file.path(comparison_dir, "General_vs_Ueno_Classification.csv")
)

write_excel(
  file.path(general_dir, "General_ClusterPrediction.xlsx"),
  list(
    Prediction = general_prediction,
    ModuleScores = general_summary$cluster_summary
  )
)

write_excel(
  file.path(ueno_dir, "Ueno_ClusterPrediction.xlsx"),
  list(
    Prediction = ueno_prediction,
    ModuleScores = ueno_summary$cluster_summary
  )
)

write_excel(
  file.path(comparison_dir, "General_vs_Ueno_Classification.xlsx"),
  list(
    Comparison = comparison,
    General = general_prediction,
    Ueno_Custom = ueno_prediction,
    Availability = availability_summary
  )
)

# ============================================================
# Score heatmaps
# ============================================================

make_score_heatmap <- function(cluster_summary, title_text) {
  ggplot(
    cluster_summary,
    aes(
      x = marker_set,
      y = factor(
        cluster_validation,
        levels = rev(cluster_levels)
      ),
      fill = mean_score
    )
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title = title_text,
      x = "Marker module",
      y = "Cluster",
      fill = "Mean score"
    ) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_text(face = "bold")
    )
}

p_general_heatmap <- make_score_heatmap(
  general_summary$cluster_summary,
  "General marker module scores by cluster"
)

save_pdf(
  p_general_heatmap,
  file.path(general_dir, "General_ModuleScore_Heatmap.pdf"),
  width = 16,
  height = 12
)

p_ueno_heatmap <- make_score_heatmap(
  ueno_summary$cluster_summary,
  "Ueno custom marker module scores by cluster"
)

save_pdf(
  p_ueno_heatmap,
  file.path(ueno_dir, "Ueno_ModuleScore_Heatmap.pdf"),
  width = 18,
  height = 12
)

# ============================================================
# Save RDS with added score metadata
# ============================================================

score_rds_file <- file.path(
  output_dir,
  "Mouse_RDS3_with_General_and_Ueno_ModuleScores.rds"
)

saveRDS(
  obj,
  score_rds_file,
  compress = FALSE
)

# ============================================================
# Logs
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

capture.output(
  warnings(),
  file = file.path(log_dir, "warnings.txt")
)

cat("\n============================================================\n")
cat("RDS3 Phase 2-2 completed\n")
cat("============================================================\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main outputs:\n")
cat("  General/General_Marker_DotPlot.pdf\n")
cat("  General/General_ModuleScore_Heatmap.pdf\n")
cat("  General/General_ClusterPrediction.xlsx\n")
cat("  Ueno_Custom/Ueno_Custom_Marker_DotPlot.pdf\n")
cat("  Ueno_Custom/Ueno_ModuleScore_Heatmap.pdf\n")
cat("  Ueno_Custom/Ueno_ClusterPrediction.xlsx\n")
cat("  Comparison/General_vs_Ueno_Classification.xlsx\n")
cat("  MarkerAvailability.xlsx\n")
cat("  Mouse_RDS3_with_General_and_Ueno_ModuleScores.rds\n")
