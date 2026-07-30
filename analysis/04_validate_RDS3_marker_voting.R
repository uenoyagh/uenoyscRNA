#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(openxlsx)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# RDS3 Phase 2-3
# Marker-voting annotation using:
#   A. General marker sets
#   B. Ueno custom marker sets
#
# Input:
#   Phase2_Markers/CSV/02_AllClusterMarkers.csv
# or
#   Phase2_Markers/CSV/04_Top30Markers_PerCluster.csv
#
# Output:
#   cluster-level marker evidence, voting scores,
#   general/Ueno predictions, concordance tables and heatmaps.
# ============================================================

base_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS3_validation"
)

marker_dir <- file.path(base_dir, "Phase2_Markers", "CSV")

all_markers_file <- file.path(
  marker_dir,
  "02_AllClusterMarkers.csv"
)

top30_file <- file.path(
  marker_dir,
  "04_Top30Markers_PerCluster.csv"
)

phase1_cluster_celltype_file <- file.path(
  base_dir,
  "CSV",
  "03_ClusterCelltypeComposition.csv"
)

output_dir <- file.path(
  base_dir,
  "Phase2_MarkerVoting"
)

general_dir <- file.path(output_dir, "General")
ueno_dir <- file.path(output_dir, "Ueno_Custom")
comparison_dir <- file.path(output_dir, "Comparison")
log_dir <- file.path(output_dir, "Logs")

for (x in c(output_dir, general_dir, ueno_dir, comparison_dir, log_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# Analysis settings
# ============================================================

# Evidence thresholds
padj_threshold <- 0.05
log2fc_threshold <- 0.25
pct1_threshold <- 0.10
detection_difference_threshold <- 0.05

# Weights
weight_log2fc <- 1.0
weight_detection_difference <- 1.0
weight_significance <- 0.25
weight_marker_specificity <- 1.0

# Confidence thresholds based on score margin
high_margin <- 2.0
moderate_margin <- 0.75

# ============================================================
# Marker definitions
# ============================================================

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

# ============================================================
# Helper functions
# ============================================================

sanitize_sheet <- function(x) {
  x <- gsub("[\\\\/:*?\\[\\]]", "_", x)
  substr(x, 1, 31)
}

write_excel <- function(path, sheets) {
  wb <- createWorkbook()

  for (nm in names(sheets)) {
    sheet_name <- sanitize_sheet(nm)
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

natural_order <- function(x) {
  z <- unique(as.character(x))
  zn <- suppressWarnings(as.numeric(z))

  if (all(!is.na(zn))) {
    as.character(sort(unique(zn)))
  } else {
    sort(z)
  }
}

detect_fc_column <- function(df) {
  candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) == 0) {
    stop(
      paste0(
        "No fold-change column found. Columns: ",
        paste(colnames(df), collapse = ", ")
      )
    )
  }

  hit[[1]]
}

marker_sets_to_table <- function(marker_sets, source_name) {
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

# ============================================================
# Load marker results
# ============================================================

if (file.exists(all_markers_file)) {
  marker_file_used <- all_markers_file
} else if (file.exists(top30_file)) {
  marker_file_used <- top30_file
} else {
  stop(
    paste0(
      "Marker input not found.\nExpected either:\n",
      all_markers_file, "\nor\n", top30_file,
      "\nRun analysis/02_validate_RDS3_markers.R first."
    )
  )
}

cat("Reading marker results:\n", marker_file_used, "\n\n")
markers <- read_csv(marker_file_used, show_col_types = FALSE)

required_columns <- c("cluster", "gene", "pct.1", "pct.2", "p_val_adj")
missing_columns <- setdiff(required_columns, colnames(markers))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Required columns missing: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

fc_col <- detect_fc_column(markers)

markers <- markers %>%
  mutate(
    cluster = as.character(cluster),
    gene = as.character(gene),
    detection_difference = if (
      "detection_difference" %in% colnames(markers)
    ) {
      detection_difference
    } else {
      pct.1 - pct.2
    }
  )

cluster_levels <- natural_order(markers$cluster)

# ============================================================
# Existing annotations from Phase 1
# ============================================================

existing_annotation <- tibble(
  cluster = cluster_levels,
  existing_celltype = NA_character_
)

if (file.exists(phase1_cluster_celltype_file)) {
  phase1 <- read_csv(
    phase1_cluster_celltype_file,
    show_col_types = FALSE
  )

  cluster_col_phase1 <- c(
    "cluster_validation",
    "cluster"
  )
  cluster_col_phase1 <- cluster_col_phase1[
    cluster_col_phase1 %in% colnames(phase1)
  ][1]

  celltype_col_phase1 <- c(
    "celltype_validation",
    "celltype"
  )
  celltype_col_phase1 <- celltype_col_phase1[
    celltype_col_phase1 %in% colnames(phase1)
  ][1]

  if (!is.na(cluster_col_phase1) && !is.na(celltype_col_phase1)) {
    existing_annotation <- phase1 %>%
      transmute(
        cluster = as.character(.data[[cluster_col_phase1]]),
        existing_celltype = as.character(.data[[celltype_col_phase1]]),
        cells = if ("cells" %in% colnames(phase1)) cells else 1
      ) %>%
      group_by(cluster, existing_celltype) %>%
      summarise(cells = sum(cells), .groups = "drop") %>%
      group_by(cluster) %>%
      slice_max(cells, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(cluster, existing_celltype)
  }
}

# ============================================================
# Marker specificity weights
# ============================================================

general_table <- marker_sets_to_table(
  general_markers,
  "General"
)

ueno_table <- marker_sets_to_table(
  ueno_markers,
  "Ueno_Custom"
)

all_marker_set_table <- bind_rows(
  general_table,
  ueno_table
)

specificity_table <- all_marker_set_table %>%
  group_by(source, gene) %>%
  summarise(
    number_of_labels = n_distinct(label),
    specificity_weight = 1 / number_of_labels,
    .groups = "drop"
  )

# ============================================================
# Marker-voting engine
# ============================================================

run_voting <- function(markers, marker_set_table, source_name) {

  evidence <- markers %>%
    inner_join(
      marker_set_table %>%
        filter(source == source_name),
      by = "gene"
    ) %>%
    left_join(
      specificity_table %>%
        filter(source == source_name),
      by = c("source", "gene")
    ) %>%
    mutate(
      passes_padj = !is.na(p_val_adj) &
        p_val_adj < padj_threshold,
      passes_fc = !is.na(.data[[fc_col]]) &
        .data[[fc_col]] >= log2fc_threshold,
      passes_pct1 = !is.na(pct.1) &
        pct.1 >= pct1_threshold,
      passes_detection_difference =
        !is.na(detection_difference) &
        detection_difference >= detection_difference_threshold,
      evidence_pass =
        passes_padj &
        passes_fc &
        passes_pct1 &
        passes_detection_difference,
      significance_component =
        if_else(
          p_val_adj > 0,
          pmin(-log10(p_val_adj), 50),
          50
        ),
      vote_score = if_else(
        evidence_pass,
        (
          weight_log2fc * pmax(.data[[fc_col]], 0) +
          weight_detection_difference *
            pmax(detection_difference, 0) +
          weight_significance *
            significance_component
        ) *
          (
            weight_marker_specificity *
              specificity_weight
          ),
        0
      )
    )

  all_clusters <- tibble(cluster = cluster_levels)
  all_labels <- marker_set_table %>%
    filter(source == source_name) %>%
    distinct(label)

  complete_grid <- tidyr::crossing(
    all_clusters,
    all_labels
  )

  score_table <- evidence %>%
    group_by(cluster, label) %>%
    summarise(
      requested_markers = n_distinct(gene),
      supporting_markers = n_distinct(gene[evidence_pass]),
      total_vote_score = sum(vote_score, na.rm = TRUE),
      mean_log2fc_support = if_else(
        supporting_markers > 0,
        mean(.data[[fc_col]][evidence_pass], na.rm = TRUE),
        NA_real_
      ),
      mean_detection_difference_support = if_else(
        supporting_markers > 0,
        mean(detection_difference[evidence_pass], na.rm = TRUE),
        NA_real_
      ),
      supporting_gene_list = paste(
        unique(gene[evidence_pass]),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    right_join(
      complete_grid,
      by = c("cluster", "label")
    ) %>%
    mutate(
      requested_markers = replace_na(requested_markers, 0L),
      supporting_markers = replace_na(supporting_markers, 0L),
      total_vote_score = replace_na(total_vote_score, 0),
      supporting_gene_list = replace_na(
        supporting_gene_list,
        ""
      )
    )

  prediction <- score_table %>%
    group_by(cluster) %>%
    arrange(
      desc(total_vote_score),
      desc(supporting_markers),
      desc(mean_log2fc_support),
      .by_group = TRUE
    ) %>%
    summarise(
      predicted_label = label[1],
      top_score = total_vote_score[1],
      top_supporting_markers = supporting_markers[1],
      top_supporting_gene_list = supporting_gene_list[1],
      second_label = ifelse(n() >= 2, label[2], NA_character_),
      second_score = ifelse(n() >= 2, total_vote_score[2], NA_real_),
      second_supporting_markers = ifelse(
        n() >= 2,
        supporting_markers[2],
        NA_integer_
      ),
      score_margin = ifelse(
        n() >= 2,
        total_vote_score[1] - total_vote_score[2],
        NA_real_
      ),
      confidence = case_when(
        top_score <= 0 ~ "No_evidence",
        top_supporting_markers < 2 ~ "Low",
        score_margin >= high_margin ~ "High",
        score_margin >= moderate_margin ~ "Moderate",
        TRUE ~ "Low"
      ),
      .groups = "drop"
    ) %>%
    left_join(
      existing_annotation,
      by = "cluster"
    ) %>%
    arrange(
      match(cluster, cluster_levels)
    )

  list(
    evidence = evidence,
    scores = score_table,
    prediction = prediction
  )
}

general_result <- run_voting(
  markers,
  general_table,
  "General"
)

ueno_result <- run_voting(
  markers,
  ueno_table,
  "Ueno_Custom"
)

# ============================================================
# Concordance
# ============================================================

normalize_label <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]", "")
}

comparison <- general_result$prediction %>%
  select(
    cluster,
    existing_celltype,
    general_prediction = predicted_label,
    general_score = top_score,
    general_supporting_markers = top_supporting_markers,
    general_supporting_genes = top_supporting_gene_list,
    general_second_label = second_label,
    general_score_margin = score_margin,
    general_confidence = confidence
  ) %>%
  left_join(
    ueno_result$prediction %>%
      select(
        cluster,
        ueno_prediction = predicted_label,
        ueno_score = top_score,
        ueno_supporting_markers = top_supporting_markers,
        ueno_supporting_genes = top_supporting_gene_list,
        ueno_second_label = second_label,
        ueno_score_margin = score_margin,
        ueno_confidence = confidence
      ),
    by = "cluster"
  ) %>%
  mutate(
    existing_normalized = normalize_label(existing_celltype),
    general_normalized = normalize_label(general_prediction),
    general_vs_existing = case_when(
      is.na(existing_celltype) ~ "Existing_missing",
      str_detect(existing_normalized, fixed(general_normalized)) |
        str_detect(general_normalized, fixed(existing_normalized)) ~
        "Concordant",
      TRUE ~ "Review"
    ),
    final_review_priority = case_when(
      general_confidence == "No_evidence" ~ "High",
      general_confidence == "Low" ~ "High",
      general_vs_existing == "Review" ~ "High",
      ueno_confidence == "Low" ~ "Moderate",
      TRUE ~ "Routine"
    )
  ) %>%
  arrange(
    match(cluster, cluster_levels)
  )

# ============================================================
# Exports
# ============================================================

write_csv(
  general_result$evidence,
  file.path(general_dir, "General_MarkerEvidence.csv")
)

write_csv(
  general_result$scores,
  file.path(general_dir, "General_VotingScores.csv")
)

write_csv(
  general_result$prediction,
  file.path(general_dir, "General_VotingPrediction.csv")
)

write_csv(
  ueno_result$evidence,
  file.path(ueno_dir, "Ueno_MarkerEvidence.csv")
)

write_csv(
  ueno_result$scores,
  file.path(ueno_dir, "Ueno_VotingScores.csv")
)

write_csv(
  ueno_result$prediction,
  file.path(ueno_dir, "Ueno_VotingPrediction.csv")
)

write_csv(
  comparison,
  file.path(comparison_dir, "General_vs_Ueno_MarkerVoting.csv")
)

write_excel(
  file.path(
    general_dir,
    "General_MarkerVoting.xlsx"
  ),
  list(
    Prediction = general_result$prediction,
    Scores = general_result$scores,
    Evidence = general_result$evidence
  )
)

write_excel(
  file.path(
    ueno_dir,
    "Ueno_MarkerVoting.xlsx"
  ),
  list(
    Prediction = ueno_result$prediction,
    Scores = ueno_result$scores,
    Evidence = ueno_result$evidence
  )
)

write_excel(
  file.path(
    comparison_dir,
    "General_vs_Ueno_MarkerVoting.xlsx"
  ),
  list(
    Comparison = comparison,
    GeneralPrediction = general_result$prediction,
    UenoPrediction = ueno_result$prediction,
    GeneralScores = general_result$scores,
    UenoScores = ueno_result$scores
  )
)

# ============================================================
# Heatmaps
# ============================================================

make_heatmap <- function(score_table, title_text) {
  ggplot(
    score_table,
    aes(
      x = label,
      y = factor(cluster, levels = rev(cluster_levels)),
      fill = total_vote_score
    )
  ) +
    geom_tile() +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#FF1A1A"
    ) +
    labs(
      title = title_text,
      x = "Marker classification",
      y = "Cluster",
      fill = "Vote score"
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

general_heatmap <- make_heatmap(
  general_result$scores,
  "General marker-voting scores by cluster"
)

ggsave(
  file.path(
    general_dir,
    "General_MarkerVoting_Heatmap.pdf"
  ),
  general_heatmap,
  width = 16,
  height = 12,
  device = cairo_pdf,
  limitsize = FALSE
)

ueno_heatmap <- make_heatmap(
  ueno_result$scores,
  "Ueno custom marker-voting scores by cluster"
)

ggsave(
  file.path(
    ueno_dir,
    "Ueno_MarkerVoting_Heatmap.pdf"
  ),
  ueno_heatmap,
  width = 18,
  height = 12,
  device = cairo_pdf,
  limitsize = FALSE
)

# ============================================================
# Review-priority plot
# ============================================================

priority_levels <- c("Routine", "Moderate", "High")

priority_summary <- comparison %>%
  mutate(
    final_review_priority = factor(
      final_review_priority,
      levels = priority_levels
    )
  ) %>%
  count(final_review_priority, name = "clusters")

write_csv(
  priority_summary,
  file.path(comparison_dir, "ReviewPriority_Summary.csv")
)

priority_plot <- ggplot(
  priority_summary,
  aes(
    x = final_review_priority,
    y = clusters
  )
) +
  geom_col() +
  labs(
    title = "Cluster review priority",
    x = "Review priority",
    y = "Number of clusters"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(
    comparison_dir,
    "ReviewPriority_Summary.pdf"
  ),
  priority_plot,
  width = 7,
  height = 5,
  device = cairo_pdf
)

# ============================================================
# Settings and logs
# ============================================================

settings <- tibble(
  item = c(
    "marker_file_used",
    "padj_threshold",
    "log2fc_threshold",
    "pct1_threshold",
    "detection_difference_threshold",
    "weight_log2fc",
    "weight_detection_difference",
    "weight_significance",
    "weight_marker_specificity",
    "high_margin",
    "moderate_margin"
  ),
  value = c(
    marker_file_used,
    padj_threshold,
    log2fc_threshold,
    pct1_threshold,
    detection_difference_threshold,
    weight_log2fc,
    weight_detection_difference,
    weight_significance,
    weight_marker_specificity,
    high_margin,
    moderate_margin
  )
)

write_csv(
  settings,
  file.path(output_dir, "MarkerVoting_Settings.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

capture.output(
  warnings(),
  file = file.path(log_dir, "warnings.txt")
)

cat("\n============================================================\n")
cat("RDS3 Phase 2-3 marker voting completed\n")
cat("============================================================\n")
cat("Marker input:\n", marker_file_used, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Main outputs:\n")
cat("  General/General_MarkerVoting.xlsx\n")
cat("  General/General_MarkerVoting_Heatmap.pdf\n")
cat("  Ueno_Custom/Ueno_MarkerVoting.xlsx\n")
cat("  Ueno_Custom/Ueno_MarkerVoting_Heatmap.pdf\n")
cat("  Comparison/General_vs_Ueno_MarkerVoting.xlsx\n")
cat("  Comparison/ReviewPriority_Summary.pdf\n")
