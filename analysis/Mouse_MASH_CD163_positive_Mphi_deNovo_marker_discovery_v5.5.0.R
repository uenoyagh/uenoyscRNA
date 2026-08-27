#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# Cd163-positive macrophage de novo marker discovery
#
# Version: v5.5.0
#
# GOAL
#   Discover genes that identify the same STD macrophage population as Cd163
#   without restricting the search to pre-selected candidate genes.
#
# DISCOVERY DESIGN
#   1. Load frozen Clean-B macrophage object (FINAL v4.14.5 / Res2.0).
#   2. Restrict discovery to STD.
#   3. Restrict to Kupffer/Macrophage cells when lineage metadata is available.
#   4. Define Cd163+ by RNA detection (>0 counts if available; otherwise data >0).
#   5. For every detected gene, calculate:
#        - expression prevalence in Cd163+ and Cd163-
#        - delta prevalence
#        - mean expression in Cd163+ and Cd163-
#        - log2 mean-expression enrichment
#        - Spearman correlation with Cd163
#        - binary overlap: sensitivity / PPV / specificity / Jaccard
#        - AUROC for Cd163+ vs Cd163-
#   6. Calculate Res2.0 cluster-level expression-pattern correlation with Cd163.
#   7. Integrate the above into an exploratory de novo score.
#   8. Export broad all-gene ranking and stringent candidate shortlist.
#
# IMPORTANT
#   - This is an exploratory single-cell discovery analysis.
#   - Cells are NOT independent biological replicates.
#   - P-values are intentionally not used for biological inference.
#   - No reintegration, reclustering, or re-UMAP is performed.
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(5500)

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

# ==============================================================================
# 1. Paths
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS_CANDIDATES <- c(
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
  ),
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
  )
)

MPHI_RDS <- MPHI_RDS_CANDIDATES[file.exists(MPHI_RDS_CANDIDATES)][1]

OUTPUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163_Positive_deNovo_Marker_Discovery_v5.5.0"
)

DIR_TABLE <- file.path(OUTPUT_DIR, "Tables")
DIR_PLOT <- file.path(OUTPUT_DIR, "Plots")
DIR_UMAP <- file.path(OUTPUT_DIR, "UMAP")
DIR_LOG <- file.path(OUTPUT_DIR, "Logs")

for (d in c(
  OUTPUT_DIR,
  DIR_TABLE,
  DIR_PLOT,
  DIR_UMAP,
  DIR_LOG
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. Analysis settings
# ==============================================================================

ASSAY_USE <- "RNA"

REFERENCE_GENE <- "Cd163"

DISCOVERY_CONDITION <- "STD"

RESTRICT_TO_KUPFFER_MACROPHAGE <- TRUE

MACROPHAGE_LINEAGE_VALUE <- "Kupffer_Macrophage"

LINEAGE_COLUMN_CANDIDATES <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_v440",
  "layer1_original",
  "celltype_for_R8plot",
  "celltype_auto_annotation"
)

CONDITION_COLUMN_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group",
  "sample_for_annotation",
  "sample"
)

CLUSTER_COLUMN_CANDIDATES <- c(
  "mphi_rpca_res_2.0",
  "mphi_rpca_cluster",
  "mphi_rpca_clusters"
)

UMAP_CANDIDATES <- c(
  "umapRPCA",
  "umap"
)

# Basic expression filter before de novo ranking.
# A gene must be detected in at least this fraction of Cd163+ OR Cd163- cells.
MIN_DETECTION_FRACTION <- 0.01

# Shortlist filters.
MIN_PCT_POS <- 10
MIN_DELTA_PCT <- 5
MIN_AUC <- 0.60
MIN_JACCARD <- 0.05
MIN_SPEARMAN <- 0.10

TOP_N_TABLE <- 100
TOP_N_PLOTS <- 24

# Previously examined markers are NOT used to define discovery.
# They are added only as reference annotations in output tables.
REFERENCE_PANEL <- c(
  "Trem2",
  "Spp1",
  "Cd163",
  "Cd68",
  "Adgre1",
  "Stat1",
  "Mrc1",
  "Mertk",
  "Mfge8",
  "Gas6",
  "Il1rn",
  "Igf1",
  "Hmox1"
)

# ==============================================================================
# 3. Utility functions
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(
  plot,
  filename,
  width,
  height
) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
}

resolve_first <- function(
  available,
  candidates
) {
  hit <- candidates[
    candidates %in% available
  ]

  if (!length(hit)) {
    return(
      NA_character_
    )
  }

  hit[[1]]
}

canonical_condition <- function(x) {

  x <- as.character(x)

  out <- case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )

  factor(
    out,
    levels = c(
      "STD",
      "CDAHFD",
      "Sham",
      "Tx"
    )
  )
}

resolve_condition <- function(obj) {

  meta_cols <- colnames(
    obj@meta.data
  )

  source <- resolve_first(
    meta_cols,
    CONDITION_COLUMN_CANDIDATES
  )

  if (is.na(source)) {
    stop(
      "No usable condition/sample metadata column was found."
    )
  }

  cond <- canonical_condition(
    obj@meta.data[[source]]
  )

  # Fallback if the first available source does not map.
  if (all(is.na(cond))) {

    for (src in CONDITION_COLUMN_CANDIDATES) {

      if (!src %in% meta_cols) {
        next
      }

      temp <- canonical_condition(
        obj@meta.data[[src]]
      )

      if (!all(is.na(temp))) {
        source <- src
        cond <- temp
        break
      }
    }
  }

  list(
    source = source,
    condition = cond
  )
}

safe_rescale01 <- function(x) {

  out <- rep(
    NA_real_,
    length(x)
  )

  ok <- is.finite(x)

  if (!any(ok)) {
    return(out)
  }

  r <- range(
    x[ok],
    na.rm = TRUE
  )

  if (
    !all(is.finite(r)) ||
    diff(r) == 0
  ) {
    out[ok] <- 0.5
    return(out)
  }

  out[ok] <- scales::rescale(
    x[ok],
    to = c(
      0,
      1
    )
  )

  out
}

safe_spearman <- function(
  x,
  y
) {

  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      NA_real_
    )
  }

  suppressWarnings(
    cor(
      x,
      y,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
}

binary_metrics <- function(
  reference_positive,
  candidate_positive
) {

  tp <- sum(
    reference_positive &
      candidate_positive
  )

  fp <- sum(
    !reference_positive &
      candidate_positive
  )

  fn <- sum(
    reference_positive &
      !candidate_positive
  )

  tn <- sum(
    !reference_positive &
      !candidate_positive
  )

  sensitivity <- if (
    tp + fn > 0
  ) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  ppv <- if (
    tp + fp > 0
  ) {
    tp / (tp + fp)
  } else {
    NA_real_
  }

  specificity <- if (
    tn + fp > 0
  ) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  jaccard <- if (
    tp + fp + fn > 0
  ) {
    tp / (tp + fp + fn)
  } else {
    NA_real_
  }

  c(
    sensitivity = sensitivity,
    ppv = ppv,
    specificity = specificity,
    jaccard = jaccard
  )
}

auc_rank <- function(
  score,
  label
) {

  label <- as.logical(label)

  n1 <- sum(label)
  n0 <- sum(!label)

  if (
    n1 == 0 ||
    n0 == 0
  ) {
    return(
      NA_real_
    )
  }

  if (
    length(unique(score)) < 2
  ) {
    return(
      0.5
    )
  }

  ranks <- rank(
    score,
    ties.method = "average"
  )

  u <- sum(
    ranks[label]
  ) - n1 * (n1 + 1) / 2

  as.numeric(
    u / (n1 * n0)
  )
}

matrix_row_mean <- function(
  mat,
  cols
) {

  if (!length(cols)) {
    return(
      rep(
        NA_real_,
        nrow(mat)
      )
    )
  }

  Matrix::rowMeans(
    mat[
      ,
      cols,
      drop = FALSE
    ]
  )
}

matrix_row_pct_positive <- function(
  mat,
  cols
) {

  if (!length(cols)) {
    return(
      rep(
        NA_real_,
        nrow(mat)
      )
    )
  }

  Matrix::rowMeans(
    mat[
      ,
      cols,
      drop = FALSE
    ] > 0
  )
}

# ==============================================================================
# 4. Input validation and load
# ==============================================================================

if (
  length(MPHI_RDS) == 0L ||
  is.na(MPHI_RDS) ||
  !file.exists(MPHI_RDS)
) {
  stop(
    paste0(
      "Clean-B macrophage RDS not found.\n",
      "Checked:\n",
      paste(
        MPHI_RDS_CANDIDATES,
        collapse = "\n"
      )
    )
  )
}

msg(
  "Loading Clean-B macrophage RDS..."
)

mphi <- readRDS(
  MPHI_RDS
)

if (
  !ASSAY_USE %in%
    Assays(mphi)
) {
  stop(
    "RNA assay is missing from macrophage object."
  )
}

DefaultAssay(
  mphi
) <- ASSAY_USE

if (
  !REFERENCE_GENE %in%
    rownames(mphi)
) {
  stop(
    "Cd163 is absent from the macrophage object."
  )
}

# ==============================================================================
# 5. Resolve metadata
# ==============================================================================

condition_info <- resolve_condition(
  mphi
)

mphi$condition_v550 <- condition_info$condition

lineage_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  LINEAGE_COLUMN_CANDIDATES
)

cluster_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  CLUSTER_COLUMN_CANDIDATES
)

umap_name <- resolve_first(
  Reductions(mphi),
  UMAP_CANDIDATES
)

msg(
  "Condition source: ",
  condition_info$source
)

msg(
  "Lineage column: ",
  ifelse(
    is.na(lineage_col),
    "<not found>",
    lineage_col
  )
)

msg(
  "Res2 cluster column: ",
  ifelse(
    is.na(cluster_col),
    "<not found>",
    cluster_col
  )
)

msg(
  "UMAP reduction: ",
  ifelse(
    is.na(umap_name),
    "<not found>",
    umap_name
  )
)

# ==============================================================================
# 6. Restrict to STD macrophages
# ==============================================================================

keep <- as.character(
  mphi$condition_v550
) == DISCOVERY_CONDITION

if (
  RESTRICT_TO_KUPFFER_MACROPHAGE &&
  !is.na(lineage_col)
) {

  lineage <- as.character(
    mphi@meta.data[[lineage_col]]
  )

  keep <- keep &
    lineage == MACROPHAGE_LINEAGE_VALUE
}

discovery_cells <- colnames(mphi)[
  keep &
    !is.na(keep)
]

if (
  length(discovery_cells) < 50
) {
  stop(
    paste0(
      "Too few discovery cells after filtering: ",
      length(discovery_cells)
    )
  )
}

std <- subset(
  mphi,
  cells = discovery_cells
)

DefaultAssay(
  std
) <- ASSAY_USE

msg(
  "STD discovery cells: ",
  ncol(std)
)

# ==============================================================================
# 7. Extract RNA matrices and define Cd163 positivity
# ==============================================================================

data_mat <- GetAssayData(
  std,
  assay = ASSAY_USE,
  layer = "data"
)

count_mat <- tryCatch(
  GetAssayData(
    std,
    assay = ASSAY_USE,
    layer = "counts"
  ),
  error = function(e) {
    NULL
  }
)

if (
  !is.null(count_mat) &&
  REFERENCE_GENE %in%
    rownames(count_mat)
) {

  positivity_source <- "RNA counts > 0"

  reference_detect <- as.numeric(
    count_mat[
      REFERENCE_GENE,
      ,
      drop = TRUE
    ]
  ) > 0

  binary_mat <- count_mat > 0

} else {

  positivity_source <- "RNA data > 0"

  reference_detect <- as.numeric(
    data_mat[
      REFERENCE_GENE,
      ,
      drop = TRUE
    ]
  ) > 0

  binary_mat <- data_mat > 0
}

reference_expr <- as.numeric(
  data_mat[
    REFERENCE_GENE,
    ,
    drop = TRUE
  ]
)

n_pos <- sum(
  reference_detect
)

n_neg <- sum(
  !reference_detect
)

pct_cd163_pos <- 100 *
  mean(
    reference_detect
  )

msg(
  "Cd163 positivity definition: ",
  positivity_source
)

msg(
  "Cd163+ cells: ",
  n_pos,
  " / ",
  ncol(std),
  " (",
  round(
    pct_cd163_pos,
    2
  ),
  "%)"
)

if (
  n_pos < 20 ||
  n_neg < 20
) {
  stop(
    paste0(
      "Cd163+ or Cd163- group is too small for de novo discovery. ",
      "Cd163+=",
      n_pos,
      ", Cd163-=",
      n_neg
    )
  )
}

# Save Cd163 status in object for plotting.
std$Cd163_status_v550 <- factor(
  ifelse(
    reference_detect,
    "Cd163+",
    "Cd163-"
  ),
  levels = c(
    "Cd163-",
    "Cd163+"
  )
)

# ==============================================================================
# 8. Global vectorized expression summaries
# ==============================================================================

genes <- rownames(
  data_mat
)

pos_idx <- which(
  reference_detect
)

neg_idx <- which(
  !reference_detect
)

mean_pos <- matrix_row_mean(
  data_mat,
  pos_idx
)

mean_neg <- matrix_row_mean(
  data_mat,
  neg_idx
)

pct_pos <- 100 *
  matrix_row_pct_positive(
    binary_mat,
    pos_idx
  )

pct_neg <- 100 *
  matrix_row_pct_positive(
    binary_mat,
    neg_idx
  )

delta_pct <- pct_pos -
  pct_neg

log2_mean_ratio <- log2(
  (
    mean_pos +
      0.01
  ) /
    (
      mean_neg +
        0.01
    )
)

max_detection_fraction <- pmax(
  pct_pos,
  pct_neg
) / 100

eligible <- max_detection_fraction >=
  MIN_DETECTION_FRACTION

genes_eligible <- genes[
  eligible
]

msg(
  "Genes passing detection filter: ",
  length(genes_eligible),
  " / ",
  length(genes)
)

# ==============================================================================
# 9. Per-gene Cd163 similarity metrics
# ==============================================================================

metric_list <- vector(
  "list",
  length(
    genes_eligible
  )
)

names(
  metric_list
) <- genes_eligible

for (
  i in seq_along(
    genes_eligible
  )
) {

  gene <- genes_eligible[[i]]

  expr <- as.numeric(
    data_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  detect <- as.logical(
    binary_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  bm <- binary_metrics(
    reference_detect,
    detect
  )

  metric_list[[i]] <- tibble(

    gene = gene,

    pct_positive_in_Cd163pos =
      pct_pos[
        match(
          gene,
          genes
        )
      ],

    pct_positive_in_Cd163neg =
      pct_neg[
        match(
          gene,
          genes
        )
      ],

    delta_pct_positive =
      delta_pct[
        match(
          gene,
          genes
        )
      ],

    mean_expression_Cd163pos =
      mean_pos[
        match(
          gene,
          genes
        )
      ],

    mean_expression_Cd163neg =
      mean_neg[
        match(
          gene,
          genes
        )
      ],

    log2_mean_expression_ratio =
      log2_mean_ratio[
        match(
          gene,
          genes
        )
      ],

    spearman_with_Cd163 =
      safe_spearman(
        reference_expr,
        expr
      ),

    sensitivity_to_Cd163pos =
      unname(
        bm[
          "sensitivity"
        ]
      ),

    ppv_for_Cd163pos =
      unname(
        bm[
          "ppv"
        ]
      ),

    specificity_for_Cd163neg =
      unname(
        bm[
          "specificity"
        ]
      ),

    jaccard_with_Cd163 =
      unname(
        bm[
          "jaccard"
        ]
      ),

    auc_Cd163pos_vs_neg =
      auc_rank(
        expr,
        reference_detect
      )
  )
}

gene_metrics <- bind_rows(
  metric_list
)

# Remove the reference itself from de novo candidates.
gene_metrics <- gene_metrics %>%
  filter(
    gene != REFERENCE_GENE
  )

# ==============================================================================
# 10. Res2.0 cluster-pattern similarity
# ==============================================================================

if (
  !is.na(cluster_col)
) {

  cluster_id <- as.character(
    std@meta.data[[cluster_col]]
  )

  cluster_levels <- sort(
    unique(cluster_id)
  )

  cluster_summary_list <- list()

  for (
    cl in cluster_levels
  ) {

    idx <- which(
      cluster_id == cl
    )

    if (!length(idx)) {
      next
    }

    cluster_summary_list[[cl]] <- tibble(
      cluster = cl,
      gene = genes_eligible,
      mean_expression =
        as.numeric(
          Matrix::rowMeans(
            data_mat[
              genes_eligible,
              idx,
              drop = FALSE
            ]
          )
        ),
      pct_positive =
        100 *
        as.numeric(
          Matrix::rowMeans(
            binary_mat[
              genes_eligible,
              idx,
              drop = FALSE
            ]
          )
        )
    )
  }

  cluster_expression_long <- bind_rows(
    cluster_summary_list
  )

  write.csv(
    cluster_expression_long,
    file.path(
      DIR_TABLE,
      "03_gene_expression_by_Res2_cluster_STD_v5.5.0.csv"
    ),
    row.names = FALSE
  )

  cd163_cluster <- cluster_expression_long %>%
    filter(
      gene == REFERENCE_GENE
    ) %>%
    select(
      cluster,
      cd163_cluster_mean =
        mean_expression,
      cd163_cluster_pct =
        pct_positive
    )

  candidate_cluster <- cluster_expression_long %>%
    filter(
      gene != REFERENCE_GENE
    ) %>%
    left_join(
      cd163_cluster,
      by = "cluster"
    )

  cluster_pattern <- candidate_cluster %>%
    group_by(
      gene
    ) %>%
    summarise(
      cluster_mean_spearman =
        safe_spearman(
          cd163_cluster_mean,
          mean_expression
        ),
      cluster_pct_spearman =
        safe_spearman(
          cd163_cluster_pct,
          pct_positive
        ),
      .groups = "drop"
    )

} else {

  cluster_expression_long <- NULL

  cluster_pattern <- tibble(
    gene = gene_metrics$gene,
    cluster_mean_spearman =
      NA_real_,
    cluster_pct_spearman =
      NA_real_
  )
}

write.csv(
  cluster_pattern,
  file.path(
    DIR_TABLE,
    "04_Cd163_Res2_cluster_pattern_similarity_v5.5.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 11. Integrated de novo ranking
# ==============================================================================

ranking <- gene_metrics %>%
  left_join(
    cluster_pattern,
    by = "gene"
  ) %>%
  mutate(

    in_previous_panel =
      gene %in%
      REFERENCE_PANEL,

    score_delta_pct =
      safe_rescale01(
        pmax(
          delta_pct_positive,
          0
        )
      ),

    score_auc =
      safe_rescale01(
        pmax(
          auc_Cd163pos_vs_neg -
            0.5,
          0
        )
      ),

    score_spearman =
      safe_rescale01(
        pmax(
          spearman_with_Cd163,
          0
        )
      ),

    score_jaccard =
      safe_rescale01(
        jaccard_with_Cd163
      ),

    score_ppv =
      safe_rescale01(
        ppv_for_Cd163pos
      ),

    score_sensitivity =
      safe_rescale01(
        sensitivity_to_Cd163pos
      ),

    score_cluster_pattern =
      safe_rescale01(
        rowMeans(
          cbind(
            pmax(
              cluster_mean_spearman,
              0
            ),
            pmax(
              cluster_pct_spearman,
              0
            )
          ),
          na.rm = TRUE
        )
      ),

    de_novo_score =
      0.20 *
        score_delta_pct +
      0.20 *
        score_auc +
      0.15 *
        score_spearman +
      0.15 *
        score_jaccard +
      0.10 *
        score_ppv +
      0.10 *
        score_sensitivity +
      0.10 *
        score_cluster_pattern
  ) %>%
  arrange(
    desc(
      de_novo_score
    )
  ) %>%
  mutate(
    de_novo_rank =
      row_number()
  )

write.csv(
  ranking,
  file.path(
    DIR_TABLE,
    "05_ALL_GENES_Cd163_deNovo_ranking_v5.5.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 12. Stringent shortlist
# ==============================================================================

shortlist <- ranking %>%
  filter(
    pct_positive_in_Cd163pos >=
      MIN_PCT_POS,
    delta_pct_positive >=
      MIN_DELTA_PCT,
    auc_Cd163pos_vs_neg >=
      MIN_AUC,
    jaccard_with_Cd163 >=
      MIN_JACCARD,
    spearman_with_Cd163 >=
      MIN_SPEARMAN
  ) %>%
  arrange(
    desc(
      de_novo_score
    )
  ) %>%
  mutate(
    shortlist_rank =
      row_number()
  )

write.csv(
  shortlist,
  file.path(
    DIR_TABLE,
    "06_STRINGENT_Cd163_deNovo_shortlist_v5.5.0.csv"
  ),
  row.names = FALSE
)

top100 <- ranking %>%
  slice_head(
    n = TOP_N_TABLE
  )

write.csv(
  top100,
  file.path(
    DIR_TABLE,
    "07_TOP100_Cd163_deNovo_candidates_v5.5.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 13. Previous-panel position in de novo search
# ==============================================================================

previous_panel_result <- ranking %>%
  filter(
    in_previous_panel
  ) %>%
  arrange(
    de_novo_rank
  )

write.csv(
  previous_panel_result,
  file.path(
    DIR_TABLE,
    "08_previous_panel_positions_in_deNovo_ranking_v5.5.0.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 14. Cd163+ / Cd163- overview plots
# ==============================================================================

if (
  !is.na(umap_name)
) {

  p_status <- DimPlot(
    std,
    reduction = umap_name,
    group.by = "Cd163_status_v550",
    pt.size = 0.45,
    raster = FALSE
  ) +
    ggtitle(
      paste0(
        "STD macrophages: Cd163 detection (",
        positivity_source,
        ")"
      )
    ) +
    theme_classic(
      base_size = 11
    )

  save_pdf(
    p_status,
    file.path(
      DIR_UMAP,
      "01_STD_Cd163_positive_negative_UMAP_v5.5.0.pdf"
    ),
    7,
    6
  )

  p_cd163 <- FeaturePlot(
    std,
    features = REFERENCE_GENE,
    reduction = umap_name,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    raster = FALSE,
    pt.size = 0.45
  ) +
    scale_colour_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    ggtitle(
      "STD macrophages: Cd163"
    ) +
    theme_classic(
      base_size = 11
    )

  save_pdf(
    p_cd163,
    file.path(
      DIR_UMAP,
      "02_STD_Cd163_FeaturePlot_v5.5.0.pdf"
    ),
    7,
    6
  )
}

# ==============================================================================
# 15. Top candidate FeaturePlots
# ==============================================================================

top_plot_genes <- ranking$gene[
  seq_len(
    min(
      TOP_N_PLOTS,
      nrow(ranking)
    )
  )
]

if (
  !is.na(umap_name) &&
  length(top_plot_genes)
) {

  plot_features <- c(
    REFERENCE_GENE,
    top_plot_genes
  )

  p_top <- FeaturePlot(
    std,
    features = plot_features,
    reduction = umap_name,
    ncol = 5,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    raster = FALSE,
    pt.size = 0.35
  ) &
    scale_colour_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) &
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_top,
    file.path(
      DIR_UMAP,
      "03_TOP24_deNovo_candidates_FeaturePlot_v5.5.0.pdf"
    ),
    16,
    14
  )
}

# ==============================================================================
# 16. Ranking plots
# ==============================================================================

top_rank_plot <- ranking %>%
  slice_head(
    n = 40
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(gene)
    )
  )

p_rank <- ggplot(
  top_rank_plot,
  aes(
    x = de_novo_score,
    y = gene
  )
) +
  geom_col() +
  labs(
    title = "Top 40 de novo genes matching Cd163+ STD macrophages",
    x = "Integrated de novo score",
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_rank,
  file.path(
    DIR_PLOT,
    "01_TOP40_deNovo_ranking_v5.5.0.pdf"
  ),
  8,
  10
)

# ------------------------------------------------------------------------------
# Metric heatmap for top 40
# ------------------------------------------------------------------------------

heat_df <- ranking %>%
  slice_head(
    n = 40
  ) %>%
  select(
    gene,
    score_delta_pct,
    score_auc,
    score_spearman,
    score_jaccard,
    score_ppv,
    score_sensitivity,
    score_cluster_pattern,
    de_novo_score
  ) %>%
  pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "score"
  )

heat_df$gene <- factor(
  heat_df$gene,
  levels = rev(
    unique(
      ranking$gene[
        seq_len(
          min(
            40,
            nrow(ranking)
          )
        )
      ]
    )
  )
)

metric_order <- c(
  "score_delta_pct",
  "score_auc",
  "score_spearman",
  "score_jaccard",
  "score_ppv",
  "score_sensitivity",
  "score_cluster_pattern",
  "de_novo_score"
)

heat_df$metric <- factor(
  heat_df$metric,
  levels = metric_order
)

p_heat <- ggplot(
  heat_df,
  aes(
    x = metric,
    y = gene,
    fill = score
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      1
    )
  ) +
  labs(
    title = "Top 40 Cd163 de novo candidates",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_heat,
  file.path(
    DIR_PLOT,
    "02_TOP40_deNovo_metric_heatmap_v5.5.0.pdf"
  ),
  10,
  11
)

# ==============================================================================
# 17. Cd163+ vs Cd163- violin for top candidates
# ==============================================================================

top_violin_genes <- ranking$gene[
  seq_len(
    min(
      12,
      nrow(ranking)
    )
  )
]

if (
  length(top_violin_genes)
) {

  p_vln <- VlnPlot(
    std,
    features = top_violin_genes,
    group.by = "Cd163_status_v550",
    assay = ASSAY_USE,
    layer = "data",
    pt.size = 0,
    ncol = 4
  ) &
    theme_classic(
      base_size = 9
    )

  save_pdf(
    p_vln,
    file.path(
      DIR_PLOT,
      "03_TOP12_Cd163pos_vs_neg_violin_v5.5.0.pdf"
    ),
    14,
    10
  )
}

# ==============================================================================
# 18. Res2 cluster Cd163 summary
# ==============================================================================

if (
  !is.na(cluster_col)
) {

  cluster_cd163_summary <- tibble(
    cluster =
      as.character(
        std@meta.data[[cluster_col]]
      ),
    Cd163_detect =
      reference_detect,
    Cd163_expression =
      reference_expr
  ) %>%
    group_by(
      cluster
    ) %>%
    summarise(
      n_cells = n(),
      n_Cd163_positive =
        sum(
          Cd163_detect
        ),
      pct_Cd163_positive =
        100 *
        mean(
          Cd163_detect
        ),
      mean_Cd163_expression =
        mean(
          Cd163_expression
        ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(
        pct_Cd163_positive
      )
    )

  write.csv(
    cluster_cd163_summary,
    file.path(
      DIR_TABLE,
      "09_Cd163_expression_by_Res2_cluster_STD_v5.5.0.csv"
    ),
    row.names = FALSE
  )

  cluster_cd163_summary$cluster <- factor(
    cluster_cd163_summary$cluster,
    levels = rev(
      cluster_cd163_summary$cluster
    )
  )

  p_cluster <- ggplot(
    cluster_cd163_summary,
    aes(
      x = pct_Cd163_positive,
      y = cluster
    )
  ) +
    geom_col() +
    labs(
      title = "Cd163-positive fraction by Res2.0 cluster in STD macrophages",
      x = "% Cd163-positive",
      y = "Res2.0 cluster"
    ) +
    theme_classic(
      base_size = 10
    )

  save_pdf(
    p_cluster,
    file.path(
      DIR_PLOT,
      "04_Cd163_positive_fraction_by_Res2_cluster_v5.5.0.pdf"
    ),
    8,
    7
  )
}

# ==============================================================================
# 19. Analysis metadata and index
# ==============================================================================

analysis_metadata <- tibble(
  parameter = c(
    "script_version",
    "input_RDS",
    "assay",
    "condition_source",
    "discovery_condition",
    "restrict_to_Kupffer_Macrophage",
    "lineage_column",
    "lineage_value",
    "cluster_column",
    "UMAP",
    "Cd163_positivity_definition",
    "n_discovery_cells",
    "n_Cd163_positive",
    "n_Cd163_negative",
    "pct_Cd163_positive",
    "min_detection_fraction",
    "min_shortlist_pct_Cd163pos",
    "min_shortlist_delta_pct",
    "min_shortlist_AUC",
    "min_shortlist_Jaccard",
    "min_shortlist_Spearman"
  ),
  value = c(
    "v5.5.0",
    MPHI_RDS,
    ASSAY_USE,
    condition_info$source,
    DISCOVERY_CONDITION,
    as.character(
      RESTRICT_TO_KUPFFER_MACROPHAGE
    ),
    ifelse(
      is.na(lineage_col),
      "NA",
      lineage_col
    ),
    MACROPHAGE_LINEAGE_VALUE,
    ifelse(
      is.na(cluster_col),
      "NA",
      cluster_col
    ),
    ifelse(
      is.na(umap_name),
      "NA",
      umap_name
    ),
    positivity_source,
    ncol(std),
    n_pos,
    n_neg,
    round(
      pct_cd163_pos,
      4
    ),
    MIN_DETECTION_FRACTION,
    MIN_PCT_POS,
    MIN_DELTA_PCT,
    MIN_AUC,
    MIN_JACCARD,
    MIN_SPEARMAN
  )
)

write.csv(
  analysis_metadata,
  file.path(
    DIR_LOG,
    "analysis_metadata_v5.5.0.csv"
  ),
  row.names = FALSE
)

index <- tibble(
  item = c(
    "All-gene de novo ranking",
    "Stringent shortlist",
    "Top 100 candidates",
    "Previous-panel positions",
    "Cluster-pattern similarity",
    "Cd163 by Res2 cluster",
    "Cd163 +/- UMAP",
    "Top candidate FeaturePlots",
    "Top 40 ranking plot",
    "Top 40 metric heatmap",
    "Top 12 Cd163+ vs Cd163- violins"
  ),
  path = c(
    file.path(
      DIR_TABLE,
      "05_ALL_GENES_Cd163_deNovo_ranking_v5.5.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "06_STRINGENT_Cd163_deNovo_shortlist_v5.5.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "07_TOP100_Cd163_deNovo_candidates_v5.5.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "08_previous_panel_positions_in_deNovo_ranking_v5.5.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "04_Cd163_Res2_cluster_pattern_similarity_v5.5.0.csv"
    ),
    file.path(
      DIR_TABLE,
      "09_Cd163_expression_by_Res2_cluster_STD_v5.5.0.csv"
    ),
    file.path(
      DIR_UMAP,
      "01_STD_Cd163_positive_negative_UMAP_v5.5.0.pdf"
    ),
    file.path(
      DIR_UMAP,
      "03_TOP24_deNovo_candidates_FeaturePlot_v5.5.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "01_TOP40_deNovo_ranking_v5.5.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "02_TOP40_deNovo_metric_heatmap_v5.5.0.pdf"
    ),
    file.path(
      DIR_PLOT,
      "03_TOP12_Cd163pos_vs_neg_violin_v5.5.0.pdf"
    )
  )
)

write.csv(
  index,
  file.path(
    OUTPUT_DIR,
    "Cd163_deNovo_MARKER_INDEX_v5.5.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_LOG,
    "sessionInfo_v5.5.0.txt"
  )
)

# ==============================================================================
# 20. Final report
# ==============================================================================

msg(
  "DONE."
)

msg(
  "Output: ",
  OUTPUT_DIR
)

cat(
  "\n============================================================\n"
)

cat(
  "Cd163-positive macrophage de novo discovery v5.5.0\n"
)

cat(
  "============================================================\n"
)

cat(
  "Discovery cells: ",
  ncol(std),
  "\n",
  sep = ""
)

cat(
  "Cd163+ cells: ",
  n_pos,
  " (",
  round(
    pct_cd163_pos,
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "Cd163 positivity: ",
  positivity_source,
  "\n\n",
  sep = ""
)

cat(
  "TOP 30 DE NOVO GENES\n"
)

print(
  ranking %>%
    select(
      de_novo_rank,
      gene,
      de_novo_score,
      pct_positive_in_Cd163pos,
      pct_positive_in_Cd163neg,
      delta_pct_positive,
      auc_Cd163pos_vs_neg,
      spearman_with_Cd163,
      jaccard_with_Cd163,
      ppv_for_Cd163pos,
      sensitivity_to_Cd163pos,
      cluster_mean_spearman,
      cluster_pct_spearman,
      in_previous_panel
    ) %>%
    slice_head(
      n = 30
    )
)

cat(
  "\nSTRINGENT SHORTLIST: ",
  nrow(shortlist),
  " genes\n",
  sep = ""
)

if (
  nrow(shortlist) > 0
) {

  print(
    shortlist %>%
      select(
        shortlist_rank,
        gene,
        de_novo_score,
        pct_positive_in_Cd163pos,
        pct_positive_in_Cd163neg,
        delta_pct_positive,
        auc_Cd163pos_vs_neg,
        spearman_with_Cd163,
        jaccard_with_Cd163,
        ppv_for_Cd163pos,
        sensitivity_to_Cd163pos,
        in_previous_panel
      ) %>%
      slice_head(
        n = 30
      )
  )
}

cat(
  "\n============================================================\n"
)
