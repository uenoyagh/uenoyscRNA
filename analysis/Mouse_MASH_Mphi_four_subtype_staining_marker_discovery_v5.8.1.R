#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(5800)

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
# Mouse MASH scRNA-seq
# Four macrophage-subtype de novo staining-marker discovery
# Version: v5.8.1
#
# Targets:
#   Inflammatory-Mphi
#   ECM-associated inflammatory-Mphi
#   Repair/Resolution-Mphi
#   Lipid-associated/TREM2-Mphi
#
# Strategy:
#   1) Search ALL expressed RNA genes; no surface-protein restriction.
#   2) Rank subtype identity independently of treatment.
#   3) Add a secondary score for agreement with the expected Sham->Tx direction.
#   4) Preserve existing Clean-B labels and umapRPCA; no re-clustering/re-UMAP.
#
# Final score = 70% subtype identity + 30% treatment-direction agreement.
# ==============================================================================

msg <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
}

save_pdf <- function(plot_obj, filename, width, height) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
}

canonical_condition <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("STD", "CDAHFD", "Sham", "Tx"))
}

safe_rescale01 <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  rr <- range(x[ok], na.rm = TRUE)
  if (!all(is.finite(rr)) || diff(rr) == 0) {
    out[ok] <- 0.5
    return(out)
  }
  out[ok] <- scales::rescale(x[ok], to = c(0, 1))
  out
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  xx <- x[ok]
  yy <- y[ok]
  if (length(unique(xx)) < 2 || length(unique(yy)) < 2) return(NA_real_)
  suppressWarnings(cor(xx, yy, method = "spearman"))
}

pairwise_direction_fraction <- function(sham_values, tx_values, direction) {
  sham_values <- sham_values[is.finite(sham_values)]
  tx_values <- tx_values[is.finite(tx_values)]
  if (!length(sham_values) || !length(tx_values)) return(NA_real_)
  d <- outer(tx_values, sham_values, FUN = "-")
  if (direction == "increase") return(mean(d > 0))
  if (direction == "decrease") return(mean(d < 0))
  NA_real_
}

make_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# ==============================================================================
# 1. Paths / metadata
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

RDS_CANDIDATES <- c(
  file.path(DATA_ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "RDS", "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"),
  file.path(DATA_ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_Res2_CleanB_FINAL_v4.14.5",
            "RDS", "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS")
)
INPUT_RDS <- RDS_CANDIDATES[file.exists(RDS_CANDIDATES)][1]

OUTPUT_DIR <- file.path(DATA_ROOT, "Mouse_MASH_Mphi_RDS",
                        "Mphi_Subtype_StainingMarker_Discovery_v5.8.1")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
UMAP_DIR <- file.path(OUTPUT_DIR, "UMAP")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")
for (d in c(OUTPUT_DIR, TAB_DIR, FIG_DIR, UMAP_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ASSAY_USE <- "RNA"
UMAP_USE <- "umapRPCA"
CLASS_COLUMN <- "macrophage_class_Res2_FINAL_v4145_char"
CONDITION_COLUMN <- "condition_FIXED2"
SAMPLE_COLUMN <- "sample_for_annotation"
CLUSTER_COLUMN <- "mphi_rpca_res_2.0"

TARGET_CLASSES <- c(
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

EXPECTED_DIRECTION <- c(
  "Inflammatory-Mphi" = "decrease",
  "ECM-associated inflammatory-Mphi" = "decrease",
  "Repair/Resolution-Mphi" = "increase",
  "Lipid-associated/TREM2-Mphi" = "increase"
)

REFERENCE_GENES <- list(
  "Inflammatory-Mphi" = c("Il1b", "Tnf", "Nfkbia", "Ccl3", "Ccl4", "Stat1", "Cxcl2", "Ptgs2"),
  "ECM-associated inflammatory-Mphi" = c("Thbs1", "Fn1", "Spp1", "Tgfb1", "Mmp14", "Mmp12", "Lgals3", "Gpnmb"),
  "Repair/Resolution-Mphi" = c("Mertk", "Mfge8", "Gas6", "Igf1", "Hmox1", "Il10ra", "Tgm2", "Mrc1"),
  "Lipid-associated/TREM2-Mphi" = c("Trem2", "Spp1", "Cd9", "Lpl", "Gpnmb", "Lgals3", "Apoc1", "Fabp5")
)

MIN_DETECTION_FRACTION <- 0.01
GRID_NX <- 35
GRID_NY <- 35
MIN_CELLS_PER_GRID_BIN <- 5
TECHNICAL_GENE_PATTERN <- "^(mt-|Rpl|Rps|Hba-|Hbb-)"
TOP_N_TABLE <- 200
TOP_N_UMAP <- 12
TOP_N_SHAMTX <- 6

# Identity weights
W_ID_SPATIAL <- 0.30
W_ID_DELTA_PCT <- 0.20
W_ID_ENRICHMENT <- 0.20
W_ID_SPECIFICITY <- 0.10
W_ID_PPV <- 0.10
W_ID_JACCARD <- 0.10

# Final integration
W_FINAL_IDENTITY <- 0.70
W_FINAL_TREATMENT <- 0.30

# Treatment direction components
W_TX_EXPR <- 0.35
W_TX_PCT <- 0.30
W_TX_PB_CONS <- 0.175
W_TX_PCT_CONS <- 0.175

# ==============================================================================
# 2. Load / validate
# ==============================================================================

if (length(INPUT_RDS) == 0L || is.na(INPUT_RDS) || !file.exists(INPUT_RDS)) {
  stop("Clean-B macrophage RDS not found.")
}

msg("Loading Clean-B macrophage RDS...")
obj <- readRDS(INPUT_RDS)

if (!ASSAY_USE %in% Assays(obj)) stop("RNA assay missing.")
DefaultAssay(obj) <- ASSAY_USE

required_meta <- c(CLASS_COLUMN, CONDITION_COLUMN, SAMPLE_COLUMN, CLUSTER_COLUMN)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta)) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
if (!UMAP_USE %in% Reductions(obj)) stop("umapRPCA missing.")

obj$class_v580 <- as.character(obj@meta.data[[CLASS_COLUMN]])
obj$condition_v580 <- canonical_condition(obj@meta.data[[CONDITION_COLUMN]])
obj$sample_v580 <- as.character(obj@meta.data[[SAMPLE_COLUMN]])

msg("Cells: ", ncol(obj))
print(table(obj$class_v580, useNA = "ifany"))

missing_classes <- setdiff(TARGET_CLASSES, unique(obj$class_v580))
if (length(missing_classes)) {
  stop("Missing target classes: ", paste(missing_classes, collapse = ", "))
}

counts_mat <- GetAssayData(obj, assay = ASSAY_USE, layer = "counts")
data_mat <- GetAssayData(obj, assay = ASSAY_USE, layer = "data")
all_genes <- rownames(counts_mat)

global_detection <- as.numeric(Matrix::rowMeans(counts_mat > 0))
eligible_genes <- all_genes[global_detection >= MIN_DETECTION_FRACTION]
msg("Genes eligible: ", length(eligible_genes), " / ", length(all_genes))

# ==============================================================================
# 3. UMAP-grid reference
# ==============================================================================

coords <- Embeddings(obj, reduction = UMAP_USE)[colnames(obj), , drop = FALSE]
x <- coords[, 1]
y <- coords[, 2]

x_breaks <- seq(min(x) - 1e-8, max(x) + 1e-8, length.out = GRID_NX + 1)
y_breaks <- seq(min(y) - 1e-8, max(y) + 1e-8, length.out = GRID_NY + 1)
x_bin <- cut(x, x_breaks, include.lowest = TRUE, labels = FALSE)
y_bin <- cut(y, y_breaks, include.lowest = TRUE, labels = FALSE)
grid_id_raw <- (y_bin - 1L) * GRID_NX + x_bin
occupied_grid <- sort(unique(grid_id_raw))
grid_map <- match(grid_id_raw, occupied_grid)

membership <- Matrix::sparseMatrix(
  i = seq_along(grid_map), j = grid_map, x = 1,
  dims = c(length(grid_map), length(occupied_grid))
)
bin_n <- as.numeric(Matrix::colSums(membership))
keep_bins <- which(bin_n >= MIN_CELLS_PER_GRID_BIN)
msg("UMAP bins retained: ", length(keep_bins))

grid_expr_sum <- data_mat %*% membership
grid_expr_mean <- sweep(
  as.matrix(grid_expr_sum[, keep_bins, drop = FALSE]),
  2, bin_n[keep_bins], FUN = "/"
)
rm(grid_expr_sum)
gc(verbose = FALSE)

grid_pos_sum <- (counts_mat > 0) %*% membership
grid_pos_frac <- sweep(
  as.matrix(grid_pos_sum[, keep_bins, drop = FALSE]),
  2, bin_n[keep_bins], FUN = "/"
)
rm(grid_pos_sum)
gc(verbose = FALSE)

# ==============================================================================
# 4. Sample-level Sham/Tx metrics in all Clean-B macrophages
# ==============================================================================

msg("Calculating sample-level Sham/Tx metrics...")

treat_idx <- which(as.character(obj$condition_v580) %in% c("Sham", "Tx"))
treat_counts <- counts_mat[, treat_idx, drop = FALSE]
treat_samples <- obj$sample_v580[treat_idx]
treat_conditions <- as.character(obj$condition_v580[treat_idx])

sample_condition <- tibble(sample = treat_samples, condition = treat_conditions) %>% distinct()
samples_present <- sort(unique(treat_samples))

pb_counts <- matrix(
  0, nrow = length(eligible_genes), ncol = length(samples_present),
  dimnames = list(eligible_genes, samples_present)
)
pct_sample <- matrix(
  NA_real_, nrow = length(eligible_genes), ncol = length(samples_present),
  dimnames = list(eligible_genes, samples_present)
)
sample_n <- setNames(integer(length(samples_present)), samples_present)

for (s in samples_present) {
  ii <- which(treat_samples == s)
  sample_n[[s]] <- length(ii)
  pb_counts[, s] <- as.numeric(Matrix::rowSums(treat_counts[eligible_genes, ii, drop = FALSE]))
  pct_sample[, s] <- 100 * as.numeric(Matrix::rowMeans(treat_counts[eligible_genes, ii, drop = FALSE] > 0))
}

lib_size <- colSums(pb_counts)
pb_cpm <- sweep(pb_counts, 2, lib_size, FUN = "/") * 1e6
pb_log2cpm <- log2(pb_cpm + 1)

sham_samples <- sample_condition %>% filter(condition == "Sham") %>% pull(sample) %>% intersect(colnames(pb_cpm))
tx_samples <- sample_condition %>% filter(condition == "Tx") %>% pull(sample) %>% intersect(colnames(pb_cpm))
if (!length(sham_samples) || !length(tx_samples)) stop("Sham/Tx samples not found.")

Sham_mean_CPM <- rowMeans(pb_cpm[, sham_samples, drop = FALSE], na.rm = TRUE)
Tx_mean_CPM <- rowMeans(pb_cpm[, tx_samples, drop = FALSE], na.rm = TRUE)
log2FC_Tx_vs_Sham <- log2((Tx_mean_CPM + 0.5) / (Sham_mean_CPM + 0.5))
Sham_mean_pct <- rowMeans(pct_sample[, sham_samples, drop = FALSE], na.rm = TRUE)
Tx_mean_pct <- rowMeans(pct_sample[, tx_samples, drop = FALSE], na.rm = TRUE)
delta_pct_Tx_minus_Sham <- Tx_mean_pct - Sham_mean_pct

sample_long <- bind_rows(lapply(samples_present, function(s) {
  tibble(
    gene = eligible_genes,
    sample = s,
    pseudobulk_CPM = pb_cpm[, s],
    pseudobulk_log2CPM = pb_log2cpm[, s],
    pct_positive = pct_sample[, s],
    n_cells = sample_n[[s]]
  )
})) %>% left_join(sample_condition, by = "sample")

write.csv(sample_long,
          file.path(TAB_DIR, "01_ALL_GENES_samplelevel_ShamTx_metrics_v5.8.1.csv"),
          row.names = FALSE)

# Subtype fractions by sample
subtype_fraction <- tibble(
  sample = obj$sample_v580,
  condition = as.character(obj$condition_v580),
  subtype = obj$class_v580
) %>%
  filter(condition %in% c("Sham", "Tx")) %>%
  count(condition, sample, subtype, name = "n_cells") %>%
  group_by(condition, sample) %>%
  mutate(total_Mphi = sum(n_cells), pct_Mphi = 100 * n_cells / total_Mphi) %>%
  ungroup()

write.csv(subtype_fraction,
          file.path(TAB_DIR, "02_subtype_fraction_by_sample_v5.8.1.csv"),
          row.names = FALSE)

# ==============================================================================
# 5. Four-subtype discovery
# ==============================================================================

all_rankings <- list()
all_top200 <- list()
all_reference <- list()

for (subtype in TARGET_CLASSES) {

  msg("Discovering: ", subtype)
  slug <- make_slug(subtype)

  # v5.8.1 fix:
  # Freeze named-list values before entering tibble()/dplyr data masking.
  # Otherwise `subtype` can be interpreted as the newly created tibble column.
  reference_genes_now <- REFERENCE_GENES[[subtype]]
  direction_now <- unname(EXPECTED_DIRECTION[[subtype]])

  target <- obj$class_v580 == subtype
  target_idx <- which(target)
  rest_idx <- which(!target)

  # Identity metrics
  pct_target <- 100 * as.numeric(Matrix::rowMeans(counts_mat[eligible_genes, target_idx, drop = FALSE] > 0))
  pct_rest <- 100 * as.numeric(Matrix::rowMeans(counts_mat[eligible_genes, rest_idx, drop = FALSE] > 0))
  mean_target <- as.numeric(Matrix::rowMeans(data_mat[eligible_genes, target_idx, drop = FALSE]))
  mean_rest <- as.numeric(Matrix::rowMeans(data_mat[eligible_genes, rest_idx, drop = FALSE]))
  delta_target <- pct_target - pct_rest
  expression_ratio <- (mean_target + 1e-4) / (mean_rest + 1e-4)

  positive_all <- counts_mat[eligible_genes, , drop = FALSE] > 0
  positive_target_n <- as.numeric(Matrix::rowSums(positive_all[, target_idx, drop = FALSE]))
  positive_total_n <- as.numeric(Matrix::rowSums(positive_all))

  ppv <- ifelse(positive_total_n > 0, positive_target_n / positive_total_n, NA_real_)
  sensitivity <- pct_target / 100
  specificity <- 1 - pct_rest / 100
  target_size <- length(target_idx)
  union_n <- positive_total_n + target_size - positive_target_n
  jaccard <- ifelse(union_n > 0, positive_target_n / union_n, NA_real_)

  # UMAP spatial similarity: subtype fraction per grid bin vs gene pattern
  subtype_bin_n <- as.numeric(Matrix::colSums(membership[target_idx, , drop = FALSE]))
  subtype_grid_fraction <- subtype_bin_n[keep_bins] / bin_n[keep_bins]

  grid_mean_cor <- rep(NA_real_, length(eligible_genes))
  grid_pct_cor <- rep(NA_real_, length(eligible_genes))
  for (i in seq_along(eligible_genes)) {
    g <- eligible_genes[[i]]
    grid_mean_cor[[i]] <- safe_spearman(
      subtype_grid_fraction,
      as.numeric(grid_expr_mean[g, , drop = TRUE])
    )
    grid_pct_cor[[i]] <- safe_spearman(
      subtype_grid_fraction,
      as.numeric(grid_pos_frac[g, , drop = TRUE])
    )
  }

  spatial_similarity <- rowMeans(
    cbind(pmax(grid_mean_cor, 0), pmax(grid_pct_cor, 0)),
    na.rm = TRUE
  )
  spatial_similarity[!is.finite(spatial_similarity)] <- NA_real_

  # Treatment direction score. Direction is intentionally secondary.
  if (direction_now == "increase") {
    signed_log2fc <- log2FC_Tx_vs_Sham
    signed_delta_tx <- delta_pct_Tx_minus_Sham
  } else {
    signed_log2fc <- -log2FC_Tx_vs_Sham
    signed_delta_tx <- -delta_pct_Tx_minus_Sham
  }

  pb_consistency <- rep(NA_real_, length(eligible_genes))
  pct_consistency <- rep(NA_real_, length(eligible_genes))
  for (i in seq_along(eligible_genes)) {
    g <- eligible_genes[[i]]
    pb_consistency[[i]] <- pairwise_direction_fraction(
      pb_log2cpm[g, sham_samples], pb_log2cpm[g, tx_samples], direction_now
    )
    pct_consistency[[i]] <- pairwise_direction_fraction(
      pct_sample[g, sham_samples], pct_sample[g, tx_samples], direction_now
    )
  }

  # Identity score
  score_spatial <- safe_rescale01(spatial_similarity)
  score_delta <- safe_rescale01(pmax(delta_target, 0))
  score_enrichment <- safe_rescale01(pmax(log2(expression_ratio), 0))
  score_specificity <- safe_rescale01(specificity)
  score_ppv <- safe_rescale01(ppv)
  score_jaccard <- safe_rescale01(jaccard)

  identity_score <-
    W_ID_SPATIAL * score_spatial +
    W_ID_DELTA_PCT * score_delta +
    W_ID_ENRICHMENT * score_enrichment +
    W_ID_SPECIFICITY * score_specificity +
    W_ID_PPV * score_ppv +
    W_ID_JACCARD * score_jaccard

  # Treatment score
  score_tx_expr <- safe_rescale01(pmax(signed_log2fc, 0))
  score_tx_pct <- safe_rescale01(pmax(signed_delta_tx, 0))
  treatment_score <-
    W_TX_EXPR * score_tx_expr +
    W_TX_PCT * score_tx_pct +
    W_TX_PB_CONS * pb_consistency +
    W_TX_PCT_CONS * pct_consistency

  final_score <- W_FINAL_IDENTITY * identity_score + W_FINAL_TREATMENT * treatment_score

  ranking <- tibble(
    subtype = subtype,
    expected_Tx_direction = direction_now,
    gene = eligible_genes,
    pct_positive_target = pct_target,
    pct_positive_other_Mphi = pct_rest,
    delta_pct_target_minus_other = delta_target,
    mean_expression_target = mean_target,
    mean_expression_other_Mphi = mean_rest,
    expression_ratio_target_vs_other = expression_ratio,
    sensitivity = sensitivity,
    specificity = specificity,
    PPV_gene_positive_cells_in_target = ppv,
    jaccard_gene_positive_with_target = jaccard,
    UMAP_grid_mean_spearman = grid_mean_cor,
    UMAP_grid_pct_spearman = grid_pct_cor,
    UMAP_spatial_similarity = spatial_similarity,
    Sham_mean_CPM = unname(Sham_mean_CPM[eligible_genes]),
    Tx_mean_CPM = unname(Tx_mean_CPM[eligible_genes]),
    allMphi_log2FC_Tx_vs_Sham = unname(log2FC_Tx_vs_Sham[eligible_genes]),
    Sham_mean_pct_positive = unname(Sham_mean_pct[eligible_genes]),
    Tx_mean_pct_positive = unname(Tx_mean_pct[eligible_genes]),
    allMphi_delta_pct_Tx_minus_Sham = unname(delta_pct_Tx_minus_Sham[eligible_genes]),
    expected_direction_signed_log2FC = unname(signed_log2fc[eligible_genes]),
    expected_direction_signed_delta_pct = unname(signed_delta_tx[eligible_genes]),
    pseudobulk_direction_consistency = pb_consistency,
    pct_positive_direction_consistency = pct_consistency,
    subtype_identity_score = identity_score,
    treatment_direction_score = treatment_score,
    final_staining_marker_score = final_score,
    technical_gene = grepl(TECHNICAL_GENE_PATTERN, eligible_genes),
    reference_gene = eligible_genes %in% reference_genes_now
  ) %>%
    arrange(desc(final_staining_marker_score)) %>%
    mutate(
      final_rank = row_number(),
      identity_rank = rank(-subtype_identity_score, ties.method = "first")
    )

  all_rankings[[subtype]] <- ranking

  write.csv(
    ranking,
    file.path(TAB_DIR, paste0("03_ALL_GENES_", slug, "_ranking_v5.8.1.csv")),
    row.names = FALSE
  )

  top200 <- ranking %>% filter(!technical_gene) %>% slice_head(n = TOP_N_TABLE)
  all_top200[[subtype]] <- top200
  write.csv(
    top200,
    file.path(TAB_DIR, paste0("04_TOP200_", slug, "_staining_candidates_v5.8.1.csv")),
    row.names = FALSE
  )

  reference_positions <- ranking %>% filter(reference_gene) %>% arrange(final_rank)
  all_reference[[subtype]] <- reference_positions
  write.csv(
    reference_positions,
    file.path(TAB_DIR, paste0("05_REFERENCE_GENES_", slug, "_positions_v5.8.1.csv")),
    row.names = FALSE
  )

  # Ranking plot
  rank_df <- top200 %>% slice_head(n = 40) %>% mutate(gene = factor(gene, levels = rev(gene)))
  p_rank <- ggplot(rank_df, aes(x = final_staining_marker_score, y = gene)) +
    geom_col() +
    labs(
      title = paste0(subtype, " | top staining-marker candidates"),
      subtitle = paste0("Expected Sham→Tx direction: ", toupper(direction_now)),
      x = "Final staining-marker score", y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  save_pdf(
    p_rank,
    file.path(FIG_DIR, paste0("01_TOP40_", slug, "_ranking_v5.8.1.pdf")),
    8.5, 10
  )

  # Component heatmap
  heat_df <- top200 %>%
    slice_head(n = 40) %>%
    transmute(
      gene,
      subtype_identity_score,
      UMAP_spatial_similarity = safe_rescale01(UMAP_spatial_similarity),
      subtype_delta_positive = safe_rescale01(pmax(delta_pct_target_minus_other, 0)),
      treatment_direction_score,
      pseudobulk_direction_consistency,
      pct_positive_direction_consistency,
      final_staining_marker_score
    )

  gene_order <- heat_df$gene
  heat_long <- heat_df %>% pivot_longer(-gene, names_to = "metric", values_to = "score")
  heat_long$gene <- factor(heat_long$gene, levels = rev(gene_order))

  p_heat <- ggplot(heat_long, aes(x = metric, y = gene, fill = score)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), oob = squish) +
    labs(title = paste0(subtype, " | ranking components"), x = NULL, y = NULL, fill = "Score") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(
    p_heat,
    file.path(FIG_DIR, paste0("02_TOP40_", slug, "_component_heatmap_v5.8.1.pdf")),
    10, 10
  )

  # Top 12 UMAPs
  top12 <- top200 %>% slice_head(n = TOP_N_UMAP) %>% pull(gene)
  p_umap <- FeaturePlot(
    obj, features = top12, reduction = UMAP_USE, ncol = 4,
    order = TRUE, min.cutoff = "q05", max.cutoff = "q98",
    raster = FALSE, pt.size = 0.28
  ) &
    scale_colour_gradientn(colours = c("#0033FF", "#FFFFFF", "#FF1A1A")) &
    theme_classic(base_size = 8)
  save_pdf(
    p_umap,
    file.path(UMAP_DIR, paste0("01_TOP12_", slug, "_FeaturePlot_v5.8.1.pdf")),
    14, 10
  )

  # Sham vs Tx top 6 UMAPs
  obj$ShamTx_v580 <- factor(
    ifelse(as.character(obj$condition_v580) %in% c("Sham", "Tx"),
           as.character(obj$condition_v580), NA_character_),
    levels = c("Sham", "Tx")
  )
  shamtx_cells <- colnames(obj)[!is.na(obj$ShamTx_v580)]
  obj_st <- subset(obj, cells = shamtx_cells)
  top6 <- top200 %>% slice_head(n = TOP_N_SHAMTX) %>% pull(gene)

  p_st <- FeaturePlot(
    obj_st, features = top6, reduction = UMAP_USE,
    split.by = "ShamTx_v580", keep.scale = "all", ncol = 2,
    order = TRUE, min.cutoff = "q05", max.cutoff = "q98",
    raster = FALSE, pt.size = 0.24
  ) &
    scale_colour_gradientn(colours = c("#0033FF", "#FFFFFF", "#FF1A1A")) &
    theme_classic(base_size = 8)
  save_pdf(
    p_st,
    file.path(UMAP_DIR, paste0("02_TOP6_", slug, "_Sham_vs_Tx_FeaturePlot_v5.8.1.pdf")),
    11, 3.4 * TOP_N_SHAMTX
  )
}

# ==============================================================================
# 6. Combined tables
# ==============================================================================

combined_ranking <- bind_rows(all_rankings)
combined_top200 <- bind_rows(all_top200)
combined_reference <- bind_rows(all_reference)

write.csv(
  combined_ranking,
  file.path(TAB_DIR, "06_COMBINED_ALL_GENES_four_subtype_ranking_v5.8.1.csv"),
  row.names = FALSE
)
write.csv(
  combined_top200,
  file.path(TAB_DIR, "07_COMBINED_TOP200_four_subtype_candidates_v5.8.1.csv"),
  row.names = FALSE
)
write.csv(
  combined_reference,
  file.path(TAB_DIR, "08_COMBINED_reference_gene_positions_v5.8.1.csv"),
  row.names = FALSE
)

# Cross-subtype comparison source
cross_subtype <- combined_ranking %>%
  select(subtype, gene, final_rank, identity_rank, subtype_identity_score,
         treatment_direction_score, final_staining_marker_score) %>%
  group_by(subtype) %>%
  slice_head(n = 500) %>%
  ungroup()
write.csv(
  cross_subtype,
  file.path(TAB_DIR, "09_TOP500_cross_subtype_specificity_source_v5.8.1.csv"),
  row.names = FALSE
)

# Reference panel UMAP for familiar markers
reference_present <- intersect(unique(unlist(REFERENCE_GENES, use.names = FALSE)), rownames(obj))
if (length(reference_present)) {
  p_ref <- FeaturePlot(
    obj, features = reference_present, reduction = UMAP_USE, ncol = 4,
    order = TRUE, min.cutoff = "q05", max.cutoff = "q98",
    raster = FALSE, pt.size = 0.25
  ) &
    scale_colour_gradientn(colours = c("#0033FF", "#FFFFFF", "#FF1A1A")) &
    theme_classic(base_size = 7)
  save_pdf(
    p_ref,
    file.path(UMAP_DIR, "03_REFERENCE_marker_panel_FeaturePlot_v5.8.1.pdf"),
    15, 3.0 * ceiling(length(reference_present) / 4)
  )
}

# ==============================================================================
# 7. Metadata / terminal summary
# ==============================================================================

metadata_out <- tibble(
  parameter = c(
    "version", "input_RDS", "class_column", "condition_column", "sample_column",
    "cluster_column", "UMAP", "assay", "weight_final_identity",
    "weight_final_treatment", "grid_nx", "grid_ny", "min_cells_per_grid_bin",
    "min_detection_fraction"
  ),
  value = c(
    "v5.8.1", INPUT_RDS, CLASS_COLUMN, CONDITION_COLUMN, SAMPLE_COLUMN,
    CLUSTER_COLUMN, UMAP_USE, ASSAY_USE, W_FINAL_IDENTITY,
    W_FINAL_TREATMENT, GRID_NX, GRID_NY, MIN_CELLS_PER_GRID_BIN,
    MIN_DETECTION_FRACTION
  )
)
write.csv(metadata_out, file.path(LOG_DIR, "analysis_metadata_v5.8.1.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(LOG_DIR, "sessionInfo_v5.8.1.txt"))

cat("\n============================================================\n")
cat("Four macrophage-subtype staining-marker discovery v5.8.1\n")
cat("============================================================\n")

for (subtype in TARGET_CLASSES) {
  cat("\n### ", subtype,
      " | expected Tx direction: ", EXPECTED_DIRECTION[[subtype]], "\n", sep = "")
  print(
    all_rankings[[subtype]] %>%
      filter(!technical_gene) %>%
      select(
        final_rank, identity_rank, gene, final_staining_marker_score,
        subtype_identity_score, treatment_direction_score,
        pct_positive_target, pct_positive_other_Mphi,
        delta_pct_target_minus_other, expression_ratio_target_vs_other,
        UMAP_spatial_similarity, allMphi_log2FC_Tx_vs_Sham,
        allMphi_delta_pct_Tx_minus_Sham,
        pseudobulk_direction_consistency, pct_positive_direction_consistency,
        reference_gene
      ) %>%
      slice_head(n = 25)
  )
}

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
