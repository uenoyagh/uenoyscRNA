# ============================================================
# Human MASH v1.3.2A.1
# Immune Manual Annotation & Cleaning
# ============================================================
# Purpose
#   - Start from the v1.3.2A pre-plot checkpoint.
#   - Keep the existing RPCA embedding and res0.6 clusters.
#   - Apply an explicit, auditable manual cluster map.
#   - Remove non-immune contaminant clusters.
#   - Generate clean immune and clean mononuclear-myeloid RDS files.
#   - Re-evaluate the IL10 / IL10R axis after cleaning.
#
# Important
#   - No RPCA reintegration.
#   - No reclustering.
#   - No FindAllMarkers rerun.
#   - Canonical whole-liver v1.3 RDS is never modified.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(scales)
})

options(stringsAsFactors = FALSE)
seed <- 20260908
set.seed(seed)

# ------------------------------------------------------------
# 0. Paths / constants
# ------------------------------------------------------------
parent_outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2A_immune_refinement"
)

input_rds <- file.path(
  parent_outdir,
  "objects",
  "Human_MASH_Immune_v1.3.2A_preplot_checkpoint.rds"
)

marker_csv <- file.path(
  parent_outdir,
  "13_immune_cluster_markers_top15.csv"
)

provisional_csv <- file.path(
  parent_outdir,
  "06_cluster_provisional_annotation.csv"
)

outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2A.1_immune_manual_annotation_cleaning"
)

object_dir <- file.path(outdir, "objects")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

sample_col <- "canonical_sample_v13"
nas_col <- "NAS_score_v13"
nasgrp_col <- "NAS_group_v13"
cluster_col <- "immune_cluster_v132A"
reduction_name <- "umapImmuneRPCA_v132A"

expected_n_cells <- 14833L
expected_clusters <- as.character(0:21)

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------
sample_order_fun <- function(x) {
  u <- unique(as.character(x))
  if (length(u) > 0 && all(grepl("^LT[0-9]+", u))) {
    n <- suppressWarnings(as.integer(sub("^LT([0-9]+).*", "\\1", u)))
    if (!anyNA(n)) return(u[order(n)])
  }
  sort(u)
}

safe_present <- function(genes, universe) {
  intersect(unique(genes), universe)
}

dedupe_feature_list <- function(x) {
  seen <- character(0)
  out <- lapply(x, function(g) {
    g <- unique(g)
    keep <- g[!g %in% seen]
    seen <<- c(seen, keep)
    keep
  })
  out[lengths(out) > 0]
}

# ------------------------------------------------------------
# 2. Input validation
# ------------------------------------------------------------
cat("============================================================\n")
cat("Human MASH v1.3.2A.1 Immune Manual Annotation & Cleaning\n")
cat("============================================================\n\n")

required_files <- c(input_rds, marker_csv, provisional_csv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

immune_all <- readRDS(input_rds)
DefaultAssay(immune_all) <- "RNA"

required_meta <- c(
  sample_col,
  nas_col,
  nasgrp_col,
  cluster_col,
  "current_celltype_v132A",
  "immune_entry_reason_v132A",
  "immune_celltype_provisional_v132A"
)
missing_meta <- setdiff(required_meta, colnames(immune_all@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

if (!reduction_name %in% Reductions(immune_all)) {
  stop("Required reduction missing: ", reduction_name)
}

if (ncol(immune_all) != expected_n_cells) {
  stop(
    "Cell-count mismatch. Expected ", expected_n_cells,
    " but found ", ncol(immune_all),
    ". Manual map NOT applied."
  )
}

cluster_vec <- as.character(immune_all[[cluster_col]][, 1])
observed_clusters <- unique(cluster_vec)
if (!setequal(observed_clusters, expected_clusters)) {
  stop(
    "Cluster-set mismatch. Expected: ", paste(expected_clusters, collapse = ","),
    " | Observed: ", paste(sort(observed_clusters), collapse = ","),
    ". Manual map NOT applied."
  )
}

cat("Input cells    :", ncol(immune_all), "\n")
cat("Input clusters :", length(observed_clusters), "\n")
cat("Reduction      :", reduction_name, "\n\n")

# ------------------------------------------------------------
# 3. Previous marker / provisional tables
# ------------------------------------------------------------
marker_tbl <- read.csv(marker_csv, check.names = FALSE)
provisional_tbl <- read.csv(provisional_csv, check.names = FALSE)

if (!all(c("cluster", "gene") %in% colnames(marker_tbl))) {
  stop("Marker CSV lacks cluster/gene columns.")
}
if (!all(c("cluster", "provisional_celltype") %in% colnames(provisional_tbl))) {
  stop("Provisional annotation CSV lacks required columns.")
}

marker_tbl$cluster <- as.character(marker_tbl$cluster)
provisional_tbl$cluster <- as.character(provisional_tbl$cluster)

top_marker_summary <- aggregate(
  gene ~ cluster,
  data = marker_tbl,
  FUN = function(x) paste(unique(as.character(x)), collapse = ";")
)
names(top_marker_summary)[2] <- "observed_top15_markers"

# ------------------------------------------------------------
# 4. Explicit manual cluster map
# ------------------------------------------------------------
manual_map <- data.frame(
  cluster = as.character(0:21),
  manual_celltype = c(
    "Resident_macrophage_TIMD4_FOLR2", # 0
    "cDC2_DC3_like",                    # 1
    "T_IL7R_TCF7_like",                 # 2
    "T_CD8_GZMK",                       # 3
    "Monocyte_inflammatory",            # 4
    "Resident_macrophage_state2",       # 5
    "NK_cytotoxic",                     # 6
    "SPP1_TREM2_macrophage",            # 7
    "NK_gamma_delta_T_like",            # 8
    "Exclude_hepatocyte_like",          # 9
    "Exclude_hepatocyte_like",          # 10
    "Plasma_cell",                      # 11
    "Cycling_immune",                   # 12
    "B_cell",                           # 13
    "Exclude_Schwann_neural_like",      # 14
    "Mast_cell",                        # 15
    "Neutrophil",                       # 16
    "Exclude_perivascular_mesenchymal", # 17
    "Exclude_endothelial",              # 18
    "Exclude_erythroid",                # 19
    "cDC1",                             # 20
    "Exclude_ductal_like"               # 21
  ),
  broad_lineage = c(
    "Macrophage", "Dendritic", "T_cell", "T_cell", "Monocyte",
    "Macrophage", "NK_cell", "Macrophage", "NK_T",
    "Excluded_nonimmune", "Excluded_nonimmune", "Plasma_cell",
    "Cycling_immune", "B_cell", "Excluded_nonimmune", "Mast_cell",
    "Neutrophil", "Excluded_nonimmune", "Excluded_nonimmune",
    "Excluded_nonimmune", "Dendritic", "Excluded_nonimmune"
  ),
  keep_clean = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE,
    FALSE, FALSE, FALSE, TRUE, FALSE
  ),
  manual_confidence = c(
    "high", "high", "high", "high", "high", "moderate", "high",
    "high", "moderate", "high", "high", "high", "high", "high",
    "high", "high", "high", "high", "high", "high", "high", "high"
  ),
  evidence_basis = c(
    "MARCO;TIMD4;CD5L;FOLR2",
    "CLEC10A;CLEC5A;F13A1;STAB1;FCGR2B",
    "IL7R;TCF7;LTB;CD28;RORC",
    "CD8A;GZMK;CCL5;TIGIT;LAG3",
    "VCAN;S100A8;S100A9;FCAR;LYZ",
    "VCAM1;LILRB5;SLC16A9;NDST3;TCF7L1",
    "GNLY;GZMB;PRF1;FGFBP2;S1PR5",
    "SPP1;TREM2;LPL;FABP5;MMP19",
    "NCR1;KLRF1;TRDC;IL2RB;KLRC1",
    "APOA1;AMBP;ORM2;SERPINA1",
    "CYP3A4;SLCO1B1;GCGR;VIL1",
    "JCHAIN;IGHA1;DERL3;IGHG1",
    "cell-cycle program",
    "MS4A1;CD19;CD79A;PAX5;CD22",
    "SOX10;PLP1;GPM6B;RXRG",
    "KIT;HDC;MS4A2;CPA3;TPSG1",
    "CEACAM8;ELANE;MPO;PRTN3;DEFA3",
    "MYH11;PTH1R;CACNA1C;BMP5",
    "KDR;EMCN;NOS3;GPR182",
    "ALAS2;HBB;HBD;SLC4A1",
    "CLEC9A;XCR1;IDO1",
    "KRT7;PROM1;SLC34A2;VTCN1"
  ),
  exclusion_reason = c(
    NA, NA, NA, NA, NA, NA, NA, NA, NA,
    "hepatocyte_like", "hepatocyte_like", NA, NA, NA,
    "Schwann_neural_like", NA, NA,
    "perivascular_mesenchymal_like", "endothelial", "erythroid",
    NA, "ductal_like"
  ),
  stringsAsFactors = FALSE
)

if (!setequal(manual_map$cluster, expected_clusters)) {
  stop("Internal manual-map cluster mismatch.")
}

manual_map$provisional_v132A <- provisional_tbl$provisional_celltype[
  match(manual_map$cluster, provisional_tbl$cluster)
]
manual_map$observed_top15_markers <- top_marker_summary$observed_top15_markers[
  match(manual_map$cluster, top_marker_summary$cluster)
]
manual_map <- manual_map[order(as.integer(manual_map$cluster)), ]

write.csv(
  manual_map,
  file.path(outdir, "01_manual_cluster_annotation_map.csv"),
  row.names = FALSE
)

cat("=== MANUAL CLUSTER MAP ===\n")
print(
  manual_map[, c("cluster", "provisional_v132A", "manual_celltype", "keep_clean", "manual_confidence")],
  row.names = FALSE
)
cat("\n")

# ------------------------------------------------------------
# 5. Apply manual annotation
# ------------------------------------------------------------
type_map <- setNames(manual_map$manual_celltype, manual_map$cluster)
broad_map <- setNames(manual_map$broad_lineage, manual_map$cluster)
keep_map <- setNames(manual_map$keep_clean, manual_map$cluster)
confidence_map <- setNames(manual_map$manual_confidence, manual_map$cluster)
exclusion_map <- setNames(manual_map$exclusion_reason, manual_map$cluster)

immune_all$immune_celltype_manual_v132A1 <- unname(type_map[cluster_vec])
immune_all$immune_broadlineage_manual_v132A1 <- unname(broad_map[cluster_vec])
immune_all$immune_keep_clean_v132A1 <- as.logical(unname(keep_map[cluster_vec]))
immune_all$immune_manual_confidence_v132A1 <- unname(confidence_map[cluster_vec])
immune_all$immune_exclusion_reason_v132A1 <- unname(exclusion_map[cluster_vec])
immune_all$immune_status_v132A1 <- ifelse(
  immune_all$immune_keep_clean_v132A1,
  "clean_immune",
  "excluded_nonimmune"
)

if (anyNA(immune_all$immune_celltype_manual_v132A1)) {
  stop("NA introduced during manual annotation.")
}

# ------------------------------------------------------------
# 6. Cleaning audit
# ------------------------------------------------------------
cluster_counts <- as.data.frame(
  table(
    cluster = immune_all$immune_cluster_v132A,
    manual_celltype = immune_all$immune_celltype_manual_v132A1,
    keep_clean = immune_all$immune_keep_clean_v132A1
  ),
  stringsAsFactors = FALSE
)
cluster_counts <- cluster_counts[cluster_counts$Freq > 0, ]
write.csv(cluster_counts, file.path(outdir, "02_cluster_manual_annotation_counts.csv"), row.names = FALSE)

n_total <- ncol(immune_all)
n_clean <- sum(immune_all$immune_keep_clean_v132A1)
n_excluded <- n_total - n_clean

cleaning_summary <- data.frame(
  metric = c(
    "input_immune_candidate_cells",
    "clean_immune_cells",
    "excluded_nonimmune_cells",
    "clean_percent",
    "excluded_percent"
  ),
  value = c(
    n_total,
    n_clean,
    n_excluded,
    100 * n_clean / n_total,
    100 * n_excluded / n_total
  ),
  stringsAsFactors = FALSE
)
write.csv(cleaning_summary, file.path(outdir, "03_cleaning_summary.csv"), row.names = FALSE)

entry_audit <- as.data.frame(
  table(
    entry_reason = immune_all$immune_entry_reason_v132A,
    status = immune_all$immune_status_v132A1
  ),
  stringsAsFactors = FALSE
)
write.csv(entry_audit, file.path(outdir, "04_entry_reason_vs_cleaning_status.csv"), row.names = FALSE)

excluded_audit <- as.data.frame(
  table(
    cluster = immune_all$immune_cluster_v132A,
    excluded_type = immune_all$immune_celltype_manual_v132A1,
    sample = immune_all[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)
excluded_clusters <- manual_map$cluster[!manual_map$keep_clean]
excluded_audit <- excluded_audit[
  excluded_audit$cluster %in% excluded_clusters & excluded_audit$Freq > 0,
]
write.csv(excluded_audit, file.path(outdir, "05_excluded_nonimmune_by_cluster_sample.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 7. Clean immune and clean mononuclear-myeloid objects
# ------------------------------------------------------------
clean_cells <- colnames(immune_all)[immune_all$immune_keep_clean_v132A1]
immune_clean <- subset(immune_all, cells = clean_cells)
DefaultAssay(immune_clean) <- "RNA"

mononuclear_myeloid_types <- c(
  "Resident_macrophage_TIMD4_FOLR2",
  "Resident_macrophage_state2",
  "SPP1_TREM2_macrophage",
  "Monocyte_inflammatory",
  "cDC2_DC3_like",
  "cDC1"
)
myeloid_cells <- colnames(immune_clean)[
  immune_clean$immune_celltype_manual_v132A1 %in% mononuclear_myeloid_types
]
myeloid_clean <- subset(immune_clean, cells = myeloid_cells)
DefaultAssay(myeloid_clean) <- "RNA"

cat("Clean immune cells :", ncol(immune_clean), "\n")
cat("Excluded cells     :", n_excluded, "\n")
cat("Myeloid cells      :", ncol(myeloid_clean), "\n\n")

# ------------------------------------------------------------
# 8. Counts / sample composition / dominance
# ------------------------------------------------------------
manual_counts <- as.data.frame(
  table(celltype = immune_clean$immune_celltype_manual_v132A1),
  stringsAsFactors = FALSE
)
names(manual_counts)[2] <- "n_cells"
manual_counts$percent <- 100 * manual_counts$n_cells / sum(manual_counts$n_cells)
manual_counts <- manual_counts[order(manual_counts$n_cells, decreasing = TRUE), ]
write.csv(manual_counts, file.path(outdir, "06_clean_manual_celltype_counts.csv"), row.names = FALSE)

sample_vec <- as.character(immune_clean[[sample_col]][, 1])
celltype_vec <- as.character(immune_clean$immune_celltype_manual_v132A1)
sample_order <- sample_order_fun(sample_vec)

type_sample <- as.data.frame(
  table(celltype = celltype_vec, sample = sample_vec),
  stringsAsFactors = FALSE
)
type_sample$celltype_total <- ave(type_sample$Freq, type_sample$celltype, FUN = sum)
type_sample$sample_total <- ave(type_sample$Freq, type_sample$sample, FUN = sum)
type_sample$within_celltype_pct <- ifelse(
  type_sample$celltype_total > 0,
  100 * type_sample$Freq / type_sample$celltype_total,
  NA_real_
)
type_sample$within_sample_pct <- ifelse(
  type_sample$sample_total > 0,
  100 * type_sample$Freq / type_sample$sample_total,
  NA_real_
)
type_sample$sample <- factor(type_sample$sample, levels = sample_order)
write.csv(type_sample, file.path(outdir, "07_clean_manual_celltype_by_sample.csv"), row.names = FALSE)

split_types <- split(seq_len(nrow(type_sample)), type_sample$celltype)
max_idx <- unlist(lapply(split_types, function(idx) {
  idx[which.max(type_sample$within_celltype_pct[idx])]
}))
sample_dominance <- type_sample[max_idx, c("celltype", "sample", "within_celltype_pct")]
names(sample_dominance)[3] <- "max_sample_contribution_pct"
sample_dominance$sample <- as.character(sample_dominance$sample)
sample_dominance <- sample_dominance[
  order(sample_dominance$max_sample_contribution_pct, decreasing = TRUE),
]
write.csv(sample_dominance, file.path(outdir, "08_clean_manual_celltype_sample_dominance.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 9. Crosswalks
# ------------------------------------------------------------
crosswalk_current <- as.data.frame(
  table(
    current_REF12 = immune_all$current_celltype_v132A,
    manual_v132A1 = immune_all$immune_celltype_manual_v132A1
  ),
  stringsAsFactors = FALSE
)
write.csv(crosswalk_current, file.path(outdir, "09_current_REF12_vs_manual_crosswalk.csv"), row.names = FALSE)

crosswalk_provisional <- as.data.frame(
  table(
    provisional_v132A = immune_all$immune_celltype_provisional_v132A,
    manual_v132A1 = immune_all$immune_celltype_manual_v132A1
  ),
  stringsAsFactors = FALSE
)
write.csv(crosswalk_provisional, file.path(outdir, "10_provisional_vs_manual_crosswalk.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 10. Palettes
# ------------------------------------------------------------
manual_palette <- c(
  "Resident_macrophage_TIMD4_FOLR2" = "#8B0000",
  "Resident_macrophage_state2" = "#C62828",
  "SPP1_TREM2_macrophage" = "#FF7F00",
  "Monocyte_inflammatory" = "#E31A1C",
  "cDC2_DC3_like" = "#00A6A6",
  "cDC1" = "#00CFEA",
  "T_IL7R_TCF7_like" = "#6A3D9A",
  "T_CD8_GZMK" = "#9C27B0",
  "NK_cytotoxic" = "#1F78B4",
  "NK_gamma_delta_T_like" = "#1976D2",
  "B_cell" = "#00A650",
  "Plasma_cell" = "#B5179E",
  "Mast_cell" = "#FF1493",
  "Neutrophil" = "#7CB342",
  "Cycling_immune" = "#FFD600",
  "Exclude_hepatocyte_like" = "#BDBDBD",
  "Exclude_Schwann_neural_like" = "#9E9E9E",
  "Exclude_perivascular_mesenchymal" = "#757575",
  "Exclude_endothelial" = "#90A4AE",
  "Exclude_erythroid" = "#8D6E63",
  "Exclude_ductal_like" = "#A5D6A7"
)

# ------------------------------------------------------------
# 11. UMAP figures
# ------------------------------------------------------------
p01 <- DimPlot(
  immune_all,
  reduction = reduction_name,
  group.by = "immune_celltype_manual_v132A1",
  pt.size = 0.18,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols = manual_palette
) +
  ggtitle("v1.3.2A.1 manual annotation - all immune candidates") +
  theme_classic(base_size = 11)

ggsave(file.path(outdir, "Fig01_ManualAnnotation_AllCandidates_UMAP.pdf"), p01, width = 13, height = 9)

p02 <- DimPlot(
  immune_clean,
  reduction = reduction_name,
  group.by = "immune_celltype_manual_v132A1",
  pt.size = 0.20,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols = manual_palette
) +
  ggtitle("v1.3.2A.1 clean immune annotation") +
  theme_classic(base_size = 11)

ggsave(file.path(outdir, "Fig02_CleanImmune_ManualAnnotation_UMAP.pdf"), p02, width = 12, height = 9)

excluded_cells <- colnames(immune_all)[!immune_all$immune_keep_clean_v132A1]
excluded_obj <- subset(immune_all, cells = excluded_cells)
DefaultAssay(excluded_obj) <- "RNA"

p03 <- DimPlot(
  excluded_obj,
  reduction = reduction_name,
  group.by = "immune_celltype_manual_v132A1",
  pt.size = 0.30,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols = manual_palette
) +
  ggtitle("Excluded non-immune / contaminant clusters") +
  theme_classic(base_size = 11)

ggsave(file.path(outdir, "Fig03_ExcludedNonimmune_UMAP.pdf"), p03, width = 11, height = 8)

p04 <- DimPlot(
  immune_clean,
  reduction = reduction_name,
  group.by = "immune_celltype_manual_v132A1",
  split.by = sample_col,
  ncol = 4,
  pt.size = 0.14,
  raster = TRUE,
  cols = manual_palette
)
ggsave(file.path(outdir, "Fig04_CleanImmune_ManualAnnotation_split_by_LT.pdf"), p04, width = 19, height = 13)

# ------------------------------------------------------------
# 12. Marker panel for clean immune types
# ------------------------------------------------------------
marker_list <- list(
  Resident_macrophage_TIMD4_FOLR2 = c("MARCO", "TIMD4", "CD5L", "FOLR2", "C1QC"),
  Resident_macrophage_state2 = c("VCAM1", "LILRB5", "SLC16A9", "NDST3", "TCF7L1"),
  SPP1_TREM2_macrophage = c("SPP1", "TREM2", "LPL", "FABP5", "MMP19"),
  Monocyte_inflammatory = c("VCAN", "S100A8", "S100A9", "FCAR", "LYZ"),
  cDC2_DC3_like = c("CLEC10A", "CLEC5A", "F13A1", "STAB1", "FCGR2B"),
  cDC1 = c("CLEC9A", "XCR1", "IDO1", "CLNK", "EGLN3"),
  T_IL7R_TCF7_like = c("IL7R", "TCF7", "LTB", "CD28", "RORC"),
  T_CD8_GZMK = c("CD8A", "GZMK", "CCL5", "TIGIT", "LAG3"),
  NK_cytotoxic = c("GNLY", "GZMB", "PRF1", "FGFBP2", "S1PR5"),
  NK_gamma_delta_T_like = c("TRDC", "NCR1", "KLRF1", "IL2RB", "KLRC1"),
  B_cell = c("MS4A1", "CD19", "CD79A", "PAX5", "CD22"),
  Plasma_cell = c("JCHAIN", "DERL3", "IGHA1", "IGHG1", "POU2AF1"),
  Mast_cell = c("KIT", "HDC", "MS4A2", "CPA3", "TPSG1"),
  Neutrophil = c("CEACAM8", "ELANE", "MPO", "PRTN3", "DEFA3"),
  Cycling_immune = c("MKI67", "TOP2A", "UBE2C", "CENPF", "RRM2")
)

marker_universe <- rownames(immune_clean)
marker_availability <- do.call(rbind, lapply(names(marker_list), function(x) {
  genes <- marker_list[[x]]
  data.frame(
    marker_group = x,
    gene = genes,
    present = genes %in% marker_universe,
    stringsAsFactors = FALSE
  )
}))
rownames(marker_availability) <- NULL
write.csv(marker_availability, file.path(outdir, "11_manual_marker_gene_availability.csv"), row.names = FALSE)

marker_plot_list <- lapply(marker_list, function(x) safe_present(x, marker_universe))
marker_plot_list <- dedupe_feature_list(marker_plot_list)
stopifnot(!anyDuplicated(unlist(marker_plot_list, use.names = FALSE)))

p05 <- DotPlot(
  immune_clean,
  features = marker_plot_list,
  assay = "RNA",
  group.by = "immune_celltype_manual_v132A1",
  dot.scale = 7
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  labs(
    title = "v1.3.2A.1 clean immune manual annotation validation",
    color = "Scaled\nexpression",
    size = "% expressed"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

ggsave(file.path(outdir, "Fig05_CleanImmune_ManualMarker_DotPlot.pdf"), p05, width = 24, height = 10)
write.csv(p05$data, file.path(outdir, "12_CleanImmune_ManualMarker_DotPlot_data.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 13. Excluded-cluster marker validation
# ------------------------------------------------------------
excluded_marker_list <- list(
  Exclude_hepatocyte_like = c("APOA1", "AMBP", "CYP3A4", "SLCO1B1", "SERPINA1"),
  Exclude_Schwann_neural_like = c("SOX10", "PLP1", "GPM6B", "RXRG"),
  Exclude_perivascular_mesenchymal = c("MYH11", "PTH1R", "CACNA1C", "BMP5"),
  Exclude_endothelial = c("KDR", "EMCN", "NOS3", "GPR182"),
  Exclude_erythroid = c("ALAS2", "HBB", "HBD", "SLC4A1"),
  Exclude_ductal_like = c("KRT7", "PROM1", "SLC34A2", "VTCN1")
)
excluded_marker_list <- lapply(excluded_marker_list, function(x) safe_present(x, rownames(excluded_obj)))
excluded_marker_list <- dedupe_feature_list(excluded_marker_list)

if (length(excluded_marker_list) > 0) {
  p06 <- DotPlot(
    excluded_obj,
    features = excluded_marker_list,
    assay = "RNA",
    group.by = "immune_celltype_manual_v132A1",
    dot.scale = 7
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    labs(
      title = "Excluded non-immune cluster marker validation",
      color = "Scaled\nexpression",
      size = "% expressed"
    ) +
    theme_bw(base_size = 10)

  ggsave(file.path(outdir, "Fig06_ExcludedNonimmune_Marker_DotPlot.pdf"), p06, width = 16, height = 7)
}

# ------------------------------------------------------------
# 14. Sample QC / composition
# ------------------------------------------------------------
p07 <- ggplot(
  type_sample,
  aes(x = sample, y = celltype, fill = within_celltype_pct)
) +
  geom_tile(color = "grey80", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(Freq > 0, sprintf("%.1f", within_celltype_pct), "")),
    size = 2.3
  ) +
  scale_fill_gradient(low = "#FFFFFF", high = "#FF1A1A") +
  labs(
    title = "Sample contribution within clean manual immune types",
    x = NULL,
    y = NULL,
    fill = "%"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outdir, "Fig07_CleanImmune_SampleContribution.pdf"), p07, width = 13, height = 9)

p08 <- ggplot(
  type_sample,
  aes(x = sample, y = within_sample_pct, fill = celltype)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = manual_palette) +
  labs(
    title = "Clean immune composition by LT sample",
    x = NULL,
    y = "% of clean immune cells within sample",
    fill = "Manual immune type"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outdir, "Fig08_CleanImmune_Composition_by_LT.pdf"), p08, width = 13, height = 9)

# ------------------------------------------------------------
# 15. IL10-axis outputs
# ------------------------------------------------------------
il10_genes <- c("IL10", "IL10RA", "IL10RB", "JAK1", "TYK2", "STAT3", "SOCS3")
il10_present <- safe_present(il10_genes, rownames(immune_clean))

il10_availability <- data.frame(
  gene = il10_genes,
  present = il10_genes %in% rownames(immune_clean),
  stringsAsFactors = FALSE
)
write.csv(il10_availability, file.path(outdir, "13_IL10_axis_gene_availability.csv"), row.names = FALSE)

if (length(il10_present) > 0) {
  p09 <- DotPlot(
    immune_clean,
    features = il10_present,
    assay = "RNA",
    group.by = "immune_celltype_manual_v132A1",
    dot.scale = 8
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    labs(
      title = "IL10 / IL10R axis after manual immune cleaning",
      subtitle = "Annotation review required before mechanistic freeze",
      color = "Scaled\nexpression",
      size = "% expressed"
    ) +
    theme_bw(base_size = 10)

  ggsave(file.path(outdir, "Fig09_CleanImmune_IL10axis_DotPlot.pdf"), p09, width = 12, height = 8)
  write.csv(p09$data, file.path(outdir, "14_CleanImmune_IL10axis_DotPlot_data.csv"), row.names = FALSE)
}

il10_myeloid_present <- safe_present(il10_genes, rownames(myeloid_clean))
if (length(il10_myeloid_present) > 0) {
  p10 <- DotPlot(
    myeloid_clean,
    features = il10_myeloid_present,
    assay = "RNA",
    group.by = "immune_celltype_manual_v132A1",
    dot.scale = 9
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    labs(
      title = "IL10 / IL10R axis - clean mononuclear myeloid compartment",
      subtitle = "Macrophage / monocyte / DC compartment",
      color = "Scaled\nexpression",
      size = "% expressed"
    ) +
    theme_bw(base_size = 10)

  ggsave(file.path(outdir, "Fig10_MononuclearMyeloid_IL10axis_DotPlot.pdf"), p10, width = 11, height = 7)
  write.csv(p10$data, file.path(outdir, "15_MononuclearMyeloid_IL10axis_DotPlot_data.csv"), row.names = FALSE)
}

feature_genes <- safe_present(c("IL10RA", "IL10RB", "STAT3", "SOCS3"), rownames(myeloid_clean))
if (length(feature_genes) > 0) {
  p11 <- FeaturePlot(
    myeloid_clean,
    features = feature_genes,
    reduction = reduction_name,
    ncol = 2,
    pt.size = 0.18,
    raster = TRUE,
    order = FALSE,
    cols = c("#F2F2F2", "#FF1A1A")
  )
  ggsave(file.path(outdir, "Fig11_MononuclearMyeloid_IL10R_STAT3_SOCS3_FeaturePlots.pdf"), p11, width = 11, height = 9)
}

# ------------------------------------------------------------
# 16. Per-cell metadata exports
# ------------------------------------------------------------
metadata_all <- data.frame(
  cell = colnames(immune_all),
  sample = immune_all[[sample_col]][, 1],
  NAS_score = immune_all[[nas_col]][, 1],
  NAS_group = immune_all[[nasgrp_col]][, 1],
  cluster_v132A = immune_all$immune_cluster_v132A,
  current_celltype_v132A = immune_all$current_celltype_v132A,
  provisional_v132A = immune_all$immune_celltype_provisional_v132A,
  manual_v132A1 = immune_all$immune_celltype_manual_v132A1,
  broad_lineage_v132A1 = immune_all$immune_broadlineage_manual_v132A1,
  keep_clean_v132A1 = immune_all$immune_keep_clean_v132A1,
  exclusion_reason_v132A1 = immune_all$immune_exclusion_reason_v132A1,
  entry_reason_v132A = immune_all$immune_entry_reason_v132A,
  stringsAsFactors = FALSE
)
write.csv(metadata_all, file.path(outdir, "16_v1.3.2A.1_all_candidate_cell_metadata.csv"), row.names = FALSE)

metadata_clean <- metadata_all[metadata_all$keep_clean_v132A1, ]
write.csv(metadata_clean, file.path(outdir, "17_v1.3.2A.1_clean_immune_cell_metadata.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 17. Review flags
# ------------------------------------------------------------
flags <- c(
  "Human MASH v1.3.2A.1 Immune Manual Annotation & Cleaning",
  "=========================================================",
  "",
  paste0("Input immune-candidate cells: ", n_total),
  paste0("Clean immune cells: ", n_clean, " (", round(100 * n_clean / n_total, 2), "%)"),
  paste0("Excluded non-immune cells: ", n_excluded, " (", round(100 * n_excluded / n_total, 2), "%)"),
  paste0("Clean manual immune types: ", length(unique(immune_clean$immune_celltype_manual_v132A1))),
  "",
  "MANUAL REVIEW NOTES:",
  "  Cluster 5 = Resident_macrophage_state2; retain, but review during macrophage-state analysis.",
  "  Cluster 8 = NK_gamma_delta_T_like; lineage remains intentionally cautious.",
  "  Cluster 12 = Cycling_immune; cycling state is clear, parent lineage is not assigned here.",
  "  No pDC population is frozen at v1.3.2A.1.",
  "",
  "EXCLUDED CLUSTERS:",
  "  9,10 = hepatocyte-like",
  "  14   = Schwann/neural-like",
  "  17   = perivascular/mesenchymal-like",
  "  18   = endothelial",
  "  19   = erythroid",
  "  21   = ductal/cholangiocyte-like",
  "",
  "DOWNSTREAM:",
  "  Clean mononuclear-myeloid object generated for future IL10 analysis.",
  "  IL10-axis plots remain validation outputs until annotation freeze.",
  "",
  "STATUS:",
  "  MANUAL_CLEANED_REVIEW_REQUIRED",
  "  Canonical whole-liver v1.3 object was NOT modified."
)

dominant_types <- sample_dominance[sample_dominance$max_sample_contribution_pct >= 70, ]
if (nrow(dominant_types) > 0) {
  flags <- c(
    flags,
    "",
    "FLAG: manual immune types with >=70% contribution from one sample:",
    paste0(
      "  ", dominant_types$celltype,
      " : ", dominant_types$sample,
      " = ", round(dominant_types$max_sample_contribution_pct, 1), "%"
    )
  )
}

writeLines(flags, file.path(outdir, "18_v1.3.2A.1_review_flags.txt"))
cat(paste(flags, collapse = "\n"), "\n\n")

# ------------------------------------------------------------
# 18. Save objects
# ------------------------------------------------------------
all_rds <- file.path(object_dir, "Human_MASH_Immune_v1.3.2A.1_manual_all.rds")
clean_rds <- file.path(object_dir, "Human_MASH_Immune_v1.3.2A.1_clean.rds")
myeloid_rds <- file.path(object_dir, "Human_MASH_MononuclearMyeloid_v1.3.2A.1_clean.rds")

saveRDS(immune_all, all_rds)
saveRDS(immune_clean, clean_rds)
saveRDS(myeloid_clean, myeloid_rds)

# ------------------------------------------------------------
# 19. Manifest / provenance
# ------------------------------------------------------------
input_md5 <- unname(tools::md5sum(input_rds))
manifest <- c(
  "version=v1.3.2A.1",
  paste0("input=", input_rds),
  paste0("input_md5=", input_md5),
  paste0("input_cells=", n_total),
  paste0("clean_immune_cells=", n_clean),
  paste0("excluded_nonimmune_cells=", n_excluded),
  "clustering_source=v1.3.2A_res0.6",
  "RPCA_recomputed=FALSE",
  "clustering_recomputed=FALSE",
  "FindAllMarkers_recomputed=FALSE",
  paste0("manual_all_object=", all_rds),
  paste0("clean_immune_object=", clean_rds),
  paste0("clean_mononuclear_myeloid_object=", myeloid_rds),
  "canonical_object_modified=FALSE",
  "annotation_status=MANUAL_CLEANED_REVIEW_REQUIRED"
)
writeLines(manifest, file.path(outdir, "98_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "99_sessionInfo.txt"))

cat("============================================================\n")
cat("v1.3.2A.1 COMPLETE\n")
cat("MANUAL IMMUNE CLEANING COMPLETE\n")
cat("ANNOTATION STATUS: REVIEW REQUIRED\n")
cat("CANONICAL v1.3 OBJECT NOT MODIFIED\n")
cat("Clean immune object:\n", clean_rds, "\n")
cat("Clean mononuclear-myeloid object:\n", myeloid_rds, "\n")
cat("Output directory:\n", outdir, "\n")
cat("============================================================\n")
