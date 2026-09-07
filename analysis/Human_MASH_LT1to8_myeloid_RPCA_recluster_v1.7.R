#!/usr/bin/env Rscript

# Human_MASH_LT1to8_myeloid_RPCA_recluster_v1.7.R
#
# PURPOSE
# -------
# Sample-aware RPCA reclustering of the 6,956-cell clean LT1-LT8 myeloid core
# created in v1.6.
#
# INPUT:
#   Human_MASH_LT1to8_MyeloidClean_v1.6.rds
#
# THIS SCRIPT DOES:
#   - create a fresh analysis object from RAW RNA counts
#   - normalize within each LT sample
#   - identify variable features
#   - PCA
#   - RPCA integration across LT1-LT8
#   - neighbors + UMAP
#   - cluster scan at res 0.2, 0.4, 0.6, 0.8, 1.0, 1.2
#   - sample-balance QC
#   - parent-v1.5 cluster composition QC
#   - marker-panel summaries and DotPlot for the working res0.6
#   - save a NEW RDS
#
# THIS SCRIPT DOES NOT:
#   - modify/overwrite the v1.6 clean object
#   - assign final human macrophage/monocyte biological labels
#   - freeze a clustering resolution
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#
# IMPORTANT:
#   res0.6 is a DISPLAY/REVIEW resolution only in v1.7.
#   The final resolution will be selected after reviewing cluster structure,
#   sample mixing, parent-cluster composition, and marker profiles.
#
# Seurat v5 convention:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.7"
SEED <- 20260907L

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_MyeloidClean_v1.6/",
  "objects/Human_MASH_LT1to8_MyeloidClean_v1.6.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Myeloid_RPCA_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_Myeloid_RPCA_", VERSION, ".rds")
)

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input v1.6 clean-myeloid RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

set.seed(SEED)

message("============================================================")
message("Human MASH LT1-LT8 clean-myeloid RPCA recluster ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("Output: ", OUTDIR)

src <- readRDS(INPUT_RDS)

required_meta <- c(
  "sample_v15",
  "NAS_v15",
  "parent_immune_cluster_res0_6_v16",
  "entry_class_v16"
)

miss <- setdiff(required_meta, colnames(src[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(src)) stop("RNA assay absent.")
if (!"counts" %in% Layers(src[["RNA"]])) stop("RNA counts layer absent.")

if (ncol(src) != 6956L) {
  stop("Unexpected v1.6 cell number: ", ncol(src), " (expected 6956)")
}
if (nrow(src[["RNA"]]) != 16424L) {
  stop("Unexpected v1.6 RNA feature number: ", nrow(src[["RNA"]]), " (expected 16424)")
}

# ------------------------------------------------------------------
# Fresh analysis object from raw counts only.
# ------------------------------------------------------------------

message("Creating fresh analysis object from v1.6 raw counts...")

raw_counts <- GetAssayData(src, assay = "RNA", layer = "counts")
meta <- src[[]]

obj <- CreateSeuratObject(
  counts = raw_counts,
  assay = "RNA",
  meta.data = meta,
  project = "Human_MASH_LT1to8_Myeloid_RPCA_v1.7",
  min.cells = 0,
  min.features = 0
)

obj$sample_v17 <- as.character(obj$sample_v15)
obj$NAS_v17 <- as.numeric(obj$NAS_v15)
obj$parent_cluster_v17 <- as.character(obj$parent_immune_cluster_res0_6_v16)
obj$entry_class_v17 <- as.character(obj$entry_class_v16)
obj$analysis_version_v17 <- VERSION

rm(raw_counts, meta, src)
gc(verbose = FALSE)

# ------------------------------------------------------------------
# Input inventory
# ------------------------------------------------------------------

sample_counts <- as.data.frame(
  table(sample = obj$sample_v17, NAS = obj$NAS_v17),
  stringsAsFactors = FALSE
)
sample_counts <- sample_counts[sample_counts$Freq > 0, , drop = FALSE]
colnames(sample_counts)[3] <- "n_cells"

write.csv(
  sample_counts,
  file.path(TABDIR, "01_input_sample_counts_v1.7.csv"),
  row.names = FALSE
)

parent_counts <- as.data.frame(
  table(parent_cluster = obj$parent_cluster_v17),
  stringsAsFactors = FALSE
)
colnames(parent_counts)[2] <- "n_cells"

write.csv(
  parent_counts,
  file.path(TABDIR, "02_input_parent_cluster_counts_v1.7.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Seurat v5 sample-aware RPCA integration.
# ------------------------------------------------------------------

message("Splitting RNA layers by biological sample...")
obj[["RNA"]] <- split(
  obj[["RNA"]],
  f = obj$sample_v17
)

message("NormalizeData...")
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

message("FindVariableFeatures...")
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

message("ScaleData...")
obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = FALSE
)

message("RunPCA...")
obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 50,
  seed.use = SEED,
  verbose = FALSE
)

message("IntegrateLayers: RPCAIntegration...")
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.anchor = 20,
  reference = NULL,
  verbose = FALSE
)

if (!"integrated.rpca" %in% Reductions(obj)) {
  stop("RPCA integration did not create integrated.rpca.")
}

message("FindNeighbors...")
obj <- FindNeighbors(
  obj,
  reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

message("RunUMAP...")
obj <- RunUMAP(
  obj,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.myeloid.rpca",
  reduction.key = "MyeloidUMAP_",
  seed.use = SEED,
  verbose = FALSE
)

# ------------------------------------------------------------------
# Resolution scan.
# ------------------------------------------------------------------

RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

for (res in RESOLUTIONS) {
  message("FindClusters resolution ", res, "...")

  obj <- FindClusters(
    obj,
    resolution = res,
    algorithm = 1,
    random.seed = SEED,
    verbose = FALSE
  )

  nm <- paste0(
    "myeloid_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )

  obj[[nm]] <- as.character(Idents(obj))
}

DISPLAY_COL <- "myeloid_RPCA_res0_6"
Idents(obj) <- obj[[DISPLAY_COL, drop = TRUE]]

# Join RNA layers for downstream marker visualization.
message("Joining RNA layers...")
obj <- JoinLayers(obj, assay = "RNA")

# ------------------------------------------------------------------
# QC tables across resolutions.
# ------------------------------------------------------------------

inventory_list <- list()
sample_list <- list()
parent_list <- list()
balance_list <- list()

for (res in RESOLUTIONS) {
  col <- paste0(
    "myeloid_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )

  cl <- as.character(obj[[col, drop = TRUE]])

  inv <- as.data.frame(table(cluster = cl), stringsAsFactors = FALSE)
  colnames(inv)[2] <- "n_cells"
  inv$resolution <- res
  inventory_list[[length(inventory_list) + 1L]] <- inv

  stb <- as.data.frame(
    table(cluster = cl, sample = obj$sample_v17),
    stringsAsFactors = FALSE
  )
  stb <- stb[stb$Freq > 0, , drop = FALSE]
  colnames(stb)[3] <- "n_cells"
  stb$resolution <- res
  sample_list[[length(sample_list) + 1L]] <- stb

  ptb <- as.data.frame(
    table(cluster = cl, parent_cluster = obj$parent_cluster_v17),
    stringsAsFactors = FALSE
  )
  ptb <- ptb[ptb$Freq > 0, , drop = FALSE]
  colnames(ptb)[3] <- "n_cells"
  ptb$resolution <- res
  parent_list[[length(parent_list) + 1L]] <- ptb

  for (cluster_id in sort(unique(cl))) {
    idx <- cl == cluster_id
    tt <- table(obj$sample_v17[idx])
    p <- as.numeric(tt) / sum(tt)
    h <- -sum(p * log(p))
    h_norm <- if (length(tt) > 1) h / log(8) else 0

    balance_list[[length(balance_list) + 1L]] <- data.frame(
      resolution = res,
      cluster = cluster_id,
      n_cells = sum(idx),
      n_samples_present = length(tt),
      dominant_sample = names(tt)[which.max(tt)],
      dominant_sample_fraction = max(p),
      normalized_sample_entropy = h_norm,
      stringsAsFactors = FALSE
    )
  }
}

inventory_all <- do.call(rbind, inventory_list)
sample_all <- do.call(rbind, sample_list)
parent_all <- do.call(rbind, parent_list)
balance_all <- do.call(rbind, balance_list)

write.csv(
  inventory_all,
  file.path(TABDIR, "03_cluster_inventory_all_resolutions_v1.7.csv"),
  row.names = FALSE
)

write.csv(
  sample_all,
  file.path(TABDIR, "04_cluster_by_sample_all_resolutions_v1.7.csv"),
  row.names = FALSE
)

write.csv(
  parent_all,
  file.path(TABDIR, "05_cluster_by_parent_all_resolutions_v1.7.csv"),
  row.names = FALSE
)

write.csv(
  balance_all,
  file.path(TABDIR, "06_sample_balance_all_resolutions_v1.7.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Marker panel for res0.6 review.
# This panel is descriptive only. No cluster labels are assigned.
# ------------------------------------------------------------------

marker_sets <- list(
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Resident_macrophage = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2",
    "TIMD4", "MARCO", "VSIG4", "CD5L"
  ),
  Monocyte_inflammatory = c(
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "FCGR3A",
    "IL1B", "VCAN"
  ),
  SPP1_TREM2_LAM = c(
    "SPP1", "TREM2", "GPNMB", "LGALS3", "LPL", "FABP5", "APOE"
  ),
  ECM_inflammatory = c(
    "F13A1", "TGFB1", "VCAN", "THBS1", "FN1", "POSTN"
  ),
  Repair_resolution = c(
    "MERTK", "ANXA1", "LGALS3", "MRC1", "IL10RA", "IL10RB"
  ),
  Lipid_handling = c(
    "CD36", "MSR1", "OLR1", "APOE", "LPL", "FABP5"
  ),
  Interferon_response = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44", "CXCL10"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  T_NK_contamination = c(
    "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"
  ),
  B_Plasma_contamination = c(
    "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"
  ),
  Mesenchymal_contamination = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  ),
  Hepatocyte_strict_contamination = c(
    "ASGR1", "CPS1", "CYP3A4", "TDO2"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "CENPF", "CCNB1"
  )
)

features <- rownames(obj[["RNA"]])

marker_manifest <- do.call(
  rbind,
  lapply(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    data.frame(
      marker_set = set_name,
      gene = genes,
      present = genes %in% features,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_manifest,
  file.path(TABDIR, "07_marker_presence_v1.7.csv"),
  row.names = FALSE
)

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% features]

rna_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")

cl06 <- as.character(obj[[DISPLAY_COL, drop = TRUE]])
clusters06 <- sort(unique(cl06))

marker_summary_list <- list()

for (cluster_id in clusters06) {
  idx <- cl06 == cluster_id

  pct <- Matrix::rowMeans(
    rna_counts[marker_order, idx, drop = FALSE] > 0
  ) * 100

  avg <- Matrix::rowMeans(
    rna_data[marker_order, idx, drop = FALSE]
  )

  marker_summary_list[[cluster_id]] <- data.frame(
    cluster_res0_6 = cluster_id,
    gene = marker_order,
    pct_detected = as.numeric(pct[marker_order]),
    mean_data = as.numeric(avg[marker_order]),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, marker_summary_list)
rownames(marker_summary) <- NULL

write.csv(
  marker_summary,
  file.path(TABDIR, "08_res0.6_marker_summary_v1.7.csv"),
  row.names = FALSE
)

# Set-level evidence at res0.6.
set_summary_list <- list()

for (cluster_id in clusters06) {
  idx <- cl06 == cluster_id

  for (set_name in names(marker_sets)) {
    genes <- intersect(marker_sets[[set_name]], features)
    if (!length(genes)) next

    sub_counts <- rna_counts[genes, idx, drop = FALSE]
    per_cell_detected <- Matrix::colSums(sub_counts > 0)

    set_summary_list[[length(set_summary_list) + 1L]] <- data.frame(
      cluster_res0_6 = cluster_id,
      marker_set = set_name,
      n_cells = sum(idx),
      n_markers_present = length(genes),
      frac_cells_ge1 = mean(per_cell_detected >= 1),
      frac_cells_ge2 = mean(per_cell_detected >= 2),
      frac_cells_ge3 = mean(per_cell_detected >= 3),
      stringsAsFactors = FALSE
    )
  }
}

set_summary <- do.call(rbind, set_summary_list)

write.csv(
  set_summary,
  file.path(TABDIR, "09_res0.6_marker_set_evidence_v1.7.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

um <- Embeddings(obj, reduction = "umap.myeloid.rpca")

base_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  sample = obj$sample_v17,
  NAS = factor(obj$NAS_v17, levels = c(0, 2, 4, 6)),
  parent_cluster = obj$parent_cluster_v17,
  stringsAsFactors = FALSE,
  row.names = colnames(obj)
)

# Cluster UMAPs for all tested resolutions.
for (res in RESOLUTIONS) {
  col <- paste0(
    "myeloid_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )

  df <- base_df
  df$cluster <- factor(obj[[col, drop = TRUE]])

  centers <- aggregate(
    cbind(UMAP_1, UMAP_2) ~ cluster,
    data = df,
    FUN = median
  )

  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.45, alpha = 0.88) +
    geom_text(
      data = centers,
      aes(label = cluster),
      color = "black",
      fontface = "bold",
      size = 4,
      show.legend = FALSE
    ) +
    coord_equal() +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    ) +
    ggtitle(
      paste0("Human MASH LT1-LT8 clean myeloid RPCA | res ", res),
      "Resolution scan only; no biological state labels assigned"
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "01_myeloid_UMAP_res",
        gsub("\\.", "_", format(res, trim = TRUE)),
        "_v1.7.pdf"
      )
    ),
    p,
    width = 9,
    height = 7.5
  )
}

# Sample UMAP.
p_sample <- ggplot(
  base_df,
  aes(UMAP_1, UMAP_2, color = sample)
) +
  geom_point(size = 0.45, alpha = 0.88) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("Human MASH LT1-LT8 clean myeloid RPCA | biological sample")

ggsave(
  file.path(FIGDIR, "02_myeloid_UMAP_by_sample_v1.7.pdf"),
  p_sample,
  width = 9.5,
  height = 7.5
)

# Parent v1.5 cluster UMAP.
p_parent <- ggplot(
  base_df,
  aes(UMAP_1, UMAP_2, color = factor(parent_cluster))
) +
  geom_point(size = 0.45, alpha = 0.88) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 clean myeloid RPCA | parent immune cluster",
    "Parent cluster identity is QC provenance, not final state annotation"
  )

ggsave(
  file.path(FIGDIR, "03_myeloid_UMAP_parent_cluster_v1.7.pdf"),
  p_parent,
  width = 9.5,
  height = 7.5
)

# res0.6 marker dot plot.
dot_df <- marker_summary
dot_df$cluster_res0_6 <- factor(
  dot_df$cluster_res0_6,
  levels = clusters06
)
dot_df$gene <- factor(
  dot_df$gene,
  levels = rev(marker_order)
)

upper <- as.numeric(
  quantile(dot_df$mean_data, 0.98, na.rm = TRUE)
)
if (!is.finite(upper) || upper <= 0) upper <- 1

p_dot <- ggplot(
  dot_df,
  aes(
    x = cluster_res0_6,
    y = gene,
    size = pct_detected,
    color = mean_data
  )
) +
  geom_point() +
  scale_size(
    range = c(0.2, 6),
    limits = c(0, 100),
    name = "% detected"
  ) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(0, upper),
    oob = scales::squish,
    name = "Mean RNA data"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7.3)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 clean myeloid RPCA | res0.6 marker review",
    "No biological macrophage/monocyte state labels assigned"
  )

ggsave(
  file.path(FIGDIR, "04_myeloid_res0.6_marker_dotplot_v1.7.pdf"),
  p_dot,
  width = max(11, 0.65 * length(clusters06) + 5),
  height = max(14, 0.21 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Save NEW analysis RDS.
# ------------------------------------------------------------------

message("Saving NEW v1.7 myeloid RPCA object...")
saveRDS(obj, OUT_RDS, compress = FALSE)

if (!file.exists(OUT_RDS)) stop("Output RDS was not created.")

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

res_summary <- aggregate(
  cluster ~ resolution,
  data = inventory_all,
  FUN = length
)
colnames(res_summary)[2] <- "n_clusters"

res06_counts <- inventory_all[
  inventory_all$resolution == 0.6,
  c("cluster", "n_cells"),
  drop = FALSE
]

res06_balance <- balance_all[
  balance_all$resolution == 0.6,
  ,
  drop = FALSE
]

summary_lines <- c(
  paste0("Human MASH LT1-LT8 clean-myeloid RPCA recluster ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Input clean-myeloid RDS: ", INPUT_RDS),
  paste0("Output RPCA RDS: ", OUT_RDS),
  "",
  paste0("Cells: ", ncol(obj)),
  paste0("RNA features: ", nrow(obj[["RNA"]])),
  paste0("Samples: ", paste(sort(unique(obj$sample_v17)), collapse = " | ")),
  paste0("Reductions: ", paste(Reductions(obj), collapse = " | ")),
  "",
  "RESOLUTION SCAN:",
  paste(
    capture.output(print(res_summary, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "RES0.6 CLUSTER COUNTS:",
  paste(
    capture.output(print(res06_counts, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "RES0.6 SAMPLE BALANCE:",
  paste(
    capture.output(print(res06_balance, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "IMPORTANT:",
  "- res0.6 is a REVIEW resolution only.",
  "- No biological macrophage/monocyte state labels were assigned.",
  "- No downstream differential analysis was performed.",
  "- Final resolution/state annotation requires review of UMAP structure, sample balance,",
  "  parent-cluster provenance, and marker profiles.",
  "",
  "NEXT DECISION GATE:",
  "Choose the working myeloid resolution and assign preliminary Human MΦ/Monocyte states.",
  "",
  "SOURCE v1.6 RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "10_run_summary_v1.7.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.7 COMPLETE")
message("Cells: ", ncol(obj))
message("RNA features: ", nrow(obj[["RNA"]]))
message("Output RDS: ", OUT_RDS)
message("No biological state labels assigned.")
message("============================================================")
