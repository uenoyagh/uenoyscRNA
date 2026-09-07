#!/usr/bin/env Rscript

# Human_MASH_LT1to8_immune_candidate_RPCA_recluster_v1.5.R
#
# PURPOSE
# -------
# Re-open the ambiguous immune compartment identified by v1.4/v1.4.1 and
# perform an LT1-LT8 sample-aware RPCA reclustering for annotation review.
#
# INPUT CELLS
# -----------
# Transferred broad labels:
#   - Inflammatory_monocyte_macrophage
#   - SPP1_TREM2_macrophage
#   - Plasma_cell
#
# WHY THESE CELLS
# ---------------
# v1.4.1 showed that "Inflammatory_monocyte_macrophage" contains substantial
# non-myeloid T/NK evidence, while the small SPP1_TREM2 group has mixed
# hepatocyte/mesenchymal/myeloid evidence. Plasma cells are retained as a
# positive immune reference population.
#
# THIS SCRIPT DOES:
#   - create a NEW immune-candidate Seurat object from canonical raw RNA counts
#   - normalize each LT sample
#   - identify variable features
#   - PCA
#   - RPCA integration across LT1-LT8
#   - neighbors + UMAP
#   - clustering at resolutions 0.2, 0.4, 0.6, 0.8, 1.0
#   - sample/source-label QC
#   - immune-lineage marker dot plots for review
#   - save a NEW RDS
#
# THIS SCRIPT DOES NOT:
#   - modify the canonical v1.3 RDS
#   - assign biological cluster names
#   - freeze immune annotations
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#   - analyze Cycling_cell (handled separately after immune review)
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.5"
SEED <- 20260907L

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/Human_MASH_LT1to8_Canonical_v1.3.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_ImmuneCandidate_RPCA_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

OUT_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_ImmuneCandidate_RPCA_", VERSION, ".rds")
)

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) stop("Input canonical RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

set.seed(SEED)

message("============================================================")
message("Human MASH LT1-LT8 immune-candidate RPCA recluster ", VERSION)
message("============================================================")
message("Input: ", INPUT_RDS)
message("Output: ", OUTDIR)

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  "celltype_ref12_v13",
  "canonical_sample_v13",
  "NAS_score_v13"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))
if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")

SOURCE_LABELS <- c(
  "Inflammatory_monocyte_macrophage",
  "SPP1_TREM2_macrophage",
  "Plasma_cell"
)

source_label <- as.character(obj$celltype_ref12_v13)
immune_cells <- colnames(obj)[source_label %in% SOURCE_LABELS]

source_counts <- as.data.frame(
  table(source_label = source_label[source_label %in% SOURCE_LABELS]),
  stringsAsFactors = FALSE
)
colnames(source_counts)[2] <- "n_cells"

write.csv(
  source_counts,
  file.path(TABDIR, "01_source_label_counts_v1.5.csv"),
  row.names = FALSE
)

message("Immune-candidate cells: ", length(immune_cells))

if (length(immune_cells) < 1000L) {
  stop("Unexpectedly small immune-candidate compartment.")
}

# ------------------------------------------------------------------
# Build a fresh object from RAW counts only.
# This intentionally avoids carrying inherited scale.data/reductions into
# the new immune analysis.
# ------------------------------------------------------------------

message("Creating fresh immune object from canonical RNA counts...")

canonical_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
immune_counts <- canonical_counts[, immune_cells, drop = FALSE]
immune_meta <- obj[[]][immune_cells, , drop = FALSE]

immune <- CreateSeuratObject(
  counts = immune_counts,
  assay = "RNA",
  meta.data = immune_meta,
  project = "Human_MASH_LT1to8_Immune_v1.5",
  min.cells = 0,
  min.features = 0
)

immune$source_broad_label_v15 <- as.character(
  immune$celltype_ref12_v13
)
immune$sample_v15 <- as.character(
  immune$canonical_sample_v13
)
immune$NAS_v15 <- as.numeric(
  immune$NAS_score_v13
)

rm(canonical_counts, immune_counts, obj)
gc(verbose = FALSE)

# Exact sample inventory before integration.
sample_counts <- as.data.frame(
  table(sample = immune$sample_v15, NAS = immune$NAS_v15),
  stringsAsFactors = FALSE
)
sample_counts <- sample_counts[sample_counts$Freq > 0, , drop = FALSE]
colnames(sample_counts)[3] <- "n_cells"

write.csv(
  sample_counts,
  file.path(TABDIR, "02_sample_counts_v1.5.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Seurat v5 RPCA integration across LT1-LT8.
# ------------------------------------------------------------------

message("Splitting RNA layers by biological sample...")
immune[["RNA"]] <- split(
  immune[["RNA"]],
  f = immune$sample_v15
)

message("NormalizeData...")
immune <- NormalizeData(
  immune,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

message("FindVariableFeatures...")
immune <- FindVariableFeatures(
  immune,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

message("ScaleData...")
immune <- ScaleData(
  immune,
  features = VariableFeatures(immune),
  verbose = FALSE
)

message("RunPCA...")
immune <- RunPCA(
  immune,
  features = VariableFeatures(immune),
  npcs = 50,
  verbose = FALSE,
  seed.use = SEED
)

message("IntegrateLayers: RPCAIntegration...")
immune <- IntegrateLayers(
  object = immune,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.anchor = 20,
  reference = NULL,
  verbose = FALSE
)

if (!"integrated.rpca" %in% Reductions(immune)) {
  stop("RPCA integration did not produce 'integrated.rpca'.")
}

message("FindNeighbors...")
immune <- FindNeighbors(
  immune,
  reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

message("RunUMAP...")
immune <- RunUMAP(
  immune,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.immune.rpca",
  reduction.key = "ImmuneUMAP_",
  seed.use = SEED,
  verbose = FALSE
)

# ------------------------------------------------------------------
# Multi-resolution clustering for INSPECTION ONLY.
# ------------------------------------------------------------------

RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0)

for (res in RESOLUTIONS) {
  message("FindClusters resolution ", res, "...")
  immune <- FindClusters(
    immune,
    resolution = res,
    algorithm = 1,
    random.seed = SEED,
    verbose = FALSE
  )

  nm <- paste0(
    "immune_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )

  immune[[nm]] <- as.character(Idents(immune))
}

# Use res0.6 only as the DISPLAY/inspection identity, not a frozen biological choice.
DISPLAY_CLUSTER_COL <- "immune_RPCA_res0_6"
Idents(immune) <- immune[[DISPLAY_CLUSTER_COL, drop = TRUE]]

# ------------------------------------------------------------------
# Join RNA layers after integration so marker visualization reads one layer.
# ------------------------------------------------------------------

message("Joining RNA layers for marker visualization...")
immune <- JoinLayers(immune, assay = "RNA")

rna_layers <- Layers(immune[["RNA"]])

# ------------------------------------------------------------------
# Cluster/sample/source-label QC tables at every resolution.
# ------------------------------------------------------------------

cluster_inventory <- list()
cluster_sample <- list()
cluster_source <- list()

for (res in RESOLUTIONS) {
  col <- paste0(
    "immune_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )

  cl <- as.character(immune[[col, drop = TRUE]])

  a <- as.data.frame(table(cluster = cl), stringsAsFactors = FALSE)
  colnames(a)[2] <- "n_cells"
  a$resolution <- res
  cluster_inventory[[length(cluster_inventory) + 1L]] <- a

  b <- as.data.frame(
    table(cluster = cl, sample = immune$sample_v15),
    stringsAsFactors = FALSE
  )
  b <- b[b$Freq > 0, , drop = FALSE]
  colnames(b)[3] <- "n_cells"
  b$resolution <- res
  cluster_sample[[length(cluster_sample) + 1L]] <- b

  ctab <- as.data.frame(
    table(cluster = cl, source_label = immune$source_broad_label_v15),
    stringsAsFactors = FALSE
  )
  ctab <- ctab[ctab$Freq > 0, , drop = FALSE]
  colnames(ctab)[3] <- "n_cells"
  ctab$resolution <- res
  cluster_source[[length(cluster_source) + 1L]] <- ctab
}

write.csv(
  do.call(rbind, cluster_inventory),
  file.path(TABDIR, "03_cluster_inventory_all_resolutions_v1.5.csv"),
  row.names = FALSE
)

write.csv(
  do.call(rbind, cluster_sample),
  file.path(TABDIR, "04_cluster_by_sample_all_resolutions_v1.5.csv"),
  row.names = FALSE
)

write.csv(
  do.call(rbind, cluster_source),
  file.path(TABDIR, "05_cluster_by_source_label_all_resolutions_v1.5.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Marker panel for immune-lineage identification.
# NO cluster annotation is assigned here.
# ------------------------------------------------------------------

marker_sets <- list(
  Myeloid_core = c(
    "LYZ", "LST1", "TYROBP", "FCER1G", "CSF1R", "CD68", "CTSS"
  ),
  Monocyte_inflammatory = c(
    "S100A8", "S100A9", "CCR2", "LILRB1", "MNDA", "FCGR3A"
  ),
  Resident_macrophage = c(
    "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "FOLR2",
    "TIMD4", "MARCO", "VSIG4", "CD5L"
  ),
  SPP1_TREM2 = c(
    "SPP1", "TREM2", "GPNMB", "LGALS3", "LPL", "FABP5", "APOE"
  ),
  T_cell = c(
    "CD3D", "CD3E", "TRAC", "CD247", "IL7R"
  ),
  NK_cell = c(
    "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRK1", "GZMB", "PRF1"
  ),
  B_cell = c(
    "MS4A1", "CD79A", "CD79B", "CD37", "CD74", "CD19"
  ),
  Plasma_cell = c(
    "MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHM"
  ),
  Dendritic = c(
    "FCER1A", "CD1C", "CLEC10A", "CLEC9A", "XCR1", "CD74"
  ),
  Hepatocyte_contamination = c(
    "ALB", "TTR", "APOA1", "ASGR1", "CPS1"
  ),
  Mesenchymal_contamination = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"
  )
)

marker_manifest <- do.call(
  rbind,
  lapply(names(marker_sets), function(set_name) {
    genes <- marker_sets[[set_name]]
    data.frame(
      marker_set = set_name,
      gene = genes,
      present = genes %in% rownames(immune[["RNA"]]),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  marker_manifest,
  file.path(TABDIR, "06_marker_presence_v1.5.csv"),
  row.names = FALSE
)

marker_order <- unique(unlist(marker_sets, use.names = FALSE))
marker_order <- marker_order[marker_order %in% rownames(immune[["RNA"]])]

rna_counts <- GetAssayData(immune, assay = "RNA", layer = "counts")
rna_data <- GetAssayData(immune, assay = "RNA", layer = "data")

cl06 <- as.character(immune[[DISPLAY_CLUSTER_COL, drop = TRUE]])
clusters06 <- sort(unique(cl06))

marker_summary <- list()

for (cl in clusters06) {
  idx <- cl06 == cl
  pct <- Matrix::rowMeans(rna_counts[marker_order, idx, drop = FALSE] > 0) * 100
  avg <- Matrix::rowMeans(rna_data[marker_order, idx, drop = FALSE])

  marker_summary[[cl]] <- data.frame(
    cluster_res0_6 = cl,
    gene = marker_order,
    pct_detected = as.numeric(pct[marker_order]),
    mean_data = as.numeric(avg[marker_order]),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}

marker_summary <- do.call(rbind, marker_summary)
rownames(marker_summary) <- NULL

write.csv(
  marker_summary,
  file.path(TABDIR, "07_res0.6_marker_summary_v1.5.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

um <- Embeddings(immune, reduction = "umap.immune.rpca")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  sample = immune$sample_v15,
  NAS = factor(immune$NAS_v15, levels = c(0, 2, 4, 6)),
  source_label = immune$source_broad_label_v15,
  stringsAsFactors = FALSE,
  row.names = colnames(immune)
)

# Cluster UMAPs for every tested resolution.
for (res in RESOLUTIONS) {
  col <- paste0(
    "immune_RPCA_res",
    gsub("\\.", "_", format(res, trim = TRUE))
  )
  plot_df$cluster <- factor(immune[[col, drop = TRUE]])

  p <- ggplot(
    plot_df,
    aes(UMAP_1, UMAP_2, color = cluster)
  ) +
    geom_point(size = 0.35, alpha = 0.85) +
    coord_equal() +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    ) +
    ggtitle(
      paste0("Human MASH LT1-LT8 immune candidates | RPCA res ", res),
      "Inspection only; no biological cluster labels assigned"
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "01_immune_UMAP_res",
        gsub("\\.", "_", format(res, trim = TRUE)),
        "_v1.5.pdf"
      )
    ),
    p,
    width = 9,
    height = 7.5
  )
}

# Source-label UMAP.
p_source <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = source_label)
) +
  geom_point(size = 0.35, alpha = 0.85) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 immune RPCA | transferred source labels",
    "Source labels are not treated as final biology"
  )

ggsave(
  file.path(FIGDIR, "02_immune_UMAP_source_labels_v1.5.pdf"),
  p_source,
  width = 9.5,
  height = 7.5
)

# Biological-sample UMAP.
p_sample <- ggplot(
  plot_df,
  aes(UMAP_1, UMAP_2, color = sample)
) +
  geom_point(size = 0.35, alpha = 0.85) +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("Human MASH LT1-LT8 immune RPCA | biological sample")

ggsave(
  file.path(FIGDIR, "03_immune_UMAP_by_sample_v1.5.pdf"),
  p_sample,
  width = 9.5,
  height = 7.5
)

# res0.6 marker dot plot.
dot_df <- marker_summary
dot_df$gene <- factor(dot_df$gene, levels = rev(marker_order))
dot_df$cluster_res0_6 <- factor(
  dot_df$cluster_res0_6,
  levels = clusters06
)

upper <- as.numeric(quantile(dot_df$mean_data, 0.98, na.rm = TRUE))
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
    axis.text.x = element_text(angle = 0, size = 9),
    axis.text.y = element_text(size = 7.5)
  ) +
  ggtitle(
    "Human MASH LT1-LT8 immune candidate RPCA | res0.6 marker review",
    "No biological cluster names assigned"
  )

ggsave(
  file.path(FIGDIR, "04_immune_res0.6_marker_dotplot_v1.5.pdf"),
  p_dot,
  width = max(10, 0.65 * length(clusters06) + 5),
  height = max(12, 0.22 * length(marker_order) + 3)
)

# ------------------------------------------------------------------
# Save NEW immune object.
# ------------------------------------------------------------------

message("Saving NEW immune RPCA object...")
saveRDS(immune, OUT_RDS, compress = FALSE)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

cl06_counts <- as.data.frame(table(cluster = cl06), stringsAsFactors = FALSE)
colnames(cl06_counts)[2] <- "n_cells"

summary_lines <- c(
  paste0("Human MASH LT1-LT8 immune-candidate RPCA recluster ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Input canonical object: ", INPUT_RDS),
  paste0("Output immune RDS: ", OUT_RDS),
  "",
  paste0("Immune-candidate cells: ", ncol(immune)),
  paste0("RNA features: ", nrow(immune[["RNA"]])),
  paste0("Samples: ", paste(sort(unique(immune$sample_v15)), collapse = " | ")),
  paste0("RNA layers after JoinLayers: ", paste(rna_layers, collapse = " | ")),
  paste0("Reductions: ", paste(Reductions(immune), collapse = " | ")),
  "",
  "SOURCE LABEL COUNTS:",
  paste(capture.output(print(source_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "RES0.6 CLUSTER COUNTS:",
  paste(capture.output(print(cl06_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "RESOLUTIONS GENERATED:",
  paste(RESOLUTIONS, collapse = " | "),
  "",
  "IMPORTANT:",
  "- No biological cluster names were assigned.",
  "- res0.6 is only the marker-review display resolution.",
  "- Review resolution structure, sample mixing, source-label distribution, and marker DotPlot",
  "  before selecting a resolution or assigning immune identities.",
  "- Cycling_cell is NOT included in this v1.5 immune object.",
  "",
  "SOURCE CANONICAL RDS MODIFICATION:",
  "NONE"
)

writeLines(
  summary_lines,
  file.path(TABDIR, "08_run_summary_v1.5.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v1.5 COMPLETE")
message("Immune-candidate cells: ", ncol(immune))
message("Output RDS: ", OUT_RDS)
message("No biological cluster names assigned.")
message("============================================================")
