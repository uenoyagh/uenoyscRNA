suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.7.1.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_v6.7.0/",
  "objects/",
  "Mouse_MASH_endothelial_parent_v6.7.0.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_",
  VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

cat("====================================================\n")
cat("Mouse MASH endothelial clean RPCA reclustering\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

cat("Reading parent object:\n")
cat(INPUT_RDS, "\n\n")

parent <- readRDS(INPUT_RDS)

cat("=== ORIGINAL PARENT ===\n")
cat("Cells:", ncol(parent), "\n")
cat("Features:", nrow(parent), "\n")
cat("Assays:", paste(Assays(parent), collapse = ", "), "\n")
cat(
  "Reductions:",
  ifelse(
    length(Reductions(parent)) == 0,
    "none",
    paste(Reductions(parent), collapse = ", ")
  ),
  "\n"
)
cat("Graphs:", length(parent@graphs), "\n\n")

required_cols <- c(
  "sample",
  "condition",
  "endothelial_parent_label_v670"
)

missing_cols <- setdiff(
  required_cols,
  colnames(parent@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_cols, collapse = ", ")
  )
}

DefaultAssay(parent) <- "RNA"

cat("=== EXTRACT RAW RNA COUNTS ===\n")

if (inherits(parent[["RNA"]], "Assay5")) {

  cat("RNA assay is Assay5.\n")
  cat(
    "RNA layers:",
    paste(Layers(parent[["RNA"]]), collapse = ", "),
    "\n"
  )

  count_layers <- grep(
    "^counts($|\\.)",
    Layers(parent[["RNA"]]),
    value = TRUE
  )

  if (length(count_layers) == 0) {
    stop("No RNA counts layer found.")
  }

  if (
    length(count_layers) == 1 &&
    identical(count_layers, "counts")
  ) {

    counts <- LayerData(
      parent,
      assay = "RNA",
      layer = "counts"
    )

  } else {

    cat(
      "Joining RNA layers temporarily for raw-count extraction...\n"
    )

    parent <- JoinLayers(
      parent,
      assay = "RNA"
    )

    counts <- LayerData(
      parent,
      assay = "RNA",
      layer = "counts"
    )
  }

} else {

  cat("RNA assay is legacy Assay.\n")

  counts <- GetAssayData(
    parent,
    assay = "RNA",
    slot = "counts"
  )
}

cat("Count matrix genes:", nrow(counts), "\n")
cat("Count matrix cells:", ncol(counts), "\n")

if (ncol(counts) != ncol(parent)) {
  stop("Count matrix cell number differs from parent object.")
}

if (!all(colnames(counts) %in% rownames(parent@meta.data))) {
  stop("Count matrix cell names do not match metadata.")
}

metadata <- parent@meta.data[
  colnames(counts),
  ,
  drop = FALSE
]

cat("\n=== CREATE COMPLETELY CLEAN OBJECT ===\n")

options(
  Seurat.object.assay.version = "v3"
)

clean <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = "Mouse_MASH_LSEC"
)

rm(parent, counts)
invisible(gc())

cat("Clean cells:", ncol(clean), "\n")
cat("Clean features:", nrow(clean), "\n")
cat(
  "Clean assays:",
  paste(Assays(clean), collapse = ", "),
  "\n"
)
cat(
  "RNA assay class:",
  class(clean[["RNA"]])[1],
  "\n"
)
cat(
  "Clean reductions:",
  length(Reductions(clean)),
  "\n"
)
cat(
  "Clean graphs:",
  length(clean@graphs),
  "\n"
)

if (!identical(Assays(clean), "RNA")) {
  stop("Clean object contains unexpected assay(s).")
}

if (length(Reductions(clean)) != 0) {
  stop("Clean object unexpectedly contains reductions.")
}

if (length(clean@graphs) != 0) {
  stop("Clean object unexpectedly contains graphs.")
}

if (!inherits(clean[["RNA"]], "Assay")) {
  stop(
    "RNA assay is not legacy Assay; clean rebuild failed."
  )
}

cat("\nCLEAN OBJECT CHECK: PASSED\n")

cat("\n=== CELLS BY SAMPLE ===\n")
print(table(clean$sample))

cat("\n=== CELLS BY CONDITION ===\n")
print(table(clean$condition))

cat("\n=== ORIGINAL ENDOTHELIAL LABEL ===\n")
print(table(clean$endothelial_parent_label_v670))

write.csv(
  as.data.frame(table(clean$sample)),
  file.path(
    TABDIR,
    "Cells_by_sample_v6.7.1.1.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(clean$condition)),
  file.path(
    TABDIR,
    "Cells_by_condition_v6.7.1.1.csv"
  ),
  row.names = FALSE
)

cat("\n=== SPLIT BY BIOLOGICAL SAMPLE ===\n")

sample_list <- SplitObject(
  clean,
  split.by = "sample"
)

for (nm in names(sample_list)) {
  cat(
    nm,
    ":",
    ncol(sample_list[[nm]]),
    "cells\n"
  )
}

cat("\n=== NORMALIZE + HVG ===\n")

sample_list <- lapply(
  sample_list,
  function(x) {

    DefaultAssay(x) <- "RNA"

    x <- NormalizeData(
      x,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )

    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 3000,
      verbose = FALSE
    )

    return(x)
  }
)

integration_features <- SelectIntegrationFeatures(
  object.list = sample_list,
  nfeatures = 3000
)

cat(
  "Integration features:",
  length(integration_features),
  "\n"
)

cat("\n=== SCALE + SAMPLE PCA ===\n")

sample_list <- lapply(
  sample_list,
  function(x) {

    x <- ScaleData(
      x,
      features = integration_features,
      verbose = FALSE
    )

    x <- RunPCA(
      x,
      features = integration_features,
      npcs = 40,
      verbose = FALSE
    )

    return(x)
  }
)

cat("\n=== FIND RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list = sample_list,
  anchor.features = integration_features,
  reduction = "rpca",
  dims = 1:30
)

cat("\n=== INTEGRATE DATA ===\n")

endo <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  new.assay.name = "integratedLSEC"
)

cat("\nIntegrated cells:", ncol(endo), "\n")
cat(
  "Assays:",
  paste(Assays(endo), collapse = ", "),
  "\n"
)

if (ncol(endo) != ncol(clean)) {
  stop(
    "Cell number changed during integration."
  )
}

DefaultAssay(endo) <- "integratedLSEC"

cat("\n=== ENDOTHELIAL PCA ===\n")

endo <- ScaleData(
  endo,
  verbose = FALSE
)

endo <- RunPCA(
  endo,
  npcs = 40,
  verbose = FALSE,
  reduction.name = "pcaLSEC",
  reduction.key = "LSECPCA_"
)

pca_sd <- Stdev(
  endo,
  reduction = "pcaLSEC"
)

pca_df <- data.frame(
  PC = seq_along(pca_sd),
  SD = pca_sd
)

pca_df$variance <- pca_df$SD^2

pca_df$variance_fraction <-
  pca_df$variance /
  sum(pca_df$variance)

pca_df$cumulative_variance <-
  cumsum(
    pca_df$variance_fraction
  )

write.csv(
  pca_df,
  file.path(
    TABDIR,
    "PCA_variance_v6.7.1.1.csv"
  ),
  row.names = FALSE
)

p_elbow <- ggplot(
  pca_df[
    1:min(40, nrow(pca_df)),
    ,
    drop = FALSE
  ],
  aes(
    x = PC,
    y = SD
  )
) +
  geom_point(size = 2) +
  geom_line() +
  theme_classic(base_size = 13) +
  labs(
    title = "Endothelial PCA",
    x = "PC",
    y = "Standard deviation"
  )

ggsave(
  file.path(
    FIGDIR,
    "PCA_elbow_v6.7.1.1.pdf"
  ),
  p_elbow,
  width = 7,
  height = 5
)

DIMS_USE <- 1:30

cat("\n=== NEIGHBORS + UMAP ===\n")

endo <- FindNeighbors(
  endo,
  reduction = "pcaLSEC",
  dims = DIMS_USE,
  graph.name = c(
    "LSEC_nn",
    "LSEC_snn"
  ),
  verbose = FALSE
)

endo <- RunUMAP(
  endo,
  reduction = "pcaLSEC",
  dims = DIMS_USE,
  reduction.name = "umapLSEC",
  reduction.key = "LSECUMAP_",
  seed.use = 20260902,
  verbose = FALSE
)

resolutions <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0,
  1.2
)

cat("\n=== RESOLUTION SCAN ===\n")

for (res in resolutions) {

  cat(
    "\n------------------------------\n"
  )
  cat("Resolution:", res, "\n")

  endo <- FindClusters(
    endo,
    graph.name = "LSEC_snn",
    resolution = res,
    algorithm = 1,
    random.seed = 20260902,
    verbose = FALSE
  )

  cluster_col <- paste0(
    "LSEC_res",
    format(
      res,
      nsmall = 1
    )
  )

  endo[[cluster_col]] <-
    as.character(Idents(endo))

  x <- sort(
    table(
      endo@meta.data[[cluster_col]]
    ),
    decreasing = TRUE
  )

  print(x)

  cat(
    "Number of clusters:",
    length(x),
    "\n"
  )
}

um <- Embeddings(
  endo,
  reduction = "umapLSEC"
)

base_df <- data.frame(
  cell = rownames(um),
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  sample = as.character(endo$sample),
  condition = as.character(endo$condition),
  original_class =
    as.character(
      endo$endothelial_parent_label_v670
    ),
  stringsAsFactors = FALSE
)

condition_colors <- c(
  "STD" = "#0057FF",
  "CDHFD" = "#FF1A1A",
  "Sham" = "#00A651",
  "Tx" = "#AA00FF",
  "Other" = "grey65"
)

p_condition <- ggplot(
  base_df,
  aes(
    UMAP_1,
    UMAP_2,
    color = condition
  )
) +
  geom_point(
    size = 0.45,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = condition_colors
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Endothelial RPCA UMAP by condition"
  )

ggsave(
  file.path(
    FIGDIR,
    "Endothelial_RPCA_UMAP_by_condition_v6.7.1.1.pdf"
  ),
  p_condition,
  width = 8,
  height = 7
)

p_sample <- ggplot(
  base_df,
  aes(
    UMAP_1,
    UMAP_2,
    color = sample
  )
) +
  geom_point(
    size = 0.45,
    alpha = 0.85
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Endothelial RPCA UMAP by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Endothelial_RPCA_UMAP_by_sample_v6.7.1.1.pdf"
  ),
  p_sample,
  width = 9,
  height = 7
)

p_original <- ggplot(
  base_df,
  aes(
    UMAP_1,
    UMAP_2,
    color = original_class
  )
) +
  geom_point(
    size = 0.45,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = c(
      "LSEC" = "#00BFC4",
      "Vascular_endothelial" = "#0066FF"
    )
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Original endothelial annotation"
  )

ggsave(
  file.path(
    FIGDIR,
    "Endothelial_RPCA_UMAP_original_annotation_v6.7.1.1.pdf"
  ),
  p_original,
  width = 8,
  height = 7
)

cat("\n=== WRITE RESOLUTION OUTPUTS ===\n")

for (res in resolutions) {

  cluster_col <- paste0(
    "LSEC_res",
    format(
      res,
      nsmall = 1
    )
  )

  tmp <- base_df

  tmp$cluster <- factor(
    endo@meta.data[[cluster_col]]
  )

  p <- ggplot(
    tmp,
    aes(
      UMAP_1,
      UMAP_2,
      color = cluster
    )
  ) +
    geom_point(
      size = 0.48,
      alpha = 0.90
    ) +
    coord_equal() +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    ) +
    ggtitle(
      paste0(
        "Endothelial RPCA - resolution ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Endothelial_RPCA_UMAP_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.pdf"
      )
    ),
    p,
    width = 8.5,
    height = 7
  )

  write.csv(
    as.data.frame(
      table(
        cluster =
          endo@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_counts_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(
      table(
        sample = endo$sample,
        cluster =
          endo@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_sample_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(
      table(
        condition = endo$condition,
        cluster =
          endo@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_condition_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(
      table(
        original_class =
          endo$endothelial_parent_label_v670,
        cluster =
          endo@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_original_annotation_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.csv"
      )
    ),
    row.names = FALSE
  )
}

cat("\n=== MARKER AUDIT ===\n")

DefaultAssay(endo) <- "RNA"

marker_panel <- c(
  "Pecam1",
  "Cdh5",
  "Kdr",
  "Klf2",
  "Rgcc",

  "Clec4g",
  "Stab1",
  "Stab2",
  "Lyve1",
  "Fcgr2b",
  "Mrc1",

  "Plvap",
  "Vwf",
  "Emcn",
  "Esm1",
  "Cd34",
  "Car4",

  "Wnt2",
  "Wnt9b",
  "Rspo3",
  "Bmp2",
  "Efnb2",
  "Dll4",
  "Ackr1",

  "Ptprc",
  "Lyz2",
  "Adgre1",

  "Col1a1",
  "Col3a1",
  "Pdgfrb",
  "Rgs5",

  "Alb",
  "Krt19"
)

marker_present <- marker_panel[
  marker_panel %in% rownames(endo)
]

marker_missing <- setdiff(
  marker_panel,
  marker_present
)

cat(
  "Present:",
  paste(marker_present, collapse = ", "),
  "\n"
)

cat(
  "Missing:",
  paste(marker_missing, collapse = ", "),
  "\n"
)

writeLines(
  c(
    paste(
      "Present:",
      paste(marker_present, collapse = ", ")
    ),
    paste(
      "Missing:",
      paste(marker_missing, collapse = ", ")
    )
  ),
  file.path(
    TABDIR,
    "Marker_availability_v6.7.1.1.txt"
  )
)

dotplot_resolutions <- c(
  0.4,
  0.6,
  0.8,
  1.0
)

for (res in dotplot_resolutions) {

  cluster_col <- paste0(
    "LSEC_res",
    format(
      res,
      nsmall = 1
    )
  )

  Idents(endo) <- factor(
    endo@meta.data[[cluster_col]]
  )

  p_dot <- DotPlot(
    endo,
    features = marker_present,
    assay = "RNA",
    dot.scale = 6
  ) +
    RotatedAxis() +
    theme_classic(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(
      paste0(
        "Endothelial markers - resolution ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Endothelial_marker_DotPlot_res",
        format(res, nsmall = 1),
        "_v6.7.1.1.pdf"
      )
    ),
    p_dot,
    width = 16,
    height = 6.5
  )
}

cat("\n=== SAVE FINAL RESOLUTION-SCAN OBJECT ===\n")

saveRDS(
  endo,
  file.path(
    OBJDIR,
    "Mouse_MASH_endothelial_clean_RPCA_resolution_scan_v6.7.1.1.rds"
  ),
  compress = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo_v6.7.1.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.1.1 COMPLETE\n")
cat("Cells:", ncol(endo), "\n")
cat(
  "Resolutions:",
  paste(resolutions, collapse = ", "),
  "\n"
)
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
