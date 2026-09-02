suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.7.4"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.3.1/objects/",
  "Mouse_MASH_LSEC_raw_clean_v6.7.3.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH clean LSEC RPCA reclustering\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

obj <- readRDS(INPUT_RDS)

required <- c("sample", "condition")

missing <- setdiff(
  required,
  colnames(obj@meta.data)
)

if (length(missing) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing, collapse=", ")
  )
}

cat("=== INPUT CHECK ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")
cat("RNA class:", class(obj[["RNA"]])[1], "\n")
cat("Reductions:", length(Reductions(obj)), "\n")
cat("Graphs:", length(obj@graphs), "\n")

if (ncol(obj) != 9190) {
  stop(
    "Unexpected input cell number: ",
    ncol(obj)
  )
}

if (!identical(Assays(obj), "RNA")) {
  stop("Unexpected assay(s) in clean LSEC object.")
}

if (length(Reductions(obj)) != 0) {
  stop("Input object unexpectedly contains reductions.")
}

if (length(obj@graphs) != 0) {
  stop("Input object unexpectedly contains graphs.")
}

cat("\nCLEAN INPUT CHECK: PASSED\n")

cat("\n=== CELLS BY SAMPLE ===\n")
print(table(obj$sample))

cat("\n=== CELLS BY CONDITION ===\n")
print(table(obj$condition))

# ---------------------------------------------------------
# Split by biological sample
# ---------------------------------------------------------

DefaultAssay(obj) <- "RNA"

sample_list <- SplitObject(
  obj,
  split.by="sample"
)

cat("\n=== SPLIT OBJECTS ===\n")

for (nm in names(sample_list)) {
  cat(
    nm, ": ",
    ncol(sample_list[[nm]]),
    " cells\n",
    sep=""
  )
}

# ---------------------------------------------------------
# Normalize and variable features
# ---------------------------------------------------------

cat("\n=== NORMALIZE + HVG ===\n")

sample_list <- lapply(
  sample_list,
  function(x) {

    DefaultAssay(x) <- "RNA"

    x <- NormalizeData(
      x,
      normalization.method="LogNormalize",
      scale.factor=10000,
      verbose=FALSE
    )

    x <- FindVariableFeatures(
      x,
      selection.method="vst",
      nfeatures=3000,
      verbose=FALSE
    )

    x
  }
)

integration_features <- SelectIntegrationFeatures(
  object.list=sample_list,
  nfeatures=3000
)

cat(
  "Integration features:",
  length(integration_features),
  "\n"
)

# ---------------------------------------------------------
# Sample PCA for RPCA anchors
# ---------------------------------------------------------

cat("\n=== SAMPLE SCALE + PCA ===\n")

sample_list <- lapply(
  sample_list,
  function(x) {

    x <- ScaleData(
      x,
      features=integration_features,
      verbose=FALSE
    )

    x <- RunPCA(
      x,
      features=integration_features,
      npcs=30,
      verbose=FALSE
    )

    x
  }
)

# ---------------------------------------------------------
# RPCA integration
# ---------------------------------------------------------

ANCHOR_DIMS <- 1:20

cat("\n=== FIND RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list=sample_list,
  anchor.features=integration_features,
  reduction="rpca",
  dims=ANCHOR_DIMS
)

cat("\n=== INTEGRATE LSEC ===\n")

lsec <- IntegrateData(
  anchorset=anchors,
  dims=ANCHOR_DIMS,
  new.assay.name="integratedLSEC"
)

cat("Integrated cells:", ncol(lsec), "\n")
cat(
  "Assays:",
  paste(Assays(lsec), collapse=", "),
  "\n"
)

if (ncol(lsec) != 9190) {
  stop(
    "Cell count changed during integration: ",
    ncol(lsec)
  )
}

# ---------------------------------------------------------
# LSEC-only PCA
# ---------------------------------------------------------

DefaultAssay(lsec) <- "integratedLSEC"

lsec <- ScaleData(
  lsec,
  verbose=FALSE
)

lsec <- RunPCA(
  lsec,
  npcs=30,
  verbose=FALSE,
  reduction.name="pcaLSECclean",
  reduction.key="LSECCLEANPC_"
)

pca_sd <- Stdev(
  lsec,
  reduction="pcaLSECclean"
)

pca_df <- data.frame(
  PC=seq_along(pca_sd),
  SD=pca_sd
)

pca_df$variance <- pca_df$SD^2
pca_df$variance_fraction <-
  pca_df$variance / sum(pca_df$variance)

pca_df$cumulative_variance <-
  cumsum(pca_df$variance_fraction)

write.csv(
  pca_df,
  file.path(
    TABDIR,
    "LSEC_clean_PCA_variance_v6.7.4.csv"
  ),
  row.names=FALSE
)

p_elbow <- ggplot(
  pca_df,
  aes(PC, SD)
) +
  geom_point(size=2) +
  geom_line() +
  theme_classic(base_size=13) +
  labs(
    title="Clean LSEC PCA",
    x="PC",
    y="Standard deviation"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_clean_PCA_elbow_v6.7.4.pdf"
  ),
  p_elbow,
  width=7,
  height=5
)

# ---------------------------------------------------------
# Baseline graph / UMAP
# ---------------------------------------------------------

DIMS_USE <- 1:20

cat("\n=== NEIGHBORS + UMAP ===\n")
cat("Dimensions used: 1:20\n")

lsec <- FindNeighbors(
  lsec,
  reduction="pcaLSECclean",
  dims=DIMS_USE,
  graph.name=c(
    "LSECclean_nn",
    "LSECclean_snn"
  ),
  verbose=FALSE
)

lsec <- RunUMAP(
  lsec,
  reduction="pcaLSECclean",
  dims=DIMS_USE,
  reduction.name="umapLSECclean",
  reduction.key="LSECCLEANUMAP_",
  seed.use=20260902,
  verbose=FALSE
)

# ---------------------------------------------------------
# Resolution scan
# ---------------------------------------------------------

resolutions <- c(
  0.1,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6
)

cat("\n=== RESOLUTION SCAN ===\n")

for (res in resolutions) {

  cat("\n------------------------------\n")
  cat("Resolution:", res, "\n")

  lsec <- FindClusters(
    lsec,
    graph.name="LSECclean_snn",
    resolution=res,
    algorithm=1,
    random.seed=20260902,
    verbose=FALSE
  )

  cluster_col <- paste0(
    "LSECclean_res",
    format(res, nsmall=1)
  )

  lsec[[cluster_col]] <-
    as.character(Idents(lsec))

  x <- sort(
    table(
      lsec@meta.data[[cluster_col]]
    ),
    decreasing=TRUE
  )

  print(x)

  cat(
    "Number of clusters:",
    length(x),
    "\n"
  )

  write.csv(
    as.data.frame(
      table(
        cluster=
          lsec@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_counts_res",
        format(res, nsmall=1),
        "_v6.7.4.csv"
      )
    ),
    row.names=FALSE
  )

  write.csv(
    as.data.frame(
      table(
        sample=lsec$sample,
        cluster=
          lsec@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_sample_res",
        format(res, nsmall=1),
        "_v6.7.4.csv"
      )
    ),
    row.names=FALSE
  )

  write.csv(
    as.data.frame(
      table(
        condition=lsec$condition,
        cluster=
          lsec@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_condition_res",
        format(res, nsmall=1),
        "_v6.7.4.csv"
      )
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# UMAP base dataframe
# ---------------------------------------------------------

um <- Embeddings(
  lsec,
  reduction="umapLSECclean"
)

base_df <- data.frame(
  cell=rownames(um),
  UMAP_1=um[,1],
  UMAP_2=um[,2],
  sample=as.character(lsec$sample),
  condition=as.character(lsec$condition),
  stringsAsFactors=FALSE
)

# ---------------------------------------------------------
# Sample UMAP
# ---------------------------------------------------------

p_sample <- ggplot(
  base_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=sample
  )
) +
  geom_point(
    size=0.50,
    alpha=0.85
  ) +
  coord_equal() +
  theme_classic(base_size=13) +
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.title=element_blank()
  ) +
  ggtitle(
    "Clean LSEC RPCA UMAP by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_clean_RPCA_UMAP_by_sample_v6.7.4.pdf"
  ),
  p_sample,
  width=9,
  height=7
)

# ---------------------------------------------------------
# Condition UMAP
# ---------------------------------------------------------

condition_colors <- c(
  "STD"="#0057FF",
  "CDHFD"="#FF1A1A",
  "Sham"="#00A651",
  "Tx"="#AA00FF"
)

p_condition <- ggplot(
  base_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=condition
  )
) +
  geom_point(
    size=0.50,
    alpha=0.85
  ) +
  scale_color_manual(
    values=condition_colors
  ) +
  coord_equal() +
  theme_classic(base_size=13) +
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.title=element_blank()
  ) +
  ggtitle(
    "Clean LSEC RPCA UMAP by condition"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_clean_RPCA_UMAP_by_condition_v6.7.4.pdf"
  ),
  p_condition,
  width=8,
  height=7
)

# ---------------------------------------------------------
# Cluster UMAPs
# ---------------------------------------------------------

for (res in resolutions) {

  cluster_col <- paste0(
    "LSECclean_res",
    format(res, nsmall=1)
  )

  tmp <- base_df
  tmp$cluster <- factor(
    lsec@meta.data[[cluster_col]]
  )

  p <- ggplot(
    tmp,
    aes(
      UMAP_1,
      UMAP_2,
      color=cluster
    )
  ) +
    geom_point(
      size=0.52,
      alpha=0.90
    ) +
    coord_equal() +
    theme_classic(base_size=13) +
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      legend.title=element_blank()
    ) +
    ggtitle(
      paste0(
        "Clean LSEC RPCA - resolution ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "LSEC_clean_RPCA_UMAP_res",
        format(res, nsmall=1),
        "_v6.7.4.pdf"
      )
    ),
    p,
    width=8.5,
    height=7
  )
}

# ---------------------------------------------------------
# RNA marker panels
# ---------------------------------------------------------

DefaultAssay(lsec) <- "RNA"

lsec <- NormalizeData(
  lsec,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

marker_panel <- c(

  # Core LSEC
  "Clec4g","Stab1","Stab2","Lyve1",
  "Fcgr2b","Mrc1","Oit3","Dnase1l3",

  # Endothelial
  "Pecam1","Cdh5","Kdr","Klf2","Rgcc",

  # Zonation / angiocrine
  "Wnt2","Wnt9b","Rspo3","Bmp2",

  # Capillarization / vascularization
  "Cd34","Vwf","Plvap","Emcn","Esm1",

  # Activation / inflammation
  "Icam1","Vcam1","Sele","Ccl2",
  "Cxcl9","Cxcl10","Nfkbia","Socs3",

  # IFN
  "Isg15","Ifit1","Ifit2","Ifit3",
  "Irf7","Stat1",

  # angiogenic
  "Apln","Aplnr","Angpt2",

  # residual contamination check
  "Ptprc","Lyz2","Adgre1",
  "Col1a1","Pdgfrb","Rgs5",
  "Krt19","Epcam",
  "Alb","Ttr"
)

marker_present <- marker_panel[
  marker_panel %in% rownames(lsec)
]

writeLines(
  marker_present,
  file.path(
    TABDIR,
    "LSEC_marker_panel_present_v6.7.4.txt"
  )
)

dotplot_resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.5
)

for (res in dotplot_resolutions) {

  cluster_col <- paste0(
    "LSECclean_res",
    format(res, nsmall=1)
  )

  Idents(lsec) <- factor(
    lsec@meta.data[[cluster_col]]
  )

  p_dot <- DotPlot(
    lsec,
    features=marker_present,
    assay="RNA",
    dot.scale=6
  ) +
    RotatedAxis() +
    theme_classic(base_size=10) +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank()
    ) +
    ggtitle(
      paste0(
        "Clean LSEC markers - resolution ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "LSEC_clean_marker_DotPlot_res",
        format(res, nsmall=1),
        "_v6.7.4.pdf"
      )
    ),
    p_dot,
    width=18,
    height=6.5
  )
}

# ---------------------------------------------------------
# Sample robustness tables
# ---------------------------------------------------------

for (res in resolutions) {

  cluster_col <- paste0(
    "LSECclean_res",
    format(res, nsmall=1)
  )

  clusters <- sort(
    unique(
      as.character(
        lsec@meta.data[[cluster_col]]
      )
    )
  )

  samples <- sort(
    unique(
      as.character(lsec$sample)
    )
  )

  robustness <- do.call(
    rbind,
    lapply(
      clusters,
      function(cl) {

        idx <-
          as.character(
            lsec@meta.data[[cluster_col]]
          ) == cl

        x <- table(
          factor(
            as.character(lsec$sample[idx]),
            levels=samples
          )
        )

        f <- as.numeric(x) / sum(x)

        data.frame(
          cluster=cl,
          n_cells=sum(x),
          dominant_sample=
            samples[which.max(f)],
          max_sample_fraction=max(f),
          n_samples_ge_5pct=
            sum(f >= 0.05),
          n_samples_ge_10pct=
            sum(f >= 0.10),
          stringsAsFactors=FALSE
        )
      }
    )
  )

  write.csv(
    robustness,
    file.path(
      TABDIR,
      paste0(
        "Cluster_sample_robustness_res",
        format(res, nsmall=1),
        "_v6.7.4.csv"
      )
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

saveRDS(
  lsec,
  file.path(
    OBJDIR,
    "Mouse_MASH_LSEC_clean_RPCA_resolution_scan_v6.7.4.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.4.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.4 COMPLETE\n")
cat("Cells:", ncol(lsec), "\n")
cat("Baseline dimensions: 1:20\n")
cat(
  "Resolutions:",
  paste(resolutions, collapse=", "),
  "\n"
)
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
