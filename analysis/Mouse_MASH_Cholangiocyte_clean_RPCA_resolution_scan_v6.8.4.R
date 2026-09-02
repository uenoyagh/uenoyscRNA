suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.4"
DIMS_USE <- 1:20

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.3.2/objects/",
  "Mouse_MASH_Cholangiocyte_lineage_clean_raw_v6.8.3.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte clean RPCA resolution scan\n")
cat("Version:", VERSION, "\n")
cat("Dims:", paste(range(DIMS_USE), collapse=":"), "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

clean <- readRDS(INPUT_RDS)
DefaultAssay(clean) <- "RNA"

cat("=== INPUT CLEAN OBJECT ===\n")
cat("Cells:", ncol(clean), "\n")
cat("Features:", nrow(clean), "\n")
cat("Assays:", paste(Assays(clean), collapse=", "), "\n")
cat("RNA assay class:", class(clean[["RNA"]])[1], "\n")
cat("Reductions:", length(Reductions(clean)), "\n")
cat("Graphs:", length(clean@graphs), "\n")

if (ncol(clean) != 15755) {
  stop(
    "Expected 15755 cells, observed ",
    ncol(clean)
  )
}

required_md <- c(
  "sample",
  "condition",
  "Chol_res03_precleanup_v6832",
  "cleanup_status_v6832"
)

missing_md <- setdiff(
  required_md,
  colnames(clean@meta.data)
)

if (length(missing_md) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing_md, collapse=", ")
  )
}

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

clean$sample <- factor(
  as.character(clean$sample),
  levels=sample_order
)

# ---------------------------------------------------------
# Track old cluster 4 explicitly
# ---------------------------------------------------------

clean$precleanup_QC_watch_v684 <- ifelse(
  as.character(
    clean$Chol_res03_precleanup_v6832
  ) == "4",
  "Precleanup_cluster4_QC_watch",
  "Other_clean_Cholangiocyte"
)

cat("\n=== PRECLEANUP QC-WATCH CELLS ===\n")
print(
  table(
    clean$precleanup_QC_watch_v684
  )
)

cat("\n=== PRECLEANUP QC-WATCH BY SAMPLE ===\n")
print(
  table(
    clean$sample,
    clean$precleanup_QC_watch_v684
  )
)

# ---------------------------------------------------------
# Extract raw RNA counts
# ---------------------------------------------------------

counts <- GetAssayData(
  clean,
  assay="RNA",
  layer="counts"
)

meta <- clean@meta.data

cat("\n=== RAW COUNTS ===\n")
cat(
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)

# ---------------------------------------------------------
# Rebuild six independent legacy-RNA sample objects
#
# This avoids inherited Assay5 split-layer complexity.
# ---------------------------------------------------------

options(
  Seurat.object.assay.version="v3"
)

obj_list <- vector(
  "list",
  length(sample_order)
)

names(obj_list) <- sample_order

for (s in sample_order) {

  cells_s <- rownames(meta)[
    as.character(meta$sample) == s
  ]

  if (length(cells_s) == 0) {
    stop("No cells for sample: ", s)
  }

  counts_s <- counts[
    ,
    cells_s,
    drop=FALSE
  ]

  meta_s <- meta[
    cells_s,
    ,
    drop=FALSE
  ]

  x <- CreateSeuratObject(
    counts=counts_s,
    assay="RNA",
    meta.data=meta_s,
    project=paste0(
      "Mouse_MASH_Chol_",
      s
    )
  )

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

  obj_list[[s]] <- x

  cat(
    "Prepared ",
    s,
    ": ",
    ncol(x),
    " cells; assay=",
    class(x[["RNA"]])[1],
    "\n",
    sep=""
  )
}

# ---------------------------------------------------------
# Integration features
# ---------------------------------------------------------

features <- SelectIntegrationFeatures(
  object.list=obj_list,
  nfeatures=3000
)

cat("\nIntegration features:", length(features), "\n")

for (s in names(obj_list)) {

  x <- obj_list[[s]]

  x <- ScaleData(
    x,
    features=features,
    verbose=FALSE
  )

  x <- RunPCA(
    x,
    features=features,
    npcs=30,
    verbose=FALSE
  )

  obj_list[[s]] <- x
}

# ---------------------------------------------------------
# RPCA integration
# ---------------------------------------------------------

cat("\n=== FIND RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list=obj_list,
  anchor.features=features,
  reduction="rpca",
  dims=DIMS_USE
)

cat("\n=== INTEGRATE DATA ===\n")

integrated <- IntegrateData(
  anchorset=anchors,
  dims=DIMS_USE
)

cat("\nIntegrated cells:", ncol(integrated), "\n")

if (ncol(integrated) != 15755) {
  stop(
    "Cell count changed during integration: ",
    ncol(integrated)
  )
}

# ---------------------------------------------------------
# Integrated PCA
# ---------------------------------------------------------

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(
  integrated,
  verbose=FALSE
)

integrated <- RunPCA(
  integrated,
  npcs=30,
  verbose=FALSE
)

p_elbow <- ElbowPlot(
  integrated,
  ndims=30
) +
  ggtitle(
    "Clean Cholangiocyte PCA elbow after RPCA"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_PCA_elbow_v6.8.4.pdf"
  ),
  p_elbow,
  width=7,
  height=5
)

# ---------------------------------------------------------
# UMAP and neighbors
# ---------------------------------------------------------

integrated <- RunUMAP(
  integrated,
  reduction="pca",
  dims=DIMS_USE,
  seed.use=20260902,
  verbose=FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction="pca",
  dims=DIMS_USE,
  verbose=FALSE
)

# ---------------------------------------------------------
# Sample and condition UMAP
# ---------------------------------------------------------

p_sample <- DimPlot(
  integrated,
  reduction="umap",
  group.by="sample",
  pt.size=0.2
) +
  ggtitle(
    "Clean Cholangiocyte RPCA UMAP by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_RPCA_UMAP_by_sample_v6.8.4.pdf"
  ),
  p_sample,
  width=8,
  height=7
)

p_condition <- DimPlot(
  integrated,
  reduction="umap",
  group.by="condition",
  pt.size=0.2
) +
  ggtitle(
    "Clean Cholangiocyte RPCA UMAP by condition"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_RPCA_UMAP_by_condition_v6.8.4.pdf"
  ),
  p_condition,
  width=8,
  height=7
)

# ---------------------------------------------------------
# Track precleanup cluster 4
# ---------------------------------------------------------

p_qcwatch <- DimPlot(
  integrated,
  reduction="umap",
  group.by="precleanup_QC_watch_v684",
  pt.size=0.25
) +
  ggtitle(
    "Clean RPCA: location of precleanup cluster 4 QC-watch cells"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_RPCA_precleanup_cluster4_QCwatch_v6.8.4.pdf"
  ),
  p_qcwatch,
  width=8,
  height=7
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
  0.6,
  0.8
)

for (res in resolutions) {

  integrated <- FindClusters(
    integrated,
    resolution=res,
    algorithm=1,
    random.seed=20260902,
    verbose=FALSE
  )

  cluster_col <- paste0(
    "integrated_snn_res.",
    res
  )

  if (!cluster_col %in% colnames(integrated@meta.data)) {
    stop(
      "Missing cluster column: ",
      cluster_col
    )
  }

  cat(
    "\n=== RESOLUTION ",
    res,
    " ===\n",
    sep=""
  )

  print(
    sort(
      table(
        integrated@meta.data[[cluster_col]]
      ),
      decreasing=TRUE
    )
  )

  # Cluster counts
  write.csv(
    as.data.frame(
      table(
        cluster=
          integrated@meta.data[[cluster_col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cholangiocyte_clean_cluster_counts_res",
        res,
        "_v6.8.4.csv"
      )
    ),
    row.names=FALSE
  )

  # Cluster x sample
  tab_sample <- table(
    cluster=
      integrated@meta.data[[cluster_col]],
    sample=
      integrated$sample
  )

  write.csv(
    as.data.frame.matrix(tab_sample),
    file.path(
      TABDIR,
      paste0(
        "Cholangiocyte_clean_cluster_by_sample_res",
        res,
        "_v6.8.4.csv"
      )
    )
  )

  # Cluster x old QC-watch
  tab_qc <- table(
    cluster=
      integrated@meta.data[[cluster_col]],
    QC_watch=
      integrated$precleanup_QC_watch_v684
  )

  write.csv(
    as.data.frame.matrix(tab_qc),
    file.path(
      TABDIR,
      paste0(
        "Cholangiocyte_clean_cluster_by_precleanup_QCwatch_res",
        res,
        "_v6.8.4.csv"
      )
    )
  )

  # UMAP
  p <- DimPlot(
    integrated,
    reduction="umap",
    group.by=cluster_col,
    label=TRUE,
    repel=TRUE,
    pt.size=0.2
  ) +
    ggtitle(
      paste0(
        "Clean Cholangiocyte RPCA, res ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Cholangiocyte_clean_RPCA_UMAP_res",
        res,
        "_v6.8.4.pdf"
      )
    ),
    p,
    width=8,
    height=7
  )
}

# ---------------------------------------------------------
# Marker panels
# ---------------------------------------------------------

DefaultAssay(integrated) <- "RNA"

marker_sets <- list(

  Biliary_identity = c(
    "Krt19","Krt7","Krt8","Krt18",
    "Epcam","Sox9","Muc1",
    "Hnf1b","Cftr","Slc4a2"
  ),

  Ductular_reactive = c(
    "Spp1","Mmp7","Krt23",
    "Tacstd2","Prom1","Klf5"
  ),

  Inflammatory = c(
    "Icam1","Vcam1",
    "Cxcl1","Cxcl2","Ccl2",
    "Il6","Nfkbia","Socs3"
  ),

  Remodeling = c(
    "Tgfb1","Tgfb2",
    "Ccn2","Jag1",
    "Serpine1","Thbs1"
  ),

  Stress_IEG = c(
    "Cdkn1a","Cdkn2a",
    "Fos","Jun","Atf3",
    "Ddit4","Gadd45a"
  ),

  Cycling = c(
    "Mki67","Top2a",
    "Birc5","Ube2c",
    "Cenpf","Pcna","Stmn1"
  ),

  Ciliated = c(
    "Cfap73","Cfap44",
    "Lrrc23","Drc1",
    "Dnali1","Mns1"
  ),

  Tuft_like = c(
    "Pou2f3","Trpm5",
    "Gnat3"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1",
    "Fabp1","Cps1","Hnf4a"
  ),

  Myeloid = c(
    "Ptprc","Lyz2",
    "Tyrobp","Aif1",
    "Adgre1"
  ),

  Endothelial = c(
    "Pecam1","Cdh5",
    "Kdr","Stab2"
  ),

  Mesenchymal = c(
    "Col1a1","Col3a1",
    "Pdgfra","Pdgfrb",
    "Lrat","Rbp1"
  ),

  Neutrophil = c(
    "S100a8","S100a9",
    "Retnlg","Mmp8"
  ),

  Lymphoid = c(
    "Cd3d","Cd3e",
    "Nkg7","Cd79a"
  )
)

marker_sets <- lapply(
  marker_sets,
  function(x) {
    intersect(
      x,
      rownames(integrated)
    )
  }
)

marker_features <- unique(
  unlist(
    marker_sets,
    use.names=FALSE
  )
)

# ---------------------------------------------------------
# DotPlots at diagnostic resolutions
# ---------------------------------------------------------

for (res in c(0.2, 0.3, 0.4, 0.5, 0.6)) {

  cluster_col <- paste0(
    "integrated_snn_res.",
    res
  )

  Idents(integrated) <-
    integrated@meta.data[[cluster_col]]

  p_dot <- DotPlot(
    integrated,
    features=marker_features,
    assay="RNA",
    dot.scale=6
  ) +
    RotatedAxis() +
    theme_classic(base_size=8) +
    theme(
      axis.title=element_blank()
    ) +
    ggtitle(
      paste0(
        "Clean Cholangiocyte marker audit, res ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Cholangiocyte_clean_marker_DotPlot_res",
        res,
        "_v6.8.4.pdf"
      )
    ),
    p_dot,
    width=23,
    height=7
  )
}

# ---------------------------------------------------------
# QC-watch enrichment across resolutions
# ---------------------------------------------------------

qcwatch_summary <- list()

for (res in resolutions) {

  cluster_col <- paste0(
    "integrated_snn_res.",
    res
  )

  clusters <- as.character(
    integrated@meta.data[[cluster_col]]
  )

  qc <- integrated$precleanup_QC_watch_v684 ==
    "Precleanup_cluster4_QC_watch"

  df <- do.call(
    rbind,
    lapply(
      sort(
        unique(
          as.numeric(clusters)
        )
      ),
      function(cl) {

        idx <- clusters == as.character(cl)

        data.frame(
          resolution=res,
          cluster=cl,
          n_cells=sum(idx),
          precleanup_cluster4_n=
            sum(qc[idx]),
          precleanup_cluster4_fraction=
            mean(qc[idx]),
          stringsAsFactors=FALSE
        )
      }
    )
  )

  qcwatch_summary[[as.character(res)]] <- df
}

qcwatch_summary <- do.call(
  rbind,
  qcwatch_summary
)

rownames(qcwatch_summary) <- NULL

write.csv(
  qcwatch_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_precleanup_cluster4_enrichment_across_resolutions_v6.8.4.csv"
  ),
  row.names=FALSE
)

cat("\n=== QC-WATCH ENRICHMENT: TOP CLUSTERS ===\n")

for (res in resolutions) {

  x <- qcwatch_summary[
    qcwatch_summary$resolution == res,
    ,
    drop=FALSE
  ]

  x <- x[
    order(
      -x$precleanup_cluster4_fraction,
      -x$precleanup_cluster4_n
    ),
    ,
    drop=FALSE
  ]

  cat(
    "\nResolution ",
    res,
    "\n",
    sep=""
  )

  print(
    head(x, 5),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Save clean integrated object
# ---------------------------------------------------------

saveRDS(
  integrated,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_clean_RPCA_resolution_scan_v6.8.4.rds"
  ),
  compress=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte clean RPCA resolution scan v6.8.4",
  "",
  paste0("- Input clean cells: ", ncol(clean)),
  paste0("- Integrated cells: ", ncol(integrated)),
  "- Source clean baseline: v6.8.3.2",
  "- Removed contaminating precleanup clusters: 3, 9, 10, 11, 12",
  "- Precleanup cluster 4 was retained as QC-watch.",
  "- Precleanup cluster 4 cells remain explicitly tracked.",
  "- RPCA dimensions: 1:20",
  "- Resolutions scanned: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8",
  "- No final Cholangiocyte state annotation assigned."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_clean_RPCA_resolution_scan_summary_v6.8.4.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.4.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.4 COMPLETE\n")
cat("Input cells:", ncol(clean), "\n")
cat("Integrated cells:", ncol(integrated), "\n")
cat("RPCA dims: 1:20\n")
cat("Precleanup cluster 4 retained and tracked\n")
cat("No final state annotation assigned\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
