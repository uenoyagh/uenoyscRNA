suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.0/objects/",
  "Mouse_MASH_Cholangiocyte_parent_raw_clean_v6.8.0.rds"
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
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

cat("=== INPUT OBJECT ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")
cat("Reductions:", length(Reductions(obj)), "\n")
cat("Graphs:", length(obj@graphs), "\n")

required_md <- c(
  "sample",
  "condition"
)

miss <- setdiff(
  required_md,
  colnames(obj@meta.data)
)

if (length(miss) > 0) {
  stop(
    "Missing metadata: ",
    paste(miss, collapse=", ")
  )
}

cat("\n=== SAMPLE COUNTS ===\n")
print(table(obj$sample))

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

obj$sample <- factor(
  obj$sample,
  levels=sample_order
)

# ---------------------------------------------------------
# Split by biological sample
# ---------------------------------------------------------

obj_list <- SplitObject(
  obj,
  split.by="sample"
)

obj_list <- obj_list[
  sample_order
]

cat("\n=== SAMPLE OBJECTS ===\n")

for (nm in names(obj_list)) {
  cat(nm, ":", ncol(obj_list[[nm]]), "cells\n")
}

# ---------------------------------------------------------
# Sample-wise preprocessing
# ---------------------------------------------------------

obj_list <- lapply(
  obj_list,
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

features <- SelectIntegrationFeatures(
  object.list=obj_list,
  nfeatures=3000
)

cat("\nIntegration features:", length(features), "\n")

obj_list <- lapply(
  obj_list,
  function(x) {

    x <- ScaleData(
      x,
      features=features,
      verbose=FALSE
    )

    x <- RunPCA(
      x,
      features=features,
      npcs=40,
      verbose=FALSE
    )

    x
  }
)

# ---------------------------------------------------------
# RPCA anchors
# ---------------------------------------------------------

cat("\n=== FIND RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list=obj_list,
  anchor.features=features,
  reduction="rpca",
  dims=1:30
)

cat("\n=== INTEGRATE DATA ===\n")

integrated <- IntegrateData(
  anchorset=anchors,
  dims=1:30
)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(
  integrated,
  verbose=FALSE
)

integrated <- RunPCA(
  integrated,
  npcs=40,
  verbose=FALSE
)

# ---------------------------------------------------------
# PCA elbow
# ---------------------------------------------------------

p_elbow <- ElbowPlot(
  integrated,
  ndims=40
) +
  ggtitle(
    "Cholangiocyte PCA elbow after RPCA integration"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_PCA_elbow_v6.8.1.pdf"
  ),
  p_elbow,
  width=7,
  height=5
)

# ---------------------------------------------------------
# UMAP / neighbors
# ---------------------------------------------------------

integrated <- RunUMAP(
  integrated,
  reduction="pca",
  dims=1:30,
  seed.use=20260902,
  verbose=FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction="pca",
  dims=1:30,
  verbose=FALSE
)

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

  col <- paste0(
    "integrated_snn_res.",
    res
  )

  if (!col %in% colnames(integrated@meta.data)) {
    stop("Missing cluster column: ", col)
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
        integrated@meta.data[[col]]
      ),
      decreasing=TRUE
    )
  )

  write.csv(
    as.data.frame(
      table(
        cluster=
          integrated@meta.data[[col]]
      )
    ),
    file.path(
      TABDIR,
      paste0(
        "Cluster_counts_res",
        res,
        "_v6.8.1.csv"
      )
    ),
    row.names=FALSE
  )

  tab_sample <- as.data.frame(
    table(
      cluster=
        integrated@meta.data[[col]],
      sample=
        integrated$sample
    )
  )

  write.csv(
    tab_sample,
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_sample_res",
        res,
        "_v6.8.1.csv"
      )
    ),
    row.names=FALSE
  )

  p <- DimPlot(
    integrated,
    reduction="umap",
    group.by=col,
    label=TRUE,
    repel=TRUE,
    pt.size=0.2
  ) +
    ggtitle(
      paste0(
        "Cholangiocyte clean RPCA, res ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Cholangiocyte_clean_RPCA_UMAP_res",
        res,
        "_v6.8.1.pdf"
      )
    ),
    p,
    width=8,
    height=7
  )
}

# ---------------------------------------------------------
# UMAP by sample / condition
# ---------------------------------------------------------

p_sample <- DimPlot(
  integrated,
  reduction="umap",
  group.by="sample",
  pt.size=0.2
) +
  ggtitle(
    "Cholangiocyte clean RPCA UMAP by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_RPCA_UMAP_by_sample_v6.8.1.pdf"
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
    "Cholangiocyte clean RPCA UMAP by condition"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_clean_RPCA_UMAP_by_condition_v6.8.1.pdf"
  ),
  p_condition,
  width=8,
  height=7
)

# ---------------------------------------------------------
# Marker audit at selected diagnostic resolutions
# ---------------------------------------------------------

DefaultAssay(integrated) <- "RNA"

marker_panel <- c(
  # biliary identity
  "Krt19","Krt7","Krt8","Krt18",
  "Epcam","Sox9","Muc1","Hnf1b",
  "Cftr","Slc4a2",

  # reactive / ductular
  "Spp1","Mmp7","Krt23",
  "Tacstd2","Prom1",

  # inflammatory
  "Icam1","Vcam1",
  "Cxcl1","Cxcl2","Ccl2",
  "Il6","Nfkbia","Socs3",

  # fibrogenic / remodeling
  "Tgfb1","Tgfb2",
  "Ccn2","Jag1",
  "Serpine1","Thbs1",

  # stress / senescence
  "Cdkn1a","Cdkn2a",
  "Fos","Jun","Atf3",

  # cycling
  "Mki67","Top2a",
  "Birc5","Ube2c","Cenpf",

  # contamination
  "Alb","Ttr","Cps1","Hnf4a",
  "Ptprc","Lyz2","Tyrobp",
  "Pecam1","Cdh5","Stab2",
  "Col1a1","Pdgfrb","Rbp1",
  "S100a8","S100a9","Nkg7"
)

marker_panel <- intersect(
  marker_panel,
  rownames(integrated)
)

for (res in c(0.2, 0.3, 0.4, 0.6)) {

  col <- paste0(
    "integrated_snn_res.",
    res
  )

  Idents(integrated) <-
    integrated@meta.data[[col]]

  p_dot <- DotPlot(
    integrated,
    features=marker_panel,
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
        "Cholangiocyte marker audit, res ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Cholangiocyte_marker_DotPlot_res",
        res,
        "_v6.8.1.pdf"
      )
    ),
    p_dot,
    width=21,
    height=7
  )
}

# ---------------------------------------------------------
# Save integrated object
# ---------------------------------------------------------

saveRDS(
  integrated,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_RPCA_resolution_scan_v6.8.1.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.1 COMPLETE\n")
cat("Input Cholangiocyte cells:", ncol(obj), "\n")
cat("Integrated cells:", ncol(integrated), "\n")
cat("Resolutions:", paste(resolutions, collapse=", "), "\n")
cat("No annotation fixed yet\n")
cat("Cycling rescue cells were NOT added\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
