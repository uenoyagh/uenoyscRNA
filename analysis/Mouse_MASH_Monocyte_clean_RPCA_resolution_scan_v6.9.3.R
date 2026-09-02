suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.9.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.2/objects/",
  "Mouse_MASH_Monocyte_lineage_clean_raw_v6.9.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte CLEAN RPCA resolution scan\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

parent <- readRDS(INPUT_RDS)
DefaultAssay(parent) <- "RNA"

cat("Input cells:", ncol(parent), "\n")
cat("Input features:", nrow(parent), "\n\n")

if (ncol(parent) != 3490) {
  warning(
    "Expected 3490 cells after v6.9.2, observed: ",
    ncol(parent)
  )
}

# =========================================================
# Sample column
# =========================================================

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "orig.ident"
)

sample_col <- sample_candidates[
  sample_candidates %in%
    colnames(parent@meta.data)
][1]

if (is.na(sample_col)) {
  stop("Could not resolve sample column.")
}

sample_vec <- as.character(
  parent@meta.data[[sample_col]]
)

sample_names <- sort(
  unique(sample_vec)
)

cat("Sample column:", sample_col, "\n")
cat(
  "Samples:",
  paste(sample_names, collapse=", "),
  "\n\n"
)

# =========================================================
# Raw counts
# =========================================================

get_raw_counts <- function(obj) {

  assay <- obj[["RNA"]]

  if (inherits(assay, "Assay5")) {

    layer_names <- Layers(assay)

    count_layers <- grep(
      "^counts",
      layer_names,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop(
        "No RNA counts layer found: ",
        paste(layer_names, collapse=", ")
      )
    }

    if (length(count_layers) > 1) {

      obj[["RNA"]] <- JoinLayers(
        obj[["RNA"]],
        layers=count_layers,
        new="counts"
      )
    }

    counts <- GetAssayData(
      obj,
      assay="RNA",
      layer="counts"
    )

  } else {

    counts <- GetAssayData(
      obj,
      assay="RNA",
      layer="counts"
    )
  }

  counts
}

counts <- get_raw_counts(parent)

if (ncol(counts) != ncol(parent)) {
  stop("Counts/cell mismatch.")
}

# =========================================================
# Rebuild legacy RNA object per biological sample
# =========================================================

old_assay_option <- getOption(
  "Seurat.object.assay.version"
)

options(
  Seurat.object.assay.version="v3"
)

obj_list <- list()

for (s in sample_names) {

  cells <- colnames(parent)[
    sample_vec == s
  ]

  cnt <- counts[
    ,
    cells,
    drop=FALSE
  ]

  md <- parent@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  x <- CreateSeuratObject(
    counts=cnt,
    meta.data=md,
    project=paste0(
      "Mouse_MASH_Monocyte_clean_",
      s
    )
  )

  x <- NormalizeData(
    x,
    normalization.method="LogNormalize",
    scale.factor=10000,
    verbose=FALSE
  )

  x <- FindVariableFeatures(
    x,
    selection.method="vst",
    nfeatures=2000,
    verbose=FALSE
  )

  obj_list[[s]] <- x

  cat(
    "Prepared:",
    s,
    "| cells:",
    ncol(x),
    "\n"
  )
}

options(
  Seurat.object.assay.version=
    old_assay_option
)

# =========================================================
# RPCA integration
# =========================================================

features <- SelectIntegrationFeatures(
  object.list=obj_list,
  nfeatures=2500
)

for (s in names(obj_list)) {

  obj_list[[s]] <- ScaleData(
    obj_list[[s]],
    features=features,
    verbose=FALSE
  )

  obj_list[[s]] <- RunPCA(
    obj_list[[s]],
    features=features,
    npcs=30,
    verbose=FALSE
  )
}

anchors <- FindIntegrationAnchors(
  object.list=obj_list,
  anchor.features=features,
  reduction="rpca",
  dims=1:20
)

integrated <- IntegrateData(
  anchorset=anchors,
  dims=1:20
)

if (ncol(integrated) != ncol(parent)) {
  stop(
    "Cell count changed during integration: ",
    ncol(parent),
    " -> ",
    ncol(integrated)
  )
}

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

integrated <- RunUMAP(
  integrated,
  reduction="pca",
  dims=1:20,
  seed.use=20260902,
  verbose=FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction="pca",
  dims=1:20,
  verbose=FALSE
)

# =========================================================
# Resolution scan
# =========================================================

resolutions <- c(
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8
)

resolution_cols <- character()

for (res in resolutions) {

  integrated <- FindClusters(
    integrated,
    resolution=res,
    algorithm=1,
    random.seed=20260902,
    verbose=FALSE
  )

  col_name <- paste0(
    "monocyte_clean_res",
    gsub(
      "\\.",
      "_",
      format(
        res,
        trim=TRUE,
        scientific=FALSE
      )
    )
  )

  integrated[[col_name]] <-
    as.character(
      Idents(integrated)
    )

  resolution_cols <- c(
    resolution_cols,
    col_name
  )

  cat(
    "Resolution",
    res,
    ":",
    length(
      unique(
        integrated@meta.data[[col_name]]
      )
    ),
    "clusters\n"
  )
}

# =========================================================
# Cluster counts
# =========================================================

resolution_counts <- do.call(
  rbind,
  lapply(
    seq_along(resolutions),
    function(i) {

      res <- resolutions[[i]]
      col_name <- resolution_cols[[i]]

      tab <- as.data.frame(
        table(
          cluster=
            integrated@meta.data[[col_name]]
        ),
        stringsAsFactors=FALSE
      )

      tab$resolution <- res

      tab[
        ,
        c(
          "resolution",
          "cluster",
          "Freq"
        )
      ]
    }
  )
)

rownames(resolution_counts) <- NULL

write.csv(
  resolution_counts,
  file.path(
    TABDIR,
    "Monocyte_clean_RPCA_resolution_cluster_counts_v6.9.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cluster x sample counts
# =========================================================

integrated_sample <- as.character(
  integrated@meta.data[[sample_col]]
)

cluster_sample_rows <- list()
counter <- 1

for (i in seq_along(resolutions)) {

  res <- resolutions[[i]]
  col_name <- resolution_cols[[i]]

  clusters <- sort(
    unique(
      integrated@meta.data[[col_name]]
    )
  )

  for (cl in clusters) {

    cells <- rownames(integrated@meta.data)[
      integrated@meta.data[[col_name]] == cl
    ]

    sample_sub <- integrated_sample[
      match(
        cells,
        rownames(integrated@meta.data)
      )
    ]

    for (s in sample_names) {

      cluster_sample_rows[[counter]] <-
        data.frame(
          resolution=res,
          cluster=cl,
          sample=s,
          n_cells=sum(sample_sub == s),
          stringsAsFactors=FALSE
        )

      counter <- counter + 1
    }
  }
}

cluster_sample_table <- do.call(
  rbind,
  cluster_sample_rows
)

write.csv(
  cluster_sample_table,
  file.path(
    TABDIR,
    "Monocyte_clean_RPCA_cluster_sample_counts_v6.9.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Marker audit at res0.4 and res0.5
#
# These are diagnostic only.
# =========================================================

DefaultAssay(integrated) <- "RNA"

integrated <- NormalizeData(
  integrated,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

marker_resolutions <- c(
  "0.4"="monocyte_clean_res0_4",
  "0.5"="monocyte_clean_res0_5"
)

top_marker_list <- list()

for (res_name in names(marker_resolutions)) {

  col_name <- marker_resolutions[[res_name]]

  if (!(col_name %in% colnames(integrated@meta.data))) {
    stop("Missing cluster column: ", col_name)
  }

  Idents(integrated) <-
    integrated@meta.data[[col_name]]

  markers <- FindAllMarkers(
    integrated,
    assay="RNA",
    only.pos=TRUE,
    min.pct=0.10,
    logfc.threshold=0.25,
    test.use="wilcox",
    verbose=FALSE
  )

  write.csv(
    markers,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_clean_res",
        res_name,
        "_markers_v6.9.3.csv"
      )
    ),
    row.names=FALSE
  )

  top10 <- do.call(
    rbind,
    lapply(
      split(
        markers,
        markers$cluster
      ),
      function(x) {

        x <- x[
          order(
            -x$avg_log2FC
          ),
          ,
          drop=FALSE
        ]

        head(x, 10)
      }
    )
  )

  rownames(top10) <- NULL

  write.csv(
    top10,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_clean_res",
        res_name,
        "_TOP10_markers_v6.9.3.csv"
      )
    ),
    row.names=FALSE
  )

  top_marker_list[[res_name]] <- top10
}

# =========================================================
# Prespecified biological marker panel
# =========================================================

marker_panel <- c(
  # Monocyte identity / classical
  "Lyz2",
  "Ccr2",
  "Fcgr3",
  "Lst1",
  "Tyrobp",
  "Ctss",
  "Sell",

  # inflammatory / neutrophil-like
  "Mmp8",
  "Chil3",
  "Cd177",
  "S100a8",
  "S100a9",
  "Il1b",
  "Tnf",

  # inflammatory activation
  "Nos2",
  "Cxcl9",
  "Cxcl10",
  "Saa3",

  # IFN response
  "Ifit1",
  "Ifit2",
  "Ifit3",
  "Rsad2",
  "Isg15",
  "Cmpk2",

  # monocyte-to-macrophage / lipid-remodeling
  "Ms4a7",
  "Lpl",
  "Mmp12",
  "Gpr84",
  "Adgre1",
  "Mertk",
  "Clec4f",
  "Timd4",
  "Vsig4",

  # antigen presentation
  "H2-Ab1",
  "H2-Aa",
  "Cd74",

  # repair / resolution
  "Mrc1",
  "Cd163",
  "Folr2",
  "Il10",
  "Tgfb1"
)

marker_panel <- intersect(
  marker_panel,
  rownames(integrated)
)

for (res_name in names(marker_resolutions)) {

  col_name <- marker_resolutions[[res_name]]

  p <- DotPlot(
    integrated,
    features=marker_panel,
    group.by=col_name,
    assay="RNA"
  ) +
    RotatedAxis() +
    ggtitle(
      paste0(
        "Monocyte clean RPCA res",
        res_name,
        " marker audit"
      )
    ) +
    scale_color_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(base_size=9)

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "Monocyte_clean_res",
        res_name,
        "_marker_DotPlot_v6.9.3.pdf"
      )
    ),
    p,
    width=15,
    height=7
  )
}

# =========================================================
# UMAPs
# =========================================================

p_sample <- DimPlot(
  integrated,
  reduction="umap",
  group.by=sample_col,
  pt.size=0.55
) +
  ggtitle(
    "Mouse MASH lineage-clean Monocyte RPCA - by sample"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_clean_RPCA_UMAP_by_sample_v6.9.3.pdf"
  ),
  p_sample,
  width=9,
  height=7
)

pdf(
  file.path(
    FIGDIR,
    "Monocyte_clean_RPCA_UMAP_resolution_scan_v6.9.3.pdf"
  ),
  width=9,
  height=7
)

for (i in seq_along(resolutions)) {

  res <- resolutions[[i]]
  col_name <- resolution_cols[[i]]

  p <- DimPlot(
    integrated,
    reduction="umap",
    group.by=col_name,
    label=TRUE,
    repel=TRUE,
    pt.size=0.55
  ) +
    ggtitle(
      paste0(
        "Lineage-clean Monocyte RPCA resolution ",
        res
      )
    ) +
    theme_classic()

  print(p)
}

dev.off()

# =========================================================
# Save object
# =========================================================

saveRDS(
  integrated,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_clean_RPCA_resolution_scan_v6.9.3.rds"
  )
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== CLEAN RESOLUTION CLUSTER COUNTS ===\n")
print(
  resolution_counts,
  row.names=FALSE
)

cat("\n=== CLEAN res0.4 CLUSTER x SAMPLE ===\n")
print(
  subset(
    cluster_sample_table,
    resolution == 0.4
  ),
  row.names=FALSE
)

cat("\n=== CLEAN res0.5 CLUSTER x SAMPLE ===\n")
print(
  subset(
    cluster_sample_table,
    resolution == 0.5
  ),
  row.names=FALSE
)

for (res_name in c("0.4", "0.5")) {

  cat(
    "\n=== CLEAN res",
    res_name,
    " TOP10 MARKERS ===\n",
    sep=""
  )

  top10 <- top_marker_list[[res_name]]

  print(
    top10[
      ,
      intersect(
        c(
          "cluster",
          "gene",
          "avg_log2FC",
          "pct.1",
          "pct.2",
          "p_val_adj"
        ),
        colnames(top10)
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )
}

cat("\n====================================================\n")
cat("v6.9.3 COMPLETE\n")
cat("Lineage-clean Monocyte RPCA resolution scan complete\n")
cat("Input cells:", ncol(parent), "\n")
cat("Integrated cells:", ncol(integrated), "\n")
cat("No cells removed\n")
cat("No annotation frozen\n")
cat("res0.4 / res0.5 are diagnostic candidates\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.3.txt"
  )
)
