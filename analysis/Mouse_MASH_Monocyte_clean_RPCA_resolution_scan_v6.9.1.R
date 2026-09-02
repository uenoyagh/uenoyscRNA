suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.9.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.0.1/objects/",
  "Mouse_MASH_Monocyte_stringent_lineage_audit_v6.9.0.1.rds"
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
cat("Mouse MASH Monocyte clean RPCA resolution scan\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

parent <- readRDS(INPUT_RDS)
DefaultAssay(parent) <- "RNA"

cat("Input cells:", ncol(parent), "\n")
cat("Input features:", nrow(parent), "\n\n")

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

cat("Sample column:", sample_col, "\n")
cat(
  "Samples:",
  paste(sort(unique(sample_vec)), collapse=", "),
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
      slot="counts"
    )
  }

  counts
}

counts <- get_raw_counts(parent)

if (ncol(counts) != ncol(parent)) {
  stop("Counts/cell mismatch.")
}

# =========================================================
# Rebuild legacy RNA assay per biological sample
#
# This avoids Assay5 layer complications during RPCA.
# =========================================================

old_assay_option <- getOption(
  "Seurat.object.assay.version"
)

options(
  Seurat.object.assay.version="v3"
)

sample_names <- sort(
  unique(sample_vec)
)

obj_list <- list()

for (s in sample_names) {

  cells <- colnames(parent)[
    sample_vec == s
  ]

  md <- parent@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  cnt <- counts[
    ,
    cells,
    drop=FALSE
  ]

  x <- CreateSeuratObject(
    counts=cnt,
    meta.data=md,
    project=paste0(
      "Mouse_MASH_Monocyte_",
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
# RPCA preparation
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

cat(
  "\nIntegration features:",
  length(features),
  "\n"
)

# =========================================================
# RPCA anchors + integration
# =========================================================

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

cat(
  "Integrated cells:",
  ncol(integrated),
  "\n"
)

if (ncol(integrated) != ncol(parent)) {
  stop(
    "Cell count changed during integration: ",
    ncol(parent),
    " -> ",
    ncol(integrated)
  )
}

# =========================================================
# Integrated PCA / UMAP / neighbors
# =========================================================

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
  0.1,
  0.2,
  0.3,
  0.4,
  0.5,
  0.6,
  0.8,
  1.0
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
    "monocyte_res",
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
# Resolution cluster counts
# =========================================================

resolution_counts <- list()

for (i in seq_along(resolutions)) {

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

  resolution_counts[[i]] <- tab[
    ,
    c(
      "resolution",
      "cluster",
      "Freq"
    )
  ]
}

resolution_counts <- do.call(
  rbind,
  resolution_counts
)

rownames(resolution_counts) <- NULL

write.csv(
  resolution_counts,
  file.path(
    TABDIR,
    "Monocyte_RPCA_resolution_cluster_counts_v6.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Sample composition at every resolution
# =========================================================

integrated_sample <- as.character(
  integrated@meta.data[[sample_col]]
)

sample_cluster_rows <- list()
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

    for (s in sample_names) {

      n <- sum(
        integrated_sample[
          match(
            cells,
            rownames(integrated@meta.data)
          )
        ] == s
      )

      sample_cluster_rows[[counter]] <-
        data.frame(
          resolution=res,
          cluster=cl,
          sample=s,
          n_cells=n,
          stringsAsFactors=FALSE
        )

      counter <- counter + 1
    }
  }
}

sample_cluster_table <- do.call(
  rbind,
  sample_cluster_rows
)

write.csv(
  sample_cluster_table,
  file.path(
    TABDIR,
    "Monocyte_RPCA_cluster_sample_counts_v6.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Strict lineage flags by cluster
# =========================================================

flag_cols <- c(
  Neutrophil=
    "Neutrophil_strict_flag_v6.9.0.1",

  Resident_Macrophage=
    "Resident_Macrophage_strict_flag_v6.9.0.1",

  DC=
    "DC_strict_flag_v6.9.0.1",

  Endothelial=
    "Endothelial_strict_flag_v6.9.0.1",

  Mesenchymal=
    "Mesenchymal_strict_flag_v6.9.0.1",

  Epithelial=
    "Epithelial_strict_flag_v6.9.0.1",

  Lymphoid=
    "Lymphoid_strict_flag_v6.9.0.1"
)

missing_flags <- setdiff(
  unname(flag_cols),
  colnames(integrated@meta.data)
)

if (length(missing_flags) > 0) {
  stop(
    "Missing lineage flags after integration: ",
    paste(missing_flags, collapse=", ")
  )
}

flag_cluster_rows <- list()
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

    md <- integrated@meta.data[
      cells,
      ,
      drop=FALSE
    ]

    for (flag_name in names(flag_cols)) {

      flag_col <- flag_cols[[flag_name]]

      flag_cluster_rows[[counter]] <-
        data.frame(
          resolution=res,
          cluster=cl,
          n_cluster=nrow(md),
          lineage_flag=flag_name,
          n_flagged=sum(
            md[[flag_col]],
            na.rm=TRUE
          ),
          fraction_flagged=mean(
            md[[flag_col]],
            na.rm=TRUE
          ),
          stringsAsFactors=FALSE
        )

      counter <- counter + 1
    }
  }
}

flag_cluster_table <- do.call(
  rbind,
  flag_cluster_rows
)

write.csv(
  flag_cluster_table,
  file.path(
    TABDIR,
    "Monocyte_RPCA_strict_lineage_flags_by_cluster_v6.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Diagnostic res0.4 marker audit
#
# res0.4 is diagnostic only, NOT frozen final resolution.
# =========================================================

diagnostic_col <- "monocyte_res0_4"

if (!(diagnostic_col %in% colnames(integrated@meta.data))) {
  stop(
    "Missing diagnostic resolution column: ",
    diagnostic_col
  )
}

DefaultAssay(integrated) <- "RNA"

integrated <- NormalizeData(
  integrated,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

Idents(integrated) <-
  integrated@meta.data[[diagnostic_col]]

markers_res04 <- FindAllMarkers(
  integrated,
  assay="RNA",
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25,
  test.use="wilcox",
  verbose=FALSE
)

if (nrow(markers_res04) > 0) {

  markers_res04 <- markers_res04[
    order(
      as.numeric(
        as.character(
          markers_res04$cluster
        )
      ),
      -markers_res04$avg_log2FC
    ),
    ,
    drop=FALSE
  ]
}

write.csv(
  markers_res04,
  file.path(
    TABDIR,
    "Monocyte_res0.4_markers_v6.9.1.csv"
  ),
  row.names=FALSE
)

top10_res04 <- do.call(
  rbind,
  lapply(
    split(
      markers_res04,
      markers_res04$cluster
    ),
    function(x) {
      head(
        x[
          order(
            -x$avg_log2FC
          ),
          ,
          drop=FALSE
        ],
        10
      )
    }
  )
)

write.csv(
  top10_res04,
  file.path(
    TABDIR,
    "Monocyte_res0.4_TOP10_markers_v6.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# UMAP figures
# =========================================================

p_sample <- DimPlot(
  integrated,
  reduction="umap",
  group.by=sample_col,
  pt.size=0.5
) +
  ggtitle(
    "Mouse MASH Monocyte RPCA - by sample"
  ) +
  theme_classic()

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_RPCA_UMAP_by_sample_v6.9.1.pdf"
  ),
  p_sample,
  width=9,
  height=7
)

pdf(
  file.path(
    FIGDIR,
    "Monocyte_RPCA_UMAP_resolution_scan_v6.9.1.pdf"
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
    pt.size=0.5
  ) +
    ggtitle(
      paste0(
        "Monocyte RPCA resolution ",
        res
      )
    ) +
    theme_classic()

  print(p)
}

dev.off()

# =========================================================
# Diagnostic res0.4 lineage flag UMAPs
# =========================================================

pdf(
  file.path(
    FIGDIR,
    "Monocyte_res0.4_strict_lineage_flag_UMAP_v6.9.1.pdf"
  ),
  width=9,
  height=7
)

for (flag_name in names(flag_cols)) {

  flag_col <- flag_cols[[flag_name]]

  p <- DimPlot(
    integrated,
    reduction="umap",
    group.by=flag_col,
    pt.size=0.6
  ) +
    ggtitle(
      paste0(
        "res0.4 diagnostic - ",
        flag_name,
        " strict flag"
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
    "Mouse_MASH_Monocyte_RPCA_resolution_scan_v6.9.1.rds"
  )
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== RESOLUTION CLUSTER COUNTS ===\n")
print(
  resolution_counts,
  row.names=FALSE
)

cat("\n=== res0.4 CLUSTER x SAMPLE ===\n")

res04_sample <- subset(
  sample_cluster_table,
  resolution == 0.4
)

print(
  res04_sample,
  row.names=FALSE
)

cat("\n=== res0.4 STRICT LINEAGE FLAG ENRICHMENT ===\n")

res04_flags <- subset(
  flag_cluster_table,
  resolution == 0.4
)

res04_flags <- res04_flags[
  order(
    res04_flags$lineage_flag,
    -res04_flags$fraction_flagged
  ),
  ,
  drop=FALSE
]

print(
  res04_flags,
  row.names=FALSE
)

cat("\n=== res0.4 TOP10 MARKERS ===\n")

if (nrow(top10_res04) > 0) {

  print(
    top10_res04[
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
        colnames(top10_res04)
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )
}

cat("\n====================================================\n")
cat("v6.9.1 COMPLETE\n")
cat("Monocyte RPCA resolution scan complete\n")
cat("Input cells:", ncol(parent), "\n")
cat("Integrated cells:", ncol(integrated), "\n")
cat("No cells removed\n")
cat("No annotation frozen\n")
cat("res0.4 is diagnostic only\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.1.txt"
  )
)
