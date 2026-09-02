suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

VERSION <- "v6.9.3.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.3/objects/",
  "Mouse_MASH_Monocyte_clean_RPCA_resolution_scan_v6.9.3.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte targeted contamination audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "monocyte_clean_res0_4"

if (!(CLUSTER_COL %in% colnames(obj@meta.data))) {
  stop("Missing cluster column: ", CLUSTER_COL)
}

DefaultAssay(obj) <- "RNA"

# =========================================================
# Raw counts
# =========================================================

assay <- obj[["RNA"]]

if (inherits(assay, "Assay5")) {

  layer_names <- Layers(assay)

  count_layers <- grep(
    "^counts",
    layer_names,
    value=TRUE
  )

  if (length(count_layers) == 0) {
    stop("No RNA counts layer.")
  }

  if (length(count_layers) > 1) {
    obj[["RNA"]] <- JoinLayers(
      obj[["RNA"]],
      layers=count_layers,
      new="counts"
    )
  }
}

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

# =========================================================
# Targeted marker sets
# =========================================================

marker_sets <- list(

  Monocyte_identity=c(
    "Lyz2",
    "Ccr2",
    "Fcgr3",
    "Lst1",
    "Tyrobp",
    "Ctss"
  ),

  Hepatocyte_strict=c(
    "Alb",
    "Ttr",
    "Apoa1",
    "Apoa2",
    "Fabp1",
    "Aldob",
    "Hao1",
    "Apom",
    "Ass1",
    "Bhmt",
    "Tdo2",
    "Hal",
    "Pck1",
    "Hsd17b13",
    "Cyp2e1",
    "Cyp3a11"
  ),

  Endothelial_strict=c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Esam",
    "Ramp2"
  ),

  Mesenchymal_strict=c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Pdgfra",
    "Des",
    "Rgs5"
  ),

  Cholangiocyte_epithelial=c(
    "Krt8",
    "Krt18",
    "Krt19",
    "Epcam",
    "Sox9"
  )
)

all_markers <- unique(
  unlist(
    marker_sets,
    use.names=FALSE
  )
)

presence <- data.frame(
  gene=all_markers,
  present=all_markers %in% rownames(counts),
  stringsAsFactors=FALSE
)

write.csv(
  presence,
  file.path(
    TABDIR,
    "Monocyte_targeted_marker_presence_v6.9.3.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Per-cell marker hits
# =========================================================

get_hits <- function(genes) {

  genes <- intersect(
    genes,
    rownames(counts)
  )

  if (length(genes) == 0) {
    return(rep(0L, ncol(obj)))
  }

  as.integer(
    Matrix::colSums(
      counts[
        genes,
        ,
        drop=FALSE
      ] > 0
    )
  )
}

for (nm in names(marker_sets)) {

  obj[[paste0(
    nm,
    "_hits_v6.9.3.1"
  )]] <- get_hits(
    marker_sets[[nm]]
  )
}

# =========================================================
# Cluster-level audit
# =========================================================

clusters <- sort(
  unique(
    as.character(
      obj@meta.data[[CLUSTER_COL]]
    )
  )
)

audit_rows <- list()

for (cl in clusters) {

  cells <- rownames(obj@meta.data)[
    as.character(
      obj@meta.data[[CLUSTER_COL]]
    ) == cl
  ]

  md <- obj@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  audit_rows[[cl]] <- data.frame(
    cluster=cl,
    n_cells=length(cells),

    median_nCount_RNA=
      median(md$nCount_RNA),

    median_nFeature_RNA=
      median(md$nFeature_RNA),

    Monocyte_median_hits=
      median(
        md$Monocyte_identity_hits_v6.9.3.1
      ),

    Monocyte_fraction_ge2=
      mean(
        md$Monocyte_identity_hits_v6.9.3.1 >= 2
      ),

    Hepatocyte_median_hits=
      median(
        md$Hepatocyte_strict_hits_v6.9.3.1
      ),

    Hepatocyte_fraction_ge2=
      mean(
        md$Hepatocyte_strict_hits_v6.9.3.1 >= 2
      ),

    Hepatocyte_fraction_ge4=
      mean(
        md$Hepatocyte_strict_hits_v6.9.3.1 >= 4
      ),

    Endothelial_median_hits=
      median(
        md$Endothelial_strict_hits_v6.9.3.1
      ),

    Endothelial_fraction_ge2=
      mean(
        md$Endothelial_strict_hits_v6.9.3.1 >= 2
      ),

    Mesenchymal_fraction_ge2=
      mean(
        md$Mesenchymal_strict_hits_v6.9.3.1 >= 2
      ),

    Epithelial_fraction_ge2=
      mean(
        md$Cholangiocyte_epithelial_hits_v6.9.3.1 >= 2
      ),

    stringsAsFactors=FALSE
  )
}

audit <- do.call(
  rbind,
  audit_rows
)

rownames(audit) <- NULL

write.csv(
  audit,
  file.path(
    TABDIR,
    "Monocyte_targeted_contamination_by_cluster_v6.9.3.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Sample composition for clusters 5 and 8
# =========================================================

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "orig.ident"
)

sample_col <- sample_candidates[
  sample_candidates %in%
    colnames(obj@meta.data)
][1]

if (is.na(sample_col)) {
  stop("Could not resolve sample column.")
}

target_clusters <- c("5", "8")

target_sample <- do.call(
  rbind,
  lapply(
    target_clusters,
    function(cl) {

      cells <- rownames(obj@meta.data)[
        as.character(
          obj@meta.data[[CLUSTER_COL]]
        ) == cl
      ]

      tab <- as.data.frame(
        table(
          sample=
            obj@meta.data[
              cells,
              sample_col
            ]
        ),
        stringsAsFactors=FALSE
      )

      tab$cluster <- cl

      tab[
        ,
        c(
          "cluster",
          "sample",
          "Freq"
        )
      ]
    }
  )
)

write.csv(
  target_sample,
  file.path(
    TABDIR,
    "Monocyte_target_clusters5_8_by_sample_v6.9.3.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# DotPlot
# =========================================================

panel <- unique(
  unlist(
    marker_sets,
    use.names=FALSE
  )
)

panel <- intersect(
  panel,
  rownames(obj)
)

p <- DotPlot(
  obj,
  features=panel,
  group.by=CLUSTER_COL,
  assay="RNA"
) +
  RotatedAxis() +
  scale_color_gradient2(
    low="#0033FF",
    mid="#FFFFFF",
    high="#FF1A1A",
    midpoint=0
  ) +
  ggtitle(
    "Monocyte res0.4 targeted lineage contamination audit"
  ) +
  theme_classic(base_size=9)

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_res0.4_targeted_contamination_DotPlot_v6.9.3.1.pdf"
  ),
  p,
  width=16,
  height=8
)

# =========================================================
# Feature-hit UMAPs
# =========================================================

hit_cols <- c(
  "Monocyte_identity_hits_v6.9.3.1",
  "Hepatocyte_strict_hits_v6.9.3.1",
  "Endothelial_strict_hits_v6.9.3.1",
  "Mesenchymal_strict_hits_v6.9.3.1",
  "Cholangiocyte_epithelial_hits_v6.9.3.1"
)

pdf(
  file.path(
    FIGDIR,
    "Monocyte_res0.4_targeted_lineage_hit_UMAP_v6.9.3.1.pdf"
  ),
  width=9,
  height=7
)

for (col in hit_cols) {

  p <- FeaturePlot(
    obj,
    features=col,
    reduction="umap",
    pt.size=0.6
  ) +
    ggtitle(col)

  print(p)
}

dev.off()

# =========================================================
# Save diagnostic object
# =========================================================

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_targeted_contamination_audit_v6.9.3.1.rds"
  )
)

# =========================================================
# Terminal summary
# =========================================================

cat("\n=== TARGETED MARKER PRESENCE ===\n")
print(
  presence,
  row.names=FALSE
)

cat("\n=== TARGETED CONTAMINATION BY CLUSTER ===\n")
print(
  audit,
  row.names=FALSE
)

cat("\n=== CLUSTERS 5 AND 8 BY SAMPLE ===\n")
print(
  target_sample,
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.3.1 COMPLETE\n")
cat("Targeted Hepatocyte / vascular contamination audit complete\n")
cat("No cells removed\n")
cat("No annotation changed\n")
cat("res0.4 remains provisional annotation scaffold\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.3.1.txt"
  )
)
