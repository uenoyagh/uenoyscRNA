suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.7.2"
RES_COL <- "LSEC_res0.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.1.1/objects/",
  "Mouse_MASH_endothelial_clean_RPCA_resolution_scan_v6.7.1.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("============================================\n")
cat("Mouse MASH LSEC expanded marker audit\n")
cat("Version:", VERSION, "\n")
cat("Fixed resolution: 0.6\n")
cat("============================================\n\n")

obj <- readRDS(INPUT_RDS)

if (!RES_COL %in% colnames(obj@meta.data)) {
  stop("Missing cluster column: ", RES_COL)
}

if (!"umapLSEC" %in% Reductions(obj)) {
  stop("Missing reduction: umapLSEC")
}

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

Idents(obj) <- factor(obj@meta.data[[RES_COL]])

cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n\n")

cat("=== CLUSTER COUNTS ===\n")
print(sort(table(Idents(obj)), decreasing=TRUE))

marker_sets <- list(

  Pan_endothelial = c(
    "Pecam1","Cdh5","Kdr","Klf2","Klf4","Eng","Esam","Rgcc"
  ),

  LSEC_identity = c(
    "Clec4g","Clec4a1","Stab1","Stab2","Lyve1",
    "Fcgr2b","Mrc1","F8","Dnase1l3","Oit3"
  ),

  Capillarization = c(
    "Cd34","Vwf","Plvap","Emcn","Esm1"
  ),

  LSEC_zonation = c(
    "Wnt2","Wnt9b","Rspo3","Bmp2"
  ),

  Portal_arterial_EC = c(
    "Efnb2","Dll4","Gja4","Gja5","Sox17","Cxcl12"
  ),

  Venous_capillary_EC = c(
    "Nr2f2","Emcn","Plvap","Car4","Aplnr"
  ),

  Angiogenic_EC = c(
    "Esm1","Apln","Aplnr","Kdr","Pgf","Angpt2"
  ),

  Activated_inflammatory_EC = c(
    "Icam1","Vcam1","Sele","Selp","Ccl2",
    "Cxcl9","Cxcl10","Nfkbia","Socs3"
  ),

  Interferon_response = c(
    "Isg15","Ifit1","Ifit2","Ifit3","Irf7","Stat1","Cxcl10"
  ),

  Myeloid = c(
    "Ptprc","Lyz2","Adgre1","Csf1r","Tyrobp",
    "Aif1","Ctss","Fcer1g"
  ),

  HSC_mesenchymal = c(
    "Col1a1","Col1a2","Col3a1","Pdgfra","Pdgfrb",
    "Des","Lrat","Rbp1","Dcn","Lum"
  ),

  Pericyte = c(
    "Rgs5","Cspg4","Pdgfrb","Des","Acta2","Mcam"
  ),

  Cholangiocyte = c(
    "Krt19","Krt8","Krt18","Epcam","Sox9","Krt7","Muc1"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1","Apoa2","Ass1","Cps1","Hnf4a","Fabp1"
  )
)

marker_sets <- lapply(
  marker_sets,
  function(x) x[x %in% rownames(obj)]
)

cat("\n=== MARKER AVAILABILITY ===\n")

for (nm in names(marker_sets)) {
  cat(
    nm, ": ",
    paste(marker_sets[[nm]], collapse=", "),
    "\n",
    sep=""
  )
}

capture.output(
  marker_sets,
  file=file.path(TABDIR, "Expanded_marker_availability_v6.7.2.txt")
)

make_dot <- function(groups, filename, width=16) {

  genes <- unique(
    unlist(marker_sets[groups], use.names=FALSE)
  )

  genes <- genes[genes %in% rownames(obj)]

  if (length(genes) == 0) return(NULL)

  p <- DotPlot(
    obj,
    features=genes,
    assay="RNA",
    dot.scale=6
  ) +
    RotatedAxis() +
    theme_classic(base_size=10) +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank()
    )

  ggsave(
    file.path(FIGDIR, filename),
    p,
    width=width,
    height=7
  )
}

make_dot(
  c("Pan_endothelial","LSEC_identity","Capillarization"),
  "DotPlot_LSEC_identity_capillarization_v6.7.2.pdf",
  16
)

make_dot(
  c("LSEC_zonation","Portal_arterial_EC",
    "Venous_capillary_EC","Angiogenic_EC"),
  "DotPlot_endothelial_zonation_vascular_v6.7.2.pdf",
  16
)

make_dot(
  c("Activated_inflammatory_EC","Interferon_response"),
  "DotPlot_endothelial_activation_inflammation_v6.7.2.pdf",
  14
)

make_dot(
  c("Myeloid","HSC_mesenchymal","Pericyte",
    "Cholangiocyte","Hepatocyte"),
  "DotPlot_lineage_contamination_audit_v6.7.2.pdf",
  19
)

p_cluster <- DimPlot(
  obj,
  reduction="umapLSEC",
  group.by=RES_COL,
  label=TRUE,
  repel=TRUE,
  pt.size=0.45
) +
  NoLegend() +
  theme_classic(base_size=13) +
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank()
  ) +
  ggtitle("Endothelial parent - res 0.6")

ggsave(
  file.path(
    FIGDIR,
    "Endothelial_res0.6_cluster_labels_v6.7.2.pdf"
  ),
  p_cluster,
  width=9,
  height=7
)

write.csv(
  as.data.frame(
    table(
      cluster=obj@meta.data[[RES_COL]],
      sample=obj$sample
    )
  ),
  file.path(TABDIR, "Cluster_by_sample_res0.6_v6.7.2.csv"),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      cluster=obj@meta.data[[RES_COL]],
      condition=obj$condition
    )
  ),
  file.path(TABDIR, "Cluster_by_condition_res0.6_v6.7.2.csv"),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      cluster=obj@meta.data[[RES_COL]],
      original_annotation=obj$endothelial_parent_label_v670
    )
  ),
  file.path(
    TABDIR,
    "Cluster_by_original_annotation_res0.6_v6.7.2.csv"
  ),
  row.names=FALSE
)

samples <- sort(unique(as.character(obj$sample)))
clusters <- sort(unique(as.character(obj@meta.data[[RES_COL]])))

robustness <- do.call(
  rbind,
  lapply(
    clusters,
    function(cl) {

      cells <- rownames(obj@meta.data)[
        as.character(obj@meta.data[[RES_COL]]) == cl
      ]

      x <- table(
        factor(
          as.character(obj$sample[cells]),
          levels=samples
        )
      )

      f <- as.numeric(x) / sum(x)

      data.frame(
        cluster=cl,
        n_cells=sum(x),
        dominant_sample=samples[which.max(f)],
        max_sample_fraction=max(f),
        n_samples_ge_5pct=sum(f >= 0.05),
        n_samples_ge_10pct=sum(f >= 0.10),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  robustness,
  file.path(
    TABDIR,
    "Cluster_sample_robustness_res0.6_v6.7.2.csv"
  ),
  row.names=FALSE
)

all_marker_genes <- unique(
  unlist(marker_sets, use.names=FALSE)
)

avg <- AverageExpression(
  obj,
  assays="RNA",
  features=all_marker_genes,
  group.by=RES_COL,
  slot="data",
  verbose=FALSE
)$RNA

write.csv(
  avg,
  file.path(
    TABDIR,
    "AverageExpression_expanded_markers_res0.6_v6.7.2.csv"
  )
)

z <- t(scale(t(avg)))
z[!is.finite(z)] <- 0

scores <- data.frame(
  cluster=colnames(avg),
  stringsAsFactors=FALSE
)

for (nm in names(marker_sets)) {

  genes <- intersect(
    marker_sets[[nm]],
    rownames(z)
  )

  scores[[nm]] <- if (length(genes) > 0) {
    colMeans(z[genes, , drop=FALSE])
  } else {
    NA_real_
  }
}

write.csv(
  scores,
  file.path(
    TABDIR,
    "Cluster_signature_Zscores_res0.6_v6.7.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== FIND ALL POSITIVE MARKERS ===\n")

Idents(obj) <- factor(obj@meta.data[[RES_COL]])

all_markers <- FindAllMarkers(
  obj,
  assay="RNA",
  only.pos=TRUE,
  test.use="wilcox",
  min.pct=0.10,
  logfc.threshold=0.25,
  return.thresh=0.05,
  verbose=TRUE
)

write.csv(
  all_markers,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.6_all_positive_v6.7.2.csv"
  ),
  row.names=FALSE
)

fc_col <- intersect(
  c("avg_log2FC","avg_logFC"),
  colnames(all_markers)
)[1]

top20 <- do.call(
  rbind,
  lapply(
    split(all_markers, all_markers$cluster),
    function(x) {

      if (!is.na(fc_col)) {
        x <- x[
          order(-x[[fc_col]], x$p_val_adj),
          ,
          drop=FALSE
        ]
      } else {
        x <- x[
          order(x$p_val_adj),
          ,
          drop=FALSE
        ]
      }

      head(x, 20)
    }
  )
)

write.csv(
  top20,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.6_TOP20_per_cluster_v6.7.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== TOP 10 MARKERS PER CLUSTER ===\n")

for (cl in unique(as.character(top20$cluster))) {

  x <- top20[
    as.character(top20$cluster) == cl,
    ,
    drop=FALSE
  ]

  cat(
    "Cluster ", cl, ": ",
    paste(head(x$gene,10), collapse=", "),
    "\n",
    sep=""
  )
}

obj$LSEC_audit_resolution_v672 <- "res0.6"

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_endothelial_res0.6_marker_audit_v6.7.2.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.2.txt"
  )
)

cat("\n============================================\n")
cat("v6.7.2 COMPLETE\n")
cat("Cells retained:", ncol(obj), "\n")
cat("No clusters removed in v6.7.2\n")
cat("Fixed audit resolution: 0.6\n")
cat("Output:", OUTDIR, "\n")
cat("============================================\n")
