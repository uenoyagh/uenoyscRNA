suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.7.4.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.4/objects/",
  "Mouse_MASH_LSEC_clean_RPCA_resolution_scan_v6.7.4.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC sample-bias / QC audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

obj <- readRDS(INPUT_RDS)

if (!"LSECclean_res0.2" %in% colnames(obj@meta.data)) {
  stop("LSECclean_res0.2 missing")
}

if (!"LSECclean_res0.3" %in% colnames(obj@meta.data)) {
  stop("LSECclean_res0.3 missing")
}

if (!"umapLSECclean" %in% Reductions(obj)) {
  stop("umapLSECclean missing")
}

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")

cat("\n=== RES0.2 COUNTS ===\n")
print(table(obj$LSECclean_res0.2))

cat("\n=== RES0.3 COUNTS ===\n")
print(table(obj$LSECclean_res0.3))

# ---------------------------------------------------------
# Marker panels
# ---------------------------------------------------------

marker_sets <- list(

  LSEC_identity = c(
    "Clec4g","Stab1","Stab2","Lyve1",
    "Fcgr2b","Mrc1","Oit3","Dnase1l3"
  ),

  Pericentral_angiocrine = c(
    "Wnt2","Wnt9b","Rspo3","Bmp2"
  ),

  Capillarization = c(
    "Cd34","Vwf","Plvap","Emcn","Esm1"
  ),

  Inflammatory = c(
    "Icam1","Vcam1","Sele","Selp",
    "Ccl2","Cxcl9","Cxcl10"
  ),

  Interferon = c(
    "Isg15","Ifit1","Ifit2","Ifit3",
    "Irf7","Stat1"
  ),

  Angiogenic = c(
    "Apln","Aplnr","Esm1","Kdr",
    "Pgf","Angpt2"
  ),

  Immediate_early_stress = c(
    "Fos","Fosb","Jun","Junb","Jund",
    "Atf3","Egr1","Ddit3",
    "Hspa1a","Hspa1b"
  ),

  Residual_myeloid = c(
    "Ptprc","Lyz2","Adgre1","Tyrobp"
  ),

  Residual_mesenchymal = c(
    "Col1a1","Col3a1","Pdgfrb","Rgs5"
  ),

  Residual_epithelial = c(
    "Krt19","Epcam","Alb","Ttr"
  )
)

marker_sets <- lapply(
  marker_sets,
  function(x) x[x %in% rownames(obj)]
)

all_markers <- unique(
  unlist(marker_sets, use.names=FALSE)
)

# ---------------------------------------------------------
# DotPlots
# ---------------------------------------------------------

for (res in c("0.2","0.3")) {

  cluster_col <- paste0(
    "LSECclean_res", res
  )

  Idents(obj) <- factor(
    obj@meta.data[[cluster_col]]
  )

  p <- DotPlot(
    obj,
    features=all_markers,
    assay="RNA",
    dot.scale=6
  ) +
    RotatedAxis() +
    theme_classic(base_size=9) +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank()
    ) +
    ggtitle(
      paste0(
        "Clean LSEC biology/QC - resolution ",
        res
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "LSEC_biology_QC_DotPlot_res",
        res,
        "_v6.7.4.1.pdf"
      )
    ),
    p,
    width=20,
    height=6.5
  )
}

# ---------------------------------------------------------
# Cluster × sample tables
# ---------------------------------------------------------

for (res in c("0.2","0.3")) {

  cluster_col <- paste0(
    "LSECclean_res", res
  )

  tab <- table(
    cluster=obj@meta.data[[cluster_col]],
    sample=obj$sample
  )

  write.csv(
    as.data.frame(tab),
    file.path(
      TABDIR,
      paste0(
        "Cluster_by_sample_res",
        res,
        "_v6.7.4.1.csv"
      )
    ),
    row.names=FALSE
  )

  # fraction of each cluster contributed by each sample
  cluster_fraction <- prop.table(
    tab,
    margin=1
  )

  write.csv(
    as.data.frame(cluster_fraction),
    file.path(
      TABDIR,
      paste0(
        "Sample_fraction_within_cluster_res",
        res,
        "_v6.7.4.1.csv"
      )
    ),
    row.names=FALSE
  )

  # fraction of each sample assigned to each cluster
  within_sample_fraction <- prop.table(
    tab,
    margin=2
  )

  write.csv(
    as.data.frame(within_sample_fraction),
    file.path(
      TABDIR,
      paste0(
        "Cluster_fraction_within_sample_res",
        res,
        "_v6.7.4.1.csv"
      )
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Stacked composition plots
# ---------------------------------------------------------

for (res in c("0.2","0.3")) {

  cluster_col <- paste0(
    "LSECclean_res", res
  )

  df <- data.frame(
    sample=as.character(obj$sample),
    cluster=as.character(
      obj@meta.data[[cluster_col]]
    ),
    stringsAsFactors=FALSE
  )

  tab <- as.data.frame(
    table(
      sample=df$sample,
      cluster=df$cluster
    )
  )

  tab$sample_total <- ave(
    tab$Freq,
    tab$sample,
    FUN=sum
  )

  tab$fraction <- tab$Freq /
    tab$sample_total

  p <- ggplot(
    tab,
    aes(
      x=sample,
      y=fraction,
      fill=cluster
    )
  ) +
    geom_col() +
    theme_classic(base_size=12) +
    theme(
      axis.text.x=element_text(
        angle=45,
        hjust=1
      )
    ) +
    labs(
      title=paste0(
        "LSEC cluster composition by sample - res ",
        res
      ),
      x=NULL,
      y="Fraction of LSEC cells",
      fill="Cluster"
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        "LSEC_cluster_fraction_by_sample_res",
        res,
        "_v6.7.4.1.pdf"
      )
    ),
    p,
    width=9,
    height=6
  )
}

# ---------------------------------------------------------
# QC summary
# ---------------------------------------------------------

qc_candidates <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.mt_for_filter"
)

qc_cols <- intersect(
  qc_candidates,
  colnames(obj@meta.data)
)

cat(
  "\nQC columns: ",
  paste(qc_cols, collapse=", "),
  "\n",
  sep=""
)

for (res in c("0.2","0.3")) {

  cluster_col <- paste0(
    "LSECclean_res", res
  )

  tmp <- obj@meta.data

  tmp$cluster_for_audit <-
    as.character(tmp[[cluster_col]])

  grouping <- interaction(
    tmp$cluster_for_audit,
    tmp$sample,
    drop=TRUE
  )

  out <- data.frame()

  for (g in levels(grouping)) {

    idx <- grouping == g

    if (!any(idx)) next

    z <- data.frame(
      group=g,
      cluster=tmp$cluster_for_audit[
        which(idx)[1]
      ],
      sample=as.character(
        tmp$sample[
          which(idx)[1]
        ]
      ),
      n_cells=sum(idx),
      stringsAsFactors=FALSE
    )

    for (q in qc_cols) {

      v <- tmp[[q]][idx]

      z[[paste0(q,"_median")]] <-
        median(
          v,
          na.rm=TRUE
        )

      z[[paste0(q,"_IQR")]] <-
        IQR(
          v,
          na.rm=TRUE
        )
    }

    out <- rbind(out, z)
  }

  write.csv(
    out,
    file.path(
      TABDIR,
      paste0(
        "QC_by_cluster_sample_res",
        res,
        "_v6.7.4.1.csv"
      )
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Feature plot: stress / inflammation
# ---------------------------------------------------------

feature_genes <- c(
  "Clec4g",
  "Stab2",
  "Wnt2",
  "Cd34",
  "Vwf",
  "Icam1",
  "Cxcl10",
  "Stat1",
  "Fos",
  "Jun",
  "Atf3",
  "Hspa1a"
)

feature_genes <- feature_genes[
  feature_genes %in% rownames(obj)
]

p_features <- FeaturePlot(
  obj,
  features=feature_genes,
  reduction="umapLSECclean",
  order=TRUE,
  min.cutoff="q05",
  max.cutoff="q95",
  ncol=4
)

ggsave(
  file.path(
    FIGDIR,
    "LSEC_key_state_FeaturePlots_v6.7.4.1.pdf"
  ),
  p_features,
  width=16,
  height=3.4 *
    ceiling(length(feature_genes)/4)
)

# ---------------------------------------------------------
# FindAllMarkers res0.2
# ---------------------------------------------------------

cat("\n=== FIND MARKERS: RES0.2 ===\n")

Idents(obj) <- factor(
  obj$LSECclean_res0.2
)

m02 <- FindAllMarkers(
  obj,
  assay="RNA",
  only.pos=TRUE,
  test.use="wilcox",
  min.pct=0.10,
  logfc.threshold=0.25,
  return.thresh=0.05,
  verbose=TRUE
)

m02 <- m02[
  is.finite(m02$p_val_adj) &
  m02$p_val_adj < 0.05,
  ,
  drop=FALSE
]

write.csv(
  m02,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.2_SIGNIFICANT_v6.7.4.1.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# FindAllMarkers res0.3
# ---------------------------------------------------------

cat("\n=== FIND MARKERS: RES0.3 ===\n")

Idents(obj) <- factor(
  obj$LSECclean_res0.3
)

m03 <- FindAllMarkers(
  obj,
  assay="RNA",
  only.pos=TRUE,
  test.use="wilcox",
  min.pct=0.10,
  logfc.threshold=0.25,
  return.thresh=0.05,
  verbose=TRUE
)

m03 <- m03[
  is.finite(m03$p_val_adj) &
  m03$p_val_adj < 0.05,
  ,
  drop=FALSE
]

write.csv(
  m03,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.3_SIGNIFICANT_v6.7.4.1.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Top20 helper
# ---------------------------------------------------------

make_top20 <- function(x) {

  fc_col <- intersect(
    c("avg_log2FC","avg_logFC"),
    colnames(x)
  )[1]

  if (is.na(fc_col)) {
    stop("Fold-change column not found.")
  }

  y <- do.call(
    rbind,
    lapply(
      split(x, x$cluster),
      function(z) {

        z <- z[
          order(
            -z[[fc_col]],
            z$p_val_adj
          ),
          ,
          drop=FALSE
        ]

        head(z, 20)
      }
    )
  )

  rownames(y) <- NULL
  y
}

top02 <- make_top20(m02)
top03 <- make_top20(m03)

write.csv(
  top02,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.2_SIGNIFICANT_TOP20_v6.7.4.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  top03,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.3_SIGNIFICANT_TOP20_v6.7.4.1.csv"
  ),
  row.names=FALSE
)

cat("\n=== RES0.2 TOP10 ===\n")

for (cl in sort(unique(top02$cluster))) {

  z <- top02[
    top02$cluster == cl,
    ,
    drop=FALSE
  ]

  cat(
    "Cluster ", cl, ": ",
    paste(
      head(z$gene,10),
      collapse=", "
    ),
    "\n",
    sep=""
  )
}

cat("\n=== RES0.3 TOP10 ===\n")

for (cl in sort(unique(top03$cluster))) {

  z <- top03[
    top03$cluster == cl,
    ,
    drop=FALSE
  ]

  cat(
    "Cluster ", cl, ": ",
    paste(
      head(z$gene,10),
      collapse=", "
    ),
    "\n",
    sep=""
  )
}

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.4.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.4.1 COMPLETE\n")
cat("Cells:", ncol(obj), "\n")
cat("No cells removed\n")
cat("No final resolution fixed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
