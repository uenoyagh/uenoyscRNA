suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.5"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.4/objects/",
  "Mouse_MASH_Cholangiocyte_clean_RPCA_resolution_scan_v6.8.4.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

RES03 <- "integrated_snn_res.0.3"
RES04 <- "integrated_snn_res.0.4"

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte res0.4 state audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

cat("=== INPUT ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")

required_md <- c(
  "sample",
  "condition",
  RES03,
  RES04,
  "precleanup_QC_watch_v684"
)

missing_md <- setdiff(
  required_md,
  colnames(obj@meta.data)
)

if (length(missing_md) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing_md, collapse=", ")
  )
}

obj$Chol_res04_v685 <- factor(
  as.character(obj@meta.data[[RES04]]),
  levels=sort(
    unique(
      as.numeric(
        as.character(obj@meta.data[[RES04]])
      )
    )
  )
)

Idents(obj) <- "Chol_res04_v685"

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

obj$sample <- factor(
  as.character(obj$sample),
  levels=sample_order
)

# ---------------------------------------------------------
# Cluster counts
# ---------------------------------------------------------

cat("\n=== RES0.4 CLUSTER COUNTS ===\n")

cluster_counts <- table(
  obj$Chol_res04_v685
)

print(cluster_counts)

write.csv(
  data.frame(
    cluster=names(cluster_counts),
    n_cells=as.integer(cluster_counts)
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_cluster_counts_v6.8.5.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# res0.3 -> res0.4 mapping
# ---------------------------------------------------------

map_03_04 <- table(
  res0.3=obj@meta.data[[RES03]],
  res0.4=obj@meta.data[[RES04]]
)

cat("\n=== RES0.3 -> RES0.4 MAPPING ===\n")
print(map_03_04)

write.csv(
  as.data.frame.matrix(map_03_04),
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_to_res0.4_mapping_v6.8.5.csv"
  )
)

# ---------------------------------------------------------
# Cluster x sample
# ---------------------------------------------------------

cluster_sample <- table(
  cluster=obj$Chol_res04_v685,
  sample=obj$sample
)

cat("\n=== RES0.4 CLUSTER x SAMPLE ===\n")
print(cluster_sample)

write.csv(
  as.data.frame.matrix(cluster_sample),
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_cluster_by_sample_v6.8.5.csv"
  )
)

# ---------------------------------------------------------
# Sample robustness
# ---------------------------------------------------------

robustness <- do.call(
  rbind,
  lapply(
    rownames(cluster_sample),
    function(cl) {

      x <- cluster_sample[
        cl,
        ,
        drop=TRUE
      ]

      total <- sum(x)

      dominant_sample <- names(
        which.max(x)
      )

      dominant_n <- max(x)

      data.frame(
        cluster=cl,
        n_cells=total,
        dominant_sample=dominant_sample,
        dominant_n=dominant_n,
        dominant_fraction=dominant_n / total,
        samples_present=sum(x > 0),
        samples_ge_5pct=sum(x / total >= 0.05),
        samples_ge_10pct=sum(x / total >= 0.10),
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(robustness) <- NULL

robustness$sample_specific_flag <-
  robustness$dominant_fraction >= 0.80

cat("\n=== SAMPLE ROBUSTNESS ===\n")
print(robustness)

write.csv(
  robustness,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_sample_robustness_v6.8.5.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Precleanup QC-watch enrichment
# ---------------------------------------------------------

qcwatch <- as.character(
  obj$precleanup_QC_watch_v684
) == "Precleanup_cluster4_QC_watch"

qc_summary <- do.call(
  rbind,
  lapply(
    levels(obj$Chol_res04_v685),
    function(cl) {

      idx <- obj$Chol_res04_v685 == cl

      data.frame(
        cluster=cl,
        n_cells=sum(idx),
        QC_watch_n=sum(qcwatch[idx]),
        QC_watch_fraction=mean(qcwatch[idx]),
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(qc_summary) <- NULL

cat("\n=== PRECLEANUP CLUSTER4 QC-WATCH ENRICHMENT ===\n")
print(
  qc_summary[
    order(
      -qc_summary$QC_watch_fraction
    ),
    ,
    drop=FALSE
  ]
)

write.csv(
  qc_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_QCwatch_enrichment_v6.8.5.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Expanded biological marker panels
# ---------------------------------------------------------

marker_sets <- list(

  Biliary_identity = c(
    "Krt19","Krt7","Krt8","Krt18",
    "Epcam","Sox9","Muc1",
    "Hnf1b","Cftr","Slc4a2"
  ),

  Ductular_reactive = c(
    "Spp1","Mmp7","Krt23",
    "Tacstd2","Prom1","Klf5",
    "Sox9","Epcam"
  ),

  Inflammatory_adhesion = c(
    "Icam1","Vcam1",
    "Cxcl1","Cxcl2","Ccl2",
    "Il6","Nfkbia","Socs3"
  ),

  Remodeling_fibrogenic = c(
    "Tgfb1","Tgfb2",
    "Ccn2","Jag1",
    "Serpine1","Thbs1",
    "Inhba"
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

  Krt20_Cdh17_epithelial = c(
    "Krt20","Cdh17",
    "Inhba","Fst"
  ),

  Dmbt1_Duox2_epithelial = c(
    "Dmbt1","Duox2","Duoxa2",
    "Slc26a9","Gcnt3",
    "Ern2","Tns4"
  ),

  Ciliated = c(
    "Cfap73","Cfap44",
    "Lrrc23","Ttc21a",
    "Drc1","Iqub",
    "Dnali1","Mns1"
  ),

  Tuft_like = c(
    "Pou2f3","Trpm5",
    "Gnat3"
  ),

  Hepatocyte_ambient = c(
    "Alb","Ttr","Apoa1",
    "Fabp1","Cps1","Hnf4a"
  ),

  Myeloid = c(
    "Ptprc","Lyz2",
    "Tyrobp","Aif1",
    "Adgre1","Csf1r"
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

marker_sets_present <- lapply(
  marker_sets,
  function(x) {
    intersect(x, rownames(obj))
  }
)

cat("\n=== MARKER AVAILABILITY ===\n")

for (nm in names(marker_sets)) {

  missing_genes <- setdiff(
    marker_sets[[nm]],
    rownames(obj)
  )

  cat(
    nm,
    ": present=",
    length(marker_sets_present[[nm]]),
    "; missing=",
    paste(missing_genes, collapse=", "),
    "\n",
    sep=""
  )
}

marker_features <- unique(
  unlist(
    marker_sets_present,
    use.names=FALSE
  )
)

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Chol_res04_v685"

# ---------------------------------------------------------
# Expanded DotPlot
# ---------------------------------------------------------

p_dot <- DotPlot(
  obj,
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
    "Clean Cholangiocyte res0.4 expanded state audit"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.4_expanded_state_DotPlot_v6.8.5.pdf"
  ),
  p_dot,
  width=26,
  height=7
)

# ---------------------------------------------------------
# Average expression
# ---------------------------------------------------------

avg <- AverageExpression(
  obj,
  assays="RNA",
  features=marker_features,
  group.by="Chol_res04_v685",
  layer="data",
  verbose=FALSE
)$RNA

write.csv(
  avg,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_marker_AverageExpression_v6.8.5.csv"
  )
)

# ---------------------------------------------------------
# FindAllMarkers
# ---------------------------------------------------------

cat("\n=== FIND ALL MARKERS ===\n")

markers <- FindAllMarkers(
  obj,
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
    "Cholangiocyte_res0.4_FindAllMarkers_all_positive_v6.8.5.csv"
  ),
  row.names=FALSE
)

fc_col <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(markers)
)[1]

if (is.na(fc_col)) {
  stop("No fold-change column detected.")
}

sig <- markers[
  is.finite(markers$p_val_adj) &
    markers$p_val_adj < 0.05,
  ,
  drop=FALSE
]

top20 <- do.call(
  rbind,
  lapply(
    split(sig, sig$cluster),
    function(x) {

      x <- x[
        order(
          -x[[fc_col]],
          x$p_val_adj
        ),
        ,
        drop=FALSE
      ]

      head(x, 20)
    }
  )
)

rownames(top20) <- NULL

write.csv(
  top20,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.4_FindAllMarkers_TOP20_v6.8.5.csv"
  ),
  row.names=FALSE
)

cat("\n=== TOP 10 MARKERS PER CLUSTER ===\n")

for (cl in levels(obj$Chol_res04_v685)) {

  x <- top20[
    as.character(top20$cluster) == cl,
    ,
    drop=FALSE
  ]

  cat(
    "Cluster ",
    cl,
    ": ",
    paste(
      head(x$gene, 10),
      collapse=", "
    ),
    "\n",
    sep=""
  )
}

# ---------------------------------------------------------
# Sample composition bar plot
# ---------------------------------------------------------

comp <- as.data.frame(
  cluster_sample
)

colnames(comp) <- c(
  "cluster",
  "sample",
  "n"
)

comp <- do.call(
  rbind,
  lapply(
    split(comp, comp$cluster),
    function(x) {
      x$fraction <- x$n / sum(x$n)
      x
    }
  )
)

rownames(comp) <- NULL

p_comp <- ggplot(
  comp,
  aes(
    x=cluster,
    y=fraction,
    fill=sample
  )
) +
  geom_col() +
  theme_classic(base_size=11) +
  labs(
    title=
      "Clean Cholangiocyte res0.4: sample composition",
    x="res0.4 cluster",
    y="Fraction within cluster",
    fill="Sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.4_sample_composition_v6.8.5.pdf"
  ),
  p_comp,
  width=10,
  height=6
)

# ---------------------------------------------------------
# Save audit object
# ---------------------------------------------------------

saveRDS(
  obj,
  file.path(
    OUTDIR,
    "Mouse_MASH_Cholangiocyte_res0.4_state_audit_v6.8.5.rds"
  ),
  compress=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte res0.4 state audit v6.8.5",
  "",
  paste0("- Cells: ", ncol(obj)),
  "- Clean baseline: v6.8.3.2",
  "- RPCA object: v6.8.4",
  "- Annotation scaffold selected: res0.4",
  "- No final biological state names assigned yet.",
  "- Precleanup cluster 4 remains explicitly tracked as QC-watch.",
  "- STD/CDHFD disease-axis interpretation remains descriptive because n=1 per group."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_res0.4_state_audit_summary_v6.8.5.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.5.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.5 COMPLETE\n")
cat("Cells:", ncol(obj), "\n")
cat("Annotation scaffold: res0.4\n")
cat("No final state names assigned yet\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
