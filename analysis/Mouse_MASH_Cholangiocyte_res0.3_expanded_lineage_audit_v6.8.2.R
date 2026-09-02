suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260902)

VERSION <- "v6.8.2"
RES <- 0.3
CLUSTER_COL <- "integrated_snn_res.0.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.1/objects/",
  "Mouse_MASH_Cholangiocyte_RPCA_resolution_scan_v6.8.1.rds"
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
cat("Mouse MASH Cholangiocyte expanded lineage audit\n")
cat("Version:", VERSION, "\n")
cat("Diagnostic resolution:", RES, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

cat("=== INPUT OBJECT ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")
cat("Reductions:", paste(Reductions(obj), collapse=", "), "\n")

if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("Missing cluster column: ", CLUSTER_COL)
}

required_md <- c(
  "sample",
  "condition",
  "nCount_RNA",
  "nFeature_RNA"
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

obj$Chol_res03_v682 <- factor(
  obj@meta.data[[CLUSTER_COL]],
  levels=sort(
    unique(
      as.numeric(
        as.character(
          obj@meta.data[[CLUSTER_COL]]
        )
      )
    )
  )
)

Idents(obj) <- "Chol_res03_v682"

cat("\n=== RES0.3 CLUSTER COUNTS ===\n")
print(
  sort(
    table(obj$Chol_res03_v682),
    decreasing=TRUE
  )
)

# =========================================================
# Helper: raw RNA counts
# =========================================================

get_counts <- function(x, assay="RNA") {

  DefaultAssay(x) <- assay

  if (inherits(x[[assay]], "Assay5")) {

    layers <- Layers(x[[assay]])

    count_layers <- grep(
      "^counts($|\\.)",
      layers,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop("No counts layer in assay ", assay)
    }

    if (
      length(count_layers) > 1 ||
      !identical(count_layers, "counts")
    ) {

      x <- JoinLayers(
        x,
        assay=assay
      )
    }

    return(
      LayerData(
        x,
        assay=assay,
        layer="counts"
      )
    )
  }

  GetAssayData(
    x,
    assay=assay,
    layer="counts"
  )
}

# =========================================================
# Marker panels
# =========================================================

marker_sets <- list(

  Biliary_core = c(
    "Krt19","Krt7","Epcam","Sox9",
    "Muc1","Hnf1b","Cftr","Slc4a2",
    "Krt8","Krt18"
  ),

  Ductular_reactive = c(
    "Spp1","Mmp7","Krt23",
    "Tacstd2","Prom1","Klf5",
    "Epcam","Sox9"
  ),

  Inflammatory_adhesion = c(
    "Icam1","Vcam1",
    "Cxcl1","Cxcl2","Ccl2",
    "Il6","Nfkbia","Socs3"
  ),

  Fibrogenic_remodeling = c(
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
    "Mki67","Top2a","Birc5",
    "Ube2c","Cenpf","Pcna","Stmn1"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1",
    "Fabp1","Cps1","Ass1",
    "Hnf4a","Asgr1","Tat"
  ),

  Myeloid = c(
    "Ptprc","Lyz2","Tyrobp",
    "Aif1","Adgre1","Csf1r",
    "Fcgr3"
  ),

  Vascular_endothelial = c(
    "Pecam1","Cdh5","Kdr",
    "Esam","Emcn","Klf2"
  ),

  LSEC = c(
    "Stab1","Stab2","Lyve1",
    "Clec4g","Fcgr2b","Kdr"
  ),

  HSC_mesenchymal = c(
    "Col1a1","Col3a1",
    "Pdgfra","Pdgfrb",
    "Lrat","Rbp1",
    "Des","Rgs5"
  ),

  Neutrophil = c(
    "S100a8","S100a9",
    "Retnlg","Mmp8",
    "Ly6g","Csf3r"
  ),

  Lymphoid = c(
    "Cd3d","Cd3e","Trbc1",
    "Nkg7","Klrb1c",
    "Ms4a1","Cd79a"
  ),

  RBC = c(
    "Hba-a1","Hba-a2",
    "Hbb-bs","Hbb-bt"
  ),

  Platelet = c(
    "Pf4","Ppbp","Gp9"
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

  cat("\n[", nm, "]\n", sep="")

  cat(
    "Present: ",
    paste(
      marker_sets_present[[nm]],
      collapse=", "
    ),
    "\n",
    sep=""
  )

  missing_genes <- setdiff(
    marker_sets[[nm]],
    rownames(obj)
  )

  cat(
    "Missing: ",
    paste(
      missing_genes,
      collapse=", "
    ),
    "\n",
    sep=""
  )
}

# =========================================================
# Cluster x sample composition
# =========================================================

cluster_sample_counts <- as.data.frame(
  table(
    cluster=obj$Chol_res03_v682,
    sample=obj$sample
  )
)

write.csv(
  cluster_sample_counts,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_cluster_by_sample_counts_v6.8.2.csv"
  ),
  row.names=FALSE
)

cluster_sample_matrix <- table(
  obj$Chol_res03_v682,
  obj$sample
)

cluster_fraction_within_cluster <-
  prop.table(
    cluster_sample_matrix,
    margin=1
  )

write.csv(
  as.data.frame.matrix(
    cluster_fraction_within_cluster
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_sample_fraction_within_cluster_v6.8.2.csv"
  )
)

sample_fraction_within_sample <-
  prop.table(
    cluster_sample_matrix,
    margin=2
  )

write.csv(
  as.data.frame.matrix(
    sample_fraction_within_sample
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_cluster_fraction_within_sample_v6.8.2.csv"
  )
)

# =========================================================
# Sample-dominance summary
# =========================================================

robustness <- do.call(
  rbind,
  lapply(
    rownames(cluster_sample_matrix),
    function(cl) {

      x <- cluster_sample_matrix[
        cl,
        ,
        drop=TRUE
      ]

      total <- sum(x)

      dominant_sample <- names(
        which.max(x)
      )

      dominant_n <- max(x)

      samples_ge_5pct <- sum(
        x / total >= 0.05
      )

      samples_ge_10pct <- sum(
        x / total >= 0.10
      )

      data.frame(
        cluster=cl,
        n_cells=total,
        dominant_sample=dominant_sample,
        dominant_n=dominant_n,
        dominant_fraction=
          dominant_n / total,
        samples_ge_5pct=
          samples_ge_5pct,
        samples_ge_10pct=
          samples_ge_10pct,
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(robustness) <- NULL

robustness$sample_specific_flag <-
  robustness$dominant_fraction >= 0.80

write.csv(
  robustness,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_cluster_sample_robustness_v6.8.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== CLUSTER SAMPLE ROBUSTNESS ===\n")
print(robustness)

# =========================================================
# QC by cluster and cluster x sample
# =========================================================

qc_cluster <- do.call(
  rbind,
  lapply(
    levels(obj$Chol_res03_v682),
    function(cl) {

      md <- obj@meta.data[
        obj$Chol_res03_v682 == cl,
        ,
        drop=FALSE
      ]

      data.frame(
        cluster=cl,
        n_cells=nrow(md),
        median_nCount_RNA=
          median(md$nCount_RNA),
        median_nFeature_RNA=
          median(md$nFeature_RNA),
        mean_nCount_RNA=
          mean(md$nCount_RNA),
        mean_nFeature_RNA=
          mean(md$nFeature_RNA),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  qc_cluster,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_QC_by_cluster_v6.8.2.csv"
  ),
  row.names=FALSE
)

qc_cluster_sample <- do.call(
  rbind,
  lapply(
    split(
      obj@meta.data,
      interaction(
        obj$Chol_res03_v682,
        obj$sample,
        drop=TRUE
      )
    ),
    function(md) {

      if (nrow(md) == 0) {
        return(NULL)
      }

      data.frame(
        cluster=
          as.character(
            md$Chol_res03_v682[1]
          ),
        sample=
          as.character(
            md$sample[1]
          ),
        n_cells=nrow(md),
        median_nCount_RNA=
          median(md$nCount_RNA),
        median_nFeature_RNA=
          median(md$nFeature_RNA),
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(qc_cluster_sample) <- NULL

write.csv(
  qc_cluster_sample,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_QC_by_cluster_sample_v6.8.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Raw-count marker hit/fraction audit
# =========================================================

cat("\n=== RAW-COUNT LINEAGE FRACTION AUDIT ===\n")

counts <- get_counts(
  obj,
  assay="RNA"
)

counts <- counts[
  ,
  colnames(obj),
  drop=FALSE
]

total_umi <- Matrix::colSums(counts)

calc_fraction <- function(genes) {

  genes <- intersect(
    genes,
    rownames(counts)
  )

  if (length(genes) == 0) {
    return(
      rep(
        0,
        ncol(counts)
      )
    )
  }

  Matrix::colSums(
    counts[
      genes,
      ,
      drop=FALSE
    ]
  ) /
    pmax(total_umi, 1)
}

calc_hits <- function(genes) {

  genes <- intersect(
    genes,
    rownames(counts)
  )

  if (length(genes) == 0) {
    return(
      rep(
        0,
        ncol(counts)
      )
    )
  }

  Matrix::colSums(
    counts[
      genes,
      ,
      drop=FALSE
    ] > 0
  )
}

audit_sets <- c(
  "Biliary_core",
  "Hepatocyte",
  "Myeloid",
  "Vascular_endothelial",
  "LSEC",
  "HSC_mesenchymal",
  "Neutrophil",
  "Lymphoid",
  "Cycling"
)

for (nm in audit_sets) {

  obj@meta.data[[paste0(nm, "_count_fraction_v682")]] <- calc_fraction(
    marker_sets_present[[nm]]
  )

  obj@meta.data[[paste0(nm, "_hits_v682")]] <- calc_hits(
    marker_sets_present[[nm]]
  )
}

# =========================================================
# Cluster-level lineage audit
# =========================================================

lineage_cluster_summary <- do.call(
  rbind,
  lapply(
    levels(obj$Chol_res03_v682),
    function(cl) {

      idx <- obj$Chol_res03_v682 == cl

      d <- data.frame(
        cluster=cl,
        n_cells=sum(idx),
        stringsAsFactors=FALSE
      )

      for (nm in audit_sets) {

        frac_col <-
          paste0(
            nm,
            "_count_fraction_v682"
          )

        hit_col <-
          paste0(
            nm,
            "_hits_v682"
          )

        d[[paste0(nm, "_fraction_median")]] <- median(
          obj@meta.data[
            idx,
            frac_col
          ]
        )

        d[[paste0(nm, "_hits_median")]] <- median(
          obj@meta.data[
            idx,
            hit_col
          ]
        )

        d[[paste0(nm, "_pct_ge2hits")]] <- mean(
          obj@meta.data[
            idx,
            hit_col
          ] >= 2
        ) * 100
      }

      d
    }
  )
)

rownames(lineage_cluster_summary) <- NULL

write.csv(
  lineage_cluster_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_lineage_fraction_hit_summary_v6.8.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cluster x sample lineage audit
# =========================================================

cluster_sample_groups <- split(
  seq_len(ncol(obj)),
  interaction(
    obj$Chol_res03_v682,
    obj$sample,
    drop=TRUE
  )
)

lineage_cluster_sample_summary <- do.call(
  rbind,
  lapply(
    cluster_sample_groups,
    function(idx) {

      md <- obj@meta.data[
        idx,
        ,
        drop=FALSE
      ]

      d <- data.frame(
        cluster=
          as.character(
            md$Chol_res03_v682[1]
          ),
        sample=
          as.character(
            md$sample[1]
          ),
        n_cells=nrow(md),
        stringsAsFactors=FALSE
      )

      for (nm in audit_sets) {

        frac_col <-
          paste0(
            nm,
            "_count_fraction_v682"
          )

        hit_col <-
          paste0(
            nm,
            "_hits_v682"
          )

        d[[paste0(nm, "_fraction_median")]] <- median(
          md[[frac_col]]
        )

        d[[paste0(nm, "_pct_ge2hits")]] <- mean(
          md[[hit_col]] >= 2
        ) * 100
      }

      d
    }
  )
)

rownames(
  lineage_cluster_sample_summary
) <- NULL

write.csv(
  lineage_cluster_sample_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_lineage_audit_by_cluster_sample_v6.8.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Expanded lineage DotPlot
# =========================================================

lineage_features <- unique(
  unlist(
    marker_sets_present[
      c(
        "Biliary_core",
        "Hepatocyte",
        "Myeloid",
        "Vascular_endothelial",
        "LSEC",
        "HSC_mesenchymal",
        "Neutrophil",
        "Lymphoid",
        "RBC",
        "Platelet"
      )
    ],
    use.names=FALSE
  )
)

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Chol_res03_v682"

p_lineage <- DotPlot(
  obj,
  features=lineage_features,
  assay="RNA",
  dot.scale=6
) +
  RotatedAxis() +
  theme_classic(base_size=8) +
  theme(
    axis.title=element_blank()
  ) +
  ggtitle(
    "Cholangiocyte res0.3 expanded lineage audit"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.3_expanded_lineage_DotPlot_v6.8.2.pdf"
  ),
  p_lineage,
  width=22,
  height=7
)

# =========================================================
# Cholangiocyte biology DotPlot
# =========================================================

biology_features <- unique(
  unlist(
    marker_sets_present[
      c(
        "Biliary_core",
        "Ductular_reactive",
        "Inflammatory_adhesion",
        "Fibrogenic_remodeling",
        "Stress_IEG",
        "Cycling"
      )
    ],
    use.names=FALSE
  )
)

p_biology <- DotPlot(
  obj,
  features=biology_features,
  assay="RNA",
  dot.scale=6
) +
  RotatedAxis() +
  theme_classic(base_size=8) +
  theme(
    axis.title=element_blank()
  ) +
  ggtitle(
    "Cholangiocyte res0.3 biological-state audit"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.3_biology_DotPlot_v6.8.2.pdf"
  ),
  p_biology,
  width=22,
  height=7
)

# =========================================================
# UMAP with res0.3 labels
# =========================================================

p_umap <- DimPlot(
  obj,
  reduction="umap",
  group.by="Chol_res03_v682",
  label=TRUE,
  repel=TRUE,
  pt.size=0.25
) +
  ggtitle(
    "Cholangiocyte res0.3 diagnostic scaffold"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.3_diagnostic_UMAP_v6.8.2.pdf"
  ),
  p_umap,
  width=8.5,
  height=7.5
)

# =========================================================
# UMAPs of lineage-count fractions
# =========================================================

fraction_features <- paste0(
  c(
    "Biliary_core",
    "Hepatocyte",
    "Myeloid",
    "Vascular_endothelial",
    "LSEC",
    "HSC_mesenchymal",
    "Neutrophil"
  ),
  "_count_fraction_v682"
)

plots <- FeaturePlot(
  obj,
  features=fraction_features,
  reduction="umap",
  cols=c(
    "grey90",
    "#0033FF",
    "#FF1A1A"
  ),
  ncol=3,
  pt.size=0.2
)

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.3_lineage_count_fraction_UMAP_v6.8.2.pdf"
  ),
  plots,
  width=15,
  height=12
)

# =========================================================
# Sample-composition bar plot
# =========================================================

comp <- as.data.frame(
  cluster_sample_matrix
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

      x$fraction <-
        x$n / sum(x$n)

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
  geom_col(
    width=0.8
  ) +
  theme_classic(base_size=11) +
  labs(
    title=
      "Cholangiocyte res0.3: sample composition within cluster",
    x="Cluster",
    y="Fraction within cluster",
    fill="Sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_res0.3_sample_composition_v6.8.2.pdf"
  ),
  p_comp,
  width=10,
  height=6
)

# =========================================================
# FindAllMarkers
# =========================================================

cat("\n=== FIND ALL MARKERS: res0.3 ===\n")

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Chol_res03_v682"

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
    "FindAllMarkers_res0.3_all_positive_v6.8.2.csv"
  ),
  row.names=FALSE
)

fc_col <- intersect(
  c(
    "avg_log2FC",
    "avg_logFC"
  ),
  colnames(markers)
)[1]

if (is.na(fc_col)) {
  stop("No fold-change column found.")
}

sig <- markers[
  is.finite(markers$p_val_adj) &
    markers$p_val_adj < 0.05,
  ,
  drop=FALSE
]

sig <- sig[
  order(
    as.numeric(
      as.character(
        sig$cluster
      )
    ),
    -sig[[fc_col]],
    sig$p_val_adj
  ),
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
    "FindAllMarkers_res0.3_SIGNIFICANT_TOP20_v6.8.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Print marker summary
# =========================================================

cat("\n=== TOP 10 SIGNIFICANT MARKERS PER CLUSTER ===\n")

for (
  cl in sort(
    unique(
      as.numeric(
        as.character(
          top20$cluster
        )
      )
    )
  )
) {

  x <- top20[
    as.character(top20$cluster) ==
      as.character(cl),
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

# =========================================================
# Preliminary automatic flags
#
# Diagnostic only.
# These are NOT deletion calls.
# =========================================================

audit <- lineage_cluster_summary

get_pct <- function(name) {
  audit[[paste0(name, "_pct_ge2hits")]]
}

audit$strong_biliary_flag <-
  get_pct("Biliary_core") >= 25

audit$myeloid_flag <-
  get_pct("Myeloid") >= 25

audit$vascular_flag <-
  get_pct("Vascular_endothelial") >= 25

audit$lsec_flag <-
  get_pct("LSEC") >= 25

audit$mesenchymal_flag <-
  get_pct("HSC_mesenchymal") >= 25

audit$neutrophil_flag <-
  get_pct("Neutrophil") >= 25

audit$lymphoid_flag <-
  get_pct("Lymphoid") >= 25

audit$hepatocyte_high_flag <-
  get_pct("Hepatocyte") >= 50

audit$preliminary_lineage_flag <- "Review"

audit$preliminary_lineage_flag[
  audit$strong_biliary_flag &
    !audit$myeloid_flag &
    !audit$vascular_flag &
    !audit$lsec_flag &
    !audit$mesenchymal_flag &
    !audit$neutrophil_flag &
    !audit$lymphoid_flag
] <- "Likely_biliary"

audit$preliminary_lineage_flag[
  audit$myeloid_flag
] <- "Possible_myeloid"

audit$preliminary_lineage_flag[
  audit$vascular_flag |
    audit$lsec_flag
] <- "Possible_endothelial"

audit$preliminary_lineage_flag[
  audit$mesenchymal_flag
] <- "Possible_mesenchymal"

audit$preliminary_lineage_flag[
  audit$neutrophil_flag
] <- "Possible_neutrophil"

audit$preliminary_lineage_flag[
  audit$lymphoid_flag
] <- "Possible_lymphoid"

write.csv(
  audit,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_preliminary_lineage_flags_v6.8.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== PRELIMINARY LINEAGE FLAGS ===\n")
print(
  audit[
    ,
    c(
      "cluster",
      "n_cells",
      "Biliary_core_pct_ge2hits",
      "Hepatocyte_pct_ge2hits",
      "Myeloid_pct_ge2hits",
      "Vascular_endothelial_pct_ge2hits",
      "LSEC_pct_ge2hits",
      "HSC_mesenchymal_pct_ge2hits",
      "Neutrophil_pct_ge2hits",
      "Lymphoid_pct_ge2hits",
      "preliminary_lineage_flag"
    )
  ]
)

# =========================================================
# Save audit object
# =========================================================

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_res0.3_audit_v6.8.2.rds"
  ),
  compress=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte res0.3 expanded lineage audit v6.8.2",
  "",
  paste0("- Input cells: ", ncol(obj)),
  "- Diagnostic clustering scaffold: RPCA res0.3",
  "- No cells removed",
  "- No final Cholangiocyte state annotation assigned",
  "- Cycling rescue candidates from v6.8.0.1 were not added",
  "- Hepatocyte RNA alone is not used as an exclusion criterion",
  "- Multiple coherent competing-lineage markers are required for removal decisions",
  "- Preliminary lineage flags are diagnostic only"
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_res0.3_expanded_lineage_audit_summary_v6.8.2.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.2.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.2 COMPLETE\n")
cat("Input cells:", ncol(obj), "\n")
cat("Diagnostic resolution: res0.3\n")
cat("No cells removed\n")
cat("No final annotation assigned\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
