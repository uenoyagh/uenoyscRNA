suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

VERSION <- "v6.9.3.2"

MON_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.3/objects/",
  "Mouse_MASH_Monocyte_clean_RPCA_resolution_scan_v6.9.3.rds"
)

WHOLE_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

ANNOTATION_COL <- "celltype_for_R8plot_FIXED2"
CLUSTER_COL <- "monocyte_clean_res0_4"

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
cat("Mouse MASH Monocyte cluster7/8 QC audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(MON_RDS)) {
  stop("Missing Monocyte RDS: ", MON_RDS)
}

if (!file.exists(WHOLE_RDS)) {
  stop("Missing whole-cell RDS: ", WHOLE_RDS)
}

mon <- readRDS(MON_RDS)

if (!(CLUSTER_COL %in% colnames(mon@meta.data))) {
  stop("Missing cluster column: ", CLUSTER_COL)
}

DefaultAssay(mon) <- "RNA"

# =========================================================
# Helpers
# =========================================================

get_counts <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  assay <- obj[["RNA"]]

  if (inherits(assay, "Assay5")) {

    count_layers <- grep(
      "^counts",
      Layers(assay),
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

  GetAssayData(
    obj,
    assay="RNA",
    layer="counts"
  )
}

resolve_sample_col <- function(obj) {

  candidates <- c(
    "sample",
    "sample_id",
    "Sample",
    "orig.ident"
  )

  out <- candidates[
    candidates %in% colnames(obj@meta.data)
  ][1]

  if (is.na(out)) {
    stop("Could not resolve sample column.")
  }

  out
}

# =========================================================
# Monocyte counts
# =========================================================

counts <- get_counts(mon)

sample_col <- resolve_sample_col(mon)

sample_vec <- as.character(
  mon@meta.data[[sample_col]]
)

cluster_vec <- as.character(
  mon@meta.data[[CLUSTER_COL]]
)

# =========================================================
# Deep Hepatocyte markers
#
# Avoid highly abundant secreted genes such as Alb/Ttr/Apoa1.
# Focus on hepatocyte metabolic identity.
# =========================================================

deep_hep <- c(
  "Fabp1",
  "Aldob",
  "Hao1",
  "Bhmt",
  "Tdo2",
  "Hal",
  "Pck1",
  "Hsd17b13",
  "Cyp2e1",
  "Cyp3a11"
)

monocyte_core <- c(
  "Lyz2",
  "Ccr2",
  "Fcgr3",
  "Lst1",
  "Tyrobp",
  "Ctss"
)

endothelial_core <- c(
  "Pecam1",
  "Cdh5",
  "Kdr",
  "Emcn",
  "Esam",
  "Ramp2"
)

deep_hep <- intersect(
  deep_hep,
  rownames(counts)
)

monocyte_core <- intersect(
  monocyte_core,
  rownames(counts)
)

endothelial_core <- intersect(
  endothelial_core,
  rownames(counts)
)

cat(
  "Deep Hepatocyte markers present:",
  length(deep_hep),
  "\n"
)

cat(
  paste(deep_hep, collapse=", "),
  "\n\n"
)

# =========================================================
# Cell-level burden
# =========================================================

deep_hep_hits <- Matrix::colSums(
  counts[
    deep_hep,
    ,
    drop=FALSE
  ] > 0
)

deep_hep_umi <- Matrix::colSums(
  counts[
    deep_hep,
    ,
    drop=FALSE
  ]
)

monocyte_hits <- Matrix::colSums(
  counts[
    monocyte_core,
    ,
    drop=FALSE
  ] > 0
)

monocyte_umi <- Matrix::colSums(
  counts[
    monocyte_core,
    ,
    drop=FALSE
  ]
)

endothelial_hits <- Matrix::colSums(
  counts[
    endothelial_core,
    ,
    drop=FALSE
  ] > 0
)

ncount <- mon@meta.data[
  colnames(mon),
  "nCount_RNA"
]

deep_hep_per10k <-
  10000 * deep_hep_umi /
  pmax(ncount, 1)

monocyte_per10k <-
  10000 * monocyte_umi /
  pmax(ncount, 1)

mon$deep_hep_hits_v6.9.3.2 <-
  deep_hep_hits[colnames(mon)]

mon$deep_hep_umi_per10k_v6.9.3.2 <-
  deep_hep_per10k[colnames(mon)]

mon$monocyte_core_hits_v6.9.3.2 <-
  monocyte_hits[colnames(mon)]

mon$monocyte_core_umi_per10k_v6.9.3.2 <-
  monocyte_per10k[colnames(mon)]

mon$endothelial_core_hits_v6.9.3.2 <-
  endothelial_hits[colnames(mon)]

# =========================================================
# Cluster-level burden
# =========================================================

cluster_summary <- do.call(
  rbind,
  lapply(
    sort(unique(cluster_vec)),
    function(cl) {

      cells <- colnames(mon)[
        cluster_vec == cl
      ]

      md <- mon@meta.data[
        cells,
        ,
        drop=FALSE
      ]

      data.frame(
        cluster=cl,
        n_cells=length(cells),

        median_nCount_RNA=
          median(md$nCount_RNA),

        median_nFeature_RNA=
          median(md$nFeature_RNA),

        deep_hep_median_hits=
          median(
            md$deep_hep_hits_v6.9.3.2
          ),

        deep_hep_fraction_ge3=
          mean(
            md$deep_hep_hits_v6.9.3.2 >= 3
          ),

        deep_hep_median_UMI_per10k=
          median(
            md$deep_hep_umi_per10k_v6.9.3.2
          ),

        monocyte_median_hits=
          median(
            md$monocyte_core_hits_v6.9.3.2
          ),

        monocyte_fraction_ge2=
          mean(
            md$monocyte_core_hits_v6.9.3.2 >= 2
          ),

        monocyte_median_UMI_per10k=
          median(
            md$monocyte_core_umi_per10k_v6.9.3.2
          ),

        endothelial_fraction_ge2=
          mean(
            md$endothelial_core_hits_v6.9.3.2 >= 2
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  cluster_summary,
  file.path(
    TABDIR,
    "Monocyte_cluster_deepHep_burden_v6.9.3.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cluster 8 vs other Monocytes by sample
# =========================================================

cluster8_sample <- do.call(
  rbind,
  lapply(
    sort(unique(sample_vec)),
    function(s) {

      cells_s <- colnames(mon)[
        sample_vec == s
      ]

      c8 <- cells_s[
        cluster_vec[
          match(cells_s, colnames(mon))
        ] == "8"
      ]

      other <- setdiff(
        cells_s,
        c8
      )

      data.frame(
        sample=s,

        cluster8_n=length(c8),
        other_n=length(other),

        cluster8_deepHep_hits_median=
          if (length(c8) > 0)
            median(
              mon@meta.data[
                c8,
                "deep_hep_hits_v6.9.3.2"
              ]
            )
          else NA_real_,

        other_deepHep_hits_median=
          if (length(other) > 0)
            median(
              mon@meta.data[
                other,
                "deep_hep_hits_v6.9.3.2"
              ]
            )
          else NA_real_,

        cluster8_deepHep_per10k_median=
          if (length(c8) > 0)
            median(
              mon@meta.data[
                c8,
                "deep_hep_umi_per10k_v6.9.3.2"
              ]
            )
          else NA_real_,

        other_deepHep_per10k_median=
          if (length(other) > 0)
            median(
              mon@meta.data[
                other,
                "deep_hep_umi_per10k_v6.9.3.2"
              ]
            )
          else NA_real_,

        cluster8_monocyte_hits_median=
          if (length(c8) > 0)
            median(
              mon@meta.data[
                c8,
                "monocyte_core_hits_v6.9.3.2"
              ]
            )
          else NA_real_,

        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  cluster8_sample,
  file.path(
    TABDIR,
    "Monocyte_cluster8_vs_other_by_sample_v6.9.3.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cluster 7 high-complexity audit
#
# Compare cluster 7 with same-sample non-cluster7 Monocytes.
# Extreme threshold = Q3 + 3 IQR.
# =========================================================

cluster7_qc <- do.call(
  rbind,
  lapply(
    sort(unique(sample_vec)),
    function(s) {

      cells_s <- colnames(mon)[
        sample_vec == s
      ]

      cl7 <- cells_s[
        cluster_vec[
          match(cells_s, colnames(mon))
        ] == "7"
      ]

      ref <- setdiff(
        cells_s,
        cl7
      )

      if (length(ref) == 0) {
        return(NULL)
      }

      nc_ref <- mon@meta.data[
        ref,
        "nCount_RNA"
      ]

      nf_ref <- mon@meta.data[
        ref,
        "nFeature_RNA"
      ]

      nc_cut <-
        as.numeric(
          quantile(nc_ref, 0.75)
        ) +
        3 * IQR(nc_ref)

      nf_cut <-
        as.numeric(
          quantile(nf_ref, 0.75)
        ) +
        3 * IQR(nf_ref)

      if (length(cl7) > 0) {

        nc7 <- mon@meta.data[
          cl7,
          "nCount_RNA"
        ]

        nf7 <- mon@meta.data[
          cl7,
          "nFeature_RNA"
        ]

        extreme_count <- nc7 > nc_cut
        extreme_feature <- nf7 > nf_cut

        data.frame(
          sample=s,
          cluster7_n=length(cl7),

          cluster7_median_nCount=
            median(nc7),

          reference_median_nCount=
            median(nc_ref),

          nCount_ratio=
            median(nc7) /
            median(nc_ref),

          cluster7_median_nFeature=
            median(nf7),

          reference_median_nFeature=
            median(nf_ref),

          nFeature_ratio=
            median(nf7) /
            median(nf_ref),

          extreme_nCount_fraction=
            mean(extreme_count),

          extreme_nFeature_fraction=
            mean(extreme_feature),

          extreme_both_fraction=
            mean(
              extreme_count &
              extreme_feature
            ),

          stringsAsFactors=FALSE
        )

      } else {

        data.frame(
          sample=s,
          cluster7_n=0,
          cluster7_median_nCount=NA_real_,
          reference_median_nCount=median(nc_ref),
          nCount_ratio=NA_real_,
          cluster7_median_nFeature=NA_real_,
          reference_median_nFeature=median(nf_ref),
          nFeature_ratio=NA_real_,
          extreme_nCount_fraction=NA_real_,
          extreme_nFeature_fraction=NA_real_,
          extreme_both_fraction=NA_real_,
          stringsAsFactors=FALSE
        )
      }
    }
  )
)

write.csv(
  cluster7_qc,
  file.path(
    TABDIR,
    "Monocyte_cluster7_high_complexity_by_sample_v6.9.3.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Whole-cell Hepatocyte reference
# =========================================================

whole <- readRDS(WHOLE_RDS)

if (!(ANNOTATION_COL %in% colnames(whole@meta.data))) {
  stop("Missing whole-cell annotation column.")
}

hep_cells <- rownames(whole@meta.data)[
  as.character(
    whole@meta.data[[ANNOTATION_COL]]
  ) == "Hepatocyte"
]

cat(
  "\nWhole-cell Hepatocyte reference cells:",
  length(hep_cells),
  "\n"
)

if (length(hep_cells) == 0) {
  stop("No Hepatocyte reference cells found.")
}

hep <- subset(
  whole,
  cells=hep_cells
)

hep_counts <- get_counts(hep)

hep_sample_col <- resolve_sample_col(hep)

hep_sample <- as.character(
  hep@meta.data[[hep_sample_col]]
)

# =========================================================
# Pseudobulk CPM for deep Hepatocyte markers
# =========================================================

aggregate_logcpm <- function(
  mat,
  samples,
  genes
) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  sample_names <- sort(
    unique(samples)
  )

  out <- matrix(
    NA_real_,
    nrow=length(genes),
    ncol=length(sample_names),
    dimnames=list(
      genes,
      sample_names
    )
  )

  for (s in sample_names) {

    cells <- colnames(mat)[
      samples == s
    ]

    pb <- Matrix::rowSums(
      mat[
        ,
        cells,
        drop=FALSE
      ]
    )

    lib <- sum(pb)

    cpm <-
      1e6 * pb /
      max(lib, 1)

    out[
      genes,
      s
    ] <- log2(
      cpm[genes] + 0.5
    )
  }

  out
}

hep_lcpm <- aggregate_logcpm(
  hep_counts,
  hep_sample,
  deep_hep
)

mon_cluster8_cells <- colnames(mon)[
  cluster_vec == "8"
]

mon_other_cells <- colnames(mon)[
  cluster_vec != "8"
]

cluster8_lcpm <- aggregate_logcpm(
  counts[
    ,
    mon_cluster8_cells,
    drop=FALSE
  ],
  sample_vec[
    match(
      mon_cluster8_cells,
      colnames(mon)
    )
  ],
  deep_hep
)

other_lcpm <- aggregate_logcpm(
  counts[
    ,
    mon_other_cells,
    drop=FALSE
  ],
  sample_vec[
    match(
      mon_other_cells,
      colnames(mon)
    )
  ],
  deep_hep
)

reference_summary <- data.frame(
  gene=deep_hep,

  Hepatocyte_median_logCPM=
    apply(
      hep_lcpm,
      1,
      median,
      na.rm=TRUE
    ),

  Cluster8_median_logCPM=
    apply(
      cluster8_lcpm,
      1,
      median,
      na.rm=TRUE
    ),

  OtherMonocyte_median_logCPM=
    apply(
      other_lcpm,
      1,
      median,
      na.rm=TRUE
    ),

  stringsAsFactors=FALSE
)

reference_summary$Cluster8_minus_OtherMonocyte <-
  reference_summary$Cluster8_median_logCPM -
  reference_summary$OtherMonocyte_median_logCPM

reference_summary$Hepatocyte_minus_Cluster8 <-
  reference_summary$Hepatocyte_median_logCPM -
  reference_summary$Cluster8_median_logCPM

write.csv(
  reference_summary,
  file.path(
    TABDIR,
    "Monocyte_cluster8_deepHep_reference_comparison_v6.9.3.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Diagnostic figures
# =========================================================

plot_df <- data.frame(
  cluster=factor(
    cluster_vec,
    levels=sort(unique(cluster_vec))
  ),
  deep_hep_per10k=
    mon@meta.data$deep_hep_umi_per10k_v6.9.3.2,
  nCount_RNA=
    mon@meta.data$nCount_RNA
)

p1 <- ggplot(
  plot_df,
  aes(
    x=cluster,
    y=deep_hep_per10k
  )
) +
  geom_violin(
    trim=TRUE
  ) +
  geom_boxplot(
    width=0.15,
    outlier.shape=NA
  ) +
  theme_classic() +
  labs(
    title="Monocyte deep-Hepatocyte RNA burden",
    x="res0.4 cluster",
    y="Deep-Hepatocyte marker UMI / 10k RNA"
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_deepHep_burden_by_cluster_v6.9.3.2.pdf"
  ),
  p1,
  width=9,
  height=6
)

p2 <- ggplot(
  plot_df,
  aes(
    x=cluster,
    y=nCount_RNA
  )
) +
  geom_violin(
    trim=TRUE
  ) +
  geom_boxplot(
    width=0.15,
    outlier.shape=NA
  ) +
  scale_y_log10() +
  theme_classic() +
  labs(
    title="Monocyte RNA complexity by cluster",
    x="res0.4 cluster",
    y="nCount_RNA (log10 scale)"
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_nCount_by_cluster_v6.9.3.2.pdf"
  ),
  p2,
  width=9,
  height=6
)

# =========================================================
# Save diagnostic object
# =========================================================

saveRDS(
  mon,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_cluster7_8_QC_audit_v6.9.3.2.rds"
  )
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== DEEP HEPATOCYTE BURDEN BY CLUSTER ===\n")
print(
  cluster_summary,
  row.names=FALSE
)

cat("\n=== CLUSTER 8 VS OTHER MONOCYTES BY SAMPLE ===\n")
print(
  cluster8_sample,
  row.names=FALSE
)

cat("\n=== CLUSTER 7 HIGH COMPLEXITY BY SAMPLE ===\n")
print(
  cluster7_qc,
  row.names=FALSE
)

cat("\n=== CLUSTER 8 DEEP HEPATOCYTE REFERENCE COMPARISON ===\n")
print(
  reference_summary,
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.3.2 COMPLETE\n")
cat("Cluster7 high-complexity audit complete\n")
cat("Cluster8 deep-Hepatocyte reference audit complete\n")
cat("No cells removed\n")
cat("No annotation changed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.3.2.txt"
  )
)
