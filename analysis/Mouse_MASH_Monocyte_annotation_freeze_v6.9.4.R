suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.9.4"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.3.2/objects/",
  "Mouse_MASH_Monocyte_cluster7_8_QC_audit_v6.9.3.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(
  OUTDIR,
  "tables"
)

FIGDIR <- file.path(
  OUTDIR,
  "figures"
)

OBJDIR <- file.path(
  OUTDIR,
  "objects"
)

dir.create(
  TABDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  FIGDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OBJDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

cat("====================================================\n")
cat("Mouse MASH Monocyte annotation freeze\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop(
    "Missing input RDS: ",
    INPUT_RDS
  )
}

obj <- readRDS(
  INPUT_RDS
)

cat(
  "Input cells:",
  ncol(obj),
  "\n"
)

cat(
  "Input features:",
  nrow(obj),
  "\n\n"
)

if (ncol(obj) != 3490) {
  warning(
    "Expected 3490 lineage-clean Monocytes, observed: ",
    ncol(obj)
  )
}

# =========================================================
# Frozen clustering scaffold
# =========================================================

CLUSTER_COL <- "monocyte_clean_res0_4"

if (!(CLUSTER_COL %in% colnames(obj@meta.data))) {
  stop(
    "Missing cluster column: ",
    CLUSTER_COL
  )
}

cluster_vec <- as.character(
  obj@meta.data[[CLUSTER_COL]]
)

observed_clusters <- sort(
  unique(
    cluster_vec
  )
)

expected_clusters <- as.character(
  0:9
)

if (!identical(
  observed_clusters,
  expected_clusters
)) {

  stop(
    "Unexpected res0.4 clusters.\nObserved: ",
    paste(
      observed_clusters,
      collapse=", "
    ),
    "\nExpected: ",
    paste(
      expected_clusters,
      collapse=", "
    )
  )
}

cat(
  "Frozen scaffold:",
  CLUSTER_COL,
  "\n"
)

cat(
  "Clusters:",
  paste(
    observed_clusters,
    collapse=", "
  ),
  "\n\n"
)

# =========================================================
# Frozen annotation map
# =========================================================

annotation_map <- data.frame(

  cluster=as.character(
    0:9
  ),

  fine_state=c(
    "S100a8_S100a9_Thbs1_stress_inflammatory_Monocyte",
    "Mmp8_Sell_Chil3_Vcan_classical_inflammatory_Monocyte",
    "Pald1_C3ar1_homeostatic_like_Monocyte",
    "Tnf_Il1rn_Olr1_Gpnmb_inflammatory_remodeling_Monocyte",
    "Cd300e_Pglyrp1_Cd36_S1pr5_activated_Monocyte",
    "Adamdec1_Pecam1_low_complexity_state",
    "Ms4a7_Mmp12_Dab2_C1q_monocyte_to_macrophage_transition",
    "Nos2_Cxcl9_Saa3_IFNg_inflammatory_Monocyte",
    "Hepatocyte_RNA_high_Monocyte_QC_watch",
    "Ifit_Rsad2_Cmpk2_IFN_responsive_Monocyte"
  ),

  broad_state=c(
    "Stress_inflammatory",
    "Classical_inflammatory",
    "Homeostatic_like",
    "Inflammatory_remodeling",
    "Activated_inflammatory",
    "Low_complexity_exploratory",
    "Monocyte_to_macrophage_transition",
    "IFNg_inflammatory",
    "Hepatocyte_RNA_high_QC_watch",
    "IFN_responsive"
  ),

  analysis_class=c(
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "sample_biased_exploratory",
    "primary",
    "disease_enriched_primary",
    "QC_watch_sensitivity",
    "primary"
  ),

  interpretation_note=c(
    "Stress/inflammatory Monocyte; modest S100a8/S100a9 and Thbs1 signal.",
    "Classical-inflammatory Monocyte with Mmp8, Sell, Chil3 and Vcan.",
    "Homeostatic-like Monocyte enriched for Pald1 and C3ar1.",
    "Inflammatory/remodeling Monocyte with Tnf, Il1rn, Olr1, Gpnmb, Lpl and Gpr84.",
    "Activated Monocyte enriched for Cd300e, Pglyrp1, Cd36 and S1pr5.",
    "Low-complexity Adamdec1/Pecam1-associated state; strongly Sham20-biased; exploratory.",
    "Monocyte-to-macrophage transition with Ms4a7, Mmp12, Dab2, C1qa and C1qb.",
    "Disease-enriched inflammatory/IFNg-response Monocyte with Nos2, Cxcl9, Saa3 and Ccl5.",
    "Monocyte identity retained but deep Hepatocyte RNA strongly enriched; retain as QC-watch and use exclusion sensitivity.",
    "Interferon-responsive Monocyte enriched for Ifit genes, Rsad2 and Cmpk2."
  ),

  stringsAsFactors=FALSE
)

write.csv(
  annotation_map,
  file.path(
    TABDIR,
    "Monocyte_frozen_annotation_map_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Apply frozen annotation
# =========================================================

idx <- match(
  cluster_vec,
  annotation_map$cluster
)

if (any(is.na(idx))) {
  stop(
    "Failed to map all res0.4 clusters."
  )
}

fine_state <- annotation_map$fine_state[
  idx
]

broad_state <- annotation_map$broad_state[
  idx
]

analysis_class <- annotation_map$analysis_class[
  idx
]

interpretation_note <- annotation_map$interpretation_note[
  idx
]

obj$Monocyte_cluster_frozen_v6.9.4 <- factor(
  cluster_vec,
  levels=annotation_map$cluster
)

obj$Monocyte_state_frozen_v6.9.4 <- factor(
  fine_state,
  levels=annotation_map$fine_state
)

obj$Monocyte_broad_state_v6.9.4 <- factor(
  broad_state,
  levels=annotation_map$broad_state
)

analysis_class_levels <- c(
  "primary",
  "disease_enriched_primary",
  "sample_biased_exploratory",
  "QC_watch_sensitivity"
)

obj$Monocyte_analysis_class_v6.9.4 <- factor(
  analysis_class,
  levels=analysis_class_levels
)

obj$Monocyte_interpretation_note_v6.9.4 <-
  interpretation_note

# =========================================================
# Analysis policy flags
# =========================================================

obj$Monocyte_primary_analysis_v6.9.4 <-
  analysis_class %in% c(
    "primary",
    "disease_enriched_primary"
  )

obj$Monocyte_QCwatch_exclude_sensitivity_v6.9.4 <-
  analysis_class !=
    "QC_watch_sensitivity"

obj$Monocyte_primary_core_v6.9.4 <-
  analysis_class %in% c(
    "primary",
    "disease_enriched_primary"
  )

# Explicit flags for special states

obj$Monocyte_cluster5_sample_biased_v6.9.4 <-
  cluster_vec == "5"

obj$Monocyte_cluster7_disease_enriched_v6.9.4 <-
  cluster_vec == "7"

obj$Monocyte_cluster8_QC_watch_v6.9.4 <-
  cluster_vec == "8"

# =========================================================
# Resolve sample column
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
  stop(
    "Could not resolve sample column."
  )
}

sample_vec <- as.character(
  obj@meta.data[[sample_col]]
)

sample_names <- sort(
  unique(
    sample_vec
  )
)

cat(
  "Sample column:",
  sample_col,
  "\n"
)

cat(
  "Samples:",
  paste(
    sample_names,
    collapse=", "
  ),
  "\n\n"
)

# =========================================================
# Condition assignment
# =========================================================

condition <- ifelse(
  grepl(
    "^STD",
    sample_vec,
    ignore.case=TRUE
  ),
  "STD",
  ifelse(
    grepl(
      "CDHFD|CDAHFD",
      sample_vec,
      ignore.case=TRUE
    ),
    "CDHFD",
    ifelse(
      grepl(
        "^Sham",
        sample_vec,
        ignore.case=TRUE
      ),
      "Sham",
      ifelse(
        grepl(
          "^Tx",
          sample_vec,
          ignore.case=TRUE
        ),
        "Tx",
        "Other"
      )
    )
  )
)

obj$condition_v6.9.4 <- factor(
  condition,
  levels=c(
    "STD",
    "CDHFD",
    "Sham",
    "Tx",
    "Other"
  )
)

# =========================================================
# Fine-state counts
# =========================================================

fine_count <- as.data.frame(
  table(
    state=
      obj$Monocyte_state_frozen_v6.9.4
  ),
  stringsAsFactors=FALSE
)

fine_count <- fine_count[
  fine_count$Freq > 0,
  ,
  drop=FALSE
]

write.csv(
  fine_count,
  file.path(
    TABDIR,
    "Monocyte_frozen_state_counts_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Broad-state counts
# =========================================================

broad_count <- as.data.frame(
  table(
    broad_state=
      obj$Monocyte_broad_state_v6.9.4
  ),
  stringsAsFactors=FALSE
)

broad_count <- broad_count[
  broad_count$Freq > 0,
  ,
  drop=FALSE
]

write.csv(
  broad_count,
  file.path(
    TABDIR,
    "Monocyte_frozen_broad_state_counts_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Analysis-class counts
# =========================================================

analysis_class_count <- as.data.frame(
  table(
    analysis_class=
      obj$Monocyte_analysis_class_v6.9.4
  ),
  stringsAsFactors=FALSE
)

analysis_class_count <-
  analysis_class_count[
    analysis_class_count$Freq > 0,
    ,
    drop=FALSE
  ]

write.csv(
  analysis_class_count,
  file.path(
    TABDIR,
    "Monocyte_analysis_class_counts_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# State x sample counts and fractions
# =========================================================

fine_levels <- annotation_map$fine_state

state_sample_rows <- list()
counter <- 1

for (state in fine_levels) {

  for (s in sample_names) {

    sample_cells <- rownames(
      obj@meta.data
    )[
      sample_vec == s
    ]

    state_cells <- sample_cells[
      as.character(
        obj@meta.data[
          sample_cells,
          "Monocyte_state_frozen_v6.9.4"
        ]
      ) == state
    ]

    n_sample <- length(
      sample_cells
    )

    n_state <- length(
      state_cells
    )

    state_sample_rows[[counter]] <-
      data.frame(
        state=state,
        sample=s,
        n_cells=n_state,
        sample_total=n_sample,
        fraction=
          if (
            n_sample > 0
          ) {
            n_state / n_sample
          } else {
            NA_real_
          },
        stringsAsFactors=FALSE
      )

    counter <- counter + 1
  }
}

state_sample <- do.call(
  rbind,
  state_sample_rows
)

write.csv(
  state_sample,
  file.path(
    TABDIR,
    "Monocyte_frozen_state_by_sample_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Broad state x sample
# =========================================================

broad_levels <- annotation_map$broad_state

broad_sample_rows <- list()
counter <- 1

for (state in broad_levels) {

  for (s in sample_names) {

    sample_cells <- rownames(
      obj@meta.data
    )[
      sample_vec == s
    ]

    state_cells <- sample_cells[
      as.character(
        obj@meta.data[
          sample_cells,
          "Monocyte_broad_state_v6.9.4"
        ]
      ) == state
    ]

    n_sample <- length(
      sample_cells
    )

    n_state <- length(
      state_cells
    )

    broad_sample_rows[[counter]] <-
      data.frame(
        broad_state=state,
        sample=s,
        n_cells=n_state,
        sample_total=n_sample,
        fraction=
          if (
            n_sample > 0
          ) {
            n_state / n_sample
          } else {
            NA_real_
          },
        stringsAsFactors=FALSE
      )

    counter <- counter + 1
  }
}

broad_sample <- do.call(
  rbind,
  broad_sample_rows
)

write.csv(
  broad_sample,
  file.path(
    TABDIR,
    "Monocyte_frozen_broad_state_by_sample_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# State x condition counts
#
# STD/CDHFD remain descriptive because n=1/group.
# =========================================================

condition_levels <- c(
  "STD",
  "CDHFD",
  "Sham",
  "Tx"
)

state_condition_rows <- list()
counter <- 1

for (state in fine_levels) {

  for (cond in condition_levels) {

    cond_cells <- rownames(
      obj@meta.data
    )[
      as.character(
        obj$condition_v6.9.4
      ) == cond
    ]

    state_cells <- cond_cells[
      as.character(
        obj@meta.data[
          cond_cells,
          "Monocyte_state_frozen_v6.9.4"
        ]
      ) == state
    ]

    state_condition_rows[[counter]] <-
      data.frame(
        state=state,
        condition=cond,
        n_cells=
          length(
            state_cells
          ),
        condition_total=
          length(
            cond_cells
          ),
        fraction=
          if (
            length(
              cond_cells
            ) > 0
          ) {
            length(
              state_cells
            ) /
              length(
                cond_cells
              )
          } else {
            NA_real_
          },
        stringsAsFactors=FALSE
      )

    counter <- counter + 1
  }
}

state_condition <- do.call(
  rbind,
  state_condition_rows
)

write.csv(
  state_condition,
  file.path(
    TABDIR,
    "Monocyte_frozen_state_by_condition_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Sample robustness / sample bias
# =========================================================

robustness_rows <- list()

for (state in fine_levels) {

  x <- state_sample[
    state_sample$state == state,
    ,
    drop=FALSE
  ]

  total_state_cells <- sum(
    x$n_cells
  )

  n_samples_present <- sum(
    x$n_cells > 0
  )

  max_sample_n <- max(
    x$n_cells
  )

  max_sample_fraction_of_state <-
    if (
      total_state_cells > 0
    ) {
      max_sample_n /
        total_state_cells
    } else {
      NA_real_
    }

  dominant_sample <-
    if (
      total_state_cells > 0
    ) {
      x$sample[
        which.max(
          x$n_cells
        )
      ]
    } else {
      NA_character_
    }

  robustness_rows[[state]] <-
    data.frame(
      state=state,
      total_cells=
        total_state_cells,
      n_samples_present=
        n_samples_present,
      dominant_sample=
        dominant_sample,
      dominant_sample_fraction=
        max_sample_fraction_of_state,
      sample_biased_ge60pct=
        !is.na(
          max_sample_fraction_of_state
        ) &&
        max_sample_fraction_of_state >= 0.60,
      sample_biased_ge80pct=
        !is.na(
          max_sample_fraction_of_state
        ) &&
        max_sample_fraction_of_state >= 0.80,
      stringsAsFactors=FALSE
    )
}

robustness <- do.call(
  rbind,
  robustness_rows
)

rownames(robustness) <- NULL

write.csv(
  robustness,
  file.path(
    TABDIR,
    "Monocyte_frozen_state_sample_robustness_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cluster-level QC summary
# =========================================================

qc_cols <- c(
  "deep_hep_hits_v6.9.3.2",
  "deep_hep_umi_per10k_v6.9.3.2",
  "monocyte_core_hits_v6.9.3.2",
  "monocyte_core_umi_per10k_v6.9.3.2",
  "endothelial_core_hits_v6.9.3.2"
)

missing_qc <- setdiff(
  qc_cols,
  colnames(
    obj@meta.data
  )
)

if (length(missing_qc) > 0) {
  warning(
    "Missing some v6.9.3.2 QC columns: ",
    paste(
      missing_qc,
      collapse=", "
    )
  )
}

qc_summary_rows <- list()

for (cl in annotation_map$cluster) {

  cells <- rownames(
    obj@meta.data
  )[
    cluster_vec == cl
  ]

  md <- obj@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  row <- data.frame(
    cluster=cl,
    state=
      annotation_map$fine_state[
        annotation_map$cluster == cl
      ],
    n_cells=
      length(
        cells
      ),
    median_nCount_RNA=
      median(
        md$nCount_RNA
      ),
    median_nFeature_RNA=
      median(
        md$nFeature_RNA
      ),
    stringsAsFactors=FALSE
  )

  if (
    "deep_hep_hits_v6.9.3.2" %in%
      colnames(md)
  ) {

    row$deep_hep_median_hits <-
      median(
        md$deep_hep_hits_v6.9.3.2
      )

    row$deep_hep_fraction_ge3 <-
      mean(
        md$deep_hep_hits_v6.9.3.2 >= 3
      )
  }

  if (
    "monocyte_core_hits_v6.9.3.2" %in%
      colnames(md)
  ) {

    row$monocyte_median_hits <-
      median(
        md$monocyte_core_hits_v6.9.3.2
      )

    row$monocyte_fraction_ge2 <-
      mean(
        md$monocyte_core_hits_v6.9.3.2 >= 2
      )
  }

  if (
    "endothelial_core_hits_v6.9.3.2" %in%
      colnames(md)
  ) {

    row$endothelial_fraction_ge2 <-
      mean(
        md$endothelial_core_hits_v6.9.3.2 >= 2
      )
  }

  qc_summary_rows[[cl]] <- row
}

qc_summary <- do.call(
  rbind,
  qc_summary_rows
)

rownames(qc_summary) <- NULL

write.csv(
  qc_summary,
  file.path(
    TABDIR,
    "Monocyte_frozen_cluster_QC_summary_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# UMAP figures
# =========================================================

fine_palette <- c(
  "#E41A1C",
  "#FF7F00",
  "#4DAF4A",
  "#A65628",
  "#F781BF",
  "#999999",
  "#377EB8",
  "#984EA3",
  "#BDBDBD",
  "#00BFC4"
)

names(
  fine_palette
) <- annotation_map$fine_state

broad_palette <- c(
  "#E41A1C",
  "#FF7F00",
  "#4DAF4A",
  "#A65628",
  "#F781BF",
  "#999999",
  "#377EB8",
  "#984EA3",
  "#BDBDBD",
  "#00BFC4"
)

names(
  broad_palette
) <- annotation_map$broad_state

p_fine <- DimPlot(
  obj,
  reduction="umap",
  group.by=
    "Monocyte_state_frozen_v6.9.4",
  label=TRUE,
  repel=TRUE,
  pt.size=0.65
) +
  scale_color_manual(
    values=fine_palette
  ) +
  ggtitle(
    "Mouse MASH Monocyte frozen annotation v6.9.4"
  ) +
  theme_classic(
    base_size=10
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_frozen_annotation_UMAP_v6.9.4.pdf"
  ),
  p_fine,
  width=12,
  height=8
)

p_broad <- DimPlot(
  obj,
  reduction="umap",
  group.by=
    "Monocyte_broad_state_v6.9.4",
  label=TRUE,
  repel=TRUE,
  pt.size=0.65
) +
  scale_color_manual(
    values=broad_palette
  ) +
  ggtitle(
    "Mouse MASH Monocyte broad states v6.9.4"
  ) +
  theme_classic(
    base_size=10
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_frozen_broad_state_UMAP_v6.9.4.pdf"
  ),
  p_broad,
  width=11,
  height=8
)

p_class <- DimPlot(
  obj,
  reduction="umap",
  group.by=
    "Monocyte_analysis_class_v6.9.4",
  pt.size=0.65
) +
  ggtitle(
    "Mouse MASH Monocyte analysis classes v6.9.4"
  ) +
  theme_classic(
    base_size=10
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_analysis_class_UMAP_v6.9.4.pdf"
  ),
  p_class,
  width=9,
  height=7
)

p_sample <- DimPlot(
  obj,
  reduction="umap",
  group.by=sample_col,
  pt.size=0.60
) +
  ggtitle(
    "Mouse MASH Monocyte frozen object - by sample"
  ) +
  theme_classic(
    base_size=10
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_frozen_UMAP_by_sample_v6.9.4.pdf"
  ),
  p_sample,
  width=9,
  height=7
)

# =========================================================
# Frozen annotation metadata summary
# =========================================================

policy <- data.frame(

  item=c(
    "Frozen resolution",
    "Lineage-clean cells",
    "Primary analysis",
    "Cluster 5 policy",
    "Cluster 7 policy",
    "Cluster 8 policy",
    "STD/CDHFD inference",
    "Sham/Tx inference"
  ),

  decision=c(
    "RPCA res0.4",
    "3490",
    "Primary + disease_enriched_primary states",
    "Retain as sample_biased_exploratory; do not use as standalone reproducible treatment state",
    "Retain as disease_enriched_primary; high complexity does not support bulk exclusion",
    "Retain as QC_watch_sensitivity; primary transcriptional analyses require exclusion sensitivity",
    "Descriptive only; n=1 biological sample/group",
    "Biological-sample comparison with n=2/group; prioritize effect size and replicate concordance"
  ),

  stringsAsFactors=FALSE
)

write.csv(
  policy,
  file.path(
    TABDIR,
    "Monocyte_annotation_analysis_policy_v6.9.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Save frozen object
# =========================================================

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_annotation_frozen_v6.9.4.rds"
  )
)

# =========================================================
# Human-readable summary
# =========================================================

summary_lines <- c(
  "# Mouse MASH Monocyte annotation freeze v6.9.4",
  "",
  "## Frozen baseline",
  "- Lineage-clean Monocyte cells: 3,490.",
  "- Frozen clustering scaffold: RPCA res0.4.",
  "- Frozen fine states: 10.",
  "",
  "## Lineage cleanup",
  "- v6.9.2 removed only res0.4 clusters 7 and 8 from the pre-cleanup object.",
  "- Those removed clusters represented T/NK-like and B-cell contamination.",
  "- No neutrophil-like, macrophage-transition, inflammatory, or IFN-responsive Monocyte states were removed.",
  "",
  "## Frozen fine states",
  "- Cluster 0: S100a8/S100a9/Thbs1 stress-inflammatory Monocyte.",
  "- Cluster 1: Mmp8/Sell/Chil3/Vcan classical-inflammatory Monocyte.",
  "- Cluster 2: Pald1/C3ar1 homeostatic-like Monocyte.",
  "- Cluster 3: Tnf/Il1rn/Olr1/Gpnmb inflammatory-remodeling Monocyte.",
  "- Cluster 4: Cd300e/Pglyrp1/Cd36/S1pr5 activated Monocyte.",
  "- Cluster 5: Adamdec1/Pecam1 low-complexity state; sample-biased exploratory.",
  "- Cluster 6: Ms4a7/Mmp12/Dab2/C1q monocyte-to-macrophage transition.",
  "- Cluster 7: Nos2/Cxcl9/Saa3 IFNg-inflammatory Monocyte; disease-enriched primary.",
  "- Cluster 8: Hepatocyte-RNA-high Monocyte QC-watch.",
  "- Cluster 9: Ifit/Rsad2/Cmpk2 IFN-responsive Monocyte.",
  "",
  "## Special-state policy",
  "- Cluster 5 is retained but treated as sample-biased exploratory.",
  "- Cluster 7 is retained as a disease-enriched inflammatory state.",
  "- Cluster 8 is retained as QC-watch because Monocyte identity is preserved despite reproducible deep-Hepatocyte RNA enrichment.",
  "- Primary transcriptional analyses should include a cluster-8 exclusion sensitivity analysis.",
  "",
  "## Statistical limits",
  "- STD vs CDHFD: n=1/group, descriptive only.",
  "- Sham vs Tx: n=2/group.",
  "- Treatment interpretation should emphasize effect size and replicate concordance rather than p values alone."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Mouse_MASH_Monocyte_annotation_freeze_summary_v6.9.4.md"
  )
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== FROZEN ANNOTATION MAP ===\n")
print(
  annotation_map[
    ,
    c(
      "cluster",
      "fine_state",
      "broad_state",
      "analysis_class"
    )
  ],
  row.names=FALSE
)

cat("\n=== FROZEN STATE COUNTS ===\n")
print(
  fine_count,
  row.names=FALSE
)

cat("\n=== ANALYSIS CLASS COUNTS ===\n")
print(
  analysis_class_count,
  row.names=FALSE
)

cat("\n=== STATE SAMPLE ROBUSTNESS ===\n")
print(
  robustness,
  row.names=FALSE
)

cat("\n=== CLUSTER QC SUMMARY ===\n")
print(
  qc_summary,
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.4 COMPLETE\n")
cat("Monocyte annotation freeze complete\n")
cat("Frozen resolution: res0.4\n")
cat("Frozen states: 10\n")
cat("Cells:", ncol(obj), "\n")
cat("No cells removed\n")
cat("Cluster 5: sample-biased exploratory\n")
cat("Cluster 7: disease-enriched primary\n")
cat("Cluster 8: QC-watch sensitivity\n")
cat(
  "Frozen object:",
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_annotation_frozen_v6.9.4.rds"
  ),
  "\n"
)
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.4.txt"
  )
)
