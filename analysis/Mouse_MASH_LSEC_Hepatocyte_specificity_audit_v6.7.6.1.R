suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(Matrix)
  library(ggplot2)
})

VERSION <- "v6.7.6.1"

WHOLE_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

LSEC_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_v6.7.5/objects/",
  "Mouse_MASH_LSEC_annotated_v6.7.5.rds"
)

V676_TABDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_v6.7.6/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC-Hepatocyte specificity audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# ---------------------------------------------------------
# Helper: robust raw-count extraction
# ---------------------------------------------------------

get_counts <- function(obj, assay="RNA") {

  DefaultAssay(obj) <- assay

  if (inherits(obj[[assay]], "Assay5")) {

    layers <- Layers(obj[[assay]])

    count_layers <- grep(
      "^counts($|\\.)",
      layers,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop("No counts layer found in ", assay)
    }

    if (
      length(count_layers) == 1 &&
      identical(count_layers, "counts")
    ) {

      return(
        LayerData(
          obj,
          assay=assay,
          layer="counts"
        )
      )

    } else {

      obj <- JoinLayers(
        obj,
        assay=assay
      )

      return(
        LayerData(
          obj,
          assay=assay,
          layer="counts"
        )
      )
    }

  } else {

    return(
      GetAssayData(
        obj,
        assay=assay,
        layer="counts"
      )
    )
  }
}

# ---------------------------------------------------------
# Read whole-liver object
# ---------------------------------------------------------

cat("=== READ WHOLE-LIVER OBJECT ===\n")

whole <- readRDS(WHOLE_RDS)

anno_col <- "celltype_for_R8plot_FIXED2"

if (!anno_col %in% colnames(whole@meta.data)) {
  stop("Missing annotation column: ", anno_col)
}

if (!"sample" %in% colnames(whole@meta.data)) {
  stop("Missing sample column in whole object")
}

anno <- as.character(
  whole@meta.data[[anno_col]]
)

hep_labels <- unique(
  anno[
    grepl(
      "Hepat",
      anno,
      ignore.case=TRUE
    )
  ]
)

cat(
  "Hepatocyte-like labels found: ",
  paste(hep_labels, collapse=", "),
  "\n",
  sep=""
)

if ("Hepatocyte" %in% hep_labels) {

  HEP_LABEL <- "Hepatocyte"

} else {

  stop(
    "Exact label 'Hepatocyte' not found. Labels were: ",
    paste(hep_labels, collapse=", ")
  )
}

hep_cells <- colnames(whole)[
  anno == HEP_LABEL
]

cat("Hepatocyte cells:", length(hep_cells), "\n")

cat("\nHepatocytes by sample:\n")
print(
  table(
    whole$sample[hep_cells]
  )
)

# ---------------------------------------------------------
# Extract whole-liver RNA counts, then Hep cells
# ---------------------------------------------------------

cat("\n=== EXTRACT HEPATOCYTE COUNTS ===\n")

whole_counts <- get_counts(
  whole,
  assay="RNA"
)

hep_cells <- intersect(
  hep_cells,
  colnames(whole_counts)
)

hep_counts <- whole_counts[
  ,
  hep_cells,
  drop=FALSE
]

hep_md <- whole@meta.data[
  hep_cells,
  ,
  drop=FALSE
]

cat(
  "Hep matrix:",
  nrow(hep_counts),
  "genes x",
  ncol(hep_counts),
  "cells\n"
)

rm(whole_counts)
invisible(gc())

# ---------------------------------------------------------
# Read clean LSEC
# ---------------------------------------------------------

cat("\n=== READ CLEAN LSEC OBJECT ===\n")

lsec <- readRDS(LSEC_RDS)

required_lsec <- c(
  "sample",
  "condition",
  "LSEC_state_v675",
  "LSEC_QC_flag_v675"
)

missing_lsec <- setdiff(
  required_lsec,
  colnames(lsec@meta.data)
)

if (length(missing_lsec) > 0) {
  stop(
    "Missing LSEC metadata: ",
    paste(missing_lsec, collapse=", ")
  )
}

lsec_counts <- get_counts(
  lsec,
  assay="RNA"
)

cat(
  "LSEC matrix:",
  nrow(lsec_counts),
  "genes x",
  ncol(lsec_counts),
  "cells\n"
)

cat("\nLSEC by sample:\n")
print(table(lsec$sample))

# ---------------------------------------------------------
# Samples
# ---------------------------------------------------------

samples <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

samples <- samples[
  samples %in%
    intersect(
      unique(as.character(hep_md$sample)),
      unique(as.character(lsec$sample))
    )
]

cat(
  "\nSamples used:",
  paste(samples, collapse=", "),
  "\n"
)

if (length(samples) != 6) {
  stop(
    "Expected 6 matched samples; found ",
    length(samples)
  )
}

# ---------------------------------------------------------
# Pseudobulk helper
# ---------------------------------------------------------

make_pb <- function(counts, metadata, sample_col, samples) {

  pb <- sapply(
    samples,
    function(s) {

      cc <- rownames(metadata)[
        as.character(metadata[[sample_col]]) == s
      ]

      cc <- intersect(
        cc,
        colnames(counts)
      )

      if (length(cc) == 0) {
        stop("Zero cells for sample: ", s)
      }

      Matrix::rowSums(
        counts[
          ,
          cc,
          drop=FALSE
        ]
      )
    }
  )

  colnames(pb) <- samples
  pb
}

# ---------------------------------------------------------
# Hep pseudobulk
# ---------------------------------------------------------

cat("\n=== HEPATOCYTE PSEUDOBULK ===\n")

hep_pb <- make_pb(
  hep_counts,
  hep_md,
  "sample",
  samples
)

# ---------------------------------------------------------
# LSEC all pseudobulk
# ---------------------------------------------------------

cat("=== ALL CLEAN LSEC PSEUDOBULK ===\n")

lsec_md_all <- lsec@meta.data

lsec_pb_all <- make_pb(
  lsec_counts,
  lsec_md_all,
  "sample",
  samples
)

# ---------------------------------------------------------
# Primary LSEC = QC-flagged state excluded
# ---------------------------------------------------------

cat("=== PRIMARY LSEC PSEUDOBULK ===\n")

primary_cells <- rownames(lsec@meta.data)[
  lsec$LSEC_QC_flag_v675 == "Pass"
]

lsec_counts_primary <- lsec_counts[
  ,
  primary_cells,
  drop=FALSE
]

lsec_md_primary <- lsec@meta.data[
  primary_cells,
  ,
  drop=FALSE
]

cat(
  "Primary LSEC cells:",
  ncol(lsec_counts_primary),
  "\n"
)

print(
  table(
    lsec_md_primary$sample
  )
)

lsec_pb_primary <- make_pb(
  lsec_counts_primary,
  lsec_md_primary,
  "sample",
  samples
)

# ---------------------------------------------------------
# Common genes
# ---------------------------------------------------------

common_genes <- Reduce(
  intersect,
  list(
    rownames(hep_pb),
    rownames(lsec_pb_all),
    rownames(lsec_pb_primary)
  )
)

hep_pb <- hep_pb[
  common_genes,
  ,
  drop=FALSE
]

lsec_pb_all <- lsec_pb_all[
  common_genes,
  ,
  drop=FALSE
]

lsec_pb_primary <- lsec_pb_primary[
  common_genes,
  ,
  drop=FALSE
]

cat(
  "Common genes:",
  length(common_genes),
  "\n"
)

# ---------------------------------------------------------
# Convert to logCPM independently by compartment
# ---------------------------------------------------------

to_logcpm <- function(pb) {

  y <- DGEList(
    counts=pb
  )

  y <- calcNormFactors(y)

  cpm(
    y,
    log=TRUE,
    prior.count=2
  )
}

hep_logcpm <- to_logcpm(
  hep_pb
)

lsec_all_logcpm <- to_logcpm(
  lsec_pb_all
)

lsec_primary_logcpm <- to_logcpm(
  lsec_pb_primary
)

# ---------------------------------------------------------
# Gene specificity table
# ---------------------------------------------------------

specificity <- data.frame(
  gene=common_genes,
  stringsAsFactors=FALSE
)

for (s in samples) {

  specificity[[paste0(
    s,
    "_Hep_logCPM"
  )]] <- hep_logcpm[,s]

  specificity[[paste0(
    s,
    "_LSEC_all_logCPM"
  )]] <- lsec_all_logcpm[,s]

  specificity[[paste0(
    s,
    "_LSEC_primary_logCPM"
  )]] <- lsec_primary_logcpm[,s]

  specificity[[paste0(
    s,
    "_Hep_minus_LSEC_primary"
  )]] <-
    hep_logcpm[,s] -
    lsec_primary_logcpm[,s]
}

delta_cols <- paste0(
  samples,
  "_Hep_minus_LSEC_primary"
)

specificity$Hep_logCPM_median <-
  apply(
    hep_logcpm,
    1,
    median
  )

specificity$LSEC_all_logCPM_median <-
  apply(
    lsec_all_logcpm,
    1,
    median
  )

specificity$LSEC_primary_logCPM_median <-
  apply(
    lsec_primary_logcpm,
    1,
    median
  )

specificity$Hep_minus_LSEC_primary_median <-
  apply(
    specificity[
      ,
      delta_cols,
      drop=FALSE
    ],
    1,
    median
  )

specificity$n_samples_Hep_gt_LSEC_by2 <-
  rowSums(
    specificity[
      ,
      delta_cols,
      drop=FALSE
    ] > 2
  )

specificity$n_samples_LSEC_gt_Hep_by1 <-
  rowSums(
    specificity[
      ,
      delta_cols,
      drop=FALSE
    ] < -1
  )

# ---------------------------------------------------------
# Exploratory specificity classification
#
# > +2 log2 units:
# roughly >4-fold Hep dominance
#
# < -1 log2 unit:
# roughly >2-fold LSEC dominance
# ---------------------------------------------------------

specificity$Specificity_class_v6761 <-
  ifelse(
    specificity$Hep_minus_LSEC_primary_median > 2,
    "Hepatocyte_dominant",
    ifelse(
      specificity$Hep_minus_LSEC_primary_median < -1,
      "LSEC_supported",
      "Ambiguous_shared"
    )
  )

# Additional robustness flag
specificity$Specificity_robustness_v6761 <-
  ifelse(
    specificity$n_samples_Hep_gt_LSEC_by2 >= 4,
    "Hep_dominant_in_4plus_samples",
    ifelse(
      specificity$n_samples_LSEC_gt_Hep_by1 >= 4,
      "LSEC_enriched_in_4plus_samples",
      "Mixed_across_samples"
    )
  )

write.csv(
  specificity,
  file.path(
    TABDIR,
    "Hepatocyte_vs_LSEC_gene_specificity_all_genes_v6.7.6.1.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Attach specificity to v6.7.6 result files
# ---------------------------------------------------------

input_files <- c(

  "Primary_no_QC_Tx_vs_Sham_edgeR_v6.7.6.csv",
  "Primary_no_QC_replicate_consistent_effect_genes_v6.7.6.csv",

  "Shared_core_Tx_vs_Sham_edgeR_v6.7.6.csv",
  "Shared_core_replicate_consistent_effect_genes_v6.7.6.csv",

  "Inflammatory_stress_high_Tx_vs_Sham_edgeR_v6.7.6.csv",
  "Inflammatory_stress_high_replicate_consistent_effect_genes_v6.7.6.csv",

  "Homeostatic_like_Tx_vs_Sham_edgeR_v6.7.6.csv",
  "Homeostatic_like_replicate_consistent_effect_genes_v6.7.6.csv",

  "Wnt_angiocrine_high_Tx_vs_Sham_edgeR_v6.7.6.csv",
  "Wnt_angiocrine_high_replicate_consistent_effect_genes_v6.7.6.csv"
)

for (fn in input_files) {

  infile <- file.path(
    V676_TABDIR,
    fn
  )

  if (!file.exists(infile)) {
    cat(
      "SKIP missing:",
      fn,
      "\n"
    )
    next
  }

  x <- read.csv(
    infile,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  if (!"gene" %in% colnames(x)) {
    cat(
      "SKIP no gene column:",
      fn,
      "\n"
    )
    next
  }

  y <- merge(
    x,
    specificity,
    by="gene",
    all.x=TRUE,
    sort=FALSE
  )

  outfile <- sub(
    "_v6\\.7\\.6\\.csv$",
    "_with_Hep_LSEC_specificity_v6.7.6.1.csv",
    fn
  )

  write.csv(
    y,
    file.path(
      TABDIR,
      outfile
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Replicate-consistent genes combined
# ---------------------------------------------------------

candidate_files <- c(

  Primary_no_QC =
    "Primary_no_QC_replicate_consistent_effect_genes_v6.7.6.csv",

  Shared_core =
    "Shared_core_replicate_consistent_effect_genes_v6.7.6.csv",

  Inflammatory_stress_high =
    "Inflammatory_stress_high_replicate_consistent_effect_genes_v6.7.6.csv",

  Wnt_angiocrine_high =
    "Wnt_angiocrine_high_replicate_consistent_effect_genes_v6.7.6.csv"
)

candidate_list <- list()

for (nm in names(candidate_files)) {

  f <- file.path(
    V676_TABDIR,
    candidate_files[[nm]]
  )

  if (!file.exists(f)) next

  x <- read.csv(
    f,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  if (nrow(x) == 0) next

  x$analysis <- nm

  candidate_list[[nm]] <- x
}

if (length(candidate_list) > 0) {

  candidates <- do.call(
    rbind,
    candidate_list
  )

  candidates <- merge(
    candidates,
    specificity[
      ,
      c(
        "gene",
        "Hep_logCPM_median",
        "LSEC_all_logCPM_median",
        "LSEC_primary_logCPM_median",
        "Hep_minus_LSEC_primary_median",
        "n_samples_Hep_gt_LSEC_by2",
        "n_samples_LSEC_gt_Hep_by1",
        "Specificity_class_v6761",
        "Specificity_robustness_v6761"
      )
    ],
    by="gene",
    all.x=TRUE
  )

  candidates <- candidates[
    order(
      candidates$analysis,
      candidates$Specificity_class_v6761,
      candidates$FDR,
      -abs(candidates$logFC)
    ),
    ,
    drop=FALSE
  ]

  write.csv(
    candidates,
    file.path(
      TABDIR,
      "LSEC_Tx_response_replicate_consistent_candidates_specificity_audit_v6.7.6.1.csv"
    ),
    row.names=FALSE
  )

  # ---------------------------------------------
  # Clean LSEC-supported candidate subset
  # ---------------------------------------------

  supported <- candidates[
    candidates$Specificity_class_v6761 !=
      "Hepatocyte_dominant",
    ,
    drop=FALSE
  ]

  supported <- supported[
    order(
      supported$FDR,
      -abs(supported$logFC)
    ),
    ,
    drop=FALSE
  ]

  write.csv(
    supported,
    file.path(
      TABDIR,
      "LSEC_Tx_response_NON_Hepatocyte_dominant_candidates_v6.7.6.1.csv"
    ),
    row.names=FALSE
  )
}

# ---------------------------------------------------------
# Plot candidate-space: Hep vs LSEC expression
# ---------------------------------------------------------

plot_genes <- unique(
  if (exists("candidates")) {
    candidates$gene
  } else {
    character()
  }
)

plot_genes <- intersect(
  plot_genes,
  specificity$gene
)

plot_df <- specificity[
  specificity$gene %in% plot_genes,
  ,
  drop=FALSE
]

if (nrow(plot_df) > 0) {

  p <- ggplot(
    plot_df,
    aes(
      x=LSEC_primary_logCPM_median,
      y=Hep_logCPM_median,
      shape=Specificity_class_v6761
    )
  ) +
    geom_abline(
      slope=1,
      intercept=0,
      linetype=2
    ) +
    geom_abline(
      slope=1,
      intercept=2,
      linetype=3
    ) +
    geom_abline(
      slope=1,
      intercept=-1,
      linetype=3
    ) +
    geom_point(
      size=2.8
    ) +
    geom_text(
      aes(label=gene),
      size=2.6,
      nudge_x=0.08,
      check_overlap=TRUE
    ) +
    theme_classic(base_size=11) +
    labs(
      title="Tx-response candidates: Hepatocyte vs LSEC expression",
      x="Primary LSEC median logCPM",
      y="Hepatocyte median logCPM",
      shape="Specificity"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Tx_response_candidates_Hepatocyte_vs_LSEC_specificity_v6.7.6.1.pdf"
    ),
    p,
    width=9,
    height=8
  )
}

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

cat("\n=== GENE SPECIFICITY SUMMARY ===\n")

print(
  table(
    specificity$Specificity_class_v6761
  )
)

if (exists("candidates")) {

  cat("\n=== REPLICATE-CONSISTENT CANDIDATES ===\n")

  print(
    table(
      candidates$analysis,
      candidates$Specificity_class_v6761
    )
  )

  cat("\n=== NON-HEPATOCYTE-DOMINANT TOP CANDIDATES ===\n")

  tmp <- candidates[
    candidates$Specificity_class_v6761 !=
      "Hepatocyte_dominant",
    ,
    drop=FALSE
  ]

  tmp <- tmp[
    order(
      tmp$FDR,
      -abs(tmp$logFC)
    ),
    ,
    drop=FALSE
  ]

  print(
    head(
      tmp[
        ,
        c(
          "gene",
          "analysis",
          "logFC",
          "FDR",
          "Hep_minus_LSEC_primary_median",
          "Specificity_class_v6761"
        )
      ],
      30
    )
  )
}

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.6.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.6.1 COMPLETE\n")
cat("Matched samples:", length(samples), "\n")
cat("Specificity primary reference: QC-pass clean LSEC\n")
cat("No cells or genes removed from source objects\n")
cat("Thresholds are exploratory, for ambient-RNA triage\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
