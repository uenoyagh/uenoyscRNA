suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.8.8.1"

WHOLE_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

CHOL_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.7/objects/",
  "Mouse_MASH_Cholangiocyte_state_module_scored_v6.8.7.rds"
)

DE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.8/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")

dir.create(
  TABDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte reference specificity audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(WHOLE_RDS)) {
  stop("Missing whole-cell RDS.")
}

if (!file.exists(CHOL_RDS)) {
  stop("Missing clean Cholangiocyte RDS.")
}

whole <- readRDS(WHOLE_RDS)
chol <- readRDS(CHOL_RDS)

ANNOTATION_COL <- "celltype_for_R8plot_FIXED2"

required_md <- c(
  ANNOTATION_COL,
  "sample"
)

missing_md <- setdiff(
  required_md,
  colnames(whole@meta.data)
)

if (length(missing_md) > 0) {
  stop(
    "Missing whole-cell metadata: ",
    paste(missing_md, collapse=", ")
  )
}

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

cat("=== INPUT OBJECTS ===\n")
cat("Whole cells:", ncol(whole), "\n")
cat("Clean Cholangiocyte cells:", ncol(chol), "\n")

# =========================================================
# Raw-count extractor
# =========================================================

get_raw_counts <- function(
  obj,
  assay="RNA"
) {

  DefaultAssay(obj) <- assay

  if (inherits(obj[[assay]], "Assay5")) {

    layers <- Layers(
      obj[[assay]]
    )

    count_layers <- grep(
      "^counts($|\\.)",
      layers,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop("No RNA count layer found.")
    }

    if (
      length(count_layers) > 1 ||
      !identical(count_layers, "counts")
    ) {

      obj <- JoinLayers(
        obj,
        assay=assay
      )
    }

    return(
      LayerData(
        obj,
        assay=assay,
        layer="counts"
      )
    )
  }

  GetAssayData(
    obj,
    assay=assay,
    layer="counts"
  )
}

cat("\n=== EXTRACT WHOLE-CELL RAW COUNTS ===\n")

counts <- get_raw_counts(
  whole,
  assay="RNA"
)

cat(
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)

# =========================================================
# Confirm clean Cholangiocytes exist in whole-cell matrix
# =========================================================

chol_cells <- intersect(
  colnames(chol),
  colnames(counts)
)

if (length(chol_cells) != ncol(chol)) {
  stop(
    "Not all clean Cholangiocyte cells found in whole-cell object."
  )
}

# =========================================================
# Reference lineages
# =========================================================

annotation <- as.character(
  whole@meta.data[[ANNOTATION_COL]]
)

reference_labels <- c(
  Hepatocyte="Hepatocyte",
  Kupffer_Macrophage="Kupffer_Macrophage",
  LSEC="LSEC",
  HSC_Mesenchymal="HSC_Mesenchymal",
  Monocyte="Monocyte",
  Neutrophil="Neutrophil"
)

reference_cells <- list(
  Cholangiocyte_clean=chol_cells
)

for (nm in names(reference_labels)) {

  label <- reference_labels[[nm]]

  reference_cells[[nm]] <-
    rownames(whole@meta.data)[
      annotation == label
    ]
}

cat("\n=== REFERENCE CELL COUNTS ===\n")

for (nm in names(reference_cells)) {

  cat(
    nm,
    ": ",
    length(reference_cells[[nm]]),
    "\n",
    sep=""
  )
}

# =========================================================
# Cell counts by reference lineage x sample
# =========================================================

sample_count_rows <- list()

for (nm in names(reference_cells)) {

  cells <- reference_cells[[nm]]

  md <- whole@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  n_by_sample <- table(
    factor(
      as.character(md$sample),
      levels=sample_order
    )
  )

  sample_count_rows[[nm]] <- data.frame(
    lineage=nm,
    sample=sample_order,
    n_cells=as.integer(n_by_sample),
    stringsAsFactors=FALSE
  )
}

sample_count_df <- do.call(
  rbind,
  sample_count_rows
)

rownames(sample_count_df) <- NULL

write.csv(
  sample_count_df,
  file.path(
    TABDIR,
    "Reference_lineage_cells_by_sample_v6.8.8.1.csv"
  ),
  row.names=FALSE
)

cat("\n=== REFERENCE CELLS BY SAMPLE ===\n")
print(
  sample_count_df,
  row.names=FALSE
)

# =========================================================
# Raw pseudobulk per lineage/sample
# =========================================================

aggregate_reference <- function(cells) {

  cells <- intersect(
    cells,
    colnames(counts)
  )

  out <- sapply(
    sample_order,
    function(s) {

      cells_s <- cells[
        as.character(
          whole@meta.data[
            cells,
            "sample"
          ]
        ) == s
      ]

      if (length(cells_s) == 0) {
        return(
          rep(
            0,
            nrow(counts)
          )
        )
      }

      Matrix::rowSums(
        counts[
          ,
          cells_s,
          drop=FALSE
        ]
      )
    }
  )

  rownames(out) <- rownames(counts)
  colnames(out) <- sample_order

  out
}

pb_list <- lapply(
  reference_cells,
  aggregate_reference
)

# =========================================================
# Library-size CPM
#
# Direct CPM is used because the purpose is lineage
# specificity/ambient-risk comparison rather than DE testing.
# =========================================================

to_logcpm <- function(pb) {

  libsize <- colSums(pb)

  if (any(libsize <= 0)) {
    stop("Zero pseudobulk library size detected.")
  }

  cpm_mat <- sweep(
    pb,
    2,
    libsize,
    "/"
  ) * 1e6

  log2(
    cpm_mat + 0.5
  )
}

logcpm_list <- lapply(
  pb_list,
  to_logcpm
)

# =========================================================
# Gene-level specificity table
# =========================================================

genes <- rownames(counts)

chol_lcpm <- logcpm_list[["Cholangiocyte_clean"]]

specificity <- data.frame(
  gene=genes,
  stringsAsFactors=FALSE
)

specificity$Cholangiocyte_median_logCPM <-
  apply(
    chol_lcpm,
    1,
    median
  )

for (nm in setdiff(
  names(logcpm_list),
  "Cholangiocyte_clean"
)) {

  ref <- logcpm_list[[nm]]

  ref_median <- apply(
    ref,
    1,
    median
  )

  delta_by_sample <-
    ref - chol_lcpm

  median_delta <- apply(
    delta_by_sample,
    1,
    median
  )

  specificity[[paste0(
    nm,
    "_median_logCPM"
  )]] <- ref_median

  specificity[[paste0(
    nm,
    "_minus_Chol_median_delta"
  )]] <- median_delta
}

# =========================================================
# Hepatocyte ambient-risk classification
#
# Same directional thresholds used in prior LSEC audit:
#
# > +2 : Hepatocyte dominant
# < -1 : Cholangiocyte supported vs Hepatocyte
# else : ambiguous/shared
# =========================================================

hep_delta <-
  specificity$Hepatocyte_minus_Chol_median_delta

specificity$Hepatocyte_specificity_class <-
  ifelse(
    hep_delta > 2,
    "Hepatocyte_dominant",
    ifelse(
      hep_delta < -1,
      "Cholangiocyte_supported_vs_Hepatocyte",
      "Ambiguous_shared_vs_Hepatocyte"
    )
  )

# =========================================================
# Multi-lineage review flag
#
# This is NOT an automatic exclusion criterion.
# It identifies genes much more abundant in another
# major liver lineage than in clean Cholangiocytes.
# =========================================================

reference_delta_cols <- c(
  "Hepatocyte_minus_Chol_median_delta",
  "Kupffer_Macrophage_minus_Chol_median_delta",
  "LSEC_minus_Chol_median_delta",
  "HSC_Mesenchymal_minus_Chol_median_delta",
  "Monocyte_minus_Chol_median_delta",
  "Neutrophil_minus_Chol_median_delta"
)

delta_matrix <- as.matrix(
  specificity[
    ,
    reference_delta_cols,
    drop=FALSE
  ]
)

max_index <- max.col(
  delta_matrix,
  ties.method="first"
)

specificity$max_reference_delta <-
  delta_matrix[
    cbind(
      seq_len(nrow(delta_matrix)),
      max_index
    )
  ]

reference_names <- sub(
  "_minus_Chol_median_delta$",
  "",
  reference_delta_cols
)

specificity$max_reference_lineage <-
  reference_names[
    max_index
  ]

specificity$multi_lineage_review_class <-
  ifelse(
    specificity$max_reference_delta > 2,
    paste0(
      specificity$max_reference_lineage,
      "_dominant_review"
    ),
    ifelse(
      specificity$max_reference_delta < -1,
      "Cholangiocyte_enriched_vs_all_references",
      "Shared_or_not_reference_dominant"
    )
  )

write.csv(
  specificity,
  file.path(
    TABDIR,
    "Cholangiocyte_reference_specificity_all_genes_v6.8.8.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Annotate v6.8.8 DE outputs
# =========================================================

de_files <- c(
  ALL=
    "Pseudobulk_DE_Sham_vs_Tx_ALL_v6.8.8.csv",

  NO_QCWATCH=
    "Pseudobulk_DE_Sham_vs_Tx_NO_QCWATCH_v6.8.8.csv",

  PRIMARY_CORE=
    "Pseudobulk_DE_Sham_vs_Tx_PRIMARY_CORE_v6.8.8.csv"
)

annotated_de <- list()

for (nm in names(de_files)) {

  path <- file.path(
    DE_DIR,
    de_files[[nm]]
  )

  if (!file.exists(path)) {
    stop("Missing DE table: ", path)
  }

  de <- read.csv(
    path,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  idx <- match(
    de$gene,
    specificity$gene
  )

  spec_cols <- setdiff(
    colnames(specificity),
    "gene"
  )

  out <- cbind(
    de,
    specificity[
      idx,
      spec_cols,
      drop=FALSE
    ]
  )

  annotated_de[[nm]] <- out

  write.csv(
    out,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_DE_",
        nm,
        "_reference_annotated_v6.8.8.1.csv"
      )
    ),
    row.names=FALSE
  )

  top50 <- head(
    out,
    50
  )

  write.csv(
    top50,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_TOP50_",
        nm,
        "_reference_annotated_v6.8.8.1.csv"
      )
    ),
    row.names=FALSE
  )
}

# =========================================================
# TOP50 specificity summary
# =========================================================

cat("\n=== TOP50 HEPATOCYTE SPECIFICITY COUNTS ===\n")

for (nm in names(annotated_de)) {

  top50 <- head(
    annotated_de[[nm]],
    50
  )

  cat("\n[", nm, "]\n", sep="")

  print(
    table(
      top50$Hepatocyte_specificity_class
    )
  )
}

# =========================================================
# PRIMARY_CORE FDR < 0.10
# =========================================================

primary <- annotated_de[["PRIMARY_CORE"]]

primary_fdr10 <- primary[
  is.finite(primary$FDR) &
    primary$FDR < 0.10,
  ,
  drop=FALSE
]

write.csv(
  primary_fdr10,
  file.path(
    TABDIR,
    "PRIMARY_CORE_FDRlt0.10_reference_annotated_v6.8.8.1.csv"
  ),
  row.names=FALSE
)

cat("\n=== PRIMARY_CORE FDR < 0.10 ANNOTATED ===\n")

if (nrow(primary_fdr10) == 0) {

  cat("No genes\n")

} else {

  print(
    primary_fdr10[
      ,
      c(
        "gene",
        "logFC",
        "FDR",
        "replicate_pattern",
        "Hepatocyte_minus_Chol_median_delta",
        "Hepatocyte_specificity_class",
        "max_reference_lineage",
        "max_reference_delta",
        "multi_lineage_review_class"
      )
    ],
    row.names=FALSE
  )
}

# =========================================================
# Primary top 30
# =========================================================

cat("\n=== PRIMARY_CORE TOP30 REFERENCE AUDIT ===\n")

primary_top30 <- head(
  primary,
  30
)

print(
  primary_top30[
    ,
    c(
      "gene",
      "logFC",
      "FDR",
      "replicate_pattern",
      "Hepatocyte_minus_Chol_median_delta",
      "Hepatocyte_specificity_class",
      "max_reference_lineage",
      "max_reference_delta",
      "multi_lineage_review_class"
    )
  ],
  row.names=FALSE
)

# =========================================================
# Candidate lists for downstream pathway analysis
#
# Primary ambient-aware universe:
# remove only Hepatocyte-dominant genes.
#
# Multi-lineage dominance remains a review flag rather
# than an automatic exclusion.
# =========================================================

for (nm in names(annotated_de)) {

  x <- annotated_de[[nm]]

  ambient_aware <- x[
    x$Hepatocyte_specificity_class !=
      "Hepatocyte_dominant",
    ,
    drop=FALSE
  ]

  write.csv(
    ambient_aware,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_DE_",
        nm,
        "_HEPATOCYTE_AMBIENT_AWARE_v6.8.8.1.csv"
      )
    ),
    row.names=FALSE
  )
}

summary_lines <- c(
  "# Mouse MASH Cholangiocyte reference specificity audit v6.8.8.1",
  "",
  "- Clean Cholangiocytes are compared against whole-cell reference lineages using raw-count pseudobulk CPM.",
  "- Reference lineages: Hepatocyte, Kupffer/Macrophage, LSEC, HSC/Mesenchymal, Monocyte, Neutrophil.",
  "- Hepatocyte ambient-risk classification follows the prior LSEC logic.",
  "- Hepatocyte minus Cholangiocyte median logCPM > +2: Hepatocyte_dominant.",
  "- Hepatocyte minus Cholangiocyte median logCPM < -1: Cholangiocyte_supported_vs_Hepatocyte.",
  "- Other values are Ambiguous_shared_vs_Hepatocyte.",
  "- Only Hepatocyte_dominant is automatically removed from the primary ambient-aware DE universe.",
  "- Other lineage dominance is retained as a review flag rather than automatically excluded.",
  "- Sham vs Tx inference remains limited by n=2 biological samples per group."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_reference_specificity_audit_summary_v6.8.8.1.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.8.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.8.1 COMPLETE\n")
cat("Hepatocyte ambient-risk audit complete\n")
cat("Multi-lineage review flags complete\n")
cat("No cells changed\n")
cat("No annotation changed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
