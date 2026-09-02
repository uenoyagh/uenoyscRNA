suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.6.3"

MON_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.5/objects/",
  "Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds"
)

WHOLE_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

DE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.6/tables"
)

ANNOTATION_COL <- "celltype_for_R8plot_FIXED2"

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte comprehensive reference specificity audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(MON_RDS)) {
  stop("Missing Monocyte RDS: ", MON_RDS)
}

if (!file.exists(WHOLE_RDS)) {
  stop("Missing whole-cell RDS: ", WHOLE_RDS)
}

# =========================================================
# Helpers
# =========================================================

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

aggregate_logcpm <- function(
  mat,
  samples,
  genes,
  sample_order
) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  out <- matrix(
    NA_real_,
    nrow=length(genes),
    ncol=length(sample_order),
    dimnames=list(
      genes,
      sample_order
    )
  )

  for (s in sample_order) {

    cells <- colnames(mat)[
      samples == s
    ]

    if (length(cells) == 0) {
      next
    }

    pb_all <- Matrix::rowSums(
      mat[
        ,
        cells,
        drop=FALSE
      ]
    )

    lib <- sum(pb_all)

    cpm <- 1e6 * pb_all /
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

row_median_na <- function(x) {

  apply(
    x,
    1,
    function(v) {

      if (all(is.na(v))) {
        return(NA_real_)
      }

      median(
        v,
        na.rm=TRUE
      )
    }
  )
}

resolve_reference_label <- function(
  available_labels,
  candidates
) {

  hit <- candidates[
    candidates %in% available_labels
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[[1]]
}

# =========================================================
# Sample definitions
# =========================================================

samples_all6 <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

samples_shamtx <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

# =========================================================
# Monocyte
# =========================================================

mon <- readRDS(
  MON_RDS
)

mon_counts <- get_counts(
  mon
)

mon_sample_col <- resolve_sample_col(
  mon
)

mon_sample <- as.character(
  mon@meta.data[[mon_sample_col]]
)

genes <- rownames(
  mon_counts
)

cat(
  "Genome-wide audit genes:",
  length(genes),
  "\n"
)

mon_all6 <- aggregate_logcpm(
  mat=mon_counts,
  samples=mon_sample,
  genes=genes,
  sample_order=samples_all6
)

mon_shamtx <- aggregate_logcpm(
  mat=mon_counts,
  samples=mon_sample,
  genes=genes,
  sample_order=samples_shamtx
)

# =========================================================
# Whole-cell reference labels
# =========================================================

whole <- readRDS(
  WHOLE_RDS
)

if (!(ANNOTATION_COL %in% colnames(whole@meta.data))) {
  stop(
    "Missing annotation column: ",
    ANNOTATION_COL
  )
}

labels <- as.character(
  whole@meta.data[[ANNOTATION_COL]]
)

available_labels <- sort(
  unique(labels)
)

label_table <- sort(
  table(labels),
  decreasing=TRUE
)

cat("\n=== AVAILABLE WHOLE-CELL LABELS ===\n")
print(label_table)

write.csv(
  data.frame(
    label=names(label_table),
    n_cells=as.integer(label_table),
    stringsAsFactors=FALSE
  ),
  file.path(
    TABDIR,
    "Whole_cell_annotation_labels_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Comprehensive competing-lineage references
#
# Monocyte itself is intentionally excluded.
# Cycling is intentionally excluded because it is a state,
# not an independent lineage reference.
# Aliases are supplied for labels that may differ between
# historical annotation versions.
# =========================================================

reference_candidates <- list(

  Hepatocyte=c(
    "Hepatocyte"
  ),

  Kupffer_Macrophage=c(
    "Kupffer_Macrophage"
  ),

  HSC_Mesenchymal=c(
    "HSC_Mesenchymal"
  ),

  LSEC=c(
    "LSEC"
  ),

  Cholangiocyte=c(
    "Cholangiocyte"
  ),

  Neutrophil=c(
    "Neutrophil"
  ),

  Dendritic_cell=c(
    "Dendritic_cell",
    "Dendritic"
  ),

  B_cell=c(
    "B_cell",
    "B"
  ),

  T_cell=c(
    "T_cell",
    "T"
  ),

  NK_cell=c(
    "NK_cell",
    "NK"
  ),

  Vascular_endothelial=c(
    "Vascular_endothelial",
    "Vascular",
    "Endothelial_Vascular",
    "Endothelial"
  ),

  Plasma_cell=c(
    "Plasma_cell",
    "Plasma"
  ),

  Mesothelial=c(
    "Mesothelial"
  ),

  RBC=c(
    "RBC",
    "Erythrocyte"
  ),

  Platelet=c(
    "Platelet"
  )
)

resolved_labels <- vapply(
  reference_candidates,
  function(x) {
    resolve_reference_label(
      available_labels,
      x
    )
  },
  character(1)
)

reference_resolution <- data.frame(
  reference=names(resolved_labels),
  resolved_annotation_label=unname(resolved_labels),
  available=!is.na(resolved_labels),
  stringsAsFactors=FALSE
)

write.csv(
  reference_resolution,
  file.path(
    TABDIR,
    "Reference_label_resolution_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

cat("\n=== REFERENCE LABEL RESOLUTION ===\n")
print(
  reference_resolution,
  row.names=FALSE
)

missing_refs <- reference_resolution$reference[
  !reference_resolution$available
]

if (length(missing_refs) > 0) {

  warning(
    "Some expected reference lineages were not found: ",
    paste(
      missing_refs,
      collapse=", "
    )
  )
}

reference_labels <- resolved_labels[
  !is.na(resolved_labels)
]

if (length(reference_labels) == 0) {
  stop("No competing reference lineages resolved.")
}

whole_sample_col <- resolve_sample_col(
  whole
)

reference_all6 <- list()
reference_shamtx <- list()
reference_counts <- list()

# =========================================================
# Process each reference one at a time to limit memory use
# =========================================================

for (ref_name in names(reference_labels)) {

  label <- unname(
    reference_labels[[ref_name]]
  )

  cells <- rownames(
    whole@meta.data
  )[
    labels == label
  ]

  if (length(cells) == 0) {
    next
  }

  ref <- subset(
    whole,
    cells=cells
  )

  ref_counts <- get_counts(
    ref
  )

  ref_sample <- as.character(
    ref@meta.data[[whole_sample_col]]
  )

  reference_all6[[ref_name]] <-
    aggregate_logcpm(
      mat=ref_counts,
      samples=ref_sample,
      genes=genes,
      sample_order=samples_all6
    )

  reference_shamtx[[ref_name]] <-
    aggregate_logcpm(
      mat=ref_counts,
      samples=ref_sample,
      genes=genes,
      sample_order=samples_shamtx
    )

  sample_count_table <- table(
    factor(
      ref_sample,
      levels=samples_all6
    )
  )

  reference_counts[[ref_name]] <-
    data.frame(
      reference=ref_name,
      annotation_label=label,
      n_cells=ncol(ref),
      STD_rep1=as.integer(sample_count_table["STD_rep1"]),
      CDHFD_rep1=as.integer(sample_count_table["CDHFD_rep1"]),
      Sham1=as.integer(sample_count_table["Sham1"]),
      Sham20=as.integer(sample_count_table["Sham20"]),
      Tx17=as.integer(sample_count_table["Tx17"]),
      Tx5=as.integer(sample_count_table["Tx5"]),
      stringsAsFactors=FALSE
    )

  cat(
    "Reference:",
    ref_name,
    "| label:",
    label,
    "| cells:",
    ncol(ref),
    "\n"
  )

  rm(
    ref,
    ref_counts
  )

  gc()
}

reference_counts <- do.call(
  rbind,
  reference_counts
)

rownames(reference_counts) <- NULL

write.csv(
  reference_counts,
  file.path(
    TABDIR,
    "Reference_lineage_cell_counts_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# Whole object is no longer needed.
rm(whole)
gc()

# =========================================================
# Genome-wide specificity comparison
# =========================================================

comparison <- data.frame(
  gene=genes,

  Monocyte_median_logCPM_ShTx=
    row_median_na(
      mon_shamtx
    )[genes],

  Monocyte_median_logCPM_all6=
    row_median_na(
      mon_all6
    )[genes],

  stringsAsFactors=FALSE
)

for (ref_name in names(reference_shamtx)) {

  ref_st <- row_median_na(
    reference_shamtx[[ref_name]]
  )

  ref_a6 <- row_median_na(
    reference_all6[[ref_name]]
  )

  comparison[[
    paste0(
      ref_name,
      "_median_logCPM_ShTx"
    )
  ]] <- ref_st[
    comparison$gene
  ]

  comparison[[
    paste0(
      ref_name,
      "_minus_Monocyte_ShTx"
    )
  ]] <-
    comparison[[
      paste0(
        ref_name,
        "_median_logCPM_ShTx"
      )
    ]] -
    comparison$Monocyte_median_logCPM_ShTx

  comparison[[
    paste0(
      ref_name,
      "_median_logCPM_all6"
    )
  ]] <- ref_a6[
    comparison$gene
  ]

  comparison[[
    paste0(
      ref_name,
      "_minus_Monocyte_all6"
    )
  ]] <-
    comparison[[
      paste0(
        ref_name,
        "_median_logCPM_all6"
      )
    ]] -
    comparison$Monocyte_median_logCPM_all6
}

# =========================================================
# Maximum competing-lineage dominance
# =========================================================

shamtx_delta_cols <- grep(
  "_minus_Monocyte_ShTx$",
  colnames(comparison),
  value=TRUE
)

all6_delta_cols <- grep(
  "_minus_Monocyte_all6$",
  colnames(comparison),
  value=TRUE
)

if (length(shamtx_delta_cols) == 0) {
  stop("No Sham/Tx reference delta columns were generated.")
}

shamtx_delta_mat <- as.matrix(
  comparison[
    ,
    shamtx_delta_cols,
    drop=FALSE
  ]
)

all6_delta_mat <- as.matrix(
  comparison[
    ,
    all6_delta_cols,
    drop=FALSE
  ]
)

safe_row_max <- function(mat) {

  apply(
    mat,
    1,
    function(v) {

      if (all(is.na(v))) {
        return(NA_real_)
      }

      max(
        v,
        na.rm=TRUE
      )
    }
  )
}

safe_which_max <- function(mat) {

  apply(
    mat,
    1,
    function(v) {

      if (all(is.na(v))) {
        return(NA_integer_)
      }

      which.max(v)
    }
  )
}

shamtx_max_idx <- safe_which_max(
  shamtx_delta_mat
)

all6_max_idx <- safe_which_max(
  all6_delta_mat
)

comparison$max_reference_delta_ShTx <-
  safe_row_max(
    shamtx_delta_mat
  )

comparison$max_reference_lineage_ShTx <-
  ifelse(
    is.na(shamtx_max_idx),
    NA_character_,
    sub(
      "_minus_Monocyte_ShTx$",
      "",
      shamtx_delta_cols[
        shamtx_max_idx
      ]
    )
  )

comparison$max_reference_delta_all6 <-
  safe_row_max(
    all6_delta_mat
  )

comparison$max_reference_lineage_all6 <-
  ifelse(
    is.na(all6_max_idx),
    NA_character_,
    sub(
      "_minus_Monocyte_all6$",
      "",
      all6_delta_cols[
        all6_max_idx
      ]
    )
  )

# =========================================================
# Specificity classes
#
# Primary decision is Sham/Tx-matched because v6.9.6 DE
# compares Sham vs Tx.
# =========================================================

comparison$specificity_class_ShTx <-
  ifelse(
    is.na(
      comparison$max_reference_delta_ShTx
    ),
    "unresolved",
    ifelse(
      comparison$max_reference_delta_ShTx > 2,
      "reference_lineage_dominant",
      ifelse(
        comparison$max_reference_delta_ShTx < -1,
        "Monocyte_enriched",
        "ambiguous_shared"
      )
    )
  )

comparison$specificity_class_all6 <-
  ifelse(
    is.na(
      comparison$max_reference_delta_all6
    ),
    "unresolved",
    ifelse(
      comparison$max_reference_delta_all6 > 2,
      "reference_lineage_dominant",
      ifelse(
        comparison$max_reference_delta_all6 < -1,
        "Monocyte_enriched",
        "ambiguous_shared"
      )
    )
  )

# =========================================================
# GSEA filter flags
# =========================================================

comparison$keep_reference_aware <-
  !is.na(
    comparison$max_reference_delta_ShTx
  ) &
  comparison$max_reference_delta_ShTx <= 2

comparison$keep_strict_reference_sensitivity <-
  !is.na(
    comparison$max_reference_delta_ShTx
  ) &
  comparison$max_reference_delta_ShTx <= 1.5

# =========================================================
# Add v6.9.6 DE information
# =========================================================

de_files <- c(
  ALL=file.path(
    DE_DIR,
    "Monocyte_pseudobulk_ALL_Sham_vs_Tx_v6.9.6.csv"
  ),
  NO_QCWATCH=file.path(
    DE_DIR,
    "Monocyte_pseudobulk_NO_QCWATCH_Sham_vs_Tx_v6.9.6.csv"
  ),
  PRIMARY_CORE=file.path(
    DE_DIR,
    "Monocyte_pseudobulk_PRIMARY_CORE_Sham_vs_Tx_v6.9.6.csv"
  )
)

for (fw in names(de_files)) {

  if (!file.exists(
    de_files[[fw]]
  )) {
    stop(
      "Missing DE file: ",
      de_files[[fw]]
    )
  }

  de <- read.csv(
    de_files[[fw]],
    stringsAsFactors=FALSE
  )

  idx <- match(
    comparison$gene,
    de$gene
  )

  comparison[[
    paste0(
      fw,
      "_logFC"
    )
  ]] <- de$logFC[idx]

  comparison[[
    paste0(
      fw,
      "_PValue"
    )
  ]] <- de$PValue[idx]

  comparison[[
    paste0(
      fw,
      "_FDR"
    )
  ]] <- de$FDR[idx]

  comparison[[
    paste0(
      fw,
      "_replicate_direction"
    )
  ]] <- de$replicate_direction[idx]
}

# =========================================================
# Save genome-wide table
# =========================================================

write.csv(
  comparison,
  file.path(
    TABDIR,
    "Monocyte_genomewide_comprehensive_reference_specificity_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Specificity summary
# =========================================================

specificity_summary <- as.data.frame(
  table(
    specificity_class_ShTx=
      comparison$specificity_class_ShTx
  ),
  stringsAsFactors=FALSE
)

write.csv(
  specificity_summary,
  file.path(
    TABDIR,
    "Monocyte_comprehensive_specificity_class_counts_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# GSEA filter counts
# =========================================================

tested <- comparison[
  is.finite(
    comparison$PRIMARY_CORE_PValue
  ),
  ,
  drop=FALSE
]

filter_summary <- data.frame(
  metric=c(
    "PRIMARY_CORE_tested_genes",
    "REFERENCE_AWARE_keep_delta_le2",
    "STRICT_keep_delta_le1.5",
    "REFERENCE_DOMINANT_delta_gt2"
  ),

  n_genes=c(
    nrow(tested),

    sum(
      tested$keep_reference_aware,
      na.rm=TRUE
    ),

    sum(
      tested$keep_strict_reference_sensitivity,
      na.rm=TRUE
    ),

    sum(
      tested$max_reference_delta_ShTx > 2,
      na.rm=TRUE
    )
  ),

  stringsAsFactors=FALSE
)

write.csv(
  filter_summary,
  file.path(
    TABDIR,
    "Monocyte_GSEA_filter_gene_counts_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Key genes
# =========================================================

key_genes <- c(
  "Gdf15",
  "Thbs1",
  "Slc25a47",
  "Rgs1",
  "Cxcr4",
  "Ccn1",
  "Ighm",
  "Cyp2c29",
  "Pdk4",
  "Krt8",
  "Gstm3",
  "Ddit4",
  "Zfp36",
  "Cebpb",
  "Cxcl2"
)

key_table <- comparison[
  comparison$gene %in%
    key_genes,
  ,
  drop=FALSE
]

write.csv(
  key_table,
  file.path(
    TABDIR,
    "Monocyte_key_gene_comprehensive_specificity_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Top30 PRIMARY_CORE DE with comprehensive specificity
# =========================================================

top30 <- comparison[
  is.finite(
    comparison$PRIMARY_CORE_PValue
  ),
  ,
  drop=FALSE
]

top30 <- top30[
  order(
    top30$PRIMARY_CORE_PValue
  ),
  ,
  drop=FALSE
]

top30 <- head(
  top30,
  30
)

write.csv(
  top30,
  file.path(
    TABDIR,
    "Monocyte_PRIMARY_CORE_top30_comprehensive_specificity_v6.9.6.3.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Reference-aware DE rank tables for next GSEA step
#
# No pathway analysis is performed here.
# These are audit-ready rank inputs only.
# =========================================================

for (fw in c(
  "ALL",
  "NO_QCWATCH",
  "PRIMARY_CORE"
)) {

  logfc_col <- paste0(
    fw,
    "_logFC"
  )

  p_col <- paste0(
    fw,
    "_PValue"
  )

  rank_base <- comparison[
    is.finite(
      comparison[[p_col]]
    ) &
    is.finite(
      comparison[[logfc_col]]
    ),
    ,
    drop=FALSE
  ]

  rank_base$rank_metric <-
    sign(
      rank_base[[logfc_col]]
    ) *
    -log10(
      pmax(
        rank_base[[p_col]],
        .Machine$double.xmin
      )
    )

  rank_reference_aware <- rank_base[
    rank_base$keep_reference_aware,
    c(
      "gene",
      logfc_col,
      p_col,
      "rank_metric",
      "max_reference_lineage_ShTx",
      "max_reference_delta_ShTx",
      "specificity_class_ShTx"
    ),
    drop=FALSE
  ]

  rank_strict <- rank_base[
    rank_base$keep_strict_reference_sensitivity,
    c(
      "gene",
      logfc_col,
      p_col,
      "rank_metric",
      "max_reference_lineage_ShTx",
      "max_reference_delta_ShTx",
      "specificity_class_ShTx"
    ),
    drop=FALSE
  ]

  rank_reference_aware <- rank_reference_aware[
    order(
      -rank_reference_aware$rank_metric
    ),
    ,
    drop=FALSE
  ]

  rank_strict <- rank_strict[
    order(
      -rank_strict$rank_metric
    ),
    ,
    drop=FALSE
  ]

  write.csv(
    rank_reference_aware,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_",
        fw,
        "_reference_aware_rank_v6.9.6.3.csv"
      )
    ),
    row.names=FALSE
  )

  write.csv(
    rank_strict,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_",
        fw,
        "_strict_reference_rank_v6.9.6.3.csv"
      )
    ),
    row.names=FALSE
  )
}

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== REFERENCE LINEAGE CELL COUNTS ===\n")
print(
  reference_counts,
  row.names=FALSE
)

cat("\n=== COMPREHENSIVE SPECIFICITY CLASS COUNTS ===\n")
print(
  specificity_summary,
  row.names=FALSE
)

cat("\n=== GSEA FILTER GENE COUNTS ===\n")
print(
  filter_summary,
  row.names=FALSE
)

cat("\n=== KEY GENE COMPREHENSIVE SPECIFICITY ===\n")

key_cols <- c(
  "gene",
  "Monocyte_median_logCPM_ShTx",
  "max_reference_lineage_ShTx",
  "max_reference_delta_ShTx",
  "specificity_class_ShTx",
  "keep_reference_aware",
  "keep_strict_reference_sensitivity",
  "PRIMARY_CORE_logFC",
  "PRIMARY_CORE_FDR",
  "PRIMARY_CORE_replicate_direction"
)

print(
  key_table[
    ,
    intersect(
      key_cols,
      colnames(key_table)
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== TOP30 PRIMARY_CORE COMPREHENSIVE SPECIFICITY ===\n")

print(
  top30[
    ,
    intersect(
      key_cols,
      colnames(top30)
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.6.3 COMPLETE\n")
cat("Comprehensive genome-wide reference specificity audit complete\n")
cat("Primary specificity comparison: Sham1/Sham20/Tx17/Tx5\n")
cat("Monocyte excluded from competing references by design\n")
cat("Cycling excluded because it is a state rather than a lineage\n")
cat("Reference-aware filter: max competing-lineage delta <= 2\n")
cat("Strict sensitivity filter: max competing-lineage delta <= 1.5\n")
cat("No genes removed from source data\n")
cat("Reference-aware and strict rank tables generated for v6.9.7\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.6.3.txt"
  )
)
