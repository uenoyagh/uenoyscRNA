suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.6.1"

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
cat("Mouse MASH Monocyte reference specificity audit\n")
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
    candidates %in%
      colnames(obj@meta.data)
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

    if (length(cells) == 0) {
      next
    }

    pb <- Matrix::rowSums(
      mat[
        ,
        cells,
        drop=FALSE
      ]
    )

    lib <- sum(pb)

    cpm <- 1e6 * pb /
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

# =========================================================
# Read DE results
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

missing_de <- de_files[
  !file.exists(de_files)
]

if (length(missing_de) > 0) {
  stop(
    "Missing DE files: ",
    paste(missing_de, collapse=", ")
  )
}

de_list <- lapply(
  de_files,
  read.csv,
  stringsAsFactors=FALSE
)

# Top50 union across frameworks
top50_union <- unique(
  unlist(
    lapply(
      de_list,
      function(x) {
        head(
          x$gene[
            order(x$PValue)
          ],
          50
        )
      }
    )
  )
)

# Always force key candidates into audit
audit_genes <- unique(
  c(
    "Gdf15",
    "Thbs1",
    "Slc25a47",
    "Rgs1",
    "Cxcr4",
    "Ccn1",
    top50_union
  )
)

# =========================================================
# Clean Monocyte reference
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

mon_samples <- as.character(
  mon@meta.data[[mon_sample_col]]
)

audit_genes <- intersect(
  audit_genes,
  rownames(mon_counts)
)

cat(
  "Audit genes present:",
  length(audit_genes),
  "\n\n"
)

mon_lcpm <- aggregate_logcpm(
  mat=mon_counts,
  samples=mon_samples,
  genes=audit_genes
)

# =========================================================
# Whole-cell reference lineages
# =========================================================

whole <- readRDS(
  WHOLE_RDS
)

if (!(ANNOTATION_COL %in%
      colnames(whole@meta.data))) {

  stop(
    "Missing annotation column: ",
    ANNOTATION_COL
  )
}

cat("=== AVAILABLE WHOLE-CELL LABELS ===\n")

label_table <- sort(
  table(
    as.character(
      whole@meta.data[[ANNOTATION_COL]]
    )
  ),
  decreasing=TRUE
)

print(label_table)

write.csv(
  data.frame(
    label=names(label_table),
    n_cells=as.integer(label_table),
    stringsAsFactors=FALSE
  ),
  file.path(
    TABDIR,
    "Whole_cell_annotation_labels_v6.9.6.1.csv"
  ),
  row.names=FALSE
)

reference_labels <- c(
  Hepatocyte="Hepatocyte",
  Kupffer_Macrophage="Kupffer_Macrophage",
  HSC_Mesenchymal="HSC_Mesenchymal",
  LSEC="LSEC",
  Cholangiocyte="Cholangiocyte",
  Neutrophil="Neutrophil",
  Dendritic="Dendritic"
)

available_labels <- unique(
  as.character(
    whole@meta.data[[ANNOTATION_COL]]
  )
)

reference_labels <- reference_labels[
  reference_labels %in% available_labels
]

if (length(reference_labels) == 0) {
  stop("No expected reference labels found.")
}

whole_sample_col <- resolve_sample_col(
  whole
)

reference_lcpm <- list()

reference_cell_counts <- list()

for (ref_name in names(reference_labels)) {

  label <- reference_labels[[ref_name]]

  cells <- rownames(
    whole@meta.data
  )[
    as.character(
      whole@meta.data[[ANNOTATION_COL]]
    ) == label
  ]

  ref <- subset(
    whole,
    cells=cells
  )

  ref_counts <- get_counts(
    ref
  )

  ref_samples <- as.character(
    ref@meta.data[[whole_sample_col]]
  )

  reference_lcpm[[ref_name]] <-
    aggregate_logcpm(
      mat=ref_counts,
      samples=ref_samples,
      genes=audit_genes
    )

  reference_cell_counts[[ref_name]] <-
    data.frame(
      reference=ref_name,
      annotation_label=label,
      n_cells=ncol(ref),
      stringsAsFactors=FALSE
    )

  cat(
    "Reference:",
    ref_name,
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

reference_cell_counts <- do.call(
  rbind,
  reference_cell_counts
)

write.csv(
  reference_cell_counts,
  file.path(
    TABDIR,
    "Reference_lineage_cell_counts_v6.9.6.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Median expression comparison
# =========================================================

mon_median <- apply(
  mon_lcpm,
  1,
  median,
  na.rm=TRUE
)

comparison <- data.frame(
  gene=audit_genes,
  Monocyte_median_logCPM=
    mon_median[audit_genes],
  stringsAsFactors=FALSE
)

for (ref_name in names(reference_lcpm)) {

  x <- reference_lcpm[[ref_name]]

  ref_median <- apply(
    x,
    1,
    median,
    na.rm=TRUE
  )

  comparison[[
    paste0(
      ref_name,
      "_median_logCPM"
    )
  ]] <- ref_median[
    comparison$gene
  ]

  comparison[[
    paste0(
      ref_name,
      "_minus_Monocyte"
    )
  ]] <-
    comparison[[
      paste0(
        ref_name,
        "_median_logCPM"
      )
    ]] -
    comparison$Monocyte_median_logCPM
}

# =========================================================
# Maximum competing-lineage dominance
# =========================================================

delta_cols <- grep(
  "_minus_Monocyte$",
  colnames(comparison),
  value=TRUE
)

delta_mat <- as.matrix(
  comparison[
    ,
    delta_cols,
    drop=FALSE
  ]
)

comparison$max_reference_delta <-
  apply(
    delta_mat,
    1,
    max,
    na.rm=TRUE
  )

max_idx <- apply(
  delta_mat,
  1,
  which.max
)

comparison$max_reference_lineage <-
  sub(
    "_minus_Monocyte$",
    "",
    delta_cols[max_idx]
  )

comparison$specificity_class <-
  ifelse(
    comparison$max_reference_delta > 2,
    "reference_lineage_dominant",
    ifelse(
      comparison$max_reference_delta < -1,
      "Monocyte_enriched",
      "ambiguous_shared"
    )
  )

# =========================================================
# Add DE information
# =========================================================

for (fw in names(de_list)) {

  de <- de_list[[fw]]

  idx <- match(
    comparison$gene,
    de$gene
  )

  comparison[[
    paste0(
      fw,
      "_logFC"
    )
  ]] <- de$logFC[
    idx
  ]

  comparison[[
    paste0(
      fw,
      "_FDR"
    )
  ]] <- de$FDR[
    idx
  ]

  comparison[[
    paste0(
      fw,
      "_replicate_direction"
    )
  ]] <- de$replicate_direction[
    idx
  ]
}

comparison <- comparison[
  order(
    comparison$PRIMARY_CORE_FDR,
    -abs(
      comparison$PRIMARY_CORE_logFC
    )
  ),
  ,
  drop=FALSE
]

write.csv(
  comparison,
  file.path(
    TABDIR,
    "Monocyte_reference_specificity_audit_v6.9.6.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Key-candidate table
# =========================================================

key_genes <- c(
  "Gdf15",
  "Thbs1",
  "Slc25a47",
  "Rgs1",
  "Cxcr4",
  "Ccn1"
)

key_table <- comparison[
  comparison$gene %in% key_genes,
  ,
  drop=FALSE
]

write.csv(
  key_table,
  file.path(
    TABDIR,
    "Monocyte_key_DE_reference_specificity_v6.9.6.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Suggested primary-use classification
#
# Reference-lineage dominance >2 log2 units is flagged.
# No gene is deleted here.
# =========================================================

comparison$primary_interpretation <-
  ifelse(
    comparison$specificity_class ==
      "reference_lineage_dominant",
    "EXCLUDE_FROM_PRIMARY_MONOCYTE_INTERPRETATION",
    "RETAIN_FOR_MONOCYTE_INTERPRETATION"
  )

write.csv(
  comparison,
  file.path(
    TABDIR,
    "Monocyte_reference_specificity_with_primary_policy_v6.9.6.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== REFERENCE LINEAGE CELL COUNTS ===\n")
print(
  reference_cell_counts,
  row.names=FALSE
)

cat("\n=== KEY DE REFERENCE SPECIFICITY ===\n")

key_cols <- c(
  "gene",
  "Monocyte_median_logCPM",
  "max_reference_lineage",
  "max_reference_delta",
  "specificity_class",
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

cat("\n=== TOP30 PRIMARY_CORE REFERENCE SPECIFICITY ===\n")

top30 <- comparison[
  order(
    comparison$PRIMARY_CORE_FDR
  ),
  ,
  drop=FALSE
]

print(
  head(
    top30[
      ,
      intersect(
        key_cols,
        colnames(top30)
      ),
      drop=FALSE
    ],
    30
  ),
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.6.1 COMPLETE\n")
cat("Monocyte reference-lineage specificity audit complete\n")
cat("No genes removed\n")
cat("Dominance flag: competing lineage > Monocyte by >2 log2 CPM\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.6.1.txt"
  )
)
