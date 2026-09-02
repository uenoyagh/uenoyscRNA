suppressPackageStartupMessages({
  library(Matrix)
})

VERSION <- "v6.9.6.4"

INPUT_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.6.3/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(
  OUTDIR,
  "tables"
)

dir.create(
  TABDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

cat("====================================================\n")
cat("Mouse MASH Monocyte support-aware reference specificity audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

SPEC_FILE <- file.path(
  INPUT_DIR,
  "Monocyte_genomewide_comprehensive_reference_specificity_v6.9.6.3.csv"
)

REFCOUNT_FILE <- file.path(
  INPUT_DIR,
  "Reference_lineage_cell_counts_v6.9.6.3.csv"
)

if (!file.exists(SPEC_FILE)) {
  stop("Missing v6.9.6.3 specificity table: ", SPEC_FILE)
}

if (!file.exists(REFCOUNT_FILE)) {
  stop("Missing v6.9.6.3 reference cell-count table: ", REFCOUNT_FILE)
}

x <- read.csv(
  SPEC_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

ref_counts <- read.csv(
  REFCOUNT_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

# =========================================================
# Primary support rule
#
# A competing lineage is eligible for the PRIMARY Sham/Tx
# specificity filter only when it has at least 20 annotated
# cells in EACH of Sham1, Sham20, Tx17 and Tx5.
#
# Sparse references are retained in the source audit but are
# not allowed to determine the primary max-reference delta.
# =========================================================

MIN_CELLS_PER_SHAMTX_SAMPLE <- 20L

required_sample_cols <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

missing_sample_cols <- setdiff(
  required_sample_cols,
  colnames(ref_counts)
)

if (length(missing_sample_cols) > 0) {
  stop(
    "Missing reference-count sample columns: ",
    paste(missing_sample_cols, collapse=", ")
  )
}

ref_counts$min_cells_ShTx <- apply(
  ref_counts[
    ,
    required_sample_cols,
    drop=FALSE
  ],
  1,
  min,
  na.rm=TRUE
)

ref_counts$primary_reference_supported <-
  ref_counts$min_cells_ShTx >=
  MIN_CELLS_PER_SHAMTX_SAMPLE

ref_counts$reference_policy <- ifelse(
  ref_counts$primary_reference_supported,
  "PRIMARY_SUPPORTED_REFERENCE",
  "SUPPLEMENTARY_SPARSE_REFERENCE_ONLY"
)

write.csv(
  ref_counts,
  file.path(
    TABDIR,
    "Reference_lineage_support_policy_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

supported_refs <- ref_counts$reference[
  ref_counts$primary_reference_supported
]

sparse_refs <- ref_counts$reference[
  !ref_counts$primary_reference_supported
]

cat(
  "Primary supported references:",
  paste(supported_refs, collapse=", "),
  "\n"
)

cat(
  "Sparse supplementary references:",
  paste(sparse_refs, collapse=", "),
  "\n\n"
)

if (length(supported_refs) == 0) {
  stop("No supported reference lineages.")
}

# =========================================================
# Resolve supported Sham/Tx delta columns
# =========================================================

supported_delta_cols <- paste0(
  supported_refs,
  "_minus_Monocyte_ShTx"
)

missing_delta_cols <- setdiff(
  supported_delta_cols,
  colnames(x)
)

if (length(missing_delta_cols) > 0) {
  stop(
    "Missing expected supported-reference delta columns: ",
    paste(missing_delta_cols, collapse=", ")
  )
}

delta_mat <- as.matrix(
  x[
    ,
    supported_delta_cols,
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

max_idx <- safe_which_max(
  delta_mat
)

x$max_supported_reference_delta_ShTx <-
  safe_row_max(
    delta_mat
  )

x$max_supported_reference_lineage_ShTx <-
  ifelse(
    is.na(max_idx),
    NA_character_,
    sub(
      "_minus_Monocyte_ShTx$",
      "",
      supported_delta_cols[
        max_idx
      ]
    )
  )

# =========================================================
# Final primary specificity classes
# =========================================================

x$specificity_class_supported_ShTx <-
  ifelse(
    is.na(
      x$max_supported_reference_delta_ShTx
    ),
    "unresolved",
    ifelse(
      x$max_supported_reference_delta_ShTx > 2,
      "reference_lineage_dominant",
      ifelse(
        x$max_supported_reference_delta_ShTx < -1,
        "Monocyte_enriched",
        "ambiguous_shared"
      )
    )
  )

x$keep_reference_aware_supported <-
  !is.na(
    x$max_supported_reference_delta_ShTx
  ) &
  x$max_supported_reference_delta_ShTx <= 2

x$keep_strict_reference_supported <-
  !is.na(
    x$max_supported_reference_delta_ShTx
  ) &
  x$max_supported_reference_delta_ShTx <= 1.5

# =========================================================
# Compare v6.9.6.3 all-reference filter vs support-aware
# =========================================================

if (!("keep_reference_aware" %in% colnames(x))) {
  stop("Missing v6.9.6.3 keep_reference_aware column.")
}

if (!("keep_strict_reference_sensitivity" %in% colnames(x))) {
  stop("Missing v6.9.6.3 strict reference column.")
}

x$reference_filter_change <-
  ifelse(
    x$keep_reference_aware ==
      x$keep_reference_aware_supported,
    "unchanged",
    ifelse(
      !x$keep_reference_aware &
      x$keep_reference_aware_supported,
      "rescued_after_sparse_reference_exclusion",
      "newly_excluded"
    )
  )

write.csv(
  x,
  file.path(
    TABDIR,
    "Monocyte_genomewide_support_aware_reference_specificity_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Summary
# =========================================================

specificity_summary <- as.data.frame(
  table(
    specificity_class_supported_ShTx=
      x$specificity_class_supported_ShTx
  ),
  stringsAsFactors=FALSE
)

write.csv(
  specificity_summary,
  file.path(
    TABDIR,
    "Monocyte_support_aware_specificity_class_counts_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

change_summary <- as.data.frame(
  table(
    reference_filter_change=
      x$reference_filter_change
  ),
  stringsAsFactors=FALSE
)

write.csv(
  change_summary,
  file.path(
    TABDIR,
    "Monocyte_reference_filter_change_summary_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# DE-tested filter counts
# =========================================================

required_de_cols <- c(
  "ALL_logFC",
  "ALL_PValue",
  "NO_QCWATCH_logFC",
  "NO_QCWATCH_PValue",
  "PRIMARY_CORE_logFC",
  "PRIMARY_CORE_PValue",
  "PRIMARY_CORE_FDR",
  "PRIMARY_CORE_replicate_direction"
)

missing_de_cols <- setdiff(
  required_de_cols,
  colnames(x)
)

if (length(missing_de_cols) > 0) {
  stop(
    "Missing DE columns: ",
    paste(missing_de_cols, collapse=", ")
  )
}

tested <- x[
  is.finite(
    x$PRIMARY_CORE_PValue
  ),
  ,
  drop=FALSE
]

filter_summary <- data.frame(
  metric=c(
    "PRIMARY_CORE_tested_genes",
    "REFERENCE_AWARE_SUPPORTED_keep_delta_le2",
    "STRICT_SUPPORTED_keep_delta_le1.5",
    "SUPPORTED_REFERENCE_DOMINANT_delta_gt2",
    "RESCUED_vs_v6.9.6.3_reference_aware",
    "NEWLY_EXCLUDED_vs_v6.9.6.3_reference_aware"
  ),
  n_genes=c(
    nrow(tested),

    sum(
      tested$keep_reference_aware_supported,
      na.rm=TRUE
    ),

    sum(
      tested$keep_strict_reference_supported,
      na.rm=TRUE
    ),

    sum(
      tested$max_supported_reference_delta_ShTx > 2,
      na.rm=TRUE
    ),

    sum(
      tested$reference_filter_change ==
        "rescued_after_sparse_reference_exclusion",
      na.rm=TRUE
    ),

    sum(
      tested$reference_filter_change ==
        "newly_excluded",
      na.rm=TRUE
    )
  ),
  stringsAsFactors=FALSE
)

write.csv(
  filter_summary,
  file.path(
    TABDIR,
    "Monocyte_GSEA_support_aware_filter_gene_counts_v6.9.6.4.csv"
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

key_table <- x[
  x$gene %in% key_genes,
  ,
  drop=FALSE
]

write.csv(
  key_table,
  file.path(
    TABDIR,
    "Monocyte_key_gene_support_aware_specificity_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Top30 PRIMARY_CORE DE
# =========================================================

top30 <- x[
  is.finite(
    x$PRIMARY_CORE_PValue
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
    "Monocyte_PRIMARY_CORE_top30_support_aware_specificity_v6.9.6.4.csv"
  ),
  row.names=FALSE
)

# =========================================================
# FINAL rank tables for v6.9.7
#
# Rank metric:
# sign(logFC) * -log10(PValue)
#
# Reference-aware is primary.
# Strict is sensitivity.
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

  rank_base <- x[
    is.finite(
      x[[logfc_col]]
    ) &
    is.finite(
      x[[p_col]]
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

  reference_aware <- rank_base[
    rank_base$keep_reference_aware_supported,
    c(
      "gene",
      logfc_col,
      p_col,
      "rank_metric",
      "max_supported_reference_lineage_ShTx",
      "max_supported_reference_delta_ShTx",
      "specificity_class_supported_ShTx"
    ),
    drop=FALSE
  ]

  strict <- rank_base[
    rank_base$keep_strict_reference_supported,
    c(
      "gene",
      logfc_col,
      p_col,
      "rank_metric",
      "max_supported_reference_lineage_ShTx",
      "max_supported_reference_delta_ShTx",
      "specificity_class_supported_ShTx"
    ),
    drop=FALSE
  ]

  reference_aware <- reference_aware[
    order(
      -reference_aware$rank_metric
    ),
    ,
    drop=FALSE
  ]

  strict <- strict[
    order(
      -strict$rank_metric
    ),
    ,
    drop=FALSE
  ]

  write.csv(
    reference_aware,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_",
        fw,
        "_FINAL_reference_aware_rank_v6.9.6.4.csv"
      )
    ),
    row.names=FALSE
  )

  write.csv(
    strict,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_",
        fw,
        "_FINAL_strict_reference_rank_v6.9.6.4.csv"
      )
    ),
    row.names=FALSE
  )
}

# =========================================================
# Terminal output
# =========================================================

cat("\n=== REFERENCE SUPPORT POLICY ===\n")
print(
  ref_counts[
    ,
    c(
      "reference",
      "Sham1",
      "Sham20",
      "Tx17",
      "Tx5",
      "min_cells_ShTx",
      "reference_policy"
    )
  ],
  row.names=FALSE
)

cat("\n=== SUPPORT-AWARE SPECIFICITY CLASS COUNTS ===\n")
print(
  specificity_summary,
  row.names=FALSE
)

cat("\n=== REFERENCE FILTER CHANGE SUMMARY ===\n")
print(
  change_summary,
  row.names=FALSE
)

cat("\n=== GSEA SUPPORT-AWARE FILTER GENE COUNTS ===\n")
print(
  filter_summary,
  row.names=FALSE
)

cat("\n=== KEY GENE SUPPORT-AWARE SPECIFICITY ===\n")

key_cols <- c(
  "gene",
  "Monocyte_median_logCPM_ShTx",
  "max_supported_reference_lineage_ShTx",
  "max_supported_reference_delta_ShTx",
  "specificity_class_supported_ShTx",
  "keep_reference_aware_supported",
  "keep_strict_reference_supported",
  "reference_filter_change",
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

cat("\n====================================================\n")
cat("v6.9.6.4 COMPLETE\n")
cat("Support-aware comprehensive reference audit complete\n")
cat("Primary competing-reference eligibility: >=20 cells in EACH Sham/Tx sample\n")
cat("Sparse references remain supplementary but do not determine primary max delta\n")
cat("Reference-aware filter: supported max delta <= 2\n")
cat("Strict sensitivity filter: supported max delta <= 1.5\n")
cat("FINAL rank tables generated for v6.9.7\n")
cat("No source RDS modified\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.6.4.txt"
  )
)
