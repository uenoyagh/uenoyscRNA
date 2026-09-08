#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

BIOLOGICAL_VERSION <- "v1.4.3A.1"
FREEZE_DATE <- "2026-09-08"

BASE <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

SOURCE_RDS <- file.path(
  BASE,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.0",
  "objects",
  "Human_MASH_12sample_Context_v1.4.0_QC2_reconciled.rds"
)

CONTEXT_META <- file.path(
  BASE,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3.1",
  "analysis",
  "v1.4.3.1_12sample_integrated_context_review",
  "tables",
  "v1.4.3_context_UMAP_metadata.tsv"
)

IL10_TABLE <- file.path(
  BASE,
  "Human_MASH_RDS",
  "Human_MASH_iPSC_derived_v1.4.2B.3.3.1",
  "analysis",
  "v1.4.2B.3.3.1_IL10_probe_integrated_audit",
  "tables",
  "IL10_cell_level_all_RDS_v3.3_with_strict_harmonized.tsv"
)

OUT_BASE <- file.path(
  BASE,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3A.1"
)

OBJ_DIR <- file.path(
  OUT_BASE,
  "objects"
)

AUDIT_DIR <- file.path(
  OUT_BASE,
  "analysis",
  "v1.4.3A.1_context_freeze"
)

OUT_RDS <- file.path(
  OBJ_DIR,
  "Human_MASH_12sample_Context_v1.4.3A.1_frozen.rds"
)

dir.create(
  OBJ_DIR,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  AUDIT_DIR,
  recursive=TRUE,
  showWarnings=FALSE
)

EXPECTED_SAMPLES <- c(
  paste0("LT", 1:8),
  "SfLB1",
  "SfLB2",
  "Tx17",
  "Tx5"
)

EXPECTED_TOTAL <- 124763L
EXPECTED_LT_TOTAL <- 81348L
EXPECTED_IPSC_TOTAL <- 43415L

EXPECTED_COUNTS <- c(
  LT1=10847L,
  LT2=12227L,
  LT3=7906L,
  LT4=9899L,
  LT5=12172L,
  LT6=10663L,
  LT7=7734L,
  LT8=9900L,
  SfLB1=19081L,
  SfLB2=9805L,
  Tx17=5958L,
  Tx5=8571L
)

# These values were independently established before the flawed v1.4.3A freeze.
# They are used here as hard validation targets for direct RNA-assay extraction.
EXPECTED_GENE_IL10_UMI <- c(
  LT1=83,
  LT2=180,
  LT3=236,
  LT4=62,
  LT5=126,
  LT6=199,
  LT7=46,
  LT8=113,
  SfLB1=0,
  SfLB2=0,
  Tx17=18,
  Tx5=4
)

EXPECTED_GENE_IL10_POSITIVE_CELLS <- c(
  LT1=63L,
  LT2=145L,
  LT3=133L,
  LT4=47L,
  LT5=101L,
  LT6=118L,
  LT7=42L,
  LT8=76L,
  SfLB1=0L,
  SfLB2=0L,
  Tx17=13L,
  Tx5=4L
)

message("============================================================")
message("Human MASH ", BIOLOGICAL_VERSION, " 12-sample Context Freeze correction")
message("============================================================")
message("Source RDS:    ", SOURCE_RDS)
message("Context meta:  ", CONTEXT_META)
message("IL10 table:    ", IL10_TABLE)
message("Frozen RDS:    ", OUT_RDS)
message("")
message("CORRECTION:")
message("- Gene-level IL10 is extracted directly from the RNA counts for ALL 12 samples.")
message("- LT1-LT8 probe-level audit fields are NA, not zero, because LT probe-level audit was not performed.")
message("- SfLB/Tx probe-level fields remain the validated v1.4.2B.3.3.1 values.")
message("")

required_files <- c(
  SOURCE_RDS,
  CONTEXT_META,
  IL10_TABLE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) missing:\n",
    paste(
      missing_files,
      collapse="\n"
    )
  )
}

normalize_sample <- function(x) {
  x <- as.character(x)
  x <- gsub("^SfLB_1$", "SfLB1", x)
  x <- gsub("^SfLB_2$", "SfLB2", x)
  x
}

to_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  y <- tolower(
    as.character(x)
  )

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in% c(
      "true",
      "t",
      "1",
      "yes"
    )
  ] <- TRUE

  out[
    y %in% c(
      "false",
      "f",
      "0",
      "no"
    )
  ] <- FALSE

  out
}

sum_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(
    x,
    na.rm=TRUE
  )
}

positive_count_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_integer_)
  }

  as.integer(
    sum(
      x > 0,
      na.rm=TRUE
    )
  )
}

extract_gene_counts_all_cells <- function(
  obj,
  assay_name="RNA",
  gene_candidates=c(
    "IL10",
    "ENSG00000136634"
  )
) {
  if (
    !assay_name %in%
    Assays(obj)
  ) {
    stop(
      "Assay '",
      assay_name,
      "' is not present. Available assays: ",
      paste(
        Assays(obj),
        collapse=", "
      )
    )
  }

  assay_obj <- obj[[assay_name]]
  all_cells <- colnames(obj)

  gene_hits <- intersect(
    gene_candidates,
    rownames(assay_obj)
  )

  if (length(gene_hits) == 0) {
    stop(
      "IL10 could not be found in RNA assay rownames. Candidates tried: ",
      paste(
        gene_candidates,
        collapse=", "
      )
    )
  }

  if ("IL10" %in% gene_hits) {
    gene_id <- "IL10"
  } else {
    gene_id <- gene_hits[1]
  }

  message(
    "Direct RNA-count extraction gene row: ",
    gene_id
  )

  layers <- tryCatch(
    Layers(
      assay_obj
    ),
    error=function(e) {
      character(0)
    }
  )

  count_layers <- layers[
    grepl(
      "^counts($|\\.)",
      layers
    )
  ]

  out <- rep(
    NA_real_,
    length(all_cells)
  )

  names(out) <- all_cells

  # Preferred Seurat v5 path: a single joined counts layer covering all cells.
  if (
    "counts" %in%
    count_layers
  ) {
    m <- LayerData(
      assay_obj,
      layer="counts"
    )

    if (
      gene_id %in%
      rownames(m) &&
      ncol(m) ==
      length(all_cells) &&
      setequal(
        colnames(m),
        all_cells
      )
    ) {
      idx <- match(
        all_cells,
        colnames(m)
      )

      vals <- as.numeric(
        m[
          gene_id,
          idx,
          drop=TRUE
        ]
      )

      names(vals) <- all_cells

      return(
        list(
          counts=vals,
          gene_id=gene_id,
          method="RNA LayerData(counts)",
          count_layers_used="counts"
        )
      )
    }
  }

  # Seurat v5 unjoined counts layers: require disjoint cell sets and full coverage.
  if (
    length(count_layers) > 0
  ) {
    seen <- rep(
      FALSE,
      length(all_cells)
    )

    names(seen) <- all_cells

    used_layers <- character(0)

    for (layer_name in count_layers) {
      m <- LayerData(
        assay_obj,
        layer=layer_name
      )

      if (
        !gene_id %in%
        rownames(m)
      ) {
        next
      }

      layer_cells <- colnames(m)

      if (
        is.null(layer_cells) ||
        length(layer_cells) == 0
      ) {
        stop(
          "Counts layer has no cell names: ",
          layer_name
        )
      }

      idx <- match(
        layer_cells,
        all_cells
      )

      if (anyNA(idx)) {
        stop(
          "Counts layer contains cells not present in the Seurat object: ",
          layer_name
        )
      }

      if (
        any(
          seen[idx]
        )
      ) {
        overlap_cells <- layer_cells[
          seen[idx]
        ]

        stop(
          "Overlapping cells detected across multiple counts layers. ",
          "Cannot safely infer a unique raw-count source. First overlaps: ",
          paste(
            head(
              overlap_cells,
              10
            ),
            collapse=", "
          )
        )
      }

      out[idx] <- as.numeric(
        m[
          gene_id,
          ,
          drop=TRUE
        ]
      )

      seen[idx] <- TRUE

      used_layers <- c(
        used_layers,
        layer_name
      )
    }

    if (
      length(used_layers) > 0 &&
      all(seen) &&
      !anyNA(out)
    ) {
      return(
        list(
          counts=out,
          gene_id=gene_id,
          method="RNA LayerData(disjoint counts layers)",
          count_layers_used=paste(
            used_layers,
            collapse=";"
          )
        )
      )
    }
  }

  # Seurat v4-compatible fallback.
  fallback <- tryCatch(
    GetAssayData(
      obj,
      assay=assay_name,
      slot="counts"
    ),
    error=function(e) {
      NULL
    }
  )

  if (!is.null(fallback)) {
    if (
      gene_id %in%
      rownames(fallback) &&
      ncol(fallback) ==
      length(all_cells) &&
      setequal(
        colnames(fallback),
        all_cells
      )
    ) {
      idx <- match(
        all_cells,
        colnames(fallback)
      )

      vals <- as.numeric(
        fallback[
          gene_id,
          idx,
          drop=TRUE
        ]
      )

      names(vals) <- all_cells

      return(
        list(
          counts=vals,
          gene_id=gene_id,
          method="GetAssayData(slot='counts')",
          count_layers_used="legacy counts"
        )
      )
    }
  }

  stop(
    "Could not extract a unique raw RNA-count vector for IL10 across all cells. ",
    "Available RNA layers: ",
    paste(
      layers,
      collapse=", "
    )
  )
}

message("Loading source 12-sample context RDS...")
obj <- readRDS(
  SOURCE_RDS
)

message("Loading validated v1.4.3.1 context metadata...")
context_meta <- read.delim(
  CONTEXT_META,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

message("Loading validated v1.4.2B.3.3.1 IL10 table...")
il10 <- read.delim(
  IL10_TABLE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

# -------------------------------------------------------------------------
# Source structure snapshot
# -------------------------------------------------------------------------

source_cells <- colnames(
  obj
)

source_assays <- Assays(
  obj
)

source_reductions <- Reductions(
  obj
)

source_idents <- as.character(
  Idents(
    obj
  )
)

if ("umapRPCA" %in% source_reductions) {
  frozen_umap_name <- "umapRPCA"
} else {
  umap_candidates <- source_reductions[
    grepl(
      "umap",
      source_reductions,
      ignore.case=TRUE
    )
  ]

  if (length(umap_candidates) == 0) {
    stop(
      "No pre-existing UMAP reduction found."
    )
  }

  frozen_umap_name <- umap_candidates[1]
}

source_umap <- Embeddings(
  obj,
  reduction=frozen_umap_name
)

# -------------------------------------------------------------------------
# Exact context-metadata reconciliation
# -------------------------------------------------------------------------

required_context_cols <- c(
  "cell",
  "sample",
  "source_group",
  "context_annotation",
  "context_state",
  "context_lineage",
  "annotation_source",
  "barcode_core"
)

missing_context_cols <- setdiff(
  required_context_cols,
  colnames(context_meta)
)

if (length(missing_context_cols) > 0) {
  stop(
    "Missing v1.4.3.1 context metadata column(s): ",
    paste(
      missing_context_cols,
      collapse=", "
    )
  )
}

if (
  anyDuplicated(
    context_meta$cell
  )
) {
  stop(
    "Duplicated cell IDs in v1.4.3.1 context metadata."
  )
}

j_context <- match(
  source_cells,
  context_meta$cell
)

if (anyNA(j_context)) {
  stop(
    "Context metadata does not map 100% to source RDS cells. Missing=",
    sum(
      is.na(
        j_context
      )
    )
  )
}

context_meta <- context_meta[
  j_context,
  ,
  drop=FALSE
]

if (
  !identical(
    as.character(
      context_meta$cell
    ),
    source_cells
  )
) {
  stop(
    "Context metadata order differs from source cell order."
  )
}

context_meta$sample <- normalize_sample(
  context_meta$sample
)

if (
  nrow(
    context_meta
  ) !=
  EXPECTED_TOTAL
) {
  stop(
    "Unexpected total cell count: ",
    nrow(context_meta),
    " expected ",
    EXPECTED_TOTAL
  )
}

if (
  !setequal(
    unique(
      context_meta$sample
    ),
    EXPECTED_SAMPLES
  )
) {
  stop(
    "Unexpected 12-sample set."
  )
}

for (S in EXPECTED_SAMPLES) {
  observed <- sum(
    context_meta$sample == S
  )

  expected <- EXPECTED_COUNTS[[S]]

  if (observed != expected) {
    stop(
      "Cell-count mismatch for ",
      S,
      ": observed=",
      observed,
      " expected=",
      expected
    )
  }
}

is_lt <- grepl(
  "^LT",
  context_meta$sample
)

is_ipsc <- context_meta$sample %in% c(
  "SfLB1",
  "SfLB2",
  "Tx17",
  "Tx5"
)

if (
  sum(is_lt) !=
  EXPECTED_LT_TOTAL
) {
  stop(
    "LT1-LT8 total mismatch."
  )
}

if (
  sum(is_ipsc) !=
  EXPECTED_IPSC_TOTAL
) {
  stop(
    "iPSC total mismatch."
  )
}

# -------------------------------------------------------------------------
# CORRECTION: direct RNA-assay gene-level IL10 for ALL 12 samples
# -------------------------------------------------------------------------

message("")
message("Extracting IL10 directly from RNA raw counts for all 12 samples...")

gene_result <- extract_gene_counts_all_cells(
  obj=obj,
  assay_name="RNA",
  gene_candidates=c(
    "IL10",
    "ENSG00000136634"
  )
)

gene_il10 <- gene_result$counts

if (
  length(
    gene_il10
  ) !=
  ncol(
    obj
  )
) {
  stop(
    "Direct RNA IL10 vector length mismatch."
  )
}

if (anyNA(gene_il10)) {
  stop(
    "NA detected in direct RNA IL10 count vector."
  )
}

if (
  any(
    gene_il10 < 0
  )
) {
  stop(
    "Negative IL10 raw count detected."
  )
}

direct_gene_summary <- do.call(
  rbind,
  lapply(
    EXPECTED_SAMPLES,
    function(S) {
      idx <- context_meta$sample == S

      data.frame(
        sample=S,
        n_cells=sum(idx),
        gene_IL10_total_UMI=sum(
          gene_il10[idx]
        ),
        gene_IL10_positive_cells=sum(
          gene_il10[idx] > 0
        ),
        max_gene_IL10_UMI=max(
          gene_il10[idx]
        ),
        stringsAsFactors=FALSE
      )
    }
  )
)

for (S in EXPECTED_SAMPLES) {
  observed_umi <- direct_gene_summary$gene_IL10_total_UMI[
    direct_gene_summary$sample == S
  ]

  expected_umi <- EXPECTED_GENE_IL10_UMI[[S]]

  if (
    length(observed_umi) != 1 ||
    observed_umi != expected_umi
  ) {
    stop(
      "Direct RNA IL10 UMI validation failed for ",
      S,
      ": observed=",
      observed_umi,
      " expected=",
      expected_umi
    )
  }

  observed_pos <- direct_gene_summary$gene_IL10_positive_cells[
    direct_gene_summary$sample == S
  ]

  expected_pos <- EXPECTED_GENE_IL10_POSITIVE_CELLS[[S]]

  if (
    length(observed_pos) != 1 ||
    observed_pos != expected_pos
  ) {
    stop(
      "Direct RNA IL10-positive cell validation failed for ",
      S,
      ": observed=",
      observed_pos,
      " expected=",
      expected_pos
    )
  }
}

write.table(
  direct_gene_summary,
  file=file.path(
    AUDIT_DIR,
    "v1.4.3A.1_direct_RNA_IL10_by_sample.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# -------------------------------------------------------------------------
# Validated iPSC probe-level audit reconciliation
# -------------------------------------------------------------------------

required_il10_cols <- c(
  "sample",
  "barcode_core",
  "gene_IL10_UMI_RDS",
  "accepted_IL10_UMI",
  "filteredout_IL10_UMI",
  "strict_harmonized_IL10_UMI",
  "gene_vs_accepted_concordant"
)

missing_il10_cols <- setdiff(
  required_il10_cols,
  colnames(il10)
)

if (length(missing_il10_cols) > 0) {
  stop(
    "Missing IL10 probe-audit column(s): ",
    paste(
      missing_il10_cols,
      collapse=", "
    )
  )
}

il10$sample <- normalize_sample(
  il10$sample
)

if (
  nrow(
    il10
  ) !=
  EXPECTED_IPSC_TOTAL
) {
  stop(
    "Unexpected iPSC IL10 audit row count: ",
    nrow(il10),
    " expected ",
    EXPECTED_IPSC_TOTAL
  )
}

il10_key <- paste(
  il10$sample,
  il10$barcode_core,
  sep="||"
)

if (
  anyDuplicated(
    il10_key
  )
) {
  stop(
    "Duplicated sample+barcode keys in IL10 probe-audit table."
  )
}

context_key <- paste(
  context_meta$sample,
  context_meta$barcode_core,
  sep="||"
)

j_il10 <- match(
  context_key,
  il10_key
)

if (
  any(
    is.na(
      j_il10[
        is_ipsc
      ]
    )
  )
) {
  stop(
    "iPSC IL10 probe-audit mapping incomplete. Missing=",
    sum(
      is.na(
        j_il10[
          is_ipsc
        ]
      )
    )
  )
}

if (
  any(
    !is.na(
      j_il10[
        is_lt
      ]
    )
  )
) {
  stop(
    "iPSC probe-audit table unexpectedly maps to LT cells."
  )
}

ipsc_map <- j_il10[
  is_ipsc
]

probe_gene_il10 <- as.numeric(
  il10$gene_IL10_UMI_RDS[
    ipsc_map
  ]
)

if (
  !identical(
    as.numeric(
      gene_il10[
        is_ipsc
      ]
    ),
    probe_gene_il10
  )
) {
  mismatch <- which(
    as.numeric(
      gene_il10[
        is_ipsc
      ]
    ) !=
    probe_gene_il10
  )

  stop(
    "Direct RNA IL10 does not match validated iPSC gene-level IL10 table cellwise. Mismatches=",
    length(mismatch)
  )
}

probe_concordance <- to_logical_safe(
  il10$gene_vs_accepted_concordant[
    ipsc_map
  ]
)

if (
  anyNA(
    probe_concordance
  ) ||
  !all(
    probe_concordance
  )
) {
  stop(
    "Validated iPSC gene-vs-accepted-probe concordance is not TRUE for all iPSC cells."
  )
}

# LT probe-level values are deliberately NA because no LT probe-level audit
# has been performed in this workflow.
accepted_il10 <- rep(
  NA_real_,
  ncol(obj)
)

filteredout_il10 <- rep(
  NA_real_,
  ncol(obj)
)

strict_il10 <- rep(
  NA_real_,
  ncol(obj)
)

gene_vs_accepted <- rep(
  NA,
  ncol(obj)
)

probe_audit_available <- rep(
  FALSE,
  ncol(obj)
)

accepted_il10[
  is_ipsc
] <- as.numeric(
  il10$accepted_IL10_UMI[
    ipsc_map
  ]
)

filteredout_il10[
  is_ipsc
] <- as.numeric(
  il10$filteredout_IL10_UMI[
    ipsc_map
  ]
)

strict_il10[
  is_ipsc
] <- as.numeric(
  il10$strict_harmonized_IL10_UMI[
    ipsc_map
  ]
)

gene_vs_accepted[
  is_ipsc
] <- probe_concordance

probe_audit_available[
  is_ipsc
] <- TRUE

accepted_positive <- rep(
  NA,
  ncol(obj)
)

strict_positive <- rep(
  NA,
  ncol(obj)
)

accepted_positive[
  is_ipsc
] <- accepted_il10[
    is_ipsc
  ] > 0

strict_positive[
  is_ipsc
] <- strict_il10[
    is_ipsc
  ] > 0

# Hard check that LT probe-level metadata are unknown, not falsely zero.
if (
  any(
    !is.na(
      accepted_il10[
        is_lt
      ]
    )
  ) ||
  any(
    !is.na(
      filteredout_il10[
        is_lt
      ]
    )
  ) ||
  any(
    !is.na(
      strict_il10[
        is_lt
      ]
    )
  ) ||
  any(
    !is.na(
      gene_vs_accepted[
        is_lt
      ]
    )
  )
) {
  stop(
    "LT probe-level metadata must remain NA."
  )
}

# iPSC validated probe totals
if (
  sum(
    accepted_il10[
      is_ipsc
    ],
    na.rm=TRUE
  ) !=
  22
) {
  stop(
    "Validated iPSC accepted IL10 total is not 22 UMI."
  )
}

if (
  sum(
    accepted_positive[
      is_ipsc
    ],
    na.rm=TRUE
  ) !=
  17
) {
  stop(
    "Validated iPSC accepted IL10-positive cell count is not 17."
  )
}

if (
  sum(
    strict_il10[
      is_ipsc
    ],
    na.rm=TRUE
  ) !=
  5
) {
  stop(
    "Validated strict b1bd8e2 total is not 5 UMI."
  )
}

if (
  sum(
    strict_positive[
      is_ipsc
    ],
    na.rm=TRUE
  ) !=
  4
) {
  stop(
    "Validated strict b1bd8e2-positive cell count is not 4."
  )
}

# -------------------------------------------------------------------------
# Add validated metadata only
# -------------------------------------------------------------------------

obj$sample_v143A1 <- context_meta$sample
obj$source_group_v143A1 <- context_meta$source_group
obj$barcode_core_v143A1 <- context_meta$barcode_core

obj$context_annotation_v143A1 <- context_meta$context_annotation
obj$context_state_v143A1 <- context_meta$context_state
obj$harmonized_lineage_v143A1 <- context_meta$context_lineage
obj$annotation_source_v143A1 <- context_meta$annotation_source

obj$gene_IL10_UMI_v143A1 <- gene_il10
obj$IL10_gene_positive_v143A1 <- gene_il10 > 0

obj$probe_IL10_audit_available_v143A1 <- probe_audit_available

obj$accepted_IL10_UMI_v143A1 <- accepted_il10
obj$filteredout_IL10_UMI_v143A1 <- filteredout_il10
obj$strict_b1bd8e2_IL10_UMI_v143A1 <- strict_il10

obj$IL10_accepted_positive_v143A1 <- accepted_positive
obj$IL10_strict_b1bd8e2_positive_v143A1 <- strict_positive
obj$IL10_gene_vs_accepted_concordant_v143A1 <- gene_vs_accepted

obj$context_source_RDS_v143A1 <- basename(
  SOURCE_RDS
)

obj$context_UMAP_reduction_v143A1 <- frozen_umap_name

obj$context_version_v143A1 <- BIOLOGICAL_VERSION
obj$context_freeze_date_v143A1 <- FREEZE_DATE
obj$context_frozen_v143A1 <- TRUE

obj$context_correction_note_v143A1 <- paste0(
  "Corrected from v1.4.3A: direct RNA gene-level IL10 is retained for all 12 samples; ",
  "LT1-LT8 probe-level IL10 fields are NA because LT probe-level audit was not performed."
)

# -------------------------------------------------------------------------
# In-memory immutability audit
# -------------------------------------------------------------------------

if (
  !identical(
    colnames(obj),
    source_cells
  )
) {
  stop(
    "Cell order changed during v1.4.3A.1 freeze."
  )
}

if (
  !identical(
    Assays(obj),
    source_assays
  )
) {
  stop(
    "Assay set changed during v1.4.3A.1 freeze."
  )
}

if (
  !identical(
    Reductions(obj),
    source_reductions
  )
) {
  stop(
    "Reduction set changed during v1.4.3A.1 freeze."
  )
}

if (
  !identical(
    as.character(
      Idents(obj)
    ),
    source_idents
  )
) {
  stop(
    "Active identities changed during v1.4.3A.1 freeze."
  )
}

if (
  !identical(
    Embeddings(
      obj,
      reduction=frozen_umap_name
    ),
    source_umap
  )
) {
  stop(
    "Existing UMAP coordinates changed during v1.4.3A.1 freeze."
  )
}

# -------------------------------------------------------------------------
# Audit tables
# -------------------------------------------------------------------------

metadata_dictionary <- data.frame(
  field=c(
    "sample_v143A1",
    "source_group_v143A1",
    "barcode_core_v143A1",
    "context_annotation_v143A1",
    "context_state_v143A1",
    "harmonized_lineage_v143A1",
    "annotation_source_v143A1",
    "gene_IL10_UMI_v143A1",
    "IL10_gene_positive_v143A1",
    "probe_IL10_audit_available_v143A1",
    "accepted_IL10_UMI_v143A1",
    "filteredout_IL10_UMI_v143A1",
    "strict_b1bd8e2_IL10_UMI_v143A1",
    "IL10_accepted_positive_v143A1",
    "IL10_strict_b1bd8e2_positive_v143A1",
    "IL10_gene_vs_accepted_concordant_v143A1",
    "context_source_RDS_v143A1",
    "context_UMAP_reduction_v143A1",
    "context_version_v143A1",
    "context_freeze_date_v143A1",
    "context_frozen_v143A1",
    "context_correction_note_v143A1"
  ),
  meaning=c(
    "Canonical 12-sample identity",
    "Native_LT / SfLB / Tx source group",
    "Normalized nucleotide barcode used for reconciliation",
    "Source-specific reconciled broad annotation",
    "Source-specific reconciled refined state",
    "Visualization-level harmonized lineage",
    "Annotation provenance",
    "Direct raw RNA-assay IL10 UMI for every one of the 12 samples",
    "Direct RNA-assay IL10-positive flag for every cell",
    "TRUE only for SfLB/Tx cells with validated probe-level audit; FALSE for LT1-LT8",
    "Sample-specific Cell Ranger accepted-probe IL10 UMI; NA for LT1-LT8",
    "Raw UMI from Cell Ranger-filtered-out IL10 probes; NA for LT1-LT8",
    "Shared accepted spliced probe b1bd8e2 UMI; NA for LT1-LT8",
    "Accepted-probe positive flag; NA for LT1-LT8",
    "Strict b1bd8e2 positive flag; NA for LT1-LT8",
    "Gene-level versus accepted-probe concordance; NA for LT1-LT8",
    "Source reconciled 12-sample RDS",
    "Existing UMAP reduction retained without recomputation",
    "Corrected context-freeze version",
    "Freeze date",
    "Frozen-context flag",
    "Explanation of v1.4.3A -> v1.4.3A.1 IL10 metadata correction"
  ),
  stringsAsFactors=FALSE
)

write.table(
  metadata_dictionary,
  file=file.path(
    AUDIT_DIR,
    "v1.4.3A.1_metadata_dictionary.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

sample_summary <- do.call(
  rbind,
  lapply(
    EXPECTED_SAMPLES,
    function(S) {
      idx <- obj$sample_v143A1 == S

      data.frame(
        sample=S,
        n_cells=sum(idx),
        gene_IL10_UMI=sum(
          obj$gene_IL10_UMI_v143A1[idx]
        ),
        gene_IL10_positive_cells=sum(
          obj$IL10_gene_positive_v143A1[idx]
        ),
        probe_audit_available=all(
          obj$probe_IL10_audit_available_v143A1[idx]
        ),
        accepted_IL10_UMI=sum_or_na(
          obj$accepted_IL10_UMI_v143A1[idx]
        ),
        accepted_IL10_positive_cells=positive_count_or_na(
          obj$accepted_IL10_UMI_v143A1[idx]
        ),
        filteredout_IL10_UMI=sum_or_na(
          obj$filteredout_IL10_UMI_v143A1[idx]
        ),
        strict_b1bd8e2_IL10_UMI=sum_or_na(
          obj$strict_b1bd8e2_IL10_UMI_v143A1[idx]
        ),
        strict_b1bd8e2_positive_cells=positive_count_or_na(
          obj$strict_b1bd8e2_IL10_UMI_v143A1[idx]
        ),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.table(
  sample_summary,
  file=file.path(
    AUDIT_DIR,
    "v1.4.3A.1_sample_summary.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

correction_audit <- data.frame(
  item=c(
    "v1.4.3A_problem",
    "v1.4.3A.1_gene_level_fix",
    "v1.4.3A.1_probe_level_semantics",
    "direct_RNA_extraction_method",
    "direct_RNA_gene_row",
    "iPSC_gene_cellwise_validation",
    "source_structure",
    "UMAP"
  ),
  result=c(
    "v1.4.3A incorrectly initialized LT1-LT8 IL10 metadata to zero because only iPSC probe-audit rows were overlaid",
    "IL10 gene-level raw counts extracted directly from RNA counts for all 124,763 cells",
    "LT1-LT8 accepted/filtered-out/strict/gene-vs-probe fields are NA because LT probe-level audit was not performed",
    gene_result$method,
    gene_result$gene_id,
    "Direct RNA gene-level IL10 exactly matches validated iPSC gene-level IL10 cellwise",
    "Assays, reductions, identities, cell order unchanged",
    paste0(frozen_umap_name, " coordinates exactly unchanged")
  ),
  stringsAsFactors=FALSE
)

write.table(
  correction_audit,
  file=file.path(
    AUDIT_DIR,
    "v1.4.3A.1_correction_audit.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# -------------------------------------------------------------------------
# Save corrected frozen RDS
# -------------------------------------------------------------------------

message("")
message("Saving corrected frozen RDS...")
message("This is a large Seurat object and may take several minutes.")

saveRDS(
  obj,
  OUT_RDS,
  compress=FALSE
)

if (
  !file.exists(
    OUT_RDS
  )
) {
  stop(
    "Corrected frozen RDS was not created."
  )
}

# -------------------------------------------------------------------------
# Reload and persisted-object validation
# -------------------------------------------------------------------------

message("")
message("Reloading corrected frozen RDS for persisted-object validation...")

check <- readRDS(
  OUT_RDS
)

if (
  ncol(
    check
  ) !=
  EXPECTED_TOTAL
) {
  stop(
    "Reloaded frozen RDS cell count mismatch."
  )
}

if (
  !identical(
    colnames(check),
    source_cells
  )
) {
  stop(
    "Reloaded frozen RDS cell order mismatch."
  )
}

if (
  !identical(
    Assays(check),
    source_assays
  )
) {
  stop(
    "Reloaded frozen RDS assay set mismatch."
  )
}

if (
  !identical(
    Reductions(check),
    source_reductions
  )
) {
  stop(
    "Reloaded frozen RDS reduction set mismatch."
  )
}

if (
  !identical(
    as.character(
      Idents(check)
    ),
    source_idents
  )
) {
  stop(
    "Reloaded frozen RDS identities mismatch."
  )
}

if (
  !identical(
    Embeddings(
      check,
      reduction=frozen_umap_name
    ),
    source_umap
  )
) {
  stop(
    "Reloaded frozen RDS UMAP coordinates changed."
  )
}

# Recheck all 12 sample gene-level IL10 totals from persisted metadata.
for (S in EXPECTED_SAMPLES) {
  idx <- check$sample_v143A1 == S

  observed_umi <- sum(
    check$gene_IL10_UMI_v143A1[idx]
  )

  expected_umi <- EXPECTED_GENE_IL10_UMI[[S]]

  if (
    observed_umi !=
    expected_umi
  ) {
    stop(
      "Reloaded gene-level IL10 mismatch for ",
      S
    )
  }

  observed_pos <- sum(
    check$IL10_gene_positive_v143A1[idx]
  )

  expected_pos <- EXPECTED_GENE_IL10_POSITIVE_CELLS[[S]]

  if (
    observed_pos !=
    expected_pos
  ) {
    stop(
      "Reloaded IL10-positive cell mismatch for ",
      S
    )
  }
}

lt_idx_check <- grepl(
  "^LT",
  check$sample_v143A1
)

ipsc_idx_check <- check$sample_v143A1 %in% c(
  "SfLB1",
  "SfLB2",
  "Tx17",
  "Tx5"
)

if (
  any(
    !is.na(
      check$accepted_IL10_UMI_v143A1[
        lt_idx_check
      ]
    )
  )
) {
  stop(
    "Reloaded LT accepted-probe IL10 must be NA."
  )
}

if (
  any(
    check$probe_IL10_audit_available_v143A1[
      lt_idx_check
    ]
  )
) {
  stop(
    "Reloaded LT probe audit availability must be FALSE."
  )
}

if (
  !all(
    check$probe_IL10_audit_available_v143A1[
      ipsc_idx_check
    ]
  )
) {
  stop(
    "Reloaded iPSC probe audit availability must be TRUE."
  )
}

if (
  sum(
    check$accepted_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ) !=
  22
) {
  stop(
    "Reloaded iPSC accepted IL10 total mismatch."
  )
}

if (
  sum(
    check$strict_b1bd8e2_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ) !=
  5
) {
  stop(
    "Reloaded iPSC strict IL10 total mismatch."
  )
}

# -------------------------------------------------------------------------
# Manifest / integrity
# -------------------------------------------------------------------------

rds_md5 <- unname(
  tools::md5sum(
    OUT_RDS
  )
)

rds_size <- file.info(
  OUT_RDS
)$size

manifest <- data.frame(
  biological_version=BIOLOGICAL_VERSION,
  freeze_date=FREEZE_DATE,
  source_RDS=SOURCE_RDS,
  frozen_RDS=OUT_RDS,
  frozen_RDS_size_bytes=rds_size,
  frozen_RDS_MD5=rds_md5,
  n_cells=ncol(check),
  n_genes=nrow(check),
  RNA_IL10_gene_row=gene_result$gene_id,
  RNA_IL10_extraction_method=gene_result$method,
  RNA_count_layers_used=gene_result$count_layers_used,
  gene_IL10_total_12sample=sum(
    check$gene_IL10_UMI_v143A1
  ),
  gene_IL10_positive_cells_12sample=sum(
    check$IL10_gene_positive_v143A1
  ),
  iPSC_accepted_IL10_total_UMI=sum(
    check$accepted_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  iPSC_accepted_IL10_positive_cells=sum(
    check$IL10_accepted_positive_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  iPSC_strict_b1bd8e2_total_UMI=sum(
    check$strict_b1bd8e2_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  iPSC_strict_b1bd8e2_positive_cells=sum(
    check$IL10_strict_b1bd8e2_positive_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  UMAP_reduction=frozen_umap_name,
  stringsAsFactors=FALSE
)

write.table(
  manifest,
  file=file.path(
    AUDIT_DIR,
    "v1.4.3A.1_freeze_manifest.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

summary_path <- file.path(
  AUDIT_DIR,
  "Human_MASH_v1.4.3A.1_12sample_ContextFreeze_summary.md"
)

con <- file(
  summary_path,
  open="wt"
)

writeLines(
  c(
    "# Human MASH v1.4.3A.1 12-sample Context Freeze",
    "",
    paste0("Freeze date: ", FREEZE_DATE),
    "",
    "## Correction from v1.4.3A",
    "",
    "- v1.4.3A incorrectly represented LT1-LT8 IL10 as zero in the newly added IL10 metadata because the metadata vector was initialized to zero and only the iPSC probe-audit table was overlaid.",
    "- v1.4.3A.1 extracts gene-level IL10 directly from the RNA raw-count matrix for all 124,763 cells.",
    "- LT1-LT8 probe-level accepted/filtered-out/strict fields are NA because LT probe-level audit was not performed.",
    "- SfLB/Tx probe-level fields retain the validated v1.4.2B.3.3.1 values.",
    "",
    "## Canonical corrected frozen object",
    "",
    paste0("- ", OUT_RDS),
    "",
    "## Direct RNA IL10 totals",
    ""
  ),
  con
)

for (i in seq_len(nrow(direct_gene_summary))) {
  x <- direct_gene_summary[i, ]

  writeLines(
    paste0(
      "- ",
      x$sample,
      ": ",
      x$gene_IL10_total_UMI,
      " UMI / ",
      x$gene_IL10_positive_cells,
      " positive cells"
    ),
    con
  )
}

writeLines(
  c(
    "",
    "## Probe-audit semantics",
    "",
    "- LT1-LT8: probe audit unavailable; probe-level metadata = NA.",
    "- SfLB1/SfLB2/Tx17/Tx5: validated probe audit available.",
    "- iPSC accepted-probe IL10 = 22 UMI / 17 cells.",
    "- iPSC strict shared spliced probe b1bd8e2 = 5 UMI / 4 cells.",
    "- Direct RNA IL10 matches the validated iPSC gene-level IL10 table cellwise.",
    "",
    "## Structural validation",
    "",
    paste0("- Total cells: ", ncol(check)),
    paste0("- Total genes: ", nrow(check)),
    paste0("- Existing UMAP retained exactly: ", frozen_umap_name),
    "- Cell order unchanged.",
    "- Assay set unchanged.",
    "- Reduction set unchanged.",
    "- Active identities unchanged.",
    "- UMAP coordinates unchanged after save/reload.",
    "",
    "## Intended canonical use",
    "",
    "- Use Human_MASH_12sample_Context_v1.4.3A.1_frozen.rds as the canonical 12-sample input for downstream Human MASH analyses.",
    "- Do not use the superseded v1.4.3A frozen RDS for IL10 metadata interpretation.",
    "",
    "## File integrity",
    "",
    paste0("- MD5: ", rds_md5),
    paste0("- Size bytes: ", rds_size),
    ""
  ),
  con
)

close(
  con
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con=file.path(
    AUDIT_DIR,
    "sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("v1.4.3A.1 CORRECTED CONTEXT FREEZE VALIDATION\n")
cat("============================================================\n")

cat(
  "Frozen RDS: ",
  OUT_RDS,
  "\n",
  sep=""
)

cat(
  "Total cells: ",
  ncol(check),
  "\n",
  sep=""
)

cat(
  "Total genes: ",
  nrow(check),
  "\n",
  sep=""
)

cat(
  "RNA IL10 extraction: ",
  gene_result$method,
  "\n",
  sep=""
)

cat(
  "RNA IL10 row: ",
  gene_result$gene_id,
  "\n",
  sep=""
)

cat(
  "12-sample gene IL10 total: ",
  sum(
    check$gene_IL10_UMI_v143A1
  ),
  " UMI / ",
  sum(
    check$IL10_gene_positive_v143A1
  ),
  " positive cells\n",
  sep=""
)

cat(
  "iPSC accepted IL10: ",
  sum(
    check$accepted_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  " UMI / ",
  sum(
    check$IL10_accepted_positive_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  " cells\n",
  sep=""
)

cat(
  "iPSC strict b1bd8e2: ",
  sum(
    check$strict_b1bd8e2_IL10_UMI_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  " UMI / ",
  sum(
    check$IL10_strict_b1bd8e2_positive_v143A1[
      ipsc_idx_check
    ],
    na.rm=TRUE
  ),
  " cells\n",
  sep=""
)

cat(
  "LT probe-level fields: NA by design\n"
)

cat(
  "UMAP reduction preserved: ",
  frozen_umap_name,
  "\n",
  sep=""
)

cat(
  "Frozen RDS size: ",
  format(
    rds_size,
    big.mark=","
  ),
  " bytes\n",
  sep=""
)

cat(
  "Frozen RDS MD5: ",
  rds_md5,
  "\n",
  sep=""
)

cat("\nCorrected sample summary:\n")

print(
  sample_summary,
  row.names=FALSE
)

cat("\n============================================================\n")
cat("PASS: DIRECT RNA IL10 VALIDATED FOR ALL 12 SAMPLES\n")
cat("PASS: LT PROBE-LEVEL FIELDS ARE NA, NOT FALSE ZERO\n")
cat("PASS: SOURCE STRUCTURE AND UMAP UNCHANGED AFTER SAVE/RELOAD\n")
cat("PASS: v1.4.3A.1 CORRECTED 12-SAMPLE CONTEXT FREEZE COMPLETE\n")
cat("============================================================\n")
