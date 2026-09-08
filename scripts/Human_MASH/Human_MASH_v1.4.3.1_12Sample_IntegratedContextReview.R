#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

VERSION <- "v1.4.3.1"
DATESTAMP <- "2026-09-08"

CONTEXT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.0/objects/",
  "Human_MASH_12sample_Context_v1.4.0_QC2_reconciled.rds"
)

LT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3.3/objects/",
  "Human_MASH_LT1to8_Canonical_v1.3.3.rds"
)

IPSC_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.1/objects/",
  "Human_MASH_iPSC_derived_v1.4.2B.1_mesenchymal_IL10_source_audit.rds"
)

IL10_TABLE <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.3.1/analysis/",
  "v1.4.2B.3.3.1_IL10_probe_integrated_audit/tables/",
  "IL10_cell_level_all_RDS_v3.3_with_strict_harmonized.tsv"
)

OUT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3.1/analysis/",
  "v1.4.3.1_12sample_integrated_context_review"
)

FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")

dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

EXPECTED_SAMPLES <- c(
  paste0("LT", 1:8),
  "SfLB1", "SfLB2",
  "Tx17", "Tx5"
)

IPSC_SAMPLES <- c(
  "SfLB1", "SfLB2",
  "Tx17", "Tx5"
)

EXPECTED_IPSC_COUNTS <- c(
  SfLB1=19081L,
  SfLB2=9805L,
  Tx17=5958L,
  Tx5=8571L
)

EXPECTED_TOTAL <- 124763L
EXPECTED_LT_TOTAL <- 81348L
EXPECTED_IPSC_TOTAL <- 43415L

SAMPLE_COLORS <- c(
  LT1="#0047FF",
  LT2="#0088FF",
  LT3="#00B8D4",
  LT4="#00C853",
  LT5="#64DD17",
  LT6="#FFD600",
  LT7="#FF9100",
  LT8="#FF3D00",
  SfLB1="#7C4DFF",
  SfLB2="#AA00FF",
  Tx17="#FF1493",
  Tx5="#C2185B"
)

SOURCE_COLORS <- c(
  Native_LT="#3F51B5",
  SfLB="#00B8D4",
  Tx="#FF1493"
)

LINEAGE_COLORS <- c(
  Hepatocyte="#008C95",
  Hepatobiliary="#00C853",
  Mesenchymal="#FF1493",
  Endothelial="#00B8D4",
  Myeloid="#FF6D00",
  Lymphoid="#7C4DFF",
  Neutrophil="#FFD600",
  Mast="#C2185B",
  Erythroid_Platelet="#D50000",
  Cycling="#795548",
  Unresolved="#616161",
  Other="#9E9E9E"
)

IL10_LINEAGE_COLORS <- c(
  Hepatocyte_like="#008C95",
  Hepatobiliary_like="#00C853",
  Mesenchymal_like="#FF1493",
  Endothelial_like="#00B8D4",
  Myeloid_like="#FF6D00",
  Unresolved="#616161"
)

message("============================================================")
message("Human MASH ", VERSION, " 12-sample integrated context review")
message("============================================================")
message("Context RDS: ", CONTEXT_RDS)
message("LT RDS:      ", LT_RDS)
message("iPSC RDS:    ", IPSC_RDS)
message("IL10 table:  ", IL10_TABLE)
message("Output:      ", OUT)
message("")

required_files <- c(
  CONTEXT_RDS,
  LT_RDS,
  IPSC_RDS,
  IL10_TABLE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) missing:\n",
    paste(missing_files, collapse="\n")
  )
}

normalize_sample <- function(x) {
  x <- as.character(x)
  x <- gsub("^SfLB_1$", "SfLB1", x)
  x <- gsub("^SfLB_2$", "SfLB2", x)
  x <- gsub("^SfLB[_-]?1$", "SfLB1", x)
  x <- gsub("^SfLB[_-]?2$", "SfLB2", x)
  x
}

extract_barcode <- function(x) {
  vapply(
    as.character(x),
    function(s) {
      if (is.na(s) || !nzchar(s)) {
        return(NA_character_)
      }

      m <- gregexpr(
        "[ACGT]{16,40}",
        s,
        perl=TRUE
      )[[1]]

      if (
        length(m) == 1 &&
        m[1] == -1
      ) {
        return(NA_character_)
      }

      ml <- attr(
        m,
        "match.length"
      )

      hits <- substring(
        s,
        first=m,
        last=m + ml - 1
      )

      hits[
        which.max(
          nchar(hits)
        )
      ]
    },
    character(1)
  )
}

sample_from_cellname <- function(x) {
  vapply(
    as.character(x),
    function(s) {
      hit <- regmatches(
        s,
        regexpr(
          "LT[1-8]|SfLB[_-]?[12]|Tx(?:17|5)",
          s,
          perl=TRUE
        )
      )

      if (
        length(hit) == 0 ||
        !nzchar(hit)
      ) {
        return(NA_character_)
      }

      normalize_sample(hit)
    },
    character(1)
  )
}

resolve_sample_column <- function(md, cell_names) {
  candidates <- c(
    "sample_v141",
    "sample",
    "orig.ident",
    "Sample",
    "sample_id",
    "sample_name"
  )

  best <- NULL
  best_score <- -Inf

  for (nm in candidates) {
    if (!nm %in% colnames(md)) {
      next
    }

    z <- normalize_sample(md[[nm]])

    score <- mean(
      z %in% EXPECTED_SAMPLES,
      na.rm=TRUE
    )

    if (
      is.finite(score) &&
      score > best_score
    ) {
      best <- z
      best_score <- score
    }
  }

  fallback <- sample_from_cellname(
    cell_names
  )

  if (is.null(best)) {
    best <- fallback
    best_score <- mean(
      best %in% EXPECTED_SAMPLES,
      na.rm=TRUE
    )
  } else {
    bad <- (
      is.na(best) |
      !best %in% EXPECTED_SAMPLES
    )

    best[bad] <- fallback[bad]
  }

  final_score <- mean(
    best %in% EXPECTED_SAMPLES,
    na.rm=TRUE
  )

  if (
    !is.finite(final_score) ||
    final_score < 0.99
  ) {
    stop(
      "Could not resolve sample identity for >=99% of cells. Coverage=",
      sprintf("%.4f", final_score)
    )
  }

  best
}

harmonize_lineage <- function(x) {
  x <- as.character(x)
  out <- rep(
    "Other",
    length(x)
  )

  out[
    grepl(
      "Hepatocyte",
      x,
      ignore.case=TRUE
    )
  ] <- "Hepatocyte"

  out[
    grepl(
      "Cholang|Duct|Hepatobiliary",
      x,
      ignore.case=TRUE
    )
  ] <- "Hepatobiliary"

  out[
    grepl(
      "HSC|Stellate|Fibro|Mesench|Pericy|Mesothel",
      x,
      ignore.case=TRUE
    )
  ] <- "Mesenchymal"

  out[
    grepl(
      "Endothel|LSEC|Vascular|Lymphatic",
      x,
      ignore.case=TRUE
    )
  ] <- "Endothelial"

  out[
    grepl(
      "Kupffer|Macroph|Monocyte|Myeloid|Dendritic",
      x,
      ignore.case=TRUE
    )
  ] <- "Myeloid"

  out[
    grepl(
      "T.?cell|NK|B.?cell|Plasma|Lymphoid",
      x,
      ignore.case=TRUE
    )
  ] <- "Lymphoid"

  out[
    grepl(
      "Neutro|Granulo",
      x,
      ignore.case=TRUE
    )
  ] <- "Neutrophil"

  out[
    grepl(
      "Mast",
      x,
      ignore.case=TRUE
    )
  ] <- "Mast"

  out[
    grepl(
      "Eryth|Platelet|Megak",
      x,
      ignore.case=TRUE
    )
  ] <- "Erythroid_Platelet"

  out[
    grepl(
      "Cycling",
      x,
      ignore.case=TRUE
    )
  ] <- "Cycling"

  out[
    grepl(
      "^Unresolved$|Ambiguous",
      x,
      ignore.case=TRUE
    )
  ] <- "Unresolved"

  out
}

save_plot_both <- function(
  p,
  stem,
  width,
  height
) {
  ggsave(
    file.path(
      FIG,
      paste0(
        stem,
        ".pdf"
      )
    ),
    p,
    width=width,
    height=height,
    limitsize=FALSE
  )

  ggsave(
    file.path(
      FIG,
      paste0(
        stem,
        ".png"
      )
    ),
    p,
    width=width,
    height=height,
    dpi=300,
    limitsize=FALSE
  )
}

message("Loading 12-sample reconciled context object...")
ctx <- readRDS(
  CONTEXT_RDS
)

message("Loading frozen LT annotation...")
lt <- readRDS(
  LT_RDS
)

message("Loading frozen iPSC lineage/state annotation...")
ipsc <- readRDS(
  IPSC_RDS
)

message("Loading validated IL10 probe table...")
il10 <- read.delim(
  IL10_TABLE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

ctx_md <- ctx@meta.data
lt_md <- lt@meta.data
ipsc_md <- ipsc@meta.data

ctx_sample <- resolve_sample_column(
  ctx_md,
  colnames(ctx)
)

lt_sample <- resolve_sample_column(
  lt_md,
  colnames(lt)
)

ipsc_sample <- resolve_sample_column(
  ipsc_md,
  colnames(ipsc)
)

ctx_barcode <- if (
  "original_cell_barcode" %in%
  colnames(ctx_md)
) {
  extract_barcode(
    ctx_md$original_cell_barcode
  )
} else {
  extract_barcode(
    colnames(ctx)
  )
}

lt_barcode <- if (
  "original_cell_barcode" %in%
  colnames(lt_md)
) {
  extract_barcode(
    lt_md$original_cell_barcode
  )
} else {
  extract_barcode(
    colnames(lt)
  )
}

ipsc_barcode <- if (
  "original_cell_barcode" %in%
  colnames(ipsc_md)
) {
  extract_barcode(
    ipsc_md$original_cell_barcode
  )
} else {
  extract_barcode(
    colnames(ipsc)
  )
}

ctx_key <- paste(
  ctx_sample,
  ctx_barcode,
  sep="||"
)

lt_key <- paste(
  lt_sample,
  lt_barcode,
  sep="||"
)

ipsc_key <- paste(
  ipsc_sample,
  ipsc_barcode,
  sep="||"
)

if (anyDuplicated(ctx_key)) {
  stop(
    "Duplicated sample+barcode keys in context object."
  )
}

if (anyDuplicated(lt_key)) {
  stop(
    "Duplicated sample+barcode keys in LT annotation object."
  )
}

if (anyDuplicated(ipsc_key)) {
  stop(
    "Duplicated sample+barcode keys in iPSC annotation object."
  )
}

if (
  !"celltype_final_v133" %in%
  colnames(lt_md)
) {
  stop(
    "LT metadata column celltype_final_v133 is missing."
  )
}

if (
  !"cellstate_final_v133" %in%
  colnames(lt_md)
) {
  stop(
    "LT metadata column cellstate_final_v133 is missing."
  )
}

if (
  !"broad_lineage_manual_v142A" %in%
  colnames(ipsc_md)
) {
  stop(
    "iPSC metadata column broad_lineage_manual_v142A is missing."
  )
}

if (
  !"manual_state_v142A" %in%
  colnames(ipsc_md)
) {
  stop(
    "iPSC metadata column manual_state_v142A is missing."
  )
}

context_annotation <- rep(
  NA_character_,
  ncol(ctx)
)

context_state <- rep(
  NA_character_,
  ncol(ctx)
)

annotation_source <- rep(
  NA_character_,
  ncol(ctx)
)

is_lt <- ctx_sample %in% paste0(
  "LT",
  1:8
)

is_ipsc <- ctx_sample %in% IPSC_SAMPLES

j_lt <- match(
  ctx_key[is_lt],
  lt_key
)

if (anyNA(j_lt)) {
  unmatched <- data.frame(
    cell=colnames(ctx)[is_lt][is.na(j_lt)],
    sample=ctx_sample[is_lt][is.na(j_lt)],
    barcode=ctx_barcode[is_lt][is.na(j_lt)],
    stringsAsFactors=FALSE
  )

  write.table(
    unmatched,
    file=file.path(
      TAB,
      "ERROR_unmatched_LT_cells.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  stop(
    "LT annotation mapping incomplete: ",
    sum(is.na(j_lt)),
    " cells unmatched."
  )
}

context_annotation[is_lt] <-
  as.character(
    lt_md$celltype_final_v133[
      j_lt
    ]
  )

context_state[is_lt] <-
  as.character(
    lt_md$cellstate_final_v133[
      j_lt
    ]
  )

annotation_source[is_lt] <-
  "LT_v1.3.3_frozen"

j_ipsc <- match(
  ctx_key[is_ipsc],
  ipsc_key
)

if (anyNA(j_ipsc)) {
  unmatched <- data.frame(
    cell=colnames(ctx)[is_ipsc][is.na(j_ipsc)],
    sample=ctx_sample[is_ipsc][is.na(j_ipsc)],
    barcode=ctx_barcode[is_ipsc][is.na(j_ipsc)],
    stringsAsFactors=FALSE
  )

  write.table(
    unmatched,
    file=file.path(
      TAB,
      "ERROR_unmatched_iPSC_cells.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  stop(
    "iPSC annotation mapping incomplete: ",
    sum(is.na(j_ipsc)),
    " cells unmatched."
  )
}

context_annotation[is_ipsc] <-
  as.character(
    ipsc_md$broad_lineage_manual_v142A[
      j_ipsc
    ]
  )

context_state[is_ipsc] <-
  as.character(
    ipsc_md$manual_state_v142A[
      j_ipsc
    ]
  )

annotation_source[is_ipsc] <-
  "iPSC_v1.4.2A_manual"

if (anyNA(context_annotation)) {
  stop(
    "Context annotation contains NA after reconciliation: ",
    sum(is.na(context_annotation))
  )
}

source_group <- ifelse(
  grepl(
    "^LT",
    ctx_sample
  ),
  "Native_LT",
  ifelse(
    grepl(
      "^SfLB",
      ctx_sample
    ),
    "SfLB",
    "Tx"
  )
)

context_lineage <- harmonize_lineage(
  context_annotation
)

reds <- Reductions(
  ctx
)

if (
  "umap" %in%
  reds
) {
  umap_name <- "umap"
} else {
  hit <- reds[
    grepl(
      "umap",
      reds,
      ignore.case=TRUE
    )
  ]

  if (length(hit) == 0) {
    stop(
      "No existing UMAP reduction found. ",
      "v1.4.3 must not recompute UMAP."
    )
  }

  umap_name <- hit[1]
}

message(
  "Using existing UMAP reduction: ",
  umap_name
)

emb <- Embeddings(
  ctx,
  reduction=umap_name
)

if (nrow(emb) != ncol(ctx)) {
  stop(
    "UMAP embedding row count does not match context cells."
  )
}

if (ncol(emb) < 2) {
  stop(
    "UMAP embedding has fewer than 2 dimensions."
  )
}

umap_df <- data.frame(
  cell=rownames(emb),
  UMAP_1=emb[, 1],
  UMAP_2=emb[, 2],
  sample=ctx_sample,
  source_group=source_group,
  context_annotation=context_annotation,
  context_state=context_state,
  context_lineage=context_lineage,
  annotation_source=annotation_source,
  barcode_core=ctx_barcode,
  stringsAsFactors=FALSE
)

umap_df$sample <- factor(
  umap_df$sample,
  levels=EXPECTED_SAMPLES
)

umap_df$source_group <- factor(
  umap_df$source_group,
  levels=c(
    "Native_LT",
    "SfLB",
    "Tx"
  )
)

lineage_levels <- names(
  LINEAGE_COLORS
)

umap_df$context_lineage <- factor(
  umap_df$context_lineage,
  levels=lineage_levels
)

sample_counts <- as.data.frame(
  table(
    umap_df$sample
  ),
  stringsAsFactors=FALSE
)

colnames(sample_counts) <- c(
  "sample",
  "n_cells"
)

sample_counts$sample <- as.character(
  sample_counts$sample
)

sample_counts <- sample_counts[
  sample_counts$sample %in%
  EXPECTED_SAMPLES,
  ,
  drop=FALSE
]

sample_counts <- sample_counts[
  match(
    EXPECTED_SAMPLES,
    sample_counts$sample
  ),
  ,
  drop=FALSE
]

if (
  nrow(umap_df) !=
  EXPECTED_TOTAL
) {
  stop(
    "Unexpected total context cells: ",
    nrow(umap_df),
    " expected ",
    EXPECTED_TOTAL
  )
}

if (
  sum(
    grepl(
      "^LT",
      as.character(
        umap_df$sample
      )
    )
  ) !=
  EXPECTED_LT_TOTAL
) {
  stop(
    "Unexpected LT1-LT8 total."
  )
}

if (
  sum(
    as.character(
      umap_df$sample
    ) %in%
    IPSC_SAMPLES
  ) !=
  EXPECTED_IPSC_TOTAL
) {
  stop(
    "Unexpected iPSC total."
  )
}

for (
  S in names(
    EXPECTED_IPSC_COUNTS
  )
) {
  observed <- sum(
    as.character(
      umap_df$sample
    ) == S
  )

  expected <- EXPECTED_IPSC_COUNTS[[S]]

  if (observed != expected) {
    stop(
      "Unexpected cell count for ",
      S,
      ": observed=",
      observed,
      " expected=",
      expected
    )
  }
}

write.table(
  sample_counts,
  file=file.path(
    TAB,
    "v1.4.3_sample_cell_counts.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

annotation_mapping_audit <- unique(
  data.frame(
    source_annotation=context_annotation,
    harmonized_lineage=as.character(
      context_lineage
    ),
    annotation_source=annotation_source,
    stringsAsFactors=FALSE
  )
)

annotation_mapping_audit <- annotation_mapping_audit[
  order(
    annotation_mapping_audit$annotation_source,
    annotation_mapping_audit$source_annotation
  ),
  ,
  drop=FALSE
]

write.table(
  annotation_mapping_audit,
  file=file.path(
    TAB,
    "v1.4.3_annotation_to_harmonized_lineage_audit.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

lineage_by_sample <- as.data.frame(
  table(
    sample=umap_df$sample,
    lineage=umap_df$context_lineage
  ),
  stringsAsFactors=FALSE
)

lineage_by_sample <- lineage_by_sample[
  lineage_by_sample$Freq > 0,
  ,
  drop=FALSE
]

colnames(lineage_by_sample)[3] <-
  "n_cells"

write.table(
  lineage_by_sample,
  file=file.path(
    TAB,
    "v1.4.3_harmonized_lineage_by_sample.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  umap_df,
  file=file.path(
    TAB,
    "v1.4.3_context_UMAP_metadata.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# ------------------------------------------------------------
# Map validated IL10 probe-audit table onto the 12-sample UMAP
# ------------------------------------------------------------

required_il10_cols <- c(
  "sample",
  "barcode_core",
  "gene_IL10_UMI_RDS",
  "accepted_IL10_UMI",
  "filteredout_IL10_UMI",
  "strict_harmonized_IL10_UMI",
  "lineage",
  "state"
)

missing_il10_cols <- setdiff(
  required_il10_cols,
  colnames(il10)
)

if (length(missing_il10_cols) > 0) {
  stop(
    "Missing IL10 table columns: ",
    paste(
      missing_il10_cols,
      collapse=", "
    )
  )
}

il10$sample <- normalize_sample(
  il10$sample
)

il10_key <- paste(
  il10$sample,
  il10$barcode_core,
  sep="||"
)

if (anyDuplicated(il10_key)) {
  stop(
    "Duplicated sample+barcode keys in IL10 table."
  )
}

j_il10 <- match(
  ctx_key,
  il10_key
)

umap_df$gene_IL10_UMI_RDS <- 0
umap_df$accepted_IL10_UMI <- 0
umap_df$filteredout_IL10_UMI <- 0
umap_df$strict_harmonized_IL10_UMI <- 0

ok_il10 <- !is.na(
  j_il10
)

umap_df$gene_IL10_UMI_RDS[
  ok_il10
] <- il10$gene_IL10_UMI_RDS[
  j_il10[
    ok_il10
  ]
]

umap_df$accepted_IL10_UMI[
  ok_il10
] <- il10$accepted_IL10_UMI[
  j_il10[
    ok_il10
  ]
]

umap_df$filteredout_IL10_UMI[
  ok_il10
] <- il10$filteredout_IL10_UMI[
  j_il10[
    ok_il10
  ]
]

umap_df$strict_harmonized_IL10_UMI[
  ok_il10
] <- il10$strict_harmonized_IL10_UMI[
  j_il10[
    ok_il10
  ]
]

accepted_cells <- umap_df[
  umap_df$accepted_IL10_UMI > 0,
  ,
  drop=FALSE
]

strict_cells <- umap_df[
  umap_df$strict_harmonized_IL10_UMI > 0,
  ,
  drop=FALSE
]

filteredout_only_cells <- umap_df[
  umap_df$filteredout_IL10_UMI > 0 &
  umap_df$accepted_IL10_UMI == 0,
  ,
  drop=FALSE
]

if (
  sum(
    accepted_cells$accepted_IL10_UMI
  ) !=
  22
) {
  stop(
    "Accepted IL10 overlay total is not 22 UMI."
  )
}

if (
  nrow(
    accepted_cells
  ) !=
  17
) {
  stop(
    "Accepted IL10 overlay cell count is not 17."
  )
}

if (
  sum(
    strict_cells$strict_harmonized_IL10_UMI
  ) !=
  5
) {
  stop(
    "Strict harmonized overlay total is not 5 UMI."
  )
}

if (
  nrow(
    strict_cells
  ) !=
  4
) {
  stop(
    "Strict harmonized overlay cell count is not 4."
  )
}

write.table(
  accepted_cells,
  file=file.path(
    TAB,
    "v1.4.3_IL10_accepted_positive_context_cells.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  strict_cells,
  file=file.path(
    TAB,
    "v1.4.3_IL10_strict_harmonized_positive_context_cells.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  filteredout_only_cells,
  file=file.path(
    TAB,
    "v1.4.3_IL10_filteredout_only_context_cells.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# ------------------------------------------------------------
# Figures
# Existing UMAP coordinates are preserved; no re-integration,
# PCA, neighbor graph, clustering, or UMAP recomputation occurs.
# ------------------------------------------------------------

theme_umap <- theme_void(
  base_size=11
) +
  theme(
    legend.position="right",
    plot.title=element_text(
      face="bold",
      size=13
    ),
    plot.subtitle=element_text(
      size=9
    ),
    strip.text=element_text(
      face="bold",
      size=9
    )
  )

p1 <- ggplot(
  umap_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=sample
  )
) +
  geom_point(
    size=0.18,
    alpha=0.78
  ) +
  scale_color_manual(
    values=SAMPLE_COLORS,
    drop=FALSE
  ) +
  coord_equal() +
  labs(
    title="12-sample integrated context UMAP",
    subtitle=paste0(
      "Existing ",
      umap_name,
      " coordinates preserved; LT1-LT8 + SfLB1/2 + Tx17/Tx5"
    ),
    color="Sample"
  ) +
  theme_umap

save_plot_both(
  p1,
  "Fig01_12sample_integrated_UMAP_by_sample",
  9.5,
  7.5
)

p2 <- ggplot(
  umap_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=source_group
  )
) +
  geom_point(
    size=0.18,
    alpha=0.8
  ) +
  scale_color_manual(
    values=SOURCE_COLORS,
    drop=FALSE
  ) +
  coord_equal() +
  labs(
    title="12-sample context UMAP by biological source",
    subtitle="Native LT versus non-transgene SfLB versus IL10-transgene-positive Tx graft",
    color="Source"
  ) +
  theme_umap

save_plot_both(
  p2,
  "Fig02_12sample_integrated_UMAP_by_source_group",
  9.5,
  7.5
)

p3 <- ggplot(
  umap_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=context_lineage
  )
) +
  geom_point(
    size=0.18,
    alpha=0.8
  ) +
  scale_color_manual(
    values=LINEAGE_COLORS,
    drop=FALSE
  ) +
  coord_equal() +
  labs(
    title="12-sample context UMAP by harmonized broad lineage",
    subtitle="LT frozen v1.3.3 annotations + iPSC manual v1.4.2A annotations; harmonized only for visualization",
    color="Broad lineage"
  ) +
  theme_umap

save_plot_both(
  p3,
  "Fig03_12sample_integrated_UMAP_harmonized_lineage",
  10.5,
  7.5
)

p4 <- ggplot(
  umap_df,
  aes(
    UMAP_1,
    UMAP_2,
    color=sample
  )
) +
  geom_point(
    size=0.18,
    alpha=0.85,
    show.legend=FALSE
  ) +
  scale_color_manual(
    values=SAMPLE_COLORS,
    drop=FALSE
  ) +
  facet_wrap(
    ~sample,
    ncol=4
  ) +
  coord_equal() +
  labs(
    title="12-sample split UMAP",
    subtitle="Same integrated coordinates and fixed axes in every panel"
  ) +
  theme_umap

save_plot_both(
  p4,
  "Fig04_12sample_split_UMAP",
  14,
  11
)

lt_bg <- umap_df[
  umap_df$source_group ==
  "Native_LT",
  ,
  drop=FALSE
]

ipsc_fg <- umap_df[
  umap_df$source_group !=
  "Native_LT",
  ,
  drop=FALSE
]

p5 <- ggplot() +
  geom_point(
    data=lt_bg,
    aes(
      UMAP_1,
      UMAP_2
    ),
    color="grey85",
    size=0.16,
    alpha=0.55
  ) +
  geom_point(
    data=ipsc_fg,
    aes(
      UMAP_1,
      UMAP_2,
      color=sample
    ),
    size=0.28,
    alpha=0.9
  ) +
  scale_color_manual(
    values=SAMPLE_COLORS[
      IPSC_SAMPLES
    ],
    drop=FALSE
  ) +
  coord_equal() +
  labs(
    title="iPSC-derived samples in native-liver context",
    subtitle="LT1-LT8 shown as grey reference background",
    color="iPSC sample"
  ) +
  theme_umap

save_plot_both(
  p5,
  "Fig05_iPSC_samples_highlight_over_LT_context",
  9.5,
  7.5
)

strongest_idx <- which.max(
  accepted_cells$accepted_IL10_UMI
)

strongest <- accepted_cells[
  strongest_idx,
  ,
  drop=FALSE
]

accepted_plot <- accepted_cells
accepted_plot$il10_lineage <- as.character(
  accepted_plot$context_annotation
)

accepted_plot$il10_lineage[
  !accepted_plot$il10_lineage %in%
  names(IL10_LINEAGE_COLORS)
] <- "Unresolved"

p6 <- ggplot() +
  geom_point(
    data=umap_df,
    aes(
      UMAP_1,
      UMAP_2
    ),
    color="grey88",
    size=0.14,
    alpha=0.5
  ) +
  geom_point(
    data=accepted_plot,
    aes(
      UMAP_1,
      UMAP_2,
      color=il10_lineage,
      size=accepted_IL10_UMI
    ),
    alpha=0.98
  ) +
  geom_point(
    data=strongest,
    aes(
      UMAP_1,
      UMAP_2
    ),
    shape=21,
    size=5,
    stroke=1.4,
    fill=NA,
    color="black"
  ) +
  scale_color_manual(
    values=IL10_LINEAGE_COLORS,
    drop=FALSE
  ) +
  scale_size_continuous(
    breaks=c(
      1,
      2,
      5
    ),
    range=c(
      1.8,
      5
    )
  ) +
  coord_equal() +
  labs(
    title="Cell Ranger-accepted IL10-positive cells in 12-sample context",
    subtitle=paste0(
      "17 Tx cells / 22 UMI; black ring = strongest cell (",
      strongest$context_state,
      ", ",
      strongest$accepted_IL10_UMI,
      " UMI)"
    ),
    color="Lineage/state class",
    size="Accepted IL10 UMI"
  ) +
  theme_umap

save_plot_both(
  p6,
  "Fig06_IL10_accepted_positive_overlay_12sample_context",
  10.5,
  7.5
)

strict_plot <- strict_cells
strict_plot$il10_lineage <- as.character(
  strict_plot$context_annotation
)

strict_plot$il10_lineage[
  !strict_plot$il10_lineage %in%
  names(IL10_LINEAGE_COLORS)
] <- "Unresolved"

strict_strongest <- strict_plot[
  which.max(
    strict_plot$strict_harmonized_IL10_UMI
  ),
  ,
  drop=FALSE
]

p7 <- ggplot() +
  geom_point(
    data=umap_df,
    aes(
      UMAP_1,
      UMAP_2
    ),
    color="grey88",
    size=0.14,
    alpha=0.5
  ) +
  geom_point(
    data=strict_plot,
    aes(
      UMAP_1,
      UMAP_2,
      color=il10_lineage,
      size=strict_harmonized_IL10_UMI
    ),
    alpha=0.98
  ) +
  geom_point(
    data=strict_strongest,
    aes(
      UMAP_1,
      UMAP_2
    ),
    shape=21,
    size=5,
    stroke=1.4,
    fill=NA,
    color="black"
  ) +
  scale_color_manual(
    values=IL10_LINEAGE_COLORS,
    drop=FALSE
  ) +
  scale_size_continuous(
    breaks=c(
      1,
      2
    ),
    range=c(
      2.2,
      5
    )
  ) +
  coord_equal() +
  labs(
    title="Strict harmonized IL10 signal in 12-sample context",
    subtitle=paste0(
      "Shared accepted spliced probe b1bd8e2 only: ",
      nrow(strict_plot),
      " Tx17 cells / ",
      sum(
        strict_plot$strict_harmonized_IL10_UMI
      ),
      " UMI"
    ),
    color="Lineage/state class",
    size="b1bd8e2 UMI"
  ) +
  theme_umap

save_plot_both(
  p7,
  "Fig07_IL10_strict_b1bd8e2_overlay_12sample_context",
  10.5,
  7.5
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_path <- file.path(
  OUT,
  "Human_MASH_v1.4.3.1_12sample_context_review_summary.md"
)

con <- file(
  summary_path,
  open="wt"
)

writeLines(
  c(
    "# Human MASH v1.4.3 12-sample integrated context review",
    "",
    paste0("Date: ", DATESTAMP),
    "",
    "## Scope",
    "",
    "- Review LT1-LT8, SfLB1/SfLB2, and Tx17/Tx5 in the pre-existing 12-sample RPCA/UMAP coordinate system.",
    "- No re-integration, PCA, neighbor graph, clustering, or UMAP recomputation was performed.",
    "- LT annotations were reconciled from the frozen v1.3.3 object.",
    "- iPSC annotations were reconciled from the frozen/manual v1.4.2A/B lineage-state object.",
    "- IL10 overlays were reconciled from the validated v1.4.2B.3.3.1 probe-audit table.",
    "",
    "## Input validation",
    "",
    paste0("- Total cells: ", nrow(umap_df)),
    paste0("- LT1-LT8 cells: ", sum(is_lt)),
    paste0("- iPSC-derived cells: ", sum(is_ipsc)),
    paste0("- Existing UMAP reduction used: ", umap_name),
    paste0("- LT annotation mapping: ", sum(is_lt), "/", sum(is_lt)),
    paste0("- iPSC annotation mapping: ", sum(is_ipsc), "/", sum(is_ipsc)),
    "",
    "## IL10 overlay validation",
    "",
    paste0(
      "- Cell Ranger-accepted IL10: ",
      nrow(accepted_cells),
      " cells / ",
      sum(accepted_cells$accepted_IL10_UMI),
      " UMI."
    ),
    paste0(
      "- Strict harmonized b1bd8e2 signal: ",
      nrow(strict_cells),
      " cells / ",
      sum(strict_cells$strict_harmonized_IL10_UMI),
      " UMI."
    ),
    paste0(
      "- Strongest accepted IL10 cell: ",
      strongest$cell,
      "; state=",
      strongest$context_state,
      "; ",
      strongest$accepted_IL10_UMI,
      " UMI."
    ),
    "",
    "## Interpretation boundary",
    "",
    "- UMAP is used only as a descriptive context map.",
    "- Geometric proximity on the integrated UMAP is not interpreted as lineage continuity, maturation distance, or quantitative biological similarity.",
    "- Sample-specific occupancy and isolated regions should be reviewed before downstream mechanistic interpretation.",
    "- IL10 overlays remain non-transgene-specific at the sequence level; they are interpreted as concordant with, not definitive proof of, transgene origin.",
    ""
  ),
  con
)

close(con)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con=file.path(
    OUT,
    "sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("v1.4.3.1 VALIDATION SUMMARY\n")
cat("============================================================\n")

cat(
  "Total cells: ",
  nrow(umap_df),
  "\n",
  sep=""
)

cat(
  "LT1-LT8 cells: ",
  sum(is_lt),
  "\n",
  sep=""
)

cat(
  "iPSC cells: ",
  sum(is_ipsc),
  "\n",
  sep=""
)

cat(
  "UMAP reduction preserved: ",
  umap_name,
  "\n",
  sep=""
)

cat(
  "Accepted IL10 overlay: ",
  nrow(accepted_cells),
  " cells / ",
  sum(
    accepted_cells$accepted_IL10_UMI
  ),
  " UMI\n",
  sep=""
)

cat(
  "Strict b1bd8e2 overlay: ",
  nrow(strict_cells),
  " cells / ",
  sum(
    strict_cells$strict_harmonized_IL10_UMI
  ),
  " UMI\n",
  sep=""
)

cat(
  "Strongest IL10 state: ",
  strongest$context_state,
  " | ",
  strongest$accepted_IL10_UMI,
  " UMI\n",
  sep=""
)

cat("\nSample counts:\n")
print(
  sample_counts,
  row.names=FALSE
)

cat("\n============================================================\n")
cat("OUTPUT FIGURES\n")
cat("============================================================\n")

print(
  list.files(
    FIG,
    pattern="\\.pdf$",
    full.names=TRUE
  )
)

cat("\n============================================================\n")
cat("v1.4.3.1 12-SAMPLE CONTEXT REVIEW COMPLETE\n")
cat("============================================================\n")
