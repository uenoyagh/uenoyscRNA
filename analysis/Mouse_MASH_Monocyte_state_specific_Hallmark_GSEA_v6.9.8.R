suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

VERSION <- "v6.9.8"

BASE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS"
)

INPUT_RDS <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.5",
  "objects",
  "Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds"
)

SPEC_FILE <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.6.4",
  "tables",
  "Monocyte_genomewide_support_aware_reference_specificity_v6.9.6.4.csv"
)

GLOBAL_GSEA_FILE <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.7.2",
  "tables",
  "Monocyte_Hallmark_GSEA_PRIMARY_v6.9.7.2.csv"
)

OUTDIR <- file.path(
  BASE_DIR,
  paste0("Mouse_MASH_Monocyte_", VERSION)
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte state-specific pseudobulk Hallmark GSEA\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing Monocyte RDS: ", INPUT_RDS)
}

if (!file.exists(SPEC_FILE)) {
  stop("Missing support-aware specificity table: ", SPEC_FILE)
}

if (!file.exists(GLOBAL_GSEA_FILE)) {
  stop("Missing global Hallmark GSEA file: ", GLOBAL_GSEA_FILE)
}

required_pkgs <- c(
  "fgsea",
  "msigdbr"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly=TRUE,
    FUN.VALUE=logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse=", ")
  )
}

# =========================================================
# Load object and metadata
# =========================================================

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

STATE_COL <- "Monocyte_state_frozen_v6.9.4"
CLUSTER_COL <- "Monocyte_cluster_frozen_v6.9.4"
CLASS_COL <- "Monocyte_analysis_class_v6.9.4"

required_meta <- c(
  STATE_COL,
  CLUSTER_COL,
  CLASS_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing_meta, collapse=", ")
  )
}

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
  stop("Could not resolve sample column.")
}

target_samples <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

sample_vec <- as.character(
  obj@meta.data[[sample_col]]
)

missing_samples <- setdiff(
  target_samples,
  unique(sample_vec)
)

if (length(missing_samples) > 0) {
  stop(
    "Missing Sham/Tx sample(s): ",
    paste(missing_samples, collapse=", ")
  )
}

# =========================================================
# Counts
# =========================================================

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

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

# =========================================================
# Determine eligible frozen states
#
# Primary rule:
# - analysis class must be primary or disease_enriched_primary
# - >=10 cells in EACH Sham/Tx biological sample
#
# This prevents cluster 7 and cluster 9 from being overinterpreted
# in state-specific pseudobulk when one or more replicates have
# too few cells.
# =========================================================

MIN_CELLS_PER_SAMPLE <- 10L

state_map <- unique(
  obj@meta.data[
    ,
    c(
      CLUSTER_COL,
      STATE_COL,
      CLASS_COL
    ),
    drop=FALSE
  ]
)

state_map <- state_map[
  order(
    as.numeric(
      as.character(
        state_map[[CLUSTER_COL]]
      )
    )
  ),
  ,
  drop=FALSE
]

state_count_rows <- list()
counter <- 1

for (i in seq_len(nrow(state_map))) {

  cl <- as.character(
    state_map[[CLUSTER_COL]][i]
  )

  state <- as.character(
    state_map[[STATE_COL]][i]
  )

  analysis_class <- as.character(
    state_map[[CLASS_COL]][i]
  )

  state_cells <- rownames(obj@meta.data)[
    as.character(
      obj@meta.data[[CLUSTER_COL]]
    ) == cl
  ]

  row <- data.frame(
    cluster=cl,
    state=state,
    analysis_class=analysis_class,
    stringsAsFactors=FALSE
  )

  for (s in target_samples) {

    row[[s]] <- sum(
      state_cells %in%
        rownames(obj@meta.data)[
          sample_vec == s
        ]
    )
  }

  row$min_cells_ShTx <- min(
    as.numeric(
      row[
        1,
        target_samples,
        drop=TRUE
      ]
    )
  )

  row$primary_class_eligible <-
    analysis_class %in% c(
      "primary",
      "disease_enriched_primary"
    )

  row$cell_support_eligible <-
    row$min_cells_ShTx >=
    MIN_CELLS_PER_SAMPLE

  row$state_specific_eligible <-
    row$primary_class_eligible &
    row$cell_support_eligible

  state_count_rows[[counter]] <- row
  counter <- counter + 1
}

state_support <- do.call(
  rbind,
  state_count_rows
)

rownames(state_support) <- NULL

write.csv(
  state_support,
  file.path(
    TABDIR,
    "Monocyte_state_specific_support_policy_v6.9.8.csv"
  ),
  row.names=FALSE
)

eligible_clusters <- state_support$cluster[
  state_support$state_specific_eligible
]

eligible_states <- state_support$state[
  state_support$state_specific_eligible
]

if (length(eligible_clusters) == 0) {
  stop("No states satisfy state-specific support criteria.")
}

cat(
  "Eligible clusters:",
  paste(eligible_clusters, collapse=", "),
  "\n"
)

cat(
  "Eligible states:",
  length(eligible_states),
  "\n\n"
)

# =========================================================
# Support-aware reference specificity
# =========================================================

spec <- read.csv(
  SPEC_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_spec_cols <- c(
  "gene",
  "keep_reference_aware_supported",
  "keep_strict_reference_supported"
)

missing_spec <- setdiff(
  required_spec_cols,
  colnames(spec)
)

if (length(missing_spec) > 0) {
  stop(
    "Missing specificity columns: ",
    paste(missing_spec, collapse=", ")
  )
}

# =========================================================
# Hallmark collection
#
# Use the SAME ortholog-mapped human Hallmark collection as v6.9.7.2
# so that global vs within-state results are directly comparable.
# =========================================================

hallmark_df <- tryCatch(
  {
    msigdbr::msigdbr(
      species="Mus musculus",
      collection="H"
    )
  },
  error=function(e) {
    msigdbr::msigdbr(
      species="Mus musculus",
      category="H"
    )
  }
)

if (!all(
  c(
    "gs_name",
    "gene_symbol"
  ) %in% colnames(hallmark_df)
)) {
  stop("Unexpected msigdbr Hallmark output.")
}

hallmark_df <- unique(
  hallmark_df[
    !is.na(hallmark_df$gene_symbol) &
    hallmark_df$gene_symbol != "",
    c(
      "gs_name",
      "gene_symbol"
    ),
    drop=FALSE
  ]
)

hallmark_pathways <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

cat(
  "Hallmark pathways loaded:",
  length(hallmark_pathways),
  "\n\n"
)

# =========================================================
# Helpers
# =========================================================

make_pseudobulk <- function(
  mat,
  cells,
  sample_vec,
  sample_order
) {

  pb <- matrix(
    0,
    nrow=nrow(mat),
    ncol=length(sample_order),
    dimnames=list(
      rownames(mat),
      sample_order
    )
  )

  cell_samples <- sample_vec[
    match(
      cells,
      colnames(obj)
    )
  ]

  for (s in sample_order) {

    scells <- cells[
      cell_samples == s
    ]

    if (length(scells) < MIN_CELLS_PER_SAMPLE) {
      stop(
        "Insufficient cells for sample ",
        s,
        ": ",
        length(scells)
      )
    }

    pb[, s] <- Matrix::rowSums(
      mat[
        ,
        scells,
        drop=FALSE
      ]
    )
  }

  pb
}

run_edger <- function(pb) {

  group <- factor(
    c(
      "Sham",
      "Sham",
      "Tx",
      "Tx"
    ),
    levels=c(
      "Sham",
      "Tx"
    )
  )

  y <- DGEList(
    counts=pb,
    group=group
  )

  keep <- filterByExpr(
    y,
    group=group,
    min.count=5
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(
    y,
    method="TMM"
  )

  design <- model.matrix(
    ~group
  )

  y <- estimateDisp(
    y,
    design=design,
    robust=TRUE
  )

  fit <- glmQLFit(
    y,
    design=design,
    robust=TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef="groupTx"
  )

  tt <- topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table

  tt$gene <- rownames(tt)
  rownames(tt) <- NULL

  tt[
    ,
    c(
      "gene",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]
}

build_rank <- function(
  de,
  spec,
  filter_col
) {

  idx <- match(
    de$gene,
    spec$gene
  )

  keep <- spec[[filter_col]][idx]
  keep[is.na(keep)] <- FALSE

  x <- de[
    keep &
    is.finite(de$logFC) &
    is.finite(de$F),
    ,
    drop=FALSE
  ]

  if (nrow(x) == 0) {
    return(NULL)
  }

  x$rank_metric <-
    sign(x$logFC) *
    sqrt(
      pmax(
        x$F,
        0
      )
    )

  x <- x[
    order(
      -abs(x$rank_metric)
    ),
    ,
    drop=FALSE
  ]

  x <- x[
    !duplicated(x$gene),
    ,
    drop=FALSE
  ]

  ranks <- x$rank_metric
  names(ranks) <- x$gene

  ranks <- sort(
    ranks,
    decreasing=TRUE
  )

  list(
    table=x,
    ranks=ranks
  )
}

run_gsea <- function(
  ranks,
  pathways,
  cluster,
  state,
  filter_name
) {

  g <- fgsea::fgseaMultilevel(
    pathways=pathways,
    stats=ranks,
    minSize=15,
    maxSize=500,
    eps=0
  )

  g <- as.data.frame(
    g,
    stringsAsFactors=FALSE
  )

  if (nrow(g) == 0) {
    return(g)
  }

  if ("leadingEdge" %in% colnames(g)) {

    g$leadingEdge <- vapply(
      g$leadingEdge,
      function(x) {
        paste(
          x,
          collapse=";"
        )
      },
      character(1)
    )
  }

  g$cluster <- cluster
  g$state <- state
  g$reference_filter <- filter_name

  g$direction <- ifelse(
    g$NES > 0,
    "Tx_enriched",
    "Sham_enriched"
  )

  g$significance_class <- ifelse(
    g$padj < 0.10,
    "FDR_lt0.10",
    ifelse(
      g$padj < 0.25,
      "exploratory_FDR_lt0.25",
      "not_significant"
    )
  )

  g[
    order(
      g$padj,
      -abs(g$NES)
    ),
    ,
    drop=FALSE
  ]
}

# =========================================================
# State-specific pseudobulk DE + Hallmark GSEA
# =========================================================

filter_defs <- c(
  REFERENCE_AWARE=
    "keep_reference_aware_supported",
  STRICT=
    "keep_strict_reference_supported"
)

de_results <- list()
gsea_results <- list()
rank_summary_rows <- list()
counter <- 1

for (cl in eligible_clusters) {

  state <- state_support$state[
    state_support$cluster == cl
  ][1]

  cat(
    "\nRunning cluster ",
    cl,
    ": ",
    state,
    "\n",
    sep=""
  )

  state_cells <- colnames(obj)[
    as.character(
      obj@meta.data[[CLUSTER_COL]]
    ) == cl &
    sample_vec %in% target_samples
  ]

  pb <- make_pseudobulk(
    mat=counts,
    cells=state_cells,
    sample_vec=sample_vec,
    sample_order=target_samples
  )

  write.csv(
    pb,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_cluster",
        cl,
        "_pseudobulk_counts_v6.9.8.csv"
      )
    )
  )

  de <- run_edger(pb)

  de$cluster <- cl
  de$state <- state

  de_results[[cl]] <- de

  write.csv(
    de,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_cluster",
        cl,
        "_Sham_vs_Tx_DE_v6.9.8.csv"
      )
    ),
    row.names=FALSE
  )

  for (filter_name in names(filter_defs)) {

    rank_obj <- build_rank(
      de=de,
      spec=spec,
      filter_col=filter_defs[[filter_name]]
    )

    if (is.null(rank_obj)) {
      next
    }

    write.csv(
      rank_obj$table,
      file.path(
        TABDIR,
        paste0(
          "Monocyte_cluster",
          cl,
          "_",
          filter_name,
          "_rank_v6.9.8.csv"
        )
      ),
      row.names=FALSE
    )

    rank_summary_rows[[counter]] <- data.frame(
      cluster=cl,
      state=state,
      reference_filter=filter_name,
      n_ranked_genes=length(rank_obj$ranks),
      max_rank=max(rank_obj$ranks),
      min_rank=min(rank_obj$ranks),
      stringsAsFactors=FALSE
    )

    counter <- counter + 1

    gsea <- run_gsea(
      ranks=rank_obj$ranks,
      pathways=hallmark_pathways,
      cluster=cl,
      state=state,
      filter_name=filter_name
    )

    key <- paste(
      cl,
      filter_name,
      sep="__"
    )

    gsea_results[[key]] <- gsea

    write.csv(
      gsea,
      file.path(
        TABDIR,
        paste0(
          "Monocyte_cluster",
          cl,
          "_Hallmark_GSEA_",
          filter_name,
          "_v6.9.8.csv"
        )
      ),
      row.names=FALSE
    )
  }
}

# =========================================================
# Combined outputs
# =========================================================

all_de <- do.call(
  rbind,
  de_results
)

rownames(all_de) <- NULL

write.csv(
  all_de,
  file.path(
    TABDIR,
    "Monocyte_state_specific_DE_all_eligible_states_v6.9.8.csv"
  ),
  row.names=FALSE
)

all_gsea <- do.call(
  rbind,
  gsea_results
)

rownames(all_gsea) <- NULL

write.csv(
  all_gsea,
  file.path(
    TABDIR,
    "Monocyte_state_specific_Hallmark_GSEA_all_v6.9.8.csv"
  ),
  row.names=FALSE
)

rank_summary <- do.call(
  rbind,
  rank_summary_rows
)

rownames(rank_summary) <- NULL

write.csv(
  rank_summary,
  file.path(
    TABDIR,
    "Monocyte_state_specific_GSEA_rank_summary_v6.9.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Pathway concordance across eligible states
# =========================================================

reference_gsea <- all_gsea[
  all_gsea$reference_filter ==
    "REFERENCE_AWARE",
  ,
  drop=FALSE
]

strict_gsea <- all_gsea[
  all_gsea$reference_filter ==
    "STRICT",
  ,
  drop=FALSE
]

summarize_pathways <- function(gsea_df, label) {

  pathways <- sort(
    unique(
      gsea_df$pathway
    )
  )

  rows <- list()

  for (pw in pathways) {

    x <- gsea_df[
      gsea_df$pathway == pw,
      ,
      drop=FALSE
    ]

    nes <- x$NES[
      is.finite(x$NES)
    ]

    rows[[pw]] <- data.frame(
      pathway=pw,
      reference_filter=label,
      n_states=nrow(x),
      n_Tx_enriched=sum(
        x$NES > 0,
        na.rm=TRUE
      ),
      n_Sham_enriched=sum(
        x$NES < 0,
        na.rm=TRUE
      ),
      n_Tx_FDRlt0.10=sum(
        x$NES > 0 &
        x$padj < 0.10,
        na.rm=TRUE
      ),
      n_Tx_FDRlt0.25=sum(
        x$NES > 0 &
        x$padj < 0.25,
        na.rm=TRUE
      ),
      n_Sham_FDRlt0.10=sum(
        x$NES < 0 &
        x$padj < 0.10,
        na.rm=TRUE
      ),
      n_Sham_FDRlt0.25=sum(
        x$NES < 0 &
        x$padj < 0.25,
        na.rm=TRUE
      ),
      median_NES=median(
        nes,
        na.rm=TRUE
      ),
      all_states_same_direction=
        length(nes) == length(eligible_clusters) &&
        (
          all(nes > 0) ||
          all(nes < 0)
        ),
      stringsAsFactors=FALSE
    )
  }

  out <- do.call(
    rbind,
    rows
  )

  rownames(out) <- NULL
  out
}

reference_summary <- summarize_pathways(
  reference_gsea,
  "REFERENCE_AWARE"
)

strict_summary <- summarize_pathways(
  strict_gsea,
  "STRICT"
)

pathway_summary <- rbind(
  reference_summary,
  strict_summary
)

write.csv(
  pathway_summary,
  file.path(
    TABDIR,
    "Monocyte_state_specific_Hallmark_concordance_v6.9.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Compare with global v6.9.7.2
# =========================================================

global_gsea <- read.csv(
  GLOBAL_GSEA_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

global_key <- global_gsea[
  global_gsea$padj < 0.25,
  c(
    "pathway",
    "NES",
    "padj",
    "direction"
  ),
  drop=FALSE
]

names(global_key)[
  names(global_key) == "NES"
] <- "global_PRIMARY_NES"

names(global_key)[
  names(global_key) == "padj"
] <- "global_PRIMARY_FDR"

names(global_key)[
  names(global_key) == "direction"
] <- "global_PRIMARY_direction"

global_state_concordance <- merge(
  global_key,
  reference_summary,
  by="pathway",
  all.x=TRUE
)

global_state_concordance <- global_state_concordance[
  order(
    global_state_concordance$global_PRIMARY_FDR
  ),
  ,
  drop=FALSE
]

write.csv(
  global_state_concordance,
  file.path(
    TABDIR,
    "Monocyte_global_vs_state_specific_Hallmark_v6.9.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Selected state-specific results
# =========================================================

selected_reference <- reference_gsea[
  reference_gsea$padj < 0.25,
  ,
  drop=FALSE
]

selected_strict <- strict_gsea[
  strict_gsea$padj < 0.25,
  ,
  drop=FALSE
]

write.csv(
  selected_reference,
  file.path(
    TABDIR,
    "Monocyte_state_specific_Hallmark_REFERENCE_AWARE_FDRlt0.25_v6.9.8.csv"
  ),
  row.names=FALSE
)

write.csv(
  selected_strict,
  file.path(
    TABDIR,
    "Monocyte_state_specific_Hallmark_STRICT_FDRlt0.25_v6.9.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Heatmap of global-primary pathways across eligible states
# =========================================================

heat_pathways <- global_key$pathway

heat_df <- reference_gsea[
  reference_gsea$pathway %in%
    heat_pathways,
  c(
    "cluster",
    "state",
    "pathway",
    "NES",
    "padj"
  ),
  drop=FALSE
]

if (nrow(heat_df) > 0) {

  heat_df$cluster <- factor(
    heat_df$cluster,
    levels=eligible_clusters
  )

  global_order <- global_key$pathway[
    order(
      global_key$global_PRIMARY_FDR
    )
  ]

  heat_df$pathway <- factor(
    heat_df$pathway,
    levels=rev(global_order)
  )

  p <- ggplot(
    heat_df,
    aes(
      x=cluster,
      y=pathway,
      fill=NES
    )
  ) +
    geom_tile() +
    geom_point(
      data=heat_df[
        heat_df$padj < 0.10,
        ,
        drop=FALSE
      ],
      shape=21,
      size=2.2
    ) +
    scale_fill_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(
      base_size=9
    ) +
    labs(
      title="Monocyte within-state Hallmark GSEA",
      subtitle="REFERENCE_AWARE; dots indicate within-state FDR < 0.10",
      x="Frozen Monocyte cluster",
      y=NULL,
      fill="NES"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Monocyte_state_specific_Hallmark_heatmap_v6.9.8.pdf"
    ),
    p,
    width=9,
    height=10
  )
}

# =========================================================
# Save a lightweight metadata checkpoint object
# =========================================================

obj$Monocyte_state_specific_GSEA_eligible_v6.9.8 <-
  as.character(
    obj@meta.data[[CLUSTER_COL]]
  ) %in%
  eligible_clusters

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_state_specific_GSEA_metadata_v6.9.8.rds"
  )
)

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== STATE-SPECIFIC SUPPORT POLICY ===\n")
print(
  state_support,
  row.names=FALSE
)

cat("\n=== STATE-SPECIFIC GSEA RANK SUMMARY ===\n")
print(
  rank_summary,
  row.names=FALSE
)

cat("\n=== GLOBAL VS STATE-SPECIFIC HALLMARK ===\n")

print_cols <- c(
  "pathway",
  "global_PRIMARY_NES",
  "global_PRIMARY_FDR",
  "n_states",
  "n_Tx_enriched",
  "n_Tx_FDRlt0.10",
  "n_Tx_FDRlt0.25",
  "n_Sham_enriched",
  "median_NES",
  "all_states_same_direction"
)

print(
  global_state_concordance[
    ,
    intersect(
      print_cols,
      colnames(global_state_concordance)
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== REFERENCE-AWARE PATHWAYS WITH >=4 TX-ENRICHED STATES ===\n")

robust_tx <- reference_summary[
  reference_summary$n_Tx_enriched >= 4,
  ,
  drop=FALSE
]

robust_tx <- robust_tx[
  order(
    -robust_tx$n_Tx_FDRlt0.10,
    -robust_tx$n_Tx_FDRlt0.25,
    -robust_tx$median_NES
  ),
  ,
  drop=FALSE
]

if (nrow(robust_tx) == 0) {
  cat("NONE\n")
} else {
  print(
    robust_tx,
    row.names=FALSE
  )
}

cat("\n=== STRICT PATHWAYS WITH >=4 TX-ENRICHED STATES ===\n")

robust_tx_strict <- strict_summary[
  strict_summary$n_Tx_enriched >= 4,
  ,
  drop=FALSE
]

robust_tx_strict <- robust_tx_strict[
  order(
    -robust_tx_strict$n_Tx_FDRlt0.10,
    -robust_tx_strict$n_Tx_FDRlt0.25,
    -robust_tx_strict$median_NES
  ),
  ,
  drop=FALSE
]

if (nrow(robust_tx_strict) == 0) {
  cat("NONE\n")
} else {
  print(
    robust_tx_strict,
    row.names=FALSE
  )
}

cat("\n====================================================\n")
cat("v6.9.8 COMPLETE\n")
cat("State-specific Monocyte pseudobulk Hallmark GSEA complete\n")
cat("Eligibility: primary class + >=10 cells in EACH Sham/Tx sample\n")
cat("Global comparison: v6.9.7.2 PRIMARY_CORE reference-aware\n")
cat("Reference filters: support-aware <=2 primary; <=1.5 strict sensitivity\n")
cat("Rank metric: sign(logFC) * sqrt(edgeR QL F)\n")
cat("Sham vs Tx biological n=2/group\n")
cat("No cells removed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.8.txt"
  )
)
