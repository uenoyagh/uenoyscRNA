suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v6.9.7.2"

BASE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS"
)

DE_DIR <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.6",
  "tables"
)

SPEC_DIR <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.6.4",
  "tables"
)

OUTDIR <- file.path(
  BASE_DIR,
  paste0("Mouse_MASH_Monocyte_", VERSION)
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte reference-aware Hallmark GSEA\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# =========================================================
# Package checks
# =========================================================

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
    paste(missing_pkgs, collapse=", "),
    "\nInstall them before running v6.9.7."
  )
}

# =========================================================
# Inputs
# =========================================================

SPEC_FILE <- file.path(
  SPEC_DIR,
  "Monocyte_genomewide_support_aware_reference_specificity_v6.9.6.4.csv"
)

if (!file.exists(SPEC_FILE)) {
  stop(
    "Missing support-aware specificity table: ",
    SPEC_FILE
  )
}

spec <- read.csv(
  SPEC_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

frameworks <- c(
  "ALL",
  "NO_QCWATCH",
  "PRIMARY_CORE"
)

de_files <- setNames(
  file.path(
    DE_DIR,
    paste0(
      "Monocyte_pseudobulk_",
      frameworks,
      "_Sham_vs_Tx_v6.9.6.csv"
    )
  ),
  frameworks
)

missing_de <- de_files[
  !file.exists(de_files)
]

if (length(missing_de) > 0) {
  stop(
    "Missing DE file(s): ",
    paste(missing_de, collapse=", ")
  )
}

# =========================================================
# Hallmark gene sets
#
# Prefer current msigdbr API; fall back to older API.
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

required_hallmark_cols <- c(
  "gs_name",
  "gene_symbol"
)

if (!all(
  required_hallmark_cols %in%
    colnames(hallmark_df)
)) {
  stop(
    "Unexpected msigdbr output; required columns not found: ",
    paste(required_hallmark_cols, collapse=", ")
  )
}

hallmark_df <- hallmark_df[
  !is.na(hallmark_df$gene_symbol) &
  hallmark_df$gene_symbol != "",
  c(
    "gs_name",
    "gene_symbol"
  ),
  drop=FALSE
]

hallmark_df <- unique(
  hallmark_df
)

hallmark_pathways <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

cat(
  "Hallmark pathways loaded:",
  length(hallmark_pathways),
  "\n"
)

# =========================================================
# Ranking method
#
# edgeR quasi-likelihood F statistic is unsigned.
# Rank = sign(logFC) * sqrt(F)
#
# This avoids the large number of ties that can occur with
# sign(logFC) * -log10(PValue), while preserving DE direction.
# =========================================================

build_rank <- function(
  de,
  spec,
  filter_col
) {

  required_de_cols <- c(
    "gene",
    "logFC",
    "F",
    "PValue"
  )

  missing_cols <- setdiff(
    required_de_cols,
    colnames(de)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing DE columns: ",
      paste(missing_cols, collapse=", ")
    )
  }

  if (!(filter_col %in% colnames(spec))) {
    stop(
      "Missing specificity filter column: ",
      filter_col
    )
  }

  idx <- match(
    de$gene,
    spec$gene
  )

  keep <- spec[[filter_col]][idx]

  keep[
    is.na(keep)
  ] <- FALSE

  x <- de[
    keep &
    is.finite(de$logFC) &
    is.finite(de$F),
    ,
    drop=FALSE
  ]

  if (nrow(x) == 0) {
    stop(
      "No genes remain after filter: ",
      filter_col
    )
  }

  x$rank_metric <-
    sign(x$logFC) *
    sqrt(
      pmax(
        x$F,
        0
      )
    )

  # One score per symbol.
  # If duplicated symbols exist, retain the largest absolute rank.
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

# =========================================================
# GSEA helper
# =========================================================

run_hallmark_gsea <- function(
  ranks,
  pathways,
  framework,
  filter_name
) {

  result <- fgsea::fgseaMultilevel(
    pathways=pathways,
    stats=ranks,
    minSize=15,
    maxSize=500,
    eps=0
  )

  result <- as.data.frame(
    result,
    stringsAsFactors=FALSE
  )

  if (nrow(result) == 0) {
    return(result)
  }

  # fgsea returns leadingEdge as a list-column.
  # Convert it to a semicolon-delimited character column so that
  # write.csv() can export the result without EncodeElement errors.
  if ("leadingEdge" %in% colnames(result)) {
    result$leadingEdge <- vapply(
      result$leadingEdge,
      function(x) paste(x, collapse=";"),
      character(1)
    )
  }

  result$framework <- framework
  result$reference_filter <- filter_name

  result$direction <- ifelse(
    result$NES > 0,
    "Tx_enriched",
    "Sham_enriched"
  )

  result$significance_class <- ifelse(
    result$padj < 0.10,
    "FDR_lt0.10",
    ifelse(
      result$padj < 0.25,
      "exploratory_FDR_lt0.25",
      "not_significant"
    )
  )

  result <- result[
    order(
      result$padj,
      -abs(result$NES)
    ),
    ,
    drop=FALSE
  ]

  result
}

# =========================================================
# Run six analyses:
# 3 frameworks x 2 reference filters
# =========================================================

filter_definitions <- c(
  REFERENCE_AWARE=
    "keep_reference_aware_supported",
  STRICT=
    "keep_strict_reference_supported"
)

all_results <- list()
rank_summary_rows <- list()
counter <- 1

for (fw in frameworks) {

  cat(
    "\nReading DE:",
    fw,
    "\n"
  )

  de <- read.csv(
    de_files[[fw]],
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  for (filter_name in names(filter_definitions)) {

    filter_col <- filter_definitions[[filter_name]]

    rank_obj <- build_rank(
      de=de,
      spec=spec,
      filter_col=filter_col
    )

    rank_table <- rank_obj$table
    ranks <- rank_obj$ranks

    rank_summary_rows[[counter]] <- data.frame(
      framework=fw,
      reference_filter=filter_name,
      n_ranked_genes=length(ranks),
      max_rank=max(ranks),
      min_rank=min(ranks),
      stringsAsFactors=FALSE
    )

    counter <- counter + 1

    write.csv(
      rank_table,
      file.path(
        TABDIR,
        paste0(
          "Monocyte_",
          fw,
          "_",
          filter_name,
          "_signed_sqrtF_rank_v6.9.7.2.csv"
        )
      ),
      row.names=FALSE
    )

    gsea <- run_hallmark_gsea(
      ranks=ranks,
      pathways=hallmark_pathways,
      framework=fw,
      filter_name=filter_name
    )

    key <- paste(
      fw,
      filter_name,
      sep="__"
    )

    all_results[[key]] <- gsea

    write.csv(
      gsea,
      file.path(
        TABDIR,
        paste0(
          "Monocyte_Hallmark_GSEA_",
          fw,
          "_",
          filter_name,
          "_v6.9.7.2.csv"
        )
      ),
      row.names=FALSE
    )
  }
}

rank_summary <- do.call(
  rbind,
  rank_summary_rows
)

write.csv(
  rank_summary,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_rank_summary_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Combined GSEA table
# =========================================================

combined <- do.call(
  rbind,
  all_results
)

rownames(combined) <- NULL

write.csv(
  combined,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_all_frameworks_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Primary analysis:
# PRIMARY_CORE + REFERENCE_AWARE
# =========================================================

primary <- all_results[["PRIMARY_CORE__REFERENCE_AWARE"]]

strict_primary <- all_results[["PRIMARY_CORE__STRICT"]]

write.csv(
  primary,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_PRIMARY_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

write.csv(
  strict_primary,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_PRIMARY_strict_sensitivity_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cross-analysis pathway concordance
# =========================================================

pathway_names <- sort(
  unique(
    combined$pathway
  )
)

concordance_rows <- list()

for (pw in pathway_names) {

  x <- combined[
    combined$pathway == pw,
    ,
    drop=FALSE
  ]

  primary_row <- x[
    x$framework == "PRIMARY_CORE" &
    x$reference_filter == "REFERENCE_AWARE",
    ,
    drop=FALSE
  ]

  strict_row <- x[
    x$framework == "PRIMARY_CORE" &
    x$reference_filter == "STRICT",
    ,
    drop=FALSE
  ]

  all_ref_row <- x[
    x$framework == "ALL" &
    x$reference_filter == "REFERENCE_AWARE",
    ,
    drop=FALSE
  ]

  noqc_ref_row <- x[
    x$framework == "NO_QCWATCH" &
    x$reference_filter == "REFERENCE_AWARE",
    ,
    drop=FALSE
  ]

  get1 <- function(df, col) {

    if (nrow(df) != 1) {
      return(NA)
    }

    df[[col]][[1]]
  }

  nes_vec <- x$NES[
    is.finite(
      x$NES
    )
  ]

  direction_concordant_all6 <-
    length(nes_vec) == 6 &&
    (
      all(nes_vec > 0) ||
      all(nes_vec < 0)
    )

  concordance_rows[[pw]] <- data.frame(
    pathway=pw,

    PRIMARY_reference_NES=
      get1(
        primary_row,
        "NES"
      ),

    PRIMARY_reference_FDR=
      get1(
        primary_row,
        "padj"
      ),

    PRIMARY_strict_NES=
      get1(
        strict_row,
        "NES"
      ),

    PRIMARY_strict_FDR=
      get1(
        strict_row,
        "padj"
      ),

    ALL_reference_NES=
      get1(
        all_ref_row,
        "NES"
      ),

    NO_QCWATCH_reference_NES=
      get1(
        noqc_ref_row,
        "NES"
      ),

    direction_concordant_all6=
      direction_concordant_all6,

    stringsAsFactors=FALSE
  )
}

concordance <- do.call(
  rbind,
  concordance_rows
)

rownames(concordance) <- NULL

concordance$primary_direction <- ifelse(
  concordance$PRIMARY_reference_NES > 0,
  "Tx_enriched",
  "Sham_enriched"
)

concordance$primary_significance <- ifelse(
  concordance$PRIMARY_reference_FDR < 0.10,
  "FDR_lt0.10",
  ifelse(
    concordance$PRIMARY_reference_FDR < 0.25,
    "exploratory_FDR_lt0.25",
    "not_significant"
  )
)

concordance$strict_direction_same <-
  sign(
    concordance$PRIMARY_reference_NES
  ) ==
  sign(
    concordance$PRIMARY_strict_NES
  )

concordance <- concordance[
  order(
    concordance$PRIMARY_reference_FDR,
    -abs(
      concordance$PRIMARY_reference_NES
    )
  ),
  ,
  drop=FALSE
]

write.csv(
  concordance,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_concordance_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Primary significant / exploratory pathways
# =========================================================

primary_selected <- primary[
  primary$padj < 0.25,
  ,
  drop=FALSE
]

write.csv(
  primary_selected,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_PRIMARY_selected_FDRlt0.25_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

strict_selected <- strict_primary[
  strict_primary$padj < 0.25,
  ,
  drop=FALSE
]

write.csv(
  strict_selected,
  file.path(
    TABDIR,
    "Monocyte_Hallmark_GSEA_PRIMARY_STRICT_selected_FDRlt0.25_v6.9.7.2.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Primary GSEA plot
# =========================================================

plot_primary <- primary[
  is.finite(primary$NES) &
  is.finite(primary$padj),
  ,
  drop=FALSE
]

plot_primary <- plot_primary[
  order(
    plot_primary$padj,
    -abs(plot_primary$NES)
  ),
  ,
  drop=FALSE
]

plot_primary <- head(
  plot_primary,
  20
)

if (nrow(plot_primary) > 0) {

  plot_primary$pathway_display <- gsub(
    "^HALLMARK_",
    "",
    plot_primary$pathway
  )

  plot_primary$pathway_display <- gsub(
    "_",
    " ",
    plot_primary$pathway_display
  )

  plot_primary$pathway_display <- factor(
    plot_primary$pathway_display,
    levels=rev(
      plot_primary$pathway_display
    )
  )

  plot_primary$minus_log10_FDR <-
    -log10(
      pmax(
        plot_primary$padj,
        .Machine$double.xmin
      )
    )

  p <- ggplot(
    plot_primary,
    aes(
      x=NES,
      y=pathway_display,
      size=minus_log10_FDR,
      color=NES
    )
  ) +
    geom_point() +
    scale_color_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(
      base_size=10
    ) +
    labs(
      title="Mouse MASH Monocyte Hallmark GSEA",
      subtitle="PRIMARY_CORE + support-aware reference filter",
      x="Normalized enrichment score (Tx vs Sham)",
      y=NULL,
      size="-log10(FDR)",
      color="NES"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Monocyte_Hallmark_GSEA_PRIMARY_v6.9.7.2.pdf"
    ),
    p,
    width=9,
    height=8
  )
}

# =========================================================
# Primary vs strict sensitivity plot
# =========================================================

compare_plot <- merge(
  primary[
    ,
    c(
      "pathway",
      "NES",
      "padj"
    )
  ],
  strict_primary[
    ,
    c(
      "pathway",
      "NES",
      "padj"
    )
  ],
  by="pathway",
  suffixes=c(
    "_reference",
    "_strict"
  )
)

compare_plot <- compare_plot[
  is.finite(compare_plot$NES_reference) &
  is.finite(compare_plot$NES_strict),
  ,
  drop=FALSE
]

if (nrow(compare_plot) > 0) {

  p2 <- ggplot(
    compare_plot,
    aes(
      x=NES_reference,
      y=NES_strict
    )
  ) +
    geom_point(
      size=2.5
    ) +
    geom_abline(
      slope=1,
      intercept=0,
      linetype=2
    ) +
    theme_classic(
      base_size=10
    ) +
    labs(
      title="Monocyte Hallmark GSEA reference-filter sensitivity",
      subtitle="PRIMARY_CORE",
      x="NES: reference-aware delta <= 2",
      y="NES: strict delta <= 1.5"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Monocyte_Hallmark_GSEA_PRIMARY_reference_vs_strict_v6.9.7.2.pdf"
    ),
    p2,
    width=7,
    height=7
  )
}

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== GSEA RANK SUMMARY ===\n")
print(
  rank_summary,
  row.names=FALSE
)

cat("\n=== PRIMARY HALLMARK GSEA TOP20 ===\n")

primary_print_cols <- c(
  "pathway",
  "NES",
  "pval",
  "padj",
  "size",
  "direction",
  "significance_class"
)

print(
  head(
    primary[
      ,
      intersect(
        primary_print_cols,
        colnames(primary)
      ),
      drop=FALSE
    ],
    20
  ),
  row.names=FALSE
)

cat("\n=== PRIMARY HALLMARK SELECTED FDR<0.25 ===\n")

if (nrow(primary_selected) == 0) {

  cat("NONE\n")

} else {

  print(
    primary_selected[
      ,
      intersect(
        primary_print_cols,
        colnames(primary_selected)
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )
}

cat("\n=== PRIMARY STRICT HALLMARK SELECTED FDR<0.25 ===\n")

if (nrow(strict_selected) == 0) {

  cat("NONE\n")

} else {

  print(
    strict_selected[
      ,
      intersect(
        primary_print_cols,
        colnames(strict_selected)
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )
}

cat("\n=== HALLMARK CONCORDANCE TOP20 ===\n")

concordance_print_cols <- c(
  "pathway",
  "PRIMARY_reference_NES",
  "PRIMARY_reference_FDR",
  "PRIMARY_strict_NES",
  "PRIMARY_strict_FDR",
  "ALL_reference_NES",
  "NO_QCWATCH_reference_NES",
  "direction_concordant_all6",
  "strict_direction_same"
)

print(
  head(
    concordance[
      ,
      intersect(
        concordance_print_cols,
        colnames(concordance)
      ),
      drop=FALSE
    ],
    20
  ),
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.7.2 COMPLETE\n")
cat("Reference-aware Hallmark GSEA complete\n")
cat("Primary analysis: PRIMARY_CORE + supported-reference delta <= 2\n")
cat("Strict sensitivity: PRIMARY_CORE + supported-reference delta <= 1.5\n")
cat("Framework sensitivity: ALL / NO_QCWATCH / PRIMARY_CORE\n")
cat("Rank metric: sign(logFC) * sqrt(edgeR QL F)\n")
cat("Sham vs Tx biological n=2/group\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.7.2.txt"
  )
)
