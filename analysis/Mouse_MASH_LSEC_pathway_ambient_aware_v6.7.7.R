suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

VERSION <- "v6.7.7"

INDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.6.1/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC ambient-aware pathway analysis\n")
cat("Version:", VERSION, "\n")
cat("Contrast: Tx vs Sham\n")
cat("====================================================\n\n")

files <- c(
  Primary_no_QC =
    "Primary_no_QC_Tx_vs_Sham_edgeR_with_Hep_LSEC_specificity_v6.7.6.1.csv",

  Shared_core =
    "Shared_core_Tx_vs_Sham_edgeR_with_Hep_LSEC_specificity_v6.7.6.1.csv",

  Inflammatory_stress_high =
    "Inflammatory_stress_high_Tx_vs_Sham_edgeR_with_Hep_LSEC_specificity_v6.7.6.1.csv",

  Wnt_angiocrine_high =
    "Wnt_angiocrine_high_Tx_vs_Sham_edgeR_with_Hep_LSEC_specificity_v6.7.6.1.csv"
)

for (f in files) {
  if (!file.exists(file.path(INDIR, f))) {
    stop("Missing input: ", f)
  }
}

# ---------------------------------------------------------
# msigdbr compatibility helper
# ---------------------------------------------------------

get_msig <- function(collection, subcollection=NULL) {

  fm <- names(formals(msigdbr::msigdbr))

  args <- list(
    species="Mus musculus"
  )

  if ("collection" %in% fm) {
    args$collection <- collection
  } else {
    args$category <- collection
  }

  if (!is.null(subcollection)) {

    if ("subcollection" %in% fm) {
      args$subcollection <- subcollection
    } else {
      args$subcategory <- subcollection
    }
  }

  do.call(
    msigdbr::msigdbr,
    args
  )
}

cat("=== LOAD MSigDB ===\n")

msig_h <- get_msig("H")

msig_bp <- get_msig(
  "C5",
  "GO:BP"
)

if (!all(
  c("gs_name","gene_symbol") %in%
    colnames(msig_h)
)) {
  stop("Unexpected msigdbr columns.")
}

hallmark_sets <- split(
  msig_h$gene_symbol,
  msig_h$gs_name
)

gobp_sets <- split(
  msig_bp$gene_symbol,
  msig_bp$gs_name
)

hallmark_sets <- lapply(
  hallmark_sets,
  unique
)

gobp_sets <- lapply(
  gobp_sets,
  unique
)

cat(
  "Hallmark pathways:",
  length(hallmark_sets),
  "\n"
)

cat(
  "GO BP pathways:",
  length(gobp_sets),
  "\n"
)

# ---------------------------------------------------------
# Rank builder
# Primary analysis excludes Hepatocyte_dominant genes
# No P-value threshold is applied before GSEA
# ---------------------------------------------------------

make_rank <- function(x, remove_hep=TRUE) {

  required <- c(
    "gene",
    "logFC",
    "PValue",
    "Specificity_class_v6761"
  )

  miss <- setdiff(
    required,
    colnames(x)
  )

  if (length(miss) > 0) {
    stop(
      "Missing columns: ",
      paste(miss, collapse=", ")
    )
  }

  x <- x[
    is.finite(x$logFC) &
    is.finite(x$PValue) &
    !is.na(x$gene),
    ,
    drop=FALSE
  ]

  if (remove_hep) {

    x <- x[
      x$Specificity_class_v6761 !=
        "Hepatocyte_dominant" |
      is.na(x$Specificity_class_v6761),
      ,
      drop=FALSE
    ]
  }

  # edgeR directional rank
  score <- sign(x$logFC) *
    -log10(
      pmax(
        x$PValue,
        1e-300
      )
    )

  names(score) <- x$gene

  # one value per gene
  o <- order(
    abs(score),
    decreasing=TRUE
  )

  score <- score[o]

  score <- score[
    !duplicated(names(score))
  ]

  score <- sort(
    score,
    decreasing=TRUE
  )

  score
}

# ---------------------------------------------------------
# fgsea runner
# ---------------------------------------------------------

run_fgsea <- function(
  stats,
  pathways,
  analysis,
  collection,
  universe
) {

  set.seed(20260902)

  fg <- fgseaMultilevel(
    pathways=pathways,
    stats=stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  fg <- as.data.frame(fg)

  fg$analysis <- analysis
  fg$collection <- collection
  fg$universe <- universe

  fg$direction <- ifelse(
    fg$NES > 0,
    "Tx_up",
    "Tx_down"
  )

  fg$leadingEdge <- vapply(
    fg$leadingEdge,
    function(x) {
      paste(x, collapse=";")
    },
    character(1)
  )

  fg <- fg[
    order(
      fg$padj,
      -abs(fg$NES)
    ),
    ,
    drop=FALSE
  ]

  fg
}

results <- list()

for (nm in names(files)) {

  cat(
    "\n========================================\n"
  )
  cat("ANALYSIS:", nm, "\n")
  cat(
    "========================================\n"
  )

  x <- read.csv(
    file.path(
      INDIR,
      files[[nm]]
    ),
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  # Primary ambient-aware universe
  rank_nonhep <- make_rank(
    x,
    remove_hep=TRUE
  )

  # Sensitivity universe
  rank_all <- make_rank(
    x,
    remove_hep=FALSE
  )

  cat(
    "Ranked genes, non-Hep:",
    length(rank_nonhep),
    "\n"
  )

  cat(
    "Ranked genes, all:",
    length(rank_all),
    "\n"
  )

  for (universe in c(
    "NonHep_primary",
    "AllGenes_sensitivity"
  )) {

    stats <- if (
      universe == "NonHep_primary"
    ) {
      rank_nonhep
    } else {
      rank_all
    }

    h <- run_fgsea(
      stats,
      hallmark_sets,
      nm,
      "Hallmark",
      universe
    )

    bp <- run_fgsea(
      stats,
      gobp_sets,
      nm,
      "GO_BP",
      universe
    )

    key_h <- paste(
      nm,
      universe,
      "Hallmark",
      sep="__"
    )

    key_bp <- paste(
      nm,
      universe,
      "GO_BP",
      sep="__"
    )

    results[[key_h]] <- h
    results[[key_bp]] <- bp

    write.csv(
      h,
      file.path(
        TABDIR,
        paste0(
          nm,
          "_",
          universe,
          "_Hallmark_fgsea_v6.7.7.csv"
        )
      ),
      row.names=FALSE
    )

    write.csv(
      bp,
      file.path(
        TABDIR,
        paste0(
          nm,
          "_",
          universe,
          "_GOBP_fgsea_v6.7.7.csv"
        )
      ),
      row.names=FALSE
    )
  }
}

# ---------------------------------------------------------
# Primary vs Shared-core pathway concordance
# ---------------------------------------------------------

make_concordance <- function(collection) {

  a <- results[[
    paste(
      "Primary_no_QC",
      "NonHep_primary",
      collection,
      sep="__"
    )
  ]]

  b <- results[[
    paste(
      "Shared_core",
      "NonHep_primary",
      collection,
      sep="__"
    )
  ]]

  out <- merge(
    a[
      ,
      c(
        "pathway",
        "NES",
        "pval",
        "padj",
        "size",
        "leadingEdge"
      )
    ],
    b[
      ,
      c(
        "pathway",
        "NES",
        "pval",
        "padj",
        "size",
        "leadingEdge"
      )
    ],
    by="pathway",
    suffixes=c(
      "_Primary",
      "_Shared"
    )
  )

  out$same_direction <-
    sign(out$NES_Primary) ==
    sign(out$NES_Shared)

  out$mean_NES <-
    (
      out$NES_Primary +
      out$NES_Shared
    ) / 2

  out$max_padj <-
    pmax(
      out$padj_Primary,
      out$padj_Shared
    )

  out$direction <- ifelse(
    out$mean_NES > 0,
    "Tx_up",
    "Tx_down"
  )

  out <- out[
    order(
      !out$same_direction,
      out$max_padj,
      -abs(out$mean_NES)
    ),
    ,
    drop=FALSE
  ]

  out
}

hallmark_conc <- make_concordance(
  "Hallmark"
)

gobp_conc <- make_concordance(
  "GO_BP"
)

write.csv(
  hallmark_conc,
  file.path(
    TABDIR,
    "Hallmark_concordance_Primary_vs_Shared_v6.7.7.csv"
  ),
  row.names=FALSE
)

write.csv(
  gobp_conc,
  file.path(
    TABDIR,
    "GOBP_concordance_Primary_vs_Shared_v6.7.7.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Exploratory consensus sets
#
# Because n=2/group, padj<0.25 is used only as a
# discovery threshold, not confirmatory significance.
# ---------------------------------------------------------

hallmark_consensus <- hallmark_conc[
  hallmark_conc$same_direction &
  hallmark_conc$max_padj < 0.25,
  ,
  drop=FALSE
]

gobp_consensus <- gobp_conc[
  gobp_conc$same_direction &
  gobp_conc$max_padj < 0.25,
  ,
  drop=FALSE
]

write.csv(
  hallmark_consensus,
  file.path(
    TABDIR,
    "Hallmark_consensus_exploratory_padj025_v6.7.7.csv"
  ),
  row.names=FALSE
)

write.csv(
  gobp_consensus,
  file.path(
    TABDIR,
    "GOBP_consensus_exploratory_padj025_v6.7.7.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Cross-analysis LSEC-supported gene candidates
# ---------------------------------------------------------

cand_file <- file.path(
  INDIR,
  "LSEC_Tx_response_NON_Hepatocyte_dominant_candidates_v6.7.6.1.csv"
)

cand <- read.csv(
  cand_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cand <- cand[
  cand$analysis %in%
    c(
      "Primary_no_QC",
      "Shared_core"
    ) &
  cand$Specificity_class_v6761 ==
    "LSEC_supported",
  ,
  drop=FALSE
]

primary <- cand[
  cand$analysis == "Primary_no_QC",
  ,
  drop=FALSE
]

shared <- cand[
  cand$analysis == "Shared_core",
  ,
  drop=FALSE
]

cross <- merge(
  primary[
    ,
    c(
      "gene",
      "logFC",
      "PValue",
      "FDR",
      "Hep_minus_LSEC_primary_median"
    )
  ],
  shared[
    ,
    c(
      "gene",
      "logFC",
      "PValue",
      "FDR"
    )
  ],
  by="gene",
  suffixes=c(
    "_Primary",
    "_Shared"
  )
)

cross$same_direction <-
  sign(cross$logFC_Primary) ==
  sign(cross$logFC_Shared)

cross$mean_logFC <-
  (
    cross$logFC_Primary +
    cross$logFC_Shared
  ) / 2

cross$max_FDR <-
  pmax(
    cross$FDR_Primary,
    cross$FDR_Shared
  )

cross$min_PValue <-
  pmin(
    cross$PValue_Primary,
    cross$PValue_Shared
  )

cross <- cross[
  cross$same_direction,
  ,
  drop=FALSE
]

cross <- cross[
  order(
    cross$max_FDR,
    -abs(cross$mean_logFC)
  ),
  ,
  drop=FALSE
]

write.csv(
  cross,
  file.path(
    TABDIR,
    "LSEC_supported_cross_analysis_candidates_v6.7.7.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Hallmark consensus plot
# ---------------------------------------------------------

if (nrow(hallmark_consensus) > 0) {

  pdat <- hallmark_consensus

  pdat <- pdat[
    order(
      abs(pdat$mean_NES),
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  pdat <- head(
    pdat,
    20
  )

  pdat$pathway <- factor(
    pdat$pathway,
    levels=rev(pdat$pathway)
  )

  p <- ggplot(
    pdat,
    aes(
      x=mean_NES,
      y=pathway
    )
  ) +
    geom_col() +
    theme_classic(base_size=10) +
    labs(
      title=
        "LSEC Tx vs Sham: concordant Hallmark pathways",
      subtitle=
        "Primary + Shared-core, non-Hepatocyte universe",
      x="Mean NES",
      y=NULL
    )

  ggsave(
    file.path(
      FIGDIR,
      "Hallmark_concordant_pathways_v6.7.7.pdf"
    ),
    p,
    width=9,
    height=7
  )
}

# ---------------------------------------------------------
# Print summary
# ---------------------------------------------------------

cat("\n=== HALLMARK CONSENSUS ===\n")

if (nrow(hallmark_consensus) == 0) {

  cat("No concordant pathways with max padj < 0.25\n")

} else {

  print(
    hallmark_consensus[
      ,
      c(
        "pathway",
        "NES_Primary",
        "padj_Primary",
        "NES_Shared",
        "padj_Shared",
        "direction"
      )
    ]
  )
}

cat("\n=== GO BP CONSENSUS TOP 30 ===\n")

if (nrow(gobp_consensus) == 0) {

  cat("No concordant pathways with max padj < 0.25\n")

} else {

  print(
    head(
      gobp_consensus[
        ,
        c(
          "pathway",
          "NES_Primary",
          "padj_Primary",
          "NES_Shared",
          "padj_Shared",
          "direction"
        )
      ],
      30
    )
  )
}

cat("\n=== CROSS-ANALYSIS LSEC-SUPPORTED GENES TOP 30 ===\n")

print(
  head(
    cross[
      ,
      c(
        "gene",
        "logFC_Primary",
        "logFC_Shared",
        "mean_logFC",
        "max_FDR",
        "Hep_minus_LSEC_primary_median"
      )
    ],
    30
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.7.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.7 COMPLETE\n")
cat("Primary pathway universe: non-Hepatocyte-dominant genes\n")
cat("Sensitivity pathway universe: all tested genes\n")
cat("Primary comparisons: Primary_no_QC + Shared_core\n")
cat("padj < 0.25 is exploratory GSEA criterion only\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
