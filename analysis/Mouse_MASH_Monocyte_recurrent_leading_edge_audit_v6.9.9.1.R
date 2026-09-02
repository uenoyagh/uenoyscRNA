suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v6.9.9.1"

BASE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS"
)

STATE_GSEA_FILE <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.8",
  "tables",
  "Monocyte_state_specific_Hallmark_GSEA_all_v6.9.8.csv"
)

SPEC_FILE <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.6.4",
  "tables",
  "Monocyte_genomewide_support_aware_reference_specificity_v6.9.6.4.csv"
)

GLOBAL_DE_FILE <- file.path(
  BASE_DIR,
  "Mouse_MASH_Monocyte_v6.9.6",
  "tables",
  "Monocyte_pseudobulk_PRIMARY_CORE_Sham_vs_Tx_v6.9.6.csv"
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
cat("Mouse MASH Monocyte recurrent Hallmark leading-edge audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

for (f in c(
  STATE_GSEA_FILE,
  SPEC_FILE,
  GLOBAL_DE_FILE
)) {
  if (!file.exists(f)) {
    stop("Missing input: ", f)
  }
}

gsea <- read.csv(
  STATE_GSEA_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

spec <- read.csv(
  SPEC_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

de <- read.csv(
  GLOBAL_DE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_gsea_cols <- c(
  "pathway",
  "NES",
  "padj",
  "leadingEdge",
  "cluster",
  "state",
  "reference_filter"
)

missing_gsea_cols <- setdiff(
  required_gsea_cols,
  colnames(gsea)
)

if (length(missing_gsea_cols) > 0) {
  stop(
    "Missing GSEA columns: ",
    paste(missing_gsea_cols, collapse=", ")
  )
}

# =========================================================
# Frozen interpretation hierarchy
#
# Core pathways:
# robust in v6.9.8 across states and retained under strict filtering.
#
# Secondary pathways:
# biologically useful but less robust or more filter-sensitive.
# =========================================================

core_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_P53_PATHWAY"
)

secondary_pathways <- c(
  "HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

target_pathways <- unique(
  c(
    core_pathways,
    secondary_pathways
  )
)

gsea <- gsea[
  gsea$pathway %in% target_pathways,
  ,
  drop=FALSE
]

if (nrow(gsea) == 0) {
  stop("No target pathways found in v6.9.8 GSEA table.")
}

# =========================================================
# Expand semicolon-delimited leadingEdge genes
# =========================================================

expand_leading_edges <- function(df) {

  rows <- list()
  k <- 1

  for (i in seq_len(nrow(df))) {

    le <- df$leadingEdge[[i]]

    if (is.na(le) || le == "") {
      next
    }

    genes <- unique(
      strsplit(
        le,
        ";",
        fixed=TRUE
      )[[1]]
    )

    genes <- genes[
      genes != ""
    ]

    if (length(genes) == 0) {
      next
    }

    rows[[k]] <- data.frame(
      pathway=df$pathway[[i]],
      cluster=as.character(df$cluster[[i]]),
      state=df$state[[i]],
      reference_filter=df$reference_filter[[i]],
      NES=df$NES[[i]],
      padj=df$padj[[i]],
      gene=genes,
      stringsAsFactors=FALSE
    )

    k <- k + 1
  }

  if (length(rows) == 0) {
    return(
      data.frame(
        pathway=character(),
        cluster=character(),
        state=character(),
        reference_filter=character(),
        NES=numeric(),
        padj=numeric(),
        gene=character(),
        stringsAsFactors=FALSE
      )
    )
  }

  out <- do.call(
    rbind,
    rows
  )

  rownames(out) <- NULL
  out
}

leading <- expand_leading_edges(
  gsea
)

write.csv(
  leading,
  file.path(
    TABDIR,
    "Monocyte_state_specific_target_pathway_leading_edges_v6.9.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Recurrent gene counts by pathway/filter
# =========================================================

count_recurrence <- function(df) {

  if (nrow(df) == 0) {
    return(data.frame())
  }

  key <- unique(
    df[
      ,
      c(
        "pathway",
        "reference_filter",
        "gene",
        "cluster"
      ),
      drop=FALSE
    ]
  )

  agg <- aggregate(
    cluster ~ pathway +
      reference_filter +
      gene,
    data=key,
    FUN=function(x) {
      length(unique(x))
    }
  )

  names(agg)[
    names(agg) == "cluster"
  ] <- "n_states_leading_edge"

  state_names <- aggregate(
    state ~ pathway +
      reference_filter +
      gene,
    data=unique(
      df[
        ,
        c(
          "pathway",
          "reference_filter",
          "gene",
          "state"
        ),
        drop=FALSE
      ]
    ),
    FUN=function(x) {
      paste(
        sort(unique(x)),
        collapse=" | "
      )
    }
  )

  out <- merge(
    agg,
    state_names,
    by=c(
      "pathway",
      "reference_filter",
      "gene"
    ),
    all.x=TRUE
  )

  out
}

recurrence <- count_recurrence(
  leading
)

# =========================================================
# Add specificity and global PRIMARY_CORE DE context
# =========================================================

spec_cols <- c(
  "gene",
  "Monocyte_median_logCPM_ShTx",
  "max_supported_reference_lineage_ShTx",
  "max_supported_reference_delta_ShTx",
  "specificity_class_supported_ShTx",
  "keep_reference_aware_supported",
  "keep_strict_reference_supported"
)

spec_small <- spec[
  ,
  intersect(
    spec_cols,
    colnames(spec)
  ),
  drop=FALSE
]

de_cols <- c(
  "gene",
  "logFC",
  "F",
  "PValue",
  "FDR",
  "replicate_direction"
)

de_small <- de[
  ,
  intersect(
    de_cols,
    colnames(de)
  ),
  drop=FALSE
]

names(de_small)[
  names(de_small) == "logFC"
] <- "PRIMARY_CORE_logFC"

names(de_small)[
  names(de_small) == "F"
] <- "PRIMARY_CORE_F"

names(de_small)[
  names(de_small) == "PValue"
] <- "PRIMARY_CORE_PValue"

names(de_small)[
  names(de_small) == "FDR"
] <- "PRIMARY_CORE_FDR"

names(de_small)[
  names(de_small) == "replicate_direction"
] <- "PRIMARY_CORE_replicate_direction"

recurrence <- merge(
  recurrence,
  spec_small,
  by="gene",
  all.x=TRUE
)

recurrence <- merge(
  recurrence,
  de_small,
  by="gene",
  all.x=TRUE
)

recurrence <- recurrence[
  order(
    recurrence$reference_filter,
    recurrence$pathway,
    -recurrence$n_states_leading_edge,
    -abs(recurrence$PRIMARY_CORE_logFC)
  ),
  ,
  drop=FALSE
]

write.csv(
  recurrence,
  file.path(
    TABDIR,
    "Monocyte_recurrent_leading_edge_gene_audit_v6.9.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Reference-aware recurrent drivers
#
# Primary recurrent driver:
# leading edge in >=4/6 eligible states.
#
# Strong recurrent driver:
# >=5/6 states.
# =========================================================

ref_rec <- recurrence[
  recurrence$reference_filter ==
    "REFERENCE_AWARE",
  ,
  drop=FALSE
]

ref_rec$driver_recurrence_class <- ifelse(
  ref_rec$n_states_leading_edge >= 5,
  "strong_recurrent_5to6_states",
  ifelse(
    ref_rec$n_states_leading_edge >= 4,
    "recurrent_4_states",
    "limited_1to3_states"
  )
)

strict_rec <- recurrence[
  recurrence$reference_filter ==
    "STRICT",
  ,
  drop=FALSE
]

strict_key <- strict_rec[
  ,
  c(
    "pathway",
    "gene",
    "n_states_leading_edge"
  ),
  drop=FALSE
]

names(strict_key)[
  names(strict_key) ==
    "n_states_leading_edge"
] <- "strict_n_states_leading_edge"

ref_rec <- merge(
  ref_rec,
  strict_key,
  by=c(
    "pathway",
    "gene"
  ),
  all.x=TRUE
)

ref_rec$strict_n_states_leading_edge[
  is.na(
    ref_rec$strict_n_states_leading_edge
  )
] <- 0

ref_rec$strict_recurrent_4plus <-
  ref_rec$strict_n_states_leading_edge >= 4

ref_rec$interpretation_tier <- ifelse(
  ref_rec$pathway %in% core_pathways &
    ref_rec$n_states_leading_edge >= 4 &
    ref_rec$strict_n_states_leading_edge >= 4,
  "CORE_ROBUST_DRIVER",
  ifelse(
    ref_rec$n_states_leading_edge >= 4,
    "SECONDARY_RECURRENT_DRIVER",
    "LIMITED_DRIVER"
  )
)

ref_rec <- ref_rec[
  order(
    factor(
      ref_rec$interpretation_tier,
      levels=c(
        "CORE_ROBUST_DRIVER",
        "SECONDARY_RECURRENT_DRIVER",
        "LIMITED_DRIVER"
      )
    ),
    ref_rec$pathway,
    -ref_rec$n_states_leading_edge,
    -ref_rec$strict_n_states_leading_edge,
    -abs(ref_rec$PRIMARY_CORE_logFC)
  ),
  ,
  drop=FALSE
]

write.csv(
  ref_rec,
  file.path(
    TABDIR,
    "Monocyte_REFERENCE_AWARE_leading_edge_driver_summary_v6.9.9.1.csv"
  ),
  row.names=FALSE
)

core_drivers <- ref_rec[
  ref_rec$interpretation_tier ==
    "CORE_ROBUST_DRIVER",
  ,
  drop=FALSE
]

write.csv(
  core_drivers,
  file.path(
    TABDIR,
    "Monocyte_CORE_ROBUST_leading_edge_drivers_v6.9.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Cross-pathway recurrent genes
# =========================================================

core_ref <- ref_rec[
  ref_rec$pathway %in% core_pathways &
    ref_rec$n_states_leading_edge >= 4,
  ,
  drop=FALSE
]

if (nrow(core_ref) > 0) {

  cross_pathway <- aggregate(
    pathway ~ gene,
    data=unique(
      core_ref[
        ,
        c(
          "gene",
          "pathway"
        ),
        drop=FALSE
      ]
    ),
    FUN=function(x) {
      length(unique(x))
    }
  )

  names(cross_pathway)[
    names(cross_pathway) ==
      "pathway"
  ] <- "n_core_pathways"

  pathway_names <- aggregate(
    pathway ~ gene,
    data=unique(
      core_ref[
        ,
        c(
          "gene",
          "pathway"
        ),
        drop=FALSE
      ]
    ),
    FUN=function(x) {
      paste(
        sort(unique(x)),
        collapse=" | "
      )
    }
  )

  cross_pathway <- merge(
    cross_pathway,
    pathway_names,
    by="gene",
    all.x=TRUE
  )

  # Maximum state recurrence across the core pathways
  max_recurrence <- aggregate(
    n_states_leading_edge ~ gene,
    data=core_ref,
    FUN=max
  )

  names(max_recurrence)[
    names(max_recurrence) ==
      "n_states_leading_edge"
  ] <- "max_states_in_core_pathway"

  cross_pathway <- merge(
    cross_pathway,
    max_recurrence,
    by="gene",
    all.x=TRUE
  )

  cross_pathway <- merge(
    cross_pathway,
    spec_small,
    by="gene",
    all.x=TRUE
  )

  cross_pathway <- merge(
    cross_pathway,
    de_small,
    by="gene",
    all.x=TRUE
  )

  cross_pathway <- cross_pathway[
    order(
      -cross_pathway$n_core_pathways,
      -cross_pathway$max_states_in_core_pathway,
      -abs(cross_pathway$PRIMARY_CORE_logFC)
    ),
    ,
    drop=FALSE
  ]

} else {

  cross_pathway <- data.frame()
}

write.csv(
  cross_pathway,
  file.path(
    TABDIR,
    "Monocyte_cross_core_pathway_recurrent_genes_v6.9.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Core driver figure
# =========================================================

plot_df <- core_drivers

if (nrow(plot_df) > 0) {

  plot_df <- plot_df[
    order(
      plot_df$pathway,
      -plot_df$n_states_leading_edge,
      -plot_df$strict_n_states_leading_edge,
      -abs(plot_df$PRIMARY_CORE_logFC)
    ),
    ,
    drop=FALSE
  ]

  plot_df <- do.call(
    rbind,
    lapply(
      split(
        plot_df,
        plot_df$pathway
      ),
      function(z) {
        head(z, 12)
      }
    )
  )

  rownames(plot_df) <- NULL

  plot_df$label <- paste0(
    plot_df$gene,
    "  [",
    plot_df$n_states_leading_edge,
    "/6; strict ",
    plot_df$strict_n_states_leading_edge,
    "/6]"
  )

  plot_df$label <- factor(
    plot_df$label,
    levels=rev(
      unique(
        plot_df$label
      )
    )
  )

  p <- ggplot(
    plot_df,
    aes(
      x=n_states_leading_edge,
      y=label,
      size=abs(PRIMARY_CORE_logFC),
      color=PRIMARY_CORE_logFC
    )
  ) +
    geom_point() +
    facet_wrap(
      ~pathway,
      scales="free_y"
    ) +
    scale_color_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    scale_x_continuous(
      breaks=1:6,
      limits=c(0.5, 6.5)
    ) +
    theme_classic(
      base_size=9
    ) +
    labs(
      title="Monocyte recurrent leading-edge drivers",
      subtitle="Core pathways; reference-aware recurrence with strict sensitivity shown in labels",
      x="Number of eligible states containing gene in leading edge",
      y=NULL,
      size="|global logFC|",
      color="Global logFC"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Monocyte_CORE_recurrent_leading_edge_drivers_v6.9.9.1.pdf"
    ),
    p,
    width=12,
    height=8
  )
}

# =========================================================
# Terminal summaries
# =========================================================

cat("\n=== CORE PATHWAY RECURRENT DRIVERS ===\n")

core_print_cols <- c(
  "pathway",
  "gene",
  "n_states_leading_edge",
  "strict_n_states_leading_edge",
  "PRIMARY_CORE_logFC",
  "PRIMARY_CORE_FDR",
  "PRIMARY_CORE_replicate_direction",
  "max_supported_reference_lineage_ShTx",
  "max_supported_reference_delta_ShTx"
)

if (nrow(core_drivers) == 0) {

  cat("NONE\n")

} else {

  print(
    core_drivers[
      ,
      intersect(
        core_print_cols,
        colnames(core_drivers)
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )
}

cat("\n=== CROSS-CORE-PATHWAY RECURRENT GENES ===\n")

if (nrow(cross_pathway) == 0) {

  cat("NONE\n")

} else {

  cross_print_cols <- c(
    "gene",
    "n_core_pathways",
    "max_states_in_core_pathway",
    "pathway",
    "PRIMARY_CORE_logFC",
    "PRIMARY_CORE_FDR",
    "PRIMARY_CORE_replicate_direction",
    "max_supported_reference_lineage_ShTx",
    "max_supported_reference_delta_ShTx"
  )

  print(
    head(
      cross_pathway[
        ,
        intersect(
          cross_print_cols,
          colnames(cross_pathway)
        ),
        drop=FALSE
      ],
      40
    ),
    row.names=FALSE
  )
}

cat("\n=== SECONDARY RECURRENT DRIVERS (TOP40) ===\n")

secondary <- ref_rec[
  ref_rec$interpretation_tier ==
    "SECONDARY_RECURRENT_DRIVER",
  ,
  drop=FALSE
]

if (nrow(secondary) == 0) {

  cat("NONE\n")

} else {

  print(
    head(
      secondary[
        ,
        intersect(
          core_print_cols,
          colnames(secondary)
        ),
        drop=FALSE
      ],
      40
    ),
    row.names=FALSE
  )
}

cat("\n====================================================\n")
cat("v6.9.9.1 COMPLETE\n")
cat("Recurrent leading-edge driver audit complete\n")
cat("Core pathways: TNFA/NFKB, HYPOXIA, P53\n")
cat("Core robust driver: leading edge in >=4 states in BOTH reference-aware and strict analyses\n")
cat("Secondary pathways retained separately\n")
cat("No new differential-expression testing performed\n")
cat("No source data modified\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.9.1.txt"
  )
)
