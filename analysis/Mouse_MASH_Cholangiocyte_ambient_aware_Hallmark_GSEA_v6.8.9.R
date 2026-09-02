suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.9"

INPUT_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.8.1/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte ambient-aware Hallmark GSEA\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# ---------------------------------------------------------
# Input DE tables
# ---------------------------------------------------------

files <- c(
  ALL=
    "Pseudobulk_DE_ALL_reference_annotated_v6.8.8.1.csv",

  NO_QCWATCH=
    "Pseudobulk_DE_NO_QCWATCH_reference_annotated_v6.8.8.1.csv",

  PRIMARY_CORE=
    "Pseudobulk_DE_PRIMARY_CORE_reference_annotated_v6.8.8.1.csv"
)

de_list <- list()

for (nm in names(files)) {

  path <- file.path(
    INPUT_DIR,
    files[[nm]]
  )

  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }

  x <- read.csv(
    path,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  required <- c(
    "gene",
    "logFC",
    "F",
    "PValue",
    "FDR",
    "Hepatocyte_specificity_class",
    "max_reference_delta"
  )

  missing_cols <- setdiff(
    required,
    colnames(x)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ",
      nm,
      ": ",
      paste(missing_cols, collapse=", ")
    )
  }

  de_list[[nm]] <- x
}

# ---------------------------------------------------------
# Hallmark gene sets
# ---------------------------------------------------------

cat("=== LOAD HALLMARK GENE SETS ===\n")

hallmark_df <- tryCatch(
  {
    msigdbr(
      species="Mus musculus",
      collection="H"
    )
  },
  error=function(e) {
    msigdbr(
      species="Mus musculus",
      category="H"
    )
  }
)

if (!all(
  c("gs_name", "gene_symbol") %in%
    colnames(hallmark_df)
)) {
  stop("Unexpected msigdbr output.")
}

hallmark <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

hallmark <- lapply(
  hallmark,
  unique
)

cat(
  "Hallmark pathways:",
  length(hallmark),
  "\n"
)

# ---------------------------------------------------------
# Ranking function
#
# signed sqrt(edgeR quasi-likelihood F)
# avoids hard DE thresholds.
# ---------------------------------------------------------

make_rank <- function(x) {

  keep <- is.finite(x$logFC) &
    is.finite(x$F) &
    !is.na(x$gene) &
    x$gene != ""

  x <- x[
    keep,
    ,
    drop=FALSE
  ]

  stat <- sign(x$logFC) *
    sqrt(
      pmax(
        x$F,
        0
      )
    )

  names(stat) <- x$gene

  # Defensive duplicate handling
  if (anyDuplicated(names(stat))) {

    genes <- unique(names(stat))

    stat <- sapply(
      genes,
      function(g) {

        z <- stat[
          names(stat) == g
        ]

        z[
          which.max(
            abs(z)
          )
        ]
      }
    )
  }

  stat <- sort(
    stat,
    decreasing=TRUE
  )

  stat
}

# ---------------------------------------------------------
# GSEA runner
# ---------------------------------------------------------

run_one_gsea <- function(
  de,
  framework,
  filter_mode
) {

  if (
    filter_mode ==
      "HEP_AMBIENT_AWARE"
  ) {

    x <- de[
      de$Hepatocyte_specificity_class !=
        "Hepatocyte_dominant",
      ,
      drop=FALSE
    ]

  } else if (
    filter_mode ==
      "STRICT_REFERENCE_SENSITIVITY"
  ) {

    x <- de[
      de$Hepatocyte_specificity_class !=
        "Hepatocyte_dominant" &
        (
          is.na(de$max_reference_delta) |
            de$max_reference_delta <= 2
        ),
      ,
      drop=FALSE
    ]

  } else {

    stop(
      "Unknown filter mode: ",
      filter_mode
    )
  }

  stats <- make_rank(x)

  cat(
    "\n=== GSEA ",
    framework,
    " / ",
    filter_mode,
    " ===\n",
    sep=""
  )

  cat(
    "Ranked genes:",
    length(stats),
    "\n"
  )

  fg <- fgseaMultilevel(
    pathways=hallmark,
    stats=stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  fg <- as.data.frame(fg)

  if (nrow(fg) == 0) {
    stop(
      "No GSEA results: ",
      framework,
      " / ",
      filter_mode
    )
  }

  fg$framework <- framework
  fg$filter_mode <- filter_mode

  fg$leadingEdge_genes <- vapply(
    fg$leadingEdge,
    function(z) {
      paste(
        z,
        collapse=";"
      )
    },
    character(1)
  )

  fg$leadingEdge <- NULL

  fg <- fg[
    order(
      fg$padj,
      -abs(fg$NES)
    ),
    ,
    drop=FALSE
  ]

  rownames(fg) <- NULL

  write.csv(
    fg,
    file.path(
      TABDIR,
      paste0(
        "Hallmark_GSEA_",
        framework,
        "_",
        filter_mode,
        "_v6.8.9.csv"
      )
    ),
    row.names=FALSE
  )

  fg
}

# ---------------------------------------------------------
# Run all combinations
# ---------------------------------------------------------

filter_modes <- c(
  "HEP_AMBIENT_AWARE",
  "STRICT_REFERENCE_SENSITIVITY"
)

results <- list()

for (framework in names(de_list)) {

  for (mode in filter_modes) {

    key <- paste(
      framework,
      mode,
      sep="__"
    )

    results[[key]] <- run_one_gsea(
      de=de_list[[framework]],
      framework=framework,
      filter_mode=mode
    )
  }
}

all_gsea <- do.call(
  rbind,
  results
)

rownames(all_gsea) <- NULL

write.csv(
  all_gsea,
  file.path(
    TABDIR,
    "Hallmark_GSEA_ALL_RESULTS_v6.8.9.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Concordance across ALL / NO_QCWATCH / PRIMARY_CORE
# ---------------------------------------------------------

make_concordance <- function(
  all_gsea,
  mode
) {

  x <- all_gsea[
    all_gsea$filter_mode == mode,
    ,
    drop=FALSE
  ]

  pathways <- sort(
    unique(x$pathway)
  )

  out <- do.call(
    rbind,
    lapply(
      pathways,
      function(pw) {

        z <- x[
          x$pathway == pw,
          ,
          drop=FALSE
        ]

        get_val <- function(
          framework,
          col
        ) {

          zz <- z[
            z$framework == framework,
            col
          ]

          if (length(zz) != 1) {
            return(NA_real_)
          }

          as.numeric(zz)
        }

        nes_all <- get_val(
          "ALL",
          "NES"
        )

        nes_noqc <- get_val(
          "NO_QCWATCH",
          "NES"
        )

        nes_primary <- get_val(
          "PRIMARY_CORE",
          "NES"
        )

        padj_all <- get_val(
          "ALL",
          "padj"
        )

        padj_noqc <- get_val(
          "NO_QCWATCH",
          "padj"
        )

        padj_primary <- get_val(
          "PRIMARY_CORE",
          "padj"
        )

        nes_vec <- c(
          nes_all,
          nes_noqc,
          nes_primary
        )

        same_direction <-
          all(
            is.finite(nes_vec)
          ) &&
          (
            all(nes_vec > 0) ||
              all(nes_vec < 0)
          )

        data.frame(
          filter_mode=mode,
          pathway=pw,

          NES_ALL=nes_all,
          NES_NO_QCWATCH=nes_noqc,
          NES_PRIMARY_CORE=nes_primary,

          padj_ALL=padj_all,
          padj_NO_QCWATCH=padj_noqc,
          padj_PRIMARY_CORE=padj_primary,

          same_direction_all3=
            same_direction,

          min_abs_NES=
            if (
              all(
                is.finite(nes_vec)
              )
            ) {
              min(
                abs(nes_vec)
              )
            } else {
              NA_real_
            },

          max_padj=
            max(
              c(
                padj_all,
                padj_noqc,
                padj_primary
              ),
              na.rm=TRUE
            ),

          stringsAsFactors=FALSE
        )
      }
    )
  )

  out <- out[
    order(
      !out$same_direction_all3,
      out$max_padj,
      -out$min_abs_NES
    ),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL
  out
}

conc_primary <- make_concordance(
  all_gsea,
  "HEP_AMBIENT_AWARE"
)

conc_strict <- make_concordance(
  all_gsea,
  "STRICT_REFERENCE_SENSITIVITY"
)

write.csv(
  conc_primary,
  file.path(
    TABDIR,
    "Hallmark_concordance_HEP_AMBIENT_AWARE_v6.8.9.csv"
  ),
  row.names=FALSE
)

write.csv(
  conc_strict,
  file.path(
    TABDIR,
    "Hallmark_concordance_STRICT_REFERENCE_SENSITIVITY_v6.8.9.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Primary-core terminal summaries
# ---------------------------------------------------------

print_primary <- function(
  mode
) {

  x <- all_gsea[
    all_gsea$framework == "PRIMARY_CORE" &
      all_gsea$filter_mode == mode,
    ,
    drop=FALSE
  ]

  x <- x[
    order(
      x$padj,
      -abs(x$NES)
    ),
    ,
    drop=FALSE
  ]

  cat(
    "\n=== PRIMARY_CORE HALLMARK: ",
    mode,
    " ===\n",
    sep=""
  )

  print(
    head(
      x[
        ,
        c(
          "pathway",
          "NES",
          "pval",
          "padj",
          "size",
          "leadingEdge_genes"
        )
      ],
      20
    ),
    row.names=FALSE
  )
}

print_primary(
  "HEP_AMBIENT_AWARE"
)

print_primary(
  "STRICT_REFERENCE_SENSITIVITY"
)

# ---------------------------------------------------------
# Robust same-direction pathways
# ---------------------------------------------------------

cat(
  "\n=== THREE-FRAMEWORK CONCORDANCE: HEP_AMBIENT_AWARE ===\n"
)

print(
  head(
    conc_primary[
      conc_primary$same_direction_all3,
      ,
      drop=FALSE
    ],
    20
  ),
  row.names=FALSE
)

cat(
  "\n=== THREE-FRAMEWORK CONCORDANCE: STRICT_REFERENCE_SENSITIVITY ===\n"
)

print(
  head(
    conc_strict[
      conc_strict$same_direction_all3,
      ,
      drop=FALSE
    ],
    20
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Primary-core GSEA plots
# ---------------------------------------------------------

plot_primary_gsea <- function(
  mode,
  filename
) {

  x <- all_gsea[
    all_gsea$framework == "PRIMARY_CORE" &
      all_gsea$filter_mode == mode,
    ,
    drop=FALSE
  ]

  x <- x[
    order(
      x$padj,
      -abs(x$NES)
    ),
    ,
    drop=FALSE
  ]

  x <- head(
    x,
    20
  )

  x$label <- sub(
    "^HALLMARK_",
    "",
    x$pathway
  )

  x$label <- gsub(
    "_",
    " ",
    x$label
  )

  x$label <- factor(
    x$label,
    levels=rev(x$label)
  )

  x$minus_log10_padj <-
    -log10(
      pmax(
        x$padj,
        1e-300
      )
    )

  p <- ggplot(
    x,
    aes(
      x=NES,
      y=label,
      size=minus_log10_padj,
      fill=NES
    )
  ) +
    geom_point(
      shape=21
    ) +
    scale_fill_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(base_size=10) +
    labs(
      title=paste0(
        "Cholangiocyte PRIMARY_CORE Hallmark GSEA: ",
        mode
      ),
      x="Normalized enrichment score",
      y=NULL,
      size="-log10 FDR",
      fill="NES"
    )

  ggsave(
    filename,
    p,
    width=10,
    height=8
  )
}

plot_primary_gsea(
  "HEP_AMBIENT_AWARE",
  file.path(
    FIGDIR,
    "Cholangiocyte_PRIMARY_CORE_Hallmark_GSEA_HEP_AMBIENT_AWARE_v6.8.9.pdf"
  )
)

plot_primary_gsea(
  "STRICT_REFERENCE_SENSITIVITY",
  file.path(
    FIGDIR,
    "Cholangiocyte_PRIMARY_CORE_Hallmark_GSEA_STRICT_REFERENCE_SENSITIVITY_v6.8.9.pdf"
  )
)

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

summary_lines <- c(
  "# Mouse MASH Cholangiocyte ambient-aware Hallmark GSEA v6.8.9",
  "",
  "- Sham vs Tx biological-sample pseudobulk source: v6.8.8.",
  "- Reference specificity source: v6.8.8.1.",
  "- GSEA ranking statistic: sign(logFC) * sqrt(edgeR QL F).",
  "- No hard DEG threshold is used for GSEA.",
  "- Primary filter: remove Hepatocyte-dominant genes only.",
  "- Strict sensitivity: additionally remove genes with max reference-lineage delta > 2.",
  "- Frameworks: ALL, NO_QCWATCH, PRIMARY_CORE.",
  "- Pathways concordant in direction across all three frameworks are prioritized.",
  "- Sham vs Tx has n=2 biological samples/group; pathway results require cautious interpretation."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_ambient_aware_Hallmark_GSEA_summary_v6.8.9.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.9.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.9 COMPLETE\n")
cat("Hallmark GSEA complete\n")
cat("Primary: Hepatocyte ambient-aware\n")
cat("Sensitivity: strict multi-lineage reference filtering\n")
cat("Frameworks: ALL / NO_QCWATCH / PRIMARY_CORE\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
