suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.9.1"

DE_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.8/tables"
)

SPEC_FILE <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.8.1/tables/",
  "Cholangiocyte_reference_specificity_all_genes_v6.8.8.1.csv"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(
  TABDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  FIGDIR,
  recursive=TRUE,
  showWarnings=FALSE
)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte state-specific Hallmark GSEA\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(SPEC_FILE)) {
  stop("Missing specificity file: ", SPEC_FILE)
}

# =========================================================
# State frameworks
# =========================================================

state_files <- c(
  Homeostatic_combined =
    "Pseudobulk_DE_Sham_vs_Tx_Homeostatic_combined_v6.8.8.csv",

  Inflammatory_reactive =
    "Pseudobulk_DE_Sham_vs_Tx_Inflammatory_reactive_v6.8.8.csv",

  Ductular_like =
    "Pseudobulk_DE_Sham_vs_Tx_Ductular_like_v6.8.8.csv",

  Krt20_Cdh17_reactive =
    "Pseudobulk_DE_Sham_vs_Tx_Krt20_Cdh17_reactive_v6.8.8.csv",

  Cycling =
    "Pseudobulk_DE_Sham_vs_Tx_Cycling_v6.8.8.csv",

  Ciliated =
    "Pseudobulk_DE_Sham_vs_Tx_Ciliated_v6.8.8.csv",

  IEG_stress =
    "Pseudobulk_DE_Sham_vs_Tx_IEG_stress_v6.8.8.csv"
)

state_analysis_class <- c(
  Homeostatic_combined="PRIMARY",
  Inflammatory_reactive="PRIMARY",
  Ductular_like="PRIMARY",
  Krt20_Cdh17_reactive="PRIMARY",
  Cycling="PRIMARY",
  Ciliated="PRIMARY",
  IEG_stress="SUPPLEMENTAL_LOW_CELL"
)

primary_states <- names(
  state_analysis_class[
    state_analysis_class == "PRIMARY"
  ]
)

# =========================================================
# Reference specificity table
# =========================================================

spec <- read.csv(
  SPEC_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_spec <- c(
  "gene",
  "Hepatocyte_specificity_class",
  "max_reference_delta",
  "max_reference_lineage",
  "multi_lineage_review_class"
)

missing_spec <- setdiff(
  required_spec,
  colnames(spec)
)

if (length(missing_spec) > 0) {
  stop(
    "Missing specificity columns: ",
    paste(missing_spec, collapse=", ")
  )
}

# =========================================================
# Read and annotate state-specific DE tables
# =========================================================

de_list <- list()

for (state in names(state_files)) {

  path <- file.path(
    DE_DIR,
    state_files[[state]]
  )

  if (!file.exists(path)) {
    stop(
      "Missing DE file for ",
      state,
      ": ",
      path
    )
  }

  de <- read.csv(
    path,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  required_de <- c(
    "gene",
    "logFC",
    "F",
    "PValue",
    "FDR",
    "replicate_pattern"
  )

  missing_de <- setdiff(
    required_de,
    colnames(de)
  )

  if (length(missing_de) > 0) {
    stop(
      "Missing DE columns for ",
      state,
      ": ",
      paste(missing_de, collapse=", ")
    )
  }

  idx <- match(
    de$gene,
    spec$gene
  )

  spec_cols <- setdiff(
    required_spec,
    "gene"
  )

  for (col in spec_cols) {
    de[[col]] <- spec[[col]][idx]
  }

  de$state <- state
  de$analysis_class <- state_analysis_class[[state]]

  de_list[[state]] <- de

  write.csv(
    de,
    file.path(
      TABDIR,
      paste0(
        "State_DE_reference_annotated_",
        state,
        "_v6.8.9.1.csv"
      )
    ),
    row.names=FALSE
  )
}

# =========================================================
# Hallmark gene sets
# =========================================================

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

# =========================================================
# Ranking statistic
#
# Same rule as v6.8.9:
# sign(logFC) * sqrt(edgeR quasi-likelihood F)
# =========================================================

make_rank <- function(x) {

  keep <-
    is.finite(x$logFC) &
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

  if (anyDuplicated(names(stat))) {

    genes <- unique(
      names(stat)
    )

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

  sort(
    stat,
    decreasing=TRUE
  )
}

# =========================================================
# Ambient/reference filters
# =========================================================

filter_de <- function(
  de,
  filter_mode
) {

  if (
    filter_mode ==
      "HEP_AMBIENT_AWARE"
  ) {

    out <- de[
      is.na(
        de$Hepatocyte_specificity_class
      ) |
        de$Hepatocyte_specificity_class !=
          "Hepatocyte_dominant",
      ,
      drop=FALSE
    ]

  } else if (
    filter_mode ==
      "STRICT_REFERENCE_SENSITIVITY"
  ) {

    out <- de[
      (
        is.na(
          de$Hepatocyte_specificity_class
        ) |
          de$Hepatocyte_specificity_class !=
            "Hepatocyte_dominant"
      ) &
        (
          is.na(
            de$max_reference_delta
          ) |
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

  out
}

# =========================================================
# Run state-specific GSEA
# =========================================================

run_state_gsea <- function(
  de,
  state,
  filter_mode
) {

  filtered <- filter_de(
    de,
    filter_mode
  )

  stats <- make_rank(
    filtered
  )

  cat(
    "\n=== STATE GSEA: ",
    state,
    " / ",
    filter_mode,
    " ===\n",
    sep=""
  )

  cat(
    "Analysis class:",
    state_analysis_class[[state]],
    "\n"
  )

  cat(
    "DE genes before filter:",
    nrow(de),
    "\n"
  )

  cat(
    "Ranked genes after filter:",
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

  fg <- as.data.frame(
    fg
  )

  if (nrow(fg) == 0) {
    stop(
      "No GSEA results for ",
      state,
      " / ",
      filter_mode
    )
  }

  fg$state <- state
  fg$analysis_class <-
    state_analysis_class[[state]]

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
        state,
        "_",
        filter_mode,
        "_v6.8.9.1.csv"
      )
    ),
    row.names=FALSE
  )

  fg
}

filter_modes <- c(
  "HEP_AMBIENT_AWARE",
  "STRICT_REFERENCE_SENSITIVITY"
)

gsea_results <- list()

for (state in names(de_list)) {

  for (mode in filter_modes) {

    key <- paste(
      state,
      mode,
      sep="__"
    )

    gsea_results[[key]] <-
      run_state_gsea(
        de=de_list[[state]],
        state=state,
        filter_mode=mode
      )
  }
}

all_gsea <- do.call(
  rbind,
  gsea_results
)

rownames(all_gsea) <- NULL

write.csv(
  all_gsea,
  file.path(
    TABDIR,
    "State_specific_Hallmark_GSEA_ALL_RESULTS_v6.8.9.1.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Terminal: top pathways by state
# =========================================================

print_top_by_state <- function(
  mode,
  n=10
) {

  cat(
    "\n====================================================\n"
  )

  cat(
    "TOP STATE-SPECIFIC HALLMARKS: ",
    mode,
    "\n",
    sep=""
  )

  cat(
    "====================================================\n"
  )

  for (state in names(de_list)) {

    x <- all_gsea[
      all_gsea$state == state &
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
      "\n--- ",
      state,
      " [",
      state_analysis_class[[state]],
      "] ---\n",
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
            "size"
          )
        ],
        n
      ),
      row.names=FALSE
    )
  }
}

print_top_by_state(
  "HEP_AMBIENT_AWARE",
  n=10
)

print_top_by_state(
  "STRICT_REFERENCE_SENSITIVITY",
  n=10
)

# =========================================================
# Selected mechanistic pathways
# =========================================================

selected_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)

make_selected_table <- function(
  mode
) {

  rows <- list()

  for (state in names(de_list)) {

    for (pw in selected_pathways) {

      x <- all_gsea[
        all_gsea$state == state &
          all_gsea$filter_mode == mode &
          all_gsea$pathway == pw,
        ,
        drop=FALSE
      ]

      if (nrow(x) == 1) {

        rows[[
          paste(
            state,
            pw,
            sep="__"
          )
        ]] <- data.frame(
          state=state,
          analysis_class=
            state_analysis_class[[state]],
          pathway=pw,
          NES=x$NES,
          pval=x$pval,
          padj=x$padj,
          stringsAsFactors=FALSE
        )
      }
    }
  }

  out <- do.call(
    rbind,
    rows
  )

  rownames(out) <- NULL
  out
}

selected_primary <- make_selected_table(
  "HEP_AMBIENT_AWARE"
)

selected_strict <- make_selected_table(
  "STRICT_REFERENCE_SENSITIVITY"
)

write.csv(
  selected_primary,
  file.path(
    TABDIR,
    "Selected_pathways_by_state_HEP_AMBIENT_AWARE_v6.8.9.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  selected_strict,
  file.path(
    TABDIR,
    "Selected_pathways_by_state_STRICT_REFERENCE_SENSITIVITY_v6.8.9.1.csv"
  ),
  row.names=FALSE
)

cat(
  "\n=== SELECTED PATHWAYS BY STATE: HEP_AMBIENT_AWARE ===\n"
)

print(
  selected_primary[
    selected_primary$state %in%
      primary_states,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n=== SELECTED PATHWAYS BY STATE: STRICT_REFERENCE_SENSITIVITY ===\n"
)

print(
  selected_strict[
    selected_strict$state %in%
      primary_states,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

# =========================================================
# Cross-state recurrence summary
# PRIMARY states only
# =========================================================

make_recurrence <- function(
  mode
) {

  x <- all_gsea[
    all_gsea$filter_mode == mode &
      all_gsea$state %in%
        primary_states,
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

        data.frame(
          filter_mode=mode,
          pathway=pw,

          states_tested=nrow(z),

          states_NES_positive=
            sum(
              z$NES > 0,
              na.rm=TRUE
            ),

          states_NES_negative=
            sum(
              z$NES < 0,
              na.rm=TRUE
            ),

          states_padj_lt_0.05=
            sum(
              z$padj < 0.05,
              na.rm=TRUE
            ),

          states_positive_padj_lt_0.05=
            sum(
              z$NES > 0 &
                z$padj < 0.05,
              na.rm=TRUE
            ),

          states_negative_padj_lt_0.05=
            sum(
              z$NES < 0 &
                z$padj < 0.05,
              na.rm=TRUE
            ),

          median_NES=
            median(
              z$NES,
              na.rm=TRUE
            ),

          stringsAsFactors=FALSE
        )
      }
    )
  )

  out <- out[
    order(
      -out$states_padj_lt_0.05,
      -abs(out$median_NES)
    ),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL
  out
}

recurrence_primary <- make_recurrence(
  "HEP_AMBIENT_AWARE"
)

recurrence_strict <- make_recurrence(
  "STRICT_REFERENCE_SENSITIVITY"
)

write.csv(
  recurrence_primary,
  file.path(
    TABDIR,
    "Hallmark_cross_state_recurrence_HEP_AMBIENT_AWARE_v6.8.9.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  recurrence_strict,
  file.path(
    TABDIR,
    "Hallmark_cross_state_recurrence_STRICT_REFERENCE_SENSITIVITY_v6.8.9.1.csv"
  ),
  row.names=FALSE
)

cat(
  "\n=== CROSS-STATE RECURRENCE: HEP_AMBIENT_AWARE ===\n"
)

print(
  head(
    recurrence_primary,
    20
  ),
  row.names=FALSE
)

cat(
  "\n=== CROSS-STATE RECURRENCE: STRICT_REFERENCE_SENSITIVITY ===\n"
)

print(
  head(
    recurrence_strict,
    20
  ),
  row.names=FALSE
)

# =========================================================
# NES heatmap helper
# =========================================================

plot_nes_heatmap <- function(
  mode,
  filename
) {

  x <- all_gsea[
    all_gsea$filter_mode == mode,
    ,
    drop=FALSE
  ]

  x <- x[
    x$pathway %in%
      selected_pathways,
    ,
    drop=FALSE
  ]

  state_levels <- c(
    primary_states,
    "IEG_stress"
  )

  pathway_levels <- rev(
    selected_pathways
  )

  x$state <- factor(
    x$state,
    levels=state_levels
  )

  x$pathway_label <- sub(
    "^HALLMARK_",
    "",
    x$pathway
  )

  x$pathway_label <- gsub(
    "_",
    " ",
    x$pathway_label
  )

  pathway_label_levels <- sub(
    "^HALLMARK_",
    "",
    pathway_levels
  )

  pathway_label_levels <- gsub(
    "_",
    " ",
    pathway_label_levels
  )

  x$pathway_label <- factor(
    x$pathway_label,
    levels=pathway_label_levels
  )

  p <- ggplot(
    x,
    aes(
      x=state,
      y=pathway_label,
      fill=NES
    )
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(base_size=9) +
    theme(
      axis.text.x=
        element_text(
          angle=45,
          hjust=1
        )
    ) +
    labs(
      title=paste0(
        "Cholangiocyte state-specific Hallmark NES: ",
        mode
      ),
      x=NULL,
      y=NULL,
      fill="NES"
    )

  ggsave(
    filename,
    p,
    width=12,
    height=8
  )
}

plot_nes_heatmap(
  "HEP_AMBIENT_AWARE",
  file.path(
    FIGDIR,
    "Cholangiocyte_state_specific_Hallmark_NES_HEP_AMBIENT_AWARE_v6.8.9.1.pdf"
  )
)

plot_nes_heatmap(
  "STRICT_REFERENCE_SENSITIVITY",
  file.path(
    FIGDIR,
    "Cholangiocyte_state_specific_Hallmark_NES_STRICT_REFERENCE_SENSITIVITY_v6.8.9.1.pdf"
  )
)

# =========================================================
# Summary
# =========================================================

summary_lines <- c(
  "# Mouse MASH Cholangiocyte state-specific Hallmark GSEA v6.8.9.1",
  "",
  "- Source state-specific pseudobulk DE: v6.8.8.",
  "- Reference specificity source: v6.8.8.1.",
  "- Ranking: sign(logFC) * sqrt(edgeR quasi-likelihood F).",
  "- Primary states: Homeostatic_combined, Inflammatory_reactive, Ductular_like, Krt20_Cdh17_reactive, Cycling, Ciliated.",
  "- IEG_stress is calculated as supplemental because each biological sample contains relatively few cells.",
  "- Primary ambient filter removes Hepatocyte-dominant genes.",
  "- Strict sensitivity additionally removes genes with max reference-lineage delta > 2.",
  "- No cells, clusters, or annotations are changed.",
  "- Sham vs Tx remains n=2 biological samples per group.",
  "- State-specific GSEA is used to distinguish within-state transcriptional reprogramming from changes in state composition."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_state_specific_Hallmark_GSEA_summary_v6.8.9.1.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.9.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.9.1 COMPLETE\n")
cat("State-specific Hallmark validation complete\n")
cat("Primary states:", length(primary_states), "\n")
cat("IEG_stress: supplemental low-cell analysis\n")
cat("No cells or annotations changed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
