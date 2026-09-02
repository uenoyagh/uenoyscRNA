suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.8.7"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.6/objects/",
  "Mouse_MASH_Cholangiocyte_annotation_frozen_v6.8.6.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte state composition + module axis\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

required_md <- c(
  "sample",
  "condition",
  "Chol_state_v686",
  "Chol_broad_state_v686",
  "Chol_analysis_class_v686",
  "Chol_QCwatch_v686"
)

missing_md <- setdiff(
  required_md,
  colnames(obj@meta.data)
)

if (length(missing_md) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing_md, collapse=", ")
  )
}

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

condition_order <- c(
  "STD",
  "CDHFD",
  "Sham",
  "Tx"
)

obj$sample <- factor(
  as.character(obj$sample),
  levels=sample_order
)

obj$condition <- factor(
  as.character(obj$condition),
  levels=condition_order
)

cat("=== INPUT ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")

cat("\n=== CELLS BY SAMPLE ===\n")
print(table(obj$sample))

# =========================================================
# PART 1
# Frozen-state composition
# =========================================================

state_levels <- levels(
  obj$Chol_state_v686
)

if (is.null(state_levels)) {
  state_levels <- unique(
    as.character(obj$Chol_state_v686)
  )
}

composition_table <- function(
  x,
  sensitivity_name
) {

  tab <- table(
    state=factor(
      as.character(x$Chol_state_v686),
      levels=state_levels
    ),
    sample=factor(
      as.character(x$sample),
      levels=sample_order
    )
  )

  frac <- prop.table(
    tab,
    margin=2
  )

  count_file <- file.path(
    TABDIR,
    paste0(
      "Cholangiocyte_state_counts_",
      sensitivity_name,
      "_v6.8.7.csv"
    )
  )

  frac_file <- file.path(
    TABDIR,
    paste0(
      "Cholangiocyte_state_fraction_",
      sensitivity_name,
      "_v6.8.7.csv"
    )
  )

  write.csv(
    as.data.frame.matrix(tab),
    count_file
  )

  write.csv(
    as.data.frame.matrix(frac),
    frac_file
  )

  list(
    counts=tab,
    fractions=frac
  )
}

# All frozen cells
comp_all <- composition_table(
  obj@meta.data,
  "ALL"
)

# Sensitivity excluding QC-watch state
md_noqc <- obj@meta.data[
  !obj$Chol_QCwatch_v686,
  ,
  drop=FALSE
]

comp_noqc <- composition_table(
  md_noqc,
  "NO_QCWATCH"
)

# ---------------------------------------------------------
# Axis summary helper
# ---------------------------------------------------------

composition_axis_summary <- function(
  frac,
  analysis_name
) {

  needed <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  if (!all(needed %in% colnames(frac))) {
    stop("Missing sample columns in fraction table.")
  }

  out <- data.frame(
    state=rownames(frac),

    STD=
      as.numeric(frac[, "STD_rep1"]),

    CDHFD=
      as.numeric(frac[, "CDHFD_rep1"]),

    disease_delta_CDHFD_minus_STD=
      as.numeric(
        frac[, "CDHFD_rep1"] -
        frac[, "STD_rep1"]
      ),

    Sham1=
      as.numeric(frac[, "Sham1"]),

    Sham20=
      as.numeric(frac[, "Sham20"]),

    Tx17=
      as.numeric(frac[, "Tx17"]),

    Tx5=
      as.numeric(frac[, "Tx5"]),

    stringsAsFactors=FALSE
  )

  out$Sham_mean <-
    rowMeans(
      out[, c("Sham1", "Sham20")]
    )

  out$Tx_mean <-
    rowMeans(
      out[, c("Tx17", "Tx5")]
    )

  out$Tx_minus_Sham <-
    out$Tx_mean -
    out$Sham_mean

  out$Tx_vs_Sham_pattern <- apply(
    out,
    1,
    function(z) {

      sham <- as.numeric(
        z[c("Sham1", "Sham20")]
      )

      tx <- as.numeric(
        z[c("Tx17", "Tx5")]
      )

      if (min(tx) > max(sham)) {
        return("Tx_all_higher")
      }

      if (max(tx) < min(sham)) {
        return("Tx_all_lower")
      }

      "Ranges_overlap"
    }
  )

  out$analysis <- analysis_name

  out <- out[
    order(
      -abs(out$Tx_minus_Sham)
    ),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL

  out
}

axis_all <- composition_axis_summary(
  comp_all$fractions,
  "ALL"
)

axis_noqc <- composition_axis_summary(
  comp_noqc$fractions,
  "NO_QCWATCH"
)

write.csv(
  axis_all,
  file.path(
    TABDIR,
    "Cholangiocyte_state_axis_summary_ALL_v6.8.7.csv"
  ),
  row.names=FALSE
)

write.csv(
  axis_noqc,
  file.path(
    TABDIR,
    "Cholangiocyte_state_axis_summary_NO_QCWATCH_v6.8.7.csv"
  ),
  row.names=FALSE
)

cat("\n=== STATE AXIS: ALL ===\n")
print(
  axis_all[
    ,
    c(
      "state",
      "STD",
      "CDHFD",
      "disease_delta_CDHFD_minus_STD",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "Tx_vs_Sham_pattern"
    )
  ],
  row.names=FALSE
)

cat("\n=== STATE AXIS: NO QC-WATCH ===\n")
print(
  axis_noqc[
    ,
    c(
      "state",
      "STD",
      "CDHFD",
      "disease_delta_CDHFD_minus_STD",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "Tx_vs_Sham_pattern"
    )
  ],
  row.names=FALSE
)

# =========================================================
# Composition heatmaps
# =========================================================

plot_fraction_heatmap <- function(
  frac,
  title_text,
  filename
) {

  df <- as.data.frame(
    as.table(frac)
  )

  colnames(df) <- c(
    "state",
    "sample",
    "fraction"
  )

  df$sample <- factor(
    as.character(df$sample),
    levels=sample_order
  )

  df$state <- factor(
    as.character(df$state),
    levels=rev(state_levels)
  )

  p <- ggplot(
    df,
    aes(
      x=sample,
      y=state,
      fill=fraction
    )
  ) +
    geom_tile() +
    scale_fill_gradient(
      low="white",
      high="#FF1A1A"
    ) +
    theme_classic(base_size=10) +
    theme(
      axis.text.x=
        element_text(
          angle=45,
          hjust=1
        )
    ) +
    labs(
      title=title_text,
      x=NULL,
      y=NULL,
      fill="Fraction"
    )

  ggsave(
    filename,
    p,
    width=10,
    height=8
  )
}

plot_fraction_heatmap(
  comp_all$fractions,
  "Cholangiocyte frozen-state fractions: all cells",
  file.path(
    FIGDIR,
    "Cholangiocyte_state_fraction_heatmap_ALL_v6.8.7.pdf"
  )
)

plot_fraction_heatmap(
  comp_noqc$fractions,
  "Cholangiocyte frozen-state fractions: QC-watch excluded",
  file.path(
    FIGDIR,
    "Cholangiocyte_state_fraction_heatmap_NO_QCWATCH_v6.8.7.pdf"
  )
)

# =========================================================
# PART 2
# Prespecified biological modules
# =========================================================

module_sets <- list(

  Biliary_identity = c(
    "Krt19","Krt7","Krt8","Krt18",
    "Epcam","Sox9","Muc1",
    "Hnf1b","Cftr","Slc4a2"
  ),

  Ductular_reactive = c(
    "Spp1","Mmp7","Krt23",
    "Tacstd2","Prom1","Klf5"
  ),

  Inflammatory_adhesion = c(
    "Icam1","Vcam1",
    "Cxcl1","Cxcl2",
    "Ccl2","Il6",
    "Nfkbia","Socs3"
  ),

  Remodeling_fibrogenic = c(
    "Tgfb1","Tgfb2",
    "Ccn2","Jag1",
    "Serpine1","Thbs1",
    "Inhba"
  ),

  Stress_IEG = c(
    "Cdkn1a","Cdkn2a",
    "Fos","Jun","Atf3",
    "Ddit4","Gadd45a"
  ),

  Cycling = c(
    "Mki67","Top2a",
    "Birc5","Ube2c",
    "Cenpf","Pcna","Stmn1"
  ),

  Krt20_Cdh17_reactive = c(
    "Krt20","Cdh17",
    "Inhba","Fst"
  ),

  Dmbt1_Duox2_reactive = c(
    "Dmbt1","Duox2",
    "Duoxa2","Slc26a9",
    "Gcnt3","Ern2","Tns4"
  ),

  Ciliated = c(
    "Cfap73","Cfap44",
    "Lrrc23","Ttc21a",
    "Drc1","Iqub",
    "Dnali1","Mns1"
  ),

  Tuft_like = c(
    "Pou2f3","Trpm5",
    "Gnat3"
  )
)

module_sets_present <- lapply(
  module_sets,
  function(x) {
    intersect(
      x,
      rownames(obj)
    )
  }
)

cat("\n=== MODULE AVAILABILITY ===\n")

for (nm in names(module_sets)) {

  cat(
    nm,
    ": ",
    length(module_sets_present[[nm]]),
    "/",
    length(module_sets[[nm]]),
    " genes",
    "\n",
    sep=""
  )
}

# ---------------------------------------------------------
# Re-normalize RNA for transparent common scoring
# ---------------------------------------------------------

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

# ---------------------------------------------------------
# AddModuleScore
# ---------------------------------------------------------

module_columns <- character()

for (nm in names(module_sets_present)) {

  genes <- module_sets_present[[nm]]

  if (length(genes) < 2) {
    warning(
      "Skipping module with <2 genes: ",
      nm
    )
    next
  }

  temp_name <- paste0(
    "TEMP_",
    nm,
    "_"
  )

  obj <- AddModuleScore(
    obj,
    features=list(genes),
    assay="RNA",
    name=temp_name,
    seed=20260902
  )

  old_col <- paste0(
    temp_name,
    "1"
  )

  new_col <- paste0(
    "CholMOD_",
    nm,
    "_v687"
  )

  obj@meta.data[[new_col]] <-
    obj@meta.data[[old_col]]

  obj@meta.data[[old_col]] <- NULL

  module_columns[nm] <- new_col
}

# =========================================================
# Sample-level module summary
# =========================================================

summarize_modules <- function(
  md,
  analysis_name
) {

  out_list <- list()

  for (s in sample_order) {

    x <- md[
      as.character(md$sample) == s,
      ,
      drop=FALSE
    ]

    if (nrow(x) == 0) {
      next
    }

    for (nm in names(module_columns)) {

      col <- module_columns[[nm]]

      out_list[[
        paste(s, nm, sep="__")
      ]] <- data.frame(
        analysis=analysis_name,
        sample=s,
        condition=
          as.character(x$condition[1]),
        module=nm,
        n_cells=nrow(x),
        mean_score=
          mean(
            x[[col]],
            na.rm=TRUE
          ),
        median_score=
          median(
            x[[col]],
            na.rm=TRUE
          ),
        stringsAsFactors=FALSE
      )
    }
  }

  out <- do.call(
    rbind,
    out_list
  )

  rownames(out) <- NULL
  out
}

module_sample_all <- summarize_modules(
  obj@meta.data,
  "ALL"
)

module_sample_noqc <- summarize_modules(
  obj@meta.data[
    !obj$Chol_QCwatch_v686,
    ,
    drop=FALSE
  ],
  "NO_QCWATCH"
)

write.csv(
  module_sample_all,
  file.path(
    TABDIR,
    "Cholangiocyte_module_scores_by_sample_ALL_v6.8.7.csv"
  ),
  row.names=FALSE
)

write.csv(
  module_sample_noqc,
  file.path(
    TABDIR,
    "Cholangiocyte_module_scores_by_sample_NO_QCWATCH_v6.8.7.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Module-axis helper
# =========================================================

module_axis_summary <- function(
  sample_df,
  analysis_name
) {

  modules <- unique(
    sample_df$module
  )

  out <- do.call(
    rbind,
    lapply(
      modules,
      function(m) {

        x <- sample_df[
          sample_df$module == m,
          ,
          drop=FALSE
        ]

        get_score <- function(s) {

          z <- x$mean_score[
            x$sample == s
          ]

          if (length(z) != 1) {
            stop(
              "Expected one score for ",
              m,
              " / ",
              s
            )
          }

          z
        }

        STD <- get_score("STD_rep1")
        CDHFD <- get_score("CDHFD_rep1")
        Sham1 <- get_score("Sham1")
        Sham20 <- get_score("Sham20")
        Tx17 <- get_score("Tx17")
        Tx5 <- get_score("Tx5")

        sham <- c(Sham1, Sham20)
        tx <- c(Tx17, Tx5)

        pattern <- if (
          min(tx) > max(sham)
        ) {
          "Tx_all_higher"
        } else if (
          max(tx) < min(sham)
        ) {
          "Tx_all_lower"
        } else {
          "Ranges_overlap"
        }

        data.frame(
          analysis=analysis_name,
          module=m,
          STD=STD,
          CDHFD=CDHFD,
          disease_delta_CDHFD_minus_STD=
            CDHFD - STD,
          Sham1=Sham1,
          Sham20=Sham20,
          Sham_mean=mean(sham),
          Tx17=Tx17,
          Tx5=Tx5,
          Tx_mean=mean(tx),
          Tx_minus_Sham=
            mean(tx) - mean(sham),
          Tx_vs_Sham_pattern=pattern,
          stringsAsFactors=FALSE
        )
      }
    )
  )

  out <- out[
    order(
      -abs(out$Tx_minus_Sham)
    ),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL
  out
}

module_axis_all <- module_axis_summary(
  module_sample_all,
  "ALL"
)

module_axis_noqc <- module_axis_summary(
  module_sample_noqc,
  "NO_QCWATCH"
)

write.csv(
  module_axis_all,
  file.path(
    TABDIR,
    "Cholangiocyte_module_axis_ALL_v6.8.7.csv"
  ),
  row.names=FALSE
)

write.csv(
  module_axis_noqc,
  file.path(
    TABDIR,
    "Cholangiocyte_module_axis_NO_QCWATCH_v6.8.7.csv"
  ),
  row.names=FALSE
)

cat("\n=== MODULE AXIS: ALL ===\n")

print(
  module_axis_all[
    ,
    c(
      "module",
      "STD",
      "CDHFD",
      "disease_delta_CDHFD_minus_STD",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "Tx_vs_Sham_pattern"
    )
  ],
  row.names=FALSE
)

cat("\n=== MODULE AXIS: NO QC-WATCH ===\n")

print(
  module_axis_noqc[
    ,
    c(
      "module",
      "STD",
      "CDHFD",
      "disease_delta_CDHFD_minus_STD",
      "Sham_mean",
      "Tx_mean",
      "Tx_minus_Sham",
      "Tx_vs_Sham_pattern"
    )
  ],
  row.names=FALSE
)

# =========================================================
# Module score heatmaps
# =========================================================

plot_module_heatmap <- function(
  df,
  title_text,
  filename
) {

  df$sample <- factor(
    df$sample,
    levels=sample_order
  )

  df$module <- factor(
    df$module,
    levels=rev(
      names(module_sets)
    )
  )

  p <- ggplot(
    df,
    aes(
      x=sample,
      y=module,
      fill=mean_score
    )
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low="#0033FF",
      mid="#FFFFFF",
      high="#FF1A1A",
      midpoint=0
    ) +
    theme_classic(base_size=10) +
    theme(
      axis.text.x=
        element_text(
          angle=45,
          hjust=1
        )
    ) +
    labs(
      title=title_text,
      x=NULL,
      y=NULL,
      fill="Module score"
    )

  ggsave(
    filename,
    p,
    width=10,
    height=7
  )
}

plot_module_heatmap(
  module_sample_all,
  "Cholangiocyte module scores by biological sample: all cells",
  file.path(
    FIGDIR,
    "Cholangiocyte_module_heatmap_ALL_v6.8.7.pdf"
  )
)

plot_module_heatmap(
  module_sample_noqc,
  "Cholangiocyte module scores by biological sample: QC-watch excluded",
  file.path(
    FIGDIR,
    "Cholangiocyte_module_heatmap_NO_QCWATCH_v6.8.7.pdf"
  )
)

# =========================================================
# Tx-minus-Sham overview
# =========================================================

plot_axis_delta <- function(
  df,
  label_col,
  value_col,
  title_text,
  filename
) {

  plot_df <- data.frame(
    label=df[[label_col]],
    delta=df[[value_col]],
    stringsAsFactors=FALSE
  )

  plot_df <- plot_df[
    order(plot_df$delta),
    ,
    drop=FALSE
  ]

  plot_df$label <- factor(
    plot_df$label,
    levels=plot_df$label
  )

  p <- ggplot(
    plot_df,
    aes(
      x=label,
      y=delta
    )
  ) +
    geom_col() +
    geom_hline(
      yintercept=0,
      linetype=2
    ) +
    coord_flip() +
    theme_classic(base_size=10) +
    labs(
      title=title_text,
      x=NULL,
      y="Tx mean - Sham mean"
    )

  ggsave(
    filename,
    p,
    width=9,
    height=7
  )
}

plot_axis_delta(
  axis_all,
  "state",
  "Tx_minus_Sham",
  "Cholangiocyte state composition: Tx - Sham",
  file.path(
    FIGDIR,
    "Cholangiocyte_state_Tx_minus_Sham_ALL_v6.8.7.pdf"
  )
)

plot_axis_delta(
  module_axis_all,
  "module",
  "Tx_minus_Sham",
  "Cholangiocyte module scores: Tx - Sham",
  file.path(
    FIGDIR,
    "Cholangiocyte_module_Tx_minus_Sham_ALL_v6.8.7.pdf"
  )
)

# =========================================================
# Save scored object
# =========================================================

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_state_module_scored_v6.8.7.rds"
  ),
  compress=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte state composition + module axis v6.8.7",
  "",
  paste0("- Cells: ", ncol(obj)),
  "- Frozen annotation source: v6.8.6",
  "- State composition analyzed using biological samples.",
  "- Primary view includes all frozen states.",
  "- Sensitivity view excludes Disease_enriched_high_complexity_QC_watch.",
  "- STD vs CDHFD is descriptive only because n=1 per condition.",
  "- Sham vs Tx uses biological-sample means with n=2 per group.",
  "- Tx/Sham outputs emphasize effect size and between-sample range rather than formal significance.",
  "- Module scores are calculated from prespecified gene sets using AddModuleScore."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_state_composition_module_axis_summary_v6.8.7.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.7.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.7 COMPLETE\n")
cat("Cells:", ncol(obj), "\n")
cat("State composition: ALL + NO_QCWATCH\n")
cat("Module analysis: ALL + NO_QCWATCH\n")
cat("STD/CDHFD: descriptive only\n")
cat("Sham/Tx: biological-sample effect-size analysis\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
