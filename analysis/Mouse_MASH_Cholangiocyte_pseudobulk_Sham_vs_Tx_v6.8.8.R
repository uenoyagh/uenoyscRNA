suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

set.seed(20260902)

VERSION <- "v6.8.8"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.7/objects/",
  "Mouse_MASH_Cholangiocyte_state_module_scored_v6.8.7.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte pseudobulk Sham vs Tx\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

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

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

md <- obj@meta.data

cat("=== INPUT ===\n")
cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(counts), "\n")

# ---------------------------------------------------------
# Framework definitions
# ---------------------------------------------------------

frameworks <- list(

  ALL =
    rownames(md),

  NO_QCWATCH =
    rownames(md)[
      !md$Chol_QCwatch_v686
    ],

  PRIMARY_CORE =
    rownames(md)[
      md$Chol_analysis_class_v686 %in%
        c(
          "primary",
          "exploratory_primary"
        )
    ],

  Homeostatic_combined =
    rownames(md)[
      as.character(
        md$Chol_broad_state_v686
      ) == "Homeostatic_like"
    ],

  Inflammatory_reactive =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "Ccl2_Vcam1_inflammatory_reactive"
    ],

  Ductular_like =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "Msln_Aqp5_ductular_like"
    ],

  IEG_stress =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "IEG_stress_response"
    ],

  Krt20_Cdh17_reactive =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "Krt20_Cdh17_reactive_epithelial"
    ],

  Cycling =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "Cycling_cholangiocyte"
    ],

  Ciliated =
    rownames(md)[
      as.character(
        md$Chol_state_v686
      ) ==
        "Ciliated_cholangiocyte"
    ]
)

PB_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

PB_GROUP <- factor(
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

# ---------------------------------------------------------
# Pseudobulk helper
# ---------------------------------------------------------

run_pseudobulk <- function(
  cells,
  framework_name
) {

  cells <- intersect(
    cells,
    colnames(counts)
  )

  sub_md <- md[
    cells,
    ,
    drop=FALSE
  ]

  cell_n <- sapply(
    PB_SAMPLES,
    function(s) {
      sum(
        as.character(
          sub_md$sample
        ) == s
      )
    }
  )

  cat(
    "\n=== FRAMEWORK: ",
    framework_name,
    " ===\n",
    sep=""
  )

  cat("Cells by sample:\n")
  print(cell_n)

  if (any(cell_n < 20)) {

    cat(
      "SKIPPED: one or more samples have <20 cells\n"
    )

    return(NULL)
  }

  pb <- do.call(
    cbind,
    lapply(
      PB_SAMPLES,
      function(s) {

        cells_s <- rownames(sub_md)[
          as.character(
            sub_md$sample
          ) == s
        ]

        Matrix::rowSums(
          counts[
            ,
            cells_s,
            drop=FALSE
          ]
        )
      }
    )
  )

  colnames(pb) <- PB_SAMPLES
  rownames(pb) <- rownames(counts)

  write.csv(
    pb,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_counts_",
        framework_name,
        "_v6.8.8.csv"
      )
    )
  )

  y <- DGEList(
    counts=pb,
    group=PB_GROUP
  )

  design <- model.matrix(
    ~ PB_GROUP
  )

  keep <- filterByExpr(
    y,
    design=design
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(y)

  y <- estimateDisp(
    y,
    design
  )

  fit <- tryCatch(
    glmQLFit(
      y,
      design,
      robust=TRUE
    ),
    error=function(e) {

      message(
        "robust=TRUE failed; retrying robust=FALSE: ",
        conditionMessage(e)
      )

      glmQLFit(
        y,
        design,
        robust=FALSE
      )
    }
  )

  test <- glmQLFTest(
    fit,
    coef=2
  )

  tab <- topTags(
    test,
    n=Inf,
    sort.by="none"
  )$table

  tab$gene <- rownames(tab)

  logcpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  idx <- match(
    tab$gene,
    rownames(logcpm)
  )

  tab$Sham1_logCPM <-
    logcpm[idx, "Sham1"]

  tab$Sham20_logCPM <-
    logcpm[idx, "Sham20"]

  tab$Tx17_logCPM <-
    logcpm[idx, "Tx17"]

  tab$Tx5_logCPM <-
    logcpm[idx, "Tx5"]

  tab$Sham_mean_logCPM <-
    rowMeans(
      tab[
        ,
        c(
          "Sham1_logCPM",
          "Sham20_logCPM"
        )
      ]
    )

  tab$Tx_mean_logCPM <-
    rowMeans(
      tab[
        ,
        c(
          "Tx17_logCPM",
          "Tx5_logCPM"
        )
      ]
    )

  tab$Tx_minus_Sham_logCPM <-
    tab$Tx_mean_logCPM -
    tab$Sham_mean_logCPM

  tab$replicate_pattern <- apply(
    tab[
      ,
      c(
        "Sham1_logCPM",
        "Sham20_logCPM",
        "Tx17_logCPM",
        "Tx5_logCPM"
      )
    ],
    1,
    function(z) {

      sham <- z[1:2]
      tx <- z[3:4]

      if (min(tx) > max(sham)) {
        return("Tx_all_higher")
      }

      if (max(tx) < min(sham)) {
        return("Tx_all_lower")
      }

      "Ranges_overlap"
    }
  )

  tab$framework <- framework_name

  tab <- tab[
    ,
    c(
      "gene",
      "framework",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR",
      "Sham1_logCPM",
      "Sham20_logCPM",
      "Tx17_logCPM",
      "Tx5_logCPM",
      "Sham_mean_logCPM",
      "Tx_mean_logCPM",
      "Tx_minus_Sham_logCPM",
      "replicate_pattern"
    )
  ]

  tab <- tab[
    order(
      tab$FDR,
      -abs(tab$logFC)
    ),
    ,
    drop=FALSE
  ]

  rownames(tab) <- NULL

  write.csv(
    tab,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_DE_Sham_vs_Tx_",
        framework_name,
        "_v6.8.8.csv"
      )
    ),
    row.names=FALSE
  )

  top50 <- head(
    tab,
    50
  )

  write.csv(
    top50,
    file.path(
      TABDIR,
      paste0(
        "Pseudobulk_DE_TOP50_",
        framework_name,
        "_v6.8.8.csv"
      )
    ),
    row.names=FALSE
  )

  data.frame(
    framework=framework_name,
    total_cells=length(cells),

    Sham1_cells=cell_n["Sham1"],
    Sham20_cells=cell_n["Sham20"],
    Tx17_cells=cell_n["Tx17"],
    Tx5_cells=cell_n["Tx5"],

    genes_tested=nrow(tab),

    FDR_lt_0.05=
      sum(
        tab$FDR < 0.05,
        na.rm=TRUE
      ),

    FDR_lt_0.10=
      sum(
        tab$FDR < 0.10,
        na.rm=TRUE
      ),

    FDR_lt_0.10_abs_logFC_ge_1=
      sum(
        tab$FDR < 0.10 &
          abs(tab$logFC) >= 1,
        na.rm=TRUE
      ),

    Tx_all_higher=
      sum(
        tab$replicate_pattern ==
          "Tx_all_higher"
      ),

    Tx_all_lower=
      sum(
        tab$replicate_pattern ==
          "Tx_all_lower"
      ),

    stringsAsFactors=FALSE
  )
}

# ---------------------------------------------------------
# Run all Sham vs Tx frameworks
# ---------------------------------------------------------

framework_summaries <- list()

for (nm in names(frameworks)) {

  ans <- run_pseudobulk(
    frameworks[[nm]],
    nm
  )

  if (!is.null(ans)) {
    framework_summaries[[nm]] <- ans
  }
}

framework_summary <- do.call(
  rbind,
  framework_summaries
)

rownames(framework_summary) <- NULL

write.csv(
  framework_summary,
  file.path(
    TABDIR,
    "Pseudobulk_framework_summary_v6.8.8.csv"
  ),
  row.names=FALSE
)

cat("\n=== PSEUDOBULK FRAMEWORK SUMMARY ===\n")
print(
  framework_summary,
  row.names=FALSE
)

# =========================================================
# Descriptive disease axis: STD vs CDHFD
# NO formal statistics
# =========================================================

run_disease_descriptive <- function(
  cells,
  framework_name
) {

  samples <- c(
    "STD_rep1",
    "CDHFD_rep1"
  )

  cells <- intersect(
    cells,
    colnames(counts)
  )

  sub_md <- md[
    cells,
    ,
    drop=FALSE
  ]

  cell_n <- sapply(
    samples,
    function(s) {
      sum(
        as.character(
          sub_md$sample
        ) == s
      )
    }
  )

  cat(
    "\n=== DESCRIPTIVE DISEASE AXIS: ",
    framework_name,
    " ===\n",
    sep=""
  )

  print(cell_n)

  pb <- do.call(
    cbind,
    lapply(
      samples,
      function(s) {

        cells_s <- rownames(sub_md)[
          as.character(
            sub_md$sample
          ) == s
        ]

        Matrix::rowSums(
          counts[
            ,
            cells_s,
            drop=FALSE
          ]
        )
      }
    )
  )

  colnames(pb) <- samples
  rownames(pb) <- rownames(counts)

  y <- DGEList(
    counts=pb
  )

  y <- calcNormFactors(y)

  lcpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  out <- data.frame(
    gene=rownames(lcpm),

    STD_logCPM=
      lcpm[, "STD_rep1"],

    CDHFD_logCPM=
      lcpm[, "CDHFD_rep1"],

    CDHFD_minus_STD=
      lcpm[, "CDHFD_rep1"] -
      lcpm[, "STD_rep1"],

    framework=framework_name,

    stringsAsFactors=FALSE
  )

  out <- out[
    order(
      -abs(
        out$CDHFD_minus_STD
      )
    ),
    ,
    drop=FALSE
  ]

  rownames(out) <- NULL

  write.csv(
    out,
    file.path(
      TABDIR,
      paste0(
        "Disease_axis_CDHFD_minus_STD_DESCRIPTIVE_",
        framework_name,
        "_v6.8.8.csv"
      )
    ),
    row.names=FALSE
  )

  invisible(out)
}

run_disease_descriptive(
  frameworks$ALL,
  "ALL"
)

run_disease_descriptive(
  frameworks$NO_QCWATCH,
  "NO_QCWATCH"
)

run_disease_descriptive(
  frameworks$PRIMARY_CORE,
  "PRIMARY_CORE"
)

# ---------------------------------------------------------
# Session info and summary
# ---------------------------------------------------------

summary_lines <- c(
  "# Mouse MASH Cholangiocyte pseudobulk Sham vs Tx v6.8.8",
  "",
  "- Biological sample is the statistical unit.",
  "- Sham: Sham1, Sham20.",
  "- Tx: Tx17, Tx5.",
  "- edgeR quasi-likelihood framework.",
  "- FDR is reported but interpreted cautiously because n=2/group.",
  "- Replicate concordance is explicitly reported.",
  "- Multiple cell frameworks are analyzed to distinguish global from state-specific effects.",
  "- QC-watch exclusion is included as a sensitivity analysis.",
  "- STD vs CDHFD is descriptive only; no formal statistical test is performed.",
  "- Hepatocyte ambient-RNA specificity has NOT yet been applied to DE interpretation."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_pseudobulk_Sham_vs_Tx_summary_v6.8.8.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.8.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.8 COMPLETE\n")
cat("Biological-sample pseudobulk complete\n")
cat("Sham vs Tx: n=2/group\n")
cat("STD vs CDHFD: descriptive only\n")
cat("No pathway interpretation performed yet\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
