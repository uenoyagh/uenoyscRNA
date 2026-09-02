suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

VERSION <- "v6.9.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.5/objects/",
  "Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte pseudobulk Sham vs Tx\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

cat("Input cells:", ncol(obj), "\n")
cat("Input features:", nrow(obj), "\n\n")

# =========================================================
# Sample column
# =========================================================

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

sample_vec <- as.character(
  obj@meta.data[[sample_col]]
)

target_samples <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

missing_samples <- setdiff(
  target_samples,
  unique(sample_vec)
)

if (length(missing_samples) > 0) {
  stop(
    "Missing samples: ",
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
# Frameworks
# =========================================================

CLASS_COL <- "Monocyte_analysis_class_v6.9.4"

if (!(CLASS_COL %in% colnames(obj@meta.data))) {
  stop("Missing analysis class column.")
}

analysis_class <- as.character(
  obj@meta.data[[CLASS_COL]]
)

framework_cells <- list(

  ALL=colnames(obj),

  NO_QCWATCH=colnames(obj)[
    analysis_class != "QC_watch_sensitivity"
  ],

  PRIMARY_CORE=colnames(obj)[
    analysis_class %in% c(
      "primary",
      "disease_enriched_primary"
    )
  ]
)

# Restrict to Sham/Tx
framework_cells <- lapply(
  framework_cells,
  function(cells) {

    cells[
      sample_vec[
        match(cells, colnames(obj))
      ] %in% target_samples
    ]
  }
)

framework_cell_summary <- do.call(
  rbind,
  lapply(
    names(framework_cells),
    function(fw) {

      cells <- framework_cells[[fw]]

      do.call(
        rbind,
        lapply(
          target_samples,
          function(s) {

            n <- sum(
              sample_vec[
                match(cells, colnames(obj))
              ] == s
            )

            data.frame(
              framework=fw,
              sample=s,
              n_cells=n,
              stringsAsFactors=FALSE
            )
          }
        )
      )
    }
  )
)

write.csv(
  framework_cell_summary,
  file.path(
    TABDIR,
    "Monocyte_pseudobulk_framework_cell_counts_v6.9.6.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Pseudobulk helper
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
    match(cells, colnames(obj))
  ]

  for (s in sample_order) {

    scells <- cells[
      cell_samples == s
    ]

    if (length(scells) == 0) {
      stop(
        "No cells for sample ",
        s
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

# =========================================================
# edgeR helper
# =========================================================

run_edger <- function(
  pb,
  framework_name
) {

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

  tt <- tt[
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

  # -------------------------------------------------------
  # Sample-level logCPM
  # -------------------------------------------------------

  lcpm <- cpm(
    y,
    log=TRUE,
    prior.count=0.5
  )

  lcpm_df <- data.frame(
    gene=rownames(lcpm),
    Sham1=lcpm[, "Sham1"],
    Sham20=lcpm[, "Sham20"],
    Tx17=lcpm[, "Tx17"],
    Tx5=lcpm[, "Tx5"],
    stringsAsFactors=FALSE
  )

  # replicate direction
  lcpm_df$Tx17_minus_Sham_mean <-
    lcpm_df$Tx17 -
    rowMeans(
      lcpm_df[
        ,
        c(
          "Sham1",
          "Sham20"
        )
      ]
    )

  lcpm_df$Tx5_minus_Sham_mean <-
    lcpm_df$Tx5 -
    rowMeans(
      lcpm_df[
        ,
        c(
          "Sham1",
          "Sham20"
        )
      ]
    )

  lcpm_df$replicate_direction <-
    ifelse(
      lcpm_df$Tx17_minus_Sham_mean > 0 &
      lcpm_df$Tx5_minus_Sham_mean > 0,
      "Tx_both_higher",
      ifelse(
        lcpm_df$Tx17_minus_Sham_mean < 0 &
        lcpm_df$Tx5_minus_Sham_mean < 0,
        "Tx_both_lower",
        "discordant"
      )
    )

  tt <- merge(
    tt,
    lcpm_df,
    by="gene",
    all.x=TRUE,
    sort=FALSE
  )

  tt <- tt[
    order(
      tt$PValue
    ),
    ,
    drop=FALSE
  ]

  # Primary reporting class
  tt$report_class <- "not_selected"

  tt$report_class[
    tt$FDR < 0.10 &
    abs(tt$logFC) >= 0.5 &
    tt$replicate_direction != "discordant"
  ] <- "FDR_lt0.10_concordant"

  tt$report_class[
    tt$report_class == "not_selected" &
    tt$FDR < 0.25 &
    abs(tt$logFC) >= 0.75 &
    tt$replicate_direction != "discordant"
  ] <- "exploratory_FDR_lt0.25_concordant"

  write.csv(
    tt,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_pseudobulk_",
        framework_name,
        "_Sham_vs_Tx_v6.9.6.csv"
      )
    ),
    row.names=FALSE
  )

  selected <- tt[
    tt$report_class != "not_selected",
    ,
    drop=FALSE
  ]

  write.csv(
    selected,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_pseudobulk_",
        framework_name,
        "_selected_v6.9.6.csv"
      )
    ),
    row.names=FALSE
  )

  list(
    full=tt,
    selected=selected,
    n_tested=nrow(tt)
  )
}

# =========================================================
# Run all frameworks
# =========================================================

results <- list()

for (fw in names(framework_cells)) {

  cat(
    "\nRunning framework:",
    fw,
    "\n"
  )

  cells <- framework_cells[[fw]]

  pb <- make_pseudobulk(
    mat=counts,
    cells=cells,
    sample_vec=sample_vec,
    sample_order=target_samples
  )

  write.csv(
    pb,
    file.path(
      TABDIR,
      paste0(
        "Monocyte_pseudobulk_raw_counts_",
        fw,
        "_v6.9.6.csv"
      )
    )
  )

  results[[fw]] <- run_edger(
    pb=pb,
    framework_name=fw
  )
}

# =========================================================
# Cross-framework concordance
# =========================================================

top_union <- unique(
  unlist(
    lapply(
      results,
      function(x) {
        head(
          x$full$gene,
          50
        )
      }
    )
  )
)

cross_rows <- list()

for (gene in top_union) {

  row <- data.frame(
    gene=gene,
    stringsAsFactors=FALSE
  )

  for (fw in names(results)) {

    tt <- results[[fw]]$full

    x <- tt[
      tt$gene == gene,
      ,
      drop=FALSE
    ]

    if (nrow(x) == 1) {

      row[[paste0(fw, "_logFC")]] <-
        x$logFC

      row[[paste0(fw, "_FDR")]] <-
        x$FDR

      row[[paste0(fw, "_direction")]] <-
        x$replicate_direction

    } else {

      row[[paste0(fw, "_logFC")]] <-
        NA_real_

      row[[paste0(fw, "_FDR")]] <-
        NA_real_

      row[[paste0(fw, "_direction")]] <-
        NA_character_
    }
  }

  cross_rows[[gene]] <- row
}

cross_framework <- do.call(
  rbind,
  cross_rows
)

rownames(cross_framework) <- NULL

write.csv(
  cross_framework,
  file.path(
    TABDIR,
    "Monocyte_pseudobulk_cross_framework_top_union_v6.9.6.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Terminal summary
# =========================================================

cat("\n=== PSEUDOBULK CELL COUNTS ===\n")
print(
  framework_cell_summary,
  row.names=FALSE
)

for (fw in names(results)) {

  cat(
    "\n=== ",
    fw,
    " TOP20 ===\n",
    sep=""
  )

  print(
    head(
      results[[fw]]$full[
        ,
        c(
          "gene",
          "logFC",
          "FDR",
          "replicate_direction",
          "report_class"
        )
      ],
      20
    ),
    row.names=FALSE
  )

  cat(
    "\n=== ",
    fw,
    " SELECTED ===\n",
    sep=""
  )

  sel <- results[[fw]]$selected

  if (nrow(sel) == 0) {

    cat("NONE\n")

  } else {

    print(
      sel[
        ,
        c(
          "gene",
          "logFC",
          "FDR",
          "replicate_direction",
          "report_class"
        )
      ],
      row.names=FALSE
    )
  }
}

cat("\n====================================================\n")
cat("v6.9.6 COMPLETE\n")
cat("Monocyte pseudobulk Sham vs Tx complete\n")
cat("Biological samples: Sham1, Sham20, Tx17, Tx5\n")
cat("n=2/group\n")
cat("Frameworks: ALL / NO_QCWATCH / PRIMARY_CORE\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.6.txt"
  )
)
