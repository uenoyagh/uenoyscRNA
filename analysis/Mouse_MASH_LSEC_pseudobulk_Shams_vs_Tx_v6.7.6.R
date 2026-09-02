suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(Matrix)
  library(ggplot2)
})

VERSION <- "v6.7.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.5/objects/",
  "Mouse_MASH_LSEC_annotated_v6.7.5.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(INPUT_RDS)

required <- c(
  "sample",
  "condition",
  "LSEC_state_v675",
  "LSEC_QC_flag_v675"
)

miss <- setdiff(required, colnames(obj@meta.data))
if (length(miss) > 0) {
  stop("Missing metadata: ", paste(miss, collapse=", "))
}

DefaultAssay(obj) <- "RNA"

counts <- GetAssayData(
  obj,
  assay="RNA",
  layer="counts"
)

# -------------------------------------------------------
# Sham / Tx only
# -------------------------------------------------------

keep_st <- colnames(obj)[
  obj$condition %in% c("Sham","Tx")
]

obj_st <- subset(
  obj,
  cells=keep_st
)

cat("============================================\n")
cat("Mouse MASH LSEC pseudobulk\n")
cat("Version:", VERSION, "\n")
cat("Contrast: Tx vs Sham\n")
cat("============================================\n\n")

cat("Cells by sample:\n")
print(table(obj_st$sample))

# -------------------------------------------------------
# Analysis definitions
# -------------------------------------------------------

analyses <- list(

  All_9190_framework = unique(
    as.character(obj_st$LSEC_state_v675)
  ),

  Primary_no_QC = setdiff(
    unique(as.character(obj_st$LSEC_state_v675)),
    "Low_quality_ambient_enriched_LSEC"
  ),

  Shared_core = c(
    "Inflammatory_stress_high_LSEC",
    "Homeostatic_like_LSEC",
    "Wnt_angiocrine_high_LSEC",
    "Cycling_LSEC"
  ),

  Inflammatory_stress_high = c(
    "Inflammatory_stress_high_LSEC"
  ),

  Homeostatic_like = c(
    "Homeostatic_like_LSEC"
  ),

  Wnt_angiocrine_high = c(
    "Wnt_angiocrine_high_LSEC"
  )
)

# -------------------------------------------------------
# Function
# -------------------------------------------------------

run_pb <- function(label, states) {

  cat("\n====================================\n")
  cat("ANALYSIS:", label, "\n")
  cat("====================================\n")

  cells <- colnames(obj_st)[
    obj_st$LSEC_state_v675 %in% states
  ]

  md <- obj_st@meta.data[
    cells,
    ,
    drop=FALSE
  ]

  sample_order <- c(
    "Sham1","Sham20","Tx17","Tx5"
  )

  sample_order <- sample_order[
    sample_order %in% unique(md$sample)
  ]

  cell_table <- table(
    factor(
      md$sample,
      levels=sample_order
    )
  )

  cat("Cell numbers:\n")
  print(cell_table)

  write.csv(
    data.frame(
      sample=names(cell_table),
      cells=as.integer(cell_table)
    ),
    file.path(
      TABDIR,
      paste0(label, "_cell_counts_v6.7.6.csv")
    ),
    row.names=FALSE
  )

  if (any(cell_table == 0)) {
    warning(label, ": zero-cell sample; skipping.")
    return(NULL)
  }

  cnt <- counts[
    ,
    cells,
    drop=FALSE
  ]

  pb <- sapply(
    sample_order,
    function(s) {

      cc <- cells[
        md$sample == s
      ]

      Matrix::rowSums(
        counts[
          ,
          cc,
          drop=FALSE
        ]
      )
    }
  )

  colnames(pb) <- sample_order

  group <- factor(
    c(
      rep("Sham", sum(grepl("^Sham", sample_order))),
      rep("Tx", sum(grepl("^Tx", sample_order)))
    ),
    levels=c("Sham","Tx")
  )

  # Safer: derive directly from metadata
  group <- factor(
    sapply(
      sample_order,
      function(s) {
        unique(
          as.character(
            md$condition[
              md$sample == s
            ]
          )
        )[1]
      }
    ),
    levels=c("Sham","Tx")
  )

  y <- DGEList(
    counts=pb,
    group=group
  )

  keep_gene <- filterByExpr(
    y,
    group=group
  )

  y <- y[
    keep_gene,
    ,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(y)

  design <- model.matrix(
    ~0 + group
  )

  colnames(design) <- levels(group)

  y <- estimateDisp(
    y,
    design,
    robust=TRUE
  )

  fit <- glmQLFit(
    y,
    design,
    robust=TRUE
  )

  contrast <- makeContrasts(
    Tx - Sham,
    levels=design
  )

  qlf <- glmQLFTest(
    fit,
    contrast=contrast
  )

  res <- topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table

  res$gene <- rownames(res)
  rownames(res) <- NULL

  # -----------------------------------------------
  # replicate-level normalized expression
  # -----------------------------------------------

  logcpm <- cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  common_genes <- intersect(
    res$gene,
    rownames(logcpm)
  )

  res <- res[
    match(common_genes, res$gene),
    ,
    drop=FALSE
  ]

  lc <- logcpm[
    common_genes,
    ,
    drop=FALSE
  ]

  res$Sham1_logCPM <- lc[, "Sham1"]
  res$Sham20_logCPM <- lc[, "Sham20"]
  res$Tx17_logCPM <- lc[, "Tx17"]
  res$Tx5_logCPM <- lc[, "Tx5"]

  # strict replicate concordance
  res$Tx_both_above_Shams <-
    pmin(
      res$Tx17_logCPM,
      res$Tx5_logCPM
    ) >
    pmax(
      res$Sham1_logCPM,
      res$Sham20_logCPM
    )

  res$Tx_both_below_Shams <-
    pmax(
      res$Tx17_logCPM,
      res$Tx5_logCPM
    ) <
    pmin(
      res$Sham1_logCPM,
      res$Sham20_logCPM
    )

  res$strict_direction_consistent <-
    ifelse(
      res$Tx_both_above_Shams,
      "Tx_up_both",
      ifelse(
        res$Tx_both_below_Shams,
        "Tx_down_both",
        "Not_strictly_consistent"
      )
    )

  res$effect_ge_0.5 <-
    abs(res$logFC) >= 0.5

  res$FDR_lt_0.10 <-
    res$FDR < 0.10

  write.csv(
    res,
    file.path(
      TABDIR,
      paste0(
        label,
        "_Tx_vs_Sham_edgeR_v6.7.6.csv"
      )
    ),
    row.names=FALSE
  )

  # -----------------------------------------------
  # ranked reproducible candidates
  # -----------------------------------------------

  candidate <- res[
    abs(res$logFC) >= 0.5 &
    res$strict_direction_consistent !=
      "Not_strictly_consistent",
    ,
    drop=FALSE
  ]

  candidate <- candidate[
    order(
      candidate$FDR,
      -abs(candidate$logFC)
    ),
    ,
    drop=FALSE
  ]

  write.csv(
    candidate,
    file.path(
      TABDIR,
      paste0(
        label,
        "_replicate_consistent_effect_genes_v6.7.6.csv"
      )
    ),
    row.names=FALSE
  )

  # -----------------------------------------------
  # pseudobulk PCA
  # -----------------------------------------------

  pca <- prcomp(
    t(logcpm),
    scale.=FALSE
  )

  pca_df <- data.frame(
    sample=rownames(pca$x),
    PC1=pca$x[,1],
    PC2=pca$x[,2],
    condition=group
  )

  var_exp <- (
    pca$sdev^2 /
    sum(pca$sdev^2)
  ) * 100

  p <- ggplot(
    pca_df,
    aes(
      PC1,
      PC2,
      label=sample,
      shape=condition
    )
  ) +
    geom_point(size=4) +
    geom_text(
      nudge_y=0.15,
      check_overlap=TRUE
    ) +
    theme_classic(base_size=12) +
    labs(
      title=paste0(
        label,
        " pseudobulk PCA"
      ),
      x=paste0(
        "PC1 (",
        round(var_exp[1],1),
        "%)"
      ),
      y=paste0(
        "PC2 (",
        round(var_exp[2],1),
        "%)"
      )
    )

  ggsave(
    file.path(
      FIGDIR,
      paste0(
        label,
        "_pseudobulk_PCA_v6.7.6.pdf"
      )
    ),
    p,
    width=7,
    height=6
  )

  # -----------------------------------------------
  # sample correlation
  # -----------------------------------------------

  cor_mat <- cor(
    logcpm,
    method="pearson"
  )

  write.csv(
    cor_mat,
    file.path(
      TABDIR,
      paste0(
        label,
        "_sample_correlation_v6.7.6.csv"
      )
    )
  )

  cat(
    "Genes tested:",
    nrow(res),
    "\n"
  )

  cat(
    "FDR < 0.10:",
    sum(res$FDR < 0.10),
    "\n"
  )

  cat(
    "|logFC| >= 0.5 and strict replicate consistency:",
    nrow(candidate),
    "\n"
  )

  invisible(res)
}

# -------------------------------------------------------
# Run all analyses
# -------------------------------------------------------

for (nm in names(analyses)) {
  run_pb(
    nm,
    analyses[[nm]]
  )
}

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.6.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.6 COMPLETE\n")
cat("Primary contrast: Tx vs Sham\n")
cat("Statistical unit: biological sample\n")
cat("n = 2 Sham + 2 Tx\n")
cat("Interpret FDR cautiously; emphasize effect size and replicate concordance\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
