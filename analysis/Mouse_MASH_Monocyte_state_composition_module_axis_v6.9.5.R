suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.9.5"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.4/objects/",
  "Mouse_MASH_Monocyte_annotation_frozen_v6.9.4.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte state composition + module axis\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

STATE_COL <- "Monocyte_state_frozen_v6.9.4"
CLASS_COL <- "Monocyte_analysis_class_v6.9.4"

required_meta <- c(
  STATE_COL,
  CLASS_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse=", ")
  )
}

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

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

observed_samples <- unique(sample_vec)

sample_order <- sample_order[
  sample_order %in% observed_samples
]

extra_samples <- setdiff(
  sort(observed_samples),
  sample_order
)

sample_order <- c(
  sample_order,
  extra_samples
)

# =========================================================
# Condition
# =========================================================

condition <- ifelse(
  grepl("^STD", sample_vec, ignore.case=TRUE),
  "STD",
  ifelse(
    grepl("CDHFD|CDAHFD", sample_vec, ignore.case=TRUE),
    "CDHFD",
    ifelse(
      grepl("^Sham", sample_vec, ignore.case=TRUE),
      "Sham",
      ifelse(
        grepl("^Tx", sample_vec, ignore.case=TRUE),
        "Tx",
        "Other"
      )
    )
  )
)

obj$condition_v6.9.5 <- factor(
  condition,
  levels=c(
    "STD",
    "CDHFD",
    "Sham",
    "Tx",
    "Other"
  )
)

# =========================================================
# Analysis frameworks
# =========================================================

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

framework_summary <- data.frame(
  framework=names(framework_cells),
  n_cells=vapply(
    framework_cells,
    length,
    integer(1)
  ),
  stringsAsFactors=FALSE
)

write.csv(
  framework_summary,
  file.path(
    TABDIR,
    "Monocyte_analysis_framework_cell_counts_v6.9.5.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Frozen state levels
# =========================================================

state_levels <- levels(
  obj@meta.data[[STATE_COL]]
)

if (is.null(state_levels)) {
  state_levels <- sort(
    unique(
      as.character(
        obj@meta.data[[STATE_COL]]
      )
    )
  )
}

# =========================================================
# State composition by sample
# =========================================================

composition_rows <- list()
counter <- 1

for (fw in names(framework_cells)) {

  fw_cells <- framework_cells[[fw]]

  md <- obj@meta.data[
    fw_cells,
    ,
    drop=FALSE
  ]

  for (s in sample_order) {

    sample_cells <- rownames(md)[
      as.character(
        md[[sample_col]]
      ) == s
    ]

    n_sample <- length(sample_cells)

    for (state in state_levels) {

      n_state <- sum(
        as.character(
          md[
            sample_cells,
            STATE_COL
          ]
        ) == state
      )

      composition_rows[[counter]] <- data.frame(
        framework=fw,
        sample=s,
        condition=ifelse(
          grepl("^STD", s),
          "STD",
          ifelse(
            grepl("CDHFD|CDAHFD", s),
            "CDHFD",
            ifelse(
              grepl("^Sham", s),
              "Sham",
              ifelse(
                grepl("^Tx", s),
                "Tx",
                "Other"
              )
            )
          )
        ),
        state=state,
        n_state=n_state,
        n_sample=n_sample,
        fraction=if (
          n_sample > 0
        ) {
          n_state / n_sample
        } else {
          NA_real_
        },
        stringsAsFactors=FALSE
      )

      counter <- counter + 1
    }
  }
}

composition <- do.call(
  rbind,
  composition_rows
)

write.csv(
  composition,
  file.path(
    TABDIR,
    "Monocyte_state_composition_by_sample_v6.9.5.csv"
  ),
  row.names=FALSE
)

# =========================================================
# State axis summary
# =========================================================

state_axis_rows <- list()
counter <- 1

for (fw in names(framework_cells)) {

  for (state in state_levels) {

    x <- composition[
      composition$framework == fw &
      composition$state == state,
      ,
      drop=FALSE
    ]

    get_fraction <- function(sample_name) {

      y <- x$fraction[
        x$sample == sample_name
      ]

      if (length(y) == 0) {
        return(NA_real_)
      }

      y[[1]]
    }

    std <- get_fraction("STD_rep1")
    cdhfd <- get_fraction("CDHFD_rep1")
    sham1 <- get_fraction("Sham1")
    sham20 <- get_fraction("Sham20")
    tx17 <- get_fraction("Tx17")
    tx5 <- get_fraction("Tx5")

    sham_vals <- c(
      sham1,
      sham20
    )

    tx_vals <- c(
      tx17,
      tx5
    )

    sham_mean <- mean(
      sham_vals,
      na.rm=TRUE
    )

    tx_mean <- mean(
      tx_vals,
      na.rm=TRUE
    )

    if (
      all(is.finite(sham_vals)) &&
      all(is.finite(tx_vals))
    ) {

      if (
        min(tx_vals) >
        max(sham_vals)
      ) {

        replicate_pattern <-
          "Tx_all_higher"

      } else if (
        max(tx_vals) <
        min(sham_vals)
      ) {

        replicate_pattern <-
          "Tx_all_lower"

      } else {

        replicate_pattern <-
          "ranges_overlap"
      }

    } else {

      replicate_pattern <-
        "insufficient"
    }

    state_axis_rows[[counter]] <- data.frame(
      framework=fw,
      state=state,

      STD_fraction=std,
      CDHFD_fraction=cdhfd,
      Disease_delta_CDHFD_minus_STD=
        cdhfd - std,

      Sham1_fraction=sham1,
      Sham20_fraction=sham20,
      Sham_mean_fraction=sham_mean,

      Tx17_fraction=tx17,
      Tx5_fraction=tx5,
      Tx_mean_fraction=tx_mean,

      Tx_minus_Sham=
        tx_mean - sham_mean,

      replicate_pattern=
        replicate_pattern,

      stringsAsFactors=FALSE
    )

    counter <- counter + 1
  }
}

state_axis <- do.call(
  rbind,
  state_axis_rows
)

write.csv(
  state_axis,
  file.path(
    TABDIR,
    "Monocyte_state_axis_summary_v6.9.5.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Functional modules
#
# Prespecified before treatment comparison.
# Ly6c2 is unavailable in this dataset and is not used.
# =========================================================

module_sets <- list(

  Classical_recruitment=c(
    "Ccr2",
    "Sell",
    "Vcan",
    "S100a8",
    "S100a9",
    "Mmp8",
    "Chil3",
    "F13a1"
  ),

  Inflammatory_activation=c(
    "Il1b",
    "Tnf",
    "Nfkbia",
    "Thbs1",
    "Ccl2",
    "Ccl3",
    "Ccl4",
    "Cxcl2",
    "Il1rn"
  ),

  IFN_IFNg_response=c(
    "Cxcl9",
    "Cxcl10",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Isg15",
    "Rsad2",
    "Cmpk2",
    "Stat1",
    "Irf7"
  ),

  Macrophage_transition_remodeling=c(
    "Ms4a7",
    "Mmp12",
    "Dab2",
    "C1qa",
    "C1qb",
    "C1qc",
    "Mertk",
    "Lpl",
    "Gpnmb"
  ),

  Lipid_scavenging=c(
    "Cd36",
    "Olr1",
    "Lpl",
    "Gpnmb",
    "Apoe",
    "Lgals3",
    "Msr1"
  )
)

all_module_genes <- unique(
  unlist(
    module_sets,
    use.names=FALSE
  )
)

present_module_genes <- intersect(
  all_module_genes,
  rownames(obj)
)

missing_module_genes <- setdiff(
  all_module_genes,
  rownames(obj)
)

module_gene_presence <- data.frame(
  gene=all_module_genes,
  present=
    all_module_genes %in%
    rownames(obj),
  stringsAsFactors=FALSE
)

write.csv(
  module_gene_presence,
  file.path(
    TABDIR,
    "Monocyte_module_gene_presence_v6.9.5.csv"
  ),
  row.names=FALSE
)

cat(
  "Module genes present:",
  length(present_module_genes),
  "/",
  length(all_module_genes),
  "\n"
)

if (length(missing_module_genes) > 0) {

  cat(
    "Missing module genes:",
    paste(
      missing_module_genes,
      collapse=", "
    ),
    "\n"
  )
}

# =========================================================
# Normalize RNA
# =========================================================

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

data_mat <- GetAssayData(
  obj,
  assay="RNA",
  layer="data"
)

# =========================================================
# Gene-wise z-score module scoring
#
# Equal weight per gene.
# No control-gene sampling.
# =========================================================

score_module <- function(
  mat,
  genes
) {

  genes <- intersect(
    genes,
    rownames(mat)
  )

  if (length(genes) == 0) {

    return(
      rep(
        NA_real_,
        ncol(mat)
      )
    )
  }

  x <- as.matrix(
    mat[
      genes,
      ,
      drop=FALSE
    ]
  )

  gene_mean <- rowMeans(
    x
  )

  gene_sd <- apply(
    x,
    1,
    sd
  )

  valid <- is.finite(gene_sd) &
    gene_sd > 0

  x <- x[
    valid,
    ,
    drop=FALSE
  ]

  gene_mean <- gene_mean[
    valid
  ]

  gene_sd <- gene_sd[
    valid
  ]

  if (nrow(x) == 0) {

    return(
      rep(
        NA_real_,
        ncol(mat)
      )
    )
  }

  z <- sweep(
    x,
    1,
    gene_mean,
    "-"
  )

  z <- sweep(
    z,
    1,
    gene_sd,
    "/"
  )

  colMeans(
    z,
    na.rm=TRUE
  )
}

for (module_name in names(module_sets)) {

  score <- score_module(
    data_mat,
    module_sets[[module_name]]
  )

  obj[[paste0(
    "MonocyteModule_",
    module_name,
    "_v6.9.5"
  )]] <- score
}

# =========================================================
# Module summary by sample/framework
# =========================================================

module_rows <- list()
counter <- 1

for (fw in names(framework_cells)) {

  fw_cells <- framework_cells[[fw]]

  for (s in sample_order) {

    cells <- intersect(
      fw_cells,
      colnames(obj)[
        sample_vec == s
      ]
    )

    for (module_name in names(module_sets)) {

      score_col <- paste0(
        "MonocyteModule_",
        module_name,
        "_v6.9.5"
      )

      vals <- obj@meta.data[
        cells,
        score_col
      ]

      module_rows[[counter]] <- data.frame(
        framework=fw,
        sample=s,
        condition=ifelse(
          grepl("^STD", s),
          "STD",
          ifelse(
            grepl("CDHFD|CDAHFD", s),
            "CDHFD",
            ifelse(
              grepl("^Sham", s),
              "Sham",
              ifelse(
                grepl("^Tx", s),
                "Tx",
                "Other"
              )
            )
          )
        ),
        module=module_name,
        n_cells=length(cells),
        mean_score=mean(
          vals,
          na.rm=TRUE
        ),
        median_score=median(
          vals,
          na.rm=TRUE
        ),
        stringsAsFactors=FALSE
      )

      counter <- counter + 1
    }
  }
}

module_by_sample <- do.call(
  rbind,
  module_rows
)

write.csv(
  module_by_sample,
  file.path(
    TABDIR,
    "Monocyte_module_scores_by_sample_v6.9.5.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Module axis summary
# =========================================================

module_axis_rows <- list()
counter <- 1

for (fw in names(framework_cells)) {

  for (module_name in names(module_sets)) {

    x <- module_by_sample[
      module_by_sample$framework == fw &
      module_by_sample$module == module_name,
      ,
      drop=FALSE
    ]

    get_score <- function(sample_name) {

      y <- x$mean_score[
        x$sample == sample_name
      ]

      if (length(y) == 0) {
        return(NA_real_)
      }

      y[[1]]
    }

    std <- get_score("STD_rep1")
    cdhfd <- get_score("CDHFD_rep1")
    sham1 <- get_score("Sham1")
    sham20 <- get_score("Sham20")
    tx17 <- get_score("Tx17")
    tx5 <- get_score("Tx5")

    sham_vals <- c(
      sham1,
      sham20
    )

    tx_vals <- c(
      tx17,
      tx5
    )

    sham_mean <- mean(
      sham_vals,
      na.rm=TRUE
    )

    tx_mean <- mean(
      tx_vals,
      na.rm=TRUE
    )

    if (
      all(is.finite(sham_vals)) &&
      all(is.finite(tx_vals))
    ) {

      if (
        min(tx_vals) >
        max(sham_vals)
      ) {

        replicate_pattern <-
          "Tx_all_higher"

      } else if (
        max(tx_vals) <
        min(sham_vals)
      ) {

        replicate_pattern <-
          "Tx_all_lower"

      } else {

        replicate_pattern <-
          "ranges_overlap"
      }

    } else {

      replicate_pattern <-
        "insufficient"
    }

    module_axis_rows[[counter]] <- data.frame(
      framework=fw,
      module=module_name,

      STD_mean=std,
      CDHFD_mean=cdhfd,
      Disease_delta_CDHFD_minus_STD=
        cdhfd - std,

      Sham1_mean=sham1,
      Sham20_mean=sham20,
      Sham_mean=sham_mean,

      Tx17_mean=tx17,
      Tx5_mean=tx5,
      Tx_mean=tx_mean,

      Tx_minus_Sham=
        tx_mean - sham_mean,

      replicate_pattern=
        replicate_pattern,

      stringsAsFactors=FALSE
    )

    counter <- counter + 1
  }
}

module_axis <- do.call(
  rbind,
  module_axis_rows
)

write.csv(
  module_axis,
  file.path(
    TABDIR,
    "Monocyte_module_axis_summary_v6.9.5.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Figures
# =========================================================

plot_composition <- composition[
  composition$framework == "ALL",
  ,
  drop=FALSE
]

plot_composition$sample <- factor(
  plot_composition$sample,
  levels=sample_order
)

p_comp <- ggplot(
  plot_composition,
  aes(
    x=sample,
    y=fraction,
    fill=state
  )
) +
  geom_col(
    width=0.82
  ) +
  theme_classic(
    base_size=9
  ) +
  theme(
    axis.text.x=
      element_text(
        angle=45,
        hjust=1
      ),
    legend.position="right"
  ) +
  labs(
    title="Mouse MASH Monocyte frozen-state composition",
    x=NULL,
    y="Fraction of Monocytes",
    fill="Frozen state"
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_state_composition_ALL_v6.9.5.pdf"
  ),
  p_comp,
  width=13,
  height=7
)

module_plot <- module_by_sample[
  module_by_sample$framework ==
    "PRIMARY_CORE",
  ,
  drop=FALSE
]

module_plot$sample <- factor(
  module_plot$sample,
  levels=sample_order
)

p_module <- ggplot(
  module_plot,
  aes(
    x=sample,
    y=mean_score,
    group=1
  )
) +
  geom_line(
    linewidth=0.6
  ) +
  geom_point(
    size=2.5
  ) +
  facet_wrap(
    ~module,
    scales="free_y",
    ncol=2
  ) +
  theme_classic(
    base_size=10
  ) +
  theme(
    axis.text.x=
      element_text(
        angle=45,
        hjust=1
      )
  ) +
  labs(
    title="Mouse MASH Monocyte functional module axis - PRIMARY_CORE",
    x=NULL,
    y="Mean gene-wise z-score"
  )

ggsave(
  file.path(
    FIGDIR,
    "Monocyte_functional_module_axis_PRIMARY_CORE_v6.9.5.pdf"
  ),
  p_module,
  width=11,
  height=10
)

# =========================================================
# Save scored object
# =========================================================

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds"
  )
)

# =========================================================
# Terminal output
# =========================================================

cat("\n=== ANALYSIS FRAMEWORKS ===\n")
print(
  framework_summary,
  row.names=FALSE
)

cat("\n=== STATE AXIS ALL ===\n")
print(
  state_axis[
    state_axis$framework == "ALL",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== STATE AXIS NO_QCWATCH ===\n")
print(
  state_axis[
    state_axis$framework == "NO_QCWATCH",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== STATE AXIS PRIMARY_CORE ===\n")
print(
  state_axis[
    state_axis$framework == "PRIMARY_CORE",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== MODULE AXIS ALL ===\n")
print(
  module_axis[
    module_axis$framework == "ALL",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== MODULE AXIS NO_QCWATCH ===\n")
print(
  module_axis[
    module_axis$framework == "NO_QCWATCH",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n=== MODULE AXIS PRIMARY_CORE ===\n")
print(
  module_axis[
    module_axis$framework == "PRIMARY_CORE",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.5 COMPLETE\n")
cat("Monocyte state composition + module axis complete\n")
cat("Cells:", ncol(obj), "\n")
cat("STD/CDHFD: descriptive only\n")
cat("Sham/Tx: n=2/group; replicate concordance reported\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.5.txt"
  )
)
