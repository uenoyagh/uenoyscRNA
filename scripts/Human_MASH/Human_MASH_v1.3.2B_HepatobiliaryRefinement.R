# ============================================================
# Human MASH v1.3.2B
# Hepatobiliary Refinement
#
# PURPOSE
#   Refine the hepatocyte / immature-hepatocyte / ductal
#   compartment of the canonical Human MASH LT1-LT8 object.
#
# SCIENTIFIC QUESTIONS
#   1. Is "Hepatoblast_or_immature_hepatocyte" a true fetal/
#      hepatoblast-like population, or mainly an adult
#      injury/transitional/ductular program?
#   2. Are bona fide cholangiocytes hidden within the current
#      immature-hepatocyte compartment?
#   3. Can hepatocyte zonation be separated from injury/
#      progenitor/ductular states?
#
# DESIGN
#   - Canonical whole-liver v1.3 object is READ ONLY.
#   - Core compartment:
#       Hepatocyte
#       Hepatoblast_or_immature_hepatocyte
#       Cholangiocyte_or_ductal
#   - Strict ductal-like cells outside the core compartment are
#     audited but NOT automatically rescued.
#   - LT1-LT8 RPCA integration is performed only within the
#     hepatobiliary core compartment.
#   - Multi-resolution clustering: 0.2-1.2.
#   - res0.6 is used as the provisional review resolution.
#   - Cluster marker discovery + curated state/signature audit.
#   - NO final hepatobiliary cell-type labels are frozen here.
#   - Manual cluster annotation is deferred to v1.3.2B.1.
#
# IMPORTANT
#   - No canonical whole-liver RDS is modified.
#   - No cluster number is interpreted across species.
#   - Integrated assay is used for clustering/UMAP.
#   - RNA assay is used for biological marker interpretation.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(scales)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 12 * 1024^3)

seed <- 20260908
set.seed(seed)

# ------------------------------------------------------------
# 0. Paths / configuration
# ------------------------------------------------------------

input_rds <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/",
  "Human_MASH_LT1to8_Canonical_v1.3.rds"
)

outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B_hepatobiliary_refinement"
)

object_dir <- file.path(
  outdir,
  "objects"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

celltype_col <- "celltype_ref12_v13"
sample_col   <- "canonical_sample_v13"
nas_col      <- "NAS_score_v13"
nasgrp_col   <- "NAS_group_v13"

core_labels <- c(
  "Hepatocyte",
  "Hepatoblast_or_immature_hepatocyte",
  "Cholangiocyte_or_ductal"
)

expected_canonical_cells <- 81348L

expected_label_counts <- c(
  Hepatocyte = 37956L,
  Hepatoblast_or_immature_hepatocyte = 9467L,
  Cholangiocyte_or_ductal = 2L
)

expected_core_cells <- sum(
  expected_label_counts
)

annotation_resolution <- 0.6

resolution_grid <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0,
  1.2
)

nfeatures_integration <- 3000
npcs <- 30

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

get_data_layer <- function(
    object,
    assay = "RNA"
) {

  tryCatch(
    GetAssayData(
      object,
      assay = assay,
      layer = "data"
    ),
    error = function(e) {

      GetAssayData(
        object,
        assay = assay,
        slot = "data"
      )
    }
  )
}

sample_order_fun <- function(x) {

  u <- unique(
    as.character(x)
  )

  if (
    length(u) > 0 &&
    all(
      grepl(
        "^LT[0-9]+",
        u
      )
    )
  ) {

    n <- suppressWarnings(
      as.integer(
        sub(
          "^LT([0-9]+).*",
          "\\1",
          u
        )
      )
    )

    if (!anyNA(n)) {

      return(
        u[
          order(n)
        ]
      )
    }
  }

  sort(u)
}

safe_present <- function(
    genes,
    universe
) {

  intersect(
    unique(genes),
    universe
  )
}

dedupe_feature_list <- function(x) {

  seen <- character(0)

  out <- lapply(
    x,
    function(g) {

      g <- unique(g)

      keep <- g[
        !g %in% seen
      ]

      seen <<- c(
        seen,
        keep
      )

      keep
    }
  )

  out[
    lengths(out) > 0
  ]
}

resolution_colname <- function(res) {

  paste0(
    "hepato_res",
    gsub(
      "\\.",
      "_",
      as.character(res)
    ),
    "_v132B"
  )
}

# ------------------------------------------------------------
# 2. Load canonical v1.3 object
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.2B Hepatobiliary Refinement\n")
cat("============================================================\n\n")

stopifnot(
  file.exists(input_rds)
)

obj <- readRDS(
  input_rds
)

DefaultAssay(obj) <- "RNA"

required_meta <- c(
  celltype_col,
  sample_col,
  nas_col,
  nasgrp_col
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    obj@meta.data
  )
)

if (
  length(
    missing_meta
  ) > 0
) {

  stop(
    "Missing required metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  ncol(obj) != expected_canonical_cells
) {

  stop(
    "Canonical-cell-count mismatch. Expected ",
    expected_canonical_cells,
    " but found ",
    ncol(obj),
    ". Do NOT continue until object provenance is checked."
  )
}

cat(
  "Canonical cells :",
  ncol(obj),
  "\n"
)

cat(
  "RNA genes       :",
  nrow(obj[["RNA"]]),
  "\n"
)

cat(
  "Samples         :",
  paste(
    sort(
      unique(
        obj[[sample_col]][, 1]
      )
    ),
    collapse = " | "
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 3. Verify current target-label counts exactly
# ------------------------------------------------------------

current_celltype <- as.character(
  obj[[celltype_col]][, 1]
)

current_sample <- as.character(
  obj[[sample_col]][, 1]
)

observed_label_counts <- table(
  factor(
    current_celltype,
    levels = names(
      expected_label_counts
    )
  )
)

observed_label_counts <- as.integer(
  observed_label_counts
)

names(
  observed_label_counts
) <- names(
  expected_label_counts
)

if (
  !identical(
    as.integer(
      observed_label_counts
    ),
    as.integer(
      expected_label_counts
    )
  )
) {

  stop(
    paste0(
      "Current target-label counts do not match the validated v1.3.1 object.\n",
      "Expected: ",
      paste(
        names(expected_label_counts),
        expected_label_counts,
        sep = "=",
        collapse = ", "
      ),
      "\nObserved: ",
      paste(
        names(observed_label_counts),
        observed_label_counts,
        sep = "=",
        collapse = ", "
      ),
      "\nDo NOT apply v1.3.2B until provenance is resolved."
    )
  )
}

target_label_counts <- data.frame(
  celltype = names(
    observed_label_counts
  ),
  n_cells = as.integer(
    observed_label_counts
  ),
  stringsAsFactors = FALSE
)

target_label_counts$percent_of_core <- 100 *
  target_label_counts$n_cells /
  sum(
    target_label_counts$n_cells
  )

write.csv(
  target_label_counts,
  file.path(
    outdir,
    "01_current_target_label_counts.csv"
  ),
  row.names = FALSE
)

cat("=== CURRENT HEPATOBILIARY CORE LABELS ===\n")
print(
  target_label_counts,
  row.names = FALSE
)
cat("\n")

# ------------------------------------------------------------
# 4. Strict ductal-like rescue AUDIT outside the core
#
# IMPORTANT:
#   This is audit-only. These cells are NOT included in the
#   RPCA object at v1.3.2B.
#
# Rule:
#   - outside the three current hepatobiliary labels
#   - PTPRC == 0 when PTPRC is available
#   - >=2 genes detected among:
#       KRT7/KRT19/SOX9/MUC1/CFTR/EPCAM/PROM1
# ------------------------------------------------------------

rna_data_whole <- get_data_layer(
  obj,
  assay = "RNA"
)

gene_universe <- rownames(
  rna_data_whole
)

strict_ductal_genes <- c(
  "KRT7",
  "KRT19",
  "SOX9",
  "MUC1",
  "CFTR",
  "EPCAM",
  "PROM1"
)

strict_ductal_present <- safe_present(
  strict_ductal_genes,
  gene_universe
)

core_mask <- current_celltype %in%
  core_labels

if (
  length(
    strict_ductal_present
  ) >= 2
) {

  ductal_hits <- Matrix::colSums(
    rna_data_whole[
      strict_ductal_present,
      ,
      drop = FALSE
    ] > 0
  )

  if (
    "PTPRC" %in%
      gene_universe
  ) {

    ptprc_negative <- as.numeric(
      rna_data_whole[
        "PTPRC",
        ,
        drop = TRUE
      ]
    ) == 0

  } else {

    ptprc_negative <- rep(
      TRUE,
      ncol(obj)
    )
  }

  strict_rescue_candidate <- (
    !core_mask &
    ptprc_negative &
    ductal_hits >= 2
  )

} else {

  warning(
    "Fewer than 2 strict ductal genes available; rescue audit skipped."
  )

  strict_rescue_candidate <- rep(
    FALSE,
    ncol(obj)
  )
}

rescue_candidate_table <- data.frame(
  cell = colnames(obj)[
    strict_rescue_candidate
  ],
  current_celltype =
    current_celltype[
      strict_rescue_candidate
    ],
  sample =
    current_sample[
      strict_rescue_candidate
    ],
  stringsAsFactors = FALSE
)

write.csv(
  rescue_candidate_table,
  file.path(
    outdir,
    "02_strict_ductal_rescue_candidate_cells_OUTSIDE_CORE.csv"
  ),
  row.names = FALSE
)

rescue_candidate_by_label <- as.data.frame(
  table(
    current_celltype =
      current_celltype[
        strict_rescue_candidate
      ]
  ),
  stringsAsFactors = FALSE
)

names(
  rescue_candidate_by_label
)[2] <- "n_candidate_cells"

write.csv(
  rescue_candidate_by_label,
  file.path(
    outdir,
    "03_strict_ductal_rescue_candidates_by_current_label.csv"
  ),
  row.names = FALSE
)

cat(
  "Strict ductal-like rescue candidates outside core:",
  sum(
    strict_rescue_candidate
  ),
  "\n"
)

cat(
  "NOTE: audit only; none are added automatically.\n\n"
)

# ------------------------------------------------------------
# 5. Create hepatobiliary core subset
# ------------------------------------------------------------

core_cells <- colnames(obj)[
  core_mask
]

if (
  length(
    core_cells
  ) != expected_core_cells
) {

  stop(
    "Core-cell-count mismatch. Expected ",
    expected_core_cells,
    " but found ",
    length(core_cells),
    "."
  )
}

hepato <- subset(
  obj,
  cells = core_cells
)

hepato$current_hepatobiliary_label_v132B <-
  as.character(
    hepato[[celltype_col]][, 1]
  )

DefaultAssay(
  hepato
) <- "RNA"

rm(
  rna_data_whole
)

gc()

cat(
  "Hepatobiliary core cells:",
  ncol(hepato),
  "\n\n"
)

# ------------------------------------------------------------
# 6. Sample-level input QC
# ------------------------------------------------------------

sample_order <- sample_order_fun(
  hepato[[sample_col]][, 1]
)

sample_counts <- as.data.frame(
  table(
    sample =
      hepato[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

names(
  sample_counts
)[2] <- "n_cells"

sample_counts$sample <- factor(
  sample_counts$sample,
  levels = sample_order
)

sample_counts <- sample_counts[
  order(
    sample_counts$sample
  ),
]

sample_counts$sample <- as.character(
  sample_counts$sample
)

write.csv(
  sample_counts,
  file.path(
    outdir,
    "04_hepatobiliary_cells_by_sample.csv"
  ),
  row.names = FALSE
)

cat("=== HEPATOBILIARY CELLS BY SAMPLE ===\n")
print(
  sample_counts,
  row.names = FALSE
)
cat("\n")

if (
  any(
    sample_counts$n_cells < 100
  )
) {

  stop(
    "At least one LT sample has <100 hepatobiliary cells; ",
    "RPCA integration aborted."
  )
}

# ------------------------------------------------------------
# 7. Split by sample and prepare RPCA integration
# ------------------------------------------------------------

cat("=== RPCA PREPARATION ===\n")

hepato_list <- SplitObject(
  hepato,
  split.by = sample_col
)

hepato_list <- hepato_list[
  sample_order
]

hepato_list <- lapply(
  hepato_list,
  function(x) {

    DefaultAssay(x) <- "RNA"

    x <- NormalizeData(
      x,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )

    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = nfeatures_integration,
      verbose = FALSE
    )

    x
  }
)

integration_features <- SelectIntegrationFeatures(
  object.list = hepato_list,
  nfeatures = nfeatures_integration
)

hepato_list <- lapply(
  hepato_list,
  function(x) {

    x <- ScaleData(
      x,
      features =
        integration_features,
      verbose = FALSE
    )

    x <- RunPCA(
      x,
      features =
        integration_features,
      npcs = npcs,
      verbose = FALSE
    )

    x
  }
)

# ------------------------------------------------------------
# 8. RPCA anchors and integration
# ------------------------------------------------------------

cat("=== RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list = hepato_list,
  anchor.features =
    integration_features,
  reduction = "rpca",
  dims = 1:npcs,
  k.anchor = 5,
  verbose = FALSE
)

cat("=== INTEGRATE DATA ===\n")

hepato_int <- IntegrateData(
  anchorset = anchors,
  dims = 1:npcs,
  new.assay.name = "integrated",
  verbose = FALSE
)

DefaultAssay(
  hepato_int
) <- "integrated"

hepato_int <- ScaleData(
  hepato_int,
  verbose = FALSE
)

hepato_int <- RunPCA(
  hepato_int,
  npcs = npcs,
  verbose = FALSE
)

hepato_int <- RunUMAP(
  hepato_int,
  reduction = "pca",
  dims = 1:npcs,
  seed.use = seed,
  reduction.name =
    "umapHepatoRPCA_v132B",
  reduction.key =
    "UMAPHEP132B_",
  verbose = FALSE
)

hepato_int <- FindNeighbors(
  hepato_int,
  reduction = "pca",
  dims = 1:npcs,
  verbose = FALSE
)

# ------------------------------------------------------------
# 9. Multi-resolution clustering
# ------------------------------------------------------------

cat("=== MULTI-RESOLUTION CLUSTERING ===\n")

for (
  res in resolution_grid
) {

  set.seed(seed)

  hepato_int <- FindClusters(
    hepato_int,
    resolution = res,
    algorithm = 1,
    random.seed = seed,
    verbose = FALSE
  )

  metadata_name <- resolution_colname(
    res
  )

  hepato_int[[metadata_name]] <- as.character(
    Idents(
      hepato_int
    )
  )
}

annotation_cluster_col <- resolution_colname(
  annotation_resolution
)

hepato_int$hepatobiliary_cluster_v132B <- as.character(
  hepato_int[[annotation_cluster_col]][, 1]
)

Idents(
  hepato_int
) <- "hepatobiliary_cluster_v132B"

# ------------------------------------------------------------
# 10. Cluster counts and resolution audit
# ------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(
    cluster =
      hepato_int$hepatobiliary_cluster_v132B
  ),
  stringsAsFactors = FALSE
)

names(
  cluster_counts
)[2] <- "n_cells"

write.csv(
  cluster_counts,
  file.path(
    outdir,
    "05_res0.6_cluster_counts.csv"
  ),
  row.names = FALSE
)

resolution_count_list <- lapply(
  resolution_grid,
  function(res) {

    col <- resolution_colname(
      res
    )

    x <- as.data.frame(
      table(
        cluster =
          hepato_int[[col]][, 1]
      ),
      stringsAsFactors = FALSE
    )

    names(x)[2] <- "n_cells"

    x$resolution <- res

    x
  }
)

resolution_counts <- do.call(
  rbind,
  resolution_count_list
)

write.csv(
  resolution_counts,
  file.path(
    outdir,
    "06_cluster_counts_all_resolutions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Current annotation vs res0.6 crosswalk
# ------------------------------------------------------------

crosswalk <- as.data.frame(
  table(
    current_label =
      hepato_int$current_hepatobiliary_label_v132B,
    cluster =
      hepato_int$hepatobiliary_cluster_v132B
  ),
  stringsAsFactors = FALSE
)

write.csv(
  crosswalk,
  file.path(
    outdir,
    "07_current_annotation_vs_res0.6_cluster_crosswalk.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Sample contribution by cluster
# ------------------------------------------------------------

sample_vec <- as.character(
  hepato_int[[sample_col]][, 1]
)

cluster_vec <- as.character(
  hepato_int$hepatobiliary_cluster_v132B
)

cluster_sample <- as.data.frame(
  table(
    cluster =
      cluster_vec,
    sample =
      sample_vec
  ),
  stringsAsFactors = FALSE
)

cluster_sample$cluster_total <- ave(
  cluster_sample$Freq,
  cluster_sample$cluster,
  FUN = sum
)

cluster_sample$sample_total <- ave(
  cluster_sample$Freq,
  cluster_sample$sample,
  FUN = sum
)

cluster_sample$within_cluster_pct <- ifelse(
  cluster_sample$cluster_total > 0,
  100 *
    cluster_sample$Freq /
    cluster_sample$cluster_total,
  NA_real_
)

cluster_sample$within_sample_pct <- ifelse(
  cluster_sample$sample_total > 0,
  100 *
    cluster_sample$Freq /
    cluster_sample$sample_total,
  NA_real_
)

cluster_sample$sample <- factor(
  cluster_sample$sample,
  levels = sample_order
)

write.csv(
  cluster_sample,
  file.path(
    outdir,
    "08_cluster_by_sample_long.csv"
  ),
  row.names = FALSE
)

split_cluster <- split(
  seq_len(
    nrow(
      cluster_sample
    )
  ),
  cluster_sample$cluster
)

max_idx <- unlist(
  lapply(
    split_cluster,
    function(idx) {

      idx[
        which.max(
          cluster_sample$within_cluster_pct[idx]
        )
      ]
    }
  )
)

cluster_dominance <- cluster_sample[
  max_idx,
  c(
    "cluster",
    "sample",
    "within_cluster_pct"
  )
]

names(
  cluster_dominance
)[3] <- "max_sample_contribution_pct"

cluster_dominance$sample <- as.character(
  cluster_dominance$sample
)

cluster_dominance <- cluster_dominance[
  order(
    cluster_dominance$max_sample_contribution_pct,
    decreasing = TRUE
  ),
]

write.csv(
  cluster_dominance,
  file.path(
    outdir,
    "09_cluster_max_sample_contribution.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Curated biological signature panel
#
# NOTE:
#   These are diagnostic state/signature scores.
#   They are NOT used to freeze cell-type labels automatically.
# ------------------------------------------------------------

signatures <- list(

  Mature_hepatocyte = c(
    "ALB",
    "APOA1",
    "TTR",
    "ASGR1",
    "HNF4A",
    "CYP3A4",
    "CYP2E1",
    "CPS1"
  ),

  Periportal_hepatocyte = c(
    "CPS1",
    "ASS1",
    "ARG1",
    "HAL",
    "PCK1",
    "GLS2",
    "SDS"
  ),

  Pericentral_hepatocyte = c(
    "GLUL",
    "CYP2E1",
    "CYP1A2",
    "OAT",
    "AXIN2"
  ),

  Fetal_hepatoblast = c(
    "AFP",
    "DLK1",
    "IGF2",
    "H19",
    "PROM1",
    "EPCAM"
  ),

  Progenitor_ductular = c(
    "EPCAM",
    "PROM1",
    "KRT19",
    "KRT7",
    "SOX9",
    "KRT23",
    "TACSTD2"
  ),

  Cholangiocyte = c(
    "KRT7",
    "KRT19",
    "SOX9",
    "MUC1",
    "CFTR",
    "EPCAM",
    "KRT17"
  ),

  Injury_transitional_hepatocyte = c(
    "KRT8",
    "KRT18",
    "KRT23",
    "LGALS3",
    "CLDN1",
    "SPP1",
    "KLF6"
  ),

  Acute_phase_stress = c(
    "SAA1",
    "SAA2",
    "CRP",
    "HP",
    "ORM1",
    "FGA",
    "FGB",
    "FGG"
  ),

  Cycling_epithelial = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "RRM2"
  ),

  Immune_contamination = c(
    "PTPRC",
    "LST1",
    "TYROBP",
    "FCER1G"
  ),

  Endothelial_contamination = c(
    "PECAM1",
    "VWF",
    "KDR",
    "EMCN"
  ),

  Mesenchymal_contamination = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "LUM"
  )
)

DefaultAssay(
  hepato_int
) <- "RNA"

if (
  "JoinLayers" %in%
    getNamespaceExports(
      "SeuratObject"
    )
) {

  joined_rna <- try(
    SeuratObject::JoinLayers(
      hepato_int[["RNA"]]
    ),
    silent = TRUE
  )

  if (
    !inherits(
      joined_rna,
      "try-error"
    )
  ) {

    hepato_int[["RNA"]] <-
      joined_rna
  }
}

rna_data <- get_data_layer(
  hepato_int,
  assay = "RNA"
)

rna_genes <- rownames(
  rna_data
)

signature_availability <- do.call(
  rbind,
  lapply(
    names(
      signatures
    ),
    function(sig) {

      genes <- signatures[[sig]]

      data.frame(
        signature = sig,
        gene = genes,
        present =
          genes %in%
            rna_genes,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  signature_availability
) <- NULL

write.csv(
  signature_availability,
  file.path(
    outdir,
    "10_signature_gene_availability.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Cell-level signature scores
# ------------------------------------------------------------

score_columns <- character(0)

for (
  sig in names(
    signatures
  )
) {

  genes_present <- safe_present(
    signatures[[sig]],
    rna_genes
  )

  score_name <- paste0(
    "score_",
    sig,
    "_v132B"
  )

  if (
    length(
      genes_present
    ) == 0
  ) {

    hepato_int[[score_name]] <- 0

  } else {

    hepato_int[[score_name]] <-
      Matrix::colMeans(
        rna_data[
          genes_present,
          ,
          drop = FALSE
        ]
      )
  }

  score_columns <- c(
    score_columns,
    score_name
  )
}

# ------------------------------------------------------------
# 15. Cluster-level signature scores
# ------------------------------------------------------------

metadata_scoring <- hepato_int@meta.data[
  ,
  c(
    "hepatobiliary_cluster_v132B",
    score_columns
  ),
  drop = FALSE
]

cluster_score_raw <- aggregate(
  metadata_scoring[
    ,
    score_columns,
    drop = FALSE
  ],
  by = list(
    cluster =
      metadata_scoring$hepatobiliary_cluster_v132B
  ),
  FUN = mean
)

rownames(
  cluster_score_raw
) <- cluster_score_raw$cluster

score_matrix <- as.matrix(
  cluster_score_raw[
    ,
    score_columns,
    drop = FALSE
  ]
)

colnames(
  score_matrix
) <- sub(
  "^score_",
  "",
  colnames(score_matrix)
)

colnames(
  score_matrix
) <- sub(
  "_v132B$",
  "",
  colnames(score_matrix)
)

score_z <- scale(
  score_matrix,
  center = TRUE,
  scale = TRUE
)

score_z[
  !is.finite(
    score_z
  )
] <- 0

score_raw_out <- data.frame(
  cluster =
    rownames(
      score_matrix
    ),
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

write.csv(
  score_raw_out,
  file.path(
    outdir,
    "11_cluster_signature_raw_scores.csv"
  ),
  row.names = FALSE
)

score_z_out <- data.frame(
  cluster =
    rownames(
      score_z
    ),
  score_z,
  check.names = FALSE,
  row.names = NULL
)

write.csv(
  score_z_out,
  file.path(
    outdir,
    "12_cluster_signature_Zscores.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Top-two diagnostic signatures per cluster
# ------------------------------------------------------------

signature_top2 <- do.call(
  rbind,
  lapply(
    seq_len(
      nrow(
        score_z
      )
    ),
    function(i) {

      x <- score_z[i, ]

      ord <- order(
        x,
        decreasing = TRUE
      )

      data.frame(
        cluster =
          rownames(
            score_z
          )[i],
        top_signature =
          names(x)[
            ord[1]
          ],
        top_z =
          round(
            x[
              ord[1]
            ],
            4
          ),
        second_signature =
          names(x)[
            ord[2]
          ],
        second_z =
          round(
            x[
              ord[2]
            ],
            4
          ),
        margin =
          round(
            x[
              ord[1]
            ] -
              x[
                ord[2]
              ],
            4
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  signature_top2,
  file.path(
    outdir,
    "13_cluster_signature_top2_DIAGNOSTIC_ONLY.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. Cluster marker discovery
# ------------------------------------------------------------

DefaultAssay(
  hepato_int
) <- "RNA"

Idents(
  hepato_int
) <- "hepatobiliary_cluster_v132B"

cat("=== FIND ALL CLUSTER MARKERS ===\n")

cluster_markers <- FindAllMarkers(
  hepato_int,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

write.csv(
  cluster_markers,
  file.path(
    outdir,
    "14_hepatobiliary_cluster_markers_all.csv"
  ),
  row.names = FALSE
)

if (
  nrow(
    cluster_markers
  ) > 0
) {

  marker_split <- split(
    cluster_markers,
    cluster_markers$cluster
  )

  top_markers <- do.call(
    rbind,
    lapply(
      marker_split,
      function(x) {

        fc_col <- intersect(
          c(
            "avg_log2FC",
            "avg_logFC"
          ),
          colnames(x)
        )

        if (
          length(
            fc_col
          ) > 0
        ) {

          x <- x[
            order(
              x[[fc_col[1]]],
              decreasing = TRUE
            ),
          ]
        }

        head(
          x,
          20
        )
      }
    )
  )

  write.csv(
    top_markers,
    file.path(
      outdir,
      "15_hepatobiliary_cluster_markers_top20.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 18. PRE-PLOT CHECKPOINT
# ------------------------------------------------------------

preplot_checkpoint_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatobiliary_v1.3.2B_preplot_checkpoint.rds"
)

saveRDS(
  hepato_int,
  preplot_checkpoint_rds
)

writeLines(
  c(
    "Human MASH v1.3.2B pre-plot checkpoint",
    paste0(
      "cells=",
      ncol(
        hepato_int
      )
    ),
    paste0(
      "annotation_resolution=",
      annotation_resolution
    ),
    paste0(
      "clusters=",
      length(
        unique(
          hepato_int$hepatobiliary_cluster_v132B
        )
      )
    ),
    paste0(
      "checkpoint_rds=",
      preplot_checkpoint_rds
    ),
    "RPCA_complete=TRUE",
    "clustering_complete=TRUE",
    "signature_scoring_complete=TRUE",
    "FindAllMarkers_complete=TRUE"
  ),
  file.path(
    outdir,
    "15A_preplot_checkpoint.txt"
  )
)

cat(
  "PRE-PLOT CHECKPOINT SAVED:\n",
  preplot_checkpoint_rds,
  "\n\n"
)

# ------------------------------------------------------------
# 19. Palettes
# ------------------------------------------------------------

current_palette <- c(
  "Hepatocyte" =
    "#00A6A6",

  "Hepatoblast_or_immature_hepatocyte" =
    "#0072B2",

  "Cholangiocyte_or_ductal" =
    "#00C853"
)

cluster_ids <- sort(
  unique(
    hepato_int$hepatobiliary_cluster_v132B
  )
)

cluster_palette <- setNames(
  grDevices::hcl.colors(
    length(
      cluster_ids
    ),
    palette = "Dark 3"
  ),
  cluster_ids
)

# ------------------------------------------------------------
# 20. Fig01 - current annotation on refined UMAP
# ------------------------------------------------------------

p01 <- DimPlot(
  hepato_int,
  reduction =
    "umapHepatoRPCA_v132B",
  group.by =
    "current_hepatobiliary_label_v132B",
  pt.size = 0.10,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    current_palette
) +
  ggtitle(
    "v1.3.2B hepatobiliary subset - current annotation"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig01_CurrentAnnotation_HepatobiliaryRPCA_UMAP.pdf"
  ),
  p01,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 21. Fig02 - res0.6 clusters
# ------------------------------------------------------------

p02 <- DimPlot(
  hepato_int,
  reduction =
    "umapHepatoRPCA_v132B",
  group.by =
    "hepatobiliary_cluster_v132B",
  pt.size = 0.10,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    cluster_palette
) +
  ggtitle(
    paste0(
      "Hepatobiliary RPCA clusters - resolution ",
      annotation_resolution
    )
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig02_HepatobiliaryClusters_res0.6_UMAP.pdf"
  ),
  p02,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 22. Fig03 - clusters split by LT sample
# ------------------------------------------------------------

p03 <- DimPlot(
  hepato_int,
  reduction =
    "umapHepatoRPCA_v132B",
  group.by =
    "hepatobiliary_cluster_v132B",
  split.by = sample_col,
  ncol = 4,
  pt.size = 0.08,
  raster = TRUE,
  cols =
    cluster_palette
)

ggsave(
  file.path(
    outdir,
    "Fig03_HepatobiliaryClusters_res0.6_split_by_LT.pdf"
  ),
  p03,
  width = 19,
  height = 13
)

# ------------------------------------------------------------
# 23. Fig04 - current annotation vs cluster crosswalk
# ------------------------------------------------------------

p04 <- ggplot(
  crosswalk,
  aes(
    x = cluster,
    y = current_label,
    fill = log1p(Freq)
  )
) +
  geom_tile(
    color = "grey85",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        Freq > 0,
        Freq,
        ""
      )
    ),
    size = 2.4
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Current hepatobiliary annotation vs res0.6 clusters",
    x =
      "res0.6 cluster",
    y =
      "Current REF12 annotation",
    fill =
      "log(1+n)"
  ) +
  theme_bw(
    base_size = 10
  )

ggsave(
  file.path(
    outdir,
    "Fig04_CurrentAnnotation_vs_res0.6_Crosswalk.pdf"
  ),
  p04,
  width = 13,
  height = 5
)

# ------------------------------------------------------------
# 24. Fig05 - cluster-level signature Z-score heatmap
# ------------------------------------------------------------

z_long <- data.frame(
  cluster =
    rep(
      rownames(
        score_z
      ),
      times =
        ncol(
          score_z
        )
    ),
  signature =
    rep(
      colnames(
        score_z
      ),
      each =
        nrow(
          score_z
        )
    ),
  z =
    as.vector(
      score_z
    ),
  stringsAsFactors = FALSE
)

p05 <- ggplot(
  z_long,
  aes(
    x = signature,
    y = cluster,
    fill = z
  )
) +
  geom_tile(
    color = "grey85",
    linewidth = 0.25
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(
      -2.5,
      2.5
    ),
    oob =
      scales::squish
  ) +
  labs(
    title =
      "Hepatobiliary cluster diagnostic signature scores",
    subtitle =
      "Z scores are diagnostic only; not final annotation",
    x = NULL,
    y =
      "res0.6 cluster",
    fill =
      "Z score"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig05_ClusterSignature_Zscore_Heatmap.pdf"
  ),
  p05,
  width = 15,
  height = 9
)

# ------------------------------------------------------------
# 25. Fig06 - curated marker DotPlot by res0.6 cluster
# ------------------------------------------------------------

marker_panel <- list(

  Mature_hepatocyte = c(
    "ALB",
    "APOA1",
    "TTR",
    "ASGR1",
    "HNF4A",
    "CYP3A4"
  ),

  Hepatocyte_zonation = c(
    "CPS1",
    "ASS1",
    "HAL",
    "PCK1",
    "GLUL",
    "CYP2E1",
    "CYP1A2",
    "OAT"
  ),

  Fetal_hepatoblast = c(
    "AFP",
    "DLK1",
    "IGF2",
    "H19"
  ),

  Progenitor_ductular = c(
    "EPCAM",
    "PROM1",
    "KRT19",
    "KRT7",
    "SOX9",
    "KRT23",
    "TACSTD2"
  ),

  Cholangiocyte = c(
    "MUC1",
    "CFTR",
    "KRT17"
  ),

  Injury_stress = c(
    "KRT8",
    "KRT18",
    "LGALS3",
    "CLDN1",
    "SPP1",
    "SAA1",
    "SAA2"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF"
  ),

  Contamination_audit = c(
    "PTPRC",
    "LST1",
    "PECAM1",
    "VWF",
    "COL1A1",
    "DCN"
  )
)

marker_plot_list <- lapply(
  marker_panel,
  function(x) {

    safe_present(
      x,
      rownames(
        hepato_int
      )
    )
  }
)

marker_plot_list <- dedupe_feature_list(
  marker_plot_list
)

stopifnot(
  !anyDuplicated(
    unlist(
      marker_plot_list,
      use.names = FALSE
    )
  )
)

p06 <- DotPlot(
  hepato_int,
  features =
    marker_plot_list,
  assay = "RNA",
  group.by =
    "hepatobiliary_cluster_v132B",
  dot.scale = 7
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(
      -2.5,
      2.5
    ),
    oob =
      scales::squish
  ) +
  labs(
    title =
      "v1.3.2B hepatobiliary cluster marker validation",
    color =
      "Scaled\nexpression",
    size =
      "% expressed"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig06_HepatobiliaryCluster_Marker_DotPlot.pdf"
  ),
  p06,
  width = 23,
  height = 11
)

write.csv(
  p06$data,
  file.path(
    outdir,
    "16_HepatobiliaryCluster_Marker_DotPlot_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 26. Fig07 - sentinel FeaturePlots
# ------------------------------------------------------------

sentinel_genes <- c(
  "ALB",
  "CYP3A4",
  "CPS1",
  "GLUL",
  "AFP",
  "DLK1",
  "EPCAM",
  "PROM1",
  "KRT19",
  "KRT7",
  "SOX9",
  "MUC1",
  "CFTR",
  "KRT23",
  "SPP1",
  "MKI67",
  "PTPRC",
  "PECAM1",
  "COL1A1"
)

sentinel_present <- safe_present(
  sentinel_genes,
  rownames(
    hepato_int
  )
)

sentinel_availability <- data.frame(
  gene =
    sentinel_genes,
  present =
    sentinel_genes %in%
      rownames(
        hepato_int
      ),
  stringsAsFactors = FALSE
)

write.csv(
  sentinel_availability,
  file.path(
    outdir,
    "17_sentinel_gene_availability.csv"
  ),
  row.names = FALSE
)

if (
  length(
    sentinel_present
  ) > 0
) {

  p07 <- FeaturePlot(
    hepato_int,
    features =
      sentinel_present,
    reduction =
      "umapHepatoRPCA_v132B",
    ncol = 4,
    pt.size = 0.08,
    raster = TRUE,
    order = FALSE,
    cols = c(
      "#F2F2F2",
      "#FF1A1A"
    )
  )

  ggsave(
    file.path(
      outdir,
      "Fig07_SentinelMarker_FeaturePlots.pdf"
    ),
    p07,
    width = 16,
    height = 20
  )
}

# ------------------------------------------------------------
# 27. Fig08 - sample contribution within each res0.6 cluster
# ------------------------------------------------------------

p08 <- ggplot(
  cluster_sample,
  aes(
    x = sample,
    y = cluster,
    fill =
      within_cluster_pct
  )
) +
  geom_tile(
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        Freq > 0,
        sprintf(
          "%.1f",
          within_cluster_pct
        ),
        ""
      )
    ),
    size = 2.2
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Sample contribution within hepatobiliary res0.6 clusters",
    x = NULL,
    y =
      "Cluster",
    fill = "%"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig08_ClusterSampleContribution_QC.pdf"
  ),
  p08,
  width = 12,
  height = 9
)

# ------------------------------------------------------------
# 28. Fig09 - hepatobiliary composition by LT sample
# ------------------------------------------------------------

current_sample_tab <- as.data.frame(
  table(
    sample =
      hepato_int[[sample_col]][, 1],
    current_label =
      hepato_int$current_hepatobiliary_label_v132B
  ),
  stringsAsFactors = FALSE
)

current_sample_tab$sample_total <- ave(
  current_sample_tab$Freq,
  current_sample_tab$sample,
  FUN = sum
)

current_sample_tab$within_sample_pct <- ifelse(
  current_sample_tab$sample_total > 0,
  100 *
    current_sample_tab$Freq /
    current_sample_tab$sample_total,
  NA_real_
)

current_sample_tab$sample <- factor(
  current_sample_tab$sample,
  levels = sample_order
)

p09 <- ggplot(
  current_sample_tab,
  aes(
    x = sample,
    y = within_sample_pct,
    fill =
      current_label
  )
) +
  geom_col(
    width = 0.85
  ) +
  scale_fill_manual(
    values =
      current_palette
  ) +
  labs(
    title =
      "Current hepatobiliary annotation composition by LT sample",
    x = NULL,
    y =
      "% of hepatobiliary-core cells within sample",
    fill =
      "Current annotation"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig09_CurrentHepatobiliaryComposition_by_LT.pdf"
  ),
  p09,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 29. Per-cell metadata export
# ------------------------------------------------------------

cell_metadata <- data.frame(
  cell =
    colnames(
      hepato_int
    ),
  sample =
    hepato_int[[sample_col]][, 1],
  NAS_score =
    hepato_int[[nas_col]][, 1],
  NAS_group =
    hepato_int[[nasgrp_col]][, 1],
  current_label =
    hepato_int$current_hepatobiliary_label_v132B,
  cluster_res0.6 =
    hepato_int$hepatobiliary_cluster_v132B,
  stringsAsFactors = FALSE
)

write.csv(
  cell_metadata,
  file.path(
    outdir,
    "18_v1.3.2B_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 30. Review flags
# ------------------------------------------------------------

dominant_clusters <- cluster_dominance[
  cluster_dominance$max_sample_contribution_pct >= 70,
  ,
  drop = FALSE
]

review_flags <- c(
  "Human MASH v1.3.2B Hepatobiliary Refinement",
  "============================================",
  "",
  paste0(
    "Canonical cells: ",
    ncol(obj)
  ),
  paste0(
    "Hepatobiliary core cells: ",
    ncol(hepato_int)
  ),
  paste0(
    "Expected hepatobiliary core cells: ",
    expected_core_cells
  ),
  paste0(
    "Strict ductal-like rescue candidates outside core (audit only): ",
    sum(
      strict_rescue_candidate
    )
  ),
  paste0(
    "Annotation review resolution: ",
    annotation_resolution
  ),
  paste0(
    "Clusters at review resolution: ",
    length(
      unique(
        hepato_int$hepatobiliary_cluster_v132B
      )
    )
  ),
  "",
  "IMPORTANT INTERPRETATION RULES:",
  "  1. No final hepatobiliary labels are frozen at v1.3.2B.",
  "  2. Hepatocyte zonation must be separated from injury/transitional biology.",
  "  3. Fetal hepatoblast identity requires AFP/DLK1/IGF2/H19 support.",
  "  4. Ductular/progenitor identity requires EPCAM/PROM1/KRT19/KRT7/SOX9 support.",
  "  5. Bona fide cholangiocyte identity requires coherent biliary markers such as KRT7/KRT19/SOX9/MUC1/CFTR.",
  "  6. Strict outside-core ductal candidates are audit-only and are not automatically rescued.",
  "",
  "NEXT DECISION:",
  "  Review cluster markers, signature heatmap, DotPlot, FeaturePlots, crosswalk, and sample bias.",
  "  Then create v1.3.2B.1 manual hepatobiliary annotation/cleaning.",
  "",
  "STATUS:",
  "  PROVISIONAL_REVIEW_REQUIRED",
  "  Canonical whole-liver v1.3 object was NOT modified."
)

if (
  nrow(
    dominant_clusters
  ) > 0
) {

  review_flags <- c(
    review_flags,
    "",
    "FLAG: clusters with >=70% contribution from one LT sample:",
    paste0(
      "  cluster ",
      dominant_clusters$cluster,
      " : ",
      dominant_clusters$sample,
      " = ",
      round(
        dominant_clusters$max_sample_contribution_pct,
        1
      ),
      "%"
    )
  )
}

writeLines(
  review_flags,
  file.path(
    outdir,
    "19_v1.3.2B_review_flags.txt"
  )
)

cat(
  paste(
    review_flags,
    collapse = "\n"
  )
)

cat("\n\n")

# ------------------------------------------------------------
# 31. Save provisional hepatobiliary object
# ------------------------------------------------------------

provisional_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatobiliary_v1.3.2B_provisional.rds"
)

saveRDS(
  hepato_int,
  provisional_rds
)

# ------------------------------------------------------------
# 32. Manifest / provenance
# ------------------------------------------------------------

input_md5 <- unname(
  tools::md5sum(
    input_rds
  )
)

manifest <- c(
  "version=v1.3.2B",
  paste0(
    "input=",
    input_rds
  ),
  paste0(
    "input_md5=",
    input_md5
  ),
  paste0(
    "canonical_cells=",
    ncol(obj)
  ),
  paste0(
    "hepatobiliary_core_cells=",
    ncol(hepato_int)
  ),
  paste0(
    "strict_outside_core_ductal_candidates_audit_only=",
    sum(
      strict_rescue_candidate
    )
  ),
  paste0(
    "annotation_resolution=",
    annotation_resolution
  ),
  paste0(
    "seed=",
    seed
  ),
  paste0(
    "preplot_checkpoint=",
    preplot_checkpoint_rds
  ),
  paste0(
    "provisional_object=",
    provisional_rds
  ),
  "RPCA_recomputed=TRUE",
  "clustering_recomputed=TRUE",
  "manual_annotation_frozen=FALSE",
  "canonical_object_modified=FALSE",
  "annotation_status=PROVISIONAL_REVIEW_REQUIRED"
)

writeLines(
  manifest,
  file.path(
    outdir,
    "98_manifest.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)

cat("============================================================\n")
cat("v1.3.2B COMPLETE\n")
cat("HEPATOBILIARY REFINEMENT COMPLETE\n")
cat("ANNOTATION STATUS: PROVISIONAL_REVIEW_REQUIRED\n")
cat("CANONICAL v1.3 OBJECT NOT MODIFIED\n")
cat("\n")
cat("Pre-plot checkpoint:\n")
cat(preplot_checkpoint_rds, "\n")
cat("\n")
cat("Provisional object:\n")
cat(provisional_rds, "\n")
cat("\n")
cat("Output directory:\n")
cat(outdir, "\n")
cat("============================================================\n")
