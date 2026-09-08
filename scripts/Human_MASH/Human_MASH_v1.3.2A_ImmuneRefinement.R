# ============================================================
# Human MASH v1.3.2A
# Immune Refinement
#
# PURPOSE
#   Refine the immune compartment of the canonical Human MASH
#   LT1-LT8 whole-liver object before disease-progression and
#   IL10 mechanistic analyses.
#
# DESIGN
#   1. Read canonical v1.3 object (READ ONLY)
#   2. Define core immune cells from current REF12 annotation
#   3. Conservatively rescue additional PTPRC+ immune-marker+
#      cells from other current labels
#   4. Re-normalize by sample
#   5. RPCA integration across LT1-LT8
#   6. Multi-resolution clustering
#   7. Curated immune-lineage scoring
#   8. Provisional cluster annotation
#   9. Marker / sample-bias / crosswalk QC
#  10. Exploratory IL10-axis visualization
#
# IMPORTANT
#   - Canonical whole-liver RDS is NEVER overwritten.
#   - v1.3.2A annotations are PROVISIONAL.
#   - Final annotation freeze occurs only after visual review.
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
options(future.globals.maxSize = 8 * 1024^3)

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
  "v1.3.2A_immune_refinement"
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

core_immune_labels <- c(
  "Inflammatory_monocyte_macrophage",
  "SPP1_TREM2_macrophage",
  "Plasma_cell"
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
# Helper functions
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

  u <- unique(as.character(x))

  if (
    length(u) > 0 &&
    all(grepl("^LT[0-9]+", u))
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
        u[order(n)]
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

# ------------------------------------------------------------
# 1. Load canonical object
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.2A Immune Refinement\n")
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
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {

  stop(
    "Missing required metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

cat("Canonical cells :", ncol(obj), "\n")
cat("RNA genes       :", nrow(obj), "\n")
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
# 2. Define core immune compartment
# ------------------------------------------------------------

current_celltype <- as.character(
  obj[[celltype_col]][, 1]
)

current_sample <- as.character(
  obj[[sample_col]][, 1]
)

core_immune <- current_celltype %in%
  core_immune_labels

cat(
  "Core immune cells from current labels:",
  sum(core_immune),
  "\n"
)

# ------------------------------------------------------------
# 3. Conservative whole-liver immune rescue audit
#
# Rescue rule:
#   outside core labels
#   AND PTPRC > 0
#   AND >=1 additional immune-lineage marker > 0
#
# This intentionally favors specificity over sensitivity.
# ------------------------------------------------------------

rescue_secondary_genes <- c(
  # T
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",

  # NK
  "NKG7",
  "GNLY",
  "KLRD1",

  # B
  "MS4A1",
  "CD79A",
  "CD79B",

  # Plasma
  "JCHAIN",
  "MZB1",
  "SDC1",

  # Myeloid
  "LST1",
  "TYROBP",
  "FCER1G",
  "CTSS",

  # DC
  "FCER1A",
  "CD1C",
  "CLEC10A",
  "LILRA4"
)

rna_data_whole <- get_data_layer(
  obj,
  assay = "RNA"
)

gene_universe <- rownames(
  rna_data_whole
)

rescue_secondary_present <- safe_present(
  rescue_secondary_genes,
  gene_universe
)

ptprc_present <- "PTPRC" %in%
  gene_universe

if (
  ptprc_present &&
  length(rescue_secondary_present) > 0
) {

  ptprc_positive <- as.numeric(
    rna_data_whole[
      "PTPRC",
      ,
      drop = TRUE
    ]
  ) > 0

  secondary_hits <- Matrix::colSums(
    rna_data_whole[
      rescue_secondary_present,
      ,
      drop = FALSE
    ] > 0
  )

  rescued_immune <- (
    !core_immune &
    ptprc_positive &
    secondary_hits >= 1
  )

} else {

  warning(
    "PTPRC or rescue marker set unavailable; ",
    "no rescue cells added."
  )

  rescued_immune <- rep(
    FALSE,
    ncol(obj)
  )
}

immune_include <- (
  core_immune |
  rescued_immune
)

entry_reason <- rep(
  "not_included",
  ncol(obj)
)

entry_reason[
  core_immune
] <- "core_current_annotation"

entry_reason[
  rescued_immune
] <- "rescued_PTPRC_plus_lineage_marker"

cat(
  "Rescued additional immune cells:",
  sum(rescued_immune),
  "\n"
)

cat(
  "Total immune cells entering v1.3.2A:",
  sum(immune_include),
  "\n\n"
)

rescue_audit <- as.data.frame(
  table(
    current_celltype = current_celltype,
    entry_reason = entry_reason
  ),
  stringsAsFactors = FALSE
)

write.csv(
  rescue_audit,
  file.path(
    outdir,
    "01_immune_entry_audit_by_current_celltype.csv"
  ),
  row.names = FALSE
)

rescued_cells_table <- data.frame(
  cell = colnames(obj)[rescued_immune],
  current_celltype =
    current_celltype[rescued_immune],
  sample =
    current_sample[rescued_immune],
  stringsAsFactors = FALSE
)

write.csv(
  rescued_cells_table,
  file.path(
    outdir,
    "02_rescued_immune_cells.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Create immune subset
# ------------------------------------------------------------

immune_cells <- colnames(obj)[
  immune_include
]

immune <- subset(
  obj,
  cells = immune_cells
)

immune$immune_entry_reason_v132A <-
  entry_reason[immune_include]

immune$current_celltype_v132A <-
  current_celltype[immune_include]

DefaultAssay(immune) <- "RNA"

rm(
  rna_data_whole
)

gc()

# ------------------------------------------------------------
# 5. Sample-level QC before integration
# ------------------------------------------------------------

sample_order <- sample_order_fun(
  immune[[sample_col]][, 1]
)

sample_counts <- as.data.frame(
  table(
    sample =
      immune[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

names(sample_counts)[2] <- "n_cells"

sample_counts$sample <- factor(
  sample_counts$sample,
  levels = sample_order
)

sample_counts <- sample_counts[
  order(sample_counts$sample),
]

sample_counts$sample <- as.character(
  sample_counts$sample
)

write.csv(
  sample_counts,
  file.path(
    outdir,
    "03_immune_cells_by_sample.csv"
  ),
  row.names = FALSE
)

cat("=== IMMUNE CELLS BY SAMPLE ===\n")
print(
  sample_counts,
  row.names = FALSE
)

cat("\n")

if (
  any(
    sample_counts$n_cells < 50
  )
) {

  stop(
    "At least one LT sample has <50 immune cells; ",
    "RPCA integration aborted."
  )
}

# ------------------------------------------------------------
# 6. Split by LT sample and prepare RPCA integration
# ------------------------------------------------------------

cat("=== RPCA PREPARATION ===\n")

immune_list <- SplitObject(
  immune,
  split.by = sample_col
)

immune_list <- immune_list[
  sample_order
]

immune_list <- lapply(
  immune_list,
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
      nfeatures =
        nfeatures_integration,
      verbose = FALSE
    )

    x
  }
)

integration_features <- SelectIntegrationFeatures(
  object.list = immune_list,
  nfeatures = nfeatures_integration
)

immune_list <- lapply(
  immune_list,
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
# 7. RPCA anchors and integration
# ------------------------------------------------------------

cat("=== RPCA ANCHORS ===\n")

anchors <- FindIntegrationAnchors(
  object.list = immune_list,
  anchor.features =
    integration_features,
  reduction = "rpca",
  dims = 1:npcs,
  k.anchor = 5,
  verbose = FALSE
)

cat("=== INTEGRATE DATA ===\n")

immune_int <- IntegrateData(
  anchorset = anchors,
  dims = 1:npcs,
  new.assay.name = "integrated",
  verbose = FALSE
)

DefaultAssay(
  immune_int
) <- "integrated"

immune_int <- ScaleData(
  immune_int,
  verbose = FALSE
)

immune_int <- RunPCA(
  immune_int,
  npcs = npcs,
  verbose = FALSE
)

immune_int <- RunUMAP(
  immune_int,
  reduction = "pca",
  dims = 1:npcs,
  seed.use = seed,
  reduction.name =
    "umapImmuneRPCA_v132A",
  reduction.key =
    "UMAPIMM132A_",
  verbose = FALSE
)

immune_int <- FindNeighbors(
  immune_int,
  reduction = "pca",
  dims = 1:npcs,
  verbose = FALSE
)

# ------------------------------------------------------------
# 8. Multi-resolution clustering
# ------------------------------------------------------------

cat("=== MULTI-RESOLUTION CLUSTERING ===\n")

for (
  res in resolution_grid
) {

  set.seed(seed)

  immune_int <- FindClusters(
    immune_int,
    resolution = res,
    algorithm = 1,
    random.seed = seed,
    verbose = FALSE
  )

  metadata_name <- paste0(
    "immune_res",
    gsub(
      "\\.",
      "_",
      as.character(res)
    ),
    "_v132A"
  )

  immune_int[[
    metadata_name
  ]] <- as.character(
    Idents(immune_int)
  )
}

annotation_cluster_col <- paste0(
  "immune_res",
  gsub(
    "\\.",
    "_",
    as.character(
      annotation_resolution
    )
  ),
  "_v132A"
)

immune_int$immune_cluster_v132A <-
  as.character(
    immune_int[[
      annotation_cluster_col
    ]][, 1]
  )

Idents(
  immune_int
) <- "immune_cluster_v132A"

cluster_counts <- as.data.frame(
  table(
    cluster =
      immune_int$immune_cluster_v132A
  ),
  stringsAsFactors = FALSE
)

names(cluster_counts)[2] <-
  "n_cells"

write.csv(
  cluster_counts,
  file.path(
    outdir,
    "04_res0.6_cluster_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Curated immune marker signatures
# ------------------------------------------------------------

immune_signatures <- list(

  T_cell = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2",
    "LCK",
    "IL7R",
    "LTB"
  ),

  NK_cell = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "PRF1",
    "GZMB",
    "GZMH",
    "TRAC"
  ),

  B_cell = c(
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD37",
    "CD74",
    "HLA-DRA",
    "CD22"
  ),

  Plasma_cell = c(
    "JCHAIN",
    "MZB1",
    "SDC1",
    "DERL3",
    "IGHG1",
    "IGHG3",
    "IGKC"
  ),

  cDC = c(
    "FCER1A",
    "CD1C",
    "CLEC10A",
    "CST3",
    "CD74",
    "HLA-DRA",
    "HLA-DPA1"
  ),

  pDC = c(
    "LILRA4",
    "GZMB",
    "TCF4",
    "IRF7",
    "SERPINF1"
  ),

  Monocyte = c(
    "FCN1",
    "S100A8",
    "S100A9",
    "S100A10",
    "CTSS",
    "LILRB1",
    "VCAN"
  ),

  Resident_macrophage = c(
    "C1QA",
    "C1QB",
    "C1QC",
    "MARCO",
    "VSIG4",
    "CD5L",
    "TIMD4",
    "FOLR2"
  ),

  SPP1_TREM2_macrophage = c(
    "SPP1",
    "TREM2",
    "GPNMB",
    "LPL",
    "APOE",
    "CD9",
    "LGALS3"
  ),

  Cycling_immune = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "STMN1"
  )
)

DefaultAssay(
  immune_int
) <- "RNA"

# Join RNA layers if needed / supported
if (
  "JoinLayers" %in%
    getNamespaceExports(
      "SeuratObject"
    )
) {

  joined_rna <- try(
    SeuratObject::JoinLayers(
      immune_int[["RNA"]]
    ),
    silent = TRUE
  )

  if (
    !inherits(
      joined_rna,
      "try-error"
    )
  ) {

    immune_int[["RNA"]] <-
      joined_rna
  }
}

rna_data <- get_data_layer(
  immune_int,
  assay = "RNA"
)

rna_genes <- rownames(
  rna_data
)

signature_availability <- do.call(
  rbind,
  lapply(
    names(
      immune_signatures
    ),
    function(sig) {

      genes <-
        immune_signatures[[sig]]

      data.frame(
        signature = sig,
        gene = genes,
        present =
          genes %in% rna_genes,
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
    "05_signature_gene_availability.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Cell-level signature scores
#
# Simple mean log-normalized expression.
# No control-gene subtraction is used.
# This avoids AddModuleScore reference-bin dependence.
# ------------------------------------------------------------

score_columns <- character(0)

for (
  sig in names(
    immune_signatures
  )
) {

  genes_present <- safe_present(
    immune_signatures[[sig]],
    rna_genes
  )

  score_name <- paste0(
    "score_",
    sig,
    "_v132A"
  )

  if (
    length(
      genes_present
    ) == 0
  ) {

    immune_int[[
      score_name
    ]] <- 0

  } else {

    immune_int[[
      score_name
    ]] <-
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
# 11. Cluster-level signature scoring
# ------------------------------------------------------------

metadata_scoring <-
  immune_int@meta.data[
    ,
    c(
      "immune_cluster_v132A",
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
      metadata_scoring$immune_cluster_v132A
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
  "_v132A$",
  "",
  colnames(score_matrix)
)

# Standardize each signature across clusters
score_z <- scale(
  score_matrix,
  center = TRUE,
  scale = TRUE
)

score_z[
  !is.finite(score_z)
] <- 0

# ------------------------------------------------------------
# 12. Provisional cluster annotation
# ------------------------------------------------------------

annotation_rows <- vector(
  "list",
  nrow(score_z)
)

for (
  i in seq_len(
    nrow(score_z)
  )
) {

  x <- score_z[i, ]

  ord <- order(
    x,
    decreasing = TRUE
  )

  top_name <- names(x)[
    ord[1]
  ]

  second_name <- if (
    length(ord) >= 2
  ) {
    names(x)[ord[2]]
  } else {
    NA_character_
  }

  top_z <- x[
    ord[1]
  ]

  second_z <- if (
    length(ord) >= 2
  ) {
    x[ord[2]]
  } else {
    NA_real_
  }

  margin <- if (
    is.finite(second_z)
  ) {
    top_z - second_z
  } else {
    NA_real_
  }

  confidence <- if (
    is.finite(margin) &&
    top_z >= 0.5 &&
    margin >= 0.75
  ) {

    "high"

  } else if (
    is.finite(margin) &&
    top_z >= 0 &&
    margin >= 0.30
  ) {

    "moderate"

  } else {

    "low"
  }

  annotation_rows[[i]] <- data.frame(
    cluster =
      rownames(score_z)[i],
    provisional_celltype =
      top_name,
    top_signature =
      top_name,
    second_signature =
      second_name,
    top_z =
      round(top_z, 4),
    second_z =
      round(second_z, 4),
    margin =
      round(margin, 4),
    confidence =
      confidence,
    stringsAsFactors = FALSE
  )
}

cluster_annotation <- do.call(
  rbind,
  annotation_rows
)

write.csv(
  cluster_annotation,
  file.path(
    outdir,
    "06_cluster_provisional_annotation.csv"
  ),
  row.names = FALSE
)

score_z_out <- data.frame(
  cluster =
    rownames(score_z),
  score_z,
  check.names = FALSE,
  row.names = NULL
)

write.csv(
  score_z_out,
  file.path(
    outdir,
    "07_cluster_signature_Zscores.csv"
  ),
  row.names = FALSE
)

score_raw_out <- data.frame(
  cluster =
    rownames(score_matrix),
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

write.csv(
  score_raw_out,
  file.path(
    outdir,
    "08_cluster_signature_raw_scores.csv"
  ),
  row.names = FALSE
)

cluster_to_type <- setNames(
  cluster_annotation$provisional_celltype,
  cluster_annotation$cluster
)

cluster_to_conf <- setNames(
  cluster_annotation$confidence,
  cluster_annotation$cluster
)

cluster_to_margin <- setNames(
  cluster_annotation$margin,
  cluster_annotation$cluster
)

immune_int$immune_celltype_provisional_v132A <-
  unname(
    cluster_to_type[
      immune_int$immune_cluster_v132A
    ]
  )

immune_int$immune_annotation_confidence_v132A <-
  unname(
    cluster_to_conf[
      immune_int$immune_cluster_v132A
    ]
  )

immune_int$immune_annotation_margin_v132A <-
  as.numeric(
    unname(
      cluster_to_margin[
        immune_int$immune_cluster_v132A
      ]
    )
  )

# ------------------------------------------------------------
# 13. Current vs refined crosswalk
# ------------------------------------------------------------

crosswalk <- as.data.frame(
  table(
    current =
      immune_int$current_celltype_v132A,
    provisional =
      immune_int$immune_celltype_provisional_v132A
  ),
  stringsAsFactors = FALSE
)

write.csv(
  crosswalk,
  file.path(
    outdir,
    "09_current_vs_provisional_crosswalk.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Sample contribution QC
# ------------------------------------------------------------

sample_vec <- as.character(
  immune_int[[sample_col]][, 1]
)

provisional_vec <- as.character(
  immune_int$immune_celltype_provisional_v132A
)

sample_order <- sample_order_fun(
  sample_vec
)

type_sample <- as.data.frame(
  table(
    celltype =
      provisional_vec,
    sample =
      sample_vec
  ),
  stringsAsFactors = FALSE
)

type_sample$celltype_total <- ave(
  type_sample$Freq,
  type_sample$celltype,
  FUN = sum
)

type_sample$sample_total <- ave(
  type_sample$Freq,
  type_sample$sample,
  FUN = sum
)

type_sample$within_celltype_pct <- ifelse(
  type_sample$celltype_total > 0,
  100 *
    type_sample$Freq /
    type_sample$celltype_total,
  NA_real_
)

type_sample$within_sample_pct <- ifelse(
  type_sample$sample_total > 0,
  100 *
    type_sample$Freq /
    type_sample$sample_total,
  NA_real_
)

type_sample$sample <- factor(
  type_sample$sample,
  levels = sample_order
)

write.csv(
  type_sample,
  file.path(
    outdir,
    "10_provisional_celltype_by_sample.csv"
  ),
  row.names = FALSE
)

# Cluster-level sample dominance
cluster_sample <- as.data.frame(
  table(
    cluster =
      immune_int$immune_cluster_v132A,
    sample =
      sample_vec
  ),
  stringsAsFactors = FALSE
)

cluster_sample$total <- ave(
  cluster_sample$Freq,
  cluster_sample$cluster,
  FUN = sum
)

cluster_sample$within_cluster_pct <- ifelse(
  cluster_sample$total > 0,
  100 *
    cluster_sample$Freq /
    cluster_sample$total,
  NA_real_
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
    "11_cluster_max_sample_contribution.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Cluster marker discovery
# ------------------------------------------------------------

DefaultAssay(
  immune_int
) <- "RNA"

Idents(
  immune_int
) <- "immune_cluster_v132A"

cat("=== FIND ALL CLUSTER MARKERS ===\n")

cluster_markers <- FindAllMarkers(
  immune_int,
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
    "12_immune_cluster_markers_all.csv"
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
          length(fc_col) > 0
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
          15
        )
      }
    )
  )

  write.csv(
    top_markers,
    file.path(
      outdir,
      "13_immune_cluster_markers_top15.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 15A. Pre-plot checkpoint
#
# Saves the completed RPCA / clustering / provisional
# annotation / marker-analysis state before figure generation.
# If a later plotting step fails, this object can be reused
# without repeating RPCA integration and FindAllMarkers().
# ------------------------------------------------------------

preplot_checkpoint_rds <- file.path(
  object_dir,
  "Human_MASH_Immune_v1.3.2A_preplot_checkpoint.rds"
)

saveRDS(
  immune_int,
  preplot_checkpoint_rds
)

writeLines(
  c(
    "Human MASH v1.3.2A pre-plot checkpoint",
    paste0(
      "cells=",
      ncol(immune_int)
    ),
    paste0(
      "annotation_resolution=",
      annotation_resolution
    ),
    paste0(
      "clusters=",
      length(
        unique(
          immune_int$immune_cluster_v132A
        )
      )
    ),
    paste0(
      "checkpoint_rds=",
      preplot_checkpoint_rds
    ),
    "RPCA_complete=TRUE",
    "clustering_complete=TRUE",
    "provisional_annotation_complete=TRUE",
    "FindAllMarkers_complete=TRUE"
  ),
  file.path(
    outdir,
    "13A_preplot_checkpoint.txt"
  )
)

cat(
  "PRE-PLOT CHECKPOINT SAVED:\n",
  preplot_checkpoint_rds,
  "\n\n"
)

# ------------------------------------------------------------
# 16. Visualization palettes
# ------------------------------------------------------------

immune_palette <- c(
  "T_cell" =
    "#6A3D9A",

  "NK_cell" =
    "#1F78B4",

  "B_cell" =
    "#00A650",

  "Plasma_cell" =
    "#B5179E",

  "cDC" =
    "#00A6A6",

  "pDC" =
    "#00CFEA",

  "Monocyte" =
    "#E31A1C",

  "Resident_macrophage" =
    "#8B0000",

  "SPP1_TREM2_macrophage" =
    "#FF7F00",

  "Cycling_immune" =
    "#FFD600"
)

cluster_palette <- grDevices::hcl.colors(
  length(
    unique(
      immune_int$immune_cluster_v132A
    )
  ),
  palette = "Dark 3"
)

# ------------------------------------------------------------
# 17. Fig01 - current annotation on refined UMAP
# ------------------------------------------------------------

p01 <- DimPlot(
  immune_int,
  reduction =
    "umapImmuneRPCA_v132A",
  group.by =
    "current_celltype_v132A",
  pt.size = 0.18,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  ggtitle(
    "v1.3.2A immune subset: current REF12 annotation"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig01_CurrentAnnotation_ImmuneRPCA_UMAP.pdf"
  ),
  p01,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 18. Fig02 - res0.6 clusters
# ------------------------------------------------------------

p02 <- DimPlot(
  immune_int,
  reduction =
    "umapImmuneRPCA_v132A",
  group.by =
    "immune_cluster_v132A",
  pt.size = 0.18,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    cluster_palette
) +
  ggtitle(
    paste0(
      "Immune RPCA clusters - resolution ",
      annotation_resolution
    )
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig02_ImmuneClusters_res0.6_UMAP.pdf"
  ),
  p02,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 19. Fig03 - provisional annotation
# ------------------------------------------------------------

p03 <- DimPlot(
  immune_int,
  reduction =
    "umapImmuneRPCA_v132A",
  group.by =
    "immune_celltype_provisional_v132A",
  pt.size = 0.18,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    immune_palette
) +
  ggtitle(
    "Provisional immune lineage annotation v1.3.2A"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig03_ProvisionalImmuneAnnotation_UMAP.pdf"
  ),
  p03,
  width = 11,
  height = 8
)

# ------------------------------------------------------------
# 20. Fig04 - provisional annotation split by LT sample
# ------------------------------------------------------------

p04 <- DimPlot(
  immune_int,
  reduction =
    "umapImmuneRPCA_v132A",
  group.by =
    "immune_celltype_provisional_v132A",
  split.by = sample_col,
  ncol = 4,
  pt.size = 0.12,
  raster = TRUE,
  cols =
    immune_palette
)

ggsave(
  file.path(
    outdir,
    "Fig04_ProvisionalImmuneAnnotation_split_by_LT.pdf"
  ),
  p04,
  width = 18,
  height = 12
)

# ------------------------------------------------------------
# 21. Fig05 - curated immune marker DotPlot
# ------------------------------------------------------------

marker_list_present <- lapply(
  immune_signatures,
  function(x) {

    safe_present(
      x,
      rownames(
        immune_int
      )
    )
  }
)

# Seurat DotPlot requires globally unique feature levels.
# Shared biological markers remain in the signature-score analysis above,
# but are shown only once in this visualization.
marker_plot_long <- data.frame(
  signature = rep(
    names(marker_list_present),
    lengths(marker_list_present)
  ),
  gene = unlist(
    marker_list_present,
    use.names = FALSE
  ),
  stringsAsFactors = FALSE
)

marker_plot_duplicates <- marker_plot_long[
  duplicated(marker_plot_long$gene) |
    duplicated(marker_plot_long$gene, fromLast = TRUE),
  ,
  drop = FALSE
]

write.csv(
  marker_plot_duplicates,
  file.path(
    outdir,
    "14A_DotPlot_duplicate_marker_audit.csv"
  ),
  row.names = FALSE
)

seen_marker_genes <- character(0)

marker_list_present <- lapply(
  marker_list_present,
  function(x) {

    keep <- x[
      !x %in% seen_marker_genes
    ]

    seen_marker_genes <<- c(
      seen_marker_genes,
      keep
    )

    keep
  }
)

marker_list_present <-
  marker_list_present[
    lengths(
      marker_list_present
    ) > 0
  ]

stopifnot(
  !anyDuplicated(
    unlist(
      marker_list_present,
      use.names = FALSE
    )
  )
)

p05 <- DotPlot(
  immune_int,
  features =
    marker_list_present,
  assay = "RNA",
  group.by =
    "immune_celltype_provisional_v132A",
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
      "v1.3.2A provisional immune annotation marker validation",
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
    "Fig05_ProvisionalImmune_Marker_DotPlot.pdf"
  ),
  p05,
  width = 22,
  height = 9
)

write.csv(
  p05$data,
  file.path(
    outdir,
    "14_ProvisionalImmune_DotPlot_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 22. Fig06 - cluster signature Z-score heatmap
# ------------------------------------------------------------

z_long <- data.frame(
  cluster =
    rep(
      rownames(score_z),
      times = ncol(score_z)
    ),
  signature =
    rep(
      colnames(score_z),
      each = nrow(score_z)
    ),
  z =
    as.vector(score_z),
  stringsAsFactors = FALSE
)

p06 <- ggplot(
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
      "Cluster-level immune lineage signature scores",
    x = NULL,
    y = "Cluster",
    fill = "Z score"
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
    "Fig06_ClusterSignature_Zscore_Heatmap.pdf"
  ),
  p06,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 23. Fig07 - sample contribution by provisional type
# ------------------------------------------------------------

p07 <- ggplot(
  type_sample,
  aes(
    x = sample,
    y = celltype,
    fill = within_celltype_pct
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
          within_celltype_pct
        ),
        ""
      )
    ),
    size = 2.5
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Sample contribution within provisional immune cell types",
    x = NULL,
    y = NULL,
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
    "Fig07_SampleContribution_ProvisionalImmuneTypes.pdf"
  ),
  p07,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 24. Fig08 - current vs provisional annotation crosswalk
# ------------------------------------------------------------

p08 <- ggplot(
  crosswalk,
  aes(
    x = provisional,
    y = current,
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
    size = 2.3
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Current REF12 annotation vs provisional v1.3.2A annotation",
    x =
      "Provisional immune type",
    y =
      "Current REF12 type",
    fill =
      "log(1+n)"
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
    "Fig08_Current_vs_Provisional_Crosswalk.pdf"
  ),
  p08,
  width = 13,
  height = 7
)

# ------------------------------------------------------------
# 25. Fig09 - exploratory IL10 axis after immune refinement
# ------------------------------------------------------------

il10_genes <- c(
  "IL10",
  "IL10RA",
  "IL10RB",
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3"
)

il10_present <- safe_present(
  il10_genes,
  rownames(
    immune_int
  )
)

il10_availability <- data.frame(
  gene = il10_genes,
  present =
    il10_genes %in%
      rownames(
        immune_int
      ),
  stringsAsFactors = FALSE
)

write.csv(
  il10_availability,
  file.path(
    outdir,
    "15_IL10_axis_gene_availability.csv"
  ),
  row.names = FALSE
)

if (
  length(
    il10_present
  ) > 0
) {

  p09 <- DotPlot(
    immune_int,
    features =
      il10_present,
    assay = "RNA",
    group.by =
      "immune_celltype_provisional_v132A",
    dot.scale = 8
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
        "Exploratory IL10 / IL10R axis after immune refinement",
      subtitle =
        "Do not interpret mechanistically until v1.3.2A annotation is reviewed",
      color =
        "Scaled\nexpression",
      size =
        "% expressed"
    ) +
    theme_bw(
      base_size = 10
    )

  ggsave(
    file.path(
      outdir,
      "Fig09_Exploratory_IL10axis_RefinedImmune_DotPlot.pdf"
    ),
    p09,
    width = 11,
    height = 7
  )

  write.csv(
    p09$data,
    file.path(
      outdir,
      "16_IL10axis_RefinedImmune_DotPlot_data.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 26. Fig10 - IL10 receptor / response FeaturePlots
# ------------------------------------------------------------

il10_feature_genes <- safe_present(
  c(
    "IL10RA",
    "IL10RB",
    "STAT3",
    "SOCS3"
  ),
  rownames(
    immune_int
  )
)

if (
  length(
    il10_feature_genes
  ) > 0
) {

  p10 <- FeaturePlot(
    immune_int,
    features =
      il10_feature_genes,
    reduction =
      "umapImmuneRPCA_v132A",
    ncol = 2,
    pt.size = 0.15,
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
      "Fig10_IL10R_STAT3_SOCS3_FeaturePlots.pdf"
    ),
    p10,
    width = 11,
    height = 9
  )
}

# ------------------------------------------------------------
# 27. Fig11 - provisional immune composition by LT sample
# ------------------------------------------------------------

p11 <- ggplot(
  type_sample,
  aes(
    x = sample,
    y = within_sample_pct,
    fill = celltype
  )
) +
  geom_col(
    width = 0.85
  ) +
  scale_fill_manual(
    values =
      immune_palette
  ) +
  labs(
    title =
      "Provisional immune composition by LT sample",
    x = NULL,
    y =
      "% of immune cells within sample",
    fill =
      "Provisional immune type"
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
    "Fig11_ProvisionalImmuneComposition_by_LT.pdf"
  ),
  p11,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 28. Cell metadata export
# ------------------------------------------------------------

cell_metadata_export <- data.frame(
  cell =
    colnames(
      immune_int
    ),

  sample =
    immune_int[[sample_col]][, 1],

  NAS_score =
    immune_int[[nas_col]][, 1],

  NAS_group =
    immune_int[[nasgrp_col]][, 1],

  current_celltype =
    immune_int$current_celltype_v132A,

  entry_reason =
    immune_int$immune_entry_reason_v132A,

  immune_cluster_v132A =
    immune_int$immune_cluster_v132A,

  immune_celltype_provisional_v132A =
    immune_int$immune_celltype_provisional_v132A,

  annotation_confidence =
    immune_int$immune_annotation_confidence_v132A,

  annotation_margin =
    immune_int$immune_annotation_margin_v132A,

  stringsAsFactors = FALSE
)

write.csv(
  cell_metadata_export,
  file.path(
    outdir,
    "17_v1.3.2A_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 29. Automated review flags
# ------------------------------------------------------------

flags <- c(
  "Human MASH v1.3.2A Immune Refinement",
  "====================================",
  "",
  paste0(
    "Canonical cells: ",
    ncol(obj)
  ),
  paste0(
    "Core immune cells: ",
    sum(core_immune)
  ),
  paste0(
    "Rescued immune cells: ",
    sum(rescued_immune)
  ),
  paste0(
    "Total immune cells analyzed: ",
    ncol(immune_int)
  ),
  paste0(
    "Annotation resolution: ",
    annotation_resolution
  ),
  paste0(
    "Clusters at annotation resolution: ",
    length(
      unique(
        immune_int$immune_cluster_v132A
      )
    )
  ),
  ""
)

low_conf <- cluster_annotation[
  cluster_annotation$confidence == "low",
]

if (
  nrow(
    low_conf
  ) > 0
) {

  flags <- c(
    flags,
    "FLAG: low-confidence provisional clusters:",
    paste0(
      "  cluster ",
      low_conf$cluster,
      " -> ",
      low_conf$provisional_celltype,
      " ; top=",
      low_conf$top_signature,
      " ; second=",
      low_conf$second_signature,
      " ; margin=",
      low_conf$margin
    ),
    ""
  )
}

dominant_clusters <- cluster_dominance[
  cluster_dominance$max_sample_contribution_pct >= 70,
]

if (
  nrow(
    dominant_clusters
  ) > 0
) {

  flags <- c(
    flags,
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
    ),
    ""
  )
}

flags <- c(
  flags,
  "IMPORTANT:",
  "Provisional cluster labels are marker-score assisted, not final.",
  "Review Fig02, Fig03, Fig05, Fig06, Fig07 and cluster marker tables.",
  "Do not merge these labels into the canonical whole-liver object yet.",
  "",
  "IL10-axis figures are exploratory at v1.3.2A.",
  "Final IL10 mechanistic testing starts only after immune annotation freeze.",
  "",
  "Canonical v1.3 RDS was not modified."
)

writeLines(
  flags,
  file.path(
    outdir,
    "18_v1.3.2A_review_flags.txt"
  )
)

cat(
  paste(
    flags,
    collapse = "\n"
  )
)

cat("\n\n")

# ------------------------------------------------------------
# 30. Save provisional immune object
# ------------------------------------------------------------

immune_rds_out <- file.path(
  object_dir,
  "Human_MASH_Immune_v1.3.2A_provisional.rds"
)

saveRDS(
  immune_int,
  immune_rds_out
)

# ------------------------------------------------------------
# 31. Manifest
# ------------------------------------------------------------

manifest <- c(
  paste0(
    "version=v1.3.2A"
  ),
  paste0(
    "input=",
    input_rds
  ),
  paste0(
    "output_object=",
    immune_rds_out
  ),
  paste0(
    "canonical_cells=",
    ncol(obj)
  ),
  paste0(
    "core_immune_cells=",
    sum(core_immune)
  ),
  paste0(
    "rescued_immune_cells=",
    sum(rescued_immune)
  ),
  paste0(
    "immune_cells_final=",
    ncol(immune_int)
  ),
  paste0(
    "annotation_resolution=",
    annotation_resolution
  ),
  paste0(
    "seed=",
    seed
  ),
  "canonical_object_modified=FALSE",
  "annotation_status=PROVISIONAL"
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
cat("v1.3.2A COMPLETE\n")
cat("PROVISIONAL IMMUNE ANNOTATION ONLY\n")
cat("CANONICAL v1.3 OBJECT NOT MODIFIED\n")
cat("Output object:\n")
cat(immune_rds_out, "\n")
cat("Output directory:\n")
cat(outdir, "\n")
cat("============================================================\n")
