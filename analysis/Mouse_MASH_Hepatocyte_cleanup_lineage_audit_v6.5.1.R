#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6510)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Hepatocyte cleanup + lineage/QC audit + clean re-clustering
#
# Version: v6.5.1
#
# INPUT:
#   v6.5.0 Hepatocyte subclustered RDS
#
# HIGH-CONFIDENCE SOURCE-CLUSTER DECISIONS FROM v6.5.0:
#
#   RETAIN:
#     0, 1, 2, 4
#
#   RETAIN BUT QC REVIEW:
#     5
#       mitochondrial-high hepatocyte candidate
#
#   RETAIN BUT LINEAGE REVIEW:
#     10, 12
#       cycling cells; lineage not fixed
#
#   REMOVE AS HIGH-CONFIDENCE NON-HEPATOCYTE CONTAMINATION:
#     3  = HSC / mesenchymal-like
#     6  = cholangiocyte / ductular-like
#     7  = macrophage / immune-like
#     8  = immune-like
#     9  = LSEC / endothelial-like
#     11 = HSC-like
#
# PURPOSE:
#   1) Remove only high-confidence contaminants.
#   2) Audit cluster 5 using nCount_RNA / nFeature_RNA / percent.mt.
#   3) Audit cycling source clusters 10/12 using lineage identity programs.
#   4) Re-run sample-aware RPCA after high-confidence cleanup.
#   5) Produce a cleaner hepatocyte map for final annotation in v6.5.2.
#
# IMPORTANT:
#   - Source clusters 10/12 are NOT automatically removed.
#   - No percent.mt threshold is used to remove cluster 5 in this version.
#   - Final hepatocyte state annotation is NOT fixed in v6.5.1.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}

save_pdf <- function(
  p,
  file,
  width,
  height
) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}

safe_join_rna <- function(object) {
  if (
    "RNA" %in%
      Assays(
        object
      ) &&
    length(
      Layers(
        object[["RNA"]]
      )
    ) >
      1
  ) {
    object[["RNA"]] <- JoinLayers(
      object[["RNA"]]
    )
  }

  object
}

present_genes <- function(
  object,
  genes
) {
  intersect(
    unique(
      genes
    ),
    rownames(
      object
    )
  )
}

add_program_score <- function(
  object,
  genes,
  name
) {

  genes_use <- present_genes(
    object,
    genes
  )

  final_name <- paste0(
    name,
    "_score"
  )

  if (
    length(
      genes_use
    ) <
      2
  ) {
    warning(
      "Program ",
      name,
      " has fewer than 2 genes present: ",
      paste(
        genes_use,
        collapse = ", "
      )
    )

    object[[
      final_name
    ]] <- NA_real_

    return(
      list(
        object = object,
        genes = genes_use
      )
    )
  }

  generated_prefix <- paste0(
    name,
    "_"
  )

  object <- AddModuleScore(
    object,
    features = list(
      genes_use
    ),
    name =
      generated_prefix,
    assay =
      "RNA",
    seed =
      6510
  )

  generated_name <- paste0(
    generated_prefix,
    "1"
  )

  object[[
    final_name
  ]] <- object[[
    generated_name
  ]][
    ,
    1
  ]

  object[[
    generated_name
  ]] <- NULL

  list(
    object = object,
    genes = genes_use
  )
}

cluster_fraction_table <- function(
  object,
  cluster_col,
  sample_col,
  samples
) {

  object@meta.data %>%
    as_tibble(
      rownames = "cell"
    ) %>%
    transmute(
      cell,
      sample =
        as.character(
          .data[[
            sample_col
          ]]
        ),
      cluster =
        as.character(
          .data[[
            cluster_col
          ]]
        )
    ) %>%
    filter(
      sample %in%
        samples
    ) %>%
    count(
      sample,
      cluster,
      name =
        "n_cells"
    ) %>%
    complete(
      sample =
        samples,
      cluster,
      fill = list(
        n_cells = 0
      )
    ) %>%
    group_by(
      sample
    ) %>%
    mutate(
      sample_total =
        sum(
          n_cells
        ),
      fraction =
        n_cells /
          sample_total
    ) %>%
    ungroup()
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_subclustering_v6.5.0",
  "RDS",
  "Mouse_MASH_Hepatocyte_subclustered_v6.5.0.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_cleanup_lineage_audit_v6.5.1"
)

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

for (
  d in c(
    OUT,
    RDS_OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

SOURCE_CLUSTER_COL <-
  "hep_rpca_res_0.4"

SAMPLE_COL <-
  "sample_hep_v650"

CONDITION_COL <-
  "condition_hep_v650"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

HIGH_CONF_REMOVE <- c(
  "3",
  "6",
  "7",
  "8",
  "9",
  "11"
)

CORE_RETAIN <- c(
  "0",
  "1",
  "2",
  "4"
)

QC_REVIEW <- c(
  "5"
)

LINEAGE_REVIEW <- c(
  "10",
  "12"
)

CLEAN_CANDIDATE_CLUSTERS <- c(
  CORE_RETAIN,
  QC_REVIEW,
  LINEAGE_REVIEW
)

NFEATURES <-
  3000

NPCS <-
  40

DIMS_USE <-
  1:30

RESOLUTIONS <- c(
  0.2,
  0.4,
  0.6
)

WORKING_RESOLUTION <-
  0.4

CLEAN_CLUSTER_COL <-
  "hepclean_rpca_res_0.4"

CLEAN_UMAP <-
  "umap.hep.clean.rpca"

CLEAN_UMAP_KEY <-
  "HEPCLEAN_"


# ==============================================================================
# 4. Identity and biological programs
# ==============================================================================

IDENTITY_PROGRAMS <- list(

  Hepatocyte_identity = c(
    "Alb",
    "Ttr",
    "Apoa1",
    "Apoa2",
    "Hnf4a",
    "Cps1",
    "Ass1",
    "G6pc",
    "Pck1"
  ),

  Cholangiocyte_identity = c(
    "Krt19",
    "Krt7",
    "Epcam",
    "Sox9",
    "Muc1",
    "Krt20",
    "Krt8",
    "Krt18"
  ),

  Endothelial_identity = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Vwf",
    "Ptprb",
    "Clec4g"
  ),

  Immune_identity = c(
    "Ptprc",
    "Lyz2",
    "Adgre1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Cd68"
  ),

  HSC_mesenchymal_identity = c(
    "Lrat",
    "Rbp1",
    "Pdgfra",
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Dcn",
    "Lum"
  )
)

BIO_PROGRAMS <- list(

  Periportal_Z1 = c(
    "Cps1",
    "Ass1",
    "Asl",
    "Arg1",
    "Pck1",
    "G6pc",
    "Gls2",
    "Hal"
  ),

  Pericentral_Z3 = c(
    "Glul",
    "Cyp2e1",
    "Cyp1a2",
    "Oat",
    "Axin2",
    "Slc1a2",
    "Cyp2a5"
  ),

  Injury_inflammatory = c(
    "Areg",
    "Ccl2",
    "Cxcl10",
    "Rsad2",
    "Lif",
    "Nupr1",
    "Sprr1a",
    "Socs3",
    "Nfkbia"
  ),

  ER_UPR_stress = c(
    "Ddit3",
    "Atf3",
    "Atf4",
    "Hspa5",
    "Xbp1",
    "Hsp90b1",
    "Herpud1",
    "Trib3"
  ),

  Oxidative_stress = c(
    "Hmox1",
    "Nqo1",
    "Gclc",
    "Gclm",
    "Txnrd1",
    "Srxn1"
  ),

  Regenerative_cycling = c(
    "Mki67",
    "Top2a",
    "Pcna",
    "Mcm2",
    "Mcm3",
    "Mcm4",
    "Mcm5",
    "Mcm6",
    "Mcm7",
    "Ccna2",
    "Ccnb1",
    "Cdk1"
  )
)


# ==============================================================================
# 5. Load v6.5.0
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input v6.5.0 RDS missing: ",
    INPUT_RDS
  )
}

msg(
  "Loading v6.5.0 hepatocyte object..."
)

hep <- readRDS(
  INPUT_RDS
)

DefaultAssay(
  hep
) <- "RNA"

hep <- safe_join_rna(
  hep
)

required_meta <- c(
  SOURCE_CLUSTER_COL,
  SAMPLE_COL,
  CONDITION_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    hep@meta.data
  )
)

if (
  length(
    missing_meta
  )
) {
  stop(
    "Missing metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

hep$source_cluster_v650 <- as.character(
  hep@meta.data[[
    SOURCE_CLUSTER_COL
  ]]
)

hep$cleanup_status_v651 <- case_when(
  hep$source_cluster_v650 %in%
    HIGH_CONF_REMOVE ~
    "Remove_high_conf_contaminant",

  hep$source_cluster_v650 %in%
    CORE_RETAIN ~
    "Retain_core_hepatocyte",

  hep$source_cluster_v650 %in%
    QC_REVIEW ~
    "Retain_QC_review",

  hep$source_cluster_v650 %in%
    LINEAGE_REVIEW ~
    "Retain_lineage_review",

  TRUE ~
    "Unexpected_source_cluster"
)

if (
  any(
    hep$cleanup_status_v651 ==
      "Unexpected_source_cluster"
  )
) {
  unexpected <- sort(
    unique(
      hep$source_cluster_v650[
        hep$cleanup_status_v651 ==
          "Unexpected_source_cluster"
      ]
    )
  )

  stop(
    "Unexpected v6.5.0 source cluster(s): ",
    paste(
      unexpected,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 6. Ensure QC metadata
# ==============================================================================

msg(
  "Computing QC metrics..."
)

if (
  !"percent.mt" %in%
    colnames(
      hep@meta.data
    )
) {
  hep[[
    "percent.mt"
  ]] <- PercentageFeatureSet(
    hep,
    pattern =
      "^mt-",
    assay =
      "RNA"
  )
}

qc_required <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt"
)

qc_missing <- setdiff(
  qc_required,
  colnames(
    hep@meta.data
  )
)

if (
  length(
    qc_missing
  )
) {
  stop(
    "QC metadata missing: ",
    paste(
      qc_missing,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 7. Source-cluster cleanup audit
# ==============================================================================

source_cleanup_counts <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    condition =
      as.character(
        .data[[
          CONDITION_COL
        ]]
      ),
    source_cluster =
      source_cluster_v650,
    cleanup_status =
      cleanup_status_v651
  ) %>%
  count(
    sample,
    condition,
    source_cluster,
    cleanup_status,
    name =
      "n_cells"
  ) %>%
  arrange(
    factor(
      sample,
      levels =
        SAMPLES
    ),
    as.numeric(
      source_cluster
    )
  )

write.csv(
  source_cleanup_counts,
  file.path(
    TAB_OUT,
    "01_source_cluster_cleanup_counts_v6.5.1.csv"
  ),
  row.names = FALSE
)

cleanup_by_sample <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    cleanup_status =
      cleanup_status_v651
  ) %>%
  count(
    sample,
    cleanup_status,
    name =
      "n_cells"
  ) %>%
  group_by(
    sample
  ) %>%
  mutate(
    total_source_hepatocyte =
      sum(
        n_cells
      ),
    fraction =
      n_cells /
        total_source_hepatocyte
  ) %>%
  ungroup()

write.csv(
  cleanup_by_sample,
  file.path(
    TAB_OUT,
    "02_cleanup_status_by_sample_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. QC audit by source cluster and sample
# ==============================================================================

qc_by_cluster_sample <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    source_cluster =
      source_cluster_v650,
    cleanup_status =
      cleanup_status_v651,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  ) %>%
  group_by(
    sample,
    source_cluster,
    cleanup_status
  ) %>%
  summarise(
    n_cells =
      n(),

    median_nCount_RNA =
      median(
        nCount_RNA,
        na.rm = TRUE
      ),

    median_nFeature_RNA =
      median(
        nFeature_RNA,
        na.rm = TRUE
      ),

    median_percent_mt =
      median(
        percent.mt,
        na.rm = TRUE
      ),

    mean_percent_mt =
      mean(
        percent.mt,
        na.rm = TRUE
      ),

    pct_cells_mt_gt_10 =
      mean(
        percent.mt >
          10,
        na.rm = TRUE
      ),

    pct_cells_mt_gt_20 =
      mean(
        percent.mt >
          20,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  qc_by_cluster_sample,
  file.path(
    TAB_OUT,
    "03_source_cluster_QC_by_sample_v6.5.1.csv"
  ),
  row.names = FALSE
)

qc_by_cluster <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  transmute(
    source_cluster =
      source_cluster_v650,
    cleanup_status =
      cleanup_status_v651,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  ) %>%
  group_by(
    source_cluster,
    cleanup_status
  ) %>%
  summarise(
    n_cells =
      n(),

    median_nCount_RNA =
      median(
        nCount_RNA,
        na.rm = TRUE
      ),

    median_nFeature_RNA =
      median(
        nFeature_RNA,
        na.rm = TRUE
      ),

    median_percent_mt =
      median(
        percent.mt,
        na.rm = TRUE
      ),

    mean_percent_mt =
      mean(
        percent.mt,
        na.rm = TRUE
      ),

    pct_cells_mt_gt_10 =
      mean(
        percent.mt >
          10,
        na.rm = TRUE
      ),

    pct_cells_mt_gt_20 =
      mean(
        percent.mt >
          20,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write.csv(
  qc_by_cluster,
  file.path(
    TAB_OUT,
    "04_source_cluster_QC_summary_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Identity scoring BEFORE cleanup
# ==============================================================================

msg(
  "Scoring lineage identity before cleanup..."
)

identity_gene_audit <- list()

for (
  nm in names(
    IDENTITY_PROGRAMS
  )
) {

  ans <- add_program_score(
    hep,
    IDENTITY_PROGRAMS[[
      nm
    ]],
    paste0(
      "v651_",
      nm
    )
  )

  hep <- ans$object

  identity_gene_audit[[
    nm
  ]] <- tibble(
    program =
      nm,
    gene =
      IDENTITY_PROGRAMS[[
        nm
      ]],
    present =
      IDENTITY_PROGRAMS[[
        nm
      ]] %in%
        ans$genes
  )
}

identity_gene_audit <- bind_rows(
  identity_gene_audit
)

write.csv(
  identity_gene_audit,
  file.path(
    TAB_OUT,
    "05_lineage_identity_gene_audit_v6.5.1.csv"
  ),
  row.names = FALSE
)

identity_score_cols <- paste0(
  "v651_",
  names(
    IDENTITY_PROGRAMS
  ),
  "_score"
)

identity_by_source_cluster <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    source_cluster =
      source_cluster_v650
  ) %>%
  group_by(
    source_cluster
  ) %>%
  summarise(
    n_cells =
      n(),
    across(
      all_of(
        identity_score_cols
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write.csv(
  identity_by_source_cluster,
  file.path(
    TAB_OUT,
    "06_lineage_identity_scores_by_source_cluster_v6.5.1.csv"
  ),
  row.names = FALSE
)

review_identity <- identity_by_source_cluster %>%
  filter(
    source_cluster %in%
      c(
        QC_REVIEW,
        LINEAGE_REVIEW
      )
  )

write.csv(
  review_identity,
  file.path(
    TAB_OUT,
    "07_review_clusters_5_10_12_identity_scores_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Direct marker audit for review clusters 5 / 10 / 12
# ==============================================================================

REVIEW_DIRECT_GENES <- c(
  # hepatocyte
  "Alb",
  "Ttr",
  "Apoa1",
  "Apoa2",
  "Hnf4a",
  "Cps1",
  "Ass1",

  # cholangiocyte / ductular
  "Krt19",
  "Krt7",
  "Epcam",
  "Sox9",
  "Muc1",
  "Krt20",
  "Krt8",
  "Krt18",

  # mesenchymal
  "Lrat",
  "Rbp1",
  "Pdgfra",
  "Col1a1",

  # endothelial
  "Pecam1",
  "Cdh5",
  "Ptprb",
  "Clec4g",

  # immune
  "Ptprc",
  "Lyz2",
  "C1qa",

  # cycling
  "Mki67",
  "Top2a",
  "Pcna",
  "Mcm5",
  "Uhrf1",
  "Dtl"
)

REVIEW_DIRECT_PRESENT <- present_genes(
  hep,
  REVIEW_DIRECT_GENES
)

review_cells <- rownames(
  hep@meta.data
)[
  hep$source_cluster_v650 %in%
    c(
      QC_REVIEW,
      LINEAGE_REVIEW
    )
]

counts_mat <- GetAssayData(
  hep,
  assay = "RNA",
  layer = "counts"
)

data_mat <- GetAssayData(
  hep,
  assay = "RNA",
  layer = "data"
)

review_direct_list <- list()

for (
  cl in c(
    QC_REVIEW,
    LINEAGE_REVIEW
  )
) {

  cells <- rownames(
    hep@meta.data
  )[
    hep$source_cluster_v650 ==
      cl
  ]

  if (
    !length(
      cells
    )
  ) {
    next
  }

  pcts <- Matrix::rowMeans(
    counts_mat[
      REVIEW_DIRECT_PRESENT,
      cells,
      drop = FALSE
    ] >
      0
  )

  means <- Matrix::rowMeans(
    data_mat[
      REVIEW_DIRECT_PRESENT,
      cells,
      drop = FALSE
    ]
  )

  review_direct_list[[
    cl
  ]] <- tibble(
    source_cluster =
      cl,
    gene =
      REVIEW_DIRECT_PRESENT,
    pct_expressed =
      as.numeric(
        pcts
      ),
    mean_logexpr =
      as.numeric(
        means
      )
  )
}

review_direct <- bind_rows(
  review_direct_list
)

write.csv(
  review_direct,
  file.path(
    TAB_OUT,
    "08_review_clusters_5_10_12_direct_marker_expression_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Create high-confidence-clean candidate object
# ==============================================================================

clean_cells <- rownames(
  hep@meta.data
)[
  hep$source_cluster_v650 %in%
    CLEAN_CANDIDATE_CLUSTERS
]

if (
  !length(
    clean_cells
  )
) {
  stop(
    "No cells remain after high-confidence cleanup."
  )
}

msg(
  "Cells before cleanup: ",
  ncol(
    hep
  )
)

msg(
  "Cells after high-confidence cleanup: ",
  length(
    clean_cells
  )
)

clean <- subset(
  hep,
  cells =
    clean_cells
)

DefaultAssay(
  clean
) <- "RNA"

clean <- safe_join_rna(
  clean
)

clean_counts <- clean@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    .data[[
      SAMPLE_COL
    ]],
    .data[[
      CONDITION_COL
    ]],
    name =
      "n_cells"
  ) %>%
  arrange(
    factor(
      .data[[
        SAMPLE_COL
      ]],
      levels =
        SAMPLES
    )
  )

write.csv(
  clean_counts,
  file.path(
    TAB_OUT,
    "09_clean_candidate_cell_counts_by_sample_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Re-run sample-aware RPCA after cleanup
# ==============================================================================

msg(
  "Splitting clean RNA assay by sample..."
)

clean[["RNA"]] <- split(
  clean[["RNA"]],
  f =
    clean@meta.data[[
      SAMPLE_COL
    ]]
)

msg(
  "NormalizeData..."
)

clean <- NormalizeData(
  clean,
  assay = "RNA",
  normalization.method =
    "LogNormalize",
  scale.factor =
    10000,
  verbose = FALSE
)

msg(
  "FindVariableFeatures..."
)

clean <- FindVariableFeatures(
  clean,
  assay = "RNA",
  selection.method =
    "vst",
  nfeatures =
    NFEATURES,
  verbose = FALSE
)

msg(
  "ScaleData..."
)

clean <- ScaleData(
  clean,
  assay = "RNA",
  features =
    VariableFeatures(
      clean
    ),
  verbose = FALSE
)

msg(
  "RunPCA..."
)

clean <- RunPCA(
  clean,
  assay = "RNA",
  features =
    VariableFeatures(
      clean
    ),
  npcs =
    NPCS,
  reduction.name =
    "pca.hepclean",
  reduction.key =
    "HEPCLEANPCA_",
  verbose = FALSE
)

msg(
  "IntegrateLayers RPCA..."
)

clean <- IntegrateLayers(
  object =
    clean,
  method =
    RPCAIntegration,
  orig.reduction =
    "pca.hepclean",
  new.reduction =
    "integrated.hepclean.rpca",
  normalization.method =
    "LogNormalize",
  dims =
    DIMS_USE,
  verbose = FALSE
)

clean <- FindNeighbors(
  clean,
  reduction =
    "integrated.hepclean.rpca",
  dims =
    DIMS_USE,
  graph.name = c(
    "hepclean_rpca_nn",
    "hepclean_rpca_snn"
  ),
  verbose = FALSE
)

for (
  res in RESOLUTIONS
) {

  nm <- paste0(
    "hepclean_rpca_res_",
    res
  )

  msg(
    "FindClusters clean resolution=",
    res
  )

  clean <- FindClusters(
    clean,
    graph.name =
      "hepclean_rpca_snn",
    resolution =
      res,
    cluster.name =
      nm,
    algorithm =
      1,
    random.seed =
      6510,
    verbose = FALSE
  )
}

msg(
  "RunUMAP clean..."
)

clean <- RunUMAP(
  clean,
  reduction =
    "integrated.hepclean.rpca",
  dims =
    DIMS_USE,
  reduction.name =
    CLEAN_UMAP,
  reduction.key =
    CLEAN_UMAP_KEY,
  seed.use =
    6510,
  verbose = FALSE
)

clean <- safe_join_rna(
  clean
)

DefaultAssay(
  clean
) <- "RNA"


# ==============================================================================
# 13. Re-score biological programs on clean object
# ==============================================================================

bio_gene_audit <- list()

for (
  nm in names(
    BIO_PROGRAMS
  )
) {

  ans <- add_program_score(
    clean,
    BIO_PROGRAMS[[
      nm
    ]],
    paste0(
      "v651_",
      nm
    )
  )

  clean <- ans$object

  bio_gene_audit[[
    nm
  ]] <- tibble(
    program =
      nm,
    gene =
      BIO_PROGRAMS[[
        nm
      ]],
    present =
      BIO_PROGRAMS[[
        nm
      ]] %in%
        ans$genes
  )
}

bio_gene_audit <- bind_rows(
  bio_gene_audit
)

write.csv(
  bio_gene_audit,
  file.path(
    TAB_OUT,
    "10_clean_biological_program_gene_audit_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Clean working-cluster markers
# ==============================================================================

if (
  !CLEAN_CLUSTER_COL %in%
    colnames(
      clean@meta.data
    )
) {
  stop(
    "Clean working cluster column missing: ",
    CLEAN_CLUSTER_COL
  )
}

Idents(
  clean
) <- clean@meta.data[[
  CLEAN_CLUSTER_COL
]]

msg(
  "FindAllMarkers clean res0.4..."
)

clean_markers <- FindAllMarkers(
  clean,
  assay = "RNA",
  only.pos =
    TRUE,
  test.use =
    "wilcox",
  min.pct =
    0.10,
  logfc.threshold =
    0.25,
  return.thresh =
    1,
  verbose = FALSE
)

write.csv(
  clean_markers,
  file.path(
    TAB_OUT,
    "11_clean_Hepatocyte_cluster_markers_res0.4_v6.5.1.csv"
  ),
  row.names = FALSE
)

clean_top20 <- clean_markers %>%
  group_by(
    cluster
  ) %>%
  arrange(
    desc(
      avg_log2FC
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write.csv(
  clean_top20,
  file.path(
    TAB_OUT,
    "12_clean_Hepatocyte_cluster_top20_markers_res0.4_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Clean cluster fractions
# ==============================================================================

clean_fraction <- cluster_fraction_table(
  clean,
  CLEAN_CLUSTER_COL,
  SAMPLE_COL,
  SAMPLES
)

write.csv(
  clean_fraction,
  file.path(
    TAB_OUT,
    "13_clean_Hepatocyte_cluster_fraction_by_sample_res0.4_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Clean cluster biological program summary
# ==============================================================================

bio_score_cols <- paste0(
  "v651_",
  names(
    BIO_PROGRAMS
  ),
  "_score"
)

clean_program_summary <- clean@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    clean_cluster =
      as.character(
        .data[[
          CLEAN_CLUSTER_COL
        ]]
      )
  ) %>%
  group_by(
    clean_cluster
  ) %>%
  summarise(
    n_cells =
      n(),
    across(
      all_of(
        bio_score_cols
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write.csv(
  clean_program_summary,
  file.path(
    TAB_OUT,
    "14_clean_cluster_program_summary_res0.4_v6.5.1.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Source-cluster 5 QC figures
# ==============================================================================

qc_plot_df <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    source_cluster =
      factor(
        source_cluster_v650,
        levels =
          as.character(
            sort(
              as.numeric(
                unique(
                  source_cluster_v650
                )
              )
            )
          )
      )
  )

p_qc_mt <- ggplot(
  qc_plot_df,
  aes(
    x =
      source_cluster,
    y =
      percent.mt
  )
) +
  geom_violin(
    scale =
      "width",
    trim =
      TRUE
  ) +
  geom_boxplot(
    width =
      0.12,
    outlier.shape =
      NA
  ) +
  geom_hline(
    yintercept =
      10,
    linetype =
      2,
    linewidth =
      0.3
  ) +
  geom_hline(
    yintercept =
      20,
    linetype =
      2,
    linewidth =
      0.3
  ) +
  labs(
    title =
      "Source-cluster mitochondrial QC audit",
    subtitle =
      "Cluster 5 is retained; thresholds are descriptive only",
    x =
      "v6.5.0 source cluster",
    y =
      "percent.mt"
  ) +
  theme_classic(
    base_size =
      9
  )

save_pdf(
  p_qc_mt,
  file.path(
    FIG_OUT,
    "01_source_cluster_percent_mt_QC_v6.5.1.pdf"
  ),
  10,
  6
)

p_qc_5 <- qc_plot_df %>%
  filter(
    source_cluster ==
      "5"
  ) %>%
  mutate(
    sample =
      factor(
        .data[[
          SAMPLE_COL
        ]],
        levels =
          SAMPLES
      )
  ) %>%
  ggplot(
    aes(
      x =
        sample,
      y =
        percent.mt,
      fill =
        .data[[
          CONDITION_COL
        ]]
    )
  ) +
  geom_violin(
    scale =
      "width",
    trim =
      TRUE
  ) +
  geom_boxplot(
    width =
      0.12,
    outlier.shape =
      NA
  ) +
  scale_fill_manual(
    values = c(
      "Sham" = "#0072B2",
      "Tx" = "#D55E00"
    )
  ) +
  labs(
    title =
      "Source cluster 5 | mitochondrial-high QC",
    x = NULL,
    y =
      "percent.mt",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size =
      9
  )

save_pdf(
  p_qc_5,
  file.path(
    FIG_OUT,
    "02_source_cluster5_percent_mt_by_sample_v6.5.1.pdf"
  ),
  8,
  5
)


# ==============================================================================
# 18. Review-cluster direct marker DotPlot
# ==============================================================================

review_obj <- subset(
  hep,
  cells =
    review_cells
)

if (
  length(
    REVIEW_DIRECT_PRESENT
  )
) {

  p_review_dot <- DotPlot(
    review_obj,
    features =
      REVIEW_DIRECT_PRESENT,
    group.by =
      "source_cluster_v650",
    assay =
      "RNA",
    dot.scale =
      8
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low =
        "#0033FF",
      mid =
        "#FFFFFF",
      high =
        "#FF1A1A",
      midpoint =
        0
    ) +
    ggtitle(
      "Source clusters 5 / 10 / 12 lineage audit"
    )

  save_pdf(
    p_review_dot,
    file.path(
      FIG_OUT,
      "03_review_clusters5_10_12_lineage_DotPlot_v6.5.1.pdf"
    ),
    14,
    5
  )
}


# ==============================================================================
# 19. Clean resolution comparison
# ==============================================================================

res_plots <- lapply(
  RESOLUTIONS,
  function(res) {

    nm <- paste0(
      "hepclean_rpca_res_",
      res
    )

    DimPlot(
      clean,
      reduction =
        CLEAN_UMAP,
      group.by =
        nm,
      label =
        TRUE,
      repel =
        TRUE,
      pt.size =
        0.22,
      raster =
        FALSE
    ) +
      ggtitle(
        paste0(
          "Clean res ",
          res
        )
      )
  }
)

p_res <- wrap_plots(
  res_plots,
  ncol =
    length(
      RESOLUTIONS
    )
) +
  plot_annotation(
    title =
      "Clean Hepatocyte RPCA resolution comparison | v6.5.1"
  )

save_pdf(
  p_res,
  file.path(
    FIG_OUT,
    "04_clean_Hepatocyte_resolution_comparison_v6.5.1.pdf"
  ),
  16,
  6
)


# ==============================================================================
# 20. Clean UMAP
# ==============================================================================

p_clean <- DimPlot(
  clean,
  reduction =
    CLEAN_UMAP,
  group.by =
    CLEAN_CLUSTER_COL,
  label =
    TRUE,
  repel =
    TRUE,
  pt.size =
    0.25,
  raster =
    FALSE
) +
  ggtitle(
    "High-confidence-clean Hepatocyte UMAP | res 0.4"
  )

save_pdf(
  p_clean,
  file.path(
    FIG_OUT,
    "05_clean_Hepatocyte_working_UMAP_res0.4_v6.5.1.pdf"
  ),
  8,
  7
)

p_clean_sample <- DimPlot(
  clean,
  reduction =
    CLEAN_UMAP,
  group.by =
    SAMPLE_COL,
  split.by =
    SAMPLE_COL,
  ncol =
    2,
  pt.size =
    0.25,
  raster =
    FALSE
) +
  plot_annotation(
    title =
      "Clean Hepatocyte UMAP by biological sample"
  )

save_pdf(
  p_clean_sample,
  file.path(
    FIG_OUT,
    "06_clean_Hepatocyte_UMAP_by_sample_v6.5.1.pdf"
  ),
  13,
  10
)


# ==============================================================================
# 21. Clean source-cluster provenance UMAP
# ==============================================================================

p_provenance <- DimPlot(
  clean,
  reduction =
    CLEAN_UMAP,
  group.by =
    "source_cluster_v650",
  label =
    TRUE,
  repel =
    TRUE,
  pt.size =
    0.25,
  raster =
    FALSE
) +
  ggtitle(
    "Clean UMAP colored by v6.5.0 source cluster"
  )

save_pdf(
  p_provenance,
  file.path(
    FIG_OUT,
    "07_clean_UMAP_source_cluster_provenance_v6.5.1.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 22. Clean canonical / state DotPlot
# ==============================================================================

CLEAN_DOT_GENES <- unique(
  c(
    "Alb",
    "Ttr",
    "Hnf4a",
    "Cps1",
    "Pck1",
    "G6pc",
    "Glul",
    "Cyp2e1",
    "Oat",

    "Areg",
    "Ccl2",
    "Cxcl10",
    "Rsad2",
    "Nupr1",
    "Sprr1a",

    "Ddit3",
    "Atf3",
    "Hspa5",
    "Hmox1",
    "Nqo1",

    "Mki67",
    "Top2a",
    "Pcna",
    "Mcm5",

    "Krt19",
    "Epcam",
    "Sox9",
    "Muc1",
    "Krt20",

    "Lrat",
    "Pdgfra",
    "Ptprc",
    "Pecam1"
  )
)

CLEAN_DOT_PRESENT <- present_genes(
  clean,
  CLEAN_DOT_GENES
)

if (
  length(
    CLEAN_DOT_PRESENT
  )
) {

  p_clean_dot <- DotPlot(
    clean,
    features =
      CLEAN_DOT_PRESENT,
    group.by =
      CLEAN_CLUSTER_COL,
    assay =
      "RNA",
    dot.scale =
      7
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low =
        "#0033FF",
      mid =
        "#FFFFFF",
      high =
        "#FF1A1A",
      midpoint =
        0
    ) +
    ggtitle(
      "Clean Hepatocyte lineage / zonation / injury audit | res 0.4"
    )

  save_pdf(
    p_clean_dot,
    file.path(
      FIG_OUT,
      "08_clean_Hepatocyte_lineage_state_DotPlot_res0.4_v6.5.1.pdf"
    ),
    15,
    7
  )
}


# ==============================================================================
# 23. Clean program score FeaturePlots
# ==============================================================================

bio_score_present <- intersect(
  bio_score_cols,
  colnames(
    clean@meta.data
  )
)

if (
  length(
    bio_score_present
  )
) {

  p_program <- FeaturePlot(
    clean,
    features =
      bio_score_present,
    reduction =
      CLEAN_UMAP,
    ncol =
      3,
    order =
      TRUE,
    pt.size =
      0.25,
    raster =
      FALSE
  ) &
    scale_color_gradient2(
      low =
        "#0033FF",
      mid =
        "#FFFFFF",
      high =
        "#FF1A1A",
      midpoint =
        0
    )

  save_pdf(
    p_program,
    file.path(
      FIG_OUT,
      "09_clean_Hepatocyte_program_FeaturePlots_v6.5.1.pdf"
    ),
    13,
    9
  )
}


# ==============================================================================
# 24. Clean cluster fraction heatmap
# ==============================================================================

clean_fraction_plot <- clean_fraction %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    cluster =
      factor(
        cluster,
        levels =
          sort(
            unique(
              cluster
            )
          )
      )
  )

p_frac <- ggplot(
  clean_fraction_plot,
  aes(
    x =
      sample,
    y =
      cluster,
    fill =
      fraction
  )
) +
  geom_tile(
    linewidth =
      0.3
  ) +
  scale_fill_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) +
  labs(
    title =
      "Clean Hepatocyte cluster fractions by sample | res 0.4",
    x = NULL,
    y =
      "Clean cluster",
    fill =
      "Fraction"
  ) +
  theme_classic(
    base_size =
      9
  )

save_pdf(
  p_frac,
  file.path(
    FIG_OUT,
    "10_clean_Hepatocyte_cluster_fraction_heatmap_v6.5.1.pdf"
  ),
  8,
  6
)


# ==============================================================================
# 25. Save RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_Hepatocyte_highconf_clean_reclustered_v6.5.1.rds"
)

saveRDS(
  clean,
  RDS_FILE,
  compress =
    FALSE
)

msg(
  "Saved clean Hepatocyte RDS: ",
  RDS_FILE
)


# ==============================================================================
# 26. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "source_cluster_column",
    "high_conf_remove_clusters",
    "core_retain_clusters",
    "QC_review_clusters",
    "lineage_review_clusters",
    "percent_mt_filter_applied",
    "cycling_review_cells_removed",
    "recluster_method",
    "resolutions",
    "working_resolution",
    "clean_cluster_column",
    "final_annotation_fixed"
  ),

  value = c(
    "v6.5.1",
    INPUT_RDS,
    SOURCE_CLUSTER_COL,
    paste(
      HIGH_CONF_REMOVE,
      collapse = ","
    ),
    paste(
      CORE_RETAIN,
      collapse = ","
    ),
    paste(
      QC_REVIEW,
      collapse = ","
    ),
    paste(
      LINEAGE_REVIEW,
      collapse = ","
    ),
    "FALSE",
    "FALSE",
    "Seurat v5 RPCAIntegration after high-confidence contaminant removal",
    paste(
      RESOLUTIONS,
      collapse = ","
    ),
    as.character(
      WORKING_RESOLUTION
    ),
    CLEAN_CLUSTER_COL,
    "FALSE"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.5.1.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.5.1.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
