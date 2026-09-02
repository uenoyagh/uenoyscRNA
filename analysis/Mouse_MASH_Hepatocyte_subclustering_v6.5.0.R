#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6500)

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
# Hepatocyte subclustering / zonation / injury-state audit
#
# Version: v6.5.0
#
# PURPOSE
#   Prepare a hepatocyte-specific map before interpreting transplantation-
#   associated hepatocyte changes.
#
# PRIMARY QUESTIONS
#   1) Are the apparent whole-cell Hepatocyte increases in Tx explained by
#      zonation/composition shifts?
#   2) Within comparable hepatocyte states, are injury/stress programs lower
#      in Tx?
#   3) Are there discrete stress-high / regenerative hepatocyte states?
#
# IMPORTANT
#   - Final hepatocyte annotation is NOT fixed in v6.5.0.
#   - This is an audit / subclustering version.
#   - Whole-cell parent annotation is preserved.
#   - No attempt is made here to recover cells labeled "Cycling" at Layer1.
#     Cycling-lineage origin will be audited separately later.
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

canonical_sample <- function(x) {
  x <- as.character(
    x
  )

  dplyr::case_when(
    grepl(
      "^Sham1$",
      x,
      ignore.case = TRUE
    ) ~
      "Sham1",

    grepl(
      "^Sham20$",
      x,
      ignore.case = TRUE
    ) ~
      "Sham20",

    grepl(
      "^Tx17$",
      x,
      ignore.case = TRUE
    ) ~
      "Tx17",

    grepl(
      "^Tx5$",
      x,
      ignore.case = TRUE
    ) ~
      "Tx5",

    TRUE ~
      x
  )
}

canonical_condition <- function(sample) {
  ifelse(
    grepl(
      "^Tx",
      sample
    ),
    "Tx",
    "Sham"
  )
}

detect_sample_column <- function(
  object,
  required_samples
) {

  candidates <- c(
    "sample_for_annotation",
    "sample_interaction_v620",
    "sample",
    "Sample",
    "orig.ident",
    "sample_id",
    "SampleID",
    "sample_name"
  )

  candidates <- intersect(
    candidates,
    colnames(
      object@meta.data
    )
  )

  if (
    !length(
      candidates
    )
  ) {
    stop(
      "No candidate sample metadata column found."
    )
  }

  for (
    nm in candidates
  ) {
    vals <- canonical_sample(
      object@meta.data[[
        nm
      ]]
    )

    if (
      all(
        required_samples %in%
          unique(
            vals
          )
      )
    ) {
      return(
        nm
      )
    }
  }

  # Wider search across character/factor metadata if canonical candidate names fail.
  for (
    nm in colnames(
      object@meta.data
    )
  ) {

    x <- object@meta.data[[
      nm
    ]]

    if (
      !is.character(
        x
      ) &&
      !is.factor(
        x
      )
    ) {
      next
    }

    vals <- canonical_sample(
      x
    )

    if (
      all(
        required_samples %in%
          unique(
            vals
          )
      )
    ) {
      return(
        nm
      )
    }
  }

  stop(
    "Could not identify a metadata column containing Sham1, Sham20, Tx17, Tx5."
  )
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
      paste0(
        name,
        "_score"
      )
    ]] <- NA_real_

    return(
      list(
        object = object,
        genes = genes_use
      )
    )
  }

  object <- AddModuleScore(
    object,
    features = list(
      genes_use
    ),
    name = paste0(
      name,
      "_"
    ),
    assay = "RNA",
    seed = 6500
  )

  generated <- paste0(
    name,
    "_1"
  )

  final_name <- paste0(
    name,
    "_score"
  )

  object[[
    final_name
  ]] <- object[[
    generated
  ]][
    ,
    1
  ]

  object[[
    generated
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
    ungroup() %>%
    mutate(
      condition =
        canonical_condition(
          sample
        )
    )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_subclustering_v6.5.0"
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

LAYER1_COL <-
  "wholecell_layer1_FINAL_v511"

TARGET_LABEL <-
  "Hepatocyte"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
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
  0.6,
  0.8
)

WORKING_RESOLUTION <-
  0.4

WORKING_CLUSTER_COL <-
  "hep_rpca_res_0.4"

UMAP_NAME <-
  "umap.hep.rpca"

UMAP_KEY <-
  "HEPUMAP_"


# ==============================================================================
# 4. A-priori hepatocyte programs
# ==============================================================================

PROGRAMS <- list(

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

  Inflammatory_stress = c(
    "Nfkbia",
    "Socs3",
    "Stat3",
    "Ccl2",
    "Cxcl1",
    "Icam1",
    "Saa1",
    "Saa2"
  ),

  Hypoxia_glycolytic = c(
    "Hif1a",
    "Vegfa",
    "Bnip3",
    "Slc2a1",
    "Ldha",
    "Pdk1"
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

CONTAMINATION_MARKERS <- list(

  Cholangiocyte = c(
    "Krt19",
    "Krt7",
    "Epcam",
    "Sox9",
    "Krt8",
    "Krt18"
  ),

  Endothelial = c(
    "Pecam1",
    "Cdh5",
    "Kdr",
    "Emcn",
    "Vwf"
  ),

  Macrophage_immune = c(
    "Ptprc",
    "Lyz2",
    "Adgre1",
    "C1qa",
    "C1qb",
    "C1qc",
    "Cd68"
  ),

  HSC_mesenchymal = c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Lrat",
    "Rbp1",
    "Dcn",
    "Lum"
  )
)


# ==============================================================================
# 5. Load frozen whole-cell parent
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Frozen whole-cell input not found: ",
    INPUT_RDS
  )
}

msg(
  "Loading frozen whole-cell parent..."
)

parent <- readRDS(
  INPUT_RDS
)

if (
  !"RNA" %in%
    Assays(
      parent
    )
) {
  stop(
    "RNA assay missing."
  )
}

if (
  !LAYER1_COL %in%
    colnames(
      parent@meta.data
    )
) {
  stop(
    "Layer1 column missing: ",
    LAYER1_COL
  )
}

SAMPLE_SOURCE_COL <- detect_sample_column(
  parent,
  SAMPLES
)

msg(
  "Detected sample metadata column: ",
  SAMPLE_SOURCE_COL
)

parent$sample_hep_v650 <- canonical_sample(
  parent@meta.data[[
    SAMPLE_SOURCE_COL
  ]]
)

parent$condition_hep_v650 <- canonical_condition(
  parent$sample_hep_v650
)

keep_cells <- rownames(
  parent@meta.data
)[
  as.character(
    parent@meta.data[[
      LAYER1_COL
    ]]
  ) ==
    TARGET_LABEL &
    parent$sample_hep_v650 %in%
      SAMPLES
]

if (
  !length(
    keep_cells
  )
) {
  stop(
    "No Hepatocyte cells found for Sham1/Sham20/Tx17/Tx5."
  )
}

msg(
  "Hepatocyte cells retained: ",
  length(
    keep_cells
  )
)

hep <- subset(
  parent,
  cells =
    keep_cells
)

rm(
  parent
)

gc()

DefaultAssay(
  hep
) <- "RNA"

hep <- safe_join_rna(
  hep
)


# ==============================================================================
# 6. Input audit
# ==============================================================================

input_counts <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  count(
    sample_hep_v650,
    condition_hep_v650,
    name =
      "n_cells"
  ) %>%
  arrange(
    factor(
      sample_hep_v650,
      levels =
        SAMPLES
    )
  )

write.csv(
  input_counts,
  file.path(
    TAB_OUT,
    "01_Hepatocyte_input_cell_counts_v6.5.0.csv"
  ),
  row.names = FALSE
)

msg(
  "Input counts: ",
  paste0(
    input_counts$sample_hep_v650,
    "=",
    input_counts$n_cells,
    collapse = "; "
  )
)


# ==============================================================================
# 7. Split RNA layers by biological sample and preprocess
# ==============================================================================

msg(
  "Splitting RNA assay by biological sample..."
)

hep[["RNA"]] <- split(
  hep[["RNA"]],
  f =
    hep$sample_hep_v650
)

msg(
  "NormalizeData..."
)

hep <- NormalizeData(
  hep,
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

hep <- FindVariableFeatures(
  hep,
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

hep <- ScaleData(
  hep,
  assay = "RNA",
  features =
    VariableFeatures(
      hep
    ),
  verbose = FALSE
)

msg(
  "RunPCA..."
)

hep <- RunPCA(
  hep,
  assay = "RNA",
  features =
    VariableFeatures(
      hep
    ),
  npcs =
    NPCS,
  reduction.name =
    "pca",
  verbose = FALSE
)


# ==============================================================================
# 8. RPCA integration
# ==============================================================================

msg(
  "IntegrateLayers using RPCAIntegration..."
)

hep <- IntegrateLayers(
  object =
    hep,
  method =
    RPCAIntegration,
  orig.reduction =
    "pca",
  new.reduction =
    "integrated.rpca",
  normalization.method =
    "LogNormalize",
  dims =
    DIMS_USE,
  verbose = FALSE
)

msg(
  "FindNeighbors..."
)

hep <- FindNeighbors(
  hep,
  reduction =
    "integrated.rpca",
  dims =
    DIMS_USE,
  graph.name = c(
    "hep_rpca_nn",
    "hep_rpca_snn"
  ),
  verbose = FALSE
)


# ==============================================================================
# 9. Resolution comparison
# ==============================================================================

for (
  res in RESOLUTIONS
) {

  col_name <- paste0(
    "hep_rpca_res_",
    res
  )

  msg(
    "FindClusters resolution=",
    res
  )

  hep <- FindClusters(
    hep,
    graph.name =
      "hep_rpca_snn",
    resolution =
      res,
    cluster.name =
      col_name,
    algorithm =
      1,
    random.seed =
      6500,
    verbose = FALSE
  )
}

msg(
  "RunUMAP..."
)

hep <- RunUMAP(
  hep,
  reduction =
    "integrated.rpca",
  dims =
    DIMS_USE,
  reduction.name =
    UMAP_NAME,
  reduction.key =
    UMAP_KEY,
  seed.use =
    6500,
  verbose = FALSE
)


# ==============================================================================
# 10. Join RNA layers again for scoring / markers
# ==============================================================================

hep <- safe_join_rna(
  hep
)

DefaultAssay(
  hep
) <- "RNA"


# ==============================================================================
# 11. Add hepatocyte program scores
# ==============================================================================

msg(
  "Adding hepatocyte program scores..."
)

program_gene_audit <- list()

for (
  nm in names(
    PROGRAMS
  )
) {

  ans <- add_program_score(
    hep,
    PROGRAMS[[
      nm
    ]],
    nm
  )

  hep <- ans$object

  program_gene_audit[[
    nm
  ]] <- tibble(
    program =
      nm,
    gene =
      PROGRAMS[[
        nm
      ]],
    present =
      PROGRAMS[[
        nm
      ]] %in%
        ans$genes
  )
}

program_gene_audit <- bind_rows(
  program_gene_audit
)

write.csv(
  program_gene_audit,
  file.path(
    TAB_OUT,
    "02_Hepatocyte_program_gene_audit_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Cluster marker discovery at working resolution
# ==============================================================================

if (
  !WORKING_CLUSTER_COL %in%
    colnames(
      hep@meta.data
    )
) {
  stop(
    "Working cluster column missing: ",
    WORKING_CLUSTER_COL
  )
}

Idents(
  hep
) <- hep@meta.data[[
  WORKING_CLUSTER_COL
]]

msg(
  "FindAllMarkers at working resolution ",
  WORKING_RESOLUTION,
  "..."
)

markers <- FindAllMarkers(
  hep,
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
  markers,
  file.path(
    TAB_OUT,
    "03_Hepatocyte_cluster_markers_res0.4_v6.5.0.csv"
  ),
  row.names = FALSE
)

top20 <- markers %>%
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
  top20,
  file.path(
    TAB_OUT,
    "04_Hepatocyte_cluster_top20_markers_res0.4_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Cluster fractions by sample
# ==============================================================================

cluster_fraction <- cluster_fraction_table(
  hep,
  WORKING_CLUSTER_COL,
  "sample_hep_v650",
  SAMPLES
)

write.csv(
  cluster_fraction,
  file.path(
    TAB_OUT,
    "05_Hepatocyte_cluster_fraction_by_sample_res0.4_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Cluster-level program summary
# ==============================================================================

score_cols <- paste0(
  names(
    PROGRAMS
  ),
  "_score"
)

cluster_program_summary <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  mutate(
    working_cluster =
      as.character(
        .data[[
          WORKING_CLUSTER_COL
        ]]
      )
  ) %>%
  group_by(
    working_cluster
  ) %>%
  summarise(
    n_cells =
      n(),
    across(
      all_of(
        score_cols
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write.csv(
  cluster_program_summary,
  file.path(
    TAB_OUT,
    "06_Hepatocyte_cluster_program_summary_res0.4_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Sample-level program summary
# ==============================================================================

sample_program_summary <- hep@meta.data %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  group_by(
    sample_hep_v650,
    condition_hep_v650
  ) %>%
  summarise(
    n_cells =
      n(),
    across(
      all_of(
        score_cols
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  arrange(
    factor(
      sample_hep_v650,
      levels =
        SAMPLES
    )
  )

write.csv(
  sample_program_summary,
  file.path(
    TAB_OUT,
    "07_Hepatocyte_program_scores_by_sample_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Contamination gene audit
# ==============================================================================

contamination_gene_table <- bind_rows(
  lapply(
    names(
      CONTAMINATION_MARKERS
    ),
    function(nm) {
      tibble(
        panel =
          nm,
        gene =
          CONTAMINATION_MARKERS[[
            nm
          ]],
        present =
          CONTAMINATION_MARKERS[[
            nm
          ]] %in%
            rownames(
              hep
            )
      )
    }
  )
)

write.csv(
  contamination_gene_table,
  file.path(
    TAB_OUT,
    "08_Hepatocyte_contamination_marker_audit_v6.5.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Figures: PCA elbow
# ==============================================================================

p_elbow <- ElbowPlot(
  hep,
  ndims =
    NPCS,
  reduction =
    "pca"
) +
  ggtitle(
    "Hepatocyte PCA elbow | v6.5.0"
  )

save_pdf(
  p_elbow,
  file.path(
    FIG_OUT,
    "01_Hepatocyte_PCA_elbow_v6.5.0.pdf"
  ),
  7,
  5
)


# ==============================================================================
# 18. Figures: resolution comparison
# ==============================================================================

res_plots <- lapply(
  RESOLUTIONS,
  function(res) {

    col_name <- paste0(
      "hep_rpca_res_",
      res
    )

    DimPlot(
      hep,
      reduction =
        UMAP_NAME,
      group.by =
        col_name,
      label =
        TRUE,
      repel =
        TRUE,
      pt.size =
        0.20,
      raster =
        FALSE
    ) +
      ggtitle(
        paste0(
          "Resolution ",
          res
        )
      )
  }
)

p_res <- wrap_plots(
  res_plots,
  ncol = 2
) +
  plot_annotation(
    title =
      "Hepatocyte RPCA clustering resolution comparison | v6.5.0"
  )

save_pdf(
  p_res,
  file.path(
    FIG_OUT,
    "02_Hepatocyte_resolution_comparison_v6.5.0.pdf"
  ),
  13,
  11
)


# ==============================================================================
# 19. Figures: working UMAP
# ==============================================================================

p_working <- DimPlot(
  hep,
  reduction =
    UMAP_NAME,
  group.by =
    WORKING_CLUSTER_COL,
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
    "Hepatocyte working clustering | RPCA res 0.4"
  )

save_pdf(
  p_working,
  file.path(
    FIG_OUT,
    "03_Hepatocyte_working_UMAP_res0.4_v6.5.0.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 20. Figures: sample and condition split
# ==============================================================================

p_sample <- DimPlot(
  hep,
  reduction =
    UMAP_NAME,
  group.by =
    "sample_hep_v650",
  split.by =
    "sample_hep_v650",
  ncol =
    2,
  pt.size =
    0.25,
  raster =
    FALSE
) +
  plot_annotation(
    title =
      "Hepatocyte UMAP by biological sample"
  )

save_pdf(
  p_sample,
  file.path(
    FIG_OUT,
    "04_Hepatocyte_UMAP_by_sample_v6.5.0.pdf"
  ),
  13,
  10
)

p_condition <- DimPlot(
  hep,
  reduction =
    UMAP_NAME,
  group.by =
    "condition_hep_v650",
  split.by =
    "condition_hep_v650",
  ncol =
    2,
  pt.size =
    0.25,
  raster =
    FALSE
) +
  plot_annotation(
    title =
      "Hepatocyte UMAP | Sham vs Tx"
  )

save_pdf(
  p_condition,
  file.path(
    FIG_OUT,
    "05_Hepatocyte_UMAP_Sham_vs_Tx_v6.5.0.pdf"
  ),
  13,
  6
)


# ==============================================================================
# 21. Figures: canonical marker FeaturePlots
# ==============================================================================

FEATURE_GENES <- unique(
  c(
    "Cps1",
    "Ass1",
    "Pck1",
    "G6pc",
    "Glul",
    "Cyp2e1",
    "Cyp1a2",
    "Oat",
    "Ddit3",
    "Atf3",
    "Hspa5",
    "Hmox1",
    "Nqo1",
    "Socs3",
    "Ccl2",
    "Mki67",
    "Top2a"
  )
)

FEATURE_PRESENT <- present_genes(
  hep,
  FEATURE_GENES
)

if (
  length(
    FEATURE_PRESENT
  )
) {

  p_features <- FeaturePlot(
    hep,
    features =
      FEATURE_PRESENT,
    reduction =
      UMAP_NAME,
    ncol =
      4,
    order =
      TRUE,
    min.cutoff =
      "q05",
    max.cutoff =
      "q95",
    pt.size =
      0.25,
    raster =
      FALSE
  ) &
    scale_color_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    )

  save_pdf(
    p_features,
    file.path(
      FIG_OUT,
      "06_Hepatocyte_canonical_marker_FeaturePlots_v6.5.0.pdf"
    ),
    14,
    max(
      8,
      3.1 *
        ceiling(
          length(
            FEATURE_PRESENT
          ) /
            4
        )
    )
  )
}


# ==============================================================================
# 22. Figures: program-score FeaturePlots
# ==============================================================================

score_present <- intersect(
  score_cols,
  colnames(
    hep@meta.data
  )
)

if (
  length(
    score_present
  )
) {

  p_scores <- FeaturePlot(
    hep,
    features =
      score_present,
    reduction =
      UMAP_NAME,
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
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    )

  save_pdf(
    p_scores,
    file.path(
      FIG_OUT,
      "07_Hepatocyte_program_score_FeaturePlots_v6.5.0.pdf"
    ),
    13,
    max(
      8,
      3.3 *
        ceiling(
          length(
            score_present
          ) /
            3
        )
    )
  )
}


# ==============================================================================
# 23. Figures: canonical marker DotPlot
# ==============================================================================

DOT_GENES <- unique(
  c(
    PROGRAMS$Periportal_Z1,
    PROGRAMS$Pericentral_Z3,
    PROGRAMS$ER_UPR_stress,
    PROGRAMS$Oxidative_stress,
    PROGRAMS$Inflammatory_stress,
    PROGRAMS$Regenerative_cycling
  )
)

DOT_PRESENT <- present_genes(
  hep,
  DOT_GENES
)

if (
  length(
    DOT_PRESENT
  )
) {

  p_dot <- DotPlot(
    hep,
    features =
      DOT_PRESENT,
    group.by =
      WORKING_CLUSTER_COL,
    assay =
      "RNA",
    dot.scale =
      7
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    ggtitle(
      "Hepatocyte canonical programs by cluster | res 0.4"
    )

  save_pdf(
    p_dot,
    file.path(
      FIG_OUT,
      "08_Hepatocyte_canonical_program_DotPlot_res0.4_v6.5.0.pdf"
    ),
    16,
    8
  )
}


# ==============================================================================
# 24. Figures: program-score violin by cluster
# ==============================================================================

if (
  length(
    score_present
  )
) {

  violin_list <- lapply(
    score_present,
    function(sc) {

      VlnPlot(
        hep,
        features =
          sc,
        group.by =
          WORKING_CLUSTER_COL,
        pt.size =
          0,
        assay =
          "RNA"
      ) +
        ggtitle(
          sc
        ) +
        NoLegend()
    }
  )

  p_violin <- wrap_plots(
    violin_list,
    ncol = 2
  ) +
    plot_annotation(
      title =
        "Hepatocyte program scores by cluster | res 0.4"
    )

  save_pdf(
    p_violin,
    file.path(
      FIG_OUT,
      "09_Hepatocyte_program_score_violin_res0.4_v6.5.0.pdf"
    ),
    14,
    max(
      10,
      4 *
        ceiling(
          length(
            violin_list
          ) /
            2
        )
    )
  )
}


# ==============================================================================
# 25. Figures: contamination DotPlot
# ==============================================================================

CONTAM_GENES <- unique(
  unlist(
    CONTAMINATION_MARKERS,
    use.names = FALSE
  )
)

CONTAM_PRESENT <- present_genes(
  hep,
  CONTAM_GENES
)

if (
  length(
    CONTAM_PRESENT
  )
) {

  p_contam <- DotPlot(
    hep,
    features =
      CONTAM_PRESENT,
    group.by =
      WORKING_CLUSTER_COL,
    assay =
      "RNA",
    dot.scale =
      7
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    ggtitle(
      "Hepatocyte contamination audit | res 0.4"
    )

  save_pdf(
    p_contam,
    file.path(
      FIG_OUT,
      "10_Hepatocyte_contamination_DotPlot_res0.4_v6.5.0.pdf"
    ),
    14,
    7
  )
}


# ==============================================================================
# 26. Figures: cluster fraction heatmap
# ==============================================================================

cluster_heat <- cluster_fraction %>%
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

p_heat <- ggplot(
  cluster_heat,
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
    linewidth = 0.25
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
      "Hepatocyte cluster fractions by biological sample | res 0.4",
    x = NULL,
    y =
      "Cluster",
    fill =
      "Fraction"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p_heat,
  file.path(
    FIG_OUT,
    "11_Hepatocyte_cluster_fraction_heatmap_res0.4_v6.5.0.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 27. Figures: cluster fraction line plot
# ==============================================================================

p_frac_line <- ggplot(
  cluster_fraction %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    ),
  aes(
    x =
      sample,
    y =
      fraction,
    group =
      cluster
  )
) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.5
  ) +
  facet_wrap(
    ~ cluster,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    title =
      "Hepatocyte cluster fractions | Sham1, Sham20, Tx17, Tx5",
    x = NULL,
    y =
      "Fraction within hepatocytes"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_frac_line,
  file.path(
    FIG_OUT,
    "12_Hepatocyte_cluster_fraction_lineplot_res0.4_v6.5.0.pdf"
  ),
  14,
  10
)


# ==============================================================================
# 28. Figures: sample-level program scores
# ==============================================================================

sample_program_long <- sample_program_summary %>%
  pivot_longer(
    cols =
      all_of(
        score_present
      ),
    names_to =
      "program",
    values_to =
      "mean_score"
  ) %>%
  mutate(
    sample_hep_v650 =
      factor(
        sample_hep_v650,
        levels =
          SAMPLES
      )
  )

if (
  nrow(
    sample_program_long
  )
) {

  p_sample_program <- ggplot(
    sample_program_long,
    aes(
      x =
        sample_hep_v650,
      y =
        mean_score,
      fill =
        condition_hep_v650
    )
  ) +
    geom_col() +
    facet_wrap(
      ~ program,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values = c(
        "Sham" = "#0072B2",
        "Tx" = "#D55E00"
      )
    ) +
    labs(
      title =
        "Hepatocyte program scores by biological sample",
      subtitle =
        "Exploratory; sample-level means",
      x = NULL,
      y =
        "Mean module score",
      fill =
        "Condition"
    ) +
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_sample_program,
    file.path(
      FIG_OUT,
      "13_Hepatocyte_program_scores_by_sample_v6.5.0.pdf"
    ),
    13,
    9
  )
}


# ==============================================================================
# 29. Save RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_Hepatocyte_subclustered_v6.5.0.rds"
)

saveRDS(
  hep,
  RDS_FILE,
  compress = FALSE
)

msg(
  "Saved Hepatocyte RDS: ",
  RDS_FILE
)


# ==============================================================================
# 30. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "parent_annotation_column",
    "target_parent_label",
    "sample_source_column",
    "samples",
    "integration",
    "nfeatures",
    "npcs",
    "dims",
    "resolutions",
    "working_resolution",
    "working_cluster_column",
    "UMAP_reduction",
    "final_annotation_fixed",
    "cycling_parent_cells_included"
  ),

  value = c(
    "v6.5.0",
    INPUT_RDS,
    LAYER1_COL,
    TARGET_LABEL,
    SAMPLE_SOURCE_COL,
    paste(
      SAMPLES,
      collapse = ","
    ),
    "Seurat v5 RPCAIntegration",
    as.character(
      NFEATURES
    ),
    as.character(
      NPCS
    ),
    paste(
      range(
        DIMS_USE
      ),
      collapse = ":"
    ),
    paste(
      RESOLUTIONS,
      collapse = ","
    ),
    as.character(
      WORKING_RESOLUTION
    ),
    WORKING_CLUSTER_COL,
    UMAP_NAME,
    "FALSE",
    "FALSE; Layer1 Cycling cells will be audited separately"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.5.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.5.0.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
