#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(5820)

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
# Whole-cell specificity audit for macrophage-subtype staining candidates
#
# Version: v5.8.6
#
# v5.8.6 CHANGE FROM v5.8.2
#   - Added Spn/CD43 back to the Inflammatory-Mphi candidate set.
#   - Added classic reference markers Cd68, Cd80, Mrc1 (CD206), Cd86.
#   - Added Clean-B Mphi reference UMAP panels (all Mphi and Sham-vs-Tx).
#   - Added these reference genes to whole-cell specificity outputs.
#   - Whole-cell specificity logic and frozen whole-cell reference remain unchanged.
#
# v5.8.6 CHANGE FROM v5.8.3
#   Error fix only:
#   - Whole-cell Layer1 macrophage/Kupffer label is resolved automatically.
#   - Monocyte label is also resolved automatically.
#   - No candidate genes, specificity formulas, UMAP coordinates, or ranking
#     definitions are changed.
#
# v5.8.6 CHANGE FROM v5.8.4
#   Error fix only:
#   - Internal helper metadata names are unified to suffix v586.
#   - Fixes the previous create/read mismatch in Layer1 helper metadata.
#   - Candidate genes, annotation source, UMAP coordinates, specificity formulas,
#     reference-marker panels, and ranking definitions are unchanged.
#
# v5.8.6 CHANGE FROM v5.8.5
#   - Adds Thbs1 and Tgfbi to ECM-associated inflammatory-Mphi candidates.
#   - Adds focused ECM-associated inflammatory-Mphi validation on Clean-B Mphi:
#       F13a1, Thbs1, Tgfbi, Vcan, Ccr2
#   - Outputs full-Mphi UMAP, Sham-vs-Tx UMAP, subtype DotPlot/Violin,
#     and sample-level Sham1/Sham20/Tx17/Tx5 quantitative metrics.
#   - Existing whole-cell specificity logic and frozen reference objects are
#     otherwise unchanged.
#
# PURPOSE
#   Evaluate candidate staining markers from the finalized Clean-B macrophage
#   subtype discovery in the FULL liver-cell population.
#
# CANDIDATE GROUPS
#   Anti-inflammatory:
#     Cd163, Clec4f, Slc40a1
#
#   Inflammatory:
#     Slc2a6, Smpdl3b
#
#   ECM-associated inflammatory:
#     F13a1, Vcan
#
#   Repair/Resolution:
#     Mmp13, Ksr2
#
#   Lipid-associated/TREM2:
#     Gpnmb, Atp6v0d2, Trem2
#
#   Reference / overlap controls:
#     Ccr2, Mmp12
#
# ANALYSES
#   1) Full-cell UMAP FeaturePlot
#   2) Layer1 DotPlot
#   3) % RNA-positive and mean expression by Layer1 cell type
#   4) Fraction of ALL gene-positive cells that are Kupffer_Macrophage
#   5) Fraction of ALL gene-positive cells that are Kupffer_Macrophage + Monocyte
#   6) Kupffer enrichment relative to all non-Kupffer cells
#   7) Top non-Kupffer cell type for each gene
#   8) Sham-vs-Tx full-cell UMAP
#   9) Candidate-by-celltype specificity heatmaps
#
# IMPORTANT
#   - Uses frozen whole-cell parent v5.1.1.
#   - No reclustering, reintegration, or re-UMAP.
#   - Existing umapRPCA and wholecell_layer1_FINAL_v511 are preserved.
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

save_pdf <- function(plot_obj, filename, width, height) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot_obj)
  grDevices::dev.off()
}

canonical_condition <- function(x) {
  x <- as.character(x)

  out <- dplyr::case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
    TRUE ~ NA_character_
  )

  factor(
    out,
    levels = c(
      "STD",
      "CDAHFD",
      "Sham",
      "Tx"
    )
  )
}

resolve_first <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

make_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

# ==============================================================================
# 2. Paths
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

RDS_PATH <- file.path(
  DATA_ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

MPHI_RDS_CANDIDATES <- c(
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
  ),
  file.path(
    DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
  )
)

MPHI_RDS_PATH <- MPHI_RDS_CANDIDATES[
  file.exists(MPHI_RDS_CANDIDATES)
][1]

OUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_MphiSubtype_StainingMarker_Specificity_v5.8.6"
)

FIG_DIR <- file.path(
  OUT_DIR,
  "Figures"
)

TAB_DIR <- file.path(
  OUT_DIR,
  "Tables"
)

LOG_DIR <- file.path(
  OUT_DIR,
  "Logs"
)

for (d in c(
  OUT_DIR,
  FIG_DIR,
  TAB_DIR,
  LOG_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 3. Settings
# ==============================================================================

ASSAY_USE <- "RNA"

ANNOTATION_COL <- "wholecell_layer1_FINAL_v511"

UMAP_CANDIDATES <- c(
  "umapRPCA",
  "umap"
)

CONDITION_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group",
  "sample_for_annotation",
  "sample"
)

TARGET_KUPFFER_LABEL <- "Kupffer_Macrophage"
TARGET_MONOCYTE_LABEL <- "Monocyte"

CANDIDATE_GROUPS <- list(

  "Anti-inflammatory" = c(
    "Cd163",
    "Clec4f",
    "Slc40a1"
  ),

  "Inflammatory" = c(
    "Spn",
    "Slc2a6",
    "Smpdl3b"
  ),

  "ECM-associated inflammatory" = c(
    "F13a1",
    "Thbs1",
    "Tgfbi",
    "Vcan"
  ),

  "Repair/Resolution" = c(
    "Mmp13",
    "Ksr2"
  ),

  "Lipid-associated/TREM2" = c(
    "Gpnmb",
    "Atp6v0d2",
    "Trem2"
  ),

  "Reference controls" = c(
    "Cd68",
    "Cd80",
    "Mrc1",
    "Cd86",
    "Ccr2",
    "Mmp12"
  )
)

# Classic / contextual reference panel for Clean-B macrophage UMAPs.
# CD206 is represented by the mouse gene Mrc1.
REFERENCE_MPHI_MARKERS <- c(
  "Cd68",
  "Cd80",
  "Mrc1",
  "Cd86",
  "Cd163",
  "Clec4f",
  "Slc40a1"
)

REFERENCE_MPHI_DISPLAY <- c(
  "Cd68" = "CD68",
  "Cd80" = "CD80",
  "Mrc1" = "CD206 (Mrc1)",
  "Cd86" = "CD86",
  "Cd163" = "CD163",
  "Clec4f" = "CLEC4F",
  "Slc40a1" = "SLC40A1"
)


# Focused ECM-associated inflammatory-Mphi validation panel.
ECM_FOCUS_MARKERS <- c(
  "F13a1",
  "Thbs1",
  "Tgfbi",
  "Vcan",
  "Ccr2"
)

ECM_FOCUS_DISPLAY <- c(
  "F13a1" = "F13A1",
  "Thbs1" = "THBS1",
  "Tgfbi" = "TGFBI",
  "Vcan" = "VCAN",
  "Ccr2" = "CCR2"
)

MPHI_SUBTYPE_COL <- "macrophage_class_Res2_FINAL_v4145_char"

MPHI_SUBTYPE_LEVELS <- c(
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Other"
)

MPHI_SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "orig.ident",
  "sample"
)

TARGET_GENES <- unique(
  unlist(
    CANDIDATE_GROUPS,
    use.names = FALSE
  )
)

# ==============================================================================
# 4. Load object
# ==============================================================================

if (!file.exists(RDS_PATH)) {
  stop(
    "Whole-cell frozen RDS not found: ",
    RDS_PATH
  )
}

msg(
  "Loading whole-cell frozen RDS..."
)

obj <- readRDS(
  RDS_PATH
)

if (!ASSAY_USE %in% Assays(obj)) {
  stop(
    "RNA assay missing."
  )
}

DefaultAssay(
  obj
) <- ASSAY_USE

if (!ANNOTATION_COL %in% colnames(obj@meta.data)) {
  stop(
    "Frozen Layer1 annotation column missing: ",
    ANNOTATION_COL
  )
}

umap_use <- resolve_first(
  Reductions(obj),
  UMAP_CANDIDATES
)

if (is.na(umap_use)) {
  stop(
    "No usable UMAP reduction found."
  )
}

condition_col <- resolve_first(
  colnames(obj@meta.data),
  CONDITION_CANDIDATES
)

if (is.na(condition_col)) {
  stop(
    "No usable condition column found."
  )
}

obj$condition_v586 <- canonical_condition(
  obj@meta.data[[condition_col]]
)

obj$layer1_v586 <- as.character(
  obj@meta.data[[ANNOTATION_COL]]
)

msg(
  "Cells: ",
  ncol(obj)
)

msg(
  "UMAP: ",
  umap_use
)

msg(
  "Condition source: ",
  condition_col
)

msg(
  "Layer1 source: ",
  ANNOTATION_COL
)

msg(
  "Layer1 counts:"
)

print(
  table(
    obj$layer1_v586,
    useNA = "ifany"
  )
)

# ==============================================================================
# 5. Gene audit
# ==============================================================================

available <- rownames(obj)
available_lower <- tolower(available)

resolved_genes <- vapply(
  TARGET_GENES,
  FUN.VALUE = character(1),
  FUN = function(g) {
    ii <- match(
      tolower(g),
      available_lower
    )

    if (is.na(ii)) {
      NA_character_
    } else {
      available[[ii]]
    }
  }
)

group_lookup <- tibble(
  requested_gene = TARGET_GENES
) %>%
  rowwise() %>%
  mutate(
    subtype_group = {
      hit <- names(
        CANDIDATE_GROUPS
      )[
        vapply(
          CANDIDATE_GROUPS,
          function(x) requested_gene %in% x,
          logical(1)
        )
      ]

      if (length(hit)) hit[[1]] else NA_character_
    }
  ) %>%
  ungroup()

gene_audit <- tibble(
  requested_gene = TARGET_GENES,
  resolved_gene = resolved_genes,
  present = !is.na(resolved_genes)
) %>%
  left_join(
    group_lookup,
    by = "requested_gene"
  )

write.csv(
  gene_audit,
  file.path(
    TAB_DIR,
    "00_gene_audit_v5.8.6.csv"
  ),
  row.names = FALSE
)

print(
  gene_audit
)

genes_use <- gene_audit %>%
  filter(present) %>%
  pull(resolved_gene) %>%
  unique()

if (!length(genes_use)) {
  stop(
    "None of the requested genes are present."
  )
}

msg(
  "Genes available: ",
  paste(
    genes_use,
    collapse = ", "
  )
)

# ==============================================================================
# 6. Matrices
# ==============================================================================

counts_mat <- GetAssayData(
  obj,
  assay = ASSAY_USE,
  layer = "counts"
)

data_mat <- GetAssayData(
  obj,
  assay = ASSAY_USE,
  layer = "data"
)

layer1 <- obj$layer1_v586

celltypes <- sort(
  unique(
    layer1
  )
)

# ------------------------------------------------------------------
# v5.8.6: robust whole-cell Layer1 label resolution
# ------------------------------------------------------------------

resolve_layer1_label <- function(
  celltypes,
  exact_label,
  regex_priority,
  role_name
) {

  if (
    exact_label %in%
      celltypes
  ) {

    msg(
      role_name,
      " label resolved by exact match: ",
      exact_label
    )

    return(
      exact_label
    )
  }

  for (
    pattern in regex_priority
  ) {

    hits <- celltypes[
      grepl(
        pattern,
        celltypes,
        ignore.case = TRUE
      )
    ]

    if (
      length(
        hits
      ) > 0
    ) {

      # If more than one label matches, choose the largest population.
      if (
        length(
          hits
        ) > 1
      ) {

        hit_counts <- vapply(
          hits,
          function(h) {
            sum(
              layer1 ==
                h,
              na.rm = TRUE
            )
          },
          numeric(1)
        )

        hits <- hits[
          order(
            hit_counts,
            decreasing = TRUE
          )
        ]
      }

      msg(
        role_name,
        " label auto-resolved: ",
        hits[[1]]
      )

      if (
        length(
          hits
        ) > 1
      ) {
        msg(
          role_name,
          " additional matching labels: ",
          paste(
            hits[-1],
            collapse = ", "
          )
        )
      }

      return(
        hits[[1]]
      )
    }
  }

  stop(
    role_name,
    " label could not be resolved. Available Layer1 labels: ",
    paste(
      celltypes,
      collapse = " | "
    )
  )
}

TARGET_KUPFFER_LABEL <- resolve_layer1_label(
  celltypes = celltypes,
  exact_label = TARGET_KUPFFER_LABEL,
  regex_priority = c(
    "^Kupffer[_ /-]*Macrophage$",
    "^Kupffer",
    "Kupffer",
    "Macrophage",
    "Mphi",
    "MΦ"
  ),
  role_name = "Kupffer/Macrophage"
)

if (
  TARGET_MONOCYTE_LABEL %in%
    celltypes
) {

  msg(
    "Monocyte label resolved by exact match: ",
    TARGET_MONOCYTE_LABEL
  )

} else {

  monocyte_hits <- celltypes[
    grepl(
      "Monocyte",
      celltypes,
      ignore.case = TRUE
    )
  ]

  if (
    length(
      monocyte_hits
    ) > 0
  ) {

    monocyte_counts <- vapply(
      monocyte_hits,
      function(h) {
        sum(
          layer1 ==
            h,
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    TARGET_MONOCYTE_LABEL <- monocyte_hits[
      order(
        monocyte_counts,
        decreasing = TRUE
      )
    ][[1]]

    msg(
      "Monocyte label auto-resolved: ",
      TARGET_MONOCYTE_LABEL
    )

  } else {

    msg(
      "No Monocyte Layer1 label found; Kupffer-only specificity metrics will still be calculated."
    )

    TARGET_MONOCYTE_LABEL <- NA_character_
  }
}

write.csv(
  tibble(
    role = c(
      "Kupffer/Macrophage",
      "Monocyte"
    ),
    resolved_label = c(
      TARGET_KUPFFER_LABEL,
      TARGET_MONOCYTE_LABEL
    )
  ),
  file.path(
    TAB_DIR,
    "00b_resolved_Layer1_labels_v5.8.6.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 7. Full-cell UMAP overview
# ==============================================================================

msg(
  "Generating full-cell UMAP FeaturePlots..."
)

p_umap <- FeaturePlot(
  obj,
  features = genes_use,
  reduction = umap_use,
  ncol = 3,
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q98",
  raster = FALSE,
  pt.size = 0.18
) &
  scale_colour_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) &
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_umap,
  file.path(
    FIG_DIR,
    "01_FULL_CELL_candidate_UMAP_v5.8.6.pdf"
  ),
  width = 14,
  height = 4.2 *
    ceiling(
      length(genes_use) / 3
    )
)

# ==============================================================================
# 8. Layer1 DotPlot
# ==============================================================================

msg(
  "Generating Layer1 DotPlot..."
)

obj$layer1_v586_factor <- factor(
  obj$layer1_v586,
  levels = celltypes
)

p_dot <- DotPlot(
  obj,
  features = genes_use,
  group.by = "layer1_v586_factor",
  assay = ASSAY_USE
) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  RotatedAxis() +
  labs(
    title =
      "Whole-liver specificity of macrophage-subtype staining candidates",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_pdf(
  p_dot,
  file.path(
    FIG_DIR,
    "02_LAYER1_DotPlot_all_candidates_v5.8.6.pdf"
  ),
  width = 12,
  height = 8
)

# ==============================================================================
# 9. Quantitative expression by Layer1
# ==============================================================================

msg(
  "Calculating expression by Layer1..."
)

rows <- list()
k <- 1L

for (gene in genes_use) {

  for (ct in celltypes) {

    idx <- which(
      layer1 ==
        ct
    )

    gene_counts <- as.numeric(
      counts_mat[
        gene,
        idx,
        drop = TRUE
      ]
    )

    gene_data <- as.numeric(
      data_mat[
        gene,
        idx,
        drop = TRUE
      ]
    )

    rows[[k]] <- tibble(
      gene = gene,
      celltype = ct,
      n_cells = length(idx),
      n_positive = sum(
        gene_counts > 0
      ),
      pct_positive = 100 *
        mean(
          gene_counts > 0
        ),
      mean_expression = mean(
        gene_data
      ),
      mean_expression_positive =
        ifelse(
          any(
            gene_counts > 0
          ),
          mean(
            gene_data[
              gene_counts > 0
            ]
          ),
          NA_real_
        )
    )

    k <- k + 1L
  }
}

by_celltype <- bind_rows(
  rows
)

group_map_resolved <- gene_audit %>%
  filter(present) %>%
  select(
    gene = resolved_gene,
    subtype_group
  )

by_celltype <- by_celltype %>%
  left_join(
    group_map_resolved,
    by = "gene"
  )

write.csv(
  by_celltype,
  file.path(
    TAB_DIR,
    "01_candidate_expression_by_Layer1_v5.8.6.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 10. Kupffer / macrophage specificity metrics
# ==============================================================================

msg(
  "Calculating Kupffer specificity metrics..."
)

specificity_rows <- list()

for (gene in genes_use) {

  gene_counts_all <- as.numeric(
    counts_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  gene_data_all <- as.numeric(
    data_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  is_pos <- gene_counts_all > 0

  is_kupffer <- layer1 ==
    TARGET_KUPFFER_LABEL

  is_monocyte <- if (
    !is.na(
      TARGET_MONOCYTE_LABEL
    ) &&
    TARGET_MONOCYTE_LABEL %in%
      celltypes
  ) {
    layer1 ==
      TARGET_MONOCYTE_LABEL
  } else {
    rep(
      FALSE,
      length(layer1)
    )
  }

  is_kupffer_or_monocyte <-
    is_kupffer |
    is_monocyte

  pct_kupffer <- 100 *
    mean(
      is_pos[
        is_kupffer
      ]
    )

  pct_non_kupffer <- 100 *
    mean(
      is_pos[
        !is_kupffer
      ]
    )

  pct_kupffer_monocyte <- 100 *
    mean(
      is_pos[
        is_kupffer_or_monocyte
      ]
    )

  pct_other <- 100 *
    mean(
      is_pos[
        !is_kupffer_or_monocyte
      ]
    )

  mean_kupffer <- mean(
    gene_data_all[
      is_kupffer
    ]
  )

  mean_non_kupffer <- mean(
    gene_data_all[
      !is_kupffer
    ]
  )

  ppv_kupffer <- if (
    sum(is_pos) > 0
  ) {
    100 *
      sum(
        is_pos &
          is_kupffer
      ) /
      sum(
        is_pos
      )
  } else {
    NA_real_
  }

  ppv_kupffer_or_monocyte <- if (
    sum(is_pos) > 0
  ) {
    100 *
      sum(
        is_pos &
          is_kupffer_or_monocyte
      ) /
      sum(
        is_pos
      )
  } else {
    NA_real_
  }

  non_kupffer_tbl <- by_celltype %>%
    filter(
      gene == !!gene,
      celltype !=
        TARGET_KUPFFER_LABEL
    ) %>%
    arrange(
      desc(
        pct_positive
      )
    )

  top_non_kupffer_type <- if (
    nrow(
      non_kupffer_tbl
    )
  ) {
    non_kupffer_tbl$celltype[[1]]
  } else {
    NA_character_
  }

  top_non_kupffer_pct <- if (
    nrow(
      non_kupffer_tbl
    )
  ) {
    non_kupffer_tbl$pct_positive[[1]]
  } else {
    NA_real_
  }

  top_non_kupffer_mean <- if (
    nrow(
      non_kupffer_tbl
    )
  ) {
    non_kupffer_tbl$mean_expression[[1]]
  } else {
    NA_real_
  }

  specificity_rows[[gene]] <- tibble(
    gene = gene,

    pct_positive_Kupffer =
      pct_kupffer,

    pct_positive_nonKupffer =
      pct_non_kupffer,

    pct_positive_enrichment_Kupffer_vs_nonKupffer =
      (
        pct_kupffer + 0.01
      ) /
        (
          pct_non_kupffer + 0.01
        ),

    mean_expression_Kupffer =
      mean_kupffer,

    mean_expression_nonKupffer =
      mean_non_kupffer,

    mean_expression_enrichment_Kupffer_vs_nonKupffer =
      (
        mean_kupffer + 1e-6
      ) /
        (
          mean_non_kupffer + 1e-6
        ),

    pct_of_ALL_positive_cells_that_are_Kupffer =
      ppv_kupffer,

    pct_positive_Kupffer_or_Monocyte =
      pct_kupffer_monocyte,

    pct_positive_all_other_cells =
      pct_other,

    pct_positive_enrichment_KupfferMonocyte_vs_rest =
      (
        pct_kupffer_monocyte + 0.01
      ) /
        (
          pct_other + 0.01
        ),

    pct_of_ALL_positive_cells_that_are_Kupffer_or_Monocyte =
      ppv_kupffer_or_monocyte,

    top_nonKupffer_celltype_by_pct =
      top_non_kupffer_type,

    top_nonKupffer_pct_positive =
      top_non_kupffer_pct,

    top_nonKupffer_mean_expression =
      top_non_kupffer_mean
  )
}

specificity <- bind_rows(
  specificity_rows
) %>%
  left_join(
    group_map_resolved,
    by = "gene"
  ) %>%
  arrange(
    subtype_group,
    desc(
      pct_of_ALL_positive_cells_that_are_Kupffer
    )
  )

write.csv(
  specificity,
  file.path(
    TAB_DIR,
    "02_Kupffer_specificity_metrics_v5.8.6.csv"
  ),
  row.names = FALSE
)

print(
  specificity
)

# ==============================================================================
# 11. Positive-cell composition
# ==============================================================================

composition <- by_celltype %>%
  group_by(
    gene
  ) %>%
  mutate(
    fraction_of_gene_positive_cells =
      ifelse(
        sum(
          n_positive
        ) > 0,
        100 *
          n_positive /
          sum(
            n_positive
          ),
        NA_real_
      )
  ) %>%
  ungroup()

write.csv(
  composition,
  file.path(
    TAB_DIR,
    "03_positive_cell_composition_by_Layer1_v5.8.6.csv"
  ),
  row.names = FALSE
)

p_comp <- composition %>%
  filter(
    n_positive > 0
  ) %>%
  ggplot(
    aes(
      x = gene,
      y =
        fraction_of_gene_positive_cells,
      fill = celltype
    )
  ) +
  geom_col(
    width = 0.75
  ) +
  labs(
    title =
      "Cell-type composition of RNA-positive cells",
    x = NULL,
    y =
      "% of all RNA-positive cells",
    fill = "Layer1"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

save_pdf(
  p_comp,
  file.path(
    FIG_DIR,
    "03_positive_cell_composition_all_candidates_v5.8.6.pdf"
  ),
  width = 12,
  height = 7
)

# ==============================================================================
# 12. Heatmap: percent positive
# ==============================================================================

heat_pct <- by_celltype %>%
  select(
    gene,
    celltype,
    pct_positive
  )

heat_pct$gene <- factor(
  heat_pct$gene,
  levels = genes_use
)

p_heat_pct <- ggplot(
  heat_pct,
  aes(
    x = gene,
    y = celltype,
    fill = pct_positive
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  labs(
    title =
      "RNA-positive fraction across whole-liver cell types",
    x = NULL,
    y = NULL,
    fill = "% positive"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

save_pdf(
  p_heat_pct,
  file.path(
    FIG_DIR,
    "04_LAYER1_positive_fraction_heatmap_v5.8.6.pdf"
  ),
  width = 12,
  height = 8
)

# ==============================================================================
# 13. Heatmap: mean expression
# ==============================================================================

heat_mean <- by_celltype %>%
  select(
    gene,
    celltype,
    mean_expression
  )

heat_mean$gene <- factor(
  heat_mean$gene,
  levels = genes_use
)

p_heat_mean <- ggplot(
  heat_mean,
  aes(
    x = gene,
    y = celltype,
    fill = mean_expression
  )
) +
  geom_tile(
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint =
      median(
        heat_mean$mean_expression,
        na.rm = TRUE
      )
  ) +
  labs(
    title =
      "Mean RNA expression across whole-liver cell types",
    x = NULL,
    y = NULL,
    fill = "Mean expression"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )

save_pdf(
  p_heat_mean,
  file.path(
    FIG_DIR,
    "05_LAYER1_mean_expression_heatmap_v5.8.6.pdf"
  ),
  width = 12,
  height = 8
)

# ==============================================================================
# 14. Sham-vs-Tx full-cell UMAPs
# ==============================================================================

msg(
  "Generating Sham-vs-Tx full-cell UMAPs..."
)

obj$ShamTx_v586 <- factor(
  ifelse(
    as.character(
      obj$condition_v586
    ) %in%
      c(
        "Sham",
        "Tx"
      ),
    as.character(
      obj$condition_v586
    ),
    NA_character_
  ),
  levels = c(
    "Sham",
    "Tx"
  )
)

shamtx_cells <- colnames(obj)[
  !is.na(
    obj$ShamTx_v586
  )
]

obj_shamtx <- subset(
  obj,
  cells =
    shamtx_cells
)

for (gene in genes_use) {

  p_shamtx <- FeaturePlot(
    obj_shamtx,
    features = gene,
    reduction = umap_use,
    split.by = "ShamTx_v586",
    keep.scale = "all",
    ncol = 2,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q98",
    raster = FALSE,
    pt.size = 0.16
  ) &
    scale_colour_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) &
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_shamtx,
    file.path(
      FIG_DIR,
      paste0(
        "06_SHAM_vs_TX_FULL_CELL_",
        gene,
        "_v5.8.6.pdf"
      )
    ),
    width = 11,
    height = 5
  )
}

# ==============================================================================
# 15. Summary ranking for staining specificity
# ==============================================================================

specificity_summary <- specificity %>%
  mutate(
    score_Kupffer_PPV =
      scales::rescale(
        pct_of_ALL_positive_cells_that_are_Kupffer,
        to = c(
          0,
          1
        )
      ),

    score_Kupffer_enrichment =
      scales::rescale(
        log2(
          pct_positive_enrichment_Kupffer_vs_nonKupffer +
            1e-6
        ),
        to = c(
          0,
          1
        )
      ),

    score_Kupffer_detection =
      scales::rescale(
        pct_positive_Kupffer,
        to = c(
          0,
          1
        )
      ),

    wholecell_Kupffer_specificity_score =
      0.45 *
        score_Kupffer_PPV +
      0.35 *
        score_Kupffer_enrichment +
      0.20 *
        score_Kupffer_detection
  ) %>%
  arrange(
    desc(
      wholecell_Kupffer_specificity_score
    )
  ) %>%
  mutate(
    wholecell_specificity_rank =
      row_number()
  )

write.csv(
  specificity_summary,
  file.path(
    TAB_DIR,
    "04_FINAL_wholecell_Kupffer_specificity_ranking_v5.8.6.csv"
  ),
  row.names = FALSE
)

cat(
  "\n============================================================\n"
)

cat(
  "Whole-cell Kupffer specificity ranking v5.8.6\n"
)

cat(
  "============================================================\n\n"
)

print(
  specificity_summary %>%
    select(
      wholecell_specificity_rank,
      subtype_group,
      gene,
      wholecell_Kupffer_specificity_score,
      pct_positive_Kupffer,
      pct_positive_nonKupffer,
      pct_positive_enrichment_Kupffer_vs_nonKupffer,
      pct_of_ALL_positive_cells_that_are_Kupffer,
      pct_of_ALL_positive_cells_that_are_Kupffer_or_Monocyte,
      top_nonKupffer_celltype_by_pct,
      top_nonKupffer_pct_positive
    )
)

# ==============================================================================
# 16. Clean-B macrophage reference-marker UMAPs
# ==============================================================================

if (
  length(MPHI_RDS_PATH) == 0L ||
  is.na(MPHI_RDS_PATH) ||
  !file.exists(MPHI_RDS_PATH)
) {
  msg(
    "Clean-B macrophage RDS not found; skipping Mphi reference UMAP panels."
  )
} else {

  msg(
    "Loading Clean-B macrophage RDS for reference-marker UMAPs..."
  )

  mphi_ref <- readRDS(
    MPHI_RDS_PATH
  )

  DefaultAssay(
    mphi_ref
  ) <- ASSAY_USE

  mphi_ref_umap <- resolve_first(
    Reductions(mphi_ref),
    c(
      "umapRPCA",
      "umap"
    )
  )

  if (is.na(mphi_ref_umap)) {
    stop(
      "No usable UMAP reduction in Clean-B macrophage RDS."
    )
  }

  if (!"condition_FIXED2" %in% colnames(mphi_ref@meta.data)) {
    stop(
      "condition_FIXED2 missing from Clean-B macrophage RDS."
    )
  }

  mphi_ref$condition_v586 <- canonical_condition(
    mphi_ref@meta.data[["condition_FIXED2"]]
  )

  ref_present <- REFERENCE_MPHI_MARKERS[
    REFERENCE_MPHI_MARKERS %in%
      rownames(mphi_ref)
  ]

  ref_missing <- setdiff(
    REFERENCE_MPHI_MARKERS,
    ref_present
  )

  write.csv(
    tibble(
      gene = REFERENCE_MPHI_MARKERS,
      display_name = unname(
        REFERENCE_MPHI_DISPLAY[
          REFERENCE_MPHI_MARKERS
        ]
      ),
      present = REFERENCE_MPHI_MARKERS %in%
        rownames(mphi_ref)
    ),
    file.path(
      TAB_DIR,
      "05_reference_Mphi_marker_gene_audit_v5.8.6.csv"
    ),
    row.names = FALSE
  )

  if (length(ref_missing)) {
    msg(
      "Missing reference genes in Clean-B Mphi RDS: ",
      paste(
        ref_missing,
        collapse = ", "
      )
    )
  }

  if (length(ref_present)) {

    # --------------------------------------------------------------------------
    # 16A. Full Clean-B macrophage UMAP
    # --------------------------------------------------------------------------

    p_ref_mphi <- FeaturePlot(
      mphi_ref,
      features = ref_present,
      reduction = mphi_ref_umap,
      ncol = 4,
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q98",
      raster = FALSE,
      pt.size = 0.28
    ) &
      scale_colour_gradientn(
        colours = c(
          "#0033FF",
          "#FFFFFF",
          "#FF1A1A"
        )
      ) &
      theme_classic(
        base_size = 8
      )

    save_pdf(
      p_ref_mphi,
      file.path(
        FIG_DIR,
        "07_REFERENCE_Cd68_Cd80_CD206_Cd86_M2markers_FULL_MPHI_UMAP_v5.8.6.pdf"
      ),
      width = 15,
      height = 7.5
    )

    # --------------------------------------------------------------------------
    # 16B. Sham vs Tx on Clean-B macrophage UMAP
    # --------------------------------------------------------------------------

    mphi_ref_shamtx_cells <- colnames(mphi_ref)[
      as.character(
        mphi_ref$condition_v586
      ) %in%
        c(
          "Sham",
          "Tx"
        )
    ]

    mphi_ref_shamtx <- subset(
      mphi_ref,
      cells = mphi_ref_shamtx_cells
    )

    mphi_ref_shamtx$ShamTx_ref_v586 <- factor(
      as.character(
        mphi_ref_shamtx$condition_v586
      ),
      levels = c(
        "Sham",
        "Tx"
      )
    )

    for (gene in ref_present) {

      p_ref_shamtx <- FeaturePlot(
        mphi_ref_shamtx,
        features = gene,
        reduction = mphi_ref_umap,
        split.by = "ShamTx_ref_v586",
        keep.scale = "all",
        ncol = 2,
        order = TRUE,
        min.cutoff = "q05",
        max.cutoff = "q98",
        raster = FALSE,
        pt.size = 0.22
      ) &
        scale_colour_gradientn(
          colours = c(
            "#0033FF",
            "#FFFFFF",
            "#FF1A1A"
          )
        ) &
        theme_classic(
          base_size = 8
        )

      save_pdf(
        p_ref_shamtx,
        file.path(
          FIG_DIR,
          paste0(
            "08_REFERENCE_MPHI_SHAM_vs_TX_",
            gene,
            "_v5.8.6.pdf"
          )
        ),
        width = 11,
        height = 5
      )
    }
  }
}


# ==============================================================================
# 17. Focused ECM-associated inflammatory-Mphi validation
# ==============================================================================

if (
  length(MPHI_RDS_PATH) == 0L ||
  is.na(MPHI_RDS_PATH) ||
  !file.exists(MPHI_RDS_PATH)
) {

  msg(
    "Clean-B macrophage RDS unavailable; skipping focused ECM validation."
  )

} else {

  # Reuse mphi_ref if section 16 already loaded it.
  if (!exists("mphi_ref")) {

    msg(
      "Loading Clean-B macrophage RDS for focused ECM validation..."
    )

    mphi_ref <- readRDS(
      MPHI_RDS_PATH
    )

    DefaultAssay(
      mphi_ref
    ) <- ASSAY_USE

    mphi_ref_umap <- resolve_first(
      Reductions(mphi_ref),
      c(
        "umapRPCA",
        "umap"
      )
    )

    if (is.na(mphi_ref_umap)) {
      stop(
        "No usable UMAP reduction in Clean-B macrophage RDS."
      )
    }

    if (!"condition_FIXED2" %in% colnames(mphi_ref@meta.data)) {
      stop(
        "condition_FIXED2 missing from Clean-B macrophage RDS."
      )
    }

    mphi_ref$condition_v586 <- canonical_condition(
      mphi_ref@meta.data[["condition_FIXED2"]]
    )
  }

  if (!MPHI_SUBTYPE_COL %in% colnames(mphi_ref@meta.data)) {
    stop(
      "Final Clean-B macrophage subtype column missing: ",
      MPHI_SUBTYPE_COL
    )
  }

  mphi_sample_col <- resolve_first(
    colnames(mphi_ref@meta.data),
    MPHI_SAMPLE_CANDIDATES
  )

  if (is.na(mphi_sample_col)) {
    stop(
      "No usable sample column found in Clean-B macrophage RDS."
    )
  }

  mphi_ref$subtype_v586 <- factor(
    as.character(
      mphi_ref@meta.data[[MPHI_SUBTYPE_COL]]
    ),
    levels = MPHI_SUBTYPE_LEVELS
  )

  mphi_ref$sample_v586 <- as.character(
    mphi_ref@meta.data[[mphi_sample_col]]
  )

  ecm_present <- ECM_FOCUS_MARKERS[
    ECM_FOCUS_MARKERS %in%
      rownames(mphi_ref)
  ]

  ecm_missing <- setdiff(
    ECM_FOCUS_MARKERS,
    ecm_present
  )

  write.csv(
    tibble(
      gene = ECM_FOCUS_MARKERS,
      display_name = unname(
        ECM_FOCUS_DISPLAY[
          ECM_FOCUS_MARKERS
        ]
      ),
      present = ECM_FOCUS_MARKERS %in%
        rownames(mphi_ref)
    ),
    file.path(
      TAB_DIR,
      "06_ECM_focus_gene_audit_v5.8.6.csv"
    ),
    row.names = FALSE
  )

  if (length(ecm_missing)) {
    msg(
      "Missing ECM focus genes in Clean-B Mphi RDS: ",
      paste(
        ecm_missing,
        collapse = ", "
      )
    )
  }

  if (!length(ecm_present)) {
    stop(
      "None of the focused ECM-associated inflammatory-Mphi genes are present."
    )
  }

  # --------------------------------------------------------------------------
  # 17A. Full Clean-B Mphi UMAP: F13A1 / THBS1 / TGFBI / VCAN / CCR2
  # --------------------------------------------------------------------------

  msg(
    "Generating focused ECM-associated inflammatory-Mphi full UMAP..."
  )

  p_ecm_full <- FeaturePlot(
    mphi_ref,
    features = ecm_present,
    reduction = mphi_ref_umap,
    ncol = 3,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q98",
    raster = FALSE,
    pt.size = 0.30
  ) &
    scale_colour_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) &
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_ecm_full,
    file.path(
      FIG_DIR,
      "09_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_FULL_MPHI_UMAP_v5.8.6.pdf"
    ),
    width = 13,
    height = 8
  )

  # --------------------------------------------------------------------------
  # 17B. Sham vs Tx split UMAP
  # --------------------------------------------------------------------------

  ecm_shamtx_cells <- colnames(mphi_ref)[
    as.character(
      mphi_ref$condition_v586
    ) %in%
      c(
        "Sham",
        "Tx"
      )
  ]

  mphi_ecm_shamtx <- subset(
    mphi_ref,
    cells = ecm_shamtx_cells
  )

  mphi_ecm_shamtx$ShamTx_ECM_v586 <- factor(
    as.character(
      mphi_ecm_shamtx$condition_v586
    ),
    levels = c(
      "Sham",
      "Tx"
    )
  )

  msg(
    "Generating focused ECM-associated inflammatory-Mphi Sham-vs-Tx UMAPs..."
  )

  for (gene in ecm_present) {

    p_ecm_split <- FeaturePlot(
      mphi_ecm_shamtx,
      features = gene,
      reduction = mphi_ref_umap,
      split.by = "ShamTx_ECM_v586",
      keep.scale = "all",
      ncol = 2,
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q98",
      raster = FALSE,
      pt.size = 0.24
    ) &
      scale_colour_gradientn(
        colours = c(
          "#0033FF",
          "#FFFFFF",
          "#FF1A1A"
        )
      ) &
      theme_classic(
        base_size = 8
      )

    save_pdf(
      p_ecm_split,
      file.path(
        FIG_DIR,
        paste0(
          "10_ECM_SHAM_vs_TX_",
          gene,
          "_v5.8.6.pdf"
        )
      ),
      width = 11,
      height = 5
    )
  }

  # --------------------------------------------------------------------------
  # 17C. Five-subtype DotPlot
  # --------------------------------------------------------------------------

  msg(
    "Generating focused ECM marker subtype DotPlot..."
  )

  p_ecm_dot <- DotPlot(
    mphi_ref,
    features = ecm_present,
    group.by = "subtype_v586",
    assay = ASSAY_USE
  ) +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    RotatedAxis() +
    labs(
      title =
        "ECM-associated inflammatory-Mphi marker specificity across Mphi subtypes",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p_ecm_dot,
    file.path(
      FIG_DIR,
      "11_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_subtype_DotPlot_v5.8.6.pdf"
    ),
    width = 9,
    height = 5.8
  )

  # --------------------------------------------------------------------------
  # 17D. Five-subtype violin plots
  # --------------------------------------------------------------------------

  msg(
    "Generating focused ECM marker subtype Violin plots..."
  )

  p_ecm_vln <- VlnPlot(
    mphi_ref,
    features = ecm_present,
    group.by = "subtype_v586",
    assay = ASSAY_USE,
    pt.size = 0,
    ncol = 2,
    combine = TRUE
  ) &
    theme_classic(
      base_size = 8
    ) &
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_pdf(
    p_ecm_vln,
    file.path(
      FIG_DIR,
      "12_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_subtype_Violin_v5.8.6.pdf"
    ),
    width = 13,
    height = 10
  )

  # --------------------------------------------------------------------------
  # 17E. Sample-level quantitative metrics:
  #      Sham1 / Sham20 / Tx17 / Tx5
  # --------------------------------------------------------------------------

  msg(
    "Calculating focused ECM marker sample-level Sham/Tx metrics..."
  )

  mphi_counts <- GetAssayData(
    mphi_ref,
    assay = ASSAY_USE,
    layer = "counts"
  )

  mphi_data <- GetAssayData(
    mphi_ref,
    assay = ASSAY_USE,
    layer = "data"
  )

  target_samples <- unique(
    mphi_ref$sample_v586[
      as.character(
        mphi_ref$condition_v586
      ) %in%
        c(
          "Sham",
          "Tx"
        )
    ]
  )

  target_samples <- target_samples[
    !is.na(target_samples)
  ]

  sample_metric_rows <- list()
  rr <- 1L

  for (sample_now in target_samples) {

    idx <- which(
      mphi_ref$sample_v586 ==
        sample_now
    )

    if (!length(idx)) {
      next
    }

    condition_now <- unique(
      as.character(
        mphi_ref$condition_v586[
          idx
        ]
      )
    )

    condition_now <- condition_now[
      !is.na(condition_now)
    ]

    condition_now <- if (
      length(condition_now)
    ) {
      condition_now[[1]]
    } else {
      NA_character_
    }

    library_size_now <- sum(
      Matrix::colSums(
        mphi_counts[
          ,
          idx,
          drop = FALSE
        ]
      )
    )

    for (gene in ecm_present) {

      gene_counts <- as.numeric(
        mphi_counts[
          gene,
          idx,
          drop = TRUE
        ]
      )

      gene_data <- as.numeric(
        mphi_data[
          gene,
          idx,
          drop = TRUE
        ]
      )

      gene_sum_counts <- sum(
        gene_counts
      )

      sample_metric_rows[[rr]] <- tibble(
        gene = gene,
        sample = sample_now,
        condition = condition_now,
        n_Mphi = length(idx),
        n_positive = sum(
          gene_counts > 0
        ),
        pct_positive = 100 *
          mean(
            gene_counts > 0
          ),
        mean_log_normalized_expression = mean(
          gene_data
        ),
        sum_gene_counts = gene_sum_counts,
        total_RNA_counts = library_size_now,
        pseudobulk_CPM = ifelse(
          library_size_now > 0,
          1e6 *
            gene_sum_counts /
            library_size_now,
          NA_real_
        )
      )

      rr <- rr + 1L
    }
  }

  ecm_sample_metrics <- bind_rows(
    sample_metric_rows
  )

  write.csv(
    ecm_sample_metrics,
    file.path(
      TAB_DIR,
      "07_ECM_focus_samplelevel_Shams_Txs_metrics_v5.8.6.csv"
    ),
    row.names = FALSE
  )

  # --------------------------------------------------------------------------
  # 17F. Sham-vs-Tx summary derived from biological samples
  # --------------------------------------------------------------------------

  ecm_condition_summary <- ecm_sample_metrics %>%
    group_by(
      gene,
      condition
    ) %>%
    summarise(
      n_samples = n(),
      mean_pct_positive = mean(
        pct_positive,
        na.rm = TRUE
      ),
      mean_expression = mean(
        mean_log_normalized_expression,
        na.rm = TRUE
      ),
      mean_pseudobulk_CPM = mean(
        pseudobulk_CPM,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  ecm_shamtx_summary <- ecm_condition_summary %>%
    filter(
      condition %in%
        c(
          "Sham",
          "Tx"
        )
    ) %>%
    select(
      gene,
      condition,
      n_samples,
      mean_pct_positive,
      mean_expression,
      mean_pseudobulk_CPM
    ) %>%
    pivot_wider(
      names_from = condition,
      values_from = c(
        n_samples,
        mean_pct_positive,
        mean_expression,
        mean_pseudobulk_CPM
      )
    ) %>%
    mutate(
      delta_pct_positive_Tx_minus_Sham =
        mean_pct_positive_Tx -
        mean_pct_positive_Sham,

      log2FC_pseudobulk_CPM_Tx_vs_Sham =
        log2(
          (
            mean_pseudobulk_CPM_Tx +
              1e-6
          ) /
            (
              mean_pseudobulk_CPM_Sham +
                1e-6
            )
        ),

      log2FC_mean_expression_Tx_vs_Sham =
        log2(
          (
            mean_expression_Tx +
              1e-6
          ) /
            (
              mean_expression_Sham +
                1e-6
            )
        )
    ) %>%
    arrange(
      log2FC_pseudobulk_CPM_Tx_vs_Sham
    )

  write.csv(
    ecm_shamtx_summary,
    file.path(
      TAB_DIR,
      "08_ECM_focus_Shams_vs_Txs_summary_v5.8.6.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\n============================================================\n"
  )

  cat(
    "ECM-associated inflammatory-Mphi focused Sham vs Tx summary v5.8.6\n"
  )

  cat(
    "============================================================\n\n"
  )

  print(
    ecm_shamtx_summary
  )

  # --------------------------------------------------------------------------
  # 17G. Sample-level positive-fraction plot
  # --------------------------------------------------------------------------

  sample_order <- c(
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )

  observed_order <- sample_order[
    sample_order %in%
      ecm_sample_metrics$sample
  ]

  remaining_samples <- setdiff(
    unique(
      ecm_sample_metrics$sample
    ),
    observed_order
  )

  ecm_sample_metrics$sample_plot <- factor(
    ecm_sample_metrics$sample,
    levels = c(
      observed_order,
      remaining_samples
    )
  )

  p_ecm_sample_pct <- ggplot(
    ecm_sample_metrics,
    aes(
      x = sample_plot,
      y = pct_positive,
      group = gene
    )
  ) +
    geom_line(
      aes(
        linetype = gene
      ),
      linewidth = 0.7
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~ gene,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title =
        "ECM-associated inflammatory-Mphi markers: RNA-positive fraction by biological sample",
      x = NULL,
      y = "% RNA-positive Mphi"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p_ecm_sample_pct,
    file.path(
      FIG_DIR,
      "13_ECM_focus_positive_fraction_by_sample_v5.8.6.pdf"
    ),
    width = 11,
    height = 7.5
  )

  # --------------------------------------------------------------------------
  # 17H. Sample-level pseudobulk CPM plot
  # --------------------------------------------------------------------------

  p_ecm_sample_cpm <- ggplot(
    ecm_sample_metrics,
    aes(
      x = sample_plot,
      y = pseudobulk_CPM,
      group = gene
    )
  ) +
    geom_line(
      aes(
        linetype = gene
      ),
      linewidth = 0.7
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~ gene,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title =
        "ECM-associated inflammatory-Mphi markers: pseudobulk CPM by biological sample",
      x = NULL,
      y = "Pseudobulk CPM"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p_ecm_sample_cpm,
    file.path(
      FIG_DIR,
      "14_ECM_focus_pseudobulk_CPM_by_sample_v5.8.6.pdf"
    ),
    width = 11,
    height = 7.5
  )
}

# ==============================================================================
# 18. Metadata / index

# ==============================================================================

metadata_out <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "assay",
    "UMAP",
    "annotation_column",
    "condition_source",
    "Kupffer_label",
    "Monocyte_label",
    "CleanB_Mphi_RDS"
  ),
  value = c(
    "v5.8.6",
    RDS_PATH,
    ASSAY_USE,
    umap_use,
    ANNOTATION_COL,
    condition_col,
    TARGET_KUPFFER_LABEL,
    TARGET_MONOCYTE_LABEL,
    MPHI_RDS_PATH
  )
)

write.csv(
  metadata_out,
  file.path(
    LOG_DIR,
    "analysis_metadata_v5.8.6.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "sessionInfo_v5.8.6.txt"
  )
)

figure_index <- tibble(
  file = c(
    "01_FULL_CELL_candidate_UMAP_v5.8.6.pdf",
    "02_LAYER1_DotPlot_all_candidates_v5.8.6.pdf",
    "03_positive_cell_composition_all_candidates_v5.8.6.pdf",
    "04_LAYER1_positive_fraction_heatmap_v5.8.6.pdf",
    "05_LAYER1_mean_expression_heatmap_v5.8.6.pdf",
    "07_REFERENCE_Cd68_Cd80_CD206_Cd86_M2markers_FULL_MPHI_UMAP_v5.8.6.pdf",
    "09_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_FULL_MPHI_UMAP_v5.8.6.pdf",
    "11_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_subtype_DotPlot_v5.8.6.pdf",
    "12_ECM_F13A1_THBS1_TGFBI_VCAN_CCR2_subtype_Violin_v5.8.6.pdf",
    "13_ECM_focus_positive_fraction_by_sample_v5.8.6.pdf",
    "14_ECM_focus_pseudobulk_CPM_by_sample_v5.8.6.pdf"
  )
)

write.csv(
  figure_index,
  file.path(
    OUT_DIR,
    "FIGURE_INDEX_v5.8.6.csv"
  ),
  row.names = FALSE
)

msg(
  "DONE."
)

msg(
  "Output: ",
  OUT_DIR
)
