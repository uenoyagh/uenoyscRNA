#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# v5.6.0-used Kupffer_Macrophage STRICT contamination re-audit
#
# Version: v5.6.4
#
# PURPOSE
#   Reproduce the exact cell-selection logic used by v5.6.0:
#       lineage == "Kupffer_Macrophage"
#       condition in c("Sham", "Tx")
#
#   Then re-audit those cells for possible:
#       - Monocyte contamination
#       - Neutrophil contamination
#       - B-cell contamination
#       - T-cell contamination
#       - NK-cell contamination
#
#   At the same time, evaluate expression of current candidate factors:
#       Marco, Aif1, Siglece, Ms4a7, Hexb, Alox5ap, Ltc4s, Ntpcr
#
# IMPORTANT
#   - No reintegration, reclustering, or re-UMAP.
#   - This script uses the frozen Clean-B macrophage RDS.
#   - The original umapRPCA embedding is preserved.
#   - Contamination calls are exploratory QC flags, not biological labels.
#   - A candidate is NOT rejected merely because it is expressed in a
#     contamination-flagged cell; the purpose is to determine whether candidate
#     signal is driven by suspicious cells or remains in clean Kupffer cells.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(5610)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

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

MPHI_RDS <- MPHI_RDS_CANDIDATES[file.exists(MPHI_RDS_CANDIDATES)][1]

V560_SAMPLE_COUNTS <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Cd163like_Sham_to_Tx_Decrease_v5.6.0",
  "Tables",
  "00_sample_cell_counts_v5.6.0.csv"
)

OUTPUT_DIR <- file.path(
  DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Kupffer_Contamination_Reaudit_v5.6.4"
)

DIR_TABLE <- file.path(OUTPUT_DIR, "Tables")
DIR_UMAP <- file.path(OUTPUT_DIR, "UMAP")
DIR_DOTPLOT <- file.path(OUTPUT_DIR, "DotPlot")
DIR_HEATMAP <- file.path(OUTPUT_DIR, "Heatmap")
DIR_LOG <- file.path(OUTPUT_DIR, "Logs")

for (d in c(
  OUTPUT_DIR,
  DIR_TABLE,
  DIR_UMAP,
  DIR_DOTPLOT,
  DIR_HEATMAP,
  DIR_LOG
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. Core settings
# ==============================================================================

ASSAY_USE <- "RNA"

TARGET_GROUPS <- c(
  "Sham",
  "Tx"
)

EXPECTED_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

MACROPHAGE_LINEAGE_VALUE <- "Kupffer_Macrophage"

LINEAGE_COLUMN_CANDIDATES <- c(
  "celltype_for_R8plot_FIXED2",
  "celltype_v440",
  "layer1_original",
  "celltype_for_R8plot",
  "celltype_auto_annotation"
)

SAMPLE_COLUMN_CANDIDATES <- c(
  "sample_for_annotation",
  "sample",
  "orig.ident"
)

CONDITION_COLUMN_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group",
  "sample_for_annotation",
  "sample"
)

CLUSTER_COLUMN_CANDIDATES <- c(
  "mphi_rpca_res_2.0",
  "mphi_rpca_cluster",
  "mphi_rpca_clusters"
)

UMAP_CANDIDATES <- c(
  "umapRPCA",
  "umap"
)

# ==============================================================================
# 3. Contamination marker sets
# ==============================================================================
#
# v5.6.4:
#   Only lineage-specific CORE markers are used for contamination calls.
#   Supportive/non-specific markers are displayed but NEVER used for flags.
#
#   This prevents macrophage MHC-II expression (Cd74/H2-Aa/H2-Ab1/Cd37)
#   from being misclassified as B-cell contamination.
# ==============================================================================

CONTAMINATION_CORE_MARKERS <- list(

  Monocyte = c(
    "Ccr2",
    "Ly6c2",
    "Plac8",
    "Sell",
    "Fcgr3",
    "Ms4a4c"
  ),

  Neutrophil = c(
    "Ly6g",
    "Csf3r",
    "Retnlg",
    "Camp",
    "Ngp",
    "Mpo",
    "Elane"
  ),

  B = c(
    "Cd79a",
    "Cd79b",
    "Ms4a1",
    "Cd19",
    "Ebf1",
    "Pax5",
    "Cd22"
  ),

  T = c(
    "Cd3d",
    "Cd3e",
    "Cd3g",
    "Trac",
    "Lck",
    "Lat",
    "Cd247",
    "Il7r"
  ),

  NK = c(
    "Ncr1",
    "Klrd1",
    "Klrk1",
    "Prf1",
    "Gzmb",
    "Ccl5"
  )
)

CONTAMINATION_SUPPORTIVE_MARKERS <- list(

  Monocyte_supportive = c(
    "S100a8",
    "S100a9"
  ),

  Neutrophil_supportive = c(
    "S100a8",
    "S100a9"
  ),

  B_supportive = c(
    "Cd74",
    "H2-Aa",
    "H2-Ab1",
    "Cd37",
    "Igkc"
  ),

  T_supportive = c(
    "Nkg7"
  ),

  NK_supportive = c(
    "Nkg7",
    "Gzma"
  )
)

# ==============================================================================
# 4. Macrophage / Kupffer identity markers for counter-check
# ==============================================================================

MACROPHAGE_IDENTITY_MARKERS <- c(
  "Adgre1",
  "Cd68",
  "Csf1r",
  "Mertk",
  "Fcgr1",
  "Lyz2",
  "Aif1",
  "Tyrobp",
  "Laptm5"
)

KUPFFER_RESIDENT_MARKERS <- c(
  "Clec4f",
  "Timd4",
  "Vsig4",
  "Marco",
  "Folr2",
  "Cd5l"
)

# ==============================================================================
# 5. Candidate factors requested/currently under consideration
# ==============================================================================

CANDIDATE_FACTORS <- c(
  "Cd163",
  "Marco",
  "Aif1",
  "Siglece",
  "Ms4a7",
  "Hexb",
  "Alox5ap",
  "Ltc4s",
  "Ntpcr"
)

# ==============================================================================
# 6. QC thresholds
# ==============================================================================
#
# Primary suspicious-cell rule:
#   lineage-specific CORE markers detected >= 2
#
# Strong suspicious-cell rule:
#   lineage-specific CORE markers detected >= 3
#
# Supportive markers never trigger contamination calls.
# ==============================================================================

PRIMARY_MIN_MARKERS <- 2
STRONG_MIN_MARKERS <- 3

# ==============================================================================
# 7. Utility functions
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
  plot,
  filename,
  width,
  height
) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
}

resolve_first <- function(
  available,
  candidates
) {
  hit <- candidates[
    candidates %in% available
  ]

  if (!length(hit)) {
    return(NA_character_)
  }

  hit[[1]]
}

canonical_condition <- function(x) {

  x <- as.character(x)

  out <- case_when(
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

resolve_condition <- function(obj) {

  meta_cols <- colnames(
    obj@meta.data
  )

  for (src in CONDITION_COLUMN_CANDIDATES) {

    if (!src %in% meta_cols) {
      next
    }

    cond <- canonical_condition(
      obj@meta.data[[src]]
    )

    if (!all(is.na(cond))) {
      return(
        list(
          source = src,
          condition = cond
        )
      )
    }
  }

  stop(
    "No usable condition metadata column found."
  )
}

# ==============================================================================
# 8. Validate and load input
# ==============================================================================

if (
  length(MPHI_RDS) == 0L ||
  is.na(MPHI_RDS) ||
  !file.exists(MPHI_RDS)
) {
  stop(
    "Clean-B macrophage RDS not found."
  )
}

msg(
  "Loading Clean-B macrophage RDS..."
)

mphi <- readRDS(
  MPHI_RDS
)

if (
  !ASSAY_USE %in%
    Assays(mphi)
) {
  stop(
    "RNA assay missing from macrophage object."
  )
}

DefaultAssay(
  mphi
) <- ASSAY_USE

# ==============================================================================
# 9. Resolve metadata exactly as in v5.6.0
# ==============================================================================

condition_info <- resolve_condition(
  mphi
)

mphi$condition_v564 <- condition_info$condition

lineage_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  LINEAGE_COLUMN_CANDIDATES
)

sample_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  SAMPLE_COLUMN_CANDIDATES
)

cluster_col <- resolve_first(
  colnames(
    mphi@meta.data
  ),
  CLUSTER_COLUMN_CANDIDATES
)

umap_name <- resolve_first(
  Reductions(mphi),
  UMAP_CANDIDATES
)

if (
  is.na(lineage_col)
) {
  stop(
    "Could not resolve lineage column."
  )
}

if (
  is.na(sample_col)
) {
  stop(
    "Could not resolve sample column."
  )
}

msg(
  "Condition source: ",
  condition_info$source
)

msg(
  "Lineage column: ",
  lineage_col
)

msg(
  "Sample column: ",
  sample_col
)

msg(
  "Cluster column: ",
  ifelse(
    is.na(cluster_col),
    "<not found>",
    cluster_col
  )
)

msg(
  "UMAP: ",
  ifelse(
    is.na(umap_name),
    "<not found>",
    umap_name
  )
)

# ==============================================================================
# 10. Reproduce v5.6.0 cell selection
# ==============================================================================

lineage <- as.character(
  mphi@meta.data[[lineage_col]]
)

condition <- as.character(
  mphi$condition_v564
)

keep <- lineage ==
  MACROPHAGE_LINEAGE_VALUE &
  condition %in%
  TARGET_GROUPS

cells_use <- colnames(mphi)[
  keep &
    !is.na(keep)
]

if (
  length(cells_use) < 100
) {
  stop(
    paste0(
      "Too few v5.6.0-like Kupffer cells: ",
      length(cells_use)
    )
  )
}

obj <- subset(
  mphi,
  cells = cells_use
)

DefaultAssay(
  obj
) <- ASSAY_USE

obj$sample_v564 <- as.character(
  obj@meta.data[[sample_col]]
)

obj$condition_v564 <- as.character(
  obj$condition_v564
)

# ==============================================================================
# 11. Verify sample counts against v5.6.0
# ==============================================================================

sample_counts_v564 <- tibble(
  sample =
    obj$sample_v564,
  condition =
    obj$condition_v564
) %>%
  count(
    condition,
    sample,
    name = "n_cells_v564"
  ) %>%
  arrange(
    condition,
    sample
  )

write.csv(
  sample_counts_v564,
  file.path(
    DIR_TABLE,
    "00_sample_cell_counts_reproduced_v5.6.4.csv"
  ),
  row.names = FALSE
)

if (
  file.exists(
    V560_SAMPLE_COUNTS
  )
) {

  sample_counts_v560 <- read.csv(
    V560_SAMPLE_COUNTS,
    stringsAsFactors = FALSE
  )

  count_check <- full_join(
    sample_counts_v560,
    sample_counts_v564,
    by = c(
      "condition",
      "sample"
    )
  ) %>%
    mutate(
      same_count =
        n_cells ==
        n_cells_v564
    )

  write.csv(
    count_check,
    file.path(
      DIR_TABLE,
      "00_v5.6.0_vs_v5.6.4_cell_count_check.csv"
    ),
    row.names = FALSE
  )

  if (
    all(
      count_check$same_count,
      na.rm = TRUE
    )
  ) {
    msg(
      "v5.6.0 sample counts reproduced successfully."
    )
  } else {
    warning(
      "v5.6.0 sample counts were NOT reproduced exactly. Inspect count-check CSV."
    )
  }
}

# Save exact cell IDs used in this audit.
cell_manifest <- tibble(
  cell = colnames(obj),
  condition = obj$condition_v564,
  sample = obj$sample_v564,
  lineage =
    as.character(
      obj@meta.data[[lineage_col]]
    )
)

if (
  !is.na(cluster_col)
) {
  cell_manifest$Res2_cluster <-
    as.character(
      obj@meta.data[[cluster_col]]
    )
}

write.csv(
  cell_manifest,
  file.path(
    DIR_TABLE,
    "01_exact_cells_used_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 12. Marker availability audit
# ==============================================================================

all_requested_genes <- unique(
  c(
    unlist(
      CONTAMINATION_CORE_MARKERS,
      use.names = FALSE
    ),
    unlist(
      CONTAMINATION_SUPPORTIVE_MARKERS,
      use.names = FALSE
    ),
    MACROPHAGE_IDENTITY_MARKERS,
    KUPFFER_RESIDENT_MARKERS,
    CANDIDATE_FACTORS
  )
)

gene_availability <- tibble(
  gene = all_requested_genes,
  present =
    all_requested_genes %in%
    rownames(obj)
)

write.csv(
  gene_availability,
  file.path(
    DIR_TABLE,
    "02_marker_availability_v5.6.4.csv"
  ),
  row.names = FALSE
)

CONTAMINATION_CORE_USE <- lapply(
  CONTAMINATION_CORE_MARKERS,
  function(x) {
    intersect(
      x,
      rownames(obj)
    )
  }
)

CONTAMINATION_SUPPORTIVE_USE <- lapply(
  CONTAMINATION_SUPPORTIVE_MARKERS,
  function(x) {
    intersect(
      x,
      rownames(obj)
    )
  }
)

MACROPHAGE_IDENTITY_USE <- intersect(
  MACROPHAGE_IDENTITY_MARKERS,
  rownames(obj)
)

KUPFFER_RESIDENT_USE <- intersect(
  KUPFFER_RESIDENT_MARKERS,
  rownames(obj)
)

CANDIDATE_FACTORS_USE <- intersect(
  CANDIDATE_FACTORS,
  rownames(obj)
)

msg(
  "Candidate factors available: ",
  paste(
    CANDIDATE_FACTORS_USE,
    collapse = ", "
  )
)

# ==============================================================================
# 13. RNA detection matrix
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

# ==============================================================================
# 14. Per-cell contamination marker counts
# ==============================================================================

msg("Calculating per-cell contamination marker counts...")

contam_count_df <- tibble(
  cell = colnames(obj)
)

for (
  lineage_name in names(
    CONTAMINATION_CORE_USE
  )
) {

  genes <- CONTAMINATION_CORE_USE[[
    lineage_name]
  ]

  if (
    length(genes) == 0
  ) {

    marker_count <- rep(
      0,
      ncol(obj)
    )

  } else {

    marker_count <- as.numeric(
      Matrix::colSums(
        counts_mat[
          genes,
          ,
          drop = FALSE
        ] > 0
      )
    )
  }

  contam_count_df[[
    paste0(
      lineage_name,
      "_marker_count"
    )]
  ] <- marker_count

  contam_count_df[[
    paste0(
      lineage_name,
      "_primary_flag"
    )]
  ] <- marker_count >=
    PRIMARY_MIN_MARKERS

  contam_count_df[[
    paste0(
      lineage_name,
      "_strong_flag"
    )]
  ] <- marker_count >=
    STRONG_MIN_MARKERS
}

primary_flag_cols <- grep(
  "_primary_flag$",
  colnames(
    contam_count_df
  ),
  value = TRUE
)

strong_flag_cols <- grep(
  "_strong_flag$",
  colnames(
    contam_count_df
  ),
  value = TRUE
)

contam_count_df$any_contamination_primary <-
  rowSums(
    contam_count_df[
      primary_flag_cols
    ]
  ) > 0

contam_count_df$any_contamination_strong <-
  rowSums(
    contam_count_df[
      strong_flag_cols
    ]
  ) > 0

contam_count_df$number_of_primary_lineages_flagged <-
  rowSums(
    contam_count_df[
      primary_flag_cols
    ]
  )

# ==============================================================================
# 15. Macrophage identity support per cell
# ==============================================================================

if (
  length(
    MACROPHAGE_IDENTITY_USE
  ) > 0
) {

  contam_count_df$macrophage_identity_marker_count <-
    as.numeric(
      Matrix::colSums(
        counts_mat[
          MACROPHAGE_IDENTITY_USE,
          ,
          drop = FALSE
        ] > 0
      )
    )

} else {

  contam_count_df$macrophage_identity_marker_count <-
    NA_real_
}

if (
  length(
    KUPFFER_RESIDENT_USE
  ) > 0
) {

  contam_count_df$Kupffer_resident_marker_count <-
    as.numeric(
      Matrix::colSums(
        counts_mat[
          KUPFFER_RESIDENT_USE,
          ,
          drop = FALSE
        ] > 0
      )
    )

} else {

  contam_count_df$Kupffer_resident_marker_count <-
    NA_real_
}

# ==============================================================================
# 16. Add QC flags to Seurat object
# ==============================================================================

msg("Attaching contamination QC vectors to Seurat metadata...")
#
# IMPORTANT v5.6.4 FIX
#   Do NOT assign row names to a tibble and do NOT pass a one-column tibble
#   into Seurat metadata.  Both can trigger unstable recursive behavior in
#   Seurat/tibble dispatch.  Match cells explicitly and assign plain vectors.
# ==============================================================================

obj_match_v564 <- match(
  colnames(obj),
  contam_count_df$cell
)

if (
  anyNA(
    obj_match_v564
  )
) {
  stop(
    "Failed to match contamination-QC table to Seurat cell names."
  )
}

if (
  !identical(
    contam_count_df$cell[
      obj_match_v564
    ],
    colnames(obj)
  )
) {
  stop(
    "Cell-order validation failed while attaching contamination metadata."
  )
}

obj$Monocyte_marker_count_v564 <-
  contam_count_df$Monocyte_marker_count[
    obj_match_v564
  ]

obj$Neutrophil_marker_count_v564 <-
  contam_count_df$Neutrophil_marker_count[
    obj_match_v564
  ]

obj$B_marker_count_v564 <-
  contam_count_df$B_marker_count[
    obj_match_v564
  ]

obj$T_marker_count_v564 <-
  contam_count_df$T_marker_count[
    obj_match_v564
  ]

obj$NK_marker_count_v564 <-
  contam_count_df$NK_marker_count[
    obj_match_v564
  ]

obj$contamination_primary_v564 <-
  factor(
    ifelse(
      contam_count_df$any_contamination_primary[
        obj_match_v564
      ],
      "Contamination_flagged",
      "Clean"
    ),
    levels = c(
      "Clean",
      "Contamination_flagged"
    )
  )

obj$contamination_strong_v564 <-
  factor(
    ifelse(
      contam_count_df$any_contamination_strong[
        obj_match_v564
      ],
      "Strong_flag",
      "No_strong_flag"
    ),
    levels = c(
      "No_strong_flag",
      "Strong_flag"
    )
  )

msg("Contamination metadata attached successfully.")

# ==============================================================================
# 17. Contamination summary
# ==============================================================================

contam_summary <- tibble(
  lineage = names(
    CONTAMINATION_CORE_USE
  )
) %>%
  rowwise() %>%
  mutate(
    n_markers_available =
      length(
        CONTAMINATION_CORE_USE[[
          lineage]
        ]
      ),
    markers_available =
      paste(
        CONTAMINATION_CORE_USE[[
          lineage]
        ],
        collapse = ";"
      ),
    n_primary_flagged =
      sum(
        contam_count_df[[
          paste0(
            lineage,
            "_primary_flag"
          )]
        ]
      ),
    pct_primary_flagged =
      100 *
      n_primary_flagged /
      ncol(obj),
    n_strong_flagged =
      sum(
        contam_count_df[[
          paste0(
            lineage,
            "_strong_flag"
          )]
        ]
      ),
    pct_strong_flagged =
      100 *
      n_strong_flagged /
      ncol(obj)
  ) %>%
  ungroup()

write.csv(
  contam_summary,
  file.path(
    DIR_TABLE,
    "03_contamination_summary_v5.6.4.csv"
  ),
  row.names = FALSE
)

overall_contam_summary <- tibble(
  metric = c(
    "total_cells",
    "primary_flagged_cells",
    "primary_flagged_pct",
    "strong_flagged_cells",
    "strong_flagged_pct"
  ),
  value = c(
    ncol(obj),
    sum(
      contam_count_df$any_contamination_primary
    ),
    100 *
      mean(
        contam_count_df$any_contamination_primary
      ),
    sum(
      contam_count_df$any_contamination_strong
    ),
    100 *
      mean(
        contam_count_df$any_contamination_strong
      )
  )
)

write.csv(
  overall_contam_summary,
  file.path(
    DIR_TABLE,
    "04_overall_contamination_summary_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 18. Contamination by sample
# ==============================================================================

contam_by_sample <- contam_count_df %>%
  mutate(
    sample =
      obj$sample_v564[
        match(
          cell,
          colnames(obj)
        )
      ],
    condition =
      obj$condition_v564[
        match(
          cell,
          colnames(obj)
        )
      ]
  ) %>%
  group_by(
    condition,
    sample
  ) %>%
  summarise(
    n_cells = n(),
    n_primary_flagged =
      sum(
        any_contamination_primary
      ),
    pct_primary_flagged =
      100 *
      mean(
        any_contamination_primary
      ),
    n_strong_flagged =
      sum(
        any_contamination_strong
      ),
    pct_strong_flagged =
      100 *
      mean(
        any_contamination_strong
      ),
    .groups = "drop"
  )

write.csv(
  contam_by_sample,
  file.path(
    DIR_TABLE,
    "05_contamination_by_sample_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 19. Candidate expression in clean vs contamination-flagged cells
# ==============================================================================

candidate_expr_list <- list()

for (
  gene in CANDIDATE_FACTORS_USE
) {

  expr <- as.numeric(
    data_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  detect <- as.numeric(
    counts_mat[
      gene,
      ,
      drop = TRUE
    ]
  ) > 0

  candidate_expr_list[[
    gene]
  ] <- tibble(
    gene = gene,

    clean_n =
      sum(
        !contam_count_df$any_contamination_primary
      ),

    flagged_n =
      sum(
        contam_count_df$any_contamination_primary
      ),

    clean_pct_positive =
      100 *
      mean(
        detect[
          !contam_count_df$any_contamination_primary
        ]
      ),

    flagged_pct_positive =
      if (
        any(
          contam_count_df$any_contamination_primary
        )
      ) {
        100 *
          mean(
            detect[
              contam_count_df$any_contamination_primary
            ]
          )
      } else {
        NA_real_
      },

    clean_mean_expression =
      mean(
        expr[
          !contam_count_df$any_contamination_primary
        ]
      ),

    flagged_mean_expression =
      if (
        any(
          contam_count_df$any_contamination_primary
        )
      ) {
        mean(
          expr[
            contam_count_df$any_contamination_primary
          ]
        )
      } else {
        NA_real_
      },

    fraction_candidate_positive_cells_that_are_flagged =
      if (
        sum(
          detect
        ) > 0
      ) {
        sum(
          detect &
            contam_count_df$any_contamination_primary
        ) /
          sum(
            detect
          )
      } else {
        NA_real_
      }
  )
}

candidate_clean_flagged <- bind_rows(
  candidate_expr_list
)

write.csv(
  candidate_clean_flagged,
  file.path(
    DIR_TABLE,
    "06_candidate_expression_clean_vs_flagged_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 20. Candidate expression by specific contamination lineage
# ==============================================================================

candidate_lineage_overlap <- list()

for (
  gene in CANDIDATE_FACTORS_USE
) {

  detect <- as.numeric(
    counts_mat[
      gene,
      ,
      drop = TRUE
    ]
  ) > 0

  for (
    lineage_name in names(
      CONTAMINATION_CORE_USE
    )
  ) {

    flag <- contam_count_df[[
      paste0(
        lineage_name,
        "_primary_flag"
      )]
    ]

    candidate_lineage_overlap[[
      paste(
        gene,
        lineage_name,
        sep = "__"
      )]
    ] <- tibble(
      gene = gene,
      contamination_lineage =
        lineage_name,
      candidate_positive_cells =
        sum(
          detect
        ),
      candidate_positive_and_flagged =
        sum(
          detect &
            flag
        ),
      fraction_candidate_positive_flagged =
        if (
          sum(
            detect
          ) > 0
        ) {
          sum(
            detect &
              flag
          ) /
            sum(
              detect
            )
        } else {
          NA_real_
        },
      pct_candidate_positive_within_flagged =
        if (
          sum(
            flag
          ) > 0
        ) {
          100 *
            mean(
              detect[
                flag
              ]
            )
        } else {
          NA_real_
        }
    )
  }
}

candidate_lineage_overlap <- bind_rows(
  candidate_lineage_overlap
)

write.csv(
  candidate_lineage_overlap,
  file.path(
    DIR_TABLE,
    "07_candidate_overlap_by_contamination_lineage_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 21. Export suspicious cells
# ==============================================================================

suspicious_cells <- contam_count_df %>%
  filter(
    any_contamination_primary
  ) %>%
  mutate(
    sample =
      obj$sample_v564[
        match(
          cell,
          colnames(obj)
        )
      ],
    condition =
      obj$condition_v564[
        match(
          cell,
          colnames(obj)
        )
      ]
  )

if (
  !is.na(cluster_col)
) {

  suspicious_cells$Res2_cluster <-
    as.character(
      obj@meta.data[
        suspicious_cells$cell,
        cluster_col
      ]
    )
}

write.csv(
  suspicious_cells,
  file.path(
    DIR_TABLE,
    "08_primary_contamination_flagged_cells_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 22. Candidate expression after excluding contamination-flagged cells
# ==============================================================================

clean_cells <- contam_count_df$cell[
  !contam_count_df$any_contamination_primary
]

clean_obj <- subset(
  obj,
  cells = clean_cells
)

candidate_by_sample_clean_list <- list()

clean_counts <- GetAssayData(
  clean_obj,
  assay = ASSAY_USE,
  layer = "counts"
)

clean_data <- GetAssayData(
  clean_obj,
  assay = ASSAY_USE,
  layer = "data"
)

for (
  gene in CANDIDATE_FACTORS_USE
) {

  for (
    s in unique(
      clean_obj$sample_v564
    )
  ) {

    idx <- which(
      clean_obj$sample_v564 == s
    )

    candidate_by_sample_clean_list[[
      paste(
        gene,
        s,
        sep = "__"
      )]
    ] <- tibble(
      gene = gene,
      sample = s,
      condition =
        unique(
          clean_obj$condition_v564[
            idx
          ]
        )[[1]],
      n_cells = length(
        idx
      ),
      pct_positive =
        100 *
        mean(
          as.numeric(
            clean_counts[
              gene,
              idx,
              drop = TRUE
            ]
          ) > 0
        ),
      mean_expression =
        mean(
          as.numeric(
            clean_data[
              gene,
              idx,
              drop = TRUE
            ]
          )
        )
    )
  }
}

candidate_by_sample_clean <- bind_rows(
  candidate_by_sample_clean_list
)

write.csv(
  candidate_by_sample_clean,
  file.path(
    DIR_TABLE,
    "09_candidate_expression_by_sample_after_QC_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 23. DotPlot: contamination markers by sample
# ==============================================================================

core_features_use <- unique(
  unlist(
    CONTAMINATION_CORE_USE,
    use.names = FALSE
  )
)

supportive_features_use <- unique(
  unlist(
    CONTAMINATION_SUPPORTIVE_USE,
    use.names = FALSE
  )
)

obj$sample_v564 <- factor(
  obj$sample_v564,
  levels = EXPECTED_SAMPLES
)

if (
  length(
    core_features_use
  ) > 0
) {

  p_contam_core <- DotPlot(
    obj,
    features = core_features_use,
    group.by = "sample_v564",
    assay = ASSAY_USE,
    dot.scale = 8
  ) +
    scale_colour_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title =
        "STRICT contamination CORE markers in v5.6.0-used Kupffer cells",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        size = 13,
        face = "bold"
      )
    )

  save_pdf(
    p_contam_core,
    file.path(
      DIR_DOTPLOT,
      "01_STRICT_contamination_CORE_markers_by_sample_v5.6.4.pdf"
    ),
    15,
    6
  )
}

if (
  length(
    supportive_features_use
  ) > 0
) {

  p_contam_support <- DotPlot(
    obj,
    features = supportive_features_use,
    group.by = "sample_v564",
    assay = ASSAY_USE,
    dot.scale = 8
  ) +
    scale_colour_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title =
        "Supportive/non-specific markers (NOT used for contamination flags)",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        size = 13,
        face = "bold"
      )
    )

  save_pdf(
    p_contam_support,
    file.path(
      DIR_DOTPLOT,
      "02_supportive_contamination_markers_by_sample_v5.6.4.pdf"
    ),
    10,
    5
  )
}

# ==============================================================================
# 24. DotPlot: macrophage identity + candidate factors
# ==============================================================================

identity_candidate_list <- list(
  Macrophage_identity =
    MACROPHAGE_IDENTITY_USE,
  Kupffer_resident =
    KUPFFER_RESIDENT_USE,
  Candidate_factors =
    CANDIDATE_FACTORS_USE
)

identity_candidate_list <- identity_candidate_list[
  lengths(
    identity_candidate_list
  ) > 0
]

# v5.6.4:
# Aif1 occurs in both Macrophage_identity and Candidate_factors.
# Marco occurs in both Kupffer_resident and Candidate_factors.
# Flatten and deduplicate before DotPlot to avoid duplicated factor levels.
identity_candidate_features_unique <- unique(
  c(
    MACROPHAGE_IDENTITY_USE,
    KUPFFER_RESIDENT_USE,
    CANDIDATE_FACTORS_USE
  )
)

p_identity_dot <- DotPlot(
  obj,
  features = identity_candidate_features_unique,
  group.by = "sample_v564",
  assay = ASSAY_USE,
  dot.scale = 8
) +
  scale_colour_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "Macrophage identity and candidate factors in v5.6.0-used cells",
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_identity_dot,
  file.path(
    DIR_DOTPLOT,
    "03_macrophage_identity_and_candidates_by_sample_v5.6.4.pdf"
  ),
  14,
  6
)

# Candidate-only DotPlot: keeps Marco / Aif1 / Siglece directly visible.
if (
  length(
    CANDIDATE_FACTORS_USE
  ) > 0
) {

  candidate_features_unique <- unique(
    CANDIDATE_FACTORS_USE
  )

  p_candidate_dot <- DotPlot(
    obj,
    features = candidate_features_unique,
    group.by = "sample_v564",
    assay = ASSAY_USE,
    dot.scale = 9
  ) +
    scale_colour_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    labs(
      title =
        "Candidate factors in v5.6.0-used Kupffer cells",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        size = 13,
        face = "bold"
      )
    )

  save_pdf(
    p_candidate_dot,
    file.path(
      DIR_DOTPLOT,
      "04_candidate_factors_by_sample_DotPlot_v5.6.4.pdf"
    ),
    10,
    5
  )
}

# ==============================================================================
# 25. UMAP: contamination flags
# ==============================================================================

if (
  !is.na(
    umap_name
  )
) {

  p_flag <- DimPlot(
    obj,
    reduction = umap_name,
    group.by = "contamination_primary_v564",
    pt.size = 0.45,
    raster = FALSE
  ) +
    ggtitle(
      "Primary contamination QC flags in v5.6.0-used Kupffer cells"
    ) +
    theme_classic(
      base_size = 10
    )

  save_pdf(
    p_flag,
    file.path(
      DIR_UMAP,
      "01_primary_contamination_flags_UMAP_v5.6.4.pdf"
    ),
    7,
    6
  )

  p_flag_split <- DimPlot(
    obj,
    reduction = umap_name,
    group.by = "contamination_primary_v564",
    split.by = "sample_v564",
    ncol = 4,
    pt.size = 0.35,
    raster = FALSE
  ) +
    theme_classic(
      base_size = 8
    )

  save_pdf(
    p_flag_split,
    file.path(
      DIR_UMAP,
      "02_primary_contamination_flags_by_sample_UMAP_v5.6.4.pdf"
    ),
    15,
    4.5
  )
}

# ==============================================================================
# 26. UMAP: each contamination lineage
# ==============================================================================

if (
  !is.na(
    umap_name
  )
) {

  for (
    lineage_name in names(
      CONTAMINATION_CORE_USE
    )
  ) {

    flag_col <- paste0(
      lineage_name,
      "_primary_flag"
    )

    lineage_flag_v564 <- contam_count_df[[flag_col]][
      obj_match_v564
    ]

    obj[[
      paste0(
        lineage_name,
        "_QC_v564"
      )
    ]] <- factor(
      ifelse(
        lineage_flag_v564,
        paste0(
          lineage_name,
          "_flag"
        ),
        "No_flag"
      ),
      levels = c(
        "No_flag",
        paste0(
          lineage_name,
          "_flag"
        )
      )
    )

    p <- DimPlot(
      obj,
      reduction = umap_name,
      group.by = paste0(
        lineage_name,
        "_QC_v564"
      ),
      pt.size = 0.45,
      raster = FALSE
    ) +
      ggtitle(
        paste0(
          lineage_name,
          " contamination-marker flag"
        )
      ) +
      theme_classic(
        base_size = 10
      )

    save_pdf(
      p,
      file.path(
        DIR_UMAP,
        paste0(
          "QC_",
          lineage_name,
          "_UMAP_v5.6.4.pdf"
        )
      ),
      7,
      6
    )
  }
}

# ==============================================================================
# 27. Candidate FeaturePlots
# ==============================================================================

if (
  !is.na(
    umap_name
  ) &&
  length(
    CANDIDATE_FACTORS_USE
  ) > 0
) {

  p_candidates <- FeaturePlot(
    obj,
    features = CANDIDATE_FACTORS_USE,
    reduction = umap_name,
    ncol = 4,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
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
    p_candidates,
    file.path(
      DIR_UMAP,
      "03_candidate_factors_FeaturePlot_v5.6.4.pdf"
    ),
    15,
    8
  )

  for (
    gene in CANDIDATE_FACTORS_USE
  ) {

    p_gene <- FeaturePlot(
      obj,
      features = gene,
      reduction = umap_name,
      split.by = "sample_v564",
      keep.scale = "all",
      ncol = 4,
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95",
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
      p_gene,
      file.path(
        DIR_UMAP,
        paste0(
          "Candidate_",
          gene,
          "_by_sample_FeaturePlot_v5.6.4.pdf"
        )
      ),
      15,
      4.5
    )
  }
}

# ==============================================================================
# 28. Heatmap of per-cell contamination-marker counts by sample summary
# ==============================================================================

count_long <- contam_count_df %>%
  select(
    cell,
    ends_with(
      "_marker_count"
    )
  ) %>%
  pivot_longer(
    cols = -cell,
    names_to = "lineage",
    values_to = "marker_count"
  ) %>%
  mutate(
    lineage =
      sub(
        "_marker_count$",
        "",
        lineage
      ),
    sample =
      obj$sample_v564[
        match(
          cell,
          colnames(obj)
        )
      ]
  ) %>%
  group_by(
    sample,
    lineage
  ) %>%
  summarise(
    mean_marker_count =
      mean(
        marker_count
      ),
    pct_cells_with_1plus =
      100 *
      mean(
        marker_count >= 1
      ),
    pct_cells_with_2plus =
      100 *
      mean(
        marker_count >= PRIMARY_MIN_MARKERS
      ),
    pct_cells_with_3plus =
      100 *
      mean(
        marker_count >= STRONG_MIN_MARKERS
      ),
    .groups = "drop"
  )

write.csv(
  count_long,
  file.path(
    DIR_TABLE,
    "10_contamination_marker_burden_by_sample_v5.6.4.csv"
  ),
  row.names = FALSE
)

count_long$sample <- factor(
  count_long$sample,
  levels = EXPECTED_SAMPLES
)

p_burden <- ggplot(
  count_long,
  aes(
    x = sample,
    y = lineage,
    fill = pct_cells_with_2plus
  )
) +
  geom_tile(
    linewidth = 0.3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  labs(
    title =
      "Cells with >=2 contamination markers",
    x = NULL,
    y = NULL,
    fill = "% cells"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold"
    )
  )

save_pdf(
  p_burden,
  file.path(
    DIR_HEATMAP,
    "01_contamination_burden_by_sample_v5.6.4.pdf"
  ),
  7,
  5
)

# ==============================================================================
# 28B. Candidate expression before vs after STRICT contamination QC
# ==============================================================================

strict_clean_mask <- !contam_count_df$any_contamination_primary

candidate_prepost_list <- list()

for (
  gene in CANDIDATE_FACTORS_USE
) {

  detect_all <- as.numeric(
    counts_mat[
      gene,
      ,
      drop = TRUE
    ]
  ) > 0

  expr_all <- as.numeric(
    data_mat[
      gene,
      ,
      drop = TRUE
    ]
  )

  all_pct <- 100 *
    mean(
      detect_all
    )

  clean_pct <- 100 *
    mean(
      detect_all[
        strict_clean_mask
      ]
    )

  all_mean <- mean(
    expr_all
  )

  clean_mean <- mean(
    expr_all[
      strict_clean_mask
    ]
  )

  candidate_prepost_list[[gene]] <- tibble(
    gene = gene,
    all_cells_n =
      length(
        detect_all
      ),
    strict_clean_n =
      sum(
        strict_clean_mask
      ),
    all_cells_pct_positive =
      all_pct,
    strict_clean_pct_positive =
      clean_pct,
    pct_positive_change_after_QC =
      clean_pct -
      all_pct,
    all_cells_mean_expression =
      all_mean,
    strict_clean_mean_expression =
      clean_mean,
    mean_expression_ratio_after_QC =
      ifelse(
        all_mean > 0,
        clean_mean /
          all_mean,
        NA_real_
      )
  )
}

candidate_prepost <- bind_rows(
  candidate_prepost_list
)

write.csv(
  candidate_prepost,
  file.path(
    DIR_TABLE,
    "11_candidate_pre_vs_post_STRICT_QC_v5.6.4.csv"
  ),
  row.names = FALSE
)

candidate_sample_prepost_list <- list()

for (
  gene in CANDIDATE_FACTORS_USE
) {

  for (
    s_id in EXPECTED_SAMPLES
  ) {

    idx_all <- which(
      as.character(
        obj$sample_v564
      ) ==
        s_id
    )

    if (
      !length(
        idx_all
      )
    ) {
      next
    }

    idx_clean <- idx_all[
      strict_clean_mask[
        idx_all
      ]
    ]

    detect_before <- as.numeric(
      counts_mat[
        gene,
        idx_all,
        drop = TRUE
      ]
    ) > 0

    expr_before <- as.numeric(
      data_mat[
        gene,
        idx_all,
        drop = TRUE
      ]
    )

    if (
      length(
        idx_clean
      ) > 0
    ) {

      detect_after <- as.numeric(
        counts_mat[
          gene,
          idx_clean,
          drop = TRUE
        ]
      ) > 0

      expr_after <- as.numeric(
        data_mat[
          gene,
          idx_clean,
          drop = TRUE
        ]
      )

      pct_after <- 100 *
        mean(
          detect_after
        )

      mean_after <- mean(
        expr_after
      )

    } else {

      pct_after <- NA_real_
      mean_after <- NA_real_
    }

    candidate_sample_prepost_list[[
      paste(
        gene,
        s_id,
        sep = "__"
      )]
    ] <- tibble(
      gene = gene,
      sample = s_id,
      condition =
        unique(
          obj$condition_v564[
            idx_all
          ]
        )[[1]],
      n_before_QC =
        length(
          idx_all
        ),
      n_after_QC =
        length(
          idx_clean
        ),
      pct_positive_before_QC =
        100 *
        mean(
          detect_before
        ),
      pct_positive_after_QC =
        pct_after,
      mean_expression_before_QC =
        mean(
          expr_before
        ),
      mean_expression_after_QC =
        mean_after
    )
  }
}

candidate_sample_prepost <- bind_rows(
  candidate_sample_prepost_list
)

write.csv(
  candidate_sample_prepost,
  file.path(
    DIR_TABLE,
    "12_candidate_by_sample_pre_vs_post_STRICT_QC_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 29. Summary verdict table for candidate robustness
# ==============================================================================

candidate_verdict <- candidate_clean_flagged %>%
  mutate(
    likely_contamination_driven =
      case_when(
        is.na(
          fraction_candidate_positive_cells_that_are_flagged
        ) ~ NA,
        fraction_candidate_positive_cells_that_are_flagged >=
          0.50 ~ TRUE,
        TRUE ~ FALSE
      ),
    retained_in_clean_Kupffer_cells =
      clean_pct_positive > 0
  ) %>%
  arrange(
    fraction_candidate_positive_cells_that_are_flagged
  )

write.csv(
  candidate_verdict,
  file.path(
    DIR_TABLE,
    "13_candidate_contamination_robustness_summary_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 30. Analysis metadata
# ==============================================================================

analysis_metadata <- tibble(
  parameter = c(
    "script_version",
    "input_RDS",
    "reference_v5.6.0_sample_counts",
    "assay",
    "condition_source",
    "sample_column",
    "lineage_column",
    "lineage_value",
    "target_groups",
    "primary_min_contamination_markers",
    "strong_min_contamination_markers",
    "candidate_factors",
    "contamination_flag_basis",
    "supportive_markers_used_for_flag"
  ),
  value = c(
    "v5.6.4",
    MPHI_RDS,
    V560_SAMPLE_COUNTS,
    ASSAY_USE,
    condition_info$source,
    sample_col,
    lineage_col,
    MACROPHAGE_LINEAGE_VALUE,
    paste(
      TARGET_GROUPS,
      collapse = ";"
    ),
    PRIMARY_MIN_MARKERS,
    STRONG_MIN_MARKERS,
    paste(
      CANDIDATE_FACTORS_USE,
      collapse = ";"
    ),
    "lineage-specific CORE markers only",
    "FALSE"
  )
)

write.csv(
  analysis_metadata,
  file.path(
    DIR_LOG,
    "analysis_metadata_v5.6.4.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_LOG,
    "sessionInfo_v5.6.4.txt"
  )
)

# ==============================================================================
# 31. Output index
# ==============================================================================

index <- tibble(
  item = c(
    "Reproduced sample counts",
    "v5.6.0 vs v5.6.4 count check",
    "Exact cells used",
    "Marker availability",
    "Strict contamination summary",
    "Overall strict contamination summary",
    "Strict contamination by sample",
    "Candidate clean-vs-flagged expression",
    "Candidate overlap by strict contamination lineage",
    "Strict contamination-flagged cells",
    "Candidate expression by sample after strict QC",
    "Strict contamination-marker burden by sample",
    "Candidate pre-vs-post strict QC",
    "Candidate by sample pre-vs-post strict QC",
    "Candidate contamination robustness summary",
    "Strict core contamination DotPlot",
    "Supportive marker DotPlot",
    "Macrophage identity/candidate DotPlot",
    "Candidate-only DotPlot",
    "Primary strict contamination UMAP",
    "Candidate FeaturePlots"
  ),
  path = c(
    file.path(DIR_TABLE, "00_sample_cell_counts_reproduced_v5.6.4.csv"),
    file.path(DIR_TABLE, "00_v5.6.0_vs_v5.6.4_cell_count_check.csv"),
    file.path(DIR_TABLE, "01_exact_cells_used_v5.6.4.csv"),
    file.path(DIR_TABLE, "02_marker_availability_v5.6.4.csv"),
    file.path(DIR_TABLE, "03_contamination_summary_v5.6.4.csv"),
    file.path(DIR_TABLE, "04_overall_contamination_summary_v5.6.4.csv"),
    file.path(DIR_TABLE, "05_contamination_by_sample_v5.6.4.csv"),
    file.path(DIR_TABLE, "06_candidate_expression_clean_vs_flagged_v5.6.4.csv"),
    file.path(DIR_TABLE, "07_candidate_overlap_by_contamination_lineage_v5.6.4.csv"),
    file.path(DIR_TABLE, "08_primary_contamination_flagged_cells_v5.6.4.csv"),
    file.path(DIR_TABLE, "09_candidate_expression_by_sample_after_QC_v5.6.4.csv"),
    file.path(DIR_TABLE, "10_contamination_marker_burden_by_sample_v5.6.4.csv"),
    file.path(DIR_TABLE, "11_candidate_pre_vs_post_STRICT_QC_v5.6.4.csv"),
    file.path(DIR_TABLE, "12_candidate_by_sample_pre_vs_post_STRICT_QC_v5.6.4.csv"),
    file.path(DIR_TABLE, "13_candidate_contamination_robustness_summary_v5.6.4.csv"),
    file.path(DIR_DOTPLOT, "01_STRICT_contamination_CORE_markers_by_sample_v5.6.4.pdf"),
    file.path(DIR_DOTPLOT, "02_supportive_contamination_markers_by_sample_v5.6.4.pdf"),
    file.path(DIR_DOTPLOT, "03_macrophage_identity_and_candidates_by_sample_v5.6.4.pdf"),
    file.path(DIR_DOTPLOT, "04_candidate_factors_by_sample_DotPlot_v5.6.4.pdf"),
    file.path(DIR_UMAP, "01_primary_contamination_flags_UMAP_v5.6.4.pdf"),
    file.path(DIR_UMAP, "03_candidate_factors_FeaturePlot_v5.6.4.pdf")
  )
)

write.csv(
  index,
  file.path(
    OUTPUT_DIR,
    "Kupffer_contamination_reaudit_INDEX_v5.6.4.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 32. Final terminal report
# ==============================================================================

msg(
  "DONE."
)

msg(
  "Output: ",
  OUTPUT_DIR
)

cat(
  "\n============================================================\n"
)

cat(
  "v5.6.0-used Kupffer_Macrophage STRICT contamination re-audit\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Cells analyzed: ",
  ncol(obj),
  "\n",
  sep = ""
)

cat(
  "Primary contamination-flagged cells: ",
  sum(
    contam_count_df$any_contamination_primary
  ),
  " (",
  round(
    100 *
      mean(
        contam_count_df$any_contamination_primary
      ),
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "Strong contamination-flagged cells: ",
  sum(
    contam_count_df$any_contamination_strong
  ),
  " (",
  round(
    100 *
      mean(
        contam_count_df$any_contamination_strong
      ),
    2
  ),
  "%)\n\n",
  sep = ""
)

cat(
  "Contamination by lineage:\n"
)

print(
  contam_summary
)

cat(
  "\nCandidate robustness:\n"
)

print(
  candidate_verdict %>%
    select(
      gene,
      clean_pct_positive,
      flagged_pct_positive,
      clean_mean_expression,
      flagged_mean_expression,
      fraction_candidate_positive_cells_that_are_flagged,
      likely_contamination_driven,
      retained_in_clean_Kupffer_cells
    )
)

cat(
  "\n============================================================\n"
)
