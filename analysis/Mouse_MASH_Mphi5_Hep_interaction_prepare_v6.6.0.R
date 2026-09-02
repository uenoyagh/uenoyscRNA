#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6600)

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
# 5 macrophage subtypes x 5 Hepatocyte states
# Interaction preparation
#
# Version: v6.6.0
#
# PURPOSE
# -------
# Build one frozen interaction-ready Seurat object for sample-wise
# macrophage <-> hepatocyte communication analyses (CellChat / NicheNet).
#
# IMPORTANT DESIGN
# ----------------
# 1) No re-clustering is performed in this script.
# 2) Final macrophage identities are transferred from the frozen Clean-B
#    macrophage object by EXACT cell-barcode match.
# 3) Final hepatocyte identities are transferred from the v6.5.2 final
#    hepatocyte object by EXACT cell-barcode match.
# 4) All transferred cells are anchored back to the frozen whole-cell parent.
# 5) Only Sham1, Sham20, Tx17, Tx5 are included.
# 6) MT_high_QC_Hepatocyte is retained in an exclusion audit but is NOT used
#    in the primary interaction-ready object.
#
# MACROPHAGE INTERACTION STATES (5)
# ---------------------------------
#   Mphi_Anti-inflammatory
#   Mphi_Inflammatory
#   Mphi_ECM-associated-inflammatory
#   Mphi_Repair-Resolution
#   Mphi_Lipid-associated-TREM2
#
# HEPATOCYTE INTERACTION STATES (5)
# ---------------------------------
# v6.5.2 detailed states are collapsed only for communication analysis:
#
#   Periportal_Hepatocyte_1
#   Periportal_Hepatocyte_2
#       -> Hep_Periportal
#
#   Pericentral_Hepatocyte
#       -> Hep_Pericentral
#
#   Injury_inflammatory_Hepatocyte
#       -> Hep_Injury-inflammatory
#
#   Intermediate_Hepatocyte
#       -> Hep_Intermediate
#
#   Cycling_G2M_Hepatocyte
#   Cycling_S_Hepatocyte
#       -> Hep_Cycling
#
#   MT_high_QC_Hepatocyte
#       -> EXCLUDED from primary interaction analysis
#
# RATIONALE
# ---------
# The detailed hepatocyte labels remain stored in metadata. The merged 5-state
# interaction annotation reduces over-fragmentation while preserving the
# biological axes needed to test macrophage <-> hepatocyte communication.
#
# PRIMARY OUTPUT
# --------------
# metadata:
#   interaction_celltype_v660
#   interaction_lineage_v660
#   sample_interaction_v660
#   condition_interaction_v660
#   mphi_state_detailed_v660
#   hepatocyte_state_detailed_v660
#
# RDS:
#   Mouse_MASH_Mphi5_Hep5_interaction_ready_v6.6.0.rds
# ==============================================================================


# ==============================================================================
# 1. Helper functions
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste0(...)
  )
}


save_pdf <- function(
  plot,
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

  print(
    plot
  )

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


canonical_condition <- function(sample_name) {

  x <- as.character(
    sample_name
  )

  case_when(
    grepl(
      "^Sham",
      x,
      ignore.case = TRUE
    ) ~
      "Sham",

    grepl(
      "^Tx",
      x,
      ignore.case = TRUE
    ) ~
      "Tx",

    grepl(
      "^STD",
      x,
      ignore.case = TRUE
    ) ~
      "STD",

    grepl(
      "CDHFD|CDAHFD",
      x,
      ignore.case = TRUE
    ) ~
      "CDAHFD",

    TRUE ~
      NA_character_
  )
}


find_sample_column <- function(
  object,
  target_samples,
  preferred = character()
) {

  meta <- object@meta.data

  candidates <- unique(
    c(
      preferred,
      "sample_for_annotation",
      "sample_interaction_v620",
      "sample_interaction_v660",
      "sample_hep_v650",
      "sample",
      "Sample",
      "orig.ident",
      "sample_id",
      "SampleID",
      "sample_name"
    )
  )

  candidates <- candidates[
    candidates %in%
      colnames(
        meta
      )
  ]

  score_candidate <- function(col) {

    x <- as.character(
      meta[[
        col
      ]]
    )

    present <- target_samples[
      target_samples %in%
        unique(
          x
        )
    ]

    tibble(
      column =
        col,
      n_target_samples_found =
        length(
          present
        ),
      target_samples_found =
        paste(
          present,
          collapse = ","
        )
    )
  }

  audit <- bind_rows(
    lapply(
      candidates,
      score_candidate
    )
  )

  # Wider search among character / factor metadata if preferred candidates fail.
  if (
    !nrow(
      audit
    ) ||
    max(
      audit$n_target_samples_found
    ) <
      length(
        target_samples
      )
  ) {

    wider <- colnames(
      meta
    )[
      vapply(
        meta,
        function(x) {
          is.character(
            x
          ) ||
            is.factor(
              x
            )
        },
        logical(
          1
        )
      )
    ]

    wider <- setdiff(
      wider,
      candidates
    )

    if (
      length(
        wider
      )
    ) {

      audit <- bind_rows(
        audit,
        bind_rows(
          lapply(
            wider,
            score_candidate
          )
        )
      )
    }
  }

  if (
    !nrow(
      audit
    )
  ) {
    stop(
      "No candidate sample metadata columns found."
    )
  }

  audit <- audit %>%
    arrange(
      desc(
        n_target_samples_found
      ),
      column
    )

  best <- audit$column[
    1
  ]

  if (
    audit$n_target_samples_found[
      1
    ] <
      length(
        target_samples
      )
  ) {

    stop(
      "Could not identify a metadata column containing all target samples. ",
      "Best column: ",
      best,
      " ; found ",
      audit$n_target_samples_found[
        1
      ],
      "/",
      length(
        target_samples
      ),
      " target samples."
    )
  }

  list(
    column =
      best,
    audit =
      audit
  )
}


exact_barcode_audit <- function(
  source_cells,
  parent_cells,
  label
) {

  source_cells <- unique(
    as.character(
      source_cells
    )
  )

  parent_cells <- unique(
    as.character(
      parent_cells
    )
  )

  matched <- source_cells[
    source_cells %in%
      parent_cells
  ]

  unmatched <- setdiff(
    source_cells,
    parent_cells
  )

  rate <- if (
    length(
      source_cells
    ) >
      0
  ) {
    length(
      matched
    ) /
      length(
        source_cells
      )
  } else {
    NA_real_
  }

  summary <- tibble(
    source =
      label,
    n_source_cells =
      length(
        source_cells
      ),
    n_exact_matched =
      length(
        matched
      ),
    n_unmatched =
      length(
        unmatched
      ),
    exact_match_rate =
      rate
  )

  list(
    summary =
      summary,
    matched =
      matched,
    unmatched =
      unmatched
  )
}


fraction_table <- function(
  meta,
  sample_col,
  group_col,
  target_samples
) {

  meta %>%
    as_tibble(
      rownames =
        "cell"
    ) %>%
    transmute(
      cell,
      sample =
        as.character(
          .data[[
            sample_col
          ]]
        ),
      group =
        as.character(
          .data[[
            group_col
          ]]
        )
    ) %>%
    filter(
      sample %in%
        target_samples
    ) %>%
    count(
      sample,
      group,
      name =
        "n_cells"
    ) %>%
    complete(
      sample =
        target_samples,
      group,
      fill =
        list(
          n_cells =
            0
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


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"


WHOLECELL_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)


MPHI_RDS_CANDIDATES <- c(
  file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
  ),
  file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS"
  )
)


HEP_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Hepatocyte",
  "Hepatocyte_FINAL_state_pseudobulk_v6.5.2",
  "RDS",
  "Mouse_MASH_Hepatocyte_FINAL_annotated_v6.5.2.rds"
)


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_interaction_ready_v6.6.0"
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
# 3. Fixed metadata columns and target samples
# ==============================================================================

TARGET_SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)


WHOLECELL_LAYER1_COL <-
  "wholecell_layer1_FINAL_v511"


MPHI_ANNOTATION_COL <-
  "macrophage_class_Res2_FINAL_v4145_char"


MPHI_SAMPLE_COL_PREFERRED <-
  "sample_for_annotation"


HEP_STATE_COL <-
  "hepatocyte_state_FINAL_v652"


HEP_SAMPLE_COL_PREFERRED <-
  "sample_hep_v650"


INTERACTION_CELLTYPE_COL <-
  "interaction_celltype_v660"


INTERACTION_LINEAGE_COL <-
  "interaction_lineage_v660"


INTERACTION_SAMPLE_COL <-
  "sample_interaction_v660"


INTERACTION_CONDITION_COL <-
  "condition_interaction_v660"


MPHI_DETAIL_COL <-
  "mphi_state_detailed_v660"


HEP_DETAIL_COL <-
  "hepatocyte_state_detailed_v660"


HEP_INTERACTION_STATE_COL <-
  "hepatocyte_state_interaction_v660"


SOURCE_COL <-
  "interaction_label_source_v660"


# ==============================================================================
# 4. Fixed state maps
# ==============================================================================

MPHI_FINAL_STATES <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)


MPHI_INTERACTION_MAP <- c(
  "Anti-inflammatory-Mphi" =
    "Mphi_Anti-inflammatory",

  "Inflammatory-Mphi" =
    "Mphi_Inflammatory",

  "ECM-associated inflammatory-Mphi" =
    "Mphi_ECM-associated-inflammatory",

  "Repair/Resolution-Mphi" =
    "Mphi_Repair-Resolution",

  "Lipid-associated/TREM2-Mphi" =
    "Mphi_Lipid-associated-TREM2"
)


HEP_DETAIL_STATES_EXPECTED <- c(
  "Periportal_Hepatocyte_1",
  "Injury_inflammatory_Hepatocyte",
  "Pericentral_Hepatocyte",
  "MT_high_QC_Hepatocyte",
  "Periportal_Hepatocyte_2",
  "Intermediate_Hepatocyte",
  "Cycling_G2M_Hepatocyte",
  "Cycling_S_Hepatocyte"
)


HEP_INTERACTION_MAP <- c(
  "Periportal_Hepatocyte_1" =
    "Hep_Periportal",

  "Periportal_Hepatocyte_2" =
    "Hep_Periportal",

  "Pericentral_Hepatocyte" =
    "Hep_Pericentral",

  "Injury_inflammatory_Hepatocyte" =
    "Hep_Injury-inflammatory",

  "Intermediate_Hepatocyte" =
    "Hep_Intermediate",

  "Cycling_G2M_Hepatocyte" =
    "Hep_Cycling",

  "Cycling_S_Hepatocyte" =
    "Hep_Cycling",

  "MT_high_QC_Hepatocyte" =
    NA_character_
)


INTERACTION_LEVELS <- c(
  "Mphi_Anti-inflammatory",
  "Mphi_Inflammatory",
  "Mphi_ECM-associated-inflammatory",
  "Mphi_Repair-Resolution",
  "Mphi_Lipid-associated-TREM2",
  "Hep_Periportal",
  "Hep_Pericentral",
  "Hep_Injury-inflammatory",
  "Hep_Intermediate",
  "Hep_Cycling"
)


LINEAGE_LEVELS <- c(
  "Macrophage",
  "Hepatocyte"
)


MIN_EXACT_MATCH_RATE <-
  0.95


MIN_CELLS_STOP <-
  10


MIN_CELLS_WARN <-
  30


# ==============================================================================
# 5. Load inputs
# ==============================================================================

if (
  !file.exists(
    WHOLECELL_RDS
  )
) {

  stop(
    "Whole-cell frozen parent RDS missing: ",
    WHOLECELL_RDS
  )
}


MPHI_RDS <- MPHI_RDS_CANDIDATES[
  file.exists(
    MPHI_RDS_CANDIDATES
  )
][
  1
]


if (
  is.na(
    MPHI_RDS
  )
) {

  stop(
    "Macrophage Clean-B RDS not found. Candidates:\n",
    paste(
      MPHI_RDS_CANDIDATES,
      collapse = "\n"
    )
  )
}


if (
  !file.exists(
    HEP_RDS
  )
) {

  stop(
    "Final Hepatocyte v6.5.2 RDS missing: ",
    HEP_RDS
  )
}


msg(
  "Loading frozen whole-cell parent..."
)


parent <- readRDS(
  WHOLECELL_RDS
)


msg(
  "Loading frozen Clean-B macrophage object..."
)


mphi <- readRDS(
  MPHI_RDS
)


msg(
  "Loading final Hepatocyte v6.5.2 object..."
)


hep <- readRDS(
  HEP_RDS
)


DefaultAssay(
  parent
) <- "RNA"


DefaultAssay(
  mphi
) <- "RNA"


DefaultAssay(
  hep
) <- "RNA"


# ==============================================================================
# 6. Validate fixed source annotation columns
# ==============================================================================

if (
  !WHOLECELL_LAYER1_COL %in%
    colnames(
      parent@meta.data
    )
) {

  stop(
    "Frozen parent Layer1 column missing: ",
    WHOLECELL_LAYER1_COL
  )
}


if (
  !MPHI_ANNOTATION_COL %in%
    colnames(
      mphi@meta.data
    )
) {

  stop(
    "Macrophage annotation column missing: ",
    MPHI_ANNOTATION_COL
  )
}


if (
  !HEP_STATE_COL %in%
    colnames(
      hep@meta.data
    )
) {

  stop(
    "Hepatocyte final state column missing: ",
    HEP_STATE_COL
  )
}


# ==============================================================================
# 7. Detect sample columns
# ==============================================================================

parent_sample_detect <- find_sample_column(
  parent,
  TARGET_SAMPLES,
  preferred = c(
    "sample_for_annotation"
  )
)


mphi_sample_detect <- find_sample_column(
  mphi,
  TARGET_SAMPLES,
  preferred = c(
    MPHI_SAMPLE_COL_PREFERRED
  )
)


hep_sample_detect <- find_sample_column(
  hep,
  TARGET_SAMPLES,
  preferred = c(
    HEP_SAMPLE_COL_PREFERRED
  )
)


PARENT_SAMPLE_COL <-
  parent_sample_detect$column


MPHI_SAMPLE_COL <-
  mphi_sample_detect$column


HEP_SAMPLE_COL <-
  hep_sample_detect$column


msg(
  "Detected parent sample column: ",
  PARENT_SAMPLE_COL
)


msg(
  "Detected Mphi sample column: ",
  MPHI_SAMPLE_COL
)


msg(
  "Detected Hepatocyte sample column: ",
  HEP_SAMPLE_COL
)


write.csv(
  parent_sample_detect$audit,
  file.path(
    TAB_OUT,
    "01_parent_sample_column_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


write.csv(
  mphi_sample_detect$audit,
  file.path(
    TAB_OUT,
    "02_Mphi_sample_column_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


write.csv(
  hep_sample_detect$audit,
  file.path(
    TAB_OUT,
    "03_Hepatocyte_sample_column_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Restrict source objects to the 4 target biological samples
# ==============================================================================

mphi_target_cells <- rownames(
  mphi@meta.data
)[
  as.character(
    mphi@meta.data[[
      MPHI_SAMPLE_COL
    ]]
  ) %in%
    TARGET_SAMPLES
]


hep_target_cells <- rownames(
  hep@meta.data
)[
  as.character(
    hep@meta.data[[
      HEP_SAMPLE_COL
    ]]
  ) %in%
    TARGET_SAMPLES
]


if (
  !length(
    mphi_target_cells
  )
) {

  stop(
    "No target Mphi cells found."
  )
}


if (
  !length(
    hep_target_cells
  )
) {

  stop(
    "No target Hepatocyte cells found."
  )
}


# ==============================================================================
# 9. Validate source state vocabularies
# ==============================================================================

mphi_states_present <- sort(
  unique(
    as.character(
      mphi@meta.data[
        mphi_target_cells,
        MPHI_ANNOTATION_COL
      ]
    )
  )
)


missing_mphi_states <- setdiff(
  MPHI_FINAL_STATES,
  mphi_states_present
)


if (
  length(
    missing_mphi_states
  )
) {

  stop(
    "Required final macrophage state(s) missing: ",
    paste(
      missing_mphi_states,
      collapse = ", "
    )
  )
}


hep_states_present <- sort(
  unique(
    as.character(
      hep@meta.data[
        hep_target_cells,
        HEP_STATE_COL
      ]
    )
  )
)


unexpected_hep_states <- setdiff(
  hep_states_present,
  HEP_DETAIL_STATES_EXPECTED
)


if (
  length(
    unexpected_hep_states
  )
) {

  stop(
    "Unexpected final Hepatocyte state(s): ",
    paste(
      unexpected_hep_states,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 10. Restrict Mphi to the final five biological states
# ==============================================================================

mphi_keep_cells <- mphi_target_cells[
  as.character(
    mphi@meta.data[
      mphi_target_cells,
      MPHI_ANNOTATION_COL
    ]
  ) %in%
    MPHI_FINAL_STATES
]


mphi_other_cells <- setdiff(
  mphi_target_cells,
  mphi_keep_cells
)


mphi_exclusion_audit <- tibble(
  reason = c(
    "Retained_final_Mphi5",
    "Excluded_nonfinal_or_Other_Mphi"
  ),
  n_cells = c(
    length(
      mphi_keep_cells
    ),
    length(
      mphi_other_cells
    )
  )
)


write.csv(
  mphi_exclusion_audit,
  file.path(
    TAB_OUT,
    "04_Mphi5_inclusion_exclusion_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Collapse detailed Hepatocyte states to 5 interaction states
# ==============================================================================

hep_detail <- as.character(
  hep@meta.data[
    hep_target_cells,
    HEP_STATE_COL
  ]
)


hep_interaction <- unname(
  HEP_INTERACTION_MAP[
    hep_detail
  ]
)


names(
  hep_interaction
) <- hep_target_cells


hep_keep_cells <- names(
  hep_interaction
)[
  !is.na(
    hep_interaction
  )
]


hep_mt_excluded_cells <- names(
  hep_interaction
)[
  is.na(
    hep_interaction
  )
]


hep_mapping_audit <- tibble(
  detailed_state =
    names(
      HEP_INTERACTION_MAP
    ),
  interaction_state =
    unname(
      HEP_INTERACTION_MAP
    ),
  primary_interaction_included =
    !is.na(
      unname(
        HEP_INTERACTION_MAP
      )
    )
)


write.csv(
  hep_mapping_audit,
  file.path(
    TAB_OUT,
    "05_Hepatocyte_state_mapping_v6.6.0.csv"
  ),
  row.names = FALSE
)


hep_exclusion_audit <- tibble(
  reason = c(
    "Retained_primary_Hep5",
    "Excluded_MT_high_QC_Hepatocyte"
  ),
  n_cells = c(
    length(
      hep_keep_cells
    ),
    length(
      hep_mt_excluded_cells
    )
  )
)


write.csv(
  hep_exclusion_audit,
  file.path(
    TAB_OUT,
    "06_Hepatocyte5_inclusion_exclusion_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Exact barcode match back to frozen whole-cell parent
# ==============================================================================

parent_cells <- colnames(
  parent
)


mphi_match <- exact_barcode_audit(
  mphi_keep_cells,
  parent_cells,
  "Mphi5"
)


hep_match <- exact_barcode_audit(
  hep_keep_cells,
  parent_cells,
  "Hepatocyte5"
)


barcode_audit <- bind_rows(
  mphi_match$summary,
  hep_match$summary
)


write.csv(
  barcode_audit,
  file.path(
    TAB_OUT,
    "07_exact_barcode_match_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


if (
  any(
    barcode_audit$exact_match_rate <
      MIN_EXACT_MATCH_RATE
  )
) {

  stop(
    "Exact barcode match rate < ",
    MIN_EXACT_MATCH_RATE,
    ". See 07_exact_barcode_match_audit_v6.6.0.csv"
  )
}


if (
  length(
    mphi_match$unmatched
  )
) {

  writeLines(
    mphi_match$unmatched,
    file.path(
      LOG_OUT,
      "unmatched_Mphi_barcodes_v6.6.0.txt"
    )
  )
}


if (
  length(
    hep_match$unmatched
  )
) {

  writeLines(
    hep_match$unmatched,
    file.path(
      LOG_OUT,
      "unmatched_Hepatocyte_barcodes_v6.6.0.txt"
    )
  )
}


# ==============================================================================
# 13. Check barcode overlap between Mphi and Hepatocyte sources
# ==============================================================================

cross_lineage_overlap <- intersect(
  mphi_match$matched,
  hep_match$matched
)


overlap_audit <- tibble(
  metric =
    "Exact barcode overlap between Mphi5 and Hepatocyte5",
  n_cells =
    length(
      cross_lineage_overlap
    )
)


write.csv(
  overlap_audit,
  file.path(
    TAB_OUT,
    "08_cross_lineage_barcode_overlap_audit_v6.6.0.csv"
  ),
  row.names = FALSE
)


if (
  length(
    cross_lineage_overlap
  )
) {

  writeLines(
    cross_lineage_overlap,
    file.path(
      LOG_OUT,
      "ERROR_cross_lineage_overlapping_barcodes_v6.6.0.txt"
    )
  )

  stop(
    "Mphi and Hepatocyte barcode sets overlap. ",
    "This must be resolved before interaction analysis."
  )
}


# ==============================================================================
# 14. Build transfer tables from source objects
# ==============================================================================

mphi_transfer <- tibble(
  cell =
    mphi_match$matched,
  mphi_state_detailed =
    as.character(
      mphi@meta.data[
        mphi_match$matched,
        MPHI_ANNOTATION_COL
      ]
    ),
  source_sample =
    as.character(
      mphi@meta.data[
        mphi_match$matched,
        MPHI_SAMPLE_COL
      ]
    )
) %>%
  mutate(
    interaction_celltype =
      unname(
        MPHI_INTERACTION_MAP[
          mphi_state_detailed
        ]
      ),
    interaction_lineage =
      "Macrophage",
    label_source =
      "Mphi_CleanB_v4.14.5"
  )


hep_transfer <- tibble(
  cell =
    hep_match$matched,
  hepatocyte_state_detailed =
    as.character(
      hep@meta.data[
        hep_match$matched,
        HEP_STATE_COL
      ]
    ),
  source_sample =
    as.character(
      hep@meta.data[
        hep_match$matched,
        HEP_SAMPLE_COL
      ]
    )
) %>%
  mutate(
    interaction_celltype =
      unname(
        HEP_INTERACTION_MAP[
          hepatocyte_state_detailed
        ]
      ),
    interaction_lineage =
      "Hepatocyte",
    label_source =
      "Hepatocyte_FINAL_v6.5.2"
  )


if (
  any(
    is.na(
      mphi_transfer$interaction_celltype
    )
  )
) {

  stop(
    "NA interaction label generated in Mphi transfer."
  )
}


if (
  any(
    is.na(
      hep_transfer$interaction_celltype
    )
  )
) {

  stop(
    "NA interaction label generated in Hepatocyte transfer."
  )
}


# ==============================================================================
# 15. Parent sample concordance audit
# ==============================================================================

parent_sample_vector <- as.character(
  parent@meta.data[
    ,
    PARENT_SAMPLE_COL
  ]
)


names(
  parent_sample_vector
) <- rownames(
  parent@meta.data
)


mphi_transfer <- mphi_transfer %>%
  mutate(
    parent_sample =
      parent_sample_vector[
        cell
      ],
    sample_concordant =
      source_sample ==
        parent_sample
  )


hep_transfer <- hep_transfer %>%
  mutate(
    parent_sample =
      parent_sample_vector[
        cell
      ],
    sample_concordant =
      source_sample ==
        parent_sample
  )


sample_concordance <- bind_rows(
  mphi_transfer %>%
    transmute(
      source =
        "Mphi5",
      cell,
      source_sample,
      parent_sample,
      sample_concordant
    ),

  hep_transfer %>%
    transmute(
      source =
        "Hepatocyte5",
      cell,
      source_sample,
      parent_sample,
      sample_concordant
    )
)


sample_concordance_summary <- sample_concordance %>%
  group_by(
    source
  ) %>%
  summarise(
    n_cells =
      n(),
    n_sample_concordant =
      sum(
        sample_concordant
      ),
    concordance_rate =
      mean(
        sample_concordant
      ),
    .groups =
      "drop"
  )


write.csv(
  sample_concordance_summary,
  file.path(
    TAB_OUT,
    "09_source_parent_sample_concordance_v6.6.0.csv"
  ),
  row.names = FALSE
)


if (
  any(
    sample_concordance_summary$concordance_rate <
      1
  )
) {

  discordant <- sample_concordance %>%
    filter(
      !sample_concordant
    )

  write.csv(
    discordant,
    file.path(
      TAB_OUT,
      "09b_source_parent_sample_DISCORDANT_cells_v6.6.0.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Source-parent sample identity mismatch detected. ",
    "See 09b_source_parent_sample_DISCORDANT_cells_v6.6.0.csv"
  )
}


# ==============================================================================
# 16. Subset frozen parent to the exact interaction cells
# ==============================================================================

interaction_cells <- c(
  mphi_transfer$cell,
  hep_transfer$cell
)


interaction_cells <- unique(
  interaction_cells
)


obj <- subset(
  parent,
  cells =
    interaction_cells
)


DefaultAssay(
  obj
) <- "RNA"


obj <- safe_join_rna(
  obj
)


# ==============================================================================
# 17. Initialize interaction metadata
# ==============================================================================

obj[[
  INTERACTION_CELLTYPE_COL
]] <- NA_character_


obj[[
  INTERACTION_LINEAGE_COL
]] <- NA_character_


obj[[
  INTERACTION_SAMPLE_COL
]] <- as.character(
  obj@meta.data[[
    PARENT_SAMPLE_COL
  ]]
)


obj[[
  INTERACTION_CONDITION_COL
]] <- canonical_condition(
  obj@meta.data[[
    INTERACTION_SAMPLE_COL
  ]]
)


obj[[
  MPHI_DETAIL_COL
]] <- NA_character_


obj[[
  HEP_DETAIL_COL
]] <- NA_character_


obj[[
  HEP_INTERACTION_STATE_COL
]] <- NA_character_


obj[[
  SOURCE_COL
]] <- NA_character_


# ==============================================================================
# 18. Transfer macrophage metadata
# ==============================================================================

mphi_idx <- match(
  mphi_transfer$cell,
  rownames(
    obj@meta.data
  )
)


obj@meta.data[
  mphi_idx,
  INTERACTION_CELLTYPE_COL
] <- mphi_transfer$interaction_celltype


obj@meta.data[
  mphi_idx,
  INTERACTION_LINEAGE_COL
] <- "Macrophage"


obj@meta.data[
  mphi_idx,
  MPHI_DETAIL_COL
] <- mphi_transfer$mphi_state_detailed


obj@meta.data[
  mphi_idx,
  SOURCE_COL
] <- mphi_transfer$label_source


# ==============================================================================
# 19. Transfer Hepatocyte metadata
# ==============================================================================

hep_idx <- match(
  hep_transfer$cell,
  rownames(
    obj@meta.data
  )
)


obj@meta.data[
  hep_idx,
  INTERACTION_CELLTYPE_COL
] <- hep_transfer$interaction_celltype


obj@meta.data[
  hep_idx,
  INTERACTION_LINEAGE_COL
] <- "Hepatocyte"


obj@meta.data[
  hep_idx,
  HEP_DETAIL_COL
] <- hep_transfer$hepatocyte_state_detailed


obj@meta.data[
  hep_idx,
  HEP_INTERACTION_STATE_COL
] <- hep_transfer$interaction_celltype


obj@meta.data[
  hep_idx,
  SOURCE_COL
] <- hep_transfer$label_source


# ==============================================================================
# 20. Validate final metadata
# ==============================================================================

required_final_cols <- c(
  INTERACTION_CELLTYPE_COL,
  INTERACTION_LINEAGE_COL,
  INTERACTION_SAMPLE_COL,
  INTERACTION_CONDITION_COL,
  SOURCE_COL
)


for (
  col in required_final_cols
) {

  if (
    any(
      is.na(
        obj@meta.data[[
          col
        ]]
      )
    )
  ) {

    stop(
      "NA values remain in required metadata column: ",
      col
    )
  }
}


if (
  !all(
    as.character(
      obj@meta.data[[
        INTERACTION_SAMPLE_COL
      ]]
    ) %in%
      TARGET_SAMPLES
  )
) {

  stop(
    "Non-target sample found in final interaction object."
  )
}


condition_values <- sort(
  unique(
    as.character(
      obj@meta.data[[
        INTERACTION_CONDITION_COL
      ]]
    )
  )
)


if (
  !setequal(
    condition_values,
    c(
      "Sham",
      "Tx"
    )
  )
) {

  stop(
    "Final interaction object does not contain exactly Sham and Tx."
  )
}


celltype_values <- sort(
  unique(
    as.character(
      obj@meta.data[[
        INTERACTION_CELLTYPE_COL
      ]]
    )
  )
)


missing_interaction_states <- setdiff(
  INTERACTION_LEVELS,
  celltype_values
)


if (
  length(
    missing_interaction_states
  )
) {

  stop(
    "Missing interaction state(s): ",
    paste(
      missing_interaction_states,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 21. Freeze factor levels / identities
# ==============================================================================

obj[[
  INTERACTION_CELLTYPE_COL
]] <- factor(
  as.character(
    obj@meta.data[[
      INTERACTION_CELLTYPE_COL
    ]]
  ),
  levels =
    INTERACTION_LEVELS
)


obj[[
  INTERACTION_LINEAGE_COL
]] <- factor(
  as.character(
    obj@meta.data[[
      INTERACTION_LINEAGE_COL
    ]]
  ),
  levels =
    LINEAGE_LEVELS
)


obj[[
  INTERACTION_SAMPLE_COL
]] <- factor(
  as.character(
    obj@meta.data[[
      INTERACTION_SAMPLE_COL
    ]]
  ),
  levels =
    TARGET_SAMPLES
)


obj[[
  INTERACTION_CONDITION_COL
]] <- factor(
  as.character(
    obj@meta.data[[
      INTERACTION_CONDITION_COL
    ]]
  ),
  levels = c(
    "Sham",
    "Tx"
  )
)


Idents(
  obj
) <- obj@meta.data[[
  INTERACTION_CELLTYPE_COL
]]


# ==============================================================================
# 22. Final cell count tables
# ==============================================================================

count_by_sample_state <- obj@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[
          INTERACTION_SAMPLE_COL
        ]]
      ),
    condition =
      as.character(
        .data[[
          INTERACTION_CONDITION_COL
        ]]
      ),
    lineage =
      as.character(
        .data[[
          INTERACTION_LINEAGE_COL
        ]]
      ),
    interaction_celltype =
      as.character(
        .data[[
          INTERACTION_CELLTYPE_COL
        ]]
      )
  ) %>%
  count(
    sample,
    condition,
    lineage,
    interaction_celltype,
    name =
      "n_cells"
  ) %>%
  complete(
    sample =
      TARGET_SAMPLES,
    interaction_celltype =
      INTERACTION_LEVELS,
    fill =
      list(
        n_cells =
          0
      )
  ) %>%
  mutate(
    condition =
      canonical_condition(
        sample
      ),
    lineage =
      ifelse(
        grepl(
          "^Mphi_",
          interaction_celltype
        ),
        "Macrophage",
        "Hepatocyte"
      )
  ) %>%
  arrange(
    factor(
      sample,
      levels =
        TARGET_SAMPLES
    ),
    factor(
      interaction_celltype,
      levels =
        INTERACTION_LEVELS
    )
  )


write.csv(
  count_by_sample_state,
  file.path(
    TAB_OUT,
    "10_interaction_cell_counts_by_sample_state_v6.6.0.csv"
  ),
  row.names = FALSE
)


count_wide <- count_by_sample_state %>%
  select(
    sample,
    interaction_celltype,
    n_cells
  ) %>%
  pivot_wider(
    names_from =
      sample,
    values_from =
      n_cells
  )


write.csv(
  count_wide,
  file.path(
    TAB_OUT,
    "11_interaction_cell_counts_wide_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 23. Cell-number suitability audit for downstream CellChat
# ==============================================================================

suitability <- count_by_sample_state %>%
  mutate(
    suitability =
      case_when(
        n_cells <
          MIN_CELLS_STOP ~
          "STOP_too_few_cells",

        n_cells <
          MIN_CELLS_WARN ~
          "WARN_low_cells",

        TRUE ~
          "OK"
      )
  )


write.csv(
  suitability,
  file.path(
    TAB_OUT,
    "12_CellChat_cell_number_suitability_v6.6.0.csv"
  ),
  row.names = FALSE
)


if (
  any(
    suitability$n_cells <
      MIN_CELLS_STOP
  )
) {

  bad <- suitability %>%
    filter(
      n_cells <
        MIN_CELLS_STOP
    )

  stop(
    "At least one sample x interaction state has fewer than ",
    MIN_CELLS_STOP,
    " cells. See 12_CellChat_cell_number_suitability_v6.6.0.csv"
  )
}


if (
  any(
    suitability$n_cells <
      MIN_CELLS_WARN
  )
) {

  warning(
    "At least one sample x interaction state has fewer than ",
    MIN_CELLS_WARN,
    " cells. Downstream CellChat should be interpreted cautiously."
  )
}


# ==============================================================================
# 24. Fraction tables within each lineage
# ==============================================================================

lineage_fraction <- obj@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          INTERACTION_SAMPLE_COL
        ]]
      ),
    condition =
      as.character(
        .data[[
          INTERACTION_CONDITION_COL
        ]]
      ),
    lineage =
      as.character(
        .data[[
          INTERACTION_LINEAGE_COL
        ]]
      ),
    interaction_celltype =
      as.character(
        .data[[
          INTERACTION_CELLTYPE_COL
        ]]
      )
  ) %>%
  count(
    sample,
    condition,
    lineage,
    interaction_celltype,
    name =
      "n_cells"
  ) %>%
  group_by(
    sample,
    lineage
  ) %>%
  mutate(
    lineage_total =
      sum(
        n_cells
      ),
    fraction_within_lineage =
      n_cells /
        lineage_total
  ) %>%
  ungroup()


write.csv(
  lineage_fraction,
  file.path(
    TAB_OUT,
    "13_interaction_state_fraction_within_lineage_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 25. Detailed Hepatocyte -> interaction-state provenance
# ==============================================================================

hep_provenance <- obj@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  filter(
    as.character(
      .data[[
        INTERACTION_LINEAGE_COL
      ]]
    ) ==
      "Hepatocyte"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          INTERACTION_SAMPLE_COL
        ]]
      ),
    detailed_state =
      as.character(
        .data[[
          HEP_DETAIL_COL
        ]]
      ),
    interaction_state =
      as.character(
        .data[[
          INTERACTION_CELLTYPE_COL
        ]]
      )
  ) %>%
  count(
    sample,
    detailed_state,
    interaction_state,
    name =
      "n_cells"
  ) %>%
  arrange(
    factor(
      sample,
      levels =
        TARGET_SAMPLES
    ),
    interaction_state,
    detailed_state
  )


write.csv(
  hep_provenance,
  file.path(
    TAB_OUT,
    "14_Hepatocyte_detailed_to_interaction_provenance_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 26. Macrophage detailed-state provenance
# ==============================================================================

mphi_provenance <- obj@meta.data %>%
  as_tibble(
    rownames =
      "cell"
  ) %>%
  filter(
    as.character(
      .data[[
        INTERACTION_LINEAGE_COL
      ]]
    ) ==
      "Macrophage"
  ) %>%
  transmute(
    sample =
      as.character(
        .data[[
          INTERACTION_SAMPLE_COL
        ]]
      ),
    detailed_state =
      as.character(
        .data[[
          MPHI_DETAIL_COL
        ]]
      ),
    interaction_state =
      as.character(
        .data[[
          INTERACTION_CELLTYPE_COL
        ]]
      )
  ) %>%
  count(
    sample,
    detailed_state,
    interaction_state,
    name =
      "n_cells"
  ) %>%
  arrange(
    factor(
      sample,
      levels =
        TARGET_SAMPLES
    ),
    interaction_state,
    detailed_state
  )


write.csv(
  mphi_provenance,
  file.path(
    TAB_OUT,
    "15_Mphi_detailed_to_interaction_provenance_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 27. Candidate communication-gene availability audit
# ==============================================================================

MPHI_TO_HEP_CANDIDATES <- c(
  "Tnf",
  "Il1a",
  "Il1b",
  "Il6",
  "Osm",
  "Tgfb1",
  "Spp1",
  "Fn1",
  "Pdgfb",
  "Gas6",
  "Areg",
  "Hbegf",
  "Igf1",
  "Ccl2",
  "Ccl3",
  "Ccl4",
  "Cxcl2",
  "Cxcl10"
)


HEP_RECEPTOR_CANDIDATES <- c(
  "Tnfrsf1a",
  "Tnfrsf1b",
  "Il1r1",
  "Il6ra",
  "Il6st",
  "Osmr",
  "Tgfbr1",
  "Tgfbr2",
  "Cd44",
  "Itgav",
  "Itgb1",
  "Egfr",
  "Igf1r",
  "Axl",
  "Mertk"
)


HEP_TO_MPHI_CANDIDATES <- c(
  "Ccl2",
  "Cxcl1",
  "Cxcl10",
  "Csf1",
  "Il6",
  "Spp1",
  "Fn1",
  "Saa1",
  "Saa2",
  "Orm1",
  "Orm2",
  "Areg",
  "Shh",
  "Ihh"
)


MPHI_RECEPTOR_CANDIDATES <- c(
  "Ccr2",
  "Ccr5",
  "Cxcr2",
  "Cxcr3",
  "Csf1r",
  "Il6ra",
  "Il6st",
  "Tnfrsf1a",
  "Tnfrsf1b",
  "Tgfbr1",
  "Tgfbr2",
  "Cd44",
  "Itgav",
  "Itgb1",
  "Axl",
  "Mertk"
)


gene_audit <- bind_rows(
  tibble(
    direction =
      "Mphi_to_Hep",
    role =
      "Candidate_ligand",
    gene =
      MPHI_TO_HEP_CANDIDATES
  ),
  tibble(
    direction =
      "Mphi_to_Hep",
    role =
      "Candidate_receptor",
    gene =
      HEP_RECEPTOR_CANDIDATES
  ),
  tibble(
    direction =
      "Hep_to_Mphi",
    role =
      "Candidate_ligand",
    gene =
      HEP_TO_MPHI_CANDIDATES
  ),
  tibble(
    direction =
      "Hep_to_Mphi",
    role =
      "Candidate_receptor",
    gene =
      MPHI_RECEPTOR_CANDIDATES
  )
) %>%
  mutate(
    present_in_interaction_object =
      gene %in%
        rownames(
          obj
        )
  )


write.csv(
  gene_audit,
  file.path(
    TAB_OUT,
    "16_candidate_communication_gene_availability_v6.6.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 28. Candidate gene expression audit by interaction state
# ==============================================================================

candidate_genes_present <- present_genes(
  obj,
  unique(
    gene_audit$gene
  )
)


if (
  length(
    candidate_genes_present
  )
) {

  counts_mat <- GetAssayData(
    obj,
    assay =
      "RNA",
    layer =
      "counts"
  )

  data_mat <- GetAssayData(
    obj,
    assay =
      "RNA",
    layer =
      "data"
  )

  expression_rows <- list()

  for (
    state_name in INTERACTION_LEVELS
  ) {

    cells <- rownames(
      obj@meta.data
    )[
      as.character(
        obj@meta.data[[
          INTERACTION_CELLTYPE_COL
        ]]
      ) ==
        state_name
    ]

    if (
      !length(
        cells
      )
    ) {
      next
    }

    pct <- Matrix::rowMeans(
      counts_mat[
        candidate_genes_present,
        cells,
        drop = FALSE
      ] >
        0
    )

    avg <- Matrix::rowMeans(
      data_mat[
        candidate_genes_present,
        cells,
        drop = FALSE
      ]
    )

    expression_rows[[
      state_name
    ]] <- tibble(
      interaction_celltype =
        state_name,
      gene =
        candidate_genes_present,
      pct_expressed =
        as.numeric(
          pct
        ),
      mean_log_normalized_expression =
        as.numeric(
          avg
        )
    )
  }

  candidate_expression <- bind_rows(
    expression_rows
  )

  write.csv(
    candidate_expression,
    file.path(
      TAB_OUT,
      "17_candidate_communication_gene_expression_by_state_v6.6.0.csv"
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 29. UMAP figures from the frozen whole-cell parent reduction
# ==============================================================================

UMAP_NAME <-
  "umapRPCA"


if (
  UMAP_NAME %in%
    Reductions(
      obj
    )
) {

  interaction_colors <- c(
    "Mphi_Anti-inflammatory" =
      "#00A651",
    "Mphi_Inflammatory" =
      "#E31A1C",
    "Mphi_ECM-associated-inflammatory" =
      "#E7298A",
    "Mphi_Repair-Resolution" =
      "#377EB8",
    "Mphi_Lipid-associated-TREM2" =
      "#FF8C00",
    "Hep_Periportal" =
      "#008C95",
    "Hep_Pericentral" =
      "#00BFC4",
    "Hep_Injury-inflammatory" =
      "#7B3294",
    "Hep_Intermediate" =
      "#66C2A5",
    "Hep_Cycling" =
      "#E6AB02"
  )


  p_umap <- DimPlot(
    obj,
    reduction =
      UMAP_NAME,
    group.by =
      INTERACTION_CELLTYPE_COL,
    cols =
      interaction_colors,
    pt.size =
      0.35,
    raster =
      FALSE,
    label =
      FALSE
  ) +
    ggtitle(
      "Mphi5 x Hep5 interaction-ready cells | v6.6.0"
    ) +
    theme_classic(
      base_size =
        9
    )


  save_pdf(
    p_umap,
    file.path(
      FIG_OUT,
      "01_Mphi5_Hep5_interaction_UMAP_v6.6.0.pdf"
    ),
    10,
    8
  )


  p_sample <- DimPlot(
    obj,
    reduction =
      UMAP_NAME,
    group.by =
      INTERACTION_CELLTYPE_COL,
    split.by =
      INTERACTION_SAMPLE_COL,
    cols =
      interaction_colors,
    pt.size =
      0.28,
    raster =
      FALSE,
    ncol =
      2
  ) +
    plot_annotation(
      title =
        "Mphi5 x Hep5 interaction-ready cells by biological sample"
    )


  save_pdf(
    p_sample,
    file.path(
      FIG_OUT,
      "02_Mphi5_Hep5_interaction_UMAP_by_sample_v6.6.0.pdf"
    ),
    15,
    11
  )


  p_lineage <- DimPlot(
    obj,
    reduction =
      UMAP_NAME,
    group.by =
      INTERACTION_LINEAGE_COL,
    pt.size =
      0.30,
    raster =
      FALSE
  ) +
    ggtitle(
      "Interaction-ready lineage overview"
    ) +
    theme_classic(
      base_size =
        9
    )


  save_pdf(
    p_lineage,
    file.path(
      FIG_OUT,
      "03_Mphi_Hep_lineage_UMAP_v6.6.0.pdf"
    ),
    8,
    7
  )

} else {

  warning(
    "Reduction ",
    UMAP_NAME,
    " not found in the frozen parent. UMAP figures skipped."
  )
}


# ==============================================================================
# 30. Cell-number heatmap
# ==============================================================================

count_heat <- count_by_sample_state %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          TARGET_SAMPLES
      ),
    interaction_celltype =
      factor(
        interaction_celltype,
        levels =
          rev(
            INTERACTION_LEVELS
          )
      )
  )


p_count_heat <- ggplot(
  count_heat,
  aes(
    x =
      sample,
    y =
      interaction_celltype,
    fill =
      n_cells
  )
) +
  geom_tile(
    linewidth =
      0.35
  ) +
  geom_text(
    aes(
      label =
        n_cells
    ),
    size =
      3
  ) +
  scale_fill_gradientn(
    colours = c(
      "#FFFFFF",
      "#00BFC4",
      "#0033FF"
    )
  ) +
  labs(
    title =
      "Interaction-ready cell counts",
    subtitle =
      "Exact barcode transfer | Sham1, Sham20, Tx17, Tx5",
    x =
      NULL,
    y =
      NULL,
    fill =
      "Cells"
  ) +
  theme_classic(
    base_size =
      9
  )


save_pdf(
  p_count_heat,
  file.path(
    FIG_OUT,
    "04_Mphi5_Hep5_cell_count_heatmap_v6.6.0.pdf"
  ),
  9,
  7
)


# ==============================================================================
# 31. Fraction plots within each lineage
# ==============================================================================

fraction_plot_df <- lineage_fraction %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          TARGET_SAMPLES
      ),
    interaction_celltype =
      factor(
        interaction_celltype,
        levels =
          INTERACTION_LEVELS
      )
  )


p_fraction <- ggplot(
  fraction_plot_df,
  aes(
    x =
      sample,
    y =
      fraction_within_lineage,
    group =
      1
  )
) +
  geom_line(
    linewidth =
      0.5
  ) +
  geom_point(
    size =
      2
  ) +
  facet_wrap(
    lineage ~ interaction_celltype,
    scales =
      "free_y",
    ncol =
      5
  ) +
  labs(
    title =
      "Interaction-state fractions within lineage",
    subtitle =
      "Exploratory abundance context for later communication analysis",
    x =
      NULL,
    y =
      "Fraction within lineage"
  ) +
  theme_classic(
    base_size =
      7
  )


save_pdf(
  p_fraction,
  file.path(
    FIG_OUT,
    "05_Mphi5_Hep5_fraction_within_lineage_v6.6.0.pdf"
  ),
  16,
  8
)


# ==============================================================================
# 32. Candidate communication DotPlot
# ==============================================================================

DOT_GENES <- unique(
  c(
    "Tnf",
    "Il1b",
    "Il6",
    "Osm",
    "Tgfb1",
    "Spp1",
    "Fn1",
    "Pdgfb",
    "Gas6",
    "Areg",
    "Igf1",
    "Ccl2",
    "Cxcl1",
    "Cxcl10",
    "Csf1",
    "Saa1",
    "Saa2",
    "Shh",
    "Ihh",
    "Tnfrsf1a",
    "Il1r1",
    "Il6st",
    "Osmr",
    "Tgfbr1",
    "Tgfbr2",
    "Cd44",
    "Egfr",
    "Igf1r",
    "Axl",
    "Mertk"
  )
)


DOT_GENES_PRESENT <- present_genes(
  obj,
  DOT_GENES
)


if (
  length(
    DOT_GENES_PRESENT
  )
) {

  p_dot <- DotPlot(
    obj,
    features =
      DOT_GENES_PRESENT,
    group.by =
      INTERACTION_CELLTYPE_COL,
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
      "Candidate Mphi <-> Hepatocyte communication genes"
    )


  save_pdf(
    p_dot,
    file.path(
      FIG_OUT,
      "06_candidate_communication_DotPlot_v6.6.0.pdf"
    ),
    16,
    8
  )
}


# ==============================================================================
# 33. Save interaction-ready RDS
# ==============================================================================

RDS_FILE <- file.path(
  RDS_OUT,
  "Mouse_MASH_Mphi5_Hep5_interaction_ready_v6.6.0.rds"
)


saveRDS(
  obj,
  RDS_FILE,
  compress =
    FALSE
)


msg(
  "Saved interaction-ready RDS: ",
  RDS_FILE
)


# ==============================================================================
# 34. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "wholecell_parent_RDS",
    "Mphi_source_RDS",
    "Hepatocyte_source_RDS",
    "target_samples",
    "Mphi_annotation_source_column",
    "Hepatocyte_annotation_source_column",
    "interaction_celltype_column",
    "interaction_lineage_column",
    "interaction_sample_column",
    "interaction_condition_column",
    "n_Mphi_interaction_states",
    "n_Hepatocyte_interaction_states",
    "Hepatocyte_MT_high_QC_in_primary_interaction",
    "barcode_transfer_method",
    "minimum_exact_match_rate",
    "reclustering_performed",
    "communication_analysis_performed"
  ),
  value = c(
    "v6.6.0",
    WHOLECELL_RDS,
    MPHI_RDS,
    HEP_RDS,
    paste(
      TARGET_SAMPLES,
      collapse = ","
    ),
    MPHI_ANNOTATION_COL,
    HEP_STATE_COL,
    INTERACTION_CELLTYPE_COL,
    INTERACTION_LINEAGE_COL,
    INTERACTION_SAMPLE_COL,
    INTERACTION_CONDITION_COL,
    "5",
    "5",
    "FALSE; retained only in source/exclusion audit",
    "Exact cell-barcode matching to frozen whole-cell parent",
    as.character(
      MIN_EXACT_MATCH_RATE
    ),
    "FALSE",
    "FALSE; this is preparation only"
  )
)


write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.0.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.0.txt"
    )
)


# ==============================================================================
# 35. Final console summary
# ==============================================================================

msg(
  "FINAL interaction cell counts:"
)


print(
  count_wide
)


msg(
  "Exact barcode match audit:"
)


print(
  barcode_audit
)


msg(
  "DONE."
)


msg(
  "Output directory: ",
  OUT
)
