# ============================================================
# Human MASH v1.3.3
# Whole-liver Annotation Freeze
#
# PURPOSE
#   Create a new canonical Human MASH LT1-LT8 object by merging:
#
#     1. v1.3 canonical whole-liver object (READ ONLY)
#     2. v1.3.2A.1 immune manual annotation / cleaning
#     3. v1.3.2B.1 hepatobiliary manual annotation / cleaning
#
#   while preserving validated structural annotations from v1.3.
#
# IMPORTANT
#   - The original canonical v1.3 RDS is NEVER modified.
#   - All 81,348 canonical cells are retained.
#   - Refined immune / hepatobiliary annotations are merged by
#     exact cell barcode.
#   - Cells excluded from a refinement are NOT deleted from the
#     whole-liver canonical object. Instead, where the excluded
#     cluster has a validated alternate lineage, that lineage is
#     used as a cross-lineage correction.
#   - Conflicting clean-clean overlap between immune and
#     hepatobiliary refinements causes a hard STOP.
#   - If BOTH refinements exclude the same cell but propose
#     DIFFERENT alternate lineages, the cell is retained as:
#         celltype_final_v133 = Ambiguous_crosslineage
#         cellstate_final_v133 = Doublet_or_crosslineage_ambiguous
#     rather than being force-assigned to either lineage.
#   - Original v1.3 annotation columns remain untouched for
#     provenance.
#
# FINAL PRIMARY METADATA
#
#   celltype_final_v133
#       final broad whole-liver cell type
#
#   cellstate_final_v133
#       refined immune/hepatobiliary state when available;
#       otherwise carried-forward structural type
#
#   annotation_source_v133
#       source of the final annotation decision
#
#   annotation_frozen_v133
#       TRUE for all cells after successful v1.3.3 freeze
#
#   state_sample_bias_ge70_v133
#       TRUE for known strongly sample-biased hepatobiliary
#       states (clusters 9, 12, 14)
#
# UMAP POLICY (v1.3.3+)
#   - publication/review PDF uses vector points: raster = FALSE
#   - cairo_pdf when available
#   - overview pt.size = 0.18
#   - split-by-sample pt.size = 0.12
#   - high-contrast saturated palette
#   - no ordered rasterization / no ggrastr dependency
#
# OUTPUT
#   New canonical directory:
#     Human_MASH_LT1to8_Canonical_v1.3.3/
#
#   New object:
#     objects/Human_MASH_LT1to8_Canonical_v1.3.3.rds
#
# STATUS
#   WHOLE_LIVER_ANNOTATION_FROZEN_v1.3.3
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

seed <- 20260908
set.seed(seed)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

base_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS"
)

canonical_v13_rds <- file.path(
  base_dir,
  "Human_MASH_LT1to8_Canonical_v1.3",
  "objects",
  "Human_MASH_LT1to8_Canonical_v1.3.rds"
)

immune_outdir <- file.path(
  base_dir,
  "Human_MASH_LT1to8_Canonical_v1.3",
  "analysis",
  "v1.3.2A.1_immune_manual_annotation_cleaning"
)

immune_metadata_csv <- file.path(
  immune_outdir,
  "16_v1.3.2A.1_all_candidate_cell_metadata.csv"
)

immune_flags_txt <- file.path(
  immune_outdir,
  "18_v1.3.2A.1_review_flags.txt"
)

hepato_outdir <- file.path(
  base_dir,
  "Human_MASH_LT1to8_Canonical_v1.3",
  "analysis",
  "v1.3.2B.1_hepatobiliary_manual_annotation_cleaning"
)

hepato_metadata_csv <- file.path(
  hepato_outdir,
  "15_v1.3.2B.1_all_candidate_cell_metadata.csv"
)

hepato_flags_txt <- file.path(
  hepato_outdir,
  "18_v1.3.2B.1_review_flags.txt"
)

outside_core_ductal_csv <- file.path(
  hepato_outdir,
  "11_outside_core_ductal_candidates_AUDIT_ONLY.csv"
)

canonical_v133_dir <- file.path(
  base_dir,
  "Human_MASH_LT1to8_Canonical_v1.3.3"
)

object_dir <- file.path(
  canonical_v133_dir,
  "objects"
)

outdir <- file.path(
  canonical_v133_dir,
  "analysis",
  "v1.3.3_whole_liver_annotation_freeze"
)

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

final_rds <- file.path(
  object_dir,
  "Human_MASH_LT1to8_Canonical_v1.3.3.rds"
)

expected_canonical_cells <- 81348L
expected_immune_candidate_cells <- 14833L
expected_hepatobiliary_candidate_cells <- 47425L

celltype_old_col <- "celltype_ref12_v13"
sample_col <- "canonical_sample_v13"
nas_col <- "NAS_score_v13"
nasgrp_col <- "NAS_group_v13"

native_umap <- "umap"
rpca_umap <- "umapRPCA_ref12_v13"

overview_pt_size <- 0.18
split_pt_size <- 0.12

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

safe_chr <- function(x) {

  out <- as.character(x)

  out[
    is.na(out) |
      out == ""
  ] <- NA_character_

  out
}

safe_logical <- function(x) {

  if (
    is.logical(x)
  ) {

    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA,
    length(y)
  )

  out[
    y %in%
      c(
        "TRUE",
        "T",
        "1"
      )
  ] <- TRUE

  out[
    y %in%
      c(
        "FALSE",
        "F",
        "0"
      )
  ] <- FALSE

  out
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

save_pdf <- function(
    filename,
    plot,
    width,
    height
) {

  if (
    capabilities(
      "cairo"
    )
  ) {

    ggsave(
      filename = filename,
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height
    )

  } else {

    warning(
      "R was built without Cairo support; falling back to standard PDF device."
    )

    ggsave(
      filename = filename,
      plot = plot,
      device = grDevices::pdf,
      width = width,
      height = height
    )
  }
}

# ------------------------------------------------------------
# 2. Validate required files
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.3 Whole-liver Annotation Freeze\n")
cat("============================================================\n\n")

required_files <- c(
  canonical_v13_rds,
  immune_metadata_csv,
  immune_flags_txt,
  hepato_metadata_csv,
  hepato_flags_txt
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing_files
  ) > 0
) {

  stop(
    "Missing required input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 3. Load canonical v1.3
# ------------------------------------------------------------

obj <- readRDS(
  canonical_v13_rds
)

DefaultAssay(
  obj
) <- "RNA"

required_meta <- c(
  celltype_old_col,
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
    "Canonical v1.3 lacks required metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  ncol(
    obj
  ) != expected_canonical_cells
) {

  stop(
    "Canonical v1.3 cell-count mismatch. Expected ",
    expected_canonical_cells,
    " but found ",
    ncol(
      obj
    ),
    "."
  )
}

if (
  anyDuplicated(
    colnames(
      obj
    )
  )
) {

  stop(
    "Canonical v1.3 contains duplicated cell names."
  )
}

for (
  red in c(
    native_umap,
    rpca_umap
  )
) {

  if (
    !red %in%
      Reductions(
        obj
      )
  ) {

    stop(
      "Required reduction missing from canonical v1.3: ",
      red
    )
  }
}

canonical_cells <- colnames(
  obj
)

old_celltype <- safe_chr(
  obj[[celltype_old_col]][, 1]
)

canonical_sample <- safe_chr(
  obj[[sample_col]][, 1]
)

cat(
  "Canonical cells:",
  ncol(
    obj
  ),
  "\n"
)

cat(
  "Old annotation types:",
  length(
    unique(
      old_celltype
    )
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 4. Load refined immune / hepatobiliary metadata
# ------------------------------------------------------------

immune <- read.csv(
  immune_metadata_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

hepato <- read.csv(
  hepato_metadata_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_immune_cols <- c(
  "cell",
  "sample",
  "cluster_v132A",
  "manual_v132A1",
  "broad_lineage_v132A1",
  "keep_clean_v132A1",
  "exclusion_reason_v132A1"
)

missing_immune_cols <- setdiff(
  required_immune_cols,
  colnames(
    immune
  )
)

if (
  length(
    missing_immune_cols
  ) > 0
) {

  stop(
    "Immune metadata lacks required columns: ",
    paste(
      missing_immune_cols,
      collapse = ", "
    )
  )
}

required_hepato_cols <- c(
  "cell",
  "sample",
  "cluster_v132B",
  "manual_state_v132B1",
  "broad_celltype_v132B1",
  "manual_confidence_v132B1",
  "keep_clean_v132B1",
  "status_v132B1"
)

missing_hepato_cols <- setdiff(
  required_hepato_cols,
  colnames(
    hepato
  )
)

if (
  length(
    missing_hepato_cols
  ) > 0
) {

  stop(
    "Hepatobiliary metadata lacks required columns: ",
    paste(
      missing_hepato_cols,
      collapse = ", "
    )
  )
}

if (
  nrow(
    immune
  ) != expected_immune_candidate_cells
) {

  stop(
    "Immune candidate-cell count mismatch. Expected ",
    expected_immune_candidate_cells,
    " but found ",
    nrow(
      immune
    ),
    "."
  )
}

if (
  nrow(
    hepato
  ) != expected_hepatobiliary_candidate_cells
) {

  stop(
    "Hepatobiliary candidate-cell count mismatch. Expected ",
    expected_hepatobiliary_candidate_cells,
    " but found ",
    nrow(
      hepato
    ),
    "."
  )
}

if (
  anyDuplicated(
    immune$cell
  )
) {

  stop(
    "Immune metadata contains duplicated cell names."
  )
}

if (
  anyDuplicated(
    hepato$cell
  )
) {

  stop(
    "Hepatobiliary metadata contains duplicated cell names."
  )
}

unknown_immune_cells <- setdiff(
  immune$cell,
  canonical_cells
)

unknown_hepato_cells <- setdiff(
  hepato$cell,
  canonical_cells
)

if (
  length(
    unknown_immune_cells
  ) > 0
) {

  stop(
    "Immune metadata contains cells not found in canonical v1.3: ",
    length(
      unknown_immune_cells
    )
  )
}

if (
  length(
    unknown_hepato_cells
  ) > 0
) {

  stop(
    "Hepatobiliary metadata contains cells not found in canonical v1.3: ",
    length(
      unknown_hepato_cells
    )
  )
}

# ------------------------------------------------------------
# 5. Match refined metadata to canonical cells
# ------------------------------------------------------------

mi <- match(
  canonical_cells,
  immune$cell
)

mh <- match(
  canonical_cells,
  hepato$cell
)

immune_member <- !is.na(
  mi
)

hepato_member <- !is.na(
  mh
)

immune_keep <- rep(
  NA,
  length(
    canonical_cells
  )
)

hepato_keep <- rep(
  NA,
  length(
    canonical_cells
  )
)

immune_keep[
  immune_member
] <- safe_logical(
  immune$keep_clean_v132A1[
    mi[
      immune_member
    ]
  ]
)

hepato_keep[
  hepato_member
] <- safe_logical(
  hepato$keep_clean_v132B1[
    mh[
      hepato_member
    ]
  ]
)

if (
  anyNA(
    immune_keep[
      immune_member
    ]
  )
) {

  stop(
    "Unable to parse immune keep_clean_v132A1."
  )
}

if (
  anyNA(
    hepato_keep[
      hepato_member
    ]
  )
) {

  stop(
    "Unable to parse hepatobiliary keep_clean_v132B1."
  )
}

# Sample consistency
immune_sample_match <- rep(
  TRUE,
  length(
    canonical_cells
  )
)

immune_sample_match[
  immune_member
] <-
  safe_chr(
    immune$sample[
      mi[
        immune_member
      ]
    ]
  ) ==
  canonical_sample[
    immune_member
  ]

hepato_sample_match <- rep(
  TRUE,
  length(
    canonical_cells
  )
)

hepato_sample_match[
  hepato_member
] <-
  safe_chr(
    hepato$sample[
      mh[
        hepato_member
      ]
    ]
  ) ==
  canonical_sample[
    hepato_member
  ]

if (
  any(
    !immune_sample_match,
    na.rm = TRUE
  )
) {

  stop(
    "Sample mismatch between canonical and immune metadata."
  )
}

if (
  any(
    !hepato_sample_match,
    na.rm = TRUE
  )
) {

  stop(
    "Sample mismatch between canonical and hepatobiliary metadata."
  )
}

# ------------------------------------------------------------
# 6. Build immune proposed final annotation
# ------------------------------------------------------------

immune_type_map_clean <- c(
  "Macrophage" =
    "Macrophage",
  "Dendritic" =
    "Dendritic_cell",
  "T_cell" =
    "T_cell",
  "Monocyte" =
    "Monocyte",
  "NK_cell" =
    "NK_cell",
  "NK_T" =
    "NK_gamma_delta_T",
  "Plasma_cell" =
    "Plasma_cell",
  "Cycling_immune" =
    "Cycling_immune",
  "B_cell" =
    "B_cell",
  "Mast_cell" =
    "Mast_cell",
  "Neutrophil" =
    "Neutrophil"
)

immune_type_map_excluded <- c(
  "Exclude_hepatocyte_like" =
    "Hepatocyte",
  "Exclude_Schwann_neural_like" =
    "Schwann_neural_like",
  "Exclude_perivascular_mesenchymal" =
    "Mesenchymal_perivascular",
  "Exclude_endothelial" =
    "Endothelial",
  "Exclude_erythroid" =
    "Erythroid",
  "Exclude_ductal_like" =
    "Cholangiocyte"
)

immune_state_map_excluded <- c(
  "Exclude_hepatocyte_like" =
    "Reassigned_hepatocyte_like",
  "Exclude_Schwann_neural_like" =
    "Schwann_neural_like",
  "Exclude_perivascular_mesenchymal" =
    "Mesenchymal_perivascular",
  "Exclude_endothelial" =
    "Reassigned_endothelial",
  "Exclude_erythroid" =
    "Erythroid",
  "Exclude_ductal_like" =
    "Reassigned_ductal_like"
)

immune_manual <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

immune_broad_raw <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

immune_cluster <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

immune_manual[
  immune_member
] <- safe_chr(
  immune$manual_v132A1[
    mi[
      immune_member
    ]
  ]
)

immune_broad_raw[
  immune_member
] <- safe_chr(
  immune$broad_lineage_v132A1[
    mi[
      immune_member
    ]
  ]
)

immune_cluster[
  immune_member
] <- safe_chr(
  immune$cluster_v132A[
    mi[
      immune_member
    ]
  ]
)

immune_proposed_type <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

immune_proposed_state <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

immune_proposal_kind <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

ii_clean <- which(
  immune_member &
    immune_keep
)

ii_excl <- which(
  immune_member &
    !immune_keep
)

immune_proposed_type[
  ii_clean
] <- unname(
  immune_type_map_clean[
    immune_broad_raw[
      ii_clean
    ]
  ]
)

immune_proposed_state[
  ii_clean
] <- immune_manual[
  ii_clean
]

immune_proposal_kind[
  ii_clean
] <- "clean"

immune_proposed_type[
  ii_excl
] <- unname(
  immune_type_map_excluded[
    immune_manual[
      ii_excl
    ]
  ]
)

immune_proposed_state[
  ii_excl
] <- unname(
  immune_state_map_excluded[
    immune_manual[
      ii_excl
    ]
  ]
)

immune_proposal_kind[
  ii_excl
] <- "excluded_correction"

if (
  anyNA(
    immune_proposed_type[
      immune_member
    ]
  )
) {

  bad <- unique(
    data.frame(
      manual =
        immune_manual[
          immune_member &
            is.na(
              immune_proposed_type
            )
        ],
      broad =
        immune_broad_raw[
          immune_member &
            is.na(
              immune_proposed_type
            )
        ],
      keep =
        immune_keep[
          immune_member &
            is.na(
              immune_proposed_type
            )
        ],
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    bad,
    file.path(
      outdir,
      "ERROR_unmapped_immune_annotations.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Unmapped immune annotation(s). See ERROR_unmapped_immune_annotations.csv"
  )
}

# ------------------------------------------------------------
# 7. Build hepatobiliary proposed final annotation
# ------------------------------------------------------------

hepato_excluded_type_map <- c(
  "Exclude_immune_T_cell" =
    "T_cell",
  "Exclude_mesenchymal" =
    "Mesenchymal_perivascular",
  "Exclude_macrophage" =
    "Macrophage"
)

hepato_excluded_state_map <- c(
  "Exclude_immune_T_cell" =
    "Reassigned_T_cell",
  "Exclude_mesenchymal" =
    "Mesenchymal_perivascular",
  "Exclude_macrophage" =
    "Reassigned_macrophage"
)

hepato_state <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_broad <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_cluster <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_confidence <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_state[
  hepato_member
] <- safe_chr(
  hepato$manual_state_v132B1[
    mh[
      hepato_member
    ]
  ]
)

hepato_broad[
  hepato_member
] <- safe_chr(
  hepato$broad_celltype_v132B1[
    mh[
      hepato_member
    ]
  ]
)

hepato_cluster[
  hepato_member
] <- safe_chr(
  hepato$cluster_v132B[
    mh[
      hepato_member
    ]
  ]
)

hepato_confidence[
  hepato_member
] <- safe_chr(
  hepato$manual_confidence_v132B1[
    mh[
      hepato_member
    ]
  ]
)

hepato_proposed_type <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_proposed_state <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

hepato_proposal_kind <- rep(
  NA_character_,
  length(
    canonical_cells
  )
)

ih_clean <- which(
  hepato_member &
    hepato_keep
)

ih_excl <- which(
  hepato_member &
    !hepato_keep
)

hepato_proposed_type[
  ih_clean
] <- hepato_broad[
  ih_clean
]

hepato_proposed_state[
  ih_clean
] <- hepato_state[
  ih_clean
]

hepato_proposal_kind[
  ih_clean
] <- "clean"

hepato_proposed_type[
  ih_excl
] <- unname(
  hepato_excluded_type_map[
    hepato_state[
      ih_excl
    ]
  ]
)

hepato_proposed_state[
  ih_excl
] <- unname(
  hepato_excluded_state_map[
    hepato_state[
      ih_excl
    ]
  ]
)

hepato_proposal_kind[
  ih_excl
] <- "excluded_correction"

if (
  anyNA(
    hepato_proposed_type[
      hepato_member
    ]
  )
) {

  bad <- unique(
    data.frame(
      state =
        hepato_state[
          hepato_member &
            is.na(
              hepato_proposed_type
            )
        ],
      broad =
        hepato_broad[
          hepato_member &
            is.na(
              hepato_proposed_type
            )
        ],
      keep =
        hepato_keep[
          hepato_member &
            is.na(
              hepato_proposed_type
            )
        ],
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    bad,
    file.path(
      outdir,
      "ERROR_unmapped_hepatobiliary_annotations.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Unmapped hepatobiliary annotation(s). See ERROR_unmapped_hepatobiliary_annotations.csv"
  )
}

# ------------------------------------------------------------
# 8. Overlap / conflict audit
# ------------------------------------------------------------

overlap <- immune_member &
  hepato_member

overlap_df <- data.frame(
  cell =
    canonical_cells[
      overlap
    ],
  sample =
    canonical_sample[
      overlap
    ],
  old_celltype =
    old_celltype[
      overlap
    ],
  immune_cluster =
    immune_cluster[
      overlap
    ],
  immune_keep =
    immune_keep[
      overlap
    ],
  immune_type =
    immune_proposed_type[
      overlap
    ],
  immune_state =
    immune_proposed_state[
      overlap
    ],
  hepato_cluster =
    hepato_cluster[
      overlap
    ],
  hepato_keep =
    hepato_keep[
      overlap
    ],
  hepato_type =
    hepato_proposed_type[
      overlap
    ],
  hepato_state =
    hepato_proposed_state[
      overlap
    ],
  stringsAsFactors = FALSE
)

write.csv(
  overlap_df,
  file.path(
    outdir,
    "01_immune_hepatobiliary_overlap_audit.csv"
  ),
  row.names = FALSE
)

clean_clean_conflict <- overlap &
  immune_keep &
  hepato_keep

if (
  any(
    clean_clean_conflict
  )
) {

  conflict <- data.frame(
    cell =
      canonical_cells[
        clean_clean_conflict
      ],
    sample =
      canonical_sample[
        clean_clean_conflict
      ],
    immune_type =
      immune_proposed_type[
        clean_clean_conflict
      ],
    immune_state =
      immune_proposed_state[
        clean_clean_conflict
      ],
    hepato_type =
      hepato_proposed_type[
        clean_clean_conflict
      ],
    hepato_state =
      hepato_proposed_state[
        clean_clean_conflict
      ],
    stringsAsFactors = FALSE
  )

  write.csv(
    conflict,
    file.path(
      outdir,
      "ERROR_clean_clean_overlap_conflict.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Clean immune and clean hepatobiliary annotations overlap. ",
    "See ERROR_clean_clean_overlap_conflict.csv"
  )
}

both_excluded <- overlap &
  !immune_keep &
  !hepato_keep

excluded_excluded_agree <- both_excluded &
  immune_proposed_type ==
    hepato_proposed_type

excluded_excluded_agree[
  is.na(
    excluded_excluded_agree
  )
] <- FALSE

excluded_excluded_conflict <- both_excluded &
  immune_proposed_type !=
    hepato_proposed_type

excluded_excluded_conflict[
  is.na(
    excluded_excluded_conflict
  )
] <- FALSE

conflict <- data.frame(
  cell =
    canonical_cells[
      excluded_excluded_conflict
    ],
  sample =
    canonical_sample[
      excluded_excluded_conflict
    ],
  old_celltype =
    old_celltype[
      excluded_excluded_conflict
    ],
  immune_cluster =
    immune_cluster[
      excluded_excluded_conflict
    ],
  immune_type =
    immune_proposed_type[
      excluded_excluded_conflict
    ],
  immune_state =
    immune_proposed_state[
      excluded_excluded_conflict
    ],
  hepato_cluster =
    hepato_cluster[
      excluded_excluded_conflict
    ],
  hepato_type =
    hepato_proposed_type[
      excluded_excluded_conflict
    ],
  hepato_state =
    hepato_proposed_state[
      excluded_excluded_conflict
    ],
  stringsAsFactors = FALSE
)

write.csv(
  conflict,
  file.path(
    outdir,
    "08_ambiguous_crosslineage_conflicts.csv"
  ),
  row.names = FALSE
)

conflict_combo <- as.data.frame(
  table(
    immune_type =
      conflict$immune_type,
    hepato_type =
      conflict$hepato_type
  ),
  stringsAsFactors = FALSE
)

conflict_combo <- conflict_combo[
  conflict_combo$Freq > 0,
  ,
  drop = FALSE
]

write.csv(
  conflict_combo,
  file.path(
    outdir,
    "09_ambiguous_crosslineage_conflict_combinations.csv"
  ),
  row.names = FALSE
)

conflict_by_sample <- as.data.frame(
  table(
    sample =
      conflict$sample
  ),
  stringsAsFactors = FALSE
)

names(
  conflict_by_sample
)[2] <- "n_cells"

conflict_by_sample <- conflict_by_sample[
  conflict_by_sample$n_cells > 0,
  ,
  drop = FALSE
]

write.csv(
  conflict_by_sample,
  file.path(
    outdir,
    "10_ambiguous_crosslineage_by_sample.csv"
  ),
  row.names = FALSE
)

cat(
  "Excluded-excluded same-lineage consensus cells:",
  sum(
    excluded_excluded_agree
  ),
  "\n"
)

cat(
  "Excluded-excluded conflicting cells retained as Ambiguous_crosslineage:",
  sum(
    excluded_excluded_conflict
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 9. Build final whole-liver annotation
# ------------------------------------------------------------

final_type <- old_celltype
final_state <- old_celltype

annotation_source <- rep(
  "v1.3_carried_forward",
  length(
    canonical_cells
  )
)

annotation_detail <- rep(
  "Validated v1.3 structural annotation carried forward",
  length(
    canonical_cells
  )
)

# Immune-only cells
idx <- which(
  immune_member &
    !hepato_member
)

final_type[
  idx
] <- immune_proposed_type[
  idx
]

final_state[
  idx
] <- immune_proposed_state[
  idx
]

annotation_source[
  idx
] <- ifelse(
  immune_keep[
    idx
  ],
  "immune_v1.3.2A.1_clean",
  "immune_v1.3.2A.1_crosslineage_correction"
)

annotation_detail[
  idx
] <- paste0(
  "immune cluster ",
  immune_cluster[
    idx
  ],
  "; ",
  immune_proposal_kind[
    idx
  ]
)

# Hepatobiliary-only cells
idx <- which(
  hepato_member &
    !immune_member
)

final_type[
  idx
] <- hepato_proposed_type[
  idx
]

final_state[
  idx
] <- hepato_proposed_state[
  idx
]

annotation_source[
  idx
] <- ifelse(
  hepato_keep[
    idx
  ],
  "hepatobiliary_v1.3.2B.1_clean",
  "hepatobiliary_v1.3.2B.1_crosslineage_correction"
)

annotation_detail[
  idx
] <- paste0(
  "hepatobiliary cluster ",
  hepato_cluster[
    idx
  ],
  "; ",
  hepato_proposal_kind[
    idx
  ]
)

# Overlap: immune clean wins over hepatobiliary excluded
idx <- which(
  overlap &
    immune_keep &
    !hepato_keep
)

final_type[
  idx
] <- immune_proposed_type[
  idx
]

final_state[
  idx
] <- immune_proposed_state[
  idx
]

annotation_source[
  idx
] <- "immune_clean_over_hepatobiliary_exclusion"

annotation_detail[
  idx
] <- paste0(
  "immune cluster ",
  immune_cluster[
    idx
  ],
  " clean; hepatobiliary cluster ",
  hepato_cluster[
    idx
  ],
  " excluded"
)

# Overlap: hepatobiliary clean wins over immune excluded
idx <- which(
  overlap &
    !immune_keep &
    hepato_keep
)

final_type[
  idx
] <- hepato_proposed_type[
  idx
]

final_state[
  idx
] <- hepato_proposed_state[
  idx
]

annotation_source[
  idx
] <- "hepatobiliary_clean_over_immune_exclusion"

annotation_detail[
  idx
] <- paste0(
  "hepatobiliary cluster ",
  hepato_cluster[
    idx
  ],
  " clean; immune cluster ",
  immune_cluster[
    idx
  ],
  " excluded"
)

# Overlap: both excluded and AGREE on alternate lineage
idx <- which(
  excluded_excluded_agree
)

if (
  length(
    idx
  ) > 0
) {

  final_type[
    idx
  ] <- immune_proposed_type[
    idx
  ]

  final_state[
    idx
  ] <- paste0(
    "Crosslineage_consensus_",
    immune_proposed_type[
      idx
    ]
  )

  annotation_source[
    idx
  ] <- "crosslineage_consensus_correction"

  annotation_detail[
    idx
  ] <- paste0(
    "immune=",
    immune_proposed_state[
      idx
    ],
    "; hepatobiliary=",
    hepato_proposed_state[
      idx
    ]
  )
}

# Overlap: both excluded but DISAGREE on alternate lineage
# Conservative freeze rule:
#   retain the cells in the canonical object, but do not force
#   them into either lineage.
idx <- which(
  excluded_excluded_conflict
)

if (
  length(
    idx
  ) > 0
) {

  final_type[
    idx
  ] <- "Ambiguous_crosslineage"

  final_state[
    idx
  ] <- "Doublet_or_crosslineage_ambiguous"

  annotation_source[
    idx
  ] <- "refinement_conflict_ambiguous"

  annotation_detail[
    idx
  ] <- paste0(
    "immune excluded -> ",
    immune_proposed_type[
      idx
    ],
    " [",
    immune_proposed_state[
      idx
    ],
    "]; hepatobiliary excluded -> ",
    hepato_proposed_type[
      idx
    ],
    " [",
    hepato_proposed_state[
      idx
    ],
    "]"
  )
}

# ------------------------------------------------------------
# 10. Sanity checks on final labels
# ------------------------------------------------------------

if (
  anyNA(
    final_type
  ) ||
    anyNA(
      final_state
    ) ||
    anyNA(
      annotation_source
    )
) {

  stop(
    "NA present in final v1.3.3 annotation."
  )
}

obsolete_labels <- c(
  "Inflammatory_monocyte_macrophage",
  "SPP1_TREM2_macrophage",
  "Hepatoblast_or_immature_hepatocyte",
  "Cholangiocyte_or_ductal"
)

obsolete_remaining <- intersect(
  unique(
    final_type
  ),
  obsolete_labels
)

if (
  length(
    obsolete_remaining
  ) > 0
) {

  stop(
    "Obsolete broad labels remain after merge: ",
    paste(
      obsolete_remaining,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 11. Carry source-specific refinement metadata into canonical
# ------------------------------------------------------------

obj$immune_candidate_v132A1 <-
  immune_member

obj$immune_cluster_v132A <-
  immune_cluster

obj$immune_state_v132A1 <-
  immune_manual

obj$immune_keep_clean_v132A1 <-
  immune_keep

obj$hepatobiliary_candidate_v132B1 <-
  hepato_member

obj$hepatobiliary_cluster_v132B <-
  hepato_cluster

obj$hepatobiliary_state_v132B1 <-
  hepato_state

obj$hepatobiliary_keep_clean_v132B1 <-
  hepato_keep

obj$hepatobiliary_manual_confidence_v132B1 <-
  hepato_confidence

# ------------------------------------------------------------
# 12. Outside-core ductal audit flag
# ------------------------------------------------------------

outside_core_ductal_flag <- rep(
  FALSE,
  length(
    canonical_cells
  )
)

outside_core_ductal_count <- NA_integer_

if (
  file.exists(
    outside_core_ductal_csv
  )
) {

  outside_core <- read.csv(
    outside_core_ductal_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (
    "cell" %in%
      colnames(
        outside_core
      )
  ) {

    outside_core_ductal_count <- nrow(
      outside_core
    )

    outside_core_ductal_flag[
      canonical_cells %in%
        outside_core$cell
    ] <- TRUE
  }
}

obj$ductal_candidate_audit_only_v133 <-
  outside_core_ductal_flag

# ------------------------------------------------------------
# 13. Strong single-sample state-bias flag
# ------------------------------------------------------------

state_sample_bias_flag <- rep(
  FALSE,
  length(
    canonical_cells
  )
)

state_sample_bias_flag[
  hepato_member &
    hepato_cluster %in%
      c(
        "9",
        "12",
        "14"
      ) &
    hepato_keep
] <- TRUE

obj$state_sample_bias_ge70_v133 <-
  state_sample_bias_flag

obj$annotation_ambiguous_crosslineage_v133 <-
  excluded_excluded_conflict

# ------------------------------------------------------------
# 14. Final frozen metadata
# ------------------------------------------------------------

obj$celltype_final_v133 <-
  final_type

obj$cellstate_final_v133 <-
  final_state

obj$annotation_source_v133 <-
  annotation_source

obj$annotation_detail_v133 <-
  annotation_detail

obj$annotation_version_v133 <-
  "v1.3.3"

obj$annotation_frozen_v133 <-
  TRUE

# Factor order for stable figures/tables
preferred_type_order <- c(
  "Hepatocyte",
  "Cholangiocyte",
  "Ductular_progenitor",
  "LSEC",
  "Endothelial",
  "HSC_quiescent",
  "HSC_fibroblast_mesenchymal",
  "Fibroblast_mesenchymal",
  "Mesenchymal_perivascular",
  "Mesothelial",
  "Macrophage",
  "Monocyte",
  "Dendritic_cell",
  "T_cell",
  "NK_cell",
  "NK_gamma_delta_T",
  "B_cell",
  "Plasma_cell",
  "Mast_cell",
  "Neutrophil",
  "Cycling_immune",
  "Cycling_cell",
  "Erythroid",
  "Schwann_neural_like",
  "Ambiguous_crosslineage"
)

observed_final_types <- unique(
  obj$celltype_final_v133
)

extra_types <- setdiff(
  observed_final_types,
  preferred_type_order
)

final_type_levels <- c(
  preferred_type_order[
    preferred_type_order %in%
      observed_final_types
  ],
  sort(
    extra_types
  )
)

obj$celltype_final_v133 <- factor(
  obj$celltype_final_v133,
  levels =
    final_type_levels
)

# ------------------------------------------------------------
# 15. Whole-liver counts / crosswalk / provenance tables
# ------------------------------------------------------------

final_counts <- as.data.frame(
  table(
    celltype_final_v133 =
      obj$celltype_final_v133
  ),
  stringsAsFactors = FALSE
)

names(
  final_counts
)[2] <- "n_cells"

final_counts$percent <- 100 *
  final_counts$n_cells /
  sum(
    final_counts$n_cells
  )

final_counts <- final_counts[
  final_counts$n_cells > 0,
  ,
  drop = FALSE
]

write.csv(
  final_counts,
  file.path(
    outdir,
    "02_final_celltype_counts.csv"
  ),
  row.names = FALSE
)

old_new_crosswalk <- as.data.frame(
  table(
    old_v13 =
      old_celltype,
    final_v133 =
      obj$celltype_final_v133
  ),
  stringsAsFactors = FALSE
)

old_new_crosswalk <- old_new_crosswalk[
  old_new_crosswalk$Freq > 0,
  ,
  drop = FALSE
]

write.csv(
  old_new_crosswalk,
  file.path(
    outdir,
    "03_v1.3_to_v1.3.3_annotation_crosswalk.csv"
  ),
  row.names = FALSE
)

source_counts <- as.data.frame(
  table(
    annotation_source =
      obj$annotation_source_v133
  ),
  stringsAsFactors = FALSE
)

names(
  source_counts
)[2] <- "n_cells"

write.csv(
  source_counts,
  file.path(
    outdir,
    "04_annotation_source_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Final type x sample QC
# ------------------------------------------------------------

sample_order <- sample_order_fun(
  obj[[sample_col]][, 1]
)

type_sample <- as.data.frame(
  table(
    celltype_final_v133 =
      obj$celltype_final_v133,
    sample =
      obj[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

type_sample$type_total <- ave(
  type_sample$Freq,
  type_sample$celltype_final_v133,
  FUN = sum
)

type_sample$sample_total <- ave(
  type_sample$Freq,
  type_sample$sample,
  FUN = sum
)

type_sample$within_type_pct <- ifelse(
  type_sample$type_total > 0,
  100 *
    type_sample$Freq /
    type_sample$type_total,
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
  levels =
    sample_order
)

write.csv(
  type_sample,
  file.path(
    outdir,
    "05_final_celltype_by_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. State counts
# ------------------------------------------------------------

state_counts <- as.data.frame(
  table(
    cellstate_final_v133 =
      obj$cellstate_final_v133
  ),
  stringsAsFactors = FALSE
)

names(
  state_counts
)[2] <- "n_cells"

state_counts <- state_counts[
  state_counts$n_cells > 0,
  ,
  drop = FALSE
]

state_counts <- state_counts[
  order(
    state_counts$n_cells,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  state_counts,
  file.path(
    outdir,
    "06_final_cellstate_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. Cell-level final annotation export
# ------------------------------------------------------------

cell_metadata <- data.frame(
  cell =
    canonical_cells,
  sample =
    obj[[sample_col]][, 1],
  NAS_score =
    obj[[nas_col]][, 1],
  NAS_group =
    obj[[nasgrp_col]][, 1],
  old_celltype_v13 =
    old_celltype,
  final_celltype_v133 =
    as.character(
      obj$celltype_final_v133
    ),
  final_cellstate_v133 =
    obj$cellstate_final_v133,
  annotation_source_v133 =
    obj$annotation_source_v133,
  annotation_detail_v133 =
    obj$annotation_detail_v133,
  state_sample_bias_ge70_v133 =
    obj$state_sample_bias_ge70_v133,
  annotation_ambiguous_crosslineage_v133 =
    obj$annotation_ambiguous_crosslineage_v133,
  ductal_candidate_audit_only_v133 =
    obj$ductal_candidate_audit_only_v133,
  stringsAsFactors = FALSE
)

write.csv(
  cell_metadata,
  file.path(
    outdir,
    "07_v1.3.3_final_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 19. Final high-contrast palette
# ------------------------------------------------------------

final_palette <- c(
  "Hepatocyte" =
    "#00A6A6",
  "Cholangiocyte" =
    "#00D34D",
  "Ductular_progenitor" =
    "#76FF03",
  "LSEC" =
    "#00D4FF",
  "Endothelial" =
    "#1E90FF",
  "HSC_quiescent" =
    "#FF66B3",
  "HSC_fibroblast_mesenchymal" =
    "#E6007E",
  "Fibroblast_mesenchymal" =
    "#B0005A",
  "Mesenchymal_perivascular" =
    "#FF3FA4",
  "Mesothelial" =
    "#FFB000",
  "Macrophage" =
    "#D7191C",
  "Monocyte" =
    "#FF5A36",
  "Dendritic_cell" =
    "#7B2CBF",
  "T_cell" =
    "#6A3D9A",
  "NK_cell" =
    "#1F78B4",
  "NK_gamma_delta_T" =
    "#4D7CFE",
  "B_cell" =
    "#00A650",
  "Plasma_cell" =
    "#B5179E",
  "Mast_cell" =
    "#FF1493",
  "Neutrophil" =
    "#7CB342",
  "Cycling_immune" =
    "#FFD600",
  "Cycling_cell" =
    "#C28A00",
  "Erythroid" =
    "#8B1A1A",
  "Schwann_neural_like" =
    "#8C6BB1",
  "Ambiguous_crosslineage" =
    "#4D4D4D"
)

missing_palette <- setdiff(
  final_type_levels,
  names(
    final_palette
  )
)

if (
  length(
    missing_palette
  ) > 0
) {

  fallback_cols <- grDevices::hcl.colors(
    length(
      missing_palette
    ),
    palette = "Dark 3"
  )

  names(
    fallback_cols
  ) <- missing_palette

  final_palette <- c(
    final_palette,
    fallback_cols
  )
}

# ------------------------------------------------------------
# 20. Sample palette
# ------------------------------------------------------------

sample_palette <- c(
  "LT1" = "#E41A1C",
  "LT2" = "#377EB8",
  "LT3" = "#4DAF4A",
  "LT4" = "#984EA3",
  "LT5" = "#FF7F00",
  "LT6" = "#00A6A6",
  "LT7" = "#A65628",
  "LT8" = "#F781BF"
)

# ------------------------------------------------------------
# 21. Improved UMAPs
#
# Final policy:
#   vector PDF, raster = FALSE, cairo_pdf when available.
# ------------------------------------------------------------

cat("=== FINAL v1.3.3 UMAP GENERATION ===\n")

p_native <- DimPlot(
  obj,
  reduction =
    native_umap,
  group.by =
    "celltype_final_v133",
  pt.size =
    overview_pt_size,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    final_palette
) +
  ggtitle(
    "Final whole-liver annotation v1.3.3 - native UMAP"
  ) +
  theme_classic(
    base_size = 11
  )

p_rpca <- DimPlot(
  obj,
  reduction =
    rpca_umap,
  group.by =
    "celltype_final_v133",
  pt.size =
    overview_pt_size,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    final_palette
) +
  ggtitle(
    "Final whole-liver annotation v1.3.3 - RPCA UMAP"
  ) +
  theme_classic(
    base_size = 11
  )

p01 <- p_native +
  p_rpca +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

save_pdf(
  file.path(
    outdir,
    "Fig01_FinalAnnotation_native_vs_RPCA_UMAP.pdf"
  ),
  p01,
  width = 18,
  height = 8.5
)

p02 <- DimPlot(
  obj,
  reduction =
    rpca_umap,
  group.by =
    sample_col,
  pt.size =
    overview_pt_size,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = FALSE,
  cols =
    sample_palette
) +
  ggtitle(
    "LT sample distribution - final canonical v1.3.3 RPCA UMAP"
  ) +
  theme_classic(
    base_size = 11
  )

save_pdf(
  file.path(
    outdir,
    "Fig02_Sample_RPCA_UMAP.pdf"
  ),
  p02,
  width = 10.5,
  height = 8.5
)

p03 <- DimPlot(
  obj,
  reduction =
    rpca_umap,
  group.by =
    "celltype_final_v133",
  split.by =
    sample_col,
  ncol = 4,
  pt.size =
    split_pt_size,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = FALSE,
  cols =
    final_palette
)

save_pdf(
  file.path(
    outdir,
    "Fig03_FinalCelltype_RPCA_UMAP_split_by_sample.pdf"
  ),
  p03,
  width = 20,
  height = 13
)

# Refined state overview (supplemental QC)
p04 <- DimPlot(
  obj,
  reduction =
    rpca_umap,
  group.by =
    "cellstate_final_v133",
  pt.size = 0.16,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = TRUE,
  repel = TRUE,
  label.size = 2.5
) +
  ggtitle(
    "Final refined cell-state annotation v1.3.3 - RPCA UMAP"
  ) +
  theme_classic(
    base_size = 10
  )

save_pdf(
  file.path(
    outdir,
    "Fig04_FinalCellState_RPCA_UMAP.pdf"
  ),
  p04,
  width = 15,
  height = 10
)

# ------------------------------------------------------------
# 22. Final type x sample QC heatmap
# ------------------------------------------------------------

p05 <- ggplot(
  type_sample,
  aes(
    x = sample,
    y =
      celltype_final_v133,
    fill =
      within_type_pct
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
          within_type_pct
        ),
        ""
      )
    ),
    size = 2.0
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Sample contribution within final v1.3.3 cell types",
    x = NULL,
    y = NULL,
    fill = "%"
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

save_pdf(
  file.path(
    outdir,
    "Fig05_FinalCelltype_SampleContribution_QC.pdf"
  ),
  p05,
  width = 14,
  height = 11
)

# ------------------------------------------------------------
# 23. Annotation-source UMAP
# ------------------------------------------------------------

p06 <- DimPlot(
  obj,
  reduction =
    rpca_umap,
  group.by =
    "annotation_source_v133",
  pt.size = 0.16,
  raster = FALSE,
  shuffle = TRUE,
  seed = seed,
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  ggtitle(
    "Annotation provenance - canonical v1.3.3"
  ) +
  theme_classic(
    base_size = 10
  )

save_pdf(
  file.path(
    outdir,
    "Fig06_AnnotationSource_RPCA_UMAP.pdf"
  ),
  p06,
  width = 12,
  height = 9
)

# ------------------------------------------------------------
# 24. Final review flags
# ------------------------------------------------------------

final_review_flags <- c(
  "Human MASH v1.3.3 Whole-liver Annotation Freeze",
  "================================================",
  "",
  paste0(
    "Canonical cells retained: ",
    ncol(
      obj
    )
  ),
  paste0(
    "Final broad cell types: ",
    length(
      unique(
        as.character(
          obj$celltype_final_v133
        )
      )
    )
  ),
  paste0(
    "Final refined states: ",
    length(
      unique(
        obj$cellstate_final_v133
      )
    )
  ),
  paste0(
    "Immune candidate cells merged: ",
    sum(
      immune_member
    )
  ),
  paste0(
    "Hepatobiliary candidate cells merged: ",
    sum(
      hepato_member
    )
  ),
  paste0(
    "Immune/hepatobiliary overlap cells audited: ",
    sum(
      overlap
    )
  ),
  paste0(
    "Excluded-excluded same-lineage consensus cells: ",
    sum(
      excluded_excluded_agree
    )
  ),
  paste0(
    "Ambiguous cross-lineage conflict cells: ",
    sum(
      excluded_excluded_conflict
    ),
    " (",
    round(
      100 *
        sum(
          excluded_excluded_conflict
        ) /
        ncol(
          obj
        ),
      3
    ),
    "% of canonical)"
  ),
  paste0(
    "Strong sample-bias state flag cells: ",
    sum(
      obj$state_sample_bias_ge70_v133
    )
  ),
  paste0(
    "Outside-core ductal audit-only cells: ",
    ifelse(
      is.na(
        outside_core_ductal_count
      ),
      "NA",
      outside_core_ductal_count
    )
  ),
  "",
  "FREEZE DECISIONS:",
  "  Immune refinement v1.3.2A.1 integrated.",
  "  Hepatobiliary refinement v1.3.2B.1 integrated.",
  "  HSC / LSEC / Endothelial and other validated structural labels carried forward unless corrected by a refinement.",
  "  Old Hepatoblast_or_immature_hepatocyte label is retired.",
  "  Old broad inflammatory/SPP1-TREM2 macrophage labels are retired in favor of refined immune lineages/states.",
  "  Cells excluded by BOTH refinements with DIFFERENT alternate-lineage calls are retained as Ambiguous_crosslineage.",
  "  These ambiguous cells are not assigned to a lineage-specific downstream analysis unless separately resolved.",
  "  All original v1.3 annotation columns are preserved for provenance.",
  "",
  "UMAP POLICY:",
  paste0(
    "  Overview point size = ",
    overview_pt_size
  ),
  paste0(
    "  Split point size = ",
    split_pt_size
  ),
  "  raster = FALSE",
  paste0(
    "  cairo_pdf available = ",
    capabilities(
      "cairo"
    )
  ),
  "  high-contrast saturated palette used",
  "",
  "CAUTION FLAGS:",
  paste0(
    "  Ambiguous_crosslineage cells retained: ",
    sum(
      excluded_excluded_conflict
    ),
    ". Exclude from lineage-specific biological inference unless separately resolved."
  ),
  "  Hepatobiliary clusters 9, 12, 14 remain annotated states but carry >=70% single-sample-bias flags.",
  "  Outside-core strict ductal candidates remain audit-only and were not automatically rescued.",
  "",
  "STATUS:",
  "  WHOLE_LIVER_ANNOTATION_FROZEN_v1.3.3",
  "  Original canonical v1.3 object was NOT modified."
)

writeLines(
  final_review_flags,
  file.path(
    outdir,
    "18_v1.3.3_review_flags.txt"
  )
)

cat(
  paste(
    final_review_flags,
    collapse = "\n"
  )
)

cat("\n\n")

# ------------------------------------------------------------
# 25. Save new canonical v1.3.3 RDS
# ------------------------------------------------------------

saveRDS(
  obj,
  final_rds
)

if (
  !file.exists(
    final_rds
  )
) {

  stop(
    "Failed to write canonical v1.3.3 RDS."
  )
}

# Re-read lightweight validation of object header / metadata
obj_check <- readRDS(
  final_rds
)

if (
  ncol(
    obj_check
  ) != expected_canonical_cells
) {

  stop(
    "Saved v1.3.3 object cell-count validation failed."
  )
}

required_final_meta <- c(
  "celltype_final_v133",
  "cellstate_final_v133",
  "annotation_source_v133",
  "annotation_version_v133",
  "annotation_frozen_v133",
  "state_sample_bias_ge70_v133",
  "annotation_ambiguous_crosslineage_v133",
  "ductal_candidate_audit_only_v133"
)

missing_final_meta <- setdiff(
  required_final_meta,
  colnames(
    obj_check@meta.data
  )
)

if (
  length(
    missing_final_meta
  ) > 0
) {

  stop(
    "Saved v1.3.3 object missing required final metadata: ",
    paste(
      missing_final_meta,
      collapse = ", "
    )
  )
}

rm(
  obj_check
)

gc()

# ------------------------------------------------------------
# 26. Manifest / provenance
# ------------------------------------------------------------

manifest <- c(
  "version=v1.3.3",
  paste0(
    "input_canonical_v1.3=",
    canonical_v13_rds
  ),
  paste0(
    "input_canonical_v1.3_md5=",
    unname(
      tools::md5sum(
        canonical_v13_rds
      )
    )
  ),
  paste0(
    "immune_metadata=",
    immune_metadata_csv
  ),
  paste0(
    "immune_metadata_md5=",
    unname(
      tools::md5sum(
        immune_metadata_csv
      )
    )
  ),
  paste0(
    "hepatobiliary_metadata=",
    hepato_metadata_csv
  ),
  paste0(
    "hepatobiliary_metadata_md5=",
    unname(
      tools::md5sum(
        hepato_metadata_csv
      )
    )
  ),
  paste0(
    "output_canonical_v1.3.3=",
    final_rds
  ),
  paste0(
    "output_canonical_v1.3.3_md5=",
    unname(
      tools::md5sum(
        final_rds
      )
    )
  ),
  paste0(
    "cells=",
    ncol(
      obj
    )
  ),
  paste0(
    "final_celltypes=",
    length(
      unique(
        as.character(
          obj$celltype_final_v133
        )
      )
    )
  ),
  paste0(
    "final_states=",
    length(
      unique(
        obj$cellstate_final_v133
      )
    )
  ),
  paste0(
    "immune_hepatobiliary_overlap_cells=",
    sum(
      overlap
    )
  ),
  paste0(
    "excluded_excluded_same_lineage_consensus_cells=",
    sum(
      excluded_excluded_agree
    )
  ),
  paste0(
    "ambiguous_crosslineage_cells=",
    sum(
      excluded_excluded_conflict
    )
  ),
  "excluded_excluded_conflict_policy=retain_as_Ambiguous_crosslineage",
  paste0(
    "cairo_pdf_available=",
    capabilities(
      "cairo"
    )
  ),
  "original_v1.3_modified=FALSE",
  "all_81348_cells_retained=TRUE",
  "annotation_status=WHOLE_LIVER_ANNOTATION_FROZEN_v1.3.3"
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

# ------------------------------------------------------------
# 27. Final terminal summary
# ------------------------------------------------------------

cat("============================================================\n")
cat("v1.3.3 COMPLETE\n")
cat("WHOLE-LIVER ANNOTATION FREEZE COMPLETE\n")
cat("ALL 81,348 CELLS RETAINED\n")
cat("IMMUNE v1.3.2A.1 MERGED\n")
cat("HEPATOBILIARY v1.3.2B.1 MERGED\n")
cat(
  "AMBIGUOUS CROSS-LINEAGE CELLS RETAINED: ",
  sum(
    excluded_excluded_conflict
  ),
  "\n",
  sep = ""
)
cat("FINAL VECTOR UMAPs GENERATED\n")
cat("ORIGINAL canonical v1.3 OBJECT NOT MODIFIED\n")
cat("\n")
cat("New canonical v1.3.3 object:\n")
cat(final_rds, "\n")
cat("\n")
cat("Output directory:\n")
cat(outdir, "\n")
cat("============================================================\n")
