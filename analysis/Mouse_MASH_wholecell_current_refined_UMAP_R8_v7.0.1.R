#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(7000)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Whole-cell UMAP update using current frozen lineage-specific annotations
#
# Version: v7.0.1
#
# PURPOSE
#   Update the mouse whole-liver UMAP before starting the human analysis.
#
#   Layer 1:
#     Frozen broad whole-cell annotation (R8 high-saturation colors)
#
#   Layer 2:
#     Transfer current frozen lineage-specific states back to the SAME
#     frozen whole-cell UMAP coordinates by exact cell-barcode matching.
#
#   Detailed lineages:
#     1) Kupffer/Macrophage
#     2) HSC/Mesenchymal
#     3) Hepatocyte
#     4) LSEC
#     5) Cholangiocyte
#     6) Monocyte
#
# IMPORTANT DESIGN
#   - NO whole-cell reclustering
#   - NO reintegration
#   - NO new UMAP
#   - NO source RDS overwrite
#   - Existing frozen whole-cell umapRPCA coordinates are reused.
#   - Lineage-specific annotations are transferred by exact barcode.
#   - Whole-cell broad lineage remains authoritative for lineage ownership.
#
#   In particular:
#   - Clean-B Mphi was originally derived from the macrophage compartment that
#     included Kupffer/Macrophage + Monocyte. For this whole-cell display,
#     Mphi subtype labels are transferred ONLY to cells whose frozen whole-cell
#     broad parent is Kupffer_Macrophage.
#   - The new v6.9.4 Monocyte annotation owns cells whose broad parent is
#     Monocyte. This prevents Mphi/Monocyte annotation collision.
#
# DISPLAY POLICY
#   - R8-like high-saturation palette.
#   - Refined states use lineage-family color sets.
#   - Cells belonging to a targeted broad lineage but absent from its clean
#     refined object remain visible as "<lineage> | parent-only/unresolved".
#   - QC-watch / sample-biased states remain visible because this is an
#     annotation-display update, not a primary statistical analysis.
#
# EXPECTED WHOLE-CELL PARENT
#   /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/
#   Mouse_MASH_RDS/WholeCell_Layer1_ParentFreeze_v5.1.1/RDS/
#   Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds
#
# OUTPUT
#   /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/
#   Mouse_MASH_RDS/WholeCell_RefinedAnnotation_UMAP_v7.0.1/
# ==============================================================================


# ==============================================================================
# 1. User settings
# ==============================================================================

VERSION <- "v7.0.1"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_Layer1_ParentFreeze_v5.1.1",
  "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "WholeCell_RefinedAnnotation_UMAP_v7.0.1"
)

FIG_OUT <- file.path(
  OUT,
  "Figures"
)

TAB_OUT <- file.path(
  OUT,
  "Tables"
)

OBJ_OUT <- file.path(
  OUT,
  "Objects"
)

LOG_OUT <- file.path(
  OUT,
  "Logs"
)

for (
  d in c(
    OUT,
    FIG_OUT,
    TAB_OUT,
    OBJ_OUT,
    LOG_OUT
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Keep FALSE for routine visualization.
# TRUE writes a full 104,588-cell Seurat RDS with the new metadata columns,
# which can be large and slow to save.
SAVE_UPDATED_WHOLE_RDS <- FALSE

WHOLE_UMAP <- "umapRPCA"

WHOLE_LAYER1_CANDIDATES <- c(
  "wholecell_layer1_FINAL_v511",
  "celltype_for_R8plot_FIXED2",
  "celltype_for_R8plot"
)

WHOLE_SAMPLE_CANDIDATES <- c(
  "sample_for_annotation",
  "sample_FIXED2",
  "sample",
  "orig.ident"
)

WHOLE_CONDITION_CANDIDATES <- c(
  "condition_FIXED2",
  "condition_v502",
  "condition",
  "sample_4group"
)

SAMPLE_LEVELS <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

CONDITION_LEVELS <- c(
  "STD",
  "CDAHFD",
  "Sham",
  "Tx"
)

# Plot settings
PT_ALL <- 0.42
PT_SPLIT <- 0.34
PT_HIGHLIGHT <- 0.38

PDF_WIDTH_ALL <- 14
PDF_HEIGHT_ALL <- 10

PDF_WIDTH_SPLIT <- 18
PDF_HEIGHT_SPLIT <- 11

LEGEND_NCOL_BROAD <- 2
LEGEND_NCOL_REFINED <- 3


# ==============================================================================
# 2. Helpers
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


save_png <- function(
  p,
  file,
  width,
  height,
  dpi = 400
) {

  ggsave(
    filename = file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )
}


first_existing <- function(
  candidates,
  available
) {

  hit <- candidates[
    candidates %in% available
  ]

  if (!length(hit)) {
    return(NA_character_)
  }

  hit[[1]]
}


canonical_sample <- function(x) {

  x <- as.character(x)

  out <- rep(
    NA_character_,
    length(x)
  )

  z <- toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      trimws(x)
    )
  )

  out[
    grepl("^STD", z)
  ] <- "STD_rep1"

  out[
    grepl("^(CDHFD|CDAHFD)", z)
  ] <- "CDHFD_rep1"

  out[
    grepl("^SHAM1$", z)
  ] <- "Sham1"

  out[
    grepl("^SHAM20$", z)
  ] <- "Sham20"

  out[
    grepl("^TX17$", z)
  ] <- "Tx17"

  out[
    grepl("^TX5$", z)
  ] <- "Tx5"

  # relaxed matching after exact canonical rules
  out[
    is.na(out) &
      grepl("SHAM.*1", z)
  ] <- "Sham1"

  out[
    is.na(out) &
      grepl("SHAM.*20", z)
  ] <- "Sham20"

  out[
    is.na(out) &
      grepl("TX.*17", z)
  ] <- "Tx17"

  out[
    is.na(out) &
      grepl("TX.*5", z)
  ] <- "Tx5"

  out
}


canonical_condition <- function(
  condition,
  sample
) {

  condition <- as.character(condition)
  sample <- as.character(sample)

  out <- rep(
    NA_character_,
    length(sample)
  )

  zc <- toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      trimws(condition)
    )
  )

  zs <- toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      trimws(sample)
    )
  )

  out[
    grepl("^STD", zc) |
      grepl("^STD", zs)
  ] <- "STD"

  out[
    grepl("^(CDHFD|CDAHFD)", zc) |
      grepl("^(CDHFD|CDAHFD)", zs)
  ] <- "CDAHFD"

  out[
    grepl("^SHAM", zc) |
      grepl("^SHAM", zs)
  ] <- "Sham"

  out[
    grepl("^TX", zc) |
      grepl("^TX", zs)
  ] <- "Tx"

  out
}


canonical_broad <- function(x) {

  x <- as.character(x)

  z <- trimws(x)

  out <- z

  out[
    z %in% c(
      "HSC/Mesenchymal",
      "HSC_Mesenchymal",
      "HSC",
      "Mesenchymal"
    )
  ] <- "HSC_Mesenchymal"

  out[
    z %in% c(
      "Kupffer/Macrophage",
      "Kupffer_Macrophage",
      "Kupffer macrophage",
      "Kupffer_macrophage"
    )
  ] <- "Kupffer_Macrophage"

  out[
    z %in% c(
      "Vascular endothelial",
      "Vascular_endothelial",
      "Endothelial (Vascular)"
    )
  ] <- "Vascular_endothelial"

  out[
    z %in% c(
      "Dendritic",
      "Dendritic_cell"
    )
  ] <- "Dendritic_cell"

  out[
    z %in% c(
      "B",
      "B cell",
      "B_cell"
    )
  ] <- "B_cell"

  out[
    z %in% c(
      "T",
      "T cell",
      "T_cell"
    )
  ] <- "T_cell"

  out[
    z %in% c(
      "NK",
      "NK cell",
      "NK_cell"
    )
  ] <- "NK_cell"

  out[
    z %in% c(
      "Plasma",
      "Plasma cell",
      "Plasma_cell"
    )
  ] <- "Plasma_cell"

  out[
    z %in% c(
      "Biliary",
      "Cholangiocyte"
    )
  ] <- "Cholangiocyte"

  out[
    z %in% c(
      "LSEC",
      "Liver sinusoidal endothelial"
    )
  ] <- "LSEC"

  out[
    z %in% c(
      "Monocyte"
    )
  ] <- "Monocyte"

  out[
    z %in% c(
      "Neutrophil"
    )
  ] <- "Neutrophil"

  out[
    z %in% c(
      "Hepatocyte"
    )
  ] <- "Hepatocyte"

  out[
    z %in% c(
      "Mesothelial"
    )
  ] <- "Mesothelial"

  out[
    z %in% c(
      "RBC",
      "Erythroid"
    )
  ] <- "RBC"

  out[
    z %in% c(
      "Platelet"
    )
  ] <- "Platelet"

  out[
    z %in% c(
      "Cycling"
    )
  ] <- "Cycling"

  out
}


resolve_state_col <- function(
  object,
  preferred,
  expected_states,
  source_name
) {

  md <- object@meta.data

  direct <- preferred[
    preferred %in%
      colnames(md)
  ]

  if (length(direct)) {

    for (col in direct) {

      vals <- as.character(
        md[[col]]
      )

      overlap <- sum(
        unique(
          vals[
            !is.na(vals)
          ]
        ) %in%
          expected_states
      )

      frac <- mean(
        vals %in%
          expected_states,
        na.rm = TRUE
      )

      if (
        overlap >=
          min(
            2L,
            length(expected_states)
          ) &&
        is.finite(frac) &&
        frac >= 0.50
      ) {
        return(col)
      }
    }
  }

  # Fallback: inspect all character/factor metadata columns.
  candidate_cols <- colnames(md)[
    vapply(
      md,
      function(x) {
        is.character(x) ||
          is.factor(x)
      },
      logical(1)
    )
  ]

  if (!length(candidate_cols)) {
    stop(
      "No character/factor metadata columns in ",
      source_name
    )
  }

  score <- lapply(
    candidate_cols,
    function(col) {

      vals <- as.character(
        md[[col]]
      )

      present <- unique(
        vals[
          !is.na(vals)
        ]
      )

      overlap <- sum(
        present %in%
          expected_states
      )

      frac <- mean(
        vals %in%
          expected_states,
        na.rm = TRUE
      )

      data.frame(
        column = col,
        overlap = overlap,
        fraction = frac,
        stringsAsFactors = FALSE
      )
    }
  )

  score <- do.call(
    rbind,
    score
  )

  score <- score[
    order(
      -score$overlap,
      -score$fraction
    ),
    ,
    drop = FALSE
  ]

  best <- score[
    1,
    ,
    drop = FALSE
  ]

  if (
    nrow(best) != 1 ||
    best$overlap <
      min(
        2L,
        length(expected_states)
      ) ||
    !is.finite(
      best$fraction
    ) ||
    best$fraction < 0.50
  ) {

    stop(
      "Could not resolve state column for ",
      source_name,
      ". Best candidate: ",
      paste(
        capture.output(
          print(
            head(
              score,
              10
            ),
            row.names = FALSE
          )
        ),
        collapse = "\n"
      )
    )
  }

  best$column[[1]]
}


candidate_rds_files <- function(
  explicit_paths,
  search_roots,
  filename_pattern
) {

  explicit_hits <- explicit_paths[
    file.exists(
      explicit_paths
    )
  ]

  search_hits <- character()

  for (root_now in search_roots) {

    if (!dir.exists(root_now)) {
      next
    }

    hits <- list.files(
      root_now,
      pattern = filename_pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    search_hits <- c(
      search_hits,
      hits
    )
  }

  unique(
    c(
      explicit_hits,
      search_hits
    )
  )
}


load_annotation_source <- function(
  source_name,
  explicit_paths,
  search_roots,
  filename_pattern,
  preferred_cols,
  expected_states
) {

  files <- candidate_rds_files(
    explicit_paths = explicit_paths,
    search_roots = search_roots,
    filename_pattern = filename_pattern
  )

  if (!length(files)) {
    stop(
      "No candidate RDS files found for ",
      source_name
    )
  }

  audit <- list()

  for (i in seq_along(files)) {

    f <- files[[i]]

    msg(
      "Inspecting ",
      source_name,
      " candidate: ",
      f
    )

    object_now <- tryCatch(
      readRDS(f),
      error = function(e) {
        NULL
      }
    )

    if (
      is.null(object_now) ||
      !inherits(
        object_now,
        "Seurat"
      )
    ) {

      audit[[i]] <- data.frame(
        file = f,
        usable = FALSE,
        state_col = NA_character_,
        matched_states = 0L,
        matched_fraction = NA_real_,
        n_cells = NA_integer_,
        stringsAsFactors = FALSE
      )

      next
    }

    col_now <- tryCatch(
      resolve_state_col(
        object = object_now,
        preferred = preferred_cols,
        expected_states = expected_states,
        source_name = source_name
      ),
      error = function(e) {
        NA_character_
      }
    )

    if (is.na(col_now)) {

      audit[[i]] <- data.frame(
        file = f,
        usable = FALSE,
        state_col = NA_character_,
        matched_states = 0L,
        matched_fraction = NA_real_,
        n_cells = ncol(object_now),
        stringsAsFactors = FALSE
      )

      rm(object_now)
      gc()

      next
    }

    vals <- as.character(
      object_now@meta.data[[col_now]]
    )

    matched_states <- sum(
      unique(
        vals[
          !is.na(vals)
        ]
      ) %in%
        expected_states
    )

    matched_fraction <- mean(
      vals %in%
        expected_states,
      na.rm = TRUE
    )

    audit[[i]] <- data.frame(
      file = f,
      usable = TRUE,
      state_col = col_now,
      matched_states = matched_states,
      matched_fraction = matched_fraction,
      n_cells = ncol(object_now),
      stringsAsFactors = FALSE
    )

    # Return the first fully usable explicit/candidate object.
    if (
      matched_states >=
        min(
          2L,
          length(expected_states)
        ) &&
      matched_fraction >= 0.50
    ) {

      audit_df <- do.call(
        rbind,
        audit
      )

      return(
        list(
          object = object_now,
          path = f,
          state_col = col_now,
          audit = audit_df
        )
      )
    }

    rm(object_now)
    gc()
  }

  audit_df <- do.call(
    rbind,
    audit
  )

  write.csv(
    audit_df,
    file.path(
      LOG_OUT,
      paste0(
        source_name,
        "_candidate_RDS_audit_",
        VERSION,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  stop(
    "No usable annotation RDS identified for ",
    source_name,
    ". See audit in ",
    LOG_OUT
  )
}


build_transfer_map <- function(
  source_object,
  state_col,
  expected_states,
  whole_cells,
  whole_broad,
  target_parent,
  source_name,
  display_map
) {

  state_raw <- as.character(
    source_object@meta.data[[state_col]]
  )

  source_cells <- colnames(
    source_object
  )

  keep <- !is.na(
    state_raw
  ) &
    state_raw %in%
      expected_states

  source_cells <- source_cells[
    keep
  ]

  state_raw <- state_raw[
    keep
  ]

  whole_idx <- match(
    source_cells,
    whole_cells
  )

  exact_match <- !is.na(
    whole_idx
  )

  parent_at_match <- rep(
    NA_character_,
    length(source_cells)
  )

  parent_at_match[
    exact_match
  ] <- whole_broad[
    whole_idx[
      exact_match
    ]
  ]

  parent_consistent <-
    exact_match &
    parent_at_match ==
      target_parent

  state_display <- unname(
    display_map[
      state_raw
    ]
  )

  if (
    any(
      is.na(
        state_display
      )
    )
  ) {

    bad_states <- unique(
      state_raw[
        is.na(
          state_display
        )
      ]
    )

    stop(
      "Missing display-map entries for ",
      source_name,
      ": ",
      paste(
        bad_states,
        collapse = ", "
      )
    )
  }

  data.frame(
    cell = source_cells,
    source = source_name,
    source_state = state_raw,
    display_state = state_display,
    exact_match_in_whole = exact_match,
    whole_parent = parent_at_match,
    target_parent = target_parent,
    parent_consistent = parent_consistent,
    stringsAsFactors = FALSE
  )
}


add_legend_style <- function(
  p,
  ncol_legend,
  base_size = 9
) {

  p +
    theme_classic(
      base_size = base_size
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 7.5
      ),
      legend.key.height = grid::unit(
        0.38,
        "cm"
      ),
      legend.key.width = grid::unit(
        0.38,
        "cm"
      )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 2.2,
          alpha = 1
        ),
        ncol = ncol_legend
      )
    )
}


# ==============================================================================
# 3. Frozen annotation definitions
# ==============================================================================

# ------------------------------------------------------------------------------
# 3A. Mphi
# ------------------------------------------------------------------------------

MPHI_STATES <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Other"
)

MPHI_DISPLAY <- c(
  "Anti-inflammatory-Mphi" =
    "Mphi | Anti-inflammatory",
  "Inflammatory-Mphi" =
    "Mphi | Inflammatory",
  "ECM-associated inflammatory-Mphi" =
    "Mphi | ECM-associated inflammatory",
  "Repair/Resolution-Mphi" =
    "Mphi | Repair/Resolution",
  "Lipid-associated/TREM2-Mphi" =
    "Mphi | Lipid-associated/TREM2",
  "Other" =
    "Mphi | Other"
)

MPHI_COL_CANDIDATES <- c(
  "macrophage_class_Res2_FINAL_v4145_char",
  "macrophage_class_Res2_FINAL_v4145",
  "macrophage_class_Res2_v484",
  "macrophage_class"
)

MPHI_EXPLICIT <- c(
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


# ------------------------------------------------------------------------------
# 3B. HSC
# ------------------------------------------------------------------------------

HSC_STATES <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

HSC_DISPLAY <- c(
  "qHSC" =
    "HSC | qHSC",
  "ECM-activated HSC" =
    "HSC | ECM-activated",
  "Contractile HSC" =
    "HSC | Contractile"
)

HSC_COL_CANDIDATES <- c(
  "HSC_state3_v612",
  "HSC_state3_v611",
  "HSC_state3",
  "HSC_state"
)

HSC_EXPLICIT <- c(
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "HSC_RefinedAnnotation_v6.1.2",
    "RDS",
    "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
  )
)


# ------------------------------------------------------------------------------
# 3C. Hepatocyte
# ------------------------------------------------------------------------------

HEP_STATES <- c(
  "Periportal_Hepatocyte_1",
  "Injury_inflammatory_Hepatocyte",
  "Pericentral_Hepatocyte",
  "MT_high_QC_Hepatocyte",
  "Periportal_Hepatocyte_2",
  "Intermediate_Hepatocyte",
  "Cycling_G2M_Hepatocyte",
  "Cycling_S_Hepatocyte"
)

HEP_DISPLAY <- c(
  "Periportal_Hepatocyte_1" =
    "Hep | Periportal-1",
  "Injury_inflammatory_Hepatocyte" =
    "Hep | Injury-inflammatory",
  "Pericentral_Hepatocyte" =
    "Hep | Pericentral",
  "MT_high_QC_Hepatocyte" =
    "Hep | MT-high/QC",
  "Periportal_Hepatocyte_2" =
    "Hep | Periportal-2",
  "Intermediate_Hepatocyte" =
    "Hep | Intermediate",
  "Cycling_G2M_Hepatocyte" =
    "Hep | Cycling-G2M",
  "Cycling_S_Hepatocyte" =
    "Hep | Cycling-S"
)

HEP_COL_CANDIDATES <- c(
  "hepatocyte_state_FINAL_v652",
  "hepatocyte_state_v652",
  "hepatocyte_state_FINAL",
  "hepatocyte_state"
)

HEP_EXPLICIT <- c(
  file.path(
    ROOT,
    "Mouse_MASH_Hepatocyte",
    "Hepatocyte_FINAL_state_pseudobulk_v6.5.2",
    "RDS",
    "Mouse_MASH_Hepatocyte_FINAL_annotated_v6.5.2.rds"
  )
)


# ------------------------------------------------------------------------------
# 3D. LSEC
# ------------------------------------------------------------------------------

LSEC_STATES <- c(
  "Homeostatic_like_LSEC",
  "Wnt_angiocrine_high_LSEC",
  "Inflammatory_stress_high_LSEC",
  "Low_quality_ambient_enriched_LSEC",
  "Tx5_enriched_LSEC_state",
  "Cycling_LSEC",
  "Tx17_enriched_Cd209b_Ctsj_LSEC"
)

LSEC_DISPLAY <- c(
  "Homeostatic_like_LSEC" =
    "LSEC | Homeostatic-like",
  "Wnt_angiocrine_high_LSEC" =
    "LSEC | Wnt/angiocrine-high",
  "Inflammatory_stress_high_LSEC" =
    "LSEC | Inflammatory/stress-high",
  "Low_quality_ambient_enriched_LSEC" =
    "LSEC | Low-quality/ambient",
  "Tx5_enriched_LSEC_state" =
    "LSEC | Tx5-enriched",
  "Cycling_LSEC" =
    "LSEC | Cycling",
  "Tx17_enriched_Cd209b_Ctsj_LSEC" =
    "LSEC | Tx17-enriched Cd209b/Ctsj"
)

LSEC_COL_CANDIDATES <- c(
  "LSEC_state_v675",
  "LSEC_state_FINAL_v675",
  "LSEC_state_frozen_v675",
  "lsec_state_v675",
  "LSEC_state",
  "lsec_state"
)

LSEC_EXPLICIT <- c(
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_LSEC_v6.7.5",
    "objects",
    "Mouse_MASH_LSEC_annotation_frozen_v6.7.5.rds"
  ),
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_LSEC_v6.7.5",
    "objects",
    "Mouse_MASH_LSEC_annotated_v6.7.5.rds"
  ),
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_LSEC_v6.7.5",
    "RDS",
    "Mouse_MASH_LSEC_annotated_v6.7.5.rds"
  )
)


# ------------------------------------------------------------------------------
# 3E. Cholangiocyte
# ------------------------------------------------------------------------------

CHOL_STATES <- c(
  "Pcp4l1_Serpina6_homeostatic_like",
  "Kptn_Pacs2_homeostatic_like",
  "Ccl2_Vcam1_inflammatory_reactive",
  "Msln_Aqp5_ductular_like",
  "Disease_enriched_high_complexity_QC_watch",
  "Cycling_cholangiocyte",
  "Ciliated_cholangiocyte",
  "IEG_stress_response",
  "Krt20_Cdh17_reactive_epithelial",
  "Dmbt1_Duox2_reactive_epithelial",
  "Tuft_like_cholangiocyte"
)

CHOL_DISPLAY <- c(
  "Pcp4l1_Serpina6_homeostatic_like" =
    "Chol | Pcp4l1/Serpina6 homeostatic",
  "Kptn_Pacs2_homeostatic_like" =
    "Chol | Kptn/Pacs2 homeostatic",
  "Ccl2_Vcam1_inflammatory_reactive" =
    "Chol | Ccl2/Vcam1 inflammatory-reactive",
  "Msln_Aqp5_ductular_like" =
    "Chol | Msln/Aqp5 ductular-like",
  "Disease_enriched_high_complexity_QC_watch" =
    "Chol | Disease-enriched QC-watch",
  "Cycling_cholangiocyte" =
    "Chol | Cycling",
  "Ciliated_cholangiocyte" =
    "Chol | Ciliated",
  "IEG_stress_response" =
    "Chol | IEG/stress-response",
  "Krt20_Cdh17_reactive_epithelial" =
    "Chol | Krt20/Cdh17 reactive",
  "Dmbt1_Duox2_reactive_epithelial" =
    "Chol | Dmbt1/Duox2 reactive",
  "Tuft_like_cholangiocyte" =
    "Chol | Tuft-like"
)

CHOL_COL_CANDIDATES <- c(
  "Cholangiocyte_state_frozen_v6.8.6",
  "Cholangiocyte_state_v6.8.6",
  "Cholangiocyte_state_v686",
  "cholangiocyte_state_v686",
  "Cholangiocyte_state",
  "cholangiocyte_state"
)

CHOL_EXPLICIT <- c(
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_Cholangiocyte_v6.8.6",
    "objects",
    "Mouse_MASH_Cholangiocyte_annotation_frozen_v6.8.6.rds"
  ),
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_Cholangiocyte_v6.8.6",
    "objects",
    "Mouse_MASH_Cholangiocyte_annotated_v6.8.6.rds"
  )
)


# ------------------------------------------------------------------------------
# 3F. Monocyte
# ------------------------------------------------------------------------------

MONO_STATES <- c(
  "S100a8_S100a9_Thbs1_stress_inflammatory_Monocyte",
  "Mmp8_Sell_Chil3_Vcan_classical_inflammatory_Monocyte",
  "Pald1_C3ar1_homeostatic_like_Monocyte",
  "Tnf_Il1rn_Olr1_Gpnmb_inflammatory_remodeling_Monocyte",
  "Cd300e_Pglyrp1_Cd36_S1pr5_activated_Monocyte",
  "Adamdec1_Pecam1_low_complexity_state",
  "Ms4a7_Mmp12_Dab2_C1q_monocyte_to_macrophage_transition",
  "Nos2_Cxcl9_Saa3_IFNg_inflammatory_Monocyte",
  "Hepatocyte_RNA_high_Monocyte_QC_watch",
  "Ifit_Rsad2_Cmpk2_IFN_responsive_Monocyte"
)

MONO_DISPLAY <- c(
  "S100a8_S100a9_Thbs1_stress_inflammatory_Monocyte" =
    "Mono | S100a8/S100a9/Thbs1 stress-inflammatory",
  "Mmp8_Sell_Chil3_Vcan_classical_inflammatory_Monocyte" =
    "Mono | Mmp8/Sell/Chil3/Vcan classical-inflammatory",
  "Pald1_C3ar1_homeostatic_like_Monocyte" =
    "Mono | Pald1/C3ar1 homeostatic-like",
  "Tnf_Il1rn_Olr1_Gpnmb_inflammatory_remodeling_Monocyte" =
    "Mono | Tnf/Il1rn/Olr1/Gpnmb inflammatory-remodeling",
  "Cd300e_Pglyrp1_Cd36_S1pr5_activated_Monocyte" =
    "Mono | Cd300e/Pglyrp1/Cd36/S1pr5 activated",
  "Adamdec1_Pecam1_low_complexity_state" =
    "Mono | Adamdec1/Pecam1 sample-biased",
  "Ms4a7_Mmp12_Dab2_C1q_monocyte_to_macrophage_transition" =
    "Mono | Ms4a7/Mmp12/Dab2/C1q transition-like",
  "Nos2_Cxcl9_Saa3_IFNg_inflammatory_Monocyte" =
    "Mono | Nos2/Cxcl9/Saa3 IFNg-inflammatory",
  "Hepatocyte_RNA_high_Monocyte_QC_watch" =
    "Mono | Hepatocyte-RNA-high QC-watch",
  "Ifit_Rsad2_Cmpk2_IFN_responsive_Monocyte" =
    "Mono | Ifit/Rsad2/Cmpk2 IFN-responsive"
)

MONO_COL_CANDIDATES <- c(
  "Monocyte_state_frozen_v6.9.4",
  "Monocyte_state_v6.9.4",
  "Monocyte_state_v694",
  "monocyte_state_v694",
  "Monocyte_state",
  "monocyte_state"
)

MONO_EXPLICIT <- c(
  file.path(
    ROOT,
    "Mouse_MASH_RDS",
    "Mouse_MASH_Monocyte_v6.9.4",
    "objects",
    "Mouse_MASH_Monocyte_annotation_frozen_v6.9.4.rds"
  )
)


# ==============================================================================
# 4. R8 broad and refined palettes
# ==============================================================================

# Broad R8-like high-saturation palette.
# Direction follows the previously preferred whole-cell scheme:
# LSEC = vivid cyan; Hepatocyte = blue-green; Cholangiocyte = vivid green;
# HSC = vivid pink; macrophage/monocyte and immune lineages strongly separated.

R8_BROAD <- c(
  "B_cell" = "#FF6B6B",
  "Cholangiocyte" = "#00C853",
  "Cycling" = "#FFB000",
  "Dendritic_cell" = "#A6D854",
  "Hepatocyte" = "#00A087",
  "HSC_Mesenchymal" = "#FF4FA3",
  "Kupffer_Macrophage" = "#F04444",
  "LSEC" = "#00C8FF",
  "Mesothelial" = "#00BFA5",
  "Monocyte" = "#2F65FF",
  "Neutrophil" = "#0066FF",
  "NK_cell" = "#5E3CFF",
  "Plasma_cell" = "#A020F0",
  "Platelet" = "#C77CFF",
  "RBC" = "#E83E9B",
  "T_cell" = "#FF1493",
  "Vascular_endothelial" = "#FF5C8A",
  "Other" = "#B8B8B8"
)


# Mphi detailed colors retain the established high-saturation state colors.
R8_MPHI <- c(
  "Mphi | Anti-inflammatory" = "#00C8FF",
  "Mphi | Inflammatory" = "#FF3131",
  "Mphi | ECM-associated inflammatory" = "#FF8C00",
  "Mphi | Repair/Resolution" = "#A020F0",
  "Mphi | Lipid-associated/TREM2" = "#0066FF",
  "Mphi | Other" = "#8A8A8A",
  "Mphi | parent-only/unresolved" = "#F7A1A1"
)


R8_HSC <- c(
  "HSC | qHSC" = "#FF76B6",
  "HSC | ECM-activated" = "#FF2D95",
  "HSC | Contractile" = "#C2185B",
  "HSC | parent-only/unresolved" = "#FFB6D9"
)


R8_HEP <- c(
  "Hep | Periportal-1" = "#006D5B",
  "Hep | Injury-inflammatory" = "#00BFA5",
  "Hep | Pericentral" = "#00A087",
  "Hep | MT-high/QC" = "#26C6DA",
  "Hep | Periportal-2" = "#008C72",
  "Hep | Intermediate" = "#00C9A7",
  "Hep | Cycling-G2M" = "#00E5FF",
  "Hep | Cycling-S" = "#1DE9B6",
  "Hep | parent-only/unresolved" = "#66CDBB"
)


R8_LSEC <- c(
  "LSEC | Homeostatic-like" = "#00BFC4",
  "LSEC | Wnt/angiocrine-high" = "#00E5FF",
  "LSEC | Inflammatory/stress-high" = "#0077FF",
  "LSEC | Low-quality/ambient" = "#7FDBFF",
  "LSEC | Tx5-enriched" = "#00A6FB",
  "LSEC | Cycling" = "#4CC9F0",
  "LSEC | Tx17-enriched Cd209b/Ctsj" = "#0066FF",
  "LSEC | parent-only/unresolved" = "#7FDBFF"
)


R8_CHOL <- c(
  "Chol | Pcp4l1/Serpina6 homeostatic" = "#007A33",
  "Chol | Kptn/Pacs2 homeostatic" = "#00A651",
  "Chol | Ccl2/Vcam1 inflammatory-reactive" = "#00C853",
  "Chol | Msln/Aqp5 ductular-like" = "#64DD17",
  "Chol | Disease-enriched QC-watch" = "#AEEA00",
  "Chol | Cycling" = "#76FF03",
  "Chol | Ciliated" = "#00E676",
  "Chol | IEG/stress-response" = "#C6FF00",
  "Chol | Krt20/Cdh17 reactive" = "#2ECC71",
  "Chol | Dmbt1/Duox2 reactive" = "#B2FF59",
  "Chol | Tuft-like" = "#00B248",
  "Chol | parent-only/unresolved" = "#8CE99A"
)


R8_MONO <- c(
  "Mono | S100a8/S100a9/Thbs1 stress-inflammatory" = "#0057FF",
  "Mono | Mmp8/Sell/Chil3/Vcan classical-inflammatory" = "#1E88E5",
  "Mono | Pald1/C3ar1 homeostatic-like" = "#00A6FB",
  "Mono | Tnf/Il1rn/Olr1/Gpnmb inflammatory-remodeling" = "#3A86FF",
  "Mono | Cd300e/Pglyrp1/Cd36/S1pr5 activated" = "#4361EE",
  "Mono | Adamdec1/Pecam1 sample-biased" = "#5E60CE",
  "Mono | Ms4a7/Mmp12/Dab2/C1q transition-like" = "#6930C3",
  "Mono | Nos2/Cxcl9/Saa3 IFNg-inflammatory" = "#4D96FF",
  "Mono | Hepatocyte-RNA-high QC-watch" = "#0096C7",
  "Mono | Ifit/Rsad2/Cmpk2 IFN-responsive" = "#2F65FF",
  "Mono | parent-only/unresolved" = "#8DB3FF"
)


R8_REFINED <- c(
  R8_MPHI,
  R8_HSC,
  R8_HEP,
  R8_LSEC,
  R8_CHOL,
  R8_MONO,

  # Broad lineages that are not refined in v7.0.0.
  R8_BROAD[
    c(
      "B_cell",
      "Cycling",
      "Dendritic_cell",
      "Mesothelial",
      "Neutrophil",
      "NK_cell",
      "Plasma_cell",
      "Platelet",
      "RBC",
      "T_cell",
      "Vascular_endothelial",
      "Other"
    )
  ]
)


# ==============================================================================
# 5. Load frozen whole-cell parent
# ==============================================================================

if (!file.exists(WHOLE_RDS)) {
  stop(
    "Whole-cell frozen RDS not found: ",
    WHOLE_RDS
  )
}

msg(
  "Loading frozen whole-cell parent..."
)

whole <- readRDS(
  WHOLE_RDS
)

if (!inherits(
  whole,
  "Seurat"
)) {
  stop(
    "WHOLE_RDS is not a Seurat object."
  )
}

if (
  !WHOLE_UMAP %in%
    Reductions(
      whole
    )
) {
  stop(
    "Required whole-cell UMAP reduction not found: ",
    WHOLE_UMAP,
    "\nAvailable reductions: ",
    paste(
      Reductions(
        whole
      ),
      collapse = ", "
    )
  )
}

layer1_col <- first_existing(
  WHOLE_LAYER1_CANDIDATES,
  colnames(
    whole@meta.data
  )
)

if (is.na(layer1_col)) {
  stop(
    "Could not resolve whole-cell Layer1 column. Available metadata:\n",
    paste(
      colnames(
        whole@meta.data
      ),
      collapse = ", "
    )
  )
}

sample_col <- first_existing(
  WHOLE_SAMPLE_CANDIDATES,
  colnames(
    whole@meta.data
  )
)

if (is.na(sample_col)) {
  stop(
    "Could not resolve whole-cell sample column."
  )
}

condition_col <- first_existing(
  WHOLE_CONDITION_CANDIDATES,
  colnames(
    whole@meta.data
  )
)

msg(
  "Whole cells: ",
  ncol(whole)
)

msg(
  "Whole features: ",
  nrow(whole)
)

msg(
  "Layer1 source: ",
  layer1_col
)

msg(
  "Sample source: ",
  sample_col
)

msg(
  "Condition source: ",
  ifelse(
    is.na(condition_col),
    "<derived from sample>",
    condition_col
  )
)


# ==============================================================================
# 6. Standardize whole-cell plotting metadata
# ==============================================================================

whole$wholecell_layer1_R8_v701 <-
  canonical_broad(
    whole@meta.data[[
      layer1_col
    ]]
  )

whole$sample_R8_v701 <-
  canonical_sample(
    whole@meta.data[[
      sample_col
    ]]
  )

condition_raw <- if (
  is.na(
    condition_col
  )
) {
  rep(
    NA_character_,
    ncol(whole)
  )
} else {
  whole@meta.data[[
    condition_col
  ]]
}

whole$condition_R8_v701 <-
  canonical_condition(
    condition = condition_raw,
    sample = whole@meta.data[[
      sample_col
    ]]
  )

whole$sample_R8_v701 <-
  factor(
    whole$sample_R8_v701,
    levels = SAMPLE_LEVELS
  )

whole$condition_R8_v701 <-
  factor(
    whole$condition_R8_v701,
    levels = CONDITION_LEVELS
  )

broad_present <- sort(
  unique(
    as.character(
      whole$wholecell_layer1_R8_v701
    )
  )
)

msg(
  "Broad labels present: ",
  paste(
    broad_present,
    collapse = ", "
  )
)

missing_broad_colors <- setdiff(
  broad_present,
  names(
    R8_BROAD
  )
)

if (length(missing_broad_colors)) {

  warning(
    "Broad labels without predefined R8 color: ",
    paste(
      missing_broad_colors,
      collapse = ", "
    ),
    ". They will be assigned fallback colors."
  )

  fallback <- grDevices::hcl.colors(
    length(
      missing_broad_colors
    ),
    palette = "Dynamic"
  )

  names(
    fallback
  ) <- missing_broad_colors

  R8_BROAD <- c(
    R8_BROAD,
    fallback
  )
}


# ==============================================================================
# 7. Load current frozen lineage-specific annotations
# ==============================================================================

source_specs <- list(

  Mphi = list(
    explicit = MPHI_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_Mphi_RDS"
      )
    ),
    pattern =
      "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4\\.14\\.5\\.[rR][dD][sS]$",
    preferred_cols = MPHI_COL_CANDIDATES,
    states = MPHI_STATES,
    parent = "Kupffer_Macrophage",
    display = MPHI_DISPLAY
  ),

  HSC = list(
    explicit = HSC_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_RDS",
        "HSC_RefinedAnnotation_v6.1.2"
      )
    ),
    pattern =
      "v6\\.1\\.2.*\\.[rR][dD][sS]$",
    preferred_cols = HSC_COL_CANDIDATES,
    states = HSC_STATES,
    parent = "HSC_Mesenchymal",
    display = HSC_DISPLAY
  ),

  Hepatocyte = list(
    explicit = HEP_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_Hepatocyte",
        "Hepatocyte_FINAL_state_pseudobulk_v6.5.2"
      )
    ),
    pattern =
      "v6\\.5\\.2.*\\.[rR][dD][sS]$",
    preferred_cols = HEP_COL_CANDIDATES,
    states = HEP_STATES,
    parent = "Hepatocyte",
    display = HEP_DISPLAY
  ),

  LSEC = list(
    explicit = LSEC_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_RDS"
      )
    ),
    pattern =
      "v6\\.7\\.5.*\\.[rR][dD][sS]$",
    preferred_cols = LSEC_COL_CANDIDATES,
    states = LSEC_STATES,
    parent = "LSEC",
    display = LSEC_DISPLAY
  ),

  Cholangiocyte = list(
    explicit = CHOL_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_RDS"
      )
    ),
    pattern =
      "v6\\.8\\.6.*\\.[rR][dD][sS]$",
    preferred_cols = CHOL_COL_CANDIDATES,
    states = CHOL_STATES,
    parent = "Cholangiocyte",
    display = CHOL_DISPLAY
  ),

  Monocyte = list(
    explicit = MONO_EXPLICIT,
    search_roots = c(
      file.path(
        ROOT,
        "Mouse_MASH_RDS",
        "Mouse_MASH_Monocyte_v6.9.4"
      )
    ),
    pattern =
      "v6\\.9\\.4.*\\.[rR][dD][sS]$",
    preferred_cols = MONO_COL_CANDIDATES,
    states = MONO_STATES,
    parent = "Monocyte",
    display = MONO_DISPLAY
  )
)


source_audit <- list()
transfer_maps <- list()

whole_cells <- colnames(
  whole
)

whole_broad <- as.character(
  whole$wholecell_layer1_R8_v701
)


for (
  source_name in names(
    source_specs
  )
) {

  spec <- source_specs[[
    source_name
  ]]

  msg(
    "Resolving annotation source: ",
    source_name
  )

  src <- load_annotation_source(
    source_name = source_name,
    explicit_paths = spec$explicit,
    search_roots = spec$search_roots,
    filename_pattern = spec$pattern,
    preferred_cols = spec$preferred_cols,
    expected_states = spec$states
  )

  state_col_now <- src$state_col

  map_now <- build_transfer_map(
    source_object = src$object,
    state_col = state_col_now,
    expected_states = spec$states,
    whole_cells = whole_cells,
    whole_broad = whole_broad,
    target_parent = spec$parent,
    source_name = source_name,
    display_map = spec$display
  )

  source_audit[[
    source_name
  ]] <- data.frame(
    source = source_name,
    RDS = src$path,
    state_column = state_col_now,
    target_parent = spec$parent,
    source_object_cells = ncol(
      src$object
    ),
    annotated_source_cells = nrow(
      map_now
    ),
    exact_whole_matches = sum(
      map_now$exact_match_in_whole
    ),
    exact_match_fraction = mean(
      map_now$exact_match_in_whole
    ),
    parent_consistent_matches = sum(
      map_now$parent_consistent
    ),
    parent_consistent_fraction_of_exact =
      sum(
        map_now$parent_consistent
      ) /
      max(
        1,
        sum(
          map_now$exact_match_in_whole
        )
      ),
    cross_parent_exact_matches = sum(
      map_now$exact_match_in_whole &
        !map_now$parent_consistent
    ),
    stringsAsFactors = FALSE
  )

  transfer_maps[[
    source_name
  ]] <- map_now

  write.csv(
    src$audit,
    file.path(
      LOG_OUT,
      paste0(
        source_name,
        "_source_resolution_audit_",
        VERSION,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  rm(src)
  gc()
}


source_audit_df <- do.call(
  rbind,
  source_audit
)

rownames(
  source_audit_df
) <- NULL

write.csv(
  source_audit_df,
  file.path(
    TAB_OUT,
    paste0(
      "01_annotation_source_and_barcode_audit_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

msg(
  "Annotation-source audit:"
)

print(
  source_audit_df
)


# Fail on poor barcode compatibility.
if (
  any(
    source_audit_df$exact_match_fraction <
      0.90,
    na.rm = TRUE
  )
) {

  stop(
    "At least one lineage-specific annotation object has <90% exact barcode ",
    "compatibility with the whole-cell parent. Review ",
    file.path(
      TAB_OUT,
      paste0(
        "01_annotation_source_and_barcode_audit_",
        VERSION,
        ".csv"
      )
    )
  )
}


# ==============================================================================
# 8. Audit raw cross-source overlap
#
# v7.0.1 policy
#   Dedicated refined lineage objects are allowed to override the older broad
#   parent label when exact cell barcodes are present.
#
#   The only expected large cross-source overlap is Mphi vs Monocyte because the
#   historical Clean-B macrophage object included the Monocyte compartment.
#
#   Ownership priority:
#     Monocyte v6.9.4 > Mphi v4.14.5 for Monocyte-source cells.
#
#   Mphi labels are therefore used only for current broad Kupffer_Macrophage
#   cells that are NOT present in the frozen Monocyte v6.9.4 object.
#
#   HSC / Hepatocyte / LSEC / Cholangiocyte / Monocyte dedicated objects own
#   their exact barcodes even if the older whole-cell broad label differs.
# ==============================================================================

exact_cells_by_source <- lapply(
  transfer_maps,
  function(x) {
    x$cell[
      x$exact_match_in_whole
    ]
  }
)

source_names <- names(
  exact_cells_by_source
)

pairwise_overlap <- list()
k_overlap <- 0L

for (i in seq_along(source_names)) {

  if (i == length(source_names)) {
    next
  }

  for (j in seq.int(i + 1L, length(source_names))) {

    a <- source_names[[i]]
    b <- source_names[[j]]

    ov <- intersect(
      exact_cells_by_source[[a]],
      exact_cells_by_source[[b]]
    )

    k_overlap <- k_overlap + 1L

    pairwise_overlap[[k_overlap]] <- data.frame(
      source_A = a,
      source_B = b,
      n_overlap = length(ov),
      stringsAsFactors = FALSE
    )
  }
}

pairwise_overlap_df <- do.call(
  rbind,
  pairwise_overlap
)

write.csv(
  pairwise_overlap_df,
  file.path(
    TAB_OUT,
    paste0(
      "02_raw_cross_source_overlap_audit_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

unexpected_overlap <- pairwise_overlap_df[
  pairwise_overlap_df$n_overlap > 0 &
    !(
      pairwise_overlap_df$source_A %in% c("Mphi", "Monocyte") &
        pairwise_overlap_df$source_B %in% c("Mphi", "Monocyte")
    ),
  ,
  drop = FALSE
]

if (nrow(unexpected_overlap) > 0) {

  stop(
    "Unexpected overlap exists between dedicated refined lineage sources. ",
    "Review 02_raw_cross_source_overlap_audit_",
    VERSION,
    ".csv"
  )
}

msg(
  "Raw cross-source overlap audit:"
)

print(
  pairwise_overlap_df
)


# ==============================================================================
# 9. Build updated whole-cell refined annotation with explicit ownership
# ==============================================================================

whole$refined_lineage_v701 <-
  NA_character_

whole$refined_state_raw_v701 <-
  NA_character_

whole$refined_source_v701 <-
  NA_character_

whole$refined_ownership_rule_v701 <-
  "broad_parent_only"

whole$annotation_refined_v701 <-
  as.character(
    whole$wholecell_layer1_R8_v701
  )


# Parent-remainder labels preserve cells excluded during lineage cleanup.
# These are shown in the same lineage color family rather than neutral gray.
parent_only_labels <- c(
  "Kupffer_Macrophage" =
    "Mphi | parent-only/unresolved",
  "HSC_Mesenchymal" =
    "HSC | parent-only/unresolved",
  "Hepatocyte" =
    "Hep | parent-only/unresolved",
  "LSEC" =
    "LSEC | parent-only/unresolved",
  "Cholangiocyte" =
    "Chol | parent-only/unresolved",
  "Monocyte" =
    "Mono | parent-only/unresolved"
)

for (parent_now in names(parent_only_labels)) {

  idx <- which(
    whole$wholecell_layer1_R8_v701 ==
      parent_now
  )

  whole$annotation_refined_v701[
    idx
  ] <- parent_only_labels[[
    parent_now
  ]]
}


# ------------------------------------------------------------------
# 9A. Dedicated non-Mphi lineage sources own their exact barcodes.
#     Monocyte is applied last and therefore has highest myeloid priority.
# ------------------------------------------------------------------

dedicated_order <- c(
  "HSC",
  "Hepatocyte",
  "LSEC",
  "Cholangiocyte"
)

for (source_name in dedicated_order) {

  map_now <- transfer_maps[[
    source_name
  ]]

  map_now <- map_now[
    map_now$exact_match_in_whole,
    ,
    drop = FALSE
  ]

  idx <- match(
    map_now$cell,
    colnames(whole)
  )

  if (anyNA(idx)) {
    stop(
      "Internal exact-barcode transfer error for ",
      source_name
    )
  }

  whole$refined_lineage_v701[
    idx
  ] <- source_name

  whole$refined_state_raw_v701[
    idx
  ] <- map_now$source_state

  whole$refined_source_v701[
    idx
  ] <- source_name

  whole$refined_ownership_rule_v701[
    idx
  ] <- "dedicated_refined_source_exact_barcode"

  whole$annotation_refined_v701[
    idx
  ] <- map_now$display_state
}


# ------------------------------------------------------------------
# 9B. Mphi ownership
#
# Historical Clean-B Mphi contains Monocyte-compartment cells.
# Use Mphi subtype labels only for current broad Kupffer_Macrophage cells
# that are not present in the frozen v6.9.4 Monocyte object.
# ------------------------------------------------------------------

mono_exact_cells <- exact_cells_by_source[[
  "Monocyte"
]]

mphi_map <- transfer_maps[[
  "Mphi"
]]

mphi_map <- mphi_map[
  mphi_map$exact_match_in_whole,
  ,
  drop = FALSE
]

mphi_idx_whole <- match(
  mphi_map$cell,
  colnames(whole)
)

mphi_keep <-
  whole_broad[
    mphi_idx_whole
  ] ==
    "Kupffer_Macrophage" &
  !mphi_map$cell %in%
    mono_exact_cells

mphi_map_use <- mphi_map[
  mphi_keep,
  ,
  drop = FALSE
]

idx <- match(
  mphi_map_use$cell,
  colnames(whole)
)

whole$refined_lineage_v701[
  idx
] <- "Mphi"

whole$refined_state_raw_v701[
  idx
] <- mphi_map_use$source_state

whole$refined_source_v701[
  idx
] <- "Mphi"

whole$refined_ownership_rule_v701[
  idx
] <- "Mphi_only_broad_Kupffer_not_in_Monocyte_v6.9.4"

whole$annotation_refined_v701[
  idx
] <- mphi_map_use$display_state


# ------------------------------------------------------------------
# 9C. Monocyte ownership
#
# The dedicated v6.9.4 Monocyte object is newer than the whole-parent
# boundary and therefore owns ALL of its exact barcodes, including cells
# whose older v5.1.1 broad label differs.
# ------------------------------------------------------------------

mono_map <- transfer_maps[[
  "Monocyte"
]]

mono_map <- mono_map[
  mono_map$exact_match_in_whole,
  ,
  drop = FALSE
]

idx <- match(
  mono_map$cell,
  colnames(whole)
)

if (anyNA(idx)) {
  stop(
    "Internal exact-barcode transfer error for Monocyte"
  )
}

whole$refined_lineage_v701[
  idx
] <- "Monocyte"

whole$refined_state_raw_v701[
  idx
] <- mono_map$source_state

whole$refined_source_v701[
  idx
] <- "Monocyte"

whole$refined_ownership_rule_v701[
  idx
] <- "Monocyte_v6.9.4_exact_barcode_priority"

whole$annotation_refined_v701[
  idx
] <- mono_map$display_state


# ------------------------------------------------------------------
# 9D. Ownership audit
# ------------------------------------------------------------------

ownership_audit <- lapply(
  names(transfer_maps),
  function(source_name) {

    map_now <- transfer_maps[[
      source_name
    ]]

    exact_now <- map_now[
      map_now$exact_match_in_whole,
      ,
      drop = FALSE
    ]

    owned_cells <- colnames(whole)[
      whole$refined_source_v701 ==
        source_name
    ]

    exact_owned <- intersect(
      exact_now$cell,
      owned_cells
    )

    exact_cross_parent_owned <- exact_now[
      exact_now$cell %in%
        exact_owned &
      !exact_now$parent_consistent,
      ,
      drop = FALSE
    ]

    data.frame(
      source = source_name,
      exact_source_cells = nrow(exact_now),
      final_owned_cells = length(exact_owned),
      cross_parent_cells_owned =
        nrow(exact_cross_parent_owned),
      stringsAsFactors = FALSE
    )
  }
)

ownership_audit <- do.call(
  rbind,
  ownership_audit
)

write.csv(
  ownership_audit,
  file.path(
    TAB_OUT,
    paste0(
      "02b_final_lineage_ownership_audit_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

msg(
  "Final lineage ownership audit:"
)

print(
  ownership_audit
)


# ==============================================================================
# 10. Refined annotation audit tables
# ==============================================================================

annotation_counts <- as.data.frame(
  table(
    annotation_refined =
      whole$annotation_refined_v701
  ),
  stringsAsFactors = FALSE
)

names(
  annotation_counts
) <- c(
  "annotation_refined",
  "n_cells"
)

annotation_counts <- annotation_counts[
  order(
    -annotation_counts$n_cells
  ),
  ,
  drop = FALSE
]

write.csv(
  annotation_counts,
  file.path(
    TAB_OUT,
    paste0(
      "03_refined_annotation_cell_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


parent_transfer_audit <- lapply(
  names(
    parent_only_labels
  ),
  function(parent_now) {

    source_now <- c(
      "Kupffer_Macrophage" = "Mphi",
      "HSC_Mesenchymal" = "HSC",
      "Hepatocyte" = "Hepatocyte",
      "LSEC" = "LSEC",
      "Cholangiocyte" = "Cholangiocyte",
      "Monocyte" = "Monocyte"
    )[[parent_now]]

    parent_cells_now <- colnames(whole)[
      whole$wholecell_layer1_R8_v701 ==
        parent_now
    ]

    unresolved_label <-
      parent_only_labels[[
        parent_now
      ]]

    unresolved_n <- sum(
      as.character(
        whole$annotation_refined_v701
      ) ==
        unresolved_label,
      na.rm = TRUE
    )

    source_owned_n <- sum(
      whole$refined_source_v701 ==
        source_now,
      na.rm = TRUE
    )

    source_owned_cross_parent_n <- sum(
      whole$refined_source_v701 ==
        source_now &
      as.character(
        whole$wholecell_layer1_R8_v701
      ) !=
        parent_now,
      na.rm = TRUE
    )

    data.frame(
      parent_lineage = parent_now,
      refined_source = source_now,
      broad_parent_cells = length(
        parent_cells_now
      ),
      final_source_owned_cells =
        source_owned_n,
      source_owned_cross_parent =
        source_owned_cross_parent_n,
      broad_parent_remainder =
        unresolved_n,
      stringsAsFactors = FALSE
    )
  }
)

parent_transfer_audit <- do.call(
  rbind,
  parent_transfer_audit
)

write.csv(
  parent_transfer_audit,
  file.path(
    TAB_OUT,
    paste0(
      "04_parent_to_refined_transfer_summary_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

msg(
  "Parent-to-refined transfer summary:"
)

print(
  parent_transfer_audit
)


# Full metadata export: small compared with a full Seurat RDS.
metadata_export <- data.frame(
  cell = colnames(
    whole
  ),
  sample =
    as.character(
      whole$sample_R8_v701
    ),
  condition =
    as.character(
      whole$condition_R8_v701
    ),
  broad_annotation =
    as.character(
      whole$wholecell_layer1_R8_v701
    ),
  refined_lineage =
    whole$refined_lineage_v701,
  refined_state_raw =
    whole$refined_state_raw_v701,
  refined_source =
    whole$refined_source_v701,
  refined_ownership_rule =
    whole$refined_ownership_rule_v701,
  refined_annotation =
    whole$annotation_refined_v701,
  stringsAsFactors = FALSE
)

write.csv(
  metadata_export,
  file.path(
    TAB_OUT,
    paste0(
      "05_wholecell_refined_metadata_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  metadata_export,
  file.path(
    OBJ_OUT,
    paste0(
      "Mouse_MASH_wholecell_refined_metadata_",
      VERSION,
      ".rds"
    )
  )
)


# Sample x refined annotation counts.
sample_counts <- as.data.frame(
  table(
    sample =
      whole$sample_R8_v701,
    refined_annotation =
      whole$annotation_refined_v701
  ),
  stringsAsFactors = FALSE
)

sample_counts <- sample_counts[
  sample_counts$Freq > 0,
  ,
  drop = FALSE
]

names(
  sample_counts
) <- c(
  "sample",
  "refined_annotation",
  "n_cells"
)

write.csv(
  sample_counts,
  file.path(
    TAB_OUT,
    paste0(
      "06_sample_by_refined_annotation_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# Condition x refined annotation counts.
condition_counts <- as.data.frame(
  table(
    condition =
      whole$condition_R8_v701,
    refined_annotation =
      whole$annotation_refined_v701
  ),
  stringsAsFactors = FALSE
)

condition_counts <- condition_counts[
  condition_counts$Freq > 0,
  ,
  drop = FALSE
]

names(
  condition_counts
) <- c(
  "condition",
  "refined_annotation",
  "n_cells"
)

write.csv(
  condition_counts,
  file.path(
    TAB_OUT,
    paste0(
      "07_condition_by_refined_annotation_counts_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# Detailed transfer maps.
for (
  source_name in names(
    transfer_maps
  )
) {

  write.csv(
    transfer_maps[[
      source_name
    ]],
    file.path(
      TAB_OUT,
      paste0(
        "transfer_map_",
        source_name,
        "_",
        VERSION,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 11. Final palette audit
# ==============================================================================

refined_present <- sort(
  unique(
    whole$annotation_refined_v701
  )
)

missing_refined_colors <- setdiff(
  refined_present,
  names(
    R8_REFINED
  )
)

if (length(missing_refined_colors)) {

  fallback <- grDevices::hcl.colors(
    length(
      missing_refined_colors
    ),
    palette = "Dynamic"
  )

  names(
    fallback
  ) <- missing_refined_colors

  R8_REFINED <- c(
    R8_REFINED,
    fallback
  )

  warning(
    "Fallback colors assigned to: ",
    paste(
      missing_refined_colors,
      collapse = ", "
    )
  )
}


palette_manifest <- data.frame(
  annotation = names(
    R8_REFINED
  ),
  hex = unname(
    R8_REFINED
  ),
  present_in_object =
    names(
      R8_REFINED
    ) %in%
      refined_present,
  stringsAsFactors = FALSE
)

write.csv(
  palette_manifest,
  file.path(
    TAB_OUT,
    paste0(
      "08_R8_refined_palette_manifest_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


broad_palette_manifest <- data.frame(
  annotation = names(
    R8_BROAD
  ),
  hex = unname(
    R8_BROAD
  ),
  present_in_object =
    names(
      R8_BROAD
    ) %in%
      broad_present,
  stringsAsFactors = FALSE
)

write.csv(
  broad_palette_manifest,
  file.path(
    TAB_OUT,
    paste0(
      "09_R8_broad_palette_manifest_",
      VERSION,
      ".csv"
    )
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Factor levels for stable legends
# ==============================================================================

broad_order <- c(
  "Hepatocyte",
  "Cholangiocyte",
  "LSEC",
  "Vascular_endothelial",
  "HSC_Mesenchymal",
  "Kupffer_Macrophage",
  "Monocyte",
  "Neutrophil",
  "Dendritic_cell",
  "NK_cell",
  "T_cell",
  "B_cell",
  "Plasma_cell",
  "Mesothelial",
  "Cycling",
  "RBC",
  "Platelet",
  "Other"
)

broad_order <- c(
  broad_order[
    broad_order %in%
      broad_present
  ],
  setdiff(
    broad_present,
    broad_order
  )
)

whole$wholecell_layer1_R8_v701 <-
  factor(
    whole$wholecell_layer1_R8_v701,
    levels = broad_order
  )


refined_order <- c(
  unname(
    HEP_DISPLAY[
      HEP_STATES
    ]
  ),
  "Hep | parent-only/unresolved",

  unname(
    CHOL_DISPLAY[
      CHOL_STATES
    ]
  ),
  "Chol | parent-only/unresolved",

  unname(
    LSEC_DISPLAY[
      LSEC_STATES
    ]
  ),
  "LSEC | parent-only/unresolved",

  "Vascular_endothelial",

  unname(
    HSC_DISPLAY[
      HSC_STATES
    ]
  ),
  "HSC | parent-only/unresolved",

  unname(
    MPHI_DISPLAY[
      MPHI_STATES
    ]
  ),
  "Mphi | parent-only/unresolved",

  unname(
    MONO_DISPLAY[
      MONO_STATES
    ]
  ),
  "Mono | parent-only/unresolved",

  "Neutrophil",
  "Dendritic_cell",
  "NK_cell",
  "T_cell",
  "B_cell",
  "Plasma_cell",
  "Mesothelial",
  "Cycling",
  "RBC",
  "Platelet",
  "Other"
)

refined_order <- unique(
  c(
    refined_order[
      refined_order %in%
        refined_present
    ],
    setdiff(
      refined_present,
      refined_order
    )
  )
)

whole$annotation_refined_v701 <-
  factor(
    whole$annotation_refined_v701,
    levels = refined_order
  )


# ==============================================================================
# 13. Broad Layer1 R8 UMAP
# ==============================================================================

msg(
  "Drawing broad Layer1 R8 UMAP..."
)

p_broad <- DimPlot(
  whole,
  reduction = WHOLE_UMAP,
  group.by =
    "wholecell_layer1_R8_v701",
  cols = R8_BROAD,
  pt.size = PT_ALL,
  raster = FALSE,
  shuffle = TRUE,
  seed = 7000
) +
  ggtitle(
    "Mouse liver whole-cell UMAP | broad annotation | R8"
  ) +
  coord_equal()

p_broad <- add_legend_style(
  p_broad,
  ncol_legend =
    LEGEND_NCOL_BROAD,
  base_size = 10
)

save_pdf(
  p_broad,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_broad_R8_",
      VERSION,
      ".pdf"
    )
  ),
  PDF_WIDTH_ALL,
  PDF_HEIGHT_ALL
)

save_png(
  p_broad,
  file.path(
    FIG_OUT,
    paste0(
      "01_wholecell_broad_R8_",
      VERSION,
      ".png"
    )
  ),
  PDF_WIDTH_ALL,
  PDF_HEIGHT_ALL
)


# ==============================================================================
# 14. Broad UMAP split by biological sample
# ==============================================================================

msg(
  "Drawing broad Layer1 R8 UMAP by sample..."
)

p_broad_sample <- DimPlot(
  whole,
  reduction = WHOLE_UMAP,
  group.by =
    "wholecell_layer1_R8_v701",
  split.by =
    "sample_R8_v701",
  cols = R8_BROAD,
  pt.size = PT_SPLIT,
  raster = FALSE,
  ncol = 3,
  shuffle = TRUE,
  seed = 7001
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | broad annotation | biological samples | R8"
  )

p_broad_sample <-
  p_broad_sample &
  theme_classic(
    base_size = 8
  ) &
  theme(
    legend.title =
      element_blank(),
    legend.text =
      element_text(
        size = 7
      )
  ) &
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 2.0,
        alpha = 1
      ),
      ncol =
        LEGEND_NCOL_BROAD
    )
  )

save_pdf(
  p_broad_sample,
  file.path(
    FIG_OUT,
    paste0(
      "02_wholecell_broad_by_sample_R8_",
      VERSION,
      ".pdf"
    )
  ),
  PDF_WIDTH_SPLIT,
  PDF_HEIGHT_SPLIT
)


# ==============================================================================
# 15. Refined all-cell R8 UMAP
# ==============================================================================

msg(
  "Drawing refined all-cell R8 UMAP..."
)

p_refined <- DimPlot(
  whole,
  reduction = WHOLE_UMAP,
  group.by =
    "annotation_refined_v701",
  cols = R8_REFINED,
  pt.size = PT_ALL,
  raster = FALSE,
  shuffle = TRUE,
  seed = 7002
) +
  ggtitle(
    "Mouse liver whole-cell UMAP | current refined annotations | R8"
  ) +
  coord_equal()

p_refined <- add_legend_style(
  p_refined,
  ncol_legend =
    LEGEND_NCOL_REFINED,
  base_size = 9
)

save_pdf(
  p_refined,
  file.path(
    FIG_OUT,
    paste0(
      "03_wholecell_refined_allstates_R8_",
      VERSION,
      ".pdf"
    )
  ),
  18,
  12
)

save_png(
  p_refined,
  file.path(
    FIG_OUT,
    paste0(
      "03_wholecell_refined_allstates_R8_",
      VERSION,
      ".png"
    )
  ),
  18,
  12
)


# ==============================================================================
# 16. Refined UMAP split by sample
# ==============================================================================

msg(
  "Drawing refined R8 UMAP by sample..."
)

p_refined_sample <- DimPlot(
  whole,
  reduction = WHOLE_UMAP,
  group.by =
    "annotation_refined_v701",
  split.by =
    "sample_R8_v701",
  cols = R8_REFINED,
  pt.size = PT_SPLIT,
  raster = FALSE,
  ncol = 3,
  shuffle = TRUE,
  seed = 7003
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | current refined annotations | biological samples | R8"
  )

p_refined_sample <-
  p_refined_sample &
  theme_classic(
    base_size = 7.5
  ) &
  theme(
    legend.title =
      element_blank(),
    legend.text =
      element_text(
        size = 6.2
      )
  ) &
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 1.8,
        alpha = 1
      ),
      ncol =
        LEGEND_NCOL_REFINED
    )
  )

save_pdf(
  p_refined_sample,
  file.path(
    FIG_OUT,
    paste0(
      "04_wholecell_refined_by_sample_R8_",
      VERSION,
      ".pdf"
    )
  ),
  22,
  13
)


# ==============================================================================
# 17. Refined UMAP split by condition
# ==============================================================================

msg(
  "Drawing refined R8 UMAP by condition..."
)

p_refined_condition <- DimPlot(
  whole,
  reduction = WHOLE_UMAP,
  group.by =
    "annotation_refined_v701",
  split.by =
    "condition_R8_v701",
  cols = R8_REFINED,
  pt.size = PT_SPLIT,
  raster = FALSE,
  ncol = 2,
  shuffle = TRUE,
  seed = 7004
) +
  plot_annotation(
    title =
      "Mouse liver whole-cell UMAP | current refined annotations | condition | R8"
  )

p_refined_condition <-
  p_refined_condition &
  theme_classic(
    base_size = 8
  ) &
  theme(
    legend.title =
      element_blank(),
    legend.text =
      element_text(
        size = 6.5
      )
  ) &
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 1.9,
        alpha = 1
      ),
      ncol =
        LEGEND_NCOL_REFINED
    )
  )

save_pdf(
  p_refined_condition,
  file.path(
    FIG_OUT,
    paste0(
      "05_wholecell_refined_by_condition_R8_",
      VERSION,
      ".pdf"
    )
  ),
  20,
  12
)


# ==============================================================================
# 18. Lineage-focused refined highlight UMAPs
# ==============================================================================

msg(
  "Drawing lineage-focused refined highlight UMAPs..."
)

highlight_specs <- list(

  Mphi = list(
    lineage = "Mphi",
    palette = R8_MPHI,
    title =
      "Whole liver | current Kupffer/Macrophage states"
  ),

  HSC = list(
    lineage = "HSC",
    palette = R8_HSC,
    title =
      "Whole liver | current HSC states"
  ),

  Hepatocyte = list(
    lineage = "Hepatocyte",
    palette = R8_HEP,
    title =
      "Whole liver | current Hepatocyte states"
  ),

  LSEC = list(
    lineage = "LSEC",
    palette = R8_LSEC,
    title =
      "Whole liver | current LSEC states"
  ),

  Cholangiocyte = list(
    lineage = "Cholangiocyte",
    palette = R8_CHOL,
    title =
      "Whole liver | current Cholangiocyte states"
  ),

  Monocyte = list(
    lineage = "Monocyte",
    palette = R8_MONO,
    title =
      "Whole liver | current Monocyte states"
  )
)


um <- Embeddings(
  whole,
  WHOLE_UMAP
)

for (
  nm in names(
    highlight_specs
  )
) {

  hs <- highlight_specs[[
    nm
  ]]

  refined_now <- as.character(
    whole$annotation_refined_v701
  )

  lineage_now <- whole$refined_lineage_v701

  # Include parent-only/unresolved cells for that lineage.
  parent_prefix <- switch(
    nm,
    Mphi = "Mphi | parent-only/unresolved",
    HSC = "HSC | parent-only/unresolved",
    Hepatocyte = "Hep | parent-only/unresolved",
    LSEC = "LSEC | parent-only/unresolved",
    Cholangiocyte = "Chol | parent-only/unresolved",
    Monocyte = "Mono | parent-only/unresolved"
  )

  highlight <- rep(
    "Other cells",
    ncol(whole)
  )

  idx_refined <- lineage_now ==
    hs$lineage

  idx_refined[
    is.na(
      idx_refined
    )
  ] <- FALSE

  idx_parent_only <-
    refined_now ==
      parent_prefix

  highlight[
    idx_refined
  ] <- refined_now[
    idx_refined
  ]

  highlight[
    idx_parent_only
  ] <- parent_prefix

  highlight_levels <- c(
    "Other cells",
    names(
      hs$palette
    )[
      names(
        hs$palette
      ) %in%
        unique(
          highlight
        )
    ]
  )

  highlight <- factor(
    highlight,
    levels =
      highlight_levels
  )

  plot_df <- data.frame(
    UMAP_1 = um[
      ,
      1
    ],
    UMAP_2 = um[
      ,
      2
    ],
    highlight = highlight,
    stringsAsFactors = FALSE
  )

  # Plot grey background first, then lineage states.
  plot_df$draw_order <-
    plot_df$highlight !=
      "Other cells"

  plot_df <- plot_df[
    order(
      plot_df$draw_order
    ),
    ,
    drop = FALSE
  ]

  palette_now <- c(
    "Other cells" =
      "#E5E5E5",
    hs$palette
  )

  p_hi <- ggplot(
    plot_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = highlight
    )
  ) +
    geom_point(
      size = PT_HIGHLIGHT,
      alpha = 0.95
    ) +
    scale_color_manual(
      values = palette_now,
      drop = TRUE
    ) +
    coord_equal() +
    labs(
      title =
        hs$title,
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
      ),
      axis.text =
        element_blank(),
      axis.ticks =
        element_blank(),
      legend.title =
        element_blank(),
      legend.text =
        element_text(
          size = 7.5
        )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 2.2,
          alpha = 1
        ),
        ncol = 2
      )
    )

  save_pdf(
    p_hi,
    file.path(
      FIG_OUT,
      paste0(
        "06_highlight_",
        nm,
        "_R8_",
        VERSION,
        ".pdf"
      )
    ),
    12,
    9
  )
}


# ==============================================================================
# 19. Optional full updated Seurat RDS
# ==============================================================================

if (SAVE_UPDATED_WHOLE_RDS) {

  updated_rds <- file.path(
    OBJ_OUT,
    paste0(
      "Mouse_MASH_wholecell_current_refined_annotations_",
      VERSION,
      ".rds"
    )
  )

  msg(
    "Saving full updated Seurat RDS: ",
    updated_rds
  )

  saveRDS(
    whole,
    updated_rds,
    compress = FALSE
  )
}


# ==============================================================================
# 20. Final README / session info
# ==============================================================================

readme <- c(

  paste0(
    "Mouse MASH whole-cell refined annotation UMAP ",
    VERSION
  ),

  "============================================================",
  "",

  paste0(
    "Whole-cell source: ",
    WHOLE_RDS
  ),

  paste0(
    "Whole-cell Layer1 metadata: ",
    layer1_col
  ),

  paste0(
    "Whole-cell UMAP reduction: ",
    WHOLE_UMAP
  ),

  paste0(
    "Whole cells: ",
    ncol(
      whole
    )
  ),

  "",

  "Design:",
  "- No reintegration",
  "- No reclustering",
  "- No new UMAP",
  "- Exact barcode transfer only",
  "- Dedicated refined lineage objects may override older broad-parent labels",
  "- Monocyte v6.9.4 has priority over historical Clean-B Mphi for overlapping cells",
  "- Mphi subtype transfer is restricted to broad Kupffer_Macrophage cells not present in Monocyte v6.9.4",
  "- Parent-only/unresolved cells remain explicitly visible in lineage-family colors",
  "- Source RDS is never overwritten",

  "",

  "Detailed frozen sources:",
  paste0(
    "- Mphi: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "Mphi"
    ]
  ),
  paste0(
    "- HSC: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "HSC"
    ]
  ),
  paste0(
    "- Hepatocyte: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "Hepatocyte"
    ]
  ),
  paste0(
    "- LSEC: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "LSEC"
    ]
  ),
  paste0(
    "- Cholangiocyte: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "Cholangiocyte"
    ]
  ),
  paste0(
    "- Monocyte: ",
    source_audit_df$RDS[
      source_audit_df$source ==
        "Monocyte"
    ]
  ),

  "",

  "Primary display outputs:",
  paste0(
    "- ",
    file.path(
      FIG_OUT,
      paste0(
        "01_wholecell_broad_R8_",
        VERSION,
        ".pdf"
      )
    )
  ),
  paste0(
    "- ",
    file.path(
      FIG_OUT,
      paste0(
        "03_wholecell_refined_allstates_R8_",
        VERSION,
        ".pdf"
      )
    )
  ),
  paste0(
    "- ",
    file.path(
      FIG_OUT,
      paste0(
        "04_wholecell_refined_by_sample_R8_",
        VERSION,
        ".pdf"
      )
    )
  ),
  paste0(
    "- ",
    file.path(
      FIG_OUT,
      paste0(
        "05_wholecell_refined_by_condition_R8_",
        VERSION,
        ".pdf"
      )
    )
  ),

  "",

  paste0(
    "SAVE_UPDATED_WHOLE_RDS: ",
    SAVE_UPDATED_WHOLE_RDS
  )
)

writeLines(
  readme,
  file.path(
    OUT,
    paste0(
      "README_",
      VERSION,
      ".txt"
    )
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_OUT,
    paste0(
      "sessionInfo_",
      VERSION,
      ".txt"
    )
  )
)


cat("\n====================================================\n")
cat("WHOLE-CELL R8 UMAP UPDATE COMPLETE\n")
cat("Version:", VERSION, "\n")
cat("Whole cells:", ncol(whole), "\n")
cat("Layer1 source:", layer1_col, "\n")
cat("UMAP:", WHOLE_UMAP, "\n")
cat("\n=== ANNOTATION SOURCE AUDIT ===\n")
print(source_audit_df)
cat("\n=== PARENT TO REFINED TRANSFER ===\n")
print(parent_transfer_audit)
cat("\n=== TOP REFINED ANNOTATIONS ===\n")
print(head(annotation_counts, 30))
cat("\nOutput:", OUT, "\n")
cat("====================================================\n")
