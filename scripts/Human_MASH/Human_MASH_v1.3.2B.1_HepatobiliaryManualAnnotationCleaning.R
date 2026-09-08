# ============================================================
# Human MASH v1.3.2B.1
# Hepatobiliary Manual Annotation & Cleaning
#
# PURPOSE
#   Apply the manually reviewed cluster map to the v1.3.2B
#   hepatobiliary RPCA object, remove clear non-hepatobiliary
#   contaminants, preserve biologically useful hepatocyte /
#   cholangiocyte / ductular states, and create clean objects
#   for downstream Human MASH progression analyses.
#
# IMPORTANT
#   - NO RPCA reintegration.
#   - NO reclustering.
#   - NO FindAllMarkers rerun.
#   - Reads the v1.3.2B pre-plot checkpoint.
#   - Requires the v1.3.2B-QC1 decision-support result and
#     validates that clusters 0 and 7 are closer to the
#     Hepatocyte reference than the LSEC reference.
#   - Canonical whole-liver v1.3 RDS is NEVER modified here.
#   - The 44 strict outside-core ductal-like candidates from
#     v1.3.2B remain audit-only and are NOT automatically rescued.
#
# MANUAL CLUSTER MAP (res0.6)
#
#   0   Mature_hepatocyte
#   1   Mature_hepatocyte_periportal_like
#   2   Hepatocyte_acute_phase_stress
#   3   Mature_hepatocyte_periportal_like
#   4   Mature_hepatocyte_pericentral
#   5   Hepatocyte_acute_phase_stress
#   6   Cholangiocyte_ductal
#   7   Hepatocyte_low_maturity_transitional
#   8   Injury_transitional_epithelial
#   9   EPCAM_positive_ductular_progenitor_like
#   10  Hepatocyte_acute_phase_stress
#   11  Exclude_immune_T_cell
#   12  Mature_hepatocyte_periportal_like
#   13  Cholangiocyte_mucinous_biliary
#   14  Injury_transitional_hepatocyte
#   15  Exclude_mesenchymal
#   16  Exclude_macrophage
#
# SAMPLE-BIAS FLAGS FROM v1.3.2B
#   cluster 9  : LT2-biased (~75%)
#   cluster 12 : LT4-biased (~91%)
#   cluster 14 : LT1-biased (~98%)
#
# INTERPRETATION
#   - Sample-biased clusters are retained as biologically
#     interpretable states but must not be treated as robust
#     cross-sample progression states without additional support.
#   - "Hepatoblast_or_immature_hepatocyte" is not retained as a
#     final manual label because fetal hepatoblast evidence was
#     not convincing (notably DLK1 absent in v1.3.2B).
#
# OUTPUT STATUS
#   MANUAL_CLEANED_REVIEW_REQUIRED
#   Clean hepatobiliary object is created.
#   Canonical whole-liver v1.3 object is not modified.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

options(stringsAsFactors = FALSE)

seed <- 20260908
set.seed(seed)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

parent_outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B_hepatobiliary_refinement"
)

checkpoint_rds <- file.path(
  parent_outdir,
  "objects",
  "Human_MASH_Hepatobiliary_v1.3.2B_preplot_checkpoint.rds"
)

top20_csv <- file.path(
  parent_outdir,
  "15_hepatobiliary_cluster_markers_top20.csv"
)

outside_core_ductal_csv <- file.path(
  parent_outdir,
  "02_strict_ductal_rescue_candidate_cells_OUTSIDE_CORE.csv"
)

qc1_outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B_QC1_LSEC_vs_Hepatocyte_audit"
)

qc1_decision_csv <- file.path(
  qc1_outdir,
  "11_target_cluster_QC_decision_support.csv"
)

qc1_flags_txt <- file.path(
  qc1_outdir,
  "18_v1.3.2B_QC1_review_flags.txt"
)

outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.2B.1_hepatobiliary_manual_annotation_cleaning"
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

cluster_col <- "hepatobiliary_cluster_v132B"
reduction_name <- "umapHepatoRPCA_v132B"

expected_cells <- 47425L
expected_clusters <- as.character(0:16)

exclude_clusters <- c(
  "11",
  "15",
  "16"
)

sample_bias_threshold <- 70

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.2B.1\n")
cat("Hepatobiliary Manual Annotation & Cleaning\n")
cat("============================================================\n\n")

required_files <- c(
  checkpoint_rds,
  qc1_decision_csv,
  qc1_flags_txt
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

hepato_all <- readRDS(
  checkpoint_rds
)

DefaultAssay(
  hepato_all
) <- "RNA"

required_meta <- c(
  celltype_col,
  sample_col,
  nas_col,
  nasgrp_col,
  cluster_col,
  "current_hepatobiliary_label_v132B"
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    hepato_all@meta.data
  )
)

if (
  length(
    missing_meta
  ) > 0
) {

  stop(
    "Missing checkpoint metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  ncol(
    hepato_all
  ) != expected_cells
) {

  stop(
    "Checkpoint cell-count mismatch. Expected ",
    expected_cells,
    " but found ",
    ncol(
      hepato_all
    ),
    ". Manual cluster map must NOT be applied."
  )
}

cluster_vec <- as.character(
  hepato_all[[cluster_col]][, 1]
)

observed_clusters <- unique(
  cluster_vec
)

if (
  !setequal(
    observed_clusters,
    expected_clusters
  )
) {

  stop(
    "Cluster-set mismatch.\nExpected: ",
    paste(
      expected_clusters,
      collapse = ","
    ),
    "\nObserved: ",
    paste(
      sort(
        observed_clusters
      ),
      collapse = ","
    ),
    "\nManual cluster map must NOT be applied."
  )
}

if (
  !reduction_name %in%
    Reductions(
      hepato_all
    )
) {

  stop(
    "Required reduction missing: ",
    reduction_name
  )
}

cat(
  "Checkpoint cells :",
  ncol(
    hepato_all
  ),
  "\n"
)

cat(
  "Clusters         :",
  length(
    observed_clusters
  ),
  "\n\n"
)

# ------------------------------------------------------------
# 3. Validate QC1 decision for clusters 0 and 7
# ------------------------------------------------------------

qc1 <- read.csv(
  qc1_decision_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_qc1_cols <- c(
  "group",
  "QC_axis_closer_to"
)

missing_qc1_cols <- setdiff(
  required_qc1_cols,
  colnames(
    qc1
  )
)

if (
  length(
    missing_qc1_cols
  ) > 0
) {

  stop(
    "QC1 decision-support CSV lacks required columns: ",
    paste(
      missing_qc1_cols,
      collapse = ", "
    )
  )
}

qc1_required_groups <- c(
  "Cluster0_target",
  "Cluster7_target"
)

qc1_target <- qc1[
  qc1$group %in%
    qc1_required_groups,
  ,
  drop = FALSE
]

if (
  !setequal(
    qc1_target$group,
    qc1_required_groups
  )
) {

  stop(
    "QC1 decision-support file does not contain both Cluster0_target and Cluster7_target."
  )
}

if (
  any(
    qc1_target$QC_axis_closer_to !=
      "Hepatocyte_reference"
  )
) {

  stop(
    "QC1 no longer supports Hepatocyte-reference assignment for clusters 0 and/or 7. ",
    "Do NOT apply this manual map."
  )
}

write.csv(
  qc1_target,
  file.path(
    outdir,
    "01_QC1_cluster0_7_validation.csv"
  ),
  row.names = FALSE
)

cat("=== QC1 VALIDATION ===\n")
print(
  qc1_target[
    ,
    intersect(
      c(
        "group",
        "n_cells",
        "median_Hepatocyte_score",
        "median_LSEC_score",
        "median_LSEC_minus_Hepatocyte",
        "pct_cells_LSEC_gt_Hepatocyte",
        "QC_axis_closer_to"
      ),
      colnames(
        qc1_target
      )
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nPASS: clusters 0 and 7 are QC1-supported as Hepatocyte-reference closer.\n\n")

# ------------------------------------------------------------
# 4. Manual cluster map
# ------------------------------------------------------------

manual_map <- data.frame(

  cluster = as.character(
    0:16
  ),

  manual_state = c(
    "Mature_hepatocyte",
    "Mature_hepatocyte_periportal_like",
    "Hepatocyte_acute_phase_stress",
    "Mature_hepatocyte_periportal_like",
    "Mature_hepatocyte_pericentral",
    "Hepatocyte_acute_phase_stress",
    "Cholangiocyte_ductal",
    "Hepatocyte_low_maturity_transitional",
    "Injury_transitional_epithelial",
    "EPCAM_positive_ductular_progenitor_like",
    "Hepatocyte_acute_phase_stress",
    "Exclude_immune_T_cell",
    "Mature_hepatocyte_periportal_like",
    "Cholangiocyte_mucinous_biliary",
    "Injury_transitional_hepatocyte",
    "Exclude_mesenchymal",
    "Exclude_macrophage"
  ),

  broad_celltype = c(
    "Hepatocyte",
    "Hepatocyte",
    "Hepatocyte",
    "Hepatocyte",
    "Hepatocyte",
    "Hepatocyte",
    "Cholangiocyte",
    "Hepatocyte",
    "Hepatocyte",
    "Ductular_progenitor",
    "Hepatocyte",
    "Excluded_nonhepatobiliary",
    "Hepatocyte",
    "Cholangiocyte",
    "Hepatocyte",
    "Excluded_nonhepatobiliary",
    "Excluded_nonhepatobiliary"
  ),

  keep_clean = c(
    TRUE,   # 0
    TRUE,   # 1
    TRUE,   # 2
    TRUE,   # 3
    TRUE,   # 4
    TRUE,   # 5
    TRUE,   # 6
    TRUE,   # 7
    TRUE,   # 8
    TRUE,   # 9
    TRUE,   # 10
    FALSE,  # 11
    TRUE,   # 12
    TRUE,   # 13
    TRUE,   # 14
    FALSE,  # 15
    FALSE   # 16
  ),

  manual_confidence = c(
    "high",      # 0
    "moderate",  # 1
    "moderate",  # 2
    "moderate",  # 3
    "high",      # 4
    "moderate",  # 5
    "high",      # 6
    "moderate",  # 7
    "moderate",  # 8
    "moderate",  # 9
    "moderate",  # 10
    "high",      # 11
    "moderate",  # 12
    "high",      # 13
    "moderate",  # 14
    "high",      # 15
    "high"       # 16
  ),

  evidence_basis = c(
    "QC1 Hepatocyte-reference; ALB/APOA1/TTR/ASGR1/HNF4A/CYP3A4",
    "hepatocyte markers with periportal-biased program",
    "hepatocyte identity with acute-phase/stress program",
    "hepatocyte markers with periportal-biased program",
    "hepatocyte identity with pericentral/GLUL/CYP2E1 program",
    "hepatocyte identity with acute-phase/stress program",
    "EPCAM/PROM1/KRT19/KRT7/SOX9/CFTR coherent biliary program",
    "QC1 excludes LSEC identity; weaker mature-hepatocyte program",
    "injury/transitional epithelial program; no fetal-hepatoblast support",
    "EPCAM/PROM1/KRT19/KRT7/SOX9 ductular/progenitor-like program",
    "hepatocyte identity with acute-phase/stress program",
    "CD2/CD8A/TRBC/ITK T-cell markers",
    "periportal-like hepatocyte program; strong single-sample bias flag",
    "MUC1/MUC5B/TFF2/CEACAM6/AGR2 biliary-mucinous program",
    "injury/transitional hepatocyte program; strong single-sample bias flag",
    "MYH11/CACNA1C/LAMA2/HGF perivascular-mesenchymal markers",
    "C1QA/C1QB/C1QC/MARCO/FOLR2/TIMD4 macrophage markers"
  ),

  predefined_review_note = c(
    "QC1-resolved: retain as Hepatocyte",
    NA,
    NA,
    NA,
    NA,
    NA,
    "bona fide cholangiocyte/ductal",
    "QC1-resolved: not LSEC; retain as transitional Hepatocyte",
    "retain as injury/transitional epithelial state",
    "LT2-biased in v1.3.2B; interpret cautiously",
    NA,
    "exclude non-hepatobiliary contaminant",
    "LT4-biased in v1.3.2B; interpret cautiously",
    "bona fide cholangiocyte/biliary state",
    "LT1-biased in v1.3.2B; interpret cautiously",
    "exclude non-hepatobiliary contaminant",
    "exclude non-hepatobiliary contaminant"
  ),

  stringsAsFactors = FALSE
)

if (
  !setequal(
    manual_map$cluster,
    expected_clusters
  )
) {

  stop(
    "Internal manual-map error."
  )
}

write.csv(
  manual_map,
  file.path(
    outdir,
    "02_manual_cluster_annotation_map.csv"
  ),
  row.names = FALSE
)

cat("=== MANUAL CLUSTER MAP ===\n")
print(
  manual_map[
    ,
    c(
      "cluster",
      "manual_state",
      "broad_celltype",
      "keep_clean",
      "manual_confidence"
    )
  ],
  row.names = FALSE
)
cat("\n")

# ------------------------------------------------------------
# 5. Apply manual annotation
# ------------------------------------------------------------

state_map <- setNames(
  manual_map$manual_state,
  manual_map$cluster
)

broad_map <- setNames(
  manual_map$broad_celltype,
  manual_map$cluster
)

keep_map <- setNames(
  manual_map$keep_clean,
  manual_map$cluster
)

confidence_map <- setNames(
  manual_map$manual_confidence,
  manual_map$cluster
)

hepato_all$hepatobiliary_state_manual_v132B1 <-
  unname(
    state_map[
      cluster_vec
    ]
  )

hepato_all$hepatobiliary_broad_manual_v132B1 <-
  unname(
    broad_map[
      cluster_vec
    ]
  )

hepato_all$hepatobiliary_keep_clean_v132B1 <-
  as.logical(
    unname(
      keep_map[
        cluster_vec
      ]
    )
  )

hepato_all$hepatobiliary_manual_confidence_v132B1 <-
  unname(
    confidence_map[
      cluster_vec
    ]
  )

hepato_all$hepatobiliary_status_v132B1 <- ifelse(
  hepato_all$hepatobiliary_keep_clean_v132B1,
  "clean_hepatobiliary",
  "excluded_nonhepatobiliary"
)

if (
  anyNA(
    hepato_all$hepatobiliary_state_manual_v132B1
  ) ||
    anyNA(
      hepato_all$hepatobiliary_broad_manual_v132B1
    ) ||
    anyNA(
      hepato_all$hepatobiliary_keep_clean_v132B1
    )
) {

  stop(
    "NA introduced while applying manual hepatobiliary map."
  )
}

# ------------------------------------------------------------
# 6. Cluster counts / sample dominance
# ------------------------------------------------------------

sample_vec_all <- as.character(
  hepato_all[[sample_col]][, 1]
)

cluster_sample <- as.data.frame(
  table(
    cluster =
      cluster_vec,
    sample =
      sample_vec_all
  ),
  stringsAsFactors = FALSE
)

cluster_sample$cluster_total <- ave(
  cluster_sample$Freq,
  cluster_sample$cluster,
  FUN = sum
)

cluster_sample$within_cluster_pct <- ifelse(
  cluster_sample$cluster_total > 0,
  100 *
    cluster_sample$Freq /
    cluster_sample$cluster_total,
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
  ),
  drop = FALSE
]

names(
  cluster_dominance
)[3] <- "max_sample_contribution_pct"

cluster_dominance <- cluster_dominance[
  order(
    as.integer(
      cluster_dominance$cluster
    )
  ),
  ,
  drop = FALSE
]

manual_map$max_sample <- cluster_dominance$sample[
  match(
    manual_map$cluster,
    cluster_dominance$cluster
  )
]

manual_map$max_sample_contribution_pct <-
  cluster_dominance$max_sample_contribution_pct[
    match(
      manual_map$cluster,
      cluster_dominance$cluster
    )
  ]

manual_map$sample_bias_ge70 <-
  manual_map$max_sample_contribution_pct >=
    sample_bias_threshold

write.csv(
  manual_map,
  file.path(
    outdir,
    "03_manual_cluster_annotation_map_with_sample_bias.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Create clean hepatobiliary object
# ------------------------------------------------------------

clean_cells <- colnames(
  hepato_all
)[
  hepato_all$hepatobiliary_keep_clean_v132B1
]

excluded_cells <- colnames(
  hepato_all
)[
  !hepato_all$hepatobiliary_keep_clean_v132B1
]

hepato_clean <- subset(
  hepato_all,
  cells = clean_cells
)

DefaultAssay(
  hepato_clean
) <- "RNA"

n_total <- ncol(
  hepato_all
)

n_clean <- ncol(
  hepato_clean
)

n_excluded <- n_total -
  n_clean

if (
  any(
    as.character(
      hepato_clean[[cluster_col]][, 1]
    ) %in%
      exclude_clusters
  )
) {

  stop(
    "Excluded clusters remain in clean hepatobiliary object."
  )
}

# ------------------------------------------------------------
# 8. Create lineage-specific clean objects
# ------------------------------------------------------------

hepatocyte_cells <- colnames(
  hepato_clean
)[
  hepato_clean$hepatobiliary_broad_manual_v132B1 ==
    "Hepatocyte"
]

cholangiocyte_ductular_cells <- colnames(
  hepato_clean
)[
  hepato_clean$hepatobiliary_broad_manual_v132B1 %in%
    c(
      "Cholangiocyte",
      "Ductular_progenitor"
    )
]

hepatocyte_clean <- subset(
  hepato_clean,
  cells = hepatocyte_cells
)

cholangiocyte_ductular_clean <- subset(
  hepato_clean,
  cells = cholangiocyte_ductular_cells
)

DefaultAssay(
  hepatocyte_clean
) <- "RNA"

DefaultAssay(
  cholangiocyte_ductular_clean
) <- "RNA"

# ------------------------------------------------------------
# 9. Cleaning summary
# ------------------------------------------------------------

cleaning_summary <- data.frame(
  metric = c(
    "input_hepatobiliary_candidate_cells",
    "clean_hepatobiliary_cells",
    "excluded_nonhepatobiliary_cells",
    "clean_percent",
    "excluded_percent",
    "clean_hepatocyte_cells",
    "clean_cholangiocyte_plus_ductular_cells"
  ),
  value = c(
    n_total,
    n_clean,
    n_excluded,
    100 *
      n_clean /
      n_total,
    100 *
      n_excluded /
      n_total,
    ncol(
      hepatocyte_clean
    ),
    ncol(
      cholangiocyte_ductular_clean
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cleaning_summary,
  file.path(
    outdir,
    "04_cleaning_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Counts by manual state / broad class
# ------------------------------------------------------------

state_counts <- as.data.frame(
  table(
    state =
      hepato_clean$hepatobiliary_state_manual_v132B1
  ),
  stringsAsFactors = FALSE
)

names(
  state_counts
)[2] <- "n_cells"

state_counts$percent <- 100 *
  state_counts$n_cells /
  sum(
    state_counts$n_cells
  )

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
    "05_clean_manual_state_counts.csv"
  ),
  row.names = FALSE
)

broad_counts <- as.data.frame(
  table(
    broad_celltype =
      hepato_clean$hepatobiliary_broad_manual_v132B1
  ),
  stringsAsFactors = FALSE
)

names(
  broad_counts
)[2] <- "n_cells"

broad_counts$percent <- 100 *
  broad_counts$n_cells /
  sum(
    broad_counts$n_cells
  )

broad_counts <- broad_counts[
  order(
    broad_counts$n_cells,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  broad_counts,
  file.path(
    outdir,
    "06_clean_broad_celltype_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Excluded-cluster audit
# ------------------------------------------------------------

excluded_audit <- as.data.frame(
  table(
    cluster =
      hepato_all[[cluster_col]][, 1],
    manual_state =
      hepato_all$hepatobiliary_state_manual_v132B1,
    sample =
      hepato_all[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

excluded_audit <- excluded_audit[
  excluded_audit$cluster %in%
    exclude_clusters &
    excluded_audit$Freq > 0,
  ,
  drop = FALSE
]

write.csv(
  excluded_audit,
  file.path(
    outdir,
    "07_excluded_nonhepatobiliary_by_cluster_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Current annotation vs manual crosswalk
# ------------------------------------------------------------

crosswalk <- as.data.frame(
  table(
    current_REF12 =
      hepato_all$current_hepatobiliary_label_v132B,
    manual_state =
      hepato_all$hepatobiliary_state_manual_v132B1,
    broad_celltype =
      hepato_all$hepatobiliary_broad_manual_v132B1
  ),
  stringsAsFactors = FALSE
)

crosswalk <- crosswalk[
  crosswalk$Freq > 0,
  ,
  drop = FALSE
]

write.csv(
  crosswalk,
  file.path(
    outdir,
    "08_current_REF12_vs_manual_crosswalk.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Clean state x sample tables
# ------------------------------------------------------------

sample_order <- sample_order_fun(
  hepato_clean[[sample_col]][, 1]
)

state_sample <- as.data.frame(
  table(
    state =
      hepato_clean$hepatobiliary_state_manual_v132B1,
    sample =
      hepato_clean[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

state_sample$state_total <- ave(
  state_sample$Freq,
  state_sample$state,
  FUN = sum
)

state_sample$sample_total <- ave(
  state_sample$Freq,
  state_sample$sample,
  FUN = sum
)

state_sample$within_state_pct <- ifelse(
  state_sample$state_total > 0,
  100 *
    state_sample$Freq /
    state_sample$state_total,
  NA_real_
)

state_sample$within_sample_pct <- ifelse(
  state_sample$sample_total > 0,
  100 *
    state_sample$Freq /
    state_sample$sample_total,
  NA_real_
)

state_sample$sample <- factor(
  state_sample$sample,
  levels =
    sample_order
)

write.csv(
  state_sample,
  file.path(
    outdir,
    "09_clean_manual_state_by_sample.csv"
  ),
  row.names = FALSE
)

broad_sample <- as.data.frame(
  table(
    broad_celltype =
      hepato_clean$hepatobiliary_broad_manual_v132B1,
    sample =
      hepato_clean[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

broad_sample$sample_total <- ave(
  broad_sample$Freq,
  broad_sample$sample,
  FUN = sum
)

broad_sample$within_sample_pct <- ifelse(
  broad_sample$sample_total > 0,
  100 *
    broad_sample$Freq /
    broad_sample$sample_total,
  NA_real_
)

broad_sample$sample <- factor(
  broad_sample$sample,
  levels =
    sample_order
)

write.csv(
  broad_sample,
  file.path(
    outdir,
    "10_clean_broad_celltype_by_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Outside-core ductal candidate status
# ------------------------------------------------------------

outside_core_count <- NA_integer_

if (
  file.exists(
    outside_core_ductal_csv
  )
) {

  outside_core_tbl <- read.csv(
    outside_core_ductal_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  outside_core_count <- nrow(
    outside_core_tbl
  )

  write.csv(
    outside_core_tbl,
    file.path(
      outdir,
      "11_outside_core_ductal_candidates_AUDIT_ONLY.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 15. Attach existing top20 markers to manual map
# ------------------------------------------------------------

if (
  file.exists(
    top20_csv
  )
) {

  top20 <- read.csv(
    top20_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (
    all(
      c(
        "cluster",
        "gene"
      ) %in%
        colnames(
          top20
        )
    )
  ) {

    top20$cluster <- as.character(
      top20$cluster
    )

    marker_summary <- aggregate(
      gene ~ cluster,
      data = top20,
      FUN = function(x) {

        paste(
          unique(
            as.character(x)
          ),
          collapse = ";"
        )
      }
    )

    names(
      marker_summary
    )[2] <-
      "observed_top20_markers"

    manual_map$observed_top20_markers <-
      marker_summary$observed_top20_markers[
        match(
          manual_map$cluster,
          marker_summary$cluster
        )
      ]

    write.csv(
      manual_map,
      file.path(
        outdir,
        "12_manual_cluster_map_with_top20_markers.csv"
      ),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------
# 16. Save pre-plot manual checkpoint
# ------------------------------------------------------------

manual_checkpoint_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatobiliary_v1.3.2B.1_preplot_manual_checkpoint.rds"
)

saveRDS(
  hepato_all,
  manual_checkpoint_rds
)

writeLines(
  c(
    "Human MASH v1.3.2B.1 pre-plot manual checkpoint",
    paste0(
      "input_cells=",
      n_total
    ),
    paste0(
      "clean_cells=",
      n_clean
    ),
    paste0(
      "excluded_cells=",
      n_excluded
    ),
    paste0(
      "manual_checkpoint_rds=",
      manual_checkpoint_rds
    ),
    "QC1_cluster0_7_validation=PASS",
    "manual_annotation_applied=TRUE",
    "cleaning_applied=TRUE",
    "canonical_modified=FALSE"
  ),
  file.path(
    outdir,
    "12A_preplot_manual_checkpoint.txt"
  )
)

cat(
  "PRE-PLOT MANUAL CHECKPOINT SAVED:\n",
  manual_checkpoint_rds,
  "\n\n"
)

# ------------------------------------------------------------
# 17. Palettes
# ------------------------------------------------------------

state_palette <- c(
  "Mature_hepatocyte" =
    "#00A6A6",
  "Mature_hepatocyte_periportal_like" =
    "#00C5C8",
  "Mature_hepatocyte_pericentral" =
    "#007F86",
  "Hepatocyte_acute_phase_stress" =
    "#FF8C00",
  "Hepatocyte_low_maturity_transitional" =
    "#3399FF",
  "Injury_transitional_epithelial" =
    "#FF5A36",
  "Injury_transitional_hepatocyte" =
    "#FF1A1A",
  "EPCAM_positive_ductular_progenitor_like" =
    "#76FF03",
  "Cholangiocyte_ductal" =
    "#00C853",
  "Cholangiocyte_mucinous_biliary" =
    "#00E676",
  "Exclude_immune_T_cell" =
    "#9E9E9E",
  "Exclude_mesenchymal" =
    "#757575",
  "Exclude_macrophage" =
    "#616161"
)

broad_palette <- c(
  "Hepatocyte" =
    "#00A6A6",
  "Cholangiocyte" =
    "#00C853",
  "Ductular_progenitor" =
    "#76FF03",
  "Excluded_nonhepatobiliary" =
    "#808080"
)

# ------------------------------------------------------------
# 18. Fig01 - manual state, all candidates
# ------------------------------------------------------------

p01 <- DimPlot(
  hepato_all,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_state_manual_v132B1",
  pt.size = 0.10,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    state_palette
) +
  ggtitle(
    "v1.3.2B.1 manual hepatobiliary annotation - all candidates"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig01_ManualAnnotation_AllCandidates_UMAP.pdf"
  ),
  p01,
  width = 13,
  height = 9
)

# ------------------------------------------------------------
# 19. Fig02 - clean manual state UMAP
# ------------------------------------------------------------

p02 <- DimPlot(
  hepato_clean,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_state_manual_v132B1",
  pt.size = 0.10,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    state_palette
) +
  ggtitle(
    "v1.3.2B.1 clean hepatobiliary states"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig02_CleanHepatobiliary_ManualState_UMAP.pdf"
  ),
  p02,
  width = 13,
  height = 9
)

# ------------------------------------------------------------
# 20. Fig03 - broad manual cell type
# ------------------------------------------------------------

p03 <- DimPlot(
  hepato_clean,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_broad_manual_v132B1",
  pt.size = 0.10,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  cols =
    broad_palette
) +
  ggtitle(
    "v1.3.2B.1 clean hepatobiliary broad cell types"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig03_CleanHepatobiliary_BroadCellType_UMAP.pdf"
  ),
  p03,
  width = 11,
  height = 8
)

# ------------------------------------------------------------
# 21. Fig04 - excluded contaminants
# ------------------------------------------------------------

excluded_obj <- subset(
  hepato_all,
  cells = excluded_cells
)

p04 <- DimPlot(
  excluded_obj,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_state_manual_v132B1",
  pt.size = 0.25,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    state_palette
) +
  ggtitle(
    "Excluded non-hepatobiliary contaminant clusters"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig04_ExcludedNonHepatobiliary_UMAP.pdf"
  ),
  p04,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 22. Fig05 - clean states split by LT
# ------------------------------------------------------------

p05 <- DimPlot(
  hepato_clean,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_state_manual_v132B1",
  split.by =
    sample_col,
  ncol = 4,
  pt.size = 0.06,
  raster = TRUE,
  cols =
    state_palette
)

ggsave(
  file.path(
    outdir,
    "Fig05_CleanHepatobiliary_ManualState_split_by_LT.pdf"
  ),
  p05,
  width = 20,
  height = 13
)

# ------------------------------------------------------------
# 23. Curated manual-validation marker panel
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

  Periportal = c(
    "CPS1",
    "ASS1",
    "ARG1",
    "HAL",
    "PCK1"
  ),

  Pericentral = c(
    "GLUL",
    "CYP2E1",
    "CYP1A2",
    "OAT"
  ),

  Acute_phase_stress = c(
    "SAA1",
    "SAA2",
    "CRP",
    "HP",
    "ORM1"
  ),

  Injury_transitional = c(
    "KRT8",
    "KRT18",
    "KRT23",
    "LGALS3",
    "CLDN1",
    "SPP1",
    "KLF6"
  ),

  Ductular_progenitor = c(
    "EPCAM",
    "PROM1",
    "KRT19",
    "KRT7",
    "SOX9",
    "TACSTD2"
  ),

  Cholangiocyte = c(
    "MUC1",
    "CFTR",
    "KRT17",
    "KRT19"
  ),

  Mucinous_biliary = c(
    "MUC5B",
    "TFF2",
    "CEACAM6",
    "AGR2"
  ),

  Fetal_hepatoblast_audit = c(
    "AFP",
    "DLK1",
    "IGF2",
    "H19"
  ),

  Contamination = c(
    "PTPRC",
    "LST1",
    "TYROBP",
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
        hepato_clean
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

# ------------------------------------------------------------
# 24. Fig06 - clean manual marker DotPlot
# ------------------------------------------------------------

p06 <- DotPlot(
  hepato_clean,
  features =
    marker_plot_list,
  assay = "RNA",
  group.by =
    "hepatobiliary_state_manual_v132B1",
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
      "v1.3.2B.1 clean hepatobiliary manual annotation validation",
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
    "Fig06_CleanHepatobiliary_ManualMarker_DotPlot.pdf"
  ),
  p06,
  width = 24,
  height = 10
)

write.csv(
  p06$data,
  file.path(
    outdir,
    "13_Fig06_ManualMarker_DotPlot_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 25. Fig07 - sample contribution within manual states
# ------------------------------------------------------------

p07 <- ggplot(
  state_sample,
  aes(
    x = sample,
    y = state,
    fill =
      within_state_pct
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
          within_state_pct
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
      "Sample contribution within clean manual hepatobiliary states",
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
    "Fig07_CleanManualState_SampleContribution.pdf"
  ),
  p07,
  width = 14,
  height = 9
)

# ------------------------------------------------------------
# 26. Fig08 - broad cell-type composition by LT
# ------------------------------------------------------------

p08 <- ggplot(
  broad_sample,
  aes(
    x = sample,
    y =
      within_sample_pct,
    fill =
      broad_celltype
  )
) +
  geom_col(
    width = 0.85
  ) +
  scale_fill_manual(
    values =
      broad_palette
  ) +
  labs(
    title =
      "Clean hepatobiliary broad cell-type composition by LT",
    x = NULL,
    y =
      "% of clean hepatobiliary cells within sample",
    fill =
      "Broad cell type"
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
    "Fig08_CleanBroadCellType_Composition_by_LT.pdf"
  ),
  p08,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 27. Fig09 - hepatocyte-state composition by LT
# ------------------------------------------------------------

hep_state_sample <- as.data.frame(
  table(
    state =
      hepatocyte_clean$hepatobiliary_state_manual_v132B1,
    sample =
      hepatocyte_clean[[sample_col]][, 1]
  ),
  stringsAsFactors = FALSE
)

hep_state_sample$sample_total <- ave(
  hep_state_sample$Freq,
  hep_state_sample$sample,
  FUN = sum
)

hep_state_sample$within_sample_pct <- ifelse(
  hep_state_sample$sample_total > 0,
  100 *
    hep_state_sample$Freq /
    hep_state_sample$sample_total,
  NA_real_
)

hep_state_sample$sample <- factor(
  hep_state_sample$sample,
  levels =
    sample_order
)

write.csv(
  hep_state_sample,
  file.path(
    outdir,
    "14_clean_hepatocyte_state_by_sample.csv"
  ),
  row.names = FALSE
)

p09 <- ggplot(
  hep_state_sample,
  aes(
    x = sample,
    y =
      within_sample_pct,
    fill = state
  )
) +
  geom_col(
    width = 0.85
  ) +
  scale_fill_manual(
    values =
      state_palette
  ) +
  labs(
    title =
      "Clean hepatocyte-state composition by LT",
    x = NULL,
    y =
      "% of clean hepatocytes within sample",
    fill =
      "Hepatocyte state"
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
    "Fig09_CleanHepatocyteState_Composition_by_LT.pdf"
  ),
  p09,
  width = 14,
  height = 9
)

# ------------------------------------------------------------
# 28. Fig10 - Cholangiocyte / ductular subset
# ------------------------------------------------------------

p10 <- DimPlot(
  cholangiocyte_ductular_clean,
  reduction =
    reduction_name,
  group.by =
    "hepatobiliary_state_manual_v132B1",
  pt.size = 0.20,
  raster = TRUE,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  cols =
    state_palette
) +
  ggtitle(
    "Clean cholangiocyte / ductular-progenitor compartment"
  ) +
  theme_classic(
    base_size = 11
  )

ggsave(
  file.path(
    outdir,
    "Fig10_Cholangiocyte_Ductular_ManualState_UMAP.pdf"
  ),
  p10,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 29. Per-cell metadata exports
# ------------------------------------------------------------

metadata_all <- data.frame(
  cell =
    colnames(
      hepato_all
    ),
  sample =
    hepato_all[[sample_col]][, 1],
  NAS_score =
    hepato_all[[nas_col]][, 1],
  NAS_group =
    hepato_all[[nasgrp_col]][, 1],
  cluster_v132B =
    hepato_all[[cluster_col]][, 1],
  current_REF12 =
    hepato_all$current_hepatobiliary_label_v132B,
  manual_state_v132B1 =
    hepato_all$hepatobiliary_state_manual_v132B1,
  broad_celltype_v132B1 =
    hepato_all$hepatobiliary_broad_manual_v132B1,
  manual_confidence_v132B1 =
    hepato_all$hepatobiliary_manual_confidence_v132B1,
  keep_clean_v132B1 =
    hepato_all$hepatobiliary_keep_clean_v132B1,
  status_v132B1 =
    hepato_all$hepatobiliary_status_v132B1,
  stringsAsFactors = FALSE
)

write.csv(
  metadata_all,
  file.path(
    outdir,
    "15_v1.3.2B.1_all_candidate_cell_metadata.csv"
  ),
  row.names = FALSE
)

metadata_clean <- metadata_all[
  metadata_all$keep_clean_v132B1,
  ,
  drop = FALSE
]

write.csv(
  metadata_clean,
  file.path(
    outdir,
    "16_v1.3.2B.1_clean_hepatobiliary_cell_metadata.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 30. Review flags
# ------------------------------------------------------------

biased_map <- manual_map[
  manual_map$sample_bias_ge70,
  ,
  drop = FALSE
]

review_flags <- c(
  "Human MASH v1.3.2B.1 Hepatobiliary Manual Annotation & Cleaning",
  "================================================================",
  "",
  paste0(
    "Input hepatobiliary-candidate cells: ",
    n_total
  ),
  paste0(
    "Clean hepatobiliary cells: ",
    n_clean,
    " (",
    round(
      100 *
        n_clean /
        n_total,
      2
    ),
    "%)"
  ),
  paste0(
    "Excluded non-hepatobiliary cells: ",
    n_excluded,
    " (",
    round(
      100 *
        n_excluded /
        n_total,
      2
    ),
    "%)"
  ),
  paste0(
    "Clean hepatocyte cells: ",
    ncol(
      hepatocyte_clean
    )
  ),
  paste0(
    "Clean cholangiocyte + ductular-progenitor cells: ",
    ncol(
      cholangiocyte_ductular_clean
    )
  ),
  "",
  "QC1:",
  "  Cluster 0 = Hepatocyte-reference closer: PASS",
  "  Cluster 7 = Hepatocyte-reference closer: PASS",
  "  Neither cluster 0 nor cluster 7 is reclassified as LSEC.",
  "",
  "FETAL HEPATOBLAST:",
  "  The previous Hepatoblast_or_immature_hepatocyte label is not retained.",
  "  v1.3.2B did not provide convincing fetal-hepatoblast support; DLK1 was absent.",
  "",
  "EXCLUDED CLUSTERS:",
  "  11 = immune/T-cell contamination",
  "  15 = mesenchymal/perivascular contamination",
  "  16 = macrophage contamination",
  "",
  "OUTSIDE-CORE DUCTAL CANDIDATES:",
  paste0(
    "  Audit-only candidates carried forward: ",
    ifelse(
      is.na(
        outside_core_count
      ),
      "not available",
      outside_core_count
    )
  ),
  "  These cells are NOT automatically rescued at v1.3.2B.1.",
  "",
  "STATUS:",
  "  MANUAL_CLEANED_REVIEW_REQUIRED",
  "  Canonical whole-liver v1.3 object was NOT modified."
)

if (
  nrow(
    biased_map
  ) > 0
) {

  review_flags <- c(
    review_flags,
    "",
    paste0(
      "FLAG: clusters with >=",
      sample_bias_threshold,
      "% contribution from one LT sample:"
    ),
    paste0(
      "  cluster ",
      biased_map$cluster,
      " -> ",
      biased_map$manual_state,
      " : ",
      biased_map$max_sample,
      " = ",
      round(
        biased_map$max_sample_contribution_pct,
        1
      ),
      "%"
    ),
    "  Retain as states, but do not treat these as robust cross-sample progression states without additional support."
  )
}

writeLines(
  review_flags,
  file.path(
    outdir,
    "18_v1.3.2B.1_review_flags.txt"
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
# 31. Save final v1.3.2B.1 objects
# ------------------------------------------------------------

manual_all_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatobiliary_v1.3.2B.1_manual_all.rds"
)

clean_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatobiliary_v1.3.2B.1_clean.rds"
)

hepatocyte_rds <- file.path(
  object_dir,
  "Human_MASH_Hepatocyte_v1.3.2B.1_clean.rds"
)

cholangiocyte_ductular_rds <- file.path(
  object_dir,
  "Human_MASH_Cholangiocyte_Ductular_v1.3.2B.1_clean.rds"
)

saveRDS(
  hepato_all,
  manual_all_rds
)

saveRDS(
  hepato_clean,
  clean_rds
)

saveRDS(
  hepatocyte_clean,
  hepatocyte_rds
)

saveRDS(
  cholangiocyte_ductular_clean,
  cholangiocyte_ductular_rds
)

# ------------------------------------------------------------
# 32. Manifest / provenance
# ------------------------------------------------------------

manifest <- c(
  "version=v1.3.2B.1",
  paste0(
    "input_checkpoint=",
    checkpoint_rds
  ),
  paste0(
    "input_checkpoint_md5=",
    unname(
      tools::md5sum(
        checkpoint_rds
      )
    )
  ),
  paste0(
    "qc1_decision_support=",
    qc1_decision_csv
  ),
  paste0(
    "qc1_decision_md5=",
    unname(
      tools::md5sum(
        qc1_decision_csv
      )
    )
  ),
  paste0(
    "input_cells=",
    n_total
  ),
  paste0(
    "clean_hepatobiliary_cells=",
    n_clean
  ),
  paste0(
    "excluded_nonhepatobiliary_cells=",
    n_excluded
  ),
  paste0(
    "clean_hepatocyte_cells=",
    ncol(
      hepatocyte_clean
    )
  ),
  paste0(
    "clean_cholangiocyte_plus_ductular_cells=",
    ncol(
      cholangiocyte_ductular_clean
    )
  ),
  paste0(
    "outside_core_ductal_candidates_audit_only=",
    ifelse(
      is.na(
        outside_core_count
      ),
      "NA",
      outside_core_count
    )
  ),
  "manual_annotation_resolution=0.6",
  "RPCA_recomputed=FALSE",
  "clustering_recomputed=FALSE",
  "FindAllMarkers_recomputed=FALSE",
  "QC1_cluster0_7_validation=PASS",
  paste0(
    "manual_all_object=",
    manual_all_rds
  ),
  paste0(
    "clean_hepatobiliary_object=",
    clean_rds
  ),
  paste0(
    "clean_hepatocyte_object=",
    hepatocyte_rds
  ),
  paste0(
    "clean_cholangiocyte_ductular_object=",
    cholangiocyte_ductular_rds
  ),
  "canonical_object_modified=FALSE",
  "annotation_status=MANUAL_CLEANED_REVIEW_REQUIRED"
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
cat("v1.3.2B.1 COMPLETE\n")
cat("HEPATOBILIARY MANUAL ANNOTATION & CLEANING COMPLETE\n")
cat("QC1 CLUSTERS 0/7 VALIDATION: PASS\n")
cat("ANNOTATION STATUS: MANUAL_CLEANED_REVIEW_REQUIRED\n")
cat("CANONICAL v1.3 OBJECT NOT MODIFIED\n")
cat("\n")
cat("Clean hepatobiliary object:\n")
cat(clean_rds, "\n")
cat("\n")
cat("Clean hepatocyte object:\n")
cat(hepatocyte_rds, "\n")
cat("\n")
cat("Clean cholangiocyte / ductular object:\n")
cat(cholangiocyte_ductular_rds, "\n")
cat("\n")
cat("Output directory:\n")
cat(outdir, "\n")
cat("============================================================\n")
