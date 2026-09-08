#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.4.3D.1"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

D0_MACHINERY <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3D.0",
  "analysis",
  "v1.4.3D.0_human_myeloid_IL10_conservation_audit",
  "tables",
  "v1.4.3D.0_06_IL10_response_machinery_by_state.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3D.1",
  "analysis",
  "v1.4.3D.1_human_myeloid_singlecell_IL10R_competence_audit"
)

TAB_OUT <- file.path(
  OUT,
  "tables"
)

LOG_OUT <- file.path(
  OUT,
  "logs"
)

dir.create(
  TAB_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOG_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

LOG_FILE <- file.path(
  LOG_OUT,
  paste0(
    VERSION,
    "_run.log"
  )
)

cat(
  "",
  file = LOG_FILE
)

msg <- function(...) {
  txt <- paste0(...)
  cat(
    txt,
    "\n"
  )
  cat(
    txt,
    "\n",
    file = LOG_FILE,
    append = TRUE
  )
}

stop2 <- function(...) {
  txt <- paste0(...)
  cat(
    "ERROR: ",
    txt,
    "\n",
    sep = "",
    file = LOG_FILE,
    append = TRUE
  )
  stop(
    txt,
    call. = FALSE
  )
}

write_tsv <- function(
  x,
  path
) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

resolve_feature <- function(
  requested_gene,
  available_features
) {
  idx_exact <- which(
    available_features ==
      requested_gene
  )

  if (
    length(
      idx_exact
    ) ==
      1L
  ) {
    return(
      available_features[
        idx_exact
      ]
    )
  }

  idx_ci <- which(
    toupper(
      available_features
    ) ==
      toupper(
        requested_gene
      )
  )

  if (
    length(
      idx_ci
    ) ==
      1L
  ) {
    return(
      available_features[
        idx_ci
      ]
    )
  }

  if (
    length(
      idx_ci
    ) >
      1L
  ) {
    stop2(
      "Ambiguous case-insensitive feature match for ",
      requested_gene,
      ": ",
      paste(
        available_features[
          idx_ci
        ],
        collapse = ", "
      )
    )
  }

  NA_character_
}

safe_mean <- function(
  x
) {
  x <- as.numeric(
    x
  )

  x <- x[
    is.finite(
      x
    )
  ]

  if (
    !length(
      x
    )
  ) {
    return(
      NA_real_
    )
  }

  mean(
    x
  )
}

safe_min <- function(
  x
) {
  x <- as.numeric(
    x
  )

  x <- x[
    is.finite(
      x
    )
  ]

  if (
    !length(
      x
    )
  ) {
    return(
      NA_real_
    )
  }

  min(
    x
  )
}

safe_max <- function(
  x
) {
  x <- as.numeric(
    x
  )

  x <- x[
    is.finite(
      x
    )
  ]

  if (
    !length(
      x
    )
  ) {
    return(
      NA_real_
    )
  }

  max(
    x
  )
}

msg(
  strrep(
    "=",
    76
  )
)

msg(
  "Human MASH ",
  VERSION,
  " | single-cell IL10R competence audit"
)

msg(
  strrep(
    "=",
    76
  )
)

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop2(
    "Final Human myeloid v1.9 RDS not found: ",
    INPUT_RDS
  )
}

suppressPackageStartupMessages({
  library(
    Seurat
  )
  library(
    SeuratObject
  )
  library(
    Matrix
  )
})

msg(
  "Seurat version: ",
  as.character(
    packageVersion(
      "Seurat"
    )
  )
)

msg(
  "SeuratObject version: ",
  as.character(
    packageVersion(
      "SeuratObject"
    )
  )
)

# ==============================================================================
# 1. Frozen design
# ==============================================================================

SAMPLE_ORDER <- c(
  "LT7",
  "LT8",
  "LT5",
  "LT6",
  "LT3",
  "LT4",
  "LT1",
  "LT2"
)

sample_info <- data.frame(
  sample =
    SAMPLE_ORDER,
  NAS = c(
    0,
    0,
    2,
    2,
    4,
    4,
    6,
    6
  ),
  fibrosis = c(
    "F0",
    "F0",
    "F0",
    "F0",
    "F0",
    "F1",
    "F4",
    "F3"
  ),
  stringsAsFactors = FALSE
)

STATE_ORDER <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

STATE_TIER <- c(
  "Resident/Kupffer-like_Mphi" =
    "PRIMARY_robust",
  "Repair/Resolution_Mphi" =
    "PRIMARY_robust",
  "Lipid-associated/TREM2_Mphi" =
    "EXPLORATORY_lower_abundance",
  "Inflammatory/Transitional_Mphi" =
    "PRIMARY_robust",
  "Inflammatory_Monocyte" =
    "EXPLORATORY_lower_abundance"
)

PRIMARY_STATES <- names(
  STATE_TIER
)[
  STATE_TIER ==
    "PRIMARY_robust"
]

AXIS_GENES <- c(
  "IL10RA",
  "IL10RB",
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3",
  "BCL3",
  "IL1RN"
)

DOWNSTREAM_GENES <- c(
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3",
  "BCL3",
  "IL1RN"
)

# ==============================================================================
# 2. Load frozen final-myeloid v1.9 object
# ==============================================================================

msg(
  "Loading final Human myeloid v1.9 RDS: ",
  INPUT_RDS
)

obj <- readRDS(
  INPUT_RDS
)

if (
  !inherits(
    obj,
    "Seurat"
  )
) {
  stop2(
    "Input object is not a Seurat object."
  )
}

required_meta <- c(
  "sample_v17",
  "NAS_v17",
  "myeloid_state_FINAL_v19"
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
  )
) {
  stop2(
    "Required frozen metadata missing: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  !"RNA" %in%
    Assays(
      obj
    )
) {
  stop2(
    "RNA assay absent."
  )
}

DefaultAssay(
  obj
) <- "RNA"

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

rna_layers <- Layers(
  obj[["RNA"]]
)

if (
  "counts" %in%
    rna_layers
) {
  counts_layer <- "counts"
} else {
  candidates <- grep(
    "^counts($|\\.)",
    rna_layers,
    value = TRUE
  )

  if (
    length(
      candidates
    ) !=
      1L
  ) {
    stop2(
      "Could not resolve unique RNA counts layer. Layers: ",
      paste(
        rna_layers,
        collapse = ", "
      )
    )
  }

  counts_layer <- candidates[
    1
  ]
}

counts <- LayerData(
  obj,
  assay = "RNA",
  layer = counts_layer
)

meta <- data.frame(
  cell =
    colnames(
      obj
    ),
  sample =
    as.character(
      obj@meta.data[["sample_v17"]]
    ),
  NAS =
    suppressWarnings(
      as.numeric(
        as.character(
          obj@meta.data[["NAS_v17"]]
        )
      )
    ),
  final_state =
    as.character(
      obj@meta.data[["myeloid_state_FINAL_v19"]]
    ),
  stringsAsFactors = FALSE
)

keep <- meta$sample %in%
  SAMPLE_ORDER &
  meta$final_state %in%
  STATE_ORDER

meta <- meta[
  keep,
  ,
  drop = FALSE
]

counts <- counts[
  ,
  keep,
  drop = FALSE
]

if (
  !identical(
    colnames(
      counts
    ),
    meta$cell
  )
) {
  stop2(
    "Expression matrix / metadata cell order mismatch."
  )
}

expected_nas <- setNames(
  sample_info$NAS,
  sample_info$sample
)

nas_check <- expected_nas[
  meta$sample
]

if (
  any(
    is.na(
      meta$NAS
    )
  ) ||
  any(
    meta$NAS !=
      as.numeric(
        nas_check
      )
  )
) {
  stop2(
    "sample_v17 and NAS_v17 are inconsistent with frozen LT1-LT8 design."
  )
}

meta$fibrosis <- sample_info$fibrosis[
  match(
    meta$sample,
    sample_info$sample
  )
]

meta$analysis_tier <- unname(
  STATE_TIER[
    meta$final_state
  ]
)

msg(
  "Retained final Human myeloid cells: ",
  nrow(
    meta
  )
)

# ==============================================================================
# 3. Feature audit
# ==============================================================================

features <- rownames(
  counts
)

feature_map <- data.frame(
  requested_gene =
    AXIS_GENES,
  actual_feature =
    vapply(
      AXIS_GENES,
      function(
        g
      ) {
        resolve_feature(
          g,
          features
        )
      },
      character(
        1
      )
    ),
  stringsAsFactors = FALSE
)

feature_map$present <-
  !is.na(
    feature_map$actual_feature
  )

write_tsv(
  feature_map,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_01_feature_audit.tsv"
    )
  )
)

if (
  !all(
    c(
      "IL10RA",
      "IL10RB"
    ) %in%
      feature_map$requested_gene[
        feature_map$present
      ]
  )
) {
  stop2(
    "IL10RA and/or IL10RB are absent from single-cell RNA counts."
  )
}

# ==============================================================================
# 4. Build cell-level binary detection matrix for IL10 axis
# ==============================================================================

detected_genes <- feature_map$requested_gene[
  feature_map$present
]

actual_genes <- feature_map$actual_feature[
  feature_map$present
]

axis_counts <- counts[
  actual_genes,
  ,
  drop = FALSE
]

rownames(
  axis_counts
) <- detected_genes

axis_detect <- axis_counts >
  0

il10ra_pos <- as.vector(
  axis_detect[
    "IL10RA",
    ,
    drop = TRUE
  ]
)

il10rb_pos <- as.vector(
  axis_detect[
    "IL10RB",
    ,
    drop = TRUE
  ]
)

meta$IL10RA_positive <-
  il10ra_pos

meta$IL10RB_positive <-
  il10rb_pos

meta$IL10R_dual_positive <-
  il10ra_pos &
  il10rb_pos

meta$IL10R_any_positive <-
  il10ra_pos |
  il10rb_pos

meta$IL10R_status <- ifelse(
  meta$IL10R_dual_positive,
  "dual_positive",
  ifelse(
    meta$IL10RA_positive,
    "IL10RA_only",
    ifelse(
      meta$IL10RB_positive,
      "IL10RB_only",
      "double_negative"
    )
  )
)

# ==============================================================================
# 5. Frozen state/sample cell counts
# ==============================================================================

cell_count_table <- as.data.frame(
  table(
    sample =
      factor(
        meta$sample,
        levels = SAMPLE_ORDER
      ),
    final_state =
      factor(
        meta$final_state,
        levels = STATE_ORDER
      )
  ),
  stringsAsFactors = FALSE
)

names(
  cell_count_table
) <- c(
  "sample",
  "final_state",
  "n_cells"
)

cell_count_table <- cell_count_table[
  cell_count_table$n_cells >
    0,
  ,
  drop = FALSE
]

cell_count_table$NAS <- sample_info$NAS[
  match(
    cell_count_table$sample,
    sample_info$sample
  )
]

cell_count_table$fibrosis <- sample_info$fibrosis[
  match(
    cell_count_table$sample,
    sample_info$sample
  )
]

cell_count_table$analysis_tier <- unname(
  STATE_TIER[
    cell_count_table$final_state
  ]
)

write_tsv(
  cell_count_table,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_02_state_sample_cell_counts.tsv"
    )
  )
)

# ==============================================================================
# 6. IL10R single-cell competence by state x biological sample
# ==============================================================================

group_key <- interaction(
  meta$sample,
  meta$final_state,
  drop = TRUE,
  sep = "|||"
)

group_idx <- split(
  seq_len(
    nrow(
      meta
    )
  ),
  group_key
)

competence_rows <- list()
cri <- 1L

for (
  idx in group_idx
) {
  z <- meta[
    idx,
    ,
    drop = FALSE
  ]

  n_cells <- nrow(
    z
  )

  dual_n <- sum(
    z$IL10R_dual_positive
  )

  any_n <- sum(
    z$IL10R_any_positive
  )

  competence_rows[[cri]] <- data.frame(
    sample =
      z$sample[
        1
      ],
    NAS =
      z$NAS[
        1
      ],
    fibrosis =
      z$fibrosis[
        1
      ],
    final_state =
      z$final_state[
        1
      ],
    analysis_tier =
      z$analysis_tier[
        1
      ],
    n_cells =
      n_cells,
    IL10RA_positive_cells =
      sum(
        z$IL10RA_positive
      ),
    IL10RA_positive_fraction =
      mean(
        z$IL10RA_positive
      ),
    IL10RB_positive_cells =
      sum(
        z$IL10RB_positive
      ),
    IL10RB_positive_fraction =
      mean(
        z$IL10RB_positive
      ),
    IL10R_dual_positive_cells =
      dual_n,
    IL10R_dual_positive_fraction =
      dual_n /
      n_cells,
    IL10R_any_positive_cells =
      any_n,
    IL10R_any_positive_fraction =
      any_n /
      n_cells,
    dual_fraction_among_receptor_detected =
      if (
        any_n >
          0
      ) {
        dual_n /
          any_n
      } else {
        NA_real_
      },
    min20_cells =
      n_cells >=
      20L,
    stringsAsFactors = FALSE
  )

  cri <- cri + 1L
}

competence_by_sample <- do.call(
  rbind,
  competence_rows
)

competence_by_sample <- competence_by_sample[
  order(
    match(
      competence_by_sample$final_state,
      STATE_ORDER
    ),
    match(
      competence_by_sample$sample,
      SAMPLE_ORDER
    )
  ),
  ,
  drop = FALSE
]

write_tsv(
  competence_by_sample,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_03_IL10R_competence_by_state_sample.tsv"
    )
  )
)

# ==============================================================================
# 7. Four-category IL10R status distribution
# ==============================================================================

status_table <- as.data.frame(
  table(
    sample =
      factor(
        meta$sample,
        levels = SAMPLE_ORDER
      ),
    final_state =
      factor(
        meta$final_state,
        levels = STATE_ORDER
      ),
    IL10R_status =
      factor(
        meta$IL10R_status,
        levels = c(
          "dual_positive",
          "IL10RA_only",
          "IL10RB_only",
          "double_negative"
        )
      )
  ),
  stringsAsFactors = FALSE
)

names(
  status_table
) <- c(
  "sample",
  "final_state",
  "IL10R_status",
  "n_cells"
)

status_table <- status_table[
  status_table$n_cells >
    0,
  ,
  drop = FALSE
]

status_total <- aggregate(
  n_cells ~ sample + final_state,
  data = status_table,
  sum
)

names(
  status_total
)[3] <- "state_sample_total"

status_table <- merge(
  status_table,
  status_total,
  by = c(
    "sample",
    "final_state"
  ),
  all.x = TRUE,
  sort = FALSE
)

status_table$fraction <-
  status_table$n_cells /
  status_table$state_sample_total

status_table$NAS <- sample_info$NAS[
  match(
    status_table$sample,
    sample_info$sample
  )
]

write_tsv(
  status_table,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_04_IL10R_status_distribution.tsv"
    )
  )
)

# ==============================================================================
# 8. Downstream machinery gene detection by state x sample
# ==============================================================================

axis_gene_rows <- list()
agi <- 1L

for (
  idx in group_idx
) {
  zmeta <- meta[
    idx,
    ,
    drop = FALSE
  ]

  n_cells <- nrow(
    zmeta
  )

  for (
    g in detected_genes
  ) {
    pos <- as.vector(
      axis_detect[
        g,
        idx,
        drop = TRUE
      ]
    )

    umi <- as.numeric(
      axis_counts[
        g,
        idx,
        drop = TRUE
      ]
    )

    axis_gene_rows[[agi]] <- data.frame(
      sample =
        zmeta$sample[
          1
        ],
      NAS =
        zmeta$NAS[
          1
        ],
      fibrosis =
        zmeta$fibrosis[
          1
        ],
      final_state =
        zmeta$final_state[
          1
        ],
      analysis_tier =
        zmeta$analysis_tier[
          1
        ],
      gene =
        g,
      n_cells =
        n_cells,
      positive_cells =
        sum(
          pos
        ),
      positive_fraction =
        mean(
          pos
        ),
      total_UMI =
        sum(
          umi
        ),
      UMI_per_1000_cells =
        sum(
          umi
        ) /
        n_cells *
        1000,
      stringsAsFactors = FALSE
    )

    agi <- agi + 1L
  }
}

axis_gene_by_sample <- do.call(
  rbind,
  axis_gene_rows
)

write_tsv(
  axis_gene_by_sample,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_05_IL10_axis_gene_detection_by_state_sample.tsv"
    )
  )
)

# ==============================================================================
# 9. Descriptive NAS-stage summaries
# ==============================================================================

stage_rows <- list()
sri <- 1L

for (
  st in STATE_ORDER
) {
  z <- competence_by_sample[
    competence_by_sample$final_state ==
      st,
    ,
    drop = FALSE
  ]

  for (
    nas_now in c(
      0,
      2,
      4,
      6
    )
  ) {
    q <- z[
      z$NAS ==
        nas_now,
      ,
      drop = FALSE
    ]

    stage_rows[[sri]] <- data.frame(
      final_state =
        st,
      analysis_tier =
        unname(
          STATE_TIER[[st]]
        ),
      NAS =
        nas_now,
      n_biological_samples =
        nrow(
          q
        ),
      mean_IL10RA_positive_fraction =
        safe_mean(
          q$IL10RA_positive_fraction
        ),
      mean_IL10RB_positive_fraction =
        safe_mean(
          q$IL10RB_positive_fraction
        ),
      mean_IL10R_dual_positive_fraction =
        safe_mean(
          q$IL10R_dual_positive_fraction
        ),
      min_IL10R_dual_positive_fraction =
        safe_min(
          q$IL10R_dual_positive_fraction
        ),
      max_IL10R_dual_positive_fraction =
        safe_max(
          q$IL10R_dual_positive_fraction
        ),
      stringsAsFactors = FALSE
    )

    sri <- sri + 1L
  }
}

stage_summary <- do.call(
  rbind,
  stage_rows
)

write_tsv(
  stage_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_06_IL10R_competence_by_state_NAS_stage.tsv"
    )
  )
)

# ==============================================================================
# 10. Advanced NAS6/F3-F4 versus NAS0/F0 descriptive bridge
# ==============================================================================

advanced_rows <- list()
advi <- 1L

for (
  st in STATE_ORDER
) {
  z <- competence_by_sample[
    competence_by_sample$final_state ==
      st,
    ,
    drop = FALSE
  ]

  b <- z[
    z$NAS ==
      0,
    ,
    drop = FALSE
  ]

  a <- z[
    z$NAS ==
      6,
    ,
    drop = FALSE
  ]

  baseline_mean <- safe_mean(
    b$IL10R_dual_positive_fraction
  )

  advanced_mean <- safe_mean(
    a$IL10R_dual_positive_fraction
  )

  advanced_rows[[advi]] <- data.frame(
    final_state =
      st,
    analysis_tier =
      unname(
        STATE_TIER[[st]]
      ),
    LT7_NAS0_F0_dual_fraction =
      z$IL10R_dual_positive_fraction[
        z$sample ==
          "LT7"
      ],
    LT8_NAS0_F0_dual_fraction =
      z$IL10R_dual_positive_fraction[
        z$sample ==
          "LT8"
      ],
    LT1_NAS6_F4_dual_fraction =
      z$IL10R_dual_positive_fraction[
        z$sample ==
          "LT1"
      ],
    LT2_NAS6_F3_dual_fraction =
      z$IL10R_dual_positive_fraction[
        z$sample ==
          "LT2"
      ],
    NAS0_mean_dual_fraction =
      baseline_mean,
    NAS6_F3F4_mean_dual_fraction =
      advanced_mean,
    advanced_minus_baseline =
      advanced_mean -
      baseline_mean,
    both_advanced_samples_detect_dual_positive_cells =
      if (
        nrow(
          a
        ) ==
          2L
      ) {
        all(
          a$IL10R_dual_positive_cells >
            0
        )
      } else {
        FALSE
      },
    both_advanced_samples_above_zero_fraction =
      if (
        nrow(
          a
        ) ==
          2L
      ) {
        all(
          a$IL10R_dual_positive_fraction >
            0
        )
      } else {
        FALSE
      },
    both_advanced_samples_above_NAS0_mean =
      if (
        nrow(
          a
        ) ==
          2L &&
        is.finite(
          baseline_mean
        )
      ) {
        all(
          a$IL10R_dual_positive_fraction >
            baseline_mean
        )
      } else {
        FALSE
      },
    interpretation_boundary =
      "Detected dual positivity in probe-based single-cell RNA counts; zero is not proof of biological absence; NAS6 is confounded with F3/F4 fibrosis",
    stringsAsFactors = FALSE
  )

  advi <- advi + 1L
}

advanced_summary <- do.call(
  rbind,
  advanced_rows
)

write_tsv(
  advanced_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_07_advanced_vs_baseline_IL10R_bridge.tsv"
    )
  )
)

# ==============================================================================
# 11. Primary-state compact translational summary
# ==============================================================================

primary_summary <- advanced_summary[
  advanced_summary$final_state %in%
    PRIMARY_STATES,
  ,
  drop = FALSE
]

primary_summary$translational_call <- ifelse(
  primary_summary$both_advanced_samples_detect_dual_positive_cells,
  "same_cell_IL10RA_IL10RB_co_detection_present_in_both_advanced_samples",
  "same_cell_dual_detection_not_replicated_in_both_advanced_samples"
)

write_tsv(
  primary_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_08_primary_state_translational_summary.tsv"
    )
  )
)

# ==============================================================================
# 12. D0 state-level receptor-availability consistency check
# ==============================================================================

d0_validation <- data.frame(
  status =
    "D0 machinery table unavailable; validation skipped",
  states_checked =
    0L,
  mismatches =
    NA_integer_,
  stringsAsFactors = FALSE
)

if (
  file.exists(
    D0_MACHINERY
  )
) {
  d0 <- read.delim(
    D0_MACHINERY,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (
    all(
      c(
        "final_state",
        "IL10RA_present",
        "IL10RB_present"
      ) %in%
        colnames(
          d0
        )
    )
  ) {
    observed <- data.frame(
      final_state =
        STATE_ORDER,
      D1_IL10RA_present =
        "IL10RA" %in%
        detected_genes,
      D1_IL10RB_present =
        "IL10RB" %in%
        detected_genes,
      stringsAsFactors = FALSE
    )

    chk <- merge(
      d0[
        ,
        c(
          "final_state",
          "IL10RA_present",
          "IL10RB_present"
        ),
        drop = FALSE
      ],
      observed,
      by = "final_state",
      all = TRUE,
      sort = FALSE
    )

    chk$mismatch <-
      chk$IL10RA_present !=
      chk$D1_IL10RA_present |
      chk$IL10RB_present !=
      chk$D1_IL10RB_present

    mismatch_n <- sum(
      chk$mismatch,
      na.rm = TRUE
    )

    d0_validation <- data.frame(
      status =
        if (
          mismatch_n ==
            0L
        ) {
          "PASS_exact_state_level_feature_availability"
        } else {
          "WARNING_feature_availability_mismatch"
        },
      states_checked =
        nrow(
          chk
        ),
      mismatches =
        mismatch_n,
      stringsAsFactors = FALSE
    )

    write_tsv(
      chk,
      file.path(
        TAB_OUT,
        paste0(
          VERSION,
          "_09_D0_feature_availability_validation.tsv"
        )
      )
    )
  }
}

write_tsv(
  d0_validation,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_10_D0_validation_summary.tsv"
    )
  )
)

# ==============================================================================
# 13. Provenance and interpretation boundaries
# ==============================================================================

provenance <- data.frame(
  key = c(
    "version",
    "input_RDS",
    "sample_metadata",
    "NAS_metadata",
    "final_state_metadata",
    "RNA_counts_layer",
    "single_cell_competence_definition",
    "assay_boundary",
    "disease_stage_boundary",
    "causal_boundary",
    "statistics_boundary"
  ),
  value = c(
    VERSION,
    INPUT_RDS,
    "sample_v17",
    "NAS_v17",
    "myeloid_state_FINAL_v19",
    counts_layer,
    "IL10R dual-positive = IL10RA count > 0 AND IL10RB count > 0 in the same cell",
    "Human data are 10x Flex/probe-based scRNA-seq; detected dual positivity supports same-cell transcript co-detection, while zero counts can reflect limited capture and are not proof of biological absence",
    "LT1/LT2 are NAS6/F4-F3; advanced disease is confounded with fibrosis and must not be described as a pure NAS effect",
    "Same-cell IL10RA/IL10RB transcript detection demonstrates receptor-transcript competence, not receptor protein, ligand exposure, receptor occupancy, JAK activation, or pSTAT3",
    "Biological samples are the interpretive unit; stage summaries are descriptive, with n=2 samples per NAS stage and no cell-level inferential p-values"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  provenance,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_11_provenance_and_boundaries.tsv"
    )
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_OUT,
    paste0(
      VERSION,
      "_sessionInfo.txt"
    )
  )
)

# ==============================================================================
# 14. Console output
# ==============================================================================

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "FEATURE AUDIT"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  feature_map,
  row.names = FALSE
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "IL10R DUAL-POSITIVE FRACTION | PRIMARY HUMAN MYELOID STATES"
)

msg(
  strrep(
    "=",
    76
  )
)

primary_competence <- competence_by_sample[
  competence_by_sample$final_state %in%
    PRIMARY_STATES,
  c(
    "final_state",
    "sample",
    "NAS",
    "fibrosis",
    "n_cells",
    "IL10RA_positive_fraction",
    "IL10RB_positive_fraction",
    "IL10R_dual_positive_fraction",
    "dual_fraction_among_receptor_detected"
  ),
  drop = FALSE
]

print(
  primary_competence,
  row.names = FALSE
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "ADVANCED NAS6/F3-F4 VS NAS0/F0 | ALL FINAL MYELOID STATES"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  advanced_summary,
  row.names = FALSE
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "PRIMARY TRANSLATIONAL SUMMARY"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  primary_summary,
  row.names = FALSE
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "D0 FEATURE-AVAILABILITY VALIDATION"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  d0_validation,
  row.names = FALSE
)

capture.output(
  print(
    feature_map,
    row.names = FALSE
  ),
  print(
    primary_competence,
    row.names = FALSE
  ),
  print(
    advanced_summary,
    row.names = FALSE
  ),
  print(
    primary_summary,
    row.names = FALSE
  ),
  print(
    d0_validation,
    row.names = FALSE
  ),
  file = LOG_FILE,
  append = TRUE
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "OUTPUT"
)

msg(
  strrep(
    "=",
    76
  )
)

msg(
  OUT
)

msg("")

msg(
  "PASS: FINAL HUMAN MYELOID v1.9 OBJECT LOADED"
)

msg(
  "PASS: SAME-CELL IL10RA/IL10RB CO-DETECTION AUDITED"
)

msg(
  "PASS: IL10 DOWNSTREAM TRANSCRIPT DETECTION AUDITED BY STATE AND BIOLOGICAL SAMPLE"
)

msg(
  "PASS: ADVANCED NAS6/F3-F4 IL10R COMPETENCE COMPARED DESCRIPTIVELY WITH NAS0/F0"
)

msg(
  "PASS: ",
  VERSION,
  " HUMAN MYELOID SINGLE-CELL IL10R COMPETENCE AUDIT COMPLETE"
)

msg(
  strrep(
    "=",
    76
  )
)
