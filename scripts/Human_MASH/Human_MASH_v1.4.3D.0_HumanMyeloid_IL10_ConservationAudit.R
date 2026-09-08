#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v1.4.3D.0"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

REPO <-
  "/Users/uenoya/Projects/uenoyscRNA"

DE_DIR <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloid_PseudobulkDE_v2.2.1"
)

LOGCPM_RDS <- file.path(
  DE_DIR,
  "objects",
  "Human_MASH_LT1to8_FinalMyeloid_TMM_logCPM_v2.2.1.rds"
)

TREND_CSV <- file.path(
  DE_DIR,
  "tables",
  "02_all_states_overall_severity_trend_v2.2.1.csv"
)

STAGE_PATTERN_CSV <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloid_StagePatternAudit_v2.3",
  "tables",
  "01_all_gene_stage_pattern_table_v2.3.csv"
)

MOUSE_PROGRAM_SOURCE <- file.path(
  REPO,
  "analysis",
  "Mouse_MASH_Mphi_FINAL_extended_figures_v4.16.1.R"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3D.0",
  "analysis",
  "v1.4.3D.0_human_myeloid_IL10_conservation_audit"
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

safe_z <- function(
  x
) {
  x <- as.numeric(
    x
  )

  s <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      s
    ) ||
      s == 0
  ) {
    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) /
    s
}

extract_top_level_assignment <- function(
  exprs,
  object_name,
  env
) {
  found <- FALSE

  for (
    e in exprs
  ) {
    if (
      !is.call(
        e
      ) ||
        length(
          e
        ) < 3L
    ) {
      next
    }

    op <- as.character(
      e[[1]]
    )

    if (
      !op %in%
        c(
          "<-",
          "="
        )
    ) {
      next
    }

    lhs <- e[[2]]

    if (
      !is.symbol(
        lhs
      )
    ) {
      next
    }

    if (
      as.character(
        lhs
      ) !=
        object_name
    ) {
      next
    }

    eval(
      e,
      envir = env
    )

    found <- TRUE
    break
  }

  if (
    !found
  ) {
    stop2(
      "Could not extract top-level assignment `",
      object_name,
      "` from ",
      MOUSE_PROGRAM_SOURCE
    )
  }
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

msg(
  strrep(
    "=",
    76
  )
)

msg(
  "Human MASH ",
  VERSION,
  " | Human myeloid IL10 conservation audit"
)

msg(
  strrep(
    "=",
    76
  )
)

required_files <- c(
  LOGCPM_RDS,
  TREND_CSV,
  STAGE_PATTERN_CSV,
  MOUSE_PROGRAM_SOURCE
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing_files
  )
) {
  stop2(
    "Required files missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ==============================================================================
# 1. Sample/state design copied from validated v2.2.1 / v2.3
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

sample_info$NAS_step <-
  sample_info$NAS /
  2

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

# ==============================================================================
# 2. Load validated human final-myeloid pseudobulk results
# ==============================================================================

logcpm_list <- readRDS(
  LOGCPM_RDS
)

if (
  !all(
    STATE_ORDER %in%
      names(
        logcpm_list
      )
  )
) {
  stop2(
    "logCPM RDS missing expected state(s): ",
    paste(
      setdiff(
        STATE_ORDER,
        names(
          logcpm_list
        )
      ),
      collapse = " | "
    )
  )
}

for (
  st in STATE_ORDER
) {
  mat <- logcpm_list[[st]]

  if (
    !all(
      SAMPLE_ORDER %in%
        colnames(
          mat
        )
    )
  ) {
    stop2(
      "State missing one or more required samples: ",
      st
    )
  }

  logcpm_list[[st]] <-
    mat[
      ,
      SAMPLE_ORDER,
      drop = FALSE
    ]
}

trend <- read.csv(
  TREND_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage_pattern <- read.csv(
  STAGE_PATTERN_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_trend_cols <- c(
  "final_state",
  "gene",
  "logFC_per_severity_step_NAS2",
  "PValue",
  "FDR"
)

missing_trend_cols <- setdiff(
  required_trend_cols,
  colnames(
    trend
  )
)

if (
  length(
    missing_trend_cols
  )
) {
  stop2(
    "Trend table missing columns: ",
    paste(
      missing_trend_cols,
      collapse = ", "
    )
  )
}

required_pattern_cols <- c(
  "final_state",
  "gene",
  "mean_NAS0_F0",
  "mean_NAS2_F0",
  "mean_NAS4_mixed_F0_F1",
  "mean_NAS6_F3_F4",
  "delta_NAS6_vs_NAS0",
  "logFC_per_severity_step_NAS2",
  "pattern_class"
)

missing_pattern_cols <- setdiff(
  required_pattern_cols,
  colnames(
    stage_pattern
  )
)

if (
  length(
    missing_pattern_cols
  )
) {
  stop2(
    "Stage-pattern table missing columns: ",
    paste(
      missing_pattern_cols,
      collapse = ", "
    )
  )
}

trend$gene_upper <- toupper(
  trend$gene
)

stage_pattern$gene_upper <- toupper(
  stage_pattern$gene
)

# ==============================================================================
# 3. Recover exact mouse v4.16.1 programs and map one-to-one symbols to human
# ==============================================================================

mouse_exprs <- parse(
  file = MOUSE_PROGRAM_SOURCE,
  keep.source = FALSE
)

mouse_env <- new.env(
  parent = baseenv()
)

extract_top_level_assignment(
  mouse_exprs,
  "PROGRAMS",
  mouse_env
)

MOUSE_PROGRAMS <- mouse_env$PROGRAMS

preferred_programs <- c(
  "IL10_STAT3",
  "Anti_inflammatory",
  "Repair_Resolution",
  "Efferocytosis",
  "Inflammatory",
  "ECM_associated_inflammatory",
  "Lipid_TREM2"
)

missing_programs <- setdiff(
  preferred_programs,
  names(
    MOUSE_PROGRAMS
  )
)

if (
  length(
    missing_programs
  )
) {
  stop2(
    "Expected mouse programs missing from v4.16.1: ",
    paste(
      missing_programs,
      collapse = ", "
    )
  )
}

HUMAN_PROGRAMS <- lapply(
  MOUSE_PROGRAMS[
    preferred_programs
  ],
  function(
    genes
  ) {
    unique(
      toupper(
        genes
      )
    )
  }
)

HUMAN_PROGRAMS[["IL10_downstream_no_receptors"]] <-
  setdiff(
    HUMAN_PROGRAMS[["IL10_STAT3"]],
    c(
      "IL10RA",
      "IL10RB"
    )
  )

PROGRAM_ORDER <- c(
  "IL10_STAT3",
  "IL10_downstream_no_receptors",
  "Anti_inflammatory",
  "Repair_Resolution",
  "Efferocytosis",
  "Inflammatory",
  "ECM_associated_inflammatory",
  "Lipid_TREM2"
)

IL10_AXIS_GENES <- c(
  "IL10RA",
  "IL10RB",
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3",
  "BCL3",
  "IL1RN"
)

IL10_DOWNSTREAM_GENES <- c(
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3",
  "BCL3",
  "IL1RN"
)

program_manifest <- do.call(
  rbind,
  lapply(
    PROGRAM_ORDER,
    function(
      program_now
    ) {
      data.frame(
        program =
          program_now,
        human_gene =
          HUMAN_PROGRAMS[[program_now]],
        mapping =
          "mouse_v4.16.1_symbol_uppercase_one_to_one_candidate",
        stringsAsFactors = FALSE
      )
    }
  )
)

write_tsv(
  program_manifest,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_01_crossspecies_program_manifest.tsv"
    )
  )
)

# ==============================================================================
# 4. Feature/program coverage by final human myeloid state
# ==============================================================================

audit_rows <- list()
ari <- 1L

for (
  st in STATE_ORDER
) {
  features <- rownames(
    logcpm_list[[st]]
  )

  requested <- unique(
    c(
      IL10_AXIS_GENES,
      unlist(
        HUMAN_PROGRAMS,
        use.names = FALSE
      )
    )
  )

  for (
    g in requested
  ) {
    actual <- resolve_feature(
      g,
      features
    )

    audit_rows[[ari]] <-
      data.frame(
        final_state =
          st,
        analysis_tier =
          unname(
            STATE_TIER[[st]]
          ),
        requested_gene =
          g,
        actual_feature =
          if (
            is.na(
              actual
            )
          ) {
            NA_character_
          } else {
            actual
          },
        present_in_filtered_pseudobulk =
          !is.na(
            actual
          ),
        stringsAsFactors = FALSE
      )

    ari <- ari + 1L
  }
}

feature_audit <- do.call(
  rbind,
  audit_rows
)

write_tsv(
  feature_audit,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_02_feature_audit_by_state.tsv"
    )
  )
)

coverage_rows <- list()
cri <- 1L

for (
  st in STATE_ORDER
) {
  features <- rownames(
    logcpm_list[[st]]
  )

  for (
    program_now in PROGRAM_ORDER
  ) {
    genes <- HUMAN_PROGRAMS[[program_now]]

    detected <- vapply(
      genes,
      function(
        g
      ) {
        !is.na(
          resolve_feature(
            g,
            features
          )
        )
      },
      logical(
        1
      )
    )

    coverage_rows[[cri]] <-
      data.frame(
        final_state =
          st,
        analysis_tier =
          unname(
            STATE_TIER[[st]]
          ),
        program =
          program_now,
        requested_genes =
          length(
            genes
          ),
        detected_genes =
          sum(
            detected
          ),
        coverage_fraction =
          mean(
            detected
          ),
        missing_genes =
          paste(
            genes[
              !detected
            ],
            collapse = ";"
          ),
        stringsAsFactors = FALSE
      )

    cri <- cri + 1L
  }
}

program_coverage <- do.call(
  rbind,
  coverage_rows
)

write_tsv(
  program_coverage,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_03_program_coverage_by_state.tsv"
    )
  )
)

# ==============================================================================
# 5. IL10-axis gene trajectories by state
# ==============================================================================

axis_sample_rows <- list()
asi <- 1L

axis_stage_rows <- list()
asti <- 1L

for (
  st in STATE_ORDER
) {
  mat <- logcpm_list[[st]]
  features <- rownames(
    mat
  )

  for (
    g in IL10_AXIS_GENES
  ) {
    actual <- resolve_feature(
      g,
      features
    )

    if (
      is.na(
        actual
      )
    ) {
      next
    }

    x <- as.numeric(
      mat[
        actual,
        SAMPLE_ORDER
      ]
    )

    names(
      x
    ) <- SAMPLE_ORDER

    for (
      i in seq_along(
        SAMPLE_ORDER
      )
    ) {
      axis_sample_rows[[asi]] <-
        data.frame(
          final_state =
            st,
          analysis_tier =
            unname(
              STATE_TIER[[st]]
            ),
          gene =
            g,
          actual_feature =
            actual,
          sample =
            SAMPLE_ORDER[i],
          NAS =
            sample_info$NAS[i],
          fibrosis =
            sample_info$fibrosis[i],
          logCPM =
            x[
              SAMPLE_ORDER[i]
            ],
          stringsAsFactors = FALSE
        )

      asi <- asi + 1L
    }

    tr <- trend[
      trend$final_state ==
        st &
      trend$gene_upper ==
        g,
      ,
      drop = FALSE
    ]

    pt <- stage_pattern[
      stage_pattern$final_state ==
        st &
      stage_pattern$gene_upper ==
        g,
      ,
      drop = FALSE
    ]

    axis_stage_rows[[asti]] <-
      data.frame(
        final_state =
          st,
        analysis_tier =
          unname(
            STATE_TIER[[st]]
          ),
        gene =
          g,
        mean_NAS0_F0 =
          safe_mean(
            x[
              c(
                "LT7",
                "LT8"
              )
            ]
          ),
        mean_NAS2_F0 =
          safe_mean(
            x[
              c(
                "LT5",
                "LT6"
              )
            ]
          ),
        mean_NAS4_mixed_F0_F1 =
          safe_mean(
            x[
              c(
                "LT3",
                "LT4"
              )
            ]
          ),
        mean_NAS6_F3_F4 =
          safe_mean(
            x[
              c(
                "LT1",
                "LT2"
              )
            ]
          ),
        delta_NAS6_vs_NAS0 =
          safe_mean(
            x[
              c(
                "LT1",
                "LT2"
              )
            ]
          ) -
          safe_mean(
            x[
              c(
                "LT7",
                "LT8"
              )
            ]
          ),
        LT1_NAS6_F4_logCPM =
          x[["LT1"]],
        LT2_NAS6_F3_logCPM =
          x[["LT2"]],
        logFC_per_severity_step_NAS2 =
          if (
            nrow(
              tr
            ) ==
              1L
          ) {
            tr$logFC_per_severity_step_NAS2
          } else {
            NA_real_
          },
        PValue =
          if (
            nrow(
              tr
            ) ==
              1L
          ) {
            tr$PValue
          } else {
            NA_real_
          },
        FDR =
          if (
            nrow(
              tr
            ) ==
              1L
          ) {
            tr$FDR
          } else {
            NA_real_
          },
        pattern_class =
          if (
            nrow(
              pt
            ) ==
              1L
          ) {
            pt$pattern_class
          } else {
            NA_character_
          },
        source_boundary =
          "TMM_logCPM_v2.2.1; advanced NAS6 samples are F3/F4; not a pure NAS effect",
        stringsAsFactors = FALSE
      )

    asti <- asti + 1L
  }
}

axis_by_sample <- do.call(
  rbind,
  axis_sample_rows
)

axis_stage_summary <- do.call(
  rbind,
  axis_stage_rows
)

write_tsv(
  axis_by_sample,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_04_IL10_axis_logCPM_by_state_sample.tsv"
    )
  )
)

write_tsv(
  axis_stage_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_05_IL10_axis_stage_summary.tsv"
    )
  )
)

# ==============================================================================
# 6. State-level IL10 response-machinery availability
# ==============================================================================

machinery_rows <- list()
mri <- 1L

for (
  st in STATE_ORDER
) {
  z <- axis_stage_summary[
    axis_stage_summary$final_state ==
      st,
    ,
    drop = FALSE
  ]

  has_gene <- function(
    g
  ) {
    g %in%
      z$gene
  }

  getv <- function(
    g,
    col
  ) {
    x <- z[
      z$gene ==
        g,
      col,
      drop = TRUE
    ]

    if (
      length(
        x
      ) ==
        1L
    ) {
      return(
        x
      )
    }

    NA_real_
  }

  downstream_present <- vapply(
    IL10_DOWNSTREAM_GENES,
    has_gene,
    logical(
      1
    )
  )

  machinery_rows[[mri]] <-
    data.frame(
      final_state =
        st,
      analysis_tier =
        unname(
          STATE_TIER[[st]]
        ),
      IL10RA_present =
        has_gene(
          "IL10RA"
        ),
      IL10RB_present =
        has_gene(
          "IL10RB"
        ),
      receptor_pair_present_at_state_level =
        has_gene(
          "IL10RA"
        ) &&
        has_gene(
          "IL10RB"
        ),
      downstream_requested_genes =
        length(
          IL10_DOWNSTREAM_GENES
        ),
      downstream_detected_genes =
        sum(
          downstream_present
        ),
      downstream_coverage_fraction =
        mean(
          downstream_present
        ),
      IL10RA_NAS0_mean_logCPM =
        getv(
          "IL10RA",
          "mean_NAS0_F0"
        ),
      IL10RA_NAS6_mean_logCPM =
        getv(
          "IL10RA",
          "mean_NAS6_F3_F4"
        ),
      IL10RA_NAS6_minus_NAS0 =
        getv(
          "IL10RA",
          "delta_NAS6_vs_NAS0"
        ),
      IL10RB_NAS0_mean_logCPM =
        getv(
          "IL10RB",
          "mean_NAS0_F0"
        ),
      IL10RB_NAS6_mean_logCPM =
        getv(
          "IL10RB",
          "mean_NAS6_F3_F4"
        ),
      IL10RB_NAS6_minus_NAS0 =
        getv(
          "IL10RB",
          "delta_NAS6_vs_NAS0"
        ),
      interpretation_boundary =
        "State-level pseudobulk transcript availability only; does not demonstrate IL10RA/IL10RB same-cell coexpression, receptor occupancy, or STAT3 phosphorylation",
      stringsAsFactors = FALSE
    )

  mri <- mri + 1L
}

machinery_summary <- do.call(
  rbind,
  machinery_rows
)

write_tsv(
  machinery_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_06_IL10_response_machinery_by_state.tsv"
    )
  )
)

# ==============================================================================
# 7. Cross-species program scores in human final-myeloid pseudobulk
# ==============================================================================

program_sample_rows <- list()
psi <- 1L

program_stage_rows <- list()
psti <- 1L

for (
  st in STATE_ORDER
) {
  mat <- logcpm_list[[st]]
  features <- rownames(
    mat
  )

  for (
    program_now in PROGRAM_ORDER
  ) {
    requested <- HUMAN_PROGRAMS[[program_now]]

    actual <- vapply(
      requested,
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
    )

    keep_gene <- !is.na(
      actual
    )

    requested_present <- requested[
      keep_gene
    ]

    actual_present <- actual[
      keep_gene
    ]

    if (
      !length(
        actual_present
      )
    ) {
      next
    }

    zmat <- matrix(
      NA_real_,
      nrow = length(
        actual_present
      ),
      ncol = length(
        SAMPLE_ORDER
      ),
      dimnames = list(
        requested_present,
        SAMPLE_ORDER
      )
    )

    for (
      j in seq_along(
        actual_present
      )
    ) {
      zmat[j, ] <- safe_z(
        mat[
          actual_present[j],
          SAMPLE_ORDER
        ]
      )
    }

    score <- colMeans(
      zmat,
      na.rm = TRUE
    )

    for (
      i in seq_along(
        SAMPLE_ORDER
      )
    ) {
      program_sample_rows[[psi]] <-
        data.frame(
          final_state =
            st,
          analysis_tier =
            unname(
              STATE_TIER[[st]]
            ),
          program =
            program_now,
          detected_genes =
            length(
              actual_present
            ),
          requested_genes =
            length(
              requested
            ),
          sample =
            SAMPLE_ORDER[i],
          NAS =
            sample_info$NAS[i],
          fibrosis =
            sample_info$fibrosis[i],
          program_mean_gene_relative_z =
            score[
              SAMPLE_ORDER[i]
            ],
          stringsAsFactors = FALSE
        )

      psi <- psi + 1L
    }

    m0 <- safe_mean(
      score[
        c(
          "LT7",
          "LT8"
        )
      ]
    )

    m2 <- safe_mean(
      score[
        c(
          "LT5",
          "LT6"
        )
      ]
    )

    m4 <- safe_mean(
      score[
        c(
          "LT3",
          "LT4"
        )
      ]
    )

    m6 <- safe_mean(
      score[
        c(
          "LT1",
          "LT2"
        )
      ]
    )

    endpoint_delta <- m6 - m0

    advanced_concordant <- if (
      is.finite(
        endpoint_delta
      ) &&
        endpoint_delta >
          0
    ) {
      score[["LT1"]] >
        m0 &&
      score[["LT2"]] >
        m0
    } else if (
      is.finite(
        endpoint_delta
      ) &&
        endpoint_delta <
          0
    ) {
      score[["LT1"]] <
        m0 &&
      score[["LT2"]] <
        m0
    } else {
      FALSE
    }

    program_stage_rows[[psti]] <-
      data.frame(
        final_state =
          st,
        analysis_tier =
          unname(
            STATE_TIER[[st]]
          ),
        program =
          program_now,
        requested_genes =
          length(
            requested
          ),
        detected_genes =
          length(
            actual_present
          ),
        coverage_fraction =
          length(
            actual_present
          ) /
          length(
            requested
          ),
        mean_NAS0_F0 =
          m0,
        mean_NAS2_F0 =
          m2,
        mean_NAS4_mixed_F0_F1 =
          m4,
        mean_NAS6_F3_F4 =
          m6,
        delta_NAS6_vs_NAS0 =
          endpoint_delta,
        LT1_NAS6_F4_score =
          score[["LT1"]],
        LT2_NAS6_F3_score =
          score[["LT2"]],
        both_advanced_samples_concordant_vs_NAS0_mean =
          advanced_concordant,
        interpretation_boundary =
          "Descriptive cross-species ortholog program score; gene-wise z-score within state across 8 samples; advanced stage confounded with fibrosis",
        stringsAsFactors = FALSE
      )

    psti <- psti + 1L
  }
}

program_by_sample <- do.call(
  rbind,
  program_sample_rows
)

program_stage_summary <- do.call(
  rbind,
  program_stage_rows
)

write_tsv(
  program_by_sample,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_07_program_score_by_state_sample.tsv"
    )
  )
)

write_tsv(
  program_stage_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_08_program_stage_summary.tsv"
    )
  )
)

# ==============================================================================
# 8. Compact translational bridge summary
# ==============================================================================

bridge_rows <- list()
bri <- 1L

for (
  st in STATE_ORDER
) {
  m <- machinery_summary[
    machinery_summary$final_state ==
      st,
    ,
    drop = FALSE
  ]

  p <- program_stage_summary[
    program_stage_summary$final_state ==
      st,
    ,
    drop = FALSE
  ]

  get_program_delta <- function(
    name
  ) {
    x <- p[
      p$program ==
        name,
      "delta_NAS6_vs_NAS0",
      drop = TRUE
    ]

    if (
      length(
        x
      ) ==
        1L
    ) {
      return(
        x
      )
    }

    NA_real_
  }

  get_program_concordance <- function(
    name
  ) {
    x <- p[
      p$program ==
        name,
      "both_advanced_samples_concordant_vs_NAS0_mean",
      drop = TRUE
    ]

    if (
      length(
        x
      ) ==
        1L
    ) {
      return(
        x
      )
    }

    FALSE
  }

  bridge_rows[[bri]] <-
    data.frame(
      final_state =
        st,
      analysis_tier =
        unname(
          STATE_TIER[[st]]
        ),
      IL10_receptor_pair_present_at_state_level =
        m$receptor_pair_present_at_state_level,
      IL10_downstream_coverage_fraction =
        m$downstream_coverage_fraction,
      IL10RA_advanced_minus_baseline =
        m$IL10RA_NAS6_minus_NAS0,
      IL10RB_advanced_minus_baseline =
        m$IL10RB_NAS6_minus_NAS0,
      IL10_downstream_program_advanced_minus_baseline =
        get_program_delta(
          "IL10_downstream_no_receptors"
        ),
      Anti_inflammatory_advanced_minus_baseline =
        get_program_delta(
          "Anti_inflammatory"
        ),
      Repair_Resolution_advanced_minus_baseline =
        get_program_delta(
          "Repair_Resolution"
        ),
      Efferocytosis_advanced_minus_baseline =
        get_program_delta(
          "Efferocytosis"
        ),
      Inflammatory_advanced_minus_baseline =
        get_program_delta(
          "Inflammatory"
        ),
      ECM_inflammatory_advanced_minus_baseline =
        get_program_delta(
          "ECM_associated_inflammatory"
        ),
      Lipid_TREM2_advanced_minus_baseline =
        get_program_delta(
          "Lipid_TREM2"
        ),
      IL10_downstream_advanced_concordant =
        get_program_concordance(
          "IL10_downstream_no_receptors"
        ),
      Anti_inflammatory_advanced_concordant =
        get_program_concordance(
          "Anti_inflammatory"
        ),
      Inflammatory_advanced_concordant =
        get_program_concordance(
          "Inflammatory"
        ),
      ECM_inflammatory_advanced_concordant =
        get_program_concordance(
          "ECM_associated_inflammatory"
        ),
      translational_boundary =
        "Supports or refutes availability of the mouse-defined transcriptional architecture in human disease states; does not establish treatment response in human subjects",
      stringsAsFactors = FALSE
    )

  bri <- bri + 1L
}

bridge_summary <- do.call(
  rbind,
  bridge_rows
)

write_tsv(
  bridge_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_09_translational_bridge_summary.tsv"
    )
  )
)

# ==============================================================================
# 9. Provenance / interpretation boundaries
# ==============================================================================

provenance <- data.frame(
  key = c(
    "version",
    "human_logCPM_RDS",
    "human_trend_CSV",
    "human_stage_pattern_CSV",
    "mouse_program_definition_source",
    "human_state_count",
    "biological_sample_count",
    "advanced_disease_boundary",
    "IL10_receptor_boundary",
    "crossspecies_program_boundary",
    "statistics_boundary"
  ),
  value = c(
    VERSION,
    LOGCPM_RDS,
    TREND_CSV,
    STAGE_PATTERN_CSV,
    MOUSE_PROGRAM_SOURCE,
    as.character(
      length(
        STATE_ORDER
      )
    ),
    as.character(
      length(
        SAMPLE_ORDER
      )
    ),
    "LT1/LT2 are NAS6 and F4/F3; advanced disease severity is confounded with fibrosis and must not be called a pure NAS effect",
    "Pseudobulk presence of IL10RA and IL10RB at state level does not prove same-cell receptor coexpression or IL10 signaling activation",
    "Mouse v4.16.1 program symbols are mapped to candidate human ortholog symbols by uppercase one-to-one nomenclature; coverage is audited state-by-state",
    "Existing v2.2.1 edgeR gene-level results are retained; new program summaries are descriptive and sample-level, not cell-level inferential tests"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  provenance,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_10_provenance_and_boundaries.tsv"
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
# 10. Console summaries
# ==============================================================================

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "IL10 RESPONSE MACHINERY BY HUMAN MYELOID STATE"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  machinery_summary,
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
  "IL10 AXIS | PRIMARY HUMAN MYELOID STATES"
)

msg(
  strrep(
    "=",
    76
  )
)

axis_primary <- axis_stage_summary[
  axis_stage_summary$final_state %in%
    PRIMARY_STATES,
  c(
    "final_state",
    "gene",
    "mean_NAS0_F0",
    "mean_NAS6_F3_F4",
    "delta_NAS6_vs_NAS0",
    "logFC_per_severity_step_NAS2",
    "FDR",
    "pattern_class"
  ),
  drop = FALSE
]

print(
  axis_primary,
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
  "CROSS-SPECIES PROGRAM CONSERVATION | PRIMARY HUMAN MYELOID STATES"
)

msg(
  strrep(
    "=",
    76
  )
)

program_primary <- program_stage_summary[
  program_stage_summary$final_state %in%
    PRIMARY_STATES,
  c(
    "final_state",
    "program",
    "coverage_fraction",
    "mean_NAS0_F0",
    "mean_NAS6_F3_F4",
    "delta_NAS6_vs_NAS0",
    "both_advanced_samples_concordant_vs_NAS0_mean"
  ),
  drop = FALSE
]

print(
  program_primary,
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
  "TRANSLATIONAL BRIDGE SUMMARY"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  bridge_summary,
  row.names = FALSE
)

capture.output(
  print(
    machinery_summary,
    row.names = FALSE
  ),
  print(
    axis_primary,
    row.names = FALSE
  ),
  print(
    program_primary,
    row.names = FALSE
  ),
  print(
    bridge_summary,
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
  "PASS: VALIDATED HUMAN FINAL-MYELOID v2.2.1/v2.3 INPUTS LOADED"
)

msg(
  "PASS: EXACT MOUSE v4.16.1 PROGRAM DEFINITIONS RECOVERED AND HUMAN ORTHOLOG COVERAGE AUDITED"
)

msg(
  "PASS: HUMAN IL10 RECEPTOR/DOWNSTREAM MACHINERY AUDITED ACROSS FIVE FINAL MYELOID STATES"
)

msg(
  "PASS: HUMAN DISEASE-SEVERITY PROGRAM TRAJECTORIES SUMMARIZED AT BIOLOGICAL-SAMPLE LEVEL"
)

msg(
  "PASS: ",
  VERSION,
  " HUMAN MYELOID IL10 CONSERVATION AUDIT COMPLETE"
)

msg(
  strrep(
    "=",
    76
  )
)
