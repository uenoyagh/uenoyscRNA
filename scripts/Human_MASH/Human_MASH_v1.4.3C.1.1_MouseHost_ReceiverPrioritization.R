#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

VERSION <- "v1.4.3C.1.1"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

C0_DIR <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3C.0",
  "analysis",
  "v1.4.3C.0_human_ligand_mouse_host_receptor_preflight",
  "tables"
)

HUMAN_B2C3_DIR <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3B.2C.3",
  "analysis",
  "v1.4.3B.2C.3_prioritized_ligand_molecule_reconstruction",
  "tables"
)

IN_COMPLEX_SAMPLE <- file.path(
  C0_DIR,
  "v1.4.3C.0_11_receiver_complex_by_sample_broad.tsv"
)

IN_IL10_AXIS <- file.path(
  C0_DIR,
  "v1.4.3C.0_12_IL10_axis_by_sample_broad.tsv"
)

IN_BRIDGE <- file.path(
  C0_DIR,
  "v1.4.3C.0_04_human_ligand_mouse_receiver_bridge.tsv"
)

IN_COVERAGE <- file.path(
  C0_DIR,
  "v1.4.3C.0_05_feature_coverage.tsv"
)

IN_HUMAN <- file.path(
  HUMAN_B2C3_DIR,
  "v1.4.3B.2C.3_crosssample_strict_shared_TRUE_summary.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3C.1.1",
  "analysis",
  "v1.4.3C.1.1_human_ligand_mouse_receiver_prioritization"
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

read_tsv <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv <- function(x, path) {
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

required_files <- c(
  IN_COMPLEX_SAMPLE,
  IN_IL10_AXIS,
  IN_BRIDGE,
  IN_COVERAGE,
  IN_HUMAN
)

missing_files <-
  required_files[
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
    "Required inputs missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
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
  " | human graft ligand -> mouse host receiver prioritization"
)

msg(
  strrep(
    "=",
    76
  )
)

# ==============================================================================
# 1. Load validated v1.4.3C.0 / v1.4.3B.2C.3 tables
# ==============================================================================

complex_sample <-
  read_tsv(
    IN_COMPLEX_SAMPLE
  )

il10_axis <-
  read_tsv(
    IN_IL10_AXIS
  )

bridge <-
  read_tsv(
    IN_BRIDGE
  )

coverage <-
  read_tsv(
    IN_COVERAGE
  )

human <-
  read_tsv(
    IN_HUMAN
  )

required_complex_cols <- c(
  "sample",
  "condition",
  "broad",
  "receiver_complex",
  "n_cells",
  "positive_cells",
  "positive_fraction"
)

missing_complex_cols <-
  setdiff(
    required_complex_cols,
    colnames(
      complex_sample
    )
  )

if (
  length(
    missing_complex_cols
  )
) {
  stop2(
    "Complex input missing columns: ",
    paste(
      missing_complex_cols,
      collapse = ", "
    )
  )
}

expected_samples <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

st <-
  complex_sample[
    complex_sample$sample %in%
      expected_samples,
    ,
    drop = FALSE
  ]

if (
  !all(
    expected_samples %in%
      unique(
        st$sample
      )
  )
) {
  stop2(
    "Sham/Tx sample set is incomplete in receiver-complex table."
  )
}

# ==============================================================================
# 2. Human ligand -> receiver-complex map
# ==============================================================================

ligand_complex_map <-
  data.frame(
    human_ligand = c(
      "IL10",
      "FGF1",
      "TNFSF12",
      "VEGFA",
      "CXCL12",
      "TGFB3",
      "BMP5",
      "BMP2",
      "HBEGF",
      "BTC"
    ),
    receiver_complex = c(
      "IL10R_Il10ra_AND_Il10rb",
      "FGFR_ANY_Fgfr1_4",
      "TWEAKR_Tnfrsf12a",
      "VEGFR_core_Flt1_OR_Kdr",
      "CXCL12R_Cxcr4_OR_Ackr3",
      "TGFBR_Tgfbr1_AND_Tgfbr2",
      "BMPR_Bmpr2_AND_typeI",
      "BMPR_Bmpr2_AND_typeI",
      "EGFR_Egfr_OR_Erbb4",
      "EGFR_Egfr_OR_Erbb4"
    ),
    stringsAsFactors = FALSE
  )

missing_complexes <-
  setdiff(
    unique(
      ligand_complex_map$receiver_complex
    ),
    unique(
      st$receiver_complex
    )
  )

if (
  length(
    missing_complexes
  )
) {
  stop2(
    "Expected receiver complexes absent: ",
    paste(
      missing_complexes,
      collapse = ", "
    )
  )
}

# ==============================================================================
# 3. Sample-level Sham vs Tx receptor competence by broad lineage
# ==============================================================================

safe_get <- function(
  z,
  sample,
  col
) {
  x <-
    z[
      z$sample == sample,
      col,
      drop = TRUE
    ]

  if (
    length(
      x
    ) == 1L
  ) {
    return(
      x
    )
  }

  NA_real_
}

complex_key <-
  interaction(
    st$receiver_complex,
    st$broad,
    drop = TRUE,
    sep = "|||"
  )

complex_split <-
  split(
    seq_len(
      nrow(
        st
      )
    ),
    complex_key
  )

complex_summary_rows <- list()

ri <- 1L

for (
  idx in complex_split
) {
  z <-
    st[
      idx,
      ,
      drop = FALSE
    ]

  receptor <-
    z$receiver_complex[
      1
    ]

  broad <-
    z$broad[
      1
    ]

  sham_pf <- c(
    safe_get(
      z,
      "Sham1",
      "positive_fraction"
    ),
    safe_get(
      z,
      "Sham20",
      "positive_fraction"
    )
  )

  tx_pf <- c(
    safe_get(
      z,
      "Tx17",
      "positive_fraction"
    ),
    safe_get(
      z,
      "Tx5",
      "positive_fraction"
    )
  )

  sham_n <- c(
    safe_get(
      z,
      "Sham1",
      "n_cells"
    ),
    safe_get(
      z,
      "Sham20",
      "n_cells"
    )
  )

  tx_n <- c(
    safe_get(
      z,
      "Tx17",
      "n_cells"
    ),
    safe_get(
      z,
      "Tx5",
      "n_cells"
    )
  )

  complete <-
    all(
      is.finite(
        c(
          sham_pf,
          tx_pf,
          sham_n,
          tx_n
        )
      )
    )

  complex_summary_rows[[ri]] <-
    data.frame(
      receiver_complex = receptor,
      broad = broad,
      Sham1_positive_fraction =
        sham_pf[
          1
        ],
      Sham20_positive_fraction =
        sham_pf[
          2
        ],
      Tx17_positive_fraction =
        tx_pf[
          1
        ],
      Tx5_positive_fraction =
        tx_pf[
          2
        ],
      Sham_mean_positive_fraction =
        if (
          all(
            is.finite(
              sham_pf
            )
          )
        ) {
          mean(
            sham_pf
          )
        } else {
          NA_real_
        },
      Tx_mean_positive_fraction =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          mean(
            tx_pf
          )
        } else {
          NA_real_
        },
      Tx_min_positive_fraction =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          min(
            tx_pf
          )
        } else {
          NA_real_
        },
      Tx_max_positive_fraction =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          max(
            tx_pf
          )
        } else {
          NA_real_
        },
      Tx_replicate_spread =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          abs(
            diff(
              tx_pf
            )
          )
        } else {
          NA_real_
        },
      Tx_minus_Sham_positive_fraction =
        if (
          complete
        ) {
          mean(
            tx_pf
          ) -
            mean(
              sham_pf
            )
        } else {
          NA_real_
        },
      replicate_complete =
        complete,
      inference =
        "descriptive_only_n2_per_condition",
      stringsAsFactors = FALSE
    )

  ri <- ri + 1L
}

complex_summary <-
  do.call(
    rbind,
    complex_summary_rows
  )

complex_summary <-
  complex_summary[
    order(
      complex_summary$receiver_complex,
      -complex_summary$Tx_mean_positive_fraction,
      -complex_summary$Tx_min_positive_fraction,
      complex_summary$broad
    ),
    ,
    drop = FALSE
  ]

complex_summary$rank_within_receiver_Tx <-
  ave(
    -complex_summary$Tx_mean_positive_fraction,
    complex_summary$receiver_complex,
    FUN = function(x) {
      rank(
        x,
        ties.method = "min",
        na.last = "keep"
      )
    }
  )

write_tsv(
  complex_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_01_receiver_complex_Sham_Tx_by_broad.tsv"
    )
  )
)

# ==============================================================================
# 4. Add human graft ligand evidence
# ==============================================================================

human_keep <- c(
  "gene",
  "SfLB1_MEX_UMI",
  "SfLB2_MEX_UMI",
  "Tx17_MEX_UMI",
  "Tx5_MEX_UMI",
  "SfLB1_shared_TRUE_UMI",
  "SfLB2_shared_TRUE_UMI",
  "Tx17_shared_TRUE_UMI",
  "Tx5_shared_TRUE_UMI",
  "strict_shared_TRUE_Tx_vs_SfLB_log2FC_per_cell_descriptive"
)

missing_human_cols <-
  setdiff(
    human_keep,
    colnames(
      human
    )
  )

if (
  length(
    missing_human_cols
  )
) {
  stop2(
    "Human ligand table missing columns: ",
    paste(
      missing_human_cols,
      collapse = ", "
    )
  )
}

human2 <-
  human[
    ,
    human_keep,
    drop = FALSE
  ]

names(
  human2
)[
  names(
    human2
  ) == "gene"
] <- "human_ligand"

ligand_receiver <-
  merge(
    ligand_complex_map,
    complex_summary,
    by = "receiver_complex",
    all.x = TRUE,
    sort = FALSE
  )

ligand_receiver <-
  merge(
    ligand_receiver,
    human2,
    by = "human_ligand",
    all.x = TRUE,
    sort = FALSE
  )

ligand_receiver <-
  ligand_receiver[
    order(
      match(
        ligand_receiver$human_ligand,
        ligand_complex_map$human_ligand
      ),
      -ligand_receiver$Tx_mean_positive_fraction,
      -ligand_receiver$Tx_min_positive_fraction,
      ligand_receiver$broad
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  ligand_receiver,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_02_human_ligand_mouse_receiver_bridge_by_broad.tsv"
    )
  )
)

# ==============================================================================
# 5. IL10-focused host receiver ranking
# ==============================================================================

il10_receiver <-
  complex_summary[
    complex_summary$receiver_complex ==
      "IL10R_Il10ra_AND_Il10rb",
    ,
    drop = FALSE
  ]

il10_receiver <-
  il10_receiver[
    order(
      -il10_receiver$Tx_mean_positive_fraction,
      -il10_receiver$Tx_min_positive_fraction,
      il10_receiver$Tx_replicate_spread,
      il10_receiver$broad
    ),
    ,
    drop = FALSE
  ]

il10_receiver$Tx_receiver_rank <-
  seq_len(
    nrow(
      il10_receiver
    )
  )

core_lineages <- c(
  "Kupffer/Macrophage",
  "Monocyte",
  "HSC/Mesenchymal",
  "Hepatocyte",
  "LSEC",
  "Cholangiocyte",
  "Vascular endothelial",
  "Dendritic cell",
  "Neutrophil"
)

il10_receiver$core_mechanistic_lineage <-
  il10_receiver$broad %in%
    core_lineages

il10_receiver$roadmap_role <-
  ifelse(
    il10_receiver$broad %in%
      c(
        "Kupffer/Macrophage",
        "Monocyte"
      ),
    "primary_myeloid_receiver_candidate",
    ifelse(
      il10_receiver$broad %in%
        c(
          "HSC/Mesenchymal",
          "Hepatocyte",
          "LSEC",
          "Cholangiocyte"
        ),
      "downstream_tissue_lineage",
      ifelse(
        il10_receiver$broad ==
          "Neutrophil",
        "secondary_high_IL10R_not_current_core_roadmap",
        "context_or_secondary"
      )
    )
  )

write_tsv(
  il10_receiver,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_03_IL10R_host_receiver_ranked.tsv"
    )
  )
)

# ==============================================================================
# 6. IL10 downstream pathway support at lineage level
# ==============================================================================

required_il10_cols <- c(
  "sample",
  "condition",
  "broad",
  "gene",
  "n_cells",
  "UMI",
  "positive_cells",
  "positive_fraction",
  "UMI_per_1000_cells"
)

missing_il10_cols <-
  setdiff(
    required_il10_cols,
    colnames(
      il10_axis
    )
  )

if (
  length(
    missing_il10_cols
  )
) {
  stop2(
    "IL10 axis table missing columns: ",
    paste(
      missing_il10_cols,
      collapse = ", "
    )
  )
}

il10_axis_st <-
  il10_axis[
    il10_axis$sample %in%
      expected_samples,
    ,
    drop = FALSE
  ]

axis_key <-
  interaction(
    il10_axis_st$gene,
    il10_axis_st$broad,
    drop = TRUE,
    sep = "|||"
  )

axis_split <-
  split(
    seq_len(
      nrow(
        il10_axis_st
      )
    ),
    axis_key
  )

axis_rows <- list()

ai <- 1L

for (
  idx in axis_split
) {
  z <-
    il10_axis_st[
      idx,
      ,
      drop = FALSE
    ]

  gene <-
    z$gene[
      1
    ]

  broad <-
    z$broad[
      1
    ]

  sham_pf <- c(
    safe_get(
      z,
      "Sham1",
      "positive_fraction"
    ),
    safe_get(
      z,
      "Sham20",
      "positive_fraction"
    )
  )

  tx_pf <- c(
    safe_get(
      z,
      "Tx17",
      "positive_fraction"
    ),
    safe_get(
      z,
      "Tx5",
      "positive_fraction"
    )
  )

  sham_u <- c(
    safe_get(
      z,
      "Sham1",
      "UMI_per_1000_cells"
    ),
    safe_get(
      z,
      "Sham20",
      "UMI_per_1000_cells"
    )
  )

  tx_u <- c(
    safe_get(
      z,
      "Tx17",
      "UMI_per_1000_cells"
    ),
    safe_get(
      z,
      "Tx5",
      "UMI_per_1000_cells"
    )
  )

  axis_rows[[ai]] <-
    data.frame(
      gene = gene,
      broad = broad,
      Sham_mean_positive_fraction =
        if (
          all(
            is.finite(
              sham_pf
            )
          )
        ) {
          mean(
            sham_pf
          )
        } else {
          NA_real_
        },
      Tx_mean_positive_fraction =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          mean(
            tx_pf
          )
        } else {
          NA_real_
        },
      Tx_min_positive_fraction =
        if (
          all(
            is.finite(
              tx_pf
            )
          )
        ) {
          min(
            tx_pf
          )
        } else {
          NA_real_
        },
      Sham_mean_UMI_per_1000_cells =
        if (
          all(
            is.finite(
              sham_u
            )
          )
        ) {
          mean(
            sham_u
          )
        } else {
          NA_real_
        },
      Tx_mean_UMI_per_1000_cells =
        if (
          all(
            is.finite(
              tx_u
            )
          )
        ) {
          mean(
            tx_u
          )
        } else {
          NA_real_
        },
      inference =
        "supporting_expression_only_not_pathway_activation",
      stringsAsFactors = FALSE
    )

  ai <- ai + 1L
}

axis_summary <-
  do.call(
    rbind,
    axis_rows
  )

axis_summary <-
  axis_summary[
    order(
      match(
        axis_summary$gene,
        c(
          "Il10ra",
          "Il10rb",
          "Jak1",
          "Tyk2",
          "Stat3",
          "Socs3"
        )
      ),
      -axis_summary$Tx_mean_positive_fraction,
      axis_summary$broad
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  axis_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_04_IL10_axis_support_by_broad.tsv"
    )
  )
)

# ==============================================================================
# 7. Core-lineage bridge matrix for manuscript/mechanistic planning
# ==============================================================================

core_bridge <-
  ligand_receiver[
    ligand_receiver$broad %in%
      core_lineages,
    ,
    drop = FALSE
  ]

core_bridge$core_lineage_role <-
  ifelse(
    core_bridge$broad %in%
      c(
        "Kupffer/Macrophage",
        "Monocyte"
      ),
    "myeloid",
    ifelse(
      core_bridge$broad ==
        "HSC/Mesenchymal",
      "fibrosis_effector",
      ifelse(
        core_bridge$broad %in%
          c(
            "LSEC",
            "Vascular endothelial"
          ),
        "endothelial",
        ifelse(
          core_bridge$broad ==
            "Hepatocyte",
          "parenchymal",
          ifelse(
            core_bridge$broad ==
              "Cholangiocyte",
            "ductal",
            "secondary_immune"
          )
        )
      )
    )
  )

core_bridge <-
  core_bridge[
    order(
      match(
        core_bridge$human_ligand,
        ligand_complex_map$human_ligand
      ),
      -core_bridge$Tx_mean_positive_fraction,
      -core_bridge$Tx_min_positive_fraction,
      core_bridge$broad
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  core_bridge,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_05_core_lineage_bridge_matrix.tsv"
    )
  )
)

# ==============================================================================
# 8. Coverage caveat
# ==============================================================================

bmpr1b_absent <-
  "Bmpr1b" %in%
    coverage$requested_gene[
      !coverage$present
    ]

coverage_note <-
  data.frame(
    item = c(
      "Bmpr1b_absent",
      "BMPR_logic",
      "IL10_activation_boundary",
      "statistics_boundary"
    ),
    value = c(
      as.character(
        bmpr1b_absent
      ),
      "BMPR competence remains computable with Bmpr2 AND ANY(Bmpr1a,Acvr1) among measured type-I receptors; Bmpr1b is unavailable in this RNA feature space",
      "Il10ra/Il10rb and Jak1/Tyk2/Stat3/Socs3 transcript availability supports response competence only; scRNA-seq does not measure pSTAT3 activation",
      "Sham1/Sham20 versus Tx17/Tx5 summaries are descriptive only (n=2 biological samples per condition)"
    ),
    stringsAsFactors = FALSE
  )

write_tsv(
  coverage_note,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_06_interpretation_boundaries.tsv"
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
# 9. Console summary
# ==============================================================================

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "IL10R DIRECT HOST RECEIVER RANKING | Tx mean dual-positive fraction"
)

msg(
  strrep(
    "=",
    76
  )
)

print_cols <- c(
  "Tx_receiver_rank",
  "broad",
  "Sham_mean_positive_fraction",
  "Tx17_positive_fraction",
  "Tx5_positive_fraction",
  "Tx_mean_positive_fraction",
  "Tx_min_positive_fraction",
  "Tx_minus_Sham_positive_fraction",
  "roadmap_role"
)

print(
  il10_receiver[
    ,
    print_cols,
    drop = FALSE
  ],
  row.names = FALSE
)

capture.output(
  print(
    il10_receiver[
      ,
      print_cols,
      drop = FALSE
    ],
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
  "TOP 3 Tx RECEIVER LINEAGES PER DIRECT HUMAN LIGAND"
)

msg(
  strrep(
    "=",
    76
  )
)

top_rows <- list()

ti <- 1L

for (
  ligand in ligand_complex_map$human_ligand
) {
  z <-
    ligand_receiver[
      ligand_receiver$human_ligand ==
        ligand,
      ,
      drop = FALSE
    ]

  z <-
    z[
      order(
        -z$Tx_mean_positive_fraction,
        -z$Tx_min_positive_fraction,
        z$broad
      ),
      ,
      drop = FALSE
    ]

  z <-
    head(
      z,
      3L
    )

  if (
    nrow(
      z
    )
  ) {
    z$rank <-
      seq_len(
        nrow(
          z
        )
      )

    top_rows[[ti]] <-
      z[
        ,
        c(
          "human_ligand",
          "receiver_complex",
          "rank",
          "broad",
          "Tx_mean_positive_fraction",
          "Tx_min_positive_fraction",
          "Tx_minus_Sham_positive_fraction",
          "strict_shared_TRUE_Tx_vs_SfLB_log2FC_per_cell_descriptive"
        ),
        drop = FALSE
      ]

    ti <- ti + 1L
  }
}

top3 <-
  do.call(
    rbind,
    top_rows
  )

print(
  top3,
  row.names = FALSE
)

capture.output(
  print(
    top3,
    row.names = FALSE
  ),
  file = LOG_FILE,
  append = TRUE
)

write_tsv(
  top3,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_07_top3_receivers_per_ligand.tsv"
    )
  )
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
  "PASS: HUMAN PROBE-VALIDATED LIGAND EVIDENCE LINKED TO MOUSE RECEIVER COMPLEXES"
)

msg(
  "PASS: SHAM/TX SAMPLE-LEVEL RECEIVER COMPETENCE SUMMARIZED"
)

msg(
  "PASS: IL10R HOST RECEIVER LINEAGES RANKED"
)

msg(
  "PASS: ",
  VERSION,
  " HUMAN-LIGAND -> MOUSE-HOST RECEIVER PRIORITIZATION COMPLETE"
)

msg(
  strrep(
    "=",
    76
  )
)
