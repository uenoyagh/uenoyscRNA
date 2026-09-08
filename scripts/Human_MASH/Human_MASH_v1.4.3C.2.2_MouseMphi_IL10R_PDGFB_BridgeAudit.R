#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

VERSION <- "v1.4.3C.2.2"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

C21_IL10R_TABLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3C.2.1",
  "analysis",
  "v1.4.3C.2.1_mouse_Mphi_IL10_responder_mechanistic_audit",
  "tables",
  "v1.4.3C.2.1_04_IL10R_competence_by_sample_class.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3C.2.2",
  "analysis",
  "v1.4.3C.2.2_mouse_Mphi_IL10R_PDGFB_bridge_audit"
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

first_existing <- function(
  x,
  candidates
) {
  y <-
    candidates[
      candidates %in%
        x
    ]

  if (
    length(
      y
    )
  ) {
    return(
      y[
        1
      ]
    )
  }

  NA_character_
}

canonical_condition <- function(
  sample
) {
  x <-
    as.character(
      sample
    )

  out <-
    rep(
      NA_character_,
      length(
        x
      )
    )

  out[
    x %in%
      c(
        "Sham1",
        "Sham20"
      )
  ] <- "Sham"

  out[
    x %in%
      c(
        "Tx17",
        "Tx5"
      )
  ] <- "Tx"

  out
}

safe_mean <- function(
  x
) {
  x <-
    x[
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
  x <-
    x[
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
  x <-
    x[
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
  " | focused mouse Mphi IL10R-PDGFB bridge audit"
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
    "Input RDS not found: ",
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
# 1. Load final Clean-B macrophage object
# ==============================================================================

msg(
  "Loading FINAL Clean-B macrophage RDS..."
)

mphi <-
  readRDS(
    INPUT_RDS
  )

if (
  !inherits(
    mphi,
    "Seurat"
  )
) {
  stop2(
    "Input object is not a Seurat object."
  )
}

if (
  !"RNA" %in%
    Assays(
      mphi
    )
) {
  stop2(
    "RNA assay is absent."
  )
}

DefaultAssay(
  mphi
) <- "RNA"

mphi <-
  JoinLayers(
    mphi,
    assay = "RNA"
  )

rna_layers <-
  Layers(
    mphi[["RNA"]]
  )

if (
  "counts" %in%
    rna_layers
) {
  counts_layer <- "counts"
} else {
  count_candidates <-
    grep(
      "^counts($|\\.)",
      rna_layers,
      value = TRUE
    )

  if (
    length(
      count_candidates
    ) != 1L
  ) {
    stop2(
      "Could not resolve a unique RNA counts layer. Layers: ",
      paste(
        rna_layers,
        collapse = ", "
      )
    )
  }

  counts_layer <-
    count_candidates[
      1
    ]
}

counts <-
  LayerData(
    mphi,
    assay = "RNA",
    layer = counts_layer
  )

if (
  !"Pdgfb" %in%
    rownames(
      counts
    )
) {
  stop2(
    "Pdgfb is absent from RNA counts."
  )
}

if (
  !all(
    c(
      "Il10ra",
      "Il10rb"
    ) %in%
      rownames(
        counts
      )
  )
) {
  stop2(
    "Il10ra and/or Il10rb are absent from RNA counts."
  )
}

meta_cols <-
  colnames(
    mphi@meta.data
  )

SAMPLE_COL <-
  first_existing(
    meta_cols,
    c(
      "sample_for_annotation",
      "sample_FIXED2",
      "sample_4group",
      "sample",
      "sample_id",
      "orig.ident"
    )
  )

CLASS_COL <-
  first_existing(
    meta_cols,
    c(
      "macrophage_class_Res2_FINAL_v4145_char",
      "macrophage_class_Res2_FINAL_v4145"
    )
  )

if (
  is.na(
    SAMPLE_COL
  )
) {
  stop2(
    "Sample metadata column not found."
  )
}

if (
  is.na(
    CLASS_COL
  )
) {
  stop2(
    "Final macrophage class column not found."
  )
}

SAMPLE_ORDER <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

CLASS_ORDER <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

PRIMARY_CLASS <-
  "Repair/Resolution-Mphi"

meta <-
  data.frame(
    cell =
      colnames(
        mphi
      ),
    sample =
      as.character(
        mphi@meta.data[[SAMPLE_COL]]
      ),
    macrophage_class =
      as.character(
        mphi@meta.data[[CLASS_COL]]
      ),
    stringsAsFactors = FALSE
  )

meta$condition <-
  canonical_condition(
    meta$sample
  )

keep <-
  meta$sample %in%
    SAMPLE_ORDER &
  meta$macrophage_class %in%
    CLASS_ORDER

meta <-
  meta[
    keep,
    ,
    drop = FALSE
  ]

counts <-
  counts[
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

if (
  anyNA(
    meta$condition
  )
) {
  stop2(
    "Unresolved Sham/Tx condition labels."
  )
}

msg(
  "Sample column: ",
  SAMPLE_COL
)

msg(
  "Final class column: ",
  CLASS_COL
)

msg(
  "Retained Sham/Tx final-class Mphi cells: ",
  nrow(
    meta
  )
)

# ==============================================================================
# 2. Cell-level exact metrics
# ==============================================================================

cell_library_UMI <-
  Matrix::colSums(
    counts
  )

pdgfb_count <-
  as.numeric(
    counts[
      "Pdgfb",
      ,
      drop = TRUE
    ]
  )

il10ra_count <-
  as.numeric(
    counts[
      "Il10ra",
      ,
      drop = TRUE
    ]
  )

il10rb_count <-
  as.numeric(
    counts[
      "Il10rb",
      ,
      drop = TRUE
    ]
  )

meta$library_UMI <-
  as.numeric(
    cell_library_UMI
  )

meta$Pdgfb_UMI <-
  pdgfb_count

meta$Pdgfb_positive <-
  pdgfb_count >
  0

meta$Pdgfb_cell_CP10k <-
  ifelse(
    meta$library_UMI >
      0,
    pdgfb_count /
      meta$library_UMI *
      1e4,
    NA_real_
  )

meta$Il10ra_positive <-
  il10ra_count >
  0

meta$Il10rb_positive <-
  il10rb_count >
  0

meta$IL10R_dual_positive <-
  meta$Il10ra_positive &
  meta$Il10rb_positive

meta$IL10R_status <-
  ifelse(
    meta$IL10R_dual_positive,
    "IL10R_dual_positive",
    "IL10R_not_dual_positive"
  )

# ==============================================================================
# 3. Reproduce v1.4.3C.2.1 IL10R fractions before focused analysis
# ==============================================================================

summarize_il10r_fraction <- function(
  meta_df
) {
  key <-
    interaction(
      meta_df$sample,
      meta_df$macrophage_class,
      drop = TRUE,
      sep = "|||"
    )

  idx_list <-
    split(
      seq_len(
        nrow(
          meta_df
        )
      ),
      key
    )

  rows <- list()
  ri <- 1L

  for (
    idx in idx_list
  ) {
    z <-
      meta_df[
        idx,
        ,
        drop = FALSE
      ]

    rows[[ri]] <-
      data.frame(
        sample =
          z$sample[
            1
          ],
        macrophage_class =
          z$macrophage_class[
            1
          ],
        n_cells =
          nrow(
            z
          ),
        IL10R_dual_positive_fraction =
          mean(
            z$IL10R_dual_positive
          ),
        stringsAsFactors = FALSE
      )

    ri <- ri + 1L
  }

  do.call(
    rbind,
    rows
  )
}

current_il10r <-
  summarize_il10r_fraction(
    meta
  )

c21_validation <-
  data.frame(
    status =
      "C2.1 table unavailable; validation skipped",
    max_abs_fraction_difference =
      NA_real_,
    stringsAsFactors = FALSE
  )

if (
  file.exists(
    C21_IL10R_TABLE
  )
) {
  old <-
    read.delim(
      C21_IL10R_TABLE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  req_old <- c(
    "sample",
    "macrophage_class",
    "IL10R_dual_positive_fraction"
  )

  if (
    all(
      req_old %in%
        colnames(
          old
        )
    )
  ) {
    chk <-
      merge(
        current_il10r,
        old[
          ,
          req_old,
          drop = FALSE
        ],
        by = c(
          "sample",
          "macrophage_class"
        ),
        suffixes = c(
          "_C22",
          "_C21"
        ),
        all = TRUE
      )

    chk$abs_difference <-
      abs(
        chk$IL10R_dual_positive_fraction_C22 -
        chk$IL10R_dual_positive_fraction_C21
      )

    max_diff <-
      max(
        chk$abs_difference,
        na.rm = TRUE
      )

    if (
      !is.finite(
        max_diff
      )
    ) {
      max_diff <- NA_real_
    }

    c21_validation <-
      data.frame(
        status =
          if (
            is.finite(
              max_diff
            ) &&
              max_diff <
                1e-12
          ) {
            "PASS_exact"
          } else {
            "WARNING_nonzero_difference"
          },
        max_abs_fraction_difference =
          max_diff,
        stringsAsFactors = FALSE
      )

    write_tsv(
      chk,
      file.path(
        TAB_OUT,
        paste0(
          VERSION,
          "_01_C21_IL10R_fraction_validation.tsv"
        )
      )
    )
  }
}

write_tsv(
  c21_validation,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_02_C21_validation_summary.tsv"
    )
  )
)

# ==============================================================================
# 4. Pdgfb summary helper
# ==============================================================================

summarize_pdgfb <- function(
  meta_df,
  group_cols
) {
  group_key <-
    do.call(
      interaction,
      c(
        as.list(
          meta_df[
            ,
            group_cols,
            drop = FALSE
          ]
        ),
        list(
          drop = TRUE,
          sep = "|||",
          lex.order = TRUE
        )
      )
    )

  idx_list <-
    split(
      seq_len(
        nrow(
          meta_df
        )
      ),
      group_key
    )

  rows <- list()
  ri <- 1L

  for (
    idx in idx_list
  ) {
    z <-
      meta_df[
        idx,
        ,
        drop = FALSE
      ]

    first <-
      z[
        1,
        group_cols,
        drop = FALSE
      ]

    n_cells <-
      nrow(
        z
      )

    total_lib <-
      sum(
        z$library_UMI
      )

    total_pdgfb <-
      sum(
        z$Pdgfb_UMI
      )

    tmp <-
      data.frame(
        n_cells =
          n_cells,
        total_library_UMI =
          total_lib,
        Pdgfb_total_UMI =
          total_pdgfb,
        Pdgfb_positive_cells =
          sum(
            z$Pdgfb_positive
          ),
        Pdgfb_positive_fraction =
          mean(
            z$Pdgfb_positive
          ),
        Pdgfb_mean_raw_UMI_per_cell =
          mean(
            z$Pdgfb_UMI
          ),
        Pdgfb_mean_per_cell_CP10k =
          safe_mean(
            z$Pdgfb_cell_CP10k
          ),
        Pdgfb_median_per_cell_CP10k =
          stats::median(
            z$Pdgfb_cell_CP10k,
            na.rm = TRUE
          ),
        Pdgfb_pseudobulk_CP10k =
          if (
            total_lib >
              0
          ) {
            total_pdgfb /
              total_lib *
              1e4
          } else {
            NA_real_
          },
        min20_cells =
          n_cells >=
          20L,
        stringsAsFactors = FALSE
      )

    rows[[ri]] <-
      cbind(
        first,
        tmp
      )

    ri <- ri + 1L
  }

  do.call(
    rbind,
    rows
  )
}

# ==============================================================================
# 5. Primary unsplit Repair/Resolution-Mphi benchmark
# ==============================================================================

rr <-
  meta[
    meta$macrophage_class ==
      PRIMARY_CLASS,
    ,
    drop = FALSE
  ]

rr_unsplit <-
  summarize_pdgfb(
    rr,
    c(
      "sample",
      "condition",
      "macrophage_class"
    )
  )

rr_unsplit <-
  rr_unsplit[
    order(
      match(
        rr_unsplit$sample,
        SAMPLE_ORDER
      )
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  rr_unsplit,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_03_RepairResolution_unsplit_Pdgfb_by_sample.tsv"
    )
  )
)

sham_rr <-
  rr_unsplit$Pdgfb_mean_per_cell_CP10k[
    rr_unsplit$condition ==
      "Sham"
  ]

tx_rr <-
  rr_unsplit$Pdgfb_mean_per_cell_CP10k[
    rr_unsplit$condition ==
      "Tx"
  ]

observed_sham_mean <-
  safe_mean(
    sham_rr
  )

observed_tx_mean <-
  safe_mean(
    tx_rr
  )

# v6.3.6 checkpoint benchmark:
# Repair/Resolution-Mphi per-cell Pdgfb CP10k
# Sham mean 5.971, Tx mean 3.568.
benchmark <-
  data.frame(
    condition = c(
      "Sham",
      "Tx"
    ),
    checkpoint_v636_mean_CP10k = c(
      5.971,
      3.568
    ),
    C22_observed_mean_CP10k = c(
      observed_sham_mean,
      observed_tx_mean
    ),
    stringsAsFactors = FALSE
  )

benchmark$absolute_difference <-
  abs(
    benchmark$C22_observed_mean_CP10k -
    benchmark$checkpoint_v636_mean_CP10k
  )

benchmark$within_0.50_CP10k <-
  benchmark$absolute_difference <=
  0.50

write_tsv(
  benchmark,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_04_v636_Pdgfb_CP10k_benchmark.tsv"
    )
  )
)

# ==============================================================================
# 6. Primary Repair/Resolution-Mphi: IL10R dual-positive vs not-dual-positive
# ==============================================================================

rr_status <-
  summarize_pdgfb(
    rr,
    c(
      "sample",
      "condition",
      "macrophage_class",
      "IL10R_status"
    )
  )

rr_status <-
  rr_status[
    order(
      match(
        rr_status$sample,
        SAMPLE_ORDER
      ),
      match(
        rr_status$IL10R_status,
        c(
          "IL10R_dual_positive",
          "IL10R_not_dual_positive"
        )
      )
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  rr_status,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_05_RepairResolution_Pdgfb_by_IL10R_status.tsv"
    )
  )
)

get_status_value <- function(
  z,
  status,
  col
) {
  x <-
    z[
      z$IL10R_status ==
        status,
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

sample_contrast_rows <- list()
sci <- 1L

for (
  sample_now in SAMPLE_ORDER
) {
  z <-
    rr_status[
      rr_status$sample ==
        sample_now,
    ,
      drop = FALSE
    ]

  dual_n <-
    get_status_value(
      z,
      "IL10R_dual_positive",
      "n_cells"
    )

  notdual_n <-
    get_status_value(
      z,
      "IL10R_not_dual_positive",
      "n_cells"
    )

  dual_cp10k <-
    get_status_value(
      z,
      "IL10R_dual_positive",
      "Pdgfb_mean_per_cell_CP10k"
    )

  notdual_cp10k <-
    get_status_value(
      z,
      "IL10R_not_dual_positive",
      "Pdgfb_mean_per_cell_CP10k"
    )

  dual_pb <-
    get_status_value(
      z,
      "IL10R_dual_positive",
      "Pdgfb_pseudobulk_CP10k"
    )

  notdual_pb <-
    get_status_value(
      z,
      "IL10R_not_dual_positive",
      "Pdgfb_pseudobulk_CP10k"
    )

  dual_pf <-
    get_status_value(
      z,
      "IL10R_dual_positive",
      "Pdgfb_positive_fraction"
    )

  notdual_pf <-
    get_status_value(
      z,
      "IL10R_not_dual_positive",
      "Pdgfb_positive_fraction"
    )

  sample_contrast_rows[[sci]] <-
    data.frame(
      sample =
        sample_now,
      condition =
        canonical_condition(
          sample_now
        ),
      IL10R_dual_positive_n_cells =
        dual_n,
      IL10R_not_dual_positive_n_cells =
        notdual_n,
      both_groups_min20 =
        is.finite(
          dual_n
        ) &&
        is.finite(
          notdual_n
        ) &&
        dual_n >=
          20 &&
        notdual_n >=
          20,
      dual_positive_mean_per_cell_CP10k =
        dual_cp10k,
      not_dual_mean_per_cell_CP10k =
        notdual_cp10k,
      dual_minus_notdual_mean_per_cell_CP10k =
        dual_cp10k -
        notdual_cp10k,
      dual_positive_pseudobulk_CP10k =
        dual_pb,
      not_dual_pseudobulk_CP10k =
        notdual_pb,
      dual_minus_notdual_pseudobulk_CP10k =
        dual_pb -
        notdual_pb,
      dual_positive_Pdgfb_positive_fraction =
        dual_pf,
      not_dual_Pdgfb_positive_fraction =
        notdual_pf,
      dual_minus_notdual_positive_fraction =
        dual_pf -
        notdual_pf,
      inference =
        "within_sample_descriptive_association_not_causal_IL10_effect",
      stringsAsFactors = FALSE
    )

  sci <- sci + 1L
}

rr_contrast <-
  do.call(
    rbind,
    sample_contrast_rows
  )

write_tsv(
  rr_contrast,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_06_RepairResolution_IL10R_Pdgfb_sample_contrast.tsv"
    )
  )
)

# ==============================================================================
# 7. Condition-level contrast summary and difference-of-differences
# ==============================================================================

sham_contrast <-
  rr_contrast$dual_minus_notdual_mean_per_cell_CP10k[
    rr_contrast$condition ==
      "Sham"
  ]

tx_contrast <-
  rr_contrast$dual_minus_notdual_mean_per_cell_CP10k[
    rr_contrast$condition ==
      "Tx"
  ]

sham_pb_contrast <-
  rr_contrast$dual_minus_notdual_pseudobulk_CP10k[
    rr_contrast$condition ==
      "Sham"
  ]

tx_pb_contrast <-
  rr_contrast$dual_minus_notdual_pseudobulk_CP10k[
    rr_contrast$condition ==
      "Tx"
  ]

sham_pf_contrast <-
  rr_contrast$dual_minus_notdual_positive_fraction[
    rr_contrast$condition ==
      "Sham"
  ]

tx_pf_contrast <-
  rr_contrast$dual_minus_notdual_positive_fraction[
    rr_contrast$condition ==
      "Tx"
  ]

tx_same_direction_cp10k <-
  all(
    is.finite(
      tx_contrast
    )
  ) &&
  sign(
    tx_contrast[
      1
    ]
  ) ==
  sign(
    tx_contrast[
      2
    ]
  )

tx_both_lower_cp10k <-
  all(
    is.finite(
      tx_contrast
    )
  ) &&
  all(
    tx_contrast <
      0
  )

tx_both_lower_pseudobulk <-
  all(
    is.finite(
      tx_pb_contrast
    )
  ) &&
  all(
    tx_pb_contrast <
      0
  )

tx_both_lower_positive_fraction <-
  all(
    is.finite(
      tx_pf_contrast
    )
  ) &&
  all(
    tx_pf_contrast <
      0
  )

bridge_summary <-
  data.frame(
    metric = c(
      "mean_per_cell_CP10k",
      "pseudobulk_CP10k",
      "Pdgfb_positive_fraction"
    ),
    Sham1_dual_minus_notdual = c(
      sham_contrast[
        1
      ],
      sham_pb_contrast[
        1
      ],
      sham_pf_contrast[
        1
      ]
    ),
    Sham20_dual_minus_notdual = c(
      sham_contrast[
        2
      ],
      sham_pb_contrast[
        2
      ],
      sham_pf_contrast[
        2
      ]
    ),
    Tx17_dual_minus_notdual = c(
      tx_contrast[
        1
      ],
      tx_pb_contrast[
        1
      ],
      tx_pf_contrast[
        1
      ]
    ),
    Tx5_dual_minus_notdual = c(
      tx_contrast[
        2
      ],
      tx_pb_contrast[
        2
      ],
      tx_pf_contrast[
        2
      ]
    ),
    Sham_mean_dual_minus_notdual = c(
      safe_mean(
        sham_contrast
      ),
      safe_mean(
        sham_pb_contrast
      ),
      safe_mean(
        sham_pf_contrast
      )
    ),
    Tx_mean_dual_minus_notdual = c(
      safe_mean(
        tx_contrast
      ),
      safe_mean(
        tx_pb_contrast
      ),
      safe_mean(
        tx_pf_contrast
      )
    ),
    Tx_minus_Sham_difference_of_differences = c(
      safe_mean(
        tx_contrast
      ) -
      safe_mean(
        sham_contrast
      ),
      safe_mean(
        tx_pb_contrast
      ) -
      safe_mean(
        sham_pb_contrast
      ),
      safe_mean(
        tx_pf_contrast
      ) -
      safe_mean(
        sham_pf_contrast
      )
    ),
    Tx_replicates_same_direction = c(
      tx_same_direction_cp10k,
      all(
        is.finite(
          tx_pb_contrast
        )
      ) &&
      sign(
        tx_pb_contrast[
          1
        ]
      ) ==
      sign(
        tx_pb_contrast[
          2
        ]
      ),
      all(
        is.finite(
          tx_pf_contrast
        )
      ) &&
      sign(
        tx_pf_contrast[
          1
        ]
      ) ==
      sign(
        tx_pf_contrast[
          2
        ]
      )
    ),
    Tx_both_dual_positive_lower = c(
      tx_both_lower_cp10k,
      tx_both_lower_pseudobulk,
      tx_both_lower_positive_fraction
    ),
    stringsAsFactors = FALSE
  )

write_tsv(
  bridge_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_07_RepairResolution_IL10R_PDGFB_bridge_summary.tsv"
    )
  )
)

# ==============================================================================
# 8. All-five-subtype sensitivity
# ==============================================================================

all_status <-
  summarize_pdgfb(
    meta,
    c(
      "sample",
      "condition",
      "macrophage_class",
      "IL10R_status"
    )
  )

write_tsv(
  all_status,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_08_all_subtypes_Pdgfb_by_IL10R_status.tsv"
    )
  )
)

all_class_contrast_rows <- list()
aci <- 1L

for (
  class_now in CLASS_ORDER
) {
  for (
    sample_now in SAMPLE_ORDER
  ) {
    z <-
      all_status[
        all_status$macrophage_class ==
          class_now &
        all_status$sample ==
          sample_now,
        ,
        drop = FALSE
      ]

    dual_n <-
      get_status_value(
        z,
        "IL10R_dual_positive",
        "n_cells"
      )

    notdual_n <-
      get_status_value(
        z,
        "IL10R_not_dual_positive",
        "n_cells"
      )

    dual_cp10k <-
      get_status_value(
        z,
        "IL10R_dual_positive",
        "Pdgfb_mean_per_cell_CP10k"
      )

    notdual_cp10k <-
      get_status_value(
        z,
        "IL10R_not_dual_positive",
        "Pdgfb_mean_per_cell_CP10k"
      )

    all_class_contrast_rows[[aci]] <-
      data.frame(
        macrophage_class =
          class_now,
        sample =
          sample_now,
        condition =
          canonical_condition(
            sample_now
          ),
        dual_n =
          dual_n,
        notdual_n =
          notdual_n,
        both_groups_min20 =
          is.finite(
            dual_n
          ) &&
          is.finite(
            notdual_n
          ) &&
          dual_n >=
            20 &&
          notdual_n >=
            20,
        dual_minus_notdual_mean_per_cell_CP10k =
          dual_cp10k -
          notdual_cp10k,
        stringsAsFactors = FALSE
      )

    aci <- aci + 1L
  }
}

all_class_contrast <-
  do.call(
    rbind,
    all_class_contrast_rows
  )

write_tsv(
  all_class_contrast,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_09_all_subtypes_IL10R_Pdgfb_contrast.tsv"
    )
  )
)

class_summary_rows <- list()
csi <- 1L

for (
  class_now in CLASS_ORDER
) {
  z <-
    all_class_contrast[
      all_class_contrast$macrophage_class ==
        class_now,
      ,
      drop = FALSE
    ]

  sh <-
    z$dual_minus_notdual_mean_per_cell_CP10k[
      z$condition ==
        "Sham"
    ]

  tx <-
    z$dual_minus_notdual_mean_per_cell_CP10k[
      z$condition ==
        "Tx"
    ]

  class_summary_rows[[csi]] <-
    data.frame(
      macrophage_class =
        class_now,
      Sham_mean_dual_minus_notdual_CP10k =
        safe_mean(
          sh
        ),
      Tx_mean_dual_minus_notdual_CP10k =
        safe_mean(
          tx
        ),
      Tx_minus_Sham_difference_of_differences =
        safe_mean(
          tx
        ) -
        safe_mean(
          sh
        ),
      Tx_replicates_same_direction =
        all(
          is.finite(
            tx
          )
        ) &&
        sign(
          tx[
            1
          ]
        ) ==
        sign(
          tx[
            2
          ]
        ),
      Tx_both_dual_positive_lower =
        all(
          is.finite(
            tx
          )
        ) &&
        all(
          tx <
            0
        ),
      stringsAsFactors = FALSE
    )

  csi <- csi + 1L
}

class_summary <-
  do.call(
    rbind,
    class_summary_rows
  )

write_tsv(
  class_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_10_all_subtypes_IL10R_Pdgfb_summary.tsv"
    )
  )
)

# ==============================================================================
# 9. Interpretation boundary
# ==============================================================================

primary_bridge_call <-
  if (
    tx_both_lower_cp10k &&
    tx_both_lower_pseudobulk
  ) {
    "Tx17_and_Tx5_both_show_lower_Pdgfb_in_IL10R_dual_positive_RepairResolution_Mphi"
  } else {
    "No_replicate_consistent_direct_IL10R_status_to_Pdgfb_bridge"
  }

tx_specificity_call <-
  if (
    tx_both_lower_cp10k &&
    is.finite(
      safe_mean(
        tx_contrast
      ) -
      safe_mean(
        sham_contrast
      )
    ) &&
    (
      safe_mean(
        tx_contrast
      ) -
      safe_mean(
        sham_contrast
      )
    ) <
      0
  ) {
    "Tx_association_is_more_negative_than_Sham_descriptively"
  } else {
    "No_clear_Tx_specific_strengthening_of_IL10R_Pdgfb_association"
  }

interpretation <-
  data.frame(
    item = c(
      "primary_bridge_call",
      "Tx_specificity_call",
      "causal_boundary",
      "benchmark_boundary",
      "primary_metric",
      "secondary_metrics"
    ),
    value = c(
      primary_bridge_call,
      tx_specificity_call,
      "IL10R dual-positive status is an observational transcript-defined cell state. Even a replicated inverse association with Pdgfb is consistent with, but does not prove, IL10-driven suppression of Pdgfb.",
      "The v6.3.6 checkpoint benchmark is Sham mean 5.971 and Tx mean 3.568 CP10k for Repair/Resolution-Mphi per-cell Pdgfb expression. C2.2 recomputes mean per-cell CP10k directly from raw counts and each cell's total RNA UMI.",
      "Pdgfb_mean_per_cell_CP10k = mean over cells of Pdgfb_UMI / total_cell_RNA_UMI * 10000",
      "Pdgfb positive fraction and aggregate pseudobulk CP10k are reported as sensitivity metrics."
    ),
    stringsAsFactors = FALSE
  )

write_tsv(
  interpretation,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_11_interpretation.tsv"
    )
  )
)

provenance <-
  data.frame(
    key = c(
      "version",
      "input_RDS",
      "C21_IL10R_table",
      "sample_column",
      "final_class_column",
      "primary_class",
      "RNA_counts_layer",
      "statistics_boundary"
    ),
    value = c(
      VERSION,
      INPUT_RDS,
      C21_IL10R_TABLE,
      SAMPLE_COL,
      CLASS_COL,
      PRIMARY_CLASS,
      counts_layer,
      "Descriptive within-sample and n=2-per-condition summaries only; no cell-level inferential p-values."
    ),
    stringsAsFactors = FALSE
  )

write_tsv(
  provenance,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_12_provenance.tsv"
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
# 10. Console summary
# ==============================================================================

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "C2.1 IL10R FRACTION VALIDATION"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  c21_validation,
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
  "v6.3.6 REPAIR/RESOLUTION-Mphi PDGFB BENCHMARK"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  benchmark,
  row.names = FALSE
)

msg("")
msg(
  "Per-sample unsplit Repair/Resolution-Mphi Pdgfb:"
)

print(
  rr_unsplit[
    ,
    c(
      "sample",
      "n_cells",
      "Pdgfb_total_UMI",
      "Pdgfb_positive_fraction",
      "Pdgfb_mean_per_cell_CP10k",
      "Pdgfb_pseudobulk_CP10k"
    ),
    drop = FALSE
  ],
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
  "PRIMARY FOCUSED BRIDGE | Repair/Resolution-Mphi"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  rr_contrast,
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
  "BRIDGE SUMMARY"
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

msg("")
msg(
  "Primary bridge call: ",
  primary_bridge_call
)

msg(
  "Tx-specificity call: ",
  tx_specificity_call
)

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "ALL-Mphi SUBTYPE SENSITIVITY | mean per-cell CP10k"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  class_summary,
  row.names = FALSE
)

capture.output(
  print(
    c21_validation,
    row.names = FALSE
  ),
  print(
    benchmark,
    row.names = FALSE
  ),
  print(
    rr_unsplit,
    row.names = FALSE
  ),
  print(
    rr_contrast,
    row.names = FALSE
  ),
  print(
    bridge_summary,
    row.names = FALSE
  ),
  print(
    class_summary,
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
  "PASS: FINAL CLEAN-B Mphi v4.14.5 LOADED"
)

msg(
  "PASS: IL10R DUAL-POSITIVE STATE RECONSTRUCTED"
)

msg(
  "PASS: REPAIR/RESOLUTION-Mphi PDGFB PER-CELL AND PSEUDOBULK METRICS COMPUTED"
)

msg(
  "PASS: TX17/TX5 REPLICATE-CONSISTENT IL10R-PDGFB BRIDGE AUDITED"
)

msg(
  "PASS: ",
  VERSION,
  " FOCUSED IL10R-PDGFB BRIDGE AUDIT COMPLETE"
)

msg(
  strrep(
    "=",
    76
  )
)
