#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

VERSION <- "v1.4.3C.2.1"

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

REPO <-
  "/Users/uenoya/Projects/uenoyscRNA"

INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

PROGRAM_SOURCE <- file.path(
  REPO,
  "analysis",
  "Mouse_MASH_Mphi_FINAL_extended_figures_v4.16.1.R"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3C.2.1",
  "analysis",
  "v1.4.3C.2.1_mouse_Mphi_IL10_responder_mechanistic_audit"
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
      candidates %in% x
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

    op <-
      as.character(
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

    lhs <-
      e[[2]]

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
      PROGRAM_SOURCE
    )
  }
}

safe_z <- function(
  x
) {
  x <-
    as.numeric(
      x
    )

  s <-
    stats::sd(
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
  " | mouse macrophage IL10 responder mechanistic audit"
)

msg(
  strrep(
    "=",
    76
  )
)

required_files <- c(
  INPUT_RDS,
  PROGRAM_SOURCE
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
    "Required files missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
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
# 1. Recover exact validated v4.16.1 program definitions
# ==============================================================================

program_exprs <-
  parse(
    file = PROGRAM_SOURCE,
    keep.source = FALSE
  )

program_env <-
  new.env(
    parent = baseenv()
  )

extract_top_level_assignment(
  program_exprs,
  "PROGRAMS",
  program_env
)

PROGRAMS_ALL <-
  program_env$PROGRAMS

if (
  !is.list(
    PROGRAMS_ALL
  ) ||
    is.null(
      names(
        PROGRAMS_ALL
      )
    )
) {
  stop2(
    "PROGRAMS extracted from v4.16.1 is not a named list."
  )
}

preferred_programs <- c(
  "IL10_STAT3",
  "Repair_Resolution",
  "Inflammatory",
  "ECM_associated_inflammatory",
  "Efferocytosis",
  "Anti_inflammatory",
  "Lipid_TREM2"
)

program_names <-
  intersect(
    preferred_programs,
    names(
      PROGRAMS_ALL
    )
  )

critical_programs <- c(
  "IL10_STAT3",
  "Repair_Resolution",
  "Inflammatory",
  "ECM_associated_inflammatory",
  "Efferocytosis"
)

missing_critical <-
  setdiff(
    critical_programs,
    program_names
  )

if (
  length(
    missing_critical
  )
) {
  stop2(
    "Critical v4.16.1 programs missing: ",
    paste(
      missing_critical,
      collapse = ", "
    )
  )
}

PROGRAMS <-
  PROGRAMS_ALL[
    program_names
  ]

# Derived IL10 downstream-only panel removes the two receptor genes.
if (
  "IL10_STAT3" %in%
    names(
      PROGRAMS
    )
) {
  PROGRAMS[["IL10_downstream_no_receptors"]] <-
    setdiff(
      PROGRAMS[["IL10_STAT3"]],
      c(
        "Il10ra",
        "Il10rb"
      )
    )
}

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

PROGRAM_ORDER <-
  PROGRAM_ORDER[
    PROGRAM_ORDER %in%
      names(
        PROGRAMS
      )
  ]

msg(
  "Programs recovered from v4.16.1: ",
  paste(
    PROGRAM_ORDER,
    collapse = " | "
  )
)

# ==============================================================================
# 2. Load final Clean-B Mphi object and resolve metadata
# ==============================================================================

msg(
  "Loading FINAL Clean-B macrophage RDS: ",
  INPUT_RDS
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

if (
  is.na(
    SAMPLE_COL
  )
) {
  stop2(
    "Sample metadata column not found."
  )
}

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

meta <-
  data.frame(
    cell = colnames(
      mphi
    ),
    sample = as.character(
      mphi@meta.data[[SAMPLE_COL]]
    ),
    macrophage_class = as.character(
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
# 3. Program-gene audit
# ==============================================================================

program_gene_audit <-
  do.call(
    rbind,
    lapply(
      PROGRAM_ORDER,
      function(
        program_now
      ) {
        genes <-
          unique(
            PROGRAMS[[program_now]]
          )

        data.frame(
          program =
            program_now,
          gene =
            genes,
          detected =
            genes %in%
              rownames(
                counts
              ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

program_coverage <-
  do.call(
    rbind,
    lapply(
      split(
        program_gene_audit,
        program_gene_audit$program
      ),
      function(
        z
      ) {
        data.frame(
          program =
            z$program[
              1
            ],
          requested_genes =
            nrow(
              z
            ),
          detected_genes =
            sum(
              z$detected
            ),
          coverage_fraction =
            mean(
              z$detected
            ),
          missing_genes =
            paste(
              z$gene[
                !z$detected
              ],
              collapse = ";"
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

write_tsv(
  program_gene_audit,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_01_program_gene_audit.tsv"
    )
  )
)

write_tsv(
  program_coverage,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_02_program_coverage.tsv"
    )
  )
)

bad_programs <-
  program_coverage$program[
    program_coverage$detected_genes <
      2L
  ]

if (
  length(
    bad_programs
  )
) {
  stop2(
    "Programs with <2 detected genes: ",
    paste(
      bad_programs,
      collapse = ", "
    )
  )
}

# ==============================================================================
# 4. Sample-level class abundance
# ==============================================================================

class_counts <-
  as.data.frame(
    table(
      sample =
        factor(
          meta$sample,
          levels = SAMPLE_ORDER
        ),
      macrophage_class =
        factor(
          meta$macrophage_class,
          levels = CLASS_ORDER
        )
    ),
    stringsAsFactors = FALSE
  )

names(
  class_counts
) <- c(
  "sample",
  "macrophage_class",
  "n_cells"
)

sample_totals <-
  aggregate(
    n_cells ~ sample,
    data = class_counts,
    sum
  )

names(
  sample_totals
)[
  2
] <- "sample_total_Mphi"

class_counts <-
  merge(
    class_counts,
    sample_totals,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )

class_counts$fraction_total_Mphi <-
  class_counts$n_cells /
  class_counts$sample_total_Mphi

class_counts$condition <-
  canonical_condition(
    class_counts$sample
  )

write_tsv(
  class_counts,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_03_class_abundance_by_sample.tsv"
    )
  )
)

# ==============================================================================
# 5. IL10 receptor competence by sample and macrophage class
# ==============================================================================

required_il10_genes <- c(
  "Il10ra",
  "Il10rb"
)

missing_il10_genes <-
  setdiff(
    required_il10_genes,
    rownames(
      counts
    )
  )

if (
  length(
    missing_il10_genes
  )
) {
  stop2(
    "Required IL10 receptor genes absent: ",
    paste(
      missing_il10_genes,
      collapse = ", "
    )
  )
}

il10ra_pos <-
  as.vector(
    counts[
      "Il10ra",
      ,
      drop = TRUE
    ] >
      0
  )

il10rb_pos <-
  as.vector(
    counts[
      "Il10rb",
      ,
      drop = TRUE
    ] >
      0
  )

meta$Il10ra_positive <-
  il10ra_pos

meta$Il10rb_positive <-
  il10rb_pos

meta$IL10R_dual_positive <-
  il10ra_pos &
  il10rb_pos

meta$IL10R_status <-
  ifelse(
    meta$IL10R_dual_positive,
    "IL10R_dual_positive",
    "IL10R_not_dual_positive"
  )

group_key <-
  interaction(
    meta$sample,
    meta$macrophage_class,
    drop = TRUE,
    sep = "|||"
  )

group_idx <-
  split(
    seq_len(
      nrow(
        meta
      )
    ),
    group_key
  )

receptor_rows <- list()
ri <- 1L

for (
  idx in group_idx
) {
  z <-
    meta[
      idx,
      ,
      drop = FALSE
    ]

  receptor_rows[[ri]] <-
    data.frame(
      sample =
        z$sample[
          1
        ],
      condition =
        z$condition[
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
      Il10ra_positive_cells =
        sum(
          z$Il10ra_positive
        ),
      Il10ra_positive_fraction =
        mean(
          z$Il10ra_positive
        ),
      Il10rb_positive_cells =
        sum(
          z$Il10rb_positive
        ),
      Il10rb_positive_fraction =
        mean(
          z$Il10rb_positive
        ),
      IL10R_dual_positive_cells =
        sum(
          z$IL10R_dual_positive
        ),
      IL10R_dual_positive_fraction =
        mean(
          z$IL10R_dual_positive
        ),
      stringsAsFactors = FALSE
    )

  ri <- ri + 1L
}

receptor_by_sample_class <-
  do.call(
    rbind,
    receptor_rows
  )

receptor_by_sample_class <-
  receptor_by_sample_class[
    order(
      match(
        receptor_by_sample_class$sample,
        SAMPLE_ORDER
      ),
      match(
        receptor_by_sample_class$macrophage_class,
        CLASS_ORDER
      )
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  receptor_by_sample_class,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_04_IL10R_competence_by_sample_class.tsv"
    )
  )
)

# ==============================================================================
# 6. Generic pseudobulk helper
# ==============================================================================

program_gene_union <-
  unique(
    unlist(
      PROGRAMS[
        PROGRAM_ORDER
      ],
      use.names = FALSE
    )
  )

program_gene_union <-
  intersect(
    program_gene_union,
    rownames(
      counts
    )
  )

if (
  !length(
    program_gene_union
  )
) {
  stop2(
    "No program genes detected."
  )
}

aggregate_pb <- function(
  counts_mat,
  meta_df,
  group_cols,
  min_cells = 1L
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
  long_gene <- list()

  ri <- 1L
  gi <- 1L

  for (
    idx in idx_list
  ) {
    first <-
      meta_df[
        idx[
          1
        ],
        group_cols,
        drop = FALSE
      ]

    n_cells <-
      length(
        idx
      )

    lib_vec <-
      Matrix::rowSums(
        counts_mat[
          ,
          idx,
          drop = FALSE
        ]
      )

    lib_size <-
      sum(
        lib_vec
      )

    pg_vec <-
      lib_vec[
        match(
          program_gene_union,
          rownames(
            counts_mat
          )
        )
      ]

    names(
      pg_vec
    ) <- program_gene_union

    cpm <-
      if (
        lib_size >
          0
      ) {
        pg_vec /
          lib_size *
          1e6
      } else {
        rep(
          0,
          length(
            pg_vec
          )
        )
      }

    logcpm <-
      log2(
        cpm +
          1
      )

    group_name <-
      paste(
        unlist(
          first,
          use.names = FALSE
        ),
        collapse = "|||"
      )

    for (
      program_now in PROGRAM_ORDER
    ) {
      genes <-
        intersect(
          PROGRAMS[[program_now]],
          names(
            logcpm
          )
        )

      rows[[ri]] <-
        cbind(
          first,
          data.frame(
            group =
              group_name,
            program =
              program_now,
            n_cells =
              n_cells,
            library_UMI =
              lib_size,
            detected_program_genes =
              length(
                genes
              ),
            program_mean_log2CPM1 =
              if (
                length(
                  genes
                )
              ) {
                mean(
                  logcpm[
                    genes
                  ]
                )
              } else {
                NA_real_
              },
            min_cells_pass =
              n_cells >=
              min_cells,
            stringsAsFactors = FALSE
          )
        )

      ri <- ri + 1L
    }

    long_gene[[gi]] <-
      cbind(
        first[
          rep(
            1,
            length(
              logcpm
            )
          ),
          ,
          drop = FALSE
        ],
        data.frame(
          group =
            group_name,
          gene =
            names(
              logcpm
            ),
          n_cells =
            n_cells,
          library_UMI =
            lib_size,
          log2CPM1 =
            as.numeric(
              logcpm
            ),
          stringsAsFactors = FALSE
        )
      )

    gi <- gi + 1L
  }

  list(
    program =
      do.call(
        rbind,
        rows
      ),
    gene =
      do.call(
        rbind,
        long_gene
      )
  )
}

# ==============================================================================
# 7. Sample x class pseudobulk program remodeling
# ==============================================================================

pb_sample_class <-
  aggregate_pb(
    counts,
    meta,
    c(
      "sample",
      "condition",
      "macrophage_class"
    )
  )

gene_sc <-
  pb_sample_class$gene

gene_sc$gene_relative_z <-
  ave(
    gene_sc$log2CPM1,
    gene_sc$gene,
    FUN = safe_z
  )

relative_rows <- list()
rzi <- 1L

unique_group_rows <-
  unique(
    gene_sc[
      ,
      c(
        "sample",
        "condition",
        "macrophage_class",
        "group",
        "n_cells",
        "library_UMI"
      ),
      drop = FALSE
    ]
  )

for (
  ii in seq_len(
    nrow(
      unique_group_rows
    )
  )
) {
  ginfo <-
    unique_group_rows[
      ii,
      ,
      drop = FALSE
    ]

  z <-
    gene_sc[
      gene_sc$group ==
        ginfo$group,
      ,
      drop = FALSE
    ]

  for (
    program_now in PROGRAM_ORDER
  ) {
    genes <-
      intersect(
        PROGRAMS[[program_now]],
        z$gene
      )

    vals <-
      z$gene_relative_z[
        match(
          genes,
          z$gene
        )
      ]

    relative_rows[[rzi]] <-
      data.frame(
        sample =
          ginfo$sample,
        condition =
          ginfo$condition,
        macrophage_class =
          ginfo$macrophage_class,
        program =
          program_now,
        n_cells =
          ginfo$n_cells,
        library_UMI =
          ginfo$library_UMI,
        detected_program_genes =
          length(
            genes
          ),
        program_mean_gene_relative_z =
          safe_mean(
            vals
          ),
        stringsAsFactors = FALSE
      )

    rzi <- rzi + 1L
  }
}

program_relative_z <-
  do.call(
    rbind,
    relative_rows
  )

write_tsv(
  pb_sample_class$program,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_05_program_logCPM_by_sample_class.tsv"
    )
  )
)

write_tsv(
  program_relative_z,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_06_program_relative_z_by_sample_class.tsv"
    )
  )
)

# ==============================================================================
# 8. Sham vs Tx subtype-level program effects
# ==============================================================================

get_value <- function(
  z,
  sample,
  col
) {
  x <-
    z[
      z$sample ==
        sample,
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

effect_key <-
  interaction(
    program_relative_z$macrophage_class,
    program_relative_z$program,
    drop = TRUE,
    sep = "|||"
  )

effect_split <-
  split(
    seq_len(
      nrow(
        program_relative_z
      )
    ),
    effect_key
  )

effect_rows <- list()
ei <- 1L

for (
  idx in effect_split
) {
  z <-
    program_relative_z[
      idx,
      ,
      drop = FALSE
    ]

  sh <- c(
    get_value(
      z,
      "Sham1",
      "program_mean_gene_relative_z"
    ),
    get_value(
      z,
      "Sham20",
      "program_mean_gene_relative_z"
    )
  )

  tx <- c(
    get_value(
      z,
      "Tx17",
      "program_mean_gene_relative_z"
    ),
    get_value(
      z,
      "Tx5",
      "program_mean_gene_relative_z"
    )
  )

  effect <-
    safe_mean(
      tx
    ) -
    safe_mean(
      sh
    )

  complete <-
    all(
      is.finite(
        c(
          sh,
          tx
        )
      )
    )

  complete_separation <-
    if (
      complete &&
        effect >
          0
    ) {
      safe_min(
        tx
      ) >
        safe_max(
          sh
        )
    } else if (
      complete &&
        effect <
          0
    ) {
      safe_max(
        tx
      ) <
        safe_min(
          sh
        )
    } else {
      FALSE
    }

  effect_rows[[ei]] <-
    data.frame(
      macrophage_class =
        z$macrophage_class[
          1
        ],
      program =
        z$program[
          1
        ],
      Sham1_relative_z =
        sh[
          1
        ],
      Sham20_relative_z =
        sh[
          2
        ],
      Tx17_relative_z =
        tx[
          1
        ],
      Tx5_relative_z =
        tx[
          2
        ],
      Sham_mean_relative_z =
        safe_mean(
          sh
        ),
      Tx_mean_relative_z =
        safe_mean(
          tx
        ),
      Tx_minus_Sham_relative_z =
        effect,
      complete_four_sample_data =
        complete,
      complete_ShTx_separation =
        complete_separation,
      inference =
        "descriptive_effect_size_n2_per_condition",
      stringsAsFactors = FALSE
    )

  ei <- ei + 1L
}

program_effects <-
  do.call(
    rbind,
    effect_rows
  )

program_effects <-
  program_effects[
    order(
      match(
        program_effects$macrophage_class,
        CLASS_ORDER
      ),
      match(
        program_effects$program,
        PROGRAM_ORDER
      )
    ),
    ,
    drop = FALSE
  ]

write_tsv(
  program_effects,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_07_Sham_Tx_program_effect_by_class.tsv"
    )
  )
)

# ==============================================================================
# 9. All-Mphi sample-level program remodeling
# ==============================================================================

pb_all <-
  aggregate_pb(
    counts,
    meta,
    c(
      "sample",
      "condition"
    )
  )

gene_all <-
  pb_all$gene

gene_all$gene_relative_z <-
  ave(
    gene_all$log2CPM1,
    gene_all$gene,
    FUN = safe_z
  )

all_rel_rows <- list()
ari <- 1L

for (
  sample_now in SAMPLE_ORDER
) {
  z <-
    gene_all[
      gene_all$sample ==
        sample_now,
      ,
      drop = FALSE
    ]

  for (
    program_now in PROGRAM_ORDER
  ) {
    genes <-
      intersect(
        PROGRAMS[[program_now]],
        z$gene
      )

    vals <-
      z$gene_relative_z[
        match(
          genes,
          z$gene
        )
      ]

    all_rel_rows[[ari]] <-
      data.frame(
        sample =
          sample_now,
        condition =
          canonical_condition(
            sample_now
          ),
        program =
          program_now,
        program_mean_gene_relative_z =
          safe_mean(
            vals
          ),
        stringsAsFactors = FALSE
      )

    ari <- ari + 1L
  }
}

all_rel <-
  do.call(
    rbind,
    all_rel_rows
  )

write_tsv(
  all_rel,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_08_all_Mphi_program_relative_z_by_sample.tsv"
    )
  )
)

all_effect_rows <- list()
aei <- 1L

for (
  program_now in PROGRAM_ORDER
) {
  z <-
    all_rel[
      all_rel$program ==
        program_now,
      ,
      drop = FALSE
    ]

  sh <- c(
    get_value(
      z,
      "Sham1",
      "program_mean_gene_relative_z"
    ),
    get_value(
      z,
      "Sham20",
      "program_mean_gene_relative_z"
    )
  )

  tx <- c(
    get_value(
      z,
      "Tx17",
      "program_mean_gene_relative_z"
    ),
    get_value(
      z,
      "Tx5",
      "program_mean_gene_relative_z"
    )
  )

  effect <-
    safe_mean(
      tx
    ) -
    safe_mean(
      sh
    )

  all_effect_rows[[aei]] <-
    data.frame(
      program =
        program_now,
      Sham1_relative_z =
        sh[
          1
        ],
      Sham20_relative_z =
        sh[
          2
        ],
      Tx17_relative_z =
        tx[
          1
        ],
      Tx5_relative_z =
        tx[
          2
        ],
      Sham_mean_relative_z =
        safe_mean(
          sh
        ),
      Tx_mean_relative_z =
        safe_mean(
          tx
        ),
      Tx_minus_Sham_relative_z =
        effect,
      complete_ShTx_separation =
        if (
          all(
            is.finite(
              c(
                sh,
                tx
              )
            )
          ) &&
            effect >
              0
        ) {
          safe_min(
            tx
          ) >
            safe_max(
              sh
            )
        } else if (
          all(
            is.finite(
              c(
                sh,
                tx
              )
            )
          ) &&
            effect <
              0
        ) {
          safe_max(
            tx
          ) <
            safe_min(
              sh
            )
        } else {
          FALSE
        },
      inference =
        "descriptive_effect_size_n2_per_condition",
      stringsAsFactors = FALSE
    )

  aei <- aei + 1L
}

all_effects <-
  do.call(
    rbind,
    all_effect_rows
  )

write_tsv(
  all_effects,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_09_all_Mphi_Sham_Tx_program_effect.tsv"
    )
  )
)

# ==============================================================================
# 10. IL10R dual-positive vs not-dual-positive within sample x subtype
# ==============================================================================

status_counts <-
  as.data.frame(
    table(
      sample =
        factor(
          meta$sample,
          levels = SAMPLE_ORDER
        ),
      macrophage_class =
        factor(
          meta$macrophage_class,
          levels = CLASS_ORDER
        ),
      IL10R_status =
        factor(
          meta$IL10R_status,
          levels = c(
            "IL10R_dual_positive",
            "IL10R_not_dual_positive"
          )
        )
    ),
    stringsAsFactors = FALSE
  )

names(
  status_counts
) <- c(
  "sample",
  "macrophage_class",
  "IL10R_status",
  "n_cells"
)

status_counts$condition <-
  canonical_condition(
    status_counts$sample
  )

status_counts$min20_cells <-
  status_counts$n_cells >=
  20L

write_tsv(
  status_counts,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_10_IL10R_status_cell_counts.tsv"
    )
  )
)

pb_status <-
  aggregate_pb(
    counts,
    meta,
    c(
      "sample",
      "condition",
      "macrophage_class",
      "IL10R_status"
    ),
    min_cells = 20L
  )

write_tsv(
  pb_status$program,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_11_IL10R_status_program_logCPM_by_sample_class.tsv"
    )
  )
)

# Within-sample positive-minus-not-dual contrast.
status_prog <-
  pb_status$program

status_key <-
  interaction(
    status_prog$sample,
    status_prog$macrophage_class,
    status_prog$program,
    drop = TRUE,
    sep = "|||"
  )

status_split <-
  split(
    seq_len(
      nrow(
        status_prog
      )
    ),
    status_key
  )

status_contrast_rows <- list()
sci <- 1L

for (
  idx in status_split
) {
  z <-
    status_prog[
      idx,
      ,
      drop = FALSE
    ]

  pos <-
    z[
      z$IL10R_status ==
        "IL10R_dual_positive",
      ,
      drop = FALSE
    ]

  neg <-
    z[
      z$IL10R_status ==
        "IL10R_not_dual_positive",
      ,
      drop = FALSE
    ]

  if (
    nrow(
      pos
    ) !=
      1L ||
      nrow(
        neg
      ) !=
        1L
  ) {
    next
  }

  both_min20 <-
    isTRUE(
      pos$min_cells_pass
    ) &&
    isTRUE(
      neg$min_cells_pass
    )

  status_contrast_rows[[sci]] <-
    data.frame(
      sample =
        pos$sample,
      condition =
        pos$condition,
      macrophage_class =
        pos$macrophage_class,
      program =
        pos$program,
      IL10R_dual_positive_n_cells =
        pos$n_cells,
      IL10R_not_dual_positive_n_cells =
        neg$n_cells,
      IL10R_dual_positive_program_mean_log2CPM1 =
        pos$program_mean_log2CPM1,
      IL10R_not_dual_positive_program_mean_log2CPM1 =
        neg$program_mean_log2CPM1,
      dual_positive_minus_not_dual_log2CPM1 =
        pos$program_mean_log2CPM1 -
        neg$program_mean_log2CPM1,
      both_groups_min20 =
        both_min20,
      inference =
        "within_sample_descriptive_association_not_IL10_activation_proof",
      stringsAsFactors = FALSE
    )

  sci <- sci + 1L
}

status_contrasts <-
  do.call(
    rbind,
    status_contrast_rows
  )

write_tsv(
  status_contrasts,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_12_IL10R_status_program_contrast_within_sample.tsv"
    )
  )
)

# Summarize the within-sample IL10R-status contrast separately for Sham and Tx.
contrast_summary_rows <- list()
csi <- 1L

summary_key <-
  interaction(
    status_contrasts$macrophage_class,
    status_contrasts$program,
    drop = TRUE,
    sep = "|||"
  )

summary_split <-
  split(
    seq_len(
      nrow(
        status_contrasts
      )
    ),
    summary_key
  )

for (
  idx in summary_split
) {
  z <-
    status_contrasts[
      idx,
      ,
      drop = FALSE
    ]

  v <- function(
    sample
  ) {
    x <-
      z$dual_positive_minus_not_dual_log2CPM1[
        z$sample ==
          sample &
        z$both_groups_min20
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

  sh <- c(
    v(
      "Sham1"
    ),
    v(
      "Sham20"
    )
  )

  tx <- c(
    v(
      "Tx17"
    ),
    v(
      "Tx5"
    )
  )

  contrast_summary_rows[[csi]] <-
    data.frame(
      macrophage_class =
        z$macrophage_class[
          1
        ],
      program =
        z$program[
          1
        ],
      Sham1_dualpos_minus_notdual =
        sh[
          1
        ],
      Sham20_dualpos_minus_notdual =
        sh[
          2
        ],
      Tx17_dualpos_minus_notdual =
        tx[
          1
        ],
      Tx5_dualpos_minus_notdual =
        tx[
          2
        ],
      Sham_mean_dualpos_minus_notdual =
        safe_mean(
          sh
        ),
      Tx_mean_dualpos_minus_notdual =
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
        if (
          all(
            is.finite(
              tx
            )
          )
        ) {
          sign(
            tx[
              1
            ]
          ) ==
            sign(
              tx[
                2
              ]
            )
        } else {
          FALSE
        },
      inference =
        "descriptive_IL10R_state_association_n2_per_condition",
      stringsAsFactors = FALSE
    )

  csi <- csi + 1L
}

status_contrast_summary <-
  do.call(
    rbind,
    contrast_summary_rows
  )

write_tsv(
  status_contrast_summary,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_13_IL10R_status_program_contrast_summary.tsv"
    )
  )
)

# ==============================================================================
# 11. Key gene expression by sample and subtype
# ==============================================================================

key_genes <- unique(
  c(
    "Il10ra",
    "Il10rb",
    "Jak1",
    "Tyk2",
    "Stat3",
    "Socs3",
    "Bcl3",
    "Il1rn",
    "Il1b",
    "Tnf",
    "Cxcl10",
    "Ccl2",
    "Stat1",
    "Thbs1",
    "Fn1",
    "Tgfb1",
    "Mertk",
    "Mfge8",
    "Gas6",
    "Mmp13",
    "Mmp14",
    "Plau",
    "Trem2",
    "Gpnmb",
    "Cd9",
    "Lpl",
    "Apoe"
  )
)

key_genes <-
  intersect(
    key_genes,
    rownames(
      counts
    )
  )

key_rows <- list()
ki <- 1L

for (
  idx in group_idx
) {
  zmeta <-
    meta[
      idx,
      ,
      drop = FALSE
    ]

  block <-
    counts[
      key_genes,
      idx,
      drop = FALSE
    ]

  umi <-
    Matrix::rowSums(
      block
    )

  pos <-
    Matrix::rowSums(
      block >
        0
    )

  n_cells <-
    length(
      idx
    )

  key_rows[[ki]] <-
    data.frame(
      sample =
        zmeta$sample[
          1
        ],
      condition =
        zmeta$condition[
          1
        ],
      macrophage_class =
        zmeta$macrophage_class[
          1
        ],
      gene =
        key_genes,
      n_cells =
        n_cells,
      UMI =
        as.numeric(
          umi
        ),
      positive_cells =
        as.integer(
          pos
        ),
      positive_fraction =
        as.numeric(
          pos
        ) /
        n_cells,
      UMI_per_1000_cells =
        as.numeric(
          umi
        ) /
        n_cells *
        1000,
      stringsAsFactors = FALSE
    )

  ki <- ki + 1L
}

key_gene_expression <-
  do.call(
    rbind,
    key_rows
  )

write_tsv(
  key_gene_expression,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_14_key_gene_expression_by_sample_class.tsv"
    )
  )
)

# ==============================================================================
# 12. Interpretation boundaries / provenance
# ==============================================================================

boundaries <-
  data.frame(
    item = c(
      "input_RDS",
      "program_definition_source",
      "sample_column",
      "final_class_column",
      "statistics_boundary",
      "IL10_activation_boundary",
      "IL10R_status_boundary",
      "primary_mechanistic_question"
    ),
    value = c(
      INPUT_RDS,
      PROGRAM_SOURCE,
      SAMPLE_COL,
      CLASS_COL,
      "Sham1/Sham20 versus Tx17/Tx5: descriptive sample-level effect sizes and replicate concordance; no cell-level inferential p-values",
      "scRNA-seq measures transcript abundance/response competence; it does not measure IL10 protein exposure, receptor occupancy, JAK phosphorylation, or pSTAT3 activation",
      "IL10R_dual_positive means Il10ra>0 AND Il10rb>0 in the RNA count matrix; not-dual includes single-positive and double-negative cells",
      "Does the Tx macrophage compartment retain IL10R competence and remodel away from inflammatory/ECM-inflammatory programs toward anti-inflammatory/repair/efferocytosis-associated programs?"
    ),
    stringsAsFactors = FALSE
  )

write_tsv(
  boundaries,
  file.path(
    TAB_OUT,
    paste0(
      VERSION,
      "_15_interpretation_boundaries.tsv"
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
# 13. Console summaries
# ==============================================================================

msg("")
msg(
  strrep(
    "=",
    76
  )
)

msg(
  "PROGRAM COVERAGE"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  program_coverage[
    order(
      match(
        program_coverage$program,
        PROGRAM_ORDER
      )
    ),
    ,
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
  "IL10R DUAL-POSITIVE FRACTION BY SAMPLE x Mphi CLASS"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  receptor_by_sample_class[
    ,
    c(
      "sample",
      "macrophage_class",
      "n_cells",
      "Il10ra_positive_fraction",
      "Il10rb_positive_fraction",
      "IL10R_dual_positive_fraction"
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
  "ALL-Mphi SHAM -> TX PROGRAM REMODELING"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  all_effects[
    order(
      match(
        all_effects$program,
        PROGRAM_ORDER
      )
    ),
    ,
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
  "SUBTYPE-SPECIFIC SHAM -> TX PROGRAM EFFECTS"
)

msg(
  strrep(
    "=",
    76
  )
)

print(
  program_effects[
    ,
    c(
      "macrophage_class",
      "program",
      "Sham_mean_relative_z",
      "Tx_mean_relative_z",
      "Tx_minus_Sham_relative_z",
      "complete_ShTx_separation"
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
  "TX IL10R-STATE ASSOCIATION | dual-positive minus not-dual-positive"
)

msg(
  strrep(
    "=",
    76
  )
)

tx_assoc <-
  status_contrast_summary[
    ,
    c(
      "macrophage_class",
      "program",
      "Tx17_dualpos_minus_notdual",
      "Tx5_dualpos_minus_notdual",
      "Tx_mean_dualpos_minus_notdual",
      "Tx_replicates_same_direction",
      "Tx_minus_Sham_difference_of_differences"
    ),
    drop = FALSE
  ]

print(
  tx_assoc,
  row.names = FALSE
)

capture.output(
  print(
    program_coverage,
    row.names = FALSE
  ),
  print(
    receptor_by_sample_class,
    row.names = FALSE
  ),
  print(
    all_effects,
    row.names = FALSE
  ),
  print(
    program_effects,
    row.names = FALSE
  ),
  print(
    tx_assoc,
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
  "PASS: EXACT v4.16.1 PROGRAM DEFINITIONS RECOVERED"
)

msg(
  "PASS: SAMPLE-LEVEL IL10R COMPETENCE AND PROGRAM REMODELING AUDITED"
)

msg(
  "PASS: WITHIN-SAMPLE IL10R-STATE ASSOCIATIONS AUDITED"
)

msg(
  "PASS: ",
  VERSION,
  " MOUSE Mphi IL10 RESPONDER MECHANISTIC AUDIT COMPLETE"
)

msg(
  strrep(
    "=",
    76
  )
)
