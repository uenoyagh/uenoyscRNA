#!/usr/bin/env Rscript

# ==============================================================================
# Human MASH graft all-ligandome
# v7.2.0.3 | Tx17 / Tx5 ligand-agnostic sender atlas
#
# PURPOSE
#   Build an unbiased human-graft ligand expression atlas from Tx17 and Tx5,
#   without privileging IL10 or the previous 15-ligand shortlist.
#
# PRIMARY PRINCIPLES
#   - Sender-first within human graft: all NicheNet-v2 human LR ligands.
#   - Tx17 and Tx5 are biological replicate graft samples.
#   - No inferential test is performed with n=2.
#   - Standard filtered MEX is used for DISCOVERY.
#   - Absence in standard MEX is NOT interpreted as biological absence if
#     probe-filter status has not been audited.
#   - Human and mouse raw UMI values are NEVER compared directly.
#   - Cross-species compatibility is assessed by strict unique human-mouse
#     orthology (babelgene, >=3 supporting databases) plus NicheNet resources.
#   - IL10 is audit-only and receives NO scoring/ranking advantage.
#
# UPSTREAM HOST TARGETS
#   Mouse Host Response Atlas v7.1.2 has already defined:
#     - contamination-clean CORE targets
#     - contamination-clean EXPANDED targets
#   Those targets are NOT used to select human ligands in this script.
#   They will be integrated in v7.2.1/v7.2.2 after receptor competence.
#
# OUTPUT
#   /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/
#     Human_MASH_RDS/Human_Graft_AllLigandome_v7.2.0.3/
# ==============================================================================

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v7.2.0.3"

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
REPO <- "/Users/uenoya/Projects/uenoyscRNA"

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  paste0("Human_Graft_AllLigandome_", VERSION)
)

TAB_OUT <- file.path(OUT, "tables")
OBJ_OUT <- file.path(OUT, "objects")
LOG_OUT <- file.path(OUT, "logs")
FIG_OUT <- file.path(OUT, "figures")

for (d in c(TAB_OUT, OBJ_OUT, LOG_OUT, FIG_OUT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  LOG_OUT,
  paste0("Human_Graft_AllLigandome_", VERSION, "_run.log")
)

cat("", file = LOG_FILE)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOG_FILE, append = TRUE)
}

stop2 <- function(...) {
  x <- paste0(...)
  cat("ERROR: ", x, "\n", sep = "", file = LOG_FILE, append = TRUE)
  stop(x, call. = FALSE)
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

first_present <- function(candidates, nms) {
  hit <- candidates[candidates %in% nms]
  if (!length(hit)) {
    return(NA_character_)
  }
  hit[[1]]
}

norm_sample <- function(x) {
  x <- trimws(as.character(x))
  x[x == "SfLB_1"] <- "SfLB1"
  x[x == "SfLB_2"] <- "SfLB2"
  x
}

extract_barcode_core_one <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }

  hits <- regmatches(
    toupper(x),
    gregexpr(
      "[ACGT]{16,40}",
      toupper(x),
      perl = TRUE
    )
  )[[1]]

  if (!length(hits)) {
    return(NA_character_)
  }

  hits[[which.max(nchar(hits))]]
}

extract_barcode_core <- function(x) {
  vapply(
    x,
    extract_barcode_core_one,
    character(1)
  )
}

safe_divide <- function(a, b) {
  out <- rep(NA_real_, length(a))
  ok <- is.finite(a) & is.finite(b) & b != 0
  out[ok] <- a[ok] / b[ok]
  out
}

# ==============================================================================
# 1. Inputs
# ==============================================================================

CONTEXT_META <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.3.1",
  "analysis",
  "v1.4.3.1_12sample_integrated_context_review",
  "tables",
  "v1.4.3_context_UMAP_metadata.tsv"
)

MEX_DIRS <- c(
  Tx17 = file.path(
    ROOT,
    "Human_MASH_RDS",
    "RH251217117_raw_MEX",
    "Tx17",
    "filtered_feature_bc_matrix"
  ),
  Tx5 = file.path(
    ROOT,
    "Human_MASH_RDS",
    "RH251217117_raw_MEX",
    "Tx5",
    "filtered_feature_bc_matrix"
  )
)

EXPECTED_RAW_MEX_CELLS <- c(
  Tx17 = 5959L,
  Tx5 = 8571L
)

EXPECTED_ANALYSIS_CELLS <- c(
  Tx17 = 5958L,
  Tx5 = 8571L
)

EXPECTED_FEATURES <- c(
  Tx17 = 18132L,
  Tx5 = 18132L
)

SAMPLES <- names(MEX_DIRS)

# NicheNet-v2 official final resources.
RESOURCE_DIR <- file.path(
  ROOT,
  "resources",
  "NicheNet_v2_Zenodo7074291"
)

dir.create(
  RESOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

RESOURCE_SPECS <- data.frame(
  name = c(
    "lr_human",
    "lr_mouse",
    "ligand_target_mouse"
  ),
  file = c(
    "lr_network_human_21122021.rds",
    "lr_network_mouse_21122021.rds",
    "ligand_target_matrix_nsga2r_final_mouse.rds"
  ),
  url = c(
    "https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds",
    "https://zenodo.org/records/7074291/files/lr_network_mouse_21122021.rds",
    "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
  ),
  md5 = c(
    "2a155f81e9ffffd5d5e709fe66bcc465",
    "cf33ee8b6bf84bdf2d11cab9c8f94b9e",
    "ac80d846fe0bfc4879a5b52ca85ffeb9"
  ),
  stringsAsFactors = FALSE
)

RESOURCE_SPECS$path <- file.path(
  RESOURCE_DIR,
  RESOURCE_SPECS$file
)

# Previous candidate panel is retained ONLY as an audit label.
LEGACY_15 <- c(
  "IL10",
  "FGF1",
  "VEGFA",
  "CXCL12",
  "TNFSF12",
  "FST",
  "IGFBP2",
  "IGFBP4",
  "INHBA",
  "BMP2",
  "BMP5",
  "TGFB3",
  "HBEGF",
  "VEGFB",
  "BTC"
)

# ==============================================================================
# 2. Preflight
# ==============================================================================

msg(strrep("=", 92))
msg("Human MASH graft all-ligandome ", VERSION)
msg(strrep("=", 92))
msg("Samples: Tx17 + Tx5")
msg("Ligand prior: NONE within NicheNet-v2 human LR ligand universe")
msg("IL10 privileged: FALSE")
msg("Inference: descriptive only (n=2 graft samples)")

required_paths <- c(
  CONTEXT_META,
  unname(MEX_DIRS)
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths)) {
  stop2(
    "Required input path(s) missing:\n",
    paste(missing_paths, collapse = "\n")
  )
}

required_pkgs <- c(
  "Seurat",
  "Matrix",
  "babelgene",
  "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs)) {
  stop2(
    "Required package(s) missing: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

for (p in required_pkgs) {
  msg(
    p,
    " version: ",
    as.character(packageVersion(p))
  )
}

# ==============================================================================
# 3. Download / validate NicheNet-v2 resources
# ==============================================================================

# Large mouse ligand-target matrix is ~182 MB.
# Local validated cache is preferred; this timeout only protects
# reproducibility if a future re-download becomes necessary.
options(
  timeout = max(
    1800,
    getOption("timeout", 60)
  )
)

download_and_validate <- function(path, url, md5_expected) {
  need_download <- !file.exists(path)

  if (!need_download) {
    md5_now <- unname(
      tools::md5sum(path)
    )

    if (!identical(
      tolower(md5_now),
      tolower(md5_expected)
    )) {
      msg(
        "Existing resource MD5 mismatch; re-downloading: ",
        basename(path)
      )
      unlink(path)
      need_download <- TRUE
    }
  }

  if (need_download) {
    msg(
      "Downloading NicheNet resource: ",
      basename(path)
    )

    tmp <- paste0(
      path,
      ".partial"
    )

    if (file.exists(tmp)) {
      unlink(tmp)
    }

    utils::download.file(
      url = url,
      destfile = tmp,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    )

    if (!file.exists(tmp)) {
      stop2(
        "Download failed: ",
        url
      )
    }

    md5_download <- unname(
      tools::md5sum(tmp)
    )

    if (!identical(
      tolower(md5_download),
      tolower(md5_expected)
    )) {
      unlink(tmp)
      stop2(
        "Downloaded resource MD5 mismatch for ",
        basename(path),
        " | expected=",
        md5_expected,
        " | observed=",
        md5_download
      )
    }

    if (!file.rename(
      tmp,
      path
    )) {
      stop2(
        "Could not move validated resource into cache: ",
        path
      )
    }
  }

  md5_final <- unname(
    tools::md5sum(path)
  )

  if (!identical(
    tolower(md5_final),
    tolower(md5_expected)
  )) {
    stop2(
      "Resource validation failed after cache step: ",
      path
    )
  }

  invisible(path)
}

for (i in seq_len(nrow(RESOURCE_SPECS))) {
  download_and_validate(
    path = RESOURCE_SPECS$path[i],
    url = RESOURCE_SPECS$url[i],
    md5_expected = RESOURCE_SPECS$md5[i]
  )
}

write_tsv(
  RESOURCE_SPECS,
  file.path(
    TAB_OUT,
    paste0(
      "00A_NicheNet_resource_manifest_",
      VERSION,
      ".tsv"
    )
  )
)

lr_human <- readRDS(
  RESOURCE_SPECS$path[
    RESOURCE_SPECS$name == "lr_human"
  ]
)

lr_mouse <- readRDS(
  RESOURCE_SPECS$path[
    RESOURCE_SPECS$name == "lr_mouse"
  ]
)

ligand_target_mouse <- readRDS(
  RESOURCE_SPECS$path[
    RESOURCE_SPECS$name == "ligand_target_mouse"
  ]
)

for (nm in c("from", "to")) {
  if (!nm %in% colnames(lr_human)) {
    stop2(
      "Human LR network missing column: ",
      nm
    )
  }

  if (!nm %in% colnames(lr_mouse)) {
    stop2(
      "Mouse LR network missing column: ",
      nm
    )
  }
}

if (is.null(
  colnames(
    ligand_target_mouse
  )
)) {
  stop2(
    "Mouse ligand-target matrix has no ligand column names."
  )
}

human_lr <- unique(
  data.frame(
    human_ligand = as.character(
      lr_human$from
    ),
    human_receptor = as.character(
      lr_human$to
    ),
    stringsAsFactors = FALSE
  )
)

mouse_lr <- unique(
  data.frame(
    mouse_ligand = as.character(
      lr_mouse$from
    ),
    mouse_receptor = as.character(
      lr_mouse$to
    ),
    stringsAsFactors = FALSE
  )
)

human_ligands <- sort(
  unique(
    human_lr$human_ligand[
      !is.na(human_lr$human_ligand) &
        nzchar(human_lr$human_ligand)
    ]
  )
)

mouse_lt_ligands <- unique(
  colnames(
    ligand_target_mouse
  )
)

msg(
  "NicheNet human ligand universe: ",
  length(human_ligands)
)

msg(
  "NicheNet mouse LR ligands: ",
  length(
    unique(
      mouse_lr$mouse_ligand
    )
  )
)

msg(
  "NicheNet mouse ligand-target ligands: ",
  length(mouse_lt_ligands)
)

# ==============================================================================
# 4. Load canonical context metadata
# ==============================================================================

context_raw <- read.delim(
  CONTEXT_META,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

sample_col <- first_present(
  c(
    "sample",
    "sample_v143A1"
  ),
  colnames(context_raw)
)

barcode_col <- first_present(
  c(
    "barcode_core",
    "barcode_core_v143A1"
  ),
  colnames(context_raw)
)

annotation_col <- first_present(
  c(
    "context_annotation",
    "context_annotation_v143A1",
    "harmonized_lineage_v143A1"
  ),
  colnames(context_raw)
)

state_col <- first_present(
  c(
    "context_state",
    "context_state_v143A1"
  ),
  colnames(context_raw)
)

resolved_cols <- c(
  sample = sample_col,
  barcode_core = barcode_col,
  context_annotation = annotation_col,
  context_state = state_col
)

if (anyNA(resolved_cols)) {
  stop2(
    "Could not resolve canonical context columns. Resolved: ",
    paste(
      names(resolved_cols),
      resolved_cols,
      sep = "=",
      collapse = "; "
    )
  )
}

context <- data.frame(
  sample = norm_sample(
    context_raw[[sample_col]]
  ),
  barcode_core = as.character(
    context_raw[[barcode_col]]
  ),
  context_annotation = as.character(
    context_raw[[annotation_col]]
  ),
  context_state = as.character(
    context_raw[[state_col]]
  ),
  stringsAsFactors = FALSE
)

context <- context[
  context$sample %in% SAMPLES,
  ,
  drop = FALSE
]

for (nm in c(
  "barcode_core",
  "context_annotation",
  "context_state"
)) {
  context[[nm]][
    is.na(context[[nm]]) |
      !nzchar(context[[nm]])
  ] <- "Unresolved"
}

context$barcode_core[
  context$barcode_core == "Unresolved"
] <- NA_character_

if (anyNA(
  context$barcode_core
)) {
  stop2(
    "Canonical context contains NA/unresolved barcode_core values."
  )
}

if (anyDuplicated(
  paste(
    context$sample,
    context$barcode_core,
    sep = "__"
  )
)) {
  stop2(
    "Duplicated sample + barcode_core pairs in canonical context."
  )
}

context_qc <- list()

for (s in SAMPLES) {
  n <- sum(
    context$sample == s
  )

  if (n != EXPECTED_ANALYSIS_CELLS[[s]]) {
    stop2(
      s,
      ": canonical context cells=",
      n,
      " expected=",
      EXPECTED_ANALYSIS_CELLS[[s]]
    )
  }

  context_qc[[s]] <- data.frame(
    sample = s,
    canonical_context_cells = n,
    expected_cells = EXPECTED_ANALYSIS_CELLS[[s]],
    stringsAsFactors = FALSE
  )
}

msg(
  "PASS: canonical Tx17/Tx5 context counts reproduced."
)

# ==============================================================================
# 5. Read Tx17 / Tx5 filtered MEX and map canonical barcodes
# ==============================================================================

read_gene_expression_matrix <- function(data_dir) {
  x <- Seurat::Read10X(
    data.dir = data_dir,
    gene.column = 2,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x <- x[["Gene Expression"]]
    } else if (length(x) == 1L) {
      x <- x[[1]]
    } else {
      stop2(
        "Read10X returned multiple modalities but no 'Gene Expression': ",
        paste(
          names(x),
          collapse = ", "
        )
      )
    }
  }

  if (!inherits(
    x,
    "Matrix"
  )) {
    x <- Matrix::Matrix(
      x,
      sparse = TRUE
    )
  }

  x
}

reconcile_sample_to_context <- function(
  mat,
  sample_name
) {

  cores <- extract_barcode_core(
    colnames(mat)
  )

  if (anyNA(cores)) {
    stop2(
      sample_name,
      ": failed barcode-core extraction for ",
      sum(is.na(cores)),
      " matrix cells."
    )
  }

  if (anyDuplicated(cores)) {
    stop2(
      sample_name,
      ": duplicated extracted barcode cores in raw MEX."
    )
  }

  meta <- context[
    context$sample == sample_name,
    ,
    drop = FALSE
  ]

  ## ----------------------------------------------------------
  ## Canonical cells must all exist in the raw filtered MEX.
  ## ----------------------------------------------------------

  missing_from_mex <- setdiff(
    meta$barcode_core,
    cores
  )

  if (length(missing_from_mex)) {
    stop2(
      sample_name,
      ": canonical cells absent from raw filtered MEX: ",
      length(missing_from_mex),
      " | first=",
      paste(
        head(missing_from_mex, 10),
        collapse = ","
      )
    )
  }

  ## ----------------------------------------------------------
  ## Raw-MEX-only cells are audited and excluded from the
  ## primary analysis cell universe.
  ## ----------------------------------------------------------

  idx_context <- match(
    cores,
    meta$barcode_core
  )

  keep <- !is.na(idx_context)

  excluded_idx <- which(!keep)

  if (length(excluded_idx)) {

    excluded_mat <- mat[
      ,
      excluded_idx,
      drop = FALSE
    ]

    excluded <- data.frame(
      sample = sample_name,
      raw_barcode = colnames(mat)[excluded_idx],
      barcode_core = cores[excluded_idx],
      total_UMI = as.numeric(
        Matrix::colSums(excluded_mat)
      ),
      nFeature = as.integer(
        Matrix::colSums(excluded_mat > 0)
      ),
      exclusion_reason =
        "raw_MEX_cell_absent_from_canonical_context",
      stringsAsFactors = FALSE
    )

    msg(
      sample_name,
      ": excluding ",
      nrow(excluded),
      " raw-MEX-only cell(s): ",
      paste(
        excluded$raw_barcode,
        collapse = ","
      )
    )

  } else {

    excluded <- data.frame(
      sample = character(0),
      raw_barcode = character(0),
      barcode_core = character(0),
      total_UMI = numeric(0),
      nFeature = integer(0),
      exclusion_reason = character(0),
      stringsAsFactors = FALSE
    )
  }

  mat_use <- mat[
    ,
    keep,
    drop = FALSE
  ]

  cores_use <- cores[keep]

  idx_use <- match(
    cores_use,
    meta$barcode_core
  )

  if (anyNA(idx_use)) {
    stop2(
      sample_name,
      ": internal canonical matching failure after filtering."
    )
  }

  aligned <- meta[
    idx_use,
    ,
    drop = FALSE
  ]

  if (!identical(
    unname(cores_use),
    unname(aligned$barcode_core)
  )) {
    stop2(
      sample_name,
      ": barcode alignment failed after canonical filtering."
    )
  }

  if (ncol(mat_use) != nrow(meta)) {
    stop2(
      sample_name,
      ": canonical matrix cells=",
      ncol(mat_use),
      " context cells=",
      nrow(meta)
    )
  }

  rownames(aligned) <- colnames(mat_use)

  list(
    mat = mat_use,
    meta = aligned,
    excluded = excluded
  )
}

mats <- list()
meta_by_sample <- list()
mex_qc_rows <- list()
excluded_qc_rows <- list()

for (s in SAMPLES) {

  msg(strrep("-", 92))

  msg(
    "Loading raw filtered MEX: ",
    s
  )

  mat_raw <- read_gene_expression_matrix(
    MEX_DIRS[[s]]
  )

  msg(
    "  raw features=",
    nrow(mat_raw),
    " raw cells=",
    ncol(mat_raw)
  )

  if (
    ncol(mat_raw) !=
      EXPECTED_RAW_MEX_CELLS[[s]]
  ) {
    stop2(
      s,
      ": raw filtered MEX cells=",
      ncol(mat_raw),
      " expected raw=",
      EXPECTED_RAW_MEX_CELLS[[s]]
    )
  }

  if (
    nrow(mat_raw) !=
      EXPECTED_FEATURES[[s]]
  ) {
    stop2(
      s,
      ": raw filtered MEX features=",
      nrow(mat_raw),
      " expected=",
      EXPECTED_FEATURES[[s]]
    )
  }

  rec <- reconcile_sample_to_context(
    mat_raw,
    s
  )

  mat <- rec$mat
  meta_aligned <- rec$meta

  if (
    ncol(mat) !=
      EXPECTED_ANALYSIS_CELLS[[s]]
  ) {
    stop2(
      s,
      ": canonical analysis cells=",
      ncol(mat),
      " expected=",
      EXPECTED_ANALYSIS_CELLS[[s]]
    )
  }

  msg(
    "  canonical analysis cells=",
    ncol(mat),
    " excluded raw-MEX-only=",
    nrow(rec$excluded)
  )

  mats[[s]] <- mat
  meta_by_sample[[s]] <- meta_aligned
  excluded_qc_rows[[s]] <- rec$excluded

  mex_qc_rows[[s]] <- data.frame(
    sample = s,

    n_features = nrow(mat_raw),

    raw_mex_cells = ncol(mat_raw),

    expected_raw_mex_cells =
      EXPECTED_RAW_MEX_CELLS[[s]],

    n_cells = ncol(mat),

    expected_cells =
      EXPECTED_ANALYSIS_CELLS[[s]],

    excluded_noncanonical_cells =
      nrow(rec$excluded),

    expected_features =
      EXPECTED_FEATURES[[s]],

    canonical_barcode_mapping =
      "PASS_context_intersection",

    total_UMI = sum(mat),

    median_UMI_per_cell =
      stats::median(
        Matrix::colSums(mat)
      ),

    stringsAsFactors = FALSE
  )
}

mex_qc <- do.call(
  rbind,
  mex_qc_rows
)

rownames(mex_qc) <- NULL

excluded_qc <- do.call(
  rbind,
  excluded_qc_rows
)

rownames(excluded_qc) <- NULL

write_tsv(
  excluded_qc,
  file.path(
    TAB_OUT,
    paste0(
      "01A_Tx_MEX_excluded_noncanonical_cells_",
      VERSION,
      ".tsv"
    )
  )
)

write_tsv(
  mex_qc,
  file.path(
    TAB_OUT,
    paste0(
      "01_Tx_MEX_context_mapping_QC_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 6. Human ligand feature-space audit
# ==============================================================================

feature_space_rows <- list()

for (gene in human_ligands) {
  present_tx17 <- gene %in% rownames(
    mats$Tx17
  )

  present_tx5 <- gene %in% rownames(
    mats$Tx5
  )

  feature_space_rows[[gene]] <- data.frame(
    human_ligand = gene,
    in_human_NicheNet_LR = TRUE,
    Tx17_feature_measurable = present_tx17,
    Tx5_feature_measurable = present_tx5,
    shared_Tx_feature_measurable =
      present_tx17 && present_tx5,
    stringsAsFactors = FALSE
  )
}

feature_space <- do.call(
  rbind,
  feature_space_rows
)

rownames(feature_space) <- NULL

feature_space$n_human_LR_receptors <- vapply(
  feature_space$human_ligand,
  function(g) {
    length(
      unique(
        human_lr$human_receptor[
          human_lr$human_ligand == g
        ]
      )
    )
  },
  integer(1)
)

feature_space$legacy_15_panel_audit_only <-
  feature_space$human_ligand %in% LEGACY_15

feature_space$is_IL10_audit_only <-
  feature_space$human_ligand == "IL10"

write_tsv(
  feature_space,
  file.path(
    TAB_OUT,
    paste0(
      "02_human_NicheNet_ligand_feature_space_",
      VERSION,
      ".tsv"
    )
  )
)

msg(
  "Human LR ligands measurable in BOTH Tx samples: ",
  sum(
    feature_space$shared_Tx_feature_measurable
  ),
  "/",
  nrow(feature_space)
)

# ==============================================================================
# 7. Per-sample all-ligand quantification
# ==============================================================================

quantify_sample <- function(
  mat,
  sample_name,
  ligand_universe
) {
  n_cells <- ncol(mat)
  total_library <- sum(mat)

  out <- data.frame(
    sample = sample_name,
    human_ligand = ligand_universe,
    feature_measurable = ligand_universe %in% rownames(mat),
    UMI = NA_real_,
    positive_cells = NA_integer_,
    positive_fraction = NA_real_,
    UMI_per_1000_cells = NA_real_,
    pseudobulk_CPM = NA_real_,
    mean_UMI_per_cell = NA_real_,
    mean_UMI_per_positive_cell = NA_real_,
    stringsAsFactors = FALSE
  )

  genes <- out$human_ligand[
    out$feature_measurable
  ]

  if (!length(genes)) {
    return(out)
  }

  x <- mat[
    genes,
    ,
    drop = FALSE
  ]

  umi <- Matrix::rowSums(x)

  pos <- Matrix::rowSums(
    x > 0
  )

  idx <- match(
    genes,
    out$human_ligand
  )

  out$UMI[idx] <- as.numeric(umi)

  out$positive_cells[idx] <- as.integer(pos)

  out$positive_fraction[idx] <-
    as.numeric(pos) / n_cells

  out$UMI_per_1000_cells[idx] <-
    as.numeric(umi) / n_cells * 1000

  out$pseudobulk_CPM[idx] <-
    as.numeric(umi) / total_library * 1e6

  out$mean_UMI_per_cell[idx] <-
    as.numeric(umi) / n_cells

  out$mean_UMI_per_positive_cell[idx] <-
    ifelse(
      pos > 0,
      as.numeric(umi) / as.numeric(pos),
      0
    )

  out
}

sample_quant_list <- lapply(
  SAMPLES,
  function(s) {
    quantify_sample(
      mats[[s]],
      s,
      human_ligands
    )
  }
)

names(sample_quant_list) <- SAMPLES

sample_quant <- do.call(
  rbind,
  sample_quant_list
)

rownames(sample_quant) <- NULL

write_tsv(
  sample_quant,
  file.path(
    TAB_OUT,
    paste0(
      "03_human_graft_all_ligands_by_sample_long_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 8. Lineage/state source atlas using FULL frozen denominators
# ==============================================================================

quantify_groups <- function(
  mat,
  meta,
  sample_name,
  ligand_genes,
  group_col,
  group_type,
  sample_gene_umi
) {
  groups <- sort(
    unique(
      meta[[group_col]]
    )
  )

  cell_library <- Matrix::colSums(
    mat
  )

  out_list <- vector(
    "list",
    length(groups)
  )

  for (i in seq_along(groups)) {
    grp <- groups[[i]]

    idx_cells <- which(
      meta[[group_col]] == grp
    )

    n_group <- length(
      idx_cells
    )

    group_library <- sum(
      cell_library[
        idx_cells
      ]
    )

    x <- mat[
      ligand_genes,
      idx_cells,
      drop = FALSE
    ]

    umi <- Matrix::rowSums(
      x
    )

    pos <- Matrix::rowSums(
      x > 0
    )

    total_gene_umi <- sample_gene_umi[
      ligand_genes
    ]

    out_list[[i]] <- data.frame(
      sample = sample_name,
      group_type = group_type,
      group = grp,
      human_ligand = ligand_genes,
      full_denominator_cells = n_group,
      UMI = as.numeric(umi),
      positive_cells = as.integer(pos),
      positive_fraction_full_denominator =
        as.numeric(pos) / n_group,
      UMI_per_1000_cells =
        as.numeric(umi) / n_group * 1000,
      pseudobulk_CPM_within_group =
        if (group_library > 0) {
          as.numeric(umi) /
            group_library *
            1e6
        } else {
          NA_real_
        },
      fraction_of_sample_ligand_UMI =
        ifelse(
          total_gene_umi > 0,
          as.numeric(umi) /
            as.numeric(total_gene_umi),
          0
        ),
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    out_list
  )
}

shared_ligands <- feature_space$human_ligand[
  feature_space$shared_Tx_feature_measurable
]

lineage_rows <- list()
state_rows <- list()

for (s in SAMPLES) {
  mat <- mats[[s]]
  meta <- meta_by_sample[[s]]

  sample_counts <- Matrix::rowSums(
    mat[
      shared_ligands,
      ,
      drop = FALSE
    ]
  )

  names(sample_counts) <- shared_ligands

  lineage_rows[[s]] <- quantify_groups(
    mat = mat,
    meta = meta,
    sample_name = s,
    ligand_genes = shared_ligands,
    group_col = "context_annotation",
    group_type = "lineage",
    sample_gene_umi = sample_counts
  )

  state_rows[[s]] <- quantify_groups(
    mat = mat,
    meta = meta,
    sample_name = s,
    ligand_genes = shared_ligands,
    group_col = "context_state",
    group_type = "state",
    sample_gene_umi = sample_counts
  )
}

lineage_atlas <- do.call(
  rbind,
  lineage_rows
)

rownames(lineage_atlas) <- NULL

state_atlas <- do.call(
  rbind,
  state_rows
)

rownames(state_atlas) <- NULL

write_tsv(
  lineage_atlas,
  file.path(
    TAB_OUT,
    paste0(
      "04_human_graft_all_ligands_lineage_source_atlas_",
      VERSION,
      ".tsv"
    )
  )
)

write_tsv(
  state_atlas,
  file.path(
    TAB_OUT,
    paste0(
      "05_human_graft_all_ligands_state_source_atlas_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 9. Dominant source summaries
# ==============================================================================

dominant_source <- function(
  atlas,
  sample_name,
  metric,
  output_name
) {
  z <- atlas[
    atlas$sample == sample_name,
    ,
    drop = FALSE
  ]

  genes <- unique(
    z$human_ligand
  )

  out <- vector(
    "list",
    length(genes)
  )

  for (i in seq_along(genes)) {
    g <- genes[[i]]

    x <- z[
      z$human_ligand == g,
      ,
      drop = FALSE
    ]

    vals <- x[[metric]]

    if (!length(vals) ||
        all(is.na(vals)) ||
        max(vals, na.rm = TRUE) <= 0) {
      grp <- NA_character_
      val <- 0
    } else {
      j <- which.max(vals)
      grp <- x$group[j]
      val <- vals[j]
    }

    out[[i]] <- data.frame(
      human_ligand = g,
      sample = sample_name,
      source_metric = output_name,
      dominant_group = grp,
      dominant_metric_value = val,
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    out
  )
}

dominant_lineage_total <- do.call(
  rbind,
  lapply(
    SAMPLES,
    function(s) {
      dominant_source(
        lineage_atlas,
        s,
        "fraction_of_sample_ligand_UMI",
        "total_UMI_fraction"
      )
    }
  )
)

dominant_lineage_per1000 <- do.call(
  rbind,
  lapply(
    SAMPLES,
    function(s) {
      dominant_source(
        lineage_atlas,
        s,
        "UMI_per_1000_cells",
        "UMI_per_1000_cells"
      )
    }
  )
)

dominant_state_total <- do.call(
  rbind,
  lapply(
    SAMPLES,
    function(s) {
      dominant_source(
        state_atlas,
        s,
        "fraction_of_sample_ligand_UMI",
        "total_UMI_fraction"
      )
    }
  )
)

dominant_sources <- rbind(
  dominant_lineage_total,
  dominant_lineage_per1000,
  dominant_state_total
)

rownames(dominant_sources) <- NULL

write_tsv(
  dominant_sources,
  file.path(
    TAB_OUT,
    paste0(
      "06_dominant_ligand_source_summary_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 10. Strict human -> mouse ligand orthology
# ==============================================================================

msg(strrep("-", 92))
msg(
  "Building babelgene human -> mouse ligand ortholog map..."
)

ortholog_raw <- babelgene::orthologs(
  genes = human_ligands,
  species = "mouse",
  human = TRUE,
  min_support = 3,
  top = FALSE
)

required_ortholog_cols <- c(
  "human_symbol",
  "symbol",
  "support_n"
)

missing_ortholog_cols <- setdiff(
  required_ortholog_cols,
  colnames(ortholog_raw)
)

if (length(missing_ortholog_cols)) {
  stop2(
    "babelgene output missing columns: ",
    paste(
      missing_ortholog_cols,
      collapse = ", "
    )
  )
}

ortholog_raw <- as.data.frame(
  ortholog_raw,
  stringsAsFactors = FALSE
)

# Uniqueness is enforced across the FULL queried human NicheNet ligand universe.
human_n <- table(
  ortholog_raw$human_symbol
)

mouse_n <- table(
  ortholog_raw$symbol
)

ortholog_raw$n_mouse_matches_for_human_ligand <-
  as.integer(
    human_n[
      ortholog_raw$human_symbol
    ]
  )

ortholog_raw$n_human_ligands_for_mouse_symbol <-
  as.integer(
    mouse_n[
      ortholog_raw$symbol
    ]
  )

ortholog_raw$strict_unique_1to1_within_ligand_universe <-
  ortholog_raw$n_mouse_matches_for_human_ligand == 1L &
  ortholog_raw$n_human_ligands_for_mouse_symbol == 1L

write_tsv(
  ortholog_raw,
  file.path(
    TAB_OUT,
    paste0(
      "07_human_mouse_ligand_ortholog_raw_babelgene_",
      VERSION,
      ".tsv"
    )
  )
)

strict_ortholog <- ortholog_raw[
  ortholog_raw$strict_unique_1to1_within_ligand_universe,
  ,
  drop = FALSE
]

strict_ortholog <- strict_ortholog[
  !duplicated(
    strict_ortholog$human_symbol
  ),
  ,
  drop = FALSE
]

strict_map <- data.frame(
  human_ligand = strict_ortholog$human_symbol,
  mouse_ortholog = strict_ortholog$symbol,
  ortholog_support_n = strict_ortholog$support_n,
  ortholog_support = if (
    "support" %in% colnames(strict_ortholog)
  ) {
    strict_ortholog$support
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)

strict_map$n_mouse_LR_receptors <- vapply(
  strict_map$mouse_ortholog,
  function(g) {
    length(
      unique(
        mouse_lr$mouse_receptor[
          mouse_lr$mouse_ligand == g
        ]
      )
    )
  },
  integer(1)
)

strict_map$in_mouse_NicheNet_LR <-
  strict_map$mouse_ortholog %in%
  mouse_lr$mouse_ligand

strict_map$in_mouse_ligand_target_matrix <-
  strict_map$mouse_ortholog %in%
  mouse_lt_ligands

strict_map$cross_species_NicheNet_compatible <-
  strict_map$in_mouse_NicheNet_LR &
  strict_map$in_mouse_ligand_target_matrix

write_tsv(
  strict_map,
  file.path(
    TAB_OUT,
    paste0(
      "08_strict_human_mouse_ligand_resource_compatibility_",
      VERSION,
      ".tsv"
    )
  )
)

msg(
  "Strict unique ligand orthologs: ",
  nrow(strict_map)
)

msg(
  "Strict ortholog + mouse LR + mouse ligand-target compatible: ",
  sum(
    strict_map$cross_species_NicheNet_compatible
  )
)

# ==============================================================================
# 11. Build paired Tx17 / Tx5 ligand expression table
# ==============================================================================

to_wide_sample <- function(
  df,
  sample_name
) {
  z <- df[
    df$sample == sample_name,
    ,
    drop = FALSE
  ]

  keep <- c(
    "human_ligand",
    "feature_measurable",
    "UMI",
    "positive_cells",
    "positive_fraction",
    "UMI_per_1000_cells",
    "pseudobulk_CPM",
    "mean_UMI_per_cell",
    "mean_UMI_per_positive_cell"
  )

  z <- z[
    ,
    keep,
    drop = FALSE
  ]

  names(z)[
    names(z) != "human_ligand"
  ] <- paste0(
    sample_name,
    "_",
    names(z)[
      names(z) != "human_ligand"
    ]
  )

  z
}

paired <- merge(
  to_wide_sample(
    sample_quant,
    "Tx17"
  ),
  to_wide_sample(
    sample_quant,
    "Tx5"
  ),
  by = "human_ligand",
  all = TRUE,
  sort = FALSE
)

paired$shared_Tx_feature_measurable <-
  paired$Tx17_feature_measurable &
  paired$Tx5_feature_measurable

paired$detected_Tx17 <-
  paired$shared_Tx_feature_measurable &
  !is.na(
    paired$Tx17_positive_cells
  ) &
  paired$Tx17_positive_cells > 0

paired$detected_Tx5 <-
  paired$shared_Tx_feature_measurable &
  !is.na(
    paired$Tx5_positive_cells
  ) &
  paired$Tx5_positive_cells > 0

paired$detected_both_Tx <-
  paired$detected_Tx17 &
  paired$detected_Tx5

paired$detected_either_Tx <-
  paired$detected_Tx17 |
  paired$detected_Tx5

paired$at_least_3_positive_cells_each_Tx <-
  paired$detected_both_Tx &
  paired$Tx17_positive_cells >= 3 &
  paired$Tx5_positive_cells >= 3

paired$min_positive_fraction_Tx <-
  ifelse(
    paired$shared_Tx_feature_measurable,
    pmin(
      paired$Tx17_positive_fraction,
      paired$Tx5_positive_fraction,
      na.rm = TRUE
    ),
    NA_real_
  )

paired$mean_positive_fraction_Tx <-
  ifelse(
    paired$shared_Tx_feature_measurable,
    rowMeans(
      cbind(
        paired$Tx17_positive_fraction,
        paired$Tx5_positive_fraction
      ),
      na.rm = TRUE
    ),
    NA_real_
  )

paired$min_pseudobulk_CPM_Tx <-
  ifelse(
    paired$shared_Tx_feature_measurable,
    pmin(
      paired$Tx17_pseudobulk_CPM,
      paired$Tx5_pseudobulk_CPM,
      na.rm = TRUE
    ),
    NA_real_
  )

paired$geometric_mean_pseudobulk_CPM_Tx <-
  ifelse(
    paired$detected_both_Tx,
    sqrt(
      pmax(
        paired$Tx17_pseudobulk_CPM,
        0
      ) *
      pmax(
        paired$Tx5_pseudobulk_CPM,
        0
      )
    ),
    0
  )

paired$replicate_balance_CPM <-
  ifelse(
    paired$detected_both_Tx &
      pmax(
        paired$Tx17_pseudobulk_CPM,
        paired$Tx5_pseudobulk_CPM
      ) > 0,
    pmin(
      paired$Tx17_pseudobulk_CPM,
      paired$Tx5_pseudobulk_CPM
    ) /
    pmax(
      paired$Tx17_pseudobulk_CPM,
      paired$Tx5_pseudobulk_CPM
    ),
    0
  )

paired$expression_reproducibility_class <-
  ifelse(
    !paired$shared_Tx_feature_measurable,
    "Not_shared_feature_measurable",
    ifelse(
      paired$at_least_3_positive_cells_each_Tx,
      "Both_Tx_at_least_3_positive_cells",
      ifelse(
        paired$detected_both_Tx,
        "Both_Tx_detected_sparse",
        ifelse(
          paired$detected_either_Tx,
          "Single_Tx_detected_only",
          "Measured_but_zero_standard_MEX"
        )
      )
    )
  )

paired$legacy_15_panel_audit_only <-
  paired$human_ligand %in% LEGACY_15

paired$is_IL10_audit_only <-
  paired$human_ligand == "IL10"

# Human LR receptor count.
human_receptor_n <- tapply(
  human_lr$human_receptor,
  human_lr$human_ligand,
  function(x) {
    length(
      unique(x)
    )
  }
)

paired$n_human_LR_receptors <-
  as.integer(
    human_receptor_n[
      paired$human_ligand
    ]
  )

# ==============================================================================
# 12. Attach dominant lineage/state source
# ==============================================================================

extract_dom <- function(
  dom_df,
  sample_name,
  metric_name,
  prefix
) {

  z <- dom_df[
    dom_df$sample == sample_name &
      dom_df$source_metric == metric_name,
    c(
      "human_ligand",
      "dominant_group",
      "dominant_metric_value"
    ),
    drop = FALSE
  ]

  if (anyDuplicated(z$human_ligand)) {

    dup <- unique(
      z$human_ligand[
        duplicated(z$human_ligand) |
          duplicated(
            z$human_ligand,
            fromLast = TRUE
          )
      ]
    )

    stop2(
      "Non-unique dominant-source table before merge: ",
      prefix,
      " | duplicated ligands=",
      length(dup),
      " | first=",
      paste(
        head(dup, 10),
        collapse = ","
      )
    )
  }

  names(z)[2:3] <- c(
    paste0(
      prefix,
      "_group"
    ),
    paste0(
      prefix,
      "_value"
    )
  )

  z
}

## IMPORTANT:
## Do not extract total_UMI_fraction from the combined dominant_sources table.
## That combined object contains BOTH lineage-total and state-total rows.
## Joining that object created a 2 x 2 Cartesian expansion per ligand.
##
## Instead, use the already separated source tables directly.

source_tables <- list(

  extract_dom(
    dominant_lineage_total,
    "Tx17",
    "total_UMI_fraction",
    "Tx17_dominant_lineage_total"
  ),

  extract_dom(
    dominant_lineage_total,
    "Tx5",
    "total_UMI_fraction",
    "Tx5_dominant_lineage_total"
  ),

  extract_dom(
    dominant_lineage_per1000,
    "Tx17",
    "UMI_per_1000_cells",
    "Tx17_dominant_lineage_per1000"
  ),

  extract_dom(
    dominant_lineage_per1000,
    "Tx5",
    "UMI_per_1000_cells",
    "Tx5_dominant_lineage_per1000"
  ),

  extract_dom(
    dominant_state_total,
    "Tx17",
    "total_UMI_fraction",
    "Tx17_dominant_state_total"
  ),

  extract_dom(
    dominant_state_total,
    "Tx5",
    "total_UMI_fraction",
    "Tx5_dominant_state_total"
  )
)

paired <- Reduce(
  function(x, y) {

    out <- merge(
      x,
      y,
      by = "human_ligand",
      all.x = TRUE,
      sort = FALSE
    )

    if (anyDuplicated(out$human_ligand)) {
      stop2(
        "Dominant-source merge created duplicated human ligands."
      )
    }

    out
  },
  c(
    list(paired),
    source_tables
  )
)

if (nrow(paired) != length(human_ligands)) {
  stop2(
    "Ligand table row count changed after source merge: ",
    nrow(paired),
    " vs expected ",
    length(human_ligands)
  )
}

if (anyDuplicated(paired$human_ligand)) {
  stop2(
    "Ligand table contains duplicates after source merge."
  )
}

paired$dominant_lineage_total_consensus <-
  !is.na(
    paired$Tx17_dominant_lineage_total_group
  ) &
  !is.na(
    paired$Tx5_dominant_lineage_total_group
  ) &
  paired$Tx17_dominant_lineage_total_group ==
    paired$Tx5_dominant_lineage_total_group

paired$dominant_lineage_per1000_consensus <-
  !is.na(
    paired$Tx17_dominant_lineage_per1000_group
  ) &
  !is.na(
    paired$Tx5_dominant_lineage_per1000_group
  ) &
  paired$Tx17_dominant_lineage_per1000_group ==
    paired$Tx5_dominant_lineage_per1000_group

paired$dominant_state_total_consensus <-
  !is.na(
    paired$Tx17_dominant_state_total_group
  ) &
  !is.na(
    paired$Tx5_dominant_state_total_group
  ) &
  paired$Tx17_dominant_state_total_group ==
    paired$Tx5_dominant_state_total_group

# ==============================================================================
# 13. Attach orthology / mouse NicheNet compatibility
# ==============================================================================

paired <- merge(
  paired,
  strict_map,
  by = "human_ligand",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(paired) != length(human_ligands)) {
  stop2(
    "Ligand table row count changed after orthology merge: ",
    nrow(paired),
    " vs expected ",
    length(human_ligands)
  )
}

if (anyDuplicated(paired$human_ligand)) {
  stop2(
    "Orthology merge created duplicated human ligands."
  )
}

paired$strict_unique_1to1_within_ligand_universe <-
  !is.na(
    paired$mouse_ortholog
  )

paired$cross_species_NicheNet_compatible[
  is.na(
    paired$cross_species_NicheNet_compatible
  )
] <- FALSE

paired$in_mouse_NicheNet_LR[
  is.na(
    paired$in_mouse_NicheNet_LR
  )
] <- FALSE

paired$in_mouse_ligand_target_matrix[
  is.na(
    paired$in_mouse_ligand_target_matrix
  )
] <- FALSE

# Primary admissibility for the NEXT stage.
# Expression ranking itself does NOT require orthology, because this table
# should preserve human-only ligands for audit.
paired$eligible_for_mouse_receiver_bridge <-
  paired$shared_Tx_feature_measurable &
  paired$detected_either_Tx &
  paired$strict_unique_1to1_within_ligand_universe &
  paired$cross_species_NicheNet_compatible

paired$eligible_for_primary_replicate_supported_bridge <-
  paired$detected_both_Tx &
  paired$strict_unique_1to1_within_ligand_universe &
  paired$cross_species_NicheNet_compatible

# No arbitrary composite biological score is used.
# Rank is lexicographic:
#   1) bridge compatibility
#   2) detected in both Tx
#   3) >=3 positive cells in both Tx
#   4) minimum detection fraction
#   5) geometric mean CPM
#   6) replicate balance
ord <- order(
  -as.integer(
    paired$cross_species_NicheNet_compatible
  ),
  -as.integer(
    paired$detected_both_Tx
  ),
  -as.integer(
    paired$at_least_3_positive_cells_each_Tx
  ),
  -ifelse(
    is.na(
      paired$min_positive_fraction_Tx
    ),
    -Inf,
    paired$min_positive_fraction_Tx
  ),
  -ifelse(
    is.na(
      paired$geometric_mean_pseudobulk_CPM_Tx
    ),
    -Inf,
    paired$geometric_mean_pseudobulk_CPM_Tx
  ),
  -ifelse(
    is.na(
      paired$replicate_balance_CPM
    ),
    -Inf,
    paired$replicate_balance_CPM
  ),
  paired$human_ligand
)

paired <- paired[
  ord,
  ,
  drop = FALSE
]

paired$descriptive_sender_rank <- seq_len(
  nrow(paired)
)

# Put key columns first.
first_cols <- c(
  "descriptive_sender_rank",
  "human_ligand",
  "expression_reproducibility_class",
  "shared_Tx_feature_measurable",
  "detected_both_Tx",
  "at_least_3_positive_cells_each_Tx",
  "detected_either_Tx",
  "min_positive_fraction_Tx",
  "mean_positive_fraction_Tx",
  "geometric_mean_pseudobulk_CPM_Tx",
  "min_pseudobulk_CPM_Tx",
  "replicate_balance_CPM",
  "Tx17_UMI",
  "Tx17_positive_cells",
  "Tx17_positive_fraction",
  "Tx17_UMI_per_1000_cells",
  "Tx17_pseudobulk_CPM",
  "Tx5_UMI",
  "Tx5_positive_cells",
  "Tx5_positive_fraction",
  "Tx5_UMI_per_1000_cells",
  "Tx5_pseudobulk_CPM",
  "Tx17_dominant_lineage_total_group",
  "Tx5_dominant_lineage_total_group",
  "dominant_lineage_total_consensus",
  "Tx17_dominant_lineage_per1000_group",
  "Tx5_dominant_lineage_per1000_group",
  "dominant_lineage_per1000_consensus",
  "Tx17_dominant_state_total_group",
  "Tx5_dominant_state_total_group",
  "dominant_state_total_consensus",
  "mouse_ortholog",
  "ortholog_support_n",
  "strict_unique_1to1_within_ligand_universe",
  "n_human_LR_receptors",
  "n_mouse_LR_receptors",
  "in_mouse_NicheNet_LR",
  "in_mouse_ligand_target_matrix",
  "cross_species_NicheNet_compatible",
  "eligible_for_mouse_receiver_bridge",
  "eligible_for_primary_replicate_supported_bridge",
  "legacy_15_panel_audit_only",
  "is_IL10_audit_only"
)

first_cols <- first_cols[
  first_cols %in% colnames(paired)
]

other_cols <- setdiff(
  colnames(paired),
  first_cols
)

paired <- paired[
  ,
  c(
    first_cols,
    other_cols
  ),
  drop = FALSE
]

write_tsv(
  paired,
  file.path(
    TAB_OUT,
    paste0(
      "09_human_graft_all_ligandome_ranked_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 14. Primary bridge candidate table + audit-only legacy table
# ==============================================================================

bridge_candidates <- paired[
  paired$eligible_for_mouse_receiver_bridge,
  ,
  drop = FALSE
]

write_tsv(
  bridge_candidates,
  file.path(
    TAB_OUT,
    paste0(
      "10_cross_species_bridge_candidate_ligands_",
      VERSION,
      ".tsv"
    )
  )
)

legacy_audit <- paired[
  paired$human_ligand %in% LEGACY_15,
  ,
  drop = FALSE
]

legacy_audit <- legacy_audit[
  match(
    LEGACY_15,
    legacy_audit$human_ligand
  ),
  ,
  drop = FALSE
]

legacy_audit <- legacy_audit[
  !is.na(
    legacy_audit$human_ligand
  ),
  ,
  drop = FALSE
]

write_tsv(
  legacy_audit,
  file.path(
    TAB_OUT,
    paste0(
      "11_legacy15_IL10_audit_only_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 15. Compact summaries
# ==============================================================================

summary_table <- data.frame(
  metric = c(
    "NicheNet_human_LR_ligands",
    "shared_Tx_feature_measurable",
    "detected_either_Tx",
    "detected_both_Tx",
    "at_least_3_positive_cells_each_Tx",
    "strict_unique_human_mouse_ligand_ortholog",
    "cross_species_NicheNet_compatible",
    "eligible_for_mouse_receiver_bridge",
    "eligible_for_primary_replicate_supported_bridge"
  ),
  n = c(
    length(human_ligands),
    sum(
      paired$shared_Tx_feature_measurable,
      na.rm = TRUE
    ),
    sum(
      paired$detected_either_Tx,
      na.rm = TRUE
    ),
    sum(
      paired$detected_both_Tx,
      na.rm = TRUE
    ),
    sum(
      paired$at_least_3_positive_cells_each_Tx,
      na.rm = TRUE
    ),
    sum(
      paired$strict_unique_1to1_within_ligand_universe,
      na.rm = TRUE
    ),
    sum(
      paired$cross_species_NicheNet_compatible,
      na.rm = TRUE
    ),
    sum(
      paired$eligible_for_mouse_receiver_bridge,
      na.rm = TRUE
    ),
    sum(
      paired$eligible_for_primary_replicate_supported_bridge,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_table,
  file.path(
    TAB_OUT,
    paste0(
      "12_all_ligandome_summary_",
      VERSION,
      ".tsv"
    )
  )
)

# ==============================================================================
# 16. Simple expression figure for top bridge-compatible ligands
# ==============================================================================

top_n <- min(
  60L,
  nrow(
    bridge_candidates
  )
)

if (top_n > 0L) {
  top_genes <- bridge_candidates$human_ligand[
    seq_len(top_n)
  ]

  fig_dat <- sample_quant[
    sample_quant$human_ligand %in% top_genes,
    ,
    drop = FALSE
  ]

  fig_dat$human_ligand <- factor(
    fig_dat$human_ligand,
    levels = rev(
      top_genes
    )
  )

  fig_dat$sample <- factor(
    fig_dat$sample,
    levels = SAMPLES
  )

  p <- ggplot(
    fig_dat,
    aes(
      x = sample,
      y = human_ligand
    )
  ) +
    geom_point(
      aes(
        size = positive_fraction,
        fill = log1p(
          pseudobulk_CPM
        )
      ),
      shape = 21
    ) +
    scale_size_continuous(
      name = "Positive fraction",
      range = c(
        0.5,
        5
      )
    ) +
    scale_fill_gradient(
      name = "log1p CPM"
    ) +
    labs(
      title = "Human graft ligandome: top cross-species-compatible ligands",
      subtitle = "Descriptive Tx17/Tx5 expression; no IL10 weighting",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.y = element_text(
        size = 7
      )
    )

  ggsave(
    filename = file.path(
      FIG_OUT,
      paste0(
        "01_top_bridge_compatible_ligands_",
        VERSION,
        ".pdf"
      )
    ),
    plot = p,
    width = 8.2,
    height = 12,
    units = "in",
    limitsize = FALSE
  )
}

# ==============================================================================
# 17. Save compact object
# ==============================================================================

saveRDS(
  list(
    version = VERSION,
    human_lr = human_lr,
    mouse_lr = mouse_lr,
    feature_space = feature_space,
    sample_quant = sample_quant,
    lineage_atlas = lineage_atlas,
    state_atlas = state_atlas,
    dominant_sources = dominant_sources,
    ortholog_raw = ortholog_raw,
    strict_map = strict_map,
    ranked_ligandome = paired,
    bridge_candidates = bridge_candidates,
    summary = summary_table
  ),
  file.path(
    OBJ_OUT,
    paste0(
      "Human_Graft_AllLigandome_",
      VERSION,
      ".rds"
    )
  ),
  compress = "xz"
)

# ==============================================================================
# 18. Manifest / session info
# ==============================================================================

manifest <- data.frame(
  item = c(
    "version",
    "analysis_principle",
    "human_sender_samples",
    "context_metadata",
    "Tx17_MEX",
    "Tx5_MEX",
    "expected_Tx17_cells",
    "expected_Tx5_cells",
    "expected_Tx17_features",
    "expected_Tx5_features",
    "ligand_universe",
    "NicheNet_resource_record",
    "orthology_method",
    "orthology_min_support",
    "human_mouse_UMI_direct_comparison",
    "IL10_priority",
    "legacy_15_priority",
    "standard_MEX_boundary",
    "next_step"
  ),
  value = c(
    VERSION,
    "ligand-agnostic human graft sender atlas",
    "Tx17;Tx5",
    CONTEXT_META,
    MEX_DIRS[["Tx17"]],
    MEX_DIRS[["Tx5"]],
    EXPECTED_ANALYSIS_CELLS[["Tx17"]],
    EXPECTED_ANALYSIS_CELLS[["Tx5"]],
    EXPECTED_FEATURES[["Tx17"]],
    EXPECTED_FEATURES[["Tx5"]],
    "all unique ligand genes in NicheNet-v2 human LR network",
    "Zenodo 7074291",
    "babelgene human->mouse; top=FALSE; unique human and mouse mapping within full human ligand universe",
    ">=3 supporting databases",
    "PROHIBITED",
    "NONE; audit annotation only",
    "NONE; audit annotation only",
    "Discovery only. A zero standard-MEX count is not proof of biological absence before probe-status audit.",
    "v7.2.1 mouse receiver receptor competence + LR admissibility; then v7.2.2 NicheNet CORE/EXPANDED target activity"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  manifest,
  file.path(
    TAB_OUT,
    paste0(
      "00_run_manifest_",
      VERSION,
      ".tsv"
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

msg(strrep("=", 92))
msg("STATUS: SUCCESS")
msg("Output: ", OUT)
msg("Key review files:")
msg("  01_Tx_MEX_context_mapping_QC_", VERSION, ".tsv")
msg("  02_human_NicheNet_ligand_feature_space_", VERSION, ".tsv")
msg("  09_human_graft_all_ligandome_ranked_", VERSION, ".tsv")
msg("  10_cross_species_bridge_candidate_ligands_", VERSION, ".tsv")
msg("  11_legacy15_IL10_audit_only_", VERSION, ".tsv")
msg("  12_all_ligandome_summary_", VERSION, ".tsv")
msg("Next: receptor competence in each v7.1.2 mouse receiver; then NicheNet against CORE/EXPANDED host targets.")
msg(strrep("=", 92))
