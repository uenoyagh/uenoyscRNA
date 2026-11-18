#!/usr/bin/env Rscript

# ==============================================================================
# Mouse_HSC_WholeLiver_DecontX_Input_Preflight_v1.0.R
#
# PURPOSE
# -------
# Identify and validate the correct whole-liver filtered-count source object
# for sample-wise DecontX sensitivity analysis of the refined Mouse HSC object.
#
# This script:
#   - does NOT run DecontX
#   - does NOT modify source objects
#   - searches known whole-liver RDS filenames under ROOT
#   - checks barcode overlap with the HSC object
#   - detects sample columns
#   - checks sample agreement for overlapping HSC barcodes
#   - inventories RNA count layers
#   - reports the best-matching whole-liver source
#
# POLICY
# ------
# Original filtered counts remain primary.
# Whole-liver filtered-count DecontX will be sensitivity analysis only.
# Raw/empty-droplet background is not assumed available.
#
# Version: v1.0
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "HSC_RefinedAnnotation_v6.1.2",
  "RDS",
  "Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds"
)

OUTDIR <- file.path(
  ROOT,
  "CrossSpecies_MASH_HSC",
  "Mouse_HSC_WholeLiver_DecontX_Input_Preflight_v1.0"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, f) {
  write.table(
    x, f,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

if (!file.exists(HSC_RDS)) stop("Missing HSC RDS: ", HSC_RDS)

message("Loading HSC object...")
hsc <- readRDS(HSC_RDS)
hmd <- hsc@meta.data

HSC_SAMPLE_COL <- if ("sample_hsc_v610" %in% colnames(hmd)) {
  "sample_hsc_v610"
} else if ("sample" %in% colnames(hmd)) {
  "sample"
} else {
  stop("No HSC sample column found.")
}

hsc_sample <- as.character(hmd[[HSC_SAMPLE_COL]])

KNOWN_FILENAMES <- c(
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds",
  "RH260519343_mo4sample_clustered_annotated.rds"
)

candidates <- unique(unlist(
  lapply(KNOWN_FILENAMES, function(fn) {
    list.files(
      ROOT,
      pattern = paste0("^", gsub("\\.", "\\\\.", fn), "$"),
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = FALSE
    )
  })
))

# Also include broad whole-liver candidates if exact filenames are absent.
if (!length(candidates)) {
  candidates <- unique(c(
    list.files(
      ROOT,
      pattern = "Mouse.*(whole|object|clustered|annotated).*\\.rds$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  ))
}

if (!length(candidates)) {
  stop("No candidate whole-liver RDS found under: ", ROOT)
}

write_tsv(
  data.frame(
    candidate_path = candidates,
    size_bytes = file.info(candidates)$size,
    stringsAsFactors = FALSE
  ),
  file.path(OUTDIR, "01_candidate_whole_liver_RDS.tsv")
)

detect_sample_col <- function(md) {
  preferred <- c(
    "sample_display_FIXED2",
    "sample",
    "orig.ident",
    "sample_display",
    "sample_hsc_v610"
  )

  for (cn in preferred) {
    if (cn %in% colnames(md)) {
      vals <- unique(as.character(md[[cn]]))
      if (any(c("Sham1", "Sham20", "Tx17", "Tx5") %in% vals)) {
        return(cn)
      }
    }
  }

  for (cn in colnames(md)) {
    vals <- unique(as.character(md[[cn]]))
    if (any(c("Sham1", "Sham20", "Tx17", "Tx5") %in% vals)) {
      return(cn)
    }
  }

  NA_character_
}

candidate_rows <- list()
sample_match_rows <- list()
layer_rows <- list()

k <- 0L
s <- 0L
l <- 0L

for (path in candidates) {
  message("Inspecting candidate: ", path)

  obj <- tryCatch(
    readRDS(path),
    error = function(e) e
  )

  if (inherits(obj, "error")) {
    k <- k + 1L
    candidate_rows[[k]] <- data.frame(
      candidate_path = path,
      load_ok = FALSE,
      n_cells = NA_integer_,
      n_features = NA_integer_,
      sample_col = NA_character_,
      overlap_n = NA_integer_,
      overlap_fraction_of_HSC = NA_real_,
      sample_match_fraction = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  md <- obj@meta.data
  sc <- detect_sample_col(md)

  overlap <- intersect(colnames(hsc), colnames(obj))
  overlap_fraction <- length(overlap) / ncol(hsc)

  sample_match_fraction <- NA_real_

  if (!is.na(sc) && length(overlap)) {
    hidx <- match(overlap, colnames(hsc))
    oidx <- match(overlap, colnames(obj))

    hs <- hsc_sample[hidx]
    os <- as.character(md[[sc]])[oidx]

    sample_match_fraction <- mean(hs == os, na.rm = TRUE)

    stab <- as.data.frame(
      table(
        HSC_sample = hs,
        whole_liver_sample = os
      ),
      stringsAsFactors = FALSE
    )
    stab <- stab[stab$Freq > 0, , drop = FALSE]
    stab$candidate_path <- path

    s <- s + 1L
    sample_match_rows[[s]] <- stab[
      ,
      c(
        "candidate_path",
        "HSC_sample",
        "whole_liver_sample",
        "Freq"
      )
    ]
  }

  assays <- Assays(obj)
  for (a in assays) {
    lyr <- tryCatch(
      Layers(obj[[a]]),
      error = function(e) character(0)
    )
    if (!length(lyr)) lyr <- NA_character_

    for (ly in lyr) {
      l <- l + 1L
      layer_rows[[l]] <- data.frame(
        candidate_path = path,
        assay = a,
        layer = ly,
        stringsAsFactors = FALSE
      )
    }
  }

  k <- k + 1L
  candidate_rows[[k]] <- data.frame(
    candidate_path = path,
    load_ok = TRUE,
    n_cells = ncol(obj),
    n_features = nrow(obj),
    sample_col = sc,
    overlap_n = length(overlap),
    overlap_fraction_of_HSC = overlap_fraction,
    sample_match_fraction = sample_match_fraction,
    stringsAsFactors = FALSE
  )

  rm(obj)
  invisible(gc())
}

candidate_summary <- do.call(rbind, candidate_rows)

candidate_summary$score <- with(
  candidate_summary,
  ifelse(
    load_ok,
    overlap_fraction_of_HSC *
      ifelse(is.na(sample_match_fraction), 0, sample_match_fraction),
    -Inf
  )
)

candidate_summary <- candidate_summary[
  order(-candidate_summary$score),
  ,
  drop = FALSE
]

write_tsv(
  candidate_summary,
  file.path(OUTDIR, "02_candidate_match_summary.tsv")
)

sample_match <- if (length(sample_match_rows)) {
  do.call(rbind, sample_match_rows)
} else {
  data.frame(
    candidate_path = character(0),
    HSC_sample = character(0),
    whole_liver_sample = character(0),
    Freq = integer(0)
  )
}

write_tsv(
  sample_match,
  file.path(OUTDIR, "03_overlap_sample_match.tsv")
)

layers <- if (length(layer_rows)) {
  do.call(rbind, layer_rows)
} else {
  data.frame(
    candidate_path = character(0),
    assay = character(0),
    layer = character(0)
  )
}

write_tsv(
  layers,
  file.path(OUTDIR, "04_candidate_assay_layers.tsv")
)

best <- candidate_summary[1, , drop = FALSE]

best_rna_layers <- layers[
  layers$candidate_path == best$candidate_path &
    layers$assay == "RNA",
  ,
  drop = FALSE
]

readiness <- data.frame(
  item = c(
    "HSC_RDS",
    "HSC_cells",
    "HSC_sample_col",
    "candidate_count",
    "best_candidate",
    "best_candidate_cells",
    "HSC_barcode_overlap_n",
    "HSC_barcode_overlap_fraction",
    "sample_match_fraction",
    "best_candidate_has_RNA_assay",
    "best_candidate_has_RNA_counts",
    "all_HSC_barcodes_overlap",
    "all_overlap_samples_match",
    "recommended_next_step"
  ),
  value = c(
    HSC_RDS,
    ncol(hsc),
    HSC_SAMPLE_COL,
    nrow(candidate_summary),
    best$candidate_path,
    best$n_cells,
    best$overlap_n,
    best$overlap_fraction_of_HSC,
    best$sample_match_fraction,
    as.character(any(best_rna_layers$assay == "RNA")),
    as.character(any(best_rna_layers$layer == "counts")),
    as.character(isTRUE(all.equal(best$overlap_fraction_of_HSC, 1))),
    as.character(isTRUE(all.equal(best$sample_match_fraction, 1))),
    "If overlap/sample mapping/counts are valid: run sample-wise whole-liver filtered-count DecontX sensitivity"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  readiness,
  file.path(OUTDIR, "05_analysis_readiness_summary.tsv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt"),
  useBytes = TRUE
)

message("============================================================")
message("STATUS: SUCCESS")
message("Best whole-liver candidate:")
message("  ", best$candidate_path)
message("HSC overlap: ", best$overlap_n, " / ", ncol(hsc))
message("Overlap fraction: ", best$overlap_fraction_of_HSC)
message("Sample match fraction: ", best$sample_match_fraction)
message("Review:")
message("  02_candidate_match_summary.tsv")
message("  03_overlap_sample_match.tsv")
message("  04_candidate_assay_layers.tsv")
message("  05_analysis_readiness_summary.tsv")
message("============================================================")
