#!/usr/bin/env Rscript

# ==============================================================================
# WholeLiver_DecontX_Input_Preflight_v1.0.R
#
# PURPOSE
# -------
# Inventory candidate whole-liver Seurat objects BEFORE any DecontX analysis.
# No decontamination is performed. No source object is modified.
#
# Outputs:
#   - object inventory
#   - likely sample / broad-celltype metadata columns
#   - sample x broad-celltype counts
#   - overlap of frozen Mphi barcodes with whole-liver candidate
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_CANDIDATES <- c(
  file.path(
    ROOT, "Human_MASH_RDS",
    "Human_MASH_LT1to8_Canonical_v1.3.3",
    "objects", "Human_MASH_LT1to8_Canonical_v1.3.3.rds"
  ),
  file.path(
    ROOT, "Human_MASH_RDS",
    "Human_MASH_LT1to8_Canonical_v1.3",
    "objects", "Human_MASH_LT1to8_Canonical_v1.3.rds"
  )
)

HUMAN_MPHI <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects", "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

mouse_hits <- list.files(
  ROOT,
  pattern = "^Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

MOUSE_MPHI <- file.path(
  ROOT, "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS", "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTDIR <- file.path(
  ROOT, "CrossSpecies_MASH_Mphi",
  "WholeLiver_DecontX_Input_Preflight_v1.0"
)
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(LOGDIR, "run.log")
cat("", file = LOG)

msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = LOG, append = TRUE)
}

write_tsv <- function(x, f) {
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

candidate_columns <- function(meta, patterns) {
  nm <- colnames(meta)
  unique(unlist(lapply(patterns, function(p) grep(p, nm, value = TRUE, ignore.case = TRUE))))
}

inventory_one <- function(path, species, role) {
  if (!file.exists(path)) {
    return(list(
      summary = data.frame(
        species = species, role = role, path = path, exists = FALSE,
        n_cells = NA_integer_, n_features = NA_integer_,
        assays = NA_character_, layers_RNA = NA_character_,
        reductions = NA_character_, graphs = NA_character_,
        stringsAsFactors = FALSE
      ),
      object = NULL
    ))
  }

  msg("Loading: ", path)
  obj <- readRDS(path)

  layers_rna <- if ("RNA" %in% Assays(obj)) {
    paste(Layers(obj[["RNA"]]), collapse = ";")
  } else {
    NA_character_
  }

  list(
    summary = data.frame(
      species = species,
      role = role,
      path = path,
      exists = TRUE,
      n_cells = ncol(obj),
      n_features = nrow(obj),
      assays = paste(Assays(obj), collapse = ";"),
      layers_RNA = layers_rna,
      reductions = paste(Reductions(obj), collapse = ";"),
      graphs = paste(names(obj@graphs), collapse = ";"),
      stringsAsFactors = FALSE
    ),
    object = obj
  )
}

human_path <- HUMAN_CANDIDATES[file.exists(HUMAN_CANDIDATES)][1]
if (length(human_path) == 0 || is.na(human_path)) human_path <- HUMAN_CANDIDATES[1]

mouse_path <- if (length(mouse_hits)) mouse_hits[1] else NA_character_

items <- list(
  inventory_one(human_path, "Human", "whole_liver_candidate"),
  inventory_one(HUMAN_MPHI, "Human", "frozen_Mphi"),
  if (!is.na(mouse_path)) inventory_one(mouse_path, "Mouse", "whole_liver_candidate")
  else list(
    summary = data.frame(
      species="Mouse", role="whole_liver_candidate", path=NA_character_,
      exists=FALSE, n_cells=NA_integer_, n_features=NA_integer_,
      assays=NA_character_, layers_RNA=NA_character_,
      reductions=NA_character_, graphs=NA_character_,
      stringsAsFactors=FALSE
    ),
    object=NULL
  ),
  inventory_one(MOUSE_MPHI, "Mouse", "frozen_Mphi")
)

inv <- do.call(rbind, lapply(items, `[[`, "summary"))
write_tsv(inv, file.path(TABDIR, "01_object_inventory.tsv"))

human_whole <- items[[1]]$object
human_mphi <- items[[2]]$object
mouse_whole <- items[[3]]$object
mouse_mphi <- items[[4]]$object

audit_meta <- function(obj, species, role) {
  if (is.null(obj)) return(data.frame())

  sample_cols <- candidate_columns(
    obj@meta.data,
    c("^sample", "orig\\.ident", "donor", "library", "condition", "group")
  )

  celltype_cols <- candidate_columns(
    obj@meta.data,
    c("celltype", "cell_type", "annotation", "lineage", "broad", "layer1", "identity")
  )

  data.frame(
    species = species,
    role = role,
    category = c(rep("sample_candidate", length(sample_cols)),
                 rep("celltype_candidate", length(celltype_cols))),
    column = c(sample_cols, celltype_cols),
    n_unique = c(
      vapply(sample_cols, function(x) length(unique(obj@meta.data[[x]])), integer(1)),
      vapply(celltype_cols, function(x) length(unique(obj@meta.data[[x]])), integer(1))
    ),
    example_values = c(
      vapply(sample_cols, function(x) paste(head(unique(as.character(obj@meta.data[[x]])), 20), collapse=";"), character(1)),
      vapply(celltype_cols, function(x) paste(head(unique(as.character(obj@meta.data[[x]])), 20), collapse=";"), character(1))
    ),
    stringsAsFactors = FALSE
  )
}

meta_audit <- rbind(
  audit_meta(human_whole, "Human", "whole_liver_candidate"),
  audit_meta(mouse_whole, "Mouse", "whole_liver_candidate")
)
write_tsv(meta_audit, file.path(TABDIR, "02_metadata_column_candidates.tsv"))

overlap_one <- function(whole, mphi, species) {
  if (is.null(whole) || is.null(mphi)) {
    return(data.frame(
      species=species, whole_cells=NA_integer_, mphi_cells=NA_integer_,
      overlap=NA_integer_, fraction_mphi_in_whole=NA_real_,
      stringsAsFactors=FALSE
    ))
  }
  ov <- intersect(colnames(whole), colnames(mphi))
  data.frame(
    species=species,
    whole_cells=ncol(whole),
    mphi_cells=ncol(mphi),
    overlap=length(ov),
    fraction_mphi_in_whole=length(ov)/ncol(mphi),
    stringsAsFactors=FALSE
  )
}

overlap <- rbind(
  overlap_one(human_whole, human_mphi, "Human"),
  overlap_one(mouse_whole, mouse_mphi, "Mouse")
)
write_tsv(overlap, file.path(TABDIR, "03_Mphi_barcode_overlap_with_whole_liver.tsv"))

# Produce frequency tables for all plausible sample / celltype pairs.
pair_tables <- list()
k <- 0L

emit_pairs <- function(obj, species) {
  if (is.null(obj)) return(list())

  sample_cols <- candidate_columns(
    obj@meta.data,
    c("^sample", "orig\\.ident", "donor", "library", "condition", "group")
  )

  celltype_cols <- candidate_columns(
    obj@meta.data,
    c("celltype", "cell_type", "annotation", "lineage", "broad", "layer1", "identity")
  )

  out <- list()
  j <- 0L

  for (scol in sample_cols) {
    for (ccol in celltype_cols) {
      tb <- as.data.frame(
        table(
          sample = as.character(obj@meta.data[[scol]]),
          celltype = as.character(obj@meta.data[[ccol]])
        ),
        stringsAsFactors = FALSE
      )
      tb <- tb[tb$Freq > 0, , drop=FALSE]
      tb$species <- species
      tb$sample_col <- scol
      tb$celltype_col <- ccol
      j <- j + 1L
      out[[j]] <- tb[, c("species","sample_col","celltype_col","sample","celltype","Freq")]
    }
  }
  out
}

pair_tables <- c(
  emit_pairs(human_whole, "Human"),
  emit_pairs(mouse_whole, "Mouse")
)

if (length(pair_tables)) {
  write_tsv(
    do.call(rbind, pair_tables),
    file.path(TABDIR, "04_sample_x_celltype_candidate_counts.tsv")
  )
}

msg("STATUS: SUCCESS")
msg("Review:")
msg("  01_object_inventory.tsv")
msg("  02_metadata_column_candidates.tsv")
msg("  03_Mphi_barcode_overlap_with_whole_liver.tsv")
msg("  04_sample_x_celltype_candidate_counts.tsv")
