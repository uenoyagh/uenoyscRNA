#!/usr/bin/env Rscript

# Human_MASH_RDS_inventory_metadata_audit_v1.0.1.R
# READ-ONLY audit only.
# v1.0.1 fix: safely combines per-object audit tables even when clinical metadata columns differ.
# NO normalization, integration, clustering, UMAP, DE, pseudobulk,
# pathway analysis, cell-cell interaction analysis, annotation change,
# saveRDS(), or RDS overwrite.

options(stringsAsFactors = FALSE, width = 200)

VERSION <- "v1.0.1"
SSD_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
OUTDIR <- file.path(SSD_ROOT, "Human_MASH_RDS",
                    paste0("Inventory_Metadata_Audit_", VERSION))

EXPLICIT_CANDIDATES <- c(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/rds/Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples.rds",
  "/Users/uenoya/Desktop/clustered2.0_20251003.rds",
  "/Users/uenoya/scRNA_MASLD_MASH_pseudobulk/base/RH24029743/Result_20251003_set1_res1.1_2.0/basic/Clustered.rds/clustered2.0_20251003.rds",
  "/Users/uenoya/scRNA_MASLD_MASH_pseudobulk/base/●20251006_scRNASeq/Result_20251003_set1/basic/clustered.rds"
)

SEARCH_DIRS <- c(
  file.path(SSD_ROOT, "Human_MASH_RDS"),
  file.path(SSD_ROOT, "Human_MASH_Mphi_RDS"),
  file.path(SSD_ROOT, "rds"),
  "/Users/uenoya/scRNA_MASLD_MASH_pseudobulk",
  "/Users/uenoya/Desktop"
)

# Recursively scan only the human-specific/internal project trees.
# Keep the generic SSD /rds and Desktop scans shallow.
RECURSIVE_SEARCH_DIRS <- c(
  file.path(SSD_ROOT, "Human_MASH_RDS"),
  file.path(SSD_ROOT, "Human_MASH_Mphi_RDS"),
  "/Users/uenoya/scRNA_MASLD_MASH_pseudobulk"
)

DISCOVERY_PATTERN <- paste(
  c("Human", "LT1", "LT2", "LT3", "LT4", "LT5", "LT6", "LT7", "LT8",
    "SfLB", "clustered2\\.0_20251003", "MASH", "MASLD"),
  collapse = "|"
)

EXPECTED_LABELS <- c(
  paste0("LT", 1:8), "SfLB1", "SfLB2", "SfLB_1", "SfLB_2", "Tx17", "Tx5"
)

SAMPLE_COL_EXACT_PRIORITY <- c(
  "sample_id", "sampleID", "sample", "Sample", "sample_name",
  "orig.ident", "patient_id", "patient", "Patient",
  "donor_id", "donor", "Donor", "subject_id", "subject",
  "specimen_id", "specimen"
)

SAMPLE_COL_REGEX <- "(^|[._])(sample|patient|donor|subject|specimen)([._]|$)|orig\\.ident"
CLINICAL_COL_REGEX <- "(condition|group|diagnos|disease|NAS|fibrosis|fibrotic|stage|steatos|balloon|inflam|sex|gender|age|BMI|body.mass)"
ANNOTATION_COL_REGEX <- "(celltype|cell_type|annotation|annot|identity|ident|cluster)"

if (!dir.exists(SSD_ROOT)) {
  stop("External SSD project root is not mounted or accessible:\n", SSD_ROOT,
       "\nAudit aborted; nothing will be written elsewhere.")
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("=== Human MASH RDS inventory + metadata audit ", VERSION, " ===")
message("READ-ONLY audit. Output: ", OUTDIR)

safe_file_info <- function(paths) {
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  fi <- file.info(paths)
  data.frame(
    path = paths,
    exists = file.exists(paths),
    size_bytes = ifelse(file.exists(paths), fi$size, NA_real_),
    size_GB = ifelse(file.exists(paths), round(fi$size / 1024^3, 3), NA_real_),
    modified = ifelse(file.exists(paths), as.character(fi$mtime), NA_character_),
    stringsAsFactors = FALSE
  )
}

discover_rds <- function() {
  out <- character(0)
  for (d in SEARCH_DIRS) {
    if (!dir.exists(d)) next
    ff <- tryCatch(
      list.files(d, pattern = "\\.rds$", full.names = TRUE,
                 recursive = d %in% RECURSIVE_SEARCH_DIRS,
                 ignore.case = TRUE),
      error = function(e) character(0)
    )
    if (length(ff)) {
      keep <- grepl(DISCOVERY_PATTERN, basename(ff),
                    ignore.case = TRUE, perl = TRUE)
      out <- c(out, ff[keep])
    }
  }
  unique(c(EXPLICIT_CANDIDATES, out))
}

collapse_unique <- function(x, max_n = 20L) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  u <- unique(x)
  if (!length(u)) return(NA_character_)
  if (length(u) > max_n) {
    paste0(paste(utils::head(u, max_n), collapse = " | "),
           " | ... [", length(u), " unique]")
  } else {
    paste(u, collapse = " | ")
  }
}

meta_col_summary <- function(md, object_name, object_path) {
  if (!ncol(md)) return(data.frame())
  do.call(rbind, lapply(names(md), function(nm) {
    x <- md[[nm]]
    data.frame(
      object = object_name,
      path = object_path,
      column = nm,
      class = paste(class(x), collapse = "/"),
      n_cells = length(x),
      n_missing = sum(is.na(x)),
      pct_missing = round(100 * mean(is.na(x)), 3),
      n_unique_nonmissing = length(unique(x[!is.na(x)])),
      example_values = collapse_unique(x, 5L),
      sample_candidate = grepl(SAMPLE_COL_REGEX, nm, ignore.case = TRUE, perl = TRUE),
      clinical_candidate = grepl(CLINICAL_COL_REGEX, nm, ignore.case = TRUE, perl = TRUE),
      annotation_candidate = grepl(ANNOTATION_COL_REGEX, nm, ignore.case = TRUE, perl = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

pick_sample_col <- function(md) {
  nms <- names(md)
  hit <- SAMPLE_COL_EXACT_PRIORITY[SAMPLE_COL_EXACT_PRIORITY %in% nms]
  if (length(hit)) return(hit[[1]])
  cand <- nms[grepl(SAMPLE_COL_REGEX, nms, ignore.case = TRUE, perl = TRUE)]
  if (!length(cand)) return(NA_character_)
  score <- vapply(cand, function(nm) {
    u <- length(unique(md[[nm]][!is.na(md[[nm]])]))
    if (u >= 2 && u <= 100) 0 else abs(u - 20)
  }, numeric(1))
  cand[[which.min(score)]]
}

make_sample_counts <- function(md, sample_col, object_name, object_path) {
  if (is.na(sample_col) || !nzchar(sample_col)) return(data.frame())
  x <- as.character(md[[sample_col]])
  x[is.na(x) | !nzchar(x)] <- "<NA/blank>"
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(
    object = object_name, path = object_path,
    sample_column = sample_col, sample_value = names(tab),
    n_cells = as.integer(tab), stringsAsFactors = FALSE
  )
}

make_crosscheck <- function(md, sample_col, object_name, object_path) {
  if (is.na(sample_col) || !nzchar(sample_col)) return(data.frame())
  clinical_cols <- names(md)[
    grepl(CLINICAL_COL_REGEX, names(md), ignore.case = TRUE, perl = TRUE)
  ]
  clinical_cols <- setdiff(clinical_cols, sample_col)
  sample_vec <- as.character(md[[sample_col]])
  sample_vec[is.na(sample_vec) | !nzchar(sample_vec)] <- "<NA/blank>"
  samples <- unique(sample_vec)

  do.call(rbind, lapply(samples, function(s) {
    idx <- sample_vec == s
    row <- data.frame(
      object = object_name, path = object_path,
      sample_column = sample_col, sample_value = s,
      n_cells = sum(idx), stringsAsFactors = FALSE
    )
    for (cc in clinical_cols) {
      row[[cc]] <- collapse_unique(md[[cc]][idx], 10L)
    }
    row
  }))
}

audit_expected_labels <- function(md, sample_col, cell_names,
                                  object_name, object_path) {
  sample_values <- character(0)
  if (!is.na(sample_col) && nzchar(sample_col)) {
    sample_values <- unique(as.character(md[[sample_col]]))
    sample_values <- sample_values[!is.na(sample_values)]
  }

  do.call(rbind, lapply(EXPECTED_LABELS, function(lbl) {
    starts <- startsWith(cell_names, paste0(lbl, "_")) |
      startsWith(cell_names, paste0(lbl, "-")) |
      startsWith(cell_names, lbl)
    data.frame(
      object = object_name,
      path = object_path,
      sample_column = ifelse(is.na(sample_col), NA_character_, sample_col),
      expected_label = lbl,
      exact_metadata_match = lbl %in% sample_values,
      normalized_metadata_match = gsub("_", "", lbl) %in%
        gsub("_", "", sample_values),
      cellname_prefix_matches = sum(starts),
      stringsAsFactors = FALSE
    )
  }))
}

safe_layers <- function(obj, assay_name) {
  tryCatch({
    if (!requireNamespace("SeuratObject", quietly = TRUE)) return(NA_character_)
    paste(SeuratObject::Layers(obj[[assay_name]]), collapse = " | ")
  }, error = function(e) paste0("<ERROR: ", conditionMessage(e), ">"))
}

audit_structure <- function(obj, object_name, object_path) {
  assays <- tryCatch(names(obj@assays), error = function(e) character(0))
  reductions <- tryCatch(names(obj@reductions), error = function(e) character(0))
  default_assay <- tryCatch({
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
      SeuratObject::DefaultAssay(obj)
    } else NA_character_
  }, error = function(e) NA_character_)

  rows <- list()
  for (a in assays) {
    rows[[length(rows) + 1L]] <- data.frame(
      object = object_name, path = object_path,
      type = "assay", name = a, detail = safe_layers(obj, a),
      stringsAsFactors = FALSE
    )
  }
  for (r in reductions) {
    dims <- tryCatch({
      if (requireNamespace("SeuratObject", quietly = TRUE)) {
        emb <- SeuratObject::Embeddings(obj[[r]])
        paste0(nrow(emb), " cells x ", ncol(emb), " dims")
      } else NA_character_
    }, error = function(e) paste0("<ERROR: ", conditionMessage(e), ">"))
    rows[[length(rows) + 1L]] <- data.frame(
      object = object_name, path = object_path,
      type = "reduction", name = r, detail = dims,
      stringsAsFactors = FALSE
    )
  }
  list(
    default_assay = default_assay,
    assays = assays,
    reductions = reductions,
    table = if (length(rows)) do.call(rbind, rows) else data.frame()
  )
}

candidate_paths <- discover_rds()
inventory <- safe_file_info(candidate_paths)
inventory$source_hint <- ifelse(
  inventory$path %in%
    normalizePath(EXPLICIT_CANDIDATES, winslash = "/", mustWork = FALSE),
  "explicit_candidate", "directory_discovery"
)
inventory <- inventory[order(!inventory$exists, inventory$path), , drop = FALSE]

write.csv(inventory, file.path(OUTDIR, "01_RDS_inventory.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

existing_paths <- inventory$path[inventory$exists]
if (!length(existing_paths)) {
  stop("No candidate human MASH RDS file found. Inventory was written to:\n",
       file.path(OUTDIR, "01_RDS_inventory.csv"))
}

object_summary_list <- list()
metadata_summary_list <- list()
sample_counts_list <- list()
crosscheck_list <- list()
structure_list <- list()
expected_labels_list <- list()
audit_log <- character(0)

for (i in seq_along(existing_paths)) {
  pth <- existing_paths[[i]]
  nm <- basename(pth)
  message("[", i, "/", length(existing_paths), "] Reading: ", pth)
  gc(verbose = FALSE)

  obj <- tryCatch(
    readRDS(pth),
    error = function(e) {
      audit_log <<- c(audit_log,
                      paste0("READ ERROR | ", pth, " | ",
                             conditionMessage(e)))
      NULL
    }
  )
  if (is.null(obj)) next

  object_class <- paste(class(obj), collapse = "/")
  n_cells <- tryCatch(ncol(obj), error = function(e) NA_integer_)
  n_features <- tryCatch(nrow(obj), error = function(e) NA_integer_)

  md <- tryCatch({
    if (inherits(obj, "Seurat")) {
      obj[[]]
    } else if ("meta.data" %in% slotNames(obj)) {
      obj@meta.data
    } else data.frame()
  }, error = function(e) data.frame())

  cell_names <- tryCatch(colnames(obj), error = function(e) character(0))
  st <- audit_structure(obj, nm, pth)
  sample_col <- if (ncol(md)) pick_sample_col(md) else NA_character_

  object_summary_list[[length(object_summary_list) + 1L]] <- data.frame(
    object = nm, path = pth, class = object_class,
    n_cells = n_cells, n_features = n_features,
    n_metadata_columns = ncol(md),
    default_assay = st$default_assay,
    assays = paste(st$assays, collapse = " | "),
    reductions = paste(st$reductions, collapse = " | "),
    selected_sample_column = sample_col,
    stringsAsFactors = FALSE
  )

  z <- meta_col_summary(md, nm, pth)
  if (nrow(z)) metadata_summary_list[[length(metadata_summary_list) + 1L]] <- z

  z <- make_sample_counts(md, sample_col, nm, pth)
  if (nrow(z)) sample_counts_list[[length(sample_counts_list) + 1L]] <- z

  z <- make_crosscheck(md, sample_col, nm, pth)
  if (nrow(z)) crosscheck_list[[length(crosscheck_list) + 1L]] <- z

  if (nrow(st$table))
    structure_list[[length(structure_list) + 1L]] <- st$table

  z <- audit_expected_labels(md, sample_col, cell_names, nm, pth)
  if (nrow(z)) expected_labels_list[[length(expected_labels_list) + 1L]] <- z

  audit_log <- c(
    audit_log,
    paste0(
      "OK | ", nm,
      " | cells=", n_cells,
      " | features=", n_features,
      " | metadata_cols=", ncol(md),
      " | sample_col=", ifelse(is.na(sample_col), "<none>", sample_col),
      " | assays=", paste(st$assays, collapse = ","),
      " | reductions=", paste(st$reductions, collapse = ",")
    )
  )

  rm(obj, md)
  gc(verbose = FALSE)
}

rbind_fill <- function(x) {
  if (!length(x)) return(data.frame())

  all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))

  x2 <- lapply(x, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols)) {
      for (nm in missing_cols) df[[nm]] <- NA_character_
    }
    df <- df[, all_cols, drop = FALSE]
    df
  })

  do.call(rbind, x2)
}

bind_or_empty <- function(x) rbind_fill(x)

object_summary <- bind_or_empty(object_summary_list)
metadata_summary <- bind_or_empty(metadata_summary_list)
sample_counts <- bind_or_empty(sample_counts_list)
crosscheck <- bind_or_empty(crosscheck_list)
structure_tbl <- bind_or_empty(structure_list)
expected_tbl <- bind_or_empty(expected_labels_list)

write.csv(object_summary, file.path(OUTDIR, "02_object_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(metadata_summary, file.path(OUTDIR, "03_metadata_columns.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sample_counts, file.path(OUTDIR, "04_sample_cell_counts.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(crosscheck, file.path(OUTDIR, "05_sample_metadata_crosscheck.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(structure_tbl, file.path(OUTDIR, "06_assays_layers_reductions.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(expected_tbl, file.path(OUTDIR, "07_expected_sample_label_audit.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

summary_lines <- c(
  paste0("Human MASH RDS inventory + metadata audit ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("SSD root: ", SSD_ROOT),
  paste0("Output directory: ", OUTDIR),
  "",
  "READ-ONLY POLICY:",
  "- No object was normalized, integrated, reclustered, or re-embedded.",
  "- No annotation was modified.",
  "- No Seurat RDS was saved or overwritten.",
  "",
  paste0("Inventory candidate paths: ", nrow(inventory)),
  paste0("Existing candidate RDS files: ", length(existing_paths)),
  paste0("Successfully audited objects: ", nrow(object_summary)),
  "",
  "AUDIT LOG:",
  audit_log,
  "",
  "OUTPUT FILES:",
  "01_RDS_inventory.csv",
  "02_object_summary.csv",
  "03_metadata_columns.csv",
  "04_sample_cell_counts.csv",
  "05_sample_metadata_crosscheck.csv",
  "06_assays_layers_reductions.csv",
  "07_expected_sample_label_audit.csv",
  "08_audit_summary.txt"
)

writeLines(summary_lines, file.path(OUTDIR, "08_audit_summary.txt"),
           useBytes = TRUE)

message("=== AUDIT COMPLETE ===")
message("Successfully audited objects: ", nrow(object_summary))
message("Output directory: ", OUTDIR)
message("No RDS object was modified or saved.")
