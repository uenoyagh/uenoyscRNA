# ============================================================
# RDS discovery and dataset registry helpers
# uenoy scRNAseq Framework v4.0
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_path_soft <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

rds_default_exclude_dirs <- c(
  ".git", ".Rproj.user", "renv", "library", "packrat",
  "analysis_results", "00_inventory", "08_logs"
)

rds_default_exclude_file_patterns <- c(
  "(?:^|_)analysis_bundle\\.rds$",
  "(?:^|_)bundle\\.rds$",
  "object_summaries\\.rds$",
  "(?:^|_)intermediate\\.rds$",
  "without_object\\.rds$"
)

#' Discover candidate RDS files under a project root
#'
#' Discovery is filename/path based and therefore does not load large objects.
#' Files in package libraries, renv, inventory output, and known bundle files are
#' excluded by default.
#'
#' @param root Directory to scan.
#' @param recursive Search recursively.
#' @param exclude_dirs Directory names to exclude.
#' @param exclude_file_patterns Regular expressions matched against basenames.
#' @return Data frame with one row per candidate RDS file.
#' @export
discover_rds_files <- function(
    root,
    recursive = TRUE,
    exclude_dirs = rds_default_exclude_dirs,
    exclude_file_patterns = rds_default_exclude_file_patterns) {
  if (!dir.exists(root)) {
    return(data.frame(
      path = character(), file = character(), relative_path = character(),
      size_gb = numeric(), modified = as.POSIXct(character()),
      species_guess = character(), scope_guess = character(),
      stringsAsFactors = FALSE
    ))
  }

  paths <- list.files(
    root, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE,
    recursive = recursive, all.files = FALSE, include.dirs = FALSE
  )
  if (!length(paths)) {
    return(data.frame(
      path = character(), file = character(), relative_path = character(),
      size_gb = numeric(), modified = as.POSIXct(character()),
      species_guess = character(), scope_guess = character(),
      stringsAsFactors = FALSE
    ))
  }

  paths <- normalize_path_soft(paths)
  root_norm <- sub("/+$", "", normalize_path_soft(root))
  rel <- ifelse(
    startsWith(paths, paste0(root_norm, "/")),
    substring(paths, nchar(root_norm) + 2L),
    basename(paths)
  )
  segments <- strsplit(rel, "/", fixed = TRUE)
  excluded_dir <- vapply(
    segments,
    function(x) any(tolower(x) %in% tolower(exclude_dirs)),
    logical(1)
  )
  basenames <- basename(paths)
  excluded_file <- Reduce(
    `|`,
    lapply(exclude_file_patterns, grepl, x = basenames, ignore.case = TRUE, perl = TRUE),
    init = rep(FALSE, length(paths))
  )
  keep <- !excluded_dir & !excluded_file
  paths <- paths[keep]
  rel <- rel[keep]
  basenames <- basenames[keep]

  if (!length(paths)) {
    return(data.frame(
      path = character(), file = character(), relative_path = character(),
      size_gb = numeric(), modified = as.POSIXct(character()),
      species_guess = character(), scope_guess = character(),
      stringsAsFactors = FALSE
    ))
  }

  info <- file.info(paths)
  search_text <- tolower(paste(rel, basenames))
  species <- ifelse(
    grepl("mouse|gse325222|mm10|grcm", search_text), "Mouse",
    ifelse(grepl("human|lt[0-9]|sflb|hg38", search_text), "Human", "Unknown")
  )
  scope <- ifelse(
    grepl("mphi|macrophage|kupffer|monocyte", search_text), "Macrophage",
    ifelse(grepl("mesenchymal|hsc", search_text), "Mesenchymal", "All_cell_or_other")
  )

  data.frame(
    path = paths,
    file = basenames,
    relative_path = rel,
    size_gb = unname(info$size) / 1024^3,
    modified = info$mtime,
    species_guess = species,
    scope_guess = scope,
    stringsAsFactors = FALSE
  )
}

pattern_hits_all <- function(text, patterns) {
  if (is.null(patterns) || !length(patterns)) return(rep(TRUE, length(text)))
  Reduce(`&`, lapply(patterns, grepl, x = text, ignore.case = TRUE, perl = TRUE))
}

pattern_hits_any <- function(text, patterns) {
  if (is.null(patterns) || !length(patterns)) return(rep(FALSE, length(text)))
  Reduce(`|`, lapply(patterns, grepl, x = text, ignore.case = TRUE, perl = TRUE))
}

#' Resolve RDS files for one registered dataset
#'
#' Preferred explicit files are used first. Missing preferred files trigger a
#' path-based fallback search under the external data root.
#'
#' @param dataset_name Registry key.
#' @param registry Dataset registry list.
#' @param external_data_root Root searched by fallback rules.
#' @param catalog Optional precomputed result from [discover_rds_files()].
#' @param strict Stop when no file resolves.
#' @return Character vector of normalized RDS paths.
#' @export
resolve_dataset_rds <- function(
    dataset_name,
    registry,
    external_data_root,
    catalog = NULL,
    strict = FALSE) {
  if (!dataset_name %in% names(registry)) {
    stop("Unknown dataset: ", dataset_name, call. = FALSE)
  }
  spec <- registry[[dataset_name]]
  preferred <- spec$preferred_files %||% character()
  preferred <- unique(normalize_path_soft(preferred))
  existing <- preferred[file.exists(preferred)]
  if (length(existing)) {
    selection <- spec$selection %||% "latest"
    if (identical(selection, "all")) return(existing)
    return(existing[[1L]])
  }

  if (is.null(catalog)) catalog <- discover_rds_files(external_data_root)
  if (!nrow(catalog)) {
    msg <- paste0("No candidate RDS files were discovered under: ", external_data_root)
    if (strict) stop(msg, call. = FALSE) else return(character())
  }

  text <- paste(catalog$relative_path, catalog$file)
  keep <- pattern_hits_all(text, spec$include_patterns %||% character())
  keep <- keep & !pattern_hits_any(text, spec$exclude_patterns %||% character())
  hits <- catalog[keep, , drop = FALSE]
  if (!nrow(hits)) {
    msg <- paste0(
      "No RDS matched dataset '", dataset_name, "'. Preferred paths were missing and fallback patterns returned no matches."
    )
    if (strict) stop(msg, call. = FALSE) else return(character())
  }

  selection <- spec$selection %||% "latest"
  if (identical(selection, "largest")) {
    hits <- hits[order(-hits$size_gb, hits$file), , drop = FALSE]
    hits <- hits[1L, , drop = FALSE]
  } else if (identical(selection, "all")) {
    hits <- hits[order(hits$relative_path), , drop = FALSE]
  } else {
    hits <- hits[order(hits$modified, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    hits <- hits[1L, , drop = FALSE]
  }
  hits$path
}

#' Summarize registered datasets without loading Seurat objects
#'
#' @param registry Dataset registry.
#' @param external_data_root Project data root.
#' @param catalog Optional RDS catalog.
#' @return Data frame describing resolved files.
#' @export
summarize_dataset_registry <- function(registry, external_data_root, catalog = NULL) {
  if (is.null(catalog)) catalog <- discover_rds_files(external_data_root)
  rows <- lapply(names(registry), function(nm) {
    paths <- resolve_dataset_rds(nm, registry, external_data_root, catalog, strict = FALSE)
    data.frame(
      dataset = nm,
      label = registry[[nm]]$label %||% nm,
      n_resolved = length(paths),
      resolved_files = paste(paths, collapse = " | "),
      status = if (length(paths)) "READY" else "NOT_FOUND",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

print_dataset_registry_summary <- function(summary_df) {
  cat("\n=== uenoyscRNA dataset registry ===\n")
  for (i in seq_len(nrow(summary_df))) {
    mark <- if (identical(summary_df$status[[i]], "READY")) "[OK]" else "[--]"
    cat(sprintf("%-4s %-28s %d file(s)\n", mark, summary_df$dataset[[i]], summary_df$n_resolved[[i]]))
  }
  invisible(summary_df)
}
