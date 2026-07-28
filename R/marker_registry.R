#' Read a marker registry
#'
#' Reads a marker registry CSV file and validates its schema.
#'
#' @param path Path to a marker registry CSV file.
#' @param validate Logical; validate the registry after reading.
#'
#' @return A data.frame containing marker definitions.
#' @export
read_marker_registry <- function(path, validate = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Marker registry file does not exist: ", path, call. = FALSE)
  }

  registry <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )

  character_columns <- c(
    "species", "tissue", "layer", "celltype", "marker_group",
    "gene", "direction", "evidence", "notes"
  )
  for (column in intersect(character_columns, names(registry))) {
    registry[[column]] <- as.character(registry[[column]])
  }

  if (isTRUE(validate)) {
    validate_marker_registry(registry)
  }

  registry
}

#' Validate a marker registry
#'
#' Checks required columns, missing values, permitted directions, and duplicate
#' marker definitions.
#'
#' @param registry A marker registry data.frame.
#'
#' @return Invisibly returns TRUE when validation succeeds.
#' @export
validate_marker_registry <- function(registry) {
  if (!is.data.frame(registry)) {
    stop("`registry` must be a data.frame.", call. = FALSE)
  }

  required_columns <- c(
    "species", "tissue", "layer", "celltype", "marker_group",
    "gene", "direction", "evidence", "notes"
  )
  missing_columns <- setdiff(required_columns, names(registry))
  if (length(missing_columns) > 0L) {
    stop(
      "Marker registry is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  required_non_missing <- c(
    "species", "tissue", "celltype", "marker_group", "gene", "direction"
  )
  for (column in required_non_missing) {
    values <- registry[[column]]
    if (any(is.na(values) | !nzchar(trimws(values)))) {
      stop("Column `", column, "` contains missing or empty values.", call. = FALSE)
    }
  }

  allowed_directions <- c("positive", "negative")
  invalid_directions <- setdiff(unique(registry$direction), allowed_directions)
  if (length(invalid_directions) > 0L) {
    stop(
      "Invalid `direction` value(s): ",
      paste(invalid_directions, collapse = ", "),
      ". Allowed values are: positive, negative.",
      call. = FALSE
    )
  }

  duplicate_key <- paste(
    registry$species,
    registry$tissue,
    registry$layer,
    registry$celltype,
    registry$marker_group,
    registry$gene,
    registry$direction,
    sep = "\r"
  )
  if (anyDuplicated(duplicate_key)) {
    stop(
      "Marker registry contains duplicate marker definitions.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Get markers from a marker registry
#'
#' Filters markers by species, tissue, layer, cell type, marker group, and
#' direction.
#'
#' @param registry A marker registry data.frame or path to a registry CSV file.
#' @param species Optional species filter.
#' @param tissue Optional tissue filter.
#' @param layer Optional annotation layer filter.
#' @param celltype Optional cell type filter.
#' @param marker_group Optional marker group filter.
#' @param direction Optional direction filter: `"positive"` or `"negative"`.
#' @param unique_genes Logical; return a unique character vector of genes.
#'
#' @return A filtered data.frame, or a character vector when
#'   `unique_genes = TRUE`.
#' @export
get_markers <- function(
    registry,
    species = NULL,
    tissue = NULL,
    layer = NULL,
    celltype = NULL,
    marker_group = NULL,
    direction = "positive",
    unique_genes = TRUE
) {
  if (is.character(registry) && length(registry) == 1L) {
    registry <- read_marker_registry(registry)
  } else {
    validate_marker_registry(registry)
  }

  filters <- list(
    species = species,
    tissue = tissue,
    layer = layer,
    celltype = celltype,
    marker_group = marker_group,
    direction = direction
  )

  out <- registry
  for (column in names(filters)) {
    value <- filters[[column]]
    if (!is.null(value)) {
      out <- out[out[[column]] %in% value, , drop = FALSE]
    }
  }

  if (isTRUE(unique_genes)) {
    return(unique(out$gene))
  }

  rownames(out) <- NULL
  out
}

#' List marker sets
#'
#' Returns distinct combinations of marker-set metadata.
#'
#' @param registry A marker registry data.frame or path to a registry CSV file.
#'
#' @return A data.frame of distinct marker-set definitions.
#' @export
list_marker_sets <- function(registry) {
  if (is.character(registry) && length(registry) == 1L) {
    registry <- read_marker_registry(registry)
  } else {
    validate_marker_registry(registry)
  }

  columns <- c(
    "species", "tissue", "layer", "celltype", "marker_group", "direction"
  )
  out <- unique(registry[, columns, drop = FALSE])
  out <- out[do.call(order, out[columns]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Write a marker registry
#'
#' Validates and writes a marker registry to CSV.
#'
#' @param registry A marker registry data.frame.
#' @param path Output CSV path.
#' @param overwrite Logical; overwrite an existing file.
#'
#' @return Invisibly returns the normalized output path.
#' @export
write_marker_registry <- function(registry, path, overwrite = FALSE) {
  validate_marker_registry(registry)

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty character string.", call. = FALSE)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("File already exists. Use `overwrite = TRUE` to replace it.", call. = FALSE)
  }

  output_dir <- dirname(path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  utils::write.csv(registry, path, row.names = FALSE, na = "")
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}
