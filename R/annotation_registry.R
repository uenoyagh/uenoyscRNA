#' Read an annotation registry
#'
#' @param path Path to a CSV file.
#' @param validate Logical; validate after reading.
#' @return A data frame.
#' @export
read_annotation_registry <- function(path, validate = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single, non-empty character string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Annotation registry file does not exist: ", path, call. = FALSE)
  }

  registry <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA"),
    colClasses = "character"
  )

  if (isTRUE(validate)) {
    validate_annotation_registry(registry)
  }
  registry
}

#' Validate an annotation registry
#'
#' The unique key is species + tissue + dataset + layer + cluster.
#'
#' @param registry A data frame.
#' @param allowed_confidence Allowed confidence labels.
#' @return The input registry, invisibly.
#' @export
validate_annotation_registry <- function(
    registry,
    allowed_confidence = c("High", "Medium", "Low")
) {
  required <- c(
    "species", "tissue", "dataset", "layer", "cluster",
    "annotation", "confidence", "markers", "evidence",
    "reviewer", "review_date"
  )

  if (!is.data.frame(registry)) {
    stop("`registry` must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required, names(registry))
  if (length(missing_columns) > 0L) {
    stop(
      "Annotation registry is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(registry) == 0L) {
    stop("Annotation registry must contain at least one row.", call. = FALSE)
  }

  essential <- c(
    "species", "tissue", "dataset", "layer",
    "cluster", "annotation", "confidence"
  )
  for (column in essential) {
    values <- trimws(as.character(registry[[column]]))
    bad <- is.na(registry[[column]]) | values == ""
    if (any(bad)) {
      stop(
        "Column `", column, "` contains missing or empty value(s) at row(s): ",
        paste(which(bad), collapse = ", "),
        call. = FALSE
      )
    }
  }

  bad_confidence <- !registry$confidence %in% allowed_confidence
  if (any(bad_confidence)) {
    stop(
      "Invalid confidence value(s) at row(s): ",
      paste(which(bad_confidence), collapse = ", "),
      ". Allowed values: ",
      paste(allowed_confidence, collapse = ", "),
      call. = FALSE
    )
  }

  key_columns <- c("species", "tissue", "dataset", "layer", "cluster")
  key <- do.call(
    paste,
    c(lapply(registry[key_columns], as.character), sep = "\r")
  )
  duplicate_key <- duplicated(key) | duplicated(key, fromLast = TRUE)
  if (any(duplicate_key)) {
    stop(
      "Duplicate annotation key(s) detected at row(s): ",
      paste(which(duplicate_key), collapse = ", "),
      call. = FALSE
    )
  }

  has_date <- !is.na(registry$review_date) &
    trimws(as.character(registry$review_date)) != ""
  if (any(has_date)) {
    date_text <- as.character(registry$review_date[has_date])
    parsed <- as.Date(date_text, format = "%Y-%m-%d")
    invalid <- is.na(parsed) | format(parsed, "%Y-%m-%d") != date_text
    if (any(invalid)) {
      rows <- which(has_date)[invalid]
      stop(
        "Invalid `review_date` value(s) at row(s): ",
        paste(rows, collapse = ", "),
        ". Use YYYY-MM-DD format.",
        call. = FALSE
      )
    }
  }

  invisible(registry)
}

#' Get cluster annotations
#'
#' @param registry A validated annotation registry.
#' @param cluster One or more cluster identifiers.
#' @param species Species label.
#' @param tissue Tissue label.
#' @param dataset Dataset identifier.
#' @param layer Annotation layer.
#' @param strict Error if requested clusters are absent.
#' @param return Either `"record"` or `"annotation"`.
#' @return A data frame or named character vector.
#' @export
get_cluster_annotation <- function(
    registry,
    cluster,
    species,
    tissue,
    dataset,
    layer,
    strict = TRUE,
    return = c("record", "annotation")
) {
  validate_annotation_registry(registry)
  return <- match.arg(return)
  requested <- as.character(cluster)

  keep <- registry$species == species &
    registry$tissue == tissue &
    registry$dataset == dataset &
    registry$layer == layer

  available <- registry[keep, , drop = FALSE]
  found <- requested %in% as.character(available$cluster)

  if (isTRUE(strict) && any(!found)) {
    stop(
      "Cluster(s) not found in the selected registry context: ",
      paste(requested[!found], collapse = ", "),
      call. = FALSE
    )
  }

  requested <- requested[found]
  result <- available[
    match(requested, as.character(available$cluster)),
    ,
    drop = FALSE
  ]
  rownames(result) <- NULL

  if (identical(return, "annotation")) {
    value <- as.character(result$annotation)
    names(value) <- requested
    return(value)
  }

  result
}

#' Write an annotation registry
#'
#' @param registry A data frame.
#' @param path Destination CSV path.
#' @param overwrite Logical; replace an existing file.
#' @return The normalized output path, invisibly.
#' @export
write_annotation_registry <- function(registry, path, overwrite = FALSE) {
  validate_annotation_registry(registry)

  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "File already exists. Set `overwrite = TRUE` to replace it: ",
      path,
      call. = FALSE
    )
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    registry,
    path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}
