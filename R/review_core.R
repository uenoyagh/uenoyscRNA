#' Check metadata required for review
#'
#' Validates that requested metadata columns are present in a Seurat object's
#' metadata and reports available values.
#'
#' @param object A Seurat object.
#' @param annotation_column Metadata column containing annotations.
#' @param sample_column Optional metadata column containing sample identifiers.
#' @param condition_column Optional metadata column containing condition labels.
#'
#' @return A list with validated column names and observed values.
#' @export
check_review_metadata <- function(
    object,
    annotation_column,
    sample_column = NULL,
    condition_column = NULL
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  metadata <- object[[]]

  required <- c(annotation_column, sample_column, condition_column)
  required <- required[!is.null(required) & !is.na(required) & nzchar(required)]

  missing_columns <- setdiff(required, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  values <- lapply(required, function(column) {
    unique(as.character(metadata[[column]]))
  })
  names(values) <- required

  list(
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column,
    values = values
  )
}

#' Resolve marker features against a Seurat object
#'
#' Compares requested markers with features present in the selected assay.
#'
#' @param object A Seurat object.
#' @param markers Character vector of marker genes.
#' @param assay Assay name. Defaults to Seurat's default assay.
#'
#' @return A list containing requested, present, missing, and coverage.
#' @export
resolve_marker_features <- function(object, markers, assay = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }
  if (!is.character(markers)) {
    stop("`markers` must be a character vector.", call. = FALSE)
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object)

  if (!assay %in% names(object@assays)) {
    stop("Assay not found in object: ", assay, call. = FALSE)
  }

  available <- rownames(object[[assay]])
  requested <- unique(markers)
  present <- requested[requested %in% available]
  missing <- requested[!requested %in% available]

  coverage <- if (length(requested) == 0L) {
    NA_real_
  } else {
    length(present) / length(requested)
  }

  list(
    assay = assay,
    requested = requested,
    present = present,
    missing = missing,
    coverage = coverage
  )
}

#' Create a review manifest
#'
#' Builds a manifest of review outputs that later plotting stages can execute.
#'
#' @param include Character vector of review components.
#' @param output_dir Output directory for the review.
#'
#' @return A data.frame describing planned outputs.
#' @export
create_review_manifest <- function(
    include = c(
      "qc",
      "annotation",
      "marker_coverage",
      "umap",
      "dotplot",
      "featureplot",
      "violin",
      "module_score"
    ),
    output_dir = "annotation_review"
) {
  supported <- c(
    "qc",
    "annotation",
    "marker_coverage",
    "umap",
    "dotplot",
    "featureplot",
    "violin",
    "module_score"
  )

  invalid <- setdiff(include, supported)
  if (length(invalid) > 0L) {
    stop(
      "Unsupported review component(s): ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  directory_map <- c(
    qc = "01_qc",
    annotation = "02_annotation",
    marker_coverage = "03_marker_coverage",
    umap = "04_umap",
    dotplot = "05_dotplot",
    featureplot = "06_featureplot",
    violin = "07_violin",
    module_score = "08_module_score"
  )

  data.frame(
    component = include,
    directory = unname(directory_map[include]),
    path = file.path(output_dir, unname(directory_map[include])),
    status = "planned",
    stringsAsFactors = FALSE
  )
}


# Internal null-coalescing helper
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
