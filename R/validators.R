#' Validate a Seurat object
#'
#' Internal helper used by plotting functions.
#'
#' @param object Object to validate.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_seurat_object <- function(object) {

  if (!inherits(object, "Seurat")) {
    stop(
      "`object` must be a Seurat object.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a character vector
#'
#' @param x Object to validate.
#' @param arg Argument name used in error messages.
#' @param allow_empty Whether an empty vector is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_character_vector <- function(
    x,
    arg,
    allow_empty = FALSE
) {

  if (!is.character(x)) {
    stop(
      sprintf("`%s` must be a character vector.", arg),
      call. = FALSE
    )
  }

  if (!allow_empty && length(x) == 0L) {
    stop(
      sprintf("`%s` must not be empty.", arg),
      call. = FALSE
    )
  }

  if (anyNA(x)) {
    stop(
      sprintf("`%s` must not contain NA values.", arg),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a single string
#'
#' @param x Object to validate.
#' @param arg Argument name used in error messages.
#' @param allow_null Whether NULL is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_single_string <- function(
    x,
    arg,
    allow_null = FALSE
) {

  if (is.null(x)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      sprintf("`%s` must be a single character string.", arg),
      call. = FALSE
    )
  }

  if (
    !is.character(x) ||
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(x)
  ) {
    stop(
      sprintf("`%s` must be a single character string.", arg),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a logical scalar
#'
#' @param x Object to validate.
#' @param arg Argument name used in error messages.
#' @param allow_null Whether NULL is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_logical <- function(
    x,
    arg,
    allow_null = FALSE
) {

  if (is.null(x)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      sprintf("`%s` must be TRUE or FALSE.", arg),
      call. = FALSE
    )
  }

  if (
    !is.logical(x) ||
    length(x) != 1L ||
    is.na(x)
  ) {
    stop(
      sprintf("`%s` must be TRUE or FALSE.", arg),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a numeric scalar
#'
#' @param x Object to validate.
#' @param arg Argument name used in error messages.
#' @param allow_null Whether NULL is allowed.
#' @param finite Whether the value must be finite.
#' @param min Minimum allowed value.
#' @param max Maximum allowed value.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_number <- function(
    x,
    arg,
    allow_null = FALSE,
    finite = TRUE,
    min = -Inf,
    max = Inf
) {

  if (is.null(x)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      sprintf("`%s` must be a single numeric value.", arg),
      call. = FALSE
    )
  }

  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    is.na(x)
  ) {
    stop(
      sprintf("`%s` must be a single numeric value.", arg),
      call. = FALSE
    )
  }

  if (finite && !is.finite(x)) {
    stop(
      sprintf("`%s` must be finite.", arg),
      call. = FALSE
    )
  }

  if (x < min || x > max) {
    stop(
      sprintf(
        "`%s` must be between %s and %s.",
        arg,
        format(min),
        format(max)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate an assay name
#'
#' @param object Seurat object.
#' @param assay Assay name.
#' @param allow_null Whether NULL is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_assay <- function(
    object,
    assay,
    allow_null = TRUE
) {

  .validate_seurat_object(object)

  if (is.null(assay)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      "`assay` must be a single character string.",
      call. = FALSE
    )
  }

  .validate_single_string(
    assay,
    arg = "assay"
  )

  available_assays <- names(object@assays)

  if (!assay %in% available_assays) {
    stop(
      sprintf(
        "Assay `%s` was not found in `object`.",
        assay
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate feature names
#'
#' @param object Seurat object.
#' @param features Character vector of feature names.
#' @param assay Assay name.
#' @param require_all Whether all features must exist.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_features <- function(
    object,
    features,
    assay = NULL,
    require_all = TRUE
) {

  .validate_seurat_object(object)

  .validate_character_vector(
    features,
    arg = "features"
  )

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  .validate_assay(
    object,
    assay,
    allow_null = FALSE
  )

  available_features <- rownames(object[[assay]])

  missing_features <- setdiff(
    features,
    available_features
  )

  if (require_all && length(missing_features) > 0L) {
    stop(
      sprintf(
        "The following features were not found in assay `%s`: %s",
        assay,
        paste(missing_features, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a metadata column
#'
#' @param object Seurat object.
#' @param column Metadata column name.
#' @param arg Argument name used in error messages.
#' @param allow_null Whether NULL is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_metadata_column <- function(
    object,
    column,
    arg = "group.by",
    allow_null = TRUE
) {

  .validate_seurat_object(object)

  if (is.null(column)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      sprintf("`%s` must be a single character string.", arg),
      call. = FALSE
    )
  }

  .validate_single_string(
    column,
    arg = arg
  )

  metadata_columns <- colnames(object[[]])

  if (!column %in% metadata_columns) {
    stop(
      sprintf(
        "Metadata column `%s` was not found in `object`.",
        column
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate group order
#'
#' Checks whether a requested ordering is valid for a metadata column.
#'
#' @param object A Seurat object.
#' @param column Metadata column name.
#' @param order Character vector specifying the desired order.
#' @param arg Argument name used in error messages.
#' @param allow_null Whether NULL is allowed.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.validate_group_order <- function(
    object,
    column,
    order,
    arg = "group_order",
    allow_null = TRUE
) {

  .validate_seurat_object(object)

  .validate_metadata_column(
    object = object,
    column = column,
    arg = "group.by",
    allow_null = FALSE
  )

  if (is.null(order)) {
    if (allow_null) {
      return(invisible(TRUE))
    }

    stop(
      sprintf("`%s` must be a character vector.", arg),
      call. = FALSE
    )
  }

  .validate_character_vector(
    x = order,
    arg = arg
  )

  if (anyDuplicated(order)) {
    stop(
      sprintf("`%s` must not contain duplicated values.", arg),
      call. = FALSE
    )
  }

  observed_values <- unique(
    as.character(object[[]][[column]])
  )

  observed_values <- observed_values[
    !is.na(observed_values)
  ]

  missing_values <- setdiff(
    order,
    observed_values
  )

  if (length(missing_values) > 0L) {
    stop(
      sprintf(
        paste0(
          "The following values in `%s` were not found in ",
          "metadata column `%s`: %s"
        ),
        arg,
        column,
        paste(missing_values, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate a single non-empty character string
#'
#' Compatibility helper used by existing plotting functions.
#'
#' @param x Object to validate.
#' @param argument Argument name used in error messages.
#'
#' @return Invisibly returns the supplied object.
#'
#' @keywords internal
validate_single_string <- function(
    x,
    argument
) {

  if (
    !is.character(x) ||
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(x)
  ) {
    stop(
      paste0(
        "`",
        argument,
        "` must be a single non-empty character string."
      ),
      call. = FALSE
    )
  }

  invisible(x)
}


#' Validate an optional single non-empty character string
#'
#' Compatibility helper used by existing plotting functions.
#'
#' @param x Object to validate.
#' @param argument Argument name used in error messages.
#'
#' @return Invisibly returns the supplied object.
#'
#' @keywords internal
validate_optional_single_string <- function(
    x,
    argument
) {

  if (is.null(x)) {
    return(invisible(x))
  }

  validate_single_string(
    x = x,
    argument = argument
  )
}


#' Validate a positive numeric scalar
#'
#' Compatibility helper used by existing plotting functions.
#'
#' @param x Object to validate.
#' @param argument Argument name used in error messages.
#'
#' @return Invisibly returns the supplied object.
#'
#' @keywords internal
validate_positive_number <- function(
    x,
    argument
) {

  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    is.na(x) ||
    x <= 0
  ) {
    stop(
      paste0(
        "`",
        argument,
        "` must be a single positive number."
      ),
      call. = FALSE
    )
  }

  invisible(x)
}
