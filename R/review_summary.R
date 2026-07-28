#' Build a review summary
#'
#' Creates a compact summary of a Seurat object or an existing
#' `uenoy_review` result. When a Seurat object is supplied, review settings
#' are detected automatically unless explicitly specified.
#'
#' @param x A Seurat object or an object returned by
#'   [review_seurat_object()].
#' @param annotation_column Optional annotation metadata column.
#' @param sample_column Optional sample metadata column.
#' @param condition_column Optional condition metadata column.
#' @param reduction Optional reduction name.
#' @param assay Optional assay name.
#' @param marker_registry Optional marker registry data.frame or CSV path.
#' @param species Optional marker-registry species filter.
#' @param tissue Optional marker-registry tissue filter.
#' @param layer Optional marker-registry layer filter.
#' @param output_dir Output directory used when `x` is a Seurat object.
#' @param marker_warning_threshold Marker-coverage threshold below which a
#'   warning is reported.
#'
#' @return An object of class `uenoy_review_summary`.
#' @export
review_summary <- function(
    x,
    annotation_column = NULL,
    sample_column = NULL,
    condition_column = NULL,
    reduction = NULL,
    assay = NULL,
    marker_registry = NULL,
    species = NULL,
    tissue = NULL,
    layer = NULL,
    output_dir = "annotation_review",
    marker_warning_threshold = 0.80
) {
  if (!is.numeric(marker_warning_threshold) ||
      length(marker_warning_threshold) != 1L ||
      is.na(marker_warning_threshold) ||
      marker_warning_threshold < 0 ||
      marker_warning_threshold > 1) {
    stop(
      "`marker_warning_threshold` must be a single number between 0 and 1.",
      call. = FALSE
    )
  }

  if (inherits(x, "Seurat")) {
    review <- review_seurat_object(
      object = x,
      annotation_column = annotation_column,
      sample_column = sample_column,
      condition_column = condition_column,
      reduction = reduction,
      assay = assay,
      marker_registry = marker_registry,
      species = species,
      tissue = tissue,
      layer = layer,
      output_dir = output_dir
    )
    object <- x
  } else if (inherits(x, "uenoy_review")) {
    review <- x
    object <- NULL
  } else {
    stop(
      "`x` must be a Seurat object or an object returned by ",
      "`review_seurat_object()`.",
      call. = FALSE
    )
  }

  annotation_levels <- .review_unique_values(
    object = object,
    metadata_table = review$metadata,
    column = review$summary$annotation_column,
    fallback_name = "annotation_values"
  )

  sample_levels <- .review_unique_values(
    object = object,
    metadata_table = review$metadata,
    column = review$summary$sample_column,
    fallback_name = "sample_values"
  )

  condition_levels <- .review_unique_values(
    object = object,
    metadata_table = review$metadata,
    column = review$summary$condition_column,
    fallback_name = "condition_values"
  )

  marker_coverage <- review$marker_coverage

  low_coverage <- NULL
  if (!is.null(marker_coverage) && nrow(marker_coverage) > 0L) {
    low_coverage <- marker_coverage[
      !is.na(marker_coverage$coverage) &
        marker_coverage$coverage < marker_warning_threshold,
      ,
      drop = FALSE
    ]
  }

  warnings <- character()

  if (is.null(review$summary$sample_column)) {
    warnings <- c(warnings, "No sample column was detected.")
  }

  if (is.null(review$summary$condition_column)) {
    warnings <- c(warnings, "No condition column was detected.")
  }

  if (!is.null(low_coverage) && nrow(low_coverage) > 0L) {
    warnings <- c(
      warnings,
      paste0(
        nrow(low_coverage),
        " cell type(s) have marker coverage below ",
        formatC(marker_warning_threshold * 100, format = "fg"),
        "%."
      )
    )
  }

  result <- list(
    object = list(
      n_cells = review$summary$n_cells,
      n_features = review$summary$n_features
    ),
    settings = list(
      assay = review$summary$assay,
      reduction = review$summary$reduction,
      annotation_column = review$summary$annotation_column,
      sample_column = review$summary$sample_column,
      condition_column = review$summary$condition_column
    ),
    counts = list(
      n_annotations = length(annotation_levels),
      n_samples = length(sample_levels),
      n_conditions = length(condition_levels)
    ),
    values = list(
      annotations = annotation_levels,
      samples = sample_levels,
      conditions = condition_levels
    ),
    marker_coverage = marker_coverage,
    low_marker_coverage = low_coverage,
    warnings = warnings,
    detected_settings = review$detected_settings,
    manifest = review$manifest
  )

  class(result) <- c("uenoy_review_summary", "list")
  result
}

#' Print a review summary
#'
#' @param x A `uenoy_review_summary` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.uenoy_review_summary <- function(x, ...) {
  cat("\n")
  cat("=================================\n")
  cat("uenoyscRNA Review Summary\n")
  cat("=================================\n\n")

  cat("Object\n")
  cat("  Cells        : ", format(x$object$n_cells, big.mark = ","), "\n", sep = "")
  cat("  Features     : ", format(x$object$n_features, big.mark = ","), "\n\n", sep = "")

  cat("Detected settings\n")
  .print_review_setting("Assay", x$settings$assay)
  .print_review_setting("Reduction", x$settings$reduction)
  .print_review_setting("Annotation", x$settings$annotation_column)
  .print_review_setting("Sample", x$settings$sample_column)
  .print_review_setting("Condition", x$settings$condition_column)
  cat("\n")

  cat("Metadata overview\n")
  cat("  Annotation levels : ", x$counts$n_annotations, "\n", sep = "")
  cat("  Samples           : ", x$counts$n_samples, "\n", sep = "")
  cat("  Conditions        : ", x$counts$n_conditions, "\n", sep = "")

  if (length(x$values$conditions) > 0L) {
    cat(
      "  Condition values  : ",
      paste(x$values$conditions, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  cat("\n")

  if (!is.null(x$marker_coverage)) {
    cat("Marker coverage\n")

    if (nrow(x$marker_coverage) == 0L) {
      cat("  No matching marker sets were found.\n")
    } else {
      coverage_table <- x$marker_coverage
      coverage_table$coverage_percent <- paste0(
        round(coverage_table$coverage * 100, 1),
        "%"
      )

      for (i in seq_len(nrow(coverage_table))) {
        cat(
          "  ",
          coverage_table$celltype[[i]],
          " : ",
          coverage_table$coverage_percent[[i]],
          "\n",
          sep = ""
        )
      }
    }
    cat("\n")
  }

  if (length(x$warnings) == 0L) {
    cat("Status\n")
    cat("  OK: no review warnings detected.\n")
  } else {
    cat("Warnings\n")
    for (warning_text in x$warnings) {
      cat("  - ", warning_text, "\n", sep = "")
    }
  }

  cat("\n")
  invisible(x)
}

#' Convert a review summary to data frames
#'
#' Returns tabular components suitable for writing to CSV files or including
#' in reports.
#'
#' @param x A `uenoy_review_summary` object.
#'
#' @return A named list of data.frames.
#' @export
as_review_tables <- function(x) {
  if (!inherits(x, "uenoy_review_summary")) {
    stop("`x` must inherit from class `uenoy_review_summary`.", call. = FALSE)
  }

  settings <- data.frame(
    item = c(
      "n_cells",
      "n_features",
      "assay",
      "reduction",
      "annotation_column",
      "sample_column",
      "condition_column",
      "n_annotations",
      "n_samples",
      "n_conditions"
    ),
    value = c(
      x$object$n_cells,
      x$object$n_features,
      x$settings$assay,
      x$settings$reduction,
      x$settings$annotation_column %||% NA_character_,
      x$settings$sample_column %||% NA_character_,
      x$settings$condition_column %||% NA_character_,
      x$counts$n_annotations,
      x$counts$n_samples,
      x$counts$n_conditions
    ),
    stringsAsFactors = FALSE
  )

  values <- do.call(
    rbind,
    list(
      .review_value_table("annotation", x$values$annotations),
      .review_value_table("sample", x$values$samples),
      .review_value_table("condition", x$values$conditions)
    )
  )
  rownames(values) <- NULL

  warnings <- data.frame(
    warning = x$warnings,
    stringsAsFactors = FALSE
  )

  list(
    settings = settings,
    values = values,
    marker_coverage = x$marker_coverage,
    low_marker_coverage = x$low_marker_coverage,
    warnings = warnings
  )
}

#' Write review summary tables
#'
#' @param x A `uenoy_review_summary` object.
#' @param output_dir Directory for CSV output.
#'
#' @return A named character vector of written file paths.
#' @export
write_review_summary <- function(
    x,
    output_dir = file.path("annotation_review", "summary")
) {
  tables <- as_review_tables(x)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- character()

  for (name in names(tables)) {
    table <- tables[[name]]

    if (is.null(table)) {
      next
    }

    path <- file.path(output_dir, paste0(name, ".csv"))
    utils::write.csv(table, path, row.names = FALSE, na = "")
    paths[[name]] <- path
  }

  paths
}

# Internal helpers ---------------------------------------------------------

.review_unique_values <- function(
    object,
    metadata_table,
    column,
    fallback_name
) {
  if (is.null(column)) {
    return(character())
  }

  if (!is.null(object)) {
    values <- object[[]][[column]]
    values <- as.character(values)
    return(sort(unique(values[!is.na(values) & nzchar(values)])))
  }

  if (!is.null(metadata_table) &&
      is.list(metadata_table) &&
      fallback_name %in% names(metadata_table)) {
    values <- as.character(metadata_table[[fallback_name]])
    return(sort(unique(values[!is.na(values) & nzchar(values)])))
  }

  character()
}

.print_review_setting <- function(label, value) {
  if (is.null(value) || length(value) == 0L || is.na(value)) {
    cat("  [!] ", label, " : not detected\n", sep = "")
  } else {
    cat("  [OK] ", label, " : ", value, "\n", sep = "")
  }
}

.review_value_table <- function(type, values) {
  if (length(values) == 0L) {
    return(data.frame(
      type = character(),
      value = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    type = rep(type, length(values)),
    value = values,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
