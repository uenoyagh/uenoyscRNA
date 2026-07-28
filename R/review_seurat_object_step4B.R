#' Review a Seurat object
#'
#' Performs core preflight checks before generating annotation-review plots.
#' Review settings are automatically detected unless explicitly supplied.
#'
#' @param object A Seurat object.
#' @param annotation_column Optional metadata column containing annotations.
#' @param sample_column Optional metadata column containing sample identifiers.
#' @param condition_column Optional metadata column containing condition labels.
#' @param reduction Optional dimensional reduction name.
#' @param assay Optional assay name.
#' @param marker_registry Optional marker registry data.frame or CSV path.
#' @param species Optional species filter for the marker registry.
#' @param tissue Optional tissue filter for the marker registry.
#' @param layer Optional layer filter for the marker registry.
#' @param output_dir Output directory used by the review manifest.
#'
#' @return An object of class `uenoy_review`.
#' @export
review_seurat_object <- function(
    object,
    annotation_column = NULL,
    sample_column = NULL,
    condition_column = NULL,
    reduction = NULL,
    assay = NULL,
    marker_registry = NULL,
    species = NULL,
    tissue = NULL,
    layer = NULL,
    output_dir = "annotation_review"
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  detected <- detect_review_settings(
    object = object,
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column,
    reduction = reduction,
    assay = assay
  )

  annotation_column <- detected$annotation_column
  sample_column <- detected$sample_column
  condition_column <- detected$condition_column
  reduction <- detected$reduction
  assay <- detected$assay

  metadata_check <- check_review_metadata(
    object = object,
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column
  )

  marker_coverage <- NULL
  if (!is.null(marker_registry)) {
    if (is.character(marker_registry) && length(marker_registry) == 1L) {
      marker_registry <- read_marker_registry(marker_registry)
    } else {
      validate_marker_registry(marker_registry)
    }

    marker_table <- get_markers(
      registry = marker_registry,
      species = species,
      tissue = tissue,
      layer = layer,
      direction = "positive",
      unique_genes = FALSE
    )

    if (nrow(marker_table) > 0L) {
      split_markers <- split(marker_table$gene, marker_table$celltype)
      coverage_rows <- lapply(names(split_markers), function(celltype) {
        resolved <- resolve_marker_features(
          object = object,
          markers = split_markers[[celltype]],
          assay = assay
        )
        data.frame(
          celltype = celltype,
          requested_n = length(resolved$requested),
          present_n = length(resolved$present),
          missing_n = length(resolved$missing),
          coverage = resolved$coverage,
          missing_genes = paste(resolved$missing, collapse = ";"),
          stringsAsFactors = FALSE
        )
      })
      marker_coverage <- do.call(rbind, coverage_rows)
      rownames(marker_coverage) <- NULL
    } else {
      marker_coverage <- data.frame(
        celltype = character(),
        requested_n = integer(),
        present_n = integer(),
        missing_n = integer(),
        coverage = numeric(),
        missing_genes = character(),
        stringsAsFactors = FALSE
      )
    }
  }

  result <- list(
    summary = list(
      n_cells = ncol(object),
      n_features = nrow(object[[assay]]),
      assay = assay,
      reduction = reduction,
      annotation_column = annotation_column,
      sample_column = sample_column,
      condition_column = condition_column
    ),
    detected_settings = detected,
    metadata = metadata_check,
    marker_coverage = marker_coverage,
    manifest = create_review_manifest(output_dir = output_dir)
  )

  class(result) <- c("uenoy_review", "list")
  result
}
