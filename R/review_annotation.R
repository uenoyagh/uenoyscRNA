#' Run an annotation review
#'
#' Creates review summaries and UMAP figures for available annotation, sample,
#' and condition metadata columns.
#'
#' @param object A Seurat object.
#' @param annotation_column Optional annotation metadata column.
#' @param sample_column Optional sample metadata column.
#' @param condition_column Optional condition metadata column.
#' @param reduction Optional reduction name.
#' @param assay Optional assay name.
#' @param marker_registry Optional marker registry data frame or CSV path.
#' @param species,tissue,layer Optional marker-registry filters.
#' @param palette Optional annotation palette.
#' @param sample_palette Optional sample palette.
#' @param condition_palette Optional condition palette.
#' @param point_size UMAP point size.
#' @param rds_file Optional source RDS path or filename.
#' @param analysis_date Analysis date printed in figure footers.
#' @param show_provenance Show provenance footers.
#' @param output_dir Output directory.
#' @param width,height PDF dimensions in inches.
#'
#' @return An object of class `uenoy_annotation_review`.
#' @export
review_annotation <- function(
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
    palette = NULL,
    sample_palette = NULL,
    condition_palette = NULL,
    point_size = 0.6,
    rds_file = NULL,
    analysis_date = Sys.Date(),
    show_provenance = TRUE,
    output_dir = "annotation_review",
    width = 8,
    height = 7
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_dir <- file.path(output_dir, "Summary")
  umap_dir <- file.path(output_dir, "UMAP")
  dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)

  summary_result <- review_summary(
    x = object,
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

  summary_files <- write_review_summary(summary_result, output_dir = summary_dir)
  settings <- summary_result$settings

  umap_result <- review_umap(
    object = object,
    annotation_column = settings$annotation_column,
    sample_column = settings$sample_column,
    condition_column = settings$condition_column,
    reduction = settings$reduction,
    assay = settings$assay,
    palette = palette,
    sample_palette = sample_palette,
    condition_palette = condition_palette,
    point_size = point_size,
    rds_file = rds_file,
    analysis_date = analysis_date,
    show_provenance = show_provenance,
    output_dir = umap_dir,
    width = width,
    height = height
  )

  umap_files <- umap_result$files

  manifest <- data.frame(
    category = c(rep("summary", length(summary_files)), rep("umap", length(umap_files))),
    name = c(names(summary_files), names(umap_files)),
    path = c(unname(summary_files), unname(umap_files)),
    stringsAsFactors = FALSE
  )

  manifest_path <- file.path(output_dir, "review_manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)

  result <- list(
    summary = summary_result,
    umap = umap_result,
    manifest = manifest,
    manifest_path = manifest_path,
    output_dir = output_dir
  )
  class(result) <- c("uenoy_annotation_review", "list")
  result
}

#' Print an annotation review result
#'
#' @param x A `uenoy_annotation_review` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.uenoy_annotation_review <- function(x, ...) {
  cat("\n")
  cat("uenoyscRNA Annotation Review\n")
  cat("-----------------------------\n")
  cat("Output directory : ", x$output_dir, "\n", sep = "")
  cat("Files generated  : ", nrow(x$manifest), "\n", sep = "")
  cat("Manifest         : ", x$manifest_path, "\n\n", sep = "")
  invisible(x)
}
