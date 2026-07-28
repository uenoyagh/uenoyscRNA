#' Detect review-relevant metadata columns
#'
#' Automatically identifies likely annotation, sample, and condition columns
#' in a Seurat object's metadata. Exact-name matches are preferred, followed by
#' suffix-aware and keyword-based matches.
#'
#' @param object A Seurat object.
#' @param annotation_column Optional explicit annotation column.
#' @param sample_column Optional explicit sample column.
#' @param condition_column Optional explicit condition column.
#'
#' @return A list containing detected columns and candidate rankings.
#' @export
detect_review_metadata <- function(
    object,
    annotation_column = NULL,
    sample_column = NULL,
    condition_column = NULL
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  metadata_columns <- colnames(object[[]])

  annotation_candidates <- .rank_review_columns(
    metadata_columns,
    preferred = c(
      "celltype_for_R8plot_FIXED2",
      "celltype_for_R8plot",
      "celltype_annotation",
      "celltype_auto_annotation",
      "celltype",
      "annotation",
      "layer2",
      "layer1"
    ),
    patterns = c(
      "^celltype_for_R8plot_FIXED[0-9]+$",
      "^celltype_for_R8plot$",
      "celltype.*annotation",
      "celltype",
      "annotation",
      "layer2",
      "layer1"
    ),
    exclude_patterns = c("confidence", "feature_annotation", "cluster")
  )

  sample_candidates <- .rank_review_columns(
    metadata_columns,
    preferred = c(
      "sample_display_FIXED2",
      "sample_for_R8plot_FIXED2",
      "sample_for_R8plot",
      "sample",
      "sample_for_annotation",
      "orig.ident"
    ),
    patterns = c(
      "^sample_display_FIXED[0-9]+$",
      "^sample_for_R8plot_FIXED[0-9]+$",
      "^sample_for_R8plot$",
      "^sample$",
      "sample",
      "^orig\\.ident$"
    )
  )

  condition_candidates <- .rank_review_columns(
    metadata_columns,
    preferred = c(
      "condition_FIXED2",
      "condition",
      "group"
    ),
    patterns = c(
      "^condition_FIXED[0-9]+$",
      "^condition$",
      "^group$",
      "condition",
      "group"
    )
  )

  annotation_column <- .resolve_review_column(
    explicit = annotation_column,
    candidates = annotation_candidates,
    metadata_columns = metadata_columns,
    label = "annotation"
  )

  sample_column <- .resolve_review_column(
    explicit = sample_column,
    candidates = sample_candidates,
    metadata_columns = metadata_columns,
    label = "sample",
    allow_null = TRUE
  )

  condition_column <- .resolve_review_column(
    explicit = condition_column,
    candidates = condition_candidates,
    metadata_columns = metadata_columns,
    label = "condition",
    allow_null = TRUE
  )

  list(
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column,
    candidates = list(
      annotation = annotation_candidates,
      sample = sample_candidates,
      condition = condition_candidates
    )
  )
}

#' Detect a review reduction
#'
#' Selects the most appropriate dimensional reduction for annotation review.
#' UMAP-like reductions are preferred, including custom names such as
#' `"umapRPCA"`.
#'
#' @param object A Seurat object.
#' @param reduction Optional explicit reduction name.
#'
#' @return A list containing the selected reduction and ranked candidates.
#' @export
detect_review_reduction <- function(object, reduction = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  reductions <- names(object@reductions)

  if (!is.null(reduction)) {
    if (!reduction %in% reductions) {
      stop("Reduction not found in object: ", reduction, call. = FALSE)
    }

    return(list(
      reduction = reduction,
      candidates = reduction
    ))
  }

  if (length(reductions) == 0L) {
    stop("No dimensional reductions were found in the Seurat object.", call. = FALSE)
  }

  preferred <- c("umap", "umapRPCA", "wnn.umap", "integrated.umap", "tsne", "pca")
  exact <- preferred[preferred %in% reductions]

  umap_like <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  tsne_like <- reductions[grepl("tsne", reductions, ignore.case = TRUE)]
  pca_like <- reductions[grepl("pca", reductions, ignore.case = TRUE)]

  candidates <- unique(c(exact, umap_like, tsne_like, pca_like, reductions))

  list(
    reduction = candidates[[1L]],
    candidates = candidates
  )
}

#' Detect the review assay
#'
#' Uses an explicit assay when supplied; otherwise uses the Seurat default
#' assay and falls back to commonly used assays if necessary.
#'
#' @param object A Seurat object.
#' @param assay Optional explicit assay name.
#'
#' @return A list containing the selected assay and candidates.
#' @export
detect_review_assay <- function(object, assay = NULL) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  assays <- names(object@assays)

  if (!is.null(assay)) {
    if (!assay %in% assays) {
      stop("Assay not found in object: ", assay, call. = FALSE)
    }

    return(list(
      assay = assay,
      candidates = assay
    ))
  }

  default_assay <- SeuratObject::DefaultAssay(object)
  preferred <- c(default_assay, "RNA", "SCT", "integrated")
  candidates <- unique(c(preferred[preferred %in% assays], assays))

  if (length(candidates) == 0L) {
    stop("No assays were found in the Seurat object.", call. = FALSE)
  }

  list(
    assay = candidates[[1L]],
    candidates = candidates
  )
}

#' Automatically detect review settings
#'
#' Detects annotation, sample, condition, reduction, and assay settings for a
#' Seurat object.
#'
#' @param object A Seurat object.
#' @param annotation_column Optional explicit annotation column.
#' @param sample_column Optional explicit sample column.
#' @param condition_column Optional explicit condition column.
#' @param reduction Optional explicit reduction.
#' @param assay Optional explicit assay.
#'
#' @return A list of detected settings and candidate rankings.
#' @export
detect_review_settings <- function(
    object,
    annotation_column = NULL,
    sample_column = NULL,
    condition_column = NULL,
    reduction = NULL,
    assay = NULL
) {
  metadata <- detect_review_metadata(
    object = object,
    annotation_column = annotation_column,
    sample_column = sample_column,
    condition_column = condition_column
  )

  reduction_result <- detect_review_reduction(
    object = object,
    reduction = reduction
  )

  assay_result <- detect_review_assay(
    object = object,
    assay = assay
  )

  list(
    annotation_column = metadata$annotation_column,
    sample_column = metadata$sample_column,
    condition_column = metadata$condition_column,
    reduction = reduction_result$reduction,
    assay = assay_result$assay,
    candidates = list(
      annotation = metadata$candidates$annotation,
      sample = metadata$candidates$sample,
      condition = metadata$candidates$condition,
      reduction = reduction_result$candidates,
      assay = assay_result$candidates
    )
  )
}

# Internal helper: rank metadata columns
.rank_review_columns <- function(
    metadata_columns,
    preferred,
    patterns,
    exclude_patterns = character()
) {
  columns <- metadata_columns

  if (length(exclude_patterns) > 0L) {
    excluded <- Reduce(
      `|`,
      lapply(
        exclude_patterns,
        function(pattern) grepl(pattern, columns, ignore.case = TRUE)
      )
    )
    columns <- columns[!excluded]
  }

  ranked <- character()

  for (name in preferred) {
    ranked <- c(ranked, columns[columns == name])
  }

  for (pattern in patterns) {
    ranked <- c(
      ranked,
      columns[grepl(pattern, columns, ignore.case = TRUE)]
    )
  }

  unique(ranked)
}

# Internal helper: resolve explicit or detected metadata column
.resolve_review_column <- function(
    explicit,
    candidates,
    metadata_columns,
    label,
    allow_null = FALSE
) {
  if (!is.null(explicit)) {
    if (!explicit %in% metadata_columns) {
      stop(
        "Explicit ", label, " column not found: ",
        explicit,
        call. = FALSE
      )
    }
    return(explicit)
  }

  if (length(candidates) == 0L) {
    if (allow_null) {
      return(NULL)
    }

    stop(
      "Could not automatically detect an ", label,
      " column. Supply it explicitly.",
      call. = FALSE
    )
  }

  candidates[[1L]]
}
