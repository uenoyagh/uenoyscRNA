#' Detect shared metadata roles in a Seurat object
#'
#' Identifies likely cell-type, cluster, condition, and sample columns. Explicit
#' overrides always take priority. The returned report records the selected
#' column, detection source, alternative candidates, and number of unique values.
#'
#' @param object A Seurat object.
#' @param overrides Named list containing any of `celltype`, `cluster`,
#'   `condition`, and `sample`.
#' @param required_roles Roles that must be detected.
#' @param strict Stop when a required role is unavailable.
#'
#' @return An object of class `uenoy_metadata_map`.
#' @export
#'
#' @examples
#' \dontrun{
#' map <- detect_metadata(object, overrides = list(
#'   celltype = "feature_annotation_added"
#' ))
#' map$report
#' }
detect_metadata <- function(
    object,
    overrides = list(),
    required_roles = c("celltype", "condition"),
    strict = TRUE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  available <- colnames(object[[]])
  rules <- list(
    celltype = list(
      exact = c(
        "feature_annotation_added", "annotation_group_final",
        "annotation_group", "cell_annotation", "celltype", "cell_type",
        "CellType", "layer2", "Layer2", "layer1", "Layer1",
        "annotation", "Annotation"
      ),
      patterns = c("feature.*annotation", "annotation", "cell.?type", "^layer[12]$")
    ),
    cluster = list(
      exact = c("seurat_clusters", "cluster", "Cluster"),
      patterns = c(
        "^integratedRPCA_snn_res\\.", "^integrated_snn_res\\.",
        "^RNA_snn_res\\.", "^SCT_snn_res\\.", "^.*_snn_res\\."
      )
    ),
    condition = list(
      exact = c(
        "NAS", "nas", "condition", "Condition", "group", "Group",
        "treatment", "Treatment", "diet", "Diet"
      ),
      patterns = c("^NAS$", "condition", "group", "treatment", "diet", "disease", "status")
    ),
    sample = list(
      exact = c("sample", "Sample", "sample_id", "SampleID", "orig.ident"),
      patterns = c("^sample$", "sample.?id", "orig\\.ident", "donor", "patient", "replicate")
    )
  )

  roles <- names(rules)
  selected <- setNames(rep(NA_character_, length(roles)), roles)
  reports <- vector("list", length(roles))

  for (i in seq_along(roles)) {
    role <- roles[[i]]
    override <- overrides[[role]]
    rule <- rules[[role]]

    if (!is.null(override)) {
      if (length(override) != 1L || !override %in% available) {
        stop(
          "Metadata override for '", role, "' was not found: ",
          paste(override, collapse = ", "),
          call. = FALSE
        )
      }
      hits <- override
      source <- "override"
    } else {
      exact_hits <- rule$exact[rule$exact %in% available]
      pattern_hits <- unique(unlist(lapply(
        rule$patterns,
        function(pattern) grep(
          pattern, available, value = TRUE,
          ignore.case = TRUE, perl = TRUE
        )
      )))
      hits <- unique(c(exact_hits, pattern_hits))
      if (identical(role, "cluster") && length(hits) > 1L) {
        resolution <- suppressWarnings(as.numeric(
          sub("^.*(?:res\\.|res_)", "", hits, perl = TRUE)
        ))
        hits <- hits[order(is.na(resolution), -resolution)]
      }
      source <- if (!length(hits)) {
        "not_found"
      } else if (hits[[1L]] %in% exact_hits) {
        "exact"
      } else {
        "pattern"
      }
    }

    if (length(hits)) selected[[role]] <- hits[[1L]]
    reports[[i]] <- data.frame(
      role = role,
      selected_column = selected[[role]],
      detection_source = source,
      alternatives = if (length(hits)) {
        paste(setdiff(hits, selected[[role]]), collapse = " | ")
      } else "",
      n_unique = if (!is.na(selected[[role]])) {
        length(unique(object[[]][[selected[[role]]]]))
      } else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  missing_required <- required_roles[is.na(selected[required_roles])]
  if (length(missing_required)) {
    message <- paste0(
      "Required metadata roles were not detected: ",
      paste(missing_required, collapse = ", "),
      ". Available columns: ", paste(available, collapse = ", ")
    )
    if (isTRUE(strict)) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  }

  result <- as.list(selected)
  result$report <- do.call(rbind, reports)
  class(result) <- c("uenoy_metadata_map", "list")
  result
}

#' Write a metadata-detection report
#'
#' @param metadata_map Result from [detect_metadata()].
#' @param path Output CSV path.
#'
#' @return The normalized output path, invisibly.
#' @export
write_metadata_detection_report <- function(metadata_map, path) {
  if (!inherits(metadata_map, "uenoy_metadata_map")) {
    stop("`metadata_map` must be returned by detect_metadata().", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(metadata_map$report, path, row.names = FALSE, na = "")
  invisible(normalizePath(path, mustWork = FALSE))
}
