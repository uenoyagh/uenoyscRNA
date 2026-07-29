# ============================================================
# Cluster-fraction calculation engine
# uenoy scRNAseq Framework v2.6
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

cf_natural_levels <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  xn <- suppressWarnings(as.numeric(x))
  if (length(x) > 0 && all(!is.na(xn))) x[order(xn)] else sort(x)
}

cf_resolve_metadata_column <- function(object, override = NULL, candidates, role) {
  available <- colnames(object[[]])
  if (!is.null(override)) {
    if (!override %in% available) {
      stop("Configured ", role, " column not found: ", override,
           "\nAvailable columns: ", paste(available, collapse = ", "))
    }
    return(override)
  }
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

cf_resolve_columns <- function(object,
                               cluster_override = NULL,
                               annotation_override = NULL,
                               condition_override = NULL,
                               sample_override = NULL) {
  list(
    cluster = cf_resolve_metadata_column(
      object, cluster_override,
      c("integratedRPCA_snn_res.3", "integratedRPCA_snn_res.3.0",
        "RNA_snn_res.3", "RNA_snn_res.3.0", "seurat_clusters",
        "cluster", "Cluster", "res3.0", "res3"), "cluster"
    ),
    annotation = cf_resolve_metadata_column(
      object, annotation_override,
      c("annotation_group", "cell_annotation", "celltype", "cell_type",
        "layer2", "Layer2", "annotation", "Annotation"), "annotation"
    ),
    condition = cf_resolve_metadata_column(
      object, condition_override,
      c("condition", "Condition", "group", "Group", "treatment",
        "Treatment", "diet", "Diet"), "condition"
    ),
    sample = cf_resolve_metadata_column(
      object, sample_override,
      c("sample", "Sample", "sample_id", "sampleID", "orig.ident",
        "library", "Library"), "sample"
    )
  )
}

cf_infer_condition <- function(sample_values, regex_map) {
  sample_values <- as.character(sample_values)
  out <- rep(NA_character_, length(sample_values))
  for (pattern in names(regex_map)) {
    idx <- is.na(out) & grepl(pattern, sample_values, ignore.case = TRUE)
    out[idx] <- unname(regex_map[[pattern]])
  }
  out[is.na(out)] <- sample_values[is.na(out)]
  out
}

cf_order_factor <- function(x, preferred = NULL) {
  present <- unique(as.character(x))
  present <- present[!is.na(present) & nzchar(present)]
  preferred_present <- preferred[preferred %in% present]
  extra <- setdiff(cf_natural_levels(present), preferred_present)
  factor(as.character(x), levels = c(preferred_present, extra))
}

cf_prepare_metadata <- function(object, columns, condition_order, regex_map) {
  md <- object[[]]
  n <- nrow(md)

  if (is.na(columns$cluster)) stop("No cluster metadata column was detected.")
  if (is.na(columns$annotation)) stop("No annotation metadata column was detected.")
  if (is.na(columns$sample)) {
    md$.cf_sample <- rep("All_cells", n)
  } else {
    md$.cf_sample <- as.character(md[[columns$sample]])
  }

  if (is.na(columns$condition)) {
    md$.cf_condition <- cf_infer_condition(md$.cf_sample, regex_map)
  } else {
    md$.cf_condition <- as.character(md[[columns$condition]])
  }

  md$.cf_cluster <- as.character(md[[columns$cluster]])
  md$.cf_annotation <- as.character(md[[columns$annotation]])

  keep <- !is.na(md$.cf_cluster) & nzchar(md$.cf_cluster) &
    !is.na(md$.cf_annotation) & nzchar(md$.cf_annotation) &
    !is.na(md$.cf_condition) & nzchar(md$.cf_condition)
  md <- md[keep, , drop = FALSE]

  md$.cf_condition <- cf_order_factor(md$.cf_condition, condition_order)
  md$.cf_cluster <- factor(md$.cf_cluster, levels = cf_natural_levels(md$.cf_cluster))
  md$.cf_annotation <- factor(md$.cf_annotation, levels = unique(as.character(md$.cf_annotation)))
  md
}

cf_calculate_tables <- function(md) {
  count <- as.data.frame(
    xtabs(~ .cf_condition + .cf_annotation + .cf_cluster, data = md),
    stringsAsFactors = FALSE
  )
  names(count) <- c("condition", "annotation", "cluster", "cell_count")

  # Preserve zero-valued combinations because transition lines should remain visible.
  total_condition <- aggregate(cell_count ~ condition, count, sum)
  names(total_condition)[2] <- "total_cells_condition"
  count <- merge(count, total_condition, by = "condition", all.x = TRUE, sort = FALSE)

  total_annotation <- aggregate(cell_count ~ condition + annotation, count, sum)
  names(total_annotation)[3] <- "total_cells_annotation"
  count <- merge(count, total_annotation,
                 by = c("condition", "annotation"), all.x = TRUE, sort = FALSE)

  count$fraction_total_percent <- ifelse(
    count$total_cells_condition > 0,
    100 * count$cell_count / count$total_cells_condition,
    NA_real_
  )
  count$fraction_within_annotation_percent <- ifelse(
    count$total_cells_annotation > 0,
    100 * count$cell_count / count$total_cells_annotation,
    NA_real_
  )

  annotation_summary <- aggregate(
    cell_count ~ condition + annotation,
    count,
    sum
  )
  annotation_summary <- merge(
    annotation_summary,
    total_condition,
    by = "condition",
    all.x = TRUE,
    sort = FALSE
  )
  annotation_summary$fraction_total_percent <- ifelse(
    annotation_summary$total_cells_condition > 0,
    100 * annotation_summary$cell_count /
      annotation_summary$total_cells_condition,
    NA_real_
  )

  list(cluster = count, annotation = annotation_summary)
}
