# ============================================================
# Generic cell-fraction calculation engine
# uenoy scRNAseq Framework v3.0
# ============================================================

cf_safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

cf_natural_levels <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values) & nzchar(values)]
  numeric_values <- suppressWarnings(as.numeric(values))
  if (length(values) > 0L && all(!is.na(numeric_values))) {
    return(values[order(numeric_values)])
  }
  sort(values)
}

cf_order_factor <- function(x, preferred = NULL) {
  values <- as.character(x)
  observed <- cf_natural_levels(values)
  if (!is.null(preferred)) {
    levels_out <- c(
      preferred[preferred %in% observed],
      setdiff(observed, preferred)
    )
  } else {
    levels_out <- observed
  }
  factor(values, levels = unique(levels_out))
}

cf_resolve_metadata_column <- function(
    object,
    override = NULL,
    candidates,
    role,
    required = TRUE) {

  available <- colnames(object[[]])

  if (!is.null(override)) {
    if (!override %in% available) {
      stop(
        "Configured ", role, " column not found: ", override,
        "\nAvailable columns: ", paste(available, collapse = ", ")
      )
    }
    return(override)
  }

  hit <- candidates[candidates %in% available]
  if (length(hit) == 0L) {
    if (isTRUE(required)) {
      stop(
        "No suitable ", role, " column was found.",
        "\nCandidates: ", paste(candidates, collapse = ", "),
        "\nAvailable columns: ", paste(available, collapse = ", ")
      )
    }
    return(NA_character_)
  }

  hit[[1L]]
}

cf_resolve_columns <- function(
    object,
    feature_override = NULL,
    parent_override = NULL,
    condition_override = NULL,
    sample_override = NULL,
    feature_candidates,
    parent_candidates,
    condition_candidates,
    sample_candidates) {

  feature <- cf_resolve_metadata_column(
    object, feature_override, feature_candidates, "feature", TRUE
  )

  parent <- cf_resolve_metadata_column(
    object, parent_override, parent_candidates, "parent", FALSE
  )

  condition <- cf_resolve_metadata_column(
    object, condition_override, condition_candidates, "condition", TRUE
  )

  sample <- cf_resolve_metadata_column(
    object, sample_override, sample_candidates, "sample", FALSE
  )

  if (is.na(parent)) parent <- feature

  list(
    feature = feature,
    parent = parent,
    condition = condition,
    sample = sample
  )
}

cf_apply_regex_map <- function(x, regex_map = NULL) {
  out <- as.character(x)
  if (is.null(regex_map) || length(regex_map) == 0L) return(out)

  original <- out
  for (pattern in names(regex_map)) {
    hit <- grepl(pattern, original, ignore.case = TRUE, perl = TRUE)
    out[hit] <- unname(regex_map[[pattern]])
  }
  out
}

cf_filter_values <- function(x, include = NULL, exclude = NULL) {
  keep <- rep(TRUE, length(x))
  values <- as.character(x)
  if (!is.null(include)) keep <- keep & values %in% include
  if (!is.null(exclude)) keep <- keep & !values %in% exclude
  keep
}

cf_prepare_metadata <- function(
    object,
    columns,
    condition_order = NULL,
    condition_regex_map = NULL,
    include_features = NULL,
    exclude_features = NULL,
    include_parents = NULL,
    exclude_parents = NULL,
    drop_na = TRUE) {

  md <- object[[]]
  out <- data.frame(
    cell = rownames(md),
    .cf_feature = as.character(md[[columns$feature]]),
    .cf_parent = as.character(md[[columns$parent]]),
    .cf_condition = as.character(md[[columns$condition]]),
    stringsAsFactors = FALSE
  )

  if (!is.na(columns$sample)) {
    out$.cf_sample <- as.character(md[[columns$sample]])
  } else {
    out$.cf_sample <- NA_character_
  }

  out$.cf_condition <- cf_apply_regex_map(
    out$.cf_condition,
    condition_regex_map
  )

  keep <- cf_filter_values(
    out$.cf_feature,
    include_features,
    exclude_features
  ) &
    cf_filter_values(
      out$.cf_parent,
      include_parents,
      exclude_parents
    )

  if (isTRUE(drop_na)) {
    keep <- keep &
      !is.na(out$.cf_feature) & nzchar(out$.cf_feature) &
      !is.na(out$.cf_parent) & nzchar(out$.cf_parent) &
      !is.na(out$.cf_condition) & nzchar(out$.cf_condition)
  }

  out <- out[keep, , drop = FALSE]
  out$.cf_condition <- cf_order_factor(
    out$.cf_condition,
    condition_order
  )
  out$.cf_feature <- factor(
    out$.cf_feature,
    levels = cf_natural_levels(out$.cf_feature)
  )
  out$.cf_parent <- factor(
    out$.cf_parent,
    levels = cf_natural_levels(out$.cf_parent)
  )
  out
}

cf_complete_counts <- function(md) {
  condition_levels <- levels(md$.cf_condition)
  parent_levels <- levels(md$.cf_parent)
  feature_levels <- levels(md$.cf_feature)

  template <- expand.grid(
    condition = condition_levels,
    parent = parent_levels,
    feature = feature_levels,
    stringsAsFactors = FALSE
  )

  observed <- as.data.frame(
    table(
      condition = md$.cf_condition,
      parent = md$.cf_parent,
      feature = md$.cf_feature
    ),
    stringsAsFactors = FALSE
  )
  names(observed)[4L] <- "cell_count"

  merged <- merge(
    template,
    observed,
    by = c("condition", "parent", "feature"),
    all.x = TRUE,
    sort = FALSE
  )
  merged$cell_count[is.na(merged$cell_count)] <- 0L
  merged
}

cf_calculate_tables <- function(md) {
  count <- cf_complete_counts(md)

  total_condition <- stats::aggregate(
    cell_count ~ condition,
    count,
    sum
  )
  names(total_condition)[2L] <- "total_cells_condition"

  total_parent <- stats::aggregate(
    cell_count ~ condition + parent,
    count,
    sum
  )
  names(total_parent)[3L] <- "total_cells_parent"

  count <- merge(
    count,
    total_condition,
    by = "condition",
    all.x = TRUE,
    sort = FALSE
  )
  count <- merge(
    count,
    total_parent,
    by = c("condition", "parent"),
    all.x = TRUE,
    sort = FALSE
  )

  count$fraction_total_percent <- ifelse(
    count$total_cells_condition > 0,
    100 * count$cell_count / count$total_cells_condition,
    NA_real_
  )
  count$fraction_within_parent_percent <- ifelse(
    count$total_cells_parent > 0,
    100 * count$cell_count / count$total_cells_parent,
    NA_real_
  )

  parent_summary <- stats::aggregate(
    cell_count ~ condition + parent,
    count,
    sum
  )
  parent_summary <- merge(
    parent_summary,
    total_condition,
    by = "condition",
    all.x = TRUE,
    sort = FALSE
  )
  parent_summary$fraction_total_percent <- ifelse(
    parent_summary$total_cells_condition > 0,
    100 * parent_summary$cell_count /
      parent_summary$total_cells_condition,
    NA_real_
  )

  list(feature = count, parent = parent_summary)
}

cf_qc_summary <- function(md, rds_name, analysis_target, profile_name) {
  data.frame(
    rds_name = rds_name,
    analysis_target = analysis_target,
    profile_name = profile_name,
    retained_cells = nrow(md),
    n_conditions = length(unique(md$.cf_condition)),
    n_parents = length(unique(md$.cf_parent)),
    n_features = length(unique(md$.cf_feature)),
    conditions = paste(levels(md$.cf_condition), collapse = " | "),
    parents = paste(levels(md$.cf_parent), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

cf_make_heatmap_table <- function(parent_df, value_column) {
  parent_df[, c("condition", "parent", value_column), drop = FALSE]
}

cf_get_profile <- function(target, profiles, default_profile) {
  if (target %in% names(profiles)) {
    profile <- profiles[[target]]
  } else {
    profile <- default_profile
    warning("No dedicated profile found for ", target, "; default profile used.")
  }

  required <- c(
    "profile_name", "feature_candidates", "parent_candidates",
    "condition_candidates", "sample_candidates",
    "condition_order", "condition_regex_map",
    "feature_label", "parent_label", "total_denominator_label"
  )
  missing <- setdiff(required, names(profile))
  if (length(missing) > 0L) {
    stop(
      "Profile for ", target, " is missing field(s): ",
      paste(missing, collapse = ", ")
    )
  }
  profile
}

cf_get_output_root <- function(
    dataset_name,
    result_folder,
    run_folder,
    create = TRUE) {

  if (!exists("get_result_root", mode = "function")) {
    stop("get_result_root() was not found after loading project_config.R")
  }

  root <- file.path(
    get_result_root(dataset_name),
    result_folder,
    run_folder
  )

  if (isTRUE(create)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  root
}
