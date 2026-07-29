#' Create marker-review DotPlots
#'
#' Creates one or more Seurat DotPlots for reviewing annotation-marker
#' relationships. Marker genes may be supplied directly as a named list or
#' resolved from a marker registry.
#'
#' @param object A Seurat object.
#' @param markers Optional named list of marker vectors. Names are used as
#'   marker-group labels.
#' @param marker_registry Optional marker registry data.frame or CSV path.
#' @param annotation_column Optional annotation metadata column.
#' @param assay Optional assay name.
#' @param species Optional species filter for the marker registry.
#' @param tissue Optional tissue filter for the marker registry.
#' @param layer Optional annotation-layer filter for the marker registry.
#' @param split_by Optional metadata column passed to `Seurat::DotPlot()`.
#' @param dot_scale Maximum dot size.
#' @param rotate_axis Logical; rotate x-axis labels.
#' @param low Low-expression color.
#' @param mid Midpoint color.
#' @param high High-expression color.
#' @param midpoint Midpoint used by the diverging color scale.
#' @param quantile_cutoff Numeric vector of length two defining lower and upper
#'   quantile clipping of scaled average expression.
#' @param max_features_per_plot Maximum number of genes in each PDF panel.
#' @param output_dir Optional output directory. If `NULL`, plots are returned
#'   without writing files.
#' @param width PDF width.
#' @param height PDF height.
#'
#' @return An object of class `uenoy_review_dotplot`.
#' @export
review_marker_dotplot <- function(
    object,
    markers = NULL,
    marker_registry = NULL,
    annotation_column = NULL,
    assay = NULL,
    species = NULL,
    tissue = NULL,
    layer = NULL,
    split_by = NULL,
    dot_scale = 6,
    rotate_axis = TRUE,
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    quantile_cutoff = c(0.02, 0.98),
    max_features_per_plot = 40L,
    output_dir = NULL,
    width = 12,
    height = 8
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  if (!is.numeric(quantile_cutoff) ||
      length(quantile_cutoff) != 2L ||
      anyNA(quantile_cutoff) ||
      quantile_cutoff[1] < 0 ||
      quantile_cutoff[2] > 1 ||
      quantile_cutoff[1] >= quantile_cutoff[2]) {
    stop(
      "`quantile_cutoff` must contain two increasing numbers between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(max_features_per_plot) ||
      length(max_features_per_plot) != 1L ||
      is.na(max_features_per_plot) ||
      max_features_per_plot < 1) {
    stop("`max_features_per_plot` must be a positive integer.", call. = FALSE)
  }

  settings <- detect_review_settings(
    object = object,
    annotation_column = annotation_column,
    assay = assay,
    require_reduction = FALSE
  )

  marker_list <- .resolve_review_markers(
    object = object,
    markers = markers,
    marker_registry = marker_registry,
    species = species,
    tissue = tissue,
    layer = layer
  )

  marker_list <- lapply(marker_list, unique)
  marker_list <- marker_list[lengths(marker_list) > 0L]

  if (length(marker_list) == 0L) {
    stop("No marker genes were available for DotPlot generation.", call. = FALSE)
  }

  all_features <- unique(unlist(marker_list, use.names = FALSE))
  present_features <- intersect(all_features, rownames(object))
  missing_features <- setdiff(all_features, present_features)

  if (length(present_features) == 0L) {
    stop("None of the supplied marker genes were found in the object.", call. = FALSE)
  }

  marker_list_present <- lapply(
    marker_list,
    function(x) intersect(x, present_features)
  )
  marker_list_present <- marker_list_present[lengths(marker_list_present) > 0L]

  chunks <- .chunk_marker_list(
    marker_list = marker_list_present,
    max_features = as.integer(max_features_per_plot)
  )

  plots <- vector("list", length(chunks))
  names(plots) <- sprintf("marker_dotplot_%02d", seq_along(chunks))

  old_assay <- SeuratObject::DefaultAssay(object)
  on.exit(
    SeuratObject::DefaultAssay(object) <- old_assay,
    add = TRUE
  )
  SeuratObject::DefaultAssay(object) <- settings$assay

  old_idents <- SeuratObject::Idents(object)
  on.exit(
    SeuratObject::Idents(object) <- old_idents,
    add = TRUE
  )
  SeuratObject::Idents(object) <- object[[]][[settings$annotation_column]]

  for (i in seq_along(chunks)) {
    plot <- Seurat::DotPlot(
      object = object,
      features = chunks[[i]],
      assay = settings$assay,
      split.by = split_by,
      dot.scale = dot_scale
    )

    plot <- .apply_review_dotplot_scale(
      plot = plot,
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      quantile_cutoff = quantile_cutoff
    ) +
      ggplot2::labs(
        title = paste0(
          "Marker DotPlot ",
          i,
          " / ",
          length(chunks)
        ),
        x = "Marker gene",
        y = settings$annotation_column,
        color = "Scaled average\nexpression",
        size = "Percent\nexpressed"
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.text.y = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(face = "bold")
      )

    if (isTRUE(rotate_axis)) {
      plot <- plot +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
        )
    }

    plots[[i]] <- plot
  }

  files <- character()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    for (name in names(plots)) {
      path <- file.path(
        output_dir,
        paste0(name, ".pdf")
      )

      grDevices::pdf(
        file = path,
        width = width,
        height = height,
        useDingbats = FALSE
      )

      tryCatch(
        {
          print(plots[[name]])
        },
        finally = {
          grDevices::dev.off()
        }
      )

      files[[name]] <- path
    }
  }

  coverage <- data.frame(
    marker_group = rep(names(marker_list), lengths(marker_list)),
    feature = unlist(marker_list, use.names = FALSE),
    present = unlist(marker_list, use.names = FALSE) %in% rownames(object),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      plots = plots,
      files = files,
      settings = settings,
      marker_groups = marker_list_present,
      coverage = coverage,
      missing_features = missing_features
    ),
    class = c("uenoy_review_dotplot", "list")
  )
}

#' Print a marker DotPlot review result
#'
#' @param x A `uenoy_review_dotplot` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.uenoy_review_dotplot <- function(x, ...) {
  cat("\n")
  cat("uenoyscRNA Marker DotPlot Review\n")
  cat("--------------------------------\n")
  cat("Plots generated   : ", length(x$plots), "\n", sep = "")
  cat("Files written     : ", length(x$files), "\n", sep = "")
  cat("Marker genes      : ", nrow(x$coverage), "\n", sep = "")
  cat("Missing genes     : ", length(x$missing_features), "\n\n", sep = "")
  invisible(x)
}

# Internal helpers ---------------------------------------------------------

.resolve_review_markers <- function(
    object,
    markers,
    marker_registry,
    species,
    tissue,
    layer
) {
  if (!is.null(markers)) {
    if (is.character(markers)) {
      markers <- list(Markers = markers)
    }

    if (!is.list(markers)) {
      stop(
        "`markers` must be a character vector or a named list.",
        call. = FALSE
      )
    }

    if (is.null(names(markers)) ||
        any(!nzchar(names(markers)))) {
      names(markers) <- paste0(
        "Marker_group_",
        seq_along(markers)
      )
    }

    return(lapply(markers, as.character))
  }

  if (is.null(marker_registry)) {
    stop(
      "Supply either `markers` or `marker_registry`.",
      call. = FALSE
    )
  }

  registry <- .read_review_marker_registry(marker_registry)

  species_column <- .find_registry_column(
    registry,
    c("species", "organism")
  )
  tissue_column <- .find_registry_column(
    registry,
    c("tissue", "organ")
  )
  layer_column <- .find_registry_column(
    registry,
    c("layer", "annotation_layer")
  )
  feature_column <- .find_registry_column(
    registry,
    c("feature", "gene", "gene_symbol", "marker")
  )
  group_column <- .find_registry_column(
    registry,
    c(
      "cell_type",
      "celltype",
      "annotation",
      "marker_group",
      "group"
    )
  )

  if (is.null(feature_column)) {
    stop(
      paste0(
        "Could not identify a marker-gene column in `marker_registry`. ",
        "Supported names include: feature, gene, gene_symbol, marker."
      ),
      call. = FALSE
    )
  }

  if (!is.null(species) && !is.null(species_column)) {
    registry <- registry[
      tolower(as.character(registry[[species_column]])) ==
        tolower(species),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(tissue) && !is.null(tissue_column)) {
    registry <- registry[
      tolower(as.character(registry[[tissue_column]])) ==
        tolower(tissue),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(layer) && !is.null(layer_column)) {
    registry <- registry[
      tolower(as.character(registry[[layer_column]])) ==
        tolower(layer),
      ,
      drop = FALSE
    ]
  }

  if (nrow(registry) == 0L) {
    stop(
      "No marker-registry rows remained after filtering.",
      call. = FALSE
    )
  }

  features <- trimws(as.character(registry[[feature_column]]))
  keep <- !is.na(features) & nzchar(features)
  registry <- registry[keep, , drop = FALSE]
  features <- features[keep]

  if (is.null(group_column)) {
    groups <- rep("Markers", length(features))
  } else {
    groups <- trimws(as.character(registry[[group_column]]))
    groups[is.na(groups) | !nzchar(groups)] <- "Unspecified"
  }

  split(features, groups)
}

.read_review_marker_registry <- function(marker_registry) {
  if (is.data.frame(marker_registry)) {
    return(marker_registry)
  }

  if (is.character(marker_registry) &&
      length(marker_registry) == 1L &&
      file.exists(marker_registry)) {
    return(
      utils::read.csv(
        marker_registry,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  stop(
    "`marker_registry` must be a data.frame or an existing CSV path.",
    call. = FALSE
  )
}

.find_registry_column <- function(data, candidates) {
  matched <- match(
    tolower(candidates),
    tolower(colnames(data)),
    nomatch = 0L
  )
  matched <- matched[matched > 0L]

  if (length(matched) == 0L) {
    return(NULL)
  }

  colnames(data)[matched[[1L]]]
}

.chunk_marker_list <- function(marker_list, max_features) {
  chunks <- list()
  current <- list()
  current_n <- 0L

  flush_current <- function() {
    if (length(current) > 0L) {
      chunks[[length(chunks) + 1L]] <<- current
      current <<- list()
      current_n <<- 0L
    }
  }

  for (group_name in names(marker_list)) {
    features <- marker_list[[group_name]]

    while (length(features) > 0L) {
      available <- max_features - current_n

      if (available == 0L) {
        flush_current()
        available <- max_features
      }

      take_n <- min(length(features), available)
      selected <- features[seq_len(take_n)]

      if (group_name %in% names(current)) {
        current[[group_name]] <- c(
          current[[group_name]],
          selected
        )
      } else {
        current[[group_name]] <- selected
      }

      current_n <- current_n + take_n
      features <- features[-seq_len(take_n)]

      if (current_n >= max_features) {
        flush_current()
      }
    }
  }

  flush_current()
  chunks
}

.apply_review_dotplot_scale <- function(
    plot,
    low,
    mid,
    high,
    midpoint,
    quantile_cutoff
) {
  values <- plot$data$avg.exp.scaled
  finite_values <- values[is.finite(values)]

  if (length(finite_values) == 0L) {
    limits <- c(-1, 1)
  } else {
    limits <- as.numeric(
      stats::quantile(
        finite_values,
        probs = quantile_cutoff,
        na.rm = TRUE,
        names = FALSE
      )
    )

    if (!all(is.finite(limits)) ||
        limits[1] == limits[2]) {
      spread <- max(abs(finite_values), na.rm = TRUE)

      if (!is.finite(spread) || spread == 0) {
        spread <- 1
      }

      limits <- c(-spread, spread)
    }
  }

  plot$data$avg.exp.scaled <- pmin(
    pmax(plot$data$avg.exp.scaled, limits[1]),
    limits[2]
  )

  plot +
    ggplot2::scale_color_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = midpoint,
      limits = limits,
      oob = scales::squish
    )
}
