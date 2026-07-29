#' Create marker-review FeaturePlots
#'
#' Creates one or more Seurat FeaturePlots for reviewing marker-expression
#' patterns on dimensional-reduction embeddings.
#'
#' @param object A Seurat object.
#' @param features Character vector of genes to plot.
#' @param assay Optional assay name.
#' @param reduction Optional dimensional-reduction name.
#' @param split_by Optional metadata column passed to `Seurat::FeaturePlot()`.
#' @param pt_size Optional point size passed to `Seurat::FeaturePlot()`.
#' @param order Logical; draw high-expression cells last.
#' @param keep_scale Scale behavior passed to `Seurat::FeaturePlot()`.
#' @param raster Optional logical value passed to `Seurat::FeaturePlot()`.
#'   When `NULL`, Seurat determines whether rasterization is used.
#' @param low Low-expression color.
#' @param mid Midpoint color.
#' @param high High-expression color.
#' @param midpoint Midpoint used by the diverging color scale.
#' @param quantile_cutoff Numeric vector of length two defining lower and upper
#'   quantile clipping of feature expression.
#' @param max_features_per_plot Maximum number of genes in each PDF panel.
#' @param ncol Number of plot columns.
#' @param rds_file Optional source RDS path or filename shown in the footer.
#' @param analysis_date Analysis date shown in the footer. Defaults to `Sys.Date()`.
#' @param show_provenance Show the RDS filename and analysis date.
#' @param provenance_size Footer text size.
#' @param provenance_color Footer text color.
#' @param output_dir Optional output directory. If `NULL`, plots are returned
#'   without writing files.
#' @param width PDF width.
#' @param height PDF height.
#'
#' @return An object of class `uenoy_review_featureplot`.
#' @export
review_featureplot <- function(
    object,
    features,
    assay = NULL,
    reduction = NULL,
    split_by = NULL,
    pt_size = NULL,
    order = TRUE,
    keep_scale = "feature",
    raster = NULL,
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    quantile_cutoff = c(0.02, 0.98),
    max_features_per_plot = 12L,
    ncol = 4L,
    rds_file = NULL,
    analysis_date = Sys.Date(),
    show_provenance = TRUE,
    provenance_size = 6,
    provenance_color = "#666666",
    output_dir = NULL,
    width = 12,
    height = 9
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from class `Seurat`.", call. = FALSE)
  }

  if (!is.character(features) ||
      length(features) == 0L ||
      anyNA(features) ||
      any(!nzchar(features))) {
    stop(
      "`features` must be a non-empty character vector without missing or empty values.",
      call. = FALSE
    )
  }

  if (!is.numeric(quantile_cutoff) ||
      length(quantile_cutoff) != 2L ||
      anyNA(quantile_cutoff) ||
      quantile_cutoff[1] < 0 ||
      quantile_cutoff[2] > 1 ||
      quantile_cutoff[1] >= quantile_cutoff[2]) {
    stop(
      paste0(
        "`quantile_cutoff` must contain two increasing ",
        "numbers between 0 and 1."
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(max_features_per_plot) ||
      length(max_features_per_plot) != 1L ||
      is.na(max_features_per_plot) ||
      max_features_per_plot < 1) {
    stop("`max_features_per_plot` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(ncol) ||
      length(ncol) != 1L ||
      is.na(ncol) ||
      ncol < 1) {
    stop("`ncol` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(split_by)) {
    if (!is.character(split_by) ||
        length(split_by) != 1L ||
        is.na(split_by) ||
        !nzchar(split_by)) {
      stop(
        "`split_by` must be `NULL` or a single non-empty character value.",
        call. = FALSE
      )
    }

    if (!split_by %in% colnames(object[[]])) {
      stop(
        paste0("`split_by` column was not found in object metadata: ", split_by),
        call. = FALSE
      )
    }
  }

  if (!is.null(raster) &&
      (!is.logical(raster) || length(raster) != 1L || is.na(raster))) {
    stop("`raster` must be `NULL`, `TRUE`, or `FALSE`.", call. = FALSE)
  }

  if (!is.logical(show_provenance) || length(show_provenance) != 1L || is.na(show_provenance)) {
    stop("`show_provenance` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(provenance_size) || length(provenance_size) != 1L ||
      is.na(provenance_size) || provenance_size <= 0) {
    stop("`provenance_size` must be one positive number.", call. = FALSE)
  }
  if (!is.character(provenance_color) || length(provenance_color) != 1L ||
      is.na(provenance_color) || !nzchar(provenance_color)) {
    stop("`provenance_color` must be one non-empty character value.", call. = FALSE)
  }
  analysis_date <- .review_normalize_analysis_date(analysis_date)

  settings <- detect_review_settings(
    object = object,
    assay = assay,
    reduction = reduction,
    require_reduction = TRUE
  )

  features <- unique(features)
  present_features <- intersect(features, rownames(object))
  missing_features <- setdiff(features, present_features)

  if (length(present_features) == 0L) {
    stop("None of the supplied features were found in the object.", call. = FALSE)
  }

  chunks <- split(
    present_features,
    ceiling(seq_along(present_features) / as.integer(max_features_per_plot))
  )

  plots <- vector("list", length(chunks))
  names(plots) <- sprintf("featureplot_%02d", seq_along(chunks))

  old_assay <- SeuratObject::DefaultAssay(object)
  on.exit(SeuratObject::DefaultAssay(object) <- old_assay, add = TRUE)
  SeuratObject::DefaultAssay(object) <- settings$assay

  for (i in seq_along(chunks)) {
    panel_plots <- list()

    for (feature in chunks[[i]]) {
      featureplot_args <- list(
        object = object,
        features = feature,
        assay = settings$assay,
        reduction = settings$reduction,
        split.by = split_by,
        pt.size = pt_size,
        order = order,
        keep.scale = keep_scale,
        combine = FALSE
      )

      if (!is.null(raster)) {
        featureplot_args$raster <- raster
      }

      feature_plots <- do.call(Seurat::FeaturePlot, featureplot_args)

      if (inherits(feature_plots, "ggplot")) {
        feature_plots <- list(feature_plots)
      }

      if (!is.list(feature_plots) ||
          length(feature_plots) == 0L ||
          !all(vapply(feature_plots, inherits, logical(1), what = "ggplot"))) {
        stop(
          paste0(
            "Seurat::FeaturePlot() returned an unexpected object for feature `",
            feature,
            "`."
          ),
          call. = FALSE
        )
      }

      feature_plots <- .apply_review_featureplot_scale(
        plots = feature_plots,
        feature = feature,
        low = low,
        mid = mid,
        high = high,
        midpoint = midpoint,
        quantile_cutoff = quantile_cutoff
      )

      panel_plots <- c(panel_plots, feature_plots)
    }

    plots[[i]] <- patchwork::wrap_plots(
      panel_plots,
      ncol = as.integer(ncol)
    ) +
      patchwork::plot_annotation(
        title = paste0("FeaturePlot Review ", i, " / ", length(chunks)),
        caption = if (isTRUE(show_provenance)) {
          .review_provenance_text(
            rds_file = rds_file,
            analysis_date = analysis_date
          )
        } else {
          NULL
        },
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          plot.caption = ggplot2::element_text(
            size = provenance_size,
            colour = provenance_color,
            hjust = 1,
            margin = ggplot2::margin(t = 5, r = 2, unit = "pt")
          )
        )
      )
  }

  files <- character()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    for (name in names(plots)) {
      path <- file.path(output_dir, paste0(name, ".pdf"))

      grDevices::pdf(
        file = path,
        width = width,
        height = height,
        useDingbats = FALSE
      )

      tryCatch(
        print(plots[[name]]),
        finally = grDevices::dev.off()
      )

      files[[name]] <- path
    }
  }

  coverage <- data.frame(
    feature = features,
    present = features %in% rownames(object),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      plots = plots,
      files = files,
      settings = settings,
      coverage = coverage,
      present_features = present_features,
      missing_features = missing_features,
      provenance = list(
        rds_file = rds_file,
        analysis_date = analysis_date,
        show = show_provenance
      )
    ),
    class = c("uenoy_review_featureplot", "list")
  )
}


#' Print a FeaturePlot review result
#'
#' @param x A `uenoy_review_featureplot` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.uenoy_review_featureplot <- function(x, ...) {
  cat("\n")
  cat("uenoyscRNA FeaturePlot Review\n")
  cat("-----------------------------\n")
  cat("Plots generated   : ", length(x$plots), "\n", sep = "")
  cat("Files written     : ", length(x$files), "\n", sep = "")
  cat("Features supplied : ", nrow(x$coverage), "\n", sep = "")
  cat("Missing features  : ", length(x$missing_features), "\n\n", sep = "")
  invisible(x)
}


# Internal helpers ---------------------------------------------------------

.apply_review_featureplot_scale <- function(
    plots,
    feature,
    low,
    mid,
    high,
    midpoint,
    quantile_cutoff
) {
  color_columns <- vapply(
    plots,
    .find_featureplot_color_column,
    character(1),
    feature = feature
  )

  values <- unlist(
    Map(
      function(plot, color_column) {
        if (is.na(color_column)) {
          return(numeric())
        }
        as.numeric(plot$data[[color_column]])
      },
      plots,
      color_columns
    ),
    use.names = FALSE
  )

  limits <- .review_featureplot_limits(
    values = values,
    quantile_cutoff = quantile_cutoff
  )

  Map(
    function(plot, color_column) {
      if (is.na(color_column)) {
        warning(
          paste0(
            "Could not identify the expression column for feature `",
            feature,
            "`. The original FeaturePlot color scale was retained."
          ),
          call. = FALSE
        )
        return(.style_review_featureplot(plot))
      }

      plot$data[[color_column]] <- pmin(
        pmax(as.numeric(plot$data[[color_column]]), limits[1]),
        limits[2]
      )

      plot$scales$remove("colour")
      plot$scales$remove("color")

      plot <- plot +
        ggplot2::scale_color_gradient2(
          low = low,
          mid = mid,
          high = high,
          midpoint = midpoint,
          limits = limits,
          oob = scales::squish,
          name = feature
        )

      .style_review_featureplot(plot)
    },
    plots,
    color_columns
  )
}


.review_featureplot_limits <- function(values, quantile_cutoff) {
  finite_values <- values[is.finite(values)]

  if (length(finite_values) == 0L) {
    return(c(-1, 1))
  }

  limits <- as.numeric(
    stats::quantile(
      finite_values,
      probs = quantile_cutoff,
      na.rm = TRUE,
      names = FALSE
    )
  )

  if (!all(is.finite(limits)) || limits[1] == limits[2]) {
    center <- unique(finite_values)

    if (length(center) == 1L) {
      spread <- max(abs(center), 1)
      limits <- c(center - spread, center + spread)
    } else {
      spread <- max(abs(finite_values), na.rm = TRUE)
      if (!is.finite(spread) || spread == 0) {
        spread <- 1
      }
      limits <- c(-spread, spread)
    }
  }

  limits
}


.style_review_featureplot <- function(plot) {
  plot +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )
}


.find_featureplot_color_column <- function(plot, feature) {
  data <- plot$data

  if (!is.data.frame(data) || nrow(data) == 0L) {
    return(NA_character_)
  }

  if (feature %in% colnames(data) && is.numeric(data[[feature]])) {
    return(feature)
  }

  feature_safe <- make.names(feature)
  if (feature_safe %in% colnames(data) && is.numeric(data[[feature_safe]])) {
    return(feature_safe)
  }

  mapped_color <- tryCatch(
    rlang::as_name(plot$mapping$colour),
    error = function(e) {
      tryCatch(
        rlang::as_name(plot$mapping$color),
        error = function(e2) NA_character_
      )
    }
  )

  if (!is.na(mapped_color) &&
      mapped_color %in% colnames(data) &&
      is.numeric(data[[mapped_color]])) {
    return(mapped_color)
  }

  numeric_columns <- names(data)[
    vapply(data, is.numeric, logical(1))
  ]

  coordinate_patterns <- c(
    "(_|\\.)1$",
    "(_|\\.)2$",
    "^UMAP",
    "^tSNE",
    "^TSNE",
    "^PC",
    "^PCA",
    "^DC",
    "^LSI",
    "^harmony"
  )

  is_coordinate <- vapply(
    numeric_columns,
    function(column) {
      any(vapply(
        coordinate_patterns,
        grepl,
        logical(1),
        x = column,
        ignore.case = TRUE
      ))
    },
    logical(1)
  )

  candidates <- numeric_columns[!is_coordinate]

  if (length(candidates) == 0L) {
    return(NA_character_)
  }

  candidates[[length(candidates)]]
}
