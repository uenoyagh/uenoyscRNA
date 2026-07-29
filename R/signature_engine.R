# ============================================================
# Unified signature analysis engine
# uenoy scRNAseq Framework v2.5
# ============================================================

match_features_ci <- function(features, genes) {
  idx <- match(toupper(genes), toupper(features))
  features[idx[!is.na(idx)]]
}

get_assay_data_safe <- function(object, assay = NULL) {

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  if (!assay %in% SeuratObject::Assays(object)) {
    stop(
      "Assay not found: ", assay,
      "\nAvailable assays: ",
      paste(SeuratObject::Assays(object), collapse = ", ")
    )
  }

  assay_object <- object[[assay]]

  available_layers <- SeuratObject::Layers(assay_object)

  # ----------------------------------------------------------
  # 1. Exact "data" layer
  # ----------------------------------------------------------
  if ("data" %in% available_layers) {
    return(
      SeuratObject::LayerData(
        object = assay_object,
        layer = "data"
      )
    )
  }

  # ----------------------------------------------------------
  # 2. Multiple normalized layers: data.1, data.2, ...
  # ----------------------------------------------------------
  data_layers <- available_layers[
    grepl("^data\\.", available_layers)
  ]

  if (length(data_layers) > 0) {

    message(
      "Joining normalized layers in assay '",
      assay,
      "': ",
      paste(data_layers, collapse = ", ")
    )

    joined_object <- SeuratObject::JoinLayers(
      object = object,
      assay = assay,
      layers = data_layers,
      new = "data"
    )

    return(
      SeuratObject::LayerData(
        object = joined_object[[assay]],
        layer = "data"
      )
    )
  }

  # ----------------------------------------------------------
  # 3. No normalized data layer: join counts and normalize copy
  # ----------------------------------------------------------
  count_layers <- available_layers[
    grepl("^counts($|\\.)", available_layers)
  ]

  if (length(count_layers) == 0) {
    stop(
      "No usable data or counts layer was found in assay: ",
      assay,
      "\nAvailable layers: ",
      paste(available_layers, collapse = ", ")
    )
  }

  temporary_object <- object

  if (
    length(count_layers) > 1 ||
    !identical(count_layers, "counts")
  ) {
    message(
      "Joining count layers in assay '",
      assay,
      "': ",
      paste(count_layers, collapse = ", ")
    )

    temporary_object <- SeuratObject::JoinLayers(
      object = temporary_object,
      assay = assay,
      layers = count_layers,
      new = "counts"
    )
  }

  message(
    "Normalized data layer was absent. ",
    "NormalizeData() is being applied to a temporary copy of assay '",
    assay,
    "'."
  )

  temporary_object <- Seurat::NormalizeData(
    object = temporary_object,
    assay = assay,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  SeuratObject::LayerData(
    object = temporary_object[[assay]],
    layer = "data"
  )
}

natural_levels_signature <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values)]
  numeric_values <- suppressWarnings(as.numeric(values))

  if (length(values) > 0 && all(!is.na(numeric_values))) {
    return(values[order(numeric_values)])
  }

  sort(values)
}

resolve_signature_group_column <- function(
    object,
    override = NULL,
    condition_column = NA_character_,
    sample_column = NA_character_,
    annotation_column = NA_character_,
    cluster_column = NA_character_) {

  available <- colnames(object[[]])

  if (!is.null(override)) {
    if (!override %in% available) {
      stop("Configured signature group column was not found: ", override)
    }
    return(override)
  }

  candidates <- c(
    condition_column,
    sample_column,
    annotation_column,
    cluster_column
  )

  candidates <- candidates[!is.na(candidates)]
  candidates <- candidates[candidates %in% available]

  if (length(candidates) == 0) {
    stop("No suitable grouping column was found.")
  }

  candidates[[1]]
}

calculate_cell_signature_scores <- function(
    object,
    signatures,
    assay = NULL,
    min_genes = 2) {

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  expr <- get_assay_data_safe(object, assay)
  features <- rownames(expr)

  score_list <- list()
  report_list <- list()

  for (signature_name in names(signatures)) {
    matched <- match_features_ci(features, signatures[[signature_name]])

    report_list[[signature_name]] <- data.frame(
      signature = signature_name,
      requested_n = length(signatures[[signature_name]]),
      matched_n = length(matched),
      genes_used = paste(matched, collapse = " | "),
      missing_genes = paste(
        setdiff(
          signatures[[signature_name]],
          signatures[[signature_name]][
            toupper(signatures[[signature_name]]) %in% toupper(matched)
          ]
        ),
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )

    if (length(matched) < min_genes) next

    score_list[[signature_name]] <- Matrix::colMeans(
      expr[matched, , drop = FALSE]
    )
  }

  if (length(score_list) == 0) {
    stop("No signature had enough matched genes.")
  }

  score_df <- as.data.frame(score_list, check.names = FALSE)
  rownames(score_df) <- colnames(object)

  list(
    scores = score_df,
    gene_report = do.call(rbind, report_list)
  )
}

calculate_cluster_signature_summary <- function(
    score_df,
    group_values,
    group_name) {

  long <- tidyr::pivot_longer(
    cbind(
      data.frame(
        cell = rownames(score_df),
        group = group_values,
        stringsAsFactors = FALSE
      ),
      score_df
    ),
    cols = -c(cell, group),
    names_to = "signature",
    values_to = "score"
  )

  summary_df <- stats::aggregate(
    score ~ group + signature,
    long,
    mean
  )

  summary_df$group <- factor(
    summary_df$group,
    levels = natural_levels_signature(summary_df$group)
  )

  summary_df$z_score <- ave(
    summary_df$score,
    summary_df$signature,
    FUN = function(x) {
      z <- as.numeric(scale(x))
      z[is.na(z)] <- 0
      z
    }
  )

  names(summary_df)[names(summary_df) == "group"] <- group_name
  summary_df
}

build_signature_footer <- function(
    rds_name,
    framework_name,
    framework_version) {

  paste0(
    "RDS: ", rds_name,
    "    |    Created: ", format(Sys.time(), "%Y-%m-%d %H:%M"),
    "    |    ", framework_name, " v", framework_version
  )
}

add_signature_footer <- function(
    plot,
    footer_text,
    footer_size = 6.5) {

  plot +
    patchwork::plot_annotation(
      caption = footer_text,
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(
          size = footer_size,
          colour = "grey30",
          hjust = 0,
          margin = ggplot2::margin(t = 8)
        )
      )
    )
}

save_plot_dual <- function(
    plot,
    base_path,
    width,
    height,
    export_pdf = TRUE,
    export_png = TRUE,
    png_dpi = 300,
    png_background = "white",
    overwrite = FALSE) {

  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(export_pdf)) {
    pdf_path <- paste0(base_path, ".pdf")

    if (!file.exists(pdf_path) || isTRUE(overwrite)) {
      grDevices::pdf(
        pdf_path,
        width = width,
        height = height,
        useDingbats = FALSE
      )
      print(plot)
      grDevices::dev.off()
    }
  }

  if (isTRUE(export_png)) {
    png_path <- paste0(base_path, ".png")

    if (!file.exists(png_path) || isTRUE(overwrite)) {
      ggplot2::ggsave(
        filename = png_path,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = png_dpi,
        bg = png_background,
        limitsize = FALSE
      )
    }
  }

  invisible(base_path)
}

make_signature_heatmap_plot <- function(
    summary_df,
    group_column,
    clip = 2,
    low_color = "#0033FF",
    mid_color = "#FFFFFF",
    high_color = "#FF1A1A",
    title = "Functional signature scores") {

  summary_df$z_plot <- pmax(
    pmin(summary_df$z_score, clip),
    -clip
  )

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data[[group_column]],
      y = signature,
      fill = z_plot
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.45
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      limits = c(-clip, clip),
      name = "Group\nZ-score"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Group-level z-scores; red = high score",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 17,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 12,
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        size = 9,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

make_signature_featureplot <- function(
    object,
    reduction,
    score_df,
    low_color = "#0033FF",
    mid_color = "#FFFFFF",
    high_color = "#FF1A1A",
    q_low = 0.02,
    q_high = 0.98,
    max_columns = 4) {

  emb <- as.data.frame(Seurat::Embeddings(object, reduction = reduction))
  emb <- emb[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$cell <- rownames(emb)

  plot_df <- cbind(
    emb,
    score_df[rownames(emb), , drop = FALSE]
  )

  long <- tidyr::pivot_longer(
    plot_df,
    cols = names(score_df),
    names_to = "signature",
    values_to = "score"
  )

  quantile_list <- lapply(
    split(long$score, long$signature),
    function(x) {
      x <- x[is.finite(x)]

      if (length(x) == 0) {
        return(c(low = 0, high = 0))
      }

      limits <- stats::quantile(
        x,
        probs = c(q_low, q_high),
        na.rm = TRUE,
        names = FALSE
      )

      if (length(limits) < 2) {
        limits <- rep(limits[[1]], 2)
      }

      c(
        low = limits[[1]],
        high = limits[[2]]
      )
    }
  )

  long$score_scaled <- mapply(
    FUN = function(value, signature_name) {
      limits <- quantile_list[[as.character(signature_name)]]

      if (
        is.null(limits) ||
        length(limits) < 2 ||
        !is.finite(value)
      ) {
        return(NA_real_)
      }

      max(
        min(value, limits[["high"]]),
        limits[["low"]]
      )
    },
    value = long$score,
    signature_name = long$signature,
    USE.NAMES = FALSE
  )

  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = score_scaled
    )
  ) +
    ggplot2::geom_point(
      size = 0.25,
      alpha = 1,
      stroke = 0
    ) +
    ggplot2::scale_colour_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      name = "Score"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(signature),
      ncol = max_columns,
      scales = "fixed"
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = "Signature score UMAP",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      strip.text = ggplot2::element_text(
        size = 9,
        face = "bold"
      )
    )
}

make_signature_violin_plot <- function(
    score_df,
    group_values,
    show_points = FALSE) {

  long <- tidyr::pivot_longer(
    cbind(
      data.frame(
        group = group_values,
        stringsAsFactors = FALSE
      ),
      score_df
    ),
    cols = -group,
    names_to = "signature",
    values_to = "score"
  )

  long$group <- factor(
    long$group,
    levels = natural_levels_signature(long$group)
  )

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = group, y = score, fill = group)
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.25
    ) +
    ggplot2::geom_boxplot(
      width = 0.10,
      outlier.shape = NA,
      alpha = 0.7,
      linewidth = 0.25
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(signature),
      scales = "free_y",
      ncol = 4
    ) +
    ggplot2::labs(
      title = "Signature scores by group",
      x = NULL,
      y = "Mean normalized expression"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "none",
      strip.text = ggplot2::element_text(
        size = 9,
        face = "bold"
      )
    )

  if (isTRUE(show_points)) {
    p <- p +
      ggplot2::geom_jitter(
        width = 0.12,
        size = 0.08,
        alpha = 0.15
      )
  }

  p
}

make_signature_gene_dotplot <- function(
    object,
    signatures,
    group_column,
    assay = NULL,
    max_genes_per_signature = 8) {

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  features <- rownames(get_assay_data_safe(object, assay))

  selected <- unlist(
    lapply(signatures, function(genes) {
      head(
        match_features_ci(features, genes),
        max_genes_per_signature
      )
    }),
    use.names = FALSE
  )

  selected <- unique(selected)

  if (length(selected) == 0) {
    stop("No genes available for DotPlot.")
  }

  Seurat::DotPlot(
    object,
    features = selected,
    group.by = group_column,
    assay = assay
  ) +
    ggplot2::scale_colour_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0
    ) +
    ggplot2::labs(
      title = "Representative signature genes",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      )
    )
}
