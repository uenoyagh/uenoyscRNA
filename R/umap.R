# ============================================================
# UMAP functions
# uenoy scRNAseq Framework v2.3
# ============================================================

resolve_from_priority <- function(
    available, priority, override = NULL,
    required = FALSE, label = "column") {

  if (!is.null(override)) {
    if (!override %in% available) {
      stop("Configured ", label, " was not found: ", override)
    }
    return(override)
  }

  hit <- priority[priority %in% available]
  if (length(hit) > 0) return(hit[[1]])

  if (isTRUE(required)) {
    stop("No suitable ", label, " was found.")
  }

  NA_character_
}

resolve_umap_reduction <- function(object, override = NULL) {
  resolve_from_priority(
    names(object@reductions),
    c("umapRPCA", "umap.rpca", "integrated_umap", "umap_integrated", "umap"),
    override, TRUE, "UMAP reduction"
  )
}

resolve_sample_column <- function(object, override = NULL) {
  resolve_from_priority(
    colnames(object[[]]),
    c("sample", "sample_for_annotation", "sample_id", "orig.ident", "replicate"),
    override, FALSE, "sample column"
  )
}

resolve_condition_column <- function(object, override = NULL) {
  resolve_from_priority(
    colnames(object[[]]),
    c("condition", "group", "treatment", "disease_status", "status", "diet"),
    override, FALSE, "condition column"
  )
}

resolve_annotation_column <- function(object, override = NULL) {
  resolve_from_priority(
    colnames(object[[]]),
    c(
      "celltype_for_R8plot_FIXED2",
      "celltype_for_R8plot",
      "celltype",
      "auto_celltype_final",
      "auto_celltype_safe",
      "celltype_auto_annotation",
      "auto_celltype"
    ),
    override, FALSE, "annotation column"
  )
}

resolve_cluster_column <- function(object, override = NULL) {
  cols <- colnames(object[[]])

  exact_priority <- c(
    "cluster_for_R8plot_FIXED2",
    "cluster_for_R8plot",
    "cluster_for_annotation",
    "integratedRPCA_snn_res.3",
    "integratedRPCA_snn_res.3.0",
    "RNA_snn_res.3",
    "RNA_snn_res.3.0",
    "integratedRPCA_snn_res.0.8",
    "RNA_snn_res.0.8",
    "seurat_clusters"
  )

  resolved <- resolve_from_priority(
    cols, exact_priority, override, FALSE, "cluster column"
  )

  if (!is.na(resolved)) return(resolved)

  hits <- grep("(snn_res|cluster)", cols, ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(hits[[1]])

  NA_character_
}

natural_level_order <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values)]
  numeric_values <- suppressWarnings(as.numeric(values))

  if (length(values) > 0 && all(!is.na(numeric_values))) {
    return(values[order(numeric_values)])
  }

  sort(values)
}

prepare_metadata_factor <- function(object, column) {
  if (is.na(column) || !column %in% colnames(object[[]])) return(object)

  x <- object[[column, drop = TRUE]]
  if (!is.factor(x)) {
    object[[column]] <- factor(x, levels = natural_level_order(x))
  }

  object
}

safe_umap_filename <- function(x) {
  x <- gsub("\\.[Rr][Dd][Ss]$", "", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

make_R8_palette <- function(n) {
  if (n <= 0) return(character(0))
  hues <- (seq_len(n) - 1) * 137.508 %% 360
  luminance <- rep(c(50, 62, 42, 56), length.out = n)
  grDevices::hcl(h = hues, c = 100, l = luminance, fixup = TRUE)
}

make_pastel_palette <- function(n) {
  if (n <= 0) return(character(0))

  hues <- (seq_len(n) - 1) * 137.508 %% 360

  grDevices::hcl(
    h = hues,
    c = 58,
    l = 74,
    fixup = TRUE
  )
}

get_group_palette <- function(object, group_by, palette_style = "R8") {
  x <- object[[group_by, drop = TRUE]]
  lv <- if (is.factor(x)) levels(x) else natural_level_order(x)

  cols <- switch(
    palette_style,
    R8 = make_R8_palette(length(lv)),
    Pastel = make_pastel_palette(length(lv)),
    stop("Unknown palette_style: ", palette_style)
  )

  names(cols) <- lv
  cols
}

build_umap_footer <- function(
    rds_file, framework_name, framework_version,
    include_rds = TRUE, include_created = TRUE,
    include_framework = TRUE) {

  parts <- character(0)

  if (isTRUE(include_rds)) {
    parts <- c(parts, paste0("RDS: ", rds_file))
  }

  if (isTRUE(include_created)) {
    parts <- c(parts, paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M")))
  }

  if (isTRUE(include_framework)) {
    parts <- c(parts, paste0(framework_name, " v", framework_version))
  }

  paste(parts, collapse = "    |    ")
}

extract_umap_dataframe <- function(
    object, reduction, metadata_columns = character(0)) {

  coords <- as.data.frame(Seurat::Embeddings(object, reduction = reduction))
  coords <- coords[, 1:2, drop = FALSE]
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords$cell <- rownames(coords)

  metadata_columns <- unique(metadata_columns[!is.na(metadata_columns)])
  metadata_columns <- metadata_columns[
    metadata_columns %in% colnames(object[[]])
  ]

  if (length(metadata_columns) > 0) {
    meta <- object[[]][rownames(coords), metadata_columns, drop = FALSE]
    coords <- cbind(coords, meta)
  }

  coords
}

choose_split_ncol <- function(n_panels, max_columns = 3) {
  if (n_panels == 4) return(2L)
  if (n_panels == 6) return(3L)
  as.integer(min(max_columns, ceiling(sqrt(n_panels))))
}

calculate_split_dimensions <- function(
    n_panels, panel_width = 5.2,
    panel_height = 4.8, max_columns = 3) {

  ncol <- choose_split_ncol(n_panels, max_columns)
  nrow <- ceiling(n_panels / ncol)

  c(
    width = panel_width * ncol,
    height = panel_height * nrow,
    ncol = ncol,
    nrow = nrow
  )
}

add_footer_annotation <- function(plot, footer_text, footer_text_size = 6.5) {
  if (is.null(footer_text) || !nzchar(footer_text)) return(plot)

  plot +
    patchwork::plot_annotation(
      caption = footer_text,
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(
          size = footer_text_size,
          colour = "grey30",
          hjust = 0,
          margin = ggplot2::margin(t = 8)
        )
      )
    )
}

make_dimplot <- function(
    object, reduction, group_by, palette_style = "R8",
    title = NULL, point_size = 0.40,
    raster = FALSE, raster_dpi = c(512, 512),
    shuffle = TRUE, label = FALSE, label_size = 4.8,
    hide_legend = FALSE, legend_point_size = 4.2,
    legend_text_size = 10.5, title_size = 14,
    axis_text_size = 10, footer_text = NULL,
    footer_text_size = 6.5) {

  cols <- get_group_palette(object, group_by, palette_style)

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = point_size,
    raster = raster,
    raster.dpi = raster_dpi,
    shuffle = shuffle,
    label = label,
    label.size = label_size,
    repel = label,
    cols = unname(cols)
  ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size, face = "bold", hjust = 0.5
      ),
      axis.text = ggplot2::element_text(
        size = axis_text_size, colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = axis_text_size + 1, colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.55, colour = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text_size),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = legend_point_size, alpha = 1)
      )
    )

  if (isTRUE(hide_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  add_footer_annotation(p, footer_text, footer_text_size)
}

make_faceted_split_umap <- function(
    umap_df, group_column, split_column,
    palette_style = "R8", title,
    point_size = 0.40, max_columns = 3,
    legend_point_size = 4.2,
    legend_text_size = 10.5,
    title_size = 14, axis_text_size = 9,
    footer_text = NULL, footer_text_size = 6.5) {

  umap_df$.group <- factor(
    umap_df[[group_column]],
    levels = natural_level_order(umap_df[[group_column]])
  )

  umap_df$.split <- factor(
    umap_df[[split_column]],
    levels = natural_level_order(umap_df[[split_column]])
  )

  palette <- switch(
    palette_style,
    R8 = make_R8_palette(nlevels(umap_df$.group)),
    Pastel = make_pastel_palette(nlevels(umap_df$.group)),
    stop("Unknown palette_style: ", palette_style)
  )
  names(palette) <- levels(umap_df$.group)

  ncol <- choose_split_ncol(nlevels(umap_df$.split), max_columns)

  p <- ggplot2::ggplot(
    umap_df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = .group)
  ) +
    ggplot2::geom_point(size = point_size, alpha = 1, stroke = 0) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::facet_wrap(
      ggplot2::vars(.split),
      ncol = ncol,
      scales = "fixed"
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title,
      x = "UMAP 1",
      y = "UMAP 2",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size, face = "bold", hjust = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "white", colour = "black", linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        size = 11, face = "bold", colour = "black"
      ),
      axis.text = ggplot2::element_text(
        size = axis_text_size, colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = axis_text_size + 1, colour = "black"
      ),
      panel.spacing = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(size = legend_text_size)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = legend_point_size, alpha = 1)
      )
    )

  add_footer_annotation(p, footer_text, footer_text_size)
}

make_single_sample_monochrome_umap <- function(
    umap_df, sample_column, sample_value,
    title, point_color = "#0057B8",
    point_size = 0.45, global_limits = TRUE,
    transparent_background = TRUE,
    title_size = 14, axis_text_size = 10,
    footer_text = NULL, footer_text_size = 6.5) {

  subset_df <- umap_df[
    as.character(umap_df[[sample_column]]) == as.character(sample_value),
    ,
    drop = FALSE
  ]

  bg_element <- if (isTRUE(transparent_background)) {
    ggplot2::element_rect(fill = "transparent", colour = NA)
  } else {
    ggplot2::element_rect(fill = "white", colour = NA)
  }

  p <- ggplot2::ggplot(
    subset_df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2)
  ) +
    ggplot2::geom_point(
      colour = point_color,
      size = point_size,
      alpha = 1,
      stroke = 0
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Cells: ", format(nrow(subset_df), big.mark = ",")),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.background = bg_element,
      panel.background = bg_element,
      legend.background = bg_element,
      legend.box.background = bg_element,
      plot.title = ggplot2::element_text(
        size = title_size, face = "bold", hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10, hjust = 0.5, colour = "grey30"
      ),
      axis.text = ggplot2::element_text(
        size = axis_text_size, colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = axis_text_size + 1, colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.55, colour = "black"
      )
    )

  if (isTRUE(global_limits)) {
    p <- p +
      ggplot2::scale_x_continuous(
        limits = range(umap_df$UMAP_1, na.rm = TRUE),
        expand = ggplot2::expansion(mult = 0.03)
      ) +
      ggplot2::scale_y_continuous(
        limits = range(umap_df$UMAP_2, na.rm = TRUE),
        expand = ggplot2::expansion(mult = 0.03)
      )
  }

  add_footer_annotation(p, footer_text, footer_text_size)
}

save_umap_pdf <- function(
    plot, path, width, height,
    overwrite = FALSE, transparent = FALSE) {

  if (file.exists(path) && !isTRUE(overwrite)) {
    message("Skipped existing file: ", path)
    return(invisible(FALSE))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(transparent) && capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = path,
      width = width,
      height = height,
      bg = "transparent",
      onefile = TRUE
    )
  } else {
    grDevices::pdf(
      file = path,
      width = width,
      height = height,
      bg = if (isTRUE(transparent)) "transparent" else "white",
      useDingbats = FALSE,
      onefile = TRUE
    )
  }

  print(plot)
  grDevices::dev.off()

  invisible(TRUE)
}
