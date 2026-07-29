# ============================================================
# UMAP resolvers and publication plotting
# uenoy scRNAseq Framework v2.1.1
# Split UMAP aspect-ratio fix
# ============================================================

resolve_from_priority <- function(
    available,
    priority,
    override = NULL,
    required = FALSE,
    label = "column") {

  if (!is.null(override)) {
    if (!override %in% available) {
      stop("Configured ", label, " was not found: ", override)
    }
    return(override)
  }

  hit <- priority[priority %in% available]

  if (length(hit) > 0) {
    return(hit[[1]])
  }

  if (isTRUE(required)) {
    stop("No suitable ", label, " was found.")
  }

  NA_character_
}

resolve_umap_reduction <- function(object, override = NULL) {
  resolve_from_priority(
    available = names(object@reductions),
    priority = c(
      "umapRPCA", "umap.rpca", "integrated_umap",
      "umap_integrated", "umap"
    ),
    override = override,
    required = TRUE,
    label = "UMAP reduction"
  )
}

resolve_sample_column <- function(object, override = NULL) {
  resolve_from_priority(
    available = colnames(object[[]]),
    priority = c(
      "sample", "sample_for_annotation", "sample_id",
      "orig.ident", "replicate"
    ),
    override = override,
    required = FALSE,
    label = "sample column"
  )
}

resolve_condition_column <- function(object, override = NULL) {
  resolve_from_priority(
    available = colnames(object[[]]),
    priority = c(
      "condition", "group", "treatment",
      "disease_status", "status", "diet"
    ),
    override = override,
    required = FALSE,
    label = "condition column"
  )
}

resolve_annotation_column <- function(object, override = NULL) {
  resolve_from_priority(
    available = colnames(object[[]]),
    priority = c(
      "celltype_for_R8plot_FIXED2",
      "celltype_for_R8plot",
      "celltype",
      "auto_celltype_final",
      "auto_celltype_safe",
      "celltype_auto_annotation",
      "auto_celltype"
    ),
    override = override,
    required = FALSE,
    label = "annotation column"
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
    available = cols,
    priority = exact_priority,
    override = override,
    required = FALSE,
    label = "cluster column"
  )

  if (!is.na(resolved)) return(resolved)

  hits <- grep(
    "(snn_res|cluster)",
    cols,
    ignore.case = TRUE,
    value = TRUE
  )

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
  if (is.na(column) || !column %in% colnames(object[[]])) {
    return(object)
  }

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

choose_split_ncol <- function(n_panels, max_columns = 3) {
  if (is.na(n_panels) || n_panels < 1) return(1L)

  # Near-square layout:
  # 2 panels -> 2 x 1
  # 4 panels -> 2 x 2
  # 6 panels -> 3 x 2
  as.integer(min(max_columns, ceiling(sqrt(n_panels))))
}

calculate_split_dimensions <- function(
    n_panels,
    panel_width = 5,
    panel_height = 4.5,
    max_columns = 3) {

  ncol <- choose_split_ncol(n_panels, max_columns)
  nrow <- ceiling(n_panels / ncol)

  c(
    width = panel_width * ncol,
    height = panel_height * nrow,
    ncol = ncol,
    nrow = nrow
  )
}

make_vivid_palette <- function(n) {
  if (n <= 0) return(character(0))

  hues <- (seq_len(n) - 1) * 137.508 %% 360
  luminance <- rep(c(52, 66, 42, 60), length.out = n)

  grDevices::hcl(
    h = hues,
    c = 95,
    l = luminance,
    fixup = TRUE
  )
}

get_group_palette <- function(object, group_by) {
  x <- object[[group_by, drop = TRUE]]
  lv <- if (is.factor(x)) levels(x) else natural_level_order(x)

  cols <- make_vivid_palette(length(lv))
  names(cols) <- lv
  cols
}

build_umap_footer <- function(
    rds_file,
    framework_name,
    framework_version,
    include_rds = TRUE,
    include_created = TRUE,
    include_framework = TRUE) {

  parts <- character(0)

  if (isTRUE(include_rds)) {
    parts <- c(parts, paste0("RDS: ", rds_file))
  }

  if (isTRUE(include_created)) {
    parts <- c(
      parts,
      paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
    )
  }

  if (isTRUE(include_framework)) {
    parts <- c(
      parts,
      paste0(framework_name, " v", framework_version)
    )
  }

  paste(parts, collapse = "    |    ")
}

apply_umap_theme <- function(
    p,
    title_size = 14,
    axis_text_size = 10,
    legend_point_size = 4.2,
    legend_text_size = 10.5,
    hide_legend = FALSE) {

  p <- p +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      ),
      axis.text = ggplot2::element_text(
        size = axis_text_size,
        colour = "black"
      ),
      axis.title = ggplot2::element_text(
        size = axis_text_size + 1,
        colour = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.55,
        colour = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = legend_text_size,
        colour = "black"
      ),
      legend.key.height = grid::unit(0.52, "cm"),
      legend.key.width = grid::unit(0.52, "cm"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(
          size = legend_point_size,
          alpha = 1
        )
      )
    )

  if (isTRUE(hide_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}

make_dimplot <- function(
    object,
    reduction,
    group_by,
    split_by = NULL,
    split_ncol = NULL,
    title = NULL,
    point_size = 0.40,
    raster = FALSE,
    raster_dpi = c(512, 512),
    shuffle = TRUE,
    label = FALSE,
    label_size = 4.8,
    hide_legend = FALSE,
    legend_point_size = 4.2,
    legend_text_size = 10.5,
    title_size = 14,
    axis_text_size = 10,
    footer_text = NULL,
    footer_text_size = 6.5) {

  cols <- get_group_palette(object, group_by)

  if (is.null(split_by)) {
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
    )

    p <- apply_umap_theme(
      p = p,
      title_size = title_size,
      axis_text_size = axis_text_size,
      legend_point_size = legend_point_size,
      legend_text_size = legend_text_size,
      hide_legend = hide_legend
    )

    p <- p + ggplot2::labs(
      title = title,
      x = "UMAP 1",
      y = "UMAP 2"
    )

  } else {
    # Seurat's default split.by output is one horizontal row.
    # combine = FALSE allows explicit multi-row patchwork layout.
    plot_list <- Seurat::DimPlot(
      object = object,
      reduction = reduction,
      group.by = group_by,
      split.by = split_by,
      pt.size = point_size,
      raster = raster,
      raster.dpi = raster_dpi,
      shuffle = shuffle,
      label = label,
      label.size = label_size,
      repel = label,
      cols = unname(cols),
      combine = FALSE
    )

    plot_list <- lapply(plot_list, function(q) {
      q <- apply_umap_theme(
        p = q,
        title_size = title_size - 1,
        axis_text_size = axis_text_size,
        legend_point_size = legend_point_size,
        legend_text_size = legend_text_size,
        hide_legend = hide_legend
      )

      q + ggplot2::labs(x = "UMAP 1", y = "UMAP 2")
    })

    if (is.null(split_ncol)) {
      split_ncol <- choose_split_ncol(length(plot_list), max_columns = 3)
    }

    p <- patchwork::wrap_plots(
      plot_list,
      ncol = split_ncol,
      guides = "collect"
    ) +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = title_size,
            face = "bold",
            hjust = 0.5,
            margin = ggplot2::margin(b = 8)
          )
        )
      )

    if (!isTRUE(hide_legend)) {
      p <- p & ggplot2::theme(legend.position = "right")
    }
  }

  if (!is.null(footer_text) && nzchar(footer_text)) {
    p <- p +
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

  p
}

save_umap_pdf <- function(
    plot,
    path,
    width,
    height,
    overwrite = FALSE) {

  if (file.exists(path) && !isTRUE(overwrite)) {
    message("Skipped existing file: ", path)
    return(invisible(FALSE))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(
    file = path,
    width = width,
    height = height,
    useDingbats = FALSE,
    onefile = TRUE
  )

  print(plot)
  grDevices::dev.off()

  invisible(TRUE)
}
