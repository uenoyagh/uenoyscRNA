# ============================================================
# uenoyscRNA Framework v4.0.3
# Split UMAP and top-3 marker violin plots
# ============================================================

sanitize_filename_v403 <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "Unknown")
}

publish_split_umap_v403 <- function(
    object,
    group_by,
    split_by,
    reduction,
    palette,
    title = NULL,
    pt_size = 0.28,
    raster = TRUE,
    source_rds = NULL,
    created_at = Sys.time(),
    ncol = NULL
) {
  metadata <- object[[]]

  needed <- c(group_by, split_by)
  missing_cols <- setdiff(needed, colnames(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!reduction %in% names(object@reductions)) {
    stop("Reduction not found: ", reduction)
  }

  group_values <- as.character(metadata[[group_by]])
  group_values[is.na(group_values) | !nzchar(group_values)] <- "Unknown"
  object[[group_by]] <- factor(group_values, levels = sort(unique(group_values)))

  split_values <- as.character(metadata[[split_by]])
  split_values[is.na(split_values) | !nzchar(split_values)] <- "Unknown"
  object[[split_by]] <- factor(split_values, levels = unique(split_values))

  active_levels <- levels(object[[group_by]][, 1])
  plot_palette <- resolve_palette_v402(
    values = active_levels,
    palette = palette,
    fallback = "#7F7F7F"
  )

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    pt.size = pt_size,
    raster = raster,
    ncol = ncol
  )

  if (length(p$layers) >= 1L) {
    for (i in seq_along(p$layers)) {
      p$layers[[i]]$aes_params$alpha <- 1
    }
  }

  p +
    ggplot2::scale_color_manual(values = plot_palette, drop = FALSE) +
    ggplot2::labs(
      title = title %||% paste0(group_by, " split by ", split_by),
      color = NULL,
      caption = make_caption_v402(
        source_rds = source_rds,
        created_at = created_at
      )
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 9
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7,
        colour = "grey25"
      ),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.height = grid::unit(0.4, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3.2, alpha = 1)
      )
    )
}

find_top3_markers_by_celltype_v403 <- function(
    object,
    identity_col = "vote_ueno_summary_v40",
    assay = NULL,
    slot = "data",
    min_pct = 0.20,
    logfc_threshold = 0.25,
    only_pos = TRUE,
    top_n = 3L
) {
  if (!identity_col %in% colnames(object[[]])) {
    stop("Metadata column not found: ", identity_col)
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  Seurat::Idents(object) <- identity_col

  markers <- Seurat::FindAllMarkers(
    object = object,
    assay = assay,
    slot = slot,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )

  if (nrow(markers) == 0L) {
    stop("No marker genes were detected.")
  }

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC", "avg_diff"),
    colnames(markers)
  )

  if (length(fc_col) == 0L) {
    stop("Fold-change column was not found in marker results.")
  }

  fc_col <- fc_col[[1]]

  markers <- markers[
    !is.na(markers[[fc_col]]) &
      !is.na(markers$p_val_adj),
    ,
    drop = FALSE
  ]

  markers <- markers[
    order(
      markers$cluster,
      -markers[[fc_col]],
      markers$p_val_adj
    ),
    ,
    drop = FALSE
  ]

  top_markers <- do.call(
    rbind,
    lapply(
      split(markers, markers$cluster),
      function(df) utils::head(df, n = top_n)
    )
  )

  rownames(top_markers) <- NULL
  attr(top_markers, "fc_col") <- fc_col
  top_markers
}

save_top3_violin_pdf_v403 <- function(
    object,
    top_markers,
    output_pdf,
    identity_col = "vote_ueno_summary_v40",
    assay = NULL,
    slot = "data",
    source_rds = NULL,
    created_at = Sys.time(),
    width = 13,
    height = 8.5
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  if (!identity_col %in% colnames(object[[]])) {
    stop("Metadata column not found: ", identity_col)
  }

  if (!all(c("gene", "cluster") %in% colnames(top_markers))) {
    stop("top_markers must contain gene and cluster columns.")
  }

  grDevices::pdf(
    output_pdf,
    width = width,
    height = height,
    onefile = TRUE,
    useDingbats = FALSE
  )

  on.exit(grDevices::dev.off(), add = TRUE)

  clusters <- unique(as.character(top_markers$cluster))

  for (cell_type in clusters) {
    genes <- unique(
      as.character(
        top_markers$gene[top_markers$cluster == cell_type]
      )
    )
    genes <- genes[genes %in% rownames(object)]

    if (length(genes) == 0L) {
      next
    }

    plots <- Seurat::VlnPlot(
      object = object,
      features = genes,
      group.by = identity_col,
      assay = assay,
      slot = slot,
      pt.size = 0,
      combine = FALSE
    )

    plots <- lapply(
      seq_along(plots),
      function(i) {
        plots[[i]] +
          ggplot2::geom_boxplot(
            width = 0.12,
            outlier.shape = NA,
            alpha = 0.35
          ) +
          ggplot2::labs(
            title = genes[[i]],
            x = NULL,
            y = "Normalized expression"
          ) +
          ggplot2::theme_classic(base_size = 10) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "bold",
              size = 12
            ),
            axis.text.x = ggplot2::element_text(
              angle = 60,
              hjust = 1,
              vjust = 1,
              size = 7
            ),
            legend.position = "none"
          )
      }
    )

    combined <- patchwork::wrap_plots(
      plots,
      ncol = length(plots)
    ) +
      patchwork::plot_annotation(
        title = paste0(
          "Top marker genes for ",
          gsub("_", " ", cell_type, fixed = TRUE)
        ),
        subtitle = paste0(
          "Cell type of interest: ",
          gsub("_", " ", cell_type, fixed = TRUE),
          " | top ",
          length(genes),
          " genes"
        ),
        caption = make_caption_v402(
          source_rds = source_rds,
          created_at = created_at
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold",
            size = 15
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5,
            size = 10
          ),
          plot.caption = ggplot2::element_text(
            hjust = 0,
            size = 7,
            colour = "grey25"
          )
        )
      )

    print(combined)
  }

  invisible(output_pdf)
}

save_each_celltype_violin_v403 <- function(
    object,
    top_markers,
    output_dir,
    identity_col = "vote_ueno_summary_v40",
    assay = NULL,
    slot = "data",
    source_rds = NULL,
    created_at = Sys.time()
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  clusters <- unique(as.character(top_markers$cluster))

  for (cell_type in clusters) {
    genes <- unique(
      as.character(
        top_markers$gene[top_markers$cluster == cell_type]
      )
    )
    genes <- genes[genes %in% rownames(object)]

    if (length(genes) == 0L) {
      next
    }

    plots <- Seurat::VlnPlot(
      object = object,
      features = genes,
      group.by = identity_col,
      assay = assay,
      slot = slot,
      pt.size = 0,
      combine = FALSE
    )

    plots <- lapply(
      seq_along(plots),
      function(i) {
        plots[[i]] +
          ggplot2::geom_boxplot(
            width = 0.12,
            outlier.shape = NA,
            alpha = 0.35
          ) +
          ggplot2::labs(
            title = genes[[i]],
            x = NULL,
            y = "Normalized expression"
          ) +
          ggplot2::theme_classic(base_size = 10) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "bold"
            ),
            axis.text.x = ggplot2::element_text(
              angle = 60,
              hjust = 1,
              vjust = 1,
              size = 7
            ),
            legend.position = "none"
          )
      }
    )

    combined <- patchwork::wrap_plots(
      plots,
      ncol = length(plots)
    ) +
      patchwork::plot_annotation(
        title = paste0(
          "Top marker genes for ",
          gsub("_", " ", cell_type, fixed = TRUE)
        ),
        caption = make_caption_v402(
          source_rds = source_rds,
          created_at = created_at
        )
      )

    file_name <- paste0(
      "Violin_top3_",
      sanitize_filename_v403(cell_type),
      ".pdf"
    )

    save_pdf(
      combined,
      file.path(output_dir, file_name),
      width = 13,
      height = 8.5
    )
  }

  invisible(output_dir)
}
