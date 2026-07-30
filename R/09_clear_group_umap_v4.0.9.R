# ============================================================
# uenoyscRNA Framework v4.0.8
# Clearer UMAPs with larger points and labels
# Integrated group-level UMAPs
# ============================================================

derive_analysis_group_v408 <- function(
    object,
    sample_col,
    output_col = "analysis_group_v408"
) {
  metadata <- object[[]]

  if (!sample_col %in% colnames(metadata)) {
    stop("Sample metadata column not found: ", sample_col)
  }

  sample_value <- as.character(metadata[[sample_col]])
  group_value <- rep("Other", length(sample_value))

  group_value[grepl("STD", sample_value, ignore.case = TRUE)] <- "STD"
  group_value[grepl("CDHFD|CDAHFD", sample_value, ignore.case = TRUE)] <- "CDHFD"
  group_value[grepl("Sham", sample_value, ignore.case = TRUE)] <- "Sham"
  group_value[grepl("^Tx|Tx17|Tx5", sample_value, ignore.case = TRUE)] <- "Tx"

  preferred_order <- c("STD", "CDHFD", "Sham", "Tx", "Other")
  present_order <- preferred_order[preferred_order %in% unique(group_value)]

  object[[output_col]] <- factor(
    group_value,
    levels = present_order
  )

  object
}

publish_split_umap_clear_v408 <- function(
    object,
    group_by,
    split_by,
    reduction,
    palette,
    title,
    pt_size = 0.68,
    raster = TRUE,
    raster_dpi = c(600, 600),
    source_rds = NULL,
    created_at = Sys.time(),
    ncol = NULL,
    legend_text_size = 9,
    strip_text_size = 11
) {
  metadata <- object[[]]

  missing_cols <- setdiff(c(group_by, split_by), colnames(metadata))
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
    fallback = "#777777"
  )

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    pt.size = pt_size,
    raster = raster,
    raster.dpi = raster_dpi,
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
      title = title,
      color = NULL,
      caption = make_caption_v402(
        source_rds = source_rds,
        created_at = created_at
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = strip_text_size
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7.5,
        colour = "grey25"
      ),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.height = grid::unit(0.44, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4.5, alpha = 1)
      )
    )
}

publish_single_group_umap_v408 <- function(
    object,
    group_value,
    analysis_group_col = "analysis_group_v408",
    annotation_col = "vote_ueno_summary_v40",
    reduction,
    palette,
    pt_size = 0.82,
    label_size = 4.0,
    source_rds = NULL,
    created_at = Sys.time()
) {
  metadata <- object[[]]

  if (!analysis_group_col %in% colnames(metadata)) {
    stop("Analysis-group column not found: ", analysis_group_col)
  }

  keep_cells <- rownames(metadata)[
    as.character(metadata[[analysis_group_col]]) == group_value
  ]

  if (length(keep_cells) == 0L) {
    stop("No cells found for group: ", group_value)
  }

  sub_object <- subset(object, cells = keep_cells)

  publish_umap_r8_v402(
    object = sub_object,
    group_by = annotation_col,
    reduction = reduction,
    palette = palette,
    title = paste0("RDS3: integrated ", group_value, " group"),
    pt_size = pt_size,
    label = TRUE,
    repel = TRUE,
    label_size = label_size,
    legend_ncol = 1L,
    raster = TRUE,
    source_rds = source_rds,
    created_at = created_at
  )
}

save_each_integrated_group_umap_v408 <- function(
    object,
    output_pdf,
    analysis_group_col = "analysis_group_v408",
    annotation_col = "vote_ueno_summary_v40",
    reduction,
    palette,
    pt_size = 0.82,
    label_size = 4.0,
    source_rds = NULL,
    created_at = Sys.time()
) {
  group_levels <- levels(object[[analysis_group_col]][, 1])

  grDevices::pdf(
    output_pdf,
    width = 12,
    height = 11,
    onefile = TRUE,
    useDingbats = FALSE
  )

  on.exit(grDevices::dev.off(), add = TRUE)

  for (group_value in group_levels) {
    p <- publish_single_group_umap_v408(
      object = object,
      group_value = group_value,
      analysis_group_col = analysis_group_col,
      annotation_col = annotation_col,
      reduction = reduction,
      palette = palette,
      pt_size = pt_size,
      label_size = label_size,
      source_rds = source_rds,
      created_at = created_at
    )
    print(p)
  }

  invisible(output_pdf)
}
