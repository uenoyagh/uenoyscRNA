# ============================================================
# uenoy scRNAseq Framework v4.0.2
# High-contrast R8 UMAP plotting
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

ueno_lineage_palette_v402 <- c(
  Biliary = "#D96F12",
  Cycling = "#9A8000",
  Endothelial = "#0B8F8F",
  Erythroid = "#A63D91",
  Hepatocyte = "#789A1F",
  Lymphoid = "#2E6FAF",
  Megakaryocytic = "#7B4FA7",
  Mesenchymal = "#4E963F",
  Myeloid = "#087A68",
  Unknown = "#7F7F7F"
)

ueno_confidence_palette_v402 <- c(
  High = "#D73027",
  Moderate = "#2C5AA0",
  Low = "#00876C",
  Unknown = "#7F7F7F",
  `NA` = "#BDBDBD"
)

ueno_celltype_palette_v402 <- c(
  Mature_hepatocyte = "#789A1F",
  Hepatic_progenitor = "#A5B520",
  Cholangiocyte = "#D96F12",
  LSEC = "#0B8F8F",
  Vascular_endothelial = "#247FA6",
  Lymphatic_endothelial = "#2E9E99",
  qHSC = "#5D9F45",
  aHSC = "#3C8F3E",
  Portal_fibroblast = "#8A7029",
  Pericyte_VSMC = "#665A2F",
  Kupffer_macrophage = "#087A68",
  Monocyte_derived_macrophage = "#005F56",
  Monocyte = "#149E8D",
  Neutrophil = "#3B86AD",
  cDC1 = "#826F20",
  cDC2 = "#A98620",
  pDC = "#5069AA",
  Mast_cell = "#843C7B",
  Basophil = "#B45186",
  B_cell = "#D9585A",
  Plasma_cell = "#6C4AA3",
  CD4_T_cell = "#D54771",
  CD8_T_cell = "#C98B16",
  Treg = "#96305E",
  Gamma_delta_T = "#6538A6",
  NK_cell = "#3A74B5",
  NKT_cell = "#2C5B93",
  ILC = "#4B865F",
  Erythroid = "#A63D91",
  Megakaryocyte = "#7B4FA7",
  Platelet = "#9B78BD",
  Cycling = "#9A8000",
  Unknown = "#7F7F7F"
)

ueno_subtype_palette_v402 <- c(
  Mature_hepatocyte = "#789A1F",
  Periportal_hepatocyte = "#9DB12F",
  Midzonal_hepatocyte = "#789A1F",
  Pericentral_hepatocyte = "#526C13",
  Stress_response_hepatocyte = "#A68B22",
  IFN_response_hepatocyte = "#627C3A",
  Cycling_hepatocyte = "#8A7100",
  Biliary = "#D96F12",
  Cholangiocyte = "#D96F12",
  Reactive_cholangiocyte = "#C85E08",
  LSEC = "#0B8F8F",
  Vascular_endothelial = "#247FA6",
  Periportal_LSEC = "#219D97",
  Pericentral_LSEC = "#0B8785",
  Capillarized_LSEC = "#006E6D",
  Angiogenic_LSEC = "#247FA6",
  Inflammatory_LSEC = "#1C608A",
  Mesenchymal = "#4E963F",
  qHSC = "#5D9F45",
  aHSC = "#3C8F3E",
  Early_activated_HSC = "#679F4E",
  Myofibroblastic_aHSC = "#3C8F3E",
  Fibrogenic_aHSC = "#237431",
  Inflammatory_aHSC = "#4E934A",
  Portal_fibroblast = "#8A7029",
  Pericyte_VSMC = "#665A2F",
  Kupffer_macrophage = "#087A68",
  Monocyte_derived_macrophage = "#005F56",
  Resident_Kupffer_like = "#118A73",
  Monocyte_like = "#006F62",
  Inflammatory_M1_like = "#BE3B54",
  Pro_resolution_M2_like = "#3E8060",
  SPP1_TREM2_MASH_associated = "#8D2F79",
  Lipid_associated_macrophage = "#613886",
  Efferocytosis_phagocytosis_high = "#347456",
  IL10_response_high_Mphi = "#286348",
  Fibrosis_associated_Mphi = "#75304F",
  IFN_response_macrophage = "#354E96",
  Classical_monocyte = "#0F907F",
  Nonclassical_monocyte = "#006F65",
  Neutrophil = "#3B86AD",
  Inflammatory_neutrophil = "#3B86AD",
  Aged_neutrophil = "#276B91",
  cDC1 = "#826F20",
  cDC2 = "#A98620",
  pDC = "#5069AA",
  Lymphoid = "#2E6FAF",
  B_cell = "#D9585A",
  Naive_B = "#D66460",
  Memory_B = "#B94646",
  Plasma_cell = "#6C4AA3",
  Plasma_cell_subtype = "#6C4AA3",
  CD4_T_cell = "#D54771",
  Naive_CD4_T = "#DA6D8C",
  Activated_CD4_T = "#C53B68",
  CD8_T_cell = "#C98B16",
  Naive_CD8_T = "#C99A26",
  Cytotoxic_CD8_T = "#A9730F",
  Exhausted_CD8_T = "#805203",
  NK_cell = "#3A74B5",
  Activated_NK = "#3A74B5",
  NKT_cell = "#2C5B93",
  Treg = "#96305E",
  Gamma_delta_T = "#6538A6",
  Basophil = "#B45186",
  Mast_cell = "#843C7B",
  Erythroid = "#A63D91",
  Cycling = "#9A8000",
  Unknown = "#7F7F7F"
)

annotation_difference_palette_v402 <- c(
  Concordant = "#AFAFAF",
  Changed = "#D73027",
  Review = "#2C5AA0",
  Unknown = "#7F7F7F"
)

resolve_palette_v402 <- function(values, palette, fallback = "#7F7F7F") {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]

  missing_values <- setdiff(values, names(palette))

  if (length(missing_values) > 0L) {
    additional_colors <- grDevices::hcl.colors(
      n = length(missing_values),
      palette = "Dark 3"
    )
    names(additional_colors) <- missing_values
    palette <- c(palette, additional_colors)
  }

  result <- palette[values]
  result[is.na(result)] <- fallback
  result
}

make_caption_v402 <- function(source_rds = NULL, created_at = Sys.time()) {
  source_name <- if (
    is.null(source_rds) ||
    length(source_rds) == 0L ||
    is.na(source_rds) ||
    !nzchar(source_rds)
  ) {
    "Not specified"
  } else {
    basename(source_rds)
  }

  paste0(
    "Source RDS: ", source_name,
    "    |    Created: ",
    format(created_at, "%Y-%m-%d %H:%M:%S %Z")
  )
}

calculate_umap_label_centres_v402 <- function(object, group_by, reduction) {
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  metadata <- object[[]]

  label_data <- data.frame(
    x = embedding[, 1],
    y = embedding[, 2],
    label = as.character(metadata[[group_by]]),
    stringsAsFactors = FALSE
  )

  label_data$label[is.na(label_data$label)] <- "Unknown"

  centres <- stats::aggregate(
    label_data[c("x", "y")],
    by = list(label = label_data$label),
    FUN = stats::median
  )

  centres$display <- gsub("_", " ", centres$label, fixed = TRUE)
  centres
}

publish_umap_r8_v402 <- function(
    object,
    group_by,
    reduction,
    palette,
    title = NULL,
    pt_size = 0.44,
    label = TRUE,
    repel = TRUE,
    label_size = 3.4,
    legend_ncol = 1L,
    raster = TRUE,
    source_rds = NULL,
    created_at = Sys.time()
) {
  metadata <- object[[]]

  if (!group_by %in% colnames(metadata)) {
    stop("Metadata field not found: ", group_by)
  }

  if (!reduction %in% names(object@reductions)) {
    stop("Reduction not found: ", reduction)
  }

  group_values <- as.character(metadata[[group_by]])
  group_values[is.na(group_values) | !nzchar(group_values)] <- "Unknown"

  object[[group_by]] <- factor(
    group_values,
    levels = sort(unique(group_values))
  )

  active_levels <- levels(object[[group_by]][, 1])

  plot_palette <- resolve_palette_v402(
    values = active_levels,
    palette = palette,
    fallback = "#7F7F7F"
  )

  caption_text <- make_caption_v402(
    source_rds = source_rds,
    created_at = created_at
  )

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = pt_size,
    label = FALSE,
    raster = raster
  )

  if (length(p$layers) >= 1L) {
    p$layers[[1]]$aes_params$alpha <- 1
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = plot_palette,
      drop = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title %||% group_by,
      x = paste0(reduction, "_1"),
      y = paste0(reduction, "_2"),
      color = NULL,
      caption = caption_text
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7.5,
        colour = "grey25",
        margin = ggplot2::margin(t = 8)
      ),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8.5),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.42, "cm"),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3.4, alpha = 1),
        ncol = legend_ncol,
        byrow = TRUE
      )
    )

  if (isTRUE(label)) {
    centres <- calculate_umap_label_centres_v402(
      object = object,
      group_by = group_by,
      reduction = reduction
    )

    if (isTRUE(repel)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = centres,
          mapping = ggplot2::aes(x = x, y = y, label = display),
          inherit.aes = FALSE,
          size = label_size,
          colour = "black",
          min.segment.length = 0,
          max.overlaps = Inf,
          box.padding = 0.35,
          point.padding = 0.15,
          segment.size = 0.25,
          segment.alpha = 0.7,
          seed = 1234
        )
    } else {
      p <- p +
        ggplot2::geom_text(
          data = centres,
          mapping = ggplot2::aes(x = x, y = y, label = display),
          inherit.aes = FALSE,
          size = label_size,
          colour = "black"
        )
    }
  }

  p
}

publish_confidence_umap_v402 <- function(
    object,
    group_by,
    reduction,
    title = NULL,
    pt_size = 0.44,
    source_rds = NULL,
    created_at = Sys.time()
) {
  publish_umap_r8_v402(
    object = object,
    group_by = group_by,
    reduction = reduction,
    palette = ueno_confidence_palette_v402,
    title = title,
    pt_size = pt_size,
    label = TRUE,
    repel = TRUE,
    label_size = 3.6,
    legend_ncol = 1L,
    raster = TRUE,
    source_rds = source_rds,
    created_at = created_at
  )
}
