# uenoy scRNAseq Framework v4.0.1
# High-contrast R8 UMAP plotting with source-RDS and creation-time captions.

`%||%` <- function(x, y) if (is.null(x)) y else x

ueno_lineage_palette_v401 <- c(
  Biliary = "#F28E2B",
  Cycling = "#B59B00",
  Endothelial = "#18A6A6",
  Erythroid = "#B85AA6",
  Hepatocyte = "#8FAF2F",
  Lymphoid = "#2F7FC1",
  Megakaryocytic = "#8E63B8",
  Mesenchymal = "#59A84F",
  Myeloid = "#169C84",
  Unknown = "#8F8F8F"
)

ueno_confidence_palette_v401 <- c(
  High = "#D73027",
  Moderate = "#2C5AA0",
  Low = "#00876C",
  Unknown = "#8F8F8F",
  `NA` = "#BDBDBD"
)

ueno_celltype_palette_v401 <- c(
  Mature_hepatocyte = "#9BBE3B",
  Hepatic_progenitor = "#C3CD2E",
  Cholangiocyte = "#F28E2B",
  LSEC = "#18A6A6",
  Vascular_endothelial = "#2C8FB8",
  Lymphatic_endothelial = "#46BEB5",
  qHSC = "#7EBE63",
  aHSC = "#4FAE52",
  Portal_fibroblast = "#A9893A",
  Pericyte_VSMC = "#7A6D3A",
  Kupffer_macrophage = "#179C83",
  Monocyte_derived_macrophage = "#087F71",
  Monocyte = "#28B7A3",
  Neutrophil = "#4C9BC0",
  cDC1 = "#9A8732",
  cDC2 = "#C29D32",
  pDC = "#617CBF",
  Mast_cell = "#9C4F93",
  Basophil = "#C96C9B",
  B_cell = "#E76F70",
  Plasma_cell = "#805CB6",
  CD4_T_cell = "#E65B83",
  CD8_T_cell = "#E4A62A",
  Treg = "#B23D73",
  Gamma_delta_T = "#7646B8",
  NK_cell = "#4C86C6",
  NKT_cell = "#376DA8",
  ILC = "#5B9E75",
  Erythroid = "#B85AA6",
  Megakaryocyte = "#8E63B8",
  Platelet = "#AD8CCB",
  Cycling = "#B59B00",
  Unknown = "#8F8F8F"
)

ueno_subtype_palette_v401 <- c(
  Periportal_hepatocyte = "#B4C64F",
  Midzonal_hepatocyte = "#94AF32",
  Pericentral_hepatocyte = "#66851E",
  Stress_response_hepatocyte = "#C0A832",
  IFN_response_hepatocyte = "#7D9850",
  Cycling_hepatocyte = "#A48700",
  Reactive_cholangiocyte = "#E9781B",
  Periportal_LSEC = "#37B9AF",
  Pericentral_LSEC = "#169C98",
  Capillarized_LSEC = "#087B79",
  Angiogenic_LSEC = "#2C8FB8",
  Inflammatory_LSEC = "#246E98",
  Early_activated_HSC = "#77B75F",
  Myofibroblastic_aHSC = "#4FAE52",
  Fibrogenic_aHSC = "#2F8C3D",
  Inflammatory_aHSC = "#62A95B",
  Resident_Kupffer_like = "#25A98F",
  Monocyte_like = "#0C8E7D",
  Inflammatory_M1_like = "#CF4B64",
  Pro_resolution_M2_like = "#4C9B73",
  SPP1_TREM2_MASH_associated = "#A83D8E",
  Lipid_associated_macrophage = "#74499A",
  Efferocytosis_phagocytosis_high = "#438C69",
  IL10_response_high_Mphi = "#347B5A",
  Fibrosis_associated_Mphi = "#8A385F",
  IFN_response_macrophage = "#405DA6",
  Classical_monocyte = "#1FAF9C",
  Nonclassical_monocyte = "#0C8679",
  Inflammatory_neutrophil = "#4C9BC0",
  Aged_neutrophil = "#347DA2",
  Naive_CD4_T = "#ED819D",
  Activated_CD4_T = "#D94D78",
  Naive_CD8_T = "#E7B33B",
  Cytotoxic_CD8_T = "#C78B16",
  Exhausted_CD8_T = "#9F6808",
  Activated_NK = "#4C86C6",
  Naive_B = "#E87B78",
  Memory_B = "#CB5A58",
  Plasma_cell_subtype = "#805CB6",
  qHSC = "#7EBE63",
  Erythroid = "#B85AA6",
  Vascular_endothelial = "#2C8FB8",
  cDC1 = "#9A8732",
  cDC2 = "#C29D32",
  pDC = "#617CBF",
  Basophil = "#C96C9B",
  Cycling = "#B59B00",
  Unknown = "#8F8F8F"
)

resolve_palette_v401 <- function(values, palette, fallback = "#8F8F8F") {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  missing <- setdiff(values, names(palette))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    palette <- c(palette, extra)
  }
  palette[values] %||% setNames(rep(fallback, length(values)), values)
}

make_caption_v401 <- function(source_rds = NULL, created_at = Sys.time()) {
  source_name <- if (is.null(source_rds) || !nzchar(source_rds)) {
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

publish_umap_r8_v401 <- function(
    object,
    group_by,
    reduction,
    palette,
    title = NULL,
    pt_size = 0.34,
    label = TRUE,
    repel = TRUE,
    label_size = 3.4,
    legend_ncol = 1L,
    raster = TRUE,
    source_rds = NULL,
    created_at = Sys.time()
) {
  if (!group_by %in% colnames(object[[]])) {
    stop("Metadata field not found: ", group_by)
  }

  vals <- as.character(object[[]][[group_by]])
  vals[is.na(vals)] <- "Unknown"
  object[[group_by]][, 1] <- factor(vals)

  pal <- resolve_palette_v401(levels(object[[group_by]][, 1]), palette)
  caption_text <- make_caption_v401(source_rds, created_at)

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = pt_size,
    label = FALSE,
    raster = raster
  )

  # Force fully opaque points, including rasterized output.
  if (length(p$layers) >= 1L) {
    p$layers[[1]]$aes_params$alpha <- 1
  }

  p <- p +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title %||% group_by,
      color = NULL,
      caption = caption_text
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 8,
        colour = "grey30",
        margin = ggplot2::margin(t = 8)
      ),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8.5),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.42, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3.3, alpha = 1),
        ncol = legend_ncol
      )
    )

  if (isTRUE(label)) {
    emb <- Seurat::Embeddings(object, reduction = reduction)
    meta <- object[[]]

    dat <- data.frame(
      x = emb[, 1],
      y = emb[, 2],
      label = as.character(meta[[group_by]]),
      stringsAsFactors = FALSE
    )
    dat$label[is.na(dat$label)] <- "Unknown"

    centers <- stats::aggregate(
      dat[c("x", "y")],
      by = list(label = dat$label),
      FUN = stats::median
    )

    # Keep full information-rich names; replace underscores only on-map.
    centers$display <- gsub("_", " ", centers$label)

    if (isTRUE(repel)) {
      p <- p + ggrepel::geom_text_repel(
        data = centers,
        ggplot2::aes(x = x, y = y, label = display),
        inherit.aes = FALSE,
        size = label_size,
        min.segment.length = 0,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.15,
        segment.size = 0.25,
        seed = 1234
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = centers,
        ggplot2::aes(x = x, y = y, label = display),
        inherit.aes = FALSE,
        size = label_size
      )
    }
  }

  p
}

publish_confidence_umap_v401 <- function(
    object,
    group_by,
    reduction,
    title = NULL,
    pt_size = 0.34,
    source_rds = NULL,
    created_at = Sys.time()
) {
  publish_umap_r8_v401(
    object = object,
    group_by = group_by,
    reduction = reduction,
    palette = ueno_confidence_palette_v401,
    title = title,
    pt_size = pt_size,
    label = TRUE,
    repel = TRUE,
    label_size = 3.6,
    source_rds = source_rds,
    created_at = created_at
  )
}
