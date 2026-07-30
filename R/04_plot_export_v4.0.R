# R8-style visualization and explanatory UMAP outputs for uenoy scRNAseq Framework v4.0.
# This file depends on Seurat, ggplot2, ggrepel and patchwork.

`%||%` <- function(x, y) if (is.null(x)) y else x

ueno_lineage_palette_v40 <- c(
  Biliary = "#F4A261",
  Cycling = "#B89B00",
  Endothelial = "#43B7B1",
  Erythroid = "#C77DBB",
  Hepatocyte = "#A8C957",
  Lymphoid = "#4EA8DE",
  Megakaryocytic = "#9D79BC",
  Mesenchymal = "#7BC96F",
  Myeloid = "#59C3B0",
  Unknown = "#BDBDBD"
)

ueno_confidence_palette_v40 <- c(
  High = "#E64B35",
  Moderate = "#3C5488",
  Low = "#00A087",
  Unknown = "#BDBDBD",
  `NA` = "#D9D9D9"
)

ueno_celltype_palette_v40 <- c(
  Mature_hepatocyte = "#B8CE63",
  Hepatic_progenitor = "#D4E157",
  Cholangiocyte = "#F4A261",
  LSEC = "#42C2B3",
  Vascular_endothelial = "#69B3D0",
  Lymphatic_endothelial = "#7AD7D0",
  qHSC = "#9BD77A",
  aHSC = "#74C476",
  Portal_fibroblast = "#C2A55F",
  Pericyte_VSMC = "#8A7D55",
  Kupffer_macrophage = "#45B8A4",
  Monocyte_derived_macrophage = "#2A9D8F",
  Monocyte = "#57C5B6",
  Neutrophil = "#89C2D9",
  cDC1 = "#B9A44C",
  cDC2 = "#D3B85A",
  pDC = "#7B9ACC",
  Mast_cell = "#B565A7",
  Basophil = "#E08DAC",
  B_cell = "#F28E8B",
  Plasma_cell = "#9A78C4",
  CD4_T_cell = "#FF8FAB",
  CD8_T_cell = "#F6BD60",
  Treg = "#D45087",
  Gamma_delta_T = "#9B5DE5",
  NK_cell = "#70A5D8",
  NKT_cell = "#5B8ECA",
  ILC = "#8EC5A4",
  Erythroid = "#C77DBB",
  Megakaryocyte = "#A07CC5",
  Platelet = "#C3A6D9",
  Cycling = "#B89B00",
  Unknown = "#BDBDBD"
)

ueno_subtype_palette_v40 <- c(
  Periportal_hepatocyte = "#C7D87A",
  Midzonal_hepatocyte = "#B5C957",
  Pericentral_hepatocyte = "#819B3A",
  Stress_response_hepatocyte = "#D6C667",
  IFN_response_hepatocyte = "#A3B86C",
  Cycling_hepatocyte = "#A68A00",
  Reactive_cholangiocyte = "#F6B26B",
  Periportal_LSEC = "#72D6C9",
  Pericentral_LSEC = "#48B8B0",
  Capillarized_LSEC = "#2F8F8A",
  Angiogenic_LSEC = "#4FA3C7",
  Inflammatory_LSEC = "#3D7EA6",
  Early_activated_HSC = "#A6D98A",
  Myofibroblastic_aHSC = "#6BBF74",
  Fibrogenic_aHSC = "#4E9E5A",
  Inflammatory_aHSC = "#8FCF86",
  Resident_Kupffer_like = "#56C5B3",
  Monocyte_like = "#2EAA9A",
  Inflammatory_M1_like = "#D95F76",
  Pro_resolution_M2_like = "#65B891",
  SPP1_TREM2_MASH_associated = "#C05A9D",
  Lipid_associated_macrophage = "#8F5FA8",
  Efferocytosis_phagocytosis_high = "#5DAE8B",
  IL10_response_high_Mphi = "#4E9F78",
  Fibrosis_associated_Mphi = "#A24A73",
  IFN_response_macrophage = "#5470B3",
  Classical_monocyte = "#4BC0B5",
  Nonclassical_monocyte = "#3B9B91",
  Inflammatory_neutrophil = "#91C4DA",
  Aged_neutrophil = "#6DA7C2",
  Naive_CD4_T = "#F8A1B7",
  Activated_CD4_T = "#E86A92",
  Naive_CD8_T = "#F6CB67",
  Cytotoxic_CD8_T = "#DDAA33",
  Exhausted_CD8_T = "#C68E17",
  Activated_NK = "#6EA8D7",
  Naive_B = "#F1948A",
  Memory_B = "#D97872",
  Plasma_cell_subtype = "#9C7BC7",
  qHSC = "#9BD77A",
  Erythroid = "#C77DBB",
  Vascular_endothelial = "#69B3D0",
  cDC1 = "#B9A44C",
  cDC2 = "#D3B85A",
  pDC = "#7B9ACC",
  Basophil = "#E08DAC",
  Cycling = "#B89B00",
  Unknown = "#BDBDBD"
)

resolve_palette_v40 <- function(values, palette, fallback = "#BDBDBD") {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  missing <- setdiff(values, names(palette))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dynamic")
    names(extra) <- missing
    palette <- c(palette, extra)
  }
  palette[values] %||% setNames(rep(fallback, length(values)), values)
}

make_two_line_label_v40 <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- gsub("SPP1 TREM2", "SPP1/TREM2", x)
  x <- gsub("MASH associated", "MASH-associated", x)
  x <- gsub("Monocyte derived", "Monocyte-derived", x)
  x <- gsub("Resident Kupffer like", "Resident Kupffer-like", x)
  x <- gsub("Periportal hepatocyte", "Periportal\nhepatocyte", x)
  x <- gsub("Pericentral hepatocyte", "Pericentral\nhepatocyte", x)
  x <- gsub("Midzonal hepatocyte", "Midzonal\nhepatocyte", x)
  x <- gsub("Reactive cholangiocyte", "Reactive\ncholangiocyte", x)
  x <- gsub("Capillarized LSEC", "Capillarized\nLSEC", x)
  x <- gsub("Pericentral LSEC", "Pericentral\nLSEC", x)
  x <- gsub("Periportal LSEC", "Periportal\nLSEC", x)
  x <- gsub("Myofibroblastic aHSC", "Myofibroblastic\naHSC", x)
  x <- gsub("Fibrogenic aHSC", "Fibrogenic\naHSC", x)
  x <- gsub("Inflammatory neutrophil", "Inflammatory\nneutrophil", x)
  x <- gsub("Cytotoxic CD8 T", "Cytotoxic\nCD8 T", x)
  x <- gsub("Naive CD8 T", "Naive\nCD8 T", x)
  x <- gsub("Activated NK", "Activated\nNK", x)
  x
}

publish_umap_r8_v40 <- function(
    object,
    group_by,
    reduction,
    palette,
    title = NULL,
    pt_size = 0.22,
    label = TRUE,
    repel = TRUE,
    label_size = 3.4,
    legend_ncol = 1L,
    raster = TRUE
) {
  if (!group_by %in% colnames(object[[]])) {
    stop("Metadata field not found: ", group_by)
  }
  vals <- as.character(object[[]][[group_by]])
  vals[is.na(vals)] <- "Unknown"
  object[[group_by]][, 1] <- factor(vals)

  pal <- resolve_palette_v40(levels(object[[group_by]][,1]), palette)

  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = pt_size,
    label = FALSE,
    raster = raster
  ) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = title %||% group_by, color = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8.5),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.42, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3.2, alpha = 1),
        ncol = legend_ncol
      )
    )

  if (isTRUE(label)) {
    emb <- Seurat::Embeddings(object, reduction = reduction)
    meta <- object[[]]
    dat <- data.frame(
      x = emb[,1],
      y = emb[,2],
      label = as.character(meta[[group_by]]),
      stringsAsFactors = FALSE
    )
    dat$label[is.na(dat$label)] <- "Unknown"

    centers <- stats::aggregate(
      dat[c("x","y")],
      by = list(label = dat$label),
      FUN = stats::median
    )
    centers$display <- make_two_line_label_v40(centers$label)

    if (isTRUE(repel)) {
      p <- p + ggrepel::geom_text_repel(
        data = centers,
        ggplot2::aes(x = x, y = y, label = display),
        inherit.aes = FALSE,
        size = label_size,
        fontface = "plain",
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

publish_confidence_umap_v40 <- function(
    object, group_by, reduction, title = NULL, pt_size = 0.22
) {
  publish_umap_r8_v40(
    object = object,
    group_by = group_by,
    reduction = reduction,
    palette = ueno_confidence_palette_v40,
    title = title,
    pt_size = pt_size,
    label = TRUE,
    repel = TRUE,
    label_size = 3.6
  )
}

derive_summary_annotation_v40 <- function(object) {
  meta <- object[[]]
  required <- c(
    "vote_ueno_lineage_v34",
    "vote_ueno_celltype_v34",
    "vote_ueno_subtype_v34",
    "vote_ueno_lineage_v34_confidence",
    "vote_ueno_celltype_v34_confidence",
    "vote_ueno_subtype_v34_confidence"
  )
  missing <- setdiff(required, colnames(meta))
  if (length(missing)) {
    stop("Required v3.4 metadata missing: ", paste(missing, collapse = ", "))
  }

  lineage <- as.character(meta$vote_ueno_lineage_v34)
  celltype <- as.character(meta$vote_ueno_celltype_v34)
  subtype <- as.character(meta$vote_ueno_subtype_v34)
  lin_conf <- as.character(meta$vote_ueno_lineage_v34_confidence)
  ct_conf <- as.character(meta$vote_ueno_celltype_v34_confidence)
  st_conf <- as.character(meta$vote_ueno_subtype_v34_confidence)

  summary <- celltype

  use_subtype <- !is.na(subtype) &
    subtype != "Unknown" &
    st_conf %in% c("High","Moderate")

  summary[use_subtype] <- subtype[use_subtype]

  use_celltype <- !use_subtype &
    !is.na(celltype) &
    celltype != "Unknown" &
    ct_conf %in% c("High","Moderate")

  summary[use_celltype] <- celltype[use_celltype]

  fallback <- (!use_subtype & !use_celltype) |
    is.na(summary) | summary == "Unknown"

  summary[fallback] <- lineage[fallback]
  summary[is.na(summary)] <- "Unknown"

  object$vote_ueno_summary_v40 <- summary
  object
}

derive_annotation_difference_v40 <- function(
    object,
    current_col = "celltype_for_R8plot_FIXED2",
    recommended_col = "vote_ueno_summary_v40"
) {
  meta <- object[[]]
  if (!all(c(current_col, recommended_col) %in% colnames(meta))) {
    stop("Difference columns not found.")
  }

  current <- as.character(meta[[current_col]])
  recommended <- as.character(meta[[recommended_col]])

  normalize <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]", "", x)
    x <- gsub("macrophage", "mphi", x)
    x
  }

  status <- ifelse(
    is.na(recommended) | recommended == "Unknown",
    "Review",
    ifelse(normalize(current) == normalize(recommended), "Concordant", "Changed")
  )

  object$annotation_difference_v40 <- status
  object
}

annotation_difference_palette_v40 <- c(
  Concordant = "#C7C7C7",
  Changed = "#E64B35",
  Review = "#3C5488"
)

make_cluster_audit_pages_v40 <- function(
    object,
    audit_table,
    cluster_col,
    reduction,
    current_col = "celltype_for_R8plot_FIXED2"
) {
  emb <- Seurat::Embeddings(object, reduction = reduction)
  meta <- object[[]]
  meta$cluster_key_v40 <- ifelse(
    startsWith(as.character(meta[[cluster_col]]), "g"),
    as.character(meta[[cluster_col]]),
    paste0("g", as.character(meta[[cluster_col]]))
  )

  clusters <- unique(audit_table$cluster)
  pages <- vector("list", length(clusters))
  names(pages) <- clusters

  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    selected <- meta$cluster_key_v40 == cl

    dat <- data.frame(
      x = emb[,1],
      y = emb[,2],
      selected = ifelse(selected, "Target cluster", "Background")
    )

    p_umap <- ggplot2::ggplot(dat, ggplot2::aes(x, y, color = selected)) +
      ggplot2::geom_point(size = 0.18, alpha = 0.9) +
      ggplot2::scale_color_manual(
        values = c("Background" = "#D9D9D9", "Target cluster" = "#D73027")
      ) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::labs(title = paste("Cluster", sub("^g","",cl)), color = NULL) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )

    rows <- audit_table[audit_table$cluster == cl, , drop = FALSE]
    rows <- rows[match(c("lineage","celltype","subtype"), rows$level), , drop = FALSE]

    fmt <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "—", as.character(x))

    text_lines <- c(
      paste0("Current annotation: ", fmt(rows$current_annotation[1])),
      "",
      paste0("Lineage: ", fmt(rows$recommended_label[rows$level == "lineage"])),
      paste0("Cell type: ", fmt(rows$recommended_label[rows$level == "celltype"])),
      paste0("Subtype: ", fmt(rows$recommended_label[rows$level == "subtype"])),
      "",
      paste0("Confidence: ",
             paste(rows$level, rows$confidence, sep = "=", collapse = "; ")),
      "",
      "Supporting markers:",
      paste(fmt(rows$supporting_markers), collapse = "\n"),
      "",
      "Contradicting markers:",
      paste(fmt(rows$contradicting_markers), collapse = "\n")
    )

    text_df <- data.frame(x = 0, y = 1, label = paste(text_lines, collapse = "\n"))

    p_text <- ggplot2::ggplot(text_df, ggplot2::aes(x, y, label = label)) +
      ggplot2::geom_text(hjust = 0, vjust = 1, size = 3.2, lineheight = 1.05) +
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
      ggplot2::theme_void()

    pages[[i]] <- p_umap + p_text + patchwork::plot_layout(widths = c(1.2, 1))
  }

  pages
}
