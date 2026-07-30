# ============================================================
# uenoyscRNA Framework v4.1.0
# Dense group UMAP + single-page target-cell violin matrix
# ============================================================

target_cell_order_v410 <- c(
  "Hepatocyte",
  "LSEC",
  "Vascular_endothelial",
  "Kupffer_M2_Mphi",
  "Monocyte_M1_Mphi",
  "SPP1_TREM2_Mphi",
  "qHSC",
  "aHSC",
  "Mesenchymal",
  "Biliary",
  "B_cell",
  "T_cell",
  "NK_cell",
  "Plasma_cell"
)

target_cell_palette_v410 <- c(
  Hepatocyte = "#009688",
  LSEC = "#00AEEF",
  Vascular_endothelial = "#B86B00",
  Kupffer_M2_Mphi = "#2F8F2F",
  Monocyte_M1_Mphi = "#00A67A",
  SPP1_TREM2_Mphi = "#B02A9B",
  qHSC = "#F34BA6",
  aHSC = "#C2187A",
  Mesenchymal = "#E950A7",
  Biliary = "#00B84A",
  B_cell = "#E83E8C",
  T_cell = "#5B5CE2",
  NK_cell = "#C89D00",
  Plasma_cell = "#2F8FD8"
)

preferred_marker_sets_v410 <- list(
  Hepatocyte = c("Alb", "Ttr", "Apoa1", "Cps1", "Ass1", "Fabp1"),
  LSEC = c("Clec4g", "Stab2", "Fcgr2b", "Lyve1", "Kdr"),
  Vascular_endothelial = c("Vwf", "Emcn", "Esam", "Pecam1", "Kdr"),
  Kupffer_M2_Mphi = c("C1qa", "C1qb", "C1qc", "Marco", "Clec4f", "Timd4", "Mrc1"),
  Monocyte_M1_Mphi = c("Lyz2", "Ccr2", "S100a8", "S100a9", "Il1b", "Ccl2", "Tyrobp"),
  SPP1_TREM2_Mphi = c("Spp1", "Trem2", "Gpnmb", "Lgals3", "Cd9"),
  qHSC = c("Lrat", "Rbp1", "Reln", "Dcn", "Col15a1"),
  aHSC = c("Col1a1", "Col1a2", "Acta2", "Tagln", "Timp1", "Lox"),
  Mesenchymal = c("Dcn", "Lum", "Col3a1", "Pdgfra", "Col1a1"),
  Biliary = c("Krt19", "Krt8", "Krt18", "Sox9", "Epcam", "Krt7"),
  B_cell = c("Cd79a", "Ms4a1", "Cd74", "Cd37", "H2-Aa"),
  T_cell = c("Cd3d", "Cd3e", "Trbc1", "Lck", "Tcf7"),
  NK_cell = c("Nkg7", "Klrd1", "Prf1", "Gzmb", "Xcl1"),
  Plasma_cell = c("Jchain", "Mzb1", "Sdc1", "Xbp1", "Igkc")
)

derive_target_cellclass_v410 <- function(
    object,
    source_col = "vote_ueno_summary_v40",
    output_col = "target_cellclass_v410"
) {
  metadata <- object[[]]

  if (!source_col %in% colnames(metadata)) {
    stop("Source annotation column not found: ", source_col)
  }

  x <- as.character(metadata[[source_col]])
  out <- rep(NA_character_, length(x))

  out[grepl("Mature_hepatocyte|Periportal_hepatocyte|Pericentral_hepatocyte|Midzonal_hepatocyte",
            x, ignore.case = TRUE)] <- "Hepatocyte"

  out[grepl("^LSEC$|Capillarized_LSEC|Periportal_LSEC|Pericentral_LSEC",
            x, ignore.case = TRUE)] <- "LSEC"

  out[grepl("Vascular_endothelial", x, ignore.case = TRUE)] <- "Vascular_endothelial"

  out[grepl("Kupffer_macrophage|Resident_Kupffer_like|Pro_resolution_M2_like|IL10_response_high_Mphi",
            x, ignore.case = TRUE)] <- "Kupffer_M2_Mphi"

  out[grepl("Monocyte_derived_macrophage|Monocyte_like|Inflammatory_M1_like|Classical_monocyte|Nonclassical_monocyte",
            x, ignore.case = TRUE)] <- "Monocyte_M1_Mphi"

  out[grepl("SPP1_TREM2_MASH_associated|Lipid_associated_macrophage|Fibrosis_associated_Mphi",
            x, ignore.case = TRUE)] <- "SPP1_TREM2_Mphi"

  out[grepl("^qHSC$", x, ignore.case = TRUE)] <- "qHSC"

  out[grepl("Myofibroblastic_aHSC|Fibrogenic_aHSC|Early_activated_HSC|Inflammatory_aHSC|^aHSC$",
            x, ignore.case = TRUE)] <- "aHSC"

  out[grepl("^Mesenchymal$|Portal_fibroblast|Pericyte_VSMC",
            x, ignore.case = TRUE)] <- "Mesenchymal"

  out[grepl("Biliary|Cholangiocyte|Reactive_cholangiocyte",
            x, ignore.case = TRUE)] <- "Biliary"

  out[grepl("^B_cell$|Naive_B|Memory_B", x, ignore.case = TRUE)] <- "B_cell"

  out[grepl("Naive_CD8_T|Cytotoxic_CD8_T|CD4_T_cell|CD8_T_cell|Treg|Gamma_delta_T|^Lymphoid$",
            x, ignore.case = TRUE)] <- "T_cell"

  out[grepl("^NK_cell$|Activated_NK|NKT_cell", x, ignore.case = TRUE)] <- "NK_cell"

  out[grepl("Plasma_cell", x, ignore.case = TRUE)] <- "Plasma_cell"

  present_levels <- target_cell_order_v410[
    target_cell_order_v410 %in% unique(out[!is.na(out)])
  ]

  object[[output_col]] <- factor(out, levels = present_levels)
  object
}

resolve_features_case_insensitive_v410 <- function(features, available_features) {
  available_lower <- tolower(available_features)
  resolved <- character(0)

  for (feature in features) {
    idx <- match(tolower(feature), available_lower)
    if (!is.na(idx)) {
      resolved <- c(resolved, available_features[[idx]])
    }
  }

  unique(resolved)
}

find_target_markers_v410 <- function(
    object,
    identity_col = "target_cellclass_v410",
    assay = NULL,
    slot = "data",
    min_pct = 0.15,
    logfc_threshold = 0.15
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  metadata <- object[[]]
  keep_cells <- rownames(metadata)[!is.na(metadata[[identity_col]])]
  target_object <- subset(object, cells = keep_cells)

  Seurat::Idents(target_object) <- identity_col

  markers <- Seurat::FindAllMarkers(
    object = target_object,
    assay = assay,
    slot = slot,
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )

  if (nrow(markers) == 0L) {
    stop("No marker genes were detected for target cell classes.")
  }

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC", "avg_diff"),
    colnames(markers)
  )

  if (length(fc_col) == 0L) {
    stop("Fold-change column was not found.")
  }

  fc_col <- fc_col[[1]]

  if (!"pct.1" %in% colnames(markers)) markers$pct.1 <- 1
  if (!"pct.2" %in% colnames(markers)) markers$pct.2 <- 0
  if (!"p_val_adj" %in% colnames(markers)) markers$p_val_adj <- 1

  markers <- markers[
    !grepl(
      "^(MT-|RPL|RPS|HBA|HBB)|^MALAT1$|^XIST$|^JUN$|^FOS$|^FOSB$|^IER2$|^IER3$",
      markers$gene,
      ignore.case = FALSE
    ),
    ,
    drop = FALSE
  ]

  markers$specificity_delta <- pmax(markers$pct.1 - markers$pct.2, 0)
  markers$characteristic_score <- (
    pmax(markers[[fc_col]], 0) *
      pmax(markers$pct.1, 0.01) *
      pmax(markers$specificity_delta, 0.01)
  )

  selected_rows <- list()
  available_features <- rownames(object)

  for (cell_class in levels(target_object[[identity_col]][, 1])) {
    class_markers <- markers[
      as.character(markers$cluster) == cell_class,
      ,
      drop = FALSE
    ]

    if (nrow(class_markers) == 0L) next

    preferred <- resolve_features_case_insensitive_v410(
      preferred_marker_sets_v410[[cell_class]],
      available_features
    )

    preferred_rows <- class_markers[
      tolower(class_markers$gene) %in% tolower(preferred),
      ,
      drop = FALSE
    ]

    candidate_rows <- if (nrow(preferred_rows) > 0L) {
      preferred_rows
    } else {
      class_markers
    }

    candidate_rows <- candidate_rows[
      order(
        -candidate_rows$characteristic_score,
        -candidate_rows[[fc_col]],
        candidate_rows$p_val_adj
      ),
      ,
      drop = FALSE
    ]

    selected_rows[[cell_class]] <- candidate_rows[1, , drop = FALSE]
  }

  selected <- do.call(rbind, selected_rows)
  rownames(selected) <- NULL
  selected$target_cellclass <- as.character(selected$cluster)
  selected$display_cellclass <- gsub("_", " ", selected$target_cellclass)
  attr(selected, "fc_col") <- fc_col
  selected
}

fetch_violin_data_v410 <- function(
    object,
    genes,
    identity_col,
    analysis_group_col,
    assay = NULL,
    slot = "data"
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  vars <- unique(c(genes, identity_col, analysis_group_col))

  data <- Seurat::FetchData(
    object = object,
    vars = vars,
    layer = slot
  )

  data$cell_id <- rownames(data)
  data
}

make_violin_panel_v410 <- function(
    data,
    gene,
    analysis_group,
    cell_order,
    palette,
    show_x_text = FALSE,
    y_max = NULL,
    show_y_title = FALSE
) {
  plot_data <- data[
    as.character(data$analysis_group_v408) == analysis_group &
      !is.na(data$target_cellclass_v410),
    ,
    drop = FALSE
  ]

  plot_data$target_cellclass_v410 <- factor(
    as.character(plot_data$target_cellclass_v410),
    levels = cell_order
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = target_cellclass_v410,
      y = .data[[gene]],
      fill = target_cellclass_v410
    )
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.22,
      colour = "black",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = palette[cell_order],
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, y_max),
      clip = "off"
    ) +
    ggplot2::labs(
      title = gene,
      x = NULL,
      y = if (show_y_title) "Expression" else NULL
    ) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 9
      ),
      axis.title.y = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 6.5),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(1, 2, 1, 2)
    )

  if (show_x_text) {
    p <- p +
      ggplot2::scale_x_discrete(
        labels = function(x) gsub("_", " ", x)
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 55,
          hjust = 1,
          vjust = 1,
          size = 6.5
        )
      )
  } else {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }

  p
}

save_singlepage_target_violin_v410 <- function(
    object,
    marker_table,
    output_pdf,
    identity_col = "target_cellclass_v410",
    analysis_group_col = "analysis_group_v408",
    assay = NULL,
    slot = "data",
    source_rds = NULL,
    created_at = Sys.time(),
    group_order = c("STD", "CDHFD", "Sham", "Tx"),
    width = 22,
    height = 12.5
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  cell_order <- target_cell_order_v410[
    target_cell_order_v410 %in%
      levels(object[[identity_col]][, 1])
  ]

  marker_table <- marker_table[
    match(cell_order, marker_table$target_cellclass),
    ,
    drop = FALSE
  ]

  genes <- as.character(marker_table$gene)

  violin_data <- fetch_violin_data_v410(
    object = object,
    genes = genes,
    identity_col = identity_col,
    analysis_group_col = analysis_group_col,
    assay = assay,
    slot = slot
  )

  colnames(violin_data)[
    colnames(violin_data) == identity_col
  ] <- "target_cellclass_v410"

  colnames(violin_data)[
    colnames(violin_data) == analysis_group_col
  ] <- "analysis_group_v408"

  panels <- list()

  for (i in seq_along(genes)) {
    gene <- genes[[i]]
    gene_values <- violin_data[[gene]]
    gene_values <- gene_values[is.finite(gene_values)]

    y_max <- if (length(gene_values) > 0L) {
      as.numeric(stats::quantile(gene_values, 0.995, na.rm = TRUE))
    } else {
      1
    }

    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    y_max <- y_max * 1.08

    for (j in seq_along(group_order)) {
      group_value <- group_order[[j]]

      panels[[length(panels) + 1L]] <- make_violin_panel_v410(
        data = violin_data,
        gene = gene,
        analysis_group = group_value,
        cell_order = cell_order,
        palette = target_cell_palette_v410,
        show_x_text = i == length(genes),
        y_max = y_max,
        show_y_title = j == 1L
      )
    }
  }

  column_headers <- lapply(
    group_order,
    function(group_value) {
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = group_value,
          fontface = "bold",
          size = 5.5
        ) +
        ggplot2::xlim(0, 1) +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
    }
  )

  header <- patchwork::wrap_plots(column_headers, nrow = 1)

  body <- patchwork::wrap_plots(
    panels,
    ncol = length(group_order),
    byrow = TRUE
  )

  final_plot <- header / body +
    patchwork::plot_layout(
      heights = c(0.35, length(genes))
    ) +
    patchwork::plot_annotation(
      title = "One characteristic marker per selected liver cell class",
      subtitle = paste0(
        "Rows: marker genes | Columns: integrated analysis groups | ",
        "x-axis: selected cell classes"
      ),
      caption = make_caption_v402(
        source_rds = source_rds,
        created_at = created_at
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 18
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

  grDevices::pdf(
    output_pdf,
    width = width,
    height = height,
    onefile = TRUE,
    useDingbats = FALSE
  )

  print(final_plot)
  grDevices::dev.off()

  invisible(output_pdf)
}
