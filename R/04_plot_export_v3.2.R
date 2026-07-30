save_pdf <- function(plot, path, width, height) {
  grDevices::cairo_pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

publish_umap_v32 <- function(object, group_by, cfg, title = NULL, label = cfg$umap_label) {
  Seurat::DimPlot(
    object,
    reduction = cfg$reduction,
    group.by = group_by,
    label = label,
    repel = TRUE,
    pt.size = cfg$umap_point_size,
    raster = TRUE
  ) +
    ggplot2::ggtitle(title %||% group_by) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.key.height = grid::unit(0.55, "cm")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1),
        ncol = 1
      )
    )
}

make_marker_dotplot <- function(object, marker_reference, group_by, cfg, source = "General") {
  genes <- unique(marker_reference$gene[
    marker_reference$source == source &
    marker_reference$direction == "positive"
  ])
  genes <- intersect(genes, rownames(object))
  Seurat::DotPlot(
    object,
    features = genes,
    group.by = group_by,
    assay = cfg$assay,
    cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
    dot.scale = 7
  ) +
    ggplot2::scale_color_gradient2(
      low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A",
      midpoint = 0,
      limits = c(cfg$dotplot_scale_min, cfg$dotplot_scale_max),
      oob = scales::squish
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_line(linewidth = 0.15)
    )
}

make_marker_violin_pages <- function(object, marker_reference, group_by, cfg, source = "Ueno") {
  genes <- unique(marker_reference$gene[
    marker_reference$source == source &
    marker_reference$direction == "positive"
  ])
  genes <- intersect(genes, rownames(object))
  chunks <- split(genes, ceiling(seq_along(genes) / 12))
  lapply(chunks, function(g) {
    Seurat::VlnPlot(
      object,
      features = g,
      group.by = group_by,
      assay = cfg$assay,
      layer = cfg$data_layer,
      pt.size = 0,
      ncol = cfg$violin_ncol
    )
  })
}

write_excel_report <- function(path, sheets) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    warning("openxlsx is unavailable; Excel report was not written.")
    return(invisible(NULL))
  }
  wb <- openxlsx::createWorkbook()
  used <- character()
  for (nm in names(sheets)) {
    sn <- substr(gsub("[\\\\/:*?\\[\\]]", "_", nm), 1, 31)
    if (sn %in% used) sn <- substr(paste0(sn, "_", length(used) + 1L), 1, 31)
    used <- c(used, sn)
    openxlsx::addWorksheet(wb, sn)
    x <- sheets[[nm]]
    if (is.matrix(x)) x <- data.frame(gene = rownames(x), x, check.names = FALSE)
    openxlsx::writeDataTable(wb, sn, x, withFilter = TRUE)
    openxlsx::freezePane(wb, sn, firstRow = TRUE)
    openxlsx::setColWidths(wb, sn, cols = seq_len(ncol(x)), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}
