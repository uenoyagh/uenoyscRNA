# ============================================================
# uenoyscRNA Framework v4.0.7
# Saturated UMAP palette + one characteristic marker per cell type
# Grouped violin panels by lineage and condition
# ============================================================

sanitize_filename_v407 <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "Unknown")
}

# ------------------------------------------------------------
# High-saturation palette modeled after the user's preferred UMAP
# ------------------------------------------------------------

ueno_lineage_palette_v407 <- c(
  Biliary      = "#00B84A",
  Cycling      = "#8A2BE2",
  Endothelial  = "#00AEEF",
  Erythroid    = "#008C2E",
  Hepatocyte   = "#009688",
  Lymphoid     = "#3D5AFE",
  Mesenchymal  = "#F34BA6",
  Myeloid      = "#00A67A",
  Unknown      = "#777777"
)

ueno_celltype_palette_v407 <- c(
  Mature_hepatocyte            = "#009688",
  Hepatic_progenitor           = "#00B39B",
  Cholangiocyte                = "#00B84A",
  LSEC                         = "#00AEEF",
  Vascular_endothelial         = "#B86B00",
  Lymphatic_endothelial        = "#00B7B7",
  qHSC                         = "#F34BA6",
  aHSC                         = "#D62C8A",
  Portal_fibroblast            = "#C4589C",
  Pericyte_VSMC                = "#8B5A2B",
  Kupffer_macrophage           = "#2F8F2F",
  Monocyte_derived_macrophage  = "#00A67A",
  Monocyte                     = "#00C08B",
  Neutrophil                   = "#D44A9A",
  cDC1                         = "#C57A00",
  cDC2                         = "#A45A00",
  pDC                          = "#7C4DFF",
  Mast_cell                    = "#AF2D8C",
  Basophil                     = "#E146A1",
  B_cell                       = "#E83E8C",
  Plasma_cell                  = "#2F8FD8",
  CD4_T_cell                   = "#5B5CE2",
  CD8_T_cell                   = "#6A5ACD",
  Treg                         = "#A0208F",
  Gamma_delta_T                = "#7B2CBF",
  NK_cell                      = "#C89D00",
  NKT_cell                     = "#9C7A00",
  ILC                          = "#3B8D5B",
  Erythroid                    = "#008C2E",
  Megakaryocyte                = "#AD2E97",
  Platelet                     = "#E95B9F",
  Cycling                      = "#8A2BE2",
  Unknown                      = "#777777"
)

ueno_subtype_palette_v407 <- c(
  Mature_hepatocyte                    = "#009688",
  Periportal_hepatocyte                = "#00A67F",
  Midzonal_hepatocyte                  = "#008D7A",
  Pericentral_hepatocyte               = "#006F63",
  Stress_response_hepatocyte           = "#8F9800",
  IFN_response_hepatocyte              = "#4B8A55",
  Cycling_hepatocyte                   = "#8A2BE2",

  Biliary                              = "#00B84A",
  Cholangiocyte                        = "#00B84A",
  Reactive_cholangiocyte               = "#00A33D",

  LSEC                                 = "#00AEEF",
  Vascular_endothelial                 = "#B86B00",
  Periportal_LSEC                      = "#00C4E0",
  Pericentral_LSEC                     = "#008FCA",
  Capillarized_LSEC                    = "#007CA8",
  Angiogenic_LSEC                      = "#1E88E5",
  Inflammatory_LSEC                    = "#3F73B7",

  Mesenchymal                          = "#F34BA6",
  qHSC                                 = "#F34BA6",
  aHSC                                 = "#D62C8A",
  Early_activated_HSC                  = "#EA3F98",
  Myofibroblastic_aHSC                 = "#C2187A",
  Fibrogenic_aHSC                      = "#A40E68",
  Inflammatory_aHSC                    = "#D83A8E",
  Portal_fibroblast                    = "#C4589C",
  Pericyte_VSMC                        = "#8B5A2B",

  Kupffer_macrophage                   = "#2F8F2F",
  Monocyte_derived_macrophage          = "#00A67A",
  Resident_Kupffer_like                = "#2F8F2F",
  Monocyte_like                        = "#00A67A",
  Inflammatory_M1_like                 = "#E53935",
  Pro_resolution_M2_like               = "#2E8B57",
  SPP1_TREM2_MASH_associated           = "#B02A9B",
  Lipid_associated_macrophage          = "#7D3CB5",
  Efferocytosis_phagocytosis_high      = "#26734D",
  IL10_response_high_Mphi              = "#17683F",
  Fibrosis_associated_Mphi             = "#8D2F56",
  IFN_response_macrophage              = "#3949AB",

  Classical_monocyte                   = "#00C08B",
  Nonclassical_monocyte                = "#008B72",
  Neutrophil                           = "#D44A9A",
  Inflammatory_neutrophil              = "#C72E8B",
  Aged_neutrophil                      = "#9A246D",

  cDC1                                 = "#C57A00",
  cDC2                                 = "#A45A00",
  pDC                                  = "#7C4DFF",

  Lymphoid                             = "#3D5AFE",
  B_cell                               = "#E83E8C",
  Naive_B                              = "#F45A9C",
  Memory_B                             = "#D42579",
  Plasma_cell                          = "#2F8FD8",

  CD4_T_cell                           = "#5B5CE2",
  Naive_CD4_T                          = "#7475F2",
  Activated_CD4_T                      = "#4A49C8",
  CD8_T_cell                           = "#6A5ACD",
  Naive_CD8_T                          = "#8B70D8",
  Cytotoxic_CD8_T                      = "#5B3DB0",
  Exhausted_CD8_T                      = "#452A88",
  NK_cell                              = "#C89D00",
  Activated_NK                        = "#B78B00",
  NKT_cell                             = "#9C7A00",
  Treg                                 = "#A0208F",
  Gamma_delta_T                        = "#7B2CBF",

  Basophil                             = "#E146A1",
  Mast_cell                            = "#AF2D8C",
  Erythroid                            = "#008C2E",
  Cycling                              = "#8A2BE2",
  Unknown                              = "#777777"
)

# ------------------------------------------------------------
# Cell-type to broader plotting group
# ------------------------------------------------------------

infer_celltype_group_v407 <- function(celltypes) {
  x <- as.character(celltypes)

  out <- rep("Other", length(x))

  out[grepl(
    "hepatocyte|Hepatic_progenitor|Biliary|Cholangiocyte",
    x, ignore.case = TRUE
  )] <- "Parenchymal_biliary"

  out[grepl(
    "LSEC|endothelial",
    x, ignore.case = TRUE
  )] <- "Endothelial"

  out[grepl(
    "HSC|Mesenchymal|fibroblast|Pericyte|VSMC",
    x, ignore.case = TRUE
  )] <- "Mesenchymal"

  out[grepl(
    "Kupffer|macrophage|Monocyte|Neutrophil|cDC|pDC|Mast|Basophil",
    x, ignore.case = TRUE
  )] <- "Myeloid"

  out[grepl(
    "B_cell|Plasma|CD4|CD8|Treg|Gamma|NK|NKT|ILC|Lymphoid",
    x, ignore.case = TRUE
  )] <- "Lymphoid"

  out[grepl(
    "Erythroid|Megakaryocyte|Platelet",
    x, ignore.case = TRUE
  )] <- "Blood_other"

  out
}

# ------------------------------------------------------------
# One most-characteristic positive marker per cell type
# ------------------------------------------------------------

find_one_characteristic_marker_v407 <- function(
    object,
    identity_col = "vote_ueno_summary_v40",
    assay = NULL,
    slot = "data",
    min_pct = 0.20,
    logfc_threshold = 0.25,
    exclude_regex = "^(MT-|RPL|RPS|HBA|HBB)|^MALAT1$|^XIST$|^JUN$|^FOS$|^FOSB$|^IER2$|^IER3$"
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
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )

  if (nrow(markers) == 0L) {
    stop("No positive markers were detected.")
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
    !grepl(exclude_regex, markers$gene, ignore.case = FALSE),
    ,
    drop = FALSE
  ]

  markers$specificity_delta <- pmax(markers$pct.1 - markers$pct.2, 0)
  markers$significance_weight <- pmin(
    -log10(pmax(markers$p_val_adj, 1e-300)),
    50
  ) / 50

  markers$characteristic_score <- (
    pmax(markers[[fc_col]], 0) *
      pmax(markers$specificity_delta, 0.01) *
      pmax(markers$pct.1, 0.01) *
      (0.5 + 0.5 * markers$significance_weight)
  )

  markers <- markers[
    order(
      markers$cluster,
      -markers$characteristic_score,
      -markers[[fc_col]],
      -markers$specificity_delta,
      markers$p_val_adj
    ),
    ,
    drop = FALSE
  ]

  selected <- do.call(
    rbind,
    lapply(
      split(markers, markers$cluster),
      function(df) df[1, , drop = FALSE]
    )
  )

  rownames(selected) <- NULL
  selected$celltype_group <- infer_celltype_group_v407(selected$cluster)
  attr(selected, "fc_col") <- fc_col
  selected
}

# ------------------------------------------------------------
# Build one-row violin plot for one marker gene
# ------------------------------------------------------------

make_marker_violin_row_v407 <- function(
    object,
    gene,
    identity_col,
    condition_col,
    condition_value,
    celltypes,
    palette,
    assay = NULL,
    slot = "data",
    show_x_text = FALSE,
    y_max = NULL
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  md <- object[[]]
  keep_cells <- rownames(md)[
    as.character(md[[condition_col]]) == condition_value &
      as.character(md[[identity_col]]) %in% celltypes
  ]

  if (length(keep_cells) == 0L) {
    return(NULL)
  }

  sub_obj <- Seurat::subset(object, cells = keep_cells)
  sub_obj[[identity_col]] <- factor(
    as.character(sub_obj[[]][[identity_col]]),
    levels = celltypes
  )

  plot_palette <- resolve_palette_v402(
    values = celltypes,
    palette = palette,
    fallback = "#777777"
  )

  p <- Seurat::VlnPlot(
    object = sub_obj,
    features = gene,
    group.by = identity_col,
    assay = assay,
    slot = slot,
    pt.size = 0,
    cols = unname(plot_palette),
    combine = TRUE
  ) +
    ggplot2::labs(
      title = gene,
      x = NULL,
      y = "Expression"
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),
      axis.title.y = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none"
    )

  if (isTRUE(show_x_text)) {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 55,
          hjust = 1,
          vjust = 1,
          size = 7
        )
      )
  } else {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank()
      )
  }

  if (!is.null(y_max) && is.finite(y_max)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
  }

  p
}

# ------------------------------------------------------------
# Multipage PDF:
# each page = one lineage group
# columns = conditions
# rows = one characteristic marker per cell type
# ------------------------------------------------------------

save_grouped_one_marker_violin_v407 <- function(
    object,
    marker_table,
    output_pdf,
    identity_col = "vote_ueno_summary_v40",
    condition_col,
    palette = ueno_subtype_palette_v407,
    assay = NULL,
    slot = "data",
    source_rds = NULL,
    created_at = Sys.time(),
    condition_order = NULL
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  if (is.null(condition_order)) {
    condition_order <- unique(
      as.character(object[[]][[condition_col]])
    )
    condition_order <- condition_order[
      !is.na(condition_order) & nzchar(condition_order)
    ]
  }

  group_order <- c(
    "Parenchymal_biliary",
    "Endothelial",
    "Mesenchymal",
    "Myeloid",
    "Lymphoid",
    "Blood_other",
    "Other"
  )

  groups <- intersect(group_order, unique(marker_table$celltype_group))

  grDevices::pdf(
    output_pdf,
    width = max(12, 6.0 * length(condition_order)),
    height = 11,
    onefile = TRUE,
    useDingbats = FALSE
  )

  on.exit(grDevices::dev.off(), add = TRUE)

  for (grp in groups) {
    tab <- marker_table[
      marker_table$celltype_group == grp,
      ,
      drop = FALSE
    ]

    if (nrow(tab) == 0L) next

    celltypes <- as.character(tab$cluster)
    genes <- as.character(tab$gene)

    condition_columns <- list()

    for (cond in condition_order) {
      rows <- list()

      for (i in seq_along(genes)) {
        rows[[i]] <- make_marker_violin_row_v407(
          object = object,
          gene = genes[[i]],
          identity_col = identity_col,
          condition_col = condition_col,
          condition_value = cond,
          celltypes = celltypes,
          palette = palette,
          assay = assay,
          slot = slot,
          show_x_text = i == length(genes)
        )
      }

      rows <- Filter(Negate(is.null), rows)
      if (length(rows) == 0L) next

      condition_columns[[cond]] <- patchwork::wrap_plots(
        rows,
        ncol = 1
      ) +
        patchwork::plot_annotation(
          title = cond,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "bold",
              size = 16
            )
          )
        )
    }

    condition_columns <- Filter(Negate(is.null), condition_columns)
    if (length(condition_columns) == 0L) next

    page <- patchwork::wrap_plots(
      condition_columns,
      nrow = 1
    ) +
      patchwork::plot_annotation(
        title = paste0(
          "One characteristic marker per cell type: ",
          gsub("_", " ", grp, fixed = TRUE)
        ),
        subtitle = paste0(
          "Rows: marker genes | x-axis: cell types | columns: ",
          condition_col
        ),
        caption = make_caption_v402(
          source_rds = source_rds,
          created_at = created_at
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold",
            size = 16
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5,
            size = 9
          ),
          plot.caption = ggplot2::element_text(
            hjust = 0,
            size = 7,
            colour = "grey25"
          )
        )
      )

    print(page)
  }

  invisible(output_pdf)
}
