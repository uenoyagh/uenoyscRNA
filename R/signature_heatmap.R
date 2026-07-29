# ============================================================
# Functional signature heatmap
# ============================================================

match_features_ci <- function(features, genes) {
  idx <- match(toupper(genes), toupper(features))
  features[idx[!is.na(idx)]]
}

get_assay_data_safe <- function(object, assay) {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) {
      SeuratObject::GetAssayData(object, assay = assay, slot = "data")
    }
  )
}

calculate_signature_cluster_scores <- function(
    object, cluster_column, signatures,
    assay = NULL, min_genes = 2) {

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  expr <- get_assay_data_safe(object, assay)
  clusters <- as.character(object[[cluster_column, drop = TRUE]])
  cluster_levels <- natural_level_order(clusters)

  score_rows <- list()
  report_rows <- list()

  for (signature_name in names(signatures)) {
    matched <- match_features_ci(
      rownames(expr),
      signatures[[signature_name]]
    )

    report_rows[[signature_name]] <- data.frame(
      signature = signature_name,
      requested_n = length(signatures[[signature_name]]),
      matched_n = length(matched),
      genes_used = paste(matched, collapse = " | "),
      stringsAsFactors = FALSE
    )

    if (length(matched) < min_genes) next

    cell_score <- Matrix::colMeans(expr[matched, , drop = FALSE])

    tmp <- data.frame(
      cluster = clusters,
      score = as.numeric(cell_score),
      stringsAsFactors = FALSE
    )

    agg <- stats::aggregate(score ~ cluster, tmp, mean)
    agg$signature <- signature_name
    score_rows[[signature_name]] <- agg
  }

  if (length(score_rows) == 0) {
    stop("No signature had enough matched genes.")
  }

  scores <- do.call(rbind, score_rows)
  scores$cluster <- factor(scores$cluster, levels = cluster_levels)

  scores$z_score <- ave(
    scores$score,
    scores$signature,
    FUN = function(x) {
      z <- as.numeric(scale(x))
      z[is.na(z)] <- 0
      z
    }
  )

  list(
    scores = scores,
    gene_report = do.call(rbind, report_rows)
  )
}

make_signature_heatmap <- function(
    score_df, clip = 2,
    title = "Functional signature scores by existing RPCA cluster",
    subtitle = "Cluster-level z-scores; red = high score") {

  score_df$z_plot <- pmax(pmin(score_df$z_score, clip), -clip)

  ggplot2::ggplot(
    score_df,
    ggplot2::aes(x = cluster, y = signature, fill = z_plot)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(-clip, clip),
      name = "Cluster\nZ-score"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 17, face = "bold", hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 12, hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 55, hjust = 1, vjust = 1,
        size = 10, face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 11, face = "bold"
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )
}
