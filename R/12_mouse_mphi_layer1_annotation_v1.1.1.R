# ============================================================
# uenoyscRNA Mouse macrophage RPCA layer-1 annotation v1.1.1
# Six-class macrophage framework
# Palette: UenoMphiLayer1_v1
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ------------------------------------------------------------
# Reproducible named palette
# ------------------------------------------------------------

UenoMphiLayer1_v1 <- c(
  "Resident Kupffer-like Mphi" = "#168A4A",
  "Monocyte-like Mphi" = "#00A67A",
  "Inflammatory M1-like Mphi" = "#E53935",
  "Pro-resolution M2-like Mphi" = "#2E8B57",
  "SPP1/TREM2 MASH-associated Mphi" = "#B02A9B",
  "Other" = "#7A7A7A"
)

# ------------------------------------------------------------
# Marker sets
# ------------------------------------------------------------

UenoMphiLayer1_markers_v1 <- list(
  "Resident Kupffer-like Mphi" = c(
    "Clec4f", "Timd4", "Vsig4", "Marco", "Cd5l",
    "C1qa", "C1qb", "C1qc", "Fcgr2b", "Slc40a1"
  ),
  "Monocyte-like Mphi" = c(
    "Lyz2", "Ccr2", "Ly6c2", "S100a8", "S100a9",
    "Plac8", "Ctss", "Fcgr3", "Ms4a7", "Tyrobp"
  ),
  "Inflammatory M1-like Mphi" = c(
    "Il1b", "Tnf", "Ccl2", "Ccl3", "Ccl4",
    "Cxcl2", "Cxcl9", "Cxcl10", "Nos2", "Ptgs2",
    "Cd80", "Cd86", "Nfkbia", "Irf1", "Stat1"
  ),
  "Pro-resolution M2-like Mphi" = c(
    "Mrc1", "Cd163", "Il10", "Maf", "Mafb",
    "Klf4", "Arg1", "Retnla", "Chil3", "Ccl24",
    "Lgals3", "Gas6", "Axl", "Mertk"
  ),
  "SPP1/TREM2 MASH-associated Mphi" = c(
    "Spp1", "Trem2", "Gpnmb", "Cd9", "Lgals3",
    "Lpl", "Fabp5", "Ctsb", "Ctsd", "Ctsk",
    "Apoe", "Cst7", "Itgax"
  )
)

UenoMphiLayer1_exclusion_v1 <- list(
  "Resident Kupffer-like Mphi" = c("S100a8", "S100a9", "Ly6c2", "Ccr2"),
  "Monocyte-like Mphi" = c("Clec4f", "Timd4", "Vsig4", "Marco"),
  "Inflammatory M1-like Mphi" = c("Mrc1", "Retnla", "Chil3"),
  "Pro-resolution M2-like Mphi" = c("Il1b", "Tnf", "Cxcl10", "Nos2"),
  "SPP1/TREM2 MASH-associated Mphi" = c("Clec4f", "Timd4")
)

resolve_features_case_insensitive_mphi_v110 <- function(
    features,
    available_features
) {
  available_lower <- tolower(available_features)
  idx <- match(tolower(features), available_lower)
  unique(available_features[idx[!is.na(idx)]])
}

detect_metadata_column_mphi_v110 <- function(
    object,
    candidates,
    required = TRUE,
    label = "metadata"
) {
  cols <- colnames(object[[]])
  hit <- candidates[candidates %in% cols]

  if (length(hit) > 0L) {
    return(hit[[1]])
  }

  if (required) {
    stop(
      label, " column was not found.\nCandidates:\n",
      paste(candidates, collapse = ", "),
      "\n\nAvailable metadata:\n",
      paste(cols, collapse = ", ")
    )
  }

  NULL
}

detect_reduction_mphi_v110 <- function(object) {
  reductions <- names(object@reductions)
  candidates <- c(
    "umapRPCA", "umap_rpca", "rpca.umap", "umap", "UMAP"
  )

  hit <- candidates[candidates %in% reductions]

  if (length(hit) > 0L) {
    return(hit[[1]])
  }

  hit <- reductions[grepl("umap", reductions, ignore.case = TRUE)]

  if (length(hit) > 0L) {
    return(hit[[1]])
  }

  stop(
    "UMAP reduction was not found. Available reductions:\n",
    paste(reductions, collapse = ", ")
  )
}


normalize_aggregate_cluster_names_mphi_v111 <- function(
    aggregate_names,
    original_cluster_values
) {
  original_cluster_values <- unique(as.character(original_cluster_values))
  original_cluster_values <- original_cluster_values[!is.na(original_cluster_values)]

  normalized <- as.character(aggregate_names)

  # Seurat v5 prefixes numeric group names with "g", e.g. 0 -> g0.
  prefixed_numeric <- grepl("^g[-+]?[0-9]+([.][0-9]+)?$", normalized)
  stripped <- sub("^g", "", normalized[prefixed_numeric])

  can_strip <- stripped %in% original_cluster_values
  normalized[prefixed_numeric][can_strip] <- stripped[can_strip]

  # General fallback: preserve names already matching original cluster IDs.
  normalized
}

compute_cluster_marker_scores_mphi_v110 <- function(
    object,
    cluster_col,
    assay = NULL,
    slot = "data",
    marker_sets = UenoMphiLayer1_markers_v1,
    exclusion_sets = UenoMphiLayer1_exclusion_v1
) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  clusters <- sort(unique(as.character(object[[]][[cluster_col]])))
  clusters <- clusters[!is.na(clusters)]

  available_features <- rownames(object)

  resolved_positive <- lapply(
    marker_sets,
    resolve_features_case_insensitive_mphi_v110,
    available_features = available_features
  )

  resolved_exclusion <- lapply(
    exclusion_sets,
    resolve_features_case_insensitive_mphi_v110,
    available_features = available_features
  )

  all_features <- unique(c(
    unlist(resolved_positive, use.names = FALSE),
    unlist(resolved_exclusion, use.names = FALSE)
  ))

  if (length(all_features) == 0L) {
    stop("None of the configured macrophage markers were found.")
  }

  avg <- Seurat::AverageExpression(
    object = object,
    assays = assay,
    features = all_features,
    group.by = cluster_col,
    slot = slot,
    verbose = FALSE
  )[[assay]]

  avg <- as.matrix(avg)

  original_cluster_values <- as.character(object[[]][[cluster_col]])
  colnames(avg) <- normalize_aggregate_cluster_names_mphi_v111(
    aggregate_names = colnames(avg),
    original_cluster_values = original_cluster_values
  )

  # z-score each gene across clusters
  z <- t(scale(t(avg)))
  z[!is.finite(z)] <- 0

  rows <- list()

  for (cluster_id in colnames(z)) {
    for (class_name in names(marker_sets)) {
      pos_genes <- resolved_positive[[class_name]]
      neg_genes <- resolved_exclusion[[class_name]]

      pos_score <- if (length(pos_genes) > 0L) {
        mean(z[pos_genes, cluster_id, drop = TRUE], na.rm = TRUE)
      } else {
        NA_real_
      }

      neg_score <- if (length(neg_genes) > 0L) {
        mean(z[neg_genes, cluster_id, drop = TRUE], na.rm = TRUE)
      } else {
        0
      }

      final_score <- pos_score - 0.5 * neg_score

      rows[[length(rows) + 1L]] <- data.frame(
        cluster = cluster_id,
        class = class_name,
        positive_score = pos_score,
        exclusion_score = neg_score,
        final_score = final_score,
        marker_coverage = length(pos_genes) / length(marker_sets[[class_name]]),
        positive_markers_found = paste(pos_genes, collapse = ";"),
        exclusion_markers_found = paste(neg_genes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }

  score_table <- do.call(rbind, rows)

  annotation_rows <- lapply(
    split(score_table, score_table$cluster),
    function(df) {
      df <- df[order(-df$final_score), , drop = FALSE]

      top1 <- df[1, , drop = FALSE]
      top2_score <- if (nrow(df) >= 2L) df$final_score[[2]] else NA_real_
      delta <- top1$final_score[[1]] - top2_score

      recommended <- if (
        is.na(top1$final_score[[1]]) ||
        top1$final_score[[1]] < 0.25 ||
        is.na(delta) ||
        delta < 0.15 ||
        top1$marker_coverage[[1]] < 0.20
      ) {
        "Other"
      } else {
        top1$class[[1]]
      }

      confidence <- if (recommended == "Other") {
        "Low"
      } else if (
        top1$final_score[[1]] >= 1.0 &&
        delta >= 0.50 &&
        top1$marker_coverage[[1]] >= 0.40
      ) {
        "High"
      } else {
        "Moderate"
      }

      data.frame(
        cluster = top1$cluster[[1]],
        mphi_layer1 = recommended,
        top_score = top1$final_score[[1]],
        score_delta = delta,
        confidence = confidence,
        marker_coverage = top1$marker_coverage[[1]],
        stringsAsFactors = FALSE
      )
    }
  )

  annotation_table <- do.call(rbind, annotation_rows)
  rownames(annotation_table) <- NULL

  list(
    score_table = score_table,
    annotation_table = annotation_table
  )
}

apply_cluster_annotation_mphi_v110 <- function(
    object,
    cluster_col,
    annotation_table,
    output_col = "mphi_layer1_v110",
    confidence_col = "mphi_layer1_confidence_v110"
) {
  metadata <- object[[]]
  cluster_values <- as.character(metadata[[cluster_col]])

  annotation_map <- setNames(
    annotation_table$mphi_layer1,
    annotation_table$cluster
  )

  confidence_map <- setNames(
    annotation_table$confidence,
    annotation_table$cluster
  )

  mapped_annotation <- unname(annotation_map[cluster_values])
  mapped_confidence <- unname(confidence_map[cluster_values])

  unmatched_clusters <- sort(unique(cluster_values[is.na(mapped_annotation)]))
  unmatched_clusters <- unmatched_clusters[!is.na(unmatched_clusters)]

  if (length(unmatched_clusters) > 0L) {
    stop(
      "Cluster annotation mapping failed for cluster IDs: ",
      paste(unmatched_clusters, collapse = ", "),
      "\nAvailable annotation IDs: ",
      paste(names(annotation_map), collapse = ", ")
    )
  }

  object[[output_col]] <- factor(
    mapped_annotation,
    levels = c(
      "Resident Kupffer-like Mphi",
      "Monocyte-like Mphi",
      "Inflammatory M1-like Mphi",
      "Pro-resolution M2-like Mphi",
      "SPP1/TREM2 MASH-associated Mphi",
      "Other"
    )
  )

  object[[confidence_col]] <- factor(
    mapped_confidence,
    levels = c("High", "Moderate", "Low")
  )

  object
}

publish_mphi_layer1_umap_v110 <- function(
    object,
    group_by,
    reduction,
    title,
    pt_size = 1.05,
    label_size = 4.6,
    source_rds = NULL,
    created_at = Sys.time()
) {
  p <- Seurat::DimPlot(
    object = object,
    reduction = reduction,
    group.by = group_by,
    pt.size = pt_size,
    label = TRUE,
    repel = TRUE,
    label.size = label_size,
    raster = TRUE,
    raster.dpi = c(700, 700)
  )

  if (length(p$layers) > 0L) {
    for (i in seq_along(p$layers)) {
      p$layers[[i]]$aes_params$alpha <- 1
    }
  }

  p +
    ggplot2::scale_color_manual(
      values = UenoMphiLayer1_v1,
      drop = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title,
      color = NULL,
      caption = paste0(
        "Source RDS: ", basename(source_rds),
        " | Palette: UenoMphiLayer1_v1",
        " | Created: ",
        format(created_at, "%Y-%m-%d %H:%M:%S %Z")
      )
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      legend.text = ggplot2::element_text(size = 9.5),
      legend.key.height = grid::unit(0.48, "cm"),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7.5,
        colour = "grey25"
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4.5, alpha = 1)
      )
    )
}
