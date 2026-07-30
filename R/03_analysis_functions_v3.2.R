prepare_seurat_for_markers <- function(object, cfg, log_file = NULL) {
  Seurat::DefaultAssay(object) <- cfg$assay
  assay_obj <- object[[cfg$assay]]

  if (inherits(assay_obj, "Assay5")) {
    layers <- SeuratObject::Layers(assay_obj)
    split_like <- sum(grepl("^data\\.", layers)) > 1L || sum(grepl("^counts\\.", layers)) > 1L
    if (split_like) {
      log_msg("Joining Seurat v5 assay layers before marker analysis.", log_file = log_file)
      object[[cfg$assay]] <- SeuratObject::JoinLayers(assay_obj)
    }
  }
  object
}

get_layer_matrix <- function(object, assay, layer) {
  mat <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    stop("Could not retrieve assay='", assay, "', layer='", layer, "'.")
  }
  mat
}

run_presto_markers <- function(object, cluster_col, cfg, log_file = NULL) {
  assert_packages(c("presto", "dplyr", "tibble"))
  mat <- get_layer_matrix(object, cfg$assay, cfg$data_layer)
  groups <- object[[]][[cluster_col]]
  names(groups) <- rownames(object[[]])

  keep_cells <- intersect(colnames(mat), names(groups))
  mat <- mat[, keep_cells, drop = FALSE]
  groups <- droplevels(factor(groups[keep_cells]))

  tab <- table(groups)
  valid <- names(tab)[tab >= cfg$min_cells_per_cluster]
  keep <- groups %in% valid
  mat <- mat[, keep, drop = FALSE]
  groups <- droplevels(groups[keep])

  log_msg(
    "Running presto::wilcoxauc on ", ncol(mat), " cells and ",
    nlevels(groups), " clusters.", log_file = log_file
  )

  res <- presto::wilcoxauc(mat, groups)
  res <- as.data.frame(res)
  names(res)[names(res) == "group"] <- "cluster"
  names(res)[names(res) == "feature"] <- "gene"

  if ("logFC" %in% names(res)) {
    res <- res[abs(res$logFC) >= cfg$logfc_threshold, , drop = FALSE]
  }
  if ("pct_in" %in% names(res)) {
    res <- res[res$pct_in >= cfg$min_pct, , drop = FALSE]
  }
  if (isTRUE(cfg$only_pos) && "logFC" %in% names(res)) {
    res <- res[res$logFC > 0, , drop = FALSE]
  }

  rank_col <- if ("auc" %in% names(res)) "auc" else if ("logFC" %in% names(res)) "logFC" else names(res)[1]
  res <- res[order(res$cluster, -res[[rank_col]]), , drop = FALSE]
  rownames(res) <- NULL
  res
}

run_seurat_markers <- function(object, cluster_col, cfg, log_file = NULL) {
  Seurat::Idents(object) <- object[[]][[cluster_col]]
  log_msg("Running Seurat::FindAllMarkers.", log_file = log_file)
  Seurat::FindAllMarkers(
    object = object,
    assay = cfg$assay,
    only.pos = cfg$only_pos,
    min.pct = cfg$min_pct,
    logfc.threshold = cfg$logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )
}

summarize_top_markers <- function(markers, n = 100L) {
  rank_col <- intersect(c("auc", "avg_log2FC", "avg_logFC", "logFC"), names(markers))[1]
  if (is.na(rank_col)) stop("No marker ranking column found.")
  split_markers <- split(markers, markers$cluster)
  do.call(rbind, lapply(split_markers, function(x) {
    x <- x[order(-x[[rank_col]]), , drop = FALSE]
    utils::head(x, n)
  }))
}

cluster_mean_expression <- function(object, cluster_col, genes, cfg) {
  genes <- intersect(unique(genes), rownames(object))
  if (!length(genes)) return(matrix(numeric(), nrow = 0, ncol = 0))

  avg <- Seurat::AverageExpression(
    object,
    assays = cfg$assay,
    features = genes,
    group.by = cluster_col,
    layer = cfg$data_layer,
    verbose = FALSE
  )[[cfg$assay]]
  as.matrix(avg)
}

score_marker_set <- function(expr, positive, negative, cfg) {
  positive <- intersect(positive, rownames(expr))
  negative <- intersect(negative, rownames(expr))

  if (length(positive) < cfg$min_present_genes) {
    return(rep(NA_real_, ncol(expr)))
  }

  pos <- colMeans(expr[positive, , drop = FALSE])
  neg <- if (length(negative)) colMeans(expr[negative, , drop = FALSE]) else rep(0, ncol(expr))
  cfg$positive_weight * pos - cfg$negative_weight * neg
}

run_marker_voting <- function(object, cluster_col, marker_reference, cfg) {
  all_genes <- unique(marker_reference$gene)
  expr <- cluster_mean_expression(object, cluster_col, all_genes, cfg)
  if (!nrow(expr)) stop("None of the marker reference genes were found in the object.")

  combos <- unique(marker_reference[c("source", "label")])
  scores <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
    src <- combos$source[i]
    lab <- combos$label[i]
    pos <- marker_reference$gene[
      marker_reference$source == src &
      marker_reference$label == lab &
      marker_reference$direction == "positive"
    ]
    neg <- marker_reference$gene[
      marker_reference$source == src &
      marker_reference$label == lab &
      marker_reference$direction == "negative"
    ]
    sc <- score_marker_set(expr, pos, neg, cfg)
    data.frame(
      source = src,
      label = lab,
      cluster = colnames(expr),
      score = as.numeric(sc),
      n_positive_present = length(intersect(pos, rownames(expr))),
      n_negative_present = length(intersect(neg, rownames(expr))),
      stringsAsFactors = FALSE
    )
  }))

  valid <- scores[is.finite(scores$score), , drop = FALSE]
  winners <- do.call(rbind, lapply(split(valid, list(valid$source, valid$cluster), drop = TRUE), function(x) {
    x <- x[order(-x$score), , drop = FALSE]
    top <- x[1, , drop = FALSE]
    second <- if (nrow(x) >= 2L) x$score[2] else NA_real_
    top$second_score <- second
    top$score_delta <- top$score - second
    top$confidence <- ifelse(
      top$score < cfg$low_confidence_score, "low",
      ifelse(!is.na(top$score_delta) && top$score_delta < cfg$ambiguity_delta, "ambiguous", "high")
    )
    top
  }))
  rownames(winners) <- NULL

  list(scores = scores, winners = winners, expression = expr)
}

attach_voting_metadata <- function(object, cluster_col, winners) {
  meta <- object[[]]
  cl <- as.character(meta[[cluster_col]])

  if (!all(startsWith(cl, "g"))) {
    cl <- paste0("g", cl)
  }
  for (src in unique(winners$source)) {
    sub <- winners[winners$source == src, , drop = FALSE]
    label_map <- setNames(sub$label, sub$cluster)
    score_map <- setNames(sub$score, sub$cluster)
    conf_map <- setNames(sub$confidence, sub$cluster)
    prefix <- paste0("vote_", tolower(src))
    meta[[paste0(prefix, "_label")]] <- unname(label_map[cl])
    meta[[paste0(prefix, "_score")]] <- unname(score_map[cl])
    meta[[paste0(prefix, "_confidence")]] <- unname(conf_map[cl])
  }
  object@meta.data <- meta
  object
}
