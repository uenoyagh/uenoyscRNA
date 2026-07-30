# Hierarchical marker-voting functions for v3.3.
# Reuses prepare_seurat_for_markers(), run_presto_markers(),
# run_seurat_markers(), summarize_top_markers(), cluster_mean_expression()
# from R/03_analysis_functions_v3.2.R.

score_marker_set_v33 <- function(expr, positive, negative, cfg) {
  positive_all <- unique(positive)
  negative_all <- unique(negative)
  positive <- intersect(positive_all, rownames(expr))
  negative <- intersect(negative_all, rownames(expr))

  n_pos_total <- length(positive_all)
  n_neg_total <- length(negative_all)
  n_pos_present <- length(positive)
  n_neg_present <- length(negative)

  if (n_pos_present < cfg$min_present_genes) {
    return(list(
      score = rep(NA_real_, ncol(expr)),
      positive_mean = rep(NA_real_, ncol(expr)),
      negative_mean = rep(NA_real_, ncol(expr)),
      n_positive_total = n_pos_total,
      n_positive_present = n_pos_present,
      positive_coverage = if (n_pos_total) n_pos_present / n_pos_total else NA_real_,
      n_negative_total = n_neg_total,
      n_negative_present = n_neg_present,
      negative_coverage = if (n_neg_total) n_neg_present / n_neg_total else NA_real_
    ))
  }

  pos_mean <- colMeans(expr[positive, , drop = FALSE])
  neg_mean <- if (n_neg_present) {
    colMeans(expr[negative, , drop = FALSE])
  } else {
    rep(0, ncol(expr))
  }

  list(
    score = cfg$positive_weight * pos_mean - cfg$negative_weight * neg_mean,
    positive_mean = pos_mean,
    negative_mean = neg_mean,
    n_positive_total = n_pos_total,
    n_positive_present = n_pos_present,
    positive_coverage = n_pos_present / n_pos_total,
    n_negative_total = n_neg_total,
    n_negative_present = n_neg_present,
    negative_coverage = if (n_neg_total) n_neg_present / n_neg_total else NA_real_
  )
}

score_reference_level_v33 <- function(expr, marker_reference, level, cfg) {
  ref <- marker_reference[marker_reference$level == level, , drop = FALSE]
  combos <- unique(ref[c(
    "source", "level", "parent_lineage", "parent_celltype", "label"
  )])

  out <- lapply(seq_len(nrow(combos)), function(i) {
    z <- combos[i, , drop = FALSE]
    keep <- ref$source == z$source &
      ref$level == z$level &
      ref$label == z$label

    pos <- ref$gene[keep & ref$direction == "positive"]
    neg <- ref$gene[keep & ref$direction == "negative"]
    sc <- score_marker_set_v33(expr, pos, neg, cfg)

    data.frame(
      source = z$source,
      level = z$level,
      parent_lineage = z$parent_lineage,
      parent_celltype = z$parent_celltype,
      label = z$label,
      cluster = colnames(expr),
      score = as.numeric(sc$score),
      positive_mean = as.numeric(sc$positive_mean),
      negative_mean = as.numeric(sc$negative_mean),
      n_positive_total = sc$n_positive_total,
      n_positive_present = sc$n_positive_present,
      positive_coverage = sc$positive_coverage,
      n_negative_total = sc$n_negative_total,
      n_negative_present = sc$n_negative_present,
      negative_coverage = sc$negative_coverage,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

classify_confidence_v33 <- function(score, delta, coverage, cfg) {
  min_cov <- cfg$ueno_min_positive_coverage %||% 0.25
  moderate_delta <- cfg$ueno_moderate_delta %||% cfg$ambiguity_delta
  high_delta <- cfg$ueno_high_delta %||% max(cfg$ambiguity_delta * 2, moderate_delta)

  if (!is.finite(score) || !is.finite(coverage) || coverage < min_cov) {
    return("Unknown")
  }
  if (score < cfg$low_confidence_score) {
    return("Low")
  }
  if (!is.finite(delta) || delta < moderate_delta) {
    return("Low")
  }
  if (delta < high_delta) {
    return("Moderate")
  }
  "High"
}

choose_winners_v33 <- function(scores, cfg, eligible = NULL) {
  x <- scores[is.finite(scores$score), , drop = FALSE]

  if (!is.null(eligible)) {
    x <- eligible(x)
  }
  if (!nrow(x)) return(x)

  keys <- unique(x[c("cluster")])
  winners <- lapply(keys$cluster, function(cl) {
    y <- x[x$cluster == cl, , drop = FALSE]
    y <- y[order(-y$score, -y$positive_coverage, y$label), , drop = FALSE]
    top <- y[1, , drop = FALSE]
    second <- if (nrow(y) >= 2L) y$score[2] else NA_real_
    top$second_score <- second
    top$score_delta <- top$score - second
    top$confidence <- classify_confidence_v33(
      top$score, top$score_delta, top$positive_coverage, cfg
    )
    top
  })

  out <- do.call(rbind, winners)
  rownames(out) <- NULL
  out
}

run_hierarchical_voting_v33 <- function(object, cluster_col, marker_reference, cfg) {
  genes <- unique(marker_reference$gene)
  expr <- cluster_mean_expression(object, cluster_col, genes, cfg)
  if (!nrow(expr)) stop("None of the v3.3 reference genes were found in the object.")

  general_scores <- score_reference_level_v33(expr, marker_reference, "general", cfg)
  lineage_scores <- score_reference_level_v33(expr, marker_reference, "lineage", cfg)
  celltype_scores <- score_reference_level_v33(expr, marker_reference, "celltype", cfg)
  subtype_scores <- score_reference_level_v33(expr, marker_reference, "subtype", cfg)

  general_winners <- choose_winners_v33(general_scores, cfg)
  lineage_winners <- choose_winners_v33(lineage_scores, cfg)

  lineage_map <- setNames(lineage_winners$label, lineage_winners$cluster)
  celltype_winners <- choose_winners_v33(
    celltype_scores,
    cfg,
    eligible = function(x) {
      expected <- unname(lineage_map[x$cluster])
      x[!is.na(expected) & x$parent_lineage == expected, , drop = FALSE]
    }
  )

  celltype_map <- setNames(celltype_winners$label, celltype_winners$cluster)
  subtype_winners <- choose_winners_v33(
    subtype_scores,
    cfg,
    eligible = function(x) {
      expected <- unname(celltype_map[x$cluster])
      x[!is.na(expected) & x$parent_celltype == expected, , drop = FALSE]
    }
  )

  list(
    expression = expr,
    scores = rbind(general_scores, lineage_scores, celltype_scores, subtype_scores),
    winners = rbind(general_winners, lineage_winners, celltype_winners, subtype_winners),
    general = list(scores = general_scores, winners = general_winners),
    lineage = list(scores = lineage_scores, winners = lineage_winners),
    celltype = list(scores = celltype_scores, winners = celltype_winners),
    subtype = list(scores = subtype_scores, winners = subtype_winners)
  )
}

normalize_cluster_keys_v33 <- function(metadata_cluster, winner_clusters) {
  cl <- as.character(metadata_cluster)
  winner_clusters <- as.character(winner_clusters)

  if (all(cl %in% winner_clusters)) return(cl)

  with_g <- ifelse(startsWith(cl, "g"), cl, paste0("g", cl))
  if (all(with_g %in% winner_clusters)) return(with_g)

  without_g <- sub("^g", "", cl)
  if (all(without_g %in% winner_clusters)) return(without_g)

  stop(
    "Could not reconcile metadata cluster IDs with voting cluster IDs.\n",
    "Metadata examples: ", paste(utils::head(unique(cl), 6), collapse = ", "), "\n",
    "Voting examples: ", paste(utils::head(unique(winner_clusters), 6), collapse = ", ")
  )
}

attach_hierarchical_voting_metadata_v33 <- function(object, cluster_col, voting) {
  meta <- object[[]]
  winner_clusters <- unique(voting$winners$cluster)
  cl <- normalize_cluster_keys_v33(meta[[cluster_col]], winner_clusters)

  attach_level <- function(meta, winners, prefix, unknown_label = "Unknown") {
    if (!nrow(winners)) {
      meta[[prefix]] <- unknown_label
      meta[[paste0(prefix, "_score")]] <- NA_real_
      meta[[paste0(prefix, "_delta")]] <- NA_real_
      meta[[paste0(prefix, "_confidence")]] <- "Unknown"
      meta[[paste0(prefix, "_coverage")]] <- NA_real_
      return(meta)
    }

    label_map <- setNames(winners$label, winners$cluster)
    score_map <- setNames(winners$score, winners$cluster)
    delta_map <- setNames(winners$score_delta, winners$cluster)
    conf_map <- setNames(winners$confidence, winners$cluster)
    cov_map <- setNames(winners$positive_coverage, winners$cluster)

    labels <- unname(label_map[cl])
    conf <- unname(conf_map[cl])
    labels[is.na(labels) | conf == "Unknown"] <- unknown_label

    meta[[prefix]] <- labels
    meta[[paste0(prefix, "_score")]] <- unname(score_map[cl])
    meta[[paste0(prefix, "_delta")]] <- unname(delta_map[cl])
    meta[[paste0(prefix, "_confidence")]] <- conf
    meta[[paste0(prefix, "_coverage")]] <- unname(cov_map[cl])
    meta
  }

  meta <- attach_level(meta, voting$general$winners, "vote_general_label")
  meta <- attach_level(meta, voting$lineage$winners, "vote_ueno_lineage")
  meta <- attach_level(meta, voting$celltype$winners, "vote_ueno_celltype")
  meta <- attach_level(meta, voting$subtype$winners, "vote_ueno_subtype")

  # Backward compatibility.
  meta$vote_ueno_label <- meta$vote_ueno_subtype

  object@meta.data <- meta
  object
}

make_hierarchical_agreement_table_v33 <- function(object, cluster_col, annotation_col) {
  meta <- object[[]]
  fields <- c(
    annotation_col,
    "vote_general_label",
    "vote_ueno_lineage",
    "vote_ueno_celltype",
    "vote_ueno_subtype"
  )
  fields <- fields[fields %in% colnames(meta)]

  out <- lapply(fields, function(field) {
    tab <- as.data.frame(table(
      cluster = as.character(meta[[cluster_col]]),
      annotation_field = rep(field, nrow(meta)),
      annotation = as.character(meta[[field]])
    ))
    tab <- tab[tab$Freq > 0, , drop = FALSE]
    tab <- do.call(rbind, lapply(split(tab, tab$cluster), function(x) {
      x <- x[order(-x$Freq), , drop = FALSE]
      x$cluster_total <- sum(x$Freq)
      x$fraction <- x$Freq / x$cluster_total
      x
    }))
    tab
  })
  do.call(rbind, out)
}
