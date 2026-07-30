# Evidence-weighted hierarchical voting for v3.4.

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_marker_clusters_v34 <- function(markers) {
  markers$cluster <- as.character(markers$cluster)
  markers$cluster_key <- ifelse(
    startsWith(markers$cluster, "g"),
    markers$cluster,
    paste0("g", markers$cluster)
  )
  markers
}

rank_markers_v34 <- function(markers) {
  split_x <- split(markers, markers$cluster_key)
  out <- lapply(split_x, function(x) {
    x <- x[order(-x$auc, -x$logFC, -(x$pct_in - x$pct_out)), , drop = FALSE]
    x$marker_rank <- seq_len(nrow(x))
    x
  })
  do.call(rbind, out)
}

gene_evidence_v34 <- function(auc, logfc, pct_in, pct_out, rank) {
  auc_component <- pmax((auc - 0.5) / 0.5, 0)
  logfc_component <- pmax(tanh(logfc / 2), 0)
  specificity_component <- pmax((pct_in - pct_out) / 100, 0)
  rank_component <- 1 / log2(rank + 1)
  0.35 * auc_component +
    0.35 * logfc_component +
    0.20 * specificity_component +
    0.10 * rank_component
}

prepare_marker_evidence_v34 <- function(markers, cfg) {
  markers <- normalize_marker_clusters_v34(markers)
  markers <- rank_markers_v34(markers)

  markers$gene_evidence <- gene_evidence_v34(
    markers$auc,
    markers$logFC,
    markers$pct_in,
    markers$pct_out,
    markers$marker_rank
  )

  min_auc <- cfg$v34_min_auc %||% 0.55
  min_logfc <- cfg$v34_min_logfc %||% 0.15
  min_pct_diff <- cfg$v34_min_pct_diff %||% 5

  markers$evidence_pass <-
    markers$auc >= min_auc &
    markers$logFC >= min_logfc &
    (markers$pct_in - markers$pct_out) >= min_pct_diff

  markers
}

annotation_prior_maps_v34 <- function() {
  lineage <- c(
    Hepatocyte = "Hepatocyte",
    Cholangiocyte = "Biliary",
    LSEC = "Endothelial",
    Vascular_endothelial = "Endothelial",
    HSC_Mesenchymal = "Mesenchymal",
    Kupffer_Macrophage = "Myeloid",
    Monocyte = "Myeloid",
    Neutrophil = "Myeloid",
    Dendritic_cell = "Myeloid",
    T_cell = "Lymphoid",
    NK_cell = "Lymphoid",
    B_cell = "Lymphoid",
    Plasma_cell = "Lymphoid",
    RBC = "Erythroid",
    Platelet = "Megakaryocytic",
    Cycling = "Cycling",
    Mesothelial = "Mesenchymal"
  )

  celltype <- c(
    Hepatocyte = "Mature_hepatocyte",
    Cholangiocyte = "Cholangiocyte",
    LSEC = "LSEC",
    Vascular_endothelial = "Vascular_endothelial",
    HSC_Mesenchymal = "qHSC",
    Kupffer_Macrophage = "Kupffer_macrophage",
    Monocyte = "Monocyte",
    Neutrophil = "Neutrophil",
    Dendritic_cell = "cDC2",
    T_cell = "CD4_T_cell",
    NK_cell = "NK_cell",
    B_cell = "B_cell",
    Plasma_cell = "Plasma_cell",
    RBC = "Erythroid",
    Platelet = "Platelet",
    Cycling = "Cycling"
  )

  list(lineage = lineage, celltype = celltype)
}

cluster_current_annotation_v34 <- function(object, cluster_col, annotation_col) {
  meta <- object[[]]
  cl <- as.character(meta[[cluster_col]])
  cl <- ifelse(startsWith(cl, "g"), cl, paste0("g", cl))
  ann <- as.character(meta[[annotation_col]])

  tab <- as.data.frame(table(cluster = cl, annotation = ann))
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  out <- do.call(rbind, lapply(split(tab, tab$cluster), function(x) {
    x[which.max(x$Freq), , drop = FALSE]
  }))
  rownames(out) <- NULL
  out
}

required_gate_v34 <- function(markers_ev, cluster, label, level, cfg) {
  req_level <- UENO_REQUIRED_MARKERS[[level]]
  req <- req_level[[label]]

  if (is.null(req) || !length(req)) {
    return(list(pass = TRUE, n_required = 0L, n_present = 0L, fraction = NA_real_))
  }

  x <- markers_ev[
    markers_ev$cluster_key == cluster &
      markers_ev$gene %in% req &
      markers_ev$evidence_pass,
    , drop = FALSE
  ]

  present <- unique(x$gene)
  min_n <- min(UENO_REQUIRED_MIN[[level]] %||% 1L, length(req))

  list(
    pass = length(present) >= min_n,
    n_required = length(req),
    n_present = length(present),
    fraction = length(present) / length(req)
  )
}

score_label_evidence_v34 <- function(markers_ev, marker_reference, cluster, level,
                                     label, cfg, prior_label = NA_character_) {
  ref <- marker_reference[
    marker_reference$level == level &
      marker_reference$label == label,
    , drop = FALSE
  ]

  pos <- unique(ref$gene[ref$direction == "positive"])
  neg <- unique(ref$gene[ref$direction == "negative"])

  cx <- markers_ev[markers_ev$cluster_key == cluster, , drop = FALSE]
  pos_x <- cx[cx$gene %in% pos & cx$evidence_pass, , drop = FALSE]
  neg_x <- cx[cx$gene %in% neg & cx$evidence_pass, , drop = FALSE]

  pos_score <- if (nrow(pos_x)) sum(pos_x$gene_evidence) else 0
  neg_score <- if (nrow(neg_x)) sum(neg_x$gene_evidence) else 0

  pos_coverage <- if (length(pos)) length(unique(pos_x$gene)) / length(pos) else NA_real_
  neg_hits <- length(unique(neg_x$gene))

  top20_hits <- length(unique(pos_x$gene[pos_x$marker_rank <= 20]))
  top50_hits <- length(unique(pos_x$gene[pos_x$marker_rank <= 50]))
  top100_hits <- length(unique(pos_x$gene[pos_x$marker_rank <= 100]))

  gate <- required_gate_v34(markers_ev, cluster, label, level, cfg)
  prior_bonus <- if (!is.na(prior_label) && identical(label, prior_label)) {
    cfg$v34_prior_bonus %||% 0.20
  } else {
    0
  }

  exclusion_weight <- cfg$v34_negative_weight %||% 1.25
  gate_penalty <- if (isTRUE(gate$pass)) 0 else (cfg$v34_gate_penalty %||% 2.0)

  score <- pos_score -
    exclusion_weight * neg_score -
    gate_penalty +
    prior_bonus

  data.frame(
    level = level,
    label = label,
    cluster = cluster,
    evidence_score = score,
    positive_evidence = pos_score,
    negative_evidence = neg_score,
    positive_coverage = pos_coverage,
    n_positive_hits = length(unique(pos_x$gene)),
    n_negative_hits = neg_hits,
    top20_hits = top20_hits,
    top50_hits = top50_hits,
    top100_hits = top100_hits,
    required_pass = gate$pass,
    n_required = gate$n_required,
    n_required_present = gate$n_present,
    required_fraction = gate$fraction,
    prior_match = !is.na(prior_label) && identical(label, prior_label),
    supporting_markers = paste(
      utils::head(pos_x$gene[order(-pos_x$gene_evidence)], 15),
      collapse = ";"
    ),
    contradicting_markers = paste(
      utils::head(neg_x$gene[order(-neg_x$gene_evidence)], 15),
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

score_level_v34 <- function(markers_ev, marker_reference, level, clusters, cfg,
                            prior_by_cluster = NULL, eligible_parent = NULL,
                            parent_column = NULL) {
  ref <- marker_reference[marker_reference$level == level, , drop = FALSE]
  labels <- unique(ref$label)

  out <- lapply(clusters, function(cluster) {
    allowed <- labels

    if (!is.null(eligible_parent) && !is.null(parent_column)) {
      parent <- unname(eligible_parent[cluster])
      allowed <- unique(ref$label[ref[[parent_column]] == parent])
    }

    prior_label <- if (is.null(prior_by_cluster)) {
      NA_character_
    } else {
      unname(prior_by_cluster[cluster])
    }

    do.call(rbind, lapply(allowed, function(label) {
      score_label_evidence_v34(
        markers_ev, marker_reference, cluster, level, label, cfg, prior_label
      )
    }))
  })

  do.call(rbind, out)
}

choose_evidence_winners_v34 <- function(scores, cfg) {
  out <- lapply(split(scores, scores$cluster), function(x) {
    x <- x[order(
      -x$evidence_score,
      -x$required_pass,
      -x$top50_hits,
      -x$positive_coverage,
      x$label
    ), , drop = FALSE]

    top <- x[1, , drop = FALSE]
    second <- if (nrow(x) >= 2) x$evidence_score[2] else NA_real_
    top$second_score <- second
    top$score_delta <- top$evidence_score - second

    min_score <- cfg$v34_unknown_score %||% 0.5
    moderate_delta <- cfg$v34_moderate_delta %||% 0.35
    high_delta <- cfg$v34_high_delta %||% 0.8

    top$confidence <- if (
      !isTRUE(top$required_pass) || top$evidence_score < min_score
    ) {
      "Unknown"
    } else if (!is.finite(top$score_delta) || top$score_delta < moderate_delta) {
      "Low"
    } else if (top$score_delta < high_delta) {
      "Moderate"
    } else {
      "High"
    }

    top
  })

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

run_evidence_voting_v34 <- function(object, cluster_col, annotation_col,
                                    markers, marker_reference, cfg) {
  markers_ev <- prepare_marker_evidence_v34(markers, cfg)
  clusters <- sort(unique(markers_ev$cluster_key))

  current <- cluster_current_annotation_v34(object, cluster_col, annotation_col)
  current_map <- setNames(current$annotation, current$cluster)
  prior_maps <- annotation_prior_maps_v34()

  prior_lineage <- setNames(
    unname(prior_maps$lineage[current_map[clusters]]),
    clusters
  )
  prior_celltype <- setNames(
    unname(prior_maps$celltype[current_map[clusters]]),
    clusters
  )

  lineage_scores <- score_level_v34(
    markers_ev, marker_reference, "lineage", clusters, cfg,
    prior_by_cluster = prior_lineage
  )
  lineage_winners <- choose_evidence_winners_v34(lineage_scores, cfg)
  lineage_map <- setNames(lineage_winners$label, lineage_winners$cluster)

  celltype_scores <- score_level_v34(
    markers_ev, marker_reference, "celltype", clusters, cfg,
    prior_by_cluster = prior_celltype,
    eligible_parent = lineage_map,
    parent_column = "parent_lineage"
  )
  celltype_winners <- choose_evidence_winners_v34(celltype_scores, cfg)
  celltype_map <- setNames(celltype_winners$label, celltype_winners$cluster)

  subtype_scores <- score_level_v34(
    markers_ev, marker_reference, "subtype", clusters, cfg,
    eligible_parent = celltype_map,
    parent_column = "parent_celltype"
  )
  subtype_winners <- choose_evidence_winners_v34(subtype_scores, cfg)

  list(
    marker_evidence = markers_ev,
    current_annotation = current,
    scores = rbind(lineage_scores, celltype_scores, subtype_scores),
    winners = rbind(lineage_winners, celltype_winners, subtype_winners),
    lineage = list(scores = lineage_scores, winners = lineage_winners),
    celltype = list(scores = celltype_scores, winners = celltype_winners),
    subtype = list(scores = subtype_scores, winners = subtype_winners)
  )
}

attach_evidence_metadata_v34 <- function(object, cluster_col, voting) {
  meta <- object[[]]
  cl <- as.character(meta[[cluster_col]])
  cl <- ifelse(startsWith(cl, "g"), cl, paste0("g", cl))

  attach_one <- function(meta, winners, prefix) {
    map <- function(field) setNames(winners[[field]], winners$cluster)

    lab <- unname(map("label")[cl])
    conf <- unname(map("confidence")[cl])
    lab[is.na(lab) | conf == "Unknown"] <- "Unknown"

    meta[[prefix]] <- lab
    meta[[paste0(prefix, "_score")]] <- unname(map("evidence_score")[cl])
    meta[[paste0(prefix, "_delta")]] <- unname(map("score_delta")[cl])
    meta[[paste0(prefix, "_confidence")]] <- conf
    meta[[paste0(prefix, "_coverage")]] <- unname(map("positive_coverage")[cl])
    meta[[paste0(prefix, "_required_pass")]] <- unname(map("required_pass")[cl])
    meta[[paste0(prefix, "_top50_hits")]] <- unname(map("top50_hits")[cl])
    meta
  }

  meta <- attach_one(meta, voting$lineage$winners, "vote_ueno_lineage_v34")
  meta <- attach_one(meta, voting$celltype$winners, "vote_ueno_celltype_v34")
  meta <- attach_one(meta, voting$subtype$winners, "vote_ueno_subtype_v34")

  # Recommended label falls back hierarchically.
  subtype <- meta$vote_ueno_subtype_v34
  celltype <- meta$vote_ueno_celltype_v34
  lineage <- meta$vote_ueno_lineage_v34

  recommended <- subtype
  recommended[is.na(recommended) | recommended == "Unknown"] <-
    celltype[is.na(recommended) | recommended == "Unknown"]
  recommended[is.na(recommended) | recommended == "Unknown"] <-
    lineage[is.na(recommended) | recommended == "Unknown"]

  meta$vote_ueno_recommended_v34 <- recommended
  object@meta.data <- meta
  object
}

make_annotation_audit_v34 <- function(voting) {
  current_map <- setNames(
    voting$current_annotation$annotation,
    voting$current_annotation$cluster
  )

  levels <- c("lineage","celltype","subtype")
  out <- lapply(levels, function(level) {
    w <- voting[[level]]$winners
    if (!nrow(w)) return(NULL)

    alternatives <- lapply(split(voting[[level]]$scores, voting[[level]]$scores$cluster), function(x) {
      x <- x[order(-x$evidence_score), , drop = FALSE]
      if (nrow(x) >= 2) x$label[2] else NA_character_
    })
    alt_map <- unlist(alternatives)

    data.frame(
      cluster = w$cluster,
      level = level,
      current_annotation = unname(current_map[w$cluster]),
      recommended_label = ifelse(w$confidence == "Unknown", "Unknown", w$label),
      alternative_label = unname(alt_map[w$cluster]),
      evidence_score = w$evidence_score,
      score_delta = w$score_delta,
      confidence = w$confidence,
      required_pass = w$required_pass,
      required_fraction = w$required_fraction,
      positive_coverage = w$positive_coverage,
      top20_hits = w$top20_hits,
      top50_hits = w$top50_hits,
      top100_hits = w$top100_hits,
      supporting_markers = w$supporting_markers,
      contradicting_markers = w$contradicting_markers,
      recommended_action = ifelse(
        w$confidence %in% c("High","Moderate"),
        "Accept after visual review",
        "Manual review required"
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), out))
}
