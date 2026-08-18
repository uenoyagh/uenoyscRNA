# ============================================================
# uenoyscRNA Mouse macrophage RPCA cluster→Layer workflow v2.0.2
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

UenoMphiLayer1_v2 <- c(
  "Resident Kupffer-like Mphi" = "#168A4A",
  "Monocyte-like Mphi" = "#00A67A",
  "Inflammatory M1-like Mphi" = "#E53935",
  "Pro-resolution M2-like Mphi" = "#2E8B57",
  "SPP1/TREM2 MASH-associated Mphi" = "#B02A9B",
  "Other" = "#7A7A7A"
)

UenoMphiLayer1_markers_v2 <- list(
  "Resident Kupffer-like Mphi" = c(
    "Clec4f","Timd4","Vsig4","Marco","Cd5l","C1qa","C1qb","C1qc","Fcgr2b","Slc40a1"
  ),
  "Monocyte-like Mphi" = c(
    "Lyz2","Ccr2","Ly6c2","S100a8","S100a9","Plac8","Ctss","Fcgr3","Ms4a7","Tyrobp"
  ),
  "Inflammatory M1-like Mphi" = c(
    "Il1b","Tnf","Ccl2","Ccl3","Ccl4","Cxcl2","Cxcl9","Cxcl10","Nos2","Ptgs2","Cd80","Cd86","Irf1","Stat1"
  ),
  "Pro-resolution M2-like Mphi" = c(
    "Mrc1","Cd163","Il10","Maf","Mafb","Klf4","Arg1","Retnla","Chil3","Ccl24","Gas6","Axl","Mertk"
  ),
  "SPP1/TREM2 MASH-associated Mphi" = c(
    "Spp1","Trem2","Gpnmb","Cd9","Lgals3","Lpl","Fabp5","Ctsb","Ctsd","Ctsk","Apoe","Cst7","Itgax"
  )
)

detect_reduction_mphi_v200 <- function(object) {
  candidates <- c("rpca_umap_mphi","umapRPCA","umap_rpca","rpca.umap","umap","UMAP")
  hit <- candidates[candidates %in% names(object@reductions)]
  if (length(hit) > 0) return(hit[[1]])
  hit <- names(object@reductions)[grepl("umap", names(object@reductions), ignore.case=TRUE)]
  if (length(hit) > 0) return(hit[[1]])
  stop("No UMAP reduction found.")
}

detect_cluster_col_mphi_v200 <- function(object) {
  candidates <- c(
    "seurat_clusters","cluster","cluster_id","RNA_snn_res.3","RNA_snn_res.3.0",
    "integrated_snn_res.3","integrated_snn_res.3.0","res3.0","cluster_for_plot"
  )
  hit <- candidates[candidates %in% colnames(object[[]])]
  if (length(hit) == 0) stop("No cluster column found.")
  hit[[1]]
}

resolve_features_ci_v200 <- function(features, available) {
  idx <- match(tolower(features), tolower(available))
  unique(available[idx[!is.na(idx)]])
}

normalize_group_names_v200 <- function(x, original) {
  y <- as.character(x)
  numeric_prefixed <- grepl("^g[-+]?[0-9]+([.][0-9]+)?$", y)
  stripped <- sub("^g", "", y[numeric_prefixed])
  ok <- stripped %in% unique(as.character(original))
  y[numeric_prefixed][ok] <- stripped[ok]
  y
}


join_assay_layers_for_de_v201 <- function(object, assay = NULL) {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  assay_object <- object[[assay]]

  if (!inherits(assay_object, "Assay5")) {
    message("Assay ", assay, " is not Assay5; JoinLayers is not required.")
    return(object)
  }

  layers_before <- SeuratObject::Layers(assay_object)

  message(
    "Assay layers before JoinLayers: ",
    paste(layers_before, collapse = ", ")
  )

  split_data_layers <- layers_before[
    grepl("^data[.]", layers_before)
  ]

  split_count_layers <- layers_before[
    grepl("^counts[.]", layers_before)
  ]

  needs_join <- (
    length(split_data_layers) > 0L ||
    length(split_count_layers) > 0L
  )

  if (needs_join) {
    object <- SeuratObject::JoinLayers(
      object = object,
      assay = assay
    )
  }

  layers_after <- SeuratObject::Layers(object[[assay]])

  message(
    "Assay layers after JoinLayers: ",
    paste(layers_after, collapse = ", ")
  )

  if (!"data" %in% layers_after) {
    message(
      "Joined assay has no normalized data layer. Running NormalizeData()."
    )
    object <- Seurat::NormalizeData(
      object = object,
      assay = assay,
      verbose = FALSE
    )
  }

  object
}

find_cluster_markers_v200 <- function(
    object, cluster_col, assay=NULL, slot="data",
    min_pct=0.10, logfc_threshold=0.10
) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  Seurat::Idents(object) <- cluster_col
  markers <- Seurat::FindAllMarkers(
    object=object, assay=assay, slot=slot, only.pos=TRUE,
    min.pct=min_pct, logfc.threshold=logfc_threshold,
    test.use="wilcox", verbose=TRUE
  )
  if (nrow(markers) == 0L) {
    stop(
      "FindAllMarkers returned zero rows after JoinLayers. ",
      "Check the assay, normalized data layer, and cluster identities."
    )
  }

  fc_col <- intersect(c("avg_log2FC","avg_logFC","avg_diff"), colnames(markers))
  if (length(fc_col)==0) stop("Fold-change column not found.")
  fc_col <- fc_col[[1]]
  markers <- markers[
    !grepl("^(mt-|Mt-|Rpl|Rps|Hba|Hbb)|^Malat1$|^Xist$|^Jun$|^Fos$", markers$gene),
    , drop=FALSE
  ]
  markers <- markers[order(markers$cluster, -markers[[fc_col]], markers$p_val_adj), , drop=FALSE]
  attr(markers, "fc_col") <- fc_col
  markers
}

top_n_by_cluster_v200 <- function(markers, n=30L) {
  out <- do.call(rbind, lapply(split(markers, markers$cluster), function(df) head(df, n)))
  rownames(out) <- NULL
  out
}

score_clusters_layer1_v200 <- function(object, cluster_col, assay=NULL, slot="data") {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  available <- rownames(object)
  resolved <- lapply(UenoMphiLayer1_markers_v2, resolve_features_ci_v200, available=available)
  genes <- unique(unlist(resolved, use.names=FALSE))
  avg <- Seurat::AverageExpression(
    object=object, assays=assay, features=genes, group.by=cluster_col, slot=slot, verbose=FALSE
  )[[assay]]
  avg <- as.matrix(avg)
  colnames(avg) <- normalize_group_names_v200(colnames(avg), object[[]][[cluster_col]])
  z <- t(scale(t(avg))); z[!is.finite(z)] <- 0

  rows <- list()
  for (cl in colnames(z)) {
    for (lay in names(resolved)) {
      g <- resolved[[lay]]
      rows[[length(rows)+1]] <- data.frame(
        cluster=cl,
        layer1=lay,
        score=if(length(g)>0) mean(z[g,cl,drop=TRUE],na.rm=TRUE) else NA_real_,
        coverage=length(g)/length(UenoMphiLayer1_markers_v2[[lay]]),
        genes_found=paste(g,collapse=";"),
        stringsAsFactors=FALSE
      )
    }
  }
  score_table <- do.call(rbind, rows)

  provisional <- do.call(rbind, lapply(split(score_table, score_table$cluster), function(df){
    df <- df[order(-df$score),,drop=FALSE]
    delta <- if(nrow(df)>=2) df$score[1]-df$score[2] else NA_real_
    label <- if(is.na(df$score[1]) || df$score[1] < 0.20 || is.na(delta) || delta < 0.10 || df$coverage[1] < 0.20) {
      "Other"
    } else df$layer1[1]
    data.frame(
      cluster=df$cluster[1],
      provisional_layer1=label,
      top_score=df$score[1],
      score_delta=delta,
      marker_coverage=df$coverage[1],
      stringsAsFactors=FALSE
    )
  }))
  rownames(provisional) <- NULL

  list(score_table=score_table, provisional=provisional)
}

apply_manual_override_v200 <- function(provisional, override_file) {
  if (!file.exists(override_file)) return(provisional)
  ov <- utils::read.csv(override_file, stringsAsFactors=FALSE, check.names=FALSE)
  needed <- c("cluster","manual_layer1")
  if (!all(needed %in% colnames(ov))) {
    stop("Override CSV must contain: cluster, manual_layer1")
  }
  valid <- c(names(UenoMphiLayer1_v2))
  ov$manual_layer1[is.na(ov$manual_layer1)] <- ""
  ov$manual_layer1 <- trimws(as.character(ov$manual_layer1))

  nonempty_manual <- ov$manual_layer1[nzchar(ov$manual_layer1)]
  bad <- setdiff(unique(nonempty_manual), valid)

  if (length(bad) > 0L) {
    stop(
      "Invalid manual_layer1 values: ",
      paste(bad, collapse = ", ")
    )
  }

  map <- setNames(ov$manual_layer1, as.character(ov$cluster))
  provisional$final_layer1 <- provisional$provisional_layer1
  hit <- provisional$cluster %in% names(map) & nzchar(map[provisional$cluster])
  provisional$final_layer1[hit] <- unname(map[provisional$cluster[hit]])
  provisional
}

apply_cluster_labels_v200 <- function(object, cluster_col, annotation_table) {
  md <- object[[]]
  cl <- as.character(md[[cluster_col]])
  map <- setNames(annotation_table$final_layer1, annotation_table$cluster)
  layer <- unname(map[cl])
  if (any(is.na(layer))) {
    stop("Unmapped clusters: ", paste(sort(unique(cl[is.na(layer)])),collapse=", "))
  }
  object[["mphi_layer1_v200"]] <- factor(layer, levels=names(UenoMphiLayer1_v2))
  object[["mphi_comprehensive_v200"]] <- factor(
    paste0("C",cl," | ",layer),
    levels=unique(paste0("C",annotation_table$cluster," | ",annotation_table$final_layer1))
  )
  object
}

cluster_palette_by_layer_v200 <- function(annotation_table) {
  base <- UenoMphiLayer1_v2[annotation_table$final_layer1]
  names(base) <- paste0("C", annotation_table$cluster, " | ", annotation_table$final_layer1)
  unname_base <- unname(base)
  names(unname_base) <- names(base)
  unname_base
}

publish_umap_v200 <- function(
    object, group_by, reduction, palette, title,
    pt_size=1.10, label_size=4.5, source_rds=NULL, created_at=Sys.time(),
    legend_text_size=8.5
) {
  p <- Seurat::DimPlot(
    object=object, reduction=reduction, group.by=group_by,
    pt.size=pt_size, label=TRUE, repel=TRUE, label.size=label_size,
    raster=TRUE, raster.dpi=c(700,700)
  )
  if(length(p$layers)>0) for(i in seq_along(p$layers)) p$layers[[i]]$aes_params$alpha <- 1
  p +
    ggplot2::scale_color_manual(values=palette, drop=FALSE) +
    ggplot2::coord_fixed(ratio=1) +
    ggplot2::labs(
      title=title, color=NULL,
      caption=paste0(
        "Source RDS: ",basename(source_rds),
        " | Palette: UenoMphiLayer1_v2",
        " | Created: ",format(created_at,"%Y-%m-%d %H:%M:%S %Z")
      )
    ) +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(hjust=.5,face="bold",size=15),
      legend.text=ggplot2::element_text(size=legend_text_size),
      legend.key.height=grid::unit(.42,"cm"),
      plot.caption=ggplot2::element_text(hjust=0,size=7.5,colour="grey25")
    ) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(size=4.2,alpha=1)))
}
