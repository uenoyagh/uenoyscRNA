# ==============================================================================
# uenoyscRNA
# Compare candidate Mouse whole-liver RDS objects
#
# Purpose:
#   1. Compare cell-type separation
#   2. Compare within-cell-type sample mixing
#   3. Review annotation-marker consistency
#   4. Select the standard Mouse_MASH_RDS object
#
# Date: 2026-07-30
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "ggplot2",
  "patchwork",
  "cluster",
  "Matrix",
  "RANN",
  "readr",
  "stringr",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing:\n",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(cluster)
  library(Matrix)
  library(RANN)
  library(readr)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

rds_dir <- paste0(
  "/Volumes/SSD990_uenoy/",
  "scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS"
)

output_dir <- file.path(
  "/Volumes/SSD990_uenoy/",
  "scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS_review_20260730"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- c(
  "Mouse_GSE325222_RH251217117_RPCA_integrated_symbolFixed.rds",
  "Mouse_object_with_celltype_annotation_R8tone_metadata.rds",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds",
  "Mouse_RH260519343_GSE325222_RPCA_integrated_celltype_annotated.rds"
)

rds_paths <- file.path(rds_dir, rds_files)

missing_files <- rds_paths[!file.exists(rds_paths)]

if (length(missing_files) > 0) {
  stop(
    "The following RDS files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

object_names <- c(
  "RDS1_symbolFixed",
  "RDS2_R8tone",
  "RDS3_FIXED2",
  "RDS4_RH260519343"
)

names(rds_paths) <- object_names

# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------

max_cells_for_metrics <- 30000L
max_cells_per_celltype <- 5000L
n_pcs_for_metrics <- 30L
knn_k <- 30L

celltype_candidates <- c(
  "celltype",
  "cell_type",
  "CellType",
  "celltype_final",
  "annotation",
  "annotation_final",
  "predicted_celltype",
  "feature_annotation_added"
)

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "orig.ident",
  "donor",
  "replicate"
)

condition_candidates <- c(
  "condition",
  "group",
  "Condition",
  "treatment",
  "diet",
  "experimental_group"
)

cluster_candidates <- c(
  "seurat_clusters",
  "RNA_snn_res.0.5",
  "integrated_snn_res.0.5",
  "SCT_snn_res.0.5",
  "cluster",
  "clusters"
)

# Broad mouse liver marker panel
marker_sets <- list(
  Hepatocyte = c(
    "Alb", "Ttr", "Apoa1", "Apoa2", "Cyp2e1",
    "Cyp3a11", "Ass1", "Asl", "Arg1"
  ),
  Cholangiocyte = c(
    "Krt19", "Krt8", "Krt18", "Sox9", "Epcam",
    "Krt7", "Muc1"
  ),
  LSEC = c(
    "Clec4g", "Stab2", "Lyve1", "Fcgr2b",
    "Kdr", "Pecam1", "Eng"
  ),
  HSC = c(
    "Lrat", "Rgs5", "Reln", "Des", "Col1a1",
    "Col1a2", "Dcn", "Rbp1", "Pdgfrb"
  ),
  Macrophage = c(
    "Adgre1", "Cd68", "C1qa", "C1qb", "C1qc",
    "Lyz2", "Tyrobp", "Aif1"
  ),
  Neutrophil = c(
    "S100a8", "S100a9", "Retnlg", "Lcn2",
    "Csf3r", "Mmp8"
  ),
  T_cell = c(
    "Cd3d", "Cd3e", "Cd3g", "Trbc1",
    "Trbc2", "Lck"
  ),
  B_cell = c(
    "Cd79a", "Cd79b", "Ms4a1", "Cd74",
    "Cd37", "Cd22"
  ),
  NK_cell = c(
    "Nkg7", "Klrd1", "Prf1", "Gzmb",
    "Klrk1", "Xcl1"
  )
)

# ------------------------------------------------------------------------------
# 3. Utility functions
# ------------------------------------------------------------------------------

first_existing_column <- function(metadata, candidates) {
  found <- candidates[candidates %in% colnames(metadata)]

  if (length(found) == 0) {
    return(NA_character_)
  }

  found[[1]]
}

safe_name <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9_.-]", "_")
}

get_assay_for_expression <- function(object) {
  assays <- Assays(object)

  preferred <- c("RNA", "SCT", "integrated")
  found <- preferred[preferred %in% assays]

  if (length(found) == 0) {
    return(DefaultAssay(object))
  }

  found[[1]]
}

get_reduction_for_metrics <- function(object) {
  reductions <- Reductions(object)

  preferred <- c("pca", "rpca", "integrated.rpca", "harmony")

  found <- preferred[preferred %in% reductions]

  if (length(found) == 0) {
    return(NA_character_)
  }

  found[[1]]
}

get_umap_reduction <- function(object) {
  reductions <- Reductions(object)

  exact_preference <- c(
    "umap",
    "integrated.umap",
    "rpca.umap",
    "harmony.umap"
  )

  found_exact <- exact_preference[exact_preference %in% reductions]

  if (length(found_exact) > 0) {
    return(found_exact[[1]])
  }

  found_partial <- reductions[
    stringr::str_detect(tolower(reductions), "umap")
  ]

  if (length(found_partial) == 0) {
    return(NA_character_)
  }

  found_partial[[1]]
}

normalized_entropy <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  tab <- table(x)
  p <- as.numeric(tab) / sum(tab)

  if (length(p) <= 1) {
    return(0)
  }

  entropy <- -sum(p * log(p))
  entropy / log(length(p))
}

inverse_simpson <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  p <- as.numeric(table(x)) / length(x)

  1 / sum(p^2)
}

scaled_inverse_simpson <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  n_groups <- length(unique(x))

  if (n_groups <= 1) {
    return(0)
  }

  score <- inverse_simpson(x)

  (score - 1) / (n_groups - 1)
}

stratified_subsample <- function(
    cells,
    labels,
    max_total = 30000L,
    max_per_group = 5000L,
    seed = 1234
) {
  set.seed(seed)

  df <- tibble(
    cell = cells,
    label = as.character(labels)
  ) |>
    filter(!is.na(label), label != "")

  sampled <- df |>
    group_by(label) |>
    group_modify(
      ~{
        n_take <- min(nrow(.x), max_per_group)
        slice_sample(.x, n = n_take)
      }
    ) |>
    ungroup()

  if (nrow(sampled) > max_total) {
    sampled <- slice_sample(sampled, n = max_total)
  }

  sampled$cell
}

calculate_silhouette <- function(
    embedding,
    labels,
    max_cells = 30000L,
    seed = 1234
) {
  labels <- as.character(labels)

  keep <- !is.na(labels) & labels != ""
  embedding <- embedding[keep, , drop = FALSE]
  labels <- labels[keep]

  if (nrow(embedding) < 10 || length(unique(labels)) < 2) {
    return(
      list(
        overall = NA_real_,
        by_group = tibble()
      )
    )
  }

  set.seed(seed)

  if (nrow(embedding) > max_cells) {
    idx <- sample(seq_len(nrow(embedding)), max_cells)
    embedding <- embedding[idx, , drop = FALSE]
    labels <- labels[idx]
  }

  label_factor <- as.integer(factor(labels))

  sil <- cluster::silhouette(
    label_factor,
    dist(embedding)
  )

  sil_df <- tibble(
    label = labels,
    silhouette = sil[, "sil_width"]
  )

  list(
    overall = mean(sil_df$silhouette, na.rm = TRUE),
    by_group = sil_df |>
      group_by(label) |>
      summarise(
        n_cells = n(),
        mean_silhouette = mean(silhouette, na.rm = TRUE),
        median_silhouette = median(silhouette, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

calculate_knn_celltype_purity <- function(
    embedding,
    labels,
    k = 30L
) {
  labels <- as.character(labels)

  keep <- !is.na(labels) & labels != ""
  embedding <- embedding[keep, , drop = FALSE]
  labels <- labels[keep]

  if (nrow(embedding) <= k + 1 || length(unique(labels)) < 2) {
    return(
      list(
        overall = NA_real_,
        by_group = tibble()
      )
    )
  }

  nn <- RANN::nn2(
    data = embedding,
    query = embedding,
    k = k + 1
  )

  neighbor_index <- nn$nn.idx[, -1, drop = FALSE]

  purity <- vapply(
    seq_len(nrow(neighbor_index)),
    function(i) {
      mean(labels[neighbor_index[i, ]] == labels[i])
    },
    numeric(1)
  )

  purity_df <- tibble(
    label = labels,
    purity = purity
  )

  list(
    overall = mean(purity_df$purity, na.rm = TRUE),
    by_group = purity_df |>
      group_by(label) |>
      summarise(
        n_cells = n(),
        mean_knn_purity = mean(purity, na.rm = TRUE),
        median_knn_purity = median(purity, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

calculate_local_sample_mixing <- function(
    embedding,
    celltypes,
    samples,
    k = 30L
) {
  celltypes <- as.character(celltypes)
  samples <- as.character(samples)

  valid <- (
    !is.na(celltypes) &
      celltypes != "" &
      !is.na(samples) &
      samples != ""
  )

  embedding <- embedding[valid, , drop = FALSE]
  celltypes <- celltypes[valid]
  samples <- samples[valid]

  if (nrow(embedding) <= k + 1) {
    return(
      list(
        overall_entropy = NA_real_,
        overall_inverse_simpson = NA_real_,
        by_celltype = tibble()
      )
    )
  }

  results <- map_dfr(
    unique(celltypes),
    function(ct) {
      idx <- which(celltypes == ct)

      if (
        length(idx) <= k + 1 ||
        length(unique(samples[idx])) < 2
      ) {
        return(
          tibble(
            celltype = ct,
            n_cells = length(idx),
            n_samples = length(unique(samples[idx])),
            local_entropy = NA_real_,
            local_inverse_simpson = NA_real_
          )
        )
      }

      emb_ct <- embedding[idx, , drop = FALSE]
      sample_ct <- samples[idx]

      local_k <- min(k, nrow(emb_ct) - 1L)

      nn <- RANN::nn2(
        data = emb_ct,
        query = emb_ct,
        k = local_k + 1L
      )

      neighbor_index <- nn$nn.idx[, -1, drop = FALSE]

      entropy_values <- apply(
        neighbor_index,
        1,
        function(neighbors) {
          normalized_entropy(sample_ct[neighbors])
        }
      )

      inverse_simpson_values <- apply(
        neighbor_index,
        1,
        function(neighbors) {
          scaled_inverse_simpson(sample_ct[neighbors])
        }
      )

      tibble(
        celltype = ct,
        n_cells = length(idx),
        n_samples = length(unique(sample_ct)),
        local_entropy = mean(entropy_values, na.rm = TRUE),
        local_inverse_simpson = mean(
          inverse_simpson_values,
          na.rm = TRUE
        )
      )
    }
  )

  weighted_mean_safe <- function(x, w) {
    valid <- is.finite(x) & is.finite(w)

    if (!any(valid)) {
      return(NA_real_)
    }

    weighted.mean(x[valid], w[valid])
  }

  list(
    overall_entropy = weighted_mean_safe(
      results$local_entropy,
      results$n_cells
    ),
    overall_inverse_simpson = weighted_mean_safe(
      results$local_inverse_simpson,
      results$n_cells
    ),
    by_celltype = results
  )
}

calculate_sample_silhouette_within_celltype <- function(
    embedding,
    celltypes,
    samples,
    max_cells_per_type = 3000L
) {
  celltypes <- as.character(celltypes)
  samples <- as.character(samples)

  valid <- (
    !is.na(celltypes) &
      celltypes != "" &
      !is.na(samples) &
      samples != ""
  )

  embedding <- embedding[valid, , drop = FALSE]
  celltypes <- celltypes[valid]
  samples <- samples[valid]

  results <- map_dfr(
    unique(celltypes),
    function(ct) {
      idx <- which(celltypes == ct)

      if (
        length(idx) < 20 ||
        length(unique(samples[idx])) < 2
      ) {
        return(
          tibble(
            celltype = ct,
            n_cells = length(idx),
            n_samples = length(unique(samples[idx])),
            sample_silhouette = NA_real_
          )
        )
      }

      set.seed(1234)

      if (length(idx) > max_cells_per_type) {
        idx <- sample(idx, max_cells_per_type)
      }

      emb_ct <- embedding[idx, , drop = FALSE]
      sample_ct <- samples[idx]

      sil <- cluster::silhouette(
        as.integer(factor(sample_ct)),
        dist(emb_ct)
      )

      tibble(
        celltype = ct,
        n_cells = length(idx),
        n_samples = length(unique(sample_ct)),
        sample_silhouette = mean(
          sil[, "sil_width"],
          na.rm = TRUE
        )
      )
    }
  )

  list(
    overall = weighted.mean(
      results$sample_silhouette[
        is.finite(results$sample_silhouette)
      ],
      results$n_cells[
        is.finite(results$sample_silhouette)
      ],
      na.rm = TRUE
    ),
    by_celltype = results
  )
}

normalize_annotation_label <- function(x) {
  x_lower <- tolower(as.character(x))

  case_when(
    str_detect(
      x_lower,
      "hepatocyte|hepatocytes|hepatic"
    ) ~ "Hepatocyte",

    str_detect(
      x_lower,
      "cholangi|bile|ductal"
    ) ~ "Cholangiocyte",

    str_detect(
      x_lower,
      "lsec|sinusoidal endothelial"
    ) ~ "LSEC",

    str_detect(
      x_lower,
      "stellate|hsc|fibroblast"
    ) ~ "HSC",

    str_detect(
      x_lower,
      "macrophage|kupffer|monocyte|mphi"
    ) ~ "Macrophage",

    str_detect(
      x_lower,
      "neutrophil"
    ) ~ "Neutrophil",

    str_detect(
      x_lower,
      "^t$|t cell|t-cell|cd4|cd8"
    ) ~ "T_cell",

    str_detect(
      x_lower,
      "^b$|b cell|b-cell"
    ) ~ "B_cell",

    str_detect(
      x_lower,
      "nk|natural killer"
    ) ~ "NK_cell",

    TRUE ~ NA_character_
  )
}

calculate_marker_scores <- function(
    object,
    celltype_column,
    marker_sets
) {
  assay_use <- get_assay_for_expression(object)
  DefaultAssay(object) <- assay_use

  available_marker_sets <- map(
    marker_sets,
    ~intersect(.x, rownames(object))
  )

  available_marker_sets <- available_marker_sets[
    lengths(available_marker_sets) >= 2
  ]

  if (length(available_marker_sets) == 0) {
    return(
      list(
        object = object,
        summary = tibble(),
        concordance = NA_real_
      )
    )
  }

  for (marker_name in names(available_marker_sets)) {
    score_name <- paste0(
      "ReviewScore_",
      marker_name,
      "_"
    )

    object <- AddModuleScore(
      object = object,
      features = list(
        available_marker_sets[[marker_name]]
      ),
      assay = assay_use,
      name = score_name,
      seed = 1234
    )

    generated_column <- paste0(score_name, "1")

    corrected_column <- paste0(
      "ReviewScore_",
      marker_name
    )

    object[[corrected_column]] <- object[[generated_column]][, 1]
    object[[generated_column]] <- NULL
  }

  score_columns <- paste0(
    "ReviewScore_",
    names(available_marker_sets)
  )

  metadata <- object[[]]

  score_matrix <- as.matrix(
    metadata[, score_columns, drop = FALSE]
  )

  predicted_index <- max.col(
    score_matrix,
    ties.method = "first"
  )

  predicted_type <- names(available_marker_sets)[predicted_index]

  annotated_type <- normalize_annotation_label(
    metadata[[celltype_column]]
  )

  comparable <- (
    !is.na(annotated_type) &
      annotated_type %in% names(available_marker_sets)
  )

  concordance <- if (sum(comparable) == 0) {
    NA_real_
  } else {
    mean(
      predicted_type[comparable] ==
        annotated_type[comparable]
    )
  }

  summary <- tibble(
    annotated_celltype = as.character(
      metadata[[celltype_column]]
    ),
    normalized_annotation = annotated_type,
    predicted_marker_type = predicted_type
  ) |>
    count(
      annotated_celltype,
      normalized_annotation,
      predicted_marker_type,
      name = "n_cells"
    ) |>
    group_by(annotated_celltype) |>
    mutate(
      fraction_within_annotation =
        n_cells / sum(n_cells)
    ) |>
    ungroup()

  list(
    object = object,
    summary = summary,
    concordance = concordance
  )
}

make_metadata_summary <- function(
    seu,
    object_name,
    file_path
) {

  metadata <- seu[[]]

  celltype_column <- first_existing_column(
    metadata,
    celltype_candidates
  )

  sample_column <- first_existing_column(
    metadata,
    sample_candidates
  )

  condition_column <- first_existing_column(
    metadata,
    condition_candidates
  )

  cluster_column <- first_existing_column(
    metadata,
    cluster_candidates
  )

  tibble(
    object_name = object_name,
    file = basename(file_path),
    file_size_gb = file.info(file_path)$size / 1024^3,
    n_cells = ncol(seu),
    n_features = nrow(seu),
    default_assay = DefaultAssay(seu),
    assays = paste(Assays(seu), collapse = "; "),
    reductions = paste(Reductions(seu), collapse = "; "),
    celltype_column = celltype_column,
    sample_column = sample_column,
    condition_column = condition_column,
    cluster_column = cluster_column,
    n_celltypes = if (is.na(celltype_column)) {
      NA_integer_
    } else {
      dplyr::n_distinct(metadata[[celltype_column]])
    },
    n_samples = if (is.na(sample_column)) {
      NA_integer_
    } else {
      dplyr::n_distinct(metadata[[sample_column]])
    },
    n_conditions = if (is.na(condition_column)) {
      NA_integer_
    } else {
      dplyr::n_distinct(metadata[[condition_column]])
    },
    metric_reduction = get_reduction_for_metrics(seu),
    umap_reduction = get_umap_reduction(seu)
  )
}

# ------------------------------------------------------------------------------
# 4. Main evaluation
# ------------------------------------------------------------------------------

all_inventory <- list()
all_cell_counts <- list()
all_sample_composition <- list()
all_silhouette_by_type <- list()
all_knn_by_type <- list()
all_mixing_by_type <- list()
all_sample_silhouette_by_type <- list()
all_marker_summary <- list()
all_scores <- list()
all_umap_plots <- list()

for (object_name in names(rds_paths)) {

  message(
    "\n============================================================"
  )
  message("Loading: ", object_name)
  message("File: ", rds_paths[[object_name]])
  message(
    "============================================================"
  )

  object <- readRDS(rds_paths[[object_name]])

  if (!inherits(object, "Seurat")) {
    stop(
      object_name,
      " is not a Seurat object."
    )
  }

  metadata <- object[[]]

  celltype_column <- first_existing_column(
    metadata,
    celltype_candidates
  )

  sample_column <- first_existing_column(
    metadata,
    sample_candidates
  )

  condition_column <- first_existing_column(
    metadata,
    condition_candidates
  )

  inventory <- make_metadata_summary(
    seu = object,
    object_name = object_name,
    file_path = rds_paths[[object_name]]
  )

  all_inventory[[object_name]] <- inventory

  # --------------------------------------------------------------------------
  # Cell counts
  # --------------------------------------------------------------------------

  if (!is.na(celltype_column)) {
    cell_counts <- metadata |>
      rownames_to_column("cell") |>
      count(
        celltype = .data[[celltype_column]],
        name = "n_cells"
      ) |>
      mutate(
        object = object_name,
        fraction = n_cells / sum(n_cells)
      )

    all_cell_counts[[object_name]] <- cell_counts
  }

  if (
    !is.na(celltype_column) &&
    !is.na(sample_column)
  ) {
    sample_composition <- metadata |>
      rownames_to_column("cell") |>
      count(
        sample = .data[[sample_column]],
        celltype = .data[[celltype_column]],
        name = "n_cells"
      ) |>
      group_by(sample) |>
      mutate(
        fraction_within_sample =
          n_cells / sum(n_cells)
      ) |>
      ungroup() |>
      mutate(object = object_name)

    all_sample_composition[[object_name]] <-
      sample_composition
  }

  # --------------------------------------------------------------------------
  # Embedding-based metrics
  # --------------------------------------------------------------------------

  metric_reduction <- get_reduction_for_metrics(object)

  if (
    !is.na(metric_reduction) &&
    !is.na(celltype_column)
  ) {
    embedding <- Embeddings(
      object,
      reduction = metric_reduction
    )

    n_dims <- min(
      ncol(embedding),
      n_pcs_for_metrics
    )

    embedding <- embedding[, seq_len(n_dims), drop = FALSE]

    evaluation_cells <- stratified_subsample(
      cells = rownames(embedding),
      labels = metadata[
        rownames(embedding),
        celltype_column
      ],
      max_total = max_cells_for_metrics,
      max_per_group = max_cells_per_celltype
    )

    embedding_eval <- embedding[
      evaluation_cells,
      ,
      drop = FALSE
    ]

    celltype_eval <- metadata[
      evaluation_cells,
      celltype_column
    ]

    celltype_silhouette <- calculate_silhouette(
      embedding = embedding_eval,
      labels = celltype_eval,
      max_cells = max_cells_for_metrics
    )

    all_silhouette_by_type[[object_name]] <-
      celltype_silhouette$by_group |>
      mutate(object = object_name)

    knn_purity <- calculate_knn_celltype_purity(
      embedding = embedding_eval,
      labels = celltype_eval,
      k = knn_k
    )

    all_knn_by_type[[object_name]] <-
      knn_purity$by_group |>
      mutate(object = object_name)

    sample_mixing_entropy <- NA_real_
    sample_mixing_inverse_simpson <- NA_real_
    sample_silhouette_value <- NA_real_

    if (!is.na(sample_column)) {
      sample_eval <- metadata[
        evaluation_cells,
        sample_column
      ]

      sample_mixing <- calculate_local_sample_mixing(
        embedding = embedding_eval,
        celltypes = celltype_eval,
        samples = sample_eval,
        k = knn_k
      )

      sample_mixing_entropy <-
        sample_mixing$overall_entropy

      sample_mixing_inverse_simpson <-
        sample_mixing$overall_inverse_simpson

      all_mixing_by_type[[object_name]] <-
        sample_mixing$by_celltype |>
        mutate(object = object_name)

      sample_silhouette <-
        calculate_sample_silhouette_within_celltype(
          embedding = embedding_eval,
          celltypes = celltype_eval,
          samples = sample_eval
        )

      sample_silhouette_value <-
        sample_silhouette$overall

      all_sample_silhouette_by_type[[object_name]] <-
        sample_silhouette$by_celltype |>
        mutate(object = object_name)
    }

  } else {
    celltype_silhouette <- list(overall = NA_real_)
    knn_purity <- list(overall = NA_real_)
    sample_mixing_entropy <- NA_real_
    sample_mixing_inverse_simpson <- NA_real_
    sample_silhouette_value <- NA_real_
  }

  # --------------------------------------------------------------------------
  # Marker consistency
  # --------------------------------------------------------------------------

  marker_concordance <- NA_real_

  if (!is.na(celltype_column)) {
    marker_result <- calculate_marker_scores(
      object = object,
      celltype_column = celltype_column,
      marker_sets = marker_sets
    )

    marker_concordance <- marker_result$concordance

    if (nrow(marker_result$summary) > 0) {
      all_marker_summary[[object_name]] <-
        marker_result$summary |>
        mutate(object = object_name)
    }
  }

  # --------------------------------------------------------------------------
  # Score summary
  # --------------------------------------------------------------------------

  all_scores[[object_name]] <- tibble(
    object = object_name,
    file = basename(rds_paths[[object_name]]),
    celltype_silhouette =
      celltype_silhouette$overall,
    knn_celltype_purity =
      knn_purity$overall,
    sample_mixing_entropy =
      sample_mixing_entropy,
    sample_mixing_inverse_simpson =
      sample_mixing_inverse_simpson,
    sample_silhouette_within_celltype =
      sample_silhouette_value,
    marker_annotation_concordance =
      marker_concordance
  )

  # --------------------------------------------------------------------------
  # UMAP plots
  # --------------------------------------------------------------------------

  umap_reduction <- get_umap_reduction(object)

  if (!is.na(umap_reduction)) {
    plot_list <- list()

    if (!is.na(celltype_column)) {
      plot_list$celltype <- DimPlot(
        object,
        reduction = umap_reduction,
        group.by = celltype_column,
        raster = ncol(object) > 100000,
        pt.size = 0.2
      ) +
        ggtitle(
          paste0(
            object_name,
            "\nCell type: ",
            celltype_column
          )
        ) +
        theme(
          plot.title = element_text(
            size = 12,
            face = "bold"
          ),
          legend.text = element_text(size = 8)
        )
    }

    if (!is.na(sample_column)) {
      plot_list$sample <- DimPlot(
        object,
        reduction = umap_reduction,
        group.by = sample_column,
        raster = ncol(object) > 100000,
        pt.size = 0.2
      ) +
        ggtitle(
          paste0(
            object_name,
            "\nSample: ",
            sample_column
          )
        ) +
        theme(
          plot.title = element_text(
            size = 12,
            face = "bold"
          ),
          legend.text = element_text(size = 8)
        )
    }

    if (!is.na(condition_column)) {
      plot_list$condition <- DimPlot(
        object,
        reduction = umap_reduction,
        group.by = condition_column,
        raster = ncol(object) > 100000,
        pt.size = 0.2
      ) +
        ggtitle(
          paste0(
            object_name,
            "\nCondition: ",
            condition_column
          )
        ) +
        theme(
          plot.title = element_text(
            size = 12,
            face = "bold"
          ),
          legend.text = element_text(size = 8)
        )
    }

    if (length(plot_list) > 0) {
      combined_plot <- wrap_plots(
        plot_list,
        ncol = 1
      )

      ggsave(
        filename = file.path(
          output_dir,
          paste0(
            "UMAP_review_",
            safe_name(object_name),
            ".pdf"
          )
        ),
        plot = combined_plot,
        width = 11,
        height = max(7, 6 * length(plot_list)),
        device = cairo_pdf
      )

      all_umap_plots[[object_name]] <- plot_list
    }
  }

  rm(object)
  invisible(gc())
}

# ------------------------------------------------------------------------------
# 5. Combine results
# ------------------------------------------------------------------------------

inventory_table <- bind_rows(all_inventory)
score_table <- bind_rows(all_scores)

cell_count_table <- bind_rows(all_cell_counts)
sample_composition_table <- bind_rows(all_sample_composition)
silhouette_by_type_table <- bind_rows(all_silhouette_by_type)
knn_by_type_table <- bind_rows(all_knn_by_type)
mixing_by_type_table <- bind_rows(all_mixing_by_type)
sample_silhouette_by_type_table <-
  bind_rows(all_sample_silhouette_by_type)
marker_summary_table <- bind_rows(all_marker_summary)

# ------------------------------------------------------------------------------
# 6. Ranking
# ------------------------------------------------------------------------------

scale_metric <- function(x, higher_is_better = TRUE) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  valid <- is.finite(x)

  if (sum(valid) <= 1 || diff(range(x[valid])) == 0) {
    result <- rep(0.5, length(x))
    result[!valid] <- NA_real_
    return(result)
  }

  result <- rep(NA_real_, length(x))
  result[valid] <- scales::rescale(
    x[valid],
    to = c(0, 1)
  )

  if (!higher_is_better) {
    result[valid] <- 1 - result[valid]
  }

  result
}

ranked_score_table <- score_table |>
  mutate(
    score_celltype_silhouette =
      scale_metric(
        celltype_silhouette,
        higher_is_better = TRUE
      ),

    score_knn_purity =
      scale_metric(
        knn_celltype_purity,
        higher_is_better = TRUE
      ),

    score_sample_entropy =
      scale_metric(
        sample_mixing_entropy,
        higher_is_better = TRUE
      ),

    score_sample_inverse_simpson =
      scale_metric(
        sample_mixing_inverse_simpson,
        higher_is_better = TRUE
      ),

    score_sample_silhouette =
      scale_metric(
        abs(sample_silhouette_within_celltype),
        higher_is_better = FALSE
      ),

    score_marker_concordance =
      scale_metric(
        marker_annotation_concordance,
        higher_is_better = TRUE
      )
  ) |>
  rowwise() |>
  mutate(
    separation_score = mean(
      c(
        score_celltype_silhouette,
        score_knn_purity,
        score_marker_concordance
      ),
      na.rm = TRUE
    ),

    batch_mixing_score = mean(
      c(
        score_sample_entropy,
        score_sample_inverse_simpson,
        score_sample_silhouette
      ),
      na.rm = TRUE
    ),

    total_score = (
      0.60 * separation_score +
        0.40 * batch_mixing_score
    )
  ) |>
  ungroup() |>
  arrange(desc(total_score)) |>
  mutate(rank = row_number())

# ------------------------------------------------------------------------------
# 7. Save CSV files
# ------------------------------------------------------------------------------

write_csv(
  inventory_table,
  file.path(output_dir, "01_RDS_inventory.csv")
)

write_csv(
  score_table,
  file.path(output_dir, "02_raw_metric_summary.csv")
)

write_csv(
  ranked_score_table,
  file.path(output_dir, "03_RDS_final_ranking.csv")
)

if (nrow(cell_count_table) > 0) {
  write_csv(
    cell_count_table,
    file.path(output_dir, "04_celltype_cell_counts.csv")
  )
}

if (nrow(sample_composition_table) > 0) {
  write_csv(
    sample_composition_table,
    file.path(output_dir, "05_sample_celltype_composition.csv")
  )
}

if (nrow(silhouette_by_type_table) > 0) {
  write_csv(
    silhouette_by_type_table,
    file.path(output_dir, "06_celltype_silhouette_by_type.csv")
  )
}

if (nrow(knn_by_type_table) > 0) {
  write_csv(
    knn_by_type_table,
    file.path(output_dir, "07_knn_purity_by_type.csv")
  )
}

if (nrow(mixing_by_type_table) > 0) {
  write_csv(
    mixing_by_type_table,
    file.path(output_dir, "08_sample_mixing_by_celltype.csv")
  )
}

if (nrow(sample_silhouette_by_type_table) > 0) {
  write_csv(
    sample_silhouette_by_type_table,
    file.path(output_dir, "09_sample_silhouette_by_celltype.csv")
  )
}

if (nrow(marker_summary_table) > 0) {
  write_csv(
    marker_summary_table,
    file.path(output_dir, "10_marker_annotation_concordance.csv")
  )
}

# ------------------------------------------------------------------------------
# 8. Comparison plots
# ------------------------------------------------------------------------------

plot_metric_long <- ranked_score_table |>
  select(
    object,
    celltype_silhouette,
    knn_celltype_purity,
    sample_mixing_entropy,
    sample_mixing_inverse_simpson,
    sample_silhouette_within_celltype,
    marker_annotation_concordance
  ) |>
  pivot_longer(
    cols = -object,
    names_to = "metric",
    values_to = "value"
  )

p_metrics <- ggplot(
  plot_metric_long,
  aes(
    x = reorder(object, value),
    y = value
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(
    ~metric,
    scales = "free_x",
    ncol = 2
  ) +
  labs(
    x = NULL,
    y = "Raw metric value",
    title = "Mouse RDS comparison: raw metrics"
  ) +
  theme_bw(base_size = 11)

ggsave(
  filename = file.path(
    output_dir,
    "11_RDS_raw_metric_comparison.pdf"
  ),
  plot = p_metrics,
  width = 12,
  height = 10,
  device = cairo_pdf
)

plot_score_long <- ranked_score_table |>
  select(
    object,
    separation_score,
    batch_mixing_score,
    total_score
  ) |>
  pivot_longer(
    cols = -object,
    names_to = "score_type",
    values_to = "score"
  )

p_scores <- ggplot(
  plot_score_long,
  aes(
    x = reorder(object, score),
    y = score
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(
    ~score_type,
    ncol = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Normalized score",
    title = "Mouse RDS comparison: final scores"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(
    output_dir,
    "12_RDS_final_score_comparison.pdf"
  ),
  plot = p_scores,
  width = 9,
  height = 10,
  device = cairo_pdf
)

# Cell-type-specific comparison
if (
  nrow(silhouette_by_type_table) > 0 &&
  nrow(knn_by_type_table) > 0
) {
  per_type_metrics <- full_join(
    silhouette_by_type_table,
    knn_by_type_table,
    by = c(
      "object",
      "label",
      "n_cells"
    )
  )

  write_csv(
    per_type_metrics,
    file.path(
      output_dir,
      "13_celltype_specific_separation_summary.csv"
    )
  )

  p_celltype_silhouette <- ggplot(
    silhouette_by_type_table,
    aes(
      x = object,
      y = mean_silhouette
    )
  ) +
    geom_point(size = 2) +
    facet_wrap(
      ~label,
      scales = "free_y"
    ) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Mean silhouette",
      title = "Cell-type-specific separation"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    filename = file.path(
      output_dir,
      "14_celltype_specific_silhouette.pdf"
    ),
    plot = p_celltype_silhouette,
    width = 14,
    height = 10,
    device = cairo_pdf
  )

  p_celltype_knn <- ggplot(
    knn_by_type_table,
    aes(
      x = object,
      y = mean_knn_purity
    )
  ) +
    geom_point(size = 2) +
    facet_wrap(
      ~label,
      scales = "free_y"
    ) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Mean kNN purity",
      title = "Cell-type-specific local purity"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    filename = file.path(
      output_dir,
      "15_celltype_specific_knn_purity.pdf"
    ),
    plot = p_celltype_knn,
    width = 14,
    height = 10,
    device = cairo_pdf
  )
}

# ------------------------------------------------------------------------------
# 9. Console report
# ------------------------------------------------------------------------------

message("\n============================================================")
message("RDS comparison completed")
message("============================================================")
message("Output directory:")
message(output_dir)

message("\nFinal ranking:")
print(
  ranked_score_table |>
    select(
      rank,
      object,
      file,
      separation_score,
      batch_mixing_score,
      total_score
    ),
  n = Inf
)

message(
  "\nImportant: the final RDS must not be selected ",
  "from total_score alone."
)

message(
  "Review the UMAP PDFs and cell-type-specific metrics, ",
  "especially hepatocyte, LSEC, HSC and macrophage clusters."
)
