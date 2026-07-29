# ============================================================
# 06_functional_signature_heatmap.R
# ============================================================

rm(list = ls())
gc()

script_root <- "/Users/uenoya/Projects/uenoyscRNA"

source(file.path(script_root, "config", "project_config.R"))
source(file.path(script_root, "config", "local_config.R"))
source(file.path(script_root, "config", "signatures_mouse_macrophage.R"))
source(file.path(script_root, "R", "io.R"))
source(file.path(script_root, "R", "utils.R"))
source(file.path(script_root, "R", "umap.R"))
source(file.path(script_root, "R", "signature_heatmap.R"))

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

target_dir <- get_dataset_dir(analysis_target)
output_dir <- file.path(
  get_result_dir(analysis_target, "module_score", create = TRUE),
  "functional_signature_heatmap"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list_rds_files(target_dir, recursive = FALSE)

if (!is.null(signature_heatmap_selected_files)) {
  rds_files <- rds_files[
    basename(rds_files) %in% signature_heatmap_selected_files
  ]
}

for (rds_path in rds_files) {
  rds_name <- basename(rds_path)
  file_stub <- safe_umap_filename(rds_name)

  object <- safe_read_rds(rds_path)
  if (inherits(object, "rds_read_error")) next
  if (!inherits(object, "Seurat")) next

  cluster_column <- resolve_cluster_column(
    object,
    signature_heatmap_cluster_column_override
  )

  if (is.na(cluster_column)) {
    warning("No cluster column: ", rds_name)
    rm(object)
    next
  }

  assay <- if (is.null(signature_heatmap_assay_override)) {
    SeuratObject::DefaultAssay(object)
  } else {
    signature_heatmap_assay_override
  }

  result <- calculate_signature_cluster_scores(
    object = object,
    cluster_column = cluster_column,
    signatures = mouse_macrophage_signatures,
    assay = assay,
    min_genes = signature_heatmap_min_genes
  )

  p <- make_signature_heatmap(
    result$scores,
    clip = signature_heatmap_clip
  )

  pdf_path <- file.path(
    output_dir,
    paste0(
      file_stub,
      "_functional_signature_heatmap_",
      cluster_column,
      ".pdf"
    )
  )

  grDevices::pdf(
    pdf_path,
    width = signature_heatmap_width,
    height = signature_heatmap_height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()

  safe_write_csv(
    result$scores,
    file.path(
      output_dir,
      paste0(file_stub, "_functional_signature_scores.csv")
    )
  )

  safe_write_csv(
    result$gene_report,
    file.path(
      output_dir,
      paste0(file_stub, "_functional_signature_gene_report.csv")
    )
  )

  rm(object)
  gc(verbose = FALSE)
}
