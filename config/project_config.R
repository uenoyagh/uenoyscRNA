# ============================================================
# uenoy scRNAseq Framework
# Project configuration
# Version 2.5.0
# ============================================================

project_config <- list(
  framework_name = "uenoy scRNAseq Framework",
  framework_version = "2.5.0",
  script_project_root = "/Users/uenoya/Projects/uenoyscRNA",
  external_data_root = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
  Mouse_MASH_RDS = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS",
  Mouse_MASH_Mphi_RDS = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS",
  Human_MASH_RDS = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS",
  Human_MASH_Mphi_RDS = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_Mphi_RDS"
)

allowed_datasets <- c(
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_RDS",
  "Human_MASH_RDS",
  "Human_MASH_Mphi_RDS"
)

get_dataset_dir <- function(dataset_name) {
  if (!dataset_name %in% allowed_datasets) {
    stop(
      "Unknown dataset_name: ", dataset_name,
      "\nAllowed values: ", paste(allowed_datasets, collapse = ", ")
    )
  }
  project_config[[dataset_name]]
}

get_result_root <- function(dataset_name) {
  file.path(get_dataset_dir(dataset_name), "analysis_results")
}

get_result_dir <- function(dataset_name, result_type = NULL, create = TRUE) {
  root <- get_result_root(dataset_name)

  result_map <- c(
    inventory = "00_inventory",
    umap = "01_UMAP",
    composition = "02_cell_composition",
    cluster_highlight = "03_cluster_highlight",
    marker_dynamics = "04_marker_dynamics",
    module_score = "05_module_score",
    de = "06_DE_analysis",
    tables = "07_tables",
    logs = "08_logs"
  )

  out <- if (is.null(result_type)) {
    root
  } else {
    if (!result_type %in% names(result_map)) {
      stop("Unknown result_type: ", result_type)
    }
    file.path(root, unname(result_map[[result_type]]))
  }

  if (isTRUE(create)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
  }

  out
}
