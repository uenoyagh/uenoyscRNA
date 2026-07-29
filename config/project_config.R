# ============================================================
# uenoy scRNAseq Framework
# Project configuration
# Version 4.0.0
# ============================================================

project_config <- list(
  framework_name = "uenoy scRNAseq Framework",
  framework_version = "4.0.0",
  script_project_root = Sys.getenv(
    "UENOY_SCRNA_ROOT",
    unset = "/Users/uenoya/Projects/uenoyscRNA"
  ),
  external_data_root = Sys.getenv(
    "UENOY_SCRNA_DATA_ROOT",
    unset = "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
  )
)

# Dataset identity is no longer tied to a directory name. Each entry may define
# one or more preferred files plus fallback search patterns. Preferred files are
# intentionally based on the current validated project structure.
.dataset_root <- project_config$external_data_root

dataset_registry <- list(
  Mouse_MASH_RDS = list(
    label = "Mouse MASH whole-cell",
    preferred_files = file.path(.dataset_root, "Mouse_MASH_RDS", c(
      "Mouse_RH260519343_GSE325222_RPCA_integrated_celltype_annotated.rds",
      "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
    )),
    include_patterns = c("mouse|gse325222", "celltype|annotated|annotation"),
    exclude_patterns = c("mphi|macrophage|bundle|intermediate|raw"),
    selection = "latest",
    result_root = file.path(.dataset_root, "Mouse_MASH_RDS", "analysis_results")
  ),
  Mouse_MASH_Mphi_RDS = list(
    label = "Mouse MASH macrophage",
    preferred_files = file.path(
      .dataset_root, "results",
      "Mouse_Mphi_RPCA_two_layer_annotation_layer2_by_layer1_color",
      "Mouse_Mphi_RPCA_two_layer_annotation_layer2_by_layer1_color.rds"
    ),
    include_patterns = c("mouse.*(?:mphi|macrophage)|(?:mphi|macrophage).*mouse"),
    exclude_patterns = c("bundle|intermediate|raw|without_object"),
    selection = "latest",
    result_root = file.path(.dataset_root, "analysis_results", "Mouse_MASH_Mphi_RDS")
  ),
  Human_MASH_RDS = list(
    label = "Human MASLD/MASH whole-cell",
    preferred_files = file.path(
      .dataset_root, "rds",
      "Human_LT1to8_SfLB1_SfLB2_Tx17_Tx5_RPCA_integrated_12samples_celltype_annotated.rds"
    ),
    include_patterns = c("human|lt1to8|12samples", "celltype|annotated"),
    exclude_patterns = c("mphi|macrophage|mesenchymal|bundle|raw"),
    selection = "latest",
    result_root = file.path(.dataset_root, "analysis_results", "Human_MASH_RDS")
  ),
  Human_MASH_Mphi_RDS = list(
    label = "Human MASLD/MASH macrophage",
    preferred_files = file.path(
      .dataset_root, "results",
      "Human_liver_macrophage_reclustering_strict_M2_definition_corrected",
      "Human_macrophage_reclustered_with_strict_M2_annotations.rds"
    ),
    include_patterns = c("human", "macrophage|mphi"),
    exclude_patterns = c("bundle|intermediate|raw|without_object"),
    selection = "latest",
    result_root = file.path(.dataset_root, "analysis_results", "Human_MASH_Mphi_RDS")
  )
)

allowed_datasets <- names(dataset_registry)

get_dataset_spec <- function(dataset_name) {
  if (!dataset_name %in% allowed_datasets) {
    stop(
      "Unknown dataset_name: ", dataset_name,
      "\nAllowed values: ", paste(allowed_datasets, collapse = ", "),
      call. = FALSE
    )
  }
  dataset_registry[[dataset_name]]
}

# Backward-compatible directory accessor. It returns the common directory of
# resolved preferred files when possible, otherwise the external data root.
get_dataset_dir <- function(dataset_name) {
  spec <- get_dataset_spec(dataset_name)
  preferred <- spec$preferred_files
  existing <- preferred[file.exists(preferred)]
  if (length(existing)) return(dirname(existing[[1L]]))
  project_config$external_data_root
}

get_dataset_files <- function(dataset_name, catalog = NULL, strict = FALSE) {
  if (!exists("resolve_dataset_rds", mode = "function")) {
    stop("R/rds_registry.R must be sourced before get_dataset_files().", call. = FALSE)
  }
  resolve_dataset_rds(
    dataset_name = dataset_name,
    registry = dataset_registry,
    external_data_root = project_config$external_data_root,
    catalog = catalog,
    strict = strict
  )
}

get_result_root <- function(dataset_name) {
  spec <- get_dataset_spec(dataset_name)
  spec$result_root
}

get_result_dir <- function(dataset_name, result_type = NULL, create = TRUE) {
  root <- get_result_root(dataset_name)
  result_map <- c(
    inventory = "00_inventory", umap = "01_UMAP",
    composition = "02_cell_composition", cluster_highlight = "03_cluster_highlight",
    marker_dynamics = "04_marker_dynamics", module_score = "05_module_score",
    de = "06_DE_analysis", tables = "07_tables", logs = "08_logs"
  )
  out <- if (is.null(result_type)) root else {
    if (!result_type %in% names(result_map)) stop("Unknown result_type: ", result_type)
    file.path(root, unname(result_map[[result_type]]))
  }
  if (isTRUE(create)) dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}
