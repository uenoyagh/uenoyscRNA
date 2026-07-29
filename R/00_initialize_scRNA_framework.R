# ============================================================
# uenoy scRNAseq Framework
# Initial directory setup
#
# Data and analysis outputs:
#   External SSD
#
# R scripts and reusable functions:
#   Mac internal SSD
# ============================================================


# ============================================================
# 0. User settings
# ============================================================

# 外付けSSD上のデータルート
external_data_root <- file.path(
  "/Volumes",
  "SSD990_uenoy",
  "scRNA_MASLD_MASH_pseudobulk"
)

# Mac内蔵SSD上のスクリプト管理場所
#
# 必要に応じて、この1行だけ変更してください。
script_project_root <- path.expand(
  "~/R/uenoy_scRNAseq_Framework"
)


# ============================================================
# 1. Check external SSD
# ============================================================

external_ssd_mount <- file.path(
  "/Volumes",
  "SSD990_uenoy"
)

if (!dir.exists(external_ssd_mount)) {
  stop(
    paste0(
      "\nExternal SSD is not mounted.\n",
      "Expected location:\n",
      external_ssd_mount,
      "\n"
    )
  )
}

if (!dir.exists(external_data_root)) {
  dir.create(
    external_data_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(external_data_root)) {
  stop(
    paste0(
      "\nFailed to create external data root:\n",
      external_data_root,
      "\n"
    )
  )
}


# ============================================================
# 2. External SSD directory structure
# ============================================================

external_dataset_names <- c(
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_RDS",
  "Human_MASH_RDS",
  "Human_MASH_Mphi_RDS"
)

external_dataset_dirs <- file.path(
  external_data_root,
  external_dataset_names
)

for (dataset_dir in external_dataset_dirs) {

  dir.create(
    dataset_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  result_dirs <- file.path(
    dataset_dir,
    c(
      "comparison_results",
      "comparison_results/00_inventory",
      "comparison_results/01_UMAP",
      "comparison_results/02_cell_composition",
      "comparison_results/03_cluster_highlight",
      "comparison_results/04_marker_dynamics",
      "comparison_results/05_module_score",
      "comparison_results/06_DE_analysis",
      "comparison_results/07_tables",
      "comparison_results/08_logs"
    )
  )

  for (result_dir in result_dirs) {
    dir.create(
      result_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}


# ============================================================
# 3. Mac internal SSD directory structure
# ============================================================

internal_dirs <- file.path(
  script_project_root,
  c(
    "analysis",
    "R",
    "config",
    "logs",
    "tests"
  )
)

dir.create(
  script_project_root,
  recursive = TRUE,
  showWarnings = FALSE
)

for (internal_dir in internal_dirs) {
  dir.create(
    internal_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 4. Define configuration paths
# ============================================================

config_dir <- file.path(
  script_project_root,
  "config"
)

config_file <- file.path(
  config_dir,
  "project_config.R"
)

local_config_file <- file.path(
  config_dir,
  "local_config.R"
)


# ============================================================
# 5. Create project_config.R
# ============================================================

project_config_lines <- c(
  "# ============================================================",
  "# uenoy scRNAseq Framework",
  "# Project configuration",
  "# ============================================================",
  "",
  "project_config <- list(",
  "",
  "  framework_name = \"uenoy scRNAseq Framework\",",
  "",
  "  framework_version = \"2.0\",",
  "",
  paste0(
    "  script_project_root = ",
    deparse(script_project_root),
    ","
  ),
  "",
  paste0(
    "  external_data_root = ",
    deparse(external_data_root),
    ","
  ),
  "",
  paste0(
    "  Mouse_MASH_RDS = ",
    deparse(file.path(external_data_root, "Mouse_MASH_RDS")),
    ","
  ),
  "",
  paste0(
    "  Mouse_MASH_Mphi_RDS = ",
    deparse(file.path(external_data_root, "Mouse_MASH_Mphi_RDS")),
    ","
  ),
  "",
  paste0(
    "  Human_MASH_RDS = ",
    deparse(file.path(external_data_root, "Human_MASH_RDS")),
    ","
  ),
  "",
  paste0(
    "  Human_MASH_Mphi_RDS = ",
    deparse(file.path(external_data_root, "Human_MASH_Mphi_RDS")),
    ","
  ),
  "",
  paste0(
    "  analysis_dir = ",
    deparse(file.path(script_project_root, "analysis")),
    ","
  ),
  "",
  paste0(
    "  functions_dir = ",
    deparse(file.path(script_project_root, "R")),
    ","
  ),
  "",
  paste0(
    "  config_dir = ",
    deparse(file.path(script_project_root, "config")),
    ","
  ),
  "",
  paste0(
    "  log_dir = ",
    deparse(file.path(script_project_root, "logs")),
    ""
  ),
  ")",
  "",
  "",
  "# ============================================================",
  "# Dataset helper",
  "# ============================================================",
  "",
  "get_dataset_dir <- function(dataset_name) {",
  "",
  "  allowed_datasets <- c(",
  "    \"Mouse_MASH_RDS\",",
  "    \"Mouse_MASH_Mphi_RDS\",",
  "    \"Human_MASH_RDS\",",
  "    \"Human_MASH_Mphi_RDS\"",
  "  )",
  "",
  "  if (!dataset_name %in% allowed_datasets) {",
  "    stop(",
  "      \"Unknown dataset_name: \",",
  "      dataset_name,",
  "      \"\\nAllowed values: \",",
  "      paste(allowed_datasets, collapse = \", \")",
  "    )",
  "  }",
  "",
  "  project_config[[dataset_name]]",
  "}",
  "",
  "",
  "# ============================================================",
  "# Result directory helper",
  "# ============================================================",
  "",
  "get_result_dir <- function(dataset_name, result_type = NULL) {",
  "",
  "  dataset_dir <- get_dataset_dir(dataset_name)",
  "",
  "  result_root <- file.path(",
  "    dataset_dir,",
  "    \"comparison_results\"",
  "  )",
  "",
  "  if (is.null(result_type)) {",
  "    return(result_root)",
  "  }",
  "",
  "  result_map <- c(",
  "    inventory = \"00_inventory\",",
  "    umap = \"01_UMAP\",",
  "    composition = \"02_cell_composition\",",
  "    cluster_highlight = \"03_cluster_highlight\",",
  "    marker_dynamics = \"04_marker_dynamics\",",
  "    module_score = \"05_module_score\",",
  "    de = \"06_DE_analysis\",",
  "    tables = \"07_tables\",",
  "    logs = \"08_logs\"",
  "  )",
  "",
  "  if (!result_type %in% names(result_map)) {",
  "    stop(",
  "      \"Unknown result_type: \",",
  "      result_type,",
  "      \"\\nAllowed values: \",",
  "      paste(names(result_map), collapse = \", \")",
  "    )",
  "  }",
  "",
  "  file.path(",
  "    result_root,",
  "    unname(result_map[[result_type]])",
  "  )",
  "}"
)

writeLines(
  project_config_lines,
  con = config_file
)


# ============================================================
# 6. Create local_config.R
# ============================================================

local_config_lines <- c(
  "# ============================================================",
  "# Local execution settings",
  "# ============================================================",
  "",
  "# 解析対象を以下の4種類から指定します。",
  "#",
  "# Mouse_MASH_RDS",
  "# Mouse_MASH_Mphi_RDS",
  "# Human_MASH_RDS",
  "# Human_MASH_Mphi_RDS",
  "",
  "analysis_target <- \"Mouse_MASH_RDS\"",
  "",
  "# FALSEの場合、既存ファイルを上書きしません。",
  "overwrite_existing <- FALSE",
  "",
  "# UMAP point size",
  "umap_point_size <- 0.5",
  "",
  "# UMAP rasterization",
  "umap_raster <- TRUE",
  "",
  "# Default PDF dimensions",
  "pdf_width <- 14",
  "pdf_height <- 10"
)

writeLines(
  local_config_lines,
  con = local_config_file
)


# ============================================================
# 7. Create starter analysis scripts
# ============================================================

analysis_dir <- file.path(
  script_project_root,
  "analysis"
)

starter_scripts <- c(
  "00_project_setup.R",
  "01_inventory_RDS.R",
  "02_publish_UMAP.R",
  "03_cell_composition.R",
  "04_cluster_highlight.R",
  "05_marker_dynamics.R",
  "06_module_score.R",
  "07_DE_analysis.R",
  "08_batch_pipeline.R"
)

for (script_name in starter_scripts) {

  script_path <- file.path(
    analysis_dir,
    script_name
  )

  if (!file.exists(script_path)) {

    script_header <- c(
      "# ============================================================",
      paste0("# ", script_name),
      "# uenoy scRNAseq Framework",
      "# ============================================================",
      "",
      "rm(list = ls())",
      "",
      paste0(
        "source(",
        deparse(config_file),
        ")"
      ),
      "",
      paste0(
        "source(",
        deparse(local_config_file),
        ")"
      ),
      "",
      "target_dir <- get_dataset_dir(analysis_target)",
      "",
      "output_dir <- get_result_dir(analysis_target)",
      "",
      "cat(\"Analysis target:\\n\", analysis_target, \"\\n\")",
      "cat(\"Input directory:\\n\", target_dir, \"\\n\")",
      "cat(\"Output directory:\\n\", output_dir, \"\\n\")",
      ""
    )

    writeLines(
      script_header,
      con = script_path
    )
  }
}


# ============================================================
# 8. Create starter reusable function files
# ============================================================

function_dir <- file.path(
  script_project_root,
  "R"
)

function_files <- c(
  "io.R",
  "metadata.R",
  "seurat_helpers.R",
  "umap.R",
  "plotting.R",
  "marker_dynamics.R",
  "module_score.R",
  "utils.R"
)

for (function_name in function_files) {

  function_path <- file.path(
    function_dir,
    function_name
  )

  if (!file.exists(function_path)) {

    function_header <- c(
      "# ============================================================",
      paste0("# ", function_name),
      "# Reusable functions",
      "# ============================================================",
      ""
    )

    writeLines(
      function_header,
      con = function_path
    )
  }
}


# ============================================================
# 9. Create README
# ============================================================

readme_file <- file.path(
  script_project_root,
  "README.md"
)

if (!file.exists(readme_file)) {

  readme_lines <- c(
    "# uenoy scRNAseq Framework",
    "",
    "## Data storage",
    "",
    "RDS files and analysis outputs are stored on the external SSD:",
    "",
    paste0("`", external_data_root, "`"),
    "",
    "## Script storage",
    "",
    "R scripts and reusable functions are stored on the Mac internal drive:",
    "",
    paste0("`", script_project_root, "`"),
    "",
    "## Dataset directories",
    "",
    "- Mouse_MASH_RDS",
    "- Mouse_MASH_Mphi_RDS",
    "- Human_MASH_RDS",
    "- Human_MASH_Mphi_RDS",
    "",
    "## Main configuration",
    "",
    "- config/project_config.R",
    "- config/local_config.R"
  )

  writeLines(
    readme_lines,
    con = readme_file
  )
}


# ============================================================
# 10. Final validation
# ============================================================

required_external_dirs <- c(
  external_data_root,
  external_dataset_dirs
)

required_internal_dirs <- c(
  script_project_root,
  internal_dirs
)

missing_external_dirs <- required_external_dirs[
  !dir.exists(required_external_dirs)
]

missing_internal_dirs <- required_internal_dirs[
  !dir.exists(required_internal_dirs)
]

if (length(missing_external_dirs) > 0) {
  stop(
    paste0(
      "\nMissing external directories:\n",
      paste(missing_external_dirs, collapse = "\n")
    )
  )
}

if (length(missing_internal_dirs) > 0) {
  stop(
    paste0(
      "\nMissing internal directories:\n",
      paste(missing_internal_dirs, collapse = "\n")
    )
  )
}


# ============================================================
# 11. Summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("uenoy scRNAseq Framework setup completed\n")
cat("============================================================\n\n")

cat("External SSD data root:\n")
cat(external_data_root, "\n\n")

cat("Mac internal script root:\n")
cat(script_project_root, "\n\n")

cat("Dataset directories:\n")
cat(
  paste0(
    "- ",
    external_dataset_dirs
  ),
  sep = "\n"
)

cat("\n\nConfiguration files:\n")
cat("- ", config_file, "\n", sep = "")
cat("- ", local_config_file, "\n", sep = "")

cat("\nAnalysis scripts:\n")
cat(
  paste0(
    "- ",
    file.path(
      analysis_dir,
      starter_scripts
    )
  ),
  sep = "\n"
)

cat("\n\nSetup completed successfully.\n")
