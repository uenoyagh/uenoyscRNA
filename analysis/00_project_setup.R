rm(list=ls())

script_root <- "/Users/uenoya/Projects/uenoyscRNA"
source(file.path(script_root, "config", "project_config.R"))

if (!dir.exists(project_config$external_data_root)) {
  stop("External SSD data root was not found: ", project_config$external_data_root)
}

for (d in file.path(script_root, c("analysis","R","config","logs","tests"))) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

for (nm in allowed_datasets) {
  dir.create(get_dataset_dir(nm), recursive=TRUE, showWarnings=FALSE)
  for (type in c("inventory","umap","composition","cluster_highlight",
                 "marker_dynamics","module_score","de","tables","logs")) {
    get_result_dir(nm, type, create=TRUE)
  }
}

cat("Framework directory validation completed.\n")
