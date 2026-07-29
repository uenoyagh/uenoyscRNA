root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(root, "R", "rds_registry.R"))
source(file.path(root, "config", "project_config.R"))
stopifnot(identical(project_config$framework_version, "4.0.0"))
stopifnot(setequal(names(dataset_registry), c(
  "Mouse_MASH_RDS", "Mouse_MASH_Mphi_RDS", "Human_MASH_RDS", "Human_MASH_Mphi_RDS"
)))
cat("Static v4 registry validation passed.\n")
