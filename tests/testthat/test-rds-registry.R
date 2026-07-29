test_that("discover_rds_files excludes bundles and renv", {
  root <- tempfile("rds_registry_")
  dir.create(root)
  dir.create(file.path(root, "rds"))
  dir.create(file.path(root, "renv", "library"), recursive = TRUE)
  saveRDS(list(a = 1), file.path(root, "rds", "Mouse_annotated.rds"))
  saveRDS(list(a = 1), file.path(root, "rds", "analysis_bundle.rds"))
  saveRDS(list(a = 1), file.path(root, "renv", "library", "package.rds"))
  out <- discover_rds_files(root)
  expect_equal(out$file, "Mouse_annotated.rds")
})

test_that("resolve_dataset_rds prefers explicit existing files", {
  root <- tempfile("rds_registry_")
  dir.create(root)
  path <- file.path(root, "preferred.rds")
  saveRDS(1, path)
  registry <- list(test = list(preferred_files = path, selection = "latest"))
  expect_equal(resolve_dataset_rds("test", registry, root), normalizePath(path))
})

test_that("resolve_dataset_rds uses fallback patterns", {
  root <- tempfile("rds_registry_")
  dir.create(root)
  path <- file.path(root, "Human_macrophage_annotated.rds")
  saveRDS(1, path)
  registry <- list(test = list(
    preferred_files = file.path(root, "missing.rds"),
    include_patterns = c("human", "macrophage"),
    exclude_patterns = "bundle",
    selection = "latest"
  ))
  expect_equal(resolve_dataset_rds("test", registry, root), normalizePath(path))
})
