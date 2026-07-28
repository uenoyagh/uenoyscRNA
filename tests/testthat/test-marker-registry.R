test_that("valid marker registry passes validation", {
  registry <- data.frame(
    species = c("Mouse", "Mouse"),
    tissue = c("Liver", "Liver"),
    layer = c("Layer2", "Layer2"),
    celltype = c("Resident Kupffer-like", "Resident Kupffer-like"),
    marker_group = c("identity", "identity"),
    gene = c("Clec4f", "Timd4"),
    direction = c("positive", "positive"),
    evidence = c("curated", "curated"),
    notes = c("", ""),
    stringsAsFactors = FALSE
  )

  expect_true(validate_marker_registry(registry))
})

test_that("missing required columns fail validation", {
  registry <- data.frame(species = "Mouse", gene = "Clec4f")
  expect_error(validate_marker_registry(registry), "missing required columns")
})

test_that("invalid directions fail validation", {
  registry <- data.frame(
    species = "Mouse",
    tissue = "Liver",
    layer = "Layer2",
    celltype = "Resident Kupffer-like",
    marker_group = "identity",
    gene = "Clec4f",
    direction = "up",
    evidence = "curated",
    notes = "",
    stringsAsFactors = FALSE
  )
  expect_error(validate_marker_registry(registry), "Invalid `direction`")
})

test_that("duplicate marker definitions fail validation", {
  registry <- data.frame(
    species = rep("Mouse", 2),
    tissue = rep("Liver", 2),
    layer = rep("Layer2", 2),
    celltype = rep("Resident Kupffer-like", 2),
    marker_group = rep("identity", 2),
    gene = rep("Clec4f", 2),
    direction = rep("positive", 2),
    evidence = rep("curated", 2),
    notes = rep("", 2),
    stringsAsFactors = FALSE
  )
  expect_error(validate_marker_registry(registry), "duplicate")
})

test_that("get_markers returns filtered unique genes", {
  registry <- data.frame(
    species = c("Mouse", "Mouse", "Human"),
    tissue = c("Liver", "Liver", "Liver"),
    layer = c("Layer2", "Layer2", "Layer2"),
    celltype = c("Resident Kupffer-like", "Monocyte-like", "Resident Kupffer-like"),
    marker_group = c("identity", "identity", "identity"),
    gene = c("Clec4f", "Ccr2", "VSIG4"),
    direction = c("positive", "positive", "positive"),
    evidence = c("curated", "curated", "curated"),
    notes = c("", "", ""),
    stringsAsFactors = FALSE
  )

  result <- get_markers(
    registry,
    species = "Mouse",
    celltype = "Resident Kupffer-like"
  )
  expect_identical(result, "Clec4f")
})

test_that("get_markers can return registry rows", {
  registry <- data.frame(
    species = "Mouse",
    tissue = "Liver",
    layer = "Layer2",
    celltype = "Resident Kupffer-like",
    marker_group = "identity",
    gene = "Clec4f",
    direction = "positive",
    evidence = "curated",
    notes = "",
    stringsAsFactors = FALSE
  )

  result <- get_markers(registry, unique_genes = FALSE)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
})

test_that("CSV round trip preserves character columns", {
  registry <- data.frame(
    species = "Mouse",
    tissue = "Liver",
    layer = "2",
    celltype = "Resident Kupffer-like",
    marker_group = "identity",
    gene = "Clec4f",
    direction = "positive",
    evidence = "curated",
    notes = "",
    stringsAsFactors = FALSE
  )

  path <- tempfile(fileext = ".csv")
  write_marker_registry(registry, path)
  actual <- read_marker_registry(path)

  expect_type(actual$layer, "character")
  expect_identical(actual$gene, "Clec4f")
})

test_that("list_marker_sets returns distinct definitions", {
  registry <- data.frame(
    species = c("Mouse", "Mouse"),
    tissue = c("Liver", "Liver"),
    layer = c("Layer2", "Layer2"),
    celltype = c("Resident Kupffer-like", "Resident Kupffer-like"),
    marker_group = c("identity", "identity"),
    gene = c("Clec4f", "Timd4"),
    direction = c("positive", "positive"),
    evidence = c("curated", "curated"),
    notes = c("", ""),
    stringsAsFactors = FALSE
  )

  result <- list_marker_sets(registry)
  expect_equal(nrow(result), 1L)
})
