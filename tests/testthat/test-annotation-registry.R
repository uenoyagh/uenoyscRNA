test_that("valid registry passes", {
  x <- data.frame(
    species = "Mouse", tissue = "Liver", dataset = "Example",
    layer = "Layer2", cluster = "1",
    annotation = "Resident Kupffer-like", confidence = "High",
    markers = "Clec4f;Timd4;Vsig4", evidence = "Marker review",
    reviewer = "YU", review_date = "2026-07-29",
    stringsAsFactors = FALSE
  )
  expect_invisible(validate_annotation_registry(x))
})

test_that("missing columns fail", {
  expect_error(
    validate_annotation_registry(data.frame(species = "Mouse")),
    "missing required column"
  )
})

test_that("duplicate keys fail", {
  x <- data.frame(
    species = c("Mouse", "Mouse"), tissue = c("Liver", "Liver"),
    dataset = c("Example", "Example"), layer = c("Layer2", "Layer2"),
    cluster = c("1", "1"), annotation = c("A", "B"),
    confidence = c("High", "Low"), markers = c("", ""),
    evidence = c("", ""), reviewer = c("", ""),
    review_date = c("", ""), stringsAsFactors = FALSE
  )
  expect_error(validate_annotation_registry(x), "Duplicate annotation key")
})

test_that("invalid confidence fails", {
  x <- data.frame(
    species = "Mouse", tissue = "Liver", dataset = "Example",
    layer = "Layer2", cluster = "1", annotation = "A",
    confidence = "Certain", markers = "", evidence = "",
    reviewer = "", review_date = "", stringsAsFactors = FALSE
  )
  expect_error(validate_annotation_registry(x), "Invalid confidence")
})

test_that("lookup returns annotations", {
  x <- data.frame(
    species = c("Mouse", "Mouse"), tissue = c("Liver", "Liver"),
    dataset = c("Example", "Example"), layer = c("Layer2", "Layer2"),
    cluster = c("1", "2"),
    annotation = c("Resident Kupffer-like", "Monocyte-like"),
    confidence = c("High", "Medium"), markers = c("", ""),
    evidence = c("", ""), reviewer = c("", ""),
    review_date = c("", ""), stringsAsFactors = FALSE
  )

  observed <- get_cluster_annotation(
    x, c("2", "1"), "Mouse", "Liver", "Example", "Layer2",
    return = "annotation"
  )

  expect_equal(
    observed,
    c("2" = "Monocyte-like", "1" = "Resident Kupffer-like")
  )
})

test_that("CSV round trip works", {
  x <- data.frame(
    species = "Mouse", tissue = "Liver", dataset = "Example",
    layer = "Layer2", cluster = "1",
    annotation = "Resident Kupffer-like", confidence = "High",
    markers = "Clec4f;Timd4;Vsig4", evidence = "Marker review",
    reviewer = "YU", review_date = "2026-07-29",
    stringsAsFactors = FALSE
  )

  path <- tempfile(fileext = ".csv")
  write_annotation_registry(x, path)
  expect_equal(read_annotation_registry(path), x)
  expect_error(write_annotation_registry(x, path), "already exists")
})
