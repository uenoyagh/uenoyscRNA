test_that("review_marker_dotplot validates object class", {
  expect_error(
    review_marker_dotplot(
      object = list(),
      markers = c("Gene1", "Gene2")
    ),
    "must inherit"
  )
})

test_that("review_marker_dotplot validates quantile cutoff", {
  object <- make_review_test_object()

  expect_error(
    review_marker_dotplot(
      object = object,
      markers = c("Gene1", "Gene2"),
      quantile_cutoff = c(0.9, 0.1)
    ),
    "quantile_cutoff"
  )
})

test_that("review_marker_dotplot creates plots", {
  object <- make_review_test_object()

  result <- review_marker_dotplot(
    object = object,
    markers = list(
      Group_A = c("Gene1", "Gene2"),
      Group_B = c("Gene3", "MissingGene")
    ),
    annotation_column = "annotation",
    assay = "RNA",
    max_features_per_plot = 3
  )

  expect_s3_class(result, "uenoy_review_dotplot")
  expect_true(length(result$plots) >= 1L)
  expect_true("MissingGene" %in% result$missing_features)
  expect_true(all(vapply(result$plots, inherits, logical(1), "ggplot")))
})

test_that("review_marker_dotplot writes PDF files", {
  object <- make_review_test_object()
  output_dir <- tempfile("dotplot-review-")

  result <- review_marker_dotplot(
    object = object,
    markers = list(
      Group_A = c("Gene1", "Gene2"),
      Group_B = "Gene3"
    ),
    annotation_column = "annotation",
    assay = "RNA",
    output_dir = output_dir
  )

  expect_true(length(result$files) >= 1L)
  expect_true(all(file.exists(result$files)))
})

test_that("marker registry columns are detected flexibly", {
  object <- make_review_test_object()

  registry <- data.frame(
    species = rep("mouse", 3),
    tissue = rep("liver", 3),
    layer = rep("layer1", 3),
    cell_type = c("Type_A", "Type_A", "Type_B"),
    gene_symbol = c("Gene1", "Gene2", "Gene3"),
    stringsAsFactors = FALSE
  )

  result <- review_marker_dotplot(
    object = object,
    marker_registry = registry,
    annotation_column = "annotation",
    assay = "RNA",
    species = "mouse",
    tissue = "liver",
    layer = "layer1"
  )

  expect_named(result$marker_groups, c("Type_A", "Type_B"))
  expect_equal(
    sort(unique(unlist(result$marker_groups))),
    c("Gene1", "Gene2", "Gene3")
  )
})
