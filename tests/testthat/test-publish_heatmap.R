test_that("publish_heatmap returns a plot object", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    3
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  plot <- publish_heatmap(
    object = object,
    features = features,
    assay = "RNA",
    raster = FALSE
  )

  expect_true(
    inherits(plot, "ggplot") ||
      inherits(plot, "patchwork")
  )
})


test_that("publish_heatmap accepts metadata grouping", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  features <- head(
    rownames(object[["RNA"]]),
    3
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_no_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      group.by = "test_group",
      group_order = c("B", "A"),
      raster = FALSE
    )
  )
})


test_that("publish_heatmap reverses feature order", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    3
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_no_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      reverse_features = TRUE,
      raster = FALSE
    )
  )
})


test_that("publish_heatmap validates the Seurat object", {
  expect_error(
    publish_heatmap(
      object = data.frame(),
      features = "GeneA"
    ),
    "`object` must be a Seurat object.",
    fixed = TRUE
  )
})


test_that("publish_heatmap validates features", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_heatmap(
      object = object,
      features = character()
    ),
    "`features` must be a non-empty character vector",
    fixed = TRUE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = NA_character_
    ),
    "`features` must be a non-empty character vector",
    fixed = TRUE
  )
})


test_that("publish_heatmap detects missing features", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_heatmap(
      object = object,
      features = "not_a_real_gene",
      assay = "RNA"
    ),
    "not_a_real_gene"
  )
})


test_that("publish_heatmap validates assay names", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_heatmap(
      object = object,
      features = rownames(object)[1],
      assay = "not_an_assay"
    ),
    "was not found"
  )
})


test_that("publish_heatmap validates slot", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_heatmap(
      object = object,
      features = rownames(object)[1],
      slot = "invalid_slot"
    ),
    "`slot` must be one of",
    fixed = TRUE
  )
})


test_that("publish_heatmap validates group.by", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    2
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      group.by = "not_a_metadata_column"
    ),
    "not_a_metadata_column"
  )
})


test_that("publish_heatmap validates group_order", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  features <- head(
    rownames(object[["RNA"]]),
    2
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      group.by = "test_group",
      group_order = c("A", "MissingGroup")
    ),
    "MissingGroup"
  )
})


test_that("publish_heatmap validates cells", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    2
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      cells = "not_a_real_cell"
    ),
    "not_a_real_cell"
  )
})


test_that("publish_heatmap validates display limits", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    2
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      disp.min = 2,
      disp.max = -2
    ),
    "`disp.max` must be greater than `disp.min`.",
    fixed = TRUE
  )
})


test_that("publish_heatmap validates logical arguments", {
  object <- SeuratObject::pbmc_small

  features <- head(
    rownames(object[["RNA"]]),
    2
  )

  object <- Seurat::ScaleData(
    object,
    features = features,
    verbose = FALSE
  )

  expect_error(
    publish_heatmap(
      object = object,
      features = features,
      assay = "RNA",
      raster = "yes"
    ),
    "`raster` must be TRUE or FALSE.",
    fixed = TRUE
  )
})
