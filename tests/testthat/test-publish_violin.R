test_that("publish_violin returns a patchwork object", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  features <- head(
    rownames(object),
    2
  )

  plot <- suppressWarnings(
    publish_violin(
      object = object,
      features = features,
      group.by = "test_group",
      point_size = 0,
      raster = FALSE,
      combine = TRUE
    )
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("publish_violin returns a list when combine is FALSE", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  features <- head(
    rownames(object),
    2
  )

  plots <- suppressWarnings(
    publish_violin(
      object = object,
      features = features,
      group.by = "test_group",
      point_size = 0,
      raster = FALSE,
      combine = FALSE
    )
  )

  expect_type(
    plots,
    "list"
  )

  expect_length(
    plots,
    length(features)
  )

  expect_true(
    all(
      vapply(
        plots,
        inherits,
        logical(1),
        what = "ggplot"
      )
    )
  )
})


test_that("publish_violin accepts group ordering", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  feature <- rownames(object)[1]

  plots <- suppressWarnings(
    publish_violin(
      object = object,
      features = feature,
      group.by = "test_group",
      group_order = c("B", "A"),
      point_size = 0,
      raster = FALSE,
      combine = FALSE
    )
  )

  observed_levels <- levels(
    plots[[1]]$data$ident
  )

  expect_equal(
    observed_levels,
    c("B", "A")
  )
})


test_that("publish_violin accepts split ordering", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  object$test_split <- rep(
    c("X", "Y"),
    length.out = ncol(object)
  )

  feature <- rownames(object)[1]

  expect_no_error(
    suppressWarnings(
      publish_violin(
        object = object,
        features = feature,
        group.by = "test_group",
        split.by = "test_split",
        split_order = c("Y", "X"),
        point_size = 0,
        raster = FALSE
      )
    )
  )
})


test_that("publish_violin validates the Seurat object", {
  expect_error(
    publish_violin(
      object = data.frame(),
      features = "GeneA"
    ),
    "`object` must be a Seurat object",
    fixed = TRUE
  )
})


test_that("publish_violin validates features", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = character()
    ),
    "`features` must be a non-empty character vector",
    fixed = TRUE
  )

  expect_error(
    publish_violin(
      object = object,
      features = NA_character_
    ),
    "`features` must be a non-empty character vector",
    fixed = TRUE
  )
})


test_that("publish_violin validates assay names", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      assay = "not_an_assay"
    ),
    "was not found"
  )
})


test_that("publish_violin validates group.by", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      group.by = "not_a_metadata_column"
    ),
    "was not found"
  )
})


test_that("publish_violin validates split.by", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      split.by = "not_a_metadata_column"
    ),
    "was not found"
  )
})


test_that("publish_violin validates group_order values", {
  object <- SeuratObject::pbmc_small

  object$test_group <- rep(
    c("A", "B"),
    length.out = ncol(object)
  )

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      group.by = "test_group",
      group_order = c("A", "MissingGroup"),
      point_size = 0
    ),
    "MissingGroup"
  )
})


test_that("publish_violin requires split.by with split_order", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      split_order = c("X", "Y"),
      point_size = 0
    ),
    "`split_order` can only be used when `split.by` is specified",
    fixed = TRUE
  )
})


test_that("publish_violin validates logical arguments", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      combine = "yes"
    ),
    "`combine` must be TRUE or FALSE",
    fixed = TRUE
  )
})


test_that("publish_violin validates point_size", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      point_size = -1
    ),
    "`point_size` must be a single non-negative number",
    fixed = TRUE
  )
})


test_that("publish_violin validates fill_by", {
  object <- SeuratObject::pbmc_small

  expect_error(
    publish_violin(
      object = object,
      features = rownames(object)[1],
      fill_by = "invalid"
    ),
    '`fill_by` must be either "feature" or "ident"',
    fixed = TRUE
  )
})
