test_that("create_review_manifest returns planned outputs", {
  manifest <- create_review_manifest(
    include = c("qc", "umap"),
    output_dir = "review"
  )

  expect_s3_class(manifest, "data.frame")
  expect_equal(manifest$component, c("qc", "umap"))
  expect_equal(manifest$status, c("planned", "planned"))
})

test_that("create_review_manifest rejects unsupported components", {
  expect_error(
    create_review_manifest(include = "heatmap"),
    "Unsupported review component"
  )
})

test_that("review functions reject non-Seurat objects", {
  expect_error(
    check_review_metadata(data.frame(), "annotation"),
    "must inherit from class `Seurat`"
  )

  expect_error(
    resolve_marker_features(data.frame(), "Clec4f"),
    "must inherit from class `Seurat`"
  )

  expect_error(
    review_seurat_object(data.frame(), "annotation"),
    "must inherit from class `Seurat`"
  )
})

test_that("review core validates a minimal Seurat object", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    sparse = TRUE,
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$annotation <- c("Resident Kupffer-like", "Monocyte-like")
  object$sample <- c("S1", "S1")
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      dimnames = list(c("cell1", "cell2"), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  result <- review_seurat_object(
    object = object,
    annotation_column = "annotation",
    sample_column = "sample"
  )

  expect_s3_class(result, "uenoy_review")
  expect_equal(result$summary$n_cells, 2L)
  expect_equal(result$summary$reduction, "umap")
})

test_that("missing metadata columns are reported", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    sparse = TRUE,
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    check_review_metadata(object, "annotation"),
    "Missing metadata column"
  )
})

test_that("resolve_marker_features reports coverage", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    sparse = TRUE,
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)

  result <- resolve_marker_features(
    object,
    c("Clec4f", "Timd4", "Ccr2")
  )

  expect_equal(result$present, c("Clec4f", "Timd4"))
  expect_equal(result$missing, "Ccr2")
  expect_equal(result$coverage, 2 / 3)
})
