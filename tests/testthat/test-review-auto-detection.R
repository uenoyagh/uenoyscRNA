test_that("detect_review_metadata prefers FIXED2 review columns", {
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
  object$sample <- c("raw1", "raw2")
  object$sample_for_R8plot_FIXED2 <- c("S1", "S2")
  object$sample_display_FIXED2 <- c("STD", "Tx")
  object$condition <- c("raw", "raw")
  object$condition_FIXED2 <- c("STD", "Tx")
  object$celltype_auto_annotation <- c("Macrophage", "Monocyte")
  object$celltype_for_R8plot_FIXED2 <- c(
    "Resident Kupffer-like",
    "Monocyte-like"
  )
  object$celltype_annotation_confidence <- c("high", "high")

  detected <- detect_review_metadata(object)

  expect_equal(
    detected$annotation_column,
    "celltype_for_R8plot_FIXED2"
  )
  expect_equal(
    detected$sample_column,
    "sample_display_FIXED2"
  )
  expect_equal(
    detected$condition_column,
    "condition_FIXED2"
  )
})

test_that("detect_review_reduction finds custom UMAP names", {
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
  object[["umapRPCA"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      dimnames = list(c("cell1", "cell2"), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  detected <- detect_review_reduction(object)

  expect_equal(detected$reduction, "umapRPCA")
})

test_that("detect_review_settings respects explicit overrides", {
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
  object$annotation_custom <- c("A", "B")
  object$sample_custom <- c("S1", "S2")
  object$condition_custom <- c("STD", "Tx")
  object[["umap_custom"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      dimnames = list(c("cell1", "cell2"), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  detected <- detect_review_settings(
    object,
    annotation_column = "annotation_custom",
    sample_column = "sample_custom",
    condition_column = "condition_custom",
    reduction = "umap_custom"
  )

  expect_equal(detected$annotation_column, "annotation_custom")
  expect_equal(detected$sample_column, "sample_custom")
  expect_equal(detected$condition_column, "condition_custom")
  expect_equal(detected$reduction, "umap_custom")
})

test_that("review_seurat_object works with automatic detection", {
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
  object$sample <- c("S1", "S2")
  object$condition_FIXED2 <- c("STD", "Tx")
  object$celltype_for_R8plot_FIXED2 <- c(
    "Resident Kupffer-like",
    "Monocyte-like"
  )
  object[["umapRPCA"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      dimnames = list(c("cell1", "cell2"), c("UMAP_1", "UMAP_2"))
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  result <- review_seurat_object(object)

  expect_s3_class(result, "uenoy_review")
  expect_equal(
    result$summary$annotation_column,
    "celltype_for_R8plot_FIXED2"
  )
  expect_equal(result$summary$condition_column, "condition_FIXED2")
  expect_equal(result$summary$reduction, "umapRPCA")
})

test_that("automatic detection reports missing annotation columns", {
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
    detect_review_metadata(object),
    "Could not automatically detect an annotation column"
  )
})
