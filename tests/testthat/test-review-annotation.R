test_that("review_umap creates annotation, sample, and condition plots", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ggplot2")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 2),
    j = c(1, 1, 2, 3, 4),
    x = c(1, 2, 3, 1, 2),
    dims = c(3, 4),
    dimnames = list(
      c("Clec4f", "Timd4", "Ccr2"),
      paste0("cell", 1:4)
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- c("S1", "S1", "S2", "S2")
  object$condition_FIXED2 <- c("STD", "STD", "Tx", "Tx")
  object$celltype_for_R8plot_FIXED2 <- c(
    "Resident Kupffer-like",
    "Resident Kupffer-like",
    "Monocyte-like",
    "Monocyte-like"
  )

  object[["umapRPCA"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(
        0, 0,
        1, 0,
        0, 1,
        1, 1
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(
        paste0("cell", 1:4),
        c("UMAP_1", "UMAP_2")
      )
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  result <- review_umap(object)

  expect_s3_class(result, "uenoy_review_umap")
  expect_true(all(c("annotation", "sample", "condition") %in% names(result$plots)))
  expect_s3_class(result$plots$annotation, "ggplot")
  expect_equal(result$settings$reduction, "umapRPCA")
})

test_that("review_umap writes PDF files", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ggplot2")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 1),
    dims = c(2, 2),
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- c("S1", "S2")
  object$condition <- c("STD", "Tx")
  object$celltype <- c("Kupffer", "Monocyte")

  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("cell1", "cell2"),
        c("UMAP_1", "UMAP_2")
      )
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  output_dir <- tempfile("review-umap-")

  result <- review_umap(
    object,
    output_dir = output_dir
  )

  expect_true(file.exists(result$files[["annotation"]]))
  expect_true(file.exists(result$files[["sample"]]))
  expect_true(file.exists(result$files[["condition"]]))
})

test_that("review_annotation creates structured output", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ggplot2")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 2),
    j = c(1, 1, 2, 2),
    x = c(1, 2, 3, 1),
    dims = c(2, 2),
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- c("S1", "S2")
  object$condition <- c("STD", "Tx")
  object$celltype <- c("Kupffer", "Monocyte")

  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("cell1", "cell2"),
        c("UMAP_1", "UMAP_2")
      )
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  output_dir <- tempfile("annotation-review-")

  result <- review_annotation(
    object,
    output_dir = output_dir
  )

  expect_s3_class(result, "uenoy_annotation_review")
  expect_true(dir.exists(file.path(output_dir, "Summary")))
  expect_true(dir.exists(file.path(output_dir, "UMAP")))
  expect_true(file.exists(result$manifest_path))
  expect_true(nrow(result$manifest) >= 3)
})

test_that("review_umap rejects incomplete named palettes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(1, 2),
    x = c(1, 1),
    dims = c(2, 2),
    dimnames = list(
      c("Clec4f", "Timd4"),
      c("cell1", "cell2")
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$celltype <- c("Kupffer", "Monocyte")

  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("cell1", "cell2"),
        c("UMAP_1", "UMAP_2")
      )
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  expect_error(
    review_umap(
      object,
      palette = c(Kupffer = "#000000")
    ),
    "missing colors"
  )
})
