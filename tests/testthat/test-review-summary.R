test_that("review_summary summarizes a Seurat object", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1),
    j = c(1, 1, 2, 3),
    x = c(1, 2, 3, 1),
    dims = c(3, 3),
    dimnames = list(
      c("Clec4f", "Timd4", "Ccr2"),
      c("cell1", "cell2", "cell3")
    )
  )

  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$sample <- c("S1", "S1", "S2")
  object$condition_FIXED2 <- c("STD", "STD", "Tx")
  object$celltype_for_R8plot_FIXED2 <- c(
    "Resident Kupffer-like",
    "Resident Kupffer-like",
    "Monocyte-like"
  )

  object[["umapRPCA"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1, 2, 2),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("cell1", "cell2", "cell3"),
        c("UMAP_1", "UMAP_2")
      )
    ),
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )

  result <- review_summary(object)

  expect_s3_class(result, "uenoy_review_summary")
  expect_equal(result$object$n_cells, 3)
  expect_equal(result$object$n_features, 3)
  expect_equal(result$settings$reduction, "umapRPCA")
  expect_equal(
    result$settings$annotation_column,
    "celltype_for_R8plot_FIXED2"
  )
  expect_equal(result$counts$n_annotations, 2)
  expect_equal(result$counts$n_samples, 2)
  expect_equal(result$counts$n_conditions, 2)
  expect_equal(result$values$conditions, c("STD", "Tx"))
})

test_that("review_summary accepts a uenoy_review object", {
  review <- structure(
    list(
      summary = list(
        n_cells = 10,
        n_features = 100,
        assay = "RNA",
        reduction = "umap",
        annotation_column = "celltype",
        sample_column = NULL,
        condition_column = NULL
      ),
      detected_settings = list(),
      metadata = list(),
      marker_coverage = NULL,
      manifest = list()
    ),
    class = c("uenoy_review", "list")
  )

  result <- review_summary(review)

  expect_s3_class(result, "uenoy_review_summary")
  expect_equal(result$object$n_cells, 10)
  expect_true("No sample column was detected." %in% result$warnings)
  expect_true("No condition column was detected." %in% result$warnings)
})

test_that("low marker coverage creates a warning", {
  review <- structure(
    list(
      summary = list(
        n_cells = 10,
        n_features = 100,
        assay = "RNA",
        reduction = "umap",
        annotation_column = "celltype",
        sample_column = "sample",
        condition_column = "condition"
      ),
      detected_settings = list(),
      metadata = list(),
      marker_coverage = data.frame(
        celltype = c("Kupffer", "Monocyte"),
        requested_n = c(10, 10),
        present_n = c(9, 5),
        missing_n = c(1, 5),
        coverage = c(0.9, 0.5),
        missing_genes = c("Timd4", "Ccr2;Ly6c2"),
        stringsAsFactors = FALSE
      ),
      manifest = list()
    ),
    class = c("uenoy_review", "list")
  )

  result <- review_summary(
    review,
    marker_warning_threshold = 0.8
  )

  expect_equal(nrow(result$low_marker_coverage), 1)
  expect_match(result$warnings[[1]], "1 cell type")
})

test_that("as_review_tables returns exportable tables", {
  summary_object <- structure(
    list(
      object = list(n_cells = 10, n_features = 100),
      settings = list(
        assay = "RNA",
        reduction = "umap",
        annotation_column = "celltype",
        sample_column = "sample",
        condition_column = "condition"
      ),
      counts = list(
        n_annotations = 2,
        n_samples = 2,
        n_conditions = 2
      ),
      values = list(
        annotations = c("A", "B"),
        samples = c("S1", "S2"),
        conditions = c("STD", "Tx")
      ),
      marker_coverage = NULL,
      low_marker_coverage = NULL,
      warnings = character(),
      detected_settings = list(),
      manifest = list()
    ),
    class = c("uenoy_review_summary", "list")
  )

  tables <- as_review_tables(summary_object)

  expect_true(all(c("settings", "values", "warnings") %in% names(tables)))
  expect_equal(nrow(tables$settings), 10)
  expect_equal(nrow(tables$values), 6)
})

test_that("write_review_summary writes CSV files", {
  summary_object <- structure(
    list(
      object = list(n_cells = 10, n_features = 100),
      settings = list(
        assay = "RNA",
        reduction = "umap",
        annotation_column = "celltype",
        sample_column = "sample",
        condition_column = "condition"
      ),
      counts = list(
        n_annotations = 2,
        n_samples = 2,
        n_conditions = 2
      ),
      values = list(
        annotations = c("A", "B"),
        samples = c("S1", "S2"),
        conditions = c("STD", "Tx")
      ),
      marker_coverage = NULL,
      low_marker_coverage = NULL,
      warnings = character(),
      detected_settings = list(),
      manifest = list()
    ),
    class = c("uenoy_review_summary", "list")
  )

  output_dir <- tempfile("review-summary-")
  paths <- write_review_summary(summary_object, output_dir)

  expect_true(file.exists(paths[["settings"]]))
  expect_true(file.exists(paths[["values"]]))
  expect_true(file.exists(paths[["warnings"]]))
})

test_that("invalid marker thresholds are rejected", {
  review <- structure(
    list(
      summary = list(
        n_cells = 10,
        n_features = 100,
        assay = "RNA",
        reduction = "umap",
        annotation_column = "celltype",
        sample_column = NULL,
        condition_column = NULL
      ),
      detected_settings = list(),
      metadata = list(),
      marker_coverage = NULL,
      manifest = list()
    ),
    class = c("uenoy_review", "list")
  )

  expect_error(
    review_summary(review, marker_warning_threshold = 1.5),
    "between 0 and 1"
  )
})
