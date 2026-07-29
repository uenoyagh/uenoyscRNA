make_review_test_object <- function() {
  counts <- matrix(
    c(
      5, 0, 1, 0, 2, 1,
      0, 4, 0, 3, 1, 0,
      2, 2, 1, 0, 0, 3,
      0, 1, 0, 1, 0, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("Gene", 1:4),
      paste0("Cell", 1:6)
    )
  )

  object <- Seurat::CreateSeuratObject(counts = counts)
  object$annotation <- factor(
    c("Type_A", "Type_A", "Type_A", "Type_B", "Type_B", "Type_B"),
    levels = c("Type_A", "Type_B")
  )
  object$sample <- c("S1", "S1", "S2", "S2", "S3", "S3")
  object$condition <- c("Control", "Control", "Control", "Treated", "Treated", "Treated")

  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object
}
