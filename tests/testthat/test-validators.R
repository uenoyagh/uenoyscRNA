test_that(".validate_seurat_object validates Seurat objects", {

  object <- SeuratObject::pbmc_small

  expect_true(
    .validate_seurat_object(object)
  )

  expect_error(
    .validate_seurat_object(data.frame()),
    "`object` must be a Seurat object.",
    fixed = TRUE
  )
})


test_that(".validate_character_vector validates character vectors", {

  expect_true(
    .validate_character_vector(
      c("A", "B"),
      arg = "features"
    )
  )

  expect_error(
    .validate_character_vector(
      1:3,
      arg = "features"
    ),
    "`features` must be a character vector.",
    fixed = TRUE
  )

  expect_error(
    .validate_character_vector(
      character(0),
      arg = "features"
    ),
    "`features` must not be empty.",
    fixed = TRUE
  )

  expect_true(
    .validate_character_vector(
      character(0),
      arg = "features",
      allow_empty = TRUE
    )
  )

  expect_error(
    .validate_character_vector(
      c("A", NA_character_),
      arg = "features"
    ),
    "`features` must not contain NA values.",
    fixed = TRUE
  )
})


test_that(".validate_single_string validates scalar strings", {

  expect_true(
    .validate_single_string(
      "RNA",
      arg = "assay"
    )
  )

  expect_true(
    .validate_single_string(
      NULL,
      arg = "assay",
      allow_null = TRUE
    )
  )

  expect_error(
    .validate_single_string(
      NULL,
      arg = "assay"
    ),
    "`assay` must be a single character string.",
    fixed = TRUE
  )

  expect_error(
    .validate_single_string(
      c("RNA", "SCT"),
      arg = "assay"
    ),
    "`assay` must be a single character string.",
    fixed = TRUE
  )

  expect_error(
    .validate_single_string(
      "",
      arg = "assay"
    ),
    "`assay` must be a single character string.",
    fixed = TRUE
  )
})


test_that(".validate_logical validates logical scalars", {

  expect_true(
    .validate_logical(
      TRUE,
      arg = "raster"
    )
  )

  expect_true(
    .validate_logical(
      NULL,
      arg = "raster",
      allow_null = TRUE
    )
  )

  expect_error(
    .validate_logical(
      "yes",
      arg = "raster"
    ),
    "`raster` must be TRUE or FALSE.",
    fixed = TRUE
  )

  expect_error(
    .validate_logical(
      c(TRUE, FALSE),
      arg = "raster"
    ),
    "`raster` must be TRUE or FALSE.",
    fixed = TRUE
  )

  expect_error(
    .validate_logical(
      NA,
      arg = "raster"
    ),
    "`raster` must be TRUE or FALSE.",
    fixed = TRUE
  )
})


test_that(".validate_number validates numeric scalars", {

  expect_true(
    .validate_number(
      1,
      arg = "pt.size"
    )
  )

  expect_true(
    .validate_number(
      NULL,
      arg = "pt.size",
      allow_null = TRUE
    )
  )

  expect_true(
    .validate_number(
      0.5,
      arg = "alpha",
      min = 0,
      max = 1
    )
  )

  expect_error(
    .validate_number(
      "1",
      arg = "pt.size"
    ),
    "`pt.size` must be a single numeric value.",
    fixed = TRUE
  )

  expect_error(
    .validate_number(
      Inf,
      arg = "pt.size"
    ),
    "`pt.size` must be finite.",
    fixed = TRUE
  )

  expect_error(
    .validate_number(
      2,
      arg = "alpha",
      min = 0,
      max = 1
    ),
    "`alpha` must be between 0 and 1.",
    fixed = TRUE
  )
})


test_that(".validate_assay validates assay names", {

  object <- SeuratObject::pbmc_small

  expect_true(
    .validate_assay(
      object,
      "RNA",
      allow_null = FALSE
    )
  )

  expect_true(
    .validate_assay(
      object,
      NULL,
      allow_null = TRUE
    )
  )

  expect_error(
    .validate_assay(
      object,
      "missing_assay",
      allow_null = FALSE
    ),
    "Assay `missing_assay` was not found in `object`.",
    fixed = TRUE
  )
})


test_that(".validate_features validates feature names", {

  object <- SeuratObject::pbmc_small
  features <- head(rownames(object[["RNA"]]), 2)

  expect_true(
    .validate_features(
      object,
      features,
      assay = "RNA"
    )
  )

  expect_error(
    .validate_features(
      object,
      c(features, "NOT_A_GENE"),
      assay = "RNA"
    ),
    "The following features were not found in assay `RNA`: NOT_A_GENE",
    fixed = TRUE
  )

  expect_true(
    .validate_features(
      object,
      c(features, "NOT_A_GENE"),
      assay = "RNA",
      require_all = FALSE
    )
  )
})


test_that(".validate_metadata_column validates metadata columns", {

  object <- SeuratObject::pbmc_small
  metadata_column <- colnames(object[[]])[1]

  expect_true(
    .validate_metadata_column(
      object,
      metadata_column,
      arg = "group.by",
      allow_null = FALSE
    )
  )

  expect_true(
    .validate_metadata_column(
      object,
      NULL,
      arg = "group.by",
      allow_null = TRUE
    )
  )

  expect_error(
    .validate_metadata_column(
      object,
      "missing_column",
      arg = "group.by",
      allow_null = FALSE
    ),
    "Metadata column `missing_column` was not found in `object`.",
    fixed = TRUE
  )
})

test_that(".validate_group_order validates metadata group order", {

  object <- SeuratObject::pbmc_small

  metadata_column <- colnames(object[[]])[1]

  observed_values <- unique(
    as.character(object[[]][[metadata_column]])
  )

  observed_values <- observed_values[
    !is.na(observed_values)
  ]

  valid_order <- head(
    observed_values,
    min(2L, length(observed_values))
  )

  expect_true(
    .validate_group_order(
      object = object,
      column = metadata_column,
      order = valid_order,
      arg = "group_order"
    )
  )

  expect_true(
    .validate_group_order(
      object = object,
      column = metadata_column,
      order = NULL,
      arg = "group_order",
      allow_null = TRUE
    )
  )

  expect_error(
    .validate_group_order(
      object = object,
      column = metadata_column,
      order = NULL,
      arg = "group_order",
      allow_null = FALSE
    ),
    "`group_order` must be a character vector.",
    fixed = TRUE
  )

  expect_error(
    .validate_group_order(
      object = object,
      column = metadata_column,
      order = c(valid_order[1], valid_order[1]),
      arg = "group_order"
    ),
    "`group_order` must not contain duplicated values.",
    fixed = TRUE
  )

  expect_error(
    .validate_group_order(
      object = object,
      column = metadata_column,
      order = "NOT_A_GROUP",
      arg = "group_order"
    ),
    paste0(
      "The following values in `group_order` were not found in metadata column `",
      metadata_column,
      "`: NOT_A_GROUP"
    ),
    fixed = TRUE
  )
})

