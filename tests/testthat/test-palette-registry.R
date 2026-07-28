test_that("mouse liver palette can be loaded", {
  palette <- get_celltype_palette(
    palette = "mouse_liver"
  )

  expect_type(palette, "character")
  expect_true(length(palette) > 0L)
  expect_false(is.null(names(palette)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", palette)))
  expect_true("Hepatocyte" %in% names(palette))
  expect_true("Kupffer cell" %in% names(palette))
  expect_true("Unclassified" %in% names(palette))
})

test_that("Unclassified can be excluded", {
  palette <- get_celltype_palette(
    palette = "mouse_liver",
    include_unclassified = FALSE
  )

  expect_false("Unclassified" %in% names(palette))
})

test_that("unsupported palette names are rejected", {
  expect_error(
    get_celltype_palette(
      palette = "unsupported_palette"
    ),
    "Unsupported palette"
  )
})

test_that("named palettes are validated", {
  palette <- c(
    Hepatocyte = "#D9A441",
    Macrophage = "#3C5488"
  )

  expect_invisible(validate_named_palette(palette))
})

test_that("unnamed palettes are rejected", {
  palette <- c("#D9A441", "#3C5488")

  expect_error(
    validate_named_palette(palette),
    "named character vector"
  )
})

test_that("invalid hexadecimal colors are rejected", {
  palette <- c(
    Hepatocyte = "gold",
    Macrophage = "#3C5488"
  )

  expect_error(
    validate_named_palette(palette),
    "hexadecimal colors"
  )
})

test_that("missing categories are detected", {
  palette <- c(Hepatocyte = "#D9A441")

  expect_error(
    validate_named_palette(
      palette = palette,
      categories = c("Hepatocyte", "Kupffer cell")
    ),
    "missing categories"
  )
})

test_that("unexpected categories can be rejected", {
  palette <- c(
    Hepatocyte = "#D9A441",
    Macrophage = "#3C5488"
  )

  expect_error(
    validate_named_palette(
      palette = palette,
      categories = "Hepatocyte",
      allow_extra = FALSE
    ),
    "unexpected categories"
  )
})
