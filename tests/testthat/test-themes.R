test_that("theme_ueno_scRNA returns a ggplot2 theme", {
  theme <- theme_ueno_scRNA()

  expect_s3_class(theme, "theme")
})

test_that("theme_ueno_scRNA validates base_size", {
  expect_error(
    theme_ueno_scRNA(base_size = 0),
    "`base_size` must be a single positive number"
  )

  expect_error(
    theme_ueno_scRNA(base_size = NA_real_),
    "`base_size` must be a single positive number"
  )
})

test_that("theme_ueno_scRNA validates base_family", {
  expect_error(
    theme_ueno_scRNA(base_family = NA_character_),
    "`base_family` must be a single character string"
  )

  expect_error(
    theme_ueno_scRNA(base_family = c("Arial", "Helvetica")),
    "`base_family` must be a single character string"
  )
})
