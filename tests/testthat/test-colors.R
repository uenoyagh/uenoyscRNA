test_that("ueno_blue_white_red returns requested number of colors", {
  colors <- ueno_blue_white_red(11)

  expect_type(colors, "character")
  expect_length(colors, 11)
  expect_identical(colors[[1]], "#0033FF")
  expect_identical(colors[[11]], "#FF1A1A")
})

test_that("ueno_blue_white_red returns white at the midpoint for odd n", {
  colors <- ueno_blue_white_red(5)

  expect_identical(colors[[3]], "#FFFFFF")
})

test_that("ueno_blue_white_red validates n", {
  expect_error(
    ueno_blue_white_red(1),
    "`n` must be a single integer"
  )

  expect_error(
    ueno_blue_white_red(5.5),
    "`n` must be a single integer"
  )

  expect_error(
    ueno_blue_white_red(NA_integer_),
    "`n` must be a single integer"
  )
})
