test_that(".apply_group_order applies the requested order", {

  x <- c(
    "Tx",
    "Sham",
    "CDAHFD",
    "STD",
    "Tx"
  )

  result <- .apply_group_order(
    x = x,
    order = c(
      "STD",
      "CDAHFD",
      "Sham",
      "Tx"
    )
  )

  expect_s3_class(
    result,
    "factor"
  )

  expect_identical(
    levels(result),
    c(
      "STD",
      "CDAHFD",
      "Sham",
      "Tx"
    )
  )
})


test_that(".apply_group_order appends unspecified groups", {

  x <- c(
    "Group_C",
    "Group_A",
    "Group_B"
  )

  result <- .apply_group_order(
    x = x,
    order = c(
      "Group_A",
      "Group_B"
    )
  )

  expect_identical(
    levels(result),
    c(
      "Group_A",
      "Group_B",
      "Group_C"
    )
  )
})


test_that(".apply_group_order preserves order when order is NULL", {

  x <- c(
    "Group_C",
    "Group_A",
    "Group_B",
    "Group_C"
  )

  result <- .apply_group_order(
    x = x,
    order = NULL
  )

  expect_identical(
    levels(result),
    c(
      "Group_C",
      "Group_A",
      "Group_B"
    )
  )
})


test_that(".apply_group_order handles factor input", {

  x <- factor(
    c("B", "A", "C"),
    levels = c("A", "B", "C", "D")
  )

  result_drop <- .apply_group_order(
    x = x,
    order = c("C", "A"),
    drop = TRUE
  )

  expect_identical(
    levels(result_drop),
    c("C", "A", "B")
  )

  result_keep <- .apply_group_order(
    x = x,
    order = c("C", "A"),
    drop = FALSE
  )

  expect_identical(
    levels(result_keep),
    c("C", "A", "B", "D")
  )
})


test_that(".apply_group_order validates input", {

  x <- c("A", "B", "C")

  expect_error(
    .apply_group_order(
      x = x,
      order = c("A", "A")
    ),
    "`order` must not contain duplicated values.",
    fixed = TRUE
  )

  expect_error(
    .apply_group_order(
      x = x,
      order = "D"
    ),
    "The following values in `order` were not found in `x`: D",
    fixed = TRUE
  )

  expect_error(
    .apply_group_order(
      x = x,
      order = c("A", "B"),
      drop = "yes"
    ),
    "`drop` must be TRUE or FALSE.",
    fixed = TRUE
  )
})


test_that(".apply_group_order preserves missing values", {

  x <- c(
    "A",
    NA_character_,
    "B"
  )

  result <- .apply_group_order(
    x = x,
    order = c("B", "A")
  )

  expect_identical(
    levels(result),
    c("B", "A")
  )

  expect_true(
    is.na(result[2])
  )
})
