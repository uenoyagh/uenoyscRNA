#' Apply an explicit order to grouping values
#'
#' Converts grouping values into a factor with a specified order.
#'
#' Values not explicitly included in `order` are appended after the
#' requested levels in their original order of appearance.
#'
#' @param x Vector of grouping values.
#' @param order Optional character vector specifying the desired order.
#' @param drop Whether unused factor levels should be removed.
#'
#' @return A factor.
#'
#' @keywords internal
.apply_group_order <- function(
    x,
    order = NULL,
    drop = TRUE
) {

  .validate_logical(
    x = drop,
    arg = "drop"
  )

  if (is.factor(x)) {
    original_levels <- levels(x)
  } else {
    original_levels <- unique(as.character(x))
  }

  original_levels <- original_levels[
    !is.na(original_levels)
  ]

  if (is.null(order)) {

    result <- factor(
      x,
      levels = original_levels
    )

    if (drop) {
      result <- droplevels(result)
    }

    return(result)
  }

  .validate_character_vector(
    x = order,
    arg = "order"
  )

  if (anyDuplicated(order)) {
    stop(
      "`order` must not contain duplicated values.",
      call. = FALSE
    )
  }

  observed_values <- unique(as.character(x))
  observed_values <- observed_values[
    !is.na(observed_values)
  ]

  missing_values <- setdiff(
    order,
    observed_values
  )

  if (length(missing_values) > 0L) {
    stop(
      sprintf(
        "The following values in `order` were not found in `x`: %s",
        paste(missing_values, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  remaining_levels <- setdiff(
    original_levels,
    order
  )

  new_levels <- c(
    order,
    remaining_levels
  )

  result <- factor(
    x,
    levels = new_levels
  )

  if (drop) {
    result <- droplevels(result)
  }

  result
}
