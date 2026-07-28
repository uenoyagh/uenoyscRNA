#' Read a packaged cell-type color palette
#'
#' Reads a named color vector from a CSV file stored in
#' `inst/extdata` within the uenoyscRNA package.
#'
#' @param palette Character scalar specifying the palette name.
#'   Currently, `"mouse_liver"` is supported.
#' @param include_unclassified Logical. If `FALSE`, the
#'   `"Unclassified"` entry is removed.
#'
#' @return A named character vector containing hexadecimal colors.
#'
#' @examples
#' \dontrun{
#' colors <- get_celltype_palette("mouse_liver")
#' colors["Hepatocyte"]
#' }
#'
#' @export
get_celltype_palette <- function(
    palette = "mouse_liver",
    include_unclassified = TRUE
) {
  if (
    length(palette) != 1L ||
    is.na(palette) ||
    !is.character(palette) ||
    !nzchar(palette)
  ) {
    stop(
      "`palette` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (
    length(include_unclassified) != 1L ||
    is.na(include_unclassified) ||
    !is.logical(include_unclassified)
  ) {
    stop(
      "`include_unclassified` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  palette_files <- c(
    mouse_liver = "mouse_liver_celltype_palette.csv"
  )

  if (!palette %in% names(palette_files)) {
    stop(
      "Unsupported palette: ",
      palette,
      ". Supported palettes: ",
      paste(names(palette_files), collapse = ", "),
      call. = FALSE
    )
  }

  palette_path <- system.file(
    "extdata",
    unname(palette_files[[palette]]),
    package = "uenoyscRNA"
  )

  if (!nzchar(palette_path)) {
    stop(
      "Palette file was not found in the installed package: ",
      unname(palette_files[[palette]]),
      call. = FALSE
    )
  }

  palette_table <- utils::read.csv(
    palette_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )

  required_columns <- c("celltype", "color")
  missing_columns <- setdiff(required_columns, colnames(palette_table))

  if (length(missing_columns) > 0L) {
    stop(
      "Palette file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  palette_table$celltype <- trimws(as.character(palette_table$celltype))
  palette_table$color <- trimws(as.character(palette_table$color))

  invalid_celltype <- is.na(palette_table$celltype) |
    !nzchar(palette_table$celltype)

  if (any(invalid_celltype)) {
    stop(
      "Palette file contains empty cell-type names.",
      call. = FALSE
    )
  }

  invalid_color <- is.na(palette_table$color) |
    !grepl("^#[0-9A-Fa-f]{6}$", palette_table$color)

  if (any(invalid_color)) {
    stop(
      "Palette file contains invalid hexadecimal colors. ",
      "Colors must use the form #RRGGBB.",
      call. = FALSE
    )
  }

  if (anyDuplicated(palette_table$celltype)) {
    duplicated_celltypes <- unique(
      palette_table$celltype[duplicated(palette_table$celltype)]
    )

    stop(
      "Palette file contains duplicated cell types: ",
      paste(duplicated_celltypes, collapse = ", "),
      call. = FALSE
    )
  }

  if (!include_unclassified) {
    palette_table <- palette_table[
      palette_table$celltype != "Unclassified",
      ,
      drop = FALSE
    ]
  }

  palette_vector <- stats::setNames(
    palette_table$color,
    palette_table$celltype
  )

  validate_named_palette(palette_vector)

  palette_vector
}


#' Validate a named color palette
#'
#' Checks that a palette is a named character vector containing valid
#' hexadecimal colors and, optionally, verifies that required categories
#' are represented.
#'
#' @param palette Named character vector of hexadecimal colors.
#' @param categories Optional character vector containing categories that
#'   must be present in the palette.
#' @param allow_extra Logical. If `FALSE`, palette entries not contained
#'   in `categories` cause an error.
#'
#' @return Invisibly returns `TRUE` when validation succeeds.
#'
#' @examples
#' validate_named_palette(
#'   c(
#'     Hepatocyte = "#D9A441",
#'     Macrophage = "#3C5488"
#'   )
#' )
#'
#' @export
validate_named_palette <- function(
    palette,
    categories = NULL,
    allow_extra = TRUE
) {
  if (!is.character(palette)) {
    stop(
      "`palette` must be a character vector.",
      call. = FALSE
    )
  }

  if (is.null(names(palette))) {
    stop(
      "`palette` must be a named character vector.",
      call. = FALSE
    )
  }

  if (length(palette) == 0L) {
    stop(
      "`palette` must contain at least one color.",
      call. = FALSE
    )
  }

  invalid_names <- is.na(names(palette)) |
    !nzchar(names(palette))

  if (any(invalid_names)) {
    stop(
      "All palette entries must have non-empty names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(names(palette))) {
    duplicated_names <- unique(
      names(palette)[duplicated(names(palette))]
    )

    stop(
      "Palette names must be unique. Duplicated names: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }

  invalid_colors <- is.na(palette) |
    !grepl("^#[0-9A-Fa-f]{6}$", palette)

  if (any(invalid_colors)) {
    stop(
      "All palette values must be hexadecimal colors ",
      "in the form #RRGGBB.",
      call. = FALSE
    )
  }

  if (
    length(allow_extra) != 1L ||
    is.na(allow_extra) ||
    !is.logical(allow_extra)
  ) {
    stop(
      "`allow_extra` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.null(categories)) {
    categories <- unique(
      trimws(as.character(categories[!is.na(categories)]))
    )
    categories <- categories[nzchar(categories)]

    missing_categories <- setdiff(categories, names(palette))

    if (length(missing_categories) > 0L) {
      stop(
        "Palette is missing categories: ",
        paste(missing_categories, collapse = ", "),
        call. = FALSE
      )
    }

    if (!allow_extra) {
      extra_categories <- setdiff(names(palette), categories)

      if (length(extra_categories) > 0L) {
        stop(
          "Palette contains unexpected categories: ",
          paste(extra_categories, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}
