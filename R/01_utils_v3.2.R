`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(..., log_file = NULL) {
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = ""))
  message(msg)
  if (!is.null(log_file)) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  invisible(msg)
}

assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop(
      "Required packages are missing: ", paste(missing, collapse = ", "),
      "\nInstall them within the active renv environment and retry."
    )
  }
}

detect_external_root <- function(candidates) {
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  mounted <- candidates[dir.exists(candidates)]
  if (!length(mounted)) {
    stop(
      "External SSD was not detected.\nChecked:\n- ",
      paste(candidates, collapse = "\n- "),
      "\nMount SSD990_uenoy and rerun."
    )
  }
  mounted[[1]]
}

resolve_project_root <- function(configured_root) {
  if (requireNamespace("here", quietly = TRUE)) {
    root <- tryCatch(here::here(), error = function(e) NULL)
    if (!is.null(root) && file.exists(file.path(root, ".here"))) return(root)
  }
  normalizePath(configured_root, winslash = "/", mustWork = FALSE)
}

ensure_dirs <- function(paths) {
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}

detect_metadata_column <- function(meta, requested = NULL, candidates, label) {
  if (!is.null(requested)) {
    if (!requested %in% colnames(meta)) stop(label, " column not found: ", requested)
    return(requested)
  }
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) {
    stop(
      "Could not auto-detect ", label, " metadata column.\nAvailable columns:\n",
      paste(colnames(meta), collapse = ", ")
    )
  }
  hit[[1]]
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

checkpoint_path <- function(ctx, name, ext = "rds") {
  file.path(ctx$checkpoint_dir, paste0(name, ".", ext))
}

checkpoint_exists <- function(ctx, name, ext = "rds") {
  file.exists(checkpoint_path(ctx, name, ext))
}

save_checkpoint <- function(object, ctx, name) {
  path <- checkpoint_path(ctx, name)
  saveRDS(object, path, compress = FALSE)
  invisible(path)
}

load_checkpoint <- function(ctx, name) {
  readRDS(checkpoint_path(ctx, name))
}

should_run <- function(ctx, name, ext = "rds") {
  if (isTRUE(ctx$cfg$overwrite)) return(TRUE)
  if (!isTRUE(ctx$cfg$resume)) return(TRUE)
  !checkpoint_exists(ctx, name, ext)
}

write_session_info <- function(path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
}
