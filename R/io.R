list_rds_files <- function(target_dir, recursive = FALSE) {
  if (!dir.exists(target_dir)) stop("Directory does not exist: ", target_dir)
  x <- list.files(target_dir, pattern="\\.[Rr][Dd][Ss]$", full.names=TRUE, recursive=recursive)
  x[order(tolower(basename(x)))]
}

safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error=function(e)
    structure(list(path=path, message=conditionMessage(e)), class="rds_read_error"))
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(x, path, row.names=FALSE, fileEncoding="UTF-8")
  invisible(path)
}
