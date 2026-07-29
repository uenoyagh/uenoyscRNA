collapse_or_na <- function(x, sep=" | ") if (!length(x)) NA_character_ else paste(x, collapse=sep)
file_size_gb <- function(path) as.numeric(file.info(path)$size) / 1024^3
timestamp_string <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
