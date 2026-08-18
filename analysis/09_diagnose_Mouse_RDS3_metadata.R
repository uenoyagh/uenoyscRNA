############################################################
## 09_diagnose_Mouse_RDS3_metadata.R
############################################################

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_metadata_diagnosis"
)

required_packages <- c("Seurat", "SeuratObject", "data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
value_dir <- file.path(output_dir, "08_candidate_annotation_value_counts")
dir.create(value_dir, recursive = TRUE, showWarnings = FALSE)

write_log <- function(x) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", x)
  cat(line, "\n")
  cat(line, "\n", file = file.path(output_dir, "diagnosis.log"), append = TRUE)
}

safe_unique_count <- function(x) length(unique(x[!is.na(x)]))

safe_preview <- function(x, n = 50L) {
  values <- unique(as.character(x[!is.na(x)]))
  values <- values[nzchar(values)]
  if (length(values) == 0) return("")
  paste(head(values, n), collapse = " | ")
}

pattern_score <- function(x, patterns) {
  values <- unique(as.character(x[!is.na(x)]))
  if (length(values) == 0) return(0L)
  sum(vapply(
    patterns,
    function(p) any(grepl(p, values, ignore.case = TRUE)),
    logical(1)
  ))
}

sanitize_filename <- function(x) {
  substr(gsub("[^A-Za-z0-9._-]+", "_", x), 1, 150)
}

if (!file.exists(rds_file)) stop("RDS not found:\n", rds_file, call. = FALSE)

write_log(paste0("Loading RDS: ", rds_file))
obj <- readRDS(rds_file)
if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object.", call. = FALSE)

write_log(paste0(
  "Cells: ", format(ncol(obj), big.mark = ","),
  " | Features: ", format(nrow(obj), big.mark = ",")
))

md <- obj@meta.data

metadata_columns <- data.frame(
  column = colnames(md),
  class = vapply(md, function(x) paste(class(x), collapse = "/"), character(1)),
  n_non_na = vapply(md, function(x) sum(!is.na(x)), integer(1)),
  n_na = vapply(md, function(x) sum(is.na(x)), integer(1)),
  n_unique_non_na = vapply(md, safe_unique_count, integer(1)),
  stringsAsFactors = FALSE
)
fwrite(metadata_columns, file.path(output_dir, "01_metadata_columns.csv"))

metadata_value_preview <- data.frame(
  column = colnames(md),
  preview = vapply(md, safe_preview, character(1)),
  stringsAsFactors = FALSE
)
fwrite(metadata_value_preview, file.path(output_dir, "02_metadata_value_preview.csv"))

layer1_patterns <- c(
  "Macrophage", "Mphi", "Kupffer", "Monocyte",
  "Hepatocyte", "LSEC", "Endothelial", "HSC",
  "Stellate", "Cholangiocyte", "Neutrophil", "T cell", "B cell"
)

layer2_patterns <- c(
  "Resident.*Kupffer", "Kupffer.*Resident", "Monocyte[- ]?like",
  "Inflammatory.*M1", "M1[- ]?like", "Pro[- ]?resolution",
  "M2[- ]?like", "SPP1", "TREM2", "MASH[- ]?associated"
)

condition_patterns <- c("^STD$", "CDAHFD", "CDHFD", "Sham", "^Tx", "Control")
sample_patterns <- c("rep[0-9]+", "STD", "CDAHFD", "CDHFD", "Sham", "^Tx")

candidate_scores <- data.frame(
  column = colnames(md),
  class = vapply(md, function(x) paste(class(x), collapse = "/"), character(1)),
  n_unique_non_na = vapply(md, safe_unique_count, integer(1)),
  layer1_score = vapply(md, pattern_score, integer(1), patterns = layer1_patterns),
  layer2_score = vapply(md, pattern_score, integer(1), patterns = layer2_patterns),
  condition_score = vapply(md, pattern_score, integer(1), patterns = condition_patterns),
  sample_score = vapply(md, pattern_score, integer(1), patterns = sample_patterns),
  stringsAsFactors = FALSE
)

candidate_scores <- candidate_scores[
  order(
    -candidate_scores$layer2_score,
    -candidate_scores$layer1_score,
    -candidate_scores$condition_score,
    -candidate_scores$sample_score,
    candidate_scores$column
  ),
]
fwrite(candidate_scores, file.path(output_dir, "03_candidate_column_scores.csv"))

reduction_names <- Reductions(obj)
reductions_table <- data.frame(
  reduction = reduction_names,
  dimensions = vapply(reduction_names, function(r) ncol(Embeddings(obj, r)), integer(1)),
  cells = vapply(reduction_names, function(r) nrow(Embeddings(obj, r)), integer(1)),
  stringsAsFactors = FALSE
)
fwrite(reductions_table, file.path(output_dir, "04_reductions.csv"))

condition_candidates <- candidate_scores$column[candidate_scores$condition_score > 0]
if (length(condition_candidates) > 0) {
  condition_long <- rbindlist(lapply(condition_candidates, function(col) {
    data.table(column = col, value = as.character(md[[col]]))[
      , .N, by = .(column, value)
    ][order(column, -N)]
  }), fill = TRUE)
  fwrite(condition_long, file.path(output_dir, "05_condition_crosscheck.csv"))
}

sample_candidates <- candidate_scores$column[candidate_scores$sample_score > 0]
if (length(sample_candidates) > 0) {
  sample_long <- rbindlist(lapply(sample_candidates, function(col) {
    data.table(column = col, value = as.character(md[[col]]))[
      , .N, by = .(column, value)
    ][order(column, -N)]
  }), fill = TRUE)
  fwrite(sample_long, file.path(output_dir, "06_sample_crosscheck.csv"))
}

annotation_candidates <- unique(c(
  candidate_scores$column[candidate_scores$layer2_score > 0],
  candidate_scores$column[candidate_scores$layer1_score > 0],
  grep("celltype|annotation|layer|mphi|macro|cluster",
       colnames(md), value = TRUE, ignore.case = TRUE)
))

for (col in annotation_candidates) {
  counts <- data.table(value = as.character(md[[col]]))[
    , .N, by = value
  ][order(-N)]
  fwrite(
    counts,
    file.path(value_dir, paste0(sanitize_filename(col), "_value_counts.csv"))
  )
}

summary_file <- file.path(output_dir, "diagnosis_summary.txt")
sink(summary_file)
cat("Mouse RDS3 metadata diagnosis\n\n")
cat("RDS:\n", rds_file, "\n\n", sep = "")
cat("Cells: ", format(ncol(obj), big.mark = ","), "\n", sep = "")
cat("Features: ", format(nrow(obj), big.mark = ","), "\n", sep = "")
cat("Metadata columns: ", ncol(md), "\n\n", sep = "")
cat("Reductions:\n")
print(reductions_table)
cat("\nTop Layer2 candidates:\n")
print(head(candidate_scores[candidate_scores$layer2_score > 0, ], 10))
cat("\nTop Layer1 candidates:\n")
print(head(candidate_scores[candidate_scores$layer1_score > 0, ], 10))
cat("\nOutput directory:\n", output_dir, "\n", sep = "")
sink()

write_log("Metadata diagnosis completed successfully")

cat("\n====================================================\n")
cat("Mouse RDS3 metadata diagnosis completed\n")
cat("====================================================\n")
cat("Output directory:\n", output_dir, "\n\n", sep = "")
cat("Review first:\n")
cat("  03_candidate_column_scores.csv\n")
cat("  diagnosis_summary.txt\n")
cat("  08_candidate_annotation_value_counts/\n")
cat("====================================================\n")
