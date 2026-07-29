metadata_column_summary <- function(meta, max_unique_values = 50) {
  do.call(rbind, lapply(names(meta), function(nm) {
    x <- meta[[nm]]
    u <- unique(as.character(x[!is.na(x)]))
    examples <- if (length(u) > 0 && length(u) <= max_unique_values)
      paste(head(u, 10), collapse=" | ") else NA_character_
    data.frame(
      column=nm, class=paste(class(x), collapse="/"),
      n_non_na=sum(!is.na(x)), n_na=sum(is.na(x)),
      n_unique=length(u), example_values=examples,
      stringsAsFactors=FALSE
    )
  }))
}

detect_metadata_candidates <- function(meta_names) {
  pats <- list(
    sample=c("^sample$","sample_id","orig.ident","donor","patient","replicate"),
    condition=c("^condition$","group","treatment","diet","disease","status"),
    cluster=c("^seurat_clusters$","snn_res","res\\.","cluster"),
    annotation=c("celltype","cell_type","annotation","annot","identity","label")
  )
  lapply(pats, function(pp) unique(unlist(lapply(pp, function(p)
    grep(p, meta_names, ignore.case=TRUE, value=TRUE)))))
}
