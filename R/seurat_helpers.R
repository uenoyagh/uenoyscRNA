is_seurat_object <- function(x) inherits(x, "Seurat")
get_assay_names <- function(x) if (is_seurat_object(x)) names(x@assays) else character(0)
get_reduction_names <- function(x) if (is_seurat_object(x)) names(x@reductions) else character(0)
get_default_assay_safe <- function(x) {
  if (!is_seurat_object(x)) return(NA_character_)
  tryCatch(SeuratObject::DefaultAssay(x), error=function(e) NA_character_)
}
get_layer_names <- function(x) {
  if (!is_seurat_object(x)) return(data.frame())
  out <- lapply(get_assay_names(x), function(a) {
    z <- tryCatch(SeuratObject::Layers(x[[a]]), error=function(e) character(0))
    if (!length(z)) z <- NA_character_
    data.frame(assay=a, layer=z, stringsAsFactors=FALSE)
  })
  do.call(rbind, out)
}
