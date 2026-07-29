# uenoyscRNA 0.8.0

- Added a dedicated single shared legend to split-panel `review_umap()` output.
- Added RDS filename and analysis-date provenance footers to `review_umap()`.
- Added the same provenance footer interface to `review_featureplot()`.
- Removed executable example scripts from `R/` so package loading never tries
  to read placeholder files.
- Consolidated duplicate `review_umap()` and `review_seurat_object()` definitions.
