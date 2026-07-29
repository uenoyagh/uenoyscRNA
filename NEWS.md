# uenoyscRNA 0.9.0 / Framework v3.1

- Added shared metadata recognition with explicit overrides.
- Added `feature_annotation_added` support.
- Added pattern-based Seurat cluster-column detection.
- Added metadata detection QC reports to Cell Fraction.
- Improved Framework-root resolution in the Cell Fraction script.
- Added regression tests for resolution and cluster-label ordering.

# uenoyscRNA 0.8.0

- Added a dedicated single shared legend to split-panel `review_umap()` output.
- Added RDS filename and analysis-date provenance footers to `review_umap()`.
- Added the same provenance footer interface to `review_featureplot()`.
- Removed executable example scripts from `R/` so package loading never tries
  to read placeholder files.
- Consolidated duplicate `review_umap()` and `review_seurat_object()` definitions.
