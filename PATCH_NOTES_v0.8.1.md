# uenoyscRNA 0.8.1

Changes made from 0.8.0:

- Restored the legacy `review_umap()` return structure and class expected by the package tests.
- Preserved the newer split-panel `review_umap()` mode through automatic dispatch when `group_by`, `split_by`, or advanced arguments are supplied.
- Added source-RDS and analysis-date footers to legacy UMAP outputs.
- Added `require_reduction` support to `detect_review_settings()`.
- Prevented `review_marker_dotplot()` from requiring a dimensional reduction.
- Preserved reduction requirements for `review_featureplot()`.
- Forced annotation-registry CSV columns to character on import, preserving cluster identifiers.
- Exported `theme_ueno_scRNA()`.
- Declared `ggrastr` and `Matrix` in Suggests.
- Expanded `.Rbuildignore` for development-only top-level files.
- Updated package version to 0.8.1.
