# Validation status for v3.1

## Completed in the build environment

- Reviewed the supplied Framework directory structure.
- Checked the edited R files for balanced parentheses, brackets, braces, quoted strings, and comments.
- Confirmed that all four Cell Fraction profiles contain pattern settings.
- Confirmed that `feature_annotation_added` is included in both Cell Fraction and review metadata detection.
- Confirmed that the new exported functions are listed in `NAMESPACE` and have `.Rd` documentation.
- Removed the stale nested patch directory and macOS `.DS_Store` files from the distributed copy.

## Not completed in the build environment

R and the user's Seurat RDS files were not available in this execution environment. Therefore the following still require local validation:

- R parser/package build check
- `testthat` execution
- Seurat-object metadata detection
- End-to-end Cell Fraction plotting

The first local run should be treated as a validation run. Inspect `08_QC/*_metadata_detection.csv` before biological interpretation.
