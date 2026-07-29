# CHANGELOG v3.1

## Changed
- Cell Fraction metadata resolver now supports exact candidates, regular-expression candidates, and explicit overrides.
- Highest numeric Seurat clustering resolution is selected when several pattern matches exist.
- Annotation candidates now include `feature_annotation_added`.
- Review metadata detection also recognizes `feature_annotation_added`.
- Cell Fraction outputs detailed metadata detection reports.
- Analysis root can be set through `UENOY_SCRNA_ROOT`.

## Added
- `R/metadata_detection.R`
- `detect_metadata()`
- `write_metadata_detection_report()`
- Metadata helper regression tests
- v3.1 README and migration guide

## Compatibility
- Existing Cell Fraction profile fields remain valid.
- New pattern fields are optional.
- Existing explicit column overrides retain highest priority.
