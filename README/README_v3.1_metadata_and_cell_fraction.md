# uenoyscRNA v3.1: metadata recognition and Cell Fraction

## Main changes

- `feature_annotation_added` is recognized as a cell-type annotation.
- Seurat clustering columns such as `integratedRPCA_snn_res.0.8` are found by pattern.
- When multiple clustering resolutions are present, the highest numeric resolution is selected unless an override is supplied.
- Cell Fraction writes `*_metadata_detection.csv` under `08_QC`.
- New package API: `detect_metadata()` and `write_metadata_detection_report()`.
- `analysis/06_cell_fraction_transition.R` accepts `UENOY_SCRNA_ROOT` and falls back to the current project directory.

## Recommended first validation

```r
setwd("/Users/uenoya/Projects/uenoyscRNA")
source("analysis/06_cell_fraction_transition.R")
```

Review the generated `08_QC/*_metadata_detection.csv` before interpreting fraction plots.

## Explicit override

Edit the relevant profile in `config/cell_fraction_config.R`:

```r
feature_column_override = "feature_annotation_added"
cluster_column_override  # not used; set feature_column_override to the desired cluster column
```

For macrophage cluster analysis, for example:

```r
feature_column_override = "integratedRPCA_snn_res.0.8"
```
