# uenoy scRNAseq Framework v3.0 — Generic Cell-Fraction Module

## Purpose

This module replaces the macrophage-specific v2.6.1 module with one generic
pipeline that can process:

- `Mouse_MASH_RDS`
- `Mouse_MASH_Mphi_RDS`
- `Human_MASH_RDS`
- `Human_MASH_Mphi_RDS`

It supports whole-cell-type composition and subset/cluster composition using
dataset-specific profiles.

## Files

```text
analysis/06_cell_fraction_transition.R
R/cell_fraction_engine.R
R/cell_fraction_plot.R
config/cell_fraction_config.R
INSTALL_v3.0.R
```

## Installation

1. Extract the ZIP.
2. Open `INSTALL_v3.0.R`.
3. Replace `source_dir` with the extracted directory.
4. Run the installer.

The installer overwrites only the four v3.0 module files. It does not modify
`project_config.R`.

## Execution

```r
source(
  "/Users/uenoya/Projects/uenoyscRNA/analysis/06_cell_fraction_transition.R"
)
```

## Output

For each dataset target:

```text
09_cell_fraction/
└── cell_fraction_transition_v3.0/
    ├── 01_line_total_fraction/
    ├── 02_line_within_parent/
    ├── 03_line_cell_count/
    ├── 04_stacked_feature/
    ├── 05_stacked_parent/
    ├── 06_heatmap/
    ├── 07_tables/
    └── 08_QC/
```

PDF and 300 dpi PNG are produced when enabled. Tables are exported as CSV.

## Profiles

Profiles are defined in `config/cell_fraction_config.R`.

- Whole-cell profile:
  - feature = cell type
  - parent = cell type
  - denominator = total cells
- Macrophage profile:
  - feature = cluster
  - parent = macrophage annotation
  - denominator = total macrophages

After the first run, inspect:

```text
07_tables/*_resolved_metadata_columns.csv
08_QC/*_QC_summary.csv
```

If automatic metadata detection is incorrect, set explicit overrides in the
corresponding profile.

## Overwriting

Set the following in `config/local_config.R`:

```r
overwrite_existing <- TRUE
```

## Important validation note

This module was constructed to fit the known framework interfaces and metadata
conventions, but it has not been executed against the user's local RDS files.
The first run should therefore be treated as validation, especially for Human
NAS and cell-type column detection.
