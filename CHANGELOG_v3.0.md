# CHANGELOG — Framework v3.0

## v3.0

- Retired the macrophage-specific v2.6.1 entry point.
- Added generic `06_cell_fraction_transition.R`.
- Added dataset-specific profiles for all four active RDS directories.
- Supports:
  - whole-cell composition,
  - macrophage cluster composition,
  - total-fraction lines,
  - within-parent fraction lines,
  - raw cell-count lines,
  - stacked feature plots,
  - stacked parent plots,
  - blue-white-red heatmaps,
  - PDF, PNG and CSV output.
- Avoids unregistered `get_result_dir(..., "cluster_fraction")`.
- Uses `get_result_root()` and creates `09_cell_fraction` directly.
- Does not modify `project_config.R`.
