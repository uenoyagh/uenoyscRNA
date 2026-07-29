# uenoyscRNA v4.0.0

## Main change: data-driven RDS registry

- Removed the assumption that each analysis class must live in a same-named directory.
- Added project-wide, non-loading RDS discovery in `R/rds_registry.R`.
- Added explicit preferred files plus pattern-based fallback rules for:
  - `Mouse_MASH_RDS`
  - `Mouse_MASH_Mphi_RDS`
  - `Human_MASH_RDS`
  - `Human_MASH_Mphi_RDS`
- Excludes `renv`, package metadata RDS files, inventory output, and common analysis bundles.
- Cell Fraction now resolves files through the dataset registry instead of calling
  `list.files()` in nonexistent fixed directories.
- Project inventory scans all candidate paths but deeply opens only the four resolved
  registered analysis objects.
- Added registry summary CSV and project RDS catalog.

## Current preferred files

The defaults reflect the project paths supplied on 2026-07-30. They can be edited in
`config/project_config.R`; environment variables can override both framework and data roots.

## Compatibility

- Existing dataset keys are preserved.
- `get_dataset_dir()`, `get_result_root()`, and `get_result_dir()` remain available.
- New code should use `get_dataset_files()`.
