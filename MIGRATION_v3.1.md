# Migration guide: v3.0 to v3.1

1. Back up the current `/Users/uenoya/Projects/uenoyscRNA` folder.
2. Replace it with the supplied `uenoyscRNA_v3.1` folder.
3. Keep or restore your machine-specific `config/local_config.R` if it contains paths not present in the supplied copy.
4. Open `uenoyscRNA.Rproj`.
5. Run `renv::restore()` only when R reports missing packages.
6. Run `source("analysis/06_cell_fraction_transition.R")`.
7. Inspect `08_QC/*_metadata_detection.csv`.

No RDS data are moved or modified by this replacement. Dataset paths remain controlled by `config/local_config.R` and `config/project_config.R`.
