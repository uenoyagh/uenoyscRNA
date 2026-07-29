# Migration to v4.0

1. Back up the existing framework directory.
2. Replace it with this distribution and keep the folder name `uenoyscRNA`.
3. Open `uenoyscRNA.Rproj`.
4. Run:

```r
source("analysis/01_inventory_RDS.R")
```

Confirm that all four registered datasets show `[OK]` and inspect:

```text
/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/
  analysis_results/00_project_inventory/dataset_registry_summary.csv
```

5. Run Cell Fraction:

```r
source("analysis/06_cell_fraction_transition.R")
```

Outputs use the run folder `cell_fraction_transition_v4.0`.

## Changing an RDS

Edit only `preferred_files` in `config/project_config.R`. The former requirement to create
`Mouse_MASH_Mphi_RDS`, `Human_MASH_RDS`, or `Human_MASH_Mphi_RDS` directories no longer applies.
