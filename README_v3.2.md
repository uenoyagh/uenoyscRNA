# uenoyscRNA RDS3 Annotation Validation v3.2

## Storage policy

Framework and R scripts are stored on the Mac internal SSD:

`/Users/uenoya/Projects/uenoyscRNA/`

RDS files, checkpoints, figures, CSV/Excel files and output RDS files are stored on the external SSD:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/`

The pipeline never copies the source RDS to the internal SSD.

## Default input

`Mouse_MASH_RDS/Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds`

Edit `config/config.R` when the actual RDS is in another data directory.

## Installation into the existing framework

Copy these folders into `/Users/uenoya/Projects/uenoyscRNA/`:

- `config/config.R`
- `R/01_utils_v3.2.R`
- `R/02_markers_v3.2.R`
- `R/03_analysis_functions_v3.2.R`
- `R/04_plot_export_v3.2.R`
- `analysis/RDS3_annotation_validation_v3.2/99_run_all.R`

Create a project-root marker once:

```r
here::i_am("analysis/RDS3_annotation_validation_v3.2/99_run_all.R")
```

Alternatively, create an empty `.here` file at the project root.

## Run

From Terminal:

```bash
cd /Users/uenoya/Projects/uenoyscRNA
Rscript analysis/RDS3_annotation_validation_v3.2/99_run_all.R
```

From RStudio:

```r
source(here::here("analysis", "RDS3_annotation_validation_v3.2", "99_run_all.R"))
```

## Resume behavior

Each major step saves an uncompressed RDS checkpoint on the external SSD.

- `01_object_prepared.rds`
- `02_cluster_markers.rds`
- `03_marker_voting.rds`
- `99_complete.rds`

With `CFG$resume = TRUE`, completed steps are restored rather than recalculated.

To rerun everything, set `CFG$overwrite = TRUE`, or delete the output checkpoint directory.

## Main outputs

- Cluster-level presto markers
- General/Ueno marker voting
- Current annotation composition
- Cluster/current/voted annotation UMAPs
- General and Ueno marker DotPlots
- Multipage Violin PDFs
- Excel workbook
- Validation metadata-added RDS
- Session information and execution log

## Important implementation detail

`future::multisession` is configured for compatible operations, but the main
`presto::wilcoxauc()` call is already vectorized and is run once on the full
sparse matrix. Launching one independent full Seurat copy per cluster would
increase memory use substantially, so v3.2 avoids that design.
