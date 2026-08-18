# uenoy scRNAseq Framework v4.0.1 visualization patch

Changes from v4.0:

- darker, higher-saturation R8-style colors
- larger points
- point alpha forced to 1
- 1:1 UMAP aspect ratio retained
- full information-rich legend names retained
- source RDS filename printed below every figure
- figure creation date and time printed below every figure

Caption example:

```text
Source RDS: RDS3_with_visualization_metadata_v4.0.rds | Created: 2026-07-30 23:15:42 JST
```

## Install

```bash
bash install_v4.0.1_visualization.sh
```

## Run

```bash
cd "/Users/uenoya/Projects/uenoyscRNA"
Rscript analysis/RDS3_annotation_visualization_v4.0.1/99_run_all.R
```

Outputs are written to:

```text
RDS3_annotation_visualization_v4.0.1/figures/
```
