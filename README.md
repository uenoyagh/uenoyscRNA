# uenoy scRNAseq Framework v4.0 visualization patch

This patch adds two complementary visualization systems:

1. Interpretation UMAPs
   - lineage
   - cell type
   - subtype
   - confidence
   - annotation difference

2. Explanatory summary UMAP
   - one UMAP containing the finest reliable annotation
   - subtype is used only at High/Moderate confidence
   - otherwise falls back to cell type, then lineage
   - two-line labels and lineage-consistent R8-style colors

All square UMAPs use `coord_fixed(ratio = 1)` and are saved at 11 × 11 inches.

## Install

```bash
cd /path/to/uenoyscRNA_v4.0_visualization_patch
bash install_v4.0_visualization.sh
```

## Run

```bash
cd "/Users/uenoya/Projects/uenoyscRNA"
Rscript analysis/RDS3_annotation_visualization_v4.0/99_run_all.R
```

## Inputs

The pipeline expects:

```text
RDS3_annotation_validation_v3.4/
├── objects/RDS3_with_annotation_validation_v3.4.rds
└── tables/annotation_audit_v3.4.csv
```

## Main outputs

```text
RDS3_annotation_visualization_v4.0/
├── figures/
│   ├── UMAP_lineage_R8_square.pdf
│   ├── UMAP_celltype_R8_square.pdf
│   ├── UMAP_subtype_R8_square.pdf
│   ├── UMAP_summary_R8_square.pdf
│   ├── UMAP_summary_R8_explanatory_wide.pdf
│   ├── UMAP_difference_R8_square.pdf
│   ├── UMAP_confidence_lineage_R8_square.pdf
│   ├── UMAP_confidence_celltype_R8_square.pdf
│   ├── UMAP_confidence_subtype_R8_square.pdf
│   └── Cluster_annotation_audit_v4.0.pdf
├── tables/
└── objects/
```

The explanatory summary label is stored in:

```text
vote_ueno_summary_v40
```

The current-versus-recommended review status is stored in:

```text
annotation_difference_v40
```
