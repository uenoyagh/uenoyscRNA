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

## MΦ Clean-B checkpoint (v4.12.1)

Primary macrophage analysis framework was fixed at this checkpoint.

- MΦ-only RPCA clustering: Res2.0
- Manual MΦ annotation: v4.8.4
- Primary cleaned dataset: Clean-B
- Clean-B exclusion:
  - Res2 cluster 24: B-cell contamination
  - Res2 cluster 27: T/NK contamination
- Original / Clean-A / Clean-B / Clean-C sensitivity analysis confirmed that the major Tx-associated increase in Anti-inflammatory-MΦ and M2/M1 balance is robust to contamination cleaning.
- Clean-B is used as the primary MΦ dataset for subsequent subtype characterization.
- Pseudobulk is calculated at the biological-sample level using explicit raw-count aggregation with Matrix::rowSums().
- Anti-inflammatory-MΦ Sham vs Tx mechanistic analysis has been completed through v4.12.1.
- Next phase:
  1. rebuild Clean-B functional-master figures
  2. subtype-specific pseudobulk DE
  3. formal compositional analysis
  4. gene-level DotPlot / heatmap
  5. Anti-inflammatory-MΦ heterogeneity
  6. MΦ output-ligand analysis
  7. subtype-specific pathway enrichment
  8. classification robustness
  9. Whole-liver Res2.0 reconstruction before MΦ-HSC interaction analysis

Checkpoint tag:
`mphi-cleanB-v4.12.1`
