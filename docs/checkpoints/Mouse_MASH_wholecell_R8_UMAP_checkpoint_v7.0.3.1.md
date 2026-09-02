# Mouse MASH whole-cell R8 UMAP checkpoint v7.0.3.1

Date: 2026-09-02

## Scope

This checkpoint updates the frozen mouse whole-cell UMAP for pre-human-analysis visualization using the current lineage annotations and R8-style high-saturation colors.

The frozen whole-cell UMAP coordinates are preserved. No reintegration, reclustering, or UMAP recomputation is performed.

## Coordinate source

- Whole-cell frozen parent:
  `/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/WholeCell_Layer1_ParentFreeze_v5.1.1/RDS/Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds`
- Reduction:
  `umapRPCA`
- Total cells:
  `104,588`

## Current refined annotation update

Current frozen lineage-specific annotations were mapped back to the whole-cell object by exact cell-barcode matching.

Lineages included:
- Kupffer/Macrophage
- HSC/Mesenchymal
- Hepatocyte
- LSEC
- Cholangiocyte
- Monocyte

Final ownership logic in v7.0.1:
- Dedicated refined lineage objects may override older whole-cell broad-parent labels by exact barcode.
- Monocyte v6.9.4 has priority over historical Clean-B Mphi for overlapping cells.
- Mphi subtypes are used for broad Kupffer/Macrophage cells not owned by Monocyte v6.9.4.
- Parent-only/unresolved cells remain visible rather than being discarded.

## Annotation-source barcode audit

All six lineage-specific source objects showed 100% exact barcode matching to the whole-cell object.

Source object sizes:
- Mphi: 18,029
- HSC: 13,166
- Hepatocyte: 10,302
- LSEC: 9,190
- Cholangiocyte: 15,755
- Monocyte: 3,490

The large historical Mphi cross-parent component reflects the original Clean-B macrophage compartment that also contained Monocyte cells.

## Biological sample recovery

Sample metadata fields in later whole-cell objects were incomplete for Sham/Tx. Biological sample identity was therefore recovered directly from the exact cell-name prefixes.

Recovered sample counts:
- STD_rep1: 12,211
- CDHFD_rep1: 11,543
- Sham1: 20,658
- Sham20: 17,724
- Tx17: 17,621
- Tx5: 24,831

Recovered condition totals:
- STD: 12,211
- CDAHFD: 11,543
- Sham: 38,382
- Tx: 42,452

All four condition totals exactly matched the validated expected whole-cell totals.

## Sham vs Tx refined UMAP

Successful script:
`analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.4.R`

Outputs include:
- full Sham vs Tx refined UMAP
- balanced Sham vs Tx refined UMAP
- Sham1 / Sham20 / Tx17 / Tx5 refined UMAP
- Sham vs Tx Mphi + Monocyte highlight UMAP

Balanced display uses equal numbers of cells from each biological sample and is for visualization only.

## Sham vs Tx balanced broad R8 UMAP

Successful script:
`analysis/Mouse_MASH_wholecell_Sham_vs_Tx_balanced_broad_R8_v7.0.3.1.R`

The current refined annotation is collapsed back to broad lineages while retaining v7.0.1 lineage ownership.

Balanced sampling:
- Sham1: 17,621
- Sham20: 17,621
- Tx17: 17,621
- Tx5: 17,621

Therefore:
- Balanced Sham total: 35,242 cells
- Balanced Tx total: 35,242 cells

The balanced broad R8 UMAP was successfully generated.

## R8 display policy

Broad R8 display retains the established high-saturation direction:
- LSEC: vivid cyan
- Hepatocyte: blue-green
- Cholangiocyte: vivid green
- HSC/Mesenchymal: vivid pink
- Kupffer/Macrophage, Monocyte, and other immune lineages use strongly separated colors

## Successful scripts included in this checkpoint

- `analysis/Mouse_MASH_wholecell_current_refined_UMAP_R8_v7.0.1.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.4.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_balanced_broad_R8_v7.0.3.1.R`

## Superseded / failed scripts intentionally excluded from checkpoint

Do not stage:
- `analysis/Mouse_MASH_wholecell_current_refined_UMAP_R8_v7.0.0.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.1.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.2.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_refined_UMAP_R8_v7.0.2.3.R`
- `analysis/Mouse_MASH_wholecell_Sham_vs_Tx_balanced_broad_R8_v7.0.3.R`

These files may remain in the working tree for provenance but are not part of the reproducible checkpoint.

## Interpretation / use

This checkpoint is a visualization and annotation-integration checkpoint before starting the human scRNA-seq analysis.

The key deliverables are:
1. updated whole-cell refined R8 UMAP
2. Sham vs Tx refined comparison
3. Sham vs Tx balanced broad R8 comparison
4. biological-replicate UMAPs for Sham1, Sham20, Tx17, and Tx5

The balanced UMAPs are descriptive visualization tools and are not used as inferential statistical tests.

## Recommended Git checkpoint

Commit message:
`Checkpoint mouse whole-cell R8 UMAP update v7.0.3.1`

Tag:
`mouse-mash-wholecell-r8-umap-v7.0.3.1`
