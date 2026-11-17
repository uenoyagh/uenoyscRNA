# Human / Mouse MΦ Final QC and Primary-Clean Checkpoint v1.0

Date: 2026-11-18

## Final QC decision

### Human MΦ
- Frozen clustering resolution: res0.8
- Source cells: 6,720
- Exclude cluster 13: 90 cells
- Primary-clean cells: 6,630
- Retained fraction: 98.66%
- Cluster 13: coherent neutrophil/non-MΦ contamination
- Cluster 3: retain
- Cluster 16: retain as cycling MΦ
- No reclustering

Primary-clean object:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_Mphi_PrimaryClean_v2.0/objects/Human_MASH_Mphi_PrimaryClean_v2.0.rds`

### Mouse MΦ
- Frozen clustering resolution: res2.0
- Source cells: 18,029
- Exclude cluster 23: 173 cells
- Primary-clean cells: 17,856
- Retained fraction: 99.04%
- Cluster 23: MΦ-LSEC heterotypic doublet / mixed-cell candidate
- Cluster 8: retain
- Cluster 15: retain
- Cluster 26: retain as M1/IFN-high MΦ
- No reclustering

Primary-clean object:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mouse_Mphi_PrimaryClean_v5.0/objects/Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`

## Ambient RNA policy

- Original filtered counts remain the primary expression matrix.
- Whole-liver sample-wise DecontX-corrected counts are sensitivity-analysis data only.
- Human cluster 3 and Mouse clusters 8/15 are retained because hepatocyte-derived signatures decrease strongly after DecontX while macrophage programs persist.
- Mouse cluster 23 is excluded from primary MΦ biology because of coherent LSEC/endothelial signal and high-complexity mixed-cell characteristics.
- Mouse cluster 23 remains available for sensitivity analyses.
- Raw/empty-droplet matrices were not available on currently mounted storage.

## Primary QC metadata

- `qc_original_frozen_cluster`
- `qc_primary_include`
- `qc_sensitivity_include`
- `qc_exclusion_reason`
- `qc_final_disposition`
- `qc_audit_version`
- `qc_resolution_policy`
- `qc_ambient_policy`

## Canonical finalization output

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/CrossSpecies_MASH_Mphi/PrimaryClean_Mphi_QC_Finalization_v1.0/`

Key tables:
- `01_final_cluster_disposition_v1.0.tsv`
- `02_cluster_counts_before_after_primary_QC_v1.0.tsv`
- `03_sample_counts_before_after_primary_QC_v1.0.tsv`
- `04_sample_cluster_counts_before_after_primary_QC_v1.0.tsv`
- `05_primary_clean_summary_v1.0.tsv`

## Downstream policy

Primary downstream analyses use:
- Human: `Human_MASH_Mphi_PrimaryClean_v2.0.rds`
- Mouse: `Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`

Priority:
1. M1 and M2 functional axes
2. M2/M1 balance
3. Resident / Recruited / LAM composition
4. Mouse Sham vs Tx
5. Mouse STD vs CDHFD
6. IL10-response program
7. Profibrotic macrophage program
8. MΦ-HSC / MΦ-Hepatocyte / MΦ-LSEC interactions

M1/M2 are functional axes and are not used as clustering identities.

## Reproducibility

- Frozen cluster labels are preserved.
- Human res0.8 and Mouse res2.0 remain fixed.
- QC decisions are independent of desired treatment/disease direction.
- Source RDS files were not overwritten.
