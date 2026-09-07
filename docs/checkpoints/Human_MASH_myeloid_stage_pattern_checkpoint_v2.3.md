# Human MASH myeloid stage-pattern checkpoint v2.3

Date: 2026-09-07  
Repository: `/Users/uenoya/Projects/uenoyscRNA`  
Branch: `main`  
Checkpoint tag target: `human-mash-myeloid-stage-pattern-v2.3`

## Scope

This checkpoint freezes the successful Human MASH LT1–LT8 analysis workflow through the stage-pattern audit v2.3, immediately before ranked pathway/GSEA analysis.

The checkpoint includes R scripts only plus this checkpoint note. RDS files and analysis outputs remain on the external SSD and are not committed to Git.

## Clinical design

| Sample | NAS | Fibrosis |
|---|---:|---:|
| LT7 | 0 | F0 |
| LT8 | 0 | F0 |
| LT5 | 2 | F0 |
| LT6 | 2 | F0 |
| LT3 | 4 | F0 |
| LT4 | 4 | F1 |
| LT1 | 6 | F4 |
| LT2 | 6 | F3 |

Critical interpretation:
- NAS and fibrosis are strongly confounded in advanced disease.
- The NAS0→2→4→6 continuous axis is interpreted as an overall disease-severity axis, not a pure NAS-specific effect.
- NAS6/F3–F4 vs NAS0/F0 is an advanced-disease contrast and cannot separate activity from fibrosis.
- NAS2/F0 vs NAS0/F0 is the cleanest fibrosis-matched early activity contrast in the current dataset.
- Additional F3/F4 samples will improve reproducibility and power; samples that break the NAS–fibrosis correlation are particularly valuable for disentangling fibrosis from activity.

## Current Human MASH workflow status

1. RDS inventory / metadata audit completed.
2. Reference object and feature provenance audited.
3. Canonical LT1–LT8 object built.
4. Broad whole-liver annotation audited.
5. Ambiguous immune compartment reopened and sample-aware RPCA reclustered.
6. Mixed/non-myeloid immune populations audited and removed from the clean myeloid entry set.
7. Clean myeloid core built.
8. Myeloid-specific sample-aware RPCA performed.
9. Working resolution selected at res0.8 after res0.6 vs res0.8 comparison.
10. Borderline clusters audited.
11. Final clean Human MΦ/Monocyte cell set fixed at 6,720 cells.
12. Final broad myeloid states frozen:
   - Resident/Kupffer-like_Mphi
   - Repair/Resolution_Mphi
   - Lipid-associated/TREM2_Mphi
   - Inflammatory/Transitional_Mphi
   - Inflammatory_Monocyte
13. NAS-stage composition quantified at biological-sample level.
14. State-specific sample-level pseudobulk QC completed.
15. Fibrosis-aware pseudobulk DE completed with edgeR QL models.
16. Stage-pattern audit completed to classify early, intermediate, advanced, monotonic, and sample-sensitive trajectories.
17. Ranked pathway/GSEA v2.4 has been prepared but is not included in this checkpoint because successful execution has not yet been confirmed.

## Final clean myeloid object

Final annotated object:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9/objects/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds`

Final clean cell count: 6,720  
RNA features: 16,424  
Frozen numeric clustering resolution: res0.8

Final broad state counts:
- Resident/Kupffer-like_Mphi: 2,692
- Repair/Resolution_Mphi: 1,662
- Inflammatory/Transitional_Mphi: 1,296
- Inflammatory_Monocyte: 569
- Lipid-associated/TREM2_Mphi: 501

## Key composition findings

Across NAS0 → NAS2 → NAS4 → NAS6:
- Lipid-associated/TREM2_Mphi shows the clearest increase toward later disease.
- Repair/Resolution_Mphi is highest at NAS0 and lower in disease groups.
- Resident/Kupffer-like_Mphi and Inflammatory/Transitional_Mphi are non-monotonic and show notable inter-sample heterogeneity.
- Total final-myeloid fraction among captured whole-liver cells does not show a strong monotonic NAS relationship.

All composition inferences are sample-level/descriptive because there are only two biological samples per NAS stage.

## Pseudobulk DE design

Statistical unit: biological sample only.

Primary:
- overall disease-severity trend across all 8 samples
- numeric effect reported per +2 NAS points
- interpreted as disease severity because advanced NAS6 samples are also F3/F4

Secondary:
- advanced disease: NAS6/F3–F4 vs NAS0/F0
- early fibrosis-matched activity: NAS2/F0 vs NAS0/F0

State analysis tiers:
- PRIMARY robust:
  - Resident/Kupffer-like_Mphi
  - Repair/Resolution_Mphi
  - Inflammatory/Transitional_Mphi
- EXPLORATORY lower abundance:
  - Lipid-associated/TREM2_Mphi
  - Inflammatory_Monocyte

## v2.2.1 gene-level results

Overall disease-severity trend, FDR < 0.05:
- Repair/Resolution_Mphi:
  - CHI3L1
- Lipid-associated/TREM2_Mphi:
  - CHI3L1
  - TM4SF19

No FDR < 0.05 genes were detected in the advanced NAS6/F3–F4 vs NAS0/F0 endpoint contrast or the early NAS2/F0 vs NAS0/F0 contrast.

## v2.3 stage-pattern findings

Tier1 FDR < 0.05:
- Repair/Resolution_Mphi / CHI3L1
  - `Early_F0_associated_sustained`
- Lipid-associated/TREM2_Mphi / CHI3L1
  - `Nonmonotonic_or_sample_sensitive`
- Lipid-associated/TREM2_Mphi / TM4SF19
  - `Gradual_monotonic_severity_shift`

Important mechanistic candidates retained for downstream review include:
- Repair/Resolution_Mphi:
  - CHI3L1, THBS1, SPP1
- Lipid-associated/TREM2_Mphi:
  - CHI3L1, TM4SF19, SPP1
- Inflammatory/Transitional_Mphi:
  - IL1R2, HMOX1, DDIT4
- Resident/Kupffer-like_Mphi:
  - FBN2, RNASE1, SIGLEC5
- Inflammatory_Monocyte:
  - DDIT4 and additional exploratory candidates

These patterns are descriptive and are not interpreted as pure fibrosis or pure NAS effects.

## Scripts included in this checkpoint

All 20 scripts below passed syntax parsing on 2026-09-07:

1. `analysis/Human_MASH_RDS_inventory_metadata_audit_v1.0.1.R`
2. `analysis/Human_MASH_reference_object_targeted_audit_v1.1.R`
3. `analysis/Human_MASH_LT_cellmap_feature_provenance_audit_v1.2.R`
4. `analysis/Human_MASH_LT1to8_build_canonical_object_v1.3.R`
5. `analysis/Human_MASH_LT1to8_broad_celltype_validation_v1.4.R`
6. `analysis/Human_MASH_LT1to8_targeted_broad_annotation_audit_v1.4.1.R`
7. `analysis/Human_MASH_LT1to8_immune_candidate_RPCA_recluster_v1.5.R`
8. `analysis/Human_MASH_LT1to8_immune_res06_identity_audit_v1.5.1.R`
9. `analysis/Human_MASH_LT1to8_immune_mixedcluster_audit_v1.5.2.R`
10. `analysis/Human_MASH_LT1to8_build_clean_myeloid_object_v1.6.R`
11. `analysis/Human_MASH_LT1to8_myeloid_RPCA_recluster_v1.7.R`
12. `analysis/Human_MASH_LT1to8_myeloid_resolution_compare_v1.7.1.R`
13. `analysis/Human_MASH_LT1to8_myeloid_res08_state_audit_v1.7.2.R`
14. `analysis/Human_MASH_LT1to8_myeloid_borderline_cluster_audit_v1.7.3.R`
15. `analysis/Human_MASH_LT1to8_build_final_myeloid_set_v1.8.R`
16. `analysis/Human_MASH_LT1to8_final_myeloid_state_annotation_v1.9.R`
17. `analysis/Human_MASH_LT1to8_final_myeloid_NAS_composition_v2.0.R`
18. `analysis/Human_MASH_LT1to8_final_myeloid_pseudobulk_QC_v2.1.R`
19. `analysis/Human_MASH_LT1to8_final_myeloid_pseudobulk_DE_v2.2.1.R`
20. `analysis/Human_MASH_LT1to8_final_myeloid_stage_pattern_audit_v2.3.R`

## Explicitly excluded from this checkpoint

- `analysis/Human_MASH_RDS_inventory_metadata_audit_v1.0.R`
  - failed/superseded (`rbind` column mismatch)
- `analysis/Human_MASH_LT1to8_final_myeloid_pseudobulk_DE_v2.2.R`
  - failed/superseded (`edgeR::makeContrasts` namespace error)
- `analysis/Human_MASH_LT1to8_final_myeloid_ranked_pathway_GSEA_v2.4.R`
  - prepared but not yet successfully executed
- all unrelated untracked Mouse MASH scripts/directories
- all `.bak` files
- all RDS/PDF/CSV/analysis outputs on the external SSD

## Next step after checkpoint

Run and validate:
`analysis/Human_MASH_LT1to8_final_myeloid_ranked_pathway_GSEA_v2.4.R`

Then review:
- FDR < 0.05 / FDR < 0.10 ranked GSEA pathways
- cross-test concordance across:
  - overall disease severity
  - advanced NAS6/F3–F4 vs NAS0/F0
  - early NAS2/F0 vs NAS0/F0

The goal is to separate:
1. early MASH-associated myeloid programs,
2. advanced F3/F4-enriched programs,
3. signals that should remain provisional until fibrosis-balanced samples are added.
