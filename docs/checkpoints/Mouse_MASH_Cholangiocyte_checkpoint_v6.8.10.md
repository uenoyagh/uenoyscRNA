# Mouse MASH Cholangiocyte mechanistic closure checkpoint v6.8.10

## Checkpoint

Final mouse Cholangiocyte mechanistic closure:

- final version: v6.8.10
- Git tag target:
  `mouse-mash-cholangiocyte-mechanistic-closure-v6.8.10`

## Input

Frozen whole-cell source:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/RDS3_with_visualization_metadata_v4.1.1.rds`

Whole-cell annotation:

`celltype_for_R8plot_FIXED2`

Samples:

- STD_rep1
- CDHFD_rep1
- Sham1
- Sham20
- Tx17
- Tx5

Conditions:

- STD
- CDHFD
- Sham
- Tx

## Lineage cleanup

Original annotated Cholangiocytes:

- 17,000 cells

Cycling-derived biliary rescue audit:

- permissive candidates: 700
- stringent rescue candidates: 3
- decision: do not add Cycling-derived cells

Initial res0.3 lineage audit identified clear non-biliary contamination.

Final excluded pre-cleanup clusters:

- cluster 3: myeloid contamination
- cluster 9: vascular/mesenchymal contamination
- cluster 10: neutrophil contamination
- cluster 11: dendritic-cell contamination
- cluster 12: B-cell contamination

Cluster 4 was initially removed in v6.8.3 but subsequently restored after sample-level QC review.

Final lineage-clean baseline:

- 15,755 cells
- 14,734 features

v6.8.3 is superseded by v6.8.3.2 and is not part of this checkpoint.

## Cluster 4 QC-watch

Pre-cleanup cluster 4:

- 420 cells
- enriched in STD/CDHFD
- unusually high RNA-count complexity in STD/CDHFD
- biliary core retained across samples
- multi-lineage signal required caution

After clean RPCA, most pre-cleanup cluster-4 cells reassembled into a discrete population.

At res0.4:

- resulting cluster 4: 565 cells
- 329/565 = 58.2% originated from pre-cleanup cluster 4

Final decision:

- retain the population
- label as `Disease_enriched_high_complexity_QC_watch`
- do not interpret as a definitive disease-specific state
- carry QC-watch exclusion as sensitivity analysis

## Clean RPCA and frozen annotation

Clean RPCA:

- dimensions: 1:20
- final annotation scaffold: res0.4
- 11 frozen states

Frozen fine states:

1. `Pcp4l1_Serpina6_homeostatic_like`
2. `Kptn_Pacs2_homeostatic_like`
3. `Ccl2_Vcam1_inflammatory_reactive`
4. `Msln_Aqp5_ductular_like`
5. `Disease_enriched_high_complexity_QC_watch`
6. `Cycling_cholangiocyte`
7. `Ciliated_cholangiocyte`
8. `IEG_stress_response`
9. `Krt20_Cdh17_reactive_epithelial`
10. `Dmbt1_Duox2_reactive_epithelial`
11. `Tuft_like_cholangiocyte`

Important limitations:

- `Dmbt1_Duox2_reactive_epithelial` is Sham20-biased and exploratory.
- `Tuft_like_cholangiocyte` contains only 15 cells and is rare exploratory.
- STD/CDHFD contain one biological sample each and are descriptive only.

## State-composition and module results

QC-watch exclusion did not materially alter the Sham-to-Tx conclusions.

Tx-associated state changes included:

- increased `Krt20_Cdh17_reactive_epithelial`
- increased `Ccl2_Vcam1_inflammatory_reactive`
- increased `Msln_Aqp5_ductular_like`
- reduced Cycling
- reduced Ciliated
- reduced `Dmbt1_Duox2_reactive_epithelial`, although this state is Sham20-biased

Module-level Tx effects included increases in:

- Stress/IEG
- Biliary identity
- Inflammatory/adhesion
- Krt20/Cdh17-reactive program

The major trends persisted after QC-watch exclusion.

## Biological-sample pseudobulk

Sham vs Tx:

- Sham1, Sham20
- Tx17, Tx5
- n=2 biological samples/group
- edgeR quasi-likelihood pseudobulk

Whole-cell and QC-watch-excluded analyses showed limited FDR-level DE.

PRIMARY_CORE identified five genes at FDR < 0.10 before reference-specificity filtering:

- Slc25a47
- Gdf15
- Ccn1
- Jun
- Actg1

## Ambient/reference specificity audit

Whole-cell reference lineages were used to identify genes dominated by:

- Hepatocyte
- Kupffer/Macrophage
- LSEC
- HSC/Mesenchymal
- Monocyte
- Neutrophil

Primary interpretation excludes Hepatocyte-dominant genes.

Among PRIMARY_CORE FDR < 0.10 genes:

Excluded as Hepatocyte-dominant:

- Slc25a47
- Gdf15

Retained as ambient-aware Cholangiocyte-associated Tx-up candidates:

- Ccn1
- Jun
- Actg1

All three were Tx-up in both Tx biological samples relative to both Sham biological samples.

## Ambient-aware Hallmark GSEA

Primary analysis:

- Hepatocyte-dominant genes excluded

Strict sensitivity:

- Hepatocyte-dominant genes excluded
- genes strongly dominated by any major reference lineage also excluded

Major Tx-enriched pathways persisted under strict reference filtering:

- TNFA signaling via NF-kB
- Hypoxia
- Epithelial-mesenchymal transition
- p53 pathway
- Apoptosis
- Inflammatory response
- TGF-beta signaling
- Apical junction
- Bile acid metabolism
- Oxidative phosphorylation

Major Tx-depleted proliferative pathways:

- E2F targets
- G2M checkpoint

The major pathway directions were concordant across:

- ALL
- NO_QCWATCH
- PRIMARY_CORE

## State-specific validation

State-specific pseudobulk GSEA showed that the global Tx response was not explained only by redistribution among Cholangiocyte states.

In particular, `Homeostatic_combined` Cholangiocytes showed Tx-associated increases in:

- TNFA/NF-kB
- Hypoxia
- EMT/remodeling
- p53
- Apoptosis
- TGF-beta signaling

and decreases in proliferative E2F/G2M programs.

TNFA/NF-kB and hypoxia responses were recurrent across multiple Cholangiocyte states.

## Mechanistic closure

The data do not support a simple model of Cholangiocyte normalization after transplantation.

The most appropriate interpretation is:

**Tx is associated with broad adaptive/reactive Cholangiocyte epithelial reprogramming.**

This includes:

- preserved or modestly increased biliary identity
- strong NF-kB/IEG response
- hypoxic/stress response
- p53/apoptotic response
- EMT/remodeling programs
- reduced proliferative E2F/G2M programs

These effects remain after:

- QC-watch exclusion
- Hepatocyte ambient-RNA filtering
- strict multi-lineage reference sensitivity
- state-specific validation

## Interpretation limits

- Sham vs Tx: n=2 biological samples/group.
- STD vs CDHFD: n=1/group; descriptive only.
- GSEA FDR values do not overcome the small biological-sample number.
- Tx-associated remodeling cannot yet be classified as beneficial or detrimental.
- Direct graft-to-Cholangiocyte signaling is not established here.
- Graft-derived ligand → host Cholangiocyte receptor mechanisms should be tested later in the planned cross-species analysis.

## Successful scripts included in this checkpoint

- `analysis/Mouse_MASH_Cholangiocyte_parent_audit_v6.8.0.R`
- `analysis/Mouse_MASH_Cholangiocyte_cycling_rescue_validation_v6.8.0.1.R`
- `analysis/Mouse_MASH_Cholangiocyte_clean_RPCA_resolution_scan_v6.8.1.R`
- `analysis/Mouse_MASH_Cholangiocyte_res0.3_expanded_lineage_audit_v6.8.2.R`
- `analysis/Mouse_MASH_Cholangiocyte_removed_cluster_sample_validation_v6.8.3.1.R`
- `analysis/Mouse_MASH_Cholangiocyte_lineage_cleanup_revised_v6.8.3.2.R`
- `analysis/Mouse_MASH_Cholangiocyte_cluster4_high_complexity_QC_v6.8.3.3.R`
- `analysis/Mouse_MASH_Cholangiocyte_clean_RPCA_resolution_scan_v6.8.4.R`
- `analysis/Mouse_MASH_Cholangiocyte_res0.4_state_audit_v6.8.5.R`
- `analysis/Mouse_MASH_Cholangiocyte_annotation_freeze_v6.8.6.R`
- `analysis/Mouse_MASH_Cholangiocyte_state_composition_module_axis_v6.8.7.R`
- `analysis/Mouse_MASH_Cholangiocyte_pseudobulk_Sham_vs_Tx_v6.8.8.R`
- `analysis/Mouse_MASH_Cholangiocyte_reference_specificity_audit_v6.8.8.1.R`
- `analysis/Mouse_MASH_Cholangiocyte_ambient_aware_Hallmark_GSEA_v6.8.9.R`
- `analysis/Mouse_MASH_Cholangiocyte_state_specific_Hallmark_GSEA_v6.8.9.1.R`
- `analysis/Mouse_MASH_Cholangiocyte_mechanistic_closure_v6.8.10.R`

## Explicitly excluded from checkpoint

Superseded:

- `analysis/Mouse_MASH_Cholangiocyte_lineage_cleanup_v6.8.3.R`

Development backups:

- `analysis/Mouse_MASH_Cholangiocyte_res0.3_expanded_lineage_audit_v6.8.2.R.bak`
- `analysis/Mouse_MASH_Cholangiocyte_lineage_cleanup_revised_v6.8.3.2.R.bak`
- `analysis/Mouse_MASH_Cholangiocyte_reference_specificity_audit_v6.8.8.1.R.bak`

All unrelated untracked analysis files remain outside this checkpoint.
