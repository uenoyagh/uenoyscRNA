# Mouse MASH Monocyte checkpoint v6.9.10

## Status
Mechanistic closure completed successfully for the Mouse MASH Monocyte analysis.

Final closure script:
- `analysis/Mouse_MASH_Monocyte_mechanistic_closure_v6.9.10.R`

Recommended Git tag:
- `mouse-mash-monocyte-mechanistic-closure-v6.9.10`

## Source object
Frozen whole-cell reference:
- `/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/RDS3_with_visualization_metadata_v4.1.1.rds`

Frozen Monocyte analysis object:
- `/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.4/objects/Mouse_MASH_Monocyte_annotation_frozen_v6.9.4.rds`

State/module-scored object:
- `/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.5/objects/Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds`

## Cell population
Frozen whole-cell annotation contained 3,556 Monocyte cells.

Lineage cleanup:
- Removed only RPCA res0.4 clusters 7 and 8 from v6.9.1, corresponding to 66 lymphoid/B-cell contaminants.
- Final lineage-clean Monocyte population: 3,490 cells.
- Hepatocyte-RNA-high Monocyte cluster was retained as QC-watch and excluded by sensitivity frameworks rather than removed.
- Disease-enriched high-complexity inflammatory Monocytes were retained because the signal was not consistent with a bulk doublet artifact.

## Frozen Monocyte states
Final RPCA scaffold: resolution 0.4.

1. `S100a8_S100a9_Thbs1_stress_inflammatory_Monocyte`
2. `Mmp8_Sell_Chil3_Vcan_classical_inflammatory_Monocyte`
3. `Pald1_C3ar1_homeostatic_like_Monocyte`
4. `Tnf_Il1rn_Olr1_Gpnmb_inflammatory_remodeling_Monocyte`
5. `Cd300e_Pglyrp1_Cd36_S1pr5_activated_Monocyte`
6. `Adamdec1_Pecam1_low_complexity_state` — sample-biased exploratory
7. `Ms4a7_Mmp12_Dab2_C1q_monocyte_to_macrophage_transition`
8. `Nos2_Cxcl9_Saa3_IFNg_inflammatory_Monocyte` — disease-enriched primary
9. `Hepatocyte_RNA_high_Monocyte_QC_watch` — QC-watch sensitivity
10. `Ifit_Rsad2_Cmpk2_IFN_responsive_Monocyte`

## Composition result
For Sham vs Tx, biological n=2/group.

The principal composition signal was a modest increase in the
`Cd300e_Pglyrp1_Cd36_S1pr5_activated_Monocyte` state in Tx.

There was no evidence for a broad anti-inflammatory normalization of the Monocyte compartment.

## Reference-specificity framework
A genome-wide competing-lineage audit was performed before pathway analysis.

Primary supported competing references required at least 20 cells in each Sham/Tx biological sample.

Supported primary references:
- Hepatocyte
- Kupffer_Macrophage
- HSC_Mesenchymal
- LSEC
- Cholangiocyte
- Neutrophil
- Dendritic_cell
- NK_cell
- Vascular_endothelial
- Plasma_cell
- Mesothelial

Sparse supplementary references excluded from determining the primary maximum reference delta:
- B_cell
- T_cell
- RBC
- Platelet

Primary reference-aware criterion:
- maximum supported competing-lineage delta <= 2 log2CPM

Strict sensitivity criterion:
- maximum supported competing-lineage delta <= 1.5 log2CPM

PRIMARY_CORE:
- tested genes: 10,227
- reference-aware retained: 7,994
- strict retained: 7,111

Important gene-level interpretation:
- `Gdf15`: reference-lineage dominant (Hepatocyte); excluded from primary Monocyte interpretation.
- `Thbs1`: retained in reference-aware analysis but excluded in strict sensitivity; shared/borderline rather than Monocyte-specific.
- `Slc25a47`, `Cyp2c29`, `Krt8`, `Gstm3`, `Ccn1`, `Rgs1`, `Ighm`: excluded from primary Monocyte interpretation because of competing-lineage dominance.
- `Cxcr4`, `Ddit4`, `Zfp36`, `Cebpb`, `Cxcl2`: retained by reference-aware filtering.

## Global Hallmark GSEA
Primary analysis:
- PRIMARY_CORE
- support-aware reference filter
- rank = sign(logFC) * sqrt(edgeR QL F)

Core Tx-enriched pathways:
- TNFA signaling via NF-kB: NES 2.4696, FDR 1.14e-11
- Hypoxia: NES 2.0279, FDR 4.74e-05
- p53 pathway: NES 1.8234, FDR 3.12e-04

All three remained significant under the strict reference filter:
- TNFA/NF-kB: NES 2.2829, FDR 1.98e-07
- Hypoxia: NES 2.0607, FDR 5.50e-05
- p53: NES 1.7313, FDR 2.13e-03

Supporting Tx-associated pathways included:
- UV/stress response
- IL2/STAT5 signaling
- KRAS signaling up
- apoptosis
- oxidative phosphorylation

Secondary/reference-filter-sensitive pathways included:
- TGF-beta signaling
- EMT-like/remodeling Hallmark
- adipogenesis
- fatty-acid/metabolic programs

`EMT` is not interpreted as literal epithelial-to-mesenchymal transition in Monocytes; it is treated as a shared remodeling/cytoskeletal/ECM-associated Hallmark program.

## State-specific Hallmark GSEA
Eligibility:
- primary analysis class
- >=10 cells in each of Sham1, Sham20, Tx17, Tx5

Eligible frozen states:
- clusters 0, 1, 2, 3, 4, 6

Core within-state result:
- TNFA/NF-kB: Tx-enriched in 6/6 states
  - reference-aware FDR <0.10 in 5/6
  - strict FDR <0.10 in 4/6
- Hypoxia: Tx-enriched in 6/6 states
  - reference-aware FDR <0.10 in 4/6
  - strict FDR <0.10 in 4/6
- p53/stress: Tx-enriched in 6/6 states
  - reference-aware FDR <0.10 in 3/6
  - strict FDR <0.10 in 3/6

Thus the Tx signal is not explained only by state-composition changes. It includes broad within-state transcriptional reprogramming.

## Recurrent leading-edge network
Core robust recurrent genes across the TNFA/NF-kB, Hypoxia, and p53/stress pathways:

Across all three core pathways:
- `Fos`
- `Atf3`

Across two core pathways:
- `Dusp1`
- `Zfp36`
- `Nfil3`
- `Sat1`
- `Ddit4`
- `Btg2`
- `Foxo3`
- `Klf6`

Additional recurrent single-core-pathway genes included:
- `Cxcr4`
- `Per1`
- `Sik1`
- `Cebpb`
- `Id2`
- `Junb`
- `Cebpa`
- `Pnrc1`
- `Plin2`
- `Nfkbia`
- `Tank`
- `Bhlhe40`
- `Zfp292`
- `Dennd5a`
- `Slc3a2`
- `Sgk1`
- `Abca1`
- `Tcn2`
- `Cdkn1b`
- `Ldha`

These genes are interpreted as a recurrent pathway-level network, not as independently significant causal driver genes. Most individual gene-level FDR values are not significant with n=2/group.

## Final mechanistic interpretation
Tx induces modest compositional redistribution together with broad within-state TNFA/NF-kB–Hypoxia–p53/stress adaptive reprogramming of the Monocyte compartment, rather than simple anti-inflammatory normalization.

The recurrent leading-edge network is consistent with integrated immediate-early, inflammatory-feedback, hypoxic, and cellular stress-adaptation responses.

In the broader myeloid context, this Monocyte response differs from the mature macrophage compartment, where Tx is associated with anti-inflammatory and repair/resolution remodeling. The combined data support heterogeneous Tx-associated myeloid reorganization rather than uniform anti-inflammatory conversion.

## Not established
The current analyses do not establish:
- direct Monocyte-to-macrophage differentiation
- causal ligand-receptor signaling
- causal activity of any individual leading-edge gene
- a beneficial or harmful functional consequence of the Monocyte reactive program by itself

## Statistical limitations
- Sham vs Tx biological n=2/group.
- STD vs CDHFD biological n=1/group and is descriptive only.
- State-specific analyses require adequate cells in every biological replicate.
- Gene-level FDR is weak for most individual genes; pathway-level and recurrent-within-state evidence is emphasized.

## Successful scripts for checkpoint
Stage only these Monocyte scripts:

- `analysis/Mouse_MASH_Monocyte_parent_audit_v6.9.0.R`
- `analysis/Mouse_MASH_Monocyte_stringent_lineage_audit_v6.9.0.1.R`
- `analysis/Mouse_MASH_Monocyte_clean_RPCA_resolution_scan_v6.9.1.R`
- `analysis/Mouse_MASH_Monocyte_lineage_cleanup_v6.9.2.R`
- `analysis/Mouse_MASH_Monocyte_clean_RPCA_resolution_scan_v6.9.3.R`
- `analysis/Mouse_MASH_Monocyte_targeted_contamination_audit_v6.9.3.1.R`
- `analysis/Mouse_MASH_Monocyte_cluster7_8_QC_audit_v6.9.3.2.R`
- `analysis/Mouse_MASH_Monocyte_annotation_freeze_v6.9.4.R`
- `analysis/Mouse_MASH_Monocyte_state_composition_module_axis_v6.9.5.R`
- `analysis/Mouse_MASH_Monocyte_pseudobulk_Sham_vs_Tx_v6.9.6.R`
- `analysis/Mouse_MASH_Monocyte_reference_specificity_audit_v6.9.6.1.R`
- `analysis/Mouse_MASH_Monocyte_genomewide_reference_specificity_v6.9.6.2.R`
- `analysis/Mouse_MASH_Monocyte_comprehensive_reference_specificity_v6.9.6.3.R`
- `analysis/Mouse_MASH_Monocyte_support_aware_reference_specificity_v6.9.6.4.R`
- `analysis/Mouse_MASH_Monocyte_Hallmark_GSEA_v6.9.7.2.R`
- `analysis/Mouse_MASH_Monocyte_state_specific_Hallmark_GSEA_v6.9.8.R`
- `analysis/Mouse_MASH_Monocyte_recurrent_leading_edge_audit_v6.9.9.1.R`
- `analysis/Mouse_MASH_Monocyte_mechanistic_closure_v6.9.10.R`

Do NOT stage failed/superseded scripts:
- `analysis/Mouse_MASH_Monocyte_Hallmark_GSEA_v6.9.7.R`
- `analysis/Mouse_MASH_Monocyte_Hallmark_GSEA_v6.9.7.1.R`
- `analysis/Mouse_MASH_Monocyte_recurrent_leading_edge_audit_v6.9.9.R`

Do not stage `.bak` files.

## Git checkpoint target
Recommended commit message:
- `Mouse MASH Monocyte mechanistic closure v6.9.10`

Recommended annotated tag:
- `mouse-mash-monocyte-mechanistic-closure-v6.9.10`
