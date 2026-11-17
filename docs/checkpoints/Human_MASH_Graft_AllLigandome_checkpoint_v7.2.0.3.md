# Human MASH Graft All-Ligandome Checkpoint v7.2.0.3

Date: 2026-09-11

## Repository

- Repository: /Users/uenoya/Projects/uenoyscRNA
- Branch: main
- Parent checkpoint: 38150e4
- Analysis script:
  - analysis/Human_MASH_Graft_AllLigandome_Tx17_Tx5_v7.2.0.3.R
- Script SHA256:
  - 117468ca69ff54f02dd1ed0480a50349e00b4a116b049125c6f009d0071c8a15

## Purpose

Ligand-agnostic characterization of the human graft secreted/signaling ligand repertoire in Tx17 and Tx5, followed by strict human-to-mouse orthology and NicheNet-v2 compatibility filtering.

IL10 was not privileged in ranking or filtering.

## Primary inputs

Human filtered raw MEX:

- Tx17:
  /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/RH251217117_raw_MEX/Tx17/filtered_feature_bc_matrix
- Tx5:
  /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/RH251217117_raw_MEX/Tx5/filtered_feature_bc_matrix

Canonical context metadata:

- /Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3.1/analysis/v1.4.3.1_12sample_integrated_context_review/tables/v1.4.3_context_UMAP_metadata.tsv

NicheNet-v2 resources:

- Zenodo record 7074291
- Human LR network
- Mouse LR network
- Mouse ligand-target matrix

## Cell-space QC

Tx17:

- Raw filtered MEX: 18,132 features / 5,959 cells
- Canonical analysis cells: 5,958
- Excluded raw-MEX-only cell:
  ATACCCCAGTCAACTTCAACTCAAAC-1

Tx5:

- Raw filtered MEX: 18,132 features / 8,571 cells
- Canonical analysis cells: 8,571
- Excluded cells: 0

## Ligand universe and bridge results

- NicheNet human LR ligand universe: 1,226
- Measurable in both Tx samples: 1,154
- Detected in either Tx: 599
- Detected in both Tx: 334
- >=3 positive cells in both Tx: 277
- Strict unique human-to-mouse ligand orthologs: 996
- Cross-species NicheNet compatible: 986
- Eligible for mouse receiver bridge: 560
- Eligible for primary replicate-supported bridge: 326

## Final uniqueness validation

Ranked ligandome:

- Rows: 1,226
- Unique human ligands: 1,226
- Duplicated human ligands: 0

Cross-species bridge candidates:

- Rows: 560
- Unique human ligands: 560
- Duplicated human ligands: 0

## Source annotation fix

Earlier v7.2.0.2 incorrectly merged lineage-total and state-total dominant-source rows through the combined dominant_sources table.

This caused a 2 x 2 Cartesian expansion for measurable ligands.

v7.2.0.3 fixes this by independently joining:

- dominant_lineage_total
- dominant_lineage_per1000
- dominant_state_total

and enforcing one human ligand per row after source and orthology joins.

## IL10 audit

IL10 remained in the unbiased ligand universe without preferential weighting.

- Descriptive sender rank: 261
- Tx17: 18 UMI / 13 positive cells
- Tx17 pseudobulk CPM: 2.125722
- Tx5: 4 UMI / 4 positive cells
- Tx5 pseudobulk CPM: 0.4450291
- Mouse ortholog: Il10
- Eligible for mouse receiver bridge: TRUE
- Eligible for primary replicate-supported bridge: TRUE

IL10 source-state interpretation remains descriptive because expression is sparse and source annotation differs between Tx17 and Tx5.

## Legacy-15 audit

14 of the 15 legacy audit genes are members of the NicheNet-v2 human LR ligand universe.

IGFBP2 is not present in the NicheNet-v2 human LR ligand universe and therefore does not appear in the ranked ligandome or NicheNet bridge table.

IGFBP2 should remain available as a NicheNet-external audit candidate where appropriate.

## Output

/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_Graft_AllLigandome_v7.2.0.3

Key outputs:

- 01_Tx_MEX_context_mapping_QC_v7.2.0.3.tsv
- 02_human_NicheNet_ligand_feature_space_v7.2.0.3.tsv
- 09_human_graft_all_ligandome_ranked_v7.2.0.3.tsv
- 10_cross_species_bridge_candidate_ligands_v7.2.0.3.tsv
- 11_legacy15_IL10_audit_only_v7.2.0.3.tsv
- 12_all_ligandome_summary_v7.2.0.3.tsv
- Human_Graft_AllLigandome_v7.2.0.3.rds

## Status

SUCCESS

## Next analysis

v7.2.1:

Mouse receiver receptor competence and cross-species ligand-receptor admissibility for:

- Kupffer/Macrophage
- Monocyte
- HSC/Mesenchymal
- Hepatocyte
- LSEC
- Cholangiocyte

Then v7.2.2:

NicheNet ligand activity against the receiver-specific v7.1.2 CORE / EXPANDED host target sets.
