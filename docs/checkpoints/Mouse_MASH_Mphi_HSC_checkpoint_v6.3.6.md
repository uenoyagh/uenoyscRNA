# Mouse MASH scRNA-seq checkpoint — Macrophage–HSC interaction

**Checkpoint date:** 2026-09-02  
**Recommended Git tag:** `mouse-mash-mphi-hsc-mechanistic-closure-v6.3.6`

## Scope

This checkpoint closes the Mouse MASH macrophage–HSC interaction analysis through v6.3.6.

Biological samples:
- Sham1
- Sham20
- Tx17
- Tx5

Biological replication:
- Sham n=2
- Tx n=2

All treatment-level inference remains exploratory because biological replication is n=2 per group.

## Confirmed upstream frozen inputs

Whole-cell parent:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/WholeCell_Layer1_ParentFreeze_v5.1.1/RDS/Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds`

Macrophage Clean-B:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mphi_Res2_CleanB_FINAL_v4.14.5/RDS/Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds`

Interaction-ready Mphi5 + HSC3 object:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Interaction/Mphi5_HSC3_interaction_ready_v6.2.0/RDS/Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds`

## Successful scripts recommended for checkpoint

### Preparation and HSC annotation
- `analysis/Mouse_MASH_Mphi5_HSC_interaction_prepare_v6.0.0.R`
- `analysis/Mouse_MASH_HSC_subclustering_v6.1.0.R`
- `analysis/Mouse_MASH_HSC_refined_annotation_v6.1.2.R`
- `analysis/Mouse_MASH_Mphi5_HSC3_interaction_prepare_v6.2.0.R`

### CellChat
- `analysis/Mouse_MASH_Mphi5_HSC3_CellChat_samplewise_v6.2.1.R`
- `analysis/Mouse_MASH_Mphi5_HSC3_CellChat_refine_v6.2.2.R`

### NicheNet and 3-way integration
- `analysis/Mouse_MASH_Mphi5_HSC3_NicheNet_v6.3.0.R`
- `analysis/Mouse_MASH_Mphi5_HSC3_3way_concordance_v6.3.1.R`

### PDGFB / FN1 validation and mechanistic closure
- `analysis/Mouse_MASH_PDGFB_axis_validation_v6.3.2.R`
- `analysis/Mouse_MASH_PDGFB_population_weighted_v6.3.3.1.R`
- `analysis/Mouse_MASH_PDGFB_independent_validation_v6.3.4.1.R`
- `analysis/Mouse_MASH_PDGF_FN1_rewiring_v6.3.5.R`
- `analysis/Mouse_MASH_PDGF_FN1_interaction_closure_v6.3.6.R`

## Superseded / non-final scripts not recommended for the final checkpoint

- `Mouse_MASH_HSC_refined_annotation_v6.1.1.R`
  - superseded by v6.1.2
- `Mouse_MASH_PDGFB_population_weighted_v6.3.3.R`
  - superseded by v6.3.3.1
- `Mouse_MASH_PDGFB_independent_validation_v6.3.4.R`
  - superseded by v6.3.4.1

Note:
The executed files named `v6.3.3.1` and `v6.3.4.1` contain a harmless header-version typo (`v6.3.3.1.1` / `v6.3.4.1.1`) in the script comment, while their output paths/manifests use the intended versions v6.3.3.1 and v6.3.4.1. Preserve the executed scripts as-is for reproducibility.

## Final HSC states

- qHSC
- ECM-activated HSC
- Contractile/myofibroblastic HSC

A two-level qHSC vs aHSC representation was also retained for supporting analyses.

## Final biological interpretation

### HSC-state remodeling

ECM-activated HSC fraction decreases after transplantation:
- Sham mean: ~0.527
- Tx mean: ~0.452
- Delta: ~-0.075
- 4/4 pairwise comparisons down
- Strong Tx-down evidence

Cycling fraction within ECM-activated HSC decreases:
- Sham mean: ~0.0560
- Tx mean: ~0.0478
- Delta: ~-0.0081
- 3/4 pairwise comparisons down
- Moderate Tx-down evidence

Cycling ECM-HSC burden within the HSC3 compartment decreases:
- Sham mean: ~0.0298
- Tx mean: ~0.0216
- approximately 28% lower
- 3/4 pairwise comparisons down
- Moderate Tx-down evidence

### PDGFB axis

Repair/Resolution-Mphi Pdgfb per-cell expression:
- Sham mean: 5.971 CP10k
- Tx mean: 3.568 CP10k
- Delta: -2.403
- 4/4 pairwise comparisons down
- Strong Tx-down evidence

CellChat PDGF communication:
Repair/Resolution-Mphi -> ECM-activated HSC
- Sham mean: ~0.1846
- Tx mean: ~0.1088
- Delta: ~-0.0758
- 4/4 pairwise comparisons down
- Strong Tx-down evidence

ECM-activated HSC PDGF mitogenic module:
- Sham mean: ~0.479
- Tx mean: ~-0.479
- Delta: ~-0.958
- 4/4 pairwise comparisons down
- Strong Tx-down evidence

Final PDGF interpretation:
`Coherent_Tx_down_mitogenic_axis`

The data support selective attenuation of a PDGF-associated mitogenic/proliferative signaling arm rather than global HSC suppression.

### Important abundance caveat

Repair/Resolution-Mphi abundance increases after Tx.

Therefore:
- per-cell Pdgfb down does not automatically mean Repair/Resolution-Mphi population-weighted PDGFB output decreases;
- per-cell phenotype, subtype abundance-weighted output, and total macrophage output must be interpreted separately.

### FN1 axis

Total Mphi5 population-weighted Fn1:
- Sham mean: 3.444
- Tx mean: 4.674
- Delta: +1.231
- 4/4 pairwise comparisons up
- Strong Tx-up evidence

CellChat FN1 communication:
all Mphi5 -> ECM-activated HSC
- Sham mean: ~0.108
- Tx mean: ~0.401
- Delta: ~+0.294
- 4/4 pairwise comparisons up
- Strong Tx-up evidence

Receiver-side FN1/remodeling transcription:
- mixed / not reproducibly changed

Final FN1 interpretation:
`Sender_and_communication_Tx_up_but_receiver_transcription_mixed`

FN1 therefore represents communication rewiring rather than a demonstrated downstream remodeling mechanism.

## Final working model

Transplantation does not cause uniform suppression of HSC biology.

Instead, it is associated with:
1. reduced ECM-activated HSC burden,
2. reduced ECM-HSC proliferative/cycling activity,
3. a coherent decrease in Repair/Resolution-Mphi PDGFB signaling toward ECM-HSC,
4. increased FN1 communication without consistent downstream FN1/remodeling transcription.

The strongest mechanistically coherent axis is PDGF.

## Main output directories

Interaction root:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Interaction/`

Important milestones:
- Mphi5/HSC preparation: v6.0.0
- HSC final annotation: v6.1.2
- Mphi5/HSC3 interaction-ready object: v6.2.0
- CellChat: v6.2.1 / v6.2.2
- NicheNet: v6.3.0
- 3-way concordance: v6.3.1
- PDGFB focused validation: v6.3.2–v6.3.4.1
- PDGF/FN1 rewiring: v6.3.5
- final interaction closure: v6.3.6

## Repository policy

- R scripts, reusable framework, checkpoint documentation: GitHub repository on Mac internal drive
- RDS, PDF, CSV, Excel and other analysis outputs: external SSD
- Use explicit `git add` only
- Never use `git add .`
