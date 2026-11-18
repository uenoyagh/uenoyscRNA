# Human MASH Mphi–HSC Mechanistic Closure Checkpoint v2.2.0
Date: 2026-11-18

## Scope
Human MASH scRNA-seq analysis of macrophage-to-HSC signaling using the canonical HSC receiver and final myeloid annotation.

## Canonical HSC receiver
- Total canonical HSC cells: 7,640
- Primary HSC labels:
  - HSC_quiescent
  - HSC_fibroblast_mesenchymal
- Sensitivity-only generic Fibroblast_mesenchymal cells were not included in the primary receiver.
- ECM-high is a functional state and must not be treated as an independent HSC subtype.

## Key HSC result
- qHSC: 4,335
- activated/fibroblast-mesenchymal HSC: 3,305
- Primary classification is lineage/state based; functional ECM activation is analyzed separately.

## Human Mphi -> HSC closure summary

### SPP1
- Primary mechanistic axis.
- NicheNet primary rank: 4.
- TREM2-Mphi SPP1 burden vs whole-HSC TGF-beta Hallmark:
  - Spearman rho = 0.7857
  - exact permutation P = 0.02793
- HSC receptor support includes ITGA9, ITGB1, ITGAV, ITGB5, ITGA4, ITGA5.
- NicheNet target support:
  - 19 targets tested
  - 19/19 positive rho
  - 6 exact P < 0.05
  - 1 target survives target-level FDR < 0.05
- Interpretation:
  - Best-supported direct mechanistic axis.
  - Not strictly TREM2-specific; inflammatory monocyte/transitional states also carry SPP1-associated signal.

### CDCP1
- Strongest novel sample-level TREM2-associated candidate.
- TREM2-Mphi CDCP1 burden vs whole-HSC TGF-beta Hallmark:
  - rho = 0.9524
  - exact P = 0.001141
  - FDR within endpoint = 0.003423
  - NAS-adjusted partial rho = 0.9483
  - partial P = 0.001135
- Leave-one-out:
  - all 8 positive
  - minimum rho = 0.9286
  - median rho = 0.9464
- Sender-state specificity:
  - TREM2 state rank = 1
  - exact P = 0.001141
  - FDR across states = 0.005704
- HSC EGFR is measurable in all 8 samples.
- Interpretation:
  - Very strong novel association.
  - Do not yet call CDCP1 a proven macrophage-to-HSC ligand axis.

### TGFB1
- Robust TREM2-state-associated signal.
- TREM2-Mphi TGFB1 burden vs whole-HSC TGF-beta Hallmark:
  - rho = 0.8810
  - exact P = 0.007242
  - FDR within endpoint = 0.01086
  - NAS-adjusted partial rho = 0.8285
  - partial P = 0.02128
- Leave-one-out:
  - all 8 positive
  - minimum rho = 0.8214
  - median rho = 0.8571
- Sender-state specificity:
  - TREM2 state rank = 1
  - exact P = 0.007242
  - FDR across states = 0.03621
- NicheNet primary rank = 206/209.
- Interpretation:
  - Strong state-associated TGF-beta signal.
  - Not supported as the dominant direct TREM2-Mphi ligand driver.

## Current mechanistic hierarchy
1. SPP1 -> HSC integrins: primary mechanistic axis.
2. CDCP1-associated signaling: strongest novel TREM2-specific candidate.
3. TGFB1/TGF-beta response: robust TREM2-state-associated auxiliary axis, but weak NicheNet ligand prioritization.

## Statistical principles
- Biological sample is the unit of inference.
- n = 8 human samples.
- No cell-level pseudoreplication for sample-level associations.
- Exact 8! permutation tests used where applicable.
- NAS-adjusted partial correlations treated as sensitivity analyses.
- No reclustering and no source RDS overwrite.

## Successful scripts to retain
- Human_MASH_HSC_ECM_Proliferation_v1.2.0.R
- Human_MASH_TREM2_TGFB1_TGFBR1_Closure_v1.4.0.R
- Human_MASH_TGFB1_Sender_Burden_Audit_v1.5.0.R
- Human_MASH_TREM2_TGFB1_Robustness_v1.6.0.R
- Human_MASH_Orthogonal_TGFbeta_Response_v1.7.0.R
- Human_MASH_Broad_Orthogonal_TGFbeta_v1.8.0.R
- Human_MASH_TGFbeta_Cell_PB_Concordance_v1.9.0.R
- Human_MASH_TREM2_HSC_NicheNet_v2.0.1.R
- Human_MASH_TREM2_HSC_Integrated_Ligand_Audit_v2.1.1.R
- Human_MASH_Focused_Ligand_Closure_v2.2.0.R

## Next analysis
Before switching to Mouse Sham vs Tx validation, re-audit HSC subfraction classification.

### HSC subfraction audit priorities
1. Confirm primary qHSC vs activated/fibroblast-mesenchymal HSC labels.
2. Audit marker concordance without reclustering.
3. Separate lineage/subfraction labels from functional programs.
4. Evaluate whether finer HSC states are supportable:
   - quiescent HSC
   - activated/myofibroblastic HSC
   - inflammatory/cytokine-responsive HSC
   - proliferating/cycling HSC
   - pericyte/VSMC-like contamination/sensitivity state
5. Retain current canonical HSC receiver unless a clearly justified relabeling is supported.

## Planned Git checkpoint
Suggested tag:
`human-mash-mphi-hsc-focused-ligand-closure-v2.2.0`
