# Human MASH HSC Classification Checkpoint v1.4
Date: 2026-11-18

## Canonical HSC receiver
- Total HSC cells: 7,640
- Identity annotation is retained:
  - qHSC: 4,335 cells
  - aHSC_mesenchymal: 3,305 cells
- No reclustering.
- No source RDS overwrite.

## Identity validation
Human_MASH_HSC_Subfraction_Validation_v1.3.2 confirmed that the current two-state identity annotation is directionally valid.

### qHSC support
Across all 8 biological samples:
- Quiescent score was higher in qHSC than aHSC_mesenchymal in 8/8 samples.
- Median qHSC-minus-aHSC delta = 0.7213.
- Wilcoxon P = 0.01427.
- FDR = 0.02378.

### aHSC support
Across all 8 biological samples:
- Composite identity delta shifted toward activated identity in aHSC_mesenchymal in 8/8 samples.
- Median aHSC-minus-qHSC identity delta = 0.6929.
- P = 0.01427.
- FDR = 0.02378.

### Important limitation
The myofibroblastic score alone was not consistently higher in aHSC_mesenchymal:
- Median delta = -0.0472.
- Positive in only 3/8 samples.
- P = 0.363.

Therefore:
- aHSC_mesenchymal should NOT be renamed "myofibroblastic HSC".
- It is better interpreted as an activated/mesenchymal identity with heterogeneous functional programs.

## Functional overlays

### IEG / initiatory activation
This is the strongest reproducible functional distinction between the two identity classes.

In v1.3.2:
- IEG overlay was higher in aHSC_mesenchymal in 8/8 samples.
- Median delta = 0.0892.
- P = 0.01427.
- FDR = 0.02378.

In v1.4.0 using sample-matched qHSC p90 thresholds:
- IEG-high fraction was higher in aHSC in 8/8 samples.
- Median delta = 0.0849.
- P = 0.01427.
- FDR = 0.03567.
- IEG score was higher in aHSC in 8/8 samples.
- Median delta = 0.0882.
- P = 0.01427.
- FDR = 0.03567.

### ECM / fibrogenic program
ECM-high is a functional state, not an identity subtype.

In v1.4.0:
- ECM-high fraction did not differ consistently between aHSC and qHSC.
- Median delta = -0.00292.
- Positive in 4/8 samples.
- P = 1.0.
- ECM score was not consistently higher in aHSC.
- Median delta = -0.0521.
- P = 0.234.

Therefore:
- ECM/fibrogenic activation must remain an independent functional axis.
- Do not equate aHSC_mesenchymal with ECM-high or myofibroblastic HSC.

### Cycling
Cycling HSC is not supported as a meaningful subfraction in this dataset.
Direct MKI67/TOP2A evidence:
- aHSC_mesenchymal: 7/3305 = 0.21%
- qHSC: 21/4335 = 0.48%
- Two-or-more cycling genes: essentially absent.

Therefore:
- Do not use Cycling as a formal HSC subtype.

### Pericyte/VSMC and portal-fibroblast programs
- Retain as sensitivity flags only.
- Do not exclude or relabel cells automatically.
- RGS5 alone is insufficient for pericyte/VSMC reassignment.
- Portal-fibroblast assignment is limited by missing PI16/SFRP2.

## Final HSC framework

### Identity layer
1. qHSC
2. aHSC_mesenchymal

### Functional overlay layer
1. IEG / initiatory activation
2. ECM / fibrogenic output

Optional sensitivity flags:
- Pericyte/VSMC-like
- Portal-fibroblast-like

Not retained as formal subtypes:
- Cycling HSC
- Myofibroblastic HSC as a replacement for all aHSC_mesenchymal cells

## Critical interpretation of v1.4 NAS correlations
The v1.4 primary functional-state calls use a sample-specific qHSC p90 threshold.

This is appropriate for:
- comparing aHSC vs qHSC within the same sample
- asking whether aHSC is more activated than its matched qHSC baseline

It is NOT appropriate as the primary method for:
- comparing absolute activation fractions across samples
- correlating high-state fractions with NAS

Reason:
- if qHSC itself becomes activated with disease severity, the sample-specific qHSC p90 threshold also rises.
- This can suppress or reverse cross-sample disease-severity associations.

Therefore:
- v1.4 NAS correlations are descriptive only and should not be interpreted as evidence that IEG/ECM activation decreases with MASH severity.
- Cross-sample disease associations should instead use fixed-reference or pseudobulk approaches.

## Scripts retained
- Human_MASH_HSC_Subfraction_Validation_v1.3.2.R
- Human_MASH_HSC_Activation_State_Matrix_v1.4.0.R

Superseded / not for checkpoint:
- Human_MASH_HSC_Subfraction_Audit_v1.3.0.R
- Human_MASH_HSC_Subfraction_Audit_v1.3.1.R

## Recommended next step
Return to the main mechanistic program:
Mouse Sham vs Tx validation of:
- Spp1 axis
- Cdcp1 axis
- Tgfb1 axis
with HSC response evaluated using the same two-layer framework:
- identity: qHSC / aHSC
- function: IEG / initiatory and ECM / fibrogenic programs
