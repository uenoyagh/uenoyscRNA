# Human–Mouse MASH MΦ→HSC cross-species hierarchy checkpoint v1.1

**Date:** 2026-11-19  
**Project:** uenoyscRNA

## Canonical distinction

### Human mechanistic priority
1. SPP1 → HSC integrins
2. CDCP1-associated → HSC EGFR-supported axis
3. TGFB1 → TGFβ receptor axis
4. PDGFA → PDGFRA

### Cross-species therapeutic-response bridge priority
1. PDGF-family → HSC PDGFR
2. CDCP1/Cdcp1-associated
3. SPP1/Spp1-integrin
4. TGFB1/Tgfb1-TGFβR

These two rankings answer different questions and must not be conflated.

## Cross-species interpretation

### 1. PDGF-family → HSC PDGFR
**Current strongest cross-species functional therapeutic-response bridge.**

- Human: PDGFA–PDGFRA is associated with the ECM-high HSC phenotype.
- Mouse: Pdgfa–Pdgfra is approximately unchanged, whereas Pdgfb–Pdgfrb is attenuated with Tx.
- Mouse Pdgfb–Pdgfrb Tx/Sham ≈ 0.746.
- Mouse host-weighted ECM-activated HSC core-collagen Tx/Sham ≈ 0.925.
- Interpretation: conserved PDGF-family → HSC function with context/species-specific paralog usage divergence.

Boundary:
- Do not claim one-to-one PDGFA–PDGFRA conservation.
- Do not claim causal PDGF dependence from these observational/descriptive data.

### 2. CDCP1/Cdcp1-associated
**Strongest directionally concordant gene-level bridge, but not a conserved dominant macrophage-source mechanism.**

Corrected Mouse full-library source decomposition:
- MΦ Cdcp1 raw-UMI share = 0.01785 (~1.78%)
- non-MΦ Cdcp1 raw-UMI share = 0.98215 (~98.22%)
- Cholangiocyte + Hepatocyte share = 0.88068 (~88.07%)
- TREM2-associated MΦ share within the MΦ Cdcp1 pool = 0.39048 (~39.05%)
- MΦ Cdcp1 CPM Tx/Sham = 0.61838
- non-MΦ Cdcp1 CPM Tx/Sham = 0.70559

Interpretation:
- Human CDCP1 is the strongest novel TREM2-linked HSC-associated candidate.
- A sparse Mouse MΦ/TREM2-associated Cdcp1 component is present and decreases with Tx.
- Mouse whole-liver Cdcp1 is overwhelmingly non-MΦ-derived.
- Therefore source conservation of a dominant macrophage-derived CDCP1/Cdcp1 axis is not supported.

Do not claim:
- Mouse Cdcp1 is MΦ-derived or MΦ-specific.
- Direct CDCP1→EGFR causality.

### 3. SPP1/Spp1-integrin
- Human: primary mechanistic/fibrotic-niche axis.
- Mouse: Spp1 is not attenuated by Tx.
- Interpretation: conserved profibrotic biology, but not supported as the main Tx-response attenuation mechanism.

### 4. TGFB1/Tgfb1-TGFβR
- Human: robust auxiliary profibrotic pathway.
- Mouse: Tgfb1 is not attenuated by Tx.
- Interpretation: conserved canonical pathway, but not supported as the main Tx-response attenuation mechanism.

## Inference boundaries

- Human: exploratory sample-level inference, n=8.
- Mouse: descriptive Sham vs Tx inference, n=2/group.
- Cross-species directional concordance does not establish causality.
- Receptor expression does not establish biochemical receptor activation.
- Functional-axis conservation does not require identical ligand paralogs.

## Integrity

- Input v1.0 rows: 4
- Output axes: 4
- Human priority preserved: TRUE
- Cdcp1 full-library denominator verified: TRUE
- Source RDS overwritten: FALSE
- Reclustering performed: FALSE
- STATUS: SUCCESS
- ANALYSIS COMPLETED

## Canonical outputs

- `Mouse_MASH_Cdcp1_WholeCell_SourceDecomposition_v6.4.11.1`
- `Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.1`

## Recommended manuscript-safe model

Human and Mouse MASH data support a conserved multi-axis myeloid-to-HSC fibrogenic network with species/context-dependent ligand usage. The PDGF-family → HSC PDGFR axis currently provides the strongest cross-species functional therapeutic-response bridge, with Human PDGFA–PDGFRA versus Mouse Pdgfb–Pdgfrb divergence. CDCP1/Cdcp1 shows strong directional concordance across species, but Mouse whole-liver Cdcp1 is predominantly non-macrophage-derived, arguing against a conserved dominant macrophage-source mechanism. SPP1 and TGFB1 remain biologically conserved profibrotic axes but are not attenuated in the Mouse Tx response.
