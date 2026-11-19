# Human MASH MΦ→HSC 4-axis mechanistic hierarchy checkpoint v2.3.1

**Date:** 2026-11-19  
**Project:** uenoyscRNA  
**Scope:** Human MASH LT1–LT8, sample-level exploratory inference (n=8)

## Canonical hierarchy

| Priority | Axis | Sample-level rho | NAS-adjusted rho | NicheNet rank | Evidence position |
|---|---|---:|---:|---:|---|
| 1 | SPP1 → HSC integrins | 0.786 | 0.649 | 4 | HIGH — primary mechanistic axis |
| 2 | CDCP1-associated → HSC EGFR-supported axis | 0.952 | 0.948 | 20 | HIGH for association / MODERATE for receptor mechanism |
| 3 | TGFB1 → TGFβ receptor axis | 0.881 | 0.828 | 206 | HIGH as auxiliary pathway |
| 4 | PDGFA → PDGFRA | 0.778 | 0.503 | NA | MODERATE — robust exploratory ECM-high axis |

## Canonical receiver display

- SPP1: `ITGA4;ITGA5;ITGA9;ITGAV;ITGB1;ITGB5`
- CDCP1: `EGFR (expression support; direct CDCP1-EGFR LR not established)`
- TGFB1: `TGFBR1;TGFBR2;TGFBR3`
- PDGFA: `PDGFRA`

Raw receptor-prior support is retained separately for provenance/audit and must not be substituted for the canonical mechanistic display.

## Interpretation

### 1. SPP1 → HSC integrins
Primary Human MASH MΦ→HSC profibrotic axis. Strong receptor support, NicheNet rank 4, and broad target support. Use the phrase **“myeloid SPP1 → HSC integrin axis”**. Do not call it a TREM2-specific direct axis and do not imply that one integrin is uniquely causal.

### 2. CDCP1-associated signaling
Strongest novel TREM2-linked candidate. It has the strongest sample-level and NAS-adjusted association among the focused ligands. HSC EGFR expression supports receiver competence, but **direct CDCP1→EGFR ligand-receptor causality is not established**.

### 3. TGFB1 → TGFβ receptors
Robust TREM2-state-associated auxiliary profibrotic pathway with broad target support. TREM2-state association does not imply that TREM2-associated macrophages are the dominant absolute TGFB1 source. Do not call TGFB1 the dominant direct ligand.

### 4. PDGFA → PDGFRA
Cross-state myeloid PDGFA-associated signaling is linked to a PDGFRA-positive ECM-high HSC phenotype. `total_PDGFA × HSC_PDGFRA → ECM_high_fraction`: rho=0.778, NAS-adjusted rho=0.503, leave-one-out sign consistency=1. TREM2-associated PDGFA is sparse, so the axis is not considered TREM2-specific. Equivalent focused NicheNet target closure has not yet been performed. Do not claim that PDGFA directly drives HSC proliferation from the current Human scRNA-seq data.

## Integrated working model

Human MASH HSC fibrogenic activation is best represented as a **multi-axis myeloid signaling network** rather than a single-ligand mechanism:

**SPP1-integrin signaling** forms the primary fibrotic-niche axis; **CDCP1-associated signaling** is the strongest novel TREM2-linked candidate; **TGFB1/TGFβ signaling** provides a canonical auxiliary profibrotic program; and **PDGFA/PDGFRA signaling** adds a cross-state myeloid axis associated with the ECM-high HSC phenotype.

All Human sample-level associations remain exploratory because `n=8`. Correlation, receptor transcript detection, and ligand-target prediction do not establish biochemical receptor activation or causal ligand dependence.

## Canonical analysis outputs

- `Human_MASH_Focused_Ligand_Closure_v2.2`
- `Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2`
- `Human_MASH_PDGF_Robustness_Adjudication_v1.1.0`
- `Human_MASH_PDGFA_SenderState_Attribution_v1.2.0`
- `Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.2`
- `Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1`

## v2.3.1 integrity checks

- Axes: `SPP1;CDCP1;TGFB1;PDGFA`
- All four axes: `n=8`
- NicheNet ranks validated: SPP1=4, CDCP1=20, TGFB1=206
- Concise TSV: 12 columns on all rows
- Detailed TSV: 34 columns on all rows
- Canonical receptor labels verified
- Source RDS overwritten: FALSE
- Reclustering performed: FALSE
- PDGFA equivalent focused NicheNet closure: NOT YET PERFORMED

## Superseded / do not use as canonical

- `Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.R`
- `Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.1.R`
- `Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.0.R`
- `Human_MASH_PDGFA_ValidatedHSC_Phenotype_v1.3.1.R`
- `Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.0.R`

## Next step

Proceed to a **Human–Mouse cross-species MΦ→HSC mechanistic hierarchy** while preserving the different evidence boundaries between species. PDGFA-focused NicheNet closure can be added later if evidence-level harmonization across all four Human axes is required.
