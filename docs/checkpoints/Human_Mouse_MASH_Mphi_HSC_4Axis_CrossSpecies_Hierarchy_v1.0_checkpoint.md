# Human–Mouse MASH MΦ→HSC 4-axis cross-species hierarchy checkpoint v1.0

**Date:** 2026-11-19
**Project:** uenoyscRNA
**Scope:** Human MASH LT1–LT8 + Mouse CDAHFD Sham vs Tx

## Canonical cross-species hierarchy

| Human priority | Functional axis | Human evidence | Human rho | Human NAS-adjusted rho | Mouse Tx axis ratio | Mouse ECM-HSC core-collagen ratio | Cross-species interpretation |
|---|---|---|---:|---:|---:|---:|---|
| 1 | SPP1/Spp1–integrin | HIGH — primary mechanistic axis | 0.786 | 0.649 | 1.034 | 0.925 | Conserved profibrotic candidate biology; Mouse Tx response is not explained by Spp1 attenuation |
| 2 | CDCP1/Cdcp1-associated | HIGH association / MODERATE receptor mechanism | 0.952 | 0.948 | 0.732 | 0.925 | Strongest directionally concordant Human–Mouse candidate bridge; Mouse source attribution remains unresolved |
| 3 | TGFB1/Tgfb1–TGFβR | HIGH — auxiliary pathway | 0.881 | 0.828 | 1.021 | 0.925 | Canonical pathway conservation; Tgfb1 attenuation does not explain Mouse Tx response |
| 4 | PDGF-family → HSC PDGFR | MODERATE — robust exploratory Human axis | 0.778 | 0.503 | 0.746 (Mouse Pdgfb–Pdgfrb) | 0.925 | Functional-axis conservation with context/species-specific paralog usage divergence |

## Mouse treatment-response details

### SPP1/Spp1
- Total-myeloid Spp1 Tx/Sham = 1.0336
- TREM2-associated Spp1 Tx/Sham = 1.0565
- Host-weighted ECM-activated HSC core-collagen Tx/Sham = 0.9252
- Interpretation: Spp1-integrin biology may be conserved, but Spp1 attenuation is not supported as the principal Mouse Tx mechanism.

### CDCP1/Cdcp1
- Total-myeloid Cdcp1 Tx/Sham = 0.7323
- TREM2-associated Cdcp1 Tx/Sham = 0.7132
- Host-weighted ECM-activated HSC core-collagen Tx/Sham = 0.9252
- Independent raw-count pseudobulk:
  - TREM2-state Cdcp1 delta log2CPM = -0.8747
  - rank_down = 218
  - percentile_down = 0.01847
- Interpretation: strongest directionally concordant Human–Mouse bridge.
- Boundary: Mouse Cdcp1 is sparse and substantially more abundant outside macrophages; an endogenous macrophage component is plausible but ambient/global contribution is not excluded.

### TGFB1/Tgfb1
- Total-myeloid Tgfb1 Tx/Sham = 1.0210
- TREM2-associated Tgfb1 Tx/Sham = 1.0272
- Host-weighted ECM-activated HSC core-collagen Tx/Sham = 0.9252
- Interpretation: canonical TGFβ profibrotic pathway is conserved, but Tgfb1 attenuation is not supported as the principal Mouse Tx mechanism.

### PDGF family
Human:
- PDGFA–PDGFRA is associated with ECM-high HSC phenotype.

Mouse:
- total Pdgfa = 1.1309
- TREM2 Pdgfa = 1.1386
- total Pdgfa×Pdgfra = 0.9898
- total Pdgfb×Pdgfrb = 0.7460
- TREM2 Pdgfb×Pdgfrb = 0.8532

Interpretation:
- Cross-species conservation is best defined at the **PDGF-family → HSC PDGFR functional-axis** level.
- Human preferentially supports PDGFA–PDGFRA.
- Mouse Tx response preferentially attenuates Pdgfb–Pdgfrb.
- This is treated as context/species-specific paralog usage divergence, not failure of axis conservation.

## Canonical receiver endpoint

Mouse treatment-response comparisons use:

`host-weighted ECM-activated HSC core-collagen burden`

Sham → Tx ratio = 0.9252.

Simple per-cell HSC intensity is not substituted for this canonical receiver endpoint.

## Inference boundaries

- Human: exploratory sample-level inference, n=8.
- Mouse: descriptive Sham vs Tx inference, n=2/group.
- Mouse n=4 correlations are descriptive only.
- Cross-species concordance does not establish causal ligand dependence or receptor activation.
- `Mouse_MASH_WholeHost_Provenance_StressAudit_v6.4.10.4` is not used quantitatively in this checkpoint because the full-library CPM denominator correction remains pending.

## Integrity audit

- Detailed cross-species TSV: 28 columns on every row.
- Concise cross-species TSV: 10 columns on every row.
- Rows = 4.
- Human priority = 1,2,3,4.
- All Mouse receiver ratios finite.
- All canonical Mouse receiver ratios < 1.
- Result: `INTEGRITY OK`.
- Source RDS overwritten: FALSE.
- Reclustering performed: FALSE.

## Canonical input analyses

Human:
- `Human_MASH_Mphi_HSC_4Axis_Hierarchy_v2.3.1`
- `Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2`

Mouse:
- `Mouse_MASH_Mphi_HSC_HumanAxis_Validation_v6.4.1`
- `Mouse_MASH_HSC_3Axis_TotalHostWeighted_v6.4.4`
- `Mouse_MASH_Mphi_HSC_IntegratedClosure_v6.4.5`
- `Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1`

Cross-species integration:
- `Human_Mouse_MASH_Mphi_HSC_4Axis_CrossSpecies_Hierarchy_v1.0`

## Working model

Human and Mouse support a conserved multi-axis myeloid-to-HSC fibrogenic network, but individual ligand usage and treatment responsiveness differ by species/context.

1. **SPP1/Spp1–integrin:** conserved primary fibrotic-niche biology; not attenuated by Mouse Tx.
2. **CDCP1/Cdcp1-associated:** strongest directionally concordant Human–Mouse candidate therapeutic bridge, with unresolved Mouse source attribution.
3. **TGFB1/Tgfb1–TGFβR:** conserved canonical auxiliary pathway; not attenuated by Mouse Tx.
4. **PDGF-family → HSC PDGFR:** conserved functional axis with Human PDGFA–PDGFRA versus Mouse Pdgfb–Pdgfrb divergence.

## Next unresolved issue

Before strengthening the Mouse CDCP1/Cdcp1 source-specific mechanistic claim, correct the whole-host provenance analysis to use the full raw RNA library size as the CPM denominator, then repeat the Cdcp1 whole-cell source decomposition.
