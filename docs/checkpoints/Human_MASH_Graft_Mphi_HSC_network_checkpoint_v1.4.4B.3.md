# Human MASH graft → host MΦ → HSC network checkpoint v1.4.4B.3

**Date:** 2026-09-09  
**Scope:** Human graft-state remodeling → mouse host macrophage receiver competence → macrophage remodeling → MΦ-to-HSC PDGF attenuation → reduced activated/cycling HSC burden.

## Status

**Network-level mechanistic integration: CLOSED WITH EXPLICIT CAUSAL BOUNDARIES**

This checkpoint integrates:
- Human graft-side ligand/state analyses through v1.4.4A.4.2
- Graft-ligand → mouse MΦ 5-state receptor bridge v1.4.4B.1
- Receiver-weighted MΦ pool → HSC network integration v1.4.4B.2
- Four-sample IL10 → MΦ → HSC directional concordance v1.4.4B.3
- Existing mouse MΦ–HSC mechanistic closure v6.3.6
- Existing direct IL10R-status → Pdgfb negative boundary v1.4.3C.2.2

## Primary biological model

Engrafted IL10 organoid grafts show graft-state remodeling and retain or enhance multiple candidate paracrine ligands. In the mouse host, transplantation is associated with expansion and favorable-state redistribution of an IL10R-transcript-competent macrophage compartment. This occurs together with reduced Repair/Resolution-MΦ Pdgfb output, selective attenuation of macrophage-to-ECM-HSC PDGF signaling, and reduced activated/cycling fibrogenic HSC burden.

The supported claim is **network-level convergence**, not a directly demonstrated intracellular IL10R → Pdgfb causal chain.

## v1.4.4B.1 — Graft ligand → mouse MΦ 5-state receiver bridge

Primary mouse MΦ states:
- Inflammatory-Mphi
- Anti-inflammatory-Mphi
- ECM-associated inflammatory-Mphi
- Repair/Resolution-Mphi
- Lipid-associated/TREM2-Mphi

`Other` is retained only as context and excluded from the primary 5-state ranking.

Frozen receptor-complex logic inherited from v1.4.3C.1.1:
- IL10 → Il10ra AND Il10rb
- FGF1 → ANY(Fgfr1,Fgfr2,Fgfr3,Fgfr4)
- TNFSF12 → Tnfrsf12a
- VEGFA → Flt1 OR Kdr
- CXCL12 → Cxcr4 OR Ackr3
- TGFB3 → Tgfbr1 AND Tgfbr2
- BMP5 / BMP2 → Bmpr2 AND ANY(Bmpr1a,Acvr1)
- HBEGF / BTC → Egfr OR Erbb4

Unmapped in frozen C.1.1 and intentionally not added without a separate receptor-definition audit:
- FST
- IGFBP2
- IGFBP4
- INHBA
- VEGFB

## v1.4.4B.2 — Receiver-weighted primary MΦ pool

B.1 denominator validation against v1.4.3C.2.1: **PASS**

Primary-5 MΦ totals:
- Sham1: 2,949
- Sham20: 2,498
- Tx17: 3,067
- Tx5: 3,718

Pooled receptor-complex-positive fractions, Sham mean → Tx mean:
- CXCL12 receiver: 0.178141 → 0.303352; Δ +0.125211; complete Tx-higher separation
- IL10R: 0.311171 → 0.364579; Δ +0.053408; complete Tx-higher separation
- TNFSF12 receiver: 0.075430 → 0.105892; Δ +0.030462; complete Tx-higher separation
- TGFB3 receiver: 0.114356 → 0.144621; Δ +0.030265; complete Tx-higher separation
- VEGFA receiver: 0.038146 → 0.054810; Δ +0.016664; complete Tx-higher separation
- BMP5/BMP2 receiver: 0.005670 → 0.009129; overlap
- HBEGF/BTC receiver: 0.037030 → 0.033998; overlap
- FGF1 receiver: 0.240339 → 0.227631; overlap

### IL10 receiver-pool state decomposition

Mean state fraction of primary-5 MΦ, Sham → Tx:
- Anti-inflammatory-Mphi: 0.288814 → 0.334347
- Lipid-associated/TREM2-Mphi: 0.323814 → 0.353128
- Repair/Resolution-Mphi: 0.027794 → 0.052948
- ECM-associated inflammatory-Mphi: 0.089785 → 0.078554
- Inflammatory-Mphi: 0.269793 → 0.181022

Mean within-state IL10R dual-positive fraction, Sham → Tx:
- Anti-inflammatory-Mphi: 0.363040 → 0.403984
- Lipid-associated/TREM2-Mphi: 0.318594 → 0.356187
- Repair/Resolution-Mphi: 0.279798 → 0.336281
- ECM-associated inflammatory-Mphi: 0.212825 → 0.273333
- Inflammatory-Mphi: 0.284874 → 0.358293

Contribution to Tx − Sham change in the receiver-weighted IL10R-competent MΦ pool:
- Anti-inflammatory-Mphi: +0.030576
- Lipid-associated/TREM2-Mphi: +0.020787
- Repair/Resolution-Mphi: +0.009794
- ECM-associated inflammatory-Mphi: +0.002136
- Inflammatory-Mphi: −0.009884

Interpretation:
Tx is associated with redistribution of the IL10R-transcript-competent macrophage pool toward anti-inflammatory / repair-associated states, rather than a uniform increase of every macrophage state.

## v1.4.4B.3 — Four-sample network directional concordance

Directional evidence ladder:

| Layer | Metric | Sham mean | Tx mean | Tx−Sham | Expected | Replicate relation |
|---|---|---:|---:|---:|---|---|
| Host receiver | IL10R-dual-positive fraction across primary five MΦ states | 0.311171 | 0.364579 | +0.053408 | Tx up | complete Tx-higher |
| MΦ remodeling | Anti-inflammatory-MΦ fraction | 0.288814 | 0.334347 | +0.045534 | Tx up | complete Tx-higher |
| MΦ remodeling | Repair/Resolution-MΦ fraction | 0.027794 | 0.052948 | +0.025154 | Tx up | complete Tx-higher |
| MΦ remodeling | Inflammatory-MΦ fraction | 0.269793 | 0.181022 | −0.088771 | Tx down | complete Tx-lower |
| MΦ effector | Repair/Resolution-MΦ Pdgfb per-cell CP10k | 5.970686 | 3.568014 | −2.402672 | Tx down | complete Tx-lower |
| MΦ→HSC | CellChat PDGF Repair/Resolution-MΦ → ECM-HSC | 0.184601 | 0.108756 | −0.075844 | Tx down | complete Tx-lower |
| HSC receiver | ECM-HSC PDGF mitogenic module | +0.478867 | −0.478867 | −0.957734 | Tx down | complete Tx-lower |
| HSC state | ECM-activated HSC fraction within HSC3 | 0.527266 | 0.452309 | −0.074956 | Tx down | complete Tx-lower |
| HSC state | Cycling-like fraction within ECM-HSC | 0.055968 | 0.047825 | −0.008143 | Tx down | overlap |
| HSC state | Cycling ECM-HSC burden within HSC3 | 0.029815 | 0.021581 | −0.008234 | Tx down | overlap |

Summary:
- Model-direction concordance: **10 / 10**
- Complete replicate-range separation: **8 / 10**
- The two cycling-related metrics retain Tx-down mean direction but show replicate-range overlap.

This is a descriptive directional concordance analysis of selected biologically linked metrics, not an independent statistical validation.

## Human graft-side findings carried into this checkpoint

### IL10 anchor
Standard filtered matrix:
- SfLB1: 2 UMI / 2 cells
- SfLB2: 0
- Tx17: 18 UMI / 13 cells
- Tx5: 4 UMI / 4 cells

Strict shared-probe sensitivity:
- SfLB1: 0
- SfLB2: 0
- Tx17: 5
- Tx5: 0

Paper-safe interpretation:
IL10 is the engineered anchor by study design and is detected in both Tx graft samples in the standard filtered gene-expression matrix with molecule/probe-level reproduction. Strict shared-probe sensitivity is heterogeneous and retains signal in Tx17 but not Tx5.

### Non-IL10 graft ligands

Tx-higher with complete replicate-range separation:
- FGF1
- TNFSF12
- VEGFA
- TGFB3
- BMP5
- BMP2
- BTC
- FST
- IGFBP2

Tx-mean higher but replicate ranges overlap:
- CXCL12
- HBEGF
- IGFBP4

No clear Tx enhancement:
- INHBA
- VEGFB

These are candidate graft-derived paracrine/tissue-remodeling factors. They must not all be labeled anti-inflammatory or anti-fibrotic.

### Hepatocyte-like graft-state remodeling

Strict vs standard direction concordance across 16 measurable hepatic genes:
- 15 / 16 = 93.75%
- Exact replicate-relation match: 14 / 16 = 87.5%

Robust Tx-associated pattern:
- HNF4A ↑
- CEBPA ↑
- APOA1 ↓
- FGB ↓
- TTR ↓
- ASS1 ↓
- OTC ↓
- ABCB11 ↓
- FABP1 ↓
- MTTP ↓
- EPCAM ↓
- SOX9 ↓
- ASL ↑
- NAGS overlap
- KRT7 overlap
- FGG standard/strict discordant

Interpretation:
Tx-associated hepatocyte-like cells show stronger hepatocyte-lineage transcription-factor identity with selective remodeling/reduction of several probe-measurable mature hepatic functional transcripts. This is not a simple immature state and must not be called global liver-function decline.

## Explicit causal / interpretation boundaries

Do **not** claim:
1. IL10 transcript proves secreted bioactive IL10 protein.
2. Il10ra+Il10rb transcript co-detection proves receptor occupancy or pathway activation.
3. Stat3 mRNA is evidence of pSTAT3.
4. Direct cell-autonomous IL10R activation in Repair/Resolution-MΦ suppresses Pdgfb.
5. PDGF is the sole mediator of HSC suppression or fibrosis improvement.
6. SfLB → Tx is an IL10-only perturbation.
7. Genes absent from the Flex feature space are biological negatives.
8. The 10/10 directional concordance is an independent inferential test.

Supported wording:
> Together, the data support a network-level model in which engrafted IL10 organoid grafts are associated with expansion and favorable-state redistribution of a receptor-transcript-competent host macrophage compartment, selective attenuation of macrophage-to-ECM-HSC PDGF signaling, and reduced activated/cycling fibrogenic HSC burden.

## Key scripts to checkpoint

- `scripts/Human_MASH/Human_MASH_v1.4.4B.1_GraftLigand_MouseMphi5_ReceiverBridge.R`
- `scripts/Human_MASH/Human_MASH_v1.4.4B.2_GraftReceiverWeighted_MphiHSC_NetworkIntegration.R`
- `scripts/Human_MASH/Human_MASH_v1.4.4B.3_IL10_Mphi_HSC_SampleLevel_NetworkConcordance.R`

Human graft-side scripts used by the integrated interpretation should already be checkpointed or tracked separately; do not stage unrelated failed/superseded scripts automatically.

## Key SSD outputs

### B.1
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.1/analysis/v1.4.4B.1_graft_ligand_mouse_Mphi5_receiver_bridge`

### B.2
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.2/analysis/v1.4.4B.2_graft_receiver_weighted_Mphi_HSC_network_integration`

### B.3
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.3/analysis/v1.4.4B.3_IL10_Mphi_HSC_sample_level_network_concordance`

Recommended summary directory:
`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.3/analysis/v1.4.4B.3_IL10_Mphi_HSC_sample_level_network_concordance/summary`

Recommended summary files:
- `Human_MASH_Graft_Mphi_HSC_network_summary_v1.4.4B.3.png`
- `Human_MASH_Graft_Mphi_HSC_network_evidence_v1.4.4B.xlsx`

## Recommended Git checkpoint

Commit message:
`Checkpoint human MASH graft-Mphi-HSC network v1.4.4B.3`

Annotated tag:
`human-mash-graft-mphi-hsc-network-v1.4.4B.3`

Only explicitly verified successful scripts and this checkpoint markdown should be staged. Do not use `git add .`.

## Summary artifact integrity

Summary PNG:
- filename: `Human_MASH_Graft_Mphi_HSC_network_summary_v1.4.4B.3.png`
- dimensions: `1672 x 941 px`
- SHA256: `90f4d4a1bcc12dd799802eb3fd7d7f23ed75f4246a7f993cf6254d5c460ae3ef`

Evidence Excel:
- filename: `Human_MASH_Graft_Mphi_HSC_network_evidence_v1.4.4B.xlsx`
- SHA256: `229206d24cc9185c1c1c6d41dae1b26fa1a291d5a71fcfbce7a4c77428ffad74`

The PNG and Excel are stored on the external SSD as analysis outputs and are not committed to Git.

