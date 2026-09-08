# Human MASH IL10 translational mechanistic closure checkpoint — v1.4.3E.0

**Date:** 2026-09-09
**Scope:** `v1.4.3B.2C.3 -> v1.4.3C.1.1 -> v1.4.3C.2.1 / v1.4.3C.2.2 -> mouse MΦ–HSC v6.3.6 -> v1.4.3D.0 / v1.4.3D.1`

## 1. Overall closure call

**MECHANISTIC CLOSURE AT THE NETWORK LEVEL: ACHIEVED**

The integrated evidence supports the following cross-species working model:

> **Human engineered graft IL10 / graft secretome -> receptor-competent mouse host myeloid compartment -> Tx-associated macrophage-state remodeling -> selective attenuation of macrophage-to-HSC PDGF mitogenic signaling -> reduced activated/cycling HSC burden -> fibrosis improvement, with the corresponding IL10 receptor-transcript-competent myeloid target compartment retained in advanced human MASH.**

The following narrower intracellular claim is **not closed**:

> **IL10R activation in an individual Repair/Resolution-MΦ directly suppresses Pdgfb.**

`v1.4.3C.2.2` specifically failed to demonstrate a Tx-specific, replicate-consistent direct IL10R-status -> `Pdgfb` relationship across all sensitivity metrics. This negative boundary must remain explicit.

## 2. Mechanistic closure table

| Link | Mechanistic question | Main result | Strength | Closure |
|---|---|---|---|---|
| L1 | Human graft -> IL10 transcript evidence | IL10 detected in Tx17 and Tx5 and reproduced at molecule/probe level; stringent shared-probe sensitivity retains Tx17 signal | High for transcript detection | CLOSED WITH BOUNDARY |
| L2 | Graft IL10 -> mouse host receiver | Host macrophage/monocyte compartments show substantial `Il10ra`/`Il10rb` competence | High | CLOSED WITH BOUNDARY |
| L3 | Tx -> macrophage remodeling | Inflammatory/ECM-inflammatory MΦ decrease; anti-inflammatory and Repair/Resolution MΦ increase | High | CLOSED WITH NEGATIVE SIGNALING BOUNDARY |
| L4 | IL10R-positive MΦ -> direct `Pdgfb` suppression | No Tx-specific, replicate-consistent direct link across per-cell/pseudobulk metrics | Negative/limiting evidence | **NOT CLOSED DIRECTLY** |
| L5 | Tx -> MΦ-to-HSC PDGF axis | Repair/Resolution-MΦ `Pdgfb` and PDGF communication to ECM-HSC decrease; ECM-HSC mitogenic module decreases | High | CLOSED, PATHWAY-SELECTIVE |
| L6 | PDGF attenuation -> HSC consequence | ECM-activated/cycling HSC burden decreases | High | CLOSED WITH BOUNDARY |
| L7 | Mouse mechanism -> human MASH target conservation | Advanced human MASH retains same-cell `IL10RA+IL10RB` myeloid cells in all PRIMARY states | High for target availability | CLOSED TRANSLATIONALLY |
| L8 | Integrated model | Network-level cross-species mechanism is coherent; direct intracellular IL10R->`Pdgfb` chain remains unproven | High at network level | **MECHANISTIC CLOSURE ACHIEVED AT NETWORK LEVEL** |

## 3. Detailed evidence ledger

### L1. Human engineered graft -> IL10 transcript evidence

- **Source:** v1.4.3B.2C.3
- **Biological unit:** Human iPSC-derived graft cells / molecules
- **Key validated evidence:** Probe/molecule-level reconstruction reproduced canonical IL10 counts in Tx17 (18 UMI / 13 cells) and Tx5 (4 UMI / 4 cells). SfLB1 had 2 singleton MEX events and SfLB2 had 0; stringent all-four shared-TRUE probe sensitivity was SfLB1=0, SfLB2=0, Tx17=5, Tx5=0.
- **Direction / result:** IL10 transcript signal is reproducibly present in both Tx graft samples; stringent shared-probe sensitivity retains Tx17 signal.
- **Replicate consistency:** Tx17 and Tx5 both positive in canonical gene-level/molecule reconstruction; strict shared-probe sensitivity is Tx17-positive / Tx5-zero.
- **Evidence strength:** High for transcript detection; Moderate for transgene-specific attribution
- **Causal boundary:** Gene-level IL10 cannot by itself distinguish endogenous from transgene-derived transcript; transcript detection does not demonstrate secreted IL10 protein.
- **Paper-safe wording:** IL10 was detected in both Tx graft samples in the standard filtered gene-expression matrix and reproduced at molecule/probe level. A stringent shared-probe sensitivity analysis retained signal in Tx17 but not Tx5.
- **Closure status:** `CLOSED_WITH_BOUNDARY`

### L2. Human graft ligands -> mouse host receiver landscape

- **Source:** v1.4.3C.1.1
- **Biological unit:** Mouse host broad liver lineages, sample-level descriptive summaries
- **Key validated evidence:** IL10R dual-positive competence was concentrated in myeloid lineages relevant to the treatment mechanism; Tx mean dual-positive fractions were ~34.4% in Kupffer/Macrophage and ~25.2% in Monocyte. Other ligand-receptor bridges showed lineage-specific patterns (e.g., VEGFA->LSEC; FGF1/TNFSF12/HBEGF/BTC->hepatobiliary/HSC compartments).
- **Direction / result:** Mouse host myeloid cells are a major receptor-competent compartment for the engineered IL10 anchor.
- **Replicate consistency:** Receiver competence retained in both Tx17 and Tx5; descriptive n=2 Sham vs n=2 Tx.
- **Evidence strength:** High for receptor-transcript competence
- **Causal boundary:** Receptor transcript co-detection/competence is not ligand binding, receptor occupancy, JAK activation, or pSTAT3.
- **Paper-safe wording:** Mouse host macrophage/monocyte compartments retain substantial Il10ra/Il10rb receptor-transcript competence after transplantation.
- **Closure status:** `CLOSED_WITH_BOUNDARY`

### L3. Tx -> macrophage state remodeling

- **Source:** v1.4.3C.2.1
- **Biological unit:** Final Clean-B mouse MΦ states; Sham1/Sham20 vs Tx17/Tx5
- **Key validated evidence:** MΦ composition shifted in a treatment-consistent direction: Inflammatory-MΦ ~27.0%->18.1%; Anti-inflammatory-MΦ ~28.9%->33.4%; ECM-associated inflammatory-MΦ ~9.0%->7.9%; Repair/Resolution-MΦ ~2.8%->5.3%. IL10R dual-positive fraction increased across all five MΦ classes.
- **Direction / result:** Tx is associated with reduced inflammatory/profibrotic MΦ abundance and increased anti-inflammatory/repair-associated MΦ abundance.
- **Replicate consistency:** Direction was concordant across Tx17 and Tx5 for the principal composition changes.
- **Evidence strength:** High
- **Causal boundary:** This is comparative single-cell association. The all-MΦ IL10_STAT3 and downstream-no-receptor programs were not increased (approximately -0.746 and -0.610 respectively), so canonical IL10->STAT3 transcriptional activation is not directly demonstrated.
- **Paper-safe wording:** Transplantation was associated with macrophage-state remodeling away from inflammatory/ECM-associated inflammatory states and toward anti-inflammatory and repair/resolution states, while transcript-level evidence of sustained canonical IL10-STAT3 activation was not observed.
- **Closure status:** `CLOSED_WITH_NEGATIVE_BOUNDARY`

### L4. IL10R-positive Repair/Resolution-MΦ -> Pdgfb suppression

- **Source:** v1.4.3C.2.2
- **Biological unit:** Repair/Resolution-MΦ, within-sample IL10R dual-positive vs not-dual-positive
- **Key validated evidence:** C2.1 IL10R fractions reproduced exactly. v6.3.6 Repair/Resolution-MΦ Pdgfb benchmark was reproduced (Sham 5.970686 vs checkpoint 5.971 CP10k; Tx 3.568014 vs checkpoint 3.568). Within Tx, mean per-cell Pdgfb was lower in IL10R-dual-positive cells in Tx17 (-2.543 CP10k) and Tx5 (-0.034), but pseudobulk contrast was Tx17 -2.214 and Tx5 +0.060; Sham also showed negative per-cell contrasts.
- **Direction / result:** No Tx-specific, replicate-consistent direct IL10R-status -> Pdgfb suppression bridge.
- **Replicate consistency:** Not sufficient across metrics; Tx5 pseudobulk direction differs and Sham shows a similar inverse association.
- **Evidence strength:** Negative / limiting evidence
- **Causal boundary:** Do not claim cell-autonomous IL10R activation directly suppresses Pdgfb from this scRNA-seq dataset.
- **Paper-safe wording:** IL10R transcript status did not provide a Tx-specific, replicate-consistent cell-autonomous explanation for Pdgfb suppression; IL10-responsive macrophage remodeling and PDGF attenuation are therefore linked at the network/treatment-phenotype level rather than by a directly demonstrated intracellular chain.
- **Closure status:** `NOT_CLOSED_DIRECTLY`

### L5. Tx -> selective MΦ-to-HSC PDGF attenuation

- **Source:** mouse HSC v6.3.6
- **Biological unit:** Repair/Resolution-MΦ sender and ECM-activated HSC receiver; Sham1/Sham20 vs Tx17/Tx5
- **Key validated evidence:** Repair/Resolution-MΦ Pdgfb per-cell: 5.971->3.568 CP10k. CellChat PDGF communication Repair/Resolution-MΦ -> ECM-activated HSC: ~0.1846->~0.1088. ECM-HSC PDGF mitogenic module: ~+0.479->~-0.479.
- **Direction / result:** Selective attenuation of a PDGF-associated mitogenic/proliferative MΦ-to-HSC signaling arm.
- **Replicate consistency:** Checkpoint classified the axis as strong/coherent Tx-down evidence with biological-sample-aware analyses.
- **Evidence strength:** High
- **Causal boundary:** Not global HSC suppression. FN1 communication increased in Tx, demonstrating pathway-selective rewiring rather than uniform suppression.
- **Paper-safe wording:** The macrophage-HSC interaction analysis supports selective attenuation of a PDGF-associated mitogenic/proliferative signaling arm after transplantation, rather than global suppression of HSC biology.
- **Closure status:** `CLOSED_WITH_PATHWAY_SELECTIVITY`

### L6. PDGF attenuation -> reduced activated/cycling HSC burden

- **Source:** mouse HSC v6.3.6
- **Biological unit:** Mouse HSC states
- **Key validated evidence:** ECM-activated HSC fraction ~0.527->~0.452; cycling fraction within ECM-activated HSC ~0.0560->~0.0478; cycling ECM-HSC burden within HSC3 ~0.0298->~0.0216.
- **Direction / result:** Activated/proliferative HSC burden decreases after Tx, consistent with reduced fibrogenic drive.
- **Replicate consistency:** Strong-to-moderate Tx-down evidence in v6.3.6 checkpoint.
- **Evidence strength:** High
- **Causal boundary:** Association is mechanistically coherent with PDGF attenuation but does not by itself prove PDGF is the sole cause of fibrosis improvement.
- **Paper-safe wording:** Transplantation reduced ECM-activated and cycling HSC burden in parallel with attenuation of the macrophage-to-HSC PDGF mitogenic axis.
- **Closure status:** `CLOSED_WITH_BOUNDARY`

### L7. Mouse-defined IL10-responsive myeloid architecture -> human MASH target conservation

- **Source:** v1.4.3D.0 + v1.4.3D.1
- **Biological unit:** Human LT1-LT8 final myeloid states; single-cell RNA counts and state-level pseudobulk
- **Key validated evidence:** All five final human myeloid states retained IL10RA and IL10RB at state level. In advanced LT1/LT2 (NAS6/F4-F3), same-cell IL10RA+IL10RB dual-positive fractions in PRIMARY states were Resident/Kupffer-like 24.7%, Repair/Resolution 28.3%, and Inflammatory/Transitional 17.9% (LT1/LT2 mean). Both advanced samples contained dual-positive cells in all three PRIMARY states.
- **Direction / result:** Advanced human MASH retains a receptor-transcript-competent myeloid target compartment corresponding to the mouse IL10-responsive architecture.
- **Replicate consistency:** Same-cell dual-positive cells detected in both LT1 and LT2 for all three PRIMARY states.
- **Evidence strength:** High for translational target availability
- **Causal boundary:** Does not show IL10 treatment response in humans; transcript co-detection is not receptor protein, ligand exposure, or pSTAT3. Advanced NAS6 samples are also F3/F4, so disease-severity effects are confounded with fibrosis.
- **Paper-safe wording:** Advanced human MASH retains same-cell IL10RA/IL10RB transcript co-detection across major myeloid states, supporting the biological availability of an IL10-responsive myeloid target compartment in human disease.
- **Closure status:** `CLOSED_TRANSLATIONAL_TARGET_AVAILABILITY`

### L8. Integrated network-level working model

- **Source:** B.2C.3 -> C.1.1 -> C.2.1/C.2.2 -> mouse HSC v6.3.6 -> D.0/D.1
- **Biological unit:** Cross-species integrated mechanism
- **Key validated evidence:** Human graft IL10 transcript evidence; mouse host myeloid IL10R competence; Tx-associated macrophage-state remodeling; selective MΦ-to-ECM-HSC PDGF attenuation; reduced activated/cycling HSC burden; advanced human MASH retains same-cell IL10RA/IL10RB-positive myeloid targets.
- **Direction / result:** Network-level closure supports: engineered graft IL10 / graft secretome -> receptor-competent host myeloid compartment -> macrophage remodeling and selective PDGF-axis rewiring -> reduced fibrogenic HSC activity -> fibrosis improvement; the corresponding IL10-responsive target architecture is retained in human MASH.
- **Replicate consistency:** Major mouse and human links are replicate-aware; direct IL10R->Pdgfb intracellular linkage is the explicit non-closed link.
- **Evidence strength:** High at network level; not a fully causal intracellular chain
- **Causal boundary:** Do not state that scRNA-seq proves IL10 protein secretion, pSTAT3 activation, or direct IL10R-mediated repression of Pdgfb. The supported model is pathway/network-level convergence.
- **Paper-safe wording:** Together, the data support a cross-species model in which the engineered graft engages a receptor-competent host myeloid compartment, is associated with macrophage-state remodeling and selective attenuation of macrophage-to-HSC PDGF signaling, and thereby converges on reduced fibrogenic HSC activity. Human MASH retains the corresponding IL10 receptor-transcript-competent myeloid target compartment.
- **Closure status:** `MECHANISTIC_CLOSURE_AT_NETWORK_LEVEL`

## 4. Quantitative anchor values to preserve

### Human graft IL10 — v1.4.3B.2C.3
- Tx17: 18 IL10 UMI / 13 positive cells in canonical MEX reconstruction.
- Tx5: 4 IL10 UMI / 4 positive cells.
- SfLB1: 2 singleton MEX events; SfLB2: 0.
- Strict all-four shared-TRUE probe sensitivity: SfLB1=0, SfLB2=0, Tx17=5, Tx5=0.
- Do **not** interpret strict shared-probe Tx5=0 as absence of IL10 from Tx5.

### Mouse macrophage remodeling — v1.4.3C.2.1
Approximate mean fractions:
- Inflammatory-MΦ: 27.0% -> 18.1%
- Anti-inflammatory-MΦ: 28.9% -> 33.4%
- ECM-associated inflammatory-MΦ: 9.0% -> 7.9%
- Repair/Resolution-MΦ: 2.8% -> 5.3%
- Lipid-associated/TREM2-MΦ: 32.4% -> 35.3%

IL10R dual-positive fractions increased across all five MΦ classes, but:
- all-MΦ `IL10_STAT3` program: ~-0.746 Tx-Sham
- `IL10_downstream_no_receptors`: ~-0.610
- `Anti_inflammatory`: ~+0.698

Therefore sustained canonical IL10->STAT3 transcriptional activation is **not** directly demonstrated.

### Focused IL10R–PDGFB audit — v1.4.3C.2.2
Repair/Resolution-MΦ:
- v6.3.6 `Pdgfb` benchmark reproduced almost exactly:
  - Sham 5.970686 CP10k vs checkpoint 5.971
  - Tx 3.568014 CP10k vs checkpoint 3.568
- IL10R dual-positive minus not-dual mean per-cell `Pdgfb`:
  - Tx17: -2.5425 CP10k
  - Tx5: -0.0340 CP10k
- Pseudobulk contrast:
  - Tx17: -2.2140 CP10k
  - Tx5: +0.0596 CP10k
- Sham also showed negative per-cell contrasts.
- Final call: **No replicate-consistent direct IL10R-status -> Pdgfb bridge; no clear Tx-specific strengthening.**

### Mouse MΦ–HSC closure — v6.3.6
- Repair/Resolution-MΦ `Pdgfb` per-cell: 5.971 -> 3.568 CP10k
- CellChat PDGF, Repair/Resolution-MΦ -> ECM-activated HSC: ~0.1846 -> ~0.1088
- ECM-HSC PDGF mitogenic module: ~+0.479 -> ~-0.479
- ECM-activated HSC fraction: ~0.527 -> ~0.452
- Cycling fraction within ECM-HSC: ~0.0560 -> ~0.0478
- Cycling ECM-HSC burden within HSC3: ~0.0298 -> ~0.0216
- FN1 communication increased rather than decreased, supporting **selective rewiring rather than global HSC suppression**.

### Human MASH IL10 target conservation — v1.4.3D.0 / D.1
Advanced LT1/LT2 mean same-cell `IL10RA + IL10RB` dual-positive fractions:
- Resident/Kupffer-like MΦ: 24.7%
- Repair/Resolution MΦ: 28.3%
- Inflammatory/Transitional MΦ: 17.9%

Both advanced samples (LT1/F4 and LT2/F3) contained dual-positive cells in all three PRIMARY states.

Do **not** claim receptor enrichment in advanced MASH:
- Resident/Kupffer-like: NAS0 mean 29.8% -> advanced 24.7%
- Repair/Resolution: 29.7% -> 28.3%
- Inflammatory/Transitional: 17.1% -> 17.9%

The supported conclusion is **retention/availability**, not disease-stage induction.

## 5. Claims explicitly permitted

1. The engineered Tx graft contains detectable human IL10 transcript signal.
2. Mouse host macrophage/monocyte populations retain substantial IL10 receptor-transcript competence.
3. Transplantation is associated with macrophage-state remodeling away from inflammatory/ECM-inflammatory states and toward anti-inflammatory/repair-associated states.
4. The macrophage-HSC network shows selective Tx-associated attenuation of a PDGF-associated mitogenic arm.
5. ECM-activated and cycling HSC burden decreases after transplantation.
6. Advanced human MASH retains same-cell `IL10RA`/`IL10RB` transcript co-detection in major myeloid states.
7. The overall cross-species mechanism is closed **at the network level**.

## 6. Claims explicitly not permitted from current data

1. IL10 transcript detection proves secretion of bioactive IL10 protein.
2. scRNA-seq proves receptor occupancy or JAK/STAT3 phosphorylation.
3. Tx causes sustained canonical IL10/STAT3 transcriptional activation across all macrophages.
4. IL10R-positive Repair/Resolution-MΦ directly suppress `Pdgfb` in a cell-autonomous manner.
5. PDGF is the sole mediator of fibrosis improvement.
6. Human MASH patients have been shown to respond therapeutically to IL10.
7. NAS6-associated human changes are pure NAS effects; LT1/LT2 are also F4/F3.

## 7. Recommended manuscript-level wording

> **Together, the data support a cross-species network model in which the engineered human graft engages a receptor-competent host myeloid compartment, is associated with macrophage-state remodeling and selective attenuation of macrophage-to-HSC PDGF mitogenic signaling, and converges on reduced fibrogenic HSC activity. Advanced human MASH retains the corresponding IL10 receptor-transcript-competent myeloid target compartment.**

More conservative version:

> **The data are consistent with a model in which graft-derived IL10 and other secreted factors contribute to host myeloid remodeling, while selective attenuation of macrophage-to-HSC PDGF signaling provides a mechanistically coherent downstream route to reduced HSC activation and fibrosis.**

## 8. Recommended checkpoint filename / Git target

Recommended repo location:

`docs/checkpoints/Human_MASH_IL10_translational_mechanistic_closure_v1.4.3E.0.md`

Recommended annotated tag after review:

`human-mash-il10-translational-mechanistic-closure-v1.4.3E.0`

Do not create the Git tag until the table and wording above have been reviewed against the saved output tables.
