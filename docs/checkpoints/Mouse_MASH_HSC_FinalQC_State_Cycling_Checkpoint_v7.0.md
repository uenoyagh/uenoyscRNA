# Mouse MASH HSC Final QC / State / Cycling Checkpoint v7.0

Date: 2026-11-18  
Project: `uenoyscRNA`  
Scope: Mouse CDAHFD MASH, Sham vs Tx, hepatic stellate cell (HSC) lineage  
Checkpoint type: HSC annotation/QC/state/cycling freeze  
Downstream MΦ→HSC mechanistic closure is intentionally excluded from this checkpoint.

---

## 1. Canonical source and final objects

### Source HSC object

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/HSC_RefinedAnnotation_v6.1.2/RDS/Mouse_MASH_HSC_clean_3state_2state_v6.1.2.rds`

Source HSC cells: **13,166**

Samples:

- Sham1
- Sham20
- Tx17
- Tx5

### QC-audited object

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_MASH_HSC_PrimaryClean_v7.0/objects/Mouse_MASH_HSC_QCAudited_v7.0.rds`

### Canonical Primary Clean object

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_MASH_HSC_PrimaryClean_v7.0/objects/Mouse_MASH_HSC_PrimaryClean_v7.0.rds`

Primary Clean cells: **12,941**

### Canonical Sensitivity Clean object

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_MASH_HSC_PrimaryClean_v7.0/objects/Mouse_MASH_HSC_SensitivityClean_v7.0.rds`

Sensitivity Clean cells: **13,026**

RDS/output data remain on the external SSD and are not committed to GitHub.

---

## 2. Frozen clustering resolution

The HSC clustering resolution was audited across res0.4, 0.6, 0.8, 1.0 and 1.2.

Final decision:

**Freeze HSC clustering at res0.8.**

Key reasons:

- res0.8 exactly reproduces the working HSC clustering (ARI = 1).
- res0.4 is overly coarse.
- res0.6 creates a tiny 4-cell cluster.
- res1.0 and res1.2 increasingly fragment the HSC continuum.
- Low silhouette is accepted because HSC activation is biologically continuous rather than composed of sharply separated discrete states.
- Resolution was not altered to remove contamination.

No reclustering is performed in final QC.

---

## 3. Frozen HSC state framework

### Three-state framework

1. `qHSC`
2. `ECM-activated HSC`
3. `Contractile HSC`

### Two-state framework

1. `qHSC`
2. `aHSC`

The state labels originate from the previously refined v6.1.2 HSC annotation and are retained through the QC freeze. QC does not re-cluster or redefine the biological state framework.

### res0.8 cluster mapping before contaminant exclusion

- Contractile HSC: c11
- ECM-activated HSC: c0, c4, c6, c8, c14, c15
- qHSC: c1, c2, c3, c5, c16

### Effective Primary Clean mapping

- Contractile HSC: c11
- ECM-activated HSC: c0, c4, c6, c8
- qHSC: c1, c2, c3, c5

c14, c15 and c16 are handled by QC as described below.

---

## 4. Final QC disposition

### c14 — HSC/B-cell mixed population

Cells: **96**

Features:

- HSC program remains detectable.
- Coherent B-cell program persists after whole-liver DecontX sensitivity analysis.
- Representative B-cell genes include `Ebf1`, `Cd19`, `Ighm`, `Cd79a`, `Cd79b`, `Pax5`, `Igkc`.

Final disposition:

- Primary Clean: **exclude**
- Sensitivity Clean: **exclude**

Interpretation: HSC–B mixed/doublet-like population; not explained by ambient RNA alone.

### c15 — HSC/LSEC mixed population

Cells: **85**

Features:

- HSC identity remains detectable.
- Strong coherent endothelial/LSEC program persists after ambient correction.
- Representative genes include `Egfl7`, `Clec4g`, `Kdr`, `Pecam1`, `Stab1`, `Stab2`, `Ptprb`, `Esam`, `Aqp1`.

Final disposition:

- Primary Clean: **exclude**
- Sensitivity Clean: **retain**

Interpretation: HSC–LSEC mixed/doublet-like population. Retained in sensitivity analysis to test whether its removal changes HSC conclusions.

### c16 — non-HSC striated/skeletal muscle-like population

Cells: **44**

Features:

- Weak HSC identity.
- Strong intrinsic striated muscle program including `Tnnt3`, `Ttn`, `Neb`, `Actn2`, `Obscn`, `Tcap`, `Atp2a1`.
- Muscle signature persists or strengthens after DecontX.

Final disposition:

- Primary Clean: **exclude**
- Sensitivity Clean: **exclude**

### c11 — bona fide Contractile HSC

Cells before QC: **235**

Retained as genuine Contractile HSC.

---

## 5. Ambient RNA / DecontX policy

Whole-liver DecontX was used **only as a sensitivity analysis**.

Whole-liver source:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds`

HSC overlap with source: **13,166 / 13,166**

Sample-level median contamination estimates were similar across Sham and Tx, arguing against major condition-specific ambient confounding.

Observed whole-HSC effects of DecontX included strong reduction of hepatocyte ambient markers and reduction of endothelial ambient signal, while coherent lineage programs in c14/c15/c16 remained sufficiently strong to support the final QC decisions.

Policy:

- Original filtered RNA counts remain primary.
- DecontX-corrected counts are not used as the canonical expression matrix.
- Whole-liver source objects are not globally overwritten or cleaned.
- Ambient hepatocyte RNA alone is not sufficient for HSC exclusion.

---

## 6. Doublet / mixed-lineage policy

Cells were **not** excluded solely because of:

- high `nCount`
- high `nFeature`
- high complexity
- cell-cycle activity
- cycling marker expression

Candidate mixed/doublet populations require:

1. a coherent foreign-lineage program,
2. retained HSC identity or a mixed bridge phenotype,
3. persistence after ambient-RNA sensitivity assessment where applicable.

Thus the final HSC QC is lineage-program based rather than complexity-based.

---

## 7. Primary Clean integrity

Finalization integrity checks:

- Source HSC cells: **13,166**
- Primary Clean cells: **12,941**
- Expected Primary Clean cells: **12,941**
- Sensitivity Clean cells: **13,026**
- Expected Sensitivity Clean cells: **13,026**
- Working clustering equals frozen res0.8: **TRUE**
- All Primary Clean cells have `qc_primary_include == TRUE`: **TRUE**
- All Sensitivity Clean cells have `qc_sensitivity_include == TRUE`: **TRUE**
- Source RDS overwritten: **FALSE**
- DecontX counts used as primary counts: **FALSE**

QC metadata carried in canonical objects include:

- `qc_original_frozen_cluster`
- `qc_primary_include`
- `qc_sensitivity_include`
- `qc_exclusion_reason`
- `qc_final_disposition`
- `qc_audit_version`
- `qc_resolution_policy`
- `qc_ambient_policy`

---

## 8. Final HSC state behavior in Sham vs Tx

Using Primary Clean HSC:

### ECM-activated HSC

Group-mean fraction:

- Sham: **0.51649**
- Tx: **0.45098**
- Δ Tx − Sham: **−0.06551**

Both Tx samples are below both Sham samples.

Interpretation: **strong Tx-associated reduction of ECM-activated HSC abundance.**

### qHSC

Group-mean fraction:

- Sham: **0.46601**
- Tx: **0.53125**
- Δ Tx − Sham: **+0.06524**

Both Tx samples are above both Sham samples.

Interpretation: **strong Tx-associated shift toward qHSC.**

### Contractile HSC

Group-mean fraction is approximately unchanged and represents a small HSC fraction.

### qHSC/aHSC balance

The qHSC/aHSC ratio increases from approximately **0.88 to 1.13**, supporting a population shift toward HSC deactivation/quiescence.

---

## 9. Cycling audit — strict retention audit

A stringent cycling audit was performed using:

`Mki67`, `Top2a`, `Cenpf`, `Ube2c`, `Ccnb1`, `Ccnb2`

Strict cycling was defined as detection of at least two of these genes in the same cell.

This definition identified only **14 / 13,166 HSC**, making it too stringent for the primary biological cycling analysis.

However, it is useful as a QC retention test:

- strict cycling cells in source: 14
- strict cycling + ECM-high Q75 cells: 6
- strict cycling + ECM-high Q90 cells: 2
- strict cycling + ECM-high Q75 removed by Primary Clean: **0**
- strict cycling + ECM-high Q90 removed by Primary Clean: **0**

Conclusion:

**Primary Clean QC did not selectively remove the strict cycling + ECM-high HSC population.**

---

## 10. Legacy-compatible Cycling definition

The previous v6.3.5 analysis used a broader and biologically more sensitive cycling-like definition.

Core genes:

`Mki67`, `Top2a`, `Pcna`, `Mcm2`, `Mcm3`, `Mcm4`, `Mcm5`, `Mcm6`, `Mcm7`, `Ccna2`, `Ccnb1`, `Cdk1`

Legacy cycling-like definition:

**Mki67 positive OR Top2a positive OR at least two core proliferation genes detected.**

This definition was reapplied without modification to Source, Primary Clean and Sensitivity Clean HSC.

### Exact reproducibility

The Source all-HSC object **exactly reproduces the v6.3.5 cycling results** for all four samples.

Source legacy cycling-like ECM-HSC counts:

- Sham1: 87 / 1,847
- Sham20: 99 / 1,527
- Tx17: 50 / 1,191
- Tx5: 98 / 1,826

Total legacy cycling-like ECM-HSC: **334**

---

## 11. Legacy Cycling results after QC

### Primary Clean

ECM-activated HSC fraction:

- Sham mean: **0.51649**
- Tx mean: **0.45098**
- Δ: **−0.06551**
- 4/4 pairwise comparisons: Tx lower
- direction: **Tx_down**

Cycling-like fraction within ECM-activated HSC:

- Sham mean: **0.04954**
- Tx mean: **0.04674**
- Δ: **−0.00280**
- 3/4 pairwise comparisons: Tx lower
- direction: **Tx_down**

Cycling ECM-HSC burden within all HSC:

- Sham mean: **0.02569**
- Tx mean: **0.02105**
- Δ: **−0.00465**
- relative reduction: approximately **18%**
- 3/4 pairwise comparisons: Tx lower
- direction: **Tx_down**

### Sensitivity Clean

Cycling-like fraction within ECM-activated HSC:

- Sham mean: **0.05111**
- Tx mean: **0.04786**
- Δ: **−0.00325**
- direction: **Tx_down**

Cycling ECM-HSC burden:

- Sham mean: **0.02665**
- Tx mean: **0.02172**
- Δ: **−0.00494**
- relative reduction: approximately **19%**
- direction: **Tx_down**

### QC loss of legacy cycling-like ECM-HSC

Of 334 source legacy cycling-like ECM-HSC:

- Primary Clean removes 33 cells (**9.88%**)
- Sensitivity Clean removes 20 cells (**5.99%**)

Despite these removals, both Primary and Sensitivity analyses preserve the same Tx-down direction.

Conclusion:

**The previously observed reduction of proliferative/cycling ECM-HSC is not an artifact of the final doublet/mixed-lineage QC.**

---

## 12. Interpretation of HSC fibrogenic burden

The HSC results should not be summarized as a uniform per-cell shutdown of ECM transcription.

The most defensible interpretation is:

1. Tx reduces the **abundance of ECM-activated HSC**.
2. Tx increases the **qHSC fraction and quiescence program**.
3. Within ECM-activated HSC, the legacy-compatible **cycling-like fraction decreases modestly**.
4. Consequently, the total **cycling ECM-HSC burden decreases**.
5. Residual HSC can retain or increase selected wound-healing / ECM-remodeling transcripts such as `Fn1` and `Serpine1`.
6. Therefore, the treatment effect is best described as **HSC state redistribution / partial deactivation with reduced proliferative fibrogenic cellular burden**, rather than uniform suppression of every ECM transcript in every HSC.

This distinction should be maintained in downstream manuscript and mechanistic figures.

---

## 13. Scripts included in this checkpoint

The following successful HSC scripts define this checkpoint:

- `analysis/Mouse_HSC_FinalQC_Preflight_v1.0.R`
- `analysis/Mouse_HSC_Resolution_Provenance_Audit_v1.0.1.R`
- `analysis/Mouse_HSC_Cluster_QC_Audit_res0.8_v1.0.R`
- `analysis/Mouse_HSC_WholeLiver_DecontX_Input_Preflight_v1.0.R`
- `analysis/Mouse_HSC_WholeLiver_DecontX_Sensitivity_v1.0.R`
- `analysis/Mouse_HSC_PrimaryClean_Finalization_v1.0.1.R`
- `analysis/Mouse_HSC_Cycling_ECM_Retention_Audit_v1.0.R`
- `analysis/Mouse_HSC_LegacyCycling_Compatibility_Audit_v1.0.R`

Downstream functional-response, ECM burden, ligand–receptor and cross-lineage mechanistic scripts are intentionally not part of this HSC state/QC freeze checkpoint.

---

## 14. Frozen downstream policy

For subsequent mouse MASH HSC analyses:

- Primary biological analyses use:
  `Mouse_MASH_HSC_PrimaryClean_v7.0.rds`
- Sensitivity analyses use:
  `Mouse_MASH_HSC_SensitivityClean_v7.0.rds`
- Frozen clustering resolution: **res0.8**
- Primary HSC state framework:
  `qHSC / ECM-activated HSC / Contractile HSC`
- Statistical replicate: **biological sample**
- Cell-level inferential p-values are not used as biological-replicate inference.
- Primary expression analyses use original filtered counts.
- c15 is the key mixed-lineage sensitivity population.
- Cycling/proliferation analyses intended to compare with v6.3.5 use the **legacy-compatible definition**, not the overly stringent six-gene ≥2-hit definition.

---

## 15. Checkpoint conclusion

The mouse HSC lineage is sufficiently audited to freeze its QC and state framework for downstream mechanistic analyses.

The principal Tx-associated HSC phenotype is:

**ECM-activated HSC ↓ + qHSC ↑ + cycling ECM-HSC burden ↓**

This phenotype remains directionally stable after Primary Clean QC and after restoring c15 in Sensitivity Clean.

The HSC QC does not support the concern that proliferative ECM-producing HSC were selectively removed by aggressive doublet filtering.

Next checkpoint should integrate the frozen HSC state framework with the already cleaned macrophage lineage to evaluate the MΦ→PDGF→HSC and matrix-remodeling mechanism.
