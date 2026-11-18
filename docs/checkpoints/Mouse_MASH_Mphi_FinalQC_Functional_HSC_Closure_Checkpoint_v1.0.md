# Mouse MASH Macrophage Final QC / Functional State / HSC Closure Checkpoint v1.0

Date: 2026-11-18  
Project: `uenoyscRNA`  
Scope: Mouse CDAHFD MASH, Sham vs Tx, macrophage lineage and macrophage–HSC biological consistency  
Checkpoint type: Mouse MΦ final QC + functional subcluster/state + HSC consistency freeze

This checkpoint builds on the previously frozen macrophage QC checkpoint:

- Git commit: `6cbe4d3`
- Git tag: `human-mouse-mphi-final-qc-primary-clean-v1.1`

and the frozen HSC checkpoint:

- Git commit: `cce1edf`
- Git tag: `mouse-mash-hsc-final-qc-state-cycling-v7.0`

The purpose of this checkpoint is to freeze the interpretation of the **Mouse macrophage lineage** for downstream mechanistic analyses.

---

## 1. Canonical Mouse macrophage source and Primary Clean object

### Final annotated Mouse macrophage source

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mphi_Res2_CleanB_FINAL_v4.14.5/RDS/Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds`

Cells: **18,029**

Frozen macrophage clustering resolution:

**res2.0**

### Canonical Primary Clean macrophage object

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mouse_Mphi_PrimaryClean_v5.0/objects/Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`

Primary Clean cells: **17,856**

Primary analyses should use this object unless a specific sensitivity analysis requires otherwise.

RDS and SSD output files are not committed to GitHub.

---

## 2. res2.0 clustering freeze

The macrophage clustering resolution has already been independently audited and frozen.

Final decision:

**Freeze Mouse macrophage clustering at res2.0.**

Key policy:

- resolution is not adjusted to remove contamination;
- macrophage QC is performed after the frozen clustering is defined;
- high complexity alone is not a reason for exclusion;
- ambient RNA alone is not a reason for exclusion;
- coherent foreign-lineage programs are required to support a mixed-lineage/doublet interpretation.

No downstream functional analysis in this checkpoint performs reclustering.

---

## 3. Final macrophage QC status

### c8

Retained.

Interpretation:

- macrophage identity remains coherent;
- apparent foreign/ambient signal is compatible with phagocytic/ambient biology;
- no sufficient evidence for exclusion after sensitivity review.

### c15

Retained.

Interpretation:

- macrophage identity remains coherent;
- ambient/phagocytic signal does not justify exclusion;
- retained in the canonical Primary Clean object.

### c23

Cells: **173**

Final disposition:

**exclude from Primary Clean**

Reason:

- coherent LSEC/endothelial program;
- high complexity / bridge phenotype;
- foreign-lineage program persists sufficiently after ambient-RNA sensitivity assessment;
- interpreted as LSEC-associated mixed/doublet-like population rather than a clean macrophage state.

### c26

Cells: **95**

Final disposition:

**retain**

Interpretation:

- coherent macrophage identity;
- M1/IFN-high phenotype;
- strong sample/condition skew toward CDAHFD requires caution;
- condition/library confounding remains a caveat;
- sample skew alone is not a sufficient exclusion criterion.

---

## 4. Doublet / ambient-RNA policy

The final macrophage QC does **not** use aggressive global doublet removal.

Cells are not excluded solely because of:

- high `nCount`
- high `nFeature`
- high RNA complexity
- inflammatory activation
- cycling state
- sample specificity

Candidate exclusion requires a coherent foreign-lineage program and biologically plausible mixed-lineage structure.

Whole-liver DecontX is used as a **sensitivity analysis**, not as the primary expression matrix.

Primary expression analyses use original filtered counts.

---

## 5. Primary functional macrophage framework

The macrophage states are interpreted as **continuous functional programs**, not rigid cell identities.

The key functional axes retained for downstream interpretation are:

### Inflammatory / M1-like programs

- `M1_core`
- `M1_IFN`
- `M1_MASH_pathogenic`

### Resolution / anti-inflammatory programs

- `M2_core`
- `M2_resolution/efferocytosis`
- anti-inflammatory / resolution-repair program

### Fibrosis / ECM-associated macrophage programs

- fibrosis-associated macrophage program
- ECM remodeling program
- LAM/TREM2-like program

### Other biologically relevant programs

- monocyte/recruited identity
- resident/Kupffer-like identity
- efferocytosis/phagocytosis
- matrix degradation
- phagolysosomal remodeling

The labels `M1` and `M2` are treated as functional axes rather than mutually exclusive terminal identities.

---

## 6. Sham vs Tx macrophage functional composition

Using the canonical Mouse Primary Clean macrophage object:

### Resolution / anti-inflammatory shift

- `M2_core`: approximately **+14.4%**
- `M2_resolution/efferocytosis`: approximately **+8.2%**
- `M2_minus_M1_core`: approximately **+12.6%**
- Anti-inflammatory class: approximately **+3.64 percentage points**
- Repair/Resolution class: approximately **+2.01 percentage points**

Both Tx samples are above both Sham samples for the strongest resolution-associated metrics.

### Inflammatory behavior

- `M1_core`: approximately **+18.2%**
- `M1_IFN`: approximately **−3.8%**
- `M1_MASH_pathogenic`: approximately **−12.8%**

The key interpretation is not a simple global M1→M2 conversion.

Instead, Tx is associated with:

**reduction of the MASH-pathogenic inflammatory program together with increased resolution/efferocytosis and remodeling programs.**

### Broad class composition

Approximate group shifts:

- Inflammatory: **−6.77 percentage points**
- ECM-associated inflammatory: **−0.81 percentage points**
- Anti-inflammatory: **+3.64 percentage points**
- Repair/Resolution: **+2.01 percentage points**
- LAM/TREM2-like: **+2.53 percentage points**

---

## 7. M2/M1 ratio interpretation

The pseudobulk M2/M1 ratio does **not** increase consistently.

Approximate sample values:

- Sham1: 3.64
- Sham20: 4.30
- Tx17: 4.28
- Tx5: 3.14

Group mean:

- Sham: 3.97
- Tx: 3.71

Therefore:

**Do not claim that Tx globally increases the M2/M1 ratio.**

The supported interpretation is a selective redistribution of inflammatory, resolution, efferocytic and remodeling programs.

---

## 8. IL10 / fibrosis-associated macrophage interpretation

Global IL10-response is not uniformly increased.

However:

- anti-inflammatory / resolution program increases;
- fibrosis-associated macrophage program increases;
- ECM-remodeling program increases;
- MASH-pathogenic inflammatory program decreases;
- efferocytosis increases.

The fibrosis-associated increase is not interpreted as simple profibrotic activation.

Gene-level decomposition shows that the increase is driven mainly by remodeling / repair-associated genes rather than canonical direct HSC-activating ligands.

Representative behavior:

- `Pdgfb`: decreases
- `Spp1`: decreases
- `Tgfb1`: decreases
- `Pdgfa`: modestly increases
- `Lgals3`: increases
- `Fn1`: increases
- `Gpnmb`: increases
- `Mmp12`: increases
- `Gdf15`: strongly increases
- `Il1rn`, `Areg`, `Apoe`, `Igf1`: increase
- `Mertk`: increases

Therefore:

**Tx macrophages are better described as matrix-turnover / phagolysosomal-remodeling / repair-associated than as uniformly profibrotic.**

---

## 9. PDGF-axis correction

A previous working hypothesis emphasized:

`macrophage PDGFA ↓ → HSC fibrogenesis ↓`

The strict Primary Clean analysis does **not** support this as the principal PDGF mechanism.

Current result:

- `Pdgfa`: modestly increases
- `Pdgfa–Pdgfra` compatibility: approximately flat
- `Pdgfb`: decreases
- HSC `Pdgfra`: decreases
- `Pdgfb–Pdgfra` compatibility decreases

Approximate ligand–receptor compatibility changes:

- all-MΦ → qHSC `Pdgfb–Pdgfra`: log2 ratio approximately **−0.33**
- all-MΦ → ECM-HSC `Pdgfb–Pdgfra`: log2 ratio approximately **−0.30**

Therefore, if the PDGF axis is retained mechanistically, the supported formulation is:

**MΦ PDGFB ↓ + HSC PDGFRA ↓**

rather than PDGFA suppression.

---

## 10. Independent HSC receiver consistency

The HSC lineage was independently frozen at:

- `res0.8`
- Primary Clean: **12,941 cells**
- Sensitivity Clean: **13,026 cells**

See:

`mouse-mash-hsc-final-qc-state-cycling-v7.0`

The following HSC findings are independent receiver-side observations.

### HSC state composition

Primary Clean:

- ECM-activated HSC: **51.65% → 45.10%**
- qHSC: **46.60% → 53.13%**
- qHSC/aHSC ratio: approximately **0.88 → 1.13**

Both Tx samples are below both Sham samples for ECM-activated HSC and above both Sham samples for qHSC.

### HSC cycling / proliferative fibrogenic burden

Using the legacy-compatible v6.3.5 cycling definition:

`Mki67+ OR Top2a+ OR >=2 core proliferation genes detected`

Primary Clean:

- cycling-like fraction within ECM-activated HSC: Tx-down
- cycling ECM-HSC burden within all HSC: Tx-down
- cycling ECM-HSC burden relative reduction: approximately **18%**

The direction is preserved in Sensitivity Clean.

Thus, HSC receiver biology is consistent with reduced activated/proliferative fibrogenic cellular burden.

---

## 11. HSC ECM transcription is not globally suppressed

The receiver-side HSC data do **not** support a claim that Tx uniformly suppresses ECM transcription in every HSC.

Primary Clean total-HSC behavior includes:

- collagen-core expression: approximately unchanged to modestly increased
- ECM structural program: approximately unchanged to modestly increased
- matrix-remodeling program: increased
- fibrogenic stress program: increased
- contractile program: decreased
- quiescence program: increased

Representative genes:

- `Col1a1`: approximately unchanged/slightly increased
- `Col1a2`: near flat
- `Col3a1`: modestly increased
- `Col5a1`: increased
- `Fn1`: increased
- `Serpine1`: increased
- `Acta2`: decreased
- `Tagln`: decreased
- `Lrat`: increased
- `Reln`: increased
- `Colec11`: increased
- `Vipr1`: increased

Therefore the supported HSC interpretation is:

**population-level deactivation/state redistribution with reduced activated/proliferative HSC burden, while residual HSC retain wound-healing / ECM-remodeling transcription.**

---

## 12. Cross-lineage MΦ–HSC consistency

The cross-lineage analysis was performed as a descriptive sample-level mechanistic closure.

### Macrophage side

Tx-associated changes include:

- matrix degradation program: approximately **+58%** by pseudobulk sum
- phagolysosomal remodeling: approximately **+25%**
- resolution/repair: increased
- matrix adhesion/remodeling: increased
- direct HSC-activating ligand panel: approximately flat to mildly decreased
- `Pdgfb`: decreased

### HSC side

Tx-associated changes include:

- qHSC/quiescence: increased
- ECM-activated HSC fraction: decreased
- cycling ECM-HSC burden: decreased
- contractile program: decreased
- HSC matrix degradation: not increased
- collagen synthesis: not globally decreased

### Integrated interpretation

The most defensible integrated model is:

**Tx shifts the macrophage compartment away from MASH-pathogenic inflammatory signaling toward resolution/efferocytosis and matrix-turnover/remodeling programs, while the HSC compartment shifts from ECM-activated/proliferative states toward qHSC.**

The PDGF-compatible component is:

**MΦ PDGFB ↓ + HSC PDGFRA ↓ → reduced HSC activation/proliferation compatibility**

This is consistent with the independent reduction of ECM-activated and cycling ECM-HSC burden.

---

## 13. What this cross-lineage analysis does and does not prove

The MΦ–HSC relationship provides **biological consistency**, not a standalone proof of macrophage subcluster identity.

Macrophage subcluster/state validity is supported by multiple independent layers:

1. frozen res2.0 clustering;
2. lineage-marker coherence;
3. ambient-RNA and mixed-lineage QC;
4. Primary Clean stability;
5. sample-level functional composition;
6. gene-level decomposition;
7. biologically coherent receiver-side HSC behavior.

Therefore:

**No major QC or biological inconsistency currently requires re-clustering or rejection of the Mouse macrophage Primary Clean functional state framework.**

However:

- M1/M2 labels remain functional continua;
- c26 remains condition/library-confounded and should be interpreted cautiously;
- cross-lineage transcriptomic compatibility is not causal proof;
- cell-level p-values must not substitute for biological-replicate inference.

---

## 14. Statistical interpretation policy

Biological replicate:

**sample**

Samples:

- Sham1
- Sham20
- Tx17
- Tx5

n = 2 biological samples per group.

Therefore:

- emphasize effect size;
- emphasize direction across biological replicates;
- use pairwise direction consistency;
- avoid over-interpreting formal p-values;
- do not use cells as independent biological replicates.

---

## 15. Final Mouse macrophage conclusion

The Mouse macrophage lineage is considered sufficiently audited for downstream mechanistic analyses.

Current frozen conclusion:

**The res2.0 / Primary Clean macrophage framework is acceptable and does not show a major subcluster-quality problem requiring re-clustering.**

The principal Tx-associated macrophage phenotype is:

**MASH-pathogenic inflammatory program ↓  
+ resolution/efferocytosis ↑  
+ matrix degradation / phagolysosomal remodeling ↑  
+ PDGFB ↓**

This macrophage state shift is biologically consistent with the independently observed HSC phenotype:

**ECM-activated HSC ↓  
+ qHSC ↑  
+ cycling ECM-HSC burden ↓**

The integrated interpretation is:

**macrophage functional rewiring toward resolution/remodeling together with reduced PDGFB–PDGFRA activation compatibility and HSC deactivation/state redistribution.**

---

## 16. Scripts included in this checkpoint

This checkpoint should include the successful Mouse Primary Clean functional and MΦ–HSC consistency scripts:

- `analysis/Mouse_PrimaryClean_ShamsTx_Preflight_v1.0.1.R`
- `analysis/Mouse_PrimaryClean_ShamsTx_FunctionalComposition_v1.0.R`
- `analysis/Mouse_PrimaryClean_ShamsTx_IL10_Profibrotic_Attribution_v1.0.R`
- `analysis/Mouse_PrimaryClean_ShamsTx_GeneLevel_Decomposition_v1.0.R`
- `analysis/Mouse_PrimaryClean_Mphi_to_HSC_LR_Compatibility_v1.0.R`
- `analysis/Mouse_PrimaryClean_HSC_ShamsTx_FunctionalResponse_v1.0.R`
- `analysis/Mouse_PrimaryClean_HSC_TotalFibrogenicBurden_v1.0.R`
- `analysis/Mouse_PrimaryClean_HSC_ECM_Synthesis_vs_Resolution_v1.0.R`
- `analysis/Mouse_CrossLineage_Mphi_HSC_ECM_MechanisticClosure_v1.0.R`

The earlier macrophage QC/freeze scripts are already captured by:

`human-mouse-mphi-final-qc-primary-clean-v1.1`

The HSC QC/freeze scripts are already captured by:

`mouse-mash-hsc-final-qc-state-cycling-v7.0`

This checkpoint therefore adds the Mouse macrophage functional interpretation and independent HSC receiver consistency layer without duplicating the earlier frozen QC checkpoints.

---

## 17. Frozen downstream policy

For subsequent Mouse MASH macrophage analyses:

- canonical MΦ object:
  `Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`
- frozen MΦ resolution:
  **res2.0**
- do not re-cluster merely to improve a downstream hypothesis;
- treat M1/M2 as functional axes;
- c23 remains excluded;
- c26 remains retained with condition/library-confounding caveat;
- use original filtered counts for primary expression analyses;
- use biological sample as replicate;
- use HSC Primary Clean v7.0 as the canonical HSC receiver object;
- preserve Primary/Sensitivity comparison where c15 HSC could matter;
- use `PDGFB–PDGFRA`, not `PDGFA↓`, when describing the PDGF-compatible mechanism;
- do not claim global HSC ECM-transcription suppression;
- do not claim causal macrophage-mediated ECM clearance from transcriptomics alone.

---

## 18. Next step after this checkpoint

After this checkpoint is frozen, the Mouse MΦ and HSC foundations should remain unchanged unless a specific new QC failure is demonstrated.

Subsequent work can focus on:

- formal composition-vs-within-state decomposition of HSC ECM burden;
- macrophage–HSC ligand/receptor prioritization;
- integration with hepatocyte/LSEC axes;
- figure-level mechanistic model construction;
- manuscript-ready claim hierarchy.

