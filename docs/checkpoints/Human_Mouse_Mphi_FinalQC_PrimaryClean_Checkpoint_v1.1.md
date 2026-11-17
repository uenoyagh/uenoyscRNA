# Human / Mouse MΦ Final QC and Primary-Clean Checkpoint v1.1

Date: 2026-11-18

## 1. Scope — this checkpoint applies to the macrophage fraction only

**The Primary Clean RDS files created in this work contain only the macrophage / myeloid fraction that had already been extracted for the MΦ analysis.**

This work did **not** create a globally cleaned whole-liver RDS.

The whole-liver objects were used only as a biologically appropriate background for the ambient-RNA sensitivity analysis:

`whole-liver filtered counts`
→ `sample-wise DecontX using broad lineage labels`
→ `extract the already-frozen MΦ barcodes`
→ `compare original vs DecontX-corrected MΦ expression`
→ `make QC inclusion/exclusion decisions`
→ `build Primary Clean MΦ RDS`

Therefore:

- whole-liver hepatocytes, HSCs, LSECs, cholangiocytes, lymphocytes, etc. were **not** globally cleaned or re-annotated by this workflow;
- the whole-liver source RDS files were not overwritten;
- DecontX-corrected counts were used as a **sensitivity analysis** and were not adopted as the primary expression matrix;
- the final canonical outputs of this checkpoint are MΦ-only Primary Clean RDS files.

## 2. Purpose

Before downstream mechanistic analysis, the Human and Mouse liver MΦ datasets were audited for:

1. clustering-resolution validity;
2. annotation / identity validity;
3. non-MΦ contamination;
4. ambient RNA / passenger RNA;
5. heterotypic doublets or mixed-cell populations;
6. robustness of MΦ functional programs after ambient-RNA correction.

The objective was to freeze an auditable and reproducible MΦ analytical foundation before proceeding to M1/M2, Resident/Recruited/LAM, IL10-response, fibrosis-related, and cell-cell interaction analyses.

## 3. Frozen clustering solutions

| Species | Frozen resolution | Policy |
|---|---:|---|
| Human MΦ | res0.8 | Retain; no reclustering |
| Mouse MΦ | res2.0 | Retain; no reclustering |

Resolution selection was based on clustering provenance, stability, marker coherence, sample reproducibility, and neighborhood structure, not on desired disease or treatment direction.

## 4. Ambient-RNA strategy

Raw / empty-droplet background matrices were not available on the currently mounted storage.

Accordingly, ambient-RNA assessment was performed as a sensitivity analysis using the whole-liver filtered-count objects:

- each sample was processed separately;
- broad lineage labels were supplied to DecontX;
- MΦ functional states such as M1/M2/LAM were not used as contamination-source labels;
- after correction, only the frozen MΦ barcodes were extracted and compared with the original counts.

**Primary expression policy:** original filtered counts remain primary.

**Sensitivity policy:** DecontX-corrected MΦ counts are retained only for robustness checks.

## 5. Final cluster decisions

### Human

- **Cluster 13 — EXCLUDE from primary MΦ analysis**
  - coherent neutrophil / non-MΦ program;
  - CSF3R, MPO, ELANE, OLFM4 positive;
  - weak macrophage-core and phagocytosis programs;
  - not explained by ambient RNA.

- **Cluster 3 — KEEP**
  - strong hepatocyte-derived signal before correction;
  - hepatocyte signal decreases markedly after DecontX;
  - macrophage, resident/Kupffer-like, M2 and phagocytosis programs remain;
  - interpreted as bona fide macrophages with ambient/passenger hepatocyte RNA.

- **Cluster 16 — KEEP**
  - cycling MΦ population.

### Mouse

- **Cluster 23 — EXCLUDE from primary MΦ analysis**
  - strong macrophage program plus coherent endothelial/LSEC program;
  - high-complexity MΦ–LSEC heterotypic doublet / mixed-cell candidate;
  - excluded from primary MΦ biology;
  - retained in flagged object for sensitivity analyses.

- **Cluster 8 — KEEP**
  - hepatocyte/HSC signal decreases strongly after DecontX;
  - macrophage, resident and phagocytosis programs persist.

- **Cluster 15 — KEEP**
  - hepatocyte signal decreases strongly after DecontX;
  - macrophage identity persists.

- **Cluster 26 — KEEP**
  - inflammatory M1/IFN-high macrophage state remains robust after correction.

## 6. Primary Clean RDS files

### Human

Source frozen MΦ set: 6,720 cells  
Excluded: cluster 13, 90 cells  
Primary Clean: **6,630 cells**  
Retained fraction: **98.66%**

Canonical RDS:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_Mphi_PrimaryClean_v2.0/objects/Human_MASH_Mphi_PrimaryClean_v2.0.rds`

Flagged source-derived object:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_Mphi_PrimaryClean_v2.0/objects/Human_MASH_Mphi_FinalQC_FLAGGED_v2.0.rds`

### Mouse

Source frozen MΦ set: 18,029 cells  
Excluded: cluster 23, 173 cells  
Primary Clean: **17,856 cells**  
Retained fraction: **99.04%**

Canonical RDS:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mouse_Mphi_PrimaryClean_v5.0/objects/Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`

Flagged source-derived object:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_Mphi_RDS/Mouse_Mphi_PrimaryClean_v5.0/objects/Mouse_MASH_Mphi_FinalQC_FLAGGED_v5.0.rds`

## 7. QC metadata added to the finalized objects

- `qc_original_frozen_cluster`
- `qc_primary_include`
- `qc_sensitivity_include`
- `qc_exclusion_reason`
- `qc_final_disposition`
- `qc_audit_version`
- `qc_resolution_policy`
- `qc_ambient_policy`

Mouse cluster 23:
- `qc_primary_include = FALSE`
- `qc_sensitivity_include = TRUE`

Human cluster 13:
- `qc_primary_include = FALSE`
- `qc_sensitivity_include = FALSE`

## 8. Canonical finalization outputs

Cross-species finalization directory:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/CrossSpecies_MASH_Mphi/PrimaryClean_Mphi_QC_Finalization_v1.0/`

Key tables:

- `01_final_cluster_disposition_v1.0.tsv`
- `02_cluster_counts_before_after_primary_QC_v1.0.tsv`
- `03_sample_counts_before_after_primary_QC_v1.0.tsv`
- `04_sample_cluster_counts_before_after_primary_QC_v1.0.tsv`
- `05_primary_clean_summary_v1.0.tsv`

Ambient / DecontX audit directories are preserved separately and are not the canonical primary expression inputs.

## 9. Git provenance

Final QC scripts and checkpoint documentation were committed and pushed.

- commit: `182cf29`
- tag: `human-mouse-mphi-final-qc-primary-clean-v1.0.1`

The earlier tag `human-mouse-mphi-final-qc-primary-clean-v1.0` is retained for provenance.

## 10. Downstream analysis policy

Primary downstream analyses must use:

- Human: `Human_MASH_Mphi_PrimaryClean_v2.0.rds`
- Mouse: `Mouse_MASH_Mphi_PrimaryClean_v5.0.rds`

Planned priorities:

1. M1 and M2 functional axes;
2. M2/M1 balance;
3. Resident / Recruited / LAM composition;
4. Mouse Sham vs Tx as the primary treatment comparison;
5. Mouse STD vs CDAHFD disease comparison;
6. IL10-response program;
7. profibrotic macrophage program;
8. MΦ–HSC / MΦ–Hepatocyte / MΦ–LSEC interaction analyses.

M1/M2 are treated as **functional axes**, not as fixed cell identities.

For Mouse cluster 23, primary analyses exclude the cluster, while selected sensitivity analyses may restore it.

## 11. Important limitation / future extension

This checkpoint does **not** certify the non-MΦ whole-liver compartments as globally cleaned.

If future analyses require a fully audited whole-liver atlas, separate lineage-specific QC should be performed for hepatocytes, HSC/mesenchymal cells, LSEC/endothelial cells, cholangiocytes, lymphoid cells, and other compartments.

This MΦ checkpoint should remain unchanged as the frozen macrophage analysis foundation unless a later audit identifies a material error.

