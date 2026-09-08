# Human MASH iPSC IL10 probe audit checkpoint v1.4.2B.3

Date: 2026-09-08

## Repository state at checkpoint creation

- Repository: /Users/uenoya/Projects/uenoyscRNA
- Branch: main
- Parent HEAD: 220e5e0df9e9feed70f745cdd2350d9d9fe82f8a
- Parent HEAD short: 220e5e0

## Biological analysis version

- Biological analysis version: v1.4.2B.3
- Validated integrated implementation: v1.4.2B.3.2
- Post-audit technical implementation: v1.4.2B.3.3.1

## Successful scripts included in this checkpoint

1. scripts/Human_MASH/Human_MASH_v1.4.2B.3.2_IL10_ProbeIntegratedAudit.sh
   - SHA256: 816285a3120620f4ffe58b0d19ce51ab2c2db7d74a63d73c186185a5a6d3bc09
   - Purpose:
     - reconstruct IL10 probe-level counts in SfLB1/SfLB2 and Tx17/Tx5
     - apply sample-specific Cell Ranger filtered_probes status
     - reconcile probe-level IL10 with current Seurat RDS gene-level IL10
     - generate sample-, cell-, lineage-, and probe-level audit outputs

2. scripts/Human_MASH/Human_MASH_v1.4.2B.3.3.1_IL10_PostAuditRepair_HarmonizedSensitivity.sh
   - SHA256: ad235e98764fbff465bba43c06fc15ab4250bae7e2b63711dc92e374647ced23
   - Purpose:
     - repair Fig03 cell-label collision using full RDS cell IDs as unique factor identities
     - retain shortened barcodes only as display labels
     - add strict harmonized sensitivity analysis using IL10 probes accepted in all four samples

## Source provenance

### SfLB

- Delivery: RH241011956
- Samples:
  - SfLB_1 -> SfLB1
  - SfLB_2 -> SfLB2
- Biological status: transgene-negative iPSC-LO controls
- Flex probe panel:
  - Chromium Human Transcriptome Probe Set v1.0.1

### Tx

- Delivery: RH251217117
- Samples:
  - Tx17
  - Tx5
- Biological status: transplanted tissues derived from IL10-transgene iPSC-LO
- Flex probe panel:
  - Chromium Human Transcriptome Probe Set v2.0.0

## IL10 probe identities

Three IL10 probe IDs were traced directly through the Cell Ranger outputs:

- 6b6892b
- b1bd8e2
- e99f652

Sample-specific Cell Ranger filtered_probes status:

| Probe | SfLB1 | SfLB2 | Tx17 | Tx5 |
|---|---:|---:|---:|---:|
| 6b6892b | FALSE | FALSE | TRUE | TRUE |
| b1bd8e2 | TRUE | TRUE | TRUE | TRUE |
| e99f652 | FALSE | FALSE | FALSE | FALSE |

Thus b1bd8e2 is the only IL10 probe accepted across all four samples.

## Primary sample-specific accepted-probe result

Cell Ranger gene-level IL10 and reconstructed sample-specific accepted-probe IL10 were exactly concordant:

| Sample | RDS cells | Gene-level IL10 UMI | Accepted-probe IL10 UMI | IL10+ cells |
|---|---:|---:|---:|---:|
| SfLB1 | 19,081 | 0 | 0 | 0 |
| SfLB2 | 9,805 | 0 | 0 | 0 |
| Tx17 | 5,958 | 18 | 18 | 13 |
| Tx5 | 8,571 | 4 | 4 | 4 |

Validation:

- exact_total_match = TRUE for all four samples
- n_cellwise_mismatches = 0 for all four samples
- cellwise gene-vs-accepted-probe concordance = 1.0000

## SfLB negative-control audit

SfLB1 had two raw single-UMI IL10-probe events in current RDS cells:

- Endothelial_like:
  - 6b6892b = 1 UMI
  - filtered_probes = FALSE

- Hepatobiliary_like:
  - e99f652 = 1 UMI
  - filtered_probes = FALSE

Both were excluded by sample-specific Cell Ranger probe filtering.

Therefore:

- SfLB1 accepted IL10 = 0
- SfLB2 accepted IL10 = 0
- SfLB1/SfLB2 gene-level IL10 = 0

## Tx17 strongest IL10 cell

Strongest IL10-positive cell:

- Cell:
  Tx17_Tx17_ATGGGTAGTACACCTGCAACTCAAAC_1
- Broad lineage:
  Mesenchymal_like
- State:
  Perivascular_mesenchymal_HSC_like
- Gene-level IL10:
  5 UMI
- Probe support:
  - 6b6892b = 3 UMI
  - b1bd8e2 = 2 UMI
  - e99f652 = 0 UMI
- Independent accepted IL10 probes:
  2

This is the only cell supported by two independent accepted IL10 probes.

A second mesenchymal cell had:

- Cycling_mesenchymal_HSC_like
- IL10 = 2 UMI
- 6b6892b = 2 UMI

All non-mesenchymal IL10-positive cells were single-UMI cells.

## Primary Tx17 lineage distribution

Accepted IL10 UMI in Tx17:

- Mesenchymal_like = 8
- Hepatocyte_like = 7
- Unresolved = 2
- Myeloid_like = 1

## Strict harmonized sensitivity analysis

To eliminate differences in sample-specific probe filtering between SfLB v1.0.1 and Tx v2.0.0, a strict analysis was performed using only the probe accepted in all four samples:

- b1bd8e2
- region: spliced

Results:

| Sample | Harmonized IL10 UMI | Positive cells |
|---|---:|---:|
| SfLB1 | 0 | 0 |
| SfLB2 | 0 | 0 |
| Tx17 | 5 | 4 |
| Tx5 | 0 | 0 |

Tx17 strict-harmonized lineage distribution:

- Mesenchymal_like = 2 UMI
- Hepatocyte_like = 1 UMI
- Unresolved = 2 UMI

The strongest strict-harmonized cell remained:

- Perivascular_mesenchymal_HSC_like
- b1bd8e2 = 2 UMI

## Figure audit

v1.4.2B.3.2 Fig03 contained a visualization-only cell-label collision because shortened barcode suffixes were used as factor identities.

This was corrected in v1.4.2B.3.3.1.

Corrected Fig03 audit:

- unique audited RDS cells = 19
- unique factor identities = 19
- unique display labels = 19

Therefore each audited RDS cell is represented exactly once.

Final figures include:

- Fig01_IL10_sample_accepted_vs_filteredout_UMI
- Fig02_IL10_cell_signal_strength
- Fig03_IL10_cell_probe_heatmap
- Fig04_IL10_accepted_UMI_by_lineage
- Fig05_IL10_strict_harmonized_shared_probe_by_sample
- Fig06_IL10_strict_harmonized_by_lineage

## Interpretation

The transgene-negative SfLB controls show no Cell Ranger-accepted IL10 signal.

The transgene-positive Tx grafts show Cell Ranger-accepted IL10 signal, and the strongest Tx17 signal occurs in a mesenchymal cell. The strongest cell is supported by two independent accepted IL10 probes, including the shared accepted spliced probe b1bd8e2.

The strict harmonized analysis using only b1bd8e2 retains IL10 signal in Tx17 while both SfLB controls remain zero.

These observations are strongly concordant with the experimental design in which the IL10 transgene was introduced into the mesenchymal component of the iPSC-LO.

However, the IL10 probes target IL10 sequence shared with endogenous IL10. Therefore this analysis does not by itself prove that an individual detected IL10 molecule originated from the transgene.

## Output directories

Validated integrated output:

/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.2/analysis/v1.4.2B.3.2_IL10_probe_integrated_audit

Post-audit corrected / harmonized output:

/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.3.1/analysis/v1.4.2B.3.3.1_IL10_probe_integrated_audit

## Explicit exclusions from checkpoint

The following failed, superseded, intermediate, or backup files are NOT included:

- Human_MASH_v1.4.2B.3.1_IL10_ProbeIntegratedAudit.sh
- Human_MASH_v1.4.2B.3.2_IL10_ProbeIntegratedAudit.sh.bak_before_metadata_rbind_fix_20260908
- Human_MASH_v1.4.2B.3.3_IL10_PostAuditRepair_HarmonizedSensitivity.sh
- Human_MASH_v1.4.2B.3_SfLB_IL10_ProbeNegativeControl.sh

No unrelated untracked Mouse/Human analysis scripts or analysis-output directories are included in this checkpoint.
