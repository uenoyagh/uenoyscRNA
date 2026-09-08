# Human MASH 12-sample Context Freeze checkpoint v1.4.3A.1

Date: 2026-09-08

## Repository state at checkpoint creation

- Repository: /Users/uenoya/Projects/uenoyscRNA
- Branch: main
- Parent HEAD: 12c0fdb746a7a8a7d2a4535b12ec7bd43da2e3c9
- Parent HEAD short: 12c0fdb

## Analysis versions

- 12-sample context review: v1.4.3
- Successful context-review implementation: v1.4.3.1
- Corrected frozen-context version: v1.4.3A.1

## Successful scripts included

1. scripts/Human_MASH/Human_MASH_v1.4.3.1_12Sample_IntegratedContextReview.R
   - SHA256: 989e14ff1a1e698feed9e615520339a61eff4add050c341d1cadb4deef55ca40
   - Purpose:
     - review LT1-LT8, SfLB1/SfLB2, and Tx17/Tx5 in the pre-existing 12-sample RPCA/UMAP coordinate system
     - preserve existing umapRPCA coordinates
     - reconcile LT v1.3.3 frozen annotations
     - reconcile iPSC v1.4.2A/B annotations
     - overlay validated IL10 probe-audit results
     - generate sample, source-group, lineage, split-UMAP, and IL10 context figures

2. scripts/Human_MASH/Human_MASH_v1.4.3A.1_12Sample_ContextFreeze_Correction.R
   - SHA256: 4ec093d3a3dc4c181897b271dae0eea6e4442c8d7f76b8089fcd92fe5db3d87c
   - Purpose:
     - create corrected canonical 12-sample frozen RDS
     - extract direct RNA raw-count IL10 for all 124,763 cells
     - retain validated SfLB/Tx probe-level IL10 audit
     - set LT1-LT8 probe-level IL10 fields to NA rather than false zero
     - verify assays, reductions, identities, cell order, and umapRPCA coordinates remain unchanged after save/reload

## Canonical frozen RDS

/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3A.1/objects/Human_MASH_12sample_Context_v1.4.3A.1_frozen.rds

Validated integrity:

- Size: 14,926,923,385 bytes
- MD5: f87e1d051f57da246474cb1995c80b71

## Canonical object dimensions

- Total cells: 124,763
- Total genes: 9,137
- LT1-LT8: 81,348 cells
- iPSC-derived samples: 43,415 cells
- Existing UMAP reduction preserved: umapRPCA

Sample cell counts:

| Sample | Cells |
|---|---:|
| LT1 | 10,847 |
| LT2 | 12,227 |
| LT3 | 7,906 |
| LT4 | 9,899 |
| LT5 | 12,172 |
| LT6 | 10,663 |
| LT7 | 7,734 |
| LT8 | 9,900 |
| SfLB1 | 19,081 |
| SfLB2 | 9,805 |
| Tx17 | 5,958 |
| Tx5 | 8,571 |

## Direct RNA IL10 validation

Gene-level IL10 was extracted directly from RNA LayerData(counts) for all 12 samples.

| Sample | IL10 UMI | IL10-positive cells |
|---|---:|---:|
| LT1 | 83 | 63 |
| LT2 | 180 | 145 |
| LT3 | 236 | 133 |
| LT4 | 62 | 47 |
| LT5 | 126 | 101 |
| LT6 | 199 | 118 |
| LT7 | 46 | 42 |
| LT8 | 113 | 76 |
| SfLB1 | 0 | 0 |
| SfLB2 | 0 | 0 |
| Tx17 | 18 | 13 |
| Tx5 | 4 | 4 |

Total:

- IL10 UMI = 1,067
- IL10-positive cells = 742

## Probe-level IL10 semantics

### LT1-LT8

No LT probe-level audit was performed in this workflow.

Therefore the following fields are deliberately NA:

- accepted_IL10_UMI_v143A1
- filteredout_IL10_UMI_v143A1
- strict_b1bd8e2_IL10_UMI_v143A1
- IL10_accepted_positive_v143A1
- IL10_strict_b1bd8e2_positive_v143A1
- IL10_gene_vs_accepted_concordant_v143A1

probe_IL10_audit_available_v143A1 = FALSE.

### SfLB / Tx

Validated probe-level audit from v1.4.2B.3.3.1 is retained.

- accepted IL10 = 22 UMI / 17 cells
- strict shared accepted spliced probe b1bd8e2 = 5 UMI / 4 cells
- direct RNA gene-level IL10 matches validated iPSC gene-level IL10 cellwise

## v1.4.3 context review

The existing 12-sample RPCA/UMAP coordinates were used without recomputation.

No:

- re-integration
- PCA recomputation
- neighbor reconstruction
- reclustering
- UMAP recomputation

was performed.

Visual audit showed:

- LT1-LT8 share the main native-liver manifold
- SfLB1/SfLB2 show replicate-consistent occupancy
- Tx17/Tx5 show replicate-consistent occupancy
- SfLB and Tx exhibit reproducible condition-associated redistribution
- harmonized lineage structure remains biologically coherent

UMAP geometry is retained as descriptive context only and is not interpreted as lineage trajectory or quantitative maturation distance.

## v1.4.3A.1 correction

The superseded v1.4.3A freeze incorrectly initialized LT1-LT8 IL10 metadata to zero because only the iPSC probe-audit table was overlaid.

v1.4.3A.1 corrects this by:

1. extracting gene-level IL10 directly from RNA counts for all 12 samples;
2. validating the previously established LT1-LT8 IL10 counts;
3. requiring direct RNA IL10 to match validated SfLB/Tx gene-level IL10 cellwise;
4. marking LT probe-level fields as NA rather than zero;
5. preserving source object structure and umapRPCA coordinates exactly.

## Freeze validation

All checks passed:

- DIRECT RNA IL10 VALIDATED FOR ALL 12 SAMPLES
- LT PROBE-LEVEL FIELDS ARE NA, NOT FALSE ZERO
- SOURCE STRUCTURE AND UMAP UNCHANGED AFTER SAVE/RELOAD
- v1.4.3A.1 CORRECTED 12-SAMPLE CONTEXT FREEZE COMPLETE

## Canonical downstream use

Use:

Human_MASH_12sample_Context_v1.4.3A.1_frozen.rds

as the canonical 12-sample input for downstream Human MASH analyses beginning with IL10 receptor competence / response analysis.

Do not use the superseded v1.4.3A frozen RDS for IL10 metadata interpretation.

## Explicit exclusions from this checkpoint

The following are not included:

- Human_MASH_v1.4.3_12Sample_IntegratedContextReview.R
  - failed syntax version
- Human_MASH_v1.4.3A_12Sample_ContextFreeze.R
  - superseded because LT IL10 metadata were incorrectly initialized to zero
- unrelated Human MASH analysis scripts
- unrelated Mouse MASH analysis scripts
- backup files
- analysis output directories
- the large RDS itself, which remains on the external SSD
