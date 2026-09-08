# Human MASH Canonical Annotation Freeze v1.3.3

Date: 2026-09-08

## Repository

- Repository: `/Users/uenoya/Projects/uenoyscRNA`
- Branch: `main`
- Previous checkpoint commit: `32c93c4`
- Previous tag: `human-mash-myeloid-stage-pattern-v2.3`

## Purpose

Freeze the refined Human MASH LT1-LT8 whole-liver annotation after targeted immune and hepatobiliary reannotation.

## Canonical v1.3.3

Canonical object:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3.3/objects/Human_MASH_LT1to8_Canonical_v1.3.3.rds`

- Canonical cells retained: 81,348 / 81,348
- Final broad cell types: 25
- Final refined states: 42
- Original canonical v1.3 object was not modified.

## Immune refinement

Integrated:

- v1.3.2A Immune Refinement
- v1.3.2A.1 Immune Manual Annotation Cleaning

Immune candidate cells merged: 14,833.

Broad immune annotation was manually refined and nonimmune contaminants were reassigned at the whole-liver freeze stage when supported.

## Hepatobiliary refinement

Integrated:

- v1.3.2B Hepatobiliary Refinement
- v1.3.2B-QC1 LSEC-vs-Hepatocyte Audit
- v1.3.2B.1 Hepatobiliary Manual Annotation Cleaning

Hepatobiliary candidate cells merged: 47,425.

Clean hepatobiliary cells before whole-liver merge:

- total clean hepatobiliary: 46,817
- Hepatocyte: 42,613
- Cholangiocyte: 3,442
- Ductular progenitor: 762

The previous broad `Hepatoblast_or_immature_hepatocyte` label was retired because the compartment resolved into adult hepatocyte, transitional, ductular/progenitor, biliary, and contaminating lineages; broad fetal hepatoblast support was not convincing.

Clusters 0 and 7 passed the targeted QC1 LSEC-vs-Hepatocyte audit and were retained on the hepatocyte side rather than reclassified as LSEC.

## Whole-liver conflict handling

Immune/hepatobiliary overlap cells audited: 434.

- excluded-excluded same-lineage consensus: 1 cell
- excluded-excluded conflicting lineage calls: 128 cells
- ambiguous fraction: 0.157% of canonical cells

The 128 conflicting cells were retained as:

- `celltype_final_v133 = Ambiguous_crosslineage`
- `cellstate_final_v133 = Doublet_or_crosslineage_ambiguous`
- `annotation_source_v133 = refinement_conflict_ambiguous`

These cells should be excluded from lineage-specific biological inference unless separately resolved.

No clean-clean overlap conflict was accepted silently.

## Caution flags

Strong sample-bias flag retained for hepatobiliary clusters:

- cluster 9: EPCAM-positive ductular/progenitor-like, LT2-biased
- cluster 12: mature hepatocyte periportal-like, LT4-biased
- cluster 14: injury/transitional hepatocyte, LT1-biased

Total cells carrying the >=70% single-sample state-bias flag: 1,249.

Outside-core strict ductal candidates: 44 cells, audit-only; not automatically rescued.

## Final UMAP policy

Beginning with canonical v1.3.3:

- final PDF uses Cairo when available
- `raster = FALSE`
- whole-liver overview `pt.size = 0.18`
- split-by-sample `pt.size = 0.12`
- high-contrast saturated palette
- vector output for final review/publication figures

Final v1.3.3 UMAPs were successfully generated and visually reviewed.

## Scripts in this checkpoint

1. `scripts/Human_MASH/Human_MASH_v1.3.1_WholeLiver_AnnotationValidation.R`
2. `scripts/Human_MASH/Human_MASH_v1.3.2A_ImmuneRefinement.R`
3. `scripts/Human_MASH/Human_MASH_v1.3.2A.1_ImmuneManualAnnotationCleaning.R`
4. `scripts/Human_MASH/Human_MASH_v1.3.2B_HepatobiliaryRefinement.R`
5. `scripts/Human_MASH/Human_MASH_v1.3.2B_QC1_LSECvsHepatocyteAudit.R`
6. `scripts/Human_MASH/Human_MASH_v1.3.2B.1_HepatobiliaryManualAnnotationCleaning.R`
7. `scripts/Human_MASH/Human_MASH_v1.3.3_WholeLiverAnnotationFreeze.R`

## Excluded from Git checkpoint

Not staged:

- RDS files
- PDFs
- CSV/analysis outputs on external SSD
- Mouse scripts
- failed/superseded Human scripts
- unexecuted downstream Human MASH scripts
- `.bak` files
- unrelated untracked directories

## Status

`WHOLE_LIVER_ANNOTATION_FROZEN_v1.3.3`

## Next step

Proceed to Human MASH disease-progression analysis using canonical v1.3.3, with sample-aware cell-composition analysis followed by lineage-specific pseudobulk DE/pathway analysis and subsequent human-mouse translational concordance analysis.
