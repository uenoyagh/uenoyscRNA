# Mouse MASH Hepatocyte + MΦ–Hepatocyte checkpoint v6.6.4.4

Date: 2026-09-02

## Checkpoint

Tag:

`mouse-mash-mphi-hep-mechanistic-closure-v6.6.4.4`

This checkpoint preserves the successful Hepatocyte and MΦ–Hepatocyte analysis series following completion of the MΦ–HSC mechanistic closure checkpoint.

Previous checkpoint:

- commit: `7868c25`
- tag: `mouse-mash-mphi-hsc-mechanistic-closure-v6.3.6`

## Hepatocyte successful versions

- `v6.5.0`
  - `analysis/Mouse_MASH_Hepatocyte_subclustering_v6.5.0.R`
- `v6.5.1`
  - `analysis/Mouse_MASH_Hepatocyte_cleanup_lineage_audit_v6.5.1.R`
- `v6.5.2`
  - `analysis/Mouse_MASH_Hepatocyte_FINAL_state_pseudobulk_v6.5.2.R`
- `v6.5.3.1`
  - `analysis/Mouse_MASH_Hepatocyte_ballooning_focused_v6.5.3.1.R`

## MΦ–Hepatocyte successful versions

- `v6.6.0`
  - `analysis/Mouse_MASH_Mphi5_Hep_interaction_prepare_v6.6.0.R`
- `v6.6.1`
  - `analysis/Mouse_MASH_Mphi5_Hep5_CellChat_samplewise_v6.6.1.R`
- `v6.6.1.1`
  - `analysis/Mouse_MASH_Mphi5_Hep5_CellChat_priority_postprocess_v6.6.1.1.R`
- `v6.6.2.1`
  - `analysis/Mouse_MASH_Mphi5_Hep5_CellChat_refine_v6.6.2.1.R`
- `v6.6.3`
  - `analysis/Mouse_MASH_Mphi5_Hep5_NicheNet_v6.6.3.R`
- `v6.6.3.1`
  - `analysis/Mouse_MASH_Mphi5_Hep5_NicheNet_concordance_v6.6.3.1.R`
- `v6.6.4.4`
  - `analysis/Mouse_MASH_Mphi5_Hep5_mechanistic_closure_correction_v6.6.4.4.R`

## Explicitly excluded / superseded versions

These versions are intentionally not included in this checkpoint:

- Hepatocyte `v6.5.3`
- MΦ–Hepatocyte `v6.6.2`
- MΦ–Hepatocyte `v6.6.4`
- MΦ–Hepatocyte `v6.6.4.1`
- MΦ–Hepatocyte `v6.6.4.2`
- MΦ–Hepatocyte `v6.6.4.3`

The duplicate repository-root file below is also intentionally excluded:

- `Mouse_MASH_Hepatocyte_ballooning_focused_v6.5.3.1.R`

The canonical `v6.5.3.1` file is the copy under `analysis/`.

## Validation

The 11 successful R scripts were syntax-checked with R `parse()` before staging.

Only the 11 successful scripts and this checkpoint markdown are intended for the checkpoint commit.

No bulk staging command such as `git add .` is used.
