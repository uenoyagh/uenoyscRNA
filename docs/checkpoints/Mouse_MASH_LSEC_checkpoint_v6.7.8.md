# Mouse MASH LSEC mechanistic closure checkpoint v6.7.8

## Checkpoint

- Tag: `mouse-mash-lsec-mechanistic-closure-v6.7.8`
- Analysis compartment: Mouse liver LSEC / endothelial cells
- Final mechanistic closure version: `v6.7.8`
- Previous major checkpoint: `mouse-mash-mphi-hep-mechanistic-closure-v6.6.4.4`

## Successful analysis series

The following scripts constitute the validated LSEC analysis series:

1. `Mouse_MASH_LSEC_endothelial_parent_audit_v6.7.0.R`
   - extracted and audited the endothelial parent population

2. `Mouse_MASH_LSEC_RPCA_clean_rebuild_resolution_scan_v6.7.1.1.R`
   - rebuilt the endothelial object from raw RNA counts
   - removed inherited reductions/graphs/integrated assay state
   - performed clean sample-aware RPCA integration
   - endothelial parent: 14,402 cells

3. `Mouse_MASH_LSEC_res0.6_expanded_marker_audit_v6.7.2.R`
   - fixed parent audit resolution at res0.6
   - expanded lineage marker audit
   - identified non-endothelial contamination

4. `Mouse_MASH_LSEC_marker_table_correction_v6.7.2.1.R`
   - regenerated significant marker tables using adjusted P < 0.05
   - correction of marker-ranking output only

5. `Mouse_MASH_LSEC_lineage_cleanup_v6.7.3.1.R`
   - removed clear non-endothelial populations
   - clean endothelial parent: 11,203 cells
   - clean LSEC candidate: 9,190 cells
   - non-endothelial cells removed: 3,199
   - rebuilt a raw-count-only clean LSEC object

6. `Mouse_MASH_LSEC_clean_RPCA_resolution_scan_v6.7.4.R`
   - sample-aware RPCA reclustering of 9,190 clean LSEC cells
   - baseline PCA dimensions: 1:20
   - resolution scan: 0.1–0.6

7. `Mouse_MASH_LSEC_sample_bias_QC_audit_v6.7.4.1.R`
   - assessed sample-specific clustering and QC effects
   - compared res0.2 and res0.3

8. `Mouse_MASH_LSEC_sample_specific_contamination_audit_v6.7.4.2.R`
   - quantified LSEC versus hepatocyte transcript fractions
   - evaluated Tx5- and Tx17-enriched states
   - no additional cells removed

9. `Mouse_MASH_LSEC_annotation_freeze_v6.7.5.R`
   - adopted res0.3 as the final LSEC annotation scaffold
   - defined 7 LSEC states

10. `Mouse_MASH_LSEC_pseudobulk_Shams_vs_Tx_v6.7.6.R`
    - sample-level edgeR pseudobulk analysis
    - biological replicate as statistical unit
    - Sham n=2, Tx n=2
    - Primary, Shared-core, and state-specific sensitivity analyses

11. `Mouse_MASH_LSEC_Hepatocyte_specificity_audit_v6.7.6.1.R`
    - classified Tx-response genes by Hepatocyte versus LSEC expression
    - separated Hepatocyte-dominant, ambiguous/shared, and LSEC-supported genes
    - ambient hepatocyte RNA was explicitly considered before mechanistic interpretation

12. `Mouse_MASH_LSEC_pathway_ambient_aware_v6.7.7.R`
    - ambient-aware ranked GSEA
    - Primary_no_QC and Shared-core analyses
    - Hallmark and GO Biological Process pathway concordance

13. `Mouse_MASH_LSEC_mechanistic_closure_v6.7.8.R`
    - integrated disease-axis and transplantation-axis findings
    - generated six-sample program heatmap
    - summarized concordant Hallmark pathways
    - generated final LSEC-supported candidate concordance output

## Final LSEC scaffold

Clean LSEC:

- 9,190 cells
- final clustering scaffold: RPCA res0.3
- PCA baseline: dimensions 1:20

Final annotated states:

1. `Inflammatory_stress_high_LSEC`
2. `Homeostatic_like_LSEC`
3. `Tx5_enriched_LSEC_state`
4. `Wnt_angiocrine_high_LSEC`
5. `Low_quality_ambient_enriched_LSEC`
6. `Tx17_enriched_Cd209b_Ctsj_LSEC`
7. `Cycling_LSEC`

Interpretation flags:

- Tx5-enriched and Tx17-enriched states are retained but treated as sample-specific exploratory states.
- Low-quality/ambient-enriched LSEC is QC-flagged.
- These populations were not simply deleted before sensitivity analyses.

## Main biological interpretation

Disease-axis comparison (CDHFD vs STD; descriptive because n=1 each):

- LSEC identity decreased.
- Wnt/angiocrine program decreased.
- Capillarization increased.
- Inflammatory/adhesion program increased.
- IFN response increased.
- Stress/immediate-early-response program increased.
- MHC-II program increased.

This is consistent with MASH-associated LSEC dysfunction/capillarization.

Transplantation-axis comparison (Tx vs Sham; n=2 each):

- no strong uniform normalization of the entire LSEC compartment was detected
- LSEC identity showed a modest recovery direction
- Wnt/angiocrine activity showed a modest recovery direction
- most prespecified module-score differences between Sham and Tx were small
- sample-to-sample heterogeneity was substantial

Ambient-aware ranked pathway analysis nevertheless detected concordant Tx-associated transcriptional remodeling in Primary and Shared-core LSEC analyses, including:

- TNFA signaling via NFKB: Tx-up
- Hypoxia: Tx-up
- Inflammatory response: Tx-up
- P53 / stress-response programs: Tx-up
- TGF-beta signaling: Tx-up
- Coagulation: Tx-down

Therefore, the transplantation response is not interpreted as simple anti-inflammatory normalization.

The preferred interpretation is:

> LSEC represents a heterogeneous secondary remodeling compartment after transplantation, with partial restoration of selected endothelial identity/angiocrine features but persistent or remodeled stress, inflammatory, hypoxic, and vascular-response programs.

The strongest mouse mechanistic conclusions remain centered on the macrophage-HSC and macrophage-hepatocyte axes. LSEC findings are supportive/exploratory and may become more important when integrated later with human MASH and human graft-to-mouse interaction analyses.

## LSEC-supported exploratory Tx-response genes

Representative Tx-up candidates shared by Primary and Shared-core analyses:

- `Kcnj8`
- `Ltbp4`
- `Flrt1`
- `Matn4`
- `Cebpd`
- `Abcc9`

Representative Tx-down candidates:

- `Car11`
- `Pi16`
- `Rab3b`
- `Gnai1`
- `Nid2`
- `Col5a2`

These gene-level signals are exploratory; conventional strong gene-level FDR significance was not obtained.

## Explicitly excluded / superseded versions

Not included in this checkpoint:

- `v6.7.1`
  - superseded by the clean raw-count rebuild in v6.7.1.1
  - previous object state generated excessive inherited assay/reduction warnings

- `v6.7.3`
  - failed because deprecated `GetAssayData(slot="counts")` was defunct under SeuratObject 5
  - corrected in v6.7.3.1 using `layer="counts"`

## Reproducibility audit

Before staging this checkpoint, all 13 selected successful R scripts were parsed successfully with R.

Result:

`RESULT: ALL 13 SCRIPTS PASSED`

No broad `git add .` is to be used for this checkpoint.
Only the explicitly listed successful scripts and this checkpoint document are staged.
