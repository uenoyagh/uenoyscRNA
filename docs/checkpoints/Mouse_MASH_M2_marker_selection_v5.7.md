# Mouse MASH M2 / anti-inflammatory macrophage marker checkpoint

Date: 2026-08-27

## Final primary marker panel

The current primary markers selected for M2 / anti-inflammatory macrophage
evaluation are:

- CD163
- CLEC4F
- SLC40A1

## Biological interpretation

### CD163
Anti-inflammatory / scavenging macrophage phenotype.

### CLEC4F
Resident / homeostatic Kupffer-cell identity.

### SLC40A1
M2-like iron-export / homeostatic metabolic phenotype.

## Experimental context

Mouse MASH model:
- NOD/scid background
- AAV8-OSM administered 1 week before transplantation
- tissue collected 4 weeks after transplantation

Because all Sham and Tx animals receive the same AAV8-OSM pretreatment,
Sham-vs-Tx comparisons are interpreted as transplantation-associated changes
within an AAV8-OSM-conditioned NOD/scid liver environment.

## Analysis checkpoints supporting marker selection

- Mouse_MASH_CD163_positive_Mphi_deNovo_marker_discovery_v5.5.0.R
- Mouse_MASH_CD163like_Sham_to_Tx_decrease_screen_v5.6.0.R
- Mouse_MASH_Kupffer_contamination_reaudit_v5.6.4.R
- Mouse_MASH_CD163_UMAPlike_Tx_increase_deNovo_screen_v5.7.1.R
- Mouse_MASH_CD163_candidate_UMAP_panel_v5.7.3.R

## Current interpretation

CD163 / CLEC4F / SLC40A1 were selected after considering:

- Cd163-like UMAP localization
- Sham-vs-Tx expression changes
- RNA-positive cell fraction
- biological replicate consistency
- contamination robustness
- macrophage / Kupffer-cell biological function

SEMA6D remains biologically interesting, particularly as a macrophage
anti-inflammatory polarization molecule, but is not included in the current
primary three-marker panel.
