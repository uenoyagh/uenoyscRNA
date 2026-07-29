# uenoy scRNAseq Framework v2.1.1

## Split UMAP correction

The previous script used Seurat's default `split.by` layout, which placed
all panels in a single horizontal row. As a result, each UMAP panel became
too narrow and appeared vertically elongated.

Version 2.1.1 changes split UMAP output as follows:

- `combine = FALSE` is used to obtain individual plots.
- Individual panels are arranged explicitly with `patchwork::wrap_plots()`.
- Six samples are arranged as 3 columns × 2 rows.
- Four conditions are arranged as 2 columns × 2 rows.
- `coord_fixed(ratio = 1)` is applied to every UMAP panel.
- A single shared legend is collected on the right.

## Files to overwrite

- `R/umap.R`
- `analysis/02_publish_UMAP.R`
- `config/project_config.R`

## Rerun

Set:

```r
overwrite_existing <- TRUE
```

in `config/local_config.R`, then run:

```r
source("/Users/uenoya/Projects/uenoyscRNA/analysis/02_publish_UMAP.R")
```

After regeneration, return:

```r
overwrite_existing <- FALSE
```
