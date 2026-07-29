# uenoy scRNAseq Framework v2.2

## Main changes

1. Split UMAP is rebuilt with `ggplot2::facet_wrap()`.
   - Six samples: 3 columns × 2 rows
   - Four conditions: 2 columns × 2 rows
   - Fixed x/y scales and `coord_fixed(ratio = 1)`

2. Separate monochrome UMAP PDF files are created for every sample.
   - One PDF per sample
   - A single blue colour by default
   - Global UMAP coordinate limits are retained across files
   - Cell count is shown below the title

## New output directory

`analysis_results/01_UMAP/07_each_sample_monochrome/`

## Files to overwrite

- `config/project_config.R`
- `config/local_config.R`
- `R/umap.R`
- `R/plotting.R`
- `analysis/02_publish_UMAP.R`

## Rerun

Set in `config/local_config.R`:

```r
overwrite_existing <- TRUE
```

Then run:

```r
source("/Users/uenoya/Projects/uenoyscRNA/analysis/02_publish_UMAP.R")
```

After regeneration:

```r
overwrite_existing <- FALSE
```
