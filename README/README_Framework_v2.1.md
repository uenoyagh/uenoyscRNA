# uenoy scRNAseq Framework v2.1

## Main changes

- Vector PDF is now the default.
- UMAP point size increased from 0.20 to 0.40.
- A vivid high-contrast palette is generated automatically.
- Cluster labels are printed directly on the UMAP.
- The cluster legend is hidden by default when labels are displayed.
- Every PDF contains a small footer with:
  - source RDS filename
  - creation date and time
  - framework name and version
- Plot titles were simplified.

## Installation

Copy the contents of this update into:

`/Users/uenoya/Projects/uenoyscRNA/`

The following files should be overwritten:

- `config/project_config.R`
- `config/local_config.R`
- `R/umap.R`
- `R/plotting.R`
- `analysis/02_publish_UMAP.R`

## Important before rerunning

Existing PDFs are not overwritten when:

```r
overwrite_existing <- FALSE
```

To regenerate the existing UMAP PDFs, temporarily set:

```r
overwrite_existing <- TRUE
```

in `config/local_config.R`.

Then run:

```r
source("/Users/uenoya/Projects/uenoyscRNA/analysis/02_publish_UMAP.R")
```

After regeneration, return the setting to:

```r
overwrite_existing <- FALSE
```
