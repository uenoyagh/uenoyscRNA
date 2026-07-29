# uenoy scRNAseq Framework v2.3

## UMAP output folders

### R8 colour version

`analysis_results/01_UMAP/R8/`

### Pastel colour version

`analysis_results/01_UMAP/Pastel/`

### Transparent monochrome sample UMAPs

`analysis_results/01_UMAP/Monochrome_transparent/07_each_sample/`

The individual sample PDFs use the same global UMAP coordinate limits and a
transparent PDF background so that they can be overlaid later.

## Functional signature heatmap

Run:

```r
source("/Users/uenoya/Projects/uenoyscRNA/analysis/06_functional_signature_heatmap.R")
```

Output:

`analysis_results/05_module_score/functional_signature_heatmap/`

The heatmap uses the standard blue-white-red scale:

- low: `#0033FF`
- midpoint: `#FFFFFF`
- high: `#FF1A1A`
