# 02_publish_UMAP.R update

Copy the included files to the corresponding locations under:

`/Users/uenoya/Projects/uenoyscRNA/`

Files:

- `analysis/02_publish_UMAP.R`
- `R/umap.R`
- `R/plotting.R`
- `config/local_config.R`

Execution:

```r
source("/Users/uenoya/Projects/uenoyscRNA/analysis/02_publish_UMAP.R")
```

Output root:

`/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Mouse_MASH_RDS/analysis_results/01_UMAP/`

The first RDS has no valid cell-type annotation column. For that RDS, annotation plots are skipped automatically and cluster-based plots are generated instead.
