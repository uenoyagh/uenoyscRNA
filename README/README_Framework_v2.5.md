# uenoy scRNAseq Framework v2.5

## Unified signature pipeline

The same signature lists are now used to generate:

- cluster-level functional heatmaps
- signature-score UMAPs
- signature-score violin plots
- representative-gene DotPlots
- score tables
- gene matching reports

Both PDF and PNG are exported.

PNG is written at 300 dpi by default for lightweight PowerPoint use.

## Run

```r
source(
  "/Users/uenoya/Projects/uenoyscRNA/analysis/05_signature_pipeline.R"
)
```

## Output structure

```text
analysis_results/
└── 05_module_score/
    └── signature_pipeline_v2.5/
        ├── Mouse_Liver_AllCell/
        │   └── <RDS name>/
        │       ├── 01_heatmap/
        │       ├── 02_featureplot/
        │       ├── 03_violin/
        │       ├── 04_dotplot/
        │       └── 05_tables/
        └── Mouse_Macrophage/
            └── <RDS name>/
                ├── 01_heatmap/
                ├── 02_featureplot/
                ├── 03_violin/
                ├── 04_dotplot/
                └── 05_tables/
```

All four RDS files are processed by default.

## PNG settings

```r
signature_export_png <- TRUE
signature_png_dpi <- 300
signature_png_background <- "white"
```

PDF output can be disabled independently:

```r
signature_export_pdf <- FALSE
```
