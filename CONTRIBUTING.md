# Contributing to uenoyscRNA

## Branches

- `master`: stable releases
- `develop`: integration branch
- `feature/<name>`: feature development

Create feature branches from `develop`.

## Before committing

```r
devtools::document()
devtools::test()
devtools::check()
```

## Annotation registry rules

- UTF-8 CSV
- One row per unique species + tissue + dataset + layer + cluster
- Confidence: High, Medium, or Low
- Marker genes separated by semicolons
- Review dates formatted as YYYY-MM-DD
- Biological uncertainty must be retained in confidence and evidence
