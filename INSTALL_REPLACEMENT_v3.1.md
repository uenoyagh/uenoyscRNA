# Replacement installation: uenoyscRNA v3.1

## 1. Back up the current framework

In Terminal:

```bash
cd /Users/uenoya/Projects
mv uenoyscRNA uenoyscRNA_backup_before_v3.1
```

## 2. Place the supplied folder

Copy the supplied `uenoyscRNA_v3.1` folder to:

```text
/Users/uenoya/Projects/uenoyscRNA
```

The final path must therefore be:

```text
/Users/uenoya/Projects/uenoyscRNA/uenoyscRNA.Rproj
```

## 3. Preserve local machine settings

Compare the following files with the backup before running analyses:

```text
config/local_config.R
config/project_config.R
```

Restore machine-specific RDS and result paths from the backup when necessary.

## 4. Open and validate

Open `uenoyscRNA.Rproj`, then run:

```r
source("analysis/06_cell_fraction_transition.R")
```

Review every generated:

```text
08_QC/*_metadata_detection.csv
```

before interpreting the plots.
