#!/bin/bash
set -euo pipefail

PROJECT="/Users/uenoya/Projects/uenoyscRNA"
PATCH_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$PROJECT/R"
mkdir -p "$PROJECT/analysis/RDS3_annotation_visualization_v4.0"

cp "$PATCH_DIR/R/04_plot_export_v4.0.R" \
   "$PROJECT/R/04_plot_export_v4.0.R"

cp "$PATCH_DIR/analysis/RDS3_annotation_visualization_v4.0/99_run_all.R" \
   "$PROJECT/analysis/RDS3_annotation_visualization_v4.0/99_run_all.R"

chmod +x \
  "$PROJECT/analysis/RDS3_annotation_visualization_v4.0/99_run_all.R"

echo "Installed v4.0 visualization patch."
echo
echo "Run:"
echo "cd \"$PROJECT\""
echo "Rscript analysis/RDS3_annotation_visualization_v4.0/99_run_all.R"
