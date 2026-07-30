#!/bin/bash
set -euo pipefail

PROJECT="/Users/uenoya/Projects/uenoyscRNA"
PATCH_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$PROJECT/R"
mkdir -p "$PROJECT/analysis/RDS3_annotation_validation_v3.3"

cp "$PATCH_DIR/R/02_markers_v3.3.R" \
   "$PROJECT/R/02_markers_v3.3.R"

cp "$PATCH_DIR/R/03_hierarchical_voting_v3.3.R" \
   "$PROJECT/R/03_hierarchical_voting_v3.3.R"

cp "$PATCH_DIR/analysis/RDS3_annotation_validation_v3.3/99_run_all.R" \
   "$PROJECT/analysis/RDS3_annotation_validation_v3.3/99_run_all.R"

chmod +x "$PROJECT/analysis/RDS3_annotation_validation_v3.3/99_run_all.R"

echo "Installed v3.3 patch into:"
echo "  $PROJECT"
echo
echo "Run:"
echo "  cd \"$PROJECT\""
echo "  Rscript analysis/RDS3_annotation_validation_v3.3/99_run_all.R"
