#!/bin/bash
# Find_Flex_Raw_Background_COMPLETE_v1.2.sh
# Search all currently mounted volumes for standard and Flex-specific raw/background files.
# Does not modify anything.

set -u

OUT="/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/CrossSpecies_MASH_Mphi/Flex_raw_background_inventory_v1.2.tsv"
mkdir -p "$(dirname "$OUT")"
printf "kind\tpath\n" > "$OUT"

echo "===== Mounted volumes ====="
ls -1 /Volumes 2>/dev/null || true
echo

find /Volumes \
  \( \
    -type f \( \
      -name 'raw_feature_bc_matrix.h5' -o \
      -name 'raw_feature_bc_matrix_h5.h5' -o \
      -name 'sample_raw_feature_bc_matrix.h5' -o \
      -name 'raw_molecule_info.h5' -o \
      -name 'sample_molecule_info.h5' -o \
      -name 'molecule_info.h5' \
    \) -o \
    -type d \( \
      -name 'raw_feature_bc_matrix' -o \
      -name 'sample_raw_feature_bc_matrix' \
    \) \
  \) \
  -print 2>/dev/null | sort -u | while IFS= read -r p; do
    b="$(basename "$p")"
    case "$b" in
      raw_feature_bc_matrix.h5|raw_feature_bc_matrix_h5.h5) kind="raw_h5" ;;
      sample_raw_feature_bc_matrix.h5) kind="sample_raw_h5" ;;
      raw_feature_bc_matrix) kind="raw_mex_dir" ;;
      sample_raw_feature_bc_matrix) kind="sample_raw_mex_dir" ;;
      raw_molecule_info.h5) kind="raw_molecule_info" ;;
      sample_molecule_info.h5) kind="sample_molecule_info" ;;
      molecule_info.h5) kind="molecule_info" ;;
      *) kind="other" ;;
    esac
    printf "%s\t%s\n" "$kind" "$p" >> "$OUT"
  done

echo "===== Raw/background candidates ====="
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$OUT" | head -400
else
  head -400 "$OUT"
fi

echo
echo "===== Number of candidates ====="
awk 'NR>1{n++} END{print n+0}' "$OUT"
echo
echo "Saved: $OUT"

echo
echo "===== Repository references to raw/Flex source files ====="
cd "/Users/uenoya/Projects/uenoyscRNA" || exit 0
grep -RniE \
'sample_raw_feature_bc_matrix|raw_feature_bc_matrix|raw_molecule_info\.h5|sample_molecule_info\.h5|molecule_info\.h5|/Volumes/RH[0-9]+' \
analysis scripts R 2>/dev/null | head -400 || true
