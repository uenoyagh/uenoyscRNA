#!/bin/bash
set -euo pipefail

VERSION="v1.4.2B.3.2"
BIOLOGICAL_VERSION="v1.4.2B.3"

RDS="/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.1/objects/Human_MASH_iPSC_derived_v1.4.2B.1_mesenchymal_IL10_source_audit.rds"

SFLB_BASE="/Volumes/RH241011956/Deliver_RH241011956/CellRanger/okmt_pool1_outs/per_sample_outs"
TX_BASE="/Volumes/RH251217117/Deliver_RH251217117/Result_human/01_CellRanger/Count_matrix"

OUT="/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2B.3.2/analysis/v1.4.2B.3.2_IL10_probe_integrated_audit"
STAGE="$OUT/stage"
FIG="$OUT/figures"
TAB="$OUT/tables"

mkdir -p "$STAGE" "$FIG" "$TAB"

export RDS_PATH="$RDS"
export OUT_DIR="$OUT"
export STAGE_DIR="$STAGE"
export FIG_DIR="$FIG"
export TABLE_DIR="$TAB"

FORCE="${FORCE:-0}"

echo "============================================================"
echo "Human MASH ${VERSION}"
echo "Integrated IL10 probe audit"
echo "Biological analysis version: ${BIOLOGICAL_VERSION}"
echo "============================================================"
echo "RDS: $RDS"
echo "OUT: $OUT"
echo

if [ ! -f "$RDS" ]; then
    echo "ERROR: RDS not found:"
    echo "$RDS"
    exit 1
fi

stage_sflb_sample () {
    S="$1"
    DONE="$STAGE/${S}.done"

    if [ -f "$DONE" ] && [ "$FORCE" != "1" ]; then
        echo "SKIP $S: cached stage exists"
        return
    fi

    DIR="$SFLB_BASE/$S/count"
    H5="$DIR/sample_raw_probe_bc_matrix.h5"
    FILTERED="$DIR/sample_filtered_barcodes.csv"
    PROBESET="$DIR/probe_set.csv"

    if [ ! -f "$H5" ] || [ ! -f "$FILTERED" ] || [ ! -f "$PROBESET" ]; then
        echo "ERROR: required SfLB files missing for $S"
        exit 1
    fi

    TMPROOT=$(mktemp -d)
    trap 'rm -rf "$TMPROOT"' EXIT

    FEATURETXT="$TMPROOT/features_id.txt"
    BARCODETXT="$TMPROOT/barcodes.txt"
    DATA_BIN="$TMPROOT/data.bin"
    INDEX_BIN="$TMPROOT/indices.bin"
    INDPTR_BIN="$TMPROOT/indptr.bin"

    POSITIVE_TSV="$STAGE/${S}_IL10_probe_positive_raw_barcodes.tsv"
    META_TSV="$STAGE/${S}_IL10_probe_metadata.tsv"

    echo
    echo "============================================================"
    echo "STAGE SfLB: $S"
    echo "============================================================"

    h5dump -d /matrix/features/id -w 1 "$H5" > "$FEATURETXT"
    h5dump -d /matrix/barcodes -w 1 "$H5" > "$BARCODETXT"
    h5dump -d /matrix/data -b LE -o "$DATA_BIN" "$H5" >/dev/null
    h5dump -d /matrix/indices -b LE -o "$INDEX_BIN" "$H5" >/dev/null
    h5dump -d /matrix/indptr -b LE -o "$INDPTR_BIN" "$H5" >/dev/null

    python3 - \
        "$S" \
        "$FEATURETXT" \
        "$BARCODETXT" \
        "$DATA_BIN" \
        "$INDEX_BIN" \
        "$INDPTR_BIN" \
        "$H5" \
        "$FILTERED" \
        "$PROBESET" \
        "$POSITIVE_TSV" \
        "$META_TSV" <<'PY'
import sys
import os
import re
import csv
import subprocess
from array import array

(
    sample,
    feature_txt,
    barcode_txt,
    data_bin,
    index_bin,
    indptr_bin,
    h5_path,
    filtered_csv,
    probe_set_csv,
    positive_tsv,
    meta_tsv
) = sys.argv[1:]

TARGETS = [
    "ENSG00000136634|IL10|6b6892b",
    "ENSG00000136634|IL10|b1bd8e2",
    "ENSG00000136634|IL10|e99f652",
]

SHORT = {
    "ENSG00000136634|IL10|6b6892b": "6b6892b",
    "ENSG00000136634|IL10|b1bd8e2": "b1bd8e2",
    "ENSG00000136634|IL10|e99f652": "e99f652",
}

def clean_h5_string(x):
    return x.split("\\000", 1)[0]

def extract_core(x):
    hits = re.findall(r"[ACGT]{16,40}", x)
    if not hits:
        return None
    return max(hits, key=len)

def read_binary(path, code, expected_itemsize):
    a = array(code)
    size = os.path.getsize(path)
    if size % expected_itemsize != 0:
        raise RuntimeError(f"Unexpected binary size for {path}: {size}")
    with open(path, "rb") as fh:
        a.fromfile(fh, size // expected_itemsize)
    if sys.byteorder != "little":
        a.byteswap()
    if a.itemsize != expected_itemsize:
        raise RuntimeError(f"Unexpected item size for {path}: {a.itemsize}")
    return a

pat = re.compile(r'\((\d+)\):\s*"([^"]*)"')

feature_index = {}
with open(feature_txt, "r", errors="ignore") as fh:
    for line in fh:
        m = pat.search(line)
        if not m:
            continue
        idx = int(m.group(1))
        value = clean_h5_string(m.group(2))
        if value in TARGETS:
            feature_index[value] = idx

missing = [x for x in TARGETS if x not in feature_index]
if missing:
    raise RuntimeError("Missing IL10 probe IDs in H5: " + ", ".join(missing))

def read_h5_enum_at(h5_path, dataset, idx):
    cmd = [
        "h5dump",
        "-d", dataset,
        "-s", str(idx),
        "-c", "1",
        h5_path
    ]
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    m = re.search(r'\(\d+\):\s*(TRUE|FALSE)', p.stdout)
    if not m:
        raise RuntimeError(
            f"Could not parse enum value for {dataset}[{idx}]\n{p.stdout}"
        )
    return m.group(1) == "TRUE"

probe_meta = {}
with open(probe_set_csv, "r", errors="ignore") as fh:
    useful = [line for line in fh if not line.startswith("#")]
reader = csv.DictReader(useful)
for row in reader:
    pid = row.get("probe_id", "")
    if pid in TARGETS:
        probe_meta[pid] = row

with open(meta_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "stage_sample", "rds_sample", "source_delivery", "panel_version",
        "feature_index", "probe_id", "short_probe", "region",
        "included", "filtered_probes", "probe_seq"
    ])
    rds_sample = "SfLB1" if sample == "SfLB_1" else "SfLB2"
    for pid in TARGETS:
        idx = feature_index[pid]
        row = probe_meta.get(pid, {})
        w.writerow([
            sample, rds_sample, "RH241011956",
            "Chromium Human Transcriptome Probe Set v1.0.1",
            idx, pid, SHORT[pid], row.get("region", ""),
            row.get("included", ""),
            read_h5_enum_at(h5_path, "/matrix/features/filtered_probes", idx),
            row.get("probe_seq", "")
        ])

barcodes = {}
with open(barcode_txt, "r", errors="ignore") as fh:
    for line in fh:
        m = pat.search(line)
        if not m:
            continue
        idx = int(m.group(1))
        value = clean_h5_string(m.group(2))
        core = extract_core(value)
        if core is not None:
            barcodes[idx] = (value, core)

filtered_cores = set()
with open(filtered_csv, "r", errors="ignore") as fh:
    for line in fh:
        core = extract_core(line)
        if core is not None:
            filtered_cores.add(core)

data = read_binary(data_bin, "i", 4)
indices = read_binary(index_bin, "q", 8)
indptr = read_binary(indptr_bin, "q", 8)

if len(data) != len(indices):
    raise RuntimeError("data / indices length mismatch")

n_barcodes = len(indptr) - 1
row_to_probe = {feature_index[pid]: SHORT[pid] for pid in TARGETS}
target_rows = set(row_to_probe.keys())
positive = {}

for col in range(n_barcodes):
    start = indptr[col]
    end = indptr[col + 1]
    for k in range(start, end):
        row = indices[k]
        if row not in target_rows:
            continue
        value = int(data[k])
        if value <= 0:
            continue
        probe = row_to_probe[row]
        if col not in positive:
            positive[col] = {"6b6892b": 0, "b1bd8e2": 0, "e99f652": 0}
        positive[col][probe] += value

with open(positive_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "stage_sample", "barcode_index", "barcode", "barcode_core",
        "is_cellranger_filtered",
        "IL10_6b6892b_UMI", "IL10_b1bd8e2_UMI", "IL10_e99f652_UMI",
        "IL10_all3_UMI", "n_positive_IL10_probes"
    ])
    for col in sorted(positive.keys()):
        if col not in barcodes:
            raise RuntimeError(f"Barcode index missing: {col}")
        barcode, core = barcodes[col]
        c = positive[col]
        all3 = c["6b6892b"] + c["b1bd8e2"] + c["e99f652"]
        n_probes = sum(x > 0 for x in [c["6b6892b"], c["b1bd8e2"], c["e99f652"]])
        w.writerow([
            sample, col, barcode, core, core in filtered_cores,
            c["6b6892b"], c["b1bd8e2"], c["e99f652"],
            all3, n_probes
        ])

print(f"{sample}: raw H5 barcodes={n_barcodes:,}")
print(f"{sample}: non-zero probe×barcode entries={len(data):,}")
print(f"{sample}: Cell Ranger filtered barcodes={len(filtered_cores):,}")
print(f"{sample}: IL10-positive raw barcodes={len(positive)}")
print(f"{sample}: IL10 all-3-probe UMI in raw barcodes={sum(sum(x.values()) for x in positive.values())}")
print(f"{sample}: metadata={meta_tsv}")
print(f"{sample}: positives={positive_tsv}")
PY

    rm -rf "$TMPROOT"
    trap - EXIT
    touch "$DONE"
}

stage_tx_sample () {
    S="$1"
    DONE="$STAGE/${S}.done"

    if [ -f "$DONE" ] && [ "$FORCE" != "1" ]; then
        echo "SKIP $S: cached stage exists"
        return
    fi

    H5="$TX_BASE/$S/molecule_info.h5"
    if [ ! -f "$H5" ]; then
        echo "ERROR: Tx molecule_info.h5 missing for $S"
        exit 1
    fi

    TMPROOT=$(mktemp -d)
    trap 'rm -rf "$TMPROOT"' EXIT

    PROBEIDTXT="$TMPROOT/probe_id.txt"
    REGIONTXT="$TMPROOT/region.txt"
    PROBETYPETXT="$TMPROOT/probe_type.txt"
    BARCODETXT="$TMPROOT/barcodes.txt"

    PROBEIDX_BIN="$TMPROOT/probe_idx.bin"
    BARCODEIDX_BIN="$TMPROOT/barcode_idx.bin"
    UMI_BIN="$TMPROOT/umi.bin"
    COUNT_BIN="$TMPROOT/count.bin"

    MOLECULE_TSV="$STAGE/${S}_IL10_molecules.tsv"
    META_TSV="$STAGE/${S}_IL10_probe_metadata.tsv"

    echo
    echo "============================================================"
    echo "STAGE Tx: $S"
    echo "============================================================"

    h5dump -d /probes/probe_id -w 1 "$H5" > "$PROBEIDTXT"
    h5dump -d /probes/region -w 1 "$H5" > "$REGIONTXT"
    h5dump -d /probes/probe_type -w 1 "$H5" > "$PROBETYPETXT"
    h5dump -d /barcodes -w 1 "$H5" > "$BARCODETXT"

    h5dump -d /probe_idx -b LE -o "$PROBEIDX_BIN" "$H5" >/dev/null
    h5dump -d /barcode_idx -b LE -o "$BARCODEIDX_BIN" "$H5" >/dev/null
    h5dump -d /umi -b LE -o "$UMI_BIN" "$H5" >/dev/null
    h5dump -d /count -b LE -o "$COUNT_BIN" "$H5" >/dev/null

    python3 - \
        "$S" \
        "$PROBEIDTXT" \
        "$REGIONTXT" \
        "$PROBETYPETXT" \
        "$BARCODETXT" \
        "$PROBEIDX_BIN" \
        "$BARCODEIDX_BIN" \
        "$UMI_BIN" \
        "$COUNT_BIN" \
        "$H5" \
        "$MOLECULE_TSV" \
        "$META_TSV" <<'PY'
import sys
import os
import re
import csv
import subprocess
from array import array

(
    sample,
    probeid_txt,
    region_txt,
    probetype_txt,
    barcode_txt,
    probeidx_bin,
    barcodeidx_bin,
    umi_bin,
    count_bin,
    h5_path,
    molecule_tsv,
    meta_tsv
) = sys.argv[1:]

TARGETS = [
    "ENSG00000136634|IL10|6b6892b",
    "ENSG00000136634|IL10|b1bd8e2",
    "ENSG00000136634|IL10|e99f652",
]

SHORT = {
    "ENSG00000136634|IL10|6b6892b": "6b6892b",
    "ENSG00000136634|IL10|b1bd8e2": "b1bd8e2",
    "ENSG00000136634|IL10|e99f652": "e99f652",
}

def clean_h5_string(x):
    return x.split("\\000", 1)[0]

def extract_core(x):
    hits = re.findall(r"[ACGT]{16,40}", x)
    if not hits:
        return None
    return max(hits, key=len)

def read_binary(path, code, expected_itemsize):
    a = array(code)
    size = os.path.getsize(path)
    if size % expected_itemsize != 0:
        raise RuntimeError(f"Unexpected binary size for {path}: {size}")
    with open(path, "rb") as fh:
        a.fromfile(fh, size // expected_itemsize)
    if sys.byteorder != "little":
        a.byteswap()
    if a.itemsize != expected_itemsize:
        raise RuntimeError(f"Unexpected item size for {path}: {a.itemsize}")
    return a

pat = re.compile(r'\((\d+)\):\s*"([^"]*)"')

def parse_indexed_strings(path):
    out = {}
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            m = pat.search(line)
            if not m:
                continue
            out[int(m.group(1))] = clean_h5_string(m.group(2))
    return out

probe_ids = parse_indexed_strings(probeid_txt)
regions = parse_indexed_strings(region_txt)
probe_types = parse_indexed_strings(probetype_txt)

target_index = {}
for idx, pid in probe_ids.items():
    if pid in TARGETS:
        target_index[pid] = idx

missing = [pid for pid in TARGETS if pid not in target_index]
if missing:
    raise RuntimeError("Missing IL10 probes in Tx molecule_info.h5: " + ", ".join(missing))

def read_h5_enum_at(h5_path, dataset, idx):
    cmd = [
        "h5dump",
        "-d", dataset,
        "-s", str(idx),
        "-c", "1",
        h5_path
    ]
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    m = re.search(r'\(\d+\):\s*(TRUE|FALSE)', p.stdout)
    if not m:
        raise RuntimeError(
            f"Could not parse enum value for {dataset}[{idx}]\n{p.stdout}"
        )
    return m.group(1) == "TRUE"

with open(meta_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "stage_sample", "rds_sample", "source_delivery", "panel_version",
        "feature_index", "probe_id", "short_probe", "region",
        "probe_type", "included", "filtered_probes", "probe_seq"
    ])
    for pid in TARGETS:
        idx = target_index[pid]
        w.writerow([
            sample, sample, "RH251217117",
            "Chromium Human Transcriptome Probe Set v2.0.0",
            idx, pid, SHORT[pid], regions.get(idx, ""),
            probe_types.get(idx, ""),
            read_h5_enum_at(h5_path, "/probes/included_probes", idx),
            read_h5_enum_at(h5_path, "/probes/filtered_probes", idx), ""
        ])

probe_idx = read_binary(probeidx_bin, "i", 4)
target_rows = []
target_probe = {}
index_to_short = {target_index[pid]: SHORT[pid] for pid in TARGETS}

for row_idx, p in enumerate(probe_idx):
    if p in index_to_short:
        target_rows.append(row_idx)
        target_probe[row_idx] = index_to_short[p]

del probe_idx

barcode_idx = read_binary(barcodeidx_bin, "Q", 8)
umi = read_binary(umi_bin, "I", 4)
count = read_binary(count_bin, "I", 4)

target_bcidx = {i: barcode_idx[i] for i in target_rows}
needed_bcidx = set(target_bcidx.values())

del barcode_idx

barcode_map = {}
with open(barcode_txt, "r", errors="ignore") as fh:
    for line in fh:
        m = pat.search(line)
        if not m:
            continue
        idx = int(m.group(1))
        if idx not in needed_bcidx:
            continue
        value = clean_h5_string(m.group(2))
        core = extract_core(value)
        barcode_map[idx] = (value, core)

missing_bc = needed_bcidx - set(barcode_map.keys())
if missing_bc:
    raise RuntimeError(f"Unresolved barcode indices: {sorted(missing_bc)[:20]}")

with open(molecule_tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "stage_sample", "molecule_row", "probe", "barcode_idx",
        "barcode", "barcode_core", "umi", "read_support"
    ])
    for row in target_rows:
        bi = target_bcidx[row]
        barcode, core = barcode_map[bi]
        w.writerow([
            sample, row, target_probe[row], bi,
            barcode, core, umi[row], count[row]
        ])

probe_counts = {"6b6892b": 0, "b1bd8e2": 0, "e99f652": 0}
for row in target_rows:
    probe_counts[target_probe[row]] += 1

print(f"{sample}: total molecule rows={len(umi):,}")
print(f"{sample}: IL10 molecule rows={len(target_rows)}")
for p in ["6b6892b", "b1bd8e2", "e99f652"]:
    print(f"{sample}: {p} molecule rows={probe_counts[p]}")
print(f"{sample}: metadata={meta_tsv}")
print(f"{sample}: molecules={molecule_tsv}")
PY

    rm -rf "$TMPROOT"
    trap - EXIT
    touch "$DONE"
}

if [ -d "$SFLB_BASE" ]; then
    stage_sflb_sample "SfLB_1"
    stage_sflb_sample "SfLB_2"
else
    echo "SfLB HDD RH241011956 not mounted: SfLB stage not run."
fi

if [ -d "$TX_BASE" ]; then
    stage_tx_sample "Tx17"
    stage_tx_sample "Tx5"
else
    echo "Tx HDD RH251217117 not mounted: Tx stage not run."
fi

REQUIRED_STAGE=(
    "$STAGE/SfLB_1.done"
    "$STAGE/SfLB_2.done"
    "$STAGE/Tx17.done"
    "$STAGE/Tx5.done"
)

READY=1
for F in "${REQUIRED_STAGE[@]}"; do
    if [ ! -f "$F" ]; then
        READY=0
    fi
done

if [ "$READY" != "1" ]; then
    echo
    echo "============================================================"
    echo "STAGING INCOMPLETE"
    echo "============================================================"
    echo "Run this same script again after mounting the other delivery HDD."
    echo "Cached stage files on the SSD will be reused."
    echo
    echo "Current stage markers:"
    ls -1 "$STAGE"/*.done 2>/dev/null || true
    exit 0
fi

echo
echo "============================================================"
echo "ALL STAGES READY: FINAL INTEGRATION"
echo "============================================================"

Rscript --vanilla - <<'RS'
suppressPackageStartupMessages({
    library(Seurat)
})

if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("R package 'ggplot2' is required for figure generation.")
}

library(ggplot2)

RDS <- Sys.getenv("RDS_PATH")
OUT <- Sys.getenv("OUT_DIR")
STAGE <- Sys.getenv("STAGE_DIR")
FIG <- Sys.getenv("FIG_DIR")
TAB <- Sys.getenv("TABLE_DIR")

dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(RDS)
md <- obj@meta.data

extract_barcode <- function(x) {
    vapply(
        as.character(x),
        function(s) {
            if (is.na(s) || !nzchar(s)) return(NA_character_)
            m <- gregexpr("[ACGT]{16,40}", s, perl=TRUE)[[1]]
            if (length(m) == 1 && m[1] == -1) return(NA_character_)
            ml <- attr(m, "match.length")
            hits <- substring(s, first=m, last=m + ml - 1)
            hits[which.max(nchar(hits))]
        },
        character(1)
    )
}

to_bool <- function(x) {
    tolower(as.character(x)) %in% c("true", "1", "t", "yes")
}

PROBES <- c("6b6892b", "b1bd8e2", "e99f652")
SAMPLES <- c("SfLB1", "SfLB2", "Tx17", "Tx5")

stage_sample <- c(
    SfLB1="SfLB_1",
    SfLB2="SfLB_2",
    Tx17="Tx17",
    Tx5="Tx5"
)

condition <- c(
    SfLB1="SfLB",
    SfLB2="SfLB",
    Tx17="Tx",
    Tx5="Tx"
)

transgene_status <- c(
    SfLB1="transgene-negative",
    SfLB2="transgene-negative",
    Tx17="IL10-transgene-positive graft",
    Tx5="IL10-transgene-positive graft"
)

source_delivery <- c(
    SfLB1="RH241011956",
    SfLB2="RH241011956",
    Tx17="RH251217117",
    Tx5="RH251217117"
)

md$barcode_core_v142B31 <- extract_barcode(md$original_cell_barcode)

META_COLS <- c(
    "stage_sample",
    "rds_sample",
    "source_delivery",
    "panel_version",
    "feature_index",
    "probe_id",
    "short_probe",
    "region",
    "probe_type",
    "included",
    "filtered_probes",
    "probe_seq"
)

meta_list <- list()

for (S in SAMPLES) {

    ss <- stage_sample[[S]]

    f <- file.path(
        STAGE,
        paste0(
            ss,
            "_IL10_probe_metadata.tsv"
        )
    )

    if (!file.exists(f)) {
        stop(
            "Missing probe metadata: ",
            f
        )
    }

    x <- read.delim(
        f,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    x$rds_sample <- S

    missing_cols <- setdiff(
        META_COLS,
        colnames(x)
    )

    if (length(missing_cols) > 0) {

        for (nm in missing_cols) {
            x[[nm]] <- NA_character_
        }
    }

    x <- x[
        ,
        META_COLS,
        drop=FALSE
    ]

    meta_list[[S]] <- x
}

probe_meta <- do.call(
    rbind,
    unname(meta_list)
)

rownames(probe_meta) <- NULL
probe_meta$filtered_probes <- to_bool(probe_meta$filtered_probes)
probe_meta$included <- to_bool(probe_meta$included)

write.table(
    probe_meta,
    file=file.path(TAB, "IL10_probe_filter_status_by_sample.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

cell_list <- list()

for (S in SAMPLES) {
    mm <- md[as.character(md$sample_v141) == S, , drop=FALSE]

    if (nrow(mm) == 0) stop("No RDS cells for sample ", S)
    if (anyDuplicated(mm$barcode_core_v142B31)) {
        stop("Duplicated barcode cores in RDS sample ", S)
    }

    cell <- data.frame(
        sample=S,
        condition=condition[[S]],
        transgene_status=transgene_status[[S]],
        source_delivery=source_delivery[[S]],
        rds_cell=rownames(mm),
        barcode_core=mm$barcode_core_v142B31,
        gene_IL10_UMI_RDS=as.numeric(mm$IL10_counts_raw_v142B),
        lineage=as.character(mm$broad_lineage_manual_v142A),
        state=as.character(mm$manual_state_v142A),
        audit_class=as.character(mm$mesenchymal_audit_class_v142B),
        IL10_6b6892b_UMI=0,
        IL10_b1bd8e2_UMI=0,
        IL10_e99f652_UMI=0,
        stringsAsFactors=FALSE
    )

    ss <- stage_sample[[S]]

    if (grepl("^SfLB", S)) {
        f <- file.path(STAGE, paste0(ss, "_IL10_probe_positive_raw_barcodes.tsv"))
        d <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE)

        if (nrow(d) > 0) {
            j <- match(d$barcode_core, cell$barcode_core)
            ok <- !is.na(j)
            cell$IL10_6b6892b_UMI[j[ok]] <- d$IL10_6b6892b_UMI[ok]
            cell$IL10_b1bd8e2_UMI[j[ok]] <- d$IL10_b1bd8e2_UMI[ok]
            cell$IL10_e99f652_UMI[j[ok]] <- d$IL10_e99f652_UMI[ok]
        }

    } else {
        f <- file.path(STAGE, paste0(ss, "_IL10_molecules.tsv"))
        d <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE)

        if (nrow(d) > 0) {
            agg <- aggregate(
                rep(1L, nrow(d)),
                by=list(barcode_core=d$barcode_core, probe=d$probe),
                FUN=sum
            )
            colnames(agg)[3] <- "UMI"

            for (P in PROBES) {
                z <- agg[agg$probe == P, , drop=FALSE]
                if (nrow(z) == 0) next
                j <- match(z$barcode_core, cell$barcode_core)
                ok <- !is.na(j)
                target_col <- paste0("IL10_", P, "_UMI")
                cell[[target_col]][j[ok]] <- z$UMI[ok]
            }
        }
    }

    accepted <- setNames(
        rep(FALSE, length(PROBES)),
        PROBES
    )

    pm <- probe_meta[probe_meta$rds_sample == S, , drop=FALSE]
    for (P in PROBES) {
        z <- pm[pm$short_probe == P, , drop=FALSE]
        if (nrow(z) != 1) stop("Probe metadata ambiguity: ", S, " ", P)
        accepted[P] <- isTRUE(z$filtered_probes[1])
    }

    for (P in PROBES) {
        cell[[paste0("accepted_", P)]] <- accepted[P]
    }

    probe_cols <- paste0("IL10_", PROBES, "_UMI")
    acc_mat <- as.matrix(cell[, probe_cols, drop=FALSE])

    acc_flag <- matrix(
        rep(as.numeric(accepted), each=nrow(cell)),
        nrow=nrow(cell),
        ncol=length(PROBES)
    )

    cell$accepted_IL10_UMI <- rowSums(acc_mat * acc_flag)
    cell$filteredout_IL10_UMI <- rowSums(acc_mat * (1 - acc_flag))
    cell$all3_IL10_probe_UMI <- rowSums(acc_mat)

    cell$n_positive_probe_all <- rowSums(acc_mat > 0)
    cell$n_positive_probe_accepted <- rowSums((acc_mat > 0) & (acc_flag > 0))

    cell$multi_probe_all <- cell$n_positive_probe_all >= 2
    cell$multi_probe_accepted <- cell$n_positive_probe_accepted >= 2

    cell$gene_vs_accepted_concordant <-
        cell$gene_IL10_UMI_RDS == cell$accepted_IL10_UMI

    cell_list[[S]] <- cell
}

cell_df <- do.call(rbind, cell_list)
rownames(cell_df) <- NULL

write.table(
    cell_df,
    file=file.path(TAB, "IL10_cell_level_all_RDS.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

positive_df <- cell_df[
    cell_df$all3_IL10_probe_UMI > 0 |
    cell_df$gene_IL10_UMI_RDS > 0,
    ,
    drop=FALSE
]

write.table(
    positive_df,
    file=file.path(TAB, "IL10_cell_level_probe_or_gene_positive.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

sample_summary <- do.call(
    rbind,
    lapply(SAMPLES, function(S) {
        x <- cell_df[cell_df$sample == S, , drop=FALSE]

        accepted_status <- paste(
            paste0(
                PROBES, "=",
                c(
                    unique(x$accepted_6b6892b),
                    unique(x$accepted_b1bd8e2),
                    unique(x$accepted_e99f652)
                )
            ),
            collapse=";"
        )

        strongest_idx <- which.max(x$accepted_IL10_UMI)
        strongest_cell <- if (max(x$accepted_IL10_UMI) > 0) x$rds_cell[strongest_idx] else NA_character_
        strongest_lineage <- if (max(x$accepted_IL10_UMI) > 0) x$lineage[strongest_idx] else NA_character_
        strongest_state <- if (max(x$accepted_IL10_UMI) > 0) x$state[strongest_idx] else NA_character_

        data.frame(
            sample=S,
            condition=unique(x$condition),
            transgene_status=unique(x$transgene_status),
            source_delivery=unique(x$source_delivery),
            n_RDS_cells=nrow(x),
            gene_IL10_total_UMI=sum(x$gene_IL10_UMI_RDS),
            gene_IL10_positive_cells=sum(x$gene_IL10_UMI_RDS > 0),
            accepted_probe_IL10_total_UMI=sum(x$accepted_IL10_UMI),
            accepted_probe_IL10_positive_cells=sum(x$accepted_IL10_UMI > 0),
            filteredout_probe_IL10_total_UMI=sum(x$filteredout_IL10_UMI),
            filteredout_only_positive_cells=sum(
                x$filteredout_IL10_UMI > 0 & x$accepted_IL10_UMI == 0
            ),
            all3_probe_IL10_total_UMI=sum(x$all3_IL10_probe_UMI),
            all3_probe_positive_cells=sum(x$all3_IL10_probe_UMI > 0),
            multi_accepted_probe_cells=sum(x$multi_probe_accepted),
            max_gene_IL10_UMI=max(x$gene_IL10_UMI_RDS),
            max_accepted_probe_IL10_UMI=max(x$accepted_IL10_UMI),
            strongest_cell=strongest_cell,
            strongest_lineage=strongest_lineage,
            strongest_state=strongest_state,
            probe_filter_status=accepted_status,
            cellwise_gene_vs_accepted_concordance=
                sprintf("%.4f", mean(x$gene_vs_accepted_concordant)),
            stringsAsFactors=FALSE
        )
    })
)

write.table(
    sample_summary,
    file=file.path(TAB, "IL10_sample_summary.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

lineage_key <- interaction(
    cell_df$sample,
    cell_df$lineage,
    drop=TRUE,
    lex.order=TRUE
)

split_cells <- split(cell_df, lineage_key)

lineage_summary <- do.call(
    rbind,
    lapply(split_cells, function(x) {
        data.frame(
            sample=unique(x$sample),
            condition=unique(x$condition),
            lineage=unique(x$lineage),
            n_cells=nrow(x),
            gene_IL10_total_UMI=sum(x$gene_IL10_UMI_RDS),
            accepted_IL10_total_UMI=sum(x$accepted_IL10_UMI),
            all3_probe_IL10_total_UMI=sum(x$all3_IL10_probe_UMI),
            gene_IL10_positive_cells=sum(x$gene_IL10_UMI_RDS > 0),
            accepted_IL10_positive_cells=sum(x$accepted_IL10_UMI > 0),
            max_accepted_IL10_UMI=max(x$accepted_IL10_UMI),
            multi_accepted_probe_cells=sum(x$multi_probe_accepted),
            stringsAsFactors=FALSE
        )
    })
)

write.table(
    lineage_summary,
    file=file.path(TAB, "IL10_lineage_summary.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

validation <- data.frame(
    sample=SAMPLES,
    gene_total=sample_summary$gene_IL10_total_UMI,
    accepted_probe_total=sample_summary$accepted_probe_IL10_total_UMI,
    exact_total_match=
        sample_summary$gene_IL10_total_UMI ==
        sample_summary$accepted_probe_IL10_total_UMI,
    n_cellwise_mismatches=sapply(
        SAMPLES,
        function(S) {
            x <- cell_df[cell_df$sample == S, , drop=FALSE]
            sum(!x$gene_vs_accepted_concordant)
        }
    ),
    stringsAsFactors=FALSE
)

write.table(
    validation,
    file=file.path(TAB, "IL10_gene_vs_probe_validation.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

if (!all(validation$exact_total_match) || any(validation$n_cellwise_mismatches != 0)) {
    stop(
        "Gene-level IL10 and sample-specific accepted probe counts are not exactly concordant. ",
        "See IL10_gene_vs_probe_validation.tsv"
    )
}

sample_summary$sample <- factor(
    sample_summary$sample,
    levels=SAMPLES
)

stack_df <- rbind(
    data.frame(
        sample=sample_summary$sample,
        signal_class="Cell Ranger accepted probe",
        UMI=sample_summary$accepted_probe_IL10_total_UMI
    ),
    data.frame(
        sample=sample_summary$sample,
        signal_class="Filtered-out probe only",
        UMI=sample_summary$filteredout_probe_IL10_total_UMI
    )
)

p1 <- ggplot(
    stack_df,
    aes(x=sample, y=UMI, fill=signal_class)
) +
    geom_col() +
    geom_text(
        data=sample_summary,
        aes(
            x=sample,
            y=accepted_probe_IL10_total_UMI +
              filteredout_probe_IL10_total_UMI,
            label=paste0("gene=", gene_IL10_total_UMI)
        ),
        inherit.aes=FALSE,
        vjust=-0.4,
        size=3.5
    ) +
    labs(
        title="IL10 probe-level signal by sample",
        subtitle="Accepted = sample-specific Cell Ranger filtered_probes == TRUE",
        x=NULL,
        y="IL10 UMI",
        fill=NULL
    ) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom")

ggsave(
    file.path(FIG, "Fig01_IL10_sample_accepted_vs_filteredout_UMI.pdf"),
    p1,
    width=8,
    height=5
)
ggsave(
    file.path(FIG, "Fig01_IL10_sample_accepted_vs_filteredout_UMI.png"),
    p1,
    width=8,
    height=5,
    dpi=300
)

plot_cells <- positive_df
plot_cells$sample <- factor(plot_cells$sample, levels=SAMPLES)
plot_cells$signal_class <- ifelse(
    plot_cells$accepted_IL10_UMI > 0,
    "Accepted IL10 signal",
    "Filtered-out-only raw signal"
)

p2 <- ggplot(
    plot_cells,
    aes(
        x=sample,
        y=all3_IL10_probe_UMI,
        shape=signal_class,
        size=pmax(n_positive_probe_accepted, 1)
    )
) +
    geom_jitter(width=0.15, height=0) +
    scale_size_continuous(
        breaks=c(1,2),
        name="Accepted probes detected"
    ) +
    labs(
        title="Cell-level IL10 probe signal strength",
        subtitle="Each point is an RDS cell with gene-level or raw probe-level IL10 signal",
        x=NULL,
        y="All-3-probe IL10 UMI",
        shape=NULL
    ) +
    theme_bw(base_size=11) +
    theme(legend.position="bottom")

ggsave(
    file.path(FIG, "Fig02_IL10_cell_signal_strength.pdf"),
    p2,
    width=8,
    height=5.5
)
ggsave(
    file.path(FIG, "Fig02_IL10_cell_signal_strength.png"),
    p2,
    width=8,
    height=5.5,
    dpi=300
)

heat <- positive_df
heat$cell_label <- paste0(
    heat$sample, " | ",
    substr(heat$barcode_core, pmax(1, nchar(heat$barcode_core)-9), nchar(heat$barcode_core)),
    " | ", heat$lineage,
    " | gene=", heat$gene_IL10_UMI_RDS
)

heat_long <- do.call(
    rbind,
    lapply(PROBES, function(P) {
        data.frame(
            sample=heat$sample,
            cell_label=heat$cell_label,
            probe=P,
            UMI=heat[[paste0("IL10_", P, "_UMI")]],
            accepted=heat[[paste0("accepted_", P)]],
            stringsAsFactors=FALSE
        )
    })
)

heat_long$probe <- factor(heat_long$probe, levels=PROBES)

cell_order <- unique(heat$cell_label)
heat_long$cell_label <- factor(
    heat_long$cell_label,
    levels=rev(cell_order)
)

p3 <- ggplot(
    heat_long,
    aes(x=probe, y=cell_label, fill=UMI, alpha=accepted)
) +
    geom_tile() +
    geom_text(
        aes(label=ifelse(UMI > 0, UMI, "")),
        size=3
    ) +
    scale_alpha_manual(
        values=c(`TRUE`=1, `FALSE`=0.35),
        name="Cell Ranger accepted"
    ) +
    facet_grid(
        sample ~ .,
        scales="free_y",
        space="free_y"
    ) +
    labs(
        title="IL10 probe composition of IL10-audited cells",
        subtitle="Tile alpha indicates sample-specific filtered_probes status",
        x="IL10 probe",
        y=NULL,
        fill="UMI"
    ) +
    theme_bw(base_size=9) +
    theme(
        legend.position="bottom",
        axis.text.y=element_text(size=7)
    )

heat_h <- max(7, 0.32 * nrow(heat) + 3)

ggsave(
    file.path(FIG, "Fig03_IL10_cell_probe_heatmap.pdf"),
    p3,
    width=11,
    height=heat_h,
    limitsize=FALSE
)
ggsave(
    file.path(FIG, "Fig03_IL10_cell_probe_heatmap.png"),
    p3,
    width=11,
    height=heat_h,
    dpi=300,
    limitsize=FALSE
)

lin_plot <- lineage_summary[
    lineage_summary$accepted_IL10_total_UMI > 0 |
    lineage_summary$all3_probe_IL10_total_UMI > 0,
    ,
    drop=FALSE
]

lin_plot$sample <- factor(lin_plot$sample, levels=SAMPLES)

p4 <- ggplot(
    lin_plot,
    aes(
        x=lineage,
        y=accepted_IL10_total_UMI,
        fill=sample
    )
) +
    geom_col(position="dodge") +
    labs(
        title="Accepted IL10 probe signal by lineage",
        subtitle="Descriptive only; no cell-level inferential testing",
        x=NULL,
        y="Accepted IL10 UMI",
        fill="Sample"
    ) +
    theme_bw(base_size=10) +
    theme(
        axis.text.x=element_text(angle=35, hjust=1),
        legend.position="bottom"
    )

ggsave(
    file.path(FIG, "Fig04_IL10_accepted_UMI_by_lineage.pdf"),
    p4,
    width=10,
    height=5.5
)
ggsave(
    file.path(FIG, "Fig04_IL10_accepted_UMI_by_lineage.png"),
    p4,
    width=10,
    height=5.5,
    dpi=300
)

summary_path <- file.path(OUT, "IL10_probe_integrated_summary.md")
con <- file(summary_path, open="wt")

writeLines(c(
    "# Human MASH v1.4.2B.3 integrated IL10 probe audit",
    "",
    "Technical implementation: v1.4.2B.3.2",
    "",
    "## Analysis definition",
    "",
    "- SfLB1/SfLB2: transgene-negative iPSC-LO controls from RH241011956.",
    "- Tx17/Tx5: IL10-transgene-positive graft samples from RH251217117.",
    "- `accepted probe` is defined strictly as `filtered_probes == TRUE` in the corresponding Cell Ranger sample output.",
    "- The analysis is sample-specific: a probe can be accepted in one delivery/sample and filtered out in another.",
    "- The IL10 probes are not transgene-specific; they cannot distinguish endogenous IL10 from IL10 transgene transcript by sequence alone.",
    "",
    "## Probe status by sample",
    ""
), con)

for (S in SAMPLES) {
    z <- probe_meta[probe_meta$rds_sample == S, , drop=FALSE]
    writeLines(paste0(
        "- ", S, ": ",
        paste(
            paste0(
                z$short_probe,
                " region=", z$region,
                " filtered_probes=", z$filtered_probes
            ),
            collapse="; "
        )
    ), con)
}

writeLines(c(
    "",
    "## Sample-level results",
    ""
), con)

for (i in seq_len(nrow(sample_summary))) {
    x <- sample_summary[i, ]
    writeLines(paste0(
        "- ", as.character(x$sample),
        ": n=", x$n_RDS_cells,
        "; gene-level IL10=", x$gene_IL10_total_UMI, " UMI / ",
        x$gene_IL10_positive_cells, " cells",
        "; accepted-probe IL10=", x$accepted_probe_IL10_total_UMI, " UMI / ",
        x$accepted_probe_IL10_positive_cells, " cells",
        "; filtered-out probe UMI=", x$filteredout_probe_IL10_total_UMI,
        "; multi-accepted-probe cells=", x$multi_accepted_probe_cells,
        "; max accepted IL10=", x$max_accepted_probe_IL10_UMI
    ), con)
}

tx17 <- cell_df[cell_df$sample == "Tx17", , drop=FALSE]
if (max(tx17$accepted_IL10_UMI) > 0) {
    z <- tx17[which.max(tx17$accepted_IL10_UMI), , drop=FALSE]
    writeLines(c(
        "",
        "## Strongest Tx17 IL10 cell",
        "",
        paste0("- Cell: ", z$rds_cell),
        paste0("- Lineage: ", z$lineage),
        paste0("- State: ", z$state),
        paste0("- Gene-level IL10: ", z$gene_IL10_UMI_RDS, " UMI"),
        paste0(
            "- Probe support: 6b6892b=", z$IL10_6b6892b_UMI,
            ", b1bd8e2=", z$IL10_b1bd8e2_UMI,
            ", e99f652=", z$IL10_e99f652_UMI
        ),
        paste0("- Accepted independent probes detected: ", z$n_positive_probe_accepted),
        ""
    ), con)
}

writeLines(c(
    "## Validation",
    "",
    "- For every sample and every RDS cell, the sum of sample-specific accepted IL10 probe UMIs must exactly equal the existing RDS gene-level IL10 count.",
    "- The script stops if this validation fails.",
    "",
    "## Interpretation boundary",
    "",
    "The comparison can support concordance between the experimental transgene design and the observed IL10 signal pattern. It does not, by itself, prove that an IL10 molecule originated from the transgene because the IL10 probes target IL10 sequence shared with endogenous IL10.",
    ""
), con)

close(con)

cat("\n============================================================\n")
cat("FINAL SAMPLE SUMMARY\n")
cat("============================================================\n")
print(sample_summary, row.names=FALSE)

cat("\n============================================================\n")
cat("GENE vs ACCEPTED-PROBE VALIDATION\n")
cat("============================================================\n")
print(validation, row.names=FALSE)

cat("\n============================================================\n")
cat("OUTPUTS\n")
cat("============================================================\n")
cat("Tables:", TAB, "\n")
cat("Figures:", FIG, "\n")
cat("Summary:", summary_path, "\n")

cat("\n============================================================\n")
cat("v1.4.2B.3.2 INTEGRATED IL10 PROBE AUDIT COMPLETE\n")
cat("============================================================\n")
RS

cat > "$OUT/VERSION.txt" <<EOF
Biological analysis version: v1.4.2B.3
Technical corrected implementation: v1.4.2B.3.2
SfLB source: RH241011956
Tx source: RH251217117
Probe acceptance definition: sample-specific Cell Ranger filtered_probes == TRUE
EOF

echo
echo "COMPLETE:"
echo "$OUT"
