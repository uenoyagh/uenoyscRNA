#!/usr/bin/env python3
# Human MASH v1.4.3B.2C
# Prioritized human graft ligand molecule-level reconstruction.
#
# No semantic assumption is made about /probes/filtered_probes.
# Literal TRUE/FALSE probe subsets are reconstructed and empirically
# calibrated against the canonical standard filtered MEX counts.
#
# No h5py/numpy dependency is required. A tiny C helper is compiled
# locally with clang/cc to scan the large molecule_info arrays efficiently.

from __future__ import annotations

import csv
import math
import mmap
import os
import re
import shutil
import struct
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "v1.4.3B.2C.3"

BASE = Path("/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk")

OUT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.3B.2C.3"
    / "analysis"
    / "v1.4.3B.2C.3_prioritized_ligand_molecule_reconstruction"
)
TABLE_DIR = OUT / "tables"
TMP_DIR = OUT / "_tmp"

CONTEXT_META = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.3.1"
    / "analysis"
    / "v1.4.3.1_12sample_integrated_context_review"
    / "tables"
    / "v1.4.3_context_UMAP_metadata.tsv"
)

INVENTORY = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.3B.1.2"
    / "analysis"
    / "v1.4.3B.1.2_human_graft_ligand_expression_inventory"
    / "tables"
    / "v1.4.3B.1.2_human_graft_ligand_inventory_ranked.tsv"
)

PROBE_STATUS = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.3B.2B.2"
    / "analysis"
    / "v1.4.3B.2B.2_prioritized_ligand_probe_status_audit"
    / "v1.4.3B.2B.2_candidate_probe_status_long.tsv"
)

H5_FILES = {
    "SfLB1": Path(
        "/Volumes/RH241011956/Deliver_RH241011956/"
        "CellRanger/okmt_pool1_outs/per_sample_outs/SfLB_1/"
        "count/sample_molecule_info.h5"
    ),
    "SfLB2": Path(
        "/Volumes/RH241011956/Deliver_RH241011956/"
        "CellRanger/okmt_pool1_outs/per_sample_outs/SfLB_2/"
        "count/sample_molecule_info.h5"
    ),
    "Tx17": Path(
        "/Volumes/RH251217117/Deliver_RH251217117/"
        "Result_human/01_CellRanger/Count_matrix/Tx17/"
        "molecule_info.h5"
    ),
    "Tx5": Path(
        "/Volumes/RH251217117/Deliver_RH251217117/"
        "Result_human/01_CellRanger/Count_matrix/Tx5/"
        "molecule_info.h5"
    ),
}

CANDIDATES = [
    "IL10",
    "FGF1",
    "VEGFA",
    "CXCL12",
    "TNFSF12",
    "FST",
    "IGFBP2",
    "IGFBP4",
    "INHBA",
    "BMP2",
    "BMP5",
    "TGFB3",
    "HBEGF",
    "VEGFB",
    "BTC",
]

SAMPLES = ["SfLB1", "SfLB2", "Tx17", "Tx5"]

LOW_FILTER_SENSITIVITY_GENES = {
    "VEGFA",
    "IGFBP2",
    "IGFBP4",
    "INHBA",
    "VEGFB",
}

EXPECTED_CANONICAL_CELLS = {
    "SfLB1": 19081,
    "SfLB2": 9805,
    "Tx17": 5958,
    "Tx5": 8571,
}

MODES = [
    "ALL_PROBES",
    "INCLUDED_TRUE",
    "FILTERED_TRUE",
    "INCLUDED_AND_FILTERED_TRUE",
    "INCLUDED_TRUE_FILTERED_FALSE",
    "SHARED_FILTERED_TRUE_ALL4",
]

MIN_FREE_BYTES = 3 * 1024**3

C_HELPER_SOURCE = r'''
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#pragma pack(push, 1)
typedef struct {
    int32_t probe_idx;       /* H5T_STD_I32LE */
    uint64_t barcode_idx;    /* H5T_STD_U64LE */
    uint32_t umi;            /* H5T_STD_U32LE */
} OutRec;
#pragma pack(pop)

#define CHUNK 1048576
#define OUT_CHUNK 65536

static void die(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(2);
}

static long long get_file_size(FILE *f) {
    if (fseek(f, 0, SEEK_END) != 0) die("fseek failed");
    long x = ftell(f);
    if (x < 0) die("ftell failed");
    if (fseek(f, 0, SEEK_SET) != 0) die("rewind failed");
    return (long long)x;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(
            stderr,
            "usage: %s probe.bin barcode.bin umi.bin target_idx.txt out.bin\n",
            argv[0]
        );
        return 2;
    }

    if (sizeof(OutRec) != 16) {
        die("unexpected OutRec size; expected 16 bytes");
    }

    FILE *fp = fopen(argv[1], "rb");
    FILE *fb = fopen(argv[2], "rb");
    FILE *fu = fopen(argv[3], "rb");
    FILE *ft = fopen(argv[4], "r");
    FILE *fo = fopen(argv[5], "wb");

    if (!fp || !fb || !fu || !ft || !fo) {
        die("could not open one or more files");
    }

    long long sp = get_file_size(fp);
    long long sb = get_file_size(fb);
    long long su = get_file_size(fu);

    if ((sp % 4) != 0) die("probe_idx binary size is not divisible by 4");
    if ((sb % 8) != 0) die("barcode_idx binary size is not divisible by 8");
    if ((su % 4) != 0) die("umi binary size is not divisible by 4");

    unsigned long long np_total = (unsigned long long)(sp / 4);
    unsigned long long nb_total = (unsigned long long)(sb / 8);
    unsigned long long nu_total = (unsigned long long)(su / 4);

    if (np_total != nb_total || np_total != nu_total) {
        die("numeric arrays have different element counts");
    }

    int value = 0;
    int max_target = -1;

    while (fscanf(ft, "%d", &value) == 1) {
        if (value > max_target) max_target = value;
    }

    if (max_target < 0) die("empty target probe list");

    unsigned char *mask = calloc(
        (size_t)max_target + 1,
        sizeof(unsigned char)
    );

    if (!mask) die("could not allocate target mask");

    rewind(ft);

    while (fscanf(ft, "%d", &value) == 1) {
        if (value >= 0 && value <= max_target) {
            mask[value] = 1;
        }
    }

    int32_t *pbuf = malloc(CHUNK * sizeof(int32_t));
    uint64_t *bbuf = malloc(CHUNK * sizeof(uint64_t));
    uint32_t *ubuf = malloc(CHUNK * sizeof(uint32_t));
    OutRec *obuf = malloc(OUT_CHUNK * sizeof(OutRec));

    if (!pbuf || !bbuf || !ubuf || !obuf) {
        die("could not allocate scan buffers");
    }

    size_t out_n = 0;
    unsigned long long scanned = 0;
    unsigned long long selected = 0;

    while (1) {
        size_t np = fread(pbuf, sizeof(int32_t), CHUNK, fp);
        size_t nb = fread(bbuf, sizeof(uint64_t), CHUNK, fb);
        size_t nu = fread(ubuf, sizeof(uint32_t), CHUNK, fu);

        if (np != nb || np != nu) {
            die("numeric arrays ended at different positions");
        }

        if (np == 0) break;

        for (size_t i = 0; i < np; ++i) {
            int32_t p = pbuf[i];

            if (p >= 0 && p <= max_target && mask[p]) {
                obuf[out_n].probe_idx = p;
                obuf[out_n].barcode_idx = bbuf[i];
                obuf[out_n].umi = ubuf[i];
                out_n++;
                selected++;

                if (out_n == OUT_CHUNK) {
                    if (fwrite(obuf, sizeof(OutRec), out_n, fo) != out_n) {
                        die("failed writing selected records");
                    }
                    out_n = 0;
                }
            }
        }

        scanned += np;
    }

    if (out_n > 0) {
        if (fwrite(obuf, sizeof(OutRec), out_n, fo) != out_n) {
            die("failed writing final records");
        }
    }

    printf(
        "SCANNED=%llu SELECTED=%llu ELEMENTS=%llu\n",
        scanned,
        selected,
        np_total
    );

    free(mask);
    free(pbuf);
    free(bbuf);
    free(ubuf);
    free(obuf);

    fclose(fp);
    fclose(fb);
    fclose(fu);
    fclose(ft);
    fclose(fo);

    return 0;
}
'''


def fail(msg: str) -> None:
    raise RuntimeError(msg)


def norm_sample(x: str) -> str:
    x = x.strip()
    if x == "SfLB_1":
        return "SfLB1"
    if x == "SfLB_2":
        return "SfLB2"
    return x


def bool_text(x: str) -> bool:
    return str(x).strip().upper() == "TRUE"


def extract_barcode_core(x: str) -> str | None:
    hits = re.findall(r"[ACGT]{16,40}", x.upper())
    if not hits:
        return None
    return max(hits, key=len)


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def find_executable(names: list[str], preferred: list[str] | None = None) -> str:
    if preferred:
        for p in preferred:
            if Path(p).exists():
                return p
    for name in names:
        x = shutil.which(name)
        if x:
            return x
    fail("Executable not found: " + "/".join(names))


def compile_helper(compiler: str) -> Path:
    c_path = TMP_DIR / "molecule_target_scan.c"
    exe_path = TMP_DIR / "molecule_target_scan"
    c_path.write_text(C_HELPER_SOURCE, encoding="utf-8")

    p = subprocess.run(
        [compiler, "-O3", "-std=c11", str(c_path), "-o", str(exe_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if p.returncode != 0:
        fail("C helper compilation failed:\n" + p.stderr)

    return exe_path


def h5dump_binary(h5dump: str, h5: Path, dataset: str, out_path: Path) -> None:
    if out_path.exists():
        out_path.unlink()

    print(f"  exporting {dataset}")

    p = subprocess.run(
        [
            h5dump,
            "-d",
            dataset,
            "-b",
            "LE",
            "-o",
            str(out_path),
            str(h5),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )

    if p.returncode != 0:
        fail(f"h5dump failed for {dataset}:\n{p.stderr}")

    if not out_path.exists() or out_path.stat().st_size == 0:
        fail(f"Empty binary export for {dataset}")


def load_context() -> dict[str, dict[str, dict]]:
    rows = read_tsv(CONTEXT_META)
    out = {s: {} for s in SAMPLES}

    required = {
        "sample",
        "barcode_core",
        "context_annotation",
        "context_state",
    }

    if not rows:
        fail("Context metadata is empty")

    missing = required - set(rows[0])
    if missing:
        fail("Context metadata missing: " + ", ".join(sorted(missing)))

    for r in rows:
        s = norm_sample(r["sample"])
        if s not in SAMPLES:
            continue

        bc = r["barcode_core"].strip()

        if bc in out[s]:
            fail(f"{s}: duplicated canonical barcode {bc}")

        out[s][bc] = {
            "annotation": r["context_annotation"],
            "state": r["context_state"],
        }

    for s in SAMPLES:
        if len(out[s]) != EXPECTED_CANONICAL_CELLS[s]:
            fail(
                f"{s}: canonical cells={len(out[s])}, "
                f"expected={EXPECTED_CANONICAL_CELLS[s]}"
            )

    return out


def load_inventory() -> dict[str, dict]:
    rows = read_tsv(INVENTORY)
    out = {}

    for r in rows:
        g = (r.get("gene") or "").strip().upper()
        if g in CANDIDATES:
            out[g] = r

    missing = [g for g in CANDIDATES if g not in out]
    if missing:
        fail("Inventory missing: " + ", ".join(missing))

    return out


def load_probe_status():
    rows = read_tsv(PROBE_STATUS)

    required = {
        "sample",
        "gene",
        "probe_index_0based",
        "probe_id",
        "included_probes",
        "filtered_probes",
    }

    if not rows:
        fail("Probe-status table is empty")

    missing = required - set(rows[0])
    if missing:
        fail("Probe-status table missing: " + ", ".join(sorted(missing)))

    by_sample = {s: {} for s in SAMPLES}
    true_sets = {
        s: {g: set() for g in CANDIDATES}
        for s in SAMPLES
    }

    for r in rows:
        s = norm_sample(r["sample"])
        g = r["gene"].strip().upper()

        if s not in SAMPLES or g not in CANDIDATES:
            continue

        idx = int(r["probe_index_0based"])

        if idx in by_sample[s]:
            fail(f"{s}: duplicate probe_idx={idx}")

        inc = bool_text(r["included_probes"])
        filt = bool_text(r["filtered_probes"])

        meta = {
            "sample": s,
            "gene": g,
            "probe_idx": idx,
            "probe_id": r["probe_id"].strip(),
            "included": inc,
            "filtered": filt,
        }

        by_sample[s][idx] = meta

        if filt:
            true_sets[s][g].add(meta["probe_id"])

    for s in SAMPLES:
        seen = {x["gene"] for x in by_sample[s].values()}
        missing_genes = [g for g in CANDIDATES if g not in seen]
        if missing_genes:
            fail(f"{s}: missing candidate probe definitions: " + ", ".join(missing_genes))

    shared_true = {
        g: set.intersection(*[true_sets[s][g] for s in SAMPLES])
        for g in CANDIDATES
    }

    return by_sample, shared_true


def export_sample_arrays(
    h5dump: str,
    sample: str,
    h5: Path,
    sample_tmp: Path,
) -> dict:
    """
    Export only the large numeric molecule arrays as raw little-endian binary.

    IMPORTANT:
    /barcodes is a fixed-width HDF5 string dataset. Some h5dump builds do not
    emit a usable raw binary file for string datasets with "-b LE". Therefore
    v1.4.3B.2C.1 deliberately does NOT binary-export /barcodes. Barcode strings
    are resolved later by streaming the normal text h5dump output and retaining
    only barcode indices actually observed among the selected ligand molecules.
    """
    sample_tmp.mkdir(parents=True, exist_ok=True)

    if shutil.disk_usage(sample_tmp).free < MIN_FREE_BYTES:
        fail(f"{sample}: less than 3 GiB free space for temporary binary export")

    paths = {
        "probe": sample_tmp / "probe_idx.bin",
        "barcode": sample_tmp / "barcode_idx.bin",
        "umi": sample_tmp / "umi.bin",
    }

    h5dump_binary(h5dump, h5, "/probe_idx", paths["probe"])
    h5dump_binary(h5dump, h5, "/barcode_idx", paths["barcode"])
    h5dump_binary(h5dump, h5, "/umi", paths["umi"])

    probe_bytes = paths["probe"].stat().st_size
    barcode_bytes = paths["barcode"].stat().st_size
    umi_bytes = paths["umi"].stat().st_size

    # Datatypes were directly audited in all four source HDF5 files:
    #   /probe_idx   H5T_STD_I32LE  -> 4 bytes/element
    #   /barcode_idx H5T_STD_U64LE  -> 8 bytes/element
    #   /umi         H5T_STD_U32LE  -> 4 bytes/element
    if probe_bytes % 4 != 0:
        fail(f"{sample}: probe_idx byte size is not divisible by 4")

    if barcode_bytes % 8 != 0:
        fail(f"{sample}: barcode_idx byte size is not divisible by 8")

    if umi_bytes % 4 != 0:
        fail(f"{sample}: umi byte size is not divisible by 4")

    n_probe = probe_bytes // 4
    n_barcode = barcode_bytes // 8
    n_umi = umi_bytes // 4

    if not (n_probe == n_barcode == n_umi):
        fail(
            f"{sample}: numeric array element counts differ: "
            f"probe_idx={n_probe}, barcode_idx={n_barcode}, umi={n_umi}"
        )

    print(
        f"  numeric array elements: {n_probe} "
        f"(probe_idx=int32, barcode_idx=uint64, umi=uint32)"
    )

    return paths


def make_target_file(probe_meta: dict[int, dict], sample_tmp: Path) -> Path:
    p = sample_tmp / "target_probe_indices.txt"
    with p.open("w", encoding="utf-8") as fh:
        for idx in sorted(probe_meta):
            fh.write(f"{idx}\n")
    return p


def scan_selected(
    helper: Path,
    paths: dict,
    target_file: Path,
    sample_tmp: Path,
) -> Path:
    out = sample_tmp / "selected_molecules.bin"

    p = subprocess.run(
        [
            str(helper),
            str(paths["probe"]),
            str(paths["barcode"]),
            str(paths["umi"]),
            str(target_file),
            str(out),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if p.returncode != 0:
        fail("C molecule scanner failed:\n" + p.stderr)

    print("  " + p.stdout.strip())

    if not out.exists() or out.stat().st_size % 16 != 0:
        fail("selected_molecules.bin was not created correctly (expected packed 16-byte records)")

    return out


def collect_selected_barcode_indices(
    selected_file: Path,
) -> set[int]:
    """
    Read the compact selected-molecule file produced by the C scanner and
    return only the barcode_idx values that actually occur in candidate rows.
    """
    wanted = set()
    rec = struct.Struct("<iQI")

    with selected_file.open("rb") as fh:
        while True:
            raw = fh.read(16 * 250000)

            if not raw:
                break

            if len(raw) % 16 != 0:
                fail("Partial selected-molecule record block while collecting barcode indices")

            for _probe_idx, barcode_idx, _umi in rec.iter_unpack(raw):
                wanted.add(barcode_idx)

    return wanted


def decode_h5dump_string(raw: str) -> str:
    """
    Decode h5dump escaped fixed-width strings such as:
      AAAC...-1\\000\\000...
    """
    try:
        x = bytes(raw, "utf-8").decode(
            "unicode_escape",
            errors="replace",
        )
    except Exception:
        x = raw

    return x.replace("\x00", "").strip()


def resolve_selected_barcodes_from_h5dump(
    h5dump: str,
    h5: Path,
    wanted_indices: set[int],
) -> dict[int, str]:
    """
    Stream indexed textual h5dump output for /barcodes and retain only the
    requested barcode indices.

    IMPORTANT:
    Earlier v1.4.3B.2C.2 used h5dump "-y", which explicitly suppresses array
    indices.  That made index-based matching impossible and yielded Found=0.
    This version keeps indices and also handles multiple quoted strings on one
    h5dump line by assigning consecutive indices from the printed base index.
    """
    if not wanted_indices:
        return {}

    print(
        f"  resolving {len(wanted_indices)} selected barcode indices "
        "from indexed textual /barcodes h5dump"
    )

    cmd = [
        h5dump,
        "-w",
        "120",
        "-d",
        "/barcodes",
        str(h5),
    ]

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    if p.stdout is None or p.stderr is None:
        fail("Could not open h5dump pipes for /barcodes")

    found = {}

    # Typical h5dump data lines may be either:
    #   (123): "BARCODE"
    # or:
    #   (123): "BARCODE1", "BARCODE2", ...
    # In the latter case, quoted values correspond to consecutive indices.
    index_pat = re.compile(r'^\s*\((\d+)\):\s*(.*)$')
    string_pat = re.compile(r'"((?:\\.|[^"\\])*)"')

    for line in p.stdout:
        m = index_pat.match(line)

        if not m:
            continue

        base_idx = int(m.group(1))
        payload = m.group(2)

        vals = string_pat.findall(payload)

        for offset, raw_value in enumerate(vals):
            idx = base_idx + offset

            if idx not in wanted_indices:
                continue

            found[idx] = decode_h5dump_string(raw_value)

        if len(found) == len(wanted_indices):
            # All required indices have been resolved; stop early.
            p.terminate()
            break

    stderr = p.stderr.read()

    try:
        rc = p.wait(timeout=10)
    except subprocess.TimeoutExpired:
        p.kill()
        rc = p.wait()

    if len(found) != len(wanted_indices):
        missing = sorted(wanted_indices - set(found))

        fail(
            f"Could not resolve all requested barcode indices from {h5}. "
            f"Found={len(found)} wanted={len(wanted_indices)} "
            f"missing_first={missing[:20]}\n"
            f"h5dump return_code={rc}\n"
            f"h5dump stderr:\n{stderr}"
        )

    print(
        f"  PASS: resolved all {len(found)} selected barcode indices"
    )

    return found


def mex_values(inventory: dict[str, dict], gene: str, sample: str) -> tuple[int, int]:
    row = inventory[gene]

    try:
        umi = int(round(float(row[f"{sample}_UMI"])))
        pos = int(round(float(row[f"{sample}_positive_cells"])))
    except Exception as e:
        fail(f"Could not parse MEX values for {sample} {gene}: {e}")

    return umi, pos


def mol_key(barcode_idx: int, umi: int) -> int:
    return (barcode_idx << 32) | umi


def summarize_mode(
    sample: str,
    gene: str,
    mode: str,
    keys: set[int],
    n_canonical: int,
) -> dict:
    cells = {k >> 32 for k in keys}

    return {
        "sample": sample,
        "gene": gene,
        "mode": mode,
        "UMI": len(keys),
        "positive_cells": len(cells),
        "positive_fraction": len(cells) / n_canonical,
        "UMI_per_1000_cells": len(keys) / n_canonical * 1000.0,
    }


def aggregate_group(
    sample: str,
    gene: str,
    mode: str,
    keys: set[int],
    canonical_by_idx: dict[int, dict],
    field: str,
) -> list[dict]:
    denom = defaultdict(set)
    umi_sets = defaultdict(set)
    pos_sets = defaultdict(set)

    for idx, meta in canonical_by_idx.items():
        denom[meta[field]].add(idx)

    for key in keys:
        idx = key >> 32
        grp = canonical_by_idx[idx][field]
        umi_sets[grp].add(key)
        pos_sets[grp].add(idx)

    rows = []

    for grp in sorted(denom):
        n_cells = len(denom[grp])
        n_umi = len(umi_sets[grp])
        n_pos = len(pos_sets[grp])

        rows.append(
            {
                "sample": sample,
                "gene": gene,
                "mode": mode,
                "group": grp,
                "n_cells": n_cells,
                "UMI": n_umi,
                "positive_cells": n_pos,
                "positive_fraction": n_pos / n_cells if n_cells else 0.0,
                "UMI_per_1000_cells": n_umi / n_cells * 1000.0 if n_cells else 0.0,
            }
        )

    return rows


def process_sample(
    sample: str,
    h5: Path,
    h5dump: str,
    helper: Path,
    probe_meta: dict[int, dict],
    shared_true: dict[str, set[str]],
    canonical_meta: dict[str, dict],
    inventory: dict[str, dict],
):
    print()
    print("-" * 60)
    print(f"Processing {sample}")
    print("-" * 60)
    print(f"H5: {h5}")

    sample_tmp = TMP_DIR / sample

    if sample_tmp.exists():
        shutil.rmtree(sample_tmp)

    sample_tmp.mkdir(parents=True, exist_ok=True)

    keep_temp = os.environ.get("KEEP_V143B2C_TEMP", "") == "1"

    try:
        paths = export_sample_arrays(h5dump, sample, h5, sample_tmp)
        target_file = make_target_file(probe_meta, sample_tmp)
        selected_file = scan_selected(helper, paths, target_file, sample_tmp)

        wanted_barcode_indices = collect_selected_barcode_indices(
            selected_file
        )

        barcode_strings = resolve_selected_barcodes_from_h5dump(
            h5dump,
            h5,
            wanted_barcode_indices,
        )

        barcode_cache = {}
        canonical_by_idx = {}

        mode_sets = {
            g: {mode: set() for mode in MODES}
            for g in CANDIDATES
        }

        probe_sets = defaultdict(set)

        canonical_selected_rows = 0
        background_selected_rows = 0

        rec = struct.Struct("<iQI")

        with selected_file.open("rb") as fh:
            while True:
                raw = fh.read(16 * 250000)

                if not raw:
                    break

                if len(raw) % 16 != 0:
                    fail(f"{sample}: partial selected-record block")

                for probe_idx, barcode_idx, umi in rec.iter_unpack(raw):
                    pm = probe_meta.get(probe_idx)

                    if pm is None:
                        fail(f"{sample}: unknown selected probe_idx={probe_idx}")

                    cached = barcode_cache.get(barcode_idx)

                    if cached is None:
                        raw_bc = barcode_strings.get(barcode_idx)

                        if raw_bc is None:
                            fail(
                                f"{sample}: selected barcode_idx={barcode_idx} "
                                "was not resolved from /barcodes"
                            )

                        core = extract_barcode_core(raw_bc)
                        meta = canonical_meta.get(core) if core else None
                        cached = (core, meta)
                        barcode_cache[barcode_idx] = cached

                        if meta is not None:
                            canonical_by_idx[barcode_idx] = meta

                    core, meta = cached

                    if meta is None:
                        background_selected_rows += 1
                        continue

                    canonical_selected_rows += 1

                    g = pm["gene"]
                    pid = pm["probe_id"]
                    key = mol_key(barcode_idx, umi)

                    mode_sets[g]["ALL_PROBES"].add(key)

                    if pm["included"]:
                        mode_sets[g]["INCLUDED_TRUE"].add(key)

                    if pm["filtered"]:
                        mode_sets[g]["FILTERED_TRUE"].add(key)

                    if pm["included"] and pm["filtered"]:
                        mode_sets[g]["INCLUDED_AND_FILTERED_TRUE"].add(key)

                    if pm["included"] and not pm["filtered"]:
                        mode_sets[g]["INCLUDED_TRUE_FILTERED_FALSE"].add(key)

                    if pid in shared_true[g]:
                        mode_sets[g]["SHARED_FILTERED_TRUE_ALL4"].add(key)

                    probe_sets[(g, pid)].add(key)

        if not canonical_by_idx:
            fail(f"{sample}: no target molecules mapped to canonical cells")

        n_canonical = len(canonical_meta)

        sample_rows = []
        probe_rows = []
        lineage_rows = []
        state_rows = []

        for g in CANDIDATES:
            mex_umi, mex_pos = mex_values(inventory, g, sample)

            mode_sum = {}

            for mode in MODES:
                sm = summarize_mode(
                    sample,
                    g,
                    mode,
                    mode_sets[g][mode],
                    n_canonical,
                )

                mode_sum[mode] = sm

                lineage_rows.extend(
                    aggregate_group(
                        sample,
                        g,
                        mode,
                        mode_sets[g][mode],
                        canonical_by_idx,
                        "annotation",
                    )
                )

                state_rows.extend(
                    aggregate_group(
                        sample,
                        g,
                        mode,
                        mode_sets[g][mode],
                        canonical_by_idx,
                        "state",
                    )
                )

            matching_umi = [
                mode
                for mode in MODES
                if mode_sum[mode]["UMI"] == mex_umi
            ]

            matching_full = [
                mode
                for mode in MODES
                if (
                    mode_sum[mode]["UMI"] == mex_umi
                    and mode_sum[mode]["positive_cells"] == mex_pos
                )
            ]

            row = {
                "sample": sample,
                "gene": g,
                "canonical_cells": n_canonical,
                "MEX_UMI": mex_umi,
                "MEX_positive_cells": mex_pos,
                "matching_UMI_modes": ";".join(matching_umi),
                "matching_UMI_and_positive_cell_modes": ";".join(matching_full),
            }

            for mode in MODES:
                row[f"{mode}_UMI"] = mode_sum[mode]["UMI"]
                row[f"{mode}_positive_cells"] = mode_sum[mode]["positive_cells"]

            sample_rows.append(row)

            gene_probes = [
                x
                for x in probe_meta.values()
                if x["gene"] == g
            ]

            for pm in sorted(gene_probes, key=lambda x: x["probe_idx"]):
                keys = probe_sets[(g, pm["probe_id"])]

                probe_rows.append(
                    {
                        "sample": sample,
                        "gene": g,
                        "probe_idx": pm["probe_idx"],
                        "probe_id": pm["probe_id"],
                        "included_probes": "TRUE" if pm["included"] else "FALSE",
                        "filtered_probes": "TRUE" if pm["filtered"] else "FALSE",
                        "is_shared_filtered_TRUE_all4": (
                            "TRUE" if pm["probe_id"] in shared_true[g] else "FALSE"
                        ),
                        "canonical_UMI": len(keys),
                        "canonical_positive_cells": len({k >> 32 for k in keys}),
                    }
                )

        audit_row = {
            "sample": sample,
            "canonical_cells_expected": n_canonical,
            "canonical_barcode_indices_with_candidate_signal": len(canonical_by_idx),
            "canonical_selected_molecule_rows": canonical_selected_rows,
            "background_selected_molecule_rows": background_selected_rows,
            "selected_barcode_indices_cached": len(barcode_cache),
        }

        return sample_rows, probe_rows, lineage_rows, state_rows, audit_row

    finally:
        if not keep_temp and sample_tmp.exists():
            print(f"  cleaning temporary binaries for {sample}")
            shutil.rmtree(sample_tmp, ignore_errors=True)


def build_cross_summary(sample_rows: list[dict]) -> list[dict]:
    by_gene = defaultdict(dict)

    for r in sample_rows:
        by_gene[r["gene"]][r["sample"]] = r

    rows = []

    for g in CANDIDATES:
        x = by_gene[g]

        if set(x) != set(SAMPLES):
            fail(f"{g}: incomplete four-sample reconstruction")

        strict = [
            x["SfLB1"]["SHARED_FILTERED_TRUE_ALL4_UMI"],
            x["SfLB2"]["SHARED_FILTERED_TRUE_ALL4_UMI"],
            x["Tx17"]["SHARED_FILTERED_TRUE_ALL4_UMI"],
            x["Tx5"]["SHARED_FILTERED_TRUE_ALL4_UMI"],
        ]

        tx_mean = (
            strict[2] / EXPECTED_CANONICAL_CELLS["Tx17"]
            + strict[3] / EXPECTED_CANONICAL_CELLS["Tx5"]
        ) / 2.0

        sflb_mean = (
            strict[0] / EXPECTED_CANONICAL_CELLS["SfLB1"]
            + strict[1] / EXPECTED_CANONICAL_CELLS["SfLB2"]
        ) / 2.0

        pseudo = 1e-6
        log2fc = math.log2((tx_mean + pseudo) / (sflb_mean + pseudo))

        rows.append(
            {
                "gene": g,
                "SfLB1_MEX_UMI": x["SfLB1"]["MEX_UMI"],
                "SfLB2_MEX_UMI": x["SfLB2"]["MEX_UMI"],
                "Tx17_MEX_UMI": x["Tx17"]["MEX_UMI"],
                "Tx5_MEX_UMI": x["Tx5"]["MEX_UMI"],
                "SfLB1_shared_TRUE_UMI": strict[0],
                "SfLB2_shared_TRUE_UMI": strict[1],
                "Tx17_shared_TRUE_UMI": strict[2],
                "Tx5_shared_TRUE_UMI": strict[3],
                "strict_shared_TRUE_Tx_mean_UMI_per_cell": tx_mean,
                "strict_shared_TRUE_SfLB_mean_UMI_per_cell": sflb_mean,
                "strict_shared_TRUE_Tx_vs_SfLB_log2FC_per_cell_descriptive": log2fc,
                "SfLB1_MEX_matching_modes": x["SfLB1"][
                    "matching_UMI_and_positive_cell_modes"
                ],
                "SfLB2_MEX_matching_modes": x["SfLB2"][
                    "matching_UMI_and_positive_cell_modes"
                ],
                "Tx17_MEX_matching_modes": x["Tx17"][
                    "matching_UMI_and_positive_cell_modes"
                ],
                "Tx5_MEX_matching_modes": x["Tx5"][
                    "matching_UMI_and_positive_cell_modes"
                ],
            }
        )

    rows.sort(
        key=lambda r: (
            0 if r["gene"] == "IL10" else 1,
            -r["strict_shared_TRUE_Tx_mean_UMI_per_cell"],
            r["gene"],
        )
    )

    return rows


def main() -> None:
    print("=" * 60)
    print(
        f"Human MASH {VERSION} prioritized ligand molecule-level reconstruction"
    )
    print("=" * 60)
    print()
    print("Literal included_probes / filtered_probes subsets will be")
    print("calibrated empirically against canonical standard MEX.")
    print("Barcode strings are resolved from indexed streamed textual h5dump output;")
    print("raw binary export is used only for numeric molecule arrays.")
    print("Validated numeric datatypes: probe_idx=int32, barcode_idx=uint64, umi=uint32.")
    print()

    required_paths = [
        CONTEXT_META,
        INVENTORY,
        PROBE_STATUS,
        *H5_FILES.values(),
    ]

    for p in required_paths:
        if not p.exists():
            fail(f"Required input missing: {p}")

    OUT.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    h5dump = find_executable(
        ["h5dump"],
        preferred=["/opt/homebrew/bin/h5dump"],
    )

    compiler = find_executable(["clang", "cc", "gcc"])

    print(f"h5dump : {h5dump}")
    print(f"compiler: {compiler}")

    helper = compile_helper(compiler)
    context = load_context()
    inventory = load_inventory()
    probe_status, shared_true = load_probe_status()

    all_sample_rows = []
    all_probe_rows = []
    all_lineage_rows = []
    all_state_rows = []
    all_audit_rows = []

    for sample in SAMPLES:
        result = process_sample(
            sample,
            H5_FILES[sample],
            h5dump,
            helper,
            probe_status[sample],
            shared_true,
            context[sample],
            inventory,
        )

        sample_rows, probe_rows, lineage_rows, state_rows, audit_row = result

        all_sample_rows.extend(sample_rows)
        all_probe_rows.extend(probe_rows)
        all_lineage_rows.extend(lineage_rows)
        all_state_rows.extend(state_rows)
        all_audit_rows.append(audit_row)

    sample_fields = [
        "sample",
        "gene",
        "canonical_cells",
        "MEX_UMI",
        "MEX_positive_cells",
        "matching_UMI_modes",
        "matching_UMI_and_positive_cell_modes",
    ]

    for mode in MODES:
        sample_fields.extend(
            [f"{mode}_UMI", f"{mode}_positive_cells"]
        )

    write_tsv(
        TABLE_DIR / f"{VERSION}_gene_reconstruction_by_sample.tsv",
        all_sample_rows,
        sample_fields,
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_probe_composition_by_sample.tsv",
        all_probe_rows,
        [
            "sample",
            "gene",
            "probe_idx",
            "probe_id",
            "included_probes",
            "filtered_probes",
            "is_shared_filtered_TRUE_all4",
            "canonical_UMI",
            "canonical_positive_cells",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_lineage_reconstruction_long.tsv",
        all_lineage_rows,
        [
            "sample",
            "gene",
            "mode",
            "group",
            "n_cells",
            "UMI",
            "positive_cells",
            "positive_fraction",
            "UMI_per_1000_cells",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_state_reconstruction_long.tsv",
        all_state_rows,
        [
            "sample",
            "gene",
            "mode",
            "group",
            "n_cells",
            "UMI",
            "positive_cells",
            "positive_fraction",
            "UMI_per_1000_cells",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_barcode_mapping_audit.tsv",
        all_audit_rows,
        [
            "sample",
            "canonical_cells_expected",
            "canonical_barcode_indices_with_candidate_signal",
            "canonical_selected_molecule_rows",
            "background_selected_molecule_rows",
            "selected_barcode_indices_cached",
        ],
    )

    cross_rows = build_cross_summary(all_sample_rows)

    write_tsv(
        TABLE_DIR / f"{VERSION}_crosssample_strict_shared_TRUE_summary.tsv",
        cross_rows,
        [
            "gene",
            "SfLB1_MEX_UMI",
            "SfLB2_MEX_UMI",
            "Tx17_MEX_UMI",
            "Tx5_MEX_UMI",
            "SfLB1_shared_TRUE_UMI",
            "SfLB2_shared_TRUE_UMI",
            "Tx17_shared_TRUE_UMI",
            "Tx5_shared_TRUE_UMI",
            "strict_shared_TRUE_Tx_mean_UMI_per_cell",
            "strict_shared_TRUE_SfLB_mean_UMI_per_cell",
            "strict_shared_TRUE_Tx_vs_SfLB_log2FC_per_cell_descriptive",
            "SfLB1_MEX_matching_modes",
            "SfLB2_MEX_matching_modes",
            "Tx17_MEX_matching_modes",
            "Tx5_MEX_matching_modes",
        ],
    )

    il10_rows = [r for r in all_sample_rows if r["gene"] == "IL10"]

    write_tsv(
        TABLE_DIR / f"{VERSION}_IL10_resolution.tsv",
        il10_rows,
        sample_fields,
    )

    calibration_failures = []

    for r in all_sample_rows:
        if r["gene"] not in LOW_FILTER_SENSITIVITY_GENES:
            continue

        if (
            r["INCLUDED_AND_FILTERED_TRUE_UMI"] != r["MEX_UMI"]
            or r["INCLUDED_AND_FILTERED_TRUE_positive_cells"]
            != r["MEX_positive_cells"]
        ):
            calibration_failures.append(r)

    print()
    print("=" * 60)
    print("LOW-FILTER-SENSITIVITY CALIBRATION")
    print("=" * 60)

    if calibration_failures:
        for r in calibration_failures:
            print(
                f"FAIL {r['sample']} {r['gene']}: "
                f"MEX={r['MEX_UMI']}/{r['MEX_positive_cells']}cells "
                f"vs INCLUDED_AND_FILTERED_TRUE="
                f"{r['INCLUDED_AND_FILTERED_TRUE_UMI']}/"
                f"{r['INCLUDED_AND_FILTERED_TRUE_positive_cells']}cells"
            )

        fail(
            "Low-filter-sensitivity calibration did not reproduce canonical MEX. "
            "Debug tables were written; do not interpret mixed-flag genes yet."
        )

    print(
        "PASS: all 5 low-filter-sensitivity genes reproduce canonical MEX "
        "in all 4 samples with INCLUDED_AND_FILTERED_TRUE."
    )

    print()
    print("=" * 60)
    print("IL10 RESOLUTION")
    print("=" * 60)

    for r in il10_rows:
        print(
            f"{r['sample']:5s} "
            f"MEX={r['MEX_UMI']:>3}/{r['MEX_positive_cells']:>3}cells  "
            f"ALL={r['ALL_PROBES_UMI']:>3}  "
            f"INCL={r['INCLUDED_TRUE_UMI']:>3}  "
            f"FILT_TRUE={r['FILTERED_TRUE_UMI']:>3}  "
            f"INCL&TRUE={r['INCLUDED_AND_FILTERED_TRUE_UMI']:>3}  "
            f"INCL&FALSE={r['INCLUDED_TRUE_FILTERED_FALSE_UMI']:>3}  "
            f"SHARED_TRUE={r['SHARED_FILTERED_TRUE_ALL4_UMI']:>3}  "
            f"MEX_MATCH={r['matching_UMI_and_positive_cell_modes']}"
        )

    print()
    print("=" * 60)
    print("STRICT SHARED-TRUE CROSS-SAMPLE SUMMARY")
    print("=" * 60)

    for r in cross_rows:
        print(
            f"{r['gene']:8s} "
            f"SfLB={r['SfLB1_shared_TRUE_UMI']:>6},"
            f"{r['SfLB2_shared_TRUE_UMI']:>6}  "
            f"Tx={r['Tx17_shared_TRUE_UMI']:>6},"
            f"{r['Tx5_shared_TRUE_UMI']:>6}  "
            f"log2FC/cell="
            f"{r['strict_shared_TRUE_Tx_vs_SfLB_log2FC_per_cell_descriptive']:+.3f}"
        )

    print()
    print("=" * 60)
    print("OUTPUT")
    print("=" * 60)
    print(OUT)
    print()
    print("PASS: CANONICAL STANDARD-MEX CALIBRATION COMPLETE")
    print("PASS: IL10 AND 14 PRIORITIZED LIGANDS RECONSTRUCTED AT MOLECULE/PROBE LEVEL")
    print(f"PASS: {VERSION} PRIORITIZED LIGAND MOLECULE RECONSTRUCTION COMPLETE")
    print("=" * 60)

    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR, ignore_errors=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
