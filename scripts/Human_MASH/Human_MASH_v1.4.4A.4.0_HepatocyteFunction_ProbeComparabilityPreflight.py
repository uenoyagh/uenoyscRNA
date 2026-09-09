#!/usr/bin/env python3
# ==============================================================================
# Human_MASH_v1.4.4A.4.0_HepatocyteFunction_ProbeComparabilityPreflight.py
#
# Purpose
# -------
# Probe-comparability preflight for the Flex-measurable hepatocyte function /
# maturity panel before any SfLB-vs-Tx scoring.
#
# Candidate measured genes (16):
#   Secretory/synthetic: TTR, APOA1, FGB, FGG
#   Urea/nitrogen:       OTC, ASS1, ASL, NAGS
#   Bile/transport:      ABCB11
#   Lipid/metabolic:     MTTP, FABP1
#   Hepatocyte identity: HNF4A, CEBPA
#   Epithelial/immature: EPCAM, KRT7, SOX9
#
# Comparison:
#   SfLB1/SfLB2 : Chromium Human Transcriptome Probe Set v1.0.1
#   Tx17/Tx5    : Chromium Human Transcriptome Probe Set v2.0.0
#
# This step DOES NOT scan molecule arrays and DOES NOT score function.
# It reports whether identical probe IDs exist and whether an
# included_probes=TRUE & filtered_probes=TRUE subset is shared across all 4.
#
# No semantic assumption is made about /probes/filtered_probes.
# ==============================================================================

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

VERSION = "v1.4.4A.4.0"

BASE = Path("/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk")

OUT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.4.0"
    / "analysis"
    / "v1.4.4A.4.0_hepatocyte_function_probe_comparability_preflight"
)

SAMPLES = ["SfLB1", "SfLB2", "Tx17", "Tx5"]

H5_FILES = {
    "SfLB1": Path(
        "/Volumes/RH241011956/Deliver_RH241011956/CellRanger/"
        "okmt_pool1_outs/per_sample_outs/SfLB_1/count/"
        "sample_molecule_info.h5"
    ),
    "SfLB2": Path(
        "/Volumes/RH241011956/Deliver_RH241011956/CellRanger/"
        "okmt_pool1_outs/per_sample_outs/SfLB_2/count/"
        "sample_molecule_info.h5"
    ),
    "Tx17": Path(
        "/Volumes/RH251217117/Deliver_RH251217117/Result_human/"
        "01_CellRanger/Count_matrix/Tx17/molecule_info.h5"
    ),
    "Tx5": Path(
        "/Volumes/RH251217117/Deliver_RH251217117/Result_human/"
        "01_CellRanger/Count_matrix/Tx5/molecule_info.h5"
    ),
}

EXPECTED_PROBE_SET_VERSION = {
    "SfLB1": "Chromium Human Transcriptome Probe Set v1.0.1",
    "SfLB2": "Chromium Human Transcriptome Probe Set v1.0.1",
    "Tx17": "Chromium Human Transcriptome Probe Set v2.0.0",
    "Tx5": "Chromium Human Transcriptome Probe Set v2.0.0",
}

PANELS = {
    "Secretory_synthetic": ["TTR", "APOA1", "FGB", "FGG"],
    "Urea_nitrogen": ["OTC", "ASS1", "ASL", "NAGS"],
    "Bile_transport": ["ABCB11"],
    "Lipid_metabolic": ["MTTP", "FABP1"],
    "Hepatocyte_identity": ["HNF4A", "CEBPA"],
    "Epithelial_immaturity_transition": ["EPCAM", "KRT7", "SOX9"],
}

CANDIDATES = sorted({g for genes in PANELS.values() for g in genes})


def fail(msg: str) -> None:
    raise RuntimeError(msg)


def find_h5dump() -> str:
    preferred = Path("/opt/homebrew/bin/h5dump")
    if preferred.exists():
        return str(preferred)
    x = shutil.which("h5dump")
    if x:
        return x
    fail("h5dump not found.")


def run_cmd(cmd: list[str]) -> str:
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if p.returncode != 0:
        fail(
            "Command failed:\n"
            + " ".join(cmd)
            + "\nSTDERR:\n"
            + p.stderr
        )
    return p.stdout


def run_h5dump(h5dump: str, h5: Path, dataset: str) -> str:
    return run_cmd(
        [h5dump, "-y", "-w", "0", "-d", dataset, str(h5)]
    )


def extract_data_block(text: str) -> str:
    m = re.search(
        r"\bDATA\s*\{\s*(.*?)\s*\}\s*\}",
        text,
        flags=re.S,
    )
    if not m:
        fail("Could not locate DATA block in h5dump output.")
    return m.group(1)


def normalize_h5_string(x: str) -> str:
    return x.replace("\x00", "").strip()


def parse_string_dataset(text: str) -> list[str]:
    data = extract_data_block(text)
    vals = re.findall(r'"((?:\\.|[^"\\])*)"', data)

    out = []
    for v in vals:
        try:
            x = bytes(v, "utf-8").decode(
                "unicode_escape",
                errors="replace",
            )
        except Exception:
            x = v
        out.append(normalize_h5_string(x))
    return out


def parse_bool_dataset(text: str) -> list[bool]:
    data = extract_data_block(text)

    tf = re.findall(
        r"\b(TRUE|FALSE)\b",
        data,
        flags=re.I,
    )
    if tf:
        return [x.upper() == "TRUE" for x in tf]

    nums = re.findall(r"(?<![\w.])-?\d+(?![\w.])", data)

    ans = []
    for x in nums:
        v = int(x)
        if v not in (0, 1):
            fail(f"Unexpected boolean enum value: {v}")
        ans.append(bool(v))
    return ans


def short_probe_id(raw_probe_id: str) -> str:
    x = normalize_h5_string(raw_probe_id)
    if "|" in x:
        return x.split("|")[-1].strip()
    return x


def bool_text(x: bool) -> str:
    return "TRUE" if x else "FALSE"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def detect_probe_set_version(h5dump: str, sample: str, h5: Path) -> str:
    header = run_cmd([h5dump, "-H", str(h5)])

    m = re.search(
        r"Chromium(?:\\ | )Human(?:\\ | )Transcriptome(?:\\ | )"
        r"Probe(?:\\ | )Set(?:\\ | )v[0-9.]+",
        header,
    )
    if not m:
        fail(f"{sample}: probe-set version not found.")

    observed = m.group(0).replace("\\ ", " ")
    expected = EXPECTED_PROBE_SET_VERSION[sample]

    if observed != expected:
        fail(
            f"{sample}: probe-set mismatch: "
            f"observed={observed!r}, expected={expected!r}"
        )

    return observed


def load_probe_rows(h5dump: str, sample: str, h5: Path) -> list[dict]:
    names = parse_string_dataset(
        run_h5dump(h5dump, h5, "/probes/feature_name")
    )
    ids = parse_string_dataset(
        run_h5dump(h5dump, h5, "/probes/probe_id")
    )
    included = parse_bool_dataset(
        run_h5dump(h5dump, h5, "/probes/included_probes")
    )
    filtered = parse_bool_dataset(
        run_h5dump(h5dump, h5, "/probes/filtered_probes")
    )
    regions = parse_string_dataset(
        run_h5dump(h5dump, h5, "/probes/region")
    )

    lens = {
        len(names),
        len(ids),
        len(included),
        len(filtered),
        len(regions),
    }
    if len(lens) != 1:
        fail(f"{sample}: inconsistent /probes array lengths.")

    rows = []

    for i, (gene, raw_pid, inc, filt, region) in enumerate(
        zip(names, ids, included, filtered, regions)
    ):
        g = normalize_h5_string(gene).upper()
        if g not in CANDIDATES:
            continue

        rows.append(
            {
                "sample": sample,
                "gene": g,
                "probe_index_0based": i,
                "probe_id": short_probe_id(raw_pid),
                "probe_id_raw": normalize_h5_string(raw_pid),
                "included_probes": bool(inc),
                "filtered_probes": bool(filt),
                "included_and_filtered_TRUE": bool(inc and filt),
                "region": normalize_h5_string(region),
            }
        )

    return rows


def main() -> None:
    print("=" * 70)
    print(
        f"Human MASH {VERSION} hepatocyte-function probe comparability preflight"
    )
    print("=" * 70)
    print()
    print(f"Candidate measured genes: {len(CANDIDATES)}")
    print(", ".join(CANDIDATES))
    print()

    OUT.mkdir(parents=True, exist_ok=True)

    h5dump = find_h5dump()
    print("h5dump:", h5dump)
    print()

    for p in H5_FILES.values():
        if not p.exists():
            fail(f"Required input missing: {p}")

    all_rows = []

    for sample in SAMPLES:
        ver = detect_probe_set_version(
            h5dump,
            sample,
            H5_FILES[sample],
        )
        print(f"PASS: {sample}: {ver}")

        rows = load_probe_rows(
            h5dump,
            sample,
            H5_FILES[sample],
        )

        genes = {r["gene"] for r in rows}
        missing = sorted(set(CANDIDATES) - genes)

        if missing:
            fail(
                f"{sample}: candidate genes missing from /probes: "
                + ", ".join(missing)
            )

        all_rows.extend(rows)

        counts = defaultdict(int)
        for r in rows:
            counts[r["gene"]] += 1

        print(
            f"  {len(rows)} probe rows; "
            f"probe rows/gene range "
            f"{min(counts.values())}-{max(counts.values())}"
        )

    long_rows = []
    for r in all_rows:
        long_rows.append(
            {
                **r,
                "included_probes": bool_text(r["included_probes"]),
                "filtered_probes": bool_text(r["filtered_probes"]),
                "included_and_filtered_TRUE":
                    bool_text(r["included_and_filtered_TRUE"]),
            }
        )

    write_tsv(
        OUT / f"{VERSION}_hepatocyte_function_probe_status_long.tsv",
        long_rows,
        [
            "sample",
            "gene",
            "probe_index_0based",
            "probe_id",
            "probe_id_raw",
            "included_probes",
            "filtered_probes",
            "included_and_filtered_TRUE",
            "region",
        ],
    )

    by_sg = defaultdict(list)
    for r in all_rows:
        by_sg[(r["sample"], r["gene"])].append(r)

    cross_rows = []

    for gene in CANDIDATES:
        all_sets = {
            s: {r["probe_id"] for r in by_sg[(s, gene)]}
            for s in SAMPLES
        }

        strict_sets = {
            s: {
                r["probe_id"]
                for r in by_sg[(s, gene)]
                if r["included_and_filtered_TRUE"]
            }
            for s in SAMPLES
        }

        shared_all = set.intersection(
            *(all_sets[s] for s in SAMPLES)
        )
        shared_strict = set.intersection(
            *(strict_sets[s] for s in SAMPLES)
        )

        if shared_strict:
            technical_class = "STRICT_SHARED_PROBE_COMPARABLE"
        elif shared_all:
            technical_class = (
                "SHARED_PROBE_IDS_BUT_NO_SHARED_INCLUDED_FILTERED_TRUE_ALL4"
            )
        else:
            technical_class = "NO_SHARED_PROBE_ID_ALL4"

        panels = [
            p for p, genes in PANELS.items()
            if gene in genes
        ]

        cross_rows.append(
            {
                "gene": gene,
                "panel": ";".join(panels),
                "technical_class": technical_class,
                "n_probe_ids_SfLB1": len(all_sets["SfLB1"]),
                "n_probe_ids_SfLB2": len(all_sets["SfLB2"]),
                "n_probe_ids_Tx17": len(all_sets["Tx17"]),
                "n_probe_ids_Tx5": len(all_sets["Tx5"]),
                "n_shared_probe_ids_all4": len(shared_all),
                "shared_probe_ids_all4": ";".join(sorted(shared_all)),
                "n_strict_shared_probe_ids_all4": len(shared_strict),
                "strict_shared_probe_ids_all4": ";".join(
                    sorted(shared_strict)
                ),
                "all_probe_sets_identical_all4": bool_text(
                    all(
                        all_sets[s] == all_sets[SAMPLES[0]]
                        for s in SAMPLES[1:]
                    )
                ),
            }
        )

    write_tsv(
        OUT / f"{VERSION}_hepatocyte_function_crosssample_probe_summary.tsv",
        cross_rows,
        [
            "gene",
            "panel",
            "technical_class",
            "n_probe_ids_SfLB1",
            "n_probe_ids_SfLB2",
            "n_probe_ids_Tx17",
            "n_probe_ids_Tx5",
            "n_shared_probe_ids_all4",
            "shared_probe_ids_all4",
            "n_strict_shared_probe_ids_all4",
            "strict_shared_probe_ids_all4",
            "all_probe_sets_identical_all4",
        ],
    )

    panel_rows = []

    by_gene = {r["gene"]: r for r in cross_rows}

    for panel, genes in PANELS.items():
        strict = [
            g for g in genes
            if by_gene[g]["technical_class"]
            == "STRICT_SHARED_PROBE_COMPARABLE"
        ]
        nonstrict = [
            g for g in genes
            if by_gene[g]["technical_class"]
            != "STRICT_SHARED_PROBE_COMPARABLE"
        ]

        panel_rows.append(
            {
                "panel": panel,
                "n_measured_genes": len(genes),
                "n_strict_comparable": len(strict),
                "fraction_strict_comparable":
                    len(strict) / len(genes),
                "strict_comparable_genes": ";".join(strict),
                "non_strict_genes": ";".join(nonstrict),
            }
        )

    write_tsv(
        OUT / f"{VERSION}_hepatocyte_function_panel_probe_coverage.tsv",
        panel_rows,
        [
            "panel",
            "n_measured_genes",
            "n_strict_comparable",
            "fraction_strict_comparable",
            "strict_comparable_genes",
            "non_strict_genes",
        ],
    )

    n_strict = sum(
        r["technical_class"] == "STRICT_SHARED_PROBE_COMPARABLE"
        for r in cross_rows
    )

    print()
    print("=" * 70)
    print("CROSS-SAMPLE PROBE COMPARABILITY")
    print("=" * 70)
    print(f"Strict-comparable genes: {n_strict}/{len(CANDIDATES)}")
    print()

    for r in panel_rows:
        print(
            f"{r['panel']}: "
            f"{r['n_strict_comparable']}/"
            f"{r['n_measured_genes']} strict-comparable"
        )

    print()
    print("=== NON-STRICT GENES ===")
    bad = [
        r for r in cross_rows
        if r["technical_class"] != "STRICT_SHARED_PROBE_COMPARABLE"
    ]
    if not bad:
        print("None")
    else:
        for r in bad:
            print(
                f"{r['gene']}: {r['technical_class']}; "
                f"shared IDs={r['n_shared_probe_ids_all4']}; "
                f"strict shared={r['n_strict_shared_probe_ids_all4']}"
            )

    print()
    print("=" * 70)
    print(f"{VERSION} COMPLETE")
    print("=" * 70)
    print("Output:", OUT)
    print("No molecule-level reconstruction performed.")
    print("No hepatocyte function score performed.")
    print(
        "Next step: reconstruct strict shared-probe molecules only for "
        "technically comparable hepatic genes."
    )


if __name__ == "__main__":
    main()
