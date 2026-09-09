#!/usr/bin/env python3
# ==============================================================================
# Human_MASH_v1.4.4A.3B.1P.1_MeasuredGene_ProbeStatusAudit.py
#
# Purpose
# -------
# Cross-sample probe-definition/status audit for the 60 unique Flex-measured genes
# used in v1.4.4A.3B.1.
#
# This step DOES NOT scan molecule-level arrays and DOES NOT make biological
# SfLB-vs-Tx claims. It asks whether each measured gene can be compared across:
#   SfLB1/SfLB2 : Chromium Human Transcriptome Probe Set v1.0.1
#   Tx17/Tx5    : Chromium Human Transcriptome Probe Set v2.0.0
#
# For each HDF5 /probes table it records:
#   - probe index
#   - gene
#   - probe ID
#   - included_probes literal flag
#   - filtered_probes literal flag
#   - included AND filtered TRUE flag
#
# Cross-sample strict comparability is defined as:
#   >=1 identical short probe ID present in ALL 4 samples AND
#   included_probes == TRUE AND filtered_probes == TRUE in ALL 4 samples.
#
# IMPORTANT
# ---------
# No semantic assumption is made about the biological meaning of
# /probes/filtered_probes. TRUE/FALSE is reported literally. The strict shared
# subset is used only as a technical sensitivity/comparability criterion.
#
# No h5py dependency is required. Uses command-line h5dump.
# ==============================================================================

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

VERSION = "v1.4.4A.3B.1P.1"

BASE = Path(
    "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
)

PREFLIGHT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.3B.0"
    / "analysis"
    / "v1.4.4A.3B.0_graft_stress_xenogeneic_pressure_feature_preflight"
    / "tables"
    / "01_candidate_gene_feature_coverage_long_v1.4.4A.3B.0.csv"
)

OUT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.3B.1P.1"
    / "analysis"
    / "v1.4.4A.3B.1P.1_measured_gene_probe_status_audit"
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

ANTIGEN_GENES = {"TAP2", "TAPBP", "PSMB10", "NLRC5"}


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
        [
            h5dump,
            "-y",       # suppress array indices
            "-w", "0",  # no wrapping
            "-d", dataset,
            str(h5),
        ]
    )


def run_h5dump_header(h5dump: str, h5: Path) -> str:
    return run_cmd([h5dump, "-H", str(h5)])


def extract_data_block(text: str) -> str:
    m = re.search(
        r"\bDATA\s*\{\s*(.*?)\s*\}\s*\}",
        text,
        flags=re.S,
    )
    if not m:
        fail(
            "Could not locate DATA block in h5dump output.\n"
            + text[:1500]
        )
    return m.group(1)


def parse_string_dataset(text: str) -> list[str]:
    data = extract_data_block(text)

    vals = re.findall(
        r'"((?:\\.|[^"\\])*)"',
        data,
    )

    decoded = []
    for v in vals:
        try:
            decoded.append(
                bytes(v, "utf-8").decode(
                    "unicode_escape",
                    errors="replace",
                )
            )
        except Exception:
            decoded.append(v)

    return [normalize_h5_string(x) for x in decoded]


def normalize_h5_string(x: str) -> str:
    return x.replace("\x00", "").strip()


def short_probe_id(raw_probe_id: str) -> str:
    x = normalize_h5_string(raw_probe_id)
    if "|" in x:
        return x.split("|")[-1].strip()
    return x


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
            fail(
                f"Unexpected numeric enum value {v} "
                "in boolean dataset."
            )
        ans.append(bool(v))

    return ans


def bool_text(x: bool) -> str:
    return "TRUE" if x else "FALSE"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_preflight() -> tuple[list[str], dict[str, list[str]]]:
    if not PREFLIGHT.exists():
        fail(f"Preflight file missing: {PREFLIGHT}")

    with PREFLIGHT.open(
        newline="",
        encoding="utf-8",
    ) as fh:
        rows = list(csv.DictReader(fh))

    required = {"panel", "gene", "present_in_Flex_RNA"}
    if not rows:
        fail("Preflight file is empty.")
    if required - set(rows[0]):
        fail(
            "Preflight missing required fields: "
            + ", ".join(sorted(required - set(rows[0])))
        )

    def as_bool(x: str) -> bool:
        return str(x).strip().upper() in {"TRUE", "T", "1", "YES"}

    measured_rows = [
        r for r in rows
        if as_bool(r["present_in_Flex_RNA"])
    ]

    genes = sorted(
        {
            r["gene"].strip().upper()
            for r in measured_rows
            if r["gene"].strip()
        }
    )

    if len(genes) != 60:
        fail(
            f"Expected 60 unique measured genes from A.3B.0, "
            f"found {len(genes)}."
        )

    panel_genes: dict[str, list[str]] = defaultdict(list)

    for r in measured_rows:
        panel = r["panel"].strip()
        gene = r["gene"].strip().upper()

        if gene not in panel_genes[panel]:
            panel_genes[panel].append(gene)

    return genes, dict(panel_genes)


def detect_probe_set_version(
    h5dump: str,
    sample: str,
    h5: Path,
) -> str:
    header = run_h5dump_header(h5dump, h5)

    hits = sorted(
        set(
            re.findall(
                r"Chromium\\? Human\\? Transcriptome\\? Probe\\? Set\\? v[0-9.]+",
                header,
            )
        )
    )

    # h5dump may escape spaces as "\ " in names. Normalize both variants.
    normalized = [
        x.replace("\\ ", " ")
        for x in hits
    ]

    if not normalized:
        # More permissive fallback on surrounding target-set text.
        m = re.search(
            r"Chromium(?:\\ | )Human(?:\\ | )Transcriptome(?:\\ | )"
            r"Probe(?:\\ | )Set(?:\\ | )v[0-9.]+",
            header,
        )
        if m:
            normalized = [m.group(0).replace("\\ ", " ")]

    if len(normalized) != 1:
        fail(
            f"{sample}: could not uniquely determine probe-set version. "
            f"Observed={normalized}"
        )

    observed = normalized[0]
    expected = EXPECTED_PROBE_SET_VERSION[sample]

    if observed != expected:
        fail(
            f"{sample}: probe-set version mismatch. "
            f"Observed={observed!r}; expected={expected!r}"
        )

    return observed


def load_probe_table(
    h5dump: str,
    sample: str,
    h5: Path,
    measured_genes: set[str],
) -> list[dict]:
    print(f"Reading /probes arrays for {sample}...")

    names = parse_string_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/feature_name",
        )
    )

    feature_ids = parse_string_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/feature_id",
        )
    )

    probe_ids = parse_string_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/probe_id",
        )
    )

    included = parse_bool_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/included_probes",
        )
    )

    filtered = parse_bool_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/filtered_probes",
        )
    )

    regions = parse_string_dataset(
        run_h5dump(
            h5dump,
            h5,
            "/probes/region",
        )
    )

    lengths = {
        "feature_name": len(names),
        "feature_id": len(feature_ids),
        "probe_id": len(probe_ids),
        "included_probes": len(included),
        "filtered_probes": len(filtered),
        "region": len(regions),
    }

    if len(set(lengths.values())) != 1:
        fail(
            f"{sample}: inconsistent /probes array lengths: "
            + repr(lengths)
        )

    rows = []

    for i, (
        gene,
        feature_id,
        raw_probe_id,
        inc,
        filt,
        region,
    ) in enumerate(
        zip(
            names,
            feature_ids,
            probe_ids,
            included,
            filtered,
            regions,
        )
    ):
        gene_u = normalize_h5_string(gene).upper()

        if gene_u not in measured_genes:
            continue

        pid = short_probe_id(raw_probe_id)

        rows.append(
            {
                "sample": sample,
                "gene": gene_u,
                "probe_index_0based": i,
                "feature_id": normalize_h5_string(feature_id),
                "probe_id": pid,
                "probe_id_raw": normalize_h5_string(raw_probe_id),
                "included_probes": bool(inc),
                "filtered_probes": bool(filt),
                "included_and_filtered_TRUE": bool(inc and filt),
                "region": normalize_h5_string(region),
            }
        )

    return rows


def main() -> None:
    print("=" * 68)
    print(
        f"Human MASH {VERSION} measured-gene probe-status audit"
    )
    print("=" * 68)
    print()
    print("This is a probe-definition/status audit only.")
    print("No molecule-array scanning and no biological SfLB-vs-Tx inference.")
    print()

    OUT.mkdir(parents=True, exist_ok=True)

    measured_genes_list, panel_genes = read_preflight()
    measured_genes = set(measured_genes_list)

    print(f"Measured unique genes from A.3B.0: {len(measured_genes_list)}")
    print()

    h5dump = find_h5dump()
    print("h5dump:", h5dump)
    print()

    for sample, path in H5_FILES.items():
        if not path.exists():
            fail(f"Required H5 missing: {path}")

    # ------------------------------------------------------------------
    # Probe-set version validation
    # ------------------------------------------------------------------

    version_rows = []

    for sample in SAMPLES:
        observed = detect_probe_set_version(
            h5dump,
            sample,
            H5_FILES[sample],
        )

        version_rows.append(
            {
                "sample": sample,
                "probe_set_version": observed,
                "molecule_info_h5": str(H5_FILES[sample]),
            }
        )

        print(f"PASS: {sample}: {observed}")

    write_tsv(
        OUT / f"{VERSION}_probe_set_version_by_sample.tsv",
        version_rows,
        [
            "sample",
            "probe_set_version",
            "molecule_info_h5",
        ],
    )

    print()

    # ------------------------------------------------------------------
    # Read /probes metadata
    # ------------------------------------------------------------------

    all_probe_rows: list[dict] = []

    for sample in SAMPLES:
        rows = load_probe_table(
            h5dump,
            sample,
            H5_FILES[sample],
            measured_genes,
        )

        all_probe_rows.extend(rows)

        observed_genes = {r["gene"] for r in rows}
        missing = sorted(measured_genes - observed_genes)

        if missing:
            fail(
                f"{sample}: measured genes missing from /probes: "
                + ", ".join(missing)
            )

        gene_counts = defaultdict(int)
        for r in rows:
            gene_counts[r["gene"]] += 1

        print(
            f"PASS: {sample}: "
            f"{len(rows)} probe rows for all 60 measured genes"
        )
        print(
            "  probe rows/gene range: "
            f"{min(gene_counts.values())}–{max(gene_counts.values())}"
        )

    print()

    # ------------------------------------------------------------------
    # Long table
    # ------------------------------------------------------------------

    long_rows = []

    for r in all_probe_rows:
        long_rows.append(
            {
                **r,
                "included_probes": bool_text(r["included_probes"]),
                "filtered_probes": bool_text(r["filtered_probes"]),
                "included_and_filtered_TRUE": bool_text(
                    r["included_and_filtered_TRUE"]
                ),
            }
        )

    long_rows.sort(
        key=lambda r: (
            r["gene"],
            SAMPLES.index(r["sample"]),
            int(r["probe_index_0based"]),
        )
    )

    write_tsv(
        OUT / f"{VERSION}_measured_gene_probe_status_long.tsv",
        long_rows,
        [
            "sample",
            "gene",
            "probe_index_0based",
            "feature_id",
            "probe_id",
            "probe_id_raw",
            "included_probes",
            "filtered_probes",
            "included_and_filtered_TRUE",
            "region",
        ],
    )

    # ------------------------------------------------------------------
    # Per-sample per-gene summary
    # ------------------------------------------------------------------

    by_sample_gene: dict[tuple[str, str], list[dict]] = defaultdict(list)

    for r in all_probe_rows:
        by_sample_gene[(r["sample"], r["gene"])].append(r)

    sample_gene_rows = []

    for gene in measured_genes_list:
        for sample in SAMPLES:
            rows = by_sample_gene[(sample, gene)]

            included_rows = [
                r for r in rows
                if r["included_probes"]
            ]

            literal_filtered_true_rows = [
                r for r in rows
                if r["filtered_probes"]
            ]

            inc_filt_true_rows = [
                r for r in rows
                if r["included_and_filtered_TRUE"]
            ]

            included_filtered_false_rows = [
                r for r in rows
                if r["included_probes"] and not r["filtered_probes"]
            ]

            sample_gene_rows.append(
                {
                    "gene": gene,
                    "sample": sample,
                    "n_probe_rows": len(rows),
                    "n_included_probes": len(included_rows),
                    "n_literal_filtered_TRUE_probes": len(
                        literal_filtered_true_rows
                    ),
                    "n_included_and_filtered_TRUE_probes": len(
                        inc_filt_true_rows
                    ),
                    "n_included_filtered_FALSE_probes": len(
                        included_filtered_false_rows
                    ),
                    "mixed_included_filtered_status": bool_text(
                        bool(inc_filt_true_rows)
                        and bool(included_filtered_false_rows)
                    ),
                    "all_probe_ids": ";".join(
                        sorted({r["probe_id"] for r in rows})
                    ),
                    "included_probe_ids": ";".join(
                        sorted({r["probe_id"] for r in included_rows})
                    ),
                    "literal_filtered_TRUE_probe_ids": ";".join(
                        sorted(
                            {
                                r["probe_id"]
                                for r in literal_filtered_true_rows
                            }
                        )
                    ),
                    "included_and_filtered_TRUE_probe_ids": ";".join(
                        sorted(
                            {
                                r["probe_id"]
                                for r in inc_filt_true_rows
                            }
                        )
                    ),
                }
            )

    write_tsv(
        OUT / f"{VERSION}_measured_gene_probe_status_by_sample.tsv",
        sample_gene_rows,
        [
            "gene",
            "sample",
            "n_probe_rows",
            "n_included_probes",
            "n_literal_filtered_TRUE_probes",
            "n_included_and_filtered_TRUE_probes",
            "n_included_filtered_FALSE_probes",
            "mixed_included_filtered_status",
            "all_probe_ids",
            "included_probe_ids",
            "literal_filtered_TRUE_probe_ids",
            "included_and_filtered_TRUE_probe_ids",
        ],
    )

    # ------------------------------------------------------------------
    # Cross-sample intersection
    # ------------------------------------------------------------------

    cross_rows = []

    for gene in measured_genes_list:
        rows_by_sample = {
            s: by_sample_gene[(s, gene)]
            for s in SAMPLES
        }

        all_sets = {
            s: {r["probe_id"] for r in rows_by_sample[s]}
            for s in SAMPLES
        }

        included_sets = {
            s: {
                r["probe_id"]
                for r in rows_by_sample[s]
                if r["included_probes"]
            }
            for s in SAMPLES
        }

        inc_filt_sets = {
            s: {
                r["probe_id"]
                for r in rows_by_sample[s]
                if r["included_and_filtered_TRUE"]
            }
            for s in SAMPLES
        }

        same_all4 = set.intersection(
            *(all_sets[s] for s in SAMPLES)
        )

        shared_included_all4 = set.intersection(
            *(included_sets[s] for s in SAMPLES)
        )

        shared_strict_all4 = set.intersection(
            *(inc_filt_sets[s] for s in SAMPLES)
        )

        all_probe_sets_identical = all(
            all_sets[s] == all_sets[SAMPLES[0]]
            for s in SAMPLES[1:]
        )

        included_probe_sets_identical = all(
            included_sets[s] == included_sets[SAMPLES[0]]
            for s in SAMPLES[1:]
        )

        strict_sets_identical = all(
            inc_filt_sets[s] == inc_filt_sets[SAMPLES[0]]
            for s in SAMPLES[1:]
        )

        if shared_strict_all4:
            technical_class = "STRICT_SHARED_PROBE_COMPARABLE"
        elif same_all4:
            technical_class = (
                "SHARED_PROBE_IDS_BUT_NO_SHARED_INCLUDED_FILTERED_TRUE_ALL4"
            )
        else:
            technical_class = "NO_SHARED_PROBE_ID_ALL4"

        cross_rows.append(
            {
                "gene": gene,
                "technical_class": technical_class,
                "n_probe_ids_SfLB1": len(all_sets["SfLB1"]),
                "n_probe_ids_SfLB2": len(all_sets["SfLB2"]),
                "n_probe_ids_Tx17": len(all_sets["Tx17"]),
                "n_probe_ids_Tx5": len(all_sets["Tx5"]),
                "n_same_probe_ids_all4": len(same_all4),
                "same_probe_ids_all4": ";".join(sorted(same_all4)),
                "n_shared_included_probe_ids_all4": len(
                    shared_included_all4
                ),
                "shared_included_probe_ids_all4": ";".join(
                    sorted(shared_included_all4)
                ),
                "n_shared_included_filtered_TRUE_probe_ids_all4": len(
                    shared_strict_all4
                ),
                "shared_included_filtered_TRUE_probe_ids_all4": ";".join(
                    sorted(shared_strict_all4)
                ),
                "all_probe_id_sets_identical_all4": bool_text(
                    all_probe_sets_identical
                ),
                "included_probe_id_sets_identical_all4": bool_text(
                    included_probe_sets_identical
                ),
                "included_filtered_TRUE_sets_identical_all4": bool_text(
                    strict_sets_identical
                ),
            }
        )

    cross_rows.sort(
        key=lambda r: (
            0 if r["technical_class"] == "STRICT_SHARED_PROBE_COMPARABLE" else 1,
            r["gene"],
        )
    )

    write_tsv(
        OUT / f"{VERSION}_measured_gene_crosssample_probe_summary.tsv",
        cross_rows,
        [
            "gene",
            "technical_class",
            "n_probe_ids_SfLB1",
            "n_probe_ids_SfLB2",
            "n_probe_ids_Tx17",
            "n_probe_ids_Tx5",
            "n_same_probe_ids_all4",
            "same_probe_ids_all4",
            "n_shared_included_probe_ids_all4",
            "shared_included_probe_ids_all4",
            "n_shared_included_filtered_TRUE_probe_ids_all4",
            "shared_included_filtered_TRUE_probe_ids_all4",
            "all_probe_id_sets_identical_all4",
            "included_probe_id_sets_identical_all4",
            "included_filtered_TRUE_sets_identical_all4",
        ],
    )

    # ------------------------------------------------------------------
    # Program-level strict-comparability coverage
    # ------------------------------------------------------------------

    cross_by_gene = {
        r["gene"]: r
        for r in cross_rows
    }

    program_rows = []

    for panel in sorted(panel_genes):
        genes = sorted(set(panel_genes[panel]))

        strict = [
            g for g in genes
            if cross_by_gene[g]["technical_class"]
            == "STRICT_SHARED_PROBE_COMPARABLE"
        ]

        shared_but_not_strict = [
            g for g in genes
            if cross_by_gene[g]["technical_class"]
            == "SHARED_PROBE_IDS_BUT_NO_SHARED_INCLUDED_FILTERED_TRUE_ALL4"
        ]

        no_shared = [
            g for g in genes
            if cross_by_gene[g]["technical_class"]
            == "NO_SHARED_PROBE_ID_ALL4"
        ]

        program_rows.append(
            {
                "panel": panel,
                "n_measured_genes": len(genes),
                "n_strict_shared_probe_comparable": len(strict),
                "fraction_strict_shared_probe_comparable": (
                    len(strict) / len(genes) if genes else 0.0
                ),
                "strict_shared_probe_comparable_genes": ";".join(strict),
                "n_shared_probe_but_not_strict": len(shared_but_not_strict),
                "shared_probe_but_not_strict_genes": ";".join(
                    shared_but_not_strict
                ),
                "n_no_shared_probe_id_all4": len(no_shared),
                "no_shared_probe_id_all4_genes": ";".join(no_shared),
            }
        )

    write_tsv(
        OUT / f"{VERSION}_program_crosssample_probe_coverage.tsv",
        program_rows,
        [
            "panel",
            "n_measured_genes",
            "n_strict_shared_probe_comparable",
            "fraction_strict_shared_probe_comparable",
            "strict_shared_probe_comparable_genes",
            "n_shared_probe_but_not_strict",
            "shared_probe_but_not_strict_genes",
            "n_no_shared_probe_id_all4",
            "no_shared_probe_id_all4_genes",
        ],
    )

    # ------------------------------------------------------------------
    # Dedicated antigen-presentation 4-gene table
    # ------------------------------------------------------------------

    antigen_rows = [
        r for r in cross_rows
        if r["gene"] in ANTIGEN_GENES
    ]

    antigen_rows.sort(key=lambda r: r["gene"])

    write_tsv(
        OUT / f"{VERSION}_antigen4_crosssample_probe_detail.tsv",
        antigen_rows,
        [
            "gene",
            "technical_class",
            "n_probe_ids_SfLB1",
            "n_probe_ids_SfLB2",
            "n_probe_ids_Tx17",
            "n_probe_ids_Tx5",
            "n_same_probe_ids_all4",
            "same_probe_ids_all4",
            "n_shared_included_probe_ids_all4",
            "shared_included_probe_ids_all4",
            "n_shared_included_filtered_TRUE_probe_ids_all4",
            "shared_included_filtered_TRUE_probe_ids_all4",
            "all_probe_id_sets_identical_all4",
            "included_probe_id_sets_identical_all4",
            "included_filtered_TRUE_sets_identical_all4",
        ],
    )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------

    n_strict = sum(
        r["technical_class"] == "STRICT_SHARED_PROBE_COMPARABLE"
        for r in cross_rows
    )
    n_shared_not_strict = sum(
        r["technical_class"]
        == "SHARED_PROBE_IDS_BUT_NO_SHARED_INCLUDED_FILTERED_TRUE_ALL4"
        for r in cross_rows
    )
    n_no_shared = sum(
        r["technical_class"] == "NO_SHARED_PROBE_ID_ALL4"
        for r in cross_rows
    )

    print()
    print("=" * 68)
    print("CROSS-SAMPLE PROBE COMPARABILITY SUMMARY")
    print("=" * 68)
    print(f"Measured unique genes:                  {len(cross_rows)}")
    print(f"STRICT shared-probe comparable:         {n_strict}")
    print(f"Shared IDs but no strict shared subset: {n_shared_not_strict}")
    print(f"No shared probe ID across all 4:        {n_no_shared}")
    print()

    print("=== PROGRAM COVERAGE ===")
    for r in program_rows:
        print(
            f"{r['panel']}: "
            f"{r['n_strict_shared_probe_comparable']}/"
            f"{r['n_measured_genes']} strict-comparable "
            f"({r['fraction_strict_shared_probe_comparable']:.3f})"
        )

    print()
    print("=== ANTIGEN-PRESENTATION 4-GENE DETAIL ===")
    for r in antigen_rows:
        print(
            f"{r['gene']}: {r['technical_class']}; "
            f"shared IDs={r['n_same_probe_ids_all4']}; "
            f"strict shared={r['n_shared_included_filtered_TRUE_probe_ids_all4']}; "
            f"all-probe-set-identical={r['all_probe_id_sets_identical_all4']}"
        )

    print()
    print("=" * 68)
    print(f"{VERSION} COMPLETE")
    print("=" * 68)
    print("Output:", OUT)
    print("No molecule-level reconstruction performed.")
    print("No biological SfLB-vs-Tx interpretation performed.")
    print(
        "Next step: use strict shared probe IDs for molecule-level "
        "reconstruction of technically comparable genes."
    )


if __name__ == "__main__":
    main()
