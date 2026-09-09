#!/usr/bin/env python3
# ==============================================================================
# Human_MASH_v1.4.4A.3B.1P.2_StrictSharedProbe_MoleculeReconstruction.py
#
# Purpose
# -------
# Molecule-level reconstruction for the 49 genes classified in v1.4.4A.3B.1P.1
# as STRICT_SHARED_PROBE_COMPARABLE across:
#
#   SfLB1 / SfLB2 : Chromium Human Transcriptome Probe Set v1.0.1
#   Tx17 / Tx5    : Chromium Human Transcriptome Probe Set v2.0.0
#
# Only probe IDs that satisfy ALL of the following are used:
#   1) same short probe ID exists in all 4 samples
#   2) included_probes == TRUE in all 4 samples
#   3) filtered_probes == TRUE in all 4 samples
#
# IMPORTANT TECHNICAL RULE
# ------------------------
# UMI / positive-cell numerators are reconstructed from molecule_info.h5.
# Denominators are ALWAYS the full frozen canonical cell universe from
# v1.4.3.1 context metadata. They are NOT based on cells with selected signal.
#
# The verified v1.4.3B.2C.3 low-level scanner / barcode-resolution functions
# are imported and reused rather than reimplemented.
#
# This step is gene-level reconstruction and QC. It does NOT yet rebuild the
# multi-gene program scores and does NOT make a final xenogeneic-rejection claim.
# ==============================================================================

from __future__ import annotations

import csv
import importlib.util
import math
import shutil
import struct
from collections import defaultdict
from pathlib import Path

VERSION = "v1.4.4A.3B.1P.2"

REPO = Path("/Users/uenoya/Projects/uenoyscRNA")

BASE_ENGINE = (
    REPO
    / "scripts"
    / "Human_MASH"
    / "Human_MASH_v1.4.3B.2C.3_PrioritizedLigand_MoleculeReconstruction.py"
)

BASE = Path("/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk")

P1_DIR = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.3B.1P.1"
    / "analysis"
    / "v1.4.4A.3B.1P.1_measured_gene_probe_status_audit"
)

P1_LONG = P1_DIR / "v1.4.4A.3B.1P.1_measured_gene_probe_status_long.tsv"
P1_CROSS = P1_DIR / "v1.4.4A.3B.1P.1_measured_gene_crosssample_probe_summary.tsv"
P1_PROGRAM = P1_DIR / "v1.4.4A.3B.1P.1_program_crosssample_probe_coverage.tsv"

CONTEXT_META = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.3.1"
    / "analysis"
    / "v1.4.3.1_12sample_integrated_context_review"
    / "tables"
    / "v1.4.3_context_UMAP_metadata.tsv"
)

A11_COUNTS = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.1.1"
    / "analysis"
    / "v1.4.4A.1.1_pretransplant_SfLB_vs_engraftedTx_lineage_composition"
    / "tables"
    / "02_broad_lineage_counts_fractions_by_sample_v1.4.4A.1.1.csv"
)

A3B1_GENE_METRICS = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.3B.1"
    / "analysis"
    / "v1.4.4A.3B.1_graft_stress_immune_interaction_measured_program_audit"
    / "tables"
    / "06_measured_gene_metrics_by_sample_lineage_v1.4.4A.3B.1.csv"
)

OUT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.3B.1P.2"
    / "analysis"
    / "v1.4.4A.3B.1P.2_strict_shared_probe_molecule_reconstruction"
)

TABLE_DIR = OUT / "tables"
TMP_DIR = OUT / "_tmp"

SAMPLES = ["SfLB1", "SfLB2", "Tx17", "Tx5"]

EXPECTED_CANONICAL_CELLS = {
    "SfLB1": 19081,
    "SfLB2": 9805,
    "Tx17": 5958,
    "Tx5": 8571,
}

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

COMPARABLE_LINEAGES = {"Hepatocyte_like", "Mesenchymal_like"}
ANTIGEN_GENES = {"TAP2", "TAPBP", "PSMB10", "NLRC5"}
MIN_FREE_BYTES = 3 * 1024**3


def fail(msg: str) -> None:
    raise RuntimeError(msg)


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


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


def norm_sample(x: str) -> str:
    x = x.strip()
    if x == "SfLB_1":
        return "SfLB1"
    if x == "SfLB_2":
        return "SfLB2"
    return x


def bool_text(x) -> bool:
    return str(x).strip().upper() == "TRUE"


def import_verified_engine():
    if not BASE_ENGINE.exists():
        fail(f"Verified B.2C.3 engine missing: {BASE_ENGINE}")

    spec = importlib.util.spec_from_file_location(
        "v143b2c3_engine",
        BASE_ENGINE,
    )
    if spec is None or spec.loader is None:
        fail("Could not construct import spec for B.2C.3 engine.")

    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    required = [
        "find_executable",
        "compile_helper",
        "export_sample_arrays",
        "make_target_file",
        "scan_selected",
        "collect_selected_barcode_indices",
        "resolve_selected_barcodes_from_h5dump",
        "extract_barcode_core",
        "mol_key",
    ]

    missing = [x for x in required if not hasattr(mod, x)]
    if missing:
        fail(
            "Verified B.2C.3 engine missing required functions: "
            + ", ".join(missing)
        )

    # compile_helper uses module-level TMP_DIR.
    mod.TMP_DIR = TMP_DIR

    return mod


def load_strict_probe_definition():
    cross = read_tsv(P1_CROSS)
    long_rows = read_tsv(P1_LONG)

    if not cross or not long_rows:
        fail("P.1 input tables are empty.")

    strict_genes = sorted(
        r["gene"].strip().upper()
        for r in cross
        if r["technical_class"].strip()
        == "STRICT_SHARED_PROBE_COMPARABLE"
    )

    if len(strict_genes) != 49:
        fail(
            f"Expected 49 strict-comparable genes from P.1, "
            f"found {len(strict_genes)}."
        )

    cross_by_gene = {
        r["gene"].strip().upper(): r
        for r in cross
    }

    shared_ids = {}

    for g in strict_genes:
        raw = cross_by_gene[g][
            "shared_included_filtered_TRUE_probe_ids_all4"
        ].strip()
        ids = {x for x in raw.split(";") if x}

        if not ids:
            fail(f"{g}: strict class but no shared strict probe IDs.")

        shared_ids[g] = ids

    probe_meta = {s: {} for s in SAMPLES}
    selected_rows = []

    for r in long_rows:
        s = norm_sample(r["sample"])
        g = r["gene"].strip().upper()

        if s not in SAMPLES or g not in shared_ids:
            continue

        pid = r["probe_id"].strip()

        if pid not in shared_ids[g]:
            continue

        if not (
            bool_text(r["included_probes"])
            and bool_text(r["filtered_probes"])
        ):
            fail(
                f"{s} {g} {pid}: P.1 shared strict probe does not "
                "have included=TRUE and filtered=TRUE in long table."
            )

        idx = int(r["probe_index_0based"])

        if idx in probe_meta[s]:
            fail(f"{s}: duplicated selected probe index {idx}")

        meta = {
            "sample": s,
            "gene": g,
            "probe_idx": idx,
            "probe_id": pid,
            "included": True,
            "filtered": True,
        }

        probe_meta[s][idx] = meta
        selected_rows.append(
            {
                "sample": s,
                "gene": g,
                "probe_index_0based": idx,
                "probe_id": pid,
                "strict_shared_all4": "TRUE",
            }
        )

    # Every strict gene must have >=1 selected probe in every sample.
    for s in SAMPLES:
        seen = {x["gene"] for x in probe_meta[s].values()}
        missing = sorted(set(strict_genes) - seen)
        if missing:
            fail(
                f"{s}: strict genes missing selected probe rows: "
                + ", ".join(missing)
            )

    selected_rows.sort(
        key=lambda r: (
            r["gene"],
            SAMPLES.index(r["sample"]),
            int(r["probe_index_0based"]),
        )
    )

    return strict_genes, shared_ids, probe_meta, selected_rows


def load_context_full():
    rows = read_tsv(CONTEXT_META)

    required = {
        "sample",
        "barcode_core",
        "context_annotation",
        "context_state",
    }

    if not rows:
        fail("Context metadata is empty.")

    missing = required - set(rows[0])
    if missing:
        fail(
            "Context metadata missing: "
            + ", ".join(sorted(missing))
        )

    out = {s: {} for s in SAMPLES}

    for r in rows:
        s = norm_sample(r["sample"])
        if s not in SAMPLES:
            continue

        bc = r["barcode_core"].strip()
        if not bc:
            fail(f"{s}: empty barcode_core")

        if bc in out[s]:
            fail(f"{s}: duplicated canonical barcode {bc}")

        out[s][bc] = {
            "annotation": r["context_annotation"].strip(),
            "state": r["context_state"].strip(),
        }

    for s in SAMPLES:
        if len(out[s]) != EXPECTED_CANONICAL_CELLS[s]:
            fail(
                f"{s}: canonical cells={len(out[s])}; "
                f"expected={EXPECTED_CANONICAL_CELLS[s]}"
            )

    return out


def build_full_denominators(context):
    sample_den = {}
    lineage_den = {}
    state_den = {}

    for s in SAMPLES:
        sample_den[s] = len(context[s])

        lc = defaultdict(int)
        sc = defaultdict(int)

        for meta in context[s].values():
            lc[meta["annotation"]] += 1
            sc[meta["state"]] += 1

        lineage_den[s] = dict(lc)
        state_den[s] = dict(sc)

    return sample_den, lineage_den, state_den


def validate_full_denominators_against_a11(lineage_den):
    rows = read_csv(A11_COUNTS)

    required = {"sample", "category", "n_cells"}
    if not rows:
        fail("A.1.1 lineage table empty.")
    if required - set(rows[0]):
        fail(
            "A.1.1 lineage table missing: "
            + ", ".join(sorted(required - set(rows[0])))
        )

    expected = {}
    for r in rows:
        s = norm_sample(r["sample"])
        if s not in SAMPLES:
            continue
        expected[(s, r["category"].strip())] = int(float(r["n_cells"]))

    observed = {}
    for s in SAMPLES:
        for lineage, n in lineage_den[s].items():
            observed[(s, lineage)] = n

    all_keys = set(expected) | set(observed)

    mismatches = []
    for k in sorted(all_keys):
        a = expected.get(k, 0)
        b = observed.get(k, 0)
        if a != b:
            mismatches.append((k, a, b))

    if mismatches:
        txt = "; ".join(
            f"{s}/{ll}:A11={a},context={b}"
            for ((s, ll), a, b) in mismatches[:30]
        )
        fail(
            "Full context lineage denominator mismatch vs A.1.1: "
            + txt
        )

    return True


def load_standard_gene_metrics():
    rows = read_csv(A3B1_GENE_METRICS)

    required = {
        "sample",
        "lineage",
        "gene",
        "n_cells",
        "positive_cells",
        "total_UMI",
        "UMI_per_1000_cells",
    }

    if not rows:
        fail("A.3B.1 gene metrics table empty.")
    if required - set(rows[0]):
        fail(
            "A.3B.1 gene metrics missing: "
            + ", ".join(sorted(required - set(rows[0])))
        )

    out = {}
    for r in rows:
        s = norm_sample(r["sample"])
        if s not in SAMPLES:
            continue

        key = (
            s,
            r["lineage"].strip(),
            r["gene"].strip().upper(),
        )

        if key in out:
            fail(f"Duplicated A.3B.1 standard metric row: {key}")

        out[key] = {
            "standard_n_cells": int(float(r["n_cells"])),
            "standard_positive_cells": int(float(r["positive_cells"])),
            "standard_total_UMI": int(round(float(r["total_UMI"]))),
            "standard_UMI_per_1000_cells": float(
                r["UMI_per_1000_cells"]
            ),
        }

    return out


def relation_2v2(sf1, sf2, tx17, tx5):
    vals = [sf1, sf2, tx17, tx5]
    if any(v is None or not math.isfinite(v) for v in vals):
        return "Not_comparable_missing_value"

    s = [sf1, sf2]
    t = [tx17, tx5]

    if min(t) > max(s):
        return "Tx_higher_complete_replicate_range_separation"
    if max(t) < min(s):
        return "Tx_lower_complete_replicate_range_separation"
    return "Replicate_ranges_overlap"


def main():
    print("=" * 72)
    print(
        f"Human MASH {VERSION} strict shared-probe molecule reconstruction"
    )
    print("=" * 72)
    print()
    print("Verified B.2C.3 low-level scanner/barcode engine will be reused.")
    print("Full frozen canonical denominators will be used.")
    print()

    required_paths = [
        BASE_ENGINE,
        P1_LONG,
        P1_CROSS,
        P1_PROGRAM,
        CONTEXT_META,
        A11_COUNTS,
        A3B1_GENE_METRICS,
        *H5_FILES.values(),
    ]

    for p in required_paths:
        if not p.exists():
            fail(f"Required input missing: {p}")

    OUT.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    free_bytes = shutil.disk_usage(OUT).free
    if free_bytes < MIN_FREE_BYTES:
        fail(
            f"Insufficient free space at {OUT}: "
            f"{free_bytes / 1024**3:.2f} GiB; "
            f"require >= {MIN_FREE_BYTES / 1024**3:.2f} GiB."
        )

    strict_genes, shared_ids, probe_meta, selected_probe_rows = (
        load_strict_probe_definition()
    )

    print(f"Strict-comparable genes: {len(strict_genes)} / 60")
    print(
        "Strict shared probes/gene range: "
        f"{min(len(shared_ids[g]) for g in strict_genes)}–"
        f"{max(len(shared_ids[g]) for g in strict_genes)}"
    )
    print()

    write_tsv(
        TABLE_DIR / f"{VERSION}_strict_probe_definition_used.tsv",
        selected_probe_rows,
        [
            "sample",
            "gene",
            "probe_index_0based",
            "probe_id",
            "strict_shared_all4",
        ],
    )

    context = load_context_full()
    sample_den, lineage_den, state_den = build_full_denominators(context)

    validate_full_denominators_against_a11(lineage_den)

    print("Canonical context cell counts: PASS")
    print("Full lineage denominators vs A.1.1: PASS")
    for s in SAMPLES:
        print(f"  {s}: {sample_den[s]} cells")
    print()

    standard_metrics = load_standard_gene_metrics()

    engine = import_verified_engine()

    h5dump = engine.find_executable(
        ["h5dump"],
        preferred=["/opt/homebrew/bin/h5dump"],
    )
    compiler = engine.find_executable(["clang", "cc", "gcc"])

    print("h5dump :", h5dump)
    print("compiler:", compiler)
    print()

    helper = engine.compile_helper(compiler)

    all_sample_rows = []
    all_lineage_rows = []
    all_probe_rows = []
    all_audit_rows = []

    for sample in SAMPLES:
        print()
        print("-" * 72)
        print(f"Processing {sample}")
        print("-" * 72)
        print("H5:", H5_FILES[sample])

        sample_tmp = TMP_DIR / sample
        if sample_tmp.exists():
            shutil.rmtree(sample_tmp)
        sample_tmp.mkdir(parents=True, exist_ok=True)

        try:
            paths = engine.export_sample_arrays(
                h5dump,
                sample,
                H5_FILES[sample],
                sample_tmp,
            )

            target_file = engine.make_target_file(
                probe_meta[sample],
                sample_tmp,
            )

            selected_file = engine.scan_selected(
                helper,
                paths,
                target_file,
                sample_tmp,
            )

            wanted_barcode_indices = (
                engine.collect_selected_barcode_indices(selected_file)
            )

            barcode_strings = (
                engine.resolve_selected_barcodes_from_h5dump(
                    h5dump,
                    H5_FILES[sample],
                    wanted_barcode_indices,
                )
            )

            barcode_cache = {}
            canonical_selected_indices = set()
            background_selected_rows = 0
            canonical_selected_rows = 0

            gene_keys = {
                g: set()
                for g in strict_genes
            }

            probe_keys = defaultdict(set)

            rec = struct.Struct("<iQI")

            with selected_file.open("rb") as fh:
                while True:
                    raw = fh.read(16 * 250000)

                    if not raw:
                        break

                    if len(raw) % 16 != 0:
                        fail(f"{sample}: partial selected-record block")

                    for probe_idx, barcode_idx, umi in rec.iter_unpack(raw):
                        pm = probe_meta[sample].get(probe_idx)

                        if pm is None:
                            fail(
                                f"{sample}: unknown selected "
                                f"probe_idx={probe_idx}"
                            )

                        cached = barcode_cache.get(barcode_idx)

                        if cached is None:
                            raw_bc = barcode_strings.get(barcode_idx)

                            if raw_bc is None:
                                fail(
                                    f"{sample}: selected barcode_idx="
                                    f"{barcode_idx} unresolved from /barcodes"
                                )

                            core = engine.extract_barcode_core(raw_bc)
                            meta = (
                                context[sample].get(core)
                                if core
                                else None
                            )

                            cached = (core, meta)
                            barcode_cache[barcode_idx] = cached

                        core, meta = cached

                        if meta is None:
                            background_selected_rows += 1
                            continue

                        canonical_selected_rows += 1
                        canonical_selected_indices.add(barcode_idx)

                        g = pm["gene"]
                        pid = pm["probe_id"]
                        key = engine.mol_key(barcode_idx, umi)

                        gene_keys[g].add(key)
                        probe_keys[(g, pid)].add(key)

            # Sample-level full denominator
            for g in strict_genes:
                keys = gene_keys[g]
                cells = {k >> 32 for k in keys}
                n_full = sample_den[sample]

                all_sample_rows.append(
                    {
                        "sample": sample,
                        "gene": g,
                        "full_canonical_cells": n_full,
                        "strict_shared_UMI": len(keys),
                        "strict_shared_positive_cells": len(cells),
                        "strict_shared_positive_fraction":
                            len(cells) / n_full if n_full else 0.0,
                        "strict_shared_UMI_per_1000_cells":
                            len(keys) / n_full * 1000.0 if n_full else 0.0,
                        "n_strict_shared_probe_ids":
                            len(shared_ids[g]),
                        "strict_shared_probe_ids":
                            ";".join(sorted(shared_ids[g])),
                    }
                )

                # Full lineage denominator, including zero-expressing cells.
                lineages = sorted(lineage_den[sample])

                # Group selected keys by mapped lineage.
                grouped_keys = defaultdict(set)

                for key in keys:
                    barcode_idx = key >> 32
                    meta = barcode_cache[barcode_idx][1]

                    if meta is None:
                        fail(
                            f"{sample} {g}: canonical key lacks metadata."
                        )

                    grouped_keys[meta["annotation"]].add(key)

                for lineage in lineages:
                    lkeys = grouped_keys.get(lineage, set())
                    lcells = {k >> 32 for k in lkeys}
                    n_full_lineage = lineage_den[sample][lineage]

                    std = standard_metrics.get(
                        (sample, lineage, g)
                    )

                    standard_umi = (
                        std["standard_total_UMI"]
                        if std is not None
                        else None
                    )

                    standard_pos = (
                        std["standard_positive_cells"]
                        if std is not None
                        else None
                    )

                    standard_rate = (
                        std["standard_UMI_per_1000_cells"]
                        if std is not None
                        else None
                    )

                    strict_rate = (
                        len(lkeys) / n_full_lineage * 1000.0
                        if n_full_lineage
                        else 0.0
                    )

                    if standard_umi is None or standard_umi == 0:
                        retained = None
                    else:
                        retained = len(lkeys) / standard_umi

                    all_lineage_rows.append(
                        {
                            "sample": sample,
                            "gene": g,
                            "lineage": lineage,
                            "full_lineage_n_cells": n_full_lineage,
                            "strict_shared_UMI": len(lkeys),
                            "strict_shared_positive_cells": len(lcells),
                            "strict_shared_positive_fraction":
                                len(lcells) / n_full_lineage
                                if n_full_lineage else 0.0,
                            "strict_shared_UMI_per_1000_full_lineage_cells":
                                strict_rate,
                            "standard_matrix_UMI":
                                "" if standard_umi is None else standard_umi,
                            "standard_matrix_positive_cells":
                                "" if standard_pos is None else standard_pos,
                            "standard_matrix_UMI_per_1000_cells":
                                "" if standard_rate is None else standard_rate,
                            "strict_to_standard_UMI_fraction":
                                "" if retained is None else retained,
                            "n_strict_shared_probe_ids":
                                len(shared_ids[g]),
                            "strict_shared_probe_ids":
                                ";".join(sorted(shared_ids[g])),
                        }
                    )

                for pid in sorted(shared_ids[g]):
                    pkeys = probe_keys[(g, pid)]
                    all_probe_rows.append(
                        {
                            "sample": sample,
                            "gene": g,
                            "probe_id": pid,
                            "canonical_strict_UMI": len(pkeys),
                            "canonical_strict_positive_cells":
                                len({k >> 32 for k in pkeys}),
                        }
                    )

            all_audit_rows.append(
                {
                    "sample": sample,
                    "full_canonical_cells": sample_den[sample],
                    "canonical_barcode_indices_with_any_strict_signal":
                        len(canonical_selected_indices),
                    "canonical_selected_molecule_rows":
                        canonical_selected_rows,
                    "background_selected_molecule_rows":
                        background_selected_rows,
                    "selected_barcode_indices_cached":
                        len(barcode_cache),
                    "n_target_probe_indices":
                        len(probe_meta[sample]),
                }
            )

            print(
                f"PASS: {sample}: "
                f"{len(canonical_selected_indices)} canonical cells "
                "with >=1 strict-target molecule"
            )

        finally:
            if sample_tmp.exists():
                print(f"  cleaning temporary binaries for {sample}")
                shutil.rmtree(sample_tmp, ignore_errors=True)

    # ------------------------------------------------------------------
    # Write reconstruction tables
    # ------------------------------------------------------------------

    write_tsv(
        TABLE_DIR / f"{VERSION}_strict_gene_reconstruction_by_sample.tsv",
        all_sample_rows,
        [
            "sample",
            "gene",
            "full_canonical_cells",
            "strict_shared_UMI",
            "strict_shared_positive_cells",
            "strict_shared_positive_fraction",
            "strict_shared_UMI_per_1000_cells",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_strict_lineage_reconstruction_FULL_DENOM.tsv",
        all_lineage_rows,
        [
            "sample",
            "gene",
            "lineage",
            "full_lineage_n_cells",
            "strict_shared_UMI",
            "strict_shared_positive_cells",
            "strict_shared_positive_fraction",
            "strict_shared_UMI_per_1000_full_lineage_cells",
            "standard_matrix_UMI",
            "standard_matrix_positive_cells",
            "standard_matrix_UMI_per_1000_cells",
            "strict_to_standard_UMI_fraction",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_strict_probe_composition_by_sample.tsv",
        all_probe_rows,
        [
            "sample",
            "gene",
            "probe_id",
            "canonical_strict_UMI",
            "canonical_strict_positive_cells",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_barcode_mapping_audit.tsv",
        all_audit_rows,
        [
            "sample",
            "full_canonical_cells",
            "canonical_barcode_indices_with_any_strict_signal",
            "canonical_selected_molecule_rows",
            "background_selected_molecule_rows",
            "selected_barcode_indices_cached",
            "n_target_probe_indices",
        ],
    )

    # ------------------------------------------------------------------
    # Comparable-lineage 2-vs-2 descriptive transitions
    # ------------------------------------------------------------------

    by_gls = {}
    for r in all_lineage_rows:
        if r["lineage"] not in COMPARABLE_LINEAGES:
            continue
        key = (r["gene"], r["lineage"], r["sample"])
        by_gls[key] = r

    transition_rows = []

    for g in strict_genes:
        for lineage in sorted(COMPARABLE_LINEAGES):
            vals = {}

            for s in SAMPLES:
                r = by_gls.get((g, lineage, s))
                if r is None:
                    vals[s] = None
                else:
                    vals[s] = float(
                        r[
                            "strict_shared_UMI_per_1000_full_lineage_cells"
                        ]
                    )

            rel = relation_2v2(
                vals["SfLB1"],
                vals["SfLB2"],
                vals["Tx17"],
                vals["Tx5"],
            )

            smean = (
                (vals["SfLB1"] + vals["SfLB2"]) / 2
                if vals["SfLB1"] is not None
                and vals["SfLB2"] is not None
                else None
            )

            tmean = (
                (vals["Tx17"] + vals["Tx5"]) / 2
                if vals["Tx17"] is not None
                and vals["Tx5"] is not None
                else None
            )

            transition_rows.append(
                {
                    "gene": g,
                    "lineage": lineage,
                    "SfLB1_strict_UMI_per_1000": vals["SfLB1"],
                    "SfLB2_strict_UMI_per_1000": vals["SfLB2"],
                    "Tx17_strict_UMI_per_1000": vals["Tx17"],
                    "Tx5_strict_UMI_per_1000": vals["Tx5"],
                    "SfLB_mean_strict_UMI_per_1000":
                        "" if smean is None else smean,
                    "Tx_mean_strict_UMI_per_1000":
                        "" if tmean is None else tmean,
                    "Tx_minus_SfLB":
                        "" if smean is None or tmean is None
                        else tmean - smean,
                    "replicate_relation": rel,
                    "n_strict_shared_probe_ids":
                        len(shared_ids[g]),
                    "strict_shared_probe_ids":
                        ";".join(sorted(shared_ids[g])),
                }
            )

    write_tsv(
        TABLE_DIR / f"{VERSION}_comparable_lineage_strict_gene_transition.tsv",
        transition_rows,
        [
            "gene",
            "lineage",
            "SfLB1_strict_UMI_per_1000",
            "SfLB2_strict_UMI_per_1000",
            "Tx17_strict_UMI_per_1000",
            "Tx5_strict_UMI_per_1000",
            "SfLB_mean_strict_UMI_per_1000",
            "Tx_mean_strict_UMI_per_1000",
            "Tx_minus_SfLB",
            "replicate_relation",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    # ------------------------------------------------------------------
    # Dedicated antigen-presentation audit.
    # TAP2 should be absent because it had no strict shared subset.
    # ------------------------------------------------------------------

    antigen_rows = [
        r for r in transition_rows
        if r["gene"] in ANTIGEN_GENES
    ]

    write_tsv(
        TABLE_DIR / f"{VERSION}_antigen_strict_gene_transition.tsv",
        antigen_rows,
        [
            "gene",
            "lineage",
            "SfLB1_strict_UMI_per_1000",
            "SfLB2_strict_UMI_per_1000",
            "Tx17_strict_UMI_per_1000",
            "Tx5_strict_UMI_per_1000",
            "SfLB_mean_strict_UMI_per_1000",
            "Tx_mean_strict_UMI_per_1000",
            "Tx_minus_SfLB",
            "replicate_relation",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------

    rel_counts = defaultdict(int)
    for r in transition_rows:
        rel_counts[r["replicate_relation"]] += 1

    print()
    print("=" * 72)
    print("STRICT SHARED-PROBE RECONSTRUCTION SUMMARY")
    print("=" * 72)
    print(f"Strict genes reconstructed: {len(strict_genes)}")
    print(f"Comparable lineages: {', '.join(sorted(COMPARABLE_LINEAGES))}")
    print()
    print("Gene x comparable-lineage replicate relations:")
    for k in sorted(rel_counts):
        print(f"  {k}: {rel_counts[k]}")

    print()
    print("=== ANTIGEN-PRESENTATION STRICT RESULTS ===")
    for r in antigen_rows:
        print(
            f"{r['gene']} | {r['lineage']} | "
            f"SfLB1={float(r['SfLB1_strict_UMI_per_1000']):.4f}, "
            f"SfLB2={float(r['SfLB2_strict_UMI_per_1000']):.4f}, "
            f"Tx17={float(r['Tx17_strict_UMI_per_1000']):.4f}, "
            f"Tx5={float(r['Tx5_strict_UMI_per_1000']):.4f} | "
            f"{r['replicate_relation']}"
        )

    if any(r["gene"] == "TAP2" for r in antigen_rows):
        fail(
            "TAP2 unexpectedly entered strict reconstruction despite "
            "P.1 class lacking a strict shared subset."
        )

    print("TAP2: correctly excluded from strict reconstruction.")

    print()
    print("=" * 72)
    print(f"{VERSION} COMPLETE")
    print("=" * 72)
    print("Output:", OUT)
    print("Full frozen sample/lineage denominators used: YES")
    print("Strict shared probe molecule reconstruction: COMPLETE")
    print("No multi-gene program score rebuilt yet.")
    print(
        "Next step: inspect strict gene transitions, then rebuild "
        "probe-comparable program scores."
    )


if __name__ == "__main__":
    main()
