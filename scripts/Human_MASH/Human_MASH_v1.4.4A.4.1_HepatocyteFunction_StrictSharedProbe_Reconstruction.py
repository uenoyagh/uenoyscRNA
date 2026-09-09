#!/usr/bin/env python3
# ==============================================================================
# Human_MASH_v1.4.4A.4.1_HepatocyteFunction_StrictSharedProbe_Reconstruction.py
#
# Purpose
# -------
# Molecule-level reconstruction of the 16 Flex-measurable hepatic
# function/maturity genes after v1.4.4A.4.0 established that ALL 16 are
# strict shared-probe comparable across SfLB1/SfLB2/Tx17/Tx5.
#
# Uses only probe IDs that are:
#   - present in all four samples,
#   - included_probes == TRUE in all four samples,
#   - filtered_probes == TRUE in all four samples.
#
# Full frozen canonical denominators are used for sample, lineage and state
# summaries. Zero-expressing canonical cells remain in denominators.
#
# This is a descriptive reconstruction/QC step only:
#   - no inferential test (n=2/group)
#   - no final aggregate functional score
#   - no claim that transcript levels equal functional protein/activity
#
# The verified v1.4.3B.2C.3 low-level scanner and barcode resolver are imported
# and reused rather than reimplemented.
# ==============================================================================

from __future__ import annotations

import csv
import importlib.util
import math
import shutil
import struct
from collections import defaultdict
from pathlib import Path

VERSION = "v1.4.4A.4.1"

REPO = Path("/Users/uenoya/Projects/uenoyscRNA")

BASE_ENGINE = (
    REPO
    / "scripts"
    / "Human_MASH"
    / "Human_MASH_v1.4.3B.2C.3_PrioritizedLigand_MoleculeReconstruction.py"
)

BASE = Path("/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk")

P0_DIR = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.4.0"
    / "analysis"
    / "v1.4.4A.4.0_hepatocyte_function_probe_comparability_preflight"
)

P0_LONG = P0_DIR / "v1.4.4A.4.0_hepatocyte_function_probe_status_long.tsv"
P0_CROSS = P0_DIR / "v1.4.4A.4.0_hepatocyte_function_crosssample_probe_summary.tsv"

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

OUT = (
    BASE
    / "Human_MASH_RDS"
    / "Human_MASH_12sample_Context_v1.4.4A.4.1"
    / "analysis"
    / "v1.4.4A.4.1_hepatocyte_function_strict_shared_probe_reconstruction"
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

PANELS = {
    "Secretory_synthetic": ["TTR", "APOA1", "FGB", "FGG"],
    "Urea_nitrogen": ["OTC", "ASS1", "ASL", "NAGS"],
    "Bile_transport": ["ABCB11"],
    "Lipid_metabolic": ["MTTP", "FABP1"],
    "Hepatocyte_identity": ["HNF4A", "CEBPA"],
    "Epithelial_immaturity_transition": ["EPCAM", "KRT7", "SOX9"],
}

CANDIDATES = sorted({g for genes in PANELS.values() for g in genes})

PRIMARY_LINEAGE = "Hepatocyte_like"
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

    mod.TMP_DIR = TMP_DIR
    return mod


def load_strict_probe_definition():
    cross = read_tsv(P0_CROSS)
    long_rows = read_tsv(P0_LONG)

    if not cross or not long_rows:
        fail("A.4.0 probe-comparability tables are empty.")

    by_gene = {
        r["gene"].strip().upper(): r
        for r in cross
    }

    missing = sorted(set(CANDIDATES) - set(by_gene))
    if missing:
        fail("A.4.0 cross table missing genes: " + ", ".join(missing))

    non_strict = [
        g for g in CANDIDATES
        if by_gene[g]["technical_class"].strip()
        != "STRICT_SHARED_PROBE_COMPARABLE"
    ]

    if non_strict:
        fail(
            "Expected all 16 hepatic genes to be strict-comparable; "
            "non-strict: " + ", ".join(non_strict)
        )

    shared_ids = {}
    for g in CANDIDATES:
        raw = by_gene[g]["strict_shared_probe_ids_all4"].strip()
        ids = {x for x in raw.split(";") if x}
        if not ids:
            fail(f"{g}: no strict shared probe IDs despite strict class.")
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
                f"{s} {g} {pid}: strict probe lacks "
                "included=TRUE & filtered=TRUE."
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
                "panel": ";".join(
                    p for p, genes in PANELS.items()
                    if g in genes
                ),
                "probe_index_0based": idx,
                "probe_id": pid,
                "strict_shared_all4": "TRUE",
            }
        )

    for s in SAMPLES:
        seen = {x["gene"] for x in probe_meta[s].values()}
        missing = sorted(set(CANDIDATES) - seen)
        if missing:
            fail(
                f"{s}: candidate genes missing strict selected probes: "
                + ", ".join(missing)
            )

    selected_rows.sort(
        key=lambda r: (
            r["gene"],
            SAMPLES.index(r["sample"]),
            int(r["probe_index_0based"]),
        )
    )

    return shared_ids, probe_meta, selected_rows


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


def validate_lineages_against_a11(lineage_den):
    rows = read_csv(A11_COUNTS)

    if not rows:
        fail("A.1.1 lineage table empty.")

    required = {"sample", "category", "n_cells"}
    missing = required - set(rows[0])
    if missing:
        fail(
            "A.1.1 lineage table missing: "
            + ", ".join(sorted(missing))
        )

    expected = {}
    for r in rows:
        s = norm_sample(r["sample"])
        if s in SAMPLES:
            expected[(s, r["category"].strip())] = int(float(r["n_cells"]))

    observed = {}
    for s in SAMPLES:
        for lineage, n in lineage_den[s].items():
            observed[(s, lineage)] = n

    all_keys = set(expected) | set(observed)
    bad = []

    for key in sorted(all_keys):
        a = expected.get(key, 0)
        b = observed.get(key, 0)
        if a != b:
            bad.append((key, a, b))

    if bad:
        txt = "; ".join(
            f"{s}/{ll}:A11={a},context={b}"
            for ((s, ll), a, b) in bad[:30]
        )
        fail("Full lineage denominator mismatch vs A.1.1: " + txt)


def summarize_group(
    sample,
    gene,
    keys,
    barcode_cache,
    group_key,
    full_denominators,
):
    grouped = defaultdict(set)

    for key in keys:
        barcode_idx = key >> 32
        cached = barcode_cache.get(barcode_idx)

        if cached is None or cached[1] is None:
            fail(
                f"{sample}/{gene}: canonical key missing cached metadata."
            )

        meta = cached[1]
        grouped[meta[group_key]].add(key)

    rows = []

    for group in sorted(full_denominators):
        gkeys = grouped.get(group, set())
        cells = {k >> 32 for k in gkeys}
        n_full = full_denominators[group]

        rows.append(
            {
                "sample": sample,
                "gene": gene,
                "group": group,
                "full_n_cells": n_full,
                "strict_shared_UMI": len(gkeys),
                "strict_shared_positive_cells": len(cells),
                "strict_shared_positive_fraction":
                    len(cells) / n_full if n_full else 0.0,
                "strict_shared_UMI_per_1000_cells":
                    len(gkeys) / n_full * 1000.0 if n_full else 0.0,
            }
        )

    return rows


def main():
    print("=" * 74)
    print(
        f"Human MASH {VERSION} hepatic-function strict shared-probe reconstruction"
    )
    print("=" * 74)
    print()
    print("All 16 A.4.0 hepatic genes are expected to be strict-comparable.")
    print("Full frozen sample / lineage / state denominators will be used.")
    print()

    required_paths = [
        BASE_ENGINE,
        P0_LONG,
        P0_CROSS,
        CONTEXT_META,
        A11_COUNTS,
        *H5_FILES.values(),
    ]

    for p in required_paths:
        if not p.exists():
            fail(f"Required input missing: {p}")

    OUT.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    if shutil.disk_usage(OUT).free < MIN_FREE_BYTES:
        fail("Insufficient free disk space (<3 GiB).")

    shared_ids, probe_meta, selected_rows = load_strict_probe_definition()

    print(f"Strict-comparable hepatic genes: {len(CANDIDATES)}/16")
    print(
        "Strict shared probes/gene range: "
        f"{min(len(shared_ids[g]) for g in CANDIDATES)}–"
        f"{max(len(shared_ids[g]) for g in CANDIDATES)}"
    )
    print()

    write_tsv(
        TABLE_DIR / f"{VERSION}_strict_probe_definition_used.tsv",
        selected_rows,
        [
            "sample",
            "gene",
            "panel",
            "probe_index_0based",
            "probe_id",
            "strict_shared_all4",
        ],
    )

    context = load_context_full()
    sample_den, lineage_den, state_den = build_full_denominators(context)

    validate_lineages_against_a11(lineage_den)

    print("Canonical context cell counts: PASS")
    print("Full lineage denominators vs A.1.1: PASS")
    print()

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
    all_state_rows = []
    all_probe_rows = []
    audit_rows = []

    for sample in SAMPLES:
        print()
        print("-" * 74)
        print(f"Processing {sample}")
        print("-" * 74)
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
            canonical_selected_rows = 0
            background_selected_rows = 0

            gene_keys = {g: set() for g in CANDIDATES}
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
                                    f"{barcode_idx} unresolved."
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

                        gene = pm["gene"]
                        pid = pm["probe_id"]
                        key = engine.mol_key(barcode_idx, umi)

                        gene_keys[gene].add(key)
                        probe_keys[(gene, pid)].add(key)

            # sample-level
            for gene in CANDIDATES:
                keys = gene_keys[gene]
                cells = {k >> 32 for k in keys}
                n_full = sample_den[sample]

                all_sample_rows.append(
                    {
                        "sample": sample,
                        "gene": gene,
                        "panel": ";".join(
                            p for p, genes in PANELS.items()
                            if gene in genes
                        ),
                        "full_canonical_cells": n_full,
                        "strict_shared_UMI": len(keys),
                        "strict_shared_positive_cells": len(cells),
                        "strict_shared_positive_fraction":
                            len(cells) / n_full if n_full else 0.0,
                        "strict_shared_UMI_per_1000_cells":
                            len(keys) / n_full * 1000.0 if n_full else 0.0,
                        "n_strict_shared_probe_ids":
                            len(shared_ids[gene]),
                        "strict_shared_probe_ids":
                            ";".join(sorted(shared_ids[gene])),
                    }
                )

                lrows = summarize_group(
                    sample,
                    gene,
                    keys,
                    barcode_cache,
                    "annotation",
                    lineage_den[sample],
                )

                for r in lrows:
                    r["panel"] = ";".join(
                        p for p, genes in PANELS.items()
                        if gene in genes
                    )
                    r["n_strict_shared_probe_ids"] = len(shared_ids[gene])
                    r["strict_shared_probe_ids"] = ";".join(
                        sorted(shared_ids[gene])
                    )
                    all_lineage_rows.append(r)

                srows = summarize_group(
                    sample,
                    gene,
                    keys,
                    barcode_cache,
                    "state",
                    state_den[sample],
                )

                for r in srows:
                    r["panel"] = ";".join(
                        p for p, genes in PANELS.items()
                        if gene in genes
                    )
                    r["n_strict_shared_probe_ids"] = len(shared_ids[gene])
                    r["strict_shared_probe_ids"] = ";".join(
                        sorted(shared_ids[gene])
                    )
                    all_state_rows.append(r)

                for pid in sorted(shared_ids[gene]):
                    pkeys = probe_keys[(gene, pid)]
                    all_probe_rows.append(
                        {
                            "sample": sample,
                            "gene": gene,
                            "probe_id": pid,
                            "canonical_strict_UMI": len(pkeys),
                            "canonical_strict_positive_cells":
                                len({k >> 32 for k in pkeys}),
                        }
                    )

            audit_rows.append(
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
                "with >=1 hepatic strict-target molecule"
            )

        finally:
            if sample_tmp.exists():
                print(f"  cleaning temporary binaries for {sample}")
                shutil.rmtree(sample_tmp, ignore_errors=True)

    # ------------------------------------------------------------------
    # Write reconstruction
    # ------------------------------------------------------------------

    write_tsv(
        TABLE_DIR / f"{VERSION}_hepatic_gene_reconstruction_by_sample.tsv",
        all_sample_rows,
        [
            "sample",
            "gene",
            "panel",
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
        TABLE_DIR / f"{VERSION}_hepatic_lineage_reconstruction_FULL_DENOM.tsv",
        all_lineage_rows,
        [
            "sample",
            "gene",
            "panel",
            "group",
            "full_n_cells",
            "strict_shared_UMI",
            "strict_shared_positive_cells",
            "strict_shared_positive_fraction",
            "strict_shared_UMI_per_1000_cells",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_hepatic_state_reconstruction_FULL_DENOM.tsv",
        all_state_rows,
        [
            "sample",
            "gene",
            "panel",
            "group",
            "full_n_cells",
            "strict_shared_UMI",
            "strict_shared_positive_cells",
            "strict_shared_positive_fraction",
            "strict_shared_UMI_per_1000_cells",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    write_tsv(
        TABLE_DIR / f"{VERSION}_hepatic_probe_composition_by_sample.tsv",
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
        audit_rows,
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
    # Primary Hepatocyte_like 2-vs-2 gene transition
    # ------------------------------------------------------------------

    by_key = {}

    for r in all_lineage_rows:
        if r["group"] != PRIMARY_LINEAGE:
            continue

        key = (r["gene"], r["sample"])

        if key in by_key:
            fail(f"Duplicated primary lineage row: {key}")

        by_key[key] = r

    transition_rows = []

    for gene in CANDIDATES:
        vals = {}

        for sample in SAMPLES:
            r = by_key.get((gene, sample))

            if r is None:
                fail(
                    f"Missing {PRIMARY_LINEAGE} row for "
                    f"{gene}/{sample}"
                )

            vals[sample] = float(
                r["strict_shared_UMI_per_1000_cells"]
            )

        sf_mean = (vals["SfLB1"] + vals["SfLB2"]) / 2
        tx_mean = (vals["Tx17"] + vals["Tx5"]) / 2

        transition_rows.append(
            {
                "gene": gene,
                "panel": ";".join(
                    p for p, genes in PANELS.items()
                    if gene in genes
                ),
                "lineage": PRIMARY_LINEAGE,
                "SfLB1_strict_UMI_per_1000": vals["SfLB1"],
                "SfLB2_strict_UMI_per_1000": vals["SfLB2"],
                "Tx17_strict_UMI_per_1000": vals["Tx17"],
                "Tx5_strict_UMI_per_1000": vals["Tx5"],
                "SfLB_mean_strict_UMI_per_1000": sf_mean,
                "Tx_mean_strict_UMI_per_1000": tx_mean,
                "Tx_minus_SfLB": tx_mean - sf_mean,
                "Tx_vs_SfLB_log2FC_pseudocount_0p1":
                    math.log2((tx_mean + 0.1) / (sf_mean + 0.1)),
                "replicate_relation": relation_2v2(
                    vals["SfLB1"],
                    vals["SfLB2"],
                    vals["Tx17"],
                    vals["Tx5"],
                ),
                "n_strict_shared_probe_ids":
                    len(shared_ids[gene]),
                "strict_shared_probe_ids":
                    ";".join(sorted(shared_ids[gene])),
            }
        )

    write_tsv(
        TABLE_DIR / f"{VERSION}_Hepatocyte_like_strict_gene_transition.tsv",
        transition_rows,
        [
            "gene",
            "panel",
            "lineage",
            "SfLB1_strict_UMI_per_1000",
            "SfLB2_strict_UMI_per_1000",
            "Tx17_strict_UMI_per_1000",
            "Tx5_strict_UMI_per_1000",
            "SfLB_mean_strict_UMI_per_1000",
            "Tx_mean_strict_UMI_per_1000",
            "Tx_minus_SfLB",
            "Tx_vs_SfLB_log2FC_pseudocount_0p1",
            "replicate_relation",
            "n_strict_shared_probe_ids",
            "strict_shared_probe_ids",
        ],
    )

    # ------------------------------------------------------------------
    # Panel-level gene-direction counts in Hepatocyte_like
    # ------------------------------------------------------------------

    panel_summary = []

    for panel, genes in PANELS.items():
        z = [r for r in transition_rows if r["gene"] in genes]

        higher = [
            r["gene"] for r in z
            if r["replicate_relation"]
            == "Tx_higher_complete_replicate_range_separation"
        ]

        lower = [
            r["gene"] for r in z
            if r["replicate_relation"]
            == "Tx_lower_complete_replicate_range_separation"
        ]

        overlap = [
            r["gene"] for r in z
            if r["replicate_relation"]
            == "Replicate_ranges_overlap"
        ]

        panel_summary.append(
            {
                "panel": panel,
                "n_genes": len(z),
                "n_Tx_higher_complete": len(higher),
                "Tx_higher_complete_genes": ";".join(sorted(higher)),
                "n_Tx_lower_complete": len(lower),
                "Tx_lower_complete_genes": ";".join(sorted(lower)),
                "n_overlap": len(overlap),
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )

    write_tsv(
        TABLE_DIR / f"{VERSION}_Hepatocyte_like_panel_direction_summary.tsv",
        panel_summary,
        [
            "panel",
            "n_genes",
            "n_Tx_higher_complete",
            "Tx_higher_complete_genes",
            "n_Tx_lower_complete",
            "Tx_lower_complete_genes",
            "n_overlap",
            "overlap_genes",
        ],
    )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------

    print()
    print("=" * 74)
    print("HEPATOCYTE-LIKE STRICT GENE TRANSITION")
    print("=" * 74)

    for r in panel_summary:
        print(
            f"{r['panel']}: "
            f"Tx-higher={r['n_Tx_higher_complete']}, "
            f"Tx-lower={r['n_Tx_lower_complete']}, "
            f"overlap={r['n_overlap']}"
        )
        print(
            "  higher: "
            + (r["Tx_higher_complete_genes"] or "(none)")
        )
        print(
            "  lower : "
            + (r["Tx_lower_complete_genes"] or "(none)")
        )
        print(
            "  overlap: "
            + (r["overlap_genes"] or "(none)")
        )

    print()
    print("=" * 74)
    print(f"{VERSION} COMPLETE")
    print("=" * 74)
    print("Output:", OUT)
    print("Full frozen denominators used: YES")
    print("Strict shared-probe molecule reconstruction: COMPLETE")
    print("Primary lineage: Hepatocyte_like")
    print("No inferential test performed (n=2/group).")
    print("No aggregate hepatic-function score created yet.")


if __name__ == "__main__":
    main()
