#!/usr/bin/env python3
"""Audit redshift/richness dependence of the approved HSC count-only preflight.

This analysis reads only one integer count per lens. It never submits a CAS
query and never reads or exports source-row, coordinate, or photo-z data.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from submit_wl_shape_counts import load_plan, read_count_csv

DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PREFLIGHT = DEFAULT_ROOT / "wl_shape_preflight_2026-10-18_01"
DEFAULT_OUTPUT = DEFAULT_ROOT / "wl_shape_count_selection_audit_2026-10-18_01"
Z_EDGES = (0.29, 0.35, 0.40, 0.45, 0.50, 0.57)
RICHNESS_EDGES = (15.0, 20.0, 30.0, 40.0, 50.0, 110.0)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect_counts(preflight):
    preflight = Path(preflight).resolve()
    plan = load_plan(preflight)
    result_dir = preflight / "cas_counts"
    state_path = result_dir / "count_status.json"
    summary_path = result_dir / "count_summary.json"
    state = json.loads(state_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (state.get("status") != "complete" or summary.get("status") != "complete"
            or state.get("source_rows_downloaded") is not False
            or summary.get("source_rows_downloaded") is not False
            or state.get("science_measurement") is not False
            or summary.get("science_measurement") is not False
            or summary.get("query_plan_sha256") != state.get("plan_sha256")):
        raise ValueError("Expected a completed count-only run with no source-row transfer")
    jobs = state.get("jobs", {})
    if len(jobs) != 32 or any(job.get("status") != "complete" for job in jobs.values()):
        raise ValueError("All 32 count-only SQL batches must be complete")

    counts = {}
    count_hashes = {}
    batch_size = int(plan["batch_size"])
    for index, item in enumerate(plan["sql_files"]):
        start = index * batch_size
        lenses = plan["lenses"][start:start + batch_size]
        ids = [str(row["lens_id"]) for row in lenses]
        path = result_dir / f"shape_counts_{index:03d}.csv"
        records = read_count_csv(path, ids)
        count_hashes[path.name] = sha256(path)
        for record in records:
            counts[str(record["lens_id"])] = int(record["n_shape_sources"])

    lenses = {str(row["lens_id"]): row for row in plan["lenses"]}
    if (len(counts) != 800 or set(counts) != set(lenses)
            or sum(counts.values()) != int(summary["sum_lens_cone_shape_counts"])):
        raise ValueError("Verified count rows do not match all 800 planned lenses")
    rows = []
    for lens_id, lens in lenses.items():
        # Richness is joined from the immutable sample manifest below.
        rows.append({"lens_id": lens_id, "z_lens": float(lens["z_lens"]),
                     "n_shape_sources": counts[lens_id]})
    return rows, state, summary, count_hashes, plan


def summarize_group(rows, variable, low, high, include_high=False):
    selected = [row for row in rows
                if row[variable] >= low and (row[variable] <= high if include_high
                                             else row[variable] < high)]
    values = np.asarray([row["n_shape_sources"] for row in selected], dtype=np.int64)
    if values.size == 0:
        raise ValueError(f"Empty {variable} group [{low}, {high}]")
    return {
        "dimension": variable,
        "lower": float(low),
        "upper": float(high),
        "n_lenses": int(values.size),
        "zero_count_lenses": int(np.count_nonzero(values == 0)),
        "zero_fraction": float(np.mean(values == 0)),
        "median_candidates_per_lens": float(np.median(values)),
        "p10_candidates_per_lens": float(np.quantile(values, 0.10)),
        "p90_candidates_per_lens": float(np.quantile(values, 0.90)),
        "lenses_ge_1000_candidates": int(np.count_nonzero(values >= 1000)),
    }


def analyze(preflight, sample_manifest, output):
    preflight = Path(preflight).resolve()
    sample_manifest = Path(sample_manifest).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit output: {output}")

    manifest = json.loads(sample_manifest.read_text(encoding="utf-8"))
    cluster_rows = {str(row["id"]): row for row in manifest["clusters"]}
    rows, state, count_summary, count_hashes, plan = collect_counts(preflight)
    if set(cluster_rows) != {row["lens_id"] for row in rows}:
        raise ValueError("Sample manifest and count-only lens IDs do not match")
    for row in rows:
        cluster = cluster_rows[row["lens_id"]]
        row["richness"] = float(cluster["richness"])
        row["sample_z_cl"] = float(cluster["z_cl"])
        if not np.isclose(row["z_lens"], float(cluster["z"]), atol=5e-7, rtol=0):
            raise ValueError(f"Geometry redshift mismatch for lens {row['lens_id']}")

    counts = np.asarray([row["n_shape_sources"] for row in rows], dtype=np.int64)
    z = np.asarray([row["z_lens"] for row in rows], dtype=float)
    richness = np.asarray([row["richness"] for row in rows], dtype=float)
    z_corr = float(spearmanr(counts, z).statistic)
    richness_corr = float(spearmanr(counts, richness).statistic)
    groups = []
    z_edges = list(Z_EDGES)
    richness_edges = list(RICHNESS_EDGES)
    for edges, variable in ((z_edges, "z_lens"), (richness_edges, "richness")):
        for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
            groups.append(summarize_group(
                rows, variable, low, high, include_high=index == len(edges) - 2,
            ))

    output.mkdir(parents=True, exist_ok=False)
    groups_path = output / "count_groups.tsv"
    fields = list(groups[0])
    with groups_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(groups)
    summary_path = output / "summary.json"
    summary = {
        "schema": "hsc-wl-count-selection-audit-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "source_rows_downloaded": False,
        "photoz_or_shape_values_returned": False,
        "n_lenses": len(rows),
        "count_rows": len(counts),
        "sum_lens_cone_candidates": int(counts.sum()),
        "zero_count_lenses": int(np.count_nonzero(counts == 0)),
        "zero_count_fraction": float(np.mean(counts == 0)),
        "candidate_count_distribution": {
            "min": int(counts.min()), "p10": float(np.quantile(counts, .10)),
            "median": float(np.median(counts)), "p90": float(np.quantile(counts, .90)),
            "max": int(counts.max()),
            "lenses_ge_1000": int(np.count_nonzero(counts >= 1000)),
            "lenses_ge_3000": int(np.count_nonzero(counts >= 3000)),
            "lenses_ge_5000": int(np.count_nonzero(counts >= 5000)),
        },
        "spearman_descriptive_only": {
            "count_vs_geometry_redshift": z_corr,
            "count_vs_CAMIRA_richness": richness_corr,
            "p_values_reported": False,
        },
        "group_table": str(groups_path),
        "group_table_sha256": sha256(groups_path),
        "inputs": {
            "sample_manifest": str(sample_manifest),
            "sample_manifest_sha256": sha256(sample_manifest),
            "preflight_manifest": str(preflight / "preflight_manifest.json"),
            "preflight_manifest_sha256": sha256(preflight / "preflight_manifest.json"),
            "count_plan_sha256": state["plan_sha256"],
            "count_summary_sha256": sha256(preflight / "cas_counts" / "count_summary.json"),
            "count_csv_sha256": count_hashes,
            "analysis_script_sha256": sha256(Path(__file__).resolve()),
        },
        "interpretation": [
            "Counts are raw shape-catalog lens-cone candidates before source-quality, tomographic-bin, magnitude, resolution, radial, and random-lens cuts.",
            "The summed count is lens-source pair support; the same source can count for multiple lenses and this is not a unique-source total.",
            "Redshift dependence is descriptive and combines angular aperture, survey footprint/depth, and catalog selection; it is not a physical trend.",
            "This audit does not estimate final effective source counts or lensing signal-to-noise.",
        ],
        "output_sha256": {},
    }
    summary["output_sha256"][str(groups_path)] = sha256(groups_path)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--sample-manifest", type=Path, default=DEFAULT_ROOT / "sample_manifest.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = analyze(args.preflight, args.sample_manifest, args.output)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "lenses": result["n_lenses"], "candidate_pairs": result["sum_lens_cone_candidates"],
        "zero_count_lenses": result["zero_count_lenses"],
        "spearman_rho_count_z": result["spearman_descriptive_only"]["count_vs_geometry_redshift"],
        "spearman_rho_count_richness": result["spearman_descriptive_only"]["count_vs_CAMIRA_richness"],
        "source_rows_downloaded": False, "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
