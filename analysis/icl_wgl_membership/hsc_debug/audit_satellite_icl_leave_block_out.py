#!/usr/bin/env python3
"""Measure connected-sky-block influence on the exploratory ICL association.

This is an influence diagnostic for an existing exploratory correlation. It
does not create a new science measurement or alter the parent profile data.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_satellite_icl_geometric_strata import (
    connected_labels,
    partial_correlation,
    read_tsv,
    sha256,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = ROOT / "satellite_icl_geometric_strata_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "satellite_icl_lobo_2026-10-19_03"


def analyze(input_dir: Path, output_dir: Path):
    input_dir, output_dir = Path(input_dir).resolve(), Path(output_dir).resolve()
    summary_path = input_dir / "summary.json"
    table_path = input_dir / "geometric_strata_pairs.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "hsc-satellite-icl-geometric-strata-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False):
        raise ValueError("unexpected or science-promoted input summary")
    expected_hash = summary.get("outputs_sha256", {}).get(table_path.name)
    if not expected_hash or sha256(table_path) != expected_hash:
        raise ValueError("geometric-strata table does not match its recorded hash")
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_dir}")

    rows = read_tsv(table_path)
    if len(rows) != int(summary.get("n_association_pairs", -1)):
        raise ValueError("geometric-strata row count differs from its summary")
    pair_keys = [(str(row["cluster_id"]), row["filter"]) for row in rows]
    if len(set(pair_keys)) != len(pair_keys):
        raise ValueError("duplicate cluster/filter row in geometric-strata input")
    outputs = []
    summaries = []
    for band in ("i", "i2"):
        for category, key in (
            ("all", "satellite_unweighted"),
            ("no_geometric_overlap", "satellite_no_geometric_overlap"),
            ("overlap_possible", "satellite_overlap_possible"),
        ):
            chosen = [row for row in rows if row["filter"] == band and row.get(key, "")]
            if len(chosen) < 5:
                summaries.append({"filter": band, "stratum": category,
                                  "status": "insufficient_pairs", "n_pairs": len(chosen)})
                continue
            point = partial_correlation(chosen, key)
            labels = connected_labels(chosen)
            unique = np.unique(labels)
            values = []
            for label in unique:
                retained = [row for row, row_label in zip(chosen, labels)
                            if row_label != label]
                status = "complete"
                try:
                    r = partial_correlation(retained, key)
                except ValueError:
                    r = None
                    status = "undefined_after_block_removal"
                delta = None if r is None else r - point
                values.append(r)
                outputs.append({
                    "filter": band,
                    "stratum": category,
                    "omitted_connected_block": str(label),
                    "n_pairs_full": len(chosen),
                    "n_blocks_full": len(unique),
                    "n_pairs_omitted": int(np.count_nonzero(labels == label)),
                    "n_pairs_remaining": len(retained),
                    "partial_r_full": point,
                    "partial_r_leave_block_out": r,
                    "delta_r_leave_block_out": delta,
                    "sign_changed": bool(r * point < 0) if r is not None else None,
                    "status": status,
                })
            valid = [value for value in values if value is not None]
            summaries.append({
                "filter": band,
                "stratum": category,
                "status": "complete_exploratory_influence_audit",
                "n_pairs": len(chosen),
                "n_connected_blocks": len(unique),
                "partial_r_full": point,
                "leave_block_out_valid": len(valid),
                "leave_block_out_undefined": len(values) - len(valid),
                "leave_block_out_min_r": min(valid) if valid else None,
                "leave_block_out_max_r": max(valid) if valid else None,
                "leave_block_out_sign_changes": sum(value * point < 0 for value in valid),
                "interpretation": "influence diagnostic only; no p-value or confirmatory claim",
            })

    output_dir.mkdir(parents=True, exist_ok=False)
    table_out = output_dir / "leave_block_out.tsv"
    fields = list(outputs[0]) if outputs else []
    with table_out.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(outputs)
    payload = {
        "schema": "hsc-satellite-icl-lobo-audit-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "input": {
            "summary_sha256": sha256(summary_path),
            "geometric_strata_pairs_sha256": sha256(table_path),
            "source_summary_bootstrap_replicas": summary["method"]["replicas"],
        },
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "method": "recompute the same z/richness-adjusted partial Pearson r after omitting each connected target/control 5-degree sky block in turn",
        "groups": summaries,
        "leave_block_out_table": str(table_out),
        "leave_block_out_table_sha256": sha256(table_out),
        "limitations": [
            "This tests block influence in the selected exploratory sample, not the full 800-cluster population.",
            "It does not include pixel/coadd covariance, PSF, sky transfer, extinction, or selection uncertainty.",
            "The strata and association were selected after exploratory analysis; no significance claim is made.",
        ],
    }
    summary_out = output_dir / "summary.json"
    summary_out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                           encoding="utf-8")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = analyze(args.input, args.output)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "groups": result["groups"], "output": str(args.output.resolve())},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
