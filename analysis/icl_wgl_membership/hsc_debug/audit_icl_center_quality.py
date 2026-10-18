#!/usr/bin/env python3
"""Test whether outer-annulus ICL excess depends on center quality flags.

This is a post-selection descriptive sensitivity analysis of the existing
K2.5/40%, 200--300 kpc pair measurements. It does not certify the automated
center flags, correct selection, or change the science-measurement gate.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def center_eligible(row):
    return row["auto_unique_center_candidate"].strip().lower() == "true" and \
        row["center_masked_by_kron"].strip().lower() == "false"


def summarize(rows, replicas=10000, seed=20261018):
    if not rows:
        return {"status": "insufficient_support", "n_clusters": 0,
                "mean_delta_flux": None, "bootstrap_ci95": [None, None],
                "n_joint_blocks": 0}
    target, control = [], []
    for row in rows:
        target.append({"ra_deg": float(row["target_ra_deg"]),
                       "dec_deg": float(row["target_dec_deg"])})
        control.append({"ra_deg": float(row["control_ra_deg"]),
                        "dec_deg": float(row["control_dec_deg"])})
    selected = np.ones(len(rows), dtype=bool)
    indices, labels, _, _, n_blocks, largest = pair_block_labels(target, control, selected)
    values = np.asarray([float(row["delta_flux"]) for row in rows], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Non-finite outer-annulus ICL value")
    unique = np.unique(labels)
    members = {block: np.flatnonzero(labels == block) for block in unique}
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    for j in range(replicas):
        chosen = rng.integers(0, len(unique), size=len(unique))
        sample = np.concatenate([members[unique[k]] for k in chosen])
        draws[j] = float(values[sample].mean())
    full_mean = float(values.mean())
    leave_one_block_out = []
    for block in unique:
        kept = labels != block
        if kept.any():
            leave_one_block_out.append(float(values[kept].mean()))
    return {
        "status": "descriptive", "n_clusters": len(rows), "n_joint_blocks": n_blocks,
        "largest_joint_block": largest, "mean_delta_flux": full_mean,
        "median_delta_flux": float(np.median(values)),
        "spatial_block_bootstrap_ci95": [float(x) for x in np.quantile(draws, [.025, .975])],
        "leave_one_block_out_mean_range": [min(leave_one_block_out), max(leave_one_block_out)],
        "leave_one_block_out_max_abs_shift": max(abs(x - full_mean) for x in leave_one_block_out),
        "replicas": replicas, "seed": seed,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Refusing to overwrite existing output: {output}")
    if args.replicas < 1000:
        parser.error("At least 1000 bootstrap replicas are required")

    all_rows = read_tsv(args.pairs)
    rows = [row for row in all_rows
            if float(row["kron_factor"]) == 2.5
            and float(row["coverage_floor"]) == .40
            and row["radial_bin_kpc"] == "200-300"
            and float(row["common_coverage"]) >= .40]
    ids = [row["id"] for row in rows]
    if len(rows) != 224 or len(set(ids)) != len(rows):
        parser.error(f"Expected 224 unique K2.5/40% outer-annulus pairs, got {len(rows)}")

    result = {
        "schema": "hsc-icl-center-quality-sensitivity-v1",
        "status": "complete_descriptive", "science_measurement": False,
        "selection_correction_applied": False,
        "definition": "center-eligible iff unique automatic center candidate and K2.5 center unmasked",
        "input_sha256": {str(args.pairs.resolve()): sha256(args.pairs)},
        "analysis_script_sha256": sha256(Path(__file__)),
        "bootstrap": {"method": "resample connected 5-degree target/control sky blocks",
                      "replicas": args.replicas, "seed": args.seed},
        "by_filter": {},
        "limitations": [
            "Center flags are automatic diagnostics, not expert BCG validation.",
            "Groups are selected from the same images and the split is post hoc.",
            "Intervals do not include PSF wings, sky-transfer, coadd pixel covariance, or mask uncertainty.",
            "No significance or physical ICL interpretation is claimed; science_measurement remains false.",
        ],
    }
    for band_i, band in enumerate(("i", "i2")):
        band_rows = [row for row in rows if row["filter"] == band]
        eligible = [row for row in band_rows if center_eligible(row)]
        ineligible = [row for row in band_rows if not center_eligible(row)]
        result["by_filter"][band] = {
            "all": summarize(band_rows, args.replicas, args.seed + band_i * 10),
            "center_eligible": summarize(eligible, args.replicas, args.seed + band_i * 10 + 1),
            "not_center_eligible": summarize(ineligible, args.replicas, args.seed + band_i * 10 + 2),
        }

    output.mkdir(parents=True, exist_ok=False)
    values_path = output / "center_quality_summary.tsv"
    with values_path.open("x", encoding="utf-8", newline="") as stream:
        fields = ["filter", "group", "n_clusters", "n_joint_blocks", "largest_joint_block",
                  "mean_delta_flux", "median_delta_flux", "ci95_low", "ci95_high",
                  "loo_mean_min", "loo_mean_max", "loo_max_abs_shift"]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for band, groups in result["by_filter"].items():
            for group, summary in groups.items():
                ci = summary.get("spatial_block_bootstrap_ci95", [None, None])
                writer.writerow({"filter": band, "group": group,
                    "n_clusters": summary["n_clusters"],
                    "n_joint_blocks": summary.get("n_joint_blocks", 0),
                    "largest_joint_block": summary.get("largest_joint_block", 0),
                    "mean_delta_flux": summary.get("mean_delta_flux"),
                    "median_delta_flux": summary.get("median_delta_flux"),
                    "ci95_low": ci[0], "ci95_high": ci[1],
                    "loo_mean_min": (summary.get("leave_one_block_out_mean_range") or [None, None])[0],
                    "loo_mean_max": (summary.get("leave_one_block_out_mean_range") or [None, None])[1],
                    "loo_max_abs_shift": summary.get("leave_one_block_out_max_abs_shift")})
    result["output_sha256"] = {values_path.name: sha256(values_path)}
    (output / "summary.json").write_text(json.dumps(result, ensure_ascii=False,
        indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_pairs": len(rows),
        "by_filter": result["by_filter"], "output": str(output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
