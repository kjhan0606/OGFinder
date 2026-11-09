#!/usr/bin/env python3
"""Joint block-bootstrap sensitivity for the 90-pair visit-split audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUT = ROOT / "icl_visit_split_covariance_sample_2026-11-06_02/paired_comparisons.tsv"
EXPECTED_SHA256 = "b41d428bf423692da0d849b8b590bf1c106fbfd0945c93cf5c3a837466638443"
DEFAULT_OUTPUT = ROOT / "icl_visit_split_covariance_joint_bootstrap_2026-11-06_01"
N_BOOT = 10000
METRICS = {
    "normalized_split_rms": ("rms_target", "rms_control"),
    "zero_lag_covariance": ("zero_lag_target", "zero_lag_control"),
    "cumulative_covariance_sum_8px": ("lag_sum_8_target", "lag_sum_8_control"),
    "clean_annulus_coverage": ("coverage_target", "coverage_control"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows() -> list[dict]:
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    digest = sha256(INPUT)
    if digest != EXPECTED_SHA256:
        raise ValueError(f"paired table hash mismatch: {digest}")
    with INPUT.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    rows = load_rows()
    report = {}
    for band_i, band in enumerate(("i", "i2")):
        required = [col for columns in METRICS.values() for col in columns]
        sample = [row for row in rows if row["filter"] == band
                  and all(row[col] not in ("", "None") for col in required)]
        blocks = sorted({row["block"] for row in sample})
        if not sample or len(blocks) < 2:
            report[band] = {"n_pairs_complete_case": len(sample), "n_blocks": len(blocks),
                            "metrics": {}, "max_t_critical_95": None}
            continue
        names = list(METRICS)
        delta = np.column_stack([
            np.asarray([float(row[t]) - float(row[c]) for row in sample], dtype=np.float64)
            for t, c in METRICS.values()
        ])
        observed = np.median(delta, axis=0)
        ix = {block: np.asarray([i for i, row in enumerate(sample) if row["block"] == block], dtype=int)
              for block in blocks}
        rng = np.random.default_rng(20261106 + band_i)
        boot = np.empty((N_BOOT, len(names)), dtype=np.float64)
        for b in range(N_BOOT):
            chosen = rng.choice(blocks, size=len(blocks), replace=True)
            index = np.concatenate([ix[x] for x in chosen])
            boot[b] = np.median(delta[index], axis=0)
        se = boot.std(axis=0, ddof=1)
        valid = se > 0
        max_t = np.max(np.abs((boot[:, valid] - observed[valid]) / se[valid]), axis=1)
        critical = float(np.quantile(max_t, 0.95))
        metric_results = {}
        for j, name in enumerate(names):
            metric_results[name] = {
                "median_target_minus_control": float(observed[j]),
                "bootstrap_se": float(se[j]),
                "pointwise_ci95": [float(x) for x in np.quantile(boot[:, j], [0.025, 0.975])],
                "simultaneous_ci95_within_filter": [float(observed[j] - critical * se[j]),
                                                    float(observed[j] + critical * se[j])],
            }
        report[band] = {
            "n_pairs_complete_case": len(sample), "n_connected_5deg_blocks": len(blocks),
            "n_replicates": N_BOOT, "max_t_critical_95": critical,
            "metrics": metric_results,
        }
    result = {
        "schema": "astrafek-icl-visit-split-joint-bootstrap-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": "complete-case paired median differences; resample connected 5-degree sky blocks; studentized max-|T| over four diagnostic metrics within each filter",
        "family_definition": "four metrics per filter, separately; filters are not combined into one family",
        "n_replicates": N_BOOT, "input_path": str(INPUT),
        "input_sha256": EXPECTED_SHA256, "results": report,
        "interpretation": [
            "The split-stack metrics remain estimator-conditional diagnostics, not calibrated ICL noise corrections.",
            "Simultaneous intervals account for four metrics within each filter but not between-filter multiplicity.",
            "The analysis uses the saved complete 90-pair visit-split products and makes no HSC request.",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "output": str(path), "results": report}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args().output_dir)
