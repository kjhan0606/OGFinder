#!/usr/bin/env python3
"""Audit frame influence and support changes in extended-PSF injection tests.

This compares per-scenario redetection with the saved fixed-baseline-mask
branch. It is a synthetic sensitivity audit, not a PSF correction or science
measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


METRICS = (
    "icl_response_to_diffuse_unit_injection",
    "icl_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def quantiles(values: list[float]) -> dict:
    a = np.asarray(values, dtype=float)
    if not a.size or not np.isfinite(a).all():
        raise ValueError("expected finite non-empty frame responses")
    return {"n": int(a.size), "median": float(np.median(a)),
            "q16": float(np.quantile(a, .16)), "q84": float(np.quantile(a, .84)),
            "min": float(np.min(a)), "max": float(np.max(a))}


def audit_batch(batch_dir: Path, expected_frames: int = 13) -> dict:
    batch_dir = Path(batch_dir).resolve()
    summary_path = batch_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete"
            or summary.get("science_measurement") is not False
            or summary.get("frame_filter") != "i"):
        raise ValueError("expected a complete exploratory i-band batch")
    records = summary.get("records", [])
    if len(records) != expected_frames:
        raise ValueError(f"expected {expected_frames} frame records, got {len(records)}")
    ids, rows = set(), []
    for record in records:
        ident = str(record["cluster_id"])
        if ident in ids or record.get("filter") != "i":
            raise ValueError("duplicate cluster or unexpected frame filter")
        ids.add(ident)
        path = Path(record["path"]).resolve()
        if not path.is_file() or sha256(path) != record.get("sha256"):
            raise ValueError(f"frame JSON missing or hash mismatch: {ident}")
        frame = json.loads(path.read_text(encoding="utf-8"))
        if (frame.get("status") != "complete" or frame.get("science_measurement") is not False
                or frame.get("filter") != "i"):
            raise ValueError(f"unexpected frame flags: {ident}")
        redetected = frame.get("transfer_response", {})
        fixed = frame.get("fixed_mask_transfer_response", {})
        delta = frame.get("redetection_minus_fixed_mask_response", {})
        for metric in METRICS:
            r, f, d = float(redetected[metric]), float(fixed[metric]), float(delta[metric])
            if not np.isclose(r - f, d, rtol=1e-10, atol=1e-10):
                raise ValueError(f"reselection arithmetic mismatch: {ident}/{metric}")
            if (not np.isclose(r, float(record["transfer_response"][metric]), rtol=1e-10, atol=1e-10)
                    or not np.isclose(f, float(record["fixed_mask_transfer_response"][metric]), rtol=1e-10, atol=1e-10)
                    or not np.isclose(d, float(record["redetection_minus_fixed_mask_response"][metric]), rtol=1e-10, atol=1e-10)):
                raise ValueError(f"summary/frame mismatch: {ident}/{metric}")
        measure = frame["measurements_by_scenario"]
        fixed_measure = frame["measurements_by_scenario_fixed_masks"]
        support = {name: measure[name]["satellite_support"] for name in ("baseline", "icl_only", "satellite_only", "both")}
        rows.append({
            "cluster_id": ident,
            "redetected": redetected,
            "fixed_mask": fixed,
            "reselection_delta": delta,
            "baseline_support": support["baseline"],
            "icl_only_support": support["icl_only"],
            "both_support": support["both"],
            "baseline_icl_coverage": float(measure["baseline"]["icl_coverage_200_300"]),
            "icl_only_icl_coverage": float(measure["icl_only"]["icl_coverage_200_300"]),
            "fixed_mask_icl_only_coverage": float(fixed_measure["icl_only"]["icl_coverage_200_300"]),
        })

    metrics = {}
    for metric in METRICS:
        red = [float(row["redetected"][metric]) for row in rows]
        fixed = [float(row["fixed_mask"][metric]) for row in rows]
        delta = [float(row["reselection_delta"][metric]) for row in rows]
        loo = []
        for index, row in enumerate(rows):
            red_without = red[:index] + red[index + 1:]
            fixed_without = fixed[:index] + fixed[index + 1:]
            loo.append({"cluster_id": row["cluster_id"],
                        "redetected_median_without": float(np.median(red_without)),
                        "redetected_median_shift": float(np.median(red_without) - np.median(red)),
                        "fixed_mask_median_without": float(np.median(fixed_without)),
                        "fixed_mask_median_shift": float(np.median(fixed_without) - np.median(fixed))})
        metrics[metric] = {
            "redetected": quantiles(red), "fixed_mask": quantiles(fixed),
            "reselection_delta": quantiles(delta),
            "largest_absolute_reselection_frames": sorted(
                ({"cluster_id": row["cluster_id"], "delta": float(row["reselection_delta"][metric])}
                 for row in rows), key=lambda item: abs(item["delta"]), reverse=True)[:3],
            "leave_one_out_median_influence": sorted(
                loo, key=lambda item: abs(item["redetected_median_shift"]), reverse=True),
        }
    rows.sort(key=lambda row: int(row["cluster_id"]))
    return {
        "schema": "astrafek-extended-psf-reselection-influence-audit-v1",
        "status": "validated_exploratory_sensitivity_only",
        "batch_dir": str(batch_dir), "summary_sha256": sha256(summary_path),
        "n_frames": len(rows), "metrics": metrics, "frames": rows,
        "science_measurement": False, "psf_approved": False,
        "limitations": [
            "Leave-one-frame-out shifts describe this support-selected 13-frame sample, not population uncertainty.",
            "This audit does not validate the public PSF, HSC sky transfer, source completeness, or coadd covariance.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", type=Path, required=True)
    parser.add_argument("--expected-frames", type=int, default=13)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite output: {output}")
    result = audit_batch(args.batch_dir, args.expected_frames)
    result["script_sha256"] = sha256(Path(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_frames": result["n_frames"],
                      "summary_sha256": result["summary_sha256"], "metrics": result["metrics"]},
                     indent=2))


if __name__ == "__main__":
    main()
