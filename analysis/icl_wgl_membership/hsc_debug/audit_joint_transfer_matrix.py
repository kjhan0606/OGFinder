#!/usr/bin/env python3
"""Audit whether fixed-mask joint-injection responses support de-mixing.

This is a conditioning and sample-overlap audit only. It deliberately does not
invert the transfer matrices or claim an astrophysical correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


METRICS = (
    "icl_response_to_diffuse_unit_injection",
    "icl_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def transfer_matrix(response: dict) -> np.ndarray:
    values = np.asarray([response[key] for key in METRICS], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("transfer response must be finite")
    return np.asarray([[values[0], values[1]], [values[2], values[3]]])


def matrix_diagnostics(response: dict) -> dict:
    matrix = transfer_matrix(response)
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    condition = (float(singular_values[0] / singular_values[-1])
                 if singular_values[-1] > 0 else math.inf)
    return {
        "matrix": matrix.tolist(),
        "determinant": float(np.linalg.det(matrix)),
        "condition_number_2": condition,
        "singular_values": singular_values.tolist(),
    }


def summarize_filter(rows: list[dict]) -> dict:
    conditions = np.asarray([r["condition_number_2"] for r in rows], dtype=float)
    determinants = np.asarray([r["determinant"] for r in rows], dtype=float)
    cross_to_diffuse = np.asarray([
        abs(r["matrix"][0][1]) / max(abs(r["matrix"][0][0]), 1e-12)
        for r in rows
    ])
    return {
        "n_pairs": len(rows),
        "n_negative_determinant": int(np.sum(determinants < 0)),
        "n_condition_gt_10": int(np.sum(conditions > 10)),
        "n_condition_gt_100": int(np.sum(conditions > 100)),
        "n_condition_gt_1000": int(np.sum(conditions > 1000)),
        "condition_q50_q90_max": np.quantile(conditions, [0.5, 0.9, 1.0]).tolist(),
        "median_icl_satellite_cross_to_diffuse_response": float(np.median(cross_to_diffuse)),
    }


def audit(batch_summary_path: Path, spatial_audit_path: Path,
          association_pairs_path: Path) -> dict:
    batch_summary = json.loads(batch_summary_path.read_text(encoding="utf-8"))
    spatial = json.loads(spatial_audit_path.read_text(encoding="utf-8"))
    if batch_summary.get("status") != "complete" or batch_summary.get("n_errors") != 0:
        raise ValueError("joint injection batch is not complete")
    if spatial.get("status") != "verified_exploratory":
        raise ValueError("spatial batch audit is not verified")
    if spatial.get("batch_summary_sha256") != sha256(batch_summary_path):
        raise ValueError("spatial audit references a different batch summary")
    if len(batch_summary.get("records", [])) != spatial.get("n_clusters"):
        raise ValueError("batch and spatial audit pair counts differ")

    records = []
    for item in batch_summary["records"]:
        path = Path(item["path"])
        if not path.is_file() or sha256(path) != item.get("sha256"):
            raise ValueError(f"frame result hash mismatch: {path}")
        frame = json.loads(path.read_text(encoding="utf-8"))
        if (str(frame.get("cluster_id")) != str(item["cluster_id"])
                or frame.get("filter") != item["filter"]
                or frame.get("script_sha256") != batch_summary.get("single_frame_script_sha256")
                or frame.get("science_measurement") is not False
                or frame.get("observational_validation_passed") is not False):
            raise ValueError(f"frame identity/hash/gate mismatch: {path}")
        diag = matrix_diagnostics(frame["transfer_response"])
        records.append({"cluster_id": str(item["cluster_id"]), "filter": item["filter"], **diag})

    association_rows = []
    with association_pairs_path.open(encoding="utf-8", newline="") as stream:
        association_rows = list(csv.DictReader(stream, delimiter="\t"))
    association_ids = [(str(row["cluster_id"]), row["filter"]) for row in association_rows]
    if len(set(association_ids)) != len(association_ids):
        raise ValueError("association pairs contain duplicate cluster/filter rows")
    injected_ids = {(row["cluster_id"], row["filter"]) for row in records}
    overlap = injected_ids & set(association_ids)
    by_filter = {}
    for filt in sorted({r["filter"] for r in records}):
        subset = [r for r in records if r["filter"] == filt]
        assoc_total = sum(row["filter"] == filt for row in association_rows)
        overlap_count = sum((row["cluster_id"], row["filter"]) in overlap
                            for row in association_rows if row["filter"] == filt)
        by_filter[filt] = {
            **summarize_filter(subset),
            "association_pairs": assoc_total,
            "association_pairs_with_injection": overlap_count,
        }
    if len(records) != spatial.get("n_clusters"):
        raise ValueError("record count does not match spatial audit")
    return {
        "schema": "hsc-joint-transfer-matrix-conditioning-audit-v1",
        "status": "verified_exploratory",
        "science_measurement": False,
        "deconvolution_performed": False,
        "reason_no_deconvolution": (
            "Some empirical response matrices have negative determinants or large condition numbers; "
            "the fixed-mask injection model also omits PSF, sky processing, detection/mask response, "
            "and image-level measurement covariance."
        ),
        "n_injection_pairs": len(records),
        "n_association_pairs": len(association_rows),
        "n_pairs_in_common": len(overlap),
        "input_sha256": {
            str(batch_summary_path): sha256(batch_summary_path),
            str(spatial_audit_path): sha256(spatial_audit_path),
            str(association_pairs_path): sha256(association_pairs_path),
        },
        "by_filter": by_filter,
        "records": records,
        "limitations": [
            "The 2x2 matrices characterize assumed fixed-mask templates, not calibrated population-level response.",
            "The association table and injection support selection overlap only partially.",
            "Condition numbers describe numerical inversion sensitivity and do not include measurement noise.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-summary", type=Path, required=True)
    parser.add_argument("--spatial-audit", type=Path, required=True)
    parser.add_argument("--association-pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = audit(args.batch_summary.resolve(), args.spatial_audit.resolve(),
                   args.association_pairs.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_injection_pairs": result["n_injection_pairs"],
                      "n_association_pairs": result["n_association_pairs"],
                      "n_pairs_in_common": result["n_pairs_in_common"],
                      "by_filter": result["by_filter"]}, indent=2))


if __name__ == "__main__":
    main()
