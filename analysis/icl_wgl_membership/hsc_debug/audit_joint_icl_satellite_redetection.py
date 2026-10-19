#!/usr/bin/env python3
"""Verify a joint ICL/satellite redetection-injection batch and its gate state."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(summary_path: Path) -> dict:
    summary_path = Path(summary_path).resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "hsc-joint-icl-satellite-fixed-mask-batch-v1":
        raise ValueError("unexpected joint-injection batch schema")
    if summary.get("science_measurement") is not False:
        raise ValueError("batch cannot be promoted to a science measurement")
    injection = summary.get("injection", {})
    if (injection.get("detection_rerun") is not True
            or injection.get("masks_and_segmentation_fixed") is not False
            or injection.get("poisson_noise_added") is not False
            or injection.get("psf_convolution_applied") is not False
            or injection.get("original_fits_modified") is not False):
        raise ValueError("batch provenance does not match the redetection injection stage")
    selection = summary.get("selection", {})
    selected = [(str(row["cluster_id"]), str(row["filter"]))
                for row in selection.get("selected_pairs", [])]
    if len(selected) != selection.get("n_selected_cluster_filter_pairs") or len(set(selected)) != len(selected):
        raise ValueError("selected pair list is missing or contains duplicates")
    records = summary.get("records", [])
    errors = summary.get("errors", [])
    record_keys = [(str(row["cluster_id"]), str(row["filter"])) for row in records]
    error_keys = [(str(row["cluster_id"]), str(row["filter"])) for row in errors]
    if len(set(record_keys)) != len(record_keys) or len(set(error_keys)) != len(error_keys):
        raise ValueError("batch records contain duplicate pair keys")
    if set(record_keys) & set(error_keys) or set(record_keys) | set(error_keys) != set(selected):
        raise ValueError("selected pairs are not partitioned exactly into results and errors")
    if (len(records) != summary.get("n_completed")
            or len(errors) != summary.get("n_errors")
            or summary.get("status") != ("complete" if not errors else "partial")):
        raise ValueError("batch totals/status do not match its records")

    per_filter = Counter()
    for record in records:
        path = Path(record["path"])
        if not path.is_file() or sha256(path) != record.get("sha256"):
            raise ValueError(f"missing or changed frame result: {path}")
        frame = json.loads(path.read_text(encoding="utf-8"))
        if (str(frame.get("cluster_id")) != str(record["cluster_id"])
                or frame.get("filter") != record["filter"]
                or frame.get("schema") != "hsc-joint-icl-satellite-redetection-injection-v1"
                or frame.get("science_measurement") is not False
                or frame.get("observational_validation_passed") is not False
                or frame.get("flux_covariance_available") is not False
                or frame.get("injection", {}).get("source_detection_and_masks_rerun") is not True
                or frame.get("injection", {}).get("original_fits_modified") is not False):
            raise ValueError(f"frame identity or science gate mismatch: {path}")
        if set(frame.get("measurements_by_scenario", {})) != {"baseline", "icl_only", "satellite_only", "both"}:
            raise ValueError(f"frame does not contain all paired injection scenarios: {path}")
        if set(frame.get("redetection_by_scenario", {})) != {"baseline", "icl_only", "satellite_only", "both"}:
            raise ValueError(f"frame is missing redetection diagnostics: {path}")
        per_filter[str(record["filter"])] += 1

    return {
        "schema": "astrafek-joint-icl-satellite-redetection-audit-v1",
        "status": "verified_complete" if not errors else "verified_partial",
        "science_measurement": False,
        "n_selected": len(selected),
        "n_completed": len(records),
        "n_errors": len(errors),
        "completed_by_filter": dict(sorted(per_filter.items())),
        "errors_by_type_and_message": {
            f"{kind}: {message}": count
            for (kind, message), count in sorted(Counter(
                (str(row.get("error_type")), str(row.get("error"))) for row in errors
            ).items())
        },
        "summary_path": str(summary_path),
        "summary_sha256": sha256(summary_path),
        "frame_sha256_verified": len(records),
        "all_frames_keep_science_measurement_false": True,
        "limitations": [
            "Redetection success is not equivalent to scientific validation; failures are not imputed.",
            "The injection stage does not add injected-source photon noise or convolve templates with measured PSFs.",
            "The delivered HSC upstream sky processing and image-level covariance are not modeled.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = audit(args.summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("status", "n_selected", "n_completed", "n_errors", "completed_by_filter",
                       "errors_by_type_and_message", "summary_sha256")}, indent=2))


if __name__ == "__main__":
    main()
