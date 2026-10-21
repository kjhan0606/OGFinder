#!/usr/bin/env python3
"""Independently verify and summarize a completed hybrid PSF grid."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from validate_extended_psf_injection_batch import audit_batch


METRICS = (
    "icl_response_to_diffuse_unit_injection",
    "icl_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize_grid(grid_root: Path) -> dict:
    grid_root = Path(grid_root).resolve()
    checkpoint_path = grid_root / "grid_status.json"
    state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    spec = state.get("spec", {})
    expected_count = int(spec.get("n_configurations", 0))
    if (state.get("status") != "complete" or expected_count < 1
            or len(state.get("records", [])) != expected_count
            or spec.get("science_measurement") is not False
            or spec.get("psf_approved") is not False):
        raise ValueError("grid is incomplete or has unexpected approval flags")

    rows, values = [], {metric: [] for metric in METRICS}
    seen = set()
    for item in state["records"]:
        key = str(item["configuration"])
        if key in seen or item.get("status") != "complete":
            raise ValueError("grid checkpoint has duplicate or incomplete configurations")
        seen.add(key)
        batch_dir = grid_root / key
        summary_path = batch_dir / "summary.json"
        if sha256(summary_path) != item.get("summary_sha256"):
            raise ValueError(f"grid summary checksum mismatch for {key}")
        batch_summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if batch_summary.get("hybrid_parameters") != item.get("parameters"):
            raise ValueError(f"grid parameter mismatch for {key}")
        validation = audit_batch(
            batch_dir, expected_filter=str(spec["filter"]), expected_frames=1)
        if (validation.get("n_frames_validated") != 1
                or validation.get("max_stored_vs_recomputed_response_error") != 0.0):
            raise ValueError(f"grid result arithmetic failed independent validation: {key}")
        record = json.loads(Path(batch_summary["records"][0]["path"]).read_text(
            encoding="utf-8"))
        response = record["transfer_response"]
        fixed = record["fixed_mask_transfer_response"]
        delta = record["redetection_minus_fixed_mask_response"]
        metrics = {}
        for metric in METRICS:
            triple = {"redetection": float(response[metric]),
                      "fixed_masks": float(fixed[metric]),
                      "redetection_minus_fixed_masks": float(delta[metric])}
            if not all(math.isfinite(value) for value in triple.values()):
                raise ValueError(f"non-finite grid response in {key}/{metric}")
            metrics[metric] = triple
            values[metric].append(triple)
        rows.append({
            "configuration": key,
            "parameters": item["parameters"],
            "frame_sha256": batch_summary["records"][0]["sha256"],
            "batch_summary_sha256": sha256(summary_path),
            "validated_frames": validation["n_frames_validated"],
            "max_arithmetic_error": validation[
                "max_stored_vs_recomputed_response_error"],
            "metrics": metrics,
        })
    if len(seen) != expected_count:
        raise ValueError("grid configuration count differs from its pinned specification")

    envelopes = {}
    for metric, entries in values.items():
        envelopes[metric] = {}
        for branch in ("redetection", "fixed_masks", "redetection_minus_fixed_masks"):
            branch_values = [entry[branch] for entry in entries]
            envelopes[metric][branch] = {
                "minimum_across_models": min(branch_values),
                "maximum_across_models": max(branch_values),
            }
    return {
        "schema": "astrafek-core-wing-hybrid-grid-audit-v1",
        "status": "validated_exploratory_sensitivity_only",
        "grid_spec": spec,
        "grid_checkpoint_sha256": sha256(checkpoint_path),
        "n_configurations_validated": len(rows),
        "configurations": rows,
        "model_sensitivity_envelopes": envelopes,
        "science_measurement": False,
        "psf_approved": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "limitations": [
            "This is a one-cluster pilot; model envelopes are not confidence intervals.",
            "The public HSC-I wing is not matched to the frame position, epoch, or seeing.",
            "No isolated point-source injection, source Poisson noise, or coadd covariance is included.",
            "The run tests synthetic-light transfer only and does not alter observed science data.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = summarize_grid(args.grid_root)
    result["audit_script_sha256"] = sha256(Path(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True,
                                 allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "configurations": result["n_configurations_validated"],
                      "science_measurement": False,
                      "output": str(output),
                      "output_sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
