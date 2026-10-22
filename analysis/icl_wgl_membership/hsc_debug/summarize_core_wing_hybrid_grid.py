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


def finite_summary(values: list[float]) -> dict:
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("summary requires a non-empty finite sample")
    ordered = sorted(values)
    mid = len(ordered) // 2
    median = (ordered[mid] if len(ordered) % 2
              else 0.5 * (ordered[mid - 1] + ordered[mid]))
    # Quantiles are deliberately not used as inferential intervals here; the
    # grid is a set of model alternatives and the sample may be small.
    import numpy as np
    return {"n_frames": len(values), "median": float(median),
            "q16": float(np.quantile(values, 0.16)),
            "q84": float(np.quantile(values, 0.84)),
            "minimum": float(min(values)), "maximum": float(max(values))}


def optional_summary(values: list[float]) -> dict:
    if values:
        return finite_summary(values)
    return {"n_frames": 0, "median": None, "q16": None, "q84": None,
            "minimum": None, "maximum": None}


def validate_frame_set(frame_ids: list[str], expected_frames: int,
                       reference_frame_ids: set[str] | None = None) -> set[str]:
    frame_id_set = set(frame_ids)
    if len(frame_id_set) != expected_frames:
        raise ValueError("grid configuration contains duplicate frame IDs")
    if reference_frame_ids is not None and frame_id_set != reference_frame_ids:
        raise ValueError("grid configurations do not use identical frame sets")
    return frame_id_set


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

    expected_frames = int(spec.get("n_frames", 0))
    if expected_frames < 1:
        raise ValueError("grid specification has no expected frame count")
    rows, values = [], {metric: {branch: [] for branch in
                                 ("redetection", "fixed_masks",
                                  "redetection_minus_fixed_masks")}
                        for metric in METRICS}
    point_values = {annulus: {name: [] for name in
                    ("mask_coverage", "masked_annulus_flux_fraction",
                     "mean_surface_brightness_njy_arcsec2_per_njy")}
                    for annulus in ("50_100_kpc", "200_300_kpc")}
    seen = set()
    expected_frame_ids = None
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
            batch_dir, expected_filter=str(spec["filter"]),
            expected_frames=expected_frames)
        if (validation.get("n_frames_validated") != expected_frames
                or validation.get("max_stored_vs_recomputed_response_error") != 0.0):
            raise ValueError(f"grid result arithmetic failed independent validation: {key}")
        if len(batch_summary.get("records", [])) != expected_frames:
            raise ValueError(f"frame count mismatch in {key}")
        frame_rows = []
        per_configuration_values = {metric: {branch: [] for branch in
                                    ("redetection", "fixed_masks",
                                     "redetection_minus_fixed_masks")}
                                    for metric in METRICS}
        per_configuration_points = {annulus: {name: [] for name in point_values[annulus]}
                                    for annulus in point_values}
        for frame_item in batch_summary["records"]:
            record_path = Path(frame_item["path"]).resolve()
            if sha256(record_path) != frame_item.get("sha256"):
                raise ValueError(f"frame checksum mismatch in {key}: {record_path}")
            record = json.loads(record_path.read_text(encoding="utf-8"))
            response = record["transfer_response"]
            fixed = record["fixed_mask_transfer_response"]
            delta = record["redetection_minus_fixed_mask_response"]
            frame_metrics = {}
            for metric in METRICS:
                triple = {"redetection": float(response[metric]),
                          "fixed_masks": float(fixed[metric]),
                          "redetection_minus_fixed_masks": float(delta[metric])}
                if not all(math.isfinite(value) for value in triple.values()):
                    raise ValueError(f"non-finite grid response in {key}/{metric}")
                frame_metrics[metric] = triple
                for branch, value in triple.items():
                    per_configuration_values[metric][branch].append(value)
                    values[metric][branch].append(value)

            point_profile = record.get("point_source_wing_profile")
            if not isinstance(point_profile, dict):
                raise ValueError(f"separate point-source mask profile is missing in {key}")
            frame_points = {}
            for annulus, annulus_values in point_profile.get("annuli", {}).items():
                if annulus not in point_values:
                    raise ValueError(f"unexpected point-source annulus: {annulus}")
                frame_points[annulus] = {}
                for name in point_values[annulus]:
                    raw = annulus_values.get(name)
                    if raw is None:
                        if (name == "mean_surface_brightness_njy_arcsec2_per_njy"
                                and int(annulus_values.get("n_valid_pixels", -1)) == 0):
                            frame_points[annulus][name] = None
                            continue
                        raise ValueError(f"missing point-source metric {annulus}/{name}")
                    if not math.isfinite(float(raw)):
                        raise ValueError(f"missing or non-finite point-source metric {annulus}/{name}")
                    value = float(raw)
                    frame_points[annulus][name] = value
                    per_configuration_points[annulus][name].append(value)
                    point_values[annulus][name].append(value)
            if set(frame_points) != set(point_values):
                raise ValueError(f"point-source annulus set is incomplete in {key}")
            frame_rows.append({"cluster_id": str(record["cluster_id"]),
                               "frame_sha256": frame_item["sha256"],
                               "metrics": frame_metrics,
                               "point_source_annuli": frame_points})

        metrics = {metric: {branch: finite_summary(branch_values)
                            for branch, branch_values in branches.items()}
                   for metric, branches in per_configuration_values.items()}
        point_summaries = {annulus: {name: optional_summary(metric_values)
                                      for name, metric_values
                                      in metric_map.items()}
                           for annulus, metric_map
                           in per_configuration_points.items()}
        frame_ids = [row["cluster_id"] for row in frame_rows]
        frame_id_set = validate_frame_set(frame_ids, expected_frames,
                                          expected_frame_ids)
        if expected_frame_ids is None:
            expected_frame_ids = frame_id_set
        rows.append({
            "configuration": key,
            "parameters": item["parameters"],
            "frame_ids": frame_ids,
            "frames": frame_rows,
            "batch_summary_sha256": sha256(summary_path),
            "validated_frames": validation["n_frames_validated"],
            "max_arithmetic_error": validation[
                "max_stored_vs_recomputed_response_error"],
            "metrics": metrics,
            "point_source_annuli": point_summaries,
        })
    if len(seen) != expected_count:
        raise ValueError("grid configuration count differs from its pinned specification")

    envelopes = {}
    for metric, branches in values.items():
        envelopes[metric] = {}
        for branch, frame_values in branches.items():
            model_medians = [row["metrics"][metric][branch]["median"]
                             for row in rows]
            envelopes[metric][branch] = {
                "minimum_model_median": min(model_medians),
                "maximum_model_median": max(model_medians),
                "minimum_frame_value": min(frame_values),
                "maximum_frame_value": max(frame_values),
            }
    point_envelopes = {}
    for annulus, metrics in point_values.items():
        point_envelopes[annulus] = {}
        for metric, frame_values in metrics.items():
            model_medians = [row["point_source_annuli"][annulus][metric]["median"]
                             for row in rows
                             if row["point_source_annuli"][annulus][metric]["median"] is not None]
            if not model_medians:
                point_envelopes[annulus][metric] = {
                    "minimum_model_median": None,
                    "maximum_model_median": None,
                    "minimum_frame_value": None,
                    "maximum_frame_value": None,
                }
                continue
            point_envelopes[annulus][metric] = {
                "minimum_model_median": min(model_medians),
                "maximum_model_median": max(model_medians),
                "minimum_frame_value": min(frame_values),
                "maximum_frame_value": max(frame_values),
            }
    return {
        "schema": "astrafek-core-wing-hybrid-grid-audit-v2",
        "status": "validated_exploratory_sensitivity_only",
        "grid_spec": spec,
        "grid_checkpoint_sha256": sha256(checkpoint_path),
        "n_configurations_validated": len(rows),
        "configurations": rows,
        "model_sensitivity_envelopes": envelopes,
        "point_source_model_sensitivity_envelopes": point_envelopes,
        "science_measurement": False,
        "psf_approved": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "limitations": [
            "The frame sample and model envelopes are sensitivity descriptions, not confidence intervals.",
            "The public HSC-I wing is not matched to the frame position, epoch, or seeing.",
            "The separate point-source profile is noiseless and measured through saved masks; it is not a source-detection recovery test.",
            "Source Poisson noise, coadd covariance, realistic sky transfer, and observed-source subtraction are not modeled.",
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
