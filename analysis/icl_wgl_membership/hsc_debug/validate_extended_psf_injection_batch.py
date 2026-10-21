#!/usr/bin/env python3
"""Independently audit extended-PSF injection outputs without rereading images."""
from __future__ import annotations

import argparse
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
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def recompute_transfer_response(record: dict, *, measurement_key: str = "measurements_by_scenario") -> dict:
    scenarios = record.get(measurement_key, {})
    required = {"baseline", "icl_only", "satellite_only", "both"}
    if set(scenarios) != required:
        raise ValueError("injection result must contain exactly four required scenarios")
    injection = record.get("injection", {})
    diffuse_amplitude = float(injection["diffuse_amplitude_at_100_njy_arcsec2"])
    satellite_amplitude = float(injection["satellite_amplitude_at_50_per_galaxy_njy_arcsec2"])
    if not (math.isfinite(diffuse_amplitude) and diffuse_amplitude > 0
            and math.isfinite(satellite_amplitude) and satellite_amplitude > 0):
        raise ValueError("injection amplitudes must be finite and positive")
    base_icl = float(scenarios["baseline"]["icl_raw_200_300_njy_arcsec2"])
    base_sat = float(scenarios["baseline"]["satellite_delta_50_75_njy_arcsec2"])
    icl_only_icl = float(scenarios["icl_only"]["icl_raw_200_300_njy_arcsec2"])
    icl_only_sat = float(scenarios["icl_only"]["satellite_delta_50_75_njy_arcsec2"])
    sat_only_icl = float(scenarios["satellite_only"]["icl_raw_200_300_njy_arcsec2"])
    sat_only_sat = float(scenarios["satellite_only"]["satellite_delta_50_75_njy_arcsec2"])
    both_icl = float(scenarios["both"]["icl_raw_200_300_njy_arcsec2"])
    both_sat = float(scenarios["both"]["satellite_delta_50_75_njy_arcsec2"])
    values = {
        "icl_response_to_diffuse_unit_injection": (icl_only_icl - base_icl) / diffuse_amplitude,
        "icl_response_to_all_satellite_envelopes": (sat_only_icl - base_icl) / satellite_amplitude,
        "satellite_response_to_smooth_diffuse": (icl_only_sat - base_sat) / diffuse_amplitude,
        "satellite_response_to_satellite_unit_injection": (sat_only_sat - base_sat) / satellite_amplitude,
        "both_minus_sum_residual_icl_njy_arcsec2":
            both_icl - base_icl - (icl_only_icl - base_icl) - (sat_only_icl - base_icl),
        "both_minus_sum_residual_satellite_njy_arcsec2":
            both_sat - base_sat - (icl_only_sat - base_sat) - (sat_only_sat - base_sat),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("injection transfer outputs must be finite")
    return values


def _validate_hybrid_kernel_metadata(metadata: dict, expected_parameters: dict,
                                    label: str) -> None:
    if (not isinstance(metadata, dict)
            or metadata.get("status") != "exploratory_core_wing_hybrid_sensitivity_only"
            or metadata.get("psf_approved") is not False
            or metadata.get("science_measurement") is not False):
        raise ValueError(f"{label} hybrid PSF metadata is missing or has unsafe science flags")
    for key in ("splice_radius_arcsec", "blend_width_arcsec", "wing_amplitude_pre_normalization"):
        expected_key = "wing_amplitude" if key == "wing_amplitude_pre_normalization" else key
        if not math.isclose(float(metadata.get(key, float("nan"))),
                            float(expected_parameters[expected_key]), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"{label} hybrid PSF parameter mismatch: {expected_key}")
    if not math.isclose(float(metadata.get("hybrid_signed_sum", float("nan"))),
                        1.0, rel_tol=0, abs_tol=2e-6):
        raise ValueError(f"{label} hybrid PSF is not signed-unit-normalized")
    if not math.isfinite(float(metadata.get("hybrid_signed_sum_before_normalization", float("nan")))):
        raise ValueError(f"{label} hybrid PSF pre-normalization sum is non-finite")


def validate_hybrid_metadata(record: dict, expected_parameters: dict) -> None:
    metadata = record.get("hybrid_psf")
    _validate_hybrid_kernel_metadata(metadata, expected_parameters, "center")
    injection = record.get("injection", {})
    if injection.get("hybrid_parameters") != expected_parameters:
        raise ValueError("frame injection metadata does not match the batch hybrid parameters")
    input_hashes = record.get("input_sha256", {})
    if not any("psf" in Path(path).name.lower() for path in input_hashes):
        raise ValueError("hybrid PSF frame record lacks a local PSF-core input checksum")
    support = record.get("support", {})
    expected_satellites = support.get("n_injected_targets_in_50_500_kpc")
    satellite_metadata = record.get("hybrid_satellite_psfs")
    if (not isinstance(expected_satellites, int)
            or not isinstance(satellite_metadata, list)
            or len(satellite_metadata) != expected_satellites):
        raise ValueError("hybrid PSF metadata must cover every injected satellite position")
    hashed_by_name = {Path(path).name: digest for path, digest in input_hashes.items()}
    for index, satellite in enumerate(satellite_metadata):
        if satellite.get("target_index") != index:
            raise ValueError("hybrid satellite PSF indices are not complete and ordered")
        _validate_hybrid_kernel_metadata(satellite, expected_parameters,
                                         f"satellite {index}")
        filename = satellite.get("local_core_filename")
        if not filename or filename not in hashed_by_name:
            raise ValueError(f"hybrid satellite {index} lacks a hashed local PSF core")
        if satellite.get("local_core_sha256") != hashed_by_name[filename]:
            raise ValueError(f"hybrid satellite {index} local PSF checksum does not match its input record")


def validate_point_source_profile(record: dict) -> None:
    profile = record.get("point_source_wing_profile")
    injection = record.get("injection", {})
    if (not isinstance(profile, dict)
            or injection.get("point_source_added_to_redetection_scenarios") is not False
            or not math.isclose(float(profile.get("source_flux_njy", float("nan"))),
                                float(injection.get("separate_unit_point_source_flux_njy", float("nan"))),
                                rel_tol=0, abs_tol=1e-12)):
        raise ValueError("separate point-source wing profile is missing or inconsistent")
    for value in (profile.get("pixel_scale_arcsec"), profile.get("kpc_per_pixel")):
        if not math.isfinite(float(value)) or float(value) <= 0:
            raise ValueError("point-source profile scales must be finite and positive")
    annuli = profile.get("annuli", {})
    if set(annuli) != {"50_100_kpc", "200_300_kpc"}:
        raise ValueError("point-source profile must include both physical annuli")
    for label, values in annuli.items():
        count = values.get("n_valid_pixels")
        expected = float(values.get("geometric_pixel_count", float("nan")))
        coverage = float(values.get("mask_coverage", float("nan")))
        fraction = float(values.get("masked_annulus_flux_fraction", float("nan")))
        mean = values.get("mean_surface_brightness_njy_arcsec2_per_njy")
        if (not isinstance(count, int) or count < 0 or not math.isfinite(expected)
                or expected <= 0 or not math.isfinite(coverage) or not 0 <= coverage <= 1
                or not math.isfinite(fraction)):
            raise ValueError(f"invalid point-source annulus metadata: {label}")
        if count == 0 and mean is not None:
            raise ValueError(f"empty point-source annulus must have a null mean: {label}")
        if count > 0 and (mean is None or not math.isfinite(float(mean))):
            raise ValueError(f"populated point-source annulus must have a finite mean: {label}")


def audit_batch(batch_dir: Path, *, expected_filter: str, expected_frames: int) -> dict:
    batch_dir = Path(batch_dir).resolve()
    summary_path = batch_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete"
            or summary.get("n_selected") != expected_frames
            or summary.get("n_completed") != expected_frames
            or summary.get("n_errors") != 0
            or summary.get("science_measurement") is not False
            or summary.get("observational_validation_passed") is not False
            or summary.get("flux_covariance_available") is not False):
        raise ValueError("batch is incomplete or has unexpected science flags")
    if summary.get("frame_filter") != expected_filter:
        raise ValueError("batch filter does not match the requested filter")
    proxy_mode = summary.get("public_reference_transfer_mode")
    expected_mode = ("exploratory_HSC_I_proxy_for_i2" if expected_filter == "i2"
                     else "HSC_I_reference_for_i_coadd")
    if proxy_mode != expected_mode:
        raise ValueError("unexpected public-reference transfer mode")
    reference = summary.get("reference", {})
    if reference.get("psf_approved") is not False:
        raise ValueError("unapproved extended PSF reference was expected")
    hybrid_parameters = summary.get("hybrid_parameters")
    hybrid_mode = summary.get("psf_model_mode") == "position_core_public_wing_hybrid_sensitivity"
    if hybrid_mode:
        required = {"splice_radius_arcsec", "blend_width_arcsec", "wing_amplitude"}
        if not isinstance(hybrid_parameters, dict) or set(hybrid_parameters) != required:
            raise ValueError("hybrid batch lacks its complete parameter record")
    elif hybrid_parameters is not None:
        raise ValueError("non-hybrid batch unexpectedly records hybrid parameters")

    records = summary.get("records", [])
    if len(records) != expected_frames:
        raise ValueError("summary record count is inconsistent")
    cluster_ids, responses = set(), {key: [] for key in METRICS}
    fixed_responses, reselection_effects = {key: [] for key in METRICS}, {key: [] for key in METRICS}
    max_arithmetic_error = 0.0
    for item in records:
        path = Path(item["path"]).resolve()
        if path.parent != batch_dir or sha256(path) != item.get("sha256"):
            raise ValueError("frame result is outside the batch or has a checksum mismatch")
        record = json.loads(path.read_text(encoding="utf-8"))
        ident = str(record.get("cluster_id"))
        if ident in cluster_ids:
            raise ValueError(f"duplicate frame result for cluster {ident}")
        cluster_ids.add(ident)
        if (record.get("schema") not in {"astrafek-extended-psf-joint-injection-v1",
                                          "astrafek-extended-psf-joint-injection-v2",
                                          "astrafek-extended-psf-joint-injection-v3"}
                or record.get("status") != "complete"
                or str(record.get("filter", "")).lower() != expected_filter
                or record.get("science_measurement") is not False
                or record.get("observational_validation_passed") is not False
                or record.get("flux_covariance_available") is not False):
            raise ValueError(f"frame {ident} has unexpected filter/status/science flags")
        if expected_filter == "i2" and record.get("injection", {}).get(
                "HSC_I_to_i2_is_proxy_sensitivity_only") is not True:
            raise ValueError(f"i2 frame {ident} lacks the proxy-sensitivity label")
        if hybrid_mode:
            validate_hybrid_metadata(record, hybrid_parameters)
        if record.get("schema") == "astrafek-extended-psf-joint-injection-v3":
            validate_point_source_profile(record)
        calculated = recompute_transfer_response(record)
        stored = record.get("transfer_response", {})
        for key, value in calculated.items():
            saved = float(stored[key])
            if not math.isclose(value, saved, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"frame {ident} stored response does not reproduce: {key}")
            if key in METRICS:
                responses[key].append(value)
                max_arithmetic_error = max(max_arithmetic_error, abs(value - saved))
        if record.get("schema") in {"astrafek-extended-psf-joint-injection-v2",
                                      "astrafek-extended-psf-joint-injection-v3"}:
            fixed = recompute_transfer_response(
                record, measurement_key="measurements_by_scenario_fixed_masks")
            saved_fixed = record.get("fixed_mask_transfer_response", {})
            saved_delta = record.get("redetection_minus_fixed_mask_response", {})
            for key, value in fixed.items():
                if key not in saved_fixed or not math.isclose(
                        value, float(saved_fixed[key]), rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"frame {ident} fixed-mask response does not reproduce: {key}")
                delta = calculated[key] - value
                if key not in saved_delta or not math.isclose(
                        delta, float(saved_delta[key]), rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"frame {ident} mask-selection effect does not reproduce: {key}")
                if key in METRICS:
                    fixed_responses[key].append(value)
                    reselection_effects[key].append(delta)

    quantiles = {
        key: {
            "n_frames": len(values),
            "median": float(np.median(values)),
            "q16": float(np.quantile(values, 0.16)),
            "q84": float(np.quantile(values, 0.84)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        }
        for key, values in responses.items()
    }
    fixed_quantiles = {
        key: {"n_frames": len(values), "median": float(np.median(values)),
              "q16": float(np.quantile(values, 0.16)), "q84": float(np.quantile(values, 0.84)),
              "min": float(np.min(values)), "max": float(np.max(values))}
        for key, values in fixed_responses.items() if values
    }
    reselection_quantiles = {
        key: {"n_frames": len(values), "median": float(np.median(values)),
              "q16": float(np.quantile(values, 0.16)), "q84": float(np.quantile(values, 0.84)),
              "min": float(np.min(values)), "max": float(np.max(values))}
        for key, values in reselection_effects.items() if values
    }
    return {
        "schema": "astrafek-extended-psf-injection-batch-validation-v1",
        "status": "validated_exploratory_transfer_only",
        "batch_dir": str(batch_dir),
        "summary_sha256": sha256(summary_path),
        "expected_filter": expected_filter,
        "n_frames_validated": len(cluster_ids),
        "cluster_ids": sorted(cluster_ids, key=int),
        "max_stored_vs_recomputed_response_error": max_arithmetic_error,
        "frame_response_summaries": quantiles,
        "fixed_mask_frame_response_summaries": fixed_quantiles,
        "redetection_minus_fixed_mask_summaries": reselection_quantiles,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "limitations": [
            "Validator checks saved JSON provenance and response arithmetic; it does not reread source FITS images.",
            "For i2, HSC-I is an exploratory proxy sensitivity, not a filter-specific PSF calibration.",
            "Frame quantiles describe selected-frame spread, not a confidence interval or population estimate.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", type=Path, required=True)
    parser.add_argument("--filter", choices=("i", "i2"), required=True)
    parser.add_argument("--expected-frames", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing validation: {output}")
    result = audit_batch(args.batch_dir, expected_filter=args.filter,
                         expected_frames=args.expected_frames)
    result["validator_sha256"] = sha256(Path(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "expected_filter", "n_frames_validated", "max_stored_vs_recomputed_response_error",
    )}, indent=2))


if __name__ == "__main__":
    main()
