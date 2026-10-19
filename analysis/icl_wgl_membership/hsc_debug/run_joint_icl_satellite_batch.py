#!/usr/bin/env python3
"""Run fixed-mask joint injection transfer tests for matched ICL/satellite IDs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import time
from pathlib import Path

import numpy as np

from joint_icl_satellite_injection import (
    DEFAULT_DATA,
    DEFAULT_ICL_MASKS,
    DEFAULT_ICL_PAIRS,
    DEFAULT_SAT_PREP,
    run_frame,
    sha256,
)
from position_psf import PositionPSFBank


DEFAULT_SATELLITE_TABLE = DEFAULT_DATA / "satellite_outer_profile_perbin_full_2026-10-18_01/satellite_profiles.tsv"
SCRIPT_SHA256_AT_IMPORT = sha256(Path(__file__))
SINGLE_FRAME_SCRIPT = Path(__file__).with_name("joint_icl_satellite_injection.py")
SINGLE_FRAME_SHA256_AT_IMPORT = sha256(SINGLE_FRAME_SCRIPT)
POSITION_PSF_SCRIPT = Path(__file__).with_name("position_psf.py")
POSITION_PSF_SHA256_AT_IMPORT = sha256(POSITION_PSF_SCRIPT)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def select_jointly_supported_pairs(icl_pairs_path: Path,
                                   satellite_profiles_path: Path,
                                   icl_coverage_floor: float = 0.40,
                                   satellite_coverage_floor: float = 0.60) -> list[tuple[str, str]]:
    """Select unique cluster/filter pairs with both annuli in saved support."""
    icl_rows = read_tsv(icl_pairs_path)
    sat_rows = read_tsv(satellite_profiles_path)
    if not icl_rows or not {"id", "filter", "coverage_200_300"}.issubset(icl_rows[0]):
        raise ValueError("ICL pair table is missing the expected common-coverage fields")
    required_sat = {"cluster_id", "filter", "radial_bin_kpc", "target_coverage", "control_coverage"}
    if not sat_rows or not required_sat.issubset(sat_rows[0]):
        raise ValueError("satellite table is missing the expected per-annulus support fields")
    icl_support = set()
    for row in icl_rows:
        if float(row["coverage_200_300"]) >= icl_coverage_floor:
            key = (row["id"], row["filter"])
            if key in icl_support:
                raise ValueError(f"duplicate ICL pair support row {key}")
            icl_support.add(key)
    satellite_support = set()
    for row in sat_rows:
        if (row["radial_bin_kpc"] == "50-75"
                and float(row["target_coverage"]) >= satellite_coverage_floor
                and float(row["control_coverage"]) >= satellite_coverage_floor):
            satellite_support.add((row["cluster_id"], row["filter"]))
    selected = sorted(icl_support & satellite_support,
                      key=lambda item: (int(item[0]), item[1]))
    if len({item[0] for item in selected}) != len(selected):
        raise ValueError("a cluster is represented by multiple filters in the joined sample")
    return selected


def summarize(records: list[dict]) -> dict:
    by_filter = {}
    for filt in sorted({r["filter"] for r in records}):
        subset = [r for r in records if r["filter"] == filt]
        metrics = {}
        for key in ("icl_response_to_diffuse_unit_injection",
                    "icl_plane_corrected_response_to_diffuse_unit_injection",
                    "icl_response_to_all_satellite_envelopes",
                    "icl_plane_corrected_response_to_all_satellite_envelopes",
                    "satellite_response_to_smooth_diffuse",
                    "satellite_response_to_satellite_unit_injection"):
            vals = [r["transfer_response"][key] for r in subset]
            metrics[key] = {
                "median": statistics.median(vals),
                "q16": float(np.quantile(vals, 0.16)),
                "q84": float(np.quantile(vals, 0.84)),
            }
        by_filter[filt] = {"n_clusters": len(subset), "responses": metrics}
    return by_filter


def run_batch(data: Path, output_dir: Path,
              icl_masks_dir: Path = DEFAULT_ICL_MASKS,
              satellite_prep_dir: Path = DEFAULT_SAT_PREP,
              icl_pairs_path: Path = DEFAULT_ICL_PAIRS,
              satellite_profiles_path: Path = DEFAULT_SATELLITE_TABLE,
              diffuse_amplitude: float = 1.0,
              satellite_amplitude: float = 1.0,
              diffuse_scale_kpc: float = 100.0,
              satellite_scale_kpc: float = 15.0,
              rerun_detection_masks: bool = False,
              unresolved_bcg_policy: str = "error",
              psf_bank: PositionPSFBank | None = None) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    selected = select_jointly_supported_pairs(icl_pairs_path, satellite_profiles_path)
    n_joint_supported = len(selected)
    if psf_bank is not None:
        selected = [pair for pair in selected if pair in psf_bank.supported_frames]
        if len(selected) != len(psf_bank.supported_frames):
            missing = sorted(psf_bank.supported_frames - set(selected))
            raise ValueError(f"approved PSF frame roles lack joint annulus support: {missing[:10]}")
    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.time()
    records, errors = [], []
    for index, (ident, filt) in enumerate(selected, 1):
        path = output_dir / f"cluster_{ident}.json"
        print(f"[{index}/{len(selected)}] joint injection cluster {ident} ({filt})", flush=True)
        try:
            result = run_frame(
                data, ident, path, icl_masks_dir, satellite_prep_dir, icl_pairs_path,
                diffuse_amplitude, satellite_amplitude, diffuse_scale_kpc,
                satellite_scale_kpc, rerun_detection_masks, unresolved_bcg_policy,
                psf_bank)
            if result["filter"] != filt:
                raise ValueError(f"image filter mismatch for cluster {ident}")
            records.append({
                "cluster_id": ident,
                "filter": filt,
                "path": str(path),
                "sha256": sha256(path),
                "support": result["support"],
                "icl_coverage_200_300": result["measurements_by_scenario"]["baseline"]["icl_coverage_200_300"],
                "transfer_response": result["transfer_response"],
            })
        except Exception as exc:
            errors.append({"cluster_id": ident, "filter": filt,
                           "error_type": type(exc).__name__, "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
    summary = {
        "schema": ("hsc-joint-icl-satellite-psf-redetection-batch-v1"
                   if psf_bank is not None else
                   "hsc-joint-icl-satellite-redetection-batch-v1"
                   if rerun_detection_masks else
                   "hsc-joint-icl-satellite-fixed-mask-batch-v1"),
        "status": "complete" if not errors else "partial",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "selection": {
            "icl_bin_kpc": [200.0, 300.0],
            "icl_common_coverage_floor": 0.40,
            "satellite_bin_kpc": [50.0, 75.0],
            "satellite_target_and_control_coverage_floor": 0.60,
            "n_selected_cluster_filter_pairs": len(selected),
            "n_joint_supported_before_psf_intersection": n_joint_supported,
            "selection_limited_to_approved_psf_positions": psf_bank is not None,
            "selected_pairs": [{"cluster_id": ident, "filter": filt} for ident, filt in selected],
        },
        "injection": {
            "diffuse_amplitude_at_100_njy_arcsec2": diffuse_amplitude,
            "diffuse_scale_kpc": diffuse_scale_kpc,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": satellite_amplitude,
            "satellite_scale_kpc": satellite_scale_kpc,
            "masks_and_segmentation_fixed": not rerun_detection_masks,
            "detection_rerun": rerun_detection_masks,
            "unresolved_bcg_policy": unresolved_bcg_policy if rerun_detection_masks else None,
            "poisson_noise_added": False,
            "psf_convolution_applied": psf_bank is not None,
            "psf_model": ("position-dependent HSC-I coadd model PSF core"
                          if psf_bank is not None else None),
            "original_fits_modified": False,
        },
        "n_completed": len(records),
        "n_errors": len(errors),
        "by_filter": summarize(records),
        "records": records,
        "errors": errors,
        "inputs": {
            str(icl_pairs_path): sha256(icl_pairs_path),
            str(satellite_profiles_path): sha256(satellite_profiles_path),
        },
        "code_sha256": {
            str(Path(__file__).resolve()): SCRIPT_SHA256_AT_IMPORT,
            str(SINGLE_FRAME_SCRIPT.resolve()): SINGLE_FRAME_SHA256_AT_IMPORT,
            str(POSITION_PSF_SCRIPT.resolve()): POSITION_PSF_SHA256_AT_IMPORT,
        },
        "elapsed_seconds": time.time() - started,
        "limitations": [
            ("This is a support-selected redetection/mask-response sensitivity, not a science measurement."
             if rerun_detection_masks else "This is a support-selected fixed-mask measurement-transfer sensitivity, not a science measurement."),
            (("SEP source detection and K2.5 masks are regenerated per injection, but injected-source photon noise and HSC upstream sky processing are not modeled."
              if rerun_detection_masks else "It does not rerun source detection, grow masks, or model upstream sky subtraction.")
             if psf_bank is not None else
             ("SEP source detection and K2.5 masks are regenerated per injection, but injected-source photon noise, PSF convolution and HSC upstream sky processing are not modeled."
              if rerun_detection_masks else "It does not rerun source detection, grow masks, convolve the injected shapes with PSFs, or model upstream sky subtraction.")),
            *( ["Only compact PSF cores are used; extended wings and source Poisson noise remain unmodeled."]
               if psf_bank is not None else [] ),
            "A missing or ambiguous central detection never establishes BCG identity; mask_all masks every detection conservatively.",
            "The artificial components and diagonal-variance estimator do not establish population completeness or flux covariance.",
            "The paired annular sample is selected by existing exploratory coverage cuts and is not an independent validation sample.",
        ],
    }
    if psf_bank is not None:
        summary["inputs"].update({
            str(psf_bank.plan_path): sha256(psf_bank.plan_path),
            str(psf_bank.mapping_path): sha256(psf_bank.mapping_path),
            str(psf_bank.audit_path): sha256(psf_bank.audit_path),
        })
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--icl-masks", type=Path, default=DEFAULT_ICL_MASKS)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_SAT_PREP)
    parser.add_argument("--icl-pairs", type=Path, default=DEFAULT_ICL_PAIRS)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_SATELLITE_TABLE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--diffuse-amplitude", type=float, default=1.0)
    parser.add_argument("--satellite-amplitude", type=float, default=1.0)
    parser.add_argument("--diffuse-scale-kpc", type=float, default=100.0)
    parser.add_argument("--satellite-scale-kpc", type=float, default=15.0)
    parser.add_argument("--rerun-detection-masks", action="store_true",
                        help="rerun SEP detection and regenerate source masks for each injection scenario")
    parser.add_argument("--unresolved-bcg-policy", choices=("error", "mask_all"), default="error",
                        help="on ambiguous/missing central detection, fail or mask every detected source")
    parser.add_argument("--psf-plan-dir", type=Path,
                        help="apply audited position-dependent HSC model PSF cores from this approved plan")
    parser.add_argument("--psf-download-dir", type=Path,
                        help="directory containing the audited PSF batch")
    parser.add_argument("--psf-test-center", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits"))
    args = parser.parse_args()
    if bool(args.psf_plan_dir) != bool(args.psf_download_dir):
        parser.error("--psf-plan-dir and --psf-download-dir must be supplied together")
    psf_bank = (PositionPSFBank(args.psf_plan_dir.resolve(), args.psf_download_dir.resolve(),
                                args.psf_test_center.resolve())
                if args.psf_plan_dir else None)
    result = run_batch(args.data.resolve(), args.output_dir.resolve(),
                       args.icl_masks.resolve(), args.satellite_prep.resolve(),
                       args.icl_pairs.resolve(), args.satellite_profiles.resolve(),
                       args.diffuse_amplitude, args.satellite_amplitude,
                       args.diffuse_scale_kpc, args.satellite_scale_kpc,
                       args.rerun_detection_masks, args.unresolved_bcg_policy, psf_bank)
    print(json.dumps({k: result[k] for k in
                      ("status", "n_completed", "n_errors", "by_filter", "elapsed_seconds")},
                     indent=2))


if __name__ == "__main__":
    main()
