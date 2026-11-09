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
    RUNTIME_CODE_SHA256_AT_IMPORT,
    run_frame,
    sha256,
)
from position_psf import PositionPSFBank, SupportedPairsPositionPSFBank


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


def select_icl_supported_pairs(icl_pairs_path: Path,
                                icl_coverage_floor: float = 0.40) -> list[tuple[str, str]]:
    """Select the parent ICL-support sample, retaining satellite support failures."""
    rows = read_tsv(icl_pairs_path)
    if not rows or not {"id", "filter", "coverage_200_300"}.issubset(rows[0]):
        raise ValueError("ICL pair table is missing the expected common-coverage fields")
    selected = sorted({
        (str(row["id"]), str(row["filter"]))
        for row in rows if float(row["coverage_200_300"]) >= icl_coverage_floor
    }, key=lambda item: (int(item[0]), item[1]))
    if len(selected) != sum(float(row["coverage_200_300"]) >= icl_coverage_floor for row in rows):
        raise ValueError("duplicate ICL-supported cluster/filter pair rows")
    if len({item[0] for item in selected}) != len(selected):
        raise ValueError("an ICL-supported cluster is represented by multiple filters")
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
            vals = [float(r["transfer_response"][key])
                    for r in subset
                    if r["transfer_response"].get(key) is not None
                    and np.isfinite(float(r["transfer_response"][key]))]
            metrics[key] = {
                "n_valid": len(vals),
                "median": statistics.median(vals) if vals else None,
                "q16": float(np.quantile(vals, 0.16)) if vals else None,
                "q84": float(np.quantile(vals, 0.84)) if vals else None,
            }
        shell = [r["transfer_response"].get("satellite_shell_excluded_icl") for r in subset
                 if r["transfer_response"].get("satellite_shell_excluded_icl") is not None]
        if shell:
            metrics["satellite_shell_excluded_icl"] = {}
            for key in ("icl_response_to_diffuse_unit_injection",
                        "icl_response_to_all_satellite_envelopes",
                        "both_minus_sum_residual_icl_njy_arcsec2"):
                vals = np.asarray([item[key] for item in shell if item.get(key) is not None], dtype=float)
                metrics["satellite_shell_excluded_icl"][key] = {
                    "n_valid": int(len(vals)),
                    "median": float(np.median(vals)) if len(vals) else None,
                    "q16": float(np.quantile(vals, 0.16)) if len(vals) else None,
                    "q84": float(np.quantile(vals, 0.84)) if len(vals) else None,
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
              psf_bank: PositionPSFBank | None = None,
              exclude_supported_satellite_shells: bool = False,
              association_pairs_path: Path | None = None,
              satellite_coverage_floor: float = 0.60,
              minimum_controls: int = 3,
              icl_supported_only: bool = False) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    if not 0. < satellite_coverage_floor <= 1. or minimum_controls < 1:
        raise ValueError("invalid satellite coverage floor or minimum control count")
    if exclude_supported_satellite_shells and rerun_detection_masks:
        raise ValueError("shell-excluded ICL injection currently requires fixed observed segmentation")
    if icl_supported_only and exclude_supported_satellite_shells:
        raise ValueError("shell exclusion requires satellite-supported rows and cannot use the ICL-only parent selection")
    if icl_supported_only and association_pairs_path is not None:
        raise ValueError("--icl-supported-only cannot be combined with the fixed association-pair subset")
    selected = (select_icl_supported_pairs(icl_pairs_path, icl_coverage_floor=0.40)
                if icl_supported_only else
                select_jointly_supported_pairs(
                    icl_pairs_path, satellite_profiles_path,
                    icl_coverage_floor=0.40,
                    satellite_coverage_floor=satellite_coverage_floor))
    joint_supported_pairs = list(selected)
    n_joint_supported = None if icl_supported_only else len(selected)
    association_scope = None
    if association_pairs_path is not None:
        association_pairs_path = association_pairs_path.resolve()
        association_summary_path = association_pairs_path.parent / "summary.json"
        association_audit_path = association_pairs_path.parent / "independent_audit.json"
        if not association_summary_path.is_file() or not association_audit_path.is_file():
            raise FileNotFoundError("association table requires its summary and independent audit")
        association_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
        association_audit = json.loads(association_audit_path.read_text(encoding="utf-8"))
        if (association_audit.get("status") != "verified"
                or association_audit.get("science_measurement") is not False
                or association_audit.get("association_rows") != 74
                or association_audit.get("association_table_sha256") != sha256(association_pairs_path)):
            raise ValueError("association table lacks the verified fixed 74-pair sample provenance")
        association_rows = read_tsv(association_pairs_path)
        association_keys = [(row["target_id"], row["filter"]) for row in association_rows]
        if len(association_keys) != 74 or len(set(association_keys)) != 74:
            raise ValueError("expected exactly 74 unique target/filter association pairs")
        before_association_filter = set(selected)
        selected = sorted(before_association_filter & set(association_keys),
                          key=lambda item: (int(item[0]), item[1]))
        association_scope = {
            "path": str(association_pairs_path),
            "sha256": sha256(association_pairs_path),
            "summary_sha256": sha256(association_summary_path),
            "independent_audit_sha256": sha256(association_audit_path),
            "n_fixed_association_pairs": 74,
            "n_association_and_joint_support": len(selected),
        }
    shell_rows_by_pair = {}
    if exclude_supported_satellite_shells:
        for row in read_tsv(satellite_profiles_path):
            if (row.get("cohort") == "full_sample" and row.get("radial_bin_kpc") == "50-75"
                    and float(row["target_coverage"]) >= satellite_coverage_floor
                    and float(row["control_coverage"]) >= satellite_coverage_floor
                    and int(row["n_controls"]) >= minimum_controls):
                shell_rows_by_pair.setdefault((row["cluster_id"], row["filter"]), []).append(row)
        missing_shells = [pair for pair in selected if pair not in shell_rows_by_pair]
        if missing_shells:
            raise ValueError(f"jointly supported pairs lack satellite-shell pixel rows: {missing_shells[:10]}")
    if psf_bank is not None:
        missing_psf_support = sorted(set(selected) - psf_bank.supported_frames)
        if missing_psf_support:
            raise ValueError(f"selected frames lack complete approved PSF role support: {missing_psf_support[:10]}")
    if association_scope is not None:
        association_scope["n_association_joint_psf_intersection"] = len(selected)
        association_scope["n_association_pairs_without_joint_or_psf_support"] = 74-len(selected)
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
                psf_bank, shell_rows_by_pair.get((ident, filt))
                if exclude_supported_satellite_shells else None,
                satellite_coverage_floor, minimum_controls)
            if result["filter"] != filt:
                raise ValueError(f"image filter mismatch for cluster {ident}")
            records.append({
                "cluster_id": ident,
                "filter": filt,
                "path": str(path),
                "sha256": sha256(path),
                "support": result["support"],
                "support_status": result["support_status"],
                "satellite_support_gaps": result["satellite_support_gaps"],
                "icl_coverage_200_300": result["measurements_by_scenario"]["baseline"]["icl_coverage_200_300"],
                "transfer_response": result["transfer_response"],
            })
        except Exception as exc:
            errors.append({"cluster_id": ident, "filter": filt,
                           "error_type": type(exc).__name__, "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
    summary = {
        "schema": ("hsc-joint-icl-satellite-shell-exclusion-batch-v1"
                   if exclude_supported_satellite_shells else
                   "hsc-joint-icl-satellite-psf-redetection-batch-v1"
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
            "satellite_target_and_control_coverage_floor": satellite_coverage_floor,
            "minimum_azimuthal_controls": minimum_controls,
            "satellite_profile_support_required_for_frame_selection": not icl_supported_only,
            "frame_selection_mode": "icl_supported_parent" if icl_supported_only else "joint_icl_satellite_supported",
            "n_selected_cluster_filter_pairs": len(selected),
            "n_frames_selected_before_psf_intersection": len(joint_supported_pairs),
            "n_joint_supported_before_psf_intersection": n_joint_supported,
            "selection_limited_to_approved_psf_positions": psf_bank is not None,
            "fixed_association_scope": association_scope,
            "selected_pairs": [{"cluster_id": ident, "filter": filt} for ident, filt in selected],
        },
        "injection": {
            "diffuse_amplitude_at_100_njy_arcsec2": diffuse_amplitude,
            "diffuse_scale_kpc": diffuse_scale_kpc,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": satellite_amplitude,
            "satellite_scale_kpc": satellite_scale_kpc,
            "satellite_profile_coverage_floor": satellite_coverage_floor,
            "minimum_azimuthal_controls": minimum_controls,
            "masks_and_segmentation_fixed": not rerun_detection_masks,
            "detection_rerun": rerun_detection_masks,
            "unresolved_bcg_policy": unresolved_bcg_policy if rerun_detection_masks else None,
            "poisson_noise_added": False,
            "psf_convolution_applied": psf_bank is not None,
            "psf_model": ("position-dependent HSC-I coadd model PSF core"
                          if psf_bank is not None else None),
            "original_fits_modified": False,
            "exclude_supported_satellite_shells_from_icl": exclude_supported_satellite_shells,
            "shell_exclusion_uses_observed_fixed_segmentation": exclude_supported_satellite_shells,
        },
        "n_completed": len(records),
        "n_errors": len(errors),
        "n_frames_with_satellite_support_gaps": sum(
            bool(row["satellite_support_gaps"]) for row in records),
        "satellite_support_gaps_by_scenario": {
            scenario: sum(scenario in row["satellite_support_gaps"] for row in records)
            for scenario in ("baseline", "icl_only", "satellite_only", "both")
        },
        "by_filter": summarize(records),
        "records": records,
        "errors": errors,
        "inputs": {
            str(icl_pairs_path): sha256(icl_pairs_path),
            str(satellite_profiles_path): sha256(satellite_profiles_path),
        },
        "code_sha256": {
            **RUNTIME_CODE_SHA256_AT_IMPORT,
            str(Path(__file__).resolve()): SCRIPT_SHA256_AT_IMPORT,
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
            str(path): sha256(path)
            for path in getattr(psf_bank, "provenance_paths", (
                psf_bank.plan_path, psf_bank.mapping_path, psf_bank.audit_path))
        })
        if getattr(psf_bank, "request_approval_verified", False):
            summary["supported_satellite_psf_request_approval"] = {
                "verified": True,
                "approved_plan_sha256": psf_bank.approved_satellite_plan_sha256,
            }
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
    parser.add_argument("--satellite-coverage-floor", type=float, default=0.60,
                        help="target and control floor used to select and measure satellite profiles")
    parser.add_argument("--minimum-controls", type=int, default=3,
                        help="minimum valid azimuthal controls per satellite profile")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--diffuse-amplitude", type=float, default=1.0)
    parser.add_argument("--satellite-amplitude", type=float, default=1.0)
    parser.add_argument("--diffuse-scale-kpc", type=float, default=100.0)
    parser.add_argument("--satellite-scale-kpc", type=float, default=15.0)
    parser.add_argument("--rerun-detection-masks", action="store_true",
                        help="rerun SEP detection and regenerate source masks for each injection scenario")
    parser.add_argument("--exclude-supported-satellite-shells-from-icl", action="store_true",
                        help="measure ICL transfer after excluding the observed supported 50--75 kpc satellite pixel union")
    parser.add_argument("--association-pairs", type=Path,
                        help="restrict processing to a verified fixed ICL--satellite association pair table")
    parser.add_argument("--icl-supported-only", action="store_true",
                        help="inject all ICL-supported frames and retain satellite-profile support failures in the output")
    parser.add_argument("--unresolved-bcg-policy", choices=("error", "mask_all"), default="error",
                        help="on ambiguous/missing central detection, fail or mask every detected source")
    parser.add_argument("--psf-plan-dir", type=Path,
                        help="apply audited position-dependent HSC model PSF cores from this approved plan")
    parser.add_argument("--psf-download-dir", type=Path,
                        help="directory containing the audited PSF batch")
    parser.add_argument("--psf-test-center", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits"))
    parser.add_argument("--target-center-psf-plan-dir", type=Path,
                        help="audited 1,542-position target/matched-control center PSF plan")
    parser.add_argument("--target-center-psf-download-dir", type=Path,
                        help="audited response for --target-center-psf-plan-dir")
    parser.add_argument("--supported-satellite-psf-plan-dir", type=Path,
                        help="approved supported-pairs plan for additional satellite PSF cores")
    parser.add_argument("--supported-satellite-psf-download-dir", type=Path,
                        help="audited response for --supported-satellite-psf-plan-dir")
    parser.add_argument("--approved-supported-satellite-psf-plan-sha256",
                        help="exact 64-character plan SHA-256 explicitly approved for the satellite PSF request")
    args = parser.parse_args()
    base_pair = (args.psf_plan_dir, args.psf_download_dir)
    supplemental = (args.target_center_psf_plan_dir, args.target_center_psf_download_dir,
                    args.supported_satellite_psf_plan_dir, args.supported_satellite_psf_download_dir)
    if bool(base_pair[0]) != bool(base_pair[1]):
        parser.error("--psf-plan-dir and --psf-download-dir must be supplied together")
    if any(value is not None for value in supplemental) and not all(supplemental):
        parser.error("the four supplemental PSF plan/download directories must be supplied together")
    if any(supplemental) and not all(base_pair):
        parser.error("supplemental PSF batches require the original --psf-plan-dir and --psf-download-dir")
    if all(supplemental) and not args.approved_supported_satellite_psf_plan_sha256:
        parser.error("supplemental satellite PSFs require --approved-supported-satellite-psf-plan-sha256")
    if args.approved_supported_satellite_psf_plan_sha256 and not all(supplemental):
        parser.error("approved satellite PSF plan hash requires all supplemental PSF directories")
    if all(supplemental):
        psf_bank = SupportedPairsPositionPSFBank(
            args.psf_plan_dir.resolve(), args.psf_download_dir.resolve(),
            args.psf_test_center.resolve(), args.target_center_psf_plan_dir.resolve(),
            args.target_center_psf_download_dir.resolve(),
            args.supported_satellite_psf_plan_dir.resolve(),
            args.supported_satellite_psf_download_dir.resolve(),
            args.approved_supported_satellite_psf_plan_sha256)
    elif args.psf_plan_dir:
        psf_bank = PositionPSFBank(args.psf_plan_dir.resolve(), args.psf_download_dir.resolve(),
                                   args.psf_test_center.resolve())
    else:
        psf_bank = None
    result = run_batch(args.data.resolve(), args.output_dir.resolve(),
                       args.icl_masks.resolve(), args.satellite_prep.resolve(),
                       args.icl_pairs.resolve(), args.satellite_profiles.resolve(),
                       args.diffuse_amplitude, args.satellite_amplitude,
                       args.diffuse_scale_kpc, args.satellite_scale_kpc,
                       args.rerun_detection_masks, args.unresolved_bcg_policy, psf_bank,
                       args.exclude_supported_satellite_shells_from_icl,
                       args.association_pairs, args.satellite_coverage_floor,
                       args.minimum_controls, args.icl_supported_only)
    print(json.dumps({k: result[k] for k in
                      ("status", "n_completed", "n_errors", "by_filter", "elapsed_seconds")},
                     indent=2))


if __name__ == "__main__":
    main()
