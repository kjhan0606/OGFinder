#!/usr/bin/env python3
"""Summarize paired diffuse-profile transfer across audited PSF injection runs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
AUDITOR = Path(__file__).with_name("audit_joint_icl_satellite_psf_redetection.py")
REFERENCE_SCALES = {100.0, 500.0}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_summary(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    if not len(array) or not np.isfinite(array).all():
        raise ValueError("summary requires finite, non-empty values")
    return {"n": int(len(array)), "median": float(np.median(array)),
            "q16": float(np.quantile(array, 0.16)),
            "q84": float(np.quantile(array, 0.84))}


def load_run(scale: float, directory: Path) -> dict:
    summary_path = directory / "summary.json"
    audit_path = directory / "psf_redetection_audit.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete" or summary.get("science_measurement") is not False:
        raise ValueError(f"run {scale:g} kpc is incomplete or marked as a science measurement")
    injection = summary.get("injection", {})
    if (not injection.get("psf_convolution_applied")
            or not injection.get("detection_rerun")
            or float(injection.get("diffuse_scale_kpc", -1)) != scale
            or float(injection.get("diffuse_amplitude_at_100_njy_arcsec2", -1)) != 1.0
            or float(injection.get("satellite_scale_kpc", -1)) != 15.0):
        raise ValueError(f"run {scale:g} kpc does not match the bracketing injection design")
    if (audit.get("status") != "verified" or audit.get("science_measurement") is not False
            or audit.get("n_frame_results_checked") != 57
            or audit.get("batch_summary_sha256") != sha256(summary_path)
            or audit.get("audit_script_sha256") != sha256(AUDITOR)):
        raise ValueError(f"run {scale:g} kpc failed its independent PSF/redetection audit")
    by_pair = {}
    for row in summary.get("records", []):
        key = (str(row["cluster_id"]), str(row["filter"]).lower())
        if key in by_pair:
            raise ValueError(f"duplicate frame in {scale:g} kpc run: {key}")
        by_pair[key] = row
    if len(by_pair) != 57 or summary.get("n_errors") != 0:
        raise ValueError(f"unexpected frame support in {scale:g} kpc run")
    code_hashes = summary.get("code_sha256")
    inputs = summary.get("inputs")
    if not code_hashes or not inputs:
        raise ValueError("run lacks code or input provenance")
    return {"scale_kpc": scale, "directory": str(directory.resolve()),
            "summary_path": str(summary_path.resolve()), "summary_sha256": sha256(summary_path),
            "audit_path": str(audit_path.resolve()), "audit_sha256": sha256(audit_path),
            "summary": summary, "records": by_pair,
            "code_sha256": code_hashes, "inputs": inputs}


def summarize_grid(runs: dict[float, dict]) -> dict:
    scales = sorted(runs)
    if not REFERENCE_SCALES.issubset(runs):
        raise ValueError("the scale grid must include the existing 100 and 500 kpc anchors")
    reference = runs[scales[0]]
    expected_pairs = set(reference["records"])
    for scale in scales[1:]:
        run = runs[scale]
        if set(run["records"]) != expected_pairs:
            raise ValueError("scale runs do not have identical cluster/filter support")
        if run["code_sha256"] != reference["code_sha256"] or run["inputs"] != reference["inputs"]:
            raise ValueError("scale runs differ in code or scientific inputs")

    per_filter = {}
    for filt in ("i", "i2"):
        keys = sorted((key for key in expected_pairs if key[1] == filt),
                      key=lambda key: int(key[0]))
        ratios: dict[float, dict[tuple[str, str], float]] = {scale: {} for scale in scales}
        metrics = {str(int(scale)): {} for scale in scales}
        nonpositive = {str(int(scale)): 0 for scale in scales}
        for scale in scales:
            raw_values, corrected_values = [], []
            for key in keys:
                response = runs[scale]["records"][key]["transfer_response"]
                raw = response.get("icl_response_to_diffuse_unit_injection")
                corrected = response.get("icl_plane_corrected_response_to_diffuse_unit_injection")
                if raw is None or corrected is None:
                    continue
                raw, corrected = float(raw), float(corrected)
                if not np.isfinite(raw) or not np.isfinite(corrected):
                    raise ValueError(f"non-finite transfer response at {scale:g} kpc/{key}")
                raw_values.append(raw)
                corrected_values.append(corrected)
                if abs(raw) > 1e-10:
                    ratios[scale][key] = corrected / raw
                    if raw <= 0:
                        nonpositive[str(int(scale))] += 1
            if not raw_values:
                raise ValueError(f"no diffuse ICL response values at {scale:g} kpc/{filt}")
            metrics[str(int(scale))] = {
                "raw_response": finite_summary(raw_values),
                "plane_corrected_response": finite_summary(corrected_values),
                "plane_retention_ratio": finite_summary(list(ratios[scale].values())),
                "nonpositive_raw_response_frames": nonpositive[str(int(scale))],
            }

        jointly_valid = set.intersection(*(set(ratios[scale]) for scale in scales))
        if not jointly_valid:
            raise ValueError(f"no common valid retention ratios across scales for {filt}")
        paired_changes = {}
        anchor = 100.0
        for scale in scales:
            changes = [ratios[scale][key] - ratios[anchor][key] for key in jointly_valid]
            paired_changes[str(int(scale))] = finite_summary(changes)
        per_filter[filt] = {
            "n_frames": len(keys), "scales": metrics,
            "paired_retention_change_from_100kpc": paired_changes,
            "n_common_frames_for_paired_retention": len(jointly_valid),
        }

    return {
        "schema": "astrafek-icl-psf-diffuse-scale-grid-v1",
        "status": "verified_exploratory_scale_grid",
        "science_measurement": False,
        "method": "same-frame paired exponential diffuse injections; local HSC-I model-PSF cores; source detection and K2.5 masks regenerated per injection; plane estimator applied to saved measurements",
        "uncertainty_note": "reported q16/q84 values describe frame-to-frame spread, not confidence intervals; full spatial covariance and upstream HSC sky processing are absent",
        "scale_lengths_kpc": scales,
        "filters": per_filter,
        "runs": {str(int(scale)): {key: runs[scale][key] for key in
                                    ("directory", "summary_sha256", "audit_path", "audit_sha256")}
                 for scale in scales},
        "common_code_sha256": reference["code_sha256"],
        "common_inputs": reference["inputs"],
        "limitations": [
            "This is a support-selected injection transfer test, not an observed ICL detection or multiplicative correction.",
            "The diffuse exponential scales bracket the 200-300 kpc measurement annulus but are not fitted ICL models.",
            "Only compact position-dependent PSF cores are used; extended wings, injected photon noise, correlated coadd noise, extinction and passband uncertainty are not modeled.",
            "The plane estimator is downstream of HSC processing; this does not recover unknown upstream sky subtraction.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, metavar="SCALE_KPC=DIR")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    runs = {}
    for item in args.run:
        try:
            scale_text, directory_text = item.split("=", 1)
            scale = float(scale_text)
        except ValueError as exc:
            parser.error(f"invalid --run {item!r}; expected SCALE_KPC=DIR")
        if scale in runs or scale <= 0:
            parser.error(f"duplicate or invalid scale: {scale}")
        runs[scale] = load_run(scale, Path(directory_text))
    result = summarize_grid(runs)
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    result["output_path"] = str(args.output.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "scale_lengths_kpc": result["scale_lengths_kpc"],
                      "filters": result["filters"], "output": str(args.output)},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
