#!/usr/bin/env python3
"""Compare matched 100/500 kpc PSF-convolved ICL injection responses.

This is an offline audit of two existing Astrafek injection runs. It makes no
HSC requests and does not estimate an observed ICL--satellite relation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


DIFFUSE_RAW = "icl_response_to_diffuse_unit_injection"
DIFFUSE_PLANE = "icl_plane_corrected_response_to_diffuse_unit_injection"
OTHER_METRICS = (
    "icl_response_to_all_satellite_envelopes",
    "icl_plane_corrected_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_summary(values: list[float]) -> dict:
    array = np.asarray([float(value) for value in values if math.isfinite(float(value))])
    if not len(array):
        return {"n": 0, "median": None, "q16": None, "q84": None}
    return {
        "n": int(len(array)),
        "median": float(np.median(array)),
        "q16": float(np.quantile(array, 0.16)),
        "q84": float(np.quantile(array, 0.84)),
    }


def load_run(run_dir: Path, expected_scale: float) -> dict:
    run_dir = run_dir.resolve()
    summary_path = run_dir / "summary.json"
    audit_path = run_dir / "psf_redetection_audit.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete" or summary.get("n_errors") != 0
            or summary.get("n_completed") != 57 or summary.get("science_measurement") is not False):
        raise ValueError(f"run is not a complete exploratory 57-frame batch: {run_dir}")
    if summary.get("injection", {}).get("diffuse_scale_kpc") != expected_scale:
        raise ValueError(f"unexpected diffuse scale in {run_dir}")
    if (audit.get("status") != "verified" or audit.get("n_approved_frames") != 57
            or audit.get("n_frame_results_checked") != 57
            or audit.get("batch_summary_sha256") != sha256(summary_path)
            or audit.get("science_measurement") is not False):
        raise ValueError(f"PSF/frame audit does not verify this run: {run_dir}")

    code_hashes = summary.get("code_sha256", {})
    if len(code_hashes) < 8:
        raise ValueError(f"runtime provenance omits imported Astrafek modules: {run_dir}")
    for source, digest in code_hashes.items():
        if sha256(Path(source)) != digest:
            raise ValueError(f"run-time source changed since execution: {source}")

    records = {}
    for row in summary.get("records", []):
        pair = (str(row["cluster_id"]), str(row["filter"]).lower())
        path = Path(row["path"]).resolve()
        if path.parent != run_dir or sha256(path) != row.get("sha256"):
            raise ValueError(f"frame output path/hash mismatch: {pair}")
        result = json.loads(path.read_text(encoding="utf-8"))
        if (result.get("science_measurement") is not False
                or result.get("injection", {}).get("diffuse_scale_kpc") != expected_scale):
            raise ValueError(f"frame scale/science flag mismatch: {pair}")
        if (row.get("satellite_support_gaps") != result.get("satellite_support_gaps")
                or row.get("support_status") != result.get("support_status")):
            raise ValueError(f"frame support provenance mismatch: {pair}")
        records[pair] = result
    expected = {(str(row["cluster_id"]), str(row["filter"]).lower())
                for row in summary["selection"]["selected_pairs"]}
    if len(records) != 57 or set(records) != expected:
        raise ValueError(f"frame records differ from selected sample: {run_dir}")
    return {
        "run_dir": str(run_dir),
        "summary": summary,
        "summary_sha256": sha256(summary_path),
        "audit_sha256": sha256(audit_path),
        "code_sha256": code_hashes,
        "records": records,
    }


def compare_runs(reference: dict, extended: dict) -> dict:
    if reference["code_sha256"] != extended["code_sha256"]:
        raise ValueError("100 and 500 kpc runs were not produced by identical code")
    if reference["summary"].get("inputs") != extended["summary"].get("inputs"):
        raise ValueError("input, PSF plan, or position mapping differs between runs")
    if set(reference["records"]) != set(extended["records"]):
        raise ValueError("the 100 and 500 kpc frame samples do not match")

    output = {"filters": {}, "science_measurement": False}
    for filt in ("i", "i2"):
        pairs = sorted((pair for pair in reference["records"] if pair[1] == filt),
                       key=lambda pair: int(pair[0]))
        if not pairs:
            raise ValueError(f"no matched frames for filter {filt}")
        metrics = {}
        for key in (DIFFUSE_RAW, DIFFUSE_PLANE) + OTHER_METRICS:
            by_scale = {}
            for label, run in (("100kpc", reference), ("500kpc", extended)):
                values = [run["records"][pair]["transfer_response"].get(key)
                          for pair in pairs]
                by_scale[label] = finite_summary([value for value in values if value is not None])
            deltas = []
            for pair in pairs:
                low = reference["records"][pair]["transfer_response"].get(key)
                high = extended["records"][pair]["transfer_response"].get(key)
                if low is not None and high is not None:
                    deltas.append(float(high) - float(low))
            metrics[key] = {
                "100kpc": by_scale["100kpc"],
                "500kpc": by_scale["500kpc"],
                "paired_500_minus_100": finite_summary(deltas),
            }

        retention = {"100kpc": [], "500kpc": []}
        raw_nonpositive = {"100kpc": 0, "500kpc": 0}
        for pair in pairs:
            for label, run in (("100kpc", reference), ("500kpc", extended)):
                response = run["records"][pair]["transfer_response"]
                raw = response.get(DIFFUSE_RAW)
                corrected = response.get(DIFFUSE_PLANE)
                if raw is None or corrected is None or abs(float(raw)) < 1e-10:
                    continue
                if float(raw) <= 0:
                    raw_nonpositive[label] += 1
                retention[label].append(float(corrected) / float(raw))
        paired_retention = []
        for pair in pairs:
            ratios = []
            valid = True
            for run in (reference, extended):
                response = run["records"][pair]["transfer_response"]
                raw, corrected = response.get(DIFFUSE_RAW), response.get(DIFFUSE_PLANE)
                if raw is None or corrected is None or abs(float(raw)) < 1e-10:
                    valid = False
                    break
                ratios.append(float(corrected) / float(raw))
            if valid:
                paired_retention.append(ratios[1] - ratios[0])
        support_gaps = {
            label: sum(bool(run["records"][pair].get("satellite_support_gaps"))
                       for pair in pairs)
            for label, run in (("100kpc", reference), ("500kpc", extended))
        }
        output["filters"][filt] = {
            "n_matched_frames": len(pairs),
            "icl_plane_retention_ratio": {
                "definition": "plane-corrected diffuse ICL response divided by raw-annulus diffuse ICL response, evaluated frame by frame",
                "100kpc": finite_summary(retention["100kpc"]),
                "500kpc": finite_summary(retention["500kpc"]),
                "paired_500_minus_100": finite_summary(paired_retention),
                "raw_response_nonpositive_frames": raw_nonpositive,
            },
            "transfer_responses": metrics,
            "frames_with_satellite_support_gaps": support_gaps,
        }
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-100kpc", type=Path, required=True)
    parser.add_argument("--run-500kpc", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    reference = load_run(args.run_100kpc, 100.0)
    extended = load_run(args.run_500kpc, 500.0)
    comparison = compare_runs(reference, extended)
    result = {
        "schema": "astrafek-icl-psf-redetection-scale-comparison-v1",
        "status": "verified_exploratory",
        "science_measurement": False,
        "analysis_code_path": str(Path(__file__).resolve()),
        "analysis_code_sha256": sha256(Path(__file__).resolve()),
        "method": "paired, same-frame 100 versus 500 kpc diffuse exponential profile; position-dependent HSC-I model-PSF cores; SEP and K2.5 masks re-detected per scenario",
        "uncertainty_note": "q16/q84 values describe frame-to-frame spread, not confidence intervals; no full spatial covariance or upstream HSC sky transfer is included",
        "n_frames": 57,
        "n_frames_by_filter": {f: int(sum(pair[1] == f for pair in reference["records"]))
                                for f in ("i", "i2")},
        "runs": {
            "100kpc": {k: reference[k] for k in ("run_dir", "summary_sha256", "audit_sha256")},
            "500kpc": {k: extended[k] for k in ("run_dir", "summary_sha256", "audit_sha256")},
        },
        "runtime_code_sha256": reference["code_sha256"],
        "filters": comparison["filters"],
        "limitations": [
            "This compares injected profiles and does not measure the observed ICL--satellite relation.",
            "The 500 kpc profile is an extended stress case, not a fitted ICL model.",
            "Only compact position-dependent PSF cores are applied; extended PSF wings are not included.",
            "Injected-source Poisson noise, correlated coadd noise, HSC upstream sky subtraction, and extinction/passband uncertainty are not modeled.",
            "Satellite cross-responses are missing in a scenario when the regenerated 50--75 kpc estimator loses its minimum support; the ICL response remains separately audited.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_frames": result["n_frames"], "n_frames_by_filter": result["n_frames_by_filter"],
                      "icl_plane_retention": {f: result["filters"][f]["icl_plane_retention_ratio"]
                                              for f in result["filters"]}}, indent=2))
    print(f"comparison_sha256={sha256(args.output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
