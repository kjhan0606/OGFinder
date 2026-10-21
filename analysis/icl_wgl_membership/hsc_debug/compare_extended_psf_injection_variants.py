#!/usr/bin/env python3
"""Pair symmetric/asymmetric extended-PSF injection outputs on identical frames."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from validate_extended_psf_injection_batch import METRICS, audit_batch


INJECTION_KEYS = (
    "diffuse_amplitude_at_100_njy_arcsec2",
    "satellite_amplitude_at_50_per_galaxy_njy_arcsec2",
    "diffuse_scale_kpc",
    "satellite_scale_kpc",
    "n_injected_satellite_centers",
    "unresolved_bcg_policy",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize_paired_differences(values) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("paired differences must be a non-empty finite vector")
    return {
        "n_frames": int(len(values)),
        "median_asymmetric_minus_symmetric": float(np.median(values)),
        "q16": float(np.quantile(values, 0.16)),
        "q84": float(np.quantile(values, 0.84)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "max_abs": float(np.max(np.abs(values))),
        "q16_q84_contains_zero": bool(np.quantile(values, 0.16) <= 0 <= np.quantile(values, 0.84)),
    }


def compare_batches(symmetric_dir: Path, asymmetric_dir: Path, *, expected_frames: int,
                    filter_name: str = "i2") -> dict:
    if filter_name not in {"i", "i2"}:
        raise ValueError("filter_name must be i or i2")
    symmetric_dir, asymmetric_dir = Path(symmetric_dir).resolve(), Path(asymmetric_dir).resolve()
    audit_batch(symmetric_dir, expected_filter=filter_name, expected_frames=expected_frames)
    audit_batch(asymmetric_dir, expected_filter=filter_name, expected_frames=expected_frames)
    symmetric_summary = json.loads((symmetric_dir / "summary.json").read_text(encoding="utf-8"))
    asymmetric_summary = json.loads((asymmetric_dir / "summary.json").read_text(encoding="utf-8"))
    if symmetric_summary["reference"]["name"] != "HSC_PSF_i.fits":
        raise ValueError("symmetric batch does not use HSC_PSF_i.fits")
    if asymmetric_summary["reference"]["name"] != "HSC_PSF_i_asymmetric.fits":
        raise ValueError("asymmetric batch does not use HSC_PSF_i_asymmetric.fits")
    if symmetric_summary["reference"]["sha256"] == asymmetric_summary["reference"]["sha256"]:
        raise ValueError("symmetric and asymmetric reference hashes unexpectedly match")

    def records_by_cluster(summary):
        rows = {}
        for item in summary["records"]:
            path = Path(item["path"]).resolve()
            record = json.loads(path.read_text(encoding="utf-8"))
            ident = str(record["cluster_id"])
            if ident in rows:
                raise ValueError(f"duplicate cluster {ident}")
            rows[ident] = record
        return rows

    sym = records_by_cluster(symmetric_summary)
    asym = records_by_cluster(asymmetric_summary)
    if set(sym) != set(asym) or len(sym) != expected_frames:
        raise ValueError("symmetric and asymmetric batches do not share exact frame support")
    for ident in sym:
        if sym[ident].get("support") != asym[ident].get("support"):
            raise ValueError(f"injection support differs for cluster {ident}")
        for key in INJECTION_KEYS:
            if sym[ident]["injection"].get(key) != asym[ident]["injection"].get(key):
                raise ValueError(f"injection parameter {key} differs for cluster {ident}")

    differences = {}
    for metric in METRICS:
        delta = [float(asym[ident]["transfer_response"][metric])
                 - float(sym[ident]["transfer_response"][metric])
                 for ident in sorted(sym, key=int)]
        differences[metric] = summarize_paired_differences(delta)
    return {
        "schema": "astrafek-extended-psf-injection-variant-comparison-v1",
        "status": "complete_paired_psf_model_sensitivity",
        "frame_filter": filter_name,
        "n_matched_frames": len(sym),
        "matched_cluster_ids": sorted(sym, key=int),
        "symmetric_summary_sha256": sha256(symmetric_dir / "summary.json"),
        "asymmetric_summary_sha256": sha256(asymmetric_dir / "summary.json"),
        "symmetric_reference_sha256": symmetric_summary["reference"]["sha256"],
        "asymmetric_reference_sha256": asymmetric_summary["reference"]["sha256"],
        "paired_asymmetric_minus_symmetric": differences,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "limitations": [
            "The PSF variants share the public star sample; this is model-form sensitivity, not independent replication.",
            ("The HSC-I model is only a proxy sensitivity for i2 and is not a filter-specific calibration."
             if filter_name == "i2" else
             "The HSC-I model is not position-, epoch-, or seeing-matched to individual i coadds."),
            "Quantiles describe selected-frame variation, not confidence intervals or a population estimate.",
            "No observed ICL correction or science measurement is produced.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symmetric-dir", type=Path, required=True)
    parser.add_argument("--asymmetric-dir", type=Path, required=True)
    parser.add_argument("--filter", choices=("i", "i2"), default="i2")
    parser.add_argument("--expected-frames", type=int, default=44)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing comparison: {output}")
    result = compare_batches(args.symmetric_dir, args.asymmetric_dir,
                             expected_frames=args.expected_frames,
                             filter_name=args.filter)
    result["comparison_script_sha256"] = sha256(Path(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "matched_frames": result["n_matched_frames"],
                      "science_measurement": False,
                      "output_sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
