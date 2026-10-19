#!/usr/bin/env python3
"""Compare public wide-PSF and position-dependent core injection responses.

This pairs existing, checksum-recorded outputs on the same i-band frames. It
does not rerun injections, download data, or promote the result to a science
measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_CORE = DEFAULT_DATA / "joint_icl_satellite_psf_redetection_2026-10-19_02"
DEFAULT_WIDE = DEFAULT_DATA / "icl_extended_psf_symmetric_2026-10-19_01"
METRICS = (
    "icl_response_to_diffuse_unit_injection",
    "icl_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)
INJECTION_KEYS = (
    "diffuse_amplitude_at_100_njy_arcsec2",
    "satellite_amplitude_at_50_per_galaxy_njy_arcsec2",
    "diffuse_scale_kpc",
    "satellite_scale_kpc",
    "n_injected_satellite_centers",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quantiles(values) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("summary values must be a non-empty finite vector")
    return {
        "median": float(np.median(values)),
        "q16": float(np.quantile(values, 0.16)),
        "q84": float(np.quantile(values, 0.84)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "q16_q84_contains_zero": bool(np.quantile(values, 0.16) <= 0 <= np.quantile(values, 0.84)),
    }


def load_batch(root: Path, expected_frames: int, kind: str, *, i_band_only=False) -> tuple[dict, dict]:
    root = Path(root).resolve()
    summary_path = root / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete"
            or summary.get("n_completed") != expected_frames
            or summary.get("n_errors") != 0
            or summary.get("science_measurement") is not False
            or summary.get("observational_validation_passed") is not False
            or summary.get("flux_covariance_available") is not False):
        raise ValueError(f"{kind} batch is incomplete or has unexpected science flags")
    records = {}
    for item in summary.get("records", []):
        path = Path(item["path"]).resolve()
        if sha256(path) != item.get("sha256"):
            raise ValueError(f"{kind} result checksum mismatch: {path}")
        record = json.loads(path.read_text(encoding="utf-8"))
        if (record.get("science_measurement") is not False
                or record.get("observational_validation_passed") is not False
                or record.get("flux_covariance_available") is not False):
            raise ValueError(f"{kind} frame has unexpected band or science flags")
        if i_band_only and str(record.get("filter")).lower() != "i":
            continue
        ident = str(record["cluster_id"])
        if ident in records:
            raise ValueError(f"duplicate {kind} cluster id {ident}")
        records[ident] = record
    return summary, records


def audit(core_dir: Path, wide_dir: Path) -> dict:
    core_summary, core = load_batch(core_dir, 57, "position-dependent core", i_band_only=True)
    wide_summary, wide = load_batch(wide_dir, 13, "public wide reference", i_band_only=True)
    if len(wide) != 13 or not set(wide).issubset(core):
        raise ValueError("wide-reference frames do not match a unique 13-frame core subset")
    ids = sorted(wide, key=int)
    for ident in ids:
        a, b = core[ident], wide[ident]
        if a.get("support") != b.get("support"):
            raise ValueError(f"matched-position support differs for cluster {ident}")
        for key in INJECTION_KEYS:
            if a["injection"].get(key) != b["injection"].get(key):
                raise ValueError(f"injection parameter {key} differs for cluster {ident}")

    by_variant, paired = {}, {}
    for metric in METRICS:
        core_values = np.asarray([core[ident]["transfer_response"][metric] for ident in ids])
        wide_values = np.asarray([wide[ident]["transfer_response"][metric] for ident in ids])
        by_variant[metric] = {
            "position_dependent_core": quantiles(core_values),
            "public_extended_reference": quantiles(wide_values),
        }
        paired[metric] = quantiles(wide_values - core_values)

    reference = wide_summary["reference"]
    return {
        "schema": "astrafek-extended-psf-core-vs-wide-audit-v1",
        "status": "complete_paired_sensitivity_audit",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "n_matched_i_band_frames": len(ids),
        "matched_cluster_ids": ids,
        "summary_sha256": {
            "position_dependent_core": sha256(Path(core_dir) / "summary.json"),
            "public_extended_reference": sha256(Path(wide_dir) / "summary.json"),
        },
        "public_reference": {
            "name": reference["name"],
            "sha256": reference["sha256"],
            "license": reference["license"],
            "psf_approved": reference["psf_approved"],
        },
        "response_units": "measured surface-brightness response per unit injected nJy/arcsec^2 amplitude",
        "response_summaries": by_variant,
        "paired_public_reference_minus_position_core": paired,
        "limitations": [
            "public wide reference is not matched to frame position, epoch, or seeing; differences are a model sensitivity, not a target-specific wing correction",
            "comparison covers only the 13 support-selected i-band frames; it does not test i2",
            "only synthetic light is convolved; injected-source photon noise, coadd covariance, and upstream sky transfer are omitted",
            "16th–84th percentiles describe frame-to-frame spread, not confidence intervals",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, default=DEFAULT_CORE)
    parser.add_argument("--wide-dir", type=Path, default=DEFAULT_WIDE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = audit(args.core_dir, args.wide_dir)
    result["audit_script_sha256"] = sha256(Path(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "matched_frames": result["n_matched_i_band_frames"],
        "science_measurement": False,
        "output": str(output),
        "output_sha256": sha256(output),
        "audit_script_sha256": result["audit_script_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
