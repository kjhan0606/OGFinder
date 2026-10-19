#!/usr/bin/env python3
"""Summarize position-dependent HSC PSF-core variation for the injection set.

All FITS files are checksum-verified against the approved 744-position audit.
Positive-flux moments and encircled-energy radii are descriptive core metrics;
the signed, unit-normalized PSF pixels remain unchanged in saved inputs and in
the existing convolution workflow.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = ROOT / "icl_psf_picker_injection_plan_2026-10-19_02"
DEFAULT_DOWNLOAD = ROOT / "icl_psf_picker_injection_2026-10-19_01"
DEFAULT_TEST = ROOT / "icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits"
TEST_SHA256 = "1276d46113d389116ded6f63bae4c9865bb42d24244b22247a84d2c6929cb202"
PIXEL_SCALE_ARCSEC = 0.168
GAUSSIAN_FWHM_FACTOR = 2.0 * math.sqrt(2.0 * math.log(2.0))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def coordinate_key(ra: float, dec: float) -> tuple[float, float]:
    return round(float(ra), 9), round(float(dec), 9)


def measure_core(psf: np.ndarray, pixel_scale_arcsec: float = PIXEL_SCALE_ARCSEC) -> dict:
    """Measure positive-flux core morphology without altering signed PSF data."""
    image = np.asarray(psf, dtype=np.float64)
    if (image.ndim != 2 or any(size % 2 != 1 for size in image.shape)
            or not np.isfinite(image).all()):
        raise ValueError("PSF must be a finite odd-sized 2-D image")
    signed_sum = float(image.sum())
    if not math.isfinite(signed_sum) or abs(signed_sum - 1.0) > 1e-3:
        raise ValueError("PSF signed pixel sum must be unit normalized")
    positive = np.clip(image, 0.0, None)
    positive_sum = float(positive.sum())
    if positive_sum <= 0 or not math.isfinite(positive_sum):
        raise ValueError("PSF has no finite positive flux")
    yy, xx = np.indices(image.shape, dtype=float)
    cx = float(np.sum(positive * xx) / positive_sum)
    cy = float(np.sum(positive * yy) / positive_sum)
    dx, dy = xx - cx, yy - cy
    cxx = float(np.sum(positive * dx * dx) / positive_sum)
    cyy = float(np.sum(positive * dy * dy) / positive_sum)
    cxy = float(np.sum(positive * dx * dy) / positive_sum)
    eigenvalues = np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]])
    if eigenvalues[0] <= 0 or not np.isfinite(eigenvalues).all():
        raise ValueError("positive-flux PSF covariance is not positive definite")
    sigma_minor, sigma_major = map(float, np.sqrt(eigenvalues))
    radius = np.hypot(dx, dy).ravel()
    order = np.argsort(radius, kind="mergesort")
    cumulative = np.cumsum(positive.ravel()[order]) / positive_sum

    def encircled_radius(fraction: float) -> float:
        index = int(np.searchsorted(cumulative, fraction, side="left"))
        return float(radius[order[min(index, len(order) - 1)]])

    return {
        "shape_yx": [int(image.shape[0]), int(image.shape[1])],
        "signed_sum": signed_sum,
        "positive_flux_fraction_over_signed_sum": positive_sum / signed_sum,
        "negative_pixel_count": int(np.count_nonzero(image < 0)),
        "centroid_xy_pixels": [cx, cy],
        "moment_fwhm_major_arcsec": GAUSSIAN_FWHM_FACTOR * sigma_major * pixel_scale_arcsec,
        "moment_fwhm_minor_arcsec": GAUSSIAN_FWHM_FACTOR * sigma_minor * pixel_scale_arcsec,
        "moment_axis_ratio_minor_over_major": sigma_minor / sigma_major,
        "encircled_energy_r50_arcsec": encircled_radius(0.50) * pixel_scale_arcsec,
        "encircled_energy_r80_arcsec": encircled_radius(0.80) * pixel_scale_arcsec,
        "encircled_energy_r90_arcsec": encircled_radius(0.90) * pixel_scale_arcsec,
    }


def distribution(records: list[dict], metric: str) -> dict:
    values = np.asarray([float(row[metric]) for row in records], dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError(f"no finite values for {metric}")
    return {
        "n": int(len(values)),
        "median": float(np.median(values)),
        "q16": float(np.quantile(values, 0.16)),
        "q84": float(np.quantile(values, 0.84)),
        "q05": float(np.quantile(values, 0.05)),
        "q95": float(np.quantile(values, 0.95)),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
    }


def read_psf(path: Path) -> np.ndarray:
    with fits.open(path, memmap=False, checksum=True) as hdus:
        if len(hdus) != 1 or hdus[0].data is None:
            raise ValueError(f"unexpected FITS structure: {path}")
        return np.asarray(hdus[0].data, dtype=np.float64)


def analyze(plan_dir: Path, download_dir: Path, test_psf: Path,
            output_path: Path) -> dict:
    plan_dir, download_dir, test_psf = map(Path, (plan_dir, download_dir, test_psf))
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    plan_path = plan_dir / "psf_picker_plan.json"
    mapping_path = plan_dir / "position_mapping.json"
    coordinate_path = plan_dir / "psf_picker_coordinates.txt"
    audit_path = download_dir / "psf_picker_audit.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-hsc-psf-picker-injection-plan-v1"
            or plan.get("user_approved_requests") != 744
            or plan.get("filter") != "HSC-I"):
        raise ValueError("unexpected or unapproved PSF request plan")
    if (audit.get("status") != "verified" or audit.get("n_requested") != 744
            or audit.get("n_returned") != 744 or audit.get("science_measurement") is not False):
        raise ValueError("744-position PSF batch audit is incomplete or science-promoted")
    plan_sha = sha256(plan_path)
    coordinate_sha = sha256(coordinate_path)
    mapping_sha = sha256(mapping_path)
    if (coordinate_sha != plan.get("coordinate_list_sha256")
            or mapping_sha != plan.get("position_mapping_sha256")
            or audit.get("coordinate_list_sha256") != coordinate_sha
            or audit.get("position_mapping_sha256") != mapping_sha):
        raise ValueError("PSF plan or audit checksum mismatch")

    expected = {}
    role_counts = {"cluster_center": 0, "satellite_injection_center": 0}
    for position in mapping:
        key = coordinate_key(position["ra_deg"], position["dec_deg"])
        if key in expected:
            raise ValueError(f"duplicate coordinate in approved mapping: {key}")
        roles = position.get("roles", [])
        if not roles:
            raise ValueError(f"approved coordinate has no role: {key}")
        for role in roles:
            name = str(role.get("role"))
            if name not in role_counts:
                raise ValueError(f"unexpected PSF request role: {name}")
            role_counts[name] += 1
        expected[key] = roles
    if len(expected) != 744 or role_counts != {
            "cluster_center": 56, "satellite_injection_center": 688}:
        raise ValueError(f"approved position-role counts changed: {role_counts}")

    measured = []
    seen = set()
    for record in audit["files"]:
        key = coordinate_key(record["ra_deg_filename"], record["dec_deg_filename"])
        if key not in expected or key in seen:
            raise ValueError(f"unexpected or repeated returned PSF coordinate: {key}")
        seen.add(key)
        path = download_dir / record["file"]
        digest = sha256(path)
        if digest != record.get("sha256"):
            raise ValueError(f"PSF file hash mismatch: {path.name}")
        metrics = measure_core(read_psf(path))
        roles = expected[key]
        metrics.update({
            "file": path.name,
            "sha256": digest,
            "ra_deg": key[0],
            "dec_deg": key[1],
            "roles": roles,
            "cluster_ids": sorted({str(role["cluster_id"]) for role in roles}),
            "role_names": sorted({str(role["role"]) for role in roles}),
            "source": "approved_744_position_batch",
        })
        measured.append(metrics)
    if len(seen) != 744 or seen != set(expected):
        raise ValueError("returned PSF set does not exactly match approved coordinates")

    if not test_psf.is_file() or sha256(test_psf) != TEST_SHA256:
        raise ValueError("single-position approved test PSF is missing or changed")
    test_metrics = measure_core(read_psf(test_psf))
    test_metrics.update({
        "file": test_psf.name,
        "sha256": sha256(test_psf),
        "ra_deg": 0.339748,
        "dec_deg": 1.030287,
        "roles": [{"cluster_id": "10", "frame_filter": "i2",
                   "role": "cluster_center", "target_index": None}],
        "cluster_ids": ["10"],
        "role_names": ["cluster_center"],
        "source": "approved_single_position_test",
    })
    if any(str(role["cluster_id"]) == "10" and role["role"] == "cluster_center"
           for roles in expected.values() for role in roles):
        raise ValueError("single-position test center unexpectedly duplicates a batch role")
    all_records = [*measured, test_metrics]
    metric_names = (
        "moment_fwhm_major_arcsec", "moment_fwhm_minor_arcsec",
        "moment_axis_ratio_minor_over_major", "encircled_energy_r50_arcsec",
        "encircled_energy_r80_arcsec", "encircled_energy_r90_arcsec",
        "positive_flux_fraction_over_signed_sum", "negative_pixel_count",
    )
    groups = {
        "all_745": all_records,
        "cluster_centers_57": [row for row in all_records if "cluster_center" in row["role_names"]],
        "satellite_positions_688": [row for row in all_records
                                     if "satellite_injection_center" in row["role_names"]],
        "batch_cluster_centers_56": [row for row in measured
                                      if "cluster_center" in row["role_names"]],
        "batch_satellite_positions_688": [row for row in measured
                                          if "satellite_injection_center" in row["role_names"]],
    }
    summaries = {
        group_name: {metric: distribution(group, metric) for metric in metric_names}
        for group_name, group in groups.items()
    }
    result = {
        "schema": "astrafek-hsc-position-psf-core-variation-v1",
        "status": "complete_exploratory_psf_core_audit",
        "science_measurement": False,
        "psf_approved": False,
        "n_batch_models": len(measured),
        "n_total_models_including_test": len(all_records),
        "request_role_counts_in_batch": role_counts,
        "filter": "HSC-I",
        "rerun": "pdr3_wide",
        "type": "coadd",
        "pixel_scale_arcsec": PIXEL_SCALE_ARCSEC,
        "metric_definition": {
            "moment_fwhm": "Gaussian-equivalent FWHM from second moments of the positive part of each signed PSF stamp; descriptive core metric, not a direct fitted seeing or wing width",
            "encircled_energy": "radius about the positive-flux centroid enclosing the stated fraction of positive stamp flux",
            "signed_psf_pixels": "original signed values are not clipped or renormalized in the FITS inputs or any convolution",
        },
        "group_summaries": summaries,
        "records": all_records,
        "inputs": {
            "plan": str(plan_path.resolve()), "plan_sha256": plan_sha,
            "mapping": str(mapping_path.resolve()), "mapping_sha256": mapping_sha,
            "coordinate_list": str(coordinate_path.resolve()), "coordinate_list_sha256": coordinate_sha,
            "batch_audit": str(audit_path.resolve()), "batch_audit_sha256": sha256(audit_path),
            "test_psf": str(test_psf.resolve()), "test_psf_sha256": TEST_SHA256,
        },
        "limitations": [
            "all PSFs are compact coadd model cores spanning only about 7 arcsec, not extended wings on ICL scales",
            "positive-only moment and encircled-energy metrics summarize the core and are not a substitute for signed-kernel convolution",
            "the 745 positions belong to the selected 57-frame injection sample, not all 800 clusters",
            "this audit does not validate sky transfer, observational noise, or scientific flux calibration",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--download-dir", type=Path, default=DEFAULT_DOWNLOAD)
    parser.add_argument("--test-psf", type=Path, default=DEFAULT_TEST)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.plan_dir, args.download_dir, args.test_psf, args.output)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "n_batch_models": result["n_batch_models"],
                      "n_total_models_including_test": result["n_total_models_including_test"],
                      "group_summaries": result["group_summaries"]}, indent=2))


if __name__ == "__main__":
    main()
