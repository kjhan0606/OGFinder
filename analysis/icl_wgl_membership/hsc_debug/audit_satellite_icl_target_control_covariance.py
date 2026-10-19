#!/usr/bin/env python3
"""Add matched-control ICL aperture noise to the shared-pixel covariance audit.

This CPU analysis uses only already downloaded HSC target/control coadds. It
estimates each side's 200--300 kpc aperture variance from its variance plane
and matched 9x9 residual-covariance kernel, then combines independent target
and control terms with the previously measured target/satellite covariance.
It is an exploratory noise diagnostic, not a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy import units as u

from audit_satellite_icl_shared_pixel_covariance import lag_weighted_covariance
from estimate_hsc_coadd_noise_correlation import (
    _matched_control_source_mask, load_matched_controls,
)
from validate_control_images import nJy_factor


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_TARGET_COVARIANCE = ROOT / "satellite_icl_lag_covariance_2026-10-19_01"
DEFAULT_LAG_DIR = ROOT / "icl_coadd_noise_correlation_2026-10-19_04"
DEFAULT_CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
DEFAULT_CONTROL_STATUS = ROOT / "control_download_01/download_status.json"
DEFAULT_METADATA = ROOT / "icl_coadd_noise_metadata_2026-10-19_01/summary.json"
DEFAULT_OUTPUT = ROOT / "satellite_icl_target_control_covariance_2026-10-19_01"
INNER_KPC, OUTER_KPC, COVERAGE_FLOOR = 200.0, 300.0, 0.40


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def kernel_bank(lag_dir: Path, role: str) -> tuple[dict, dict]:
    if role not in {"target", "matched_control"}:
        raise ValueError("role must be target or matched_control")
    summary_path, table_path = lag_dir / "summary.json", lag_dir / "lag_correlation.tsv"
    summary = json.loads(summary_path.read_text())
    if (summary.get("schema") != "astrafek-hsc-coadd-noise-correlation-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("n_target_frames") != 57
            or summary.get("n_matched_control_frames") != 57):
        raise ValueError("expected the audited matched 57-frame target/control lag-kernel run")
    if sha256(table_path) != summary.get("lag_results_sha256"):
        raise ValueError("lag-correlation table checksum mismatch")
    for raw_path, expected in summary.get("input_sha256", {}).items():
        if sha256(Path(raw_path)) != expected:
            raise ValueError(f"lag-kernel source checksum mismatch: {raw_path}")
    source_script = Path(__file__).with_name("estimate_hsc_coadd_noise_correlation.py")
    if sha256(source_script) != summary.get("analysis_script_sha256"):
        raise ValueError("lag-kernel source code differs from its audited version")

    width = int(summary["max_lag_pixels"])
    grouped: dict[tuple[str, str], dict[tuple[int, int], float]] = {}
    for row in read_tsv(table_path):
        if row["role"] != role:
            continue
        key = str(row["cluster_id"]), str(row["filter"])
        lag = int(row["dy"]), int(row["dx"])
        value = float(row["median_variance_normalized_covariance"])
        if not math.isfinite(value) or lag in grouped.setdefault(key, {}):
            raise ValueError(f"invalid or duplicate {role} lag for {key}: {lag}")
        grouped[key][lag] = value
    if len(grouped) != 57:
        raise ValueError(f"expected 57 {role} lag kernels, got {len(grouped)}")
    expected_lags = {(dy, dx) for dy in range(-width, width + 1)
                     for dx in range(-width, width + 1)}
    kernels = {}
    for key, values in grouped.items():
        if set(values) != expected_lags:
            raise ValueError(f"incomplete {role} lag grid for {key}")
        kernel = np.empty((2 * width + 1, 2 * width + 1), dtype=np.float64)
        for dy, dx in expected_lags:
            kernel[dy + width, dx + width] = 0.5 * (
                values[(dy, dx)] + values[(-dy, -dx)])
        kernels[key] = kernel
    return kernels, {
        "summary_path": str(summary_path.resolve()), "summary_sha256": sha256(summary_path),
        "lag_table_path": str(table_path.resolve()), "lag_table_sha256": sha256(table_path),
        "role": role,
    }


def annulus_weights(shape: tuple[int, int], center_xy: tuple[float, float],
                    kpc_per_pixel: float, valid: np.ndarray,
                    inner_kpc: float = INNER_KPC, outer_kpc: float = OUTER_KPC
                    ) -> tuple[np.ndarray, float, int]:
    """Return uniform valid-pixel weights and geometric annulus coverage."""
    if valid.shape != shape or valid.ndim != 2 or not np.isfinite(valid).all():
        raise ValueError("valid-pixel mask must be a finite 2D image mask")
    if not (math.isfinite(kpc_per_pixel) and kpc_per_pixel > 0
            and math.isfinite(inner_kpc) and math.isfinite(outer_kpc)
            and 0 <= inner_kpc < outer_kpc):
        raise ValueError("invalid pixel scale or annulus radii")
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    cx, cy = map(float, center_xy)
    radius = np.hypot((xx - cx) * kpc_per_pixel,
                      (yy - cy) * kpc_per_pixel)
    annulus = (radius >= inner_kpc) & (radius < outer_kpc)
    geometric = int(annulus.sum())
    selected = annulus & valid
    count = int(selected.sum())
    if geometric == 0 or count == 0:
        raise ValueError("ICL annulus has no valid pixels")
    weights = np.zeros(shape, dtype=np.float64)
    weights[selected] = 1.0 / count
    return weights, count / geometric, count


def combine_variances(target_diag: float, target_lag: float,
                      control_diag: float, control_lag: float,
                      satellite_diag: float, satellite_lag: float,
                      covariance_diag: float, covariance_lag: float) -> dict:
    values = (target_diag, target_lag, control_diag, control_lag,
              satellite_diag, satellite_lag, covariance_diag, covariance_lag)
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError("all covariance terms must be finite")
    if min(target_diag, target_lag, control_diag, control_lag,
           satellite_diag, satellite_lag) <= 0:
        raise ValueError("all diagonal and lag variances must be positive")
    icl_diag, icl_lag = target_diag + control_diag, target_lag + control_lag
    return {
        "icl_delta_diagonal_variance": icl_diag,
        "icl_delta_lag_variance": icl_lag,
        "icl_delta_sigma_inflation": math.sqrt(icl_lag / icl_diag),
        "control_to_target_lag_variance_ratio": control_lag / target_lag,
        "icl_satellite_diagonal_noise_correlation": covariance_diag / math.sqrt(icl_diag * satellite_diag),
        "icl_satellite_lag_noise_correlation": covariance_lag / math.sqrt(icl_lag * satellite_lag),
    }


def diagonal_weighted_variance(weights: np.ndarray, variance: np.ndarray,
                               conversion: float) -> float:
    """Compute diagonal aperture variance without evaluating 0*inf pixels."""
    weights, variance = np.asarray(weights, dtype=np.float64), np.asarray(variance, dtype=np.float64)
    if weights.shape != variance.shape or weights.ndim != 2:
        raise ValueError("weights and variance must be matching 2D arrays")
    selected = weights != 0
    if (not np.any(selected) or not np.isfinite(weights[selected]).all()
            or not np.isfinite(variance[selected]).all() or np.any(variance[selected] <= 0)
            or not math.isfinite(float(conversion)) or conversion <= 0):
        raise ValueError("weighted pixels and conversion must be finite with positive variance")
    return float(np.sum(weights[selected]**2 * variance[selected], dtype=np.float64)
                 * float(conversion)**2)


def run(target_covariance_dir: Path, lag_dir: Path, metadata_path: Path,
        control_manifest_path: Path, control_status_path: Path,
        output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    target_summary_path = target_covariance_dir / "summary.json"
    target_table_path = target_covariance_dir / "pair_pixel_covariance.tsv"
    target_summary = json.loads(target_summary_path.read_text())
    if (target_summary.get("schema") != "hsc-icl-satellite-shared-pixel-lag-covariance-v1"
            or target_summary.get("status") != "complete_exploratory"
            or target_summary.get("science_measurement") is not False
            or sha256(target_table_path) != target_summary.get("output_sha256", {}).get(str(target_table_path))):
        raise ValueError("target/satellite covariance input is missing, changed, or not exploratory")
    for raw_path, expected in target_summary.get("input_sha256", {}).items():
        if sha256(Path(raw_path)) != expected:
            raise ValueError(f"target/satellite covariance source checksum mismatch: {raw_path}")
    target_rows = read_tsv(target_table_path)
    if len(target_rows) != 57 or len({(r["cluster_id"], r["filter"]) for r in target_rows}) != 57:
        raise ValueError("expected exactly 57 unique target/satellite covariance pairs")
    target_summary_sha = sha256(target_summary_path)

    metadata = json.loads(metadata_path.read_text())
    if (metadata.get("schema") != "astrafek-icl-coadd-noise-metadata-audit-v1"
            or metadata.get("status") != "verified_metadata_only"
            or metadata.get("science_measurement") is not False
            or metadata.get("n_frames") != 57):
        raise ValueError("expected the checksum-audited 57-frame target metadata")
    controls = load_matched_controls(metadata, control_manifest_path, control_status_path)
    control_by_key = {(r["cluster_id"], r["filter"]): r for r in controls}
    if len(control_by_key) != 57:
        raise ValueError("matched-control image mapping is not one-to-one")
    control_manifest = json.loads(control_manifest_path.read_text())
    control_rows = {(str(r["target_cluster_id"]), str(r["coadd_header_filter"])): r
                    for r in control_manifest["clusters"]}
    if len(control_rows) != len(control_manifest["clusters"]):
        raise ValueError("matched-control manifest has duplicate target/filter mappings")
    control_status = json.loads(control_status_path.read_text())
    completed = {str(r["id"]): r for r in control_status["completed"]}
    sample_manifest_path = ROOT / "sample_manifest.json"
    sample_manifest = json.loads(sample_manifest_path.read_text())
    frames = {str(r["id"]): r for r in sample_manifest["clusters"]}

    target_kernels, target_kernel_provenance = kernel_bank(lag_dir, "target")
    control_kernels, control_kernel_provenance = kernel_bank(lag_dir, "matched_control")
    if set(target_kernels) != set(control_kernels):
        raise ValueError("target and matched-control lag-kernel pairs differ")

    lag_rows = read_tsv(lag_dir / "lag_correlation.tsv")
    max_lag = int(json.loads((lag_dir / "summary.json").read_text())["max_lag_pixels"])
    expected_kernel_rows = (2 * max_lag + 1) ** 2
    control_patch_counts = {}
    for ident, filt in control_kernels:
        counts = [int(r["n_patches"]) for r in lag_rows
                  if r["role"] == "matched_control" and r["cluster_id"] == ident
                  and r["filter"] == filt]
        if len(counts) != expected_kernel_rows:
            raise ValueError(f"incomplete matched-control patch counts for {ident}/{filt}")
        control_patch_counts[(ident, filt)] = float(np.median(counts))

    h0 = float(sample_manifest.get("h0", 70.0))
    omega_m = float(sample_manifest.get("omega_m", 0.3))
    cosmo = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0 * u.K)
    records, errors = [], []
    for index, row in enumerate(sorted(target_rows, key=lambda x: (int(x["cluster_id"]), x["filter"])), 1):
        ident, filt = str(row["cluster_id"]), str(row["filter"])
        key = ident, filt
        print(f"[{index}/57] target+control ICL aperture noise {ident} {filt}", flush=True)
        try:
            control = control_by_key.get(key)
            planned = control_rows.get(key)
            if control is None or planned is None:
                raise ValueError("no unique matched control image for target/filter")
            saved = completed[str(planned["id"])]
            review = saved["review"]
            if sha256(Path(control["path"])) != control["source_image_sha256"]:
                raise ValueError("matched-control FITS SHA-256 changed")
            frame = frames[ident]
            if str(planned["target_cluster_id"]) != ident or not math.isclose(
                    float(planned["z"]), float(frame["z"]), rel_tol=0, abs_tol=1e-12):
                raise ValueError("matched-control target identity/redshift mismatch")
            if str(planned["coadd_header_filter"]) != filt:
                raise ValueError("matched-control coadd filter differs from target")

            with fits.open(control["path"], memmap=True, checksum=True) as hdus:
                image = np.asarray(hdus[1].data, dtype=np.float32)
                image_mask = np.asarray(hdus[2].data, dtype=np.int64)
                variance = np.asarray(hdus[3].data, dtype=np.float32)
                primary = hdus[0].header
            if image.shape != variance.shape or image.shape != image_mask.shape:
                raise ValueError("matched-control IMAGE/MASK/VARIANCE dimensions differ")
            if str(primary.get("FILTER", "")).strip().lower() != filt:
                raise ValueError("matched-control FITS FILTER header mismatch")
            if not math.isclose(float(review["pixel_scale_arcsec"]),
                                float(control["pixel_scale_arcsec"]), rel_tol=0, abs_tol=1e-9):
                raise ValueError("matched-control pixel-scale review mismatch")
            valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0)
            excluded = _matched_control_source_mask(
                image, variance, image_mask, float(review["pixel_scale_arcsec"]))
            if excluded.shape != image.shape:
                raise ValueError("derived matched-control source mask shape mismatch")
            valid &= ~excluded
            scale = float(review["pixel_scale_arcsec"])
            da_mpc = float(cosmo.angular_diameter_distance(float(frame["z"])).value)
            kpc_per_pixel = float(cosmo.kpc_proper_per_arcmin(float(frame["z"])).value * scale / 60.)
            dra = math.radians(float(planned["ra_deg"]) - float(frame["ra_deg"]))
            dec1, dec2 = math.radians(float(frame["dec_deg"])), math.radians(float(planned["dec_deg"]))
            hav = (math.sin((dec2-dec1)/2)**2
                   + math.cos(dec1)*math.cos(dec2)*math.sin(dra/2)**2)
            separation_arcsec = 2 * math.asin(min(1.0, math.sqrt(max(0.0, hav)))) * 206264.806247
            nonoverlap_min_arcsec = math.degrees(2 * OUTER_KPC / (da_mpc * 1000.0)) * 3600.0
            if separation_arcsec <= nonoverlap_min_arcsec:
                raise ValueError("target/control 200--300 kpc annuli overlap on the sky")
            weights, coverage, n_pixels = annulus_weights(
                image.shape, tuple(review["center_pixel_xy"]), kpc_per_pixel, valid)
            if coverage < COVERAGE_FLOOR:
                raise ValueError(f"matched-control ICL annulus coverage {coverage:.4f} < {COVERAGE_FLOOR}")
            conversion = float(nJy_factor({
                "input_unit": review["input_unit"], "fluxmag0": review["fluxmag0"],
                "pixel_scale_arcsec": scale,
            }))
            kernel = control_kernels[key]
            control_diag = diagonal_weighted_variance(weights, variance, conversion)
            control_lag = lag_weighted_covariance(
                weights, weights, variance, kernel, conversion)
            if control_diag <= 0 or control_lag <= 0:
                raise ValueError("matched-control lag model returned nonpositive aperture variance")

            target_diag = float(row["diagonal_var_icl_target_njy2_arcsec4"])
            target_lag = float(row["lag_kernel_var_icl_target_njy2_arcsec4"])
            sat_diag = float(row["diagonal_var_satellite_unweighted_njy2_arcsec4"])
            sat_lag = float(row["lag_kernel_var_satellite_unweighted_njy2_arcsec4"])
            cov_diag = float(row["diagonal_covariance_njy2_arcsec4"])
            cov_lag = float(row["lag_kernel_covariance_njy2_arcsec4"])
            combined = combine_variances(target_diag, target_lag, control_diag, control_lag,
                                         sat_diag, sat_lag, cov_diag, cov_lag)
            control_record = completed[str(planned["id"])]
            records.append({
                "cluster_id": ident, "filter": filt, "control_id": str(planned["id"]),
                "control_image_sha256": str(control["source_image_sha256"]),
                "control_source_mask": "matched SEP/Kron-2.5 mask regenerated from downloaded control coadd",
                "target_ra_deg": float(frame["ra_deg"]), "target_dec_deg": float(frame["dec_deg"]),
                "control_ra_deg": float(planned["ra_deg"]), "control_dec_deg": float(planned["dec_deg"]),
                "target_control_center_separation_arcsec": separation_arcsec,
                "minimum_separation_for_disjoint_300kpc_annuli_arcsec": nonoverlap_min_arcsec,
                "target_control_annuli_nonoverlapping": True,
                "control_center_pixel_x": float(review["center_pixel_xy"][0]),
                "control_center_pixel_y": float(review["center_pixel_xy"][1]),
                "control_icl_coverage_200_300": coverage, "control_icl_valid_pixels": n_pixels,
                "control_pixel_scale_arcsec": scale, "control_ab_zeropoint": float(review["ab_zeropoint"]),
                "control_diagonal_variance_njy2_arcsec4": control_diag,
                "control_lag_variance_njy2_arcsec4": control_lag,
                "control_lag_sigma_inflation": math.sqrt(control_lag / control_diag),
                "control_lag_kernel_n_patches_median": control_patch_counts[key],
                "control_download_bytes": int(control_record["bytes"]),
                **combined,
                "target_icl_diagonal_variance_njy2_arcsec4": target_diag,
                "target_icl_lag_variance_njy2_arcsec4": target_lag,
                "target_satellite_diagonal_variance_njy2_arcsec4": sat_diag,
                "target_satellite_lag_variance_njy2_arcsec4": sat_lag,
                "target_icl_satellite_diagonal_covariance_njy2_arcsec4": cov_diag,
                "target_icl_satellite_lag_covariance_njy2_arcsec4": cov_lag,
                "assumed_target_control_noise_independent": True,
            })
        except Exception as exc:
            errors.append({"cluster_id": ident, "filter": filt,
                           "error_type": type(exc).__name__, "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)

    if not records:
        raise RuntimeError("no matched target/control aperture covariance rows completed")
    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "target_control_pixel_covariance.tsv"
    with table_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(records)
    by_filter = {}
    for filt in sorted({r["filter"] for r in records}):
        subset = [r for r in records if r["filter"] == filt]
        by_filter[filt] = {
            "n_pairs": len(subset),
            "median_control_to_target_lag_variance_ratio": float(np.median([
                r["control_to_target_lag_variance_ratio"] for r in subset])),
            "median_total_icl_sigma_inflation_vs_diagonal": float(np.median([
                r["icl_delta_sigma_inflation"] for r in subset])),
            "median_diagonal_icl_satellite_noise_correlation": float(np.median([
                r["icl_satellite_diagonal_noise_correlation"] for r in subset])),
            "median_lag_icl_satellite_noise_correlation": float(np.median([
                r["icl_satellite_lag_noise_correlation"] for r in subset])),
            "q16_q84_lag_icl_satellite_noise_correlation": np.quantile([
                r["icl_satellite_lag_noise_correlation"] for r in subset], [.16, .84]).tolist(),
            "control_coverage_q16_q84": np.quantile([
                r["control_icl_coverage_200_300"] for r in subset], [.16, .84]).tolist(),
        }
    source_paths = [target_summary_path, target_table_path, metadata_path,
                    control_manifest_path, control_status_path,
                    sample_manifest_path, Path(__file__).resolve(),
                    Path(__file__).with_name("estimate_hsc_coadd_noise_correlation.py"),
                    Path(__file__).with_name("audit_satellite_icl_shared_pixel_covariance.py"),
                    Path(__file__).with_name("validate_control_images.py")]
    summary = {
        "schema": "hsc-icl-satellite-target-control-pixel-covariance-v1",
        "status": "complete_exploratory" if not errors else "partial_exploratory",
        "science_measurement": False,
        "target_covariance_summary_sha256": target_summary_sha,
        "target_covariance_table_sha256": sha256(target_table_path),
        "lag_kernel_provenance": [target_kernel_provenance, control_kernel_provenance],
        "selection": {"n_requested_pairs": len(target_rows), "n_completed_pairs": len(records),
                      "n_errors": len(errors), "icl_annulus_kpc": [INNER_KPC, OUTER_KPC],
                      "control_coverage_floor": COVERAGE_FLOOR},
        "assumptions": [
            "Target and matched-control coadd noise are treated as independent; their verified sky centers and non-overlapping annuli are recorded in the row table.",
            "The ICL delta estimator variance is target plus matched-control aperture variance.",
            "Cross covariance with the satellite estimator comes from its shared target-image pixels; matched-control image cross covariance is assumed zero.",
            "Target and control apertures use their saved/regenerated K2.5 masks and delivered variance planes.",
        ],
        "by_filter": by_filter, "errors": errors,
        "input_sha256": {str(p): sha256(p) for p in source_paths},
        "output_sha256": {str(table_path): sha256(table_path)},
        "limitations": [
            "The residual covariance kernel is stationary and measured on masked patches; it is not the exact position-dependent covariance of each aperture.",
            "This calculation excludes source Poisson noise, sky-transfer uncertainty, PSF and photometric-calibration uncertainty, extinction, and astrophysical covariance.",
            "The ICL–satellite covariance is conditional on the existing detections, masks, estimator supports, and equal weighting.",
            "A per-frame noise correlation is not a population-level association test and must not be interpreted as a science detection.",
        ],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-covariance", type=Path, default=DEFAULT_TARGET_COVARIANCE)
    parser.add_argument("--lag-dir", type=Path, default=DEFAULT_LAG_DIR)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--control-manifest", type=Path, default=DEFAULT_CONTROL_MANIFEST)
    parser.add_argument("--control-status", type=Path, default=DEFAULT_CONTROL_STATUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run(args.target_covariance.resolve(), args.lag_dir.resolve(), args.metadata.resolve(),
                   args.control_manifest.resolve(), args.control_status.resolve(), args.output.resolve())
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "selection": summary["selection"], "by_filter": summary["by_filter"],
                      "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    main()
