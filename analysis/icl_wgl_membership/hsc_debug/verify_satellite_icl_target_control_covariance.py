#!/usr/bin/env python3
"""Independently verify target/control ICL covariance outputs and provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_RUN = ROOT / "satellite_icl_target_control_covariance_2026-10-19_03"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def close(a: float, b: float) -> bool:
    return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-14)


def max_lag_distance_arcsec(kernel_shape: tuple[int, int], pixel_scale: float) -> float:
    """Largest pixel-center separation represented by a finite lag kernel."""
    if (len(kernel_shape) != 2 or any(size < 1 or size % 2 == 0 for size in kernel_shape)
            or not math.isfinite(pixel_scale) or pixel_scale <= 0):
        raise ValueError("an odd 2-D kernel and positive finite pixel scale are required")
    half_y, half_x = ((size - 1) / 2 for size in kernel_shape)
    return math.hypot(half_x, half_y) * pixel_scale


def verify(run_dir: Path) -> dict:
    summary_path = run_dir / "summary.json"
    table_path = run_dir / "target_control_pixel_covariance.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "hsc-icl-satellite-target-control-pixel-covariance-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("selection", {}).get("n_requested_pairs") != 57
            or summary.get("selection", {}).get("n_completed_pairs") != 57
            or summary.get("selection", {}).get("n_errors") != 0):
        raise ValueError("run summary is incomplete or has an unexpected science status")
    if sha256(table_path) != summary.get("output_sha256", {}).get(str(table_path)):
        raise ValueError("output table checksum differs from the summary")
    for raw_path, expected in summary.get("input_sha256", {}).items():
        if sha256(Path(raw_path)) != expected:
            raise ValueError(f"source checksum mismatch: {raw_path}")

    target_summary_path = ROOT / "satellite_icl_lag_covariance_2026-10-19_01/summary.json"
    target_table_path = ROOT / "satellite_icl_lag_covariance_2026-10-19_01/pair_pixel_covariance.tsv"
    target_summary = json.loads(target_summary_path.read_text(encoding="utf-8"))
    if sha256(target_summary_path) != summary.get("target_covariance_summary_sha256"):
        raise ValueError("target covariance summary checksum mismatch")
    if sha256(target_table_path) != summary.get("target_covariance_table_sha256"):
        raise ValueError("target covariance table checksum mismatch")
    target_rows = {(r["cluster_id"], r["filter"]): r for r in read_tsv(target_table_path)}
    rows = read_tsv(table_path)
    keys = [(r["cluster_id"], r["filter"]) for r in rows]
    if len(rows) != 57 or len(set(keys)) != 57 or set(keys) != set(target_rows):
        raise ValueError("combined rows do not match the 57 unique target covariance pairs")
    if sum(r["filter"] == "i" for r in rows) != 13 or sum(r["filter"] == "i2" for r in rows) != 44:
        raise ValueError("filter counts differ from the audited support sample")

    control_manifest_path = ROOT / "control_download_01/download_manifest.json"
    control_status_path = ROOT / "control_download_01/download_status.json"
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    planned = {(str(r["target_cluster_id"]), str(r["coadd_header_filter"])): r
               for r in control_manifest["clusters"]}
    completed = {str(r["id"]): r for r in control_status["completed"]}
    short_range_gaps = []
    short_range_supports = []
    for row in rows:
        key = row["cluster_id"], row["filter"]
        plan = planned.get(key)
        if plan is None or str(plan["id"]) != row["control_id"]:
            raise ValueError(f"wrong matched-control mapping for {key}")
        saved = completed.get(row["control_id"])
        if (saved is None or saved.get("file") != plan.get("image")
                or saved.get("sha256") != row["control_image_sha256"]):
            raise ValueError(f"control image provenance mismatch for {key}")
        if (float(row["control_icl_coverage_200_300"]) < 0.4
                or row["target_control_annuli_nonoverlapping"] != "True"
                or float(row["target_control_center_separation_arcsec"])
                <= float(row["minimum_separation_for_disjoint_300kpc_annuli_arcsec"])):
            raise ValueError(f"control annulus fails support/separation gate for {key}")
        gap = (float(row["target_control_center_separation_arcsec"])
               - float(row["minimum_separation_for_disjoint_300kpc_annuli_arcsec"]))
        kernel_support = max_lag_distance_arcsec(
            (9, 9), float(row["control_pixel_scale_arcsec"]))
        if gap <= kernel_support:
            raise ValueError(f"target/control apertures approach within 9x9 lag support for {key}")
        short_range_gaps.append(gap)
        short_range_supports.append(kernel_support)

        target = target_rows[key]
        target_diag = float(target["diagonal_var_icl_target_njy2_arcsec4"])
        target_lag = float(target["lag_kernel_var_icl_target_njy2_arcsec4"])
        sat_diag = float(target["diagonal_var_satellite_unweighted_njy2_arcsec4"])
        sat_lag = float(target["lag_kernel_var_satellite_unweighted_njy2_arcsec4"])
        cov_diag = float(target["diagonal_covariance_njy2_arcsec4"])
        cov_lag = float(target["lag_kernel_covariance_njy2_arcsec4"])
        ctrl_diag = float(row["control_diagonal_variance_njy2_arcsec4"])
        ctrl_lag = float(row["control_lag_variance_njy2_arcsec4"])
        combined_diag, combined_lag = target_diag + ctrl_diag, target_lag + ctrl_lag
        expected = {
            "icl_delta_diagonal_variance": combined_diag,
            "icl_delta_lag_variance": combined_lag,
            "icl_delta_sigma_inflation": math.sqrt(combined_lag / combined_diag),
            "control_to_target_lag_variance_ratio": ctrl_lag / target_lag,
            "icl_satellite_diagonal_noise_correlation": cov_diag / math.sqrt(combined_diag * sat_diag),
            "icl_satellite_lag_noise_correlation": cov_lag / math.sqrt(combined_lag * sat_lag),
        }
        for name, value in expected.items():
            if not close(row[name], value):
                raise ValueError(f"combined formula mismatch for {key}/{name}")

    summary_checks = {}
    for filt, claimed in summary["by_filter"].items():
        subset = [r for r in rows if r["filter"] == filt]
        if len(subset) != claimed["n_pairs"]:
            raise ValueError(f"filter row count mismatch for {filt}")
        columns = {
            "median_control_to_target_lag_variance_ratio": "control_to_target_lag_variance_ratio",
            "median_total_icl_sigma_inflation_vs_diagonal": "icl_delta_sigma_inflation",
            "median_diagonal_icl_satellite_noise_correlation": "icl_satellite_diagonal_noise_correlation",
            "median_lag_icl_satellite_noise_correlation": "icl_satellite_lag_noise_correlation",
        }
        recomputed = {}
        for output_name, column in columns.items():
            value = float(np.median([float(r[column]) for r in subset]))
            if not close(value, claimed[output_name]):
                raise ValueError(f"summary median mismatch for {filt}/{output_name}")
            recomputed[output_name] = value
        q16_q84 = np.quantile([float(r["icl_satellite_lag_noise_correlation"])
                               for r in subset], [0.16, 0.84]).tolist()
        if not np.allclose(q16_q84, claimed["q16_q84_lag_icl_satellite_noise_correlation"],
                           rtol=1e-12, atol=1e-14):
            raise ValueError(f"summary quantiles mismatch for {filt}")
        summary_checks[filt] = {**recomputed, "q16_q84_lag_icl_satellite_noise_correlation": q16_q84,
                                "n_pairs": len(subset)}
    return {
        "schema": "astrafek-icl-satellite-target-control-covariance-verification-v1",
        "status": "verified", "science_measurement": False,
        "run_summary_sha256": sha256(summary_path), "run_table_sha256": sha256(table_path),
        "n_source_hashes_verified": len(summary["input_sha256"]),
        "n_rows_formula_verified": len(rows), "n_control_mappings_verified": len(rows),
        "n_disjoint_covered_apertures_verified": len(rows),
        "finite_lag_target_control_separation": {
            "kernel_shape_pixels": [9, 9],
            "pairs_outside_finite_lag_covariance_support": len(short_range_gaps),
            "minimum_annulus_boundary_gap_arcsec": min(short_range_gaps),
            "median_annulus_boundary_gap_arcsec": float(np.median(short_range_gaps)),
            "maximum_9x9_kernel_pixel_center_separation_arcsec": max(short_range_supports),
            "minimum_gap_to_kernel_support_ratio": min(
                gap / support for gap, support in zip(short_range_gaps, short_range_supports)),
            "scope": "excludes cross-aperture pixel covariance within the measured 9x9 finite-lag kernel only; large-scale sky/calibration covariance remains unmeasured",
        },
        "filter_summary_recomputed": summary_checks,
        "verifier_sha256": sha256(Path(__file__).resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--output", type=Path,
                        help="optional new JSON file; defaults to stdout only")
    args = parser.parse_args()
    result = verify(args.run_dir.resolve())
    if args.output:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
