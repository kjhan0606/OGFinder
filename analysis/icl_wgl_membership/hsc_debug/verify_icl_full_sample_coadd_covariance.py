#!/usr/bin/env python3
"""Independently audit the saved full-sample HSC coadd covariance tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RUNNER = Path(__file__).with_name("audit_icl_full_sample_coadd_covariance.py")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fval(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("non-finite output value")
    return result


def close(a: float, b: float, atol: float = 1e-10, rtol: float = 1e-10) -> bool:
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def audit(work: Path) -> dict:
    work = work.resolve()
    summary_path = work / "summary.json"
    frame_path = work / "frame_noise_correlation.tsv"
    lag_path = work / "lag_correlation.tsv"
    launch_path = work / "launch_provenance.json"
    for path in (summary_path, frame_path, lag_path, launch_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    launch = json.loads(launch_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-icl-full-sample-coadd-residual-covariance-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("method", {}).get("image_content_rehashed") is not True
            or launch.get("slurm_job_id") != "448446"
            or launch.get("hsc_query_submitted") is not False
            or launch.get("new_data_downloaded") is not False
            or launch.get("runner_sha256_at_submission") != sha256(RUNNER)):
        raise ValueError("run status, provenance, or analysis-code hash is unexpected")
    frame_rows, lag_rows = rows(frame_path), rows(lag_path)
    if len(frame_rows) != 1566 or len(lag_rows) != 1566 * 81:
        raise ValueError("frame or 9x9 lag row count does not match the full sample")
    if (sha256(frame_path) != summary.get("output_sha256", {}).get(str(frame_path))
            or sha256(lag_path) != summary.get("output_sha256", {}).get(str(lag_path))):
        raise ValueError("output table SHA-256 mismatch")

    target_status_path = SAMPLE / "download_01/download_status.json"
    control_status_path = SAMPLE / "control_download_01/download_status.json"
    target_status = json.loads(target_status_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    target_files = {str(r["id"]): r for r in target_status["completed"]}
    control_files = {str(r["id"]): r for r in control_status["completed"]}
    lag_by_key: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in lag_rows:
        key = (str(row["cluster_id"]), str(row["role"]), str(row["filter"]))
        lag_by_key.setdefault(key, []).append(row)
    frame_by_key = {}
    for frame in frame_rows:
        key = (str(frame["cluster_id"]), str(frame["role"]), str(frame["filter"]))
        if key in frame_by_key:
            raise ValueError(f"duplicate frame row {key}")
        frame_by_key[key] = frame
        expected = target_files.get(key[0]) if key[1] == "target" else control_files.get(str(frame["image_id"]))
        if expected is None or frame.get("image_sha256") != expected.get("sha256"):
            raise ValueError(f"image hash does not match the acquisition ledger for {key}")
        path = (SAMPLE / "download_01/images" / str(expected["file"])) if key[1] == "target" else (
            SAMPLE / "control_download_01/images" / str(expected["file"]))
        if not path.is_file() or path.stat().st_size != expected.get("bytes"):
            raise ValueError(f"image path/size differs from its acquisition ledger for {key}")
        kernel_rows = lag_by_key.get(key, [])
        lag_values = {(int(r["dy"]), int(r["dx"])): r for r in kernel_rows}
        expected_lags = {(dy, dx) for dy in range(-4, 5) for dx in range(-4, 5)}
        if set(lag_values) != expected_lags:
            raise ValueError(f"lag grid incomplete for {key}")
        factor = sum(fval(row["median_variance_normalized_covariance"])
                     for row in kernel_rows if row["median_variance_normalized_covariance"] != "")
        center = fval(lag_values[(0, 0)]["median_variance_normalized_covariance"])
        if (not close(factor, fval(frame["top_hat_variance_inflation_approx_square_lag_window"]))
                or not close(center, fval(frame["variance_plane_ratio_at_zero_lag"]))
                or not close(math.sqrt(factor), fval(frame["top_hat_sigma_inflation_approx_square_lag_window"]))):
            raise ValueError(f"saved covariance summary arithmetic does not reproduce for {key}")
    if set(lag_by_key) != set(frame_by_key):
        raise ValueError("lag-table frame keys differ from the frame summary")

    by_filter = {}
    for band, expected_n in (("i", 302), ("i2", 481)):
        target = {key[0]: frame for key, frame in frame_by_key.items()
                  if key[1:] == ("target", band)}
        control = {key[0]: frame for key, frame in frame_by_key.items()
                   if key[1:] == ("matched_control", band)}
        if len(target) != expected_n or set(target) != set(control):
            raise ValueError(f"unexpected paired-frame support in {band}")
        by_filter[band] = {}
        for role, data in (("target", list(target.values())),
                           ("matched_control", list(control.values()))):
            factors = np.asarray([fval(r["top_hat_variance_inflation_approx_square_lag_window"])
                                  for r in data])
            center = np.asarray([fval(r["variance_plane_ratio_at_zero_lag"]) for r in data])
            actual = summary["by_filter_and_role"][band][role]
            if (not close(float(np.median(factors)), float(actual["median_variance_factor"]))
                    or not np.allclose(np.quantile(factors, [0.16, 0.84]), actual["q16_q84_variance_factor"],
                                       rtol=1e-10, atol=1e-10)
                    or not close(float(np.median(center)), float(actual["median_zero_lag_variance_ratio"]))):
                raise ValueError(f"filter summary does not reproduce for {band}/{role}")
        delta = np.asarray([
            fval(target[ident]["top_hat_variance_inflation_approx_square_lag_window"])
            - fval(control[ident]["top_hat_variance_inflation_approx_square_lag_window"])
            for ident in target], dtype=float)
        saved = summary["by_filter_and_role"][band]["paired_target_minus_control"]
        if (not close(float(np.median(delta)), float(saved["median"]))
                or not np.allclose(np.quantile(delta, [0.16, 0.84]), saved["q16_q84"],
                                   rtol=1e-10, atol=1e-10)):
            raise ValueError(f"paired target-control summary does not reproduce for {band}")
        by_filter[band] = {"pairs": expected_n, "target_control_ids_match": True,
                           "frame_and_pair_summaries_reproduced": True}

    result = {
        "schema": "astrafek-icl-full-sample-coadd-covariance-independent-audit-v1",
        "status": "verified", "science_measurement": False,
        "n_frame_rows": len(frame_rows), "n_lag_rows": len(lag_rows),
        "by_filter": by_filter,
        "frame_and_lag_arithmetic_reproduced": True,
        "input_acquisition_hashes_matched": True,
        "frame_table_sha256": sha256(frame_path), "lag_table_sha256": sha256(lag_path),
        "summary_sha256": sha256(summary_path),
        "limitations": [
            "The finite lag window is a local covariance diagnostic, not the exact radial-estimator covariance.",
            "This table-level validation does not independently remeasure image pixels.",
            "The result remains exploratory and is not a correction to spatial-block errors.",
        ],
    }
    output = work / "independent_audit.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result["audit_sha256"] = sha256(output)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    audit(args.work)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
