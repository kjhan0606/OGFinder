#!/usr/bin/env python3
"""Audit the saved 5-degree/2.5-degree ICL-split WL sensitivity summaries."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_5DEG = DATA / "wl_random_subtracted_icl_split_2026-11-05_02/summary.json"
DEFAULT_2P5DEG = DATA / "wl_random_subtracted_icl_split_2026-11-05_03/summary.json"
DEFAULT_OUTPUT = DATA / "wl_random_subtracted_icl_split_blockwidth_2026-11-05_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> dict:
    doc = json.loads(path.read_text(encoding="utf-8"))
    if (doc.get("schema") != "astrafek-hsc-random-subtracted-icl-split-v1"
            or doc.get("status") != "complete_exploratory"
            or doc.get("science_measurement") is not False
            or doc.get("filter") != "i2"
            or doc.get("source_tomographic_bins") != [2, 3, 4]):
        raise ValueError(f"unexpected ICL-split scope/status: {path}")
    if len(doc.get("matched_lens_pairs", [])) < 10 or len(doc.get("profile", [])) != 7:
        raise ValueError(f"incomplete ICL-split sample/profile: {path}")
    if doc.get("analysis_script_sha256") != "705d6389c70f1a7d6dce6ff2c07bc1517ede631b0a507757b81ef39417bb4a3f":
        raise ValueError(f"ICL-split runner provenance mismatch: {path}")
    return doc


def run(path_5deg: Path, path_2p5deg: Path, output: Path) -> dict:
    path_5deg, path_2p5deg, output = map(Path, (path_5deg, path_2p5deg, output))
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite block-width comparison: {output}")
    wide, narrow = load(path_5deg), load(path_2p5deg)
    if (wide.get("block_width_deg") != 5.0 or narrow.get("block_width_deg") != 2.5
            or wide.get("n_shared_delete_one_blocks") != 26
            or narrow.get("n_shared_delete_one_blocks") != 54):
        raise ValueError("expected the exact 5-degree/2.5-degree comparison scope")
    if (wide["matched_lens_pairs"] != narrow["matched_lens_pairs"]
            or wide["n_random_centers_in_selected_pair_sample"]
            != narrow["n_random_centers_in_selected_pair_sample"]
            or wide["icl_split_and_matching"] != narrow["icl_split_and_matching"]):
        raise ValueError("matched sample, ICL split, or random-center scope changed")
    point_5 = np.asarray([
        row["high_minus_low_delta_sigma_msun_mpc2"] for row in wide["profile"]
    ], dtype=float)
    point_2p5 = np.asarray([
        row["high_minus_low_delta_sigma_msun_mpc2"] for row in narrow["profile"]
    ], dtype=float)
    error_5 = np.asarray([
        row["high_minus_low_delta_sigma_error_msun_mpc2"] for row in wide["profile"]
    ], dtype=float)
    error_2p5 = np.asarray([
        row["high_minus_low_delta_sigma_error_msun_mpc2"] for row in narrow["profile"]
    ], dtype=float)
    if not np.array_equal(point_5, point_2p5):
        raise ValueError("tangential point estimates differ between block widths")
    ratios = error_2p5 / error_5
    bins = []
    for index, (left, right) in enumerate(zip(wide["profile"], narrow["profile"])):
        bins.append({
            "r_min_mpc": left["r_min_mpc"],
            "r_max_mpc": left["r_max_mpc"],
            "high_minus_low_delta_sigma_msun_mpc2": float(point_5[index]),
            "error_5deg_msun_mpc2": float(error_5[index]),
            "error_2p5deg_msun_mpc2": float(error_2p5[index]),
            "error_ratio_2p5_over_5": float(ratios[index]),
        })
    result = {
        "schema": "astrafek-wl-icl-split-blockwidth-comparison-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "source_rows_read": False,
        "new_data_requested": False,
        "method": "compare independently saved Astrafek ICL-split summaries; assert exact matched pairs, split definition, source bins, and identical point estimates before comparing block covariance",
        "scope": {
            "filter": "i2",
            "source_tomographic_bins": [2, 3, 4],
            "n_matched_lens_pairs": len(wide["matched_lens_pairs"]),
            "n_random_centers": wide["n_random_centers_in_selected_pair_sample"],
            "block_widths_deg": [5.0, 2.5],
            "occupied_blocks": [26, 54],
            "same_matched_pairs": True,
            "same_icl_split_and_calipers": True,
            "point_estimates_bitwise_equal": True,
        },
        "profile": bins,
        "tangential_null": {
            "5deg": wide["conditional_tangential_wald"],
            "2p5deg": narrow["conditional_tangential_wald"],
        },
        "cross_null": {
            "5deg": wide["conditional_cross_null"],
            "2p5deg": narrow["conditional_cross_null"],
        },
        "error_ratio_range_2p5_over_5": [float(ratios.min()), float(ratios.max())],
        "interpretation": "nominal tangential contrast changes from p=2.42e-4 at 5 degrees to p=0.099 at 2.5 degrees; it does not pass this block-width sensitivity and neither p value is calibrated because the regions are approximate equal-area cells",
        "inputs": {
            "5deg_summary": str(path_5deg.resolve()),
            "5deg_summary_sha256": sha256(path_5deg),
            "2p5deg_summary": str(path_2p5deg.resolve()),
            "2p5deg_summary_sha256": sha256(path_2p5deg),
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    result["comparison_summary"] = str(summary_path)
    result["comparison_summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "point_estimates_bitwise_equal": True,
        "error_ratio_range": result["error_ratio_range_2p5_over_5"],
        "p_tangential_5deg": wide["conditional_tangential_wald"]["nominal_p"],
        "p_tangential_2p5deg": narrow["conditional_tangential_wald"]["nominal_p"],
        "comparison_summary_sha256": result["comparison_summary_sha256"],
    }, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--five-degree", type=Path, default=DEFAULT_5DEG)
    parser.add_argument("--two-point-five-degree", type=Path, default=DEFAULT_2P5DEG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.five_degree, args.two_point_five_degree, args.output)


if __name__ == "__main__":
    main()
