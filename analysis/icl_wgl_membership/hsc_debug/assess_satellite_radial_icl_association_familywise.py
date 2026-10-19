#!/usr/bin/env python3
"""Family-wise max-T block-bootstrap intervals for seven radial ICL associations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = DATA / "satellite_outer_icl_baseline_assoc_2026-10-20_07"
RADIAL_BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def simultaneous_intervals(point: np.ndarray, draws: np.ndarray,
                           alpha: float = 0.05) -> dict:
    point, draws = np.asarray(point, float), np.asarray(draws, float)
    if point.ndim != 1 or draws.ndim != 2 or draws.shape[1] != len(point):
        raise ValueError("point/draw arrays have incompatible dimensions")
    valid = np.isfinite(draws).all(axis=1)
    sample = draws[valid]
    if len(sample) < 0.8 * len(draws):
        raise ValueError(f"too few valid joint replicates: {len(sample)}/{len(draws)}")
    scale = np.std(sample, axis=0, ddof=1)
    if np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("max-T scale is non-finite or zero")
    max_t = np.max(np.abs((sample - point) / scale), axis=1)
    critical = float(np.quantile(max_t, 1.0 - alpha))
    return {
        "familywise_confidence": 1.0 - alpha,
        "max_t_critical_value": critical,
        "standard_errors": scale.tolist(),
        "lower": (point - critical * scale).tolist(),
        "upper": (point + critical * scale).tolist(),
        "valid_joint_replicates": int(len(sample)),
        "replicates_requested": int(len(draws)),
    }


def radial_matrices(rows: list[dict[str, str]], band: str):
    by_bin = {radial: {} for radial in RADIAL_BINS}
    for row in rows:
        if row["filter"] != band:
            continue
        radial = row["radial_bin_kpc"]
        cluster_id = str(row["cluster_id"])
        if radial not in by_bin or cluster_id in by_bin[radial]:
            raise ValueError(f"unexpected or duplicate row: {band}/{radial}/{cluster_id}")
        by_bin[radial][cluster_id] = row
    cluster_ids = sorted(by_bin[RADIAL_BINS[0]])
    if not cluster_ids or any(set(by_bin[r]) != set(cluster_ids) for r in RADIAL_BINS):
        raise ValueError(f"cluster sample is not identical in all radial bins for {band}")
    matrices = {name: np.empty((len(cluster_ids), len(RADIAL_BINS)), dtype=float)
                for name in ("icl_raw", "icl_outer_subtracted", "satellite")}
    targets, controls = [], []
    for i, cluster_id in enumerate(cluster_ids):
        reference = by_bin[RADIAL_BINS[0]][cluster_id]
        targets.append({"ra_deg": float(reference["ra_deg"]), "dec_deg": float(reference["dec_deg"])})
        controls.append({"ra_deg": float(reference["control_ra_deg"]),
                         "dec_deg": float(reference["control_dec_deg"])})
        for j, radial in enumerate(RADIAL_BINS):
            row = by_bin[radial][cluster_id]
            if row["control_id"] != reference["control_id"]:
                raise ValueError(f"control mapping varies across radial bins for {cluster_id}/{band}")
            matrices["icl_raw"][i, j] = float(row["icl_raw"])
            matrices["icl_outer_subtracted"][i, j] = float(row["icl_outer_subtracted"])
            matrices["satellite"][i, j] = float(row["satellite"])
    return cluster_ids, targets, controls, matrices


def estimate(rows: list[dict[str, str]], band: str, replicas: int, seed: int) -> dict:
    cluster_ids, targets, controls, values = radial_matrices(rows, band)
    n = len(cluster_ids)
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        targets, controls, np.ones(n, dtype=bool))
    if n_blocks < 3:
        raise ValueError(f"insufficient connected blocks for {band}: {n_blocks}")
    members = [np.flatnonzero(labels == label) for label in np.unique(labels)]
    point = {
        kind: np.asarray([spearmanr(values[kind][:, j], values["satellite"][:, j]).statistic
                          for j in range(len(RADIAL_BINS))], dtype=float)
        for kind in ("icl_raw", "icl_outer_subtracted")
    }
    point["corrected_minus_raw"] = point["icl_outer_subtracted"] - point["icl_raw"]
    rng = np.random.default_rng(seed)
    draws = {kind: np.full((replicas, len(RADIAL_BINS)), np.nan)
             for kind in point}
    for b in range(replicas):
        chosen = rng.integers(0, len(members), size=len(members))
        indices = np.concatenate([members[int(index)] for index in chosen])
        sampled_satellite = values["satellite"][indices]
        sampled_raw = values["icl_raw"][indices]
        sampled_corrected = values["icl_outer_subtracted"][indices]
        for j in range(len(RADIAL_BINS)):
            raw = float(spearmanr(sampled_raw[:, j], sampled_satellite[:, j]).statistic)
            corrected = float(spearmanr(sampled_corrected[:, j], sampled_satellite[:, j]).statistic)
            if math.isfinite(raw) and math.isfinite(corrected):
                draws["icl_raw"][b, j] = raw
                draws["icl_outer_subtracted"][b, j] = corrected
                draws["corrected_minus_raw"][b, j] = corrected - raw
    result = {
        "filter": band, "n_clusters": n, "n_target_cells": int(n_target),
        "n_control_cells": int(n_control), "n_connected_blocks": int(n_blocks),
        "largest_block_clusters": int(largest), "radial_bins_kpc": list(RADIAL_BINS),
        "point_estimates": {k: v.tolist() for k, v in point.items()},
        "familywise_max_t_95": {k: simultaneous_intervals(point[k], draws[k]) for k in point},
        "bootstrap": {"replicates": replicas, "seed": seed, "unit": "connected target/control 5-degree block"},
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 bootstrap replicates are required")
    if args.output_dir.exists():
        parser.error(f"refusing to overwrite existing output: {args.output_dir}")
    summary_path = args.input_dir / "summary.json"
    detail_path = args.input_dir / "matched_radial_cluster_values.tsv"
    source = json.loads(summary_path.read_text(encoding="utf-8"))
    if (source.get("status") != "complete_exploratory" or source.get("science_measurement") is not False
            or len(source.get("radial_common_support_results", [])) != 14):
        parser.error("unexpected or incomplete radial-association source")
    expected = source.get("outputs_sha256", {}).get(detail_path.name)
    if not expected or sha256(detail_path) != expected:
        parser.error("radial association input table checksum mismatch")
    rows = read_tsv(detail_path)
    analyses = [estimate(rows, band, args.replicas, args.seed + 10000 * i)
                for i, band in enumerate(FILTERS)]
    args.output_dir.mkdir(parents=True, exist_ok=False)
    result = {
        "schema": "astrafek-satellite-radial-icl-familywise-max-t-v1",
        "status": "complete_exploratory_familywise_sensitivity",
        "science_measurement": False,
        "method": {
            "estimand": "Spearman cluster-level association of satellite SB with raw or same-pair outer-subtracted ICL",
            "family": "seven radial bins per filter; raw, corrected, and corrected-minus-raw are separate families",
            "interval": "95% simultaneous max-T interval standardized by the connected-block bootstrap standard error",
            "sampling": "same connected 5-degree target/control blocks resampled jointly across all seven radial bins within each filter",
            "p_values": "none; exploratory intervals only",
        },
        "limitations": [
            "the 400-500 kpc annulus is not verified blank sky",
            "max-T intervals depend on the chosen connected-block partition and omit pixel/coadd covariance",
            "radial satellite member lists differ by bin despite fixed cluster support",
            "PSF wings, sky transfer, source completeness, extinction, and photometric systematics remain unresolved",
        ],
        "source_summary": {"path": str(summary_path.resolve()), "sha256": sha256(summary_path)},
        "source_radial_table": {"path": str(detail_path.resolve()), "sha256": sha256(detail_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "filters": analyses,
    }
    output = args.output_dir / "summary.json"
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "filters": [{"filter": item["filter"], "n_clusters": item["n_clusters"],
                     "raw_95_simultaneous": item["familywise_max_t_95"]["icl_raw"],
                     "outer_subtracted_95_simultaneous": item["familywise_max_t_95"]["icl_outer_subtracted"],
                     "delta_95_simultaneous": item["familywise_max_t_95"]["corrected_minus_raw"]}
                    for item in analyses],
        "output": str(args.output_dir),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
