#!/usr/bin/env python3
"""Recompute 74/94-pair ICL--satellite association intervals by block width."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from analyze_icl_satellite_supported_94 import (
    ASSOCIATION_TABLE,
    CONTROL_MANIFEST,
    SAMPLE_MANIFEST,
    SUPPORT_TABLE,
    partial,
    read_tsv,
    sha256,
)
from prestack_core import spatial_labels


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PARENT = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
PAIRS = PARENT / "supported_pairs.tsv"
PARENT_SUMMARY = PARENT / "summary.json"
DEFAULT_OUTPUT = ROOT / "icl_satellite_supported_94_blockwidth_2026-11-05_01"
WIDTHS = (2.5, 5.0, 10.0)


def connected_blocks(rows: list[dict], width_deg: float) -> np.ndarray:
    """Join pairs sharing a target or control equal-area sky cell."""
    n = len(rows)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    owner: dict[tuple[str, str], int] = {}
    target_ra = np.asarray([float(row["target_ra_deg"]) for row in rows])
    target_dec = np.asarray([float(row["target_dec_deg"]) for row in rows])
    control_ra = np.asarray([float(row["control_ra_deg"]) for row in rows])
    control_dec = np.asarray([float(row["control_dec_deg"]) for row in rows])
    target_cell = spatial_labels(target_ra, target_dec, width_deg)
    control_cell = spatial_labels(control_ra, control_dec, width_deg)
    for i in range(n):
        for key in (("target", str(target_cell[i])), ("control", str(control_cell[i]))):
            if key in owner:
                union(i, owner[key])
            else:
                owner[key] = i
    roots = sorted({find(i) for i in range(n)})
    labels = {root: f"block_{j:03d}" for j, root in enumerate(roots)}
    return np.asarray([labels[find(i)] for i in range(n)], dtype=str)


def partial_bootstrap(rows: list[dict], labels: np.ndarray,
                      replicas: int, seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    if len(groups) < 6:
        raise ValueError(f"at least six connected blocks required; found {len(groups)}")
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        selected = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_block[group] for group in selected])
        estimate = partial([rows[int(i)] for i in indices])
        if np.isfinite(estimate):
            draws.append(float(estimate))
    if len(draws) < replicas * 0.9:
        raise ValueError(f"valid bootstrap fraction too low: {len(draws)}/{replicas}")
    return {
        "partial_r_z_richness": float(partial(rows)),
        "block_bootstrap_ci95": np.percentile(draws, [2.5, 97.5]).tolist(),
        "bootstrap_fraction_positive": float(np.mean(np.asarray(draws) > 0)),
        "n_pairs": len(rows),
        "n_connected_blocks": int(len(groups)),
        "largest_block_pairs": int(max(map(len, by_block.values()))),
        "replicas_requested": replicas,
        "replicas_valid": len(draws),
        "seed": int(seed),
    }


def paired_scope_bootstrap(rows: list[dict], labels: np.ndarray,
                           keep_narrow: np.ndarray, replicas: int,
                           seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    keep_narrow = np.asarray(keep_narrow, dtype=bool)
    groups = np.unique(labels)
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    narrow_draws, full_draws, delta_draws = [], [], []
    for _ in range(replicas):
        selected = rng.choice(groups, size=len(groups), replace=True)
        all_idx = np.concatenate([by_block[group] for group in selected])
        narrow_idx = all_idx[keep_narrow[all_idx]]
        if len(narrow_idx) < 6:
            continue
        narrow_r = partial([rows[int(i)] for i in narrow_idx])
        full_r = partial([rows[int(i)] for i in all_idx])
        if np.isfinite(narrow_r) and np.isfinite(full_r):
            narrow_draws.append(float(narrow_r))
            full_draws.append(float(full_r))
            delta_draws.append(float(full_r - narrow_r))
    if len(delta_draws) < replicas * 0.8:
        raise ValueError(f"valid paired scope fraction too low: {len(delta_draws)}/{replicas}")
    narrow_rows = [row for row, include in zip(rows, keep_narrow) if include]
    return {
        "scope_delta_94_minus_74": float(partial(rows) - partial(narrow_rows)),
        "paired_delta_block_bootstrap_ci95": np.percentile(delta_draws, [2.5, 97.5]).tolist(),
        "support94_ci95_on_same_resamples": np.percentile(full_draws, [2.5, 97.5]).tolist(),
        "support74_ci95_on_same_resamples": np.percentile(narrow_draws, [2.5, 97.5]).tolist(),
        "replicas_requested": replicas,
        "replicas_valid": len(delta_draws),
        "seed": int(seed),
    }


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite ICL--satellite block-width audit: {output}")
    parent_summary = json.loads(PARENT_SUMMARY.read_text(encoding="utf-8"))
    if (parent_summary.get("status") != "complete_diagnostic"
            or parent_summary.get("science_measurement") is not False
            or sha256(PAIRS) != "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"):
        raise ValueError("the verified 94-pair parent output changed or is not diagnostic")
    rows = read_tsv(PAIRS)
    if len(rows) != 94 or any(row.get("science_measurement") != "False" for row in rows):
        raise ValueError("expected 94 exploratory pairs with science_measurement=false")
    narrow_keys = {
        (row["cluster_id"], row["filter"])
        for row in read_tsv(ASSOCIATION_TABLE)
    }
    keys = {(row["cluster_id"], row["filter"]) for row in rows}
    if len(narrow_keys) != 74 or not narrow_keys < keys:
        raise ValueError("the 74-pair diagnostic sample is not a strict subset of the 94-pair support table")
    keep_narrow = np.asarray([
        (row["cluster_id"], row["filter"]) in narrow_keys for row in rows
    ], dtype=bool)
    results = {}
    for width_index, width in enumerate(WIDTHS):
        results[str(width)] = {}
        for band_index, band in enumerate(("i", "i2")):
            band_rows = [row for row in rows if row["filter"] == band]
            band_narrow = np.asarray([
                (row["cluster_id"], row["filter"]) in narrow_keys for row in band_rows
            ], dtype=bool)
            if len(band_rows) != (26 if band == "i" else 68):
                raise ValueError(f"unexpected {band} support count: {len(band_rows)}")
            labels = connected_blocks(band_rows, width)
            narrow_rows = [row for row, keep in zip(band_rows, band_narrow) if keep]
            narrow_labels = connected_blocks(narrow_rows, width)
            results[str(width)][band] = {
                "support94": partial_bootstrap(
                    band_rows, labels, replicas, seed + width_index * 100 + band_index * 10,
                ),
                "previous74": partial_bootstrap(
                    narrow_rows, narrow_labels, replicas,
                    seed + width_index * 100 + band_index * 10 + 1,
                ),
                "paired_scope_change": paired_scope_bootstrap(
                    band_rows, labels, band_narrow, replicas,
                    seed + width_index * 100 + band_index * 10 + 2,
                ),
                "block_definition": "connected components when either target or control endpoint shares an equal-area sky cell of the stated width",
            }
    summary = {
        "schema": "astrafek-icl-satellite-supported-94-blockwidth-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "n_support_pairs": len(rows),
        "n_previous_diagnostic_pairs": len(narrow_keys),
        "block_widths_deg": list(WIDTHS),
        "replicas": replicas,
        "seed_base": seed,
        "estimator": "raw target-minus-control 200-300 kpc ICL and 50-75 kpc satellite partial Pearson r, controlling redshift and richness",
        "results_by_width_and_filter": results,
        "input_sha256": {
            str(PAIRS): sha256(PAIRS),
            str(PARENT_SUMMARY): sha256(PARENT_SUMMARY),
            str(SUPPORT_TABLE): sha256(SUPPORT_TABLE),
            str(ASSOCIATION_TABLE): sha256(ASSOCIATION_TABLE),
            str(SAMPLE_MANIFEST): sha256(SAMPLE_MANIFEST),
            str(CONTROL_MANIFEST): sha256(CONTROL_MANIFEST),
        },
        "limitations": [
            "block-width intervals describe the fixed, selected 94-pair support subset and its nested 74-pair diagnostic subset; neither represents the 800-cluster parent",
            "the block widths are approximate equal-area cells, not validated contiguous HSC footprint regions",
            "the estimator uses raw saved profile contrasts and does not correct image-level noise, PSF wings, sky transfer, or completeness",
            "this is a sensitivity diagnostic, not a calibrated detection or causal estimate",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "results": {
            width: {band: {
                "support94_r": result["support94"]["partial_r_z_richness"],
                "support94_ci95": result["support94"]["block_bootstrap_ci95"],
                "support94_blocks": result["support94"]["n_connected_blocks"],
            } for band, result in by_filter.items()}
            for width, by_filter in results.items()
        },
        "summary_sha256": summary["summary_sha256"],
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicates are required")
    run(args.output, args.replicas, args.seed)


if __name__ == "__main__":
    main()
