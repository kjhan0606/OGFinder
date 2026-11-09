#!/usr/bin/env python3
"""Reassess ICL--satellite association on the saved 94-pair support sample.

This is an offline selection-scope diagnostic. It reads saved pair values and
coordinate manifests only; it does not read image pixels or submit HSC queries.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_TABLE = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03/matched_cluster_values.tsv"
ASSOCIATION_TABLE = ROOT / "icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv"
SAMPLE_MANIFEST = ROOT / "sample_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial(rows: list[dict]) -> float:
    return partial_pearson(
        [float(row["icl_delta_200_300"]) for row in rows],
        [float(row["satellite_delta_50_75"]) for row in rows],
        [[float(row["redshift"]) for row in rows],
         [float(row["richness"]) for row in rows]],
    )


def bootstrap(rows: list[dict], labels: np.ndarray, replicas: int,
              seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    if len(groups) < 6:
        raise ValueError(f"at least six spatial blocks are required; found {len(groups)}")
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        idx = np.concatenate([by_block[group] for group in chosen])
        value = partial([rows[index] for index in idx])
        if np.isfinite(value):
            draws.append(value)
    if len(draws) < replicas * 0.9:
        raise ValueError(f"too few valid block-bootstrap draws: {len(draws)}/{replicas}")
    return {
        "partial_r_z_richness": partial(rows),
        "block_bootstrap_ci95": np.percentile(draws, [2.5, 97.5]).tolist(),
        "bootstrap_fraction_positive": float(np.mean(np.asarray(draws) > 0)),
        "n_pairs": len(rows), "n_connected_5deg_blocks": int(len(groups)),
        "largest_block_pairs": int(max(len(index) for index in by_block.values())),
        "replicas_requested": replicas, "replicas_valid": len(draws),
        "seed": seed,
    }


def paired_scope_bootstrap(rows: list[dict], labels: np.ndarray,
                           is_narrow: np.ndarray, replicas: int,
                           seed: int) -> dict:
    """Compare nested 74/94 samples using the same 94-pair block resamples."""
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    narrow_draws, full_draws, delta_draws = [], [], []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        indices = [by_block[group] for group in chosen]
        all_idx = np.concatenate(indices)
        narrow_idx = all_idx[is_narrow[all_idx]]
        if len(narrow_idx) < 6:
            continue
        narrow_r = partial([rows[index] for index in narrow_idx])
        full_r = partial([rows[index] for index in all_idx])
        if np.isfinite(narrow_r) and np.isfinite(full_r):
            narrow_draws.append(narrow_r)
            full_draws.append(full_r)
            delta_draws.append(full_r - narrow_r)
    if len(delta_draws) < replicas * 0.8:
        raise ValueError(f"too few valid paired scope draws: {len(delta_draws)}/{replicas}")
    return {
        "scope_delta_full94_minus_narrow74": partial(rows) - partial(
            [row for row, keep in zip(rows, is_narrow) if keep]),
        "paired_delta_block_bootstrap_ci95": np.percentile(delta_draws, [2.5, 97.5]).tolist(),
        "narrow74_paired_bootstrap_ci95": np.percentile(narrow_draws, [2.5, 97.5]).tolist(),
        "full94_paired_bootstrap_ci95": np.percentile(full_draws, [2.5, 97.5]).tolist(),
        "replicas_requested": replicas, "replicas_valid": len(delta_draws), "seed": seed,
        "resampling_scope": "94-pair connected 5-degree blocks, applied identically to nested 74- and 94-pair samples",
    }


def run(output: Path, replicas: int, seed: int) -> dict:
    for path in (SUPPORT_TABLE, ASSOCIATION_TABLE, SAMPLE_MANIFEST, CONTROL_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")

    support = read_tsv(SUPPORT_TABLE)
    association = read_tsv(ASSOCIATION_TABLE)
    support_by_key = {(row["cluster_id"], row["filter"]): row for row in support}
    if len(support) != 94 or len(support_by_key) != 94:
        raise ValueError(f"expected 94 unique support pairs; got {len(support_by_key)}")
    association_by_key = {(row["cluster_id"], row["filter"]): row for row in association}
    if len(association_by_key) != len(association) or len(association) != 74:
        raise ValueError("the comparison association table must contain 74 unique pairs")
    if not set(association_by_key).issubset(support_by_key):
        raise ValueError("the 74-pair association sample is not nested in the 94-pair support sample")

    overlap_diffs = {"icl_raw_max_abs": 0.0, "satellite_raw_max_abs": 0.0}
    for key, old in association_by_key.items():
        new = support_by_key[key]
        overlap_diffs["icl_raw_max_abs"] = max(
            overlap_diffs["icl_raw_max_abs"],
            abs(float(old["icl_delta_raw_njy_arcsec2"]) - float(new["icl_delta_200_300"])),
        )
        overlap_diffs["satellite_raw_max_abs"] = max(
            overlap_diffs["satellite_raw_max_abs"],
            abs(float(old["satellite_delta_raw_njy_arcsec2"]) - float(new["satellite_delta_50_75"])),
        )
    if overlap_diffs["icl_raw_max_abs"] > 5e-7 or overlap_diffs["satellite_raw_max_abs"] > 1e-12:
        raise ValueError(f"saved pair values disagree on the 74-pair overlap: {overlap_diffs}")

    sample = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in sample["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798:
        raise ValueError("unexpected target/control manifest sizes")

    pair_rows = []
    by_filter = {}
    blocks_meta = {}
    for band in ("i", "i2"):
        rows = []
        for source in support:
            if source["filter"] != band:
                continue
            ident = source["cluster_id"]
            if ident not in targets or ident not in controls:
                raise ValueError(f"missing target or matched control coordinates for {ident}")
            target, control = targets[ident], controls[ident]
            row = {
                "cluster_id": ident, "filter": band,
                "control_id": str(controls[ident]["id"]),
                "icl_delta_200_300": float(source["icl_delta_200_300"]),
                "satellite_delta_50_75": float(source["satellite_delta_50_75"]),
                "redshift": float(source["redshift"]),
                "richness": float(source["richness"]),
                "target_ra_deg": float(target["ra_deg"]),
                "target_dec_deg": float(target["dec_deg"]),
                "control_ra_deg": float(control["ra_deg"]),
                "control_dec_deg": float(control["dec_deg"]),
                "in_previous_74_pair_sample": (ident, band) in association_by_key,
                "science_measurement": False,
            }
            rows.append(row)
        if len(rows) != (26 if band == "i" else 68):
            raise ValueError(f"unexpected {band}-band support count: {len(rows)}")
        target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in rows]
        control_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in rows]
        _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
            target_coords, control_coords, np.ones(len(rows), dtype=bool))
        for row, label in zip(rows, labels):
            row["connected_5deg_block"] = str(label)
        keep = np.asarray([row["in_previous_74_pair_sample"] for row in rows], dtype=bool)
        narrow_rows = [row for row, use in zip(rows, keep) if use]
        narrow_target = [coord for coord, use in zip(target_coords, keep) if use]
        narrow_control = [coord for coord, use in zip(control_coords, keep) if use]
        _, narrow_labels, nt_narrow, nc_narrow, nblocks_narrow, largest_narrow = pair_block_labels(
            narrow_target, narrow_control, np.ones(len(narrow_rows), dtype=bool))
        by_filter[band] = {
            "support94": bootstrap(rows, labels, replicas, seed + (0 if band == "i" else 10)),
            "previous74": bootstrap(narrow_rows, narrow_labels, replicas,
                                    seed + (1 if band == "i" else 11)),
            "paired_scope_change": paired_scope_bootstrap(
                rows, labels, keep, replicas, seed + (2 if band == "i" else 12)),
        }
        blocks_meta[band] = {
            "support94": {"target_cells": n_target_cells, "control_cells": n_control_cells,
                          "connected_blocks": n_blocks, "largest_block_pairs": largest},
            "previous74": {"target_cells": nt_narrow, "control_cells": nc_narrow,
                           "connected_blocks": nblocks_narrow, "largest_block_pairs": largest_narrow},
            "rule": "connected pairs sharing either endpoint's equal-area 5-degree sky cell; components recomputed within each sample/filter",
        }
        rows.sort(key=lambda row: int(row["cluster_id"]))
        pair_rows.extend(rows)

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "supported_pairs.tsv"
    with pair_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pair_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pair_rows)
    summary = {
        "schema": "astrafek-icl-satellite-supported-94-association-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "n_pairs": len(pair_rows),
        "filter_counts": {band: sum(row["filter"] == band for row in pair_rows) for band in ("i", "i2")},
        "n_previous_74_pairs": len(association_by_key),
        "n_added_supported_pairs": len(set(support_by_key) - set(association_by_key)),
        "overlap_pair_value_max_abs_difference": overlap_diffs,
        "estimator": "unweighted partial Pearson correlation of raw target-minus-control ICL 200-300 kpc surface brightness and raw satellite 50-75 kpc mean contrast, controlling linearly for redshift and richness",
        "support_selection": "the fixed 94-pair ICL plus per-satellite-profile support table; this is a selected support subset, not representative of all 800 parents",
        "block_bootstrap": "5-degree connected target/control sky-cell blocks, recomputed by filter and sample; percentile 95% intervals",
        "replicas": replicas,
        "seed_base": seed,
        "spatial_blocks_by_filter": blocks_meta,
        "results_by_filter": by_filter,
        "input_sha256": {str(path): sha256(path) for path in
                         (SUPPORT_TABLE, ASSOCIATION_TABLE, SAMPLE_MANIFEST, CONTROL_MANIFEST)},
        "output_sha256": {str(pair_path): sha256(pair_path)},
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "limitations": [
            "This measures sensitivity to expanding from the earlier 74-pair diagnostic input to the pre-existing 94-pair ICL-plus-satellite support intersection; it does not remove support-selection bias.",
            "Raw surface-brightness contrasts are not corrected here for extinction, pixel covariance, PSF wings, sky modes, completeness, or image-level systematics.",
            "The block bootstrap intervals describe this finite selected sample under the chosen spatial partition and are not a calibrated detection significance.",
            "No inference is made about causality or the full 800-cluster parent sample.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicas are required")
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({
        "status": result["status"], "n_pairs": result["n_pairs"],
        "filter_counts": result["filter_counts"],
        "result": {band: {
            "94_r": result["results_by_filter"][band]["support94"]["partial_r_z_richness"],
            "74_r": result["results_by_filter"][band]["previous74"]["partial_r_z_richness"],
            "94_ci95": result["results_by_filter"][band]["support94"]["block_bootstrap_ci95"],
            "94_blocks": result["results_by_filter"][band]["support94"]["n_connected_5deg_blocks"],
            "added_pair_delta_ci95": result["results_by_filter"][band]["paired_scope_change"]["paired_delta_block_bootstrap_ci95"],
        } for band in ("i", "i2")},
        "output_sha256": result["output_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
