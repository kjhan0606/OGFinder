#!/usr/bin/env python3
"""Recompute ICL--satellite sky-residual sensitivity with sample-scoped blocks.

This reuses the completed 783-pair outer-residual table. For each filter, it
builds connected target/control sky blocks only from the fixed 74-pair
ICL--satellite association subset. It does not reread images or request data.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_icl_target_control_outer_sky import (
    ASSOCIATION_TABLE,
    DATA,
    block_bootstrap_association,
    pair_block_labels,
    read_json,
    read_tsv,
    sha256,
)


DEFAULT_INPUT = DATA / "icl_target_control_outer_sky_2026-11-03_01"
DEFAULT_OUTPUT = DATA / "icl_outer_sky_association_reblocked_2026-11-03_01"


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty association sample")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def replay(input_dir: Path, output_dir: Path, replicas: int = 5000,
           seed: int = 20261103) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    source_summary_path = input_dir / "summary.json"
    pair_path = input_dir / "paired_outer_sky_diagnostics.tsv"
    source_summary = read_json(source_summary_path)
    if (source_summary.get("status") != "complete_diagnostic"
            or source_summary.get("n_target_control_pairs") != 783
            or source_summary.get("science_measurement") is not False):
        raise ValueError("source outer-sky result is incomplete or has an unexpected scope")
    pair_rows = read_tsv(pair_path)
    if len(pair_rows) != 783:
        raise ValueError(f"expected 783 paired outer-sky rows, found {len(pair_rows)}")
    association_rows = read_tsv(ASSOCIATION_TABLE)
    association_by_key = {(str(row["cluster_id"]), str(row["filter"])): row
                          for row in association_rows}

    results, saved_rows, counts = {}, [], {}
    for band in ("i2", "i"):
        selected = []
        for pair in pair_rows:
            if pair["filter"] != band:
                continue
            assoc = association_by_key.get((pair["target_id"], band))
            if assoc is None:
                continue
            selected.append({
                "target_id": pair["target_id"], "control_id": pair["control_id"], "filter": band,
                "icl_delta_pixel_sfd_njy_arcsec2": float(assoc["icl_delta_pixel_sfd_njy_arcsec2"]),
                "satellite_delta_center_sfd_njy_arcsec2": float(assoc["satellite_delta_center_sfd_njy_arcsec2"]),
                "z": float(assoc["z"]), "richness": float(assoc["richness"]),
                "outer_plane_delta_njy_arcsec2": float(pair["target_minus_control_plane_constant_njy_arcsec2"]),
                "target_ra_deg": float(pair["target_ra_deg"]), "target_dec_deg": float(pair["target_dec_deg"]),
                "control_ra_deg": float(pair["control_ra_deg"]), "control_dec_deg": float(pair["control_dec_deg"]),
            })
        if len(selected) != (53 if band == "i2" else 21):
            raise ValueError(f"unexpected fixed association sample size for {band}: {len(selected)}")
        target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]}
                         for row in selected]
        control_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]}
                          for row in selected]
        _, labels, n_target_cells, n_control_cells, n_blocks, max_block = pair_block_labels(
            target_coords, control_coords, np.ones(len(selected), dtype=bool))
        for row, label in zip(selected, labels):
            row["connected_5deg_block"] = str(label)
            saved_rows.append(row)
        result = block_bootstrap_association(
            selected, np.asarray([row["connected_5deg_block"] for row in selected]), replicas,
            seed + (0 if band == "i2" else 1))
        block_meta = {"target_cells": n_target_cells, "control_cells": n_control_cells,
                      "connected_blocks": n_blocks, "largest_block_pairs": max_block,
                      "scope": "association pairs only; blocks are recomputed independently within each filter"}
        result["spatial_block_assignment"] = block_meta
        results[band] = result
        counts[band] = len(selected)

    output_dir.mkdir(parents=True, exist_ok=False)
    profile_path = output_dir / "association_rows.tsv"
    write_tsv(profile_path, saved_rows)
    summary = {
        "schema": "astrafek-icl-outer-sky-association-reblocked-v1",
        "status": "complete_exploratory_replay",
        "science_measurement": False,
        "method_correction": "recompute connected target/control 5-degree blocks from only the fixed ICL-satellite association sample within each filter; do not inherit labels from the full 783-pair ICL stack",
        "supersedes_source_summary_association_bootstrap": True,
        "reuses_outer_sky_frame_measurements_without_rereading_images": True,
        "n_association_pairs_by_filter": counts,
        "association_outer_residual_sensitivity_by_filter": results,
        "input_sha256": {
            str(source_summary_path.resolve()): sha256(source_summary_path),
            str(pair_path.resolve()): sha256(pair_path),
            str(ASSOCIATION_TABLE.resolve()): sha256(ASSOCIATION_TABLE),
            str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
        },
        "hsc_network_used": False,
        "data_requested": False,
        "sky_subtracted": False,
        "limitations": [
            "The 600--900 kpc target/control intercept difference is a residual indicator, not pure sky; it was not subtracted from measured fluxes.",
            "The 74-pair sample and ICL/satellite estimators remain exploratory and selected; blocks do not establish independent random fields.",
            "This replay corrects only the spatial-block assignment for the association bootstrap; it does not add PSF-wing, upstream-sky, completeness, or full image-pixel covariance corrections.",
        ],
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    hash_manifest = {
        "association_rows_sha256": sha256(profile_path),
        "summary_sha256": sha256(summary_path),
    }
    with (output_dir / "output_hashes.json").open("x", encoding="utf-8") as stream:
        json.dump(hash_manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": summary["status"], "pairs": counts,
                      "blocks": {band: results[band]["spatial_block_assignment"]["connected_blocks"]
                                 for band in results},
                      "results": results, "science_measurement": False}, ensure_ascii=False), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicas are required")
    try:
        replay(args.input.resolve(), args.output.resolve(), args.replicas, args.seed)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
