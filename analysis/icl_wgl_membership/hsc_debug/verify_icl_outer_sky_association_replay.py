#!/usr/bin/env python3
"""Verify the corrected, association-sample-scoped spatial bootstrap replay."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from audit_icl_target_control_outer_sky import (
    ASSOCIATION_TABLE,
    block_bootstrap_association,
    pair_block_labels,
    read_json,
    read_tsv,
    sha256,
)


def close(a, b, tol=2e-10):
    return math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol)


def audit(input_dir: Path, replay_dir: Path) -> dict:
    pair_path = input_dir / "paired_outer_sky_diagnostics.tsv"
    source_summary = input_dir / "summary.json"
    summary_path = replay_dir / "summary.json"
    rows_path = replay_dir / "association_rows.tsv"
    hashes_path = replay_dir / "output_hashes.json"
    summary = read_json(summary_path)
    if (summary.get("schema") != "astrafek-icl-outer-sky-association-reblocked-v1"
            or summary.get("status") != "complete_exploratory_replay"
            or summary.get("science_measurement") is not False
            or summary.get("supersedes_source_summary_association_bootstrap") is not True
            or summary.get("sky_subtracted") is not False):
        raise ValueError("reblocked replay has an invalid schema, status, or gate")
    expected_hashes = summary.get("input_sha256", {})
    for path in (source_summary, pair_path, ASSOCIATION_TABLE):
        recorded = expected_hashes.get(str(path.resolve()))
        if not recorded or recorded != sha256(path):
            raise ValueError(f"replay source hash does not match: {path}")
    rows = read_tsv(rows_path)
    source_pairs = read_tsv(pair_path)
    associations = read_tsv(ASSOCIATION_TABLE)
    assoc_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in associations}
    source_by_key = {(row["target_id"], row["filter"]): row for row in source_pairs}
    if len(source_by_key) != 783 or len(rows) != 74:
        raise ValueError("unexpected source or replay row count")

    replay_results = {}
    for band, expected_count in (("i2", 53), ("i", 21)):
        selected = []
        for pair in source_pairs:
            if pair["filter"] != band:
                continue
            assoc = assoc_by_key.get((pair["target_id"], band))
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
        if len(selected) != expected_count:
            raise ValueError(f"wrong association subset size in {band}: {len(selected)}")
        t = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in selected]
        c = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in selected]
        _, labels, nt, nc, nb, largest = pair_block_labels(t, c, np.ones(len(selected), dtype=bool))
        for row, label in zip(selected, labels):
            row["connected_5deg_block"] = str(label)
        saved = [row for row in rows if row["filter"] == band]
        saved_lookup = {row["target_id"]: row for row in saved}
        if len(saved_lookup) != len(saved) or set(saved_lookup) != {r["target_id"] for r in selected}:
            raise ValueError(f"replay table membership mismatch for {band}")
        for row in selected:
            actual = saved_lookup[row["target_id"]]
            for key in ("control_id", "connected_5deg_block"):
                if actual[key] != row[key]:
                    raise ValueError(f"saved {key} mismatch for {band}/{row['target_id']}")
            for key in ("icl_delta_pixel_sfd_njy_arcsec2", "satellite_delta_center_sfd_njy_arcsec2",
                        "z", "richness", "outer_plane_delta_njy_arcsec2"):
                if not close(actual[key], row[key]):
                    raise ValueError(f"saved {key} mismatch for {band}/{row['target_id']}")
        expected = summary["association_outer_residual_sensitivity_by_filter"][band]
        replay = block_bootstrap_association(
            selected, np.asarray([r["connected_5deg_block"] for r in selected]),
            int(expected["bootstrap_replicas"]), int(expected["seed"]))
        meta = expected["spatial_block_assignment"]
        if (replay["n_connected_blocks"] != nb or meta["connected_blocks"] != nb
                or meta["target_cells"] != nt or meta["control_cells"] != nc
                or meta["largest_block_pairs"] != largest):
            raise ValueError(f"association sample blocks failed to reproduce for {band}")
        for key in ("baseline_partial_r_z_richness", "outer_residual_adjusted_partial_r_z_richness_sky",
                    "paired_adjusted_minus_baseline_r"):
            if not close(replay[key], expected[key]):
                raise ValueError(f"bootstrap replay mismatch for {band}/{key}")
        for key in ("baseline_block_bootstrap_ci95", "outer_residual_adjusted_block_bootstrap_ci95", "paired_delta_ci95"):
            if any(not close(a, b) for a, b in zip(replay[key], expected[key])):
                raise ValueError(f"bootstrap interval mismatch for {band}/{key}")
        replay_results[band] = {"pairs": len(selected), "blocks": nb, "bootstrap_reproduced": True}

    hash_manifest = read_json(hashes_path)
    if hash_manifest.get("association_rows_sha256") != sha256(rows_path):
        raise ValueError("association-row output hash mismatch")
    if hash_manifest.get("summary_sha256") != sha256(summary_path):
        raise ValueError("replay-summary output hash mismatch")
    audit_result = {
        "schema": "astrafek-icl-outer-sky-association-replay-verification-v1",
        "status": "verified", "science_measurement": False,
        "association_blocks_recomputed_within_filter_sample": replay_results,
        "source_outer_pair_rows": len(source_pairs), "association_rows": len(rows),
        "association_table_sha256": sha256(rows_path),
        "summary_sha256": sha256(summary_path),
        "limitation": "This verifies deterministic table-level resampling and sample membership; it does not independently remeasure image pixels or establish pure-sky transfer.",
    }
    output = replay_dir / "independent_audit.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    with output.open("x", encoding="utf-8") as stream:
        json.dump(audit_result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(audit_result, ensure_ascii=False))
    return audit_result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--replay", type=Path, required=True)
    args = parser.parse_args()
    audit(args.input.resolve(), args.replay.resolve())


if __name__ == "__main__":
    main()
