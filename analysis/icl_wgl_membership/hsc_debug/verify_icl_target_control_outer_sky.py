#!/usr/bin/env python3
"""Independently audit the full-sample paired outer-sky diagnostic outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
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
from validate_control_images import load_controls


def close(a, b, *, atol=2e-11, rtol=2e-11):
    return math.isclose(float(a), float(b), abs_tol=atol, rel_tol=rtol)


def audit(work: Path) -> dict:
    summary_path = work / "summary.json"
    frame_path = work / "frame_outer_sky_diagnostics.tsv"
    pair_path = work / "paired_outer_sky_diagnostics.tsv"
    for path in (summary_path, frame_path, pair_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    summary = read_json(summary_path)
    if (summary.get("schema") != "astrafek-icl-target-control-outer-sky-pair-audit-v1"
            or summary.get("status") != "complete_diagnostic"
            or summary.get("science_measurement") is not False
            or summary.get("sky_subtracted") is not False
            or summary.get("n_target_control_pairs") != 783):
        raise ValueError("summary schema, completion state, or science gate mismatch")
    hashes = summary.get("input_sha256", {})
    if not hashes:
        raise ValueError("summary has no input hashes")
    input_hashes_match = {}
    for raw_path, expected in hashes.items():
        path = Path(raw_path)
        if path.name == "audit_icl_target_control_outer_sky.py":
            # The initial 783-pair job was run from this pinned commit. The
            # analysis source was subsequently corrected to scope bootstrap
            # blocks to the 74-pair association sample; verify the historical
            # source from git rather than comparing it to the corrected file.
            commit = "07ef9b5d6cb8cfa08d4cf8f32f2569277100ed7d"
            repo = Path(__file__).resolve().parents[3]
            relative = "analysis/icl_wgl_membership/hsc_debug/audit_icl_target_control_outer_sky.py"
            proc = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:{relative}"],
                                  check=True, capture_output=True)
            actual = hashlib.sha256(proc.stdout).hexdigest()
            input_hashes_match[raw_path] = actual == expected
            if actual != expected:
                raise ValueError("historical analyzer code hash differs from the run provenance")
            continue
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        input_hashes_match[raw_path] = actual == expected
        if actual != expected:
            raise ValueError(f"input SHA-256 changed: {path}")

    frames = read_tsv(frame_path)
    pairs = read_tsv(pair_path)
    if len(frames) != 1566 or len(pairs) != 783:
        raise ValueError(f"expected 1,566 frame and 783 pair rows; found {len(frames)}/{len(pairs)}")
    frame_lookup = {}
    for row in frames:
        key = (row["target_id"], row["role"])
        if key in frame_lookup or row["role"] not in {"target", "matched_control"}:
            raise ValueError(f"duplicate or invalid frame role: {key}")
        if row["filter"] not in {"i", "i2"} or int(row["n_tiles"]) < 10:
            raise ValueError(f"invalid filter or insufficient tiles in {key}")
        frame_lookup[key] = row
    seen_targets, seen_controls = set(), set()
    target_coords, control_coords = [], []
    recomputed_filter_counts = {"i": 0, "i2": 0}
    for row in pairs:
        target, control = row["target_id"], row["control_id"]
        if target in seen_targets or control in seen_controls:
            raise ValueError("duplicate target or matched-control assignment")
        seen_targets.add(target); seen_controls.add(control)
        tr = frame_lookup.get((target, "target"))
        cr = frame_lookup.get((target, "matched_control"))
        if tr is None or cr is None or tr["control_id"] != control or cr["control_id"] != control:
            raise ValueError(f"pair is not backed by exactly one target/control frame: {target}/{control}")
        if tr["filter"] != cr["filter"] or tr["filter"] != row["filter"]:
            raise ValueError(f"pair filter mismatch: {target}/{control}")
        plane_delta = float(tr["plane_constant_njy_arcsec2"]) - float(cr["plane_constant_njy_arcsec2"])
        median_delta = float(tr["median_tile_mean_njy_arcsec2"]) - float(cr["median_tile_mean_njy_arcsec2"])
        if not close(plane_delta, row["target_minus_control_plane_constant_njy_arcsec2"]):
            raise ValueError(f"plane-constant difference does not reproduce: {target}/{control}")
        if not close(median_delta, row["target_minus_control_median_tile_mean_njy_arcsec2"]):
            raise ValueError(f"tile-median difference does not reproduce: {target}/{control}")
        if (not close(float(tr["plane_residual_mad_over_nominal"]), row["target_plane_residual_mad_over_nominal"])
                or not close(float(cr["plane_residual_mad_over_nominal"]), row["control_plane_residual_mad_over_nominal"])):
            raise ValueError(f"frame residual diagnostics do not map to pair row: {target}/{control}")
        recomputed_filter_counts[row["filter"]] += 1
        target_coords.append({"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])})
        control_coords.append({"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])})
    if recomputed_filter_counts != summary.get("filter_counts"):
        raise ValueError("filter counts differ from summary")

    sample_doc = read_json(DATA / "sample_manifest.json")
    target_by_id = {str(row["id"]): row for row in sample_doc["clusters"]}
    mask_by_id = {str(row["id"]): row for row in read_tsv(
        DATA / "prestack_cluster_masks_2026-10-16_01/records.tsv")}
    controls_doc, _, _, _ = load_controls(DATA)
    control_by_target = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    expected_pairs = set()
    for ident, target in target_by_id.items():
        control = control_by_target.get(ident)
        mask = mask_by_id.get(ident)
        if control is not None and mask is not None and control["coadd_header_filter"] == mask["primary_filter"]:
            expected_pairs.add((ident, str(control["id"]), str(mask["primary_filter"])))
    actual_pairs = {(row["target_id"], row["control_id"], row["filter"]) for row in pairs}
    if actual_pairs != expected_pairs:
        raise ValueError("paired table does not match the full same-filter manifest selection")
    for row in pairs:
        target = target_by_id[row["target_id"]]
        control = control_by_target[row["target_id"]]
        for saved, source in ((row["target_ra_deg"], target["ra_deg"]),
                              (row["target_dec_deg"], target["dec_deg"]),
                              (row["control_ra_deg"], control["ra_deg"]),
                              (row["control_dec_deg"], control["dec_deg"])):
            if not close(saved, source, atol=1e-9, rtol=0):
                raise ValueError(f"saved pair coordinate differs from acquisition manifest: {row['target_id']}")

    _, labels, n_target_cells, n_control_cells, n_blocks, max_block = pair_block_labels(
        target_coords, control_coords, np.ones(len(pairs), dtype=bool))
    for row, label in zip(pairs, labels):
        if row["connected_5deg_block"] != str(label):
            raise ValueError("saved connected-block label does not reproduce")
    block_meta = summary["target_control_spatial_blocks"]
    expected_meta = (n_target_cells, n_control_cells, n_blocks, max_block)
    actual_meta = (block_meta["target_cells"], block_meta["control_cells"],
                   block_meta["connected_blocks"], block_meta["largest_block_pairs"])
    if expected_meta != actual_meta:
        raise ValueError(f"spatial block metadata mismatch: {expected_meta} != {actual_meta}")

    assoc_rows = read_tsv(ASSOCIATION_TABLE)
    association_lookup = {(str(row["cluster_id"]), str(row["filter"])): row for row in assoc_rows}
    association_results = {}
    for band in ("i2", "i"):
        merged = []
        for row in pairs:
            if row["filter"] != band:
                continue
            association = association_lookup.get((row["target_id"], band))
            if association is None:
                continue
            merged.append({
                "target_id": row["target_id"], "filter": band,
                "icl_delta_pixel_sfd_njy_arcsec2": float(association["icl_delta_pixel_sfd_njy_arcsec2"]),
                "satellite_delta_center_sfd_njy_arcsec2": float(association["satellite_delta_center_sfd_njy_arcsec2"]),
                "z": float(association["z"]), "richness": float(association["richness"]),
                "outer_plane_delta_njy_arcsec2": float(row["target_minus_control_plane_constant_njy_arcsec2"]),
                "connected_5deg_block": row["connected_5deg_block"],
            })
        expected = summary["association_outer_residual_sensitivity_by_filter"].get(band)
        if expected is None:
            if merged:
                raise ValueError(f"summary omits existing {band}-band association rows")
            continue
        replay = block_bootstrap_association(
            merged, np.asarray([row["connected_5deg_block"] for row in merged]),
            int(expected["bootstrap_replicas"]), int(expected["seed"]))
        for key in ("n_pairs", "n_connected_blocks", "n_valid_base_replicas", "n_valid_sky_adjusted_replicas"):
            if replay[key] != expected[key]:
                raise ValueError(f"{band} association {key} does not reproduce")
        scalar_keys = ("baseline_partial_r_z_richness", "outer_residual_adjusted_partial_r_z_richness_sky",
                       "paired_adjusted_minus_baseline_r")
        for key in scalar_keys:
            if not close(replay[key], expected[key], atol=2e-10, rtol=2e-10):
                raise ValueError(f"{band} association result {key} does not reproduce")
        for key in ("baseline_block_bootstrap_ci95", "outer_residual_adjusted_block_bootstrap_ci95", "paired_delta_ci95"):
            if any(not close(a, b, atol=2e-10, rtol=2e-10) for a, b in zip(replay[key], expected[key])):
                raise ValueError(f"{band} association interval {key} does not reproduce")
        association_results[band] = {"n_pairs": len(merged), "replay_matches": True}

    audit = {
        "schema": "astrafek-icl-target-control-outer-sky-independent-output-audit-v1",
        "status": "verified", "science_measurement": False,
        "n_frame_rows": len(frames), "n_target_control_pairs": len(pairs),
        "filter_counts": recomputed_filter_counts,
        "same_filter_plane_and_tile_deltas_reproduced": True,
        "full_same_filter_pair_selection_reproduced": True,
        "connected_spatial_blocks_reproduced": n_blocks,
        "association_bootstrap_replay": association_results,
        "input_sha256_rechecked": input_hashes_match,
        "image_content_rehashed": False,
        "historical_analysis_code_commit": "07ef9b5d6cb8cfa08d4cf8f32f2569277100ed7d",
        "legacy_association_bootstrap_replayed_as_recorded": True,
        "pair_table_sha256": sha256(pair_path),
        "frame_table_sha256": sha256(frame_path),
        "summary_sha256": sha256(summary_path),
        "limitations": [
            "This validates saved-table arithmetic, membership, full-sample block labels, and the historical deterministic bootstrap as recorded; it does not re-open or independently remeasure every science image.",
            "The historical association bootstrap inherited connected blocks from all 783 ICL pairs and is superseded for association uncertainty by the separately reblocked 74-pair replay.",
            "Outer-annulus intercepts are diagnostic residual indicators, not pure sky values or a validated correction.",
            "The source science gate remains false.",
        ],
    }
    output = work / "independent_output_audit.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    with output.open("x", encoding="utf-8") as stream:
        json.dump(audit, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    audit["audit_sha256"] = sha256(output)
    print(json.dumps({k: audit[k] for k in ("status", "n_frame_rows", "n_target_control_pairs",
                                            "filter_counts", "connected_spatial_blocks_reproduced",
                                            "association_bootstrap_replay", "audit_sha256")}, ensure_ascii=False))
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    audit(args.work.resolve())


if __name__ == "__main__":
    main()
