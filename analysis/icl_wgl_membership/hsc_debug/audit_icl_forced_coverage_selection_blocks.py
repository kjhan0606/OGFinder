#!/usr/bin/env python3
"""Audit spatial positivity of the 20/20 forced-profile support selection.

Uses saved Astrafek tables and target/control manifests only. Reports which
connected 5-degree paired-sky blocks contain both supported and unsupported
ICL/satellite pairs, and descriptive within-block ICL differences. It does not
fit a missing-at-random correction or request HSC data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCE = DATA / "icl_satellite_forced_coverage_both_sides_2026-11-06_03"
SUMMARY = SOURCE / "summary.json"
PAIRS = SOURCE / "paired_cluster_comparison.tsv"
TARGETS = DATA / "download_01" / "download_manifest.json"
CONTROLS = DATA / "control_download_01" / "download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_forced_coverage_selection_blocks_2026-11-06_02"
BOOTSTRAP_REPLICATES = 10_000
SEED = 20261106
EXPECTED = {
    SUMMARY: "b8b79a65e1c838e757fc4ee6c1e1fdf61340d72e989bf035ab29bc0a7a6a314a",
    PAIRS: "b613b53457dc9a85aa1be2b3ed0de0026a47d86fa32897af392859c34f5fabd2",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    for path in (SUMMARY, PAIRS, TARGETS, CONTROLS):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path, expected in EXPECTED.items():
        if sha256(path) != expected:
            raise ValueError(f"pinned input hash mismatch: {path}")

    source_summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if (source_summary.get("schema") != "astrafek-icl-satellite-forced-coverage-paired-comparison-v1"
            or source_summary.get("status") != "complete_exploratory"
            or source_summary.get("science_measurement") is not False):
        raise ValueError("pinned 20/20 paired comparison has unexpected state")
    for path_text, expected in source_summary.get("inputs_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"source-comparison input hash mismatch: {path}")
    target_manifest = json.loads(TARGETS.read_text(encoding="utf-8"))
    control_manifest = json.loads(CONTROLS.read_text(encoding="utf-8"))
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    control_by_id = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    rows = read_tsv(PAIRS)
    if len(rows) != 224 or len({(row["cluster_id"], row["filter"]) for row in rows}) != 224:
        raise ValueError("expected one row for each of 224 ICL-supported cluster/filter pairs")

    block_rows: list[dict] = []
    contrasts: list[dict] = []
    band_summary: dict[str, dict] = {}
    for band in ("i", "i2"):
        selected_rows = [row for row in rows if row["filter"] == band]
        selected_mask = np.asarray([row["selected_20_20"] == "True" for row in selected_rows])
        target_coords, control_coords = [], []
        for row in selected_rows:
            cluster_id = row["cluster_id"]
            if cluster_id not in target_by_id or cluster_id not in control_by_id:
                raise ValueError(f"missing paired target/control coordinates for {cluster_id}")
            target = target_by_id[cluster_id]
            control = control_by_id[cluster_id]
            target_coords.append({"ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"])})
            control_coords.append({"ra_deg": float(control["ra_deg"]), "dec_deg": float(control["dec_deg"])})
        _, labels, n_target_cells, n_control_cells, _, _ = pair_block_labels(
            target_coords, control_coords, np.ones(len(selected_rows), dtype=bool))
        blocks: dict[str, list[int]] = {}
        for index, label in enumerate(labels):
            blocks.setdefault(str(label), []).append(index)
        selected_total = int(selected_mask.sum())
        excluded_total = int((~selected_mask).sum())
        mixed_block_ids = []
        selected_in_mixed = excluded_in_mixed = 0
        block_differences = []
        for block, indexes in sorted(blocks.items()):
            take = selected_mask[indexes]
            n_selected = int(take.sum())
            n_excluded = len(indexes) - n_selected
            mixed = n_selected > 0 and n_excluded > 0
            block_rows.append({
                "filter": band, "block": block, "n_pairs": len(indexes),
                "n_selected_20_20": n_selected, "n_excluded": n_excluded,
                "mixed_selection": mixed,
            })
            if not mixed:
                continue
            mixed_block_ids.append(block)
            selected_in_mixed += n_selected
            excluded_in_mixed += n_excluded
            y = np.asarray([float(selected_rows[index]["icl_delta_200_300"])
                            for index in indexes])
            selected_y = y[take]
            excluded_y = y[~take]
            difference = float(np.mean(selected_y) - np.mean(excluded_y))
            block_differences.append(difference)
            contrasts.append({
                "filter": band, "block": block,
                "n_selected_20_20": n_selected, "n_excluded": n_excluded,
                "selected_mean_icl_delta_200_300": float(np.mean(selected_y)),
                "excluded_mean_icl_delta_200_300": float(np.mean(excluded_y)),
                "selected_minus_excluded": difference,
            })

        rng = np.random.default_rng(SEED + (0 if band == "i" else 1))
        if block_differences:
            draws = rng.choice(np.asarray(block_differences),
                               size=(BOOTSTRAP_REPLICATES, len(block_differences)),
                               replace=True).mean(axis=1)
            block_mean_ci = np.quantile(draws, [0.025, 0.975]).tolist()
        else:
            block_mean_ci = None

        band_summary[band] = {
            "icl_supported_pairs": len(selected_rows),
            "selected_20_20_pairs": selected_total,
            "excluded_pairs": excluded_total,
            "target_5deg_cells": n_target_cells,
            "control_5deg_cells": n_control_cells,
            "connected_paired_sky_blocks": len(blocks),
            "mixed_selection_blocks": len(mixed_block_ids),
            "selected_pairs_in_mixed_blocks": selected_in_mixed,
            "excluded_pairs_in_mixed_blocks": excluded_in_mixed,
            "selected_pairs_in_single_status_blocks": selected_total - selected_in_mixed,
            "excluded_pairs_in_single_status_blocks": excluded_total - excluded_in_mixed,
            "mixed_block_selected_minus_excluded_icl_delta_mean": float(np.mean(block_differences)) if block_differences else None,
            "mixed_block_selected_minus_excluded_icl_delta_median": float(np.median(block_differences)) if block_differences else None,
            "mixed_block_difference_signs_span_zero": bool(min(block_differences) <= 0 <= max(block_differences)) if block_differences else None,
            "equal_block_mean_difference_bootstrap_ci95": block_mean_ci,
            "block_bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "block_bootstrap_seed": SEED + (0 if band == "i" else 1),
            "interpretation": "descriptive support-overlap diagnostic; block contrasts are not an inferential correction",
        }

    output.mkdir(parents=True, exist_ok=False)
    block_path = output / "spatial_block_support.tsv"
    contrast_path = output / "mixed_block_icl_contrasts.tsv"
    write_tsv(block_path, block_rows, ["filter", "block", "n_pairs", "n_selected_20_20", "n_excluded", "mixed_selection"])
    write_tsv(contrast_path, contrasts, ["filter", "block", "n_selected_20_20", "n_excluded", "selected_mean_icl_delta_200_300", "excluded_mean_icl_delta_200_300", "selected_minus_excluded"])
    inputs = (SUMMARY, PAIRS, TARGETS, CONTROLS)
    result = {
        "schema": "astrafek-icl-forced-coverage-spatial-selection-audit-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "image_pixels_read": False,
        "scope": "spatial support overlap for 224 existing ICL-supported cluster/filter pairs under the 20/20 target/control coverage rule",
        "method": "connected paired target/control sky blocks on a 5-degree grid; summarize selected/excluded pair counts and descriptive within-mixed-block ICL residual contrasts",
        "by_filter": band_summary,
        "inputs_sha256": {str(path): sha256(path) for path in inputs},
        "output_sha256": {block_path.name: sha256(block_path), contrast_path.name: sha256(contrast_path)},
        "limitations": [
            "selected/excluded support is deterministic and was not randomized; mixed-block contrasts are descriptive",
            "block construction approximates shared target/control spatial dependence and does not validate a population transport model",
            "no adjustment for redshift/richness is applied within the mixed-block contrasts",
            "the 20% coverage threshold is post-hoc; no multiplicity correction across thresholds is included",
            "forced CAMIRA positions without unique SEP detections are not confirmed isolated galaxies",
            "no correction is applied and science_measurement remains false",
        ],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "by_filter": result["by_filter"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
