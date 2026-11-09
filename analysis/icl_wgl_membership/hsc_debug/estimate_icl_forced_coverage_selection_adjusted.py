#!/usr/bin/env python3
"""Estimate ICL residual differences by 20/20 profile support within mixed sky blocks.

Fits a block-fixed-effect linear model with redshift and log richness on the
existing 224 ICL-supported pairs. The estimate is restricted to connected
5-degree blocks containing both selected and excluded pairs. A block bootstrap
describes uncertainty; this is an exploratory support audit, not a correction.
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
DEFAULT_OUTPUT = DATA / "icl_forced_coverage_selection_adjusted_2026-11-06_01"
REPLICATES = 10_000
SEED = 20261109
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


def fit_within_block(block_data: dict[str, dict[str, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    """Return coefficients and predictor matrix for block-demeaned OLS."""
    x_blocks, y_blocks = [], []
    for values in block_data.values():
        x_blocks.append(values["x"])
        y_blocks.append(values["y"])
    x = np.vstack(x_blocks)
    y = np.concatenate(y_blocks)
    gram = x.T @ x
    if np.linalg.matrix_rank(gram) < x.shape[1]:
        raise ValueError("within-block design matrix is rank deficient")
    return np.linalg.solve(gram, x.T @ y), x


def run(output: Path, replicates: int = REPLICATES) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicates < 1000:
        raise ValueError("at least 1,000 block-bootstrap replicates are required")
    for path in (SUMMARY, PAIRS, TARGETS, CONTROLS):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path, expected in EXPECTED.items():
        if sha256(path) != expected:
            raise ValueError(f"pinned input hash mismatch: {path}")
    source_summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if (source_summary.get("status") != "complete_exploratory"
            or source_summary.get("science_measurement") is not False):
        raise ValueError("pinned forced-coverage comparison has unexpected status")
    for path_text, expected in source_summary.get("inputs_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"forced-coverage source input hash mismatch: {path}")

    targets = json.loads(TARGETS.read_text(encoding="utf-8"))
    controls = json.loads(CONTROLS.read_text(encoding="utf-8"))
    target_by_id = {str(row["id"]): row for row in targets["clusters"]}
    control_by_id = {str(row["target_cluster_id"]): row for row in controls["clusters"]}
    all_rows = read_tsv(PAIRS)
    band_results: dict[str, dict] = {}
    block_table: list[dict] = []
    coefficient_draws_by_band: dict[str, np.ndarray] = {}

    for band_index, band in enumerate(("i", "i2")):
        rows = [row for row in all_rows if row["filter"] == band]
        tcoords, ccoords = [], []
        for row in rows:
            cid = row["cluster_id"]
            if cid not in target_by_id or cid not in control_by_id:
                raise ValueError(f"missing target/control coordinate for cluster {cid}")
            target, control = target_by_id[cid], control_by_id[cid]
            tcoords.append({"ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"])})
            ccoords.append({"ra_deg": float(control["ra_deg"]), "dec_deg": float(control["dec_deg"])})
        _, labels, *_ = pair_block_labels(tcoords, ccoords, np.ones(len(rows), dtype=bool))
        grouped: dict[str, list[dict]] = {}
        for row, label in zip(rows, labels):
            grouped.setdefault(str(label), []).append(row)
        mixed = {
            label: group for label, group in grouped.items()
            if any(row["selected_20_20"] == "True" for row in group)
            and any(row["selected_20_20"] != "True" for row in group)
        }
        if len(mixed) < 5:
            raise ValueError(f"too few mixed spatial blocks in {band}: {len(mixed)}")

        block_data: dict[str, dict[str, np.ndarray]] = {}
        for label, group in sorted(mixed.items()):
            y = np.asarray([float(row["icl_delta_200_300"]) for row in group])
            selected = np.asarray([row["selected_20_20"] == "True" for row in group], dtype=float)
            redshift = np.asarray([float(row["redshift"]) for row in group])
            log_richness = np.log([float(row["richness"]) for row in group])
            x = np.column_stack((selected, redshift, log_richness))
            # Within transformation is algebraically equivalent to including
            # one fixed intercept per connected spatial block.
            x -= x.mean(axis=0, keepdims=True)
            y -= y.mean()
            block_data[label] = {"x": x, "y": y}
            nsel = int(selected.sum())
            block_table.append({"filter": band, "block": label, "n_pairs": len(group),
                                "n_selected_20_20": nsel, "n_excluded": len(group) - nsel,
                                "mean_icl_delta": float(np.mean([float(row["icl_delta_200_300"]) for row in group]))})

        coef, x = fit_within_block(block_data)
        # The first regressor is selected_20_20; the remaining regressors are
        # redshift and log richness, all within-block centered.
        observed = float(coef[0])
        block_names = list(block_data)
        rng = np.random.default_rng(SEED + band_index)
        draws = np.empty(replicates, dtype=float)
        valid = 0
        for draw_index in range(replicates):
            chosen = rng.choice(block_names, size=len(block_names), replace=True)
            x_parts = [block_data[name]["x"] for name in chosen]
            y_parts = [block_data[name]["y"] for name in chosen]
            x_boot = np.vstack(x_parts)
            y_boot = np.concatenate(y_parts)
            gram = x_boot.T @ x_boot
            if np.linalg.matrix_rank(gram) < x_boot.shape[1]:
                continue
            draws[valid] = np.linalg.solve(gram, x_boot.T @ y_boot)[0]
            valid += 1
        draws = draws[:valid]
        if valid < 0.9 * replicates:
            raise ValueError(f"only {valid}/{replicates} valid block draws for {band}")
        ci = np.quantile(draws, [0.025, 0.975]).tolist()
        coefficient_draws_by_band[band] = draws
        selected_count = sum(row["selected_20_20"] == "True" for group in mixed.values() for row in group)
        excluded_count = sum(row["selected_20_20"] != "True" for group in mixed.values() for row in group)
        band_results[band] = {
            "estimand": "pair-weighted selected-minus-excluded ICL residual coefficient within mixed blocks, adjusted for redshift and log richness",
            "n_mixed_blocks": len(mixed), "n_pairs_in_mixed_blocks": sum(map(len, mixed.values())),
            "n_selected_in_mixed_blocks": selected_count, "n_excluded_in_mixed_blocks": excluded_count,
            "adjusted_selected_minus_excluded": observed,
            "block_bootstrap_ci95": ci,
            "block_bootstrap_replicates_requested": replicates,
            "block_bootstrap_replicates_valid": valid,
            "block_bootstrap_seed": SEED + band_index,
            "adjustment_terms": ["redshift", "log_richness", "connected_5deg_block_fixed_effects"],
        }

    output.mkdir(parents=True, exist_ok=False)
    block_path = output / "mixed_block_model_support.tsv"
    write_tsv(block_path, block_table,
              ["filter", "block", "n_pairs", "n_selected_20_20", "n_excluded", "mean_icl_delta"])
    inputs = (SUMMARY, PAIRS, TARGETS, CONTROLS)
    result = {
        "schema": "astrafek-icl-forced-coverage-adjusted-selection-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "image_pixels_read": False,
        "scope": "mixed connected 5-degree blocks among 224 ICL-supported pairs under the post-hoc 20/20 coverage rule",
        "method": "within-block OLS of ICL residual on support selection, redshift, and log richness; percentile bootstrap resampling of mixed spatial blocks",
        "by_filter": band_results,
        "inputs_sha256": {str(path): sha256(path) for path in inputs},
        "output_sha256": {block_path.name: sha256(block_path)},
        "limitations": [
            "coverage selection is deterministic; the adjusted coefficient is descriptive and not a causal effect",
            "analysis is restricted to mixed blocks and cannot transport to single-status blocks",
            "the post-hoc coverage floor and model form are exploratory; no multiplicity adjustment across thresholds is included",
            "measurement uncertainty in ICL residuals is not propagated",
            "forced CAMIRA positions without unique SEP detections are not confirmed isolated galaxies",
            "no correction is applied and science_measurement remains false",
        ],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    result = run(args.output, args.replicates)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "by_filter": result["by_filter"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
