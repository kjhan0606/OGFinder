#!/usr/bin/env python3
"""Block-bootstrap the 17x17 to 33x33 ICL aperture-noise sensitivity."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_icl_supported_long_lag_noise import DEFAULT_SAMPLE, build_records
from evaluate_prestack_pair_support import pair_block_labels

ROOT = DEFAULT_SAMPLE
PREVIOUS = ROOT / "icl_aperture_covariance_200_300_2026-11-06_01"
EXTENDED = ROOT / "icl_aperture_covariance_200_300_16px_2026-11-06_01"
OUTPUT = ROOT / "icl_aperture_covariance_16px_block_bootstrap_2026-11-06_01"
REPLICAS = 5000
SEED = 20261106


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_pair_variances(directory: Path, width: int) -> tuple[dict, dict]:
    summary_path = directory / "summary.json"
    progress_path = directory / "progress.json"
    table_path = directory / "pair_aperture_covariance.tsv"
    summary, progress = read_json(summary_path), read_json(progress_path)
    if (summary.get("status") != "complete_exploratory_diagnostic"
            or summary.get("science_measurement") is not False
            or summary.get("n_supported_pairs") != 224
            or progress.get("status") != "complete_exploratory_diagnostic"
            or progress.get("output_hashes", {}).get("summary_sha256") != sha256(summary_path)
            or progress.get("output_hashes", {}).get("pair_aperture_covariance_sha256") != sha256(table_path)):
        raise ValueError(f"incomplete or changed {width}x{width} aperture covariance output")
    rows = read_tsv(table_path)
    if len(rows) != 224:
        raise ValueError(f"expected 224 {width}x{width} pair rows")
    key = f"target_finite_{width}x{width}_lag_variance"
    control = f"control_finite_{width}x{width}_lag_variance"
    result = {}
    for row in rows:
        pair = (str(row["cluster_id"]), str(row["filter"]))
        values = np.asarray([float(row[key]), float(row[control])], dtype=float)
        if pair in result or not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError(f"invalid/duplicate aperture variance for {pair}")
        result[pair] = float(values.sum())
    return result, {"summary": summary_path, "progress": progress_path, "table": table_path}


def bootstrap_band(band: str, ids: list[str], previous: dict, extended: dict,
                   pair_lookup: dict, replicas: int, seed: int) -> dict:
    target_rows = [{"ra_deg": pair_lookup[(ident, band)]["target_ra_deg"],
                    "dec_deg": pair_lookup[(ident, band)]["target_dec_deg"]}
                   for ident in ids]
    control_rows = [{"ra_deg": pair_lookup[(ident, band)]["control_ra_deg"],
                     "dec_deg": pair_lookup[(ident, band)]["control_dec_deg"]}
                    for ident in ids]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_rows, control_rows, np.ones(len(ids), dtype=bool),
    )
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    members = [np.flatnonzero(labels == group) for group in groups]
    old_var = np.asarray([previous[(ident, band)] for ident in ids], dtype=float)
    new_var = np.asarray([extended[(ident, band)] for ident in ids], dtype=float)
    old_sigma = np.sqrt(old_var.sum() / len(ids) ** 2)
    new_sigma = np.sqrt(new_var.sum() / len(ids) ** 2)
    rng = np.random.default_rng(seed)
    old_draw = np.empty(replicas)
    new_draw = np.empty(replicas)
    for index in range(replicas):
        selection = np.concatenate([members[int(j)] for j in rng.integers(0, len(members), size=len(members))])
        n = len(selection)
        old_draw[index] = np.sqrt(old_var[selection].sum() / n**2)
        new_draw[index] = np.sqrt(new_var[selection].sum() / n**2)
    difference = new_draw - old_draw
    ratio = new_draw / old_draw
    return {
        "n_pairs": len(ids), "n_target_cells": n_target_cells,
        "n_control_cells": n_control_cells, "n_connected_5deg_blocks": n_blocks,
        "largest_block_pairs": largest,
        "previous_17x17_stack_sigma_njy_arcsec2": float(old_sigma),
        "extended_33x33_stack_sigma_njy_arcsec2": float(new_sigma),
        "point_estimate_sigma_change_percent": float(100 * (new_sigma / old_sigma - 1)),
        "paired_block_bootstrap_sigma_change_percent_ci95": np.quantile(100 * (ratio - 1), [.025, .975]).tolist(),
        "paired_block_bootstrap_sigma_difference_njy_arcsec2_ci95": np.quantile(difference, [.025, .975]).tolist(),
        "bootstrap_replicas": replicas, "seed": seed,
        "interpretation": "conditional uncertainty for the stationarity/patch-median covariance estimators; not uncertainty of total ICL flux noise",
    }


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    old, old_paths = load_pair_variances(PREVIOUS, 17)
    new, new_paths = load_pair_variances(EXTENDED, 33)
    records, setup = build_records(DEFAULT_SAMPLE.resolve())
    pair_lookup = {(str(row["cluster_id"]), str(row["filter"])): row
                   for row in setup["pair_rows"]}
    filters = {}
    for band in ("i", "i2"):
        ids = sorted(ident for ident, filt in new if filt == band)
        if len(ids) != (71 if band == "i" else 153):
            raise ValueError(f"unexpected supported pair count for {band}")
        filters[band] = bootstrap_band(band, ids, old, new, pair_lookup,
                                       REPLICAS, SEED + (0 if band == "i" else 1))
    input_paths = [*old_paths.values(), *new_paths.values()]
    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    summary = {
        "schema": "astrafek-icl-aperture-covariance-16px-block-bootstrap-v1",
        "status": "complete_exploratory_diagnostic", "science_measurement": False,
        "n_supported_pairs": len(new), "filters": filters,
        "method": "paired 5-degree connected-block resampling of exact per-pair target+control annulus variances under the 17x17 and 33x33 stationary median-kernel estimates",
        "limitations": [
            "This bootstrap quantifies spatial-block sampling variation conditional on each estimated covariance kernel; it does not propagate uncertainty in lag estimates within patches.",
            "The 33x33 lag kernel is patch-median covariance after affine-plane subtraction; negative long-lag shells may reflect estimator filtering and are not treated as a physical correction.",
            "Source Poisson, sky transfer, PSF wings, masks, nonstationarity, longer lags, and cross-frame covariance remain unmodeled.",
            "No measured ICL profile or science flag is changed.",
        ],
        "input_sha256": input_hashes,
        "records_rebuilt": len(records),
        "hsc_network_used": False, "data_requested": False,
    }
    out = OUTPUT / "summary.json"
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "summary_sha256": sha256(out),
                      "filters": filters, "science_measurement": False}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
