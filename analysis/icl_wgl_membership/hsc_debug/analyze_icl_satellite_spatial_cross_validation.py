#!/usr/bin/env python3
"""Spatial leave-one-block-out prediction check for the ICL/satellite relation.

Uses only the pinned 94-pair Astrafek exploratory table. Within each filter,
it compares out-of-block prediction of satellite 50--75 kpc excess from
redshift and richness alone against the same model plus 200--300 kpc ICL
contrast. Spatial blocks are recomputed jointly from target and matched-control
coordinates. This is a predictive robustness diagnostic, not a causal or
physical ICL measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = DATA / "icl_satellite_supported_94_association_2026-11-05_01"
DEFAULT_OUTPUT = DATA / "icl_satellite_spatial_cv_2026-11-06_01"
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def fit_predict(train: list[dict[str, str]], test: list[dict[str, str]],
                include_icl: bool) -> np.ndarray:
    z_train = np.asarray([float(row["redshift"]) for row in train])
    richness_train = np.log(np.asarray([float(row["richness"]) for row in train]))
    icl_train = np.asarray([float(row["icl_delta_200_300"]) for row in train])
    z_mean, z_sd = float(z_train.mean()), float(z_train.std(ddof=0))
    r_mean, r_sd = float(richness_train.mean()), float(richness_train.std(ddof=0))
    i_mean, i_sd = float(icl_train.mean()), float(icl_train.std(ddof=0))
    if min(z_sd, r_sd, i_sd) <= 0:
        raise ValueError("training fold has a zero-variance covariate")

    def design(rows: list[dict[str, str]]) -> np.ndarray:
        z = (np.asarray([float(row["redshift"]) for row in rows]) - z_mean) / z_sd
        richness = (np.log(np.asarray([float(row["richness"]) for row in rows])) - r_mean) / r_sd
        columns = [np.ones(len(rows)), z, richness]
        if include_icl:
            icl = (np.asarray([float(row["icl_delta_200_300"]) for row in rows]) - i_mean) / i_sd
            columns.append(icl)
        return np.column_stack(columns)

    x_train, x_test = design(train), design(test)
    y_train = np.asarray([float(row["satellite_delta_50_75"]) for row in train])
    beta, _, rank, _ = np.linalg.lstsq(x_train, y_train, rcond=None)
    if rank != x_train.shape[1]:
        raise ValueError("training model is rank deficient")
    return x_test @ beta


def run(args: argparse.Namespace) -> None:
    input_dir, output = args.input_dir.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output: {output}")
    table = input_dir / "supported_pairs.tsv"
    summary_path = input_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete_diagnostic"
            or summary.get("science_measurement") is not False
            or summary.get("n_pairs") != 94
            or summary.get("output_sha256", {}).get(str(table)) != sha256(table)):
        raise ValueError("pinned 94-pair input failed status/hash validation")
    rows = read_tsv(table)
    if len(rows) != 94 or len({(row["cluster_id"], row["filter"]) for row in rows}) != 94:
        raise ValueError("expected 94 unique cluster/filter rows")

    from evaluate_prestack_pair_support import pair_block_labels

    targets = [{"ra_deg": float(row["target_ra_deg"]),
                "dec_deg": float(row["target_dec_deg"])} for row in rows]
    controls = [{"ra_deg": float(row["control_ra_deg"]),
                "dec_deg": float(row["control_dec_deg"])} for row in rows]
    _, labels, *_ = pair_block_labels(targets, controls, np.ones(len(rows), dtype=bool))
    for row, label in zip(rows, labels):
        row["spatial_block"] = str(label)

    filter_summaries, fold_rows = {}, []
    for band in FILTERS:
        selected = [row for row in rows if row["filter"] == band]
        blocks = sorted({row["spatial_block"] for row in selected})
        if len(blocks) < args.minimum_blocks or len(selected) < 10:
            raise ValueError(f"insufficient sample/block count in {band}")
        oof_null = np.full(len(selected), np.nan)
        oof_icl = np.full(len(selected), np.nan)
        for block in blocks:
            test_idx = [i for i, row in enumerate(selected) if row["spatial_block"] == block]
            train = [row for row in selected if row["spatial_block"] != block]
            test = [selected[i] for i in test_idx]
            if len(train) <= 4:
                raise ValueError(f"too few training rows after holding out {band}/{block}")
            oof_null[test_idx] = fit_predict(train, test, include_icl=False)
            oof_icl[test_idx] = fit_predict(train, test, include_icl=True)
            fold_rows.append({
                "filter": band, "held_out_block": block,
                "n_train": len(train), "n_test": len(test),
                "mse_nuisance_only": float(np.mean((np.asarray(
                    [float(row["satellite_delta_50_75"]) for row in test])
                    - oof_null[test_idx]) ** 2)),
                "mse_nuisance_plus_icl": float(np.mean((np.asarray(
                    [float(row["satellite_delta_50_75"]) for row in test])
                    - oof_icl[test_idx]) ** 2)),
            })
        y = np.asarray([float(row["satellite_delta_50_75"]) for row in selected])
        if not (np.isfinite(oof_null).all() and np.isfinite(oof_icl).all()):
            raise ValueError(f"incomplete out-of-block predictions for {band}")
        pointwise_gain = (y - oof_null) ** 2 - (y - oof_icl) ** 2
        mse_null = float(np.mean((y - oof_null) ** 2))
        mse_icl = float(np.mean((y - oof_icl) ** 2))
        rng = np.random.default_rng(args.seed + (0 if band == "i" else 1))
        by_block = {block: np.asarray([i for i, row in enumerate(selected)
                                       if row["spatial_block"] == block], dtype=int)
                    for block in blocks}
        gains = np.empty(args.replicas)
        for b in range(args.replicas):
            sampled = rng.choice(blocks, size=len(blocks), replace=True)
            idx = np.concatenate([by_block[block] for block in sampled])
            gains[b] = float(np.mean(pointwise_gain[idx]))
        filter_summaries[band] = {
            "n_pairs": len(selected), "n_spatial_blocks": len(blocks),
            "nuisance_only_oof_mse": mse_null,
            "nuisance_plus_icl_oof_mse": mse_icl,
            "relative_oof_mse_reduction_from_icl": (mse_null - mse_icl) / mse_null,
            "mean_paired_squared_error_gain": float(np.mean(pointwise_gain)),
            "block_bootstrap_ci95_mean_squared_error_gain": np.quantile(
                gains, [0.025, 0.975]).tolist(),
            "bootstrap_replicates": args.replicas,
            "bootstrap_seed": args.seed + (0 if band == "i" else 1),
        }

    output.mkdir(parents=True, exist_ok=False)
    fold_path = output / "leave_one_spatial_block_out.tsv"
    write_tsv(fold_path, fold_rows)
    result = {
        "schema": "astrafek-icl-satellite-spatial-cv-v1",
        "status": "complete_exploratory_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "input": str(table), "input_sha256": sha256(table),
        "input_summary_sha256": sha256(summary_path),
        "runner_sha256": sha256(Path(__file__)),
        "method": "leave-one-joint-connected-5-degree-target/control-block-out prediction; training-fold standardization; compare redshift+log-richness OLS against same OLS plus ICL; block bootstrap of paired out-of-block squared-error gain",
        "filter_results": filter_summaries,
        "fold_table": str(fold_path), "fold_table_sha256": sha256(fold_path),
        "limitations": [
            "uses the selected 94-pair support intersection, not all 800 clusters",
            "connected equal-area 5-degree cells approximate HSC field dependence",
            "cross-validation assesses prediction, not causality or a physical ICL/satellite decomposition",
            "input photometry omits source Poisson, full coadd weights, sky-transfer, and extended-PSF uncertainty",
            "the diagnostic was selected after earlier exploratory analyses; intervals are not a global research-family correction",
        ],
    }
    result_path = output / "summary.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({"summary": str(result_path), "results": filter_summaries}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261109)
    parser.add_argument("--minimum-blocks", type=int, default=6)
    args = parser.parse_args()
    if args.replicas < 1000 or args.minimum_blocks < 6:
        raise SystemExit("replicas must be >=1000 and minimum blocks >=6")
    run(args)


if __name__ == "__main__":
    main()
