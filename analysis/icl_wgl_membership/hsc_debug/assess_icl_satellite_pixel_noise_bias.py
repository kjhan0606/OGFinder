#!/usr/bin/env python3
"""Estimate pixel-noise impact on the fixed ICL--satellite association sample.

Uses uncorrected surface-brightness measurements so per-pair variances and
covariances stay in the same native measurement system. It subtracts the
known finite-lag pixel-noise contribution from the residual cross-products
after controlling for redshift and richness. Spatial-block bootstrap intervals
remain exploratory and do not include large-scale sky or calibration errors.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def corrected_partial_correlation(x, y, z, richness, variance_x, variance_y,
                                  covariance_xy, weights=None):
    """Return observed and noise-corrected weighted partial correlations."""
    x, y, z, richness, variance_x, variance_y, covariance_xy = [
        np.asarray(values, dtype=float) for values in
        (x, y, z, richness, variance_x, variance_y, covariance_xy)
    ]
    n = len(x)
    if any(values.shape != (n,) for values in
           (y, z, richness, variance_x, variance_y, covariance_xy)):
        raise ValueError("all input vectors must have matching one-dimensional shape")
    w = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != (n,) or not np.all(np.isfinite(w)) or np.any(w < 0) or np.count_nonzero(w) < 4:
        return None
    active = w > 0
    if np.linalg.matrix_rank(np.column_stack((np.ones(np.count_nonzero(active)), z[active], richness[active]))) < 3:
        return None
    design = np.column_stack((np.ones(n), z, richness))
    gram = design.T @ (w[:, None] * design)
    if np.linalg.matrix_rank(gram) < 3:
        return None
    inverse = np.linalg.inv(gram)
    # Symmetric weighted residual-maker for the weighted residual cross-products.
    projector = np.diag(w) - (w[:, None] * design) @ inverse @ (design.T * w[None, :])
    xx = float(x @ projector @ x)
    yy = float(y @ projector @ y)
    xy = float(x @ projector @ y)
    raw = xy / np.sqrt(xx * yy) if xx > 0 and yy > 0 else None
    diagonal = np.diag(projector)
    corrected_xx = xx - float(np.sum(diagonal * variance_x))
    corrected_yy = yy - float(np.sum(diagonal * variance_y))
    corrected_xy = xy - float(np.sum(diagonal * covariance_xy))
    corrected = (corrected_xy / np.sqrt(corrected_xx * corrected_yy)
                 if corrected_xx > 0 and corrected_yy > 0 else None)
    return {
        "raw_partial_r": raw,
        "pixel_noise_corrected_partial_r": corrected,
        "raw_residual_ss_x": xx,
        "raw_residual_ss_y": yy,
        "pixel_noise_residual_ss_x": corrected_xx,
        "pixel_noise_residual_ss_y": corrected_yy,
        "raw_residual_cross_product": xy,
        "pixel_noise_corrected_residual_cross_product": corrected_xy,
        "pixel_noise_fraction_residual_ss_x": float(np.sum(diagonal * variance_x) / xx) if xx > 0 else None,
        "pixel_noise_fraction_residual_ss_y": float(np.sum(diagonal * variance_y) / yy) if yy > 0 else None,
        "pixel_noise_cross_product_fraction": float(np.sum(diagonal * covariance_xy) / xy) if abs(xy) > 0 else None,
    }


def run(args):
    association = read_tsv(args.association_rows)
    block_rows = read_tsv(args.reblocked_association_rows)
    covariance = read_tsv(args.target_control_covariance / "target_control_pixel_covariance.tsv")
    scope = json.loads(args.scope_verification.read_text(encoding="utf-8"))
    cov_summary_path = args.target_control_covariance / "summary.json"
    cov_summary = json.loads(cov_summary_path.read_text(encoding="utf-8"))
    if scope.get("status") != "verified" or scope.get("science_measurement") is not False:
        raise ValueError("independent covariance support-scope verification did not pass")
    if scope.get("n_association_pairs_covered") != len(association) or scope.get("missing_association_pairs"):
        raise ValueError("the fixed association sample is not fully covered by verified covariance")
    if sha256(args.target_control_covariance / "target_control_pixel_covariance.tsv") != cov_summary.get("output_sha256", {}).get(
            str((args.target_control_covariance / "target_control_pixel_covariance.tsv").resolve())):
        raise ValueError("target/control covariance table checksum mismatch")
    block_by_key = {(str(row["target_id"]), str(row["filter"])): str(row["connected_5deg_block"])
                    for row in block_rows}
    covariance_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in covariance}
    if len(block_by_key) != len(block_rows) or len(covariance_by_key) != len(covariance):
        raise ValueError("duplicate pair keys in association or covariance input")

    rows_by_band = {band: [] for band in ("i", "i2")}
    for row in association:
        key = (str(row["cluster_id"]), str(row["filter"]))
        block = block_by_key.get(key)
        cov = covariance_by_key.get(key)
        if block is None or cov is None:
            raise ValueError(f"association pair lacks block label or covariance: {key}")
        values = {
            "cluster_id": key[0], "filter": key[1], "block": block,
            "x": float(row["icl_delta_raw_njy_arcsec2"]),
            "y": float(row["satellite_delta_raw_njy_arcsec2"]),
            "z": float(row["z"]), "richness": float(row["richness"]),
            "var_x": float(cov["icl_delta_lag_variance"]),
            "var_y": float(cov["target_satellite_lag_variance_njy2_arcsec4"]),
            "cov_xy": float(cov["target_icl_satellite_lag_covariance_njy2_arcsec4"]),
        }
        if not all(np.isfinite(value) for name, value in values.items()
                   if name not in {"cluster_id", "filter", "block"}):
            raise ValueError(f"non-finite association/noise input for {key}")
        rows_by_band[key[1]].append(values)

    output = {}
    for band_index, band in enumerate(("i", "i2")):
        rows = rows_by_band[band]
        if len(rows) < 8:
            raise ValueError(f"too few association pairs for {band}")
        arrays = {name: np.asarray([row[key] for row in rows], dtype=float)
                  for name, key in (("x", "x"), ("y", "y"), ("z", "z"), ("richness", "richness"),
                                    ("variance_x", "var_x"), ("variance_y", "var_y"),
                                    ("covariance_xy", "cov_xy"))}
        arrays["z"] = (arrays["z"] - arrays["z"].mean()) / arrays["z"].std(ddof=1)
        arrays["richness"] = (arrays["richness"] - arrays["richness"].mean()) / arrays["richness"].std(ddof=1)
        point = corrected_partial_correlation(**arrays)
        block_labels = np.asarray([row["block"] for row in rows], dtype=str)
        unique_blocks = np.unique(block_labels)
        index_by_block = {block: np.flatnonzero(block_labels == block) for block in unique_blocks}
        rng = np.random.default_rng(args.seed + band_index)
        boot_raw, boot_corrected = [], []
        for _ in range(args.replicas):
            draws = rng.choice(unique_blocks, size=len(unique_blocks), replace=True)
            weight = np.zeros(len(rows), dtype=float)
            for block in draws:
                weight[index_by_block[block]] += 1.0
            result = corrected_partial_correlation(**arrays, weights=weight)
            if result is None:
                continue
            if result["raw_partial_r"] is not None:
                boot_raw.append(result["raw_partial_r"])
            if result["pixel_noise_corrected_partial_r"] is not None:
                boot_corrected.append(result["pixel_noise_corrected_partial_r"])

        def interval(values):
            return np.quantile(np.asarray(values, dtype=float), [0.025, 0.5, 0.975]).tolist() if values else None

        output[band] = {
            "n_pairs": len(rows),
            "n_connected_5deg_blocks": int(len(unique_blocks)),
            "point_estimate": point,
            "raw_partial_r_bootstrap_2p5_50_97p5": interval(boot_raw),
            "noise_corrected_partial_r_bootstrap_2p5_50_97p5": interval(boot_corrected),
            "bootstrap_replicas_requested": args.replicas,
            "bootstrap_raw_valid": len(boot_raw),
            "bootstrap_noise_corrected_valid": len(boot_corrected),
        }

    return {
        "schema": "astrafek-icl-satellite-pixel-noise-bias-assessment-v1",
        "status": "exploratory",
        "science_measurement": False,
        "estimator": "partial Pearson correlation of raw ICL 200-300 kpc and satellite 50-75 kpc surface brightness, residualized against z and richness; subtract diagonal per-pair finite-lag pixel-noise covariance from residual cross-products",
        "measurement_units": "nJy arcsec^-2; covariance terms in (nJy arcsec^-2)^2",
        "bootstrap": {
            "method": "resample connected target/control 5-degree blocks within each filter; use block multiplicities as weights",
            "seed": args.seed,
            "replicas": args.replicas,
        },
        "filters": output,
        "input_sha256": {
            str(args.association_rows.resolve()): sha256(args.association_rows),
            str(args.reblocked_association_rows.resolve()): sha256(args.reblocked_association_rows),
            str((args.target_control_covariance / "target_control_pixel_covariance.tsv").resolve()): sha256(args.target_control_covariance / "target_control_pixel_covariance.tsv"),
            str(cov_summary_path.resolve()): sha256(cov_summary_path),
            str(args.scope_verification.resolve()): sha256(args.scope_verification),
            str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
        },
        "limitations": [
            "pixel-noise errors use the stationary finite 9x9 lag kernel and delivered variance plane; larger-scale sky covariance is excluded",
            "noise is assumed independent between different cluster/filter pairs; repeated survey-level modes are not represented",
            "the correction assumes the supplied covariance describes measurement noise in the raw, unextinguished estimators",
            "the spatial blocks are an approximate equal-area partition, not validated independent HSC fields",
            "redshift, richness, selection, source-mask, sky-transfer, extended-PSF, extinction, and calibration uncertainties remain",
            "the noise-corrected correlation is a diagnostic, not a causal association or science detection",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-rows", type=Path, required=True)
    parser.add_argument("--reblocked-association-rows", type=Path, required=True)
    parser.add_argument("--target-control-covariance", type=Path, required=True)
    parser.add_argument("--scope-verification", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 block resamples are required")
    if args.output.exists():
        parser.error(f"refusing to overwrite {args.output}")
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=False)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "filter_estimates": {band: values["point_estimate"] for band, values in result["filters"].items()},
        "output": str(args.output),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
