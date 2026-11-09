#!/usr/bin/env python3
"""Familywise wild-cluster max-T audit of the saved 29 triple-overlap relations.

Uses only the hash-pinned Astrafek triple-overlap table. No HSC access occurs.
Inference is exploratory and is never promoted to a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCE_DIR = DATA / "icl_satellite_wl_triple_overlap_2026-11-06_01"
LENSES = SOURCE_DIR / "triple_overlap_lenses.tsv"
METRICS = SOURCE_DIR / "metrics.tsv"
SUMMARY = SOURCE_DIR / "summary.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_wl_triple_overlap_maxt_2026-11-06_01"
EXPECTED = {
    LENSES: "2eedc5d349c4ad2ecc4b95c32c18925c8dab18380ef6e77d3d08bbfe34c94c25",
    METRICS: "afe6fc396119bcb2b6c23af1afda3ae79d8080ecea5b75cc71800c1c6df34b8d",
    SUMMARY: "105f010d936850dbc137408c357e42c96dfc3f8e36e2e0a67fdf321bc722e469",
}
REPLICATES = 10000
SEED = 20261107
WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5),
                 math.sqrt(0.5), 1.0, math.sqrt(1.5)])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_r(x: np.ndarray, y: np.ndarray, c: np.ndarray) -> float:
    bx = np.linalg.lstsq(c, x, rcond=None)[0]
    by = np.linalg.lstsq(c, y, rcond=None)[0]
    rx, ry = x - c @ bx, y - c @ by
    den = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / den) if den > 0 else float("nan")


def fit_t(design: np.ndarray, y: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """Return last-column OLS coefficient and CR1 cluster-robust t statistic."""
    n, p = design.shape
    if n <= p or np.linalg.matrix_rank(design) != p:
        return float("nan"), float("nan")
    xtx_inv = np.linalg.inv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    residual = y - design @ beta
    blocks = np.unique(labels)
    g = len(blocks)
    if g < 6:
        return float(beta[-1]), float("nan")
    meat = np.zeros((p, p), dtype=float)
    for block in blocks:
        take = labels == block
        score = design[take].T @ residual[take]
        meat += np.outer(score, score)
    correction = (g / (g - 1.0)) * ((n - 1.0) / (n - p))
    variance = correction * xtx_inv @ meat @ xtx_inv
    se = math.sqrt(max(float(variance[-1, -1]), 0.0))
    return float(beta[-1]), float(beta[-1] / se) if se > 0 else float("nan")


def metric_specs(rows: list[dict[str, str]]) -> list[dict]:
    specs = [{
        "relationship": "icl_vs_satellite", "component": "na", "radial_bin": "na",
        "x": "icl_delta_200_300", "y": "satellite_delta_50_75",
        "covars": ["redshift", "log_richness"],
        "include": np.ones(len(rows), dtype=bool),
    }]
    for radial_bin in range(7):
        for component, ykey in (("tangential", f"wl_tangential_bin{radial_bin}_msun_mpc2"),
                                ("cross", f"wl_cross_bin{radial_bin}_msun_mpc2")):
            target_key, random_key = f"wl_target_pairs_bin{radial_bin}", f"wl_random_pairs_bin{radial_bin}"
            use = np.array([
                all(math.isfinite(float(r[k])) for k in
                    [ykey, "icl_delta_200_300", "satellite_delta_50_75", "redshift", "log_richness", target_key, random_key])
                and float(r[target_key]) > 0 and float(r[random_key]) > 0
                for r in rows
            ], dtype=bool)
            for predictor, xkey, other in (
                ("icl", "icl_delta_200_300", "satellite_delta_50_75"),
                ("satellite", "satellite_delta_50_75", "icl_delta_200_300"),
            ):
                specs.append({
                    "relationship": f"wl_{predictor}", "component": component,
                    "radial_bin": radial_bin, "x": xkey, "y": ykey,
                    "covars": [other, "redshift", "log_richness"], "include": use,
                })
    return specs


def prepare_model(rows: list[dict[str, str]], spec: dict, block_ids: np.ndarray) -> dict:
    indices = np.flatnonzero(spec["include"])
    chosen = [rows[i] for i in indices]
    x = np.asarray([float(r[spec["x"]]) for r in chosen])
    y = np.asarray([float(r[spec["y"]]) for r in chosen])
    covars = np.column_stack([[float(r[k]) for r in chosen] for k in spec["covars"]])
    reduced = np.column_stack((np.ones(len(indices)), covars))
    full = np.column_stack((reduced, x))
    labels = block_ids[indices]
    point_r = partial_r(x, y, reduced)
    beta, observed_t = fit_t(full, y, labels)
    reduced_fit = reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]
    restricted_resid = y - reduced_fit
    return {**spec, "indices": indices, "reduced_fit": reduced_fit,
            "restricted_resid": restricted_resid, "full": full, "labels": labels,
            "n": len(indices), "n_blocks": len(np.unique(labels)),
            "partial_r": point_r, "beta": beta, "observed_t": observed_t}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lenses", type=Path, default=LENSES)
    parser.add_argument("--metrics", type=Path, default=METRICS)
    parser.add_argument("--source-summary", type=Path, default=SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"Refusing to overwrite existing result: {args.output}")
    if args.replicates < 1000:
        parser.error("At least 1,000 Webb wild-cluster replicates are required")

    sources = {"lenses": args.lenses, "metrics": args.metrics, "source_summary": args.source_summary}
    hashes = {name: sha256(path) for name, path in sources.items()}
    for name, path in sources.items():
        pinned = {"lenses": LENSES, "metrics": METRICS, "source_summary": SUMMARY}[name]
        if path.resolve() == pinned.resolve() and hashes[name] != EXPECTED[path.resolve()]:
            raise ValueError(f"Pinned input hash mismatch: {path}")
    source_summary = json.loads(args.source_summary.read_text())
    if source_summary.get("science_measurement") is not False or source_summary.get("n_reported_relationships") != 29:
        raise ValueError("Source summary scope/status changed")
    rows = read_tsv(args.lenses)
    if len(rows) != 27 or len({r["cluster_id"] for r in rows}) != 27:
        raise ValueError("Expected 27 unique triple-overlap lenses")

    coords_t = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])} for r in rows]
    coords_c = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])} for r in rows]
    _, block_ids, _, _, n_blocks, largest = pair_block_labels(
        coords_t, coords_c, np.ones(len(rows), dtype=bool))
    models = [prepare_model(rows, spec, block_ids) for spec in metric_specs(rows)]
    if len(models) != 29 or any(m["n_blocks"] < 6 for m in models):
        raise ValueError("Expected 29 estimable relationships with at least six shared spatial blocks")
    if any(not math.isfinite(m["observed_t"]) for m in models):
        raise ValueError("A point estimate has a non-finite cluster-robust t statistic")

    rng = np.random.default_rng(SEED)
    block_values = np.unique(block_ids)
    t_draws = np.full((args.replicates, len(models)), np.nan)
    # Restricted wild-cluster bootstrap under each predictor's zero-slope null.
    # One shared Webb multiplier per spatial block per draw preserves dependence
    # across the 29 relationships, including radial-bin subsets.
    for b in range(args.replicates):
        multiplier = {block: WEBB[rng.integers(0, len(WEBB))] for block in block_values}
        for j, model in enumerate(models):
            y_star = model["reduced_fit"] + model["restricted_resid"] * np.asarray(
                [multiplier[label] for label in model["labels"]])
            _, t_value = fit_t(model["full"], y_star, model["labels"])
            t_draws[b, j] = t_value
    valid = np.isfinite(t_draws).all(axis=1)
    valid_draws = t_draws[valid]
    if valid_draws.shape[0] < int(0.99 * args.replicates):
        raise RuntimeError(f"Only {valid_draws.shape[0]}/{args.replicates} complete max-T draws")
    max_abs = np.max(np.abs(valid_draws), axis=1)

    output_rows = []
    for j, model in enumerate(models):
        abs_draw = np.abs(valid_draws[:, j])
        obs = abs(model["observed_t"])
        p_raw = (1 + int(np.count_nonzero(abs_draw >= obs))) / (1 + len(abs_draw))
        p_max = (1 + int(np.count_nonzero(max_abs >= obs))) / (1 + len(max_abs))
        output_rows.append({
            "relationship": model["relationship"], "component": model["component"],
            "radial_bin": model["radial_bin"], "n_lenses": model["n"],
            "n_shared_blocks": model["n_blocks"], "partial_r": model["partial_r"],
            "ols_slope": model["beta"], "cr1_t": model["observed_t"],
            "wild_cluster_p_unadjusted": p_raw, "wild_cluster_p_max_t_29": p_max,
            "bootstrap_replicates": len(valid_draws),
        })

    args.output.mkdir(parents=True)
    out_table = args.output / "max_t_metrics.tsv"
    with out_table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    report = {
        "schema": "astrafek-icl-satellite-wl-triple-overlap-max-t-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "network_used": False, "hsc_query_submitted": False,
        "family_size": len(models), "n_triple_overlap_lenses": len(rows),
        "shared_connected_5deg_blocks": n_blocks, "largest_block_n": largest,
        "method": "restricted Webb six-point wild-cluster bootstrap under each zero-slope null; CR1 studentized statistics; one shared multiplier per 5-degree block and replicate; max absolute statistic across 29 relationships",
        "replicates_requested": args.replicates, "replicates_valid": int(valid.sum()),
        "seed": SEED,
        "selection": "same audited 27-lens i2 common-support sample and per-bin target/random-pair support rules as the source triple-overlap analysis",
        "source_input_sha256": {str(path): hashes[name] for name, path in sources.items()},
        "output_sha256": {out_table.name: sha256(out_table)},
        "limitations": [
            "selected 27-lens common-support sample; not population representative",
            "5-degree equal-area connected components are approximate blocks, not validated HSC survey fields",
            "small numbers of blocks and familywise inference remain exploratory",
            "inherits unresolved weak-lensing source redshift, shear calibration, magnification, masking, PSF, and covariance systematics",
            "does not correct ICL fluxes for sky transfer, source Poisson, PSF wings, or image covariance",
        ],
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    smallest = sorted(output_rows, key=lambda r: float(r["wild_cluster_p_max_t_29"]))[:5]
    print(json.dumps({"status": report["status"], "science_measurement": False,
                      "valid_replicates": int(valid.sum()),
                      "smallest_adjusted_p": [{k: r[k] for k in ("relationship", "component", "radial_bin", "partial_r", "wild_cluster_p_unadjusted", "wild_cluster_p_max_t_29")} for r in smallest]}, indent=2))


if __name__ == "__main__":
    main()
