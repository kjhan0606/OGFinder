#!/usr/bin/env python3
"""Test ICL--satellite correlation sensitivity to nonlinear z/richness controls.

This offline, post-hoc sensitivity keeps the established 94-pair sample and
adds quadratic/interaction terms to its linear redshift and CAMIRA-richness
adjustment. It does not select a preferred model or promote a detection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_DIR = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
SUPPORT = SUPPORT_DIR / "supported_pairs.tsv"
SUPPORT_SHA256 = "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"
DEFAULT_OUTPUT = ROOT / "icl_satellite_nonlinear_adjustment_2026-11-06_01"
FILTERS = ("i", "i2")
MODELS = ("quadratic_redshift", "quadratic_richness", "interaction", "full_quadratic")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def prepare(rows: list[dict[str, str]]) -> list[dict]:
    output = []
    for source in rows:
        if source.get("science_measurement") != "False":
            raise ValueError("support input unexpectedly promotes science_measurement")
        row = {key: source[key] for key in (
            "cluster_id", "control_id", "filter", "icl_delta_200_300",
            "satellite_delta_50_75", "redshift", "richness",
            "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg")}
        row.update({key: float(row[key]) for key in (
            "icl_delta_200_300", "satellite_delta_50_75", "redshift", "richness",
            "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg")})
        output.append(row)
    if len(output) != 94 or {band: sum(row["filter"] == band for row in output)
                             for band in FILTERS} != {"i": 26, "i2": 68}:
        raise ValueError("expected the fixed 94-pair support (26 i, 68 i2)")
    return output


def make_design(rows: list[dict], model: str, indices: np.ndarray | None = None,
                centers: tuple[float, float, float, float] | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = rows if indices is None else [rows[int(i)] for i in indices]
    z = np.asarray([row["redshift"] for row in selected], dtype=float)
    richness = np.asarray([row["richness"] for row in selected], dtype=float)
    if centers is None:
        z0, zs = float(np.mean(z)), float(np.std(z, ddof=0))
        r0, rs = float(np.mean(richness)), float(np.std(richness, ddof=0))
    else:
        z0, zs, r0, rs = centers
    if min(zs, rs) <= 0:
        raise ValueError("redshift or richness has zero spread")
    zc, rc = (z - z0) / zs, (richness - r0) / rs
    controls = [z, richness]
    if model in ("quadratic_redshift", "full_quadratic"):
        controls.append(zc * zc)
    if model in ("quadratic_richness", "full_quadratic"):
        controls.append(rc * rc)
    if model in ("interaction", "full_quadratic"):
        controls.append(zc * rc)
    design = np.column_stack([np.ones(len(selected)), *controls])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError(f"rank-deficient adjustment design: {model}")
    x = np.asarray([row["icl_delta_200_300"] for row in selected], dtype=float)
    y = np.asarray([row["satellite_delta_50_75"] for row in selected], dtype=float)
    return x, y, design


def estimate(rows: list[dict], model: str | None,
             indices: np.ndarray | None = None,
             centers: tuple[float, float, float, float] | None = None) -> float:
    selected = rows if indices is None else [rows[int(i)] for i in indices]
    z = np.asarray([row["redshift"] for row in selected], dtype=float)
    richness = np.asarray([row["richness"] for row in selected], dtype=float)
    if centers is None:
        centers = (float(np.mean(z)), float(np.std(z, ddof=0)),
                   float(np.mean(richness)), float(np.std(richness, ddof=0)))
    effective_model = "baseline" if model is None else model
    x, y, design = make_design(selected, effective_model, centers=centers)
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / denom) if denom > 0 else float("nan")


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if not SUPPORT.is_file() or sha256(SUPPORT) != SUPPORT_SHA256:
        raise ValueError("pinned 94-pair support table is missing or changed")
    parent_summary = json.loads((SUPPORT_DIR / "summary.json").read_text(encoding="utf-8"))
    if parent_summary.get("science_measurement") is not False:
        raise ValueError("support parent must remain diagnostic-only")
    rows = prepare(read_rows(SUPPORT))

    results, draws, rows_by_filter, blocks_by_filter = {}, {}, {}, {}
    observed = {}
    for band_index, band in enumerate(FILTERS):
        selected = [row for row in rows if row["filter"] == band]
        target = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in selected]
        control = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in selected]
        _, labels, nt, nc, nblocks, largest = pair_block_labels(
            target, control, np.ones(len(selected), dtype=bool))
        labels = np.asarray(labels, dtype=str)
        groups = np.unique(labels)
        if len(groups) < 6:
            raise ValueError(f"at least six spatial blocks required; {band} has {len(groups)}")
        by_block = {group: np.flatnonzero(labels == group) for group in groups}
        z = np.asarray([row["redshift"] for row in selected])
        richness = np.asarray([row["richness"] for row in selected])
        centers = (float(np.mean(z)), float(np.std(z)),
                   float(np.mean(richness)), float(np.std(richness)))
        rows_by_filter[band] = (selected, centers)
        blocks_by_filter[band] = by_block
        results[band] = {"n_pairs": len(selected), "n_connected_5deg_blocks": len(groups),
                         "n_target_cells": int(nt), "n_control_cells": int(nc),
                         "largest_block_pairs": int(largest)}
        models = ("baseline", *MODELS)
        for model in models:
            observed[(band, model)] = estimate(selected, None if model == "baseline" else model,
                                               centers=centers)
            draws[(band, model)] = []
        rng = np.random.default_rng(seed + band_index * 1000)
        for _ in range(replicas):
            selected_groups = rng.choice(groups, size=len(groups), replace=True)
            indices = np.concatenate([by_block[group] for group in selected_groups])
            estimates = {model: estimate(selected, None if model == "baseline" else model,
                                         indices, centers) for model in models}
            if all(np.isfinite(value) for value in estimates.values()):
                for model, value in estimates.items():
                    draws[(band, model)].append(float(value))

    comparisons = [(band, model) for band in FILTERS for model in MODELS]
    delta_bootstrap, ses = {}, {}
    for key in comparisons:
        band, model = key
        base = np.asarray(draws[(band, "baseline")])
        adjusted = np.asarray(draws[(band, model)])
        valid = min(len(base), len(adjusted))
        if valid < replicas * 0.9:
            raise ValueError(f"too few valid bootstrap draws for {band}/{model}: {valid}/{replicas}")
        delta_bootstrap[key] = adjusted[:valid] - base[:valid]
        ses[key] = float(np.std(delta_bootstrap[key], ddof=1))
        if not np.isfinite(ses[key]) or ses[key] <= 0:
            raise ValueError(f"degenerate bootstrap change for {band}/{model}")

    max_t_draws = []
    for index in range(replicas):
        values = []
        if any(index >= len(delta_bootstrap[key]) for key in comparisons):
            continue
        for key in comparisons:
            band, model = key
            delta = observed[(band, model)] - observed[(band, "baseline")]
            values.append(abs((delta_bootstrap[key][index] - delta) / ses[key]))
        max_t_draws.append(max(values))
    if len(max_t_draws) < replicas * 0.9:
        raise ValueError(f"too few complete max-T resamples: {len(max_t_draws)}/{replicas}")
    critical = float(np.quantile(max_t_draws, 0.95))
    max_t_draws = np.asarray(max_t_draws)
    for band, model in comparisons:
        delta = observed[(band, model)] - observed[(band, "baseline")]
        sample = delta_bootstrap[(band, model)]
        results[band].setdefault("models", {})[model] = {
            "baseline_partial_r": observed[(band, "baseline")],
            "adjusted_partial_r": observed[(band, model)],
            "paired_change_in_r": delta,
            "paired_change_ci95_pointwise": np.percentile(sample, [2.5, 97.5]).tolist(),
            "paired_change_ci95_maxT": [delta - critical * ses[(band, model)],
                                         delta + critical * ses[(band, model)]],
            "maxT_adjusted_p": float((1 + np.count_nonzero(
                max_t_draws >= abs(delta / ses[(band, model)]))) / (len(max_t_draws) + 1)),
            "replicas_valid": len(sample),
        }

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "joined_supported_pairs.tsv"
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema": "astrafek-icl-satellite-nonlinear-adjustment-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "n_pairs": len(rows),
        "filter_counts": {band: sum(row["filter"] == band for row in rows) for band in FILTERS},
        "baseline_controls": ["redshift", "linear CAMIRA richness"],
        "additional_models": {
            "quadratic_redshift": "baseline plus centered/scaled redshift squared",
            "quadratic_richness": "baseline plus centered/scaled richness squared",
            "interaction": "baseline plus redshift-by-richness interaction",
            "full_quadratic": "baseline plus both squares and the interaction",
        },
        "resampling": {
            "method": "connected 5-degree target/control endpoint blocks, resampled within filter",
            "replicas_requested": replicas, "seed_base": seed,
            "family": "four nonlinear adjustments in each filter, eight paired changes",
            "maxT_critical_95": critical, "replicas_valid_joint": len(max_t_draws),
        },
        "results": results,
        "inputs_sha256": {str(SUPPORT): sha256(SUPPORT),
                           str(SUPPORT_DIR / "summary.json"): sha256(SUPPORT_DIR / "summary.json")},
        "runner_sha256": sha256(Path(__file__).resolve()),
        "joined_pairs_path": str(joined_path),
        "joined_pairs_sha256": sha256(joined_path),
        "limitations": [
            "This is a post-hoc model-form sensitivity on 26 i and 68 i2 selected pairs.",
            "The adjusted covariates are descriptive controls, not an identified causal model.",
            "The bootstrap blocks are approximate spatial groups, not validated HSC fields.",
            "The analysis does not address pixel noise, sky transfer, PSF, source Poisson, or selection.",
            "No profile correction or science_measurement promotion is made.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 spatial-block bootstrap replicates are required")
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "results": result["results"], "output": str(args.output.resolve()),
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
