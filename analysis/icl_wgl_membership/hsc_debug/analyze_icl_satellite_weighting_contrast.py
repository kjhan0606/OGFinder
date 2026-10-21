#!/usr/bin/env python3
"""Paired spatial-block bootstrap for CAMIRA-weight vs equal-satellite ICL association."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from analyze_icl_satellite_radial_matrix_familywise import (
    DEFAULT_EXTRACTION, DEFAULT_ICL, DEFAULT_SATELLITE, load_inputs,
)
from analyze_icl_satellite_weighted_radial_matrix_familywise import (
    DEFAULT_MEMBER, load_weighted_profiles,
)
from audit_icl_satellite_radial_support_matrix import sha256


PRIMARY_CELLS = (("i", "200-300", "50-75"), ("i2", "200-300", "50-75"))


def partial_r(rows: list[dict], satellite_field: str) -> float:
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[satellite_field]) for row in rows])
    z = np.asarray([float(row["z"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    if len(rows) < 5 or not all(np.isfinite(v).all() for v in (x, y, z, richness)):
        raise ValueError("partial correlation requires at least five finite rows")
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("redshift/richness design is rank deficient")
    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(xr) * np.linalg.norm(yr))
    if denominator <= 0:
        raise ValueError("partial correlation is undefined")
    return float(xr @ yr / denominator)


def paired_bootstrap(equal_rows: list[dict], weighted_rows: list[dict],
                     blocks: list[str], *, replicas: int, seed: int) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 bootstrap replicates are required")
    if len(equal_rows) != len(weighted_rows) or len(equal_rows) != len(blocks):
        raise ValueError("equal/weighted rows and blocks must be paired")
    for a, b in zip(equal_rows, weighted_rows):
        for field in ("cluster_id", "filter", "icl_radial_bin_kpc", "satellite_radial_bin_kpc"):
            if a[field] != b[field]:
                raise ValueError(f"equal/weighted pair key differs: {field}")
    groups = [np.flatnonzero(np.asarray(blocks) == label)
              for label in sorted(set(blocks))]
    if len(groups) < 3:
        raise ValueError("at least three spatial blocks are required")
    point_equal = partial_r(equal_rows, "satellite_delta_50_75")
    point_weighted = partial_r(weighted_rows, "satellite_delta_50_75")
    draws_equal = np.full(replicas, np.nan)
    draws_weighted = np.full(replicas, np.nan)
    draws_delta = np.full(replicas, np.nan)
    rng = np.random.default_rng(seed)
    for index in range(replicas):
        chosen = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[int(choice)] for choice in chosen])
        eq = [equal_rows[int(i)] for i in indices]
        wt = [weighted_rows[int(i)] for i in indices]
        try:
            draws_equal[index] = partial_r(eq, "satellite_delta_50_75")
            draws_weighted[index] = partial_r(wt, "satellite_delta_50_75")
            draws_delta[index] = draws_weighted[index] - draws_equal[index]
        except ValueError:
            continue
    valid = np.isfinite(draws_delta)
    if valid.sum() < 0.8 * replicas:
        raise ValueError(f"too few paired bootstrap draws: {valid.sum()}/{replicas}")

    def interval(values):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]

    return {
        "n_clusters": len(equal_rows), "n_spatial_blocks": len(groups),
        "equal_satellite_partial_r": point_equal,
        "equal_satellite_ci95_block": interval(draws_equal),
        "camira_weighted_partial_r": point_weighted,
        "camira_weighted_ci95_block": interval(draws_weighted),
        "weighted_minus_equal_partial_r": point_weighted - point_equal,
        "paired_delta_ci95_block": interval(draws_delta),
        "bootstrap_fraction_delta_below_zero": float(np.mean(draws_delta[valid] < 0)),
        "valid_replicates": int(valid.sum()), "replicas_requested": replicas,
        "seed": seed,
    }


def analyze(member_dir: Path, icl_dir: Path, satellite_dir: Path,
            extraction_path: Path, output_dir: Path, *, replicas: int, seed: int):
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output: {output_dir}")
    weighted, weight_provenance = load_weighted_profiles(member_dir)
    groups, support, base_provenance = load_inputs(icl_dir, satellite_dir, extraction_path)
    results = []
    family_input = []
    for band, icl_bin, satellite_bin in PRIMARY_CELLS:
        key = (band, icl_bin, satellite_bin)
        group = groups[key]
        equal_rows = group["rows"]
        weighted_rows = []
        for row in equal_rows:
            item = weighted.get((str(row["cluster_id"]), band, satellite_bin))
            if item is None or item["n_satellites"] != int(row["n_satellites"]):
                raise ValueError(f"weighted profile mismatch for {key}/{row['cluster_id']}")
            weighted_rows.append({**row, "satellite_delta_50_75": item["weighted_mean"]})
        labels = list(map(str, group["blocks"]))
        result = paired_bootstrap(equal_rows, weighted_rows, labels,
                                  replicas=replicas, seed=seed + (0 if band == "i" else 1))
        result.update({"filter": band, "icl_radial_bin_kpc": icl_bin,
                       "satellite_radial_bin_kpc": satellite_bin})
        results.append(result)
        family_input.append({
            "filter": band, "icl_radial_bin_kpc": icl_bin,
            "satellite_radial_bin_kpc": satellite_bin,
            "n_support_rows": next(r["n_common"] for r in support
                                    if (r["filter"], r["icl_radial_bin_kpc"],
                                        r["satellite_radial_bin_kpc"]) == key),
            "n_blocks": len(set(labels)),
        })
    output_dir.mkdir(parents=True, exist_ok=False)
    result_path = output_dir / "paired_weight_contrast.tsv"
    fields = list(results[0])
    with result_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    summary = {
        "schema": "astrafek-icl-satellite-camira-weight-paired-contrast-v1",
        "status": "complete_exploratory_paired_sensitivity",
        "science_measurement": False,
        "method": {
            "endpoint": "ICL 200-300 kpc vs satellite 50-75 kpc",
            "estimand": "CAMIRA-weighted partial Pearson r minus equal-satellite partial Pearson r, controlling redshift and richness",
            "resampling": "paired resampling of the same connected 5-degree target/control blocks; duplicate blocks retained",
            "interval": "central percentile 95% paired block-bootstrap interval; no null-test p-value",
            "interpretation": "descriptive estimator sensitivity only; endpoint was selected in prior exploratory work",
            "replicas": replicas, "seed": seed,
        },
        "provenance": {"icl_satellite": base_provenance,
                       "camira_weighted_profiles": weight_provenance},
        "support": family_input,
        "results": results,
        "output_sha256": {str(result_path.resolve()): sha256(result_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "two primary-filter contrasts are not a new confirmatory test and are not multiplicity-adjusted across prior endpoint selection",
            "CAMIRA weight includes its mass/spatial selection and is not an independent unbiased membership probability",
            "5-degree spatial blocks approximate dependence and do not include pixel/coadd covariance",
            "PSF wings, sky transfer, extinction, satellite completeness, and selection bias remain unresolved",
            "science_measurement remains false",
        ],
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--member-dir", type=Path, default=DEFAULT_MEMBER)
    parser.add_argument("--icl-dir", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261023)
    args = parser.parse_args()
    result = analyze(args.member_dir, args.icl_dir, args.satellite_dir,
                     args.extraction_summary, args.output_dir,
                     replicas=args.replicas, seed=args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "results": result["results"], "output": str(args.output_dir.resolve())},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
