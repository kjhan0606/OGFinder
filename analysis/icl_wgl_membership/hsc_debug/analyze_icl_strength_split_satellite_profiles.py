#!/usr/bin/env python3
"""Stack satellite-centered profiles by cluster ICL strength.

This exploratory analysis uses the existing 94 paired ICL/satellite sample.
It compares high- and low-ICL cluster-level satellite excess profiles at each
radial bin, adjusts the high-low contrast for redshift and richness, and uses
a common restricted Webb wild-cluster bootstrap across filters and radii.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import sha
from evaluate_prestack_pair_support import pair_block_labels
from familywise_block_wild_association import WEBB_WEIGHTS, holm_adjust


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOCIATION = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
DEFAULT_SATELLITE = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "icl_satellite_icl_strength_stacked_profiles_2026-11-06_01"
RADIAL_BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
FILTERS = ("i", "i2")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def cluster_t(beta: np.ndarray, residual: np.ndarray, design: np.ndarray,
              bread: np.ndarray, labels: np.ndarray, coefficient: int = 3) -> float:
    groups = np.unique(labels)
    n, k = design.shape
    if len(groups) < 2 or n <= k:
        return math.nan
    scores = np.vstack([design[labels == group].T @ residual[labels == group]
                        for group in groups])
    correction = (len(groups) / (len(groups) - 1)) * ((n - 1) / (n - k))
    covariance = correction * (bread @ (scores.T @ scores) @ bread)
    variance = covariance[coefficient, coefficient]
    if not np.isfinite(variance) or variance <= 0:
        return math.nan
    return float(beta[coefficient] / math.sqrt(variance))


def prepare(rows: list[dict], icl_cut: float) -> dict:
    y = np.asarray([float(row["satellite_sb"]) for row in rows])
    x_icl = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    high = (x_icl >= icl_cut).astype(float)
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    zsd = float(z.std(ddof=0))
    rsd = float(np.log(richness).std(ddof=0))
    if zsd <= 0 or rsd <= 0:
        raise ValueError("redshift or richness has no variance in a split profile")
    z = (z - z.mean()) / zsd
    richness = (np.log(richness) - np.log(richness).mean()) / rsd
    nuisance = np.column_stack((np.ones(len(rows)), z, richness))
    design = np.column_stack((nuisance, high))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("high-low profile regression is rank deficient")
    bread = np.linalg.inv(design.T @ design)
    beta = bread @ design.T @ y
    residual = y - design @ beta
    null_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_fit = nuisance @ (null_bread @ nuisance.T @ y)
    null_resid = y - null_fit
    leverage = np.sum((nuisance @ null_bread) * nuisance, axis=1)
    adjusted = null_resid / np.sqrt(np.maximum(1.0 - leverage, 0.05))
    blocks = np.asarray([row["connected_5deg_block"] for row in rows], dtype=str)
    t_obs = cluster_t(beta, residual, design, bread, blocks)
    if not math.isfinite(t_obs):
        raise ValueError("cluster-robust high-low statistic is invalid")
    return {
        "rows": rows, "y": y, "high": high, "design": design, "bread": bread,
        "beta": beta, "residual": residual, "blocks": blocks,
        "null_fit": null_fit, "adjusted_null_residual": adjusted,
        "t_observed": t_obs, "icl_cut": icl_cut,
        "low_mean": float(np.mean(y[high == 0])),
        "high_mean": float(np.mean(y[high == 1])),
        "low_n": int(np.count_nonzero(high == 0)),
        "high_n": int(np.count_nonzero(high == 1)),
        "n_blocks": int(len(np.unique(blocks))),
        "largest_block": int(max(np.count_nonzero(blocks == b) for b in np.unique(blocks))),
        "n_clusters": len(rows),
    }


def run(args):
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    association_summary_path = args.association_dir / "summary.json"
    association_path = args.association_dir / "supported_pairs.tsv"
    satellite_summary_path = args.satellite_dir / "summary.json"
    satellite_path = args.satellite_dir / "cluster_profiles.tsv"
    association_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    satellite_summary = json.loads(satellite_summary_path.read_text(encoding="utf-8"))
    if (association_summary.get("status") != "complete_diagnostic"
            or association_summary.get("science_measurement") is not False
            or association_summary.get("n_pairs") != 94
            or association_summary.get("output_sha256", {}).get(str(association_path)) != sha(association_path)
            or satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("output_sha256", {}).get(str(satellite_path)) != sha(satellite_path)):
        raise ValueError("association or satellite-profile input fails its saved hash/status audit")

    association = read_tsv(association_path)
    satellite = read_tsv(satellite_path)
    if len(association) != 94 or len({(r["cluster_id"], r["filter"]) for r in association}) != 94:
        raise ValueError("expected 94 unique cluster/filter association pairs")
    targets = [{"ra_deg": float(row["target_ra_deg"]),
                "dec_deg": float(row["target_dec_deg"])} for row in association]
    controls = [{"ra_deg": float(row["control_ra_deg"]),
                 "dec_deg": float(row["control_dec_deg"])} for row in association]
    pair_indices, global_labels, _, _, _, _ = pair_block_labels(
        targets, controls, np.ones(len(association), dtype=bool))
    if len(pair_indices) != len(association):
        raise ValueError("global spatial-block construction changed pair order")
    for row, label in zip(association, global_labels):
        row["global_connected_5deg_block"] = str(label)
    profile_lookup = {}
    for row in satellite:
        if row["cohort"] != "full_sample":
            continue
        key = (str(row["cluster_id"]), row["filter"], row["radial_bin_kpc"])
        if key in profile_lookup:
            raise ValueError(f"duplicate per-bin satellite profile: {key}")
        profile_lookup[key] = row

    by_filter = {band: [row for row in association if row["filter"] == band]
                 for band in FILTERS}
    cuts = {
        band: float(np.median([float(row["icl_delta_200_300"]) for row in by_filter[band]]))
        for band in FILTERS
    }
    records = []
    prepared = {}
    for band in FILTERS:
        sample = by_filter[band]
        for radial_bin in RADIAL_BINS:
            rows = []
            for pair in sample:
                profile = profile_lookup.get((pair["cluster_id"], band, radial_bin))
                if profile is None:
                    continue
                record = dict(pair)
                record["satellite_sb"] = float(profile["mean_delta_sb_njy_arcsec2"])
                record["n_satellites"] = int(profile["n_satellites"])
                rows.append(record)
            if not rows:
                raise ValueError(f"no common ICL/satellite support for {band} {radial_bin}")
            group = prepare(rows, cuts[band])
            if (min(group["high_n"], group["low_n"]) < 3
                    or group["n_blocks"] < args.minimum_blocks):
                raise ValueError(f"insufficient support for {band} {radial_bin}: {group['n_clusters']} rows, {group['n_blocks']} blocks")
            key = (band, radial_bin)
            prepared[key] = group

    global_blocks = sorted({str(row["global_connected_5deg_block"])
                            for row in association})
    rng = np.random.default_rng(args.seed)
    t_draws = {key: np.empty(args.replicas) for key in prepared}
    max_t_draws = np.empty(args.replicas)
    for replica in range(args.replicas):
        block_weights = dict(zip(global_blocks, rng.choice(
            WEBB_WEIGHTS, size=len(global_blocks), replace=True,
        )))
        maximum = 0.0
        for key, group in prepared.items():
            weight = np.asarray([block_weights[block] for block in group["blocks"]])
            y_star = group["null_fit"] + group["adjusted_null_residual"] * weight
            design, bread = group["design"], group["bread"]
            beta_star = bread @ design.T @ y_star
            residual_star = y_star - design @ beta_star
            value = cluster_t(beta_star, residual_star, design, bread, group["blocks"])
            if not math.isfinite(value):
                raise ValueError(f"non-finite restricted-bootstrap draw for {key}")
            t_draws[key][replica] = value
            maximum = max(maximum, abs(value))
        max_t_draws[replica] = maximum

    raw_p = {}
    max_t_p = {}
    for key, group in prepared.items():
        raw_p[key] = float((1 + np.count_nonzero(
            np.abs(t_draws[key]) >= abs(group["t_observed"]))) / (args.replicas + 1))
        max_t_p[key] = float((1 + np.count_nonzero(
            max_t_draws >= abs(group["t_observed"]))) / (args.replicas + 1))
    holm_values = holm_adjust(np.asarray([raw_p[key] for key in prepared]))
    holm_p = dict(zip(prepared, map(float, holm_values)))

    for key, group in prepared.items():
        beta = group["beta"]
        t_observed = group["t_observed"]
        blocks = group["n_blocks"]
        # Design order is intercept, redshift, log richness, high-ICL.
        high_low = float(beta[3])
        se = float(high_low / t_observed) if t_observed else None
        critical = float(student_t.ppf(0.975, blocks - 1))
        records.append({
            "filter": key[0], "radial_bin_kpc": key[1],
            "icl_split_median_njy_arcsec2": group["icl_cut"],
            "n_clusters": group["n_clusters"], "n_spatial_blocks": blocks,
            "largest_spatial_block": group["largest_block"],
            "low_icl_n": group["low_n"], "high_icl_n": group["high_n"],
            "low_icl_mean_satellite_excess_njy_arcsec2": group["low_mean"],
            "high_icl_mean_satellite_excess_njy_arcsec2": group["high_mean"],
            "adjusted_high_minus_low_njy_arcsec2": high_low,
            "cr1_standard_error_njy_arcsec2": abs(se) if se is not None else None,
            "cr1_t": t_observed,
            "cr1_t95_low": float(high_low - critical * abs(se)),
            "cr1_t95_high": float(high_low + critical * abs(se)),
            "p_restricted_webb_wild_raw": raw_p[key],
            "p_holm_across_14_profiles": holm_p[key],
            "p_maxT_across_14_profiles": max_t_p[key],
        })

    output.mkdir(parents=True, exist_ok=False)
    profile_path = output / "icl_strength_split_satellite_profile.tsv"
    fields = list(records[0])
    write_tsv(profile_path, records, fields)
    null_path = output / "maxT_null.tsv"
    write_tsv(null_path,
              [{"replica": i + 1, "max_abs_t": float(value)}
               for i, value in enumerate(max_t_draws)],
              ["replica", "max_abs_t"])
    summary = {
        "schema": "astrafek-icl-strength-split-satellite-profile-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "method": {
            "estimand": "equal-cluster mean target-minus-azimuthal-control satellite excess at each available radial bin",
            "split": "within-filter median ICL 200-300 kpc signal from the fixed 94-pair association table; high group includes ties",
            "adjustment": "OLS high-minus-low contrast adjusted for standardized cluster redshift and log richness",
            "resampling": "restricted HC2 Webb six-point wild-cluster bootstrap-t with common 5-degree connected pair blocks recomputed jointly across both filters, shared across all 14 filter/radial-bin tests",
            "multiplicity": "Holm and single-step max-T across all 14 profile points",
            "profile_selection": "available per-annulus cluster/filter satellite profiles joined to the fixed 94-pair ICL sample; N varies by radius and no profile is imputed",
        },
        "n_association_pairs": len(association),
        "n_by_filter": {band: len(by_filter[band]) for band in FILTERS},
        "icl_split_median_by_filter_njy_arcsec2": cuts,
        "n_profile_tests": len(records),
        "replicas": args.replicas,
        "seed": args.seed,
        "profile_table": str(profile_path),
        "input_sha256": {
            str(association_path): sha(association_path),
            str(association_summary_path): sha(association_summary_path),
            str(satellite_path): sha(satellite_path),
            str(satellite_summary_path): sha(satellite_summary_path),
        },
        "output_sha256": {str(path): sha(path) for path in (profile_path, null_path)},
        "runner_sha256": sha(Path(__file__).resolve()),
        "limitations": [
            "the 94-pair sample was selected by existing ICL, satellite, mask, and control support gates",
            "the high/low ICL split is a descriptive sensitivity rather than a prespecified independent validation",
            "the 5-degree blocks approximate spatial dependence and do not represent validated HSC survey fields",
            "bootstrap intervals omit photometric covariance, source Poisson, sky transfer, PSF-wing, completeness, and transfer-matrix uncertainty",
            "multiplicity adjustment covers the 14 filter-by-radius profile contrasts only and does not undo prior analysis selection",
        ],
    }
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "profile_tests": len(records), "replicas": args.replicas,
        "output": str(output / "summary.json"),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-dir", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--minimum-blocks", type=int, default=6)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("--replicas must be at least 1,000")
    if args.minimum_blocks < 6:
        parser.error("--minimum-blocks must be at least 6")
    run(args)


if __name__ == "__main__":
    main()
