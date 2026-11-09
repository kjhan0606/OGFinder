#!/usr/bin/env python3
"""Compare continuous and median-split ICL relations across satellite radii.

Both association families are evaluated jointly on the existing supported
cluster sample. Restricted Webb wild-cluster bootstrap weights are shared
across all 28 filter/radius/method contrasts for familywise inference.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from familywise_block_wild_association import WEBB_WEIGHTS, holm_adjust
from fetch_ancillary_catalogues import sha


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOCIATION = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
DEFAULT_SATELLITE = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "icl_satellite_radial_relations_joint_2026-11-06_01"
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
              bread: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(labels)
    n, k = design.shape
    if n <= k or len(unique) < 2:
        return math.nan
    scores = np.vstack([design[labels == group].T @ residual[labels == group]
                        for group in unique])
    correction = (len(unique) / (len(unique) - 1)) * ((n - 1) / (n - k))
    covariance = correction * bread @ (scores.T @ scores) @ bread
    variance = covariance[-1, -1]
    if not np.isfinite(variance) or variance <= 0:
        return math.nan
    return float(beta[-1] / math.sqrt(variance))


def prepare(rows: list[dict], mode: str, icl_cut: float) -> dict:
    y = np.asarray([float(row["satellite_sb"]) for row in rows])
    icl = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    zsd, rlog = float(z.std()), np.log(richness)
    rsd, xsd = float(rlog.std()), float(icl.std())
    if min(zsd, rsd, xsd) <= 0:
        raise ValueError("predictor or adjustment covariate has zero variance")
    z = (z - z.mean()) / zsd
    rlog = (rlog - rlog.mean()) / rsd
    nuisance = np.column_stack((np.ones(len(rows)), z, rlog))
    if mode == "continuous_icl":
        signal = (icl - icl.mean()) / xsd
    elif mode == "median_split":
        signal = (icl >= icl_cut).astype(float)
    else:
        raise ValueError(f"unknown relation mode: {mode}")
    design = np.column_stack((nuisance, signal))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("ICL relation model is rank deficient")
    bread = np.linalg.inv(design.T @ design)
    beta = bread @ design.T @ y
    residual = y - design @ beta
    null_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_fit = nuisance @ (null_bread @ nuisance.T @ y)
    null_residual = y - null_fit
    leverage = np.sum((nuisance @ null_bread) * nuisance, axis=1)
    adjusted = null_residual / np.sqrt(np.maximum(1.0 - leverage, 0.05))
    blocks = np.asarray([str(row["global_connected_5deg_block"]) for row in rows])
    observed_t = cluster_t(beta, residual, design, bread, blocks)
    if not math.isfinite(observed_t):
        raise ValueError("invalid observed spatial-cluster statistic")
    groups = (icl >= icl_cut) if mode == "median_split" else None
    return {
        "rows": rows, "mode": mode, "design": design, "bread": bread,
        "beta": beta, "blocks": blocks, "null_fit": null_fit,
        "adjusted_null_residual": adjusted, "t_observed": observed_t,
        "estimate": float(beta[-1]),
        "se": float(abs(beta[-1] / observed_t)) if observed_t else math.nan,
        "n_clusters": len(rows), "n_blocks": len(np.unique(blocks)),
        "largest_block": int(max(np.count_nonzero(blocks == b) for b in np.unique(blocks))),
        "low_n": int(np.count_nonzero(~groups)) if groups is not None else None,
        "high_n": int(np.count_nonzero(groups)) if groups is not None else None,
        "low_mean": float(np.mean(y[~groups])) if groups is not None else None,
        "high_mean": float(np.mean(y[groups])) if groups is not None else None,
    }


def run(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output: {output}")
    association_dir, satellite_dir = args.association_dir.resolve(), args.satellite_dir.resolve()
    association_path = association_dir / "supported_pairs.tsv"
    association_summary_path = association_dir / "summary.json"
    satellite_path = satellite_dir / "cluster_profiles.tsv"
    satellite_summary_path = satellite_dir / "summary.json"
    a_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    s_summary = json.loads(satellite_summary_path.read_text(encoding="utf-8"))
    if (a_summary.get("status") != "complete_diagnostic"
            or a_summary.get("science_measurement") is not False
            or a_summary.get("n_pairs") != 94
            or a_summary.get("output_sha256", {}).get(str(association_path)) != sha(association_path)
            or s_summary.get("status") != "complete"
            or s_summary.get("science_measurement") is not False
            or s_summary.get("output_sha256", {}).get(str(satellite_path)) != sha(satellite_path)):
        raise ValueError("input summary status or pinned hashes do not match")
    association, profiles = read_tsv(association_path), read_tsv(satellite_path)
    if len(association) != 94 or len({(r["cluster_id"], r["filter"]) for r in association}) != 94:
        raise ValueError("expected 94 unique cluster/filter association pairs")
    targets = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])}
               for r in association]
    controls = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])}
                for r in association]
    indices, labels, _, _, _, _ = pair_block_labels(
        targets, controls, np.ones(len(association), dtype=bool))
    if len(indices) != len(association):
        raise ValueError("joint spatial block labeling changed row order")
    for row, label in zip(association, labels):
        row["global_connected_5deg_block"] = str(label)

    profile_lookup = {}
    for row in profiles:
        if row["cohort"] != "full_sample":
            continue
        key = (row["cluster_id"], row["filter"], row["radial_bin_kpc"])
        if key in profile_lookup:
            raise ValueError(f"duplicate cluster profile: {key}")
        profile_lookup[key] = row
    cuts = {band: float(np.median([float(r["icl_delta_200_300"])
                                   for r in association if r["filter"] == band]))
            for band in FILTERS}
    prepared = {}
    for band in FILTERS:
        band_pairs = [r for r in association if r["filter"] == band]
        for radial in RADIAL_BINS:
            rows = []
            for pair in band_pairs:
                profile = profile_lookup.get((pair["cluster_id"], band, radial))
                if profile is None:
                    continue
                row = dict(pair)
                row["satellite_sb"] = float(profile["mean_delta_sb_njy_arcsec2"])
                rows.append(row)
            if not rows:
                raise ValueError(f"no supported profile at {band}/{radial}")
            for mode in ("median_split", "continuous_icl"):
                group = prepare(rows, mode, cuts[band])
                if (group["n_blocks"] < args.minimum_blocks
                        or (mode == "median_split" and min(group["high_n"], group["low_n"]) < 3)):
                    raise ValueError(f"insufficient cluster/block support for {mode}, {band}/{radial}")
                prepared[(mode, band, radial)] = group

    global_blocks = sorted({str(row["global_connected_5deg_block"]) for row in association})
    rng = np.random.default_rng(args.seed)
    bootstrap = {key: np.empty(args.replicas) for key in prepared}
    max_t = np.empty(args.replicas)
    for b in range(args.replicas):
        weights = dict(zip(global_blocks, rng.choice(WEBB_WEIGHTS, len(global_blocks), replace=True)))
        maximum = 0.0
        for key, group in prepared.items():
            w = np.asarray([weights[block] for block in group["blocks"]])
            y_star = group["null_fit"] + group["adjusted_null_residual"] * w
            design, bread = group["design"], group["bread"]
            beta_star = bread @ design.T @ y_star
            residual_star = y_star - design @ beta_star
            statistic = cluster_t(beta_star, residual_star, design, bread, group["blocks"])
            if not math.isfinite(statistic):
                raise ValueError(f"invalid bootstrap statistic for {key}")
            bootstrap[key][b] = statistic
            maximum = max(maximum, abs(statistic))
        max_t[b] = maximum

    raw_p = {key: float((1 + np.count_nonzero(np.abs(draw) >= abs(prepared[key]["t_observed"])))
                        / (args.replicas + 1)) for key, draw in bootstrap.items()}
    holm = dict(zip(prepared, map(float, holm_adjust(np.asarray(list(raw_p.values()))))))
    records = []
    for key, group in prepared.items():
        mode, band, radial = key
        critical = float(student_t.ppf(0.975, group["n_blocks"] - 1))
        lo, hi = group["estimate"] - critical * group["se"], group["estimate"] + critical * group["se"]
        records.append({
            "relation": mode, "filter": band, "radial_bin_kpc": radial,
            "n_clusters": group["n_clusters"], "n_spatial_blocks": group["n_blocks"],
            "largest_spatial_block": group["largest_block"],
            "low_icl_n": group["low_n"], "high_icl_n": group["high_n"],
            "low_icl_mean_satellite_excess_njy_arcsec2": group["low_mean"],
            "high_icl_mean_satellite_excess_njy_arcsec2": group["high_mean"],
            "adjusted_effect_njy_arcsec2": group["estimate"],
            "effect_definition": ("high-minus-low at median ICL" if mode == "median_split"
                                  else "per 1 SD higher ICL signal"),
            "cr1_se_njy_arcsec2": group["se"], "cr1_t95_low": lo, "cr1_t95_high": hi,
            "cr1_t": group["t_observed"], "p_webb_raw": raw_p[key],
            "p_holm_across_28": holm[key],
            "p_maxT_across_28": float((1 + np.count_nonzero(max_t >= abs(group["t_observed"])))
                                       / (args.replicas + 1)),
        })
    output.mkdir(parents=True, exist_ok=False)
    table = output / "icl_satellite_radial_relation_profiles.tsv"
    write_tsv(table, records, list(records[0]))
    null_table = output / "maxT_null.tsv"
    write_tsv(null_table, [{"replica": i + 1, "max_abs_t_across_28": float(v)}
                           for i, v in enumerate(max_t)], ["replica", "max_abs_t_across_28"])
    summary = {
        "schema": "astrafek-icl-satellite-radial-relations-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "hsc_network_used": False, "image_pixels_read": False,
        "n_association_pairs": len(association), "n_by_filter": {
            band: sum(r["filter"] == band for r in association) for band in FILTERS},
        "n_relation_tests": len(records), "replicas": args.replicas, "seed": args.seed,
        "methods": {
            "continuous": "per-1-SD ICL signal slope at each radius, adjusted for standardized redshift and log richness",
            "median_split": "high-minus-low contrast at filter-specific median ICL signal, adjusted for standardized redshift and log richness",
            "covariance": "CR1 cluster-robust covariance with approximate 5-degree connected target/control pair blocks",
            "bootstrap": "restricted HC2 Webb six-point wild-cluster bootstrap-t; same joint spatial-block weights across all 28 tests",
            "familywise": "Holm and single-step max-T correction across both methods, both filters, and all seven radii (28 tests)",
        },
        "icl_split_median_njy_arcsec2": cuts,
        "input_sha256": {str(p): sha(p) for p in
            (association_path, association_summary_path, satellite_path, satellite_summary_path)},
        "output_sha256": {str(p): sha(p) for p in (table, null_table)},
        "runner_sha256": sha(Path(__file__).resolve()),
        "limitations": [
            "94 pairs are a selected ICL/satellite support intersection and are not representative validation of all 800 cutouts",
            "both relation families are exploratory and post-hoc; 28-way correction does not undo prior research selection",
            "connected equal-area 5-degree cells approximate spatial dependence, not validated HSC survey fields",
            "intervals omit photometric covariance, source Poisson, coadd-weight, sky-transfer, PSF-wing, completeness, and injection-response uncertainty",
            "per-radius support varies and missing cluster profiles are not imputed",
        ],
    }
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "tests": len(records), "replicas": args.replicas,
                      "output": str(output / "summary.json")}, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-dir", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--minimum-blocks", type=int, default=6)
    args = parser.parse_args()
    if args.replicas < 1_000:
        parser.error("--replicas must be at least 1,000")
    if args.minimum_blocks < 6:
        parser.error("--minimum-blocks must be at least 6")
    run(args)


if __name__ == "__main__":
    main()
