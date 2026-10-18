#!/usr/bin/env python3
"""Exploratory cluster-level link between satellite envelopes and ICL light.

Joins equal-satellite-weighted 50--75 kpc satellite outer-light profiles to the
existing 50--300 kpc K2.5/c40 cluster-light pair measurements. Filters remain
separate. This is a selected-sample association diagnostic, not a science
measurement or a physical-causality test.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def partial_pearson(x, y, redshift, richness):
    covariates = np.column_stack((np.asarray(redshift, float), np.asarray(richness, float)))
    scales = np.std(covariates, axis=0)
    if (not np.isfinite(covariates).all() or np.any(scales <= 0)
            or np.std(x) <= 0 or np.std(y) <= 0):
        return math.nan
    design = np.column_stack((np.ones(len(x)),
                              (covariates - np.mean(covariates, axis=0)) / scales))
    x_residual = np.asarray(x, float) - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_residual = np.asarray(y, float) - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(x_residual) <= 0 or np.std(y_residual) <= 0:
        return math.nan
    return float(pearsonr(x_residual, y_residual).statistic)


def analyze_group(rows, selection, band, replicas, seed):
    n = len(rows)
    result = {"selection": selection, "filter": band, "n_clusters": n}
    if n < 5:
        result.update(status="insufficient_support", spearman_rho=None,
                      spearman_p_nominal=None, partial_pearson_r=None)
        return result
    x = np.asarray([row["icl_delta_200_300"] for row in rows], dtype=float)
    y = np.asarray([row["satellite_delta_50_75"] for row in rows], dtype=float)
    redshift = np.asarray([row["redshift"] for row in rows], dtype=float)
    richness = np.asarray([row["richness"] for row in rows], dtype=float)
    rho = float(spearmanr(x, y).statistic)
    p_value = float(spearmanr(x, y).pvalue)
    partial = partial_pearson(x, y, redshift, richness)
    rng = np.random.default_rng(seed)
    rho_boot = np.full(replicas, np.nan)
    partial_boot = np.full(replicas, np.nan)
    for index in range(replicas):
        draw = rng.integers(0, n, size=n)
        bx, by = x[draw], y[draw]
        if np.std(bx) > 0 and np.std(by) > 0:
            rho_boot[index] = spearmanr(bx, by).statistic
            partial_boot[index] = partial_pearson(bx, by, redshift[draw], richness[draw])

    def interval(values):
        valid = values[np.isfinite(values)]
        return ([float(value) for value in np.percentile(valid, (2.5, 97.5))]
                if len(valid) else [None, None])

    valid_rho = rho_boot[np.isfinite(rho_boot)]
    result.update({
        "status": "exploratory",
        "mean_icl_delta_200_300_njy_arcsec2": float(np.mean(x)),
        "mean_satellite_delta_50_75_njy_arcsec2": float(np.mean(y)),
        "spearman_rho": rho,
        "spearman_p_nominal": p_value,
        "spearman_cluster_bootstrap_ci95": interval(rho_boot),
        "spearman_bootstrap_fraction_positive": (float(np.mean(valid_rho > 0))
                                                   if len(valid_rho) else None),
        "partial_pearson_r_control_z_richness": partial if math.isfinite(partial) else None,
        "partial_cluster_bootstrap_ci95": interval(partial_boot),
        "bootstrap_valid_replicas": int(len(valid_rho)),
    })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--stack", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
                                     "prestack_flux_profiles_2026-10-18_02"))
    parser.add_argument("--satellite", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
                                     "satellite_outer_profile_full_2026-10-18_02"))
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error("replicas must be at least 100")
    data, stack, satellite, work = (args.data.resolve(), args.stack.resolve(),
                                    args.satellite.resolve(), args.work.resolve())
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    stack_summary_path = stack / "summary.json"
    pair_path = stack / "pair_metrics_k2p5_c40_r50_300.tsv"
    satellite_summary_path = satellite / "summary.json"
    satellite_cluster_path = satellite / "cluster_profiles.tsv"
    required = (stack_summary_path, pair_path, satellite_summary_path, satellite_cluster_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        parser.error("Missing input(s): " + ", ".join(missing))
    stack_summary, satellite_summary = read_json(stack_summary_path), read_json(satellite_summary_path)
    if (stack_summary.get("status") != "complete"
            or stack_summary.get("science_gates", {}).get("science_measurement") is not False
            or satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("n_frames_completed") != 800):
        parser.error("Expected complete exploratory inputs with science_measurement=false")
    for path, summary in ((pair_path, stack_summary),
                          (satellite_cluster_path, satellite_summary)):
        if summary.get("output_sha256", {}).get(str(path)) != sha256(path):
            parser.error(f"Input checksum mismatch: {path}")

    pairs = read_tsv(pair_path)
    satellite_rows = read_tsv(satellite_cluster_path)
    satellite_by_key = {}
    for row in satellite_rows:
        if row["radial_bin_kpc"] == "50-75":
            key = (str(row["cluster_id"]), str(row["filter"]))
            if key in satellite_by_key:
                parser.error(f"Duplicate satellite cluster profile: {key}")
            satellite_by_key[key] = row
    pair_keys = [(str(row["id"]), str(row["filter"])) for row in pairs]
    if len(set(pair_keys)) != len(pair_keys):
        parser.error("Pair table contains duplicate cluster/filter rows")

    joined_by_filter = {"i": [], "i2": []}
    for pair in pairs:
        band = str(pair["filter"])
        key = (str(pair["id"]), band)
        profile = satellite_by_key.get(key)
        if profile is None:
            continue
        joined = {
            "cluster_id": key[0], "filter": band,
            "icl_delta_200_300": float(pair["delta_sb_200_300"]),
            "satellite_delta_50_75": float(profile["mean_delta_sb_njy_arcsec2"]),
            "n_satellites": int(profile["n_satellites"]),
            "redshift": float(pair["z"]), "richness": float(pair["richness"]),
            "center_eligible": (pair["auto_unique_center_candidate"] == "True"
                                 and pair["center_masked_k2p5"] == "False"),
        }
        if not all(math.isfinite(joined[key]) for key in
                   ("icl_delta_200_300", "satellite_delta_50_75", "redshift", "richness")):
            parser.error(f"Non-finite joined value for {joined['cluster_id']}/{band}")
        joined_by_filter[band].append(joined)

    results = []
    for band_index, band in enumerate(("i", "i2")):
        all_supported = joined_by_filter[band]
        center_eligible = [row for row in all_supported if row["center_eligible"]]
        results.append(analyze_group(all_supported, "all_supported_pairs", band,
                                     args.replicas, args.seed + band_index * 10000))
        results.append(analyze_group(center_eligible, "center_eligible", band,
                                     args.replicas, args.seed + band_index * 10000 + 1))

    work.mkdir(parents=True, exist_ok=False)
    detail_rows = []
    for band in ("i", "i2"):
        detail_rows.extend(joined_by_filter[band])
    fields = ["cluster_id", "filter", "icl_delta_200_300", "satellite_delta_50_75",
              "n_satellites", "redshift", "richness", "center_eligible"]
    write_tsv(work / "matched_cluster_values.tsv", detail_rows, fields)
    result_fields = sorted({key for row in results for key in row})
    write_tsv(work / "correlation_results.tsv", results, result_fields)

    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), squeeze=False)
    for ax, band in zip(axes[0], ("i", "i2")):
        rows = joined_by_filter[band]
        ordinary = [row for row in rows if not row["center_eligible"]]
        selected = [row for row in rows if row["center_eligible"]]
        ax.scatter([r["icl_delta_200_300"] for r in ordinary],
                   [r["satellite_delta_50_75"] for r in ordinary],
                   s=25, facecolors="none", edgecolors="#89939d", label="center uncertain/masked")
        ax.scatter([r["icl_delta_200_300"] for r in selected],
                   [r["satellite_delta_50_75"] for r in selected],
                   s=30, color="#245a87", label="unique center; center unmasked")
        result = next(r for r in results
                      if r["filter"] == band and r["selection"] == "all_supported_pairs")
        rho = result.get("spearman_rho")
        ax.set_title(f"{band}: n={result['n_clusters']}, Spearman ρ="
                     f"{rho:.2f}" if rho is not None else f"{band}: insufficient overlap")
        ax.set_xlabel("Cluster ICL excess at 200–300 kpc (nJy arcsec⁻²)")
        ax.set_ylabel("Mean satellite excess at 50–75 kpc (nJy arcsec⁻²)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Selected cluster-level association; exploratory, not science-approved")
    fig.tight_layout()
    fig.savefig(work / "satellite_icl_correlation.png", dpi=160)
    plt.close(fig)

    satellite_support_rule = satellite_summary.get("method", {}).get(
        "support", "all seven radial bins required for each satellite profile")
    limitations = [
        "The ICL pair sample has central-mask/BCG limitations; center-eligible sensitivity subsets are smaller.",
        "The retained satellite profiles are an area/mask/control-selected subset of the 800-cluster images.",
        "Neither input has calibrated PSF wings, sky transfer, correlated-noise covariance, extinction, or complete member mass/weight information.",
        "Nominal p-values and this selected-sample correlation do not establish a physical ICL-satellite relation.",
    ]
    if satellite_summary.get("schema") == "hsc-satellite-outer-profile-per-annulus-v1":
        limitations.append(
            "Satellite annuli were selected independently; the 50-75 kpc ICL correlation uses only clusters supported in that outer bin and cannot establish a radial trend on a common cluster sample.")
    summary = {
        "schema": "hsc-satellite-outer-icl-correlation-v1",
        "status": "complete", "science_measurement": False,
        "method": {
            "icl_variable": "raw K2.5, common-coverage >=40%, 50-300 kpc target-control pair delta; 200-300 kpc bin",
            "satellite_variable": "equal-satellite-weighted per-cluster target-control mean; 50-75 kpc bin",
            "satellite_support_rule": satellite_support_rule,
            "matching": "same cluster ID and exact filter; one row per cluster/filter",
            "filters": "i and i2 analyzed independently",
            "covariates": "redshift and richness linear residuals for partial Pearson diagnostic",
            "uncertainty": "ordinary cluster bootstrap; not spatial-block bootstrap and not corrected for coadd covariance",
            "center_eligible": "auto_unique_center_candidate=True and center_masked_k2p5=False",
        },
        "limitations": limitations,
        "replicas": args.replicas, "seed": args.seed,
        "stack_summary_sha256": sha256(stack_summary_path),
        "icl_pair_table_sha256": sha256(pair_path),
        "satellite_summary_sha256": sha256(satellite_summary_path),
        "satellite_cluster_table_sha256": sha256(satellite_cluster_path),
        "results": results,
        "n_matched_clusters": {band: len(joined_by_filter[band]) for band in ("i", "i2")},
        "outputs_sha256": {},
    }
    summary_path = work / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    for path in sorted(work.iterdir()):
        if path.is_file() and path.name != "summary.json":
            summary["outputs_sha256"][path.name] = sha256(path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_matched_clusters": summary["n_matched_clusters"],
                      "results": results, "output": str(work)},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
