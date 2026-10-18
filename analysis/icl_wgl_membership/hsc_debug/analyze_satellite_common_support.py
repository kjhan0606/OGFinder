#!/usr/bin/env python3
"""Recompute satellite profiles on clusters supported in every radial bin.

This is an exploratory common-cluster sensitivity analysis. It uses only the
existing satellite cluster table, sample manifest, and annulus-supported ICL
pair table; it does not read images or request new data.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr


BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")


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


def sky_cells(ra, dec, width=5.0):
    if not 0 < width <= 30:
        raise ValueError("block width must be in (0, 30] degrees")
    return np.asarray([f"{int(float(r) // width)}_{int((float(d) + 90) // width)}"
                       for r, d in zip(ra, dec)], dtype=object)


def block_bootstrap_mean(values, labels, replicas, seed):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=object)
    unique = np.unique(labels)
    if values.ndim != 2 or len(values) != len(labels) or len(unique) < 2:
        raise ValueError("Need aligned profile rows in at least two spatial blocks")
    members = {label: values[labels == label] for label in unique}
    rng = np.random.default_rng(seed)
    draws = np.empty((replicas, values.shape[1]), dtype=float)
    for index in range(replicas):
        selected = rng.choice(unique, size=len(unique), replace=True)
        sample = np.concatenate([members[label] for label in selected], axis=0)
        draws[index] = sample.mean(axis=0)
    covariance = np.atleast_2d(np.cov(draws, rowvar=False, ddof=1))
    low, high = np.percentile(draws, (2.5, 97.5), axis=0)
    return draws.mean(axis=0), low, high, covariance, len(unique), max(
        len(members[label]) for label in unique)


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


def association_result(rows, selection, band, replicas, seed):
    n = len(rows)
    base = {"selection": selection, "filter": band, "n_clusters": n}
    if n < 5:
        return {**base, "status": "insufficient_support", "spearman_rho": None,
                "spearman_p_nominal": None, "partial_pearson_r": None,
                "spearman_block_bootstrap_ci95": [None, None],
                "partial_block_bootstrap_ci95": [None, None]}
    x = np.asarray([row["icl_delta_200_300"] for row in rows], dtype=float)
    y = np.asarray([row["satellite_delta_50_75"] for row in rows], dtype=float)
    z = np.asarray([row["redshift"] for row in rows], dtype=float)
    richness = np.asarray([row["richness"] for row in rows], dtype=float)
    labels = np.asarray([row["sky_cell"] for row in rows], dtype=object)
    unique = np.unique(labels)
    if len(unique) < 2:
        return {**base, "status": "insufficient_spatial_blocks", "n_spatial_blocks": len(unique),
                "spearman_rho": None, "spearman_p_nominal": None,
                "partial_pearson_r": None,
                "spearman_block_bootstrap_ci95": [None, None],
                "partial_block_bootstrap_ci95": [None, None]}

    rho, p_value = spearmanr(x, y)
    partial = partial_pearson(x, y, z, richness)
    members = {label: np.flatnonzero(labels == label) for label in unique}
    rng = np.random.default_rng(seed)
    rho_draws, partial_draws = [], []
    for _ in range(replicas):
        selected = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([members[label] for label in selected])
        bx, by = x[indices], y[indices]
        if np.std(bx) > 0 and np.std(by) > 0:
            rho_draws.append(float(spearmanr(bx, by).statistic))
            value = partial_pearson(bx, by, z[indices], richness[indices])
            if math.isfinite(value):
                partial_draws.append(value)

    def interval(values):
        return ([float(value) for value in np.percentile(values, (2.5, 97.5))]
                if len(values) else [None, None])

    return {
        **base, "status": "exploratory", "n_spatial_blocks": len(unique),
        "mean_icl_delta_200_300_njy_arcsec2": float(x.mean()),
        "mean_satellite_delta_50_75_njy_arcsec2": float(y.mean()),
        "spearman_rho": float(rho), "spearman_p_nominal": float(p_value),
        "spearman_block_bootstrap_ci95": interval(rho_draws),
        "partial_pearson_r": partial if math.isfinite(partial) else None,
        "partial_block_bootstrap_ci95": interval(partial_draws),
        "valid_spearman_replicas": len(rho_draws),
        "valid_partial_replicas": len(partial_draws),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite", type=Path, required=True)
    parser.add_argument("--icl", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--block-width-deg", type=float, default=5.0)
    args = parser.parse_args()
    if args.replicas < 100 or not 0 < args.block_width_deg <= 30:
        parser.error("Invalid bootstrap or block-width configuration")

    satellite, icl, work = (args.satellite.resolve(), args.icl.resolve(),
                            args.work.resolve())
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    sat_summary_path = satellite / "summary.json"
    sat_profiles_path = satellite / "cluster_profiles.tsv"
    icl_summary_path = icl / "summary.json"
    icl_pairs_path = icl / "pair_metrics_k2p5_c40_r50_300.tsv"
    sat_summary, icl_summary = read_json(sat_summary_path), read_json(icl_summary_path)
    if (sat_summary.get("schema") != "hsc-satellite-outer-profile-per-annulus-v1"
            or sat_summary.get("status") != "complete"
            or sat_summary.get("science_measurement") is not False
            or sat_summary.get("n_frames_completed") != 800
            or icl_summary.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_summary.get("status") != "complete"
            or icl_summary.get("science_gates", {}).get("science_measurement") is not False):
        parser.error("Expected complete exploratory 800-frame inputs")
    if (sat_summary.get("output_sha256", {}).get(str(sat_profiles_path)) != sha256(sat_profiles_path)
            or icl_summary.get("output_sha256", {}).get(str(icl_pairs_path)) != sha256(icl_pairs_path)):
        parser.error("Satellite or ICL input checksum mismatch")

    manifest_path = Path(sat_summary["input_manifest"])
    if sha256(manifest_path) != sat_summary.get("manifest_sha256"):
        parser.error("Sample-manifest checksum mismatch")
    manifest = read_json(manifest_path)
    frames = manifest.get("clusters", [])
    if len(frames) != 800:
        parser.error("Expected an 800-cluster sample manifest")
    frame_by_id = {str(row["id"]): row for row in frames}
    if len(frame_by_id) != 800:
        parser.error("Sample manifest contains duplicate cluster IDs")

    satellite_rows = read_tsv(sat_profiles_path)
    grouped = {band: {radial_bin: {} for radial_bin in BINS} for band in ("i", "i2")}
    for row in satellite_rows:
        band, radial_bin, ident = row["filter"], row["radial_bin_kpc"], str(row["cluster_id"])
        if band not in grouped or radial_bin not in BINS:
            parser.error(f"Unexpected satellite profile key: {band}/{radial_bin}")
        records = grouped[band][radial_bin]
        if ident in records:
            parser.error(f"Duplicate satellite cluster profile: {ident}/{band}/{radial_bin}")
        frame = frame_by_id.get(ident)
        if frame is None or frame.get("coadd_header_filter") != band:
            parser.error(f"Satellite/manifest ID or filter mismatch: {ident}/{band}")
        if ("redshift" in row and not math.isclose(
                float(row["redshift"]), float(frame["z"]), rel_tol=0, abs_tol=1e-10)):
            parser.error(f"Satellite/manifest redshift mismatch: {ident}/{band}")
        if ("richness" in row and not math.isclose(
                float(row["richness"]), float(frame["richness"]), rel_tol=0, abs_tol=1e-8)):
            parser.error(f"Satellite/manifest richness mismatch: {ident}/{band}")
        records[ident] = row

    profile_rows, covariance_rows, common_rows = [], [], []
    covariance_archive = {}
    common_by_band = {}
    for band_index, band in enumerate(("i", "i2")):
        common = set.intersection(*(set(grouped[band][radial_bin]) for radial_bin in BINS))
        if len(common) < 5:
            parser.error(f"Insufficient all-annulus cluster support in {band}: {len(common)}")
        ids = sorted(common)
        rows_by_bin = [[grouped[band][radial_bin][ident] for ident in ids]
                       for radial_bin in BINS]
        values = np.column_stack([[float(row["mean_delta_sb_njy_arcsec2"])
                                   for row in rows] for rows in rows_by_bin])
        ra = [frame_by_id[ident]["ra_deg"] for ident in ids]
        dec = [frame_by_id[ident]["dec_deg"] for ident in ids]
        labels = sky_cells(ra, dec, args.block_width_deg)
        boot_mean, low, high, covariance, n_blocks, largest = block_bootstrap_mean(
            values, labels, args.replicas, args.seed + band_index)
        covariance_archive[band] = covariance
        common_by_band[band] = set(ids)
        for index, radial_bin in enumerate(BINS):
            n_satellites = sum(int(row["n_satellites"]) for row in rows_by_bin[index])
            profile_rows.append({
                "filter": band, "radial_bin_kpc": radial_bin,
                "n_common_clusters": len(ids), "n_satellites": n_satellites,
                "mean_equal_cluster_delta_sb_njy_arcsec2": float(values[:, index].mean()),
                "spatial_block_bootstrap_mean_njy_arcsec2": float(boot_mean[index]),
                "spatial_block_bootstrap_95_low": float(low[index]),
                "spatial_block_bootstrap_95_high": float(high[index]),
                "n_spatial_blocks": n_blocks, "largest_block_clusters": largest,
            })
            for other_index, other_bin in enumerate(BINS):
                covariance_rows.append({
                    "filter": band, "radial_bin_kpc": radial_bin,
                    "other_radial_bin_kpc": other_bin,
                    "bootstrap_covariance_njy2_arcsec4": float(covariance[index, other_index]),
                })
        for ident in ids:
            frame = frame_by_id[ident]
            common_rows.append({
                "filter": band, "cluster_id": ident, "ra_deg": float(frame["ra_deg"]),
                "dec_deg": float(frame["dec_deg"]), "sky_cell": str(labels[ids.index(ident)]),
                "redshift": float(frame["z"]), "richness": float(frame["richness"]),
            })

    icl_rows = read_tsv(icl_pairs_path)
    icl_by_key = {}
    for row in icl_rows:
        key = (str(row["id"]), row["filter"])
        if key in icl_by_key:
            parser.error(f"Duplicate ICL pair row: {key}")
        icl_by_key[key] = row

    association_rows = []
    results = []
    for band_index, band in enumerate(("i", "i2")):
        band_rows = []
        for ident in sorted(common_by_band[band]):
            pair = icl_by_key.get((ident, band))
            if pair is None:
                continue
            sat = grouped[band]["50-75"][ident]
            frame = frame_by_id[ident]
            redshift, richness = float(frame["z"]), float(frame["richness"])
            if (not math.isclose(float(pair["z"]), redshift, rel_tol=0, abs_tol=1e-10)
                    or not math.isclose(float(pair["richness"]), richness, rel_tol=0, abs_tol=1e-8)):
                parser.error(f"ICL/satellite z or richness mismatch: {ident}/{band}")
            row = {
                "cluster_id": ident, "filter": band,
                "icl_delta_200_300": float(pair["delta_sb_200_300"]),
                "satellite_delta_50_75": float(sat["mean_delta_sb_njy_arcsec2"]),
                "redshift": redshift, "richness": richness,
                "sky_cell": str(sky_cells([frame["ra_deg"]], [frame["dec_deg"]],
                                          args.block_width_deg)[0]),
                "center_eligible": (pair["auto_unique_center_candidate"] == "True"
                                    and pair["center_masked_k2p5"] == "False"),
            }
            band_rows.append(row)
            association_rows.append({**row, "sky_cell": row["sky_cell"]})
        results.append(association_result(
            band_rows, "all_common_annulus_clusters", band, args.replicas,
            args.seed + 10000 + band_index * 100))
        eligible = [row for row in band_rows if row["center_eligible"]]
        results.append(association_result(
            eligible, "center_eligible", band, args.replicas,
            args.seed + 10001 + band_index * 100))

    work.mkdir(parents=True, exist_ok=False)
    write_tsv(work / "common_support_profile.tsv", profile_rows,
              ["filter", "radial_bin_kpc", "n_common_clusters", "n_satellites",
               "mean_equal_cluster_delta_sb_njy_arcsec2",
               "spatial_block_bootstrap_mean_njy_arcsec2", "spatial_block_bootstrap_95_low",
               "spatial_block_bootstrap_95_high", "n_spatial_blocks", "largest_block_clusters"])
    write_tsv(work / "profile_covariance.tsv", covariance_rows,
              ["filter", "radial_bin_kpc", "other_radial_bin_kpc",
               "bootstrap_covariance_njy2_arcsec4"])
    write_tsv(work / "common_cluster_ids.tsv", common_rows,
              ["filter", "cluster_id", "ra_deg", "dec_deg", "sky_cell", "redshift", "richness"])
    write_tsv(work / "icl_satellite_common_support.tsv", association_rows,
              ["cluster_id", "filter", "icl_delta_200_300", "satellite_delta_50_75",
               "redshift", "richness", "sky_cell", "center_eligible"])
    result_fields = sorted({key for row in results for key in row})
    write_tsv(work / "association_results.tsv", results, result_fields)
    np.savez_compressed(work / "profile_bootstrap_covariance.npz", **covariance_archive)

    outputs = sorted(path for path in work.iterdir() if path.is_file())
    summary = {
        "schema": "hsc-satellite-common-cluster-support-v1", "status": "complete",
        "science_measurement": False,
        "method": {
            "sample": "intersection of clusters with at least one retained satellite annulus profile in all seven radial bins, separately by filter",
            "weighting": "equal cluster weights; satellite counts may differ by cluster and radial bin",
            "uncertainty": f"{args.block_width_deg:g}-degree target-sky-cell block bootstrap; covariance preserves all seven bins jointly",
            "icl_link": "same cluster ID and filter; per-annulus ICL 200-300 kpc pair; common satellite support only",
            "selection": "no CAMIRA member mass or weight; ICL center eligibility assessed as a sensitivity subset",
        },
        "n_satellite_cluster_bin_rows": len(satellite_rows),
        "n_common_clusters": {band: len(common_by_band[band]) for band in ("i", "i2")},
        "n_spatial_blocks": {band: int(next(row["n_spatial_blocks"] for row in profile_rows
                                             if row["filter"] == band)) for band in ("i", "i2")},
        "n_matched_icl_clusters": {band: sum(row["filter"] == band for row in association_rows)
                                   for band in ("i", "i2")},
        "replicas": args.replicas, "seed": args.seed,
        "block_width_deg": args.block_width_deg,
        "input_sha256": {
            str(sat_summary_path): sha256(sat_summary_path),
            str(sat_profiles_path): sha256(sat_profiles_path),
            str(manifest_path): sha256(manifest_path),
            str(icl_summary_path): sha256(icl_summary_path),
            str(icl_pairs_path): sha256(icl_pairs_path),
        },
        "program_sha256": sha256(Path(__file__).resolve()),
        "results": results,
        "output_sha256": {str(path): sha256(path) for path in outputs},
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                    allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_common_clusters": summary["n_common_clusters"],
                      "n_matched_icl_clusters": summary["n_matched_icl_clusters"],
                      "results": results, "output": str(work)},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
