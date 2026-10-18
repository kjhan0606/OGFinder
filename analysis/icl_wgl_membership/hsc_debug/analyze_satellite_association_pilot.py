#!/usr/bin/env python3
"""Explore CAMIRA member density versus matched-control diffuse-light profiles.

This is deliberately limited to the existing 50-pair validation cohort and
the saved K2.5 / 40% / 50--300 kpc area-support selection. It is a descriptive
pilot, not the satellite outer-light stack or a science measurement.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import pearsonr, spearmanr


EDGES_KPC = np.array([50.0, 100.0, 150.0, 200.0, 300.0])
BIN_NAMES = ["50-100", "100-150", "150-200", "200-300"]
FACTOR_INDEX = 0  # K2.5
SUPPORT_KEY = "k2.5_c40_r50_300__selected_indices"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value):
    value = float(value)
    return value if math.isfinite(value) else None


def association_stats(rows, band):
    selected = [row for row in rows if band == "all" or row["filter"] == band]
    result = {"filter": band, "n_pairs": len(selected), "bins": []}
    for index, name in enumerate(BIN_NAMES):
        x = np.asarray([row["satellite_density_mpc2"][index] for row in selected], float)
        y = np.asarray([row["delta_sb_njy_arcsec2"][index] for row in selected], float)
        item = {
            "radial_bin_kpc": name,
            "mean_satellite_density_mpc2": finite_float(np.mean(x)) if len(x) else None,
            "mean_delta_sb_njy_arcsec2": finite_float(np.mean(y)) if len(y) else None,
            "pearson_r": None,
            "pearson_p_nominal": None,
            "spearman_rho": None,
            "spearman_p_nominal": None,
            "partial_pearson_r_control_z_richness": None,
        }
        if len(x) >= 4 and np.std(x) > 0 and np.std(y) > 0:
            pr = pearsonr(x, y)
            sr = spearmanr(x, y)
            item.update(pearson_r=float(pr.statistic), pearson_p_nominal=float(pr.pvalue),
                        spearman_rho=float(sr.statistic), spearman_p_nominal=float(sr.pvalue))
            if len(x) >= 6:
                covariates = np.asarray([[row["richness"], row["z"]] for row in selected], float)
                scale = np.std(covariates, axis=0)
                if np.all(scale > 0):
                    design = np.column_stack([np.ones(len(x)), (covariates - np.mean(covariates, axis=0)) / scale])
                    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
                    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
                    if np.std(xr) > 0 and np.std(yr) > 0:
                        item["partial_pearson_r_control_z_richness"] = float(pearsonr(xr, yr).statistic)
        result["bins"].append(item)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    data = args.data.resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")

    validation_path = data / "validation_sample_2026-10-13_01/validation_manifest.json"
    sample_path = data / "sample_manifest.json"
    cluster_dir = data / "prestack_clusters_2026-10-15_01"
    cluster_summary_path = cluster_dir / "summary.json"
    cluster_profiles_path = cluster_dir / "profiles.npz"
    extraction_dir = data / "prestack_extraction_2026-10-15_01"
    extraction_summary_path = extraction_dir / "extraction_summary.json"
    mask_dir = data / "prestack_cluster_masks_2026-10-16_01"
    mask_summary_path = mask_dir / "summary.json"
    mask_records_path = mask_dir / "records.tsv"
    pair_dir = data / "prestack_pair_support_2026-10-16_02"
    pair_summary_path = pair_dir / "summary.json"
    support_path = pair_dir / "support_counts.tsv"
    covariance_path = pair_dir / "pilot_covariances.npz"
    input_paths = [validation_path, sample_path, cluster_summary_path, cluster_profiles_path,
                   extraction_summary_path, mask_summary_path, mask_records_path,
                   pair_summary_path, support_path, covariance_path]
    missing = [str(path) for path in input_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input(s): " + ", ".join(missing))

    validation = read_json(validation_path)
    sample = read_json(sample_path)
    cluster_summary = read_json(cluster_summary_path)
    extraction_summary = read_json(extraction_summary_path)
    mask_summary = read_json(mask_summary_path)
    pair_summary = read_json(pair_summary_path)
    if len(validation.get("clusters", [])) != 50 or cluster_summary.get("n_clusters") != 50:
        raise ValueError("Expected the existing 50-cluster validation cohort")
    if cluster_summary.get("status") != "complete" or sha256(cluster_profiles_path) != cluster_summary.get("output_sha256"):
        raise ValueError("Cluster profile archive failed its recorded SHA-256 check")
    if extraction_summary.get("status") != "complete" or mask_summary.get("status") != "complete":
        raise ValueError("An upstream extraction or mask-support summary is not complete")
    if pair_summary.get("status") != "complete" or pair_summary.get("science_measurement") is not False:
        raise ValueError("Unexpected paired-support provenance")
    for output in (support_path, covariance_path, pair_dir / "injection_recovery_sensitivity.tsv"):
        expected = pair_summary.get("output_sha256", {}).get(str(output))
        if expected is None or sha256(output) != expected:
            raise ValueError(f"Paired-support artifact failed its SHA-256 check: {output}")

    sample_rows = sample.get("clusters", [])
    sample_by_id = {str(row["id"]): row for row in sample_rows}
    if len(sample_by_id) != 800:
        raise ValueError(f"Expected 800 unique sample clusters; found {len(sample_by_id)}")
    validation_rows = validation["clusters"]
    if [str(row["id"]) for row in validation_rows] != [str(row["id"]) for row in cluster_summary["records"]]:
        raise ValueError("Cluster profile order differs from the validation manifest")

    mask_records = {str(row["id"]): row for row in read_tsv(mask_records_path)}
    if len(mask_records) != 800:
        raise ValueError(f"Expected 800 mask-support records; found {len(mask_records)}")
    extraction_by_target = {}
    for row in extraction_summary.get("records", []):
        key = str(row["target_cluster_id"])
        if key in extraction_by_target:
            raise ValueError(f"Duplicate matched controls for target {key}")
        extraction_by_target[key] = row

    with np.load(cluster_profiles_path) as archive:
        cluster_raw = archive["raw"].copy()
        cluster_coverage = archive["coverage"].copy()
    if cluster_raw.shape != (50, 3, 6) or cluster_coverage.shape != cluster_raw.shape:
        raise ValueError("Unexpected cluster-profile array shape")
    with np.load(covariance_path) as archive:
        indices = archive[SUPPORT_KEY].astype(int)
    if len(indices) != 24 or len(set(indices.tolist())) != len(indices):
        raise ValueError("Saved exploratory support selection differs from the reviewed 24-pair result")
    cosmology = FlatLambdaCDM(H0=float(validation["h0"]), Om0=float(validation["omega_m"]), Tcmb0=0)
    source_hashes = {str(path): sha256(path) for path in input_paths}
    output_rows = []
    for index in indices:
        frame = validation_rows[int(index)]
        ident = str(frame["id"])
        source = sample_by_id.get(ident)
        mask = mask_records.get(ident)
        control = extraction_by_target.get(ident)
        if source is None or mask is None or control is None:
            raise ValueError(f"Missing member, mask or matched-control record for {ident}")
        control_id = control["control_id"]
        baseline_path = extraction_dir / control_id / "baseline.npz"
        expected_baseline = control.get("baseline_sha256")
        if not baseline_path.is_file() or sha256(baseline_path) != expected_baseline:
            raise ValueError(f"Matched-control profile failed its SHA-256 check: {control_id}")
        source_hashes[str(baseline_path)] = expected_baseline
        with np.load(baseline_path) as archive:
            control_raw = archive["raw"][FACTOR_INDEX, :4].copy()
            control_coverage = archive["coverage"][FACTOR_INDEX, :4].copy()
        target_raw = cluster_raw[int(index), FACTOR_INDEX, :4]
        target_coverage = cluster_coverage[int(index), FACTOR_INDEX, :4]
        if np.any(np.minimum(target_coverage, control_coverage) < 0.4):
            raise ValueError(f"Saved support index fails the K2.5/40% coverage definition: {ident}")
        if mask["primary_filter"] != control["header_filter"]:
            continue

        positions = np.asarray(source.get("members_radec_deg", []), dtype=float)
        if positions.ndim != 2 or positions.shape[1] != 2 or not len(positions):
            raise ValueError(f"Member positions are missing or malformed for {ident}")
        centre = SkyCoord(float(source["ra_deg"]) * u.deg, float(source["dec_deg"]) * u.deg)
        members = SkyCoord(positions[:, 0] * u.deg, positions[:, 1] * u.deg)
        radii_kpc = (centre.separation(members).rad *
                     cosmology.angular_diameter_distance(float(source["z"])).to_value(u.kpc))
        counts = [int(np.count_nonzero((radii_kpc >= lo) & (radii_kpc < hi)))
                  for lo, hi in zip(EDGES_KPC[:-1], EDGES_KPC[1:])]
        area_mpc2 = math.pi * (EDGES_KPC[1:] ** 2 - EDGES_KPC[:-1] ** 2) / 1.0e6
        densities = (np.asarray(counts, float) / area_mpc2).tolist()
        delta_sb = (target_raw - control_raw).tolist()
        common_coverage = np.minimum(target_coverage, control_coverage).tolist()
        centre_masked = ast.literal_eval(mask["center_masked_by_kron"])
        output_rows.append({
            "id": ident,
            "filter": mask["primary_filter"],
            "control_id": control_id,
            "richness": float(source["richness"]),
            "z": float(source["z"]),
            "n_camira_member_positions": int(len(positions)),
            "auto_unique_center_candidate": mask["auto_unique_center_candidate"] == "True",
            "bcg_match_count": int(mask["bcg_match_count"]),
            "center_masked_k2p5": bool(centre_masked["2.5"]),
            "common_coverage": common_coverage,
            "satellite_counts": counts,
            "satellite_density_mpc2": densities,
            "delta_sb_njy_arcsec2": delta_sb,
        })

    if len(output_rows) != 24:
        raise ValueError(f"Expected 24 exact-filter matched pilot pairs; found {len(output_rows)}")
    by_filter = {band: association_stats(output_rows, band) for band in ("i", "i2", "all")}
    profile_summary = {}
    for band in ("i", "i2", "all"):
        selected = [row for row in output_rows if band == "all" or row["filter"] == band]
        profile_summary[band] = {
            "n_pairs": len(selected),
            "mean_delta_sb_njy_arcsec2": [finite_float(np.mean([row["delta_sb_njy_arcsec2"][i] for row in selected]))
                                           for i in range(4)],
            "mean_satellite_density_mpc2": [finite_float(np.mean([row["satellite_density_mpc2"][i] for row in selected]))
                                             for i in range(4)],
            "n_auto_unique_center_candidates": sum(row["auto_unique_center_candidate"] for row in selected),
            "n_center_masked_k2p5": sum(row["center_masked_k2p5"] for row in selected),
        }

    work.mkdir(parents=True, exist_ok=False)
    fields = ["id", "filter", "control_id", "richness", "z", "n_camira_member_positions",
              "auto_unique_center_candidate", "bcg_match_count", "center_masked_k2p5",
              "coverage_50_100", "coverage_100_150", "coverage_150_200", "coverage_200_300",
              "members_50_100", "members_100_150", "members_150_200", "members_200_300",
              "satellite_density_mpc2_50_100", "satellite_density_mpc2_100_150",
              "satellite_density_mpc2_150_200", "satellite_density_mpc2_200_300",
              "delta_sb_njy_arcsec2_50_100", "delta_sb_njy_arcsec2_100_150",
              "delta_sb_njy_arcsec2_150_200", "delta_sb_njy_arcsec2_200_300"]
    flat_rows = []
    for row in output_rows:
        flat = {key: row[key] for key in ["id", "filter", "control_id", "richness", "z",
                                          "n_camira_member_positions", "auto_unique_center_candidate",
                                          "bcg_match_count", "center_masked_k2p5"]}
        for index, name in enumerate(BIN_NAMES):
            suffix = name.replace("-", "_")
            flat[f"coverage_{suffix}"] = row["common_coverage"][index]
            flat[f"members_{suffix}"] = row["satellite_counts"][index]
            flat[f"satellite_density_mpc2_{suffix}"] = row["satellite_density_mpc2"][index]
            flat[f"delta_sb_njy_arcsec2_{suffix}"] = row["delta_sb_njy_arcsec2"][index]
        flat_rows.append(flat)
    write_tsv(work / "cluster_metrics.tsv", flat_rows, fields)
    summary = {
        "schema": "hsc-icl-camira-association-pilot-v1",
        "status": "complete",
        "support_definition": "saved K2.5, 40% common cluster/control area, 50-300 kpc support from the 50-pair validation cohort",
        "n_validation_pairs_supported": 24,
        "n_exact_filter_pairs": len(output_rows),
        "filter_counts": {band: sum(row["filter"] == band for row in output_rows) for band in ("i", "i2")},
        "radial_bins_kpc": BIN_NAMES,
        "profile_summary": profile_summary,
        "association_tests": by_filter,
        "statistical_caveats": [
            "Correlation p-values are nominal, descriptive, and uncorrected for four radial-bin comparisons.",
            "No spatial-block resampling was applied to the member-light correlation.",
            "The 24 pairs are a coverage-selected subset of a 50-pair engineering cohort, not the 800-cluster science sample.",
            "CAMIRA member positions are counts only here: no member photometry, stellar mass, membership-probability weighting, or contamination subtraction is used.",
            "The light profile is raw cluster minus matched-control native annular SB; Galactic extinction, i/i2 cross-band calibration, PSF-wing transfer, and upstream sky transfer remain unresolved.",
            "Only 8 of these 24 targets are automatic unique-center candidates; center/BCG confirmation is incomplete.",
            "This does not measure satellite-centered outer-light profiles and cannot establish the requested ICL--satellite outer-structure relation.",
        ],
        "science_measurement": False,
        "weak_lensing_measurement": False,
        "source_sha256": source_hashes,
        "script_sha256": {str(Path(__file__).resolve()): sha256(__file__)},
    }
    with (work / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": summary["status"], "output": str(work),
                      "n_pairs": len(output_rows), "filters": summary["filter_counts"],
                      "science_measurement": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
