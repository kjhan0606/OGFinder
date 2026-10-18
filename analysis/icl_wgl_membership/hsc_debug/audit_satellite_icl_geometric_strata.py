#!/usr/bin/env python3
"""Compare ICL association for satellite outskirts with/without annulus overlap.

The geometric split is fixed by the ICL 200--300 kpc annulus and the satellite
50--75 kpc annulus. This is an exploratory diagnostic, not a flux correction.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_weighted_satellite_icl_association import partial_pearson
from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = ROOT / "weighted_satellite_icl_association_2026-10-18_02"
DEFAULT_PROFILES = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "satellite_icl_geometric_strata_2026-10-18_01"
ICL_INNER_KPC = 200.0
ICL_OUTER_KPC = 300.0
SAT_INNER_KPC = 50.0
SAT_OUTER_KPC = 75.0
NO_OVERLAP_INNER_EDGE = ICL_INNER_KPC - SAT_OUTER_KPC
NO_OVERLAP_OUTER_EDGE = ICL_OUTER_KPC + SAT_OUTER_KPC


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"cannot write an empty table: {path}")
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def geometric_class(cluster_radius_kpc):
    radius = float(cluster_radius_kpc)
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("cluster-centric radius must be finite and non-negative")
    if radius < NO_OVERLAP_INNER_EDGE or radius > NO_OVERLAP_OUTER_EDGE:
        return "no_geometric_overlap"
    return "overlap_possible"


def connected_labels(rows):
    target = [{"ra_deg": float(row["target_ra_deg"]),
               "dec_deg": float(row["target_dec_deg"])} for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]),
                "dec_deg": float(row["control_dec_deg"])} for row in rows]
    indices, labels, *_ = pair_block_labels(target, control, np.ones(len(rows), dtype=bool))
    if len(indices) != len(rows):
        raise ValueError("connected-block helper changed the row count")
    return labels


def partial_correlation(rows, satellite_key):
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[satellite_key]) for row in rows])
    z = np.asarray([float(row["z"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    return partial_pearson(x, y, z, richness)


def block_bootstrap(rows, satellite_key, replicas, rng):
    labels = connected_labels(rows)
    blocks = np.unique(labels)
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five cluster-filter pairs")
    point = partial_correlation(rows, satellite_key)
    if len(blocks) < 6:
        return {"n_blocks": int(len(blocks)), "point_r": point,
                "replicas_requested": replicas, "replicas_valid": 0,
                "ci95": None, "status": "insufficient_spatial_blocks"}
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    draws = []
    invalid = 0
    for _ in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in selected])
        sampled = [rows[index] for index in indices]
        try:
            draws.append(partial_correlation(sampled, satellite_key))
        except ValueError:
            invalid += 1
    return {
        "n_blocks": int(len(blocks)), "point_r": point,
        "replicas_requested": replicas, "replicas_valid": len(draws),
        "replicas_invalid": invalid,
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist() if draws else None,
        "status": "complete_exploratory" if draws else "no_valid_replicates",
    }


def paired_difference_bootstrap(rows, replicas, rng):
    labels = connected_labels(rows)
    blocks = np.unique(labels)
    if len(rows) < 5:
        raise ValueError("paired difference requires at least five common cluster-filter pairs")
    nonoverlap = partial_correlation(rows, "satellite_no_overlap")
    overlap = partial_correlation(rows, "satellite_overlap_possible")
    point = overlap - nonoverlap
    if len(blocks) < 6:
        return {"n_blocks": int(len(blocks)), "point_delta_r": point,
                "r_overlap": overlap, "r_no_overlap": nonoverlap,
                "replicas_requested": replicas, "replicas_valid": 0,
                "ci95": None, "status": "insufficient_spatial_blocks"}
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    draws = []
    invalid = 0
    for _ in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in selected])
        sampled = [rows[index] for index in indices]
        try:
            r_overlap = partial_correlation(sampled, "satellite_overlap_possible")
            r_no_overlap = partial_correlation(sampled, "satellite_no_overlap")
            draws.append(r_overlap - r_no_overlap)
        except ValueError:
            invalid += 1
    return {
        "n_blocks": int(len(blocks)), "point_delta_r": point,
        "r_overlap": overlap, "r_no_overlap": nonoverlap,
        "replicas_requested": replicas, "replicas_valid": len(draws),
        "replicas_invalid": invalid,
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist() if draws else None,
        "status": "complete_exploratory" if draws else "no_valid_replicates",
    }


def load_verified_inputs(association_dir, profile_dir):
    association_dir, profile_dir = Path(association_dir), Path(profile_dir)
    association_summary_path = association_dir / "summary.json"
    association_table = association_dir / "weighted_satellite_icl_pairs.tsv"
    profile_summary_path = profile_dir / "summary.json"
    profile_table = profile_dir / "satellite_profiles.tsv"
    for path in (association_summary_path, association_table,
                 profile_summary_path, profile_table):
        if not path.is_file():
            raise ValueError(f"missing required input: {path}")
    association_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    if (association_summary.get("schema") != "hsc-camira-weighted-satellite-icl-association-v1"
            or association_summary.get("status") != "complete_exploratory"
            or association_summary.get("science_measurement") is not False):
        raise ValueError("unexpected association input or science gate")
    if (association_summary.get("outputs_sha256", {}).get(association_table.name)
            != sha256(association_table)):
        raise ValueError("association table checksum differs from its summary")
    if (profile_summary.get("schema") != "hsc-satellite-outer-profile-per-annulus-v1"
            or profile_summary.get("status") != "complete"
            or profile_summary.get("science_measurement") is not False):
        raise ValueError("unexpected per-annulus profile input or science gate")
    if profile_summary.get("output_sha256", {}).get(str(profile_table)) != sha256(profile_table):
        raise ValueError("satellite profile checksum differs from its summary")
    return association_summary_path, association_table, profile_summary_path, profile_table


def analyze(association_dir, profile_dir, output, replicas=5000, seed=20261020):
    association_summary, association_table, profile_summary, profile_table = load_verified_inputs(
        association_dir, profile_dir)
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    association_rows = read_tsv(association_table)
    satellite_rows = read_tsv(profile_table)
    satellite_groups = {}
    for row in satellite_rows:
        if row["cohort"] != "full_sample" or row["radial_bin_kpc"] != "50-75":
            continue
        key = (str(row["cluster_id"]), row["filter"])
        satellite_groups.setdefault(key, []).append(row)
    strata = {"all": [], "no_geometric_overlap": [], "overlap_possible": []}
    satellite_records = []
    for association in association_rows:
        key = (str(association["cluster_id"]), association["filter"])
        profiles = satellite_groups.get(key, [])
        if not profiles:
            raise ValueError(f"association pair lacks 50--75 kpc satellite data: {key}")
        all_delta = [float(row["delta_sb_njy_arcsec2"]) for row in profiles]
        if not math.isclose(float(association["satellite_unweighted"]),
                            float(np.mean(all_delta)), rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"full satellite mean does not reproduce for {key}")
        selected = {"all": profiles, "no_geometric_overlap": [], "overlap_possible": []}
        for profile in profiles:
            category = geometric_class(profile["cluster_radius_kpc"])
            selected[category].append(profile)
        base = dict(association)
        base["n_satellites_all"] = str(len(profiles))
        for category in ("no_geometric_overlap", "overlap_possible"):
            if selected[category]:
                base[f"satellite_{category}"] = str(float(np.mean([
                    float(row["delta_sb_njy_arcsec2"]) for row in selected[category]])))
                base[f"n_satellites_{category}"] = str(len(selected[category]))
            else:
                base[f"satellite_{category}"] = ""
                base[f"n_satellites_{category}"] = "0"
        satellite_records.append(base)
        for category in ("all", "no_geometric_overlap", "overlap_possible"):
            if not selected[category]:
                continue
            record = dict(association)
            record["satellite_delta"] = str(float(np.mean([
                float(row["delta_sb_njy_arcsec2"]) for row in selected[category]])))
            record["n_satellites_selected"] = str(len(selected[category]))
            record["geometric_stratum"] = category
            strata[category].append(record)

    rng = np.random.default_rng(seed)
    metrics = {}
    for band in ("i", "i2"):
        metrics[band] = {}
        for category, satellite_key in (
                ("all", "satellite_delta"),
                ("no_geometric_overlap", "satellite_delta"),
                ("overlap_possible", "satellite_delta")):
            rows = [row for row in strata[category] if row["filter"] == band]
            if len(rows) < 5:
                metrics[band][category] = {"status": "insufficient_cluster_filter_pairs",
                                           "n_cluster_filter_pairs": len(rows)}
                continue
            metrics[band][category] = {
                "n_cluster_filter_pairs": len(rows),
                "n_satellite_profiles": sum(int(row["n_satellites_selected"]) for row in rows),
                **block_bootstrap(rows, satellite_key, replicas, rng),
            }

        common = []
        no_overlap = {row["cluster_id"]: row for row in strata["no_geometric_overlap"]
                      if row["filter"] == band}
        overlap = {row["cluster_id"]: row for row in strata["overlap_possible"]
                   if row["filter"] == band}
        for cluster_id in sorted(set(no_overlap) & set(overlap)):
            left, right = no_overlap[cluster_id], overlap[cluster_id]
            if left["filter"] != right["filter"]:
                raise ValueError("paired geometric strata have mismatched filters")
            common.append({
                **left,
                "satellite_no_overlap": left["satellite_delta"],
                "satellite_overlap_possible": right["satellite_delta"],
            })
        if len(common) >= 5:
            metrics[band]["paired_common_sample_difference"] = {
                "n_cluster_filter_pairs": len(common),
                **paired_difference_bootstrap(common, replicas, rng),
            }
        else:
            metrics[band]["paired_common_sample_difference"] = {
                "status": "insufficient_common_pairs", "n_cluster_filter_pairs": len(common),
            }

    output.mkdir(parents=True, exist_ok=False)
    rows_path = output / "geometric_strata_pairs.tsv"
    write_tsv(rows_path, satellite_records)
    summary_path = output / "summary.json"
    summary = {
        "schema": "hsc-satellite-icl-geometric-strata-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "geometry": {
            "icl_annulus_kpc": [ICL_INNER_KPC, ICL_OUTER_KPC],
            "satellite_outer_annulus_kpc": [SAT_INNER_KPC, SAT_OUTER_KPC],
            "overlap_possible_cluster_radius_kpc": [NO_OVERLAP_INNER_EDGE,
                                                     NO_OVERLAP_OUTER_EDGE],
            "overlap_rule": "annuli can intersect for 125 <= projected cluster-centric satellite radius <= 375 kpc",
            "no_overlap_rule": "satellite annulus is geometrically disjoint from the ICL annulus outside that interval",
        },
        "method": {
            "association": "partial Pearson correlation after linear residualization on z and richness",
            "uncertainty": "5-degree connected target/control sky-block bootstrap; paired common-sample difference uses identical block draws",
            "replicas": replicas, "seed": seed,
            "profile_aggregation": "equal mean of satellite 50-75 kpc target-minus-azimuthal-control profiles within each cluster/filter and geometric stratum",
            "p_values": "not computed",
        },
        "n_association_pairs": len(association_rows),
        "metrics_by_filter": metrics,
        "input_sha256": {
            "association_summary": sha256(association_summary),
            "association_pairs": sha256(association_table),
            "satellite_profile_summary": sha256(profile_summary),
            "satellite_profiles": sha256(profile_table),
            "analysis_script": sha256(Path(__file__).resolve()),
        },
        "outputs_sha256": {rows_path.name: sha256(rows_path)},
        "limitations": [
            "Geometric non-overlap does not remove PSF wings, sky errors, source-mask leakage, or residual correlated noise.",
            "The association sample and geometric strata are exploratory and were not selected as a preregistered confirmatory analysis.",
            "The block bootstrap does not include per-pixel/coadd covariance, profile measurement uncertainty, or sample-selection correction.",
            "An overlapping-versus-nonoverlapping difference is not a causal decomposition and must not be used to correct the full-sample correlation.",
            "No new HSC data were queried or downloaded; all outputs retain science_measurement=false.",
        ],
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--profiles", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("replicas must be at least 1000")
    result = analyze(args.association, args.profiles, args.output, args.replicas, args.seed)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "pairs": result["n_association_pairs"],
        "metrics_by_filter": result["metrics_by_filter"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
