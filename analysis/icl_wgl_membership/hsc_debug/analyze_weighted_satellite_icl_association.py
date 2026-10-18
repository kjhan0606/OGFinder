#!/usr/bin/env python3
"""Exploratory CAMIRA-weighted satellite-light versus ICL association.

This pairs the existing 50--75 kpc satellite-centered profile with the
200--300 kpc target-minus-control ICL profile for the fixed K2.5/40% support
branch. It compares equal-satellite and CAMIRA-weighted cluster summaries and
uses the same connected 5-degree sky blocks for a paired bootstrap. The result
is a sensitivity analysis, not a science measurement or a confirmatory test.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ICL = DEFAULT_ROOT / "icl_satellite_cut_sensitivity_2026-10-18_01"
DEFAULT_MEMBER = DEFAULT_ROOT / "satellite_member_weight_sensitivity_2026-10-18_04"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def cluster_satellite_summaries(member_rows):
    """Collapse satellite annuli to per-cluster unweighted and w-weighted means."""
    grouped = {}
    seen = set()
    for row in member_rows:
        if row["radial_bin_kpc"] != "50-75":
            continue
        cluster_id, band = str(row["cluster_id"]), row["filter"]
        key = (cluster_id, band)
        if "member_index" not in row:
            raise ValueError("satellite profile table must preserve member_index")
        satellite_key = (*key, str(row["member_index"]))
        if satellite_key in seen:
            raise ValueError(f"duplicate satellite 50--75 kpc profile row: {satellite_key}")
        seen.add(satellite_key)
        delta, weight = float(row["delta"]), float(row["camira_weight"])
        if not math.isfinite(delta) or not math.isfinite(weight) or not 0 <= weight <= 1:
            raise ValueError(f"invalid satellite delta or CAMIRA weight: {key}")
        grouped.setdefault(key, []).append((delta, weight))

    output = {}
    for key, values in grouped.items():
        delta = np.asarray([v[0] for v in values], dtype=float)
        weight = np.asarray([v[1] for v in values], dtype=float)
        if weight.sum() <= 0:
            continue
        output[key] = {
            "satellite_unweighted": float(delta.mean()),
            "satellite_weighted": float(np.average(delta, weights=weight)),
            "n_satellites": int(len(delta)),
            "n_eff_weight": float(weight.sum() ** 2 / np.sum(weight ** 2)),
            "weight_sum": float(weight.sum()),
        }
    return output


def partial_pearson(x, y, redshift, richness):
    """Pearson correlation after linear residualization on z and richness."""
    arrays = [np.asarray(value, dtype=float) for value in (x, y, redshift, richness)]
    if any(value.ndim != 1 for value in arrays) or len({len(value) for value in arrays}) != 1:
        raise ValueError("partial-correlation inputs must be equal-length vectors")
    if len(arrays[0]) < 5 or any(not np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("partial-correlation inputs need >=5 finite rows")
    design = np.column_stack((np.ones(len(arrays[0])), arrays[2], arrays[3]))
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("redshift/richness nuisance design is rank deficient")
    x_resid = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    y_resid = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    denominator = float(np.linalg.norm(x_resid) * np.linalg.norm(y_resid))
    if denominator <= 0:
        raise ValueError("partial correlation is undefined for a constant residual")
    return float(x_resid @ y_resid / denominator)


def block_bootstrap_correlations(rows, labels, replicas=10000, seed=20261018):
    """Paired block bootstrap for weighted/unweighted partial correlations."""
    if len(rows) != len(labels):
        raise ValueError("one spatial block label is required per row")
    labels = np.asarray(labels, dtype=str)
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("at least six connected sky blocks are required")
    rng = np.random.default_rng(seed)
    unweighted = np.asarray([float(r["satellite_unweighted"]) for r in rows])
    weighted = np.asarray([float(r["satellite_weighted"]) for r in rows])
    icl = np.asarray([float(r["icl_delta_200_300"]) for r in rows])
    z = np.asarray([float(r["z"]) for r in rows])
    richness = np.asarray([float(r["richness"]) for r in rows])
    point_unweighted = partial_pearson(icl, unweighted, z, richness)
    point_weighted = partial_pearson(icl, weighted, z, richness)
    draws = np.empty((replicas, 3), dtype=float)
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    for replica in range(replicas):
        selected_blocks = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in selected_blocks])
        r_unweighted = partial_pearson(icl[indices], unweighted[indices], z[indices], richness[indices])
        r_weighted = partial_pearson(icl[indices], weighted[indices], z[indices], richness[indices])
        draws[replica] = (r_unweighted, r_weighted, r_weighted - r_unweighted)
    quantiles = np.quantile(draws, [0.025, 0.5, 0.975], axis=0)
    return {
        "n_blocks": int(len(blocks)),
        "replicas": int(replicas),
        "point_r_unweighted": point_unweighted,
        "point_r_weighted": point_weighted,
        "point_delta_r_weighted_minus_unweighted": point_weighted - point_unweighted,
        "bootstrap_ci95_r_unweighted": quantiles[[0, 2], 0].tolist(),
        "bootstrap_ci95_r_weighted": quantiles[[0, 2], 1].tolist(),
        "bootstrap_ci95_delta_r": quantiles[[0, 2], 2].tolist(),
        "bootstrap_median_delta_r": float(quantiles[1, 2]),
        "science_measurement": False,
    }


def load_inputs(icl_dir, member_dir):
    icl_dir, member_dir = Path(icl_dir), Path(member_dir)
    icl_summary_path = icl_dir / "summary.json"
    icl_table = icl_dir / "matched_cluster_values.tsv"
    member_summary_path = member_dir / "summary.json"
    member_table = member_dir / "profile_member_properties.tsv"
    for path in (icl_summary_path, icl_table, member_summary_path, member_table):
        if not path.is_file():
            raise ValueError(f"missing required input: {path}")
    icl_summary = json.loads(icl_summary_path.read_text(encoding="utf-8"))
    member_summary = json.loads(member_summary_path.read_text(encoding="utf-8"))
    if icl_summary.get("schema") != "hsc-icl-satellite-cut-sensitivity-v1" \
            or icl_summary.get("status") != "complete" \
            or icl_summary.get("science_measurement") is not False:
        raise ValueError("unexpected ICL sensitivity input or science gate")
    if icl_summary.get("output_sha256", {}).get(str(icl_table)) != sha256(icl_table):
        raise ValueError("ICL matched-cluster table hash mismatch")
    if member_summary.get("status") != "complete_exploratory" \
            or member_summary.get("science_measurement") is not False:
        raise ValueError("unexpected member-weighted profile input or science gate")
    if member_summary.get("column_mapping", {}).get("column_8", "").startswith("weight,") is False:
        raise ValueError("member weight column is not identified by acquisition metadata")
    if member_summary.get("outputs_sha256", {}).get(member_table.name) != sha256(member_table):
        raise ValueError("member profile table hash mismatch")
    return (icl_summary_path, icl_table, member_summary_path, member_table,
            icl_summary, member_summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--members", type=Path, default=DEFAULT_MEMBER)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("replicas must be at least 1000")
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    icl_summary_path, icl_table_path, member_summary_path, member_table_path, icl_summary, member_summary = load_inputs(args.icl.resolve(), args.members.resolve())
    cluster_satellites = cluster_satellite_summaries(read_tsv(member_table_path))
    raw_rows = read_tsv(icl_table_path)
    rows = []
    seen = set()
    for row in raw_rows:
        if float(row["kron_factor"]) != 2.5 or float(row["coverage_floor"]) != 0.40:
            continue
        key = (str(row["cluster_id"]), row["filter"])
        if key in seen:
            raise ValueError(f"duplicate fixed-cut ICL cluster row: {key}")
        seen.add(key)
        satellite = cluster_satellites.get(key)
        if satellite is None:
            continue
        equal_reference = float(row["satellite_delta_50_75"])
        if not math.isclose(satellite["satellite_unweighted"], equal_reference,
                            rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"equal-satellite profile does not reproduce ICL table value: {key}")
        rows.append({
            "cluster_id": key[0], "filter": key[1],
            "icl_delta_200_300": float(row["icl_delta_200_300"]),
            "satellite_unweighted": satellite["satellite_unweighted"],
            "satellite_weighted": satellite["satellite_weighted"],
            "n_satellites": satellite["n_satellites"],
            "n_eff_weight": satellite["n_eff_weight"],
            "weight_sum": satellite["weight_sum"],
            "z": float(row["z"]), "richness": float(row["richness"]),
            "target_ra_deg": float(row["target_ra_deg"]),
            "target_dec_deg": float(row["target_dec_deg"]),
            "control_ra_deg": float(row["control_ra_deg"]),
            "control_dec_deg": float(row["control_dec_deg"]),
            "center_eligible": row["center_eligible"],
        })
    if len({(row["cluster_id"], row["filter"]) for row in rows}) != len(rows):
        raise ValueError("duplicate joined cluster/filter rows")
    if not rows:
        parser.error("no cluster/filter rows overlap the fixed ICL and satellite samples")
    from evaluate_prestack_pair_support import pair_block_labels
    target = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in rows]
    control = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in rows]
    indices, labels, *_ = pair_block_labels(target, control, np.ones(len(rows), dtype=bool))
    if len(indices) != len(rows):
        raise ValueError("spatial block builder changed joined row order")
    work.mkdir(parents=True, exist_ok=False)
    outputs = []
    stats = {}
    for band in ("i", "i2"):
        chosen = [row for row in rows if row["filter"] == band]
        chosen_indices = [i for i, row in enumerate(rows) if row["filter"] == band]
        if len(chosen) < 5:
            stats[band] = {"status": "insufficient_clusters", "n_clusters": len(chosen)}
            continue
        result = block_bootstrap_correlations(
            chosen, [labels[i] for i in chosen_indices], args.replicas,
            args.seed + (0 if band == "i" else 1),
        )
        result.update(status="exploratory_sensitivity", filter=band,
                      n_clusters=len(chosen), n_satellite_profile_sum=sum(r["n_satellites"] for r in chosen),
                      sum_cluster_effective_satellite_counts=float(
                          sum(r["n_eff_weight"] for r in chosen)))
        stats[band] = result
    table_path = work / "weighted_satellite_icl_pairs.tsv"
    fields = ["cluster_id", "filter", "icl_delta_200_300", "satellite_unweighted",
              "satellite_weighted", "n_satellites", "n_eff_weight", "weight_sum",
              "z", "richness", "target_ra_deg", "target_dec_deg",
              "control_ra_deg", "control_dec_deg", "center_eligible"]
    write_tsv(table_path, rows, fields)
    summary_path = work / "summary.json"
    summary = {
        "schema": "hsc-camira-weighted-satellite-icl-association-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": "per-cluster CAMIRA-weighted 50-75 kpc satellite-centered profile versus K2.5/40% 200-300 kpc ICL profile; z/richness adjusted partial correlation with paired 5-degree connected-block bootstrap",
        "selection": {"kron_factor": 2.5, "coverage_floor": 0.40,
                      "icl_annulus_kpc": [200, 300], "satellite_annulus_kpc": [50, 75]},
        "n_joined_cluster_filter_rows": len(rows), "results_by_filter": stats,
        "limitations": [
            "This is a post hoc sensitivity analysis, not a confirmatory or familywise-corrected test.",
            "CAMIRA weights include the algorithm's mass and spatial filters and do not give an unbiased satellite population average.",
            "Column 7 mass interpretation is unresolved and is not used here.",
            "ICL and satellite profiles remain exploratory; masks, sky/PSF transfer, completeness, and correlated coadd errors are not calibrated.",
            "Bootstrap intervals account for connected 5-degree sky blocks only and do not propagate per-satellite measurement error.",
        ],
        "input_sha256": {str(path): sha256(path) for path in
                          (icl_summary_path, icl_table_path, member_summary_path, member_table_path)},
        "member_acquisition_summary_sha256": member_summary.get("acquisition_schema", {}).get("source_summary_sha256"),
        "analysis_script_sha256": sha256(Path(__file__)),
        "bootstrap_replicas": args.replicas, "seed": args.seed,
        "outputs_sha256": {},
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    summary["outputs_sha256"] = {table_path.name: sha256(table_path)}
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "joined_rows": len(rows), "results_by_filter": stats}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
