#!/usr/bin/env python3
"""Measure support and association sensitivity to satellite coverage floor.

Compares the existing 60% per-annulus profile to a local CPU rerun at 40%.
Both products are exploratory. The relaxed-floor result is a sensitivity, not
an adopted science estimator.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_icl_target_control_outer_sky import partial_pearson
from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ICL_PAIRS = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
ICL_PROFILES = DATA / "icl_annulus_support_2026-10-18_01/annulus_pair_profiles.tsv"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
SATELLITE_60 = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01/cluster_profiles.tsv"
SATELLITE_60_SUMMARY = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01/summary.json"
SATELLITE_40 = DATA / "satellite_outer_profile_perbin_cov040_2026-11-05_01/cluster_profiles.tsv"
SATELLITE_40_SUMMARY = DATA / "satellite_outer_profile_perbin_cov040_2026-11-05_01/summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_data() -> tuple[dict[str, dict], dict[str, dict]]:
    for path in (ICL_PAIRS, ICL_PROFILES, CONTROL_MANIFEST, SATELLITE_60,
                 SATELLITE_60_SUMMARY, SATELLITE_40, SATELLITE_40_SUMMARY):
        if not path.is_file():
            raise FileNotFoundError(path)
    s60 = json.loads(SATELLITE_60_SUMMARY.read_text(encoding="utf-8"))
    s40 = json.loads(SATELLITE_40_SUMMARY.read_text(encoding="utf-8"))
    if (s60.get("status") != "complete" or s40.get("status") != "complete"
            or s60.get("science_measurement") is not False
            or s40.get("science_measurement") is not False
            or s60.get("method", {}).get("coverage_floor") != 0.6
            or s40.get("method", {}).get("coverage_floor") != 0.4
            or s60.get("n_frames_completed") != 800 or s40.get("n_frames_completed") != 800):
        raise ValueError("expected matching full-sample 60% and 40% exploratory profiles")

    icl_metric = {}
    for row in read_tsv(ICL_PAIRS):
        key = (row["id"], row["filter"])
        if key in icl_metric:
            raise ValueError(f"duplicate ICL support key: {key}")
        icl_metric[key] = {
            "delta": float(row["delta_sb_200_300"]),
            "center_eligible": row["center_masked_k2p5"] == "False",
        }
    icl_meta = {}
    for row in read_tsv(ICL_PROFILES):
        if row["radial_bin_kpc"] != "200-300":
            continue
        key = (row["id"], row["filter"])
        if key in icl_meta:
            raise ValueError(f"duplicate ICL coordinate/support key: {key}")
        icl_meta[key] = {
            "z": float(row["z"]), "richness": float(row["richness"]),
            "target_ra_deg": float(row["ra_deg"]), "target_dec_deg": float(row["dec_deg"]),
        }
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    controls = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    if len(controls) != 798:
        raise ValueError("matched-control manifest has an unexpected number of rows")

    def load_satellite(path: Path) -> dict[tuple[str, str], dict]:
        selected = [r for r in read_tsv(path) if r["radial_bin_kpc"] == "50-75"]
        rows = {}
        for row in selected:
            key = (row["cluster_id"], row["filter"])
            if key in rows:
                raise ValueError(f"duplicate 50-75 kpc satellite profile: {key}")
            rows[key] = row
        return rows

    sat60, sat40 = load_satellite(SATELLITE_60), load_satellite(SATELLITE_40)
    if not set(sat60).issubset(sat40):
        raise ValueError("relaxed 40% coverage support is not a superset of 60% support")

    def samples(satellite: dict[tuple[str, str], dict]) -> dict[str, list[dict]]:
        result = {"i": [], "i2": []}
        for (ident, band), sat in satellite.items():
            key = (ident, band)
            if key not in icl_metric or key not in icl_meta or ident not in controls:
                continue
            ctrl = controls[ident]
            result[band].append({
                "cluster_id": ident, "filter": band,
                "icl_delta": icl_metric[key]["delta"],
                "satellite_delta": float(sat["mean_delta_sb_njy_arcsec2"]),
                "n_satellites": int(sat["n_satellites"]),
                "center_eligible": icl_metric[key]["center_eligible"],
                **icl_meta[key],
                "control_ra_deg": float(ctrl["ra_deg"]),
                "control_dec_deg": float(ctrl["dec_deg"]),
            })
        for band in result:
            result[band].sort(key=lambda r: int(r["cluster_id"]))
        return result

    return samples(sat60), samples(sat40)


def labels_for(rows: list[dict]) -> np.ndarray:
    target = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    control = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, _, _, _, _ = pair_block_labels(target, control, np.ones(len(rows), dtype=bool))
    return np.asarray(labels, dtype=str)


def partial(rows: list[dict]) -> float:
    return partial_pearson(
        [r["icl_delta"] for r in rows], [r["satellite_delta"] for r in rows],
        [[r["z"] for r in rows], [r["richness"] for r in rows]],
    )


def ci(values: list[float] | np.ndarray) -> list[float]:
    return [float(x) for x in np.percentile(np.asarray(values, dtype=float), [2.5, 97.5])]


def bootstrap_partial(rows: list[dict], labels: np.ndarray, replicas: int, seed: int) -> list[float]:
    groups = np.unique(labels)
    by_block = {block: np.flatnonzero(labels == block) for block in groups}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_block[block] for block in chosen])
        value = partial([rows[i] for i in indices])
        if np.isfinite(value):
            values.append(float(value))
    if len(values) < 0.9 * replicas:
        raise ValueError(f"insufficient valid block bootstrap draws: {len(values)}/{replicas}")
    return values


def paired_nested_bootstrap(rows40: list[dict], rows60: list[dict], labels40: np.ndarray,
                            replicas: int, seed: int) -> dict:
    keys60 = {(r["cluster_id"], r["filter"]) for r in rows60}
    keep60 = np.asarray([(r["cluster_id"], r["filter"]) in keys60 for r in rows40], dtype=bool)
    if int(keep60.sum()) != len(rows60):
        raise ValueError("60% sample does not exactly match a subset of 40% rows")
    row60_by_key = {(r["cluster_id"], r["filter"]): r for r in rows60}
    row40_subset = [row for row, keep in zip(rows40, keep60) if keep]
    if { (r["cluster_id"], r["filter"]) for r in row40_subset } != set(row60_by_key):
        raise ValueError("nested support mismatch")
    blocks = np.unique(labels40)
    by_block = {block: np.flatnonzero(labels40 == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = {"r40": [], "r60": [], "delta": [], "mean_satellite_shift_40_minus_60": []}
    index_to_key = [(r["cluster_id"], r["filter"]) for r in rows40]
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        idx40 = np.concatenate([by_block[block] for block in chosen])
        idx60 = [i for i in idx40 if keep60[i]]
        if len(idx60) < 6:
            continue
        draw40 = [rows40[i] for i in idx40]
        draw60 = [row60_by_key[index_to_key[i]] for i in idx60]
        r40, r60 = partial(draw40), partial(draw60)
        if np.isfinite(r40) and np.isfinite(r60):
            draws["r40"].append(float(r40))
            draws["r60"].append(float(r60))
            draws["delta"].append(float(r40-r60))
            draws["mean_satellite_shift_40_minus_60"].append(float(np.mean(
                [draw40[j]["satellite_delta"] - row60_by_key[index_to_key[i]]["satellite_delta"]
                 for j, i in enumerate(idx40) if keep60[i]])))
    if len(draws["delta"]) < 0.8 * replicas:
        raise ValueError(f"insufficient paired bootstrap draws: {len(draws['delta'])}/{replicas}")
    return {key: {"ci95": ci(value), "valid": len(value)} for key, value in draws.items()}


def analyze(output: Path, replicas: int, seed: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    rows60_by_filter, rows40_by_filter = load_data()
    summaries, pair_rows, center_summaries = {}, [], {}
    for band in ("i", "i2"):
        rows60, rows40 = rows60_by_filter[band], rows40_by_filter[band]
        if not set((r["cluster_id"], r["filter"]) for r in rows60).issubset(
                {(r["cluster_id"], r["filter"]) for r in rows40}):
            raise ValueError(f"60% sample is not nested in the 40% sample for {band}")
        labels60, labels40 = labels_for(rows60), labels_for(rows40)
        boot60 = bootstrap_partial(rows60, labels60, replicas, seed + (1 if band == "i" else 10))
        boot40 = bootstrap_partial(rows40, labels40, replicas, seed + (2 if band == "i" else 20))
        paired = paired_nested_bootstrap(rows40, rows60, labels40, replicas,
                                         seed + (3 if band == "i" else 30))
        paired_keys = {(r["cluster_id"], r["filter"]) for r in rows60}
        sat40_map = {(r["cluster_id"], r["filter"]): r for r in rows40}
        shift = np.asarray([sat40_map[k]["satellite_delta"] - r["satellite_delta"]
                            for r in rows60 if (k := (r["cluster_id"], r["filter"])) in paired_keys])
        summaries[band] = {
            "n_60pct_support": len(rows60), "n_40pct_support": len(rows40),
            "n_added_at_40pct": len(rows40)-len(rows60),
            "n_spatial_blocks_60pct": int(len(np.unique(labels60)),),
            "n_spatial_blocks_40pct": int(len(np.unique(labels40)),),
            "partial_r_60pct": partial(rows60), "partial_r_60pct_block_bootstrap_ci95": ci(boot60),
            "partial_r_40pct": partial(rows40), "partial_r_40pct_block_bootstrap_ci95": ci(boot40),
            "paired_40_minus_60_block_bootstrap": paired,
            "same_60pct_pairs_satellite_delta_shift_40_minus_60": {
                "mean": float(np.mean(shift)), "median": float(np.median(shift)),
                "q16_q84": [float(x) for x in np.percentile(shift, [16, 84])],
            },
            "center_eligible_pairs_60pct": sum(r["center_eligible"] for r in rows60),
            "center_eligible_pairs_40pct": sum(r["center_eligible"] for r in rows40),
            "science_measurement": False,
        }
        center_summaries[band] = {}
        for floor_name, selected in (("60pct", [r for r in rows60 if r["center_eligible"]]),
                                     ("40pct", [r for r in rows40 if r["center_eligible"]])):
            if len(selected) < 6:
                center_summaries[band][floor_name] = {
                    "n_pairs": len(selected), "partial_r": None,
                    "block_bootstrap_ci95": None, "n_spatial_blocks": None,
                    "status": "too_few_pairs_for_diagnostic",
                }
                continue
            selected_labels = labels_for(selected)
            selected_bootstrap = bootstrap_partial(
                selected, selected_labels, replicas,
                seed + (100 if band == "i" else 200) + (0 if floor_name == "60pct" else 1))
            center_summaries[band][floor_name] = {
                "n_pairs": len(selected), "partial_r": partial(selected),
                "block_bootstrap_ci95": ci(selected_bootstrap),
                "n_spatial_blocks": int(len(np.unique(selected_labels))),
                "status": "exploratory_diagnostic",
            }
        m60 = {(r["cluster_id"], r["filter"]): r for r in rows60}
        for row40 in rows40:
            key = (row40["cluster_id"], row40["filter"])
            row60 = m60.get(key)
            pair_rows.append({
                "cluster_id": row40["cluster_id"], "filter": band,
                "in_60pct_support": row60 is not None,
                "icl_delta_200_300_njy_arcsec2": row40["icl_delta"],
                "satellite_delta_50_75_40pct_njy_arcsec2": row40["satellite_delta"],
                "satellite_delta_50_75_60pct_njy_arcsec2": row60["satellite_delta"] if row60 else "",
                "satellite_delta_shift_40_minus_60_njy_arcsec2": row40["satellite_delta"]-row60["satellite_delta"] if row60 else "",
                "n_satellites_40pct": row40["n_satellites"],
                "n_satellites_60pct": row60["n_satellites"] if row60 else "",
                "science_measurement": False,
            })
    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "coverage_threshold_pair_values.tsv"
    write_tsv(pair_path, pair_rows)
    inputs = (ICL_PAIRS, ICL_PROFILES, CONTROL_MANIFEST, SATELLITE_60,
              SATELLITE_60_SUMMARY, SATELLITE_40, SATELLITE_40_SUMMARY)
    summary = {
        "schema": "astrafek-icl-satellite-coverage-threshold-sensitivity-v1",
        "status": "complete_exploratory_sensitivity", "science_measurement": False,
        "hsc_network_used": False, "images_reread": False, "data_requested": False,
        "comparison": "per-annulus 60% target/control coverage floor vs same 800-parent images remeasured at 40%, minimum 3 controls retained",
        "association_estimator": "unweighted partial Pearson of ICL 200-300 kpc target-minus-control surface brightness and satellite 50-75 kpc target-minus-azimuthal-control mean, controlling linearly for redshift and richness",
        "bootstrap": {"method": "5-degree connected blocks; 60-vs-40 nested comparison resamples 40% support blocks and applies identical block multiplicities to the nested 60% sample", "replicas": replicas, "seed": seed},
        "results_by_filter": summaries,
        "center_eligible_sensitivity_by_filter": center_summaries,
        "limitations": [
            "40% support is an estimator sensitivity and may introduce coverage and mask-edge bias; it is not an adopted threshold.",
            "The 60% and 40% values for overlapping pairs can differ because the relaxed threshold admits additional azimuthal controls.",
            "Block intervals describe the selected finite sample and do not correct support-selection bias, PSF wings, sky transfer, completeness, extinction, or source-Poisson variance.",
            "No causal or parent-population inference is made; science_measurement remains false.",
        ],
        "input_sha256": {str(p): sha256(p) for p in inputs},
        "output_sha256": {str(pair_path): sha256(pair_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({"output": str(output.resolve()), "science_measurement": False,
                      "results_by_filter": summaries}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()
    if args.replicas < 500:
        parser.error("at least 500 bootstrap replicas are required")
    analyze(args.output.resolve(), args.replicas, args.seed)


if __name__ == "__main__":
    main()
