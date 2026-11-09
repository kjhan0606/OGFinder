#!/usr/bin/env python3
"""Compare forced CAMIRA satellite stacks across coverage/mask settings.

This offline analysis uses the existing forced-position outputs at 60% and
20% annular coverage. It measures (a) whether masking other unmatched CAMIRA
positions changed the 60% result, and (b) the paired spatial-block change in
the ICL--satellite correlation when the 20% measurement is added.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
MASKED_60 = DATA / "satellite_forced_outer_profile_2026-11-06_01"
UNMASKED_60 = DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_2026-11-06_01"
UNMASKED_20 = DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_2026-11-06_01"
TARGETS = DATA / "download_01/download_manifest.json"
CONTROLS = DATA / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_forced_coverage_sensitivity_2026-11-06_02"
REPLICATES = 10000
SEED = 20261106


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value: str) -> float | None:
    return None if value in ("", "None", "null") else float(value)


def partial_r(rows: list[dict], y_key: str) -> float:
    if len(rows) < 6:
        return float("nan")
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[y_key]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.log(np.asarray([float(row["richness"]) for row in rows]))
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if np.linalg.matrix_rank(design) < 3:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.linalg.norm(rx) == 0 or np.linalg.norm(ry) == 0:
        return float("nan")
    return float(np.dot(rx, ry) / (np.linalg.norm(rx) * np.linalg.norm(ry)))


def bootstrap_pairwise_difference(
    low_full: list[dict], low_common: list[dict], high_common: list[dict],
    band: str, target_manifest: dict, control_manifest: dict,
    replicas: int, seed: int,
) -> dict:
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    control_by_target = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    low_full = sorted(low_full, key=lambda row: int(row["cluster_id"]))
    common_keys = {(row["cluster_id"], row["filter"]) for row in low_common}
    low_common_map = {(row["cluster_id"], row["filter"]): row for row in low_common}
    high_common_map = {(row["cluster_id"], row["filter"]): row for row in high_common}
    low_common = [low_common_map[key] for key in sorted(common_keys, key=lambda k: int(k[0]))]
    high_common = [high_common_map[key] for key in sorted(common_keys, key=lambda k: int(k[0]))]
    if len(low_common) != len(high_common) or len(low_common) != len(common_keys):
        raise ValueError(f"{band}: paired 60/20 common sample mismatch")

    target_coords, control_coords = [], []
    for row in low_full:
        target = target_by_id[row["cluster_id"]]
        control = control_by_target.get(row["cluster_id"])
        if control is None:
            raise ValueError(f"{band}: missing matched control for {row['cluster_id']}")
        target_coords.append({"ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"])})
        control_coords.append({"ra_deg": float(control["ra_deg"]), "dec_deg": float(control["dec_deg"])})
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(low_full), dtype=bool))
    groups = sorted(set(labels))
    by_block_full = {group: np.flatnonzero(labels == group) for group in groups}
    key_to_full = {(row["cluster_id"], row["filter"]): i for i, row in enumerate(low_full)}
    by_block_common = {
        group: np.asarray([
            i for i, row in enumerate(low_common)
            if labels[key_to_full[(row["cluster_id"], row["filter"])]] == group
        ], dtype=int)
        for group in groups
    }

    full_rows_by_index = low_full
    low_rows_by_index = low_common
    high_rows_by_index = high_common
    rng = np.random.default_rng(seed)
    changes_expansion, changes_measurement, changes_total = [], [], []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        full_idx = np.concatenate([by_block_full[group] for group in chosen])
        common_idx = np.concatenate([by_block_common[group] for group in chosen])
        if len(common_idx) < 6:
            continue
        r20_full = partial_r([full_rows_by_index[i] for i in full_idx], "satellite_delta_combined_50_75")
        r20_common = partial_r([low_rows_by_index[i] for i in common_idx], "satellite_delta_combined_50_75")
        r60_common = partial_r([high_rows_by_index[i] for i in common_idx], "satellite_delta_combined_50_75")
        if np.isfinite([r20_full, r20_common, r60_common]).all():
            changes_expansion.append(r20_full - r20_common)
            changes_measurement.append(r20_common - r60_common)
            changes_total.append(r20_full - r60_common)
    if len(changes_total) < 0.9 * replicas:
        raise ValueError(f"{band}: insufficient paired bootstrap replicates {len(changes_total)}/{replicas}")

    def summarize(values: list[float]) -> dict:
        array = np.asarray(values)
        return {"median": float(np.median(array)),
                "ci95": np.quantile(array, [0.025, 0.975]).tolist(),
                "fraction_positive": float(np.mean(array > 0))}

    return {
        "n_20pct_pairs": len(low_full), "n_common_pairs": len(low_common),
        "n_connected_5deg_blocks_union": n_blocks,
        "n_target_cells_union": n_target, "n_control_cells_union": n_control,
        "largest_union_block_pairs": largest,
        "partial_r_20pct_full": partial_r(low_full, "satellite_delta_combined_50_75"),
        "partial_r_20pct_common": partial_r(low_common, "satellite_delta_combined_50_75"),
        "partial_r_60pct_common": partial_r(high_common, "satellite_delta_combined_50_75"),
        "paired_change_from_sample_expansion_20pct_full_minus_common": summarize(changes_expansion),
        "paired_change_from_measurement_20pct_minus_60pct_common": summarize(changes_measurement),
        "paired_change_total_20pct_full_minus_60pct_common": summarize(changes_total),
        "replicates_requested": replicas, "replicates_valid": len(changes_total), "seed": seed,
        "bootstrap_note": "shared connected 5-degree target/control blocks defined on the 20-percent sample; exploratory paired uncertainty",
    }


def run(output: Path, replicas: int = REPLICATES) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing existing output directory: {output}")
    summaries = {name: read_json(path / "summary.json")
                 for name, path in (("masked_60", MASKED_60),
                                    ("unmasked_60", UNMASKED_60),
                                    ("unmasked_20", UNMASKED_20))}
    if any(summary.get("status") != "complete_exploratory"
           or summary.get("science_measurement") is not False for summary in summaries.values()):
        raise ValueError("one forced-position run is incomplete or has an unexpected science flag")
    if (summaries["unmasked_60"].get("method", {}).get("mask_other_unmatched_camira") is not False
            or summaries["unmasked_20"].get("method", {}).get("mask_other_unmatched_camira") is not False
            or summaries["unmasked_20"].get("coverage_floor") != 0.2):
        raise ValueError("sensitivity run settings do not match the intended scope")

    target_manifest = read_json(TARGETS)
    control_manifest = read_json(CONTROLS)
    catalogs = {}
    for name, path in (("masked_60", MASKED_60), ("unmasked_60", UNMASKED_60), ("unmasked_20", UNMASKED_20)):
        diagnostics = read_tsv(path / "candidate_diagnostics.tsv")
        associations = read_tsv(path / "icl_satellite_association.tsv")
        profiles = read_tsv(path / "forced_satellite_profiles.tsv")
        catalogs[name] = {"diagnostics": diagnostics, "associations": associations, "profiles": profiles}
        if len(diagnostics) != 3197:
            raise ValueError(f"{name}: forced candidate diagnostics are incomplete")

    masked_rows = {(row["cluster_id"], row["filter"], row["member_index"], row["radial_bin_kpc"]): row
                   for row in catalogs["masked_60"]["profiles"]}
    unmasked_rows = {(row["cluster_id"], row["filter"], row["member_index"], row["radial_bin_kpc"]): row
                     for row in catalogs["unmasked_60"]["profiles"]}
    common_profile_keys = set(masked_rows) & set(unmasked_rows)
    mask_profile_diffs = {field: [] for field in
                          ("delta_sb_njy_arcsec2", "target_coverage", "control_coverage")}
    for key in sorted(common_profile_keys):
        a, b = masked_rows[key], unmasked_rows[key]
        for field in ("delta_sb_njy_arcsec2", "target_coverage", "control_coverage"):
            mask_profile_diffs[field].append(abs(float(a[field]) - float(b[field])))
    coverage_diffs = []
    diagnostics_by = []
    for a, b in zip(catalogs["masked_60"]["diagnostics"], catalogs["unmasked_60"]["diagnostics"]):
        if (a["cluster_id"], a["filter"], a["member_index"]) != (b["cluster_id"], b["filter"], b["member_index"]):
            raise ValueError("mask comparison candidate order/key mismatch")
        ca = np.asarray([float(v) for v in a["target_coverage_by_bin"].split(",")])
        cb = np.asarray([float(v) for v in b["target_coverage_by_bin"].split(",")])
        coverage_diffs.append(float(np.max(np.abs(ca - cb))))
        diagnostics_by.append({"filter": a["filter"], "cluster_id": a["cluster_id"],
                               "member_index": a["member_index"],
                               "masked_60_retained": int(a["retained_annuli"]),
                               "unmasked_60_retained": int(b["retained_annuli"])})
    if len(common_profile_keys) != len(masked_rows) or len(common_profile_keys) != len(unmasked_rows):
        raise ValueError("mask setting changed the 60-percent retained profile key set")

    association_by_name = {
        name: {(row["cluster_id"], row["filter"]): row for row in data["associations"]}
        for name, data in catalogs.items()
    }
    association_results = {}
    for band_index, band in enumerate(("i", "i2")):
        maps = {name: {key: row for key, row in table.items()
                       if key[1] == band and number(row["satellite_delta_combined_50_75"]) is not None}
                for name, table in association_by_name.items()}
        full20_keys = sorted(maps["unmasked_20"], key=lambda key: int(key[0]))
        common_keys = sorted(set(maps["unmasked_20"]) & set(maps["unmasked_60"])
                             & set(maps["masked_60"]), key=lambda key: int(key[0]))
        low_full = [maps["unmasked_20"][key] for key in full20_keys]
        low_common = [maps["unmasked_20"][key] for key in common_keys]
        high_common = [maps["unmasked_60"][key] for key in common_keys]
        association_results[band] = bootstrap_pairwise_difference(
            low_full, low_common, high_common, band, target_manifest, control_manifest,
            replicas, 20261106 + band_index)
        association_results[band]["n_60pct_common_profiles"] = len(maps["unmasked_60"])
        association_results[band]["n_20pct_full_profiles"] = len(maps["unmasked_20"])
        association_results[band]["partial_r_60pct"] = partial_r(
            [maps["unmasked_60"][key] for key in sorted(maps["unmasked_60"], key=lambda k: int(k[0]))],
            "satellite_delta_combined_50_75")

    output.mkdir(parents=True, exist_ok=False)
    diagnostic = {
        "status": "complete_exploratory", "science_measurement": False,
        "scope": "existing HSC cutouts and masks only; no HSC request/download",
        "settings": {name: {"mask_other_unmatched_camira": summaries[name].get("method", {}).get("mask_other_unmatched_camira", True),
                            "coverage_floor": summaries[name].get("coverage_floor", 0.6)} for name in summaries},
        "mask_sensitivity_60pct": {
            "n_common_position_bin_keys": len(common_profile_keys),
            "max_abs_differences": {field: max(values, default=0.0)
                                     for field, values in mask_profile_diffs.items()},
            "max_abs_candidate_coverage_vector_difference": max(coverage_diffs, default=0.0),
            "note": "the 60-percent support cutoff is unchanged; tests whether fallback masking of neighboring unmatched CAMIRA positions explains coverage loss",
        },
        "association_by_filter": association_results,
        "limitations": [
            "The 20-percent coverage threshold is a post-hoc exploratory sensitivity, not a validated science cut.",
            "CAMIRA positions without SEP detections are not confirmed galaxies; the signal is mean light around selected catalog positions.",
            "Selection and missing-at-random are unresolved; no multiplicity correction across analysis choices is applied.",
            "The profile still masks detected sources and uses approximate 5-degree spatial blocks.",
            "PSF wings, sky transfer, source-Poisson noise, coadd covariance, and full-aperture depth/mask matching remain uncalibrated.",
        ],
        "input_sha256": {str(path): sha256(path) for directory in (MASKED_60, UNMASKED_60, UNMASKED_20)
                          for path in (directory / "summary.json", directory / "candidate_diagnostics.tsv",
                                       directory / "forced_satellite_profiles.tsv", directory / "icl_satellite_association.tsv")}
                        | {str(path): sha256(path) for path in (TARGETS, CONTROLS)},
        "replicates": replicas,
    }
    (output / "summary.json").write_text(json.dumps(diagnostic, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return diagnostic


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    if args.replicates < 1000:
        parser.error("at least 1,000 replicates are required")
    result = run(args.output, args.replicates)
    print(json.dumps({"status": result["status"],
                      "mask_sensitivity_60pct": result["mask_sensitivity_60pct"],
                      "association_by_filter": result["association_by_filter"],
                      "science_measurement": result["science_measurement"]}, indent=2))


if __name__ == "__main__":
    main()
