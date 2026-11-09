#!/usr/bin/env python3
"""Decompose the saved i-band ICL-split cross signal by HSC tract component.

Uses previously approved target and random shape-row exports. It assigns each
lens/center and source to the existing approximate connected tract components,
reapplies the calibrated estimator cuts, and reports regional point profiles
and boundary-crossing fractions. It is descriptive only: no region covariance
or calibrated significance is inferred.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_wl_random_subtracted_icl_split import (
    NZ_SACC, RANDOM_PLAN, RANDOM_ROWS, SOURCE_BINS, active_lens_ids, subset,
)
from fit_wl_random_stack_nfw_mass import DATA, RANDOM_SUMMARY, load_pairs_and_provenance
from ogfkit.wlcluster import Cosmology
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels, tract_centers_from_patch_qa, tract_region_labels,
)
from run_wl_random_shape_stack import verify_and_load_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc, pair_geometry, stack_s19a_tomographic_pairs,
    tomographic_nz_lensing_efficiency,
)


REGION_SUMMARY = DATA / "wl_hsc_tract_region_influence_2026-10-20_01/summary.json"
PATCH_QA = DATA / "catalogues_2026-10-13_01/patch_qa.csv"
SPLIT_SUMMARY = DATA / "wl_random_subtracted_icl_split_i_2026-11-06_01/summary.json"
DEFAULT_OUTPUT = DATA / "wl_icl_split_source_region_cross_null_2026-11-06_01"
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
MASKED_B_MODE_VALUE = 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def labels_for_ids(pairs: dict, selected_ids: set[str], id_column: str,
                   ra_column: str, dec_column: str, centers: np.ndarray,
                   tract_regions: np.ndarray) -> dict[str, int]:
    ids = np.asarray(pairs[id_column]).astype(str)
    selected = np.isin(ids, list(selected_ids))
    if not np.any(selected):
        raise ValueError(f"no rows found for selected IDs in {id_column}")
    locations: dict[str, tuple[float, float]] = {}
    for lens_id, ra, dec in zip(ids[selected],
                                np.asarray(pairs[ra_column])[selected],
                                np.asarray(pairs[dec_column])[selected]):
        point = (float(ra), float(dec))
        previous = locations.setdefault(str(lens_id), point)
        if abs(previous[0] - point[0]) > 1e-9 or abs(previous[1] - point[1]) > 1e-9:
            raise ValueError(f"endpoint coordinates vary across rows for {lens_id}")
    found = nearest_region_labels(
        [locations[key][0] for key in locations],
        [locations[key][1] for key in locations], centers, tract_regions,
    )
    return dict(zip(locations, map(int, found)))


def usable_mask(pairs: dict, nz_z: np.ndarray, nz_bins: np.ndarray,
                cosmology: Cosmology) -> tuple[np.ndarray, dict]:
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    beta = tomographic_nz_lensing_efficiency(
        pairs["lens_z"], source_bin, nz_z, nz_bins, cosmology,
    )
    radius, _phi = pair_geometry(pairs, cosmology)
    shape = np.hypot(pairs["e1"], pairs["e2"])
    keep = (
        np.asarray(pairs["b_mode_pass"], dtype=bool)
        & (pairs["sigma_e"] >= 0.0) & (pairs["sigma_e"] <= 0.4)
        & (shape < 2.0) & (pairs["mag_a10"] <= 25.5)
        & (pairs["resolution"] >= 0.3) & (beta > 0.0)
        & (radius >= RADIAL_EDGES_MPC[0])
        & (radius <= RADIAL_EDGES_MPC[-1])
    )
    return keep, {
        "input_pairs": int(len(keep)),
        "used_pairs": int(np.count_nonzero(keep)),
        "source_zbin_counts": {
            str(k): int(np.count_nonzero(keep & (source_bin == k)))
            for k in SOURCE_BINS
        },
    }


def profile_rows(pairs: dict, nz_z: np.ndarray, nz_bins: np.ndarray,
                 cosmology: Cosmology) -> tuple[list[dict], dict]:
    stacked = stack_s19a_tomographic_pairs(
        pairs, nz_z, nz_bins, RADIAL_EDGES_MPC, cosmology,
    )
    return stacked["profile"], stacked["counts"]


def run(output: Path = DEFAULT_OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite source-region audit: {output}")

    split = json.loads(SPLIT_SUMMARY.read_text(encoding="utf-8"))
    region_doc = json.loads(REGION_SUMMARY.read_text(encoding="utf-8"))
    random_doc = json.loads(RANDOM_SUMMARY.read_text(encoding="utf-8"))
    if (split.get("status") != "complete_exploratory"
            or split.get("science_measurement") is not False
            or split.get("filter") != "i"
            or split.get("source_tomographic_bins") != list(SOURCE_BINS)):
        raise ValueError("unexpected saved i-band ICL/WL split")
    if (region_doc.get("status") != "exploratory_field_influence"
            or random_doc.get("status") not in {"exploratory", "complete_exploratory"}):
        raise ValueError("parent region or random stack summary is not complete")

    target_pairs, random_stack_doc, _tomo = load_pairs_and_provenance()
    if sha256(RANDOM_SUMMARY) != "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8":
        raise ValueError("random-subtracted source-row stack provenance changed")
    if random_stack_doc.get("schema") != random_doc.get("schema"):
        raise ValueError("random-row and random-stack summary schema mismatch")
    random_pairs, random_provenance = verify_and_load_pairs(
        RANDOM_ROWS, RANDOM_PLAN, masked_value=MASKED_B_MODE_VALUE,
    )
    target_pairs = subset(target_pairs, np.isin(target_pairs["source_zbin"], SOURCE_BINS))
    random_pairs = subset(random_pairs, np.isin(random_pairs["source_zbin"], SOURCE_BINS))

    matched = split["matched_lens_pairs"]
    low_ids = {str(row["lower_lens_id"]) for row in matched}
    high_ids = {str(row["upper_lens_id"]) for row in matched}
    if len(matched) != split["n_matched_lens_pairs"] or low_ids & high_ids:
        raise ValueError("saved split matched-lens list is invalid")
    selected_ids = low_ids | high_ids

    center_to_lens = {
        str(k): str(v) for k, v in
        random_doc["input"]["random_center_matched_lens_ids"].items()
    }
    centers_for_lenses = {
        center for center, lens_id in center_to_lens.items()
        if lens_id in selected_ids
    }
    count_by_lens: dict[str, int] = {}
    for center in centers_for_lenses:
        lens_id = center_to_lens[center]
        count_by_lens[lens_id] = count_by_lens.get(lens_id, 0) + 1
    if set(count_by_lens.values()) != {5} or set(count_by_lens) != selected_ids:
        raise ValueError("selected lenses do not have exactly five approved random centers")

    nz_z, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    active_lenses = active_lens_ids(target_pairs, nz_z, nz_bins, cosmology)
    active_random_centers = active_lens_ids(random_pairs, nz_z, nz_bins, cosmology)
    active_random_lenses = {center_to_lens[c] for c in active_random_centers if c in center_to_lens}
    if not selected_ids <= active_lenses & active_random_lenses:
        raise ValueError("matched split contains a lens without active target and random support")

    tract_keys, centers = tract_centers_from_patch_qa(PATCH_QA)
    components = tract_region_labels(centers, adjacency_deg=2.0)
    target_region = labels_for_ids(
        target_pairs, selected_ids, "lens_id", "lens_ra_deg", "lens_dec_deg",
        centers, components,
    )
    random_region = labels_for_ids(
        random_pairs, centers_for_lenses, "lens_id", "lens_ra_deg", "lens_dec_deg",
        centers, components,
    )
    source_region_target = nearest_region_labels(
        target_pairs["source_ra_deg"], target_pairs["source_dec_deg"], centers, components,
    )
    source_region_random = nearest_region_labels(
        random_pairs["source_ra_deg"], random_pairs["source_dec_deg"], centers, components,
    )

    parent_lens_region = {}
    for item in region_doc["region_deletions"]:
        for lens_id in item["deleted_lens_ids"]:
            parent_lens_region[str(lens_id)] = int(item["region_id"])
    if any(parent_lens_region.get(lens_id) != region for lens_id, region in target_region.items()):
        raise ValueError("recomputed i-lens region assignment differs from the parent audit")

    id_target = np.asarray(target_pairs["lens_id"]).astype(str)
    id_random = np.asarray(random_pairs["lens_id"]).astype(str)
    random_lens_for_center = np.asarray([center_to_lens.get(x, "") for x in id_random])
    group_ids = {"low_icl": low_ids, "high_icl": high_ids}
    all_used = {}
    selected_target_rows = np.isin(id_target, list(selected_ids))
    selected_random_rows = np.isin(id_random, list(centers_for_lenses))
    endpoint_arrays = {
        "target": subset(target_pairs, selected_target_rows),
        "random": subset(random_pairs, selected_random_rows),
    }
    endpoint_ids = {
        "target": np.asarray(endpoint_arrays["target"]["lens_id"]).astype(str),
        "random": np.asarray(endpoint_arrays["random"]["lens_id"]).astype(str),
    }
    endpoint_lens_region = {"target": target_region, "random": random_region}
    endpoint_source_region = {
        "target": source_region_target[selected_target_rows],
        "random": source_region_random[selected_random_rows],
    }
    endpoint_group_ids = {
        "target": group_ids,
        "random": {
            name: {center for center, lens_id in center_to_lens.items() if lens_id in ids}
            for name, ids in group_ids.items()
        },
    }

    for endpoint, arrays in endpoint_arrays.items():
        all_used[endpoint] = {}
        keep_used, counts = usable_mask(arrays, nz_z, nz_bins, cosmology)
        if counts["used_pairs"] != profile_rows(arrays, nz_z, nz_bins, cosmology)[1]["used_pairs"]:
            raise ValueError(f"manual selection mask disagrees with stack estimator for {endpoint}")
        all_used[endpoint]["all"] = {
            "selection_counts": counts,
            "lens_source_same_component_fraction": float(np.mean(
                np.asarray([endpoint_lens_region[endpoint][str(x)] for x in endpoint_ids[endpoint][keep_used]])
                == endpoint_source_region[endpoint][keep_used]
            )),
        }
        all_used[endpoint]["by_group_region"] = {}
        for group_name, ids in endpoint_group_ids[endpoint].items():
            if endpoint == "random":
                selected_row = np.isin(endpoint_ids[endpoint], list(ids))
            else:
                selected_row = np.isin(endpoint_ids[endpoint], list(ids))
            all_used[endpoint]["by_group_region"][group_name] = {}
            for region in sorted(set(target_region.values())):
                if endpoint == "random":
                    use_endpoint_id = endpoint_ids[endpoint]
                else:
                    use_endpoint_id = endpoint_ids[endpoint]
                endpoint_regions = np.asarray([
                    endpoint_lens_region[endpoint].get(str(value), -1)
                    for value in use_endpoint_id
                ], dtype=int)
                selected = selected_row & (endpoint_regions == region)
                regional = subset(arrays, selected)
                profile, stack_counts = profile_rows(regional, nz_z, nz_bins, cosmology)
                regional_used, selection_counts = usable_mask(regional, nz_z, nz_bins, cosmology)
                if selection_counts["used_pairs"] != stack_counts["used_pairs"]:
                    raise ValueError(f"manual selection mask disagrees for {endpoint}/{group_name}/{region}")
                lens_regions = np.asarray([
                    endpoint_lens_region[endpoint][str(value)]
                    for value in np.asarray(regional["lens_id"]).astype(str)[regional_used]
                ], dtype=int)
                source_regions = endpoint_source_region[endpoint][selected][regional_used]
                same_fraction = (float(np.mean(lens_regions == source_regions))
                                 if len(lens_regions) else None)
                all_used[endpoint]["by_group_region"][group_name][str(region)] = {
                    "selected_endpoint_count": len(set(np.asarray(regional["lens_id"]).astype(str))),
                    "selection_counts": selection_counts,
                    "lens_source_same_component_fraction": same_fraction,
                    "profile": [{
                        "r_min_mpc": row["r_min_mpc"],
                        "r_max_mpc": row["r_max_mpc"],
                        "tangential_msun_mpc2": row["delta_sigma_weak_msun_mpc2"],
                        "cross_msun_mpc2": row["cross_delta_sigma_weak_msun_mpc2"],
                        "used_pairs": row["n_pairs"],
                    } for row in profile],
                }

    regional_contrasts = {}
    regions_present = sorted(set(target_region.values()))
    for region in regions_present:
        low_target = all_used["target"]["by_group_region"]["low_icl"][str(region)]["profile"]
        high_target = all_used["target"]["by_group_region"]["high_icl"][str(region)]["profile"]
        low_random = all_used["random"]["by_group_region"]["low_icl"][str(region)]["profile"]
        high_random = all_used["random"]["by_group_region"]["high_icl"][str(region)]["profile"]
        regional_contrasts[str(region)] = []
        for rows in zip(low_target, high_target, low_random, high_random):
            lt, ht, lr, hr = rows
            regional_contrasts[str(region)].append({
                "r_min_mpc": lt["r_min_mpc"],
                "r_max_mpc": lt["r_max_mpc"],
                "high_minus_low_tangential_random_subtracted_msun_mpc2": (
                    ht["tangential_msun_mpc2"] - hr["tangential_msun_mpc2"]
                    - lt["tangential_msun_mpc2"] + lr["tangential_msun_mpc2"]
                ),
                "high_minus_low_cross_random_subtracted_msun_mpc2": (
                    ht["cross_msun_mpc2"] - hr["cross_msun_mpc2"]
                    - lt["cross_msun_mpc2"] + lr["cross_msun_mpc2"]
                ),
            })

    result = {
        "schema": "astrafek-wl-icl-split-source-region-cross-null-v1",
        "status": "complete_descriptive_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "new_data_requested": False,
        "images_read": False,
        "shape_rows_read": True,
        "filter": "i",
        "source_tomographic_bins": list(SOURCE_BINS),
        "method": "reapply the existing HSC Y3 estimator cuts, assign lens/random-center and source coordinates to the same approximate nearest-tract connected component, then stack target and five matched random centers separately by group and endpoint component",
        "n_pairs_in_high_low_split": len(matched),
        "lens_counts_by_component": {
            "low_icl": {str(region): sum(target_region[lens_id] == region for lens_id in low_ids)
                        for region in regions_present},
            "high_icl": {str(region): sum(target_region[lens_id] == region for lens_id in high_ids)
                         for region in regions_present},
        },
        "matched_pair_region_crossings": sum(
            target_region[str(row["lower_lens_id"])] != target_region[str(row["upper_lens_id"])]
            for row in matched
        ),
        "estimator_counts_and_region_profiles": all_used,
        "regional_random_subtracted_high_minus_low": regional_contrasts,
        "covariance": "not estimated; regional profiles are descriptive, source reuse and field covariance are not propagated",
        "limitations": [
            "tract components are approximate nearest-center groups, not official HSC field masks",
            "regional decomposition conditions on the same 10-pair ICL/WL matched sample",
            "regional profiles have no covariance or uncertainty estimate",
            "cross-shear structure does not by itself identify a physical or instrumental cause",
        ],
        "provenance": {
            "split_summary_sha256": sha256(SPLIT_SUMMARY),
            "parent_region_summary_sha256": sha256(REGION_SUMMARY),
            "random_stack_summary_sha256": sha256(RANDOM_SUMMARY),
            "random_export_plan_sha256": random_provenance["random_row_plan_sha256"],
            "random_export_batch_hashes": random_provenance["random_row_batch_sha256"],
            "patch_qa_sha256": sha256(PATCH_QA),
            "hsc_y3_nz_sha256": sha256(NZ_SACC),
        },
    }
    output.mkdir(parents=True)
    path = output / "summary.json"
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(path)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "lens_counts_by_component": result["lens_counts_by_component"],
        "matched_pair_region_crossings": result["matched_pair_region_crossings"],
        "regional_cross_profiles": {
            region: [round(row["high_minus_low_cross_random_subtracted_msun_mpc2"] / 1e12, 3)
                     for row in rows]
            for region, rows in regional_contrasts.items()
        },
        "summary": str(path),
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
