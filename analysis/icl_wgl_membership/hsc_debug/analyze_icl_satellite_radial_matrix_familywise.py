#!/usr/bin/env python3
"""Family-wise spatial wild-bootstrap audit of the ICL x satellite radial matrix.

All 6x7 radial pairings in both filters are tested as one exploratory family.
The script uses the exact target/control pairs from the existing image stack,
one common five-degree spatial partition, and a restricted Webb wild-cluster
bootstrap. It does not promote the science-measurement gate.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_icl_satellite_radial_support_matrix import (
    DATA, FILTERS, ICL_BINS, SATELLITE_BINS, read_tsv, sha256, validate_summary,
)
from evaluate_prestack_pair_support import pair_block_labels
from familywise_block_wild_association import (
    WEBB_WEIGHTS, bootstrap_t, holm_adjust, prepare_group,
)


DEFAULT_ICL = DATA / "icl_annulus_support_2026-10-18_01"
DEFAULT_SATELLITE = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_EXTRACTION = DATA / "prestack_extraction_2026-10-15_01/extraction_summary.json"


def load_inputs(icl_dir: Path, satellite_dir: Path, extraction_path: Path):
    icl_dir, satellite_dir, extraction_path = map(Path, (icl_dir, satellite_dir, extraction_path))
    icl_summary_path, satellite_summary_path = icl_dir / "summary.json", satellite_dir / "summary.json"
    icl_table, satellite_table = icl_dir / "annulus_pair_profiles.tsv", satellite_dir / "cluster_profiles.tsv"
    icl_summary = validate_summary(icl_summary_path, "hsc-icl-per-annulus-support-v1", icl_table,
                                   require_exploratory=False)
    satellite_summary = validate_summary(satellite_summary_path,
                                         "hsc-satellite-outer-profile-per-annulus-v1",
                                         satellite_table, require_exploratory=True)
    if satellite_summary.get("n_frames_completed") != 800:
        raise ValueError("satellite input is not the approved full 800-frame profile set")
    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    if (extraction.get("status") != "complete"
            or str(extraction.get("science_measurement", "")).lower() != "false"
            or len(extraction.get("records", [])) != int(extraction.get("n_measured", -1))):
        raise ValueError("matched-control extraction summary is incomplete or promoted")
    controls = {}
    for record in extraction["records"]:
        if record.get("frame_role") != "matched_control":
            continue
        control_id = str(record["control_id"])
        if control_id in controls:
            raise ValueError(f"duplicate matched-control coordinate: {control_id}")
        ra, dec = float(record["ra_deg"]), float(record["dec_deg"])
        if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
            raise ValueError(f"invalid coordinate for matched control {control_id}")
        controls[control_id] = record

    icl_rows = [row for row in read_tsv(icl_table) if row["radial_bin_kpc"] in ICL_BINS]
    satellite_rows = [row for row in read_tsv(satellite_table)
                      if row.get("cohort") == "full_sample" and row["radial_bin_kpc"] in SATELLITE_BINS]
    satellite_by_key = {}
    for row in satellite_rows:
        key = (str(row["cluster_id"]), row["filter"], row["radial_bin_kpc"])
        if key in satellite_by_key:
            raise ValueError(f"duplicate satellite profile {key}")
        satellite_by_key[key] = row

    pair_info = {}
    icl_by_key = {}
    for row in icl_rows:
        band, radial, cluster = row["filter"], row["radial_bin_kpc"], str(row["id"])
        if band not in FILTERS:
            raise ValueError(f"unexpected ICL filter {band}")
        key = (cluster, band, radial)
        if key in icl_by_key:
            raise ValueError(f"duplicate ICL profile {key}")
        icl_by_key[key] = row
        pair_key = (cluster, band)
        control_id = str(row["control_id"])
        control = controls.get(control_id)
        if control is None:
            raise ValueError(f"missing validated control coordinate for {control_id}")
        if control.get("header_filter") != band:
            raise ValueError(f"target/control filter mismatch for {pair_key}")
        metadata = {
            "cluster_id": cluster, "filter": band, "control_id": control_id,
            "target_ra_deg": float(row["ra_deg"]), "target_dec_deg": float(row["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]), "control_dec_deg": float(control["dec_deg"]),
            "redshift": float(row["z"]), "richness": float(row["richness"]),
        }
        if pair_key in pair_info:
            old = pair_info[pair_key]
            for field in metadata:
                if field not in ("cluster_id", "filter") and old[field] != metadata[field]:
                    raise ValueError(f"target/control metadata changes across ICL radii for {pair_key}/{field}")
        else:
            pair_info[pair_key] = metadata

    pair_keys = sorted(pair_info)
    targets = [{"ra_deg": pair_info[key]["target_ra_deg"],
                "dec_deg": pair_info[key]["target_dec_deg"]} for key in pair_keys]
    controls_for_pairs = [{"ra_deg": pair_info[key]["control_ra_deg"],
                           "dec_deg": pair_info[key]["control_dec_deg"]} for key in pair_keys]
    indices, labels, n_target_cells, n_control_cells, n_blocks, largest_block = pair_block_labels(
        targets, controls_for_pairs, np.ones(len(pair_keys), dtype=bool))
    if len(indices) != len(pair_keys) or len(labels) != len(pair_keys):
        raise ValueError("spatial block construction changed pair order")
    block_by_pair = dict(zip(pair_keys, map(str, labels)))

    groups, records = {}, []
    for band in FILTERS:
        for icl_bin in ICL_BINS:
            for satellite_bin in SATELLITE_BINS:
                key = (band, icl_bin, satellite_bin)
                rows = []
                for (cluster, current_band), info in pair_info.items():
                    if current_band != band:
                        continue
                    icl = icl_by_key.get((cluster, band, icl_bin))
                    sat = satellite_by_key.get((cluster, band, satellite_bin))
                    if icl is None or sat is None:
                        continue
                    row = {
                        "icl_delta_200_300": float(icl["delta_sb_njy_arcsec2"]),
                        "satellite_delta_50_75": float(sat["mean_delta_sb_njy_arcsec2"]),
                        "z": info["redshift"], "richness": info["richness"],
                        "n_satellites": int(sat["n_satellites"]),
                        "cluster_id": cluster, "filter": band,
                        "icl_radial_bin_kpc": icl_bin,
                        "satellite_radial_bin_kpc": satellite_bin,
                    }
                    if not all(math.isfinite(float(row[field])) for field in (
                            "icl_delta_200_300", "satellite_delta_50_75", "z", "richness")):
                        raise ValueError(f"non-finite analysis value for {key}/{cluster}")
                    if row["n_satellites"] < 1:
                        raise ValueError(f"satellite radial row has no member for {key}/{cluster}")
                    rows.append(row)
                if not rows:
                    raise ValueError(f"no common support for {key}")
                block_by_row = {id(row): block_by_pair[(row["cluster_id"], band)] for row in rows}
                group, support = prepare_group(rows, block_by_row, min_blocks=8)
                records.append({
                    "filter": band, "icl_radial_bin_kpc": icl_bin,
                    "satellite_radial_bin_kpc": satellite_bin,
                    "n_common": len(rows), **support,
                })
                if group is not None:
                    groups[key] = group
    return groups, records, {
        "icl_summary": {"path": str(icl_summary_path.resolve()), "sha256": sha256(icl_summary_path)},
        "icl_table": {"path": str(icl_table.resolve()), "sha256": sha256(icl_table)},
        "satellite_summary": {"path": str(satellite_summary_path.resolve()), "sha256": sha256(satellite_summary_path)},
        "satellite_table": {"path": str(satellite_table.resolve()), "sha256": sha256(satellite_table)},
        "control_extraction_summary": {"path": str(extraction_path.resolve()), "sha256": sha256(extraction_path)},
        "n_unique_target_control_pairs": len(pair_keys),
        "n_target_5deg_cells": n_target_cells,
        "n_control_5deg_cells": n_control_cells,
        "n_global_connected_blocks": n_blocks,
        "largest_global_block_pairs": largest_block,
    }


def run_bootstrap(groups: dict, records: list[dict], *, replicas: int, seed: int) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 wild-bootstrap replicates are required")
    if len(groups) < 2:
        raise ValueError("fewer than two radial cells meet the cluster/block support rule")
    global_blocks = sorted({str(block) for group in groups.values() for block in group["blocks"]})
    rng = np.random.default_rng(seed)
    max_t_draws = np.empty(replicas, dtype=float)
    group_draws = {key: np.empty(replicas, dtype=float) for key in groups}
    for replica in range(replicas):
        weights = dict(zip(global_blocks, rng.choice(WEBB_WEIGHTS, size=len(global_blocks), replace=True)))
        max_abs = 0.0
        for key, group in groups.items():
            value = bootstrap_t(group, weights)
            group_draws[key][replica] = value
            if not math.isfinite(value):
                raise ValueError(f"non-finite bootstrap t for {key}, replica {replica}")
            max_abs = max(max_abs, abs(value))
        max_t_draws[replica] = max_abs
    keys = list(groups)
    observed = np.asarray([groups[key]["t_observed"] for key in keys])
    raw_p = np.asarray([
        (1 + np.count_nonzero(np.abs(group_draws[key]) >= abs(groups[key]["t_observed"])))
        / (replicas + 1) for key in keys
    ])
    max_t_p = np.asarray([
        (1 + np.count_nonzero(max_t_draws >= abs(groups[key]["t_observed"])))
        / (replicas + 1) for key in keys
    ])
    adjusted = holm_adjust(raw_p)
    results = []
    for index, key in enumerate(keys):
        group = groups[key]
        results.append({
            "filter": key[0], "icl_radial_bin_kpc": key[1],
            "satellite_radial_bin_kpc": key[2], "status": "exploratory_inference",
            "n_clusters": group["n_clusters"],
            "n_spatial_blocks": group["n_spatial_blocks"],
            "largest_spatial_block": group["largest_spatial_block"],
            "partial_pearson_r_z_richness": group["partial_r"],
            "beta_icl_njy_arcsec2_per_njy_arcsec2": group["beta_icl"],
            "t_cluster_cr1": group["t_observed"],
            "p_wild_cluster_raw": raw_p[index],
            "p_holm_all_eligible_cells": adjusted[index],
            "p_maxT_single_step_complete_null": max_t_p[index],
        })
    return {"results": results, "max_t_draws": max_t_draws,
            "n_eligible_cells": len(groups), "n_all_cells": len(records)}


def analyze(icl_dir: Path, satellite_dir: Path, extraction_path: Path,
            output_dir: Path, *, replicas: int, seed: int):
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output: {output_dir}")
    groups, support_records, provenance = load_inputs(icl_dir, satellite_dir, extraction_path)
    bootstrap = run_bootstrap(groups, support_records, replicas=replicas, seed=seed)
    output_dir.mkdir(parents=True, exist_ok=False)
    association_path = output_dir / "radial_matrix_association.tsv"
    fields = list(bootstrap["results"][0])
    with association_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(bootstrap["results"])
    support_path = output_dir / "radial_matrix_support.tsv"
    support_fields = list(support_records[0])
    with support_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=support_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(support_records)
    null_path = output_dir / "maxT_null.tsv"
    with null_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["replica", "max_abs_t"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({"replica": i + 1, "max_abs_t": float(value)}
                         for i, value in enumerate(bootstrap["max_t_draws"]))
    result = {
        "schema": "astrafek-icl-satellite-radial-matrix-familywise-wild-v1",
        "status": "complete_exploratory_familywise_inference",
        "science_measurement": False,
        "method": {
            "outcome": "satellite target-minus-control mean surface brightness per radial annulus",
            "predictor": "ICL target-minus-control surface brightness per radial annulus",
            "model": "satellite excess ~ ICL excess + redshift + CAMIRA richness; cluster-level OLS with CR1 spatial-block standard errors",
            "null_bootstrap": "restricted null without ICL predictor; HC2 residual adjustment; Webb six-point wild bootstrap-t",
            "spatial_blocks": "common 5-degree target/control sky-cell connected components across all 84 radial/filter tests",
            "family": "all eligible 6 ICL annuli x 7 satellite annuli x 2 filters as one family",
            "multiplicity": "single-step max absolute t under complete null and Holm-adjusted raw bootstrap p-values",
            "support_rule": "at least 10 cluster/filter rows and at least 8 connected spatial blocks",
            "replicas": replicas, "seed": seed,
        },
        "n_expected_cells": len(support_records),
        "n_eligible_cells": bootstrap["n_eligible_cells"],
        "provenance": provenance,
        "output_sha256": {str(path.resolve()): sha256(path)
                           for path in (association_path, support_path, null_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "the 5-degree connected cells approximate spatial dependence; fewer than 10 blocks is a weak basis for wild-cluster inference",
            "max-T calibration is a complete-null diagnostic; strong family-wise control is not claimed because subset pivotality is unverified",
            "Holm correction assumes valid marginal p-values; finite-block Webb bootstrap remains approximate for the smallest support groups",
            "radial/filter cells have different selected support and are not corrected for missing-not-at-random selection",
            "unresolved PSF wings, sky transfer, satellite completeness, extinction, passband response, and image/coadd covariance keep the result exploratory",
            "correlation is not causal evidence for tidal stripping or ICL production",
        ],
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl-dir", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    result = analyze(args.icl_dir, args.satellite_dir, args.extraction_summary,
                     args.output_dir, replicas=args.replicas, seed=args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "eligible_cells": result["n_eligible_cells"],
                      "expected_cells": result["n_expected_cells"],
                      "output": str(args.output_dir.resolve())}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
