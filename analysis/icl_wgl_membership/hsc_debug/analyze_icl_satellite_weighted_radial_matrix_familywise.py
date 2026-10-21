#!/usr/bin/env python3
"""Joint family-wise test of equal-satellite and CAMIRA-weighted radial profiles."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_icl_satellite_radial_matrix_familywise import (
    DEFAULT_EXTRACTION, DEFAULT_ICL, DEFAULT_SATELLITE, run_bootstrap,
)
from audit_icl_satellite_radial_support_matrix import DATA, sha256


DEFAULT_MEMBER = DATA / "satellite_member_weight_sensitivity_2026-10-18_04"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def aggregate_member_profiles(detail_rows: list[dict[str, str]]) -> dict:
    accumulators = defaultdict(lambda: {"sum_w_delta": 0.0, "sum_w": 0.0,
                                       "sum_delta": 0.0, "n": 0})
    seen = set()
    for row in detail_rows:
        key = (str(row["cluster_id"]), row["filter"], row["radial_bin_kpc"],
               int(row["member_index"]))
        if key in seen:
            raise ValueError(f"duplicate CAMIRA-weighted member profile {key}")
        seen.add(key)
        delta, weight = float(row["delta"]), float(row["camira_weight"])
        if not math.isfinite(delta) or not math.isfinite(weight) or not 0 <= weight <= 1:
            raise ValueError(f"invalid member brightness/weight {key}")
        group_key = key[:3]
        value = accumulators[group_key]
        value["sum_w_delta"] += weight * delta
        value["sum_w"] += weight
        value["sum_delta"] += delta
        value["n"] += 1
    profiles = {}
    for key, item in accumulators.items():
        if item["sum_w"] <= 0:
            raise ValueError(f"zero total CAMIRA weight for profile {key}")
        profiles[key] = {
            "weighted_mean": item["sum_w_delta"] / item["sum_w"],
            "equal_mean": item["sum_delta"] / item["n"],
            "n_satellites": item["n"],
            "sum_weight": item["sum_w"],
        }
    return profiles


def load_weighted_profiles(directory: Path):
    directory = Path(directory).resolve()
    summary_path = directory / "summary.json"
    detail_path = directory / "profile_member_properties.tsv"
    aggregate_path = directory / "member_weighted_profiles.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "hsc-camira-weighted-satellite-profile-audit-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("column_mapping", {}).get("certainty") !=
            "column 8 is identified as weight; column 7 mass interpretation remains provisional"):
        raise ValueError("unexpected CAMIRA member-weighted source or schema certainty")
    if summary.get("acquisition_schema", {}).get("weight_field") != \
            "column 8, explicitly labeled weight by the original fetch/parser":
        raise ValueError("CAMIRA weight source is not explicitly identified")
    for path in (detail_path, aggregate_path):
        if summary.get("outputs_sha256", {}).get(path.name) != sha256(path):
            raise ValueError(f"member-weighted profile hash mismatch: {path}")
    for source, digest in summary.get("input_sha256", {}).items():
        if not Path(source).is_file() or sha256(Path(source)) != digest:
            raise ValueError(f"member-weighted input hash mismatch: {source}")
    member_summary = Path(summary["acquisition_schema"]["source_summary"])
    if sha256(member_summary) != summary["acquisition_schema"]["source_summary_sha256"]:
        raise ValueError("CAMIRA member acquisition summary hash mismatch")

    detail_rows = read_tsv(detail_path)
    profiles = aggregate_member_profiles(detail_rows)

    # Independently reproduce the previously saved estimator-level aggregate.
    aggregate_rows = [row for row in read_tsv(aggregate_path)
                      if row["mass_selection"] == "all"]
    for row in aggregate_rows:
        band, radial = row["filter"], row["radial_bin_kpc"]
        selected = [v for (cluster, f, r), v in profiles.items() if f == band and r == radial]
        if not selected:
            raise ValueError(f"no member-weighted profile aggregate for {band}/{radial}")
        eq = float(np.mean([v["equal_mean"] for v in selected]))
        wt = float(np.mean([v["weighted_mean"] for v in selected]))
        if (int(row["n_clusters"]) != len(selected)
                or int(row["n_satellites"]) != sum(v["n_satellites"] for v in selected)
                or not math.isclose(float(row["equal_cluster_mean_delta_njy_arcsec2"]), eq,
                                    rel_tol=1e-9, abs_tol=1e-8)
                or not math.isclose(float(row["camira_weighted_equal_cluster_mean_delta_njy_arcsec2"]), wt,
                                    rel_tol=1e-9, abs_tol=1e-8)):
            raise ValueError(f"saved aggregate does not reproduce from member weights: {band}/{radial}")
    return profiles, {
        "summary": {"path": str(summary_path), "sha256": sha256(summary_path)},
        "detail_table": {"path": str(detail_path), "sha256": sha256(detail_path)},
        "aggregate_table": {"path": str(aggregate_path), "sha256": sha256(aggregate_path)},
        "source_inputs": {
            source: {"path": source, "sha256": digest}
            for source, digest in summary.get("input_sha256", {}).items()
        },
        "n_member_profile_rows": len(detail_rows),
        "n_cluster_filter_radial_profiles": len(profiles),
        "n_zero_weight_profiles": 0,
        "weight_field": "CAMIRA column 8 named weight by original hash-matched acquisition parser",
    }


def paired_estimator_groups(weighted_profiles: dict, icl_dir: Path,
                            satellite_dir: Path, extraction_path: Path):
    equal_groups, support_records, base_provenance = __import__(
        "analyze_icl_satellite_radial_matrix_familywise", fromlist=["load_inputs"]
    ).load_inputs(icl_dir, satellite_dir, extraction_path)
    combined, support = {}, []
    for key, equal_group in equal_groups.items():
        band, icl_bin, satellite_bin = key
        equal_group_key = (*key, "equal_satellite")
        combined[equal_group_key] = equal_group
        support.append({"estimator": "equal_satellite", **next(
            row for row in support_records if (row["filter"], row["icl_radial_bin_kpc"],
                                                row["satellite_radial_bin_kpc"]) == key)})

        weighted_rows = []
        for base in equal_group["rows"]:
            w = weighted_profiles.get((str(base["cluster_id"]), band, satellite_bin))
            if w is None:
                raise ValueError(f"CAMIRA-weighted profile missing for common support {key}/{base['cluster_id']}")
            if w["n_satellites"] != int(base["n_satellites"]):
                raise ValueError(f"member count changed between equal/weighted profiles {key}/{base['cluster_id']}")
            weighted_rows.append({**base, "satellite_delta_50_75": w["weighted_mean"]})
        block_by_row = {id(row): str(label)
                        for row, label in zip(weighted_rows, equal_group["blocks"])}
        from familywise_block_wild_association import prepare_group
        weighted_group, weighted_support = prepare_group(weighted_rows, block_by_row, min_blocks=8)
        if weighted_group is None:
            raise ValueError(f"CAMIRA-weighted cell lost support: {key}/{weighted_support}")
        weighted_key = (*key, "camira_weighted")
        combined[weighted_key] = weighted_group
        support.append({"estimator": "camira_weighted", "filter": band,
                        "icl_radial_bin_kpc": icl_bin,
                        "satellite_radial_bin_kpc": satellite_bin,
                        "n_common": weighted_support["n_clusters"],
                        **weighted_support})
    if len(combined) != 168 or len(support) != 168:
        raise ValueError(f"expected 168 estimator/radial/filter cells, got {len(combined)}")
    return combined, support, base_provenance


def analyze(member_dir: Path, icl_dir: Path, satellite_dir: Path,
            extraction_path: Path, output_dir: Path, *, replicas: int, seed: int):
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output: {output_dir}")
    weighted_profiles, weight_provenance = load_weighted_profiles(member_dir)
    groups, support, base_provenance = paired_estimator_groups(
        weighted_profiles, icl_dir, satellite_dir, extraction_path)
    bootstrap = run_bootstrap(groups, support, replicas=replicas, seed=seed)
    output_dir.mkdir(parents=True, exist_ok=False)
    result_rows = [
        {"satellite_estimator": key[3], **row}
        for key, row in zip(groups, bootstrap["results"])
    ]
    association_path = output_dir / "weighted_radial_matrix_association.tsv"
    fields = ["satellite_estimator", *[key for key in result_rows[0] if key != "satellite_estimator"]]
    with association_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result_rows)
    support_path = output_dir / "weighted_radial_matrix_support.tsv"
    support_fields = list(support[0])
    with support_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=support_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(support)
    null_path = output_dir / "joint_maxT_null.tsv"
    with null_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["replica", "max_abs_t"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({"replica": i + 1, "max_abs_t": float(value)}
                         for i, value in enumerate(bootstrap["max_t_draws"]))
    result = {
        "schema": "astrafek-icl-satellite-weighted-radial-matrix-familywise-v1",
        "status": "complete_exploratory_familywise_sensitivity",
        "science_measurement": False,
        "method": {
            "estimators": ["equal satellite weight within cluster", "CAMIRA column-8 weight within cluster"],
            "cluster_weighting": "equal weight across clusters after within-cluster member aggregation",
            "model": "satellite excess ~ ICL excess + redshift + CAMIRA richness; cluster-level OLS with CR1 spatial-block standard errors",
            "null_bootstrap": "restricted null without ICL predictor; HC2 residual adjustment; Webb six-point wild bootstrap-t",
            "spatial_blocks": "common 5-degree target/control sky-cell connected components",
            "family": "all eligible 6 ICL x 7 satellite annuli x 2 filters x 2 satellite estimators, 168 tests",
            "multiplicity": "single-step max absolute t under complete null and Holm-adjusted raw bootstrap p-values across all 168 cells",
            "support_rule": "at least 10 cluster/filter rows and at least 8 connected spatial blocks",
            "replicas": replicas, "seed": seed,
        },
        "provenance": {"icl_satellite": base_provenance, "camira_weighted_profiles": weight_provenance},
        "n_expected_cells": 168, "n_eligible_cells": bootstrap["n_eligible_cells"],
        "output_sha256": {str(path.resolve()): sha256(path)
                           for path in (association_path, support_path, null_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "CAMIRA weight includes the algorithm's stellar-mass and spatial filters and is not an independent unbiased membership probability",
            "all 168 cells reuse the same selected images and have different radial support; missing-not-at-random selection is not corrected",
            "5-degree connected blocks approximate spatial dependence; some i-band tests have only 8 blocks",
            "max-T is calibrated under the complete null; strong family-wise control is not claimed without subset pivotality",
            "Holm correction uses approximate few-cluster wild-bootstrap marginal p-values",
            "PSF wings, sky transfer, satellite completeness, extinction, passband response, and full image/coadd covariance remain unresolved",
            "CAMIRA weight and HSC image profiles do not establish a causal tidal-stripping measurement",
        ],
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--member-dir", type=Path, default=DEFAULT_MEMBER)
    parser.add_argument("--icl-dir", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261022)
    args = parser.parse_args()
    result = analyze(args.member_dir, args.icl_dir, args.satellite_dir,
                     args.extraction_summary, args.output_dir,
                     replicas=args.replicas, seed=args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "cells": result["n_eligible_cells"], "output": str(args.output_dir.resolve())},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
