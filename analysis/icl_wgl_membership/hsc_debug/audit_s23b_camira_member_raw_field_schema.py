#!/usr/bin/env python3
"""Check whether raw CAMIRA member fields reproduce published richness.

This is an offline consistency audit of the already approved and imported
S23B member file. It performs no network access, download, image read, or
science measurement. Empirical consistency cannot replace an authoritative
file-format definition.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE = DATA / "sample_800_2026-10-13"
MEMBER_FILE = SAMPLE / "s23b_camira_member_catalog_2026-11-07_01/camira_s23b_wide_sm_member.dat.gz"
IMPORT_SUMMARY = SAMPLE / "s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/summary.json"
CATALOG = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
OUTPUT = SAMPLE / "s23b_camira_member_schema_audit_2026-11-07_04"
EXPECTED_MEMBER_SHA256 = "66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a"
EXPECTED_CATALOG_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_ROWS = 637_146
EXPECTED_GROUPS = 12_162
RICHNESS_ROUNDING_TOL = 0.00051
MSTAR_LOG10_THRESHOLD = 10.2


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_catalog() -> tuple[dict, dict]:
    by_id, by_coord = {}, {}
    with CATALOG.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            values = line.strip()
            if not values or values.startswith("#"):
                continue
            fields = values.split()
            if len(fields) != 7:
                raise ValueError(f"unexpected catalog row width at line {line_number}")
            row = {
                "catalog_id": fields[0], "name": fields[1],
                "ra_deg": float(fields[2]), "dec_deg": float(fields[3]),
                "z_cl": float(fields[4]), "richness": float(fields[5]),
                "z_spec_bcg": float(fields[6]),
            }
            key = (round(row["ra_deg"], 6), round(row["dec_deg"], 6))
            by_id[row["catalog_id"]] = row
            by_coord.setdefault(key, []).append(row)
    if len(by_id) != EXPECTED_GROUPS:
        raise ValueError(f"catalog has {len(by_id)} rows, expected {EXPECTED_GROUPS}")
    return by_id, by_coord


def qstats(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    q = np.percentile(array, [0, 1, 16, 50, 84, 99, 100])
    return {
        "min": float(q[0]), "q01": float(q[1]), "q16": float(q[2]),
        "median": float(q[3]), "q84": float(q[4]), "q99": float(q[5]),
        "max": float(q[6]),
    }


def main() -> dict:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite schema audit: {OUTPUT}")
    for path in (MEMBER_FILE, IMPORT_SUMMARY, CATALOG):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(MEMBER_FILE) != EXPECTED_MEMBER_SHA256:
        raise ValueError("approved member-file checksum changed")
    if sha256(CATALOG) != EXPECTED_CATALOG_SHA256:
        raise ValueError("pinned cluster-catalog checksum changed")
    importer = json.loads(IMPORT_SUMMARY.read_text(encoding="utf-8"))
    if (importer.get("status") != "complete_catalog_join_and_image_footprint_audit"
            or importer.get("member_file", {}).get("header_present") is not False
            or importer.get("member_file", {}).get("column_count") != 8):
        raise ValueError("approved member-file import record is not the expected verified product")

    catalog_by_id, catalog_by_coord = read_catalog()
    group_weights: dict[str, list[float]] = {}
    group_weights_above_mass_threshold: dict[str, list[float]] = {}
    group_field7: dict[str, list[float]] = {}
    group_rows: dict[str, int] = {}
    raw7, raw8 = [], []
    unmatched_examples = []
    n_rows = 0
    n_weight_outside_unit = 0
    n_weight_nonpositive = 0
    n_field7_below_10p2 = 0

    with gzip.open(MEMBER_FILE, "rt", encoding="ascii", newline="") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"member row {line_number} has {len(fields)} fields")
            values = [float(value) for value in fields]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"member row {line_number} has a nonfinite field")
            ra_cluster, dec_cluster, row_richness, row_z = values[:4]
            field7, field8 = values[6], values[7]
            if (not 0 <= ra_cluster < 360 or not -90 <= dec_cluster <= 90
                    or not 0 <= values[4] < 360 or not -90 <= values[5] <= 90):
                raise ValueError(f"member row {line_number} has invalid coordinates")
            candidates = catalog_by_coord.get((round(ra_cluster, 6), round(dec_cluster, 6)), [])
            matches = [row for row in candidates
                       if abs(row["richness"] - row_richness) <= RICHNESS_ROUNDING_TOL
                       and abs(row["z_cl"] - row_z) <= RICHNESS_ROUNDING_TOL]
            if len(matches) != 1:
                if len(unmatched_examples) < 10:
                    unmatched_examples.append({
                        "line": line_number, "coordinate_candidates": len(candidates),
                        "signature_matches": len(matches),
                    })
                continue
            catalog_id = matches[0]["catalog_id"]
            group_weights.setdefault(catalog_id, []).append(field8)
            if field7 > MSTAR_LOG10_THRESHOLD:
                group_weights_above_mass_threshold.setdefault(catalog_id, []).append(field8)
            group_field7.setdefault(catalog_id, []).append(field7)
            group_rows[catalog_id] = group_rows.get(catalog_id, 0) + 1
            raw7.append(field7)
            raw8.append(field8)
            n_rows += 1
            n_weight_outside_unit += int(not 0.0 <= field8 <= 1.0)
            n_weight_nonpositive += int(field8 <= 0.0)
            n_field7_below_10p2 += int(field7 < 10.2)

    if unmatched_examples:
        raise ValueError(f"member file has ambiguous/unmatched rows: {unmatched_examples}")
    if n_rows != EXPECTED_ROWS or len(group_rows) != EXPECTED_GROUPS:
        raise ValueError(f"matched {n_rows} rows / {len(group_rows)} groups; expected {EXPECTED_ROWS}/{EXPECTED_GROUPS}")

    records = []
    diffs = []
    threshold_diffs = []
    within_rounding = 0
    threshold_within_rounding = 0
    for catalog_id in sorted(catalog_by_id, key=lambda value: int(value)):
        cat = catalog_by_id[catalog_id]
        weights = group_weights[catalog_id]
        threshold_weights = group_weights_above_mass_threshold.get(catalog_id, [])
        field7s = group_field7[catalog_id]
        weight_sum = math.fsum(weights)
        delta = weight_sum - cat["richness"]
        threshold_weight_sum = math.fsum(threshold_weights)
        threshold_delta = threshold_weight_sum - cat["richness"]
        diffs.append(delta)
        threshold_diffs.append(threshold_delta)
        within_rounding += int(abs(delta) <= RICHNESS_ROUNDING_TOL)
        threshold_within_rounding += int(abs(threshold_delta) <= RICHNESS_ROUNDING_TOL)
        records.append({
            "catalog_id": catalog_id,
            "name": cat["name"],
            "richness": cat["richness"],
            "z_cl": cat["z_cl"],
            "n_member_rows": group_rows[catalog_id],
            "field8_sum": weight_sum,
            "field8_sum_minus_richness": delta,
            "abs_delta_within_0p00051": abs(delta) <= RICHNESS_ROUNDING_TOL,
            "n_member_rows_field7_gt_10p2": len(threshold_weights),
            "field8_sum_field7_gt_10p2": threshold_weight_sum,
            "field8_sum_field7_gt_10p2_minus_richness": threshold_delta,
            "field8_threshold_sum_abs_delta_within_0p00051": abs(threshold_delta) <= RICHNESS_ROUNDING_TOL,
            "field8_min": min(weights),
            "field8_median": statistics.median(weights),
            "field8_max": max(weights),
            "field7_min": min(field7s),
            "field7_median": statistics.median(field7s),
            "field7_max": max(field7s),
        })

    OUTPUT.mkdir(parents=True, exist_ok=False)
    rows_path = OUTPUT / "cluster_field_sums.tsv"
    import csv
    with rows_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    abs_diffs = [abs(value) for value in diffs]
    rel_diffs = [abs(value) / max(1e-12, catalog_by_id[record["catalog_id"]]["richness"])
                 for value, record in zip(diffs, records)]
    sum_count = sum(record["n_member_rows"] for record in records)
    summary = {
        "schema": "astrafek-s23b-camira-member-raw-field-schema-audit-v1",
        "status": "complete_empirical_consistency_audit",
        "science_measurement": False,
        "network_used": False,
        "images_read": False,
        "inputs_sha256": {
            str(MEMBER_FILE): sha256(MEMBER_FILE),
            str(IMPORT_SUMMARY): sha256(IMPORT_SUMMARY),
            str(CATALOG): sha256(CATALOG),
        },
        "published_sources": [
            {
                "title": "Oguri et al. 2026, HSC S23B CAMIRA catalog paper",
                "url": "https://arxiv.org/html/2512.13954v1",
                "relevant_claims": [
                    "Wide star-mask CAMIRA catalog has 12,162 clusters.",
                    "The member catalog includes more than 600,000 galaxies with w_mem >= 1e-3.",
                    "Richness applies stellar-mass and spatial filters and a compensated background filter; it is not the unqualified sum of exported weights.",
                ],
            },
            {
                "title": "Official Oguri cluster_catalogs S23B README",
                "url": "https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md",
                "relevant_claims": ["Documents the cluster catalog files; does not give the member-file column order."],
            },
        ],
        "n_member_rows_joined": n_rows,
        "n_catalog_clusters_joined": len(group_rows),
        "cluster_groups_missing": len(set(catalog_by_id) - set(group_rows)),
        "field7_raw_distribution": qstats(raw7),
        "field8_raw_distribution": qstats(raw8),
        "field8_rows_outside_0_1": n_weight_outside_unit,
        "field8_nonpositive_rows": n_weight_nonpositive,
        "field8_rows_at_or_above_1e_3": sum(value >= 1.0e-3 for value in raw8),
        "field7_rows_below_10p2": n_field7_below_10p2,
        "published_member_definition_consistency": {
            "published_wide_star_mask_cluster_count": 12_162,
            "observed_member_file_rows": n_rows,
            "published_member_sample_description": "more than 600,000 member galaxies with w_mem >= 10^-3",
            "observed_field8_min": min(raw8),
            "observed_field8_max": max(raw8),
            "all_field8_values_at_least_1e_3": all(value >= 1.0e-3 for value in raw8),
            "interpretation": "The row count and field8 range match the published w_mem >= 1e-3 member sample description; this strongly supports, but does not authoritatively prove, that raw field 8 is w_mem.",
        },
        "field8_sum_vs_catalog_richness": {
            "catalog_richness_definition_in_publication": "membership-probability weighted number of galaxies above approximately 10^10.2 solar masses",
            "threshold_test": "sum raw field 8 only where raw field 7 > 10.2; field 7 is tested as log10 stellar mass in solar units, not assigned that meaning from this consistency check alone",
            "clusters_with_field8_sum_within_0p00051_of_catalog_richness": within_rounding,
            "fraction_within_rounding": within_rounding / len(records),
            "absolute_difference_quantiles": qstats(abs_diffs),
            "relative_absolute_difference_quantiles": qstats(rel_diffs),
            "max_absolute_difference": max(abs_diffs),
            "clusters_with_thresholded_field8_sum_within_0p00051_of_catalog_richness": threshold_within_rounding,
            "thresholded_fraction_within_rounding": threshold_within_rounding / len(records),
            "thresholded_absolute_difference_quantiles": qstats([abs(value) for value in threshold_diffs]),
            "thresholded_relative_absolute_difference_quantiles": qstats([
                abs(value) / max(1e-12, catalog_by_id[record["catalog_id"]]["richness"])
                for value, record in zip(threshold_diffs, records)
            ]),
            "thresholded_max_absolute_difference": max(abs(value) for value in threshold_diffs),
            "n_member_rows_field7_gt_10p2": sum(record["n_member_rows_field7_gt_10p2"] for record in records),
            "sum_member_row_counts": sum_count,
        },
        "interpretation": (
            "Field 8 has the published member-weight range and threshold behavior, so it is strongly consistent with w_mem. Its sum is not equal to catalog richness because published N_mem is computed through stellar-mass and spatial filters plus compensated background subtraction, not as a simple sum of the exported membership weights. "
            "The public paper and catalog README do not specify this headerless file's column order, so field 8 remains a probable w_mem candidate pending direct file-format confirmation."
        ),
        "limitations": [
            "The published paper defines w_mem and the member-sample threshold, but the public catalog README lists cluster catalog files and does not map the member-file column order.",
            "The richness comparison is a diagnostic only; the CAMIRA richness estimator also applies mass/spatial filters and compensated background subtraction.",
            "The file remains headerless; authoritative member-file schema or source-code confirmation is still needed before scientific labeling.",
        ],
        "outputs": {"cluster_field_sums_tsv": str(rows_path)},
        "output_sha256": {rows_path.name: sha256(rows_path)},
    }
    summary_path = OUTPUT / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "rows": result["n_member_rows_joined"],
        "clusters": result["n_catalog_clusters_joined"],
        "field8_sum_matches_richness": result["field8_sum_vs_catalog_richness"]["clusters_with_field8_sum_within_0p00051_of_catalog_richness"],
        "field8_outside_unit": result["field8_rows_outside_0_1"],
        "output": result["outputs"]["cluster_field_sums_tsv"],
    }, ensure_ascii=False))
