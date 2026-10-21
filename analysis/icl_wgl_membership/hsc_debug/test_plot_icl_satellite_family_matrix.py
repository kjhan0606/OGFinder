import csv
import json

import pytest

from plot_icl_satellite_family_matrix import load_family, sha256


def _write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_rejects_unverified_or_incomplete_family_inputs(tmp_path):
    summary = tmp_path / "summary.json"
    association = tmp_path / "association.tsv"
    support = tmp_path / "support.tsv"
    summary.write_text(json.dumps({"schema": "wrong", "science_measurement": False}))
    association.write_text("x\n1\n")
    support.write_text("x\n1\n")
    with pytest.raises(ValueError, match="summary schema"):
        load_family(summary, association, support)


def test_requires_exploratory_gate_before_plotting(tmp_path):
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "schema": "astrafek-icl-satellite-weighted-radial-matrix-familywise-v1",
        "science_measurement": True,
    }))
    with pytest.raises(ValueError, match="remain exploratory"):
        load_family(summary, tmp_path / "missing.tsv", tmp_path / "also-missing.tsv")


def test_loads_hash_pinned_complete_radial_family(tmp_path):
    summary_path = tmp_path / "summary.json"
    association_path = tmp_path / "association.tsv"
    support_path = tmp_path / "support.tsv"
    assoc_fields = [
        "satellite_estimator", "filter", "icl_radial_bin_kpc", "satellite_radial_bin_kpc",
        "status", "n_clusters", "n_spatial_blocks", "partial_pearson_r_z_richness",
        "p_holm_all_eligible_cells", "p_maxT_single_step_complete_null",
    ]
    support_fields = [
        "estimator", "filter", "icl_radial_bin_kpc", "satellite_radial_bin_kpc",
        "n_clusters", "status",
    ]
    icl_bins = ("50-100", "100-150", "150-200", "200-300", "300-400", "400-500")
    sat_bins = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
    association, support = [], []
    for estimator in ("equal_satellite", "camira_weighted"):
        for filter_name in ("i", "i2"):
            for icl_bin in icl_bins:
                for sat_bin in sat_bins:
                    n = 20
                    association.append({
                        "satellite_estimator": estimator, "filter": filter_name,
                        "icl_radial_bin_kpc": icl_bin, "satellite_radial_bin_kpc": sat_bin,
                        "status": "exploratory_inference", "n_clusters": n,
                        "n_spatial_blocks": 8, "partial_pearson_r_z_richness": 0.1,
                        "p_holm_all_eligible_cells": 1.0,
                        "p_maxT_single_step_complete_null": 0.8,
                    })
                    support.append({
                        "estimator": estimator, "filter": filter_name,
                        "icl_radial_bin_kpc": icl_bin, "satellite_radial_bin_kpc": sat_bin,
                        "n_clusters": n, "status": "eligible",
                    })
    _write_tsv(association_path, assoc_fields, association)
    _write_tsv(support_path, support_fields, support)
    summary_path.write_text(json.dumps({
        "schema": "astrafek-icl-satellite-weighted-radial-matrix-familywise-v1",
        "science_measurement": False, "n_expected_cells": 168,
        "output_sha256": {
            str(association_path): sha256(association_path),
            str(support_path): sha256(support_path),
        },
    }))
    summary, assoc, sup, actual_icl_bins, actual_sat_bins = load_family(
        summary_path, association_path, support_path)
    assert summary["science_measurement"] is False
    assert len(assoc) == len(sup) == 168
    assert actual_icl_bins == list(icl_bins)
    assert actual_sat_bins == list(sat_bins)
