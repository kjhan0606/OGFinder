import csv

import numpy as np
import pytest

from audit_satellite_icl_pixel_overlap import selected_pairs, shared_pixel_counts


def test_shared_pixel_counts_distinguishes_geometric_from_valid_overlap():
    satellite = np.array([[1, 1, 1, 0]], dtype=bool)
    icl = np.array([[0, 1, 1, 1]], dtype=bool)
    satellite_valid = np.array([[1, 1, 0, 1]], dtype=bool)
    icl_valid = np.array([[1, 0, 1, 1]], dtype=bool)

    result = shared_pixel_counts(satellite, icl, satellite_valid, icl_valid)

    assert result == {
        "n_satellite_annulus_geometric_pixels": 3,
        "n_satellite_annulus_valid_pixels": 2,
        "n_icl_annulus_valid_pixels_local": 2,
        "n_geometric_intersection_pixels": 2,
        "n_shared_valid_measurement_pixels": 0,
    }


def test_shared_pixel_counts_requires_aligned_crops():
    with pytest.raises(ValueError, match="identical shapes"):
        shared_pixel_counts(np.ones((2, 2), bool), np.ones((2, 1), bool),
                            np.ones((2, 2), bool), np.ones((2, 2), bool))


def test_selected_pairs_applies_both_annulus_support_gates(tmp_path):
    icl = tmp_path / "icl.tsv"
    sat = tmp_path / "sat.tsv"
    with icl.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t",
                                fieldnames=["id", "filter", "coverage_200_300"])
        writer.writeheader()
        writer.writerows([
            {"id": "1", "filter": "i2", "coverage_200_300": 0.5},
            {"id": "2", "filter": "i2", "coverage_200_300": 0.39},
        ])
    fields = ["cluster_id", "filter", "member_index", "detection_number", "cohort",
              "radial_bin_kpc", "target_coverage", "control_coverage", "n_controls"]
    rows = [
        {"cluster_id": "1", "filter": "i2", "member_index": "4", "detection_number": "10",
         "cohort": "full_sample", "radial_bin_kpc": "50-75", "target_coverage": "0.8",
         "control_coverage": "0.7", "n_controls": "4"},
        {"cluster_id": "1", "filter": "i2", "member_index": "5", "detection_number": "11",
         "cohort": "full_sample", "radial_bin_kpc": "50-75", "target_coverage": "0.5",
         "control_coverage": "0.9", "n_controls": "4"},
        {"cluster_id": "2", "filter": "i2", "member_index": "6", "detection_number": "12",
         "cohort": "full_sample", "radial_bin_kpc": "50-75", "target_coverage": "0.9",
         "control_coverage": "0.9", "n_controls": "4"},
    ]
    with sat.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)

    result = selected_pairs(icl, sat)

    assert list(result) == ["1"]
    assert len(result["1"]) == 1
    assert result["1"][0]["member_index"] == "4"
