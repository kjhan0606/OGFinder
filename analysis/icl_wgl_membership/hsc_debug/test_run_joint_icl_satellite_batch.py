import csv

import pytest

from run_joint_icl_satellite_batch import (
    select_jointly_supported_pairs,
    summarize,
)


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_selection_intersects_common_icl_and_satellite_support(tmp_path):
    icl = tmp_path / "icl.tsv"
    satellite = tmp_path / "satellite.tsv"
    write_tsv(icl, ["id", "filter", "coverage_200_300"], [
        {"id": "1", "filter": "i", "coverage_200_300": ".5"},
        {"id": "2", "filter": "i2", "coverage_200_300": ".39"},
        {"id": "3", "filter": "i2", "coverage_200_300": ".8"},
    ])
    fields = ["cluster_id", "filter", "radial_bin_kpc", "target_coverage", "control_coverage"]
    write_tsv(satellite, fields, [
        {"cluster_id": "1", "filter": "i", "radial_bin_kpc": "50-75", "target_coverage": ".6", "control_coverage": ".6"},
        {"cluster_id": "2", "filter": "i2", "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"},
        {"cluster_id": "3", "filter": "i2", "radial_bin_kpc": "50-75", "target_coverage": ".6", "control_coverage": ".59"},
        {"cluster_id": "3", "filter": "i2", "radial_bin_kpc": "40-50", "target_coverage": "1", "control_coverage": "1"},
    ])
    assert select_jointly_supported_pairs(icl, satellite) == [("1", "i")]


def test_duplicate_icl_pair_is_rejected(tmp_path):
    icl = tmp_path / "icl.tsv"
    satellite = tmp_path / "satellite.tsv"
    row = {"id": "1", "filter": "i", "coverage_200_300": ".5"}
    write_tsv(icl, list(row), [row, row])
    write_tsv(satellite, ["cluster_id", "filter", "radial_bin_kpc", "target_coverage", "control_coverage"], [
        {"cluster_id": "1", "filter": "i", "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"}
    ])
    with pytest.raises(ValueError, match="duplicate ICL"):
        select_jointly_supported_pairs(icl, satellite)


def test_summary_reports_filter_level_medians_and_keeps_exploratory_gate():
    records = [{"filter": "i", "transfer_response": {
        "icl_response_to_diffuse_unit_injection": x,
        "icl_plane_corrected_response_to_diffuse_unit_injection": x / 2,
        "icl_response_to_all_satellite_envelopes": 2*x,
        "icl_plane_corrected_response_to_all_satellite_envelopes": x,
        "satellite_response_to_smooth_diffuse": 0.1*x,
        "satellite_response_to_satellite_unit_injection": 0.5*x,
    }} for x in (1.0, 3.0, 5.0)]
    result = summarize(records)
    assert result["i"]["n_clusters"] == 3
    response = result["i"]["responses"]["icl_response_to_diffuse_unit_injection"]
    assert response["median"] == 3.0
    assert response["q16"] == pytest.approx(1.64)
    assert response["q84"] == pytest.approx(4.36)
