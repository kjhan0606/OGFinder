import csv

from audit_satellite_icl_annulus_overlap import (
    annulus_intersection_area,
    disk_intersection_area,
    analyze,
)


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_disk_intersection_disjoint_contained_and_identical():
    assert disk_intersection_area(1, 2, 4) == 0
    assert disk_intersection_area(1, 3, 1) == 3.141592653589793
    assert abs(disk_intersection_area(2, 2, 0) - 4 * 3.141592653589793) < 1e-12


def test_annulus_overlap_only_when_center_distance_brackets_shell():
    # The annuli can touch for 125 <= R <= 375 kpc.
    assert annulus_intersection_area(124.999) == 0
    assert annulus_intersection_area(375.001) == 0
    assert annulus_intersection_area(200) > 0
    assert annulus_intersection_area(300) > 0


def test_annulus_overlap_is_bounded_by_both_annulus_areas():
    area = annulus_intersection_area(250)
    satellite_area = 3.141592653589793 * (75**2 - 50**2)
    icl_area = 3.141592653589793 * (300**2 - 200**2)
    assert 0 <= area <= min(satellite_area, icl_area)


def test_analysis_joins_only_same_cluster_filter_and_reports_exploratory(tmp_path):
    icl = tmp_path / "icl.tsv"
    sat = tmp_path / "sat.tsv"
    write_tsv(icl, ["id", "filter", "coverage_200_300"], [
        {"id": "a", "filter": "i", "coverage_200_300": "0.5"},
        {"id": "b", "filter": "i", "coverage_200_300": "0.3"},
    ])
    write_tsv(sat, ["cluster_id", "filter", "member_index", "cluster_radius_kpc",
                    "radial_bin_kpc", "target_coverage", "control_coverage"], [
        {"cluster_id": "a", "filter": "i", "member_index": "1", "cluster_radius_kpc": "250",
         "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"},
        {"cluster_id": "a", "filter": "i2", "member_index": "2", "cluster_radius_kpc": "250",
         "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"},
        {"cluster_id": "b", "filter": "i", "member_index": "3", "cluster_radius_kpc": "250",
         "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"},
    ])
    out = analyze(icl, sat)
    assert out["n_icl_cluster_filter_pairs"] == 1
    assert out["n_matched_satellite_profiles"] == 1
    assert out["science_measurement"] is False
    assert out["flux_covariance_available"] is False


def test_duplicate_satellite_bin_is_rejected(tmp_path):
    icl = tmp_path / "icl.tsv"
    sat = tmp_path / "sat.tsv"
    write_tsv(icl, ["id", "filter", "coverage_200_300"], [
        {"id": "a", "filter": "i", "coverage_200_300": "0.5"},
    ])
    row = {"cluster_id": "a", "filter": "i", "member_index": "1", "cluster_radius_kpc": "250",
           "radial_bin_kpc": "50-75", "target_coverage": "1", "control_coverage": "1"}
    write_tsv(sat, list(row), [row, row])
    try:
        analyze(icl, sat)
    except ValueError as exc:
        assert "duplicate satellite" in str(exc)
    else:
        raise AssertionError("duplicate profile row was not rejected")
