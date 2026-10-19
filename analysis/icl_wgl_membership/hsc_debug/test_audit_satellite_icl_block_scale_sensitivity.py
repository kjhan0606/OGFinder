import numpy as np
import pytest

from audit_satellite_icl_block_scale_sensitivity import (
    _bootstrap,
    BANDS,
    CATEGORIES,
    compare_five_degree_baseline,
    connected_pair_blocks,
    output_fields,
)
from audit_satellite_icl_geometric_strata import connected_labels


def pair(index, ra, control_ra, satellite):
    return {
        "cluster_id": str(index), "filter": "i2",
        "icl_delta_200_300": str(index * 0.1 + 1.0),
        "satellite_delta": str(satellite),
        "z": str(0.3 + index * 0.002 + index * index * 0.00003),
        "richness": str(10 + index + 0.07 * (index % 3)), "target_ra_deg": str(ra),
        "target_dec_deg": "0", "control_ra_deg": str(control_ra),
        "control_dec_deg": "0",
    }


def test_five_degree_components_match_existing_definition():
    rows = [pair(0, 1, 101, 1.2), pair(1, 4, 109, 0.8),
            pair(2, 12, 121, 1.5)]
    assert np.array_equal(connected_pair_blocks(rows, 5), connected_labels(rows))
    assert len(np.unique(connected_pair_blocks(rows, 5))) == 2
    assert len(np.unique(connected_pair_blocks(rows, 2.5))) == 3


def test_block_bootstrap_is_seeded_and_reports_support():
    rows = [pair(i, i * 10.0, 180 + i * 10.0, 0.5 * i + (-1) ** i)
            for i in range(12)]
    first = _bootstrap(rows, "satellite_delta", 5.0, 1000,
                       np.random.default_rng(19), 6)
    second = _bootstrap(rows, "satellite_delta", 5.0, 1000,
                        np.random.default_rng(19), 6)
    assert first == second
    assert first["n_spatial_blocks"] == 12
    assert first["replicas_valid"] == 1000
    assert first["status"] == "complete_exploratory"


def test_invalid_width_and_sky_coordinates_are_rejected():
    rows = [pair(0, 1, 101, 1.0)]
    with pytest.raises(ValueError, match="block width"):
        connected_pair_blocks(rows, 0)
    rows[0]["target_ra_deg"] = "360"
    with pytest.raises(ValueError, match="RA"):
        connected_pair_blocks(rows, 5)


def test_output_schema_preserves_paired_only_columns():
    fields = output_fields([{"block_width_deg": 5, "partial_r": 0.1},
                            {"block_width_deg": 5, "r_overlap": 0.2,
                             "r_no_overlap": 0.1}])
    assert fields == ["block_width_deg", "partial_r", "r_overlap", "r_no_overlap"]


def test_five_degree_baseline_comparison_requires_exact_statistics():
    actual, baseline = {}, {}
    for band in BANDS:
        actual[band], baseline[band] = {}, {}
        for category in CATEGORIES:
            actual[band][category] = {
                "n_pairs": 12, "n_spatial_blocks": 8, "partial_r": 0.2,
                "replicas_valid": 1000, "replicas_invalid": 0,
                "ci95": [-0.1, 0.4],
            }
            baseline[band][category] = {
                "n_cluster_filter_pairs": 12, "n_blocks": 8, "point_r": 0.2,
                "replicas_valid": 1000, "replicas_invalid": 0,
                "ci95": [-0.1, 0.4],
            }
        actual[band]["paired_common_sample_difference"] = {
            "n_pairs": 12, "n_spatial_blocks": 8, "r_overlap": 0.3,
            "r_no_overlap": 0.1, "delta_r_overlap_minus_no_overlap": 0.2,
            "replicas_valid": 1000, "replicas_invalid": 0, "ci95": [-0.1, 0.4],
        }
        baseline[band]["paired_common_sample_difference"] = {
            "n_cluster_filter_pairs": 12, "n_blocks": 8, "r_overlap": 0.3,
            "r_no_overlap": 0.1, "point_delta_r": 0.2,
            "replicas_valid": 1000, "replicas_invalid": 0, "ci95": [-0.1, 0.4],
        }
    assert compare_five_degree_baseline(actual, baseline)["status"] == "exact_match_within_1e-12"
    baseline["i2"]["all"]["ci95"][0] = -0.2
    with pytest.raises(ValueError, match="baseline mismatch"):
        compare_five_degree_baseline(actual, baseline)
