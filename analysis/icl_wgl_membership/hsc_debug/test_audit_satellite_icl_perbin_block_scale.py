import numpy as np
import pytest

from audit_satellite_icl_perbin_block_scale import (
    BANDS,
    SELECTIONS,
    bootstrap,
    compare_baseline,
)
from audit_satellite_icl_spatial_blocks import coordinate_labels
from audit_satellite_icl_block_scale_sensitivity import connected_pair_blocks


def make_row(index):
    return {
        "cluster_id": str(index), "filter": "i2",
        "icl_delta_200_300": str(np.sin(index) + index * 0.07),
        "satellite_delta": str(np.cos(index * 0.7) + index * 0.03),
        "redshift": str(0.3 + index * 0.003 + index * index * 0.00002),
        "richness": str(8 + index + 0.1 * (index % 4)),
        "target_ra_deg": str(index * 8.0 + 0.01), "target_dec_deg": "2.0",
        "control_ra_deg": str(180 + index * 7.0 + 0.02), "control_dec_deg": "-2.0",
    }


def test_five_degree_blocks_match_archived_coordinate_definition():
    rows = [make_row(i) for i in range(12)]
    assert np.array_equal(connected_pair_blocks(rows, 5.0), coordinate_labels(rows))


def test_bootstrap_reproducibility_and_support_counts():
    rows = [make_row(i) for i in range(16)]
    first = bootstrap(rows, 5.0, 1000, 77)
    second = bootstrap(rows, 5.0, 1000, 77)
    assert first == second
    assert first["n_pairs"] == 16
    assert first["replicas_valid"] == 1000


def test_archived_baseline_comparison_rejects_mismatches():
    metrics, baseline = {}, {"results": []}
    for band in BANDS:
        metrics[band] = {}
        for selection in SELECTIONS:
            metrics[band][selection] = {
                "n_pairs": 12, "n_spatial_blocks": 8, "partial_r": 0.2,
                "replicas_valid": 1000, "replicas_invalid": 0,
                "ci95": [-0.1, 0.4],
            }
            baseline["results"].append({
                "filter": band, "selection": selection,
                "partial_spatial_block_bootstrap": {
                    "n_pairs": 12, "n_blocks": 8, "partial_r_point": 0.2,
                    "replicas_valid": 1000, "replicas_invalid": 0,
                    "ci95_partial_r": [-0.1, 0.4],
                },
            })
    assert compare_baseline(metrics, baseline)["status"] == "exact_match_within_1e-12"
    metrics["i2"]["all_supported_pairs"]["ci95"][1] = 0.5
    with pytest.raises(ValueError, match="interval mismatch"):
        compare_baseline(metrics, baseline)
