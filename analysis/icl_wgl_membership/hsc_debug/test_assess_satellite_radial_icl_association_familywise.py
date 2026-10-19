import numpy as np
import pytest

from assess_satellite_radial_icl_association_familywise import (
    RADIAL_BINS,
    radial_matrices,
    simultaneous_intervals,
)


def test_simultaneous_interval_is_wider_than_pointwise_scale():
    point = np.zeros(7)
    draws = np.random.default_rng(7).normal(size=(4000, 7))

    result = simultaneous_intervals(point, draws)

    assert result["valid_joint_replicates"] == 4000
    assert result["max_t_critical_value"] > 1.96
    assert all(low < 0 < high for low, high in zip(result["lower"], result["upper"]))
    assert np.allclose(np.asarray(result["upper"]), -np.asarray(result["lower"]))


def test_joint_interval_rejects_insufficient_valid_replicates():
    draws = np.ones((20, 2))
    draws[:5, 1] = np.nan

    with pytest.raises(ValueError, match="too few valid joint replicates"):
        simultaneous_intervals(np.zeros(2), draws)


def test_radial_matrices_require_identical_cluster_and_control_support():
    rows = []
    for radial in RADIAL_BINS:
        for cluster in ("a", "b", "c"):
            rows.append({
                "filter": "i", "radial_bin_kpc": radial, "cluster_id": cluster,
                "control_id": f"control-{cluster}", "ra_deg": "1", "dec_deg": "2",
                "control_ra_deg": "3", "control_dec_deg": "4",
                "icl_raw": "1", "icl_outer_subtracted": "0.5", "satellite": "2",
            })

    ids, targets, controls, values = radial_matrices(rows, "i")

    assert ids == ["a", "b", "c"]
    assert len(targets) == len(controls) == 3
    assert values["icl_raw"].shape == (3, len(RADIAL_BINS))

    rows.pop()
    with pytest.raises(ValueError, match="cluster sample is not identical"):
        radial_matrices(rows, "i")
