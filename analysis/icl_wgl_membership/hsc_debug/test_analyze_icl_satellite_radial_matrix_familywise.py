import numpy as np

from analyze_icl_satellite_radial_matrix_familywise import (
    WEBB_WEIGHTS, bootstrap_t, holm_adjust, prepare_group,
)


def _rows(n=24):
    rows = []
    for i in range(n):
        z = 0.2 + 0.015 * i
        richness = 10 + (i * 7) % 19
        x = 2.0 + 0.3 * i + 0.04 * (i % 5)
        y = 5.0 + 0.45 * i + 0.02 * ((i * 3) % 7)
        rows.append({"icl_delta_200_300": x, "satellite_delta_50_75": y,
                     "z": z, "richness": richness})
    return rows


def test_group_preparation_and_restricted_wild_bootstrap_are_finite():
    rows = _rows()
    labels = [f"block_{i // 3}" for i in range(len(rows))]
    block_by_row = {id(row): label for row, label in zip(rows, labels)}

    group, summary = prepare_group(rows, block_by_row, min_blocks=8)

    assert summary["status"] == "eligible"
    assert group["n_clusters"] == 24
    assert group["n_spatial_blocks"] == 8
    weights = {label: float(WEBB_WEIGHTS[i % len(WEBB_WEIGHTS)])
               for i, label in enumerate(sorted(set(labels)))}
    assert np.isfinite(bootstrap_t(group, weights))


def test_group_is_excluded_below_spatial_block_rule():
    rows = _rows(18)
    labels = [f"block_{i // 3}" for i in range(len(rows))]
    block_by_row = {id(row): label for row, label in zip(rows, labels)}
    group, summary = prepare_group(rows, block_by_row, min_blocks=8)
    assert group is None
    assert summary["status"] == "insufficient_spatial_blocks"
    assert summary["n_spatial_blocks"] == 6


def test_holm_adjustment_is_monotone_and_bounded():
    adjusted = holm_adjust(np.asarray([0.01, 0.03, 0.2]))
    assert np.allclose(adjusted, [0.03, 0.06, 0.2])
    assert np.all(adjusted >= np.asarray([0.01, 0.03, 0.2]))
    assert np.all((adjusted >= 0) & (adjusted <= 1))
