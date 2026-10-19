import numpy as np
import pytest

from run_wl_cluster_spatial_jackknife import delete_block_covariance, equal_area_block_ids


def test_equal_area_blocks_wrap_ra_and_bound_poles():
    blocks = equal_area_block_ids(
        [0.1, 359.9, 45.0, 45.0], [-90.0, -90.0, 0.0, 90.0],
        cell_width_deg=5.0,
    )
    assert blocks[0] == blocks[1]
    assert len(set(blocks)) >= 3


def test_equal_area_blocks_validate_positions_and_width():
    with pytest.raises(ValueError, match="sky coordinates"):
        equal_area_block_ids([360.0], [0.0])
    with pytest.raises(ValueError, match="cell_width_deg"):
        equal_area_block_ids([0.0], [0.0], cell_width_deg=0.0)


def test_delete_block_covariance_matches_standard_formula():
    replicates = np.array([[1.0, 0.0], [2.0, 1.0], [3.0, -1.0]])
    covariance = delete_block_covariance(replicates)
    expected = (2.0 / 3.0) * np.cov(replicates, rowvar=False, ddof=0) * len(replicates)
    np.testing.assert_allclose(covariance, expected)


def test_delete_block_covariance_requires_finite_replicates():
    with pytest.raises(ValueError, match="at least 3"):
        delete_block_covariance([[1.0], [2.0]])
    with pytest.raises(ValueError, match="finite"):
        delete_block_covariance([[1.0], [np.nan], [3.0]])
