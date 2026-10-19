import numpy as np
import pytest

from audit_icl_outer_annulus_baseline import block_bootstrap, weighted_annulus


def _bin(r0, r1, value, valid, geometric):
    return {
        "r_min_kpc": r0,
        "r_max_kpc": r1,
        "delta_mean_njy_arcsec2": value,
        "n_valid_pixels": valid,
        "geometric_pixels": geometric,
    }


def test_weighted_annulus_uses_valid_pixel_counts():
    profile = [_bin(400, 425, 2.0, 1, 10), _bin(425, 450, 8.0, 3, 10)]

    mean, n_pixels, coverage = weighted_annulus(profile, 400, 450)

    assert mean == pytest.approx(6.5)
    assert n_pixels == 4
    assert coverage == pytest.approx(0.2)


def test_weighted_annulus_retains_no_support_as_missing():
    profile = [_bin(400, 425, None, 0, 10)]

    assert weighted_annulus(profile, 400, 425) == (None, 0, 0.0)


def test_connected_block_bootstrap_keeps_paired_outer_shift():
    records = [
        {"ra_deg": 130.0, "dec_deg": 1.0, "control_ra_deg": 130.2,
         "control_dec_deg": 1.1, "inner_delta": 11.0, "outer_delta": 1.0,
         "inner_coverage": 0.2, "outer_coverage": 0.1},
        {"ra_deg": 131.0, "dec_deg": 1.0, "control_ra_deg": 130.8,
         "control_dec_deg": 1.1, "inner_delta": 21.0, "outer_delta": 1.0,
         "inner_coverage": 0.3, "outer_coverage": 0.1},
        {"ra_deg": 220.0, "dec_deg": 2.0, "control_ra_deg": 220.2,
         "control_dec_deg": 2.1, "inner_delta": 31.0, "outer_delta": 1.0,
         "inner_coverage": 0.4, "outer_coverage": 0.1},
        {"ra_deg": 221.0, "dec_deg": 2.0, "control_ra_deg": 220.8,
         "control_dec_deg": 2.1, "inner_delta": 41.0, "outer_delta": 1.0,
         "inner_coverage": 0.5, "outer_coverage": 0.1},
    ]

    result = block_bootstrap(records, replicas=200, seed=11)

    assert result["n_connected_blocks"] == 2
    assert result["raw_mean_njy_arcsec2"] == pytest.approx(26.0)
    assert result["outer_400_500_mean_njy_arcsec2"] == pytest.approx(1.0)
    assert result["outer_subtracted_mean_njy_arcsec2"] == pytest.approx(25.0)
    assert result["paired_outer_subtraction_shift_njy_arcsec2"] == pytest.approx(-1.0)
    assert result["median_inner_coverage"] == pytest.approx(0.35)
    assert np.isfinite(result["outer_subtracted_block_bootstrap_95_interval"]).all()
