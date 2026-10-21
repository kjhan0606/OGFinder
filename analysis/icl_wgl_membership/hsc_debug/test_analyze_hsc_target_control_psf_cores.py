from __future__ import annotations

import pytest

from analyze_hsc_target_control_psf_cores import _pair_rows, distribution


def _record(file: str, fwhm: float, r80: float) -> dict:
    return {
        "file": file,
        "control_id": "C1",
        "moment_fwhm_major_arcsec": fwhm,
        "moment_fwhm_minor_arcsec": fwhm * 0.9,
        "moment_axis_ratio_minor_over_major": 0.9,
        "encircled_energy_r50_arcsec": r80 * 0.5,
        "encircled_energy_r80_arcsec": r80,
        "encircled_energy_r90_arcsec": r80 * 1.2,
    }


def test_pairs_by_target_id_and_keeps_signed_difference():
    targets = {"T1": _record("target.fits", 0.8, 0.4)}
    controls = {"T1": _record("control.fits", 0.7, 0.5)}
    row, = _pair_rows(targets, controls)
    assert row["target_cluster_id"] == "T1"
    assert row["target_minus_control_moment_fwhm_major_arcsec"] == pytest.approx(0.1)
    assert row["absolute_difference_moment_fwhm_major_arcsec"] == pytest.approx(0.1)
    assert row["target_minus_control_encircled_energy_r80_arcsec"] == pytest.approx(-0.1)


def test_rejects_control_without_target_match():
    with pytest.raises(ValueError, match="no target PSF"):
        _pair_rows({}, {"T-missing": _record("control.fits", 0.7, 0.5)})


def test_distribution_is_descriptive_and_uses_expected_quantiles():
    result = distribution([1.0, 2.0, 3.0, 4.0])
    assert result["n"] == 4
    assert result["median"] == pytest.approx(2.5)
    assert result["minimum"] == 1.0
    assert result["maximum"] == 4.0
