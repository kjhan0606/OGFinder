import numpy as np
import pytest

from summarize_hsc_psf_core_variation import measure_core


def test_measure_core_uses_positive_flux_for_descriptive_moments_and_keeps_input():
    axis = np.arange(-10, 11, dtype=float)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    positive = np.exp(-0.5 * ((xx / 2.0) ** 2 + (yy / 1.5) ** 2))
    positive /= positive.sum()
    psf = positive.copy()
    psf[0, 0] -= 1e-4
    psf[10, 10] += 1e-4
    original = psf.copy()

    result = measure_core(psf, pixel_scale_arcsec=1.0)

    assert np.array_equal(psf, original)
    assert result["negative_pixel_count"] == 1
    assert result["signed_sum"] == pytest.approx(1.0)
    assert result["centroid_xy_pixels"] == pytest.approx([10.0, 10.0], abs=0.02)
    assert result["moment_fwhm_major_arcsec"] > result["moment_fwhm_minor_arcsec"]
    assert 0 < result["moment_axis_ratio_minor_over_major"] < 1
    assert (result["encircled_energy_r50_arcsec"]
            < result["encircled_energy_r80_arcsec"]
            < result["encircled_energy_r90_arcsec"])


@pytest.mark.parametrize("psf", [
    np.ones((4, 4)),
    np.full((3, 3), np.nan),
    np.zeros((3, 3)),
])
def test_measure_core_rejects_invalid_psf_arrays(psf):
    with pytest.raises(ValueError):
        measure_core(psf)
