import numpy as np
import pytest

from audit_icl_psf_region_fixed_effects import partial_correlation


def test_region_effect_removes_between_region_correlation():
    rows = []
    for region in range(3):
        for offset in (-1.0, 0.0, 1.0, 2.0):
            rows.append({
                "icl_delta": 10.0 * region + offset,
                "psf_fwhm_major_arcsec": 2.0 * region - offset,
                "target_region": region,
            })

    raw = partial_correlation(rows, (), region_effect=False)
    adjusted = partial_correlation(rows, (), region_effect=True)

    assert raw["partial_r"] > 0.5
    assert adjusted["partial_r"] == pytest.approx(-1.0)
    assert adjusted["design_rank"] == 3
    assert adjusted["residual_degrees_of_freedom"] == 9


def test_region_adjustment_requires_residual_degrees_of_freedom():
    rows = [{"icl_delta": float(i), "psf_fwhm_major_arcsec": float(i * i),
             "target_region": i} for i in range(5)]
    with pytest.raises(ValueError, match="degrees of freedom"):
        partial_correlation(rows, (), region_effect=True)


def test_region_model_requires_a_single_label_for_every_row():
    rows = [{"icl_delta": float(i), "psf_fwhm_major_arcsec": float(i * i),
             "target_region": region} for i, region in enumerate([0, 0, 1, 1, 1, 0])]
    result = partial_correlation(rows, (), region_effect=True)
    assert result["n_frames"] == 6
    assert result["residual_degrees_of_freedom"] == 4
