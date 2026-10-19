import numpy as np
import pytest

from audit_extended_psf_core_vs_wide import quantiles


def test_quantiles_are_descriptive_and_mark_zero_overlap():
    result = quantiles([-2.0, -1.0, 0.0, 1.0, 4.0])
    assert result["median"] == pytest.approx(0.0)
    assert result["q16"] == pytest.approx(-1.36)
    assert result["q84"] == pytest.approx(2.08)
    assert result["q16_q84_contains_zero"] is True


def test_quantiles_reject_empty_or_non_finite_input():
    with pytest.raises(ValueError):
        quantiles([])
    with pytest.raises(ValueError):
        quantiles([1.0, np.nan])
