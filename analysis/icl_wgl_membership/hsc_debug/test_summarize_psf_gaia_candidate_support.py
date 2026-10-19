import pytest

from summarize_psf_gaia_candidate_support import quantiles


def test_quantiles_report_finite_sixteenth_median_and_eighty_fourth_percentiles():
    result = quantiles([0, 1, 2, 3, 4])

    assert result["q16"] == pytest.approx(0.64)
    assert result["median"] == pytest.approx(2.0)
    assert result["q84"] == pytest.approx(3.36)


def test_quantiles_mark_empty_samples_explicitly():
    assert quantiles([]) == {"q16": None, "median": None, "q84": None}
