from validate_icl_satellite_weighted_radial_matrix_familywise import (
    ESTIMATORS, FILTERS, ICL_BINS, SATELLITE_BINS, holm_adjust,
)


def test_expected_family_has_168_unique_hypotheses():
    keys = {(estimator, band, icl, satellite)
            for estimator in ESTIMATORS for band in FILTERS
            for icl in ICL_BINS for satellite in SATELLITE_BINS}
    assert len(keys) == 168


def test_holm_adjustment_reproduces_monotone_familywise_values():
    values = [0.01, 0.03, 0.2]
    adjusted = holm_adjust(values)
    assert adjusted == [0.03, 0.06, 0.2]
    assert all(b >= a for a, b in zip(values, adjusted))
