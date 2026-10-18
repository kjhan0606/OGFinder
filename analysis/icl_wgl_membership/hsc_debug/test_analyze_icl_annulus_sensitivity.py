import numpy as np

from analyze_icl_annulus_sensitivity import pair_is_eligible


def test_pair_eligibility_requires_matching_filter_and_both_area_cuts():
    target = {"filter": "i"}
    control = {"header_filter": "i"}
    target_coverage = np.full((3, 6), 0.8)
    control_coverage = np.full((3, 6), 0.6)
    target_raw = np.full((3, 6), 10.0)
    control_raw = np.full((3, 6), 4.0)
    assert pair_is_eligible(target, target_coverage, target_raw, control,
                            control_coverage, control_raw, 0, 0.6)
    assert not pair_is_eligible(target, target_coverage, target_raw, control,
                                control_coverage, control_raw, 0, 0.61)
    assert not pair_is_eligible(target, target_coverage, target_raw,
                                {"header_filter": "i2"}, control_coverage,
                                control_raw, 0, 0.25)


def test_pair_eligibility_rejects_nonfinite_flux():
    target = {"filter": "i2"}
    control = {"header_filter": "i2"}
    coverage = np.ones((3, 6))
    target_raw = np.ones((3, 6))
    target_raw[2, 3] = np.nan
    assert not pair_is_eligible(target, coverage, target_raw, control,
                                coverage, np.ones((3, 6)), 2, 0.25)
