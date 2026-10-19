import pytest

from audit_icl_coadd_noise_metadata import summarize_keyword_presence


def test_keyword_summary_distinguishes_missing_gain_from_variance_and_calibration():
    records = [
        {"headers": [{"FLUXMAG0": 1.0}, {}, {"BUNIT": "adu2"}]},
        {"headers": [{"FLUXMAG0": 1.0, "GAIN": 2.0}, {"EXPTIME": 300.0}]},
    ]
    summary = summarize_keyword_presence(records)
    assert summary["FLUXMAG0"]["frames_with_keyword"] == 2
    assert summary["GAIN"]["frames_with_keyword"] == 1
    assert summary["EXPTIME"]["frames_with_keyword"] == 1
    assert summary["BUNIT"]["frames_with_keyword"] == 1


def test_keyword_summary_requires_frames():
    with pytest.raises(ValueError, match="at least one"):
        summarize_keyword_presence([])
