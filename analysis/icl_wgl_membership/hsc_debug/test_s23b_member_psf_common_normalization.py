import math
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from analyze_s23b_member_psf_common_normalization import components


def sample(nominal, flux):
    return {'nominal_unmasked_50_75_mean_njy_arcsec2': str(nominal),
            'recovered_delta_increment_njy_arcsec2': str(flux),
            'recovered_increment_over_nominal_input': str(flux / nominal)}


def test_common_normalization_decomposes_ratio_change():
    result = components(sample(2, 1), sample(4, 3))
    assert result['nominal_fractional_change'] == pytest.approx(1)
    assert result['common_denominator_response_change'] == pytest.approx(1)
    assert result['denominator_contribution'] == pytest.approx(-.75)
    assert result['original_ratio_change'] == pytest.approx(.25)
    assert result['common_denominator_response_change'] + result['denominator_contribution'] == pytest.approx(result['original_ratio_change'])


def test_invalid_nominal_and_inconsistent_stored_ratio_are_rejected():
    good = sample(2, 1)
    bad = dict(good, nominal_unmasked_50_75_mean_njy_arcsec2='0')
    with pytest.raises(ValueError, match='nonpositive'):
        components(bad, good)
    bad = dict(good, recovered_increment_over_nominal_input='2')
    with pytest.raises(ValueError, match='stored response'):
        components(good, bad)
    bad = dict(good, recovered_delta_increment_njy_arcsec2=str(math.nan))
    with pytest.raises(ValueError, match='nonfinite'):
        components(good, bad)
