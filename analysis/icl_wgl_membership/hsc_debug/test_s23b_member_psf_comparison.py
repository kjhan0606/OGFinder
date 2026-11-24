import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location('psf_comparison', Path(__file__).with_name('compare_s23b_member_psf_injection.py'))
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def row(target='one', status='supported', response='0.5', increment='0.8'):
    return dict(target_id=target, filter='i2', member_index='7', source_row='12',
                amplitude_at_50_kpc_njy_arcsec2='1.0', baseline_delta_sb_njy_arcsec2='2.0',
                injected_status=status, recovered_increment_over_nominal_input=response,
                recovered_delta_increment_njy_arcsec2=increment)


def convolved(old, **kwargs):
    return dict(old, psf_convolution_applied='True', member_psf_sha256='a' * 64, **kwargs)


def test_pairing_uses_same_member_and_reports_both_response_and_flux():
    old = row()
    psf = convolved(old)
    psf['recovered_increment_over_nominal_input'] = '0.7'
    psf['recovered_delta_increment_njy_arcsec2'] = '0.9'
    paired = comparison.paired_rows([old], [psf])[0]
    assert paired['both_supported']
    assert paired['delta_response'] == pytest.approx(.2)
    assert paired['delta_increment'] == pytest.approx(.1)
    psf['injected_status'] = 'insufficient_coverage'
    paired = comparison.paired_rows([old], [psf])[0]
    assert not paired['both_supported'] and paired['delta_response'] is None


def test_pairing_rejects_mismatch_replay_changes_and_missing_psf():
    old = row()
    psf = convolved(old)
    wrong = dict(psf, source_row='13')
    with pytest.raises(ValueError, match='rows do not match'):
        comparison.paired_rows([old], [wrong])
    wrong = dict(psf, baseline_delta_sb_njy_arcsec2='2.1')
    with pytest.raises(ValueError, match='baseline replay'):
        comparison.paired_rows([old], [wrong])
    wrong = dict(psf, member_psf_sha256='')
    with pytest.raises(ValueError, match='PSF row'):
        comparison.paired_rows([old], [wrong])


def test_sky_cell_bootstrap_resamples_clusters_reproducibly():
    rows = [dict(target_id=str(i), delta_response=float(i)) for i in range(5)]
    frames = {str(i): {'ra_deg': 84.056247280625 + 5 * i, 'dec_deg': 0} for i in range(5)}
    first = comparison.median_block_interval(rows, 'delta_response', frames, 500, 42)
    second = comparison.median_block_interval(rows, 'delta_response', frames, 500, 42)
    assert first == second
    assert first['n_cells'] == 5 and first['median'] == 2.0
    assert first['ci_low'] <= first['median'] <= first['ci_high']
    few = comparison.median_block_interval(rows[:4], 'delta_response', frames, 500, 42)
    assert few['status'] == 'insufficient_sky_cells' and few['ci_low'] is None
