"""Behavioral checks of native profiles, area transport and paired errors."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from icl.pilot import (area_warp_metric, calibrate_flux, detection_mask,
                       native_profile, paired_summary)
from ogfkit.wlcluster import ARCSEC, Cosmology


def scale_for_kpc(kpc=1.0, z=0.425):
    return kpc / (float(Cosmology().angular_diameter(z)) * ARCSEC * 1000)


def test_calibration_uses_stored_flux_zero_point_and_explicit_extinction():
    calibration = dict(input_unit='instrumental_flux/pixel', ab_zeropoint=27.0,
                       evidence='fixture header calibration', a_i_mag=0.2,
                       extinction_evidence='fixture dust map')
    converted, factor = calibrate_flux(np.array([[4.0, -1.0]]), calibration)
    expected = 3631e9 * 10 ** (-0.4 * 27) * 10 ** (0.4 * 0.2)
    assert factor == pytest.approx(expected)
    assert converted == pytest.approx(np.array([[4, -1]]) * expected)
    with pytest.raises(ValueError, match='evidence'):
        calibrate_flux(np.ones((2, 2)), dict(calibration, evidence=None))
    with pytest.raises(ValueError, match='per pixel'):
        calibrate_flux(np.ones((2, 2)), dict(calibration, input_unit='MJy/sr'))


def test_area_transport_averages_pixels_instead_of_aliasing_a_checkerboard():
    scale = scale_for_kpc()
    image = np.tile(np.arange(65) % 2 * 2.0, (65, 1)) * scale ** 2
    mask = np.zeros(image.shape, bool)
    result = area_warp_metric(image, mask, (32, 32), scale, .425, kpc_max=16, n_side=17)
    assert result['sb'] == pytest.approx(np.ones((17, 17)))
    assert result['coverage'] == pytest.approx(np.ones((17, 17)))
    # Sampling just output centres would recover zero everywhere instead.
    assert image[32, 32] == 0


def test_area_transport_conserves_an_impulse_and_excludes_hidden_light():
    scale = scale_for_kpc()
    image = np.zeros((65, 65))
    image[32, 32] = 7.0
    mask = np.zeros(image.shape, bool)
    result = area_warp_metric(image, mask, (32, 32), scale, .425, kpc_max=16, n_side=17)
    transported = np.nansum(result['sb'] * result['coverage'] * result['output_area_arcsec2'])
    assert transported == pytest.approx(7.0)
    clean = np.full(image.shape, scale ** 2)
    clean[32, 32] = 1e12
    mask[32, 32] = True
    hidden = area_warp_metric(clean, mask, (32, 32), scale, .425, kpc_max=16, n_side=17)
    assert hidden['sb'] == pytest.approx(np.ones((17, 17)))
    assert hidden['coverage'][8, 8] == pytest.approx(.75)


def test_native_profile_counts_cropped_area_and_retains_negative_flux():
    scale = scale_for_kpc()
    image = np.full((31, 31), -2.0 * scale ** 2)
    rows = native_profile(image, np.zeros(image.shape, bool), (15, 15), scale, .425,
                          [0, 10, 30, 50], min_coverage=.5)
    assert rows[0]['sb_njy_arcsec2'] == pytest.approx(-2.0)
    assert rows[1]['sb_njy_arcsec2'] is None
    assert rows[1]['coverage'] < .5
    assert rows[2]['sb_njy_arcsec2'] is None


def profile(values):
    return [dict(sb_njy_arcsec2=v) for v in values]


def test_paired_jackknife_removes_common_region_fluctuation_in_both_inputs():
    clusters = [profile([10+x, 20+2*x]) for x in range(6)]
    randoms = [[profile([3+x, 5+2*x])] for x in range(6)]
    result = paired_summary(clusters, randoms, ['a', 'a', 'b', 'b', 'c', 'c'])
    assert result['corrected_sb'] == pytest.approx([7, 15])
    assert result['covariance'] == pytest.approx(np.zeros((2, 2)), abs=1e-28)
    assert result['n_regions'] == 3
    assert result['science_measurement'] is False


def test_covariance_is_missing_when_regions_or_random_fields_are_missing():
    clusters = [profile([4, 6]), profile([8, None])]
    result = paired_summary(clusters, [[], []], [None, None])
    assert result['unpaired_raw_sb'] == pytest.approx([6, 6])
    assert result['unpaired_n_clusters'].tolist() == [2, 1]
    assert np.isnan(result['corrected_sb']).all()
    assert np.isnan(result['covariance']).all()
    assert result['covariance_method'] == 'not estimated'


def test_all_non_bcg_detections_are_masked_without_camira_matching():
    rng = np.random.default_rng(10)
    yy, xx = np.indices((96, 96))
    image = rng.normal(0, .03, (96, 96))
    image += 50 * np.exp(-((xx-48)**2 + (yy-48)**2)/8)
    image += 30 * np.exp(-((xx-75)**2 + (yy-22)**2)/8)
    original = image.copy()
    blank = np.zeros(image.shape, bool)
    mask, info = detection_mask(image, (48, 48), .2, blank, blank)
    assert mask[22, 75]
    assert not mask[48, 48]
    assert info['n_masked_sources'] >= 1
    assert info['detection_background_subtracted_from_science'] is False
    random_mask, _ = detection_mask(image, (48, 48), .2, blank, blank, bcg_policy='none')
    assert random_mask[48, 48] and random_mask[22, 75]
    assert np.array_equal(image, original)


def test_unmatched_bcg_stops_and_mask_dilation_is_monotonic():
    rng = np.random.default_rng(12)
    yy, xx = np.indices((96, 96))
    image = rng.normal(0, .05, (96, 96))
    image += 30 * np.exp(-((xx-48)**2 + (yy-48)**2)/12)
    image += 15 * np.exp(-((xx-73)**2 + (yy-25)**2)/16)
    blank = np.zeros(image.shape, bool)
    small, _ = detection_mask(image, (48, 48), .2, blank, blank, kron_factor=2.5)
    large, _ = detection_mask(image, (48, 48), .2, blank, blank, kron_factor=5)
    assert np.all(~small | large)
    with pytest.raises(ValueError, match='BCG detection'):
        detection_mask(image, (5, 5), .2, blank, blank)


def pilot_cli():
    path = ROOT / 'analysis/icl_wgl_membership/hsc_pilot.py'
    spec = importlib.util.spec_from_file_location('hsc_pilot_cli', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fits_frame(tmp_path):
    from astropy.io import fits
    from astropy.wcs import WCS
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [24, 24]
    wcs.wcs.crval = [180, 2]
    wcs.wcs.cdelt = [-1/3600, 1/3600]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    fits.PrimaryHDU(np.ones((48, 48)), wcs.to_header()).writeto(tmp_path / 'science.fits')
    fits.PrimaryHDU(np.zeros((48, 48))).writeto(tmp_path / 'mask.fits')
    fits.PrimaryHDU(np.ones((48, 48))).writeto(tmp_path / 'variance.fits')
    mask = dict(file='mask.fits', hdu=0, evidence='fixture mask', same_grid_reviewed=True, kind='boolean')
    return dict(id='fixture', ra_deg=180, dec_deg=2, z=.425, image='science.fits', science_hdu=0,
                sky='global', product='coadd/bg', product_reviewed=True, product_evidence='fixture provenance',
                calibration=dict(input_unit='nJy/pixel', evidence='fixture calibration', a_i_mag=0,
                                 extinction_evidence='fixture extinction'), bad_mask=mask.copy(),
                bright_star_mask=mask.copy(), variance=dict(file='variance.fits', hdu=0,
                                                          evidence='fixture variance', same_grid_reviewed=True))


def test_header_inspection_preserves_product_and_zero_based_wcs(tmp_path):
    module = pilot_cli()
    row = fits_frame(tmp_path)
    result = module.frame_metadata(row, tmp_path)
    assert result['center'] == pytest.approx([23, 23], abs=1e-8)
    assert result['pixel_scale_arcsec'] == pytest.approx(1)
    assert result['product'] == 'coadd/bg'
    assert result['science_measurement'] is False
    with pytest.raises(ValueError, match='global-sky'):
        module.frame_metadata(dict(row, sky='local'), tmp_path)
    with pytest.raises(ValueError, match='science_hdu'):
        module.frame_metadata(dict(row, science_hdu=None), tmp_path)


def test_mask_bit_values_and_json_booleans_are_explicit(tmp_path):
    from astropy.io import fits
    module = pilot_cli()
    fits.PrimaryHDU(np.array([[0, 1, 4, 5]], dtype=np.int16)).writeto(tmp_path / 'bits.fits')
    mask = module.load_bool_mask(tmp_path, dict(file='bits.fits', hdu=0, kind='bits', exclude_bits=[4]))
    assert mask.tolist() == [[False, False, True, True]]
    with pytest.raises(ValueError, match='powers of two'):
        module.load_bool_mask(tmp_path, dict(file='bits.fits', hdu=0, kind='bits', exclude_bits=[3]))
    assert module.plain({'science_measurement': False}) == {'science_measurement': False}
    assert isinstance(module.plain(False), bool)


def test_measurement_refuses_the_login_node(monkeypatch):
    module = pilot_cli()
    monkeypatch.delenv('SLURM_JOB_ID', raising=False)
    with pytest.raises(SystemExit) as exc:
        module.main(['measure', '--manifest', 'missing.json', '--work', 'unused'])
    assert exc.value.code == 2


def test_primary_calibration_is_inherited_and_conflicts_stop_inspection(tmp_path):
    from astropy.io import fits
    module = pilot_cli()
    row = fits_frame(tmp_path)
    with fits.open(tmp_path / 'science.fits') as hdul:
        image = hdul[0].data.copy()
        header = hdul[0].header.copy()
    primary = fits.PrimaryHDU()
    primary.header['FLUXMAG0'] = 10 ** (27 / 2.5)
    header['INHERIT'] = True
    fits.HDUList([primary, fits.ImageHDU(image, header, name='IMAGE')]).writeto(
        tmp_path / 'inherited.fits')
    row.update(image='inherited.fits', science_hdu='IMAGE')
    row['calibration'].update(input_unit='instrumental_flux/pixel', ab_zeropoint=27)
    result = module.frame_metadata(row, tmp_path)
    assert result['center'] == pytest.approx([23, 23], abs=1e-8)
    assert result['header_fluxmag0'] == pytest.approx(10 ** (27 / 2.5))
    row['calibration']['ab_zeropoint'] = 26
    with pytest.raises(ValueError, match='conflicts with inherited FLUXMAG0'):
        module.frame_metadata(row, tmp_path)


@pytest.mark.skipif(not __import__('os').environ.get('SLURM_JOB_ID'),
                    reason='complete measurement requires an actual Slurm allocation')
def test_cli_measurement_retains_sky_and_writes_paired_profiles(tmp_path):
    import copy
    import csv
    import json
    from astropy.io import fits
    module = pilot_cli()
    row = fits_frame(tmp_path)
    rng = np.random.default_rng(123)
    yy, xx = np.indices((48, 48))
    with fits.open(tmp_path / 'science.fits') as hdul:
        header = hdul[0].header.copy()
    sky = 5 + rng.normal(0, .03, (48, 48))
    cluster = sky + 50 * np.exp(-((xx-23)**2 + (yy-23)**2)/2)
    fits.PrimaryHDU(cluster, header).writeto(tmp_path / 'cluster.fits')
    fits.PrimaryHDU(sky, header).writeto(tmp_path / 'random.fits')
    row.update(image='cluster.fits', members_radec_deg=[])
    clusters = []
    for i in range(3):
        frame = copy.deepcopy(row)
        frame.update(id='cluster%d' % i, region='region%d' % i,
                     random_matching_evidence='synthetic matched sky fixture')
        random = copy.deepcopy(row)
        random.update(id='random%d' % i, image='random.fits')
        frame['randoms'] = [random]
        clusters.append(frame)
    doc = dict(schema='hsc-bcg-diffuse-pilot-v1', band='HSC-I', rest_frame='observed',
               bcg_policy='keep', h0=70, omega_m=.3, kpc_max=60, n_side=25,
               edges_kpc=[0, 30, 60], min_coverage=.5, kron_factors=[2.5],
               checks={'injection_recovery': {'status': 'pending'}}, clusters=clusters)
    module.write_json(tmp_path / 'manifest.json', doc)
    work = tmp_path / 'measurement'
    assert module.main(['measure', '--manifest', str(tmp_path / 'manifest.json'),
                        '--work', str(work)]) == 0
    status = json.loads((work / 'pilot_status.json').read_text())
    assert status['science_measurement'] is False and status['mass_used'] is False
    summary = json.loads((work / 'kron_2p5/covariance.json').read_text())
    assert summary['n_clusters'] == [3, 3]
    assert summary['n_regions'] == 3
    assert summary['corrected_sb'][0] > 1
    assert summary['raw_sb'][1] == pytest.approx(5, abs=.02)
    assert summary['corrected_sb'][1] == pytest.approx(0, abs=.02)
    with fits.open(work / 'kron_2p5/kpc_stack.fits') as hdul:
        assert hdul[0].header['SCIENCE'] is False
        assert hdul[0].header['BUNIT'] == 'nJy/arcsec2'
        assert hdul[0].data.shape == (25, 25)
    with (work / 'kron_2p5/profile.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert float(rows[0]['CORRECTED_SB']) == pytest.approx(summary['corrected_sb'][0])
    assert float(rows[0]['PAIRED_RAW_SB']) - float(rows[0]['RANDOM_SB']) == pytest.approx(
        float(rows[0]['CORRECTED_SB']))
    assert int(rows[1]['N_PAIRED']) == 3
    # Re-running may not replace even a successful diagnostic output.
    with pytest.raises(SystemExit) as exc:
        module.main(['measure', '--manifest', str(tmp_path / 'manifest.json'),
                     '--work', str(work)])
    assert exc.value.code == 2


def test_global_probe_handles_literal_slash_and_records_all_headers(tmp_path, monkeypatch):
    import io
    import json
    import types
    import urllib.parse
    from astropy.io import fits
    path = ROOT / 'analysis/icl_wgl_membership/hsc_debug/probe_login.py'
    spec = importlib.util.spec_from_file_location('hsc_interactive_probe_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    primary = fits.PrimaryHDU()
    primary.header['FLUXMAG0'] = 10 ** (27 / 2.5)
    body = io.BytesIO()
    fits.HDUList([primary, fits.ImageHDU(np.ones((2, 2)), name='IMAGE'),
                  fits.ImageHDU(np.zeros((2, 2), dtype=np.int16), name='MASK')]).writeto(body)
    manual = b'<html>Fixture: coadd/bg</html>'
    class Opener:
        def open(self, url, timeout):
            return io.BytesIO(manual)
    monkeypatch.setattr(module.urllib.request, 'build_opener', lambda handler: Opener())
    password = 'fixture-secret/only'
    def fetch(url, user, supplied_password):
        assert user == 'fixture-user' and supplied_password == password
        query = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        assert query['type'] == ['coadd/bg'] and query['mask'] == ['on']
        return 200, 'image/fits', body.getvalue(), 40000000
    probe = types.SimpleNamespace(ROOT=str(tmp_path), BASE='https://example.test/cutout?',
                                  first_selected=lambda: ('fixture', 180, 2), probe=fetch)
    assert module.global_probe(probe, 'fixture-user', password) == 0
    logs = tmp_path / 'logs'
    status = json.loads((logs / 'cutout_global_status.json').read_text())
    assert status['product_verified'] is False and status['science_measurement'] is False
    headers = (logs / 'cutout_global_headers.txt').read_text()
    assert 'HDU 0' in headers and 'HDU 2' in headers and 'FLUXMAG0' in headers
    for filename in ('cutout_global_status.json', 'cutout_global_headers.txt', 'cutout_global_manual.html'):
        assert password not in (logs / filename).read_text()
    with pytest.raises(ValueError, match='Previous record exists'):
        module.global_probe(probe, 'fixture-user', password)
