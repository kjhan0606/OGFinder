#!/usr/bin/env python3
"""Independent raw-pixel replay of control QA, nested masks, and paired tiles."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import platform

from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from download_pilot import file_sha256 as sha
import prepare_detection_psf as prep


def read(path):
    return json.loads(path.read_text())


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('sample', 'quality', 'comparison', 'work'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False)
    qs, cs = read(a.quality/'summary.json'), read(a.comparison/'summary.json')
    assert qs['status'] == cs['status'] == 'complete'
    assert qs['science_measurement'] is cs['science_measurement'] is False
    assert qs['controls_ready_for_science'] is cs['sky_correction_approved'] is False
    hashes = dict(qs['source_sha256'], **cs['source_sha256'], **cs['output_masks_sha256'])
    assert qs['script_sha256'] == cs['script_sha256']
    for path, digest in cs['script_sha256'].items():
        assert sha(path) == digest, path
    for path, digest in hashes.items():
        assert sha(path) == digest, path
    n_verified = len(hashes)
    for directory, names in [(a.quality, ['summary.json', 'quality_frames.json', 'radial_coverage.tsv']),
                              (a.comparison, ['summary.json', 'frame_stats.tsv', 'paired_stats.tsv',
                                              'tile_stats.tsv', 'lag_products.tsv'])]:
        for name in names:
            path = directory/name; hashes[str(path)] = sha(path)
    manifest = read(a.sample/'control_download_01/download_manifest.json')
    frames = {f['id']: f for f in manifest['clusters']}
    downloads = {r['id']: r for r in read(a.sample/'control_download_01/download_status.json')['completed']}
    quality = {r['id']: r for r in read(a.quality/'quality_frames.json')}
    assert set(frames) == set(quality) == set(downloads)
    assert len(quality) == qs['n_reported'] == qs['n_full_file_hashes_matched']
    assert sum(r['bytes'] for r in quality.values()) == qs['bytes_rechecked']
    for ident, row in quality.items():
        assert row['sha256_rechecked'] == downloads[ident]['sha256']
        assert row['target_cluster_id'] == frames[ident]['target_cluster_id']
        assert row['frame_role'] == 'matched_control' and row['bcg_policy'] == 'none'
        assert not row['science_measurement']
    radial = np.array([r['coverage_by_25kpc_annulus_0_500'] for r in quality.values()])
    assert int(np.all(radial >= .5, axis=1).sum()) == qs['n_ge_50pct_in_every_0_500kpc_bin']
    for i, row in enumerate(table(a.quality/'radial_coverage.tsv')):
        assert int(row['n_controls_ge_50pct']) == int((radial[:, i] >= .5).sum())
    stats = table(a.comparison/'frame_stats.tsv'); paired = table(a.comparison/'paired_stats.tsv')
    cluster_stats = {(r['cluster_id'], r['kron_factor']): r for r in table(
        a.sample/'outer_sky_diagnostics_2026-10-14_02/frame_stats.tsv')}
    cluster_downloads = {r['id']: r for r in read(a.sample/'download_01/download_status.json')['completed']}
    control_stats = {(r['control_id'], r['kron_factor']): r for r in stats}
    tiles = table(a.comparison/'tile_stats.tsv'); lags = table(a.comparison/'lag_products.tsv')
    assert len(stats) == cs['n_method_rows'] == 3*cs['n_control_frames']
    assert len(paired) == len(stats) and len(tiles) == cs['n_tiles'] and len(lags) == 8*len(stats)
    counts = Counter((r['control_id'], r['kron_factor']) for r in tiles)
    groups = defaultdict(list); tile_groups = defaultdict(list)
    for row in stats:
        assert int(row['n_tiles']) == counts[(row['control_id'], row['kron_factor'])]
        groups[row['control_id']].append(row)
    for row in tiles:
        tile_groups[(row['control_id'], row['kron_factor'])].append(row)
    cosmo = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'], Tcmb0=0)
    replays = []; raw_checks = 0; mask_checks = 0; wcs_checks = 0
    bad_value = sum(1 << bit for bit in (0, 1, 2, 3, 4, 7, 8, 9, 10, 12, 13, 14, 15, 16))
    for ident, methods in groups.items():
        frame = frames[ident]; report = quality[ident]
        source = Path(report['path']); masks = a.comparison/ident/'masks.fits.gz'
        with fits.open(source, memmap=False) as h, fits.open(masks, memmap=False) as mh:
            image = h[frame['science_hdu']].data; variance = h[frame['variance']['hdu']].data
            bits = h[frame['bad_mask']['hdu']].data
            valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ((bits & bad_value) == 0)
            np.testing.assert_allclose(valid.sum()/valid.size, report['quality_unmasked_fraction'], atol=1e-15)
            np.testing.assert_array_equal(~valid, mh['BASE_MASK'].data.astype(bool))
            wcs = WCS(h[frame['science_hdu']].header).celestial
            cx, cy = wcs.all_world2pix(frame['ra_deg'], frame['dec_deg'], 0)
            # Independent coordinate formula; same declared cosmology and nominal pixel scale.
            scale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
            kpp = cosmo.kpc_proper_per_arcmin(frame['z']).value*scale/60
            yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
            radius = np.hypot(xx-cx, yy-cy)*kpp
            for i in range(20):
                annulus = (radius >= i*25) & (radius < (i+1)*25)
                coverage = int((annulus & valid).sum())/int(annulus.sum())
                np.testing.assert_allclose(coverage, report['coverage_by_25kpc_annulus_0_500'][i], atol=1e-12)
                raw_checks += 1
            outer = (radius >= 600) & (radius < 900)
            last = ~valid
            points = np.array([[0, 0], [image.shape[1]-1, 0], [0, image.shape[0]-1], [cx, cy]])
            world = wcs.all_pix2world(points, 0)
            for method in sorted(methods, key=lambda r: float(r['kron_factor'])):
                factor = method['kron_factor']; hdu = mh['MASK_K'+factor.replace('.', 'P')]
                mask = hdu.data.astype(bool)
                assert not np.any(last & ~mask) and np.all(mask[mh['SEGMENTATION'].data > 0])
                assert hdu.header['BCGPOL'] == 'none' and not hdu.header['BGSUB']
                np.testing.assert_allclose(WCS(hdu.header).celestial.all_world2pix(world, 0), points, atol=.01, rtol=0)
                assert int((outer & ~mask).sum()) == int(method['n_outer_valid'])
                assert int(outer.sum()) == int(method['n_outer_region'])
                last = mask; mask_checks += 1; wcs_checks += len(points)
                selected = tile_groups[(ident, factor)]
                for row in ([selected[0], selected[len(selected)//2], selected[-1]] if selected else []):
                    x, y = int(row['tile_x']), int(row['tile_y']); sl = np.s_[y:y+64, x:x+64]
                    good = outer[sl] & ~mask[sl] & valid[sl]
                    pixels = image[sl].astype(float)[good]; var = variance[sl].astype(float)[good]
                    assert len(pixels) == int(row['n_valid'])
                    assert int(outer[sl].sum()) == int(row['n_region'])
                    delta = abs(pixels.mean()-float(row['raw_mean']))
                    np.testing.assert_allclose(pixels.mean(), float(row['raw_mean']), atol=1e-12)
                    np.testing.assert_allclose(np.median(pixels), float(row['raw_median']), atol=1e-12)
                    np.testing.assert_allclose(np.sqrt(var.sum())/len(var), float(row['nominal_mean_sigma']), rtol=1e-12)
                    replays.append(dict(control_id=ident, kron_factor=factor, tile_x=x, tile_y=y, mean_difference=delta))
    for row in paired:
        if row['control_status'] == row['cluster_status'] == 'diagnostic_only':
            for role, source, download in (
                    ('cluster', cluster_stats[(row['target_cluster_id'], row['kron_factor'])], cluster_downloads[row['target_cluster_id']]),
                    ('control', control_stats[(row['control_id'], row['kron_factor'])], downloads[row['control_id']])):
                info = download['review']
                assert info['input_unit'] == 'instrumental_flux/pixel'
                converted = float(source['median_tile_mean'])*3.631e12/info['fluxmag0']/info['pixel_scale_arcsec']**2
                np.testing.assert_allclose(converted, float(row[role+'_median_tile_mean_njy_arcsec2']), rtol=1e-12, atol=1e-12)
            expected = float(row['cluster_median_tile_mean_njy_arcsec2'])-float(row['control_median_tile_mean_njy_arcsec2'])
            np.testing.assert_allclose(float(row['cluster_minus_control_median_tile_mean_njy_arcsec2']), expected, atol=1e-12)
    result = dict(passed=True, input_hashes_verified=n_verified, script_hashes_verified=len(cs['script_sha256']),
                  n_quality_records=len(quality), n_pair_conversions_checked=len(paired),
                  n_independent_radial_replays=raw_checks, n_full_mask_replays=mask_checks,
                  n_mask_wcs_points=wcs_checks, n_raw_tile_replays=len(replays),
                  maximum_tile_mean_difference=max(r['mean_difference'] for r in replays),
                  source_sha256=hashes, script_sha256=sha(Path(__file__)),
                  science_measurement=False, controls_ready_for_science=False, sky_correction_approved=False)
    prep.write_json(a.work/'audit.json', result)
    prep.write_json(a.work/'tile_replays.json', replays)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
