#!/usr/bin/env python3
"""Rehash all controls and independently replay both sides of the sky comparison.

Does not call tile_diagnostics: weighted normal equations verify the plane fit;
whole-image pair arrays with explicit tile boundaries verify the lag products.
No download, sky correction, or scientific acceptance is performed.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS


def read(path):
    return json.loads(path.read_text())


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def unique(rows, fields):
    result = {}
    for row in rows:
        key = tuple(str(row[f]) for f in fields)
        if key in result:
            raise ValueError('Duplicate row: '+str(key))
        result[key] = row
    return result


def fit_plane(tiles):
    design = np.array([[1., r['x_kpc']/1000, r['y_kpc']/1000] for r in tiles])
    y = np.array([r['raw_mean'] for r in tiles])
    sigma = np.array([r['nominal_mean_sigma'] for r in tiles])
    weighted = design/sigma[:, None]
    normal = weighted.T@weighted
    if len(tiles) < 10 or np.linalg.matrix_rank(normal) != 3:
        raise ValueError('Not enough independent tile positions for a plane')
    coef = np.linalg.solve(normal, weighted.T@(y/sigma))
    residual = y-design@coef
    mad = 1.4826*np.median(np.abs(residual-np.median(residual)))
    return dict(plane_constant=coef[0], plane_x_per_mpc=coef[1], plane_y_per_mpc=coef[2],
                median_tile_mean=np.median(y), median_nominal_mean_sigma=np.median(sigma),
                plane_residual_mad=mad, plane_residual_mad_over_nominal=mad/np.median(sigma)), residual


def lag_statistic(normalized, accepted, lag, axis, tile_size=64):
    # An image-wide pair is legal only inside ONE accepted tile, with both ends valid.
    first = [slice(None), slice(None)]; second = [slice(None), slice(None)]
    first[axis] = slice(lag, None); second[axis] = slice(None, -lag)
    first, second = tuple(first), tuple(second)
    coords = np.arange(normalized.shape[axis]-lag)
    same_tile = (coords//tile_size) == ((coords+lag)//tile_size)
    boundary = same_tile[:, None] if axis == 0 else same_tile[None, :]
    keep = accepted[first] & accepted[second] & boundary
    x = normalized[first][keep]; y = normalized[second][keep]
    if not len(x):
        return dict(n_pairs=0, variance_normalized_covariance=None, correlation=None)
    # Center pair arrays explicitly, independently of the production sum-products formula.
    x = x-x.mean(); y = y-y.mean()
    covariance = np.mean(x*y); vx = np.mean(x*x); vy = np.mean(y*y)
    return dict(n_pairs=len(x), variance_normalized_covariance=covariance,
                correlation=covariance/np.sqrt(vx*vy) if vx > 0 and vy > 0 else None)


def near(actual, expected, name, errors, rtol=2e-10, atol=2e-11):
    if actual is None or expected in (None, ''):
        if actual is not None or expected not in (None, ''):
            raise AssertionError('Missing value mismatch: '+name)
        return
    actual, expected = float(actual), float(expected)
    if not np.isfinite(actual) or not np.isfinite(expected) or not np.isclose(actual, expected, rtol=rtol, atol=atol):
        raise AssertionError(f'{name}: recomputed={actual}, saved={expected}')
    errors[name] = max(errors.get(name, 0.), abs(actual-expected))


def replay(image, variance, excluded, region, center, kpp, saved_tiles, saved_stat, saved_lags, errors):
    valid = region & ~excluded & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    normalized = np.zeros(image.shape, dtype=float); accepted = np.zeros(image.shape, dtype=bool)
    rows = []; positions = set()
    for y in range(0, image.shape[0], 64):
        for x in range(0, image.shape[1], 64):
            block = np.s_[y:y+64, x:x+64]; mask = valid[block]
            n = int(mask.sum()); area = int(region[block].sum())
            selected = area >= 128 and n >= 128 and 2*n >= area
            key = (str(x), str(y))
            if selected != (key in saved_tiles):
                raise AssertionError('Accepted tile set changed at '+str(key))
            if not selected:
                continue
            positions.add(key); saved = saved_tiles[key]
            pixels = image[block].astype(float)[mask]; var = variance[block].astype(float)[mask]
            gy, gx = np.nonzero(mask)
            row = dict(n_valid=n, n_region=area, valid_fraction=n/area,
                       raw_mean=pixels.mean(), raw_median=np.median(pixels),
                       nominal_mean_sigma=np.sqrt(var.sum())/n,
                       x_kpc=(x+gx.mean()-center[0])*kpp, y_kpc=(y+gy.mean()-center[1])*kpp)
            for field, value in row.items():
                near(value, saved[field], field, errors)
            normalized[block][mask] = (pixels-row['raw_median'])/np.sqrt(var)
            accepted[block] = mask
            row['saved_plane_residual'] = saved['plane_residual']; rows.append(row)
    if positions != set(saved_tiles):
        raise AssertionError('Saved tiles include an unknown grid location')
    for key, value in dict(n_tiles=len(rows), n_outer_valid=int(valid.sum()),
                           n_outer_region=int(region.sum()), outer_valid_fraction=valid.sum()/region.sum()).items():
        near(value, saved_stat[key], key, errors)
    fit, residual = fit_plane(rows)
    for key, value in fit.items():
        near(value, saved_stat[key], key, errors)
    for row, value in zip(rows, residual):
        near(value, row['saved_plane_residual'], 'tile_plane_residual', errors)
    expected_keys = {(str(lag), axis) for lag in (1, 2, 4, 8) for axis in ('x', 'y')}
    if set(saved_lags) != expected_keys:
        raise AssertionError('Missing or extra lag products')
    for lag in (1, 2, 4, 8):
        for axis, name in [(0, 'y'), (1, 'x')]:
            result = lag_statistic(normalized, accepted, lag, axis)
            for key, value in result.items():
                near(value, saved_lags[(str(lag), name)][key], key, errors)
    return len(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False); start = time.time(); hashes = {}
    script_digest = sha(Path(__file__))
    sample = a.sample.resolve()
    quality_dir = sample/'control_quality_2026-10-14_01'
    control_dir = sample/'control_outer_2026-10-14_01'
    cluster_dir = sample/'outer_sky_diagnostics_2026-10-14_02'
    # The preceding audit pins every result table from the original run.
    prior_path = sample/'control_image_audit_2026-10-14_01/audit.json'
    prior = read(prior_path)
    if prior['passed'] is not True:
        raise ValueError('Prior audit did not pass')
    for path, digest in prior['source_sha256'].items():
        actual = sha(path)
        if actual != digest:
            raise ValueError('Changed recorded artifact: '+path)
        hashes[path] = actual
    hashes[str(prior_path)] = sha(prior_path)
    manifest = read(sample/'control_download_01/download_manifest.json')
    frames = unique(manifest['clusters'], ['id'])
    download = unique(read(sample/'control_download_01/download_status.json')['completed'], ['id'])
    qa = unique(read(quality_dir/'quality_frames.json'), ['id'])
    if set(frames) != set(download) or set(qa) != set(frames):
        raise ValueError('Manifest, completed downloads, and QA IDs disagree')
    total_bytes = 0
    for i, (key, frame) in enumerate(frames.items(), 1):
        path = Path(qa[key]['path']); digest = sha(path)
        if digest != download[key]['sha256'] or digest != qa[key]['sha256_rechecked']:
            raise ValueError('Control image hash mismatch: '+str(key))
        size = path.stat().st_size
        if size != download[key]['bytes'] or size != qa[key]['bytes']:
            raise ValueError('Control image size mismatch: '+str(key))
        hashes[str(path)] = digest; total_bytes += size
        if i % 100 == 0 or i == len(frames):
            print(f'Full control hashes {i}/{len(frames)}', flush=True)
    cosmo = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'], Tcmb0=0)
    target_path = sample/'validation_sample_2026-10-13_01/validation_manifest.json'
    targets = unique(read(target_path)['clusters'], ['id'])
    parent = read(cluster_dir/'summary.json')
    for path, digest in parent['source_sha256'].items():
        if path not in hashes:
            hashes[path] = sha(path)
        if hashes[path] != digest:
            raise ValueError('Changed cluster dependency: '+path)
    records = []; errors = {}; tile_count = 0
    for role, directory, id_field in [('control', control_dir, 'control_id'), ('cluster', cluster_dir, 'cluster_id')]:
        stats = unique(table(directory/'frame_stats.tsv'), [id_field, 'kron_factor'])
        tile_groups = defaultdict(list); lag_groups = defaultdict(list)
        for row in table(directory/'tile_stats.tsv'):
            tile_groups[(row[id_field], row['kron_factor'])].append(row)
        for row in table(directory/'lag_products.tsv'):
            lag_groups[(row[id_field], row['kron_factor'])].append(row)
        ids = sorted({key[0] for key in stats})
        if len(ids) != 50 or len(stats) != 150 or set(tile_groups) != set(stats) or set(lag_groups) != set(stats):
            raise ValueError('Unexpected diagnostic frame/method coverage')
        for index, ident in enumerate(ids, 1):
            frame = frames[(ident,)] if role == 'control' else targets[(ident,)]
            path = Path(qa[(ident,)]['path']) if role == 'control' else Path(frame['image'])
            mask_path = (directory/ident/'masks.fits.gz' if role == 'control' else
                         sample/'bcg_masks_2026-10-13_01'/('cluster_'+ident)/'proposed_masks.fits.gz')
            with fits.open(path, memmap=False) as h, fits.open(mask_path, memmap=False) as mh:
                image = h[frame['science_hdu']].data; variance = h[frame['variance']['hdu']].data
                wcs = WCS(h[frame['science_hdu']].header).celestial
                center = np.array(wcs.all_world2pix(frame['ra_deg'], frame['dec_deg'], 0))
                scale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
                kpp = cosmo.kpc_proper_per_arcmin(frame['z']).value*scale/60
                yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
                radius = np.hypot(xx-center[0], yy-center[1])*kpp
                region = (radius >= 600) & (radius < 900)
                for factor in ('2.5', '3.5', '5.0'):
                    key = (ident, factor)
                    count = replay(image, variance, mh['MASK_K'+factor.replace('.', 'P')].data.astype(bool),
                                   region, center, kpp, unique(tile_groups[key], ['tile_x', 'tile_y']),
                                   stats[key], unique(lag_groups[key], ['lag_pixels', 'axis']), errors)
                    tile_count += count
                    records.append(dict(role=role, id=ident, kron_factor=factor, tiles_replayed=count, lags_replayed=8))
            if index % 10 == 0:
                print(f'{role} full statistical replay {index}/{len(ids)}', flush=True)
        for name in ('summary.json', 'frame_stats.tsv', 'tile_stats.tsv', 'lag_products.tsv'):
            path = directory/name; hashes[str(path)] = sha(path)
    result = dict(passed=True, n_control_files_rehashed=len(frames), bytes_rehashed=total_bytes,
                  n_method_fits_replayed=len(records), n_tiles_replayed=tile_count,
                  n_lag_products_replayed=len(records)*8, max_absolute_differences=errors,
                  methods=records, source_sha256=hashes, script_sha256=script_digest,
                  elapsed_seconds=time.time()-start, host=platform.node(),
                  science_measurement=False, controls_ready_for_science=False, sky_correction_approved=False,
                  limitations=['Reproduces the existing diagnostic definition; no injection/recovery or sky covariance calibration.',
                               'Finite floating point agreement is not observational acceptance or sample completeness.'])
    with (a.work/'verification.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'methods', 'limitations')}))


if __name__ == '__main__':
    main()
