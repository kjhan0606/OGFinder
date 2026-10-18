#!/usr/bin/env python3
"""Prepare traceable SFD estimates and HSC catalogue requests, without credentials.

No raw science image or existing manifest is changed. Nominal i/i2 extinction
estimates are not certified as PDR3 effective-passband catalogue corrections.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import re
import urllib.request

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates

MAP_BASE = 'https://portal.nersc.gov/project/cosmo/data/dust/v0_1/maps/'
COEFFICIENT_SOURCE = 'https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7_e/colorterms.html'
COEFFICIENTS = {'i': 1.641, 'i2': 1.617}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def acquire_map(cache, pole):
    name = f'SFD_dust_4096_{pole}.fits'
    path = cache / name
    url = MAP_BASE + name
    if not path.exists():
        temporary = path.with_suffix('.fits.part')
        with urllib.request.urlopen(url, timeout=90) as response, temporary.open('xb') as out:
            for block in iter(lambda: response.read(1024 * 1024), b''):
                out.write(block)
        with fits.open(temporary) as hdul:
            if hdul[0].shape != (4096, 4096):
                raise ValueError('Unexpected SFD map dimensions')
        temporary.rename(path)
    return path, dict(url=url, file=str(path), bytes=path.stat().st_size, sha256=digest(path))


def query_sfd(coords, maps):
    """Unscaled SFD E(B-V), bilinear interpolation of the native ZEA maps."""
    gal = coords.galactic
    result = np.full(gal.shape, np.nan)
    for pole, sign in [('ngp', 1), ('sgp', -1)]:
        selected = gal.b.deg >= 0 if sign == 1 else gal.b.deg < 0
        if not np.any(selected):
            continue
        with fits.open(maps[pole], memmap=True) as hdul:
            header, data = hdul[0].header, hdul[0].data
            x, y = WCS(header).all_world2pix(gal.l.deg[selected], gal.b.deg[selected], 0)
            if np.any((x < 0) | (y < 0) | (x > data.shape[1]-1) | (y > data.shape[0]-1)):
                raise ValueError('SFD coordinate outside map')
            result[selected] = map_coordinates(data, [y, x], order=1, prefilter=False)
    if not np.all(np.isfinite(result)) or np.any(result < 0):
        raise ValueError('Invalid SFD reddening')
    return result


def catalogue_queries(frame):
    ident = str(frame['id'])
    if not re.fullmatch(r'[A-Za-z0-9_-]+', ident):
        raise ValueError('Unsafe cluster identifier')
    ra, dec = float(frame['ra_deg']), float(frame['dec_deg'])
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError('Invalid sky coordinate')
    # Cone-search radius is in arcsec, as specified by the official PDR3 API.
    # Keep all nearby forced objects; never silently choose a BCG counterpart.
    extinction = (f"SELECT '{ident}' AS target_cluster_id, object_id, ra, dec, a_i "
                  f"FROM pdr3_wide.forced WHERE isprimary "
                  f"AND conesearch(coord, {ra:.10f}, {dec:.10f}, 3.0);")
    random = (f"SELECT '{ident}' AS target_cluster_id, r.* FROM pdr3_wide.random AS r "
              f"WHERE conesearch(coord, {ra:.10f}, {dec:.10f}, 1800.0) "
              "AND i_pixelflags_bright_objectcenter IS FALSE AND i_pixelflags_edge IS FALSE "
              "ORDER BY md5(object_id::text || ':20261013') LIMIT 200;")
    return dict(cluster_id=ident, extinction_sql=extinction, random_candidates_sql=random)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True, type=Path)
    p.add_argument('--work', required=True, type=Path)
    p.add_argument('--cache', required=True, type=Path)
    args = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    args.work.mkdir(parents=True, exist_ok=False)
    args.cache.mkdir(parents=True, exist_ok=True)
    manifest_path = args.sample / 'download_01/download_manifest.json'
    status_path = args.sample / 'download_01/download_status.json'
    manifest = json.loads(manifest_path.read_text())
    status = json.loads(status_path.read_text())
    frames = manifest['clusters']
    done = {str(row['id']): row for row in status['completed']}
    if status['status'] != 'complete' or set(done) != {str(f['id']) for f in frames}:
        raise ValueError('Incomplete or inconsistent download manifest')
    maps, provenance = {}, {}
    for pole in ['ngp', 'sgp']:
        maps[pole], provenance[pole] = acquire_map(args.cache, pole)
        print('SFD map ready:', pole, flush=True)
    cosmology = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'])
    centres = SkyCoord([f['ra_deg'] for f in frames]*u.deg, [f['dec_deg'] for f in frames]*u.deg)
    # Nine samples are a gradient diagnostic, not a pixel-level extinction map.
    halfwidth = (manifest.get('requested_outer_radius_kpc', 1000) /
                 cosmology.kpc_proper_per_arcmin([f['z'] for f in frames]).value)*u.arcmin
    coordinates = [centres]
    for dx, dy in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]:
        coordinates.append(centres.spherical_offsets_by(dx*halfwidth, dy*halfwidth))
    samples = np.array([query_sfd(c, maps) for c in coordinates])
    rows = []
    for i, frame in enumerate(frames):
        filt = done[str(frame['id'])]['review']['coadd_header_filter']
        coefficient = COEFFICIENTS.get(filt)
        if coefficient is None:
            raise ValueError('Unreviewed filter: '+str(filt))
        ebv = float(samples[0,i])
        rows.append(dict(cluster_id=str(frame['id']), ra_deg=frame['ra_deg'], dec_deg=frame['dec_deg'],
                         header_filter=filt, ebv_sfd=ebv, ebv_9point_min=float(samples[:,i].min()),
                         ebv_9point_max=float(samples[:,i].max()), coefficient=coefficient,
                         a_i_nominal_mag=ebv*coefficient, flux_factor_nominal=10**(.4*ebv*coefficient),
                         pdr3_catalogue_crosscheck='pending'))
    with (args.work/'extinction_candidates.tsv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    requests = [catalogue_queries(frame) for frame in frames]
    write_json(args.work/'hsc_catalogue_requests.json', requests)
    (args.work/'patch_qa.sql').write_text('SELECT * FROM pdr3_wide.patch_qa;\n')
    write_json(args.work/'control_selection.json', dict(
        status='awaiting_official_random_and_patch_qa', preliminary_target=800,
        candidates_per_cluster=200, search_radius_arcmin=30, seed=20261013,
        query_validation='generated from documented tables/functions; not executed against CAS',
        requirements=['patch_inner and tract_inner; HSC-I coverage across entire cutout',
                      'match cluster sky region, depth, PSF, extinction and usable area',
                      'same angular cutout, assigned cluster redshift, masks and calibration',
                      'exclude overlap with all CAMIRA clusters, not only selected 800',
                      'inspect cluster-exclusion radius sensitivity and resulting selection bias',
                      'unique controls, with independence/overlap audit and regional jackknife',
                      'select final count from measured background covariance'],
        controls_selected=0, controls_downloaded=0, science_measurement=False))
    write_json(args.work/'preparation_status.json', dict(
        status='extinction_candidates_ready_controls_await_catalogues', n_clusters=len(frames),
        source_manifest=str(manifest_path), source_manifest_sha256=digest(manifest_path),
        source_status_sha256=digest(status_path), maps=provenance,
        sfd_scaling=1.0, coefficient_source=COEFFICIENT_SOURCE, coefficients=COEFFICIENTS,
        calibration_policy='Unscaled SFD times documented nominal HSC i/i2 coefficients; no extra 0.86 factor. '
                           'PDR3 effective passband/catalogue comparison remains required.',
        ebv_min=float(samples[0].min()), ebv_median=float(np.median(samples[0])), ebv_max=float(samples[0].max()),
        a_i_nominal_min=min(r['a_i_nominal_mag'] for r in rows),
        a_i_nominal_max=max(r['a_i_nominal_mag'] for r in rows),
        raw_images_modified=False, active_calibration_modified=False,
        catalogue_authentication_available=False, science_measurement=False))
    print('Prepared', len(rows), 'extinction candidates and paired catalogue requests', flush=True)


if __name__ == '__main__':
    main()
