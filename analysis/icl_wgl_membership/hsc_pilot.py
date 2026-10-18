#!/usr/bin/env python3
"""Prepare, inspect, and measure a small BCG+diffuse kpc pilot.

No downloader, credentials, mass conversion, or scientific certification.
Inspection reads headers. Measurement must run inside a Slurm allocation.
Existing catalog/debug/synthetic outputs are never overwritten.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from icl.pilot import (area_warp_metric, calibrate_flux, detection_mask,
                       native_profile, paired_summary)
from ogfkit.wlcluster import ARCSEC, Cosmology


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(plain(value), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')


def new_work(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def template_frame(row):
    return {'id': row['id'], 'name': row.get('name'), 'ra_deg': float(row['ra_deg']),
            'dec_deg': float(row['dec_deg']), 'z': float(row['z_spec_bcg']),
            'richness': float(row['N_mem']), 'region': None,
            'image': None, 'science_hdu': None, 'sky': None, 'product': None,
            'product_evidence': None, 'product_reviewed': False,
            'calibration': {'input_unit': None, 'ab_zeropoint': None, 'evidence': None,
                            'a_i_mag': None, 'extinction_evidence': None},
            'bad_mask': None, 'bright_star_mask': None, 'variance': None,
            'members_radec_deg': [], 'randoms': []}


def prepare(args):
    data = Path(args.data).resolve()
    with (data / 'logs/selected_clusters.tsv').open(encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    all_selected = getattr(args, 'all_selected', False)
    if not all_selected and (not 10 <= args.n_clusters <= 20 or len(rows) < args.n_clusters):
        raise ValueError('pilot size must be 10..20 and available in the selected list')
    if not rows:
        raise ValueError('selected cluster list is empty')
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('duplicate selected cluster ids')
    ordered = sorted(rows, key=lambda r: (float(r['N_mem']), r['id']))
    indices = (range(len(ordered)) if all_selected else
               np.rint(np.linspace(0, len(ordered) - 1, args.n_clusters)).astype(int))
    frames = [template_frame(ordered[i]) for i in indices]
    by_id = {f['id']: f for f in frames}
    with (data / 'logs/members.tsv').open(encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['id'] in by_id:
                by_id[row['id']]['members_radec_deg'].append([float(row['ra_deg']), float(row['dec_deg'])])
    cosmo = Cosmology(H0=70.0, Om=0.3)
    for frame in frames:
        kpa = float(cosmo.angular_diameter(frame['z'])) * ARCSEC * 1000.0
        frame['requested_half_width_arcsec'] = 1000.0 / kpa
    work = new_work(args.work)
    doc = {
        'schema': 'hsc-bcg-diffuse-pilot-v1', 'product': 'bcg+diffuse',
        'science_measurement': False, 'band': 'HSC-I', 'rest_frame': 'observed',
        'luminosity_unit': 'nJy/arcsec2, corrected for explicitly supplied Galactic extinction',
        'bcg_policy': 'keep', 'mass_used': False, 'h0': 70.0, 'omega_m': 0.3,
        'kpc_max': 500.0, 'n_side': 201, 'edges_kpc': list(np.linspace(0, 500, 21)),
        'min_coverage': 0.5, 'kron_factors': [2.5, 3.5, 5.0],
        'selection': '12 (or caller-selected 10..20) equally spaced ranks in richness; engineering pilot, not a representative sample',
        'parent_selection': str(data / 'logs/selected_clusters.tsv'),
        'requested_outer_radius_kpc': 1000.0,
        'outer_region_note': 'Field size allowance, not an assertion that the outer region has no cluster light.',
        'checks': {key: {'status': 'pending', 'evidence': None} for key in
                   ('sky_transfer', 'random_matching', 'mask_completeness', 'psf_wings',
                    'injection_recovery', 'profile_agreement', 'sample_selection')},
        'clusters': frames,
    }
    manifest_name = 'pilot_manifest.json'
    if all_selected:
        manifest_name = 'sample_manifest.json'
        doc.update(schema='hsc-bcg-diffuse-sample-v1',
                   selection='All rows of the existing selected_clusters.tsv; quality selection pending',
                   selection_limits='BCG spectroscopic matching, a narrow z_cl window and 6-arcmin separation; not a complete CAMIRA population')
        import hashlib
        doc['parent_selection_sha256'] = hashlib.sha256(
            (data / 'logs/selected_clusters.tsv').read_bytes()).hexdigest()
    write_json(work / manifest_name, doc)
    write_json(work / 'preparation.json', {'status': 'prepared', 'n_clusters': len(frames),
                                         'credentials_requested': False, 'downloaded': False,
                                         'stack_executed': False, 'science_measurement': False})
    print('Prepared %d clusters: %s' % (len(frames), work / manifest_name))


def resolve(base, filename):
    if not filename:
        raise ValueError('input filename is missing')
    path = Path(filename)
    return path if path.is_absolute() else base / path


def frame_metadata(row, base):
    """Require reviewed provenance and explicit HDUs; never classify sky from guesses."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as units
    if row.get('sky') not in ('global', 'global-sky') or row.get('product') not in ('deepCoadd', 'coadd/bg'):
        raise ValueError('reviewed global-sky deepCoadd or coadd/bg is required')
    if row.get('product_reviewed') is not True or not row.get('product_evidence'):
        raise ValueError('source product review/evidence is missing')
    if row.get('science_hdu') is None:
        raise ValueError('science_hdu must explicitly identify the science plane')
    calibrate_flux(np.zeros((1, 1)), row.get('calibration') or {})
    image_path = resolve(base, row.get('image'))
    with fits.open(image_path, memmap=True) as hdul:
        hdu = hdul[row['science_hdu']]
        header = hdu.header.copy()
        if header.get('INHERIT') is True:
            # Keep extension dimensions/WCS and inherit missing primary metadata.
            header.extend(hdul[0].header, unique=True)
        shape = tuple(hdu.shape)
    if len(shape) != 2:
        raise ValueError('science HDU is not a 2-D image')
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError('celestial WCS is required')
    centre = SkyCoord(float(row['ra_deg']) * units.deg, float(row['dec_deg']) * units.deg)
    cx, cy = wcs.world_to_pixel(centre)
    if not np.isfinite([cx, cy]).all() or not (0 <= cx < shape[1] and 0 <= cy < shape[0]):
        raise ValueError('catalog centre is outside the science image')
    matrix = np.asarray(wcs.pixel_scale_matrix, float)
    norms = np.sqrt(np.sum(matrix ** 2, axis=0))
    orthogonality = abs(np.dot(matrix[:, 0], matrix[:, 1]) / np.prod(norms))
    if not np.isfinite(norms).all() or np.any(norms <= 0) or abs(norms[0] / norms[1] - 1) > 1e-3 or orthogonality > 1e-3:
        raise ValueError('area diagnostic requires locally square, orthogonal WCS pixels')
    scale = 3600.0 * math_sqrt_det(matrix)
    probes = np.array([[cx, cy], [0, 0], [shape[1]-2, 0], [0, shape[0]-2], [shape[1]-2, shape[0]-2]])
    a = wcs.pixel_to_world(probes[:, 0], probes[:, 1])
    sx = a.separation(wcs.pixel_to_world(probes[:, 0] + 1, probes[:, 1])).arcsec
    sy = a.separation(wcs.pixel_to_world(probes[:, 0], probes[:, 1] + 1)).arcsec
    variation = float(np.max(np.abs(np.concatenate([sx, sy]) / scale - 1)))
    if not np.isfinite(variation) or variation > 1e-3:
        raise ValueError('WCS scale varies by >0.1%; a full geometric reprojection is required')
    bunit = str(header.get('BUNIT', '')).strip()
    if bunit and ('mag' in bunit.lower() or 'sr' in bunit.lower()):
        raise ValueError('science plane must contain per-pixel linear flux')
    if 'njy' in bunit.lower() and row['calibration']['input_unit'] != 'nJy/pixel':
        raise ValueError('BUNIT nJy conflicts with instrumental-flux calibration')
    fluxmag0 = header.get('FLUXMAG0')
    if fluxmag0 is not None and row['calibration']['input_unit'] == 'instrumental_flux/pixel':
        fluxmag0 = float(fluxmag0)
        if not np.isfinite(fluxmag0) or fluxmag0 <= 0:
            raise ValueError('invalid inherited FLUXMAG0')
        if abs(2.5 * np.log10(fluxmag0) - float(row['calibration']['ab_zeropoint'])) > 1e-6:
            raise ValueError('explicit zero point conflicts with inherited FLUXMAG0')
    for name in ('bad_mask', 'bright_star_mask', 'variance'):
        spec = row.get(name)
        if not isinstance(spec, dict) or spec.get('hdu') is None or not spec.get('evidence') or spec.get('same_grid_reviewed') is not True:
            raise ValueError('%s requires filename, HDU, evidence, and same_grid_reviewed' % name)
        with fits.open(resolve(base, spec.get('file')), memmap=True) as hdul:
            if tuple(hdul[spec['hdu']].shape) != shape:
                raise ValueError('%s shape differs from science image' % name)
    return {'id': str(row['id']), 'shape': shape, 'center': [float(cx), float(cy)],
            'pixel_scale_arcsec': scale, 'wcs_scale_variation': variation,
            'bunit': bunit, 'product': row['product'], 'sky': row['sky'],
            'header_fluxmag0': header.get('FLUXMAG0'), 'science_measurement': False}


def math_sqrt_det(matrix):
    return float(np.sqrt(abs(np.linalg.det(matrix))))


def load_manifest(path):
    path = Path(path).resolve()
    doc = json.loads(path.read_text(encoding='utf-8'))
    if doc.get('schema') != 'hsc-bcg-diffuse-pilot-v1' or doc.get('band') != 'HSC-I' or doc.get('rest_frame') != 'observed' or doc.get('bcg_policy') != 'keep':
        raise ValueError('unexpected pilot schema/measurement definition')
    clusters = doc.get('clusters') or []
    if not clusters or len({str(r['id']) for r in clusters}) != len(clusters):
        raise ValueError('pilot clusters must have unique ids')
    edges = np.asarray(doc.get('edges_kpc'), dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or edges[0] < 0 or np.any(np.diff(edges) <= 0) or edges[-1] > float(doc['kpc_max']):
        raise ValueError('radial bins must increase within the diagnostic image radius')
    factors = doc.get('kron_factors') or []
    if not factors or len(set(factors)) != len(factors) or not all(np.isfinite(v) and v > 0 for v in factors):
        raise ValueError('Kron factors must be unique and positive')
    if not 0 < float(doc['min_coverage']) <= 1:
        raise ValueError('min_coverage must lie in (0,1]')
    if int(doc['n_side']) != doc['n_side'] or int(doc['n_side']) < 3:
        raise ValueError('n_side must be an integer >= 3')
    for row in clusters:
        if row.get('randoms') and not row.get('random_matching_evidence'):
            raise ValueError('random-field matching evidence is required')
        random_ids = [str(r['id']) for r in row.get('randoms', [])]
        if len(set(random_ids)) != len(random_ids):
            raise ValueError('duplicate random fields in a matched set')
        for random in row.get('randoms', []):
            if float(random['z']) != float(row['z']):
                raise ValueError('random fields must use their lens redshift for the metric grid')
    return doc, path.parent


def inspect(args):
    doc, base = load_manifest(args.manifest)
    results = []
    for cluster in doc['clusters']:
        for row in [cluster] + cluster.get('randoms', []):
            try:
                results.append(dict(frame_metadata(row, base), status='headers_checked'))
            except (ValueError, KeyError, TypeError, OSError, IndexError) as exc:
                results.append({'id': str(row['id']), 'status': 'stopped', 'reason': str(exc)})
    work = new_work(args.work)
    stopped = sum(r['status'] == 'stopped' for r in results)
    write_json(work / 'inspection.json', {'frames': results, 'n_stopped': stopped,
                                         'stack_executed': False, 'science_measurement': False})
    print('Inspected %d frames; %d stopped. %s' % (len(results), stopped, work))
    return 2 if stopped else 0


def load_plane(base, spec):
    from astropy.io import fits
    with fits.open(resolve(base, spec['file']), memmap=False) as hdul:
        return np.array(hdul[spec['hdu']].data, copy=True)


def load_bool_mask(base, spec):
    data = load_plane(base, spec)
    if not np.isfinite(data).all():
        raise ValueError('mask contains non-finite values')
    if spec.get('kind') == 'bits':
        bits = spec.get('exclude_bits')
        if not bits or any(not isinstance(v, int) or v <= 0 or v & (v - 1) for v in bits):
            raise ValueError('exclude_bits must name positive powers of two, not bit indices')
        if not np.equal(data, np.floor(data)).all() or np.any(data < 0):
            raise ValueError('bit mask must contain non-negative integers')
        flag = 0
        for bit in bits:
            flag |= bit
        return (data.astype(np.int64) & flag) != 0
    if spec.get('kind') != 'boolean' or not np.isin(data, [0, 1]).all():
        raise ValueError('boolean mask must explicitly contain only 0 and 1')
    return data.astype(bool)


def measure_frame(row, base, doc, factor, random=False):
    from astropy.io import fits
    from astropy.wcs import WCS
    meta = frame_metadata(row, base)
    with fits.open(resolve(base, row['image']), memmap=False) as hdul:
        image = np.array(hdul[row['science_hdu']].data, dtype=float, copy=True)
        header = hdul[row['science_hdu']].header.copy()
        if header.get('INHERIT') is True:
            header.extend(hdul[0].header, unique=True)
        wcs = WCS(header).celestial
    bad = load_bool_mask(base, row['bad_mask'])
    bright = load_bool_mask(base, row['bright_star_mask'])
    variance = load_plane(base, row['variance'])
    bad |= ~np.isfinite(variance) | (variance <= 0)
    positions = row.get('members_radec_deg', [])
    members = np.array(wcs.all_world2pix(np.asarray(positions, float), 0)).tolist() if positions else []
    mask, detection = detection_mask(image, meta['center'], meta['pixel_scale_arcsec'], bad, bright,
                                     kron_factor=factor, bcg_policy='none' if random else 'keep',
                                     member_xy=members)
    calibrated, conversion = calibrate_flux(image, row['calibration'])
    cosmo = Cosmology(H0=doc['h0'], Om=doc['omega_m'])
    profile = native_profile(calibrated, mask, meta['center'], meta['pixel_scale_arcsec'], row['z'],
                             doc['edges_kpc'], cosmo, doc['min_coverage'])
    warped = area_warp_metric(calibrated, mask, meta['center'], meta['pixel_scale_arcsec'], row['z'],
                              kpc_max=doc['kpc_max'], n_side=doc['n_side'], cosmo=cosmo)
    meta.update(njy_per_input_unit=conversion, detection=detection, profile=profile)
    return meta, warped


def measure(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('measure must run in a Slurm allocation; use pilot_job.sh')
    doc, base = load_manifest(args.manifest)
    # Validate every header before creating output or reading large data planes.
    for cluster in doc['clusters']:
        for row in [cluster] + cluster.get('randoms', []):
            frame_metadata(row, base)
    work = new_work(args.work)
    from astropy.io import fits
    all_results = []
    for factor in doc['kron_factors']:
        records, frames, random_profiles = [], [], []
        for row in doc['clusters']:
            meta, frame = measure_frame(row, base, doc, factor)
            records.append(meta); frames.append(frame)
            random_records = []
            for random in row.get('randoms', []):
                rmeta, _ = measure_frame(random, base, doc, factor, random=True)
                random_records.append(rmeta)
            random_profiles.append([r['profile'] for r in random_records])
            meta['random_frames'] = random_records
        summary = paired_summary([r['profile'] for r in records], random_profiles,
                                 [r.get('region') for r in doc['clusters']])
        good = np.array([(f['coverage'] >= doc['min_coverage']) & np.isfinite(f['sb']) for f in frames])
        values = np.array([f['sb'] for f in frames])
        count = np.sum(good, axis=0)
        image = np.full(count.shape, np.nan)
        np.divide(np.sum(np.where(good, values, 0), axis=0), count, out=image, where=count > 0)
        axis = frames[0]['kpc']
        yy, xx = np.meshgrid(axis, axis, indexing='ij')
        radius = np.hypot(xx, yy)
        stack_profile = []
        for lo, hi in zip(doc['edges_kpc'][:-1], doc['edges_kpc'][1:]):
            use = (radius >= lo) & (radius < hi) & np.isfinite(image)
            stack_profile.append(float(np.mean(image[use])) if np.any(use) else None)
        tag = ('kron_%g' % factor).replace('.', 'p')
        out = work / tag
        out.mkdir()
        header = fits.Header({'BUNIT': 'nJy/arcsec2', 'SCIENCE': False, 'PRODUCT': 'bcg+diffuse',
                              'AXISNOTE': 'input image axes; radial diagnostic only', 'KPCMAX': doc['kpc_max'],
                              'CRPIX1': (len(axis)+1)/2, 'CRPIX2': (len(axis)+1)/2,
                              'CRVAL1': 0.0, 'CRVAL2': 0.0,
                              'CDELT1': axis[1]-axis[0], 'CDELT2': axis[1]-axis[0],
                              'CTYPE1': 'LINEAR', 'CTYPE2': 'LINEAR', 'CUNIT1': 'kpc', 'CUNIT2': 'kpc'})
        fits.PrimaryHDU(image, header).writeto(out / 'kpc_stack.fits')
        coverage_header = header.copy(); coverage_header['BUNIT'] = 'clusters'
        fits.PrimaryHDU(count.astype(float), coverage_header).writeto(out / 'cluster_count.fits')
        write_json(out / 'frames.json', records)
        write_json(out / 'covariance.json', summary)
        with (out / 'profile.tsv').open('x') as handle:
            handle.write('R_MIN_KPC\tR_MAX_KPC\tRAW_SB\tPAIRED_RAW_SB\tRANDOM_SB\tCORRECTED_SB\tN_RAW\tN_PAIRED\tJK_SIGMA\tIMAGE_SB\tIMAGE_MINUS_NATIVE\n')
            for i, (lo, hi) in enumerate(zip(doc['edges_kpc'][:-1], doc['edges_kpc'][1:])):
                sigma = np.sqrt(max(0, summary['covariance'][i, i])) if np.isfinite(summary['covariance'][i, i]) else None
                diff = stack_profile[i] - summary['unpaired_raw_sb'][i] if stack_profile[i] is not None else None
                values = [lo, hi, summary['unpaired_raw_sb'][i], summary['raw_sb'][i], summary['random_sb'][i], summary['corrected_sb'][i],
                          summary['unpaired_n_clusters'][i], summary['n_clusters'][i], sigma, stack_profile[i], diff]
                handle.write('\t'.join('NA' if plain(v) is None else str(plain(v)) for v in values) + '\n')
        all_results.append({'kron_factor': factor, 'path': tag, 'n_regions': summary['n_regions'],
                            'covariance_bins': summary['covariance_bins'], 'science_measurement': False})
    write_json(work / 'manifest_used.json', doc)
    write_json(work / 'pilot_status.json', {'status': 'diagnostics_written', 'variants': all_results,
                                           'mass_used': False, 'science_measurement': False,
                                           'checks': doc.get('checks'),
                                           'note': 'Raw diagnostic images; sky-transfer, PSF and injection validation remain required.'})
    print('Pilot diagnostics written: %s; science_measurement false' % work)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    tasks = parser.add_subparsers(dest='task', required=True)
    prep = tasks.add_parser('prepare')
    prep.add_argument('--data', default='/gpfs/kjhan/Astrafek/data')
    prep.add_argument('--n-clusters', type=int, default=12)
    prep.add_argument('--all-selected', action='store_true',
                      help='prepare every selected candidate for a separate sample download')
    prep.add_argument('--work', required=True)
    prep.set_defaults(func=prepare)
    for name, func in [('inspect', inspect), ('measure', measure)]:
        task = tasks.add_parser(name)
        task.add_argument('--manifest', required=True)
        task.add_argument('--work', required=True)
        task.set_defaults(func=func)
    args = parser.parse_args(argv)
    try:
        return args.func(args) or 0
    except (ValueError, KeyError, TypeError, OSError, IndexError) as exc:
        parser.exit(2, 'hsc_pilot: %s\n' % exc)


if __name__ == '__main__':
    sys.exit(main())
