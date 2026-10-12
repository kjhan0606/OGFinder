#!/usr/bin/env python3
"""Scaled BCG+diffuse stack for the ICL menu.

    icl_stack.py --clusters clusters.json --work DIR

The one-frame profile remains ds9/library/ds9_icl.py. This program reads a JSON
list of frames and calls icl.wgl.stack_clusters. Each frame needs an image, a
centre, a pixel scale, a mass and a redshift. A frame whose sky is local, or
whose named product is not deepCoadd, is refused. science_measurement follows
the survey checklist in the JSON and stays false while that checklist is empty.
The frame count is not a sky measurement.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from icl.wgl import stack_clusters  # noqa: E402
from ogfkit.imageio import load_image, load_mask, save_fits  # noqa: E402
from ogfkit.wlcluster import Cosmology  # noqa: E402

_LOCAL_SKY = ('local', 'local-sky')
_GLOBAL_SKY = ('global', 'global-sky')
_NOTE = (
    'The stack is the mean of the supplied frames. '
    'science_measurement follows the survey checklist. '
    'An empty checklist is not a named-data measurement.')


def _fail(message):
    raise SystemExit('icl_stack: %s' % message)


def _finite(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        _fail('%s must be finite' % name)
    if not np.isfinite(number):
        _fail('%s must be finite' % name)
    return number


def _centre(value):
    if isinstance(value, dict) and 'x' in value and 'y' in value:
        pair = (value['x'], value['y'])
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        pair = value
    else:
        _fail('center must be [x, y]')
    return [_finite(pair[0], 'center x'), _finite(pair[1], 'center y')]


def _resolve(base, path):
    path = str(path)
    if os.path.isabs(path):
        return path
    return os.path.join(base, path)


def _sky_ok(record):
    if 'sky' not in record or record['sky'] is None or str(record['sky']).strip() == '':
        return
    sky = str(record['sky']).strip().lower()
    if sky in _LOCAL_SKY:
        _fail('local-sky is not an ICL input')
    if sky not in _GLOBAL_SKY:
        _fail('sky was not named global')


def _product_ok(record):
    if 'product' not in record or record['product'] is None or str(record['product']).strip() == '':
        return
    if str(record['product']).strip() != 'deepCoadd':
        _fail('product is not deepCoadd')


def _load_records(path):
    if not path or not os.path.isfile(path):
        _fail('clusters file was not found: %s' % (path or ''))
    try:
        with open(path, encoding='utf-8') as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        _fail('clusters JSON could not be read (%s)' % exc)
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict) and isinstance(payload.get('clusters'), list):
        rows = payload['clusters']
    else:
        _fail('clusters JSON needs a list of frames')
    if not rows:
        _fail('clusters JSON has no frames')
    base = os.path.dirname(os.path.abspath(path))
    survey = payload.get('survey') if isinstance(payload, dict) else None
    rest_frame = payload.get('rest_frame') if isinstance(payload, dict) else None
    records = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            _fail('each frame must be a record')
        _sky_ok(row)
        _product_ok(row)
        for key in ('image', 'center', 'pixel_scale', 'mass', 'z'):
            if key not in row:
                _fail('frame %s needs %s' % (index, key))
        image_path = _resolve(base, row['image'])
        if not os.path.isfile(image_path):
            _fail('image was not found: %s' % image_path)
        image, _header = load_image(image_path)
        record = {
            'image': image,
            'center': _centre(row['center']),
            'pixel_scale': _finite(row['pixel_scale'], 'pixel_scale'),
            'mass': _finite(row['mass'], 'mass'),
            'z': _finite(row['z'], 'z'),
            'id': str(row['id']) if row.get('id') is not None else str(index),
        }
        if row.get('mask'):
            record['mask'] = load_mask(_resolve(base, row['mask']), image.shape)
        if row.get('variance'):
            variance, _header = load_image(_resolve(base, row['variance']))
            if variance.shape != image.shape:
                _fail('variance shape does not match the image')
            record['variance'] = variance
        records.append(record)
    return records, survey, rest_frame


def _plain_bins(radial):
    rows = []
    for item in radial.get('bins') or []:
        row = {}
        for key, value in item.items():
            if key == 'ids':
                row[key] = [str(name) for name in value]
            elif isinstance(value, (bool, np.bool_)):
                row[key] = bool(value)
            elif isinstance(value, (int, np.integer)):
                row[key] = int(value)
            else:
                row[key] = float(value)
        rows.append(row)
    return {
        'bins': rows,
        'sample_changes_with_radius': bool(radial.get('sample_changes_with_radius')),
        'note': radial.get('note'),
    }


def _write_profile(path, bins, low_key, high_key, low_name, high_name):
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('%s\t%s\tN_CLUSTERS\tIDS\n' % (low_name, high_name))
        for row in bins:
            handle.write('%s\t%s\t%s\t%s\n' % (
                row[low_key], row[high_key], row['n_clusters'], ','.join(row['ids'])))


def _summary(stacked, inverse_variance):
    survey = stacked['survey']
    return {
        'product': stacked['product'],
        'description': stacked['description'],
        'delta': stacked['delta'],
        'n': int(stacked['n']),
        'ids': list(stacked['ids']),
        'science_measurement': bool(survey['science_measurement']),
        'survey_empty': list(survey['empty']),
        'surface_brightness_note': stacked['surface_brightness_note'],
        'radial': _plain_bins(stacked['radial']),
        'kpc_max': float(stacked['kpc']['kpc_max']),
        'kpc_max_source': stacked['kpc']['kpc_max_source'],
        'kpc_radial': _plain_bins(stacked['kpc']['radial']),
        'inverse_variance': bool(inverse_variance),
        'independent_evidence': False,
        'named_data_mass': False,
        'note': _NOTE,
    }


def run(args):
    records, survey, rest_frame = _load_records(args.clusters)
    u_max = _finite(args.u_max, 'u-max')
    n_side = int(args.n_side)
    delta = _finite(args.delta, 'delta')
    if n_side < 3:
        _fail('n-side must be at least 3')
    if u_max <= 0 or delta <= 0:
        _fail('u-max and delta must be positive')
    h0 = _finite(args.h0, 'h0')
    omega_m = _finite(args.omega_m, 'omega-m')
    if h0 <= 0 or omega_m <= 0:
        _fail('h0 and omega-m must be positive')
    kpc_max = None
    if str(args.kpc_max).strip() != '':
        kpc_max = _finite(args.kpc_max, 'kpc-max')
        if kpc_max <= 0:
            _fail('kpc-max must be positive')
    cosmo = Cosmology(H0=h0, Om=omega_m)
    try:
        stacked = stack_clusters(
            records, cosmo=cosmo, delta=delta, u_max=u_max, n_side=n_side, kpc_max=kpc_max,
            survey=survey, rest_frame=rest_frame)
    except ValueError as exc:
        _fail(str(exc))
    os.makedirs(args.work, exist_ok=True)
    inverse = 'inverse_variance_sb' in stacked
    save_fits(os.path.join(args.work, 'icl_stack_u.fits'), stacked['sb'],
              extra={'PRODUCT': 'bcg+diffuse', 'STACK': 'u'})
    save_fits(os.path.join(args.work, 'icl_stack_u_weight.fits'), stacked['weight'],
              extra={'PRODUCT': 'bcg+diffuse', 'STACK': 'u-coverage'})
    save_fits(os.path.join(args.work, 'icl_stack_kpc.fits'), stacked['kpc']['sb'],
              extra={'PRODUCT': 'bcg+diffuse', 'STACK': 'kpc'})
    save_fits(os.path.join(args.work, 'icl_stack_kpc_weight.fits'), stacked['kpc']['weight'],
              extra={'PRODUCT': 'bcg+diffuse', 'STACK': 'kpc-coverage'})
    if inverse:
        save_fits(os.path.join(args.work, 'icl_stack_u_ivar.fits'), stacked['inverse_variance_sb'],
                  extra={'PRODUCT': 'bcg+diffuse', 'STACK': 'u-inverse-variance'})
    _write_profile(
        os.path.join(args.work, 'icl_stack_profile.tsv'),
        _plain_bins(stacked['radial'])['bins'], 'u_min', 'u_max', 'U_MIN', 'U_MAX')
    _write_profile(
        os.path.join(args.work, 'icl_stack_kpc_profile.tsv'),
        _plain_bins(stacked['kpc']['radial'])['bins'], 'kpc_min', 'kpc_max', 'KPC_MIN', 'KPC_MAX')
    summary = _summary(stacked, inverse)
    with open(os.path.join(args.work, 'icl_stack.json'), 'w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=1)
        handle.write('\n')
    print('clusters %d product %s science_measurement %s' % (
        summary['n'], summary['product'], 'true' if summary['science_measurement'] else 'false'))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='Scaled BCG+diffuse stack of cluster frames.')
    parser.add_argument('--clusters', required=True)
    parser.add_argument('--work', required=True)
    parser.add_argument('--u-max', type=float, default=2.0)
    parser.add_argument('--n-side', type=int, default=201)
    parser.add_argument('--delta', type=float, default=200.0)
    parser.add_argument('--h0', type=float, default=70.0)
    parser.add_argument('--omega-m', type=float, default=0.3)
    parser.add_argument('--kpc-max', default='')
    args = parser.parse_args(argv)
    run(args)


if __name__ == '__main__':
    main()
