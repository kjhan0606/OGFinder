#!/usr/bin/env python3
"""Write clusters.json and run the ICL menu stack when the inputs exist.

masses.tsv columns: id, m200c_msun, mass_definition, citation.
mass_definition must be 200c.

cutouts/frames.tsv columns: id, filename, sky, product.
sky must be global or global-sky. product is empty or deepCoadd.
A local sky or any other product name refuses the stack and writes no JSON.

The center is the astropy world_to_pixel result, already 0-based.
The pixel scale is the CD or CDELT plate scale in arcsec per pixel.
The survey checklist is left empty.
"""
import json
import math
import os
import subprocess
import sys

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')
OGF = os.environ.get('ASTRAFEX_OGF_ROOT', '/home/kjhan/BACKUP/Astrafek/OGFinder')
PYTHON = os.environ.get('ASTRAFEK_PYTHON', sys.executable)
STACK = os.path.join(OGF, 'plugins', 'icl', 'icl_stack.py')


def write_status(payload, code):
    path = os.path.join(ROOT, 'logs', 'stack_status.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=1)
        handle.write('\n')
    print('stack %s' % payload['status'])
    return code


def read_tsv(path):
    with open(path, encoding='utf-8') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        rows = []
        for line in handle:
            if line.strip():
                rows.append(line.rstrip('\n').split('\t'))
    return header, rows


def selected_rows():
    header, rows = read_tsv(os.path.join(ROOT, 'logs', 'selected_clusters.tsv'))
    if header[0] != 'id':
        raise SystemExit('selected_clusters.tsv has no id column')
    return rows


def pixel_scale_arcsec(header):
    if 'CD1_1' in header and 'CD2_1' in header:
        return 3600.0 * math.hypot(float(header['CD1_1']), float(header['CD2_1']))
    if 'CDELT1' in header:
        return 3600.0 * abs(float(header['CDELT1']))
    raise SystemExit('pixel scale is missing from the cutout header')


def center_xy(header, ra, dec):
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    import astropy.units as u
    x, y = WCS(header).world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    return [float(x), float(y)]


def science_header(path):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, 'ndim', 0) == 2:
                return hdu.header
    raise SystemExit('no 2-D image in %s' % path)


def main():
    masses_path = os.path.join(ROOT, 'masses.tsv')
    frames_path = os.path.join(ROOT, 'cutouts', 'frames.tsv')
    base = {
        'clusters_json_written': False,
        'stack_executed': False,
        'science_measurement': False,
    }
    if not os.path.isfile(masses_path) or not os.path.isfile(frames_path):
        base.update(status='stopped', reason='masses.tsv or cutouts/frames.tsv is absent')
        return write_status(base, 0)
    mass_header, mass_rows = read_tsv(masses_path)
    frame_header, frame_rows = read_tsv(frames_path)
    if mass_header != ['id', 'm200c_msun', 'mass_definition', 'citation']:
        base.update(status='refused', reason='masses.tsv header is unexpected')
        return write_status(base, 1)
    if frame_header != ['id', 'filename', 'sky', 'product']:
        base.update(status='refused', reason='frames.tsv header is unexpected')
        return write_status(base, 1)
    masses = {}
    for row in mass_rows:
        if len(row) != 4 or row[2] != '200c':
            base.update(status='refused', reason='a mass row is not 200c')
            return write_status(base, 1)
        masses[row[0]] = row
    frames = {}
    for row in frame_rows:
        if len(row) != 4:
            base.update(status='refused', reason='a frame row does not have 4 columns')
            return write_status(base, 1)
        sky = row[2].strip().lower()
        product = row[3].strip()
        if sky not in ('global', 'global-sky'):
            base.update(status='refused', reason='sky is not global for id %s' % row[0])
            return write_status(base, 1)
        if product and product != 'deepCoadd':
            base.update(status='refused', reason='product is not deepCoadd for id %s' % row[0])
            return write_status(base, 1)
        frames[row[0]] = row
    records = []
    for row in selected_rows():
        ident, _name, ra, dec, _zcl, _nmem, zspec = row[:7]
        if ident not in masses or ident not in frames:
            base.update(status='stopped', reason='inputs missing for id %s' % ident)
            return write_status(base, 0)
        filename = frames[ident][1]
        path = filename if os.path.isabs(filename) else os.path.join(ROOT, 'cutouts', filename)
        if not os.path.isfile(path):
            base.update(status='stopped', reason='cutout missing for id %s' % ident)
            return write_status(base, 0)
        header = science_header(path)
        record = {
            'id': ident,
            'image': path,
            'center': center_xy(header, float(ra), float(dec)),
            'pixel_scale': pixel_scale_arcsec(header),
            'mass': float(masses[ident][1]),
            'z': float(zspec),
            'sky': frames[ident][2].strip(),
        }
        if frames[ident][3].strip():
            record['product'] = frames[ident][3].strip()
        records.append(record)
    out_dir = os.path.join(ROOT, 'work')
    os.makedirs(out_dir, exist_ok=True)
    clusters = os.path.join(out_dir, 'clusters.json')
    with open(clusters, 'w', encoding='utf-8') as handle:
        json.dump({'clusters': records}, handle, indent=1)
        handle.write('\n')
    work = os.path.join(out_dir, 'iclstack')
    proc = subprocess.run(
        [PYTHON, STACK, '--clusters', clusters, '--work', work,
         '--u-max', '2.0', '--n-side', '201', '--delta', '200',
         '--h0', '70', '--omega-m', '0.3'],
        capture_output=True, text=True)
    payload = {
        'status': 'ran' if proc.returncode == 0 else 'failed',
        'returncode': proc.returncode,
        'stdout': proc.stdout.strip(),
        'stderr': proc.stderr.strip(),
        'clusters_json_written': True,
        'stack_executed': proc.returncode == 0,
        'science_measurement': False,
        'command': (
            'python plugins/icl/icl_stack.py --clusters clusters.json --work iclstack '
            '--u-max 2.0 --n-side 201 --delta 200 --h0 70 --omega-m 0.3'),
    }
    code = write_status(payload, 0 if proc.returncode == 0 else 1)
    return code


if __name__ == '__main__':
    sys.exit(main())
