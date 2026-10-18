#!/usr/bin/env python3
"""Download prepared HSC-I pilot or sample frames, including image/mask/variance.

Run from an interactive syntax SSH terminal. Credentials stay in this process
memory. Each FITS is streamed to a new, exclusive work directory, validated,
then renamed into place. The script never downloads catalogs or random fields.
"""
import argparse
import csv
import getpass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np


ROOT = Path('/gpfs/kjhan/Astrafek/data')
API = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout?'
MANUAL = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/manual.html'
MAX_FILE_BYTES = 512 * 1024 * 1024
CHUNK_BYTES = 1024 * 1024

# MP_* values in the HSC mask header are bit positions; these are bit values.
# DETECTED is left to the source-mask stage so the BCG can remain in the image.
BAD_BIT_POSITIONS = (0, 1, 2, 3, 4, 7, 8, 10, 12, 13, 14, 15, 16)
BRIGHT_OBJECT_BIT = 1 << 9
EXPECTED_MASK_PLANES = {
    'MP_BAD': 0, 'MP_SAT': 1, 'MP_INTRP': 2, 'MP_CR': 3, 'MP_EDGE': 4,
    'MP_DETECTED': 5, 'MP_DETECTED_NEGATIVE': 6, 'MP_SUSPECT': 7,
    'MP_NO_DATA': 8, 'MP_BRIGHT_OBJECT': 9, 'MP_CROSSTALK': 10,
    'MP_NOT_DEBLENDED': 11, 'MP_UNMASKEDNAN': 12, 'MP_REJECTED': 13,
    'MP_CLIPPED': 14, 'MP_SENSOR_EDGE': 15, 'MP_INEXACT_PSF': 16,
}


def write_json(path, value):
    temp = Path(str(path) + '.tmp')
    with temp.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while True:
            block = handle.read(CHUNK_BYTES)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def estimate_bytes(frames):
    total = 0
    rows = []
    for frame in frames:
        half = float(frame['requested_half_width_arcsec'])
        if not math.isfinite(half) or not 30 <= half <= 600:
            raise ValueError('invalid requested half-width for cluster %s' % frame['id'])
        n = math.ceil(2 * half / 0.168)
        estimate = n * n * 12  # float32 image + int32 mask + float32 variance
        total += estimate
        rows.append({'id': str(frame['id']), 'estimated_pixels_per_side': n,
                     'estimated_bytes': estimate})
    return total, rows


def image_hdu_map(hdul):
    found = {}
    for index, hdu in enumerate(hdul):
        kind = str(hdu.header.get('EXTTYPE', '')).strip().upper()
        if kind in ('IMAGE', 'MASK', 'VARIANCE'):
            if kind in found:
                raise ValueError('duplicate %s HDUs' % kind)
            found[kind] = index
    missing = {'IMAGE', 'MASK', 'VARIANCE'} - set(found)
    if missing:
        raise ValueError('missing FITS plane(s): %s' % ', '.join(sorted(missing)))
    return found


def validate_fits(path, frame):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as units

    with fits.open(path, memmap=True, checksum=True) as hdul:
        size = Path(path).stat().st_size
        for hdu in hdul:
            info = hdu.fileinfo()
            if info and info['datLoc'] + info['datSpan'] > size:
                raise ValueError('truncated FITS data plane')
        primary = hdul[0].header
        planes = image_hdu_map(hdul)
        image_index = planes['IMAGE']
        mask_index = planes['MASK']
        variance_index = planes['VARIANCE']
        image_hdu = hdul[image_index]
        image_header = image_hdu.header.copy()
        if image_header.get('INHERIT') is True:
            image_header.extend(primary, unique=True)
        shape = tuple(image_hdu.shape)
        if len(shape) != 2 or min(shape) < 3:
            raise ValueError('science plane must be a 2-D image')
        source_wcs = WCS(image_header).celestial
        if not source_wcs.has_celestial:
            raise ValueError('science plane lacks celestial WCS')
        sample_x = np.array([0, shape[1] - 1, 0, shape[1] - 1, (shape[1] - 1) / 2])
        sample_y = np.array([0, 0, shape[0] - 1, shape[0] - 1, (shape[0] - 1) / 2])
        sample_sky = source_wcs.pixel_to_world(sample_x, sample_y)
        for label, index in [('MASK', mask_index), ('VARIANCE', variance_index)]:
            if tuple(hdul[index].shape) != shape:
                raise ValueError('%s plane shape differs from science plane' % label)
            test_header = hdul[index].header.copy()
            if test_header.get('INHERIT') is True:
                test_header.extend(primary, unique=True)
            test_wcs = WCS(test_header).celestial
            tx, ty = test_wcs.world_to_pixel(sample_sky)
            sx, sy = sample_x, sample_y
            if (not np.isfinite([*tx, *ty]).all()
                    or max(math.hypot(float(dx), float(dy))
                           for dx, dy in zip(tx - sx, ty - sy)) > 0.01):
                raise ValueError('%s WCS is not aligned with science plane' % label)
        primary_filter = str(primary.get('FILTER', '')).strip().lower()
        if primary_filter not in ('i', 'i2'):
            raise ValueError('primary FILTER %r is not an HSC-I-family filter' % primary_filter)
        fluxmag0 = float(primary.get('FLUXMAG0', float('nan')))
        if not math.isfinite(fluxmag0) or fluxmag0 <= 0:
            raise ValueError('primary FLUXMAG0 is missing or invalid')
        wcs = source_wcs
        x, y = wcs.world_to_pixel(SkyCoord(float(frame['ra_deg']) * units.deg,
                                           float(frame['dec_deg']) * units.deg))
        if not np.isfinite([x, y]).all() or not (0 <= x < shape[1] and 0 <= y < shape[0]):
            raise ValueError('catalog center falls outside the cutout')
        matrix = np.asarray(wcs.pixel_scale_matrix, float)
        scales = np.sqrt(np.sum(matrix ** 2, axis=0))
        orth = abs(float(np.dot(matrix[:, 0], matrix[:, 1]) / np.prod(scales)))
        scale = 3600 * math.sqrt(abs(float(np.linalg.det(matrix))))
        if (not np.isfinite(scales).all() or np.any(scales <= 0)
                or abs(scales[0] / scales[1] - 1) > 1e-3 or orth > 1e-3
                or abs(scale / 0.168 - 1) > 1e-3):
            raise ValueError('WCS pixels do not match the reviewed HSC square grid')
        bunit = str(image_header.get('BUNIT', '')).strip()
        if bunit and 'mag' in bunit.lower():
            raise ValueError('science plane is not linear flux')
        if bunit and 'njy' in bunit.lower():
            unit = 'nJy/pixel'
        elif not bunit:
            unit = 'instrumental_flux/pixel'
        else:
            raise ValueError('unreviewed science BUNIT: %s' % bunit)
        mask_header = hdul[mask_index].header
        for plane in ('MP_BAD', 'MP_BRIGHT_OBJECT', 'MP_DETECTED', 'MP_NO_DATA'):
            if plane not in mask_header:
                raise ValueError('MASK header omits %s' % plane)
        mask_bits = {key: int(mask_header[key]) for key in mask_header if key.startswith('MP_')}
        if mask_bits != EXPECTED_MASK_PLANES:
            raise ValueError('HSC mask-plane mapping differs from the reviewed probe')
        pixel_width = shape[1] * scale
        pixel_height = shape[0] * scale
        required = 2 * float(frame['requested_half_width_arcsec'])
        if min(pixel_width, pixel_height) < 0.98 * required:
            raise ValueError('cutout is smaller than the requested 1000-kpc half-width')
    return {'shape': shape, 'science_hdu': image_index, 'mask_hdu': mask_index,
            'variance_hdu': variance_index, 'pixel_scale_arcsec': scale,
            'center_pixel_xy': [float(x), float(y)], 'filter': 'HSC-I',
            'coadd_header_filter': primary_filter,
            'fluxmag0': fluxmag0, 'ab_zeropoint': 2.5 * math.log10(fluxmag0),
            'bunit': bunit, 'input_unit': unit, 'mask_planes': mask_bits}


def build_manifest_row(frame, filename, review):
    evidence = ('HSC PDR3 cutout manual maps coadd/bg to deepCoadd before aggressive '
                'background subtraction; one-probe header review is logs/cutout_global_review.json')
    common = {'file': filename, 'hdu': review['mask_hdu'],
              'evidence': 'HSC MASK HDU, coadd/bg header MP_* definitions; same WCS and shape reviewed',
              'same_grid_reviewed': True, 'kind': 'bits'}
    result = dict(frame)
    result.update(image=filename, science_hdu=review['science_hdu'], sky='global-sky',
                  product='coadd/bg', product_reviewed=True, product_evidence=evidence,
                  band='HSC-I', coadd_header_filter=review['coadd_header_filter'],
                  filter_evidence=('HSC PDR3 coadds combine i and i2 variants; preserve the '
                                   'actual primary FILTER value for audit: '
                                   'https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/'),
                  calibration={'input_unit': review['input_unit'],
                               'ab_zeropoint': review['ab_zeropoint'] if review['input_unit'] == 'instrumental_flux/pixel' else None,
                               'evidence': 'HSC PDR3 cutout manual flux calibration and inherited primary FLUXMAG0',
                               'a_i_mag': None, 'extinction_evidence': None},
                  bad_mask=dict(common, exclude_bits=[1 << bit for bit in BAD_BIT_POSITIONS]),
                  bright_star_mask=dict(common, exclude_bits=[BRIGHT_OBJECT_BIT]),
                  variance={'file': filename, 'hdu': review['variance_hdu'],
                            'evidence': 'HSC VARIANCE HDU; extension grid matched to IMAGE',
                            'same_grid_reviewed': True},
                  randoms=[])
    return result


def fetch_one(opener, frame, outdir, index, count):
    half = float(frame['requested_half_width_arcsec'])
    params = urllib.parse.urlencode({
        'ra': format(float(frame['ra_deg']), '.9f'),
        'dec': format(float(frame['dec_deg']), '.9f'),
        'sw': '%.6farcsec' % half, 'sh': '%.6farcsec' % half,
        'type': 'coadd/bg', 'image': 'on', 'mask': 'on', 'variance': 'on',
        'filter': 'HSC-I', 'rerun': 'pdr3_wide',
    })
    request = urllib.request.Request(API + params)
    safe_id = re.sub(r'[^A-Za-z0-9_-]', '_', str(frame['id']))
    final_path = outdir / ('cluster_%s.fits' % safe_id)
    part_path = outdir / ('cluster_%s.fits.part' % safe_id)
    if final_path.exists() or part_path.exists():
        raise ValueError('output for cluster %s already exists' % frame['id'])
    digest = hashlib.sha256()
    total = 0
    first = b''
    try:
        with opener.open(request, timeout=180) as response:
            content_type = response.headers.get('Content-Type', '')
            with part_path.open('xb') as handle:
                while True:
                    block = response.read(CHUNK_BYTES)
                    if not block:
                        break
                    if not first:
                        first = block[:80]
                    total += len(block)
                    if total > MAX_FILE_BYTES:
                        raise ValueError('cluster %s exceeds the 512-MiB per-file safety limit' % frame['id'])
                    digest.update(block)
                    handle.write(block)
                handle.flush()
                os.fsync(handle.fileno())
        if not first.startswith(b'SIMPLE'):
            raise ValueError('cluster %s returned a non-FITS body (%s)' % (frame['id'], content_type))
        review = validate_fits(part_path, frame)
        os.replace(part_path, final_path)
        print('[%d/%d] cluster %s: saved %d bytes; %s' %
              (index, count, frame['id'], total, final_path.name), flush=True)
        return {'id': str(frame['id']), 'file': final_path.name, 'bytes': total,
                'sha256': digest.hexdigest(), 'content_type': content_type,
                'review': review}, build_manifest_row(frame, final_path.name, review)
    except urllib.error.HTTPError as exc:
        raise ValueError('cluster %s request failed with HTTP %d' % (frame['id'], exc.code)) from None
    except urllib.error.URLError as exc:
        raise ValueError('cluster %s request failed (%s)' %
                         (frame['id'], type(exc.reason).__name__)) from None


def recover_existing(work, frames, status, status_path):
    images = work / 'images'
    prior = {}
    for artifact in status.get('completed', []):
        ident = str(artifact.get('id'))
        if ident not in {str(frame['id']) for frame in frames} or ident in prior:
            raise ValueError('resume status has an unknown or duplicate completed cluster')
        prior[ident] = artifact

    downloaded, output_frames = {}, {}
    for frame in frames:
        ident = str(frame['id'])
        safe_id = re.sub(r'[^A-Za-z0-9_-]', '_', ident)
        final_path = images / ('cluster_%s.fits' % safe_id)
        part_path = images / ('cluster_%s.fits.part' % safe_id)
        if final_path.exists() and part_path.exists():
            raise ValueError('both final and partial files exist for cluster %s; preserved for review' % ident)
        known = prior.get(ident)
        path = final_path if final_path.exists() else part_path if part_path.exists() else None
        if known is not None and path != final_path:
            raise ValueError('recorded completed FITS is missing for cluster %s' % ident)
        if path is None:
            continue
        try:
            review = validate_fits(path, frame)
        except (ValueError, OSError):
            if path != part_path:
                raise
            saved = part_path.with_name(part_path.name + '.incomplete-%d' % time.time_ns())
            part_path.rename(saved)
            status.setdefault('preserved_incomplete_files', []).append(str(saved))
            print('Preserved incomplete cluster %s response; it will be downloaded again.' % ident,
                  flush=True)
            continue
        size = path.stat().st_size
        digest = file_sha256(path)
        if known is not None and (size != known.get('bytes') or digest != known.get('sha256')):
            raise ValueError('recorded FITS size or SHA-256 differs for cluster %s' % ident)
        if path == part_path:
            os.replace(part_path, final_path)
        artifact = dict(known or {})
        artifact.update(id=ident, file=final_path.name, bytes=size, sha256=digest,
                        content_type=artifact.get('content_type', 'image/fits'), review=review)
        downloaded[ident] = artifact
        output_frames[ident] = build_manifest_row(frame, final_path.name, review)
        if known is None:
            status['completed'] = [downloaded[str(old['id'])]
                                   for old in frames if str(old['id']) in downloaded]
            status.update(status='downloading', bulk_download=False)
            status.pop('error', None)
            write_json(status_path, status)
            print('Recovered and verified cluster %s from the saved FITS response.' % ident,
                  flush=True)
    status['completed'] = [downloaded[str(frame['id'])]
                           for frame in frames if str(frame['id']) in downloaded]
    status.update(status='downloading', bulk_download=False)
    status.pop('error', None)
    write_json(status_path, status)
    return downloaded, output_frames


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', default=str(ROOT / 'pilot_2026-10-12/pilot_manifest.json'))
    parser.add_argument('--user', default=os.environ.get('HSC_USER'),
                        help='HSC account name (default: HSC_USER or interactive prompt)')
    parser.add_argument('--resume', action='store_true',
                        help='resume an incomplete work directory after revalidating saved FITS files')
    parser.add_argument('--work', required=True,
                        help='output directory; it must be new unless --resume is set')
    args = parser.parse_args(argv)
    work = Path(args.work).resolve()
    try:
        with Path(args.manifest).open(encoding='utf-8') as handle:
            doc = json.load(handle)
        frames = doc.get('clusters') or []
        pilot = doc.get('schema') == 'hsc-bcg-diffuse-pilot-v1' and 10 <= len(frames) <= 20
        sample = doc.get('schema') == 'hsc-bcg-diffuse-sample-v1' and 1 <= len(frames) <= 1000
        if not (pilot or sample):
            raise ValueError('manifest must be a prepared HSC pilot or selected sample (up to 1000 clusters)')
        if any(frame.get('image') for frame in frames):
            raise ValueError('manifest already contains image paths')
        ids = [str(frame.get('id')) for frame in frames]
        if len(set(ids)) != len(ids) or any(not re.fullmatch(r'[A-Za-z0-9_-]+', ident) for ident in ids):
            raise ValueError('cluster ids must be unique safe filename components')
        estimate, estimate_rows = estimate_bytes(frames)
        work.parent.mkdir(parents=True, exist_ok=True)
        usage = shutil.disk_usage(work.parent)
        if usage.free < int(estimate * 1.3):
            raise ValueError('free space is below the 30%% estimate margin (estimated %.2f GiB)' %
                             (estimate / 1024**3))
        status_path = work / 'download_status.json'
        images = work / 'images'
        manifest_digest = file_sha256(args.manifest)
        if args.resume:
            if not work.is_dir() or not images.is_dir() or not status_path.is_file():
                raise ValueError('--resume requires an existing work directory, images/, and download_status.json')
            if (work / 'download_manifest.json').exists():
                raise ValueError('download manifest already exists; completed work is not resumable')
            with status_path.open(encoding='utf-8') as handle:
                status = json.load(handle)
            if status.get('status') not in ('starting', 'downloading', 'stopped', 'interrupted'):
                raise ValueError('download status is not an incomplete run')
            if (status.get('n_expected') != len(frames)
                    or status.get('estimated_bytes') != estimate
                    or status.get('science_measurement') is not False):
                raise ValueError('resume manifest does not match the recorded download')
            if status.get('manifest_sha256') not in (None, manifest_digest):
                raise ValueError('resume manifest SHA-256 differs from the recorded input')
            if status.get('expected_ids') not in (None, ids):
                raise ValueError('resume cluster list differs from the recorded input')
            status.setdefault('manifest_sha256', manifest_digest)
            status.setdefault('expected_ids', ids)
            downloaded, output_frames = recover_existing(work, frames, status, status_path)
        else:
            work.mkdir(exist_ok=False)
            images.mkdir()
            status = {'status': 'starting', 'n_expected': len(frames), 'completed': [],
                      'estimated_bytes': estimate, 'manifest_sha256': manifest_digest,
                      'expected_ids': ids, 'science_measurement': False,
                      'credentials_saved': False, 'bulk_download': False}
            write_json(status_path, status)
            downloaded, output_frames = {}, {}

        pending = [frame for frame in frames if str(frame['id']) not in downloaded]
        if pending:
            manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
            username = (args.user or input('HSC account name: ')).strip()
            if not username:
                raise ValueError('HSC account name is required')
            password = getpass.getpass('HSC password (hidden): ')
            manager.add_password(None, 'https://hsc-release.mtk.nao.ac.jp/',
                                 username, password)
            if not password:
                raise ValueError('account name and password are required')
            opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(manager))
            print('Downloading %d remaining HSC-I coadd/bg frames (%d/%d already verified).' %
                  (len(pending), len(downloaded), len(frames)), flush=True)
        for i, frame in enumerate(frames, 1):
            if str(frame['id']) in downloaded:
                continue
            artifact, outframe = fetch_one(opener, frame, images, i, len(frames))
            ident = str(frame['id'])
            downloaded[ident] = artifact
            output_frames[ident] = outframe
            status['completed'] = [downloaded[str(old['id'])]
                                   for old in frames if str(old['id']) in downloaded]
            status['status'] = 'downloading'
            write_json(status_path, status)
        output_manifest = dict(doc)
        output_manifest.update(clusters=[output_frames[str(frame['id'])] for frame in frames],
                               downloaded_from=MANUAL,
                               downloaded_images=str(images), bulk_download=True,
                               science_measurement=False)
        with (work / 'download_manifest.json').open('x', encoding='utf-8') as handle:
            json.dump(output_manifest, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write('\n')
        status.update(status='complete', downloaded=len(downloaded), bulk_download=True,
                      frames_manifest='download_manifest.json',
                      extinction_values='pending', random_fields='pending',
                      science_measurement=False)
        write_json(status_path, status)
        print('Downloaded %d/%d verified FITS files to %s' % (len(downloaded), len(frames), work), flush=True)
        print('Next: review download_manifest.json; explicit Galactic extinction and matched random fields remain required.', flush=True)
        return 0
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as exc:
        if 'status_path' in locals() and 'status' in locals():
            status.update(status='stopped', error=str(exc), science_measurement=False)
            try:
                write_json(status_path, status)
            except OSError:
                pass
        parser.exit(2, 'download_pilot: %s\n' % exc)
    except (EOFError, KeyboardInterrupt):
        if 'status_path' in locals() and 'status' in locals():
            status.update(status='interrupted', science_measurement=False)
            try:
                write_json(status_path, status)
            except OSError:
                pass
        parser.exit(1, 'download_pilot: interrupted; completed files/status are preserved\n')


if __name__ == '__main__':
    sys.exit(main())
