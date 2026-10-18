#!/usr/bin/env python3
"""Pixel-level technical QA for the downloaded HSC-I cluster cutouts.

Reports structural/data/mask/area diagnostics. It does not claim sky,
PSF-wing, random-field, completeness, or injection/recovery validation.
Run this over the full sample inside Slurm with hsc_quality_job.sh.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, validate_fits

REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from ogfkit.wlcluster import ARCSEC, Cosmology

MASK_BITS = {'BAD': 0, 'SAT': 1, 'INTRP': 2, 'CR': 3, 'EDGE': 4,
             'DETECTED': 5, 'NO_DATA': 8, 'BRIGHT_OBJECT': 9,
             'UNMASKEDNAN': 12, 'REJECTED': 13, 'CLIPPED': 14,
             'SENSOR_EDGE': 15, 'INEXACT_PSF': 16}
BAD_VALUE = sum(1 << bit for bit in BAD_BIT_POSITIONS)
RADIAL_EDGES_KPC = np.arange(0.0, 525.0, 25.0)


def write_json(path, value):
    temp = Path(str(path) + '.tmp-%d' % os.getpid())
    with temp.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


def summarize_frame(frame, base, recorded, cosmology):
    ident = str(frame['id'])
    image_name = frame.get('image')
    if not image_name:
        raise ValueError('manifest omits image filename')
    path = Path(image_name)
    if not path.is_absolute():
        path = Path(base) / path
    if not path.is_file():
        raise ValueError('image file is missing')
    actual_bytes = path.stat().st_size
    if recorded and actual_bytes != int(recorded['bytes']):
        raise ValueError('file size differs from download SHA record')

    header_review = validate_fits(path, frame)
    with fits.open(path, memmap=True) as hdul:
        sci = hdul[header_review['science_hdu']].data
        mask = hdul[header_review['mask_hdu']].data
        variance = hdul[header_review['variance_hdu']].data
        hdu_shapes = [list(hdu.shape) if hdu.data is not None else None for hdu in hdul]
        size = path.stat().st_size
        for hdu in hdul:
            info = hdu.fileinfo()
            if info and info['datLoc'] + info['datSpan'] > size:
                raise ValueError('FITS data plane extends past end of file')
            if 'CHECKSUM' in hdu.header and hdu.verify_checksum() != 1:
                raise ValueError('FITS CHECKSUM failed')
            if 'DATASUM' in hdu.header and hdu.verify_datasum() != 1:
                raise ValueError('FITS DATASUM failed')
        header = hdul[header_review['science_hdu']].header.copy()
        if header.get('INHERIT') is True:
            header.extend(hdul[0].header, unique=True)
        wcs = WCS(header).celestial

        image_finite = np.isfinite(sci)
        var_finite = np.isfinite(variance)
        var_positive = var_finite & (variance > 0)
        mask_finite = np.isfinite(mask)
        if not mask_finite.all():
            raise ValueError('MASK plane contains nonfinite values')
        integer_mask = mask.astype(np.int64, copy=False)
        bad = (integer_mask & BAD_VALUE) != 0
        bright = (integer_mask & BRIGHT_OBJECT_BIT) != 0
        detected = (integer_mask & (1 << MASK_BITS['DETECTED'])) != 0
        analysis_valid = image_finite & var_positive & ~bad & ~bright
        ny, nx = sci.shape
        cx, cy = header_review['center_pixel_xy']
        kpc_per_arcsec = float(cosmology.angular_diameter(float(frame['z']))) * ARCSEC * 1000.0
        pix_kpc = header_review['pixel_scale_arcsec'] * kpc_per_arcsec
        yy, xx = np.ogrid[:ny, :nx]
        radius = np.hypot((xx - cx) * pix_kpc, (yy - cy) * pix_kpc)
        radial_total, _ = np.histogram(radius.ravel(), RADIAL_EDGES_KPC)
        radial_valid, _ = np.histogram(radius[analysis_valid].ravel(), RADIAL_EDGES_KPC)
        coverage = np.divide(radial_valid, radial_total,
                             out=np.zeros_like(radial_valid, dtype=float),
                             where=radial_total > 0)

        outer = (radius >= 750) & (radius < 950) & analysis_valid & ~detected
        outer_values = np.asarray(sci[outer], dtype=np.float64)
        if outer_values.size >= 100:
            outer_median = float(np.median(outer_values))
            outer_mad = float(1.4826 * np.median(np.abs(outer_values - outer_median)))
            outer_sigma_variance = float(np.sqrt(np.median(variance[outer])))
            if outer_sigma_variance > 0:
                variance_ratio = outer_mad / outer_sigma_variance
            else:
                variance_ratio = None
        else:
            outer_median = outer_mad = outer_sigma_variance = variance_ratio = None

        # Outer-cell medians expose large spatial residuals. They are not a
        # blank-sky estimate: cluster light can extend into these radii.
        cell_medians = []
        for iy in range(6):
            y0, y1 = round(iy * ny / 6), round((iy + 1) * ny / 6)
            for ix in range(6):
                x0, x1 = round(ix * nx / 6), round((ix + 1) * nx / 6)
                block = sci[y0:y1, x0:x1]
                good = analysis_valid[y0:y1, x0:x1] & ~detected[y0:y1, x0:x1]
                if int(np.count_nonzero(good)) >= 50:
                    cell_medians.append(float(np.median(block[good])))
        if len(cell_medians) >= 4:
            cell_median = float(np.median(cell_medians))
            cell_span = float(np.percentile(cell_medians, 90) - np.percentile(cell_medians, 10))
            cell_mad = float(1.4826 * np.median(np.abs(np.asarray(cell_medians) - cell_median)))
        else:
            cell_median = cell_span = cell_mad = None

        core_radius = max(1, int(round(10.0 / pix_kpc)))
        core_cx, core_cy = int(round(cx)), int(round(cy))
        core_y0, core_y1 = max(0, core_cy - core_radius), min(ny, core_cy + core_radius + 1)
        core_x0, core_x1 = max(0, core_cx - core_radius), min(nx, core_cx + core_radius + 1)
        core_yy, core_xx = np.ogrid[core_y0:core_y1, core_x0:core_x1]
        core_circle = ((core_xx - core_cx) ** 2 + (core_yy - core_cy) ** 2
                       <= core_radius ** 2)
        core_mask = integer_mask[core_y0:core_y1, core_x0:core_x1]
        core_bad = bad[core_y0:core_y1, core_x0:core_x1]
        core_bright = bright[core_y0:core_y1, core_x0:core_x1]
        core_valid = analysis_valid[core_y0:core_y1, core_x0:core_x1]
        core_area = int(np.count_nonzero(core_circle))

        total = int(sci.size)
        fluxmag0 = float(hdul[0].header['FLUXMAG0'])
        area_arcsec2 = header_review['pixel_scale_arcsec'] ** 2
        nJy_per_pixel = 3.631e12 / fluxmag0
        to_njy_arcsec2 = nJy_per_pixel / area_arcsec2
        bit_fractions = {name: float(np.count_nonzero(integer_mask & (1 << bit)) / total)
                         for name, bit in MASK_BITS.items()}
        record = {
            'id': ident, 'name': frame.get('name'), 'path': str(path),
            'bytes': actual_bytes, 'sha256_from_download_record': recorded.get('sha256') if recorded else None,
            'fits_hdu_shapes': hdu_shapes, 'image_shape': [ny, nx],
            'header_filter': header_review['coadd_header_filter'],
            'center_pixel_xy': [float(cx), float(cy)],
            'pixel_scale_arcsec': header_review['pixel_scale_arcsec'],
            'wcs_scale_variation': header_review.get('wcs_scale_variation'),
            'pixel_scale_kpc': pix_kpc,
            'image_nonfinite_fraction': float(1 - np.count_nonzero(image_finite) / total),
            'variance_nonfinite_fraction': float(1 - np.count_nonzero(var_finite) / total),
            'variance_nonpositive_fraction': float(1 - np.count_nonzero(var_positive) / total),
            'mask_bit_fractions': bit_fractions,
            'quality_unmasked_fraction': float(np.count_nonzero(analysis_valid) / total),
            'center_pixel_has_detected_bit': bool(int(integer_mask[core_cy, core_cx]) & (1 << 5)),
            'center_pixel_has_bright_object_bit': bool(int(integer_mask[core_cy, core_cx]) & BRIGHT_OBJECT_BIT),
            'center_10kpc_radius_pixel_count': core_area,
            'center_10kpc_bright_object_mask_fraction': float(np.count_nonzero(core_circle & core_bright) / core_area),
            'center_10kpc_other_bad_mask_fraction': float(np.count_nonzero(core_circle & core_bad) / core_area),
            'center_10kpc_usable_fraction': float(np.count_nonzero(core_circle & core_valid) / core_area),
            'coverage_by_25kpc_annulus_0_500': coverage.tolist(),
            'first_annulus_below_50pct_coverage_kpc': (
                float(RADIAL_EDGES_KPC[np.flatnonzero(coverage < 0.5)[0]])
                if np.any(coverage < 0.5) else None),
            'outer_750_950_kpc_pixel_count': int(outer_values.size),
            'outer_750_950_kpc_median_instrumental_flux_per_pixel': outer_median,
            'outer_750_950_kpc_mad_sigma_instrumental_flux_per_pixel': outer_mad,
            'outer_750_950_kpc_median_variance_sigma': outer_sigma_variance,
            'outer_mad_to_variance_sigma_ratio': variance_ratio,
            'outer_grid_cell_count': len(cell_medians),
            'grid_cell_median_span_90_10_instrumental_flux_per_pixel': cell_span,
            'grid_cell_median_mad_sigma_instrumental_flux_per_pixel': cell_mad,
            'instrumental_flux_to_njy_arcsec2_factor': to_njy_arcsec2,
            'outer_median_njy_arcsec2_uncorrected_for_extinction': (
                outer_median * to_njy_arcsec2 if outer_median is not None else None),
            'science_measurement': False,
        }
        return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--status', required=True)
    parser.add_argument('--work', required=True)
    parser.add_argument('--allow-cpu-host',
                        help='permit a direct single-process run only when the hostname matches')
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        if not args.allow_cpu_host or not platform.node().lower().startswith(args.allow_cpu_host.lower()):
            parser.error('outside Slurm, direct execution requires an explicitly allowed matching host')
    manifest_path = Path(args.manifest).resolve()
    status_path = Path(args.status).resolve()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=False)
    doc = json.loads(manifest_path.read_text(encoding='utf-8'))
    state = json.loads(status_path.read_text(encoding='utf-8'))
    frames = doc.get('clusters') or []
    expected = {str(f['id']) for f in frames}
    if (doc.get('schema') != 'hsc-bcg-diffuse-sample-v1' or len(expected) != len(frames)
            or state.get('status') != 'complete' or state.get('n_expected') != len(frames)
            or state.get('downloaded') != len(frames) or doc.get('science_measurement') is not False):
        raise ValueError('download completion, sample schema, or IDs do not agree')
    records = {str(a['id']): a for a in state.get('completed', [])}
    if set(records) != expected:
        raise ValueError('download records do not cover all manifest IDs')
    image_base = Path(doc.get('downloaded_images', manifest_path.parent / 'images'))
    if not image_base.is_absolute():
        image_base = manifest_path.parent / image_base
    cosmology = Cosmology(H0=float(doc['h0']), Om=float(doc['omega_m']))

    details, failures = [], []
    for index, frame in enumerate(frames, 1):
        ident = str(frame['id'])
        print('[%d/%d] QA cluster %s' % (index, len(frames), ident), flush=True)
        try:
            details.append(summarize_frame(frame, image_base, records[ident], cosmology))
        except Exception as exc:
            failures.append({'id': ident, 'error_type': type(exc).__name__, 'reason': str(exc)})
    write_json(work / 'quality_frames.json', details)
    with (work / 'quality_frames.tsv').open('x', encoding='utf-8', newline='') as stream:
        fields = ['id', 'name', 'bytes', 'header_filter', 'image_shape',
                  'image_nonfinite_fraction', 'variance_nonfinite_fraction',
                  'variance_nonpositive_fraction', 'quality_unmasked_fraction',
                  'first_annulus_below_50pct_coverage_kpc',
                  'outer_750_950_kpc_pixel_count',
                  'outer_median_njy_arcsec2_uncorrected_for_extinction',
                  'outer_mad_to_variance_sigma_ratio', 'grid_cell_count',
                  'grid_cell_median_span_90_10_instrumental_flux_per_pixel']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in details:
            slim = dict(row)
            slim['image_shape'] = 'x'.join(map(str, row['image_shape']))
            slim['grid_cell_count'] = row['outer_grid_cell_count']
            writer.writerow(slim)

    with (work / 'quality_bcg_centres.tsv').open('x', encoding='utf-8', newline='') as stream:
        fields = ['id', 'name', 'center_pixel_has_detected_bit',
                  'center_pixel_has_bright_object_bit', 'center_10kpc_radius_pixel_count',
                  'center_10kpc_bright_object_mask_fraction',
                  'center_10kpc_other_bad_mask_fraction', 'center_10kpc_usable_fraction']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(details)

    radial_rows = []
    for i in range(len(RADIAL_EDGES_KPC) - 1):
        values = sorted(float(row['coverage_by_25kpc_annulus_0_500'][i]) for row in details)
        radial_rows.append({
            'r_min_kpc': float(RADIAL_EDGES_KPC[i]),
            'r_max_kpc': float(RADIAL_EDGES_KPC[i + 1]),
            'n_clusters_coverage_ge_0p50': sum(v >= 0.5 for v in values),
            'median_usable_fraction': float(np.median(values)),
            'p10_usable_fraction': values[int(0.1 * (len(values) - 1))],
        })
    with (work / 'quality_radial_coverage.tsv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(radial_rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(radial_rows)

    coverage = [r for r in details if r.get('first_annulus_below_50pct_coverage_kpc') is None
                or r['first_annulus_below_50pct_coverage_kpc'] >= 500]
    summary = {
        'status': 'complete' if not failures else 'partial',
        'n_requested': len(frames), 'n_frames_reported': len(details),
        'n_frame_errors': len(failures),
        'n_with_ge_50pct_usable_area_in_every_0_500kpc_bin': len(coverage),
        'n_clusters_coverage_ge_0p50_by_25kpc_bin_0_500': [r['n_clusters_coverage_ge_0p50'] for r in radial_rows],
        'minimum_n_clusters_coverage_ge_0p50_in_any_25kpc_bin_0_500': min(
            r['n_clusters_coverage_ge_0p50'] for r in radial_rows),
        'radial_coverage_table': 'quality_radial_coverage.tsv',
        'n_center_pixels_with_detected_bit': sum(r['center_pixel_has_detected_bit'] for r in details),
        'n_center_pixels_with_bright_object_bit': sum(r['center_pixel_has_bright_object_bit'] for r in details),
        'n_with_bright_object_mask_over_25pct_of_center_10kpc': sum(
            r['center_10kpc_bright_object_mask_fraction'] > 0.25 for r in details),
        'median_usable_area_fraction_in_center_10kpc': float(np.median(
            [r['center_10kpc_usable_fraction'] for r in details])),
        'bcg_center_table': 'quality_bcg_centres.tsv',
        'n_with_any_nonfinite_science_pixels': sum(r['image_nonfinite_fraction'] > 0 for r in details),
        'n_with_any_nonfinite_or_nonpositive_variance_pixels': sum(
            r['variance_nonpositive_fraction'] > 0 for r in details),
        'n_with_under_50pct_total_unmasked_area': sum(r['quality_unmasked_fraction'] < 0.5 for r in details),
        'filter_header_counts': {key: sum(r['header_filter'] == key for r in details)
                                 for key in sorted({r['header_filter'] for r in details})},
        'download_bytes': sum(r['bytes'] for r in details),
        'checks': {
            'fits_structure_and_wcs': 'reported_per_frame; errors indicate hard technical failures',
            'finite_pixels_variance_and_masks': 'measured_per_frame; inspect quantitative results; no unvalidated global rejection cut imposed',
            'radial_usable_area': 'measured per 25 kpc bin; 50% is the existing profile engineering floor, not a scientific completeness claim',
            'outer_background': 'spatial diagnostics only; 750-950 kpc may contain real cluster light and is not treated as blank sky',
            'variance_noise_agreement': 'reported as a pixel-level diagnostic; coadd noise correlation remains unmodeled',
            'galactic_extinction': 'pending: per-cluster extinction has not been attached',
            'matched_random_fields': 'pending: no random frames are in the download manifest',
            'mask_completeness_and_growth': 'pending: raw HSC flags measured; source masks and injected-source recovery not run',
            'psf_wings': 'pending: PSF-star/profile data not in these cutout HDUs',
            'sky_transfer': 'pending: requires blank/random-field and injection-based transfer test',
            'profile_agreement': 'pending: no science stack executed',
            'sample_selection': 'pending: candidate list selection/completeness not independently validated',
        },
        'science_measurement': False,
        'notes': ['Outer-cell gradients are not interpreted as a pure sky error.',
                  'Technical QA metrics do not certify scientific quality.'],
        'frame_errors': failures,
    }
    write_json(work / 'quality_summary.json', summary)
    print('QA %s: %d/%d frames; %d errors; %d satisfy the 50%% per-bin area floor to 500 kpc.' %
          (summary['status'], len(details), len(frames), len(failures), len(coverage)), flush=True)
    print('Report: %s' % work, flush=True)
    return 0 if not failures else 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print('quality QA interrupted; partial frame report is preserved.', file=sys.stderr)
        sys.exit(130)
    except Exception as exc:
        print('hsc_quality: %s: %s' % (type(exc).__name__, exc), file=sys.stderr)
        sys.exit(2)
