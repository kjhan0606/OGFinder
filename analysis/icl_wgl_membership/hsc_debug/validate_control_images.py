#!/usr/bin/env python3
"""Offline whole-image control QA and paired outer-field diagnostics on lageunha.

The legacy download manifest calls its list 'clusters'; these are control fields,
with the target redshift used only for angular scale. No BCG is kept in controls.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import platform
import time

import astropy
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sep

import hsc_quality as quality
import prepare_detection_psf as prep
from diagnose_outer_sky import tile_diagnostics
from download_pilot import file_sha256 as sha, validate_fits


class ExactCosmology:
    """Keep the existing pixel-QA interface, using the current Astropy geometry."""
    def __init__(self, h0, omega_m):
        self.model = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)

    def angular_diameter(self, z):
        return self.model.angular_diameter_distance(z).value


def read_json(path):
    return json.loads(path.read_text())


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_table(path, rows):
    if rows:
        prep.table(path, rows, list(dict.fromkeys(k for row in rows for k in row)))


def sources(paths):
    return {str(path.resolve()): sha(path) for path in paths}


def script_hashes():
    here = Path(__file__).resolve().parent
    return sources([Path(__file__), here/'hsc_quality.py', here/'download_pilot.py',
                    here/'prepare_detection_psf.py', here/'diagnose_outer_sky.py'])


def load_controls(sample):
    mp = sample/'control_download_01/download_manifest.json'
    sp = sample/'control_download_01/download_status.json'
    doc, state = read_json(mp), read_json(sp)
    frames = doc['clusters']
    ids = [str(f['id']) for f in frames]
    targets = [str(f['target_cluster_id']) for f in frames]
    records = {str(r['id']): r for r in state['completed']}
    if (not frames or len(set(ids)) != len(ids) or len(set(targets)) != len(targets)
            or set(records) != set(ids) or len(state['completed']) != len(frames)
            or state.get('status') != 'complete' or state.get('downloaded') != len(frames)
            or state.get('n_expected') != len(frames)
            or doc.get('science_measurement') is not False
            or doc.get('frame_role') != 'matched_control' or doc.get('bcg_policy') != 'none'):
        raise ValueError('Control identity/download completion mismatch')
    for f in frames:
        if (f.get('frame_role') != 'matched_control' or f.get('bcg_policy') != 'none'
                or f.get('members_radec_deg') != [] or f.get('randoms') != []
                or f.get('science_measurement') is not False):
            raise ValueError('Control frame must have no BCG/member/science activation')
    base = Path(doc.get('downloaded_images', mp.parent/'images'))
    if not base.is_absolute():
        base = mp.parent/base
    return doc, records, base, sources([mp, sp])


def image_path(base, frame):
    p = Path(frame['image'])
    return p if p.is_absolute() else base/p


def nJy_factor(info):
    """Native flux/pixel to nJy/arcsec^2, before extinction/passband correction."""
    unit = info['input_unit']
    if unit == 'instrumental_flux/pixel':
        per_pixel = 3.631e12/info['fluxmag0']
    elif unit == 'nJy/pixel':
        per_pixel = 1.
    else:
        raise ValueError('Unknown input unit: '+unit)
    return per_pixel/info['pixel_scale_arcsec']**2


def run_quality(sample, work):
    started = time.time()
    doc, download, base, hashes = load_controls(sample)
    cosmo = ExactCosmology(doc['h0'], doc['omega_m'])
    rows, errors = [], []
    with (work/'quality_frames.jsonl').open('x') as stream:
        for index, frame in enumerate(doc['clusters'], 1):
            ident = str(frame['id'])
            try:
                path = image_path(base, frame)
                digest = sha(path)
                if digest != download[ident]['sha256']:
                    raise ValueError('Full-file SHA-256 differs from download record')
                # Existing QA reports fluxmag0-scaled native flux. Reject a different
                # unit here rather than silently mislabeling that legacy diagnostic.
                if download[ident]['review']['input_unit'] != 'instrumental_flux/pixel':
                    raise ValueError('Legacy pixel QA requires instrumental_flux/pixel')
                row = quality.summarize_frame(frame, base, download[ident], cosmo)
                row.update(frame_role='matched_control', bcg_policy='none',
                           target_cluster_id=str(frame['target_cluster_id']),
                           random_object_id=str(frame['random_object_id']),
                           assigned_target_redshift=frame['z'],
                           redshift_policy=frame['redshift_policy'], sha256_rechecked=digest)
                rows.append(row)
                stream.write(json.dumps(row, allow_nan=False)+'\n'); stream.flush()
            except Exception as exc:
                errors.append(dict(id=ident, error_type=type(exc).__name__, error=str(exc)))
            if index % 20 == 0 or index == len(doc['clusters']):
                print(f'QA {index}/{len(doc["clusters"])}; errors={len(errors)}', flush=True)
    prep.write_json(work/'quality_frames.json', rows)
    radial = []
    for i in range(20):
        values = [r['coverage_by_25kpc_annulus_0_500'][i] for r in rows]
        radial.append(dict(r_min_kpc=i*25, r_max_kpc=(i+1)*25,
                           n_controls_ge_50pct=sum(v >= .5 for v in values),
                           median_valid_fraction=float(np.median(values)) if values else None))
    write_table(work/'radial_coverage.tsv', radial)
    summary = dict(status='complete' if not errors else 'partial', n_requested=len(doc['clusters']),
                   n_reported=len(rows), n_errors=len(errors), errors=errors,
                   n_full_file_hashes_matched=len(rows), bytes_rechecked=sum(r['bytes'] for r in rows),
                   filter_header_counts=dict(Counter(r['header_filter'] for r in rows)),
                   n_nonfinite_science_frames=sum(r['image_nonfinite_fraction'] > 0 for r in rows),
                   n_invalid_variance_frames=sum(r['variance_nonpositive_fraction'] > 0 for r in rows),
                   n_bright_mask_at_random_center=sum(r['center_pixel_has_bright_object_bit'] for r in rows),
                   n_under_50pct_total_quality_area=sum(r['quality_unmasked_fraction'] < .5 for r in rows),
                   n_ge_50pct_in_every_0_500kpc_bin=sum(min(r['coverage_by_25kpc_annulus_0_500']) >= .5 for r in rows),
                   source_sha256=hashes, script_sha256=script_hashes(), elapsed_seconds=time.time()-started,
                   software=dict(numpy=np.__version__, astropy=astropy.__version__, sep=sep.__version__),
                   cosmology=dict(H0=doc['h0'], Om0=doc['omega_m'], Tcmb0=0),
                   science_measurement=False, controls_ready_for_science=False,
                   limitations=[
                       'Quality-area masks retain all native bad bits and BRIGHT_OBJECT; source-mask growth is a separate stage.',
                       'Control centers are random positions, never BCGs. Redshift is assigned from the paired target for angular scale only.',
                       'The 50% annulus area floor is diagnostic, not sample completeness or scientific acceptance.',
                       'Legacy 750-950 kpc and grid medians exclude native DETECTED pixels only; residual sources and PSF wings remain.',
                       'No sky subtraction, extinction correction, PSF approval, or science stack.'])
    prep.write_json(work/'summary.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('source_sha256','script_sha256','limitations')}), flush=True)
    return 1 if errors else 0


def control_masks(image, variance, bits, center, scale):
    bad = ((bits & quality.BAD_VALUE) != 0) | ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    base = bad | ((bits & prep.BRIGHT_OBJECT_BIT) != 0)
    empty = np.empty((0, 2))
    objects, seg, kron, flags, _, snr, extra = prep.detect(image, variance, bad, center, empty, scale)
    # A source coinciding with the random center is still a contaminant, not a BCG.
    masks, fallback = prep.source_masks(objects, seg, kron, flags, base, -1, empty, scale)
    if fallback:
        raise ValueError('Control has unexpected member fallback')
    last = base
    for mask in masks:
        if np.any(last & ~mask) or np.any((seg > 0) & ~mask):
            raise ValueError('Control mask lost a defect or retained a detected source')
        last = mask
    detections = [dict(number=i+1, x=o['x'], y=o['y'], a=o['a'], b=o['b'], theta=o['theta'],
                       kron=kron[i], kron_flags=flags[i], detection_flags=o['flag'], nominal_snr=snr[i])
                  for i, o in enumerate(objects)]
    return base, masks, seg, detections, dict(
        background_global_for_detection_only=extra['background_global'],
        background_rms_for_detection_only=extra['background_rms'],
        n_detections=len(objects), center_segment=int(seg[round(center[1]), round(center[0])]))


def paired_row(control, cluster, control_factor, cluster_factor, same_filter):
    def value(r, key):
        x = r.get(key)
        return float(x) if x not in (None, '') else None
    cmean, tmean = value(control, 'median_tile_mean'), value(cluster, 'median_tile_mean')
    return dict(target_cluster_id=str(cluster['cluster_id']), control_id=control['control_id'],
                cohort=cluster['cohort'], kron_factor=float(cluster['kron_factor']), same_header_filter=same_filter,
                control_status=control['status'], cluster_status=cluster['status'],
                control_n_tiles=control['n_tiles'], cluster_n_tiles=int(cluster['n_tiles']),
                control_outer_valid_fraction=control['outer_valid_fraction'],
                cluster_outer_valid_fraction=float(cluster['outer_valid_fraction']),
                control_residual_ratio=value(control, 'plane_residual_mad_over_nominal'),
                cluster_residual_ratio=value(cluster, 'plane_residual_mad_over_nominal'),
                control_median_tile_mean_njy_arcsec2=cmean*control_factor if cmean is not None else None,
                cluster_median_tile_mean_njy_arcsec2=tmean*cluster_factor if tmean is not None else None,
                cluster_minus_control_median_tile_mean_njy_arcsec2=(
                    tmean*cluster_factor-cmean*control_factor if cmean is not None and tmean is not None else None))


def run_compare(sample, quality_dir, work):
    started = time.time()
    doc, download, base, hashes = load_controls(sample)
    qp = quality_dir/'summary.json'; qfp = quality_dir/'quality_frames.json'
    qs = read_json(qp); qa = {r['id']: r for r in read_json(qfp)}
    if qs['status'] != 'complete' or qs['n_reported'] != len(doc['clusters']) or qs['source_sha256'] != hashes:
        raise ValueError('Whole-image QA is incomplete or belongs to different inputs')
    vp = sample/'validation_sample_2026-10-13_01/validation_manifest.json'
    cp = sample/'download_01/download_status.json'
    outer = sample/'outer_sky_diagnostics_2026-10-14_02'
    op = outer/'summary.json'; fp = outer/'frame_stats.tsv'
    validation, old = read_json(vp), read_json(op)
    if (old['status'] != 'complete' or old['source_sha256'][str(vp)] != sha(vp)
            or old['source_sha256'][str(cp)] != sha(cp)
            or old['settings']['annulus_kpc'] != [600, 900]
            or old['settings']['tile_pixels'] != 64
            or old['settings']['H0'] != doc['h0'] or old['settings']['Om0'] != doc['omega_m']):
        raise ValueError('Cluster diagnostics provenance/settings mismatch')
    cluster_download = {str(r['id']): r for r in read_json(cp)['completed']}
    cluster_rows = read_table(fp)
    expected = {(str(f['id']), k) for f in validation['clusters'] for k in prep.FACTORS}
    lookup = {(r['cluster_id'], float(r['kron_factor'])): r for r in cluster_rows}
    if set(lookup) != expected or len(lookup) != len(cluster_rows):
        raise ValueError('Cluster diagnostic rows are incomplete or duplicated')
    hashes.update(sources([qp, qfp, vp, cp, op, fp]))
    by_target = {str(f['target_cluster_id']): f for f in doc['clusters']}
    cosmo = ExactCosmology(doc['h0'], doc['omega_m']).model
    records, tiles_all, lags_all, paired, masks_hashes, errors, missing = [], [], [], [], {}, [], []
    for index, target in enumerate(validation['clusters'], 1):
        cid = str(target['id'])
        if cid not in by_target:
            missing.append(cid)
            continue
        frame = by_target[cid]; ident = str(frame['id'])
        try:
            if (frame['z'] != target['z'] or
                    frame['requested_half_width_arcsec'] != target['requested_half_width_arcsec']):
                raise ValueError('Paired angular scale mismatch')
            path = image_path(base, frame); digest = sha(path)
            if digest != download[ident]['sha256'] or digest != qa[ident]['sha256_rechecked']:
                raise ValueError('Control changed since QA')
            hashes[str(path)] = digest
            info = validate_fits(path, frame)
            with fits.open(path, memmap=False) as hdus:
                image = np.asarray(hdus[info['science_hdu']].data, dtype='f4')
                variance = np.asarray(hdus[info['variance_hdu']].data, dtype='f4')
                bits = hdus[info['mask_hdu']].data
                header = WCS(hdus[info['science_hdu']].header).celestial.to_header()
            center = info['center_pixel_xy']; scale = info['pixel_scale_arcsec']
            kpp = float(cosmo.kpc_proper_per_arcmin(frame['z']).value)*scale/60
            native, masks, seg, det, detection = control_masks(image, variance, bits, center, scale)
            out = work/ident; out.mkdir()
            write_table(out/'detections.tsv', det)
            header['SCIVAL'] = False; header['BCGPOL'] = 'none'; header['BGSUB'] = False
            hdus = [fits.PrimaryHDU(), fits.ImageHDU(seg.astype('i4'), header=header, name='SEGMENTATION'),
                    fits.ImageHDU(native.astype('u1'), header=header, name='BASE_MASK')]
            for factor, mask in zip(prep.FACTORS, masks):
                hdus.append(fits.ImageHDU(mask.astype('u1'), header=header, name='MASK_K'+str(factor).replace('.', 'P')))
            mask_path = out/'masks.fits.gz'; fits.HDUList(hdus).writeto(mask_path)
            masks_hashes[str(mask_path)] = sha(mask_path)
            prep.write_json(out/'preparation.json', dict(detection, source_image=str(path), source_sha256=digest,
                control_id=ident, target_cluster_id=cid, bcg_policy='none', science_measurement=False,
                detection_background_subtracted_from_science=False, bright_mask_retained=True))
            yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
            radius = np.hypot(xx-center[0], yy-center[1])*kpp
            region = (radius >= 600) & (radius < 900)
            if not region.any():
                raise ValueError('Outer annulus has no pixels')
            target_info = cluster_download[cid]['review']
            for factor, mask in zip(prep.FACTORS, masks):
                tiles, lags, result = tile_diagnostics(image, variance, mask, region, center, kpp)
                common = dict(control_id=ident, target_cluster_id=cid, cohort=target['preparation_cohort'], kron_factor=factor)
                record = dict(common, **result, kpc_per_pixel=kpp, **detection)
                records.append(record)
                tiles_all.extend(dict(common, **r) for r in tiles)
                lags_all.extend(dict(common, **r) for r in lags)
                cluster = lookup[(cid, factor)]
                if cluster['cohort'] != target['preparation_cohort']:
                    raise ValueError('Paired cohort mismatch')
                paired.append(paired_row(record, cluster, nJy_factor(info), nJy_factor(target_info),
                                         info['coadd_header_filter'] == target_info['coadd_header_filter']))
            print(f'Paired outer diagnostics {index}/{len(validation["clusters"])}: {cid}', flush=True)
        except Exception as exc:
            errors.append(dict(control_id=ident, error_type=type(exc).__name__, error=str(exc)))
            print(f'ERROR {cid}: {exc}', flush=True)
    for name, rows in [('frame_stats', records), ('tile_stats', tiles_all), ('lag_products', lags_all), ('paired_stats', paired)]:
        write_table(work/(name+'.tsv'), rows)
    groups = []
    for cohort in sorted({p['cohort'] for p in paired}):
        for factor in prep.FACTORS:
            for same_only in (False, True):
                subset = [p for p in paired if p['cohort'] == cohort and p['kron_factor'] == factor
                          and (not same_only or p['same_header_filter'])]
                good = [p for p in subset if p['control_status'] == p['cluster_status'] == 'diagnostic_only']
                def median(key):
                    return float(np.median([p[key] for p in good])) if good else None
                groups.append(dict(cohort=cohort, kron_factor=factor, same_header_filter_only=same_only,
                    n_pairs=len(subset), n_pairs_with_ten_tiles_each=len(good),
                    median_cluster_coverage=median('cluster_outer_valid_fraction'),
                    median_control_coverage=median('control_outer_valid_fraction'),
                    median_cluster_residual_ratio=median('cluster_residual_ratio'),
                    median_control_residual_ratio=median('control_residual_ratio'),
                    median_paired_mean_difference_njy_arcsec2=median('cluster_minus_control_median_tile_mean_njy_arcsec2')))
    prep.write_json(work/'summary.json', dict(status='complete' if not errors and not missing else 'partial',
        n_target_frames=len(validation['clusters']), n_control_frames=len({r['control_id'] for r in records}),
        n_method_rows=len(records), n_paired_rows=len(paired), n_tiles=len(tiles_all), n_lag_rows=len(lags_all),
        n_pairs_with_header_filter_difference=len({r['control_id'] for r in paired if not r['same_header_filter']}),
        missing_target_ids=missing, errors=errors, groups=groups, source_sha256=hashes,
        output_masks_sha256=masks_hashes, script_sha256=script_hashes(), elapsed_seconds=time.time()-started,
        software=dict(numpy=np.__version__, astropy=astropy.__version__, sep=sep.__version__),
        settings=old['settings'], science_measurement=False, sky_correction_approved=False,
        controls_ready_for_science=False,
        limitations=[
            'Paired 40 area-qualified and 10 mask-stress fields are a selected validation sample, not representative of all 800 targets.',
            'All control detections are masked, including sources at the random center. No target members or BCG are invented.',
            'Cluster masks additionally protect the associated BCG and mask known member positions; this difference is retained explicitly.',
            'Cluster statistics reuse the hashed, previously audited outer_sky_diagnostics_2026-10-14_02 tables.',
            'Kron growth changes retained pixels/tiles; this is not a fixed-footprint causal comparison.',
            'Residual means and correlations include unmasked sources, PSF wings, cirrus and real cluster light; not a pure noise kernel.',
            'nJy/arcsec2 uses native FLUXMAG0 and pixel scale only, without extinction or effective passband correction.',
            'Header i/i2 agreement does not establish an identical effective coadd passband.',
            'No control subtraction, uncertainty rescaling, PSF approval, injection recovery or science stack.']))
    print(json.dumps(groups), flush=True)
    return 1 if errors or missing else 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['quality', 'compare'])
    p.add_argument('--sample', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--quality', type=Path)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    if a.stage == 'compare' and a.quality is None:
        p.error('compare requires --quality')
    a.work.mkdir(parents=True, exist_ok=False)
    return run_quality(a.sample.resolve(), a.work.resolve()) if a.stage == 'quality' else run_compare(
        a.sample.resolve(), a.quality.resolve(), a.work.resolve())


if __name__ == '__main__':
    raise SystemExit(main())
