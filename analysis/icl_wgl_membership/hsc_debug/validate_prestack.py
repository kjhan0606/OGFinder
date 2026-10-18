#!/usr/bin/env python3
"""Coadd-level diffuse injections with re-detection; full control profile inputs."""
import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import prepare_detection_psf as prep
from download_pilot import file_sha256 as sha, validate_fits
from validate_control_images import load_controls, image_path, nJy_factor, control_masks
from prestack_core import EDGES, FACTORS, MODELS, geometry, measure, diffuse_model


def surface_result(result, factor):
    out = dict(result)
    for key in ('raw', 'corrected', 'plane'):
        out[key] = result[key]*factor
    for key in ('nominal_raw', 'nominal_corrected'):
        out[key] = result[key]*factor**2
    return out


def vectors(results):
    return {key:np.array([r[key] for r in results]) for key in
            ('raw', 'corrected', 'coverage', 'count', 'nominal_raw', 'nominal_corrected', 'plane', 'n_tiles')}


def process(frame, target, record, base, cosmo, out):
    ident = frame['id']; path = image_path(base, frame); digest = sha(path)
    if digest != record['sha256']:
        raise ValueError('Control image changed: '+ident)
    info = validate_fits(path, frame); factor = nJy_factor(info)
    with fits.open(path, memmap=False) as h:
        image = h[info['science_hdu']].data.astype('f4')
        variance = h[info['variance_hdu']].data.astype('f4')
        bits = h[info['mask_hdu']].data
    center = info['center_pixel_xy']; scale = info['pixel_scale_arcsec']
    kpp = cosmo.kpc_proper_per_arcmin(frame['z']).value*scale/60
    geom = geometry(image.shape, center, kpp)
    _, masks, _, detections, detection = control_masks(image, variance, bits, center, scale)
    measured = [surface_result(measure(image, variance, mask, geom, center, kpp), factor) for mask in masks]
    payload = vectors(measured)
    payload.update(mask_bits=np.packbits(np.array(masks), axis=-1, bitorder='little'),
                   image_shape=np.array(image.shape), center=np.array(center), kpp=kpp, native_to_sb=factor)
    out.mkdir()
    np.savez_compressed(out/'baseline.npz', **payload)
    rows = []; source_count = []; mask_hashes = {}
    if target is not None:
        for index, (hscale, q, amplitude) in enumerate(MODELS):
            model = diffuse_model(geom[0], geom[1], hscale, q, amplitude)
            native_model = model/factor
            added = image.astype(float)+native_model
            _, after_masks, _, _, after_info = control_masks(added, variance, bits, center, scale)
            # Geometry truth has no source exclusions; no cosmological dimming/extinction applied.
            bins = geom[3]; good = bins >= 0
            truth = np.bincount(bins[good], weights=model[good], minlength=6)/np.bincount(bins[good], minlength=6)
            case = {'truth':truth, 'parameters':np.array([hscale, q, amplitude]),
                    'after_mask_bits':np.packbits(np.array(after_masks), axis=-1, bitorder='little')}
            for j, (kron, before_mask, after_mask) in enumerate(zip(FACTORS, masks, after_masks)):
                before = measured[j]
                post = surface_result(measure(added, variance, after_mask, geom, center, kpp), factor)
                fixed_signal = surface_result(measure(native_model, variance, before_mask, geom, center, kpp), factor)
                common_signal = surface_result(measure(native_model, variance, before_mask|after_mask, geom, center, kpp), factor)
                for name, result in [('after',post), ('fixed_signal',fixed_signal), ('common_signal',common_signal)]:
                    for key in ('raw','corrected','coverage','count','plane','n_tiles'):
                        case[f'{name}_{key}_{j}'] = result[key]
                for mode in ('raw', 'corrected'):
                    increment = post[mode]-before[mode]
                    for b in range(6):
                        eligible = bool(before['coverage'][b] >= .5 and post['coverage'][b] >= .5 and np.isfinite(increment[b]))
                        rows.append(dict(control_id=ident, target_cluster_id=frame['target_cluster_id'],
                            cohort=target['preparation_cohort'], model_index=index, scale_kpc=hscale,
                            axis_ratio=q, amplitude_at_100kpc=amplitude, kron_factor=kron, extraction=mode,
                            r_min_kpc=EDGES[b], r_max_kpc=EDGES[b+1], truth_njy_arcsec2=truth[b],
                            recovered_increment=increment[b], recovery_fraction=increment[b]/truth[b],
                            fixed_support_fraction=fixed_signal[mode][b]/truth[b],
                            common_support_fraction=common_signal[mode][b]/truth[b],
                            before_coverage=before['coverage'][b], after_coverage=post['coverage'][b],
                            area_50pct_eligible=eligible,
                            detections_before=len(detections), detections_after=after_info['n_detections']))
            case_path = out/f'injection_{index:02d}.npz'
            np.savez_compressed(case_path, **case); mask_hashes[str(case_path)] = sha(case_path)
            source_count.append(after_info['n_detections'])
            del model, native_model, added, after_masks
    metadata = dict(control_id=ident, target_cluster_id=frame['target_cluster_id'],
        source_image=str(path), source_sha256=digest, baseline_sha256=sha(out/'baseline.npz'),
        ra_deg=frame['ra_deg'], dec_deg=frame['dec_deg'], assigned_z=frame['z'],
        cohort=target['preparation_cohort'] if target else 'covariance_extension',
        frame_role='matched_control', bcg_policy='none', header_filter=info['coadd_header_filter'],
        n_detections=detection['n_detections'], injection_detections=source_count,
        injection_sha256=mask_hashes, n_injection_models=len(source_count),
        detection_background_subtracted_from_science=False, noise_added=False,
        science_measurement=False, sky_correction_approved=False)
    prep.write_json(out/'metadata.json', metadata)
    if rows:
        prep.table(out/'recovery.tsv', rows, list(rows[0]))
    return metadata


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', type=Path, required=True); p.add_argument('--work', type=Path, required=True)
    p.add_argument('--limit', type=int)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False); start = time.time()
    doc, download, base, hashes = load_controls(a.sample)
    vp = a.sample/'validation_sample_2026-10-13_01/validation_manifest.json'
    targets = {str(f['id']):f for f in json.loads(vp.read_text())['clusters']}
    hashes[str(vp)] = sha(vp)
    # Keep an immutable record of the exact source used, before processing starts.
    paths = [Path(__file__), Path(__file__).with_name('prestack_core.py'),
             Path(prep.__file__), Path(__file__).with_name('validate_control_images.py'),
             Path(__file__).with_name('download_pilot.py'), Path(__file__).with_name('hsc_quality.py')]
    scripts = {str(p):sha(p) for p in paths}
    prep.write_json(a.work/'protocol.json', dict(source_sha256=hashes, script_sha256=scripts,
        radial_edges_kpc=EDGES, kron_factors=FACTORS, models=MODELS, min_area_fraction=.5,
        outer_annulus_kpc=[600,900], covariance_seed=20261015, bootstrap_replicas=2000,
        engineering_limits=dict(max_absolute_median_recovery_bias=.10, min_eligible_injections=10,
                                min_covariance_n_over_p=5, min_sky_cells_over_p=2),
        injection_scope='coadd-level observed SB; deterministic diffuse additions; re-detection and source masks rerun',
        missing_validations=['upstream HSC sky transfer','PSF wing calibration','satellite/ICL joint transfer',
                             'source photon noise (gain/exposure mapping unavailable)','extinction/effective passband'],
        science_measurement=False, sky_correction_approved=False))
    frames = (sorted(doc['clusters'], key=lambda f: str(f['target_cluster_id']) not in targets)[:a.limit]
              if a.limit else doc['clusters'])
    cosmo = FlatLambdaCDM(H0=doc['h0'], Om0=doc['omega_m'], Tcmb0=0)
    records = []; failures = []
    for index, frame in enumerate(frames, 1):
        try:
            record = process(frame, targets.get(str(frame['target_cluster_id'])), download[frame['id']],
                             base, cosmo, a.work/frame['id'])
            records.append(record)
        except Exception as exc:
            failures.append(dict(control_id=frame['id'], error_type=type(exc).__name__, error=str(exc)))
            print('ERROR', frame['id'], type(exc).__name__, str(exc), flush=True)
        if index % 10 == 0 or index == len(frames):
            print(f'Processed {index}/{len(frames)}; injections={sum(r["n_injection_models"] for r in records)}; errors={len(failures)}', flush=True)
    prep.write_json(a.work/'extraction_summary.json', dict(status='complete' if not failures else 'partial',
        n_expected=len(frames), n_measured=len(records), n_injection_fields=sum(r['n_injection_models']>0 for r in records),
        n_injection_images=sum(r['n_injection_models'] for r in records), records=records, errors=failures,
        source_sha256=hashes, script_sha256=scripts, protocol_sha256=sha(a.work/'protocol.json'),
        elapsed_seconds=time.time()-start, science_measurement=False, sky_correction_approved=False))
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
