#!/usr/bin/env python3
"""Observed outer-field residual diagnostics, not clean-sky or ICL measurement."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import platform
import time

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM

from fetch_ancillary_catalogues import sha
import prepare_detection_psf as prep

LAGS = (1, 2, 4, 8)


def tile_diagnostics(image, variance, excluded, region, center, kpp, tile_size=64):
    valid = region & ~excluded & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    tiles = []; sums = {(lag, axis):np.zeros(6) for lag in LAGS for axis in (0, 1)}
    for y0 in range(0, image.shape[0], tile_size):
        for x0 in range(0, image.shape[1], tile_size):
            sl = np.s_[y0:y0+tile_size, x0:x0+tile_size]
            good = valid[sl]; area = int(region[sl].sum()); n = int(good.sum())
            if area < 128 or n < 128 or n < .5*area:
                continue
            pix = image[sl].astype(float); var = variance[sl].astype(float)
            median = float(np.median(pix[good])); mean = float(pix[good].mean())
            yy, xx = np.nonzero(good)
            tiles.append(dict(tile_x=x0, tile_y=y0, n_valid=n, n_region=area, valid_fraction=n/area,
                              x_kpc=float((x0+xx.mean()-center[0])*kpp),
                              y_kpc=float((y0+yy.mean()-center[1])*kpp),
                              raw_mean=mean, raw_median=median, nominal_mean_sigma=float(np.sqrt(var[good].sum())/n)))
            residual = np.zeros(pix.shape)
            residual[good] = (pix[good]-median)/np.sqrt(var[good])
            # Pairs stay within the same accepted tile; no padded/periodic boundary pairs.
            for lag in LAGS:
                for axis in (0, 1):
                    a, b = ((np.s_[lag:, :], np.s_[:-lag, :]) if axis == 0 else
                            (np.s_[:, lag:], np.s_[:, :-lag]))
                    keep = good[a] & good[b]; x = residual[a][keep]; y = residual[b][keep]
                    sums[(lag, axis)] += (len(x), x.sum(), y.sum(), x@x, y@y, x@y)
    lag_rows = []
    for (lag, axis), (n, sx, sy, sxx, syy, sxy) in sums.items():
        vx = sxx/n-(sx/n)**2 if n else 0.; vy = syy/n-(sy/n)**2 if n else 0.
        covariance = sxy/n-sx*sy/n**2 if n else None
        lag_rows.append(dict(lag_pixels=lag, axis='y' if axis==0 else 'x', n_pairs=int(n),
                             variance_normalized_covariance=covariance,
                             correlation=float(covariance/np.sqrt(vx*vy)) if vx>0 and vy>0 else None))
    result = dict(status='insufficient_tiles', n_tiles=len(tiles), outer_valid_fraction=float(valid.sum()/region.sum()),
                  n_outer_valid=int(valid.sum()), n_outer_region=int(region.sum()),
                  median_tile_mean=None, median_nominal_mean_sigma=None,
                  plane_residual_mad=None, plane_residual_mad_over_nominal=None)
    if len(tiles) >= 10:
        matrix = np.array([[1., t['x_kpc']/1000, t['y_kpc']/1000] for t in tiles])
        means = np.array([t['raw_mean'] for t in tiles]); sigma = np.array([t['nominal_mean_sigma'] for t in tiles])
        coef, _, rank, _ = np.linalg.lstsq(matrix/sigma[:, None], means/sigma, rcond=None)
        if rank == 3:
            residual = means-matrix@coef
            mad = float(1.4826*np.median(abs(residual-np.median(residual))))
            result.update(status='diagnostic_only', median_tile_mean=float(np.median(means)),
                          median_nominal_mean_sigma=float(np.median(sigma)),
                          plane_residual_mad=mad, plane_residual_mad_over_nominal=mad/float(np.median(sigma)),
                          plane_constant=float(coef[0]), plane_x_per_mpc=float(coef[1]), plane_y_per_mpc=float(coef[2]))
            for t, v in zip(tiles, residual):
                t['plane_residual'] = float(v)
    return tiles, lag_rows, result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('sample', 'work'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--limit', type=int)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False); started = time.time()
    manifest_path = a.sample/'validation_sample_2026-10-13_01/validation_manifest.json'
    manifest = json.loads(manifest_path.read_text()); frames = manifest['clusters'][:a.limit]
    status_path = a.sample/'download_01/download_status.json'
    download = {str(r['id']):r for r in json.loads(status_path.read_text())['completed']}
    masks_path = a.sample/'bcg_masks_2026-10-13_01'
    input_hashes = {str(manifest_path):sha(manifest_path), str(status_path):sha(status_path),
                    str(masks_path/'summary.json'):sha(masks_path/'summary.json')}
    cosmo = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'])
    tile_rows = []; lag_rows = []; records = []
    for frame in frames:
        cid = str(frame['id']); image_path = Path(frame['image']); digest = sha(image_path)
        if digest != download[cid]['sha256']:
            raise ValueError('Source image changed: '+cid)
        input_hashes[str(image_path)] = digest
        mask_path = masks_path/('cluster_'+cid)/'proposed_masks.fits.gz'
        input_hashes[str(mask_path)] = sha(mask_path)
        info = prep.validate_fits(image_path, frame)
        with fits.open(image_path, memmap=False) as h:
            image = np.asarray(h[info['science_hdu']].data, dtype='f4')
            variance = np.asarray(h[info['variance_hdu']].data, dtype='f4')
            bits = h[info['mask_hdu']].data
            bad = ((bits & (sum(1<<b for b in prep.BAD_BIT_POSITIONS)|prep.BRIGHT_OBJECT_BIT)) != 0)
            image_wcs = WCS(h[info['science_hdu']].header).celestial
        center = info['center_pixel_xy']; kpp = float(cosmo.kpc_proper_per_arcmin(frame['z']).value)*info['pixel_scale_arcsec']/60
        yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
        radius = np.hypot(xx-center[0], yy-center[1])*kpp
        region = (radius >= 600) & (radius < 900)
        native = region & ~bad & np.isfinite(variance) & (variance>0) & np.isfinite(image)
        pixel_sigma = float(np.sqrt(np.median(variance[native]))) if native.any() else None
        with fits.open(mask_path, memmap=False) as mh:
            base = mh['BASE_MASK'].data.astype(bool)
            if base.shape != image.shape or np.any(bad & ~base):
                raise ValueError('Quality mask not retained')
            points = np.array([[0.,0.],[image.shape[1]-1,0.],[0.,image.shape[0]-1],center])
            world = image_wcs.all_pix2world(points, 0)
            last = base
            for factor in prep.FACTORS:
                hdu = mh['MASK_K'+str(factor).replace('.', 'P')]
                projected = WCS(hdu.header).celestial.all_world2pix(world, 0)
                if not np.allclose(points, projected, atol=.01, rtol=0):
                    raise ValueError('Mask WCS mismatch')
                mask = hdu.data.astype(bool)
                if mask.shape != image.shape or np.any(last & ~mask):
                    raise ValueError('Mask nesting/base preservation failed')
                last = mask
                tiles, lags, result = tile_diagnostics(image, variance, mask, region, center, kpp)
                common = dict(cluster_id=cid, cohort=frame['preparation_cohort'], kron_factor=factor)
                tile_rows.extend(dict(common, **r) for r in tiles); lag_rows.extend(dict(common, **r) for r in lags)
                records.append(dict(common, **result, kpc_per_pixel=kpp, outer_pixel_noise_sigma=pixel_sigma,
                                    median_tile_mean_over_pixel_sigma=result['median_tile_mean']/pixel_sigma if pixel_sigma and result['median_tile_mean'] is not None else None,
                                    plane_residual_mad_over_pixel_sigma=result['plane_residual_mad']/pixel_sigma if pixel_sigma and result['plane_residual_mad'] is not None else None,
                                    science_measurement=False, sky_correction_approved=False))
        print(cid, 'outer residual diagnostics complete', flush=True)
    for filename, rows in (('tile_stats.tsv', tile_rows), ('lag_products.tsv', lag_rows), ('frame_stats.tsv', records)):
        if rows:
            prep.table(a.work/filename, rows, list(dict.fromkeys(k for r in rows for k in r)))
    grouped = defaultdict(list)
    for r in records:
        grouped[(r['cohort'], r['kron_factor'])].append(r)
    cohorts = []
    for (cohort, factor), rows in grouped.items():
        good = [r for r in rows if r['status']=='diagnostic_only']
        chosen_lags = [r['correlation'] for r in lag_rows if r['cohort']==cohort and r['kron_factor']==factor and r['lag_pixels']==1 and r['correlation'] is not None]
        cohorts.append(dict(cohort=cohort, kron_factor=factor, n_frames=len(rows), n_with_ten_tiles=len(good),
                            median_outer_valid_fraction=float(np.median([r['outer_valid_fraction'] for r in rows])),
                            median_plane_residual_mad_over_nominal=float(np.median([r['plane_residual_mad_over_nominal'] for r in good])) if good else None,
                            median_lag1_correlation=float(np.median(chosen_lags)) if chosen_lags else None))
    prep.write_json(a.work/'summary.json', dict(status='complete', n_frames=len(frames), n_method_rows=len(records),
        n_tiles=len(tile_rows), n_lag_rows=len(lag_rows), elapsed_seconds=time.time()-started,
        source_sha256=input_hashes, script_sha256=sha(Path(__file__)), cohorts=cohorts,
        settings=dict(annulus_kpc=[600,900], tile_pixels=64, tile_min_pixels=128, tile_min_fraction=.5, lags_pixels=LAGS,
                      cosmology='astropy FlatLambdaCDM, Tcmb0=0', H0=manifest['h0'], Om0=manifest['omega_m']),
        science_measurement=False, sky_correction_approved=False, observational_validation_passed=False,
        limitations=['Cluster-centered outer fields can contain cluster light, PSF wings, cirrus and unmasked sources.',
                     'Residual means/plane fits are diagnostics, never subtracted from science images or activated as sky priors.',
                     'Correlation and aperture excess include astrophysical residuals and masks; not a calibrated coadd noise kernel or uncertainty multiplier.',
                     'Lag products use within-tile pairs after subtracting each tile median and dividing by per-pixel nominal sigma.',
                     'Threshold-selected tiles and retained frames are reported with counts; no representative 800-cluster claim.',
                     'Patch-level random sky statistics are not substituted for local image residuals.']))
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for cohort, color in (('area_qualified','#286a9d'), ('mask_stress','#c06b20')):
        xs = []; ys = []; cov = []
        for factor in prep.FACTORS:
            rows = [r for r in records if r['cohort']==cohort and r['kron_factor']==factor and r['status']=='diagnostic_only']
            xs.extend([factor]*len(rows)); ys.extend([r['plane_residual_mad_over_nominal'] for r in rows]); cov.extend([r['outer_valid_fraction'] for r in rows])
        axes[0].scatter(xs, ys, alpha=.6, s=15, label=cohort, color=color)
        axes[1].scatter(cov, ys, alpha=.6, s=15, label=cohort, color=color)
    axes[0].set_xlabel('Kron mask factor'); axes[1].set_xlabel('Outer annulus retained fraction')
    for ax in axes:
        ax.set_yscale('log'); ax.set_ylabel('Tile residual MAD / nominal mean sigma'); ax.axhline(1, color='k', lw=.5)
    axes[0].legend(fontsize=8)
    fig.suptitle('Observed 600-900 kpc residuals; cluster light may remain\nDiagnostic only, no sky subtraction or calibrated noise factor', fontsize=10)
    fig.tight_layout(); fig.savefig(a.work/'outer_residuals.png', dpi=140); plt.close(fig)
    print(json.dumps(cohorts))


if __name__ == '__main__':
    main()
