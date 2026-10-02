#!/usr/bin/env python3
"""Spatially varying PSF model from field stars (PSFEx-like) - CLI of the `psfex` plugin; engine: ogfkit/psfmodel.py.

    psfex.py IMAGE --work DIR [--catalog CAT] [--mask M] [--kind empirical|gaussian|moffat|lorentz|penny] [--order 2] [--oversample auto|1|2|3]
             [--size 0] [--snr-min 20] [--max-stars 300] [--saturation 0] [--fwhm 3] [--star-list FILE] [--grid 7]

Writes in --work: psfex_model.json / psfex_model.fits (the model), psfex_center.fits (PSF stamp at the field centre: the plain PSF image the
other steps take), psfex_stars.tsv (per star: position, flux, measured FWHM / ellipticity / PA, residual rms fraction, used flag), psfex_maps.tsv
(FWHM / ellipticity / PA on a grid), psfex_info.json, psfex_diag.png.  Few stars: constant model (>= 4 stars), Moffat fit (1-3), Gaussian of
the prior FWHM (none); the mode is printed and stored in the JSON.
With --catalog: `NUMBER PSFM_FWHM PSFM_E PSFM_PA PSFM_NSTAR` = the model's FWHM / ellipticity / PA at every catalog object's position.
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.abspath(os.path.join(_here, '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio, psfmodel as pm  # noqa: E402

COLUMNS = ['PSFM_FWHM', 'PSFM_E', 'PSFM_PA', 'PSFM_NSTAR']


def build(data, mask=None, kind='empirical', order=2, oversample='auto', size=0, snr_min=20.0, max_stars=300, saturation=None, fwhm=3.0, xy=None, neighbour_iter=0):
    """Server-ready entry point: returns (PSFModel, info dict, star rows, grid rows).

    neighbour_iter > 0: after the first model the other sources of the field are fitted with it (ALLSTAR-like group fit, `ogfkit.daophot.allstar`) and subtracted from the
    PSF-star stamps before the model is rebuilt (instead of masking the neighbours), `neighbour_iter` times."""
    nb_fit = None
    if neighbour_iter > 0:
        from ogfkit import daophot as dp
        bkg, rms = pm.background(data, mask)
        f = dp.find(data, bkg, rms, fwhm=fwhm, thresh=4.0, mask=mask)

        def nb_fit(mdl, _s):
            r, _ = dp.allstar(data, mdl, bkg=bkg, rms=rms, mask=mask, init_stars=f, fwhm=fwhm, thresh=4.0, n_iter=1, refind=False, fit_sky=False)
            q = [t for t in r['stars'] if t['flux'] > 0]
            return [t['x'] for t in q], [t['y'] for t in q], [t['flux'] for t in q]
    mdl, info = pm.build_psf_model(data, mask=mask, xy=xy, fwhm_prior=fwhm, size=size or None, oversample=oversample, degree=order, snr_min=snr_min,
                                   saturation=saturation, max_stars=max_stars, kind=kind, neighbour_fit=nb_fit, n_iter=max(2, neighbour_iter + 1))
    return mdl, info, pm.star_table(info), pm.grid_table(mdl, 7)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--work', default='.')
    ap.add_argument('--catalog', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--kind', default='empirical', choices=['empirical'] + list(pm.ANALYTIC_KINDS))
    ap.add_argument('--order', type=int, default=2)
    ap.add_argument('--oversample', default='auto')
    ap.add_argument('--size', type=int, default=0)
    ap.add_argument('--snr-min', type=float, default=20.0)
    ap.add_argument('--max-stars', type=int, default=300)
    ap.add_argument('--saturation', type=float, default=0.0)
    ap.add_argument('--fwhm', type=float, default=3.0)
    ap.add_argument('--star-list', default='', help='TSV with X_IMAGE Y_IMAGE (1-based) of the PSF stars; default: automatic selection')
    ap.add_argument('--neighbour-iter', type=int, default=0, help='iterations of neighbour subtraction (group fit of the other sources) before the model is rebuilt; 0 = mask neighbours')
    ap.add_argument('--grid', type=int, default=7)
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, 'psfex_' + n)
    data, hdr = imageio.load_image(a.image)
    mask = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    xy = None
    if a.star_list and os.path.isfile(a.star_list):
        _, rows = tsvio.read_catalog(a.star_list)
        xy = np.array([[tsvio.fnum(r['X_IMAGE']) - 1, tsvio.fnum(r['Y_IMAGE']) - 1] for r in rows])
    mdl, info, stars, grid = build(data, mask, a.kind, a.order, a.oversample if a.oversample == 'auto' else int(a.oversample), a.size, a.snr_min, a.max_stars,
                                   a.saturation if a.saturation > 0 else None, a.fwhm, xy, a.neighbour_iter)
    mdl.save(W('model.json'))
    mdl.save(W('model.fits'))
    imageio.save_fits(W('center.fits'), mdl.stamp(0.5 * sum(mdl.xr), 0.5 * sum(mdl.yr), 0, 0).astype(np.float32))
    tsvio.write_table(W('stars.tsv'), ['ID', 'X_IMAGE', 'Y_IMAGE', 'FLUX', 'FWHM', 'E', 'PA', 'RESID_FRAC', 'USED'], stars)
    tsvio.write_table(W('maps.tsv'), ['X', 'Y', 'FWHM', 'E', 'PA'], grid)
    with open(W('info.json'), 'w') as fh:
        json.dump(pm._jsonable(dict(info={k: v for k, v in info.items() if k != 'stars'}, summary=mdl.summary(), diag=pm.diagnostics(mdl, a.grid))), fh)
    try:
        pm.plot_maps(mdl, info, W('diag.png'))
    except Exception as e:  # plotting must never fail the step
        sys.stderr.write('psfex: plot failed: %s\n' % e)
    sm = mdl.summary()
    lines = ['PSF model: mode %s, %s stars (%d candidates, %d on the stellar locus), polynomial order %d, oversampling %d, stamp %d px' % (
        mdl.meta.get('mode'), info.get('n_used', info.get('n_stars')), info.get('n_candidates', 0), info.get('n_locus', info.get('n_stars', 0)), mdl.degree, mdl.oversample, mdl.size),
        'FWHM %.2f - %.2f px (field mean %.2f), ellipticity %.3f - %.3f' % (sm['fwhm_pix_min'], sm['fwhm_pix_max'], sm['constant_fwhm_pix'], sm['e_min'], sm['e_max']),
        'files: %s, %s, %s' % (W('model.json'), W('center.fits'), W('stars.tsv'))]
    if a.catalog:
        cols, rows = tsvio.read_catalog(a.catalog)
        x = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows])
        y = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
        ok = np.isfinite(x) & np.isfinite(y)
        f, e, pa = pm.shape_at(mdl, np.where(ok, x, 0), np.where(ok, y, 0))
        out = [(r['NUMBER'], dict(PSFM_FWHM=f[i] if ok[i] else None, PSFM_E=e[i] if ok[i] else None, PSFM_PA=pa[i] if ok[i] else None, PSFM_NSTAR=mdl.meta.get('n_stars'))) for i, r in enumerate(rows)]
        for t in lines:
            sys.stderr.write(t + '\n')
        tsvio.write_columns(COLUMNS, out)
    else:
        print('\n'.join(lines))
    return 0


if __name__ == '__main__':
    sys.exit(main())
