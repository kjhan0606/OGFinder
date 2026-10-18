#!/usr/bin/env python3
"""Paired synthetic halo-size/sky experiment. No observational science approval."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

import inject_neighbor_recovery as previous
import prepare_detection_psf as prep

SCALES = (4., 6., 8., 10., 12., 14., 16.)
SCENARIOS = previous.SCENARIOS + (
    'extended_curved', 'off_grid_curved', 'heterogeneous_curved', 'unmodeled_sky')
METHODS = ('fixed8_plane', 'grid_plane', 'grid_quadratic',
           'grid_sky_anchor', 'grid_biased_anchor')
ANCHOR_SIGMA = .01
ANCHOR_CURVATURE_BIAS = .02


def sky_basis(plane):
    """Full quadratic, centred using the whole stamp, never the surviving pixels."""
    _, x, y = plane
    return plane + [v - v.mean() for v in (x*x, x*y, y*y)]


def make_columns(data, kpp, scale, quadratic):
    core, diffuse, cores, _, plane = previous.design(
        data['base'].shape, data['center'], data['members'], kpp)
    halos = [previous.profile(data['base'].shape, xy, scale/kpp, q=.8, theta=.5*i)
             for i, xy in enumerate(data['members'])]
    return [core, diffuse] + cores + halos + (sky_basis(plane) if quadratic else plane)


def operator(columns, variance, mask, n_sky, anchor_sigma=None):
    """Scaled WLS, optionally augmented by independent sky-coefficient measurements."""
    valid = ~mask & np.isfinite(variance) & (variance > 0)
    if valid.sum() < 10*len(columns):
        return None, 'too_few_pixels'
    weight = 1/np.sqrt(variance[valid])
    matrix = np.column_stack([c[valid]*weight for c in columns])
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms == 0):
        return None, 'zero_column'
    scaled = matrix/norms
    prior = np.zeros((0, len(columns)))
    if anchor_sigma is not None:
        if not np.isfinite(anchor_sigma) or anchor_sigma <= 0:
            raise ValueError('anchor_sigma must be finite and positive')
        prior = np.zeros((n_sky, len(columns)))
        prior[:, -n_sky:] = np.eye(n_sky)/anchor_sigma
        prior /= norms
    singular = np.linalg.svd(np.vstack([scaled, prior]), compute_uv=False)
    condition = float(singular[0]/singular[-1]) if singular[-1] else float('inf')
    if condition > 1e6:
        return None, 'ill_conditioned'
    gram = scaled.T@scaled + prior.T@prior
    inv = np.linalg.inv(gram)
    # Linear map into scaled coefficient units; avoid raw-unit rank decisions.
    projection = (inv@scaled.T)*weight[None, :]
    prior_projection = inv@prior.T
    return dict(valid=valid, norms=norms, gram=gram, projection=projection,
                prior_projection=prior_projection, anchor_sigma=anchor_sigma,
                covariance=inv/norms[:, None]/norms[None, :], condition=condition,
                min_visible_fraction=min((float(c[valid].sum()) for c in columns[:-n_sky]), default=1.)), 'ok'


def solve(op, image, variance, anchor=None):
    """Joint image + external-measurement chi square, including anchor residuals."""
    pixels = image[op['valid']]
    scaled_coef = op['projection']@pixels
    norm2 = float(np.sum(pixels.astype(float)**2/variance[op['valid']]))
    if op['anchor_sigma'] is not None:
        if anchor is None or not np.all(np.isfinite(anchor)):
            raise ValueError('Finite sky measurement required')
        whitened = np.asarray(anchor)/op['anchor_sigma']
        scaled_coef += op['prior_projection']@whitened
        norm2 += float(whitened@whitened)
    chi2 = norm2 - float(scaled_coef@op['gram']@scaled_coef)
    if chi2 < -1e-5 or not np.isfinite(chi2):
        raise ValueError('Invalid residual norm')
    return dict(coef=scaled_coef/op['norms'], covariance=op['covariance'],
                chi2=max(0., chi2), condition=op['condition'],
                min_visible_fraction=op['min_visible_fraction'])


def profile_result(fits, n_neighbors):
    """Profile over the declared discrete common-scale grid; no truth-based choice.

    The hull of Delta-chi2<=1 intervals includes scale ambiguity. It is an
    approximate diagnostic set, not a calibrated 68% confidence interval.
    """
    if not fits:
        return dict(status='no_valid_scale')
    best = min(fits, key=lambda v: v['chi2'])
    supported = [(f, max(0., 1 - (f['chi2']-best['chi2']))) for f in fits
                 if f['chi2']-best['chi2'] <= 1. + 1e-9]
    low = min(f['coef'][1] - np.sqrt(delta*f['covariance'][1, 1]) for f, delta in supported)
    high = max(f['coef'][1] + np.sqrt(delta*f['covariance'][1, 1]) for f, delta in supported)
    c = best['covariance']; h = slice(2+n_neighbors, 2+2*n_neighbors)
    hsigma = float(np.sqrt(c[h, h].sum()))
    dsigma = float(np.sqrt(c[1, 1]))
    return dict(status='ok', best_scale_kpc=best['scale'], diffuse_estimate=float(best['coef'][1]),
                halo_estimate=float(best['coef'][h].sum()), conditional_diffuse_sigma=dsigma,
                conditional_halo_sigma=hsigma, profile_low=float(low), profile_high=float(high),
                conditional_diffuse_halo_correlation=float(c[1, h].sum()/dsigma/hsigma),
                supported_scale_count=len(supported), valid_scale_count=len(fits),
                supported_scale_min=min(f['scale'] for f, _ in supported),
                supported_scale_max=max(f['scale'] for f, _ in supported),
                scale_at_boundary=len(fits)>1 and best['scale'] in (SCALES[0], SCALES[-1]),
                chi2=best['chi2'], condition_number=best['condition'],
                min_visible_fraction=best['min_visible_fraction'],
                information_status=('blocked_hidden_component' if best['min_visible_fraction'] < .01
                                    else 'engineering_only'))


def truth_templates(data, kpp, scenario):
    shape = data['base'].shape; center = data['center']; members = data['members']
    wing = scenario == 'psf_wing_mismatch'
    core = previous.profile(shape, center, 5/kpp, q=.85, theta=.3, wing=wing)
    diffuse = previous.profile(shape, center, 35/kpp, wing=wing)
    cores = [previous.profile(shape, xy, 2/kpp, q=.75, theta=.5*i, wing=wing)
             for i, xy in enumerate(members)]
    length = 12 if scenario in ('extended_neighbor_halo', 'extended_curved') else (
        10.5 if scenario == 'off_grid_curved' else 8)
    lengths = [float((6, 11, 16)[i]) if scenario == 'heterogeneous_curved' else float(length)
               for i in range(len(members))]
    halos = [previous.profile(shape, xy, lengths[i]/kpp, q=.8, theta=.5*i, wing=wing)
             for i, xy in enumerate(members)]
    *_, plane = previous.design(shape, center, members, kpp)
    basis = sky_basis(plane)
    theta = np.array([.03, .02, -.01, 0., 0., 0.])
    if scenario in ('curved_sky', 'extended_curved', 'off_grid_curved',
                    'heterogeneous_curved', 'unmodeled_sky'):
        theta[[3, 5]] = .05
    sky = sum(t*c for t, c in zip(theta, basis))
    if scenario == 'unmodeled_sky':
        # Local sky feature absent from control-field coefficients and polynomial family.
        # This failure case deliberately resembles broad target light.
        x, y = plane[1:]; bump = np.exp(-((x-.15)**2+(y+.1)**2)/(2*.3**2))
        sky += .06*(bump-bump.mean())
    return core, diffuse, cores, halos, sky, theta, lengths


def synthetic_image(templates, noise, diffuse_flux, halo_flux):
    core, diffuse, cores, halos, sky, _, _ = templates
    injected = (4000*core + diffuse_flux*diffuse
                + sum((2000+1000*i)*v for i, v in enumerate(cores))
                + sum(halo_flux*v for v in halos))
    return np.asarray(injected + sky + noise, dtype='f4')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--previous', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--case-limit', type=int, default=6)
    p.add_argument('--repeats', type=int, default=3)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    if not 1 <= a.case_limit <= 6 or not 1 <= a.repeats <= 3:
        p.error('Bounded experiment: cases 1..6, repeats 1..3')
    old = json.loads((a.previous/'summary.json').read_text())
    cases = old['cases'][:a.case_limit]
    a.work.mkdir(parents=True, exist_ok=False)
    rows = []; details = []; started = time.time()
    for case in cases:
        cid = case['id']; kpp = case['kpc_per_pixel']
        path = a.previous/'inputs'/(cid+'.npz')
        if hashlib.sha256(path.read_bytes()).hexdigest() != case['saved_input_sha256']:
            raise ValueError('Input hash changed: '+cid)
        with np.load(path) as saved:
            data = {key: saved[key] for key in saved.files}
        var = data['variance']; shape = data['base'].shape; n = len(data['members'])
        banks = {}
        for family in ('plane', 'quadratic', 'anchor'):
            banks[family] = []
            for scale in SCALES:
                quadratic = family != 'plane'
                cols = make_columns(data, kpp, scale, quadratic)
                op, status = operator(cols, var, data['base'], 6 if quadratic else 3,
                                      ANCHOR_SIGMA if family == 'anchor' else None)
                banks[family].append((scale, op, status))
        for scenario in SCENARIOS:
            templates = truth_templates(data, kpp, scenario)
            for repeat in range(a.repeats):
                seed = 2026101300+int(cid)*10+repeat
                sky_seed = 3026101300+int(cid)*10+repeat
                noise = np.random.default_rng(seed).normal(size=shape)*np.sqrt(var)
                anchor = templates[5] + ANCHOR_SIGMA*np.random.default_rng(sky_seed).normal(size=6)
                biased = anchor.copy(); biased[[3, 5]] += ANCHOR_CURVATURE_BIAS
                for dflux in (1500., 4500.):
                    for hflux in (0., 500., 1000.):
                        image = synthetic_image(templates, noise, dflux, hflux)
                        for method in METHODS:
                            family = ('plane' if method.endswith('plane') else
                                      'quadratic' if method == 'grid_quadratic' else 'anchor')
                            measurement = biased if method == 'grid_biased_anchor' else anchor
                            common = dict(cluster_id=cid, cohort=case['cohort'], scenario=scenario,
                                          repeat=repeat, seed=seed, sky_seed=sky_seed, method=method,
                                          true_diffuse_flux=dflux, true_neighbor_halo_flux=hflux*n)
                            fits = []; failed = []
                            for scale, op, status in banks[family]:
                                if method == 'fixed8_plane' and scale != 8:
                                    continue
                                if op is None:
                                    failed.append(dict(scale=scale, reason=status)); continue
                                result = solve(op, image, var, measurement)
                                result['scale'] = scale; fits.append(result)
                            record = dict(common, **profile_result(fits, n), failed_scales=json.dumps(failed))
                            if record['status'] == 'ok':
                                error = record['diffuse_estimate']-dflux
                                record.update(diffuse_bias_fraction=error/dflux,
                                              halo_error=record['halo_estimate']-hflux*n,
                                              profile_contains_truth=record['profile_low'] <= dflux <= record['profile_high'],
                                              conditional_contains_truth=abs(error) <= record['conditional_diffuse_sigma'])
                            rows.append(record)
                            for fit in fits:
                                details.append(dict(common, scale_kpc=fit['scale'], chi2=fit['chi2'],
                                                    diffuse_estimate=float(fit['coef'][1]),
                                                    conditional_diffuse_sigma=float(np.sqrt(fit['covariance'][1, 1]))))
            print(cid, scenario, 'complete', flush=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    prep.table(a.work/'recovery.tsv', rows, fields)
    prep.table(a.work/'scale_profiles.tsv', details, list(details[0]))
    prep.write_json(a.work/'summary.json', dict(status='complete', cases=cases, repeats=a.repeats,
        n_scenes=len(rows)//len(METHODS), n_method_rows=len(rows), n_scale_fits=len(details),
        previous_work=str(a.previous.resolve()), scenarios=SCENARIOS, methods=METHODS,
        scales_kpc=SCALES, sky_anchor_sigma=ANCHOR_SIGMA, biased_anchor_curvature_offset=ANCHOR_CURVATURE_BIAS,
        elapsed_seconds=time.time()-started, science_measurement=False, observational_validation_passed=False,
        previous_summary_sha256=hashlib.sha256((a.previous/'summary.json').read_bytes()).hexdigest(),
        code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                     (Path(__file__), Path(previous.__file__), Path(prep.__file__))},
        software=dict(numpy=np.__version__),
        limitations=[
            'Synthetic images only; real saved variance, masks, geometry from previous injection run.',
            'Fixed known source positions, shapes, BCG/diffuse sizes, and assumed narrow Gaussian PSF.',
            'One shared discrete halo scale for all neighbours, signed linear amplitudes.',
            'Sky anchor is a simulated independent coefficient measurement, not observed control images.',
            'Anchor errors are independent of image noise; paired seeds reused across truth cells and methods.',
            'Biased anchor offsets both x^2 and y^2 by +0.02 (two declared sigma each).',
            'Unmodeled sky is a local broad bump absent from simulated external anchor.',
            'Delta chi-square 1 profile hull is approximate and restricted to the declared grid; no calibrated coverage claim.',
            'One-percent template visibility is a diagnostic veto, not a science acceptance rule.',
            'No correlated noise, source shot noise, empirical PSF, missing sources, or survey sky-transfer validation.',
            'Area cohort five cases and one mask-stress case remain separate; no representative 800-cluster claim.']))
    print(json.dumps(dict(work=str(a.work), scenes=len(rows)//5, rows=len(rows), elapsed=time.time()-started)))


if __name__ == '__main__':
    main()
