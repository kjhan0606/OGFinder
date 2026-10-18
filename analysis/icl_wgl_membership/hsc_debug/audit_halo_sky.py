#!/usr/bin/env python3
"""Audit paired halo/sky artifacts, independently refit samples, report by cohort."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

import validate_halo_sky as experiment
import prepare_detection_psf as prep


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def key(row):
    return tuple(row[x] for x in ('cluster_id', 'scenario', 'repeat', 'true_diffuse_flux',
                                 'true_neighbor_halo_flux', 'method'))


def independent_fit(cols, data, image, anchor):
    valid = ~data['base'] & np.isfinite(data['variance']) & (data['variance'] > 0)
    weight = 1/np.sqrt(data['variance'][valid])
    matrix = np.column_stack([c[valid]*weight for c in cols])
    y = image[valid]*weight
    if anchor is not None:
        constraint = np.zeros((6, len(cols)))
        constraint[:, -6:] = np.eye(6)/experiment.ANCHOR_SIGMA
        matrix = np.vstack([matrix, constraint]); y = np.r_[y, anchor/experiment.ANCHOR_SIGMA]
    norms = np.linalg.norm(matrix, axis=0)
    coef, _, rank, _ = np.linalg.lstsq(matrix/norms, y, rcond=None)
    assert rank == len(cols)
    coef /= norms
    _, singular, vt = np.linalg.svd(matrix/norms, full_matrices=False)
    covariance = ((vt.T/singular**2)@vt)/norms[:, None]/norms[None, :]
    return coef, covariance, float(np.sum((y-matrix@coef)**2))


def cohort_metrics(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row['cohort'], row['scenario'], row['method'])].append(row)
    records = []
    for (cohort, scenario, method), subset in groups.items():
        good = [r for r in subset if r['status'] == 'ok']
        result = dict(cohort=cohort, scenario=scenario, method=method, n_total=len(subset),
                      n_numerical_success=len(good), statuses=dict(Counter(r['status'] for r in subset)),
                      n_information_veto=sum(r.get('information_status') == 'blocked_hidden_component' for r in good))
        if good:
            median = lambda name: float(np.median([float(r[name]) for r in good]))
            result.update(median_diffuse_bias_fraction=median('diffuse_bias_fraction'),
                          median_absolute_diffuse_bias_fraction=float(np.median([abs(float(r['diffuse_bias_fraction'])) for r in good])),
                          median_conditional_sigma=median('conditional_diffuse_sigma'),
                          median_profile_half_width=float(np.median([(float(r['profile_high'])-float(r['profile_low']))/2 for r in good])),
                          profile_truth_containment=float(np.mean([r['profile_contains_truth']=='True' for r in good])),
                          conditional_truth_containment=float(np.mean([r['conditional_contains_truth']=='True' for r in good])),
                          median_supported_scales=median('supported_scale_count'),
                          n_grid_boundary=sum(r['scale_at_boundary']=='True' for r in good),
                          n_multiple_supported_scales=sum(int(r['supported_scale_count'])>1 for r in good),
                          median_conditional_diffuse_halo_correlation=median('conditional_diffuse_halo_correlation'))
            for label, selected in (
                    ('zero_halo', [r for r in good if float(r['true_neighbor_halo_flux']) == 0]),
                    ('halo_present', [r for r in good if float(r['true_neighbor_halo_flux']) > 0])):
                result[label] = dict(n=len(selected),
                    n_grid_boundary=sum(r['scale_at_boundary']=='True' for r in selected),
                    n_multiple_supported_scales=sum(int(r['supported_scale_count'])>1 for r in selected))
            paired = defaultdict(list)
            cells = defaultdict(list)
            for row in subset:
                paired[(row['cluster_id'], row['repeat'])].append(row)
                cells[(row['cluster_id'], row['true_diffuse_flux'], row['true_neighbor_halo_flux'])].append(row)
            transfers = []
            for group in paired.values():
                grid = {(float(r['true_diffuse_flux']), float(r['true_neighbor_halo_flux'])) for r in group}
                if len(group) != 6 or len(grid) != 6 or any(r['status'] != 'ok' for r in group):
                    continue
                x = np.array([[1., float(r['true_diffuse_flux']), float(r['true_neighbor_halo_flux'])] for r in group])
                y = np.array([[float(r['diffuse_estimate']), float(r['halo_estimate'])] for r in group])
                coef = np.linalg.lstsq(x, y, rcond=None)[0]
                transfers.append(coef)
            result.update(n_complete_amplitude_grids=len(transfers), n_amplitude_grids=len(paired))
            if transfers:
                middle = np.median(transfers, axis=0)
                result.update(diffuse_gain=float(middle[1, 0]), halo_to_diffuse_response=float(middle[2, 0]),
                              diffuse_to_halo_response=float(middle[1, 1]), halo_gain=float(middle[2, 1]))
            residuals = []
            for group in cells.values():
                if len(group) < 2 or any(r['status'] != 'ok' for r in group):
                    continue
                e = np.array([[float(r['diffuse_estimate'])-float(r['true_diffuse_flux']), float(r['halo_error'])] for r in group])
                residuals.extend(e-e.mean(axis=0))
            if residuals:
                result['demeaned_noise_error_correlation'] = float(np.corrcoef(np.asarray(residuals).T)[0, 1])
        records.append(result)
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', type=Path, required=True)
    p.add_argument('--output', type=Path, help='New directory for an additional audit; preserve existing reports')
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    destination = a.output or a.work
    if a.output is not None:
        destination.mkdir(parents=True, exist_ok=False)
    elif any((destination/name).exists() for name in ('cohort_report.json', 'artifact_audit.json', 'halo_sky_comparison.png')):
        p.error('Reports already exist; use --output with a new directory')
    summary = json.loads((a.work/'summary.json').read_text()); old = Path(summary['previous_work'])
    assert not summary['science_measurement'] and not summary['observational_validation_passed']
    assert hashlib.sha256((old/'summary.json').read_bytes()).hexdigest() == summary['previous_summary_sha256']
    for name, digest in summary['code_sha256'].items():
        assert hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() == digest, name
    rows = read_table(a.work/'recovery.tsv'); profiles = read_table(a.work/'scale_profiles.tsv')
    by_key = {key(r): r for r in rows}; scale_groups = defaultdict(list)
    assert len(rows) == len(by_key) == summary['n_method_rows']
    expected_keys = set()
    for case in summary['cases']:
        for scenario in summary['scenarios']:
            for repeat in range(summary['repeats']):
                for dflux in (1500., 4500.):
                    for hflux in (0., 500., 1000.):
                        for method in summary['methods']:
                            expected_keys.add((case['id'], scenario, str(repeat), str(dflux), str(hflux*case['n_neighbors']), method))
    assert set(by_key) == expected_keys
    assert len(profiles) == summary['n_scale_fits']
    for r in profiles:
        assert key(r) in by_key
        scale_groups[key(r)].append(r)
    for row in rows:
        cid = int(row['cluster_id']); repeat = int(row['repeat'])
        assert int(row['seed']) == 2026101300+cid*10+repeat
        assert int(row['sky_seed']) == 3026101300+cid*10+repeat
        group = scale_groups[key(row)]; failed = json.loads(row['failed_scales'])
        expected_scales = {8.} if row['method'] == 'fixed8_plane' else set(summary['scales_kpc'])
        actual = {float(r['scale_kpc']) for r in group}; failed_set = {r['scale'] for r in failed}
        assert len(actual) == len(group) and not actual & failed_set and actual | failed_set == expected_scales
        if row['status'] != 'ok':
            assert not group
            continue
        best = min(group, key=lambda r: float(r['chi2']))
        assert float(best['scale_kpc']) == float(row['best_scale_kpc'])
        support = [r for r in group if float(r['chi2'])-float(best['chi2']) <= 1.+1e-9]
        bounds = [(float(r['diffuse_estimate']), np.sqrt(max(0., 1.-float(r['chi2'])+float(best['chi2'])))*float(r['conditional_diffuse_sigma'])) for r in support]
        np.testing.assert_allclose([float(row['profile_low']), float(row['profile_high'])],
                                   [min(x-s for x, s in bounds), max(x+s for x, s in bounds)], atol=1e-5)
        assert len(support) == int(row['supported_scale_count'])
        assert float(row['profile_low']) <= float(row['diffuse_estimate']) <= float(row['profile_high'])
        assert (row['information_status'] == 'blocked_hidden_component') == (float(row['min_visible_fraction']) < .01)
    old_rows = read_table(old/'recovery.tsv'); baseline_diffs = []
    for row in old_rows:
        if row['method'] != 'neighbor_core_halo':
            continue
        k = key(dict(row, method='fixed8_plane'))
        if k not in by_key:
            continue
        newer = by_key[k]; difference = abs(float(newer['diffuse_estimate'])-float(row['diffuse_estimate']))
        assert difference < .001, (k, difference)
        baseline_diffs.append(difference)
    expected_baselines = len(summary['cases'])*4*summary['repeats']*6
    assert len(baseline_diffs) == expected_baselines
    audits = []
    for case in summary['cases']:
        cid = case['id']; path = old/'inputs'/(cid+'.npz')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == case['saved_input_sha256']
        with np.load(path) as saved:
            data = {k: saved[k] for k in saved.files}
        repeat = summary['repeats']-1; kpp = case['kpc_per_pixel']; n = len(data['members'])
        noise = np.random.default_rng(2026101300+int(cid)*10+repeat).normal(size=data['base'].shape)*np.sqrt(data['variance'])
        for scenario in ('matched', 'extended_curved'):
            templates = experiment.truth_templates(data, kpp, scenario)
            image = experiment.synthetic_image(templates, noise, 4500., 1000.)
            anchor = templates[5]+experiment.ANCHOR_SIGMA*np.random.default_rng(3026101300+int(cid)*10+repeat).normal(size=6)
            for method in summary['methods']:
                row = by_key[(cid, scenario, str(repeat), '4500.0', str(float(1000*n)), method)]
                if row['status'] != 'ok':
                    continue
                quadratic = not method.endswith('plane')
                measurement = anchor.copy() if method.endswith('anchor') else None
                if method == 'grid_biased_anchor':
                    measurement[[3, 5]] += experiment.ANCHOR_CURVATURE_BIAS
                cols = experiment.make_columns(data, kpp, float(row['best_scale_kpc']), quadratic)
                coef, covariance, chi2 = independent_fit(cols, data, image, measurement)
                np.testing.assert_allclose(coef[1], float(row['diffuse_estimate']), atol=1e-5, rtol=1e-8)
                np.testing.assert_allclose(coef[2+n:2+2*n].sum(), float(row['halo_estimate']), atol=1e-4, rtol=1e-8)
                np.testing.assert_allclose(np.sqrt(covariance[1, 1]), float(row['conditional_diffuse_sigma']), atol=1e-5, rtol=1e-8)
                np.testing.assert_allclose(chi2, float(row['chi2']), atol=1e-5, rtol=1e-10)
                audits.append(dict(cluster_id=cid, scenario=scenario, method=method,
                                   diffuse_absolute_difference=float(abs(coef[1]-float(row['diffuse_estimate'])))))
    metrics = cohort_metrics(rows)
    prep.write_json(destination/'cohort_report.json', dict(science_measurement=False, records=metrics,
        note='Descriptive paired-cell statistics. Only three independent seeds per case. Profile containment is not calibrated coverage; linear response slopes summarize potentially nonlinear grid selection. Stress case separate.'))
    prep.write_json(destination/'artifact_audit.json', dict(passed=True, rows=len(rows), scale_fits=len(profiles),
        paired_previous_baselines=len(baseline_diffs), max_previous_diffuse_difference=max(baseline_diffs),
        independent_svd_replays=len(audits), max_svd_diffuse_difference=max(r['diffuse_absolute_difference'] for r in audits),
        replay_records=audits, science_measurement=False, observational_validation_passed=False,
        audit_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    from matplotlib import pyplot as plt
    selected = ('extended_neighbor_halo', 'curved_sky', 'extended_curved', 'off_grid_curved', 'unmodeled_sky')
    colors = ('#777777', '#df8b20', '#4469a8', '#268268', '#b64c71')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    area = [r for r in metrics if r['cohort'] == 'area_qualified']
    lookup = {(r['scenario'], r['method']): r for r in area}
    for i, method in enumerate(summary['methods']):
        x = np.arange(len(selected)) + (i-2)*.15
        axes[0].bar(x, [100*lookup[(s, method)]['median_diffuse_bias_fraction'] for s in selected], width=.14, color=colors[i], label=method)
        axes[1].bar(x, [lookup[(s, method)]['median_profile_half_width'] for s in selected], width=.14, color=colors[i])
    axes[0].axhline(0, color='k', lw=.6); axes[0].set_ylabel('Median diffuse error (%)')
    axes[1].set_ylabel('Median profile half-width (synthetic flux units)')
    for ax in axes:
        ax.set_xticks(np.arange(len(selected)), [s.replace('_', '\n') for s in selected], fontsize=8)
    axes[0].legend(fontsize=8)
    fig.suptitle('Five area-cohort geometries, paired synthetic scenes; 6940 reported separately\nKnown shapes/positions, assumed PSF, simulated sky anchors; no observational approval', fontsize=10)
    fig.tight_layout(); fig.savefig(destination/'halo_sky_comparison.png', dpi=140); plt.close(fig)
    print(json.dumps(dict(passed=True, baseline_checks=len(baseline_diffs), independent_svd=len(audits),
                          max_baseline_difference=max(baseline_diffs), max_svd_difference=max(r['diffuse_absolute_difference'] for r in audits))))
    for r in area:
        if r['scenario'] in ('extended_curved', 'unmodeled_sky', 'extended_neighbor_halo', 'curved_sky'):
            print(r['scenario'], r['method'], 'bias%', round(100*r['median_diffuse_bias_fraction'], 2),
                  'profile half-width', round(r['median_profile_half_width'], 1),
                  'containment', round(r['profile_truth_containment'], 3),
                  'halo response', round(r['halo_to_diffuse_response'], 3))


if __name__ == '__main__':
    main()
