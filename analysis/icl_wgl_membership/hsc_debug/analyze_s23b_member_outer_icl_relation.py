#!/usr/bin/env python3
"""Fit spatial-block uncertainty for source-aware S23B member outer profiles."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from analyze_s23b_icl_spatial_blocks import connected_pair_blocks

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
TARGET_MANIFEST = DATA / 's23b_icl_target_image_plan_2026-11-07_01/target_manifest.json'
CONTROL_MANIFEST = DATA / 's23b_icl_control_image_plan_2026-11-07_02/control_manifest.json'
PAIR_PROFILES = DATA / 's23b_icl_target_control_stack_2026-11-07_01/paired_radial_profiles.tsv'
RELATIVE_TABLE = DATA / 's23b_icl_relative_radial_contrasts_2026-11-07_02/paired_relative_contrasts.tsv'
STACK_DIR = DATA / 's23b_camira_member_outer_profiles_full_2026-11-07_01'
STACK_SUMMARY = STACK_DIR / 'summary.json'
SELECTED = STACK_DIR / 'selected_target_frames.tsv'
MEMBER_BINS = STACK_DIR / 'member_profile_bins.tsv'
MEMBER_DIAGNOSTICS = STACK_DIR / 'member_profile_diagnostics.tsv'
CLUSTER_RELATION = STACK_DIR / 'cluster_profile_relation.tsv'
DEFAULT_OUTPUT = DATA / 's23b_camira_member_outer_icl_relation_2026-11-07_01'
BLOCK_HELPER = Path(__file__).resolve().with_name('analyze_s23b_icl_spatial_blocks.py')

EXPECTED_HASHES = {
    str(TARGET_MANIFEST): 'e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44',
    str(CONTROL_MANIFEST): 'ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df',
    str(PAIR_PROFILES): '285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43',
    str(RELATIVE_TABLE): 'c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800',
    str(STACK_SUMMARY): 'e6b91fd7fe07caf89a6012053a46bba38be17fdb5298eb265b931ccc5e0c4cdd',
    str(SELECTED): '3fe2f2240a750be2ea597d7ceeeefef3593da2cd9cf9cb97f0cb193b1723457e',
    str(MEMBER_BINS): '56ff7feaaf2c2a7dd138c62f919df3d631d7798c9ff44cdc9482564045b963e1',
    str(MEMBER_DIAGNOSTICS): 'ea4982c568b124500dc83fc54fa5c0cefd98e5622d37169779ae0fe966694883',
    str(CLUSTER_RELATION): 'b00515ced0a3599cb3f5845e01c577ac1856411b24aead67aa60b3d87757b661',
}
EXPECTED_BLOCK_HELPER_SHA256 = 'd7a21f68988fb782e5bcda147ff9a9040ded1b8c3f6a12ac883598beb2cd0db5'
WIDTHS = (2.5, 5.0, 10.0)
PHASES = (0.0, 0.5)
BAND_ORDER = ('i', 'i2')
BOOTSTRAP_REPLICAS = 5000
SEED = 20261107
OUTER_BIN = '50-75'
CONTRAST = '50-300_minus_300-800'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8'))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def annular_support(profile_json: str, low: int, high: int) -> float:
    rows = json.loads(profile_json)
    selected = [row for row in rows
                if float(row['r_min_kpc']) >= low and float(row['r_max_kpc']) <= high]
    if not selected:
        raise ValueError(f'missing common-support bins for {low}-{high} kpc')
    valid = sum(int(row['n_valid_pixels']) for row in selected)
    geometric = sum(int(row['geometric_pixels']) for row in selected)
    if geometric <= 0:
        raise ValueError('non-positive geometric area in ICL profile')
    return valid / geometric


def assemble_rows(target_manifest: dict, control_manifest: dict,
                  pair_profiles: list[dict], relative_rows: list[dict],
                  selected_frames: list[dict], member_bins: list[dict],
                  member_diagnostics: list[dict], cluster_rows: list[dict]) -> list[dict]:
    targets = {str(row['id']): row for row in target_manifest['clusters']}
    controls = {str(row['target_cluster_id']): row for row in control_manifest['clusters']}
    frame_lookup = {(row['target_id'], row['filter']): row for row in selected_frames}
    contrast_lookup = {}
    for row in relative_rows:
        if row['contrast'] != CONTRAST:
            continue
        key = (row['target_id'], row['filter'])
        if key in contrast_lookup:
            raise ValueError(f'duplicate relative ICL row {key}')
        try:
            value = float(row['contrast_njy_arcsec2'])
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            contrast_lookup[key] = value

    support_lookup = {}
    for pair in pair_profiles:
        key = (str(pair['target_id']), str(pair['filter']))
        support_lookup[key] = min(
            annular_support(pair['profile_json'], 50, 300),
            annular_support(pair['profile_json'], 300, 800),
        )

    coverage_by_target = {}
    for row in member_bins:
        if row['radial_bin_kpc'] != OUTER_BIN:
            continue
        key = (row['target_id'], row['filter'])
        coverage_by_target.setdefault(key, {'target': [], 'control': []})
        coverage_by_target[key]['target'].append(float(row['target_coverage']))
        coverage_by_target[key]['control'].append(float(row['control_coverage_mean']))

    diag_by_target = {}
    for row in member_diagnostics:
        key = (row['target_id'], row['filter'])
        diag_by_target.setdefault(key, {'candidate_count': 0, 'sep_matched': 0, 'any_profile': 0})
        diag = diag_by_target[key]
        diag['candidate_count'] += 1
        if row['status'] != 'no_unique_sep_match_within_1arcsec':
            diag['sep_matched'] += 1
        if int(row['profile_bins_retained']) > 0:
            diag['any_profile'] += 1

    outer_clusters = {}
    for row in cluster_rows:
        if row['radial_bin_kpc'] != OUTER_BIN:
            continue
        key = (row['target_id'], row['filter'])
        if key in outer_clusters:
            raise ValueError(f'duplicate cluster outer profile {key}')
        outer_clusters[key] = row

    output = []
    for key, outer in sorted(outer_clusters.items()):
        target_id, band = key
        if key not in frame_lookup or key not in contrast_lookup or key not in support_lookup:
            raise ValueError(f'outer profile lacks pinned ICL/sample row {key}')
        if target_id not in targets or target_id not in controls:
            raise ValueError(f'missing target/control coordinates {target_id}')
        target = targets[target_id]
        control = controls[target_id]
        if control['target_cluster_id'] != target_id:
            raise ValueError(f'matched control ID mismatch for {target_id}')
        cov = coverage_by_target.get(key)
        diag = diag_by_target.get(key)
        if not cov or not diag:
            raise ValueError(f'member profile support/selection metadata missing for {key}')
        output.append({
            'target_id': target_id,
            'filter': band,
            'ra_deg': float(target['ra_deg']),
            'dec_deg': float(target['dec_deg']),
            'control_ra_deg': float(control['ra_deg']),
            'control_dec_deg': float(control['dec_deg']),
            'z_cl': float(target['z_cl']),
            'richness': float(target['richness']),
            'icl_contrast_njy_arcsec2': contrast_lookup[key],
            'member_outer_delta_njy_arcsec2': float(outer['mean_delta_sb_njy_arcsec2']),
            'n_member_positions_50_75': int(outer['n_member_positions']),
            'min_icl_common_support_fraction': support_lookup[key],
            'median_member_target_coverage_50_75': float(np.median(cov['target'])),
            'median_member_control_coverage_50_75': float(np.median(cov['control'])),
            'n_member_candidates_50_500': int(diag['candidate_count']),
            'n_sep_matched_50_500': int(diag['sep_matched']),
            'n_with_any_profile_bin_50_500': int(diag['any_profile']),
            'sep_match_fraction_50_500': (diag['sep_matched'] / diag['candidate_count']
                                          if diag['candidate_count'] else 0.0),
            'science_measurement': False,
        })
    if len(output) != len(outer_clusters):
        raise ValueError('cluster-level output row count mismatch')
    return output


def model_matrix(rows: list[dict]) -> tuple[np.ndarray, np.ndarray, dict]:
    y = np.asarray([row['member_outer_delta_njy_arcsec2'] for row in rows], dtype=float)
    x = np.asarray([row['icl_contrast_njy_arcsec2'] for row in rows], dtype=float)
    richness = np.log10(np.asarray([row['richness'] for row in rows], dtype=float))
    redshift = np.asarray([row['z_cl'] for row in rows], dtype=float)
    support = np.asarray([row['min_icl_common_support_fraction'] for row in rows], dtype=float)
    target_cov = np.asarray([row['median_member_target_coverage_50_75'] for row in rows], dtype=float)
    control_cov = np.asarray([row['median_member_control_coverage_50_75'] for row in rows], dtype=float)
    cols = [np.ones(len(rows)), x - x.mean()]
    names = ['intercept', 'icl_relative_contrast']
    centers = {'icl_contrast_center': float(x.mean())}
    nuisance = [
        ('log10_richness', richness),
        ('redshift_per_0p1', (redshift - redshift.mean()) / 0.1),
        ('minimum_icl_common_support_fraction', support),
        ('median_member_target_coverage_50_75', target_cov),
        ('median_member_control_coverage_50_75', control_cov),
    ]
    centers.update({
        'log10_richness_center': float(richness.mean()),
        'redshift_center': float(redshift.mean()),
        'minimum_icl_common_support_center': float(support.mean()),
        'median_member_target_coverage_center': float(target_cov.mean()),
        'median_member_control_coverage_center': float(control_cov.mean()),
    })
    for name, values in nuisance:
        centered = values - values.mean()
        if float(np.ptp(centered)) > 1e-10:
            cols.append(centered)
            names.append(name)
    X = np.column_stack(cols)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError(f'rank-deficient adjusted model: {names}')
    centers['model_columns'] = names
    return X, y, centers


def fit_block_slope(X: np.ndarray, y: np.ndarray, labels: np.ndarray,
                    replicas: int, seed: int) -> dict:
    groups = np.unique(labels)
    if len(groups) < X.shape[1] + 2:
        return {'status': 'insufficient_connected_blocks', 'n_blocks': len(groups),
                'slope': None, 'ci_low': None, 'ci_high': None, 'valid_draws': 0}
    full, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1]:
        raise ValueError('full-sample model is rank deficient')
    index_by_label = {label: index for index, label in enumerate(groups)}
    row_group = np.asarray([index_by_label[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)), size=replicas)
    weights = draws[:, row_group].astype(float, copy=False)
    n_param = X.shape[1]
    gram = np.empty((replicas, n_param, n_param), dtype=float)
    for i in range(n_param):
        for j in range(i, n_param):
            values = weights @ (X[:, i] * X[:, j])
            gram[:, i, j] = values
            gram[:, j, i] = values
    rhs = np.column_stack([weights @ (X[:, i] * y) for i in range(n_param)])
    cond = np.linalg.cond(gram)
    good = np.isfinite(cond) & (cond < 1e10)
    beta = np.full((replicas, n_param), np.nan)
    if good.any():
        beta[good] = np.linalg.solve(gram[good], rhs[good, :, None])[..., 0]
    slopes = beta[:, 1]
    slopes = slopes[np.isfinite(slopes)]
    if slopes.size < replicas * 0.9:
        raise ValueError(f'only {slopes.size}/{replicas} block-bootstrap fits are valid')
    low, high = np.percentile(slopes, [2.5, 97.5])
    return {
        'status': 'complete', 'n_blocks': len(groups),
        'slope': float(full[1]), 'ci_low': float(low), 'ci_high': float(high),
        'valid_draws': int(slopes.size), 'ci_excludes_zero': bool(low > 0 or high < 0),
    }


def run_layouts(rows_by_band: dict[str, list[dict]], replicas: int) -> tuple[list[dict], dict]:
    output, primary = [], {}
    origin = 84.056247280625
    for wi, width in enumerate(WIDTHS):
        for pi, phase in enumerate(PHASES):
            for bi, band in enumerate(BAND_ORDER):
                rows = rows_by_band[band]
                controls = {row['target_id']: {'ra_deg': row['control_ra_deg'],
                                                'dec_deg': row['control_dec_deg']}
                            for row in rows}
                labels, block_rows = connected_pair_blocks(rows, controls, width, phase, origin)
                X, y, metadata = model_matrix(rows)
                fit = fit_block_slope(X, y, labels, replicas,
                                      SEED + wi * 1000 + pi * 100 + bi)
                row = {
                    'filter': band, 'block_width_deg': width,
                    'grid_phase_fraction': phase, 'n_target_clusters': len(rows),
                    'n_connected_blocks': fit['n_blocks'],
                    'slope_outer_sb_per_icl_contrast': fit['slope'],
                    'block_bootstrap_ci95_low': fit['ci_low'],
                    'block_bootstrap_ci95_high': fit['ci_high'],
                    'ci_excludes_zero': fit.get('ci_excludes_zero'),
                    'n_bootstrap_requested': replicas,
                    'n_bootstrap_valid': fit['valid_draws'],
                    'status': fit['status'],
                    'model_columns': ','.join(metadata['model_columns']),
                    'science_measurement': False,
                }
                output.append(row)
                if width == 5.0 and phase == 0.0:
                    primary[band] = (rows, X, y, metadata, labels, fit)
    return output, primary


def support_sensitivity(rows_by_band: dict[str, list[dict]], replicas: int) -> list[dict]:
    output = []
    origin = 84.056247280625
    for floor_index, floor in enumerate((0.05, 0.10, 0.20)):
        for wi, width in enumerate(WIDTHS):
            for pi, phase in enumerate(PHASES):
                for bi, band in enumerate(BAND_ORDER):
                    rows = [row for row in rows_by_band[band]
                            if row['min_icl_common_support_fraction'] >= floor]
                    record = {
                        'filter': band, 'minimum_icl_common_support_fraction': floor,
                        'block_width_deg': width, 'grid_phase_fraction': phase,
                        'n_target_clusters': len(rows), 'n_connected_blocks': None,
                        'slope_outer_sb_per_icl_contrast': None,
                        'block_bootstrap_ci95_low': None, 'block_bootstrap_ci95_high': None,
                        'n_bootstrap_requested': replicas, 'n_bootstrap_valid': 0,
                        'status': 'insufficient_sample', 'reason': '',
                        'science_measurement': False,
                    }
                    try:
                        controls = {row['target_id']: {'ra_deg': row['control_ra_deg'],
                                                       'dec_deg': row['control_dec_deg']}
                                    for row in rows}
                        labels, _ = connected_pair_blocks(rows, controls, width, phase, origin)
                        X, y, _ = model_matrix(rows)
                        fit = fit_block_slope(X, y, labels, replicas,
                                              SEED + 10000 + floor_index * 2000 + wi * 1000 + pi * 100 + bi)
                        record.update({
                            'n_connected_blocks': fit['n_blocks'],
                            'slope_outer_sb_per_icl_contrast': fit['slope'],
                            'block_bootstrap_ci95_low': fit['ci_low'],
                            'block_bootstrap_ci95_high': fit['ci_high'],
                            'n_bootstrap_valid': fit['valid_draws'],
                            'status': fit['status'],
                        })
                    except (ValueError, np.linalg.LinAlgError) as exc:
                        record['reason'] = str(exc)
                    output.append(record)
    return output


def leave_one_block_out(primary: dict) -> list[dict]:
    output = []
    for band in BAND_ORDER:
        rows, X, y, _, labels, fit = primary[band]
        for block in np.unique(labels):
            keep = labels != block
            beta, _, rank, _ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            if rank != X.shape[1]:
                raise ValueError(f'leave-one-block fit rank deficient: {band}/{block}')
            output.append({
                'filter': band, 'deleted_block_id': str(block),
                'n_target_clusters_retained': int(keep.sum()),
                'full_sample_slope': fit['slope'],
                'deleted_block_slope': float(beta[1]),
                'slope_shift': float(beta[1] - fit['slope']),
                'sign_changed': bool(np.sign(beta[1]) != np.sign(fit['slope'])),
            })
    return output


def render_plot(path: Path, primary: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    colors = {'i': '#3568a8', 'i2': '#b65b30'}
    for ax, band in zip(axes, BAND_ORDER):
        rows, X, y, _, _, fit = primary[band]
        x = X[:, 1]
        nuisance = X[:, [0] + list(range(2, X.shape[1]))]
        xr = x - nuisance @ np.linalg.lstsq(nuisance, x, rcond=None)[0]
        yr = y - nuisance @ np.linalg.lstsq(nuisance, y, rcond=None)[0]
        ax.scatter(xr, yr, s=20, alpha=.45, color=colors[band], linewidths=0)
        line_x = np.linspace(float(xr.min()), float(xr.max()), 100)
        ax.plot(line_x, fit['slope'] * line_x, color='#bb4e2a', linewidth=1.8)
        ax.axhline(0, color='0.35', linewidth=.8)
        ax.axvline(0, color='0.35', linewidth=.8)
        ax.text(.04, .96,
                f'HSC-{band}  N={len(rows)}\n'
                f'adjusted slope={fit["slope"]:.3f} '
                f'[{fit["ci_low"]:.3f}, {fit["ci_high"]:.3f}]',
                transform=ax.transAxes, ha='left', va='top', fontsize=9)
        ax.set_xlabel('Residual relative ICL contrast [nJy arcsec⁻²]')
        ax.grid(alpha=.15)
    axes[0].set_ylabel('Residual member-profile excess at 50–75 kpc [nJy arcsec⁻²]')
    fig.suptitle('S23B ICL–member outer-light relation (5° connected-block model)')
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replicas', type=int, default=BOOTSTRAP_REPLICAS)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error('replicas must be at least 100')
    output = args.output.resolve()
    if output.exists():
        parser.error(f'refusing to overwrite existing output directory: {output}')
    actual = {path: sha256(Path(path)) for path in EXPECTED_HASHES}
    if actual != EXPECTED_HASHES:
        changed = [path for path, expected in EXPECTED_HASHES.items() if actual.get(path) != expected]
        raise ValueError(f'pinned input hashes differ: {changed}')
    if sha256(BLOCK_HELPER) != EXPECTED_BLOCK_HELPER_SHA256:
        raise ValueError('connected sky-block helper checksum differs')

    target = read_json(TARGET_MANIFEST)
    control = read_json(CONTROL_MANIFEST)
    stack_summary = read_json(STACK_SUMMARY)
    if (stack_summary.get('science_measurement') is not False
            or stack_summary.get('status') != 'complete_exploratory_source_aware_profile_stack'):
        raise ValueError('member-profile stack status/gate differs from the pinned exploratory run')
    if len(target.get('clusters', [])) != 1043 or len(control.get('clusters', [])) != 1041:
        raise ValueError('target/control sample counts changed')

    rows = assemble_rows(
        target, control, read_tsv(PAIR_PROFILES), read_tsv(RELATIVE_TABLE),
        read_tsv(SELECTED), read_tsv(MEMBER_BINS), read_tsv(MEMBER_DIAGNOSTICS),
        read_tsv(CLUSTER_RELATION),
    )
    rows_by_band = {band: [row for row in rows if row['filter'] == band]
                    for band in BAND_ORDER}
    if len(rows_by_band['i']) != 62 or len(rows_by_band['i2']) != 365:
        raise ValueError('outer member-profile cluster support differs from the stack summary')

    layouts, primary = run_layouts(rows_by_band, args.replicas)
    sensitivity = support_sensitivity(rows_by_band, args.replicas)
    loo = leave_one_block_out(primary)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name + '.tmp.', dir=output.parent))
    fields = list(rows[0]) + ['science_measurement'] if 'science_measurement' not in rows[0] else list(rows[0])
    write_tsv(tmp / 'cluster_outer_relation_inputs.tsv', rows, fields)
    write_tsv(tmp / 'connected_block_slope_sensitivity.tsv', layouts, list(layouts[0]))
    write_tsv(tmp / 'icl_support_floor_sensitivity.tsv', sensitivity, list(sensitivity[0]))
    write_tsv(tmp / 'leave_one_block_out.tsv', loo, list(loo[0]))
    render_plot(tmp / 'adjusted_outer_relation.png', primary)

    summary = {
        'schema': 'astrafek-s23b-camira-member-outer-icl-relation-v1',
        'status': 'complete_exploratory_spatial_block_relation',
        'science_measurement': False,
        'images_ready_for_science': False,
        'network_used': False,
        'science_image_pixels_read': False,
        'sample': {band: {'target_clusters': len(rows_by_band[band]),
                          'member_profile_bin': OUTER_BIN}
                   for band in BAND_ORDER},
        'model': {
            'outcome': 'cluster mean candidate-minus-rotated-control surface brightness at 50-75 proper kpc',
            'predictor': 'paired target-control relative ICL contrast at 50-300 minus 300-800 proper kpc',
            'adjustments': 'log richness, redshift, minimum inner/outer ICL common support, median member target/control profile coverage',
            'uncertainty': '5,000 connected target/control sky-block bootstrap draws per layout; six 2.5/5/10-degree layouts and two phases',
            'primary_layout': '5 degree, phase 0',
        },
        'primary_results': {
            band: {
                'slope': primary[band][5]['slope'],
                'ci95': [primary[band][5]['ci_low'], primary[band][5]['ci_high']],
                'n_clusters': len(rows_by_band[band]),
                'n_blocks': primary[band][5]['n_blocks'],
                'model_columns': primary[band][3]['model_columns'],
                'leave_one_block_sign_changes': sum(row['sign_changed'] for row in loo if row['filter'] == band),
            } for band in BAND_ORDER
        },
        'layout_ci_excludes_zero': {
            band: sum(bool(row['ci_excludes_zero']) for row in layouts
                      if row['filter'] == band and row['status'] == 'complete')
            for band in BAND_ORDER
        },
        'limitations': [
            'The profile sample includes SEP-matched CAMIRA photometric candidates only; unmatched, ambiguous, or low-support members are excluded.',
            'The cluster-level mean is conditioned on passing a 60% annular coverage cut and at least three rotated controls.',
            'Rotated controls are within the same cluster image and can contain diffuse cluster light.',
            'No CAMIRA membership weight is applied because the local member-file columns remain unverified.',
            'PSF wings, residual sky transfer, member purity/completeness, extinction, source-Poisson noise, and coadd pixel covariance remain uncalibrated.',
            'The fitted relation is exploratory and is not a physical detection or causal estimate.',
        ],
        'input_sha256': actual,
        'block_helper_sha256': EXPECTED_BLOCK_HELPER_SHA256,
        'analysis_script_sha256': sha256(Path(__file__).resolve()),
        'bootstrap_replicas_per_layout': args.replicas,
    }
    (tmp / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for path in tmp.iterdir():
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    os.replace(tmp, output)
    print(json.dumps({'status': summary['status'], 'output': str(output),
                      'primary_results': summary['primary_results'],
                      'layout_ci_excludes_zero': summary['layout_ci_excludes_zero'],
                      'summary_sha256': sha256(output / 'summary.json')}, indent=2))


if __name__ == '__main__':
    main()
