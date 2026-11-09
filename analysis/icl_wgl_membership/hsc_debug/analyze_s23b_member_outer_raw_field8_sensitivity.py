#!/usr/bin/env python3
"""Compare unweighted and raw member-field-8 weighted S23B outer profiles.

Field 8 has not been mapped to a documented CAMIRA column. This runner keeps
it explicitly labeled as a raw-field sensitivity, never as w_mem or a
membership probability. It uses only already-imported catalog positions and
already-measured member profiles; it makes no network request and reads no
science-image pixels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import analyze_s23b_member_outer_icl_relation as relation
from analyze_s23b_icl_spatial_blocks import connected_pair_blocks


DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
MEMBERS = DATA / 's23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/target_member_positions.tsv'
PROFILE_DIR = DATA / 's23b_camira_member_outer_profiles_full_2026-11-07_01'
PROFILE_SUMMARY = PROFILE_DIR / 'summary.json'
PROFILE_BINS = PROFILE_DIR / 'member_profile_bins.tsv'
RELATION_DIR = DATA / 's23b_camira_member_outer_icl_relation_2026-11-07_01'
RELATION_SUMMARY = RELATION_DIR / 'summary.json'
RELATION_INPUTS = RELATION_DIR / 'cluster_outer_relation_inputs.tsv'
DEFAULT_OUTPUT = DATA / 's23b_camira_member_outer_field8_sensitivity_2026-11-07_01'

EXPECTED_HASHES = {
    str(MEMBERS): 'd1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae',
    str(PROFILE_SUMMARY): 'e6b91fd7fe07caf89a6012053a46bba38be17fdb5298eb265b931ccc5e0c4cdd',
    str(PROFILE_BINS): '56ff7feaaf2c2a7dd138c62f919df3d631d7798c9ff44cdc9482564045b963e1',
    str(RELATION_SUMMARY): '584082c4907be8f9d37cd361a5bdb95887cc5163a4cf18b928fa86786149a884',
    str(RELATION_INPUTS): 'b44c375ed003a460e289589472b2fddbe8281a079b336fb408c9509028b28674',
}
EXPECTED_RELATION_RUNNER_SHA256 = 'ce0c56c0d234f805071880144e70a3ac5238d93c4664e9998374ceb20a02bc8f'
EXPECTED_BLOCK_RUNNER_SHA256 = 'd7a21f68988fb782e5bcda147ff9a9040ded1b8c3f6a12ac883598beb2cd0db5'
EXPECTED_PROFILE_RUNNER_SHA256 = 'cc91779a92b3e6cfa86a564d9b9072f17f519292f95fca17c7186f1cbe3f4c87'
RADIAL_BINS = ('5-10', '10-15', '15-20', '20-30', '30-40', '40-50', '50-75')
OUTER_BIN = '50-75'
SUPPORT_FLOORS = (0.05, 0.10, 0.20)
RA_GRID_ORIGIN_DEG = 84.056247280625


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def join_member_weights(member_rows: list[dict[str, str]],
                        profile_rows: list[dict[str, str]]) -> tuple[list[dict], dict]:
    weight_by_source: dict[tuple[str, str], float] = {}
    raw_weights = []
    for row in member_rows:
        key = (row['target_cluster_id'], row['source_row'])
        if key in weight_by_source:
            raise ValueError(f'duplicate member-position identity: {key}')
        weight = float(row['member_field8_raw'])
        if not math.isfinite(weight) or not 0.0 <= weight <= 1.0:
            raise ValueError(f'raw field 8 is outside the observed [0,1] contract: {key}')
        weight_by_source[key] = weight
        raw_weights.append(weight)

    joined = []
    for row in profile_rows:
        key = (row['target_id'], row['source_row'])
        if key not in weight_by_source:
            raise ValueError(f'profile row has no imported raw field-8 match: {key}')
        item = dict(row)
        item['raw_field8_weight'] = weight_by_source[key]
        joined.append(item)
    summary = {
        'n_imported_member_positions': len(member_rows),
        'n_unique_member_identities': len(weight_by_source),
        'n_profile_bin_rows': len(joined),
        'n_unmatched_profile_rows': 0,
        'raw_field8_min': float(np.min(raw_weights)),
        'raw_field8_median': float(np.median(raw_weights)),
        'raw_field8_max': float(np.max(raw_weights)),
        'raw_field8_zero_count': int(np.count_nonzero(np.asarray(raw_weights) == 0)),
    }
    return joined, summary


def aggregate_cluster_bins(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    grouped: dict[tuple[str, str, str], list[dict]] = {}
    for row in rows:
        key = (row['target_id'], row['filter'], row['radial_bin_kpc'])
        grouped.setdefault(key, []).append(row)
    clusters = []
    radial_summary = []
    for (target_id, band, radial_bin), members in sorted(grouped.items()):
        values = np.asarray([float(row['delta_sb_njy_arcsec2']) for row in members])
        weights = np.asarray([float(row['raw_field8_weight']) for row in members])
        weight_sum = float(weights.sum())
        if not np.isfinite(weight_sum) or weight_sum <= 0:
            raise ValueError(f'non-positive field-8 weight sum for {target_id}/{band}/{radial_bin}')
        weighted_mean = float(np.average(values, weights=weights))
        unweighted_mean = float(values.mean())
        effective_n = float(weight_sum * weight_sum / np.square(weights).sum())
        clusters.append({
            'target_id': target_id,
            'filter': band,
            'radial_bin_kpc': radial_bin,
            'n_profiled_candidates': int(values.size),
            'raw_field8_weight_sum': weight_sum,
            'raw_field8_weight_effective_n': effective_n,
            'unweighted_mean_delta_sb_njy_arcsec2': unweighted_mean,
            'raw_field8_weighted_mean_delta_sb_njy_arcsec2': weighted_mean,
            'science_measurement': False,
        })
    for band in ('i', 'i2'):
        for radial_bin in RADIAL_BINS:
            band_rows = [row for row in clusters
                         if row['filter'] == band and row['radial_bin_kpc'] == radial_bin]
            if not band_rows:
                continue
            raw = np.asarray([row['unweighted_mean_delta_sb_njy_arcsec2'] for row in band_rows])
            weighted = np.asarray([row['raw_field8_weighted_mean_delta_sb_njy_arcsec2']
                                   for row in band_rows])
            radial_summary.append({
                'filter': band,
                'radial_bin_kpc': radial_bin,
                'n_clusters': int(len(band_rows)),
                'unweighted_cluster_median_njy_arcsec2': float(np.median(raw)),
                'unweighted_cluster_q16_njy_arcsec2': float(np.quantile(raw, 0.16)),
                'unweighted_cluster_q84_njy_arcsec2': float(np.quantile(raw, 0.84)),
                'raw_field8_weighted_cluster_median_njy_arcsec2': float(np.median(weighted)),
                'raw_field8_weighted_cluster_q16_njy_arcsec2': float(np.quantile(weighted, 0.16)),
                'raw_field8_weighted_cluster_q84_njy_arcsec2': float(np.quantile(weighted, 0.84)),
                'science_measurement': False,
            })
    return clusters, radial_summary


def pair_model(X: np.ndarray, y_unweighted: np.ndarray, y_weighted: np.ndarray,
               labels: np.ndarray, replicas: int, seed: int) -> dict:
    groups = np.unique(labels)
    n_param = X.shape[1]
    if len(groups) < n_param + 2:
        return {'status': 'insufficient_connected_blocks', 'n_blocks': len(groups),
                'unweighted_slope': None, 'weighted_slope': None,
                'delta_slope': None, 'valid_draws': 0}
    beta_u, _, rank_u, _ = np.linalg.lstsq(X, y_unweighted, rcond=None)
    beta_w, _, rank_w, _ = np.linalg.lstsq(X, y_weighted, rcond=None)
    if rank_u != n_param or rank_w != n_param:
        raise ValueError('full-sample paired model is rank deficient')

    group_index = {label: index for index, label in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)),
                            size=replicas)
    sample_weights = draws[:, row_group].astype(float, copy=False)
    gram = np.empty((replicas, n_param, n_param), dtype=float)
    for i in range(n_param):
        for j in range(i, n_param):
            values = sample_weights @ (X[:, i] * X[:, j])
            gram[:, i, j] = values
            gram[:, j, i] = values
    condition = np.linalg.cond(gram)
    good = np.isfinite(condition) & (condition < 1e10)
    draw_u = np.full((replicas, n_param), np.nan)
    draw_w = np.full((replicas, n_param), np.nan)
    if good.any():
        rhs_u = np.column_stack([sample_weights @ (X[:, i] * y_unweighted)
                                 for i in range(n_param)])
        rhs_w = np.column_stack([sample_weights @ (X[:, i] * y_weighted)
                                 for i in range(n_param)])
        draw_u[good] = np.linalg.solve(gram[good], rhs_u[good, :, None])[..., 0]
        draw_w[good] = np.linalg.solve(gram[good], rhs_w[good, :, None])[..., 0]
    valid = np.isfinite(draw_u[:, 1]) & np.isfinite(draw_w[:, 1])
    if np.count_nonzero(valid) < replicas * 0.9:
        raise ValueError(f'only {np.count_nonzero(valid)}/{replicas} paired bootstrap fits are valid')
    delta = draw_w[valid, 1] - draw_u[valid, 1]

    def interval(values: np.ndarray) -> tuple[float, float]:
        return tuple(float(value) for value in np.percentile(values, [2.5, 97.5]))

    ci_u = interval(draw_u[valid, 1])
    ci_w = interval(draw_w[valid, 1])
    ci_delta = interval(delta)
    return {
        'status': 'complete',
        'n_blocks': int(len(groups)),
        'unweighted_slope': float(beta_u[1]),
        'unweighted_ci_low': ci_u[0],
        'unweighted_ci_high': ci_u[1],
        'raw_field8_weighted_slope': float(beta_w[1]),
        'raw_field8_weighted_ci_low': ci_w[0],
        'raw_field8_weighted_ci_high': ci_w[1],
        'weighted_minus_unweighted_slope': float(beta_w[1] - beta_u[1]),
        'paired_delta_ci_low': ci_delta[0],
        'paired_delta_ci_high': ci_delta[1],
        'paired_delta_ci_excludes_zero': bool(ci_delta[0] > 0 or ci_delta[1] < 0),
        'valid_draws': int(np.count_nonzero(valid)),
    }


def get_controls(rows: list[dict]) -> dict:
    return {row['target_id']: {'ra_deg': row['control_ra_deg'],
                               'dec_deg': row['control_dec_deg']}
            for row in rows}


def evaluate_rows(rows: list[dict], replicas: int, seed_offset: int,
                  support_floor: float | None = None) -> list[dict]:
    result = []
    for wi, width in enumerate(relation.WIDTHS):
        for pi, phase in enumerate(relation.PHASES):
            for bi, band in enumerate(('i', 'i2')):
                band_rows = [row for row in rows if row['filter'] == band]
                controls = get_controls(band_rows)
                labels, _ = connected_pair_blocks(
                    band_rows, controls, width, phase, RA_GRID_ORIGIN_DEG)
                unweighted = [dict(row) for row in band_rows]
                weighted = [dict(row) for row in band_rows]
                for row, source in zip(unweighted, band_rows):
                    row['member_outer_delta_njy_arcsec2'] = source['unweighted_outer_delta_njy_arcsec2']
                for row, source in zip(weighted, band_rows):
                    row['member_outer_delta_njy_arcsec2'] = source['raw_field8_weighted_outer_delta_njy_arcsec2']
                X_u, y_u, _ = relation.model_matrix(unweighted)
                X_w, y_w, _ = relation.model_matrix(weighted)
                if not np.allclose(X_u, X_w, rtol=0, atol=1e-12):
                    raise ValueError('weighted and unweighted covariate matrices differ')
                fit = pair_model(X_u, y_u, y_w, labels, replicas,
                                 relation.SEED + seed_offset + wi * 1000 + pi * 100 + bi)
                result.append({
                    'filter': band,
                    'minimum_icl_common_support_fraction': support_floor,
                    'block_width_deg': width,
                    'grid_phase_fraction': phase,
                    'n_target_clusters': len(band_rows),
                    **fit,
                    'n_bootstrap_requested': replicas,
                    'science_measurement': False,
                })
    return result


def evaluate_support_floors(rows_by_band: dict[str, list[dict]], replicas: int) -> list[dict]:
    result = []
    for floor_index, floor in enumerate(SUPPORT_FLOORS):
        subset = [row for band in ('i', 'i2') for row in rows_by_band[band]
                  if float(row['min_icl_common_support_fraction']) >= floor]
        evaluated = evaluate_rows(subset, replicas, 10000 + floor_index * 2000,
                                  support_floor=floor)
        result.extend(evaluated)
    return result


def plot_results(path: Path, layouts: list[dict], radial: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    colors = {'i': '#3568a8', 'i2': '#b65b30'}
    x = np.arange(7)
    labels = ['5–10', '10–15', '15–20', '20–30', '30–40', '40–50', '50–75']
    for ax, band in zip(axes[:2], ('i', 'i2')):
        rows = [row for row in radial if row['filter'] == band]
        rows.sort(key=lambda row: RADIAL_BINS.index(row['radial_bin_kpc']))
        y0 = np.asarray([row['unweighted_cluster_median_njy_arcsec2'] for row in rows])
        lo0 = np.asarray([row['unweighted_cluster_q16_njy_arcsec2'] for row in rows])
        hi0 = np.asarray([row['unweighted_cluster_q84_njy_arcsec2'] for row in rows])
        y1 = np.asarray([row['raw_field8_weighted_cluster_median_njy_arcsec2'] for row in rows])
        lo1 = np.asarray([row['raw_field8_weighted_cluster_q16_njy_arcsec2'] for row in rows])
        hi1 = np.asarray([row['raw_field8_weighted_cluster_q84_njy_arcsec2'] for row in rows])
        ax.errorbar(x - 0.08, y0, yerr=[y0 - lo0, hi0 - y0], fmt='o-',
                    color=colors[band], alpha=0.7, label='Unweighted')
        ax.errorbar(x + 0.08, y1, yerr=[y1 - lo1, hi1 - y1], fmt='s--',
                    color='#bd7c21', alpha=0.9, label='Raw field 8 weighted')
        ax.axhline(0, color='0.35', linewidth=0.8)
        ax.set_yscale('symlog', linthresh=2.0)
        ax.set_xticks(x, labels, rotation=35, ha='right')
        ax.set_title(f'HSC-{band}: member-centered profile', loc='left')
        ax.set_xlabel('Member-centered radius [proper kpc]')
        ax.grid(alpha=0.18)
    axes[0].set_ylabel('Candidate minus rotated control [nJy arcsec⁻²]')
    axes[0].legend(frameon=False, loc='upper right')

    ax = axes[2]
    primary = [row for row in layouts if row['block_width_deg'] == 5.0
               and row['grid_phase_fraction'] == 0.0]
    for band_index, band in enumerate(('i', 'i2')):
        row = next(item for item in primary if item['filter'] == band)
        ax.errorbar(band_index - 0.12, row['unweighted_slope'],
                    yerr=[[row['unweighted_slope'] - row['unweighted_ci_low']],
                          [row['unweighted_ci_high'] - row['unweighted_slope']]],
                    fmt='o', color=colors[band], label='Unweighted' if band_index == 0 else None)
        ax.errorbar(band_index + 0.12, row['raw_field8_weighted_slope'],
                    yerr=[[row['raw_field8_weighted_slope'] - row['raw_field8_weighted_ci_low']],
                          [row['raw_field8_weighted_ci_high'] - row['raw_field8_weighted_slope']]],
                    fmt='s', color='#bd7c21', label='Raw field 8 weighted' if band_index == 0 else None)
    ax.axhline(0, color='0.35', linewidth=0.8)
    ax.set_xticks([0, 1], ['HSC-i', 'HSC-i2'])
    ax.set_title('5° layout: ICL–outer-light slope', loc='left')
    ax.set_ylabel('Adjusted slope [nJy arcsec⁻² per nJy arcsec⁻²]')
    ax.legend(frameon=False, loc='upper right')
    ax.grid(alpha=0.18)
    fig.suptitle('Raw member-field-8 weighting sensitivity (field meaning unverified)', y=0.99)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replicas', type=int, default=5000)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error('replicas must be at least 100')
    output = args.output.resolve()
    if output.exists():
        parser.error(f'refusing to overwrite existing output directory: {output}')

    actual = {path: sha256(Path(path)) for path in EXPECTED_HASHES}
    if actual != EXPECTED_HASHES:
        changed = [path for path, expected in EXPECTED_HASHES.items()
                   if actual.get(path) != expected]
        raise ValueError(f'pinned input hashes differ: {changed}')
    helper_dir = Path(__file__).resolve().parent
    for name, expected in (
        ('analyze_s23b_member_outer_icl_relation.py', EXPECTED_RELATION_RUNNER_SHA256),
        ('analyze_s23b_icl_spatial_blocks.py', EXPECTED_BLOCK_RUNNER_SHA256),
        ('measure_s23b_camira_member_outer_profiles.py', EXPECTED_PROFILE_RUNNER_SHA256),
    ):
        if sha256(helper_dir / name) != expected:
            raise ValueError(f'pinned helper source changed: {name}')

    profile_summary = json.loads(PROFILE_SUMMARY.read_text(encoding='utf-8'))
    relation_summary = json.loads(RELATION_SUMMARY.read_text(encoding='utf-8'))
    if profile_summary.get('science_measurement') is not False or relation_summary.get('science_measurement') is not False:
        raise ValueError('input science gate is not held false')
    if profile_summary.get('status') != 'complete_exploratory_source_aware_profile_stack':
        raise ValueError('pinned member profile stack is not complete')
    if relation_summary.get('status') != 'complete_exploratory_spatial_block_relation':
        raise ValueError('pinned base relation is not complete')

    member_rows = read_tsv(MEMBERS)
    profile_rows_all = read_tsv(PROFILE_BINS)
    profile_rows, weight_summary = join_member_weights(member_rows, profile_rows_all)
    candidate_bins, radial_summary = aggregate_cluster_bins(profile_rows)
    relation_rows = read_tsv(RELATION_INPUTS)

    outer_lookup = {(row['target_id'], row['filter']): row for row in candidate_bins
                    if row['radial_bin_kpc'] == OUTER_BIN}
    relation_cluster_rows = []
    weighted_outer_rows = []
    for source in relation_rows:
        key = (source['target_id'], source['filter'])
        candidate = outer_lookup.get(key)
        if candidate is None:
            raise ValueError(f'outer relation target has no raw-field-8 profile: {key}')
        n_source = int(source['n_member_positions_50_75'])
        if candidate['n_profiled_candidates'] != n_source:
            raise ValueError(f'outer candidate count differs from pinned unweighted result: {key}')
        unweighted = float(candidate['unweighted_mean_delta_sb_njy_arcsec2'])
        pinned = float(source['member_outer_delta_njy_arcsec2'])
        tolerance = max(1e-8, 1e-10 * abs(pinned))
        if abs(unweighted - pinned) > tolerance:
            raise ValueError(f'unweighted profile did not reproduce pinned base result: {key}')
        item = dict(source)
        item['unweighted_outer_delta_njy_arcsec2'] = unweighted
        item['raw_field8_weighted_outer_delta_njy_arcsec2'] = float(
            candidate['raw_field8_weighted_mean_delta_sb_njy_arcsec2'])
        item['raw_field8_weight_sum_50_75'] = float(candidate['raw_field8_weight_sum'])
        item['raw_field8_effective_member_n_50_75'] = float(
            candidate['raw_field8_weight_effective_n'])
        relation_cluster_rows.append(item)
        weighted_outer_rows.append({
            'target_id': key[0], 'filter': key[1],
            'n_profiled_candidates_50_75': candidate['n_profiled_candidates'],
            'raw_field8_weight_sum_50_75': candidate['raw_field8_weight_sum'],
            'raw_field8_effective_member_n_50_75': candidate['raw_field8_weight_effective_n'],
            'unweighted_outer_delta_njy_arcsec2': unweighted,
            'raw_field8_weighted_outer_delta_njy_arcsec2': item['raw_field8_weighted_outer_delta_njy_arcsec2'],
            'science_measurement': False,
        })
    if len(relation_cluster_rows) != len(outer_lookup):
        raise ValueError('outer weighted profile and base association samples differ')
    rows_by_band = {band: [row for row in relation_cluster_rows if row['filter'] == band]
                    for band in ('i', 'i2')}
    if len(rows_by_band['i']) != 62 or len(rows_by_band['i2']) != 365:
        raise ValueError('raw-field sensitivity support differs from the primary profile sample')

    layouts = evaluate_rows(relation_cluster_rows, args.replicas, 0)
    floor_sensitivity = evaluate_support_floors(rows_by_band, args.replicas)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name + '.tmp.', dir=output.parent))
    write_tsv(tmp / 'cluster_outer_raw_field8_inputs.tsv', relation_cluster_rows,
              list(relation_cluster_rows[0]))
    write_tsv(tmp / 'cluster_outer_raw_field8_summary.tsv', weighted_outer_rows,
              list(weighted_outer_rows[0]))
    write_tsv(tmp / 'raw_field8_weighted_radial_profiles.tsv', radial_summary,
              list(radial_summary[0]))
    write_tsv(tmp / 'raw_field8_relation_layout_sensitivity.tsv', layouts,
              list(layouts[0]))
    write_tsv(tmp / 'raw_field8_support_floor_sensitivity.tsv', floor_sensitivity,
              list(floor_sensitivity[0]))
    plot_results(tmp / 'raw_field8_weight_sensitivity.png', layouts, radial_summary)

    primary = [row for row in layouts if row['block_width_deg'] == 5.0
               and row['grid_phase_fraction'] == 0.0]
    summary = {
        'schema': 'astrafek-s23b-member-outer-raw-field8-sensitivity-v1',
        'status': 'complete_exploratory_unverified_raw_field8_sensitivity',
        'science_measurement': False,
        'images_ready_for_science': False,
        'network_used': False,
        'science_image_pixels_read': False,
        'data_downloaded': False,
        'weight_definition': 'unverified eighth raw field in the headerless S23B member file; never labeled w_mem or membership probability',
        'measurement': 'within each target/filter/radial bin, calculate the raw-field-8 weighted mean of candidate-minus-rotated-control profiles; retain the original SEP-match and 60% target/control coverage cuts',
        'relation': 'repeat the existing richness/redshift/common-support/paired-coverage adjusted 50-75 kpc model with paired connected-block draws for weighted-minus-unweighted slope changes',
        'sample': {
            **weight_summary,
            'outer_profile_clusters': {band: len(rows_by_band[band]) for band in ('i', 'i2')},
            'outer_profile_candidate_measurements': int(sum(
                row['n_profiled_candidates_50_75'] for row in weighted_outer_rows)),
        },
        'primary_5deg_phase0': {
            band: next(row for row in primary if row['filter'] == band)
            for band in ('i', 'i2')
        },
        'support_floors': [float(value) for value in SUPPORT_FLOORS],
        'limitations': [
            'The source meaning of raw field 8 has not been verified against a file schema; this is a sensitivity analysis only.',
            'CAMIRA candidates are photometric red-sequence candidates; member purity and completeness are not calibrated.',
            'The SEP-match and annular-coverage cuts condition the profile sample; weighting does not correct missing or undetected members.',
            'Rotated controls can contain diffuse cluster light and do not remove PSF-wing, residual-sky, extinction, source-Poisson, or coadd-covariance effects.',
            'The relation is exploratory and is not a physical detection or causal estimate.',
        ],
        'input_sha256': actual,
        'helper_sha256': {
            'relation_runner': EXPECTED_RELATION_RUNNER_SHA256,
            'spatial_block_runner': EXPECTED_BLOCK_RUNNER_SHA256,
            'profile_runner': EXPECTED_PROFILE_RUNNER_SHA256,
        },
        'analysis_script_sha256': sha256(Path(__file__).resolve()),
        'bootstrap_replicas_per_layout': args.replicas,
        'science_measurement_gate': False,
    }
    (tmp / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n',
                                      encoding='utf-8')
    for path in tmp.iterdir():
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    os.replace(tmp, output)
    print(json.dumps({
        'status': summary['status'], 'output': str(output),
        'primary_5deg_phase0': summary['primary_5deg_phase0'],
        'summary_sha256': sha256(output / 'summary.json'),
    }, indent=2))


if __name__ == '__main__':
    main()
