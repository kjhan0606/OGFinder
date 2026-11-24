#!/usr/bin/env python3
"""Decompose the adjusted ICL–satellite slope on the fixed matched sample."""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_s23b_icl_satellite_relation import (
    BLOCK_WIDTHS, GRID_PHASES, connected_pair_blocks, model_matrix,
    read_tsv, sha256, write_tsv,
)

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
RELATION = DATA/'s23b_icl_satellite_relation_2026-11-07_03'
RESIDUAL = DATA/'s23b_icl_control_radial_residual_2026-11-18_01'
REPLICAS = 5000


def joint_slopes(design, response, labels, seed):
    groups, inverse = np.unique(labels, return_inverse=True)
    if len(groups) < design.shape[1] + 2:
        raise ValueError('too few connected blocks')
    # Aggregate sufficient statistics by block before drawing multiplicities.
    gram = np.array([design[inverse == i].T @ design[inverse == i] for i in range(len(groups))])
    rhs = np.array([design[inverse == i].T @ response[inverse == i] for i in range(len(groups))])
    rng = np.random.default_rng(seed)
    multiplicities = rng.multinomial(len(groups), np.full(len(groups), 1/len(groups)), size=REPLICAS)
    draw_gram = np.einsum('bi,ijk->bjk', multiplicities, gram)
    draw_rhs = np.einsum('bi,ijk->bjk', multiplicities, rhs)
    good = np.linalg.cond(draw_gram) < 1e10
    if good.sum() < .9 * REPLICAS:
        raise ValueError('too few valid bootstrap fits')
    slopes = np.linalg.solve(draw_gram[good], draw_rhs[good])[:, 1, :]
    point = np.linalg.lstsq(design, response, rcond=None)[0][1]
    return point, slopes, len(groups)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    paths = {}
    for folder, expected, name in (
        (RELATION, 'fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39', 'paired_icl_satellite_profiles.tsv'),
        (RESIDUAL, '24bdf0f8e4cac0e35fc592193a5f7c901fd26d8c4059c88ff2998a759288fd5f', 'pair_contrasts.tsv'),
    ):
        summary = folder/'summary.json'
        if sha256(summary) != expected:
            raise ValueError(f'changed saved summary: {summary}')
        table = folder/name
        digest = json.loads(summary.read_text())['output_sha256'][name]
        if sha256(table) != digest:
            raise ValueError(f'changed saved input: {table}')
        paths.update({str(summary): expected, str(table): digest})
    relation_rows = read_tsv(RELATION/'paired_icl_satellite_profiles.tsv')
    residual_rows = read_tsv(RESIDUAL/'pair_contrasts.tsv')

    def key(row):
        return row['target_id'], row['control_id'], row['filter']

    pair_lookup = {key(row): row for row in residual_rows}
    relation_keys = {key(row) for row in relation_rows}
    residual_supported = {key(row) for row in residual_rows if row['delta_contrast']}
    if (len(pair_lookup) != len(residual_rows) or len(relation_keys) != len(relation_rows)
            or relation_keys != residual_supported):
        raise ValueError('pair ledgers differ')
    fits = []
    for band in ('i', 'i2'):
        supported = [r for r in relation_rows if r['filter'] == band and r['icl_relative_contrast_njy_arcsec2']]
        for floor in (0., .1):
            rows = [r for r in supported if min(float(r['inner_common_support_fraction']),
                                                float(r['outer_common_support_fraction'])) >= floor]
            design, delta, _ = model_matrix(rows)
            response = np.array([[float(pair_lookup[key(r)][f'{role}_contrast']) for role in ('target', 'control')] for r in rows])
            if not np.allclose(response[:, 0]-response[:, 1], delta, rtol=0, atol=1e-9):
                raise ValueError('saved target minus control differs from relation outcome')
            controls = {r['target_id']: {'ra_deg': r['control_ra_deg'], 'dec_deg': r['control_dec_deg']} for r in rows}
            for width in BLOCK_WIDTHS:
                for phase in GRID_PHASES:
                    labels, _ = connected_pair_blocks(rows, controls, width, phase, 84.056247280625)
                    point, slopes, n_blocks = joint_slopes(design, response, labels, 2026111803)
                    delta_slope = point[0]-point[1]
                    if not np.isclose(delta_slope, np.linalg.lstsq(design, delta, rcond=None)[0][1], atol=1e-9, rtol=0):
                        raise ValueError('point slope decomposition differs from direct fit')
                    result = dict(filter=band, support_floor=floor, block_width_deg=width,
                                  grid_phase=phase, n_pairs=len(rows), n_blocks=n_blocks, valid_bootstrap=len(slopes))
                    for role, value, draws in zip(('target','control','delta'),
                            (*point, delta_slope), (*slopes.T, slopes[:, 0]-slopes[:, 1]), strict=True):
                        low, high = np.quantile(draws, [.025,.975])
                        result.update({f'{role}_slope': float(value), f'{role}_ci_low': float(low), f'{role}_ci_high': float(high)})
                    fits.append(result)
    output.mkdir(parents=True, exist_ok=False)
    table = output/'slope_decomposition.tsv'
    write_tsv(table, fits)
    primary = [r for r in fits if r['block_width_deg'] == 5 and r['grid_phase'] == 0]
    summary = dict(science_measurement=False, units='nJy arcsec^-2 dex^-1',
                   primary=primary, n_layouts=len(fits), replicas=REPLICAS,
                   limitation='Control outcome is evaluated against the matched target member-density predictor on common pixels. The decomposition diagnoses control coupling but cannot isolate sky or mask selection.',
                   input_sha256=paths, output_sha256={table.name:sha256(table)}, runner_sha256=sha256(Path(__file__)))
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(primary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
