#!/usr/bin/env python3
"""Compare weighted point-redshift distributions of S16A and COSMOS."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from ogfkit.wlcluster import Cosmology


def weighted_cdf_at(z, weight, grid):
    order = np.argsort(z)
    z, weight = z[order], weight[order]
    if weight.sum() <= 0:
        raise ValueError('No positive weight')
    cumulative = np.concatenate(([0.], np.cumsum(weight) / weight.sum()))
    return cumulative[np.searchsorted(z, grid, side='right')]


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_path = INTAKE / 'summary.json'
    intake = json.loads(summary_path.read_text())
    table = INTAKE / 'cosmos_join.tsv'
    if sha(table) != intake['output_sha256'][table.name]:
        raise ValueError('Changed COSMOS table')
    with table.open() as stream:
        data = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    if len(data) != intake['eligible_rows']:
        raise ValueError('Changed COSMOS eligibility')
    cz, cw = np.asarray([[float(r['photo_z']), float(r['shape_weight']) * float(r['som_weight'])] for r in data]).T
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    redshifts = finite_point_redshifts(available)
    sources = defaultdict(list)
    for pair in pairs:
        if pair['object_id'] in selected and pair['object_id'] in redshifts:
            sources[pair['center_id']].append((redshifts[pair['object_id']], float(pair['weight'])))
    cosmo = Cosmology()
    output_rows = []
    for center in centers:
        array = np.asarray(sources[center['center_id']], dtype=float).reshape(-1, 2)
        zl = float(center['z_lens'])
        obs = array[array[:, 0] > zl + 0.2]
        cal = cz > zl + 0.2
        if not len(obs) or not cal.any():
            continue
        z, w = obs.T
        chi = float(cosmo.comoving(zl))
        obs_w = w * (1 - chi / cosmo.comoving(z))**2
        cal_w = cw[cal] * (1 - chi / cosmo.comoving(cz[cal]))**2
        if obs_w.sum() <= 0 or cal_w.sum() <= 0:
            continue
        grid = np.unique(np.concatenate((z, cz[cal])))
        cdf1 = weighted_cdf_at(z, obs_w, grid)
        cdf2 = weighted_cdf_at(cz[cal], cal_w, grid)
        diff = cdf1-cdf2
        index = np.argmax(np.abs(diff))
        output_rows.append(dict(center_id=center['center_id'], lens_id=center['matched_lens_id'],
                                role=center['role'], lens_z=zl, n_pairs=len(z),
                                weighted_neff=float(obs_w.sum()**2 / (obs_w @ obs_w)),
                                weighted_cdf_max_abs_difference=float(abs(diff[index])),
                                max_difference_at_z=float(grid[index])))
    output.mkdir(parents=True)
    result = output/'center_photoz_distribution.tsv'
    write_table(result, output_rows)
    report = dict(science_measurement=False, calibration_applied=False, cosmos_intake_sha256=sha(summary_path),
                  runner_sha256=sha(Path(__file__)), output_sha256={result.name:sha(result)}, **provenance,
                  n_centers=len(output_rows), gap=0.2,
                  limitation='Weighted empirical CDF distances are descriptive, not a statistical test. '
                             'Reused objects, missing photo-z, selection mismatch, and cosmic variance are not corrected.')
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
