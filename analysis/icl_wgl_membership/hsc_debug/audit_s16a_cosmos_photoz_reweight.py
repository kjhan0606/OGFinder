#!/usr/bin/env python3
"""Diagnostic reweighting of COSMOS by observed estimated redshift."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from ogfkit.wlcluster import Cosmology
from ogfkit.wl_photoz import point_redshift_response


def run(output, bin_width=0.2):
    if output.exists():
        raise FileExistsError(output)
    if bin_width not in (0.1, 0.2, 0.4):
        raise ValueError('bin width must be 0.1, 0.2, or 0.4')
    intake_path = INTAKE/'summary.json'
    intake = json.loads(intake_path.read_text())
    table = INTAKE/'cosmos_join.tsv'
    if sha(table) != intake['output_sha256'][table.name]:
        raise ValueError('COSMOS table changed')
    with table.open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t')
                if r['eligible_for_selection_review'] == 'True']
    if len(rows) != intake['eligible_rows']:
        raise ValueError('COSMOS eligibility changed')
    photo, reference, shape, som = np.array([
        [float(r[k]) for k in ('photo_z', 'reference_z', 'shape_weight', 'som_weight')]
        for r in rows]).T
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    point_z = finite_point_redshifts(available)
    observed = defaultdict(list)
    for pair in pairs:
        if pair['object_id'] in selected and pair['object_id'] in point_z:
            observed[pair['center_id']].append((point_z[pair['object_id']], float(pair['weight'])))
    cosmology = Cosmology()
    edges = np.r_[np.linspace(0, 4, round(4/bin_width)+1), 8.000001]
    result = []
    support_rows = []
    for center in centers:
        if center['role'] != 'lens':
            continue
        zl = float(center['z_lens'])
        sample = np.asarray(observed[center['center_id']], float).reshape(-1, 2)
        sample = sample[sample[:, 0] > zl+0.2]
        keep = photo > zl+0.2
        p, ref, sh, so = photo[keep], reference[keep], shape[keep], som[keep]
        z, w = sample.T
        chi = float(cosmology.comoving(zl))
        observed_weight = w * (1-chi/cosmology.comoving(z))**2
        cal_weight = sh*so*(1-chi/cosmology.comoving(p))**2
        if not np.all((z >= edges[0]) & (z < edges[-1])) or not np.all((p >= edges[0]) & (p < edges[-1])):
            raise ValueError('Photo-z outside diagnostic bins')
        obs_bins = np.histogram(z, bins=edges, weights=observed_weight)[0]
        cal_bins = np.histogram(p, bins=edges, weights=cal_weight)[0]
        obs_bins /= obs_bins.sum()
        cal_bins /= cal_bins.sum()
        cal_indices = np.searchsorted(edges, p, side='right')-1
        for bin_index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
            weights = cal_weight[cal_indices == bin_index]
            total = float(weights.sum())
            support_rows.append(dict(lens_id=center['matched_lens_id'], z_min=low, z_max=high,
                                     observed_weight_fraction=float(obs_bins[bin_index]),
                                     cosmos_weight_fraction=float(cal_bins[bin_index]),
                                     cosmos_count=len(weights),
                                     cosmos_neff=total**2/float(weights @ weights) if total else 0.))
        unsupported = float(obs_bins[cal_bins == 0].sum())
        baseline = point_redshift_response(zl, p, ref, sh, so, cosmology)
        row = dict(lens_id=center['matched_lens_id'], lens_z=zl, obs_pairs=len(z),
                   unsupported_fraction=unsupported, baseline=baseline['delta_sigma_multiplier'],
                   reweighted=None, fractional_change=None, reweighted_neff=None, max_bin_factor=None)
        if unsupported == 0:
            factors = np.divide(obs_bins, cal_bins, out=np.zeros_like(obs_bins), where=cal_bins > 0)
            index = np.searchsorted(edges, p, side='right')-1
            weighted = point_redshift_response(zl, p, ref, sh, so*factors[index], cosmology)
            row.update(reweighted=weighted['delta_sigma_multiplier'],
                       fractional_change=weighted['delta_sigma_multiplier']/baseline['delta_sigma_multiplier']-1,
                       reweighted_neff=weighted['neff_calibration'],
                       max_bin_factor=float(factors.max()))
        result.append(row)
    if len(result) != 11:
        raise ValueError('Unexpected pilot lens count')
    output.mkdir(parents=True)
    dest = output/'lens_photoz_reweight.tsv'
    write_table(dest, result)
    support = output/'photoz_bin_support.tsv'
    write_table(support, support_rows)
    report = dict(science_measurement=False, calibration_applied=False, intake_sha256=sha(intake_path),
                  runner_sha256=sha(Path(__file__)), output_sha256={dest.name:sha(dest), support.name:sha(support)}, **provenance,
                  edges=edges.tolist(), n_lenses=len(result),
                  n_full_support=sum(r['reweighted'] is not None for r in result),
                  limitation='Post-hoc one-dimensional estimated-z matching with no factor clipping; '
                             'finite COSMOS sample, unmeasured covariates, redshift errors and field variance '
                             'are not controlled. Not a validated calibration.')
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bin-width', type=float, choices=(0.1, 0.2, 0.4), default=0.2)
    args = parser.parse_args()
    run(args.output.resolve(), args.bin_width)
