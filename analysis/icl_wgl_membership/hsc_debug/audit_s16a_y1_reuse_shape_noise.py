#!/usr/bin/env python3
"""Propagate shared-source shape noise through the S16A pilot geometry."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from audit_s16a_photoz_radial_support import OBJECT_QA, PLAN, RADIAL_EDGES, load_inputs, sha, write_table

GEOMETRY = PLAN.parent/'s16a_y1_pair_geometry_2026-11-18_01'
GEOMETRY_SHA = '3a9a8aff8aaabcf4adcb56facc9655f001de422ae5369a92f4b59cfc5c2b963e'


def covariance(output_rows, origin_ids, phi, amplitudes, n_outputs):
    """Coefficients already include the square root of per-component variance."""
    _, inverse = np.unique(origin_ids, return_inverse=True)
    cosine, sine = np.cos(2*phi), np.sin(2*phi)
    # Adjacent rows are tangential/cross. The two columns per origin are
    # independent shape components in the source's local sky basis.
    cols = np.tile(np.r_[2*inverse, 2*inverse+1], 2)
    rows = np.r_[output_rows, output_rows, output_rows+1, output_rows+1]
    values = np.tile(amplitudes, 4)*np.r_[-cosine, -sine, sine, -cosine]
    operator = coo_matrix((values, (rows, cols)), shape=(n_outputs, 2*(inverse.max()+1))).tocsr()
    return (operator @ operator.T).toarray()


def estimate(by_center, groups, membership):
    n_bins = len(RADIAL_EDGES)-1
    rows, objects, angles, amplitudes, pair_keys, labels = [], [], [], [], [], []
    for index, mode in enumerate(('random', 'lens_minus_random')):
        for k, labs in enumerate(membership):
            output_row = 2*(index*n_bins+k)
            labels.extend(f'{mode}:{k}:{component}' for component in ('t','x'))
            for lens_id in labs:
                for center in groups[lens_id]:
                    cid, role = center['center_id'], center['role']
                    if role == 'lens' and mode == 'random':
                        continue
                    pairs = by_center[cid, k]
                    weights = np.array([r[2] for r in pairs])
                    coefficient = (1 if role == 'lens' else (.2 if mode == 'random' else -.2))/len(labs)
                    rows.extend([output_row]*len(pairs))
                    objects.extend(r[0] for r in pairs)
                    pair_keys.extend(f'{cid}:{r[0]}' for r in pairs)
                    angles.extend(r[1] for r in pairs)
                    amplitudes.extend(coefficient*np.sqrt(weights)/weights.sum())
    args = (np.array(rows), np.array(angles), np.array(amplitudes), len(labels))
    actual = covariance(args[0], objects, *args[1:])
    independent = covariance(args[0], pair_keys, *args[1:])
    for matrix in (actual, independent):
        if np.linalg.eigvalsh(matrix).min() < -1e-12 or np.any(np.diag(matrix) <= 0):
            raise ValueError('invalid shape-noise covariance')
    return actual, independent, labels


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    geometry_summary = GEOMETRY/'summary.json'
    geometry_table = GEOMETRY/'pair_geometry.tsv'
    if sha(geometry_summary) != GEOMETRY_SHA:
        raise ValueError('geometry summary differs from verified pilot')
    existing = json.loads(geometry_summary.read_text())
    if any(existing[key] != value for key,value in provenance.items()) or sha(geometry_table) != existing['output_sha256'][geometry_table.name]:
        raise ValueError('geometry table or source lineage changed')
    by_center, weight_by_id = defaultdict(list), {}
    with geometry_table.open(newline='') as stream:
        for pair, geo in zip(pairs, csv.DictReader(stream, delimiter='\t'), strict=True):
            obj = pair['object_id']
            if (pair['center_id'], obj) != (geo['center_id'], geo['object_id']):
                raise ValueError('geometry and pilot pair identity differ')
            if obj in selected and obj in available:
                weight = float(pair['weight'])
                if not np.isfinite(weight) or weight <= 0 or weight_by_id.setdefault(obj, weight) != weight:
                    raise ValueError('same source must have the same positive catalog shape weight')
                by_center[pair['center_id'], int(geo['radial_bin'])].append((obj, float(geo['source_angle_rad']), weight))
    groups = defaultdict(list)
    for center in centers:
        groups[center['matched_lens_id']].append(center)
    if any(sorted(row['role'] for row in g) != ['lens']+['random']*5 for g in groups.values()):
        raise ValueError('fixed matched-control membership differs')
    membership = [sorted(lens_id for lens_id, g in groups.items() if all(by_center[c['center_id'], k] for c in g))
                  for k in range(len(RADIAL_EDGES)-1)]
    fixed = sorted(set.intersection(*map(set, membership)))
    if fixed != ['S23B_3012','S23B_9131'] or list(map(len, membership)) != [2,3,2,3,4,4,4,4]:
        raise ValueError('common support differs from prior availability audit')
    matrices, diagnostics = {}, []
    for support, group_ids in [('fixed', [fixed]*len(membership)), ('per_annulus', membership)]:
        actual, independent, labels = estimate(by_center, groups, group_ids)
        matrices.update({support+'_reuse':actual, support+'_independent':independent})
        for label, a, b in zip(labels, np.diag(actual), np.diag(independent), strict=True):
            diagnostics.append(dict(support=support, statistic=label, reuse_sd=np.sqrt(a),
                                    independent_pair_sd=np.sqrt(b), sd_ratio=np.sqrt(a/b)))
    output.mkdir(parents=True, exist_ok=False)
    table = output/'shape_noise_diagonal.tsv'
    matrix = output/'shape_noise_covariance.npz'
    write_table(table, diagnostics)
    np.savez_compressed(matrix, labels=np.asarray(labels), **matrices)
    report = dict(**provenance, geometry_summary_sha256=GEOMETRY_SHA, fixed_lens_ids=fixed,
                  per_annulus_lens_ids=membership, science_measurement=False, observed_shape_used=False,
                  model='Independent Gaussian sky components per object with variance 1/catalog shape weight. Center means use shape weights; equal groups and five equal controls per group.',
                  limitation='Only shared-object shape noise under a fixed geometry and source-availability subset. No background efficiency, full survey covariance or shear calibration.',
                  output_sha256={p.name:sha(p) for p in (table,matrix)}, runner_sha256=sha(Path(__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    for support in matrices:
        print(support, matrices[support].shape)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
