#!/usr/bin/env python3
"""Check the fixed-52 joint relations after omitting lenses by footprint field.

The retained lenses keep their saved ICL, member, and WL profiles, including
their original matched control contributions. This is a partial-correlation
influence diagnostic with no source endpoint deletion or survey covariance.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np

from adjust_s23b_icl_satellite_wl_fixed52_max_t import (
    EXPECTED_SHA256, MIN_LENSES, load_and_check, make_models, pair_block_labels,
    partial_r, sha256, write_tsv,
)
from analyze_s23b_wl_field_component_influence import FIELD, FIELD_SHA, TARGET, TARGET_SHA
from hsc_pdr3_tract_patch_geometry import load_index


def assign_lens_fields(ids):
    for path, digest in ((FIELD / 'summary.json', FIELD_SHA), (TARGET, TARGET_SHA)):
        if sha256(path) != digest:
            raise ValueError(f'Pinned footprint input differs: {path}')
    prior = json.loads((FIELD / 'summary.json').read_text())
    if prior['status'] != 'complete_exact_field_influence_only' or prior['science_measurement'] is not False:
        raise ValueError('Unexpected footprint audit status')
    geometry = Path(prior['inputs']['official_geometry_directory'])
    geometry_digests = {name: info['sha256'] for name, info in prior['official_geometry_summary']['files'].items()}
    for name, digest in geometry_digests.items():
        if sha256(geometry / name) != digest:
            raise ValueError(f'Pinned footprint file differs: {name}')
    centers = {str(row['id']): row for row in json.loads(TARGET.read_text())['clusters']}
    index, _ = load_index(geometry, level='tract')
    assigned = index.assign([float(centers[i]['ra_deg']) for i in ids],
                            [float(centers[i]['dec_deg']) for i in ids], level='tract')
    if any(f is None for f in assigned.field) or np.any(assigned.cross_field_ambiguous) or np.any(assigned.tract < 0):
        raise ValueError('Unassigned or ambiguous official field')
    fields = dict(zip(ids, map(str, assigned.field), strict=True))
    expected = {row['region_label']: row['lens_centers_in_region'] for row in prior['region_deletions']}
    if Counter(fields.values()) != expected:
        raise ValueError('Lens field totals differ from the exact endpoint-deletion audit')
    return fields, geometry_digests


def omit_field(model, lens_fields, field):
    cluster_ids = model['cluster_ids']
    keep = np.array([lens_fields[i] != field for i in cluster_ids])
    retained = [i for i, use in zip(cluster_ids, keep, strict=True) if use]
    design = model['design'][keep]
    n = len(retained)
    rank = int(np.linalg.matrix_rank(design)) if n else 0
    estimate = None
    if n == len(cluster_ids):
        status = 'field_absent'
        estimate = model['partial_r']
    elif n < MIN_LENSES:
        status = 'below_original_lens_floor'
    elif rank != design.shape[1]:
        status = 'rank_deficient'
    else:
        outcome = (model['reduced_fit'] + model['restricted_resid'])[keep]
        value = partial_r(design[:, -1], outcome, design[:, :-1])
        status = 'complete' if math.isfinite(value) else 'undefined_correlation'
        estimate = value if math.isfinite(value) else None
    reference = model['partial_r']
    return dict(
        **{key: model[key] for key in ('filter', 'support_group', 'relationship', 'predictor', 'outcome', 'component', 'radial_bin', 'covariates')},
        omitted_field=field, status=status, full_n_lenses=len(cluster_ids), n_lenses=n,
        n_omitted=len(cluster_ids)-n, design_rank=rank, design_columns=design.shape[1],
        retained_lens_ids='|'.join(retained), full_partial_r=reference, omission_partial_r=estimate,
        delta_partial_r=None if estimate is None else estimate-reference,
        sign_flip=None if estimate is None or reference == 0 or estimate == 0 else np.signbit(estimate).item() != np.signbit(reference).item(),
    )


def run(output):
    if output.exists():
        raise FileExistsError(output)
    _, metrics, joined, icl_by_cluster = load_and_check()
    ids = sorted(icl_by_cluster)
    lens_fields, geometry_digests = assign_lens_fields(ids)

    # Reproduce the original models and partial r before restricting any rows.
    targets = [{f'{angle}_deg': float(icl_by_cluster[i][f'{angle}_deg']) for angle in ('ra', 'dec')} for i in ids]
    controls = [{f'{angle}_deg': float(icl_by_cluster[i][f'control_{angle}_deg']) for angle in ('ra', 'dec')} for i in ids]
    _, labels, _, _, _, _ = pair_block_labels(targets, controls, np.ones(len(ids), dtype=bool))
    models = make_models(metrics, joined, icl_by_cluster, dict(zip(ids, map(str, labels), strict=True)))
    fields = sorted(set(lens_fields.values()))
    rows = [omit_field(model, lens_fields, field) for model in models for field in fields]
    if len(models) != 116 or len(rows) != 348:
        raise ValueError('Omission coverage differs from the fixed 116-test family and three fields')
    changed = [row for row in rows if row['status'] == 'complete']
    if len(changed) == 0:
        raise ValueError('No evaluable field omission')

    samples = {}
    for model in models:
        sample_key = f"{model['filter']}/{model['support_group']}"
        field_counts = dict(Counter(lens_fields[i] for i in model['cluster_ids']))
        if samples.setdefault(sample_key, field_counts) != field_counts:
            raise ValueError('Sample field totals vary by relation')

    output.mkdir(parents=True, exist_ok=False)
    table = output / 'field_omission_relations.tsv'
    membership = output / 'lens_field_membership.tsv'
    write_tsv(table, rows)
    write_tsv(membership, [dict(cluster_id=i, field=lens_fields[i], icl_filter=icl_by_cluster[i]['filter']) for i in ids])
    summary = dict(
        status='complete_lens_field_omission_influence', science_measurement=False,
        definition='Omit lenses by official parent footprint field and refit residual partial r on saved profiles',
        lens_fields=dict(Counter(lens_fields.values())), sample_fields=samples,
        full_family_size=len(models), omission_rows=len(rows), original_model_lens_floor=MIN_LENSES,
        omission_status_counts=dict(Counter(row['status'] for row in rows)),
        sign_reversals=sum(bool(row['sign_flip']) for row in changed),
        largest_changes=sorted(changed, key=lambda row: abs(row['delta_partial_r']), reverse=True)[:5],
        direct_icl_satellite_omissions=[row for row in rows if row['relationship'] == 'icl_vs_satellite'],
        limitations=['Retained lenses use their original source and matched-control endpoints.',
                     'Three fields give an influence check without a survey covariance or detection test.',
                     'Filtered subsets follow the original minimum of 12 lenses and require full model rank.'],
        input_sha256={**EXPECTED_SHA256, 'official_field_summary': FIELD_SHA, 'target_manifest': TARGET_SHA,
                      'official_geometry_files': geometry_digests},
        runner_sha256=sha256(Path(__file__)),
        output_sha256={path.name: sha256(path) for path in (table, membership)},
    )
    summary_path = output / 'summary.json'
    with summary_path.open('x') as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({key: summary[key] for key in ('lens_fields', 'sample_fields', 'omission_status_counts', 'sign_reversals')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
