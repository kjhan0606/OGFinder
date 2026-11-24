#!/usr/bin/env python3
"""Compare S16A pilot and COSMOS point-z selections and lens responses offline."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_pointz_background_support import (
    GAPS, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table,
)
from ogfkit.wl_photoz import point_redshift_response
from ogfkit.wlcluster import Cosmology

INTAKE = PLAN.parent / 's16a_ephor_ab_cosmos_intake_2026-11-19_01'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_path = INTAKE / 'summary.json'
    intake = json.loads(summary_path.read_text())
    table = INTAKE / 'cosmos_join.tsv'
    if sha(table) != intake['output_sha256'][table.name]:
        raise ValueError('COSMOS intake table hash differs')
    with table.open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t')
                if r['eligible_for_selection_review'] == 'True']
    if len(rows) != intake['eligible_rows']:
        raise ValueError('COSMOS eligibility count differs')
    columns = ['photo_z', 'reference_z', 'shape_weight', 'som_weight', 'photo_z_risk']
    photo_z, reference, shape, som, risk = np.array([[float(r[k]) for k in columns] for r in rows]).T
    tracts = np.array([r['tract'] for r in rows])
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    point_z = finite_point_redshifts(available)
    observed = defaultdict(list)
    for pair in pairs:
        obj = pair['object_id']
        if obj in selected and obj in point_z:
            observed[pair['center_id']].append((point_z[obj], float(pair['weight'])))
    cosmology = Cosmology()
    results = []
    lenses = [c for c in centers if c['role'] == 'lens']
    if len(lenses) != 11:
        raise ValueError('S16A pilot lens count differs')
    for lens in lenses:
        z_lens = float(lens['z_lens'])
        chi_lens = cosmology.comoving(z_lens)
        obs_z, obs_shape = np.asarray(observed[lens['center_id']], dtype=float).reshape(-1, 2).T
        obs_background = obs_z > z_lens
        observed_weights = np.zeros_like(obs_z)
        observed_weights[obs_background] = obs_shape[obs_background] * (
            1 - chi_lens / cosmology.comoving(obs_z[obs_background])
        ) ** 2
        for gap in GAPS:
            obs_keep = obs_z > z_lens + gap
            for footprint in ('combined', *sorted(set(tracts))):
                region = np.ones(len(tracts), dtype=bool) if footprint == 'combined' else tracts == footprint
                base = region & np.isfinite(risk) & (photo_z > z_lens)
                keep = base & (photo_z > z_lens + gap)
                response = point_redshift_response(
                    z_lens, photo_z[keep], reference[keep], shape[keep], som[keep], cosmology,
                )
                if gap == 0:
                    baseline = response['weight_sum']
                else:
                    baseline = point_redshift_response(
                        z_lens, photo_z[base], reference[base], shape[base], som[base], cosmology,
                    )['weight_sum']
                results.append(dict(
                    lens_id=lens['matched_lens_id'], lens_z=z_lens, gap=gap, cosmos_tract=footprint,
                    observed_pairs=int(obs_keep.sum()),
                    observed_retained_weight=float(observed_weights[obs_keep].sum()/observed_weights.sum())
                    if observed_weights.sum() else None,
                    cosmos_retained_weight=response['weight_sum']/baseline, **response,
                ))
    output.mkdir()
    path = output / 'lens_response.tsv'
    write_table(path, results)
    report = dict(
        science_measurement=False, calibration_applied=False, compatibility='partial_only',
        rule='Finite Ephor_AB best/risk followed by strict z_best > z_lens + gap, no risk threshold.',
        limitations=[
            'HSC COSMOS reference weights target baseline S16A shapes. Exact local shape selection, '
            'accepted pixels and missing-photo-z selection cannot be reproduced on the joined COSMOS rows.',
            'The S17A Ephor_AB reference uses the same named method, without a proven release-level '
            'equivalence of its outputs to S16A.',
            'Tract splits in the same COSMOS field test sensitivity and do not give independent-field errors.',
            'Responses use reference redshifts as a proxy, and do not model cluster-member dilution or a PDF cut.',
        ],
        n_lenses=len(lenses), n_rows=len(results), cosmology={'H0':70, 'Om0':0.3},
        runner_sha256=sha(Path(__file__)), cosmos_intake_sha256=sha(summary_path),
        output_sha256={path.name:sha(path)}, **provenance,
    )
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
