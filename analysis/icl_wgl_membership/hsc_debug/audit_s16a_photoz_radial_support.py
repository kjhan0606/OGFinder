#!/usr/bin/env python3
"""Audit finite Ephor_AB input coverage in the fixed S16A pilot annuli.

Use the verified receipt and the existing Table 4 diagnostic selection.
Finite best/risk values describe input availability. No redshift cut, lensing
efficiency, or science measurement is computed.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from audit_s23b_s16a_quality_supplement import OUTPUT as OBJECT_QA, common_support
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES, radial_bins, separations_radians, weight_diagnostics
from audit_s23b_s16a_wl_random_support import write_table
from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_s23b_s16a_quality_supplement import load_pilot
from s16a_y1_source_contract import finite
from submit_s23b_s16a_quality_supplement import APPROVED_SHA, OUTPUT as PLAN, load_plan, verify_complete


def load_inputs(plan_dir, qa_dir):
    plan = load_plan(plan_dir)
    results = plan_dir / 'results'
    receipt_file = results / 'download_summary.json'
    status = json.loads((results / 'quality_supplement_status.json').read_text())
    receipt_sha = sha(receipt_file)
    if status['status'] != 'complete' or status['summary_sha256'] != receipt_sha:
        raise ValueError('supplement receipt is incomplete or changed')
    receipt = verify_complete(plan, plan_dir, status)
    summary_file = qa_dir / 'summary.json'
    summary = json.loads(summary_file.read_text())
    if summary['receipt_sha256'] != receipt_sha or summary['approved_plan_sha256'] != APPROVED_SHA:
        raise ValueError('object QA refers to a different receipt or plan')

    def qa_rows(name):
        path = qa_dir / name
        if sha(path) != summary['output_sha256'][name]:
            raise ValueError(f'changed object QA file: {name}')
        with path.open(newline='') as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    selected = set()
    selection_ids = set()
    for row in qa_rows('unique_object_selection.tsv'):
        ident = row['object_id']
        if ident in selection_ids or row['selected_diagnostic'] not in ('True', 'False'):
            raise ValueError('invalid object selection ledger')
        selection_ids.add(ident)
        if row['selected_diagnostic'] == 'True':
            selected.add(ident)

    available = set()
    received = set()
    for name in receipt['output_sha256']:
        with (results / name).open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                ident = row['object_id']
                received.add(ident)
                if finite(row['photoz_best']) is not None and finite(row['photoz_risk_best']) is not None:
                    available.add(ident)

    pilot, pairs, object_ids, _ = load_pilot(Path(plan['pilot_directory']))
    if received != selection_ids or received != set(object_ids):
        raise ValueError('object IDs differ between pilot, supplement and selection ledger')
    return pilot['centers'], pairs, selected, available, qa_rows('per_center_radial_support.tsv'), {
        'plan_sha256': APPROVED_SHA, 'receipt_sha256': receipt_sha, 'object_qa_sha256': sha(summary_file),
    }


def coverage(centers, pairs, selected, available):
    by_center = defaultdict(list)
    for row in pairs:
        if row['object_id'] in selected:
            by_center[row['center_id']].append(row)

    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    profiles = []
    for center in centers:
        center_id = center['center_id']
        rows = by_center[center_id]
        ra = [float(row['source_ra']) for row in rows]
        dec = [float(row['source_dec']) for row in rows]
        weights = np.array([float(row['weight']) for row in rows])
        has_photoz = np.array([row['object_id'] in available for row in rows], dtype=bool)
        radius = separations_radians(ra, dec, center['ra_deg'], center['dec_deg'])
        bins = radial_bins(radius * cosmo.angular_diameter_distance(center['z_lens']).value)

        for radial_bin, (low, high) in enumerate(zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])):
            in_bin = bins == radial_bin
            complete = in_bin & has_photoz
            total = int(in_bin.sum())
            count = int(complete.sum())
            before = weight_diagnostics(weights[in_bin])
            after = weight_diagnostics(weights[complete])
            profiles.append(dict(
                center_id=center_id, lens_id=center['matched_lens_id'], role=center['role'],
                icl_filter=center['icl_filter'], radial_bin=radial_bin, r_min_mpc=float(low), r_max_mpc=float(high),
                selected_pairs=total, finite_photoz_pairs=count, missing_photoz_pairs=total - count,
                pair_fraction=count / total if total else None,
                shape_weight_fraction=after['sum_shape_weight'] / before['sum_shape_weight'] if total else None,
                selected_shape_weight=before['sum_shape_weight'], finite_photoz_shape_weight=after['sum_shape_weight'],
                selected_neff=before['effective_number_shape_weight'], finite_photoz_neff=after['effective_number_shape_weight'],
            ))
    return profiles


def group_support(profiles, centers):
    # The existing support reducer compares a lens and its five fixed randoms.
    records = [dict(center_id=row['center_id'], radial_bin=row['radial_bin'],
                    raw_count=row['selected_pairs'], selected_count=row['finite_photoz_pairs']) for row in profiles]
    return [dict(lens_id=r['matched_lens_id'], icl_filter=r['icl_filter'], radial_bin=r['radial_bin'],
                 r_min_mpc=r['r_min_mpc'], r_max_mpc=r['r_max_mpc'],
                 target_selected_pairs=r['target_raw_count'], target_finite_photoz_pairs=r['target_selected_count'],
                 randoms_with_selected=r['randoms_with_raw_rows'], randoms_with_finite_photoz=r['randoms_with_selected_rows'],
                 common_selected=r['raw_common_support'], common_finite_photoz=r['selected_common_support'])
            for r in common_support(records, centers)]


def plot_coverage(profiles, centers, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 10})
    labels = []
    matrix = []
    values = {(r['center_id'], r['radial_bin']): r['shape_weight_fraction'] for r in profiles}
    lens_ids = sorted({c['matched_lens_id'] for c in centers})
    for lens_id in lens_ids:
        group = [c for c in centers if c['matched_lens_id'] == lens_id]
        group.sort(key=lambda c: c['role'] != 'lens')
        for index, center in enumerate(group):
            labels.append(lens_id.replace('S23B_', '') + ' lens' if index == 0 else '')
            matrix.append([values[center['center_id'], k] for k in range(len(RADIAL_EDGES)-1)])

    fig, ax = plt.subplots(figsize=(9.1, 8), layout='constrained')
    colormap = plt.colormaps['viridis'].copy()
    colormap.set_bad('#e3e3e3')
    image = ax.imshow(np.array(matrix, dtype=float), aspect='auto', interpolation='nearest',
                      cmap=colormap, vmin=0, vmax=1)
    for boundary in range(6, len(matrix), 6):
        ax.axhline(boundary-.5, color='white', linewidth=.8)
    ax.set(xticks=range(len(RADIAL_EDGES)-1),
           xticklabels=[f'{lo:g}–{hi:g}' for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])],
           yticks=range(len(labels)), yticklabels=labels,
           xlabel='Proper radius [Mpc]', ylabel='Lens and five fixed randoms')
    ax.tick_params(axis='y', length=0)
    fig.colorbar(image, ax=ax, label='Shape weight fraction with finite Ephor_AB inputs')
    # Reserve an empty strip inside the panel for the title.
    ax.set_ylim(len(labels)-.5, -4.5)
    ax.text(.02, .98, 'Ephor_AB input availability by center and annulus', transform=ax.transAxes, va='top')
    fig.savefig(output / 'photoz_radial_coverage.png', dpi=160)
    fig.savefig(output / 'photoz_radial_coverage.pdf')
    plt.close(fig)


def run(plan_dir, qa_dir, output):
    if output.exists():
        raise ValueError('output already exists')
    centers, pairs, selected, available, prior, provenance = load_inputs(plan_dir, qa_dir)
    profiles = coverage(centers, pairs, selected, available)
    expected = {(r['center_id'], int(r['radial_bin'])): int(r['selected_count']) for r in prior}
    observed = {(r['center_id'], r['radial_bin']): r['selected_pairs'] for r in profiles}
    if observed != expected or len(profiles) != len(prior):
        raise ValueError('annular counts differ from verified object QA')

    common = group_support(profiles, centers)
    finite_pairs = sum(r['finite_photoz_pairs'] for r in profiles)
    selected_pairs = sum(r['selected_pairs'] for r in profiles)
    if selected_pairs != sum(r['object_id'] in selected for r in pairs):
        raise ValueError('selected source pair totals differ')
    if finite_pairs != sum(r['object_id'] in selected and r['object_id'] in available for r in pairs):
        raise ValueError('finite photo-z pair totals differ')

    annuli = []
    for k, (low, high) in enumerate(zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])):
        rows = [r for r in common if r['radial_bin'] == k]
        annuli.append(dict(r_min_mpc=float(low), r_max_mpc=float(high),
                           common_selected=sum(r['common_selected'] for r in rows),
                           common_finite_photoz=sum(r['common_finite_photoz'] for r in rows)))
    intersection = {}
    for minimum in (0., .1, .7):
        comparison = defaultdict(list)
        for row in common:
            if row['r_min_mpc'] >= minimum:
                comparison[row['lens_id']].append(row)
        intersection[str(minimum)] = {key: sorted(lens for lens, rows in comparison.items() if all(r[key] for r in rows))
                                     for key in ('common_selected', 'common_finite_photoz')}

    positive = [r for r in profiles if r['selected_pairs']]
    absent = [r for r in positive if not r['finite_photoz_pairs']]
    subtotal = {}
    for role in ('lens', 'random'):
        rows = [r for r in profiles if r['role'] == role]
        total_weight = sum(r['selected_shape_weight'] for r in rows)
        known_weight = sum(r['finite_photoz_shape_weight'] for r in rows)
        subtotal[role] = dict(selected_pairs=sum(r['selected_pairs'] for r in rows),
                             finite_photoz_pairs=sum(r['finite_photoz_pairs'] for r in rows),
                             shape_weight_fraction=known_weight / total_weight)

    output.mkdir(parents=True, exist_ok=False)
    write_table(output / 'center_annulus_coverage.tsv', profiles)
    write_table(output / 'lens_annulus_common_support.tsv', common)
    plot_coverage(profiles, centers, output)
    result = dict(schema='astrafek-s16a-photoz-radial-support-v1', **provenance,
                  source_code_sha256=sha(Path(__file__)), center_annuli=len(profiles),
                  selected_unique_objects=len(selected), finite_photoz_unique_objects=len(selected & available),
                  selected_pairs=selected_pairs, finite_photoz_pairs=finite_pairs, role_pair_diagnostics=subtotal,
                  positive_center_annuli=len(positive), zero_photoz_positive_annuli=absent,
                  lowest_weight_fraction_annuli=sorted(positive, key=lambda r: r['shape_weight_fraction'])[:5],
                  common_support=annuli, common_group_intersection=intersection,
                  photoz_used_as_selection=False, science_measurement=False,
                  output_sha256={p.name: sha(p) for p in sorted(output.iterdir())})
    with (output / 'summary.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('selected_pairs', 'finite_photoz_pairs', 'role_pair_diagnostics',
                                            'common_support', 'common_group_intersection')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-dir', type=Path, default=PLAN)
    parser.add_argument('--object-qa', type=Path, default=OBJECT_QA)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.plan_dir.resolve(), args.object_qa.resolve(), args.output.resolve())
