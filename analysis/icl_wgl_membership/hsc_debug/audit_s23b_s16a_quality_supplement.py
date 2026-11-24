#!/usr/bin/env python3
"""Offline object-cut, tail, radial-support and point-photo-z QA of S16A.

Table-4 cuts are diagnostics, not a calibrated background-source selection.
All original centers including prior-count zeros are retained. No shear, mass,
boost, source-z correction or covariance is computed.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from fetch_ancillary_catalogues import catalogue_reader, sha
from s16a_y1_source_contract import select_row, merge_supplement, finite
from submit_s23b_s16a_quality_supplement import load_plan, verify_complete, OUTPUT as PLAN
from prepare_s23b_s16a_quality_supplement import load_pilot
from audit_s23b_s16a_shape_qa_pilot import radial_bins, RADIAL_EDGES, separations_radians, weight_diagnostics
from audit_s23b_s16a_wl_random_support import describe, write_table

OUTPUT = PLAN.parent/'s16a_quality_supplement_qa_2026-11-18_01'
DIAGNOSTIC = 'Table-4 object cuts and finite calibration/positive shape weight only; no calibrated background-source selection'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def common_support(profiles, centers):
    groups = defaultdict(list)
    for center in centers:
        groups[center['matched_lens_id']].append(center)
    lookup = {(r['center_id'], r['radial_bin']):r for r in profiles}
    require(len(lookup) == len(profiles), 'duplicate radial profile')
    common = []
    for lens, group in sorted(groups.items()):
        targets = [c for c in group if c['role'] == 'lens']
        randoms = [c for c in group if c['role'] == 'random']
        require(len(targets) == 1 and len(randoms) == 5, 'expected one lens and five fixed randoms')
        for k in range(len(RADIAL_EDGES)-1):
            target = lookup[(targets[0]['center_id'], k)]
            rs = [lookup[(c['center_id'], k)] for c in randoms]
            common.append(dict(matched_lens_id=lens, icl_filter=group[0]['icl_filter'], radial_bin=k,
                r_min_mpc=float(RADIAL_EDGES[k]), r_max_mpc=float(RADIAL_EDGES[k+1]),
                target_raw_count=target['raw_count'], target_selected_count=target['selected_count'],
                randoms_with_raw_rows=sum(r['raw_count'] > 0 for r in rs),
                randoms_with_selected_rows=sum(r['selected_count'] > 0 for r in rs),
                raw_common_support=target['raw_count'] > 0 and all(r['raw_count'] > 0 for r in rs),
                selected_common_support=target['selected_count'] > 0 and all(r['selected_count'] > 0 for r in rs)))
    return common


def numeric_calibration(row):
    values = [finite(row.get(k)) for k in ('sigma_e', 'rms_e', 'weight', 'm', 'c1', 'c2')]
    return all(v is not None for v in values) and values[2] > 0


def analyze(centers, records):
    center_map = {c['center_id']:c for c in centers}
    require(len(center_map) == len(centers), 'duplicate center')
    by_center, objects = defaultdict(list), {}
    pairs, pair_status = set(), Counter()
    for row in records:
        key = (row['center_id'], row['object_id'])
        require(key not in pairs and row['center_id'] in center_map, 'duplicate pair or unknown center')
        pairs.add(key); by_center[row['center_id']].append(row)
        obj = row['object_id']
        if obj in objects:
            require(all(row[k] == objects[obj]['row'][k] for k in row if k != 'center_id'), 'source properties differ between reused pairs')
        else:
            selection = select_row(row)
            objects[obj] = dict(row=row, selection=selection, numeric=numeric_calibration(row), pairs=0)
        objects[obj]['pairs'] += 1
        pair_status[objects[obj]['selection']['table4_diagnostic']] += 1
    unique_status = Counter(v['selection']['table4_diagnostic'] for v in objects.values())
    criteria = []
    if objects:
        for name in next(iter(objects.values()))['selection']['criteria']:
            unique, paired = Counter(), Counter()
            for item in objects.values():
                status = item['selection']['criteria'][name]
                unique[status] += 1; paired[status] += item['pairs']
            criteria.append(dict(criterion=name, **{f'unique_{s}':unique[s] for s in ('pass','fail','unknown')},
                                 **{f'pair_{s}':paired[s] for s in ('pass','fail','unknown')}))
    selection_rows, tails = [], []
    for obj, item in sorted(objects.items(), key=lambda x:int(x[0])):
        row, sel = item['row'], item['selection']
        failures = '|'.join(k for k,v in sel['criteria'].items() if v == 'fail')
        unknowns = '|'.join(k for k,v in sel['criteria'].items() if v == 'unknown')
        selected = sel['table4_diagnostic'] == 'pass' and item['numeric']
        selection_rows.append(dict(object_id=obj, table4_diagnostic=sel['table4_diagnostic'],
            numeric_calibration=item['numeric'], selected_diagnostic=selected,
            pair_multiplicity=item['pairs'], failed_criteria=failures, unknown_criteria=unknowns))
        m, sigma = finite(row.get('m')), finite(row.get('sigma_e'))
        if (m is not None and m <= -1) or (sigma is not None and sigma > .4):
            tails.append(dict(object_id=obj, pair_multiplicity=item['pairs'], m=m,
                derived_sigma_e=sigma, raw_sigma=finite(row.get('raw_sigma')),
                resolution=finite(row.get('resolution')), m_le_minus_one=m is not None and m <= -1,
                derived_sigma_gt_0p4=sigma is not None and sigma > .4,
                table4_diagnostic=sel['table4_diagnostic'], failed_criteria=failures, unknown_criteria=unknowns))
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0)
    center_rows, radial_rows = [], []
    for center in centers:
        rows = by_center[center['center_id']]
        require(len(rows) == center['expected_raw_pairs'], 'center raw count differs')
        selected = np.array([objects[r['object_id']]['selection']['table4_diagnostic'] == 'pass'
                             and objects[r['object_id']]['numeric'] for r in rows], dtype=bool)
        radius = separations_radians([float(r['source_ra']) for r in rows],
            [float(r['source_dec']) for r in rows], center['ra_deg'], center['dec_deg'])
        radius *= float(cosmo.angular_diameter_distance(center['z_lens']).value)
        bins = radial_bins(radius)
        weights = np.array([float(r['weight']) for r in rows])
        status = Counter(objects[r['object_id']]['selection']['table4_diagnostic'] for r in rows)
        common = dict(center_id=center['center_id'], matched_lens_id=center['matched_lens_id'], role=center['role'],
                      icl_filter=center['icl_filter'], source_origin='received_rows' if center['query_planned'] else 'prior_count_only_zero_not_requeried')
        center_rows.append(dict(**common, raw_count=len(rows), table4_pass=status['pass'],
            table4_fail=status['fail'], table4_unknown=status['unknown'], selected_count=int(selected.sum()),
            **weight_diagnostics(weights[selected])))
        for k in range(len(RADIAL_EDGES)-1):
            inside = bins == k; keep = inside & selected
            radial_rows.append(dict(**common, radial_bin=k, r_min_mpc=float(RADIAL_EDGES[k]),
                r_max_mpc=float(RADIAL_EDGES[k+1]), raw_count=int(inside.sum()),
                selected_count=int(keep.sum()), **weight_diagnostics(weights[keep])))
    common = common_support(radial_rows, centers)
    require(sum(r['raw_count'] for r in radial_rows) == len(records), 'radial raw accounting differs')
    require(sum(r['selected_count'] for r in radial_rows) == sum(v['pairs'] for v in objects.values()
        if v['selection']['table4_diagnostic'] == 'pass' and v['numeric']), 'radial selected accounting differs')
    photoz = {}
    for name in ('photoz_best', 'photoz_risk_best'):
        for scope, items in (('all_unique', list(objects.values())),
                             ('table4_pass_unique', [v for v in objects.values() if v['selection']['table4_diagnostic'] == 'pass'])):
            values = [finite(item['row'].get(name)) for item in items]
            good = [v for v in values if v is not None]
            photoz[name+'_'+scope] = dict(total_objects=len(values), missing_or_nonfinite=len(values)-len(good),
                finite_distribution=describe(good), negative=sum(v < 0 for v in good),
                above_one=sum(v > 1 for v in good) if name == 'photoz_risk_best' else None)
    intersection = {}
    for start, label in ((0., 'all_0_to_3'), (.1, 'all_0p1_to_3'), (.7, 'all_0p7_to_3')):
        relevant = [r for r in common if r['r_min_mpc'] >= start]
        groups = defaultdict(list)
        for row in relevant: groups[row['matched_lens_id']].append(row)
        intersection[label] = [lens for lens, group in sorted(groups.items()) if all(r['selected_common_support'] for r in group)]
    support = [dict(radial_bin=k, r_min_mpc=float(RADIAL_EDGES[k]),r_max_mpc=float(RADIAL_EDGES[k+1]),
                    raw_common_groups=sum(r['raw_common_support'] for r in common if r['radial_bin']==k),
                    selected_common_groups=sum(r['selected_common_support'] for r in common if r['radial_bin']==k))
               for k in range(len(RADIAL_EDGES)-1)]
    return dict(unique_status=dict(unique_status), pair_status=dict(pair_status), criteria=criteria,
                selection_rows=selection_rows, tails=tails, centers=center_rows, radial=radial_rows,
                common=common, photoz=photoz, radial_common_support=support, common_group_intersection=intersection)


def make_plot(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    counts = [result['unique_status'].get(k,0) for k in ('pass','fail','unknown')]
    axes[0].bar(['Pass', 'Fail', 'Unknown'], counts, color=['#2a8f77','#be6a4d','#8a8a8a'])
    axes[0].set_ylabel('Unique objects'); axes[0].set_ylim(0,max(counts)*1.25 if max(counts) else 1)
    axes[0].text(.03,.97,'Table-4 object diagnostics',va='top',transform=axes[0].transAxes)
    for x,n in enumerate(counts): axes[0].text(x,n+max(counts)*.025,str(n),ha='center')
    support = result['radial_common_support']
    x = [(r['r_min_mpc']+r['r_max_mpc'])/2 for r in support]
    axes[1].step(x,[r['raw_common_groups'] for r in support],where='mid',label='Raw support')
    axes[1].step(x,[r['selected_common_groups'] for r in support],where='mid',label='Object cuts + numeric calibration')
    axes[1].set_xlabel('Proper radius [Mpc]');axes[1].set_ylabel('Lens groups with all 5 randoms')
    axes[1].set_ylim(0,12);axes[1].text(.03,.97,'Diagnostic radial support',va='top',transform=axes[1].transAxes)
    axes[1].legend(loc='upper right',bbox_to_anchor=(1,.86),fontsize=11)
    fig.tight_layout();fig.savefig(output/'selection_qa.png',dpi=160);fig.savefig(output/'selection_qa.pdf');plt.close(fig)


def main(output):
    require(not output.exists(), 'refusing to overwrite QA output')
    plan = load_plan(PLAN)
    state = json.loads((PLAN/'results/quality_supplement_status.json').read_text())
    require(state['status'] == 'complete' and sha(PLAN/'results/download_summary.json') == state['summary_sha256'], 'supplement not complete/verified')
    receipt = verify_complete(plan, PLAN, state)
    pilot, pilot_rows, _, _ = load_pilot(Path(plan['pilot_directory']))
    add = []
    for batch in plan['batches']:
        with (PLAN/'results'/(Path(batch['name']).stem+'.csv')).open() as f:
            add.extend(dict(r) for r in catalogue_reader(f))
    records = merge_supplement(pilot_rows, add)
    result = analyze(pilot['centers'], records)
    output.mkdir(parents=True, exist_ok=False)
    for key, name in [('criteria','criterion_attrition.tsv'),('selection_rows','unique_object_selection.tsv'),
                      ('centers','per_center_selection.tsv'),('radial','per_center_radial_support.tsv'),
                      ('common','per_lens_common_support.tsv')]:
        write_table(output/name, result[key])
    if result['tails']: write_table(output/'tail_selection_diagnostics.tsv', result['tails'])
    (output/'point_photoz_qa.json').write_text(json.dumps(result['photoz'],indent=2,allow_nan=False)+'\n')
    make_plot(result, output)
    summary = dict(schema='astrafek-s16a-quality-supplement-qa-v1', status='local_diagnostics_complete',
        diagnostic_definition=DIAGNOSTIC, n_unique_objects=len(add), n_source_pairs=len(records), n_centers=66,
        known_zero_controls_retained=23, unique_object_table4_status=result['unique_status'],
        source_pair_table4_status=result['pair_status'], radial_common_support=result['radial_common_support'],
        common_group_intersection=result['common_group_intersection'],
        tail_objects=len(result['tails']), tail_objects_with_table4_pass=sum(r['table4_diagnostic']=='pass' for r in result['tails']),
        verified_csv_bytes=receipt['verified_csv_bytes'], receipt_sha256=sha(PLAN/'results/download_summary.json'),
        approved_plan_sha256=plan['query_plan_sha256'], source_code_sha256=sha(Path(__file__)),
        output_sha256={f.name:sha(f) for f in sorted(output.iterdir()) if f.is_file()},
        new_query_submitted=False, point_photoz_used_for_selection=False, calibrated_background_selection=False,
        science_measurement=False, limitations=['Deliberately stratified 11-lens pilot, not representative retention.',
            'Object cuts do not establish usable area, background calibration, PSF null tests or covariance.',
            'Shape-weight effective number is descriptive; random source overlap remains.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    main(parser.parse_args().output.resolve())
