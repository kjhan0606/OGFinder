#!/usr/bin/env python3
"""Audit saved S16A target/random counts without source rows or network use.

The original 114-lens/570-center selection is preserved. Count-floor and
count-ratio checks are descriptive sensitivities, not new scientific cuts.
PDR3 polygon labels describe the parent imaging field, not S16A usable area.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from prepare_s23b_s16a_wl_random_count_plan import (
    OUTPUT as RANDOM_PLAN, RANDOM_CENTERS, TARGET_SELECTION, TARGET_MANIFEST,
    SELECTION_SUMMARY,
)
from submit_s23b_s16a_wl_random_counts import validate_plan, read_count_csv
from hsc_pdr3_tract_patch_geometry import load_index

DATA = RANDOM_PLAN.parent
TARGET_RESULTS = DATA/'wl_s23b_s16a_shape_count_plan_2026-11-11_01/results'
VERIFICATION = DATA/'wl_s23b_s16a_random_count_verification_2026-11-17_01/verification.json'
GEOMETRY = DATA/'control_patches_2026-10-14_02/public_geometry'
DEFAULT_OUTPUT = DATA/'s23b_s16a_wl_random_support_audit_2026-11-17_01'
PINNED = {
    TARGET_RESULTS/'count_results.tsv': '0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12',
    TARGET_RESULTS/'count_summary.json': '0cf90929c2a69eda0e62f3a9350359fd430b37bf869b2ead53246129ac51bb8b',
    RANDOM_PLAN/'results/count_results.tsv': '91f7521d49dc8c734fadf375f80fe58d610ff72b851569c2c38e3a5dc6eb2b77',
    RANDOM_PLAN/'results/count_summary.json': 'c3a92b0420c5c8af8cb6191ec1264305e3fc38e503365b456f5dcad3f07332ef',
    VERIFICATION: '8394aeeda9b1d4f885b2a884bd7be14fde9278ed2dcf750aebe8d601d937578f',
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_tsv(path, columns=None):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        require(reader.fieldnames is not None and len(set(reader.fieldnames)) == len(reader.fieldnames),
                'missing or duplicate table header')
        if columns is not None:
            require(reader.fieldnames == columns, 'table columns changed')
        rows = list(reader)
    require(all(None not in r and all(v is not None for v in r.values()) for r in rows), 'malformed table row')
    return rows


def count_map(rows, key):
    result = {}
    for r in rows:
        value = int(r['n_shape_sources'])
        require(value >= 0 and r[key] and r[key] not in result, 'negative count or duplicate ID')
        result[r[key]] = value
    return result


def describe(values):
    a = np.asarray(values, dtype=float)
    require(a.ndim == 1 and np.isfinite(a).all(), 'nonfinite descriptive input')
    if not len(a):
        return dict(n=0, min=None, q16=None, median=None, q84=None, max=None)
    return dict(n=len(a), **dict(zip(('min', 'q16', 'median', 'q84', 'max'),
                                    map(float, np.quantile(a, [0, .16, .5, .84, 1])))))


def angular_separation_arcmin(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = np.deg2rad([ra1, dec1, ra2, dec2])
    h = np.sin((dec2-dec1)/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin((ra2-ra1)/2)**2
    return float(np.rad2deg(2*np.arcsin(np.sqrt(np.clip(h, 0, 1))))*60)


def join_centers(selected, counts, target_counts, targets):
    ids = [r['random_center_id'] for r in selected]
    require(len(ids) == len(set(ids)) and set(ids) == set(counts), 'random count IDs differ from selected centers')
    objects, positions, slots = set(), set(), defaultdict(set)
    output = []
    for r in selected:
        ident, lens = r['random_center_id'], r['matched_lens_id']
        require(lens in target_counts and lens in targets, 'unknown lens')
        target_count = target_counts[lens]
        require(target_count >= 1000 and int(r['target_s16a_raw_count']) == target_count, 'target count disagrees with original query')
        count = counts[ident]
        require(isinstance(count, int) and not isinstance(count, bool) and count >= 0, 'invalid random count')
        z, ra, dec = (float(r[k]) for k in ('z_lens', 'ra_deg', 'dec_deg'))
        require(all(math.isfinite(x) for x in (z, ra, dec)) and 0 <= ra < 360 and -90 <= dec <= 90, 'invalid coordinates')
        require(math.isclose(z, float(targets[lens]['z_cl']), rel_tol=0, abs_tol=1e-12), 'lens redshift differs')
        require(r['icl_filter'] in ('i', 'i2'), 'unknown ICL filter')
        obj, pos, slot = r['random_object_id'], (round(ra, 12), round(dec, 12)), int(r['control_slot'])
        require(obj not in objects and pos not in positions, 'random position/object reused')
        require(slot in range(1, 6) and slot not in slots[lens], 'invalid or duplicate control slot')
        objects.add(obj); positions.add(pos); slots[lens].add(slot)
        meta = {k: float(r[k]) for k in ('offset_from_target_mpc', 'offset_from_target_arcmin',
                    'match_score', 'delta_depth_mag', 'seeing_ratio', 'delta_ebv_sfd')}
        require(all(math.isfinite(x) for x in meta.values()) and meta['seeing_ratio'] > 0, 'nonfinite matching metadata')
        require(meta['offset_from_target_mpc'] > 6 and int(r['cluster_overlap_2mpc']) == 0, 'saved candidate cuts differ')
        sep = angular_separation_arcmin(float(targets[lens]['ra_deg']), float(targets[lens]['dec_deg']), ra, dec)
        require(math.isclose(sep, meta['offset_from_target_arcmin'], rel_tol=1e-6, abs_tol=1e-6), 'saved angular separation disagrees with target coordinates')
        output.append(dict(random_center_id=ident, matched_lens_id=lens, control_slot=slot,
            random_object_id=obj, icl_filter=r['icl_filter'], ra_deg=ra, dec_deg=dec, z_lens=z,
            target_raw_count=target_count, random_raw_count=count, raw_count_ratio=count/target_count,
            zero_raw_count=count == 0, candidate_pdr3_tract=int(r['tract']), candidate_pdr3_patch=r['patch'],
            **meta))
    require(all(v == {1, 2, 3, 4, 5} for v in slots.values()), 'each lens needs exactly five controls')
    return sorted(output, key=lambda r: (r['matched_lens_id'], r['control_slot']))


def lens_rows(centers):
    groups = defaultdict(list)
    for r in centers:
        groups[r['matched_lens_id']].append(r)
    result = []
    for lens, group in sorted(groups.items()):
        require(len(group) == 5 and {r['control_slot'] for r in group} == {1, 2, 3, 4, 5}, 'invalid five-center group')
        require(len({(r['target_raw_count'], r['icl_filter'], r['z_lens']) for r in group}) == 1, 'inconsistent lens metadata')
        counts = [r['random_raw_count'] for r in group]
        ratios = [r['raw_count_ratio'] for r in group]
        result.append(dict(matched_lens_id=lens, icl_filter=group[0]['icl_filter'],
            lens_pdr3_field=group[0].get('lens_pdr3_field', 'UNASSIGNED'), z_lens=group[0]['z_lens'],
            target_raw_count=group[0]['target_raw_count'], n_random_centers=5,
            n_positive=sum(n > 0 for n in counts), n_ge1000=sum(n >= 1000 for n in counts),
            n_ratio_0p5_to_2=sum(.5 <= r <= 2 for r in ratios),
            min_random_count=min(counts), median_random_count=float(np.median(counts)),
            mean_random_count=sum(counts)/5, max_random_count=max(counts),
            mean_random_to_target=sum(ratios)/5, median_random_to_target=float(np.median(ratios)),
            min_random_to_target=min(ratios), max_random_to_target=max(ratios),
            pdr3_random_fields=','.join(sorted({r.get('pdr3_field', 'UNASSIGNED') for r in group})),
            n_other_or_unknown_pdr3_fields=sum(r.get('pdr3_field', 'UNASSIGNED') != group[0].get('lens_pdr3_field', 'UNASSIGNED')
                or r.get('pdr3_field', 'UNASSIGNED') == 'UNASSIGNED' for r in group)))
    return result


def summarize(centers, lenses):
    result = []
    for band in ('all', 'i', 'i2'):
        c = [r for r in centers if band == 'all' or r['icl_filter'] == band]
        l = [r for r in lenses if band == 'all' or r['icl_filter'] == band]
        result.append(dict(icl_filter=band, n_lenses=len(l), n_centers=len(c),
            zero_centers=sum(r['zero_raw_count'] for r in c),
            lenses_with_five_positive=sum(r['n_positive'] == 5 for r in l),
            lenses_by_positive_center_count=dict(sorted(Counter(r['n_positive'] for r in l).items())),
            target_count=describe([r['target_raw_count'] for r in l]),
            random_count=describe([r['random_raw_count'] for r in c]),
            center_ratio=describe([r['raw_count_ratio'] for r in c]),
            per_lens_mean_ratio=describe([r['mean_random_to_target'] for r in l]),
            zero_vs_positive_metadata={str(zero): {k: describe([r[k] for r in c if r['zero_raw_count'] == zero])
                for k in ('offset_from_target_mpc', 'match_score', 'delta_depth_mag', 'seeing_ratio', 'delta_ebv_sfd')}
                for zero in (False, True)}))
    return result


def sensitivities(lenses):
    result = []
    criteria = {'five_positive': lambda r: r['n_positive'] == 5,
                'five_ge1000': lambda r: r['n_ge1000'] == 5,
                'five_ratios_between_0p5_and_2': lambda r: r['n_ratio_0p5_to_2'] == 5,
                'five_ge1000_and_ratios_between_0p5_and_2': lambda r: r['n_ge1000'] == r['n_ratio_0p5_to_2'] == 5}
    for band in ('all', 'i', 'i2'):
        group = [r for r in lenses if band == 'all' or r['icl_filter'] == band]
        for name, criterion in criteria.items():
            result.append(dict(icl_filter=band, diagnostic=name, n_tested=len(group),
                n_matching=sum(criterion(r) for r in group), diagnostic_only=True))
    return result


def verify_inputs():
    plan, plan_hash = validate_plan(RANDOM_PLAN)
    for p, expected in PINNED.items():
        require(sha(p) == expected, f'pinned input changed: {p}')
    verified = json.loads(VERIFICATION.read_text())
    require(verified['status'] == 'verified' and verified['approved_plan_sha256'] == plan_hash, 'previous count verification differs')
    random_counts = count_map(read_tsv(RANDOM_PLAN/'results/count_results.tsv', ['random_center_id', 'n_shape_sources']), 'random_center_id')
    state = json.loads((RANDOM_PLAN/'results/count_status.json').read_text())
    summary = json.loads((RANDOM_PLAN/'results/count_summary.json').read_text())
    require(state['status'] == 'complete' and sha(RANDOM_PLAN/'results/count_summary.json') == state['summary_sha256'], 'random query not complete')
    require(state['plan_sha256'] == plan_hash and summary['query_plan_sha256'] == plan_hash, 'random query plan differs')
    original = {}
    for i, b in enumerate(plan['batches']):
        f = RANDOM_PLAN/'results'/(Path(b['name']).stem+'.csv')
        require(sha(f) == state['jobs'][f.stem]['sha256'], 'random CSV changed')
        ids = [str(r['random_center_id']) for r in plan['centers'][25*i:25*(i+1)]]
        for row in read_count_csv(f, ids):
            require(row['random_center_id'] not in original, 'duplicate batch ID')
            original[row['random_center_id']] = row['n_shape_sources']
    require(original == random_counts and len(original) == 570, 'combined random counts differ from original CSVs')
    target_counts = count_map(read_tsv(TARGET_RESULTS/'count_results.tsv', ['lens_id', 'n_shape_sources']), 'lens_id')
    target_summary = json.loads((TARGET_RESULTS/'count_summary.json').read_text())
    require(target_summary['status'] == 'complete' and sum(target_counts.values()) == target_summary['sum_lens_cone_shape_counts'], 'target count summary differs')
    targets = {str(r['id']): r for r in json.loads(TARGET_MANIFEST.read_text())['clusters']}
    require(set(target_counts) == set(targets) and len(targets) == 1043, 'target ID set differs')
    selection = read_tsv(TARGET_SELECTION)
    require(len(selection) == 116, 'original pre-query ledger differs')
    for r in selection:
        require(int(r['target_s16a_raw_count']) == target_counts[r['matched_lens_id']], 'selected target count differs')
    return plan, read_tsv(RANDOM_CENTERS), random_counts, target_counts, targets, selection


def assign_fields(centers, targets):
    ids = sorted({r['matched_lens_id'] for r in centers})
    points = [(float(targets[i]['ra_deg']), float(targets[i]['dec_deg'])) for i in ids]
    points.extend((r['ra_deg'], r['dec_deg']) for r in centers)
    index, geometry = load_index(GEOMETRY, level='tract')
    assigned = index.assign([p[0] for p in points], [p[1] for p in points], level='tract')
    fields = [str(x) if x is not None else 'UNASSIGNED' for x in assigned.field]
    lens_fields = dict(zip(ids, fields[:len(ids)]))
    for i, r in enumerate(centers, start=len(ids)):
        r.update(pdr3_field=fields[i], lens_pdr3_field=lens_fields[r['matched_lens_id']],
                 assigned_pdr3_tract=int(assigned.tract[i]),
                 containing_pdr3_polygons=int(assigned.containing_polygon_count[i]),
                 pdr3_cross_field_ambiguous=bool(assigned.cross_field_ambiguous[i]),
                 candidate_tract_equals_assigned=int(assigned.tract[i]) == r['candidate_pdr3_tract'])
    return dict(source_directory=str(GEOMETRY), summary=geometry,
        unmatched=assigned.unmatched_count, cross_field_ambiguous=assigned.cross_field_ambiguous_count,
        limitation='PDR3 imaging polygon labels are not S16A shear-catalog masks or aperture-area measurements.')


def write_table(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def plot(centers, lenses, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4))
    ax, sky = axes
    for band, color in (('i', '#2563eb'), ('i2', '#c2410c')):
        for complete, marker in ((True, 'o'), (False, 'x')):
            group = [r for r in lenses if r['icl_filter'] == band and (r['n_positive'] == 5) == complete]
            ax.scatter([r['target_raw_count'] for r in group], [r['mean_random_count'] for r in group],
                s=27, marker=marker, color=color, alpha=.8,
                label=f'{band}: '+('5 positive' if complete else '<5 positive'))
    limit = math.ceil(max(r['target_raw_count'] for r in lenses)/1000)*1000
    ax.plot([0, limit], [0, limit], color='.35', lw=1, ls='--', zorder=0)
    ax.set(xlim=(-100, limit), ylim=(-100, limit), xlabel=r'$N_\mathrm{lens}$', ylabel=r'$\overline{N}_\mathrm{random}$')
    ax.text(.04, .96, 'Raw counts; fixed five-center mean', transform=ax.transAxes, va='top')
    ax.legend(loc='upper left', bbox_to_anchor=(.02, .87), framealpha=.9, fontsize=11)
    for zero, color, marker, size in ((False, '#9ca3af', '.', 15), (True, '#b91c1c', 'x', 29)):
        group = [r for r in centers if r['zero_raw_count'] == zero]
        ras = [r['ra_deg']-360 if r['ra_deg'] > 300 else r['ra_deg'] for r in group]
        sky.scatter(ras, [r['dec_deg'] for r in group], color=color, marker=marker, s=size,
            alpha=.85, label=f'{"Zero" if zero else "Positive"}: {len(group)} centers')
    sky.set(xlabel='RA [deg; wrapped]', ylabel='Dec [deg]', ylim=(-12, 52))
    sky.text(.04, .96, 'S16A cone availability', transform=sky.transAxes, va='top')
    sky.legend(loc='upper left', bbox_to_anchor=(.02, .87), framealpha=.9, fontsize=11)
    fig.subplots_adjust(left=.08, right=.985, top=.97, bottom=.14, wspace=.28)
    for suffix in ('.png', '.pdf'):
        fig.savefig(destination.with_suffix(suffix), dpi=180)
    plt.close(fig)


def run(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite output directory')
    plan, selected, counts, target_counts, targets, original_selection = verify_inputs()
    centers = join_centers(selected, counts, target_counts, targets)
    require(len(centers) == 570, 'center count differs')
    # Check the saved physical separation using the same explicitly pinned cosmology.
    from astropy.cosmology import FlatLambdaCDM
    cosmology = FlatLambdaCDM(H0=plan['cosmology']['H0_km_s_mpc'], Om0=plan['cosmology']['Omega_m'], Tcmb0=0)
    for r in centers:
        physical = float(cosmology.angular_diameter_distance(r['z_lens']).value)*math.radians(r['offset_from_target_arcmin']/60)
        require(math.isclose(physical, r['offset_from_target_mpc'], rel_tol=1e-9, abs_tol=1e-9), 'physical separation differs')
    geometry = assign_fields(centers, targets)
    lenses = lens_rows(centers)
    require(len(lenses) == 114, 'lens count differs')
    incomplete = [r for r in lenses if r['n_positive'] < 5]
    require(len(incomplete) == 18 and sum(r['zero_raw_count'] for r in centers) == 45, 'count verification not reproduced')
    sensitivity = sensitivities(lenses)
    field_groups = []
    for field in sorted({r['lens_pdr3_field'] for r in lenses}):
        l = [r for r in lenses if r['lens_pdr3_field'] == field]
        c = [r for r in centers if r['lens_pdr3_field'] == field]
        field_groups.append(dict(lens_pdr3_field=field, n_lenses=len(l), n_centers=len(c),
            zero_centers=sum(r['zero_raw_count'] for r in c), n_five_positive=sum(r['n_positive'] == 5 for r in l),
            center_ratio=describe([r['raw_count_ratio'] for r in c])))
    output.mkdir(parents=True, exist_ok=False)
    write_table(output/'per_random_center.tsv', centers)
    write_table(output/'per_lens_support.tsv', lenses)
    write_table(output/'incomplete_lens_groups.tsv', incomplete, list(lenses[0]))
    write_table(output/'diagnostic_sensitivities.tsv', sensitivity)
    write_table(output/'pre_query_exclusions.tsv', [r for r in original_selection if r['selection_status'] != 'selected'])
    plot(centers, lenses, output/'raw_count_support')
    result = dict(schema='astrafek-s16a-wl-raw-random-support-audit-v1',
        status='complete_local_count_support_diagnosis', science_measurement=False,
        network_used=False, cas_query_submitted=False, source_rows_read=False, images_read=False,
        original_selection_modified=False, n_lenses=114, n_random_centers=570,
        incomplete_lens_ids=[r['matched_lens_id'] for r in incomplete],
        by_icl_filter=summarize(centers, lenses), by_pdr3_lens_field=field_groups,
        diagnostic_sensitivities=sensitivity, geometry=geometry,
        original_count_plan_sha256=plan['query_plan_sha256'],
        input_sha256={str(p): sha(p) for p in [*PINNED, RANDOM_CENTERS, TARGET_SELECTION, TARGET_MANIFEST, SELECTION_SUMMARY]},
        script_sha256=sha(Path(__file__)),
        geometry_parser_sha256=sha(Path(__file__).with_name('hsc_pdr3_tract_patch_geometry.py')),
        output_sha256={f.name: sha(f) for f in sorted(output.iterdir()) if f.is_file()},
        limitations=['Counts precede project-specific background-redshift, radial, finite-value and weight checks.',
            'A zero count is catalog unavailability in a cone, not proof of an imaging failure or its cause.',
            'Positive counts and PDR3 field containment do not establish S16A usable aperture area or common radial support.',
            'The 1000-count floor and 0.5-to-2 count-ratio interval are descriptive sensitivities, not validated source/sample cuts.',
            'Mean raw random counts include all five centers, including zeros; no shear is assigned to an empty center.',
            'The catalog-count ratio is not a boost correction or a weak-lensing measurement.',
            'Controls in a lens group and overlapping source cones are dependent; no significance or covariance is estimated.',
            'No S16A/S19A combination or calibrated S16A source-redshift distribution is established.'])
    (output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], output=str(output), summary_sha256=sha(output/'summary.json'),
        by_filter=[{k:r[k] for k in ('icl_filter','n_lenses','n_centers','zero_centers','lenses_with_five_positive','center_ratio','per_lens_mean_ratio')} for r in result['by_icl_filter']],
        by_field=field_groups, sensitivities=sensitivity), indent=2))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    run(p.parse_args().output)
