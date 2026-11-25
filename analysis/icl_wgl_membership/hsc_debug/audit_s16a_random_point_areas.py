#!/usr/bin/env python3
"""S16A random-point areas and the bright-object null fraction (2026-11-25 follow-up of the three approved queries).

A. f_rand (bright-object flag plan null expectation) from the thinned random candidate rows: object_id-unique,
   inner-deduplicated, i-band processed non-edge points inside the union of the non-empty 3-Mpc source discs;
   flagged fraction of iflags_pixel_bright_object_center/_any, pooled and source-weighted, and the pre-registered
   decision rule re-evaluated with it.
B. Inner areas (0-0.1, 0.1-0.2, 0.2-0.4 proper Mpc) of the 66 pilot centres from full-density random points:
   area = geometric annulus x (deduplicated points passing the hole rule)/(deduplicated points present). Hole rule
   (bright-object flag outcome): random iflags_pixel_bright_object_center, any of the five band edge flags, or
   i-unprocessed (NULL i flags, icountinputs = 0). Pre-registered: NULL rows count as holes; a sensitivity treats
   NULL rows lying under full inner density as overlap duplicates. Centre-point status is reported only.
C. 0-0.4 Mpc boosts (expected-count stack over the supported lens groups) with label areas (the recorded 1-arcsec,
   30-arcsec convention, own-source rule; reproduced exactly) and with random-point areas.

Local verified inputs only; descriptive; no selection change, correction or science promotion.
"""
import argparse
import copy
import csv
import glob
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_afterburner_mask_effects import tangent_arcsec
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
INNER = DATA/'s16a_random_inner_rows_plan_2026-11-25_01'
CAND = DATA/'s16a_random_candidate_rows_plan_2026-11-25_02'
FLAGS = DATA/'s16a_bright_object_flag_rows_plan_2026-11-25_03'
PLANS = {INNER: 'e8b05d6e8ed2235930b797a0e21332b2915f9df116d7918186a194d31f0d5ff6',
         CAND: 'cdd04062b79415328e7bbdae08a938a5224bf19daa2717d2d9665f81ed21b356',
         FLAGS: 'e73c0f9027c5f75162516430703ee70b7eb6e87d4d454ca5f33b6ea5c2c7062b'}
PRIOR = Path('/home/kjhan/BACKUP/Astrafek/s16a_boost_mask_field_tests_20261125_03')
PRIOR_SUMMARY_SHA = '55f94a3d'  # prefix check; full hash recorded in the output
DENSITY = 100.0  # documented full-density random points per arcmin^2
STEP = 1.0
FOOTPRINT = 30.0
INNER_BINS = (0, 1, 2)
CENTRE_ARCSEC = 6.0
DUP_PROBE_ARCSEC = 15.0
DUP_MIN_FRACTION = 0.5  # NULL row is an overlap duplicate if inner density within 15" >= half nominal
DEFICIT_SIGMA = 3.0
BANDS = 'grizy'
EDGES = [f'{b}flags_pixel_edge' for b in BANDS]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def truthy(v):
    return v in ('True', 'true', 't', '1')


def load_rows(plan_dir):
    summary = json.loads((plan_dir/'results'/'download_summary.json').read_text())
    if summary['approved_scope_plan_sha256'] != PLANS[plan_dir] or summary.get('credentials_saved') is not False:
        raise ValueError('receipt differs '+str(plan_dir))
    rows = []
    for name, digest in sorted(summary['output_sha256'].items()):
        path = plan_dir/'results'/name
        if sha(path) != digest:
            raise ValueError('csv hash differs '+name)
        with path.open() as stream:
            reader = csv.reader(stream)
            head = [h.lstrip('# ').strip() for h in next(reader)]
            rows.extend(dict(zip(head, r)) for r in reader)
    return summary, rows


def inner_ok(r):
    return truthy(r['idetect_is_patch_inner']) and truthy(r['idetect_is_tract_inner'])


def i_null(r):
    return r['icountinputs'] == '' or r['idetect_is_patch_inner'] == ''


def hole(r):
    if i_null(r) or r['icountinputs'] == '0':
        return True
    return truthy(r['iflags_pixel_bright_object_center']) or any(truthy(r[e]) for e in EDGES)


def unit(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.c_[np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)]


def sep_arcsec(u, v):
    return np.degrees(2*np.arcsin(np.clip(np.linalg.norm(u-v, axis=-1)/2, 0, 1)))*3600


def part_f_rand(centers, data, flag_summary):
    _, rows = load_rows(CAND)
    seen, pts = set(), []
    for r in rows:
        if r['object_id'] in seen or i_null(r) or not inner_ok(r) or truthy(r['iflags_pixel_edge']) or r['icountinputs'] == '0':
            continue
        seen.add(r['object_id'])
        pts.append(r)
    u = unit([float(r['ra']) for r in pts], [float(r['dec']) for r in pts])
    centre = np.array([truthy(r['iflags_pixel_bright_object_center']) for r in pts])
    anyf = np.array([truthy(r['iflags_pixel_bright_object_any']) for r in pts])
    lens_pos = {c['matched_lens_id']: unit([float(c['ra_deg'])], [float(c['dec_deg'])])[0] for c in centers if c['role'] == 'lens'}
    union = np.zeros(len(pts), bool)
    discs = []
    for c in centers:
        e = data[c['center_id']]
        if not e['n']:
            continue
        cu = unit([float(c['ra_deg'])], [float(c['dec_deg'])])[0]
        inside = sep_arcsec(u, cu) <= e['radius']
        union |= inside
        reach = sep_arcsec(cu[None], lens_pos[c['matched_lens_id']][None])[0]+e['radius']
        discs.append(dict(center_id=c['center_id'], role=c['role'], lens_id=c['matched_lens_id'], n_sources=e['n'], n_points=int(inside.sum()),
                          f_center=float(centre[inside].mean()) if inside.any() else None,
                          f_any=float(anyf[inside].mean()) if inside.any() else None,
                          disc_inside_query_cone=bool(reach <= 90*60)))
    n = int(union.sum())
    f_c, f_a = float(centre[union].mean()), float(anyf[union].mean())
    w = [(d['n_sources'], d['f_center'], d['f_any']) for d in discs if d['f_center'] is not None]
    sw_c = sum(a*b for a, b, _ in w)/sum(a for a, _, _ in w)
    sw_a = sum(a*c for a, _, c in w)/sum(a for a, _, _ in w)
    plan = json.loads((FLAGS/'plan.json').read_text())['decision_rule']
    true = flag_summary['true_by_column']
    N = flag_summary['n_unique_objects']
    T_c, T_a = true.get('meas_iflags_pixel_bright_object_center', 0), true.get('meas_iflags_pixel_bright_object_any', 0)
    amax, frac, mint = plan['applied_max_true'], plan['not_applied_random_fraction'], plan['min_testable_expected']

    def centre_rule(f):
        if T_c <= amax:
            return 'APPLIED'
        if frac*N*f < mint:
            return 'UNTESTABLE'
        return 'NOT APPLIED' if T_c/N >= frac*f else 'AMBIGUOUS'
    return dict(n_disc_points=n, n_discs=len(discs), discs_partly_outside_cone=sum(not d['disc_inside_query_cone'] for d in discs),
                f_rand_center=f_c, f_rand_center_binomial_se=math.sqrt(f_c*(1-f_c)/n), f_rand_any=f_a,
                f_rand_center_source_weighted=sw_c, f_rand_any_source_weighted=sw_a,
                N_sources=N, T_meas_center=T_c, T_meas_any=T_a, source_any_fraction=T_a/N,
                expected_center_if_not_applied=f_c*N, testability_expected=frac*N*f_c, min_testable_expected=mint,
                rule_center_outcome=centre_rule(f_c), rule_center_outcome_source_weighted=centre_rule(sw_c),
                rule_any_outcome='APPLIED' if T_a <= amax else 'NOT APPLIED',
                any_source_over_random=(T_a/N)/f_a if f_a > 0 else None), discs


def part_areas(centers, data):
    summary, rows = load_rows(INNER)
    by_c = {}
    for r in rows:
        by_c.setdefault(r['center_id'], []).append(r)
    # duplicate probe for NULL rows: inner, i-processed points pooled over all centres (object_id unique)
    pool, seen = [], set()
    for r in rows:
        if not i_null(r) and inner_ok(r) and r['object_id'] not in seen:
            seen.add(r['object_id']); pool.append((float(r['ra']), float(r['dec'])))
    pu = unit(*zip(*pool))
    tree = cKDTree(pu)
    chord = 2*math.sin(math.radians(DUP_PROBE_ARCSEC/3600)/2)
    nominal = DENSITY*math.pi*(DUP_PROBE_ARCSEC/60)**2
    out, centre_status, null_info = {}, [], []
    for c in centers:
        cid, e = c['center_id'], data[c['center_id']]
        rs = by_c.get(cid, [])
        x, y = tangent_arcsec([float(r['ra']) for r in rs], [float(r['dec']) for r in rs], float(c['ra_deg']), float(c['dec_deg']))
        rad = np.hypot(x, y)
        mpc = rad*3.0/e['radius']
        nulls = np.array([i_null(r) for r in rs], bool)
        dup = np.zeros(len(rs), bool)
        if nulls.any():
            nu = unit([float(r['ra']) for r, k in zip(rs, nulls) if k], [float(r['dec']) for r, k in zip(rs, nulls) if k])
            local = np.array([len(v) for v in tree.query_ball_point(nu, chord)])
            dup[np.flatnonzero(nulls)] = local >= DUP_MIN_FRACTION*nominal
            null_info.append(dict(center_id=cid, null_rows=int(nulls.sum()), null_rows_under_full_inner_density=int(dup.sum()),
                                  median_local_inner_over_nominal=float(np.median(local)/nominal)))
        keep = np.array([inner_ok(r) for r in rs], bool) | nulls
        holes = np.array([hole(r) for r in rs], bool)
        res = {}
        for variant, present in (('preregistered', keep), ('null_duplicates', keep & ~dup)):
            vals = []
            for k in INNER_BINS:
                lo, hi = RADIAL_EDGES[k], RADIAL_EDGES[k+1]
                sel = present & (mpc >= lo) & (mpc < hi)
                a_geo = math.pi*(hi**2-lo**2)*(e['radius']/3.0/60)**2
                n_pres, n_pass = int(sel.sum()), int((sel & ~holes).sum())
                expd = a_geo*DENSITY
                vals.append(dict(annulus=k, a_geo=a_geo, n_expected=expd, n_present=n_pres, n_pass=n_pass,
                                 n_bright=int(sum(truthy(rs[i]['iflags_pixel_bright_object_center']) for i in np.flatnonzero(sel))),
                                 n_edge=int(sum(any(truthy(rs[i][ed]) for ed in EDGES) for i in np.flatnonzero(sel))),
                                 n_unprocessed=int(sum(i_null(rs[i]) or rs[i]['icountinputs'] == '0' for i in np.flatnonzero(sel))),
                                 deficit_sigma=(expd-n_pres)/math.sqrt(expd),
                                 pass_fraction=n_pass/n_pres if n_pres else None,
                                 area=a_geo*n_pass/n_pres if n_pres else 0.0,
                                 area_binomial_se=a_geo*math.sqrt(n_pass*(n_pres-n_pass)/n_pres**3) if n_pres else None))
            res[variant] = vals
        out[cid] = res
        near = keep & ~dup & (rad <= CENTRE_ARCSEC)
        nh = int((near & holes).sum())
        centre_status.append(dict(center_id=cid, role=c['role'], lens_id=c['matched_lens_id'], n_sources=e['n'], points_within_6arcsec=int(near.sum()),
                                  holes_within_6arcsec=nh, status=('no_points' if not near.any() else 'masked' if 2*nh > near.sum() else 'unmasked')))
    dens = [v['n_present']/v['a_geo'] for r in out.values() for v in r['null_duplicates'] if v['n_present'] > 0 and v['annulus'] == 2]
    return summary, out, centre_status, null_info, float(np.median(dens))


def label_parts(entry):
    """Recorded 1-arcsec, 30-arcsec label parts restricted to the inner lattice (same lattice as the full-disc grid)."""
    xy = bt.grid(entry['radius']*RADIAL_EDGES[3]/3.0, STEP)
    bins = bt.annulus(np.hypot(xy[:, 0], xy[:, 1])*3.0/entry['radius'])
    if entry['n']:
        dist, index = entry['tree'].query(xy, workers=bt.WORKERS)
    else:
        dist, index = np.full(len(xy), np.inf), np.zeros(len(xy), int)
    return bt.parts(entry, dict(xy=xy, bins=bins, dist=dist, index=index, step=STEP), FOOTPRINT, entry.get('missing'))


def with_area(p, vals):
    q = copy.deepcopy(p)
    for v in vals:
        q['a_unm'][v['annulus']] = v['area']
        q['a_fp'][v['annulus']] = v['a_geo']*min(1.0, v['n_present']/v['n_expected']) if v['n_expected'] else 0.0
    return q


def run(output):
    if output.exists():
        raise FileExistsError(output)
    prior_sha = sha(PRIOR/'summary.json')
    if not prior_sha.startswith(PRIOR_SUMMARY_SHA):
        raise ValueError('prior boost summary changed')
    prior = json.loads((PRIOR/'summary.json').read_text())
    if sha(PRIOR/'normalisation_summary.tsv') != prior['output_sha256']['normalisation_summary.tsv']:
        raise ValueError('prior normalisation table changed')
    centers, data, _, _ = bt.load_centres()
    flag_summary = json.loads((FLAGS/'results'/'download_summary.json').read_text())
    if flag_summary['approved_scope_plan_sha256'] != PLANS[FLAGS]:
        raise ValueError('flag receipt differs')
    f_rand, discs = part_f_rand(centers, data, flag_summary)
    inner_summary, areas, centre_status, null_info, median_density = part_areas(centers, data)

    labels = {c['center_id']: label_parts(data[c['center_id']]) for c in centers}
    lenses = {c['matched_lens_id']: c['center_id'] for c in centers if c['role'] == 'lens'}
    rand = {}
    for c in centers:
        if c['role'] == 'random':
            rand.setdefault(c['matched_lens_id'], []).append(c['center_id'])
    nonempty = {g: sorted(r for r in rs if data[r]['n']) for g, rs in rand.items()}
    supported = [g for g in lenses if sum(labels[r]['a_fp'].sum() for r in rand[g]) > 0]
    supported = [g for g in sorted(lenses) if g in supported]

    # Exact reproduction of the recorded expected-count, area-state lens stacks (own-source rule) in 0-0.4 Mpc.
    recorded = {}
    with (PRIOR/'normalisation_summary.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['definition'] == 'expected' and row['state'] in ('pre', 'area') and int(row['annulus']) in INNER_BINS:
                recorded[(row['state'], row['weighting'], int(row['annulus']))] = float(row['lens_stack_supported_minus_1'])
    prior_groups = sorted(supported)

    def items(source, groups, randoms=nonempty):
        return [(source[lenses[g]], [source[r] for r in randoms[g]]) for g in groups]
    max_diff = 0.0
    for (state, kind, k), value in recorded.items():
        mine = bt.stack('expected', items(labels, prior_groups), k, kind, state, (), True)[0]-1
        max_diff = max(max_diff, abs(mine-value))
    if max_diff > 1e-12:
        raise ValueError(f'recorded inner expected-count stacks not reproduced ({max_diff})')

    variants = {}
    for name in ('preregistered', 'null_duplicates'):
        variants[name] = {cid: with_area(labels[cid], areas[cid][name]) for cid in labels}
    all_rand = {g: sorted(rs) for g, rs in rand.items()}
    boost_rows = []
    configs = [('label_area_own_rule (recorded)', labels, True, nonempty, 'area'),
               ('label_area_no_own_rule', labels, False, nonempty, 'area'),
               ('label_pre_mask', labels, True, nonempty, 'pre'),
               ('random_point_area (preregistered NULL rule)', variants['preregistered'], False, nonempty, 'area'),
               ('random_point_area (NULL overlap duplicates)', variants['null_duplicates'], False, nonempty, 'area'),
               ('random_point_area + empty randoms', variants['null_duplicates'], False, all_rand, 'area')]
    for label, source, own, randoms, state in configs:
        for kind in ('n', 'w'):
            for k in INNER_BINS:
                value, used, se, _ = bt.stack_jackknife('expected', items(source, prior_groups, randoms), k, kind, state, (), None, own)
                boost_rows.append(dict(area=label, state=state, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                       lens_groups=len(used), boost=value, group_jackknife_se=se))

    centre_rows = []
    for c in centers:
        cid = c['center_id']
        for v in areas[cid]['null_duplicates']:
            k = v['annulus']
            pre = areas[cid]['preregistered'][k]
            lp = labels[cid]
            centre_rows.append(dict(center_id=cid, role=c['role'], lens_id=c['matched_lens_id'], n_sources_disc=data[cid]['n'], annulus=k,
                                    r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], **{kk: v[kk] for kk in v if kk != 'annulus'},
                                    area_preregistered=pre['area'], label_area_fp=lp['a_fp'][k], label_area_unmasked=lp['a_unm'][k],
                                    n_pre=lp['n_pre'][k], n_post=lp['n_post'][k],
                                    random_point_over_label_unmasked=v['area']/lp['a_unm'][k] if lp['a_unm'][k] > 0 else None,
                                    missing_patch_flag=v['deficit_sigma'] > DEFICIT_SIGMA))
    output.mkdir(parents=True)
    tables = {'f_rand_discs.tsv': discs, 'centre_annulus_areas.tsv': centre_rows, 'centre_point_status.tsv': centre_status,
              'inner_boosts.tsv': boost_rows, 'null_rows.tsv': null_info}
    for name, rows in tables.items():
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    focus = [r for r in centre_rows if r['center_id'] in ('RC-6065306c1c2aaca36993', 'RC-1d54db20ff04be924a4f')]
    summary = dict(
        schema='astrafek-s16a-random-point-areas-v1', science_measurement=False, selection_change=False,
        inputs=dict(plans=dict((str(k), v) for k, v in PLANS.items()), prior_summary_sha256=prior_sha,
                    inner_csv_sha256=inner_summary['output_sha256']),
        constants=dict(density=DENSITY, step_arcsec=STEP, footprint_arcsec=FOOTPRINT, centre_arcsec=CENTRE_ARCSEC,
                       duplicate_probe_arcsec=DUP_PROBE_ARCSEC, duplicate_min_fraction=DUP_MIN_FRACTION, deficit_sigma=DEFICIT_SIGMA,
                       hole_rule='i bright_object_center OR any grizy edge OR i-unprocessed (NULL, icountinputs=0)'),
        f_rand=f_rand, recorded_inner_stack_reproduction_max_abs_diff=max_diff, supported_groups=prior_groups,
        median_dedup_density_0p2_0p4=median_density, null_rows=null_info,
        missing_patch_flags=sorted({r['center_id'] for r in centre_rows if r['missing_patch_flag']}),
        centre_status_counts={s: sum(x['status'] == s for x in centre_status) for s in ('masked', 'unmasked', 'no_points')},
        masked_centres=[x for x in centre_status if x['status'] != 'unmasked'],
        focus_centres=focus,
        limitations=['Random points carry no wl_fulldepth_fullcolor flag: FDFC cuts are not in the random-point areas.',
                     'Binomial area noise (area_binomial_se) is not propagated into the group jackknife.',
                     'The hole rule follows the bright-object flag plan outcome; sources are masked only by S16A afterburner absence.'],
        runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in tables})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_random_point_areas_20261125_01'))
    s = run(parser.parse_args().output)
    print(json.dumps({k: s[k] for k in ('f_rand', 'recorded_inner_stack_reproduction_max_abs_diff', 'null_rows', 'missing_patch_flags',
                                        'centre_status_counts', 'median_dedup_density_0p2_0p4')}, indent=1, default=float))
