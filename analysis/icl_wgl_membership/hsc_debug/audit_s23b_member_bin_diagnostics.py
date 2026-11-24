#!/usr/bin/env python3
"""Validate S23B per-bin diagnostics and compare a rerun with its saved reference.

No network or image reads. Threshold summaries are conditional descriptive
support counts; no completeness correction or science measurement is made.
"""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import tempfile

from analyze_s23b_member_outer_icl_relation import (
    MEMBER_BINS, MEMBER_DIAGNOSTICS, EXPECTED_HASHES, sha256, read_tsv, write_tsv,
)


def key(row):
    return tuple(row[f] for f in ('target_id', 'filter', 'member_index', 'radial_bin_kpc'))


def audit(measurement, output):
    if output.exists():
        raise FileExistsError(f'refusing to overwrite {output}')
    summary = json.loads((measurement/'summary.json').read_text())
    if summary.get('schema') != 'astrafek-s23b-camira-member-outer-profile-stack-v2':
        raise ValueError('expected diagnostic measurement schema v2')
    if summary.get('science_measurement') is not False or summary.get('network_used') is not False:
        raise ValueError('unexpected science/network status')
    for name, expected in summary['output_sha256'].items():
        if Path(name).name != name or sha256(measurement/name) != expected:
            raise ValueError('measurement output hash mismatch')
    rows = read_tsv(measurement/'member_profile_bin_diagnostics.tsv')
    bins = read_tsv(measurement/'member_profile_bins.tsv')
    diags = read_tsv(measurement/'member_profile_diagnostics.tsv')
    selected = read_tsv(measurement/'selected_target_frames.tsv')
    ids = {(r['target_id'], r['filter']) for r in selected}
    diagnostic_keys = {key(r) for r in rows}
    if len(diagnostic_keys) != len(rows) or len(rows) != len(diags)*7:
        raise ValueError('expected exactly seven unique diagnostics per candidate')
    bin_map = {key(r):r for r in bins}
    if len(bin_map) != len(bins):
        raise ValueError('duplicate retained bin')
    dmap = {(r['target_id'],r['filter'],r['member_index']):r for r in diags}
    if len(dmap) != len(diags):
        raise ValueError('duplicate candidate')
    radial = {'5-10','10-15','15-20','20-30','30-40','40-50','50-75'}
    from_diags = {k+(b,) for k in dmap for b in radial}
    if diagnostic_keys != from_diags:
        raise ValueError('candidate/bin diagnostic coverage changed')
    cfg = summary['measurement']
    tcut, ccut, minimum = cfg['coverage_floor'], cfg['control_coverage_floor'], cfg['minimum_valid_controls_per_bin']
    retained_keys, causes = set(), Counter()
    threshold = Counter()
    threshold_clusters = {}
    max_reconstruction_difference = 0.
    for r in rows:
        k = key(r)
        if k[:2] not in ids or r['source_row'] != dmap[k[:3]]['source_row']:
            raise ValueError('diagnostic candidate/frame join mismatch')
        matched = r['sep_matched'] == 'True'
        controls = json.loads(r['control_bin_measurements_json'])
        if matched and len(controls) != int(r['n_control_profiles']):
            raise ValueError('raw control count mismatch')
        available = r['target_profile_available'] == 'True'
        coverage = float(r['target_coverage']) if available else None
        passing = [c for c in controls if c['coverage'] >= ccut]
        expected = bool(matched and available and coverage >= tcut and len(passing) >= minimum)
        if expected != (r['retained'] == 'True'):
            raise ValueError('saved diagnostic decision fails reconstruction')
        if expected:
            retained_keys.add(k)
            actual = float(r['target_sb_njy_arcsec2']) - sum(c['sb_njy_arcsec2'] for c in passing)/len(passing)
            delta = abs(actual-float(bin_map[k]['delta_sb_njy_arcsec2']))
            max_reconstruction_difference = max(max_reconstruction_difference, delta)
            if delta > 1e-9:
                raise ValueError('saved flux differs from raw-control reconstruction')
        causes[(r['filter'],r['radial_bin_kpc'],r['failure_reason'])] += 1
        for floor in (.4, .6, .8):
            tk = (r['filter'],r['radial_bin_kpc'],floor)
            threshold.setdefault(tk, 0)
            threshold_clusters.setdefault(tk, set())
            if matched and available and coverage >= floor and sum(c['coverage'] >= floor for c in controls) >= minimum:
                threshold[tk] += 1
                threshold_clusters[tk].add(r['target_id'])
    if retained_keys != set(bin_map):
        raise ValueError('retained diagnostics and output profile keys differ')
    if sha256(MEMBER_BINS) != EXPECTED_HASHES[str(MEMBER_BINS)] or sha256(MEMBER_DIAGNOSTICS) != EXPECTED_HASHES[str(MEMBER_DIAGNOSTICS)]:
        raise ValueError('pinned historical reference changed')
    # A historical comparison is valid only for the original measurement cuts.
    if (tcut,ccut,minimum) != (.6,.6,3):
        raise ValueError('historical comparison requires original .6/.6/3 cuts')
    old = {key(r):r for r in read_tsv(MEMBER_BINS) if (r['target_id'],r['filter']) in ids}
    if old != bin_map:
        raise ValueError('rerun profile rows differ from the historical selected subset')
    old_diags = {(r['target_id'],r['filter'],r['member_index']):r for r in read_tsv(MEMBER_DIAGNOSTICS) if (r['target_id'],r['filter']) in ids}
    if old_diags != dmap:
        raise ValueError('rerun candidate diagnostics differ from historical subset')
    cause_rows = [dict(filter=k[0],radial_bin_kpc=k[1],failure_reason=k[2],n_candidate_bins=n)
                  for k,n in sorted(causes.items())]
    support_rows = [dict(filter=k[0],radial_bin_kpc=k[1],coverage_floor=k[2],
                        min_controls=minimum,n_retained_candidates=n,
                        n_retained_clusters=len(threshold_clusters[k]))
                    for k,n in sorted(threshold.items())]
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name+'.partial.',dir=output.parent))
    write_tsv(tmp/'failure_combinations.tsv',cause_rows,list(cause_rows[0]))
    write_tsv(tmp/'coverage_threshold_support.tsv',support_rows,list(support_rows[0]))
    result = dict(schema='astrafek-s23b-member-bin-diagnostic-audit-v1',
        status='complete_local_reconstruction_and_reference_audit',
        measurement=str(measurement),measurement_summary_sha256=sha256(measurement/'summary.json'),
        script_sha256=sha256(Path(__file__)),
        reference_sha256={str(p):sha256(p) for p in (MEMBER_BINS,MEMBER_DIAGNOSTICS)},
        n_frames=len(selected),n_candidates=len(diags),n_candidate_bin_diagnostics=len(rows),
        n_retained_bins=len(bins),historical_profile_rows_identical=True,
        historical_candidate_diagnostics_identical=True,
        max_flux_reconstruction_difference=max_reconstruction_difference,
        output_sha256={p.name:sha256(p) for p in tmp.iterdir()},
        science_measurement=False,network_used=False,images_read=False,
        limitations=['Threshold support counts are descriptive and conditional on upstream SEP selection.',
                    'Different thresholds select different candidates; no flux trend or association inference is made.',
                    'Pilot results cannot be generalized to the full 979-frame sample.'])
    (tmp/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    os.rename(tmp,output)
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--measurement',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(audit(args.measurement.resolve(),args.output.resolve()),indent=2))
