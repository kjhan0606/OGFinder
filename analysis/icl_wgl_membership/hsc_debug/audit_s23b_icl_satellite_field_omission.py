#!/usr/bin/env python3
"""Omit ICL target/control endpoints by official HSC footprint field."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from analyze_s23b_icl_satellite_relation import (
    block_bootstrap_slope, connected_pair_blocks, model_matrix, read_tsv, sha256, write_tsv,
)
from analyze_s23b_wl_field_component_influence import DATA, FIELD, FIELD_SHA
from hsc_pdr3_tract_patch_geometry import load_index

RELATION = DATA/'s23b_icl_satellite_relation_2026-11-07_03'
RELATION_SHA = 'fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39'
FLOORS = (0., .1, .2)
N_BOOTSTRAP = 5000
ORIGIN = 84.056247280625


def fit(rows, seed):
    estimate = dict(n_pairs=len(rows), connected_blocks=0, slope=None, ci_low=None, ci_high=None,
                    status='not_estimable', reason='')
    if len(rows) < 12:
        estimate['reason'] = 'fewer than 12 pairs'
        return estimate
    controls = {r['target_id']:dict(ra_deg=r['control_ra_deg'],dec_deg=r['control_dec_deg']) for r in rows}
    labels, blocks = connected_pair_blocks(rows,controls,5.,0.,ORIGIN)
    estimate['connected_blocks'] = len(blocks)
    try:
        design,y,_ = model_matrix(rows)
    except ValueError as error:
        estimate['reason'] = str(error)
        return estimate
    estimate.update(status='point_only',slope=float(np.linalg.lstsq(design,y,rcond=None)[0][1]))
    try:
        slope, low, high, _, n = block_bootstrap_slope(design,y,labels,N_BOOTSTRAP,seed)
    except ValueError as error:
        estimate['reason'] = str(error)
    else:
        estimate.update(status='complete',slope=slope,ci_low=low,ci_high=high,valid_bootstrap=n)
    return estimate


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_file = RELATION/'summary.json'
    field_file = FIELD/'summary.json'
    if sha256(summary_file)!=RELATION_SHA or sha256(field_file)!=FIELD_SHA:
        raise ValueError('changed ICL or field provenance')
    prior, geometry = (json.loads(path.read_text()) for path in (summary_file,field_file))
    inputs = {}
    for name in ('paired_icl_satellite_profiles.tsv','connected_block_slope_sensitivity.tsv','common_support_floor_sensitivity.tsv'):
        path = RELATION/name
        inputs[name] = sha256(path)
        if inputs[name]!=prior['output_sha256'][name]:
            raise ValueError(f'changed relation table: {name}')
    rows = [row for row in read_tsv(RELATION/'paired_icl_satellite_profiles.tsv') if row['icl_relative_contrast_njy_arcsec2']]
    if Counter(r['filter'] for r in rows) != {'i':232,'i2':756} or len({r['target_id'] for r in rows})!=988:
        raise ValueError('supported ICL membership differs')
    directory = Path(geometry['inputs']['official_geometry_directory'])
    for name,meta in geometry['official_geometry_summary']['files'].items():
        if sha256(directory/name)!=meta['sha256']:
            raise ValueError(f'changed footprint geometry: {name}')
    index,_ = load_index(directory,level='tract')
    for member, keys in (('target',('ra_deg','dec_deg')),('control',('control_ra_deg','control_dec_deg'))):
        assigned = index.assign(*[[float(r[k]) for r in rows] for k in keys],level='tract')
        if any(f is None for f in assigned.field) or np.any(assigned.cross_field_ambiguous) or np.any(assigned.tract<0):
            raise ValueError(f'{member} endpoint has no unique official survey field')
        for row,field in zip(rows,assigned.field,strict=True):
            row[member+'_field'] = str(field)
    fields = sorted({r[k+'_field'] for r in rows for k in ('target','control')})
    previous = {}
    for floor, name in ((0.,'connected_block_slope_sensitivity.tsv'),(.1,'common_support_floor_sensitivity.tsv'),(.2,'common_support_floor_sensitivity.tsv')):
        for r in read_tsv(RELATION/name):
            if float(r['block_width_deg'])==5 and float(r['grid_phase_fraction'])==0 and float(r.get('minimum_common_support_fraction_both_annuli',0))==floor:
                previous[r['filter'],floor] = r
    fits = []
    for j,band in enumerate(('i','i2')):
        for k,floor in enumerate(FLOORS):
            selected = [r for r in rows if r['filter']==band and min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=floor]
            original = fit(selected, 2026111801+j*10+k)
            baseline = previous[band,floor]
            if original['n_pairs']!=int(baseline['n_pairs']) or not np.isclose(original['slope'],float(baseline['density_ratio_slope_njy_arcsec2_per_dex']),atol=1e-10,rtol=0):
                raise ValueError('full-sample slope replay differs')
            fits.append(dict(filter=band,support_floor=floor,omitted_field='none',mode='full',**original))
            for i,field in enumerate(fields):
                # Match seeds when the two omission rules select identical rows.
                for mode in ('target_only','either_endpoint'):
                    after = [r for r in selected if r['target_field']!=field and
                             (mode=='target_only' or r['control_field']!=field)]
                    result = fit(after,2026111901+j*100+k*10+i)
                    fits.append(dict(filter=band,support_floor=floor,omitted_field=field,mode=mode,**result))
    output.mkdir(parents=True,exist_ok=False)
    table, ledger = output/'field_slope_sensitivity.tsv', output/'pair_field_membership.tsv'
    write_tsv(table,fits)
    write_tsv(ledger,[{key:r[key] for key in ('target_id','control_id','filter','target_field','control_field')} for r in rows])
    matches = Counter((r['filter'],r['target_field'],r['control_field']) for r in rows)
    summary = dict(science_measurement=False, relation_summary_sha256=RELATION_SHA, field_summary_sha256=FIELD_SHA,
                   layouts='5 degrees, phase 0. Intervals are exploratory connected-block percentile intervals.',
                   field_pair_counts=[dict(filter=b,target_field=t,control_field=c,n_pairs=n) for (b,t,c),n in sorted(matches.items())],
                   cross_field_pairs=sum(n for (_,a,b),n in matches.items() if a!=b),
                   n_fits=len(fits), bootstrap_replicas=N_BOOTSTRAP,
                   output_sha256={p.name:sha256(p) for p in (table,ledger)},runner_sha256=sha256(Path(__file__)),
                   limitation='Field influence of saved ICL profiles. Pairs with an endpoint in the field are removed. No pixel-level masking, selection refit, or survey covariance.')
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
