#!/usr/bin/env python3
"""Freeze 40 area-qualified engineering targets plus 10 mask stress targets."""
import argparse
from collections import Counter, defaultdict
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import random


def strata(row):
    return (int(math.floor((float(row['z_cl'])-.3+1e-8)/.05)),
            next((i for i,v in enumerate((20,30,50)) if float(row['N_mem'])<v),3))


def select(rows, n, seed):
    if len(rows)<n: raise ValueError('Insufficient eligible rows')
    groups=defaultdict(list)
    for row in rows: groups[strata(row)].append(row)
    quotas={k:0 for k in groups}
    if len(groups)<=n: quotas={k:1 for k in groups}
    while sum(quotas.values())<n:
        available=[k for k in groups if quotas[k]<len(groups[k])]
        k=max(sorted(available),key=lambda k:n*len(groups[k])/len(rows)-quotas[k])
        quotas[k]+=1
    rng=random.Random(seed);chosen=[]
    for k in sorted(groups):
        group=sorted(groups[k],key=lambda r:str(r['id']))
        chosen.extend(rng.sample(group,quotas[k]))
    return chosen, [dict(z_bin=k[0],richness_bin=k[1],eligible=len(groups[k]),selected=quotas[k],
                         conditional_inclusion_probability=quotas[k]/len(groups[k])) for k in sorted(groups)]


def write(path, value):
    with path.open('x') as f: json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    args=p.parse_args()
    if not platform.node().lower().startswith('lageunha'): p.error('Run on lageunha CPU')
    inputs={k:args.sample/v for k,v in {
        'selection':'selected_clusters.tsv','manifest':'download_01/download_manifest.json',
        'quality':'quality_2026-10-13/quality_frames.json',
        'requests':'ancillary_2026-10-13_01/hsc_catalogue_requests.json'}.items()}
    rows=list(csv.DictReader(inputs['selection'].open(),delimiter='\t'))
    quality={str(r['id']):r for r in json.loads(inputs['quality'].read_text())}
    if len(quality)!=len(rows) or set(quality)!={str(r['id']) for r in rows}: raise ValueError('Input IDs differ')
    good=[];stress=[]
    for row in rows:
        q=quality[row['id']];coverage=q['coverage_by_25kpc_annulus_0_500']
        if len(coverage)!=20 or not all(math.isfinite(v) and 0<=v<=1 for v in coverage):
            raise ValueError('Invalid radial coverage')
        if min(coverage)>=.5: good.append(row)
        elif q['center_10kpc_bright_object_mask_fraction']>.25: stress.append(row)
    main_rows,main_strata=select(good,40,20261014)
    stress_rows,stress_strata=select(stress,10,20261015)
    membership={r['id']:group for group,rs in [('area_qualified',main_rows),('mask_stress',stress_rows)] for r in rs}
    assert len(membership)==50
    doc=json.loads(inputs['manifest'].read_text());out=copy.deepcopy(doc)
    out['clusters']=[]
    for frame in doc['clusters']:
        if str(frame['id']) not in membership: continue
        frame=copy.deepcopy(frame);frame['preparation_cohort']=membership[str(frame['id'])]
        frame['image']=str((args.sample/'download_01/images'/frame['image']).resolve())
        for key in ('bad_mask','bright_star_mask','variance'):
            if frame.get(key) and frame[key].get('file'):
                frame[key]['file']=str((args.sample/'download_01/images'/frame[key]['file']).resolve())
        if not Path(frame['image']).is_file(): raise ValueError('Missing science image')
        out['clusters'].append(frame)
    out.update(science_measurement=False,selection='40 stratified area-qualified + 10 disjoint mask-stress engineering targets; not representative of all 800')
    out['preparation_gates']={k:'pending' for k in ('pdr3_extinction_crosscheck','matched_random_images','psf_wings','source_mask_injection_recovery','sky_transfer')}
    args.work.mkdir(parents=True,exist_ok=False)
    write(args.work/'validation_manifest.json',out)
    with (args.work/'targets.tsv').open('x',newline='') as f:
        fields=['id','cohort','z_cl','z_spec_bcg','N_mem','minimum_radial_coverage','center_bright_mask_fraction','header_filter']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for row in main_rows+stress_rows:
            q=quality[row['id']]
            w.writerow(dict(id=row['id'],cohort=membership[row['id']],z_cl=row['z_cl'],z_spec_bcg=row['z_spec_bcg'],N_mem=row['N_mem'],
                minimum_radial_coverage=min(q['coverage_by_25kpc_annulus_0_500']),
                center_bright_mask_fraction=q['center_10kpc_bright_object_mask_fraction'],header_filter=q['header_filter']))
    requests=json.loads(inputs['requests'].read_text())
    write(args.work/'hsc_catalogue_requests.json',[r for r in requests if r['cluster_id'] in membership])
    (args.work/'patch_qa.sql').write_text('SELECT * FROM pdr3_wide.patch_qa;\n')
    write(args.work/'selection_audit.json',dict(parent=len(rows),eligible=len(good),mask_stress_pool=len(stress),
        main_count=40,stress_count=10,main_seed=20261014,stress_seed=20261015,
        main_strata=main_strata,stress_strata=stress_strata,
        source_sha256={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in inputs.items()},
        representative_of_parent=False,science_measurement=False,
        notes=['z_cl controls stratification; z_spec_bcg controls physical scale',
               'No brightness/profile outcome used in target choice',
               'Conditional probabilities apply only to each engineering eligibility pool',
               'Raw images, extinction acceptance and original sample selection unchanged']))
    print(f'Prepared 40/{len(good)} area-qualified and 10/{len(stress)} mask-stress targets')


if __name__=='__main__': main()
