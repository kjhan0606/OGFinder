#!/usr/bin/env python3
"""Declare a deterministic sky-stratified subset before reading image pixels."""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import measure_s23b_camira_member_outer_profiles as p


def plan(output,n,seed):
    if output.exists(): raise FileExistsError(output)
    paths=(p.TARGET_MANIFEST,p.PAIR_CONTRASTS,p.MEMBER_TABLE,p.MASK_RECORDS)
    hashes={str(x):p.sha256(x) for x in paths}
    expected=(p.EXPECTED_TARGET_MANIFEST_SHA256,p.EXPECTED_PAIR_CONTRASTS_SHA256,
              p.EXPECTED_MEMBER_TABLE_SHA256,p.EXPECTED_MASK_RECORDS_SHA256)
    if list(hashes.values())!=list(expected):raise ValueError('pinned selection inputs changed')
    manifest=json.loads(p.TARGET_MANIFEST.read_text())
    frames={str(r['id']):r for r in manifest['clusters']}
    records={r['id']:r for r in p.read_tsv(p.MASK_RECORDS)}
    cosmo=FlatLambdaCDM(H0=float(manifest['h0']),Om0=float(manifest['omega_m']))
    pool,_=p.select_sample(p.read_tsv(p.PAIR_CONTRASTS),p.read_tsv(p.MEMBER_TABLE),frames,records,cosmo,None,seed)
    rng=np.random.default_rng(seed)
    chosen=[];stats={}
    for band in ('i','i2'):
        grouped=defaultdict(list)
        for item in pool:
            if item['filter']!=band:continue
            f=item['frame']
            cell=(math.floor(((float(f['ra_deg'])-84.056247280625)%360)/5),
                  math.floor((float(f['dec_deg'])+90)/5))
            grouped[cell].append(item)
        if sum(map(len,grouped.values()))<n:raise ValueError('insufficient eligible frames')
        for cell in sorted(grouped):
            grouped[cell].sort(key=lambda x:x['target_id'])
            order=rng.permutation(len(grouped[cell]))
            grouped[cell]=[grouped[cell][int(i)] for i in order]
        remaining={k:list(v) for k,v in grouped.items()}
        selected=[]
        while len(selected)<n:
            cells=sorted(k for k,v in remaining.items() if v)
            for index in rng.permutation(len(cells)):
                cell=cells[int(index)]
                item=remaining[cell].pop()
                f=item['frame']
                selected.append(dict(target_id=item['target_id'],filter=band,
                    ra_deg=float(f['ra_deg']),dec_deg=float(f['dec_deg']),sky_cell=list(cell)))
                if len(selected)==n:break
        chosen.extend(selected)
        stats[band]=dict(eligible=sum(map(len,grouped.values())),eligible_cells=len(grouped),
                         selected=n,selected_cells=len({tuple(x['sky_cell']) for x in selected}))
    payload=dict(schema='astrafek-s23b-member-mask-spatial-sample-v1',status='prepared_local_only',
        seed=seed,cell_width_deg=5,ra_origin_deg=84.056247280625,selection='random order within occupied cells, round-robin random cell order',
        input_sha256=hashes,selection_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sample=sorted(chosen,key=lambda r:(r['filter'],r['target_id'])),stats=stats,
        images_read=False,network_used=False,science_measurement=False)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:json.dump(payload,f,indent=2,allow_nan=False);f.write('\n')
    return payload

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--n-per-filter',type=int,default=48)
    parser.add_argument('--seed',type=int,default=20261117)
    args=parser.parse_args()
    if args.n_per_filter<1:parser.error('positive sample size required')
    print(json.dumps(plan(args.output.resolve(),args.n_per_filter,args.seed)['stats'],indent=2))
