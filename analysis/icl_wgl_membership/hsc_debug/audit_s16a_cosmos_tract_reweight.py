#!/usr/bin/env python3
"""Compare separately reweighted COSMOS tracts with explicit support accounting."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, sha, write_table
from ogfkit.wlcluster import Cosmology

SUPPORT = INTAKE.parent/'s16a_cosmos_bin_support_2026-11-21_01'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    intake_path, support_path = INTAKE/'summary.json', SUPPORT/'summary.json'
    intake, support = (json.loads(p.read_text()) for p in (intake_path, support_path))
    for folder, report in ((INTAKE, intake), (SUPPORT, support)):
        for filename, expected in report['output_sha256'].items():
            if sha(folder/filename) != expected:
                raise ValueError('Changed input file')
    if support['intake_sha256'] != sha(intake_path):
        raise ValueError('Different COSMOS input lineage')
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        data = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review']=='True']
    if len(data) != intake['eligible_rows']:
        raise ValueError('COSMOS eligibility mismatch')
    photo, reference, shape, som = np.asarray([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in data]).T
    tract=np.asarray([r['tract'] for r in data])
    with (SUPPORT/'lens_photoz_reweight.tsv').open() as stream:
        lenses=list(csv.DictReader(stream,delimiter='\t'))
    with (SUPPORT/'photoz_bin_support.tsv').open() as stream:
        bins=list(csv.DictReader(stream,delimiter='\t'))
    edges=np.asarray(support['edges'])
    model=Cosmology()
    results=[]
    common_results=[]
    for lens in lenses:
        zl=float(lens['lens_z'])
        target=np.asarray([float(r['observed_weight_fraction']) for r in bins if r['lens_id']==lens['lens_id']])
        if len(target)!=len(edges)-1 or not np.isclose(target.sum(),1):
            raise ValueError('Invalid target histogram')
        region_bins={}
        for region in ('9812','9813'):
            keep=(photo>zl+0.2)&(tract==region)
            p,ref,sh,so=photo[keep],reference[keep],shape[keep],som[keep]
            chi=float(model.comoving(zl))
            est=1-chi/model.comoving(p)
            true=np.zeros_like(ref)
            behind=ref>zl
            true[behind]=1-chi/model.comoving(ref[behind])
            weights=sh*so*est**2
            index=np.searchsorted(edges,p,side='right')-1
            mass=np.bincount(index,weights=weights,minlength=len(target))
            num=np.bincount(index,weights=weights*true/est,minlength=len(target))
            region_bins[region]=(mass,num)
            missing=float(target[mass<=0].sum())
            response=None
            if missing==0:
                means=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
                response=float(1/(target@means))
            results.append(dict(lens_id=lens['lens_id'],tract=region,n_objects=int(keep.sum()),
                                unsupported_weight_fraction=missing,
                                baseline_combined=float(lens['reweighted']),tract_multiplier=response,
                                fractional_change=response/float(lens['reweighted'])-1 if response is not None else None))
        common=(region_bins['9812'][0]>0)&(region_bins['9813'][0]>0)
        retained=float(target[common].sum())
        if retained <= 0:
            raise ValueError('No common weighted bins')
        common_target=np.where(common,target/retained,0.)
        multipliers={}
        for region,(mass,num) in region_bins.items():
            means=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
            multipliers[region]=float(1/(common_target@means))
        common_results.append(dict(lens_id=lens['lens_id'],retained_weight_fraction=retained,
                                   multiplier_9812=multipliers['9812'],multiplier_9813=multipliers['9813'],
                                   relative_difference=multipliers['9812']/multipliers['9813']-1))
    output.mkdir(parents=True)
    table=output/'tract_reweight.tsv'
    write_table(table,results)
    common_table=output/'tract_common_support.tsv'
    write_table(common_table,common_results)
    report=dict(science_measurement=False, calibration_applied=False,
                input_intake_sha256=sha(intake_path),input_support_sha256=sha(support_path),
                runner_sha256=sha(Path(__file__)),output_sha256={table.name:sha(table),common_table.name:sha(common_table)},
                n_results=len(results), n_fully_supported=sum(r['tract_multiplier'] is not None for r in results),
                limitation='Common-support estimates renormalize the observed histogram after excluding unsupported bins. '
                           'Tracts are subsets of the same COSMOS field and are not independent cosmic-variance realizations. '
                           'Missing target bins are not extrapolated; no area corrections, sampling errors, '
                           'reference-z errors or full selection compatibility are included.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
