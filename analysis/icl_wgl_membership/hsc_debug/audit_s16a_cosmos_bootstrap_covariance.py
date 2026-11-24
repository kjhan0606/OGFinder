#!/usr/bin/env python3
"""Conditional shared-COSMOS object bootstrap across pilot lenses."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, sha, write_table
from ogfkit.wlcluster import Cosmology

SUPPORT = INTAKE.parent/'s16a_cosmos_bin_support_2026-11-21_01'
REPLICAS = 2000
SEED = 2026112101


def run(output):
    if output.exists():
        raise FileExistsError(output)
    intake_path, support_path = INTAKE/'summary.json', SUPPORT/'summary.json'
    intake, support = (json.loads(p.read_text()) for p in (intake_path, support_path))
    for directory, report in ((INTAKE, intake), (SUPPORT, support)):
        for filename, expected in report['output_sha256'].items():
            if sha(directory/filename) != expected:
                raise ValueError('Changed diagnostic input')
    if support['intake_sha256'] != sha(intake_path):
        raise ValueError('Input lineage differs')
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        data = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    if len(data) != intake['eligible_rows']:
        raise ValueError('Unexpected eligible row count')
    photo, ref, shape, som = np.asarray([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in data]).T
    with (SUPPORT/'lens_photoz_reweight.tsv').open() as stream:
        lenses = list(csv.DictReader(stream, delimiter='\t'))
    with (SUPPORT/'photoz_bin_support.tsv').open() as stream:
        bins = list(csv.DictReader(stream, delimiter='\t'))
    edges = np.asarray(support['edges'])
    model = Cosmology()
    prepared = []
    for lens in lenses:
        zl=float(lens['lens_z'])
        active=photo > zl+0.2
        idx=np.searchsorted(edges, photo[active], side='right')-1
        chi=float(model.comoving(zl))
        est=1-chi/model.comoving(photo[active])
        true=np.zeros(active.sum())
        behind=ref[active]>zl
        true[behind]=1-chi/model.comoving(ref[active][behind])
        w=shape[active]*som[active]*est**2
        target=np.asarray([float(r['observed_weight_fraction']) for r in bins if r['lens_id']==lens['lens_id']])
        if len(target)!=len(edges)-1 or not np.isclose(target.sum(),1):
            raise ValueError('Invalid target histogram')
        mass=np.bincount(idx,weights=w,minlength=len(target))
        num=np.bincount(idx,weights=w*true/est,minlength=len(target))
        avg=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
        if np.any((mass==0)&(target>0)) or not np.isclose(1/(target@avg),float(lens['reweighted']),rtol=1e-12):
            raise ValueError('Baseline mismatch')
        prepared.append((active,idx,w,w*true/est,target))
    if len(prepared)!=11:
        raise ValueError('Expected 11 lenses')
    rng=np.random.default_rng(SEED)
    draws=np.full((REPLICAS,len(lenses)),np.nan)
    for i in range(REPLICAS):
        counts=np.bincount(rng.integers(0,len(data),size=len(data)),minlength=len(data))
        for j,(active,idx,w,numerator,target) in enumerate(prepared):
            mult=counts[active]
            mass=np.bincount(idx,weights=w*mult,minlength=len(target))
            if np.any((mass==0)&(target>0)):
                continue
            num=np.bincount(idx,weights=numerator*mult,minlength=len(target))
            avg=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
            response=target@avg
            if response>0:
                draws[i,j]=1/response
    valid=np.isfinite(draws).all(axis=1)
    if valid.sum()<REPLICAS*0.9:
        raise ValueError('Insufficient jointly valid replicates')
    x=draws[valid]
    covariance=np.cov(x,rowvar=False,ddof=1)
    correlation=np.corrcoef(x,rowvar=False)
    output.mkdir(parents=True)
    lens_ids=[r['lens_id'] for r in lenses]
    per_lens=[dict(lens_id=lid, baseline=float(lenses[i]['reweighted']), sd=float(np.std(x[:,i],ddof=1)),
                   low=float(np.percentile(x[:,i],2.5)), high=float(np.percentile(x[:,i],97.5)))
              for i,lid in enumerate(lens_ids)]
    matrix=[dict(lens_a=a,lens_b=b,covariance=float(covariance[i,j]),correlation=float(correlation[i,j]))
            for i,a in enumerate(lens_ids) for j,b in enumerate(lens_ids)]
    a,b=output/'lens_bootstrap.tsv',output/'lens_bootstrap_covariance.tsv'
    write_table(a,per_lens)
    write_table(b,matrix)
    report=dict(science_measurement=False, calibration_applied=False, n_replicates=REPLICAS,
                valid_joint=int(valid.sum()), seed=SEED, n_lenses=len(lenses),
                input_intake_sha256=sha(intake_path),input_support_sha256=sha(support_path),
                runner_sha256=sha(Path(__file__)),output_sha256={a.name:sha(a),b.name:sha(b)},
                method='Shared fixed-size object bootstrap of the joined COSMOS reference rows, '
                       'with fixed observed lens histograms and re-fitted bin weights.',
                limitation='Conditional IID resampling; does not include spatial covariance, cosmic variance, '
                           'reference-z errors, observed histogram uncertainty or selection mismatch.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
