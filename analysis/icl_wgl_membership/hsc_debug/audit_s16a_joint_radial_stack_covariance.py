#!/usr/bin/env python3
"""Propagate joint lens/annulus calibration covariance to diagnostic stacks."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import sha, write_table

DATA=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
JOINT=DATA/'s16a_radial_joint_bootstrap_2026-11-21_01'
RADIAL=DATA/'s16a_radial_stack_photoz_2026-11-21_01'
SUPPORT=DATA/'s16a_pointz_background_support_2026-11-18_01'
OLD=DATA/'s16a_calibration_stack_covariance_2026-11-21_01'


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    reports={}
    for folder in (JOINT,RADIAL,SUPPORT,OLD):
        report=json.loads((folder/'summary.json').read_text())
        for name,digest in report['output_sha256'].items():
            if sha(folder/name)!=digest:
                raise ValueError('Changed input artifact')
        reports[folder]=report
    if reports[JOINT]['inputs'][str(SUPPORT/'summary.json')]!=sha(SUPPORT/'summary.json') or reports[OLD]['support_summary_sha256']!=sha(SUPPORT/'summary.json'):
        raise ValueError('Support provenance differs')
    lens=read(JOINT/'joint_bootstrap_lens_annulus.tsv')
    keys=[(r['lens_id'],int(r['annulus'])) for r in lens]
    if len(keys)!=24 or len(set(keys))!=24:
        raise ValueError('Unexpected joint cases')
    index={key:i for i,key in enumerate(keys)}
    covariance=np.full((24,24),np.nan)
    for row in read(JOINT/'joint_bootstrap_covariance.tsv'):
        i=index[row['lens_a'],int(row['annulus_a'])]
        j=index[row['lens_b'],int(row['annulus_b'])]
        if np.isfinite(covariance[i,j]):
            raise ValueError('Duplicate covariance entry')
        covariance[i,j]=float(row['covariance'])
    if not np.all(np.isfinite(covariance)) or not np.allclose(covariance,covariance.T) or np.linalg.eigvalsh(covariance).min() < -1e-13:
        raise ValueError('Invalid joint covariance')
    profile=read(SUPPORT/'center_annulus_pointz_support.tsv')
    weights={(r['lens_id'],int(r['annulus'])):float(r['pair_weight_sum'])
             for r in profile if r['role']=='lens' and float(r['gap'])==0.2}
    old={int(r['annulus']):r for r in read(OLD/'radial_stack_calibration.tsv')}
    results=[]
    gradients=[]
    for row in read(RADIAL/'radial_stack_sensitivity.tsv'):
        annulus=int(row['annulus'])
        members=row['lens_ids'].split('|')
        positions=[index[lid,annulus] for lid in members]
        w=np.asarray([weights[lid,annulus] for lid in members])
        if np.any(w<=0):
            raise ValueError('Nonpositive lens weight')
        w/=w.sum()
        mult=np.asarray([float(lens[i]['baseline']) for i in positions])
        stack=float(1/(w@(1/mult)))
        if not np.isclose(stack,float(row['annular_multiplier']),rtol=1e-12):
            raise ValueError('Previous radial stack differs')
        gradient=np.zeros(24)
        gradient[positions]=w*stack**2/mult**2
        var=float(gradient@covariance@gradient)
        diag=float(np.sum(gradient**2*np.diag(covariance)))
        if var<0 or diag<=0:
            raise ValueError('Invalid propagated variance')
        gradients.append(gradient)
        results.append(dict(annulus=annulus,r_min_mpc=row['r_min_mpc'],r_max_mpc=row['r_max_mpc'],
                            lens_ids=row['lens_ids'],annular_multiplier=stack,
                            conditional_sd=float(np.sqrt(var)),independent_sd=float(np.sqrt(diag)),
                            previous_full_radius_sd=float(old[annulus]['conditional_sd'])))
    gradients=np.asarray(gradients)
    propagated=gradients@covariance@gradients.T
    output.mkdir(parents=True)
    a,b=output/'joint_stack_uncertainty.tsv',output/'joint_stack_covariance.tsv'
    write_table(a,results)
    write_table(b,[dict(annulus_a=i,annulus_b=j,covariance=float(propagated[i,j]))
                   for i in range(len(results)) for j in range(len(results))])
    summary=dict(science_measurement=False,calibration_applied=False,
                 inputs={str(folder/'summary.json'):sha(folder/'summary.json') for folder in reports},
                 runner_sha256=sha(Path(__file__)),output_sha256={a.name:sha(a),b.name:sha(b)},
                 limitation='Linear propagation of conditional source/reference Poisson-bootstrap multiplier '
                            'covariance with lens pair weights held fixed. Weight fluctuations and their '
                            'covariance with multipliers are excluded. Spatial clustering, selection '
                            'uncertainty, cosmic variance and measured shear are not included.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
