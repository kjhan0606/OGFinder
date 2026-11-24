#!/usr/bin/env python3
"""Propagate conditional COSMOS multiplier covariance to diagnostic stacks."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import sha, write_table

DATA=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SUPPORT=DATA/'s16a_pointz_background_support_2026-11-18_01'
COV=DATA/'s16a_cosmos_bootstrap_covariance_2026-11-21_01'


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summaries={}
    for folder in (SUPPORT,COV):
        payload=json.loads((folder/'summary.json').read_text())
        for name, expected in payload['output_sha256'].items():
            if sha(folder/name)!=expected:
                raise ValueError('Changed diagnostic output')
        summaries[folder]=payload
    profile=read(SUPPORT/'center_annulus_pointz_support.tsv')
    lenses=read(COV/'lens_bootstrap.tsv')
    matrix=read(COV/'lens_bootstrap_covariance.tsv')
    ids=[r['lens_id'] for r in lenses]
    if len(ids)!=11 or len(set(ids))!=11:
        raise ValueError('Unexpected calibration lens IDs')
    indices={name:i for i,name in enumerate(ids)}
    cov=np.full((len(ids),len(ids)),np.nan)
    for row in matrix:
        i,j=indices[row['lens_a']],indices[row['lens_b']]
        if np.isfinite(cov[i,j]):
            raise ValueError('Duplicate covariance entry')
        cov[i,j]=float(row['covariance'])
    if not np.all(np.isfinite(cov)) or not np.allclose(cov,cov.T) or np.linalg.eigvalsh(cov).min() < -1e-14:
        raise ValueError('Invalid covariance matrix')
    multipliers=np.asarray([float(r['baseline']) for r in lenses])
    common=summaries[SUPPORT]['per_annulus_common_lens_ids']['0.2']
    edges=summaries[SUPPORT]['radial_edges_proper_mpc']
    results=[]
    gradients=[]
    for k,members in enumerate(common):
        weights=np.zeros(len(ids))
        for row in profile:
            if row['role']=='lens' and int(row['annulus'])==k and float(row['gap'])==0.2 and row['lens_id'] in members:
                weights[indices[row['lens_id']]]+=float(row['pair_weight_sum'])
        if len(members)!=np.count_nonzero(weights):
            raise ValueError('Missing lens weight')
        weights/=weights.sum()
        correction=1/np.sum(weights/multipliers)
        gradient=weights*correction**2/multipliers**2
        variance=float(gradient@cov@gradient)
        independent=float(np.sum(gradient**2*np.diag(cov)))
        if variance<0 or independent<=0:
            raise ValueError('Invalid propagated variance')
        gradients.append(gradient)
        results.append(dict(annulus=k,r_min_mpc=edges[k],r_max_mpc=edges[k+1],
                            lens_ids='|'.join(members),n_lenses=len(members),
                            diagnostic_stack_multiplier=correction,
                            conditional_sd=float(np.sqrt(variance)),
                            independent_sd=float(np.sqrt(independent)),
                            sd_ratio=float(np.sqrt(variance/independent))))
    cross=np.asarray(gradients)@cov@np.asarray(gradients).T
    output.mkdir(parents=True)
    a,b=output/'radial_stack_calibration.tsv',output/'radial_calibration_covariance.tsv'
    write_table(a,results)
    write_table(b,[dict(annulus_a=i,annulus_b=j,covariance=float(cross[i,j]))
                   for i in range(len(results)) for j in range(len(results))])
    report=dict(science_measurement=False,calibration_applied=False,
                support_summary_sha256=sha(SUPPORT/'summary.json'),covariance_summary_sha256=sha(COV/'summary.json'),
                runner_sha256=sha(Path(__file__)),output_sha256={a.name:sha(a),b.name:sha(b)},
                method='Linear propagation of conditional calibration multiplier covariance to '
                       'inverse-response weighted effective multipliers using diagnostic lens-pair weights.',
                limitation='Uses per-lens overall photo-z responses for all annuli, fixed diagnostic weights '
                           'and matched-support sets. It is not a calibrated shear stack or full covariance; '
                           'selection, radial redshift dependence, reference errors and cosmic variance remain.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
