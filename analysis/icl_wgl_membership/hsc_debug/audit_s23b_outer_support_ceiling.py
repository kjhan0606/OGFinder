#!/usr/bin/env python3
"""Audit target and third-best control coverage using saved 50–75 kpc bins."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(root, out):
    if out.exists():
        raise FileExistsError(out)
    m = root/'measurement'
    ms = json.loads((m/'summary.json').read_text())
    audit = json.loads((root/'audit/summary.json').read_text())
    path = m/'member_profile_bin_diagnostics.tsv'
    if audit.get('historical_profile_rows_identical') is not True or audit.get('measurement_summary_sha256') != sha(m/'summary.json') or sha(path) != ms['output_sha256'][path.name]:
        raise ValueError('audited input mismatch')
    bands = {'i': [], 'i2': []}
    with path.open(newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['radial_bin_kpc'] != '50-75' or row['sep_matched'] != 'True':
                continue
            controls = sorted((float(c['coverage']) for c in json.loads(row['control_bin_measurements_json'])), reverse=True)
            target = float(row['target_coverage']) if row['target_profile_available'] == 'True' else 0.
            third = controls[2] if len(controls) >= 3 else 0.
            bands[row['filter']].append((target, third))
    result = {}
    for band,values in bands.items():
        arr=np.asarray(values)
        t,c=arr[:,0],arr[:,1]
        result[band]={'n_matched':len(arr),'target_at_least_0_8':int((t>=.8).sum()),
            'third_control_at_least_0_8':int((c>=.8).sum()),
            'both_at_least_0_8':int(((t>=.8)&(c>=.8)).sum()),
            'target_q16_median_q84_max':[float(x) for x in np.quantile(t,[.16,.5,.84,1])],
            'third_control_q16_median_q84_max':[float(x) for x in np.quantile(c,[.16,.5,.84,1])],
            'joint_ceiling_q16_median_q84_max':[float(x) for x in np.quantile(np.minimum(t,c),[.16,.5,.84,1])]}
    out.mkdir(parents=True,exist_ok=False)
    payload={'status':'complete_descriptive_support_ceiling','science_measurement':False,'images_read':False,'network_used':False,
             'diagnostics_sha256':sha(path),'audit_summary_sha256':sha(root/'audit/summary.json'),
             'script_sha256':sha(Path(__file__)),'bands':result,
             'limitation':'Conditional on uniquely SEP-matched candidates; coverage does not identify mask components or completeness.'}
    (out/'summary.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    return payload


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-root',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.run_root.resolve(),a.output.resolve()),indent=2))
