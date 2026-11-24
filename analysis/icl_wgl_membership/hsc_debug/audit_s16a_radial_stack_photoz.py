#!/usr/bin/env python3
"""Compare radial and all-radius estimated-z calibration in diagnostic stacks."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import sha, write_table

DATA=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SUPPORT=DATA/'s16a_pointz_background_support_2026-11-18_01'
RADIAL=DATA/'s16a_radial_photoz_response_2026-11-21_01'
STACK=DATA/'s16a_calibration_stack_covariance_2026-11-21_01'


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def run(output, fixed_common=False):
    if output.exists():
        raise FileExistsError(output)
    summaries={}
    for folder in (SUPPORT,RADIAL,STACK):
        summary=json.loads((folder/'summary.json').read_text())
        for name,digest in summary['output_sha256'].items():
            if sha(folder/name)!=digest:
                raise ValueError('Changed input artifact')
        summaries[folder]=summary
    if summaries[RADIAL]['support_sha256']!=sha(SUPPORT/'summary.json') or summaries[STACK]['support_summary_sha256']!=sha(SUPPORT/'summary.json'):
        raise ValueError('Support lineage mismatch')
    radial=read(RADIAL/'lens_radial_photoz_response.tsv')
    profile=read(SUPPORT/'center_annulus_pointz_support.tsv')
    old=read(STACK/'radial_stack_calibration.tsv')
    lookup={(r['lens_id'],int(r['annulus'])):r for r in radial}
    if len(lookup)!=len(radial):
        raise ValueError('Duplicate radial results')
    weights={(r['lens_id'],int(r['annulus'])):float(r['pair_weight_sum'])
             for r in profile if r['role']=='lens' and float(r['gap'])==0.2}
    result=[]
    fixed=summaries[SUPPORT]['common_0p1_to_3']['0.2']
    for baseline in old:
        annulus=int(baseline['annulus'])
        if fixed_common and annulus==0:
            continue
        members=fixed if fixed_common else baseline['lens_ids'].split('|')
        if not members or any((x,annulus) not in lookup or (x,annulus) not in weights for x in members):
            raise ValueError('Missing supported lens')
        rows=[lookup[x,annulus] for x in members]
        if any(not r['radial_multiplier'] for r in rows):
            raise ValueError('Unsupported calibration input')
        w=np.array([weights[x,annulus] for x in members])
        if np.any(w<=0):
            raise ValueError('Invalid pair weights')
        w/=w.sum()
        full=np.array([float(r['all_radius_multiplier']) for r in rows])
        radial_m=np.array([float(r['radial_multiplier']) for r in rows])
        full_stack=float(1/(w@(1/full)))
        annular_stack=float(1/(w@(1/radial_m)))
        if members==baseline['lens_ids'].split('|') and not np.isclose(full_stack,float(baseline['diagnostic_stack_multiplier']),rtol=1e-12):
            raise ValueError('Previous stack multiplier differs')
        result.append(dict(annulus=annulus,r_min_mpc=baseline['r_min_mpc'],r_max_mpc=baseline['r_max_mpc'],
                           n_lenses=len(members),lens_ids='|'.join(members),
                           full_radius_multiplier=full_stack,annular_multiplier=annular_stack,
                           fractional_difference=annular_stack/full_stack-1,
                           min_selected_pairs=min(int(r['selected_pairs']) for r in rows)))
    output.mkdir(parents=True)
    table=output/'radial_stack_sensitivity.tsv'
    write_table(table,result)
    summary=dict(science_measurement=False,calibration_applied=False,fixed_common=fixed_common,
                 inputs={str(folder/'summary.json'):sha(folder/'summary.json') for folder in summaries},
                 runner_sha256=sha(Path(__file__)),output_sha256={table.name:sha(table)},
                 limitation='Fixed diagnostic pair weights and exploratory annulus-specific point-z responses. '
                            'No measured shear or statistical uncertainty is propagated. Small annuli, '
                            'selection mismatch, reference-z errors and cosmic variance remain.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--fixed-common',action='store_true')
    args=parser.parse_args()
    run(args.output.resolve(),args.fixed_common)
