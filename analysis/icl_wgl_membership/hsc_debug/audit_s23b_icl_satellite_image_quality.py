#!/usr/bin/env python3
"""Check ICL–satellite slope sensitivity to depth and seeing covariates."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_field_omission import (
    DATA, FLOORS, N_BOOTSTRAP, ORIGIN, RELATION, RELATION_SHA,
    block_bootstrap_slope, connected_pair_blocks, model_matrix, read_tsv, sha256, write_tsv,
)

CONTROL = DATA/'s23b_icl_control_image_plan_2026-11-07_02/control_manifest.json'
CONTROL_SHA = 'ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    relation_path=RELATION/'summary.json'
    if sha256(relation_path)!=RELATION_SHA or sha256(CONTROL)!=CONTROL_SHA:
        raise ValueError('Changed pinned input')
    relation=json.loads(relation_path.read_text())
    table=RELATION/'paired_icl_satellite_profiles.tsv'
    if sha256(table)!=relation['output_sha256'][table.name]:
        raise ValueError('Changed relation table')
    manifest=json.loads(CONTROL.read_text())
    controls={r['id']:r for r in manifest['clusters']}
    if len(controls)!=len(manifest['clusters']):
        raise ValueError('Duplicate control IDs')
    rows=[r for r in read_tsv(table) if r['icl_relative_contrast_njy_arcsec2']]
    if len(rows)!=988:
        raise ValueError('Unexpected supported sample')
    for row in rows:
        c=controls[row['control_id']]
        if c['target_cluster_id']!=row['target_id']:
            raise ValueError('Control target mismatch')
        row['target_depth']=float(c['target_depth_mag'])
        row['depth_difference']=float(c['delta_depth_mag'])
        row['target_seeing']=float(c['target_seeing_arcsec'])
        row['seeing_ratio']=float(c['seeing_ratio'])
    results=[]
    for band_index,band in enumerate(('i','i2')):
        for floor_index,floor in enumerate(FLOORS):
            sample=[r for r in rows if r['filter']==band and min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=floor]
            X,y,_=model_matrix(sample)
            quality=np.asarray([[r[k] for k in ('target_depth','depth_difference','target_seeing','seeing_ratio')] for r in sample])
            if not np.all(np.isfinite(quality)):
                raise ValueError('Missing image quality')
            quality-=quality.mean(axis=0)
            quality/=quality.std(axis=0)
            expanded=np.column_stack((X,quality))
            controls_by_target={r['target_id']:dict(ra_deg=r['control_ra_deg'],dec_deg=r['control_dec_deg']) for r in sample}
            labels,blocks=connected_pair_blocks(sample,controls_by_target,5.,0.,ORIGIN)
            for name,design in (('baseline',X),('image_quality',expanded)):
                record=dict(filter=band,support_floor=floor,model=name,n_pairs=len(sample),
                            n_blocks=len(blocks),n_columns=design.shape[1],slope=None,
                            ci_low=None,ci_high=None,valid_bootstrap=0,status='not_estimable')
                try:
                    slope,low,high,_,count=block_bootstrap_slope(design,y,labels,N_BOOTSTRAP,2026112100+band_index*10+floor_index)
                except ValueError as error:
                    record['reason']=str(error)
                else:
                    record.update(slope=slope,ci_low=low,ci_high=high,valid_bootstrap=count,status='complete',reason='')
                results.append(record)
    output.mkdir(parents=True)
    dest=output/'quality_covariate_sensitivity.tsv'
    write_tsv(dest,results)
    summary=dict(science_measurement=False,n_models=len(results),bootstrap_replicas=N_BOOTSTRAP,
                 relation_sha256=RELATION_SHA,control_manifest_sha256=CONTROL_SHA,
                 runner_sha256=sha256(Path(__file__)),output_sha256={dest.name:sha256(dest)},
                 limitation='Exploratory adjustment for target depth, depth difference, target seeing, '
                            'and seeing ratio; errors are 5-degree connected-block percentiles with '
                            'one grid phase, and do not include selection, PSF-wing or pixel covariance systematics.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
