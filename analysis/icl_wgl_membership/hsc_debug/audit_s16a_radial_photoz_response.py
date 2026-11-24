#!/usr/bin/env python3
"""Assess radial sensitivity of the COSMOS photo-z response in the S16A pilot."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from audit_s16a_y1_reuse_shape_noise import GEOMETRY, GEOMETRY_SHA
from ogfkit.wlcluster import Cosmology

DATA=INTAKE.parent
SUPPORT=DATA/'s16a_pointz_background_support_2026-11-18_01'
CALIB=DATA/'s16a_cosmos_bin_support_2026-11-21_01'


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    intake_file, support_file, cal_file = INTAKE/'summary.json', SUPPORT/'summary.json', CALIB/'summary.json'
    intake,support,cal=(json.loads(p.read_text()) for p in (intake_file,support_file,cal_file))
    for folder,report in ((INTAKE,intake),(SUPPORT,support),(CALIB,cal)):
        for name,expected in report['output_sha256'].items():
            if sha(folder/name)!=expected:
                raise ValueError('Changed input file')
    if cal['intake_sha256']!=sha(intake_file):
        raise ValueError('COSMOS lineage mismatch')
    centers,pairs,selected,available,_,provenance=load_inputs(PLAN,OBJECT_QA)
    for key,value in provenance.items():
        if support[key]!=value or cal[key]!=value:
            raise ValueError('Pilot lineage mismatch')
    point_z=finite_point_redshifts(available)
    geometry_path=GEOMETRY/'summary.json'
    if sha(geometry_path)!=GEOMETRY_SHA:
        raise ValueError('Geometry summary changed')
    geometry=json.loads(geometry_path.read_text())
    geometry_table=GEOMETRY/'pair_geometry.tsv'
    if sha(geometry_table)!=geometry['output_sha256'][geometry_table.name]:
        raise ValueError('Geometry table changed')
    groups=defaultdict(list)
    with geometry_table.open() as stream:
        for pair,geo in zip(pairs,csv.DictReader(stream,delimiter='\t'),strict=True):
            if (pair['center_id'],pair['object_id'])!=(geo['center_id'],geo['object_id']):
                raise ValueError('Pair/geometry mismatch')
            if pair['object_id'] in selected and pair['object_id'] in point_z:
                groups[pair['center_id'],int(geo['radial_bin'])].append((point_z[pair['object_id']],float(pair['weight'])))
    cat=[r for r in read(INTAKE/'cosmos_join.tsv') if r['eligible_for_selection_review']=='True']
    if len(cat)!=intake['eligible_rows']:
        raise ValueError('Eligibility mismatch')
    photo,ref,shape,som=np.asarray([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in cat]).T
    baselines={r['lens_id']:float(r['reweighted']) for r in read(CALIB/'lens_photoz_reweight.tsv')}
    common=support['per_annulus_common_lens_ids']['0.2']
    edges=np.asarray(cal['edges'])
    cosmo=Cosmology()
    rng=np.random.default_rng(2026112103)
    replicas=2000
    result=[]
    for center in centers:
        if center['role']!='lens':
            continue
        lid=center['matched_lens_id']
        zl=float(center['z_lens'])
        chi=float(cosmo.comoving(zl))
        keep=photo>zl+0.2
        p,r,sh,so=photo[keep],ref[keep],shape[keep],som[keep]
        est=1-chi/cosmo.comoving(p)
        true=np.zeros_like(r)
        back=r>zl
        true[back]=1-chi/cosmo.comoving(r[back])
        w=sh*so*est**2
        index=np.searchsorted(edges,p,side='right')-1
        mass=np.bincount(index,weights=w,minlength=len(edges)-1)
        num=np.bincount(index,weights=w*true/est,minlength=len(edges)-1)
        means=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
        for k,members in enumerate(common):
            if lid not in members:
                continue
            values=np.asarray(groups[center['center_id'],k],float).reshape(-1,2)
            values=values[values[:,0]>zl+0.2]
            z,obs_weight=values.T
            weights=obs_weight*(1-chi/cosmo.comoving(z))**2
            histogram=np.histogram(z,bins=edges,weights=weights)[0]
            if weights.sum()<=0 or not np.isclose(histogram.sum(),weights.sum()):
                raise ValueError('Missing or out-of-range observed weights')
            target=histogram/histogram.sum()
            missing=float(target[mass<=0].sum())
            multiplier=float(1/(target@means)) if missing==0 else None
            selected_bins=np.searchsorted(edges,z,side='right')-1
            if np.any((mass[selected_bins]<=0)):
                raise ValueError('Observed redshift has no COSMOS support')
            # Each draw resamples observed source pairs within this annulus.
            # COSMOS calibration means are fixed throughout this diagnostic.
            selections=rng.integers(0,len(z),size=(replicas,len(z)))
            sampled_weights=weights[selections]
            weighted_response=sampled_weights*means[selected_bins[selections]]
            response=weighted_response.sum(axis=1)/sampled_weights.sum(axis=1)
            if np.any(response<=0):
                raise ValueError('A resample has zero reference response')
            draws=1/response
            result.append(dict(lens_id=lid,annulus=k,selected_pairs=len(z),
                               pair_neff=float(weights.sum()**2/(weights@weights)),
                               unsupported_weight_fraction=missing,all_radius_multiplier=baselines[lid],
                               radial_multiplier=multiplier,
                               fractional_change=multiplier/baselines[lid]-1 if multiplier else None,
                               observed_bootstrap_sd=float(np.std(draws,ddof=1)),
                               observed_bootstrap_low=float(np.percentile(draws,2.5)),
                               observed_bootstrap_high=float(np.percentile(draws,97.5))))
    if len(result)!=sum(map(len,common)):
        raise ValueError('Common lens/annulus membership mismatch')
    output.mkdir(parents=True)
    table=output/'lens_radial_photoz_response.tsv'
    write_table(table,result)
    report=dict(science_measurement=False,calibration_applied=False,n_rows=len(result),
                bootstrap_replicas=replicas,bootstrap_seed=2026112103,
                n_supported=sum(r['radial_multiplier'] is not None for r in result),
                intake_sha256=sha(intake_file),support_sha256=sha(support_file),calibration_sha256=sha(cal_file),
                geometry_sha256=GEOMETRY_SHA,runner_sha256=sha(Path(__file__)),
                output_sha256={table.name:sha(table)},
                limitation='Exploratory point-z reweighting with fixed 0.1 bins. Per-annulus bootstrap '
                           'resamples observed pairs independently, fixes COSMOS bin responses, and is not '
                           'a cross-annulus covariance; reused sources and spatial clustering are ignored. Full selection compatibility, '
                           'reference redshift errors, spatial uncertainty, and shear covariance remain unresolved.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
