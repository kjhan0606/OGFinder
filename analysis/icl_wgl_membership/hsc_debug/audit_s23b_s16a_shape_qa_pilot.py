#!/usr/bin/env python3
"""Offline intake QA of the approved S16A pilot, never a shear measurement.

Keeps all 66 ledger centers. The 23 zero-count centers are explicitly marked
as prior count-only zeros, not as newly queried or assigned zero shear.
Numeric-only support is distinct from a validated background-source sample.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from submit_s23b_s16a_shape_qa_pilot import validate_plan, verify_complete, numeric, FLOAT_COLUMNS
from fetch_ancillary_catalogues import catalogue_reader, sha
from audit_s23b_s16a_wl_random_support import require, describe, write_table
from prepare_s23b_s16a_shape_qa_schema_fix import OUTPUT as PLAN_DIRECTORY

OUTPUT = PLAN_DIRECTORY.parent/'wl_s23b_s16a_shape_qa_pilot_intake_2026-11-17_01'
# Seven existing WL-scale annuli plus a separate central raw-support bin.
RADIAL_EDGES = np.array([0., .1, .2, .4, .7, 1., 1.5, 2., 3.])
CALIBRATION_COLUMNS = FLOAT_COLUMNS[2:]


def radial_bins(radius, edges=RADIAL_EDGES):
    radius = np.asarray(radius, dtype=float)
    edges = np.asarray(edges, dtype=float)
    require(radius.ndim == 1 and np.isfinite(radius).all() and np.all(radius >= edges[0]), 'invalid radial coordinates')
    require(len(edges) >= 2 and np.isfinite(edges).all() and np.all(np.diff(edges) > 0), 'invalid radial edges')
    require(np.all(radius <= edges[-1]+1e-4), 'source lies outside proper aperture')
    # Include the final cone edge, allowing only the checked coordinate tolerance.
    return np.minimum(np.searchsorted(edges, radius, side='right')-1, len(edges)-2)


def weight_diagnostics(weights):
    w = np.asarray(weights, dtype=float)
    require(w.ndim == 1 and np.isfinite(w).all() and np.all(w > 0), 'diagnostic weights must be finite positive')
    if not len(w):
        return dict(sum_shape_weight=0., effective_number_shape_weight=0.)
    # Scaling avoids overflow; this is not SigmaCrit-weighted lensing support.
    scaled=w/np.max(w)
    effective=float(np.sum(scaled)**2/np.sum(scaled**2))
    require(1-1e-10 <= effective <= len(w)+1e-8, 'effective-number arithmetic failed')
    return dict(sum_shape_weight=float(np.sum(w)), effective_number_shape_weight=effective)


def weight_identity(sigma, rms, weights):
    """Descriptive check of the public Y1 inverse-variance identity, no cut."""
    sigma,rms,weights=(np.asarray(x,dtype=float) for x in (sigma,rms,weights))
    require(sigma.shape==rms.shape==weights.shape and sigma.ndim==1, 'weight identity input shapes differ')
    denominator=sigma*sigma+rms*rms
    eligible=np.isfinite(sigma)&np.isfinite(rms)&np.isfinite(weights)&np.isfinite(denominator)&(denominator>0)&(weights>0)
    residual=np.abs(weights[eligible]*denominator[eligible]-1)
    return dict(checked_rows=int(eligible.sum()),excluded_rows=int((~eligible).sum()),
        diagnostic_relative_tolerance=5e-6,above_tolerance=int((residual>5e-6).sum()),
        absolute_relative_residual=describe(residual),used_as_selection_cut=False)


def separations_radians(ra, dec, center_ra, center_dec):
    ra,dec=np.deg2rad(np.asarray(ra,float)),np.deg2rad(np.asarray(dec,float))
    ra0,dec0=np.deg2rad([center_ra,center_dec])
    h=np.sin((dec-dec0)/2)**2+np.cos(dec)*np.cos(dec0)*np.sin((ra-ra0)/2)**2
    return 2*np.arcsin(np.sqrt(np.clip(h,0,1)))


def common_support(center_profiles, centers):
    groups=defaultdict(list)
    for c in centers:groups[c['matched_lens_id']].append(c)
    profiles={(r['center_id'],r['radial_bin']):r for r in center_profiles}
    require(len(profiles)==len(center_profiles), 'duplicate radial profile')
    result=[]
    for lens,group in sorted(groups.items()):
        target=[c for c in group if c['role']=='lens'];randoms=[c for c in group if c['role']=='random']
        require(len(target)==1 and len(randoms)==5, 'expected one lens and five random centers')
        for k in range(len(RADIAL_EDGES)-1):
            t=profiles[(target[0]['center_id'],k)]
            rs=[profiles[(r['center_id'],k)] for r in randoms]
            result.append(dict(matched_lens_id=lens,icl_filter=group[0]['icl_filter'],radial_bin=k,
                r_min_mpc=float(RADIAL_EDGES[k]),r_max_mpc=float(RADIAL_EDGES[k+1]),
                lens_raw_count=t['raw_count'],n_randoms_with_raw_rows=sum(r['raw_count']>0 for r in rs),
                lens_numeric_positive_weight=t['numeric_positive_weight_count'],
                n_randoms_with_numeric_positive_weight=sum(r['numeric_positive_weight_count']>0 for r in rs),
                raw_common_support=t['raw_count']>0 and all(r['raw_count']>0 for r in rs),
                numeric_weight_common_support=t['numeric_positive_weight_count']>0 and all(r['numeric_positive_weight_count']>0 for r in rs)))
    return result


def reuse_audit(center_ids, object_ids, values, center_map):
    groups=defaultdict(list)
    for i,obj in enumerate(object_ids):groups[obj].append(i)
    require(len(set(zip(center_ids,object_ids))) == len(object_ids), 'duplicate center/object pair')
    inconsistent=across_lenses=across_roles=own_target_random=0
    for indices in groups.values():
        if len(indices)<2:continue
        first=values[indices[0]]
        if any(not np.array_equal(first,values[i],equal_nan=True) for i in indices[1:]):inconsistent+=1
        cs=[center_map[center_ids[i]] for i in indices]
        across_lenses+=len({c['matched_lens_id'] for c in cs})>1
        across_roles+=len({c['role'] for c in cs})>1
        own=defaultdict(set)
        for c in cs:own[c['matched_lens_id']].add(c['role'])
        own_target_random+=any(roles=={'lens','random'} for roles in own.values())
    sizes=Counter(len(indices) for indices in groups.values())
    return dict(source_pair_rows=len(object_ids),unique_objects=len(groups),
        pairs_minus_unique=len(object_ids)-len(groups),objects_in_multiple_centers=sum(n for size,n in sizes.items() if size>1),
        maximum_centers_per_object=max(sizes,default=0),multiplicity_histogram=dict(sorted(sizes.items())),
        objects_reused_across_lens_groups=across_lenses,objects_seen_in_lens_and_random_roles=across_roles,
        objects_in_own_lens_and_random=own_target_random,inconsistent_source_properties=inconsistent)


def analyze(plan, records):
    centers=plan['centers'];center_map={r['center_id']:r for r in centers}
    require(len(center_map)==66, 'original center ledger differs')
    ids=[r[0] for r in records];objects=[r[1] for r in records]
    values=np.asarray([r[2] for r in records],dtype=float)
    require(values.shape==(len(records),10), 'source column matrix differs')
    per_center=defaultdict(list)
    for i,ident in enumerate(ids):
        require(ident in center_map and center_map[ident]['query_planned'], 'unexpected/known-zero source center')
        per_center[ident].append(i)
    require(len(records)==plan['expected_source_pairs_from_prior_counts'], 'source-pair total differs')
    numeric_mask=np.isfinite(values[:,2:]).all(axis=1)&(values[:,6]>0)
    cosmology=FlatLambdaCDM(H0=70.,Om0=.3,Tcmb0=0)
    center_table=[];radial_table=[]
    for c in centers:
        indices=np.asarray(per_center[c['center_id']],dtype=int)
        require(len(indices)==c['expected_raw_pairs'], 'center row count differs')
        ra,dec=values[indices,0],values[indices,1]
        radius=separations_radians(ra,dec,c['ra_deg'],c['dec_deg'])*float(cosmology.angular_diameter_distance(c['z_lens']).value)
        bins=radial_bins(radius)
        valid=numeric_mask[indices]
        local_weights=values[indices,6]
        origin='received_source_rows' if c['query_planned'] else 'prior_count_only_zero_not_requeried'
        center_table.append(dict(center_id=c['center_id'],matched_lens_id=c['matched_lens_id'],role=c['role'],
            icl_filter=c['icl_filter'],source_origin=origin,expected_raw_count=c['expected_raw_pairs'],
            received_raw_count=len(indices),numeric_positive_weight_count=int(valid.sum()),
            **weight_diagnostics(local_weights[valid])))
        for k in range(len(RADIAL_EDGES)-1):
            selected=bins==k;good=selected&valid
            radial_table.append(dict(center_id=c['center_id'],matched_lens_id=c['matched_lens_id'],role=c['role'],
                icl_filter=c['icl_filter'],source_origin=origin,radial_bin=k,
                r_min_mpc=float(RADIAL_EDGES[k]),r_max_mpc=float(RADIAL_EDGES[k+1]),
                raw_count=int(selected.sum()),numeric_positive_weight_count=int(good.sum()),
                **weight_diagnostics(local_weights[good])))
    require(sum(r['raw_count'] for r in radial_table)==len(records), 'radial bin accounting differs')
    distributions={}
    for j,name in enumerate(FLOAT_COLUMNS):
        a=values[:,j];finite=np.isfinite(a)
        distributions[name]=dict(nonfinite=int((~finite).sum()),negative=int((a[finite]<0).sum()),
                                 zero=int((a[finite]==0).sum()),finite=describe(a[finite]))
    return dict(center_table=center_table,radial_table=radial_table,common=common_support(radial_table,centers),
        distributions=distributions,numeric_positive_weight_rows=int(numeric_mask.sum()),
        weight_identity=weight_identity(values[:,4],values[:,5],values[:,6]),
        reuse=reuse_audit(ids,objects,values,center_map),plot_data=(values[:,6],np.asarray([center_map[i]['icl_filter'] for i in ids])))


def make_plot(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11})
    fig,(weight_ax,support_ax)=plt.subplots(1,2,figsize=(12.5,5.2))
    weights,bands=result['plot_data'];valid=np.isfinite(weights)&(weights>0)
    require(valid.any(), 'no finite positive shape weights to visualize')
    low,high=np.min(np.log10(weights[valid])),np.max(np.log10(weights[valid]))
    if low==high:low-=.1;high+=.1
    bins=np.linspace(low,high,35)
    max_height=0.
    for band,color in (('i','#2563eb'),('i2','#c2410c')):
        values=np.log10(weights[valid&(bands==band)])
        counts,_,_=weight_ax.hist(values,bins=bins,histtype='step',lw=1.6,color=color,label='Lens ICL '+band)
        max_height=max(max_height,float(counts.max(initial=0)))
    weight_ax.set(xlabel=r'$\log_{10} w_\mathrm{shape}$',ylabel='Source pairs',ylim=(0,max_height*1.45))
    weight_ax.text(.04,.96,'Positive catalog shape weight',transform=weight_ax.transAxes,va='top')
    weight_ax.legend(loc='upper right',bbox_to_anchor=(.98,.86),framealpha=.9,fontsize=11)
    x=np.sqrt(RADIAL_EDGES[:-1]*RADIAL_EDGES[1:]);x[0]=.05
    for field,color,marker,label in (('raw_common_support','#2563eb','o','Raw rows'),
                                    ('numeric_weight_common_support','#c2410c','x','Finite calibration + w>0')):
        n=[sum(r[field] for r in result['common'] if r['radial_bin']==k) for k in range(len(x))]
        support_ax.plot(x,n,color=color,marker=marker,ms=5,lw=1.3,label=label)
    support_ax.set(xscale='log',xlabel='R [proper Mpc]',ylabel='Lens groups with all five controls',ylim=(-.3,11.8))
    support_ax.text(.04,.96,'Diagnostic radial support (11 groups)',transform=support_ax.transAxes,va='top')
    support_ax.legend(loc='upper right',bbox_to_anchor=(.98,.86),framealpha=.9,fontsize=11)
    fig.subplots_adjust(left=.08,right=.985,bottom=.15,top=.97,wspace=.28)
    for suffix in ('.png','.pdf'):fig.savefig(output/('pilot_intake_qa'+suffix),dpi=180)
    plt.close(fig)


def run(plan_directory, output):
    plan_directory=Path(plan_directory).resolve();output=Path(output).resolve()
    require(not output.exists(),'refusing to overwrite QA output')
    plan=validate_plan(plan_directory);results=plan_directory/'results'
    state=json.loads((results/'shape_qa_status.json').read_text())
    require(state['status']=='complete' and sha(results/'download_summary.json')==state['summary_sha256'], 'download is not complete/verified')
    receipt=verify_complete(plan,results,state)
    records=[]
    for batch in plan['batches']:
        file=results/(Path(batch['name']).stem+'.csv')
        with file.open(newline='',encoding='utf-8-sig') as f:
            for row in catalogue_reader(f):
                records.append((row['center_id'],str(int(row['object_id'])),[numeric(row[k]) for k in FLOAT_COLUMNS]))
    result=analyze(plan,records)
    output.mkdir(parents=True,exist_ok=False)
    write_table(output/'per_center_intake.tsv',result['center_table'])
    write_table(output/'per_center_radial_support.tsv',result['radial_table'])
    write_table(output/'per_lens_common_radial_support.tsv',result['common'])
    make_plot(result,output)
    radial_summary=[dict(radial_bin=k,r_min_mpc=float(RADIAL_EDGES[k]),r_max_mpc=float(RADIAL_EDGES[k+1]),
        n_lens_groups=11,raw_common_groups=sum(r['raw_common_support'] for r in result['common'] if r['radial_bin']==k),
        numeric_weight_common_groups=sum(r['numeric_weight_common_support'] for r in result['common'] if r['radial_bin']==k))
        for k in range(len(RADIAL_EDGES)-1)]
    summary=dict(schema='astrafek-s16a-shape-pilot-intake-qa-v1',status='complete_diagnostic_intake',
        science_measurement=False,network_used=False,new_query_submitted=False,
        authorized_scope_plan_sha256=receipt['approved_scope_plan_sha256'],executed_plan_sha256=plan['query_plan_sha256'],
        original_ledger_centers=66,queried_centers=43,prior_count_only_zero_centers=23,
        source_pair_rows=len(records),numeric_positive_weight_rows=result['numeric_positive_weight_rows'],
        source_column_distributions=result['distributions'],source_reuse=result['reuse'],radial_common_support=radial_summary,
        catalog_weight_identity=result['weight_identity'],
        weight_identity_reference=dict(url='https://github.com/PrincetonUniversity/hsc-y1-shear-calib/blob/main/utilities.py',
            git_blob_sha='55a640a6b3327b7fe0b7e3e2feb82dfecd3224b9',expression='weight * (sigma_e^2 + rms_e^2) = 1',external_code_imported=False),
        coordinates='Raw c.ira/c.idec interpreted as degrees only after public-schema unit and all source-cone checks; no conversion used.',
        numeric_support_definition='All eight received shape/calibration fields finite and catalog shape weight >0; no background/quality/selection claim.',
        byte_receipt={k:receipt[k] for k in ('verified_csv_bytes','confirmed_bytes_received','budget_bytes_charged','maximum_total_received_bytes')},
        input_sha256=dict(plan_json=sha(plan_directory/'shape_qa_plan.json'),download_summary=sha(results/'download_summary.json')),
        source_csv_sha256=receipt['output_sha256'],script_sha256=sha(Path(__file__)),
        runtime=dict(python=sys.version,numpy=np.__version__),
        output_sha256={f.name:sha(f) for f in sorted(output.iterdir()) if f.is_file()},
        limitations=['Diagnostic pilot chosen by raw-support strata, not a representative sample or final scientific selection.',
            'Known zero centers are inherited from earlier count-only queries, not freshly queried empty CSVs.',
            'Shape-weight effective number is a numerical concentration diagnostic, not a calibrated background-source/lensing effective number.',
            'Total and radial raw counts do not measure usable unmasked area or validate matched observational selection.',
            'No redshift, photo-z, PDF, full Y1 flags/resolution/photometry or PSF-moment fields were received.',
            'No shear, responsivity, mass, boost, covariance, n(z), source-redshift calibration or ICL-WL correlation was measured.',
            'No S19A/Y3 cut/selection coefficient or redshift calibration was imposed on S16A.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'],output=str(output),summary_sha256=sha(output/'summary.json'),
        rows=len(records),numeric_positive_weight_rows=result['numeric_positive_weight_rows'],source_reuse=result['reuse'],
        radial_common_support=radial_summary,bytes=summary['byte_receipt']),indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-dir',type=Path,default=PLAN_DIRECTORY)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();run(args.plan_dir,args.output)
