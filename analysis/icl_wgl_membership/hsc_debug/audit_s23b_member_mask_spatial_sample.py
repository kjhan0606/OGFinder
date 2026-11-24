#!/usr/bin/env python3
"""Audit base-mask reasons on a predeclared sky-stratified image sample.

Masks overlap; total footprints and exclusive pixels are reported separately.
The original source mask and saved retained profile rows must be reproduced.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import measure_s23b_camira_member_outer_profiles as pipeline
import measure_satellite_outer_profiles as core
from analyze_s23b_member_outer_icl_relation import MEMBER_BINS, EXPECTED_HASHES
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT


def run(out, plan_path, expected_plan_hash):
    if out.exists(): raise FileExistsError(out)
    for path, expected in ((pipeline.TARGET_MANIFEST,pipeline.EXPECTED_TARGET_MANIFEST_SHA256),
        (pipeline.PAIR_CONTRASTS,pipeline.EXPECTED_PAIR_CONTRASTS_SHA256),
        (pipeline.MEMBER_TABLE,pipeline.EXPECTED_MEMBER_TABLE_SHA256),
        (pipeline.MASK_RECORDS,pipeline.EXPECTED_MASK_RECORDS_SHA256),
        (MEMBER_BINS,EXPECTED_HASHES[str(MEMBER_BINS)])):
        if pipeline.sha256(path)!=expected: raise ValueError('pinned input changed')
    manifest=json.loads(pipeline.TARGET_MANIFEST.read_text())
    cosmo=FlatLambdaCDM(H0=float(manifest['h0']),Om0=float(manifest['omega_m']))
    frames={str(r['id']):r for r in manifest['clusters']}
    records={r['id']:r for r in pipeline.read_tsv(pipeline.MASK_RECORDS)}
    if pipeline.sha256(plan_path)!=expected_plan_hash: raise ValueError('selection plan hash changed')
    plan=json.loads(plan_path.read_text())
    if plan.get('schema')!='astrafek-s23b-member-mask-spatial-sample-v1' or plan.get('science_measurement') is not False or plan.get('images_read') is not False:
        raise ValueError('unexpected predeclared selection plan')
    for path,digest in plan['input_sha256'].items():
        if pipeline.sha256(Path(path))!=digest: raise ValueError('selection input hash changed')
    pool,_=pipeline.select_sample(pipeline.read_tsv(pipeline.PAIR_CONTRASTS),
        pipeline.read_tsv(pipeline.MEMBER_TABLE),frames,records,cosmo,None,plan['seed'])
    by_key={(r['target_id'],r['filter']):r for r in pool}
    requested=[(r['target_id'],r['filter']) for r in plan['sample']]
    if len(requested)!=96 or len(set(requested))!=96 or any(k not in by_key for k in requested):
        raise ValueError('selected sample does not match 96 unique eligible frames')
    selected=[by_key[k] for k in requested]
    if {b:sum(item['filter']==b for item in selected) for b in ('i','i2')}!={'i':48,'i2':48}:
        raise ValueError('sample must contain 48 frames per filter')
    for row in plan['sample']:
        f=frames[row['target_id']]
        if float(f['ra_deg'])!=row['ra_deg'] or float(f['dec_deg'])!=row['dec_deg']:
            raise ValueError('selected coordinates changed')
    old_mask,old_profile=core.make_stamp_mask,core.profile_at
    rows=[]; state={}; zeros={}
    def instrument_profile(image,variance,base,segmentation,detections,fallback_tree,
                           center_xy,pixel_scale,kpc_per_pixel,conversion,target_number,
                           coverage_floor,kron_factor=core.KRON_FACTOR):
        state.update(kpc_per_pixel=kpc_per_pixel, role='target' if target_number is not None else 'control', image=image, variance=variance)
        return old_profile(image,variance,base,segmentation,detections,fallback_tree,
                           center_xy,pixel_scale,kpc_per_pixel,conversion,target_number,
                           coverage_floor,kron_factor)
    def instrument_mask(base,segmentation,detections,fallback_tree,center_xy,x0,y0,
                        target_number,pixel_scale,half_size,kron_factor):
        full=old_mask(base,segmentation,detections,fallback_tree,center_xy,x0,y0,
                      target_number,pixel_scale,half_size,kron_factor)
        shape=base.shape
        if shape not in zeros:
            zeros.clear()
            zeros[shape]=(np.zeros(shape,bool),np.zeros(shape,segmentation.dtype))
        zero_base,zero_seg=zeros[shape]
        empty=core.detection_arrays([],pixel_scale)
        h,w=full.shape
        component_base=base[y0:y0+h,x0:x0+w]
        labels=segmentation[y0:y0+h,x0:x0+w]
        component_seg=(labels>0) if target_number is None else ((labels>0)&(labels!=target_number))
        component_source=old_mask(zero_base,zero_seg,detections,None,center_xy,x0,y0,
                                  target_number,pixel_scale,half_size,kron_factor)
        component_fallback=old_mask(zero_base,zero_seg,empty,fallback_tree,center_xy,x0,y0,
                                    target_number,pixel_scale,half_size,kron_factor)
        comps=[component_base,component_seg,component_source,component_fallback]
        if not np.array_equal(full,np.logical_or.reduce(comps)):
            raise ValueError('mask component union differs from original mask')
        yy,xx=np.indices(full.shape,dtype=np.float32)
        radius=np.hypot(xx+x0-center_xy[0],yy+y0-center_xy[1])*state['kpc_per_pixel']
        ring=(radius>=50)&(radius<75)
        expected=math.pi*(75**2-50**2)/state['kpc_per_pixel']**2
        cx,cy=map(float,center_xy)
        ux0,ux1=int(math.floor(cx-half_size)),int(math.ceil(cx+half_size))+1
        uy0,uy1=int(math.floor(cy-half_size)),int(math.ceil(cy+half_size))+1
        uy,ux=np.indices((uy1-uy0,ux1-ux0),dtype=np.float32)
        ur=np.hypot(ux+ux0-cx,uy+uy0-cy)*state['kpc_per_pixel']
        full_ring=(ur>=50)&(ur<75)
        inside=(ux+ux0>=0)&(ux+ux0<shape[1])&(uy+uy0>=0)&(uy+uy0<shape[0])
        if int((full_ring&inside).sum())!=int(ring.sum()):
            raise ValueError('discrete in-image annulus does not reconcile')
        bitstamp=state['bits'][y0:y0+h,x0:x0+w]
        img=state['image'][y0:y0+h,x0:x0+w]
        var=state['variance'][y0:y0+h,x0:x0+w]
        quality=(bitstamp & sum(1<<b for b in BAD_BIT_POSITIONS))!=0
        bright=(bitstamp & BRIGHT_OBJECT_BIT)!=0
        bad_image=~np.isfinite(img)
        bad_variance=(~np.isfinite(var)) | (var<=0)
        base_parts=[quality,bright,bad_image,bad_variance]
        if not np.array_equal(component_base,np.logical_or.reduce(base_parts)):
            raise ValueError('base-mask causes do not reconstruct original base')
        row={'target_id':state['target_id'],'filter':state['filter'],'role':state['role'],
             'geometric_expected_pixels':expected,'in_image_ring_pixels':int(ring.sum()),
             'full_discrete_annulus_pixels':int(full_ring.sum()),
             'outside_image_annulus_pixels':int((full_ring&~inside).sum()),
             'unmasked_pixels':int((ring&~full).sum()),'masked_union_pixels':int((ring&full).sum())}
        for i,name in enumerate(('quality_invalid','segmentation','neighbor_shapes','unmatched_circles')):
            others=np.logical_or.reduce([c for j,c in enumerate(comps) if j!=i])
            row[name+'_pixels']=int((ring&comps[i]).sum())
            row[name+'_exclusive_pixels']=int((ring&comps[i]&~others).sum())
        for i,name in enumerate(('quality_bits','bright_object','nonfinite_image','invalid_variance')):
            others=np.logical_or.reduce([c for j,c in enumerate(base_parts) if j!=i])
            row[name+'_pixels']=int((ring&base_parts[i]).sum())
            row[name+'_base_exclusive_pixels']=int((ring&base_parts[i]&~others).sum())
        for bit in BAD_BIT_POSITIONS:
            name=state['bit_names'][bit]
            row[name+'_pixels']=int((ring&((bitstamp & (1<<bit))!=0)).sum())
        rows.append(row)
        return full
    core.make_stamp_mask=instrument_mask;core.profile_at=instrument_profile
    results=[]
    try:
        for item in selected:
            state.update(target_id=item['target_id'],filter=item['filter'])
            review=pipeline.validate_fits(Path(item['record']['image']),item['frame'])
            state['bit_names']={int(v):k for k,v in review['mask_planes'].items()}
            with core.fits.open(item['record']['image'],memmap=False,checksum=True) as hdus:
                state['bits']=np.asarray(hdus[review['mask_hdu']].data,dtype=np.int32)
            print('processing',item['filter'],item['target_id'],flush=True)
            profiles,diagnostics,frame=pipeline.process_frame(item,cosmo,.6,3,.6)
            results.extend(profiles)
    finally:
        core.make_stamp_mask=old_mask;core.profile_at=old_profile
    keys={(x['target_id'],x['filter']) for x in selected}
    old=[r for r in pipeline.read_tsv(MEMBER_BINS) if (r['target_id'],r['filter']) in keys]
    canonical=lambda r:(str(r['target_id']),str(r['filter']),str(r['member_index']),str(r['radial_bin_kpc']))
    oldmap={canonical(r):r for r in old};newmap={canonical(r):r for r in results}
    if set(oldmap)!=set(newmap): raise ValueError('retained keys differ from saved result')
    for key,r in newmap.items():
        for column in ('target_sb_njy_arcsec2','control_sb_njy_arcsec2','delta_sb_njy_arcsec2','target_coverage','control_coverage_mean'):
            if not np.isclose(float(r[column]),float(oldmap[key][column]),rtol=0,atol=1e-9):
                raise ValueError('retained photometry differs from saved result')
    if not rows: raise ValueError('no measured stamps')
    out.mkdir(parents=True,exist_ok=False)
    pipeline.write_tsv(out/'stamp_component_counts.tsv',rows,list(rows[0]))
    summary={'status':'complete_exploratory_spatial_mask_sample','science_measurement':False,
        'network_used':False,'selection_plan_sha256':pipeline.sha256(plan_path),'n_frames':len(selected),'frame_ids':sorted(list(keys)),
        'n_profile_stamps':len(rows),'retained_profile_rows_reproduced':len(results),
        'mask_union_exact':True,'base_mask_union_exact':True,'discrete_annulus_reconciled':True,'analysis_script_sha256':pipeline.sha256(Path(__file__)),
        'input_profile_sha256':pipeline.sha256(MEMBER_BINS),
        'output_sha256':pipeline.sha256(out/'stamp_component_counts.tsv'),
        'limitations':['Spatially balanced selection is not a population-proportional random sample.',
            'Components overlap: footprints cannot be added, and exclusive areas are not causal impacts.',
            'Control positions may contain cluster light; science gates remain unresolved.']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--selection-plan',type=Path,required=True)
    p.add_argument('--plan-sha256',required=True)
    a=p.parse_args()
    run(a.output.resolve(),a.selection_plan.resolve(),a.plan_sha256)
