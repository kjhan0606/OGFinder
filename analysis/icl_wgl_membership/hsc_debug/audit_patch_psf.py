#!/usr/bin/env python3
"""Independent projected-polygon, spherical-exclusion and native-pixel replay."""
import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import platform

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import astropy.units as u

from fetch_ancillary_catalogues import sha,catalogue_reader
from connect_control_patches import parse_patches,read_table
from download_pilot import estimate_bytes,BAD_BIT_POSITIONS,BRIGHT_OBJECT_BIT,build_manifest_row
from summarize_psf_gaia import classify
import prepare_detection_psf as prep


def projected_contains(poly,ra,dec):
    """Gnomonic projection + planar edge signs, independent of spherical edge normals."""
    lon,lat=np.deg2rad(poly['center_radec'])
    east=np.array([-np.sin(lon),np.cos(lon),0.])
    north=np.array([-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)])
    center=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    coords=poly['corners_radec']+[[float(ra),float(dec)]]
    v=SkyCoord(np.array(coords)[:,0]*u.deg,np.array(coords)[:,1]*u.deg).cartesian.xyz.value.T
    if np.any(v@center<=0):return False
    xy=np.column_stack((v@east,v@north))/(v@center)[:,None]
    vertices=xy[:4];point=xy[4];edge=np.roll(vertices,-1,axis=0)-vertices;relative=point-vertices
    cross=edge[:,0]*relative[:,1]-edge[:,1]*relative[:,0]
    return bool(np.all(cross>=-1e-14) or np.all(cross<=1e-14))


def controls(a,report):
    polygons={}
    for source in report['official_geometry']:
        path=Path(source['file']);assert sha(path)==source['sha256']
        polygons.update(parse_patches(path.read_text()))
    parent_path=next(Path(p) for p in report['source_sha256'] if p.endswith('/download_manifest.json'))
    parent=json.loads(parent_path.read_text());frames={str(f['id']):f for f in parent['clusters']}
    target=read_table(a.controls/'target_patch_associations.tsv')
    polygon_checks=0
    for row in target:
        f=frames[row['target_cluster_id']]
        assert row['selected_primary_tract_verified']=='False'
        for pid in filter(None,row['polygon_ids'].split(',')):
            assert projected_contains(polygons[pid],f['ra_deg'],f['dec_deg']);polygon_checks+=1
    with (a.sample/'catalogues_2026-10-13_01/patch_qa.csv').open() as f:qa={r['skymap_id']:r for r in catalogue_reader(f)}
    previous=json.loads((a.sample/'control_candidates_2026-10-14_01/summary.json').read_text())
    cat_path=next(Path(p) for p in previous['source_sha256'] if p.endswith('camira_s20a_wide_sm.tbl'))
    assert sha(cat_path)==previous['source_sha256'][str(cat_path)]
    cat=np.loadtxt(cat_path,usecols=(0,2,3,4));cp=SkyCoord(cat[:,1]*u.deg,cat[:,2]*u.deg)
    cosmo=FlatLambdaCDM(H0=parent['h0'],Om0=parent['omega_m'])
    radii=(1*u.Mpc/cosmo.angular_diameter_distance(cat[:,3])).value*180/np.pi
    rows=read_table(a.controls/'selected_controls.tsv');plan=json.loads((a.controls/'control_download_manifest.json').read_text())
    assert len(rows)==report['n_selected']==len(plan['clusters'])
    assert len({r['object_id'] for r in rows})==len(rows)==len({r['target_cluster_id'] for r in rows})
    pos=SkyCoord([float(r['ra']) for r in rows]*u.deg,[float(r['dec']) for r in rows]*u.deg)
    caps=np.array([float(r['footprint_cap_deg']) for r in rows]);min_cluster=180.;min_pair=180.
    for i,(row,transport) in enumerate(zip(rows,plan['clusters'])):
        cid=row['target_cluster_id'];f=frames[cid];q=qa[row['skymap_id']]
        assert projected_contains(polygons[row['skymap_id']],row['ra'],row['dec']);polygon_checks+=1
        assert row['center_flags_pass']==row['metadata_match']==row['random_patch_contains_point']=='True'
        for pid in row['target_patch_ids'].split(','):
            assert projected_contains(polygons[pid],f['ra_deg'],f['dec_deg']);polygon_checks+=1
            assert q['icalexp_exists']==qa[pid]['icalexp_exists']=='True'
            assert abs(float(q['imag_psf_depth'])-float(qa[pid]['imag_psf_depth']))<=.3
            assert 1/1.2<=float(q['iseeing'])/float(qa[pid]['iseeing'])<=1.2
        assert abs(float(row['ebv_sfd'])-float(row['target_ebv_sfd']))<=.02
        gap=pos[i].separation(cp).deg-caps[i]-radii
        assert np.all(gap>0);min_cluster=min(min_cluster,float(gap.min()))
        if i:
            gap=pos[i].separation(pos[:i]).deg-caps[i]-caps[:i]
            assert np.all(gap>0);min_pair=min(min_pair,float(gap.min()))
        np.testing.assert_allclose(caps[i],np.sqrt(2)*f['requested_half_width_arcsec']/3600+1/3600,atol=1e-12)
        assert transport['frame_role']=='matched_control' and transport['bcg_policy']=='none'
        assert transport['target_cluster_id']==cid and transport['random_object_id']==row['object_id']
        assert transport['ra_deg']==float(row['ra']) and transport['dec_deg']==float(row['dec']) and transport['z']==f['z']
        assert not transport['science_measurement'] and not transport['product_reviewed'] and 'image' not in transport
        # Original downloader output path preserves purpose and calibration gates.
        review=dict(mask_hdu=2,science_hdu=1,variance_hdu=3,coadd_header_filter='i',input_unit='instrumental_flux/pixel',ab_zeropoint=27.)
        saved=build_manifest_row(transport,'mock.fits',review)
        assert saved['frame_role']=='matched_control' and saved['bcg_policy']=='none' and saved['calibration']['a_i_mag'] is None
    estimate,_=estimate_bytes(plan['clusters'])
    np.testing.assert_allclose(estimate/1024**3,report['estimated_download_gib'])
    candidates=read_table(a.controls/'candidate_patch_audit.tsv');missing=[]
    chosen={r['target_cluster_id'] for r in rows}
    for cid in sorted(set(frames)-chosen):
        subset=[r for r in candidates if r['target_cluster_id']==cid]
        missing.append(dict(target_cluster_id=cid,**{k:v for k,v in next(r for r in target if r['target_cluster_id']==cid).items() if k!='target_cluster_id'},
            n_metadata_matches=sum(r['metadata_match']=='True' for r in subset),
            n_exclusion_clear_matches=sum(r['metadata_match']=='True' and r['overlap_1.0']=='0' for r in subset)))
    return dict(n_controls=len(rows),n_camira=len(cat),n_projected_polygon_checks=polygon_checks,
        min_cluster_cap_gap_arcsec=min_cluster*3600,min_control_cap_gap_arcsec=min_pair*3600,
        unassigned_targets=missing,download_manifest_sha256=sha(a.controls/'control_download_manifest.json'))


def psf(a,report):
    stars=read_table(a.psf/'star_ledger.tsv');methods=read_table(a.psf/'measurements.tsv');profiles=read_table(a.psf/'radial_profiles.tsv')
    assert len(stars)==report['n_candidates'] and len(methods)==len(stars)*4==report['n_method_rows']
    assert len(profiles)==report['n_profile_rows']
    assert dict(Counter(r['status'] for r in methods))==report['method_status_counts']
    ed=a.sample/'psf_cds_2026-10-13_01/evidence';gaia=defaultdict(list)
    with (ed/'gaia_matches.csv').open() as f:
        for r in csv.DictReader(f):gaia[int(r['candidate_id'])].append(r)
    manifest=json.loads((a.sample/'validation_sample_2026-10-13_01/validation_manifest.json').read_text())
    frames={str(f['id']):f for f in manifest['clusters']}
    download={str(f['id']):f for f in json.loads((a.sample/'download_01/download_status.json').read_text())['completed']}
    replay=0;mask_checks=0;native_stamps=0;max_mean_difference=0.;failures=Counter()
    for star in stars:
        sid=star['candidate_id'];cid=star['cluster_id'];frame=frames[cid];number=int(star['number'])
        assert classify(gaia[int(sid)])[0]=='stellar_evidence_requires_epoch_and_psf_review'
        group=[r for r in methods if r['candidate_id']==sid]
        assert {(r['kron_factor'],r['background_r_min_pix']) for r in group}=={('2.5','120'),('2.5','160'),('3.5','120'),('3.5','160')}
        source=a.sample/'detection_psf_2026-10-13_01'/('cluster_'+cid)
        detections=read_table(source/'detections.tsv');d=detections[number-1]
        assert float(star['native_x'])==float(d['X_IMAGE'])-1 and float(star['native_y'])==float(d['Y_IMAGE'])-1
        x0,y0=int(star['x0']),int(star['y0'])
        with fits.open(frame['image'],memmap=False) as raw:
            if not star['stamp']:
                edge=x0<0 or y0<0 or x0+401>raw[frame['science_hdu']].data.shape[1] or y0+401>raw[frame['science_hdu']].data.shape[0]
                assert (edge and all(r['status']=='blocked_stamp_edge' for r in group)) or star['gaia_variability']=='VARIABLE'
                failures.update(r['status'] for r in group);continue
            assert sha(Path(star['stamp']))==star['stamp_sha256']
            with fits.open(star['stamp'],memmap=False) as h,fits.open(source/'masks.fits.gz',memmap=False) as original:
                sl=np.s_[y0:y0+401,x0:x0+401];image=h['SCIENCE'].data;variance=h['VARIANCE'].data
                assert image.shape==(401,401)
                np.testing.assert_array_equal(image,raw[frame['science_hdu']].data[sl])
                np.testing.assert_array_equal(variance,raw[frame['variance']['hdu']].data[sl])
                bits=raw[download[cid]['review']['mask_hdu']].data[sl]
                base=(bits&(sum(1<<b for b in BAD_BIT_POSITIONS)|BRIGHT_OBJECT_BIT))!=0
                base|=~np.isfinite(image)|~np.isfinite(variance)|(variance<=0)
                np.testing.assert_array_equal(base,h['BASE_MASK'].data.astype(bool))
                np.testing.assert_array_equal(h['SEGMENTATION'].data,original['SEGMENTATION'].data[sl])
                masks={k:h['MASK_K'+k.replace('.','P')].data.astype(bool) for k in ('2.5','3.5')}
                assert all(np.all(m[base]) for m in masks.values()) and np.all(~masks['2.5']|masks['3.5'])
                # Independent rotated-ellipse containment at a fixed 9 x 9 native grid.
                other=[r for r in detections if int(r['NUMBER'])!=number]
                ox=np.array([float(r['X_IMAGE'])-1-x0 for r in other]);oy=np.array([float(r['Y_IMAGE'])-1-y0 for r in other])
                aa=np.maximum(.5,[float(r['A_IMAGE']) for r in other]);bb=np.maximum(.5,[float(r['B_IMAGE']) for r in other])
                theta=np.array([float(r['THETA_RAD']) for r in other]);kr=np.array([float(r['KRON_RADIUS']) for r in other])
                flags=np.array([int(r['KRON_FLAGS']) for r in other]);kr=np.where((flags==0)&np.isfinite(kr)&(kr>0),kr,1.2)
                for py in range(0,401,50):
                    for px in range(0,401,50):
                        dx,dy=px-ox,py-oy
                        rho=((dx*np.cos(theta)+dy*np.sin(theta))/aa)**2+((-dx*np.sin(theta)+dy*np.cos(theta))/bb)**2
                        label=h['SEGMENTATION'].data[py,px]
                        for factor,mask in masks.items():
                            expected=base[py,px] or (label>0 and label!=number) or np.any(rho<=(float(factor)*kr)**2)
                            assert bool(mask[py,px])==bool(expected);mask_checks+=1
                native_stamps+=1
                yy,xx=np.indices(image.shape);dx=xx-(float(star['native_x'])-x0);dy=yy-(float(star['native_y'])-y0)
                radius=np.sqrt(dx*dx+dy*dy);sector=np.minimum(7,np.floor((np.angle(dx+1j*dy)+np.pi)*4/np.pi).astype(int))
                for method in group:
                    good=~masks[method['kron_factor']]&np.isfinite(image)&np.isfinite(variance)&(variance>0)
                    ann=(radius>=int(method['background_r_min_pix']))&(radius<int(method['background_r_max_pix']))
                    sky=ann&good;counts=[np.count_nonzero(sky&(sector==k)) for k in range(8)]
                    fractions=[counts[k]/max(1,np.count_nonzero(ann&(sector==k))) for k in range(8)]
                    assert sky.sum()==int(method['sky_n_valid'])
                    if sky.sum()<300 or sky.sum()/ann.sum()<.3 or min(counts)<20 or min(fractions)<.25:
                        assert method['status']=='blocked_background_coverage';failures[method['status']]+=1;continue
                    background=np.mean([image[sky&(sector==k)].mean() for k in range(8)])
                    np.testing.assert_allclose(background,float(method['background']),atol=1e-12)
                    core=radius<12;valid=core&good;flux=np.sum(image[valid]-background)
                    if valid.sum()/core.sum()<.99:
                        assert method['status']=='blocked_core_mask';failures[method['status']]+=1;continue
                    assert method['status']=='diagnostic_only'
                    np.testing.assert_allclose(flux,float(method['core_flux']),rtol=1e-12)
                    subset=[r for r in profiles if r['candidate_id']==sid and r['kron_factor']==method['kron_factor'] and r['background_r_min_pix']==method['background_r_min_pix']]
                    assert len(subset)==11
                    for row in subset:
                        ring=(radius>=int(row['r_min_pix']))&(radius<int(row['r_max_pix']));valid=ring&good;n=valid.sum()
                        assert n==int(row['n_valid']);raw_mean=image[valid].mean()
                        difference=abs(raw_mean-float(row['raw_mean']));max_mean_difference=max(max_mean_difference,difference)
                        assert difference<1e-12
                        np.testing.assert_allclose((raw_mean-background)/flux,float(row['mean_per_core_flux']),rtol=1e-10,atol=1e-15)
                        np.testing.assert_allclose(np.sqrt(variance[valid].sum())/n,float(row['nominal_raw_mean_sigma']),rtol=1e-12)
                        replay+=1
    assert replay==len(profiles) and native_stamps==report['n_stamps']
    return dict(native_stamps_replayed=native_stamps,radial_bins_replayed=replay,ellipse_grid_checks=mask_checks,
        max_raw_mean_difference=max_mean_difference,blocked_method_counts=dict(failures))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sample','controls','psf','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    c=json.loads((a.controls/'summary.json').read_text());s=json.loads((a.psf/'summary.json').read_text())
    assert not c['science_measurement'] and not s['science_measurement'] and not s['psf_approved']
    hashes=dict(c['source_sha256'],**s['source_sha256'])
    for path,digest in hashes.items():assert sha(Path(path))==digest,path
    for name,report in (('connect_control_patches.py',c),('diagnose_empirical_psf.py',s)):
        assert sha(Path(__file__).parent/name)==report['script_sha256']
    cr=controls(a,c);print(json.dumps(dict(control_audit=cr)),flush=True)
    sr=psf(a,s)
    result=dict(passed=True,n_input_hashes=len(hashes),controls=cr,psf=sr,science_measurement=False,psf_approved=False,
        script_sha256=sha(Path(__file__)),report_sha256={str(d/'summary.json'):sha(d/'summary.json') for d in (a.controls,a.psf)})
    prep.write_json(a.work/'audit.json',result);print(json.dumps(result))


if __name__=='__main__':main()
