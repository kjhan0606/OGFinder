#!/usr/bin/env python3
"""Radial/azimuthal area and deblend sensitivity; no flux inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree
import prepare_detection_psf as prep
from refine_bcg_member_masks import satellite_members
from review_bcg_association import peak_candidate


def coverage(mask, radius, angle, edges):
    rows=[]
    for low,high in zip(edges[:-1],edges[1:]):
        ring=(radius>=low)&(radius<high);total=int(ring.sum());valid=int(np.count_nonzero(ring&~mask))
        sectors=[]
        for sector in range(8):
            region=ring&(angle>=sector*np.pi/4)&(angle<(sector+1)*np.pi/4)
            n=int(region.sum());sectors.append(float(np.count_nonzero(region&~mask)/n) if n else None)
        rows.append(dict(r_min_kpc=low,r_max_kpc=high,n_grid=total,n_valid=valid,
            grid_coverage=valid/total if total else None,
            min_sector_coverage=min(x for x in sectors if x is not None) if any(x is not None for x in sectors) else None,
            sectors_below_half=sum(x is not None and x<.5 for x in sectors)))
    return rows


def center_causes(objects,seg,kron,kflags,base,bcg,members,center,scale,factor):
    cx,cy=center;ix,iy=round(float(cx)),round(float(cy))
    dx=cx-objects['x'];dy=cy-objects['y'];ct=np.cos(objects['theta']);st=np.sin(objects['theta'])
    size=factor*np.where((kflags==0)&np.isfinite(kron)&(kron>0),kron,1.2)
    q=((dx*ct+dy*st)/np.maximum(objects['a'],.5))**2+((-dx*st+dy*ct)/np.maximum(objects['b'],.5))**2
    ids=np.flatnonzero((q<=size**2)&(np.arange(len(objects))!=bcg))
    dist=cKDTree(np.column_stack([objects['x'],objects['y']])).query(members)[0] if len(objects) and len(members) else np.full(len(members),np.inf)
    fallback=members[dist>max(1/scale,1.5)]
    fallnear=np.hypot(fallback[:,0]-cx,fallback[:,1]-cy)<=max(2/scale,1.5)
    return dict(base_center=bool(base[iy,ix]),source_segment_center=bool(seg[iy,ix]>0 and seg[iy,ix]!=bcg+1),
        neighbor_ellipse_numbers=(ids+1).tolist(),fallback_members_covering_center=int(fallnear.sum()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','prepared','refined','work'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    doc=json.loads(a.manifest.read_text());ref=json.loads((a.refined/'summary.json').read_text())
    assert ref['manifest_sha256']==hashlib.sha256(a.manifest.read_bytes()).hexdigest()
    lookup={r['id']:r for r in ref['records']};a.work.mkdir(parents=True,exist_ok=False)
    cosmo=prep.Cosmology(H0=doc['h0'],Om=doc['omega_m']);edges=np.arange(0,501,25).tolist()
    records=[];all_rows=[]
    for n,frame in enumerate(doc['clusters'],1):
        cid=str(frame['id']);r=lookup[cid];sub='cluster_'+cid;d=a.prepared/sub
        with fits.open(a.refined/sub/'proposed_masks.fits.gz') as h:
            base=h['BASE_MASK'].data.astype(bool);header=h[1].header.copy();original=[x.data.astype(bool) for x in h[2:]]
        with fits.open(d/'masks.fits.gz') as h:seg=h['SEGMENTATION'].data.copy()
        wcs=WCS(header).celestial;center=tuple(float(x) for x in wcs.all_world2pix(frame['ra_deg'],frame['dec_deg'],0));scale=3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
        with (d/'detections.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        objects=np.zeros(len(rows),dtype=[(k,'f8') for k in ('x','y','a','b','theta')])
        for k,col in [('x','X_IMAGE'),('y','Y_IMAGE'),('a','A_IMAGE'),('b','B_IMAGE'),('theta','THETA_RAD')]:objects[k]=[float(row[col])-(1 if k in ('x','y') else 0) for row in rows]
        kron=np.array([float(row['KRON_RADIUS']) for row in rows]);flags=np.array([int(row['KRON_FLAGS']) for row in rows])
        with (d/'members.tsv').open() as f:members=np.array([[float(x['X_IMAGE'])-1,float(x['Y_IMAGE'])-1] for x in csv.DictReader(f,delimiter='\t') if x['IN_FRAME']=='True']).reshape(-1,2)
        members,_=satellite_members(members,center,scale,r['stable_peak_candidate']);bcg=r['candidate_labels'][0]-1 if r['stable_peak_candidate'] else -1
        causes=center_causes(objects,seg,kron,flags,base,bcg,members,center,scale,2.5)
        info=prep.validate_fits(Path(frame['image']),frame)
        with fits.open(frame['image']) as h:
            image=np.array(h[info['science_hdu']].data,dtype='f4');var=np.array(h[info['variance_hdu']].data,dtype='f4');bits=h[info['mask_hdu']].data.copy()
        bad=((bits&sum(1<<b for b in prep.BAD_BIT_POSITIONS))!=0)|~np.isfinite(image)|~np.isfinite(var)|(var<=0)
        alt=prep.detect(image,var,bad,center,members,scale,deblend_cont=.001)
        ai,reason=peak_candidate(alt[0],alt[1],bad,center,scale)
        if not r['stable_peak_candidate']:ai=-1
        if r['stable_peak_candidate'] and ai+1!=r['candidate_labels'][1]:raise ValueError('Alternative candidate no longer reproduces')
        amasks,_=prep.source_masks(alt[0],alt[1],alt[2],alt[3],base,ai,members,scale)
        yy,xx=np.indices(base.shape,dtype='f4');dx=xx-float(center[0]);dy=yy-float(center[1])
        kpp=float(cosmo.angular_diameter(frame['z']))*prep.ARCSEC*1000*scale
        radius=np.hypot(dx,dy)*kpp;angle=np.mod(np.arctan2(dy,dx),2*np.pi)
        products=[('base',0,base)]+[('deblend_005',k,m) for k,m in zip(prep.FACTORS,original)]+[('deblend_001',k,m) for k,m in zip(prep.FACTORS,amasks)]
        local=[]
        for method,k,mask in products:
            for row in coverage(mask,radius,angle,edges):
                geometric=np.pi*(row['r_max_kpc']**2-row['r_min_kpc']**2)/kpp**2
                local.append(dict(cluster_id=cid,cohort=r['cohort'],method=method,kron_factor=k,
                    geometric_coverage=min(1.,row['n_valid']/geometric),**row))
        all_rows.extend(local)
        delta=max(abs(x['grid_coverage']-y['grid_coverage']) for x,y in zip(local[20:80],local[80:140]))
        record=dict(id=cid,cohort=r['cohort'],stable_peak_candidate=r['stable_peak_candidate'],center_causes_kron2p5=causes,
            max_radial_coverage_change_deblend=delta,
            n_alternative_detections=len(alt[0]),science_measurement=False,mask_ready_for_science=False)
        records.append(record);print(f'[{n}/50] {cid}: max area change {delta:.3f}',flush=True)
    prep.table(a.work/'radial_coverage.tsv',all_rows,list(all_rows[0]))
    prep.write_json(a.work/'summary.json',dict(status='complete',records=records,n_radial_rows=len(all_rows),science_measurement=False,
        mask_ready_for_science=False,sector_count=8,edges_kpc=edges,manifest_sha256=hashlib.sha256(a.manifest.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note='Area only; no flux comparison or mask-completeness/injection approval. Deblend comparison retains blocked BCG gates. Center ellipse attribution uses continuous catalogue position; pixel mask uses integer sample.'))


if __name__=='__main__':main()
