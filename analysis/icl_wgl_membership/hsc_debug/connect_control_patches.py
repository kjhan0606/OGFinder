#!/usr/bin/env python3
"""Public HSC patch polygons -> coordinate matching -> control acquisition input.

No HSC credentials, image retrieval, science approval or active-manifest edits.
"""
import argparse
from collections import defaultdict, Counter
import csv
import json
from pathlib import Path
import platform
import re
import shutil
import time
import urllib.request

import numpy as np
from scipy.spatial import cKDTree

from fetch_ancillary_catalogues import sha, catalogue_reader
from prepare_control_candidates import xyz, angle, qa_ok, select_nonoverlap
from download_pilot import estimate_bytes
import prepare_detection_psf as prep

BASE = 'https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/tract_patches/info/'
FIELDS = ('AEGIS', 'hectomap', 'autumn', 'spring')
LINE = re.compile(r'Tract:\s*(\d+)\s+Patch:\s*(\d+),(\d+)\s+(Center|Corner[0-4])\s+\(RA, Dec\):\s*\(([^,]+),\s*([^\)]+)\)')


def parse_patches(text):
    records = defaultdict(dict)
    for line in text.splitlines():
        match = LINE.fullmatch(line.strip())
        if not match:
            continue
        tract, x, y, label, ra, dec = match.groups(); key = int(tract)*10000+int(x)*100+int(y)
        if label in records[key]:
            raise ValueError('Duplicate polygon vertex')
        records[key][label] = (float(ra), float(dec))
    result = {}
    for key, record in records.items():
        if set(record) != {'Center', 'Corner0', 'Corner1', 'Corner2', 'Corner3', 'Corner4'}:
            raise ValueError('Incomplete official polygon')
        if record['Corner0'] != record['Corner4']:
            raise ValueError('Polygon is not closed')
        coords = np.array([record['Corner'+str(i)] for i in range(4)])
        vertices = xyz(coords[:,0], coords[:,1]); center = xyz([record['Center'][0]], [record['Center'][1]])[0]
        normals = np.cross(vertices, np.roll(vertices, -1, axis=0))
        norms = np.linalg.norm(normals, axis=1)
        if np.any(norms < 1e-12):
            raise ValueError('Degenerate polygon')
        normals /= norms[:,None]; normals *= np.where(normals@center >= 0, 1., -1.)[:,None]
        result[str(key)] = dict(center=center, normals=normals, vertices=vertices,
                                center_radec=record['Center'], corners_radec=coords.tolist())
    if not result:
        raise ValueError('No official patch polygons found')
    return result


def contains(poly, point, tolerance=1e-12):
    point = np.asarray(point)
    return (point@poly['center'] > 0) & np.all(point@poly['normals'].T >= -tolerance, axis=-1)


class Geometry:
    def __init__(self, polygons):
        self.polygons = polygons; self.ids = list(polygons)
        self.centers = np.array([polygons[i]['center'] for i in self.ids])
        self.tree = cKDTree(self.centers)
        self.radius_deg = max(float(np.max(angle(p['center'],p['vertices']))) for p in polygons.values())+1e-5

    def locate(self, ra, dec):
        points = xyz(ra, dec)
        nearby = self.tree.query_ball_point(points, 2*np.sin(np.deg2rad(self.radius_deg)/2))
        return [[self.ids[j] for j in ids if contains(self.polygons[self.ids[j]], point)]
                for point, ids in zip(points, nearby)]


def read_table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def transport_manifest(selected, parent):
    frames = []
    for row in selected:
        frames.append(dict(id='control_'+row['target_cluster_id']+'_'+row['object_id'],
            frame_role='matched_control', target_cluster_id=row['target_cluster_id'], random_object_id=row['object_id'],
            ra_deg=row['ra'], dec_deg=row['dec'], z=row['assigned_z'],
            redshift_policy='angular scale assigned from paired target, not a physical control redshift',
            requested_half_width_arcsec=row['requested_half_width_arcsec'], members_radec_deg=[],
            bcg_policy='none', sky='global-sky', product='coadd/bg', product_reviewed=False,
            target_patch_ids=row['target_patch_ids'], random_patch_id=row['skymap_id'],
            center_geometry_verified=True, full_footprint_qa='pending_download',
            calibration=dict(a_i_mag=None, extinction_evidence=None), ebv_sfd=row['ebv_sfd'],
            science_measurement=False))
    return dict(schema='hsc-bcg-diffuse-sample-v1', purpose='matched_control_image_validation',
        frame_role='matched_control', bcg_policy='none', band='HSC-I', requested_product='coadd/bg',
        h0=parent['h0'], omega_m=parent['omega_m'], requested_outer_radius_kpc=parent['requested_outer_radius_kpc'],
        clusters=frames, science_measurement=False, controls_ready_for_science=False,
        note='clusters is the existing download_pilot transport field; every entry is a control. '
             'This acquisition manifest is not the active science sample. Post-download whole-image/PSF/sky checks remain.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('sample','prepared','work'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--geometry-cache',type=Path,help='Reuse the four previously retrieved public text files')
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False); (a.work/'public_geometry').mkdir(); started=time.time()
    polygons = {}; provenance = []
    for field in FIELDS:
        name = 'tracts_patches_W-'+field+'.txt'; url=BASE+name; path=a.work/'public_geometry'/name
        if a.geometry_cache:
            source=a.geometry_cache/name
            count=source.stat().st_size
            if count>64*1024**2: raise ValueError('Geometry file exceeds bound')
            with source.open('rb') as inp,path.open('xb') as out: shutil.copyfileobj(inp,out)
        else:
            with urllib.request.urlopen(url,timeout=90) as response, path.open('xb') as out:
                count=0
                for block in iter(lambda:response.read(1024*1024),b''):
                    count += len(block)
                    if count>64*1024**2: raise ValueError('Geometry file exceeds bound')
                    out.write(block)
        parsed = parse_patches(path.read_text())
        for key, value in parsed.items():
            if key in polygons and not np.allclose(value['vertices'], polygons[key]['vertices'], atol=1e-12, rtol=0):
                raise ValueError('Conflicting public geometry')
            polygons[key]=value
        provenance.append(dict(url=url,file=str(path),bytes=count,sha256=sha(path),n_patches=len(parsed)))
        print(field,len(parsed),'patch polygons',flush=True)
    geo = Geometry(polygons)
    mpath=a.sample/'download_01/download_manifest.json'; parent=json.loads(mpath.read_text())
    frames={str(f['id']):f for f in parent['clusters']}
    qpath=a.sample/'catalogues_2026-10-13_01/patch_qa.csv'
    previous=json.loads((a.prepared/'summary.json').read_text())
    if sha(qpath)!=previous['source_sha256'][str(qpath)] or sha(mpath)!=previous['source_sha256'][str(mpath)]:
        raise ValueError('Source metadata changed')
    with qpath.open(encoding='utf-8-sig',newline='') as f:
        qa={r['skymap_id']:r for r in catalogue_reader(f)}
    target_ids=list(frames)
    found=geo.locate([frames[i]['ra_deg'] for i in target_ids],[frames[i]['dec_deg'] for i in target_ids])
    target_records=[]; refs={}
    for cid,ids in zip(target_ids,found):
        good=[i for i in ids if qa_ok(qa.get(i))]; refs[cid]=good
        target_records.append(dict(target_cluster_id=cid, polygon_ids=','.join(sorted(ids)),
            available_i_patch_ids=','.join(sorted(good)),n_geometric_patches=len(ids),n_available_i_patches=len(good),
            status='unique_available_patch' if len(good)==1 else 'overlap_preserved' if good else 'no_available_patch',
            selected_primary_tract_verified=False))
    raw=read_table(a.prepared/'candidate_audit.tsv'); rows=[]; patch_checks=Counter(); offsets=[]
    for source in raw:
        row=dict(source)
        for key in ('ra','dec','assigned_z','requested_half_width_arcsec','footprint_cap_deg','ebv_sfd','target_ebv_sfd'):
            row[key]=float(row[key])
        for radius in (.5,1.,2.): row['overlap_'+str(radius)]=int(row['overlap_'+str(radius)])
        pid=row['skymap_id']; cid=row['target_cluster_id']; ids=refs[cid]
        geometric=pid in polygons and bool(contains(polygons[pid],xyz([row['ra']],[row['dec']])[0]))
        row.update(random_patch_contains_point=geometric,target_patch_ids=','.join(sorted(ids)),
                   metadata_match=False,match_score=None,delta_depth_max=None,seeing_log_ratio_max=None,
                   qa_proxy_object_id=None,qa_proxy_distance_arcmin=None,delta_depth=None,seeing_ratio=None,delta_ebv=None)
        if row['center_flags_pass']=='True': patch_checks['flag_clean_rows']+=1;patch_checks['flag_clean_geometry_consistent']+=int(geometric)
        if geometric and ids and row['center_flags_pass']=='True' and qa_ok(qa.get(pid)):
            q=qa[pid]; delta=[float(q['imag_psf_depth'])-float(qa[i]['imag_psf_depth']) for i in ids]
            seeing=[np.log(float(q['iseeing'])/float(qa[i]['iseeing'])) for i in ids]
            depth=max(abs(v) for v in delta); seeing=max(abs(v) for v in seeing); de=row['ebv_sfd']-row['target_ebv_sfd']
            row.update(delta_depth_max=depth,seeing_log_ratio_max=seeing,delta_ebv=de,
                       match_score=(depth/.3)**2+(seeing/np.log(1.2))**2+(de/.02)**2,
                       metadata_match=depth<=.3 and seeing<=np.log(1.2) and abs(de)<=.02)
        rows.append(row)
    missing_centers=0
    for pid in polygons.keys() & qa.keys():
        q=qa[pid]
        try: ra,dec=float(q['ra']),float(q['dec'])
        except (ValueError,TypeError): missing_centers+=1;continue
        if not np.isfinite([ra,dec]).all(): missing_centers+=1;continue
        offsets.append(float(angle(polygons[pid]['center'],xyz([ra],[dec])[0])*3600))
    chosen=select_nonoverlap(rows,1.)
    old=json.loads((a.prepared/'request_plan.json').read_text())['controls']; old_by={r['target_cluster_id']:r['object_id'] for r in old}
    now={r['target_cluster_id']:r['object_id'] for r in chosen}
    plan=transport_manifest(chosen,parent); estimate,_=estimate_bytes(plan['clusters'])
    prep.table(a.work/'target_patch_associations.tsv',target_records,list(target_records[0]))
    prep.table(a.work/'candidate_patch_audit.tsv',rows,list(dict.fromkeys(k for r in rows for k in r)))
    if chosen: prep.table(a.work/'selected_controls.tsv',chosen,list(dict.fromkeys(k for r in chosen for k in r)))
    prep.write_json(a.work/'control_download_manifest.json',plan)
    prep.write_json(a.work/'summary.json',dict(status='complete',official_geometry=provenance,
        script_sha256=sha(Path(__file__)),elapsed_seconds=time.time()-started,
        source_sha256={str(p):sha(p) for p in (mpath,qpath,a.prepared/'candidate_audit.tsv',a.prepared/'request_plan.json',a.prepared/'summary.json')},
        n_polygons=len(polygons),n_targets=len(frames),target_status_counts=dict(Counter(r['status'] for r in target_records)),
        random_geometry_checks=dict(patch_checks),qa_centers_compared=len(offsets),qa_centers_missing=missing_centers,
        qa_center_max_separation_arcsec=max(offsets) if offsets else None,
        n_selected=len(chosen),newly_covered_targets=sorted(set(now)-set(old_by)),lost_targets=sorted(set(old_by)-set(now)),
        changed_existing_assignments=sum(now.get(i)!=oid for i,oid in old_by.items()),
        estimated_download_gib=estimate/1024**3,science_measurement=False,controls_downloaded=0,controls_ready_for_science=0,
        limitations=['Public polygons locate catalogue patch cells; tract overlap is retained, not an inferred primary-tract choice.',
                     'Candidate quality must match every available target patch when geometry overlaps.',
                     'Patch centers and boundaries do not prove whole-cutout i-band coverage, usable area, PSF or effective passband.',
                     'Fixed 1 Mpc CAMIRA exclusion and greedy nonoverlap assignment retain earlier selection limits.',
                     'Acquisition input reuses the original download_pilot container and transfer code; no image request made.']))
    print(json.dumps(dict(selected=len(chosen),target_statuses=dict(Counter(r['status'] for r in target_records)),
                          geometry_checks=dict(patch_checks),estimated_gib=estimate/1024**3)))


if __name__=='__main__':main()
