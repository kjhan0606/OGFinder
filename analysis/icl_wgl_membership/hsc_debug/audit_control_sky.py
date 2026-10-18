#!/usr/bin/env python3
"""Independent sky-coordinate exclusion audit and observed-tile numerical replay."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import platform

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from fetch_ancillary_catalogues import sha
import prepare_detection_psf as prep


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('controls', 'sky', 'work'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False)
    control = json.loads((a.controls/'summary.json').read_text())
    sky = json.loads((a.sky/'summary.json').read_text())
    assert not control['science_measurement'] and control['controls_ready_for_science'] == 0
    assert not sky['sky_correction_approved'] and not sky['observational_validation_passed']
    hashes = dict(control['source_sha256'], **sky['source_sha256'])
    for path, digest in hashes.items():
        assert sha(Path(path)) == digest, path
    for name, report in [('prepare_control_candidates.py', control), ('diagnose_outer_sky.py', sky)]:
        assert sha(Path(__file__).parent/name) == report['script_sha256']
    manifest_path = next(Path(x) for x in control['source_sha256'] if x.endswith('/download_manifest.json'))
    manifest = json.loads(manifest_path.read_text())
    cosmo = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'])
    frames = {str(r['id']):r for r in manifest['clusters']}
    cat_path = next(Path(x) for x in control['source_sha256'] if x.endswith('camira_s20a_wide_sm.tbl'))
    cat = np.loadtxt(cat_path, usecols=(0,2,3,4)); cc = SkyCoord(cat[:,1]*u.deg, cat[:,2]*u.deg)
    per_mpc = (1*u.Mpc/cosmo.angular_diameter_distance(cat[:,3])).value*180/np.pi
    candidates = table(a.controls/'candidate_audit.tsv'); plans = table(a.controls/'provisional_controls.tsv')
    assert len(candidates) == control['n_random_rows']
    assert len({(r['target_cluster_id'],r['object_id']) for r in candidates}) == len(candidates)
    lookup = {(r['target_cluster_id'],r['object_id']):r for r in candidates}
    for row in candidates:
        assert int(row['overlap_0.5']) <= int(row['overlap_1.0']) <= int(row['overlap_2.0'])
        if row['metadata_match'] == 'True':
            assert row['center_flags_pass'] == row['patch_qa_pass'] == 'True'
            assert abs(float(row['delta_depth'])) <= .3 and 1/1.2 <= float(row['seeing_ratio']) <= 1.2
            assert abs(float(row['delta_ebv'])) <= .02 and float(row['qa_proxy_distance_arcmin']) <= 5
    checks = []; allocation = {}
    for radius in (.5,1.,2.):
        selected = [r for r in plans if float(r['exclusion_mpc'])==radius]
        assert len({r['object_id'] for r in selected}) == len(selected)
        assert len({r['target_cluster_id'] for r in selected}) == len(selected)
        pos = SkyCoord([float(r['ra']) for r in selected]*u.deg, [float(r['dec']) for r in selected]*u.deg)
        caps = np.array([float(r['footprint_cap_deg']) for r in selected])
        pair_min_gap = 180.; cluster_min_gap = 180.
        for i,row in enumerate(selected):
            original = lookup[(row['target_cluster_id'],row['object_id'])]
            assert original['metadata_match'] == 'True' and int(original['overlap_'+str(radius)]) == 0
            assert float(row['ra'])==float(original['ra']) and float(row['dec'])==float(original['dec'])
            f = frames[row['target_cluster_id']]
            assert float(row['assigned_z']) == f['z']
            np.testing.assert_allclose(caps[i], np.sqrt(2)*f['requested_half_width_arcsec']/3600+1/3600, atol=1e-12)
            gap = pos[i].separation(cc).deg - caps[i] - radius*per_mpc
            assert np.all(gap > 0), (radius, row['target_cluster_id'])
            cluster_min_gap = min(cluster_min_gap, float(gap.min()))
            if i:
                gaps = pos[i].separation(pos[:i]).deg-caps[i]-caps[:i]
                assert np.all(gaps > 0)
                pair_min_gap = min(pair_min_gap, float(gaps.min()))
        expected = next(r for r in control['sensitivity'] if r['exclusion_mpc']==radius)
        assert len(selected) == expected['assigned_nonoverlap']
        checks.append(dict(exclusion_mpc=radius, independent_all_camira_checks=len(selected),
                           min_cluster_cap_gap_arcsec=cluster_min_gap*3600,
                           min_control_cap_gap_arcsec=pair_min_gap*3600))
        allocation[radius] = {r['target_cluster_id']:r['object_id'] for r in selected}
    changed = {str(radius):sum(allocation[radius].get(cid)!=oid for cid,oid in allocation[1.].items()) for radius in (.5,2.)}
    measured = table(a.sky/'frame_stats.tsv'); tiles = table(a.sky/'tile_stats.tsv'); lags = table(a.sky/'lag_products.tsv')
    assert len(measured)==sky['n_method_rows']==sky['n_frames']*3
    assert len(tiles)==sky['n_tiles'] and len(lags)==sky['n_lag_rows']==len(measured)*8
    counts = Counter((r['cluster_id'],r['kron_factor']) for r in tiles)
    cases = defaultdict(list)
    for row in measured:
        assert not row['science_measurement']=='True' and row['sky_correction_approved']=='False'
        assert int(row['n_tiles']) == counts[(row['cluster_id'],row['kron_factor'])]
        assert (row['status']=='diagnostic_only') == (int(row['n_tiles'])>=10)
        cases[row['cluster_id']].append(row)
    for cid, rows in cases.items():
        rows.sort(key=lambda r:float(r['kron_factor']))
        assert [float(r['kron_factor']) for r in rows] == [2.5,3.5,5.]
        assert int(rows[0]['n_outer_valid'])>=int(rows[1]['n_outer_valid'])>=int(rows[2]['n_outer_valid'])
    # Reconstruct raw masks and pixels for six previous injection geometries, all three masks.
    validation_path = next(Path(x) for x in sky['source_sha256'] if x.endswith('/validation_manifest.json'))
    validation = {str(f['id']):f for f in json.loads(validation_path.read_text())['clusters']}
    replays = []
    for cid in ('5724','619','5814','5104','3505','6940'):
        frame = validation[cid]
        mask_path = next(Path(x) for x in sky['source_sha256'] if x.endswith('/cluster_'+cid+'/proposed_masks.fits.gz'))
        with fits.open(frame['image'],memmap=False) as h, fits.open(mask_path,memmap=False) as mh:
            wcs = WCS(h[frame['science_hdu']].header).celestial
            cx,cy = wcs.all_world2pix(frame['ra_deg'],frame['dec_deg'],0)
            scale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
            kpp = cosmo.kpc_proper_per_arcmin(frame['z']).value*scale/60
            for factor in ('2.5','3.5','5.0'):
                subset = [r for r in tiles if r['cluster_id']==cid and r['kron_factor']==factor]
                if not subset:
                    continue
                for row in (subset[0],subset[len(subset)//2],subset[-1]):
                    x,y = int(row['tile_x']),int(row['tile_y']); sl = np.s_[y:y+64,x:x+64]
                    image = np.asarray(h[frame['science_hdu']].data[sl],dtype=float)
                    var = np.asarray(h[frame['variance']['hdu']].data[sl],dtype=float)
                    mask = mh['MASK_K'+factor.replace('.','P')].data[sl].astype(bool)
                    yy,xx = np.indices(image.shape); radius = np.hypot(xx+x-cx,yy+y-cy)*kpp
                    region = (radius>=600)&(radius<900)
                    valid = region & ~mask & np.isfinite(image)&np.isfinite(var)&(var>0)
                    n = int(valid.sum())
                    assert n==int(row['n_valid']) and int(region.sum())==int(row['n_region']), (cid,factor,x,y,n,row['n_valid'],int(region.sum()),row['n_region'])
                    diff = abs(float(image[valid].mean())-float(row['raw_mean']))
                    assert diff<1e-9
                    np.testing.assert_allclose(np.sqrt(var[valid].sum())/n,float(row['nominal_mean_sigma']),rtol=1e-10)
                    replays.append(dict(cluster_id=cid,kron_factor=factor,tile_x=x,tile_y=y,mean_absolute_difference=diff))
    prep.write_json(a.work/'audit.json',dict(passed=True,input_hashes_verified=len(hashes),
        control_geometry=checks, changed_assignment_vs_1mpc=changed, n_independent_tile_replays=len(replays),
        max_tile_mean_difference=max(r['mean_absolute_difference'] for r in replays), tile_replays=replays,
        science_measurement=False, controls_ready_for_science=0, sky_correction_approved=False,
        source_reports=dict(controls=str(a.controls),sky=str(a.sky)),script_sha256=sha(Path(__file__))))
    print(json.dumps(dict(passed=True,geometry=checks,changed_assignment_vs_1mpc=changed,
                          independent_tile_replays=len(replays),input_hashes_verified=len(hashes))))


if __name__ == '__main__':
    main()
