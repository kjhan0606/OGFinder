#!/usr/bin/env python3
"""Offline provisional controls from saved official randoms; no science activation."""
import argparse
from collections import defaultdict, Counter
import json
from pathlib import Path
import platform
import time

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.spatial import cKDTree

from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_ancillary import query_sfd
import prepare_detection_psf as prep

RADII_MPC = (.5, 1., 2.)
CAMIRA_SHA = 'd929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa'
BAD_FLAGS = ('inputcount_flag', 'inputcount_flag_noinputs', 'inputcount_flag_badcentroid',
             'pixelflags_offimage', 'pixelflags_edge', 'pixelflags_bad',
             'pixelflags_interpolatedcenter', 'pixelflags_saturatedcenter',
             'pixelflags_crcenter', 'pixelflags_suspectcenter', 'pixelflags_bright_objectcenter')


def xyz(ra, dec):
    ra = np.deg2rad(np.asarray(ra)); dec = np.deg2rad(np.asarray(dec))
    return np.column_stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])


def angle(a, b):
    """Stable great-circle degrees, including RA wrap and small separations."""
    return np.rad2deg(2*np.arcsin(np.clip(np.linalg.norm(a-b, axis=-1)/2, 0, 1)))


def number(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def center_ok(row):
    # Unknown/blank flags fail closed; SQL's bright/edge filter alone is insufficient.
    return (all(row.get(k) == 'True' for k in ('isprimary', 'detect_ispatchinner', 'detect_istractinner'))
            and all(row.get('i_'+k) == 'False' for k in BAD_FLAGS)
            and (number(row.get('i_inputcount_value')) or 0) > 0
            and (number(row.get('i_pix_variance')) or 0) > 0)


def qa_ok(row):
    return (row is not None and row.get('icalexp_exists') == 'True'
            and (number(row.get('imag_psf_depth')) or 0) > 0
            and (number(row.get('iseeing')) or 0) > 0)


def exclusion_counts(points, footprint_deg, clusters, radius_per_mpc_deg, radius_mpc):
    """Conservative spherical cap encloses whole square cutout, not centre-only veto."""
    tree = cKDTree(clusters)
    query_radius = np.deg2rad(footprint_deg + radius_mpc*np.max(radius_per_mpc_deg))
    nearby = tree.query_ball_point(points, 2*np.sin(query_radius/2))
    count = np.zeros(len(points), int)
    for i, ids in enumerate(nearby):
        if ids:
            sep = angle(points[i], clusters[ids])
            count[i] = np.count_nonzero(sep <= footprint_deg[i]+radius_mpc*radius_per_mpc_deg[ids])
    return count


def select_nonoverlap(rows, radius_mpc):
    """Greedy scarce-target-first assignment, reproducible; not a survey random sample."""
    groups = defaultdict(list)
    for row in rows:
        if row['metadata_match'] and row['overlap_'+str(radius_mpc)] == 0:
            groups[row['target_cluster_id']].append(row)
    chosen = []; used = set(); positions = []; radii = []
    for cid in sorted(groups, key=lambda c: (len(groups[c]), c)):
        for row in sorted(groups[cid], key=lambda r: (r['match_score'], r['object_id'])):
            if row['object_id'] in used:
                continue
            p = xyz([row['ra']], [row['dec']])[0]
            if positions and np.any(angle(p, np.asarray(positions)) <= row['footprint_cap_deg']+np.asarray(radii)):
                continue
            chosen.append(row); positions.append(p); radii.append(row['footprint_cap_deg']); used.add(row['object_id'])
            break
    return chosen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', type=Path, required=True); p.add_argument('--work', type=Path, required=True)
    p.add_argument('--catalog', type=Path, default=Path('/gpfs/kjhan/Astrafek/data/catalog/camira_s20a_wide_sm.tbl'))
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True, exist_ok=False); started = time.time(); inputs = {}
    def read_json(path):
        inputs[str(path)] = sha(path)
        return json.loads(path.read_text())
    manifest = read_json(a.sample/'download_01/download_manifest.json')
    frames = {str(f['id']): f for f in manifest['clusters']}
    cats = a.sample/'catalogues_2026-10-13_01'
    state = read_json(cats/'catalogue_status.json')
    if state['status'] != 'complete':
        raise ValueError('Official catalogue retrieval not complete')
    def read_csv(name):
        path = cats/(name+'.csv'); digest = sha(path); item = state['jobs'][name]
        if item['status'] != 'complete' or digest != item['sha256']:
            raise ValueError('Unverified catalogue '+name)
        inputs[str(path)] = digest
        with path.open(encoding='utf-8-sig', newline='') as f:
            rows = list(catalogue_reader(f))
        if len(rows) != item['rows']:
            raise ValueError('Catalogue row count changed')
        return rows
    qa_rows = read_csv('patch_qa'); qa = {r['skymap_id']: r for r in qa_rows}
    if len(qa) != len(qa_rows):
        raise ValueError('Ambiguous patch QA key')
    raw = []
    for name in sorted(state['jobs']):
        if name.startswith('random_'):
            for row in read_csv(name):
                row['_source'] = name; raw.append(row)
    if set(r['target_cluster_id'] for r in raw) != set(frames):
        raise ValueError('Random/target IDs differ')
    unique = {}; groups = defaultdict(list); records = []
    for r in raw:
        cid = r['target_cluster_id']; ident = r['object_id']
        ra, dec = number(r['ra']), number(r['dec'])
        if ra is None or dec is None or not 0 <= ra < 360 or not -90 <= dec <= 90:
            raise ValueError('Invalid coordinates')
        value = (ra, dec, r['skymap_id'])
        if ident in unique and unique[ident] != value:
            raise ValueError('Duplicate random ID with inconsistent coordinates or patch')
        unique[ident] = value
        q = qa.get(r['skymap_id']); f = frames[cid]
        separation = float(angle(xyz([ra], [dec])[0], xyz([f['ra_deg']], [f['dec_deg']])[0])*60)
        if separation > 30.00001:
            raise ValueError('Random outside prepared query cone')
        row = dict(target_cluster_id=cid, object_id=ident, ra=ra, dec=dec, skymap_id=r['skymap_id'],
                   tract=r['tract'], source_batch=r['_source'], separation_arcmin=separation,
                   center_flags_pass=center_ok(r), patch_qa_pass=qa_ok(q),
                   depth=number(q.get('imag_psf_depth')) if q else None,
                   seeing_native=number(q.get('iseeing')) if q else None,
                   requested_half_width_arcsec=f['requested_half_width_arcsec'], assigned_z=f['z'],
                   footprint_cap_deg=np.sqrt(2)*f['requested_half_width_arcsec']/3600 + 1/3600)
        records.append(row); groups[cid].append(row)
    del raw
    # SFD is a matching covariate; no passband-dependent flux correction is activated.
    provenance = read_json(a.sample/'ancillary_2026-10-13_01/preparation_status.json')
    maps = {}
    for pole, item in provenance['maps'].items():
        path = Path(item['file']); digest = sha(path)
        if digest != item['sha256']:
            raise ValueError('SFD map changed')
        inputs[str(path)] = digest; maps[pole] = path
    ids = list(unique)
    ebv = query_sfd(SkyCoord([unique[i][0] for i in ids]*u.deg, [unique[i][1] for i in ids]*u.deg), maps)
    ebv_by_id = dict(zip(ids, ebv))
    target_ids = list(frames)
    target_ebv = query_sfd(SkyCoord([frames[i]['ra_deg'] for i in target_ids]*u.deg,
                                  [frames[i]['dec_deg'] for i in target_ids]*u.deg), maps)
    target_ebv = dict(zip(target_ids, target_ebv)); references = []
    for cid, group in groups.items():
        local = [r for r in group if r['center_flags_pass'] and r['patch_qa_pass'] and r['separation_arcmin'] <= 5]
        ref = min(local, key=lambda r: (r['separation_arcmin'], r['object_id'])) if local else None
        references.append(dict(target_cluster_id=cid, proxy_object_id=ref['object_id'] if ref else None,
                               proxy_distance_arcmin=ref['separation_arcmin'] if ref else None,
                               proxy_skymap_id=ref['skymap_id'] if ref else None,
                               target_patch_identity_verified=False, ebv_sfd=target_ebv[cid]))
        for row in group:
            row.update(ebv_sfd=ebv_by_id[row['object_id']], target_ebv_sfd=target_ebv[cid],
                       qa_proxy_object_id=ref['object_id'] if ref else None,
                       qa_proxy_distance_arcmin=ref['separation_arcmin'] if ref else None,
                       delta_depth=None, seeing_ratio=None, delta_ebv=None, match_score=None, metadata_match=False)
            if ref and row['center_flags_pass'] and row['patch_qa_pass']:
                dd = row['depth']-ref['depth']; ratio = row['seeing_native']/ref['seeing_native']
                de = row['ebv_sfd']-target_ebv[cid]
                row.update(delta_depth=dd, seeing_ratio=ratio, delta_ebv=de,
                           match_score=(dd/.3)**2+(np.log(ratio)/np.log(1.2))**2+(de/.02)**2,
                           metadata_match=abs(dd) <= .3 and 1/1.2 <= ratio <= 1.2 and abs(de) <= .02)
    if sha(a.catalog) != CAMIRA_SHA:
        raise ValueError('CAMIRA source hash changed')
    inputs[str(a.catalog)] = CAMIRA_SHA
    cat = np.loadtxt(a.catalog, usecols=(0, 2, 3, 4))
    if not np.all(np.isfinite(cat)) or np.any(cat[:, 3] <= 0) or len(np.unique(cat[:, 0])) != len(cat):
        raise ValueError('Invalid full CAMIRA catalogue')
    cosmo = FlatLambdaCDM(H0=manifest['h0'], Om0=manifest['omega_m'])
    one_mpc = np.rad2deg(1/cosmo.angular_diameter_distance(cat[:, 3]).value)
    clusters = xyz(cat[:, 1], cat[:, 2]); points = xyz([r['ra'] for r in records], [r['dec'] for r in records])
    caps = np.array([r['footprint_cap_deg'] for r in records]); sensitivity = []; plans = []
    for radius in RADII_MPC:
        counts = exclusion_counts(points, caps, clusters, one_mpc, radius)
        for row, n in zip(records, counts):
            row['overlap_'+str(radius)] = int(n)
        chosen = select_nonoverlap(records, radius)
        eligible = [r for r in records if r['metadata_match'] and r['overlap_'+str(radius)] == 0]
        sensitivity.append(dict(exclusion_mpc=radius, candidate_rows=len(eligible),
                                targets_with_candidates=len({r['target_cluster_id'] for r in eligible}),
                                assigned_nonoverlap=len(chosen)))
        plans.extend(dict(r, exclusion_mpc=radius, science_ready=False, full_footprint_qa='pending_image',
                          target_patch_match='local_proxy_only', status='provisional_geometry_queue') for r in chosen)
        print('CAMIRA exclusion', radius, 'Mpc:', sensitivity[-1], flush=True)
    prep.table(a.work/'candidate_audit.tsv', records, list(records[0]))
    prep.table(a.work/'target_qa_proxies.tsv', references, list(references[0]))
    if plans:
        prep.table(a.work/'provisional_controls.tsv', plans, list(dict.fromkeys(k for r in plans for k in r)))
    default = [r for r in plans if r['exclusion_mpc'] == 1.]
    prep.write_json(a.work/'request_plan.json', dict(status='provisional_not_science_ready',
        controls=default, requested_product='coadd/bg', requested_band='HSC-I', controls_downloaded=0,
        science_measurement=False, authentication_requested=False,
        pending=['Exact target patch identification and metadata matching', 'Full image footprint/filter/mask/variance QA',
                 'Per-control extinction/passband validation', 'Empirical sky covariance and PSF matching',
                 'Region/environment selection bias and missing-cluster catalogue coverage']))
    prep.write_json(a.work/'summary.json', dict(status='complete', source_sha256=inputs,
        script_sha256=sha(Path(__file__)), elapsed_seconds=time.time()-started,
        n_targets=len(frames), n_random_rows=len(records), n_unique_random_ids=len(unique), n_camira=len(cat),
        n_center_flags_pass=sum(r['center_flags_pass'] for r in records),
        n_patch_qa_pass=sum(r['patch_qa_pass'] for r in records),
        n_metadata_match=sum(r['metadata_match'] for r in records),
        n_targets_missing_local_proxy=sum(r['proxy_object_id'] is None for r in references), sensitivity=sensitivity,
        settings=dict(depth_tolerance_mag=.3, seeing_ratio_max=1.2, ebv_tolerance=.02, qa_proxy_max_arcmin=5,
                      footprint='sqrt(2)*requested_half_width + 1 arcsec conservative cap', default_exclusion_mpc=1.),
        science_measurement=False, controls_ready_for_science=0, controls_downloaded=0,
        limitations=['Thresholds and radii are declared engineering choices, not fitted science acceptance criteria.',
                     'Target patch is a nearest valid random within five arcminutes proxy, not a verified WCS patch assignment.',
                     'All entries of the supplied bright-star-masked S20A CAMIRA catalogue at all redshifts are excluded; catalogue is not all clusters in the Universe.',
                     'Photometric cluster redshift sets angular exclusion; no R200 or mass claim.',
                     'Conservative cap can reject disjoint square footprints; selection biases and spatial covariance remain.',
                     'Greedy one-per-target assignment maximizes neither global completeness nor independent regions.',
                     'Original hash-ordered truncated random queries are not an area-unbiased random sample.',
                     'Random sky_mean/std are patch-level two-arcsecond aperture statistics, not local coadd/bg sky estimates.'],
        official_sources=['https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/',
                          'https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/']))


if __name__ == '__main__':
    main()
