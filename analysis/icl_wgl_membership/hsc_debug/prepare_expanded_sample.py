#!/usr/bin/env python3
"""Prepare an 800-cluster extension, preserving existing selection and files."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import random
import sys
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from download_pilot import estimate_bytes, file_sha256, validate_fits
from hsc_pilot import template_frame
from ogfkit.wlcluster import ARCSEC, Cosmology


def write_json(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def strata(row):
    zbin = int(math.floor((row[4] - .30 + 1e-8) / .05))
    richness = next((i for i, upper in enumerate((20, 30, 50)) if row[5] < upper), 3)
    return zbin, richness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', default='/gpfs/kjhan/Astrafek/data')
    parser.add_argument('--previous', default='sample_2026-10-12')
    parser.add_argument('--work', required=True)
    parser.add_argument('--target', type=int, default=800)
    parser.add_argument('--seed', type=int, default=20261013)
    args = parser.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        parser.error('prepare and rehash existing files on lageunha')
    data, work = Path(args.data).resolve(), Path(args.work).resolve()
    if work.exists():
        raise ValueError('refusing to overwrite existing output directory')
    spec = importlib.util.spec_from_file_location('original_selection', HERE / '02_select_clusters.py')
    selector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(selector)
    catalog_path = data / 'catalog/camira_s20a_wide_sm.tbl'
    catalog = selector.read_catalog(catalog_path)
    by_id = {str(row[0]): row for row in catalog}
    previous = data / args.previous
    old = json.loads((previous / 'sample_manifest.json').read_text())
    old_status = json.loads((previous / 'download_01/download_status.json').read_text())
    old_frames = {str(f['id']): f for f in old['clusters']}
    old_records = {str(r['id']): r for r in old_status['completed']}
    if old_status['status'] != 'complete' or set(old_records) != set(old_frames):
        raise ValueError('previous download is not complete')
    if not len(old_frames) <= args.target <= 1000:
        raise ValueError('target must preserve existing sample and be at most 1000')
    chosen = [by_id[i] for i in sorted(old_frames)]
    eligible = [row for row in catalog if .30 <= row[4] < .55 and row[6] > 0
                and abs(row[6] - row[4]) / (1 + row[4]) < .01
                and str(row[0]) not in old_frames]
    separated = []
    for row in sorted(eligible, key=lambda r: (-r[5], r[0])):
        if all(selector.angular_sep_deg(row, other) > .1 for other in chosen):
            chosen.append(row)
            separated.append(row)
    pending_count = args.target - len(old_frames)
    if len(separated) < pending_count:
        raise ValueError('insufficient separated candidates')
    groups = defaultdict(list)
    for row in separated:
        groups[strata(row)].append(row)
    exact = {key: pending_count * len(values) / len(separated) for key, values in groups.items()}
    allocation = {key: int(value) for key, value in exact.items()}
    remainder = pending_count - sum(allocation.values())
    for key in sorted(groups, key=lambda k: (-(exact[k] - allocation[k]), k))[:remainder]:
        allocation[key] += 1
    rng = random.Random(args.seed)
    added = []
    for key in sorted(groups):
        added.extend(rng.sample(sorted(groups[key], key=lambda r: r[0]), allocation[key]))
    selected = sorted([by_id[i] for i in old_frames] + added, key=lambda r: (r[4], r[0]))
    min_sep = min(selector.angular_sep_deg(a, b) for i, a in enumerate(selected) for b in selected[i+1:])
    if min_sep <= .1 or len(selected) != args.target:
        raise ValueError('expanded sample separation or count failed')
    frames = []
    cosmology = Cosmology(H0=old['h0'], Om=old['omega_m'])
    for row in selected:
        ident = str(row[0])
        if ident in old_frames:
            frame = dict(old_frames[ident])
        else:
            frame = template_frame({'id': ident, 'name': row[1], 'ra_deg': row[2],
                                    'dec_deg': row[3], 'z_spec_bcg': row[6], 'N_mem': row[5]})
            frame['requested_half_width_arcsec'] = 1000 / (float(cosmology.angular_diameter(frame['z'])) * ARCSEC * 1000)
        zbin, rbin = strata(row)
        frame.update(z_cl=row[4], z_spec_bcg=row[6],
                     redshift_bin='%0.2f-%0.2f' % (.30 + .05*zbin, .35 + .05*zbin),
                     richness_bin=('15-20', '20-30', '30-50', '50+')[rbin],
                     reused_from_previous=ident in old_frames)
        frames.append(frame)
    new_by_coord = {(round(f['ra_deg'], 5), round(f['dec_deg'], 5)): f
                    for f in frames if f['id'] not in old_frames}
    member_path = data / 'catalog/camira_s20a_wide_sm_member.dat.gz'
    member_sha = file_sha256(member_path)
    if member_sha != 'bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e':
        raise ValueError('member catalogue digest mismatch')
    with gzip.open(member_path, 'rt') as stream:
        for line in stream:
            if not line.strip() or line.startswith('#'):
                continue
            values = list(map(float, line.split()))
            if len(values) != 8:
                raise ValueError('unexpected member catalogue columns')
            frame = new_by_coord.get((round(values[0], 5), round(values[1], 5)))
            if frame is None:
                continue
            if abs(values[2] - frame['richness']) > .001 or abs(values[3] - frame['z_cl']) > .001:
                raise ValueError('member association mismatch for %s' % frame['id'])
            frame['members_radec_deg'].append(values[4:6])
    if any(not f['members_radec_deg'] for f in frames):
        raise ValueError('one or more clusters have no matched members')
    verified = []
    for ident, rec in old_records.items():
        path = (previous / 'download_01/images' / rec['file']).resolve()
        if path.stat().st_size != rec['bytes'] or file_sha256(path) != rec['sha256']:
            raise ValueError('existing file digest mismatch: %s' % ident)
        validate_fits(path, old_frames[ident])
        verified.append((ident, rec, path))
    work.mkdir(parents=True, exist_ok=False)
    table = work / 'selected_clusters.tsv'
    with table.open('x', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['id','name','ra_deg','dec_deg','z_cl','N_mem','z_spec_bcg','frac_dz'])
        for row in selected:
            writer.writerow([*row, abs(row[6]-row[4])/(1+row[4])])
    manifest = {key: value for key, value in old.items() if key not in ('clusters','parent_selection','parent_selection_sha256')}
    manifest.update(clusters=frames, parent_selection=str(table),
                    parent_selection_sha256=file_sha256(table),
                    selection='Existing 213 retained; additions proportionally sampled by z_cl x richness strata from richness-priority separated candidates',
                    selection_limits='0.30 <= z_cl < 0.55; BCG spectroscopy frac_dz < 0.01; >6 arcmin between all selected centers; not population complete',
                    selection_seed=args.seed, geometry_redshift='z_spec_bcg',
                    catalog_sha256=file_sha256(catalog_path), member_catalog_sha256=member_sha)
    manifest_path = work / 'sample_manifest.json'
    write_json(manifest_path, manifest)
    download = work / 'download_01'
    images = download / 'images'
    images.mkdir(parents=True)
    rec_by_id = {ident: rec for ident, rec, path in verified}
    for ident, rec, path in verified:
        (images / rec['file']).symlink_to(path)
    estimated, estimates = estimate_bytes(frames)
    added_ids = {str(row[0]) for row in added}
    estimated_new = sum(row['estimated_bytes'] for row in estimates if row['id'] in added_ids)
    state = {'status':'starting', 'n_expected':len(frames), 'downloaded':len(verified),
             'completed':[rec_by_id[f['id']] for f in frames if f['id'] in rec_by_id],
             'estimated_bytes':estimated, 'manifest_sha256':file_sha256(manifest_path),
             'expected_ids':[f['id'] for f in frames], 'science_measurement':False,
             'credentials_saved':False, 'bulk_download':False,
             'reused_from':str(previous / 'download_01')}
    write_json(download / 'download_status.json', state)
    write_json(work / 'random_controls_plan.json', {
        'status':'pending_survey_random_catalog_and_observing_metadata',
        'preliminary_target_controls':len(frames), 'download_started':False,
        'matching':['survey footprint and usable area','sky region','HSC-I depth and PSF','same angular cutout and assigned cluster redshift','identical masking and calibration'],
        'note':'No arbitrary coordinate offsets are certified as matched controls. Final control count follows measured background covariance.'})
    summary = {'status':'prepared', 'n_total':len(frames), 'n_reused':len(verified),
               'n_new_downloads':len(added), 'n_separated_pool':len(separated)+len(old_frames),
               'min_separation_arcmin':60*min_sep,
               'redshift_bin_counts':dict(Counter(f['redshift_bin'] for f in frames)),
               'richness_bin_counts':dict(Counter(f['richness_bin'] for f in frames)),
               'estimated_total_gib':estimated/1024**3,
               'estimated_additional_gib':estimated_new/1024**3,
               'n_member_rows':sum(len(f['members_radec_deg']) for f in frames),
               'download_started':False, 'science_measurement':False,
               'manifest':str(manifest_path)}
    write_json(work / 'preparation.json', summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
