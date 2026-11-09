#!/usr/bin/env python3
"""Run a small source-aware CAMIRA-member outer-light profile pilot.

Uses only the already-downloaded S23B target images and imported member
positions. SEP segmentation is regenerated in memory so each matched member's
own detection can be preserved while neighboring sources and quality pixels
are masked. Rotated source-free positions at the same cluster-centric radius
provide a local baseline. This is a method/support pilot, not a science result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree

import measure_satellite_outer_profiles as core
import prepare_detection_psf as prep
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, validate_fits

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
TARGET_MANIFEST = DATA / 's23b_icl_target_image_plan_2026-11-07_01/target_manifest.json'
PAIR_CONTRASTS = DATA / 's23b_icl_relative_radial_contrasts_2026-11-07_02/paired_relative_contrasts.tsv'
MEMBER_TABLE = DATA / 's23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/target_member_positions.tsv'
MASK_RECORDS = DATA / 's23b_icl_mask_support_2026-11-07_01/target/records.tsv'
DEFAULT_OUTPUT = DATA / 's23b_camira_member_outer_profile_pilot_2026-11-07_01'
EXPECTED_TARGET_MANIFEST_SHA256 = 'e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44'
EXPECTED_PAIR_CONTRASTS_SHA256 = 'c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800'
EXPECTED_MEMBER_TABLE_SHA256 = 'd1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae'
EXPECTED_MASK_RECORDS_SHA256 = '609354bfa7777088e995a487d0c2eaaed717b80e94589f1c5e0c3758f3e32ff8'
MIN_CLUSTER_RADIUS_KPC = 50.0
MAX_CLUSTER_RADIUS_KPC = 500.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        if fields is None:
            raise ValueError(f'cannot infer columns for empty table {path}')
        columns = fields
    else:
        columns = fields or list(rows[0])
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value):
    return float(value) if np.isfinite(value) else None


def detection_rows(objects, kron, kron_flags):
    rows = []
    for i, obj in enumerate(objects):
        rows.append({
            'NUMBER': str(i + 1),
            'X_IMAGE': str(float(obj['x']) + 1.0),
            'Y_IMAGE': str(float(obj['y']) + 1.0),
            'A_IMAGE': str(float(obj['a'])),
            'B_IMAGE': str(float(obj['b'])),
            'THETA_RAD': str(float(obj['theta'])),
            'KRON_RADIUS': '' if not np.isfinite(kron[i]) else str(float(kron[i])),
            'KRON_FLAGS': str(int(kron_flags[i])),
        })
    return core.parse_detections(rows)


def build_member_rows(group: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    for index, row in enumerate(group):
        rows.append({
            'INDEX': index,
            'X_IMAGE': str(float(row['image_x_pixel']) + 1.0),
            'Y_IMAGE': str(float(row['image_y_pixel']) + 1.0),
            'IN_FRAME': str(row['inside_image_pixel_centers'] == 'True'),
            'source_row': row['source_row'],
            'member_ra_deg': row['member_ra_deg'],
            'member_dec_deg': row['member_dec_deg'],
            'member_field8_raw': row['member_field8_raw'],
        })
    return rows


def select_sample(pair_rows, member_rows, frame_by_id, record_by_id, cosmo,
                  n_per_filter: int, seed: int) -> tuple[list[dict], dict[str, list[dict]]]:
    members_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in member_rows:
        members_by_target[row['target_cluster_id']].append(row)

    candidate_info = {}
    eligible_by_filter: dict[str, list[dict]] = defaultdict(list)
    for row in pair_rows:
        if row['contrast'] != '50-300_minus_300-800':
            continue
        try:
            icl_contrast = float(row['contrast_njy_arcsec2'])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(icl_contrast):
            continue
        ident, band = row['target_id'], row['filter']
        frame = frame_by_id.get(ident)
        mask_record = record_by_id.get(ident)
        if frame is None or mask_record is None or band not in ('i', 'i2'):
            continue
        if mask_record['primary_filter'] != band:
            raise ValueError(f'profile/mask filter mismatch for {ident}')
        z = float(frame['z'])
        kpc_per_arcsec = float(cosmo.kpc_proper_per_arcmin(z).value / 60.0)
        group = members_by_target.get(ident, [])
        eligible = []
        for index, member in enumerate(group):
            if member['inside_image_pixel_centers'] != 'True':
                continue
            radius_kpc = float(member['cluster_offset_arcsec']) * kpc_per_arcsec
            if MIN_CLUSTER_RADIUS_KPC <= radius_kpc < MAX_CLUSTER_RADIUS_KPC:
                eligible.append((index, radius_kpc))
        if len(eligible) < 4:
            continue
        candidate_info[ident] = {
            'target_id': ident, 'filter': band, 'z': z,
            'richness': float(frame['richness']),
            'icl_contrast_njy_arcsec2': icl_contrast,
            'n_member_positions_50_500_kpc_in_frame': len(eligible),
            'members': group,
            'eligible_members': eligible,
            'frame': frame,
            'record': mask_record,
        }
        eligible_by_filter[band].append(candidate_info[ident])

    rng = np.random.default_rng(seed)
    selected = []
    for band in ('i', 'i2'):
        pool = sorted(eligible_by_filter[band], key=lambda item: item['target_id'])
        if len(pool) < n_per_filter:
            raise ValueError(f'{band}: only {len(pool)} eligible clusters for requested {n_per_filter}')
        indices = np.sort(rng.choice(len(pool), size=n_per_filter, replace=False))
        selected.extend(pool[int(index)] for index in indices)
    return sorted(selected, key=lambda item: (item['filter'], item['target_id'])), candidate_info


def process_frame(item: dict, cosmo, coverage_floor: float,
                  min_controls: int, control_coverage_floor: float):
    frame = dict(item['frame'])
    record = item['record']
    ident = item['target_id']
    image_path = Path(record['image'])
    if not image_path.is_file() or image_path.stat().st_size != int(record['image_bytes']):
        raise ValueError(f'image path/size disagrees with pinned mask record for {ident}')
    image_hash = sha256(image_path)
    if image_hash != record['image_sha256']:
        raise ValueError(f'image checksum differs from pinned target audit for {ident}')

    info = validate_fits(image_path, frame)
    if info['coadd_header_filter'] != item['filter']:
        raise ValueError(f'primary filter changed for {ident}')
    with core.fits.open(image_path, mode='readonly', memmap=False, checksum=True) as hdus:
        image = np.asarray(hdus[info['science_hdu']].data, dtype=np.float32)
        variance = np.asarray(hdus[info['variance_hdu']].data, dtype=np.float32)
        bits = np.asarray(hdus[info['mask_hdu']].data, dtype=np.int32)
    if image.shape != variance.shape or image.shape != bits.shape:
        raise ValueError(f'image/variance/mask shapes differ for {ident}')

    bad_bitmask = sum(1 << bit for bit in BAD_BIT_POSITIONS)
    bad = (bits & bad_bitmask) != 0
    base = ((bits & (bad_bitmask | BRIGHT_OBJECT_BIT)) != 0)
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    base |= invalid

    group = item['members']
    member_rows = build_member_rows(group)
    member_xy_all = np.asarray([
        [float(row['X_IMAGE']) - 1.0, float(row['Y_IMAGE']) - 1.0]
        for row in member_rows if row['IN_FRAME'] == 'True'
    ], dtype=float)
    center_xy = np.asarray(info['center_pixel_xy'], dtype=float)
    members_for_detector = member_xy_all if len(member_xy_all) else np.empty((0, 2))
    objects, segmentation, kron, kron_flags, bcg, snr, detection_audit = prep.detect(
        image, variance, bad, center_xy, members_for_detector, float(info['pixel_scale_arcsec'])
    )
    detections = detection_rows(objects, kron, kron_flags)
    if not detections:
        raise ValueError(f'no SEP detections for {ident}')
    arrays = core.detection_arrays(detections, float(info['pixel_scale_arcsec']))
    matched, match_stats = core.match_members(member_rows, detections, float(info['pixel_scale_arcsec']))
    matched_by_member = {int(row['member_index']): row for row in matched}

    in_frame_indices = [int(row['INDEX']) for row in member_rows if row['IN_FRAME'] == 'True']
    matched_indices = set(matched_by_member)
    fallback_xy = np.asarray([
        [float(member_rows[index]['X_IMAGE']) - 1.0,
         float(member_rows[index]['Y_IMAGE']) - 1.0]
        for index in in_frame_indices if index not in matched_indices
    ], dtype=float)
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(member_xy_all) if len(member_xy_all) else None
    kpc_per_pixel = float(cosmo.kpc_proper_per_arcmin(float(frame['z'])).value
                          * float(info['pixel_scale_arcsec']) / 60.0)
    conversion = float(core.nJy_factor(info))

    selected_members = item['eligible_members']
    profile_rows, diagnostics = [], []
    for member_index, radius_kpc in selected_members:
        match = matched_by_member.get(member_index)
        if match is None:
            diagnostics.append({
                'target_id': ident, 'filter': item['filter'],
                'member_index': member_index,
                'source_row': member_rows[member_index]['source_row'],
                'cluster_radius_kpc': radius_kpc,
                'sep_match_arcsec': None,
                'profile_bins_retained': 0,
                'n_valid_controls': 0,
                'status': 'no_unique_sep_match_within_1arcsec',
            })
            continue
        xy = np.asarray([match['member_x'], match['member_y']], dtype=float)
        number = int(match['detection_number'])
        target_profile = core.profile_at(
            image, variance, base, segmentation, arrays, fallback_tree, xy,
            float(info['pixel_scale_arcsec']), kpc_per_pixel, conversion,
            number, coverage_floor,
        )
        cluster_radius_pix = radius_kpc / kpc_per_pixel
        control_xy = core.control_centers(
            xy, center_xy, cluster_radius_pix, arrays, member_tree,
            float(info['pixel_scale_arcsec']), kpc_per_pixel, image.shape,
        )
        control_profiles = []
        for point in control_xy:
            result = core.profile_at(
                image, variance, base, segmentation, arrays, fallback_tree, point,
                float(info['pixel_scale_arcsec']), kpc_per_pixel, conversion,
                None, control_coverage_floor,
            )
            if result is not None:
                control_profiles.append(result)
        kept = 0
        for bin_index, (lo, hi) in enumerate(zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])):
            if target_profile is None or target_profile['coverage'][bin_index] < coverage_floor:
                continue
            controls = [row for row in control_profiles
                        if row['coverage'][bin_index] >= control_coverage_floor]
            if len(controls) < min_controls:
                continue
            control_sb = float(np.mean([row['sb'][bin_index] for row in controls]))
            target_sb = float(target_profile['sb'][bin_index])
            profile_rows.append({
                'target_id': ident, 'filter': item['filter'],
                'member_index': member_index,
                'source_row': member_rows[member_index]['source_row'],
                'detection_number': number,
                'member_match_arcsec': float(match['match_arcsec']),
                'cluster_radius_kpc': radius_kpc,
                'radial_bin_kpc': f'{lo:g}-{hi:g}',
                'target_sb_njy_arcsec2': target_sb,
                'control_sb_njy_arcsec2': control_sb,
                'delta_sb_njy_arcsec2': target_sb - control_sb,
                'target_coverage': float(target_profile['coverage'][bin_index]),
                'control_coverage_mean': float(np.mean([r['coverage'][bin_index] for r in controls])),
                'n_valid_controls': len(controls),
                'target_err_nominal_njy_arcsec2': float(target_profile['err'][bin_index]),
                'science_measurement': False,
            })
            kept += 1
        diagnostics.append({
            'target_id': ident, 'filter': item['filter'],
            'member_index': member_index,
            'source_row': member_rows[member_index]['source_row'],
            'cluster_radius_kpc': radius_kpc,
            'sep_match_arcsec': float(match['match_arcsec']),
            'profile_bins_retained': kept,
            'n_valid_controls': len(control_profiles),
            'status': 'profiled' if kept else 'no_bin_passed_coverage',
        })

    frame_record = {
        'target_id': ident,
        'filter': item['filter'],
        'image_sha256': image_hash,
        'n_detections': len(detections),
        'n_member_positions_in_frame': len(member_xy_all),
        'n_member_positions_50_500_kpc': len(selected_members),
        'n_unique_sep_matches_all_in_frame': len(matched),
        'match_audit': match_stats,
        'bcg_match_count': detection_audit['bcg_matches'],
        'bcg_detection_number': bcg + 1 if bcg >= 0 else None,
        'n_member_profile_bins': len(profile_rows),
        'science_measurement': False,
    }
    return profile_rows, diagnostics, frame_record


def summarize(profile_rows, selected, output_plot: Path):
    cluster_rows = []
    by_target_bin: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    icl_lookup = {(item['target_id'], item['filter']): item['icl_contrast_njy_arcsec2']
                  for item in selected}
    for row in profile_rows:
        by_target_bin[(row['target_id'], row['filter'], row['radial_bin_kpc'])].append(
            float(row['delta_sb_njy_arcsec2'])
        )
    cluster_bins = []
    for (target_id, band, radial_bin), values in sorted(by_target_bin.items()):
        cluster_bins.append({
            'target_id': target_id, 'filter': band, 'radial_bin_kpc': radial_bin,
            'n_member_positions': len(values),
            'mean_delta_sb_njy_arcsec2': float(np.mean(values)),
            'median_delta_sb_njy_arcsec2': float(np.median(values)),
            'icl_contrast_50_300_minus_300_800_njy_arcsec2': icl_lookup.get((target_id, band)),
            'science_measurement': False,
        })
    outer = [row for row in cluster_bins if row['radial_bin_kpc'] == '50-75']
    by_filter = {}
    for band in ('i', 'i2'):
        selected_band = [item for item in selected if item['filter'] == band]
        outer_band = [row for row in outer if row['filter'] == band]
        profiles = [row for row in cluster_bins if row['filter'] == band]
        bin_summaries = {}
        for radial_bin in (f'{lo:g}-{hi:g}' for lo, hi in zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])):
            vals = [row['mean_delta_sb_njy_arcsec2'] for row in profiles
                    if row['radial_bin_kpc'] == radial_bin]
            bin_summaries[radial_bin] = {
                'n_target_clusters': len(vals),
                'median_cluster_delta_njy_arcsec2': float(np.median(vals)) if vals else None,
                'q16_q84_cluster_delta_njy_arcsec2': [float(x) for x in np.quantile(vals, [0.16, 0.84])] if vals else None,
            }
        if len(outer_band) >= 3:
            x = np.asarray([row['icl_contrast_50_300_minus_300_800_njy_arcsec2'] for row in outer_band], dtype=float)
            y = np.asarray([row['mean_delta_sb_njy_arcsec2'] for row in outer_band], dtype=float)
            r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) and np.std(y) else None
        else:
            r = None
        by_filter[band] = {
            'selected_clusters': len(selected_band),
            'clusters_with_50_75_kpc_profile': len(outer_band),
            'median_camira_candidates_50_500_kpc_in_frame': float(np.median([
                item['n_member_positions_50_500_kpc_in_frame'] for item in selected_band
            ])) if selected_band else None,
            'descriptive_pearson_r_outer_profile_vs_icl_contrast': r,
            'radial_profile_summary': bin_summaries,
            'interpretation': 'pilot descriptive statistic; no spatial-block uncertainty or detection claim',
        }

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2))
    colors = {'i': '#3568a8', 'i2': '#b65b30'}
    for band, label in (('i', 'HSC-i'), ('i2', 'HSC-i2')):
        rows = [row for row in cluster_bins if row['filter'] == band]
        xs, ys, los, his = [], [], [], []
        for radial_bin in (f'{lo:g}-{hi:g}' for lo, hi in zip(core.R_EDGES_KPC[:-1], core.R_EDGES_KPC[1:])):
            values = [row['mean_delta_sb_njy_arcsec2'] for row in rows if row['radial_bin_kpc'] == radial_bin]
            if not values:
                continue
            center = float(radial_bin.split('-')[0])
            edge = float(radial_bin.split('-')[1])
            xs.append((center + edge) / 2)
            ys.append(float(np.median(values)))
            q16, q84 = np.quantile(values, [0.16, 0.84])
            los.append(float(ys[-1] - q16)); his.append(float(q84 - ys[-1]))
        if xs:
            axes[0].errorbar(xs, ys, yerr=np.asarray([los, his]), marker='o', capsize=2,
                             color=colors[band], label=f'{label} (cluster q16–q84)')
        points = [row for row in outer if row['filter'] == band]
        if points:
            axes[1].scatter(
                [row['icl_contrast_50_300_minus_300_800_njy_arcsec2'] for row in points],
                [row['mean_delta_sb_njy_arcsec2'] for row in points],
                s=28, alpha=.72, color=colors[band], label=f'{label}, N={len(points)}'
            )
    axes[0].axhline(0, color='0.35', lw=.8)
    axes[0].set_xlabel('Member-centered radius [proper kpc]')
    axes[0].set_ylabel('Candidate minus rotated-control SB [nJy arcsec⁻²]')
    axes[0].set_title('Source-aware outer-light profiles')
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].axhline(0, color='0.35', lw=.8)
    axes[1].axvline(0, color='0.35', lw=.8)
    axes[1].set_xlabel('Relative ICL contrast, 50–300 minus 300–800 kpc [nJy arcsec⁻²]')
    axes[1].set_ylabel('Member profile excess, 50–75 kpc [nJy arcsec⁻²]')
    axes[1].set_title('Exploratory cluster-level relation')
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle('S23B CAMIRA member-profile method pilot; exploratory, no science promotion')
    fig.tight_layout()
    fig.savefig(output_plot, dpi=160)
    plt.close(fig)
    return cluster_bins, by_filter


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n-per-filter', type=int, default=24)
    parser.add_argument('--seed', type=int, default=20261107)
    parser.add_argument('--coverage-floor', type=float, default=0.60)
    parser.add_argument('--control-coverage-floor', type=float, default=0.60)
    parser.add_argument('--min-controls', type=int, default=3)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.n_per_filter < 1 or args.min_controls < 1:
        parser.error('n-per-filter and min-controls must be positive')
    output = args.output.resolve()
    if output.exists():
        parser.error(f'refusing to overwrite existing output: {output}')
    for path in (TARGET_MANIFEST, PAIR_CONTRASTS, MEMBER_TABLE, MASK_RECORDS):
        if not path.is_file():
            parser.error(f'missing input: {path}')
    hashes = {
        str(TARGET_MANIFEST): sha256(TARGET_MANIFEST),
        str(PAIR_CONTRASTS): sha256(PAIR_CONTRASTS),
        str(MEMBER_TABLE): sha256(MEMBER_TABLE),
        str(MASK_RECORDS): sha256(MASK_RECORDS),
    }
    expected = {
        str(TARGET_MANIFEST): EXPECTED_TARGET_MANIFEST_SHA256,
        str(PAIR_CONTRASTS): EXPECTED_PAIR_CONTRASTS_SHA256,
        str(MEMBER_TABLE): EXPECTED_MEMBER_TABLE_SHA256,
        str(MASK_RECORDS): EXPECTED_MASK_RECORDS_SHA256,
    }
    if hashes != expected:
        raise ValueError('one or more pinned input SHA-256 values changed')

    manifest = json.loads(TARGET_MANIFEST.read_text(encoding='utf-8'))
    frame_by_id = {str(row['id']): row for row in manifest['clusters']}
    record_by_id = {row['id']: row for row in read_tsv(MASK_RECORDS)}
    pair_rows = read_tsv(PAIR_CONTRASTS)
    member_rows = read_tsv(MEMBER_TABLE)
    cosmo = FlatLambdaCDM(H0=float(manifest['h0']), Om0=float(manifest['omega_m']))
    selected, _ = select_sample(pair_rows, member_rows, frame_by_id, record_by_id,
                                cosmo, args.n_per_filter, args.seed)
    if len(selected) != 2 * args.n_per_filter:
        raise ValueError('selected pilot sample does not contain the requested two-band counts')

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name + '.partial.', dir=output.parent))
    selection_rows = []
    all_profiles, all_diagnostics, frame_records = [], [], []
    for i, item in enumerate(selected, 1):
        print(f'[{i}/{len(selected)}] {item["filter"]} {item["target_id"]}: '
              f'{item["n_member_positions_50_500_kpc_in_frame"]} candidate positions', flush=True)
        profiles, diagnostics, frame_record = process_frame(
            item, cosmo, args.coverage_floor, args.min_controls,
            args.control_coverage_floor,
        )
        all_profiles.extend(profiles)
        all_diagnostics.extend(diagnostics)
        frame_records.append(frame_record)
        selection_rows.append({
            'target_id': item['target_id'], 'filter': item['filter'],
            'z': item['z'], 'richness': item['richness'],
            'icl_contrast_50_300_minus_300_800_njy_arcsec2': item['icl_contrast_njy_arcsec2'],
            'n_member_positions_50_500_kpc_in_frame': item['n_member_positions_50_500_kpc_in_frame'],
            'image_sha256': frame_record['image_sha256'],
        })
    write_tsv(tmp / 'selected_target_frames.tsv', selection_rows)
    write_tsv(tmp / 'member_profile_bins.tsv', all_profiles, [
        'target_id', 'filter', 'member_index', 'source_row', 'detection_number',
        'member_match_arcsec', 'cluster_radius_kpc', 'radial_bin_kpc',
        'target_sb_njy_arcsec2', 'control_sb_njy_arcsec2', 'delta_sb_njy_arcsec2',
        'target_coverage', 'control_coverage_mean', 'n_valid_controls',
        'target_err_nominal_njy_arcsec2', 'science_measurement',
    ])
    write_tsv(tmp / 'member_profile_diagnostics.tsv', all_diagnostics, [
        'target_id', 'filter', 'member_index', 'source_row', 'cluster_radius_kpc',
        'sep_match_arcsec', 'profile_bins_retained', 'n_valid_controls', 'status',
    ])
    cluster_bins, by_filter = summarize(all_profiles, selected, tmp / 'pilot_profiles_and_relation.png')
    write_tsv(tmp / 'cluster_profile_relation.tsv', cluster_bins, [
        'target_id', 'filter', 'radial_bin_kpc', 'n_member_positions',
        'mean_delta_sb_njy_arcsec2', 'median_delta_sb_njy_arcsec2',
        'icl_contrast_50_300_minus_300_800_njy_arcsec2', 'science_measurement',
    ])
    summary = {
        'schema': 'astrafek-s23b-camira-member-outer-profile-pilot-v1',
        'status': 'complete_exploratory_method_pilot',
        'science_measurement': False,
        'images_ready_for_science': False,
        'data_downloaded': False,
        'network_used': False,
        'selection': {
            'targets_per_filter': args.n_per_filter,
            'seed': args.seed,
            'parent': 'targets with saved 50-300 minus 300-800 kpc ICL contrast and at least four in-frame CAMIRA positions at 50-500 proper kpc',
            'member_cut': [MIN_CLUSTER_RADIUS_KPC, MAX_CLUSTER_RADIUS_KPC],
            'detected_member_match_radius_arcsec': core.MATCH_ARCSEC,
        },
        'measurement': {
            'radial_edges_proper_kpc': core.R_EDGES_KPC.tolist(),
            'source_mask': 'retain the matched SEP detection label for each member; mask other SEP segments and Kron ellipses, bad/bright-object pixels, and unmatched member fallback circles',
            'controls': f'up to {core.N_CONTROLS_MAX} azimuthal rotations at the same cluster-centric radius, rejecting detections within {core.CONTROL_SOURCE_CLEARANCE_ARCSEC} arcsec and CAMIRA positions within {core.CONTROL_MEMBER_CLEARANCE_KPC} kpc',
            'coverage_floor': args.coverage_floor,
            'control_coverage_floor': args.control_coverage_floor,
            'minimum_valid_controls_per_bin': args.min_controls,
            'uncertainty': 'nominal independent-pixel variance only; no spatial-block inference in this method pilot',
        },
        'sample': {
            'n_frames': len(frame_records),
            'n_frames_by_filter': {band: sum(row['filter'] == band for row in frame_records) for band in ('i', 'i2')},
            'member_positions_in_selected_cluster_radius': int(sum(row['n_member_positions_50_500_kpc_in_frame'] for row in selection_rows)),
            'unique_sep_matches_all_in_frame': int(sum(row['n_unique_sep_matches_all_in_frame'] for row in frame_records)),
            'member_candidate_diagnostics': len(all_diagnostics),
            'member_profile_bin_rows_retained': len(all_profiles),
            'filter_results': by_filter,
        },
        'limitations': [
            'CAMIRA members are photometric candidates, not spectroscopically confirmed satellites.',
            'Only SEP-matched members are profiled; faint and ambiguous candidates are excluded from the measured profile sample.',
            'Rotated controls are within the same target image and can contain cluster light; they estimate a local baseline, not a blank field.',
            'The pipeline preserves a detected candidate source but does not model its PSF wings or sky-transfer bias.',
            'Coadd pixel covariance, source Poisson noise, Galactic-extinction correction, and spatially blocked uncertainty are not calibrated here.',
            'The random two-band subset is a method pilot and is not a representative full-sample estimate.',
        ],
        'input_sha256': hashes,
        'analysis_script_sha256': sha256(Path(__file__).resolve()),
        'elapsed_seconds': time.time() - started,
        'frame_provenance': frame_records,
    }
    (tmp / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for path in tmp.iterdir():
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    os.replace(tmp, output)
    print(json.dumps({
        'status': summary['status'], 'output': str(output),
        'sample': summary['sample'],
        'summary_sha256': sha256(output / 'summary.json'),
    }, indent=2))


if __name__ == '__main__':
    started = time.time()
    main()
