#!/usr/bin/env python3
"""Audit S23B CAMIRA candidate positions against existing target mask products.

This is a mask-coordinate/coverage audit only. It reads no science-image pixels,
performs no network access, and cannot separate a candidate's own source mask
from neighboring-source masks because the saved products have no segmentation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
MEMBER_TABLE = DATA / 's23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/target_member_positions.tsv'
MASK_RECORDS = DATA / 's23b_icl_mask_support_2026-11-07_01/target/records.tsv'
MASK_SUMMARY = DATA / 's23b_icl_mask_support_2026-11-07_01/target/summary.json'
OUTPUT_DEFAULT = DATA / 's23b_icl_satellite_member_mask_audit_2026-11-07_02'
EXPECTED_MEMBER_SHA256 = 'd1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae'
EXPECTED_MASK_RECORDS_SHA256 = '609354bfa7777088e995a487d0c2eaaed717b80e94589f1c5e0c3758f3e32ff8'
EXPECTED_MASK_SUMMARY_SHA256 = '3ce9f384b98a3c1a56cba291309ef3abe80c68dfc8e74c36b5b74e166a38bc80'
EXPECTED_MEMBER_ROWS = 45634
EXPECTED_MASK_RECORDS = 1039
MASK_EXTENSIONS = ('BASE_MASK', 'MASK_K2P5', 'MASK_K3P5', 'MASK_K5P0')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def to_bool(value: str) -> bool:
    return value.strip().lower() == 'true'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f'refusing to overwrite existing output directory: {output}')
    for path in (MEMBER_TABLE, MASK_RECORDS, MASK_SUMMARY):
        if not path.is_file():
            parser.error(f'missing pinned input: {path}')
    member_sha = sha256(MEMBER_TABLE)
    record_sha = sha256(MASK_RECORDS)
    summary_sha = sha256(MASK_SUMMARY)
    if member_sha != EXPECTED_MEMBER_SHA256:
        parser.error('member-position table checksum differs from the pinned import')
    if record_sha != EXPECTED_MASK_RECORDS_SHA256:
        parser.error('mask record checksum differs from the audited S23B mask set')
    if summary_sha != EXPECTED_MASK_SUMMARY_SHA256:
        parser.error('mask summary checksum differs from the audited S23B mask set')

    members = read_tsv(MEMBER_TABLE)
    records = read_tsv(MASK_RECORDS)
    if len(members) != EXPECTED_MEMBER_ROWS or len(records) != EXPECTED_MASK_RECORDS:
        parser.error('member or mask-record row count differs from the pinned audit')
    by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in members:
        by_target[row['target_cluster_id']].append(row)
    record_by_id = {row['id']: row for row in records}
    if len(record_by_id) != len(records):
        parser.error('duplicate target ID in mask records')

    mask_sha256_by_path: dict[str, str] = {}
    per_member: list[dict] = []
    target_rows: list[dict] = []
    counts = Counter()
    by_filter: dict[str, Counter] = defaultdict(Counter)
    area_fraction_by_filter: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {ext: [] for ext in MASK_EXTENSIONS}
    )
    max_wcs_delta = 0.0
    n_mask_center_samples = 0

    for target_id in sorted(by_target, key=lambda value: int(value.removeprefix('S23B_'))):
        group = by_target[target_id]
        record = record_by_id.get(target_id)
        if record is None:
            counts['targets_without_mask_record'] += 1
            for row in group:
                counts['members_without_mask_record'] += 1
                counts['members_inside_without_mask_record'] += int(
                    to_bool(row['inside_image_pixel_centers'])
                )
                per_member.append({
                    'target_cluster_id': target_id,
                    'source_row': row['source_row'],
                    'primary_filter': '',
                    'inside_image_pixel_centers': row['inside_image_pixel_centers'],
                    'mask_record_status': 'no_successful_mask_record',
                    'mask_center_pixel_x': '', 'mask_center_pixel_y': '',
                    'base_masked_at_nearest_pixel': '',
                    'mask_k2p5_at_nearest_pixel': '',
                    'mask_k3p5_at_nearest_pixel': '',
                    'mask_k5p0_at_nearest_pixel': '',
                    'wcs_delta_pixel_max': '',
                })
            continue

        image_hashes = {row['image_sha256_from_pinned_audit'] for row in group}
        if len(image_hashes) != 1 or image_hashes.pop() != record['image_sha256']:
            raise ValueError(f'image hash join failed for {target_id}')
        mask_path = Path(record['mask_file'])
        if not mask_path.is_file():
            raise FileNotFoundError(mask_path)
        actual_mask_sha = sha256(mask_path)
        if actual_mask_sha != record['mask_sha256']:
            raise ValueError(f'mask checksum mismatch for {target_id}')
        mask_sha256_by_path[str(mask_path)] = actual_mask_sha

        with fits.open(mask_path, mode='readonly', memmap=False) as hdus:
            planes = {}
            for ext in MASK_EXTENSIONS:
                if ext not in hdus or hdus[ext].data is None:
                    raise ValueError(f'{target_id}: missing mask extension {ext}')
                planes[ext] = np.asarray(hdus[ext].data, dtype=bool)
            shape = planes['BASE_MASK'].shape
            if any(plane.shape != shape for plane in planes.values()):
                raise ValueError(f'{target_id}: mask extension shapes disagree')
            wcs = WCS(hdus['BASE_MASK'].header).celestial
            if not wcs.has_celestial:
                raise ValueError(f'{target_id}: saved mask has no celestial WCS')
            ras = np.asarray([float(row['member_ra_deg']) for row in group])
            decs = np.asarray([float(row['member_dec_deg']) for row in group])
            x_calc, y_calc = wcs.world_to_pixel_values(ras, decs)
            x_saved = np.asarray([float(row['image_x_pixel']) for row in group])
            y_saved = np.asarray([float(row['image_y_pixel']) for row in group])
            deltas = np.maximum(np.abs(x_calc - x_saved), np.abs(y_calc - y_saved))
            finite = np.isfinite(deltas)
            if not finite.all():
                raise ValueError(f'{target_id}: nonfinite WCS transform in mask audit')
            target_max_delta = float(deltas.max(initial=0.0))
            max_wcs_delta = max(max_wcs_delta, target_max_delta)
            if target_max_delta > 1e-3:
                raise ValueError(f'{target_id}: mask WCS differs from imported image WCS')

            ny, nx = shape
            filter_name = record['primary_filter']
            for ext in MASK_EXTENSIONS:
                area_fraction_by_filter[filter_name][ext].append(
                    float(np.mean(planes[ext]))
                )
            inside = (
                np.isfinite(x_saved) & np.isfinite(y_saved)
                & (x_saved >= 0) & (x_saved <= nx - 1)
                & (y_saved >= 0) & (y_saved <= ny - 1)
            )
            if not np.array_equal(inside, np.asarray([
                to_bool(row['inside_image_pixel_centers']) for row in group
            ])):
                raise ValueError(f'{target_id}: saved footprint flags disagree with mask geometry')
            valid_indices = np.flatnonzero(inside)
            ix = np.floor(x_saved[valid_indices] + 0.5).astype(int)
            iy = np.floor(y_saved[valid_indices] + 0.5).astype(int)
            if np.any(ix < 0) or np.any(ix >= nx) or np.any(iy < 0) or np.any(iy >= ny):
                raise ValueError(f'{target_id}: nearest-pixel indices are out of bounds')
            point_masks = {ext: planes[ext][iy, ix] for ext in MASK_EXTENSIONS}
            for local_index, row in enumerate(group):
                result = {
                    'target_cluster_id': target_id,
                    'source_row': row['source_row'],
                    'primary_filter': filter_name,
                    'inside_image_pixel_centers': bool(inside[local_index]),
                    'mask_record_status': 'mask_available',
                    'mask_center_pixel_x': float(x_saved[local_index]),
                    'mask_center_pixel_y': float(y_saved[local_index]),
                    'base_masked_at_nearest_pixel': '',
                    'mask_k2p5_at_nearest_pixel': '',
                    'mask_k3p5_at_nearest_pixel': '',
                    'mask_k5p0_at_nearest_pixel': '',
                    'wcs_delta_pixel_max': float(deltas[local_index]),
                }
                if inside[local_index]:
                    member_index = int(np.searchsorted(valid_indices, local_index))
                    result.update({
                        'base_masked_at_nearest_pixel': bool(point_masks['BASE_MASK'][member_index]),
                        'mask_k2p5_at_nearest_pixel': bool(point_masks['MASK_K2P5'][member_index]),
                        'mask_k3p5_at_nearest_pixel': bool(point_masks['MASK_K3P5'][member_index]),
                        'mask_k5p0_at_nearest_pixel': bool(point_masks['MASK_K5P0'][member_index]),
                    })
                    n_mask_center_samples += 1
                    counts['in_footprint_with_mask'] += 1
                    for ext, key in (
                        ('BASE_MASK', 'base_masked_at_nearest_pixel'),
                        ('MASK_K2P5', 'mask_k2p5_at_nearest_pixel'),
                        ('MASK_K3P5', 'mask_k3p5_at_nearest_pixel'),
                        ('MASK_K5P0', 'mask_k5p0_at_nearest_pixel'),
                    ):
                        by_filter[filter_name][f'{key}_n'] += 1
                        if result[key]:
                            by_filter[filter_name][f'{key}_yes'] += 1
                else:
                    counts['outside_image_pixel_centers'] += 1
                per_member.append(result)

        target_rows.append({
            'target_cluster_id': target_id,
            'primary_filter': record['primary_filter'],
            'n_member_rows': len(group),
            'n_inside_image_pixel_centers': sum(to_bool(r['inside_image_pixel_centers']) for r in group),
            'n_masked_by_base_at_center': sum(
                bool(per_member[-len(group) + idx]['base_masked_at_nearest_pixel'])
                for idx in range(len(group))
                if per_member[-len(group) + idx]['base_masked_at_nearest_pixel'] != ''
            ),
            'image_sha256_matches_mask_record': True,
            'mask_sha256_verified': True,
            'mask_wcs_max_delta_pixel': target_max_delta,
            'base_mask_image_area_fraction': float(np.mean(planes['BASE_MASK'])),
            'mask_k2p5_image_area_fraction': float(np.mean(planes['MASK_K2P5'])),
            'mask_k3p5_image_area_fraction': float(np.mean(planes['MASK_K3P5'])),
            'mask_k5p0_image_area_fraction': float(np.mean(planes['MASK_K5P0'])),
        })
        counts['targets_with_mask_record'] += 1

    if len(per_member) != len(members):
        raise ValueError('member-level output row count mismatch')
    counts['outside_image_pixel_centers_total'] = sum(
        not to_bool(row['inside_image_pixel_centers']) for row in members
    )
    counts['outside_image_pixel_centers_with_mask'] = counts['outside_image_pixel_centers']
    counts['members_inside_with_mask'] = n_mask_center_samples

    summary_by_filter = {}
    for filter_name, item in sorted(by_filter.items()):
        summary_by_filter[filter_name] = {}
        for ext in MASK_EXTENSIONS:
            stem = {
                'BASE_MASK': 'base_masked_at_nearest_pixel',
                'MASK_K2P5': 'mask_k2p5_at_nearest_pixel',
                'MASK_K3P5': 'mask_k3p5_at_nearest_pixel',
                'MASK_K5P0': 'mask_k5p0_at_nearest_pixel',
            }[ext]
            n = int(item[f'{stem}_n'])
            yes = int(item[f'{stem}_yes'])
            summary_by_filter[filter_name][stem] = {
                'masked': yes, 'sampled': n, 'fraction': yes / n if n else None,
            }
        summary_by_filter[filter_name]['median_target_image_area_fraction'] = {
            ext: float(np.median(area_fraction_by_filter[filter_name][ext]))
            for ext in MASK_EXTENSIONS
        }
        summary_by_filter[filter_name]['target_image_area_fraction_q10_q50_q90'] = {
            ext: [float(value) for value in np.quantile(
                area_fraction_by_filter[filter_name][ext], [0.1, 0.5, 0.9]
            )]
            for ext in MASK_EXTENSIONS
        }

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name + '.tmp.', dir=output.parent))
    member_fields = [
        'target_cluster_id', 'source_row', 'primary_filter', 'inside_image_pixel_centers',
        'mask_record_status', 'mask_center_pixel_x', 'mask_center_pixel_y',
        'base_masked_at_nearest_pixel', 'mask_k2p5_at_nearest_pixel',
        'mask_k3p5_at_nearest_pixel', 'mask_k5p0_at_nearest_pixel', 'wcs_delta_pixel_max',
    ]
    target_fields = [
        'target_cluster_id', 'primary_filter', 'n_member_rows',
        'n_inside_image_pixel_centers', 'n_masked_by_base_at_center',
        'image_sha256_matches_mask_record', 'mask_sha256_verified', 'mask_wcs_max_delta_pixel',
        'base_mask_image_area_fraction', 'mask_k2p5_image_area_fraction',
        'mask_k3p5_image_area_fraction', 'mask_k5p0_image_area_fraction',
    ]
    write_tsv(tmp / 'member_mask_overlap.tsv', per_member, member_fields)
    write_tsv(tmp / 'target_mask_overlap_summary.tsv', target_rows, target_fields)
    document = {
        'schema': 'astrafek-s23b-camira-member-mask-overlap-v1',
        'status': 'complete_exploratory_mask_coordinate_audit',
        'science_measurement': False,
        'images_ready_for_science': False,
        'mask_semantics': {
            'BASE_MASK': 'bad-quality, nonfinite/nonpositive-variance, and bright-object pixels from the existing mask preparation',
            'MASK_K2P5_MASK_K3P5_MASK_K5P0': 'BASE_MASK plus SEP segmentation/Kron source exclusions at the named Kron factor',
            'limitation': 'saved mask products contain no segmentation labels, so a candidate own-source mask cannot be distinguished from a neighboring-source mask',
        },
        'scope': {
            'member_positions': len(members),
            'successful_target_masks': len(records),
            'members_with_center_mask_sample': n_mask_center_samples,
            'targets_with_mask_record': counts['targets_with_mask_record'],
            'targets_without_mask_record': counts['targets_without_mask_record'],
            'members_inside_without_mask_record': counts['members_inside_without_mask_record'],
            'members_without_mask_record': counts['members_without_mask_record'],
            'outside_image_pixel_centers_total': counts['outside_image_pixel_centers_total'],
            'outside_image_pixel_centers_with_mask': counts['outside_image_pixel_centers_with_mask'],
            'max_mask_wcs_delta_pixel': max_wcs_delta,
            'science_image_pixels_read': False,
            'mask_pixels_read': True,
            'network_used': False,
        },
        'masked_fraction_by_filter_and_plane': summary_by_filter,
        'input_sha256': {
            str(MEMBER_TABLE): member_sha,
            str(MASK_RECORDS): record_sha,
            str(MASK_SUMMARY): summary_sha,
        },
        'per_mask_sha256_verified_count': len(mask_sha256_by_path),
        'per_mask_sha256_verified': mask_sha256_by_path,
        'analysis_script_sha256': sha256(Path(__file__).resolve()),
        'interpretation': [
            'A center marked in BASE_MASK has no valid center pixel under this mask definition.',
            'Overlap with a Kron source mask does not by itself mean the candidate is unusable; it may be the candidate galaxy itself.',
            'The current mask product cannot support member-centered profile extraction without source segmentation labels or a regenerated candidate-excluded mask.',
            'These are CAMIRA candidate positions, not spectroscopically confirmed satellites.',
        ],
    }
    (tmp / 'summary.json').write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    os.replace(tmp, output)
    print(json.dumps({
        'status': document['status'],
        'output': str(output),
        'scope': document['scope'],
        'masked_fraction_by_filter_and_plane': summary_by_filter,
        'summary_sha256': sha256(output / 'summary.json'),
        'member_mask_overlap_sha256': sha256(output / 'member_mask_overlap.tsv'),
        'target_summary_sha256': sha256(output / 'target_mask_overlap_summary.tsv'),
    }, indent=2))


if __name__ == '__main__':
    main()
