#!/usr/bin/env python3
"""Prepare an offline request for missing position PSFs in the 224-pair ICL sample.

Only derives coordinates from already-downloaded images/catalogue products. It
never contacts HSC, submits a request, or downloads data.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from prepare_hsc_psf_picker_plan import frame_positions

DEFAULT_DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
DEFAULT_ICL_SUPPORT = DEFAULT_DATA / 'icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv'
DEFAULT_PROFILE_SUMMARY = DEFAULT_DATA / 'satellite_outer_profile_perbin_cov040_2026-11-05_01/summary.json'
DEFAULT_MASKS = DEFAULT_DATA / 'prestack_cluster_masks_2026-10-16_01/records.tsv'
DEFAULT_PREP = DEFAULT_DATA / 'satellite_detection_psf_2026-10-18_01'
DEFAULT_OLD_INJECTION_PLAN = DEFAULT_DATA / 'icl_psf_picker_injection_plan_2026-10-19_02'
DEFAULT_OLD_INJECTION_DOWNLOAD = DEFAULT_DATA / 'icl_psf_picker_injection_2026-10-19_01'
DEFAULT_CENTER_PLAN = DEFAULT_DATA / 'icl_psf_picker_target_control_centers_plan_2026-10-21_01'
DEFAULT_CENTER_DOWNLOAD = DEFAULT_DATA / 'icl_psf_picker_target_control_centers_download_2026-10-21_01'
DEFAULT_OUTPUT = DEFAULT_DATA / 'icl_psf_picker_satellite_support224_plan_2026-11-05_01'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def key(row):
    return (str(row['cluster_id']), str(row['filter']))


def coord_key(row):
    ra, dec = float(row['ra_deg']), float(row['dec_deg'])
    if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError(f'invalid sky position {ra}, {dec}')
    return round(ra, 9), round(dec, 9)


def fits_audit(mapping_path: Path, download: Path, audit_path: Path, expected_mapping_sha: str):
    mapping = json.loads(mapping_path.read_text())
    audit = json.loads(audit_path.read_text())
    files = sorted(download.glob('*.fits'))
    if (sha256(mapping_path) != expected_mapping_sha
            or audit.get('status') != 'verified'
            or audit.get('n_returned') != len(mapping)
            or audit.get('n_matched_coordinates') != len(mapping)
            or len(files) != len(mapping)
            or sum(p.stat().st_size for p in files) != audit.get('total_bytes')):
        raise ValueError(f'existing PSF bank failed mapping/file/audit consistency: {download}')
    return mapping, audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, default=DEFAULT_DATA)
    p.add_argument('--icl-support', type=Path, default=DEFAULT_ICL_SUPPORT)
    p.add_argument('--profile-summary', type=Path, default=DEFAULT_PROFILE_SUMMARY)
    p.add_argument('--masks', type=Path, default=DEFAULT_MASKS)
    p.add_argument('--satellite-prep', type=Path, default=DEFAULT_PREP)
    p.add_argument('--old-injection-plan', type=Path, default=DEFAULT_OLD_INJECTION_PLAN)
    p.add_argument('--old-injection-download', type=Path, default=DEFAULT_OLD_INJECTION_DOWNLOAD)
    p.add_argument('--center-plan', type=Path, default=DEFAULT_CENTER_PLAN)
    p.add_argument('--center-download', type=Path, default=DEFAULT_CENTER_DOWNLOAD)
    p.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    a = p.parse_args()
    output = a.output.resolve()
    if output.exists():
        raise FileExistsError(f'refusing to overwrite {output}')

    manifest_path = a.data / 'sample_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    clusters = {str(r['id']): r for r in manifest['clusters']}
    if len(clusters) != 800:
        raise ValueError(f'expected 800 unique image records, got {len(clusters)}')
    icl_rows = read_tsv(a.icl_support)
    icl_keys = {(str(r['id']), str(r['filter'])) for r in icl_rows}
    if len(icl_keys) != 224:
        raise ValueError(f'expected 224 supported ICL frame/filter pairs, got {len(icl_keys)}')
    profile = json.loads(a.profile_summary.read_text())
    if (profile.get('status') != 'complete' or profile.get('science_measurement') is not False
            or profile.get('method', {}).get('coverage_floor') != 0.4
            or profile.get('n_frames_completed') != 800):
        raise ValueError('40% profile summary status/scope does not match expected audited input')
    records = {(str(r['cluster_id']), str(r['filter'])): r for r in profile['frame_records']}
    if not icl_keys.issubset(records):
        raise ValueError('ICL-supported pairs are missing satellite profile summary rows')
    mask_rows = read_tsv(a.masks)
    masks = {str(r['id']): r for r in mask_rows}
    if len(masks) != len(mask_rows):
        raise ValueError('duplicate mask frame IDs')
    frames = []
    for ident, band in sorted(icl_keys):
        record = records[(ident, band)]
        frame = clusters[ident]
        if str(frame.get('filter', frame.get('primary_filter', ''))).lower() not in ('', band.lower()):
            raise ValueError(f'manifest filter mismatch for {ident}/{band}')
        if ident not in masks:
            raise ValueError(f'missing mask/image for {ident}')
        frame = dict(frame)
        frame['image'] = masks[ident]['image']
        expected_n = int(record['n_matches_50_500kpc'])
        fake_injection_record = {
            'cluster_id': ident, 'filter': band,
            'support': {'n_injected_targets_in_50_500_kpc': expected_n},
        }
        frames.extend(frame_positions(frame, fake_injection_record, masks[ident],
                                      a.satellite_prep, manifest.get('h0', 70.0),
                                      manifest.get('omega_m', 0.3)))

    positions = {}
    for row in frames:
        ck = coord_key(row)
        if ck not in positions:
            positions[ck] = {**row, 'roles': []}
        positions[ck]['roles'].append({
            'cluster_id': row['cluster_id'], 'frame_filter': row['frame_filter'],
            'role': row['role'], 'target_index': row.get('target_index'),
            'projected_radius_kpc': row.get('projected_radius_kpc'),
        })
    if sum(1 for x in frames if x['role'] == 'satellite_injection_center') != sum(
            int(records[k]['n_matches_50_500kpc']) for k in icl_keys):
        raise ValueError('member coordinate count failed frame-summary consistency')

    old_map_path = a.old_injection_plan / 'position_mapping.json'
    old_plan = json.loads((a.old_injection_plan / 'psf_picker_plan.json').read_text())
    old_audit_path = a.old_injection_download / 'psf_picker_audit.json'
    old_mapping, old_audit = fits_audit(old_map_path, a.old_injection_download,
                                         old_audit_path, old_plan['position_mapping_sha256'])
    center_map_path = a.center_plan / 'position_mapping.json'
    center_plan = json.loads((a.center_plan / 'psf_picker_plan.json').read_text())
    center_audit_path = a.center_download / 'psf_picker_audit_reproducible.json'
    center_mapping, center_audit = fits_audit(
        center_map_path, a.center_download, center_audit_path,
        center_plan['position_mapping_sha256'])
    existing = {coord_key(x) for x in old_mapping + center_mapping}
    covered = set(positions) & existing
    for ck in covered:
        positions[ck]['roles'] = positions[ck]['roles']
    missing = [positions[k] for k in sorted(set(positions) - existing)]
    missing_satellite_positions = sum(
        any(role['role'] == 'satellite_injection_center' for role in row['roles'])
        for row in missing)
    missing_target_centers = sum(
        any(role['role'] == 'cluster_center' for role in row['roles'])
        for row in missing)
    if missing_target_centers:
        raise ValueError(f'{missing_target_centers} of 224 required cluster-center PSFs are missing')
    total_sat_positions = sum(
        any(role['role'] == 'satellite_injection_center' for role in row['roles'])
        for row in positions.values())
    old_satellite_positions = sum(
        any(role.get('role') == 'satellite_injection_center' for role in row.get('roles', []))
        for row in old_mapping)
    if missing_satellite_positions + sum(
        any(role['role'] == 'satellite_injection_center' for role in positions[k]['roles'])
        for k in covered) != total_sat_positions:
        raise ValueError('satellite PSF position coverage arithmetic failed')

    sizes = [p.stat().st_size for p in list(a.old_injection_download.glob('*.fits'))
             + list(a.center_download.glob('*.fits'))]
    median_size = statistics.median(sizes)
    output.mkdir(parents=True, exist_ok=False)
    coordinates = output / 'psf_picker_coordinates.txt'
    with coordinates.open('x', encoding='ascii', newline='') as f:
        f.write('#? ra dec\n')
        for row in missing:
            f.write(f"{row['ra_deg']:.9f} {row['dec_deg']:+.9f}\n")
    mapping_path = output / 'position_mapping.json'
    mapping_path.write_text(json.dumps(missing, indent=2, ensure_ascii=False) + '\n')
    plan = {
        'schema': 'astrafek-icl-supported-pairs-psf-request-plan-v1',
        'status': 'prepared_not_submitted',
        'purpose': 'position-dependent HSC-I coadd PSF cores for satellite injection/recovery in 224 ICL-supported cluster/filter pairs',
        'n_supported_icl_pairs': len(icl_keys),
        'n_unique_required_positions': len(positions),
        'n_unique_satellite_injection_positions': total_sat_positions,
        'n_unique_positions_already_covered': len(covered),
        'n_satellite_positions_already_covered': total_sat_positions - missing_satellite_positions,
        'n_satellite_positions_missing': missing_satellite_positions,
        'n_cluster_center_positions_missing': missing_target_centers,
        'n_requests_remaining': len(missing),
        'n_previous_injection_psfs': old_audit['n_returned'],
        'n_previous_target_control_center_psfs': center_audit['n_returned'],
        'old_injection_psf_audit_sha256': sha256(old_audit_path),
        'target_control_psf_audit_sha256': sha256(center_audit_path),
        'existing_mapping_sha256': {
            str(old_map_path): sha256(old_map_path), str(center_map_path): sha256(center_map_path)},
        'estimated_bytes_median_per_psf': median_size,
        'estimated_total_bytes_median': round(len(missing) * median_size),
        'estimated_total_mib_median': round(len(missing) * median_size / 1024**2, 2),
        'request': {'rerun': 'pdr3_wide', 'type': 'coadd', 'filter': 'HSC-I', 'centered': True},
        'coordinate_list': str(coordinates), 'coordinate_list_sha256': sha256(coordinates),
        'position_mapping': str(mapping_path), 'position_mapping_sha256': sha256(mapping_path),
        'inputs': {
            'sample_manifest': str(manifest_path), 'sample_manifest_sha256': sha256(manifest_path),
            'icl_support': str(a.icl_support), 'icl_support_sha256': sha256(a.icl_support),
            'satellite_profile_summary': str(a.profile_summary), 'satellite_profile_summary_sha256': sha256(a.profile_summary),
            'satellite_preparation': str(a.satellite_prep),
        },
        'images_requested': False, 'shape_rows_requested': False,
        'photoz_requested': False, 'credentials_saved': False,
        'network_used': False, 'request_submitted': False,
        'requires_explicit_user_approval': True,
        'science_measurement': False,
        'limitations': [
            'PSF cores support the injection/recovery validation only; they do not expand the cluster image sample.',
            'The 224-pair sample is selected by the existing ICL 200-300 kpc support table.',
            'A completed request must pass exact coordinate, file-count, FITS, and position-bank audits before analysis.',
        ],
    }
    plan_path = output / 'psf_request_plan.json'
    plan_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: plan[k] for k in (
        'status', 'n_supported_icl_pairs', 'n_unique_required_positions',
        'n_unique_satellite_injection_positions', 'n_satellite_positions_already_covered',
        'n_satellite_positions_missing', 'n_cluster_center_positions_missing',
        'n_requests_remaining', 'estimated_total_mib_median',
        'coordinate_list_sha256', 'position_mapping_sha256')}, indent=2))
    print(f"plan_sha256={sha256(plan_path)}")

if __name__ == '__main__':
    main()
