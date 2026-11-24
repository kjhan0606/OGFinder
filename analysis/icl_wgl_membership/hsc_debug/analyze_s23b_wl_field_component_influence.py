#!/usr/bin/env python3
"""Decompose fixed-52 WL leave-field influence into saved lens/control terms.

No source rows, image reads, CAS queries, or downloads. Official polygon field
labels are applied only to 52 saved lens-center coordinates. The exact earlier
field-deletion result additionally removes source/control endpoints and may
change support; its residual versus lens-only omission is explicitly retained.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math

import numpy as np
from hsc_pdr3_tract_patch_geometry import load_index

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
FIELD = DATA / 'wl_s23b_fixed52_exact_field_influence_2026-11-17_01'
FIELD_SHA = '422afe04d11cd10f5bb047f02dae35d0990a1890615ab7b5c49e3abb716c7d33'
REF = DATA / 'wl_s23b_equal_random_center_spatial_covariance_2026-11-11_02'
REF_SHA = '00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033'
TARGET = DATA / 's23b_icl_target_image_plan_2026-11-07_01/target_manifest.json'
TARGET_SHA = 'e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44'


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def decompose(rows, removed, component):
    lens_key, control_key = (('lens_delta_sigma_msun_mpc2', 'equal_center_random_delta_sigma_msun_mpc2')
                             if component == 'tangential' else ('lens_cross_msun_mpc2', 'equal_center_random_cross_msun_mpc2'))
    all_lens = np.asarray([float(r[lens_key]) for r in rows])
    all_control = np.asarray([float(r[control_key]) for r in rows])
    keep = np.asarray([r['field'] != removed for r in rows])
    if not keep.any() or keep.all() or not np.all(np.isfinite(all_lens)) or not np.all(np.isfinite(all_control)):
        raise ValueError('invalid field selection or nonfinite values')
    dl = float(all_lens[keep].mean() - all_lens.mean())
    dc = float(all_control[keep].mean() - all_control.mean())
    full = float((all_lens - all_control).mean())
    leave = float((all_lens[keep] - all_control[keep]).mean())
    if not math.isclose(leave - full, dl - dc, rel_tol=1e-12, abs_tol=1.0):
        raise ValueError('lens and control decomposition mismatch')
    return dict(full_lens=float(all_lens.mean()), full_control=float(all_control.mean()),
                lens_only_omission_profile=leave, full_profile=full,
                lens_component_shift=dl, minus_control_component_shift=-dc,
                lens_only_omission_shift=leave-full, n_retained=int(keep.sum()))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    for path, expected in ((FIELD/'summary.json', FIELD_SHA), (REF/'summary.json', REF_SHA), (TARGET, TARGET_SHA)):
        if sha(path) != expected:
            raise ValueError(f'pinned input differs: {path}')
    field = json.loads((FIELD/'summary.json').read_text())
    ref = json.loads((REF/'summary.json').read_text())
    if field['status'] != 'complete_exact_field_influence_only' or field['science_measurement'] is not False:
        raise ValueError('unexpected official-field result')
    if sha(FIELD/'region_profile_shifts.tsv') != field['outputs']['region_profile_shifts_sha256']:
        raise ValueError('official field shift file differs')
    profiles_path = REF/'per_lens_set_profiles.tsv'
    if sha(profiles_path) != ref['output_sha256'][profiles_path.name]:
        raise ValueError('per-lens file differs')
    geometry_dir = Path(field['inputs']['official_geometry_directory'])
    for name, meta in field['official_geometry_summary']['files'].items():
        if sha(geometry_dir/name) != meta['sha256']:
            raise ValueError('official geometry file changed')
    with profiles_path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    ids = sorted({r['matched_lens_id'] for r in rows})
    if len(ids) != 52 or len(rows) != 52 * 7 or len({(r['matched_lens_id'], r['radial_bin']) for r in rows}) != len(rows):
        raise ValueError('unexpected per-lens radial coverage')
    target = json.loads(TARGET.read_text())
    centers = {str(r['id']): r for r in target['clusters']}
    if not set(ids).issubset(centers):
        raise ValueError('lens ID missing in pinned target catalog')
    index, _ = load_index(geometry_dir, level='tract')
    result = index.assign([float(centers[i]['ra_deg']) for i in ids],
                          [float(centers[i]['dec_deg']) for i in ids], level='tract')
    if any(x is None for x in result.field) or np.any(result.cross_field_ambiguous) or np.any(result.tract < 0):
        raise ValueError('unassigned or cross-field ambiguous lens center')
    lens_fields = dict(zip(ids, result.field))
    counts = {r['region_label']: r['lens_centers_in_region'] for r in field['region_deletions']}
    if {f: list(lens_fields.values()).count(f) for f in counts} != counts:
        raise ValueError('lens-center field counts disagree with prior exact run')
    for r in rows:
        r['field'] = lens_fields[r['matched_lens_id']]
    with (FIELD/'region_profile_shifts.tsv').open(newline='') as f:
        exact = {(r['region_label'], r['radial_bin']): r for r in csv.DictReader(f, delimiter='\t')
                 if r['estimator'] == 'equal_control_centers'}
    output_rows = []
    for band in range(7):
        radial = [r for r in rows if int(r['radial_bin']) == band]
        if len(radial) != 52:
            raise ValueError('missing lens in radius bin')
        for region in sorted(counts):
            saved = exact[(region, str(band))]
            for component in ('tangential', 'cross'):
                stats = decompose(radial, region, component)
                full_key = 'full_delta_sigma_msun_mpc2' if component == 'tangential' else 'full_cross_msun_mpc2'
                shift_key = 'delta_sigma_shift_msun_mpc2' if component == 'tangential' else 'cross_shift_msun_mpc2'
                if not math.isclose(stats['full_profile'], float(saved[full_key]), rel_tol=1e-12, abs_tol=1.0):
                    raise ValueError('saved full profile is not reproduced')
                exact_shift = float(saved[shift_key])
                output_rows.append(dict(field=region, radial_bin=band, r_min_mpc=float(saved['r_min_mpc']),
                    r_max_mpc=float(saved['r_max_mpc']), component=component,
                    n_exact_retained=int(saved['common_lens_sets']), **stats,
                    exact_endpoint_omission_shift=exact_shift,
                    endpoint_and_support_residual=exact_shift-stats['lens_only_omission_shift']))
    output.mkdir(parents=True, exist_ok=False)
    ledger = output/'field_component_decomposition.tsv'
    with ledger.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(output_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(output_rows)
    summary = dict(status='complete_exploratory_field_component_decomposition', science_measurement=False,
        network_used=False, source_rows_read=False, images_read=False, n_lens=52,
        lens_field_counts=counts, output_rows=len(output_rows), input_sha256={'field_summary': FIELD_SHA,
        'reference_summary': REF_SHA, 'target_manifest': TARGET_SHA},
        script_sha256=sha(Path(__file__)), geometry_parser_sha256=sha(Path(__file__).with_name('hsc_pdr3_tract_patch_geometry.py')),
        output_sha256={ledger.name: sha(ledger)}, results=output_rows,
        limitations=['Field grouping by lens-center is descriptive and not a causal decomposition.',
        'Residual between endpoint deletion and lens-only exclusion combines removed source/control endpoints, estimator changes and support changes.',
        'Three fields and sample selection do not support a stable survey covariance or significance claim.',
        'Shape calibration, n(z), source reuse, PSF, masks and contamination remain imperfect.'])
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':summary['status'], 'n_lens':len(ids), 'output':str(output),
                      'summary_sha256':sha(output/'summary.json')}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    run(args.output.resolve())
