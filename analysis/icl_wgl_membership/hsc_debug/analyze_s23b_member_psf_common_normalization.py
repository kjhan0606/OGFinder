#!/usr/bin/env python3
"""Separate PSF recovery change from changed nominal-input normalization.

Uses completed, hash-pinned research ledgers only. No images, catalog queries or
new downloads; support is conditional on both experiments passing selection.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import compare_s23b_member_psf_injection as paired

BASE = paired.BASE
NEW = BASE / 's23b_member_psf_injection_run_2026-11-17_01/injection'
NEW_SHA = '7cdd34b10f3e976a1ecd94518c338e286aaa59fd4b42c310c3607da4c27431c4'
COMPARE = BASE / 's23b_member_psf_injection_run_2026-11-17_01/comparison/summary.json'
COMPARE_SHA = '3ecb1ad6900ac6255b4143b461782a3978e020e3dbed9e22b9209e534bec19cf'


def components(old, new):
    d0 = float(old['nominal_unmasked_50_75_mean_njy_arcsec2'])
    d1 = float(new['nominal_unmasked_50_75_mean_njy_arcsec2'])
    f0 = float(old['recovered_delta_increment_njy_arcsec2'])
    f1 = float(new['recovered_delta_increment_njy_arcsec2'])
    if not all(map(math.isfinite, (d0, d1, f0, f1))) or min(d0, d1) <= 0:
        raise ValueError('nonfinite flux or nonpositive nominal input')
    common = (f1 - f0) / d0
    denominator = f1 / d1 - f1 / d0
    observed = f1 / d1 - f0 / d0
    if not math.isclose(common + denominator, observed, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError('normalization decomposition failed')
    for source, flux, denom in ((old, f0, d0), (new, f1, d1)):
        if not math.isclose(float(source['recovered_increment_over_nominal_input']), flux / denom,
                            rel_tol=1e-8, abs_tol=1e-10):
            raise ValueError('stored response does not match flux / nominal input')
    return dict(unconvolved_nominal=d0, convolved_nominal=d1,
                nominal_fractional_change=d1 / d0 - 1,
                common_denominator_response_change=common,
                denominator_contribution=denominator, original_ratio_change=observed)


def run(output, replicas):
    if output.exists():
        raise FileExistsError(output)
    for path, expected in ((paired.OLD / 'summary.json', paired.OLD_SUMMARY_SHA),
                           (NEW / 'summary.json', NEW_SHA), (COMPARE, COMPARE_SHA),
                           (paired.TARGET, paired.TARGET_SHA)):
        if paired.sha(path) != expected:
            raise ValueError(f'pinned input changed: {path}')
    comparison = json.loads(COMPARE.read_text())
    if comparison['science_measurement'] is not False or comparison['status'] != 'complete_exploratory_psf_core_paired_comparison':
        raise ValueError('unexpected paired comparison status')
    for d in (paired.OLD, NEW):
        summary = json.loads((d / 'summary.json').read_text())
        if summary['science_measurement'] is not False or summary['n_frames_complete'] != 427 or summary['n_frames_failed'] != 0:
            raise ValueError('incomplete injection input')
        for filename in ('injection_responses.tsv', 'selected_members.tsv'):
            if paired.sha(d / filename) != summary['output_sha256'][filename]:
                raise ValueError('injection product hash mismatch')
    old = paired.read(paired.OLD / 'injection_responses.tsv')
    new = paired.read(NEW / 'injection_responses.tsv')
    matched = paired.paired_rows(old, new)
    if len(matched) != 1281:
        raise ValueError('expected exactly 427 members at three amplitudes')
    index_old = {paired.key(r): r for r in old}
    index_new = {paired.key(r): r for r in new}
    rows = []
    for item in matched:
        if not item['both_supported']:
            continue
        key = (item['target_id'], item['filter'], item['member_index'], item['source_row'], item['amplitude'])
        rows.append(dict(target_id=item['target_id'], filter=item['filter'], amplitude=item['amplitude'],
                         member_index=item['member_index'], source_row=item['source_row'],
                         **components(index_old[key], index_new[key])))
    frames = {str(r['id']): r for r in json.loads(paired.TARGET.read_text())['clusters']}
    results = []
    metrics = ('nominal_fractional_change', 'common_denominator_response_change',
               'denominator_contribution', 'original_ratio_change')
    for band_no, band in enumerate(('i', 'i2')):
        for amp_no, amp in enumerate((.5, 1., 2.)):
            subset = [r for r in rows if r['filter'] == band and r['amplitude'] == amp]
            expected = next(x for x in comparison['comparisons'] if x['filter'] == band and x['amplitude'] == amp and x['metric'] == 'delta_response')
            if len(subset) != expected['n_both']:
                raise ValueError('paired support count differs from saved comparison')
            for metric in metrics:
                stats = paired.median_block_interval(subset, metric, frames, replicas, 20261117 + band_no * 10 + amp_no)
                results.append(dict(filter=band, amplitude=amp, n_both=len(subset), metric=metric, **stats))
    output.mkdir(parents=True, exist_ok=False)
    ledger = output / 'common_normalization_paired.tsv'
    with ledger.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    result = dict(schema='astrafek-member-psf-common-normalization-v1', status='complete_exploratory',
                  science_measurement=False, network_used=False, images_read=False,
                  input_sha256={'old_summary': paired.OLD_SUMMARY_SHA, 'new_summary': NEW_SHA,
                                'paired_summary': COMPARE_SHA, 'target_manifest': paired.TARGET_SHA},
                  script_sha256=paired.sha(Path(__file__)), helper_sha256=paired.sha(Path(paired.__file__)),
                  output_sha256={ledger.name: paired.sha(ledger)}, n_paired_supported=len(rows),
                  bootstrap_replicas=replicas, results=results,
                  limitations=['Selected on support in both experiments; not a completeness estimate.',
                    'Individual decomposition terms sum exactly, but their medians need not.',
                    'PSF core only; no source Poisson, extended wings, coadd weights, or calibrated sky/covariance.',
                    'Intervals are marginal exploratory sky-cell bootstraps, not multiplicity-adjusted detection tests.'])
    (output / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replicas', type=int, default=5000)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error('at least 100 bootstrap replicas required')
    run(args.output.resolve(), args.replicas)
