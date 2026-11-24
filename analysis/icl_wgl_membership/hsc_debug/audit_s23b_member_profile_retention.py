#!/usr/bin/env python3
"""Audit saved S23B member-profile attrition without image reads or network IO.

Retention is conditional on the upstream selected frames and photometric member
candidates. Fractions are not detection completeness. Balance contrasts are
purely descriptive; no independence-based significance or missing-data correction
is applied. Undocumented member fields are not interpreted.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyze_s23b_member_outer_icl_relation import (
    EXPECTED_HASHES, SELECTED, MEMBER_BINS, MEMBER_DIAGNOSTICS, CLUSTER_RELATION,
    STACK_SUMMARY, read_json, read_tsv, sha256, write_tsv,
)

DEFAULT_OUTPUT = SELECTED.parent.parent / 's23b_member_profile_retention_2026-11-17_01'
BINS = ('5-10', '10-15', '15-20', '20-30', '30-40', '40-50', '50-75')
VARIABLES = {
    'redshift': 'z',
    'log10_richness': 'richness',
    'icl_contrast': 'icl_contrast_50_300_minus_300_800_njy_arcsec2',
    'candidate_count': 'n_member_positions_50_500_kpc_in_frame',
}


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f'refusing to overwrite {output}')
    paths = (SELECTED, MEMBER_BINS, MEMBER_DIAGNOSTICS, CLUSTER_RELATION, STACK_SUMMARY)
    hashes = {str(p): sha256(p) for p in paths}
    if any(hashes[str(p)] != EXPECTED_HASHES[str(p)] for p in paths):
        raise ValueError('pinned upstream inputs changed')
    if read_json(STACK_SUMMARY).get('science_measurement') is not False:
        raise ValueError('upstream science status changed')
    frames, bins, diagnostics, clusters = map(read_tsv, paths[:4])
    frame_map = {(r['target_id'], r['filter']): r for r in frames}
    if len(frames) != 979 or len(frame_map) != len(frames):
        raise ValueError('expected 979 unique selected target/filter frames')
    diag_map = {}
    by_frame = defaultdict(list)
    for r in diagnostics:
        key = (r['target_id'], r['filter'], r['member_index'])
        if key in diag_map or key[:2] not in frame_map:
            raise ValueError('duplicate or unknown member diagnostic')
        if r['status'] not in ('no_unique_sep_match_within_1arcsec', 'profiled', 'no_bin_passed_coverage'):
            raise ValueError('unknown diagnostic status')
        diag_map[key] = r
        by_frame[key[:2]].append(r)
    observed = Counter()
    counts = Counter()
    member_bin_keys = set()
    for r in bins:
        key = (r['target_id'], r['filter'], r['member_index'])
        bk = key + (r['radial_bin_kpc'],)
        if bk in member_bin_keys or key not in diag_map or bk[-1] not in BINS:
            raise ValueError('duplicate or unrecognized profile bin')
        if diag_map[key]['status'] != 'profiled' or r['source_row'] != diag_map[key]['source_row']:
            raise ValueError('profile/diagnostic join changed')
        if min(float(r['target_coverage']), float(r['control_coverage_mean'])) < .6 or int(r['n_valid_controls']) < 3:
            raise ValueError('retained profile fails original coverage/control gate')
        member_bin_keys.add(bk)
        observed[key] += 1
        counts[(key[0], key[1], bk[-1])] += 1
    for key, r in diag_map.items():
        if observed[key] != int(r['profile_bins_retained']):
            raise ValueError('diagnostic retained-bin count mismatch')
    cluster_map = {(r['target_id'], r['filter'], r['radial_bin_kpc']): r for r in clusters}
    if len(cluster_map) != len(clusters) or set(cluster_map) != set(counts):
        raise ValueError('saved cluster rows disagree with member bin support')
    if any(int(cluster_map[k]['n_member_positions']) != n for k, n in counts.items()):
        raise ValueError('cluster member count mismatch')

    ledger = []
    for key, frame in sorted(frame_map.items()):
        ds = by_frame[key]
        n = int(frame['n_member_positions_50_500_kpc_in_frame'])
        if len(ds) != n:
            raise ValueError('selected frame candidate count mismatch')
        matched = sum(r['status'] != 'no_unique_sep_match_within_1arcsec' for r in ds)
        for radial in BINS:
            retained = counts[key + (radial,)]
            ledger.append(dict(target_id=key[0], filter=key[1], radial_bin_kpc=radial,
                n_candidates=n, n_unmatched=n-matched, n_matched=matched,
                n_matched_without_bin=matched-retained, n_retained=retained,
                cluster_retained=retained > 0, candidate_retention_fraction=retained/n,
                redshift=float(frame['z']), richness=float(frame['richness']),
                icl_contrast=float(frame[VARIABLES['icl_contrast']]), science_measurement=False))
            if not 0 <= retained <= matched <= n:
                raise ValueError('attrition partition is inconsistent')

    radial_rows, balance_rows, quartile_rows = [], [], []
    for band in ('i', 'i2'):
        fs = [r for r in frames if r['filter'] == band]
        for radial in BINS:
            ls = [r for r in ledger if r['filter'] == band and r['radial_bin_kpc'] == radial]
            total = sum(r['n_candidates'] for r in ls)
            kept = sum(r['n_retained'] for r in ls)
            nr = sum(r['cluster_retained'] for r in ls)
            radial_rows.append(dict(filter=band, radial_bin_kpc=radial, n_selected_clusters=len(ls),
                n_retained_clusters=nr, cluster_retention_fraction=nr/len(ls),
                n_candidates=total, n_unmatched=sum(r['n_unmatched'] for r in ls),
                n_matched_without_bin=sum(r['n_matched_without_bin'] for r in ls),
                n_retained_candidates=kept, candidate_retention_fraction=kept/total))
        outer = {r['target_id']: r for r in ledger if r['filter'] == band and r['radial_bin_kpc'] == '50-75'}
        keep = np.array([outer[r['target_id']]['cluster_retained'] for r in fs], dtype=bool)
        for label, field in VARIABLES.items():
            x = np.array([float(r[field]) for r in fs])
            if label == 'log10_richness':
                x = np.log10(x)
            if not np.all(np.isfinite(x)):
                raise ValueError('nonfinite balance covariate')
            retained, lost = x[keep], x[~keep]
            scale = float(np.sqrt((np.var(retained, ddof=1) + np.var(lost, ddof=1))/2))
            balance_rows.append(dict(filter=band, variable=label, n_retained=len(retained),
                n_lost=len(lost), retained_mean=float(np.mean(retained)), lost_mean=float(np.mean(lost)),
                retained_median=float(np.median(retained)), lost_median=float(np.median(lost)),
                standardized_mean_difference=float((np.mean(retained)-np.mean(lost))/scale) if scale > 0 else None))
            edges = np.quantile(x, [.25, .5, .75])
            strata = np.searchsorted(edges, x, side='right')
            for q in range(4):
                group = strata == q
                size = int(np.sum(group))
                nret = int(np.sum(group & keep))
                quartile_rows.append(dict(filter=band, variable=label, quartile=q+1,
                    lower_bound=None if q == 0 else float(edges[q-1]),
                    upper_bound=None if q == 3 else float(edges[q]),
                    n_selected=size, n_retained=nret,
                    retention_fraction=nret/size if size else None))

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name + '.partial.', dir=output.parent))
    products = {'cluster_radial_retention.tsv': ledger, 'radial_retention.tsv': radial_rows,
                'outer_retention_balance.tsv': balance_rows, 'outer_retention_quartiles.tsv': quartile_rows}
    for name, rows in products.items():
        write_tsv(tmp/name, rows, list(rows[0]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    midpoints = [7.5, 12.5, 17.5, 25, 35, 45, 62.5]
    for band, color in [('i', '#2864aa'), ('i2', '#c65a1d')]:
        rows = [r for r in radial_rows if r['filter'] == band]
        axes[0].plot(midpoints, [r['cluster_retention_fraction'] for r in rows], 'o-', label=band, color=color)
        rows = [r for r in balance_rows if r['filter'] == band]
        axes[1].plot([r['standardized_mean_difference'] for r in rows], np.arange(4), 'o', color=color, label=band)
    axes[0].set(xlabel='Member radius [proper kpc]', ylabel='Cluster retention', ylim=(0, 1.12))
    axes[0].text(.04, .95, 'Conditional on selected frames', transform=axes[0].transAxes, va='top', fontsize=10)
    axes[0].legend(loc='upper right', fontsize=10)
    axes[1].set_yticks(np.arange(4), ['Redshift', 'log richness', 'ICL contrast', 'Candidate count'])
    axes[1].axvline(0, color='gray', lw=.7)
    axes[1].set_xlabel('Retained − lost / pooled SD')
    axes[1].text(.04, .95, '50–75 kpc selection', transform=axes[1].transAxes, va='top', fontsize=10)
    axes[1].set_ylim(3.8, -.8)
    fig.tight_layout()
    fig.savefig(tmp/'retention.png', dpi=160)
    plt.close(fig)
    summary = dict(schema='astrafek-s23b-member-profile-retention-v1',
        status='complete_local_descriptive_retention_audit',
        created_utc=datetime.now(timezone.utc).isoformat(), n_selected_clusters=len(frames),
        n_candidate_positions=len(diagnostics), n_retained_candidate_bins=len(bins),
        input_sha256=hashes, script_sha256=sha256(Path(__file__)),
        output_sha256={p.name:sha256(p) for p in tmp.iterdir()},
        outer_bin=[r for r in radial_rows if r['radial_bin_kpc']=='50-75'],
        outer_balance=balance_rows, network_used=False, images_read=False,
        cas_query_submitted=False, source_rows_downloaded=False, science_measurement=False,
        limitations=['Conditional selected-frame/member-candidate support; not detection completeness.',
            'Matched-without-bin combines target-area and rotated-control failures; saved tables cannot separate them.',
            'Standardized differences and quartile fractions are descriptive, without spatial confidence intervals.',
            'No inverse selection weighting or missing-profile imputation.',
            'Membership purity, PSF, source Poisson, sky, coadd covariance and other science gates remain unresolved.'])
    (tmp/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    os.rename(tmp, output)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output.resolve())
    print(json.dumps({'status':result['status'], 'output':str(args.output.resolve()),
        'outer_bin':result['outer_bin'], 'outer_balance':result['outer_balance']}, indent=2))
