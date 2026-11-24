#!/usr/bin/env python3
"""Measure the saved paired control radial residual and joint block variance."""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_s23b_icl_relative_radial_contrasts import (
    CONTROL_MANIFEST, EXPECTED_CONTROL_MANIFEST_SHA256, EXPECTED_PAIR_PROFILE_SHA256,
    EXPECTED_SPATIAL_SUMMARY_SHA256, EXPECTED_TARGET_MANIFEST_SHA256, PAIRED_PROFILES,
    PHASES, SPATIAL_SUMMARY, TARGET_MANIFEST, WIDTHS_DEG, broad_profile,
    connected_pair_blocks, make_pair_contrasts, read_json, read_tsv, sha256, write_tsv,
)


def contrast(profile, role):
    bins = [dict(item, delta_mean_njy_arcsec2=item[f'{role}_mean_njy_arcsec2']) for item in profile]
    inner, n_inner = broad_profile(bins, 50, 300)
    outer, n_outer = broad_profile(bins, 300, 800)
    if not n_inner or not n_outer:
        return None
    return inner - outer


def run(output):
    if output.exists():
        raise FileExistsError(output)
    inputs = {TARGET_MANIFEST: EXPECTED_TARGET_MANIFEST_SHA256,
              CONTROL_MANIFEST: EXPECTED_CONTROL_MANIFEST_SHA256,
              PAIRED_PROFILES: EXPECTED_PAIR_PROFILE_SHA256,
              SPATIAL_SUMMARY: EXPECTED_SPATIAL_SUMMARY_SHA256}
    for path, digest in inputs.items():
        if sha256(path) != digest:
            raise ValueError(f'Input hash mismatch: {path}')
    spatial = read_json(SPATIAL_SUMMARY)
    targets = {str(row['id']): row for row in read_json(TARGET_MANIFEST)['clusters']}
    controls = {str(row['target_cluster_id']): row for row in read_json(CONTROL_MANIFEST)['clusters']}
    pair_rows = read_tsv(PAIRED_PROFILES)
    by_band = make_pair_contrasts(pair_rows, targets, controls)
    profiles = {str(row['target_id']): json.loads(row['profile_json']) for row in pair_rows}
    if len(profiles) != len(pair_rows):
        raise ValueError('Expected one matched profile per target')
    summaries, pairs = [], []
    for band, rows in by_band.items():
        valid_rows = []
        values = []
        for row in rows:
            ident = row['target_id']
            target, control = (contrast(profiles[ident], role) for role in ('target', 'control'))
            delta = row['contrasts']['50-300_minus_300-800']['contrast_njy_arcsec2']
            if delta is not None:
                if target is None or control is None or not np.isclose(target-control, delta, rtol=1e-12, atol=1e-9):
                    raise ValueError(f'Contrast identity fails for {ident}')
                valid_rows.append(row)
                values.append((target, control))
            pairs.append(dict(target_id=ident, control_id=row['control_id'], filter=band,
                              target_contrast=target, control_contrast=control, delta_contrast=delta))
        values = np.asarray(values)
        rng = np.random.default_rng(2026111801 if band == 'i' else 2026111802)
        # Derive blocks from the full sample to preserve connections through any unsupported pair.
        supported = {row['target_id'] for row in valid_rows}
        for width in WIDTHS_DEG:
            for phase in PHASES:
                labels, _ = connected_pair_blocks(rows, controls, width, phase, spatial['ra_grid_origin_deg'])
                kept = labels[[row['target_id'] in supported for row in rows]]
                groups, group_index = np.unique(kept, return_inverse=True)
                if len(groups) < 3:
                    raise ValueError(f'Insufficient blocks for {band}/{width}/{phase}')
                counts = np.bincount(group_index)
                sums = np.column_stack([np.bincount(group_index, weights=values[:, i]) for i in range(2)])
                draws = rng.multinomial(len(groups), np.full(len(groups), 1/len(groups)), size=20000)
                means = (draws @ sums) / (draws @ counts)[:, None]
                paired = means[:, 0] - means[:, 1]
                covariance = np.cov(means, rowvar=False)
                variance = paired.var(ddof=1)
                independent = np.trace(covariance)
                if not np.isclose(variance, independent-2*covariance[0, 1]):
                    raise ValueError('Bootstrap contrast variance does not reconcile')
                ci_control = np.quantile(means[:, 1], [.025, .975])
                ci_delta = np.quantile(paired, [.025, .975])
                summaries.append(dict(filter=band, block_width_deg=width, grid_phase=phase,
                    supported_pairs=len(values), blocks=len(groups), largest_block_pairs=int(max(counts)),
                    target_mean=float(values[:, 0].mean()), control_mean=float(values[:, 1].mean()),
                    control_ci_low=float(ci_control[0]), control_ci_high=float(ci_control[1]),
                    delta_mean=float((values[:, 0]-values[:, 1]).mean()),
                    delta_ci_low=float(ci_delta[0]), delta_ci_high=float(ci_delta[1]),
                    target_control_covariance=float(covariance[0, 1]),
                    paired_variance=float(variance), independence_variance=float(independent),
                    paired_over_independence_sd=float(np.sqrt(variance/independent))))
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output/'pair_contrasts.tsv', pairs)
    write_tsv(output/'block_diagnostics.tsv', summaries)
    summary = dict(status='complete_saved_control_radial_residual_diagnostic', science_measurement=False,
                   contrast_kpc=[[50, 300], [300, 800]], surface_brightness_unit='nJy arcsec^-2',
                   bootstrap_replicas=20000, summaries=summaries,
                   limitation='Matched control residual includes sky, sources, masking and resampling; it does not isolate their causes.',
                   input_sha256={str(path): digest for path, digest in inputs.items()},
                   runner_sha256=sha256(Path(__file__)),
                   output_sha256={f.name: sha256(f) for f in output.iterdir()})
    with (output/'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps([r for r in summaries if r['block_width_deg'] == 5 and r['grid_phase'] == 0], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
