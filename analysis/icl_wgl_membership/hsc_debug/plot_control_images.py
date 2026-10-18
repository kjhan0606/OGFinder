#!/usr/bin/env python3
"""Plot saved matched-field diagnostics without altering any image or table."""
import argparse
import csv
from pathlib import Path
import platform

from matplotlib import pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison', type=Path, required=True)
    a = p.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        p.error('Run on lageunha CPU')
    output = a.comparison/'paired_outer_diagnostics.png'
    with (a.comparison/'paired_stats.tsv').open() as f:
        paired = list(csv.DictReader(f, delimiter='\t'))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for cohort, color in [('area_qualified', '#286a9d'), ('mask_stress', '#c06b20')]:
        subset = [r for r in paired if r['cohort'] == cohort and float(r['kron_factor']) == 3.5
                  and r['control_status'] == r['cluster_status'] == 'diagnostic_only']
        axes[0].scatter([float(r['control_residual_ratio']) for r in subset],
                        [float(r['cluster_residual_ratio']) for r in subset], color=color, s=22, label=cohort)
        axes[1].scatter([float(r['control_outer_valid_fraction']) for r in subset],
                        [float(r['cluster_outer_valid_fraction']) for r in subset], color=color, s=22)
    axes[0].set_xlabel('Control tile residual MAD / nominal sigma')
    axes[0].set_ylabel('Cluster tile residual MAD / nominal sigma')
    axes[1].set_xlabel('Control outer retained fraction')
    axes[1].set_ylabel('Cluster outer retained fraction')
    for ax in axes:
        lo = min(ax.get_xlim()[0], ax.get_ylim()[0]); hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], 'k--', linewidth=.7)
    axes[0].legend(fontsize=8)
    fig.suptitle('Matched 600-900 kpc fields, Kron 3.5 masks\nDiagnostic only: selected pairs; residual sources and sky may remain', fontsize=10)
    fig.tight_layout()
    with output.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
