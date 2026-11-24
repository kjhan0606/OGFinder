#!/usr/bin/env python3
"""Summarize radial gap differences by source and control support."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics

BASE = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SOURCE = BASE / 's16a_radial_gap_comparison_2026-11-20_01'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((SOURCE / 'summary.json').read_text())
    path = SOURCE / 'matched_annulus_gap.tsv'
    if sha(path) != manifest['output_sha256'][path.name]:
        raise ValueError('Changed input')
    with path.open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t') if float(r['gap']) == 0.2]
    if len(rows) != 88:
        raise ValueError('Expected 88 lens-annulus rows')
    valid = [r for r in rows if r['lens_minus_random_retained']]
    conditions = {
        'all_evaluable': valid,
        'all_randoms_supported': [r for r in valid if int(r['random_supported_centers']) == 5],
        'all_randoms_and_20_lens_pairs': [r for r in valid if int(r['random_supported_centers']) == 5 and int(r['lens_base_pairs']) >= 20],
        'outside_inner_0p1_mpc': [r for r in valid if int(r['annulus']) > 0],
    }
    summaries = []
    for label, sample in conditions.items():
        values = [float(r['lens_minus_random_retained']) for r in sample]
        summaries.append(dict(group=label, n=len(values), median=statistics.median(values) if values else None,
                              minimum=min(values) if values else None, maximum=max(values) if values else None))
    worst = min(valid, key=lambda r: float(r['lens_minus_random_retained']))
    output.mkdir(parents=True)
    result = dict(science_measurement=False, input_sha256=sha(path), input_summary_sha256=sha(SOURCE/'summary.json'),
                  runner_sha256=sha(Path(__file__)), summaries=summaries, largest_negative=worst,
                  note='Descriptive post-hoc subsets. Pair counts include reused sources. No uncertainties or survey inference.')
    (output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())
