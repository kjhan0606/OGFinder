#!/usr/bin/env python3
"""Compare COSMOS CDF distances between lenses and their random centers."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import statistics

BASE = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SOURCE = BASE/'s16a_cosmos_photoz_distribution_2026-11-20_01'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_path = SOURCE/'summary.json'
    summary = json.loads(summary_path.read_text())
    source = SOURCE/'center_photoz_distribution.tsv'
    if digest(source) != summary['output_sha256'][source.name]:
        raise ValueError('Input digest differs')
    with source.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if len(rows) != summary['n_centers'] or len({r['center_id'] for r in rows}) != len(rows):
        raise ValueError('Unexpected center table')
    groups = defaultdict(list)
    for row in rows:
        groups[row['lens_id']].append(row)
    results = []
    for lens_id, members in sorted(groups.items()):
        lens = [r for r in members if r['role'] == 'lens']
        random = [r for r in members if r['role'] == 'random']
        if len(lens) != 1 or len(random) > 5:
            raise ValueError('Invalid group membership')
        lens = lens[0]
        values = [float(r['weighted_cdf_max_abs_difference']) for r in random]
        median = statistics.median(values) if values else None
        results.append(dict(lens_id=lens_id, lens_z=float(lens['lens_z']),
                            lens_distance=float(lens['weighted_cdf_max_abs_difference']),
                            lens_neff=float(lens['weighted_neff']), random_count=len(random),
                            random_median=median, random_min=min(values) if values else None,
                            random_max=max(values) if values else None,
                            lens_minus_random_median=float(lens['weighted_cdf_max_abs_difference'])-median if values else None))
    output.mkdir(parents=True)
    dest = output/'paired_cdf_distances.tsv'
    with dest.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    comparable = [r for r in results if r['random_count']]
    report = dict(science_measurement=False, n_lenses=len(results), n_comparable=len(comparable),
                  lens_above_random_median=sum(r['lens_minus_random_median']>0 for r in comparable),
                  median_difference=statistics.median(r['lens_minus_random_median'] for r in comparable),
                  n_all_five_randoms=sum(r['random_count']==5 for r in results),
                  input_summary_sha256=digest(summary_path), input_sha256=digest(source),
                  runner_sha256=digest(Path(__file__)), output_sha256={dest.name:digest(dest)},
                  limitation='Differences compare nonnegative CDF distance magnitudes, not signed CDF '
                             'curves. Random centers have varying sample sizes and missing coverage. '
                             'Reused sources prevent independent-sample inference.')
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
