#!/usr/bin/env python3
"""Describe paired missing-input fractions; empty random centers stay explicit."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

BASE = Path('/home/kjhan/BACKUP/Astrafek/s16a_selection_missingness_20261121_01')
EXPECTED = '755e1dbf6ff721da772d6765f96434efa318eb3f64d7d672141272f3a2cb4d9b'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(source, output):
    if output.exists():
        raise FileExistsError(output)
    summary_path = source/'summary.json'
    if digest(summary_path) != EXPECTED:
        raise ValueError('input summary differs from audited result')
    summary = json.loads(summary_path.read_text())
    table = source/'center_missingness.tsv'
    if digest(table) != summary['output_sha256'][table.name]:
        raise ValueError('center ledger digest differs')
    with table.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) != 66 or len({r['center_id'] for r in rows}) != 66:
        raise ValueError('expected 66 distinct centers')
    groups = defaultdict(list)
    for row in rows:
        n, missing = int(row['selected_pairs']), int(row['missing_pairs'])
        if not 0 <= missing <= n:
            raise ValueError('invalid counts')
        fraction = row['missing_pair_fraction']
        if n and abs(float(fraction)-missing/n) > 1e-12:
            raise ValueError('fraction disagrees with counts')
        if not n and fraction:
            raise ValueError('empty center must not have a fraction')
        groups[row['lens_id']].append(row)
    results = []
    for lens_id, group in sorted(groups.items()):
        lens = [r for r in group if r['role']=='lens']
        random = [r for r in group if r['role']=='random']
        if len(lens) != 1 or len(random) != 5:
            raise ValueError('expected one lens and five fixed randoms')
        lens = lens[0]
        n_l = int(lens['selected_pairs'])
        n_r = sum(int(r['selected_pairs']) for r in random)
        f_l = int(lens['missing_pairs'])/n_l if n_l else None
        f_r = sum(int(r['missing_pairs']) for r in random)/n_r if n_r else None
        fractions = [float(r['missing_pair_fraction']) for r in random if int(r['selected_pairs'])]
        results.append(dict(lens_id=lens_id,lens_pairs=n_l,random_pairs=n_r,
            nonempty_randoms=len(fractions),empty_randoms=5-len(fractions),
            lens_missing_fraction=f_l,pooled_random_missing_fraction=f_r,
            lens_minus_pooled_random=f_l-f_r if f_l is not None and f_r is not None else None,
            random_center_fraction_min=min(fractions) if fractions else None,
            random_center_fraction_max=max(fractions) if fractions else None,
            all_five_randoms_nonempty=len(fractions)==5))
    output.mkdir(parents=True)
    target = output/'matched_missingness.tsv'
    with target.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t')
        writer.writeheader();writer.writerows(results)
    differences=[r for r in results if r['lens_minus_pooled_random'] is not None]
    report=dict(science_measurement=False,calibration_applied=False,n_lenses=len(results),
        empty_randoms=sum(r['empty_randoms'] for r in results),
        all_five_supported=sum(r['all_five_randoms_nonempty'] for r in results),
        minimum_difference=min(differences,key=lambda r:r['lens_minus_pooled_random']),
        maximum_difference=max(differences,key=lambda r:r['lens_minus_pooled_random']),
        definition='Missing best/risk input pair counts divided by selected counts; pooled randoms use count denominators, not an equal-center average.',
        limitations=['Pairs reuse sources; no binomial error or significance claimed.',
            'Empty randoms retain undefined fractions and cannot validate matched spatial support.',
            'No lensing-efficiency weights or calibration correction can be inferred from these fractions.'],
        input_summary_sha256=EXPECTED,input_table_sha256=digest(table),
        runner_sha256=digest(Path(__file__)),output_sha256={target.name:digest(target)})
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=BASE)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.input.resolve(),args.output.resolve())
