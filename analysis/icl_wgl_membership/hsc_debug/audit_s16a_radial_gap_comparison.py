#!/usr/bin/env python3
"""Compare matched lens and random point-z selection by annulus."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
INPUT = DATA / 's16a_pointz_background_support_2026-11-18_01'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_path = INPUT / 'summary.json'
    summary = json.loads(summary_path.read_text())
    table = INPUT / 'center_annulus_pointz_support.tsv'
    if digest(table) != summary['output_sha256'][table.name]:
        raise ValueError('Input table hash mismatch')
    with table.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    keys = {(r['center_id'], int(r['annulus']), float(r['gap'])): r for r in rows}
    if len(keys) != 1584 or len(rows) != len(keys):
        raise ValueError('Unexpected or duplicate center records')
    groups = {}
    for r in rows:
        groups.setdefault(r['lens_id'], {})[r['center_id']] = r['role']
    if len(groups) != 11 or any(list(g.values()).count('lens') != 1 or
                                 list(g.values()).count('random') != 5 for g in groups.values()):
        raise ValueError('Unexpected matched-center membership')
    result = []
    for lens, centers in sorted(groups.items()):
        for annulus in range(8):
            for gap in (0.1, 0.2):
                record = dict(lens_id=lens, annulus=annulus, gap=gap)
                for role in ('lens', 'random'):
                    ids = [c for c, value in centers.items() if value == role]
                    base = [keys[c, annulus, 0.] for c in ids]
                    selected = [keys[c, annulus, gap] for c in ids]
                    before = sum(float(r['pair_weight_sum']) for r in base)
                    after = sum(float(r['pair_weight_sum']) for r in selected)
                    record[role + '_base_weight'] = before
                    record[role + '_retained_fraction'] = after / before if before else ''
                    record[role + '_base_pairs'] = sum(int(r['candidate_pairs']) for r in base)
                    record[role + '_retained_pairs'] = sum(int(r['candidate_pairs']) for r in selected)
                    record[role + '_supported_centers'] = sum(float(r['pair_weight_sum']) > 0 for r in selected)
                if record['lens_retained_fraction'] != '' and record['random_retained_fraction'] != '':
                    record['lens_minus_random_retained'] = record['lens_retained_fraction'] - record['random_retained_fraction']
                else:
                    record['lens_minus_random_retained'] = ''
                result.append(record)
    output.mkdir(parents=True)
    target = output / 'matched_annulus_gap.tsv'
    with target.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=result[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(result)
    metadata = dict(science_measurement=False, input_summary_sha256=digest(summary_path),
                    input_table_sha256=digest(table), runner_sha256=digest(Path(__file__)),
                    output_sha256={target.name:digest(target)}, count=len(result),
                    method='Per lens and annulus, compare retained inverse-critical-density-squared '
                           'weighted fractions after the same z_best gap cut. Random fractions pool '
                           'the five matched random centers. Empty denominators are missing.',
                    limitation='Diagnostic retention fractions cannot prove selection or calibration '
                               'equivalence. Sources are reused, random coverage is incomplete, '
                               'and uncertainties and area corrections are not calculated.')
    (output/'summary.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
