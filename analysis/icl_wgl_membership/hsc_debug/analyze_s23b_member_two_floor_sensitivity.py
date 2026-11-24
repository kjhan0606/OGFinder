#!/usr/bin/env python3
"""Exploratory 0.4/0.6 common-support comparison after empty 0.8 i support.

Reuses the pinned three-floor estimator, preserving the 0.6 reference model
and connected-block bootstrap. The two-floor choice is post hoc and its
intervals are descriptive, not a new detection test.
"""
import argparse
import functools
import hashlib
import json
from pathlib import Path
import analyze_s23b_member_coverage_sensitivity as method


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--replicas', type=int, default=5000)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error('at least 100 bootstrap draws required')
    method.FLOORS = (.4, .6)
    method.candidate_profiles = functools.partial(method.candidate_profiles, floors=method.FLOORS)
    result = method.run(args.run_root.resolve(), args.output.resolve(), args.replicas)
    metadata = {'schema': 'astrafek-member-two-floor-invocation-v1',
                'selection': 'post hoc: three-floor i intersection empty',
                'floors': list(method.FLOORS), 'science_measurement': False,
                'wrapper_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'method_sha256': result['script_sha256']}
    with (args.output / 'two_floor_invocation.json').open('x') as f:
        json.dump(metadata, f, indent=2)
        f.write('\n')
    print(json.dumps({'status': result['status'], 'primary_layout': result['primary_layout'],
                      'primary_paired_changes': result['primary_paired_changes']}, indent=2))
