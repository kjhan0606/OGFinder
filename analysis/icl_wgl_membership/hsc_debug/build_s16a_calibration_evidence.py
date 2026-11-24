#!/usr/bin/env python3
"""Validate saved S16A/COSMOS diagnostics and write an evidence checkpoint."""
import argparse
import hashlib
import json
from pathlib import Path

BASE = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
PRODUCTS = {
    'intake': 's16a_ephor_ab_cosmos_intake_2026-11-19_01',
    'response': 's16a_cosmos_lens_response_2026-11-19_01',
    'support': 's16a_pointz_background_support_2026-11-18_01',
    'radial': 's16a_radial_gap_comparison_2026-11-20_01',
    'extremes': 's16a_radial_gap_extremes_2026-11-20_01',
    'distribution': 's16a_cosmos_photoz_distribution_2026-11-20_01',
    'paired': 's16a_cosmos_paired_distances_2026-11-20_01',
}
LATEST_PRODUCTS = {
    'precision': 's16a_cosmos_distance_precision_2026-11-20_01',
    'bin_support': 's16a_cosmos_bin_support_2026-11-21_01',
    'influence': 's16a_cosmos_object_influence_2026-11-21_01',
    'covariance': 's16a_cosmos_bootstrap_covariance_2026-11-21_01',
    'tract': 's16a_cosmos_tract_common_support_2026-11-21_01',
    'tract_bootstrap': 's16a_cosmos_tract_bootstrap_2026-11-21_01',
}
SELECTION_PRODUCTS = {
    'missingness': 's16a_selection_missingness_20261121_01',
    'matched_missingness': 's16a_missingness_matched_centers_20261124_01',
    'selection_fields': 's16a_cosmos_selection_fields_20261124_01',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selection_evidence(root, intake_path, pilot):
    """Require all selected diagnostics and verify their cross-product lineage."""
    payloads, paths = {}, {}
    for name, directory in SELECTION_PRODUCTS.items():
        path = root/directory/'summary.json'
        payload = json.loads(path.read_text())
        if payload.get('science_measurement') is not False or payload.get('calibration_applied') is not False:
            raise ValueError('Unexpected selection science status: '+name)
        outputs = payload.get('output_sha256')
        if not isinstance(outputs, dict) or not outputs:
            raise ValueError('Missing selection artifact manifest: '+name)
        for filename, expected in outputs.items():
            if Path(filename).name != filename:
                raise ValueError('Selection artifact must be a local basename')
            if sha(path.parent/filename) != expected:
                raise ValueError('Changed selection artifact: '+filename)
        paths[name], payloads[name] = path, payload
    missing, matched, fields = (payloads[k] for k in SELECTION_PRODUCTS)
    if matched['input_summary_sha256'] != sha(paths['missingness']):
        raise ValueError('Matched missingness summary lineage mismatch')
    if matched['input_table_sha256'] != missing['output_sha256']['center_missingness.tsv']:
        raise ValueError('Matched missingness table lineage mismatch')
    if fields['intake_summary_sha256'] != sha(intake_path):
        raise ValueError('Selection field intake lineage mismatch')
    for key in ('plan_sha256', 'receipt_sha256', 'object_qa_sha256'):
        if missing['provenance'][key] != pilot[key]:
            raise ValueError('Selection pilot lineage mismatch: '+key)
    return dict(status='verified',inputs={name:dict(path=str(path),sha256=sha(path)) for name,path in paths.items()},
        selected_objects=missing['selected_objects'],missing_objects=missing['source_object_counts']['missing']['n'],
        n_lenses=matched['n_lenses'],all_five_randoms_nonempty=matched['all_five_supported'],
        n_literal_quality_fields_present=fields['n_literal_names_present'],
        n_required_quality_fields=fields['n_required_fields'],
        limitations=missing['limitations']+matched['limitations']+fields['limitations'])


def run(output, data_root=BASE, selection_root=None):
    if output.exists():
        raise FileExistsError(output)
    latest_present = {name for name, directory in LATEST_PRODUCTS.items()
                      if (data_root/directory/'summary.json').is_file()}
    if latest_present and latest_present != set(LATEST_PRODUCTS):
        raise ValueError('Incomplete additional S16A diagnostics')
    products = {**PRODUCTS, **(LATEST_PRODUCTS if latest_present else {})}
    summaries, paths = {}, {}
    for name, directory in products.items():
        path = data_root/directory/'summary.json'
        payload = json.loads(path.read_text())
        if payload.get('science_measurement') is not False:
            raise ValueError(f'Unexpected science status: {name}')
        for filename, expected in payload.get('output_sha256', {}).items():
            artifact = path.parent/filename
            if sha(artifact) != expected:
                raise ValueError(f'Changed output: {artifact}')
        paths[name] = path
        summaries[name] = payload
    checks = [
        ('response', 'cosmos_intake_sha256', 'intake'),
        ('distribution', 'cosmos_intake_sha256', 'intake'),
        ('radial', 'input_summary_sha256', 'support'),
        ('extremes', 'input_summary_sha256', 'radial'),
        ('paired', 'input_summary_sha256', 'distribution'),
    ]
    if latest_present:
        checks.extend([
            ('precision', 'intake_summary_sha256', 'intake'),
            ('precision', 'response_summary_sha256', 'response'),
            ('bin_support', 'intake_sha256', 'intake'),
            ('influence', 'input_intake_sha256', 'intake'),
            ('influence', 'input_support_sha256', 'bin_support'),
            ('covariance', 'input_intake_sha256', 'intake'),
            ('covariance', 'input_support_sha256', 'bin_support'),
            ('tract', 'input_intake_sha256', 'intake'),
            ('tract', 'input_support_sha256', 'bin_support'),
            ('tract_bootstrap', 'input_intake_sha256', 'intake'),
            ('tract_bootstrap', 'input_support_sha256', 'bin_support'),
            ('tract_bootstrap', 'input_prior_sha256', 'tract'),
        ])
    for child, field, parent in checks:
        if summaries[child][field] != sha(paths[parent]):
            raise ValueError(f'Lineage mismatch: {child} and {parent}')
    for key in ('plan_sha256', 'receipt_sha256', 'object_qa_sha256'):
        if summaries['response'][key] != summaries['distribution'][key] or summaries['response'][key] != summaries['support'][key]:
            raise ValueError(f'Pilot lineage mismatch: {key}')
    if summaries['response']['calibration_applied'] or summaries['response']['compatibility'] != 'partial_only':
        raise ValueError('Unexpected calibration state')
    selection = selection_evidence(selection_root, paths['intake'], summaries['response']) if selection_root else dict(status='not_requested')
    report = dict(schema='astrafek-s16a-calibration-evidence-v1', science_measurement=False,
                  additional_diagnostics='verified' if latest_present else 'not_available',
                  calibration_applied=False, inputs={name:dict(path=str(path),sha256=sha(path))
                                                     for name,path in paths.items()},
                  runner_sha256=sha(Path(__file__)),
                  findings=dict(reference_joined=summaries['intake']['joined_rows'],
                                reference_total=summaries['intake']['reference_rows'],
                                lenses=summaries['response']['n_lenses'],
                                matched_comparable=summaries['paired']['n_comparable'],
                                matched_lens_above_random=summaries['paired']['lens_above_random_median']),
                  selection_diagnostics=selection,
                  limitations=summaries['response']['limitations']+[summaries['paired']['limitation']])
    output.mkdir(parents=True)
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['findings'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', '--output-dir', dest='output', type=Path, required=True)
    parser.add_argument('--data-root', type=Path, default=BASE)
    parser.add_argument('--selection-root', default='', help='Optional directory holding all three saved selection diagnostics')
    args = parser.parse_args()
    run(args.output.resolve(), args.data_root.resolve(), Path(args.selection_root).resolve() if args.selection_root else None)
