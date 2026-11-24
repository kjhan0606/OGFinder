#!/usr/bin/env python3
"""Prepare an unsubmitted, bounded S16A shape/calibration QA pilot.

This creates only local SQL and provenance. It never submits a query, reads
source rows, downloads calibration data, or changes the scientific selection.
Only the twelve output columns documented in the plan may be requested.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM

from audit_s23b_s16a_wl_random_support import sha, require, read_tsv, write_table, TARGET_MANIFEST, DATA
from prepare_s23b_s16a_shape_count_plan import canonical_sha, sql_literal, sql_number

AUDIT = DATA/'s23b_s16a_wl_random_support_audit_2026-11-17_01'
AUDIT_SHA = '56c463d81b7c099dec6380aa1feb84efacd2625d6fbf16e8c189eec7719dea79'
OUTPUT = DATA/'wl_s23b_s16a_shape_qa_pilot_plan_2026-11-17_01'
BATCH_SIZE = 10
COLUMN_EXPRESSIONS = {
    'object_id': 'b.object_id',
    'source_ra': 'c.ira',
    'source_dec': 'c.idec',
    'e1': 'a.ishape_hsm_regauss_e1',
    'e2': 'a.ishape_hsm_regauss_e2',
    'sigma_e': 'b.ishape_hsm_regauss_derived_sigma_e',
    'rms_e': 'b.ishape_hsm_regauss_derived_rms_e',
    'weight': 'b.ishape_hsm_regauss_derived_weight',
    'm': 'b.ishape_hsm_regauss_derived_shear_bias_m',
    'c1': 'b.ishape_hsm_regauss_derived_shear_bias_c1',
    'c2': 'b.ishape_hsm_regauss_derived_shear_bias_c2',
}
OUTPUT_COLUMNS = ['center_id', *COLUMN_EXPRESSIONS]


def support_class(row):
    n = int(row['n_positive'])
    require(0 <= n <= 5, 'invalid positive-control count')
    return 'five_positive' if n == 5 else 'all_zero' if n == 0 else 'partial_positive'


def pilot_lenses(rows):
    require(len({r['matched_lens_id'] for r in rows}) == len(rows), 'duplicate input lens')
    require(all(r['icl_filter'] in ('i', 'i2') and int(r['target_raw_count']) >= 1000 for r in rows), 'invalid lens metadata')
    output = []
    for band in ('i', 'i2'):
        for label in ('five_positive', 'partial_positive', 'all_zero'):
            group = sorted((r for r in rows if r['icl_filter'] == band and support_class(r) == label),
                           key=lambda r: (int(r['target_raw_count']), r['matched_lens_id']))
            if not group:
                continue
            # At most two deterministic count-rank positions per stratum.
            indices = sorted({len(group)//4, min(len(group)-1, 3*len(group)//4)})
            for rank in indices:
                output.append(dict(group[rank], pilot_stratum=f'{band}:{label}',
                                   stratum_size=len(group), stratum_rank_zero_based=rank))
    return output


def make_centers(pilot, random_rows, targets):
    selected_ids = {r['matched_lens_id'] for r in pilot}
    require(len(selected_ids) == len(pilot), 'pilot repeats a lens')
    groups = defaultdict(list)
    for r in random_rows:
        if r['matched_lens_id'] in selected_ids:
            groups[r['matched_lens_id']].append(r)
    require(set(groups) == selected_ids, 'missing pilot random group')
    cosmology = FlatLambdaCDM(H0=70, Om0=.3, Tcmb0=0)
    result = []
    for row in pilot:
        lens = row['matched_lens_id']
        require(lens in targets and len(groups[lens]) == 5, 'invalid pilot lens/group')
        require({int(r['control_slot']) for r in groups[lens]} == {1, 2, 3, 4, 5}, 'invalid random slots')
        t = targets[lens]
        z = float(t['z_cl'])
        require(math.isfinite(z) and z > 0 and math.isclose(z, float(row['z_lens']), rel_tol=0, abs_tol=1e-12), 'pilot redshift differs')
        radius = math.degrees(3/float(cosmology.angular_diameter_distance(z).value))*3600
        common = dict(matched_lens_id=lens, icl_filter=row['icl_filter'],
                      pdr3_lens_field=row['lens_pdr3_field'], z_lens=z,
                      radius_mpc_proper=3., radius_arcsec=radius)
        result.append(dict(center_id='LC-'+lens, role='lens', ra_deg=float(t['ra_deg']), dec_deg=float(t['dec_deg']),
                           expected_raw_pairs=int(row['target_raw_count']), query_planned=True, **common))
        for r in sorted(groups[lens], key=lambda r: int(r['control_slot'])):
            count = int(r['random_raw_count'])
            require(count >= 0 and int(r['target_raw_count']) == int(row['target_raw_count']), 'pilot count mismatch')
            result.append(dict(center_id=r['random_center_id'], role='random', ra_deg=float(r['ra_deg']), dec_deg=float(r['dec_deg']),
                               expected_raw_pairs=count, query_planned=count > 0, **common))
    require(len({r['center_id'] for r in result}) == len(result), 'duplicate center ID')
    require(len({(round(r['ra_deg'], 12), round(r['dec_deg'], 12)) for r in result}) == len(result), 'duplicate center position')
    require(all(math.isfinite(r['ra_deg']) and 0 <= r['ra_deg'] < 360 and math.isfinite(r['dec_deg']) and -90 <= r['dec_deg'] <= 90 for r in result), 'invalid pilot coordinates')
    return result


def query(center):
    fields = [f"{sql_literal(center['center_id'])} AS center_id"]
    fields.extend(f'{expr} AS {alias}' for alias, expr in COLUMN_EXPRESSIONS.items())
    return ('SELECT DISTINCT '+', '.join(fields)+'\n'
            'FROM s16a_wide.meas2 AS a\n'
            'INNER JOIN s16a_wide.weaklensing_hsm_regauss AS b USING (object_id)\n'
            'INNER JOIN s16a_wide.meas AS c USING (object_id)\n'
            f"WHERE coneSearch(c.icoord, {sql_number(center['ra_deg'])}, {sql_number(center['dec_deg'])}, {sql_number(center['radius_arcsec'])})")


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite pilot plan')
    require(sha(AUDIT/'summary.json') == AUDIT_SHA, 'support audit summary changed')
    audit = json.loads((AUDIT/'summary.json').read_text())
    require(audit['status'] == 'complete_local_count_support_diagnosis'
            and audit['science_measurement'] is False and audit['original_selection_modified'] is False, 'support audit is inconsistent')
    for name in ('per_lens_support.tsv', 'per_random_center.tsv'):
        require(sha(AUDIT/name) == audit['output_sha256'][name], 'support ledger changed')
    require(sha(TARGET_MANIFEST) == audit['input_sha256'][str(TARGET_MANIFEST)], 'target manifest changed')
    lenses = read_tsv(AUDIT/'per_lens_support.tsv')
    require(len(lenses) == 114, 'original lens count differs')
    pilot = pilot_lenses(lenses)
    require(len(pilot) == 11, 'reviewed pilot size differs')
    targets = {str(r['id']): r for r in json.loads(TARGET_MANIFEST.read_text())['clusters']}
    centers = make_centers(pilot, read_tsv(AUDIT/'per_random_center.tsv'), targets)
    positive = [r for r in centers if r['query_planned']]
    require(len(centers) == 66 and len(positive) == 43, 'reviewed pilot center scope differs')
    require({r['lens_pdr3_field'] for r in pilot} == {'AUTUMN', 'SPRING', 'HECTOMAP'}, 'pilot field coverage differs')
    expected_pairs = sum(r['expected_raw_pairs'] for r in positive)
    output.mkdir(parents=True, exist_ok=False)
    write_table(output/'pilot_lenses.tsv', pilot)
    write_table(output/'center_ledger.tsv', centers)
    batches = []
    for start in range(0, len(positive), BATCH_SIZE):
        group = positive[start:start+BATCH_SIZE]
        sql = '\nUNION ALL\n'.join(query(r) for r in group)+';\n'
        path = output/f's16a_shape_qa_{start//BATCH_SIZE:03d}.sql'
        path.write_text(sql)
        batches.append(dict(name=path.name, sha256=sha(path), sql_bytes=path.stat().st_size,
                            center_ids=[r['center_id'] for r in group],
                            expected_raw_pairs=sum(r['expected_raw_pairs'] for r in group)))
    estimate = [expected_pairs*160, expected_pairs*360]
    budget = 2**math.ceil(math.log2(estimate[1]*1.3))
    plan = dict(schema='astrafek-s16a-diagnostic-shape-qa-pilot-plan-v1',
        status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, network_used=False, source_rows_downloaded=False,
        images_returned=False, psfs_returned=False, photoz_returned=False, science_measurement=False,
        purpose='Diagnostic catalog/schema/calibration-value and radial raw-support QA only; no shear, mass, boost, n(z) or ICL-WL correlation measurement.',
        cas_release_version='pdr1', source_catalog='HSC S16A Wide', radius_mpc_proper=3.,
        cosmology=dict(H0_km_s_mpc=70., Omega_m=.3),
        n_lenses=11, n_random_centers_in_ledger=55, n_centers_in_ledger=66,
        n_positive_centers_to_query=43, n_known_zero_centers_not_queried=23,
        n_query_batches=len(batches), batch_size=BATCH_SIZE,
        expected_source_pairs_from_prior_counts=expected_pairs,
        estimated_uncompressed_csv_bytes_range=estimate,
        estimated_bytes_per_pair_range=[160, 360], estimated_bytes_are_not_a_guarantee=True,
        maximum_total_received_bytes=budget,
        output_columns=OUTPUT_COLUMNS, column_expressions=COLUMN_EXPRESSIONS,
        position_contract='c.ira/c.idec are the documented S16A sky-coordinate columns. Verify their units against known query cones before interpreting radial coordinates; do not silently relabel or convert.',
        query_tables=['s16a_wide.meas2', 's16a_wide.weaklensing_hsm_regauss', 's16a_wide.meas'],
        pilot_selection='For each ICL band i/i2 and raw support class five_positive/partial_positive/all_zero, sort by (target count, lens ID), take ranks floor(n/4) and floor(3n/4), deduplicating singleton ranks. Diagnostic sample only; original 114-lens/570-random ledger is unchanged.',
        centers=centers, batches=batches,
        input_sha256={str(AUDIT/'summary.json'):AUDIT_SHA, str(TARGET_MANIFEST):sha(TARGET_MANIFEST)},
        script_sha256=sha(Path(__file__)),
        output_sha256={f.name: sha(f) for f in sorted(output.iterdir()) if f.is_file()},
        primary_schema_source='https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/',
        new_sql_server_validated=False,
        before_any_submission=['Obtain explicit approval for this exact plan hash and its row/byte/column scope.',
            'Use a guarded serial submitter; verify all plan/SQL/input hashes, never persist credentials and preserve resumable job IDs.',
            'Stop on schema/HTTP/query errors or byte-budget exhaustion; do not silently amend the SQL scope.'],
        intake_checks=['Expected center IDs and unique (center_id, object_id) pairs versus the pinned count ledger; preserve explicit zero centers without SQL.',
            'Column names/types, finite calibration fields, shape/weight distributions, coordinate-unit/cone consistency and per-radius support.',
            'Cross-cone source reuse and consistent source properties; no S19A/Y3 source cuts, masks, z-bins or n(z) imposed.'],
        limitations=['No source redshift, photo-z PDF, image, PSF or auxiliary selection/PSF-shape fields are requested.',
            'These twelve fields do not independently reproduce all Y1 catalog flag/selection cuts or selection-bias corrections.',
            'The existing S19A/Y3 estimator and redshift calibration are not validated for S16A.',
            'Raw-control availability strata are deliberately selected for diagnosis, not for a representative or science sample.'])
    plan['query_plan_sha256'] = canonical_sha(plan)
    (output/'shape_qa_plan.json').write_text(json.dumps(plan, indent=2, allow_nan=False)+'\n')
    (output/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  shape_qa_plan.json\n')
    print(json.dumps({k:plan[k] for k in ('status','query_plan_sha256','n_lenses',
        'n_positive_centers_to_query','n_known_zero_centers_not_queried','n_query_batches',
        'expected_source_pairs_from_prior_counts','estimated_uncompressed_csv_bytes_range',
        'maximum_total_received_bytes','output_columns')}, indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    build(parser.parse_args().output)
