#!/usr/bin/env python3
"""Correct only the documented weight identifier; retain approved data scope.

Reads public schema documentation, never HSC catalog rows. Preserves the
rejected original plan and proves every SQL change is the same substitution.
"""
import argparse
import copy
import json
from pathlib import Path
import urllib.request

from prepare_s23b_s16a_shape_count_plan import canonical_sha
from prepare_s23b_s16a_shape_qa_pilot import OUTPUT as PREDECESSOR, DATA
from audit_s23b_s16a_wl_random_support import require, sha

PREDECESSOR_SHA = '7af70c2be81f83e61bc6f74caf98b66a6309ae1b0caed0b9b7d1a32f42bd3f10'
OLD = 'b.ishape_hsm_regauss_derived_weight'
NEW = 'b.ishape_hsm_regauss_derived_shape_weight'
OUTPUT = DATA/'wl_s23b_s16a_shape_qa_pilot_plan_2026-11-17_02'


def schema_column(document, name):
    columns = [dict(zip(document['headers'], row)) for row in document['columns']]
    found = [c for c in columns if c['name'].replace('(g|r|i|z|y)', 'i') == name]
    require(len(found) == 1, 'public schema does not identify requested column: '+name)
    return found[0]


def same_scope(original, amended):
    # Both entire documents are hash-pinned by the receiver; this independently
    # verifies semantic scope rather than treating a new SQL byte hash as consent.
    for key in original:
        if key not in {'column_expressions','output_sha256','batches','query_plan_sha256','before_any_submission'}:
            require(amended[key] == original[key], 'identifier correction changed scope: '+key)
    expressions = dict(original['column_expressions']);expressions['weight'] = NEW
    require(original['column_expressions']['weight'] == OLD and amended['column_expressions'] == expressions,
            'correction changed another column')
    require(len(amended['batches']) == len(original['batches']), 'batch count changed')
    for old, new in zip(original['batches'], amended['batches']):
        for key in ('name','center_ids','expected_raw_pairs'):
            require(old[key] == new[key], 'batch membership changed')


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite schema-fix plan')
    original = json.loads((PREDECESSOR/'shape_qa_plan.json').read_text())
    require(canonical_sha(original) == original['query_plan_sha256'] == PREDECESSOR_SHA, 'original approved plan changed')
    failure = json.loads((PREDECESSOR/'results/shape_qa_status.json').read_text())
    require(failure['confirmed_bytes_received'] == failure['budget_bytes_charged'] == 0
            and failure['query_submitted'] is False
            and all('job_id' not in j for j in failure['jobs'].values()), 'original attempt has transfers/jobs; review before amendment')
    require(failure['jobs']['s16a_shape_qa_000']['http_code'] == 406, 'unexpected failure to amend')
    schemas, source_bytes = {}, {}
    for table in ('meas','meas2','weaklensing_hsm_regauss'):
        url='https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.'+table+'.json'
        with urllib.request.urlopen(url, timeout=30) as f: raw=f.read(2*1024*1024)
        document=json.loads(raw)
        require(document['db'] == 'pdr1' and document['schema'] == 's16a_wide' and document['table'] == table, 'public schema identity differs')
        schemas[table]=document;source_bytes[table]=(url,raw)
    expressions=dict(original['column_expressions']);expressions['weight']=NEW
    aliases={'a':'meas2','b':'weaklensing_hsm_regauss','c':'meas'}
    evidence=[]
    for alias, expression in expressions.items():
        table_key, name=expression.split('.',1)
        col=schema_column(schemas[aliases[table_key]],name)
        evidence.append(dict(alias=alias, expression=expression, type=col['type'], unit=col['unit']))
    for table in schemas: schema_column(schemas[table],'object_id')
    schema_column(schemas['meas'],'icoord')
    require(schema_column(schemas['meas'],'ira')['unit'] == schema_column(schemas['meas'],'idec')['unit'] == 'degree', 'coordinate unit not confirmed')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('pilot_lenses.tsv','center_ledger.tsv'):
        require(sha(PREDECESSOR/name) == original['output_sha256'][name], 'original center ledger changed')
        (output/name).write_bytes((PREDECESSOR/name).read_bytes())
    schema_sources=[]
    for table,(url,raw) in source_bytes.items():
        file=output/('schema_'+table+'.json');file.write_bytes(raw)
        schema_sources.append(dict(url=url, filename=file.name, bytes=len(raw), sha256=sha(file)))
    (output/'schema_column_check.json').write_text(json.dumps(dict(
        source='Public schema-browser documentation, not a CAS metadata query or astronomical row transfer',
        checked_columns=evidence, source_documents=schema_sources),indent=2)+'\n')
    plan=copy.deepcopy(original);plan['column_expressions']=expressions
    for batch in plan['batches']:
        prior=(PREDECESSOR/batch['name']).read_text()
        require(prior.count(OLD) == len(batch['center_ids']), 'original SQL has unexpected weight references')
        sql=prior.replace(OLD,NEW)
        require(sql.replace(NEW,OLD) == prior, 'SQL differs beyond identifier substitution')
        file=output/batch['name'];file.write_text(sql)
        batch.update(sha256=sha(file),sql_bytes=file.stat().st_size)
    plan['before_any_submission'][0]='User approval of predecessor scope '+PREDECESSOR_SHA+' covers this strictly verified weight-identifier correction only; any data-scope change requires separate approval.'
    plan['identifier_correction']=dict(predecessor_plan_sha256=PREDECESSOR_SHA,
        predecessor_directory=str(PREDECESSOR), approved_predecessor_csv_bytes_received=0,
        output_alias='weight', original_expression=OLD, corrected_expression=NEW,
        semantic_data_scope_changed=False, explanation='Official overview used an absent shorthand; public schema names the same inverse-variance shape weight explicitly.',
        source_documents=schema_sources)
    plan['schema_fix_script_sha256']=sha(Path(__file__))
    plan['output_sha256']={f.name:sha(f) for f in sorted(output.iterdir()) if f.is_file()}
    same_scope(original,plan)
    plan['query_plan_sha256']=canonical_sha(plan)
    (output/'shape_qa_plan.json').write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    (output/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  shape_qa_plan.json\n')
    print(json.dumps(dict(status='prepared_schema_identifier_fix', executed_plan_sha256=plan['query_plan_sha256'],
        authorization_scope_plan_sha256=PREDECESSOR_SHA, changed_expression={OLD:NEW},
        centers=43, columns=12, expected_rows=169535, maximum_csv_bytes=134217728),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    build(parser.parse_args().output)
