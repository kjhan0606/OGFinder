import csv
import io
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_random_id_supplement as p
import submit_s16a_random_id_supplement as s

needs_plan = pytest.mark.skipif(not all((p.output(k)/'plan.json').is_file() for k in (1, 2, 3)), reason='pinned plans are local data')


def test_parts_split_sorted_ids_contiguously():
    ids = list(range(10, 17))
    parts = [p.part_ids(k, ids) for k in (1, 2, 3)]
    assert parts == [['10', '11', '12'], ['13', '14', '15'], ['16']]


def test_kinds_use_pilot_validators_and_queries():
    assert set(s.VALIDATORS) == set(p.KINDS)
    for kind, (query, columns, _) in p.KINDS.items():
        sql = query(['5', '7'])
        assert 'FROM (VALUES\n  (5::bigint),\n  (7::bigint)\n) AS r(object_id)' in sql and 'LEFT JOIN' in sql and 'LIMIT' not in sql


def test_flag_validator_requires_every_id(tmp_path):
    cols = p.KINDS['flags'][1]
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(cols)
    w.writerow(['5' if c.endswith('object_id') else 'f' for c in cols])
    f = tmp_path/'f.csv'; f.write_text(stream.getvalue())
    assert s.VALIDATORS['flags'](f, ['5'])['rows'] == 1
    with pytest.raises(ValueError):
        s.VALIDATORS['flags'](f, ['5', '7'])


@needs_plan
def test_plans_verify_and_cover_all_ids_once():
    ids = p.all_ids()
    seen = []
    for k in (1, 2, 3):
        checked = p.validate(p.output(k), k)
        assert checked['batches'] == 75 and checked['estimate']['csv_bytes'] < p.CAP
        assert s.APPROVED_SHA[k] in (None, checked['plan_sha256'])
        seen += p.part_ids(k, ids)
    assert seen == [str(i) for i in ids] and len(set(seen)) == p.N_OBJECTS


@needs_plan
def test_submit_refuses_without_approval():
    had = (p.output(1)/'results').exists()
    r = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_random_id_supplement.py'), '--submit', '--part', '1'],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert r.returncode == 1
    assert ('not the exact approved' if s.APPROVED_SHA[1] is None else 'explicit user approval') in r.stderr
    assert (p.output(1)/'results').exists() == had
