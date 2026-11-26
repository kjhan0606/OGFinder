import csv
import io
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_local_background_rows as p
import submit_s16a_local_background_rows as s

needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


def payload(rows):
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(list(p.EXPRESSIONS))
    for ident, joined in rows:
        vals = {'meas_object_id': joined, 'forced_object_id': ident}
        for k in list(p.EXPRESSIONS)[3:]:
            vals[k] = 'f' if k in s.BOOLEANS else ('3' if k in s.INTEGERS else '0.5')
        w.writerow([ident] + [vals[k] for k in list(p.EXPRESSIONS)[1:]])
    return stream.getvalue()


def test_query_form():
    sql = p.query(['5', '7'])
    assert 'FROM (VALUES\n  (5::bigint),\n  (7::bigint)\n) AS r(object_id)' in sql
    assert 'LEFT JOIN s16a_wide.meas AS m ON m.object_id = r.object_id' in sql and 'LIMIT' not in sql
    assert 'LEFT JOIN s16a_wide.forced AS f ON f.object_id = r.object_id' in sql
    assert 'f.ivariance AS i_variance' in sql and 'm.ideblend_masked AS i_deblend_masked' in sql


def test_validator(tmp_path):
    f = tmp_path/'a.csv'
    f.write_text(payload([('5', '5'), ('7', '')]))
    assert s.validate_csv(f, ['5', '7'])['null_by_column']['meas_object_id'] == 1
    for rows, ids in (([('5', '5')], ['5', '7']), ([('5', '6'), ('7', '7')], ['5', '7']), ([('5', '5'), ('5', '5')], ['5'])):
        f.write_text(payload(rows))
        with pytest.raises(ValueError):
            s.validate_csv(f, ids)
    f.write_text(payload([('5', '5')]).replace(',f', ',maybe', 1))
    with pytest.raises(ValueError):
        s.validate_csv(f, ['5'])


@needs_plan
def test_plan_verifies():
    checked = p.validate()
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['estimated_csv_bytes'] < p.CAP


@needs_plan
def test_submit_refuses_without_approval():
    had = (p.OUTPUT/'results').exists()
    r = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_local_background_rows.py'), '--submit'],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert r.returncode == 1
    assert ('not the exact approved' if s.APPROVED_SHA is None else 'explicit user approval') in r.stderr
    assert (p.OUTPUT/'results').exists() == had
