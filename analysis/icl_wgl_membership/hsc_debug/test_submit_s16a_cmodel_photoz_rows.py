import csv
import io
from pathlib import Path
import re
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_cmodel_photoz_rows as p
import submit_s16a_cmodel_photoz_rows as s

ID = '9007199254740993'


def payload(ids=(ID,), columns=None, overrides=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or ['# object_id', *list(p.COLUMNS)[1:]])
    for ident in ids:
        row = {k: (ident if k in ('object_id', *s.JOINED_IDS) else '0.85') for k in p.COLUMNS}
        row.update(overrides or {})
        writer.writerow(row.values())
    return stream.getvalue().encode()


class Response(io.BytesIO):
    def __init__(self, data):
        super().__init__(data); self.headers = {'Content-Length': str(len(data))}


class Client:
    def __init__(self, data=None):
        self.calls = []; self.credential = {'password': 'test-secret'}
        self.data = payload() if data is None else data
    def json(self, action, **fields):
        self.calls.append(action)
        return {'id': 42} if action == 'submit' else {'status': 'done'}
    def post(self, action, **fields):
        self.calls.append(action); return Response(self.data)


def state():
    return dict(jobs={}, confirmed_bytes_received=0, budget_bytes_charged=0, source_rows_downloaded=False)


def test_query_joins_two_photoz_tables_without_filters():
    sql = p.query([ID, '9007199254740995'])
    assert len(p.COLUMNS) == 13 and sql.count('LEFT JOIN s16a_wide.') == 2
    assert 'm.reduced_chisq AS mizuki_reduced_chisq' in sql and 'm.prob_star AS mizuki_prob_star' in sql
    assert 'e.object_id AS ephor_object_id' in sql and 'm.object_id AS mizuki_object_id' in sql
    assert 'photoz_ephor_ab' not in sql
    assert not re.search(r'\b(WHERE|LIMIT|DISTINCT|COALESCE)\b', sql)


def test_null_rows_and_nan_are_counted_separately(tmp_path):
    file = tmp_path/'a.csv'
    mizuki = [c for c in p.COLUMNS if c.startswith('mizuki_')]
    file.write_bytes(payload(overrides={**{c: '' for c in mizuki}, 'ephor_photoz_std_best': 'nan'}))
    result = s.validate_csv(file, [ID])
    assert result['null_by_column'] == {c: 1 for c in mizuki}
    assert result['nonfinite_by_column'] == {'ephor_photoz_std_best': 1}


@pytest.mark.parametrize('kind', ['missing', 'extra', 'duplicate', 'float', 'columns', 'joined', 'text'])
def test_refuses_scope_or_structural_corruption(tmp_path, kind):
    ids, columns, overrides = [ID], None, None
    if kind == 'missing': ids = []
    elif kind == 'extra': ids.append('9007199254740992')
    elif kind == 'duplicate': ids += ids
    elif kind == 'float': ids = ['9007199254740992.0']
    elif kind == 'columns': columns = ['object_id', *list(p.COLUMNS)[1:-1], 'unapproved']
    elif kind == 'joined': overrides = {'ephor_object_id': '9007199254740992'}
    elif kind == 'text': overrides = {'ephor_photoz_best': 'missing'}
    file = tmp_path/'a.csv'; file.write_bytes(payload(ids, columns, overrides))
    with pytest.raises(ValueError):
        s.validate_csv(file, [ID])


def test_resume_and_malformed_payload_gates(tmp_path):
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    client = Client(); checkpoint = tmp_path/'state.json'
    s.fetch_batch(client, {'name': 'b.sql'}, 'SELECT', [ID], tmp_path, st, checkpoint)
    s.fetch_batch(client, {'name': 'b.sql'}, 'SELECT', [ID], tmp_path, st, checkpoint)
    assert client.calls == ['status', 'download']
    bad = Client(payload()+b'"'+b'x'*140000+b'\n'); st2 = state(); cp2 = tmp_path/'state2.json'
    with pytest.raises(csv.Error):
        s.fetch_batch(bad, {'name': 'c.sql'}, '', [ID], tmp_path, st2, cp2)
    with pytest.raises(ValueError, match='requires review'):
        s.fetch_batch(bad, {'name': 'c.sql'}, '', [ID], tmp_path, st2, cp2)
    assert st2['jobs']['c']['status'] == 'downloaded_invalid' and bad.calls == ['submit', 'status', 'download']


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_reuses_the_executed_supplement_id_lists_exactly():
    assert p.validate(p.OUTPUT)['plan_sha256'] == s.APPROVED_SHA
    values = lambda text: text[text.index('FROM (VALUES'):text.index(') AS r(object_id)')]
    for i in range(24):
        new = (p.OUTPUT/f's16a_cmodel_photoz_rows_{i:03d}.sql').read_text()
        old = (p.SUPPLEMENT/f's16a_quality_supplement_{i:03d}.sql').read_text()
        assert values(new) == values(old)


@needs_plan
def test_submit_requires_explicit_plan_hash():
    # Stops at the approval-hash gate, before any results directory, credential or network use.
    checkpoint = p.OUTPUT/'results'/'status.json'
    had_results = (p.OUTPUT/'results').exists()
    before = s.sha(checkpoint) if checkpoint.exists() else None
    result = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_cmodel_photoz_rows.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1 and 'explicit user approval' in result.stderr
    assert (s.sha(checkpoint) if checkpoint.exists() else None) == before
    assert (p.OUTPUT/'results').exists() == had_results
