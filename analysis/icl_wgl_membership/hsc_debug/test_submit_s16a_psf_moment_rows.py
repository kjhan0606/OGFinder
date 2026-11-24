import csv
import io
from pathlib import Path
import re
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_psf_moment_rows as p
import submit_s16a_psf_moment_rows as s

ID = '9007199254740993'


def payload(ids=(ID,), columns=None, overrides=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or ['# object_id', *list(p.COLUMNS)[1:]])
    for ident in ids:
        row = {k: (ident if k in ('object_id', *s.JOINED_IDS) else '0.12') for k in p.COLUMNS}
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


def test_query_joins_forced_psf_moments_only_without_filters():
    sql = p.query([ID, '9007199254740995'])
    assert len(p.COLUMNS) == 17 and sql.count('LEFT JOIN s16a_wide.') == 1 and 'LEFT JOIN s16a_wide.forced AS f' in sql
    assert all(f'f.{b}shape_sdss_psf_{m} AS {b}shape_sdss_psf_{m}' in sql for b in 'grizy' for m in ('11', '22', '12'))
    assert 'f.object_id AS forced_object_id' in sql
    assert not re.search(r'\b(WHERE|LIMIT|DISTINCT|COALESCE)\b', sql)


def test_null_rows_and_nan_are_counted_separately(tmp_path):
    file = tmp_path/'a.csv'
    moments = [c for c in p.COLUMNS if c.startswith('y')]
    file.write_bytes(payload(overrides={**{c: '' for c in moments}, 'rshape_sdss_psf_12': 'nan'}))
    result = s.validate_csv(file, [ID])
    assert result['null_by_column'] == {c: 1 for c in moments}
    assert result['nonfinite_by_column'] == {'rshape_sdss_psf_12': 1}


def test_absent_forced_row_is_counted_but_moments_without_it_are_refused(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload(overrides={c: '' for c in list(p.COLUMNS)[1:]}))
    assert s.validate_csv(file, [ID])['null_by_column'] == {c: 1 for c in list(p.COLUMNS)[1:]}
    file.write_bytes(payload(overrides={'forced_object_id': ''}))
    with pytest.raises(ValueError, match='without a forced row'):
        s.validate_csv(file, [ID])


@pytest.mark.parametrize('kind', ['missing', 'extra', 'duplicate', 'float', 'columns', 'joined', 'text'])
def test_refuses_scope_or_structural_corruption(tmp_path, kind):
    ids, columns, overrides = [ID], None, None
    if kind == 'missing': ids = []
    elif kind == 'extra': ids.append('9007199254740992')
    elif kind == 'duplicate': ids += ids
    elif kind == 'float': ids = ['9007199254740992.0']
    elif kind == 'columns': columns = ['object_id', *list(p.COLUMNS)[1:-1], 'unapproved']
    elif kind == 'joined': overrides = {'forced_object_id': '9007199254740992'}
    elif kind == 'text': overrides = {'rshape_sdss_psf_11': 'missing'}
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
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256'])
    values = lambda text: text[text.index('FROM (VALUES'):text.index(') AS r(object_id)')]
    for i in range(24):
        new = (p.OUTPUT/f's16a_psf_moment_rows_{i:03d}.sql').read_text()
        old = (p.SUPPLEMENT/f's16a_quality_supplement_{i:03d}.sql').read_text()
        assert values(new) == values(old)


@needs_plan
def test_submit_refuses_without_approval():
    # Stops before any results directory, credential or network use: unapproved plan, or no explicit hash.
    checkpoint = p.OUTPUT/'results'/'status.json'
    had_results = (p.OUTPUT/'results').exists()
    before = s.sha(checkpoint) if checkpoint.exists() else None
    result = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_psf_moment_rows.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1
    assert ('not the exact approved PSF-moment plan' if s.APPROVED_SHA is None else 'explicit user approval') in result.stderr
    assert (s.sha(checkpoint) if checkpoint.exists() else None) == before
    assert (p.OUTPUT/'results').exists() == had_results
