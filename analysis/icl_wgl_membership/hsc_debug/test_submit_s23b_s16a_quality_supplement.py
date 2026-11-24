import csv
import io
import json
from pathlib import Path
import sys
import urllib.error
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import submit_s23b_s16a_quality_supplement as s


def payload(ids=('9007199254740993',), columns=None, overrides=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or ['# object_id', *list(s.EXPRESSIONS)[1:]])
    for ident in ids:
        row = {k: (ident if k == 'object_id' else ('f' if k in s.BAD_FLAGS else '1')) for k in s.EXPRESSIONS}
        row['is_primary'] = 't'; row.update(overrides or {})
        writer.writerow(row.values())
    return stream.getvalue().encode()


class Response(io.BytesIO):
    def __init__(self, data):
        super().__init__(data); self.headers = {'Content-Length': str(len(data))}


class Client:
    def __init__(self, data=None):
        self.calls = []; self.credential = {'password':'test-secret'}
        self.data = payload() if data is None else data
    def json(self, action, **fields):
        self.calls.append(action)
        return {'id': 42} if action == 'submit' else {'status':'done'}
    def post(self, action, **fields):
        self.calls.append(action); return Response(self.data)


def state():
    return dict(jobs={}, confirmed_bytes_received=0, budget_bytes_charged=0, source_rows_downloaded=False)


def test_exact_bigint_ids_preserve_null_cells_for_later_selection(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload(overrides={'photoz_best':'', 'photoz_risk_best':'nan', 'pixel_edge':'null'}))
    result = s.validate_csv(file, ['9007199254740993'])
    assert result['rows'] == 1
    assert result['missing_or_nonfinite_by_column'] == {'photoz_best':1, 'photoz_risk_best':1, 'pixel_edge':1}


@pytest.mark.parametrize('kind', ['missing','extra','duplicate','float','columns','boolean','short_row'])
def test_refuses_scope_or_structural_corruption(tmp_path, kind):
    ids = ['9007199254740993']; columns = None; overrides = None
    if kind == 'missing': ids = []
    elif kind == 'extra': ids.append('9007199254740992')
    elif kind == 'duplicate': ids += ids
    elif kind == 'float': ids = ['9007199254740992.0']
    elif kind == 'columns': columns = ['object_id', *list(s.EXPRESSIONS)[1:-1], 'unapproved']
    elif kind == 'boolean': overrides = {'pixel_edge':'yes'}
    file = tmp_path/'a.csv'; file.write_bytes(payload(ids, columns, overrides))
    if kind == 'short_row': file.write_bytes(file.read_bytes().rsplit(b',',1)[0]+b'\n')
    with pytest.raises(ValueError): s.validate_csv(file, ['9007199254740993'])


def test_resume_reuses_job_and_complete_data_never_calls_network(tmp_path):
    st = state(); st['jobs']['b'] = {'status':'submitted','job_id':'42'}
    client = Client(); checkpoint = tmp_path/'state.json'
    s.fetch_batch(client, {'name':'b.sql'}, 'SELECT...', ['9007199254740993'], tmp_path, st, checkpoint)
    assert client.calls == ['status','download']
    assert st['budget_bytes_charged'] == st['confirmed_bytes_received'] == len(payload())
    s.fetch_batch(client, {'name':'b.sql'}, 'SELECT...', ['9007199254740993'], tmp_path, st, checkpoint)
    assert client.calls == ['status','download']
    (tmp_path/'b.csv').write_text('changed')
    with pytest.raises(ValueError): s.fetch_batch(client, {'name':'b.sql'}, '', ['9007199254740993'], tmp_path, st, checkpoint)


def test_ambiguous_submit_is_not_repeated(tmp_path):
    st = state(); st['jobs']['b'] = {'status':'submitting'}; client = Client()
    with pytest.raises(ValueError, match='unknown'):
        s.fetch_batch(client, {'name':'b.sql'}, '', ['9007199254740993'], tmp_path, st, tmp_path/'state.json')
    assert client.calls == []


def test_invalid_csv_stops_and_is_preserved(tmp_path):
    client = Client(payload(['7'])); st = state(); checkpoint = tmp_path/'state.json'
    for _ in range(2):
        with pytest.raises(ValueError):
            s.fetch_batch(client, {'name':'b.sql'}, '', ['9007199254740993'], tmp_path, st, checkpoint)
    assert client.calls == ['submit','status','download']
    assert (tmp_path/'b.csv.part').exists() and not (tmp_path/'b.csv').exists()
    assert st['jobs']['b']['status'] == 'downloaded_invalid'


def test_schema_rejection_saves_redacted_diagnostic_and_no_job(tmp_path):
    class Rejected(Client):
        def json(self, action, **fields):
            raise urllib.error.HTTPError('https://example.invalid', 406, '', {}, io.BytesIO(b'bad SQL test-secret'))
    st = state(); checkpoint = tmp_path/'state.json'
    with pytest.raises(urllib.error.HTTPError):
        s.fetch_batch(Rejected(), {'name':'b.sql'}, '', ['9007199254740993'], tmp_path, st, checkpoint)
    assert st['jobs']['b']['http_code'] == 406 and 'job_id' not in st['jobs']['b']
    assert st['confirmed_bytes_received'] == st['budget_bytes_charged'] == 0
    assert 'test-secret' not in checkpoint.read_text()


def test_shared_download_cap_reserves_before_interrupted_read(tmp_path, monkeypatch):
    import submit_s23b_s16a_shape_qa_pilot as original
    monkeypatch.setattr(original, 'MAX_BYTES', 5)
    class Broken(Response):
        def read(self, amount): raise OSError('interrupted')
    st = state(); checkpoint = tmp_path/'state.json'
    with pytest.raises(OSError): s.budgeted_download(Broken(b'12345'), tmp_path/'part', st, checkpoint)
    assert json.loads(checkpoint.read_text())['budget_bytes_charged'] == 5
    assert st['confirmed_bytes_received'] == 0
