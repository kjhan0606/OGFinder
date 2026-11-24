import csv
import io
from pathlib import Path
import re
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_afterburner_photoz_rows as p
import submit_s16a_afterburner_photoz_rows as s

ID = '9007199254740993'


def payload(ids=(ID,), columns=None, overrides=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or ['# object_id', *list(p.COLUMNS)[1:]])
    for ident in ids:
        row = {k: (ident if k in ('object_id', *s.JOINED_IDS) else 'False' if k in s.BOOLEANS
                   else '3' if k in s.INTEGERS else '1.5e-30') for k in p.COLUMNS}
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


def test_query_joins_three_tables_without_filters():
    sql = p.query([ID, '9007199254740995'])
    assert len(p.COLUMNS) == 35 and sql.count('LEFT JOIN s16a_wide.') == 3
    assert 'b.yparent_flux_convolved_2_deconv AS' in sql and 'b.gparent_flux_convolved_seeing AS' in sql
    assert 'b.object_id AS afterburner_object_id' in sql and 'z.object_id AS photoz_object_id' in sql
    assert not re.search(r'\b(WHERE|LIMIT|DISTINCT|COALESCE)\b', sql)
    with pytest.raises(ValueError):
        p.query([ID, ID])


def test_null_row_absence_and_nan_are_counted_separately(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload(overrides={'photoz_object_id': '', 'photoz_best': '', 'photoz_risk_best': '',
                                        'yparent_flux_convolved_2_1': 'nan', 'yparent_flux_convolved_2_1_flags': 'True'}))
    result = s.validate_csv(file, [ID])
    assert result['null_by_column'] == {'photoz_object_id': 1, 'photoz_best': 1, 'photoz_risk_best': 1}
    assert result['nonfinite_by_column'] == {'yparent_flux_convolved_2_1': 1}
    assert result['true_by_column'] == {'yparent_flux_convolved_2_1_flags': 1}


@pytest.mark.parametrize('kind', ['missing', 'extra', 'duplicate', 'float', 'columns', 'boolean', 'joined', 'count', 'text'])
def test_refuses_scope_or_structural_corruption(tmp_path, kind):
    ids, columns, overrides = [ID], None, None
    if kind == 'missing': ids = []
    elif kind == 'extra': ids.append('9007199254740992')
    elif kind == 'duplicate': ids += ids
    elif kind == 'float': ids = ['9007199254740992.0']
    elif kind == 'columns': columns = ['object_id', *list(p.COLUMNS)[1:-1], 'unapproved']
    elif kind == 'boolean': overrides = {'gparent_flux_convolved_flag': 'yes'}
    elif kind == 'joined': overrides = {'afterburner_object_id': '9007199254740992'}
    elif kind == 'count': overrides = {'icountinputs': '2.5'}
    elif kind == 'text': overrides = {'iparent_flux_convolved_2_1': 'missing'}
    file = tmp_path/'a.csv'; file.write_bytes(payload(ids, columns, overrides))
    with pytest.raises(ValueError):
        s.validate_csv(file, [ID])


def test_resume_reuses_job_and_complete_data_never_calls_network(tmp_path):
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    client = Client(); checkpoint = tmp_path/'state.json'
    s.fetch_batch(client, {'name': 'b.sql'}, 'SELECT', [ID], tmp_path, st, checkpoint)
    assert client.calls == ['status', 'download']
    assert st['budget_bytes_charged'] == st['confirmed_bytes_received'] == len(payload())
    s.fetch_batch(client, {'name': 'b.sql'}, 'SELECT', [ID], tmp_path, st, checkpoint)
    assert client.calls == ['status', 'download']


def test_invalid_csv_stops_and_is_preserved(tmp_path):
    client = Client(payload(['7'])); st = state(); checkpoint = tmp_path/'state.json'
    for _ in range(2):
        with pytest.raises(ValueError):
            s.fetch_batch(client, {'name': 'b.sql'}, '', [ID], tmp_path, st, checkpoint)
    assert client.calls == ['submit', 'status', 'download']
    assert st['jobs']['b']['status'] == 'downloaded_invalid' and not (tmp_path/'b.csv').exists()


def test_malformed_csv_field_is_marked_invalid_and_not_downloaded_again(tmp_path):
    client = Client(payload()+b'"'+b'x'*140000+b'\n'); st = state(); checkpoint = tmp_path/'state.json'
    with pytest.raises(csv.Error):
        s.fetch_batch(client, {'name': 'b.sql'}, '', [ID], tmp_path, st, checkpoint)
    charged = st['budget_bytes_charged']
    with pytest.raises(ValueError, match='requires review'):
        s.fetch_batch(client, {'name': 'b.sql'}, '', [ID], tmp_path, st, checkpoint)
    assert st['jobs']['b']['status'] == 'downloaded_invalid' and st['budget_bytes_charged'] == charged
    assert client.calls == ['submit', 'status', 'download']


@pytest.mark.parametrize('secret', ['p"w\\x', 'pässwörd', 'a<b>&c', 'sl/ash'])
def test_redaction_covers_escaped_echoes(secret):
    class C:
        credential = {'account_name': 'someone', 'password': secret}
    import html, json
    body = ' '.join([secret, json.dumps(secret), json.dumps(secret, ensure_ascii=False), html.escape(secret), repr(secret),
                     secret.replace('/', '\\/')])
    assert secret not in s.redact(body, C) and json.dumps(secret)[1:-1] not in s.redact(body, C)


def test_sql_changed_after_verification_is_refused(tmp_path):
    (tmp_path/'b.sql').write_text('SELECT (9007199254740993::bigint)')
    batch = {'name': 'b.sql', 'sha256': s.sha(tmp_path/'b.sql')}
    assert s.batch_ids(tmp_path, batch) == [ID]
    (tmp_path/'b.sql').write_text('SELECT (9007199254740995::bigint)')
    with pytest.raises(ValueError, match='SQL changed'):
        s.batch_text(tmp_path, batch)


def test_environment_credential_requires_password(monkeypatch):
    class Args:
        use_env_auth, user = True, 'someone'
    monkeypatch.delenv('HSC_PASSWORD', raising=False)
    with pytest.raises(ValueError, match='not set'):
        s.credentials(Args)
    monkeypatch.setenv('HSC_PASSWORD', 'fixture'); monkeypatch.delenv('HSC_USER', raising=False)
    assert s.credentials(Args) == ('someone', 'fixture')


def test_prompt_credential_refuses_without_terminal():
    code = ('import sys; sys.path.insert(0, %r)\n'
            'from submit_s16a_afterburner_photoz_rows import credentials\n'
            'class A: use_env_auth=False; user="u"\n'
            'credentials(A)') % str(HERE)
    result = subprocess.run([sys.executable, '-c', code], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    assert result.returncode != 0 and 'interactive terminal' in result.stderr


def test_analysis_separates_absent_rows_null_and_nan():
    import audit_s16a_afterburner_photoz_rows as a
    row = {k: ('1' if k in ('object_id', *s.JOINED_IDS) else 'False' if k in s.BOOLEANS else '3' if k in s.INTEGERS
               else '1e-30') for k in p.COLUMNS}
    assert a.classify(row) == ('flux_finite_all_bands', 'photoz_row_finite')
    nan = dict(row, rparent_flux_convolved_2_1='nan', photoz_object_id='', photoz_best='', photoz_risk_best='')
    assert a.classify(nan) == ('flux_nonfinite_any_band', 'photoz_row_absent') and a.band_state(nan, 'r') == 'nonfinite'
    assert a.classify(dict(row, gparent_flux_convolved_2_1=''))[0] == 'flux_null_any_band'
    assert a.classify(dict(row, afterburner_object_id=''))[0] == 'afterburner_row_absent'
    assert a.classify(dict(row, photoz_best='nan'))[1] == 'photoz_row_value_unavailable'


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_reuses_the_executed_supplement_id_lists_exactly():
    assert p.validate(p.OUTPUT)['plan_sha256'] == s.APPROVED_SHA
    values = lambda text: text[text.index('FROM (VALUES'):text.index(') AS r(object_id)')]
    for i in range(24):
        new = (p.OUTPUT/f's16a_afterburner_photoz_rows_{i:03d}.sql').read_text()
        old = (p.SUPPLEMENT/f's16a_quality_supplement_{i:03d}.sql').read_text()
        assert values(new) == values(old)


@needs_plan
def test_submit_requires_explicit_plan_hash():
    checkpoint = p.OUTPUT/'results'/'status.json'
    before = s.sha(checkpoint) if checkpoint.exists() else None
    result = subprocess.run([sys.executable, str(HERE/'submit_s16a_afterburner_photoz_rows.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1 and 'explicit user approval' in result.stderr
    assert (s.sha(checkpoint) if checkpoint.exists() else None) == before
