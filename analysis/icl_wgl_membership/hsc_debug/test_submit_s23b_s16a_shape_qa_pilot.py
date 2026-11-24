import csv
import copy
import io
import json
from pathlib import Path
import sys
import urllib.error

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import submit_s23b_s16a_shape_qa_pilot as submit


def center(ident='L'):
    return dict(center_id=ident, expected_raw_pairs=1, ra_deg=135., dec_deg=1., radius_arcsec=600.)


def payload(rows=None, columns=None):
    out = io.StringIO(newline='')
    writer = csv.writer(out)
    writer.writerow(columns or ['# '+submit.OUTPUT_COLUMNS[0], *submit.OUTPUT_COLUMNS[1:]])
    writer.writerows(rows or [['L','9007199254740993',135.01,1.01,.1,.2,.1,.3,10,-.02,.001,.002]])
    return out.getvalue().encode()


class Response(io.BytesIO):
    def __init__(self, data, length=True):
        super().__init__(data)
        self.headers = {'Content-Length':str(len(data))} if length else {}


def state():
    return dict(jobs={}, budget_bytes_charged=0, confirmed_bytes_received=0,
                source_rows_downloaded=False, credentials_saved=False)


class Client:
    def __init__(self, data=None):
        self.data = data or payload()
        self.calls = []
        self.credential = {'password':'test-secret'}
    def json(self, action, **fields):
        self.calls.append(action)
        return {'id':42} if action == 'submit' else {'status':'done'}
    def post(self, action, **fields):
        self.calls.append(action)
        return Response(self.data)


def test_csv_keeps_large_object_id_precise_and_accepts_nonfinite_calibration_for_qa(tmp_path):
    p = tmp_path/'data.csv'
    p.write_bytes(payload([['L','9007199254740993',135.01,1.01,.1,.2,'nan',.3,10,-.02,.001,.002]]))
    result = submit.validate_csv(p, [center()])
    assert result['rows'] == 1
    assert result['nonfinite_by_column'] == {'sigma_e':1}
    assert result['coordinates_pass_degree_cones'] is True


@pytest.mark.parametrize('kind', ['extra_column','unknown_center','duplicate','wrong_count','wrong_units','negative_id'])
def test_bad_csv_scope_or_coordinates_are_refused(tmp_path, kind):
    row = ['L','123',135.01,1.01,.1,.2,.1,.3,10,-.02,.001,.002]
    rows = [row]
    columns = None
    if kind == 'extra_column': columns = [*submit.OUTPUT_COLUMNS,'photoz'];rows=[row+['.8']]
    elif kind == 'unknown_center': row[0] = 'outside'
    elif kind == 'duplicate': rows = [row, row]
    elif kind == 'wrong_count': row = list(row);row[1]='124';rows.append(row)
    elif kind == 'wrong_units': row[2] = 2.35
    elif kind == 'negative_id': row[1] = '-1'
    p = tmp_path/'data.csv';p.write_bytes(payload(rows,columns))
    with pytest.raises(ValueError): submit.validate_csv(p, [center()])


def test_same_object_in_distinct_centers_is_allowed(tmp_path):
    a = ['L','123',135.01,1.01,.1,.2,.1,.3,10,-.02,.001,.002]
    b = list(a);b[0]='R'
    p=tmp_path/'data.csv';p.write_bytes(payload([a,b]))
    assert submit.validate_csv(p,[center('L'),center('R')])['rows'] == 2


def test_download_budget_applies_to_all_attempts_and_refuses_oversized_response(tmp_path, monkeypatch):
    monkeypatch.setattr(submit,'MAX_BYTES',10)
    s=state();checkpoint=tmp_path/'state.json'
    submit.budgeted_download(Response(b'123456'),tmp_path/'first',s,checkpoint)
    assert s['confirmed_bytes_received'] == s['budget_bytes_charged'] == 6
    with pytest.raises(ValueError,match='budget'):
        submit.budgeted_download(Response(b'56789'),tmp_path/'second',s,checkpoint)
    assert not (tmp_path/'second').exists()
    assert s['source_rows_downloaded'] is True


def test_interrupted_read_reserves_budget_before_io(tmp_path, monkeypatch):
    monkeypatch.setattr(submit,'MAX_BYTES',10)
    class Failed(Response):
        def read(self,n): raise OSError('interrupted')
    s=state();checkpoint=tmp_path/'state.json'
    with pytest.raises(OSError):
        submit.budgeted_download(Failed(b'123456'),tmp_path/'partial',s,checkpoint)
    recorded=json.loads(checkpoint.read_text())
    assert recorded['budget_bytes_charged'] == 6
    assert recorded['confirmed_bytes_received'] == 0
    assert (tmp_path/'partial').exists()


def test_unknown_length_is_stopped_at_byte_cap_without_reading_extra(tmp_path, monkeypatch):
    monkeypatch.setattr(submit,'MAX_BYTES',5)
    response=Response(b'123456789',length=False);s=state()
    with pytest.raises(ValueError,match='budget'):
        submit.budgeted_download(response,tmp_path/'partial',s,tmp_path/'state.json')
    assert response.tell() == 5
    assert (tmp_path/'partial').read_bytes() == b'12345'
    assert s['budget_bytes_charged'] == 5


def test_resume_uses_existing_job_and_completed_files_need_no_network(tmp_path):
    s=state();s['jobs']['b']={'status':'submitted','job_id':'42'}
    batch={'name':'b.sql'}; client=Client();checkpoint=tmp_path/'state.json'
    submit.fetch_batch(client,batch,'SELECT ...',[center()],tmp_path,s,checkpoint)
    assert client.calls == ['status','download']
    assert s['jobs']['b']['status'] == 'complete'
    assert (tmp_path/'b.csv').is_file()
    submit.fetch_batch(client,batch,'SELECT ...',[center()],tmp_path,s,checkpoint)
    assert client.calls == ['status','download']
    (tmp_path/'b.csv').write_bytes(b'changed')
    with pytest.raises(ValueError,match='changed'):
        submit.fetch_batch(client,batch,'SELECT ...',[center()],tmp_path,s,checkpoint)


def test_ambiguous_submit_is_never_repeated(tmp_path):
    s=state();s['jobs']['b']={'status':'submitting'};client=Client()
    with pytest.raises(ValueError,match='unknown'):
        submit.fetch_batch(client,{'name':'b.sql'},'SELECT ...',[center()],tmp_path,s,tmp_path/'state.json')
    assert client.calls == []


def test_invalid_download_is_preserved_and_stops_resume(tmp_path):
    client=Client(payload(columns=[*submit.OUTPUT_COLUMNS[:-1],'photoz']))
    s=state();batch={'name':'b.sql'};checkpoint=tmp_path/'state.json'
    with pytest.raises(ValueError,match='columns'):
        submit.fetch_batch(client,batch,'SELECT ...',[center()],tmp_path,s,checkpoint)
    assert (tmp_path/'b.csv.part').is_file() and not (tmp_path/'b.csv').exists()
    assert s['jobs']['b']['status'] == 'downloaded_invalid'
    calls=list(client.calls)
    with pytest.raises(ValueError,match='review'):
        submit.fetch_batch(client,batch,'SELECT ...',[center()],tmp_path,s,checkpoint)
    assert client.calls == calls


def test_http_rejection_keeps_redacted_diagnostic_without_job_or_retry(tmp_path):
    class Rejected(Client):
        def json(self, action, **fields):
            self.calls.append(action)
            raise urllib.error.HTTPError('https://example.invalid',406,'schema',{},
                io.BytesIO(b'column absent; test-secret'))
    client=Rejected();s=state();checkpoint=tmp_path/'state.json'
    with pytest.raises(urllib.error.HTTPError):
        submit.fetch_batch(client,{'name':'b.sql'},'SELECT ...',[center()],tmp_path,s,checkpoint)
    saved=json.loads(checkpoint.read_text())
    assert client.calls == ['submit']
    assert 'job_id' not in saved['jobs']['b']
    assert saved['jobs']['b']['http_code'] == 406
    assert 'test-secret' not in checkpoint.read_text()
    assert 'column absent' in saved['jobs']['b']['redacted_schema_error']


@pytest.mark.parametrize('changed', [None,'centers','maximum_total_received_bytes','other_column','batch'])
def test_identifier_amendment_cannot_expand_approved_semantic_scope(changed):
    original=dict(centers=[{'center_id':'L'}], maximum_total_received_bytes=128,
        column_expressions={'weight':submit.OLD_WEIGHT,'e1':'a.e1'},
        batches=[{'name':'b.sql','center_ids':['L'],'expected_raw_pairs':1}])
    amended=copy.deepcopy(original);amended['column_expressions']['weight']=submit.NEW_WEIGHT
    if changed == 'centers': amended['centers'].append({'center_id':'extra'})
    elif changed == 'maximum_total_received_bytes': amended['maximum_total_received_bytes']=129
    elif changed == 'other_column': amended['column_expressions']['e1']='a.other'
    elif changed == 'batch': amended['batches'][0]['center_ids']=['extra']
    if changed is None: submit.same_scope(original,amended)
    else:
        with pytest.raises(ValueError): submit.same_scope(original,amended)
