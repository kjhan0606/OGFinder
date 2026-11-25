import csv
import io
from pathlib import Path
import re
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_fdfc_grid as p
import submit_s16a_fdfc_grid as s

LENS = dict(lens_id='S23B_9999', center_id='LC-S23B_9999', ra_deg='150.0', dec_deg='2.0', z_lens='0.4', radius_arcsec='550.0')
BATCH = dict(name='b.sql', lens=LENS, row_limit=5)


def payload(cells=((0, 0, 30, 29, 0),), lens='S23B_9999', columns=None, n_cells=None, n_objects=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or p.OUTPUT_COLUMNS)
    total = sum(c[2] for c in cells)
    for x, y, n, f, null in cells:
        row = dict(lens_id=lens, cell_x=x, cell_y=y, n_primary=n, n_fdfc=f, n_fdfc_null=null,
                   n_cells=len(cells) if n_cells is None else n_cells, n_objects=total if n_objects is None else n_objects)
        writer.writerow([row.get(k, '1') for k in (columns or p.OUTPUT_COLUMNS)])
    return stream.getvalue().encode()


def test_query_is_aggregated_primary_cone_only():
    sql = p.query(LENS)
    where = sql[sql.index('WHERE'):sql.index('GROUP BY')]
    assert 'coneSearch(f.coord, 150.0000000000, 2.0000000000, 5400.0)' in where and where.count('\n  AND ') == 1
    assert 'f.detect_is_primary' in where and 'FROM s16a_wide.forced AS f' in sql
    assert 'count(*) OVER () AS n_cells' in sql and 'sum(count(*)) OVER () AS n_objects' in sql
    assert not re.search(r'object_id|JOIN|cmodel|shape', sql) and sql.rstrip().endswith(f'LIMIT {p.LIMIT_PER_LENS};')
    assert len(p.OUTPUT_COLUMNS) == 8


def test_valid_batch(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload([(0, 0, 30, 29, 0), (-90, 5, 3, 0, 1)]))
    assert s.validate_csv(file, BATCH) == dict(rows=2, n_primary=33, n_fdfc=29, n_fdfc_null=1)


@pytest.mark.parametrize('kind', ['lens', 'duplicate', 'outside', 'counts', 'columns', 'limit', 'cut_cells', 'cut_objects', 'float'])
def test_refuses_corruption_and_truncation(tmp_path, kind):
    cells, lens, columns, nc, no = [(0, 0, 30, 29, 0)], 'S23B_9999', None, None, None
    if kind == 'lens': lens = 'S23B_1'
    elif kind == 'duplicate': cells = [(0, 0, 30, 29, 0), (0, 0, 5, 1, 0)]
    elif kind == 'outside': cells = [(95, 0, 30, 29, 0)]
    elif kind == 'counts': cells = [(0, 0, 30, 29, 2)]
    elif kind == 'columns': columns = [*p.OUTPUT_COLUMNS[:-1], 'unapproved']
    elif kind == 'limit': cells = [(i, 0, 1, 0, 0) for i in range(5)]
    elif kind == 'cut_cells': nc = 2
    elif kind == 'cut_objects': no = 31
    file = tmp_path/'a.csv'
    data = payload(cells, lens, columns, nc, no)
    if kind == 'float': data = data.replace(b',30,29,', b',30.5,29,')
    file.write_bytes(data)
    with pytest.raises(ValueError):
        s.validate_csv(file, BATCH)


class Response(io.BytesIO):
    def __init__(self, data, headers=None):
        super().__init__(data); self.headers = {'Content-Length': str(len(data))} if headers is None else headers


class Client:
    def __init__(self, data=None, headers=None):
        self.calls = []; self.credential = {'password': 'test-secret'}
        self.data = payload() if data is None else data
        self.headers = headers
    def json(self, action, **fields):
        self.calls.append(action)
        return {'id': 42} if action == 'submit' else {'status': 'done'}
    def post(self, action, **fields):
        self.calls.append(action); return Response(self.data, self.headers)


def state():
    return dict(jobs={}, confirmed_bytes_received=0, budget_bytes_charged=0, source_rows_downloaded=False)


def test_resume_transport_and_line_end_gates(tmp_path):
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    client = Client(); checkpoint = tmp_path/'state.json'
    s.fetch_batch(client, BATCH, 'SELECT', tmp_path, st, checkpoint)
    s.fetch_batch(client, BATCH, 'SELECT', tmp_path, st, checkpoint)
    assert client.calls == ['status', 'download'] and st['jobs']['b']['rows'] == 1
    st = state(); st['jobs']['c'] = {'status': 'submitted', 'job_id': '42'}
    with pytest.raises(ValueError, match='neither Content-Length nor chunked'):
        s.fetch_batch(Client(headers={}), dict(BATCH, name='c.sql'), 'SELECT', tmp_path, st, tmp_path/'s2.json')
    st = state(); st['jobs']['d'] = {'status': 'submitted', 'job_id': '42'}
    with pytest.raises(ValueError, match='line break'):
        s.fetch_batch(Client(payload()[:-1], headers={'Transfer-Encoding': 'chunked'}), dict(BATCH, name='d.sql'), 'SELECT', tmp_path, st, tmp_path/'s3.json')
    assert st['jobs']['d']['status'] == 'downloaded_invalid'


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_verifies_and_covers_the_eleven_lenses():
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['batches'] == 11


@needs_plan
def test_submit_refuses_without_approval():
    checkpoint = p.OUTPUT/'results'/'status.json'
    had_results = (p.OUTPUT/'results').exists()
    result = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_fdfc_grid.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1
    assert ('not the exact approved FDFC grid plan' if s.APPROVED_SHA is None else 'explicit user approval') in result.stderr
    assert (p.OUTPUT/'results').exists() == had_results and not checkpoint.exists() or had_results
