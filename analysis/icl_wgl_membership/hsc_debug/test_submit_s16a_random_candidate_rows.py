import csv
import io
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_random_candidate_rows as p
import submit_s16a_random_candidate_rows as s

LENS = dict(lens_id='S23B_9999', center_id='LC-S23B_9999', ra_deg='150.0', dec_deg='2.0', z_lens='0.4', radius_arcsec='550.0')
BATCH = dict(name='b.sql', lens=LENS, row_limit=5)


def offset(ra, dec, arcmin, angle=0.3):
    """A point at the given angular distance from (ra, dec), built on the unit sphere."""
    c = np.array(p.unit_vector(ra, dec))
    east = np.array([-math.sin(math.radians(ra)), math.cos(math.radians(ra)), 0.0])
    north = np.cross(c, east)
    t = math.radians(arcmin/60)
    v = math.cos(t)*c+math.sin(t)*(math.cos(angle)*east+math.sin(angle)*north)
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2]))


def payload(points=((1, 10.0),), lens='S23B_9999', columns=None, overrides=None, n_matched=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or p.OUTPUT_COLUMNS)
    for ident, arcmin in points:
        ra, dec = offset(150.0, 2.0, arcmin)
        row = {k: ('False' if k in p.BOOLEANS else '7' if k in p.INTEGERS else '0.5') for k in p.OUTPUT_COLUMNS}
        row.update(lens_id=lens, object_id=str(ident), ra=repr(ra), dec=repr(dec), iadjust_density='0.001', skymap_id='93750107',
                   n_matched=str(len(points) if n_matched is None else n_matched))
        row.update(overrides or {})
        writer.writerow([row.get(k, '0.5') for k in (columns or p.OUTPUT_COLUMNS)])
    return stream.getvalue().encode()


def test_query_has_only_box_cone_and_thinning_filters():
    sql = p.query(LENS)
    where = sql[sql.index('WHERE'):sql.index('ORDER BY')]
    assert len(p.OUTPUT_COLUMNS) == len(p.COLUMNS)+2 == 47 and "'S23B_9999'::text AS lens_id" in sql and 'count(*) OVER () AS n_matched' in sql and 'FROM s16a_wide.random AS r' in sql
    assert where.count('\n  AND ') == 3 and where.count('BETWEEN') == 2 and 'r.dec BETWEEN' in where and 'r.ra BETWEEN' in where and 'r.cx*(' in where
    assert f'r.iadjust_density < {p.ADJUST_MAX}' in where and not re.search(r'flags|detect|countinputs|JOIN|DISTINCT', where+sql[sql.index('FROM'):])
    assert sql.rstrip().endswith(f'ORDER BY r.object_id\nLIMIT {p.LIMIT_PER_LENS};')


@pytest.mark.parametrize('dec', [-1.5, 2.0, 42.8])
def test_box_contains_the_cone(dec):
    (dlo, dhi), (alo, ahi) = p.box(150.0, dec)
    for angle in np.linspace(0, 2*np.pi, 721):
        ra, d = offset(150.0, dec, p.RADIUS_ARCMIN, angle)
        assert alo <= ra <= ahi and dlo <= d <= dhi


def test_box_refuses_wrap_and_pole():
    with pytest.raises(ValueError):
        p.box(0.3, 0.0)
    with pytest.raises(ValueError):
        p.box(150.0, 88.5)


def test_valid_batch_counts_flags_and_nulls(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload([(1, 10.0), (2, p.RADIUS_ARCMIN-0.1)], overrides={'iflags_pixel_bright_object_any': 't', 'isky_std': ''}))
    out = s.validate_csv(file, BATCH)
    assert out['rows'] == 2 and out['true_by_column'] == {'iflags_pixel_bright_object_any': 2} and out['null_by_column'] == {'isky_std': 2}


@pytest.mark.parametrize('kind', ['lens', 'duplicate', 'outside', 'density', 'boolean', 'columns', 'limit', 'null_ra', 'empty', 'integer',
                                  'short', 'count_text'])
def test_refuses_corruption_and_truncation(tmp_path, kind):
    points, lens, columns, overrides, n_matched = [(1, 10.0)], 'S23B_9999', None, None, None
    if kind == 'lens': lens = 'S23B_0001'
    elif kind == 'duplicate': points = [(1, 10.0), (1, 20.0)]
    elif kind == 'outside': points = [(1, p.RADIUS_ARCMIN+0.05)]
    elif kind == 'density': overrides = {'iadjust_density': repr(p.ADJUST_MAX)}
    elif kind == 'boolean': overrides = {'gflags_pixel_edge': 'maybe'}
    elif kind == 'columns': columns = [*p.OUTPUT_COLUMNS[:-1], 'unapproved']
    elif kind == 'limit': points = [(i, 10.0) for i in range(1, 6)]
    elif kind == 'null_ra': overrides = {'ra': ''}
    elif kind == 'empty': points = []
    elif kind == 'integer': overrides = {'icountinputs': '4.5'}
    elif kind == 'short': points, n_matched = [(1, 10.0), (2, 20.0)], 3
    elif kind == 'count_text': overrides = {'n_matched': 'x'}
    file = tmp_path/'a.csv'; file.write_bytes(payload(points, lens, columns, overrides, n_matched))
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


def test_resume_and_malformed_payload_gates(tmp_path):
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    client = Client(); checkpoint = tmp_path/'state.json'
    s.fetch_batch(client, BATCH, 'SELECT', tmp_path, st, checkpoint)
    s.fetch_batch(client, BATCH, 'SELECT', tmp_path, st, checkpoint)
    assert client.calls == ['status', 'download'] and st['jobs']['b']['rows'] == 1
    bad = Client(payload(overrides={'iadjust_density': '0.5'})); st2 = state(); cp2 = tmp_path/'state2.json'
    batch = dict(BATCH, name='c.sql')
    with pytest.raises(ValueError, match='thinning'):
        s.fetch_batch(bad, batch, '', tmp_path, st2, cp2)
    with pytest.raises(ValueError, match='requires review'):
        s.fetch_batch(bad, batch, '', tmp_path, st2, cp2)
    assert st2['jobs']['c']['status'] == 'downloaded_invalid' and bad.calls == ['submit', 'status', 'download']


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_verifies_and_covers_the_eleven_lenses():
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['batches'] == 11
    lenses = [p.query(lens) for _, lens in p.ledger()]
    assert len(set(lenses)) == 11


@needs_plan
def test_submit_refuses_without_approval():
    checkpoint = p.OUTPUT/'results'/'status.json'
    had_results = (p.OUTPUT/'results').exists()
    before = s.sha(checkpoint) if checkpoint.exists() else None
    result = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_random_candidate_rows.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1
    assert ('not the exact approved random-candidate plan' if s.APPROVED_SHA is None else 'explicit user approval') in result.stderr
    assert (s.sha(checkpoint) if checkpoint.exists() else None) == before
    assert (p.OUTPUT/'results').exists() == had_results


def test_download_without_length_or_chunking_and_cut_line_are_refused(tmp_path):
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    with pytest.raises(ValueError, match='neither Content-Length nor chunked'):
        s.fetch_batch(Client(headers={}), BATCH, 'SELECT', tmp_path, st, tmp_path/'s.json')
    st = state(); st['jobs']['b'] = {'status': 'submitted', 'job_id': '42'}
    body = payload([(1, 10.0), (2, 20.0)])[:-2]
    with pytest.raises(ValueError, match='line break'):
        s.fetch_batch(Client(body, headers={'Transfer-Encoding': 'chunked'}), BATCH, 'SELECT', tmp_path, st, tmp_path/'s2.json')
    assert st['jobs']['b']['status'] == 'downloaded_invalid' and st['jobs']['b']['transport'] == 'chunked'


def test_poll_budget_and_extended_redaction():
    import cas_transfer_checks as c
    with pytest.raises(ValueError, match='polling budget'):
        c.require_poll_budget(c.MAX_POLL_SECONDS)
    client = Client(); client.credential = {'password': "pw'<x/ y"}
    text = "a pw&#39;&lt;x/ y b pw%27%3Cx%2F%20y c pw'<x/ y"
    assert 'pw' not in c.redact(text, client)


def test_transport_uses_http_client_attributes():
    import cas_transfer_checks as c
    class Real:
        headers = {'Transfer-Encoding': 'gzip, chunked'}
        def __init__(self, chunked, length): self.chunked, self.length = chunked, length
    assert c.transport(Real(True, None)) == 'chunked' and c.transport(Real(False, 10)) == 'content-length'
    with pytest.raises(ValueError, match='neither'):
        c.transport(Real(False, None))
