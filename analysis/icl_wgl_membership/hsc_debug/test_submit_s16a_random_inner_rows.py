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
import prepare_s16a_random_inner_rows as p
import submit_s16a_random_inner_rows as s
from prepare_s16a_random_candidate_rows import unit_vector

CENTRES = [dict(center_id='LC-S23B_9999', role='lens', ra_deg='150.0', dec_deg='2.0', z_lens='0.4', radius_arcsec='550.0',
                query_radius_arcsec=550.0*0.4/3+30),
           dict(center_id='RC-abc', role='random', ra_deg='150.5', dec_deg='2.2', z_lens='0.4', radius_arcsec='550.0',
                query_radius_arcsec=550.0*0.4/3+30)]
BATCH = dict(name='b.sql', lens_id='S23B_9999', centres=CENTRES, row_limit_per_centre=4)


def offset(ra, dec, arcsec, angle=0.3):
    c = np.array(unit_vector(ra, dec))
    east = np.array([-math.sin(math.radians(ra)), math.cos(math.radians(ra)), 0.0])
    north = np.cross(c, east)
    t = math.radians(arcsec/3600)
    v = math.cos(t)*c+math.sin(t)*(math.cos(angle)*east+math.sin(angle)*north)
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(v[2]))


def payload(points=(('LC-S23B_9999', 1, 10.0),), columns=None, overrides=None, counts=None, total=None):
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(columns or p.OUTPUT_COLUMNS)
    centres = {c['center_id']: c for c in CENTRES}
    for cid, ident, arcsec in points:
        c = centres.get(cid, CENTRES[0])
        ra, dec = offset(float(c['ra_deg']), float(c['dec_deg']), arcsec)
        row = {k: ('False' if k in p.BOOLEANS else '7' if k in p.INTEGERS else '0.5') for k in p.OUTPUT_COLUMNS}
        per = sum(1 for q in points if q[0] == cid)
        row.update(center_id=cid, object_id=str(ident), ra=repr(ra), dec=repr(dec), skymap_id='93750107',
                   n_matched=str((counts or {}).get(cid, per)), n_batch=str(len(points) if total is None else total))
        row.update(overrides or {})
        writer.writerow([row.get(k, '0.5') for k in (columns or p.OUTPUT_COLUMNS)])
    return stream.getvalue().encode()


def test_query_is_six_unfiltered_unthinned_cones():
    group = [dict(c, center_id=f'RC-{i}') for i, c in enumerate(CENTRES*3)]
    sql = p.query(group)
    assert sql.count('UNION ALL') == 5 and sql.count('FROM s16a_wide.random AS r') == 6 and sql.rstrip().endswith(') AS c)\n) AS u;')
    assert sql.startswith('SELECT u.*, count(*) OVER () AS n_batch FROM (') and sql.count('count(*) OVER () AS n_matched') == 6
    assert sql.count(f'LIMIT {p.LIMIT_PER_CENTRE}') == 6 and sql.count("::text AS center_id") == 6
    assert not re.search(r'adjust_density|JOIN|DISTINCT', sql)
    for block in sql.split('UNION ALL'):
        where = block[block.index('WHERE'):block.index('ORDER BY')]
        assert where.count('\n    AND ') == 2 and 'flags' not in where and 'detect' not in where
    assert len(p.OUTPUT_COLUMNS) == 35


def test_radius_is_inner_disc_plus_margin():
    assert abs(p.radius_arcsec({'radius_arcsec': '600'})-(80+30)) < 1e-12


def test_valid_batch_allows_an_empty_random_centre(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload([('LC-S23B_9999', 1, 5.0), ('LC-S23B_9999', 2, 90.0)], overrides={'iflags_pixel_bright_object_center': 't'}))
    out = s.validate_csv(file, BATCH)
    assert out['rows'] == 2 and out['rows_by_centre'] == {'LC-S23B_9999': 2, 'RC-abc': 0}
    assert out['true_by_column'] == {'iflags_pixel_bright_object_center': 2}


@pytest.mark.parametrize('kind', ['centre', 'duplicate', 'outside', 'boolean', 'columns', 'limit', 'no_lens', 'integer', 'short_centre',
                                  'missing_centre'])
def test_refuses_corruption_and_truncation(tmp_path, kind):
    points, columns, overrides, counts, total = [('LC-S23B_9999', 1, 10.0)], None, None, None, None
    if kind == 'centre': points = [('RC-other', 1, 10.0)]
    elif kind == 'duplicate': points = [('LC-S23B_9999', 1, 10.0), ('LC-S23B_9999', 1, 20.0)]
    elif kind == 'outside': points = [('LC-S23B_9999', 1, CENTRES[0]['query_radius_arcsec']+0.5)]
    elif kind == 'boolean': overrides = {'gflags_pixel_edge': 'maybe'}
    elif kind == 'columns': columns = [*p.OUTPUT_COLUMNS[:-1], 'unapproved']
    elif kind == 'limit': points = [('LC-S23B_9999', i, 10.0) for i in range(1, 5)]
    elif kind == 'no_lens': points = [('RC-abc', 1, 10.0)]
    elif kind == 'integer': overrides = {'icountinputs': '4.5'}
    elif kind == 'short_centre': counts = {'LC-S23B_9999': 2}
    elif kind == 'missing_centre': total = 2
    file = tmp_path/'a.csv'; file.write_bytes(payload(points, columns, overrides, counts, total))
    with pytest.raises(ValueError):
        s.validate_csv(file, BATCH)


def test_same_point_may_serve_two_centres(tmp_path):
    file = tmp_path/'a.csv'
    file.write_bytes(payload([('LC-S23B_9999', 7, 5.0), ('RC-abc', 7, 5.0)]))
    assert s.validate_csv(file, BATCH)['rows'] == 2


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
    assert client.calls == ['status', 'download'] and st['jobs']['b']['transport'] == 'content-length'
    st = state(); st['jobs']['c'] = {'status': 'submitted', 'job_id': '42'}
    with pytest.raises(ValueError, match='neither Content-Length nor chunked'):
        s.fetch_batch(Client(headers={}), dict(BATCH, name='c.sql'), 'SELECT', tmp_path, st, tmp_path/'s2.json')
    st = state(); st['jobs']['d'] = {'status': 'submitted', 'job_id': '42'}
    with pytest.raises(ValueError, match='line break'):
        s.fetch_batch(Client(payload()[:-1], headers={'Transfer-Encoding': 'chunked'}), dict(BATCH, name='d.sql'), 'SELECT', tmp_path, st, tmp_path/'s3.json')
    assert st['jobs']['d']['status'] == 'downloaded_invalid'


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_verifies_and_covers_66_centres():
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['centres'] == 66 and checked['batches'] == 11


@needs_plan
def test_submit_refuses_without_approval():
    checkpoint = p.OUTPUT/'results'/'status.json'
    had_results = (p.OUTPUT/'results').exists()
    before = s.sha(checkpoint) if checkpoint.exists() else None
    result = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_random_inner_rows.py'), '--submit'],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert result.returncode == 1
    assert ('not the exact approved inner random-point plan' if s.APPROVED_SHA is None else 'explicit user approval') in result.stderr
    assert (s.sha(checkpoint) if checkpoint.exists() else None) == before
    assert (p.OUTPUT/'results').exists() == had_results
