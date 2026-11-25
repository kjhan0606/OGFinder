import csv
import io
import math
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_random_source_rows as p
import submit_s16a_random_source_rows as s

C = [dict(center_id='RC-aaaa', lens_id='S23B_9999', ra_deg=150.0, dec_deg=2.0, z_lens=0.4, radius_arcsec=550.0, fdfc_fraction=1.0),
     dict(center_id='RC-bbbb', lens_id='S23B_9999', ra_deg=151.0, dec_deg=2.0, z_lens=0.4, radius_arcsec=550.0, fdfc_fraction=1.0)]


def at(c, offset_arcsec):
    return c['ra_deg']+offset_arcsec/3600/math.cos(math.radians(c['dec_deg'])), c['dec_deg']


def pair_payload(rows, n_batch=None, n_matched=None):
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(p.PAIR_COLUMNS)
    per = {}
    for cid, _, _ in rows:
        per[cid] = per.get(cid, 0)+1
    for i, (cid, off, obj) in enumerate(rows):
        ra, dec = at(next(c for c in C if c['center_id'] == cid), off)
        w.writerow([cid, obj, repr(ra), repr(dec), 0.1, -0.1, 0.3, 0.4, 2.0, -0.1, 0.0, 0.0,
                    per[cid] if n_matched is None else n_matched, len(rows) if n_batch is None else n_batch])
    return stream.getvalue().encode()


def inner_payload(rows):
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(p.INNER_COLUMNS)
    per = {}
    for cid, _, _ in rows:
        per[cid] = per.get(cid, 0)+1
    for cid, off, obj in rows:
        ra, dec = at(next(c for c in C if c['center_id'] == cid), off)
        row = {k: '' for k in p.INNER_COLUMNS}
        row.update(center_id=cid, object_id=str(obj), ra=repr(ra), dec=repr(dec), n_matched=str(per[cid]), n_batch=str(len(rows)))
        w.writerow([row[k] for k in p.INNER_COLUMNS])
    return stream.getvalue().encode()


def pair_batch(limit=5):
    return dict(centres=C, row_limit_per_centre=limit)


def inner_batch(limit=5):
    return dict(centres=[p.inner_centre(c) for c in C], row_limit_per_centre=limit)


GOOD = [('RC-aaaa', 10.0, 11), ('RC-aaaa', 500.0, 12), ('RC-bbbb', 0.0, 11)]


def test_pair_query_wraps_pilot_form():
    sql = p.pair_query(C)
    assert sql.count('SELECT DISTINCT') == 2 and sql.count('count(*) OVER () AS n_matched') == 2
    assert sql.count(f'LIMIT {p.PAIR_LIMIT_PER_CENTRE}') == 2 and 'AS n_batch' in sql and p.NEW_WEIGHT in sql and p.OLD_WEIGHT+' ' not in sql
    assert 'coneSearch(c.icoord, 150, 2, 550)' in sql


def test_valid_batches(tmp_path):
    f = tmp_path/'p.csv'; f.write_bytes(pair_payload(GOOD))
    assert s.validate_pairs(f, pair_batch())['rows_by_centre'] == {'RC-aaaa': 2, 'RC-bbbb': 1}
    g = tmp_path/'i.csv'; g.write_bytes(inner_payload([('RC-aaaa', 1.0, 5), ('RC-bbbb', 2.0, 5)]))
    assert s.validate_inner(g, inner_batch())['rows'] == 2


@pytest.mark.parametrize('kind', ['outside', 'duplicate', 'empty_centre', 'limit', 'n_batch', 'n_matched', 'foreign'])
def test_pair_refusals(tmp_path, kind):
    rows, extra = list(GOOD), {}
    if kind == 'outside':
        rows[1] = ('RC-aaaa', 551.0, 12)
    elif kind == 'duplicate':
        rows[1] = ('RC-aaaa', 20.0, 11)
    elif kind == 'empty_centre':
        rows = rows[:2]
    elif kind == 'limit':
        rows += [('RC-bbbb', 1.0, 20+i) for i in range(5)]
    elif kind == 'n_batch':
        extra = dict(n_batch=99)
    elif kind == 'n_matched':
        extra = dict(n_matched=7)
    else:
        rows[2] = ('RC-cccc', 0.0, 11)
    f = tmp_path/'p.csv'
    if kind == 'foreign':
        C.append(dict(C[1], center_id='RC-cccc'))
        try:
            f.write_bytes(pair_payload(rows))
        finally:
            C.pop()
    else:
        f.write_bytes(pair_payload(rows, **extra))
    with pytest.raises(ValueError):
        s.validate_pairs(f, pair_batch())


def test_inner_refuses_outside_and_empty(tmp_path):
    radius = p.inner_centre(C[0])['query_radius_arcsec']
    for rows in ([('RC-aaaa', radius+1, 5), ('RC-bbbb', 0.0, 6)], [('RC-aaaa', 1.0, 5)]):
        f = tmp_path/'i.csv'; f.write_bytes(inner_payload(rows))
        with pytest.raises(ValueError):
            s.validate_inner(f, inner_batch())


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_verifies():
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['batches'] == 22
    assert checked['estimate']['csv_bytes'] < p.CAP


@needs_plan
def test_submit_refuses_without_approval():
    had = (p.OUTPUT/'results').exists()
    r = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_random_source_rows.py'), '--submit'],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert r.returncode == 1
    assert ('not the exact approved' if s.APPROVED_SHA is None else 'explicit user approval') in r.stderr
    assert (p.OUTPUT/'results').exists() == had
