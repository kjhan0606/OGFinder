import csv
import io
import math
from pathlib import Path
import re
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_s16a_fdfc_grid as pf
import prepare_s16a_random_candidate_rows as pr
import prepare_s16a_wide_extension as p
import submit_s16a_wide_extension as s

LENS = dict(lens_id='S23B_9999', center_id='LC-S23B_9999', ra_deg='150.0', dec_deg='2.0', z_lens='0.4', radius_arcsec='550.0')


def at(offset_arcmin):
    return 150.0+offset_arcmin/60/math.cos(math.radians(2.0)), 2.0


def random_payload(offsets=(120.0,), density=0.001, n_matched=None):
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(pr.OUTPUT_COLUMNS)
    for i, off in enumerate(offsets):
        ra, dec = at(off)
        row = {k: '' for k in pr.OUTPUT_COLUMNS}
        row.update(lens_id=LENS['lens_id'], object_id=str(1000+i), ra=repr(ra), dec=repr(dec), iadjust_density=str(density),
                   n_matched=str(len(offsets) if n_matched is None else n_matched))
        w.writerow([row[k] for k in pr.OUTPUT_COLUMNS])
    return stream.getvalue().encode()


def fdfc_payload(cells=((150, 0, 30, 29, 0),)):
    stream = io.StringIO(); w = csv.writer(stream); w.writerow(pf.OUTPUT_COLUMNS)
    total = sum(c[2] for c in cells)
    for x, y, n, f, null in cells:
        w.writerow([LENS['lens_id'], x, y, n, f, null, len(cells), total])
    return stream.getvalue().encode()


def test_queries_are_annulus_and_wide_cone():
    sql = p.random_query(LENS)
    assert sql.count('r.cx*') == 2 and f'< {math.cos(math.radians(1.5)):.17g}' in sql and f'>= {math.cos(math.radians(3.0)):.17g}' in sql
    assert sql.rstrip().endswith(f'LIMIT {p.RANDOM_LIMIT};') and 'r.iadjust_density < 0.003' in sql
    f = p.fdfc_query(LENS)
    assert 'coneSearch(f.coord, 150.0000000000, 2.0000000000, 10800.0)' in f and f.rstrip().endswith(f'LIMIT {p.FDFC_LIMIT};')
    assert not re.search(r'object_id|JOIN|shape', f)


def test_valid_batches(tmp_path):
    a = tmp_path/'r.csv'; a.write_bytes(random_payload((95.0, 179.0)))
    assert s.validate_random(a, dict(lens=LENS, row_limit=5))['rows'] == 2
    b = tmp_path/'f.csv'; b.write_bytes(fdfc_payload([(150, 0, 30, 29, 0), (-170, 3, 5, 0, 1)]))
    assert s.validate_fdfc(b, dict(lens=LENS, row_limit=5)) == dict(rows=2, n_primary=35, n_fdfc=29, n_fdfc_null=1)


@pytest.mark.parametrize('kind', ['inside_90', 'outside_180', 'density', 'cut', 'limit'])
def test_random_refusals(tmp_path, kind):
    args = dict(inside_90=dict(offsets=(80.0,)), outside_180=dict(offsets=(185.0,)), density=dict(density=0.004),
                cut=dict(n_matched=5), limit=dict(offsets=tuple(100.0+i for i in range(5))))[kind]
    f = tmp_path/'r.csv'; f.write_bytes(random_payload(**args))
    with pytest.raises(ValueError):
        s.validate_random(f, dict(lens=LENS, row_limit=5))


def test_fdfc_refuses_outside_cells(tmp_path):
    f = tmp_path/'f.csv'; f.write_bytes(fdfc_payload([(190, 0, 30, 29, 0)]))
    with pytest.raises(ValueError):
        s.validate_fdfc(f, dict(lens=LENS, row_limit=5))


needs_plan = pytest.mark.skipif(not (p.OUTPUT/'plan.json').is_file(), reason='pinned plan is local data')


@needs_plan
def test_plan_verifies():
    checked = p.validate(p.OUTPUT)
    assert s.APPROVED_SHA in (None, checked['plan_sha256']) and checked['batches'] == 8


@needs_plan
def test_submit_refuses_without_approval():
    had = (p.OUTPUT/'results').exists()
    r = subprocess.run([sys.executable, '-B', str(HERE/'submit_s16a_wide_extension.py'), '--submit'],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=HERE)
    assert r.returncode == 1
    assert ('not the exact approved' if s.APPROVED_SHA is None else 'explicit user approval') in r.stderr
    assert (p.OUTPUT/'results').exists() == had
