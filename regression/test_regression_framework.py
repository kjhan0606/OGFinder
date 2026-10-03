"""Offline tests of the regression-set machinery (no data download): range comparison, baseline coverage, skip-on-unavailable behaviour."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from regression import run_regression as R, cases as C, datasets as D  # noqa: E402


def test_compare_ranges():
    ok, bad = R.compare(dict(a=1.0, b=5.0), dict(a=[0.5, 1.5], b=[0, 4]))
    assert not ok and len(bad) == 1 and 'b' in bad[0]
    ok, bad = R.compare(dict(a=float('nan')), dict(a=[0, 1]))
    assert not ok and 'missing' in bad[0]
    ok, bad = R.compare(dict(), dict(a=[0, 1]))
    assert not ok
    assert R.compare(dict(a=0.5), dict(a=[0, 1]))[0]


def test_every_case_has_baseline_with_ranges():
    base = json.load(open(R.BASE))
    assert set(C.CASES) == set(base), set(C.CASES) ^ set(base)
    for n, b in base.items():
        assert b['ranges'], n
        for k, (lo, hi) in b['ranges'].items():
            assert lo <= hi, (n, k)


def test_covered_tools_documented():
    txt = open(os.path.join(os.path.dirname(HERE), 'docs', 'testing.md')).read()
    for n in C.CASES:
        assert n in txt, n


def test_unavailable_data_is_skip(tmp_path, monkeypatch):
    monkeypatch.setenv('OGF_DATA_CACHE', str(tmp_path))
    monkeypatch.setattr(D, 'CACHE', str(tmp_path))
    monkeypatch.setattr(D, 'LOCAL_HUDF', [])

    def boom(*a, **k):
        raise D.Unavailable('offline')
    monkeypatch.setattr(D, '_download', boom)
    rc = R.main(['--only', 'extract_hudf', '--report', str(tmp_path / 'r.json')])
    rep = json.load(open(tmp_path / 'r.json'))
    assert rc == 0 and rep['extract_hudf']['status'] == 'SKIP'
