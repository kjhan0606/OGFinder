"""The GALFIT binary is optional and is not bundled.

Gates fixed before the run:
- an empty or missing path raises FileNotFoundError and does not create
  a g/ directory
- a stub executable is run, and its galfit.01 table is stored under the
  key galfit, with the Sersic type and Chi^2/nu taken from that file
- the stub result has no ogfmeas key
"""
import os
import stat
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
VAL = os.path.join(ROOT, "plugins", "multifit", "validation")
sys.path.insert(0, VAL)

import galfit_real_compare as GRC  # noqa: E402

FEED = """\
A) data.fits
B) out.fits
J) 26.0

0) sersic
1) 10.0 12.0 1 1
3) 20.0 1
4) 4.0 1
5) 1.5 1
9) 0.8 1
10) 20.0 1
0) sky
1) 1.0 1
2) 0.0 0
3) 0.0 0
Chi^2/nu = 1.23
"""


class _Args:
    def __init__(self, galfit, ld=""):
        self.galfit = galfit
        self.ld = ld


def test_missing_binary_writes_nothing(tmp_path):
    md = tmp_path / "cut"
    md.mkdir()
    (md / "start.feedme").write_text("A) data.fits\n", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="not found"):
        GRC.run_galfit_optional(_Args(""), str(md), 30)
    assert GRC._galfit_bin(_Args(str(tmp_path / "missing-galfit"))) == ""
    assert not (md / "g").exists()


def test_stub_binary_is_called_and_labeled_galfit(tmp_path):
    md = tmp_path / "cut"
    md.mkdir()
    (md / "start.feedme").write_text("A) data.fits\n", encoding="utf-8")
    stub = tmp_path / "galfit"
    stub.write_text(
        "#!%s\n" % sys.executable
        + "import pathlib, sys\n"
        "here = pathlib.Path.cwd()\n"
        "(here / 'called.txt').write_text(sys.argv[0], encoding='utf-8')\n"
        "(here / 'galfit.01').write_text(%r, encoding='utf-8')\n" % FEED,
        encoding="utf-8",
    )
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    out = GRC.run_galfit_optional(_Args(str(stub)), str(md), 30)
    assert out["galfit"]["comps"][0]["type"] == "sersic"
    assert out["galfit"]["chi2nu"] == pytest.approx(1.23)
    assert "ogfmeas" not in out
    assert (md / "g" / "called.txt").read_text(encoding="utf-8") == str(stub)
