"""Explicit optional libraries. The default path stays in-tree.

Gates fixed before the run:
- measurement_library() is ogfmeas unless OGF_USE_SEP=1
- OGF_USE_SEP=1 with sep blocked raises "not found" and does not return ogfmeas
- a default galaxy stamp matches _make_numpy even when galsim imports
- backend='galsim' with the package blocked raises "not found"
- backend='galsim' with a fake module returns that module's array
- get_backend('auto') stays analytic when a fake fsps is importable
- an explicit fsps backend returns magnitudes from that fake
- orbit.py has no rebound or assist import; the worker does
- propagate_assist raises "not found" and returns no ASSIST-labeled states
- engine=package with eazy missing raises "not found"
"""
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ORBIT = os.path.join(ROOT, "moving", "orbit.py")
WORKER = os.path.join(ROOT, "moving", "assist_worker.py")
PROFILES = os.path.join(ROOT, "ai_merge", "simulation", "profiles.py")
FORCE_MODEL = (
    "ogfmeas(builtin planets and Moon, solar Schwarzschild, "
    "Sun J2, Earth J2/J3/J4)"
)


def _run(argv, stdin=None):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    env["PYTHONPATH"] = ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        argv, input=stdin, capture_output=True, text=True, cwd=ROOT, env=env,
    )


def test_measurement_library_defaults_to_ogfmeas(monkeypatch):
    monkeypatch.delenv("OGF_USE_SEP", raising=False)
    sys.modules.pop("sep", None)
    import ogfmeas
    got = ogfmeas.measurement_library()
    assert got is ogfmeas
    assert "sep" not in sys.modules


def test_measurement_library_uses_installed_sep_when_asked(monkeypatch):
    monkeypatch.setenv("OGF_USE_SEP", "1")
    import ogfmeas
    got = ogfmeas.measurement_library()
    assert got.__name__ == "sep"


def test_measurement_library_errors_when_sep_is_blocked(monkeypatch):
    monkeypatch.setenv("OGF_USE_SEP", "1")
    monkeypatch.setitem(sys.modules, "sep", None)
    import ogfmeas
    with pytest.raises(RuntimeError, match="not found") as exc:
        got = ogfmeas.measurement_library()
        assert got is not ogfmeas
    assert "ogfmeas" in str(exc.value)


def test_default_stamp_matches_numpy():
    text = open(PROFILES, encoding="utf-8").read()
    assert "import galsim\n" not in text.split("def _make_galsim")[0]
    from ai_merge.simulation.profiles import _make_numpy, make_galaxy_stamp
    kw = dict(
        stamp_size=15, flux=100.0, re=2.0, sersic_n=1.5,
        ellip=0.2, theta_rad=0.3, psf_fwhm=1.5,
        x_offset=0.0, y_offset=0.0,
    )
    assert np.allclose(make_galaxy_stamp(**kw), _make_numpy(**kw))


def test_explicit_galsim_missing_is_not_found(monkeypatch):
    monkeypatch.setitem(sys.modules, "galsim", None)
    from ai_merge.simulation.profiles import make_galaxy_stamp
    with pytest.raises(RuntimeError, match="not found"):
        make_galaxy_stamp(
            15, 100.0, 2.0, 1.5, 0.2, 0.3, 1.5, backend="galsim",
        )


def test_explicit_galsim_returns_the_fake_array(monkeypatch):
    import types
    sentinel = np.arange(9, dtype=float).reshape(3, 3)

    class _Obj:
        def shear(self, **kwargs):
            return self

        def shift(self, *args):
            return self

        def drawImage(self, image=None, method=None):
            image.array[:] = sentinel
            return image

    class _Img:
        def __init__(self, *args, **kwargs):
            self.array = np.zeros((3, 3))

    fake = types.ModuleType("galsim")
    fake.Sersic = lambda **kwargs: _Obj()
    fake.Gaussian = lambda **kwargs: _Obj()
    fake.Convolve = lambda items: _Obj()
    fake.Image = _Img
    monkeypatch.setitem(sys.modules, "galsim", fake)
    from ai_merge.simulation.profiles import make_galaxy_stamp
    out = make_galaxy_stamp(3, 10.0, 1.0, 1.0, 0.0, 0.0, 1.0, backend="galsim")
    assert np.array_equal(out, sentinel)


def test_explicit_fsps_uses_the_fake_and_auto_stays_analytic(tmp_path, monkeypatch):
    (tmp_path / "fsps.py").write_text(
        "import numpy as np\n"
        "class StellarPopulation:\n"
        "    def __init__(self, **kwargs):\n"
        "        self.params = {}\n"
        "    def get_mags(self, tage=None, redshift=None, bands=None):\n"
        "        return np.zeros(len(bands)) + 20.0\n",
        encoding="utf-8",
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.delenv("SPS_HOME", raising=False)
    monkeypatch.delenv("FSPS_DIR", raising=False)
    sys.modules.pop("fsps", None)
    from sed_fit.backends import get_backend
    from sed_fit.backends.fsps_backend import FSPSBackend
    FSPSBackend._sp = None
    try:
        assert FSPSBackend.is_available() is True
        assert get_backend("auto").name() == "analytic"
        backend = get_backend("fsps")
        assert backend.name() == "fsps"
        mags = backend.generate_sed(0.4, 10.0, 9.0, 0.0, 0.2, 9.0, ["F435W"])
        assert mags[0] == pytest.approx(-5.0)
    finally:
        FSPSBackend._sp = None
        sys.modules.pop("fsps", None)


def test_parent_orbit_does_not_import_rebound():
    text = open(ORBIT, encoding="utf-8").read()
    assert "import rebound" not in text
    assert "import assist" not in text
    worker = open(WORKER, encoding="utf-8").read()
    assert "import rebound" in worker
    assert "import assist" in worker


def test_propagate_assist_reports_not_found():
    from moving.orbit import Propagator, assist_available, propagate_assist
    assert assist_available() is True
    prop = Propagator()
    assert prop.force_model == FORCE_MODEL
    with pytest.raises(RuntimeError, match="not found"):
        propagate_assist([[1.0, 0.0, 0.0, 0.0, 0.01, 0.0]], 2450000.0, [2450001.0])
    proc = _run([sys.executable, "-m", "moving.assist_worker"], "{}")
    assert proc.returncode == 2
    assert "not found" in proc.stderr
    assert "ASSIST" not in proc.stdout
    assert "Traceback" not in proc.stderr


def test_eazy_package_engine_errors_when_missing():
    from sed_adapters import eazy_adapter
    bands = ["F435W", "F606W", "F775W", "F850LP"]
    recs = [{
        "id": "1",
        "mags": {b: 23.0 for b in bands},
        "mag_errs": {b: 0.05 for b in bands},
    }]
    with pytest.raises(RuntimeError, match="not found"):
        eazy_adapter.process(recs, {"engine": "package"}, "photoz")
    native = open(
        os.path.join(ROOT, "sed_adapters", "native", "eazy_native.py"),
        encoding="utf-8",
    ).read()
    assert "import eazy" not in native
