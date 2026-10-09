"""Explicit SED names may import a user-installed package. auto stays analytic.

The spectrum gate was fixed before the run. It is the same Gaussian as
test_template_redshift_and_circular_orbit: rest 5000 A, z = 0.4, tolerance 0.02.
Pivot-wavelength samples are a coarse stand-in for a filter integral. One
scale on a few bands does not pick out a unique redshift, so that path is
checked only for the label and for a redshift on the requested grid.
"""
import json
import os
import re
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HELPER = os.path.join(ROOT, "sed_adapters", "native", "eazy_native.py")
BAGPIPES_FIT = os.path.join(ROOT, "sed_adapters", "native", "bagpipes_fit.py")
PROSPECTOR_FIT = os.path.join(ROOT, "sed_adapters", "native", "prospector_fit.py")

# Spectrum: rest-frame Gaussian at 5000 A, planted at z = 0.4.
Z_TRUE = 0.4
Z_TOL = 0.02
# Native band path: the caller asks for this grid. The returned redshift has
# to be one of its nodes.
Z_MIN = 0.0
Z_MAX = 0.6
Z_STEP = 0.02

COLUMN0_RE = re.compile(
    r"(?m)^(?:import|from)[ \t]+(?:pcigale|fsps|bagpipes|eazy|prospect|sedpy|dense_basis)(?:\s|\.|$)"
)
BANNED_ANY_RE = re.compile(
    r"(?m)^[ \t]*(?:import|from)[ \t]+(?:pcigale|dense_basis)(?:\s|\.|$)"
)
INDENTED_RE = re.compile(
    r"(?m)^[ \t]+(?:import|from)[ \t]+(?:fsps|bagpipes|eazy|prospect|sedpy)(?:\s|\.|$)"
)
INDENTED_OK = {
    "sed_fit/backends/fsps_backend.py",
    "sed_fit/backends/bagpipes_backend.py",
    "sed_fit/backends/prospector_backend.py",
    "sed_adapters/native/bagpipes_fit.py",
    "sed_adapters/native/prospector_fit.py",
    "sed_adapters/native/eazy_package.py",
}
BANNED = ("pcigale", "fsps", "bagpipes", "eazy", "prospect", "sedpy", "dense_basis")


def _run(argv, stdin=None):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    return subprocess.run(
        argv, input=stdin, capture_output=True, text=True, cwd=ROOT, env=env,
    )


def test_sources_do_not_import_the_external_packages():
    roots = [
        os.path.join(ROOT, "sed_fit", "backends"),
        os.path.join(ROOT, "sed_adapters"),
        os.path.join(ROOT, "plugins", "sedcodes"),
    ]
    bad = []
    for base in roots:
        for dirpath, _, files in os.walk(base):
            for name in files:
                if not name.endswith(".py") or name == "test_real_codes.py":
                    continue
                path = os.path.join(dirpath, name)
                text = open(path, encoding="utf-8").read()
                rel = os.path.relpath(path, ROOT)
                if BANNED_ANY_RE.search(text) or COLUMN0_RE.search(text):
                    bad.append(rel)
                elif INDENTED_RE.search(text) and rel not in INDENTED_OK:
                    bad.append(rel)
    assert bad == []
    adapter = open(os.path.join(ROOT, "sed_adapters", "eazy_adapter.py"), encoding="utf-8").read()
    script = open(os.path.join(ROOT, "sed_adapters", "script_adapter.py"), encoding="utf-8").read()
    assert "import eazy" not in adapter
    assert "import %s" not in script


def test_auto_backend_is_analytic_and_does_not_import():
    script = r"""
import os
import sys
os.environ["PATH"] = "/usr/bin:/bin"
os.environ.pop("SPS_HOME", None)
os.environ.pop("FSPS_DIR", None)
sys.path.insert(0, %r)
for name in %r:
    sys.modules.pop(name, None)
import numpy as np
from sed_fit.backends import get_backend, list_available
backend = get_backend("auto")
assert backend.name() == "analytic", backend.name()
mags = backend.generate_sed(0.4, 10.0, 9.0, -0.3, 0.2, 9.0, ["F435W", "F606W", "F814W"])
assert len(mags) == 3 and np.all(np.isfinite(mags))
loaded = [name for name in %r if name in sys.modules]
assert not loaded, loaded
from sed_fit.backends.fsps_backend import FSPSBackend
from sed_fit.backends.bagpipes_backend import BagpipesBackend
from sed_fit.backends.prospector_backend import ProspectorBackend
from sed_fit.backends.cigale_backend import CIGALEBackend
from sed_fit.backends.dense_basis_backend import DenseBasisBackend
for cls in (FSPSBackend, BagpipesBackend, ProspectorBackend, CIGALEBackend, DenseBasisBackend):
    assert cls.is_available() is False, cls.name()
    try:
        cls().generate_sed(0.4, 10.0, 9.0, 0.0, 0.2, 9.0, ["F606W"])
    except RuntimeError as exc:
        text = str(exc)
        assert "analytic" in text
        assert "Install" not in text
    else:
        raise SystemExit("generate_sed returned for " + cls.name())
    try:
        cls().fit_sed(np.array([[22.0, 21.0]]), np.array([[0.05, 0.05]]), np.array([0.4]), ["F435W", "F606W"])
    except RuntimeError as exc:
        text = str(exc)
        assert "analytic" in text
        assert "Install" not in text
    else:
        raise SystemExit("fit_sed returned for " + cls.name())
avail = list_available()
assert "analytic" in avail
for name in ("fsps", "bagpipes", "prospector", "cigale", "dense_basis"):
    assert name not in avail
    try:
        get_backend(name)
    except ValueError as exc:
        text = str(exc)
        assert "analytic" in text
        assert "Install the required package" not in text
        if name == "dense_basis":
            assert "does not call" in text
        else:
            assert "not found" in text
    else:
        raise SystemExit(name)
from sed_adapters import eazy_adapter, bagpipes_adapter, prospector_adapter
assert eazy_adapter.check()["engine_native"] is False
assert eazy_adapter.check()["version"] == "ogfmeas.photoz"
assert bagpipes_adapter.check()["engine_native"] is False
assert prospector_adapter.check()["engine_native"] is False
loaded = [name for name in %r if name in sys.modules]
assert not loaded, loaded
print("ok")
""" % (ROOT, BANNED, BANNED, BANNED)
    proc = _run([sys.executable, "-c", script])
    assert proc.returncode == 0, proc.stderr[-800:]
    assert proc.stdout.strip().splitlines()[-1] == "ok"


def test_listed_backends_mark_optional_names():
    env = dict(os.environ)
    env["PATH"] = "/usr/bin:/bin"
    env.pop("SPS_HOME", None)
    env.pop("FSPS_DIR", None)
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    proc = subprocess.run(
        [
            sys.executable,
            os.path.join(ROOT, "ds9", "library", "ds9_sed_fit.py"),
            "--mode", "list-backends",
        ],
        capture_output=True, text=True, cwd=ROOT, env=env,
    )
    assert proc.returncode == 0, proc.stderr[-400:]
    rows = dict(line.split("\t") for line in proc.stdout.splitlines() if "\t" in line)
    assert rows["analytic"] == "available"
    for name in ("fsps", "bagpipes", "prospector", "cigale"):
        assert rows[name] == "not_found", (name, rows[name])
    assert rows["dense_basis"] == "not_called"


def test_template_helper_recovers_the_planted_line():
    wave = np.linspace(3000.0, 9000.0, 240)
    rest = np.exp(-0.5 * ((wave - 5000.0) / 80.0) ** 2)
    observed = np.interp(wave, wave * (1.0 + Z_TRUE), rest, left=0.0, right=0.0)
    req = {
        "ids": ["1"],
        "wavelength": wave.tolist(),
        "flux": observed.tolist(),
        "err": [0.02] * int(wave.size),
        "templates": [[wave.tolist(), rest.tolist()]],
        "redshifts": np.linspace(0.0, 1.0, 51).tolist(),
    }
    proc = _run([sys.executable, HELPER], json.dumps(req))
    assert proc.returncode == 0, proc.stderr[-400:]
    out = json.loads(proc.stdout)
    assert out["version"] == "ogfmeas.photoz"
    assert abs(out["rows"][0]["z"] - Z_TRUE) <= Z_TOL


def test_native_engine_uses_pivot_samples_when_templates_are_supplied():
    from sed_adapters import eazy_adapter
    wave = np.linspace(2000.0, 12000.0, 2000)
    rest = np.exp(-0.5 * ((wave - 5000.0) / 40.0) ** 2)
    bands = ["F435W", "F606W", "F775W", "F850LP"]
    recs = [{
        "id": "1",
        "mags": {"F435W": 28.0, "F606W": 22.0, "F775W": 28.0, "F850LP": 28.0},
        "mag_errs": {b: 0.05 for b in bands},
    }]
    res, model = eazy_adapter.process(recs, {
        "engine": "native",
        "templates": [[wave.tolist(), rest.tolist()]],
        "z_min": Z_MIN,
        "z_max": Z_MAX,
        "z_step": Z_STEP,
        "min_mag_err": 0.0,
        "n_min_colors": 3,
    }, "photoz")
    assert model == "ogfmeas.photoz"
    assert "EAZY" not in model
    assert "error" not in res[0]
    z = res[0]["PHOTOZ"]
    n = int(round((Z_MAX - Z_MIN) / Z_STEP))
    assert abs(z - (Z_MIN + round((z - Z_MIN) / Z_STEP) * Z_STEP)) < 1e-8
    assert Z_MIN - 1e-8 <= z <= Z_MIN + n * Z_STEP + 1e-8


def test_native_engine_stops_when_templates_are_absent():
    from sed_adapters import eazy_adapter
    bands = ["F435W", "F606W", "F775W", "F850LP"]
    recs = [{
        "id": "1",
        "mags": {b: 23.0 for b in bands},
        "mag_errs": {b: 0.05 for b in bands},
    }]
    with pytest.raises(RuntimeError, match="not called"):
        eazy_adapter.process(recs, {
            "engine": "native",
            "templates_file": "/no/such/templates.param",
            "filters_res": "/no/such/FILTER.RES",
            "z_min": 0.0,
            "z_max": 1.0,
            "z_step": 0.05,
        }, "photoz")


def test_cigale_program_is_optional(tmp_path, monkeypatch):
    """pcigale is a separate program. It is not imported."""
    import numpy as np
    from sed_fit.backends import get_backend
    from sed_fit.backends.cigale_backend import CIGALEBackend

    source = open(
        os.path.join(ROOT, "sed_fit", "backends", "cigale_backend.py"),
        encoding="utf-8",
    ).read()
    assert "import pcigale" not in source
    assert "from pcigale" not in source

    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    assert CIGALEBackend.is_available() is False
    with pytest.raises(RuntimeError, match="not found"):
        CIGALEBackend().fit_sed(
            np.array([[22.0, 21.0]]),
            np.array([[0.05, 0.05]]),
            np.array([0.4]),
            ["F435W", "F606W"],
        )
    with pytest.raises(ValueError, match="not found"):
        get_backend("cigale")
    assert get_backend("auto").name() == "analytic"

    stub = tmp_path / "pcigale"
    stub.write_text(
        "#!%s\n" % sys.executable
        + "import os, sys\n"
        "open('called.txt', 'w', encoding='utf-8').write(' '.join(sys.argv))\n"
        "os.makedirs('out', exist_ok=True)\n"
        "open('out/results.txt', 'w', encoding='utf-8').write(\n"
        "    'id stellar.m_star attenuation.Av best.chi2\\n'\n"
        "    '0 1e10 0.2 1.5\\n')\n"
    )
    stub.chmod(stub.stat().st_mode | 0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + "/usr/bin:/bin")
    assert CIGALEBackend.is_available() is True
    assert get_backend("auto").name() == "analytic"
    assert get_backend("cigale").name() == "cigale"
    fitted = CIGALEBackend().fit_sed(
        np.array([[22.0, 21.0]]),
        np.array([[0.05, 0.05]]),
        np.array([0.4]),
        ["F435W", "F606W"],
    )
    assert fitted.log_mass[0] == pytest.approx(10.0)
    with pytest.raises(RuntimeError, match="in-process"):
        CIGALEBackend().generate_sed(0.4, 10.0, 9.0, 0.0, 0.2, 9.0, ["F606W"])


@pytest.mark.parametrize("path", [BAGPIPES_FIT, PROSPECTOR_FIT])
def test_shipped_sed_scripts_refuse_without_importing(path):
    proc = _run([sys.executable, path], "{}")
    assert proc.returncode == 2
    assert "not called" in proc.stderr
    assert "not found" in proc.stderr
    assert "analytic" in proc.stderr
    assert "Traceback" not in proc.stderr
    assert "No module named" not in proc.stderr
