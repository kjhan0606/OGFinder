"""TinyTim and WebbPSF are optional and explicit.

Gates fixed before the run:
- no module-level ``import webbpsf`` and no ``pip install webbpsf``
- a missing webbpsf import raises "not found", writes no FITS, and does
  not leave the name in sys.modules
- a fake webbpsf module is what the JWST path runs, and the header says
  PSFMETH=webbpsf
- with tiny1/tiny2/tiny3 absent from PATH, check_tinytim_available()
  is False, generate_tinytim raises before writing, and the HST CLI
  exits non-zero without an output FITS
- with stub programs on PATH, the check is True, tiny1 and tiny2 are
  executed, and the written FITS header says PSFMETH=tinytim
- check_sim prints WEBBPSF=0 or 1 from the import, and TINYTIM from PATH
"""
import os
import re
import stat
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CLI = os.path.join(ROOT, "ds9", "library", "ds9_psf_deconv.py")
SOURCES = (
    os.path.join(ROOT, "psf_deconv", "psf", "sim_psf.py"),
    os.path.join(ROOT, "ds9", "library", "ds9_psf_deconv.py"),
    os.path.join(ROOT, "plugins", "star_psf", "star_psf.tcl"),
)
CLEAN_PATH = "/usr/bin:/bin"


def _env(path):
    env = dict(os.environ)
    env["PATH"] = path
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    return env


def _run(argv, path):
    return subprocess.run(
        argv, capture_output=True, text=True, cwd=ROOT, env=_env(path),
    )


def _install_stubs(directory):
    directory.mkdir()
    log = directory / "calls.txt"
    log.write_text("", encoding="utf-8")
    shebang = "#!%s\n" % sys.executable
    tiny2 = (
        shebang
        + "import os, sys\n"
        "open(os.environ['TINY_LOG'], 'a', encoding='utf-8').write('tiny2\\n')\n"
        "arg = sys.argv[1]\n"
        "base = arg[:-4] if arg.endswith('.tt3') else arg\n"
        "import numpy as np\n"
        "from astropy.io import fits\n"
        "fits.PrimaryHDU(np.ones((5, 5), np.float32)).writeto(\n"
        "    base + '00.fits', overwrite=True)\n"
    )
    bodies = {
        "tiny1": shebang + "import os\nopen(os.environ['TINY_LOG'], 'a', encoding='utf-8').write('tiny1\\n')\n",
        "tiny2": tiny2,
        "tiny3": shebang + "import os\nopen(os.environ['TINY_LOG'], 'a', encoding='utf-8').write('tiny3\\n')\n",
    }
    for name, body in bodies.items():
        path = directory / name
        path.write_text(body, encoding="utf-8")
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return log


def test_webbpsf_import_is_inside_the_function():
    level0 = re.compile(r"(?m)^import webbpsf\b|^from webbpsf\b")
    for path in SOURCES:
        text = open(path, encoding="utf-8").read()
        assert "pip install webbpsf" not in text
        assert level0.search(text) is None
    sim = open(SOURCES[0], encoding="utf-8").read()
    assert "import webbpsf" in sim
    tcl = open(SOURCES[2], encoding="utf-8").read()
    assert "TinyTim: Available" in tcl
    assert "TinyTim: Not found (tiny1/tiny2/tiny3 not on PATH)" in tcl
    assert "WebbPSF: Available" in tcl
    assert "WebbPSF: Not found" in tcl


def test_webbpsf_stays_disconnected_and_tinytim_needs_programs(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", CLEAN_PATH)
    sys.path.insert(0, ROOT)
    from psf_deconv.psf.sim_psf import (
        check_tinytim_available,
        check_webbpsf_available,
        generate_tinytim,
        generate_webbpsf,
    )

    assert check_webbpsf_available() is False
    assert check_tinytim_available() is False
    webb = tmp_path / "webb.fits"
    tiny = tmp_path / "tiny.fits"
    with pytest.raises(RuntimeError, match="not found") as webb_err:
        generate_webbpsf("NIRCAM", "F150W", output=str(webb))
    with pytest.raises(RuntimeError, match="not found on PATH") as tiny_err:
        generate_tinytim("ACS_WFC", "F814W", output=str(tiny))
    assert "Install" not in str(webb_err.value)
    assert not webb.exists()
    assert not tiny.exists()
    assert "produced" in str(tiny_err.value)


def test_stub_tinytim_is_executed_and_labeled(tmp_path, monkeypatch):
    log = _install_stubs(tmp_path / "bin")
    monkeypatch.setenv("PATH", str(tmp_path / "bin") + os.pathsep + CLEAN_PATH)
    monkeypatch.setenv("TINY_LOG", str(log))
    sys.path.insert(0, ROOT)
    from psf_deconv.psf.sim_psf import check_tinytim_available, generate_tinytim
    from astropy.io import fits

    assert check_tinytim_available() is True
    out = tmp_path / "hst.fits"
    psf = generate_tinytim("ACS_WFC", "F814W", psf_size=5, output=str(out))
    calls = log.read_text(encoding="utf-8").split()
    assert calls == ["tiny1", "tiny2"]
    assert psf.shape == (5, 5)
    assert abs(float(psf.sum()) - 1.0) < 1e-6
    assert out.is_file()
    header = fits.getheader(str(out))
    assert header["PSFMETH"] == "tinytim"
    assert header["TELESCOP"] == "HST"
    assert header["INSTRUME"] == "ACS_WFC"
    assert header["FILTER"] == "F814W"


def test_fresh_interpreter_does_not_import_webbpsf():
    script = (
        "import sys\n"
        "sys.path.insert(0, %r)\n"
        "import psf_deconv.psf.sim_psf as sim\n"
        "assert sim.check_webbpsf_available() is False\n"
        "assert 'webbpsf' not in sys.modules\n"
        "print('ok')\n"
    ) % ROOT
    proc = _run([sys.executable, "-c", script], CLEAN_PATH)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "ok"


def test_check_sim_and_cli_follow_the_program(tmp_path):
    proc = _run([sys.executable, CLI, "dummy.fits", "--mode", "check_sim"], CLEAN_PATH)
    assert proc.returncode == 0, proc.stderr
    assert "#SIM_STATUS\tWEBBPSF=0\tTINYTIM=0" in proc.stdout
    assert "WebbPSF: Not found" in proc.stderr
    assert "not on PATH" in proc.stderr
    assert "pip install" not in proc.stderr

    missing = tmp_path / "missing.fits"
    hst = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "hst",
        "--sim-instrument", "ACS_WFC",
        "--sim-filter", "F814W",
        "--psf-output", str(missing),
    ], CLEAN_PATH)
    assert hst.returncode != 0
    assert "not found on PATH" in hst.stderr
    assert not missing.exists()

    jwst = tmp_path / "jwst.fits"
    refused = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "jwst",
        "--sim-instrument", "NIRCAM",
        "--sim-filter", "F150W",
        "--psf-output", str(jwst),
    ], CLEAN_PATH)
    assert refused.returncode != 0
    assert "not found" in refused.stderr
    assert not jwst.exists()

    log = _install_stubs(tmp_path / "bin")
    found = _run(
        [sys.executable, CLI, "dummy.fits", "--mode", "check_sim"],
        str(tmp_path / "bin") + os.pathsep + CLEAN_PATH,
    )
    assert found.returncode == 0, found.stderr
    assert "#SIM_STATUS\tWEBBPSF=0\tTINYTIM=1" in found.stdout
    assert "WebbPSF: Not found" in found.stderr
    assert "TinyTim: available" in found.stderr
    assert "pip install" not in found.stderr
    assert log.read_text(encoding="utf-8") == ""


def _fake_webbpsf(directory):
    directory.mkdir()
    (directory / "webbpsf.py").write_text(
        "import os\n"
        "import numpy as np\n"
        "class _HDU:\n"
        "    def __init__(self, data):\n"
        "        self.data = data\n"
        "class NIRCam:\n"
        "    def __init__(self):\n"
        "        self.filter = None\n"
        "        self.options = {}\n"
        "        open(os.environ['WEBB_LOG'], 'w', encoding='utf-8').write('NIRCam')\n"
        "    def calc_psf(self, oversample=1, fov_pixels=201):\n"
        "        return [_HDU(np.full((int(fov_pixels), int(fov_pixels)), 4.0))]\n",
        encoding="utf-8",
    )


def test_fake_webbpsf_is_labeled_when_requested(tmp_path):
    fake = tmp_path / "fake"
    _fake_webbpsf(fake)
    log = tmp_path / "webb.log"
    env = _env(CLEAN_PATH)
    env["PYTHONPATH"] = str(fake) + os.pathsep + ROOT
    env["WEBB_LOG"] = str(log)
    out = tmp_path / "jwst.fits"
    proc = subprocess.run(
        [
            sys.executable, CLI, "dummy.fits",
            "--mode", "sim_psf",
            "--sim-telescope", "jwst",
            "--sim-instrument", "NIRCAM",
            "--sim-filter", "F150W",
            "--sim-psf-size", "5",
            "--psf-output", str(out),
        ],
        capture_output=True, text=True, cwd=ROOT, env=env,
    )
    assert proc.returncode == 0, proc.stderr[-800:]
    assert log.read_text(encoding="utf-8") == "NIRCam"
    assert out.is_file()
    from astropy.io import fits
    header = fits.getheader(str(out))
    assert header["PSFMETH"] == "webbpsf"
    assert header["TELESCOP"] == "JWST"
    data = fits.getdata(str(out))
    assert data.shape == (5, 5)
    assert abs(float(data.sum()) - 1.0) < 1e-6
    status = subprocess.run(
        [sys.executable, CLI, "dummy.fits", "--mode", "check_sim"],
        capture_output=True, text=True, cwd=ROOT, env=env,
    )
    assert status.returncode == 0, status.stderr
    assert "#SIM_STATUS\tWEBBPSF=1\tTINYTIM=0" in status.stdout
    assert "WebbPSF: Available" in status.stderr


def test_webbpsf_runtime_error_on_import_is_not_found(tmp_path):
    fake = tmp_path / "boom"
    fake.mkdir()
    (fake / "webbpsf.py").write_text("raise RuntimeError('boom')\n", encoding="utf-8")
    script = (
        "import sys\n"
        "sys.path.insert(0, %r)\n"
        "sys.path.insert(0, %r)\n"
        "import psf_deconv.psf.sim_psf as sim\n"
        "assert sim.check_webbpsf_available() is False\n"
        "print('ok')\n"
    ) % (str(fake), ROOT)
    proc = _run([sys.executable, "-c", script], CLEAN_PATH)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "ok"
