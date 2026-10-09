"""Bundle 3 stays disconnected even when the package is installed.

MultiNest's non-commercial sampler is not selected. PyMuPDF (AGPL) and
pymultinest are not imported. SWarp is not invoked. Unverified stellar
libraries and photo-z templates are not downloaded into this tree.
"""
import os
import re
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BAGPIPES_FIT = os.path.join(ROOT, "sed_adapters", "native", "bagpipes_fit.py")
SKIP_DIRS = {".git", "__pycache__", "tcl", "tk", "tcllib", "tklib", "tkimg"}
FORBIDDEN = re.compile(
    r"(?m)^[ \t]*(?:import|from)[ \t]+(?:fitz|pymupdf|pymultinest)\b"
    r"|importlib\.import_module\(\s*['\"](?:fitz|pymupdf|pymultinest)['\"]"
    r"|sampler\s*=\s*['\"]multinest['\"]"
    r"|which\(\s*['\"]swarp['\"]"
    r"|\[\s*['\"]swarp['\"]"
)
DATA_FILES = (
    "sed_fit/backends/fsps_backend.py",
    "sed_fit/backends/bagpipes_backend.py",
    "sed_adapters/native/eazy_package.py",
    "sed_adapters/native/prospector_fit.py",
    "sed_adapters/native/bagpipes_fit.py",
    "moving/assist_worker.py",
)


def _py_files():
    for dirpath, dirnames, files in os.walk(ROOT):
        dirnames[:] = [name for name in dirnames if name not in SKIP_DIRS]
        for name in files:
            if name.endswith(".py") and name != "test_bundle3_disconnected.py":
                yield os.path.join(dirpath, name)


def test_forbidden_imports_and_programs_are_absent():
    bad = []
    for path in _py_files():
        text = open(path, encoding="utf-8", errors="replace").read()
        if FORBIDDEN.search(text):
            bad.append(os.path.relpath(path, ROOT))
    assert bad == []


def test_unverified_libraries_are_not_downloaded():
    for rel in DATA_FILES:
        text = open(os.path.join(ROOT, rel), encoding="utf-8").read()
        assert "urlretrieve" not in text
        assert "git clone" not in text


def test_bagpipes_sampler_is_nautilus_and_multinest_raises():
    from sed_fit.backends.bagpipes_backend import nautilus_sampler

    assert nautilus_sampler(None) == "nautilus"
    assert nautilus_sampler("nautilus") == "nautilus"
    for name in ("multinest", "MultiNest", "pymultinest"):
        with pytest.raises(RuntimeError, match="MultiNest") as exc:
            nautilus_sampler(name)
        assert "not called" in str(exc.value)
        assert "analytic" in str(exc.value)
    backend = open(
        os.path.join(ROOT, "sed_fit", "backends", "bagpipes_backend.py"),
        encoding="utf-8",
    ).read()
    script = open(BAGPIPES_FIT, encoding="utf-8").read()
    assert 'sampler="multinest"' not in backend
    assert "sampler='multinest'" not in backend
    assert "sampler=sampler" in backend
    assert "sampler=sampler" in script


def test_fit_sed_rejects_multinest_before_bagpipes_runs(tmp_path, monkeypatch):
    (tmp_path / "bagpipes.py").write_text(
        "def galaxy(*args, **kwargs):\n"
        "    raise AssertionError('galaxy built')\n"
        "def fit(*args, **kwargs):\n"
        "    raise AssertionError('fit built')\n",
        encoding="utf-8",
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    sys.modules.pop("bagpipes", None)
    from sed_fit.backends.bagpipes_backend import BagpipesBackend
    try:
        with pytest.raises(RuntimeError, match="MultiNest"):
            BagpipesBackend().fit_sed(
                np.ones((1, 2)),
                np.full((1, 2), 0.05),
                np.array([0.4]),
                ["F435W", "F606W"],
                sampler="multinest",
            )
    finally:
        sys.modules.pop("bagpipes", None)


def test_script_rejects_multinest_before_fitting(tmp_path):
    fake = tmp_path / "fake"
    fake.mkdir()
    (fake / "bagpipes.py").write_text("x = 1\n", encoding="utf-8")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(fake) + os.pathsep + ROOT
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    proc = subprocess.run(
        [sys.executable, BAGPIPES_FIT],
        input='{"params": {"sampler": "multinest"}}',
        capture_output=True,
        text=True,
        cwd=ROOT,
        env=env,
    )
    assert proc.returncode == 2
    assert "MultiNest" in proc.stderr
    assert "not called" in proc.stderr
    assert "analytic" in proc.stderr
    assert "Traceback" not in proc.stderr
