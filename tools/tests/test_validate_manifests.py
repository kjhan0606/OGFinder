"""tools/validate_manifests.py: the shipped manifests are clean and each kind of mistake in a cli template is caught."""
import copy, json, os, shutil, subprocess, sys
import pytest
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import validate_manifests as V


def test_shipped_manifests_are_clean():
    probs, st = V.validate(os.path.join(ROOT, "plugins"))
    assert probs == []
    assert st["cli"] >= 6 and st["plugins"] >= 18


def _tree(tmp_path, mutate):
    dst = tmp_path / "plugins"; shutil.copytree(os.path.join(ROOT, "plugins"), dst, ignore=shutil.ignore_patterns("__pycache__", "tests"))
    p = dst / "morphology" / "plugin.json"; d = json.load(open(p))
    mutate(d); json.dump(d, open(p, "w"))
    return str(dst)


def _step(d, sid):
    return next(s for s in d["steps"] if s["id"] == sid)


@pytest.mark.parametrize("name,mut,needle", [
    ("unknown token", lambda d: _step(d, "sersic")["cli"].append("{no-such-param}"), "unknown token"),
    ("other plugin param missing", lambda d: _step(d, "sersic")["cli"].append("{extract:nope}"), "has no parameter"),
    ("no such plugin", lambda d: _step(d, "sersic")["cli"].append("{ghost:x}"), "no plugin ghost"),
    ("missing script", lambda d: _step(d, "sersic")["cli"].__setitem__(1, "{script:ds9_missing.py}"), "not found"),
    ("typo in flag", lambda d: _step(d, "sersic")["cli"].__setitem__(5, "--mag-zeropointt"), "does not occur"),
    ("catalog without needs", lambda d: _step(d, "sersic").__setitem__("needs", ["image"]), "{catalog} without needs"),
    ("image without needs", lambda d: _step(d, "sersic").__setitem__("needs", ["catalog"]), "{image} without needs"),
    ("cli and proc", lambda d: _step(d, "sersic").__setitem__("proc", "CatalogPanelSersicFit"), "use one"),
    ("condition on missing param", lambda d: _step(d, "bulge_disk")["cli"].append({"if": "ghost", "argv": ["--x"]}), "does not exist"),
    ("bad output mode", lambda d: _step(d, "sersic")["output"].__setitem__("mode", "bogus"), "output.mode"),
    ("undefined proc", lambda d: _step(d, "morphometry").update({"proc": "NoSuchProc", "cli": None}) or _step(d, "morphometry").pop("cli"), "not defined"),
    ("duplicate param", lambda d: d["params"].append(copy.deepcopy(d["params"][0])), "duplicate parameter"),
    ("default outside range", lambda d: d["params"][0].update(default=10 ** 9), "> max"),
])
def test_mistakes_are_reported(tmp_path, name, mut, needle):
    probs, _ = V.validate(_tree(tmp_path, mut))
    assert any(needle in p for p in probs), (name, probs)


def test_cli_exit_status(tmp_path):
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "validate_manifests.py")], capture_output=True, text=True)
    assert r.returncode == 0 and "0 problem(s)" in r.stdout
    bad = _tree(tmp_path, lambda d: _step(d, "sersic")["cli"].append("{zzz}"))
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "validate_manifests.py"), bad], capture_output=True, text=True)
    assert r.returncode == 1 and "PROBLEM" in r.stdout
