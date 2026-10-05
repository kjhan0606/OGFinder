"""Stand-alone shell of the shared Astrafex Core: manifest/licence tags, launcher, and the command line against independent numpy references.
The compute core (package ds10core of the Astrafex Web repository) is optional here: tests that need it are skipped when it cannot be found.
Nothing in this file imports ds10core: it is only ever run as a separate process, exactly like the plugin does."""
import csv
import hashlib
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

import ds10 as launcher

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(PLUGIN))
CORE, WHY = launcher.find_core()
needs_core = pytest.mark.skipif(CORE is None, reason="ds10core (ds10-web/server) not found: " + str(WHY))


def run(*args, check=True, env=None):
    e = dict(os.environ, **(env or {}))
    r = subprocess.run([sys.executable, os.path.join(PLUGIN, "ds10.py"), *map(str, args)], capture_output=True, text=True, env=e)
    if check:
        assert r.returncode == 0, (r.returncode, r.stdout[-800:], r.stderr[-1500:])
    return r


# ------------------------------------------------------------------------------------------------ manifest + launcher (no core needed)
def test_manifest_has_licence_tags_and_is_stand_alone_only():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json")))
    lic = m["license"]
    assert lic["ds9_code_in_core"] is False and "separate process" in lic["boundary"] and "NOT GPL" in lic["core"]
    assert m["web"] is False                                                       # the web app does not list its own stand-alone shell as a tool
    assert {s["id"] for s in m["steps"]} >= {"info", "calib", "maps", "regions-stats", "regions-mask", "regions-convert", "pixtab", "forced", "run-script"}
    src = open(os.path.join(PLUGIN, "ds10.py")).read()
    assert "import ds10core" not in src and "from ds10core" not in src                # process boundary: the launcher never imports the core


def test_no_ds9_ast_sep_code_is_imported_by_the_plugin():
    for f in os.listdir(PLUGIN):
        if f.endswith((".py", ".tcl")):
            txt = open(os.path.join(PLUGIN, f)).read().lower()
            for bad in ("import sep", "from sep ", "import ogfkit", "from ogfkit", "libast", "starlink"):
                assert bad not in txt, (f, bad)


def test_launcher_finds_core_by_env_and_reports_missing_core(tmp_path):
    fake = tmp_path / "core"
    (fake / "ds10core").mkdir(parents=True)
    (fake / "ds10core" / "__init__.py").write_text("")
    r = run("--where", env={"DS10_CORE": str(fake)})
    assert json.loads(r.stdout)["core_dir"] == str(fake)
    r = run("--where", "--core-dir", str(tmp_path / "nope"), env={"DS10_CORE": str(tmp_path / "nope2"), "HOME": str(tmp_path)}, check=False)
    # explicit directory that does not exist falls through to the next candidate; with none found the exit status is 3 and the message lists what was tried
    if CORE is None:
        assert r.returncode == 3 and "ds10core not found" in r.stderr


def test_split_words_expand():
    out = subprocess.run([sys.executable, "-c", "import ds10,sys; sys.argv=['x']; print(ds10.__doc__ is not None)"], cwd=PLUGIN, capture_output=True, text=True)
    assert out.stdout.strip() == "True"


# ------------------------------------------------------------------------------------------------ synthetic data
def write_fits(path, data, **cards):
    from astropy.io import fits
    h = fits.PrimaryHDU(np.asarray(data, np.float32))
    h.header.update({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRPIX1": data.shape[1] / 2 + 0.5, "CRPIX2": data.shape[0] / 2 + 0.5, "CRVAL1": 53.1, "CRVAL2": -27.8,
                     "CD1_1": -0.06 / 3600, "CD1_2": 0.0, "CD2_1": 0.0, "CD2_2": 0.06 / 3600, "CUNIT1": "deg", "CUNIT2": "deg"})
    h.header.update(cards)
    h.writeto(path, overwrite=True)


@pytest.fixture(scope="module")
def img(tmp_path_factory):
    d = tmp_path_factory.mktemp("ds10img")
    rng = np.random.default_rng(5)
    a = rng.normal(0.0, 0.01, (200, 300)).astype(np.float32)
    yy, xx = np.mgrid[0:200, 0:300]
    for cx, cy, f in ((100.0, 80.0, 40.0), (200.5, 120.25, 15.0)):
        a += (f / (2 * math.pi * 2.0 ** 2) * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 2.0 ** 2))).astype(np.float32)
    p = d / "one_f160w.fits"
    write_fits(p, a, FILTER="F160W", BUNIT="ELECTRONS/S", PHOTFLAM=1.9275707e-20, PHOTPLAM=15369.161, PHOTZPT=-21.1, EXPTIME=1000.0, INSTRUME="WFC3")
    return p, a


@needs_core
def test_calib_header_zero_point_matches_the_formula(img):
    p, _ = img
    r = json.loads(run("calib", p).stdout)
    zp = -2.5 * math.log10(1.9275707e-20) + (-21.1) - 5 * math.log10(15369.161) + 18.6921
    assert r["effective"]["origin"] == "header" and r["effective"]["zeropoint"] == pytest.approx(zp, abs=1e-5)
    t = run("calib", p, "--text").stdout
    assert "AB zero point: %.4f" % zp in t and "F160W" in t
    o = json.loads(run("calib", p, "--zp", 30.5).stdout)
    assert o["effective"]["origin"] != "header" and o["effective"]["zeropoint"] == 30.5


@needs_core
def test_regions_stats_mask_and_convert_against_numpy(img, tmp_path):
    p, a = img
    reg = tmp_path / "r.reg"
    reg.write_text("# Region file format: DS9 version 4.1\nimage\ncircle(100,80,10.3)\nbox(200.3,120.2,20,12,0)\nannulus(100,80,15.3,25.3)\n")
    st = json.loads(run("regions", "stats", reg, "--image", p, "--background", "none").stdout)
    ys, xs = np.mgrid[1:201, 1:301]                                                 # FITS pixel numbers
    ref = {"circle": (xs - 100) ** 2 + (ys - 80) ** 2 <= 10.3 ** 2, "box": (abs(xs - 200.3) <= 10) & (abs(ys - 120.2) <= 6),
           "annulus": ((xs - 100) ** 2 + (ys - 80) ** 2 > 15.3 ** 2) & ((xs - 100) ** 2 + (ys - 80) ** 2 <= 25.3 ** 2)}
    # (regions whose boundary passes exactly through pixel centres can flip a few pixels when an image-frame region is stored in sky coordinates:
    #  docs/regions_and_pixel_table.md of Astrafex Web; the test therefore uses boundaries between pixel centres)
    for row in st["stats"]:
        m = ref[row["shape"]]
        assert row["n_pix"] == int(m.sum()), (row["shape"], row["n_pix"], int(m.sum()))
        assert row["sum"] == pytest.approx(float(a[m].sum()), rel=1e-5, abs=1e-4)
    out = tmp_path / "mask.fits"
    run("regions", "mask", reg, "--image", p, "-o", out)
    from astropy.io import fits
    mk = fits.getdata(out)
    assert mk.shape == a.shape and set(np.unique(mk)) <= {0, 1}
    union = ref["circle"] | ref["box"] | ref["annulus"]
    assert (mk.astype(bool) != union).sum() == 0
    # DS9 -> JSON -> DS9 keeps the shapes
    js = tmp_path / "r.json"
    run("regions", "convert", reg, "--image", p, "--to", "json", "-o", js)
    back = run("regions", "convert", js, "--image", p, "--to", "ds9").stdout
    assert back.count("circle(") == 1 and back.count("box(") == 1 and back.count("annulus(") == 1


@needs_core
def test_pixel_table_csv_equals_the_pixels(img, tmp_path):
    p, a = img
    out = tmp_path / "pt.csv"
    run("pixtab", p, "--xy", 100, 80, "--size", 5, "--csv", out)
    vals = {(int(r["x"]), int(r["y"])): float(r["value"]) for r in csv.DictReader(open(out))}
    assert len(vals) == 25
    for (x, y), v in vals.items():
        assert v == pytest.approx(float(a[y - 1, x - 1]), rel=1e-6, abs=1e-7)


# ------------------------------------------------------------------------------------------------ multi-band chain
@pytest.fixture(scope="module")
def bands(tmp_path_factory):
    d = tmp_path_factory.mktemp("ds10bands")
    rng = np.random.default_rng(11)
    n = 220
    yy, xx = np.mgrid[0:n, 0:n]
    src = [(40 + 25 * i + 3 * (i % 2), 40 + 30 * (i % 5) + 5 * (i // 5), 8.0 + 4.0 * i) for i in range(7)]
    spec = {"F105W": (1.0, 26.2687, 1.0551e-19, 10551.0), "F125W": (1.4, 26.2303, 2.2486e-20, 12486.0), "F160W": (2.0, 25.9463, 1.9276e-20, 15369.0)}
    paths = []
    for lab, (scale, zp, flam, plam) in spec.items():
        a = rng.normal(0.0, 0.02, (n, n))
        for cx, cy, f in src:
            a += scale * f / (2 * math.pi * 1.8 ** 2) * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.8 ** 2)) * 0.05
        p = d / f"syn_{lab.lower()}.fits"
        write_fits(p, a, FILTER=lab, BUNIT="ELECTRONS/S", PHOTFLAM=flam, PHOTPLAM=plam, PHOTZPT=-21.1, EXPTIME=1000.0, INSTRUME="WFC3")
        paths.append(p)
    return paths, src


@needs_core
def test_forced_chain_runs_and_script_run_reproduces_it_exactly(bands, tmp_path):
    paths, src = bands
    out = tmp_path / "out"
    r = run("forced", *paths, "--work", tmp_path / "w1", "--out", out, "--thresh", 3.0, "--smooth-fwhm", 2.0, "--minarea", 4, "--aperture-radii", "0.2,0.4", "--cog-radius", 0.5)
    res = json.loads(r.stdout)
    assert res["ok"] and {"forced_summary.json", "forced_catalog.tsv", "detect_catalog.tsv", "script.json", "detection.fits", "detection_rms.fits"} <= set(res["files"])
    summ = json.loads((out / "forced_summary.json").read_text())
    assert [b["label"] for b in summ["bands"]] == ["F105W", "F125W", "F160W"] and all(b["zeropoint_origin"] == "header" for b in summ["bands"])
    assert summ["n_sources"] >= 5
    header = open(out / "forced_catalog.tsv").readline().rstrip("\n").split("\t")
    assert any(c.startswith("MAG_AUTO_F105W") for c in header) and any(c.startswith("COLOR_F105W_F160W") for c in header)
    # the same step list run through `astrafex run` gives byte-identical catalogues
    run("run", out / "script.json", "--files", *paths, "--work", tmp_path / "w2", "--out", tmp_path / "out2")
    h = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
    assert h(tmp_path / "w1" / "steps" / "s3" / "work" / "catalog.tsv") == h(tmp_path / "w2" / "steps" / "s3" / "work" / "catalog.tsv")
    j = lambda q: {k: v for k, v in json.loads(q.read_text()).items() if k != "seconds"}          # only the wall-clock field differs
    assert j(out / "forced_summary.json") == j(tmp_path / "w2" / "steps" / "s3" / "work" / "forced_summary.json")


@needs_core
def test_run_validates_before_running_and_reports_bad_scripts(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": "astrafex-script/1", "steps": [{"id": "s1", "plugin": "nosuch", "step": "x"}]}))
    r = run("run", bad, "--work", tmp_path / "w", check=False)
    assert r.returncode == 2 and "not valid" in r.stderr
    drop = tmp_path / "drop.json"
    drop.write_text(json.dumps({"schema": "astrafex-script/1", "steps": [{"id": "s1", "ui": "delete_file"}]}))
    assert run("run", drop, "--work", tmp_path / "w3", check=False).returncode == 2


# ------------------------------------------------------------------------------------------------ offline bundle of the web app (astrafex open) + paths with spaces
BUNDLE = os.path.join(HERE, "data", "offline_bundle_small.zip")           # made by Astrafex Web tools/make_offline_fixture.py: one 96x96 image, one region, one detect job


def test_split_words_separators_and_spaces_in_names():
    assert launcher.split_words("a.fits b.fits c.fits") == ["a.fits", "b.fits", "c.fits"]
    assert launcher.split_words("/my data/a.fits|/my data/b.fits") == ["/my data/a.fits", "/my data/b.fits"]
    assert launcher.split_words("/my data/a.fits\n /my data/b.fits \n") == ["/my data/a.fits", "/my data/b.fits"]
    assert launcher.split_words('"/my data/a.fits" /b.fits') == ["/my data/a.fits", "/b.fits"]
    assert launcher.split_words("/o'brien/a.fits /b.fits") == ["/o'brien/a.fits", "/b.fits"]          # unbalanced quote: plain spaces


def test_open_without_dest_goes_to_home_astrafex_offline():
    home = os.path.expanduser("~")
    w = launcher.default_open_dest(["open", "/x/y/session one.zip", "--force", "--text"])
    assert w[-2:] == ["--dest", os.path.join(home, "astrafex-offline", "session one")]
    assert launcher.default_open_dest(["open", "b.zip", "--dest", "", "--force"]) == ["open", "b.zip", "--force", "--dest", os.path.join(home, "astrafex-offline", "b")]
    assert launcher.default_open_dest(["open", "b.zip", "--dest", "--force"]) == ["open", "b.zip", "--force", "--dest", os.path.join(home, "astrafex-offline", "b")]   # empty value dropped by the template
    assert launcher.default_open_dest(["open", "b.zip", "--dest", "/w"]) == ["open", "b.zip", "--dest", "/w"]
    assert launcher.default_open_dest(["info", "i.fits"]) == ["info", "i.fits"]


def test_manifest_has_open_bundle_step_with_a_tcl_hook():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json")))
    st = {s["id"]: s for s in m["steps"]}["open-bundle"]
    assert st["after"] == "OGFDs10coreOpenAfter" and m["tcl"] == "ds10core.tcl" and "open" in st["cli"] and "--force" in st["cli"]
    assert "proc OGFDs10coreOpenAfter" in open(os.path.join(PLUGIN, "ds10core.tcl")).read()
    assert {p["name"] for p in m["params"]} >= {"bundle", "open-dest"}


@needs_core
def test_open_restores_the_web_session_and_continues(tmp_path):
    ws = tmp_path / "my workspace"                                                      # a space in the folder name on purpose
    r = run("open", BUNDLE, "--dest", ws, "--text")
    assert "field96.fits" in r.stdout and "regions:" in r.stdout and "steps with stored results: s1" in r.stdout
    assert (ws / "files" / "field96.fits").is_file() and (ws / "regions" / "regions.reg").is_file() and (ws / "session_script.json").is_file()
    wsj = json.loads((ws / "workspace.json").read_text())
    assert wsj["schema"] == "astrafex-workspace/1" and wsj["files"][0]["kind"] == "original" and wsj["steps_with_results"] == ["s1"]
    again = run("open", BUNDLE, "--dest", ws, check=False)
    assert again.returncode != 0 and "--force" in again.stderr                          # never overwrites silently
    run("open", BUNDLE, "--dest", ws, "--force", "--text")
    # the stored result is used (status restored), nothing is recomputed
    run("run", ws / "session_script.json", "--text")
    rep = json.loads(next((ws / "runs").glob("r*/out/run_report.json")).read_text())
    assert rep["ok"] and [s["status"] for s in rep["steps"]] == ["restored"]
    st = run("status", ws, "--text")
    assert "s1:restored" in st.stdout


@needs_core
def test_open_default_destination_under_home(tmp_path):
    r = run("open", BUNDLE, "--text", env={"HOME": str(tmp_path)})
    assert (tmp_path / "astrafex-offline" / "offline_bundle_small" / "workspace.json").is_file() and str(tmp_path / "astrafex-offline") in r.stdout


@needs_core
def test_forced_with_spaces_in_the_band_image_paths(bands, tmp_path):
    import shutil
    paths, _ = bands
    d = tmp_path / "my band images"
    d.mkdir()
    new = [str(shutil.copy(p, d / os.path.basename(p))) for p in paths]
    r = run("forced", "split:" + "|".join(new), "--work", tmp_path / "w", "--out", tmp_path / "o", "--thresh", 3.0, "--smooth-fwhm", 2.0, "--minarea", 4, "--text")
    assert "forced photometry:" in r.stdout and (tmp_path / "o" / "forced_catalog.tsv").is_file()


def test_default_pack_drops_an_empty_output_name():
    assert launcher.default_pack(["pack", "/w", "-o", "--text"]) == ["pack", "/w", "--text"]
    assert launcher.default_pack(["pack", "/w", "-o", "", "--text"]) == ["pack", "/w", "--text"]
    assert launcher.default_pack(["pack", "/w", "-o", "/x/r.zip", "--text"]) == ["pack", "/w", "-o", "/x/r.zip", "--text"]
    assert launcher.default_pack(["info", "i.fits"]) == ["info", "i.fits"]


def test_manifest_has_pack_return_step():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json")))
    st = {s["id"]: s for s in m["steps"]}["pack-return"]
    assert "pack" in st["cli"] and "{pack-workspace}" in st["cli"] and {p["name"] for p in m["params"]} >= {"pack-workspace", "pack-out"}


@needs_core
def test_pack_writes_the_return_bundle_of_the_workspace(tmp_path):
    ws = tmp_path / "my workspace"
    run("open", BUNDLE, "--dest", ws, "--text")
    script = json.loads((ws / "session_script.json").read_text())
    script["steps"].append(dict(script["steps"][0], id="s2", params=dict(script["steps"][0]["params"], **{"detect-thresh": 4.0})))
    (ws / "session_script.json").write_text(json.dumps(script))
    run("run", ws / "session_script.json", "--text")
    out = tmp_path / "back" / "my return.zip"
    r = run("pack", ws, "-o", out, "--text")
    assert "s1: web" in r.stdout and "s2: local" in r.stdout and out.is_file()
    import zipfile
    names = zipfile.ZipFile(out).namelist()
    assert "return.json" in names and "results/s2/catalog.tsv" in names and not any(n.startswith("results/s1/") for n in names)
    d = json.loads(run("pack", ws).stdout)                                          # default name: in the workspace folder
    assert d["path"].startswith(str(ws))


# ------------------------------------------------------------------------------------------------ Astrafex names + backward compatibility
def test_user_visible_labels_say_astrafex_core_and_ids_are_kept():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json")))
    assert m["short"] == "Astrafex Core" and m["name"].startswith("Astrafex Core")
    assert m["id"] == "ds10core"                                                   # internal: parameter file, session records and the web manifests depend on it
    for s in m["steps"]:
        assert "ds10" not in s.get("label", "") + s.get("title", "") + s.get("done_status", ""), s["id"]
        assert s["record"].startswith("ds10core.")


def test_scripts_astrafex_and_the_former_ds10_name_both_work():
    outs = []
    for name in ("astrafex", "ds10"):
        r = subprocess.run([os.path.join(ROOT, "scripts", name), "--where"], capture_output=True, text=True)
        assert r.returncode in (0, 3), (name, r.stderr)
        outs.append(r.stdout)
    assert outs[0] == outs[1]


@needs_core
def test_astrafex_env_aliases_and_old_ones(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DS10_", "ASTRAFEX_"))}
    for var in ("ASTRAFEX_CORE", "DS10_CORE"):
        r = subprocess.run([sys.executable, os.path.join(PLUGIN, "ds10.py"), "--where"], capture_output=True, text=True, env=dict(env, **{var: CORE}))
        d = json.loads(r.stdout)
        assert d["core_dir"] == CORE and d["found_by"] in ("DS10_CORE", "ASTRAFEX_CORE")


@needs_core
def test_the_fixture_is_a_real_pre_rename_bundle_and_still_opens(tmp_path):
    import zipfile
    with zipfile.ZipFile(BUNDLE) as z:
        assert json.loads(z.read("offline/offline.json"))["schema"] == "ds10-offline/1"
        assert json.loads(z.read("offline/session_script.json"))["schema"] == "ds10-script/1"
    ws = tmp_path / "w"
    r = run("open", BUNDLE, "--dest", ws, "--text", env={"HOME": str(tmp_path)})
    assert (ws / "workspace.json").is_file()
    out = run("run", ws / "session_script.json", env={"HOME": str(tmp_path)})
    assert json.loads(out.stdout)["ok"]


def test_manifest_exposes_psf_match_params():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json")))
    names = {q["name"] for q in m["params"]}
    assert {"psf-match", "psf-target", "psf-method", "psf-size"} <= names
    forced = next(s for s in m["steps"] if s["id"] == "forced")
    assert any(isinstance(x, dict) and x.get("if") == "psf-match" for x in forced["cli"])
