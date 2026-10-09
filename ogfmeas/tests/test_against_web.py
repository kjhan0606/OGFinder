"""Compare ogfmeas with the web core on shared synthetic scenes.

The web process imports only ds10core. This process imports only ogfmeas.
Coordinates in the report are 0-based pixel centres.

Gates, fixed before looking at the run:
- pure-sky background within 0.15 of the planted level
- each isolated source matched within 0.6 px
- isophotal flux of those matches within 25 percent
- circular aperture: each code within 5 percent of the analytic Gaussian
  integral, and the two codes within 3 percent of that integral
- two blended peaks: each code places a source within 1.2 px of both centres
"""
import json
import math
import os
import subprocess
import sys
import textwrap

import numpy as np

from ogfmeas import Background, extract, flux_radius, sum_circle

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
WEB_SERVER = os.path.join(ROOT, "ds10-web", "server")

WEB = textwrap.dedent(
    r'''
    import json, sys
    import numpy as np
    from ds10core.tools.detect import detect
    from ds10core.aperphot import aperture_flux
    from ds10core.morph import _one

    def quiet(*_a, **_k):
        return None

    z = np.load(sys.argv[1])

    def pack(rows):
        return [{"x": float(r["X_IMAGE"]) - 1.0, "y": float(r["Y_IMAGE"]) - 1.0,
                 "flux": float(r["FLUX_ISO"]), "a": float(r["A_IMAGE"]), "b": float(r["B_IMAGE"]),
                 "flags": int(r["FLAGS"]), "snr": float(r["SNR"])} for r in rows]

    def run_detect(img, deblend):
        rows, _seg, info = detect(
            img, thresh=3.0, minarea=5, smooth=False, back_size=32,
            noise_correction=False, deblend=deblend, deblend_nlevels=32,
            deblend_mincont=0.01, report=quiet, max_sources=500)
        return {"back": float(info["background_median"]), "rms": float(info["rms_median"]),
                "sources": pack(rows)}

    gauss = z["gauss"]
    aper = aperture_flux(gauss, 32.0, 32.0, 8.0, 0.0, 0.0, 100)
    morph = _one(gauss, 32.0, 32.0, None, None, None, [], sky=0.0, rms=0.01,
                 eta=0.2, r_max=20.0, kron_scale=2.5, kron_min=1.0,
                 neighbour_radius=3.0, zeropoint=25.0, mask_neighbours=False,
                 circular=True, seed=1)
    out = {
        "sep_loaded": ("sep" in sys.modules) or ("sep_pjw" in sys.modules),
        "sky": run_detect(z["sky"], False),
        "isolated": run_detect(z["isolated"], False),
        "blended": run_detect(z["blended"], True),
        "aperture": None if aper is None else float(aper["flux"]),
        "morph_r50": morph.get("MORPH_R50"),
        "morph_kron": morph.get("MORPH_KRON_R"),
    }
    json.dump(out, sys.stdout)
    '''
)


def _gauss(shape, x0, y0, amp, sig):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amp * np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / sig ** 2)


def _scene():
    rng = np.random.default_rng(7)
    sky = 8.0 + rng.normal(0.0, 0.4, size=(96, 96))
    isolated = 5.0 + rng.normal(0.0, 0.35, size=(96, 96))
    planted = [(30.0, 40.0, 90.0, 1.8), (62.0, 28.0, 70.0, 2.0), (70.0, 68.0, 55.0, 1.6)]
    for x0, y0, amp, sig in planted:
        isolated += _gauss(isolated.shape, x0, y0, amp, sig)
    blended = 1.0 + rng.normal(0.0, 0.25, size=(80, 80))
    blend_at = [(34.0, 40.0, 220.0, 2.2), (44.0, 40.0, 180.0, 2.2)]
    for x0, y0, amp, sig in blend_at:
        blended += _gauss(blended.shape, x0, y0, amp, sig)
    gauss = _gauss((64, 64), 32.0, 32.0, 100.0, 2.0)
    return sky, isolated, blended, gauss, planted, blend_at


def _ogf_detect(image, deblend):
    work = np.asarray(image, dtype=np.float64)
    bad = ~np.isfinite(work)
    filled = np.where(bad, 0.0, work)
    bkg = Background(filled, mask=bad, bw=32, bh=32)
    sub = filled - bkg.back()
    rms = np.maximum(bkg.rms(), 1e-8)
    cat = extract(
        sub, 3.0, err=rms, mask=bad, minarea=5, filter_kernel=None,
        deblend_nthresh=32 if deblend else 1, deblend_cont=0.01)
    sources = [{"x": float(o["x"]), "y": float(o["y"]), "flux": float(o["flux"]),
                "a": float(o["a"]), "b": float(o["b"]), "flags": int(o["flag"])} for o in cat]
    return {"back": float(bkg.globalback), "rms": float(bkg.globalrms), "sources": sources}


def _match(found, truth, radius):
    used = set()
    pairs = []
    for x0, y0, *_rest in truth:
        best, bestd = None, radius
        for i, src in enumerate(found):
            if i in used:
                continue
            dist = math.hypot(src["x"] - x0, src["y"] - y0)
            if dist < bestd:
                best, bestd = i, dist
        if best is None:
            pairs.append({"truth": (x0, y0), "found": None, "dist": None})
        else:
            used.add(best)
            pairs.append({"truth": (x0, y0), "found": found[best], "dist": bestd})
    return pairs


def _web(path):
    env = os.environ.copy()
    env["PYTHONPATH"] = WEB_SERVER + os.pathsep + env.get("PYTHONPATH", "")
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    proc = subprocess.run(
        [sys.executable, "-c", WEB, path],
        cwd=WEB_SERVER, env=env, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise AssertionError(proc.stderr[-2000:] or proc.stdout[-2000:])
    return json.loads(proc.stdout)


def test_web_and_ogfmeas_agree_on_synthetic_scenes(tmp_path):
    sky, isolated, blended, gauss, planted, blend_at = _scene()
    npz = tmp_path / "scenes.npz"
    np.savez(npz, sky=sky, isolated=isolated, blended=blended, gauss=gauss)
    web = _web(str(npz))
    ogf = {
        "sky": _ogf_detect(sky, False),
        "isolated": _ogf_detect(isolated, False),
        "blended": _ogf_detect(blended, True),
    }
    flux, _, _ = sum_circle(gauss, 32.0, 32.0, 8.0, subpix=5)
    half, _ = flux_radius(gauss, 32.0, 32.0, 15.0, 0.5, subpix=5)
    analytic = 100.0 * 2.0 * math.pi * 4.0
    half_true = 2.0 * math.sqrt(2.0 * math.log(2.0))
    report = {
        "sky": {"truth": 8.0, "web": web["sky"]["back"], "ogfmeas": ogf["sky"]["back"],
                "web_n": len(web["sky"]["sources"]), "ogf_n": len(ogf["sky"]["sources"])},
        "isolated": {
            "web": _match(web["isolated"]["sources"], planted, 2.0),
            "ogfmeas": _match(ogf["isolated"]["sources"], planted, 2.0),
            "web_n": len(web["isolated"]["sources"]),
            "ogf_n": len(ogf["isolated"]["sources"]),
        },
        "blended": {
            "web": _match(web["blended"]["sources"], blend_at, 2.0),
            "ogfmeas": _match(ogf["blended"]["sources"], blend_at, 2.0),
            "web_n": len(web["blended"]["sources"]),
            "ogf_n": len(ogf["blended"]["sources"]),
        },
        "aperture": {"analytic": analytic, "web": web["aperture"], "ogfmeas": float(flux[0])},
        "half_light": {"analytic": half_true, "ogfmeas_flux_radius": float(half[0]),
                       "web_morph_r50": web["morph_r50"], "web_morph_kron": web["morph_kron"]},
        "web_loaded_sep": web["sep_loaded"],
    }
    # Make the failure message the whole report.
    problems = []
    if web["sep_loaded"]:
        problems.append("web process loaded sep")
    for name in ("web", "ogfmeas"):
        if abs(report["sky"][name] - 8.0) > 0.15:
            problems.append(f"sky {name} background {report['sky'][name]}")
    for side in ("web", "ogfmeas"):
        for pair in report["isolated"][side]:
            if pair["found"] is None or pair["dist"] > 0.6:
                problems.append(f"isolated {side} {pair}")
    web_by_truth = {tuple(p["truth"]): p for p in report["isolated"]["web"]}
    ogf_by_truth = {tuple(p["truth"]): p for p in report["isolated"]["ogfmeas"]}
    for key in web_by_truth:
        a, b = web_by_truth[key]["found"], ogf_by_truth[key]["found"]
        if a is None or b is None:
            continue
        scale = max(abs(a["flux"]), abs(b["flux"]), 1.0)
        if abs(a["flux"] - b["flux"]) / scale > 0.25:
            problems.append(f"isophotal flux {key}: web {a['flux']:.4g} ogfmeas {b['flux']:.4g}")
    aper = report["aperture"]
    for name in ("web", "ogfmeas"):
        if abs(aper[name] - analytic) / analytic > 0.05:
            problems.append(f"aperture {name} {aper[name]:.4g} vs analytic {analytic:.4g}")
    if abs(aper["web"] - aper["ogfmeas"]) / analytic > 0.03:
        problems.append(f"aperture web-ogfmeas {aper['web']:.4g} vs {aper['ogfmeas']:.4g}")
    if abs(float(half[0]) - half_true) > 0.45:
        problems.append(f"ogfmeas half-light {float(half[0]):.3f} analytic {half_true:.3f}")
    for side in ("web", "ogfmeas"):
        for pair in report["blended"][side]:
            if pair["found"] is None or pair["dist"] > 1.2:
                problems.append(f"blended {side} {pair}")
    assert not problems, json.dumps({"problems": problems, "report": report}, indent=2, default=str)
