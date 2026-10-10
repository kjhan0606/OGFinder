"""Write synthetic supernova spectra and record sn_analysis plus --task sntype.

Five vacuum files, all at host redshift 0.05.  A sixth call repeats the
Type II spectrum with no host redshift.  Flux is the generator's file unit.
The desktop menu applies one parameter set to the catalogue; this script
does that, in one CLI call, because these files share a frame.
"""
import json
import os
import subprocess
import sys

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _ROOT)

from ogfkit import sntype as sn

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CLI = os.path.join(_ROOT, "plugins", "spectra", "spectra.py")
Z = 0.05

# name, kind, velocity, gaussian sigma, seed, absorption depths, emission (name, sigma_v, amp)
OBJECTS = [
    ("1", "Ia", 10000.0, 2500.0, 1,
     {"SiII6355": 0.55, "SII5454": 0.25, "SII5640": 0.25, "OI7774": 0.15}, None),
    ("2", "II", 8000.0, 2500.0, 2,
     {"Ha": 0.45, "Hb": 0.30}, None),
    ("3", "Ib", 9000.0, 2500.0, 3,
     {"HeI5876": 0.40, "HeI6678": 0.30, "HeI7065": 0.28}, None),
    ("4", "Ic-BL", 32000.0, 8000.0, 7,
     {"OI7774": 0.25, "SiII6355": 0.25}, None),
    ("5", "IIn", 0.0, 2500.0, 6,
     {}, [("Ha", 400.0, 2.0)]),
]


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _spectrum(depths, z, v, sigma_v, seed, emission):
    w = np.arange(4500.0, 9200.0, 2.0)
    f = np.ones(w.shape, float)
    for name, depth in depths.items():
        lam = sn.FEAT[name][0]
        cen = lam * (1.0 + z) * (1.0 - v / sn.C_KMS)
        sig = max(cen * sigma_v / sn.C_KMS, 1.0)
        f = f - depth * np.exp(-0.5 * ((w - cen) / sig) ** 2)
    if emission:
        for name, sig_v, amp in emission:
            lam = sn.FEAT[name][0]
            cen = lam * (1.0 + z)
            sig = max(cen * sig_v / sn.C_KMS, 1.0)
            f = f + amp * np.exp(-0.5 * ((w - cen) / sig) ** 2)
    rng = np.random.default_rng(seed)
    f = f + rng.normal(0.0, 0.015, size=w.shape)
    err = np.full(w.shape, 0.015)
    return w, f, err


def _write(path, wave, flux, err, header):
    with open(path, "w") as fh:
        fh.write(header)
        for w, f, e in zip(wave, flux, err):
            fh.write("%.8e %.8e %.8e\n" % (w, f, e))


def _slim(an):
    narrow = an.get("narrow") or {}
    ha = narrow.get("Ha")
    return dict(
        type=an["type"], z=an["z"], quality=an["quality"], v=an["v"], v_line=an["v_line"],
        a_si=an["a_si"], a_ha=an["a_ha"], a_he=an["a_he"], a_oi=an["a_oi"],
        hydrogen=an["hydrogen"], helium_n=an["helium_n"], v_grid=an["v_grid"],
        narrow_ha=None if not ha else dict(ok=ha["ok"], fwhm=ha["fwhm"], dv=ha["dv"], snr=ha["snr"]),
        note=an["note"],
    )


def _cli(catalog_text, work):
    os.makedirs(work, exist_ok=True)
    catalog = os.path.join(work, "catalog.tsv")
    with open(catalog, "w") as fh:
        fh.write(catalog_text)
    cmd = [sys.executable, CLI, "--task", "sntype", "--catalog", catalog, "--work", work,
           "--spec-dir", DATA, "--spec-pattern", "sn_{NUMBER}.txt",
           "--z-column", "Z", "--wave-frame", "vacuum", "--ebv-mw", "0"]
    proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(proc.stderr[-800:])
    lines = [ln for ln in proc.stdout.splitlines() if ln]
    cols = lines[0].split("\t")
    rows = [dict(zip(cols, ln.split("\t"))) for ln in lines[1:]]
    return dict(command=cmd[1:], columns=rows, stderr=proc.stderr.strip())


def main():
    os.makedirs(DATA, exist_ok=True)
    arrays = {}
    out = dict(calibration=sn.CALIBRATION, objects={})
    for number, kind, vel, sigma_v, seed, depths, emission in OBJECTS:
        w, f, e = _spectrum(depths, Z, vel, sigma_v, seed, emission)
        arrays[number] = (w, f, e)
        path = os.path.join(DATA, "sn_%s.txt" % number)
        _write(path, w, f, e, "# synthetic supernova %s host_z=0.05 vacuum file-unit flux\n" % kind)
        an = sn.sn_analysis(w, f, e, z_host=Z, frame="vacuum")
        out["objects"][number] = dict(
            truth=dict(kind=kind, z=Z, v=vel, sigma_v=sigma_v, seed=seed, frame="vacuum"),
            analysis=_slim(an),
        )
        print(number, kind, an["type"], an["v"], an["v_line"], an["a_si"], an["a_ha"])
    w2, f2, e2 = arrays["2"]
    missing = sn.sn_analysis(w2, f2, e2, z_host=None, frame="vacuum")
    catalog = "NUMBER\tZ\n" + "".join("%s\t0.05\n" % n for n, *_rest in OBJECTS)
    cli = _cli(catalog, os.path.join(HERE, "work", "sntype"))
    by_num = {row["NUMBER"]: row for row in cli["columns"]}
    for number in out["objects"]:
        out["objects"][number]["cli"] = by_num[number]
    miss_cli = _cli("NUMBER\tZ\n2\t\n", os.path.join(HERE, "work", "sntype_noz"))
    out["missing_host_redshift"] = dict(
        file="sn_2.txt",
        analysis=_slim(missing),
        cli=miss_cli["columns"][0],
        stderr=miss_cli["stderr"],
    )
    out["cli_stderr"] = cli["stderr"]
    with open(os.path.join(HERE, "sntype.json"), "w") as fh:
        json.dump(_clean(out), fh, indent=2)
        fh.write("\n")
    print("missing", missing["type"], missing["quality"])


if __name__ == "__main__":
    main()
