"""Measure the three example spectra and write science.json.

Each file uses its own redshift limit and detection threshold.  The desktop
menu applies one parameter set to the whole catalogue, so this script also
runs plugins/spectra/spectra.py --task science once per file.
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

from ogfkit import spectra as sp
from ogfkit import specscience as sci

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(_ROOT, "plugins", "spectra", "spectra.py")
PARAMS = {
    "1": dict(frame="air", ebv_mw=0.0, inst_fwhm_a=2.0, zmax=0.6, snr_min=4.0, h0=70.0, om0=0.3),
    "2": dict(frame="vacuum", ebv_mw=0.0, inst_fwhm_a=2.0, zmax=0.5, snr_min=5.0, h0=70.0, om0=0.3),
    "3": dict(frame="vacuum", ebv_mw=0.0, inst_fwhm_a=2.0, zmax=2.0, snr_min=4.0, h0=70.0, om0=0.3),
}


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


def _slim(an):
    z = an["redshift"]
    disp = an["dispersion"]
    sfr = an["sfr"]
    balmer = an["balmer"] or {}
    bpt = an["bpt"]
    uv = an["uv"]
    return dict(
        type=an["type"],
        redshift=dict(z=z.get("z"), z_err=z.get("z_err"), quality=z.get("quality"), source=z.get("source")),
        emission=dict(z=an["emission"]["z"], quality=an["emission"]["quality"], n_lines=an["emission"]["n_lines"]),
        absorption=dict(z=an["absorption"]["z"], quality=an["absorption"]["quality"], n_lines=an["absorption"]["n_lines"]),
        uv=None if not uv else dict(z=uv["z"], z_err=uv.get("z_err"), quality=uv["quality"], n_lines=uv["n_lines"]),
        dispersion=dict(sigma=disp.get("sigma"), sigma_err=disp.get("sigma_err"), kind=disp.get("kind"),
                        corrected=disp.get("corrected"), n=disp.get("n"), note=disp.get("note")),
        bpt=None if bpt is None else bpt.get("class"),
        broad=dict(flag=an["broad"]["flag"], line=an["broad"].get("line"), fwhm_kms=an["broad"].get("fwhm_kms"),
                   cut_kms=an["broad"].get("cut_kms"), reference=an["broad"].get("reference")),
        balmer=dict(used=bool(balmer.get("used")), ebv=balmer.get("ebv"), ratio=balmer.get("ratio")),
        ha_flux=an.get("ha_flux"),
        sfr=dict(sfr=sfr.get("sfr"), sfr_err=sfr.get("sfr_err"), applies=sfr.get("applies"),
                 flux_unit=sfr.get("flux_unit"), error_includes_balmer=sfr.get("error_includes_balmer"),
                 calibration=sfr.get("calibration")),
        preprocess=an["preprocess"],
        calibrations=an["calibrations"],
    )


def _cli(number, params):
    work = os.path.join(HERE, "work", number)
    os.makedirs(work, exist_ok=True)
    catalog = os.path.join(work, "catalog.tsv")
    with open(catalog, "w") as fh:
        fh.write("NUMBER\tX_IMAGE\tY_IMAGE\n%s\t1\t1\n" % number)
    cmd = [sys.executable, CLI, "--task", "science", "--catalog", catalog, "--work", work,
           "--spec-dir", os.path.join(HERE, "data"), "--spec-pattern", "spec_{NUMBER}.txt",
           "--wave-frame", params["frame"], "--ebv-mw", str(params["ebv_mw"]),
           "--inst-fwhm", str(params["inst_fwhm_a"]), "--zmax", str(params["zmax"]),
           "--snr-min", str(params["snr_min"]), "--h0", str(params["h0"]), "--omega-m", str(params["om0"])]
    proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(proc.stderr[-800:])
    lines = [ln for ln in proc.stdout.splitlines() if ln]
    cols = lines[0].split("\t")
    row = dict(zip(cols, lines[1].split("\t")))
    return dict(command=cmd[1:], columns=row, stderr=proc.stderr.strip())


def main():
    pre_path = os.path.join(HERE, "preprocess.json")
    if not os.path.isfile(pre_path):
        raise SystemExit("run 01_preprocess.py first")
    pre = json.load(open(pre_path))
    out = dict(objects={})
    for number, params in PARAMS.items():
        path = os.path.join(HERE, "data", "spec_%s.txt" % number)
        raw = sp.read_spectrum(path)
        an = sci.galaxy_analysis(raw["wave"], raw["flux"], raw["err"], **params)
        slim = _slim(an)
        cli = _cli(number, params)
        out["objects"][number] = dict(
            truth=dict(z=pre["objects"][number]["z_truth"], kind=pre["objects"][number]["kind"]),
            parameters=params,
            analysis=slim,
            cli=cli,
        )
        print(number, slim["type"], slim["redshift"]["source"], slim["redshift"]["z"],
              "cli", cli["columns"].get("SP_TYPE"), cli["columns"].get("SP_SCI_ZSRC"))
    with open(os.path.join(HERE, "science.json"), "w") as fh:
        json.dump(_clean(out), fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
