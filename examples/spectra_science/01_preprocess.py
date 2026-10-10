"""Write three synthetic 1D spectra and record prepare_spectrum.

Object 1 is a star-forming galaxy at z = 0.2, stored in air wavelengths,
with the rows reversed and two bad samples appended.  Object 2 is a star
at z = 0.0002 plus one emission line at 4500 Angstrom.  Object 3 is a
broad-line spectrum at z = 1.6 (C IV, C III], Mg II).  Flux is the file
unit of the generator.  It is not scaled to erg/s/cm^2.
"""
import json
import math
import os
import sys

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _ROOT)

from ogfkit import spectra as sp
from ogfkit import specscience as sci
from ogfkit import spectrasynth as ss

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FWHM = sp.FWHM


def _write_columns(path, wave, flux, err, header):
    with open(path, "w") as fh:
        fh.write(header)
        for w, f, e in zip(wave, flux, err):
            fh.write("%.8e %.8e %.8e\n" % (w, f, e))


def _star_forming():
    lines = {"Ha": 100.0, "Hb": 100.0 / 2.86, "[OIII]5008": 10.0, "[NII]6585": 10.0}
    w, f, e, _t = ss.spectrum_1d(z=0.2, seed=4, wave=(5400.0, 8400.0), dlam=1.0, cont=1.0, slope=0.0,
                                 noise=0.03, lines=lines, fwhm_kms=200.0, res_fwhm_A=2.0)
    air = np.asarray(sci.vacuum_to_air(w), float)
    return air[::-1].copy(), f[::-1].copy(), e[::-1].copy()


def _star():
    rng = np.random.default_rng(3)
    w = np.arange(3800.0, 9200.0, 1.0)
    z = 0.0002
    cont = 40.0
    f = np.full_like(w, cont)
    for _name, lam0, _wt in sci.ABS_LINES:
        lam = lam0 * (1.0 + z)
        if w[0] + 8 < lam < w[-1] - 8:
            sig = math.hypot(180.0 / sci.C_KMS * lam, 2.0) / FWHM
            f = f - cont * 0.55 * np.exp(-0.5 * ((w - lam) / sig) ** 2)
    err = np.full_like(w, 0.04)
    f = f + rng.normal(0.0, 0.04, len(w))
    lam = 4500.0
    sig = 2.0
    f = f + 30.0 / (sig * math.sqrt(2.0 * math.pi)) * np.exp(-0.5 * ((w - lam) / sig) ** 2)
    return w, f, err


def _qso():
    lines = {"CIV": 80.0, "CIII]": 40.0, "MgII": 60.0}
    w, f, e, _t = ss.spectrum_1d(z=1.6, seed=9, wave=(3700.0, 9000.0), dlam=1.0, cont=1.0, slope=0.0,
                                 noise=0.05, lines=lines, fwhm_kms=4000.0, res_fwhm_A=2.0)
    return w, f, e


def _prep_record(path, frame, ebv_mw):
    raw = sp.read_spectrum(path)
    prep = sci.prepare_spectrum(raw["wave"], raw["flux"], raw["err"], frame=frame, ebv_mw=ebv_mw)
    return dict(n_file_rows=int(len(raw["wave"])), n_dropped=prep["n_dropped"],
                n_kept=int(len(prep["wave"])), n_no_extinction=prep["n_no_extinction"],
                input_frame=prep["input_frame"], frame=prep["frame"], ebv_mw=prep["ebv_mw"],
                file_wave_first=float(raw["wave"][0]), file_wave_last=float(raw["wave"][-1]),
                wave_first=float(prep["wave"][0]), wave_last=float(prep["wave"][-1]))


def main():
    os.makedirs(DATA, exist_ok=True)
    sf_w, sf_f, sf_e = _star_forming()
    # Two numeric rows that prepare_spectrum must drop: a non-finite flux, and err <= 0.
    sf_w = np.concatenate([sf_w, [np.nan, 6100.0]])
    sf_f = np.concatenate([sf_f, [1.0, 1.0]])
    sf_e = np.concatenate([sf_e, [0.03, 0.0]])
    _write_columns(os.path.join(DATA, "spec_1.txt"), sf_w, sf_f, sf_e,
                   "# object 1: star-forming, air Angstrom, rows reversed; last two rows are bad\n"
                   "# wavelength flux err\n")
    st_w, st_f, st_e = _star()
    _write_columns(os.path.join(DATA, "spec_2.txt"), st_w, st_f, st_e,
                   "# object 2: absorption star plus one emission line at 4500 Angstrom, vacuum\n"
                   "# wavelength flux err\n")
    q_w, q_f, q_e = _qso()
    _write_columns(os.path.join(DATA, "spec_3.txt"), q_w, q_f, q_e,
                   "# object 3: CIV, C III] and Mg II at z=1.6, vacuum Angstrom\n"
                   "# wavelength flux err\n")
    with open(os.path.join(HERE, "catalog.tsv"), "w") as fh:
        fh.write("NUMBER\tX_IMAGE\tY_IMAGE\n1\t1\t1\n2\t1\t1\n3\t1\t1\n")

    air_ha = 6562.801
    vac = sci.air_to_vacuum(air_ha)
    ccm_w = np.array([5000.0, 6000.0, 7000.0])
    ccm_f = np.array([2.0, 3.0, 4.0])
    ccm_e = np.array([0.2, 0.2, 0.2])
    ebv = 0.1
    ccm = sci.prepare_spectrum(ccm_w, ccm_f, ccm_e, frame="vacuum", ebv_mw=ebv)
    report = dict(
        air_ha_input=air_ha,
        air_ha_vacuum=vac,
        air_ha_back=sci.vacuum_to_air(vac),
        ccm=dict(ebv_mw=ebv, wave=ccm_w.tolist(), flux_in=ccm_f.tolist(), err_in=ccm_e.tolist(),
                 flux_out=[float(v) for v in ccm["flux"]], err_out=[float(v) for v in ccm["err"]],
                 n_dropped=ccm["n_dropped"], n_no_extinction=ccm["n_no_extinction"]),
        objects={
            "1": dict(file="data/spec_1.txt", frame="air", ebv_mw=0.0, z_truth=0.2, kind="star-forming",
                      seed=4, lines={"Ha": 100.0, "Hb": 100.0 / 2.86, "[OIII]5008": 10.0, "[NII]6585": 10.0},
                      fwhm_kms=200.0, prepare=_prep_record(os.path.join(DATA, "spec_1.txt"), "air", 0.0)),
            "2": dict(file="data/spec_2.txt", frame="vacuum", ebv_mw=0.0, z_truth=0.0002, kind="star",
                      seed=3, false_emission_a=4500.0,
                      prepare=_prep_record(os.path.join(DATA, "spec_2.txt"), "vacuum", 0.0)),
            "3": dict(file="data/spec_3.txt", frame="vacuum", ebv_mw=0.0, z_truth=1.6, kind="broad-line",
                      seed=9, lines={"CIV": 80.0, "CIII]": 40.0, "MgII": 60.0}, fwhm_kms=4000.0,
                      prepare=_prep_record(os.path.join(DATA, "spec_3.txt"), "vacuum", 0.0)),
        },
    )
    with open(os.path.join(HERE, "preprocess.json"), "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    p1 = report["objects"]["1"]["prepare"]
    print("air 6562.801 -> %.6f" % vac)
    print("object 1 dropped %d kept %d first %.3f last %.3f" % (
        p1["n_dropped"], p1["n_kept"], p1["wave_first"], p1["wave_last"]))


if __name__ == "__main__":
    main()
