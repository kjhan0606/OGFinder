"""Shared pieces of the SED-code adapters: photometry conversion, record -> photometry, p(z) summary, result schema, subprocess helper.

Everything is a pure function of JSON-serialisable input (records are the ai_bridge input records), so the adapters can run in a service.
"""
import json
import math
import os
import subprocess
import sys

import numpy as np

from . import filters

AB_ZP_UJY = 23.9                      # AB magnitude of 1 microJansky
PHOTOZ_COLS = ['PHOTOZ', 'PHOTOZ_ERR', 'PHOTOZ_P16', 'PHOTOZ_P50', 'PHOTOZ_P84', 'PHOTOZ_CHI2']
SED_COLS = ['LOG_MASS', 'LOG_MASS_ERR', 'LOG_AGE', 'LOG_AGE_ERR', 'LOG_Z', 'AV', 'SFR', 'SED_CHI2']
TASK_COLS = {'photoz': PHOTOZ_COLS, 'sed_fit': SED_COLS}


def mag_to_ujy(m):
    return 10 ** (-0.4 * (np.asarray(m, float) - AB_ZP_UJY))


def ujy_to_mag(f):
    f = np.asarray(f, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(f > 0, AB_ZP_UJY - 2.5 * np.log10(f), np.nan)


def magerr_to_ujy_err(m, e):
    """Flux error (microJy) of a magnitude m +- e (first order: sigma_f = f sigma_m ln(10)/2.5)."""
    return mag_to_ujy(m) * np.asarray(e, float) * math.log(10) / 2.5


def photometry(record, min_mag_err=0.02, default_mag_err=0.1, bands=None):
    """Record -> (list of band entries, warnings).  Each entry: dict(band, entry (filters registry), mag, magerr, flux (uJy), fluxerr (uJy)).
    The error is combined in quadrature with `min_mag_err` (calibration floor); a missing error becomes `default_mag_err`.  Bands the registry does not know,
    and bands without a measurement, are left out (the warnings list says which)."""
    out, warn = [], []
    for b, m in (record.get('mags') or {}).items():
        if bands and b not in bands:
            continue
        ent = filters.lookup(b)
        if ent is None:
            if b.upper() not in ('AUTO', 'ISO', 'BEST'):
                warn.append('band %s unknown (see sed_adapters.filters)' % b)
            continue
        if m is None:
            continue
        e = (record.get('mag_errs') or {}).get(b)
        e = default_mag_err if e is None else e
        e = math.hypot(e, min_mag_err)
        out.append(dict(band=ent['name'], entry=ent, mag=float(m), magerr=float(e), flux=float(mag_to_ujy(m)), fluxerr=float(magerr_to_ujy_err(m, e))))
    out.sort(key=lambda d: d['entry']['lam'])
    return out, warn


def pdf_summary(z, pz):
    """Percentiles / mean / sigma of a p(z) sampled on grid z (unnormalised ok)."""
    z = np.asarray(z, float)
    p = np.clip(np.asarray(pz, float), 0, None)
    s = np.trapezoid(p, z) if hasattr(np, 'trapezoid') else np.trapz(p, z)
    if not (s > 0):
        return dict(z_mean=float('nan'), z_sigma=float('nan'), p16=float('nan'), p50=float('nan'), p84=float('nan'), z_peak=float('nan'))
    p = p / s
    cdf = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(z))])
    cdf /= cdf[-1]
    q = lambda f: float(np.interp(f, cdf, z))
    mean = float(np.trapezoid(z * p, z) if hasattr(np, 'trapezoid') else np.trapz(z * p, z))
    var = float((np.trapezoid((z - mean) ** 2 * p, z)) if hasattr(np, 'trapezoid') else np.trapz((z - mean) ** 2 * p, z))
    return dict(z_mean=mean, z_sigma=math.sqrt(max(var, 0)), p16=q(0.1587), p50=q(0.5), p84=q(0.8413), z_peak=float(z[int(np.argmax(p))]))


def photoz_row(zbest, p16, p50, p84, chi2=None, sigma=None):
    sig = sigma if sigma is not None else (0.5 * (p84 - p16) if (p84 is not None and p16 is not None) else None)
    return {'PHOTOZ': zbest, 'PHOTOZ_ERR': sig, 'PHOTOZ_P16': p16, 'PHOTOZ_P50': p50, 'PHOTOZ_P84': p84, 'PHOTOZ_CHI2': chi2}


def clean(d):
    """JSON-safe copy: NaN/inf -> None, numpy scalars -> python."""
    if isinstance(d, dict):
        return {k: clean(v) for k, v in d.items()}
    if isinstance(d, (list, tuple)):
        return [clean(v) for v in d]
    if isinstance(d, (np.floating, float)):
        v = float(d)
        return v if math.isfinite(v) else None
    if isinstance(d, np.integer):
        return int(d)
    return d


def run_command(argv, cwd=None, stdin_text=None, timeout=600, env_extra=None):
    """Run an argv list without a shell, minimal environment (+ env_extra); returns (rc, stdout, stderr)."""
    env = {k: os.environ[k] for k in ('PATH', 'HOME', 'LANG', 'TMPDIR', 'VIRTUAL_ENV', 'PYTHONPATH', 'EAZYCODE', 'SPS_HOME', 'OMP_NUM_THREADS') if k in os.environ}
    env.update(env_extra or {})
    p = subprocess.run(argv, input=stdin_text, capture_output=True, text=True, cwd=cwd, timeout=timeout, env=env)
    return p.returncode, p.stdout, p.stderr


def adapter_main(process, argv=None):
    """local_command entry: stdin = bridge JSON {task, params, records}, stdout = {"results": [...], "model": ...}."""
    req = json.load(sys.stdin)
    task = req.get('task') or 'photoz'
    res, model = process(req.get('records') or [], req.get('params') or {}, task)
    json.dump(clean({'results': res, 'model': model}), sys.stdout)
    return 0


def read_filter_res(path):
    """Parse an EAZY FILTER.RES file -> list of (name, wave_A array, throughput array)."""
    out, cur = [], None
    with open(path) as f:
        for ln in f:
            t = ln.split()
            if not t:
                continue
            if len(t) >= 2 and t[0].isdigit() and not (len(t) >= 3 and t[1].replace('.', '').replace('e', '').replace('-', '').replace('+', '').isdigit() and cur and len(cur[1]) < cur[3]):
                if cur is not None and len(cur[1]) >= cur[3]:
                    out.append((cur[0], np.array(cur[1]), np.array(cur[2])))
                    cur = None
                if cur is None:
                    cur = [' '.join(t[1:]), [], [], int(t[0])]
                    continue
            if cur is not None and len(t) >= 3:
                cur[1].append(float(t[1])); cur[2].append(float(t[2]))
                if len(cur[1]) >= cur[3]:
                    out.append((cur[0], np.array(cur[1]), np.array(cur[2])))
                    cur = None
    return out


def export_filter_files(bands, filters_res, outdir, kind='bagpipes'):
    """Write a two-column (Angstrom, throughput) file per band from an EAZY FILTER.RES; returns {band: path}."""
    os.makedirs(outdir, exist_ok=True)
    flt = read_filter_res(filters_res)
    names = [f[0] for f in flt]
    res = {}
    for b in bands:
        k = filters.eazy_filter_index(filters.B[b], names)
        if k is None:
            continue
        _, w, t = flt[k - 1]
        p = os.path.join(outdir, '%s.dat' % b)
        np.savetxt(p, np.c_[w, t])
        res[b] = p
    return res


def command_of(params):
    """params['command'] as an argv list (a string is split like a shell would, without running a shell)."""
    import shlex
    cmd = params.get('command')
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    return list(cmd) if cmd else None
