# External photo-z and SED-fitting codes (plugin `sedcodes`, package `sed_adapters/`)

Adapters for **EAZY** (photo-z), **CIGALE**, **Bagpipes** and **Prospector** (SED fitting), written in the style of `ai_bridge/`: every adapter is a pure function
`process(records, params, task) -> (results, model)` over the bridge's standard input records (so a service could call it), reachable as

* `python_callable`: `sed_adapters.<code>_adapter:run(records, params, context)`,
* `local_command`: `python -m sed_adapters.<code>_adapter` (bridge JSON on stdin, `{"results": [...], "model": ...}` on stdout),
* profile files: `plugins/sedcodes/sedcodes.py --write-profiles FILE` writes `ai_bridge` profiles `sedcodes_eazy|cigale|bagpipes|prospector` (disabled by default; run them with
  `ds9_ai_bridge.py --mode run --services-file FILE --service sedcodes_eazy --task photoz --param python=... --param eazy_data=...` and `PYTHONPATH` = the checkout),
* the plugin (chip "SED codes", Measure tab) and its CLI `plugins/sedcodes/sedcodes.py` (what the recorder/session script replays).

Photometry comes from `MAG_<band>` / `MAGERR_<band>` catalog columns. Bands known to `sed_adapters/filters.py`: HST ACS (F435W ... F850LP), WFC3 IR/UVIS, JWST NIRCam, `SDSS_U..Z`, `J H KS`, `IRAC1 IRAC2`;
`COLUMN:BAND` in the *Magnitude columns* parameter assigns a band to any column (`MAG_AUTO:F160W`). Unknown bands are skipped with a warning, never guessed. Magnitudes (AB) become micro-Jansky
(`f = 10^(-0.4 (m - 23.9))`, `sigma_f = f sigma_m ln10/2.5`); a calibration floor (`min-mag-err`, default 0.02 mag) is added in quadrature.

## Codes and engines

| code | task | native engine | external engine |
| --- | --- | --- | --- |
| EAZY | `photoz` -> `EZ_Z EZ_ZERR EZ_Z16 EZ_Z50 EZ_Z84 EZ_CHI2` | eazy-py (`pip install eazy` + the `eazy-photoz` repository for filters/templates) in the interpreter `python`; p(z) percentiles from the p(z) grid | classic EAZY files: `catalog.cat` (`F_<band>`/`E_<band>`, uJy), `zphot.translate` (filter numbers found by name in `FILTER.RES.info`), `zphot.param`; runs `command zphot.param`, parses `photz.zout` (`z_peak`/`z_m1`, `l68`/`u68`) |
| CIGALE | `sed_fit` -> `SC_LOGM SC_LOGM_ERR SC_LOGAGE SC_LOGAGE_ERR SC_LOGZ SC_AV SC_SFR SC_CHI2` | - | writes `observations.txt` (mJy, columns `id redshift <filter> <filter>_err`), `pcigale.ini` (pdf_analysis), runs `command` (default `pcigale run`), parses `out/results.txt` (`bayes.stellar.m_star`, `bayes.sfh.sfr`, `best.reduced_chi_square`, ...) |
| Bagpipes | `sed_fit` | `native/bagpipes_fit.py` (delayed-exp SFH, Calzetti dust, fixed redshift, Nautilus/MultiNest) in the interpreter `python`; filter curves exported from an EAZY FILTER.RES or `filter_dir/<BAND>.dat` | `command` = any program with the JSON protocol below |
| Prospector | `sed_fit` | `native/prospector_fit.py` (v1 API: parametric SFH, `SedModel`, `CSPSpecBasis`, dynesty; needs python-fsps + `SPS_HOME`) | same protocol |

JSON protocol of the Bagpipes/Prospector fit scripts: stdin `{wd, ids, z[], bands[], flux[][] (uJy, -99 = missing), err[][], filter_files{}, params{}}`, stdout `{"version": ..., "rows": [{LOG_MASS, LOG_MASS_ERR, LOG_AGE, LOG_AGE_ERR, LOG_Z, AV, SFR, SED_CHI2} | {"error": ...}]}`.
Per-object failures (no redshift, fewer than 3 (EAZY) / 4 (SED fit) usable bands) become `error` rows -> empty result columns and a reason in `EZ_NOTE` (photo-z) / `SC_NOTE` (SED fit); a failing external program raises.

The SED-fit steps take the redshift from the catalog column *Redshift column* (default `Z_SPEC`; set `EZ_Z` to chain after the EAZY step). Results/metadata: `<work>/sedcodes/sedcodes_results.json`, `catalog_meta.json` key `sedcodes`, catalog key `sedcodes,results_file`.

## Mock executables and the TOY library (what the tests do and do not show)

`sed_adapters/mocks/` holds test doubles that read/write the *native file formats* (`mock_eazy.py`: EAZY catalog/translate/param -> `.zout`; `mock_cigale.py`: `pcigale.ini` + observations -> `out/results.txt`;
`mock_fit.py`: the JSON protocol) and fit `sed_adapters/toy.py` - five made-up SED shapes, **not physical**. A passing mock test proves unit conversion, band-to-filter mapping, file layout, parsing and entry points, not the accuracy of any code.

## Validation (`plugins/sedcodes/tests`, 12 mock tests ~13 s + real-code tests)

| Test | Result |
| --- | --- |
| EAZY via the external (classic file) interface + mock, 80 TOY galaxies z 0.2-3.2, 7 HST bands | sigma_NMAD 0.0071, outliers 1.3 %, 68 % interval holds truth for 82 %; adapter result equals the direct TOY fit; files: `F_F160W` = `mag_to_ujy`, translate -> filter named `.../f160w.dat`, `Z_STEP`/`APPLY_PRIOR` in `zphot.param` |
| CIGALE (mock), 40 galaxies | log M rms 0.033 dex (median 0.000), Av rms 0.078; observations in mJy (checked), ini keys, redshift column |
| Bagpipes / Prospector protocol (mock), 30 galaxies | log M rms 0.025 dex, median -0.001 |
| entry points | python_callable == local_command result (identical), `ds9_ai_bridge.py --mode run` with a generated profile (25 objects, sigma_NMAD < 0.04), plugin CLI columns/exit codes |
| **real eazy-py 0.8.7** (`pip install eazy`, eazy-photoz templates), 60 galaxies built from EAZY's own templates z 0.2-3.0, mags 21-23.5, 4 % errors, 7 bands | sigma_NMAD 0.0167, outliers (|dz|/(1+z) > 0.15) 6.7 %, 68 % interval coverage 0.77, median reduced chi2 0.30 |
| **real Bagpipes 1.3.6** (Nautilus, n_live 300), 3 Bagpipes model galaxies z 0.4-1.2 | log M(fit) - log M(truth) = +0.018, +0.099, -0.073 dex (pull +0.17, +0.9, -0.66) |
| real Prospector (round 2, item 5) | **exercised**: `plugins/sedcodes/tests/test_real_codes.py::test_real_prospector` (3 synthetic FSPS galaxies from `make_prospector_truth.py`, 7 HST bands, nested sampling, `n_live=250`, 992 s on 1 core): log M - truth = -0.021, +0.035, +0.089 dex (pulls -0.2, 0.28, 0.78); A_V - truth = +0.16, -0.13, +0.24; log age - truth = -0.15, +0.20, +0.27.  Needs (a) FSPS **master** data (`git clone --depth 1 https://github.com/cconroy20/fsps`, `SPS_HOME=<clone>`): the earlier failure was a data/library mismatch (python-fsps 0.5.0 expects 13 MIST metallicities; the v3.2 tag has 12) and (b) **dynesty < 3** (dynesty 3.1 changed the return tuple of its sampler and Prospector 1.4.1's `run_dynesty_sampler` fails with `too many values to unpack (expected 15)`; 2.1.5 works).  3 galaxies is a smoke test, not a calibration |
| real CIGALE | **not exercised**: not on PyPI and `gitlab.lam.fr` was unreachable (TLS handshake fails with `unexpected eof`, also forced TLS 1.2 / IPv4, round 2 item 5 retry; `cigale.lam.fr` answers but the tarballs live on gitlab; the only copy reachable, `JohannesBuchner/cigale`, is a 2016 v0.9 snapshot with a different module set and does not run on Python 3.13); file formats follow the CIGALE documentation; the filter names in `filters.py` (`hst.wfc3.F160W` style) are UNVERIFIED and can be overridden with `params.filter_names` |

## Limits

* The real-code accuracy numbers come from small, easy simulations (templates of the code itself, 3-60 objects); they show the adapters work, not how good the codes are on real data.
* EAZY native engine: priors, template error function (`TEMPLATE_ERROR.eazy_v1.0`), zero-points and extinction use eazy-py defaults of the shipped helper; `PHOTOZ_P50` of the external engine is `z_m1` (mean), not a median.
* Bagpipes/Prospector run at a fixed redshift with a fixed simple parametrisation; `SFR` is the current SFR of that model, `SC_AV` for Prospector is `1.086 * dust2`.
* CIGALE: the module list/grid defaults in `cigale_adapter.py` are a small example, tune via `module_params`; `SC_AV` = `E_BV_lines * av_per_ebv` (4.05 by default).
* Native runs of Bagpipes cost minutes per object; no caching between sessions beyond Bagpipes' own posterior files in the work directory.
