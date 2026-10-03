# Test infrastructure

`scripts/run_all_checks.sh [--quick] [--only a,b] [--list] [--display :77]` runs every harness and prints a table
(CHECK / PASS|FAIL|SKIP / seconds / last log line); exit status 1 if anything FAILs (SKIP is not a failure).  Logs: `$OGF_CHECK_OUT`
(default `/tmp/ogf_checks.<pid>/NAME.log`).  GUI checks start an Xvfb on the display if none runs there.

| check | what | needs | in `--quick` |
|---|---|---|---|
| `tools_syntax` | `py_compile` of `tools/*.py`, `scripts/*.py`; `tclprocs.tcl`; `compare_menus.py` on the stored dumps | python, tclsh | yes |
| `ai_bridge_tests` | `pytest ai_bridge/tests` | venv | yes |
| `moving_tests` | `pytest moving/tests` | venv | yes |
| `ai_gui` | `scripts/verify_ai_gui.py` (AI panels under Xvfb) | Xvfb | yes |
| `agent_gui` | `scripts/verify_agent_gui.py`: Backend dropdown, detection, confirmation text, consent flags, recorder, replay and pipeline mode of the agent-CLI backends, **with fake CLIs on PATH** (34 checks) | Xvfb | yes |
| `agent_real` | `scripts/verify_agent_cli_real.py`: runs the *real* `codex` / `claude` / `agy` / `gemini` / `grok` binaries if installed, in an empty HOME without any login: every flag is accepted, the not-logged-in error is parsed. **No model answer is obtained.** SKIP (77) when none is installed | any installed CLI | yes |
| `icl_export` | `scripts/verify_icl_export.tcl` smoke on m51 | Xvfb, `/workspace/fits/m51.fits` | yes |
| `click_chooser` | `scripts/verify_click_chooser.tcl` (click procs called directly, layout geometry) | Xvfb | yes |
| `click_xevent` | `scripts/verify_click_xevent.tcl` (real X events with xdotool; SKIP if xdotool absent) | Xvfb, xdotool | yes |
| `cat_api` | `scripts/verify_cat_api.tcl` (unit test of the `::ogf::cat` accessor: get/set/trace/registry) | Xvfb | yes |
| `cat_behavior` | `scripts/verify_cat_behavior.sh` (82 features: exec argv, `catpanel` keys, `.prf` files, session steps and status text, diffed against `scripts/golden/cat_behavior.golden`; `--update` only from deliberately reviewed code) | Xvfb | yes |
| `report_tests` | `plugins/report/tests` (pytest: HTML report contents on m51 and a HUDF F160W crop, cut-outs, escaping, optional PDF) | python deps, `bin/ds9_sextract` | yes |
| `geometry` | `scripts/verify_geometry.tcl` (table y/h and info-area height 181/769/154 in 21 panel states: 6 tabs, extract, search, all time-domain views, review filter, selection, tile/single, detached/reattached; main width 1300 -> 736 -> 1300) | Xvfb | yes |
| `mouse` | `scripts/verify_mouse.tcl` (real xdotool events: notebook tab, chip menu + entry, cascade submenu, table row click, header-click sort, wheel, Tools-menu detach/reattach; 25 checks) | Xvfb, xdotool | yes |
| `newplugins` | `scripts/verify_newplugins.sh`: GUI check (isophote, completeness, daophot, psfex, multifit, extended-morphology, noise-model, sedcodes, cluster, spectra, xmatch, lightcurves, batch, repro sections, chip-row overflow into the "More" menu) of the analysis plugins added after the restructuring (`scripts/verify_newplugins.tcl`: real m51 extraction, each step through `::ogf::step::run`, outputs/frames/windows, layout invariants) **and** headless replay of the exported session script, requiring identical stdout and catalog CRCs per step | Xvfb | yes |
| `isophote_tests` | `pytest plugins/isophote/tests` (synthetic Sersic e/PA/intensity/growth curve/boxy-disky/mask, CLI, M51 smoke) | venv + photutils | yes |
| `completeness_tests` | `pytest plugins/completeness/tests` (analytic SEP recovery curve, false positives, flux conservation, mask, callback loading, logistic fit, FITS metadata, CLI) | venv + sep | yes |
| `daophot_tests` | `pytest plugins/daophot/tests` (18 tests: PSF model spatial variation / fallbacks / round trip, FIND cuts, PHOT apertures + sky modes, PICKPSF, ALLSTAR flux/position/completeness/chi/pulls, order-2 vs constant PSF, flags, SUBSTAR/ADDSTAR, aperture correction, DAOMATCH, artificial-star test, CLI) | venv + photutils | yes |
| `psfex_tests` | `pytest plugins/psfex/tests` (9 tests: recovery of a spatially varying PSF (FWHM / ellipticity vs truth at 15 field positions), unit-sum stamps, order 0 vs 2, fallbacks (no stars / 4 stars), JSON/FITS round trip, CLI + catalog columns + tables + star list, PSF photometry and crowded photometry with the model vs one constant PSF, tiled ZOGY with model fields vs constant PSF) |
| `multifit_tests` | `pytest plugins/multifit/tests` (11 tests: single Sersic n=1/2/4 recovery, pull calibration of the parameter errors over 40 noise realisations, blends vs separation (simultaneous vs masked vs ignored neighbour), bulge+disk B/T, PSF+host, fixed / bounds / tie / sky plane / mask, spatially varying PSF model at the field edge, auto model selection, Poisson weights, catalog CLI on a synthetic field with a blend pair and a star, config mode) |
| `morphext_tests` | `pytest plugins/morph_ext/tests` (13 tests: Petrosian radii (eta 0.1/0.2/0.3), flux fractions R20-R90, concentration and Petrosian flux vs the analytic Sersic functions for 4 profiles / axis ratios, Kron flux vs analytic aperture fraction, noise robustness and smoothness noise correction, smooth vs clumpy, Petrosian-segmented Gini/M20, masked companion, edge flags, driver on a synthetic field) |
| `noisemodel_tests` | `pytest plugins/noisemodel/tests` (8 tests: pixel ACF vs exact discrete truth, blank-aperture sigma_N(r) vs analytic sum over the ACF, white-noise law, flux-error pulls corrected vs naive on 400 faint stars, gradient + scattered-light model error (mesh / poly / mesh+poly), flat-field residual chi2, non-uniform rms map, real HUDF F160W crop) |
| `sedcodes_tests` | `pytest plugins/sedcodes/tests` (12 mock tests: photometry conversion, filter registry vs EAZY filter lists, EAZY/CIGALE/Bagpipes/Prospector adapters against mock executables with a TOY truth, native file formats, python_callable / local_command / ai_bridge profile / plugin CLI; plus real eazy-py and Bagpipes accuracy tests that are skipped when those are not installed) |
| `cluster_tests` | `pytest plugins/cluster/tests` (12 tests: red-sequence recovery over 40 synthetic clusters, membership purity/completeness and p_mem calibration, Poisson calibration of the density significance map, peak recovery with masked area, arc length/radius/alignment vs truth with decoys, faint-arc limits, real HUDF crop null test and injected arcs; ~22 s) |
| `spectra_tests` | `pytest plugins/spectra/tests` (9 tests: line-fit pulls, continuum, blind redshift vs strength, 2D boxcar/optimal extraction, IFU cube moments, CLI, real SDSS plate when present; ~42 s) | venv | yes |
| `xmatch_tests` | `pytest plugins/xmatch/tests` (10 tests: completeness/purity vs radius, chance-match estimate vs truth, systematic shift recovery, ambiguity flags, pixel mode, table readers, CLI file mode, fake local TAP server, live Gaia DR3 (skipped offline); ~20 s) | venv | yes |
| `lightcurves_tests` | `pytest plugins/lightcurves/tests` (9 tests: held-out 7-class accuracy and SN selection, calibration, S/N dependence, host offset, injection through difference-image photometry, readers, CLI link modes, plot, reproducible training; ~60 s) | venv | yes |
| `batch_tests` | `pytest plugins/batch/tests` (6 tests: cli template expansion, serial vs parallel, resume, failure isolation, recipe validation, 60-field throughput; ~20 s) | venv | yes |
| `repro_tests` | `pytest plugins/repro/tests` (8 tests: bundle contents, batch round trip verify, tampered input / bundle / processing, table comparison tolerance, environment diff, catalog bundle with included inputs; ~15 s) | venv | yes |
| `review_td` | `scripts/verify_review_td.tcl` (review on moving / transient / detection rows and the All view: columns, filters, tint, steps, save/load, re-run pruning, report from a time-domain view; 55 checks) | Xvfb | yes |
| `review_gui` | `scripts/verify_review_gui.tcl` (review columns, table filter and tint, session steps, save/load round trip, export job, layout invariants; 51 checks) | Xvfb | yes |
| `session_replay` | `scripts/verify_session_replay.sh`: record GUI sessions, export scripts, replay, compare (78 checks) | Xvfb, venv, ~10+ min | no |
| `link_bench` | `moving/validation/link_bench.py` on the injected sets `/workspace/work/inj1-6.json.pkl` (SKIP when missing); when the regenerated sets `/workspace/work/i1/sets/inj1..12.json.pkl` exist (incl. CR-heavy 9-12, held-out 7-8) all twelve run and per-group totals are printed | venv | no |
| `moving_options` | `scripts/verify_moving_options.tcl` (Xvfb): ZOGY parameters of the moving plugin -> CLI flags, default argv unchanged, recorded argv, layout invariants | X server | yes |
| `manifests` | `tools/validate_manifests.py` (static check of all plugin manifests and `cli` templates) + `pytest tools/tests` (13 mistake classes are detected) | venv | yes |
| `cli_templates` | `scripts/verify_cli_templates.tcl` (Xvfb, fake interpreter `scripts/fake_python_for_templates.sh`): the five converted steps vs. their legacy procs - argv, recorder record, resulting catalog, status | X server | yes |
| `cli_headless` | `scripts/verify_cli_headless.tcl` (Xvfb, fake interpreter): legacy GUI proc vs. the step's `headless` block for 14 steps (extract, dual, multiband, crossmatch, segmap, completeness, photo-z, SED, deconvolution, ICL ×3, mask.auto, LSBG) - argv, recorder record, catalog and status identical. The real-driver counterpart `scripts/verify_headless_real.py` (m51, ~2 min, not in `run_all_checks.sh`) compares catalogs/images of the GUI path and the batch path byte for byte. | X server | yes |
| `regression_tests` | `pytest regression` (4 offline tests of the regression-set machinery: range comparison, every case has a baseline, every case is documented here, unavailable data -> SKIP) |
| `regression_data` | (long) `regression/run_regression.py`: 11 per-tool cases on small public data compared with `regression/baselines.json` (see "Regression-validation set" below); a case whose data cannot be fetched is SKIP, the check is skipped when nothing could run |
| `moving_session` | self-contained by default (`scripts/verify_moving_session_auto.sh`): a real ds9 runs the Moving Objects GUI steps on the 4 cached BB89 exposures (`scripts/verify_moving_record.tcl`), saves the session script, and `verify_moving_session.sh` replays it and compares detections / movers / tracklets / identified / orbit / transients / light curves / export byte for byte; with `OGF_MOVING_SESSION`, `OGF_MOVING_REF`, `OGF_MOVING_FIELD` set it replays that session instead. SKIPs (77) without the cached exposures, an X server, or (network or `~/.ds9/moving_cache`); ~10 min | cache (+ network when uncached) | no (long) |

Test data location: `OGF_TEST_FITS` (default `/workspace/fits`).  `docs/windows_macos_build.md` lists what is *not* verified on
other platforms.  The generator/compare scripts are in `tools/` (see `tools/README.md`).

## Notes on the analysis-plugin checks (items A-N)
* `newplugins` runs one GUI section per plugin (each step through `::ogf::step::run`, outputs, windows, layout invariant), exports the session script, replays it headless (stdout and catalog CRC per step identical) and, since the `repro` section,
  re-runs the GUI-made reproducibility bundle (`plugins/repro/repro.py verify`, session replay) and requires the re-run catalog to match the GUI catalog.  `OGF_NP_ONLY=lightcurves,batch,repro` limits it to some sections.
* `batch_tests` and `repro_tests` need `bin/ds9_sextract` (skipped otherwise); `xmatch_tests` contains one live Gaia DR3 TAP test (skipped when offline); `spectra_tests` has a real SDSS test (skipped when `/workspace/fits/sdss` is absent);
  `sedcodes_tests` uses mock codes plus the real eazy-py / Bagpipes when installed.
* Plugin tests print their measured numbers (`pytest -s`); the values quoted in `docs/<plugin>.md` and `docs/progress_log.md` come from those runs.
* Do not edit repository files or plugins while `run_all_checks.sh` is running: the stages read them live.


## Regression-validation set (`regression/`, checks `regression_tests` and `regression_data`)

One or two cases per analysis tool, run on small **real public data** (not synthetic images), each producing a few headline numbers that must stay inside the
ranges of `regression/baselines.json`.  `python3 regression/run_regression.py [--list] [--only a,b] [--tool t] [--report out.json] [--record]`
(`--record` prints the observed metrics instead of comparing - use it to write or refresh a baseline).  Exit 1 = a metric left its range or a case crashed.

**Data and data path.**  Everything lives in `$OGF_DATA_CACHE` (default `~/.cache/ogfinder_regression`) and is fetched on demand, once:
* HUDF12 WFC3/IR F160W mosaic (HLSP `hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits`, https://archive.stsci.edu/hlsps/hudf12/, ~0.5 GB; `/workspace/fits/hudf_f160w.fits` or `$OGF_FITS_DIR/hudf_f160w.fits` is used instead when present) -> a 1400 x 1400 crop `hudf12_f160w_crop.fits` (pixel box [1000:2400,1000:2400]) is cut and cached;
* SDSS J0946+1006 (SLACS), HST ACS/WFC F814W, MAST product `j9op14010_drc.fits` (program 10886; 215 MB) -> 500 x 500 cutout in `lens/`;
* SDSS Stripe 82 run 94 frames and the SkyBoT truth for camcol 3 in `sdss_run94/` (written by `moving/validation/sdss_known_asteroids.py`; NOT downloaded automatically - the case is SKIP until that script has been run once);
* `photo_z/data/sdss_specphoto.h5` + the MDN checkpoint in the repository (the photo-z case needs torch).
Unreachable data (offline box, MAST down) give SKIP, never FAIL.

| case | tool | what is measured (headline metrics) |
|---|---|---|
| `extract_hudf` | `bin/ds9_sextract` | source count, count brighter than 24, brightest and median MAG_AUTO on the crop |
| `noise_blank_apertures` | noisemodel | blank-aperture pull std with the local noise law (ideal 1) and with the naive law; pixel sigma |
| `photerr_injection` | noisemodel / `ogfkit.photerr` | injected-source pull std (isolated) and naive pull std, fraction measured |
| `depth_hudf` | completeness | median, p10, p90 5-sigma limiting magnitude in a 3-pixel-radius aperture |
| `psf_heldout_residual` | psfex | residual of held-out stars for PSF degree 0/1/2/3 and rank 1 |
| `stack_null_calibration` | stacking | z-score mean/std of 17 null stacks (calibration of the quoted errors) |
| `isophote_bright_galaxy` | isophote | median over 3 bright galaxies of number of isophotes, total magnitude, outer ellipticity / PA |
| `multifit_sersic` | multifit | 6 faint extended galaxies, Sersic fit: median Re, n, offset from MAG_AUTO |
| `moving_sdss_run94` | moving | SDSS run 94 camcol 3, orbit-prior linking: true tracklets, number of tracklets, precision lower bound, recall |
| `photoz_closure_sdss` | photoz_sed | PIT KS p, 68 % coverage and CRPS before / after the recalibration map and cross-fitted; injected-width recovery |
| `lens_j0946_evidence` | lensmodel | HST J0946+1006: chi2_red of the stored SIE+shear model and the evidence loss when theta_E +3 % or q -0.05 |

**What these baselines are.**  They are *self-baselines*: the ranges were set around the values the current code produces on this data (generous, 3-10 sigma of the
case-to-case spread seen between runs or wide where few objects enter), so they detect regressions and silent behaviour changes, **not** accuracy against truth.  Accuracy
against truth is established by the dedicated validations (`docs/progress_log.md` R1-R22; GALFIT cross-check, SkyBoT comparison, SDSS photo-z truth, literature lens model).
Several cases deliberately record uncomfortable values (e.g. `stack_null_calibration` z mean about +0.7, `photerr_injection` only 7 evaluable injections) - the ranges are wide there and the
docs of the tool say why.

**Not covered** (no public-data regression case yet): ai_bridge, cluster (red-sequence), spectra, sedcodes (needs the EAZY/FSPS environments; the photo-z case disables EAZY), lightcurves, xmatch, daophot PSF photometry,
morph_ext, batch/repro, the GUI plugins themselves (covered by the Xvfb checks above) and the GALFIT wrapper (GALFIT binary not redistributable).  Typical run time: about 4-5 minutes with the data cached on a quiet box
(`moving_sdss_run94` 1-3 min, `multifit_sersic` 1 min, `photoz_closure_sdss` 1 min, the HUDF cases a few minutes).
