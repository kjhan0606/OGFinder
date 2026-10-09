# Astrafex Core in Astrafex Desktop - the shared compute core, stand-alone shell (plugin `ds10core`)

Names: the product family is **Astrafex** (Astrafex Web, **Astrafex Desktop** = this OGFinder tool set with its plugins, **Astrafex Core**). The plugin id `ds10core`, its
folder, the parameter file `~/.ds9/ogf_params/ds10core.json`, the session records `ds10core.*` and the Python package `ds10core` keep their original names (existing
settings, sessions and recorded command lines depend on them; the package is also importable as `astrafex_core`); the menu/chip label, titles and help texts say
"Astrafex Core". Old names still work: `scripts/ds10`, `DS10_*` variables, `ds10-script/1|2`, `ds10-offline/1`, `ds10-return/1` files. See Astrafex Web `docs/naming.md`.

Design rule (docs/phase2_web_design.md, chapter 12): one compute core, one plugin manifest format, two shells, shared parity tests. The web version
(`Astrafex Web`, private repository) gained regions with exact pixel statistics and masks, a pixel table, weight/rms/exposure/mask map association,
header-driven AB zero points, multi-band forced photometry and `astrafex-script/1` step lists / session replay. These now live in a Python package
**`ds10core`** (numpy + scipy only, original code) in the Astrafex Web repository (`server/ds10core`; full description: Astrafex Web `docs/shared_core.md`).
OGFinder is the second shell.

## Licence boundary
* `ds10core` contains no ds9, AST, SEP or other GPL/LGPL code and imports nothing from OGFinder (checked by Astrafex Web `test_core_independent.py`).
* OGFinder does not import or copy it: `plugins/ds10core/ds10.py` locates the package and starts `python -m ds10core ...` (= `python -m astrafex_core`) as a **separate process**
  (`plugins/ds10core/tests` checks that the launcher never imports it). The manifest carries the tags (`"license": {...}`, `"web": false`).
* Where the core is looked for: `--core-dir` (plugin parameter *Core directory*), `DS10_CORE`, `~/.ds9/ogf_params/ds10core.json` (`{"core-dir": ...}`),
  `../ds10-web/server`, `/workspace/ds10-web/server`. Interpreter: `DS10_PYTHON`, else the OGFinder venv (numpy + scipy needed). Without the core the
  plugin still loads; its steps fail with a clear message, the tests and the `ds10core_*` checks SKIP.

## Command line (`scripts/astrafex ...` = `python3 plugins/ds10core/ds10.py ...`; `scripts/ds10` is the former name and forwards)
| command | what |
|---|---|
| `astrafex info IMAGE [--text]` / `calib IMAGE [--zp X]` | header summary, WCS, AB zero point with origin (header PHOTFLAM/PHOTPLAM/PHOTZPT, manual override, flagged fallback) |
| `astrafex maps IMAGE suggest\|auto\|associate\|show\|clear [--dir D --map F --kind K]` | weight / rms / exposure / mask map association (sidecar `IMAGE.ds10maps.json`) |
| `astrafex regions convert\|stats\|mask REGFILE --image IMG [...]` | DS9 `.reg` / JSON / CSV, frames image/fk5, exact pixel membership statistics (background region, median/mean), mask FITS |
| `astrafex pixtab IMAGE --xy X Y [--size N --mode value\|ra\|dec --csv F --json]` | pixel table (float32 data, WCS per cell) |
| `astrafex forced IMG IMG [IMG...] --out DIR [--psf-match]` | detection image -> detect -> forced photometry (aperture / Kron-like / isophotal, correlated-noise errors, upper limits, aperture corrections) -> colours; optional `--psf-match` homogenises bands for colours (Aniano/Boucaud kernels; see Astrafex Web `docs/psf_matching.md`); writes `forced_catalog.*`, `detect_*`, `detection*.fits`, `script.json` |
| `astrafex run SCRIPT [--files DIR\|FITS...]` | run a `astrafex-script/1` step list, a web-exported `replay.py` (`--inputs`), or a bundle `.zip` (replay local mode) |
| `astrafex open BUNDLE [--dest DIR --force --text]`, `astrafex status [WORKSPACE]` | restore an **offline bundle** of the web app (Offline export button) into a workspace folder and show it; see below |
| `astrafex validate SCRIPT`, `export-script SOURCE`, `compare A B`, `list` | validate against the manifests; replay/bundle -> step list; sha256 + tolerance comparison; plugins/steps |

## GUI (what is done)
The plugin is a normal OGFinder plugin: manifest `plugins/ds10core/plugin.json`, tab **Measure**, chip "Astrafex Core" (it sits in the chip row's **More** menu
because the Measure row is full), one menu entry per step, parameter dialog generated from the manifest (Settings...). Steps: image info, AB zero point
(▶ of the chip), maps, region statistics, regions -> mask FITS, convert regions, pixel table at (X, Y), forced photometry (band list typed in the dialog),
run script, **open offline bundle**, LSBG finder, **PSF photometry / CMD** (`psf-phot`: `run_psf_phot.py` -> `python -m ds10core.tools.psf_phot`;
second band, single exposures and short exposures are given in the Settings group *PSF photometry*; the catalogue fills the table), **Sky catalogue query** (`catalog-query`: Gaia DR3 / SDSS / 2MASS / DESI / LAMOST / Rubin around the open image; the catalogue fills the table and `markers.reg` is drawn on the frame), **Cross-match catalogue** (`CX_*` columns added to the table), **MAST search** / **MAST download**, **Run user plugin** (`astrafex-user-plugin/1`). Results open in the text window; every run is recorded by the session recorder (`ds10core.*`) with its exact argument vector.
**Not done:** no interactive region editing/drawing beyond ds9's own (regions are used from `.reg` files - ds9's own Region menu saves them), no table view of the forced
catalogue inside the catalogue panel (the files are written and the path is shown), no band picking from the loaded frames, no mask hand-over to the shared
mask manager, no GUI progress strip for long runs beyond the job runner's.
GUI check: `scripts/verify_ds10core.sh` (real `bin/ds9`, 41 checks; the catalogue query runs against a local fake TAP service, `scripts/fake_sky_tap.py`,
built on Astrafex Web's `server/tests/catmock.py`, so no network is needed); screenshots `docs/shots/astrafex_*.png`:
`astrafex_01_chip_menu_measure_tab`, `astrafex_01b_plugin_menu`, `astrafex_02_calibration_result`, `astrafex_03_region_statistics`,
`astrafex_04_pixel_table`, `astrafex_05_forced_photometry_result`, and `sky_catalog_desktop_1440x900`, `sky_catalog_desktop_params_1440x900` (catalogue table + markers, Settings group *Catalogue*).

## Continue a web session here (offline bundle)
When the network to the Astrafex Web server drops, a bundle made with the web app's **Offline export** button (or kept up to date by its automatic refresh) restores the
session on this computer; the full workflow, limits and the way back are in Astrafex Web `docs/offline_workflow.md`.
* Command line: `astrafex open session-offline.zip --dest ~/work/field` (checksums and the event hash chain are verified first; an existing non-empty folder is only replaced
  with `--force`). It creates `files/` (the images: **originals** for pro and above, float32 tool copies for free, as the bundle says), `regions/regions.reg`,
  `session_script.json`, `results/<job>/` (stored outputs), map/calibration sidecars `<image>.ds10maps.json`, `workspace.json`. `astrafex run field/session_script.json`
  then **reuses the stored result of every step whose text and inputs did not change** (status `restored`) and runs the steps you add at the end; each run is written to
  `runs/rNNNN/`. `astrafex status field` lists the runs.
* GUI: Measure tab, chip **Astrafex Core** menu, **Open offline bundle (continue a web session here)**. Parameters *Offline bundle* and *Folder to restore into* (empty =
  `~/astrafex-offline/<bundle name>`). After restoring, the first image of the session is loaded in ds9 with the regions of the session (`regions/regions.reg`), the
  summary opens in the text window and the status line names the workspace; then **Run ds10-script** with `<workspace>/session_script.json` and *Folder with the images* =
  `<workspace>/files` continues the work. Screenshots: `docs/shots/astrafex_offline_open_bundle_ds9.png`, `astrafex_offline_open_bundle_script_run.png` (1440x900).
* The way back: **`astrafex pack field --text`** (GUI: Measure tab, chip **Astrafex Core** menu, **Pack results for the web server (return bundle)**, parameters *Workspace folder* and
  *Return bundle file*) writes `<name>-return.zip` with the script as it stands and the results (`results/<step>/...`) of the steps that were run **here**, each with the
  checksums of its inputs and outputs. Give that file to the web app (Offline export panel, *Import a return bundle...*): the server checks everything again and opens a
  **new** session "... (continued locally)" in which those steps are marked *local*; the original session is never changed and the server does not recompute the steps.
  A step without a result from this computer is reported and stays out. Images used by a step must be files of the web account (same content); files that exist only here
  must be uploaded first. Details and limits: Astrafex Web `docs/offline_workflow.md`.
* Paths with spaces: the *Band images* list of the forced-photometry dialog is split at `|` or new lines when the text contains one (names may then contain spaces), else at spaces;
  a name with spaces can be put in double quotes.
* Limits: the bundle is a snapshot from the moment it was made; for multi-extension originals only the first image HDU is restored; bit-identical new results need the same
  Python/numpy/scipy versions as the server (the bundle records them).

## Parity with the web version
Same functions, same results. Test list (Astrafex Web `server/tests/test_shared_core.py`, `test_core_independent.py`; OGFinder `plugins/ds10core/tests`) and the
HUDF F105W/F125W/F160W web-versus-stand-alone batch (script exported from a web session, run on both): see Astrafex Web `docs/shared_core.md`. Summary: all
catalogues, detection images, segmentation, zero points, region statistics/masks and pixel tables have the **same sha256**; run-time fields (`seconds`) are
the only differences (identical Python/numpy/scipy versions).

### Star aperture photometry (steps `aper-phot`, `aper-series`)
Menu entries **Aperture photometry (click stars)** and **Aperture photometry light curve (frames)**: `run_aper_phot.py` ->
`python -m ds10core.tools.aper_phot` (the same tool as the web panel *Photometry*; research notes, feature list and specification in the Astrafex Web
repository, `docs/aperture_photometry.md`).
1. Click on the stars in ds9 (Region edit mode, circle or point regions; the click need not hit the centre: the centroid finds the star).
   The `before` hook (`OGFDs10coreAperBefore`) saves the frame's regions in image coordinates to `~/.ds9/aper_picks.reg`.
   *Stars* = `regions` (default) uses them, `catalog` uses the catalogue table (X_IMAGE/Y_IMAGE or RA/Dec, e.g. a Gaia query), `detect` finds stars.
2. Settings group *Aperture photometry*: radii (default 1,1.5,2,3) and unit (fwhm / px / arcsec; also the unit of the sky annulus, default 4-6),
   FWHM (auto), sky estimator (median / mode / mean, sigma-clipped, neighbours masked), centroid (com / gauss / none) and recentring limit,
   gain (0 = header), saturation, aperture correction (cog / infinity / none), zero point (auto = header, a number, or `fit` with a reference column).
3. Result: the catalogue table (`NUMBER NAME X_IMAGE Y_IMAGE ... FLUX_APER_k MAG_APER_k ... FLUX FLUXERR MAG MAGERR SNR APCOR FLAGS`), the picks
   replaced by the apertures (green) and sky annuli (red) (`after` hook `OGFDs10coreAperAfter`, `~/.ds9/aper_apertures.reg`), a status line with
   the number of stars, FWHM and zero point. Files: `~/.ds9/aperphot_{catalog,summary,cog}.*`.
4. Light curve: Settings group *Light curve*: frames (comma-separated paths; the open image is frame 1), target / comparison / check star numbers
   (order of the picks), tracking (wcs refined by the bright-star shift / offset / none), AAVSO star name and observer code. The light curve
   (target - ensemble, check star) opens in its own window; `~/.ds9/aperseries_{lightcurve,series}.tsv` and `aperseries_aavso.txt` (AAVSO Extended).

GUI check: `scripts/verify_aper_phot.sh` (real `bin/ds9`; six stars picked with **real X clicks** through xdotool when available, otherwise point
regions loaded into the frame; flux within 3 % of the injected value, centroids < 0.2 px, table, 6 annuli + 24 apertures on the frame, a
6-frame light curve with an 8 % variable recovered to 0.2 %, a flat check star, AAVSO file, light-curve window); `run_all_checks.sh` check
`aperphot_gui`; driver tests `plugins/ds10core/tests/test_aper_phot_driver.py` (in `ds10core_tests`). Screenshots
`docs/shots/aper_desktop_table_1440x900.png`, `docs/shots/aper_desktop_lightcurve_1440x900.png` (all four, including the parameter tab and the
overlay zoom, are in the Astrafex Web user guide, `docs/user_guide/figures/aper_desktop_*`).

### Saturation mask, varying PSF, photometric solution, star-galaxy

Four more menu entries on the same chip. Each driver locates the core and starts it as a separate process. Settings groups use the same names.

| Menu | Driver | Core module | What the step writes |
|---|---|---|---|
| Saturation mask | `run_pixmask.py` | `ds10core.tools.pixmask` | `{work}/pixmask.fits` (uint8: saturation 1, bleed 2, persistence 4, DQ 8) and `{work}/pixmask_summary.json`. The status line opens in the text window. |
| Varying PSF | `run_psfvar.py` | `ds10core.tools.psfvar` | `{work}/psf_model.json`, `{work}/psf_center.fits`, `{work}/psfvar_summary.json`. An open catalogue supplies the stars. With no catalogue, or an empty catalogue file, the core uses compact detections. |
| Photometric solution | `run_photcal.py` | `ds10core.tools.photcal` | Fits `m_std = m_inst + zp + k X + c C` on the open catalogue and replaces the table. Files: `{work}/photcal.tsv`, `{work}/photcal_summary.json`. |
| Star-galaxy | `run_stargal.py` | `ds10core.tools.stargal` | Half-light radius and concentration against the PSF. Adds `CLASS_SG`, `STELLARITY`, `R50`, `R80`, `CONCENTRATION` and leaves the existing `CLASS` column as it is. Files: `{work}/stargal.tsv`, `{work}/stargal_summary.json`. |

Bleed and persistence are on unless the Settings booleans are cleared. An empty Saturation value uses the `SATURATE` card, then `SATLEVEL`. An empty PSF FWHM estimates the PSF from the bright round sources and needs at least eight of them.

### Original-code tools on the same menu

Ten more menu entries. Each driver locates the core and starts it as a separate process. CCD reduction is already the `ccdreduce` plugin and is not repeated here. Plugins that still need a licence review stay out of this menu.

| Menu | Driver | Core module | What the step writes |
|---|---|---|---|
| Solve WCS | `run_wcs.py` | `ds10core.tools.wcs_solve` | `{work}/wcs_solved.fits`, `{work}/wcs_summary.json`. The status line opens in the text window. |
| Dither stack | `run_dither.py` | `ds10core.tools.dither_stack` | `{work}/dither_stack.fits`, `{work}/dither_summary.json`. Frame paths are joined by a vertical bar. |
| Resample | `run_regrid.py` | `ds10core.tools.regrid` | `{work}/resampled.fits`, `{work}/regrid_summary.json`. Bilinear sampling onto the reference shape and WCS. |
| Cosmic-ray mask | `run_crmask.py` | `ds10core.tools.cr_mask` | `{work}/crmask.fits`, `{work}/crmask_summary.json`. The step id is `crmask`. |
| Detect sources | `run_detect.py` | `ds10core.tools.detect` | Replaces the table with `{work}/detect_catalog.tsv`. Also `{work}/detect_segmap.fits` and `{work}/detect_summary.json`. |
| Galactic extinction | `run_extinct.py` | `ds10core.tools.galactic_ext` | Replaces the table. Files: `{work}/extinct.tsv`, `{work}/extinct_summary.json`. |
| Periodogram | `run_period.py` | `ds10core.tools.periodogram` | Replaces the table with `{work}/periodogram.tsv`. Summary: `{work}/period_summary.json`. |
| Lyman-break | `run_dropout.py` | `ds10core.tools.dropout` | Replaces the table with `{work}/dropout.tsv`. Summary: `{work}/dropout_summary.json`. |
| Completeness and counts | `run_counts.py` | `ds10core.tools.counts` | Replaces the table with `{work}/counts.tsv`. The step id is `counts`. Summary: `{work}/counts_summary.json`. |
| Luminosity function | `run_lf.py` | `ds10core.tools.lumfunc` | Replaces the table with `{work}/lf.tsv`. Summary: `{work}/lf_summary.json`. |

The WCS catalogue file is a reference list of ra and dec. The default survey is `none`, which does not query the network. Choosing `gaia_dr3` does. Deblend stays off unless that Settings boolean is turned on. Completeness keeps the short defaults of 4 magnitude bins and 6 objects per bin. An empty catalogue file on that step is ignored. The luminosity function skips the Schechter fit below 30 galaxies. No dust map is read.

## Checks
`scripts/run_all_checks.sh --only ds10core_tests,ds10core_gui,aperphot_gui`; `tools/validate_manifests.py` (0 problems).
Driver tests: `plugins/ds10core/tests/test_image_tools_driver.py`, `plugins/ds10core/tests/test_core_tools_driver.py`.
