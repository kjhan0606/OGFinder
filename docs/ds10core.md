# ds10core - the shared compute core, stand-alone shell (plugin `ds10core`)

Design rule (docs/phase2_web_design.md, chapter 12): one compute core, one plugin manifest format, two shells, shared parity tests. The web version
(`ds10-web`, private repository) gained regions with exact pixel statistics and masks, a pixel table, weight/rms/exposure/mask map association,
header-driven AB zero points, multi-band forced photometry and `ds10-script/1` step lists / session replay. These now live in a Python package
**`ds10core`** (numpy + scipy only, original code) in the ds10-web repository (`server/ds10core`; full description: ds10-web `docs/shared_core.md`).
OGFinder is the second shell.

## Licence boundary
* `ds10core` contains no ds9, AST, SEP or other GPL/LGPL code and imports nothing from OGFinder (checked by ds10-web `test_core_independent.py`).
* OGFinder does not import or copy it: `plugins/ds10core/ds10.py` locates the package and starts `python -m ds10core ...` as a **separate process**
  (`plugins/ds10core/tests` checks that the launcher never imports it). The manifest carries the tags (`"license": {...}`, `"web": false`).
* Where the core is looked for: `--core-dir` (plugin parameter *Core directory*), `DS10_CORE`, `~/.ds9/ogf_params/ds10core.json` (`{"core-dir": ...}`),
  `../ds10-web/server`, `/workspace/ds10-web/server`. Interpreter: `DS10_PYTHON`, else the OGFinder venv (numpy + scipy needed). Without the core the
  plugin still loads; its steps fail with a clear message, the tests and the `ds10core_*` checks SKIP.

## Command line (`scripts/ds10 ...` = `python3 plugins/ds10core/ds10.py ...`)
| command | what |
|---|---|
| `ds10 info IMAGE [--text]` / `calib IMAGE [--zp X]` | header summary, WCS, AB zero point with origin (header PHOTFLAM/PHOTPLAM/PHOTZPT, manual override, flagged fallback) |
| `ds10 maps IMAGE suggest\|auto\|associate\|show\|clear [--dir D --map F --kind K]` | weight / rms / exposure / mask map association (sidecar `IMAGE.ds10maps.json`) |
| `ds10 regions convert\|stats\|mask REGFILE --image IMG [...]` | DS9 `.reg` / JSON / CSV, frames image/fk5, exact pixel membership statistics (background region, median/mean), mask FITS |
| `ds10 pixtab IMAGE --xy X Y [--size N --mode value\|ra\|dec --csv F --json]` | pixel table (float32 data, WCS per cell) |
| `ds10 forced IMG IMG [IMG...] --out DIR` | detection image -> detect -> forced photometry (aperture / Kron-like / isophotal, correlated-noise errors, upper limits, aperture corrections) -> colours; writes `forced_catalog.*`, `detect_*`, `detection*.fits`, `script.json` |
| `ds10 run SCRIPT [--files DIR\|FITS...]` | run a `ds10-script/1` step list, a web-exported `replay.py` (`--inputs`), or a bundle `.zip` (replay local mode) |
| `ds10 open BUNDLE [--dest DIR --force --text]`, `ds10 status [WORKSPACE]` | restore an **offline bundle** of the web app (Offline export button) into a workspace folder and show it; see below |
| `ds10 validate SCRIPT`, `export-script SOURCE`, `compare A B`, `list` | validate against the manifests; replay/bundle -> step list; sha256 + tolerance comparison; plugins/steps |

## GUI (what is done)
The plugin is a normal OGFinder plugin: manifest `plugins/ds10core/plugin.json`, tab **Measure**, chip "ds10 core" (it sits in the chip row's **More** menu
because the Measure row is full), one menu entry per step, parameter dialog generated from the manifest (Settings...). Steps: image info, AB zero point
(▶ of the chip), maps, region statistics, regions -> mask FITS, convert regions, pixel table at (X, Y), forced photometry (band list typed in the dialog),
run script, **open offline bundle**. Results open in the text window; every run is recorded by the session recorder (`ds10core.*`) with its exact argument vector.
**Not done:** no interactive region editing/drawing beyond ds9's own (regions are used from `.reg` files - ds9's own Region menu saves them), no table view of the forced
catalogue inside the catalogue panel (the files are written and the path is shown), no band picking from the loaded frames, no mask hand-over to the shared
mask manager, no GUI progress strip for long runs beyond the job runner's.
GUI check: `scripts/verify_ds10core.sh` (real `bin/ds9`, 29 checks); screenshots `docs/shots/standalone_*.png`:
`standalone_01_chip_menu_measure_tab`, `standalone_01b_plugin_menu`, `standalone_02_calibration_result`, `standalone_03_region_statistics`,
`standalone_04_pixel_table`, `standalone_05_forced_photometry_result`.

## Continue a web session here (offline bundle)
When the network to the ds10 web server drops, a bundle made with the web app's **Offline export** button (or kept up to date by its automatic refresh) restores the
session on this computer; the full workflow, limits and the way back are in ds10-web `docs/offline_workflow.md`.
* Command line: `ds10 open session-offline.zip --dest ~/work/field` (checksums and the event hash chain are verified first; an existing non-empty folder is only replaced
  with `--force`). It creates `files/` (the images: **originals** for pro and above, float32 tool copies for free, as the bundle says), `regions/regions.reg`,
  `session_script.json`, `results/<job>/` (stored outputs), map/calibration sidecars `<image>.ds10maps.json`, `workspace.json`. `ds10 run field/session_script.json`
  then **reuses the stored result of every step whose text and inputs did not change** (status `restored`) and runs the steps you add at the end; each run is written to
  `runs/rNNNN/`. `ds10 status field` lists the runs.
* GUI: Measure tab, chip **ds10 core** menu, **Open offline bundle (continue a web session here)**. Parameters *Offline bundle* and *Folder to restore into* (empty =
  `~/ds10-offline/<bundle name>`). After restoring, the first image of the session is loaded in ds9 with the regions of the session (`regions/regions.reg`), the
  summary opens in the text window and the status line names the workspace; then **Run ds10-script** with `<workspace>/session_script.json` and *Folder with the images* =
  `<workspace>/files` continues the work. Screenshots: `docs/shots/offline_open_bundle_ds9.png`, `offline_open_bundle_script_run.png` (1440x900).
* Paths with spaces: the *Band images* list of the forced-photometry dialog is split at `|` or new lines when the text contains one (names may then contain spaces), else at spaces;
  a name with spaces can be put in double quotes.
* Limits: the bundle is a snapshot from the moment it was made; for multi-extension originals only the first image HDU is restored; bit-identical new results need the same
  Python/numpy/scipy versions as the server (the bundle records them).

## Parity with the web version
Same functions, same results. Test list (ds10-web `server/tests/test_shared_core.py`, `test_core_independent.py`; OGFinder `plugins/ds10core/tests`) and the
HUDF F105W/F125W/F160W web-versus-stand-alone batch (script exported from a web session, run on both): see ds10-web `docs/shared_core.md`. Summary: all
catalogues, detection images, segmentation, zero points, region statistics/masks and pixel tables have the **same sha256**; run-time fields (`seconds`) are
the only differences (identical Python/numpy/scipy versions).

## Checks
`scripts/run_all_checks.sh --only ds10core_tests,ds10core_gui`; `tools/validate_manifests.py` (0 problems).
