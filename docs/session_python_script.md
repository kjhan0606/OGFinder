# Session recorder and "Save Session as Python Script..."

The catalog panel records every computational step you run in the GUI and can export
them as one stand-alone Python 3 script (standard library only). The script is an
**automated pipeline**: run on new FITS files it repeats every *automatic* step and
skips the steps you did by hand. In `--mode replay` it repeats everything on the
original files and reproduces the GUI result.

Menu (Analysis menu of the catalog panel, below the separator):

| label | action |
|---|---|
| `Save Session as Python Script...` | write `ogfinder_session.py` |
| `Show Session Log` | list the recorded steps (class, time, status) |
| `Reset Session Log` | forget all recorded steps (start a new session) |

The existing ICL/LSBG "Export CLI Script" (bash) is unchanged and independent.

## What is recorded

Every step stores the exact argv it ran, its parameters, its class, the wall time and
a CRC32 of the tool's stdout / of the catalog after the step (used by replay checks).

| class | steps |
|---|---|
| **automatic** (run in pipeline mode) | Extract Sources; forced photometry in bands (`bands.measure`); Auto Mask, Grow/Shrink/Invert/Clear/Masked/Export mask; Copy mask to bands; catalog Trim, Sort, Save, Export FITS; Sersic, morphometry, segmap, PSF photometry, crowded-field photometry, multiband, dual-image extraction, completeness, bulge+disk, photo-z, SED fit, galaxy morphology, star classification, find PSF stars, build PSF (also extended / WebbPSF / TinyTim), deconvolution; ICL background / profile / measure / multi-threshold / decompose / colour profile; LSBG steps; cross-match (needs `--allow-network`); external AI service runs on the whole catalogue (`ai.run`, network services need `--allow-network`, see `docs/ai_services.md`) |
| **config** (applied silently) | Register band, Detection band, Remove band. In pipeline mode the band registry is built from the supplied images instead. |
| **manual** (skipped in pipeline mode, listed with a message) | hand-drawn mask Add/Erase regions (image pixel coordinates), mask Undo/Redo/Import, catalog Merge sources, Delete selected, Separate selected, Add object at position, Load catalog, Clear catalog, "show visible only" subsets (`catalog.unrecorded`), AI-merge accept/reject (recorded as manual merges), external AI service run on a hand-selected row subset (`ai.run` with `--numbers`) |

Deliberately **not** recorded: view actions (zoom, pan, colormap, scale, frames, hover,
selection),
Mask Stats, Show Overlay, tile view, dialogs that only display data. The AI-merge
*detection* call and `psf.load` (Load PSF) are not recorded either.

Automatic substitutes used in pipeline mode for image-specific choices:

* ICL centre (`--center`) = the brightest (smallest `MAG_AUTO`) source of the extraction catalogue of the new field.
* PSF stars (`--star-indices`) = the sources chosen by the automatic star finder step in the same run.
* Band registry: band name from the FITS `FILTER` keyword (or `NAME=path`), AB zeropoint from `ZP_AB`/`PHOTFLAM`+`PHOTPLAM`
  (via `ds9_sextract --info`), pixel scale from `PIXSCALE`/WCS; overridable with `--zp`, `--pixel-scale`, `--fwhm`.
  If the header has no zeropoint the value recorded in the GUI session is used and a warning is logged.
* Catalog Trim limits recorded in the GUI are reused (override with `--trim COL=MIN:MAX`); if they remove every source on the new data the catalog is kept unchanged and a warning is logged.
* Steps that need two or more bands (`bands.measure`, `mask.copy_to_bands`, colour steps) are skipped with a message when fewer bands are supplied.

The header of the generated script lists all recorded steps and marks each AUTO / CONFIG / MANUAL;
`ogfinder_session.py --list-steps` prints the same table.

## Using the script

    # one image
    python3 ogfinder_session.py image.fits --outdir out
    # several bands of one field (band names from FITS FILTER or NAME=path)
    python3 ogfinder_session.py --bands F105W=a.fits F125W=b.fits F160W=c.fits --detection-band F160W --outdir out
    # many fields, 4 in parallel, resumable; one field per line:  NAME path1 [path2 ...]  or  NAME band=path ...
    python3 ogfinder_session.py --input-list fields.txt --jobs 4 --resume --outdir out
    python3 ogfinder_session.py --field A:F105W=a1.fits,F160W=a2.fits --field B:b.fits --outdir out
    # exact replay of the GUI session on the original files
    python3 ogfinder_session.py --mode replay --outdir out

Common options: `--mode pipeline|replay` (default pipeline), `--include-manual`
(pipeline mode plus the manual steps; only meaningful on the original data),
`--python PATH` (interpreter with numpy/scipy/astropy/sep for the OGFinder tools; default
the current one), `--ogfinder-root DIR` (OGFinder install; recorded default), `--zp`,
`--pixel-scale`, `--fwhm` (`BAND=value` or one value for all), `--trim COL=MIN:MAX`,
`--override STEP:--flag=VALUE` (STEP is a glob, e.g. `mask.auto:--detect-thresh=4.0`,
`mask.*:--n-workers=1`), `--allow-network`, `--dry-run`, `--list-steps`.

JSON config (`--config cfg.json`; standard library `json`):

    {"options": {"jobs": 2},
     "fields": {"m51":  {"image": "m51.fits", "zp": {"*": "25.0"}},
                "hudf3": {"bands": {"F105W": "a.fits", "F160W": "c.fits"}, "detection_band": "F160W",
                          "zp": {"F105W": "26.27"}, "pixel_scale": {"*": "0.06"}, "fwhm": {"F160W": "0.18"},
                          "overrides": [["mask.auto", "--detect-thresh", "4.0"]]}}}

Outputs: `<outdir>/<field>/{work/, outputs/, manifest.json, catalog_final.tsv, session.log, tool_stderr.log}` and
`<outdir>/summary.tsv` + `summary.json` (objects, columns, masked fraction, steps ran/cached/skipped, seconds,
failed step, message). Exit code is non-zero if any field failed. No files are created in `~/.ds9` by the script itself (the verification checks this).

`--resume`: every step stores a hash of its argv, inputs (FITS sha256) and parameters; steps whose hash is
unchanged are skipped. Because mask steps edit files in place, a changed step makes the field restart from the
beginning (no partial reuse).

`manifest.json` contains: inputs (path, size, sha256, FITS data sha256, header values), parameters, tool versions
(OGFinder git head, `ds9_sextract` sha256 and source blob, numpy/scipy/astropy/sep versions), per-step status, seconds,
output file sha256 (undo/redo stacks and logs are marked `volatile`) and the final catalog hash.

## What "identical" means (replay mode)

Replay on the original inputs with the recorded literal values; the script itself checks the CRC32 of each
tool's stdout and of the catalog after every step against the GUI values. The verification harness additionally
compares (see below) catalogs byte-for-byte and column-by-column with tolerance 0, and FITS products by pixel
array bit equality and header equality. Absolute paths in tool stdout are normalised before the CRC.

## Verification

    scripts/verify_session_replay.sh        # needs Xvfb, bin/ds9, OGFINDER_PYTHON with numpy/astropy/sep

It runs a scripted GUI session (`scripts/verify_gui_session.tcl`, same procs the menus call) and exports the script, then
T1 HUDF 3-band replay vs. GUI, T2 m51 replay incl. Sersic, morphometry and ICL chain, T3 determinism (n-workers 1/3/default),
T4 pipeline mode on new data (3-band, 2-band and 1-band HUDF cut-outs and m51 in one batch with `--jobs 3`, `--resume`
all-cached, `--resume` with a changed parameter). `scripts/verify_icl_export.tcl` is the ICL bash-export smoke test.

## Known limits

* Tcl `exec` is wrapped globally while the recorder is active; catalogs are per frame in the GUI but the recorder is global
  (one catalog state per session).
* Steps that need PyTorch (photo-z, SED fit, galaxy morphology, star classify, AI merge) are recorded but were not
  executable in the verification environment.
* LSBG steps are recorded but not covered by the harness; the ICL chain is verified on m51 only.
* A mask file left in `~/.ds9` from an earlier session makes the GUI start from it; the script starts from an empty mask (a note is recorded).
* Manual catalog edits that cannot be expressed as parameters (`catalog.unrecorded`) are skipped even in replay (with a warning),
  so such sessions cannot be reproduced exactly.
* Results are bit-identical only with the same tool binaries / package versions; versions are written to the manifest.
