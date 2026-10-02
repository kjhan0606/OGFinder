# Reproducibility bundle (plugin `repro`, tab Measure, chip "Repro")

A zip that records **what was run, with which software, on which inputs, giving which outputs**, and a `verify` command that re-checks and re-runs it. Engine `ogfkit/repro.py`, CLI `plugins/repro/repro.py`.
GUI: *Create Bundle*, *Verify Bundle…*, *Bundle Info…*.

## Bundle contents
`manifest.json` (schema 1: creation time local + UTC, OGFinder git HEAD / branch / dirty flag / list of modified files / SHA-256 of the local diff, platform, Python, versions of the key packages
numpy scipy astropy sep matplotlib scikit-learn photutils eazy bagpipes torch Pillow numba, SHA-256 of `bin/ds9_sextract` and `bin/ds9`, **input files** (path, size, SHA-256), **outputs** (SHA-256, rows, columns), kind, verify command),
`params.json` (parameter store of every plugin), `steps.json` (the recorded session steps with argv, parameters and stdout CRCs, or the per-field step status of a batch run), `requirements.lock` (`name==version` of every installed distribution, from `importlib.metadata`, no network),
`catalog_final.tsv` (or one catalog per field), `README.txt`, `verify.sh`, and by kind: `session.py` (the exported Python session script, kind `session`) or `recipe.json` + `fields.txt` (kind `batch`); `inputs/` only with *Include copies of the input images*.

## Commands
```
python3 plugins/repro/repro.py bundle --out B.zip --kind session --catalog C.tsv --inputs a.fits --params P.json --steps S.json --session-script S.py [--include-inputs] [--note TXT]
python3 plugins/repro/repro.py bundle-batch --out B.zip --batch-out out --fields fields.txt --recipe recipe.json
python3 plugins/repro/repro.py verify B.zip [--inputs-dir DIR] [--rtol 0] [--strict-env] [--no-rerun] [--jobs N] [--json R.json] [-- extra args for the session script, e.g. the image path]
python3 plugins/repro/repro.py info B.zip
```
`verify` checks, in order: bundle readable / schema; **environment** (package and Python versions vs the manifest: WARN, FAIL with `--strict-env`); **code version** (git HEAD and local diff vs the manifest: WARN); **binaries** (SHA-256);
**inputs** (SHA-256 of the file at the recorded path, in `--inputs-dir` by name, or the bundled copy: FAIL when missing or different, and the re-run is then skipped); **bundled outputs** against the manifest checksums (detects an edited bundle);
**re-run** (kind `batch`: the recipe through `ogfkit.batchrun`; kind `session`: `session.py --mode replay`; kind `catalog`: none) and **compare** every output with the re-run output
(`identical` = byte-identical; `numeric` = all differences within `--rtol/--atol`; otherwise FAIL with the number of differing cells, the columns and the largest relative difference). Exit status 1 on any FAIL.
In the GUI the re-run is optional (*Re-run the processing when verifying*), because a session replay takes minutes.

## Validation (tests in `plugins/repro/tests`, measured)
* Round trip, 3 synthetic fields × (detect + noisemodel + hello), batch bundle (~140 KB without inputs): `verify` → 15 checks, 0 failures; the 3 re-run catalogs are **byte-identical** to the bundled ones.
* Tampered input (one pixel changed by 1.0): FAIL `input img1.fits`, re-run skipped; the same input restored in another directory (`--inputs-dir`) verifies.
* Bundled catalog edited after the fact: FAIL (checksum vs manifest). Recipe edited (detection threshold 3 → 5): FAIL on all 3 catalogs (≈ 590 cells in 38 columns differ, max relative difference 1.7–2.0).
* Table comparison: identical / within tolerance (1e-7 difference accepted at rtol 1e-6, rejected at 0) / text column different / different column lists. A changed numpy version in the manifest is reported as WARN (FAIL with `--strict-env`).
* Bundle with `--include-inputs`: verifies after the original image was deleted.
* GUI (real m51 session of the check harness): bundle written with session script, parameters of ≥20 plugins and ≥10 recorded steps; integrity verification; the harness additionally re-runs the session bundle (`verify` with replay) — see docs/testing.md.

## Limits
* Bit-identical re-runs need the same code, binaries and package versions; the bundle records them and `verify` warns about differences but cannot reinstall them (use `requirements.lock` to rebuild an environment).
* The environment lock is the list of installed Python distributions of the interpreter that created the bundle; it does not capture system libraries, the C compiler, GPU/torch determinism or the OS beyond the platform string.
* Network-dependent steps (TAP queries, reference fetches) are recorded but re-run live, so their results can legitimately differ; the data returned are not stored.
* Session-kind verification re-runs the exported script, which skips manual steps (hand-drawn masks etc.); sessions containing such steps cannot be reproduced exactly (see docs/session_python_script.md).
* Input checksums cover whole files; FITS files whose data are unchanged but header/metadata differ count as different.
