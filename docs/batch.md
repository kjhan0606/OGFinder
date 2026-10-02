# Batch processing (plugin `batch`, tab Measure, chip "Batch")

Runs a **recipe** of plugin steps over many fields with parallel workers, resume, per-field logs and a summary table. Engine `ogfkit/batchrun.py`, headless template expansion `ogfkit/cliexpand.py`
(the Python twin of `::ogf::step::build_argv`), CLI `plugins/batch/batch.py`, GUI: *Run Batch…*, *Progress…*, *Summary Table…*, *Write Example Recipe*.
It complements the exported session script (`--input-list/--jobs/--resume`, docs/session_python_script.md), which replays recorded GUI sessions; the batch plugin instead drives the **manifest cli templates**, so every plugin step with a `cli` block (33 at present) can be applied to hundreds of fields without the GUI.

## Inputs
* `fields.txt`: `NAME path/to/image.fits` per line (`#` comments).
* `recipe.json`:
```
{"detect":  {"args": ["--detect-thresh", "3.0"]},          # bin/ds9_sextract <image> + args  (or  "catalog": "/data/{name}.tsv")
 "steps": [{"plugin": "noisemodel", "step": "model", "params": {"bw": 32}},
           {"plugin": "cluster", "step": "density"}]}
```
Parameters not given take the manifest defaults. The recipe is validated before anything runs (unknown plugin/step/parameter, step without cli template → error).

## Run
```
python3 plugins/batch/batch.py --fields fields.txt --recipe recipe.json --outdir out --jobs 6 --resume [--retries 1] [--timeout 600] [--only A,B] [--dry-run]
python3 plugins/batch/batch.py --example-recipe recipe.json
python3 plugins/batch/batch.py --expand xmatch.match --params-json '{"radius-arcsec": 2}'    # print the expanded argv
```
Per field `out/<NAME>/`: `catalog.tsv` (detection catalog + all `add_columns` outputs merged by NUMBER; a same-named column is replaced), `state.json` (per step: status, key, seconds), `field.log` (command, return code, stderr tail of every step),
`work/` (step work files, `snap_<i>.tsv` catalog after step i). `out/summary.tsv|json`: field, status, objects, columns, steps ok / cached / failed, failed step, seconds, message. `out/status.json` is updated after every step (the progress window reads it).
Exit status 1 when a field failed; a failing field never stops the others (its later steps are skipped and logged).

## Resume
Each step has a key = SHA-256 of (previous step's key + catalog hash, image checksum, the exact argv). With `--resume` a step is skipped when its key is unchanged and its snapshot exists; the first changed step and everything after it are recomputed.
Changing the detection arguments recomputes the whole field; changing a late step's parameter recomputes only that step and its successors. Interrupted runs (Cancel / SIGTERM / crash) continue where each field stopped.

## GUI
*Run Batch…* starts the runner in the background (the GUI stays usable), opens the progress window (one row per field: state, current step, steps done; progress bar; live runner output; Cancel, Resume, Summary table). Status is polled from `status.json` every 0.7 s.

## Validation (tests in `plugins/batch/tests`, measured)
* 33 cli templates of all plugins expand headlessly; in the GUI check the Python expansion is argument-for-argument identical to the Tcl `build_argv` for 8 steps with the live parameter stores.
* 6 synthetic fields × (detect + noisemodel + hello): `--jobs 1` 8.2 s, `--jobs 3` 2.8 s (×2.9); catalogs byte-identical between serial and parallel runs.
* Resume: second run 18/18 steps cached; changing noisemodel `bw` → 12 steps recomputed (2 per field), detection cached; changing the detection threshold → all 18 recomputed.
* Failure isolation: one missing image fails only its field (`detect: image not found`), the others finish; step timeout → rc 124 failure; recipe validation errors; dry run.
* Scale: 60 catalog-only fields × 1 step: 1.7 s serial (35 fields/s), 0.4 s with 6 workers (143 fields/s). This measures the orchestration overhead only (process start per step), not science throughput.
* GUI: 3-field run through the dialog (progress rows ok/ok/ok, summary, logs), resume all cached, summary window.

## Limits
* Only steps with a manifest `cli` template run (33 of 195 steps, see `tools/validate_manifests.py`); display-only and legacy Tcl-dialog steps (extract dialog, photometry/photo-z dialogs, deconvolution, galaxy morphology, masks, ICL/LSBG, Moving objects) are not available; for those use the exported session script. Steps needing a PSF/mask file take them from their parameters.
* Detection uses `bin/ds9_sextract` with the arguments you give (no per-field zero point / pixel scale lookup — put them in `args`, or use a catalog per field).
* One machine: workers are processes on this host; there is no distributed queue, no priority scheduling, and Cancel stops the runner, the step subprocesses of the fields in flight finish on their own.
* Steps with network access are run as is (the recipe must carry the `allow-network` parameter explicitly).
* Resume keys use the first MB + size of the image, not a full hash, and ignore changes of installed package versions.
