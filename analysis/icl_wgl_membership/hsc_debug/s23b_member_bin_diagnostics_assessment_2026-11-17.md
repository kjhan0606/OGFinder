# S23B member-profile bin diagnostics — 2026-11-17

## Pipeline change

`measure_s23b_camira_member_outer_profiles.py` now emits a diagnostic row for
all seven radial bins of every selected candidate, including unmatched and
failed candidates. Retained profile columns and the three-item `process_frame`
return interface are preserved. The optional last argument collects bin
records. The measurement summary schema is v2; existing v1 products are left
intact. A new diagnostic table records:

- Target profile availability, valid pixel count, area fraction and raw flux.
- Number of accepted rotated centers, available profiles and per-bin controls
  passing the coverage floor.
- Separate target-area and control-support pass flags, retention and failure
  reason(s). Overlapping failures are joined with `|`.
- Raw per-control coverage, pixel count and flux for subsequent offline
  threshold comparisons. Zero-pixel/nonfinite flux is JSON null.

The mutually distinguished control failures are too few accepted centers,
too few profiles available, or too few controls with adequate per-bin area.
Upstream rotation rejections (neighbor/edge/member proximity) are not
individually classified. Unmatched candidates have unmeasured support fields
rather than artificial failed measurements. No photometric selection,
masking, flux estimator or threshold default was changed. Coverage floors
must be finite and within (0,1]; the requested control count must lie within
the implemented rotation limit. New summaries include output/dependency hashes.

## Verification

Six new unit tests cover exact threshold boundaries, bins independent of
other bins, overlapping target/control failures, absent profiles, null
unmatched values and zero-pixel JSON. Nine existing shared-profile regression
tests also pass. `git diff --check` passed.

A real CPU pilot remeasured the already approved i image `S23B_10882` and i2
image `S23B_7961`: 24 candidates and 168 diagnostic rows, with 19 retained
candidate/bin measurements. The automated audit proves that all retained
profile rows and the original candidate diagnostics are **identical** to the
corresponding historical full-sample subset. The raw-control reconstruction
of all retained target-minus-control fluxes has maximum difference 0.0.
The pilot's image-analysis elapsed time was 3.75 seconds; it does not predict
full-run performance for denser or larger images.

At 50–75 kpc the pilot's failure combinations are:

| Combination | i candidates | i2 candidates |
| --- | ---: | ---: |
| No unique SEP match | 4 | 4 |
| Target area alone fails | 1 | 6 |
| Control-bin area alone fails | 4 | 1 |
| Both target and control-bin areas fail | 2 | 0 |
| Retained | 0 | 2 |

These two frames are an implementation check, not a population estimate.
Offline coverage floors 0.4/0.6/0.8 retain 5/0/0 i candidates and 5/2/0 i2
candidates at this radius (>=3 controls). No threshold is recommended from
these counts; a changing candidate sample cannot establish a flux trend.

Pilot measurement directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_pilot_2026-11-17_01/`

Pilot measurement summary SHA-256:
`9c1e5b28ed7c456269e62e95137413fdeb46e0c8a2306eeb9ee6a2b01487e1de`.

Pilot audit directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_pilot_audit_2026-11-17_01/`

Pilot audit summary SHA-256:
`27b2f75d0cc352ffe17098b44b6166729f6ca7938f514309c798f2b3a2cd0e13`.

## Full-sample run

After the pilot passed, the 979 already selected target frames were started
on LagEunha as one serial CPU run in a persistent execution session. No Slurm, GPU, catalog query,
authentication or data transfer is involved. The runner shell was verified
live as PID 3746049 (image child PID 3746134; exec session 30421) with a child image-measurement process and advancing frame
logs. This is a launch record; do not infer completion from it or from a
status file alone. Recheck the PID/command and log, and inspect `audit/summary.json`
before claiming the full comparison passed.

The first detached attempt (`...full_2026-11-17_01/`, PID 3742503) stopped
after four frame-start logs: neither its runner nor its image child remained
alive. Its stale RUNNING status was corrected to INTERRUPTED and all partial
files were preserved. The cause was not established. No prior process was
alive when the new run was started. Do not use detached shell launch success
as evidence of continuing execution.

Historical second-attempt run root (subsequently interrupted):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_02/`

The root holds `runner.pid`, `status`, `provenance.txt` and `run.log`.
On success it gains `measurement/` and `audit/`. The wrapper records COMPLETE
only after the measurement and exact historical comparison both succeed.
Existing output paths are refused; an interrupted run may leave a partial
measurement directory and must be inspected before any new run is launched.

From the repository, a future independently named rerun uses:

```bash
bash analysis/icl_wgl_membership/hsc_debug/run_s23b_member_bin_diagnostics.sh \
  /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/NEW_UNIQUE_RUN_DIRECTORY
```

The offline audit can be run separately against completed v2 products:

```bash
/home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_s23b_member_bin_diagnostics.py \
  --measurement /ABSOLUTE/COMPLETED/measurement \
  --output /ABSOLUTE/NEW/audit
```

## Scientific status

All products retain `science_measurement=false`. The full failure census is
needed before selecting common candidate/cluster support for threshold
sensitivity. The new diagnostics alone do not solve missing-member
completeness, membership purity, PSF wings, source Poisson, sky transfer or
coadd covariance. The S16A 570-center WL random count plan remains a separate
approval-gated action and was not submitted by this work.


## Checkpoint/restart correction and current full run

Persistent exec session 30421 subsequently ended with exit 143. Runner and
image child were both confirmed absent; `_02` had reached frame-start log
662/979 but produced no completed measurement. Its status was corrected to
INTERRUPTED and all partial files retained. The exact termination cause is
unconfirmed. Neither earlier attempt is a completed dataset.

The measurement now saves a JSON checkpoint after each completed image.
Each checkpoint is bound to the input hashes, selected sample, thresholds,
measurement/helper/dependency source hashes, Python and NumPy versions. A
payload checksum catches corruption; an exclusive file lock blocks concurrent
writers. Files are written through a temporary file, fsynced and atomically
renamed. Explicit `--resume` is required for existing checkpoints, and
changed configurations are rejected. `--max-new-frames 100` stops cleanly
without publishing an incomplete final measurement. A resumed invocation
loads completed frame payloads and computes only missing frames. A complete
measurement is published after all selected frames are available.

Three checkpoint tests passed, with ResourceWarning promoted to an error:
resume/configuration binding and no overwrite; corrupted payload rejection;
explicit resume and exclusive-writer protection. A real two-frame pilot was
stopped after one new frame, then resumed to finish. Its 19 retained profile
rows and candidate diagnostics match the historical subset exactly; the
168 new diagnostic rows reconstruct retained fluxes with zero difference.

Resume pilot:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_resume_pilot_2026-11-17_01/`

Resume-pilot audit:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_resume_pilot_audit_2026-11-17_01/`

The revised wrapper runs 100-frame invocations in tmux session
`astrafek_member_diag_20261117`. Its current root is:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_03/`.
The runner was confirmed live as PID 3767575, with an active image child and
completed `measurement.checkpoints/S23B_*_FILTER.json` files. Each batch has
a new image child PID; the tmux session and runner identify the whole run.
The full audit and scientific sensitivity remain unfinished at this record.

If execution is interrupted, first inspect tmux, runner and child state.
After confirming no writer remains, resume a bounded invocation from the repo:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/measure_s23b_camira_member_outer_profiles.py \
  --all-eligible --resume --max-new-frames 100 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_03/measurement
```

Repeat bounded invocations until the final measurement is published, then
run the existing offline historical audit into the root's new `audit/`
directory. Preserve any interrupted run's existing artifacts. Do not edit
measurement/helper code mid-run: the next resume rejects changed code hashes.
