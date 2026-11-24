# S23B member outer-light coverage sensitivity — methods, 2026-11-17

## Question and estimator

The original 50–75 proper-kpc member-light association retains only a small
subset of selected CAMIRA candidates. Changing coverage cuts changes both
which candidates contribute and the mean flux of accepted rotated controls.
The new `analyze_s23b_member_coverage_sensitivity.py` separates these effects
using the v2 per-bin diagnostics, after the historical reconstruction audit
has passed. It reads saved tables only, without image access or network IO.

Coverage floors are 0.4, 0.6 and 0.8 for both target and control annuli; at
least three rotated controls must pass. Each candidate flux is target mean
minus the mean of the passing controls. A cluster value is the equal-weight
mean over its retained candidates. Three sample definitions are evaluated:

1. **Native:** all candidates passing each floor and their supported clusters.
2. **Common clusters:** the same clusters at all floors, allowing their retained
   candidate lists to vary.
3. **Common candidates:** the same individual candidates at all floors; only
   the passing-control mean can change the candidate flux.

Membership supports must be nested. The 0.6 reconstructed candidate keys and
fluxes must match the saved outer profile table within 1e-9. Duplicate
candidate rows or changed input/product hashes stop the analysis.

The existing adjusted outer-light model is reused: outer light versus ICL
relative contrast with log richness, redshift, minimum ICL common-area
fraction, median target coverage and median control coverage as nuisance
covariates. Constant nuisance columns are omitted by the existing helper.
In both fixed-sample modes, the two coverage covariates are frozen at their
0.6 values so that changing the threshold cannot change the design matrix.
Native models use their own coverage covariates and therefore differ in both
sample and design. Interpret native and fixed-mode slopes accordingly.

The six established connected target/control sky-block layouts (2.5°, 5°,
10°, each with two grid phases) are evaluated with 5,000 block-bootstrap
draws per fit. A 5° phase-zero layout is reported as the primary descriptive
layout. In fixed modes, the paired outcome difference relative to floor 0.6
is fitted with the same design matrix and blocks. This produces an interval
for the slope change rather than subtracting independently estimated interval
endpoints. Unusable samples/rank/bootstrap support are recorded explicitly;
missing slope estimates are null, not zero.

## Verification so far

Four tests passed: nested support and recomputed control means; common
clusters versus common members; duplicate rejection; and a known linear
signal whose paired slope change is recovered with the shared spatial design.
A complete table-only pilot workflow was also executed using the audited two
real frames. Its threshold rows, cluster tables and paired-change tables were
written successfully. All pilot association estimates correctly report
insufficient sample; no scientific regression result is inferred from two
frames. The full 979-frame diagnostics were still processing when this method
was prepared, so no full-sample threshold association is reported here.

Pilot output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_coverage_sensitivity_pilot_2026-11-17_01/`.

## Full-run command after its audit completes

From `/home/kjhan/BACKUP/Astrafek/OGFinder`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/analyze_s23b_member_coverage_sensitivity.py \
  --run-root /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_03 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_coverage_sensitivity_full_2026-11-17_01
```

Explicit `--measurement` and `--audit-summary` paths support separately stored
pilot products. The runner refuses existing output paths. Outputs include
cluster-level threshold profiles, all layout fits, paired slope changes and
summary/input/output/script/helper hashes. Partial output is staged under a
new directory and published only on successful completion.

## Limits and decision rules

These are sensitivity estimates conditional on detection and upstream
selection. Even common candidates remain selected by the strictest coverage
cut. A nonzero interval does not establish physical association: these are
multiple marginal comparisons, with unresolved PSF, source-Poisson, sky,
coadd covariance and membership systematics. Native sample changes cannot
be interpreted as pure photometric response. Do not use threshold selection
to maximize apparent significance or convert retention to completeness.
Keep `science_measurement=false` for every mode and for all paired changes.


## Standalone figure

`plot_s23b_member_coverage_sensitivity.py --analysis COMPLETED_ANALYSIS --output NEW_FIGURE_DIRECTORY`
verifies the analysis output hashes and exports PNG/PDF panels showing the
number of supported clusters and the adjusted slope at each cut for i/i2.
It draws only complete 5° phase-zero fits, with marginal spatial-block
intervals, and labels unsupported panels explicitly. Panel titles/legend are
inside the panels; adjacent panels share row ranges and abut. The renderer
was executed and visually checked on the completed pilot, which correctly
shows no supported slope estimates. This rendering does not substitute for
the full sensitivity calculation or its six-layout robustness assessment.


The original non-checkpointed full attempts stopped before completion. The
current full measurement uses tmux session `astrafek_member_diag_20261117`,
100-frame invocations and per-frame checkpoints under the `_03` run root.
Only use that root after its historical audit completes. Previous `_01`/`_02`
partial directories are preserved and do not authorize this sensitivity run.

## Scheduled local follow-up for the current full run

`run_s23b_member_coverage_after_diagnostics.py` is running in tmux session
`astrafek_member_followup_20261117` (PID 3776762), waiting for the exact runner
PID 3767575 plus its kernel start-time identity. The current Python build
exposes no `os.pidfd_open`; the supported fallback checks `/proc/PID/stat`
start-time ticks every five seconds and does not follow a reused PID. Four
tests passed for live identity, PID reuse, missing processes and wrong-runner
rejection before creating work products.

The follow-up verifies COMPLETE status, all 979 frames, successful historical
comparison and the matching measurement-summary hash before executing the
analysis. Analysis/helper/plot source hashes are pinned at scheduling; changed
sources abort execution. Do not manually start a duplicate full sensitivity
calculation while this scheduled follow-up remains live.

Scheduled outputs under the current `_03` root:

- `coverage_followup/plan.json`: runner identity and pinned source hashes.
- `coverage_followup/status`: waiting, running, complete or failed.
- `coverage_followup/run.log`: downstream command output after prerequisites pass.
- `coverage_followup/analysis/`: full 5,000-draw sensitivity tables and summary.
- `coverage_followup/figures/`: verified PNG/PDF exports and provenance.

The waiting process and original diagnostic runner were both confirmed live;
the status was WAITING_FOR_DIAGNOSTIC_RUNNER and no premature analysis output
existed. That verification establishes a live wait, not completed science.
