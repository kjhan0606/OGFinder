# Run guided workflow steps from a Unix shell

The guided workflow catalogue and desktop plugin manifests define the mapping.
Inspect the mapping before choosing steps:

```sh
cd /path/to/OGFinder
python3 -m ogfkit.workflow_cli list
python3 -m ogfkit.workflow_cli show icl
```

`direct_cli` means the desktop manifest exposes a command that can be attempted
without GUI state; it does not imply that dependencies are installed or that
the result is scientifically valid. GUI hooks, catalog merges, manual steps,
missing CLI templates and network steps require separate support. The guided
recipes contain scientific review and validation stages that this runner does
not automatically carry out. Never treat a completed selected subset as a
complete or scientifically validated recipe.

Run one supported step in a **new** directory. For example, prepare a JSON
object containing the `core-dir` parameter if the shared core is not found by
the normal discovery mechanism, then run:

```sh
python3 -m ogfkit.workflow_cli run-step icl audit \
  --image /data/image.fits --params-json /data/core-params.json \
  --work-dir /data/output/audit
```

The optional ICL step **기존 HSC S23B 검증 재현** uses the saved S23B research
files. Set `s23b-evidence-root` in the desktop ICL settings or pass a JSON
object with that key through `--params-json`. The default points to the
research archive on GPFS. Both entry points call the same evidence builder and
record `s23b_icl_evidence/summary.json` and `assessment.md`. Use a fresh work
directory for a new replay. The web catalogue shows the desktop/shell step.

```sh
python3 -m ogfkit.workflow_cli run-step icl s23b_evidence --work-dir /data/output/icl-replay
```

For a sequence, write an explicit JSON plan. Paths in it are resolved relative
to the shell's working directory, so absolute paths are recommended:

```json
{
  "schema": "astrafek-workflow-plan-v1",
  "steps": [
    {"id": "first", "recipe": "icl", "step": "audit",
     "image": "/data/first.fits", "params": {"core-dir": "/path/to/ds10-web/server"}},
    {"id": "second", "recipe": "icl", "step": "audit",
     "image": "/data/second.fits", "params": {"core-dir": "/path/to/ds10-web/server"}}
  ]
}
```

```sh
python3 -m ogfkit.workflow_cli run-plan /data/plan.json --work-dir /data/output/batch
```

An input (`image` or `catalog`) may refer to a **declared output** of a prior
step using `{"step": "first", "output": "filename"}`. Use an output of the
appropriate format for the next tool; this syntax does not convert data. Before
execution, the runner checks the sequence, dependencies, and plugin parameter
names. It verifies the preceding output hash before passing it to another
step, stops at the first error, and never overwrites an existing work directory.
It writes the selected plan, plan hash, a plan result, and per-step receipts,
logs, exact command arguments, and input/output SHA-256 hashes. A result of
`completed` applies only to those selected commands and their declared files;
`science_validated` remains false. The existing GUI session exporter provides
a separate recorded-session replay and multi-field interface; see
[session_python_script.md](session_python_script.md).

## Masked ICL surface-brightness profile

`icl/profile` uses the same `ds9_icl.py` measurement as the desktop ICL
profile action. The shell adapter requires an image, a separately prepared mask
and explicit **zero-based** BCG coordinates (`center-x`, `center-y`) in the
ICL plugin parameter JSON. Nonzero mask pixels are excluded; the mask must be
finite, have the same pixel dimensions, and already be registered to the image.
The program verifies dimensions and finiteness, not celestial WCS alignment.
The desktop profile action also passes its shared mask if one exists and warns
in its status when no mask is available.

```sh
python3 -m ogfkit.workflow_cli run-step icl profile \
  --image /data/sky_checked.fits --mask /data/source_mask.fits \
  --params-json /data/icl-parameters.json --work-dir /data/output/profile
```

Example parameter JSON (values here are illustrative, not recommended
calibration values):

```json
{"center-x": 100, "center-y": 100, "rmin": 5, "rmax": 80,
 "nsteps": 20, "spacing": "linear", "pixel-scale": 0.168,
 "mag-zeropoint": 27.0}
```

The product is `icl_profile.tsv`. Its scatter-based errors assume independent
pixels and do not include correlated background noise, mask-selection effects,
PSF wings, residual sky or cluster-level covariance. A mask of the BCG or its
extended light can bias the profile. Establish these choices and validate the
input before any science interpretation. This single-image profile is distinct
from the research-specific HSC ICL and satellite injection/stacking pipeline.

When a desktop background run displays its result in a new frame, the app keeps
an in-memory association between that result, its original image, and the mask
used by the profile (including the iterative refined mask, when available).
This association is not a persistent project manifest: after a fresh app
session or an independently opened background-subtracted FITS, reconnect and
verify the corresponding mask and other inputs before measurement.

### Automatic source-mask to profile example

The `icl/mask` shell adapter runs the desktop's `ds9_mask.py --mode auto`
engine with the ICL preset and a new empty mask. It saves `mask.fits` (flags)
and `mask_bool.fits` (effective mask). Unlike the GUI path that can reuse an
existing mask and manual edits, this adapter starts fresh; its receipt records
that limitation. The mask may cover the BCG or genuine diffuse light: review
it and the input's sky history before scientific use.

The following plan connects the boolean output to the profile by recorded
hash; it does not perform sky modelling, PSF correction or ICL validation:

```json
{
  "schema": "astrafek-workflow-plan-v1",
  "steps": [
    {"id": "mask", "recipe": "icl", "step": "mask", "image": "/data/image.fits",
     "params": {"detect-thresh": 5, "mag-zeropoint": 27}},
    {"id": "profile", "recipe": "icl", "step": "profile", "image": "/data/image.fits",
     "mask": {"step": "mask", "output": "mask_bool.fits"},
     "params": {"center-x": 100, "center-y": 100, "pixel-scale": 0.168,
                "mag-zeropoint": 27, "rmin": 5, "rmax": 80, "nsteps": 20}}
  ]
}
```

All parameter values above are examples and must be chosen for the dataset.

### Explicit ICL background modelling

The `icl/sky` adapter runs the desktop `ds9_icl.py --mode background` engine
on the supplied image and mask. Set `bkg-method` explicitly (`polynomial`,
`chebyshev`, or `sep_large`); the runner does not assume that a local fit is
appropriate for your field. Its declared products are `icl_background.fits`,
`icl_bgsub.fits`, and `background_mask.fits`. The last file copies the input
mask in the standard fit or contains the expanded mask from an iterative fit.
A following profile may reference both products:

```json
{"id": "sky", "recipe": "icl", "step": "sky", "image": "/data/image.fits",
 "mask": {"step": "mask", "output": "mask_bool.fits"},
 "params": {"bkg-method": "polynomial", "bkg-order": 1}}
```

```json
{"id": "profile", "recipe": "icl", "step": "profile",
 "image": {"step": "sky", "output": "icl_bgsub.fits"},
 "mask": {"step": "sky", "output": "background_mask.fits"},
 "params": {"center-x": 100, "center-y": 100, "pixel-scale": 0.168}}
```

Put these objects after a `mask` step inside the plan's `steps` array, using
values appropriate to the data. This adapter fits and subtracts from the
supplied image directly; GUI background fitting may first interpolate masked
pixels. Although the same fitting engine and mask are used, byte-identical
GUI replay is therefore not promised. Local sky fits can remove diffuse ICL;
iterative source detection can also mask the BCG or the whole measurement
annulus. The runner records the warning and fails if the profile has no valid
pixels. Neither successful execution nor masked-pixel exclusion validates sky
transfer, PSF, WCS alignment, or covariance.

### Save and continue a plan

For plans created with the current runner, `--stop-after` deliberately stops
after a successful step so its products can be inspected:

```sh
python3 -m ogfkit.workflow_cli run-plan /data/plan.json \
  --work-dir /data/output/batch --stop-after mask
python3 -m ogfkit.workflow_cli run-plan /data/plan.json \
  --work-dir /data/output/batch --resume
```

Resume checks the exact saved plan, completed receipts, input and output file
hashes and sizes, plus a fingerprint of the workflow sources, known shared-core
sources, Python executable, installed Python distribution versions and selected
environment settings. It uses an exclusive run lock and preserves completed
steps rather than executing them again. A changed input, result, receipt, code
or recorded environment is refused. If a step failed or an interruption left
an incomplete step directory, preserve it, inspect its logs and start a new
work directory; do not run a second process against a still-running step.
Older plan results without the fingerprint cannot be resumed. To edit a
mask, save a reviewed copy separately and use it as input to a new plan.

The fingerprint is an integrity check over the recorded scope, not a complete
capture of arbitrary external programs, shared libraries, machine state or
scientific assumptions. A resumed `completed` plan remains
`science_validated=false` and its sky, mask, PSF and noise must still be
reviewed.

### Export a Python launcher

A validated explicit plan can be embedded into a Python script for Unix-shell
execution. The inputs must exist when the script is exported; the output script
is created without overwriting an existing file.

```sh
python3 -m ogfkit.workflow_cli export-python /data/plan.json --output /data/pipeline.py
python3 /data/pipeline.py --work-dir /data/output/run --stop-after mask
python3 /data/pipeline.py --work-dir /data/output/run --resume
```

The launcher embeds only the selected offline steps and parameters; it imports
Astr afek/OGFinder at runtime. If the installation has moved, pass
`--ogfinder-root /path/to/OGFinder` and ensure its shared core and Python
dependencies are installed. It uses the same pause, fingerprint and receipt
checks as `run-plan`. Omitted stages, manual review, user-specific inputs, sky
validation and science approval are not synthesized by the launcher. The
standalone science workflow window also offers **분석 세션 Python 저장** for
exporting the entire recorded GUI session via the existing session recorder;
its batch behavior and manual-step limits are described in
[session_python_script.md](session_python_script.md).
