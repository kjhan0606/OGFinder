# HSC-I BCG+diffuse pilot

The approved 2026-10-12 workflow measures calibrated native-pixel annular
profiles on a mass-independent kpc grid. It retains the BCG and masks all
other detections, defects and reviewed bright-star regions. The output is
BCG+diffuse light. A separate area-overlap image supplies a radial diagnostic.

From the OGFinder repository, use
`analysis/icl_wgl_membership/hsc_pilot.py prepare` to create a new input
manifest, `inspect` to check explicit HDUs, global-sky provenance, calibration
and WCS, then submit `analysis/icl_wgl_membership/pilot_job.sh` for measurement.
Missing inputs stop inspection. Measurement requires a Slurm allocation.

The full input schema, commands, output columns and limitations are in
[pilot_README.md](../analysis/icl_wgl_membership/pilot_README.md).
The sequence and acceptance criteria are in
[implementation_plan.md](../analysis/icl_wgl_membership/implementation_plan.md).

Authentication is tested with
`python3 analysis/icl_wgl_membership/hsc_debug/probe_login.py` in an interactive
terminal. Enter the password only at its hidden prompt. A successful default
probe authenticates a `coadd` request; it does not classify that product as
global-sky. `--global-sky` retains the authenticated official manual and all
HDU headers of a single `coadd/bg` probe for review. The reviewed manual maps
that type to `deepCoadd` immediately after coaddition and before aggressive
background subtraction. The reviewed header identifies IMAGE/MASK/VARIANCE in
HDUs 1/2/3 and inherited primary `FLUXMAG0`, equivalent to AB zero point 27.0.
The source mapping passes the pilot product gate; this header-only probe does
not test sky residuals or retain image pixels. Details are in the GPFS review
JSON. The remaining twelve frames, PSF, masks, random fields and
injection/recovery are pending.

Native profiles retain negative linear flux and the original sky-subtracted
science plane. SEP background subtraction is used for source detection only.
Matched random profiles enter the corrected profile with the same cluster
weights; spatial delete-one jackknife removes both terms together. Without
sufficient regions there is no covariance estimate. PSF, sky transfer and
realistic injection/recovery remain pending; all pilot products retain
`science_measurement=false` and use no mass.

The existing scaled menu remains a separate path. Its survey metadata
checklist is not scientific acceptance, and its profile TSV records radial
participation rather than a calibrated surface-brightness profile. Supplied
same-lens covariance is recorded, but not propagated through a joint fit.
