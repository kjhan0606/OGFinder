# Position-dependent HSC PSF injection assessment — 2026-10-19

## Purpose and approved inputs

This run closes one missing transfer test in the joint ICL–satellite injection
pipeline: it convolves synthetic cluster-centered diffuse light with the
cluster-center PSF and each synthetic satellite envelope with its own local
PSF. The signed HSC model values, including negative interpolation pixels, are
preserved. No HSC query or download was made in this calculation.

The run uses the user-approved 744-position PDR3 PSF batch plus the previously
approved cluster-10 center test. It covers 57 cluster/filter frames and 688
satellite positions, or 745 PSF inputs including the test center. Six other
jointly supported candidate frames have no PSF position in that approval and
were excluded. The PSFs are `HSC-I`, `pdr3_wide`, `coadd` model cores; each
stamp is about 7 arcsec across. The calibration inputs and all 744 file hashes
are recorded in the PSF audit linked from
[`psf_picker_batch_2026-10-19.md`](psf_picker_batch_2026-10-19.md).

## Execution and independent audit

The versioned runner is [`run_joint_icl_satellite_psf_redetection.sbatch`](run_joint_icl_satellite_psf_redetection.sbatch),
which calls [`run_joint_icl_satellite_batch.py`](run_joint_icl_satellite_batch.py),
[`joint_icl_satellite_injection.py`](joint_icl_satellite_injection.py), and
[`position_psf.py`](position_psf.py). Slurm job **428480** completed with exit
code 0 in 4 min 54 sec; peak RSS was 720,972 KiB. The calculation uses CPU
code. SEP detection and K2.5 masks were regenerated for baseline, ICL-only,
satellite-only, and joint-injection scenarios, using the explicit conservative
`mask_all` policy when a unique BCG candidate was unavailable. The six frames
without an approved PSF mapping were not run. No original FITS image was
modified.

The frozen-code output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_2026-10-19_02/`.
All 57 selected frames completed, with 13 `i` and 44 `i2` frames, and zero
errors. The independent audit is `psf_redetection_audit.json`; it verified all
57 result hashes, the exact approved center/satellite position mapping, all
745 PSF uses, and the run-time source hashes. Audit SHA-256:
`ce377060cfa2143936b36a400f5bfdb33fb355956bf3513d7a2fa4964595ac1e`.
The batch summary SHA-256 is
`df53035c6066bb4215c0d7d2aec9fef423e62a52bb6c60b89eb5be0e9d17ccc1`.

The first complete calculation in `_01` is preserved, but not used as the
authoritative run: its source file changed during execution, producing mixed
metadata hashes and a stale batch schema label. The `_02` run was made after
the code was committed and frozen, and is the audited result reported below.

## Paired PSF-core sensitivity

Each row is the per-frame transfer response with PSF convolution minus the
matching no-PSF redetection run on the same 57 frames. The bracket gives the
16th–84th percentiles across frames, not a confidence interval. Responses are
output surface brightness per unit injected amplitude.

| Transfer response change | i (n=13) | i2 (n=44) |
|---|---:|---:|
| Diffuse injection → ICL 200–300 kpc | +0.000032 [−0.000024, +0.000063] | +0.000040 [+0.000012, +0.000066] |
| Satellite envelopes → ICL 200–300 kpc | +0.000307 [−0.028870, +0.027425] | +0.003651 [−0.006927, +0.012599] |
| Diffuse injection → satellite 50–75 kpc | −0.000006 [−0.000240, +0.000068] | −0.000015 [−0.000102, +0.000165] |
| Satellite envelopes → satellite 50–75 kpc | −0.003545 [−0.035981, +0.041034] | +0.008016 [−0.036573, +0.034880] |

The table values are paired descriptive differences, not a detection test. The
median PSF-core effect is small in these fixed injection shapes, while the
frame-to-frame response is broader and includes zero in every displayed
16th–84th percentile interval. This characterizes compact-core sensitivity
for the chosen annuli; it does not show that PSF systematics are negligible.

Joint injection and separate injections are not exactly additive after mask
re-detection. Median `both − (ICL-only + satellite-only)` residuals are
`+0.0334` nJy/arcsec² (i) and `−0.00121` nJy/arcsec² (i2) for the ICL metric;
for the satellite metric they are `−0.00660` and `+0.00452` nJy/arcsec².
This nonlinearity is part of the detection/masking response, not an estimate
of the observed ICL–satellite flux covariance.

## Limits and next analysis gate

This is an injection/mask sensitivity on the existing selected support sample,
not an observational detection, population completeness estimate, or
calibration correction. It omits the extended PSF wings, injected-source
photon noise, correlated coadd noise, HSC upstream sky subtraction, extinction
and passband uncertainties. The compact PSF models do not validate light on
the much larger angular scales of ICL wings. `science_measurement=false`,
`observational_validation_passed=false`, and
`flux_covariance_available=false` remain unchanged.

Next, add extended-PSF-wing and noise-transfer sensitivity only after their
inputs are calibrated and matched to these frames. Then assess observed-pixel
injection recovery and spatial covariance; keep the i and i2 filters separate
and retain the full mask/redetection failure accounting. The current result
does not authorize further HSC data requests.
