# Position-dependent PSF and mask-redetection scale comparison — 2026-10-20

## Question and design

The fixed-mask local-plane audit showed that broader diffuse profiles lose more
signal after Astrafek's downstream plane estimator. This paired test asks
whether that scale dependence remains when the approved position-dependent
HSC-I coadd model-PSF cores are applied and SEP detections/K2.5 masks are
regenerated for every injection scenario.

Astrafek ran the same 57 cluster/filter frames at two diffuse exponential
scales, 100 and 500 kpc, with a surface brightness of 1 nJy/arcsec² at
R=100 kpc. Both runs use the same satellite-envelope injections (15 kpc
scale), local PSF mapping, `mask_all` policy for ambiguous BCG detections, and
runtime source files. The filters remain separate: i has 13 frames and i2 has
44. This is a same-frame injection comparison, not an observed target-minus-
control measurement.

For each frame the reported retention ratio is the ICL response in the
200–300 kpc annulus after local-plane correction divided by its raw-annulus
response. The table gives the median and 16th–84th percentiles across frames.
The final column is the same-frame 500-minus-100 kpc ratio difference. These
quantiles describe frame-to-frame spread; they are not confidence intervals.

| Filter | Frames | 100 kpc retention | 500 kpc retention | Paired difference |
|---|---:|---:|---:|---:|
| i | 13 | 1.002 [0.954, 1.066] | 0.644 [0.611, 0.679] | −0.341 [−0.431, −0.304] |
| i2 | 44 | 0.997 [0.957, 1.055] | 0.633 [0.588, 0.658] | −0.362 [−0.423, −0.316] |

The raw and plane-corrected ICL responses per unit input amplitude were:

| Filter | Scale | Raw response | Plane-corrected response |
|---|---:|---:|---:|
| i | 100 kpc | 0.215 [0.170, 0.232] | 0.206 [0.163, 0.239] |
| i | 500 kpc | 0.746 [0.682, 0.758] | 0.476 [0.411, 0.498] |
| i2 | 100 kpc | 0.224 [0.183, 0.252] | 0.227 [0.174, 0.256] |
| i2 | 500 kpc | 0.739 [0.710, 0.770] | 0.465 [0.423, 0.510] |

This confirms scale-dependent attenuation on the selected target frames after
PSF convolution and mask re-detection: the broad 500 kpc profile retains about
63–64% of its raw-annulus response, while the compact 100 kpc profile is near
unity in the median. The direction and size agree with the separate fixed-mask
74-pair plane-injection audit. They do not establish a universal correction
for observed ICL. In i, one frame at each scale has a non-positive raw
response, so its frame-wise ratio is unstable; the tabulated median and
quantiles should remain descriptive.

## Support, provenance, and validation

Both batch summaries report 57/57 frames and zero execution errors. The
position/PSF audit verified every selected frame against the approved 744-file
PSF plan. The two runs use identical code hashes for nine Astrafek source
modules and identical data/PSF input hashes. The support-aware audit also
preserved the valid ICL response for all 57 frames.

One i2 frame (cluster 2984) lost support for the satellite 50–75 kpc metric in
the `icl_only` scenario at 500 kpc. That satellite cross-response is omitted
for this one scenario, while the ICL response remains available. The
500 kpc i2 smooth-diffuse-to-satellite summary therefore has 43 frames; other
reported ICL metrics use all 44 i2 frames. The satellite-only injection
response is supported in all 57 frames and is unchanged between diffuse
scales.

Slurm job **428959** ran both batches sequentially and completed in 9 min 48 s;
peak RSS was 911,448 KiB. The outputs and independent audits are outside Git:

| Run | Summary SHA-256 | Audit SHA-256 |
|---|---|---|
| 100 kpc | `6145e85df3c3fca9f520b07df63abc03dc5749c338194f6941d06fd9e2cb35b8` | `9fb3c0b9249abcbe6cdf8cba653795ccf7291a830fdac416cd7a301a4e44ada3` |
| 500 kpc | `86d751be9e101dd003666b61410273570c3bb8ad685bdf014fac451cf64f0934` | `5ba44d44dce85db86a2eca4c247ce4550ea8a9093fbf65c24e959003c8eaec45` |

The paired comparison JSON is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_redetection_scale_comparison_2026-10-20_02/comparison.json`
(SHA-256 `7cb45820786021aeae4345d6b26521aa371c3835f91c37d3881b194e08de079e`).
It records the nine runtime code hashes, all input summary/audit hashes, and
frame-level support counts. Reproduce it with
[`summarize_joint_icl_satellite_psf_scale_comparison.py`](summarize_joint_icl_satellite_psf_scale_comparison.py)
after checking both run audits. The batch launcher is
[`run_joint_icl_satellite_psf_redetection_scale_pair.sbatch`](run_joint_icl_satellite_psf_redetection_scale_pair.sbatch).

## Limits and interpretation

This is a method-transfer stress test, not an ICL detection or an
ICL–satellite association measurement. The 500 kpc profile is an extended
stress case, not a fitted universal ICL model. The calculation uses only the
compact position-dependent PSF cores; it omits extended PSF wings, injected
source Poisson noise, correlated coadd noise, and HSC's upstream sky
subtraction. It also does not apply controls, extinction, or passband
corrections. The 16th–84th percentiles are not confidence intervals and do
not account for spatial covariance. Keep `science_measurement=false`.

The next science gate is to apply a validated target-specific PSF-wing and
sky-transfer model symmetrically to target and matched control coadds, then
estimate the ICL and satellite response covariance with the final fixed masks
and independent spatial blocks. These injection responses must not be used as
multiplicative corrections to the observed profiles.
