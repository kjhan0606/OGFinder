# Visit-split covariance pilot for ICL noise characterization

## Scope and method

The pilot uses only the already acquired, approved PDR3 Wide HSC-I warp
archives and corresponding coadds for the supported `i2` pair with cluster ID
101 and its matched control. The same seven expected visit IDs were assigned to
halves by alternating the sorted visit list, with the identical assignment for
target and control. Each half-stack uses inverse-variance weights after mapping
warp pixels to the delivered coadd grid with the affine WCS fit used by the
existing warp audit and nearest-neighbor sampling. The normalized difference
is `(A-B)/sqrt(VAR_A+VAR_B)`.

For the 200–300 proper-kpc annulus, the pilot measures normalized residual
covariance in 64×64 patches after a local affine-plane fit and iterative 4σ
clipping. The center mask uses finite coadd image/variance and excludes the
standard bad, no-data, bright-object, and related mask planes. Patch windows
overlap at 32 pixels, so their lag estimates are not independent observations.

## Results

Both centers have seven planned visits, but one visit contributes no valid warp
pixels in each center. The usable split counts are 3/2 visits for the target
and 4/2 for the control. Common clean support covers 67.5% of the geometric
annulus for the target (139,132 pixels) and 58.3% for the control (120,024
pixels).

Before the local plane fit, the normalized split-difference RMS is 1.093 for
the target and 1.166 for the control. After plane fitting and clipping, the
measured normalized zero-lag covariance is 1.182 for both. The symmetric
cumulative lag-kernel sum through ±8 pixels is 0.948 for the target and 1.187
for the control. These one-pair diagnostics show that the variance-normalized
split residual and its long-lag behavior are not yet explained by the simple
inverse-variance model. They do not identify which part comes from residual
sky, unmasked source structure, visit-dependent PSF/masks, variance calibration,
or the approximate warp-to-coadd mapping.

## Interpretation and limits

This is a method pilot on one target/control pair, not a population result.
The nearest-neighbor mapping does not reproduce the HSC coadd resampling
kernel; duplicate source-pixel mappings can induce covariance. No white-noise
simulation has yet been passed through the same WCS mapping and plane-fitting
estimator. The two split halves also have unequal effective visit support, and
the measured lag kernel is conditional on the chosen mask and estimator.
Therefore no noise correction was applied to the ICL stack, no association
claim is made, and `science_measurement=false` remains.

## Reproducibility

- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_visit_split_covariance_pilot.py`
- Runner SHA-256: `1f28f66ed9477abbef03bb53188f0cb14cd9354498a94eaf43989802be68c929`
- Output summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_visit_split_covariance_pilot_2026-11-06_03/summary.json`
- Output summary SHA-256: `09e557de604d4161cd0281065a59f32669744ab480b180029d192f79ba192359`
- Parent warp plan SHA-256: `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`
- Parent download summary SHA-256: `0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4`

No HSC request/download, GPU, or Slurm job was used for this analysis.
