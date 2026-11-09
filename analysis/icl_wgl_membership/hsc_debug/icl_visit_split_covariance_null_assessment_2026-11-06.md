# Heteroscedastic white-noise propagation through the visit-split estimator

## Setup

This offline Monte Carlo uses the supported `i2` cluster 101 and its matched
control, reusing only the already approved and downloaded HSC-I warps. For each
center it draws independent Gaussian noise on every native warp pixel using
that visit's supplied variance plane. It then applies the saved valid-pixel
masks, the exact target/control visit split, affine WCS mapping with
nearest-neighbor sampling, inverse-variance half-stacks, the `(A-B)/sqrt(VAR_A+
VAR_B)` normalization, local affine-plane subtraction, iterative 4σ clipping,
and the same overlapping-patch ±8 covariance estimator. One hundred
realizations were run separately for each center with fixed seeds.

The null reproduces the project's approximate nearest-neighbor mapping. It
therefore includes covariance when adjacent coadd pixels reuse a native warp
pixel, as well as mask geometry, heteroscedasticity, and the plane/clipping
estimator. It does not include covariance already present in the observed
warps, residual sky, astronomical source residuals, PSF differences, or HSC's
production resampling kernel.

## Results

| Center | Quantity | Observed split | White-noise mean | White-noise 16th–84th percentile |
|---|---|---:|---:|---:|
| Target | Raw normalized RMS | 1.0934 | 1.0000 | 0.9981–1.0019 |
| Target | Plane-fitted zero-lag covariance | 1.1819 | 0.9976 | 0.9940–1.0014 |
| Target | Cumulative covariance through ±8 px | 0.9481 | 0.7156 | 0.6582–0.7675 |
| Control | Raw normalized RMS | 1.1664 | 0.9997 | 0.9976–1.0015 |
| Control | Plane-fitted zero-lag covariance | 1.1815 | 0.9972 | 0.9937–1.0008 |
| Control | Cumulative covariance through ±8 px | 1.1867 | 0.7147 | 0.6502–0.7859 |

The pure-noise simulation reproduces a substantial reduction in the cumulative
lag sum after local plane removal, confirming that this estimator itself
creates negative long-lag tails. However, all three observed quantities lie
above the simulation's central 68% range at both centers. The independent
white-noise model with delivered variance planes and nearest-neighbor mapping
is therefore insufficient to explain the observed split residuals. This is a
model-check result, not a p-value: there are only two centers and the null omits
several real image components.

## Interpretation and limits

The excess can arise from correlated noise already present in warp images,
residual sky structure, unmasked astronomical sources, variance-plane
miscalibration, or other reduction effects; this pilot does not separate them.
In particular, it does not validate a scalar covariance correction. Do not
apply the simulated kernel to the ICL stack or to weak-lensing measurements.
`science_measurement=false` remains.

The next calibration step is to estimate spatial covariance empirically in
source-masked sky annuli and matched controls, then test whether it transfers to
the ICL aperture. A broader simulation can then perturb the supplied variance
normalization within limits supported by independent sky measurements.

## Reproducibility

- Runner: `analysis/icl_wgl_membership/hsc_debug/simulate_icl_visit_split_covariance_null.py`
- Runner SHA-256: `7c4dcb02f7d6b0b31cd3a1c85f47923c2e6e9ffbef75520f0687143b64d57410`
- Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_visit_split_covariance_null_2026-11-06_01/summary.json`
- Summary SHA-256: `d9dc7090b1cee14fdb88a020adc5164a6da610bcaa76d9d190142b8a2ad86890`
- Realizations per center: 100
- Parent warp plan SHA-256: `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`
- Split-pilot summary SHA-256: `09e557de604d4161cd0281065a59f32669744ab480b180029d192f79ba192359`

No HSC request/download, GPU, or Slurm job was used.
