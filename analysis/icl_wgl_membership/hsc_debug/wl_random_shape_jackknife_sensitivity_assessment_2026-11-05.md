# Random-subtracted weak-lensing jackknife block-width sensitivity (2026-11-05)

## Scope

This offline sensitivity uses the already completed HSC S19A shape-row export and
changes only the equal-area jackknife cell width. No HSC query or data download
was performed. The 5-degree reference has 51 occupied cells; the 2.5-degree run
has 154. Both use the same 435-lens common-support target sample, 2,175 associated
random centers, 8,905,892 random source pairs, 7 radial bins, HSC Y3 n(z), B-mode
bit choice, cosmology, and estimator.

The 2.5-degree run used the existing Astrafek stacker at low priority on CPU. It
completed in approximately 90 minutes wall time. The machine summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_02p5/summary.json`,
SHA-256 `653d27ae479175d6031daceddbc4fe5ff68526bac9e8fd1bc0a4e88786eab396`.
The 5-degree comparison summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json`,
SHA-256 `83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8`.

## Results

The profile point estimates are identical between the two partitions to the
stored precision. The tangential ΔΣ diagonal jackknife errors and cross-component
errors are:

| Proper radius (Mpc) | ΔΣ error, 5° | ΔΣ error, 2.5° | Cross error, 5° | Cross error, 2.5° |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | 1.524 | 1.420 | 1.525 | 1.641 |
| 0.2–0.4 | 0.764 | 0.815 | 0.770 | 0.761 |
| 0.4–0.7 | 0.597 | 0.579 | 0.450 | 0.508 |
| 0.7–1.0 | 0.413 | 0.423 | 0.372 | 0.356 |
| 1.0–1.5 | 0.310 | 0.276 | 0.198 | 0.253 |
| 1.5–2.0 | 0.311 | 0.255 | 0.210 | 0.225 |
| 2.0–3.0 | 0.172 | 0.170 | 0.148 | 0.133 |

Errors are in units of `10^13 Msun Mpc^-2`. The 2.5°/5° ΔΣ error ratio spans
0.820–1.067, with the largest decrease (18%) at 1.5–2 Mpc. The 7-component
covariance-weighted cross null statistic is χ²=10.594 (p=0.157) for 5° and
χ²=7.294 (p=0.399) for 2.5°. These are conditional diagnostic comparisons; the
p-value difference is not evidence that either partition is the correct survey
covariance.

Both 14×14 matrices are finite, full rank, symmetric to relative floating-point
roundoff (maximum relative asymmetry below 3.2×10⁻¹⁷), and positive definite.
The smallest eigenvalues are 9.29×10²³ (5°) and 1.31×10²⁴ (2.5°). This numerical
property does not establish calibrated uncertainty coverage.

The 2.5-degree radial profile figure is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_02p5/random_subtracted_profile.png`,
SHA-256 `9428ac5efddf8c01a0cc716c5da0ad3d831d5fe40f23aa3949884fb10dc9f503`.

## Interpretation and limitations

The positive declining ΔΣ point profile is unchanged, but its diagonal uncertainty
is modestly partition-sensitive; the cross null statistic also shifts. The result
therefore remains exploratory and `science_measurement=false`. The 2.5-degree
check is only one alternative partition, not a complete block-scale study.
Equal-area cells can still combine disconnected HSC regions, and random-center
full-aperture depth, mask, PSF, and sky matching remain unverified. Photo-z
calibration uncertainty and reduced-shear correction are also absent.

The existing implementation reruns the full target and random estimator for every
occupied block. The 154-block run took about 90 minutes, so further grid sweeps
should follow an Astrafek implementation that reuses per-pair sufficient
statistics while preserving the exact delete-either-endpoint rule. Do not treat
that optimization as scientifically equivalent until it reproduces the current
full-recompute results to a documented tolerance.

## Reproduction

Run from the Astrafek repository with the same input paths recorded in
`run_wl_random_shape_stack_after_export.py`, changing the output to a fresh path
and passing `--block-width-deg 2.5`. The saved summary contains the exact input
batch hashes and estimator provenance. Plot it with
`plot_wl_random_shape_stack.py`; that script reads the occupied block count from
the summary rather than assuming 51 cells.
