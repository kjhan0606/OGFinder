# HSC cluster-only weak-lensing spatial-jackknife assessment — 2026-10-19

## Run and input audit

Slurm job **428563** completed successfully (12 min 29 sec; peak RSS 2.36 GiB;
exit 0). The Astrafek runner deletes one occupied 5-degree equal-area cell at a
time, removes every lens–source pair whose lens or source falls in that cell,
and recomputes the tomographic calibration for each replicate. It found 52
occupied cells.

The stack input comprised 3,069,544 planned lens–source rows for 800 lenses.
After the HSC quality cuts, B-mode mask, foreground/zero-efficiency and radial
cuts, 2,153,004 pairs remained, involving 609 lenses and 2,051,991 unique
sources. The stack uses HSC Y3 fiducial tomographic `n(z)` with raw B-mode mask
value 0 excluded. The result is cluster-only: no random-lens subtraction,
matched-field variance, `n(z)` uncertainty, shear-calibration uncertainty, or
reduced-shear correction is included.

## Radial profile

Values below are weak-shear `DeltaSigma` in `10^12 M_sun Mpc^-2`. The quoted
error is the diagonal spatial-jackknife error; S/N is shown only as a compact
description of this provisional profile and is not a detection significance.

| R (Mpc) | DeltaSigma | Jackknife error | Descriptive S/N | Cross component / jackknife error |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 161.48 | 12.09 | 13.36 | +0.25 |
| 0.2–0.4 | 116.27 | 5.28 | 22.02 | −0.39 |
| 0.4–0.7 | 72.40 | 3.96 | 18.30 | −0.85 |
| 0.7–1.0 | 54.64 | 2.68 | 20.42 | −1.45 |
| 1.0–1.5 | 38.48 | 2.49 | 15.42 | −1.04 |
| 1.5–2.0 | 27.92 | 1.86 | 14.99 | +1.80 |
| 2.0–3.0 | 15.27 | 1.01 | 15.13 | +0.05 |

The tangential profile is positive and decreases with radius. The joint
cross-component diagnostic gives `chi2=9.30` for 7 covariance modes
(`p=0.232` under the stated Gaussian approximation). It finds no cross-shear
inconsistency in this diagnostic; it does not establish a lensing detection or
validate the tangential signal. The 14×14 stored covariance is symmetric
positive definite in a numerical eigenvalue check (smallest eigenvalue
`5.19e23`, largest `1.92e26` in its stored units; condition number about 370).

## Interpretation and next gate

This is useful evidence that the 800-lens shape-row pipeline produces a
well-populated, declining cluster-centered profile and a non-anomalous
cross-component diagnostic under a spatial deletion scheme. The covariance is
more informative than the earlier independent-pair shape-noise estimate, but
the 52 grid cells are only an approximate fallback: a cell can combine
disconnected HSC regions, and cells are not guaranteed to be independent
survey fields. The displayed bin S/N therefore must not be read as a
population-level detection.

The next weak-lensing gate is a matched random-point subtraction with the same
shape selection and a paired spatial covariance. That needs the already
approved export of random-point coordinates followed by footprint, mask,
depth, PSF, extinction, and aperture-overlap screening before any random-source
shape rows are requested. At the time of this report, the coordinate CAS
export checkpoint has two of 32 jobs complete and one submitted; no shape rows,
photo-z rows, or images were part of that export. Keep the existing submitter
running while it polls the submitted CAS job; do not resubmit that batch.

## Reproducibility

- Runner: `run_wl_cluster_spatial_jackknife.py`
- Slurm wrapper: `run_wl_cluster_spatial_jackknife.sbatch`
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_spatial_jackknife_2026-10-19_01/summary.json`
- Output SHA-256: `16579be08d02ee29f5ad912f9274f64b728bafc408ddca0746de3ad69494a6c3`
- Runner SHA-256: `b444c47b6c2f2dd700513541e65860a96a0bba8c5571517c22ce5067ba09c70a`
- The machine-readable result retains the full joint covariance, profile,
  selection counts, per-cell deletion counts, and input checksums.

`science_measurement=false` remains in force. This is an exploratory
cluster-only weak-lensing stack, not a final mass measurement.
