# Exploratory NFW mass scale for the random-subtracted HSC stack (2026-11-05)

## Scope and estimator

Astrafek fitted a spherical NFW weak-shear model to the saved random-subtracted
HSC S19A stack after excluding source tomographic bin 1. The model averages

\[
\Delta\Sigma_{\rm NFW}(R_i,z_l;M_{200c},c)
\]

over the exact 435 matched cluster-lens source rows, using each pair's observed
proper projected radius, lens redshift, HSC Y3 source-bin \(n(z)\) lensing
efficiency, B-mode/source cuts, and the same per-bin responsivity, multiplicative
and selection-bias weights as the measured weak-shear stack. The fit minimizes
the full-covariance \(\chi^2\) using the tangential submatrix of the saved
paired spatial-jackknife covariance. The standard analytic projected NFW
profile follows [Wright & Brainerd](https://arxiv.org/abs/astro-ph/9908213).

The radial fit uses 0.2–1.5 proper Mpc. It excludes 0.1–0.2 Mpc, where the
tomographic analysis found the largest low-source-bin sensitivity, and 1.5–3
Mpc, where a one-halo-only model is inadequate. Concentration is fixed at 3, 4,
or 5; the run is repeated with the saved 5° and 2.5° endpoint-jackknife
covariances. The row-level forward average avoids replacing the 435-lens
redshift and pair-weight distribution with one effective redshift.

## Results

Mass is an effective single-halo stack scale, not the arithmetic mean mass of
the member clusters. The intervals below are conditional \(\Delta\chi^2=1\)
profile intervals for fixed concentration and saved covariance.

| Fixed concentration | 5° jackknife: \(M_{200c}\) | 5° \(\chi^2/3\) | 2.5° jackknife: \(M_{200c}\) | 2.5° \(\chi^2/3\) |
|---:|---:|---:|---:|---:|
| 3 | \(1.78\times10^{14}\ M_\odot\) [log interval 14.218–14.282] | 1.39/3 | \(1.80\times10^{14}\ M_\odot\) [14.224–14.283] | 2.11/3 |
| 4 | \(1.67\times10^{14}\ M_\odot\) [14.191–14.252] | 7.43/3 | \(1.59\times10^{14}\ M_\odot\) [14.173–14.230] | 6.15/3 |
| 5 | \(1.54\times10^{14}\ M_\odot\) [14.155–14.216] | 23.59/3 | \(1.44\times10^{14}\ M_\odot\) [14.130–14.185] | 14.95/3 |

The two block widths give similar fits at fixed concentration. Across the
fixed-\(c\) grid, the fitted mass changes by roughly 20% from \(c=3\) to 5.
The \(c=5\) profile has a large residual \(\chi^2\) for either covariance;
this is a model-shape warning under the selected radial range, not a
concentration measurement. The much narrower conditional mass intervals do
not include the concentration or other systematic uncertainties.

## Interpretation limits

This is a model screen that confirms a positive stacked signal can be expressed
as an effective halo scale. It is not a science-grade mass calibration and
does not promote the lensing or ICL result. The fit assumes one mass and a
spherical NFW one-halo profile for a richness-selected population with a mass
distribution. It omits concentration variation, reduced shear, magnification,
miscentering, two-halo structure, photometric-redshift systematic uncertainty,
random-aperture depth/mask/PSF/sky mismatch, and covariance uncertainty. The
nominal fit probabilities are conditional on the selected profile bins and
jackknife covariance. Keep `science_measurement=false`.

## Reproducibility

Astrafek runner: `fit_wl_random_stack_nfw_mass.py` (SHA-256
`b8999e9e534def7eeeff9f29921e9ed479d57eb31e935d21c1deaba2c271c3e2`).
Machine result: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_nfw_mass_2026-11-05_01/summary.json`
(SHA-256 `d3ca78a24fc030f64694247fd45d33f607d89c18fb102b6b5a6879d53c52cd4b`).
The summary hash-binds the stack, row-integrity audit, tomography profile,
full target stack, shape-row plan, and HSC Y3 `n(z)` SACC inputs.
