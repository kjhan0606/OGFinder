# Expanded matched-control ICL–satellite pixel-noise covariance — 2026-11-03

## Sample and method

Astrafek expanded the target/control pixel-noise covariance audit from 57 to 63
cluster/filter pairs. The expanded sample contains 17 `i` pairs and 46 `i2`
pairs, with 162 accepted satellite profiles. The analysis retained the saved
support rules: at least 40% valid area in the cluster-centered 200–300 kpc ICL
annulus, at least 60% satellite-annulus coverage in the target and controls,
and at least three accepted azimuthal controls per satellite.

Astrafek reconstructed the actual signed pixel weights for the ICL and
satellite estimators. The variance-plane maps and the independently verified
full-sample 9×9 lag kernels supplied the pixel-noise covariance. Matched-control
ICL variance was added to target ICL variance. Target and control coadds were
treated as independent after checking the angular gaps between apertures.
Astrafek made no HSC query or download.

## Results

| Filter | Pairs | Control/target ICL lag-variance median | Total ICL sigma factor over diagonal-only variance | Median ICL–satellite noise correlation [16th, 84th percentile] |
|---|---:|---:|---:|---:|
| i | 17 | 1.027 | 1.250 | −0.00721 [−0.01578, +0.00554] |
| i2 | 46 | 0.990 | 1.241 | −0.00145 [−0.01641, +0.00473] |

The local lag model raises the median ICL aperture noise standard deviation by
about 24–25% relative to the delivered variance plane's diagonal term. The
matched-control contribution is comparable to the target contribution in both
filters. The per-pair ICL–satellite pixel-noise correlations remain close to
zero. The percentile ranges describe the selected pair distribution, not
confidence intervals.

The independent verifier confirmed all 63 target/control mappings and
aperture formulas. Every target and control annulus lies outside the measured
9×9 cross-aperture lag support. The minimum boundary gap is 297.497 arcsec,
compared with a maximum represented pixel-center separation of 0.9504 arcsec.

## Interpretation and limits

The finite-lag noise model predicts a moderate increase in per-pair ICL
aperture noise relative to independent-pixel propagation. The model predicts
little direct shared-pixel noise covariance between the cluster ICL annulus
and the satellite 50–75 kpc estimator. The near-zero pixel-noise covariance
does not identify the observed ICL–satellite relation as astrophysical. Sky
transfer, PSF wings, masks, completeness, and astrophysical covariance remain
outside the model.

The lag kernel is a stationary median over sampled image patches. The kernel
does not represent a position-dependent full covariance, and large-scale
target/control sky covariance remains unmeasured. Source Poisson, sky-transfer,
PSF, calibration, extinction, and mask-growth uncertainties remain excluded.
The analysis retains `science_measurement=false` and does not correct the
observed ICL–satellite profile.

## Reproducibility

- Target/satellite pixel covariance: [`audit_satellite_icl_shared_pixel_covariance.py`](audit_satellite_icl_shared_pixel_covariance.py)
- Target/control covariance: [`audit_satellite_icl_target_control_covariance.py`](audit_satellite_icl_target_control_covariance.py)
- Independent verifier: [`verify_satellite_icl_target_control_covariance.py`](verify_satellite_icl_target_control_covariance.py)
- Target/satellite output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lag_covariance_2026-11-03_01/`
- Combined output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_target_control_covariance_2026-11-03_01/`
- Combined summary SHA-256: `be808c38407ddffeb63ba59f7072371affbbe6adf41169391490b4422d71158d`
- Combined table SHA-256: `64d9123322d9d2d670908faa6c07de9ca812bb1ade00a17823f64f7092b0a0f4`
- Verification status: `verified`; `science_measurement=false`.
