# Finite-lag pixel covariance for the ICL–satellite estimators — 2026-10-19

## Purpose and scope

This calculation contracts the actual signed pixel weights of the target-side
ICL 200–300 kpc annulus and the satellite 50–75 kpc target-minus-azimuthal-
control estimator with the measured, symmetrized 9×9 short-range residual
covariance kernel for each HSC coadd. It extends the earlier independent-pixel
diagonal calculation to include local pixel correlations. The analysis uses
the existing 57 matched cluster/filter frames (13 `i`, 44 `i2`) and 156
satellite profiles. Six otherwise supported frames have no measured kernel and
are excluded. No new HSC query or download was made.

The implementation is in
[`audit_satellite_icl_shared_pixel_covariance.py`](audit_satellite_icl_shared_pixel_covariance.py),
with tests in
[`test_audit_satellite_icl_shared_pixel_covariance.py`](test_audit_satellite_icl_shared_pixel_covariance.py).
The CPU Slurm entry point is
[`run_satellite_icl_lag_covariance.sbatch`](run_satellite_icl_lag_covariance.sbatch).

## Results

The median absolute-noise standard deviation rises by about 25% when the
measured finite-lag kernel is applied. The normalized pixel-noise correlation
between the two estimators remains small and close to the diagonal estimate.
The intervals below are the 16th–84th percentiles across frames, not confidence
intervals.

| Filter | Frames | Median ICL sigma factor | Median satellite sigma factor | Median ICL–satellite noise correlation [q16, q84] | Diagonal correlation median |
|---|---:|---:|---:|---:|---:|
| i | 13 | 1.251 | 1.244 | −0.0109 [−0.0225, +0.0099] | −0.0112 |
| i2 | 44 | 1.249 | 1.244 | −0.0014 [−0.0227, +0.0080] | −0.0016 |

This indicates that independent-pixel variance understates the target-side
aperture noise for these weights, while the short-range kernel changes the
estimated cross-estimator correlation little. It does not imply that the
observed ICL–satellite association is absent or significant: the correlation
here is a noise-only covariance coefficient, not the astrophysical profile
correlation.

## Reproducibility and limits

Slurm job **428612** completed on CPU in 1 min 16 sec (peak RSS 1.31 GiB,
exit 0). Output is stored outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lag_covariance_2026-10-19_01/`.
The audited `summary.json` SHA-256 is
`291870ac5469151452052ecd56aca3a94bd4cc392ea406b1417345b741fa720f`; the
`pair_pixel_covariance.tsv` SHA-256 is
`d0f4632426f5db691ac9a8466730c29833eadf32da1c5f3678d17c477ea13b2c`. An
independent pass rechecked all nine input hashes, the 57 unique frame rows,
the 13/44 filter counts, the output hash, and recomputed the reported medians
and quantiles from the TSV.

This is a stationary finite-lag patch model, not the full position-dependent
coadd covariance. It covers only the target-side ICL term; the matched random-
control ICL term is not included. Fixed masks, detections, support cuts, and
equal cluster/satellite weighting are assumed. Source Poisson noise, upstream
sky transfer, PSF, calibration, and astrophysical covariance are excluded.
Some spatial-block scatter may already contain part of this noise, so applying
the factor on top of that scatter could double-count uncertainty. The output
remains `science_measurement=false` and must not be used as a final error-bar
correction.

The next useful covariance step is a paired spatial-block resampling of the
target and matched-control ICL estimators with this local kernel sensitivity
reported separately, rather than multiplied into the existing block errors.
