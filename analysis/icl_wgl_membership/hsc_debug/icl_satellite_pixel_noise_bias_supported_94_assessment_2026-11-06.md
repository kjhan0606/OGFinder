# ICL–satellite pixel-noise correction on the 94-pair support (2026-11-06)

Applied the saved finite-lag pixel-noise covariance to all 94 existing
ICL/satellite support pairs (26 `i`, 68 `i2`). This is the same residual
cross-product correction used for the earlier 74-pair diagnostic. An
independent support-scope audit verified that the covariance table covers
exactly the 94 selected pairs. The raw values for the nested 74-pair subset
match the earlier association table. No HSC query, download, or image read was
performed in this step.

| Filter | Pairs | Corrected partial *r* | 2.5° blocks (95% interval) | 5° blocks (95% interval) | 10° blocks (95% interval) |
|---|---:|---:|---:|---:|---:|
| `i` | 26 | 0.081 | 20: [−0.359, 0.610] | 13: [−0.195, 0.511] | 10: [−0.501, 0.524] |
| `i2` | 68 | 0.306 | 47: [−0.004, 0.523] | 33: [0.050, 0.519] | 23: [0.062, 0.512] |

The estimates partial out redshift and richness. The correction changes the
`i2` point estimate by only +0.013 and the `i` estimate by +0.002. The `i2`
point estimate is stable across block widths; its 2.5° interval reaches just
below zero, while the 5° and 10° intervals remain positive. The `i` intervals
are broad and include zero at every width. This indicates that the estimated
stationary finite-lag pixel noise is small relative to the observed
association, but uncertainty from spatial sampling is consequential near the
zero boundary.

The main 5° interval from the saved block labels is `[-0.203,0.494]` for `i`
and `[0.045,0.516]` for `i2`. The width-sensitivity run independently rebuilds
connected blocks from target/control coordinates; its 5° intervals are close
but not identical (`[-0.195,0.511]` and `[0.050,0.519]`) because the saved and
reconstructed block assignments differ slightly. Both implementations agree
on the scientific reading.

These intervals use 5,000 resamples of approximate connected 5° equal-area
blocks; they are not calibrated detection significances. The previous raw
block-width sensitivity showed borderline intervals at 2.5°, and those
smaller blocks are not independent validated HSC fields. The calculation is
conditional on fixed masks, detections, and support cuts. It does not include
support-selection bias, source Poisson noise, sky-transfer uncertainty,
extended PSF wings, extinction/calibration uncertainty, or astrophysical
covariance. It does not establish a population-level or causal relation;
`science_measurement=false` remains in force.

Reproduce all three block widths with
[`assess_icl_satellite_pixel_noise_bias_supported_94.py`](assess_icl_satellite_pixel_noise_bias_supported_94.py).
The independent covariance-scope check and machine-readable summary are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_pixel_noise_bias_supported_94_2026-11-06_04/`.
Summary SHA-256:
`965d0e81e78ad8d27ebfc7726a9ae170828b824fab41d29292398fbdf2940725`.
Runner SHA-256:
`4bd020428ba8ff2dc8ccb0a5bfa961e044f71985bdb79e31ff011da7138af7da`.
