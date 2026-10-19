# Matched-control contribution to ICL–satellite pixel covariance — 2026-10-19

## Purpose and inputs

The preceding finite-lag audit measured the ICL 200–300 kpc aperture noise on
the cluster target image only. This follow-up estimates the same annular
variance on each already-downloaded matched-control coadd, then adds the
target and control variances for the target-minus-control ICL estimator. The
ICL–satellite cross-covariance remains the shared target-image term. The run
uses the same 57 cluster/filter pairs and 156 satellite profiles (13 `i`, 44
`i2`) as the prior audit. It made no HSC query or download.

Control masks were regenerated with the existing SEP/Kron-2.5 method from the
matched-control `coadd/bg` images. Each control annulus passed the 40% coverage
floor. Target and control 200–300 kpc annuli are non-overlapping on the sky in
all 57 pairs. The variance calculation uses each image's delivered variance
plane and its matched 9×9 empirical lag kernel.

The reusable implementation is
[`audit_satellite_icl_target_control_covariance.py`](audit_satellite_icl_target_control_covariance.py),
with unit tests in
[`test_audit_satellite_icl_target_control_covariance.py`](test_audit_satellite_icl_target_control_covariance.py)
and the CPU Slurm wrapper
[`run_satellite_icl_target_control_covariance.sbatch`](run_satellite_icl_target_control_covariance.sbatch).

## Results

Adding the matched-control variance gives total ICL sigma inflation of about
24% relative to independent-pixel diagonal variance. The normalized
noise-only covariance with the satellite estimator remains small. Intervals
are frame-to-frame 16th–84th percentiles, not confidence intervals.

| Filter | Pairs | Median control/target ICL variance | Median total ICL sigma factor | Median total ICL–satellite noise correlation [q16, q84] |
|---|---:|---:|---:|---:|
| i | 13 | 0.825 | 1.234 | −0.00731 [−0.01812, +0.00703] |
| i2 | 44 | 1.005 | 1.240 | −0.00091 [−0.01476, +0.00604] |

Thus, for this aperture and mask selection, the matched-control image supplies
a variance contribution comparable to the target term in i2 and somewhat
smaller in i. Including it lowers the absolute normalized noise correlation
relative to the target-only estimate. These coefficients describe shared
pixel noise, not the astrophysical ICL–satellite profile correlation.

## Execution and independent audit

Slurm job **428622** completed on CPU in 1 min 19 sec (peak RSS 536 MiB,
exit 0). The authoritative output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_target_control_covariance_2026-10-19_03/`.
Independent validation rechecked all ten source checksums and the output TSV
checksum, confirmed 57 unique rows (13 `i`, 44 `i2`), verified control
coverage and disjoint apertures, and recomputed the combined diagonal/lag
variances, sigma factors, correlations, medians, and quantiles from the row
table.

- `summary.json` SHA-256: `961614b00d8c73d3464396c7f384b25599710f4b6b4115cf0131f04ebfab4bc2`
- `target_control_pixel_covariance.tsv` SHA-256: `35f59684196029f4bc1525aab079edf2c1616dbda9194d3ccd14ababe3dfb479`
- Analysis script SHA-256: `ac1d3fd474290ad44e8671979f0c18a5ac540d0ee7d9247abd6aadd2bc511df7`

## Interpretation and limits

The target and matched-control noise terms are added as independent. Their
annuli do not overlap, but the cross-image covariance is not directly measured;
large-scale sky-subtraction residuals could still couple separated coadds.
The lag kernel is a stationary local patch estimate rather than the exact
position-dependent full covariance. The masks/detections and selected sample
are fixed. Source Poisson noise, upstream sky transfer, PSF and calibration
uncertainty, extinction, and astrophysical covariance remain excluded. Spatial
block scatter may contain some of this noise, so multiplying the sigma factor
into block errors could double-count it. This remains an exploratory
noise-only audit with `science_measurement=false`, not a correction to the
observed ICL–satellite relation or a final uncertainty estimate.
