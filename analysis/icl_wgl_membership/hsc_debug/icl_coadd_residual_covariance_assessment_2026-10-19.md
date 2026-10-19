# Matched HSC coadd residual-covariance assessment — 2026-10-19

## Scope and execution

Astrafek compared 57 target coadds with the 57 already-downloaded random
control coadds matched to the same cluster and filter. The target masks use the
saved SEP/Kron 2.5 mask with detected segmentation cores excluded; the control
masks were regenerated with the same SEP/Kron 2.5 procedure. No HSC query or
download was made. Slurm job **428608** completed in 2 min 2 sec (peak RSS
359 MiB; exit 0).

For each frame, the runner selected up to 96 non-overlapping 64×64 pixel
patches with at least 80% valid pixels. It fitted a local affine plane with
iterative 4-sigma clipping, divided residuals by the delivered variance
standard deviation, and measured median covariance at each of 81 lags within
±4 pixels. The target and control patch residuals both include unresolved
sources and any remaining diffuse-light contributions.

## Results

The sum of the 9×9 lag-covariance grid, relative to the delivered diagonal
variance, is an approximate large-top-hat residual-variance factor. It is not
the exact variance of the ICL or satellite estimator.

The q16–q84 intervals below describe frame-to-frame spread, not confidence
intervals.

| Filter | Frames | Target factor (q16–q84) | Matched-control factor (q16–q84) | Paired target−control median (q16–q84) |
|---|---:|---:|---:|---:|
| i | 13 | 1.594 (1.516–1.687) | 1.545 (1.519–1.607) | −0.012 (−0.083–0.095) |
| i2 | 44 | 1.588 (1.518–1.681) | 1.541 (1.470–1.609) | +0.034 (−0.019–0.110) |

The median zero-lag residual variance is 1.145 times the delivered variance in
targets and 1.143 in controls. The lag pattern is similar in both: median
variance-normalized axial covariance is about +0.06 at one pixel, −0.04 at two
pixels, +0.04 at three pixels, and −0.02 at four pixels. Paired target-minus-
control factor intervals span zero in both filters. This is consistent with a
broadly shared short-range coadd/residual covariance component; it does not
separate instrumental resampling covariance from all astrophysical or
unresolved-source contributions.

If this approximate factor were applied to a large top-hat sum, it would imply
about 1.24–1.26 times the nominal diagonal-noise standard deviation. Do not
apply it as a correction to the ICL–satellite result: the estimate uses a
finite lag window, is not propagated through the exact aperture weights, and
does not include the external ICL-control term, source-photon noise, or
upstream sky-transfer uncertainty. Spatial-block uncertainty already
contains field-to-field variation, so adding this factor directly could
double-count some fluctuations.

## Reproducibility and status

- Astrafek runner: `estimate_hsc_coadd_noise_correlation.py`
- Source-mask implementation: `prepare_detection_psf.py`
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_noise_correlation_2026-10-19_04/`
- `summary.json` SHA-256: `68a50489629f20da927acf761b6925125c2188cca1777b2531136f11cd9a3be4`
- `frame_noise_correlation.tsv` SHA-256: `30d5a000d35f31445cd11c3bbe4a0a22c660ebe3e4a2ab29ef08a00f342affd5`
- `lag_correlation.tsv` SHA-256: `a72d335d763b15599dec70b745bad7fe6bdb6f9a1176507953dc3094a6069015`
- Runner SHA-256: `cad0cd3bd619737d6aad3a716789c9b360194618664600a33139e5094a562b39`

An independent pass verified the input, code, and output checksums; all 114
frame rows; all 9,234 lag rows; each frame's factor recomputed from its lag
grid; and the target/control filter summaries. Earlier target-only or
unequal-mask trials are retained as diagnostics but are not the reported
matched result. `science_measurement=false` remains in force.

The next useful ICL step is to propagate these per-frame lag kernels through
the actual signed ICL and satellite pixel weights, then compare that noise-only
covariance with the paired spatial-block bootstrap. Keep it a sensitivity
analysis until source-Poisson and sky-transfer terms are covered.
