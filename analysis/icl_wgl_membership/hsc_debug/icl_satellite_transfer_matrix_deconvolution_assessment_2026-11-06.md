# Pairwise ICL–satellite transfer-matrix inversion (2026-11-06)

## Estimator

The measured response matrix for each supported cluster/filter pair is

\[
\begin{pmatrix} I_{\rm obs} \\ S_{\rm obs} \end{pmatrix}
=
\begin{pmatrix} a & b \\ c & d \end{pmatrix}
\begin{pmatrix} A_{\rm diffuse} \\ A_{\rm satellite} \end{pmatrix}.
\]

The matrix columns are measured responses to diffuse and satellite injections.
Target-only injections measure the incremental response of the
target-minus-control ICL contrast when the matched random control remains fixed.
The analysis solves the 2×2 system for every pair and calculates the partial
correlation between the two inferred amplitudes after removing redshift and log
richness. The connected 5-degree block bootstrap resamples clusters while
holding all measured response matrices fixed.

## Results

| Filter | Pairs / blocks | Observed partial *r* | Pairwise inverted partial *r* | Conditional inverted 95% block interval | Median matrix condition number (maximum) |
|---|---:|---:|---:|---:|---:|
| `i` | 14 / 8 | 0.188 | −0.961 | [−0.999, 0.802] | 4.10 (17.17) |
| `i2` | 44 / 22 | 0.399 | −0.381 | [−0.971, −0.281] | 6.79 (126.72) |

The inverted amplitudes are unstable. Six of 14 `i` diffuse amplitudes and
seven of 14 `i` satellite amplitudes are negative. In `i2`, 21 of 44 diffuse
amplitudes and 18 of 44 satellite amplitudes are negative. The 16th, 50th, and
84th percentiles are [−26.19, 10.97, 50.50] for the `i2` diffuse amplitude and
[−8.74, 1.07, 12.54] for the `i2` satellite amplitude. Maximum absolute
amplitudes reach 1,839.55 and 2,211.23, respectively, in the adopted
nJy arcsec⁻² injection-amplitude units.

## Interpretation

The inversion demonstrates that pair-dependent cross-talk can change the
estimated association. The negative `i2` interval is conditional on the
observed response matrices and does not represent a calibrated physical
anti-correlation. The inversion produces many negative and extremely large
amplitudes, revealing instability or failure of the assumed zero-intercept,
linear two-component model for the current selected sample. The interval also
omits injection-response uncertainty, measurement covariance, support-selection
uncertainty, source Poisson noise, extended PSF wings, and upstream sky transfer.
The 58-pair overlap is selected from the 94-pair satellite-supported sample.

The current injection set therefore does not provide a science correction.
Improved calibration requires response measurements across representative
support, uncertainty estimates for every response-matrix element, and a model
that includes the observed background and measurement covariance. The
target-minus-control estimator should retain actual random controls fixed while
cluster light is injected into target images. A separate control-field null
tests false-positive response without inserting cluster light into random
fields. `science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_icl_satellite_transfer_matrix_deconvolution.py`
- Runner SHA-256: `042f94f3c90b01ee915516d8db271991beea7caa7f73887b1892e6d272ff2c4c`
- Input transfer summary SHA-256: `3e48028bef2d5f58a9815bb3568197ba54d578fd4036c99ad9764448996b667b`
- Output summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_transfer_matrix_deconvolution_2026-11-06_01/summary.json`
- Output summary SHA-256: `aa3566d8bfe537bb249a06c0f0b53db0bf7f70ac33d9f55c195e0e3df3066230`
- Output hash file SHA-256: `8e4e899aee23bf683c44c9665a5a526d69a6c43cff06ef51768d5cf7695e0c8f`
- HSC requests and downloads: none
