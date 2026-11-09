# Target-only injection cross-talk in the ICL–satellite measurement (2026-11-06)

## Scope

Astrafek joined the observed ICL/satellite association sample to the audited
PSF-convolved injection responses. The analysis tests whether injected satellite
light can enter the 200–300 kpc ICL annulus. The injection responses came from
target cluster images only. The target-minus-control ICL contrast subtracts a
matched random control held fixed during source injection. The incremental
target response therefore measures the signal response of the target-minus-
control contrast. Injecting the same cluster source into a random control field
would remove part of the intended signal and is not required for this response.

The parent association sample contains 94 pairs, 26 in `i` and 68 in `i2`.
Valid injection responses overlap 58 pairs, 14 in `i` and 44 in `i2`. Thirty-six
pairs lack a valid transfer response. Thirty-one lack an approved PSF injection
record, and five lack satellite-annulus support. The overlap is selected by PSF
availability and image coverage.

## Measured target-side response

| Filter | Valid pairs | ICL response to diffuse injection | ICL response to satellite injection | Satellite response to diffuse injection | Satellite response to satellite injection | Median-response ratio `b/d` |
|---|---:|---:|---:|---:|---:|---:|
| `i` | 14 | 0.211 (q16–q84: 0.181–0.230) | 0.357 (0.125–0.579) | −0.035 (−0.086–0.093) | 0.349 (−0.306–0.955) | 1.022 |
| `i2` | 44 | 0.224 (0.183–0.252) | 0.490 (0.195–0.837) | 0.016 (−0.184–0.184) | 0.409 (−0.384–1.263) | 1.198 |

The target-side median response places injected satellite light in the ICL
annulus at a level comparable to the diffuse-injection response. Satellite
profile responses vary broadly and cross zero. Per-pair `b/d` ratios are
therefore unstable. The two-component response matrix has median condition
numbers 4.10 in `i` and 6.79 in `i2`, with broad q16–q84 ranges of 1.72–16.03
and 3.11–17.53, respectively.

The observed redshift- and richness-adjusted partial correlation in the full
94-pair sample is 0.082 in `i` and 0.285 in `i2`. In the selected injection
overlap, the corresponding raw partial correlations are 0.188 and 0.399.

## Leakage stress test and interpretation

A deliberately simplified stress test subtracted the measured satellite
profile multiplied by the ratio of median target-side responses. The stress
test changed the partial correlation to −0.743 in `i` and −0.567 in `i2`.
The paired 5-degree block-bootstrap change intervals were [−1.660, 0.020] and
[−1.464, −0.501]. The adjusted-correlation intervals were [−0.996, 0.910] and
[−0.963, −0.078].

The stress test is not a correction and does not establish a physical
anti-correlation. Only 58 selected pairs have transfer responses, and broad
response distributions make the median ratio an uncertain de-mixing
coefficient. The test also omits injected source-Poisson noise, coadd weights,
extended PSF wings, and upstream HSC sky transfer. The stress test shows that a
single median leakage coefficient can strongly change the association estimate.
The analysis does not establish a corrected population relation because the
pairwise inversion is unstable and omits response uncertainty, measurement
covariance, and support-selection uncertainty.

A separate pairwise matrix inversion was subsequently run on the 58-pair overlap.
The inversion produced many negative and extremely large amplitudes. The
conditional `i2` correlation interval excludes zero, but the interval omits
response uncertainty and the inverted amplitudes violate the expected physical
range. The inversion is therefore an instability diagnostic, not an association
correction. See the [matrix-inversion assessment](icl_satellite_transfer_matrix_deconvolution_assessment_2026-11-06.md).

## Required next calibration

The next calibration must apply independent diffuse and satellite injections
to the target image, keep the actual matched random control fixed, and rerun
the exact target-minus-control ICL and satellite estimators. The calibration
must span the 94-pair support or establish representative transfer weights
from the selected 58-pair subset. The analysis must propagate pair-dependent
matrix uncertainty and spatial-block uncertainty through a matrix-based
association estimate. A separate null run on random control fields can measure
false-positive response without adding cluster light to the controls. Injected
source Poisson noise, coadd pixel weights, extended PSF wings, and upstream sky
transfer remain separate calibration requirements. Until those gates pass,
`science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_icl_satellite_injection_cross_talk.py`
- Runner SHA-256: `174c1c3a868a7549f7347db8458a1addf12e0f588979ecc060eb1264689d5028`
- Corrected machine summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_cross_talk_2026-11-06_02/summary.json`
- Corrected summary SHA-256: `3e48028bef2d5f58a9815bb3568197ba54d578fd4036c99ad9764448996b667b`
- Corrected output-hash file SHA-256: `39f902cf62e2ea45ef48f824b06e25cc446122dbcd55efd374bb8ccd143eaba5`
- Corrected runner script SHA-256: `174c1c3a868a7549f7347db8458a1addf12e0f588979ecc060eb1264689d5028`
- Original v1 machine summary remains preserved in the `_01` output directory as historical provenance; use the `_02` interpretation.
- Output hashes: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_cross_talk_2026-11-06_01/output_hashes.json`
- Hash-file SHA-256: `286ad88d2d3168ca266ff53d921f446fd44f21ab55ade937643eb57eaf30c365`
- HSC queries and downloads: none
- Science gate: `science_measurement=false`
