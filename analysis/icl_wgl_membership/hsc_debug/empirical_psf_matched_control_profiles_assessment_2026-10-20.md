# Matched-control profile diagnostic for Gaia PSF stamps — 2026-10-20

## Result

An exploratory paired diagnostic compared the existing Gaia-star stamps on target
HSC coadds against source-masked profiles at the same angular offsets in the
already downloaded matched-control coadds. It tests whether local sky and mask
residuals at those positions could affect the extended-wing measurement. A
control position contains no matched star, while the target profile can retain
ICL curvature, so target minus control is not a PSF-only estimate.

The support audit admitted 319 same-subfilter stamps in 217 target frames. All
401×401 control footprints fit in their images. Three predeclared diagnostic
variants were run on those existing images:

| Run | Control masks and central gate | Complete paired radial bins / 14,036 | Best-supported outer bin |
|---|---|---:|---|
| `_01` | HSC DETECTED/NEGATIVE/NOT_DEBLENDED masks; require at least 99% open pixels within 12 px for every ring | 287 | i2, 64–96 px: 7 frames |
| `_02` | Astrafek SEP detections with matching Kron 2.5/3.5 masks; same strict central gate | 149 | i2, 64–96 px: 4 frames |
| `_03` | Same Astrafek SEP masks; central openness gates only rings inside 12 px, while outer rings use their own masked area and 8-sector coverage | 279 | i2, 64–96 px: 12 frames |

The less restrictive outer-ring gate recovers additional support, but the best
outer bin still has only 12 frames, and other filter/mask/sky choices have 1–6
frames. In `_03`, frame-equal target-minus-control means in the i2 64–96 px bin
are `6.46e-7` and `1.48e-6` per target-core pixel for the 120–150 and 160–190 px
sky annuli, respectively. Their descriptive frame distributions include zero.
The corresponding i-band means are `2.49e-6` and `1.48e-6`, from only 5 and 3
frames. These small, sparse differences are consistent with residual background
or target-profile curvature and cannot establish the presence, shape, or
amplitude of extended PSF wings.

The strict variants show that support is highly sensitive to source masks and
to whether central-source openness is incorrectly required for outer annuli.
Even the recovered `_03` sample is too small and too selection-sensitive for a
calibration. The q16/q84 values are descriptive distributions across frames,
not confidence intervals; no bootstrap, coadd covariance, photon-noise
injection, or Gaia epoch correction was included.

## Reproducibility

The runner is [`analyze_empirical_psf_matched_control_profiles.py`](analyze_empirical_psf_matched_control_profiles.py).
The source data were the existing HSC target/control coadds, support table, and
hash-verified Gaia stamps; no new HSC query or download occurred. The test suite
for the support and profile diagnostic contains 7 tests.

Outputs are outside Git under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/`:

| Run | Output directory | `summary.json` SHA-256 |
|---|---|---|
| `_01` | `empirical_psf_matched_control_profiles_2026-10-20_01` | `17951c4736e1e4f350e704e0c85db592435995b231bf28d6a6f2ab546f8f0e4c` |
| `_02` | `empirical_psf_matched_control_profiles_2026-10-20_02` | `f7733a8201705ed883eca7a3e8e29fe649ef93eafc81ead45abcc048738056b6` |
| `_03` | `empirical_psf_matched_control_profiles_2026-10-20_03` | `7478c905a4e46c0c82f3fe569e4ee1b4cfb53280c9acb635267bc4188a98cd5b` |

Independent checks confirmed `_03` paired-table SHA-256
`5d57b6b4faa0b8003e0b807e91f001ad8129582ee33907aaf796cd9f9a1e8178` and
frame-weighted table SHA-256
`84c36fd7b6f47ae6d73dbbae405819f4b409c127d8c4b45307e5d7b06d095b8f` against
the run summary. Its status counts are 279 complete profile pairs, 9,449 blocked
by target background coverage, 1,449 blocked by control source/mask at the
center, 1,323 blocked by control background coverage, 748 blocked by target
core masks, 595 by target radial coverage, 105 by too few control quality pixels,
and 88 by radial coverage.

## Interpretation and next gate

This diagnostic shows that matched controls can help screen local background
and masking systematics, but the available same-filter stellar sample does not
support a stable empirical wing calibration. Preserve the published reference
PSF only as a model sensitivity input. Before using any wing correction in ICL
measurements, a larger same-filter isolated-star sample, epoch-aware Gaia
positions, spatial/coadd uncertainty, and diffuse-light plus star injection
through the full detection and masking pipeline are still required. Keep
`psf_approved=false` and `science_measurement=false`.
