# Global weak-lensing profile null audit

## Scope

Evaluated the saved random-subtracted HSC S19A ΔΣ profile against a zero
signal using the complete paired endpoint-jackknife covariance. The stack uses
435 common cluster lenses and 2,175 matched random centers; the two saved
covariance partitions have 51 occupied 5-degree blocks and 154 occupied
2.5-degree blocks. The profile point estimates are identical between
partitions. The analysis uses no new HSC query or download and preserves
`science_measurement=false`.

## Results

The covariance-aware tangential profile is inconsistent with zero in every
radial-range check below, conditional on the stored covariance and the
chi-square approximation. Dropping inner bins reduces concerns about the
innermost annulus but leaves a large positive profile.

| Jackknife blocks | Radius included (Mpc) | χ² vs zero | dof | nominal p |
|---|---|---:|---:|---:|
| 51, 5° | 0.1–3.0 | 479.81 | 7 | 1.76×10⁻⁹⁹ |
| 51, 5° | 0.4–3.0 | 460.37 | 5 | 2.85×10⁻⁹⁷ |
| 51, 5° | 0.7–3.0 | 331.37 | 4 | 1.84×10⁻⁷⁰ |
| 154, 2.5° | 0.1–3.0 | 642.83 | 7 | 1.45×10⁻¹³⁴ |
| 154, 2.5° | 0.4–3.0 | 537.23 | 5 | 7.32×10⁻¹¹⁴ |
| 154, 2.5° | 0.7–3.0 | 406.04 | 4 | 1.38×10⁻⁸⁶ |

The cross-component null remains acceptable: χ²=10.59 for 7 bins
(“5°”, p=0.157) and χ²=7.29 (2.5°, p=0.399). Both 14×14
covariances are full rank and positive definite. Their condition numbers are
270 (5°) and 216 (2.5°). The change in nominal p-value between partitions
shows that uncertainty remains partition-sensitive even though the point
profile is fixed.

## Interpretation and limits

This is a strong conditional tangential-shear profile with a null-consistent
cross component in the current common-support sample. The quoted probabilities
assume the jackknife covariance is adequate and the chi-square approximation
applies; they do not include uncertainty in the covariance estimate or
calibration systematics. Equal-area blocks are not verified contiguous HSC
survey regions. Random-center full-aperture depth/mask/PSF/sky matching,
photo-z dilution calibration, B-mode-mask polarity, reduced shear, and
correlated coadd noise remain open. Consequently this result supports the
weak-lensing pipeline and motivates those remaining gates, but it is not yet a
publication-grade detection claim. The machine status remains false.

## Reproduction

- Astrafek analysis: [`analyze_wl_random_shape_stack_detection.py`](analyze_wl_random_shape_stack_detection.py), SHA-256 `e0934796e204acaec678a2a7dde5ae971920d098dae24fc1539659c1dc93bed1`.
- Saved result: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_detection_audit_2026-11-05_01/summary.json`, SHA-256 `ed6040a61157df74532fccf3e67e23c6f8231a33fe953045fd19ae35acc1cffb`.
- The result references the exact 5-degree and 2.5-degree stack summary hashes and evaluates radial cuts from the same 7-bin ΔΣ vector.
