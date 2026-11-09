# ICL-supported random-subtracted weak-lensing stacks (2026-11-05)

## Method and samples

Astrafek [`run_icl_supported_wl_subset_stack.py`](run_icl_supported_wl_subset_stack.py)
reuses the verified HSC S19A shape rows, calibrated tomographic estimator,
matched random centers, and paired endpoint jackknife from the full random-
subtracted stack. It reproduced the saved 800-lens cluster profile before
subsetting and verifies that each filter's IDs are nested in the exact
435-lens common target/random sample. The filters and selections are kept
separate.

| Selection | `i` lenses / shape pairs | `i2` lenses / shape pairs |
|---|---:|---:|
| ICL 200–300 kpc profile support | 38 / 194,044 | 70 / 327,773 |
| ICL plus saved 50–75 kpc satellite profile | 12 / 60,624 | 27 / 109,336 |

Each lens retains five matched random centers. The random side contributes
822,845 / 1,361,433 shape pairs to the ICL-supported `i` / `i2` stacks, and
235,676 / 538,909 pairs to the narrower joint satellite/WL stacks.

## Results and null checks

The tangential profiles are positive and decline with projected radius. On the
ICL-supported set, the seven-bin tangential zero-profile test gives
`χ²=307.3` (`i`) and `502.4` (`i2`) under the 5° covariance; under 2.5° it gives
`125.9` and `180.0`. These numbers are conditional on the jackknife
partitions, not publication detections.

The cross component provides a critical diagnostic. For `i` it rejects zero
for both the ICL-supported and joint satellite/WL samples under both block
widths:

| Selection | 5° cross-null p | 2.5° cross-null p |
|---|---:|---:|
| ICL-supported `i` (38 lenses) | `3.4×10⁻¹¹` | `8.0×10⁻⁴` |
| Joint ICL/satellite `i` (12 lenses) | `2.5×10⁻²⁷⁴` | `3.5×10⁻²²` |

This persistent cross signal is a systematic/covariance failure warning. The
`i` tangential profile must not be interpreted as a physical lensing detection
until its cause is isolated.

For `i2`, the cross component is consistent with zero at both widths:

| Selection | 5° cross-null p | 2.5° cross-null p |
|---|---:|---:|
| ICL-supported `i2` (70 lenses) | `0.310` | `0.421` |
| Joint ICL/satellite `i2` (27 lenses) | `0.282` | `0.575` |

The `i2` tangential zero-profile null is rejected conditionally for both
selections and widths. This is a promising exploratory mass-profile result,
but it still lacks validated contiguous HSC jackknife regions, full-aperture
random matching in depth/mask/PSF/sky, photo-z uncertainty, source-density boost
and reduced-shear corrections, and coadd covariance validation. None of these
subset stacks tests a continuous cluster-by-cluster ICL–satellite correlation;
the joint 12/27 sample is too small for that claim. Every output retains
`science_measurement=false`.

Tangential point estimates are identical across the 5° and 2.5° runs. The
median diagonal-error ratios (2.5°/5°) are 1.05 / 1.12 for the 38/70-lens
ICL-supported sample, and 0.99 / 1.21 for the 12/27-lens joint sample. The
block comparison is partial: equal-area cells can cross disconnected survey
regions, and two widths do not validate the covariance model.

## Astrafek artifacts

- Reproducible runner: [`run_icl_supported_wl_subset_stack.py`](run_icl_supported_wl_subset_stack.py).
- Full-covariance null audit and profile plot: [`analyze_icl_supported_wl_subset_stack.py`](analyze_icl_supported_wl_subset_stack.py).
- 5°/2.5° comparison: [`compare_icl_supported_wl_blockwidth.py`](compare_icl_supported_wl_blockwidth.py).
- ICL-supported stack JSON: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_wl_subset_stack_2026-11-05_01/summary.json` (SHA-256 `c88067c23f237dbad3b2c337f587988140f69fa49a4ac29a6de335107c703752`).
- ICL-supported 2.5° stack JSON: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_wl_subset_stack_2026-11-05_02p5/summary.json` (SHA-256 `23395661a3104b1b8c3945831bae96301a986137b82d3b6c7e9ffce9211a09df`).
- Joint satellite/WL 5° and 2.5° stack JSONs: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_wl_subset_stack_2026-11-05_01/summary.json` and `_02p5/summary.json` (SHA-256 `5a4112876829601c3222288467685bf37c31e8f091248e6f173cd14c5cfd24a5`, `b2178b6212f48c2c0956e1096b66d76f64063891833c310fd8f6afc59aa26f1d`).
- Both comparison JSONs and filter-separated figures are in the corresponding `*_blockwidth_comparison_*` and `*_assessment_*` output directories.
