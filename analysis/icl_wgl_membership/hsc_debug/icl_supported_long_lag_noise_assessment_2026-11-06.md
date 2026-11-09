# Supported ICL sample: long-lag covariance audit

## Result

The residual covariance kernel was extended from ±4 to ±8 pixels for the 224
target/control cluster-filter pairs that support the saved 200–300 proper-kpc
ICL contrast: 71 pairs in `i` and 153 in `i2` (448 images total). This reused
only locally held HSC images, masks, variance planes, and metadata; it made no
HSC request. The HSC frame-metadata query approved under plan SHA-256
`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3` had
already completed all 16 batches (67,071 rows, no nulls in the eight requested
calibration fields). Those CCD-frame summaries still do not contain pixel
coadd weights or pixel-to-amplifier maps.

The approximate square-kernel variance factors at radius 4 and radius 8 were:

| Filter | Radius | Target median | Control median | Median paired target−control | 95% connected-block interval |
|---|---:|---:|---:|---:|---:|
| `i` | 4 px | 1.630 | 1.585 | 0.025 | −0.002 to 0.047 |
| `i` | 8 px | 1.655 | 1.559 | 0.081 | −0.001 to 0.132 |
| `i2` | 4 px | 1.626 | 1.564 | 0.059 | 0.026 to 0.090 |
| `i2` | 8 px | 1.661 | 1.611 | 0.077 | 0.030 to 0.129 |

The cumulative sums are not monotonic with radius because some covariance
shells are negative. This means the ±4 kernel is not demonstrated to be
converged. At radius 8, the target-control paired median difference is
positive in both filters, but its interval includes zero in `i`; the `i2`
interval excludes zero in this exploratory filter-specific analysis. There
are only 18 independent 5-degree blocks for `i` and 40 for `i2`, and the
filter-wise intervals are not multiplicity-adjusted. This supports treating
the observed target/control covariance mismatch as a systematic to model, not
as a corrected ICL signal. A follow-up contraction through the actual annulus
mask found that extending the exact-aperture calculation from 9×9 to 17×17
increased the independent-frame stack pixel-noise sigma by only 0.8% in `i`
and 1.3% in `i2`; see the [aperture propagation assessment](icl_long_lag_aperture_covariance_assessment_2026-11-06.md).

## Method and limits

For each image, the run measured median variance-normalized residual
covariances in deterministic, non-overlapping 64×64 patches (up to 32 patches
per image), after an affine-plane subtraction and 4-sigma clipping. The
aperture proxy sums a symmetrized stationary 17×17 median kernel over square
lag windows. Uncertainty uses 5,000 paired resamples of connected 5-degree
sky blocks. Pair support is exactly the saved 200–300 kpc set, and targets and
their matched controls remain paired.

This estimate is noisier than the earlier 96-patch calculation and a square
window only approximates radial-aperture covariance. It does not model spatial
variation in the covariance map, lags beyond 8 pixels, source Poisson noise,
coadd weights, sky transfer, PSF wings, or residual source masks. The result
therefore does not repair the ICL flux uncertainty or qualify the stack as a
science measurement. Keep `science_measurement=false`.

## Reproducibility

- Runner: `hsc_debug/audit_icl_supported_long_lag_noise.py`
- Runner SHA-256: `7ba16e7ce7ba49d0953666c646aa18d34422396bab5e2f13bf9e5787316023c1`
- Output summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_long_lag_noise_2026-11-06_01/summary.json`
- Summary SHA-256: `c64191c18bdab3cfe5c330e5f385988e8a3c752f4b209eb562b9dfb391a14d6c`
- Per-image table SHA-256: `9256abd57b0b41bbcd3636abc0dbd8d6dba028f1f6469138339725c887554096`
- Lag table SHA-256: `47e46f866d40901601472e0237d5c87946a3b00dc0c6b514db88c71d6b1475b1`
- Frame metadata audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_02/summary.json`
- Frame metadata audit SHA-256: `983458e302e1c27a2a4cb34243dd65a99c96c7d1cb00e31ef9d33093f7de0061`
