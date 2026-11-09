# 17×17 covariance propagated through the ICL annuli

## Finding

The new ±8-pixel residual kernels were contracted with the exact valid
200–300 proper-kpc annulus masks for all 224 saved target/control pairs. The
reconstructed annular means and coverage matched the archived profile values
for every pair. Compared with the previous exact-aperture 9×9 calculation,
the independent-frame pixel-noise stack sigma rose only 0.8% in `i` and 1.3%
in `i2`:

| Filter | Pairs | Diagonal pixel-noise sigma | 9×9 sigma | 17×17 sigma | 17×17 / spatial 95% half-width |
|---|---:|---:|---:|---:|---:|
| `i` | 71 | 0.0627 | 0.0787 | 0.0793 | 7.6% |
| `i2` | 153 | 0.0601 | 0.0736 | 0.0745 | 7.0% |

Units for the three sigma columns are nJy arcsec⁻². The current stack means are
3.226 nJy arcsec⁻² in `i` and 3.787 nJy arcsec⁻² in `i2`; their saved spatial
block 95% intervals have half-widths 1.041 and 1.060 nJy arcsec⁻². Under the
independent-frame approximation, finite-lag pixel covariance contributes a
small part of that empirical spatial scatter. The increase from a 9×9 to a
17×17 kernel is also small for the actual annulus mean, despite visible
target/control differences in some per-pixel lag sums.

This is a finite-lag pixel-noise floor, not a total uncertainty estimate. The
17×17 kernels use patch medians from 32 patches per image, fewer than the 96
patches used in the earlier 9×9 run. The comparison is therefore a stability
check, not proof that the true covariance is converged. It omits source
Poisson noise, covariance beyond 8 pixels, sky transfer, PSF wings, mask
uncertainty, and spatial covariance among cluster sightlines. The calculation
does not adjust the profile or turn the current signal into a validated ICL
measurement. Keep `science_measurement=false`.

## Reproducibility

- Runner: `hsc_debug/analyze_icl_long_lag_aperture_covariance.py`
- Runner SHA-256: `973a787f5a8869772311231227addcb9c2d821f1e491d26b03d6e91e13869623`
- Output summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-06_01/summary.json`
- Summary SHA-256: `ebb78c5202ff60c5484b5b6aff69d61050bf57d11b005f548e10e8fdc6fb1be4`
- Frame aperture table SHA-256: `c6ada9453d33572caf72755a5ddee739830dd926dad1cf79988a7e3519d3e32b`
- Pair aperture table SHA-256: `658045f4513f2526a94797d16303a1103f018bed125dc8d01752a8bde6fadb8e`
- Input long-lag summary SHA-256: `c64191c18bdab3cfe5c330e5f385988e8a3c752f4b209eb562b9dfb391a14d6c`
- `hsc_network_used=false`, `data_requested=false`, `covariance_correction_applied=false`, `science_measurement=false`.
