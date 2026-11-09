# HSC ICL warp pilot: coadd VARIANCE consistency assessment (2026-11-06)

## Scope and provenance

The user explicitly approved one pilot request: target cluster 6966, HSC PDR3 Wide, HSC-I, tract 9703 / patch 306, 2 arcmin square, and the image, mask, and variance planes for the nine visits listed in the pinned request. The approved plan SHA-256 is `a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e`. No additional positions or products were downloaded.

The archive is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_download_2026-11-06_01/warp_pilot.tar`, 55,511,040 bytes (52.92 MiB), SHA-256 `cdc9b642c5c5d0361b58a0bf7c0e733a60cb78caad5ed496ae00cac89b23953c`. The offline FITS audit found exactly the nine approved visits. Their image/mask/variance arrays have consistent 715×715 shapes and aligned WCS; the maximum warp-to-coadd pixel-grid residual is below `7e-10` pixel and all warp footprints map to coadd offset (661, 661). Two visits (77296 and 77304) have zero finite image/variance coverage; visit 77330 covers 59.5%; remaining visits cover 86.9–100%.

## Method

`analyze_icl_warp_pilot_coadd_variance.py` reads only the pinned local archive, target coadd, sample manifest, and approval plan. For each coadd pixel it predicts the variance as `1 / sum(1 / VARIANCE_visit)` for finite, positive warp variances not flagged `MP_NO_DATA`; it then compares the delivered coadd VARIANCE to that prediction. The primary profile requires at least three valid visits. A mask-clean sensitivity additionally rejects the mask planes shared by all warps and the coadd: `MP_BAD`, `MP_BRIGHT_OBJECT`, `MP_CR`, `MP_CROSSTALK`, `MP_EDGE`, `MP_INTRP`, `MP_NO_DATA`, `MP_SAT`, `MP_SUSPECT`, and `MP_UNMASKEDNAN`.

The annuli use the pinned BCG spectroscopic redshift 0.46215, flat cosmology H0=70 km s⁻¹ Mpc⁻¹ and Ωm=0.3, and 0.168 arcsec pixels. The 200–300 proper-kpc annulus spans 34.21–51.31 arcsec. Pixel quantiles below are descriptive: neighboring pixels are correlated and these ranges are not confidence intervals.

## Result

Within the 200–300 kpc annulus, 162,352 pixels (99.71% of its geometric area) meet the primary requirement. The coadd/predicted-variance ratio has q16/median/q84 `1.00372 / 1.00748 / 1.01506`; 67.84% of those pixels lie within 1% of unity and 87.63% within 2%. The mask-clean sensitivity retains 79,509 pixels (48.83% of annulus area), with ratio q16/median/q84 `1.00312 / 1.00656 / 1.04093`; 72.55% lie within 1% of unity.

Across all five annuli from 0–400 kpc, the primary median ratio ranges from 1.00630 to 1.00748. In the warp footprint, median valid-visit count is six, 93.85% of pixels have at least six valid visits, and 99.89% overlap finite coadd pixels.

## Interpretation and limitations

The delivered coadd VARIANCE is locally consistent to roughly 0.6–0.8% in its median with an inverse-variance combination of the delivered visit VARIANCE planes. This supports the internal bookkeeping of the single downloaded coadd position. It does not independently calibrate detector noise or demonstrate the noise model for the 800-cluster sample.

The comparison assumes unit relative weights in the variance-only reconstruction. It does not measure source-Poisson variance, effective gain, correlated resampling noise, sky-subtraction transfer/covariance, extended-PSF residuals, or covariance of the ICL estimator. It is one cluster and one tract/patch; the pixel quantiles are not spatially independent. The strict mask cut also removes over half of the science annulus. These results do not support an ICL flux correction or a science claim. Keep `science_measurement=false`.

## Reproduction and outputs

- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_warp_pilot_coadd_variance.py` (SHA-256 `dc67da5b7e6c88a6202f3d07de67e52bd7765cafe34339e3b0eb705ad48f3695`)
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_coadd_variance_2026-11-06_02/`
- Summary SHA-256: `ad523613f9efdd5226759d13bef1973cb2add1990614a80b96b27adca0753dae`
- Radial TSV SHA-256: `c15e9cea8dc9004502e714da83eca5086f2d7856d17dee4d2e57dd73603ed785`
- No HSC request occurred during the analysis; no GPU or Slurm job was used.

An earlier same-day output ending `_01` was generated during development before the strict-mask overlap and footprint denominator were corrected. It is superseded; use `_02` as the recorded result.
