# ICL annulus covariance extended to 16 pixels

## Result

Extended the variance-normalized coadd residual-lag measurement from ±8 to
±16 pixels for the same 224 target/control pairs supporting the 200–300 proper-
kpc ICL annulus (71 `i`, 153 `i2`). All 448 target/control frames completed.
The shared lag cells from −8 to +8 pixels match the prior run exactly (129,472
cells; zero differences). The 33×33 lag table contains 487,872 frame-lag rows.

Contracting the new kernels with each pair's exact saved annulus mask gives:

| Filter | Pairs | ±8-pixel sigma (nJy arcsec⁻²) | ±16-pixel sigma | Change | Paired 5° block 95% interval for change |
|---|---:|---:|---:|---:|---:|
| `i` | 71 | 0.07931 | 0.06311 | −20.4% | −21.3% to −19.1% |
| `i2` | 153 | 0.07449 | 0.06290 | −15.6% | −21.5% to −9.4% |

The bootstrap resamples the 18 `i` and 40 `i2` connected 5° blocks for 5,000
replicates. These intervals are conditional on the estimated kernels; they do
not propagate patch-to-patch uncertainty in each lag covariance. The extended
finite-lag sigma is 6.1% of the saved spatial-interval half-width in `i` and
5.9% in `i2` under the existing independent-frame approximation.

## Interpretation and limits

The square-window cumulative lag sum falls sharply between radii 8 and 16
pixels: target/control medians change from 1.656/1.559 to 0.825/0.755 in `i`,
and 1.661/1.611 to 0.827/0.740 in `i2`. Negative long-lag shells therefore
dominate the 16-pixel change. The per-lag estimator takes a median across 32
patches after fitting and subtracting a local affine plane; assembling those
separate medians does not guarantee a positive-semidefinite covariance kernel.
The decrease is consequently an estimator-sensitivity result, not evidence
that the physical image noise is lower or that a correction should be applied.

The actual-mask aperture contraction stayed positive and reproduced the saved
annular means and coverages for all 224 pairs. Even so, stationarity, patch
medians, affine-plane removal, unmodeled lags beyond 16 pixels, source Poisson,
sky transfer, PSF wings, and cross-frame field covariance limit its use. No
profile was corrected and `science_measurement=false` remains.

## Reproduction

- Kernel runner: `audit_icl_supported_covariance_to_16px.py` (SHA-256
  `2fc70a3bd0c8f1b64e4cb4338c8f4c1019120081fa20a387078966968aead948`).
- Exact-aperture runner: `analyze_icl_covariance_to_16px_aperture.py` (SHA-256
  `5fee4230b19d4451cc59d8cb764f69cb8478aa2da4046fb994659f814f15f014`).
- Block-bootstrap runner: `assess_icl_covariance_16px_aperture_blocks.py`
  (SHA-256 `04d4097dfaa8095e5a27f385b508795c4412660de9c02b1af238a5ad15100e38`).
- 33×33 lag summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_covariance_to_16px_2026-11-06_01/summary.json`
  (SHA-256 `7fcc002718cbac88c04fac0fb6b00b02b176f9f4813d75f8d4f4ae03403be80a`).
- Exact-aperture summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_16px_2026-11-06_01/summary.json`
  (SHA-256 `474489c8b052188753dc4b70c954ac1b983a4ccd60ab25475bf814158922bf8c`).
- Paired block-bootstrap summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_16px_block_bootstrap_2026-11-06_01/summary.json`
  (SHA-256 `8c2acc5266bee2497418ab3ebfd13141a789c164a9d8d2b9086ffac07770344a`).

The runs used only existing local HSC coadds, masks, profiles, and warp-audit
outputs. No HSC request, download, GPU, or Slurm job was used.
