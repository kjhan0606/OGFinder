# PSF model plugin (`psfex`, PSFEx-like)

Spatially varying PSF model from the field stars, usable by every PSF consumer of the program.

## Model
`ogfkit/psfmodel.py` (`PSFModel`, pure functions, JSON/FITS serialisable):

    PSF(x, y) = sum_k phi_k(x, y) C_k     (phi_k: polynomial basis of order 0-3 in the normalised image coordinates)

* `C_k` are pixel-grid coefficient images (sampling 1/oversample px); a stamp at a sub-pixel offset is made by cubic-convolution interpolation, every stamp has unit sum.
* **Star selection**: DAOFIND-like candidates, isolation (neighbour distance ≥ 2.5 FWHM), peak S/N ≥ `snr-min`, no saturation / masked pixels, the stellar locus in
  (flux, half-light radius) is selected by sigma clipping, up to `max-stars` brightest.
* **Order** is reduced automatically with few stars: < 20 stars → constant, < 60 → linear, < 150 → quadratic (`order` is the maximum).
* **Kinds**: `empirical` (pixel table only) or an analytic base (`gaussian`, `moffat`, `lorentz`, `penny`) with the spatially varying empirical residual table on top (DAOPHOT style).
* **Iteration**: the stars are refitted with the model (flux, sub-pixel offset), outliers (chi) are rejected, the coefficients are re-solved; with `neighbour-iter > 0` the other sources
  are fitted with the model (ALLSTAR group fit) and subtracted from the star stamps instead of being masked.
* **Oversampling** `auto`: 2 only for undersampled PSFs (FWHM prior < 2.5 px, ≥ 25 stars; strong smoothing), else 1 — measured on synthetic fields, s = 2 with a well sampled PSF
  is worse than s = 1 (FWHM bias +4 % vs −2 %, ellipticity error 0.03 vs 0.007) because the sub-pixel phase fit is noisier.
* **Fallbacks**: 4-19 stars constant empirical, 1-3 stars a Moffat fit (`moffat_fit`), 0 stars a Gaussian of the prior FWHM (`gaussian_prior`); `model.meta['mode']` tells which.
* `model.at(x, y)` (centre stamp) makes the model a drop-in replacement for the per-tile `PSFField` of the tiled ZOGY.

## GUI (Measure tab, "PSF model" chip) and CLI
Steps: **Build PSF Model** (adds `PSFM_FWHM`, `PSFM_E`, `PSFM_PA`, `PSFM_NSTAR` = model FWHM / ellipticity / PA at every catalog row; records `analysis.psfex`), **Use Model as
Catalog PSF** (sets `psf,file` to the central stamp and `psf,model` to the model JSON), **PSF Maps Plot…**, **PSF Star Table…**, **PSF Map Table…**.
Files in `<work>/psfex/`: `psfex_model.json/.fits`, `psfex_center.fits`, `psfex_stars.tsv`, `psfex_maps.tsv`, `psfex_info.json`, `psfex_diag.png`. Catalog keys `psfex,*`.

    python plugins/psfex/psfex.py IMAGE --work DIR [--catalog TSV] [--mask FITS] [--kind empirical|gaussian|moffat|lorentz|penny] [--order 2]
        [--oversample auto|1|2|3] [--size 0] [--snr-min 20] [--max-stars 300] [--saturation 0] [--fwhm 3] [--neighbour-iter 0] [--star-list TSV] [--grid 7]

## Consumers
* **PSF photometry / crowded photometry** (Photometry plugin): when `psf,model` is set (after *Use Model as Catalog PSF*) the steps append `--psf-model FILE`; every source (crowded:
  every group) uses the model PSF at its position. The argv is unchanged when the key is unset.
* **Tiled ZOGY** (Moving plugin): parameter *PSF source* = `model` (with *PSF tile* > 0, CLI `--psf-source model`) builds one model per chip and evaluates it at every tile centre; a model
  narrower than the hot-pixel floor or built from < 4 stars falls back to the default PSF. Default `tiles` (unchanged argv).
* DAOPHOT plugin (own model), completeness injection (`PSFModel` accepted), and the multi-component fitting plugin use the same class.

## Fix made on the way
`psf_phot/photometry.py`: the single-source fit compared the PSF stamp with the *top-left corner* of the (larger) cutout, so the fitted flux was ≈ 0 for any star. The stamp is now centred
in the cutout (recovered flux −1 %, scatter 10 %). argv and columns are unchanged.

## Validation (synthetic 600² field, Moffat PSF FWHM 2.8→3.8 px along x, e 0.04→0.14 along y; `plugins/psfex/tests`)
| test | result |
|---|---|
| model vs truth at 15 positions (90 stars, order 2) | FWHM rel. error rms 0.017 (max 0.040, bias −0.017), ellipticity error rms 0.007 (max 0.017) |
| order 0 vs order 2 FWHM error rms | 0.095 vs 0.017 |
| 4 stars / no stars | constant empirical model / Gaussian prior, valid unit-sum stamps |
| PSF photometry, 245 isolated stars mag 16-19.5 | constant PSF: bias −0.012, rms 0.099; model: bias −0.011, rms 0.026 |
| crowded photometry, 92 stars in a crowded 300² field | constant: median −0.081, rms 0.115; model: −0.079, rms 0.059 (the legacy fitter's −8 % bias is unchanged) |
| tiled ZOGY, seeing FWHM 2.4→4.2 vs reference 2.2→2.6, 20 injected transients | residual at bright stars (rms fraction of flux): 0.121 constant PSF → 0.030 model; transient flux ratio median 0.949 → 0.932, scatter 0.123 → 0.096 |

## Limitations
Stars are selected on a stellar locus (a field with < ~4 usable stars gets the fallback); the polynomial order is capped by the star count so the corners of a sparse field are
extrapolated (corner FWHM error up to 4-8 % with 22 stars, order 1); no colour dependence, no chip-to-chip models in one call (the Moving step builds one per chip); the model is built
from one image (no multi-exposure stacking of the star table).
