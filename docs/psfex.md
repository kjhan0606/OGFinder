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

## Extensions (precision batch item 2): reduced rank, wings, extended PSF, error budget, stacked-star PSF
Code: `ogfkit/psfext.py`; CLI options of `plugins/psfex/psfex.py` (GUI group *Extensions*); tests `plugins/psfex/tests/test_psfext.py` (5); validation `plugins/psfex/validation/psf_validate.py synth|real`, reports `psf_synth_report.json`, `psf_real_report.json`.
* **Reduced rank** (`--rank K`): the polynomial variation of the model (coefficient images of the x/y terms) is projected onto its K leading singular images, giving a smoother model with fewer degrees of freedom. The result is an ordinary `PSFModel` (`meta.rank`, `meta.rank_variance_kept`), so all consumers work unchanged. The step is post-hoc (SVD of the fitted coefficients), not a constrained fit.
* **Wings** (`--wings`, `--ext-size N`): azimuthal profiles of bright isolated stars beyond the core with the sky *given* (a ring sky would absorb the halo), a power-law fit I = A (r/r0)^-gamma (accepted only for gamma > 2.1, else `valid=false`), the flux beyond the stamp, and an extended PSF stamp = model core joined to the power-law tail (`psfex_wings.json`, `psfex_extended.fits`, `info.tail_fraction`). Use the FITS through `--psf` of multifit/daophot.
* **Error budget** (`--assess`, `psfex_budget.json`): for Sersic profiles (n, R_e) the maximum |model difference|/peak, flux change and the bias of a free (mag, R_e, n) fit when the PSF used differs from the true one; the file reports the model evaluated at the field centre used at a corner and vice versa. Not propagated into catalogues per object.
* **Stacked-star PSF in the Stacking plugin** (`--psf-stamp`): median stack of the selected sources with recentring, local sky, negatives clipped, azimuthal-median wings beyond 4 px (otherwise noise gives negative wings and spurious chi2 differences, seen in the GALFIT comparison) -> `stacking_psf.fits`.

### Numbers
Synthetic (spatially varying Moffat FWHM 3 -> 4 px, e 0.02 -> 0.14, 150 stars + 15 % star-like galaxies + 10 % blends, 2 seeds), median max|model - truth|/peak: degree 0 0.150, degree 1 0.044, degree 2/3 0.041, degree 3 rank 1 0.040, rank 2 0.042 (rank reduction costs nothing here; the polynomial variation cuts the error 3.7x). Star selection used 81 % of the true stars, 2/33 star-like galaxies, 0/1 blends. Halo test (8 % second Moffat, true sky supplied): gamma 3.07, rms 0.05 dex, flux beyond 8 px 0.125 fitted vs 0.111 true (+13 %). The extended PSF does *not* beat the core-only PSF everywhere: enclosed flux at R = 6 px 0.814 vs truth 0.852 (core-only 0.899), at R = 55 px 0.981 vs 0.9988 (a pure power law is flatter than the Moffat tail), but it recovers the halo beyond 20 px that the core stamp misses. Error budget with big galaxies (R_e 8 px, true halo PSF): core-only 31 px stamp max|d|/peak 0.040 (n=1) / 0.058 (n=4), extended 0.012 / 0.013; but the fitted parameter biases are not clearly smaller (R_e +3.1 %/+9.6 % -> -1.5 %/-6.7 %, n -2 %/-6.5 %): the image error shrinks 3-4x, the parameter bias only changes. For a wider-than-true PSF (+15 % FWHM) an n=1 galaxy gets R_e -9 % (n = 4: +2 % through the n degeneracy, n +62 %).
Real data (held-out stars, |residual|/flux median; fit on half the stars, test on the other half): HUDF F160W 2000^2 crop, 172 stars: degree 0 0.350, degree 1 0.243, degree 2/3 0.252 (degree capped by the star number), rank 1 0.250, rank 2 0.248 - a 30 % gain from spatial variation, none from rank. M51 (46 stars, undersampled, saturated cores): ~0.80 for every model (no gain), wings not measurable (gamma < 2, flagged invalid). The HUDF wings fit gamma 4.3 beyond 12 px (negligible wing, 5.9 % of the core flux beyond the 31 px stamp). No ground-based halo data was fitted: the halo recovery is verified only on synthetic data.

### Limitations / not done
Rank selection is by hand (`--rank`), no automatic criterion; local-sky builders absorb part of the halo and the wings are noise-limited beyond ~25 px; no iterative star re-selection beyond the existing `reject_sigma` / `n_iter` (these already reject ~94 % of star-like galaxies on synthetic data); no colour or chip-to-chip dependence; PSF-error propagation is a diagnostic table, not per-object error inflation.

## Limitations
Stars are selected on a stellar locus (a field with < ~4 usable stars gets the fallback); the polynomial order is capped by the star count so the corners of a sparse field are
extrapolated (corner FWHM error up to 4-8 % with 22 stars, order 1); no colour dependence, no chip-to-chip models in one call (the Moving step builds one per chip); the model is built
from one image (no multi-exposure stacking of the star table).
