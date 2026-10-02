# Background and noise model (plugin `noisemodel`, tab Measure, chip "Noise")

Measures what the other plugins assume: a background model that follows large-scale structure (gradients, scattered light), a local rms map, how correlated the pixel noise is
(drizzled / resampled images), and the resulting flux and magnitude errors. `sqrt(N) * sigma_pix` underestimates the error of any aperture larger than a few pixels
in a correlated image; this plugin measures the real law from blank apertures and applies it to the catalog.

## What is computed (`plugins/noisemodel/noisemodel.py`, engine `ogfkit/noise.py`, pure functions)

1. **Source mask** - `sep` detection at `mask-thresh` sigma (default 2.0, min. area 8, grown by 3 px) plus bad pixels (`{mask}` flags, zeros, NaN). Lower thresholds mask noise peaks of correlated noise and bias the sky low (measured: 1.5 sigma / 5 px masks 17 % of a pure correlated-noise image and gives sigma_pix 4.79 instead of 5.00).
2. **Background model** (`model`): `mesh` = clipped mesh (`bw` px cells, `fw` median filter); `poly` = polynomial surface of total degree `poly-order` fitted to block means with clipping; `mesh+poly` = polynomial first, mesh on the residual. Outputs `noisemodel_bkg.fits`, `noisemodel_rms.fits` (local pixel rms) and, with `correct`, `noisemodel_sub.fits`.
3. **Flatness**: block medians on a 6x6 grid before and after subtraction: peak-to-peak, chi2 of the block medians against a constant (expected ~1), fitted plane. A flat-field residual or an unremoved gradient shows as chi2 >> 1.
4. **Pixel autocorrelation** (FFT, masked pixels excluded): rho(1,0), rho(0,1), rho(1,1) and the FWHM of the equivalent Gaussian kernel.
5. **Noise law**: `n-aper` source-free circular apertures per radius (positions avoid the mask), sigma of their sums sigma_N(N); fit `sigma_N = sigma_pix * alpha * N^beta` (white noise: alpha 1, beta 0.5). The columns below use this law.
6. **Catalog columns**: `NM_SKY` (model at the object), `NM_RMS` (local pixel rms), `NM_NAPER` (aperture area), `NM_CORR` = sigma_N / (sigma_pix sqrt(N)) from the law, `NM_FLUXERR` = sqrt(sigma_N,local^2 + flux/gain), `NM_MAGERR`, `NM_SNR` (for `aperture = auto`: FLUX_AUTO and a Kron-like ellipse `kron-fact` x KRON_RADIUS x (A, B)). With `aperture = fixed` the circular-aperture flux on the background-subtracted image is measured too: `NM_FLUX_AP`, `NM_MAG_AP`, `NM_MAGERR_AP`. Metadata (`catalog_meta.json`, key `noisemodel`): sigma_pix, alpha, beta, rho1, kernel FWHM.

Windows: *Noise Curve / Autocorrelation Plot* (sigma_N vs N with sqrt(N) and the fitted law, ACF image, block medians), *Noise Summary*, *Show Background / RMS Frames*. Catalog keys: `noisemodel,bkg_file|rms_file|sub_file|json_file|plot_file|curve_file`.

CLI (also what the session script replays): `python3 plugins/noisemodel/noisemodel.py IMAGE --work DIR [--catalog TSV --meta-out JSON] [--mask FITS] [--model mesh|poly|mesh+poly] [--bw 64] [--radii 1,2,3,4,6,8] [--n-aper 500] [--aperture auto|fixed --aper-radius 3] [--gain G] [--correct]`. Without `--catalog` it prints the summary text. `noisemodel.analyse(data, bad, args, radii)` returns (JSON-serialisable summary, arrays) and is what a server would call.

## Validation (`plugins/noisemodel/tests/test_noisemodel.py`, 8 tests, ~7 s; truth known by construction)

| Test | Result |
| --- | --- |
| correlation recovered (Gaussian kernel 0.8 / 1.2 px, exact discrete ACF) | rho(1,0) 0.6709 (truth 0.6718), 0.8369 (0.8406); sigma_pix within 3 % of 5.00 |
| sigma_N(r) of blank apertures vs analytic sum over the exact ACF (r = 1.5, 3, 5, 8) | -2.1 %, -0.1 %, -2.1 %, -2.8 %; naive sqrt(N) is low by x2.1, 2.8, 3.1, 3.3 |
| white noise | law matches sigma sqrt(sum w^2) (fractional-pixel edge weights) within 8 %; beta 0.54, alpha 0.81 for r 1.5-8 (not exactly 0.5 / 1 because of edge pixels), 0.5 +- 0.07 for r >= 4 |
| flux-error pulls, 400 faint stars (flux 400, aperture r 6, correlated noise) | corrected pull std 1.14 (mean 0.16), naive sqrt(N) pull std 3.65 |
| gradient (+-2 counts) + scattered-light bump (6 counts, sigma 150 px) | rms error of the background model 0.25 counts (mesh), 0.27 (mesh+poly), 1.05 (quadratic only); block-median peak-to-peak 7.0 -> 0.64 (mesh); a plane is recovered to < 0.05 |
| flat-field residual (6 counts, 600-700 px scale) | flatness chi2 1033 vs 0.84 for a clean image; also reported by the CLI |
| non-uniform rms (4 -> 8) | rms map / truth 0.969 ... 0.994 |
| real drizzled HUDF F160W crop (800x800) | sigma_pix 4.4e-4, beta 0.71, alpha 1.20, rho(1,0) 0.48, kernel FWHM 1.4 px; naive errors too small by x1.3 (r 1 px), x2.1 (2), x2.7 (4), x3.6 (8) |

## Limits

* The law is a two-parameter power law fitted over the radii you give; for noise with a different correlation structure (e.g. strong large-scale 1/f pattern) it is only good inside the fitted range (extrapolating past the largest radius is not checked). Over r 1.5-8 on the synthetic image it follows the exact curve within 12 %.
* Blank apertures avoid masked pixels, so very crowded fields give few positions (< 30 positions: the radius is skipped). Sources fainter than the mask threshold still add a little to sigma_N.
* The noise is assumed stationary over the image (one alpha, beta); only sigma_pix varies via the rms map. Variable correlation (e.g. mosaics with differently drizzled tiles) is not modelled.
* Source masking at low threshold biases the sky low and sigma_pix slightly (see above); the default 2 sigma / 8 px keeps that below 0.1 counts (0.02 sigma) in the tests, but real extended faint emission is still absorbed into the mesh at `bw` px.
* `NM_FLUXERR` is a sky + Poisson error of a given aperture; it does not include background-model error (~0.2 counts per pixel in the tests, correlated over the mesh), PSF/aperture-loss errors or the error of the zero point. `gain` is only meaningful for images in electron-like units.
* The mean flux pull in the star test is +0.16 sigma (background low by ~0.06 counts/pixel from masked noise peaks).
