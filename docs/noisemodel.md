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

## Unified photometric error model (`ogfkit/photerr.py`, precision batch item 3)
Several tools estimated aperture errors separately (noise plugin: sky + correlation law; DAOPHOT `phot`: `N sigma^2 (1 + N/nsky)`; stacking: bootstrap; depth maps: own copy of the blank-aperture sampler). `photerr.aperture_error` is now the one place where the terms are added:

| term | meaning |
|---|---|
| sky | `sigma_pix sqrt(N)` times the correlation factor of the blank-aperture law (`alpha N^(beta-1/2)`), or the measured law of the *local-sky-subtracted* aperture (`NoiseModel.add_local`, includes the next term) |
| poisson | `flux / gain` |
| sky_est | error of the annulus sky level times N (`N^2 sigma^2 kappa^2 / nsky`; kappa from the law or measured with `annulus_factor`: median of white noise 1.25, correlated noise 3-5) |
| sky_sys | per-pixel background-model error times N |
| contam | `neighbour_contamination`: flux the PSF wings of the other sources put into the aperture, net of the over-subtracted annulus sky (a bias, positive or negative) and its uncertainty (`c_ij sigma_fj` in quadrature + 10 % of the contamination) |
| apcorr | `aperture_fraction` (net of the source's own wings in the sky annulus), `frac_err` of the correction enters `total = f_ap/frac` |

Replaced duplicates: `noise.flux_error` (wrapper, identical numbers), `depth.blank_circle_positions` (now `noise.blank_positions(box=...)`, bit-identical positions), `daophot.phot` variance (identical numbers without a noise model, tested against the old formula). New optional inputs: Noise plugin parameters *Local sky annulus* (`--sky-annulus 12,18`: columns `NM_SKYERR`, `NM_FLUXERR_LOC`, and the measured local law is stored in `catalog_meta.json`) and *PSF for neighbour contamination* (`--psf`: `NM_CONTAM`, `NM_FLUXERR_TOT`); DAOPHOT PHOT reads the meta (`--noise-meta`, used only when the stored annulus/sky mode match). Defaults leave all existing columns unchanged. Not replaced: stacking bootstrap errors (they measure source-to-source scatter, a different quantity), multifit covariance errors, moving-object detection S/N (`sqrt(npix) rms`, uncorrelated; changing it would move detection thresholds).

### Validation (`validation/photerr_validate.py`, `photerr_report.json`; tests `tests/test_photerr.py`, 7)
Blank apertures with local annulus sky (12-18 px, median), positions independent of those used for the law; pull = flux/error, ideal std 1:
| image | r (px) | naive `sigma sqrt(N)` | legacy DAOPHOT | global law + annulus | measured local law |
|---|---|---|---|---|---|
| synthetic, kernel 1 px | 2 / 3 / 5 / 8 | 2.5 / 3.0 / 3.4 / 3.7 | 2.5 / 2.9 / 3.2 / 3.2 | 1.05 / 1.11 / 1.09 / 0.99 | 1.06 / 1.07 / 1.02 / 0.94 |
| real HUDF F160W (1400^2) | 2 / 3 / 5 / 8 | 1.7 / 2.0 / 2.3 / 2.7 | 1.7 / 1.9 / 2.2 / 2.3 | 0.72 / 0.64 / 0.54 / 0.46 | 1.00 / 1.00 / 0.98 / 0.97 |
On the real drizzled image the global-sky law *over*-estimates the local-sky error by 1.4-2.2x (large-scale noise is removed by the local sky); the measured local law is correct to 3 %. Use `--sky-annulus` whenever the photometry uses a local sky.
Injection-recovery (Moffat FWHM 3 px stars, flux 100-10^4 x sigma/5, source Poisson noise, aperture r 4, annulus 8-14, aperture correction from a stacked-star PSF with bootstrap error; pull of the total-flux estimate, std / mean):
| sample | naive | legacy | unified, no contamination | unified + contamination |
|---|---|---|---|---|
| synthetic, isolated (n 138 / 281 / 316 in sparse / moderate / crowded fields) | 3.4 / 3.6 / 3.3 | 2.7 / 2.7 / 2.8 | 1.07 / 1.07 / 0.95 | same |
| synthetic, a neighbour within 22 px (n 10 / 127 / 804) | 1.8 / 14.8 / 18.9 | 1.5 / 11.1 / 13.3 | 0.57 / 4.6 / 4.3 | 0.58 / 1.84 / 1.97 (mean pull +0.4 / +1.0 / +0.75) |
The contamination correction brings the pull of close pairs from 4.3-4.6 to 1.8-2.0, not to 1: neighbours below the detection cut (S/N 3) and the 10 % PSF-model uncertainty are not modelled, and the neighbours' fluxes are themselves measured. Real HUDF injection (inside the real, source-masked image; only 20 isolated + 52 close samples): pull std 0.65-0.96 isolated, 1.8 close pairs with +0.9...+1.4 mean pull; there the aperture correction from 40 injected stars measured in the real background was off by 2-8 % (0.757-0.803 vs 0.820 true), which dominates the bright-star pulls - the aperture-correction term is only as good as its `frac_err`. Verified on synthetic data for the contamination; on real data only the blank-aperture test is solid.

### Limits
Neighbours are point sources of the catalogue flux (extended neighbours' wings are wrong); the annulus factor/local law is stationary over the image; the local law is only valid for apertures within the measured radii (outside half/twice the area range `NM_*_LOC` is left empty); the Poisson term assumes a known gain; no covariance between aperture and annulus pixels; the 10 % PSF-model uncertainty of the contamination is an assumption, not measured.
