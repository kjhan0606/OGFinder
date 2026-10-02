# Isophote ellipse fitting (plugin `isophote`)

IRAF-`ellipse`-style isophote fitting for the galaxies of the catalog.  Engine: `photutils.isophote` (Jedrzejewski 1987 algorithm);
OGFinder adds the catalog driver, mask handling, derived quantities, output tables, the frames and the plot window.
Measure tab, chip **Isophote** (`[Isophote v] [>] [gear]`).

## Steps

| step | what |
|---|---|
| **Fit Isophotes** (`isophote.fit`, recorder name `analysis.isophote`, AUTO) | CLI step `plugins/isophote/isophote.py`; adds columns to the catalog; writes `~/.ds9/isophote_{model,resid}.fits`, `isophote_profiles.tsv/json`, `isophote_plot.png`; the model and the residual open in two new frames (parameter `show-frames`), the original frame stays current |
| Profile Plot... | window with the 6-panel PNG: surface brightness vs `sma^0.25`, ellipticity, PA, B4, growth curve, centre drift |
| Profile Table... | the long-form table (one row per isophote and galaxy) in a text window |
| Show Model and Residual | reopen the frames of the last run |

## What is measured

Per isophote (`isophote_profiles.tsv`; PA in degrees counter-clockwise from +x like `THETA_IMAGE`, centres 1-based): `sma`, `intens` (+err),
`mu` (+err, `zp - 2.5 log10(I / pixel_scale^2)`), `growth_flux`, `ellipticity` (+err), `pa` (+err), `x0`, `y0`, `a3 b3 a4 b4` (+err), `tflux_e`, `npix_e`,
`rms`, `n_data`, `n_flag`, `stop_code`.

* **Growth curve.** `growth_flux` integrates the fitted mean intensity over the elliptical annuli (area `pi a^2 (1-e)`).  Masked pixels do not bias it,
  because the isophote intensity is the mean of the *unmasked* pixels; photutils' own `tflux_e` (a sum of the unmasked pixels) is kept and *is* biased low
  where pixels are masked.  `ISO_R50`, `ISO_R80` and `ISO_MAG_TOT` use `growth_flux`; **they are relative to the flux inside the last isophote**
  (no extrapolation to infinity), so `ISO_MAG_TOT` is a truncated magnitude when the outer isophote is shallower than the galaxy.
* **Harmonics.** `B4 > +0.01` = disky, `< -0.01` = boxy (`ISO_SHAPE`).  In photutils' convention the cos(4θ) term that carries the boxy/disky signal is **B4**
  (checked on generalised-ellipse models, see tests); `A4`, `A3`, `B3` are reported as well.
* **Catalog columns**: `ISO_NISO ISO_RMAX ISO_X0 ISO_Y0 ISO_R50 ISO_R80 ISO_MAG_TOT ISO_EPS_HL ISO_PA_HL ISO_EPS_OUT ISO_PA_OUT ISO_B4_HL ISO_A4_HL ISO_SHAPE ISO_RCONV ISO_FLAG`
  (`_HL` = mean over the *converged* isophotes (photutils stop code 0-2) at 0.8-1.25 x R50, widened to 0.5-2 x R50 once; empty when none converged;
  `ISO_RCONV` = largest converged SMA; `ISO_FLAG` 0 ok, 1 fewer than 6 isophotes, 2 failed).
* **Model / residual.** `build_ellipse_model` (with the harmonic terms unless `no-harmonics`), added up over the fitted galaxies; residual = data - model.

## Mask awareness

* The shared mask manager's effective mask (`{mask}` token = `~/.ds9/mask_<base>_bool.fits`) is used automatically when a mask exists; a bit-flag file
  (values > 1) is reduced with the manager's rule `(flags & 23) != 0 and (flags & 8) == 0`.
* The target itself is never masked (a radius of max(2 px, A/2) around its centre is cleared - the auto mask contains the target as a detected source).
* `mask-neighbours` (default 2.5): the other catalog objects are masked inside `factor x (A_IMAGE, B_IMAGE, THETA_IMAGE)` ellipses.
* Non-finite pixels are masked.

## CLI (also what the exported session script runs)

    python plugins/isophote/isophote.py IMAGE --catalog CAT.tsv [--numbers 1,5-9 | --max-objects N] [--mask M.fits] [--mag-zeropoint 25] [--pixel-scale 1]
        [--maxsma-scale 3] [--step 0.1] [--linear] [--fix-center] [--no-harmonics] [--mask-neighbours 2.5] [--n-workers N]
        [--model-out F] [--resid-out F] [--table-out F] [--json-out F] [--plot-out F]
    python plugins/isophote/isophote.py IMAGE --xy X,Y [--sma0 --eps0 --pa0 --maxsma]        # one galaxy, no catalog

Python API (server-ready: pure functions, JSON-serialisable results): `fit_galaxy(data, x0, y0, mask=...) -> (result_dict, model_array)`,
`fit_source`, `run_catalog`, `summarize`, `growth_curve`, `growth_radii`.  A failed start is retried with two other (sma0, e, PA) starts.

## Validation (measured)

Synthetic Sersic, 301x301, Gaussian noise sigma 0.05 (peak/noise 2e2-4e5), `plugins/isophote/tests/test_isophote.py`:

| test | result |
|---|---|
| n=1 and n=4, e=0.30, PA=40 deg, re=15 px, fit between 0.4 and 3 re | median e within 0.01 of the input, max deviation < 0.03, median PA within 1 deg, centre within 0.1 px, intensity within 3 % of the analytic profile at 0.3, 1, 2 re |
| growth curve, exponential | `growth_flux` and `tflux_e` within 2 % of the analytic flux inside 3 re; r50 within 5 % and total magnitude within 0.03 mag of the analytic truncated values |
| boxy (c=+0.5) / disky (c=-0.4) / elliptical generalised ellipses | classified boxy (B4 < -0.01) / disky (B4 > +0.01) / elliptical |
| extended companion on the major axis (amplitude 1.5, sigma 5 px) | unmasked: max abs(dI/I) 0.39, max abs(de) 0.066; masked: 0.020 and 0.0021; mask-aware `growth_flux` within 3 % of the analytic, photutils `tflux_e` biased low |
| residual, annulus 4 < r < 60 px | mean < 0.02, rms < 0.18 (noise 0.05; the n=4 core inside 4 px is not reproduced by a sampled isophote model) |
| CLI, two galaxies in one image (catalog + bit-flag mask) | e, PA of both within 0.02 / 1.5 deg, model + residual = image to 1e-4 |

Real data (`plugins/isophote/tests/validate_real.py`, HUDF F160W crop 1000x1000, SEP detections; independent reference = SEP elliptical aperture with the
same ellipse as the last converged isophote): growth-curve flux vs SEP, 5 galaxies: median dF/F +0.002; four of them within 2.2 %; the fifth
(-45 %) has a bright neighbour inside the aperture which the isophote fit masks and SEP sums (expected).  Ellipticity at R50 vs second-moment
ellipticity of SEP: median difference 0.000, rms 0.09; PA differences are meaningless for e < 0.2 (rms 40 deg overall; 2 deg for the e > 0.3 galaxies).  M51: the three catalog
blobs of the core region give 2 of 3 fits (the third fails: saturated core), photutils flags the outermost isophotes stop code 4 (geometry fixed), which
is why the `_HL` columns use converged isophotes only.

## Limitations

* No extrapolation to infinity; `ISO_MAG_TOT` is truncated at the last isophote.
* Overlapping galaxies are fitted separately; their models are *added* in the model image, so a shared region is double counted.
* The isophote model is a sampled interpolation: steep cusps (n=4 inside ~4 px) are underestimated; cmodel residuals there are not noise.
* Time: 5-15 s per galaxy (photutils is pure Python); use `n-workers` for several galaxies.
* Edge: isophotes that leave the image stop the fit (large `maxsma-scale` on a small image).
