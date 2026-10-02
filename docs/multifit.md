# Multi-component fitting (`multifit`, GALFIT-like)

PSF-convolved model fitting of catalog objects with several components, a shared sky and simultaneous fitting of neighbours.
Engine `ogfkit/multifit.py` (pure functions, JSON-serialisable inputs/outputs), CLI/driver `plugins/multifit/multifit.py`, GUI chip "Multi-fit" (Measure tab).

## Model
* Components: `sersic` (free n), `exp` (n=1), `dev` (n=4), `psf` (point source). Parameters: centre (x, y), total flux (reported as magnitude with the zero point), effective radius `re`
  (major axis, pixels), `n`, axis ratio `q`, position angle `pa` (degrees, counter-clockwise from +x). Sersic images use the sub-pixel integrated renderer of `ogfkit.models`
  (finer grid in the cusp), are convolved with the PSF (FFT) and summed with the sky.
* **PSF**: `PSFModel` from the PSF-model plugin (`psf,model`; the stamp at the position of the component, i.e. the spatially varying PSF), a PSF image (`psf,file`), a Gaussian
  (`PSF FWHM fallback`) or - by default - a model built from the field stars of the image.
* **Sky**: one sky shared by all components: fitted constant, constant + gradient plane, or fixed.
* **Noise**: sigma = background rms (measured, or given), plus with `gain` the Poisson noise of the model. Pixels are weighted, masked pixels (shared mask `{mask}`, neighbour masks) ignored.
* **Fit**: bounded trust-region least squares (`scipy.optimize.least_squares`), parameter scaling, errors from the weighted Jacobian (`(JᵀJ)⁻¹`, valid because the weights are the true noise), `chi2_red`,
  BIC. Parameters can be **fixed**, **bounded** (defaults: re 0.3-…, n 0.3-8, q 0.05-1) and **tied** (common centre of bulge and disk).
* **Presets** for catalog objects: `sersic`, `exp`, `dev`, `psf`, `bulge+disk` (dev + exp, common centre, B/T), `psf+sersic` (AGN + host, common centre) and `auto` (lowest BIC of
  psf / sersic / bulge+disk, a more complex model must win by the BIC margin).
* **Neighbours**: `fit` (the brightest neighbours inside the cutout are added as free PSF / Sersic components, the others masked), `mask` (elliptical masks of `neighbour-radius` x A_IMAGE x B_IMAGE),
  `ignore`.

## GUI and CLI
Steps: **Fit Components** (all / the brightest `max-objects` rows; adds `GF_X GF_Y GF_MAG GF_MAGERR GF_RE GF_REERR GF_N GF_NERR GF_Q GF_PA GF_BT GF_MAG2 GF_SKY GF_CHI2 GF_NCOMP GF_NNEIGH GF_RESFRAC GF_FLAG`;
`GF_MAG` is the total of the target's components, `GF_RE/N/Q/PA` those of its brightest extended component, `GF_BT` the flux fraction of the first component of a two-component model),
**Fit Selected Rows** (NUMBERs of the catalog selection go into the recorded parameter `objects`), **Fit From Config File** (text output), **Show Model / Residual Frames**, **Montage…**, **Component Table…**.
Flag bits: 1 parameter at a bound, 2 not converged, 4 chi2_red > limit, 8 neighbours fitted, 16 cutout cut by the image edge, 32 > 50 % masked, 64 singular covariance, 1024 no usable start values, 2048 failed.
Files in `<work>/multifit/`: `multifit_results.tsv` (one row per component), `multifit_model.fits` / `multifit_residual.fits` (full frame, sky not removed), `multifit_montage.png`, `multifit_psf.json`.
Catalog metadata key `multifit`. Catalog keys `multifit,*`.

    python plugins/multifit/multifit.py IMAGE --catalog TSV --work DIR [--mask FITS] [--psf-model JSON|--psf FITS|--psf-fwhm PX] [--model auto] [--neighbours fit] [--objects "3;7"] ...
    python plugins/multifit/multifit.py IMAGE --config cfg.json --work DIR

Config file (explicit components, 1-based image coordinates unless `"one_based": false`):

    {"bbox": [x0, x1, y0, y1], "sky": "const|plane|fixed", "sky_value": 0, "gain": 0, "rms": 1.2, "zp": 25.0, "tie": [["1.x", "0.x"], ["1.y", "0.y"]],
     "components": [{"kind": "sersic", "x": 130, "y": 61, "mag": 16.8, "re": 3.0, "n": 2.0, "q": 0.7, "pa": 60, "fixed": ["n"], "bounds": {"q": [0.3, 1.0]}}]}

Python: `from ogfkit import multifit as MF; MF.fit(cutout, comps, psf=..., rms=..., mask=..., sky='const', tie=[...])`.

## Validation (synthetic images with known truth; `plugins/multifit/tests`)
| test | result |
|---|---|
| single Sersic, n = 1 / 2 / 4 (re 6 / 4 / 5, PSF FWHM 3 px, noise 1) | dmag -0.010 / -0.018 / +0.081, dre/re -0.4 % / -0.4 % / -13 %, dn/n +4.5 % / -0.3 % / -8.5 %, dq ≤ 0.02, dpa ≤ 4°, dpos ≤ 0.06 px, chi2_red 0.99 |
| n = 4 note | the fitted sky absorbs part of the wings: with the sky fixed to its true value the same fit gives dmag -0.001, dre 0.0 %, n 4.20; in a 161 px frame with fitted sky dmag +0.07 |
| error calibration, 40 noise realisations (mag 18, re 5, n 2) | pull mean / std: mag -0.09 / 0.79, re -0.12 / 0.99, n -0.13 / 0.95, q -0.25 / 1.05, x -0.18 / 1.06 |
| blend (target mag 17, re 5, n 2; neighbour mag 17.6, re 4, n 1) vs separation 1.0 / 1.6 / 2.8 / 4.8 re | mean abs dmag simultaneous **0.113 / 0.079 / 0.050 / 0.035**; masked neighbour 0.741 / 0.248 / 0.191 / 0.074; neighbour ignored 0.407 / 0.352 / 1.506 / 1.394; abs dre/re simultaneous 0.127 / 0.088 / 0.067 / 0.037 |
| bulge + disk (B/T true 0.344) | B/T 0.389, disk re 7.01 (7.0), q 0.65 (0.6), bulge re 2.84 (2.5), centres tied exactly, chi2 1.01 |
| PSF + host | PSF mag 19.10 (19.0), host mag 17.503 (17.5), re 5.64 (6.0), n 2.31 (2.0) |
| sky plane | gradient 0.0503 / -0.0295 (true 0.05 / -0.03), sky 99.98 (100); constant sky on the same image: chi2_red worse by > 0.2 |
| spatially varying PSF (FWHM 4.21 px at the field edge) | re error with the PSF model at the galaxy -4.5 %, with the central PSF +18.4 % |
| Poisson weights (gain 4) | dmag +0.028, dre/re -0.1 %, chi2_red 0.96 |
| catalog CLI, 7 galaxies (one blend pair 16 px apart, n 1-4) + 1 star | median dmag -0.019, max 0.069; max dre/re 0.107; chi2_red median 0.96, max 0.99; blend pair |dmag| sum 0.071 simultaneous vs 2.496 neighbour ignored; star (psf model) dmag -0.010, position error 0.003 px; residual rms in an n=3 galaxy region 4.30 → 0.97 (noise 1.0) |
| auto selection | star → psf, galaxy → sersic |

## Limitations
No analytic Moffat/Gaussian/Ferrers/Nuker/Fourier-mode/truncation components, no spiral/bar modes; the PSF of a component is evaluated at its start position (not re-evaluated while fitting);
the sky is local to the cutout (n ≳ 4 profiles and large galaxies in small cutouts have the usual sky-wing degeneracy: fix the sky for those); errors are formal (no covariance between components in
`GF_MAGERR`, which adds the component errors in quadrature); at most `max-neighbours` neighbours are fitted per cutout and only if they have > 1 % of the target flux; neighbours are started from the catalog
shape parameters (A_IMAGE, B_IMAGE, THETA_IMAGE, FLUX_RADIUS, FLUX_AUTO, CLASS_STAR); the model image sums target-only components over objects.
