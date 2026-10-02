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
| catalog CLI, 7 galaxies (one blend pair 16 px apart, n 1-4) + 1 star | median dmag -0.019, max 0.069; max dre/re 0.107; chi2_red median 0.96, max 0.99; blend pair |dmag| sum 0.071 simultaneous vs 2.408 neighbour ignored; star (psf model) dmag -0.010, position error 0.003 px; residual rms in an n=3 galaxy region 4.30 → 0.97 (noise 1.0) |
| auto selection | star → psf, galaxy → sersic |

## GALFIT interoperability (`ogfkit/galfitio.py`)
* **Import**: `--config FILE` (step "Fit From Config File") accepts a **GALFIT feedme** as well as the JSON config (detected by content). With image `-` the data, PSF (D), mask (F), sigma image (C) and fitting region (H) named in
  the feedme are used (paths relative to the feedme). Objects: `sersic`, `devauc`, `expdisk` (R_e = 1.67835 R_s), `psf`, `gaussian` (= Sersic n = 0.5, R_e = FWHM/2, exact) and `sky` (value, dsky/dx, dsky/dy with the free/fixed
  flags -> constant / plane / fixed sky). Conventions are converted: GALFIT PA (up = 0, left = 90) = multifit PA - 90; positions are 1-based; magnitudes include the exposure time (`EXPTIME` of the image: m_multifit = m_GALFIT - 2.5 log10 t_exp);
  the sky value is moved from the centre of the GALFIT region to multifit's array centre. Fix flags become `fixed`; a **constraints file** (G) gives bounds (`N param lo hi` relative, `N param lo to hi` absolute, for
  x y mag re rs n q pa) and ties (`N_M param offset`, accepted when the start values are equal, i.e. common centres / PAs). Anything else (moffat, nuker, ferrer, king, edgedisk, Fourier modes, ratio / unequal-offset constraints, PSF fine
  sampling E > 1, diffusion kernel) is reported (`ValueError` for objects, a warning for constraint lines) and not silently ignored. Config / feedme fits start from the given values and from `--restarts` (default 2) alternative
  starts (R_e x 0.6 / x 1.6, n x 1.3 / x 0.75); the lowest chi2 wins (`starts_chi2` in the result JSON).
* **Export**: `--export-feedme DIR` (parameter "Also write GALFIT feedme files") writes the *fitted* model as a feedme: config mode -> `galfit.feedme` (+ `galfit.constraints` when bounds / ties exist), catalog mode ->
  `galfit_<NUMBER>.feedme` with the target and the fitted neighbours (masked neighbours are not exported), the PSF stamp at the object (`galfit_<NUMBER>_psf.fits`), region = cutout, sigma none (GALFIT builds its own), mask = your mask.
  The file starts GALFIT from the multifit solution (set `P) 1` to render it, `P) 0` to refine it). `galfitio.parse_feedme` / `write_feedme` / `write_constraints` can be used from Python.
* x / y **bounds** in a config are now given in image coordinates like `x` / `y` (before they were silently interpreted in cutout coordinates).

## Head-to-head with GALFIT (`plugins/multifit/validation/`, GALFIT 3.0.5 Linux binary downloaded from the author's page for these tests only; it is not in the repository and may not be redistributed)
Run with `GALFIT_BIN=/path/galfit [GALFIT_LD=dir with libncurses.so.5/libtinfo.so.5]`; without a binary the GALFIT-specific tests skip and everything else runs.
1. **Rendering** (`tests/test_galfit.py::test_rendering_matches_galfit_pixelwise`; GALFIT P=1 vs `render_model` for the same feedme, 81x93, Gaussian PSF sigma 1.4 px, sky with gradient): max |difference| / peak =
   0.9 % (Sersic n 2.7), 0.15 % (expdisk), 1.4 % (devauc, the cusp pixel), 0.17 % (gaussian), 0.37 % (psf), 0.47 % (all five + sky); total flux differs by 0.0006-0.36 %. Position, PA, R_s, exposure-time and sky-gradient conventions are therefore right.
   **A bug was found and fixed this way**: `render_sersic` used an odd sub-pixel grid in the cusp that sampled r = 0 exactly; for n = 5 centred on a pixel the peak was +16 % and the flux +4 % (regression test added; golden
   behaviour of other plugins unchanged: isophote 8, morph_ext 13, completeness 9, multifit tests pass). Against a brute-force 120x120 sub-sampled reference (pixel-integrated Sersic convolved with the PSF image) the multifit model is now accurate to
   0.15 % of the peak (n = 5, R_e 2.9, real PSF, off-centre; 0.04 % for n = 4, R_e 6) while GALFIT's own rendering differs from that reference by 7.5 % (4.4 %) of the peak off-centre and 2.5 % on a pixel centre - GALFIT's
   cusp sampling is the coarser one for steep profiles (the reference assumes the PSF image already includes the pixel response, as both programs do).
2. **Fits on synthetic images** (`galfit_compare.py`, report `galfit_compare_report.json`; 30 random objects per case, 101x101, Moffat PSF FWHM 3.2 px, constant sigma 2 (sigma image given to both), sky 50, **images rendered by GALFIT itself**, both programs
   fit the same start feedme (start values perturbed by up to +-1 px, +-0.4 mag, R_e +-30 %, n x 0.6-1.5, +-30 deg) with the same constraints file; chi^2 of every result is evaluated with GALFIT's renderer on the same data):

   | case | chi2 difference multifit - GALFIT (median; fraction within 0.5 / multifit lower / GALFIT lower) | multifit - GALFIT, 90th percentile of the absolute difference | median time GALFIT / multifit CLI |
   |---|---|---|---|
   | A single Sersic (n 0.8-4, R_e 3.5-9, q 0.4-0.95) | +0.006; 0.97 / 0 / 0.03 (max 0.66) | x 0.025 px, mag 0.021, R_e 2.8 %, n 6.8 %, q 0.004, PA 0.23 deg | 0.24 s / 1.6 s |
   | B devauc + expdisk, tied centre | +0.07; 0.50 / 0.17 / 0.33 (-14.3 ... +2.7) | bulge: mag 0.41, R_e up to 23x (degenerate); disc: mag 0.08, R_e 7 % | 5.0 s / 6.9 s |
   | C psf + Sersic host, tied centre | -0.08; 0.63 / 0.30 / 0.07 (-1.6 ... +0.6) | host mag 0.17, R_e 9.8 %, n 48 %; psf mag 0.34 | 0.38 s / 1.6 s |
   | D Sersic + Sersic neighbour (15-21 px away) | +0.11; 0.80 / 0 / 0.20 (max 1.8) | mag 0.024, R_e 3.8 %, n 5.5 %; neighbour mag 0.035 | 0.50 s / 3.4 s |
   | E Sersic + free sky gradient | +0.007; 1.00 / 0 / 0 (max 0.29) | x 0.018 px, mag 0.017, R_e 2.5 %, n 3.6 % | 0.25 s / 1.6 s |

   chi2 per degree of freedom: 0.994-1.006 for both programs in all cases. Equal chi2 means both reach the same minimum; the parameter differences in B and C are flat directions of chi2 (bulge/disc and PSF/host degeneracies), not different fit quality. Against the *truth*
   both programs have the same, and large, errors in these deliberately faint/extended/free-n/free-sky cases (case A 90th percentile |error|: mag 0.44 GALFIT / 0.42 multifit, R_e 100 % / 97 %, n 70 % / 69 %), i.e. the comparison is about equivalence, not about
   absolute accuracy. Timing: GALFIT is faster per object (multifit's number includes ~0.6 s Python/scipy start-up and a 3-point numerical Jacobian).
3. **pysersic 0.1.5** (JAX, maximum a posteriori, flat sky; `pysersic_map.py`, `compare_pysersic.py`; case A only, priors centred on the same start values with widths 2 px / 50 % flux / 50 % R_e, uniform n 0.5-6): chi2 difference to GALFIT median +0.12 (90th
   percentile +0.66, max +2.1; pysersic worse by more than 1 in 7 % of the objects); its priors suppress the large-R_e outliers (90th percentile |dR_e/R_e| 24 % vs 100 %), at the price of being prior-dependent; 3.3 s per object on CPU. imfit and galfitm have no pip package and were not tried.
4. **GALFIT's own example** (`galfit_example.py`, the gal.fits/psf.fits/galfit.feedme of the GALFIT distribution; a real image, with a PSF image that has a 0.1 px off-centre centroid and negative noise pixels): GALFIT (own Poisson sigma) m = 22.976, R_e = 2.623, n = 4.110, q = 0.857, PA = -71.75;
   multifit (constant sigma, same feedme): m = 22.931, R_e = 2.865, n = 4.993, q = 0.819, PA = -68.86 (positions agree to 0.015 / 0.056 px). With one constant sigma image given to both, chi2 (rendered by GALFIT) is 13663 for GALFIT and 13868 for multifit
   (+204 on 8649 pixels, 1.5 %), whereas rendered by multifit's model the multifit solution has chi2 13630 and GALFIT's parameters 13742: on this steep profile (n = 5) with an undersampled, off-centre PSF the two programs' renderings differ at the 5-8 % level in the two brightest pixels
   (see 1.), and each fit prefers its own model; which one is closer to the sky cannot be decided from this image. This example also exposed a weakness of the optimiser (2-point finite-difference Jacobian: 150 evaluations, not converged, wrong minimum n = 7.5, R_e = 23); the Jacobian is now 3-point (20 evaluations, converged) and config fits use multi-start.

Not done: no comparison on real survey images with an independent truth; multi-component fits to real galaxies (bulge + disc + bar) not compared; GALFIT features not implemented above; the GALFIT result is a closed-source binary run on Linux only; the synthetic cases cover 150 fits each for two programs, not a statistically exhaustive survey of the parameter space.


## Advanced components (GALFIT 3 compatible; `ogfkit/profiles.py`)
Used from `--config` (JSON or GALFIT feedme), `multifit.py IMAGE --config cfg --work DIR`; the dialog step "Fit From Config File" has the same parameters (`config`, `mask-catalog`, `mask-exclude`, `neighbour-radius`).
All keys are flat entries of a component dict, so they can be fixed, bounded and tied like any other parameter.

| feature | keys (GALFIT feedme keyword) |
|---|---|
| generalised ellipse (boxy/disky) | `c0` (`C0`) |
| Fourier modes (azimuthal) | `f<m>a`, `f<m>p` amplitude and phase [deg] (`F<m>) amp phase`), any m |
| bending modes | `b<m>` (`B<m>`), any m (B1 = banana, B2 = S-shape) |
| coordinate rotation (spiral) | `rot_func` `'power'`/`'log'`, `rot_in`, `rot_out`, `rot_theta` [deg], `rot_alpha` (power) or `rot_ws` (log), `rot_incl`, `rot_pa` (`R0`-`R4`, `R9`, `R10`) |
| truncation | separate component `kind='trunc'` with `rbreak`, `dsoft` (optional own `x,y,q,pa`, else the parent's), referenced by `trunc_in`/`trunc_out` = list of component indices of the parent (`T0) radial`, `T4`, `T5`, `Ti)`, `To)`) |
| flux normalisation | `norm`: `total` (default, `flux`), `center`, `re`, `break` (GALFIT `sersic1/2/3`, also `expdisk1/2/3`, `devauc1/2/3`; then the amplitude is the surface brightness `i0` in counts/pixel) |
| other profiles | `moffat` (flux, fwhm, beta), `ferrer` (i0, rout, alpha, beta), `king` (i0, rc, rt, alpha), `nuker` (ib, rb, alpha, beta, gamma), `edgedisk` (i0, hs, rs, pa), `brokenexp` (i0, h1, h2, rbreak, alpha; multifit-only, no GALFIT counterpart; `alpha` fixed unless freed) |

Ties (`tie` list): `["1.x","0.x"]` equal; `["1.re","0.re","offset"]`, `["1.flux","0.flux","ratio"]`, optional 4th element = value. Without a value the offset/ratio is that of the start values, as in GALFIT
(a GALFIT `mag offset` is a flux ratio). Neighbour masking in config fits: `--mask-catalog TSV --neighbour-radius K --mask-exclude PX` masks the ellipses (K x A_IMAGE, B_IMAGE, THETA_IMAGE) of the catalog objects in
the fitted box except those within PX pixels of a fitted component (the catalog mode already had `--neighbours mask`). Results: `multifit_results.tsv` has the usual columns plus `EXTRA` (`name=value;...` of all other parameters); `--export-feedme` writes all of it back as feedme + constraints.

Feedme/constraints import-export (`ogfkit/galfitio.py`): C0, F1..Fn, B1..Bn, R0-R4/R9/R10, T0-T5/T9/T10 truncation objects, `Ti)/To)` object references, `sersic1/2/3`, moffat, ferrer, king, nuker, edgedisk; surface-brightness
objects are converted with the plate scale and exposure time (`i0 = dx dy t_exp 10^(-0.4 (mu - ZP))`). The constraints file follows what GALFIT 3.0.5 actually accepts (measured): parameters are named by position -
`x y mag q pa c0`, `f<m>` (Fourier amplitude), `f<m>p` (phase), `re/rs/rb` = parameter 4, `n/alpha` = 5, `beta` = 6, `gamma` = 7 (for surface-brightness objects `mag` constrains the surface brightness; any unknown name crashes GALFIT);
bending and rotation parameters cannot be constrained in GALFIT (written as comments / warned on import); `radial2`, length/height truncation and other rotation functions are rejected on import.

Conventions measured against the GALFIT 3.0.5 binary (not from the paper where they differ): Fourier phase uses the mirrored angle `theta = atan2(y'/q, -x')` with `r = r0 (1 + sum a_m cos(m (theta + phi_m)))`; bending `y' = y + sum (-1)^(m+1) b_m (x/R_scale)^m`
(R_scale = R_e; R_s for exponentials; FWHM-radius for Moffat, R_out for Ferrer, ...); rotation uses the tanh edge function with `A = 2 (20 deg)/|theta_out| - 1.00001` (paper: 0.23), power law `[0.5 (r/r_out + 1)]^alpha`, log `log(r/ws + 1)/log(r_out/ws + 1)`,
inclination deprojection about `rot_pa` before the spiral rotation; truncation `P = 0.5 (tanh((2 - B) r/r_b + B) + 1)`, `B = 2.65 - 4.95 r_b/(r_b - r_s)` (binary 4.95, paper 4.98); small-radius sharp truncation edges are 5x5 sub-sampled.

Validation (`plugins/multifit/validation/galfit_advanced.py parity|h2h|recovery`; reports `galfit_advanced_*.json`; GALFIT 3.0.5 only for the local comparison, tests skip without `GALFIT_BIN`/`GALFIT_LD`):
* **Render parity** with GALFIT (P=1, 101x101, 25 model definitions each without and with a Moffat PSF, models written to GALFIT through `write_feedme`, so the exporter is tested too): all 50 cases agree, median max|d|/peak 0.32 %,
  48/50 within 1 % of the peak and 49/50 within 2 %, worst 2.6 % (n = 4, R_e = 4 px, PSF-convolved cusp pixel), total flux within 0.21 % (worst: Moffat, 0.21 %). By component (max|d|/peak, no PSF): C0 0.07-0.19 %, Fourier 0.06-0.22 %, bending 0.06 %,
  rotation power/log 0.08-0.29 %, outer truncation 0.06 %, inner truncation 0.8 % (sharp edge, sub-sampling), two truncations 1.6 %, Moffat 0.2 %, Ferrer 0.13 %, King 0.12 %, Nuker 0.54 %, edge-on disc 0.79 %.
* **Head-to-head fits** on GALFIT-rendered data (101x101, Moffat PSF, constant sigma, same start feedme/constraints for both, 12 realisations per case, chi2 re-evaluated with GALFIT's renderer): chi2(multifit) - chi2(GALFIT), median / range / fraction within 0.5:
  Fourier F1+F3 +0.003 / [-0.25, +40.9] / 0.92; bending B1+B2 +0.010 / [-1.6, 5.2] / 0.73 (GALFIT crashed once); rotation (power, theta free) +0.077 / [-6.7, 0.5] / 0.83; outer truncation (R_break free) +0.014 / [-1.0, 0.09] / 0.92;
  boxy Sersic + Moffat +0.007 / [-2.6, 0.02] / 0.92; tied bulge + disc (x/y offset, flux ratio) +0.27 / [0.006, 0.87] / 0.67 (GALFIT lower in 4/12 by < 1); masked bright neighbour with F2 -2.6 / [-5.5, -0.09] / 0.17 (multifit lower in 10/12:
  GALFIT stops in a q/F2 degeneracy). Parameter differences multifit - GALFIT (median / 90th percentile |.|): Fourier mag 0.0001/0.002, R_e 0.07 %/0.25 %, F3 amplitude 0.0004/0.0009; truncation R_break 0.0004/0.027 px (relative 0.04 %/2.7 %); Moffat FWHM 0.2 %/0.3 %;
  well-constrained parameters agree to better than 1 %, the degenerate ones (bending B1/B2, rotation theta with PA, F2 with q) differ by the same amounts as either differs from the truth. Wall time per fit: multifit CLI 1.8-4.7 s, GALFIT 0.5-5 s.
* **Synthetic truth recovery** (ogfkit renderer, Moffat PSF, sigma 2, sky 50, 12 realisations per case, `fit_multistart` with 2 restarts, 600 evaluations): reduced chi2 0.99-1.01 in all 14 cases, converged 12/12; parameter errors (median / 90th percentile |error|):
  Fourier (m = 1, 3) mag 0.0001 / 0.12, R_e -1.4 % / 12 %, F3 amplitude -0.14 / 0.19 (strongly degenerate with q at S/N ~ 100; the sign of an amplitude is degenerate with a phase shift of 180/m); truncation R_break 0.0003 / 0.36 px, R_e 1.7 % / 4.4 %;
  King r_c +4 % / 13 %, r_t -6 % / 21 %; Ferrer r_out 0.5 % / 8.5 %; Nuker r_b 1.6 % / 32 %, gamma -0.017 / 0.11; edge-on disc h_s 1.4 % / 4.7 %, r_s 0.5 % / 4 %, PA 0.2 / 0.4 deg, centre y scatter 0.3 px;
  broken exponential h_1 5 % (90th percentile 190 %: the break and h_1 are degenerate when the break is near the edge), R_break -4 % / 50 %; inner truncation R_break 0.0 / 2.6 %; spiral (log) rotation theta_out median -16 deg, 90th percentile 109 deg (the winding angle is only weakly constrained);
  bending B1 90th percentile 3.5 (unbounded fit, degenerate with PA; fits needing B1/B2 should be bounded); boxy C0 absolute error 0.6 / 1.8 in a joint boxy Sersic + Moffat fit (C0 is weakly constrained).
Reproduce: `GALFIT_BIN=... GALFIT_LD=... python plugins/multifit/validation/galfit_advanced.py parity|h2h|recovery OUT.json [--n 12]` (recovery needs no GALFIT).

Known differences and limitations: GALFIT 3.0.5 normalises the flux of a profile that has **both** C0 and Fourier modes differently from Peng et al. 2010 (factor >= 2, depends on amplitudes and phases; `C0 = 0` with Fourier modes crashes it) -
multifit uses the paper's normalisation, so the total magnitudes of such objects differ from GALFIT's (a warning is printed on feedme import; shapes agree to 0.2-0.6 % after scaling); GALFIT returns zeros for a truncated plain `sersic` (needs `sersic1/2/3`) and multifit rejects that feedme;
no PSF fine sampling (`E) != 1` is refused); only `radial` truncation; rotation only power/log; the broken exponential, `rot_*`/`b*` bounds and amplitudes beyond GALFIT's constrainable set are multifit-only extensions; the catalog mode (`--catalog`) does not use the advanced components;
the advanced components are numerically sub-sampled and therefore slower than the plain Sersic path (not benchmarked separately).

## Limitations
The GALFIT advanced components (next section) are available only in config / feedme fits, not in the catalog mode (`--model` choices unchanged); the PSF of a component is evaluated at its start position (not re-evaluated while fitting);
the sky is local to the cutout (n ≳ 4 profiles and large galaxies in small cutouts have the usual sky-wing degeneracy: fix the sky for those); errors are formal (no covariance between components in
`GF_MAGERR`, which adds the component errors in quadrature); at most `max-neighbours` neighbours are fitted per cutout and only if they have > 1 % of the target flux; neighbours are started from the catalog
shape parameters (A_IMAGE, B_IMAGE, THETA_IMAGE, FLUX_RADIUS, FLUX_AUTO, CLASS_STAR); the model image sums target-only components over objects.
