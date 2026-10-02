# DAOPHOT-style crowded-field photometry (plugin `daophot`)

The classic DAOPHOT / ALLSTAR / DAOMASTER workflow as an OGFinder plugin: Measure tab, chip **DAOPHOT** (it sits in the chip row's **More** menu when the row is full).
The existing *Photometry* plugin (`psf_phot`, `crowded`) is **unchanged** (same argv, same columns); `daophot` is a sibling plugin that shares the PSF engine with
the PSF-modelling item (`ogfkit/psfmodel.py`).  The settings dialog is a notebook with the tabs **Find, Phot, PSF, Fit, Subtract/Add, Match**.

| DAOPHOT task | step (recorder name `analysis.daophot_<id>`) | implementation |
|---|---|---|
| FIND | `find` | `photutils.detection.DAOStarFinder` (Gaussian matched filter), threshold in sigma of the sky noise (height of the convolved peak, the DAOFIND convention), **sharpness and roundness cuts** |
| PHOT | `phot` | `dp.phot`: up to 3 circular apertures (exact pixel overlap), sky annulus, **sky modes** `median`, `mode` (3 median - 2 mean), `mean` (3-sigma clipped), `global` (background map); errors from sky noise + annulus uncertainty (+ Poisson with `gain`) |
| PICKPSF | `pickpsf` | brightest unsaturated stars with S/N, sharp/round cuts, no neighbour brighter than 10 % of their flux within 1.5 stamp half-widths, away from the edge; **manual add/remove** (`psf-add`, `psf-remove` = "x,y;x,y", 1-based) and the menu entries *Selected Rows -> PSF Stars* / *-> Not PSF Stars* that edit these two recorded parameters from the catalog selection |
| PSF | `psf` | `psfmodel.build_psf_model`: analytic **gaussian / moffat / lorentz / penny** base fitted to the recentred star stack, plus the **lookup table** of residuals (sparse weighted least squares on the pixel grid, cubic-convolution sub-pixel shifts, oversampling 2 with >= 25 stars), **spatial variation order 0-3** (reduced automatically: >= 150 stars order 3, >= 60 order 2, >= 20 order 1, else constant; 1-3 stars a Moffat fit; no star a Gaussian of the prior FWHM; `info['mode']` says which), **iterative neighbour subtraction** (`neighbour-iter` passes: fit all stars with the current model, subtract the neighbours of every PSF star, refit the PSF) |
| NSTAR / ALLSTAR | `fit` (and `run` = the whole chain) | `dp.allstar`: groups of overlapping stars (connected components of pairs closer than 2 fit radii, chunked above `max-group`), simultaneous Levenberg-Marquardt fit of flux, x, y per star and a common sky offset per group (**sky re-fit**), fit radius (`fit-radius`, default max(2.5, 1.5 FWHM)), **iterative re-find on the residual image** (`n-iter` passes; DAOFIND again, stars closer than 0.9 FWHM to a fitted star are skipped), parallel over groups (fork pool, `n-workers`) |
| SUBSTAR / ADDSTAR | `substar`, `addstar` | subtract the fitted stars (keep `keep-ids`) -> `daophot_sub.fits`; add `add-n` artificial PSF stars at random positions -> `daophot_add.fits` + truth table |
| residual frame | after `fit`/`run` | `daophot_resid.fits` opens in a new frame (parameter `show-frames`), the original frame stays current |
| DAOMATCH / DAOMASTER | `match` | second image (`second`) gets the same pipeline; **triangle matching** of the 20 brightest stars of both lists gives a similarity transform, then iterated mutual-nearest-neighbour matching with growing transform order (`shift`, `similarity`, `affine`, `poly2`, `poly3`, 3-sigma clipping); `dp.master_list` builds a master list over any number of frames with per-frame magnitude offsets |
| aperture correction | `apcorr` (also inside `fit`) | growth curve of the PSF stars with **all other fitted stars subtracted**: `corr(R) = median(m_ap(R) - m_psf)`, compared with the enclosed-flux curve of the PSF model; `DAO_MAG_APC = DAO_MAG + corr(apcorr-radius)` |
| CMD-ready output | `match` | `daophot_cmd.tsv` (matched pairs: positions, both magnitudes, `COLOR = mag1 - mag2`, errors), `daophot_master.tsv`, `daophot_cmd.png`; catalog columns `DAO_MATCH DAO_MX DAO_MY DAO_MSEP DAO_MMAG2 DAO_COLOR` |
| quality flags and cuts | `fit` | `DAO_FLAG` bit field + `DAO_GOOD` (see below), cuts `chi-max`, `snr-min`, `fit-sharp-lo/hi`, `fit-round-lo/hi` |
| artificial-star test | `arttest` | `plugins/completeness` (item B) with the DAOPHOT fitter (`dp.DaophotDetector`) as the detection callback and the **spatially varying PSF model as the injected PSF** (`run_completeness(psf=PSFModel)`); injected on top of the real stars (`avoid_detected=False`, minimum separation 2.5 FWHM between injected stars); catalog columns `DAO_ART_LIM50/90` |

Extra windows: *Star Table...* (all fitted stars), *Diagnostic Plots...* (chi, sharp, magnitude error vs magnitude; PSF FWHM and ellipticity maps; aperture-correction curve; `daophot_diag.png`), *CMD (matched)...*.

## Outputs

Everything is written to `~/.ds9/daophot/` (`--work`, recorded as `@{WORK}/daophot`; a replayed session uses its own work directory): `daophot_find.tsv`, `_phot.tsv`, `_psfstars.tsv`,
`_psf.json` + `_psf.fits` (the PSF model, JSON-serialisable; the FITS cube is a PSFEx-like layout), `_psf_diag.json`, `_stars.tsv`, `_resid.fits`, `_sub.fits`, `_add.fits`, `_add_truth.tsv`,
`_result.json`, `_diag.png`, `_cmd.tsv`, `_master.tsv`, `_cmd.png`, `_art.json`.

Catalog columns after *fit* / *run* (nearest fitted star within `match-radius`): `DAO_FLUX DAO_FLUXERR DAO_MAG DAO_MAGERR DAO_MAG_APC DAO_X DAO_Y DAO_CHI DAO_SHARP DAO_ROUND DAO_SNR DAO_SKY DAO_NGRP DAO_FLAG DAO_GOOD`.
Stars found by the re-find passes that have no catalog row are only in `daophot_stars.tsv` (the catalog is not extended).

Definitions: `CHI` = sqrt(reduced chi^2) of the group fit; `SHARP` = residual in the PSF core (r < 0.7 FWHM) projected on the PSF, relative to the star flux (0 = perfect, > 0 = more compact than the PSF
(cosmic ray), < 0 = broader (galaxy/blend)); `ROUND` is the DAOFIND `roundness1` of the star at its first detection (not recomputed after the fit); `SNR` = flux / error.
`DAO_FLAG` bits: 1 not converged, 2 chi > `chi-max`, 4 sharp outside range, 8 S/N < `snr-min`, 16 fit radius reaches the image edge, 32 group chunked, 64 flux <= 0 / no error,
128 found on a residual pass (not in the first list), 256 round outside range.  `DAO_GOOD` = no bit of 1|2|4|8|64|256.

## CLI (this is what the exported session script replays)

    python plugins/daophot/daophot.py IMAGE --mode run|find|phot|pickpsf|psf|fit|substar|addstar|apcorr|match|arttest --work DIR
        [--catalog CAT] [--mask M] [--mag-zeropoint 25] [--n-workers N] [--fwhm 3] [--find-thresh 4] [--psf-kind empirical|gaussian|moffat|lorentz|penny] [--psf-order 2] ...
    (`--help` lists all options; every GUI parameter has the same name)

Python API (pure functions, JSON-serialisable results; `ogfkit/daophot.py`, `ogfkit/psfmodel.py`): `find`, `phot`, `pick_psf`, `build_psf_model`, `allstar`, `run_pipeline`, `substar`, `addstar`,
`aperture_correction`, `match_lists`, `master_list`, `apply_transform`, `DaophotDetector`.

## Validation (measured, `plugins/daophot/tests`, 18 tests, ~65 s; numbers from `tests/validate_numbers.py`)

Synthetic crowded field 600 x 600, 1200 stars (power-law luminosity function, mag 15.5-23, zero point 25, sky 50, noise 0.5), **spatially varying Moffat PSF** (FWHM 2.8 -> 3.8 px
along x, ellipticity 0.04 -> 0.14 along y, beta 3), truth known.  Fit radius 4.9 px, 3 passes.

| PSF model | residual rms / noise | flux scatter dF/F (mag 15.5-17) | (18-19) | (19-20) |
|---|---|---|---|---|
| order 0 (constant) | 1.184 | 0.094 | 0.186 | 0.158 |
| order 1 | 1.057 | 0.011 | 0.028 | 0.036 |
| order 2 | 1.082 | 0.013 | 0.041 | 0.035 |

(the `test_spatially_varying_psf_beats_constant` test requires the order-2 scatter < 0.8 x the constant-PSF scatter.)  Order-2 run, per truth magnitude bin:

| mag | stars | recovered (< 1 px) | completeness | median dF/F (PSF mag) | dpos median (px) | chi median |
|---|---|---|---|---|---|---|
| 15.5-17 | 10 | 10 | 1.00 | -0.019 | 0.029 | 0.49 |
| 17-18 | 17 | 17 | 1.00 | -0.027 | 0.069 | 0.76 |
| 18-19 | 41 | 41 | 1.00 | -0.028 | 0.072 | 0.96 |
| 19-20 | 77 | 74 | 0.96 | -0.029 | 0.086 | 0.99 |
| 20-21 | 151 | 99 | 0.66 | -0.022 | 0.205 | 0.97 |
| 21-22 | 296 | 17 | 0.06 | (selection-biased) | 0.39 | 1.00 |

* **Flux bias.** The PSF flux is that inside the PSF stamp (the wings of a beta = 3 Moffat beyond 10 px are 1-3 %), so PSF magnitudes are 2-3 % too faint; the aperture correction measures this:
  after `DAO_MAG_APC` the median bias of stars brighter than mag 19 is within 0.03 mag (test), 0.014 / 0.009 / 0.000 mag in the 15.5-18 / 18-20 / 20-21 bins in the default pipeline run, scatter (MAD) 0.03-0.08 mag.
  The error estimates are honest: the pull (F_fit - F_true)/sigma of mag 19.5-20.5 stars has median |pull| < 0.5 and robust sigma 0.6-1.6 (test).
* **Completeness.** Limited by the FIND threshold (4 sigma of the *peak height*: 50 % at mag ~20.5, analytic peak S/N 5.4 for flux 60 at FWHM 3.5); with a threshold of 3 sigma the 50 % point moves ~0.4 mag fainter at the price of more false stars.
  Artificial-star test (DAOPHOT fitter as the callback, 80 injected stars into the crowded field): recovered fraction 1.00 / 0.75 / 0.19 / 0 / 0 at mag 19.4 / 20.2 / 21.0 / 21.8 / 22.6, 50 % limit 20.55, 90 % limit 19.92; bias at mag 19.4 +0.003 mag, scatter 0.027.
* **Positions.** Median error 0.03-0.09 px for mag < 20, 0.2 px at mag 20-21 (S/N 5-10).  **False stars:** 5-10 of 250 fitted stars have no truth within 1.5 px.
* **PSF model.** FWHM of the model within 5 % and ellipticity within 0.02 (median) / 0.07 (field corners, polynomial extrapolation) of the truth on a 3 x 3 grid; oversampling 2 is accurate to ~1 % in FWHM (oversampling 1 with cubic interpolation sharpens the model by ~3 %, so it is only used with < 25 stars).
* **DAOMATCH.** 120 stars, 25 % missing in frame 2, 15 spurious stars, rotation 0.8 deg, scale 1.01, shift (13.4, -7.2), position noise 0.05 px, magnitude offset 0.5: at least all but 3 of the true pairs found (asserted), rms < 0.12 px, offset recovered to 0.02 mag; a 90-degree rotation + shift is found without a guess (similarity, rms < 0.05 px).
  End-to-end (two synthetic frames, both through the pipeline): 173 pairs, rms 0.149 px, offset 0.492 mag (scatter 0.067).
* **Analytic functions.** A constant Moffat PSF (FWHM 3.2, beta 3) is recovered by `moffat` to 3 % in FWHM and 0.5 in beta (test); `penny` within 6 %; `gaussian` / `lorentz` only give the core roughly (25 % / 50 % in FWHM) - which is what the lookup table fixes.
* **Real data.** M51 (`/workspace/fits/m51.fits`, undersampled, galaxy light, saturated cores; FWHM prior 3.5): 456 candidates, 40 PSF stars (37 used, order 1, model FWHM 1.89 px), 529 fitted stars (373 pass the cuts), aperture correction +0.119 mag at 10 px.  There is no truth: the many "new" stars on the residual passes include HII regions and clumps of the galaxy.

## Limitations

* Pixel-based coordinates only (DAOMATCH fits polynomial transforms between the star lists; the WCS is not used).  The catalog is not extended with the stars that have no catalog row.
* The fit weights are 1/sigma_sky (+ Poisson of the model with `gain`); the background is the sep mesh map plus a per-group offset, which is biased high in very dense fields (sky annulus / `global` modes differ by the same amount).
* PSF stars with neighbours rely on the neighbour subtraction; with fewer than ~20 PSF stars the model is a constant lookup table.  Saturated stars must be excluded with `saturation`.
* `ROUND` is not recomputed after the fit; `SHARP` is residual-based (not DAOPHOT's definition); chi uses the whole fit region.
* DAOMASTER here is the matching + magnitude-offset part; it does not iterate proper-motion/variability solutions.  The CLI `match` mode runs the full pipeline on the second image (no PSF reuse across images).
