# Precision batch (items 1-7): final report (2026-10-03, all commits local, not pushed)

**Final full run** (`scripts/run_all_checks.sh`, log `/workspace/final_checks4.log`, out `/workspace/work/final4`, HEAD 5fd061d1c):
**45 checks: 45 PASS, 0 FAIL, 0 SKIP** (incl. `newplugins` GUI 1402 s, `lensmodel_gui`, `regression_tests`, `regression_data` 11/11 cases in 271 s,
`moving_session`, `session_replay`).  `ai_bridge_tests` ran with the live-agent tests skipped (`OGF_LIVE=1` includes them; they are flaky on transient service errors,
fixed in 43d9fc8ca; the previous run's one failure was exactly that live test).  `agent_real` passes although the grok CLI there reports "not signed in".

## 1  GALFIT-advanced / real-galaxy cross-check, bar-ring-spiral (5c7aec59e, 364a86acf; earlier 0cbc7148c 4dbbdba66 937fa314d 087801bb9)
* 93 real-galaxy cutouts vs GALFIT: single-Sersic mag within 0.05 in 91 %, NMAD 0.0024; R_e NMAD 0.23 %; chi2 within 0.5 % in 89 %.
* Bulge+disc is degenerate: GALFIT failed on 9/93, multifit on 0.  `structfit.py`: bar/ring/spiral found in 77/80 (high S/N) and 72/80 (low S/N) synthetic, 0/20 false features.
* Real M95: bar + ring selected, bar PA 108 deg (literature ~112), Ferrers bar length overestimated.  Limitation: spiral/ring detection validated mostly on synthetic images; one real galaxy.

## 2  PSF modelling (512e45fce, d96cecfd0; `ogfkit/psfext.py`, psfex `--rank --wings --ext-size --assess`)
* Synthetic max|model-truth|/peak: degree 0 0.150 -> 1 0.044 -> 2/3 0.041, rank-1 0.040; halo gamma 3.07 recovered (wing flux +13 %).
* Real HUDF F160W held-out stars: residual 0.350 -> 0.243 (degree 1); rank reduction gave no gain; wings not measurable in this field.
* Limitations: halo recovery verified only on synthetic data; rank reduction is post hoc; no new star re-selection.  In the 1300^2 regression crop degree 1+ is *not* better than degree 0 (0.279 vs 0.258) - recorded as is.

## 3  Photometric error model (fc36a0781, dc80b2cde; `ogfkit/photerr.py`, noisemodel `--sky-annulus --psf`, daophot `--noise-meta`)
* Blank-aperture pull std (ideal 1): naive 1.7-3.7, legacy 1.7-3.2, measured local law 0.94-1.07 synthetic and 0.97-1.00 real HUDF (global-sky law 0.46-0.72).
* Injection: isolated 0.95-1.07; close pairs 4.3-4.6 -> 1.8-2.0 with the contamination model (still not 1); real-data injection aperture correction off by 2-8 %.
* Not replaced: stacking bootstrap errors, multifit covariance, moving S/N.  Real-data injection test has only 7 evaluable sources in the regression crop.

## 4  Moving objects: orbit-population prior (5ca4ca372; `moving/orbitlink.py`, `ds9_moving.py --orbit-prior`, SDSS-vs-SkyBoT validation)
* SDSS run 94 (Stripe 82) vs SkyBoT: camcol 3 baseline 342 tracklets / 49 true (precision lower bound 0.14, recall 0.82, ceiling 60) -> 99 % prior 83 / 49 (0.59, 0.82).
  Camcol 4 (nothing refit): 350 / 67 (0.19, 0.97) -> 110 / 67 (0.61, 0.97).  68 moving tests pass.
* Limitations: this is a rate-prior vetting from an orbit population, **not** cross-night orbit-fit linking (not implemented); a 1.5-5.5 AU osculating two-body population, H<=19.5, no NEO/comets/TNOs;
  two fields, one night; precision is a lower bound (SkyBoT-incomplete truth); 3-exposure tracklets are still dominated by false links.

## 5  Photo-z closure (fe8ff423d; `ogfkit/pz_closure.py`, `photoz_quality.py --recal-out/--recal-in/--mode consistency`)
* SDSS MDN: map fitted on 1500, tested on an independent 1500: PIT KS p 5e-4 -> 0.18, CRPS unchanged 0.0261 (MDN already near calibrated); 5-fold cross-fit KS p 4e-7 -> 1.0.
  Injected width errors x0.5 / x2 recovered (68 % coverage 0.463 / 0.932 -> 0.677 / 0.680).
* EAZY (eazy-py 0.8.7, no prior, ugriz, z<=1) vs MDN: sigma_NMAD 0.041 vs 0.0235, outliers 7.7 % vs 1.3 %; 2.9 % flagged at 3c, 95 % of flagged are EAZY outliers; inverse-variance merge (0.0262) is worse than the MDN alone.
* Limitations: one survey (SDSS low-z), one global map, no GUI step; EAZY comparison needs the external eazy environment (disabled in the regression case).

## 6  Strong lenses (42737073c, f3c403134, 3e8c1409c, 43d9fc8ca for tests; `ogfkit/lensextra.py`, `--task multiplane`, `--task source --source-method inversion`)
* Pixel inversion, 18 **synthetic** lenses: chi2_red 1.02, source correlation 0.9993, centroid 2.3 mas, flux ratio 1.001.  Pixel-level lens fit from an offset start: 14/18 within 1 % in theta_E
  (robust sigma 0.2 %, q 0.024, phi 2 deg), **4/18 doubles end 3-13 % off** (local optimiser) - use several starts.
* Multi-source-plane, 24 **synthetic** Jackpot-like lenses: joint fit theta_E scatter 0.16 % vs 0.52 % (plane 1 alone), q 0.009 vs 0.054, phi 1.8 vs 17 deg; free plane-2 weight 1.0006 +- 0.0028.
* **Real lens SDSS J0946+1006** (HST ACS F814W via MAST, reachable): theta_E 1.389", q 0.94, gamma 0.066; literature 1.397 (Collett & Auger) / 1.43 +- 0.01 (Gavazzi) -> agrees to 0.6-3 %.
  **chi2_red = 9.1**: fit is poor (Gaussian PSF stand-in, correlated drizzle noise treated as white, lens-light residuals); the outer ring (z=2.035) is not modelled/tested.
* Limitations: linear Tikhonov inversion (no positivity), SIE+shear only, multi-plane on image positions only, one real system/band, no GUI step.

## 7  Regression set on public data (5fd061d1c; `regression/`, docs/testing.md, checks `regression_tests` + `regression_data`)
* 11 cases: extract, noisemodel (blank apertures, injection), depth, psfex, stacking, isophote, multifit, moving (SDSS run 94), photoz_sed, lensmodel (J0946).  Data fetched on demand into `$OGF_DATA_CACHE`
  (HUDF12 F160W crop, SDSS run 94 frames [must be cached once by `sdss_known_asteroids.py`], MAST J0946); unreachable data -> SKIP.  271 s with data cached; 11 pass, 0 fail.
* Baselines are **self-baselines** (current code on these data), not truth; stack null z-mean +0.73 and the 7-source injection set are recorded as is with wide ranges.
* Not covered: ai_bridge, cluster, spectra, sedcodes, lightcurves, xmatch, daophot, morph_ext, batch/repro, GALFIT wrapper.

## Verified only on synthetic data
Item 2 halo/wings; item 3 close-pair contamination model; item 6 pixel/multi-plane accuracy numbers; item 1 spiral/ring detection rates.
## Not done
Cross-night orbit-fit linking (4); iterative PSF star re-selection (2); positivity/adaptive-mesh/sub-halo lens modelling and outer-ring fit on J0946 (6); GUI steps for the new lens and calibration tasks (5, 6);
regression coverage of the plugins listed under item 7; `ast/config.h.in` left uncommitted (not mine).
