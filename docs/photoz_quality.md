# Photo-z distribution tools (plugin `photoz_sed`: steps *Photo-z Calibration*, *Spec-z Sample Representativeness*)

Library `ogfkit/photoz_stats.py` (numpy / scipy), headless CLI `plugins/photoz_sed/photoz_quality.py`, tests `plugins/photoz_sed/tests`
(`photoz_tests` in `scripts/run_all_checks.sh`), validation `plugins/photoz_sed/validation/photoz_validate.py` (+ `photoz_*_report.json`).
They work on any photo-z catalog with a spectroscopic redshift for part of the objects: the built-in photo-z (`ds9_photo_z.py`, empirical or MDN), the EAZY/SED-code
steps (`EZ_Z`, `EZ_Z16`, `EZ_Z84`), the AI bridge (`PHOTOZ`, `PHOTOZ_P16/P84`) or an external catalog.

## Predictive distributions

| kind | needs | description |
|---|---|---|
| `gauss` | `--zphot-col`, `--zerr-col` | N(zphot, err) |
| `split` | `--zphot-col` (median), `--q16-col`, `--q84-col` | asymmetric Gaussian with sigma_lo = p50 - p16, sigma_hi = p84 - p50 |
| `mixture` | `--mixture-file` (npz: number, pi, mu, sigma) | Gaussian mixture; written by `ds9_photo_z.py --mixture-out` (MDN) |
| `grid` | `--pdf-file` (npz: number, zgrid, pdf) | tabulated p(z) |

`--pdf-kind auto` takes the mixture file, else the percentile columns, else the Gaussian.  All are truncated at `--zmin` (default 0) and renormalised.

## Step *Photo-z Calibration (PIT, accuracy)* (`--mode quality`)
Uses the objects with a spec-z (`--zspec-col`, or a TSV `--zspec-file` with NUMBER + the column).  Output in `{work}/photoz_quality/`:
* **Percentile accuracy** with dz = (zphot - zspec)/(1 + zspec): bias (mean, median), sigma_NMAD, sigma_std, outlier fraction (|dz| > `--outlier`, default 0.15, and > 3 sigma_NMAD),
  the 68/90/95/99th percentiles of |dz|, bootstrap errors; the same in bins of z_spec, z_phot and magnitude (`pzq_binned.tsv`), weighted when `--weights-col` is given.
* **Calibration of the PDFs**: PIT = F_i(zspec_i) (uniform if calibrated): histogram (`pzq_pit.tsv`), KS, Cramer-von Mises and chi^2 tests, mean / variance and a verdict
  (U-shaped = overconfident or outliers, hump = underconfident, tilted = bias), coverage of the 50...99 % equal-tail credible intervals with Wilson errors (`pzq_coverage.tsv`),
  CRPS and the logarithmic score, the z-score (zspec - mean)/std, the optimal width scale (minimum log score; gauss / split / mixture) and the PIT-recalibration check (2-fold cross-fit KS p),
  stacked p(z) against the histogram of z_spec (`pzq_stacked_nz.tsv`: KS distance of the CDFs, means), and `pzq_plot.png` (PIT, coverage, Q-Q, z_phot vs z_spec).  Everything is in `pzq_report.json`.
* Catalog columns (spec-z objects only): `PZ_PIT`, `PZ_CRPS`, `PZ_ZSCORE`, `PZ_DZ`, `PZ_INCI68` (z_spec inside the 68 % interval), `PZ_OUTLIER`; headline numbers go to the catalog metadata (`photoz_quality`).

## Step *Spec-z Sample Representativeness* (`--mode repr`)
Is the spectroscopic (training / calibration) sample representative of the target sample (all catalog objects, or `--target nospec`) in the features `--features`
(columns or colours `A-B`; default MAG_AUTO and the colours between the MAG_* columns; robust-scaled by the target median / MAD)?
* per feature: KS, standardised mean difference, quantiles, fraction of target objects outside the spec range (`pzr_features.tsv`);
* multi-dimensional k-NN two-sample test: AUC of the leave-one-out k-NN "spec vs target" score with a permutation p-value (AUC 0.5 = indistinguishable);
* support: distance of every target object to its nearest spec object in units of the p95 (`--coverage-q`) of the spec nearest-neighbour distances (expected outside fraction 5 % if representative);
* k-NN density-ratio weights (Lima et al. 2008, `--knn`) of the spec objects, the effective sample size fraction, weighted feature moments / KS; with a photo-z column the **weighted accuracy**
  (unweighted vs weighted bias / NMAD / outliers) is the estimate for the target sample;
* a verdict (`representative` / `NOT representative`), `pzr_report.json`, `pzr_plot.png`; columns `PZ_SPECDIST` (distance ratio, all objects; > 1 = outside the spec coverage), `PZ_INSPEC`, `PZ_SPECW` (weights, spec objects).
Feed `PZ_SPECW` back through `--weights-col` of the calibration step for weighted accuracy.

CLI: `python plugins/photoz_sed/photoz_quality.py --catalog CAT --work DIR --mode quality|repr [...]` (prints the new columns; `--meta-out catalog_meta.json`).
GUI: *Photo-z and SED* tab (`pq-*` parameters), viewer steps *Photo-z Calibration Plot...* / *Report...*; session export replays the steps headless.
`ds9_photo_z.py --mixture-out FILE` (parameter `photoz-mixture-out`) writes the MDN mixtures.

## Fixes made on the way
The MDN branch of `ds9_photo_z.py` never worked (it passed an array to `extract_features`, which wants dictionaries, and the feature order was per band while the checkpoint was trained on blocks
`mag_ugriz, magerr_ugriz, flag_ugriz, colours, petroR50_r`): every run silently fell back to the empirical relation.  Now `photo_z.features.extract_features_block` builds the training layout and
`predict_photoz` pads the columns the catalog cannot supply (petroR50_r) with the training mean (`return_mixture=True` returns pi / mu / sigma).

## Validation (`photoz_synth_report.json`, `photoz_real_report.json`)
* **Synthetic truth** (Gaussian errors with object-dependent sigma, 200 trials of 1000 objects): calibrated PDFs are rejected at the nominal rate (KS 4 %, CvM 5 %, chi2 1.5 % at p < 0.05, 68 % coverage 0.681);
  widths x0.8 / x1.25 and a bias of 0.01 (1 + z) are rejected in 100 % of the trials (coverage 0.573 / 0.789 / 0.646); 3 % catastrophic outliers are *not* caught by the PIT tests at n = 1000
  (8 % KS) but raise the mean CRPS by 35 % and show in the outlier fraction (expected 0.0300, measured 0.0295); the optimal width scale recovers the applied scale to 0.3 % (0.6 ... 1.5).
* **Real SDSS (10 000 objects, repository MDN, held-out test split of the training script, 1500 objects)**: sigma_NMAD 0.0235, bias +0.0018, outliers 1.27 % (p95 |dz| 0.076),
  PIT KS p 5e-4, CvM 9e-4 (mean 0.482, var 0.0755 vs 0.0833: slightly under-dispersed on the tails -> 68 % coverage 0.726, 95 % 0.956), CRPS 0.0261; optimal width scale 0.94 (log score gain tiny), PIT
  recalibration cross-fit KS p 0.999; deliberately scaled PDFs (x0.5 / x2) are rejected (KS p 7e-26 / 3e-43, coverage 0.46 / 0.93); Gaussian N(z_point, z_err) instead of the mixtures: KS p 3e-7.
  NMAD grows from 0.014 (r < 16.5) to 0.030 (r > 18.5).  End-to-end through `ds9_photo_z.py --mixture-out` from raw SDSS magnitudes (petroR50 padded): median |z_pipeline - z_direct| 0.003 (max 0.059),
  NMAD 0.0253, KS p 0.049.
* **Representativeness**, truth = the whole 3000-object held-out set (NMAD 0.0235, outliers 2.0 %, p95 0.083): random half -> AUC 0.496, p = 1.0, 2.7 % outside the support (expected 5 %); soft brightness selection
  (p = clip(10^(-0.3 (r - 17.2)), 0.05, 1), 1193 objects): AUC 0.73, 19 % outside support, unweighted NMAD 0.0170 / outliers 0.59 % -> kNN weights k = 5 / 10 / 20: 0.0193 / 0.0203 / 0.0204 and 0.89 / 1.29 / 1.39 %
  (ESS 0.67 / 0.48 / 0.38); the second soft selection (827 objects) 0.0162 -> 0.0187 / 0.0198 / 0.0200.  A hard cut r < 17.2 (470 objects): AUC 0.91, **74 % of the target outside the spec support**, the weights cannot help
  (NMAD 0.0133 -> 0.0136): the step says so (verdict NOT representative, `PZ_INSPEC = 0`) instead of fixing it.
* GUI (M51 catalog, toy EAZY photo-z against a synthetic spec-z file): `verify_newplugins.tcl` section `photozq`, 24 checks, session replay identical.

## Limitations
* The PIT tests are insensitive to a few per cent of catastrophic outliers at a few hundred objects (use the outlier fraction, CRPS, `frac_extreme`); PIT uniformity is necessary, not sufficient (a PDF that ignores the data can be PIT-calibrated; use CRPS / log score to compare models).
* The calibration is only as good as the spec-z sample: spec-z are biased to bright, low-z, emission-line / red-sequence objects, and a spec-z failure is an "outlier"; check the representativeness first.
* k-NN weights cannot extrapolate beyond the spec support; the effective sample size drops quickly with the strength of the selection; weights are of the features chosen (not of z).
* Gaussian mixtures are used untruncated for CRPS / z-score (mass below `--zmin` is renormalised in the PIT and log score only).
* The petroR50_r input of the shipped MDN is not in a catalog; it is replaced by the training mean (the pipeline NMAD is close to the direct one, the calibration is not identical).
* One predictive distribution family per run; no per-object calibration (PIT recalibration is global, optional and only reported, not applied to the catalog PDFs).
