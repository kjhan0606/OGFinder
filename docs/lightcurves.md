# Transient light curves (plugin `lightcurves`, tab Measure, chip "Light curves")

Classifies transient light curves and gives a calibrated supernova score; sits on top of the Moving-objects pipeline (difference images → Transient Candidates → Light Curve) but also reads any long-form table.
Engine `ogfkit/lcclass.py` (features, Bazin fit, pure-numpy random-forest predictor), synthetic library `ogfkit/lcsynth.py`, CLI `plugins/lightcurves/lightcurves.py`, reproducible training `plugins/lightcurves/train.py`,
shipped models `model_lc.json` / `model_host.json` and `training_report.json`, demo data `make_demo.py`.

## Input
* **file**: TSV/CSV, one row per epoch: `id` (or name/object), `mjd`, `flux`+`flux_err` or `mag`+`mag_err` (AB); optional `number` (catalog NUMBER), `ra`/`dec`, `host_offset_re` (offset from the host in effective radii). Fluxes are put on the zero point 25 scale (parameter *zp*).
* **moving**: the Moving-objects work directory (`lightcurves.json` from *Light Curve*, `transients.tsv` from *Transient Candidates*: host id, positions, `offset_re`). Per-epoch zero points are recovered from the detected points (mag + 2.5 log10 flux).
* Linking to catalog rows: `number` column / `host_id`, else numeric id, else nearest catalog object to ra/dec within the link radius. Several light curves on one object: the highest SN score is written. Light curves without a catalog counterpart stay in `lc_results.tsv`.

## Output
Step **Classify Light Curves** adds `LC_ID LC_N LC_CLASS LC_PCLASS LC_PSN LC_PEAK_MJD LC_PEAK_MAG LC_RISE_D LC_DECL LC_CHI2R`. `LC_CLASS` ∈ SNIa, SNIbc, SNII, fast, AGN, variable, static, or `unclassified`
(fewer than *min points*, peak S/N below *min S/N*, or top probability below *min probability*; the best guess stays in `lc_results.tsv`). `LC_PSN` = summed probability of the three SN classes.
`LC_RISE_D` is given only when the first epoch is below half the peak; `LC_DECL` is the decline rate (mag/day) within 25 d after the peak. `lc_results.tsv/json` have one row per light curve with all class probabilities. Windows: *Light-curve Viewer…* (flux plot with Bazin fit, Prev/Next) and *Result Table…*.

## Method
21 features per light curve: baseline, number of points, peak S/N, χ² against constant flux, von Neumann η, Bazin fit (rise/fall time, peak epoch, χ², improvement over constant), half-maximum duration, rise/decline rates, skew, kurtosis, fraction beyond 1σ, percentile ratio, rank trend, Lomb–Scargle power and period, maximum slope.
Random forest (50 trees) trained on the synthetic library; the host model adds the log host offset in R_e. Training: `python3 plugins/lightcurves/train.py --per-class 600 --seed 20260101` (4200 light curves; the same seed gives byte-identical models — tested).
The library (`lcsynth.py`) is a set of parametric caricatures: SN Ia / Ibc (Bazin shapes, time-dilated), SN II (plateau + drop, or linear decline), fast transients (rise ~1 d, fade 1–4 d), AGN (damped random walk), periodic variables, static residuals; random cadence (6–40 points over 15–220 d), depths 23.8–26.3 AB.

## Validation (tests in `plugins/lightcurves/tests`, numbers measured on held-out light curves, seed different from the training seeds)
* 280 held-out synthetic light curves: 7-class accuracy 0.746 (light curve only) / 0.807 (with host offset); per class (lc only) SNIa 0.47, SNIbc 0.68, SNII 0.78, fast 0.78, AGN 0.60, variable 0.97, static 0.95 — the three SN types are only partly separable (SN Ia is confused with Ibc/II).
* SN vs non-SN (P_SN ≥ 0.5): completeness 0.917, purity 0.948 (lc only); 0.933 / 0.957 with host. Calibration of P_SN: expected calibration error 0.042 (5 bins).
* Accuracy grows with S/N (0.70 below peak S/N 15, 0.83 above 30). The host offset lifts AGN recall from 0.60 to 0.93.
* **Injection through the existing difference-image photometry**: 70 transients (10 per class) injected as PSF-convolved sources into 16 simulated difference images (σ_PSF 1.5 px, depth 25.3 mag), light curves measured with `moving.transients.lightcurve` (3-px aperture; recovered/injected flux 0.848, aperture fraction 0.865), written as Moving-objects files and classified through the CLI: 7-class accuracy 0.77, SN selection completeness 0.93, purity 0.88.
* Shipped models: held-out accuracy 0.783 (lc) / 0.846 (host) on 1050 light curves of the training script's own validation split.
* Readers (flux/mag, zero point), link modes (number, position), plot PNG (640×340), reproducible training: tested.

## Limits
* The shipped default `model_lc.json` is trained on synthetic caricatures, **not on real survey light curves** (real-data validation and an optional ZTF-trained model: section "Real ZTF validation" below): the numbers above are self-consistent but real SNe have colour evolution, diversity (91bg, IIn, IIb, SLSN, TDE, novae, …) and non-Gaussian errors that are not represented; expect lower real-world accuracy. Retrain with `train.py` after replacing `lcsynth.py` by real templates (e.g. SNANA/PLAsTiCC-style) or labelled data.
* Single band only (no colour); no redshift prior or absolute-magnitude consistency check; host information is only the offset in R_e; AGN vs SN confusion remains for sparse light curves.
* Non-detections are not used (only the forced-photometry fluxes); errors from the Moving-objects step contain background scatter only (no Poisson term) — fine for difference images, optimistic otherwise.
* P_SN is calibrated on the training library; with a different cadence/depth mix the calibration will drift.

## Real ZTF validation (round 2, item 3)
Data: `fetch_ztf.py` (network; public services only) downloads real ZTF alert-detection light curves (ALeRCE API, difference-image PSF photometry, one band per
object = the one with more detections, cut to a 220-day window so that the duration is not a class feature).  Labels: **SNe = SPECTROSCOPIC types** of the
ZTF Bright Transient Survey (VizieR J/ApJ/895/32, Fremling et al. 2020; Ia/91T/91bg/02cx... -> SNIa, II/IIb/IIn/87A -> SNII, Ib/Ic/Ib-c/Ic-BL/Ibn -> SNIbc; SLSN and
"ambiguous" skipped).  **AGN and periodic variables are labelled by the ALeRCE lc_classifier at p >= 0.9 (a machine classifier, not spectroscopy)**, so those two
classes measure agreement with ALeRCE, not truth.  1186 objects requested, 941 light curves with >= 8 detections, 770 with peak S/N >= 3 (SNIa 307, SNII 99, SNIbc 20, AGN 67, variable 277).
`train_real.py ZTF_LCS.json` evaluates on those objects only (5-fold stratified CV, `real_validation_report.json`):

| model | accuracy | balanced accuracy | SN-vs-rest AUC | SN completeness / purity (P_SN >= 0.5) |
|---|---|---|---|---|
| A shipped synthetic-trained `model_lc.json` (no retraining) | **0.604** | 0.541 | 0.964 | 0.967 / 0.920 |
| B random forest on real training folds only | 0.848 | 0.709 | 0.992 | 0.984 / 0.931 |
| C synthetic library (400/class) + real folds (weight 5) | 0.851 | 0.643 | 0.989 | 0.969 / 0.941 |
| B / C with `n_obs`, `log_baseline` removed (cadence/duration confound check) | 0.832 / 0.838 | 0.692 / 0.633 | 0.986 / 0.986 | 0.979 / 0.925, 0.967 / 0.934 |

Reading: the synthetic model transfers well for **SN vs. not-SN** (AUC 0.96, similar to its synthetic 0.92/0.95 completeness/purity) but the 7-class accuracy drops
from 0.78 (synthetic hold-out) to 0.60 on real data: SNIa -> SNIbc (127 of 307), SNII -> SNIa/SNIbc (64 of 99), AGN with low amplitude -> "static" (17 of 67);
accuracy grows with peak S/N (0.41 below 10, 0.64 for 10-30, 0.70 above 30).  Training on real objects raises accuracy to 0.85, but the gain is mostly the majority classes: SNIbc
(20 objects) is recognised 3 times of 20 even with real training data, and the SNIa-vs-SNII split (B: 28 SNIa called SNII, 25 SNII called SNIa) stays uncertain.
Limits: one band, ALeRCE detections only (no non-detection / forced photometry), BTS is magnitude-limited (r < 18.5 peak) and therefore biased to bright, low-redshift
SNe (the plugin's targets are fainter), class labels for non-SNe are ALeRCE's, no `fast`/`static` real examples (their real recall is untested), small SNIbc sample.
`ztf/model_lc.json` (model C trained on all 770 + synthetic library, 7 classes kept; use `lightcurves.py ... --model-dir plugins/lightcurves/ztf`; no host model) is
**not the default**: its numbers above are cross-validated, but BTS-like bright SNe are not the typical input of the Moving/Transient steps; the default stays the synthetic model so
existing golden outputs do not change.  `tests/data/ztf_sample.json` holds 20 of the real objects for the loader tests.
Reproduce: `python3 plugins/lightcurves/fetch_ztf.py DIR && python3 plugins/lightcurves/train_real.py DIR/ztf_lcs.json [--write-model M.json]`.  Other surveys: any JSON list with
`label, t, mag, err[, sign]` or `label, t, flux, err` is accepted by `train_real.py`; PLAsTiCC-style tables must be converted to that layout (not done: PLAsTiCC is simulated data).
