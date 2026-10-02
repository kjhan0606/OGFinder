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
* Trained on synthetic caricatures, **not on real survey light curves**: the numbers above are self-consistent but real SNe have colour evolution, diversity (91bg, IIn, IIb, SLSN, TDE, novae, …) and non-Gaussian errors that are not represented; expect lower real-world accuracy. Retrain with `train.py` after replacing `lcsynth.py` by real templates (e.g. SNANA/PLAsTiCC-style) or labelled data.
* Single band only (no colour); no redshift prior or absolute-magnitude consistency check; host information is only the offset in R_e; AGN vs SN confusion remains for sparse light curves.
* Non-detections are not used (only the forced-photometry fluxes); errors from the Moving-objects step contain background scatter only (no Poisson term) — fine for difference images, optimistic otherwise.
* P_SN is calibrated on the training library; with a different cadence/depth mix the calibration will drift.
