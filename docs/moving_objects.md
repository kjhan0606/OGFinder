# Moving Objects: asteroids, transients and reference images in OGFinder

*Moving Objects* is a plugin (`plugins/moving/`, chip **Moving** on the **Time-domain** tab of the workflow panel) that runs a multi-epoch pipeline on HST/ACS-type
`*_flc.fits` exposures: fetch references, align, difference, detect, link tracklets, identify known objects, fit an orbit,
find static transients, make light curves and export. All numerics are in the `moving/` package and are driven through one
CLI, `ds9/library/ds9_moving.py`; `plugins/moving/moving.tcl` is only the plugin's dialogs, kind registration for the shared
table, markers and the details window.

**What existed before.** There was no asteroid mode in this repository. It had stock DS9 SkyBoT access (`catskybot.tcl`) and a
wish-list section about difference imaging in `docs/manual/main.tex`. Everything described here is new.

## Menu (exact labels) and the shared table

Chip menu **Moving** on the Time-domain tab (it replaced the old top-level *Moving Objects* menu) → `Fetch Reference...`,
`Align`, `Difference`, `Detect`, `Link Tracklets`, `Identify Known Objects`, `Orbit Fit...`, `Transient Candidates...`,
`Light Curve`, `Export...`, (separator) `Select Exposures...`, `Time-domain details...`, `Work Directory...`.

There is **no separate results window** any more (the old `Show Results Table` entry is gone).  Results are shown in the
one catalogue table, which has a `kind` filter (radio buttons *Galaxies / Moving / Transients / Detections / All* under the
tab chips).  Sorting, filtering, row selection (table row <-> marker) and *Save* behave as for galaxies.  Moving objects
get an `overlap` / `overlap_id` column (a galaxy of the current catalogue within the tracklet's footprint), transients get
`host_id` / `host_sep` (nearest catalogue galaxy within 5 arcsec).  Marker colours: moving = magenta, transients = red,
detections = white.  Markers are drawn in every tile in tile mode.

**Time-domain details** is ONE detachable window (`.ogftd`) with the tabs *Orbit* and *Light curve*; it is created once
and reused.  Selecting a moving row shows the Orbit tab, a transient row the Light curve tab.  The CLI argv of each step is
unchanged.

*Fetch Reference...* has a provider radio button **MAST / LSST / both**; the choice is written to the CLI as
`--provider mast|lsst|both` and therefore appears in the session-recorder step arguments.

## Methods: implemented vs simplified

| step | implemented | simplification / caveat |
|---|---|---|
| Fetch | MAST (astroquery, public), Pan-STARRS cutouts, SkyBoT, Horizons | see LSST table below |
| Align | Gaia DR3 (TAP) matching when stars exist, else relative chain with a polynomial/shift+rotation model between exposures | in the test field no Gaia stars matched; result is relative only (4–19 mas rms), the anchor exposure has no absolute tie |
| Difference | ZOGY (Zackay, Ofek & Gal-Yam 2016) in Fourier space; template = median of the other aligned exposures. `zogy.zogy` also returns `alpha_new`/`sigma_alpha_new` = flux in the *target exposure's* flux scale (alpha * Fn / Fr), which `difference_chip` uses. **Optional, off by default** (see "Difference-image quality" below): source-noise variance maps `Vn`/`Vr` (eqs. 26-33 of the paper), astrometric-registration variance `astrom_n`/`astrom_r` (`S_corr`, `V_S`), per-tile empirical PSF (`imaging.measure_psf_field`, `zogy.zogy_tiled`), separate target/template PSFs | default call = sky-noise-only normalisation exactly as before; one PSF per chip unless a `PSFField` is passed; on the BB89 data the PSF field falls back to a constant PSF (see below). A PSF 20% too narrow gives ~22% low flux (`moving/tests/test_zogy.py`) |
| Detect | thresholding on the S/N image (since 2026-10-02 + L.A.Cosmic features, trailed-Gaussian centroid refit, real/bogus score - see "Detection stage" below), trail (elongated) channel, classes `point trail artefact_cr artefact_edge artefact_static artefact_dipole negative faint` | the archive CR mask is dense and mislabels real detections as `artefact_cr`; `artefact_static` = non-trail detection with template S/N > 8 |
| Link | triplet seeds (first / middle / last exposure) with an analytic parallax factor k from the middle exposure (k grid when the baseline is degenerate), KD-tree extension to the other exposures, vectorised 5-parameter weighted least squares (position, rate, k), rate and k gates, chi2-like score (stage 1). **Precision stage (2026-10):** detection veto masks before the pool cap (`pipeline.detection_veto`: stationary in all other exposures, template residual, chip edge), a tracklet-level logistic score (`tracklet.LLR_MODEL`: log rate, members, archive-CR fraction, sharpness, trail members, min S/N, chi2) that becomes `score` (= minus logit) with `prob`, a Sun-bound-orbit cut (`bound_orbit_ratio`), and clustering of near-identical (position, apparent velocity) candidates; a candidate sharing >= 2 detections with a better one is still dropped. The previous pair-based linker is kept as `tracklet.link_exposures_legacy`; `rescore=False` reproduces the stage-1 ranking. **Still no multi-night heliocentric clustering** (`helio_linc` exists but is not wired into the linker) | see the injection benchmark below |
| Identify | SkyBoT candidates re-checked with Horizons using the true HST observer (`500@-48`) | |
| Orbit | statistical ranging over the admissible region (Granvik et al. 2009 / Muinonen et al. 2010 style, simplified), 2-body polish, then ASSIST/REBOUND differential correction (DE440 + 16 massive asteroids + GR) with Carpino et al. (2003) outlier rejection, Eggl et al. (2020) star-catalogue debiasing, default station sigmas | no classical Gauss IOD (not implemented); station weights are defaults, not the Veres et al. (2017) table; solar light deflection and stellar aberration neglected; `mcmc_posterior` (emcee) exists but is untested; arcs shorter than ~1 day are flagged `determined=False`: elements are blanked, only class probabilities and ranging quantiles are shown; **classification is arc-aware** (below) |
| Transients | static grouping across ≥2 files (CRs cannot repeat), host association with any OGFinder/SExtractor TSV (separation, in R_e, z), heuristic labels | heuristic only, not a classifier; TNS crossmatch is optional via an environment variable and not tested (the site returned 403) |
| Light curve | forced aperture photometry on difference images, AB magnitude from the header zero point | aperture only |
| Photometry | HG phase function (Bowell et al. 1989), absolute magnitude H, diameter with an *assumed* albedo | |
| Export | MPC 80-column optical lines (`C` and `s` satellite lines) and JSON, **local files only; nothing is ever submitted to the MPC** | |

## LSST / Rubin provider: what is actually public (probed 2026-10-01)

| endpoint | result |
|---|---|
| `https://data.lsst.cloud/api/tap/availability` | 200, 292 B |
| `https://data.lsst.cloud/api/hips/v2/dp1/list` | 200, 13301 B (metadata only) |
| `https://data.lsst.cloud/api/hips/v2/dp2/list` | 200, 7982 B (metadata only) |
| `.../api/sia/dp1/query`, `.../api/tap/sync`, `.../api/datalink`, DP1 HiPS properties/tiles | **401** (RSP login required) |
| CDS HiPS `Rubin/CDS_P_Rubin_FirstLook` (`Moc.fits`) | 200, 17280 B; colour PNG press art, not science data |
| MAST | no LSST collection (0 rows) |

DP1 is available through the Rubin Science Platform to data-rights holders only, and covers seven ~1 deg² fields (47 Tuc,
Low Ecliptic Latitude, Fornax dSph, ECDFS, EDFS, Low Galactic Latitude, Seagull). The HUDF position is inside ECDFS but
login-gated; the COSMOS asteroid test field (150.14, +2.36) is outside DP1. The provider therefore reports
`provider_unavailable` (with the reason; example output in `moving/validation/lsst_probe.txt`) and falls back to
MAST/Pan-STARRS (`MAST(fallback)` label). No credential is requested or stored. A user with data rights can download DP1 FITS
themselves and pass `--lsst-local-dir`; that path was **not** tested with real DP1 files.

## Files

* `moving/` Python package (`tests/` pytest suite, `validation/` the scripts used for the numbers below)
* `ds9/library/ds9_moving.py` CLI (`--mode setup|fetch|align|difference|link|nightlink|identify|orbit|transients|lightcurve|export`);
  prints `#MOVING {json}` status lines; `plugins/moving/moving.tcl` GUI (plugin; manifest `plugins/moving/plugin.json`)
* `scripts/run_moving_pipeline.py` standard-library driver that runs align → difference → link → identify and writes a
  manifest with the command lines, package versions and sha256 of inputs and outputs
* Each GUI step is logged with `OGFSessLog` (one argv, declared outputs) when the session recorder is present. Table selection
  and region loading are display-only and are not recorded. Replay through the recorder was not tested.

## Validation (real data; numbers from this box)

Test field: HST ACS/WFC F814W, COSMOS, 2004-04-19, `j8pu38c7q/caq/ceq/ciq` (MAST obsid 24820553 for `j8pu38010`), containing
2015 BB89 (558706, V 21.5, 14.5 arcsec/h geocentric-equivalent).

* Detection: the mover is detected in all 4 exposures, 0.05/0.16/0.27/0.07 arcsec from the JPL Horizons position (observer
  `500@-48`); S/N 760/777/591 in the trail channel for the first three exposures.
* Linking: after raising the pool to 900 detections per exposure the true tracklet is found (rms 0.03 arcsec, rate 9.2 arcsec/h
  in the parallax model, not the geocentric 14.5) as tracklet id 19, i.e. **rank 20 of 400** by score (an earlier run with a smaller pool had it at rank 9; ranking is not stable), and `Identify`
  marks it known (2015 BB89, 0.27 arcsec). It is not top-ranked and a second, partially wrong tracklet (id 37) also matched.
  Rates differ from the Horizons geocentric rate because the fit includes the parallax hypothesis; this was not reconciled.
* Orbit from this 4-exposure, 0.024 d arc: not determined (as expected). Ranging gives class probabilities (Centaur/JFC 0.82,
  Hilda 0.17; ranging a 16/50/84% = 5.2 / 7.8 / 10.4 AU) while the true object is a main-belt asteroid (a = 3.17 AU, e = 0.106,
  i = 19.0 deg).  Before: `orbit_class` was the wrong top class.  Now (`kepler.classify_arc`): `orbit_class = "unclassified (arc 0.0244 d too short)"`,
  `class_status = unclassified`, the probabilities are reported only as `orbit_class_probs` and the best-fit class as `orbit_class_bestfit`.
* Injection–recovery, original linker (30 synthetic Gaussian-PSF movers in real chips, mag 21–26.5, 1.5–40 arcsec/h geocentric): 30/30 detected in ≥3
  exposures as single detections, but only 8/30 linked correctly (all with mag 21.4–22.8; none at mag ≥ 23.0). Fitted rates match the
  injected ones (e.g. 5.3/5.3, 11.7/11.7, 8.3/6.8 arcsec/h; one object at 32.8 arcsec/h was fitted as 74.9). The 13 correct tracklets (covering the 8 linked objects; `rank4.py`) ranked at positions 16–284 (list order, sorted by n then score) of the 400 returned tracklets (score 6.7–11.8); most returned tracklets are false. Injected sources use a Gaussian PSF, not a real ACS PSF.
* Injection benchmark of the new linker (`python moving/validation/link_bench.py injN.json.pkl [--legacy] [--max-tracklets 400]`;
  six injected sets inj1..inj6 = 10 + 5x30 movers in the real BB89 chips, tolerance 1.0 arcsec, recovered = a returned tracklet has
  >= 3 members within 1 arcsec of the truth; "ceiling" = injected objects with >= 3 detections in the linker's input pool, i.e. the
  most any linker can recover):

  | set | injected | legacy (top 400) | new (top 400) | new (top 3000) | ceiling | legacy / new time |
  |---|---|---|---|---|---|---|
  | inj1 | 10 | 1 | 4 | 5 | 7 | 95 s / 7-9 s |
  | inj2 | 30 | 5 | 15 | 17 | 24 | 54 s / 3-5 s |
  | inj3 | 30 | 2 | 6 | 7 | 15 | 136 s / 11 s |
  | inj4 | 30 | 8 | 8 | 8 | 20 | 134 s / 10 s |
  | inj5 | 30 | 7 | 10 | 12 | 21 | 136 s / 13 s |
  | inj6 | 30 | 10 | 10 | 10 | 18 | 131 s / 11 s |
  | total | 160 | **33** | **53** | 59 | 105 | ~10-15x faster |

  The "8/30" case above is inj4: the new linker also gets 8 there (no gain; its ceiling is 20, the remaining 12 are lost at the
  detection/class selection or in the 400-tracklet cap).  Honest limits: recovery is bounded by detection; precision is low (about
  390 of the 400 returned tracklets are false in this crowded pool of ~2100 detections per exposure) so ranking matters
  (top 50: new 2/14/4/5/8/8 vs legacy 1/3/1/7/5/7); 3-member candidates cannot be validated because the parallax factor is free.
  Unit tests: `moving/tests/test_tracklet.py` (3 detections, missing exposure, no duplicates, unphysical rate, ranking, a 10-mover
  mini-injection with >= 8 recovered).
* **Precision / ranking improvement (item 1, 2026-10-01).** Same six sets, tolerance 1.0 arcsec, `max_tracklets=400`, same box.  "Stage 1" =
  the previous linker (`link_bench.py ... --no-rescore --no-veto`, identical to the table above), "new" = defaults of this commit.
  `top-N` = distinct injected objects recovered among the first N returned tracklets; `P>=0.5` = tracklets with logistic probability >= 0.5.

  | set | injected | ceiling | recovered stage1 -> new | top-10 | top-25 | top-50 | top-100 | tracklets with P>=0.5 (true/returned) | time stage1 / new |
  |---|---|---|---|---|---|---|---|---|---|
  | inj1 | 10 | 7 | 4 -> 5 | 1 -> 4 | 2 -> 4 | 2 -> 5 | 2 -> 5 | 1/2 | 9.1 / 10.0 s |
  | inj2 | 30 | 24 | 15 -> 19 | 8 -> 9 | 14 -> 15 | 14 -> 16 | 14 -> 17 | 15/21 | 3.5 / 4.3 s |
  | inj3 | 30 | 15 | 6 -> 7 | 3 -> 3 | 4 -> 7 | 4 -> 7 | 6 -> 7 | 2/4 | 10.7 / 12.3 s |
  | inj4 | 30 | 20 | 8 -> 8 | 3 -> 4 | 4 -> 8 | 5 -> 8 | 7 -> 8 | 4/5 | 10.0 / 11.4 s |
  | inj5 | 30 | 21 | 10 -> 15 | 4 -> 6 | 7 -> 12 | 8 -> 13 | 9 -> 14 | 8/16 | 10.1 / 10.8 s |
  | inj6 | 30 | 18 | 10 -> 12 | 7 -> 9 | 8 -> 10 | 8 -> 11 | 8 -> 11 | 9/10 | 10.0 / 10.7 s |
  | **total 1-6** | 160 | 105 | **53 -> 66** | 26 -> 35 | 39 -> 56 | 41 -> 60 | 46 -> 62 | 39/58 (precision 0.67) | 53 / 60 s |
  | inj7 (held out) | 30 | 20 | 11 -> 10 | 7 -> 5 | 8 -> 7 | 8 -> 8 | 9 -> 8 | 6/13 | 10.7 / 11.0 s |
  | inj8 (held out) | 30 | 18 | 7 -> 10 | 6 -> 7 | 7 -> 8 | 7 -> 9 | 7 -> 9 | 6/6 | 9.9 / 11.5 s |

  Held-out sets inj7/inj8 (new injections with other random seeds, **not** used to fit the logistic score; sets 1-6 were) give 18 -> 20 of 60: the
  gain on sets 1-6 (53 -> 66) is partly in-sample.  The leave-one-set-out estimate of the fitted model on candidate lists is in
  `validation/fit_link_score.py` output (sets 1-6: 7/7, 18/22, 9/11, 11/18, 15/19, 11/15 reachable objects in the top 400).  Honest summary of the
  generalisation: top-50 on the held-out sets 8 -> 8 and 7 -> 9, i.e. about +10 % there, not the +25 % seen in-sample.
  **Precision (true / returned) over all 400 returned tracklets is still 1.5-4.8 %** (390 of 400 false) because 400 tracklets are always returned and
  the crowded pool of ~2100 detections per exposure contains thousands of chance 4-point alignments (95320 candidates survive the gates in inj4 with a
  400-cap).  What improved is the *ranking*: at `prob >= 0.5` the list is 58 tracklets with 39 true (precision 0.67; per set 0.5-1.0, 0.46 on inj7),
  recovering 39 of 66 recovered objects; at `prob >= 0.2` and 0.05 precision is 0.2-0.4 and 0.05-0.14 (`validation/` sweep, the numbers are
  a calibration of this instrument/field type only).  The pipeline returns the same `max_tracklets` as before (CLI argv unchanged); the `prob`
  column is in `tracklets.json`.
  Ingredients and what each was measured to do (all numbers on sets 1-8, `moving/validation/`):
  * **Veto masks** (`pipeline.detection_veto`): *stationary* (same sky position <= 0.2 arcsec, S/N >= 8, in all other exposures) removes
    37 false and 1 injected pool detections of 62052; *template residual* (non-trail, template S/N > 4, a > 1.5 px) removes 386 false and 5 injected
    of 62052 (0.6 % of false pool, 0.3 % of true); *edge* (< 10 px from the chip border; in the CLI only, benchmark pickles have no chip
    shapes) removes 1807 false / 3 injected on sets 1-6.  Bright-galaxy cores are covered by the template rule; the `neg_frac` dipole test was
    not added (217 false vs 10 true detections hit by neg_frac > 0.03, no net gain).  The effect of vetoes alone on recovery is nil on this benchmark
    (66 -> 66 with `--no-veto`), they mostly cut the 600-per-exposure pool wasted on static residuals and help real crowded data; the "ceiling" is
    unchanged because the injected objects are rarely vetoed.
  * **Bound-orbit hypothesis**: the Sun-relative speed of a candidate at k = 1/Delta (k within +-2 sigma) must not exceed 1.25 x the local escape
    speed (only the sky-plane velocity is known, so this is the *minimum* possible speed).  On 4-exposure HST tracklets k is only known to +-1 /AU, so this
    cut is weak: it removes 434 of the 50000 rescored candidates in inj4, and keeps 90-97 % of the true candidates and 76-83 % of the false ones (a
    test with a plain threshold on the best-fit k).  It is physically right but costs one benchmark object: inj4 object 23 (mag 21.5, 32.8 arcsec/h at
    Delta 4.4 AU) was *injected* unbound (ratio 2.3); the linker finds it at rank 20 without the cut and drops it with the cut.  Without the cut
    inj4 recovers 9 and the total is 67 instead of 66 (`link_bench.py --no-bound`).  The injection generator draws rate and Delta independently, so 6 of the 30
    objects of inj4 are unbound by construction; real small bodies are not.  Not a clustering by heliocentric distance: with one HST orbit (0.025 d) there is
    no baseline for that (`helio_linc` needs >= 2 nights and is not wired into the single-exposure linker; the multi-night case is handled by `moving/nightlink.py`, below).
  * **Tracklet-level logistic score**: features and coefficients in `tracklet.LLR_MODEL`; the dominant terms are the fraction of members on
    archive-CR-flagged pixels (-10.2 per unit; 91 % of all S/N >= 8 detections are CR-flagged, but only 0.4 % of injected-object candidates have all
    members CR-flagged and no trail), the number of members (+4.8; 4 vs 3), a negative log-rate term (-0.83; slow objects are rarer to align by chance) and
    mean sharpness (-5.0).  Linearity (rms) and flux scatter did **not** improve it in cross-validation (tested, removed).  Cross-validated recovery with the
    chi2-like score alone was 55 (no dedupe, candidate list) vs 70-73 with the logistic features.
  * **Clustering in (position, apparent velocity)**: candidates within 1 arcsec and 3 arcsec/h of a better one are merged.  On the leave-one-out
    candidate lists: 71 -> 74 recovered, top-50 of the *kept* list 138 -> 125 true entries (duplicates removed); this is a small effect.
  * Dropped idea: de-duplication by "no shared detection" (instead of "no shared pair") loses objects (57 vs 71).
  * **inj4 (30 injected, ceiling 20): still 8 recovered (9 with `--no-bound`).**  Per-object diagnosis (`validation` scripts, truth positions): 10 of the 30
    objects have <= 2 detections in the linker pool (this *is* the ceiling gap: 20 reachable).  Of the 20 reachable, 8 are recovered in the top 400.
    The other 12 have their true tracklet at rank 52-93000 of the ~95000-candidate list (no-dedupe, stage-1 order; 3 more are not even in it): 18 of the 30 injected objects have at
    least one detected member 0.65-0.95 arcsec from the true position (trail centroids of faint/long trails and CR-contaminated point detections, 1 arcsec
    match radius), so the track has an rms of 0.25-0.5 arcsec and loses against thousands of chance alignments with a smaller rms; 6 objects are unbound by
    construction and 4 of them are at mag >= 25.  The logistic score reorders (top-10 3 -> 4, top-25 4 -> 8) but cannot lift tracks whose members are missing or
    mis-centred; the count does not change (8 -> 8).  Only detection-side work (better trail centroids, CR handling for faint movers) would help; out of scope here.
  * Real data: with `max_per_exposure=900`, tol 1.0, the true 2015 BB89 tracklet is rank 0 of 400 both with the stage-1 ranking and with the new
    score (prob 0.24); the earlier docs value (rank 20 with other settings) was not reproduced here, so no improvement is claimed on the real field.
  * Runtime: +0.5-1.7 s per set (the logistic features are vectorised over all gated candidates; the Sun-bound check runs on at most 50000).
  * Tests: `moving/tests/test_tracklet.py` (stage-1 reproducibility with `rescore=False`, probability monotone along the list, CR-flagged vs clean
    track, bound-orbit ratio and cut, cluster merge, veto rules and counts).  `test_chance_alignments_rank_below_a_real_mover` now asserts the strict
    first rank only for `rescore=False`; with the log-rate prior the mover must be within the top 3 in that uniform-clutter toy.

* Alignment: relative rms 19.4 / 4.4 / 3.7 mas (exposures 2–4 w.r.t. the anchor), no Gaia stars matched.
* Orbit fit pipeline check on an asteroid with a long observed arc (25153, public MPC astrometry, ASSIST force model, Eggl debiasing):
  * 90-day arc, 277 obs (100 used), rms 0.51 arcsec, χ²_red 1.15. Elements vs JPL SBDB, (fit − JPL)/σ_fit: a −1.4, e −0.9, i +0.5,
    Ω −0.6, ω −1.5, M +1.6.
  * 4-day arc, 39 thinned obs (13 used after rejection, 26 rejected), rms 0.44 arcsec, χ²_red 1.15. Deviations +1.2 σ (a), +1.2 (e),
    −1.5 (i), +1.4 (Ω), +1.7 (ω), −1.4 (M); acceptable but at the edge.
  * Same 4-day window through the shipped CLI (`--mode orbit --designation 25153 --mjd-min 58495 --mjd-max 58499`, all 115 observations, 37 used, no thinning): rms 0.56 arcsec, χ²_red 1.46; deviations +1.2 σ (a), +1.3 (e), −2.3 (i), +1.7 (Ω), +1.9 (ω), +2.4 (M), i.e. 2σ-level offsets with 78 observations rejected; this short arc is marginal.
  * A bug was found and fixed in the ASSIST propagator wrapper (one simulation reused for forward and backward sweeps): before the fix
    a differential correction started from the JPL state gave 57 arcsec rms, after the fix 1.76 arcsec (χ²_red 0.51) with default weights.

## Detection stage (item 1, 2026-10-02): cosmic rays, trail centroids, real/bogus

What changed (all in `moving/`, all ON by default in `detect.difference_chip`; `cr_reject=False, trail_fit=False, realbogus=False`
reproduces the previous detections; CLI: `--no-cr-reject --no-trail-fit --no-realbogus`, default argv unchanged):

* `crrej.py` - L.A.Cosmic (van Dokkum 2001) through `astroscrappy` (optional dependency, `pip install astroscrappy`; without it the numpy
  `imaging.lacosmic_mask` is used).  The result is a **feature** (`lac3`, `lac_n7`, trail `lac_frac`), not a veto: on a BB89 chip the LAC mask covers 1.1 %
  of the pixels (archive flag 1.13 %), 87 % of LAC pixels coincide with archive flags and 84-87 % of archive-flagged pixels are found by LAC.
  On inj4, 124 of 191 injected non-trail detections lie on LAC pixels, so a hard veto would delete most faint movers.
* `centroid.py` - trailed-Gaussian model (closed form, soft-L1 least squares, CR / bad-pixel mask) re-fits the position and length of trail
  detections up to `trail_fit_max=40` per chip.  Median centroid error against the injected truth (inj4): trails < 25 px
  0.28 -> 0.08 px; 25-50 px 5.8 -> 1.6 px; **> 50 px (> 2.5 arcsec) 18.9 -> 17.7 px: no improvement, 5 of 10 are still off by > 14 px**.
* `realbogus.py` + `data/realbogus_model.json` (150 depth-3 trees, evaluated in pure numpy, `export_sklearn` for inspection) - a
  gradient-boosted real/bogus score `rb` from 32 shape / ring / PSF-chi2 / LAC / trail features.  **Trained only on eight non-BB89 COSMOS
  visits with injected movers** (169,840 detections: 1,493 real, 168,347 bogus; the first version used four fields and 60 trees); leave-one-field-out
  AUC 0.975 / 0.986 / 0.981 / 0.994 / 0.997 / 0.990 / 0.982 / 0.982, recall at 0.2 % bogus kept 0.72-0.84.  Top features: `lac_n7`, `trail_peak_frac`, `n_hi`, `snr`, `b_pix`.
  **Held-out BB89** (injected movers in the real BB89 chips, never trained on; `validation/eval_realbogus.py`): plain sets AUC 0.991, recall 0.77 at 0.2 % / 0.93 at 2 % bogus kept;
  CR-heavy sets AUC 0.948 (0.41 / 0.65) - the old 4-field model gave 0.990 (0.78 / 0.91) and 0.942 (0.39 / 0.69), i.e. the gain from more fields is small.
  Hyper-parameter sweep (trees, depth): (60,3) LOFO AUC 0.983 | (150,3) 0.986 | (100,5) 0.984 | (200,4) 0.982 - chosen on the LOFO numbers, not on BB89.  Optional ZOGY
  improvements (source-noise, astrometric variance) gave no gain in this classifier and stay OFF.
  **Caveats:** labels are *injections into real images*, not real asteroids; the only real mover available is 2015 BB89.  `moving/tests/data/rb_heldout_sample.json`
  (60 real + 240 bogus held-out detections) backs a regression test (AUC 0.989, recall 0.93 at 2 %).  The previous unit test ordered two hand-made feature dicts,
  which are outside the training distribution (the 150-tree model scores them 0.025 vs 0.073); it was replaced by the held-out-sample test.
  Reproduce: `moving/validation/train_realbogus.py --make-sets DIR` then fit (`--trees N --depth D --seed S --out F`; `--oof` writes out-of-fold scores with the same settings).
* `pipeline.link_detections(use_rb, rb_min, rb_snr_min=6, use_lac_cr)`: with `rb` present the per-exposure pool is ordered by `rb` then S/N and the
  pool S/N floor drops to 6; without it the old S/N order is used.  `tracklet.LLR_MODEL_RB` (base features + mean / min member real/bogus
  logit) replaces the base logistic score when `rb` is present; fitted on out-of-fold `rb` of the eight non-BB89 training sets (the constants in `LLR_MODEL_RB` were refitted with the 150-tree model).
  A pinv fallback prevents a crash on degenerate (identical-epoch) fits.
* `detections.tsv` gained the columns `rb lac3 lac_n7 theta_pix trail_fit trail_len_pix` (appended; existing columns unchanged).
* `validation/inject.py` now writes exact `real` / `inj_box` labels and trail truth, supports `--nodq --cr-extra N` (CR-heavy sets: no archive CR flag
  and N synthetic CR tracks per chip), `--chips`/`--save-chips`; with no options its output is bit-identical to the old script for inj3-8
  (inj1/inj2 regenerate slightly differently, so the "before" numbers of those two differ from the older table above).
  `link_bench.py` prints per-group totals.

Benchmark on the regenerated sets (same box, `link_bench.py`, top 400, tolerance 1.0 arcsec; "before" = the unchanged pipeline on the same
chips, "after" = this commit).  inj1-6 are the sets the **base** tracklet score was fitted on (in-sample for that part), inj7-8 are held out from it, inj9-12 are new
CR-heavy sets (inj9-11 `--nodq --cr-extra 3000`, inj12 `--cr-extra 3000`).  The real/bogus model and `LLR_MODEL_RB` never saw any of them.

| set | injected | ceiling | recovered | top-10 | top-50 | precision | time before / after |
|---|---|---|---|---|---|---|---|
| inj1 | 10 | 6 -> 7 | 5 -> 5 | 4 -> 5 | 5 -> 5 | .013 -> .013 | 13.7 / 12.1 s |
| inj2 | 30 | 19 -> 22 | 10 -> 13 | 7 -> 7 | 9 -> 10 | .025 -> .033 | 15.9 / 11.6 s |
| inj3 | 30 | 16 -> 22 | 7 -> 10 | 4 -> 6 | 7 -> 9 | .018 -> .028 | 14.5 / 12.8 s |
| inj4 | 30 | 20 -> 21 | 8 -> 13 | 4 -> 8 | 8 -> 13 | .025 -> .035 | 14.6 / 11.2 s |
| inj5 | 30 | 21 -> 21 | 15 -> 19 | 6 -> 9 | 13 -> 19 | .040 -> .050 | 14.8 / 12.0 s |
| inj6 | 30 | 18 -> 21 | 12 -> 17 | 9 -> 9 | 11 -> 16 | .033 -> .045 | 14.4 / 14.4 s |
| inj7 | 30 | 20 -> 22 | 10 -> 16 | 5 -> 9 | 8 -> 15 | .030 -> .045 | 12.7 / 14.6 s |
| inj8 | 30 | 18 -> 25 | 10 -> 17 | 7 -> 9 | 9 -> 14 | .025 -> .045 | 13.5 / 11.5 s |
| inj9 | 30 | 6 -> 14 | 1 -> 5 | 0 -> 1 | 0 -> 2 | .003 -> .013 | 20.1 / 21.2 s |
| inj10 | 30 | 7 -> 14 | 0 -> 1 | 0 -> 0 | 0 -> 1 | .000 -> .003 | 18.5 / 22.6 s |
| inj11 | 30 | 7 -> 19 | 0 -> 4 | 0 -> 3 | 0 -> 3 | .000 -> .010 | 18.5 / 24.5 s |
| inj12 | 30 | 9 -> 16 | 0 -> 3 | 0 -> 1 | 0 -> 1 | .000 -> .007 | 18.9 / 21.4 s |
| **inj1-6** (160) | | 100 -> 114 | **57 -> 77** | 34 -> 44 | 53 -> 72 | | |
| **inj7-8** held out (60) | | 38 -> 47 | **20 -> 33** | 12 -> 18 | 17 -> 29 | | |
| **inj9-12** CR-heavy (120) | | 29 -> 63 | **1 -> 13** | 0 -> 5 | 0 -> 7 | | |

Honest reading: ranking and recovery improve on every group, but precision stays low (3-5 % of the 400 returned tracklets are true) and the CR-heavy
sets remain poor (13 / 120).  On those sets the per-exposure cap of 600 limits what is reachable (6 reachable at cap 600, 24 at 1200, 27 at 5000 in inj9)
and the median real/bogus score of real detections is only 0.09 there, i.e. the classifier is not reliable when the archive CR flag is absent *and*
thousands of CR tracks are added.  inj1-6 gains are partly shared with the in-sample fit of the base score; inj7-8 and inj9-12 are the cleaner evidence.
On inj1-8, 215 of 220 injected objects had >= 3 pool-eligible detections within 1 arcsec when the cap is removed: the remaining loss is cap and ranking.

Run time (4 chips of 1240 x 1240): detection 76 s -> 82 s with the CR features -> about 100 s with trail fits (320 fits, ~59 ms each); linking 11-15 s
(20-25 s on the CR-heavy sets).

**Real BB89 check (the known mover, 2015 BB89, CLI as above, `--max-per-exposure 900`):** before: rank 0 of 400, prob 0.24; after (4-field model): rank 0 of 400, prob 0.98 (8-field/150-tree model: rank 0 of 400, prob 0.993, identified as known, max separation 1.58 arcsec; link benchmark over injected sets 122 -> 123 recovered, per set within +-1, top-10 sum 67 -> 66) (3 members within 1 arcsec,
0.056-0.095 arcsec from the ephemeris positions; the 4th exposure's detection (0.07 arcsec away) has rb 0.016 and lies on a LAC pixel, so only 3 of
4 exposures enter the pool).  Nearest detection to the ephemeris position in exposures 0-3, before -> after: 0.049 -> 0.081, 0.166 -> 0.056, 0.268 -> 0.095, 0.071 -> 0.071 arcsec;
S/N unchanged (751 / 655 / 591 / 194).  One real object: this shows no regression, not a statistically meaningful gain.  Detection stage 164 s, link 38 s, identify 22 s;
the linker still reports `seed limit 4000000 reached` on this crowded field.

## Difference-image quality (ZOGY extensions)

Code: `moving/zogy.py` (`zogy(..., Vn, Vr, astrom_n, astrom_r)`, `zogy_tiled`), `moving/imaging.py` (`PSFField`,
`measure_psf_field`, `constant_psf_field`, `collect_star_cutouts`), `moving/detect.py` (`difference_chip(..., source_noise=,
astrom_sigma=, psf_field_t=, psf_field_r=, template_psf=, tile=, astrom_template_sigma=)`).  Every new argument defaults to the old
behaviour: `zogy()` called as before returns exactly the old keys and values.  With variance maps / astrometric sigmas it adds
`S_corr` (the corrected score image), `V_S`, `sigma_alpha_map`, `sigma_alpha_new_map`.  `astrom_sigma=True` takes the per-chip
registration rms from the alignment sidecar (`chip.align['rms_mas']`).  Since item 2 (2026-10-02) the options are wired into `pipeline.detect_in_region`, the CLI and the GUI (below); all default OFF, default argv unchanged.

Measured on synthetic Poisson star fields (bright stars + 25 injected 300-count transients, 4 seeds; `moving/tests/test_zogy.py` has the
corresponding assertions):

| quantity | classic (sky-only) | new |
|---|---|---|
| false positives > 5 sigma per image (no real variability) | 90-99 | 0-0.5 with source-noise terms |
| injected transients recovered > 5 sigma | ~100% | 98-100% |
| flux bias (the noise terms do not change alpha) | 0 to -12% depending on a 0-0.3 px astrometric offset, scatter ~8% | same |
| PSF FWHM_x varies 1.6 -> 3.6 px across the chip: flux bias left / right | +4.9% / -7.8% (constant PSF) | -0.6% / -3.2% (tiled PSF) |
| FWHM_x 2.0 -> 3.2 px: flux bias left / right | +1.9% / -2.3% | -1.4% / +0.5% |
| 0.3 px mis-registered 2e5-count star: dipole significance | 1 | < 0.2 (astrometric term) |

Real HST ACS BB89 data (cached `j8pu38{c7,a,e,i}q_flc.fits`, 4 chips of 1240x1240 px around 150.1375, +2.361; template = median of the
other exposures; script not in the repo, numbers from this box).  "pos" = positive detections with S/N >= 8 other than the known mover,
a *proxy* for false positives (the field is full of cosmic-ray residuals, moving objects and unmasked stars, so most are not noise
fluctuations); "mover" = number of difference-image detections of the known asteroid (7 with the classic setup):

| configuration | pos (S/N>=8) | negative >= 5 | mover hits | static-star residuals |
|---|---|---|---|---|
| classic | 7639 | 20 | 7 | 148 |
| source noise | 7870 | 7 | 6 | 143 |
| source noise + astrometry 0.5 px | 2021 | 2 | 4 | 36 |
| source noise + astrometry from alignment sidecar (19.4/4.4/3.7 mas) | 6442 | 4 | 4 | 109 |
| source noise + astrometry 0.3 px + tiled PSF | 2660 | 4 | 4 | 48 |
| tiled PSF only | 7648 | 20 | 7 | 147 |
| template PSF measured separately | 7639 | 20 | 7 | 148 |

Honest reading: (a) the astrometric term removes most of the residuals around bright static stars (static residuals 148 -> 36-48
with a 0.3-0.5 px assumption, 148 -> 109 with the sidecar rms, which at 0.05"/px is 0.39 / 0.09 / 0.07 px for exposures 2-4 and absent for the anchor exposure), which is its purpose, but
the large drop in "pos" for 0.3-0.5 px is mostly because it down-weights *everything* near sources in a crowded field, and it also
costs 3 of the 7 mover detections (the asteroid passes near stars in some exposures); (b) the tiled-PSF and separate-template-PSF
options changed nothing on this data because **no tile had enough field stars** (9-16 usable stars per chip with `min_stars` per tile
not reached in any of the 25 tiles: `tiles_from_stars = 0`) so every tile used the constant fallback, and the pooled empirical PSF has
FWHM 1.13 px, below the 0.085" (~1.7 px) floor, i.e. it is probably contaminated by cosmic-ray/hot-pixel "stars"; (c) the source-noise
term helps negatives (20 -> 7) but not positives here.  The synthetic gains for a varying PSF are therefore *not* demonstrated on real
data.

### ZOGY options in the CLI / GUI and the measured registration error (item 2, 2026-10-02)

* CLI (`--mode difference`; none of these appear in the default argv): `--source-noise`, `--astrom off|sidecar|measure|<sigma_pix>`, `--psf-tile N`,
  `--template-psf target|measure`.  GUI: *Tools > Plugin settings > Moving* (group "Difference (ZOGY)": Source-noise term, Registration error,
  PSF tile size and Template PSF as expert options); `OGFMovDifference` appends only non-default values, so the argv recorded by the session recorder
  is identical to the old one unless an option is changed (checked by `scripts/verify_moving_options.tcl`, 17 checks, part of `run_all_checks.sh` as `moving_options`).
* `--astrom measure` = `regerr.measure_registration`: sources bright in both the target and the template (stars and compact galaxies, CRs and
  movers drop out) are re-centred with the same windowed centroid on both images; sigma_x/y^2 = (median offset)^2 + (1.4826 MAD)^2 - centroid noise^2
  (centroid noise from size / S/N, template included).  It replaces the typed or sidecar number by a measurement for **each exposure/template pair**
  (the anchor exposure, which has no sidecar rms, gets one too).  If fewer than 8 usable sources it falls back to the sidecar value, else to no
  astrometric term (`info["registration"]` records what happened).  Synthetic test: known 0.30 px x-offset recovered (`moving/tests/test_regerr.py`, 6 tests).
* `--psf-tile N`: `pipeline.detect_in_region(psf_tile=N)` measures the PSF per tile from the stars of all exposures of that detector and runs `zogy_tiled`.

Real BB89 data (same four chips, 4 x 1240 px cutouts, template = median of the other three; script `/workspace/work/i2_eval.py`, not in the repo; the
detection-stage extras of item 1 are off in this table so the rows are comparable with the table above; the Gaussian PSF of FWHM 1.70 px is used because the
empirical one is unusable on this field).  Registration error measured per exposure (pixels, 0.05"/px): exposure 1: 0.23 (x) / 0.07 (y); exposure 2: 0.21 / 0.25;
exposure 3: no estimate (too few sources -> sidecar value 4.4 mas = 0.09 px); exposure 4: 0.0 / 0.14 (8 sources).  The sidecar says 19.4 / 4.4 / 3.7 mas =
0.39 / 0.09 / 0.07 px, i.e. the measured error on the *template pair* is larger than the alignment-fit rms for two exposures.

| configuration | pos (S/N>=8, not the mover) | negative >= 5 | mover hits (of 7) | static-star residuals |
|---|---|---|---|---|
| classic | 7707 | 26 | 7 | 151 |
| source noise | 7972 | 10 | 7 (S/N 760 -> 741 ... 193 -> 81) | 148 |
| source noise + astrometry from sidecar | 6654 | 6 | 4 | 116 |
| **source noise + astrometry measured** | 7427 | 6 | **5** | 133 |
| source noise + 0.15 px | 6965 | 4 | 4 | 125 |
| tiled PSF (256 px tiles; 0 tiles had >= 5 stars -> constant PSF) | 7714 | 25 | 7 | 148 |

The measured registration error does what it was meant to: it keeps 5 of 7 mover hits where the sidecar value keeps 4 (and a typed 0.15 px keeps 4) while still
removing the negative residuals (26 -> 6) and some of the static residuals (151 -> 133).  It does **not** recover all 7: the mover passes near stars and any
astrometric term that is large enough to reduce residuals also lowers S/N near sources.  The mover S/N drops from 193/193/187/89 to 70/57 in the faint exposures.
The pos(S/N>=8) proxy changes by only 4 % (it is dominated by CR residuals and unmasked sources, not by registration).
Linking benchmark with the options at detection time (`link_bench.py`, regenerated sets, defaults of item 1 otherwise; recovered / ceiling, top-50):

| set | default | source noise | source noise + measured astrometry |
|---|---|---|---|
| inj2 | 13/22, 10 | 15/22, 13 | 14/19, 14 |
| inj4 | 13/21, 13 | 11/20, 10 | 15/19, 15 |
| inj7 | 16/22, 15 | 13/19, 12 | 15/20, 13 |
| total | 42/65 | 39/61 | 44/58 |

Honest reading: 44 vs 42 recovered is within the set-to-set noise (3 sets, single noise realisation); the ceiling drops (65 -> 58) because the
astrometric term lowers S/N of real faint movers near sources, while top-50 is similar.  **Defaults therefore stay OFF**: the evidence supports "measured
beats sidecar / typed" for the astrometric term, not "on beats off" for movers.  Source noise alone does not help the linking benchmark (39 vs 42).
The tiled-PSF and separate template-PSF options again changed nothing on BB89 (no tile has enough stars); they stay untested on real data with a varying PSF.
Real BB89 end-to-end with `--source-noise --astrom measure` (CLI, as the earlier real run): the known 2015 BB89 tracklet is still rank 0 of 400, prob 0.995 (default: 0.983),
3 members 0.056-0.095 arcsec from the ephemeris; the fourth exposure's detection falls to S/N 44 (default 194) and is not linked; 8906 detections instead of 9808.

### Remaining limits

* Spatial PSF variation needs enough isolated stars per tile (default `min_stars` per tile); sparse fields fall back to one PSF.
* The astrometric term uses one scalar (or x/y) sigma per image, not a measured per-star residual map; it assumes the template and
  target noise are independent and stationary within a tile.
* Source-noise maps need a gain/readnoise-consistent variance estimate; `source_noise=True` uses the Poisson variance of the
  (smoothed) image, which is biased for the template of a median stack (correlated across exposures).
* Per-tile ZOGY (`zogy_tiled`) is seam-free only within the `margin` overlap; very strong PSF gradients across one tile are unmodelled.
* No colour term, no differential chromatic refraction, no handling of correlated noise from drizzle/resampling.
* The options are exposed in the CLI/GUI (next section) and are recorded in the session argv only when they differ from the default.

## Linking precision: orbit-population prior and known-asteroid comparison (precision batch item 4)

`moving/orbitlink.py` adds a precision stage to the tracklet linker for **fixed-observer** (ground-based) data: a synthetic population of bound two-body orbits
(JPL SBDB elements, `a` 1.5-5.5 AU, H <= 19.5, 1.42 M objects; random mean anomaly, real Earth state at the epoch, HG magnitudes, V < 22) is propagated to the
field to give the density of apparent motion vectors (east, north; arcsec/h).  A Gaussian kernel density (bandwidth 3 arcsec/h) defines a highest-density
region; tracklets whose rate vector lies outside the 99 % region are dropped.  No data of the field are used and the known objects are not in the prior.
Use: `P.link_detections(..., orbit_prior='auto'|elements.npz|RatePrior, orbit_prior_fraction=0.99)` or `ds9_moving.py --mode link --orbit-prior auto`
(default off, the HST path keeps its parallax-based `bound_orbit` cut).  `fetch_elements` downloads the table once (~145 MB query) to
`~/.cache/ogfinder_regression/sbdb_elements_H19.5.npz`.  Tests: `moving/tests/test_orbitlink.py` (synthetic elements, no network).

**Real-field comparison** (`moving/validation/sdss_known_asteroids.py`, reports `ska_c3.json`, `ska_c4.json`): SDSS run 94 (Stripe 82, 1998-09-19, near opposition),
fields 136-167, camcols 3 and 4, five exposures (r, i, u, z, g) 72 s apart, detections from the existing detector, linker `link_exposures` tile by tile in RA with
fixed observer, rates 10-80 arcsec/h, tolerance 0.8 arcsec.  Truth: IMCCE SkyBoT (observer 645) objects, matched by *predicted motion* (offset <= 3 arcsec,
tolerance 0.6 arcsec, moves >= 1.2 arcsec); the **ceiling** is the number of known objects detected in >= 3 bands.  Precision is a **lower bound** (an unmatched tracklet
can be an unknown real asteroid or a SkyBoT-faint object).

| field | ceiling (>=4 / 5 bands) | stage | tracklets | true | precision (lower bound) | recall |
|---|---|---|---|---|---|---|
| camcol 3 (39,525 det, 231 known) | 60 (33 / 19) | baseline linker | 342 | 49 | 0.14 | 0.82 |
| | | orbit prior 95 % | 70 | 47 | 0.67 | 0.78 |
| | | **orbit prior 99 %** | 83 | 49 | 0.59 | 0.82 |
| | | 99.9 % | 108 | 49 | 0.45 | 0.82 |
| camcol 4 (38,869 det, 222 known), nothing fitted | 69 (23 / 14) | baseline linker | 350 | 67 | 0.19 | 0.97 |
| | | 95 % | 88 | 63 | 0.72 | 0.91 |
| | | **99 %** | 110 | 67 | 0.61 | 0.97 |
| | | 99.9 % | 129 | 67 | 0.52 | 0.97 |

Kept tracklets by number of members (true / total, 99 %): camcol 3: 5 members 18/21, 4: 14/21, 3: 17/41; camcol 4: 5: 14/16, 4: 10/15, 3: 43/79 - the >= 4-exposure
tracklets are > 85 % real, the 3-exposure ones are where the remaining false links are.  Link time ~65 s per camcol.  The 99 % level was chosen on camcol 3 and applied
unchanged to camcol 4.  At the 99 % level the prior removes 76 % / 69 % of the baseline tracklets (camcol 3 / 4) without losing any matched one (49 -> 49, 67 -> 67); the population median rate at the field
(-31, -13) arcsec/h agrees with the known objects (-31, -14).

**Limitations / not done**
* This is *rate-prior vetting from an orbit population*: SDSS gives five exposures within 5 min, so no orbit can be determined there.  Cross-night orbit-fit linking now exists as a separate stage
  (`moving/nightlink.py`, next section) for data with tracklets on several nights; it is not applied to the SDSS fields (single night) and the rate prior remains the only vetting for single-night tracklets.
* Osculating two-body population with a in 1.5-5.5 AU and H <= 19.5: NEOs, comets, Centaurs/TNOs and objects fainter in H (but near enough to be detected) are not in the prior and
  would be dropped when their motion lies outside the region (SkyBoT shows none of the matched ones lost at 99 % in these two fields, but both are near opposition at low ecliptic latitude).
  The prior is specific to the field position and epoch (recomputed per call, ~tens of seconds with the 1.4 M table).
* Two fields of one night and one survey; truth from SkyBoT (MPC-known only); precision is a lower bound; the 3-exposure subset is still dominated by false links.

## Cross-night orbit-fit linking (`moving/nightlink.py`, R24)

**What it does.**  Takes tracklets (>= 2 detections within a night/visit) from different nights (and optionally different stations) and returns groups of >= 3 tracklets that are consistent with ONE Sun-bound two-body orbit, with the fitted
state, a, e, i, chi2/dof and rms.  CLI: `ds9_moving.py --mode nightlink --tracklet-files night1/tracklets.json night2/tracklets.json ... [--obs-code 500|250|...] [--nl-chi2 4] [--nl-floor 0.3] [--nl-min 3] [--nl-nbody]` -> `nightlinks.json`, `nightlinks.tsv`
(kind = linked | pair).  Python: `nightlink.Tracklet`, `link_nights`, `vet_with_links`; synthetic generator `nightlink_sim.py`; validation `moving/validation/nightlink_validate.py`; tests `moving/tests/test_nightlink.py` (run_all_checks: `nightlink_tests`, also inside `moving_tests`).
**Method** (linear + gravity, Herget/HelioLinC-flavoured): (1) tracklet -> attributable (ra, dec, rates, with covariance from the linear tangent-plane fit); (2) pair gate: a grid of 20 x 13 (log rho, rho_dot/v_esc) heliocentric states, Sun-bound only (E < 0, r > 0.15 AU), is propagated with universal-variable two-body motion and
light-time to the later tracklet and compared with its position and rate, using both attributable uncertainties grown over the gap plus a model floor; (3) the best 4 hypotheses of each surviving pair start a 6-parameter least-squares orbit fit (Levenberg-Marquardt; soft prior v < 0.98 v_esc) on all observations, with
sigma_eff^2 = sigma^2 + 0.3"^2 + (0.05"/day x |t - t_ref|)^2 (this floor absorbs planetary perturbations and the geocentre/parallax approximations); a pair is accepted at chi2/dof <= 4 and max residual <= 6 sigma_eff; (4) accepted pairs are grown greedily: a tracklet is added if the orbit predicts it within 6 sigma_eff (or it
formed an accepted pair with a member) and the refit still passes; one tracklet per night and station per group.  Observer = Earth-Moon barycentre from the astropy built-in ephemeris plus the station offset from the MPC parallax constants when `--obs-code` is a ground code.  `vet_with_links` marks tracklets `linked` / `pair` / `single`: a linked
tracklet has passed an orbit fit over >= 2 nights, which replaces the population rate prior for it (the rate prior in `orbitlink.py` is still what vets single-night tracklets).
**Why >= 3 tracklets.**  Two short-arc tracklets give 8 attributable numbers for 6 orbital parameters, i.e. only 2 degrees of freedom of constraint; in a dense field chance pairings pass the fit.  Measured earlier with pair links allowed (Horizons set, long baselines): 3-10 wrong groups per run; with >= 3 tracklets required there were 0 wrong groups in every test below, and 2-tracklet links are returned separately as `pairs` (candidates only; at 1" noise and nights {0,1,7,14} 9 of 13 pairs were wrong, see below).
**Validation** (numbers from this box; truth labels known in all of them):
* SYNTHETIC two-body truth (exact dynamics, so optimistic), 30 objects with p_detect 0.9 per night + 30 random single-night decoy tracklets per night, one 6x6 deg field (RA 100-106, Dec 0-6), 3 exposures 0.5 h apart per night, 2 seeds per row, near-circular bound main-belt-like states: recall of objects with >= 3 tracklets 1.000 at 0.1" and 0.3" noise for nights {0,1,3}, {0,2,6}, {0,1,7,14}, {0,3,10} (47-56 groups per row), 0.979-1.000 at 1.0"; **0 wrong groups in all 12 configurations** (24 fields); 2-tracklet pairs: 49 of 49 pure (two-night sets) at 0.1", 1 wrong of 49 at 0.3" and 3 of 51 at 1.0" for a 7-day gap.  Run time 14-42 s per 170-230 tracklets (321 s at 1.0" with a 14-day span).
* SEMI-SYNTHETIC with REAL dynamics: JPL Horizons n-body positions of 70 real numbered asteroids in one 8x8 deg window (SBDB elements picked the window, RA 190 Dec -4; MJD 60950 = 2025-10-02), Gaussian noise added, decoys = real tracklets copied from other epochs and shifted to random positions in the window (realistic rates; 70 per night), 420 tracklets: noise 0.3": nights {0,1,3}: 67 of 70 objects recovered (recall 0.957), 0 wrong groups; nights {0,2,6}: 64 of 70 (0.914), 0 wrong groups, 6 pure and 4 wrong pairs.
  The 3-6 missed objects were not analysed individually (candidate causes: model floor too small for them, coarse hypothesis grid, a tracklet lost to the one-per-night rule).  The {0,1,7,14} set and 1.0" noise were started but not finished (runtime), so there is no semi-synthetic number for baselines > 6 days; the real-MPC set below covers 5-14 days.
* REAL astrometry (MPC `get-obs` archive of 48 numbered near-Earth asteroids; tracklets = >= 2 observations from one station within 3 h; the 14-day window starting MJD 59246 with the most objects on >= 3 station-nights; real astrometric errors, parallax for each station via the MPC constants, sigma from `obs.station_sigma`; truth = designation, objects pre-selected as known and well observed):
  5 days: 27 tracklets of 11 objects, recall 1.000 (3 of 3 objects with >= 3 tracklets), 0 wrong groups; 10 days: 112 tracklets of 18 objects, recall 0.933 (14 of 15), 0 wrong, 1 pure pair; 14 days: 176 tracklets of 21 objects, recall 0.944 (17 of 18), 0 wrong, 1 pure pair.  Sky density is low (objects are spread over the sky, not in one field), so this tests model
  errors and real astrometric noise, **not** confusion; the confusion test is the semi-synthetic one above.  The network test `test_real_mpc_astrometry_links_known_asteroids` (12 of these objects, 10 days) passes (about 2 min with a cold cache).
**Limits.**  The search itself is still two-body and Sun-only: planetary perturbations and non-gravitational forces are covered only by the error floor (0.3" + 0.05"/day), so arcs of more than ~2-4 weeks, close approaches to planets and fast NEOs at small distances can fail the gate.  `--nl-nbody` is an optional second step, not part of the gate.  It re-fits each group the two-body linker already accepted.  It does not go back and recover a group the gate rejected.  The old Horizons and MPC recall numbers above were measured without this step.

**N-body refinement** (`moving/nbody_refine.py`, `--nl-nbody`).  The linker process does not import rebound, assist, or spiceypy.  It writes the group's observations and the two-body state to a child process, `python -m moving.nbody_worker`.  If ASSIST and the DE440 kernels answer `orbit.assist_available()`, the child calls `orbit.differential_correction`.  Otherwise, if `OGF_CODES_ROOT` or `~/BACKUP/3.5ST` contains `neo_orbit_calculator` with `de440s.bsp`, the child fits the same observations with that external Fortran integrator (Dormand-Prince in real128, JPL DE440s, SB441-N16, full multi-body 1PN, J2/J4/J6; area-to-mass is 0, so radiation and solar-wind terms are off).  CODES is not shipped in this repository and is not imported by the web.  The N-body residual uses the observation sigma plus a flat 0.05" floor, not the two-body 0.3" + 0.05"/day floor.  A fit is kept when chi2/dof <= 4 and the largest residual is <= 6 sigma, including a fit that missed the formal step tolerance but already met those cuts.  A converged fit that fails the cuts is removed from `groups` and listed in `nbody_rejected`.  An integrator error leaves the two-body group in place.  Tests: `moving/tests/test_nbody_refine.py` (`run_all_checks.sh`: `nbody_refine_tests`).  The CODES recovery tests skip when that tree is absent.  The ASSIST branch was not executed on this machine (assist is not installed); `moving/tests/test_orbit.py` still covers `differential_correction` when the kernels exist.

Measured here with CODES only, one massless particle, ecliptic elements a=2.7 AU, e=0.15, i=8 deg (node 40, argument 20, mean anomaly 30), geocentre, three observations per night 0.5 h apart, Gaussian noise 0.2".  Truth positions come from the same integrator.  Element errors are osculating elements at the N-body epoch (0.1 day before the first observation; MJD 60999.9 when the first night is MJD 61000), two-body state moved to that epoch with the Kepler propagator.

| nights | seed | N-body status | iterations | raw rms two-body | raw rms N-body | chi2/dof | Δa two-body | Δa N-body | Δi two-body | Δi N-body |
|---|---|---|---|---|---|---|---|---|---|---|
| 0, 7, 21 | 7 | converged | 3 | 0.1110" | 0.1125" | 0.446 | +0.002718 AU | +0.002741 AU | +0.000798 deg | +0.001000 deg |
| 0, 20, 45 | 9 | not_converged, kept | 12 | 0.1836" | 0.1837" | 1.191 | +0.000085 AU | −0.000266 AU | −0.000581 deg | −0.000103 deg |
| 0, 40, 120 | 11 | converged | 10 | 0.1371" | 0.1370" | 0.662 | −0.000383 AU | +0.000002 AU | −0.001387 deg | +0.000022 deg |

The in-arc residual does not improve on these arcs: both models adjust six parameters to nine observations, and the planetary signal sits under the 0.2" noise.  At 21 days the two semimajor-axis errors agree at 0.00002 AU.  At 45 days the two-body semimajor axis was closer.  At 120 days the N-body semimajor axis and inclination were closer (Δe was +0.000010 and +0.000005).  One particle and one noise seed per row, not a recall on the Horizons or MPC sets.  Close approaches and fast NEOs were not part of this measurement.
Tracklets need >= 2 detections and short arcs (hours); the (rho, rho_dot) grid is coarse (the LM fit from the 4 best hypotheses does the real work); hypothesis generation and gate cost O(pairs x 160 states): ~10-60 s for ~200 tracklets, it has not been run on 10^4-10^5 tracklets and would need a spatial pre-index (HelioLinC-style) for that;
`max_motion_deg_day` (2.5) rejects faster pairs; the group assignment is greedy; 2-tracklet candidates have a high false rate in dense fields; no covariance-based (Mahalanobis) acceptance, only chi2/dof and max-residual thresholds that were chosen on the synthetic and Horizons sets (the MPC set was not used to choose them; the Horizons set was used during development);
a missed night splits an object into pairs or singletons; moving-object tracklets with a wrong within-night association poison the group.  Not applied automatically by the pipeline `run` mode: the user supplies several `tracklets.json`.

## Session recorder: Moving Objects steps

The eight steps (`moving.align/difference/link/identify/orbit/transients/lightcurve/export`) are recorded and exported by
"Save Session as Python Script...".  The exposures after `--files` count as session images (also from `~/.ds9/mast_cache`),
so the exported script contains `@{IMG:key}` placeholders; `--mode replay --field NAME:main=a.fits,img2=b.fits,...` re-runs on the
original exposures, `--mode pipeline` (default) on any set of exposures supplied for the field (all images given are used).
Ephemeris files: the exported `SESSION_META['ephem_dir']` is used as `OGF_EPHEM_DIR` when the script runs under a different HOME.
Verified on the real cached HST ACS 2015 BB89 exposures (GUI run: align 34 s, difference 136 s, link 35 s, identify 21 s, orbit 37 s):
replay in an empty HOME with `env -i` ran all 8 steps (321 s) and produced files byte-identical to the GUI run (detections.tsv,
movers.tsv, tracklets.json, identified.tsv, orbit_tracklet0.json/.kv, transients.tsv, lightcurves.json, elements_export.json,
tracklet0.obs80).  `scripts/verify_moving_session.sh` repeats that comparison (needs network and the cached data; exits 77 when
they are missing).  Pipeline mode on the same four files in reversed order ran all 8 steps in 316 s, but the outputs differ
(9824 vs 9809 detections; the template/anchor exposure changes), as expected.  Caveat: the recorded `moving.orbit --tracklet 0`
is "whichever tracklet is rank 0" - on new data that is not necessarily a known object.

## Limitations (summary)

* Cross-night orbit-fit linking still searches with a two-body model.  `--nl-nbody` re-fits accepted groups in a separate process (ASSIST if installed, otherwise the external CODES integrator).  The CODES measurement in the nightlink section is one main-belt-like particle, not a new recall on the Horizons or MPC sets.
* LSST provider unavailable without an RSP login; test field outside DP1. No real DP1 FITS tested.
* No absolute astrometric tie for the test field; Gaussian PSF; dense archive CR flag; recovery limited to mag ≲ 23 in the first injection test (new linker: see the benchmark, ceiling set by detection);
  the true BB89 tracklet is not the best-scored candidate.
* Tracklets from one HST orbit cannot give orbits.  Arc-aware classification (`kepler.classify_arc`, `class_probs_from_covariance`):
  statuses `classified` / `ambiguous` / `unclassified`; arcs < 1 d, degenerate fits (e >= 1, a <= 0, not determined) and class
  probabilities below 0.70 are *not* given a class label; the uncertainty (Monte-Carlo from the covariance, or ranging samples for
  short arcs) is reported as `orbit_class_probs`, `class_note`.  Tests: `moving/tests/test_classify_arc.py`.
* Station sigmas are defaults; Carpino rejection is aggressive on short arcs.
* Transient / light-curve / export modes were exercised on one field only; light curves have no PSF photometry; host association needs a user catalogue.

## References

Zackay, Ofek & Gal-Yam 2016, ApJ 830, 27 (ZOGY); Holman et al. 2018, AJ 156, 135 and Heinze et al. 2022, PSJ 3, 28 (HelioLinC);
Holman et al. 2023 (ASSIST, REBOUND); Milani et al. 2004; Virtanen et al. 2001; Granvik et al. 2009, Icarus 202, 447;
Muinonen et al. 2010; Carpino et al. 2003, Icarus 166, 248; Veres et al. 2017, Icarus 296, 139; Eggl et al. 2020, Icarus 339, 113417;
Bowell et al. 1989, in Asteroids II.
