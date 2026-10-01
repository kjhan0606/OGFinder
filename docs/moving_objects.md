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
| Difference | ZOGY (Zackay, Ofek & Gal-Yam 2016) in Fourier space; template = median of the other aligned exposures. `zogy.zogy` also returns `alpha_new`/`sigma_alpha_new` = flux in the *target exposure's* flux scale (alpha * Fn / Fr), which `difference_chip` uses | sky-noise-only normalisation (no source-noise term, no astrometric-error terms); one PSF per chip (empirical from field stars, Gaussian fallback with FWHM floor 0.085 arcsec), **no spatial PSF variation**; template PSF = target PSF. A PSF 20% too narrow gives ~22% low flux (`moving/tests/test_zogy.py`) |
| Detect | thresholding on the S/N image, trail (elongated) channel, classes `point trail artefact_cr artefact_edge artefact_static artefact_dipole negative faint` | the archive CR mask is dense and mislabels real detections as `artefact_cr`; `artefact_static` = non-trail detection with template S/N > 8 |
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
* `ds9/library/ds9_moving.py` CLI (`--mode setup|fetch|align|difference|link|identify|orbit|transients|lightcurve|export`);
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
    no baseline for that (`helio_linc` needs >= 2 nights and is not wired into the linker).
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
