# Moving Objects: asteroids, transients and reference images in OGFinder

*Moving Objects* is a menu (next to *Analysis* in the catalogue panel) that runs a multi-epoch pipeline on HST/ACS-type
`*_flc.fits` exposures: fetch references, align, difference, detect, link tracklets, identify known objects, fit an orbit,
find static transients, make light curves and export. All numerics are in the `moving/` package and are driven through one
CLI, `ds9/library/ds9_moving.py`; `ds9/library/ogf_moving.tcl` is only menu, dialogs, table, markers and the orbit window.

**What existed before.** There was no asteroid mode in this repository. It had stock DS9 SkyBoT access (`catskybot.tcl`) and a
wish-list section about difference imaging in `docs/manual/main.tex`. Everything described here is new.

## Menu (exact labels)

`Moving Objects` → `Fetch Reference...`, `Align`, `Difference`, `Detect`, `Link Tracklets`, `Identify Known Objects`,
`Orbit Fit...`, `Transient Candidates...`, `Light Curve`, `Export...`, (separator) `Select Exposures...`,
`Show Results Table`, `Work Directory...`.

*Fetch Reference...* has a provider radio button **MAST / LSST / both**; the choice is written to the CLI as
`--provider mast|lsst|both` and therefore appears in the session-recorder step arguments.

## Methods: implemented vs simplified

| step | implemented | simplification / caveat |
|---|---|---|
| Fetch | MAST (astroquery, public), Pan-STARRS cutouts, SkyBoT, Horizons | see LSST table below |
| Align | Gaia DR3 (TAP) matching when stars exist, else relative chain with a polynomial/shift+rotation model between exposures | in the test field no Gaia stars matched; result is relative only (4–19 mas rms), the anchor exposure has no absolute tie |
| Difference | ZOGY (Zackay, Ofek & Gal-Yam 2016) in Fourier space; template = median of the other aligned exposures | sky-noise-only normalisation; template PSF approximated by the target PSF; PSF is Gaussian with FWHM = max(measured, 0.085 arcsec) because the empirical pooled PSF was unusable here |
| Detect | thresholding on the S/N image, trail (elongated) channel, classes `point trail artefact_cr artefact_edge artefact_static artefact_dipole negative faint` | the archive CR mask is dense and mislabels real detections as `artefact_cr`; `artefact_static` = non-trail detection with template S/N > 8 |
| Link | observer-parallax hypothesis grid over 1/Δ, pair + intermediate search with a KD-tree, trail PA/length/flux consistency in the score (HelioLinC-style, Holman et al. 2018 / Heinze et al. 2022; **simplified, no heliocentric clustering**) | see recovery numbers |
| Identify | SkyBoT candidates re-checked with Horizons using the true HST observer (`500@-48`) | |
| Orbit | statistical ranging over the admissible region (Granvik et al. 2009 / Muinonen et al. 2010 style, simplified), 2-body polish, then ASSIST/REBOUND differential correction (DE440 + 16 massive asteroids + GR) with Carpino et al. (2003) outlier rejection, Eggl et al. (2020) star-catalogue debiasing, default station sigmas | no classical Gauss IOD (not implemented); station weights are defaults, not the Veres et al. (2017) table; solar light deflection and stellar aberration neglected; `mcmc_posterior` (emcee) exists but is untested; arcs shorter than ~1 day are flagged `determined=False`: elements are blanked, only class probabilities and ranging quantiles are shown |
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
  prints `#MOVING {json}` status lines; `ds9/library/ogf_moving.tcl` GUI
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
  Hilda 0.17) while the true object is a main-belt asteroid (a = 3.17 AU, e = 0.106, i = 19.0 deg); this classification is **wrong**
  and must not be used for such short arcs.
* Injection–recovery (30 synthetic Gaussian-PSF movers in real chips, mag 21–26.5, 1.5–40 arcsec/h geocentric): 30/30 detected in ≥3
  exposures as single detections, but only 8/30 linked correctly (all with mag 21.4–22.8; none at mag ≥ 23.0). Fitted rates match the
  injected ones (e.g. 5.3/5.3, 11.7/11.7, 8.3/6.8 arcsec/h; one object at 32.8 arcsec/h was fitted as 74.9). The 13 correct tracklets (covering the 8 linked objects; `rank4.py`) ranked at positions 16–284 (list order, sorted by n then score) of the 400 returned tracklets (score 6.7–11.8); most returned tracklets are false. Injected sources use a Gaussian PSF, not a real ACS PSF.
* Alignment: relative rms 19.4 / 4.4 / 3.7 mas (exposures 2–4 w.r.t. the anchor), no Gaia stars matched.
* Orbit fit pipeline check on an asteroid with a long observed arc (25153, public MPC astrometry, ASSIST force model, Eggl debiasing):
  * 90-day arc, 277 obs (100 used), rms 0.51 arcsec, χ²_red 1.15. Elements vs JPL SBDB, (fit − JPL)/σ_fit: a −1.4, e −0.9, i +0.5,
    Ω −0.6, ω −1.5, M +1.6.
  * 4-day arc, 39 thinned obs (13 used after rejection, 26 rejected), rms 0.44 arcsec, χ²_red 1.15. Deviations +1.2 σ (a), +1.2 (e),
    −1.5 (i), +1.4 (Ω), +1.7 (ω), −1.4 (M); acceptable but at the edge.
  * Same 4-day window through the shipped CLI (`--mode orbit --designation 25153 --mjd-min 58495 --mjd-max 58499`, all 115 observations, 37 used, no thinning): rms 0.56 arcsec, χ²_red 1.46; deviations +1.2 σ (a), +1.3 (e), −2.3 (i), +1.7 (Ω), +1.9 (ω), +2.4 (M), i.e. 2σ-level offsets with 78 observations rejected; this short arc is marginal.
  * A bug was found and fixed in the ASSIST propagator wrapper (one simulation reused for forward and backward sweeps): before the fix
    a differential correction started from the JPL state gave 57 arcsec rms, after the fix 1.76 arcsec (χ²_red 0.51) with default weights.

## Limitations (summary)

* LSST provider unavailable without an RSP login; test field outside DP1. No real DP1 FITS tested.
* No absolute astrometric tie for the test field; Gaussian PSF; dense archive CR flag; recovery limited to mag ≲ 23 in the injection test;
  the true BB89 tracklet is not the best-scored candidate.
* Tracklets from one HST orbit cannot give orbits; class probabilities from such arcs are unreliable.
* Station sigmas are defaults; Carpino rejection is aggressive on short arcs.
* Transient / light-curve / export modes were exercised on one field only; light curves have no PSF photometry; host association needs a user catalogue.

## References

Zackay, Ofek & Gal-Yam 2016, ApJ 830, 27 (ZOGY); Holman et al. 2018, AJ 156, 135 and Heinze et al. 2022, PSJ 3, 28 (HelioLinC);
Holman et al. 2023 (ASSIST, REBOUND); Milani et al. 2004; Virtanen et al. 2001; Granvik et al. 2009, Icarus 202, 447;
Muinonen et al. 2010; Carpino et al. 2003, Icarus 166, 248; Veres et al. 2017, Icarus 296, 139; Eggl et al. 2020, Icarus 339, 113417;
Bowell et al. 1989, in Asteroids II.
