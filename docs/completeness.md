# Completeness and limiting magnitude (plugin `completeness`)

Artificial-source injection tests for the *existing extraction* (or any other detection step).  Measure tab, chip **Complete.**
(`[Complete. v] [>] [gear]`).  Backend `plugins/completeness/completeness.py` (pure functions, JSON in/out, usable on a server).

## Method

1. A work image is prepared (optionally a `--crop-size` square around the centre; the shared mask is respected: sources are never put on masked pixels).
2. `n_bins` magnitude bins between `--mag-min` and `--mag-max`; per bin `--per-bin` stars (PSF stamp from `--psf` or a Gaussian of `--psf-fwhm`, sub-pixel shifted)
   or PSF-convolved Sersic galaxies (n in `--n-sersic`, axis ratio, PA, half-light radius drawn from `--re-bins`) are added at random positions, `--per-image` per realisation,
   with a minimum separation.  Positions on baseline detections are avoided unless `--allow-blends` (so recovery is *not* reduced by blending with real sources by construction).
3. The detection callback runs on each realisation; a recovered source is the nearest detection within `--match-radius` px (greedy one-to-one).
4. Per bin: recovery fraction with Wilson 68 % interval, photometric bias (`median(m_det - m_in)`) and scatter (1.4826 MAD), and (galaxies) the same per size bin (`by_size`).
5. Limits: logistic `f(m) = fmax / (1 + exp((m - mh) / w))` fit (`fit.m50`, `fit.m90`); `lim50`, `lim90` are the direct interpolation of the measured fractions
   (NaN/`null` when not reached in the magnitude range - no extrapolation).
6. False positives: the detector is also run on the *negated* image (sky noise only, real sources become troughs); detections per 1e6 pixels and per arcmin^2 are reported,
   plus unmatched detections per injected image.  This is only an estimate: it assumes symmetric noise.

## Detection callback API

    detect(image: float32 2-D array, mask: bool array | None) -> dict-of-arrays | list-of-dicts   # keys x, y (0-based), optional mag, a, b
    run_completeness(image, detect, kind='star', mag_range=(22, 29), n_bins=14, per_bin=60, zp=25.0, psf=None, ...) -> dict

Built-in detectors: `sextract` (`bin/ds9_sextract`, the application's extraction, with `--detect-thresh/--detect-minarea` taken from the extract parameters),
`sep` (options `thresh`, `minarea`, `filter`), or `package.module:factory` (factory called with the detector options as keywords; picklable for `--n-workers > 1`).

## Catalog metadata

The result is written into the catalog metadata so downstream code can find the depth: `~/.ds9/catalog_meta.json` (`ogfkit/meta.py`; valid only while the catalog row
count matches), saved as `<name>.meta.json` next to *Save Catalog* output (also in the exported session script), and copied into the FITS-table header as
`HIERARCH OGF COMPLETENESS_LIM50` etc. by `ds9_fits_export.py`.  Catalog columns: `COMPL_FRAC` (completeness at the object's MAG_AUTO from the fitted curve), `COMPL_LIM50`, `COMPL_LIM90`.
Cataloged objects with `COMPL_FRAC < 0.5` are therefore flagged by value rather than removed.

## CLI

    python plugins/completeness/completeness.py IMAGE [--detector sextract|sep|pkg.mod:fn] [--kind star|galaxy] [--mag-min --mag-max --n-bins --per-bin --per-image]
        [--mag-zeropoint 25] [--psf FILE | --psf-fwhm 3] [--match-radius 3] [--re-bins 2,4,8,16] [--catalog CAT] [--mask M] [--crop-size N]
        [--json-out --curve-out --plot-out --meta-out] [--n-workers N] [--seed 1]

## Validation (measured; `plugins/completeness/tests`, 9 tests, ~6 s)

* **Analytic check.** SEP with filter off, minarea 1, 4 sigma threshold: a star is detected iff any pixel exceeds 4 sigma.  The expected recovery averaged over sub-pixel positions was
  computed analytically (Gaussian PSF, Gaussian noise): all 10 bins |z| < 3.5, chi2 < 35, lim50 within 0.08 mag and lim90 within 0.12 mag of the analytic limits.
* False positives on pure noise match N(1 - Phi(4)) within 4 sqrt(N).  Injected flux conserved (< 1 % n=1, < 5 % n=4 compact).  Mask respected.  Logistic fit recovers (26.0, 0.4, 0.97) within 0.03/0.05/0.01.
* **HUDF F160W 800x800 crop** (zp 25.94, FWHM 3 px): stars sextract lim50 30.45 / lim90 29.48, sep 30.59 / 29.63; galaxies sextract lim50 29.07 / lim90 27.35
  (Re 2-4 px 29.66 / 28.71; 4-8 px 29.24 / 27.77; 8-16 px 28.14 / not reached); negative-image false positives about 330 per 1e6 px at 1.5 sigma / minarea 5;
  star bias ~0 mag, scatter 0.16 mag (sextract) vs 0.02 (sep) at 26.8; galaxy MAG_AUTO bias +0.5...+0.9 mag at the faint end.
* GUI on M51 (zp 25, mags 12.5-17.5, 60 stars): 25/60 recovered, lim50 = 15.0, lim90 not reached (null); replay writes the identical catalog.

## Limitations

* The catalog/metadata link is by row count only.  Injected sources get no extra photon noise unless `--gain` is given.  The recovered fraction is for sources that are *isolated* from
  already-detected objects unless `--allow-blends`; crowded-field completeness should use `--allow-blends` (and the crowded-photometry artificial-star test).
* Galaxy magnitudes of the sextract detector are biased (MAG_AUTO) - the bias column documents it; the threshold is not a surface-brightness completeness.
* The false-positive estimate comes from the negative image (not a genuine "no-source" image) and counts no spurious deblended pieces of real sources.
