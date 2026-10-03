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

## Depth and completeness maps (steps *Depth Maps* and *Completeness Map by Region*; `plugins/completeness/depthmap.py`, library `ogfkit/depth.py`)

**Depth (per image and per position).** The pixel rms varies over an image (exposure, dither pattern, scattered light) and the noise of an aperture sum is not sqrt(N) sigma
for correlated (drizzled / resampled) noise.  The step measures

* the local pixel rms map (sep mesh with sources masked, or a user rms / weight map rescaled to the measured noise: `--rms-file`, `--weight-file`),
* the factor `f_ap = sigma_N(r) / rms` from blank apertures of radius `--aper-radius` (each aperture sum divided by the rms at its centre, 3000 apertures, clipped std; the
  fit `f = alpha N^beta` through 0.5-3 r is reported as `law`), and likewise `f_box` from blank `--box-arcsec` boxes (a power law through 0.2-1.0 x box when there are fewer than
  300 independent blank boxes - `factor_direct` and `factor_law` are both in the summary),

and builds `depth_mag_map.fits` (`zp - 2.5 log10(nsigma rms f_ap) - apcorr`: the `--nsigma` point-source limit in the aperture, with the aperture correction `--apcorr`),
`depth_sb_map.fits` (`zp - 2.5 log10(box-nsigma rms f_box / area_arcsec2)`, the usual "3 sigma in 10 x 10 arcsec" surface-brightness limit; needs the pixel scale), `depth_rms_map.fits`,
`depth_tiles.tsv` (per tile of `--tile` px: blank-aperture scatter measured in the tile versus the model, ratio ~ 1 +- 1/sqrt(2 n)), `depth_area.tsv` (fraction of the valid area deeper than a magnitude),
`depth_summary.json` (median / 5-95 % / min / max limits, factors, tile check) and `depth_plot.png`.  Catalog columns: `DEPTH_RMS`, `DEPTH_LIM`, `DEPTH_SBLIM` at each object's position.
Pixels with a non-finite value or an rms below 5 % of the median (empty margins) are invalid (NaN).

**Completeness by region.** The existing injection / recovery (`run_completeness`, now with `return_records`) runs once over the whole image with `--per-bin` injections per magnitude bin
(use a few hundred per region wanted); the recovered fractions are binned by region (`--regions grid --grid 3x3`, `rms` = classes of equal area of the depth map, `file` = label image) giving
`lim50` / `lim90` per region (`compmap_regions.tsv`), and all injections together give a **universal curve** in `x = m - depth_map(x, y)` (`compmap_universal.tsv`, logistic fit; its
50 % / 90 % points are `delta50` / `delta90`).  Hence `compmap_lim50_map.fits` / `compmap_lim90_map.fits` = depth map + delta (limit at each pixel) and
`compmap_frac_mXX.fits` (`--maps-at 26,27`: completeness at magnitude XX at each pixel).  The summary reports the collapse: standard deviation of `lim50 - depth` over the regions (small = depth
map explains the regional differences) and the regression slope of lim50 on depth (1 expected).  Catalog columns: `COMPL_REGION`, `COMPL_LOC` (completeness at the object's MAG_AUTO), `COMPL_LIM50_LOC`, `COMPL_LIM90_LOC`.
The `sep` detector gets `local_rms` (threshold relative to the local rms map, as SExtractor does) in this mode (`--global-threshold` for the old behaviour); with a global threshold the limits do *not* follow the noise.

CLI: `python plugins/completeness/depthmap.py IMAGE --work DIR --mode depth|compmap [--catalog CAT] --mag-zeropoint ZP --pixel-scale S --aper-radius 3 --nsigma 5 ...`.

### Validation (`plugins/completeness/validation/depth_validate.py`, `depth_*_report.json`; tests `tests/test_depth.py`, 9 + the 9 above)
* **Truth = brute force** (20000 blank apertures, clipped std, no mask) on noise-only images; depth map median minus truth: white noise (sigma 0.5 and 2) +0.002 mag; correlated noise (Gaussian 1.0 / 1.5 px, sigma_pix fixed; the aperture factor grows from 5.1 to 14.7 / 19.3) +0.038 / +0.069 mag;
  white-noise SB limit (50 px boxes) -0.043 mag.  Spatially varying noise (gradient 1 -> 2.5, four patches 1 / 1.5 / 1.4 / 2.1) with 700 injected sources and a NaN region, 8 blocks of 400^2 each: mean +0.046 / +0.057 mag, scatter 0.02, worst block 0.07 / 0.09 mag
  (map is ~0.05 mag too deep: the model's source mask removes a little of the positive noise); with an exact weight map as input +0.06 (gradient) and +0.19 +- 0.11 (patches; blocks that straddle a patch edge, where the 64 px mesh rms and the exact weights differ).
* **Regional completeness, truth = stand-alone images** (four noise levels sigma 1 / 1.4 / 2 / 2.8 as four bands of one 1000^2 image, sep local threshold, ~1000 injections per band): lim50 minus the lim50 measured on an image of that noise alone:
  -0.10 / -0.07 / -0.01 / +0.03 mag; lim90 +0.03 / -0.26 / -0.24 / +0.26 (the 90 % point is poorly determined with ~80 injections per bin); lim50 - depth = 0.77 / 0.79 / 0.78 / 0.72, std 0.03 mag across a 1.05 mag depth range.
* **HUDF F160W** (1800^2 crop, zp 25.94, 0.06"/px): 5 sigma limit in r = 5.8 px (0.35") median 28.80 [28.67, 28.85] (the rms class regions differ by only 0.1 mag - the mosaic is nearly uniform here), 3 sigma 10" x 10" limit 30.22 (power-law factor; no blank 10" box exists);
  tile check ratio median 0.91 (80 tiles; one edge tile at 8.3); injection (3500 stars, sep local threshold): lim50 30.15, lim90 26.6 (poorly constrained), lim50 - depth = 1.36 +- 0.14 over four depth classes; by 3 x 3 grid the lim50 of the cells ranges 29.76-30.33 although the depth map is flat to 0.1 mag - crowding and blending decide the regional completeness there, not the noise.
* **M51** (600^2, zp 25, 1.8"/px): 5 sigma limit 15.43, 3 sigma 2' x 2' limit 21.6 mag/arcsec^2, lim50 14.96 (galaxy light and crowding dominate: lim50 - depth = -0.69 +- 0.35, no collapse - the maps document the noise only).
* GUI (M51, session replay identical): `verify_newplugins.tcl` section `depth`, 29 checks.

Limitations: the depth is the *noise* depth of an isolated source; crowding, blending and extended-source surface brightness are measured only by the injection (regional lim50 beyond the noise depth); the correlated-noise factor is global (a patchy noise correlation, e.g. mosaics of different drizzle kernels, needs separate regions); the surface-brightness limit of boxes larger than ~1/5 of the image is extrapolated; rms/weight maps are rescaled by a single median ratio; no PSF-dependent optimal-aperture depth; injections are not given Poisson noise.

## Limitations

* The catalog/metadata link is by row count only.  Injected sources get no extra photon noise unless `--gain` is given.  The recovered fraction is for sources that are *isolated* from
  already-detected objects unless `--allow-blends`; crowded-field completeness should use `--allow-blends` (and the crowded-photometry artificial-star test).
* Galaxy magnitudes of the sextract detector are biased (MAG_AUTO) - the bias column documents it; the threshold is not a surface-brightness completeness.
* The false-positive estimate comes from the negative image (not a genuine "no-source" image) and counts no spurious deblended pieces of real sources.
