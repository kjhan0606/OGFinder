# Satellite / aircraft trail removal (`plugins/trails`, `ogfkit/trails.py`)

One click in the ds9 GUI: **Remove Trails** (Detect tab, chip *Remove Trails*, also in the plugin menu) detects straight satellite and aircraft
trails, builds a mask whose width comes from the measured cross-profile plus a margin, writes it into the **shared mask manager** (the same flag mask
that Auto Mask / ICL / LSBG use; trail bit 32, `N_TRAIL` in the mask statistics), turns the mask overlay on, flags the catalogue objects that touch a
trail (`TRAIL_FLAG`, `TRAIL_ID`, `TRAIL_DIST` columns in the catalogue panel) and opens a small result panel with **Undo**, *Settings* and *Close*.
Undo is the mask manager's own snapshot undo (`ds9_mask.py --mode undo`); the catalogue columns are reset together with it.  The image file is never modified.

## Menu / steps (all declared in `plugin.json`, so they appear in the recorder, session export/replay and the headless batch runner)
| step | what it does |
|---|---|
| `remove` (primary button) | detect + mask + overlay + catalogue flags (`trails.py IMAGE --work DIR --mask MASK [--catalog TSV]`) |
| `preview` | detect only: line regions and a text report, nothing is written to the mask |
| `undo` / `clear` | mask-manager undo; remove the trail bit entirely (`--mode trails-clear`) |
| `stack` | dialog: combine registered frames (median / mean / sigma-clip) excluding each frame's trail (`--task stack`) |

## Fill options (parameter `fill`)
* `mask` (default, right for photometry): trail pixels are only masked (bit 32); nothing is changed.  Aperture photometry that honours the mask excludes them.
* `interpolate`: `<base>_trailfree.fits` in which the band is replaced by a linear interpolation across the trail between the two flank medians
  (`fill-noise` adds Gaussian noise of the sky sigma so that the noise level is preserved).  For display / morphology; **not** for photometry of objects on the trail.
* `stack`: `<base>_trails_mef.fits` with `SCI` + `TRAILMASK` extensions (and `trails_mask.fits`) for stackers that read mask extensions, plus the `stack` step
  that combines frames with per-frame trail masks.

## Parameters (declarative, `plugin.json`; CLI flag = `--name`)
threshold 8 (robust z of the best run against trail-free parallel lines), min-length 0 (= 12 % of the image diagonal), max-trails 8, smooth 1 px,
min-aspect 12 (length / FWHM), max-fwhm 40 px, catalog-check (accept z >= 0.7 x threshold when >= 2 elongated aligned catalogue objects lie on the line),
margin 2 px and edge-sigma 2 (mask half-width = box half-width + edge-sigma x blur + margin), end-extend 3 px, fill, fill-noise, flag-catalog, touch-k 2.5
(footprint ellipse scale for the touch test), show-overlay, confirm-panel, show-filled.  CLI only: `--keep-bleeds`, `--pixel-scale/--mag-zeropoint`
(adds peak surface brightness and width in arcsec to `trails.json`), `--spike-fac` (default 30, see step 6), `--curved` (also a manifest checkbox "Also look for curved / segmented trails"), `--register {none,wcs,auto}` for the stack step (the stack dialog has a "Register by WCS" field).

## Algorithm
1. *Standardise*: background mesh (32 px) subtraction, local high-pass-MAD noise map (floored at 0.7 x the global value so that galaxies are down-weighted),
   Gaussian smoothing, clipping at +-3 sigma (a bright galaxy cannot dominate a line integral).
2. *Candidates*: validity-weighted Radon transform of the clipped map binned to <= 768 px (0.5 deg steps), per-angle MAD-normalised, peaks > 4.
3. *Verification at full resolution*: for each candidate, neighbouring (theta, rho) are scanned with strips of width 3/6/12/24 px; the statistic is the
   contrast of the strip against its two flanking strips, the best contiguous run along the line >= min-length is kept, and its z-score is measured against
   100 random null lines of the same image processed identically.  Partial and faint trails therefore need no a-priori length; crossing galaxies only
   interrupt a run when they are inconsistent with a constant-amplitude line (the contrast is clipped).
4. *Greedy acceptance with deflation*: after a trail is accepted its band is blanked and the remaining candidates are re-verified (no ghost trails along the wings).
5. *Profile*: median cross-profile along the run fitted with box (x) Gaussian -> centre, width, blur, amplitude.  Mask half-width = box half-width + 2 blur + margin.
6. *Rejections*: broad short ridges (galaxy chains: length < 12 x FWHM or FWHM > 40 px) and **linear features of very bright stars** (CCD bleed columns near the pixel
   axes: core > 300 x the run amplitude; diffraction spikes at any angle: run through a core > `--spike-fac` (30) x the run amplitude whose brightness falls off along the line by a
   factor > 2.5).  30 / 2.5 replaced 300 / 3.0 after the SDSS held-out check (below) showed that 53 of 56 straight detections in 160 SDSS frames were saturated-star spikes; the price on HUDF injections is
   a lower recall (4 sigma: 35/40 -> 32/40, 2 sigma: 26/40 -> 22/40) because trails that cross a bright core are sometimes read as spikes; `--spike-fac 300` restores most of the old recall.
7. *Curved / segmented trails* (`--curved`, `ogfkit/trails_ext.py:detect_curved`): the straight finder is run tile-wise (tile = half of the smaller image side, >= 200 px, 50 % overlap) with a low
   threshold (z 6) and minimum segment length 0.45 tile; segments are linked into a chain when the gap between their ends is <= 0.6 tile and the direction changes by <= 30 deg per link; a chain (2-8 segments) is accepted when
   sum(z)/sqrt(n) >= 8.  Accepted segments are added as trails with a `group` number (same group = one curved trail) and `curved` = total turn >= 3 deg; the mask is the union of the segment masks (piecewise linear).  Straight trails
   already found are not duplicated.  *Flicker*: `trails.json` gets `duty` (fraction of 12 px bins along the track above half the median amplitude of the brightest half), `n_gaps`, `cv` and a `flicker` flag (duty < 0.8) from `along_profile`;
   the straight finder tolerates gaps through its best-contiguous-run statistic, there is no extra gap-bridging rule.
8. *Contaminating flux* (`trail_flux`, catalogue columns `TRAIL_FLUX`, `TRAIL_FRAC`): the fitted box x Gaussian cross-profile (centre, width, blur, amplitude) is evaluated over the object's 2.5A x 2.5B ellipse
   (from A_IMAGE/B_IMAGE/THETA_IMAGE, scale k = 2.5) and summed = the light the trail adds to the object's aperture (image flux units); `TRAIL_FRAC` = that / `FLUX_AUTO` (or `FLUX_ISO`) when present.  Objects whose ellipse does not reach the trail get 0.
9. *Stack registration* (`--task stack --register auto|wcs`): a frame whose WCS differs from the first frame's (or whose shape differs) is resampled onto the first frame's grid with `wcs_resample` (bilinear), and its trail mask
   with `register_masks` (bilinear mask > 0.05, grown by 1 px; pixels without coverage are excluded from the stack); `auto` leaves already registered frames untouched, `wcs` forces it, `none` = old behaviour.
   The output line reports `registered=N`.
10. *Catalogue flags*: bit 1 = footprint (touch-k x A x B ellipse) reaches the mask band, bit 2 = centre inside the band, bit 4 = elongated (A/B >= 3) and aligned
   with the trail (a fragment of the trail itself); `TRAIL_DIST` = distance of the centre to the band edge (negative inside).

## Validation (code rev 58f69ec69 for the straight-trail tables - the detection-rate numbers below are from before the spike-rule change, see "Held-out" for the effect; `plugins/trails/validation/`, tables in `trails_validation.md`, raw numbers reproducible with `trails_validate.py`)
Synthetic trails (box FWHM 3/6/12/24 px convolved with a 1.2 px Gaussian, random angle/offset, half full-chord and half partial = 40-60 % of the chord)
are added to **real pixels**: 40 random 700x700 windows of the HUDF12 F160W mosaic (sky sigma 4.2e-4 e/s/px, 0.06"/px, AB zp 25.94) and the 600x600 M51 image.
Amplitude = peak per pixel in units of the pixel sigma measured by the pipeline; HUDF mag/arcsec^2 uses the peak per pixel.  Default parameters, 10 trials per cell.
* **Detection rate vs surface brightness** (all widths, full+partial): HUDF 0/200 at peak <= 0.45 sigma (29.1-30.8 mag/arcsec^2), 2/40 at 0.7 sigma (28.6), 6/40 at 1 sigma (28.3),
  26/40 at 2 sigma (27.5), 35/40 at 4 sigma (26.8).  M51: 0/200 at <= 0.45 sigma, 1/40 at 0.7, 9/40 at 1, 22/40 at 2, 25/40 at 4.  Full-length trails at >= 2 sigma: HUDF 38/40, M51 36/40; partial: 23/40 and 11/40.
  The ~1-2 sigma/pixel limit is far above the white-noise limit (0.15 sigma was found on synthetic white noise) because the drizzled HUDF noise is correlated and faint galaxies fill the strips;
  M51 trails inside the bright disc lose because the local noise there is the galaxy's own structure.  Misses at high amplitude: partial wide trails rejected by the aspect rule (length < 12 x FWHM) and partial trails inside the M51 disc.
* **Trails through bright galaxies** (24 lines through the brightest HUDF object per window / 24 through the M51 core, w = 6 px): HUDF 0/24, 0/24, 9/24 at 0.3, 0.6, 1.2 sigma; M51 0/24, 2/24, 21/24.
* **False trails**: 0 detections in 88 trail-free images (40 HUDF windows - overlapping, so not independent -, 8 flipped/transposed M51 variants - one image -, 40 synthetic galaxy fields incl. correlated noise).
  Best z of the trail-free images: median 3.5 / 6.2 / 5.1 (HUDF / M51 / synthetic), 95th percentile 7.1 / 10.0 / 7.0, maximum 15.3 / 10.2 / 7.8 (the 15.3 and 11.0 are diffraction spikes of a bright HUDF star, now rejected by the spike rule;
  that rule was designed after seeing them, so the HUDF number is not a blind test).  Fraction of images with best z above t: t=7 8 %, t=8 4.5 % (those are the rejected shapes/spikes).
  Injected images additionally produced 11 extra (unmatched) detections in 720 runs (9 of 11 for trails 12-24 px wide).  Upper 95 % limit on the false-trail rate per image: ~3.4 % (0/88).
* **Accuracy** (126 matched detections; peak >= 0.5 sigma for the per-width numbers): angle error median -0.003 deg, RMS 0.16 / 0.18 / 0.58 / 0.52 deg for w = 3 / 6 / 12 / 24 px; offset RMS 0.35 / 0.42 / 1.6 / 1.0 px;
  fitted FWHM / true FWHM median 1.01 (16-84 %: 0.95-1.22; 1.25 for w = 3 because of the 1.2 px blur); along-track extent IoU 0.99.  The mask covers 100 % of the pixels where the trail exceeds 0.25 sigma; mask area is 1.45x that visible area.
* **Photometry of objects near a trail** (HUDF, 40 windows x amplitude 0.3/1/3 sigma, w = 6 px through a bright object, r = 5 px aperture, local sky annulus 8-14 px, S/N > 8; bias in units of the nominal aperture noise; median / RMS).
  At 3 sigma (trail found 39/40): *unmasked*: centre on the trail +8.2 / 9.7, aperture touching the band -4.4 / 5.4 (sky annulus contaminated), clear of it 0.00 / 0.77.
  *Masked* (unmasked pixels rescaled, annulus without masked pixels): touching 2.2 / 13.9 and clear 0.00 / 2.65, identical to the same treatment applied to the trail-free image (2.2 / 13.8, 0.00 / 2.64): the trail itself no longer biases the flux, the remaining scatter is the cost of the missing aperture/annulus pixels; on-trail objects lose > 60 % of the aperture (175/207 not measurable).
  *Interpolated* (`fill=interpolate`): on-trail objects are erased (-49 / 361, their light is replaced by the flank values), touching -1.2 / 23, clear 0.00 / 7.4, again equal to the trail-free control: interpolation is for display, not for photometry.
  At 1 sigma (found 22/40) the unmasked bias is +2.6 / 3.2 (on), -1.9 / 2.0 (touching).  At 0.3 sigma the trail is not detected and not harmful (<1 sigma_ap).
* **5-frame stacks** (HUDF windows, frame noise 4x the pixel noise, a different trail per frame, trail amplitude 3 frame-sigma; trails found 117/120; objects within 5 px of a trail; bias vs the stack of the same frames without trails, in stack aperture sigma, median / RMS):
  plain mean +0.93 / 4.0, plain median +0.68 / 2.5, plain sigma-clip +0.63 / 2.8; with the trail masks excluded: mean -0.04 / 1.0, median +0.02 / 0.96, sigma-clip -0.06 / 0.98.
  Mean residual on the trail pixels: 0.61 / 0.43 / 0.46 sigma -> 0.03 / 0.02 / 0.02 sigma.  At 1 frame-sigma: +0.23 / 1.30 -> +0.01 / 0.77 (mean).
* **Real trails: HST ACS/WFC exposure `jc8m32j5q_flc.fits`** (MACS J0717, F606W, 1207 s; the acstools `findsat_mrt` example; `real_acs.py`, DQ flags masked): chip 1: a narrow bright trail (FWHM 2.6 px, peak 332 e-, 2014 px, z 94) and a wide faint one
  (FWHM 21 px, peak 8.5 e- = 0.54 pixel sigma, 3468 px, z 18); chip 2: the continuation of the narrow trail (FWHM 2.5, 2603 px, z 96).  Three near-vertical candidates (z 10-21) are the bleed/diffraction spikes of saturated stars and were rejected (before the bleed/spike rules they were reported as trails).
  Against acstools 3.8.2 MRT masks (chip 1): our masks lie 98 % (narrow) / 100 % (wide) inside the MRT masks, IoU 0.48 / 0.78 (MRT masks are about twice as wide); residual excess in the 4 px just outside our mask: +0.002 / -0.03 of the trail peak (0.04 pixel sigma for the narrow trail).
  There is no labelled truth for this exposure, so no recall/precision number.

### Extensions and held-out checks (rev 24bd44d51; tables in `trails_validation.md`, per-detection list `heldout_sdss_prefix.md`)
* **Curved arcs** (SYNTHETIC arcs, w = 6 px, injected into 700x700 HUDF windows, 8 per cell, coverage = fraction of the visible trail pixels inside the mask): straight trails (curvature 0): median coverage 0.99-1.00.  Curvature 4 deg/100 px (heading change ~28 deg over 700 px):
  median coverage of the straight finder alone 0.00 / 0.55 / 0.39 at 1.5 / 3 / 6 sigma, with `--curved` 0.00 / 0.89 / 0.70 (runs with coverage > 0.5: 0 / 4 / 3 of 8 -> 1 / 8 / 6 of 8 - the table lists them as "straight / both", note that the counts are over all 8 runs); curvature 8: median coverage 0.00 / 0.00 / 0.29 -> 0.00 / 0.44 / 0.73 (runs > 0.5: 0 / 1 / 3 -> 1 / 4 / 8).
  Faint (1.5 sigma) curved trails are not recovered.  At 6 sigma and curvature 8 the chain finder produced spurious extra segments in some runs (mean 492 stray masked pixels per run, i.e. a mask that is partly wrong - inspect before using on crowded fields).
* **Flickering trails** (SYNTHETIC, period 60 px, duty = fraction on): duty 0.7: covered at all amplitudes (7/8 found, coverage 0.98-0.99); duty 0.5: 3 and 6 sigma 7/8 (coverage 0.95-0.96), 1.5 sigma 2/8; duty 0.3: only 3/8 at 6 sigma and none below.  Measured duty (median) 0.41 / 0.51-0.59 / 0.61-0.68 for true 0.3 / 0.5 / 0.7: it is biased high because the bins are 12 px wide.
* **Contaminating flux** (`TRAIL_FLUX`; SYNTHETIC trails of 0.3-3 sigma, HUDF windows, 721 catalogue objects touched by a detected trail): for the 282 objects with true contamination (injected light inside the 2.5A x 2.5B ellipse) > 3 aperture sigma, model/truth median 1.05 (16-84 %: 0.90-1.33);
  the RMS of (model - truth) is 24.6 aperture sigma (dominated by large bright objects, where a few-percent mismatch of the profile is many sigma), 11.3 sigma for the others.  Median `TRAIL_FRAC` of the contaminated objects 0.062 (90th percentile 0.28).  The estimate is the trail model's light, not a measured quantity of the object itself, and ignores the local sky subtraction of the catalogue.
* **Dithered stack registration**: tested with three synthetic frames dithered by (0,0), (7.3,-4.6), (-5.1,9.4) px through their WCS (`test_stack_registers_dithered_frames_through_the_wcs`): `--register auto` reports `registered=2` and the stacked star peak is > 0.7 of the ideal single-frame peak (and at least 0.8 of the `none` stack), whereas `none` smears it.  Only integer-free shifts, no rotation, noise-free stars: a functional test, not a photometric accuracy measurement.  Superseded by the real-data run below.
* **Real dithered data: HST ACS/WFC F606W association `jc8m32010`** (GO-13498, MACS J0717.5+3745, 2014-11-11; four flc exposures j3q/j5q/j9q/jcq of 1207-1307 s, POSTARG dithers up to 2.4", Gaia-aligned WCS `IDC_4bb1536cj-FIT_REL_GAIAeDR3`), compared with the MAST drizzled `jc8m32010_drc.fits` and Gaia DR3; scripts and numbers in Astrafex Web `examples/hst_dither_stack/` (README).  Bugs found and fixed:
  1. HST flc/flt WCS have lookup-table distortions (D2IMARR / WCSDVARR): `WCS(header)` raises without the HDUList, `_wcs` returned None and the stack was written **unregistered** (`registered=0`) without a warning.  Now read with `fobj`; a warning is printed when the reference has no WCS.
  2. `wcs_resample` replaced masked (NaN) pixels by 0 inside the bilinear sum, pulling every neighbour of a masked pixel down (half the value for a half-pixel shift); now normalised interpolation (masked pixels excluded, output NaN where the valid weight < 0.5).
  3. The stack header kept the CPDIS/D2IM keywords of the reference chip but not the tables: the stack's WCS was unreadable; the tables are now appended to `trails_stack.fits`.
  4. MEF handling: only the first 2-D HDU was read (no DQ, no second chip, counts in electrons with different EXPTIME and sky levels mixed).  New `--all-chips`, `--dq-bits`, `--scale exptime`, `--sky median`; `trails_stack.json` records the per-chip sky, exposure time and trails.
  5. `sigclip` uses the MAD of the 3-4 frame values: on 4 pure-noise frames it rejects ~22 % of good pixels and it cannot reject a pixel hit by cosmic rays in 2 of 4 exposures (frequent in 1300 s ACS exposures).  New `--stack-method crrej`: AstroDrizzle driz_cr-style test against the ERR-array noise model with the derivative term, a minmed reference (minimum where the median is far above it, tolerance including the minimum image's derivative so star cores keep the median), neighbour growth, exposure-time weights.
  Tests: `plugins/trails/tests/test_dither_stack.py` (5).  Results of the crrej stack against the drc are in the Astrafex example README and user guide.
* **Held-out check 1 - SDSS run 94 camcol 3, 160 REAL frames (32 fields x ugriz)** with the straight detector as frozen at rev 58f69ec69 and the (untuned, first-guess) curved finder; multi-band consistency through the WCS of each frame (a feature seen in another band at SNR >= 6 is "static").
  Result: 37 of 160 frames had >= 1 straight detection (56 detections) and 10 had curved-chain segments (22): 60 of the 78 are present in other bands (static), 18 are single-band candidates.  **53 of the 56 straight detections were within 3 deg of 45/135 deg**: diffraction spikes of saturated stars (SDSS spike directions).
  That is, the straight finder's false-detection rate on SDSS was **23 % of frames (37/160)** before the change, not zero.  Looking at single-band candidates through an image-description tool (not by eye) also gave bright stars on the line.  After changing the spike rule (30 x / 2.5; this was done after seeing these frames, so it is no longer blind):
  11 straight detections in 11 frames (7 %) and 18 curved segments in 8 frames remain, 9 of the 11 straight at spike angles; I could not confirm any of them as a satellite trail.  No SDSS satellite-trail truth exists for these frames.
  **Visual inspection (2026-10-08)**: the same 160 frames were run again on master `8984e74` (detector identical to `ad876106`), `--curved`, 6 workers, 4593 s.  The count is again 29 (straight 11, curved 18) in 17 frames, and 9 of the 11 straight are within 1.2 deg of 45/135.  Straightening each segment and removing the brightest 3 % of columns (`plugins/trails/validation/sdss_ridge.py`; table in `trails_validation.md`) gives **16 diffraction spikes, 5 stars sitting on the segment, 8 with no continuous streak, and no confirmed satellite trail**.  The long straight detection at 159 deg (i field 165, 2155 px) has a ridge of 0.19 pixel-sigma.  The one faint coherent ridge (i field 158, 44.4 deg, 0.87 pixel-sigma, FWHM 6 px) is at a spike angle and passes through a star; z-band collapsed SNR along it is 5.0, under the static cut of 6.  The old heldout column headed `amp [sigma_pix]` is profile SNR, not per-pixel amplitude.
* **Held-out check 2 - HUDF F105W / F125W** (REAL, never used for tuning, same sky as F160W but independent pixels): 8 non-overlapping 700x700 windows (the 40 F160W tuning windows cover all fully covered F160W area, so no unseen F160W window exists), 0 detections (straight or curved).  This is a small sample (upper 95 % limit ~31 % per window).
* **Cost of the spike-rule change**: HUDF injection recall 35/40 -> 32/40 at 4 sigma and 26/40 -> 22/40 at 2 sigma; false trails stay 0/88.

## Limitations

## Limitations
* Curved trails are only handled as piecewise-linear chains (`--curved`, opt-in): recall 44-89 % at 3 sigma for 4-8 deg/100 px, none at 1.5 sigma, with occasional spurious segments; low-duty (<= 0.3) flickering trails are mostly missed; trails wider than ~40 px FWHM or shorter than 12 x their width are rejected on purpose.
* On real wide-field frames with saturated stars (SDSS) the false-detection rate was 23 % of frames (before the spike rule change) and 7 % after.  The 2026-10-08 inspection of the 29 that remain classified 16 as diffraction spikes, 5 as a star on the segment with no narrow ray, and 8 as having no continuous streak (table in `plugins/trails/validation/trails_validation.md`).  Trails crossing a bright star core may be read as spikes (lower recall).
* Sensitivity is set by the correlated noise and the galaxies: about 1-2 pixel-sigma peak amplitude (HUDF F160W: ~28 mag/arcsec^2) for 50 %; faint partial trails inside a galaxy disc (M51) are often missed.
* The 32 px background mesh subtracts part of very wide faint trails (> ~30 px) before detection; the mask width then comes from the fitted profile, not the full trail.
* Bleed/spike rejection is a heuristic; a genuine trail running along a column or through a saturated star core with falling brightness would be rejected (`--keep-bleeds`).  Spikes that stay below 30 x the run amplitude (`--spike-fac`) are reported as trails.
* Interpolation destroys objects on the trail and perturbs objects within ~15 px of it; use masks for photometry.  The stack step resamples frames and masks through the WCS (bilinear, so it correlates the noise slightly and smooths narrow features); frames without a usable WCS are assumed registered.
* The masks are binary; `TRAIL_FLUX`/`TRAIL_FRAC` model the trail's light in the object's aperture (median +5 % bias, 16-84 % spread 0.9-1.33 of the truth) but are not a measurement of the object's own corrected flux.
