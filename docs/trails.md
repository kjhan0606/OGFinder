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
(adds peak surface brightness and width in arcsec to `trails.json`).

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
   axes, and diffraction spikes at any angle: run through a core > 300 x the run amplitude whose brightness falls off along the line).
7. *Catalogue flags*: bit 1 = footprint (touch-k x A x B ellipse) reaches the mask band, bit 2 = centre inside the band, bit 4 = elongated (A/B >= 3) and aligned
   with the trail (a fragment of the trail itself); `TRAIL_DIST` = distance of the centre to the band edge (negative inside).

## Validation (code rev 58f69ec69; `plugins/trails/validation/`, tables in `trails_validation.md`, raw numbers reproducible with `trails_validate.py`)
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
  There is no labelled truth for this exposure, so no recall/precision number; an SDSS r-band frame set (run 94) gave detections in 2 of 6 frames that were **not inspected** (unknown whether real satellite/asteroid trails or artefacts).

## Limitations
* Straight trails only (curved, tumbling/flickering or strongly varying amplitude trails are partly found or cut into pieces); trails wider than ~40 px FWHM or shorter than 12 x their width are rejected on purpose.
* Sensitivity is set by the correlated noise and the galaxies: about 1-2 pixel-sigma peak amplitude (HUDF F160W: ~28 mag/arcsec^2) for 50 %; faint partial trails inside a galaxy disc (M51) are often missed.
* The 32 px background mesh subtracts part of very wide faint trails (> ~30 px) before detection; the mask width then comes from the fitted profile, not the full trail.
* Bleed/spike rejection is a heuristic; a genuine trail running along a column or through a saturated star core with falling brightness would be rejected (`--keep-bleeds`).  Spikes that stay below 300 x the run amplitude are reported as trails.
* Interpolation destroys objects on the trail and perturbs objects within ~15 px of it; use masks for photometry.  Stack uses registered frames; no resampling of masks between frames (use the mask reproject mode of the mask manager first).
* The masks are binary, the catalogue flag is geometric (it does not estimate the contaminating flux).
