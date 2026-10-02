# Extended morphology (plugin `morph_ext`, step "Extended Morphology")

Adds what `morphometry/` (CAS, Gini, M20, one Petrosian radius) lacked, in the new module `morphometry/extended.py` (pure functions) and the driver `ds9/library/ds9_morph_ext.py`.
The existing step and its columns (`CONC ASYM GINI M20 R_PETRO`) are unchanged.

## Measurements (columns `MX_*`)
| column | definition |
|---|---|
| `MX_RP`, `MX_RP_LO`, `MX_RP_HI` | Petrosian radii (semi-major axis, px) for eta = 0.2 / 0.1 / 0.3 (parameters): eta(r) = I(r) / <I>(<r) with I(r) from the 0.8 r - 1.25 r annulus (Bershady et al. 2000), fine log-spaced growth curve in elliptical apertures, linear interpolation of the crossing |
| `MX_FLUX_P`, `MX_MAG_P` | Petrosian flux = flux inside 2 R_P(eta) and its magnitude |
| `MX_R20 MX_R50 MX_R80 MX_R90` | radii containing 20/50/80/90 % of the Petrosian flux (from the growth curve) |
| `MX_CONC` | 5 log10(R80 / R20) |
| `MX_KRON_R`, `MX_KRON_MAG` | Kron radius (first moment, iterated once) and the magnitude in `kron-scale` x r_k (>= `kron-min`) elliptical aperture |
| `MX_SMOOTH` | smoothness / clumpiness S (Conselice 2003): sum(abs(I - I_boxcar) - abs(B - B_boxcar)) / sum(I) inside 1.5 R_P, boxcar width 0.25 R_P, central box excluded, B = Gaussian noise patch of the measured rms (noise correction) |
| `MX_GINI_P`, `MX_M20_P` | Gini and M20 on the Lotz et al. (2004) segmentation: pixels inside 1.5 R_P brighter than the surface brightness at R_P |
| `MX_FLAG` | 1 object too close to the image edge / too small, 2 no Petrosian radius found, 4 2 R_P beyond the search radius or image, 8 neighbours masked, 1024 no position |

Apertures follow the catalog ellipse (A_IMAGE, B_IMAGE, THETA_IMAGE) unless *Circular apertures*. Other catalog objects (> 3 % of the target flux) are masked with ellipses of `neighbour-radius` x A_IMAGE;
masked pixels are replaced by the mean of the unmasked pixels at the same elliptical radius (sep's area-ratio correction is not used: it is wrong when the mask sits in the outskirts).
Growth curves (`morph_ext_growth.tsv`: NUMBER R FLUX ETA) and a growth-curve / eta plot of the brightest objects (`morph_ext_curves.png`) are written to `<work>/morph_ext/`
(menu: *Growth Curves Plot…*, *Growth Curve Table…*). Parameters live in the dialog of the "Morph+" chip (Measure tab); catalog keys `morphext,*`.

    python ds9/library/ds9_morph_ext.py IMAGE --catalog TSV [--mask FITS] [--work DIR] [--eta 0.2] [--kron-scale 2.5] [--mask-neighbours] [--circular] ...

## Validation (analytic Sersic references: exact Petrosian eta from the incomplete gamma function, noise-free renderings, 401² images, F_tot = 1e5; `plugins/morph_ext/tests`)
| case | R_P (eta 0.1/0.2/0.3) rel. error | R20/R50/R80/R90 max rel. error | C error | F_P/F_tot (true) |
|---|---|---|---|---|
| n=1, re 15, q 1 | +0.0001 / 0.0000 / -0.0004 | 0.0010 | 0.000 | 0.9931 (0.9933) |
| n=1, re 15, q 0.6, PA 40 | +0.0004 / +0.0002 / +0.0001 | 0.0029 | -0.005 | 0.9930 (0.9933) |
| n=4, re 12, q 1 | +0.0004 / +0.0017 / +0.0024 | 0.0072 | -0.013 | 0.8169 (0.8167) |
| n=2, re 12, q 0.7, PA 100 | +0.0005 / +0.0009 / +0.0026 | 0.0062 | -0.012 | 0.9427 (0.9431) |

* Kron: the flux in 2.5 r_k equals the analytic Sersic aperture fraction at the measured r_k to 0.0005 (n=1: 0.9214 vs 0.9216; n=4: 0.7429 vs 0.7430; q 0.6 and 0.7 cases the same); r_k = 14.99 px (n=1, re 15).
* Noise (12 realisations, n=1, re 15, q 0.7): R_P bias 0.000 / scatter 0.001 at rms 0.5, bias 0.000 / scatter 0.013 at rms 8; smoothness with the noise correction 0.022 (rms 0.5) and 0.005 ± 0.006 (rms 8), without it 0.045 and 0.376.
* Smoothness: smooth disk 0.007, the same disk with 10 Gaussian clumps (20 % of the flux) 0.208.
* Gini_P / M20_P: n=1 0.437 / -1.77, n=4 0.575 / -2.35, n=4 with noise 0.577 / -2.40 (peak S/N 989).
* Companion in the outskirts (30 % of the flux, 14 px mask): R_P error +77.5 % unmasked, +0.1 % masked; F_P/F_tot 1.299 → 0.994 (true 0.993).
* Driver on a synthetic 400² field (4 galaxies, noise 1, F = 4e4): R_P errors +0.5 / -1.4 / +0.2 / +0.1 %, F_P -0.1…-1.5 %, R50 ≤ 2.2 %, S 0.01-0.04, no flags.

## Merger / interaction indicators (`morphometry/features.py`; step "Merger / Tidal Features", CLI `ds9_morph_ext.py --features-only`)
Added on top of the Petrosian-segmented Gini/M20 of the extended step (the step runs the same measurement, always masks the catalog neighbours, and returns only these columns; `--features` adds them to the `MX_*` columns of the normal run):

| column | meaning |
|---|---|
| `MX_ASYM` | rotational asymmetry (Conselice 2003) in 1.5 R_P, minimised over centres within +-2 px, minus the same quantity of a noise image (minimised over the same centres, normalised by the galaxy light) |
| `MX_NPEAK` | peaks of the Gaussian-smoothed light (sigma max(1.5, 0.08 R_P) px) above 5 sigma and 15 % of the brightest, >= 2.5 sigma apart (double nuclei / clumps) |
| `MX_LOTZ` | Lotz et al. 2008 class from the segmented G and M20: `merger` (G > -0.14 M20 + 0.33), `E/S0/Sa` (G > 0.14 M20 + 0.80), `Sb-Irr` |
| `MX_NSHELL`, `MX_SHELL_FRAC` | number of shell/arc components and their flux relative to F_P |
| `MX_NTAIL`, `MX_TAIL_FRAC`, `MX_TAIL_LEN` | tidal-tail components, flux fraction, length in R_P |
| `MX_TIDAL_FRAC`, `MX_RES_FRAC` | all positive significant residual flux / sum of the absolute residual in 0.4-4 R_P relative to F_P |
| `MX_NPAIR`, `MX_PAIR_SEP`, `MX_PAIR_RATIO`, `MX_PAIR_KIND` | catalog pairs: companions within `--pair-sep` (1.5) x (R_P + R_P,comp) and flux ratio within 1:10; separation (px) / flux ratio of the nearest; kind 0 / 1 minor (1:10) / 2 major (1:4) |
| `MX_MORPH_FLAGS` | bits: 1 G-M20 merger, 2 A >= 0.35, 4 several peaks, 8 shells, 16 tails, 32 minor pair, 64 major pair |

Shell / tail detection: the galaxy is modelled by the azimuthal median in elliptical annuli (catalog centre, q, PA; 1 px bins); the residual is smoothed (sigma max(1.5, 0.08 R_P)) and thresholded at `--feat-nsig` (3) times the larger of the analytic and the empirical noise
(rms of the *negative* half of the smoothed residual, so that correlated noise, clumps and arms - which produce negative as well as positive residuals - raise the threshold); connected positive components (>= 25 px) between 0.4 and 4 R_P are classified from their geometry about the centre:
**shell** = angular extent >= 45 deg, tangential direction (> 55 deg from radial), tangential/radial size >= 1.8, r in 0.5-4 R_P; **tail** = elongation >= 2.5, direction within 50 deg of radial, reaching > 1.1 R_P, length >= 0.8 R_P, starting inside 1.6 R_P;
components with > 25 % of the galaxy flux or within sigma + 3 px of masked neighbours are companions / halos and never shells or tails.

Validation (`plugins/morph_ext/validation/features_validate.py`, reports `features_report.json`, `features_real_injection_report.json`; tests `tests/test_features.py`):
* **Synthetic**, Sersic hosts (n 2-4, R_e 9-14 px, q 0.55-1, random PA, F = 1e5, noise sigma 1 / 4 / 12 per pixel), 60 realisations per cell; injected shells (two arcs at 1.3 and 1.9 R_e) and tails (3 R_e) with a fraction f of the host flux.
  False positives on smooth controls (shell or tail) 0 / 60 at every noise level; asymmetry of the controls 0.009 / 0.013 / 0.023 (95th percentile 0.017 / 0.021 / 0.030). Detection rate (shell / tail):
  sigma 1: f = 1 % 0.32 / 0.65, 2 % 0.93 / 0.98, 5 % 0.92 / 0.97, 10 % 0.85 / 0.93; sigma 4: 1 % 0 / 0.33, 2 % 0.93 / 0.98, 5 % 0.92 / 0.98; sigma 12: 2 % 0.15 / 0.53, 5 % 0.92 / 1.00, 10 % 0.92 / 0.97
  (the ~8 % misses are shells split into pieces or merged with the neighbouring component; recovered flux fraction of the shells 0.7-0.9 of the injected one, tails 0.9-0.98 at f >= 5 %). Nothing is detected at 0.5 %.
* **Mergers** (host + companion, sigma 4): equal mass at 1.6 and 2.5 R_e separation: 100 % flagged (A 0.42-0.46; G-M20 merger and A > 0.35), at 1.0 R_e (overlapping) 57 %; 1:2 at >= 1.6 R_e 100 % (A 0.31, flagged by G-M20 only: the A cut misses it);
  1:4 87-93 % (A 0.19-0.20, G-M20 only), at 1.0 R_e 3 %; 1:10 0 % at <= 1.6 R_e and 27 % at 2.5 R_e. The multiple-peak test finds the companion in 70-100 % of the cases at every ratio. A alone separates only mass ratios >= 1:2.
* **Real noise and crowding**: smooth Sersic hosts injected into the real HUDF F160W mosaic (random empty positions, ~1000 real catalog objects around them masked as neighbours; 30 per cell, noise correlated by drizzling): false positives 0 / 30 for shells and tails at 3, 4 and 5 sigma;
  detections at 3 sigma (shell / tail): 2 %: 0.73 / 0.70, 5 %: 0.87 / 0.70, 10 %: 0.87 / 0.63 (at 5 sigma: 2 %: 0.30 / 0.37, 5 %: 0.83 / 0.80), 1 %: 0.13 / 0.20.
* **Real galaxies** (HUDF F160W 1800x1800 crop, 133 measurable catalog objects with R_flux >= 6 px): G-M20 merger 9, A >= 0.35 2, several peaks 23, shells 29, tails 3, minor / major pairs 28 / 55 (default separation 1.5 x sum of R_P in a crowded field). The shell rate of 22 % is the arcs / clumps of the
  star-forming and irregular galaxies themselves (the smooth injected hosts give 0 %): the shell flag is meaningful for smooth (early-type) hosts, use it together with `MX_SMOOTH` / `MX_LOTZ`. On the 600x600 M51 field (NGC 5194 + 5195) the main galaxy is flagged for shells (its spiral arms), several peaks and a major pair.

## Limitations
Merger indicators: tails curving by more than ~50 deg from the radial direction are missed, spiral arms and star-forming clumps are shell candidates (see above), the pair criterion uses projected pixel separations (no redshift information: line-of-sight projections are not removed), the asymmetry cut alone misses mass ratios below 1:2, thresholds are not calibrated on observed merger samples. Apertures use the catalog ellipse (no iterative centre / shape refinement); the asymmetry index of the existing step is not changed (its background correction and centre search are the old ones);
smoothness uses a simulated Gaussian noise patch (correlated noise, e.g. drizzled data, under-corrects: use a real empty-sky patch via the Python API `noise_patch=`); R_P needs the profile to
reach eta < level inside the cutout (flag 2/4 otherwise); the growth curve is not PSF-corrected.
